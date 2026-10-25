package cli

// Verification of keyed-seed fingerprints by regeneration (#27).
//
// A keyed-seed run is an ordinary run whose seed is derived from the key,
// customer, dataset and salt (fingerprinting.DeriveSeed). To check whether
// suspect files come from a run, verify indexes up to MaxReads reads of each
// file, regenerates each candidate run from its configuration and derived
// seed one file at a time (nothing is written), and counts the suspect reads
// it reproduces: exactly, or by k-mer containment for trimmed, re-oriented
// or edited reads (fingerprinting.Matcher). Qualities and read names are not
// used, so rewritten qualities and renamed reads still match.

import (
	"bufio"
	"compress/gzip"
	"crypto/rand"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// DefaultVerifyMaxReads is the default number of reads examined per
// suspect file.
const DefaultVerifyMaxReads = 10000

// regenerationRequested reports whether opts select verification by
// regeneration rather than legacy quality-score extraction.
func (o VerifyOptions) regenerationRequested() bool {
	return len(o.Manifests) > 0 || o.ConfigPath != "" || o.Candidates != ""
}

// candidate is one run the suspect files may come from.
type candidate struct {
	label string // for reports
	// keyed runs derive their seed from the key; otherwise seed is the
	// run's recorded seed.
	keyed             bool
	customer, dataset string
	salt, seed        uint64
	config            types.SequencerConfig
	files             []genomics.RunFile // nil: the configured plan
	// fromManifest marks a candidate whose configuration comes from a
	// manifest (directly, or as the base of --candidates); inputs are the
	// input files that manifest records (see verify_inputs.go). A --config
	// candidate's configuration names its files.
	fromManifest bool
	inputs       runInputs
}

// candidateResult is the outcome of matching against one candidate.
type candidateResult struct {
	c       *candidate
	files   []fingerprinting.FileMatch
	regen   int   // files regenerated
	reads   int64 // reads regenerated
	elapsed time.Duration
}

// verifyByRegeneration implements VerifyFiles for keyed-seed fingerprints.
func verifyByRegeneration(w io.Writer, paths []string, opts VerifyOptions) error {
	cands, err := loadCandidates(opts)
	if err != nil {
		return err
	}
	if err := resolveCandidateInputs(w, cands, opts); err != nil {
		return err
	}
	var key []byte
	for _, c := range cands {
		if c.keyed {
			if key, err = fingerprinting.ResolveKey(opts.KeyFile); err != nil {
				return err
			}
			if key == nil {
				return fmt.Errorf("verifying a keyed-seed fingerprint needs the fingerprint key: pass --key-file or set %s", fingerprinting.KeyEnvVar)
			}
			break
		}
	}
	for _, c := range cands {
		if c.keyed {
			if c.seed, err = fingerprinting.DeriveSeed(key, c.customer, c.dataset, c.salt); err != nil {
				return fmt.Errorf("%s: %w", c.label, err)
			}
		}
	}

	maxReads := opts.MaxReads
	if maxReads <= 0 {
		maxReads = DefaultVerifyMaxReads
	}
	m := fingerprinting.NewMatcher()
	for _, p := range paths {
		if err := indexSuspect(m, p, maxReads); err != nil {
			return err
		}
	}
	if m.Reads() == 0 {
		return fmt.Errorf("no reads in %s", strings.Join(paths, ", "))
	}

	regen := genomics.RegenOptions{MaxFiles: opts.MaxFiles, MaxReads: opts.MaxRegenReads}
	results := make([]candidateResult, 0, len(cands))
	for _, c := range cands {
		r, err := matchCandidate(m, c, c.seed, regen)
		if err != nil {
			return err
		}
		fmt.Fprintf(w, "   regenerated %s: %d files, %d reads in %s\n", c.label, r.regen, r.reads, r.elapsed.Round(time.Millisecond))
		results = append(results, r)
	}
	var decoy *candidateResult
	if opts.Decoy {
		seed, err := randomSeed()
		if err != nil {
			return err
		}
		r, err := matchCandidate(m, cands[0], seed, regen)
		if err != nil {
			return err
		}
		fmt.Fprintf(w, "   decoy run (%s with a random seed): %d files, %d reads in %s\n", cands[0].label, r.regen, r.reads, r.elapsed.Round(time.Millisecond))
		decoy = &r
	}

	failed := 0
	for i, p := range paths {
		if !reportFile(w, p, i, results, decoy) {
			failed++
		}
	}
	if failed > 0 {
		return fmt.Errorf("%d of %d files match no candidate", failed, len(paths))
	}
	return nil
}

// matchCandidate regenerates c with seed and matches it against m.
func matchCandidate(m *fingerprinting.Matcher, c *candidate, seed uint64, regen genomics.RegenOptions) (candidateResult, error) {
	cfg := c.config
	cfg.Seed = seed
	cfg.Fingerprinting = types.FingerprintConfig{}
	cfg.Enterprise.Fingerprinting = types.FingerprintConfig{}
	cfg.OutputDir, cfg.Workers = "", 0
	regen.Files = c.files
	// resolveCandidateInputs has matched the other input files to the
	// manifest; checking them again would hash them for every candidate.
	regen.Variants = c.inputs.variants
	start := time.Now()
	m.Reset()
	e := genomics.NewEngineLog(cfg, io.Discard)
	defer e.Close()
	n, err := e.Regenerate(regen, func(_ genomics.RunFile, seq []byte) error {
		m.Observe(seq)
		return nil
	})
	if err != nil {
		return candidateResult{}, fmt.Errorf("%s: %w", c.label, err)
	}
	reads, _ := m.Regenerated()
	return candidateResult{c: c, files: m.Results(), regen: n, reads: reads, elapsed: time.Since(start)}, nil
}

// minMatches is the matched reads a verdict needs (all of them for files
// with fewer reads).
const minMatches = 3

// reportFile writes the verdict for suspect file i and reports whether a
// candidate matches it. A candidate matches when at least minMatches reads
// (or every read of a smaller file) match and the matches exceed ten times
// both the uniform-random chance bound and the decoy run's matches.
func reportFile(w io.Writer, path string, i int, results []candidateResult, decoy *candidateResult) bool {
	type scored struct {
		r  candidateResult
		f  fingerprinting.FileMatch
		ok bool
	}
	var all []scored
	best := -1
	for _, r := range results {
		f := r.files[i]
		background := f.ExpectedChance
		if decoy != nil {
			background = max(background, float64(decoy.files[i].Matched()))
		}
		ok := f.Checked > 0 && f.Matched() >= min(minMatches, f.Checked) && float64(f.Matched()) > 10*background
		all = append(all, scored{r, f, ok})
		if best < 0 || f.Matched() > all[best].f.Matched() {
			best = len(all) - 1
		}
	}
	b := all[best]
	decoyNote := ""
	if decoy != nil {
		decoyNote = fmt.Sprintf("; decoy run %d", decoy.files[i].Matched())
	}
	if b.ok {
		fmt.Fprintf(w, "✅ %s: matches %s: %s%s\n", path, b.r.c.label, describe(b.f), decoyNote)
	} else if len(all) == 1 {
		fmt.Fprintf(w, "❌ %s: does not match %s: %s%s\n", path, b.r.c.label, describe(b.f), decoyNote)
	} else {
		fmt.Fprintf(w, "❌ %s: matches none of %d candidates (best %s: %s%s)\n", path, len(all), b.r.c.label, describe(b.f), decoyNote)
	}
	for j, s := range all {
		if j == best || len(all) == 1 {
			continue
		}
		mark := "  "
		if s.ok {
			mark = "also matches"
		}
		fmt.Fprintf(w, "    %s %s: %s\n", mark, s.r.c.label, describe(s.f))
	}
	return b.ok
}

// describe summarises a file's matches.
func describe(f fingerprinting.FileMatch) string {
	return fmt.Sprintf("%d of %d reads (%d exact, %d by %d-mer containment); expected by chance ≤ %.1g",
		f.Matched(), f.Checked, f.Exact, f.Contained, fingerprinting.MatchK, f.ExpectedChance)
}

// randomSeed draws a non-zero seed for the decoy run.
func randomSeed() (uint64, error) {
	var b [8]byte
	for {
		if _, err := rand.Read(b[:]); err != nil {
			return 0, fmt.Errorf("drawing a decoy seed: %w", err)
		}
		if s := binary.LittleEndian.Uint64(b[:]); s != 0 {
			return s, nil
		}
	}
}

// loadCandidates builds the candidate runs from opts.
func loadCandidates(opts VerifyOptions) ([]*candidate, error) {
	var cands []*candidate
	if len(opts.Manifests) > 0 && (opts.Customer != "" || opts.Dataset != "" || opts.SaltSet) {
		return nil, errors.New("--customer, --dataset and --salt describe a candidate for --config; a --manifest records its own")
	}
	var base *types.SequencerConfig
	var baseManifest *candidate // the manifest base is taken from
	for _, p := range opts.Manifests {
		c, err := manifestCandidate(p)
		if err != nil {
			return nil, err
		}
		if base == nil {
			cfg := c.config
			base, baseManifest = &cfg, c
		}
		if opts.Candidates == "" {
			cands = append(cands, c)
		}
	}
	single := opts.Customer != "" || opts.Dataset != "" || opts.SaltSet
	switch {
	case opts.Candidates != "" && single:
		return nil, errors.New("--candidates lists the customers, datasets and salts: do not also pass --customer, --dataset or --salt")
	case opts.Candidates == "" && opts.ConfigPath != "" && len(opts.Manifests) > 0:
		return nil, errors.New("a --manifest carries its run's configuration: pass --config only with --customer/--dataset/--salt or --candidates")
	case opts.Candidates == "" && opts.ConfigPath != "" && (opts.Customer == "" || opts.Dataset == "" || !opts.SaltSet):
		return nil, errors.New("--config needs --customer, --dataset and --salt (from the run's manifest or log), or --candidates")
	}
	if opts.ConfigPath != "" {
		cfg, err := NewConfigLoader(false).LoadSequencerConfig(opts.ConfigPath)
		if err != nil {
			return nil, err
		}
		base, baseManifest = cfg, nil
	}
	switch {
	case opts.Candidates != "":
		if base == nil {
			return nil, errors.New("--candidates needs the run configuration: pass --config (or a --manifest of the same configuration)")
		}
		list, err := readCandidates(opts.Candidates, *base)
		if err != nil {
			return nil, err
		}
		if baseManifest != nil {
			// The manifest's configuration carries no input paths: the
			// candidates read the input files it records.
			for _, c := range list {
				c.fromManifest, c.inputs = true, baseManifest.inputs
			}
		}
		cands = append(cands, list...)
	case opts.ConfigPath != "":
		cands = append(cands, &candidate{
			label: candidateLabel(opts.Customer, opts.Dataset, opts.Salt), keyed: true,
			customer: opts.Customer, dataset: opts.Dataset, salt: opts.Salt, config: *base,
		})
	}
	if len(cands) == 0 {
		return nil, errors.New("no candidate runs: pass --manifest, or --config with --customer, --dataset and --salt, or --candidates")
	}
	return cands, nil
}

// candidateLabel names a keyed-seed candidate in reports.
func candidateLabel(customer, dataset string, salt uint64) string {
	return fmt.Sprintf("customer=%s dataset=%s salt=%d", customer, dataset, salt)
}

// manifestCandidate reads the run described by the manifest at path.
func manifestCandidate(path string) (*candidate, error) {
	m, err := manifest.Read(path)
	if err != nil {
		return nil, fmt.Errorf("reading manifest: %w", err)
	}
	if m.Signal != nil {
		return nil, fmt.Errorf("%s: a signal run's BLOW5 output cannot be verified by read regeneration", path)
	}
	c := &candidate{config: m.Config, fromManifest: true, inputs: runInputs{
		variants: m.Variants, reference: m.Reference, model: m.HMMModel, readLengths: m.ReadLengths,
	}}
	switch fp := m.Fingerprint; {
	case fp != nil && fp.Scheme == fingerprinting.KeyedSeedScheme:
		c.keyed, c.customer, c.dataset, c.salt = true, fp.Customer, fp.Dataset, fp.Salt
		c.label = candidateLabel(fp.Customer, fp.Dataset, fp.Salt)
	case fp != nil:
		return nil, fmt.Errorf("%s: unknown fingerprint scheme %q", path, fp.Scheme)
	case m.Seed != 0:
		c.seed = m.Seed
		c.label = fmt.Sprintf("seed=%d (not fingerprinted)", m.Seed)
	default:
		return nil, fmt.Errorf("%s: the manifest records neither a fingerprint nor a seed", path)
	}
	c.label = filepath.Base(filepath.Dir(path)) + "/" + filepath.Base(path) + " " + c.label
	mates := 1
	if m.Config.PairedEnd {
		mates = 2
	}
	unit := map[int]int{}
	for _, o := range m.Outputs {
		if o.Role != manifest.RoleFASTQ {
			continue
		}
		fc := max(o.FlowCell, 1)
		n := unit[fc]
		unit[fc]++
		c.files = append(c.files, genomics.RunFile{FlowCell: fc, Batch: n/mates + 1, Mate: n%mates + 1, Reads: int(o.Reads)})
	}
	if len(c.files) == 0 {
		return nil, fmt.Errorf("%s: the manifest lists no FASTQ outputs", path)
	}
	return c, nil
}

// readCandidates parses a candidates file: one candidate per line,
// "customer<TAB>dataset<TAB>salt"; blank lines, lines starting with '#' and a
// "customer dataset salt" header line are skipped.
func readCandidates(path string, base types.SequencerConfig) ([]*candidate, error) {
	f, err := os.Open(path) // #nosec G304 -- operator-supplied candidates file
	if err != nil {
		return nil, fmt.Errorf("candidates: %w", err)
	}
	defer f.Close()
	var out []*candidate
	sc := bufio.NewScanner(f)
	for line := 1; sc.Scan(); line++ {
		text := strings.TrimRight(sc.Text(), "\r")
		if strings.TrimSpace(text) == "" || strings.HasPrefix(text, "#") {
			continue
		}
		fields := strings.Split(text, "\t")
		if len(fields) != 3 {
			return nil, fmt.Errorf("%s:%d: want customer<TAB>dataset<TAB>salt, got %d fields", path, line, len(fields))
		}
		if len(out) == 0 && strings.EqualFold(fields[0], "customer") && strings.EqualFold(fields[2], "salt") {
			continue
		}
		salt, err := strconv.ParseUint(strings.TrimSpace(fields[2]), 10, 64)
		if err != nil {
			return nil, fmt.Errorf("%s:%d: salt %q: %w", path, line, fields[2], err)
		}
		cust, ds := strings.TrimSpace(fields[0]), strings.TrimSpace(fields[1])
		out = append(out, &candidate{
			label: candidateLabel(cust, ds, salt), keyed: true,
			customer: cust, dataset: ds, salt: salt, config: base,
		})
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("%s: no candidates", path)
	}
	return out, nil
}

// indexSuspect adds the first maxReads reads of the FASTQ file at path
// (plain or .gz) to m.
func indexSuspect(m *fingerprinting.Matcher, path string, maxReads int) error {
	f, err := os.Open(path) // #nosec G304 -- user-supplied input files are the point of this command
	if err != nil {
		return err
	}
	defer f.Close()
	var r io.Reader = f
	if isGzip(path) {
		zr, err := gzip.NewReader(f)
		if err != nil {
			return fmt.Errorf("%s: %w", path, err)
		}
		defer zr.Close()
		r = zr
	}
	id := m.AddFile(path)
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 1<<20), 1<<30)
	n := 0
	for line := 0; sc.Scan() && n < maxReads; line++ {
		text := sc.Bytes()
		switch line % 4 {
		case 0:
			if len(text) == 0 || text[0] != '@' {
				return fmt.Errorf("%s: line %d: not a FASTQ header", path, line+1)
			}
		case 1:
			m.AddRead(id, text)
			n++
		}
	}
	if err := sc.Err(); err != nil {
		return fmt.Errorf("%s: %w", path, err)
	}
	return nil
}
