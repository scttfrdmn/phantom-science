package cli

// Verification of keyed-seed fingerprints by regeneration (#27, #31).
//
// A keyed-seed run is an ordinary run whose seed is derived from the key,
// customer, dataset and salt (fingerprinting.DeriveSeed). To check whether
// suspect files come from a run, verify indexes up to MaxReads reads of each
// file, regenerates each candidate run from its records and derived seed
// (nothing is written), and counts the suspect reads it reproduces.
//
// Suspect FASTQ files (from generate, scenario or transcriptome runs) are
// matched by bases: exactly, or by k-mer containment for trimmed,
// re-oriented or edited reads (fingerprinting.Matcher). Qualities and read
// names are not used, so rewritten qualities and renamed reads still match.
// Suspect BLOW5 or SLOW5 files (from signal runs) are matched by raw signal:
// exactly, or by windows of consecutive samples
// (fingerprinting.SignalMatcher). Read IDs, calibration and auxiliary
// fields are not used, so renamed, re-packed or trimmed reads still match.
// A signal verdict is always measured against a decoy run (the first signal
// candidate regenerated with a random seed), because how much of the seed a
// window of samples carries depends on the profile's noise.

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"context"
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
	"github.com/scttfrdmn/phantom-science/internal/genomics/signalrun"
	"github.com/scttfrdmn/phantom-science/internal/genomics/slow5"
	"github.com/scttfrdmn/phantom-science/internal/genomics/transcriptome"
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

// runKind is the command a candidate run was made by.
type runKind int

const (
	kindGenerate runKind = iota
	kindScenario
	kindSignal
	kindTranscriptome
)

// signal reports whether the run's output is raw signal (BLOW5) rather
// than reads (FASTQ).
func (k runKind) signal() bool { return k == kindSignal }

// alwaysDecoy reports whether verdicts on the kind's output are always
// measured against a decoy run. Scenario and transcriptome reads are
// sampled from sequence every run of the configuration shares (the
// scenario window or reference, the transcripts), so their k-mers recur
// between runs and the uniform-random chance bound does not hold; raw
// signal has no chance bound at all.
func (k runKind) alwaysDecoy() bool { return k != kindGenerate }

// candidate is one run the suspect files may come from.
type candidate struct {
	label string // for reports
	kind  runKind
	// keyed runs derive their seed from the key; otherwise seed is the
	// run's recorded seed.
	keyed             bool
	customer, dataset string
	salt, seed        uint64
	// config is the sequencer configuration of a generate or scenario run.
	config types.SequencerConfig
	files  []genomics.RunFile // generate: nil for the configured plan
	// fromManifest marks a candidate whose configuration comes from a
	// manifest (directly, or as the base of --candidates); inputs are the
	// input files that manifest records (see verify_inputs.go). A --config
	// candidate's configuration names its files.
	fromManifest bool
	inputs       runInputs
	// The records of scenario, signal and transcriptome runs, and the
	// input files verify resolved for them (verify_kinds.go).
	scenario *manifest.Scenario
	sig      *manifest.Signal
	sigSpec  signalrun.RegenSpec
	tx       *manifest.Transcriptome
	txSpec   transcriptome.RunSpec
}

// candidateResult is the outcome of matching against one candidate.
type candidateResult struct {
	c       *candidate
	files   []fingerprinting.FileMatch // per suspect file of the candidate's kind
	regen   int                        // files regenerated
	reads   int64                      // reads regenerated
	elapsed time.Duration
}

// suspect is one suspect file: FASTQ reads or raw signal, and its index in
// the matcher of its kind.
type suspect struct {
	path   string
	signal bool
	index  int
}

// matchers hold the indexed suspect reads and signal.
type matchers struct {
	reads  *fingerprinting.Matcher
	signal *fingerprinting.SignalMatcher
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
	mm := matchers{reads: fingerprinting.NewMatcher(), signal: fingerprinting.NewSignalMatcher()}
	suspects := make([]suspect, 0, len(paths))
	for _, p := range paths {
		s, err := indexSuspectFile(mm, p, maxReads)
		if err != nil {
			return err
		}
		suspects = append(suspects, s)
	}
	haveReads, haveSignal := mm.reads.Reads() > 0, mm.signal.Reads() > 0
	if !haveReads && !haveSignal {
		return fmt.Errorf("no reads in %s", strings.Join(paths, ", "))
	}
	readCands, signalCands := 0, 0
	for _, c := range cands {
		if c.kind.signal() {
			signalCands++
		} else {
			readCands++
		}
	}
	for _, s := range suspects {
		switch {
		case s.signal && signalCands == 0:
			return fmt.Errorf("%s holds raw signal (SLOW5/BLOW5), but no candidate is a signal run: pass the manifest of the signal run", s.path)
		case !s.signal && readCands == 0:
			return fmt.Errorf("%s holds reads (FASTQ), but every candidate is a signal run, whose output is BLOW5: pass the manifest of the run that wrote the reads", s.path)
		}
	}

	regen := genomics.RegenOptions{MaxFiles: opts.MaxFiles, MaxReads: opts.MaxRegenReads}
	results := make([]candidateResult, 0, len(cands))
	decoyBase := map[runKind]*candidate{}
	for _, c := range cands {
		if (c.kind.signal() && !haveSignal) || (!c.kind.signal() && !haveReads) {
			fmt.Fprintf(w, "   skipped %s: no suspect file holds its kind of output\n", c.label)
			continue
		}
		if decoyBase[c.kind] == nil && (opts.Decoy || c.kind.alwaysDecoy()) {
			decoyBase[c.kind] = c
		}
		r, err := matchCandidate(mm, c, c.seed, regen)
		if err != nil {
			return err
		}
		fmt.Fprintf(w, "   regenerated %s: %d files, %d reads in %s\n", c.label, r.regen, r.reads, r.elapsed.Round(time.Millisecond))
		results = append(results, r)
	}
	// Decoy runs: the first candidate of each kind with a random seed, for
	// every kind with --decoy and always for the kinds whose output is
	// sampled from shared sequence or whose chance of matching has no
	// bound (see alwaysDecoy). A candidate is compared with the decoy of its
	// kind.
	decoys := map[runKind]*candidateResult{}
	for _, k := range []runKind{kindGenerate, kindScenario, kindTranscriptome, kindSignal} {
		base := decoyBase[k]
		if base == nil {
			continue
		}
		seed, err := randomSeed()
		if err != nil {
			return err
		}
		r, err := matchCandidate(mm, base, seed, regen)
		if err != nil {
			return err
		}
		fmt.Fprintf(w, "   decoy run (%s with a random seed): %d files, %d reads in %s\n", base.label, r.regen, r.reads, r.elapsed.Round(time.Millisecond))
		decoys[k] = &r
	}

	failed := 0
	for _, s := range suspects {
		var rs []candidateResult
		for _, r := range results {
			if r.c.kind.signal() == s.signal {
				rs = append(rs, r)
			}
		}
		if !reportFile(w, s, rs, decoys) {
			failed++
		}
	}
	if failed > 0 {
		return fmt.Errorf("%d of %d files match no candidate", failed, len(paths))
	}
	return nil
}

// errRegenStop ends a regeneration at the MaxReads limit.
var errRegenStop = errors.New("regeneration limit")

// matchCandidate regenerates c with seed and matches it against the
// matcher of its kind.
func matchCandidate(mm matchers, c *candidate, seed uint64, regen genomics.RegenOptions) (candidateResult, error) {
	start := time.Now()
	var files int
	var err error
	if c.kind.signal() {
		sm := mm.signal
		sm.Reset()
		spec := c.sigSpec
		spec.MaxReads = regen.MaxReads
		_, err = signalrun.Regenerate(spec, seed, func(_ string, raw []int16) error {
			sm.Observe(raw)
			return nil
		})
		files = len(spec.Inputs)
		if err != nil {
			return candidateResult{}, fmt.Errorf("%s: %w", c.label, err)
		}
		reads, _ := sm.Regenerated()
		return candidateResult{c: c, files: sm.Results(), regen: files, reads: reads, elapsed: time.Since(start)}, nil
	}

	m := mm.reads
	m.Reset()
	var observed int64
	observe := func(seq []byte) error {
		m.Observe(seq)
		if observed++; regen.MaxReads > 0 && observed >= regen.MaxReads {
			return errRegenStop
		}
		return nil
	}
	switch c.kind {
	case kindGenerate:
		cfg := regenConfig(c.config, seed)
		regen.Files = c.files
		// resolveCandidateInputs has matched the other input files to the
		// manifest; checking them again would hash them for every candidate.
		regen.Variants = c.inputs.variants
		e := genomics.NewEngineLog(cfg, io.Discard)
		files, err = e.Regenerate(regen, func(_ genomics.RunFile, seq []byte) error {
			m.Observe(seq)
			return nil
		})
		e.Close()
	case kindScenario:
		_, err = scenarioRegenerate(context.Background(), regenConfig(c.config, seed), c.scenario, func(_ string, seq []byte) error { return observe(seq) })
		files = len(c.scenario.Members)
	case kindTranscriptome:
		_, err = transcriptome.Regenerate(c.txSpec, seed, observe)
		files = 1
	}
	if errors.Is(err, errRegenStop) {
		err = nil
	}
	if err != nil {
		return candidateResult{}, fmt.Errorf("%s: %w", c.label, err)
	}
	reads, _ := m.Regenerated()
	return candidateResult{c: c, files: m.Results(), regen: files, reads: reads, elapsed: time.Since(start)}, nil
}

// regenConfig returns c set up to regenerate a run with seed: fingerprinting
// off (the seed is already derived) and no output settings.
func regenConfig(c types.SequencerConfig, seed uint64) types.SequencerConfig {
	c.Seed = seed
	c.Fingerprinting = types.FingerprintConfig{}
	c.Enterprise.Fingerprinting = types.FingerprintConfig{}
	c.OutputDir, c.Workers = "", 0
	return c
}

// minMatches is the matched reads a verdict needs (all of them for files
// with fewer reads).
const minMatches = 3

// attributed reports whether f, a suspect file's matches against a
// candidate, attributes the file to it given the decoy run's matches d (nil
// without a decoy). Each test compares a count of matched reads with a
// background and needs at least minMatches reads (or every read of a
// smaller file) and more than ten times the background; any test suffices:
//
//   - all matches, against the larger of the uniform-random chance bound
//     and the decoy run's matches;
//   - with a decoy run, exact matches against the decoy run's exact
//     matches, and whole matches (every sampled k-mer or window found;
//     FileMatch.Whole) against the decoy run's whole matches, each
//     background at least 1 (so at least 11 reads).
//
// The decoy tests are needed because reads sampled from shared sequence (a
// reference, transcripts) are contained in any run of the configuration,
// so half their k-mers are found in any run, while reproducing a read
// exactly, or every sampled k-mer of a trimmed read with its sequencing
// errors, takes the run that made it. The floor of 1 allows for the few
// reads two runs share by chance: the decoy run measures that rate from
// one sample.
func attributed(f fingerprinting.FileMatch, d *fingerprinting.FileMatch) bool {
	if f.Checked == 0 {
		return false
	}
	test := func(n int, background float64, need int) bool {
		return n >= need && float64(n) > 10*background
	}
	background := f.ExpectedChance
	if d != nil {
		background = max(background, float64(d.Matched()))
	}
	if test(f.Matched(), background, min(minMatches, f.Checked)) {
		return true
	}
	if d == nil {
		return false
	}
	return test(f.Exact, max(f.ExpectedChance, float64(d.Exact), 1), minMatches) ||
		test(f.Whole, max(f.ExpectedChance, float64(d.Whole), 1), minMatches)
}

// reportFile writes the verdict for suspect file s and reports whether a
// candidate matches it (see attributed). decoys holds the decoy run of each
// kind of candidate.
func reportFile(w io.Writer, s suspect, results []candidateResult, decoys map[runKind]*candidateResult) bool {
	type scored struct {
		r  candidateResult
		f  fingerprinting.FileMatch
		d  *fingerprinting.FileMatch
		ok bool
	}
	var all []scored
	best := -1
	for _, r := range results {
		f := r.files[s.index]
		var d *fingerprinting.FileMatch
		if dr := decoys[r.c.kind]; dr != nil {
			d = &dr.files[s.index]
		}
		all = append(all, scored{r, f, d, attributed(f, d)})
		if b := len(all) - 1; best < 0 || all[b].ok && !all[best].ok ||
			all[b].ok == all[best].ok && f.Matched() > all[best].f.Matched() {
			best = b
		}
	}
	path := s.path
	describe := describeReads
	if s.signal {
		describe = describeSignal
	}
	b := all[best]
	decoyNote := ""
	if b.d != nil {
		decoyNote = fmt.Sprintf("; decoy run %d (%d exact, %d whole)", b.d.Matched(), b.d.Exact, b.d.Whole)
	}
	if b.ok {
		fmt.Fprintf(w, "✅ %s: matches %s: %s%s\n", path, b.r.c.label, describe(b.f), decoyNote)
	} else if len(all) == 1 {
		fmt.Fprintf(w, "❌ %s: does not match %s: %s%s\n", path, b.r.c.label, describe(b.f), decoyNote)
	} else {
		fmt.Fprintf(w, "❌ %s: matches none of %d candidates (best %s: %s%s)\n", path, len(all), b.r.c.label, describe(b.f), decoyNote)
	}
	for j, sc := range all {
		if j == best || len(all) == 1 {
			continue
		}
		mark := "  "
		if sc.ok {
			mark = "also matches"
		}
		fmt.Fprintf(w, "    %s %s: %s\n", mark, sc.r.c.label, describe(sc.f))
	}
	return b.ok
}

// describeReads summarises a FASTQ file's matches.
func describeReads(f fingerprinting.FileMatch) string {
	return fmt.Sprintf("%d of %d reads (%d exact, %d by %d-mer containment); %d whole; expected by chance ≤ %.1g",
		f.Matched(), f.Checked, f.Exact, f.Contained, fingerprinting.MatchK, f.Whole, f.ExpectedChance)
}

// describeSignal summarises a signal file's matches.
func describeSignal(f fingerprinting.FileMatch) string {
	return fmt.Sprintf("%d of %d reads by raw signal (%d exact, %d by %d-sample windows); %d whole",
		f.Matched(), f.Checked, f.Exact, f.Contained, fingerprinting.SignalWindow, f.Whole)
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
	var base *candidate // the run --candidates are variants of
	for _, p := range opts.Manifests {
		c, err := manifestCandidate(p)
		if err != nil {
			return nil, err
		}
		if base == nil {
			base = c
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
		base = &candidate{kind: kindGenerate, config: *cfg}
	}
	switch {
	case opts.Candidates != "":
		if base == nil {
			return nil, errors.New("--candidates needs the run configuration: pass --config (or a --manifest of the same configuration)")
		}
		list, err := readCandidates(opts.Candidates, base)
		if err != nil {
			return nil, err
		}
		cands = append(cands, list...)
	case opts.ConfigPath != "":
		c := *base
		c.label, c.keyed = candidateLabel(opts.Customer, opts.Dataset, opts.Salt), true
		c.customer, c.dataset, c.salt = opts.Customer, opts.Dataset, opts.Salt
		cands = append(cands, &c)
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
	switch {
	case m.Signal != nil:
		c.kind, c.sig = kindSignal, m.Signal
		c.fromManifest = false // its inputs are not sequencer-config inputs
		if len(m.Signal.Inputs) == 0 {
			return nil, fmt.Errorf("%s: the signal manifest lists no input reads", path)
		}
		return c, nil
	case m.Transcriptome != nil:
		c.kind, c.tx = kindTranscriptome, m.Transcriptome
		c.fromManifest = false
		return c, nil
	case m.Scenario != nil:
		c.kind, c.scenario = kindScenario, m.Scenario
		return c, nil
	}
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
// "customer<TAB>dataset<TAB>salt", each a copy of base with that keyed
// seed; blank lines, lines starting with '#' and a "customer dataset salt"
// header line are skipped.
func readCandidates(path string, base *candidate) ([]*candidate, error) {
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
		c := *base
		c.label, c.keyed = candidateLabel(cust, ds, salt), true
		c.customer, c.dataset, c.salt, c.seed = cust, ds, salt, 0
		out = append(out, &c)
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("%s: no candidates", path)
	}
	return out, nil
}

// indexSuspectFile adds the first maxReads reads of the file at path to the
// matcher of its kind: raw signal for SLOW5 and BLOW5 files (recognised by
// their first bytes), reads for FASTQ (plain or .gz).
func indexSuspectFile(mm matchers, path string, maxReads int) (suspect, error) {
	f, err := os.Open(path) // #nosec G304 -- user-supplied input files are the point of this command
	if err != nil {
		return suspect{}, err
	}
	defer f.Close()
	br := bufio.NewReaderSize(f, 1<<20)
	head, _ := br.Peek(6)
	if bytes.HasPrefix(head, []byte("BLOW5")) || bytes.HasPrefix(head, []byte("#slow5")) {
		return indexSignal(mm.signal, path, br, maxReads)
	}
	s := suspect{path: path, index: mm.reads.AddFile(path)}
	return s, indexReads(mm.reads, s.index, path, br, maxReads)
}

// indexSignal adds the raw signal of the first maxReads reads of a SLOW5 or
// BLOW5 stream to sm.
func indexSignal(sm *fingerprinting.SignalMatcher, path string, r io.Reader, maxReads int) (suspect, error) {
	rd, err := slow5.NewReader(r)
	if err != nil {
		return suspect{}, fmt.Errorf("%s: %w", path, err)
	}
	defer rd.Close()
	s := suspect{path: path, signal: true, index: sm.AddFile(path)}
	for n := 0; n < maxReads; n++ {
		rec, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return suspect{}, fmt.Errorf("%s: %w", path, err)
		}
		sm.AddRead(s.index, rec.RawSignal)
	}
	return s, nil
}

// indexReads adds the first maxReads reads of the FASTQ stream r (gzip when
// path ends in .gz) to m as file id.
func indexReads(m *fingerprinting.Matcher, id int, path string, r io.Reader, maxReads int) error {
	if isGzip(path) {
		zr, err := gzip.NewReader(r)
		if err != nil {
			return fmt.Errorf("%s: %w", path, err)
		}
		defer zr.Close()
		r = zr
	}
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
