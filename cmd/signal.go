package cmd

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"crypto/sha1" // #nosec G505 -- UUIDv5 name hashing (RFC 4122), not security
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"fmt"
	"hash"
	"hash/fnv"
	"io"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore"
	"github.com/scttfrdmn/phantom-science/internal/genomics/signal"
	"github.com/scttfrdmn/phantom-science/internal/genomics/slow5"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"github.com/spf13/cobra"
)

// SignalTruthName is the truth FASTA the signal command writes: one record
// per read, named by its BLOW5 read ID, holding the sequence its signal was
// synthesised from.
const SignalTruthName = "truth.fa"

var signalCmd = &cobra.Command{
	Use:   "signal",
	Short: "Synthesise raw nanopore signal (BLOW5) from reads",
	Long: `Turn reads (FASTQ or FASTA, plain or gzip) into raw nanopore current traces
written as BLOW5 (#14), for basecaller benchmarking: basecall the BLOW5 and
align the calls to truth.fa.

Signal is synthesised from a k-mer pore model (--pore-model) with the
--profile chemistry (squigulator's algorithm, internal/genomics/signal). ONT
tables are not shipped: fetch one with scripts/fetch-pore-model.sh. Two-column
ONT kmer_models tables are in standardised units and are scaled to pA with
the profile's scaling, every k-mer getting the profile's model_stdv; tables
with a level_stdv column are taken as pA and keep their per-k-mer stdv
(--model-units overrides the detection). Level smoothing between k-mers is
--profile-set level_tau=SAMPLES (0, a step function, by default).

Source sequence. Each read's signal is synthesised from its sequence as
written in the input. With --truth FASTA (or FASTQ), a read whose ID (or
mapped UUID) is in the truth file is synthesised from the truth sequence
instead, e.g. the error-free molecule, so that basecalling identity measures
the basecaller rather than the input's simulated errors. truth.fa in --out
always holds the sequence each signal was made from.

Read IDs. UUID read IDs (as written by generate --platform ont) are kept;
any other ID is mapped to a name-based UUID (version 5, SHA-1 of the ID in a
fixed namespace), and the original ID is recorded in truth.fa as
"source=<id>". BLOW5, pod5 and Dorado require UUID read IDs.

Run metadata. MinKNOW header fields (runid, read, ch, start_time,
flow_cell_id, sample_id, protocol_group_id) are carried into the BLOW5 read
group and the ONT auxiliary fields (channel_number, read_number, start_time
in samples since exp_start_time, median_before, start_mux). Reads without
them get a run ID, flow cell and channel placement derived from --seed
(nanopore.ChannelModel). Each run is written as its own series of files,
<flowcell>_<runid>_<n>.blow5, of --reads-per-file reads.

Outputs in --out: the BLOW5 files, truth.fa and phantom-manifest.json
(roles "signal" and "truth"; seed, profile, inputs by name and SHA-256).
Output is a pure function of the inputs, flags and --seed: byte-identical
for any --workers.

Example (from the repository root; the last step needs Dorado, blue-crab,
minimap2 and samtools, see scripts/acceptance/README.md):
  phantom-science generate --platform ont --config configs/sequencers/ont_promethion24.yaml \
      --total-bases 5M --seed 1 --output reads
  model=$(scripts/fetch-pore-model.sh dna_r10.4.1_e8.2_400bps)
  phantom-science signal --input reads/*.fastq.gz --out sig --pore-model "$model" --seed 1
  scripts/acceptance/dorado_identity.sh --blow5 sig/*_0.blow5 --truth sig/truth.fa`,
	RunE: func(cmd *cobra.Command, args []string) error {
		return runSignal(cmd, cmd.OutOrStdout(), cmd.ErrOrStderr())
	},
}

func init() {
	rootCmd.AddCommand(signalCmd)
	signalCmd.SilenceUsage = true
	f := signalCmd.Flags()
	f.StringArray("input", nil, "Input reads, FASTQ or FASTA, plain or gzip (repeatable; required)")
	f.String("out", "", "Output directory (required)")
	f.String("pore-model", "", "k-mer pore model table (required; see scripts/fetch-pore-model.sh)")
	f.String("model-units", "auto", "Pore model units: auto, pa or standardised")
	f.String("profile", signal.PresetDNAR10, "Signal profile: "+strings.Join(signal.PresetNames(), ", "))
	f.StringArray("profile-set", nil, "Override a profile parameter, NAME=VALUE (repeatable): "+strings.Join(signal.ParameterNames(), ", "))
	f.String("truth", "", "Optional FASTA/FASTQ of true sequences by read ID; matching reads are synthesised from it")
	f.Uint64("seed", 0, "Random seed (0 or unset: draw one and log it to stderr)")
	f.Int("reads-per-file", slow5.DefaultReadsPerFile, "Reads per BLOW5 file")
	f.String("compression", "zlib", "BLOW5 record compression: zlib, zstd or none (PyPI blue-crab cannot read zstd)")
	f.String("signal-compression", "svb-zd", "BLOW5 signal compression: svb-zd or none")
	f.Int("workers", 0, "Synthesis goroutines (0: number of CPUs); output does not depend on it")
	f.Bool("overwrite", false, "Replace existing outputs in --out")
}

// signalOptions are the parsed flags of the signal command.
type signalOptions struct {
	inputs      []string
	outDir      string
	modelPath   string
	modelUnits  string
	profileName string
	profileSet  map[string]float64
	truthPath   string
	seed        uint64
	perFile     int
	rc          slow5.RecordCompression
	sc          slow5.SignalCompression
	workers     int
	overwrite   bool
}

func runSignal(cmd *cobra.Command, out, logw io.Writer) error {
	f := cmd.Flags()
	var o signalOptions
	o.inputs, _ = f.GetStringArray("input")
	o.outDir, _ = f.GetString("out")
	o.modelPath, _ = f.GetString("pore-model")
	o.modelUnits, _ = f.GetString("model-units")
	o.profileName, _ = f.GetString("profile")
	sets, _ := f.GetStringArray("profile-set")
	for _, kv := range sets {
		name, val, ok := strings.Cut(kv, "=")
		v, err := strconv.ParseFloat(strings.TrimSpace(val), 64)
		if !ok || err != nil {
			return fmt.Errorf("--profile-set %q: want NAME=VALUE with a number", kv)
		}
		if o.profileSet == nil {
			o.profileSet = map[string]float64{}
		}
		o.profileSet[strings.TrimSpace(name)] = v
	}
	o.truthPath, _ = f.GetString("truth")
	o.seed, _ = f.GetUint64("seed")
	o.perFile, _ = f.GetInt("reads-per-file")
	o.workers, _ = f.GetInt("workers")
	o.overwrite, _ = f.GetBool("overwrite")
	rc, _ := f.GetString("compression")
	sc, _ := f.GetString("signal-compression")

	switch {
	case len(o.inputs) == 0:
		return fmt.Errorf("--input is required")
	case o.outDir == "":
		return fmt.Errorf("--out is required")
	case o.modelPath == "":
		return fmt.Errorf("--pore-model is required (fetch one with scripts/fetch-pore-model.sh)")
	case o.perFile <= 0:
		return fmt.Errorf("--reads-per-file must be positive")
	case o.workers < 0:
		return fmt.Errorf("--workers must not be negative")
	}
	switch rc {
	case "zlib":
		o.rc = slow5.RecordZlib
	case "zstd":
		o.rc = slow5.RecordZstd
	case "none":
		o.rc = slow5.RecordNone
	default:
		return fmt.Errorf("unknown --compression %q: want zlib, zstd or none", rc)
	}
	switch sc {
	case "svb-zd":
		o.sc = slow5.SignalSvbZd
	case "none":
		o.sc = slow5.SignalNone
	default:
		return fmt.Errorf("unknown --signal-compression %q: want svb-zd or none", sc)
	}
	switch o.modelUnits {
	case "auto", "pa", "standardised":
	default:
		return fmt.Errorf("unknown --model-units %q: want auto, pa or standardised", o.modelUnits)
	}
	if o.workers == 0 {
		o.workers = runtime.NumCPU()
	}
	seed, drawn, err := genomics.ResolveSeed(o.seed)
	if err != nil {
		return err
	}
	if drawn {
		fmt.Fprintf(logw, "   [info] seed %d (pass --seed %d to reproduce this run)\n", seed, seed)
	}
	o.seed = seed

	res, err := signalRun(o, logw)
	if err != nil {
		return err
	}
	fmt.Fprintf(out, "Wrote %d reads (%d samples, %.1f s of signal) in %d BLOW5 file(s) for %d run(s) to %s\n",
		res.reads, res.samples, res.seconds, res.files, res.runs, filepath.Clean(o.outDir))
	if res.fromTruth > 0 {
		fmt.Fprintf(out, "  %d reads synthesised from --truth sequences, %d from the input sequences\n", res.fromTruth, res.reads-res.fromTruth)
	}
	if res.remapped > 0 {
		fmt.Fprintf(out, "  %d non-UUID read IDs mapped to UUIDs (original IDs in %s)\n", res.remapped, SignalTruthName)
	}
	if res.skipped > 0 {
		fmt.Fprintf(logw, "   [warn] skipped %d reads shorter than the pore model k\n", res.skipped)
	}
	return nil
}

// signalResult summarises a signal run.
type signalResult struct {
	reads, samples, skipped int64
	fromTruth, remapped     int64
	files, runs             int
	seconds                 float64
}

// sigRun is one sequencing run (read group) of the output: the reads of one
// MinKNOW run ID, or the run synthesised from the seed for reads without
// one. Each run has its own series of BLOW5 files.
type sigRun struct {
	key      string // run ID from the headers, "" for the synthesised run
	index    int    // 1-based, in order of first appearance
	info     nanopore.RunInfo
	minStart time.Time // earliest header start_time (zero if none)
	expStart time.Time // exp_start_time: start_time origin
	place    int       // reads needing a synthesised placement
	maxLen   int
	channels *nanopore.ChannelModel
	chanRand *rand.Rand
	w        *slow5.ChunkedWriter
	fcID     string // flow cell and run ID as used in file names
	runName  string
	reads    int64
	bases    int64
	// Per-file read and base counts, in file order.
	fileReads, fileBases []int64
}

// ontFields are the MinKNOW header fields of one read.
type ontFields struct {
	runID, flowCell, sampleID, protocolGroup string
	read, ch                                 int
	start                                    time.Time
	hasRead, hasCh, hasStart                 bool
}

func (f ontFields) placed() bool { return f.hasRead && f.hasCh && f.hasStart }

// parseONTComment extracts MinKNOW key=value fields from a FASTQ comment.
func parseONTComment(c string) ontFields {
	var f ontFields
	for _, tok := range strings.Fields(c) {
		k, v, ok := strings.Cut(tok, "=")
		if !ok || v == "" {
			continue
		}
		switch k {
		case "runid":
			f.runID = v
		case "flow_cell_id":
			f.flowCell = v
		case "sample_id":
			f.sampleID = v
		case "protocol_group_id":
			f.protocolGroup = v
		case "read":
			if n, err := strconv.Atoi(v); err == nil && n >= 0 && n <= math.MaxInt32 {
				f.read, f.hasRead = n, true
			}
		case "ch":
			if n, err := strconv.Atoi(v); err == nil && n > 0 {
				f.ch, f.hasCh = n, true
			}
		case "start_time":
			if t, err := time.Parse(time.RFC3339Nano, v); err == nil {
				f.start, f.hasStart = t.UTC(), true
			}
		}
	}
	return f
}

// signalUUIDNamespace is the UUIDv5 namespace of mapped read IDs: the
// RFC 4122 URL namespace applied to the project URL.
var signalUUIDNamespace = uuidV5(uuidNamespaceURL, "https://github.com/scttfrdmn/phantom-science/signal")

// uuidNamespaceURL is the RFC 4122 URL namespace, 6ba7b811-9dad-11d1-80b4-00c04fd430c8.
var uuidNamespaceURL = [16]byte{0x6b, 0xa7, 0xb8, 0x11, 0x9d, 0xad, 0x11, 0xd1, 0x80, 0xb4, 0x00, 0xc0, 0x4f, 0xd4, 0x30, 0xc8}

// uuidV5 is the RFC 4122 name-based (SHA-1) UUID of name in namespace ns.
func uuidV5(ns [16]byte, name string) [16]byte {
	h := sha1.New() // #nosec G401 -- RFC 4122 version 5 is defined over SHA-1
	h.Write(ns[:])
	h.Write([]byte(name))
	var u [16]byte
	copy(u[:], h.Sum(nil))
	u[6] = u[6]&0x0f | 0x50 // version 5
	u[8] = u[8]&0x3f | 0x80 // RFC 4122 variant
	return u
}

func formatUUID(u [16]byte) string {
	var b [36]byte
	hex.Encode(b[0:8], u[0:4])
	b[8] = '-'
	hex.Encode(b[9:13], u[4:6])
	b[13] = '-'
	hex.Encode(b[14:18], u[6:8])
	b[18] = '-'
	hex.Encode(b[19:23], u[8:10])
	b[23] = '-'
	hex.Encode(b[24:36], u[10:16])
	return string(b[:])
}

// isUUID reports whether s has the 8-4-4-4-12 hex form of a UUID.
func isUUID(s string) bool {
	if len(s) != 36 {
		return false
	}
	for i := 0; i < 36; i++ {
		c := s[i]
		if i == 8 || i == 13 || i == 18 || i == 23 {
			if c != '-' {
				return false
			}
			continue
		}
		if !(c >= '0' && c <= '9' || c >= 'a' && c <= 'f' || c >= 'A' && c <= 'F') {
			return false
		}
	}
	return true
}

// signalReadID returns the BLOW5 read ID for an input read ID: the ID itself
// (lowercased) if it is a UUID, otherwise its UUIDv5. mapped reports whether
// the result differs from name.
func signalReadID(name string) (id string, mapped bool) {
	if isUUID(name) {
		id = strings.ToLower(name)
		return id, id != name
	}
	return formatUUID(uuidV5(signalUUIDNamespace, name)), true
}

// seqFile reads FASTQ (4-line records) or FASTA (multi-line) records from a
// plain or gzip file, hashing the file's bytes as they are read.
type seqFile struct {
	path  string
	f     *os.File
	raw   io.Reader // file bytes (through the hash)
	br    *bufio.Reader
	gz    *gzip.Reader
	h     hash.Hash
	n     int64
	fasta bool
	began bool
	line  int
}

type countWriter struct{ n *int64 }

func (c countWriter) Write(p []byte) (int, error) { *c.n += int64(len(p)); return len(p), nil }

func openSeqFile(path string) (*seqFile, error) {
	f, err := os.Open(path) // #nosec G304 -- user-supplied input path
	if err != nil {
		return nil, err
	}
	s := &seqFile{path: path, f: f, h: sha256.New()}
	rawBuf := bufio.NewReaderSize(io.TeeReader(f, io.MultiWriter(s.h, countWriter{&s.n})), 1<<20)
	s.raw = rawBuf
	if magic, _ := rawBuf.Peek(2); len(magic) == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
		gz, err := gzip.NewReader(rawBuf)
		if err != nil {
			_ = f.Close()
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		s.gz = gz
		s.br = bufio.NewReaderSize(gz, 1<<20)
	} else {
		s.br = rawBuf
	}
	return s, nil
}

// finish reads the rest of the file (so the hash covers all of it) and
// closes it, returning the size and hex SHA-256.
func (s *seqFile) finish() (manifest.Input, error) {
	_, err := io.Copy(io.Discard, s.raw)
	if s.gz != nil {
		err = errors.Join(err, s.gz.Close())
	}
	err = errors.Join(err, s.f.Close())
	return manifest.Input{Name: filepath.Base(s.path), Bytes: s.n, SHA256: hex.EncodeToString(s.h.Sum(nil))}, err
}

func (s *seqFile) close() { _ = s.f.Close() }

func (s *seqFile) readLine() ([]byte, error) {
	b, err := s.br.ReadBytes('\n')
	if err == io.EOF && len(b) > 0 {
		err = nil
	}
	if err != nil {
		return nil, err
	}
	s.line++
	return bytes.TrimRight(b, "\r\n"), nil
}

func (s *seqFile) errorf(format string, args ...any) error {
	return fmt.Errorf("%s:%d: %s", s.path, s.line, fmt.Sprintf(format, args...))
}

// next returns the next record's ID, header comment and sequence, or io.EOF.
func (s *seqFile) next() (name, comment string, seq []byte, err error) {
	var head []byte
	for {
		if head, err = s.readLine(); err != nil {
			return "", "", nil, err
		}
		if len(head) > 0 {
			break
		}
	}
	if !s.began {
		s.began = true
		switch head[0] {
		case '@':
		case '>':
			s.fasta = true
		default:
			return "", "", nil, s.errorf("not FASTQ or FASTA (line starts with %q)", head[0])
		}
	}
	marker := byte('@')
	if s.fasta {
		marker = '>'
	}
	if head[0] != marker {
		return "", "", nil, s.errorf("expected a %q header line", marker)
	}
	name, comment, _ = strings.Cut(string(head[1:]), " ")
	if i := strings.IndexByte(name, '\t'); i >= 0 {
		name, comment = name[:i], name[i+1:]+" "+comment
	}
	if name == "" {
		return "", "", nil, s.errorf("empty read ID")
	}
	if s.fasta {
		for {
			p, perr := s.br.Peek(1)
			if perr != nil || p[0] == '>' {
				break
			}
			l, lerr := s.readLine()
			if lerr != nil {
				return "", "", nil, lerr
			}
			seq = append(seq, l...)
		}
		return name, comment, seq, nil
	}
	if seq, err = s.readLine(); err != nil {
		return "", "", nil, s.errorf("truncated FASTQ record %s", name)
	}
	plus, err := s.readLine()
	if err != nil || len(plus) == 0 || plus[0] != '+' {
		return "", "", nil, s.errorf("FASTQ record %s: missing '+' line", name)
	}
	qual, err := s.readLine()
	if err != nil || len(qual) != len(seq) {
		return "", "", nil, s.errorf("FASTQ record %s: quality length differs from sequence length", name)
	}
	return name, comment, append([]byte(nil), seq...), nil
}

// loadTruth reads a FASTA/FASTQ of true sequences keyed by record ID.
func loadTruth(path string) (map[string][]byte, manifest.Input, error) {
	s, err := openSeqFile(path)
	if err != nil {
		return nil, manifest.Input{}, err
	}
	defer s.close()
	truth := map[string][]byte{}
	for {
		name, _, seq, err := s.next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, manifest.Input{}, err
		}
		if _, dup := truth[name]; dup {
			return nil, manifest.Input{}, fmt.Errorf("%s: duplicate record %s", path, name)
		}
		truth[name] = seq
	}
	in, err := s.finish()
	return truth, in, err
}

// detectModelUnits returns "pa" for a squigulator/nanopolish table (header
// row or level_stdv column) and "standardised" for a two-column ONT
// kmer_models table.
func detectModelUnits(path string) (string, error) {
	f, err := os.Open(path) // #nosec G304 -- user-supplied model path
	if err != nil {
		return "", fmt.Errorf("pore model: %w", err)
	}
	defer f.Close()
	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return "", fmt.Errorf("pore model %s: %w", path, err)
		}
		defer gz.Close()
		r = gz
	}
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		fields := strings.Fields(sc.Text())
		if len(fields) == 0 || strings.HasPrefix(fields[0], "#") {
			continue
		}
		if strings.EqualFold(fields[0], "kmer") || len(fields) >= 3 {
			return "pa", nil
		}
		return "standardised", nil
	}
	if err := sc.Err(); err != nil {
		return "", fmt.Errorf("pore model %s: %w", path, err)
	}
	return "", fmt.Errorf("pore model %s: no k-mer rows", path)
}

// Run metadata defaults derived per profile: MinKNOW read-group attributes of
// a PromethION R10.4.1 (Kit 14) or RNA004 run.
func profileAttrs(p signal.Profile) (experiment, kit, product string) {
	if p.Reverse {
		return "rna", "sqk-rna004", "FLO-PRO004RA"
	}
	return "genomic_dna", "sqk-lsk114", "FLO-PRO114M"
}

// fileNamePart makes s usable in a chunk file name (letters, digits, '-'
// and '.'; see slow5.ChunkFileName).
func fileNamePart(s string) string {
	b := []byte(s)
	for i, c := range b {
		if !(c == '-' || c == '.' || c >= '0' && c <= '9' || c >= 'a' && c <= 'z' || c >= 'A' && c <= 'Z') {
			b[i] = '-'
		}
	}
	if len(b) == 0 || b[0] == '.' {
		b = append([]byte("x"), b...)
	}
	return string(b)
}

func fnv64(s string) uint64 {
	h := fnv.New64a()
	_, _ = h.Write([]byte(s))
	return h.Sum64()
}

// Generator streams of the signal command, independent of the per-read
// signal streams (signal.ReadRand).
const (
	sigChannelStream = 0x7369675f6368   // "sig_ch": synthesised placements
	sigMuxSalt       = 0x7369675f6d7578 // "sig_mux": start_mux draws
)

// header returns the BLOW5 header of a run.
func (r *sigRun) header(p signal.Profile) slow5.Header {
	experiment, kit, product := profileAttrs(p)
	return slow5.Header{
		ReadGroups: []slow5.ReadGroup{{
			"run_id":                 r.info.RunID,
			"flow_cell_id":           r.info.FlowCellID,
			"flow_cell_product_code": product,
			"asic_id":                strings.ToUpper(strconv.FormatUint(fnv64("asic/"+r.info.RunID), 16)),
			"device_type":            "promethion",
			"exp_start_time":         r.expStart.Format(time.RFC3339),
			"experiment_type":        experiment,
			"sample_frequency":       strconv.FormatFloat(p.SamplingRate, 'f', -1, 64),
			"sample_id":              r.info.SampleID,
			"sequencing_kit":         kit,
			"protocol_group_id":      r.info.ProtocolGroupID,
			"protocol_run_id":        formatUUID(uuidV5(signalUUIDNamespace, "protocol_run/"+r.info.RunID)),
			"exp_script_purpose":     "sequencing_run",
			"signal_profile":         p.Name,
		}},
		Aux: slow5.ONTAuxFields,
	}
}

// sigJob is one read on its way through synthesis.
type sigJob struct {
	id     string
	source string // input read ID when it differs from id
	seq    []byte
	run    *sigRun
	aux    slow5.ONTAux
	done   chan sigDone
}

type sigDone struct {
	rec *slow5.Record
	err error
}

// sigPlan is the result of the first pass over the inputs.
type sigPlan struct {
	runs    []*sigRun
	byKey   map[string]*sigRun
	inputs  []manifest.Input
	reads   int64
	skipped int64
}

// sigSource resolves an input record to its output read ID, the sequence to
// synthesise and whether that came from --truth.
func sigSource(name string, seq []byte, truth map[string][]byte) (id string, mapped bool, src []byte, fromTruth bool) {
	id, mapped = signalReadID(name)
	if t, ok := truth[name]; ok {
		return id, mapped, t, true
	}
	if t, ok := truth[id]; ok {
		return id, mapped, t, true
	}
	return id, mapped, seq, false
}

// scanInputs is the first pass: it checks read IDs are unique, groups reads
// by run, finds each run's earliest start time and hashes the inputs.
func scanInputs(o signalOptions, truth map[string][]byte, k int) (*sigPlan, error) {
	pl := &sigPlan{byKey: map[string]*sigRun{}}
	seen := map[[16]byte]struct{}{} // 16 bytes per read, not the ID strings
	for _, path := range o.inputs {
		s, err := openSeqFile(path)
		if err != nil {
			return nil, err
		}
		for {
			name, comment, seq, err := s.next()
			if err == io.EOF {
				break
			}
			if err != nil {
				s.close()
				return nil, err
			}
			id, _, src, _ := sigSource(name, seq, truth)
			var key [16]byte
			_, _ = hex.Decode(key[:], []byte(strings.ReplaceAll(id, "-", "")))
			if _, dup := seen[key]; dup {
				s.close()
				return nil, fmt.Errorf("%s: read %s maps to read ID %s, which an earlier read has; read IDs must be unique", path, name, id)
			}
			seen[key] = struct{}{}
			if len(src) < k {
				pl.skipped++
				continue
			}
			f := parseONTComment(comment)
			r := pl.byKey[f.runID]
			if r == nil {
				r = &sigRun{key: f.runID, index: len(pl.runs) + 1}
				r.info.RunID = f.runID
				pl.byKey[f.runID] = r
				pl.runs = append(pl.runs, r)
			}
			if r.info.FlowCellID == "" {
				r.info.FlowCellID = f.flowCell
			}
			if r.info.SampleID == "" {
				r.info.SampleID = f.sampleID
			}
			if r.info.ProtocolGroupID == "" {
				r.info.ProtocolGroupID = f.protocolGroup
			}
			if f.placed() {
				if r.minStart.IsZero() || f.start.Before(r.minStart) {
					r.minStart = f.start
				}
			} else {
				r.place++
				r.maxLen = max(r.maxLen, len(src))
			}
			pl.reads++
		}
		in, err := s.finish()
		if err != nil {
			return nil, err
		}
		pl.inputs = append(pl.inputs, in)
	}
	if pl.reads == 0 {
		return nil, fmt.Errorf("no reads of at least %d bases in the input", k)
	}
	return pl, nil
}

// setupRuns fills in each run's metadata and opens its writer.
func setupRuns(o signalOptions, pl *sigPlan, p signal.Profile) error {
	names := map[string]bool{}
	for _, r := range pl.runs {
		runSeed := o.seed
		if r.key != "" {
			runSeed ^= fnv64("run/" + r.key)
		}
		r.info.SpeedBPS = p.SpeedBPS
		r.info = r.info.WithDefaults(runSeed)
		if err := r.info.Validate(); err != nil {
			return fmt.Errorf("run %s: %w", r.info.RunID, err)
		}
		r.expStart = r.info.Start
		if !r.minStart.IsZero() {
			r.expStart = r.minStart.Truncate(time.Second)
		}
		if r.place > 0 {
			r.channels = nanopore.NewChannelModel(r.info, 0, r.place, r.maxLen)
			// #nosec G404 -- simulation metadata, reproducible by design
			r.chanRand = rand.New(rand.NewPCG(runSeed, sigChannelStream))
		}
		r.fcID, r.runName = fileNamePart(r.info.FlowCellID), fileNamePart(r.info.RunID)
		prefix := r.fcID + "_" + r.runName
		if names[prefix] {
			return fmt.Errorf("runs %s and another share the file prefix %s", r.info.RunID, prefix)
		}
		names[prefix] = true
		first := filepath.Join(o.outDir, slow5.ChunkFileName(r.fcID, r.runName, 0))
		if _, err := os.Stat(first); err == nil && !o.overwrite {
			return fmt.Errorf("%s exists (use --overwrite)", first)
		}
		w, err := slow5.NewChunkedWriter(o.outDir, r.header(p), slow5.ChunkOptions{
			ReadsPerFile: o.perFile,
			BLOW5:        slow5.BLOW5Options{RecordCompression: o.rc, SignalCompression: o.sc},
			FlowCellID:   r.fcID, RunID: r.runName,
			Overwrite: o.overwrite,
		})
		if err != nil {
			return err
		}
		r.w = w
	}
	return nil
}

// signalRun synthesises every read of the inputs and writes the outputs.
func signalRun(o signalOptions, logw io.Writer) (*signalResult, error) {
	p, err := signal.Preset(o.profileName)
	if err != nil {
		return nil, err
	}
	for name, v := range o.profileSet {
		if err := p.Set(name, v); err != nil {
			return nil, err
		}
	}
	if err := p.Validate(); err != nil {
		return nil, err
	}
	units := o.modelUnits
	if units == "auto" {
		if units, err = detectModelUnits(o.modelPath); err != nil {
			return nil, err
		}
	}
	var lopts []signal.LoadOption
	if units == "standardised" {
		lopts = append(lopts, signal.WithScaling(p.ModelScale, p.ModelShift), signal.WithDefaultStdv(p.ModelStdv))
	}
	model, err := signal.LoadPoreModel(o.modelPath, lopts...)
	if err != nil {
		return nil, err
	}
	if p.ModelK > 0 && model.K != p.ModelK {
		fmt.Fprintf(logw, "   [warn] pore model k=%d, profile %s expects k=%d\n", model.K, p.Name, p.ModelK)
	}
	_, modelSum, err := manifest.HashFile(o.modelPath)
	if err != nil {
		return nil, err
	}
	modelSize, _ := fileSize(o.modelPath)

	var truth map[string][]byte
	var truthIn *manifest.Input
	if o.truthPath != "" {
		t, in, err := loadTruth(o.truthPath)
		if err != nil {
			return nil, err
		}
		truth, truthIn = t, &in
	}

	truthPath := filepath.Join(o.outDir, SignalTruthName)
	manPath := filepath.Join(o.outDir, manifest.FileName)
	if !o.overwrite {
		for _, f := range []string{truthPath, manPath} {
			if _, err := os.Stat(f); err == nil {
				return nil, fmt.Errorf("%s exists (use --overwrite)", f)
			}
		}
	}

	pl, err := scanInputs(o, truth, model.K)
	if err != nil {
		return nil, err
	}
	if err := setupRuns(o, pl, p); err != nil {
		return nil, err
	}
	tf, err := os.OpenFile(filepath.Clean(truthPath), os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o600)
	if err != nil {
		return nil, err
	}
	tw := bufio.NewWriterSize(tf, 1<<20)

	res := &signalResult{runs: len(pl.runs), skipped: pl.skipped}
	synthErr := synthesizeAll(o, pl, truth, model, p, func(j *sigJob, rec *slow5.Record) error {
		if err := j.run.w.Write(rec); err != nil {
			return err
		}
		r := j.run
		if r.reads%int64(o.perFile) == 0 {
			r.fileReads, r.fileBases = append(r.fileReads, 0), append(r.fileBases, 0)
		}
		r.reads++
		r.bases += int64(len(j.seq))
		r.fileReads[len(r.fileReads)-1]++
		r.fileBases[len(r.fileBases)-1] += int64(len(j.seq))
		res.reads++
		res.samples += int64(len(rec.RawSignal))
		if j.source != "" {
			res.remapped++
		}
		tw.WriteByte('>')
		tw.WriteString(j.id)
		if j.source != "" {
			tw.WriteString(" source=")
			tw.WriteString(j.source)
		}
		tw.WriteByte('\n')
		tw.Write(j.seq)
		_, err := tw.WriteString("\n")
		return err
	}, &res.fromTruth)
	closeErr := errors.Join(tw.Flush(), tf.Close())
	for _, r := range pl.runs {
		closeErr = errors.Join(closeErr, r.w.Close())
	}
	if err := errors.Join(synthErr, closeErr); err != nil {
		return nil, err
	}
	res.seconds = float64(res.samples) / p.SamplingRate

	m := manifest.New(Version, o.seed, pl.runs[0].expStart, types.SequencerConfig{})
	m.Signal = &manifest.Signal{
		Profile:           p.Name,
		ProfileOverrides:  o.profileSet,
		PoreModel:         manifest.Input{Name: filepath.Base(o.modelPath), Bytes: modelSize, SHA256: modelSum},
		ModelUnits:        units,
		RecordCompression: o.rc.String(),
		SignalCompression: o.sc.String(),
		ReadsPerFile:      o.perFile,
		Inputs:            pl.inputs,
		Truth:             truthIn,
		Samples:           res.samples,
		Skipped:           pl.skipped,
	}
	for _, r := range pl.runs {
		files := r.w.Files()
		if o.overwrite {
			if err := removeStaleChunks(o.outDir, r, len(files)); err != nil {
				return nil, err
			}
		}
		m.FlowCells = append(m.FlowCells, manifest.FlowCell{
			Index: r.index, Seed: o.seed, RunID: r.info.RunID, FlowCellID: r.info.FlowCellID,
			Files: len(files), Reads: r.reads, Bases: r.bases,
		})
		for i, path := range files {
			n, sum, err := manifest.HashFile(path)
			if err != nil {
				return nil, err
			}
			m.Add(manifest.Output{
				Path: filepath.Base(path), Role: manifest.RoleSignal, FlowCell: r.index, Bytes: n, SHA256: sum,
				Reads: r.fileReads[i], Bases: r.fileBases[i],
			})
		}
		res.files += len(files)
	}
	n, sum, err := manifest.HashFile(truthPath)
	if err != nil {
		return nil, err
	}
	m.Add(manifest.Output{Path: SignalTruthName, Role: manifest.RoleTruth, Bytes: n, SHA256: sum})
	data, err := m.Marshal()
	if err != nil {
		return nil, err
	}
	if err := os.WriteFile(manPath, data, 0o600); err != nil {
		return nil, fmt.Errorf("writing manifest: %w", err)
	}
	return res, nil
}

func fileSize(path string) (int64, error) {
	st, err := os.Stat(path)
	if err != nil {
		return 0, err
	}
	return st.Size(), nil
}

// removeStaleChunks deletes chunk files of run r numbered from n up, left by
// an earlier, larger run into the same directory, so the directory holds
// exactly this run's files.
func removeStaleChunks(dir string, r *sigRun, n int) error {
	for i := n; ; i++ {
		path := filepath.Join(dir, slow5.ChunkFileName(r.fcID, r.runName, i))
		if err := os.Remove(path); err != nil {
			if errors.Is(err, os.ErrNotExist) {
				return nil
			}
			return err
		}
	}
}

// synthesizeAll is the second pass. A reader goroutine assigns each read its
// ID, run and placement in input order (placements are drawn sequentially
// from each run's channel model); o.workers goroutines synthesise signal,
// which depends only on (seed, read ID); commit is called for every read in
// input order. The output is therefore the same for any worker count.
func synthesizeAll(o signalOptions, pl *sigPlan, truth map[string][]byte, model *signal.PoreModel, p signal.Profile,
	commit func(*sigJob, *slow5.Record) error, fromTruth *int64) error {
	order := make(chan *sigJob, 2*o.workers)
	work := make(chan *sigJob, o.workers)
	quit := make(chan struct{})
	var wg sync.WaitGroup
	var readErr error

	wg.Add(1)
	go func() {
		defer wg.Done()
		defer close(work)
		defer close(order)
		readErr = produceJobs(o, pl, truth, model.K, p, quit, order, work, fromTruth)
	}()
	for i := 0; i < o.workers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range work {
				select {
				case <-quit:
					j.done <- sigDone{err: errors.New("cancelled")}
					continue
				default:
				}
				rd, err := signal.SynthesizeRead(o.seed, j.id, j.seq, p, model)
				if err != nil {
					j.done <- sigDone{err: fmt.Errorf("read %s: %w", j.id, err)}
					continue
				}
				a := j.aux
				a.MedianBefore = rd.MedianBefore
				j.done <- sigDone{rec: &slow5.Record{
					ReadID: rd.ID, Digitisation: rd.Digitisation, Offset: rd.Offset, Range: rd.Range,
					SamplingRate: rd.SamplingRate, RawSignal: rd.Raw, Aux: a.Map(),
				}}
			}
		}()
	}

	var err error
	for j := range order {
		if err != nil {
			continue // drain
		}
		d := <-j.done
		if d.err != nil {
			err = d.err
		} else {
			err = commit(j, d.rec)
		}
		if err != nil {
			close(quit)
		}
	}
	wg.Wait()
	return errors.Join(err, readErr)
}

// produceJobs reads the inputs again and queues every read for synthesis.
func produceJobs(o signalOptions, pl *sigPlan, truth map[string][]byte, k int, p signal.Profile,
	quit <-chan struct{}, order, work chan<- *sigJob, fromTruth *int64) error {
	for _, path := range o.inputs {
		s, err := openSeqFile(path)
		if err != nil {
			return err
		}
		for {
			name, comment, seq, err := s.next()
			if err == io.EOF {
				break
			}
			if err != nil {
				s.close()
				return err
			}
			id, mapped, src, ft := sigSource(name, seq, truth)
			if len(src) < k {
				continue
			}
			if ft {
				*fromTruth++
			}
			f := parseONTComment(comment)
			r := pl.byKey[f.runID]
			j := &sigJob{id: id, seq: src, run: r, done: make(chan sigDone, 1)}
			if mapped {
				j.source = name
			}
			if f.placed() {
				j.aux.ChannelNumber = strconv.Itoa(f.ch)
				j.aux.ReadNumber = int32(f.read) // #nosec G115 -- range checked in parseONTComment
				j.aux.StartTime = uint64(math.Round(f.start.Sub(r.expStart).Seconds() * p.SamplingRate))
			} else {
				pc := r.channels.Next(r.chanRand, len(src))
				j.aux.ChannelNumber = strconv.Itoa(pc.Channel)
				j.aux.ReadNumber = int32(min(pc.Read, math.MaxInt32)) // #nosec G115 -- clamped
				j.aux.StartTime = uint64(math.Round(pc.Start.Seconds() * p.SamplingRate))
			}
			// MinKNOW cycles wells (mux 1-4) over a run; draw the starting
			// well from a stream of its own so it never perturbs the signal.
			j.aux.StartMux = uint8(1 + signal.ReadRand(o.seed^sigMuxSalt, id).IntN(4)) // #nosec G115 -- 1..4
			select {
			case order <- j:
			case <-quit:
				s.close()
				return nil
			}
			select {
			case work <- j:
			case <-quit:
				// The committer has stopped; it drains order without
				// waiting on this job.
				s.close()
				return nil
			}
		}
		s.close()
	}
	return nil
}
