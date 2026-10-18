//go:build !release

package cmd

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/signal"
	"github.com/scttfrdmn/phantom-science/internal/genomics/slow5"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"github.com/spf13/pflag"
)

var signalFlagNames = []string{
	"input", "out", "pore-model", "model-units", "profile", "truth", "seed", "reads-per-file",
	"compression", "signal-compression", "workers", "overwrite",
}

// runSignalCLI runs "phantom-science signal args..." through the root
// command and resets the subcommand's flags afterwards.
func runSignalCLI(t *testing.T, args ...string) (stdout, stderr string, err error) {
	t.Helper()
	var out, errb bytes.Buffer
	rootCmd.SetOut(&out)
	rootCmd.SetErr(&errb)
	rootCmd.SetArgs(append([]string{"signal"}, args...))
	defer func() {
		rootCmd.SetArgs([]string{})
		rootCmd.SetOut(nil)
		rootCmd.SetErr(nil)
		for _, name := range signalFlagNames {
			f := signalCmd.Flags().Lookup(name)
			if sv, ok := f.Value.(pflag.SliceValue); ok {
				_ = sv.Replace(nil)
			} else {
				_ = f.Value.Set(f.DefValue)
			}
			f.Changed = false
		}
	}()
	err = rootCmd.Execute()
	return out.String(), errb.String(), err
}

// writeTestPoreModel writes a synthetic k=5 table in pA (squigulator format:
// header row, level_mean and level_stdv columns): levels 45-125 pA, stdv
// 1.0-2.2 pA. ONT's real tables are never used in tests.
func writeTestPoreModel(t *testing.T) string {
	t.Helper()
	var b strings.Builder
	b.WriteString("#k\t5\nkmer\tlevel_mean\tlevel_stdv\n")
	for i := 0; i < 1<<10; i++ {
		level := 4*math.Mod(float64(i)*0.6180339887, 1) - 2
		fmt.Fprintf(&b, "%s\t%.6f\t%.6f\n", signal.KmerString(i, 5), 85+20*level, 1+0.3*float64(i%5))
	}
	p := filepath.Join(t.TempDir(), "test_5mer.tsv")
	if err := os.WriteFile(p, []byte(b.String()), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

// testModelAndProfile loads the test model as the command does and returns
// the R10 preset.
func testModelAndProfile(t *testing.T, path string) (*signal.PoreModel, signal.Profile) {
	t.Helper()
	m, err := signal.LoadPoreModel(path)
	if err != nil {
		t.Fatal(err)
	}
	p, err := signal.Preset(signal.PresetDNAR10)
	if err != nil {
		t.Fatal(err)
	}
	return m, p
}

// generateONTReads runs the engine (generate --platform ont) and returns the
// FASTQ files: 2 files of 15 reads of 600 bases with MinKNOW headers.
func generateONTReads(t *testing.T) []string {
	t.Helper()
	dir := t.TempDir()
	e := genomics.NewEngine(types.SequencerConfig{
		Platform: "ont", SampleName: "sigS", Lane: 1, OutputDir: dir,
		ReadLength: 600, ReadsPerBatch: 15, TotalBatches: 2,
		BaseQuality: 40, QualityDecay: 0.01, Seed: 41,
		ONT: types.ONTSpec{Channels: 16, ProtocolGroupID: "pg1"},
	})
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	files := e.GeneratedFiles()
	if len(files) != 2 {
		t.Fatalf("generate wrote %v", files)
	}
	return files
}

type fqRead struct {
	name string
	f    ontFields
	seq  string
}

// readFASTQs parses FASTQ (or FASTA) files with the command's own reader.
func readFASTQs(t *testing.T, paths ...string) []fqRead {
	t.Helper()
	var out []fqRead
	for _, p := range paths {
		s, err := openSeqFile(p)
		if err != nil {
			t.Fatal(err)
		}
		for {
			name, comment, seq, err := s.next()
			if err != nil {
				break
			}
			out = append(out, fqRead{name, parseONTComment(comment), string(seq)})
		}
		if _, err := s.finish(); err != nil {
			t.Fatal(err)
		}
	}
	return out
}

// readTruthFA returns truth.fa as ID -> (sequence, source).
func readTruthFA(t *testing.T, path string) (seqs, sources map[string]string, order []string) {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	seqs, sources = map[string]string{}, map[string]string{}
	lines := strings.Split(strings.TrimSuffix(string(data), "\n"), "\n")
	if len(lines)%2 != 0 {
		t.Fatalf("truth.fa has %d lines, want 2 per read", len(lines))
	}
	for i := 0; i < len(lines); i += 2 {
		id, rest, _ := strings.Cut(strings.TrimPrefix(lines[i], ">"), " ")
		seqs[id] = lines[i+1]
		if src, ok := strings.CutPrefix(rest, "source="); ok {
			sources[id] = src
		}
		order = append(order, id)
	}
	return seqs, sources, order
}

// readBLOW5Dir reads every .blow5 file of dir, in name order.
func readBLOW5Dir(t *testing.T, dir string) (files []string, hdrs []slow5.Header, recs []*slow5.Record) {
	t.Helper()
	files, _ = filepath.Glob(filepath.Join(dir, "*.blow5"))
	sort.Slice(files, func(i, j int) bool { return chunkIndex(files[i]) < chunkIndex(files[j]) })
	for _, f := range files {
		fh, err := os.Open(f)
		if err != nil {
			t.Fatal(err)
		}
		rd, err := slow5.NewReader(fh)
		if err != nil {
			t.Fatalf("%s: %v", filepath.Base(f), err)
		}
		if rd.Format() != slow5.FormatBLOW5 {
			t.Errorf("%s: format %v", filepath.Base(f), rd.Format())
		}
		hdrs = append(hdrs, rd.Header())
		for {
			r, err := rd.Read()
			if err != nil {
				break
			}
			recs = append(recs, r)
		}
		_ = rd.Close()
		_ = fh.Close()
	}
	return files, hdrs, recs
}

func chunkIndex(path string) int {
	base := strings.TrimSuffix(filepath.Base(path), ".blow5")
	var n int
	fmt.Sscanf(base[strings.LastIndexByte(base, '_')+1:], "%d", &n)
	return n
}

func readDirFiles(t *testing.T, dir string) map[string][]byte {
	t.Helper()
	ents, err := os.ReadDir(dir)
	if err != nil {
		t.Fatal(err)
	}
	out := map[string][]byte{}
	for _, e := range ents {
		b, err := os.ReadFile(filepath.Join(dir, e.Name()))
		if err != nil {
			t.Fatal(err)
		}
		out[e.Name()] = b
	}
	return out
}

// TestSignalCommand_GenerateRoundTrip is the end-to-end check: ONT reads
// from the engine -> signal -> BLOW5 read back with the slow5 reader.
func TestSignalCommand_GenerateRoundTrip(t *testing.T) {
	fastqs := generateONTReads(t)
	model := writeTestPoreModel(t)
	pm, prof := testModelAndProfile(t, model)
	in := readFASTQs(t, fastqs...)
	if len(in) != 30 {
		t.Fatalf("generate wrote %d reads, want 30", len(in))
	}

	out := filepath.Join(t.TempDir(), "sig")
	args := []string{"--input", fastqs[0], "--input", fastqs[1], "--out", out, "--pore-model", model,
		"--seed", "9", "--reads-per-file", "12", "--workers", "3"}
	stdout, stderr, err := runSignalCLI(t, args...)
	if err != nil {
		t.Fatalf("%v\n%s", err, stderr)
	}
	if !strings.Contains(stdout, "Wrote 30 reads") || !strings.Contains(stdout, "3 BLOW5 file(s) for 1 run(s)") {
		t.Errorf("stdout %q", stdout)
	}
	if !strings.Contains(stderr, "pore model k=5, profile "+prof.Name+" expects k=9") {
		t.Errorf("no k mismatch warning: %q", stderr)
	}

	files, hdrs, recs := readBLOW5Dir(t, out)
	if len(files) != 3 || len(recs) != 30 {
		t.Fatalf("%d files, %d records; want 3 files (12+12+6), 30 records", len(files), len(recs))
	}
	runID, flowCell := in[0].f.runID, in[0].f.flowCell
	for i, f := range files {
		if want := slow5.ChunkFileName(flowCell, runID, i); filepath.Base(f) != want {
			t.Errorf("file %d is %s, want %s", i, filepath.Base(f), want)
		}
	}
	h := hdrs[0]
	rg := h.ReadGroups[0]
	if len(h.ReadGroups) != 1 || rg["run_id"] != runID || rg["flow_cell_id"] != flowCell ||
		rg["sample_id"] != in[0].f.sampleID || rg["protocol_group_id"] != "pg1" || rg["sample_frequency"] != "5000" {
		t.Errorf("read group %v; want run %s flow cell %s", rg, runID, flowCell)
	}
	expStart, err := time.Parse(time.RFC3339, rg["exp_start_time"])
	if err != nil {
		t.Fatalf("exp_start_time %q: %v", rg["exp_start_time"], err)
	}
	if !isUUID(rg["protocol_run_id"]) {
		t.Errorf("protocol_run_id %q is not a UUID", rg["protocol_run_id"])
	}

	truth, sources, order := readTruthFA(t, filepath.Join(out, SignalTruthName))
	var samples, kmers int
	for i, r := range recs {
		src := in[i]
		// Generated ONT read IDs are UUID4s and are kept, in input order.
		if r.ReadID != src.name || order[i] != src.name {
			t.Fatalf("record %d is %s (truth.fa %s), want input read %s", i, r.ReadID, order[i], src.name)
		}
		if truth[r.ReadID] != src.seq || sources[r.ReadID] != "" {
			t.Errorf("truth.fa %s does not hold the input sequence", r.ReadID)
		}
		a, ok := slow5.ONTAuxOf(r)
		if !ok {
			t.Fatalf("record %s lacks the ONT aux fields: %v", r.ReadID, r.Aux)
		}
		// start_time is in samples since exp_start_time (slow5 spec,
		// Table 5); the header times are microsecond-exact.
		wantStart := math.Round(src.f.start.Sub(expStart).Seconds() * prof.SamplingRate)
		if a.ChannelNumber != fmt.Sprint(src.f.ch) || int(a.ReadNumber) != src.f.read || float64(a.StartTime) != wantStart {
			t.Errorf("%s aux ch=%s read=%d start=%d; header ch=%d read=%d start=%.0f samples",
				r.ReadID, a.ChannelNumber, a.ReadNumber, a.StartTime, src.f.ch, src.f.read, wantStart)
		}
		if a.StartMux < 1 || a.StartMux > 4 {
			t.Errorf("%s start_mux %d outside 1..4", r.ReadID, a.StartMux)
		}
		if d := math.Abs(a.MedianBefore-prof.MedianBeforeMean) / prof.MedianBeforeStd; d > 6 {
			t.Errorf("%s median_before %.1f pA is %.1f SD from the profile mean", r.ReadID, a.MedianBefore, d)
		}
		// The stored record is exactly what the synthesiser produces for
		// (seed, read ID): calibration and every sample survive the
		// BLOW5 round trip.
		want, err := signal.SynthesizeRead(9, r.ReadID, []byte(src.seq), prof, pm)
		if err != nil {
			t.Fatal(err)
		}
		if r.Digitisation != want.Digitisation || r.Offset != want.Offset || r.Range != want.Range ||
			r.SamplingRate != want.SamplingRate || a.MedianBefore != want.MedianBefore {
			t.Errorf("%s calibration %g/%g/%g/%g, want %g/%g/%g/%g", r.ReadID, r.Digitisation, r.Offset, r.Range,
				r.SamplingRate, want.Digitisation, want.Offset, want.Range, want.SamplingRate)
		}
		if len(r.RawSignal) != want.NumSamples || !slicesEqual(r.RawSignal, want.Raw) {
			t.Fatalf("%s: %d samples differ from the synthesised %d", r.ReadID, len(r.RawSignal), want.NumSamples)
		}
		// pA = (raw + offset) * range / digitisation: the read body (after
		// the adapter) lies around the model's 45-125 pA levels.
		var sum float64
		for j := want.SeqStart; j < want.SeqEnd; j++ {
			sum += r.PicoAmps(j)
		}
		if mean := sum / float64(want.SeqEnd-want.SeqStart); mean < 60 || mean > 110 {
			t.Errorf("%s mean read current %.1f pA, want within the model's level range", r.ReadID, mean)
		}
		samples += len(r.RawSignal) - want.AdapterSamples
		kmers += len(src.seq) - pm.K + 1
	}
	// Mean samples per base is sampling rate / speed = 12.5 (plus a few
	// percent from stalls: 0.4% of k-mers at 5-20x dwell).
	if perBase := float64(samples) / float64(kmers); perBase < 12.5 || perBase > 14.5 {
		t.Errorf("%.2f samples per base, want about %.1f", perBase, prof.MeanDwell())
	}

	m, err := manifest.Read(filepath.Join(out, manifest.FileName))
	if err != nil {
		t.Fatal(err)
	}
	if err := m.VerifyDir(out); err != nil {
		t.Errorf("VerifyDir: %v", err)
	}
	roles := map[string]int{}
	var fileReads int64
	for _, o := range m.Outputs {
		roles[o.Role]++
		fileReads += o.Reads
	}
	if roles[manifest.RoleSignal] != 3 || roles[manifest.RoleTruth] != 1 || len(m.Outputs) != 4 || fileReads != 30 {
		t.Errorf("manifest outputs %+v", m.Outputs)
	}
	if m.Seed != 9 || m.Totals.Reads != 30 || m.Totals.Files != 3 || m.Totals.Bases != 30*600 ||
		len(m.FlowCells) != 1 || m.FlowCells[0].RunID != runID || !m.RunStart.Equal(expStart) {
		t.Errorf("manifest seed %d totals %+v flow cells %+v run start %v", m.Seed, m.Totals, m.FlowCells, m.RunStart)
	}
	s := m.Signal
	if s == nil || s.Profile != prof.Name || s.ModelUnits != "pa" || s.RecordCompression != "zlib" ||
		s.SignalCompression != "svb-zd" || s.ReadsPerFile != 12 || len(s.Inputs) != 2 || s.Samples == 0 {
		t.Fatalf("manifest signal section %+v", s)
	}
	for i, f := range fastqs {
		n, sum, _ := manifest.HashFile(f)
		if s.Inputs[i] != (manifest.Input{Name: filepath.Base(f), Bytes: n, SHA256: sum}) {
			t.Errorf("input %d recorded as %+v", i, s.Inputs[i])
		}
	}
	raw, _ := os.ReadFile(filepath.Join(out, manifest.FileName))
	if bytes.Contains(raw, []byte(`"config"`)) || bytes.Contains(raw, []byte(filepath.Dir(model))) {
		t.Error("manifest records a sequencer config or a local path")
	}

	// Determinism: one worker reproduces every byte; another seed changes
	// the signal but not the sequences.
	again := filepath.Join(t.TempDir(), "sig")
	if _, stderr, err := runSignalCLI(t, append(args[:len(args)-1:len(args)-1], "1", "--out", again)...); err != nil {
		t.Fatalf("%v\n%s", err, stderr)
	}
	a, b := readDirFiles(t, out), readDirFiles(t, again)
	if len(a) != 5 || len(a) != len(b) {
		t.Fatalf("outputs %d and %d files, want 5", len(a), len(b))
	}
	for name := range a {
		if !bytes.Equal(a[name], b[name]) {
			t.Errorf("%s differs between --workers 3 and --workers 1", name)
		}
	}
	other := filepath.Join(t.TempDir(), "sig")
	if _, stderr, err := runSignalCLI(t, "--input", fastqs[0], "--input", fastqs[1], "--out", other,
		"--pore-model", model, "--seed", "10", "--reads-per-file", "12"); err != nil {
		t.Fatalf("%v\n%s", err, stderr)
	}
	c := readDirFiles(t, other)
	if !bytes.Equal(a[SignalTruthName], c[SignalTruthName]) {
		t.Error("truth.fa depends on the seed")
	}
	for _, f := range files {
		if bytes.Equal(a[filepath.Base(f)], c[filepath.Base(f)]) {
			t.Errorf("%s is identical for seeds 9 and 10", filepath.Base(f))
		}
	}
}

func slicesEqual(a, b []int16) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func writeGzip(t *testing.T, path, content string) {
	t.Helper()
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	_, _ = zw.Write([]byte(content))
	if err := zw.Close(); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, buf.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
}

// TestSignalCommand_PlainFASTQWithTruth covers reads without MinKNOW
// headers or UUID IDs (e.g. Illumina-style names), gzip input, --truth, a
// read shorter than k, and zstd output.
func TestSignalCommand_PlainFASTQWithTruth(t *testing.T) {
	dir := t.TempDir()
	model := writeTestPoreModel(t)
	trueSeq := strings.Repeat("ACGTTGCAAGGCTTAC", 20)
	errSeq := trueSeq[:100] + "T" + trueSeq[101:] // one substitution
	other := strings.Repeat("GATTACA", 30)
	upper := "0F8FAD5B-D9CB-469F-A165-70867728950E"
	fq := fmt.Sprintf("@readA extra\n%s\n+\n%s\n@readB\n%s\n+\n%s\n@tiny\nACG\n+\nIII\n@%s\n%s\n+\n%s\n",
		errSeq, strings.Repeat("I", len(errSeq)), other, strings.Repeat("I", len(other)),
		upper, other, strings.Repeat("I", len(other)))
	input := filepath.Join(dir, "reads.fq.gz")
	writeGzip(t, input, fq)
	truthIn := filepath.Join(dir, "true.fa")
	if err := os.WriteFile(truthIn, []byte(">readA\n"+trueSeq[:150]+"\n"+trueSeq[150:]+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}

	out := filepath.Join(dir, "sig")
	stdout, stderr, err := runSignalCLI(t, "--input", input, "--truth", truthIn, "--out", out,
		"--pore-model", model, "--seed", "3", "--compression", "zstd")
	if err != nil {
		t.Fatalf("%v\n%s", err, stderr)
	}
	if !strings.Contains(stdout, "Wrote 3 reads") || !strings.Contains(stdout, "1 reads synthesised from --truth") ||
		!strings.Contains(stdout, "3 non-UUID read IDs mapped") || !strings.Contains(stderr, "skipped 1 reads shorter") {
		t.Errorf("stdout %q stderr %q", stdout, stderr)
	}
	_, hdrs, recs := readBLOW5Dir(t, out)
	if len(recs) != 3 {
		t.Fatalf("%d records, want 3", len(recs))
	}
	truth, sources, _ := readTruthFA(t, filepath.Join(out, SignalTruthName))
	wantSrc := []string{"readA", "readB", upper}
	wantSeq := []string{trueSeq, other, other}
	for i, r := range recs {
		if sources[r.ReadID] != wantSrc[i] || truth[r.ReadID] != wantSeq[i] {
			t.Errorf("record %d %s: source %q, truth seq ok=%v", i, r.ReadID, sources[r.ReadID], truth[r.ReadID] == wantSeq[i])
		}
	}
	// Non-UUID IDs get RFC 4122 version-5 UUIDs; a UUID is kept,
	// lowercased.
	if v := recs[0].ReadID; !isUUID(v) || v[14] != '5' || !strings.ContainsRune("89ab", rune(v[19])) {
		t.Errorf("mapped ID %s is not a version-5 UUID", v)
	}
	if recs[2].ReadID != strings.ToLower(upper) {
		t.Errorf("UUID read ID became %s", recs[2].ReadID)
	}
	// Run metadata synthesised from the seed: a 40-hex MinKNOW run ID and
	// channels on a PromethION flow cell.
	rg := hdrs[0].ReadGroups[0]
	if len(rg["run_id"]) != 40 || strings.Trim(rg["run_id"], "0123456789abcdef") != "" {
		t.Errorf("synthesised run_id %q", rg["run_id"])
	}
	for _, r := range recs {
		a, ok := slow5.ONTAuxOf(r)
		var ch int
		fmt.Sscan(a.ChannelNumber, &ch)
		if !ok || ch < 1 || ch > 2675 || a.ReadNumber < 1 {
			t.Errorf("%s aux %+v", r.ReadID, a)
		}
	}
	f, _ := os.Open(filepath.Join(out, slow5.ChunkFileName(rg["flow_cell_id"], rg["run_id"], 0)))
	defer f.Close()
	rd, err := slow5.NewReader(f)
	if err != nil || rd.RecordCompression() != slow5.RecordZstd || rd.SignalCompression() != slow5.SignalSvbZd {
		t.Errorf("reader %v: compression %v/%v, want zstd/svb-zd", err, rd.RecordCompression(), rd.SignalCompression())
	}
	m, err := manifest.Read(filepath.Join(out, manifest.FileName))
	if err != nil {
		t.Fatal(err)
	}
	if err := m.VerifyDir(out); err != nil {
		t.Error(err)
	}
	if m.Signal.Truth == nil || m.Signal.Truth.Name != "true.fa" || m.Signal.Skipped != 1 {
		t.Errorf("manifest signal %+v", m.Signal)
	}
}

func TestSignalCommand_OverwriteAndErrors(t *testing.T) {
	dir := t.TempDir()
	model := writeTestPoreModel(t)
	var fq strings.Builder
	for i := 0; i < 12; i++ {
		s := strings.Repeat("ACGGTCA", 10+i)
		fmt.Fprintf(&fq, "@r%d\n%s\n+\n%s\n", i, s, strings.Repeat("5", len(s)))
	}
	input := filepath.Join(dir, "r.fastq")
	if err := os.WriteFile(input, []byte(fq.String()), 0o600); err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(dir, "sig")
	base := []string{"--input", input, "--out", out, "--pore-model", model, "--seed", "5"}
	if _, stderr, err := runSignalCLI(t, append(base, "--reads-per-file", "2")...); err != nil {
		t.Fatalf("%v\n%s", err, stderr)
	}
	if files, _ := filepath.Glob(filepath.Join(out, "*.blow5")); len(files) != 6 {
		t.Fatalf("%d files, want 6", len(files))
	}
	if _, _, err := runSignalCLI(t, append(base, "--reads-per-file", "5")...); err == nil || !strings.Contains(err.Error(), "--overwrite") {
		t.Errorf("rerun into a used directory: %v", err)
	}
	// --overwrite replaces the run, leaving no stale chunks from the
	// larger file count.
	if _, stderr, err := runSignalCLI(t, append(base, "--reads-per-file", "5", "--overwrite")...); err != nil {
		t.Fatalf("%v\n%s", err, stderr)
	}
	files, _, recs := readBLOW5Dir(t, out)
	if len(files) != 3 || len(recs) != 12 {
		t.Errorf("after --overwrite: %d files, %d records; want 3, 12", len(files), len(recs))
	}
	m, err := manifest.Read(filepath.Join(out, manifest.FileName))
	if err != nil {
		t.Fatal(err)
	}
	if err := m.VerifyDir(out); err != nil {
		t.Error(err)
	}

	dup := filepath.Join(dir, "dup.fastq")
	if err := os.WriteFile(dup, []byte("@x\nACGTACGT\n+\nIIIIIIII\n@x\nACGTACGT\n+\nIIIIIIII\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	bad := filepath.Join(dir, "bad.fastq")
	if err := os.WriteFile(bad, []byte("@x\nACGTACGT\n+\nIII\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	cases := map[string]struct {
		args []string
		want string
	}{
		"no input":         {[]string{"--out", out, "--pore-model", model}, "--input"},
		"no model":         {[]string{"--input", input, "--out", out}, "--pore-model"},
		"bad compression":  {[]string{"--input", input, "--out", out, "--pore-model", model, "--compression", "lz4"}, "--compression"},
		"bad signal comp":  {[]string{"--input", input, "--out", out, "--pore-model", model, "--signal-compression", "vbz"}, "--signal-compression"},
		"bad profile":      {[]string{"--input", input, "--out", filepath.Join(dir, "p"), "--pore-model", model, "--profile", "r9"}, "unknown signal profile"},
		"duplicate IDs":    {[]string{"--input", dup, "--out", filepath.Join(dir, "d"), "--pore-model", model}, "must be unique"},
		"bad FASTQ":        {[]string{"--input", bad, "--out", filepath.Join(dir, "b"), "--pore-model", model}, "quality length"},
		"model units":      {[]string{"--input", input, "--out", out, "--pore-model", model, "--model-units", "mV"}, "--model-units"},
		"standardised k=5": {[]string{"--input", input, "--out", filepath.Join(dir, "s"), "--pore-model", model, "--model-units", "standardised", "--seed", "1"}, ""},
	}
	for name, c := range cases {
		_, _, err := runSignalCLI(t, c.args...)
		if c.want == "" {
			if err != nil {
				t.Errorf("%s: %v", name, err)
			}
			continue
		}
		if err == nil || !strings.Contains(err.Error(), c.want) {
			t.Errorf("%s: error %v, want %q", name, err, c.want)
		}
	}
}

// TestUUIDv5KnownVector checks the read-ID mapping against the RFC 4122
// name-based UUID of a published example: Python's
// uuid.uuid5(uuid.NAMESPACE_DNS, "python.org").
func TestUUIDv5KnownVector(t *testing.T) {
	dns := [16]byte{0x6b, 0xa7, 0xb8, 0x10, 0x9d, 0xad, 0x11, 0xd1, 0x80, 0xb4, 0x00, 0xc0, 0x4f, 0xd4, 0x30, 0xc8}
	if got := formatUUID(uuidV5(dns, "python.org")); got != "886313e1-3b8a-5372-9b90-0c9aee199e5d" {
		t.Errorf("uuid5(DNS, python.org) = %s", got)
	}
	if id, mapped := signalReadID("A00123:8:H7:1:1101:1000:2000"); !mapped || !isUUID(id) {
		t.Errorf("signalReadID mapped=%v id=%s", mapped, id)
	}
	if a, _ := signalReadID("r1"); a != formatUUID(uuidV5(signalUUIDNamespace, "r1")) {
		t.Errorf("signalReadID(r1) = %s", a)
	}
}

func TestParseONTComment(t *testing.T) {
	f := parseONTComment("runid=abc read=7 ch=12 start_time=2024-03-01T12:00:00.250000+00:00 flow_cell_id=PAW1 protocol_group_id=pg sample_id=s1 junk")
	want := time.Date(2024, 3, 1, 12, 0, 0, 250000000, time.UTC)
	if !f.placed() || f.runID != "abc" || f.read != 7 || f.ch != 12 || !f.start.Equal(want) || f.flowCell != "PAW1" || f.sampleID != "s1" || f.protocolGroup != "pg" {
		t.Errorf("%+v", f)
	}
	if g := parseONTComment("1:N:0:1"); g.placed() || g.runID != "" {
		t.Errorf("Illumina comment parsed as %+v", g)
	}
}
