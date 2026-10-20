//go:build !release

package genomics

// ONT headers and sequencing summary, gzip output and file patterns (#12).

import (
	"bytes"
	"compress/gzip"
	"io"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// ontHeader matches a MinKNOW header; the sample_id field may carry
// read-ID fingerprint tags appended to the line.
var ontHeader = regexp.MustCompile(`^@([0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}) runid=([0-9a-f]{40}) read=(\d+) ch=(\d+) start_time=(\S+) flow_cell_id=(\S+) protocol_group_id=(\S+) sample_id=(\S+)(\tRG:Z:\S+)?(\tXC:Z:\S+)?(\tXD:Z:\S+)?$`)

func ontConfig(dir string) types.SequencerConfig {
	return types.SequencerConfig{
		Platform:      "ont",
		SampleName:    "ontS",
		Lane:          1,
		OutputDir:     dir,
		ReadLength:    3000,
		ReadsPerBatch: 40,
		TotalBatches:  3,
		BaseQuality:   40,
		QualityDecay:  0.01,
		Seed:          77,
		ONT:           types.ONTSpec{Channels: 16, ProtocolGroupID: "pg1"},
	}
}

// readMaybeGzip returns a file's content, decompressing .gz files and
// checking that they really are gzip.
func readMaybeGzip(t *testing.T, path string) []byte {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.HasSuffix(path, ".gz") {
		return data
	}
	zr, err := gzip.NewReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("%s is not gzip: %v", filepath.Base(path), err)
	}
	plain, err := io.ReadAll(zr)
	if err != nil {
		t.Fatalf("%s: %v", filepath.Base(path), err)
	}
	return plain
}

type summaryRow struct {
	file, id, run string
	ch, length    int
	start, dur    time.Duration
	q             float64
}

func readSummary(t *testing.T, path string) []summaryRow {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("reading summary: %v", err)
	}
	lines := strings.Split(strings.TrimSuffix(string(data), "\n"), "\n")
	if lines[0] != strings.Join(nanopore.SummaryColumns, "\t") {
		t.Fatalf("summary header %q", lines[0])
	}
	secs := func(s string) time.Duration {
		f, err := strconv.ParseFloat(s, 64)
		if err != nil {
			t.Fatalf("bad seconds %q", s)
		}
		return time.Duration(math.Round(f*1e6)) * time.Microsecond
	}
	var rows []summaryRow
	for _, l := range lines[1:] {
		c := strings.Split(l, "\t")
		if len(c) != len(nanopore.SummaryColumns) {
			t.Fatalf("summary row %q has %d columns", l, len(c))
		}
		ch, _ := strconv.Atoi(c[3])
		n, _ := strconv.Atoi(c[6])
		q, _ := strconv.ParseFloat(c[7], 64)
		rows = append(rows, summaryRow{c[0], c[1], c[2], ch, n, secs(c[4]), secs(c[5]), q})
	}
	return rows
}

// ontMeanQ recomputes the ONT mean qscore from a FASTQ quality string.
func ontMeanQ(qual string) float64 {
	var sum float64
	for i := 0; i < len(qual); i++ {
		sum += math.Pow(10, -float64(qual[i]-'!')/10)
	}
	return -10 * math.Log10(sum/float64(len(qual)))
}

// checkONTRun cross-checks every FASTQ file of a run against its summary:
// one row per read, same IDs, channels, start times and lengths, mean
// qscores recomputed from the written qualities, and no two reads on a
// channel overlapping in time across the whole run.
func checkONTRun(t *testing.T, e *Engine, c types.SequencerConfig, files []string, summary string) {
	t.Helper()
	info := e.ontRunInfo()
	rows := readSummary(t, summary)
	byID := map[string]summaryRow{}
	for _, r := range rows {
		if _, dup := byID[r.id]; dup {
			t.Fatalf("read %s has two summary rows", r.id)
		}
		byID[r.id] = r
	}
	type span struct {
		read       int
		start, end time.Duration
	}
	byCh := map[int][]span{}
	reads := 0
	for _, f := range files {
		recs := parseFASTQ(t, string(readMaybeGzip(t, f)))
		if len(recs) != c.ReadsPerBatch {
			t.Fatalf("%s: %d reads, want %d", filepath.Base(f), len(recs), c.ReadsPerBatch)
		}
		for _, rec := range recs {
			reads++
			m := ontHeader.FindStringSubmatch(rec.header)
			if m == nil {
				t.Fatalf("header %q is not a MinKNOW header", rec.header)
			}
			if m[2] != info.RunID || m[6] != info.FlowCellID || m[7] != info.ProtocolGroupID || !strings.HasPrefix(m[8], c.SampleName) {
				t.Errorf("header %q: run fields differ from %+v", rec.header, info)
			}
			row, ok := byID[m[1]]
			if !ok {
				t.Fatalf("read %s is missing from the summary", m[1])
			}
			st, err := time.Parse(time.RFC3339Nano, m[5])
			if err != nil {
				t.Fatalf("start_time %q: %v", m[5], err)
			}
			if row.file != filepath.Base(f) || row.run != info.RunID || strconv.Itoa(row.ch) != m[4] ||
				!info.Start.Add(row.start).Equal(st) || row.length != len(rec.seq) || row.length != c.ReadLength {
				t.Errorf("read %s: summary row %+v disagrees with %s (file %s, %d bases)", m[1], row, rec.header, filepath.Base(f), len(rec.seq))
			}
			if want := ontMeanQ(rec.qual); math.Abs(row.q-want) > 1e-5 {
				t.Errorf("read %s: mean_qscore %v, recomputed %v", m[1], row.q, want)
			}
			n, _ := strconv.Atoi(m[3])
			byCh[row.ch] = append(byCh[row.ch], span{n, row.start, row.start + row.dur})
		}
	}
	if reads != len(rows) {
		t.Errorf("%d reads in FASTQ, %d summary rows", reads, len(rows))
	}
	for ch, s := range byCh {
		if ch < 1 || ch > info.Channels {
			t.Errorf("channel %d outside 1..%d", ch, info.Channels)
		}
		sort.Slice(s, func(i, j int) bool { return s[i].start < s[j].start })
		for i := 1; i < len(s); i++ {
			if s[i].start < s[i-1].end || s[i].read <= s[i-1].read {
				t.Fatalf("channel %d: read %d at %v overlaps or precedes read %d ending %v", ch, s[i].read, s[i].start, s[i-1].read, s[i-1].end)
			}
		}
	}
}

func TestONT_FilesAndSummaryAgree(t *testing.T) {
	c := ontConfig(t.TempDir())
	e := NewEngine(c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	files := e.GeneratedFiles()
	if len(files) != c.TotalBatches {
		t.Fatalf("GeneratedFiles = %v", files)
	}
	info := e.ontRunInfo()
	for i, f := range files {
		want := "ontS_" + info.FlowCellID + "_" + info.RunID[:8] + "_" + strconv.Itoa(i) + ".fastq"
		if filepath.Base(f) != want {
			t.Errorf("file %d is %s, want %s", i, filepath.Base(f), want)
		}
	}
	checkONTRun(t, e, c, files, filepath.Join(c.OutputDir, nanopore.SummaryFileName))
}

func TestONT_GzipPatternMatchesPlain(t *testing.T) {
	plainCfg := ontConfig(t.TempDir())
	plain := NewEngine(plainCfg)
	defer plain.Close()
	if err := plain.Generate(); err != nil {
		t.Fatal(err)
	}

	c := ontConfig(t.TempDir())
	c.Output = types.OutputSpec{
		Compression: "gzip",
		FilePattern: "{sample_name}_{barcode}_{flowcell_id}_{file_number}.fastq.gz",
		SummaryFile: "summary.tsv",
	}
	c.ONT.FlowCellID = "PAW00001"
	e := NewEngine(c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	files := e.GeneratedFiles()
	entries, _ := os.ReadDir(c.OutputDir)
	var names []string
	for _, de := range entries {
		names = append(names, de.Name())
	}
	want := []string{"ontS_barcode01_PAW00001_0.fastq.gz", "ontS_barcode01_PAW00001_1.fastq.gz", "ontS_barcode01_PAW00001_2.fastq.gz", "phantom-manifest.json", "summary.tsv"}
	if strings.Join(names, ",") != strings.Join(want, ",") {
		t.Fatalf("output files %v, want %v (no plain FASTQ left behind)", names, want)
	}
	for i, f := range files {
		if filepath.Base(f) != want[i] {
			t.Errorf("GeneratedFiles[%d] = %s, want %s", i, f, want[i])
		}
		// Only the flow-cell ID differs, and it is not drawn from the read streams.
		got := strings.ReplaceAll(string(readMaybeGzip(t, f)), "flow_cell_id=PAW00001", "flow_cell_id="+plain.ontRunInfo().FlowCellID)
		if ref := string(readMaybeGzip(t, plain.GeneratedFiles()[i])); got != ref {
			t.Errorf("%s decompressed differs from the plain output of the same seed", filepath.Base(f))
		}
	}
	checkONTRun(t, e, c, files, filepath.Join(c.OutputDir, "summary.tsv"))
}

// seededGolden is goldenConfig with a fixed seed, so runs can be compared.
func seededGolden(dir string) types.SequencerConfig {
	c := goldenConfig(dir)
	c.Seed = 2024
	return c
}

func TestIllumina_GzipRoundTrip(t *testing.T) {
	ref := fileOutputs(t, seededGolden(t.TempDir()))
	c := seededGolden(t.TempDir())
	c.Output.Compression = "gzip"
	e := NewEngine(c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	files := e.GeneratedFiles()
	if len(files) != len(ref) {
		t.Fatalf("GeneratedFiles = %v, reference %d files", files, len(ref))
	}
	for _, f := range files {
		plainName := strings.TrimSuffix(filepath.Base(f), ".gz")
		if plainName == filepath.Base(f) {
			t.Errorf("%s lacks .gz", f)
		}
		if _, err := os.Stat(filepath.Join(c.OutputDir, plainName)); err == nil {
			t.Errorf("plain %s left behind", plainName)
		}
		if !bytes.Equal(readMaybeGzip(t, f), ref[plainName]) {
			t.Errorf("%s decompressed differs from the uncompressed output", filepath.Base(f))
		}
	}
	if _, err := os.Stat(filepath.Join(c.OutputDir, nanopore.SummaryFileName)); err == nil {
		t.Error("Illumina run wrote a sequencing summary")
	}
}

func TestGzip_StreamAndCallerNamedFile(t *testing.T) {
	c := seededGolden(t.TempDir())
	ref := streamOutput(t, c, GenerateOpts{SampleID: "rg"})
	c.Output.Compression = "gzip"
	gz := streamOutput(t, c, GenerateOpts{SampleID: "rg"})
	zr, err := gzip.NewReader(bytes.NewReader(gz))
	if err != nil {
		t.Fatalf("streamed output is not gzip: %v", err)
	}
	if plain, _ := io.ReadAll(zr); !bytes.Equal(plain, ref) {
		t.Error("decompressed stream differs from the plain stream")
	}

	// GenerateFASTQFile with a caller-chosen plain name (the TUI) ends up as .gz.
	e := NewEngine(c)
	defer e.Close()
	name := filepath.Join(c.OutputDir, "tui.fastq")
	if err := e.GenerateFASTQFile(name, 1, 1); err != nil {
		t.Fatal(err)
	}
	if got := e.GeneratedFiles(); len(got) != 1 || got[0] != name+".gz" {
		t.Fatalf("GeneratedFiles = %v", got)
	}
	if _, err := os.Stat(name); err == nil {
		t.Error("plain file left behind")
	}
	readMaybeGzip(t, name+".gz")
}

func TestONT_StreamSummaryOnlyWhenRequested(t *testing.T) {
	c := ontConfig(t.TempDir())
	c.TotalBatches = 2
	out := streamOutput(t, c, GenerateOpts{})
	recs := parseFASTQ(t, string(out))
	if len(recs) != 2*c.ReadsPerBatch {
		t.Fatalf("%d reads", len(recs))
	}
	for _, r := range recs {
		if !ontHeader.MatchString(r.header) {
			t.Fatalf("streamed header %q", r.header)
		}
	}
	if entries, _ := os.ReadDir(c.OutputDir); len(entries) != 0 {
		t.Errorf("streamed run without summary_file wrote %d files", len(entries))
	}

	c.Output.SummaryFile = "stream_summary.txt"
	again := streamOutput(t, c, GenerateOpts{})
	if !bytes.Equal(again, out) {
		t.Error("requesting a summary changed the reads")
	}
	rows := readSummary(t, filepath.Join(c.OutputDir, "stream_summary.txt"))
	if len(rows) != len(recs) {
		t.Fatalf("%d summary rows for %d reads", len(rows), len(recs))
	}
	for i, r := range recs {
		if m := ontHeader.FindStringSubmatch(r.header); m[1] != rows[i].id || rows[i].file != "" || rows[i].length != len(r.seq) {
			t.Errorf("row %d %+v vs %q", i, rows[i], r.header)
		}
	}
}

func TestONT_ReadIDFingerprintKeepsReadID(t *testing.T) {
	c := ontConfig(t.TempDir())
	c.TotalBatches = 1
	c.DatasetID = "ds9"
	c.Enterprise.CustomerID = "cust7"
	c.Fingerprinting = types.FingerprintConfig{Enable: true, CustomerFingerprint: true, DatasetFingerprint: true}
	e := NewEngine(c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	recs := parseFASTQ(t, string(readMaybeGzip(t, e.GeneratedFiles()[0])))
	// The tags are SAM-style fields after the MinKNOW comment (#27): the
	// read ID and sample_id are unchanged.
	for _, r := range recs {
		if !strings.HasSuffix(r.header, " sample_id=ontS\tXC:Z:cust7\tXD:Z:ds9") {
			t.Fatalf("header %q: want sample_id=ontS then the XC/XD tags", r.header)
		}
	}
	checkONTRun(t, e, c, e.GeneratedFiles(), filepath.Join(c.OutputDir, nanopore.SummaryFileName))
}

func TestONT_Reproducible(t *testing.T) {
	run := func() map[string][]byte {
		c := ontConfig(t.TempDir())
		e := NewEngine(c)
		defer e.Close()
		if err := e.Generate(); err != nil {
			t.Fatal(err)
		}
		out := map[string][]byte{}
		for _, f := range append(e.GeneratedFiles(), filepath.Join(c.OutputDir, nanopore.SummaryFileName)) {
			out[filepath.Base(f)] = readMaybeGzip(t, f)
		}
		return out
	}
	a, b := run(), run()
	if len(a) != len(b) || len(a) != 4 { // three FASTQ files and the summary
		t.Fatalf("outputs %d and %d files", len(a), len(b))
	}
	for name := range a {
		if !bytes.Equal(a[name], b[name]) {
			t.Errorf("%s differs between runs with the same seed", name)
		}
	}
}

func TestOutputConfig_Rejected(t *testing.T) {
	cases := map[string]func(*types.SequencerConfig){
		"ont paired":          func(c *types.SequencerConfig) { c.Platform = "ont"; c.PairedEnd = true },
		"ont bad speed":       func(c *types.SequencerConfig) { c.Platform = "ont"; c.ONT.SpeedBPS = -1 },
		"ont bad channels":    func(c *types.SequencerConfig) { c.Platform = "ont"; c.ONT.Channels = -5 },
		"ont spaced run id":   func(c *types.SequencerConfig) { c.Platform = "ont"; c.ONT.RunID = "a b" },
		"ont spaced sample":   func(c *types.SequencerConfig) { c.Platform = "ont"; c.SampleName = "my sample" },
		"summary on illumina": func(c *types.SequencerConfig) { c.Output.SummaryFile = "s.txt" },
		"summary with dir":    func(c *types.SequencerConfig) { c.Platform = "ont"; c.Output.SummaryFile = "../s.txt" },
		"unknown placeholder": func(c *types.SequencerConfig) { c.Output.FilePattern = "{sample}_{file_number}.fastq" },
		"pattern with dir":    func(c *types.SequencerConfig) { c.Output.FilePattern = "sub/{file_number}.fastq" },
		"no file_number":      func(c *types.SequencerConfig) { c.Output.FilePattern = "{sample_name}.fastq" },
		"paired no mate":      func(c *types.SequencerConfig) { c.PairedEnd = true; c.Output.FilePattern = "{file_number}.fastq" },
		"gz without gzip":     func(c *types.SequencerConfig) { c.Output.FilePattern = "{file_number}.fastq.gz" },
		"unbalanced brace":    func(c *types.SequencerConfig) { c.Output.FilePattern = "{file_number}{.fastq" },
		"empty placeholder":   func(c *types.SequencerConfig) { c.Output.FilePattern = "{run_id}_{file_number}.fastq" },
	}
	for name, edit := range cases {
		t.Run(name, func(t *testing.T) {
			c := goldenConfig(t.TempDir())
			edit(&c)
			e := NewEngine(c)
			defer e.Close()
			if err := e.Generate(); err == nil {
				t.Error("Generate accepted the configuration")
			}
			if entries, _ := os.ReadDir(c.OutputDir); len(entries) != 0 {
				t.Errorf("output written despite the error: %d files", len(entries))
			}
		})
	}
}

func TestIllumina_FilePattern(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.PairedEnd = true
	c.Output.FilePattern = "{sample_name}_L{lane}_R{mate}_{file_number}.fq"
	e := NewEngine(c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	var names []string
	for _, f := range e.GeneratedFiles() {
		names = append(names, filepath.Base(f))
	}
	if got := strings.Join(names, ","); got != "Gold_L3_R1_0.fq,Gold_L3_R2_0.fq,Gold_L3_R1_1.fq,Gold_L3_R2_1.fq" {
		t.Errorf("files %s", got)
	}
}
