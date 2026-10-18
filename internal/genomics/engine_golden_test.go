//go:build !release

package genomics

// Golden/structure tests that pin the shape of Illumina output: file names,
// record counts, header formats, read and quality lengths, the quality string
// of the linear-decay model, and the RG tag. They were written before the
// read pipeline refactor so the refactor provably preserves the output
// (only the randomness source changes).

import (
	"bytes"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// illuminaHeader matches @<sample>:<lane>:HXXXXXX:<batch>:1101:<x>:<y> <mate>:N:0:1[\tRG:Z:<id>]
var illuminaHeader = regexp.MustCompile(`^@([^:]+):(\d+):HXXXXXX:(\d+):1101:(\d+):(\d+) ([12]):N:0:1(\tRG:Z:(\S+))?$`)

var acgt = regexp.MustCompile(`^[ACGT]+$`)

func goldenConfig(dir string) types.SequencerConfig {
	return types.SequencerConfig{
		SampleName:    "Gold",
		Lane:          3,
		OutputDir:     dir,
		ReadLength:    60,
		ReadsPerBatch: 7,
		TotalBatches:  2,
		BaseQuality:   38,
		QualityDecay:  0.25,
	}
}

// decayQuality is the expected ASCII quality string of the linear-decay model.
func decayQuality(c types.SequencerConfig) string {
	const chars = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
	var b strings.Builder
	for i := 0; i < c.ReadLength; i++ {
		q := c.BaseQuality - int(float64(i)*c.QualityDecay)
		if q < 2 {
			q = 2
		}
		if q > len(chars)-1 {
			q = len(chars) - 1
		}
		b.WriteByte(chars[q])
	}
	return b.String()
}

type fastqRec struct{ header, seq, plus, qual string }

func parseFASTQ(t *testing.T, data string) []fastqRec {
	t.Helper()
	lines := strings.Split(strings.TrimSuffix(data, "\n"), "\n")
	if len(lines)%4 != 0 {
		t.Fatalf("FASTQ has %d lines, not a multiple of 4", len(lines))
	}
	var recs []fastqRec
	for i := 0; i < len(lines); i += 4 {
		recs = append(recs, fastqRec{lines[i], lines[i+1], lines[i+2], lines[i+3]})
	}
	return recs
}

func checkRecord(t *testing.T, c types.SequencerConfig, r fastqRec, batch, mate int, x, y int, rg string) {
	t.Helper()
	m := illuminaHeader.FindStringSubmatch(r.header)
	if m == nil {
		t.Fatalf("header %q does not match the Illumina format", r.header)
	}
	if m[1] != c.SampleName || m[2] != strconv.Itoa(c.Lane) || m[3] != strconv.Itoa(batch) || m[6] != strconv.Itoa(mate) {
		t.Errorf("header %q: want sample %s lane %d batch %d mate %d", r.header, c.SampleName, c.Lane, batch, mate)
	}
	if x >= 0 && (m[4] != strconv.Itoa(x) || m[5] != strconv.Itoa(y)) {
		t.Errorf("header %q: want x=%d y=%d", r.header, x, y)
	}
	if m[8] != rg {
		t.Errorf("header %q: RG = %q, want %q", r.header, m[8], rg)
	}
	if len(r.seq) != c.ReadLength || !acgt.MatchString(r.seq) {
		t.Errorf("sequence %q: want %d ACGT bases", r.seq, c.ReadLength)
	}
	if r.plus != "+" {
		t.Errorf("separator %q, want +", r.plus)
	}
	if want := decayQuality(c); r.qual != want {
		t.Errorf("quality %q, want %q", r.qual, want)
	}
}

func TestGolden_FileOutputSingleEnd(t *testing.T) {
	c := goldenConfig(t.TempDir())
	e := NewEngine(c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	var names []string
	entries, _ := os.ReadDir(c.OutputDir)
	for _, de := range entries {
		names = append(names, de.Name())
	}
	sort.Strings(names)
	want := []string{"Gold_S1_L003_001.fastq", "Gold_S1_L003_002.fastq"}
	if strings.Join(names, ",") != strings.Join(want, ",") {
		t.Fatalf("files %v, want %v", names, want)
	}
	if got := e.GeneratedFiles(); len(got) != 2 || filepath.Base(got[0]) != want[0] || filepath.Base(got[1]) != want[1] {
		t.Errorf("GeneratedFiles = %v", got)
	}
	for bi, name := range want {
		data, err := os.ReadFile(filepath.Join(c.OutputDir, name))
		if err != nil {
			t.Fatal(err)
		}
		recs := parseFASTQ(t, string(data))
		if len(recs) != c.ReadsPerBatch {
			t.Fatalf("%s: %d records, want %d", name, len(recs), c.ReadsPerBatch)
		}
		for i, r := range recs {
			n := i + 1
			checkRecord(t, c, r, bi+1, 1, n*1000, n*2000, "")
		}
	}
}

func TestGolden_FileOutputPairedEnd(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.PairedEnd = true
	c.TotalBatches = 1
	e := NewEngine(c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	for mate, name := range map[int]string{1: "Gold_S1_L003_R1_001.fastq", 2: "Gold_S1_L003_R2_001.fastq"} {
		data, err := os.ReadFile(filepath.Join(c.OutputDir, name))
		if err != nil {
			t.Fatal(err)
		}
		recs := parseFASTQ(t, string(data))
		if len(recs) != c.ReadsPerBatch {
			t.Fatalf("%s: %d records", name, len(recs))
		}
		for i, r := range recs {
			checkRecord(t, c, r, 1, mate, (i+1)*1000, (i+1)*2000, "")
		}
	}
}

func TestGolden_GenerateFASTQFile(t *testing.T) {
	c := goldenConfig(t.TempDir())
	e := NewEngine(c)
	defer e.Close()
	path := filepath.Join(c.OutputDir, "custom.fastq")
	if err := e.GenerateFASTQFile(path, 5, 2); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	recs := parseFASTQ(t, string(data))
	if len(recs) != c.ReadsPerBatch {
		t.Fatalf("%d records", len(recs))
	}
	for i, r := range recs {
		checkRecord(t, c, r, 5, 2, (i+1)*1000, (i+1)*2000, "")
	}
}

func TestGolden_WriterInterleavedPairedWithRG(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.PairedEnd = true
	e := NewEngine(c)
	defer e.Close()
	var buf bytes.Buffer
	if err := e.GenerateToWriter(&buf, GenerateOpts{SampleID: "tumor"}); err != nil {
		t.Fatal(err)
	}
	recs := parseFASTQ(t, buf.String())
	if len(recs) != c.TotalBatches*c.ReadsPerBatch*2 {
		t.Fatalf("%d records", len(recs))
	}
	i := 0
	for batch := 1; batch <= c.TotalBatches; batch++ {
		for n := 1; n <= c.ReadsPerBatch; n++ {
			for mate := 1; mate <= 2; mate++ {
				checkRecord(t, c, recs[i], batch, mate, n*1000, n*2000, "tumor")
				i++
			}
		}
	}
}

func TestGolden_WriterRegionHeaderEncodesStart(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.ReadsPerBatch = 200
	c.TotalBatches = 1
	e := NewEngine(c)
	defer e.Close()
	region := &GenomicRegion{Chromosome: "chr7", Start: 123_400_000, End: 123_400_500}
	var buf bytes.Buffer
	if err := e.GenerateToWriter(&buf, GenerateOpts{RefRegion: region}); err != nil {
		t.Fatal(err)
	}
	recs := parseFASTQ(t, buf.String())
	if len(recs) != c.ReadsPerBatch {
		t.Fatalf("%d records", len(recs))
	}
	for _, r := range recs {
		checkRecord(t, c, r, 1, 1, -1, -1, "")
		m := illuminaHeader.FindStringSubmatch(r.header)
		x, _ := strconv.ParseInt(m[4], 10, 64)
		y, _ := strconv.ParseInt(m[5], 10, 64)
		// Header encodes readStart%1e5 and (readStart/1e5)%1e5.
		start := (region.Start/10_000_000_000)*10_000_000_000 + y*100000 + x
		if start < region.Start || start >= region.End-int64(c.ReadLength) {
			t.Errorf("decoded read start %d outside [%d,%d)", start, region.Start, region.End-int64(c.ReadLength))
		}
	}
}

func TestGolden_WriterRegionTooSmall(t *testing.T) {
	c := goldenConfig(t.TempDir())
	e := NewEngine(c)
	defer e.Close()
	err := e.GenerateToWriter(&bytes.Buffer{}, GenerateOpts{RefRegion: &GenomicRegion{Chromosome: "chr1", Start: 0, End: 60}})
	if err == nil || !strings.Contains(err.Error(), "too small") {
		t.Errorf("want region too small error, got %v", err)
	}
}
