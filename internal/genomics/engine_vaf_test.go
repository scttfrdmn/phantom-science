//go:build !release

package genomics

// Tests for issue #17: exact Bernoulli allele draws, indel-safe positioned
// injection and no 'N' padding. They use only the public API.

import (
	"bytes"
	"math"
	"strconv"
	"strings"
	"testing"
)

// fixedStartRun generates n reads that all start at regionStart (the region is
// exactly one base longer than the read) and returns the parsed records.
func fixedStartRun(t *testing.T, n, readLen int, regionStart int64, variants []VariantRecord) []fastqRec {
	t.Helper()
	c := goldenConfig(t.TempDir())
	c.ReadLength = readLen
	c.ReadsPerBatch = n
	c.TotalBatches = 1
	e := NewEngine(c)
	defer e.Close()
	region := &GenomicRegion{Chromosome: "chr1", Start: regionStart, End: regionStart + int64(readLen) + 1}
	var buf bytes.Buffer
	if err := e.GenerateToWriter(&buf, GenerateOpts{Variants: variants, RefRegion: region}); err != nil {
		t.Fatal(err)
	}
	recs := parseFASTQ(t, buf.String())
	if len(recs) != n {
		t.Fatalf("%d records, want %d", len(recs), n)
	}
	for _, r := range recs {
		if strings.ContainsRune(r.seq, 'N') {
			t.Fatalf("read contains N: %q", r.seq)
		}
		if len(r.seq) != readLen || len(r.qual) != readLen {
			t.Fatalf("read/quality length %d/%d, want %d", len(r.seq), len(r.qual), readLen)
		}
		if got := headerReadStart(t, r.header); got != regionStart {
			t.Fatalf("read start %d, want %d", got, regionStart)
		}
	}
	return recs
}

// headerReadStart decodes the read start encoded in a region-mode header
// (x = start % 1e5, y = (start / 1e5) % 1e5); valid for starts below 1e10.
func headerReadStart(t *testing.T, header string) int64 {
	t.Helper()
	m := illuminaHeader.FindStringSubmatch(header)
	if m == nil {
		t.Fatalf("bad header %q", header)
	}
	x, _ := strconv.ParseInt(m[4], 10, 64)
	y, _ := strconv.ParseInt(m[5], 10, 64)
	return y*100000 + x
}

// altCount counts reads carrying alt at the given read offset.
func altCount(recs []fastqRec, off int, alt string) int {
	n := 0
	for _, r := range recs {
		if strings.HasPrefix(r.seq[off:], alt) {
			n++
		}
	}
	return n
}

func TestInjectAtPosition_FreqOneAlwaysAlt(t *testing.T) {
	const n = 10000
	// SNV at 1-based 1101 → read offset 25 for reads starting at 0-based 1075.
	recs := fixedStartRun(t, n, 50, 1075, []VariantRecord{{Chromosome: "chr1", Position: 1101, Ref: "A", Alt: "T", Freq: 1.0}})
	refs := 0
	for _, r := range recs {
		if r.seq[25] != 'T' {
			refs++
		}
	}
	if refs != 0 {
		t.Errorf("Freq=1.0: %d of %d reads carry REF, want 0", refs, n)
	}
}

func TestInjectAtPosition_FreqZeroNeverAlt(t *testing.T) {
	const n = 10000
	recs := fixedStartRun(t, n, 50, 1075, []VariantRecord{{Chromosome: "chr1", Position: 1101, Ref: "A", Alt: "TGCATGCA", Freq: 0}})
	if got := altCount(recs, 25, "TGCATGCA"); got != 0 {
		t.Errorf("Freq=0: %d reads carry ALT, want 0", got)
	}
	for _, r := range recs {
		if r.seq[25] != 'A' {
			t.Fatalf("Freq=0: REF base not enforced at offset 25: %q", r.seq)
		}
	}
}

func TestInjectAtPosition_LowVAFBinomial(t *testing.T) {
	const n = 10000
	for _, p := range []float64{0.3, 0.02} {
		recs := fixedStartRun(t, n, 50, 1075, []VariantRecord{{Chromosome: "chr1", Position: 1101, Ref: "A", Alt: "T", Freq: p}})
		alt := 0
		for _, r := range recs {
			if r.seq[25] == 'T' {
				alt++
			}
		}
		got := float64(alt) / n
		tol := 4 * math.Sqrt(p*(1-p)/n)
		if math.Abs(got-p) > tol {
			t.Errorf("Freq=%v: observed VAF %.4f outside ±%.4f", p, got, tol)
		}
	}
}

// Two deletions and one insertion in one read, all at Freq=1: every allele
// must sit at its expected offset after the indels to its left.
func TestInjectAtPosition_MultipleIndelsInOneRead(t *testing.T) {
	const start = 2000 // 0-based read start
	variants := []VariantRecord{
		// Out of order on purpose; positions are 1-based.
		{Chromosome: "chr1", Position: start + 61, Ref: "G", Alt: "GTTTTTT", Freq: 1}, // insertion +6
		{Chromosome: "chr1", Position: start + 11, Ref: "CAGTA", Alt: "C", Freq: 1},   // deletion -4
		{Chromosome: "chr1", Position: start + 31, Ref: "AGG", Alt: "A", Freq: 1},     // deletion -2
	}
	recs := fixedStartRun(t, 500, 100, start, variants)
	// Expected offsets in the output read: original offset plus the net
	// length change of all variants to its left.
	type want struct {
		off    int
		allele string
	}
	wants := []want{
		{10, "C"},
		{30 - 4, "A"},
		{60 - 4 - 2, "GTTTTTT"},
	}
	for _, r := range recs {
		for _, w := range wants {
			if !strings.HasPrefix(r.seq[w.off:], w.allele) {
				t.Fatalf("allele %q not at offset %d in %q", w.allele, w.off, r.seq)
			}
		}
	}
}

// Overlapping variants: the rightmost is applied, the overlapping one to its
// left is skipped and counted, and the read is unharmed.
func TestInjectAtPosition_OverlapSkippedAndCounted(t *testing.T) {
	seq := []byte(strings.Repeat("A", 40))
	variants := []VariantRecord{ // sorted by descending position
		{Chromosome: "chr1", Position: 21, Ref: "AAAA", Alt: "GGGG", Freq: 1},
		{Chromosome: "chr1", Position: 19, Ref: "AAA", Alt: "C", Freq: 1}, // span 18..21 overlaps 20..24
		{Chromosome: "chr1", Position: 5, Ref: "A", Alt: "T", Freq: 1},
	}
	got, q, skipped := injectAtPosition(testRand(), seq, make([]byte, 40), 0, variants)
	if skipped != 1 {
		t.Errorf("skipped = %d, want 1", skipped)
	}
	want := strings.Repeat("A", 4) + "T" + strings.Repeat("A", 15) + "GGGG" + strings.Repeat("A", 16)
	if string(got) != want || len(q) != len(got) {
		t.Errorf("got %q, want %q", got, want)
	}
}

func TestGenerateToWriter_OtherChromosomeIgnored(t *testing.T) {
	recs := fixedStartRun(t, 200, 50, 1075, []VariantRecord{
		{Chromosome: "chr2", Position: 1101, Ref: "A", Alt: "TTTTTTTTTT", Freq: 1},
		{Chromosome: "1", Position: 1091, Ref: "A", Alt: "CCCCCCCCCC", Freq: 1}, // "1" matches "chr1"
	})
	if altCount(recs, 25, "TTTTTTTTTT") != 0 {
		t.Error("variant on another chromosome was injected")
	}
	if altCount(recs, 15, "CCCCCCCCCC") != len(recs) {
		t.Error("variant on the region chromosome (without chr prefix) was not injected")
	}
}

// A read whose window contains a deletion is extended from generated bases,
// never padded with N, and the position-blind path behaves the same way.
func TestGenerateToWriter_NoNPadding(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.ReadsPerBatch = 2000
	c.TotalBatches = 1
	e := NewEngine(c)
	defer e.Close()
	var buf bytes.Buffer
	// Position-blind: delete any 2-mer "AC" → "A".
	if err := e.GenerateToWriter(&buf, GenerateOpts{Variants: []VariantRecord{{Ref: "AC", Alt: "A", Freq: 1}, {Ref: "GT", Alt: "G", Freq: 1}}}); err != nil {
		t.Fatal(err)
	}
	for _, r := range parseFASTQ(t, buf.String()) {
		if strings.ContainsRune(r.seq, 'N') || len(r.seq) != c.ReadLength || len(r.qual) != c.ReadLength {
			t.Fatalf("bad record %q / %q", r.seq, r.qual)
		}
	}
}
