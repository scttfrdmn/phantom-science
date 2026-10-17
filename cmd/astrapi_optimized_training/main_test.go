//go:build !release
// +build !release

package main

import (
	"math"
	"os"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func TestPassFailIcon_Pass(t *testing.T) {
	if got := passFailIcon(true); got != "✅" {
		t.Errorf("passFailIcon(true) = %q, want ✅", got)
	}
}

func TestPassFailIcon_Fail(t *testing.T) {
	if got := passFailIcon(false); got != "❌" {
		t.Errorf("passFailIcon(false) = %q, want ❌", got)
	}
}

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	if got := calculateSequenceComplexity([]byte{}, 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(empty) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_TooShort(t *testing.T) {
	if got := calculateSequenceComplexity([]byte("ATG"), 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(short) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_Empty(t *testing.T) {
	if got := calculateGCContent([][]byte{}); got != 0.0 {
		t.Errorf("calculateGCContent([]) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	if got := calculateGCContent(seqs); got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestCalculateGCContent_HalfGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCATGCAT")}
	got := calculateGCContent(seqs)
	if math.Abs(got-0.5) > 0.01 {
		t.Errorf("calculateGCContent(half GC) = %f, want ~0.5", got)
	}
}

func TestFindMaxHomopolymerRun_Empty(t *testing.T) {
	if got := findMaxHomopolymerRun([][]byte{}); got != 0 {
		t.Errorf("findMaxHomopolymerRun([]) = %d, want 0", got)
	}
}

func TestFindMaxHomopolymerRun_WithRun(t *testing.T) {
	seqs := [][]byte{[]byte("ATAAAACG")}
	if got := findMaxHomopolymerRun(seqs); got != 4 {
		t.Errorf("findMaxHomopolymerRun(AAAA run) = %d, want 4", got)
	}
}

func TestExtractKmerCounts_Empty(t *testing.T) {
	counts := extractKmerCounts([][]byte{}, 4)
	if len(counts) != 0 {
		t.Errorf("expected empty counts, got %v", counts)
	}
}

func TestNormalizeCounts_SumsToOne(t *testing.T) {
	counts := map[string]int{"AAAA": 3, "TTTT": 1, "GCGC": 2}
	probs := normalizeCounts(counts)
	total := 0.0
	for _, p := range probs {
		total += p
	}
	if math.Abs(total-1.0) > 1e-9 {
		t.Errorf("probs sum = %f, want 1.0", total)
	}
}

func TestJensenShannonDivergence_Identical(t *testing.T) {
	p := map[string]float64{"AAAA": 0.5, "TTTT": 0.5}
	q := map[string]float64{"AAAA": 0.5, "TTTT": 0.5}
	got := jensenShannonDivergence(p, q)
	if math.Abs(got) > 1e-9 {
		t.Errorf("JSD(identical) = %f, want ~0", got)
	}
}

func TestJensenShannonDivergence_Disjoint(t *testing.T) {
	p := map[string]float64{"AAAA": 1.0}
	q := map[string]float64{"CCCC": 1.0}
	got := jensenShannonDivergence(p, q)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("JSD(disjoint) = %f, want 1.0", got)
	}
}

// ---------- calculateSequenceComplexity (with actual kmer loop) ----------

func TestCalculateSequenceComplexity_RepetitiveLoop(t *testing.T) {
	// AAAAAAA -> 1 unique kmer of 4 out of 4 total
	seq := []byte("AAAAAAA") // 7 bases, k=4: 4 kmers all "AAAA"
	got := calculateSequenceComplexity(seq, 4)
	want := 1.0 / 4.0
	if math.Abs(got-want) > 1e-9 {
		t.Errorf("complexity(repetitive) = %f, want %f", got, want)
	}
}

func TestCalculateSequenceComplexity_HighDiversity(t *testing.T) {
	seq := []byte("ATCGATCGATCGATCGATCGATCG")
	got := calculateSequenceComplexity(seq, 4)
	if got <= 0 || got > 1.0 {
		t.Errorf("complexity = %f, want (0,1]", got)
	}
}

// ---------- calculateCpGDepletion ----------

func TestCalculateCpGDepletion_Empty(t *testing.T) {
	got := calculateCpGDepletion([][]byte{})
	if got != 0.0 {
		t.Errorf("cpgDepletion(empty) = %f, want 0.0", got)
	}
}

func TestCalculateCpGDepletion_NoCpG(t *testing.T) {
	seqs := [][]byte{[]byte("AAAATTTTAAAATTTT")}
	got := calculateCpGDepletion(seqs)
	if got != 0.0 {
		t.Errorf("cpgDepletion(no CpG) = %f, want 0.0", got)
	}
}

func TestCalculateCpGDepletion_WithCpG(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCGATCG")}
	got := calculateCpGDepletion(seqs)
	if got < 0 {
		t.Errorf("cpgDepletion should be non-negative, got %f", got)
	}
}

// ---------- calculateJensenShannonDivergence ----------

func TestCalculateJensenShannonDivergence_Identical(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCGATCG")}
	got := calculateJensenShannonDivergence(seqs, seqs, 4)
	if got < 0 {
		t.Errorf("JSD(identical) should be non-negative, got %f", got)
	}
}

func TestCalculateJensenShannonDivergence_Empty(t *testing.T) {
	got := calculateJensenShannonDivergence([][]byte{}, [][]byte{}, 4)
	if got != 0.0 {
		t.Errorf("JSD(empty,empty) = %f, want 0.0", got)
	}
}

func TestCalculateJensenShannonDivergence_DifferentSeqs(t *testing.T) {
	p := [][]byte{[]byte("AAAAAAAAAAAAAAAA")}
	q := [][]byte{[]byte("CCCCCCCCCCCCCCCC")}
	got := calculateJensenShannonDivergence(p, q, 4)
	if got < 0 || got > 1.0+1e-9 {
		t.Errorf("JSD = %f, should be in [0,1]", got)
	}
}

// ---------- extractKmerCounts (more coverage) ----------

func TestExtractKmerCounts_WithSeqs(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCG")}
	counts := extractKmerCounts(seqs, 4)
	// ATCG, TCGA, CGAT, GATC, ATCG -> 4 unique, 5 total
	if len(counts) == 0 {
		t.Error("expected non-empty counts")
	}
	if counts["ATCG"] == 0 {
		t.Error("expected ATCG kmer to be present")
	}
}

// ---------- loadReferenceSequences ----------

func TestLoadReferenceSequences_EmptyDir_ReturnsError(t *testing.T) {
	// With empty tmpDir, LoadChromosome("chr1") will fail
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	_, err := loadReferenceSequences(refGenome, 5, 1000)
	if err == nil {
		t.Error("loadReferenceSequences(empty refGenome) expected error, got nil")
	}
}

func TestLoadReferenceSequences_WithChr1FASTA(t *testing.T) {
	tmpDir := t.TempDir()
	seq := make([]byte, 500)
	for i := range seq {
		seq[i] = []byte("ATCG")[i%4]
	}
	fastaContent := ">chr1 Homo sapiens chromosome 1\n" + string(seq) + "\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	seqs, err := loadReferenceSequences(refGenome, 3, 100)
	if err != nil {
		t.Fatalf("loadReferenceSequences error: %v", err)
	}
	if len(seqs) == 0 {
		t.Error("expected at least 1 sequence")
	}
}

// Ensure math and os imports are used
var _ = math.Abs
var _ = os.Stat
