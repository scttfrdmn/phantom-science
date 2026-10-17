//go:build !release
// +build !release

package main

import (
	"math"
	"os"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// Note: calculateGCContent returns NaN for empty sequences (0/0) in this package.

func TestIsValidKmer_Valid(t *testing.T) {
	if !isValidKmer("ATCG") {
		t.Error("ATCG should be valid kmer")
	}
	if !isValidKmer("GCGCGCGC") {
		t.Error("GCGCGCGC should be valid kmer")
	}
}

func TestIsValidKmer_Empty(t *testing.T) {
	if !isValidKmer("") {
		t.Error("empty string should be valid kmer")
	}
}

func TestIsValidKmer_WithN(t *testing.T) {
	if isValidKmer("ATCGN") {
		t.Error("ATCGN should not be valid kmer (N is invalid)")
	}
}

func TestIsValidKmer_Lowercase(t *testing.T) {
	if isValidKmer("atcg") {
		t.Error("lowercase should not be valid kmer")
	}
}

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	if got := calculateSequenceComplexity([]byte{}, 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(empty) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_TooShort(t *testing.T) {
	if got := calculateSequenceComplexity([]byte("ATG"), 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(too short) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_Empty(t *testing.T) {
	// This implementation returns NaN for empty sequences (0/0)
	got := calculateGCContent([][]byte{})
	if !math.IsNaN(got) && got != 0.0 {
		t.Errorf("calculateGCContent([]) = %f, want NaN or 0.0", got)
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	if got := calculateGCContent(seqs); got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %f, want 1.0", got)
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

func TestCalculateJSDivergenceFromProbs_Identical(t *testing.T) {
	p := map[string]float64{"AAAA": 0.5, "TTTT": 0.5}
	q := map[string]float64{"AAAA": 0.5, "TTTT": 0.5}
	got := calculateJSDivergenceFromProbs(p, q)
	if math.Abs(got) > 1e-9 {
		t.Errorf("JSD(identical) = %f, want ~0", got)
	}
}

func TestCalculateJSDivergenceFromProbs_Disjoint(t *testing.T) {
	p := map[string]float64{"AAAA": 1.0}
	q := map[string]float64{"CCCC": 1.0}
	got := calculateJSDivergenceFromProbs(p, q)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("JSD(disjoint) = %f, want 1.0", got)
	}
}

func TestCalculateCpGDepletion_Empty(t *testing.T) {
	if got := calculateCpGDepletion([][]byte{}); got != 0.0 {
		t.Errorf("calculateCpGDepletion([]) = %f, want 0.0", got)
	}
}

func TestCalculateCpGDepletion_NoCpG(t *testing.T) {
	seqs := [][]byte{[]byte("AAAATTTTAAAATTTT")}
	if got := calculateCpGDepletion(seqs); got != 0.0 {
		t.Errorf("calculateCpGDepletion(no CpG) = %f, want 0.0", got)
	}
}

func TestCalculateCpGDepletion_WithCpG(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCGATCG")}
	got := calculateCpGDepletion(seqs)
	if got < 0 {
		t.Errorf("calculateCpGDepletion should be non-negative, got %f", got)
	}
}

func TestCalculateSequenceComplexity_WithSeq(t *testing.T) {
	seq := []byte("ATCGATCGATCGATCG")
	got := calculateSequenceComplexity(seq, 4)
	if got <= 0 || got > 1.0 {
		t.Errorf("complexity = %f, want (0,1]", got)
	}
}

func TestCalculateSequenceComplexity_Repetitive(t *testing.T) {
	seq := make([]byte, 100)
	for i := range seq {
		seq[i] = 'A'
	}
	got := calculateSequenceComplexity(seq, 4)
	if got >= 0.5 {
		t.Errorf("repetitive complexity = %f, want < 0.5", got)
	}
}

func TestCalculateJSDivergence_EmptySeqsEmptyRef(t *testing.T) {
	got := calculateJSDivergence([][]byte{}, []string{}, 4)
	if got < 0 {
		t.Errorf("JSD should be non-negative, got %f", got)
	}
}

func TestCalculateJSDivergence_WithSeqsAndRef(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCGATCG")}
	refKmers := []string{"ATCG", "TCGA", "CGAT", "GATC", "ATCG"}
	got := calculateJSDivergence(seqs, refKmers, 4)
	if got < 0 || got > 1.0+1e-9 {
		t.Errorf("JSD = %f, should be in [0,1]", got)
	}
}

func TestCalculateJSDivergence_DisjointKmers(t *testing.T) {
	seqs := [][]byte{[]byte("AAAAAAAAAAAAAAAA")}
	refKmers := []string{"CCCC", "CCCC", "CCCC"}
	got := calculateJSDivergence(seqs, refKmers, 4)
	if got < 0 || got > 1.0+1e-9 {
		t.Errorf("JSD = %f, should be in [0,1]", got)
	}
}

// ---------- extractLimitedReferenceKmers ----------

func TestExtractLimitedReferenceKmers_EmptyDir_ReturnsEmpty(t *testing.T) {
	// With empty tmpDir, chr22 cannot be loaded; the function should return empty slice
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	kmers := extractLimitedReferenceKmers(refGenome, 4, 100)
	if kmers == nil {
		kmers = []string{}
	}
	if len(kmers) != 0 {
		t.Errorf("extractLimitedReferenceKmers(empty dir) = %d kmers, want 0", len(kmers))
	}
}

func TestExtractLimitedReferenceKmers_WithChr22FASTA(t *testing.T) {
	tmpDir := t.TempDir()
	fastaContent := ">chr22 Homo sapiens chromosome 22\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	kmers := extractLimitedReferenceKmers(refGenome, 4, 10)
	if len(kmers) == 0 {
		t.Error("extractLimitedReferenceKmers should return kmers when chr22 loads successfully")
	}
	if len(kmers) > 10 {
		t.Errorf("extractLimitedReferenceKmers returned %d kmers, want at most 10", len(kmers))
	}
}

// Ensure os import is used
var _ = os.Stat
var _ = math.Abs
