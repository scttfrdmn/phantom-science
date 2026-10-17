//go:build !release
// +build !release

package main

import (
	"math"
	"os"
	"runtime"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// ---------- calculateGCContent ([][]byte) ----------

func TestCalculateGCContent_Empty(t *testing.T) {
	if got := calculateGCContent([][]byte{}); got != 0.0 {
		t.Errorf("calculateGCContent(empty) = %v, want 0", got)
	}
}

func TestCalculateGCContent_AllAT(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATATAT")}
	if got := calculateGCContent(seqs); got != 0.0 {
		t.Errorf("calculateGCContent(all AT) = %v, want 0", got)
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGCGC")}
	if got := calculateGCContent(seqs); got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %v, want 1.0", got)
	}
}

func TestCalculateGCContent_Mixed(t *testing.T) {
	seqs := [][]byte{[]byte("AATTGGCC")} // 4 GC out of 8 = 0.5
	got := calculateGCContent(seqs)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("calculateGCContent(AATTGGCC) = %v, want 0.5", got)
	}
}

func TestCalculateGCContent_MultipleSeqs(t *testing.T) {
	seqs := [][]byte{
		[]byte("GGGG"), // all GC
		[]byte("AAAA"), // no GC
	} // 4 GC out of 8 total = 0.5
	got := calculateGCContent(seqs)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("calculateGCContent(two seqs) = %v, want 0.5", got)
	}
}

// ---------- calculateComplexity ----------

func TestCalculateComplexity_Empty(t *testing.T) {
	if got := calculateComplexity([][]byte{}); got != 0.0 {
		t.Errorf("calculateComplexity(empty) = %v, want 0", got)
	}
}

func TestCalculateComplexity_ShortSeqs(t *testing.T) {
	// Sequences shorter than k=4 should contribute 0 each
	seqs := [][]byte{[]byte("AT"), []byte("GC")}
	got := calculateComplexity(seqs)
	if got != 0.0 {
		t.Errorf("calculateComplexity(short seqs) = %v, want 0", got)
	}
}

func TestCalculateComplexity_SingleKmer(t *testing.T) {
	// Exactly k=4 bases -> 1 unique kmer / 1 total = 1.0
	seqs := [][]byte{[]byte("ATGC")}
	got := calculateComplexity(seqs)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("calculateComplexity(ATGC) = %v, want 1.0", got)
	}
}

func TestCalculateComplexity_RepetitiveSeq(t *testing.T) {
	// AAAAAAAA k=4: 5 total kmers, 1 unique -> 1/5 = 0.2
	seqs := [][]byte{[]byte("AAAAAAAA")}
	got := calculateComplexity(seqs)
	want := 1.0 / 5.0
	if math.Abs(got-want) > 1e-9 {
		t.Errorf("calculateComplexity(repetitive) = %v, want %v", got, want)
	}
}

func TestCalculateComplexity_HighComplexity(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCGATCG")}
	got := calculateComplexity(seqs)
	if got <= 0 || got > 1.0 {
		t.Errorf("calculateComplexity should be in (0,1], got %v", got)
	}
}

func TestCalculateComplexity_MultipleSeqs(t *testing.T) {
	seqs := [][]byte{
		[]byte("ATGCATGC"),
		[]byte("GCTAGCTA"),
	}
	got := calculateComplexity(seqs)
	if got <= 0 || got > 1.0 {
		t.Errorf("calculateComplexity(multiple seqs) = %v, out of range (0,1]", got)
	}
}

// ---------- detectAppleSiliconCores ----------

func TestDetectAppleSiliconCores_NonZero(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %v, want > 0", got)
	}
}

func TestDetectAppleSiliconCores_NotExceedTotal(t *testing.T) {
	got := detectAppleSiliconCores()
	total := runtime.NumCPU()
	if got > total {
		t.Errorf("detectAppleSiliconCores() = %v, exceeds NumCPU = %v", got, total)
	}
}

// TestCalculateGCContent_EmptySequenceInSlice covers the inner len(seq)==0
// branch, which skips an empty sequence entry but still divides by len(sequences).
func TestCalculateGCContent_EmptySequenceInSlice(t *testing.T) {
	seqs := [][]byte{{}, []byte("GCGC")} // empty + 100% GC
	// totalGC = 1.0 (GCGC only), divisor = 2 → 0.5
	got := calculateGCContent(seqs)
	want := 0.5
	if got != want {
		t.Errorf("calculateGCContent([{}, GCGC]) = %v, want %v", got, want)
	}
}

// ---------- trainIntensively ----------

func TestTrainIntensively_MissingFile(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	err := trainIntensively(model, "/nonexistent/path/file.fasta", 1)
	if err == nil {
		t.Error("trainIntensively(missing file) expected error, got nil")
	}
}

func TestTrainIntensively_WithFASTAFile(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := tmpDir + "/test.fasta"
	// Create a minimal FASTA file with a simple sequence
	content := ">chr1\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n>chr2\nGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n"
	if err := os.WriteFile(fastaPath, []byte(content), 0600); err != nil {
		t.Fatalf("failed to create FASTA file: %v", err)
	}
	model := authenticity.NewAdaptiveHMMModel(4)
	err := trainIntensively(model, fastaPath, 1)
	if err != nil {
		t.Errorf("trainIntensively(valid FASTA) error: %v", err)
	}
}

func TestTrainIntensively_MultipleWorkers(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := tmpDir + "/multi.fasta"
	content := ">seq1\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	if err := os.WriteFile(fastaPath, []byte(content), 0600); err != nil {
		t.Fatalf("failed to create FASTA file: %v", err)
	}
	model := authenticity.NewAdaptiveHMMModel(4)
	err := trainIntensively(model, fastaPath, runtime.NumCPU())
	if err != nil {
		t.Errorf("trainIntensively(multiple workers) error: %v", err)
	}
}
