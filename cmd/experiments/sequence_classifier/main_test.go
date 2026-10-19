//go:build !release
// +build !release

package main

import (
	"os"
	"testing"
)

func TestIsValidDNA_AllValid(t *testing.T) {
	if !isValidDNA("ATCG") {
		t.Error("ATCG should be valid DNA")
	}
	if !isValidDNA("GCGCGCGC") {
		t.Error("GCGCGCGC should be valid DNA")
	}
}

func TestIsValidDNA_Empty(t *testing.T) {
	if !isValidDNA("") {
		t.Error("empty string should be valid DNA")
	}
}

func TestIsValidDNA_WithN(t *testing.T) {
	if isValidDNA("ATCGN") {
		t.Error("ATCGN should not be valid DNA (N is invalid)")
	}
}

func TestIsValidDNA_Lowercase(t *testing.T) {
	if isValidDNA("atcg") {
		t.Error("lowercase should not be valid DNA")
	}
}

func TestIsValidDNA_WithDigit(t *testing.T) {
	if isValidDNA("ATC1") {
		t.Error("ATC1 should not be valid DNA (digit)")
	}
}

func TestAbs_Positive(t *testing.T) {
	if got := abs(5.0); got != 5.0 {
		t.Errorf("abs(5.0) = %f, want 5.0", got)
	}
}

func TestAbs_Negative(t *testing.T) {
	if got := abs(-5.0); got != 5.0 {
		t.Errorf("abs(-5.0) = %f, want 5.0", got)
	}
}

func TestAbs_Zero(t *testing.T) {
	if got := abs(0.0); got != 0.0 {
		t.Errorf("abs(0.0) = %f, want 0.0", got)
	}
}

func TestGetStatusFromConfidence_HighConfidence(t *testing.T) {
	if got := getStatusFromConfidence(0.9); got != "High Confidence" {
		t.Errorf("getStatusFromConfidence(0.9) = %q, want \"High Confidence\"", got)
	}
}

func TestGetStatusFromConfidence_MediumConfidence(t *testing.T) {
	if got := getStatusFromConfidence(0.6); got != "Medium Confidence" {
		t.Errorf("getStatusFromConfidence(0.6) = %q, want \"Medium Confidence\"", got)
	}
}

func TestGetStatusFromConfidence_LowConfidence(t *testing.T) {
	if got := getStatusFromConfidence(0.4); got != "Low Confidence" {
		t.Errorf("getStatusFromConfidence(0.4) = %q, want \"Low Confidence\"", got)
	}
}

func TestGetStatusFromConfidence_Unclassified(t *testing.T) {
	if got := getStatusFromConfidence(0.1); got != "Unclassified" {
		t.Errorf("getStatusFromConfidence(0.1) = %q, want \"Unclassified\"", got)
	}
}

func TestNewPopulationClassifier_NotNil(t *testing.T) {
	pc := NewPopulationClassifier()
	if pc == nil {
		t.Fatal("NewPopulationClassifier() returned nil")
	}
}

func TestNewPopulationClassifier_HasPopulationGCTargets(t *testing.T) {
	pc := NewPopulationClassifier()
	if len(pc.PopulationGCTargets) == 0 {
		t.Error("Expected non-empty PopulationGCTargets map")
	}
}

func TestPopulationClassifier_CalculateGCContent_Empty(t *testing.T) {
	pc := NewPopulationClassifier()
	if got := pc.calculateGCContent(""); got != 0.0 {
		t.Errorf("calculateGCContent(\"\") = %f, want 0.0", got)
	}
}

func TestPopulationClassifier_CalculateGCContent_AllGC(t *testing.T) {
	pc := NewPopulationClassifier()
	if got := pc.calculateGCContent("GCGCGCGC"); got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestPopulationClassifier_CalculateGCContent_AllAT(t *testing.T) {
	pc := NewPopulationClassifier()
	if got := pc.calculateGCContent("ATATATAT"); got != 0.0 {
		t.Errorf("calculateGCContent(all AT) = %f, want 0.0", got)
	}
}

func TestPopulationClassifier_CalculateGCContent_HalfGC(t *testing.T) {
	pc := NewPopulationClassifier()
	got := pc.calculateGCContent("GCATGCAT")
	if got < 0.49 || got > 0.51 {
		t.Errorf("calculateGCContent(half GC) = %f, want ~0.5", got)
	}
}

func TestPopulationClassifier_CalculateComplexity_Empty(t *testing.T) {
	pc := NewPopulationClassifier()
	if got := pc.calculateComplexity(""); got != 0.0 {
		t.Errorf("calculateComplexity(\"\") = %f, want 0.0", got)
	}
}

func TestPopulationClassifier_CalculateComplexity_TooShort(t *testing.T) {
	pc := NewPopulationClassifier()
	if got := pc.calculateComplexity("ATG"); got != 0.0 {
		t.Errorf("calculateComplexity(3-char) = %f, want 0.0", got)
	}
}

func TestPopulationClassifier_CalculateComplexity_LowComplexity(t *testing.T) {
	pc := NewPopulationClassifier()
	got := pc.calculateComplexity("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
	if got >= 0.5 {
		t.Errorf("calculateComplexity(all-A) = %f, expected < 0.5", got)
	}
}

func TestPopulationClassifier_CalculateComplexity_HighComplexity(t *testing.T) {
	pc := NewPopulationClassifier()
	got := pc.calculateComplexity("ATCGATCGATCGATCGATCGATCGATCGATCG")
	if got <= 0 {
		t.Errorf("calculateComplexity(diverse) = %f, expected > 0", got)
	}
}

func TestPopulationClassifier_CalculateGCContent_AllN(t *testing.T) {
	// Sequence with no A/T/G/C chars → totalCount == 0 → returns 0.0
	pc := NewPopulationClassifier()
	if got := pc.calculateGCContent("NNNNNNNN"); got != 0.0 {
		t.Errorf("calculateGCContent(all-N) = %f, want 0.0", got)
	}
}

func TestClassifySequence_NonExistentFile(t *testing.T) {
	pc := NewPopulationClassifier()
	_, _, err := pc.ClassifySequence("/nonexistent/path/to/file.fasta")
	if err == nil {
		t.Error("ClassifySequence with non-existent file should return an error")
	}
}

func TestClassifySequence_LowScoringSequence(t *testing.T) {
	// All-N sequence → gcContent=0, complexity=0 → all pop scores < 0 → bestScore stays 0.0
	// → triggers "if bestScore < 0.3" branch → returns "Unknown"
	f, err := os.CreateTemp("", "test_lowscore_*.fasta")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(f.Name())
	if _, err := f.WriteString(">seq1\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"); err != nil {
		t.Fatalf("failed to write: %v", err)
	}
	f.Close()

	pc := NewPopulationClassifier()
	pop, score, err := pc.ClassifySequence(f.Name())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if pop != "Unknown" {
		t.Errorf("expected 'Unknown' for all-N sequence, got %q", pop)
	}
	if score != 0.0 {
		t.Errorf("expected score 0.0 for unknown, got %f", score)
	}
}

func TestClassifySequence_ValidFasta(t *testing.T) {
	f, err := os.CreateTemp("", "test.fasta")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(f.Name())

	if _, err := f.WriteString(">seq1\nACGTACGTGCGCGCGCGCATATATGCGCATATACGTACGTGCGCGCGCGCATATATGCGCATAT\n"); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	f.Close()

	pc := NewPopulationClassifier()
	pop, score, err := pc.ClassifySequence(f.Name())
	if err != nil {
		t.Errorf("ClassifySequence with valid FASTA returned unexpected error: %v", err)
	}
	if pop == "" {
		t.Error("ClassifySequence returned empty population string")
	}
	if score < 0.0 || score > 1.0 {
		t.Errorf("ClassifySequence returned score %f, want value in [0, 1]", score)
	}
}

// TestClassifySequence_Deterministic: all GC targets are equal, so every
// population ties; the reported population must not depend on map
// iteration order.
func TestClassifySequence_Deterministic(t *testing.T) {
	p := t.TempDir() + "/seq.fasta"
	seq := ">seq1\nACGTTGCAAGCTTCGAGGATCCATGCAAGTCCGATAGCTAGGCTAACGTTGCAAGCTTCGAGG\n"
	if err := os.WriteFile(p, []byte(seq), 0o644); err != nil {
		t.Fatal(err)
	}
	pc := NewPopulationClassifier()
	first, _, err := pc.ClassifySequence(p)
	if err != nil {
		t.Fatal(err)
	}
	if first == "Unknown" {
		t.Fatal("fixture scored below the threshold; choose a sequence that scores")
	}
	for i := 0; i < 50; i++ {
		if got, _, _ := pc.ClassifySequence(p); got != first {
			t.Fatalf("call %d returned %s, first call %s", i, got, first)
		}
	}
}
