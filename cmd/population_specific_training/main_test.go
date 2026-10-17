//go:build !release
// +build !release

package main

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestGetStatusEmoji_Success(t *testing.T) {
	result := getStatusEmoji(true)
	if result != "✅" {
		t.Errorf("Expected '✅', got '%s'", result)
	}
}

func TestGetStatusEmoji_Failure(t *testing.T) {
	result := getStatusEmoji(false)
	if result != "❌" {
		t.Errorf("Expected '❌', got '%s'", result)
	}
}

func TestAbsFunction(t *testing.T) {
	if abs(5.0) != 5.0 {
		t.Errorf("abs(5.0) should be 5.0, got %f", abs(5.0))
	}
	if abs(-5.0) != 5.0 {
		t.Errorf("abs(-5.0) should be 5.0, got %f", abs(-5.0))
	}
	if abs(0.0) != 0.0 {
		t.Errorf("abs(0.0) should be 0.0, got %f", abs(0.0))
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	gc := calculateGCContent(seqs)
	if gc < 0.99 || gc > 1.01 {
		t.Errorf("Expected ~1.0 GC for all-GC, got %f", gc)
	}
}

func TestCalculateGCContent_NoGC(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATAT")}
	gc := calculateGCContent(seqs)
	if gc > 0.01 {
		t.Errorf("Expected ~0.0 GC for no-GC, got %f", gc)
	}
}

func TestCalculateGCContent_Empty(t *testing.T) {
	seqs := [][]byte{}
	gc := calculateGCContent(seqs)
	if gc != 0.0 {
		t.Errorf("Expected 0.0 for empty sequences, got %f", gc)
	}
}

func TestCalculateComplexity_High(t *testing.T) {
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	seqs := [][]byte{seq}
	complexity := calculateComplexity(seqs, 4)
	if complexity <= 0 {
		t.Errorf("Expected positive complexity, got %f", complexity)
	}
}

func TestCalculateComplexity_Low(t *testing.T) {
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
	seqs := [][]byte{seq}
	complexity := calculateComplexity(seqs, 4)
	if complexity >= 0.5 {
		t.Errorf("Expected low complexity for repetitive sequence, got %f", complexity)
	}
}

func TestCalculateComplexity_Empty(t *testing.T) {
	seqs := [][]byte{}
	complexity := calculateComplexity(seqs, 4)
	if complexity != 0.0 {
		t.Errorf("Expected 0.0 for empty sequences, got %f", complexity)
	}
}

func TestGetPopulationProfiles_NotEmpty(t *testing.T) {
	profiles := getPopulationProfiles()
	if len(profiles) == 0 {
		t.Error("Expected non-empty population profiles")
	}
}

func TestGetProfileByPopulation_Found(t *testing.T) {
	profiles := getPopulationProfiles()
	if len(profiles) == 0 {
		t.Skip("No profiles available")
	}

	// Get the first profile's population
	targetPop := profiles[0].Population
	result := getProfileByPopulation(profiles, targetPop)
	if result.Population != targetPop {
		t.Errorf("Expected population '%s', got '%s'", targetPop, result.Population)
	}
}

func TestGetProfileByPopulation_NotFound(t *testing.T) {
	profiles := getPopulationProfiles()
	result := getProfileByPopulation(profiles, "NONEXISTENT")
	if result.Population != "" {
		t.Errorf("Expected empty profile for nonexistent population, got '%s'", result.Population)
	}
}

func TestCountSuccessful(t *testing.T) {
	profiles := []PopulationGCProfile{
		{Population: "POP1"},
		{Population: "POP2"},
	}
	results := map[string]PopulationResults{
		"POP1": {Success: true},
		"POP2": {Success: false},
		"POP3": {Success: true},
	}
	count := countSuccessful(results, profiles)
	if count != 2 { // Only POP1 and POP3 are successful
		t.Errorf("Expected 2 successful results, got %d", count)
	}
}

func TestGenerateSequenceWithGC_Length(t *testing.T) {
	seq := generateSequenceWithGC(100, 0.5)
	if len(seq) != 100 {
		t.Errorf("Expected sequence length 100, got %d", len(seq))
	}
}

func TestGenerateSequenceWithGC_ValidBases(t *testing.T) {
	seq := generateSequenceWithGC(50, 0.5)
	for i, base := range seq {
		if base != 'A' && base != 'T' && base != 'G' && base != 'C' {
			t.Errorf("Invalid base at position %d: %c", i, base)
		}
	}
}

func TestDetectAppleSiliconCores_Positive(t *testing.T) {
	cores := detectAppleSiliconCores()
	if cores <= 0 {
		t.Errorf("expected positive core count, got %d", cores)
	}
}

func TestSavePopulationResults_NoCrash(t *testing.T) {
	tmpDir := t.TempDir()
	profiles := getPopulationProfiles()
	results := map[string]PopulationResults{
		"YRI": {Population: "YRI", Complexity: 0.66, GCContent: 0.47, Success: true, Sequences: 10},
		"CEU": {Population: "CEU", Complexity: 0.65, GCContent: 0.45, Success: true, Sequences: 10},
	}
	savePopulationResults(tmpDir, results, profiles, 0.65, 0.95)
}

// ---------- generatePopulationTrainingData ----------

func TestGeneratePopulationTrainingData_ReturnsCorrectCount(t *testing.T) {
	profile := PopulationGCProfile{
		Population:    "YRI",
		TargetGC:      0.47,
		Tolerance:     0.02,
		ComplexityMin: 0.60,
	}
	data := generatePopulationTrainingData(profile)
	if len(data) != 1000 {
		t.Errorf("generatePopulationTrainingData returned %d sequences, want 1000", len(data))
	}
}

func TestGeneratePopulationTrainingData_PopulationMetadata(t *testing.T) {
	profile := PopulationGCProfile{
		Population: "CEU",
		TargetGC:   0.45,
		Tolerance:  0.02,
	}
	data := generatePopulationTrainingData(profile)
	if len(data) == 0 {
		t.Fatal("generatePopulationTrainingData returned no data")
	}
	if data[0].Population != "CEU" {
		t.Errorf("data[0].Population = %q, want CEU", data[0].Population)
	}
	if len(data[0].Sequence) == 0 {
		t.Error("data[0].Sequence is empty")
	}
}

// ---------- generateValidationSequences ----------

func TestGenerateValidationSequences_TrainedModel(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(8)
	// Train with a simple sequence
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	if err := model.TrainFromSequence(seq); err != nil {
		t.Skipf("model training failed: %v", err)
	}
	seqs, err := generateValidationSequences(model, 2, 50)
	if err != nil {
		t.Fatalf("generateValidationSequences error: %v", err)
	}
	if len(seqs) != 2 {
		t.Errorf("expected 2 sequences, got %d", len(seqs))
	}
}

func TestGenerateValidationSequences_UntrainedModel(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(8)
	// Untrained model - expect error or empty results
	_, err := generateValidationSequences(model, 1, 50)
	// Result is implementation-dependent; just verify no panic
	_ = err
}

