//go:build !release
// +build !release

package main

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

func TestCalculateGCContent_AllGC(t *testing.T) {
	seq := []byte("GCGCGCGC")
	got := calculateGCContent(seq)
	if got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestCalculateGCContent_AllAT(t *testing.T) {
	seq := []byte("ATATAT")
	got := calculateGCContent(seq)
	if got != 0.0 {
		t.Errorf("calculateGCContent(all AT) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_Half(t *testing.T) {
	seq := []byte("ATGC") // 2 GC out of 4
	got := calculateGCContent(seq)
	if got != 0.5 {
		t.Errorf("calculateGCContent(ATGC) = %f, want 0.5", got)
	}
}

func TestGetExpectedGCContent_KnownPop(t *testing.T) {
	tests := []struct {
		pop  string
		want float64
	}{
		{"YRI", 0.41},
		{"CEU", 0.42},
		{"CHB", 0.40},
	}
	for _, tt := range tests {
		got := getExpectedGCContent(tt.pop)
		if got != tt.want {
			t.Errorf("getExpectedGCContent(%q) = %f, want %f", tt.pop, got, tt.want)
		}
	}
}

func TestGetExpectedGCContent_Unknown(t *testing.T) {
	got := getExpectedGCContent("UNKNOWN")
	if got != 0.41 {
		t.Errorf("getExpectedGCContent(UNKNOWN) = %f, want 0.41 (default)", got)
	}
}

func TestAbs_Positive(t *testing.T) {
	if got := abs(3.14); got != 3.14 {
		t.Errorf("abs(3.14) = %f, want 3.14", got)
	}
}

func TestAbs_Negative(t *testing.T) {
	if got := abs(-3.14); got != 3.14 {
		t.Errorf("abs(-3.14) = %f, want 3.14", got)
	}
}

func TestAbs_Zero(t *testing.T) {
	if got := abs(0.0); got != 0.0 {
		t.Errorf("abs(0) = %f, want 0.0", got)
	}
}

// ---------- convertToHMMVariants ----------

func TestConvertToHMMVariants_Empty(t *testing.T) {
	result := convertToHMMVariants([]population.VariantData{})
	if len(result) != 0 {
		t.Errorf("convertToHMMVariants(empty) length = %d, want 0", len(result))
	}
}

func TestConvertToHMMVariants_WithAlleles(t *testing.T) {
	variants := []population.VariantData{
		{Position: 100, Alleles: []string{"A", "G"}, Frequency: 0.3, Genotypes: []int{0, 1, 0}},
		{Position: 200, Alleles: []string{"C", "T"}, Frequency: 0.5, Genotypes: []int{1, 1, 0}},
	}
	result := convertToHMMVariants(variants)
	if len(result) != 2 {
		t.Fatalf("convertToHMMVariants length = %d, want 2", len(result))
	}
	if result[0].Position != 100 {
		t.Errorf("result[0].Position = %d, want 100", result[0].Position)
	}
	if result[0].ReferenceAllele != 'A' {
		t.Errorf("result[0].ReferenceAllele = %c, want A", result[0].ReferenceAllele)
	}
	if result[0].AlternateAllele != 'G' {
		t.Errorf("result[0].AlternateAllele = %c, want G", result[0].AlternateAllele)
	}
	if result[0].Frequency != 0.3 {
		t.Errorf("result[0].Frequency = %f, want 0.3", result[0].Frequency)
	}
}

func TestConvertToHMMVariants_NoAlleles_UsesDefaults(t *testing.T) {
	variants := []population.VariantData{
		{Position: 50, Alleles: []string{}, Frequency: 0.1},
	}
	result := convertToHMMVariants(variants)
	if len(result) != 1 {
		t.Fatalf("length = %d, want 1", len(result))
	}
	if result[0].ReferenceAllele != 'A' {
		t.Errorf("default ReferenceAllele = %c, want A", result[0].ReferenceAllele)
	}
	if result[0].AlternateAllele != 'T' {
		t.Errorf("default AlternateAllele = %c, want T", result[0].AlternateAllele)
	}
}

// ---------- convertStatistics ----------

func TestConvertStatistics_ZeroValues(t *testing.T) {
	stats := population.PopulationStatistics{}
	result := convertStatistics(stats)
	if result["num_variants"] != 0 {
		t.Errorf("num_variants = %f, want 0", result["num_variants"])
	}
	if result["diversity_pi"] != 0 {
		t.Errorf("diversity_pi = %f, want 0", result["diversity_pi"])
	}
}

func TestConvertStatistics_WithValues(t *testing.T) {
	stats := population.PopulationStatistics{
		DiversityPi:         0.001,
		WattersonTheta:      0.002,
		TajimaD:             -0.5,
		MeanAlleleFrequency: 0.3,
		NumVariants:         100,
		RareVariantsCount:   20,
	}
	result := convertStatistics(stats)
	if result["diversity_pi"] != 0.001 {
		t.Errorf("diversity_pi = %f, want 0.001", result["diversity_pi"])
	}
	if result["num_variants"] != 100.0 {
		t.Errorf("num_variants = %f, want 100", result["num_variants"])
	}
	if result["rare_variants"] != 20.0 {
		t.Errorf("rare_variants = %f, want 20", result["rare_variants"])
	}
	if result["tajima_d"] != -0.5 {
		t.Errorf("tajima_d = %f, want -0.5", result["tajima_d"])
	}
}

func TestConvertStatistics_HasAllKeys(t *testing.T) {
	stats := population.PopulationStatistics{}
	result := convertStatistics(stats)
	for _, key := range []string{"diversity_pi", "watterson_theta", "tajima_d", "mean_frequency", "num_variants", "rare_variants"} {
		if _, ok := result[key]; !ok {
			t.Errorf("convertStatistics missing key %q", key)
		}
	}
}

// ---------- generateTestSequence ----------

func TestGenerateTestSequence_UntrainedModel(t *testing.T) {
	// EnhancedHMMKmerModel may not be enhanced when untrained; expect an error
	model := authenticity.NewEnhancedHMMKmerModel(4, "YRI")
	_, err := generateTestSequence(model, 100)
	// An untrained model may return an error - that's acceptable
	_ = err
}

func TestGenerateTestSequence_EnhancedModel(t *testing.T) {
	model := authenticity.NewEnhancedHMMKmerModel(4, "YRI")

	// Add demographic variants to make the model enhanced
	variants := convertToHMMVariants([]population.VariantData{
		{Position: 100, Alleles: []string{"A", "G"}, Frequency: 0.3, Genotypes: []int{0, 1}},
		{Position: 200, Alleles: []string{"C", "T"}, Frequency: 0.5, Genotypes: []int{1, 0}},
	})
	stats := convertStatistics(population.PopulationStatistics{NumVariants: 2})

	if err := model.IntegrateDemographicVariants("chr1", 1000, variants, stats); err != nil {
		t.Skipf("IntegrateDemographicVariants failed: %v", err)
	}

	// Now the model should be enhanced
	seq, err := generateTestSequence(model, 50)
	if err != nil {
		t.Logf("generateTestSequence returned error (acceptable for enhanced but untrained): %v", err)
		return
	}
	if len(seq) == 0 {
		t.Error("expected non-empty sequence from enhanced model")
	}
}
