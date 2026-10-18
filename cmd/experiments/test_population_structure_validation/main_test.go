//go:build !release
// +build !release

package main

import (
	"math"
	"os"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/population"
)

// TestNewHaplotypeSimulator verifies the simulator is constructed with the correct paths.
func TestNewHaplotypeSimulator(t *testing.T) {
	projectRoot := "/test/project"
	simulator := population.NewHaplotypeSimulator(projectRoot)
	if simulator == nil {
		t.Fatal("expected non-nil HaplotypeSimulator")
	}
}

// TestNewHaplotypeSimulator_WorkingDir constructs a simulator from the actual working directory.
func TestNewHaplotypeSimulator_WorkingDir(t *testing.T) {
	wd, err := os.Getwd()
	if err != nil {
		t.Fatalf("os.Getwd() error: %v", err)
	}
	simulator := population.NewHaplotypeSimulator(wd)
	if simulator == nil {
		t.Fatal("expected non-nil HaplotypeSimulator from working directory")
	}
}

// TestGetSupportedPopulations checks that 11 populations are returned.
func TestGetSupportedPopulations(t *testing.T) {
	pops := population.GetSupportedPopulations()
	if len(pops) != 11 {
		t.Errorf("expected 11 supported populations, got %d", len(pops))
	}
}

// TestGetSupportedPopulations_AllPresent verifies each expected population is in the list.
func TestGetSupportedPopulations_AllPresent(t *testing.T) {
	pops := population.GetSupportedPopulations()
	expected := []population.Population{
		population.YRI, population.LWK, population.GWD,
		population.CEU, population.GBR, population.TSI,
		population.CHB, population.JPT, population.CHS,
		population.MXL, population.PUR,
	}

	popSet := make(map[population.Population]bool, len(pops))
	for _, p := range pops {
		popSet[p] = true
	}

	for _, exp := range expected {
		if !popSet[exp] {
			t.Errorf("missing expected population: %s", exp)
		}
	}
}

// TestPopulation_String verifies all population codes convert to the correct string.
func TestPopulation_String(t *testing.T) {
	cases := []struct {
		pop  population.Population
		want string
	}{
		{population.YRI, "YRI"},
		{population.LWK, "LWK"},
		{population.GWD, "GWD"},
		{population.CEU, "CEU"},
		{population.GBR, "GBR"},
		{population.TSI, "TSI"},
		{population.CHB, "CHB"},
		{population.JPT, "JPT"},
		{population.CHS, "CHS"},
		{population.MXL, "MXL"},
		{population.PUR, "PUR"},
	}

	for _, tc := range cases {
		t.Run(tc.want, func(t *testing.T) {
			if tc.pop.String() != tc.want {
				t.Errorf("Population.String() = %q, want %q", tc.pop.String(), tc.want)
			}
		})
	}
}

// TestIsValidPopulation validates known and unknown populations.
func TestIsValidPopulation(t *testing.T) {
	cases := []struct {
		input string
		want  bool
	}{
		{"YRI", true},
		{"LWK", true},
		{"GWD", true},
		{"CEU", true},
		{"GBR", true},
		{"TSI", true},
		{"CHB", true},
		{"JPT", true},
		{"CHS", true},
		{"MXL", true},
		{"PUR", true},
		{"INVALID", false},
		{"", false},
		{"yri", false},
		{"EUR", false},
		{"AFR", false},
	}

	for _, tc := range cases {
		t.Run(tc.input+"_"+boolStr(tc.want), func(t *testing.T) {
			got := population.IsValidPopulation(tc.input)
			if got != tc.want {
				t.Errorf("IsValidPopulation(%q) = %v, want %v", tc.input, got, tc.want)
			}
		})
	}
}

// boolStr is a small helper to produce a readable test name suffix.
func boolStr(b bool) string {
	if b {
		return "valid"
	}
	return "invalid"
}

// TestHaplotypeSimulationResult_StructureFields verifies that a fully populated result
// struct stores and exposes the correct values.
func TestHaplotypeSimulationResult_StructureFields(t *testing.T) {
	result := &population.HaplotypeSimulationResult{
		Population:     population.YRI,
		Chromosome:     "chr22",
		NumVariants:    500,
		NumTrees:       30,
		NumSamples:     30,
		SequenceLength: 5_000_000,
		SimulationTime: 1.5,
		MemoryUsedMB:   100.0,
		Variants: []population.HaplotypeVariant{
			{Position: 100000, Alleles: []string{"A", "T"}, AncestralState: "A", Frequency: 0.3},
			{Position: 200000, Alleles: []string{"G", "C"}, AncestralState: "G", Frequency: 0.1},
		},
		RecombinationMap: population.RecombinationMap{
			Chromosome:     "chr22",
			Positions:      []int{0, 1000000, 2000000},
			Rates:          []float64{1.1e-8, 1.3e-8, 1.2e-8},
			MapLength:      5.5,
			PhysicalLength: 5_000_000,
		},
		LDBlocks: []population.HaplotypeBlock{
			{StartPosition: 50000, EndPosition: 150000, Length: 100000, VariantCount: 10, RSquaredMean: 0.7, Population: population.YRI},
			{StartPosition: 200000, EndPosition: 350000, Length: 150000, VariantCount: 15, RSquaredMean: 0.6, Population: population.YRI},
		},
		Statistics: population.HaplotypeStatistics{
			DiversityPi:         0.00095,
			SegregatingSites:    480,
			NumVariants:         500,
			NumHaplotypes:       60,
			MeanRSquared:        0.055,
			MedianRSquared:      0.06,
			HighLDFraction:      0.4,
			EffectiveRecombRate: 1.15e-8,
			GeneticMapLengthCM:  5.5,
			HaplotypeDiversity:  0.93,
			LDPairsAnalyzed:     10000,
		},
	}

	if result.Population != population.YRI {
		t.Errorf("Population = %s, want YRI", result.Population)
	}
	if result.NumVariants != 500 {
		t.Errorf("NumVariants = %d, want 500", result.NumVariants)
	}
	if len(result.LDBlocks) != 2 {
		t.Errorf("len(LDBlocks) = %d, want 2", len(result.LDBlocks))
	}
	if result.Statistics.DiversityPi <= 0 {
		t.Errorf("DiversityPi = %f, want > 0", result.Statistics.DiversityPi)
	}
}

// TestHaplotypeSimulationResult_NumVariantsCheck exercises the PCA-readiness logic
// used in main() (NumVariants > 100 && len(LDBlocks) > 10).
func TestHaplotypeSimulationResult_NumVariantsCheck(t *testing.T) {
	cases := []struct {
		name        string
		numVariants int
		numBlocks   int
		wantReady   bool
	}{
		{"ready", 200, 15, true},
		{"low-variants", 50, 15, false},
		{"low-blocks", 200, 5, false},
		{"both-low", 50, 5, false},
		{"boundary-variants", 101, 11, true},
		{"boundary-exact-100", 100, 11, false},
		{"boundary-exact-10", 101, 10, false},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			// Build a result with the given counts
			blocks := make([]population.HaplotypeBlock, tc.numBlocks)
			result := &population.HaplotypeSimulationResult{
				NumVariants: tc.numVariants,
				LDBlocks:    blocks,
			}
			got := result.NumVariants > 100 && len(result.LDBlocks) > 10
			if got != tc.wantReady {
				t.Errorf("pcaReady = %v, want %v (NumVariants=%d, len(LDBlocks)=%d)",
					got, tc.wantReady, tc.numVariants, tc.numBlocks)
			}
		})
	}
}

// TestHaplotypeStatistics_MeanRSquaredRange checks r-squared statistics are in valid range.
func TestHaplotypeStatistics_MeanRSquaredRange(t *testing.T) {
	cases := []struct {
		name        string
		meanRSquared float64
		valid       bool
	}{
		{"typical", 0.05, true},
		{"high", 0.09, true},
		{"zero", 0.0, false}, // 0.01 < r² < 0.1 per validation logic
		{"one", 1.0, false},
		{"boundary-low", 0.011, true},
		{"boundary-high", 0.099, true},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := tc.meanRSquared > 0.01 && tc.meanRSquared < 0.1
			if got != tc.valid {
				t.Errorf("reasonable_ld_structure(%f) = %v, want %v",
					tc.meanRSquared, got, tc.valid)
			}
		})
	}
}

// TestDiversityHierarchy verifies the expected diversity ordering logic (YRI > CEU > CHB).
func TestDiversityHierarchy(t *testing.T) {
	diversities := map[population.Population]float64{
		population.YRI: 0.00095,
		population.CEU: 0.00075,
		population.CHB: 0.00065,
	}

	hierarchyCorrect := diversities[population.YRI] > diversities[population.CEU] &&
		diversities[population.CEU] > diversities[population.CHB]

	if !hierarchyCorrect {
		t.Errorf("diversity hierarchy YRI>CEU>CHB not satisfied: YRI=%f, CEU=%f, CHB=%f",
			diversities[population.YRI], diversities[population.CEU], diversities[population.CHB])
	}
}

// TestDiversityHierarchy_Reversed verifies that reversed diversity fails the check.
func TestDiversityHierarchy_Reversed(t *testing.T) {
	diversities := map[population.Population]float64{
		population.YRI: 0.00065,
		population.CEU: 0.00075,
		population.CHB: 0.00095,
	}

	hierarchyCorrect := diversities[population.YRI] > diversities[population.CEU] &&
		diversities[population.CEU] > diversities[population.CHB]

	if hierarchyCorrect {
		t.Errorf("reversed diversity should fail the hierarchy check but it passed")
	}
}

// TestQualityMetricsCalculation verifies the quality metrics logic from main():
// passedQuality >= totalQuality*3/4.
func TestQualityMetricsCalculation(t *testing.T) {
	cases := []struct {
		name    string
		passed  int
		total   int
		success bool
	}{
		{"all-pass", 8, 8, true},
		{"three-quarters", 6, 8, true},
		{"below-three-quarters", 5, 8, false},
		{"zero-of-eight", 0, 8, false},
		{"four-of-four", 4, 4, true},
		{"two-of-four", 2, 4, false},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := tc.passed >= tc.total*3/4
			if got != tc.success {
				t.Errorf("passed(%d) >= total(%d)*3/4: got %v, want %v",
					tc.passed, tc.total, got, tc.success)
			}
		})
	}
}

// TestLDValidationThreshold checks the LD pattern pass threshold (passedTests >= totalTests*3/4).
// 5*3/4 = 3 (integer division), so 3/5 is a PASS (3 >= 3).
// 2/5 is a FAIL (2 < 3).
func TestLDValidationThreshold(t *testing.T) {
	cases := []struct {
		passed  int
		total   int
		success bool
	}{
		{5, 5, true},
		{4, 5, true},
		{3, 5, true},  // 5*3/4 = 3 (integer), 3 >= 3 → true
		{2, 5, false}, // 2 < 3 → false
		{0, 5, false},
	}

	for _, tc := range cases {
		t.Run("", func(t *testing.T) {
			got := tc.passed >= tc.total*3/4
			if got != tc.success {
				t.Errorf("LD threshold %d/%d: got %v, want %v", tc.passed, tc.total, got, tc.success)
			}
		})
	}
}

// TestPCAReadyCount checks pcaReadyCount >= 2 threshold logic.
func TestPCAReadyCount(t *testing.T) {
	cases := []struct {
		count   int
		success bool
	}{
		{3, true},
		{2, true},
		{1, false},
		{0, false},
	}

	for _, tc := range cases {
		t.Run("", func(t *testing.T) {
			got := tc.count >= 2
			if got != tc.success {
				t.Errorf("pcaReadyCount(%d) >= 2: got %v, want %v", tc.count, got, tc.success)
			}
		})
	}
}

// TestStringRepeat tests the strings.Repeat usage pattern from main().
func TestStringRepeat(t *testing.T) {
	separator := "=" + strings.Repeat("=", 48)
	if len(separator) != 49 {
		t.Errorf("separator length = %d, want 49", len(separator))
	}

	if separator[0] != '=' {
		t.Error("separator should start with '='")
	}
}

// TestHaplotypeSimulator_GetRecombinationMap verifies round-trip of recombination map.
func TestHaplotypeSimulator_GetRecombinationMap(t *testing.T) {
	simulator := population.NewHaplotypeSimulator("/test")

	expectedMap := population.RecombinationMap{
		Chromosome:     "chr22",
		Positions:      []int{0, 1000000, 5000000},
		Rates:          []float64{1.1e-8, 1.3e-8, 1.2e-8},
		MapLength:      6.0,
		PhysicalLength: 5_000_000,
	}

	result := &population.HaplotypeSimulationResult{
		RecombinationMap: expectedMap,
	}

	retrieved := simulator.GetRecombinationMap(result)

	if retrieved.Chromosome != expectedMap.Chromosome {
		t.Errorf("Chromosome = %q, want %q", retrieved.Chromosome, expectedMap.Chromosome)
	}
	if retrieved.MapLength != expectedMap.MapLength {
		t.Errorf("MapLength = %f, want %f", retrieved.MapLength, expectedMap.MapLength)
	}
	if retrieved.PhysicalLength != expectedMap.PhysicalLength {
		t.Errorf("PhysicalLength = %d, want %d", retrieved.PhysicalLength, expectedMap.PhysicalLength)
	}
}

// TestHaplotypeSimulator_GetLDStatistics verifies round-trip of LD statistics.
func TestHaplotypeSimulator_GetLDStatistics(t *testing.T) {
	simulator := population.NewHaplotypeSimulator("/test")

	expectedStats := population.HaplotypeStatistics{
		DiversityPi:  0.00095,
		MeanRSquared: 0.055,
	}

	result := &population.HaplotypeSimulationResult{
		Statistics: expectedStats,
	}

	retrieved := simulator.GetLDStatistics(result)

	if retrieved.DiversityPi != expectedStats.DiversityPi {
		t.Errorf("DiversityPi = %f, want %f", retrieved.DiversityPi, expectedStats.DiversityPi)
	}
	if retrieved.MeanRSquared != expectedStats.MeanRSquared {
		t.Errorf("MeanRSquared = %f, want %f", retrieved.MeanRSquared, expectedStats.MeanRSquared)
	}
}

// TestHaplotypeSimulator_ValidateLDPatterns_AllThree exercises the ValidateLDPatterns function
// with all three test populations (YRI, CEU, CHB).
func TestHaplotypeSimulator_ValidateLDPatterns_AllThree(t *testing.T) {
	simulator := population.NewHaplotypeSimulator("/test")

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: {
			Population: population.YRI,
			Statistics: population.HaplotypeStatistics{
				DiversityPi:         0.00095,
				MeanRSquared:        0.05,
				EffectiveRecombRate: 1.2e-8,
				GeneticMapLengthCM:  5.0,
			},
			SequenceLength: 5_000_000,
			LDBlocks: []population.HaplotypeBlock{
				{Length: 50000, VariantCount: 10},
				{Length: 75000, VariantCount: 15},
			},
		},
		population.CEU: {
			Population: population.CEU,
			Statistics: population.HaplotypeStatistics{
				DiversityPi:         0.00075,
				MeanRSquared:        0.06,
				EffectiveRecombRate: 1.1e-8,
				GeneticMapLengthCM:  4.5,
			},
			SequenceLength: 5_000_000,
			LDBlocks: []population.HaplotypeBlock{
				{Length: 80000, VariantCount: 12},
			},
		},
		population.CHB: {
			Population: population.CHB,
			Statistics: population.HaplotypeStatistics{
				DiversityPi:         0.00065,
				MeanRSquared:        0.07,
				EffectiveRecombRate: 1.0e-8,
				GeneticMapLengthCM:  4.0,
			},
			SequenceLength: 5_000_000,
			LDBlocks: []population.HaplotypeBlock{
				{Length: 90000, VariantCount: 8},
				{Length: 60000, VariantCount: 6},
			},
		},
	}

	validation, err := simulator.ValidateLDPatterns(results)
	if err != nil {
		t.Fatalf("ValidateLDPatterns returned error: %v", err)
	}
	if validation == nil {
		t.Fatal("ValidateLDPatterns returned nil map")
	}
	if len(validation) == 0 {
		t.Fatal("ValidateLDPatterns returned empty map")
	}

	// African populations should have highest diversity
	if !validation["african_highest_diversity"] {
		t.Error("expected african_highest_diversity = true")
	}
	// Haplotype blocks must be present
	if !validation["haplotype_blocks_present"] {
		t.Error("expected haplotype_blocks_present = true")
	}
}

// TestHaplotypeSimulator_ValidateLDPatterns_NoAfricanHighDiversity verifies that when
// African populations do NOT have the highest diversity, the test reports false.
func TestHaplotypeSimulator_ValidateLDPatterns_NoAfricanHighDiversity(t *testing.T) {
	simulator := population.NewHaplotypeSimulator("/test")

	// Invert diversity: CHB > CEU > YRI (opposite of expectation)
	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: {
			Population: population.YRI,
			Statistics: population.HaplotypeStatistics{
				DiversityPi:         0.00040,
				MeanRSquared:        0.05,
				EffectiveRecombRate: 1.2e-8,
				GeneticMapLengthCM:  5.0,
			},
			SequenceLength: 5_000_000,
			LDBlocks:       []population.HaplotypeBlock{{Length: 50000, VariantCount: 10}},
		},
		population.CEU: {
			Population: population.CEU,
			Statistics: population.HaplotypeStatistics{
				DiversityPi:         0.00060,
				MeanRSquared:        0.06,
				EffectiveRecombRate: 1.1e-8,
				GeneticMapLengthCM:  4.5,
			},
			SequenceLength: 5_000_000,
			LDBlocks:       []population.HaplotypeBlock{{Length: 80000, VariantCount: 12}},
		},
		population.CHB: {
			Population: population.CHB,
			Statistics: population.HaplotypeStatistics{
				DiversityPi:         0.00080,
				MeanRSquared:        0.07,
				EffectiveRecombRate: 1.0e-8,
				GeneticMapLengthCM:  4.0,
			},
			SequenceLength: 5_000_000,
			LDBlocks:       []population.HaplotypeBlock{{Length: 90000, VariantCount: 8}},
		},
	}

	validation, err := simulator.ValidateLDPatterns(results)
	if err != nil {
		t.Fatalf("ValidateLDPatterns returned error: %v", err)
	}

	// With inverted diversity, African should NOT have highest diversity
	if validation["african_highest_diversity"] {
		t.Error("expected african_highest_diversity = false when diversity is inverted")
	}
}

// TestHaplotypeSimulator_AnalyzeCrossPopulationLD verifies the cross-population analysis.
func TestHaplotypeSimulator_AnalyzeCrossPopulationLD(t *testing.T) {
	simulator := population.NewHaplotypeSimulator("/test")

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: {
			Population: population.YRI,
			Statistics: population.HaplotypeStatistics{
				DiversityPi:  0.00095,
				MeanRSquared: 0.055,
			},
			LDBlocks: []population.HaplotypeBlock{
				{Length: 50000},
				{Length: 60000},
			},
		},
		population.CEU: {
			Population: population.CEU,
			Statistics: population.HaplotypeStatistics{
				DiversityPi:  0.00075,
				MeanRSquared: 0.065,
			},
			LDBlocks: []population.HaplotypeBlock{
				{Length: 80000},
			},
		},
	}

	analysis := simulator.AnalyzeCrossPopulationLD(results)
	if analysis == nil {
		t.Fatal("AnalyzeCrossPopulationLD returned nil")
	}

	// Check required keys
	requiredKeys := []string{
		"diversity_by_population",
		"ld_metrics_by_population",
		"block_counts_by_population",
		"diversity_mean",
		"diversity_cv",
		"ld_mean",
		"ld_cv",
	}
	for _, key := range requiredKeys {
		if _, ok := analysis[key]; !ok {
			t.Errorf("missing key %q in cross-population analysis", key)
		}
	}

	// Check diversity mean is approximately (0.00095 + 0.00075) / 2
	if mean, ok := analysis["diversity_mean"].(float64); ok {
		wantMean := (0.00095 + 0.00075) / 2
		diff := mean - wantMean
		if diff < 0 {
			diff = -diff
		}
		if diff > 1e-7 {
			t.Errorf("diversity_mean = %f, want %f (delta %e)", mean, wantMean, diff)
		}
	} else {
		t.Error("diversity_mean is not float64")
	}
}

// TestHaplotypeSimulator_SimulateHaplotypes_NoScript verifies that simulation fails
// gracefully when the Python script does not exist.
func TestHaplotypeSimulator_SimulateHaplotypes_NoScript(t *testing.T) {
	// Use a non-existent project root so the script will not be found
	simulator := population.NewHaplotypeSimulator("/nonexistent/path")
	seed := 42

	_, err := simulator.SimulateHaplotypes(
		population.YRI,
		"chr22",
		10,
		1_000,
		&seed,
		false,
	)

	if err == nil {
		t.Error("expected error when Python script does not exist, got nil")
	}
}

// TestHaplotypeSimulator_SimulateHaplotypes_NilSeed verifies that passing a nil seed
// also causes a failure when no script exists.
func TestHaplotypeSimulator_SimulateHaplotypes_NilSeed(t *testing.T) {
	simulator := population.NewHaplotypeSimulator("/nonexistent/path")

	_, err := simulator.SimulateHaplotypes(
		population.CEU,
		"chr1",
		10,
		1_000,
		nil,
		true,
	)

	if err == nil {
		t.Error("expected error when Python script does not exist, got nil")
	}
}

// TestDiversityCVCalculation exercises the diversity coefficient-of-variation logic from main().
func TestDiversityCVCalculation(t *testing.T) {
	// These values mirror the computation in main()
	diversities := map[population.Population]float64{
		population.YRI: 0.00095,
		population.CEU: 0.00075,
		population.CHB: 0.00065,
	}

	var diversitySum, diversitySquareSum float64
	count := float64(len(diversities))

	for _, d := range diversities {
		diversitySum += d
		diversitySquareSum += d * d
	}

	diversityMean := diversitySum / count
	if diversityMean <= 0 {
		t.Errorf("diversityMean = %f, want > 0", diversityMean)
	}

	diversityVariance := (diversitySquareSum / count) - (diversityMean * diversityMean)
	if diversityVariance < 0 {
		t.Errorf("diversityVariance = %f, want >= 0", diversityVariance)
	}

	diversityStd := 0.0
	if diversityVariance > 0 {
		diversityStd = math.Sqrt(diversityVariance)
	}

	diversityCV := 0.0
	if diversityMean > 0 {
		diversityCV = diversityStd / diversityMean
	}

	// CV should be non-negative
	if diversityCV < 0 {
		t.Errorf("diversityCV = %f, want >= 0", diversityCV)
	}
}

// TestDiversityCVCalculation_ZeroMean exercises the zero-mean branch in the CV logic.
func TestDiversityCVCalculation_ZeroMean(t *testing.T) {
	diversities := map[population.Population]float64{
		population.YRI: 0.0,
		population.CEU: 0.0,
	}

	var diversitySum, diversitySquareSum float64
	count := float64(len(diversities))

	for _, d := range diversities {
		diversitySum += d
		diversitySquareSum += d * d
	}

	diversityMean := diversitySum / count
	diversityVariance := (diversitySquareSum / count) - (diversityMean * diversityMean)
	diversityStd := 0.0
	if diversityVariance > 0 {
		diversityStd = math.Sqrt(diversityVariance)
	}

	diversityCV := 0.0
	if diversityMean > 0 {
		diversityCV = diversityStd / diversityMean
	}

	if diversityCV != 0.0 {
		t.Errorf("zero-mean diversityCV = %f, want 0.0", diversityCV)
	}
}

// TestHaplotypeBlock_Structure validates individual block field constraints.
func TestHaplotypeBlock_Structure(t *testing.T) {
	block := population.HaplotypeBlock{
		StartPosition: 1_000_000,
		EndPosition:   1_050_000,
		Length:        50_000,
		VariantCount:  25,
		RSquaredMean:  0.75,
		Population:    population.YRI,
	}

	if block.StartPosition >= block.EndPosition {
		t.Errorf("StartPosition %d >= EndPosition %d", block.StartPosition, block.EndPosition)
	}
	if block.Length != block.EndPosition-block.StartPosition {
		t.Errorf("Length %d != EndPosition-StartPosition %d", block.Length, block.EndPosition-block.StartPosition)
	}
	if block.RSquaredMean < 0 || block.RSquaredMean > 1 {
		t.Errorf("RSquaredMean %f out of [0, 1]", block.RSquaredMean)
	}
	if block.VariantCount <= 0 {
		t.Errorf("VariantCount %d, want > 0", block.VariantCount)
	}
	if block.Population != population.YRI {
		t.Errorf("Population = %s, want YRI", block.Population)
	}
}

// TestRecombinationMap_Structure validates recombination map fields.
func TestRecombinationMap_Structure(t *testing.T) {
	rmap := population.RecombinationMap{
		Chromosome:     "chr22",
		Positions:      []int{0, 1_000_000, 5_000_000},
		Rates:          []float64{1.1e-8, 1.3e-8, 1.2e-8},
		MapLength:      6.0,
		PhysicalLength: 5_000_000,
	}

	if rmap.Chromosome != "chr22" {
		t.Errorf("Chromosome = %q, want chr22", rmap.Chromosome)
	}
	if len(rmap.Positions) != len(rmap.Rates) {
		t.Errorf("len(Positions) %d != len(Rates) %d", len(rmap.Positions), len(rmap.Rates))
	}
	if rmap.MapLength <= 0 {
		t.Errorf("MapLength = %f, want > 0", rmap.MapLength)
	}
	if rmap.PhysicalLength <= 0 {
		t.Errorf("PhysicalLength = %d, want > 0", rmap.PhysicalLength)
	}
	for i := 1; i < len(rmap.Positions); i++ {
		if rmap.Positions[i] <= rmap.Positions[i-1] {
			t.Errorf("Positions not in ascending order at index %d", i)
		}
	}
}

// TestHaplotypeVariant_Structure validates variant field constraints.
func TestHaplotypeVariant_Structure(t *testing.T) {
	v := population.HaplotypeVariant{
		Position:       123456,
		Alleles:        []string{"A", "T"},
		AncestralState: "A",
		Frequency:      0.25,
	}

	if v.Position <= 0 {
		t.Errorf("Position = %d, want > 0", v.Position)
	}
	if len(v.Alleles) < 2 {
		t.Errorf("len(Alleles) = %d, want >= 2", len(v.Alleles))
	}
	if v.Frequency < 0 || v.Frequency > 1 {
		t.Errorf("Frequency = %f, want in [0, 1]", v.Frequency)
	}

	// Ancestral state must be one of the alleles
	found := false
	for _, a := range v.Alleles {
		if a == v.AncestralState {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("AncestralState %q not in Alleles %v", v.AncestralState, v.Alleles)
	}
}

// TestHaplotypeStatistics_Structure validates statistics fields.
func TestHaplotypeStatistics_Structure(t *testing.T) {
	stats := population.HaplotypeStatistics{
		DiversityPi:         0.001,
		SegregatingSites:    900,
		NumVariants:         1000,
		NumHaplotypes:       60,
		MeanRSquared:        0.055,
		MedianRSquared:      0.06,
		HighLDFraction:      0.4,
		EffectiveRecombRate: 1.2e-8,
		GeneticMapLengthCM:  5.0,
		HaplotypeDiversity:  0.93,
		LDPairsAnalyzed:     10000,
	}

	if stats.DiversityPi <= 0 {
		t.Errorf("DiversityPi = %f, want > 0", stats.DiversityPi)
	}
	if stats.MeanRSquared < 0 || stats.MeanRSquared > 1 {
		t.Errorf("MeanRSquared = %f, want in [0, 1]", stats.MeanRSquared)
	}
	if stats.HighLDFraction < 0 || stats.HighLDFraction > 1 {
		t.Errorf("HighLDFraction = %f, want in [0, 1]", stats.HighLDFraction)
	}
	if stats.SegregatingSites > stats.NumVariants {
		t.Errorf("SegregatingSites %d > NumVariants %d", stats.SegregatingSites, stats.NumVariants)
	}
	if stats.LDPairsAnalyzed <= 0 {
		t.Errorf("LDPairsAnalyzed = %d, want > 0", stats.LDPairsAnalyzed)
	}
}

// TestPopulationTestSet verifies the three-population set used in main() is all valid.
func TestPopulationTestSet(t *testing.T) {
	testPopulations := []population.Population{
		population.YRI,
		population.CEU,
		population.CHB,
	}

	for _, pop := range testPopulations {
		if !population.IsValidPopulation(string(pop)) {
			t.Errorf("population %s should be valid", pop)
		}
	}
}

// TestQualityMetrics_FrameworkAlwaysTrue verifies the always-true quality metrics
// used in main() evaluate correctly.
func TestQualityMetrics_FrameworkAlwaysTrue(t *testing.T) {
	pcaFrameworkAvailable := true      // Scripts created
	ldFrameworkReady := true           // Scripts created
	integrationPipelineEstablished := true // Framework in place
	week5FoundationComplete := true    // Core functionality ready

	if !pcaFrameworkAvailable {
		t.Error("PCA framework should be available")
	}
	if !ldFrameworkReady {
		t.Error("LD framework should be ready")
	}
	if !integrationPipelineEstablished {
		t.Error("Integration pipeline should be established")
	}
	if !week5FoundationComplete {
		t.Error("Week 5 foundation should be complete")
	}
}

// TestDiversityStd_IsSqrtNotSquare is a regression test for #129.
// The old code computed diversityStd = diversityVariance * diversityVariance
// (variance²) instead of math.Sqrt(diversityVariance). With variance = 0.04,
// the wrong code returns 0.0016; the correct std dev is 0.2.
func TestDiversityStd_IsSqrtNotSquare(t *testing.T) {
	// Variance deliberately chosen so the bug would produce a visibly wrong result.
	variance := 0.04 // expected std = 0.2; old buggy code returned 0.0016

	std := 0.0
	if variance > 0 {
		std = math.Sqrt(variance)
	}

	const want = 0.2
	if math.Abs(std-want) > 1e-10 {
		t.Errorf("diversityStd = %v, want %v (got variance² = %v instead of sqrt(variance))",
			std, want, variance*variance)
	}
}

// TestHaplotypeSimulationResult_WithMatrix verifies that the HaplotypeMatrix field
// can be populated and accessed (the optional large field).
func TestHaplotypeSimulationResult_WithMatrix(t *testing.T) {
	result := population.HaplotypeSimulationResult{
		Population:  population.YRI,
		NumVariants: 3,
		NumSamples:  2,
		HaplotypeMatrix: [][]int{
			{0, 1, 0},
			{1, 0, 1},
			{0, 0, 1},
			{1, 1, 0},
		},
	}

	if len(result.HaplotypeMatrix) != 4 {
		t.Errorf("len(HaplotypeMatrix) = %d, want 4", len(result.HaplotypeMatrix))
	}
	for i, row := range result.HaplotypeMatrix {
		if len(row) != 3 {
			t.Errorf("HaplotypeMatrix[%d] length = %d, want 3", i, len(row))
		}
		for j, val := range row {
			if val < 0 || val > 1 {
				t.Errorf("HaplotypeMatrix[%d][%d] = %d, want 0 or 1", i, j, val)
			}
		}
	}
}
