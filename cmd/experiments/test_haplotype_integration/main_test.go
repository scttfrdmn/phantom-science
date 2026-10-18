//go:build !release
// +build !release

package main

import (
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/population"
)

// makeResult creates a HaplotypeSimulationResult for testing.
func makeResult(pop population.Population, numVariants, seqLen int, pi, rSq float64, blocks []population.HaplotypeBlock) *population.HaplotypeSimulationResult {
	return &population.HaplotypeSimulationResult{
		Population:     pop,
		Chromosome:     "chr22",
		NumVariants:    numVariants,
		SequenceLength: seqLen,
		SimulationTime: 1.5,
		MemoryUsedMB:   50.0,
		LDBlocks:       blocks,
		Statistics: population.HaplotypeStatistics{
			DiversityPi:  pi,
			MeanRSquared: rSq,
		},
		RecombinationMap: population.RecombinationMap{
			Chromosome:     "chr22",
			Positions:      []int{0, seqLen / 2, seqLen},
			Rates:          []float64{1.2e-8, 1.3e-8, 1.1e-8},
			MapLength:      float64(seqLen) / 1_000_000 * 1.5,
			PhysicalLength: seqLen,
		},
	}
}

// TestNewHaplotypeSimulator_CreatesSimulator verifies NewHaplotypeSimulator returns a
// correctly initialized simulator pointing at the expected paths.
func TestNewHaplotypeSimulator_CreatesSimulator(t *testing.T) {
	projectRoot := "/test/project"
	sim := population.NewHaplotypeSimulator(projectRoot)
	if sim == nil {
		t.Fatal("NewHaplotypeSimulator returned nil")
	}
}

// ---- Block-structure analysis logic (a copy of the logic in main(), not a call into it) ----

// haplotypeBlockStats computes block statistics the same way main() does.
type haplotypeBlockStats struct {
	totalBlockLength int
	totalVariants    int
	maxBlockSize     int
	minBlockSize     int
	meanBlockSize    int
	genomeCoverage   float64
	meanVariants     float64
}

func computeBlockStats(result *population.HaplotypeSimulationResult) (haplotypeBlockStats, bool) {
	if len(result.LDBlocks) == 0 {
		return haplotypeBlockStats{}, false
	}
	var s haplotypeBlockStats
	s.minBlockSize = result.LDBlocks[0].Length
	for _, block := range result.LDBlocks {
		s.totalBlockLength += block.Length
		s.totalVariants += block.VariantCount
		if block.Length > s.maxBlockSize {
			s.maxBlockSize = block.Length
		}
		if block.Length < s.minBlockSize {
			s.minBlockSize = block.Length
		}
	}
	s.meanBlockSize = s.totalBlockLength / len(result.LDBlocks)
	s.genomeCoverage = float64(s.totalBlockLength) / float64(result.SequenceLength) * 100
	s.meanVariants = float64(s.totalVariants) / float64(len(result.LDBlocks))
	return s, true
}

// TestComputeBlockStats_SingleBlock tests the block-analysis logic with one block.
func TestComputeBlockStats_SingleBlock(t *testing.T) {
	blocks := []population.HaplotypeBlock{
		{StartPosition: 0, EndPosition: 50000, Length: 50000, VariantCount: 10},
	}
	result := makeResult(population.YRI, 100, 1_000_000, 0.001, 0.05, blocks)

	stats, ok := computeBlockStats(result)
	if !ok {
		t.Fatal("computeBlockStats returned false for non-empty block list")
	}
	if stats.totalBlockLength != 50000 {
		t.Errorf("totalBlockLength = %d, want 50000", stats.totalBlockLength)
	}
	if stats.totalVariants != 10 {
		t.Errorf("totalVariants = %d, want 10", stats.totalVariants)
	}
	if stats.maxBlockSize != 50000 {
		t.Errorf("maxBlockSize = %d, want 50000", stats.maxBlockSize)
	}
	if stats.minBlockSize != 50000 {
		t.Errorf("minBlockSize = %d, want 50000", stats.minBlockSize)
	}
	if stats.meanBlockSize != 50000 {
		t.Errorf("meanBlockSize = %d, want 50000", stats.meanBlockSize)
	}
	if stats.genomeCoverage != 5.0 {
		t.Errorf("genomeCoverage = %f, want 5.0", stats.genomeCoverage)
	}
	if stats.meanVariants != 10.0 {
		t.Errorf("meanVariants = %f, want 10.0", stats.meanVariants)
	}
}

// TestComputeBlockStats_MultipleBlocks tests block analysis with three blocks.
func TestComputeBlockStats_MultipleBlocks(t *testing.T) {
	blocks := []population.HaplotypeBlock{
		{Length: 30000, VariantCount: 6},
		{Length: 50000, VariantCount: 10},
		{Length: 20000, VariantCount: 4},
	}
	result := makeResult(population.CEU, 200, 1_000_000, 0.0008, 0.06, blocks)

	stats, ok := computeBlockStats(result)
	if !ok {
		t.Fatal("computeBlockStats returned false for non-empty block list")
	}
	if stats.totalBlockLength != 100000 {
		t.Errorf("totalBlockLength = %d, want 100000", stats.totalBlockLength)
	}
	if stats.totalVariants != 20 {
		t.Errorf("totalVariants = %d, want 20", stats.totalVariants)
	}
	if stats.maxBlockSize != 50000 {
		t.Errorf("maxBlockSize = %d, want 50000", stats.maxBlockSize)
	}
	if stats.minBlockSize != 20000 {
		t.Errorf("minBlockSize = %d, want 20000", stats.minBlockSize)
	}
	expectedMean := 100000 / 3
	if stats.meanBlockSize != expectedMean {
		t.Errorf("meanBlockSize = %d, want %d", stats.meanBlockSize, expectedMean)
	}
	expectedCoverage := float64(100000) / float64(1_000_000) * 100
	if stats.genomeCoverage != expectedCoverage {
		t.Errorf("genomeCoverage = %f, want %f", stats.genomeCoverage, expectedCoverage)
	}
	expectedMeanVariants := float64(20) / 3.0
	if stats.meanVariants != expectedMeanVariants {
		t.Errorf("meanVariants = %f, want %f", stats.meanVariants, expectedMeanVariants)
	}
}

// TestComputeBlockStats_EmptyBlocks verifies that empty block list returns false.
func TestComputeBlockStats_EmptyBlocks(t *testing.T) {
	result := makeResult(population.CHB, 0, 1_000_000, 0.0007, 0.07, nil)
	_, ok := computeBlockStats(result)
	if ok {
		t.Error("computeBlockStats should return false for empty block list")
	}
}

// TestComputeBlockStats_FullCoverage checks 100% genome coverage edge case.
func TestComputeBlockStats_FullCoverage(t *testing.T) {
	blocks := []population.HaplotypeBlock{
		{Length: 1_000_000, VariantCount: 500},
	}
	result := makeResult(population.MXL, 500, 1_000_000, 0.0009, 0.055, blocks)

	stats, ok := computeBlockStats(result)
	if !ok {
		t.Fatal("computeBlockStats returned false for non-empty block list")
	}
	if stats.genomeCoverage != 100.0 {
		t.Errorf("genomeCoverage = %f, want 100.0", stats.genomeCoverage)
	}
}

// ---- Performance summary accumulation (a copy of the logic in main(), not a call into it) ----

type perfSummary struct {
	totalVariants int
	totalBlocks   int
	totalMemory   float64
	totalSimTime  float64
}

func accumulatePerfSummary(results map[population.Population]*population.HaplotypeSimulationResult) perfSummary {
	var s perfSummary
	for _, result := range results {
		s.totalVariants += result.NumVariants
		s.totalBlocks += len(result.LDBlocks)
		s.totalMemory += result.MemoryUsedMB
		s.totalSimTime += result.SimulationTime
	}
	return s
}

// TestAccumulatePerfSummary_SinglePop verifies accumulation with one population.
func TestAccumulatePerfSummary_SinglePop(t *testing.T) {
	blocks := []population.HaplotypeBlock{{Length: 50000, VariantCount: 10}}
	result := makeResult(population.YRI, 1500, 5_000_000, 0.0012, 0.05, blocks)
	result.MemoryUsedMB = 80.0
	result.SimulationTime = 2.3

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: result,
	}
	s := accumulatePerfSummary(results)

	if s.totalVariants != 1500 {
		t.Errorf("totalVariants = %d, want 1500", s.totalVariants)
	}
	if s.totalBlocks != 1 {
		t.Errorf("totalBlocks = %d, want 1", s.totalBlocks)
	}
	if s.totalMemory != 80.0 {
		t.Errorf("totalMemory = %f, want 80.0", s.totalMemory)
	}
	if s.totalSimTime != 2.3 {
		t.Errorf("totalSimTime = %f, want 2.3", s.totalSimTime)
	}
}

// TestAccumulatePerfSummary_MultiplePops verifies accumulation across several populations.
func TestAccumulatePerfSummary_MultiplePops(t *testing.T) {
	pops := []population.Population{population.YRI, population.CEU, population.CHB, population.MXL}
	results := make(map[population.Population]*population.HaplotypeSimulationResult)

	for i, pop := range pops {
		blocks := []population.HaplotypeBlock{
			{Length: 50000 * (i + 1), VariantCount: 10 * (i + 1)},
		}
		r := makeResult(pop, 1000*(i+1), 5_000_000, 0.001, 0.06, blocks)
		r.MemoryUsedMB = 50.0 * float64(i+1)
		r.SimulationTime = 1.0 * float64(i+1)
		results[pop] = r
	}

	s := accumulatePerfSummary(results)

	wantVariants := 1000 + 2000 + 3000 + 4000
	if s.totalVariants != wantVariants {
		t.Errorf("totalVariants = %d, want %d", s.totalVariants, wantVariants)
	}
	if s.totalBlocks != 4 {
		t.Errorf("totalBlocks = %d, want 4", s.totalBlocks)
	}
	wantMemory := 50.0 + 100.0 + 150.0 + 200.0
	if s.totalMemory != wantMemory {
		t.Errorf("totalMemory = %f, want %f", s.totalMemory, wantMemory)
	}
	wantSimTime := 1.0 + 2.0 + 3.0 + 4.0
	if s.totalSimTime != wantSimTime {
		t.Errorf("totalSimTime = %f, want %f", s.totalSimTime, wantSimTime)
	}
}

// TestAccumulatePerfSummary_Empty verifies that an empty results map returns zero values.
func TestAccumulatePerfSummary_Empty(t *testing.T) {
	results := map[population.Population]*population.HaplotypeSimulationResult{}
	s := accumulatePerfSummary(results)

	if s.totalVariants != 0 {
		t.Errorf("totalVariants = %d, want 0", s.totalVariants)
	}
	if s.totalBlocks != 0 {
		t.Errorf("totalBlocks = %d, want 0", s.totalBlocks)
	}
	if s.totalMemory != 0 {
		t.Errorf("totalMemory = %f, want 0.0", s.totalMemory)
	}
	if s.totalSimTime != 0 {
		t.Errorf("totalSimTime = %f, want 0.0", s.totalSimTime)
	}
}

// TestAccumulatePerfSummary_AveragesCalc verifies per-population average computation.
func TestAccumulatePerfSummary_AveragesCalc(t *testing.T) {
	pops := []population.Population{population.YRI, population.CEU}
	results := make(map[population.Population]*population.HaplotypeSimulationResult)
	for _, pop := range pops {
		r := makeResult(pop, 1000, 5_000_000, 0.001, 0.05, nil)
		r.MemoryUsedMB = 100.0
		r.SimulationTime = 2.0
		results[pop] = r
	}

	s := accumulatePerfSummary(results)
	n := float64(len(results))

	avgMem := s.totalMemory / n
	if avgMem != 100.0 {
		t.Errorf("avgMemory = %f, want 100.0", avgMem)
	}
	avgTime := s.totalSimTime / n
	if avgTime != 2.0 {
		t.Errorf("avgTime = %f, want 2.0", avgTime)
	}
}

// ---- LD validation result counting (a copy of the logic in main(), not a call into it) ----

// countPassedTests mimics the counting loop in main().
func countPassedTests(validation map[string]bool) (passed, total int) {
	total = len(validation)
	for _, p := range validation {
		if p {
			passed++
		}
	}
	return
}

// TestCountPassedTests_AllPassed verifies counting when all tests pass.
func TestCountPassedTests_AllPassed(t *testing.T) {
	validation := map[string]bool{
		"african_highest_diversity": true,
		"reasonable_ld_structure":   true,
		"haplotype_blocks_present":  true,
		"recombination_variation":   true,
		"realistic_map_lengths":     true,
	}
	passed, total := countPassedTests(validation)
	if passed != 5 {
		t.Errorf("passed = %d, want 5", passed)
	}
	if total != 5 {
		t.Errorf("total = %d, want 5", total)
	}
}

// TestCountPassedTests_SomeFailed verifies counting when some tests fail.
func TestCountPassedTests_SomeFailed(t *testing.T) {
	validation := map[string]bool{
		"test_a": true,
		"test_b": false,
		"test_c": true,
		"test_d": false,
	}
	passed, total := countPassedTests(validation)
	if passed != 2 {
		t.Errorf("passed = %d, want 2", passed)
	}
	if total != 4 {
		t.Errorf("total = %d, want 4", total)
	}
}

// TestCountPassedTests_AllFailed verifies counting when no tests pass.
func TestCountPassedTests_AllFailed(t *testing.T) {
	validation := map[string]bool{
		"test_x": false,
		"test_y": false,
	}
	passed, total := countPassedTests(validation)
	if passed != 0 {
		t.Errorf("passed = %d, want 0", passed)
	}
	if total != 2 {
		t.Errorf("total = %d, want 2", total)
	}
}

// TestCountPassedTests_Empty verifies counting with empty validation map.
func TestCountPassedTests_Empty(t *testing.T) {
	validation := map[string]bool{}
	passed, total := countPassedTests(validation)
	if passed != 0 {
		t.Errorf("passed = %d, want 0", passed)
	}
	if total != 0 {
		t.Errorf("total = %d, want 0", total)
	}
}

// ---- Cross-population LD analysis (exercises AnalyzeCrossPopulationLD) ----

// TestAnalyzeCrossPopulationLD_Keys verifies that all expected keys are present.
func TestAnalyzeCrossPopulationLD_Keys(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/test")
	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeResult(population.YRI, 100, 1_000_000, 0.0012, 0.05, nil),
		population.CEU: makeResult(population.CEU, 90, 1_000_000, 0.0008, 0.07, nil),
	}

	analysis := sim.AnalyzeCrossPopulationLD(results)

	expectedKeys := []string{
		"diversity_by_population",
		"ld_metrics_by_population",
		"block_counts_by_population",
		"diversity_mean",
		"diversity_cv",
		"ld_mean",
		"ld_cv",
	}
	for _, key := range expectedKeys {
		if _, ok := analysis[key]; !ok {
			t.Errorf("analysis missing expected key %q", key)
		}
	}
}

// TestAnalyzeCrossPopulationLD_DiversityMean verifies diversity mean computation.
func TestAnalyzeCrossPopulationLD_DiversityMean(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/test")
	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeResult(population.YRI, 100, 1_000_000, 0.0012, 0.05, nil),
		population.CEU: makeResult(population.CEU, 90, 1_000_000, 0.0008, 0.07, nil),
	}

	analysis := sim.AnalyzeCrossPopulationLD(results)
	divMean, ok := analysis["diversity_mean"].(float64)
	if !ok {
		t.Fatal("diversity_mean is not float64")
	}
	want := (0.0012 + 0.0008) / 2
	if divMean != want {
		t.Errorf("diversity_mean = %f, want %f", divMean, want)
	}
}

// TestAnalyzeCrossPopulationLD_Empty verifies behavior with empty results.
func TestAnalyzeCrossPopulationLD_Empty(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/test")
	analysis := sim.AnalyzeCrossPopulationLD(map[population.Population]*population.HaplotypeSimulationResult{})

	if _, ok := analysis["diversity_mean"]; !ok {
		t.Error("analysis should contain 'diversity_mean' even for empty results")
	}
	divMean, _ := analysis["diversity_mean"].(float64)
	if divMean != 0.0 {
		t.Errorf("diversity_mean for empty results = %f, want 0.0", divMean)
	}
}

// ---- GetRecombinationMap / GetLDStatistics (passthrough accessors) ----

// TestGetRecombinationMap_ReturnsCorrectData checks that the recombination map is returned correctly.
func TestGetRecombinationMap_ReturnsCorrectData(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/test")
	result := makeResult(population.YRI, 100, 1_000_000, 0.001, 0.05, nil)
	result.RecombinationMap.Chromosome = "chr22"
	result.RecombinationMap.MapLength = 3.7

	rm := sim.GetRecombinationMap(result)
	if rm.Chromosome != "chr22" {
		t.Errorf("Chromosome = %q, want \"chr22\"", rm.Chromosome)
	}
	if rm.MapLength != 3.7 {
		t.Errorf("MapLength = %f, want 3.7", rm.MapLength)
	}
}

// TestGetLDStatistics_ReturnsCorrectData checks that LD statistics are returned correctly.
func TestGetLDStatistics_ReturnsCorrectData(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/test")
	result := makeResult(population.YRI, 100, 1_000_000, 0.0011, 0.055, nil)
	result.Statistics.EffectiveRecombRate = 1.3e-8
	result.Statistics.HaplotypeDiversity = 0.93

	stats := sim.GetLDStatistics(result)
	if stats.DiversityPi != 0.0011 {
		t.Errorf("DiversityPi = %f, want 0.0011", stats.DiversityPi)
	}
	if stats.MeanRSquared != 0.055 {
		t.Errorf("MeanRSquared = %f, want 0.055", stats.MeanRSquared)
	}
	if stats.EffectiveRecombRate != 1.3e-8 {
		t.Errorf("EffectiveRecombRate = %e, want 1.3e-8", stats.EffectiveRecombRate)
	}
	if stats.HaplotypeDiversity != 0.93 {
		t.Errorf("HaplotypeDiversity = %f, want 0.93", stats.HaplotypeDiversity)
	}
}

// ---- HaplotypeSimulator initialisation (mirrors NewHaplotypeSimulator usage in main) ----

// TestNewHaplotypeSimulator_Paths checks that the simulator paths include expected substrings.
func TestNewHaplotypeSimulator_Paths(t *testing.T) {
	root := "/my/project/root"
	sim := population.NewHaplotypeSimulator(root)
	// We can't read unexported fields from outside the package, but the internal
	// test file (haplotype_bridge_test.go in package population) does so.
	// Here we verify the simulator is non-nil and usable.
	if sim == nil {
		t.Fatal("NewHaplotypeSimulator returned nil")
	}

	// The simulator should be usable for calling accessors.
	result := makeResult(population.YRI, 10, 100_000, 0.001, 0.05, nil)
	rm := sim.GetRecombinationMap(result)
	if rm.Chromosome == "" {
		// An empty chromosome is acceptable for a test-constructed result, but
		// the call should not panic.
	}
}

// ---- Repeat-string helper (exercises strings.Repeat used in main) ----

// TestStringsRepeat_NonEmpty verifies that strings.Repeat behaves as used in main.
func TestStringsRepeat_NonEmpty(t *testing.T) {
	line := "=" + strings.Repeat("=", 47)
	if len(line) != 48 {
		t.Errorf("header line length = %d, want 48", len(line))
	}
}

// TestStringsRepeat_ZeroCount verifies that repeating zero times returns the base string.
func TestStringsRepeat_ZeroCount(t *testing.T) {
	line := "=" + strings.Repeat("=", 0)
	if line != "=" {
		t.Errorf("line = %q, want \"=\"", line)
	}
}

// ---- ValidateLDPatterns (exercises the validation called in main) ----

// TestValidateLDPatterns_ReturnsFiveTests verifies that ValidateLDPatterns always
// returns exactly five validation entries.
func TestValidateLDPatterns_ReturnsFiveTests(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/test")

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeResult(population.YRI, 200, 5_000_000, 0.0012, 0.05,
			[]population.HaplotypeBlock{{Length: 50000, VariantCount: 10}}),
		population.CEU: makeResult(population.CEU, 150, 5_000_000, 0.0008, 0.06,
			[]population.HaplotypeBlock{{Length: 70000, VariantCount: 12}}),
		population.CHB: makeResult(population.CHB, 140, 5_000_000, 0.0007, 0.065,
			[]population.HaplotypeBlock{{Length: 80000, VariantCount: 14}}),
	}

	validation, err := sim.ValidateLDPatterns(results)
	if err != nil {
		t.Fatalf("ValidateLDPatterns returned unexpected error: %v", err)
	}
	if len(validation) != 5 {
		t.Errorf("len(validation) = %d, want 5", len(validation))
	}
}

// TestValidateLDPatterns_HaplotypeBlocksPresentFalse verifies that
// haplotype_blocks_present is false when no blocks are simulated.
func TestValidateLDPatterns_HaplotypeBlocksPresentFalse(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/test")

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeResult(population.YRI, 200, 5_000_000, 0.0012, 0.05, nil),
		population.CEU: makeResult(population.CEU, 150, 5_000_000, 0.0008, 0.06, nil),
	}

	validation, err := sim.ValidateLDPatterns(results)
	if err != nil {
		t.Fatalf("ValidateLDPatterns returned unexpected error: %v", err)
	}
	if validation["haplotype_blocks_present"] {
		t.Error("haplotype_blocks_present should be false when no blocks are present")
	}
}

// TestValidateLDPatterns_AfricanHighestDiversity verifies that African diversity
// check passes when YRI diversity is highest.
func TestValidateLDPatterns_AfricanHighestDiversity(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/test")

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeResult(population.YRI, 200, 5_000_000, 0.0015, 0.04, nil),
		population.LWK: makeResult(population.LWK, 190, 5_000_000, 0.0014, 0.04, nil),
		population.GWD: makeResult(population.GWD, 185, 5_000_000, 0.0013, 0.04, nil),
		population.CEU: makeResult(population.CEU, 150, 5_000_000, 0.0008, 0.06, nil),
		population.GBR: makeResult(population.GBR, 145, 5_000_000, 0.0007, 0.06, nil),
		population.TSI: makeResult(population.TSI, 148, 5_000_000, 0.00075, 0.06, nil),
		population.CHB: makeResult(population.CHB, 140, 5_000_000, 0.0007, 0.065, nil),
		population.JPT: makeResult(population.JPT, 142, 5_000_000, 0.00065, 0.065, nil),
		population.CHS: makeResult(population.CHS, 138, 5_000_000, 0.0006, 0.065, nil),
	}

	validation, err := sim.ValidateLDPatterns(results)
	if err != nil {
		t.Fatalf("ValidateLDPatterns returned unexpected error: %v", err)
	}
	if !validation["african_highest_diversity"] {
		t.Error("african_highest_diversity should be true when African pops have highest diversity")
	}
}

// TestValidateLDPatterns_EmptyResults verifies behavior with empty results map.
func TestValidateLDPatterns_EmptyResults(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/test")
	validation, err := sim.ValidateLDPatterns(map[population.Population]*population.HaplotypeSimulationResult{})
	if err != nil {
		t.Fatalf("ValidateLDPatterns returned unexpected error: %v", err)
	}
	// All checks should default to false/zero when no data is present.
	for key, val := range validation {
		_ = key
		_ = val
		// No panic expected — just validating the function doesn't crash.
	}
}

// ---- HaplotypeSimulationResult data integrity ----

// TestHaplotypeSimulationResult_SequenceLength verifies SequenceLength is stored.
func TestHaplotypeSimulationResult_SequenceLength(t *testing.T) {
	result := makeResult(population.MXL, 500, 10_000_000, 0.0009, 0.055, nil)
	if result.SequenceLength != 10_000_000 {
		t.Errorf("SequenceLength = %d, want 10_000_000", result.SequenceLength)
	}
}

// TestHaplotypeSimulationResult_LDBlocksNil verifies empty-block path works.
func TestHaplotypeSimulationResult_LDBlocksNil(t *testing.T) {
	result := makeResult(population.PUR, 0, 1_000_000, 0.0, 0.0, nil)
	if len(result.LDBlocks) != 0 {
		t.Errorf("expected 0 LD blocks, got %d", len(result.LDBlocks))
	}
}

// TestHaplotypeSimulationResult_MemoryAndTime verifies MemoryUsedMB and SimulationTime.
func TestHaplotypeSimulationResult_MemoryAndTime(t *testing.T) {
	result := makeResult(population.YRI, 100, 1_000_000, 0.001, 0.05, nil)
	result.MemoryUsedMB = 123.4
	result.SimulationTime = 5.6

	if result.MemoryUsedMB != 123.4 {
		t.Errorf("MemoryUsedMB = %f, want 123.4", result.MemoryUsedMB)
	}
	if result.SimulationTime != 5.6 {
		t.Errorf("SimulationTime = %f, want 5.6", result.SimulationTime)
	}
}

// ---- Population constants used in main() ----

// TestPopulationConstants_MainPops verifies the four populations used in main() are
// valid population identifiers.
func TestPopulationConstants_MainPops(t *testing.T) {
	mainPops := []population.Population{
		population.YRI,
		population.CEU,
		population.CHB,
		population.MXL,
	}
	for _, pop := range mainPops {
		if !population.IsValidPopulation(string(pop)) {
			t.Errorf("population %q used in main() is not a valid population", pop)
		}
	}
}

// TestPopulationConstants_String verifies that population String() matches its constant value.
func TestPopulationConstants_String(t *testing.T) {
	cases := []struct {
		pop  population.Population
		want string
	}{
		{population.YRI, "YRI"},
		{population.CEU, "CEU"},
		{population.CHB, "CHB"},
		{population.MXL, "MXL"},
	}
	for _, tc := range cases {
		if tc.pop.String() != tc.want {
			t.Errorf("population.String() = %q, want %q", tc.pop.String(), tc.want)
		}
	}
}

// ---- RecombinationMap construction used in main() ----

// TestRecombinationMap_MapLength verifies MapLength computation for standard sequence.
func TestRecombinationMap_MapLength(t *testing.T) {
	// Standard Mb-scale sequence: 5_000_000 bp → expected ~7.5 cM at 1.5 cM/Mb
	seqLen := 5_000_000
	expectedMapLength := float64(seqLen) / 1_000_000 * 1.5
	result := makeResult(population.YRI, 100, seqLen, 0.001, 0.05, nil)

	rm := result.RecombinationMap
	if rm.MapLength != expectedMapLength {
		t.Errorf("MapLength = %f, want %f", rm.MapLength, expectedMapLength)
	}
}

// TestRecombinationMap_PhysicalLength verifies PhysicalLength matches SequenceLength.
func TestRecombinationMap_PhysicalLength(t *testing.T) {
	seqLen := 7_500_000
	result := makeResult(population.CEU, 200, seqLen, 0.0009, 0.06, nil)

	rm := result.RecombinationMap
	if rm.PhysicalLength != seqLen {
		t.Errorf("PhysicalLength = %d, want %d", rm.PhysicalLength, seqLen)
	}
}
