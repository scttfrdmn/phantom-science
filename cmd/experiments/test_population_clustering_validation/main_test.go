//go:build !release
// +build !release

package main

import (
	"math"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/population"
)

// ---------------------------------------------------------------------------
// Copies of logic in main(), tested in isolation. main() does not call
// these, so the tests do not exercise main() and the copies can drift.
// ---------------------------------------------------------------------------

// computeDiversityStats mirrors the diversity statistics block in main().
func computeDiversityStats(diversities map[population.Population]float64) (mean, rangeVal, cv float64) {
	if len(diversities) == 0 {
		return 0, 0, 0
	}

	var sum float64
	var values []float64
	for _, d := range diversities {
		sum += d
		values = append(values, d)
	}

	mean = sum / float64(len(diversities))

	minD, maxD := values[0], values[0]
	for _, d := range values {
		if d < minD {
			minD = d
		}
		if d > maxD {
			maxD = d
		}
	}

	rangeVal = maxD - minD
	if mean > 0 {
		cv = rangeVal / mean
	}
	return
}

// clusteringReadiness mirrors the readiness scoring in main().
func clusteringReadiness(numVariants int, numLDBlocks int, diversityPi float64) int {
	score := 0
	if numVariants > 200 {
		score++
	}
	if numLDBlocks > 5 {
		score++
	}
	if diversityPi > 0.00001 {
		score++
	}
	return score
}

// hierarchyPassed mirrors the hierarchy-test logic from main().
func hierarchyPassed(
	yri, ceu, chb float64,
	hasEuropean, hasEastAsian bool,
) (passed, total int) {
	if hasEuropean {
		total++
		if yri > ceu {
			passed++
		}
	}
	if hasEastAsian {
		total++
		if yri > chb {
			passed++
		}
	}
	if hasEuropean && hasEastAsian {
		total++
		if ceu > chb {
			passed++
		}
	}
	return
}

// qualityMetrics mirrors the final assessment block from main().
func buildQualityMetrics(
	numHaplotypeResults int,
	clusteringReadyCount int,
	continentsRepresented int,
) map[string]bool {
	return map[string]bool{
		"PCA clustering validation framework available": true,
		"Sequence generation operational":               numHaplotypeResults >= 3,
		"Population differentiation detected":           numHaplotypeResults >= 3,
		"LD structure analysis functional":              clusteringReadyCount >= 2,
		"Continental groups represented":                continentsRepresented >= 3,
		"Clustering readiness achieved":                 clusteringReadyCount >= numHaplotypeResults*2/3,
		"Population diversity hierarchy established":    true,
		"Phase 2 Week 6 foundation complete":            true,
	}
}

// ---------------------------------------------------------------------------
// Tests for diversity statistics
// ---------------------------------------------------------------------------

func TestComputeDiversityStats_BasicValues(t *testing.T) {
	diversities := map[population.Population]float64{
		population.YRI: 0.0009,
		population.CEU: 0.0006,
		population.CHB: 0.0004,
	}

	mean, rangeVal, cv := computeDiversityStats(diversities)

	expectedMean := (0.0009 + 0.0006 + 0.0004) / 3.0
	if math.Abs(mean-expectedMean) > 1e-12 {
		t.Errorf("mean = %v, want %v", mean, expectedMean)
	}

	expectedRange := 0.0009 - 0.0004
	if math.Abs(rangeVal-expectedRange) > 1e-12 {
		t.Errorf("range = %v, want %v", rangeVal, expectedRange)
	}

	expectedCV := expectedRange / expectedMean
	if math.Abs(cv-expectedCV) > 1e-10 {
		t.Errorf("cv = %v, want %v", cv, expectedCV)
	}
}

func TestComputeDiversityStats_SingleValue(t *testing.T) {
	diversities := map[population.Population]float64{
		population.YRI: 0.001,
	}

	mean, rangeVal, cv := computeDiversityStats(diversities)

	if math.Abs(mean-0.001) > 1e-12 {
		t.Errorf("mean = %v, want 0.001", mean)
	}
	if rangeVal != 0 {
		t.Errorf("range = %v, want 0", rangeVal)
	}
	if cv != 0 {
		t.Errorf("cv = %v, want 0", cv)
	}
}

func TestComputeDiversityStats_EmptyMap(t *testing.T) {
	mean, rangeVal, cv := computeDiversityStats(nil)
	if mean != 0 || rangeVal != 0 || cv != 0 {
		t.Errorf("expected all zeros for empty input, got mean=%v range=%v cv=%v", mean, rangeVal, cv)
	}
}

func TestComputeDiversityStats_ZeroMean(t *testing.T) {
	diversities := map[population.Population]float64{
		population.YRI: 0,
		population.CEU: 0,
	}

	_, _, cv := computeDiversityStats(diversities)
	if cv != 0 {
		t.Errorf("cv with zero mean should be 0, got %v", cv)
	}
}

func TestComputeDiversityStats_AllSameValues(t *testing.T) {
	diversities := map[population.Population]float64{
		population.YRI: 0.005,
		population.CEU: 0.005,
		population.CHB: 0.005,
		population.MXL: 0.005,
	}

	mean, rangeVal, cv := computeDiversityStats(diversities)

	if math.Abs(mean-0.005) > 1e-12 {
		t.Errorf("mean = %v, want 0.005", mean)
	}
	if rangeVal != 0 {
		t.Errorf("range = %v, want 0", rangeVal)
	}
	if cv != 0 {
		t.Errorf("cv = %v, want 0", cv)
	}
}

// ---------------------------------------------------------------------------
// Tests for clustering readiness scoring
// ---------------------------------------------------------------------------

func TestClusteringReadiness_AllCriteriaMet(t *testing.T) {
	score := clusteringReadiness(500, 10, 0.001)
	if score != 3 {
		t.Errorf("score = %d, want 3", score)
	}
}

func TestClusteringReadiness_NoVariants(t *testing.T) {
	score := clusteringReadiness(50, 10, 0.001)
	if score != 2 {
		t.Errorf("score = %d, want 2 (no variants threshold)", score)
	}
}

func TestClusteringReadiness_NoLDBlocks(t *testing.T) {
	score := clusteringReadiness(500, 2, 0.001)
	if score != 2 {
		t.Errorf("score = %d, want 2 (no LD blocks threshold)", score)
	}
}

func TestClusteringReadiness_NoDiversity(t *testing.T) {
	score := clusteringReadiness(500, 10, 0.000001)
	if score != 2 {
		t.Errorf("score = %d, want 2 (no diversity threshold)", score)
	}
}

func TestClusteringReadiness_NoCriteriaMet(t *testing.T) {
	score := clusteringReadiness(10, 2, 0.000001)
	if score != 0 {
		t.Errorf("score = %d, want 0", score)
	}
}

func TestClusteringReadiness_ExactBoundaryVariants(t *testing.T) {
	// 200 is NOT > 200, so it should not count
	score200 := clusteringReadiness(200, 10, 0.001)
	// 201 IS > 200
	score201 := clusteringReadiness(201, 10, 0.001)

	if score200 == score201 {
		t.Errorf("score at 200 (%d) should be less than score at 201 (%d)", score200, score201)
	}
}

func TestClusteringReadiness_ExactBoundaryLDBlocks(t *testing.T) {
	// 5 is NOT > 5
	score5 := clusteringReadiness(500, 5, 0.001)
	// 6 IS > 5
	score6 := clusteringReadiness(500, 6, 0.001)

	if score5 == score6 {
		t.Errorf("score at 5 blocks (%d) should be less than score at 6 blocks (%d)", score5, score6)
	}
}

// ---------------------------------------------------------------------------
// Tests for diversity hierarchy checking
// ---------------------------------------------------------------------------

func TestHierarchyPassed_AllExpected(t *testing.T) {
	// Expected hierarchy: YRI > CEU > CHB
	passed, total := hierarchyPassed(0.0009, 0.0006, 0.0004, true, true)
	if total != 3 {
		t.Errorf("total = %d, want 3", total)
	}
	if passed != 3 {
		t.Errorf("passed = %d, want 3", passed)
	}
}

func TestHierarchyPassed_NoneExpected(t *testing.T) {
	// Inverted hierarchy
	passed, total := hierarchyPassed(0.0003, 0.0006, 0.0009, true, true)
	if total != 3 {
		t.Errorf("total = %d, want 3", total)
	}
	if passed != 0 {
		t.Errorf("passed = %d, want 0", passed)
	}
}

func TestHierarchyPassed_OnlyEuropean(t *testing.T) {
	passed, total := hierarchyPassed(0.0009, 0.0006, 0.0004, true, false)
	if total != 1 {
		t.Errorf("total = %d, want 1", total)
	}
	if passed != 1 {
		t.Errorf("passed = %d, want 1", passed)
	}
}

func TestHierarchyPassed_OnlyEastAsian(t *testing.T) {
	passed, total := hierarchyPassed(0.0009, 0.0006, 0.0004, false, true)
	if total != 1 {
		t.Errorf("total = %d, want 1", total)
	}
	if passed != 1 {
		t.Errorf("passed = %d, want 1", passed)
	}
}

func TestHierarchyPassed_NeitherPresent(t *testing.T) {
	passed, total := hierarchyPassed(0.0009, 0.0006, 0.0004, false, false)
	if total != 0 {
		t.Errorf("total = %d, want 0", total)
	}
	if passed != 0 {
		t.Errorf("passed = %d, want 0", passed)
	}
}

func TestHierarchyPassed_YRILessThanCEU(t *testing.T) {
	// YRI=0.0003 < CEU=0.0009, CHB=0.0004
	// YRI > CEU? No → fail
	// YRI > CHB? No (0.0003 < 0.0004) → fail
	// CEU > CHB? Yes (0.0009 > 0.0004) → pass
	// => 1 passed out of 3
	passed, total := hierarchyPassed(0.0003, 0.0009, 0.0004, true, true)
	if total != 3 {
		t.Errorf("total = %d, want 3", total)
	}
	if passed != 1 {
		t.Errorf("passed = %d, want 1", passed)
	}
}

// ---------------------------------------------------------------------------
// Tests for quality metrics computation
// ---------------------------------------------------------------------------

func TestBuildQualityMetrics_AllPassing(t *testing.T) {
	metrics := buildQualityMetrics(4, 4, 4)

	totalQuality := len(metrics)
	passedQuality := 0
	for _, passed := range metrics {
		if passed {
			passedQuality++
		}
	}

	if totalQuality != 8 {
		t.Errorf("expected 8 quality metrics, got %d", totalQuality)
	}
	if passedQuality != 8 {
		t.Errorf("expected all 8 quality metrics to pass, got %d", passedQuality)
	}
}

func TestBuildQualityMetrics_PartialPassing(t *testing.T) {
	// Only 1 haplotype result, 0 clustering-ready, 1 continent
	metrics := buildQualityMetrics(1, 0, 1)

	// These should be false:
	// "Sequence generation operational": 1 >= 3 → false
	// "Population differentiation detected": 1 >= 3 → false
	// "LD structure analysis functional": 0 >= 2 → false
	// "Continental groups represented": 1 >= 3 → false
	// "Clustering readiness achieved": 0 >= 1*2/3 = 0 → true (0 >= 0 is true in integer division)

	if metrics["PCA clustering validation framework available"] != true {
		t.Error("PCA framework should always be true")
	}
	if metrics["Sequence generation operational"] != false {
		t.Error("sequence generation should be false for 1 result")
	}
	if metrics["Population differentiation detected"] != false {
		t.Error("population differentiation should be false for 1 result")
	}
	if metrics["LD structure analysis functional"] != false {
		t.Error("LD structure should be false for 0 clustering-ready")
	}
	if metrics["Continental groups represented"] != false {
		t.Error("continental groups should be false for 1 continent")
	}
}

func TestBuildQualityMetrics_ThresholdBoundary(t *testing.T) {
	// Exactly 3 results, 2 clustering-ready, 3 continents
	metrics := buildQualityMetrics(3, 2, 3)

	if !metrics["Sequence generation operational"] {
		t.Error("sequence generation operational: 3 >= 3 should be true")
	}
	if !metrics["LD structure analysis functional"] {
		t.Error("LD structure functional: 2 >= 2 should be true")
	}
	if !metrics["Continental groups represented"] {
		t.Error("continental groups: 3 >= 3 should be true")
	}
}

// ---------------------------------------------------------------------------
// Tests for population constants (referenced in main())
// ---------------------------------------------------------------------------

func TestPopulationConstants(t *testing.T) {
	tests := []struct {
		pop      population.Population
		expected string
	}{
		{population.YRI, "YRI"},
		{population.CEU, "CEU"},
		{population.CHB, "CHB"},
		{population.MXL, "MXL"},
	}

	for _, tt := range tests {
		if string(tt.pop) != tt.expected {
			t.Errorf("population %v string = %q, want %q", tt.pop, string(tt.pop), tt.expected)
		}
	}
}

func TestGetSupportedPopulations(t *testing.T) {
	pops := population.GetSupportedPopulations()
	if len(pops) != 11 {
		t.Errorf("expected 11 supported populations, got %d", len(pops))
	}

	required := []population.Population{
		population.YRI, population.CEU, population.CHB, population.MXL,
	}
	for _, req := range required {
		found := false
		for _, p := range pops {
			if p == req {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("required population %s not in supported list", req)
		}
	}
}

func TestIsValidPopulation(t *testing.T) {
	validPops := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"}
	for _, p := range validPops {
		if !population.IsValidPopulation(p) {
			t.Errorf("expected %s to be valid", p)
		}
	}

	invalidPops := []string{"INVALID", "", "abc", "EUR", "AFR"}
	for _, p := range invalidPops {
		if population.IsValidPopulation(p) {
			t.Errorf("expected %s to be invalid", p)
		}
	}
}

func TestPopulationString(t *testing.T) {
	if population.YRI.String() != "YRI" {
		t.Errorf("YRI.String() = %q, want %q", population.YRI.String(), "YRI")
	}
	if population.CEU.String() != "CEU" {
		t.Errorf("CEU.String() = %q, want %q", population.CEU.String(), "CEU")
	}
}

// ---------------------------------------------------------------------------
// Tests for HaplotypeSimulator construction (NewHaplotypeSimulator)
// ---------------------------------------------------------------------------

func TestNewHaplotypeSimulator_NotNil(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")
	if sim == nil {
		t.Error("expected non-nil HaplotypeSimulator")
	}
}

func TestNewHaplotypeSimulator_WithDifferentRoots(t *testing.T) {
	roots := []string{"/tmp", "/home/user", "/opt/project", ""}
	for _, root := range roots {
		sim := population.NewHaplotypeSimulator(root)
		if sim == nil {
			t.Errorf("expected non-nil simulator for root %q", root)
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for clustering continental grouping logic
// ---------------------------------------------------------------------------

func TestContinentalGrouping_AllPresent(t *testing.T) {
	haplotypeResults := map[population.Population]bool{
		population.YRI: true,
		population.CEU: true,
		population.CHB: true,
		population.MXL: true,
	}

	continentalGroups := map[string][]population.Population{
		"African":    {population.YRI},
		"European":   {population.CEU},
		"East_Asian": {population.CHB},
		"Americas":   {population.MXL},
	}

	continentsRepresented := 0
	for _, pops := range continentalGroups {
		for _, pop := range pops {
			if haplotypeResults[pop] {
				continentsRepresented++
				break
			}
		}
	}

	if continentsRepresented != 4 {
		t.Errorf("expected 4 continents represented, got %d", continentsRepresented)
	}
}

func TestContinentalGrouping_SomeAbsent(t *testing.T) {
	haplotypeResults := map[population.Population]bool{
		population.YRI: true,
		population.CEU: true,
		// CHB and MXL absent
	}

	continentalGroups := map[string][]population.Population{
		"African":    {population.YRI},
		"European":   {population.CEU},
		"East_Asian": {population.CHB},
		"Americas":   {population.MXL},
	}

	continentsRepresented := 0
	for _, pops := range continentalGroups {
		for _, pop := range pops {
			if haplotypeResults[pop] {
				continentsRepresented++
				break
			}
		}
	}

	if continentsRepresented != 2 {
		t.Errorf("expected 2 continents represented, got %d", continentsRepresented)
	}
}

func TestContinentalGrouping_ThresholdCheck(t *testing.T) {
	// At least 3 continents needed for "sufficient diversity"
	testCases := []struct {
		represented int
		wantSuff    bool
	}{
		{4, true},
		{3, true},
		{2, false},
		{1, false},
		{0, false},
	}

	for _, tc := range testCases {
		gotSuff := tc.represented >= 3
		if gotSuff != tc.wantSuff {
			t.Errorf("represented=%d: sufficient=%v, want %v", tc.represented, gotSuff, tc.wantSuff)
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for phase assessment thresholds
// ---------------------------------------------------------------------------

func TestPhaseAssessmentThreshold(t *testing.T) {
	// main() checks: passedQuality >= totalQuality*3/4
	testCases := []struct {
		passed   int
		total    int
		wantPass bool
	}{
		{8, 8, true},
		{6, 8, true},
		{5, 8, false},
		{3, 4, true},
		{2, 4, false},
		{0, 8, false},
	}

	for _, tc := range testCases {
		got := tc.passed >= tc.total*3/4
		if got != tc.wantPass {
			t.Errorf("passed=%d total=%d: got %v, want %v", tc.passed, tc.total, got, tc.wantPass)
		}
	}
}

func TestClusteringReadyThreshold(t *testing.T) {
	// main() checks: clusteringReadyCount >= len(haplotypeResults)*3/4
	testCases := []struct {
		ready    int
		total    int
		wantPass bool
	}{
		{4, 4, true},
		{3, 4, true},
		{2, 4, false},
		{1, 4, false},
	}

	for _, tc := range testCases {
		got := tc.ready >= tc.total*3/4
		if got != tc.wantPass {
			t.Errorf("ready=%d total=%d: got %v, want %v", tc.ready, tc.total, got, tc.wantPass)
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for HaplotypeSimulator exported methods using constructed results
// ---------------------------------------------------------------------------

func makeTestHaplotypeResult(pop population.Population, numVariants int, ldBlocks []population.HaplotypeBlock, stats population.HaplotypeStatistics) *population.HaplotypeSimulationResult {
	return &population.HaplotypeSimulationResult{
		Population:  pop,
		Chromosome:  "chr22",
		NumVariants: numVariants,
		NumTrees:    10,
		NumSamples:  25,
		LDBlocks:    ldBlocks,
		Statistics:  stats,
		RecombinationMap: population.RecombinationMap{
			Chromosome:     "chr22",
			Positions:      []int{1000, 2000, 3000},
			Rates:          []float64{1e-8, 1.2e-8, 0.9e-8},
			MapLength:      3.0,
			PhysicalLength: 3000000,
		},
	}
}

func TestGetRecombinationMap(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")
	blocks := []population.HaplotypeBlock{
		{StartPosition: 0, EndPosition: 100000, Length: 100000, VariantCount: 50, RSquaredMean: 0.4},
	}
	stats := population.HaplotypeStatistics{
		DiversityPi:        0.0009,
		MeanRSquared:       0.04,
		GeneticMapLengthCM: 3.0,
	}
	result := makeTestHaplotypeResult(population.YRI, 500, blocks, stats)

	recombMap := sim.GetRecombinationMap(result)
	if recombMap.Chromosome != "chr22" {
		t.Errorf("chromosome = %q, want %q", recombMap.Chromosome, "chr22")
	}
	if recombMap.MapLength != 3.0 {
		t.Errorf("map length = %v, want 3.0", recombMap.MapLength)
	}
	if len(recombMap.Positions) != 3 {
		t.Errorf("positions count = %d, want 3", len(recombMap.Positions))
	}
}

func TestGetLDStatistics(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")
	blocks := []population.HaplotypeBlock{}
	stats := population.HaplotypeStatistics{
		DiversityPi:      0.0007,
		MeanRSquared:     0.03,
		SegregatingSites: 350,
	}
	result := makeTestHaplotypeResult(population.CEU, 350, blocks, stats)

	ldStats := sim.GetLDStatistics(result)
	if math.Abs(ldStats.DiversityPi-0.0007) > 1e-10 {
		t.Errorf("DiversityPi = %v, want 0.0007", ldStats.DiversityPi)
	}
	if math.Abs(ldStats.MeanRSquared-0.03) > 1e-10 {
		t.Errorf("MeanRSquared = %v, want 0.03", ldStats.MeanRSquared)
	}
}

func TestAnalyzeCrossPopulationLD_BasicStructure(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeTestHaplotypeResult(population.YRI, 600, []population.HaplotypeBlock{
			{StartPosition: 0, EndPosition: 100000},
		}, population.HaplotypeStatistics{DiversityPi: 0.0009, MeanRSquared: 0.03}),
		population.CEU: makeTestHaplotypeResult(population.CEU, 450, []population.HaplotypeBlock{
			{StartPosition: 0, EndPosition: 200000},
			{StartPosition: 200000, EndPosition: 400000},
		}, population.HaplotypeStatistics{DiversityPi: 0.0006, MeanRSquared: 0.05}),
		population.CHB: makeTestHaplotypeResult(population.CHB, 400, []population.HaplotypeBlock{
			{StartPosition: 0, EndPosition: 150000},
		}, population.HaplotypeStatistics{DiversityPi: 0.0005, MeanRSquared: 0.06}),
	}

	analysis := sim.AnalyzeCrossPopulationLD(results)

	if analysis == nil {
		t.Fatal("expected non-nil analysis")
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
			t.Errorf("analysis missing key %q", key)
		}
	}
}

func TestAnalyzeCrossPopulationLD_EmptyResults(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")
	results := map[population.Population]*population.HaplotypeSimulationResult{}
	analysis := sim.AnalyzeCrossPopulationLD(results)

	if analysis == nil {
		t.Fatal("expected non-nil analysis even for empty results")
	}

	if _, ok := analysis["diversity_mean"]; !ok {
		t.Error("analysis should have diversity_mean key even for empty results")
	}
}

func TestAnalyzeCrossPopulationLD_DiversityValues(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeTestHaplotypeResult(population.YRI, 600, nil,
			population.HaplotypeStatistics{DiversityPi: 0.0009, MeanRSquared: 0.03}),
		population.CEU: makeTestHaplotypeResult(population.CEU, 450, nil,
			population.HaplotypeStatistics{DiversityPi: 0.0006, MeanRSquared: 0.05}),
	}

	analysis := sim.AnalyzeCrossPopulationLD(results)

	diversityMean, ok := analysis["diversity_mean"].(float64)
	if !ok {
		t.Fatal("diversity_mean should be float64")
	}
	// Mean of 0.0009 and 0.0006 = 0.00075
	expectedMean := (0.0009 + 0.0006) / 2.0
	if math.Abs(diversityMean-expectedMean) > 1e-10 {
		t.Errorf("diversity_mean = %v, want %v", diversityMean, expectedMean)
	}
}

func TestValidateLDPatterns_WithSufficientData(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")

	// Create results with African > European/EastAsian diversity
	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeTestHaplotypeResult(population.YRI, 600,
			[]population.HaplotypeBlock{
				{StartPosition: 0, EndPosition: 500000},
				{StartPosition: 500000, EndPosition: 1000000},
			},
			population.HaplotypeStatistics{
				DiversityPi:         0.0009,
				MeanRSquared:        0.03,
				EffectiveRecombRate: 1.0e-8,
				GeneticMapLengthCM:  4.5,
			}),
		population.CEU: makeTestHaplotypeResult(population.CEU, 450,
			[]population.HaplotypeBlock{
				{StartPosition: 0, EndPosition: 700000},
			},
			population.HaplotypeStatistics{
				DiversityPi:         0.0006,
				MeanRSquared:        0.05,
				EffectiveRecombRate: 1.2e-8,
				GeneticMapLengthCM:  4.0,
			}),
		population.CHB: makeTestHaplotypeResult(population.CHB, 400,
			[]population.HaplotypeBlock{
				{StartPosition: 0, EndPosition: 800000},
			},
			population.HaplotypeStatistics{
				DiversityPi:         0.0005,
				MeanRSquared:        0.06,
				EffectiveRecombRate: 0.9e-8,
				GeneticMapLengthCM:  3.5,
			}),
	}

	validation, err := sim.ValidateLDPatterns(results)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if validation == nil {
		t.Fatal("expected non-nil validation map")
	}

	// African diversity is highest (0.0009 > 0.0006 and 0.0005)
	if !validation["african_highest_diversity"] {
		t.Error("african_highest_diversity should be true when YRI > CEU and CHB")
	}
}

func TestValidateLDPatterns_HaplotypeBlocksPresent(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeTestHaplotypeResult(population.YRI, 600,
			[]population.HaplotypeBlock{
				{StartPosition: 0, EndPosition: 500000},
			},
			population.HaplotypeStatistics{
				DiversityPi:         0.0009,
				MeanRSquared:        0.04,
				EffectiveRecombRate: 1.0e-8,
				GeneticMapLengthCM:  4.5,
			}),
	}

	validation, err := sim.ValidateLDPatterns(results)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if !validation["haplotype_blocks_present"] {
		t.Error("haplotype_blocks_present should be true when blocks exist")
	}
}

func TestValidateLDPatterns_EmptyResults(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")

	results := map[population.Population]*population.HaplotypeSimulationResult{}

	validation, err := sim.ValidateLDPatterns(results)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if validation == nil {
		t.Fatal("expected non-nil validation map")
	}

	// No blocks present
	if validation["haplotype_blocks_present"] {
		t.Error("haplotype_blocks_present should be false for empty results")
	}
}

func TestValidateLDPatterns_ReturnsMultipleKeys(t *testing.T) {
	sim := population.NewHaplotypeSimulator("/tmp/test-project")

	results := map[population.Population]*population.HaplotypeSimulationResult{
		population.YRI: makeTestHaplotypeResult(population.YRI, 600, nil,
			population.HaplotypeStatistics{
				DiversityPi:         0.0009,
				MeanRSquared:        0.04,
				EffectiveRecombRate: 1.0e-8,
				GeneticMapLengthCM:  4.5,
			}),
	}

	validation, err := sim.ValidateLDPatterns(results)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	expectedKeys := []string{
		"african_highest_diversity",
		"reasonable_ld_structure",
		"haplotype_blocks_present",
		"recombination_variation",
		"realistic_map_lengths",
	}
	for _, key := range expectedKeys {
		if _, ok := validation[key]; !ok {
			t.Errorf("validation missing key %q", key)
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for HaplotypeBlock and HaplotypeVariant types (used in main())
// ---------------------------------------------------------------------------

func TestHaplotypeBlockFields(t *testing.T) {
	block := population.HaplotypeBlock{
		StartPosition: 100000,
		EndPosition:   200000,
		Length:        100000,
		VariantCount:  50,
		RSquaredMean:  0.35,
		Population:    population.YRI,
	}

	if block.StartPosition != 100000 {
		t.Errorf("StartPosition = %d, want 100000", block.StartPosition)
	}
	if block.EndPosition != 200000 {
		t.Errorf("EndPosition = %d, want 200000", block.EndPosition)
	}
	if block.Population != population.YRI {
		t.Errorf("Population = %v, want YRI", block.Population)
	}
}

func TestHaplotypeVariantFields(t *testing.T) {
	variant := population.HaplotypeVariant{
		Position:       12345,
		Alleles:        []string{"A", "G"},
		AncestralState: "A",
		Frequency:      0.15,
	}

	if variant.Position != 12345 {
		t.Errorf("Position = %d, want 12345", variant.Position)
	}
	if len(variant.Alleles) != 2 {
		t.Errorf("Alleles length = %d, want 2", len(variant.Alleles))
	}
	if variant.Frequency != 0.15 {
		t.Errorf("Frequency = %v, want 0.15", variant.Frequency)
	}
}

func TestHaplotypeStatisticsFields(t *testing.T) {
	stats := population.HaplotypeStatistics{
		DiversityPi:      0.0008,
		SegregatingSites: 400,
		NumVariants:      400,
		NumHaplotypes:    50,
		MeanRSquared:     0.04,
		MedianRSquared:   0.03,
		HighLDFraction:   0.15,
	}

	if math.Abs(stats.DiversityPi-0.0008) > 1e-12 {
		t.Errorf("DiversityPi = %v, want 0.0008", stats.DiversityPi)
	}
	if stats.SegregatingSites != 400 {
		t.Errorf("SegregatingSites = %d, want 400", stats.SegregatingSites)
	}
}

// ---------------------------------------------------------------------------
// Tests for GetPopulationParameters (bridge.go)
// ---------------------------------------------------------------------------

func TestGetPopulationParameters_KnownPops(t *testing.T) {
	sim := population.NewDemographicSimulator("/tmp/test-project")

	knownPops := []population.Population{
		population.YRI, population.CEU, population.CHB, population.MXL,
	}

	for _, pop := range knownPops {
		params, err := sim.GetPopulationParameters(pop)
		if err != nil {
			t.Errorf("unexpected error for %s: %v", pop, err)
			continue
		}
		if params == nil {
			t.Errorf("expected non-nil params for %s", pop)
			continue
		}
		if params.Population != pop {
			t.Errorf("params.Population = %v, want %v", params.Population, pop)
		}
		if params.EffectivePopulation <= 0 {
			t.Errorf("%s effective population = %d, expected > 0", pop, params.EffectivePopulation)
		}
	}
}

func TestGetPopulationParameters_UnknownPop(t *testing.T) {
	sim := population.NewDemographicSimulator("/tmp/test-project")
	_, err := sim.GetPopulationParameters(population.Population("UNKNOWN"))
	if err == nil {
		t.Error("expected error for unknown population, got nil")
	}
}

func TestGetPopulationParameters_YRIHasHigherNe(t *testing.T) {
	sim := population.NewDemographicSimulator("/tmp/test-project")

	yriParams, _ := sim.GetPopulationParameters(population.YRI)
	chbParams, _ := sim.GetPopulationParameters(population.CHB)

	if yriParams == nil || chbParams == nil {
		t.Skip("params not available")
	}

	// YRI (12300) > CHB (9025)
	if yriParams.EffectivePopulation <= chbParams.EffectivePopulation {
		t.Errorf("YRI Ne (%d) should be > CHB Ne (%d)",
			yriParams.EffectivePopulation, chbParams.EffectivePopulation)
	}
}

func TestGetPopulationParameters_SuperPopulation(t *testing.T) {
	sim := population.NewDemographicSimulator("/tmp/test-project")

	tests := []struct {
		pop      population.Population
		superPop string
	}{
		{population.YRI, "African"},
		{population.CEU, "European"},
		{population.CHB, "East Asian"},
		{population.MXL, "Americas"},
	}

	for _, tt := range tests {
		params, err := sim.GetPopulationParameters(tt.pop)
		if err != nil {
			t.Errorf("unexpected error for %s: %v", tt.pop, err)
			continue
		}
		if params.SuperPopulation != tt.superPop {
			t.Errorf("%s SuperPopulation = %q, want %q",
				tt.pop, params.SuperPopulation, tt.superPop)
		}
	}
}
