//go:build !release
// +build !release

package main

import (
	"math"
	"testing"
	"time"

	phantomtesting "github.com/scttfrdmn/phantom-science/internal/testing"
)

func TestCalculateMean_Empty(t *testing.T) {
	if got := calculateMean([]float64{}); got != 0.0 {
		t.Errorf("calculateMean([]) = %f, want 0.0", got)
	}
}

func TestCalculateMean_Single(t *testing.T) {
	if got := calculateMean([]float64{5.0}); got != 5.0 {
		t.Errorf("calculateMean([5.0]) = %f, want 5.0", got)
	}
}

func TestCalculateMean_Multiple(t *testing.T) {
	got := calculateMean([]float64{1.0, 2.0, 3.0, 4.0})
	if math.Abs(got-2.5) > 1e-9 {
		t.Errorf("calculateMean([1,2,3,4]) = %f, want 2.5", got)
	}
}

func TestCalculateStdDev_Empty(t *testing.T) {
	if got := calculateStdDev([]float64{}, 0.0); got != 0.0 {
		t.Errorf("calculateStdDev([], 0) = %f, want 0.0", got)
	}
}

func TestCalculateStdDev_AllSame(t *testing.T) {
	values := []float64{5.0, 5.0, 5.0}
	got := calculateStdDev(values, 5.0)
	if got != 0.0 {
		t.Errorf("calculateStdDev([5,5,5], 5.0) = %f, want 0.0", got)
	}
}

func TestParsePlatformList_All(t *testing.T) {
	platforms := parsePlatformList("all")
	if len(platforms) == 0 {
		t.Error("parsePlatformList(all) returned empty list")
	}
}

func TestParsePlatformList_Illumina(t *testing.T) {
	platforms := parsePlatformList("illumina")
	if len(platforms) != 3 {
		t.Errorf("parsePlatformList(illumina) = %d platforms, want 3", len(platforms))
	}
}

func TestParsePlatformList_PacBio(t *testing.T) {
	platforms := parsePlatformList("pacbio")
	if len(platforms) != 1 || platforms[0] != phantomtesting.PacBioRevio {
		t.Errorf("parsePlatformList(pacbio) = %v", platforms)
	}
}

func TestParsePlatformList_NoDuplicates(t *testing.T) {
	platforms := parsePlatformList("all,novaseq")
	seen := make(map[phantomtesting.SequencingPlatform]bool)
	for _, p := range platforms {
		if seen[p] {
			t.Errorf("duplicate platform %q in result", p)
		}
		seen[p] = true
	}
}

func TestParsePopulationList_All(t *testing.T) {
	pops := parsePopulationList("all")
	if len(pops) == 0 {
		t.Error("parsePopulationList(all) returned empty list")
	}
}

func TestParsePopulationList_African(t *testing.T) {
	pops := parsePopulationList("african")
	if len(pops) != 3 {
		t.Errorf("parsePopulationList(african) = %d populations, want 3", len(pops))
	}
}

func TestParsePopulationList_SpecificCode(t *testing.T) {
	pops := parsePopulationList("yri")
	if len(pops) != 1 || pops[0] != "YRI" {
		t.Errorf("parsePopulationList(yri) = %v, want [YRI]", pops)
	}
}

func TestParsePopulationList_NoDuplicates(t *testing.T) {
	pops := parsePopulationList("all,african")
	seen := make(map[string]bool)
	for _, p := range pops {
		if seen[p] {
			t.Errorf("duplicate population %q in result", p)
		}
		seen[p] = true
	}
}

func TestCalculateConsistencyScore_SinglePlatform(t *testing.T) {
	scores := map[phantomtesting.SequencingPlatform][]float64{
		phantomtesting.IlluminaNovaSeq6000: {80.0, 82.0, 81.0},
	}
	got := calculateConsistencyScore(scores)
	if got != 100.0 {
		t.Errorf("calculateConsistencyScore(single platform) = %f, want 100.0", got)
	}
}

func TestCalculateConsistencyScore_Empty(t *testing.T) {
	got := calculateConsistencyScore(map[phantomtesting.SequencingPlatform][]float64{})
	if got != 100.0 {
		t.Errorf("calculateConsistencyScore(empty) = %f, want 100.0", got)
	}
}

func TestCalculateAuthenticityVariance_Empty(t *testing.T) {
	got := calculateAuthenticityVariance(map[phantomtesting.SequencingPlatform][]float64{})
	if got != 0.0 {
		t.Errorf("calculateAuthenticityVariance(empty) = %f, want 0.0", got)
	}
}

func TestCalculateAuthenticityVariance_SingleScore(t *testing.T) {
	scores := map[phantomtesting.SequencingPlatform][]float64{
		phantomtesting.IlluminaNovaSeq6000: {80.0},
	}
	got := calculateAuthenticityVariance(scores)
	if got != 0.0 {
		t.Errorf("calculateAuthenticityVariance(single score) = %f, want 0.0", got)
	}
}

func TestValidateArgs_EmptyOutputDir(t *testing.T) {
	args := &CommandLineArgs{
		outputDir:   "",
		platforms:   "all",
		populations: "all",
	}
	if err := validateArgs(args); err == nil {
		t.Error("validateArgs with empty outputDir should return an error")
	}
}

func TestValidateArgs_EmptyPlatforms(t *testing.T) {
	args := &CommandLineArgs{
		outputDir:   "/tmp",
		platforms:   "invalid_platform",
		populations: "all",
	}
	if err := validateArgs(args); err == nil {
		t.Error("validateArgs with invalid platform should return an error")
	}
}

func TestValidateArgs_Valid(t *testing.T) {
	args := &CommandLineArgs{
		outputDir:   "/tmp",
		platforms:   "all",
		populations: "all",
	}
	if err := validateArgs(args); err != nil {
		t.Errorf("validateArgs with valid args should not return an error, got: %v", err)
	}
}

func TestCalculateOverallScore_Empty(t *testing.T) {
	results := &TestResults{
		PlatformResults: map[phantomtesting.SequencingPlatform]*PlatformTestResult{},
	}
	got := calculateOverallScore(results)
	if got != 0.0 {
		t.Errorf("calculateOverallScore(empty) = %f, want 0.0", got)
	}
}

func TestCalculateOverallScore_WithResult(t *testing.T) {
	results := &TestResults{
		PlatformResults: map[phantomtesting.SequencingPlatform]*PlatformTestResult{
			phantomtesting.IlluminaNovaSeq6000: {
				SuccessRate:         1.0,
				AverageAuthenticity: 80.0,
			},
		},
	}
	got := calculateOverallScore(results)
	if got <= 0 {
		t.Errorf("calculateOverallScore with one result = %f, want > 0", got)
	}
}

func TestCalculateConsistencyScore_TwoPlatforms(t *testing.T) {
	scores := map[phantomtesting.SequencingPlatform][]float64{
		phantomtesting.IlluminaNovaSeq6000: {80.0, 81.0, 80.5},
		phantomtesting.IlluminaNextSeq2000: {80.0, 81.0, 80.5},
	}
	got := calculateConsistencyScore(scores)
	if got < 90.0 {
		t.Errorf("calculateConsistencyScore with similar scores = %f, want >= 90.0", got)
	}
}

func TestCalculateAuthenticityVariance_MultipleScores(t *testing.T) {
	scores := map[phantomtesting.SequencingPlatform][]float64{
		phantomtesting.IlluminaNovaSeq6000: {80.0, 82.0},
		phantomtesting.IlluminaNextSeq2000: {78.0, 81.0},
	}
	got := calculateAuthenticityVariance(scores)
	if got < 0 {
		t.Errorf("calculateAuthenticityVariance = %f, want >= 0", got)
	}
}

func TestPrintBanner_NoCrash(t *testing.T) {
	printBanner()
}

func TestPrintUsage_NoCrash(t *testing.T) {
	printUsage()
}

func makeTestResults(overallScore float64) *TestResults {
	return &TestResults{
		Platforms:       []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000},
		Populations:     []string{"YRI", "CEU"},
		PlatformResults: map[phantomtesting.SequencingPlatform]*PlatformTestResult{
			phantomtesting.IlluminaNovaSeq6000: {
				SuccessRate:         0.95,
				AverageAuthenticity: 80.0,
				TestDuration:        time.Second * 5,
				ErrorCount:          1,
			},
		},
		BenchmarkResults: map[phantomtesting.SequencingPlatform]*phantomtesting.PlatformBenchmark{},
		OverallScore:     overallScore,
		TotalDuration:    time.Second * 30,
	}
}

func TestGenerateJSONReport_NoCrash(t *testing.T) {
	results := makeTestResults(85.0)
	args := &CommandLineArgs{outputDir: "/tmp/e2e_test"}
	generateJSONReport(results, args)
}

func TestGenerateHTMLReport_NoCrash(t *testing.T) {
	results := makeTestResults(85.0)
	args := &CommandLineArgs{outputDir: "/tmp/e2e_test"}
	generateHTMLReport(results, args)
}

func TestGenerateCSVSummary_NoCrash(t *testing.T) {
	results := makeTestResults(85.0)
	args := &CommandLineArgs{outputDir: "/tmp/e2e_test"}
	generateCSVSummary(results, args)
}

func TestGenerateReports_NoCrash(t *testing.T) {
	results := makeTestResults(85.0)
	args := &CommandLineArgs{outputDir: "/tmp/e2e_test"}
	generateReports(results, args)
}

func TestPrintTestSummary_Excellent(t *testing.T) {
	results := makeTestResults(95.0) // >= 90 → EXCELLENT
	args := &CommandLineArgs{outputDir: "/tmp", report: true}
	printTestSummary(results, args)
}

func TestPrintTestSummary_Good(t *testing.T) {
	results := makeTestResults(78.0) // >= 75 → GOOD
	args := &CommandLineArgs{outputDir: "/tmp"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_Fair(t *testing.T) {
	results := makeTestResults(65.0) // >= 60 → FAIR
	args := &CommandLineArgs{outputDir: "/tmp"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_Poor(t *testing.T) {
	results := makeTestResults(45.0) // < 60 → POOR
	args := &CommandLineArgs{outputDir: "/tmp"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_WithAuthenticityResults(t *testing.T) {
	results := makeTestResults(85.0)
	results.AuthenticityResults = &AuthenticityResult{
		OverallAccuracy: 0.503, // ~50.3% → EXCELLENT
	}
	results.CrossPlatformResults = &CrossPlatformResult{
		ConsistencyScore:     95.0,
		AuthenticityVariance: 0.001,
	}
	args := &CommandLineArgs{outputDir: "/tmp"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_AuthenticityGood(t *testing.T) {
	results := makeTestResults(75.0)
	results.AuthenticityResults = &AuthenticityResult{
		OverallAccuracy: 0.55, // abs = 0.047 < 0.1 → GOOD
	}
	args := &CommandLineArgs{}
	printTestSummary(results, args)
}

func TestPrintTestSummary_AuthenticityNeedsImprovement(t *testing.T) {
	results := makeTestResults(50.0)
	results.AuthenticityResults = &AuthenticityResult{
		OverallAccuracy: 0.7, // abs = 0.197 > 0.1 → NEEDS IMPROVEMENT
	}
	args := &CommandLineArgs{}
	printTestSummary(results, args)
}

// Ensure math and time imports are used
var _ = math.Abs
var _ = time.Second
