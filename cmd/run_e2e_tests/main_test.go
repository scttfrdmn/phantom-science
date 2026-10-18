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
		Platforms:   []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000},
		Populations: []string{"YRI", "CEU"},
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

func TestPrintTestSummary_WithCalibrationResults(t *testing.T) {
	results := makeTestResults(85.0)
	results.AuthenticityResults = &AuthenticityResult{
		Platforms: map[phantomtesting.SequencingPlatform]*PlatformCalibration{
			phantomtesting.IlluminaNovaSeq6000: gradeCalibration(0.5, 0.99, authHoldoutPerClass),
		},
	}
	results.CrossPlatformResults = &CrossPlatformResult{
		ConsistencyScore:     95.0,
		AuthenticityVariance: 0.001,
	}
	args := &CommandLineArgs{outputDir: "/tmp"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_CalibrationError(t *testing.T) {
	results := makeTestResults(50.0)
	results.AuthenticityResults = &AuthenticityResult{
		Platforms: map[phantomtesting.SequencingPlatform]*PlatformCalibration{
			phantomtesting.IlluminaNovaSeq6000: {Err: "generate training data: boom"},
		},
	}
	printTestSummary(results, &CommandLineArgs{})
}

// TestNullBand checks the band is 4 binomial standard deviations of the
// holdout accuracy under H0: 4*sqrt(0.25/2000) ~= 0.0447 for 1000 reads/class.
func TestNullBand(t *testing.T) {
	got := nullBand(1000)
	want := 4 * math.Sqrt(0.25/2000)
	if math.Abs(got-want) > 1e-12 {
		t.Errorf("nullBand(1000) = %f, want %f", got, want)
	}
	if nullBand(100) <= got {
		t.Errorf("band should shrink as the holdout grows: nullBand(100)=%f, nullBand(1000)=%f", nullBand(100), got)
	}
}

func TestGradeCalibration(t *testing.T) {
	band := nullBand(authHoldoutPerClass)
	tests := []struct {
		name           string
		null, control  float64
		wantNull       bool
		wantControl    bool
		wantPlatformOK bool
	}{
		{"chance and detected control", 0.5, 0.99, true, true, true},
		{"just inside upper band", 0.5 + band*0.99, 0.95, true, true, true},
		{"just inside lower band", 0.5 - band*0.99, 0.95, true, true, true},
		// 0.06 from chance is more than 5 sigma out at n=2000.
		{"outside band high", 0.56, 0.99, false, true, false},
		{"outside band low", 0.44, 0.99, false, true, false},
		// 0.55 is 0.05 from chance, beyond the 0.0447 band at n=2000.
		{"just outside band", 0.55, 0.99, false, true, false},
		{"control below threshold", 0.5, 0.85, true, false, false},
		{"control exactly at threshold fails", 0.5, positiveControlMin, true, false, false},
		{"both fail", 0.7, 0.5, false, false, false},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			c := gradeCalibration(tt.null, tt.control, authHoldoutPerClass)
			if c.NullPass != tt.wantNull {
				t.Errorf("NullPass = %v, want %v (null=%.4f, band=%.4f)", c.NullPass, tt.wantNull, tt.null, c.Band)
			}
			if c.ControlPass != tt.wantControl {
				t.Errorf("ControlPass = %v, want %v (control=%.4f)", c.ControlPass, tt.wantControl, tt.control)
			}
			if c.Passed() != tt.wantPlatformOK {
				t.Errorf("Passed() = %v, want %v", c.Passed(), tt.wantPlatformOK)
			}
		})
	}
}

func TestAuthenticityResult_Passed(t *testing.T) {
	ok := gradeCalibration(0.5, 0.99, authHoldoutPerClass)
	bad := gradeCalibration(0.6, 0.99, authHoldoutPerClass)
	errored := &PlatformCalibration{Err: "boom", NullPass: true, ControlPass: true}

	if (&AuthenticityResult{}).Passed() {
		t.Error("an empty result must not pass")
	}
	all := &AuthenticityResult{Platforms: map[phantomtesting.SequencingPlatform]*PlatformCalibration{
		phantomtesting.IlluminaNovaSeq6000: ok,
		phantomtesting.IlluminaMiSeq:       ok,
	}}
	if !all.Passed() || all.FailedCount() != 0 {
		t.Errorf("all-pass result: Passed=%v FailedCount=%d", all.Passed(), all.FailedCount())
	}
	mixed := &AuthenticityResult{Platforms: map[phantomtesting.SequencingPlatform]*PlatformCalibration{
		phantomtesting.IlluminaNovaSeq6000: ok,
		phantomtesting.IlluminaMiSeq:       bad,
		phantomtesting.PacBioRevio:         errored,
	}}
	if mixed.Passed() || mixed.FailedCount() != 2 {
		t.Errorf("mixed result: Passed=%v FailedCount=%d, want false/2", mixed.Passed(), mixed.FailedCount())
	}
}

// TestRunAuthenticityValidation runs the real calibration on one short-read
// platform: the null comparison must land in the band and the GC-skewed
// positive control must be detected.
func TestRunAuthenticityValidation(t *testing.T) {
	if testing.Short() {
		t.Skip("generates 4000 reads per platform")
	}
	framework := phantomtesting.NewE2EFramework()
	platform := phantomtesting.IlluminaNovaSeq6000
	res := runAuthenticityValidation(framework, []phantomtesting.SequencingPlatform{platform}, nil, &CommandLineArgs{})

	c, ok := res.Platforms[platform]
	if !ok {
		t.Fatalf("no calibration result for %s", platform)
	}
	if c.Err != "" {
		t.Fatalf("calibration error: %s", c.Err)
	}
	if !c.NullPass {
		t.Errorf("null accuracy %.3f outside 0.5 ± %.3f", c.NullAccuracy, c.Band)
	}
	if !c.ControlPass {
		t.Errorf("positive-control accuracy %.3f, want > %.2f", c.ControlAccuracy, positiveControlMin)
	}
	if !res.Passed() {
		t.Error("expected the calibration to pass")
	}
	t.Logf("null %.3f ± %.3f, control %.3f", c.NullAccuracy, c.Band, c.ControlAccuracy)
}

// Ensure math and time imports are used
var _ = math.Abs
var _ = time.Second
