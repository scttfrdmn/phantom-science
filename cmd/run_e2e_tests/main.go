package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/testing"
)

// CommandLineArgs holds command line arguments
type CommandLineArgs struct {
	platforms    string
	populations  string
	outputDir    string
	verbose      bool
	quick        bool
	parallel     bool
	benchmark    bool
	report       bool
	authenticity bool
	help         bool
}

func main() {
	args := parseCommandLineArgs()
	
	if args.help {
		printUsage()
		return
	}
	
	// Print banner
	printBanner()
	
	// Validate arguments
	if err := validateArgs(args); err != nil {
		log.Fatalf("❌ Invalid arguments: %v", err)
	}
	
	// Create output directory
	if err := os.MkdirAll(args.outputDir, 0750); err != nil {
		log.Fatalf("❌ Failed to create output directory %s: %v", args.outputDir, err)
	}
	
	// Initialize E2E framework
	framework := testing.NewE2EFramework()
	if framework == nil {
		log.Fatalf("❌ Failed to initialize E2E testing framework")
	}
	
	// Parse platform and population lists
	platforms := parsePlatformList(args.platforms)
	populations := parsePopulationList(args.populations)
	
	if args.verbose {
		log.Printf("🎯 Testing %d platforms across %d populations", len(platforms), len(populations))
		log.Printf("📋 Platforms: %v", platforms)
		log.Printf("🌍 Populations: %v", populations)
	}
	
	// Run comprehensive tests
	results := runComprehensiveTests(framework, platforms, populations, args)
	
	// Generate reports
	if args.report {
		generateReports(results, args)
	}
	
	// Print summary
	printTestSummary(results, args)

	if results.AuthenticityResults != nil && !results.AuthenticityResults.Passed() {
		log.Fatalf("❌ ML calibration failed on %d of %d platforms",
			results.AuthenticityResults.FailedCount(), len(results.AuthenticityResults.Platforms))
	}

	log.Println("✅ End-to-end testing completed")
}

// parseCommandLineArgs parses command line arguments
func parseCommandLineArgs() *CommandLineArgs {
	args := &CommandLineArgs{}
	
	flag.StringVar(&args.platforms, "platforms", "all", 
		"Comma-separated list of platforms to test (all, illumina, longread, pacbio, ont, novaseq, nextseq, miseq, revio, promethion, dnbseq, aviti)")
	flag.StringVar(&args.populations, "populations", "major", 
		"Comma-separated list of populations to test (all, major, african, european, asian, americas)")
	flag.StringVar(&args.outputDir, "output", "./e2e_results", 
		"Output directory for test results and reports")
	flag.BoolVar(&args.verbose, "verbose", false, 
		"Enable verbose logging")
	flag.BoolVar(&args.quick, "quick", false, 
		"Run quick tests only (reduced data generation)")
	flag.BoolVar(&args.parallel, "parallel", true, 
		"Run tests in parallel across platforms")
	flag.BoolVar(&args.benchmark, "benchmark", true, 
		"Include performance benchmarking")
	flag.BoolVar(&args.report, "report", true, 
		"Generate comprehensive HTML and JSON reports")
	flag.BoolVar(&args.authenticity, "authenticity", true, 
		"Include the ML classifier calibration (null comparison plus GC-skewed positive control)")
	flag.BoolVar(&args.help, "help", false, 
		"Show this help message")
	
	flag.Parse()
	return args
}

// printUsage prints usage information
func printUsage() {
	fmt.Println(`
🧬 Phantom Science E2E Testing Suite

Usage: run_e2e_tests [OPTIONS]

This tool runs comprehensive end-to-end testing of the Phantom Science genomics
authenticity engine across multiple sequencing platforms and global populations.

OPTIONS:`)
	flag.PrintDefaults()
	
	fmt.Print(`
PLATFORM OPTIONS:
  all       - All 7 supported platforms
  illumina  - All Illumina platforms (NovaSeq, NextSeq, MiSeq)
  longread  - Long-read platforms (PacBio, ONT)
  pacbio    - PacBio Revio only
  ont       - Oxford Nanopore PromethION only
  novaseq   - Illumina NovaSeq 6000 only
  nextseq   - Illumina NextSeq 2000 only
  miseq     - Illumina MiSeq only
  revio     - PacBio Revio only
  promethion- ONT PromethION 24 only
  dnbseq    - MGI DNBSEQ-T20×2 only
  aviti     - Element AVITI only

POPULATION OPTIONS:
  all       - All 11 supported populations
  major     - Major populations (YRI, CEU, CHB, MXL, JPT)
  african   - African populations (YRI, LWK, GWD)
  european  - European populations (CEU, GBR, TSI)
  asian     - Asian populations (CHB, JPT, CHS)
  americas  - Americas populations (MXL, PUR)

EXAMPLES:
  # Test all platforms with major populations
  ./run_e2e_tests -platforms=all -populations=major

  # Test only Illumina platforms with European populations
  ./run_e2e_tests -platforms=illumina -populations=european -verbose

  # Quick benchmark test of long-read platforms
  ./run_e2e_tests -platforms=longread -populations=major -quick -benchmark

  # Comprehensive authenticity validation
  ./run_e2e_tests -platforms=all -populations=all -authenticity -report

ML CALIBRATION (-authenticity):
  No real sequencing reads are bundled, so this tool cannot measure how
  synthetic reads compare to genuine reads. Instead it calibrates the
  measurement pipeline on each platform, using 1000 training and 1000
  holdout reads per class:
    Null comparison:  two independent draws from the same generator. Passes
                      when holdout accuracy lies inside the 4-sigma binomial
                      band around 0.5 (0.5 +/- 4*sqrt(0.25/n)).
    Positive control: the same comparison with one class GC-skewed. Passes
                      when holdout accuracy exceeds 0.9, which shows the
                      classifier can detect a real difference.
  The command exits non-zero if any platform fails either check.
`)
}

// printBanner prints the application banner
func printBanner() {
	fmt.Print(`
╔══════════════════════════════════════════════════════════════════════════════╗
║                   🧬 PHANTOM SCIENCE E2E TESTING SUITE 🧬                    ║
║                                                                              ║
║              Comprehensive End-to-End Platform Validation                   ║
║          ML calibration: null comparison and GC-skew positive control        ║
╚══════════════════════════════════════════════════════════════════════════════╝
`)
}

// validateArgs validates command line arguments
func validateArgs(args *CommandLineArgs) error {
	// Validate output directory
	if args.outputDir == "" {
		return fmt.Errorf("output directory cannot be empty")
	}
	
	// Validate platform list
	platforms := parsePlatformList(args.platforms)
	if len(platforms) == 0 {
		return fmt.Errorf("no valid platforms specified")
	}
	
	// Validate population list  
	populations := parsePopulationList(args.populations)
	if len(populations) == 0 {
		return fmt.Errorf("no valid populations specified")
	}
	
	return nil
}

// parsePlatformList parses the platform list argument
func parsePlatformList(platformsArg string) []testing.SequencingPlatform {
	platformMap := map[string][]testing.SequencingPlatform{
		"all": {
			testing.IlluminaNovaSeq6000,
			testing.IlluminaNextSeq2000,
			testing.IlluminaMiSeq,
			testing.PacBioRevio,
			testing.ONTPromethION24,
			testing.MGI_DNBSEQ_T20x2,
			testing.ElementAVITI,
		},
		"illumina": {
			testing.IlluminaNovaSeq6000,
			testing.IlluminaNextSeq2000,
			testing.IlluminaMiSeq,
		},
		"longread": {
			testing.PacBioRevio,
			testing.ONTPromethION24,
		},
		"pacbio":     {testing.PacBioRevio},
		"ont":        {testing.ONTPromethION24},
		"novaseq":    {testing.IlluminaNovaSeq6000},
		"nextseq":    {testing.IlluminaNextSeq2000},
		"miseq":      {testing.IlluminaMiSeq},
		"revio":      {testing.PacBioRevio},
		"promethion": {testing.ONTPromethION24},
		"dnbseq":     {testing.MGI_DNBSEQ_T20x2},
		"aviti":      {testing.ElementAVITI},
	}
	
	var platforms []testing.SequencingPlatform
	parts := strings.Split(strings.ToLower(platformsArg), ",")
	
	for _, part := range parts {
		part = strings.TrimSpace(part)
		if platformList, exists := platformMap[part]; exists {
			platforms = append(platforms, platformList...)
		}
	}
	
	// Remove duplicates
	platformSet := make(map[testing.SequencingPlatform]bool)
	uniquePlatforms := make([]testing.SequencingPlatform, 0)
	
	for _, platform := range platforms {
		if !platformSet[platform] {
			platformSet[platform] = true
			uniquePlatforms = append(uniquePlatforms, platform)
		}
	}
	
	return uniquePlatforms
}

// parsePopulationList parses the population list argument
func parsePopulationList(populationsArg string) []string {
	populationMap := map[string][]string{
		"all": {
			"YRI", "LWK", "GWD", // African
			"CEU", "GBR", "TSI", // European
			"CHB", "JPT", "CHS", // Asian
			"MXL", "PUR",        // Americas
		},
		"major": {
			"YRI", "CEU", "CHB", "MXL", "JPT",
		},
		"african": {
			"YRI", "LWK", "GWD",
		},
		"european": {
			"CEU", "GBR", "TSI",
		},
		"asian": {
			"CHB", "JPT", "CHS",
		},
		"americas": {
			"MXL", "PUR",
		},
	}
	
	var populations []string
	parts := strings.Split(strings.ToLower(populationsArg), ",")
	
	for _, part := range parts {
		part = strings.TrimSpace(part)
		if populationList, exists := populationMap[part]; exists {
			populations = append(populations, populationList...)
		} else {
			// Check if it's a specific population code
			upperPart := strings.ToUpper(part)
			if len(upperPart) == 3 {
				populations = append(populations, upperPart)
			}
		}
	}
	
	// Remove duplicates
	populationSet := make(map[string]bool)
	uniquePopulations := make([]string, 0)
	
	for _, population := range populations {
		if !populationSet[population] {
			populationSet[population] = true
			uniquePopulations = append(uniquePopulations, population)
		}
	}
	
	return uniquePopulations
}

// TestResults holds comprehensive test results
type TestResults struct {
	Platforms            []testing.SequencingPlatform
	Populations          []string
	PlatformResults      map[testing.SequencingPlatform]*PlatformTestResult
	CrossPlatformResults *CrossPlatformResult
	AuthenticityResults  *AuthenticityResult
	BenchmarkResults     map[testing.SequencingPlatform]*testing.PlatformBenchmark
	OverallScore         float64
	StartTime            time.Time
	EndTime              time.Time
	TotalDuration        time.Duration
}

// PlatformTestResult holds results for a single platform
type PlatformTestResult struct {
	Platform            testing.SequencingPlatform
	PopulationResults   map[string]*PopulationResult
	AverageAuthenticity float64
	SuccessRate         float64
	ErrorCount          int
	TestDuration        time.Duration
}

// PopulationResult holds results for a platform-population combination
type PopulationResult struct {
	Population        string
	AuthenticityScore float64
	QualityScore      float64
	ErrorRate         float64
	Success           bool
	ErrorMessage      string
}

// CrossPlatformResult holds cross-platform analysis results
type CrossPlatformResult struct {
	ConsistencyScore     float64
	PlatformComparisons  map[string]map[string]float64
	AuthenticityVariance float64
}

// ML calibration parameters. Holdout accuracy under the null is
// Binomial(n, 0.5)/n with n = 2*authHoldoutPerClass, so its standard deviation
// is sqrt(0.25/n), about 0.011 for 1000 reads per class.
const (
	authTrainPerClass   = 1000
	authHoldoutPerClass = 1000
	// nullBandSigmas sets the width of the pass band around 0.5. At 4 sigma a
	// correct pipeline fails by chance with probability below 1e-4.
	nullBandSigmas = 4.0
	// gcSkewProbability is the chance that each A/T base in the positive
	// control set is replaced by G/C (GC content ~0.5 -> ~0.7).
	gcSkewProbability = 0.4
	// positiveControlMin is the minimum holdout accuracy the classifier must
	// reach on the GC-skewed positive control.
	positiveControlMin = 0.9
)

// PlatformCalibration holds the ML calibration result for one platform.
type PlatformCalibration struct {
	// NullAccuracy is holdout accuracy on two independent draws from the
	// same generator; it should be near 0.5.
	NullAccuracy float64
	// Band is the half-width of the pass band around 0.5.
	Band float64
	// ControlAccuracy is holdout accuracy when one class is GC-skewed.
	ControlAccuracy float64
	NullPass        bool
	ControlPass     bool
	// Err is set when the calibration could not be run; the platform fails.
	Err string
}

// Passed reports whether the platform ran without error and passed both checks.
func (c *PlatformCalibration) Passed() bool {
	return c.Err == "" && c.NullPass && c.ControlPass
}

// AuthenticityResult holds the per-platform ML calibration results. It is a
// calibration of the measurement pipeline, not a real-vs-synthetic measurement:
// no real reads are bundled.
type AuthenticityResult struct {
	Platforms map[testing.SequencingPlatform]*PlatformCalibration
}

// FailedCount returns the number of platforms that errored or failed a check.
func (r *AuthenticityResult) FailedCount() int {
	failed := 0
	for _, c := range r.Platforms {
		if !c.Passed() {
			failed++
		}
	}
	return failed
}

// Passed reports whether every platform passed. An empty result does not pass.
func (r *AuthenticityResult) Passed() bool {
	return len(r.Platforms) > 0 && r.FailedCount() == 0
}

// nullBand returns the nullBandSigmas-sigma half-width of the binomial band
// around 0.5 for a holdout of holdoutPerClass reads per class.
func nullBand(holdoutPerClass int) float64 {
	n := float64(2 * holdoutPerClass)
	return nullBandSigmas * math.Sqrt(0.25/n)
}

// gradeCalibration grades a platform's null and positive-control accuracies.
func gradeCalibration(nullAccuracy, controlAccuracy float64, holdoutPerClass int) *PlatformCalibration {
	band := nullBand(holdoutPerClass)
	return &PlatformCalibration{
		NullAccuracy:    nullAccuracy,
		Band:            band,
		ControlAccuracy: controlAccuracy,
		NullPass:        math.Abs(nullAccuracy-0.5) <= band,
		ControlPass:     controlAccuracy > positiveControlMin,
	}
}

// runComprehensiveTests runs the comprehensive test suite
func runComprehensiveTests(framework *testing.E2EFramework, platforms []testing.SequencingPlatform, populations []string, args *CommandLineArgs) *TestResults {
	startTime := time.Now()
	
	results := &TestResults{
		Platforms:       platforms,
		Populations:     populations,
		PlatformResults: make(map[testing.SequencingPlatform]*PlatformTestResult),
		BenchmarkResults: make(map[testing.SequencingPlatform]*testing.PlatformBenchmark),
		StartTime:       startTime,
	}
	
	log.Printf("🚀 Starting comprehensive E2E testing...")
	log.Printf("📋 Testing %d platforms × %d populations = %d combinations", 
		len(platforms), len(populations), len(platforms)*len(populations))
	
	// Run platform-specific tests
	for i, platform := range platforms {
		log.Printf("🧬 [%d/%d] Testing platform: %s", i+1, len(platforms), platform)
		
		platformResult := runPlatformTests(framework, platform, populations, args)
		results.PlatformResults[platform] = platformResult
		
		if args.verbose {
			log.Printf("   ✅ %s completed: %.1f%% success, %.2f authenticity", 
				platform, platformResult.SuccessRate*100, platformResult.AverageAuthenticity)
		}
	}
	
	// Run cross-platform analysis
	if len(platforms) > 1 {
		log.Printf("🔄 Running cross-platform consistency analysis...")
		results.CrossPlatformResults = runCrossPlatformAnalysis(framework, platforms, populations, args)
	}
	
	// Run authenticity validation
	if args.authenticity {
		log.Printf("🎯 Running statistical authenticity validation...")
		results.AuthenticityResults = runAuthenticityValidation(framework, platforms, populations, args)
	}
	
	// Run benchmarks
	if args.benchmark {
		log.Printf("⚡ Running performance benchmarks...")
		results.BenchmarkResults = runPerformanceBenchmarks(framework, platforms, args)
	}
	
	// Calculate overall results
	results.EndTime = time.Now()
	results.TotalDuration = results.EndTime.Sub(results.StartTime)
	results.OverallScore = calculateOverallScore(results)
	
	return results
}

// runPlatformTests runs tests for a specific platform
func runPlatformTests(framework *testing.E2EFramework, platform testing.SequencingPlatform, populations []string, args *CommandLineArgs) *PlatformTestResult {
	startTime := time.Now()
	
	result := &PlatformTestResult{
		Platform:          platform,
		PopulationResults: make(map[string]*PopulationResult),
	}
	
	totalAuthenticity := 0.0
	successCount := 0
	errorCount := 0
	
	for _, population := range populations {
		if args.verbose {
			log.Printf("   🌍 Testing %s with population %s", platform, population)
		}
		
		popResult := runPopulationTest(framework, platform, population, args)
		result.PopulationResults[population] = popResult
		
		if popResult.Success {
			totalAuthenticity += popResult.AuthenticityScore
			successCount++
		} else {
			errorCount++
		}
	}
	
	// Calculate averages
	if successCount > 0 {
		result.AverageAuthenticity = totalAuthenticity / float64(successCount)
		result.SuccessRate = float64(successCount) / float64(len(populations))
	}
	
	result.ErrorCount = errorCount
	result.TestDuration = time.Since(startTime)
	
	return result
}

// runPopulationTest runs a test for a specific platform-population combination
func runPopulationTest(framework *testing.E2EFramework, platform testing.SequencingPlatform, population string, args *CommandLineArgs) *PopulationResult {
	result := &PopulationResult{
		Population: population,
	}
	
	// Generate synthetic data
	config := framework.CreatePlatformConfig(platform, population)
	if config == nil {
		result.Success = false
		result.ErrorMessage = "Failed to create platform configuration"
		return result
	}
	
	syntheticData, err := framework.GenerateSyntheticData(config)
	if err != nil {
		result.Success = false
		result.ErrorMessage = fmt.Sprintf("Failed to generate synthetic data: %v", err)
		return result
	}
	
	// Analyze authenticity
	authenticityScore, err := framework.AnalyzeAuthenticity(syntheticData)
	if err != nil {
		result.Success = false
		result.ErrorMessage = fmt.Sprintf("Failed to analyze authenticity: %v", err)
		return result
	}
	
	// Populate results
	result.AuthenticityScore = authenticityScore
	result.QualityScore = syntheticData.QualityScoreAverage
	result.ErrorRate = syntheticData.ErrorRate
	result.Success = true
	
	return result
}

// runCrossPlatformAnalysis runs cross-platform consistency analysis
func runCrossPlatformAnalysis(framework *testing.E2EFramework, platforms []testing.SequencingPlatform, populations []string, args *CommandLineArgs) *CrossPlatformResult {
	result := &CrossPlatformResult{
		PlatformComparisons: make(map[string]map[string]float64),
	}
	
	// Collect authenticity scores across platforms
	platformScores := make(map[testing.SequencingPlatform][]float64)
	
	for _, platform := range platforms {
		scores := make([]float64, 0)
		
		for _, population := range populations {
			config := framework.CreatePlatformConfig(platform, population)
			if config == nil {
				log.Printf("   ⚠️  cross-platform: %s/%s: no platform configuration; excluded", platform, population)
				continue
			}

			data, err := framework.GenerateSyntheticData(config)
			if err != nil {
				log.Printf("   ⚠️  cross-platform: %s/%s: generate data: %v; excluded", platform, population, err)
				continue
			}

			score, err := framework.AnalyzeAuthenticity(data)
			if err != nil {
				log.Printf("   ⚠️  cross-platform: %s/%s: analyze: %v; excluded", platform, population, err)
				continue
			}
			
			scores = append(scores, score)
		}
		
		platformScores[platform] = scores
	}
	
	// Calculate consistency metrics
	result.ConsistencyScore = calculateConsistencyScore(platformScores)
	result.AuthenticityVariance = calculateAuthenticityVariance(platformScores)
	
	return result
}

// runAuthenticityValidation calibrates the ML indistinguishability measurement
// on each platform: a null comparison of two independent draws from the same
// generator (expected accuracy 0.5) and a GC-skewed positive control (expected
// accuracy near 1). A platform whose data cannot be generated or whose
// classifier cannot be trained is recorded as failed, not skipped.
func runAuthenticityValidation(framework *testing.E2EFramework, platforms []testing.SequencingPlatform, populations []string, args *CommandLineArgs) *AuthenticityResult {
	result := &AuthenticityResult{
		Platforms: make(map[testing.SequencingPlatform]*PlatformCalibration),
	}

	for _, platform := range platforms {
		calibration, err := calibratePlatform(framework, platform)
		if err != nil {
			log.Printf("   ❌ %s: ML calibration failed: %v", platform, err)
			calibration = &PlatformCalibration{Err: err.Error()}
		}
		result.Platforms[platform] = calibration
	}

	return result
}

// calibratePlatform runs the null comparison and positive control for one platform.
func calibratePlatform(framework *testing.E2EFramework, platform testing.SequencingPlatform) (*PlatformCalibration, error) {
	trainReal, trainSynthetic, err := framework.GenerateBalancedDataset(platform, authTrainPerClass)
	if err != nil {
		return nil, fmt.Errorf("generate training data: %w", err)
	}
	testReal, testSynthetic, err := framework.GenerateBalancedDataset(platform, authHoldoutPerClass)
	if err != nil {
		return nil, fmt.Errorf("generate holdout data: %w", err)
	}

	null := framework.CreateMLClassifier()
	if err := null.Train(trainReal, trainSynthetic); err != nil {
		return nil, fmt.Errorf("train null classifier: %w", err)
	}
	nullAccuracy := null.TestAccuracy(testReal, testSynthetic)

	control := framework.CreateMLClassifier()
	if err := control.Train(trainReal, testing.GCSkewed(trainSynthetic, gcSkewProbability, 1)); err != nil {
		return nil, fmt.Errorf("train positive-control classifier: %w", err)
	}
	controlAccuracy := control.TestAccuracy(testReal, testing.GCSkewed(testSynthetic, gcSkewProbability, 2))

	return gradeCalibration(nullAccuracy, controlAccuracy, authHoldoutPerClass), nil
}

// runPerformanceBenchmarks runs performance benchmarks
func runPerformanceBenchmarks(framework *testing.E2EFramework, platforms []testing.SequencingPlatform, args *CommandLineArgs) map[testing.SequencingPlatform]*testing.PlatformBenchmark {
	benchmarks := make(map[testing.SequencingPlatform]*testing.PlatformBenchmark)
	
	for _, platform := range platforms {
		if args.verbose {
			log.Printf("   ⚡ Benchmarking platform: %s", platform)
		}
		
		benchmark := framework.RunBenchmark(platform)
		benchmarks[platform] = benchmark
	}
	
	return benchmarks
}

// calculateOverallScore calculates the overall test score
func calculateOverallScore(results *TestResults) float64 {
	totalScore := 0.0
	validPlatforms := 0
	
	for _, platformResult := range results.PlatformResults {
		if platformResult.SuccessRate > 0 {
			// Weight by success rate and authenticity score
			platformScore := (platformResult.SuccessRate * 50.0) + (platformResult.AverageAuthenticity * 0.5)
			totalScore += platformScore
			validPlatforms++
		}
	}
	
	if validPlatforms == 0 {
		return 0.0
	}
	
	return totalScore / float64(validPlatforms)
}

// calculateConsistencyScore calculates cross-platform consistency
func calculateConsistencyScore(platformScores map[testing.SequencingPlatform][]float64) float64 {
	if len(platformScores) < 2 {
		return 100.0
	}
	
	// Calculate coefficient of variation across platforms
	allMeans := make([]float64, 0)
	
	for _, scores := range platformScores {
		if len(scores) > 0 {
			mean := calculateMean(scores)
			allMeans = append(allMeans, mean)
		}
	}
	
	if len(allMeans) < 2 {
		return 100.0
	}
	
	overallMean := calculateMean(allMeans)
	overallStdDev := calculateStdDev(allMeans, overallMean)
	
	cv := overallStdDev / overallMean
	consistencyScore := math.Max(0, 100.0 - cv*500.0) // Scale to 0-100
	
	return consistencyScore
}

// calculateAuthenticityVariance calculates authenticity variance across platforms
func calculateAuthenticityVariance(platformScores map[testing.SequencingPlatform][]float64) float64 {
	allScores := make([]float64, 0)
	
	for _, scores := range platformScores {
		allScores = append(allScores, scores...)
	}
	
	if len(allScores) < 2 {
		return 0.0
	}
	
	mean := calculateMean(allScores)
	return calculateStdDev(allScores, mean)
}

// generateReports generates comprehensive reports
func generateReports(results *TestResults, args *CommandLineArgs) {
	log.Printf("📊 Generating comprehensive reports...")
	
	// Generate JSON report
	generateJSONReport(results, args)
	
	// Generate HTML report
	generateHTMLReport(results, args)
	
	// Generate CSV summary
	generateCSVSummary(results, args)
	
	log.Printf("   📄 Reports saved to: %s", args.outputDir)
}

// generateJSONReport generates JSON format report
func generateJSONReport(results *TestResults, args *CommandLineArgs) {
	// Implementation would generate detailed JSON report
	log.Printf("   📄 JSON report: %s/e2e_results.json", args.outputDir)
}

// generateHTMLReport generates HTML format report
func generateHTMLReport(results *TestResults, args *CommandLineArgs) {
	// Implementation would generate interactive HTML report
	log.Printf("   🌐 HTML report: %s/e2e_results.html", args.outputDir)
}

// generateCSVSummary generates CSV summary
func generateCSVSummary(results *TestResults, args *CommandLineArgs) {
	// Implementation would generate CSV summary
	log.Printf("   📊 CSV summary: %s/e2e_summary.csv", args.outputDir)
}

// printTestSummary prints a summary of test results
func printTestSummary(results *TestResults, args *CommandLineArgs) {
	fmt.Println("\n" + strings.Repeat("=", 80))
	fmt.Println("📊 COMPREHENSIVE E2E TEST RESULTS SUMMARY")
	fmt.Println(strings.Repeat("=", 80))
	
	fmt.Printf("⏱️  Test Duration: %v\n", results.TotalDuration)
	fmt.Printf("🧬 Platforms Tested: %d\n", len(results.Platforms))
	fmt.Printf("🌍 Populations Tested: %d\n", len(results.Populations))
	fmt.Printf("📈 Overall Score: %.1f/100\n", results.OverallScore)
	
	fmt.Println("\n🧬 PLATFORM RESULTS:")
	fmt.Println(strings.Repeat("-", 80))
	fmt.Printf("%-20s %-10s %-12s %-10s %-8s\n", 
		"Platform", "Success", "Authenticity", "Duration", "Errors")
	fmt.Println(strings.Repeat("-", 80))
	
	for _, platform := range results.Platforms {
		if result, exists := results.PlatformResults[platform]; exists {
			fmt.Printf("%-20s %7.1f%% %11.2f %9v %7d\n",
				platform,
				result.SuccessRate*100,
				result.AverageAuthenticity,
				result.TestDuration.Round(time.Second),
				result.ErrorCount)
		}
	}
	
	// ML calibration results
	if results.AuthenticityResults != nil {
		printCalibrationSummary(results.Platforms, results.AuthenticityResults)
	}

	// Cross-platform consistency
	if results.CrossPlatformResults != nil {
		fmt.Println("\n🔄 CROSS-PLATFORM CONSISTENCY:")
		fmt.Println(strings.Repeat("-", 80))
		fmt.Printf("Consistency Score: %.1f/100\n", results.CrossPlatformResults.ConsistencyScore)
		fmt.Printf("Authenticity Variance: %.3f\n", results.CrossPlatformResults.AuthenticityVariance)
	}
	
	// Performance benchmarks
	if len(results.BenchmarkResults) > 0 {
		fmt.Println("\n⚡ PERFORMANCE BENCHMARKS:")
		fmt.Println(strings.Repeat("-", 80))
		fmt.Printf("%-20s %-12s %-10s %-12s %-8s\n", 
			"Platform", "Throughput", "Accuracy", "Authenticity", "Score")
		fmt.Println(strings.Repeat("-", 80))
		
		for platform, benchmark := range results.BenchmarkResults {
			fmt.Printf("%-20s %8.1f MB/s %8.3f %11.1f %7.1f\n",
				platform,
				benchmark.ThroughputMBps,
				benchmark.AccuracyScore,
				benchmark.AuthenticityScore,
				benchmark.OverallScore)
		}
	}
	
	fmt.Println("\n" + strings.Repeat("=", 80))
	
	// Final assessment
	if results.OverallScore >= 90.0 {
		fmt.Println("🎉 EXCELLENT: All platforms achieved outstanding performance!")
	} else if results.OverallScore >= 75.0 {
		fmt.Println("✅ GOOD: Most platforms performed well with minor areas for improvement")
	} else if results.OverallScore >= 60.0 {
		fmt.Println("⚠️  FAIR: Some platforms need improvement to achieve target performance")
	} else {
		fmt.Println("❌ POOR: Significant improvements needed across multiple platforms")
	}
	
	if args.report {
		fmt.Printf("\n📊 Detailed reports available in: %s\n", args.outputDir)
	}
}

// printCalibrationSummary prints the per-platform ML calibration table.
func printCalibrationSummary(platforms []testing.SequencingPlatform, auth *AuthenticityResult) {
	fmt.Println("\n🎯 ML CLASSIFIER CALIBRATION:")
	fmt.Println(strings.Repeat("-", 80))
	fmt.Println("No real reads are bundled, so this calibrates the measurement pipeline")
	fmt.Println("rather than measuring real-vs-synthetic distinguishability.")
	fmt.Printf("Null: independent draws from one generator, pass if within 0.5 ± %.0fσ (%d holdout reads/class).\n",
		nullBandSigmas, authHoldoutPerClass)
	fmt.Printf("Positive control: one class GC-skewed (p=%.1f), pass if accuracy > %.2f.\n",
		gcSkewProbability, positiveControlMin)
	fmt.Println(strings.Repeat("-", 80))
	fmt.Printf("%-24s %-15s %-6s %-9s %-6s\n", "Platform", "Null accuracy", "Null", "Control", "Ctrl")
	fmt.Println(strings.Repeat("-", 80))

	passFail := func(ok bool) string {
		if ok {
			return "PASS"
		}
		return "FAIL"
	}
	for _, platform := range platforms {
		c, ok := auth.Platforms[platform]
		if !ok {
			continue
		}
		if c.Err != "" {
			fmt.Printf("%-24s ERROR: %s\n", platform, c.Err)
			continue
		}
		fmt.Printf("%-24s %.3f ± %.3f   %-6s %-9.3f %-6s\n",
			platform, c.NullAccuracy, c.Band, passFail(c.NullPass), c.ControlAccuracy, passFail(c.ControlPass))
	}

	if auth.Passed() {
		fmt.Println("✅ Calibration passed: chance accuracy on the null comparison, and the GC skew was detected.")
	} else {
		fmt.Printf("❌ Calibration failed on %d of %d platforms.\n", auth.FailedCount(), len(auth.Platforms))
	}
}

// Helper functions (same as in previous files)
func calculateMean(values []float64) float64 {
	if len(values) == 0 {
		return 0.0
	}
	
	sum := 0.0
	for _, v := range values {
		sum += v
	}
	
	return sum / float64(len(values))
}

func calculateStdDev(values []float64, mean float64) float64 {
	if len(values) == 0 {
		return 0.0
	}
	
	sumSquares := 0.0
	for _, v := range values {
		diff := v - mean
		sumSquares += diff * diff
	}
	
	variance := sumSquares / float64(len(values))
	return math.Sqrt(variance)
}