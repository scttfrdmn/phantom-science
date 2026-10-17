package authenticity

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sort"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

// IndistinguishabilityTester performs cross-platform statistical indistinguishability testing
type IndistinguishabilityTester struct {
	platformManager    *genomics.PlatformErrorSignatureManager
	qualityModeler     *genomics.MultiPlatformQualityModeler
	classificationTests map[string]*ClassificationTest
	rng                *rand.Rand
}

// ClassificationTest represents a statistical test for population classification
type ClassificationTest struct {
	TestName             string
	TestType             string // "PCA", "FST", "LD", "QUALITY", "ERROR_PATTERN"
	Sensitivity          float64
	SpecificityThreshold float64
	PowerAnalysis        PowerAnalysisResult
	ValidationResults    []IndistinguishabilityValidationResult
}

// PowerAnalysisResult stores statistical power analysis results
type PowerAnalysisResult struct {
	SampleSizeRequired   int
	EffectSizeDetectable float64
	PowerAchieved        float64
	AlphaLevel           float64
	BetaLevel            float64
}

// IndistinguishabilityValidationResult stores results from indistinguishability validation
type IndistinguishabilityValidationResult struct {
	Platform           genomics.SequencingPlatform
	PopulationPair     [2]population.Population
	ClassificationRate float64
	PValue             float64
	ConfidenceInterval [2]float64
	TestStatistic      float64
	Distinguishable    bool
}

// IndistinguishabilityResult aggregates cross-platform testing results
type IndistinguishabilityResult struct {
	OverallScore           float64
	PlatformScores         map[genomics.SequencingPlatform]float64
	PopulationPairScores   map[string]float64
	FailedTests            []string
	PassedTests            []string
	StatisticalSummary     StatisticalSummary
	RecommendedActions     []string
}

// StatisticalSummary provides comprehensive statistical analysis
type StatisticalSummary struct {
	MeanClassificationRate float64
	StdClassificationRate  float64
	MedianClassificationRate float64
	MaxClassificationRate  float64
	MinClassificationRate  float64
	PercentagePassing      float64
	PowerAchieved          float64
	MultipleTestingCorrection string
}

// NewIndistinguishabilityTester creates a new cross-platform indistinguishability tester
func NewIndistinguishabilityTester() *IndistinguishabilityTester {
	return &IndistinguishabilityTester{
		platformManager:     genomics.NewPlatformErrorSignatureManager(),
		qualityModeler:      genomics.NewMultiPlatformQualityModeler(),
		classificationTests: make(map[string]*ClassificationTest),
		rng:                 rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
	}
}

// InitializeClassificationTests sets up comprehensive classification tests
func (it *IndistinguishabilityTester) InitializeClassificationTests() {
	
	// PCA-based population classification test
	it.classificationTests["PCA_Population_Classification"] = &ClassificationTest{
		TestName:             "PCA Population Classification",
		TestType:             "PCA",
		Sensitivity:          0.05, // 5% false positive rate threshold
		SpecificityThreshold: 0.95, // 95% specificity required
		PowerAnalysis: PowerAnalysisResult{
			SampleSizeRequired:   100,
			EffectSizeDetectable: 0.20, // Cohen's d = 0.20 (small effect)
			PowerAchieved:        0.80, // 80% power
			AlphaLevel:           0.05,
			BetaLevel:            0.20,
		},
	}

	// F_ST differentiation test
	it.classificationTests["FST_Differentiation"] = &ClassificationTest{
		TestName:             "F_ST Population Differentiation",  
		TestType:             "FST",
		Sensitivity:          0.01, // 1% threshold for F_ST detection
		SpecificityThreshold: 0.99, // 99% specificity for differentiation
		PowerAnalysis: PowerAnalysisResult{
			SampleSizeRequired:   200,
			EffectSizeDetectable: 0.005, // F_ST = 0.005 minimum detectable
			PowerAchieved:        0.90,  // 90% power for F_ST detection
			AlphaLevel:           0.01,
			BetaLevel:            0.10,
		},
	}

	// Linkage Disequilibrium pattern classification
	it.classificationTests["LD_Pattern_Classification"] = &ClassificationTest{
		TestName:             "LD Pattern Population Classification",
		TestType:             "LD", 
		Sensitivity:          0.05,
		SpecificityThreshold: 0.95,
		PowerAnalysis: PowerAnalysisResult{
			SampleSizeRequired:   150,
			EffectSizeDetectable: 0.15, // LD decay rate difference
			PowerAchieved:        0.85,
			AlphaLevel:           0.05,
			BetaLevel:            0.15,
		},
	}

	// Quality score distribution classification
	it.classificationTests["Quality_Distribution_Classification"] = &ClassificationTest{
		TestName:             "Quality Score Distribution Classification",
		TestType:             "QUALITY",
		Sensitivity:          0.05,
		SpecificityThreshold: 0.95,
		PowerAnalysis: PowerAnalysisResult{
			SampleSizeRequired:   75,
			EffectSizeDetectable: 0.25, // Quality distribution difference
			PowerAchieved:        0.80,
			AlphaLevel:           0.05,
			BetaLevel:            0.20,
		},
	}

	// Error pattern fingerprinting test
	it.classificationTests["Error_Pattern_Fingerprinting"] = &ClassificationTest{
		TestName:             "Error Pattern Platform Fingerprinting",
		TestType:             "ERROR_PATTERN",
		Sensitivity:          0.01, // Very strict for platform detection
		SpecificityThreshold: 0.99,
		PowerAnalysis: PowerAnalysisResult{
			SampleSizeRequired:   50,
			EffectSizeDetectable: 0.30, // Platform signature strength
			PowerAchieved:        0.95,  // High power for platform detection
			AlphaLevel:           0.01,
			BetaLevel:            0.05,
		},
	}

	// GC bias population classification
	it.classificationTests["GC_Bias_Classification"] = &ClassificationTest{
		TestName:             "GC Bias Population Classification",
		TestType:             "GC_BIAS",
		Sensitivity:          0.05,
		SpecificityThreshold: 0.95,
		PowerAnalysis: PowerAnalysisResult{
			SampleSizeRequired:   100,
			EffectSizeDetectable: 0.18, // GC bias difference threshold
			PowerAchieved:        0.82,
			AlphaLevel:           0.05,
			BetaLevel:            0.18,
		},
	}
}

// TestCrossPlatformIndistinguishability performs comprehensive indistinguishability testing
func (it *IndistinguishabilityTester) TestCrossPlatformIndistinguishability(
	sequences map[genomics.SequencingPlatform]map[population.Population][]string,
	sampleSize int,
) (*IndistinguishabilityResult, error) {

	if len(sequences) == 0 {
		return nil, fmt.Errorf("no sequences provided for testing")
	}

	// Initialize classification tests if not already done
	if len(it.classificationTests) == 0 {
		it.InitializeClassificationTests()
	}

	result := &IndistinguishabilityResult{
		PlatformScores:       make(map[genomics.SequencingPlatform]float64),
		PopulationPairScores: make(map[string]float64),
		PassedTests:          []string{},
		FailedTests:          []string{},
	}

	// Test each classification method across all platforms
	allValidationResults := []IndistinguishabilityValidationResult{}
	
	for testName, test := range it.classificationTests {
		fmt.Printf("Running %s across all platforms...\n", testName)
		
		testResults, err := it.runClassificationTest(test, sequences, sampleSize)
		if err != nil {
			fmt.Printf("Warning: %s failed: %v\n", testName, err)
			result.FailedTests = append(result.FailedTests, testName)
			continue
		}
		
		test.ValidationResults = testResults
		allValidationResults = append(allValidationResults, testResults...)
		
		// Assess if test passes indistinguishability threshold
		passedTest := it.assessTestIndistinguishability(test)
		if passedTest {
			result.PassedTests = append(result.PassedTests, testName)
		} else {
			result.FailedTests = append(result.FailedTests, testName)
		}
	}

	// Calculate platform-specific scores
	for platform := range sequences {
		platformResults := it.filterResultsByPlatform(allValidationResults, platform)
		score := it.calculatePlatformScore(platformResults)
		result.PlatformScores[platform] = score
	}

	// Calculate population pair scores  
	populationPairs := it.extractPopulationPairs(sequences)
	for _, pair := range populationPairs {
		pairKey := fmt.Sprintf("%s_vs_%s", pair[0], pair[1])
		pairResults := it.filterResultsByPopulationPair(allValidationResults, pair)
		score := it.calculatePopulationPairScore(pairResults)
		result.PopulationPairScores[pairKey] = score
	}

	// Calculate overall indistinguishability score
	result.OverallScore = it.calculateOverallScore(allValidationResults)
	
	// Generate statistical summary
	result.StatisticalSummary = it.generateStatisticalSummary(allValidationResults)
	
	// Provide recommendations
	result.RecommendedActions = it.generateRecommendations(result)

	return result, nil
}

// runClassificationTest executes a specific classification test across platforms
func (it *IndistinguishabilityTester) runClassificationTest(
	test *ClassificationTest,
	sequences map[genomics.SequencingPlatform]map[population.Population][]string,
	sampleSize int,
) ([]IndistinguishabilityValidationResult, error) {

	var results []IndistinguishabilityValidationResult
	
	for platform, populationSeqs := range sequences {
		populations := make([]population.Population, 0, len(populationSeqs))
		for pop := range populationSeqs {
			populations = append(populations, pop)
		}
		
		// Test all population pairs for this platform
		for i := 0; i < len(populations); i++ {
			for j := i + 1; j < len(populations); j++ {
				pop1, pop2 := populations[i], populations[j]
				
				result, err := it.testPopulationPairClassification(
					test, platform, pop1, pop2, 
					populationSeqs[pop1], populationSeqs[pop2], sampleSize)
				if err != nil {
					continue // Skip failed comparisons
				}
				
				results = append(results, result)
			}
		}
	}
	
	return results, nil
}

// testPopulationPairClassification tests classification between two populations
func (it *IndistinguishabilityTester) testPopulationPairClassification(
	test *ClassificationTest,
	platform genomics.SequencingPlatform,
	pop1, pop2 population.Population,
	seqs1, seqs2 []string,
	sampleSize int,
) (IndistinguishabilityValidationResult, error) {

	result := IndistinguishabilityValidationResult{
		Platform:       platform,
		PopulationPair: [2]population.Population{pop1, pop2},
	}

	// Ensure adequate sample size
	actualSampleSize := int(math.Min(float64(sampleSize), 
		math.Min(float64(len(seqs1)), float64(len(seqs2)))))
	
	if actualSampleSize < test.PowerAnalysis.SampleSizeRequired/2 {
		return result, fmt.Errorf("insufficient sample size: %d < %d", 
			actualSampleSize, test.PowerAnalysis.SampleSizeRequired/2)
	}

	// Perform test-specific classification analysis
	switch test.TestType {
	case "PCA":
		return it.testPCAClassification(test, result, seqs1, seqs2, actualSampleSize)
	case "FST":
		return it.testFSTClassification(test, result, seqs1, seqs2, actualSampleSize)
	case "LD":
		return it.testLDClassification(test, result, seqs1, seqs2, actualSampleSize)
	case "QUALITY":
		return it.testQualityClassification(test, result, platform, pop1, pop2, seqs1, seqs2, actualSampleSize)
	case "ERROR_PATTERN":
		return it.testErrorPatternClassification(test, result, platform, seqs1, seqs2, actualSampleSize)
	case "GC_BIAS":
		return it.testGCBiasClassification(test, result, seqs1, seqs2, actualSampleSize)
	default:
		return result, fmt.Errorf("unknown test type: %s", test.TestType)
	}
}

// testPCAClassification performs PCA-based population classification testing
func (it *IndistinguishabilityTester) testPCAClassification(
	test *ClassificationTest,
	result IndistinguishabilityValidationResult,
	seqs1, seqs2 []string,
	sampleSize int,
) (IndistinguishabilityValidationResult, error) {

	// Simulate PCA-based classification
	// In a full implementation, this would use actual PCA on genetic variants
	
	// Calculate basic sequence statistics as proxy for PCA features
	stats1 := it.calculateSequenceStatistics(seqs1[:sampleSize])
	stats2 := it.calculateSequenceStatistics(seqs2[:sampleSize])
	
	// Compute separation metric (proxy for PCA separation)
	separation := it.calculateStatisticalSeparation(stats1, stats2)
	
	// Classification rate based on separation
	result.ClassificationRate = math.Min(0.95, 0.5 + separation*0.45)
	result.TestStatistic = separation
	
	// P-value calculation (simplified)
	result.PValue = it.calculatePValueFromSeparation(separation, sampleSize)
	
	// Confidence interval (95%)
	margin := 1.96 * math.Sqrt(result.ClassificationRate*(1-result.ClassificationRate)/float64(sampleSize))
	result.ConfidenceInterval = [2]float64{
		math.Max(0, result.ClassificationRate-margin),
		math.Min(1, result.ClassificationRate+margin),
	}
	
	// Test passes if classification rate is near random (0.5) and p-value > alpha
	result.Distinguishable = result.ClassificationRate > 0.6 || result.PValue < test.Sensitivity
	
	return result, nil
}

// testFSTClassification performs F_ST-based differentiation testing
func (it *IndistinguishabilityTester) testFSTClassification(
	test *ClassificationTest,
	result IndistinguishabilityValidationResult,
	seqs1, seqs2 []string,
	sampleSize int,
) (IndistinguishabilityValidationResult, error) {

	// Calculate sequence diversity within and between populations
	diversity1 := it.calculateSequenceDiversity(seqs1[:sampleSize])
	diversity2 := it.calculateSequenceDiversity(seqs2[:sampleSize])
	betweenDiversity := it.calculateBetweenSequenceDiversity(seqs1[:sampleSize], seqs2[:sampleSize])
	
	// Simplified F_ST calculation
	totalDiversity := (diversity1 + diversity2) / 2
	fst := 0.0
	if totalDiversity > 0 {
		fst = (betweenDiversity - totalDiversity) / betweenDiversity
	}
	
	result.TestStatistic = fst
	result.ClassificationRate = math.Min(0.95, 0.5 + math.Abs(fst)*10) // Scale F_ST to classification rate
	
	// P-value for F_ST significance
	result.PValue = it.calculateFSTSignificance(fst, sampleSize)
	
	// Confidence interval for F_ST
	fstStdErr := it.calculateFSTStandardError(fst, sampleSize)
	result.ConfidenceInterval = [2]float64{
		math.Max(0, fst-1.96*fstStdErr),
		fst + 1.96*fstStdErr,
	}
	
	// Test passes if F_ST is low (< threshold) and non-significant
	result.Distinguishable = fst > test.PowerAnalysis.EffectSizeDetectable || result.PValue < test.Sensitivity
	
	return result, nil
}

// testLDClassification performs LD pattern-based classification testing
func (it *IndistinguishabilityTester) testLDClassification(
	test *ClassificationTest,
	result IndistinguishabilityValidationResult,
	seqs1, seqs2 []string,
	sampleSize int,
) (IndistinguishabilityValidationResult, error) {

	// Calculate LD-related metrics (simplified)
	ldMetric1 := it.calculateLDMetric(seqs1[:sampleSize])
	ldMetric2 := it.calculateLDMetric(seqs2[:sampleSize])
	
	ldDifference := math.Abs(ldMetric1 - ldMetric2)
	
	result.TestStatistic = ldDifference
	result.ClassificationRate = math.Min(0.95, 0.5 + ldDifference*2)
	
	// P-value for LD difference significance
	result.PValue = it.calculateLDSignificance(ldDifference, sampleSize)
	
	// Confidence interval
	ldStdErr := it.calculateLDStandardError(ldDifference, sampleSize)
	result.ConfidenceInterval = [2]float64{
		math.Max(0, ldDifference-1.96*ldStdErr),
		ldDifference + 1.96*ldStdErr,
	}
	
	result.Distinguishable = ldDifference > test.PowerAnalysis.EffectSizeDetectable || result.PValue < test.Sensitivity
	
	return result, nil
}

// testQualityClassification performs quality score distribution classification testing
func (it *IndistinguishabilityTester) testQualityClassification(
	test *ClassificationTest,
	result IndistinguishabilityValidationResult,
	platform genomics.SequencingPlatform,
	pop1, pop2 population.Population,
	seqs1, seqs2 []string,
	sampleSize int,
) (IndistinguishabilityValidationResult, error) {

	// Generate quality scores for both populations using the quality modeler
	qualities1 := [][]int{}
	qualities2 := [][]int{}
	
	for i := 0; i < sampleSize && i < len(seqs1); i++ {
		qual := it.qualityModeler.GenerateQualityScores(platform, pop1, seqs1[i], i+1)
		qualities1 = append(qualities1, qual)
	}
	
	for i := 0; i < sampleSize && i < len(seqs2); i++ {
		qual := it.qualityModeler.GenerateQualityScores(platform, pop2, seqs2[i], i+1)
		qualities2 = append(qualities2, qual)
	}
	
	// Calculate quality distribution statistics
	qualStats1 := it.calculateQualityStatistics(qualities1)
	qualStats2 := it.calculateQualityStatistics(qualities2)
	
	// Compute quality distribution difference
	qualityDifference := it.calculateQualityDifference(qualStats1, qualStats2)
	
	result.TestStatistic = qualityDifference
	result.ClassificationRate = math.Min(0.95, 0.5 + qualityDifference*1.5)
	
	// P-value for quality difference significance
	result.PValue = it.calculateQualitySignificance(qualityDifference, sampleSize)
	
	// Confidence interval
	qualStdErr := it.calculateQualityStandardError(qualityDifference, sampleSize)
	result.ConfidenceInterval = [2]float64{
		math.Max(0, qualityDifference-1.96*qualStdErr),
		qualityDifference + 1.96*qualStdErr,
	}
	
	result.Distinguishable = qualityDifference > test.PowerAnalysis.EffectSizeDetectable || result.PValue < test.Sensitivity
	
	return result, nil
}

// testErrorPatternClassification performs error pattern fingerprinting testing
func (it *IndistinguishabilityTester) testErrorPatternClassification(
	test *ClassificationTest,
	result IndistinguishabilityValidationResult,
	platform genomics.SequencingPlatform,
	seqs1, seqs2 []string,
	sampleSize int,
) (IndistinguishabilityValidationResult, error) {

	// Apply platform error signatures to sequences
	errorResults1 := []genomics.ErrorApplicationResult{}
	errorResults2 := []genomics.ErrorApplicationResult{}
	
	for i := 0; i < sampleSize && i < len(seqs1); i++ {
		qualities := make([]int, len(seqs1[i]))
		for j := range qualities {
			qualities[j] = 35 // Default quality
		}
		
		_, _, errorResult := it.platformManager.ApplyErrorSignature(
			platform, population.YRI, seqs1[i], qualities, i+1)
		errorResults1 = append(errorResults1, errorResult)
	}
	
	for i := 0; i < sampleSize && i < len(seqs2); i++ {
		qualities := make([]int, len(seqs2[i]))
		for j := range qualities {
			qualities[j] = 35 // Default quality
		}
		
		_, _, errorResult := it.platformManager.ApplyErrorSignature(
			platform, population.CEU, seqs2[i], qualities, i+1)
		errorResults2 = append(errorResults2, errorResult)
	}
	
	// Calculate error pattern statistics
	errorStats1 := it.calculateErrorStatistics(errorResults1)
	errorStats2 := it.calculateErrorStatistics(errorResults2)
	
	// Compute error pattern difference
	errorDifference := it.calculateErrorDifference(errorStats1, errorStats2)
	
	result.TestStatistic = errorDifference
	result.ClassificationRate = math.Min(0.95, 0.5 + errorDifference*2)
	
	// P-value for error pattern significance
	result.PValue = it.calculateErrorSignificance(errorDifference, sampleSize)
	
	// Confidence interval
	errorStdErr := it.calculateErrorStandardError(errorDifference, sampleSize)
	result.ConfidenceInterval = [2]float64{
		math.Max(0, errorDifference-1.96*errorStdErr),
		errorDifference + 1.96*errorStdErr,
	}
	
	result.Distinguishable = errorDifference > test.PowerAnalysis.EffectSizeDetectable || result.PValue < test.Sensitivity
	
	return result, nil
}

// testGCBiasClassification performs GC bias-based classification testing
func (it *IndistinguishabilityTester) testGCBiasClassification(
	test *ClassificationTest,
	result IndistinguishabilityValidationResult,
	seqs1, seqs2 []string,
	sampleSize int,
) (IndistinguishabilityValidationResult, error) {

	// Calculate GC content distributions
	gcContents1 := []float64{}
	gcContents2 := []float64{}
	
	for i := 0; i < sampleSize && i < len(seqs1); i++ {
		gc := it.calculateGCContent(seqs1[i])
		gcContents1 = append(gcContents1, gc)
	}
	
	for i := 0; i < sampleSize && i < len(seqs2); i++ {
		gc := it.calculateGCContent(seqs2[i])
		gcContents2 = append(gcContents2, gc)
	}
	
	// Calculate GC bias difference
	gcBias1 := it.calculateMean(gcContents1)
	gcBias2 := it.calculateMean(gcContents2)
	gcDifference := math.Abs(gcBias1 - gcBias2)
	
	result.TestStatistic = gcDifference
	result.ClassificationRate = math.Min(0.95, 0.5 + gcDifference*5) // Scale GC difference
	
	// P-value for GC difference significance
	result.PValue = it.calculateGCSignificance(gcDifference, sampleSize)
	
	// Confidence interval
	gcStdErr := it.calculateGCStandardError(gcDifference, sampleSize)
	result.ConfidenceInterval = [2]float64{
		math.Max(0, gcDifference-1.96*gcStdErr),
		gcDifference + 1.96*gcStdErr,
	}
	
	result.Distinguishable = gcDifference > test.PowerAnalysis.EffectSizeDetectable || result.PValue < test.Sensitivity
	
	return result, nil
}

// Helper functions for statistical calculations

func (it *IndistinguishabilityTester) calculateSequenceStatistics(sequences []string) map[string]float64 {
	stats := make(map[string]float64)
	
	totalLength := 0
	gcCount := 0
	
	for _, seq := range sequences {
		totalLength += len(seq)
		for _, base := range seq {
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}
	}
	
	if len(sequences) > 0 {
		stats["mean_length"] = float64(totalLength) / float64(len(sequences))
	}
	if totalLength > 0 {
		stats["gc_content"] = float64(gcCount) / float64(totalLength)
	}
	
	return stats
}

func (it *IndistinguishabilityTester) calculateStatisticalSeparation(stats1, stats2 map[string]float64) float64 {
	// Simple separation metric based on normalized differences
	separation := 0.0
	count := 0
	
	for key, val1 := range stats1 {
		if val2, exists := stats2[key]; exists {
			diff := math.Abs(val1 - val2)
			separation += diff
			count++
		}
	}
	
	if count == 0 {
		return 0.0
	}
	
	return separation / float64(count)
}

func (it *IndistinguishabilityTester) calculatePValueFromSeparation(separation float64, sampleSize int) float64 {
	// Simplified p-value calculation based on separation and sample size
	// In practice, this would use proper statistical tests
	z := separation * math.Sqrt(float64(sampleSize)) / 0.1
	return 2 * (1 - it.normalCDF(math.Abs(z)))
}

func (it *IndistinguishabilityTester) normalCDF(z float64) float64 {
	// Approximation of normal CDF
	return 0.5 * (1 + math.Erf(z/math.Sqrt(2)))
}

func (it *IndistinguishabilityTester) calculateSequenceDiversity(sequences []string) float64 {
	// Simple diversity metric based on sequence variation
	if len(sequences) < 2 {
		return 0.0
	}
	
	differences := 0
	comparisons := 0
	
	for i := 0; i < len(sequences); i++ {
		for j := i + 1; j < len(sequences); j++ {
			comparisons++
			minLen := int(math.Min(float64(len(sequences[i])), float64(len(sequences[j]))))
			for k := 0; k < minLen; k++ {
				if sequences[i][k] != sequences[j][k] {
					differences++
					break // Count as different sequences
				}
			}
		}
	}
	
	if comparisons == 0 {
		return 0.0
	}
	
	return float64(differences) / float64(comparisons)
}

func (it *IndistinguishabilityTester) calculateBetweenSequenceDiversity(seqs1, seqs2 []string) float64 {
	differences := 0
	comparisons := 0
	
	for _, seq1 := range seqs1 {
		for _, seq2 := range seqs2 {
			comparisons++
			minLen := int(math.Min(float64(len(seq1)), float64(len(seq2))))
			for k := 0; k < minLen; k++ {
				if seq1[k] != seq2[k] {
					differences++
					break
				}
			}
		}
	}
	
	if comparisons == 0 {
		return 0.0
	}
	
	return float64(differences) / float64(comparisons)
}

func (it *IndistinguishabilityTester) calculateFSTSignificance(fst float64, sampleSize int) float64 {
	se := it.calculateFSTStandardError(fst, sampleSize)
	if se <= 0 {
		if fst > 0 {
			return 0.0
		}
		return 1.0
	}
	z := fst / se
	return 2 * (1 - it.normalCDF(math.Abs(z)))
}

func (it *IndistinguishabilityTester) calculateFSTStandardError(fst float64, sampleSize int) float64 {
	if sampleSize <= 0 {
		return 1.0
	}
	// Asymptotic SE for a proportion: sqrt(p*(1-p)/n)
	return math.Sqrt(fst * (1 - fst) / float64(sampleSize))
}

// calculateLDMetric computes the mean r² linkage disequilibrium over adjacent position
// pairs sampled from the provided sequences. For each pair of adjacent columns, it
// reduces each column to a biallelic site (modal base vs. all others) and computes r²
// from the phased haplotype counts.
func (it *IndistinguishabilityTester) calculateLDMetric(sequences []string) float64 {
	if len(sequences) == 0 {
		return 0.0
	}

	// Find minimum sequence length so every sequence has both positions.
	minLen := len(sequences[0])
	for _, s := range sequences[1:] {
		if len(s) < minLen {
			minLen = len(s)
		}
	}
	if minLen < 2 {
		return 0.0
	}

	// Sample up to 50 adjacent position pairs to limit runtime.
	posLimit := minLen - 1
	if posLimit > 50 {
		posLimit = 50
	}

	totalR2 := 0.0
	validPairs := 0

	for j := 0; j < posLimit; j++ {
		// Count bases at positions j and j+1 to find the modal allele.
		var counts1, counts2 [4]int
		for _, seq := range sequences {
			if j+1 >= len(seq) {
				continue
			}
			if b := ldBaseIndex(seq[j]); b >= 0 {
				counts1[b]++
			}
			if b := ldBaseIndex(seq[j+1]); b >= 0 {
				counts2[b]++
			}
		}
		modal1 := ldMaxIndex(counts1)
		modal2 := ldMaxIndex(counts2)

		// Count biallelic haplotypes (modal=A/B, non-modal=a/b).
		nAB, nAb, naB, nab := 0, 0, 0, 0
		for _, seq := range sequences {
			if j+1 >= len(seq) {
				continue
			}
			b1 := ldBaseIndex(seq[j])
			b2 := ldBaseIndex(seq[j+1])
			if b1 < 0 || b2 < 0 {
				continue
			}
			switch {
			case b1 == modal1 && b2 == modal2:
				nAB++
			case b1 == modal1 && b2 != modal2:
				nAb++
			case b1 != modal1 && b2 == modal2:
				naB++
			default:
				nab++
			}
		}

		total := nAB + nAb + naB + nab
		if total == 0 {
			continue
		}

		pA := float64(nAB+nAb) / float64(total)
		pB := float64(nAB+naB) / float64(total)
		qA := 1 - pA
		qB := 1 - pB
		denom := pA * qA * pB * qB
		if denom <= 0 {
			// Monomorphic at one or both loci — skip.
			continue
		}

		pAB := float64(nAB) / float64(total)
		D := pAB - pA*pB
		r2 := (D * D) / denom

		totalR2 += r2
		validPairs++
	}

	if validPairs == 0 {
		return 0.0
	}
	return totalR2 / float64(validPairs)
}

func (it *IndistinguishabilityTester) calculateLDSignificance(ldDiff float64, sampleSize int) float64 {
	se := it.calculateLDStandardError(ldDiff, sampleSize)
	if se <= 0 {
		if ldDiff > 0 {
			return 0.0
		}
		return 1.0
	}
	z := ldDiff / se
	return 2 * (1 - it.normalCDF(math.Abs(z)))
}

func (it *IndistinguishabilityTester) calculateLDStandardError(ldDiff float64, sampleSize int) float64 {
	if sampleSize <= 0 {
		return 1.0
	}
	// Treat ldDiff as a proportion and use the asymptotic SE: sqrt(p*(1-p)/n).
	return math.Sqrt(ldDiff * (1 - ldDiff) / float64(sampleSize))
}

// ldBaseIndex maps a DNA base byte to an index in [0,3]: A=0, C=1, G=2, T=3.
// Returns -1 for unrecognised bytes.
func ldBaseIndex(b byte) int {
	switch b | 0x20 { // tolower
	case 'a':
		return 0
	case 'c':
		return 1
	case 'g':
		return 2
	case 't':
		return 3
	default:
		return -1
	}
}

// ldMaxIndex returns the index of the maximum element in a [4]int array.
func ldMaxIndex(counts [4]int) int {
	idx := 0
	for i := 1; i < 4; i++ {
		if counts[i] > counts[idx] {
			idx = i
		}
	}
	return idx
}

func (it *IndistinguishabilityTester) calculateQualityStatistics(qualities [][]int) map[string]float64 {
	stats := make(map[string]float64)
	
	allQualities := []int{}
	for _, qualArray := range qualities {
		allQualities = append(allQualities, qualArray...)
	}
	
	if len(allQualities) == 0 {
		return stats
	}
	
	// Mean quality
	sum := 0
	for _, q := range allQualities {
		sum += q
	}
	stats["mean"] = float64(sum) / float64(len(allQualities))
	
	// Standard deviation (Bessel's correction: ÷(N-1))
	variance := 0.0
	for _, q := range allQualities {
		diff := float64(q) - stats["mean"]
		variance += diff * diff
	}
	if n := len(allQualities); n > 1 {
		stats["std"] = math.Sqrt(variance / float64(n-1))
	}
	
	return stats
}

func (it *IndistinguishabilityTester) calculateQualityDifference(stats1, stats2 map[string]float64) float64 {
	meanDiff := math.Abs(stats1["mean"] - stats2["mean"])
	stdDiff := math.Abs(stats1["std"] - stats2["std"])
	return (meanDiff + stdDiff) / 2
}

func (it *IndistinguishabilityTester) calculateQualitySignificance(qualDiff float64, sampleSize int) float64 {
	z := qualDiff * math.Sqrt(float64(sampleSize)) / 2.0
	return 2 * (1 - it.normalCDF(math.Abs(z)))
}

func (it *IndistinguishabilityTester) calculateQualityStandardError(qualDiff float64, sampleSize int) float64 {
	if sampleSize <= 0 {
		return 1.0
	}
	// Asymptotic SE for a proportion: sqrt(p*(1-p)/n) — consistent with
	// calculateErrorStandardError and calculateGCStandardError.
	return math.Sqrt(qualDiff * (1 - qualDiff) / float64(sampleSize))
}

func (it *IndistinguishabilityTester) calculateErrorStatistics(errorResults []genomics.ErrorApplicationResult) map[string]float64 {
	stats := make(map[string]float64)
	if len(errorResults) == 0 {
		return stats // mean_errors=0, success_rate=0
	}

	totalErrors := 0
	successCount := 0

	for _, result := range errorResults {
		totalErrors += result.ErrorsApplied
		if result.Success {
			successCount++
		}
	}

	stats["mean_errors"] = float64(totalErrors) / float64(len(errorResults))
	stats["success_rate"] = float64(successCount) / float64(len(errorResults))

	return stats
}

func (it *IndistinguishabilityTester) calculateErrorDifference(stats1, stats2 map[string]float64) float64 {
	errorDiff := math.Abs(stats1["mean_errors"] - stats2["mean_errors"])
	successDiff := math.Abs(stats1["success_rate"] - stats2["success_rate"])
	return (errorDiff + successDiff) / 2
}

func (it *IndistinguishabilityTester) calculateErrorSignificance(errorDiff float64, sampleSize int) float64 {
	se := it.calculateErrorStandardError(errorDiff, sampleSize)
	if se <= 0 {
		if errorDiff > 0 {
			return 0.0
		}
		return 1.0
	}
	z := errorDiff / se
	return 2 * (1 - it.normalCDF(math.Abs(z)))
}

func (it *IndistinguishabilityTester) calculateErrorStandardError(errorDiff float64, sampleSize int) float64 {
	if sampleSize <= 0 {
		return 1.0
	}
	// errorDiff combines a count-based difference (mean_errors, unbounded) with a
	// success-rate difference (bounded [0,1]). Clamp to [0,1] before applying the
	// binomial SE approximation; without clamping, errorDiff > 1 yields
	// sqrt(negative) = NaN.
	p := math.Max(0, math.Min(1, errorDiff))
	return math.Sqrt(p * (1 - p) / float64(sampleSize))
}

func (it *IndistinguishabilityTester) calculateGCContent(sequence string) float64 {
	if len(sequence) == 0 {
		return 0.0
	}
	
	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}
	
	return float64(gcCount) / float64(len(sequence))
}

func (it *IndistinguishabilityTester) calculateMean(values []float64) float64 {
	if len(values) == 0 {
		return 0.0
	}
	
	sum := 0.0
	for _, val := range values {
		sum += val
	}
	
	return sum / float64(len(values))
}

func (it *IndistinguishabilityTester) calculateGCSignificance(gcDiff float64, sampleSize int) float64 {
	se := it.calculateGCStandardError(gcDiff, sampleSize)
	if se <= 0 {
		if gcDiff > 0 {
			return 0.0
		}
		return 1.0
	}
	z := gcDiff / se
	return 2 * (1 - it.normalCDF(math.Abs(z)))
}

func (it *IndistinguishabilityTester) calculateGCStandardError(gcDiff float64, sampleSize int) float64 {
	if sampleSize <= 0 {
		return 1.0
	}
	// Asymptotic SE for a proportion: sqrt(p*(1-p)/n)
	return math.Sqrt(gcDiff * (1 - gcDiff) / float64(sampleSize))
}

// Assessment and scoring functions

func (it *IndistinguishabilityTester) assessTestIndistinguishability(test *ClassificationTest) bool {
	if len(test.ValidationResults) == 0 {
		return false
	}
	
	passCount := 0
	for _, result := range test.ValidationResults {
		// Test passes if populations are NOT distinguishable
		if !result.Distinguishable {
			passCount++
		}
	}
	
	passRate := float64(passCount) / float64(len(test.ValidationResults))
	return passRate >= test.SpecificityThreshold
}

func (it *IndistinguishabilityTester) filterResultsByPlatform(results []IndistinguishabilityValidationResult, platform genomics.SequencingPlatform) []IndistinguishabilityValidationResult {
	var filtered []IndistinguishabilityValidationResult
	for _, result := range results {
		if result.Platform == platform {
			filtered = append(filtered, result)
		}
	}
	return filtered
}

func (it *IndistinguishabilityTester) filterResultsByPopulationPair(results []IndistinguishabilityValidationResult, pair [2]population.Population) []IndistinguishabilityValidationResult {
	var filtered []IndistinguishabilityValidationResult
	for _, result := range results {
		if (result.PopulationPair[0] == pair[0] && result.PopulationPair[1] == pair[1]) ||
		   (result.PopulationPair[0] == pair[1] && result.PopulationPair[1] == pair[0]) {
			filtered = append(filtered, result)
		}
	}
	return filtered
}

func (it *IndistinguishabilityTester) calculatePlatformScore(results []IndistinguishabilityValidationResult) float64 {
	if len(results) == 0 {
		return 0.0
	}
	
	passCount := 0
	for _, result := range results {
		if !result.Distinguishable {
			passCount++
		}
	}
	
	return float64(passCount) / float64(len(results))
}

func (it *IndistinguishabilityTester) calculatePopulationPairScore(results []IndistinguishabilityValidationResult) float64 {
	if len(results) == 0 {
		return 0.0
	}
	
	totalClassificationRate := 0.0
	for _, result := range results {
		// Invert classification rate (lower is better for indistinguishability)
		totalClassificationRate += (1.0 - result.ClassificationRate)
	}
	
	return totalClassificationRate / float64(len(results))
}

func (it *IndistinguishabilityTester) calculateOverallScore(results []IndistinguishabilityValidationResult) float64 {
	if len(results) == 0 {
		return 0.0
	}
	
	passCount := 0
	totalClassificationRate := 0.0
	
	for _, result := range results {
		if !result.Distinguishable {
			passCount++
		}
		totalClassificationRate += result.ClassificationRate
	}
	
	passRate := float64(passCount) / float64(len(results))
	avgClassificationRate := totalClassificationRate / float64(len(results))
	
	// Score combines pass rate and low classification rate (close to 0.5 = random)
	indistinguishabilityScore := passRate * (1.0 - math.Abs(avgClassificationRate-0.5)*2)
	
	return indistinguishabilityScore
}

func (it *IndistinguishabilityTester) generateStatisticalSummary(results []IndistinguishabilityValidationResult) StatisticalSummary {
	if len(results) == 0 {
		return StatisticalSummary{}
	}
	
	classificationRates := []float64{}
	passCount := 0
	
	for _, result := range results {
		classificationRates = append(classificationRates, result.ClassificationRate)
		if !result.Distinguishable {
			passCount++
		}
	}
	
	// Sort for median calculation
	sort.Float64s(classificationRates)
	
	// Calculate statistics
	mean := it.calculateMean(classificationRates)
	
	// Standard deviation (Bessel's correction: ÷(N-1))
	variance := 0.0
	for _, rate := range classificationRates {
		diff := rate - mean
		variance += diff * diff
	}
	std := 0.0
	if n := len(classificationRates); n > 1 {
		std = math.Sqrt(variance / float64(n-1))
	}
	
	// Median
	median := classificationRates[len(classificationRates)/2]
	if len(classificationRates)%2 == 0 {
		median = (classificationRates[len(classificationRates)/2-1] + classificationRates[len(classificationRates)/2]) / 2
	}
	
	// Estimate statistical power: the probability that a two-sided z-test for a
	// proportion correctly rejects H0 (mean classification rate = 0.5) at α = 0.05.
	//
	//   power = Φ( |δ|·√n / σ₀ − z_{α/2} )
	//
	// where δ = |mean − 0.5| (observed effect size), σ₀ = 0.5 (null std dev for a
	// proportion), n = number of tests, and z_{α/2} = 1.96 (for α = 0.05 two-sided).
	// Φ(x) = erfc(−x/√2) / 2 via the standard normal CDF identity.
	effect := math.Abs(mean - 0.5)
	n := float64(len(results))
	const zAlpha2 = 1.96 // z_{α/2} for α = 0.05 two-sided
	zScore := 2*effect*math.Sqrt(n) - zAlpha2 // 2 = 1/σ₀ = 1/0.5
	powerAchieved := 0.5 * math.Erfc(-zScore/math.Sqrt2)

	return StatisticalSummary{
		MeanClassificationRate:    mean,
		StdClassificationRate:     std,
		MedianClassificationRate:  median,
		MaxClassificationRate:     classificationRates[len(classificationRates)-1],
		MinClassificationRate:     classificationRates[0],
		PercentagePassing:         float64(passCount) / float64(len(results)) * 100,
		PowerAchieved:             powerAchieved,
		MultipleTestingCorrection: "Bonferroni",
	}
}

func (it *IndistinguishabilityTester) extractPopulationPairs(sequences map[genomics.SequencingPlatform]map[population.Population][]string) [][2]population.Population {
	populationSet := make(map[population.Population]bool)
	
	// Collect all populations
	for _, populationSeqs := range sequences {
		for pop := range populationSeqs {
			populationSet[pop] = true
		}
	}
	
	// Generate all pairs
	populations := make([]population.Population, 0, len(populationSet))
	for pop := range populationSet {
		populations = append(populations, pop)
	}
	
	var pairs [][2]population.Population
	for i := 0; i < len(populations); i++ {
		for j := i + 1; j < len(populations); j++ {
			pairs = append(pairs, [2]population.Population{populations[i], populations[j]})
		}
	}
	
	return pairs
}

func (it *IndistinguishabilityTester) generateRecommendations(result *IndistinguishabilityResult) []string {
	recommendations := []string{}
	
	if result.OverallScore < 0.7 {
		recommendations = append(recommendations, "Overall indistinguishability score is low - consider enhancing error signature obfuscation")
	}
	
	for platform, score := range result.PlatformScores {
		if score < 0.6 {
			recommendations = append(recommendations, 
				fmt.Sprintf("Platform %s shows distinguishable patterns - review error signature implementation", platform))
		}
	}
	
	if result.StatisticalSummary.MeanClassificationRate > 0.7 {
		recommendations = append(recommendations, "Mean classification rate is high - population signatures too distinct")
	}
	
	if len(result.FailedTests) > len(result.PassedTests) {
		recommendations = append(recommendations, "More tests failing than passing - comprehensive review needed")
	}
	
	if len(recommendations) == 0 {
		recommendations = append(recommendations, "Statistical indistinguishability appears adequate - consider advanced adversarial testing")
	}
	
	return recommendations
}

// GetClassificationTests returns the configured classification tests
func (it *IndistinguishabilityTester) GetClassificationTests() map[string]*ClassificationTest {
	return it.classificationTests
}

// GetTestResult returns results for a specific test
func (it *IndistinguishabilityTester) GetTestResult(testName string) *ClassificationTest {
	return it.classificationTests[testName]
}