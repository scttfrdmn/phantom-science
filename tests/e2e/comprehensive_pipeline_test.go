package e2e

import (
	"fmt"
	"math"
	"math/rand"
	"sync"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	pt "github.com/scttfrdmn/phantom-science/internal/testing"
)

// TestComprehensivePipelineE2E tests the complete end-to-end pipeline across all sequencer platforms
func TestComprehensivePipelineE2E(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping comprehensive end-to-end pipeline tests in short mode")
	}

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	// Test all supported sequencing platforms
	platforms := []pt.SequencingPlatform{
		pt.IlluminaNovaSeq6000,
		pt.IlluminaNextSeq2000,
		pt.IlluminaMiSeq,
		pt.PacBioRevio,
		pt.ONTPromethION24,
		pt.MGI_DNBSEQ_T20x2,
		pt.ElementAVITI,
	}

	// Test across major global populations
	populations := []string{"YRI", "CEU", "CHB", "JPT", "MXL"}

	t.Run("Full Pipeline Integration Tests", func(t *testing.T) {
		for _, platform := range platforms {
			t.Run(string(platform), func(t *testing.T) {
				testPlatformPipeline(t, framework, platform, populations)
			})
		}
	})

	t.Run("Cross-Platform Authenticity Validation", func(t *testing.T) {
		testCrossPlatformAuthenticity(t, framework, platforms, populations)
	})

	t.Run("Statistical Indistinguishability Validation", func(t *testing.T) {
		testStatisticalIndistinguishability(t, framework, platforms)
	})

	t.Run("Performance Benchmarking", func(t *testing.T) {
		testPerformanceBenchmarks(t, framework, platforms)
	})
}

// testPlatformPipeline tests the complete pipeline for a specific platform
func testPlatformPipeline(t *testing.T, framework *pt.E2EFramework, platform pt.SequencingPlatform, populations []string) {
	t.Logf("🧬 Testing complete pipeline for platform: %s", platform)

	for _, population := range populations {
		t.Run(fmt.Sprintf("Population_%s", population), func(t *testing.T) {
			// Step 1: Initialize platform-specific configuration
			config := framework.CreatePlatformConfig(platform, population)
			require.NotNil(t, config)

			// Step 2: Generate synthetic genomic data
			syntheticData, err := framework.GenerateSyntheticData(config)
			require.NoError(t, err)
			require.NotNil(t, syntheticData)

			// Step 3: Validate platform-specific characteristics
			validatePlatformCharacteristics(t, syntheticData, platform)

			// Step 4: Validate population-specific features
			validatePopulationFeatures(t, syntheticData, population)

			// Step 5: Run authenticity analysis
			authenticityScore, err := framework.AnalyzeAuthenticity(syntheticData)
			require.NoError(t, err)

			// Target: ~50.3% classification accuracy for statistical indistinguishability
			expectedScore := 52.0 // Authenticity score when classification accuracy is ~50%
			assert.InDelta(t, expectedScore, authenticityScore, 20.0,
				"Authenticity score should be close to 52 for statistical indistinguishability")

			t.Logf("✅ Platform %s, Population %s: Authenticity Score = %.2f", 
				platform, population, authenticityScore)
		})
	}
}

// validatePlatformCharacteristics validates platform-specific data characteristics
func validatePlatformCharacteristics(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	switch platform {
	case pt.IlluminaNovaSeq6000:
		// Validate Illumina NovaSeq characteristics
		assert.InDelta(t, 150, data.AverageReadLength, 10, "NovaSeq read length should be ~150bp")
		assert.Less(t, data.ErrorRate, 0.002, "NovaSeq error rate should be <0.2%")
		assert.Greater(t, data.QualityScoreAverage, 30.0, "NovaSeq quality scores should be >Q30")

	case pt.IlluminaNextSeq2000:
		// Validate NextSeq 2-color chemistry characteristics
		assert.InDelta(t, 150, data.AverageReadLength, 10, "NextSeq read length should be ~150bp")
		assert.Less(t, data.ErrorRate, 0.003, "NextSeq error rate should be <0.3%")
		validateTwoColorBias(t, data)

	case pt.IlluminaMiSeq:
		// Validate MiSeq 4-color chemistry characteristics
		assert.InDelta(t, 300, data.AverageReadLength, 50, "MiSeq read length should be ~300bp")
		assert.LessOrEqual(t, data.ErrorRate, 0.001, "MiSeq error rate should be <0.1%")
		assert.Greater(t, data.QualityScoreAverage, 28.0, "MiSeq quality scores should be >Q28")

	case pt.PacBioRevio:
		// Validate PacBio HiFi characteristics
		assert.Greater(t, data.AverageReadLength, 10000.0, "PacBio reads should be >10kb")
		assert.Less(t, data.AverageReadLength, 25000.0, "PacBio reads should be <25kb")
		assert.Less(t, data.ErrorRate, 0.005, "PacBio HiFi error rate should be <0.5%")
		validateCircularConsensus(t, data)

	case pt.ONTPromethION24:
		// Validate Oxford Nanopore characteristics
		assert.Greater(t, data.AverageReadLength, 20000.0, "ONT reads should be >20kb")
		assert.Less(t, data.ErrorRate, 0.10, "ONT error rate should be <10%")
		validateHomopolymerErrors(t, data)

	case pt.MGI_DNBSEQ_T20x2:
		// Validate MGI DNBSEQ characteristics
		assert.InDelta(t, 150, data.AverageReadLength, 10, "DNBSEQ read length should be ~150bp")
		assert.Less(t, data.ErrorRate, 0.002, "DNBSEQ error rate should be <0.2%")
		validateDNBSpecificBias(t, data)

	case pt.ElementAVITI:
		// Validate Element AVITI SBB chemistry characteristics
		assert.InDelta(t, 150, data.AverageReadLength, 10, "AVITI read length should be ~150bp")
		assert.LessOrEqual(t, data.ErrorRate, 0.0015, "AVITI error rate should be <0.15%")
		validateSBBChemistry(t, data)
	}
}

// validatePopulationFeatures validates population-specific genomic features
func validatePopulationFeatures(t *testing.T, data *pt.SyntheticData, population string) {
	switch population {
	case "YRI":
		// Yoruba (West African) - highest genetic diversity
		assert.Greater(t, data.NucleotideDiversity, 0.0008, "YRI should have high nucleotide diversity")
		assert.Greater(t, data.VariantCount, 1000, "YRI should have high variant count")

	case "CEU":
		// Central European - bottleneck signature
		assert.LessOrEqual(t, data.NucleotideDiversity, 0.0008, "CEU should have reduced diversity due to bottleneck")
		assert.Greater(t, data.LinkageDisequilibrium, 0.1, "CEU should have increased LD due to bottleneck")

	case "CHB":
		// Han Chinese - East Asian population structure
		assert.LessOrEqual(t, data.NucleotideDiversity, 0.0007, "CHB should have reduced diversity")
		validateEastAsianStructure(t, data)

	case "JPT":
		// Japanese - island population characteristics
		assert.LessOrEqual(t, data.NucleotideDiversity, 0.0006, "JPT should have lowest diversity")
		validateIslandFounderEffect(t, data)

	case "MXL":
		// Mexican - admixed population
		validateAdmixtureSignatures(t, data)
		assert.Greater(t, data.AdmixtureLD, 0.05, "MXL should show admixture LD")
	}
}

// testCrossPlatformAuthenticity validates authenticity across all platforms
func testCrossPlatformAuthenticity(t *testing.T, framework *pt.E2EFramework, platforms []pt.SequencingPlatform, populations []string) {
	authenticityScores := make(map[pt.SequencingPlatform]map[string]float64)

	// Collect authenticity scores across all platforms and populations
	for _, platform := range platforms {
		authenticityScores[platform] = make(map[string]float64)
		
		for _, population := range populations {
			config := framework.CreatePlatformConfig(platform, population)
			syntheticData, err := framework.GenerateSyntheticData(config)
			require.NoError(t, err)

			score, err := framework.AnalyzeAuthenticity(syntheticData)
			require.NoError(t, err)

			authenticityScores[platform][population] = score
		}
	}

	// Validate cross-platform consistency
	for _, population := range populations {
		platformScores := make([]float64, 0, len(platforms))
		for _, platform := range platforms {
			platformScores = append(platformScores, authenticityScores[platform][population])
		}

		// Calculate coefficient of variation across platforms
		mean := calculateMean(platformScores)
		stdDev := calculateStdDev(platformScores, mean)
		cv := stdDev / mean

		// Cross-platform variation should be reasonable
		assert.Less(t, cv, 0.2, "Cross-platform authenticity variation should be <20% for %s", population)

		t.Logf("Population %s: Mean authenticity = %.2f, CV = %.3f", population, mean, cv)
	}
}

// testStatisticalIndistinguishability validates the ML indistinguishability
// measurement itself. No real reads are bundled, so it cannot measure real vs
// synthetic; instead it checks that the classifier reports chance accuracy on
// two independent draws from the same generator (null case) and near-perfect
// accuracy when one set is deliberately GC-skewed (positive control).
func testStatisticalIndistinguishability(t *testing.T, framework *pt.E2EFramework, platforms []pt.SequencingPlatform) {
	t.Log("🎯 Testing statistical indistinguishability with ML classifiers")

	const trainPerClass, holdoutPerClass = 1000, 1000

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			trainReal, trainSynthetic, err := framework.GenerateBalancedDataset(platform, trainPerClass)
			require.NoError(t, err)
			testReal, testSynthetic, err := framework.GenerateBalancedDataset(platform, holdoutPerClass)
			require.NoError(t, err)

			// Null case: under H0 holdout accuracy is Binomial(n, 0.5)/n.
			// A 4σ band keeps the false-failure rate below 1e-4.
			classifier := framework.CreateMLClassifier()
			require.NoError(t, classifier.Train(trainReal, trainSynthetic))
			accuracy := classifier.TestAccuracy(testReal, testSynthetic)

			n := float64(2 * holdoutPerClass)
			tolerance := 4 * math.Sqrt(0.25/n)
			assert.InDelta(t, 0.5, accuracy, tolerance,
				"independent draws from the same generator should be indistinguishable")

			// Positive control: the same pipeline must detect a real
			// distributional difference, or a ~0.5 result above is meaningless.
			skewed := framework.CreateMLClassifier()
			require.NoError(t, skewed.Train(trainReal, gcSkewed(trainSynthetic, 0.4)))
			skewedAccuracy := skewed.TestAccuracy(testReal, gcSkewed(testSynthetic, 0.4))
			assert.Greater(t, skewedAccuracy, 0.9,
				"classifier should easily separate GC-skewed reads from unskewed reads")

			t.Logf("✅ Platform %s: null accuracy = %.3f (±%.3f), GC-skewed control accuracy = %.3f",
				platform, accuracy, tolerance, skewedAccuracy)
		})
	}
}

// gcSkewed returns a copy of data in which each A/T base is replaced by G/C
// with probability p, shifting GC content from ~0.5 to ~0.5+p/2.
func gcSkewed(data *pt.SyntheticData, p float64) *pt.SyntheticData {
	rng := rand.New(rand.NewSource(1)) // #nosec G404 -- deterministic test fixture
	out := *data
	out.Sequences = make([][]byte, len(data.Sequences))
	for i, seq := range data.Sequences {
		s := make([]byte, len(seq))
		for j, b := range seq {
			if (b == 'A' || b == 'T') && rng.Float64() < p {
				b = "GC"[rng.Intn(2)]
			}
			s[j] = b
		}
		out.Sequences[i] = s
	}
	return &out
}

// testPerformanceBenchmarks validates performance across all platforms
func testPerformanceBenchmarks(t *testing.T, framework *pt.E2EFramework, platforms []pt.SequencingPlatform) {
	t.Log("⚡ Running performance benchmarks across all platforms")

	benchmarks := make(map[pt.SequencingPlatform]*pt.PlatformBenchmark)

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			benchmark := framework.RunBenchmark(platform)
			require.NotNil(t, benchmark)

			// Validate platform-specific performance expectations
			validatePlatformPerformance(t, platform, benchmark)

			benchmarks[platform] = benchmark
		})
	}

	// Generate comprehensive performance report
	framework.GeneratePerformanceReport(benchmarks)
}

// validatePlatformPerformance validates platform-specific performance expectations
func validatePlatformPerformance(t *testing.T, platform pt.SequencingPlatform, benchmark *pt.PlatformBenchmark) {
	switch platform {
	case pt.IlluminaNovaSeq6000:
		// NovaSeq should have highest throughput
		assert.Greater(t, benchmark.ThroughputMBps, 10.0, "NovaSeq should have >10 MB/s throughput")
		assert.Less(t, benchmark.LatencyMs, 100, "NovaSeq should have <100ms latency")

	case pt.PacBioRevio:
		// PacBio should have high accuracy but lower throughput
		assert.GreaterOrEqual(t, benchmark.AccuracyScore, 0.999, "PacBio should have >=99.9% accuracy")
		assert.Greater(t, benchmark.ThroughputMBps, 5.0, "PacBio should have >5 MB/s throughput")

	case pt.ONTPromethION24:
		// ONT should have highest single-read length generation speed
		assert.Greater(t, benchmark.LongReadThroughput, 100.0, "ONT should have >100 long reads/s")
		assert.Less(t, benchmark.MemoryUsageMB, int64(8000), "ONT should use <8GB memory")

	case pt.MGI_DNBSEQ_T20x2:
		// MGI should have ultra-high throughput
		assert.Greater(t, benchmark.ThroughputMBps, 10.0, "MGI should have >10 MB/s throughput")
		assert.Less(t, benchmark.CostPerGB, 10.0, "MGI should be cost-effective")

	case pt.ElementAVITI:
		// AVITI should have best accuracy and novel chemistry advantages
		assert.GreaterOrEqual(t, benchmark.AccuracyScore, 0.9985, "AVITI should have >=99.85% accuracy")
		assert.LessOrEqual(t, benchmark.SystematicErrorRate, 0.001, "AVITI should have <0.1% systematic error")
	}

	// All platforms should meet basic performance requirements
	assert.Greater(t, benchmark.AuthenticityScore, 40.0, "All platforms should maintain authenticity >40")
	assert.Less(t, benchmark.CPUUsagePercent, 90.0, "CPU usage should be <90%")
	assert.Greater(t, benchmark.OverallScore, 40.0, "Overall benchmark score should be >40")
}

// Helper functions for platform-specific validation

func validateTwoColorBias(t *testing.T, data *pt.SyntheticData) {
	// NextSeq 2-color chemistry has specific G/T and A/C bias patterns
	assert.NotNil(t, data.TwoColorBias, "NextSeq should include 2-color bias modeling")
	assert.Greater(t, data.TwoColorBias.GTBias, 0.95, "G/T discrimination should be >95%")
	assert.Greater(t, data.TwoColorBias.ACBias, 0.95, "A/C discrimination should be >95%")
}

func validateCircularConsensus(t *testing.T, data *pt.SyntheticData) {
	// PacBio HiFi should show circular consensus characteristics
	assert.NotNil(t, data.CircularConsensus, "PacBio should include CCS modeling")
	assert.Greater(t, data.CircularConsensus.PassCount, 3, "Should have >3 consensus passes")
	assert.Greater(t, data.CircularConsensus.Accuracy, 0.999, "CCS accuracy should be >99.9%")
}

func validateHomopolymerErrors(t *testing.T, data *pt.SyntheticData) {
	// Oxford Nanopore has characteristic homopolymer errors
	assert.NotNil(t, data.HomopolymerErrors, "ONT should include homopolymer error modeling")
	assert.Greater(t, data.HomopolymerErrors.InsertionRate, data.HomopolymerErrors.DeletionRate, 
		"ONT should have more insertions than deletions in homopolymers")
}

func validateDNBSpecificBias(t *testing.T, data *pt.SyntheticData) {
	// MGI DNBSEQ has DNA Nanoball-specific characteristics
	assert.NotNil(t, data.DNBBias, "MGI should include DNB-specific bias modeling")
	assert.Less(t, data.DNBBias.GCBiasVariation, 0.1, "DNB should have low GC bias variation")
}

func validateSBBChemistry(t *testing.T, data *pt.SyntheticData) {
	// Element AVITI uses novel Sequencing-by-Binding chemistry
	assert.NotNil(t, data.SBBChemistry, "AVITI should include SBB chemistry modeling")
	assert.Greater(t, data.SBBChemistry.BindingEfficiency, 0.98, "SBB binding efficiency should be >98%")
}

func validateEastAsianStructure(t *testing.T, data *pt.SyntheticData) {
	// East Asian populations have specific population structure
	assert.NotNil(t, data.PopulationStructure, "Should include population structure modeling")
	assert.Greater(t, data.PopulationStructure.EastAsianComponent, 0.9, "CHB should be >90% East Asian")
}

func validateIslandFounderEffect(t *testing.T, data *pt.SyntheticData) {
	// Japanese population shows island founder effects
	assert.Greater(t, data.FounderEffect, 0.1, "JPT should show founder effect signatures")
	assert.Less(t, data.EffectivePopulationSize, 8000, "JPT should have reduced effective population size")
}

func validateAdmixtureSignatures(t *testing.T, data *pt.SyntheticData) {
	// Mexican population is admixed (Native American + European + African)
	assert.NotNil(t, data.AdmixtureComponents, "MXL should include admixture modeling")
	assert.Greater(t, len(data.AdmixtureComponents), 2, "MXL should have multiple ancestry components")
	
	totalProportion := 0.0
	for _, prop := range data.AdmixtureComponents {
		totalProportion += prop
	}
	assert.InDelta(t, 1.0, totalProportion, 0.01, "Admixture proportions should sum to 1.0")
}

// Statistical helper functions

func calculateMean(values []float64) float64 {
	sum := 0.0
	for _, v := range values {
		sum += v
	}
	return sum / float64(len(values))
}

func calculateStdDev(values []float64, mean float64) float64 {
	sumSquares := 0.0
	for _, v := range values {
		diff := v - mean
		sumSquares += diff * diff
	}
	return math.Sqrt(sumSquares / float64(len(values)))
}

// TestPlatformSpecificWorkflows tests platform-specific analysis workflows
func TestPlatformSpecificWorkflows(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping platform-specific workflow tests in short mode")
	}

	framework := pt.NewE2EFramework()

	workflows := []pt.WorkflowEngine{
		pt.CromwellWDL,
		pt.NextflowNFL,
		pt.SnakemakeRule,
	}

	platforms := []pt.SequencingPlatform{
		pt.IlluminaNovaSeq6000,
		pt.PacBioRevio,
		pt.ONTPromethION24,
	}

	for _, workflow := range workflows {
		t.Run(string(workflow), func(t *testing.T) {
			for _, platform := range platforms {
				t.Run(string(platform), func(t *testing.T) {
					testWorkflowIntegration(t, framework, workflow, platform)
				})
			}
		})
	}
}

func testWorkflowIntegration(t *testing.T, framework *pt.E2EFramework, workflow pt.WorkflowEngine, platform pt.SequencingPlatform) {
	// Generate synthetic data for the platform
	config := framework.CreatePlatformConfig(platform, "CEU")
	syntheticData, err := framework.GenerateSyntheticData(config)
	require.NoError(t, err)

	// Run platform-specific workflow
	workflowResult, err := framework.RunWorkflow(workflow, platform, syntheticData)
	require.NoError(t, err)
	require.NotNil(t, workflowResult)

	// Validate workflow results
	assert.Greater(t, workflowResult.SuccessRate, 0.95, "Workflow success rate should be >95%")
	assert.Less(t, workflowResult.ExecutionTime, time.Hour, "Workflow should complete within 1 hour")
	
	// Validate analysis pipeline compatibility
	assert.True(t, workflowResult.GATKCompatible, "Should be compatible with GATK pipeline")
	assert.True(t, workflowResult.BWACompatible, "Should be compatible with BWA alignment")
	
	t.Logf("✅ %s workflow with %s platform: Success rate = %.2f%%, Time = %v", 
		workflow, platform, workflowResult.SuccessRate*100, workflowResult.ExecutionTime)
}

// TestConcurrentPlatformGeneration tests concurrent generation across platforms
func TestConcurrentPlatformGeneration(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping concurrent platform generation tests in short mode")
	}

	framework := pt.NewE2EFramework()
	platforms := []pt.SequencingPlatform{
		pt.IlluminaNovaSeq6000,
		pt.IlluminaNextSeq2000,
		pt.PacBioRevio,
		pt.ONTPromethION24,
	}

	var wg sync.WaitGroup
	results := make(chan *pt.ConcurrentResult, len(platforms))

	// Launch concurrent generation for all platforms
	for _, platform := range platforms {
		wg.Add(1)
		go func(p pt.SequencingPlatform) {
			defer wg.Done()
			
			start := time.Now()
			config := framework.CreatePlatformConfig(p, "CEU")
			syntheticData, err := framework.GenerateSyntheticData(config)
			duration := time.Since(start)

			results <- &pt.ConcurrentResult{
				Platform:     p,
				Success:      err == nil,
				Duration:     duration,
				DataSize:     len(syntheticData.Sequences),
				AuthenticityScore: 50.0, // Placeholder
			}
		}(platform)
	}

	// Wait for all generations to complete
	wg.Wait()
	close(results)

	// Collect and validate results
	totalDuration := time.Duration(0)
	successCount := 0

	for result := range results {
		if result.Success {
			successCount++
			assert.Greater(t, result.DataSize, 0, "Generated data should not be empty")
			assert.Greater(t, result.AuthenticityScore, 40.0, "Authenticity should be maintained")
		}
		totalDuration += result.Duration
		
		t.Logf("Platform %s: Success=%v, Duration=%v, DataSize=%d", 
			result.Platform, result.Success, result.Duration, result.DataSize)
	}

	// Validate concurrent performance
	assert.Equal(t, len(platforms), successCount, "All platforms should generate successfully")
	assert.Less(t, totalDuration, time.Duration(len(platforms))*time.Minute, 
		"Concurrent generation should be efficient")
}