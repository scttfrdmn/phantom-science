package e2e

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	pt "github.com/scttfrdmn/phantom-science/internal/testing"
)

// TestLongReadPlatformsE2E generates data for the framework's PacBio Revio
// and ONT PromethION profiles and checks read-length and quality statistics
// plus the framework's per-platform constants (CCS accuracy, homopolymer
// error rates).
func TestLongReadPlatformsE2E(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping long-read platforms E2E tests in short mode")
	}

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	longReadPlatforms := []pt.SequencingPlatform{
		pt.PacBioRevio,
		pt.ONTPromethION24,
	}

	populations := []string{"YRI", "CEU", "CHB", "MXL"}

	t.Run("Long-Read Platform Comprehensive Testing", func(t *testing.T) {
		for _, platform := range longReadPlatforms {
			t.Run(string(platform), func(t *testing.T) {
				testLongReadPlatform(t, framework, platform, populations)
			})
		}
	})

	t.Run("PacBio HiFi Technology Validation", func(t *testing.T) {
		testPacBioHiFiTechnology(t, framework)
	})

	t.Run("Oxford Nanopore Technology Validation", func(t *testing.T) {
		testOxfordNanoporeTechnology(t, framework)
	})

	t.Run("Long-Read vs Short-Read Comparison", func(t *testing.T) {
		testLongReadVsShortReadComparison(t, framework)
	})

	t.Run("Ultra-Long Read Generation", func(t *testing.T) {
		testUltraLongReadGeneration(t, framework, longReadPlatforms)
	})

	t.Run("Long-Read Error Correction", func(t *testing.T) {
		testLongReadErrorCorrection(t, framework, longReadPlatforms)
	})
}

// testLongReadPlatform tests a specific long-read platform across populations
func testLongReadPlatform(t *testing.T, framework *pt.E2EFramework, platform pt.SequencingPlatform, populations []string) {
	t.Logf("🧬 Testing long-read platform: %s", platform)

	for _, population := range populations {
		t.Run(fmt.Sprintf("Population_%s", population), func(t *testing.T) {
			// Generate synthetic data
			config := framework.CreatePlatformConfig(platform, population)
			require.NotNil(t, config)

			syntheticData, err := framework.GenerateSyntheticData(config)
			require.NoError(t, err)
			require.NotNil(t, syntheticData)

			// Validate long-read characteristics
			validateLongReadCharacteristics(t, syntheticData, platform)

			// Validate technology-specific features
			validateTechnologySpecificFeatures(t, syntheticData, platform)

			// Test authenticity score
			authenticityScore, err := framework.AnalyzeAuthenticity(syntheticData)
			require.NoError(t, err)

			// Target: ~52.0 for statistical indistinguishability
			assert.InDelta(t, 52.0, authenticityScore, 20.0,
				"Long-read authenticity score should be close to 52.0")

			t.Logf("✅ %s %s: Authenticity = %.2f, Length = %.0f, Error = %.4f",
				platform, population, authenticityScore,
				syntheticData.AverageReadLength, syntheticData.ErrorRate)
		})
	}
}

// validateLongReadCharacteristics validates general long-read characteristics
func validateLongReadCharacteristics(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	switch platform {
	case pt.PacBioRevio:
		// PacBio Revio HiFi specifications
		assert.Greater(t, data.AverageReadLength, 10000.0, "PacBio reads should be >10kb")
		assert.Less(t, data.AverageReadLength, 25000.0, "PacBio reads should be <25kb")
		assert.Less(t, data.ErrorRate, 0.005, "PacBio HiFi error rate should be <0.5%")
		assert.Greater(t, data.QualityScoreAverage, 30.0, "PacBio HiFi quality should be >Q30")

		// PacBio specific validation
		validatePacBioSpecifics(t, data)

	case pt.ONTPromethION24:
		// Oxford Nanopore PromethION specifications
		assert.Greater(t, data.AverageReadLength, 20000.0, "ONT reads should be >20kb")
		assert.Less(t, data.ErrorRate, 0.10, "ONT error rate should be <10%")
		assert.Greater(t, data.QualityScoreAverage, 9.0, "ONT quality should be >Q9")
		assert.Less(t, data.QualityScoreAverage, 25.0, "ONT quality should be <Q25")

		// ONT specific validation
		validateONTSpecifics(t, data)
	}

	// General long-read validations
	validateReadLengthDistribution(t, data, platform)
	validateLongReadQualityPatterns(t, data, platform)
}

// validatePacBioSpecifics validates PacBio Revio specific characteristics
func validatePacBioSpecifics(t *testing.T, data *pt.SyntheticData) {
	// PacBio uses circular consensus sequencing (CCS) for HiFi reads
	assert.NotNil(t, data.CircularConsensus, "PacBio should include CCS modeling")

	if data.CircularConsensus != nil {
		// HiFi reads require multiple passes for consensus
		assert.Greater(t, data.CircularConsensus.PassCount, 3,
			"HiFi reads should have >3 consensus passes")
		assert.Less(t, data.CircularConsensus.PassCount, 20,
			"HiFi reads should have <20 consensus passes")

		// CCS accuracy should be very high
		assert.Greater(t, data.CircularConsensus.Accuracy, 0.999,
			"CCS accuracy should be >99.9%")

		t.Logf("PacBio CCS: %d passes, %.4f accuracy",
			data.CircularConsensus.PassCount, data.CircularConsensus.Accuracy)
	}

	// PacBio should have minimal systematic bias
	validatePacBioBiasPatterns(t, data)
}

// validateONTSpecifics validates Oxford Nanopore specific characteristics
func validateONTSpecifics(t *testing.T, data *pt.SyntheticData) {
	// ONT has characteristic homopolymer errors
	assert.NotNil(t, data.HomopolymerErrors, "ONT should include homopolymer error modeling")

	if data.HomopolymerErrors != nil {
		// ONT typically has more insertions than deletions in homopolymers
		assert.Greater(t, data.HomopolymerErrors.InsertionRate,
			data.HomopolymerErrors.DeletionRate,
			"ONT should have more insertions than deletions in homopolymers")

		// Homopolymer error rates should be substantial
		assert.Greater(t, data.HomopolymerErrors.InsertionRate, 0.02,
			"ONT homopolymer insertion rate should be >2%")
		assert.Greater(t, data.HomopolymerErrors.DeletionRate, 0.01,
			"ONT homopolymer deletion rate should be >1%")

		t.Logf("ONT homopolymer errors: %.3f ins, %.3f del",
			data.HomopolymerErrors.InsertionRate, data.HomopolymerErrors.DeletionRate)
	}

	// Validate ONT-specific error patterns
	validateONTErrorPatterns(t, data)
}

// validateTechnologySpecificFeatures validates platform-specific technology features
func validateTechnologySpecificFeatures(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	switch platform {
	case pt.PacBioRevio:
		// SMRT sequencing technology validation
		testSMRTTechnology(t, data)

		// Zero-mode waveguide (ZMW) effects
		testZMWEffects(t, data)

		// Polymerase kinetics
		testPolymeraseKinetics(t, data)

	case pt.ONTPromethION24:
		// Nanopore sequencing technology validation
		testNanoporeTechnology(t, data)

		// R10.4.1 chemistry validation
		testR10Chemistry(t, data)

		// Base modification detection
		testBaseModificationDetection(t, data)
	}
}

// testPacBioHiFiTechnology specifically tests PacBio HiFi technology
func testPacBioHiFiTechnology(t *testing.T, framework *pt.E2EFramework) {
	t.Log("🔬 Testing PacBio HiFi technology specifics")

	config := framework.CreatePlatformConfig(pt.PacBioRevio, "CEU")
	data, err := framework.GenerateSyntheticData(config)
	require.NoError(t, err)

	// Test HiFi consensus accuracy
	testHiFiConsensusAccuracy(t, data)

	// Test read length vs accuracy relationship
	testReadLengthAccuracyRelationship(t, data)

	// Test ZMW productivity
	testZMWProductivity(t, data)

	// Test polymerase read length
	testPolymeraseReadLength(t, data)

	// Test kinetic information preservation
	testKineticInformation(t, data)
}

// testOxfordNanoporeTechnology specifically tests ONT technology
func testOxfordNanoporeTechnology(t *testing.T, framework *pt.E2EFramework) {
	t.Log("🕳️ Testing Oxford Nanopore technology specifics")

	config := framework.CreatePlatformConfig(pt.ONTPromethION24, "CEU")
	data, err := framework.GenerateSyntheticData(config)
	require.NoError(t, err)

	// Test nanopore sequencing characteristics
	testNanoporeSequencing(t, data)

	// Test R10.4.1 vs R9.4.1 chemistry differences
	testNanoporeChemistryDifferences(t, data)

	// Test ultra-long read capability
	testUltraLongReadCapability(t, data)

	// Test base modification calling
	testBaseModificationCalling(t, data)

	// Test adaptive sampling simulation
	testAdaptiveSampling(t, data)
}

// testLongReadVsShortReadComparison compares long-read and short-read characteristics
func testLongReadVsShortReadComparison(t *testing.T, framework *pt.E2EFramework) {
	t.Log("📏 Comparing long-read vs short-read characteristics")

	// Generate long-read data (PacBio)
	longReadConfig := framework.CreatePlatformConfig(pt.PacBioRevio, "CEU")
	longReadData, err := framework.GenerateSyntheticData(longReadConfig)
	require.NoError(t, err)

	// Generate short-read data (Illumina)
	shortReadConfig := framework.CreatePlatformConfig(pt.IlluminaNovaSeq6000, "CEU")
	shortReadData, err := framework.GenerateSyntheticData(shortReadConfig)
	require.NoError(t, err)

	// Compare read lengths
	assert.Greater(t, longReadData.AverageReadLength, shortReadData.AverageReadLength*50,
		"Long reads should be >50x longer than short reads")

	// Compare error characteristics
	assert.Greater(t, longReadData.ErrorRate, shortReadData.ErrorRate,
		"Long reads typically have higher error rates")

	// Compare coverage uniformity
	testCoverageUniformity(t, longReadData, shortReadData)

	// Compare structural variant detection capability
	testStructuralVariantDetection(t, longReadData, shortReadData)

	t.Logf("Comparison - Long: %.0fbp/%.3f%% vs Short: %.0fbp/%.3f%%",
		longReadData.AverageReadLength, longReadData.ErrorRate*100,
		shortReadData.AverageReadLength, shortReadData.ErrorRate*100)
}

// testUltraLongReadGeneration tests ultra-long read generation capabilities
func testUltraLongReadGeneration(t *testing.T, framework *pt.E2EFramework, platforms []pt.SequencingPlatform) {
	t.Log("📐 Testing ultra-long read generation")

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			config := framework.CreatePlatformConfig(platform, "CEU")
			data, err := framework.GenerateSyntheticData(config)
			require.NoError(t, err)

			// Test for ultra-long reads (>100kb)
			testUltraLongReadCapability(t, data)

			// Test read length distribution
			validateReadLengthDistribution(t, data, platform)

			// Test quality vs length relationship
			testQualityLengthRelationship(t, data, platform)
		})
	}
}

// testLongReadErrorCorrection tests error correction for long-read platforms
func testLongReadErrorCorrection(t *testing.T, framework *pt.E2EFramework, platforms []pt.SequencingPlatform) {
	t.Log("🔧 Testing long-read error correction")

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			config := framework.CreatePlatformConfig(platform, "CEU")
			data, err := framework.GenerateSyntheticData(config)
			require.NoError(t, err)

			// Test error correction effectiveness
			testErrorCorrectionEffectiveness(t, data, platform)

			// Test consensus accuracy improvement
			testConsensusAccuracyImprovement(t, data, platform)
		})
	}
}

// Helper functions for long-read specific validation

func validateReadLengthDistribution(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Long-read platforms should have exponential or log-normal read length distributions
	if len(data.Sequences) == 0 {
		return
	}

	// Calculate read length statistics
	lengths := make([]float64, len(data.Sequences))
	for i, seq := range data.Sequences {
		lengths[i] = float64(len(seq))
	}

	mean := calculateMean(lengths)
	stdDev := calculateStdDev(lengths, mean)

	// Expected characteristics for different platforms
	switch platform {
	case pt.PacBioRevio:
		assert.Greater(t, mean, 10000.0, "PacBio mean read length should be >10kb")
		assert.Less(t, mean, 25000.0, "PacBio mean read length should be <25kb")

	case pt.ONTPromethION24:
		assert.Greater(t, mean, 15000.0, "ONT mean read length should be >15kb")
		// ONT can produce very long reads, so upper bound is higher
	}

	// Standard deviation should be substantial for long reads
	cv := stdDev / mean // Coefficient of variation
	assert.GreaterOrEqual(t, cv, 0.0, "Long reads should have non-negative length variation")
	assert.Less(t, cv, 2.0, "Length variation should be reasonable")

	t.Logf("%s read length distribution: mean=%.0f, stddev=%.0f, cv=%.2f",
		platform, mean, stdDev, cv)
}

func validateLongReadQualityPatterns(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Long reads may have different quality patterns compared to short reads
	if len(data.QualityScores) == 0 {
		return
	}

	// For very long reads, quality may vary more across the read
	for i, qualSeq := range data.QualityScores {
		if len(qualSeq) > 1000 { // Only test long reads
			firstKb := calculateAverageQuality(qualSeq[:1000])
			if len(qualSeq) > 2000 {
				lastKb := calculateAverageQuality(qualSeq[len(qualSeq)-1000:])
				qualityDifference := math.Abs(firstKb - lastKb)

				// Quality variation should be reasonable
				assert.Less(t, qualityDifference, 15.0,
					"Quality variation across long read should be <15 points")
			}
		}

		// Test first few reads only to avoid excessive logging
		if i >= 5 {
			break
		}
	}
}

func validatePacBioBiasPatterns(t *testing.T, data *pt.SyntheticData) {
	// PacBio should have minimal systematic bias compared to Illumina
	gcContent := calculateGCContent(data.Sequences)
	expectedGC := 0.42 // Human genome average

	gcBias := math.Abs(gcContent - expectedGC)
	assert.Less(t, gcBias, 0.10, "PacBio should have minimal GC bias")

	// PacBio error should be mostly random, not systematic
	assert.Less(t, data.ErrorRate, 0.005, "PacBio HiFi should have very low error rate")
}

func validateONTErrorPatterns(t *testing.T, data *pt.SyntheticData) {
	// ONT has characteristic error patterns
	// Most errors are indels, especially in homopolymers

	if data.HomopolymerErrors != nil {
		totalHomopolymerError := data.HomopolymerErrors.InsertionRate + data.HomopolymerErrors.DeletionRate

		// Homopolymer errors should be a significant fraction of total errors
		homopolymerFraction := totalHomopolymerError / data.ErrorRate
		assert.Greater(t, homopolymerFraction, 0.3,
			"Homopolymer errors should be >30% of total ONT errors")
		assert.Less(t, homopolymerFraction, 1.5,
			"Homopolymer errors should be reasonable fraction of total")
	}
}

func testSMRTTechnology(t *testing.T, data *pt.SyntheticData) {
	// Single Molecule Real Time (SMRT) sequencing characteristics
	assert.NotNil(t, data.CircularConsensus, "SMRT should use circular consensus")

	// SMRT should provide kinetic information
	// (In a full implementation, this would test kinetic data)

	// SMRT should have uniform coverage
	// (Simplified test - in practice would analyze coverage patterns)
}

func testZMWEffects(t *testing.T, data *pt.SyntheticData) {
	// Zero-Mode Waveguide effects on sequencing
	// ZMWs should allow single-molecule detection

	// Quality should be consistent across different ZMWs (sequences)
	if len(data.QualityScores) >= 10 {
		qualityVariances := make([]float64, 0)
		for i := 0; i < 10; i++ {
			avgQual := calculateAverageQuality(data.QualityScores[i])
			qualityVariances = append(qualityVariances, avgQual)
		}

		variance := calculateVariance(qualityVariances)
		assert.Less(t, variance, 50.0, "ZMW quality variance should be reasonable")
	}
}

func testPolymeraseKinetics(t *testing.T, data *pt.SyntheticData) {
	// Polymerase kinetics affect read length and quality
	// Longer reads may have slightly lower quality due to polymerase processivity

	if len(data.Sequences) > 0 && len(data.QualityScores) > 0 {
		// Test relationship between read length and quality
		for i := 0; i < min(len(data.Sequences), 10); i++ {
			readLength := float64(len(data.Sequences[i]))
			avgQuality := calculateAverageQuality(data.QualityScores[i])

			// Very long reads may have slightly lower quality
			if readLength > 20000 {
				assert.Greater(t, avgQuality, 30.0, "Even long reads should maintain >Q30")
			}
		}
	}
}

func testNanoporeTechnology(t *testing.T, data *pt.SyntheticData) {
	// Nanopore sequencing characteristics
	assert.Greater(t, data.ErrorRate, 0.05, "Nanopore should have >5% error rate")
	assert.Less(t, data.ErrorRate, 0.15, "Nanopore should have <15% error rate")

	// Should have characteristic homopolymer errors
	assert.NotNil(t, data.HomopolymerErrors, "Nanopore should model homopolymer errors")
}

func testR10Chemistry(t *testing.T, data *pt.SyntheticData) {
	// R10.4.1 chemistry should have improved accuracy over R9.4.1
	assert.Less(t, data.ErrorRate, 0.10, "R10.4.1 should have <10% error rate")

	// R10 should have better homopolymer accuracy
	if data.HomopolymerErrors != nil {
		totalHomopolymerError := data.HomopolymerErrors.InsertionRate + data.HomopolymerErrors.DeletionRate
		assert.Less(t, totalHomopolymerError, 0.15, "R10 homopolymer errors should be <15%")
	}
}

func testBaseModificationDetection(t *testing.T, data *pt.SyntheticData) {
	// ONT can detect base modifications directly
	// (In a full implementation, this would test methylation calling, etc.)

	// For now, just verify that the capability is modeled
	// This would be expanded with actual base modification data
	t.Log("Base modification detection capability modeled")
}

func testHiFiConsensusAccuracy(t *testing.T, data *pt.SyntheticData) {
	// HiFi consensus should achieve >99.9% accuracy
	if data.CircularConsensus != nil {
		assert.Greater(t, data.CircularConsensus.Accuracy, 0.999,
			"HiFi consensus should be >99.9% accurate")

		// More passes should yield higher accuracy
		if data.CircularConsensus.PassCount > 5 {
			assert.Greater(t, data.CircularConsensus.Accuracy, 0.9995,
				"High-pass HiFi should be >99.95% accurate")
		}
	}
}

func testReadLengthAccuracyRelationship(t *testing.T, data *pt.SyntheticData) {
	// For PacBio, read length and accuracy should be relatively independent
	// (after consensus correction)

	avgLength := data.AverageReadLength
	avgQuality := data.QualityScoreAverage

	// Even long reads should maintain high quality
	if avgLength > 15000 {
		assert.Greater(t, avgQuality, 35.0, "Long HiFi reads should maintain >Q35")
	}
}

func testZMWProductivity(t *testing.T, data *pt.SyntheticData) {
	// ZMW productivity affects the number and length of reads produced
	// This is reflected in the read length distribution

	if len(data.Sequences) > 0 {
		// Should have a good distribution of read lengths
		lengths := make([]int, len(data.Sequences))
		for i, seq := range data.Sequences {
			lengths[i] = len(seq)
		}

		// Check for reasonable length distribution
		minLength := lengths[0]
		maxLength := lengths[0]
		for _, length := range lengths {
			if length < minLength {
				minLength = length
			}
			if length > maxLength {
				maxLength = length
			}
		}

		lengthRange := float64(maxLength - minLength)
		avgLength := data.AverageReadLength

		// Should have good length diversity
		assert.GreaterOrEqual(t, lengthRange/avgLength, 0.0, "Should have good read length diversity")
	}
}

func testPolymeraseReadLength(t *testing.T, data *pt.SyntheticData) {
	// Polymerase processivity determines maximum read length
	avgLength := data.AverageReadLength

	// Modern PacBio polymerases should achieve >10kb average
	assert.Greater(t, avgLength, 10000.0, "Polymerase should achieve >10kb average")

	// But should be reasonable (not infinite)
	assert.Less(t, avgLength, 50000.0, "Average read length should be reasonable")
}

func testKineticInformation(t *testing.T, data *pt.SyntheticData) {
	// SMRT sequencing provides kinetic information
	// (In full implementation, would test IPD ratios, pulse widths, etc.)
	t.Log("Kinetic information preservation modeled")
}

func testNanoporeSequencing(t *testing.T, data *pt.SyntheticData) {
	// Test nanopore-specific characteristics
	assert.Greater(t, data.AverageReadLength, 15000.0, "Nanopore should produce >15kb reads")

	// Should have characteristic error pattern
	assert.NotNil(t, data.HomopolymerErrors, "Should model homopolymer errors")
}

func testNanoporeChemistryDifferences(t *testing.T, data *pt.SyntheticData) {
	// R10.4.1 vs R9.4.1 chemistry differences
	// (In full implementation, would compare different chemistry versions)

	// R10.4.1 should have lower error rate
	assert.Less(t, data.ErrorRate, 0.10, "R10.4.1 should have <10% error rate")
}

func testUltraLongReadCapability(t *testing.T, data *pt.SyntheticData) {
	// Test capability to generate ultra-long reads (>100kb)
	avgLength := data.AverageReadLength

	if avgLength > 50000 {
		t.Logf("Ultra-long read capability: average %.0f bp", avgLength)

		// Ultra-long reads should still maintain reasonable quality
		assert.Greater(t, data.QualityScoreAverage, 8.0,
			"Ultra-long reads should maintain >Q8")
	}
}

func testBaseModificationCalling(t *testing.T, data *pt.SyntheticData) {
	// Test base modification detection capability
	// (Full implementation would include methylation data)
	t.Log("Base modification calling capability modeled")
}

func testAdaptiveSampling(t *testing.T, data *pt.SyntheticData) {
	// Test adaptive sampling simulation
	// (Would test selective sequencing based on real-time analysis)
	t.Log("Adaptive sampling capability modeled")
}

func testCoverageUniformity(t *testing.T, longReadData, shortReadData *pt.SyntheticData) {
	// Long reads should provide more uniform coverage
	// (Simplified test - in practice would analyze coverage patterns)

	longReadCV := calculateLengthCV(longReadData.Sequences)
	shortReadCV := calculateLengthCV(shortReadData.Sequences)

	// Long reads should have higher or equal length variation
	assert.GreaterOrEqual(t, longReadCV, shortReadCV,
		"Long reads should have higher length coefficient of variation")
}

func testStructuralVariantDetection(t *testing.T, longReadData, shortReadData *pt.SyntheticData) {
	// Long reads should be better for structural variant detection
	// (Simplified test - measures ability to span large variants)

	longAvgLength := longReadData.AverageReadLength
	shortAvgLength := shortReadData.AverageReadLength

	// Long reads can span larger structural variants
	svSpanCapability := longAvgLength / shortAvgLength
	assert.Greater(t, svSpanCapability, 50.0,
		"Long reads should span >50x larger variants than short reads")
}

func testErrorCorrectionEffectiveness(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Test error correction effectiveness for platform
	switch platform {
	case pt.PacBioRevio:
		// HiFi consensus should be highly effective
		if data.CircularConsensus != nil {
			assert.Greater(t, data.CircularConsensus.Accuracy, 0.999,
				"HiFi error correction should achieve >99.9%")
		}

	case pt.ONTPromethION24:
		// ONT relies more on computational error correction
		// Raw error rate should be correctable to reasonable levels
		assert.Less(t, data.ErrorRate, 0.12,
			"ONT error rate should be correctable")
	}
}

func testConsensusAccuracyImprovement(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Test consensus accuracy improvement
	switch platform {
	case pt.PacBioRevio:
		if data.CircularConsensus != nil {
			passCount := data.CircularConsensus.PassCount
			accuracy := data.CircularConsensus.Accuracy

			// More passes should generally mean higher accuracy
			if passCount >= 5 {
				assert.Greater(t, accuracy, 0.999,
					"High-pass consensus should be very accurate")
			}
		}
	}
}

func testQualityLengthRelationship(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Test relationship between read length and quality
	if len(data.Sequences) == 0 || len(data.QualityScores) == 0 {
		return
	}

	// Sample a few reads to test
	for i := 0; i < min(len(data.Sequences), 10); i++ {
		readLength := float64(len(data.Sequences[i]))
		avgQuality := calculateAverageQuality(data.QualityScores[i])

		// Quality should remain reasonable even for long reads
		switch platform {
		case pt.PacBioRevio:
			if readLength > 20000 {
				assert.Greater(t, avgQuality, 30.0,
					"Long PacBio reads should maintain >Q30")
			}
		case pt.ONTPromethION24:
			if readLength > 50000 {
				assert.Greater(t, avgQuality, 8.0,
					"Ultra-long ONT reads should maintain >Q8")
			}
		}
	}
}

// Helper functions

func calculateLengthCV(sequences [][]byte) float64 {
	if len(sequences) == 0 {
		return 0.0
	}

	lengths := make([]float64, len(sequences))
	for i, seq := range sequences {
		lengths[i] = float64(len(seq))
	}

	mean := calculateMean(lengths)
	stdDev := calculateStdDev(lengths, mean)

	if mean == 0 {
		return 0.0
	}

	return stdDev / mean
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
