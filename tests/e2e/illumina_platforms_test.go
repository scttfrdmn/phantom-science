package e2e

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	pt "github.com/scttfrdmn/phantom-science/internal/testing"
)

// TestIlluminaPlatformsE2E generates data for the framework's Illumina
// platforms (NovaSeq 6000, NextSeq 2000, MiSeq) and checks read lengths,
// quality and GC statistics plus the framework's per-platform constants.
func TestIlluminaPlatformsE2E(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping Illumina platforms E2E tests in short mode")
	}

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	illuminaPlatforms := []pt.SequencingPlatform{
		pt.IlluminaNovaSeq6000,
		pt.IlluminaNextSeq2000,
		pt.IlluminaMiSeq,
	}

	populations := []string{"YRI", "CEU", "CHB", "MXL"}

	t.Run("Illumina Platform Comprehensive Testing", func(t *testing.T) {
		for _, platform := range illuminaPlatforms {
			t.Run(string(platform), func(t *testing.T) {
				testIlluminaPlatform(t, framework, platform, populations)
			})
		}
	})

	t.Run("SBS Chemistry Validation", func(t *testing.T) {
		testSBSChemistryValidation(t, framework, illuminaPlatforms)
	})

	t.Run("Cluster Density and Phasing Effects", func(t *testing.T) {
		testClusterDensityAndPhasing(t, framework, illuminaPlatforms)
	})

	t.Run("Two-Color vs Four-Color Chemistry", func(t *testing.T) {
		testColorChemistryDifferences(t, framework)
	})

	t.Run("Illumina Quality Score Calibration", func(t *testing.T) {
		testIlluminaQualityCalibration(t, framework, illuminaPlatforms)
	})
}

// testIlluminaPlatform tests a specific Illumina platform across populations
func testIlluminaPlatform(t *testing.T, framework *pt.E2EFramework, platform pt.SequencingPlatform, populations []string) {
	t.Logf("🧬 Testing Illumina platform: %s", platform)

	for _, population := range populations {
		t.Run(fmt.Sprintf("Population_%s", population), func(t *testing.T) {
			// Generate synthetic data
			config := framework.CreatePlatformConfig(platform, population)
			require.NotNil(t, config)

			syntheticData, err := framework.GenerateSyntheticData(config)
			require.NoError(t, err)
			require.NotNil(t, syntheticData)

			// Validate Illumina-specific characteristics
			validateIlluminaCharacteristics(t, syntheticData, platform)

			// Validate SBS chemistry effects
			validateSBSChemistry(t, syntheticData, platform)

			// Test authenticity score
			authenticityScore, err := framework.AnalyzeAuthenticity(syntheticData)
			require.NoError(t, err)

			// Target: ~52.0 for statistical indistinguishability
			assert.InDelta(t, 52.0, authenticityScore, 8.0,
				"Illumina authenticity score should be close to 52.0")

			t.Logf("✅ %s %s: Authenticity = %.2f, Quality = %.1f, Error = %.4f",
				platform, population, authenticityScore,
				syntheticData.QualityScoreAverage, syntheticData.ErrorRate)
		})
	}
}

// validateIlluminaCharacteristics validates Illumina-specific data characteristics
func validateIlluminaCharacteristics(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	switch platform {
	case pt.IlluminaNovaSeq6000:
		// NovaSeq 6000 specifications
		assert.InDelta(t, 150, data.AverageReadLength, 10, "NovaSeq read length should be ~150bp")
		assert.Less(t, data.ErrorRate, 0.002, "NovaSeq error rate should be <0.2%")
		assert.Greater(t, data.QualityScoreAverage, 30.0, "NovaSeq quality should be >Q30")
		assert.Less(t, data.QualityScoreAverage, 42.0, "NovaSeq quality should be <Q42")

		// NovaSeq specific validation
		validateNovaSeqSpecifics(t, data)

	case pt.IlluminaNextSeq2000:
		// NextSeq 2000 specifications
		assert.InDelta(t, 150, data.AverageReadLength, 10, "NextSeq read length should be ~150bp")
		assert.Less(t, data.ErrorRate, 0.003, "NextSeq error rate should be <0.3%")
		assert.Greater(t, data.QualityScoreAverage, 28.0, "NextSeq quality should be >Q28")

		// NextSeq 2-color chemistry validation
		validateNextSeq2ColorChemistry(t, data)

	case pt.IlluminaMiSeq:
		// MiSeq specifications
		assert.InDelta(t, 300, data.AverageReadLength, 50, "MiSeq read length should be ~300bp")
		assert.LessOrEqual(t, data.ErrorRate, 0.001, "MiSeq error rate should be <0.1%")
		assert.Greater(t, data.QualityScoreAverage, 28.0, "MiSeq quality should be >Q28")

		// MiSeq 4-color chemistry validation
		validateMiSeq4ColorChemistry(t, data)
	}
}

// validateNovaSeqSpecifics validates NovaSeq 6000 specific characteristics
func validateNovaSeqSpecifics(t *testing.T, data *pt.SyntheticData) {
	// NovaSeq uses exclusion amplification (ExAmp) for cluster generation
	// Should have high cluster density with minimal cross-talk

	// Validate quality drop pattern (typical for bridge amplification)
	if len(data.QualityScores) > 0 && len(data.QualityScores[0]) > 100 {
		firstHalfQuality := calculateAverageQuality(data.QualityScores[0][:len(data.QualityScores[0])/2])
		secondHalfQuality := calculateAverageQuality(data.QualityScores[0][len(data.QualityScores[0])/2:])

		// Quality should decline slightly towards the end
		assert.Greater(t, firstHalfQuality, secondHalfQuality-5.0,
			"NovaSeq quality decline should be minimal")
	}

	// NovaSeq should have consistent error patterns
	assert.LessOrEqual(t, data.ErrorRate, 0.0015, "NovaSeq should have very low error rate")
}

// validateNextSeq2ColorChemistry validates NextSeq 2000 2-color chemistry
func validateNextSeq2ColorChemistry(t *testing.T, data *pt.SyntheticData) {
	// NextSeq uses 2-color chemistry: Green (G/T) and Red (A/C)
	assert.NotNil(t, data.TwoColorBias, "NextSeq should have 2-color bias modeling")

	if data.TwoColorBias != nil {
		// G/T discrimination should be high
		assert.Greater(t, data.TwoColorBias.GTBias, 0.95,
			"G/T discrimination should be >95%")

		// A/C discrimination should be high
		assert.Greater(t, data.TwoColorBias.ACBias, 0.95,
			"A/C discrimination should be >95%")

		// Validate 2-color specific error patterns
		validateTwoColorErrorPatterns(t, data)
	}
}

// validateMiSeq4ColorChemistry validates MiSeq 4-color chemistry
func validateMiSeq4ColorChemistry(t *testing.T, data *pt.SyntheticData) {
	// MiSeq uses 4-color chemistry with individual nucleotide detection
	// Should have the highest accuracy due to distinct signals

	assert.Greater(t, data.QualityScoreAverage, 28.0,
		"MiSeq 4-color chemistry should achieve >Q28 average")

	// 4-color chemistry should have minimal systematic bias
	gcContent := calculateGCContent(data.Sequences)
	assert.InDelta(t, 0.42, gcContent, 0.10,
		"MiSeq should have minimal GC bias")
}

// validateSBSChemistry validates Sequencing by Synthesis chemistry effects
func validateSBSChemistry(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// All Illumina platforms use SBS chemistry
	// Validate typical SBS characteristics

	// SBS should show position-dependent quality decline
	if len(data.QualityScores) > 0 && len(data.QualityScores[0]) > 50 {
		early := calculateAverageQuality(data.QualityScores[0][:25])
		late := calculateAverageQuality(data.QualityScores[0][len(data.QualityScores[0])-25:])

		// Quality decline should be present but moderate
		qualityDecline := early - late
		assert.Greater(t, qualityDecline, -2.0, "Quality decline should be minimal")
		assert.Less(t, qualityDecline, 16.0, "Quality decline should be reasonable")

		t.Logf("%s SBS quality decline: %.1f (early=%.1f, late=%.1f)",
			platform, qualityDecline, early, late)
	}

	// SBS should have characteristic error patterns
	validateSBSErrorPatterns(t, data, platform)
}

// testSBSChemistryValidation tests SBS chemistry across Illumina platforms
func testSBSChemistryValidation(t *testing.T, framework *pt.E2EFramework, platforms []pt.SequencingPlatform) {
	t.Log("🧪 Testing SBS chemistry validation across Illumina platforms")

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			config := framework.CreatePlatformConfig(platform, "CEU")
			data, err := framework.GenerateSyntheticData(config)
			require.NoError(t, err)

			// Validate SBS-specific characteristics
			validateSBSChemistry(t, data, platform)

			// Test incorporation efficiency
			testIncorporationEfficiency(t, data, platform)

			// Test phasing and prephasing effects
			testPhasingEffects(t, data, platform)
		})
	}
}

// testClusterDensityAndPhasing tests cluster density and phasing effects
func testClusterDensityAndPhasing(t *testing.T, framework *pt.E2EFramework, platforms []pt.SequencingPlatform) {
	t.Log("🎯 Testing cluster density and phasing effects")

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			config := framework.CreatePlatformConfig(platform, "CEU")
			data, err := framework.GenerateSyntheticData(config)
			require.NoError(t, err)

			// Test cluster density effects
			testClusterDensityEffects(t, data, platform)

			// Test phasing and prephasing
			testPhasingAndPrephasing(t, data, platform)

			// Test signal intensity consistency
			testSignalIntensity(t, data, platform)
		})
	}
}

// testColorChemistryDifferences tests differences between 2-color and 4-color chemistry
func testColorChemistryDifferences(t *testing.T, framework *pt.E2EFramework) {
	t.Log("🌈 Comparing 2-color vs 4-color chemistry differences")

	// Generate data for NextSeq (2-color) and MiSeq (4-color)
	nextseqConfig := framework.CreatePlatformConfig(pt.IlluminaNextSeq2000, "CEU")
	nextseqData, err := framework.GenerateSyntheticData(nextseqConfig)
	require.NoError(t, err)

	miseqConfig := framework.CreatePlatformConfig(pt.IlluminaMiSeq, "CEU")
	miseqData, err := framework.GenerateSyntheticData(miseqConfig)
	require.NoError(t, err)

	// Compare accuracy
	assert.Greater(t, miseqData.QualityScoreAverage, nextseqData.QualityScoreAverage,
		"4-color chemistry (MiSeq) should have higher quality than 2-color (NextSeq)")

	// Compare error rates
	assert.Less(t, miseqData.ErrorRate, nextseqData.ErrorRate,
		"4-color chemistry should have lower error rate")

	// Test chemistry-specific biases
	testChemistrySpecificBiases(t, nextseqData, miseqData)

	t.Logf("Chemistry comparison - NextSeq Q%.1f/%.4f%% vs MiSeq Q%.1f/%.4f%%",
		nextseqData.QualityScoreAverage, nextseqData.ErrorRate*100,
		miseqData.QualityScoreAverage, miseqData.ErrorRate*100)
}

// testIlluminaQualityCalibration tests Illumina quality score calibration
func testIlluminaQualityCalibration(t *testing.T, framework *pt.E2EFramework, platforms []pt.SequencingPlatform) {
	t.Log("📊 Testing Illumina quality score calibration")

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			config := framework.CreatePlatformConfig(platform, "CEU")
			data, err := framework.GenerateSyntheticData(config)
			require.NoError(t, err)

			// Test quality score distribution
			testQualityScoreDistribution(t, data, platform)

			// Test quality vs error correlation
			testQualityErrorCorrelation(t, data, platform)

			// Test Phred score accuracy
			testPhredScoreAccuracy(t, data, platform)
		})
	}
}

// Helper functions for Illumina-specific validation

func calculateAverageQuality(qualityBytes []byte) float64 {
	if len(qualityBytes) == 0 {
		return 0.0
	}

	total := 0.0
	for _, q := range qualityBytes {
		total += float64(q - 33) // Convert from Phred+33
	}

	return total / float64(len(qualityBytes))
}

func calculateGCContent(sequences [][]byte) float64 {
	totalBases := 0
	gcBases := 0

	for _, seq := range sequences {
		for _, base := range seq {
			totalBases++
			if base == 'G' || base == 'C' {
				gcBases++
			}
		}
	}

	if totalBases == 0 {
		return 0.0
	}

	return float64(gcBases) / float64(totalBases)
}

func validateTwoColorErrorPatterns(t *testing.T, data *pt.SyntheticData) {
	// 2-color chemistry has specific error patterns
	// Should have slightly higher substitution errors between G/T and A/C
	assert.NotNil(t, data.TwoColorBias, "2-color bias should be modeled")

	// Validate that bias values are realistic
	if data.TwoColorBias != nil {
		assert.Greater(t, data.TwoColorBias.GTBias, 0.90, "G/T bias should be >90%")
		assert.Greater(t, data.TwoColorBias.ACBias, 0.90, "A/C bias should be >90%")
		assert.Less(t, data.TwoColorBias.GTBias, 0.99, "G/T bias should be <99%")
		assert.Less(t, data.TwoColorBias.ACBias, 0.99, "A/C bias should be <99%")
	}
}

func validateSBSErrorPatterns(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// SBS chemistry has characteristic error patterns
	// Most errors are substitutions, with some indels

	// Error rate should be within platform specifications
	switch platform {
	case pt.IlluminaNovaSeq6000:
		assert.Less(t, data.ErrorRate, 0.002, "NovaSeq SBS error rate")
	case pt.IlluminaNextSeq2000:
		assert.Less(t, data.ErrorRate, 0.003, "NextSeq SBS error rate")
	case pt.IlluminaMiSeq:
		assert.LessOrEqual(t, data.ErrorRate, 0.001, "MiSeq SBS error rate")
	}
}

func testIncorporationEfficiency(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// SBS incorporation efficiency affects quality scores
	// Higher efficiency = higher quality scores

	expectedMinQuality := map[pt.SequencingPlatform]float64{
		pt.IlluminaNovaSeq6000: 30.0,
		pt.IlluminaNextSeq2000: 28.0,
		pt.IlluminaMiSeq:       28.0,
	}

	if minQuality, exists := expectedMinQuality[platform]; exists {
		assert.Greater(t, data.QualityScoreAverage, minQuality,
			"Platform %s should achieve minimum quality", platform)
	}
}

func testPhasingEffects(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Phasing and prephasing cause quality decline over read length
	// Should be minimal for modern Illumina platforms

	if len(data.QualityScores) > 0 && len(data.QualityScores[0]) > 100 {
		firstQuarter := calculateAverageQuality(data.QualityScores[0][:len(data.QualityScores[0])/4])
		lastQuarter := calculateAverageQuality(data.QualityScores[0][3*len(data.QualityScores[0])/4:])

		phasingEffect := firstQuarter - lastQuarter

		// Modern platforms should have minimal phasing
		assert.Less(t, phasingEffect, 15.0, "Phasing effect should be <15 quality points")
		assert.Greater(t, phasingEffect, -2.0, "Phasing effect should be positive or minimal")

		t.Logf("%s phasing effect: %.1f quality points", platform, phasingEffect)
	}
}

func testClusterDensityEffects(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Higher cluster density can affect signal quality
	// Should be optimized for each platform

	// Quality should be consistent across sequences (no systematic cluster effects)
	if len(data.QualityScores) >= 100 {
		qualityVariances := make([]float64, 0)

		for i := 0; i < 100; i++ {
			avgQuality := calculateAverageQuality(data.QualityScores[i])
			qualityVariances = append(qualityVariances, avgQuality)
		}

		variance := calculateVariance(qualityVariances)
		assert.Less(t, variance, 25.0, "Quality variance between clusters should be low")
	}
}

func testPhasingAndPrephasing(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Test for phasing (lagging) and prephasing (leading) effects
	// Modern Illumina platforms should have <0.3% phasing and <0.1% prephasing

	expectedPhasingRates := map[pt.SequencingPlatform]float64{
		pt.IlluminaNovaSeq6000: 0.003, // 0.3%
		pt.IlluminaNextSeq2000: 0.004, // 0.4%
		pt.IlluminaMiSeq:       0.002, // 0.2%
	}

	if expectedRate, exists := expectedPhasingRates[platform]; exists {
		// Simulate phasing rate calculation from quality decline
		estimatedPhasingRate := estimatePhasingRate(data)
		assert.Less(t, estimatedPhasingRate, expectedRate*10.0,
			"Estimated phasing rate should be within expected range")
	}
}

func testSignalIntensity(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Signal intensity should be consistent and strong
	// Reflected in quality score consistency

	if len(data.QualityScores) > 0 {
		avgQuality := data.QualityScoreAverage

		// Signal intensity should be sufficient for reliable base calling
		minimumSignalQuality := map[pt.SequencingPlatform]float64{
			pt.IlluminaNovaSeq6000: 25.0,
			pt.IlluminaNextSeq2000: 23.0,
			pt.IlluminaMiSeq:       30.0,
		}

		if minQuality, exists := minimumSignalQuality[platform]; exists {
			assert.Greater(t, avgQuality, minQuality,
				"Signal intensity should be sufficient for quality >Q%.0f", minQuality)
		}
	}
}

func testChemistrySpecificBiases(t *testing.T, twoColorData, fourColorData *pt.SyntheticData) {
	// 2-color chemistry should show specific bias patterns
	assert.NotNil(t, twoColorData.TwoColorBias, "2-color data should have bias modeling")

	// 4-color chemistry should have minimal bias
	gcContent2Color := calculateGCContent(twoColorData.Sequences)
	gcContent4Color := calculateGCContent(fourColorData.Sequences)

	// Both should be close to expected GC content, but 4-color should be closer
	expectedGC := 0.42
	bias2Color := math.Abs(gcContent2Color - expectedGC)
	bias4Color := math.Abs(gcContent4Color - expectedGC)

	assert.Less(t, bias4Color, bias2Color*1.2,
		"4-color chemistry should have equal or less GC bias")
}

func testQualityScoreDistribution(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Quality scores should follow expected distribution for platform
	if len(data.QualityScores) == 0 {
		return
	}

	// Collect all quality scores
	allQualities := make([]float64, 0)
	for _, qualSeq := range data.QualityScores {
		for _, q := range qualSeq {
			allQualities = append(allQualities, float64(q-33))
		}
	}

	if len(allQualities) == 0 {
		return
	}

	// Calculate distribution statistics
	mean := calculateMean(allQualities)
	stdDev := calculateStdDev(allQualities, mean)

	// Expected ranges for different platforms
	expectedRanges := map[pt.SequencingPlatform][2]float64{
		pt.IlluminaNovaSeq6000: {30.0, 40.0},
		pt.IlluminaNextSeq2000: {28.0, 38.0},
		pt.IlluminaMiSeq:       {28.0, 42.0},
	}

	if expectedRange, exists := expectedRanges[platform]; exists {
		assert.Greater(t, mean, expectedRange[0], "Mean quality should be above minimum")
		assert.Less(t, mean, expectedRange[1], "Mean quality should be below maximum")
		assert.Less(t, stdDev, 15.0, "Quality standard deviation should be reasonable")
	}

	t.Logf("%s quality distribution: mean=%.1f, stddev=%.1f", platform, mean, stdDev)
}

func testQualityErrorCorrelation(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Quality scores should correlate with actual error rates
	// Higher quality scores should correspond to lower error rates

	avgQuality := data.QualityScoreAverage
	errorRate := data.ErrorRate

	// Convert quality to expected error rate: P_error = 10^(-Q/10)
	expectedErrorRate := math.Pow(10, -avgQuality/10)

	// Observed error rate should be reasonably close to expected
	errorRatio := errorRate / expectedErrorRate
	assert.Greater(t, errorRatio, 0.1, "Error rate shouldn't be much lower than expected")
	assert.Less(t, errorRatio, 10.0, "Error rate shouldn't be much higher than expected")

	t.Logf("%s quality-error correlation: Q%.1f expects %.1e%%, observed %.1e%% (ratio=%.2f)",
		platform, avgQuality, expectedErrorRate*100, errorRate*100, errorRatio)
}

func testPhredScoreAccuracy(t *testing.T, data *pt.SyntheticData, platform pt.SequencingPlatform) {
	// Phred calibration is not measured: there are no true errors to compare
	// against. This only checks a value derived from the configured
	// ErrorRate against arbitrary per-platform thresholds.
	expectedAccuracy := map[pt.SequencingPlatform]float64{
		pt.IlluminaNovaSeq6000: 0.95,
		pt.IlluminaNextSeq2000: 0.92,
		pt.IlluminaMiSeq:       0.97,
	}

	if accuracy, exists := expectedAccuracy[platform]; exists {
		estimatedAccuracy := 1.0 - (data.ErrorRate * 5.0)
		assert.Greater(t, estimatedAccuracy, accuracy*0.8,
			"Phred score accuracy should be reasonable for %s", platform)
	}
}

// calculateVariance must return the variance, not the standard deviation:
// its callers compare it with variance thresholds.
func TestCalculateVariance(t *testing.T) {
	values := []float64{2, 4, 4, 4, 5, 5, 7, 9} // sum of squared deviations 32
	if got, want := calculateVariance(values), 32.0/7.0; math.Abs(got-want) > 1e-12 {
		t.Errorf("calculateVariance = %v, want sample variance %v", got, want)
	}
	if got := calculateVariance([]float64{3}); got != 0 {
		t.Errorf("calculateVariance of one value = %v, want 0", got)
	}
}

// calculateVariance returns the sample variance (n-1 denominator) of values,
// or 0 for fewer than two values.
func calculateVariance(values []float64) float64 {
	if len(values) < 2 {
		return 0
	}
	mean := calculateMean(values)
	sumSquares := 0.0
	for _, v := range values {
		d := v - mean
		sumSquares += d * d
	}
	return sumSquares / float64(len(values)-1)
}

func estimatePhasingRate(data *pt.SyntheticData) float64 {
	// Simplified phasing rate estimation based on quality decline
	if len(data.QualityScores) == 0 || len(data.QualityScores[0]) < 50 {
		return 0.0
	}

	firstTenth := calculateAverageQuality(data.QualityScores[0][:len(data.QualityScores[0])/10])
	lastTenth := calculateAverageQuality(data.QualityScores[0][9*len(data.QualityScores[0])/10:])

	qualityDecline := firstTenth - lastTenth

	// Convert quality decline to estimated phasing rate (empirical formula)
	estimatedRate := qualityDecline / 1000.0
	if estimatedRate < 0 {
		estimatedRate = 0.0
	}

	return estimatedRate
}
