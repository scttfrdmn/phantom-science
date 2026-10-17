package quality

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewGCbiasQualityModel tests model creation for different platforms
func TestNewGCbiasQualityModel(t *testing.T) {
	testCases := []struct {
		platform string
		expectedOptimalGC float64
		expectedBiasSlope float64
		expectedSaturation string
		expectedAmplification bool
		expectedWindowSize int
	}{
		{"illumina", 0.45, 2.0, "exponential", true, 50},
		{"novaseq", 0.45, 2.0, "exponential", true, 50},
		{"nextseq", 0.45, 2.0, "exponential", true, 50},
		{"miseq", 0.45, 2.0, "exponential", true, 50},
		{"pacbio", 0.50, 0.5, "linear", false, 100},
		{"pacbio_revio", 0.50, 0.5, "linear", false, 100},
		{"nanopore", 0.50, 0.3, "linear", false, 200},
		{"promethion", 0.50, 0.3, "linear", false, 200},
		{"unknown", 0.45, 2.0, "exponential", true, 50}, // Default to Illumina
	}

	for _, tc := range testCases {
		t.Run(tc.platform, func(t *testing.T) {
			model := NewGCbiasQualityModel(tc.platform)

			assert.NotNil(t, model)
			assert.Equal(t, tc.platform, model.platform)
			assert.Equal(t, tc.expectedOptimalGC, model.optimalGCContent)
			assert.Equal(t, tc.expectedBiasSlope, model.biasSlope)
			assert.Equal(t, tc.expectedSaturation, model.saturationType)
			assert.Equal(t, tc.expectedAmplification, model.amplificationBias)
			assert.Equal(t, tc.expectedWindowSize, model.windowSize)
			assert.Equal(t, 30.0, model.baseQuality)
			assert.Equal(t, 1.0, model.biasStrength)
			assert.False(t, model.modelTrained)
			assert.NotNil(t, model.rng)
		})
	}
}

// TestGCbiasQualityModel_PredictQuality tests quality prediction
func TestGCbiasQualityModel_PredictQuality(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	model.SetSeed(12345) // For reproducible results

	t.Run("empty sequence", func(t *testing.T) {
		result, err := model.PredictQuality([]byte{}, []byte{})
		assert.Error(t, err)
		assert.Nil(t, result)
		assert.Contains(t, err.Error(), "sequence cannot be empty")
	})

	t.Run("mismatched lengths", func(t *testing.T) {
		sequence := []byte("ATCG")
		qualities := []byte{30, 30} // Wrong length

		result, err := model.PredictQuality(sequence, qualities)
		assert.Error(t, err)
		assert.Nil(t, result)
		assert.Contains(t, err.Error(), "sequence and quality lengths must match")
	})

	t.Run("normal sequence", func(t *testing.T) {
		sequence := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
		qualities := make([]byte, len(sequence))
		for i := range qualities {
			qualities[i] = 30
		}

		result, err := model.PredictQuality(sequence, qualities)
		require.NoError(t, err)
		require.NotNil(t, result)

		// Validate result structure
		assert.Equal(t, len(sequence), len(result.OriginalQualities))
		assert.Equal(t, len(sequence), len(result.AdjustedQualities))
		assert.Equal(t, len(sequence), len(result.GCContentProfile))
		assert.Equal(t, len(sequence), len(result.QualityProfile))
		assert.Equal(t, len(sequence), len(result.BiasCorrection))

		// Check that adjustments were made
		hasAdjustments := false
		for i := 0; i < len(sequence); i++ {
			if result.AdjustedQualities[i] != result.OriginalQualities[i] {
				hasAdjustments = true
				break
			}
		}
		assert.True(t, hasAdjustments, "Quality adjustments should be made")

		// Validate quality range
		for _, quality := range result.AdjustedQualities {
			assert.GreaterOrEqual(t, int(quality), 0)
			assert.LessOrEqual(t, int(quality), 41)
		}
	})

	t.Run("GC-rich sequence", func(t *testing.T) {
		// Create a GC-rich sequence
		sequence := []byte("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")
		qualities := make([]byte, len(sequence))
		for i := range qualities {
			qualities[i] = 30
		}

		result, err := model.PredictQuality(sequence, qualities)
		require.NoError(t, err)

		// GC-rich regions should have different qualities due to bias
		avgOriginal := float64(0)
		avgAdjusted := float64(0)
		for i := 0; i < len(sequence); i++ {
			avgOriginal += float64(result.OriginalQualities[i])
			avgAdjusted += float64(result.AdjustedQualities[i])
		}
		avgOriginal /= float64(len(sequence))
		avgAdjusted /= float64(len(sequence))

		// Should see some bias effect
		assert.NotEqual(t, avgOriginal, avgAdjusted, "GC bias should affect quality scores")
	})

	t.Run("AT-rich sequence", func(t *testing.T) {
		// Create an AT-rich sequence
		sequence := []byte("ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT")
		qualities := make([]byte, len(sequence))
		for i := range qualities {
			qualities[i] = 30
		}

		result, err := model.PredictQuality(sequence, qualities)
		require.NoError(t, err)

		// Check GC content profile
		for _, gcContent := range result.GCContentProfile {
			assert.LessOrEqual(t, gcContent, 0.1, "AT-rich sequence should have low GC content")
		}
	})
}

// TestGCbiasQualityModel_calculateLocalGCContent tests GC content calculation
func TestGCbiasQualityModel_calculateLocalGCContent(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")

	testCases := []struct {
		sequence []byte
		position int
		expected float64
	}{
		{[]byte("GGCCGGCCGGCCGGCC"), 8, 1.0}, // Pure GC
		{[]byte("ATATATATATATATAT"), 8, 0.0}, // Pure AT
		{[]byte("ATGCATGCATGCATGC"), 8, 0.5}, // Balanced
		{[]byte("ATCG"), 1, 0.5}, // Short sequence
		{[]byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), 30, 0.0}, // Long AT-rich
	}

	for _, tc := range testCases {
		t.Run(string(tc.sequence), func(t *testing.T) {
			gcContent := model.calculateLocalGCContent(tc.sequence, tc.position)
			assert.InDelta(t, tc.expected, gcContent, 0.1, "GC content should be close to expected")
		})
	}
}

// TestGCbiasQualityModel_calculateGCBiasAdjustment tests bias adjustment calculation
func TestGCbiasQualityModel_calculateGCBiasAdjustment(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	model.SetSeed(12345)

	testCases := []struct {
		gcContent float64
		description string
	}{
		{0.45, "optimal GC content"}, // Should have minimal bias
		{0.65, "high GC content"},    // Should have negative bias
		{0.25, "low GC content"},     // Should have negative bias
		{0.0, "no GC content"},       // Should have strong negative bias
		{1.0, "pure GC content"},     // Should have strong negative bias
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			adjustment := model.calculateGCBiasAdjustment(tc.gcContent)

			// Adjustment should be a reasonable value
			assert.GreaterOrEqual(t, adjustment, -20.0, "Adjustment should not be extremely negative")
			assert.LessOrEqual(t, adjustment, 5.0, "Adjustment should not be extremely positive")

			// Optimal GC content should have minimal bias
			if tc.gcContent == model.optimalGCContent {
				assert.LessOrEqual(t, math.Abs(adjustment), 2.0, "Optimal GC should have minimal bias")
			}
		})
	}
}

// TestGCbiasQualityModel_calculatePlatformEffects tests platform-specific effects
func TestGCbiasQualityModel_calculatePlatformEffects(t *testing.T) {
	testCases := []struct {
		platform string
		sequence []byte
		position int
	}{
		{"illumina", []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"), 25},
		{"pacbio", []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"), 25},
		{"nanopore", []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"), 25},
	}

	for _, tc := range testCases {
		t.Run(tc.platform, func(t *testing.T) {
			model := NewGCbiasQualityModel(tc.platform)
			model.SetSeed(12345)

			effect := model.calculatePlatformEffects(tc.sequence, tc.position)

			// Effect should be reasonable
			assert.GreaterOrEqual(t, effect, -10.0, "Platform effect should not be extremely negative")
			assert.LessOrEqual(t, effect, 5.0, "Platform effect should not be extremely positive")
		})
	}
}

// TestGCbiasQualityModel_TrainFromRealData tests model training
func TestGCbiasQualityModel_TrainFromRealData(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")

	t.Run("empty training data", func(t *testing.T) {
		err := model.TrainFromRealData(nil)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "no training data provided")
	})

	t.Run("valid training data", func(t *testing.T) {
		// Create synthetic training data
		trainingData := []GCQualityDataPoint{
			{GCContent: 0.2, QualityScore: 25.0, Platform: "illumina"},
			{GCContent: 0.3, QualityScore: 28.0, Platform: "illumina"},
			{GCContent: 0.4, QualityScore: 32.0, Platform: "illumina"},
			{GCContent: 0.45, QualityScore: 35.0, Platform: "illumina"}, // Optimal
			{GCContent: 0.5, QualityScore: 33.0, Platform: "illumina"},
			{GCContent: 0.6, QualityScore: 28.0, Platform: "illumina"},
			{GCContent: 0.7, QualityScore: 22.0, Platform: "illumina"},
			{GCContent: 0.8, QualityScore: 18.0, Platform: "illumina"},
		}

		err := model.TrainFromRealData(trainingData)
		require.NoError(t, err)

		assert.True(t, model.modelTrained)
		assert.Len(t, model.trainingData, len(trainingData))
		assert.Greater(t, len(model.qualityAdjustment), 0)

		// Model should identify optimal GC content
		assert.InDelta(t, 0.45, model.optimalGCContent, 0.05, "Should identify optimal GC content")
	})
}

// TestGCbiasQualityModel_ValidateModel tests model validation
func TestGCbiasQualityModel_ValidateModel(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")

	t.Run("untrained model", func(t *testing.T) {
		err := model.ValidateModel()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "model has not been trained")
	})

	t.Run("trained model", func(t *testing.T) {
		// Train with minimal data
		trainingData := []GCQualityDataPoint{
			{GCContent: 0.45, QualityScore: 30.0, Platform: "illumina"},
		}

		err := model.TrainFromRealData(trainingData)
		require.NoError(t, err)

		err = model.ValidateModel()
		assert.NoError(t, err)
	})
}

// TestGCbiasQualityModel_SettersAndGetters tests configuration methods
func TestGCbiasQualityModel_SettersAndGetters(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")

	t.Run("SetBiasStrength", func(t *testing.T) {
		model.SetBiasStrength(1.5)
		assert.Equal(t, 1.5, model.biasStrength)

		// Negative values should be ignored
		model.SetBiasStrength(-0.5)
		assert.Equal(t, 1.5, model.biasStrength) // Should remain unchanged
	})

	t.Run("SetWindowSize", func(t *testing.T) {
		model.SetWindowSize(75)
		assert.Equal(t, 75, model.windowSize)

		// Non-positive values should be ignored
		model.SetWindowSize(0)
		assert.Equal(t, 75, model.windowSize) // Should remain unchanged
	})

	t.Run("GetModelStats", func(t *testing.T) {
		stats := model.GetModelStats()

		assert.Equal(t, "illumina", stats["platform"])
		assert.Equal(t, false, stats["model_trained"])
		assert.Equal(t, 0.45, stats["optimal_gc_content"])
		assert.Equal(t, 2.0, stats["bias_slope"])
		assert.Equal(t, "exponential", stats["saturation_type"])
		assert.Equal(t, 1.5, stats["bias_strength"]) // From previous test
		assert.Equal(t, 75, stats["window_size"]) // From previous test
		assert.Equal(t, true, stats["amplification_bias"])
		assert.Equal(t, 0, stats["training_data_points"])
		assert.Equal(t, 0, stats["gc_bins"])
	})

	t.Run("GetPlatformParameters", func(t *testing.T) {
		params := model.GetPlatformParameters()

		assert.Equal(t, "illumina", params["platform"])
		assert.Equal(t, 0.45, params["optimal_gc_content"])
		assert.Equal(t, 2.0, params["bias_slope"])
		assert.Equal(t, "exponential", params["saturation_type"])
		assert.Equal(t, true, params["amplification_bias"])
		assert.Equal(t, 75, params["window_size"])
	})
}

// TestGCbiasQualityModel_SimulateQualityProfile tests quality simulation
func TestGCbiasQualityModel_SimulateQualityProfile(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	model.SetSeed(12345)

	t.Run("empty sequence", func(t *testing.T) {
		qualities, err := model.SimulateQualityProfile([]byte{}, 30.0)
		assert.Error(t, err)
		assert.Nil(t, qualities)
	})

	t.Run("normal sequence", func(t *testing.T) {
		sequence := []byte("ATCGATCGATCGATCGATCGATCGATCGATCG")
		qualities, err := model.SimulateQualityProfile(sequence, 30.0)
		require.NoError(t, err)

		assert.Equal(t, len(sequence), len(qualities))

		// Qualities should be in valid range
		for _, quality := range qualities {
			assert.GreaterOrEqual(t, int(quality), 0)
			assert.LessOrEqual(t, int(quality), 41)
		}
	})
}

// TestGCbiasQualityModel_calculateSequenceComplexity tests complexity calculation
func TestGCbiasQualityModel_calculateSequenceComplexity(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")

	testCases := []struct {
		sequence []byte
		expected float64
		description string
	}{
		{[]byte{}, 0.0, "empty sequence"},
		{[]byte("AAA"), 0.0, "too short sequence"},
		{[]byte("AAAA"), 0.5, "low complexity"}, // Only "AAA" k-mer
		{[]byte("ATCG"), 1.0, "high complexity"}, // "ATC", "TCG" k-mers
		{[]byte("ATCGATCG"), 1.0, "varied sequence"}, // Multiple unique k-mers
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			complexity := model.calculateSequenceComplexity(tc.sequence)
			if len(tc.sequence) >= 3 {
				assert.GreaterOrEqual(t, complexity, 0.0)
				assert.LessOrEqual(t, complexity, 1.0)
			} else {
				assert.Equal(t, tc.expected, complexity)
			}
		})
	}
}

// TestGCbiasQualityModel_calculateEntropy tests entropy calculation
func TestGCbiasQualityModel_calculateEntropy(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")

	testCases := []struct {
		sequence []byte
		minEntropy float64
		maxEntropy float64
		description string
	}{
		{[]byte{}, 0.0, 0.0, "empty sequence"},
		{[]byte("AAAA"), 0.0, 0.0, "no entropy"},
		{[]byte("ATCG"), 1.8, 2.0, "high entropy"},
		{[]byte("ATAT"), 0.9, 1.1, "medium entropy"},
		{[]byte("AAAAAAAAT"), 0.4, 0.6, "low entropy"},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			entropy := model.calculateEntropy(tc.sequence)
			assert.GreaterOrEqual(t, entropy, tc.minEntropy)
			assert.LessOrEqual(t, entropy, tc.maxEntropy)
		})
	}
}

// TestGCQualityDataPoint_Validation tests data point structure
func TestGCQualityDataPoint_Validation(t *testing.T) {
	dataPoint := GCQualityDataPoint{
		GCContent:      0.45,
		QualityScore:   32.5,
		SequenceLength: 150,
		Platform:       "illumina",
		Cycle:          75,
		ReadPosition:   50,
	}

	assert.GreaterOrEqual(t, dataPoint.GCContent, 0.0)
	assert.LessOrEqual(t, dataPoint.GCContent, 1.0)
	assert.Greater(t, dataPoint.QualityScore, 0.0)
	assert.Greater(t, dataPoint.SequenceLength, 0)
	assert.NotEmpty(t, dataPoint.Platform)
	assert.GreaterOrEqual(t, dataPoint.Cycle, 0)
	assert.GreaterOrEqual(t, dataPoint.ReadPosition, 0)
}

// TestQualityPredictionResult_Validation tests result structure
func TestQualityPredictionResult_Validation(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	model.SetSeed(12345)

	sequence := []byte("ATCGATCGATCGATCGATCGATCGATCG")
	qualities := make([]byte, len(sequence))
	for i := range qualities {
		qualities[i] = 30
	}

	result, err := model.PredictQuality(sequence, qualities)
	require.NoError(t, err)

	// Validate all arrays have correct length
	assert.Equal(t, len(sequence), len(result.OriginalQualities))
	assert.Equal(t, len(sequence), len(result.AdjustedQualities))
	assert.Equal(t, len(sequence), len(result.GCContentProfile))
	assert.Equal(t, len(sequence), len(result.QualityProfile))
	assert.Equal(t, len(sequence), len(result.BiasCorrection))

	// Validate GC content values
	for _, gcContent := range result.GCContentProfile {
		assert.GreaterOrEqual(t, gcContent, 0.0)
		assert.LessOrEqual(t, gcContent, 1.0)
	}

	// Validate quality values
	for _, quality := range result.AdjustedQualities {
		assert.GreaterOrEqual(t, int(quality), 0)
		assert.LessOrEqual(t, int(quality), 41)
	}

	// Validate platform effects
	for _, effect := range result.PlatformEffects {
		assert.GreaterOrEqual(t, effect.Position, 0)
		assert.Less(t, effect.Position, len(sequence))
		assert.NotEmpty(t, effect.EffectType)
		assert.NotEmpty(t, effect.Description)
		assert.GreaterOrEqual(t, math.Abs(effect.Magnitude), 2.0) // Should be significant
	}
}

// TestPlatformEffect_Validation tests platform effect structure
func TestPlatformEffect_Validation(t *testing.T) {
	effect := PlatformEffect{
		Position:    150,
		EffectType:  "gc_bias",
		Magnitude:   -2.5,
		Description: "High GC content reduces polymerase efficiency",
	}

	assert.GreaterOrEqual(t, effect.Position, 0)
	assert.Contains(t, []string{"gc_bias", "amplification", "priming", "cluster_density"}, effect.EffectType)
	assert.NotEqual(t, 0.0, effect.Magnitude)
	assert.NotEmpty(t, effect.Description)
}

func TestGetEffectType_HighGC(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	// All G/C sequence = high GC -> "gc_bias"
	seq := []byte("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")
	got := model.getEffectType(seq, 15)
	if got != "gc_bias" {
		t.Errorf("getEffectType(high GC) = %q, want gc_bias", got)
	}
}

func TestGetEffectType_LowGC(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	// All A/T sequence = low GC -> "gc_bias"
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
	got := model.getEffectType(seq, 15)
	if got != "gc_bias" {
		t.Errorf("getEffectType(low GC) = %q, want gc_bias", got)
	}
}

func TestGetEffectType_PrimingContext(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	// Include AAAA at position-6 to trigger "priming"
	seq := []byte("AAAAATCGATCGATCGATCGATCGATCGAT")
	// Middle GC content should be fine; hexamer check at position>=6
	got := model.getEffectType(seq, 6)
	// Either "priming" (if hexamer triggers) or "cluster_density" (fallback)
	if got == "" {
		t.Errorf("getEffectType(priming context) returned empty string")
	}
}

func TestGetEffectDescription_HighGC(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	seq := []byte("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")
	got := model.getEffectDescription(seq, 15)
	if got == "" {
		t.Error("getEffectDescription(high GC) returned empty string")
	}
}

func TestGetEffectDescription_LowGC(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
	got := model.getEffectDescription(seq, 15)
	if got == "" {
		t.Error("getEffectDescription(low GC) returned empty string")
	}
}

func TestGetEffectDescription_Normal(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	seq := []byte("ATGCATGCATGCATGCATGCATGCATGCAT")
	got := model.getEffectDescription(seq, 15)
	if got == "" {
		t.Error("getEffectDescription(normal GC) returned empty string")
	}
}

func TestStringIndex_Found(t *testing.T) {
	got := stringIndex("ATCGATCG", "TCG")
	if got != 1 {
		t.Errorf("stringIndex found = %d, want 1", got)
	}
}

func TestStringIndex_NotFound(t *testing.T) {
	got := stringIndex("ATCG", "NNN")
	if got != -1 {
		t.Errorf("stringIndex not found = %d, want -1", got)
	}
}

func TestGetModelStats_NotNil(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	stats := model.GetModelStats()
	if stats == nil {
		t.Error("GetModelStats returned nil")
	}
	if _, ok := stats["platform"]; !ok {
		t.Error("GetModelStats missing 'platform' key")
	}
}

func TestGetModelStats_WithTraining(t *testing.T) {
	model := NewGCbiasQualityModel("illumina")
	dataPoints := []GCQualityDataPoint{
		{GCContent: 0.42, QualityScore: 35, ReadPosition: 1},
		{GCContent: 0.50, QualityScore: 32, ReadPosition: 2},
	}
	_ = model.TrainFromRealData(dataPoints)
	stats := model.GetModelStats()
	if stats["training_data_points"].(int) != 2 {
		t.Errorf("GetModelStats training_data_points = %v, want 2", stats["training_data_points"])
	}
}

// BenchmarkGCbiasQualityModel_PredictQuality benchmarks quality prediction
func BenchmarkGCbiasQualityModel_PredictQuality(b *testing.B) {
	model := NewGCbiasQualityModel("illumina")
	model.SetSeed(12345)

	// Create a realistic sequence
	sequence := make([]byte, 1000)
	for i := 0; i < len(sequence); i += 4 {
		if i+3 < len(sequence) {
			sequence[i] = 'A'
			sequence[i+1] = 'T'
			sequence[i+2] = 'C'
			sequence[i+3] = 'G'
		}
	}

	qualities := make([]byte, len(sequence))
	for i := range qualities {
		qualities[i] = 30
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = model.PredictQuality(sequence, qualities)
	}
}

// BenchmarkGCbiasQualityModel_calculateLocalGCContent benchmarks GC calculation
func BenchmarkGCbiasQualityModel_calculateLocalGCContent(b *testing.B) {
	model := NewGCbiasQualityModel("illumina")

	sequence := make([]byte, 10000)
	for i := 0; i < len(sequence); i += 4 {
		if i+3 < len(sequence) {
			sequence[i] = 'A'
			sequence[i+1] = 'T'
			sequence[i+2] = 'C'
			sequence[i+3] = 'G'
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for pos := 0; pos < len(sequence); pos += 100 {
			model.calculateLocalGCContent(sequence, pos)
		}
	}
}