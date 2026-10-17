//go:build !release

package authenticity

import (
	"testing"

	genomics "github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// helpers ----------------------------------------------------------------

func newFOSWithStrategies() *FingerprintObfuscationSystem {
	fos := NewFingerprintObfuscationSystem()
	fos.InitializeObfuscationStrategies()
	fos.InitializeFingerprintAnalyzers()
	fos.InitializeCrossPlatformModels()
	return fos
}

var testSeqs = []string{"ACGTACGTACGT", "TTTTCCCCGGGG", "ACGACGACG"}
var testPlatform = genomics.NovaSeq6000
var testPop = population.CEU

// strategyWithout returns a copy of the named strategy with a specific key removed.
func strategyWithout(fos *FingerprintObfuscationSystem, name, removeKey string) *ObfuscationStrategy {
	orig := fos.obfuscationStrategies[name]
	params := make(map[string]interface{}, len(orig.ObfuscationParameters))
	for k, v := range orig.ObfuscationParameters {
		params[k] = v
	}
	delete(params, removeKey)
	return &ObfuscationStrategy{
		StrategyName:          orig.StrategyName,
		ObfuscationType:       orig.ObfuscationType,
		ObfuscationParameters: params,
	}
}

// strategyWithWrongType returns a copy with a specific key set to the wrong type.
func strategyWithWrongType(fos *FingerprintObfuscationSystem, name, key string, wrongVal interface{}) *ObfuscationStrategy {
	s := strategyWithout(fos, name, key)
	s.ObfuscationParameters[key] = wrongVal
	return s
}

// applyErrorPatternMasking -----------------------------------------------

func TestApplyErrorPatternMasking_MissingParams(t *testing.T) {
	fos := newFOSWithStrategies()

	for _, missingKey := range []string{"error_rate_normalization", "pattern_randomization", "signature_morphing"} {
		t.Run("missing_"+missingKey, func(t *testing.T) {
			s := strategyWithout(fos, "Error_Pattern_Masking", missingKey)
			_, _, _, err := fos.applyErrorPatternMasking(s, testPlatform, testPop, testSeqs)
			require.Error(t, err)
			assert.Contains(t, err.Error(), missingKey)
		})
	}
}

func TestApplyErrorPatternMasking_WrongType(t *testing.T) {
	fos := newFOSWithStrategies()
	s := strategyWithWrongType(fos, "Error_Pattern_Masking", "error_rate_normalization", "not-a-float")
	_, _, _, err := fos.applyErrorPatternMasking(s, testPlatform, testPop, testSeqs)
	require.Error(t, err)
	assert.Contains(t, err.Error(), "error_rate_normalization")
}

func TestApplyErrorPatternMasking_Valid(t *testing.T) {
	fos := newFOSWithStrategies()
	s := fos.obfuscationStrategies["Error_Pattern_Masking"]
	seqs, quals, preservation, err := fos.applyErrorPatternMasking(s, testPlatform, testPop, testSeqs)
	require.NoError(t, err)
	assert.Len(t, seqs, len(testSeqs))
	assert.Len(t, quals, len(testSeqs))
	assert.GreaterOrEqual(t, preservation, 0.0)
	assert.LessOrEqual(t, preservation, 1.0)
}

// applyQualityHomogenization ---------------------------------------------

func TestApplyQualityHomogenization_MissingParams(t *testing.T) {
	fos := newFOSWithStrategies()

	for _, missingKey := range []string{"target_mean_quality", "target_std_quality", "quality_range"} {
		t.Run("missing_"+missingKey, func(t *testing.T) {
			s := strategyWithout(fos, "Quality_Homogenization", missingKey)
			_, _, _, err := fos.applyQualityHomogenization(s, testPlatform, testPop, testSeqs)
			require.Error(t, err)
			assert.Contains(t, err.Error(), missingKey)
		})
	}
}

func TestApplyQualityHomogenization_Valid(t *testing.T) {
	fos := newFOSWithStrategies()
	s := fos.obfuscationStrategies["Quality_Homogenization"]
	seqs, quals, preservation, err := fos.applyQualityHomogenization(s, testPlatform, testPop, testSeqs)
	require.NoError(t, err)
	assert.Len(t, seqs, len(testSeqs))
	assert.Len(t, quals, len(testSeqs))
	assert.Equal(t, 1.0, preservation, "quality homogenization must preserve all sequences")
}

// applySignatureMorphing -------------------------------------------------

func TestApplySignatureMorphing_MissingParams(t *testing.T) {
	fos := newFOSWithStrategies()

	for _, missingKey := range []string{"morphing_complexity", "signature_scrambling", "cross_platform_mimicry"} {
		t.Run("missing_"+missingKey, func(t *testing.T) {
			s := strategyWithout(fos, "Advanced_Signature_Morphing", missingKey)
			_, _, _, err := fos.applySignatureMorphing(s, testPlatform, testPop, testSeqs)
			require.Error(t, err)
			assert.Contains(t, err.Error(), missingKey)
		})
	}
}

func TestApplySignatureMorphing_WrongType(t *testing.T) {
	fos := newFOSWithStrategies()
	// morphing_complexity must be int; passing float64 should fail
	s := strategyWithWrongType(fos, "Advanced_Signature_Morphing", "morphing_complexity", 5.0)
	_, _, _, err := fos.applySignatureMorphing(s, testPlatform, testPop, testSeqs)
	require.Error(t, err)
	assert.Contains(t, err.Error(), "morphing_complexity")
}

func TestApplySignatureMorphing_Valid(t *testing.T) {
	fos := newFOSWithStrategies()
	s := fos.obfuscationStrategies["Advanced_Signature_Morphing"]
	seqs, quals, _, err := fos.applySignatureMorphing(s, testPlatform, testPop, testSeqs)
	require.NoError(t, err)
	assert.Len(t, seqs, len(testSeqs))
	assert.Len(t, quals, len(testSeqs))
}

// applyPlatformMimicry ---------------------------------------------------

func TestApplyPlatformMimicry_MissingParams(t *testing.T) {
	fos := newFOSWithStrategies()

	for _, missingKey := range []string{"target_platform_rotation", "mimicry_fidelity"} {
		t.Run("missing_"+missingKey, func(t *testing.T) {
			s := strategyWithout(fos, "Cross_Platform_Mimicry", missingKey)
			_, _, _, err := fos.applyPlatformMimicry(s, testPlatform, testPop, testSeqs)
			require.Error(t, err)
			assert.Contains(t, err.Error(), missingKey)
		})
	}
}

func TestApplyPlatformMimicry_Valid(t *testing.T) {
	fos := newFOSWithStrategies()
	s := fos.obfuscationStrategies["Cross_Platform_Mimicry"]
	seqs, quals, _, err := fos.applyPlatformMimicry(s, testPlatform, testPop, testSeqs)
	require.NoError(t, err)
	assert.Len(t, seqs, len(testSeqs))
	assert.Len(t, quals, len(testSeqs))
}

// applyGCBiasNeutralization ----------------------------------------------

func TestApplyGCBiasNeutralization_MissingParams(t *testing.T) {
	fos := newFOSWithStrategies()

	for _, missingKey := range []string{"target_gc_content", "bias_strength_reduction"} {
		t.Run("missing_"+missingKey, func(t *testing.T) {
			s := strategyWithout(fos, "GC_Bias_Neutralization", missingKey)
			_, _, _, err := fos.applyGCBiasNeutralization(s, testPlatform, testPop, testSeqs)
			require.Error(t, err)
			assert.Contains(t, err.Error(), missingKey)
		})
	}
}

func TestApplyGCBiasNeutralization_Valid(t *testing.T) {
	fos := newFOSWithStrategies()
	s := fos.obfuscationStrategies["GC_Bias_Neutralization"]
	seqs, quals, _, err := fos.applyGCBiasNeutralization(s, testPlatform, testPop, testSeqs)
	require.NoError(t, err)
	assert.Len(t, seqs, len(testSeqs))
	assert.Len(t, quals, len(testSeqs))
}

// applyHomopolymerNormalization ------------------------------------------

func TestApplyHomopolymerNormalization_MissingParams(t *testing.T) {
	fos := newFOSWithStrategies()

	for _, missingKey := range []string{"max_run_length", "error_rate_normalization"} {
		t.Run("missing_"+missingKey, func(t *testing.T) {
			s := strategyWithout(fos, "Homopolymer_Normalization", missingKey)
			_, _, _, err := fos.applyHomopolymerNormalization(s, testPlatform, testPop, testSeqs)
			require.Error(t, err)
			assert.Contains(t, err.Error(), missingKey)
		})
	}
}

func TestApplyHomopolymerNormalization_WrongType(t *testing.T) {
	fos := newFOSWithStrategies()
	// max_run_length must be int; passing float64 should fail
	s := strategyWithWrongType(fos, "Homopolymer_Normalization", "max_run_length", 8.0)
	_, _, _, err := fos.applyHomopolymerNormalization(s, testPlatform, testPop, testSeqs)
	require.Error(t, err)
	assert.Contains(t, err.Error(), "max_run_length")
}

func TestApplyHomopolymerNormalization_Valid(t *testing.T) {
	fos := newFOSWithStrategies()
	s := fos.obfuscationStrategies["Homopolymer_Normalization"]
	seqs, quals, _, err := fos.applyHomopolymerNormalization(s, testPlatform, testPop, testSeqs)
	require.NoError(t, err)
	assert.Len(t, seqs, len(testSeqs))
	assert.Len(t, quals, len(testSeqs))
}

// applyObfuscationStrategy (dispatcher) ----------------------------------

func TestApplyObfuscationStrategy_UnknownType_NoError(t *testing.T) {
	fos := newFOSWithStrategies()
	s := &ObfuscationStrategy{
		ObfuscationType:       "Unknown_Type",
		ObfuscationParameters: map[string]interface{}{},
	}
	seqs, quals, preservation, err := fos.applyObfuscationStrategy(s, testPlatform, testPop, testSeqs)
	require.NoError(t, err)
	assert.Equal(t, testSeqs, seqs)
	assert.Len(t, quals, len(testSeqs))
	assert.Equal(t, 1.0, preservation)
}

func TestApplyObfuscationStrategy_PropagatesParamError(t *testing.T) {
	fos := newFOSWithStrategies()
	// Error_Pattern_Masking strategy with a missing required parameter.
	s := strategyWithout(fos, "Error_Pattern_Masking", "error_rate_normalization")
	_, _, _, err := fos.applyObfuscationStrategy(s, testPlatform, testPop, testSeqs)
	require.Error(t, err)
	assert.Contains(t, err.Error(), "error_rate_normalization")
}
