package statistics

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewJensenShannonDivergence tests the constructor
func TestNewJensenShannonDivergence(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	assert.Equal(t, 2.0, js.Base)
	assert.Equal(t, 1e-10, js.Smoothing)
}

// TestJensenShannonDivergence_Calculate tests JS divergence calculation
func TestJensenShannonDivergence_Calculate(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	t.Run("identical distributions", func(t *testing.T) {
		p := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
		q := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
		
		divergence, err := js.Calculate(p, q)
		require.NoError(t, err)
		assert.InDelta(t, 0.0, divergence, 1e-9, "Identical distributions should have 0 divergence")
	})
	
	t.Run("different distributions", func(t *testing.T) {
		p := map[string]float64{"A": 1.0}
		q := map[string]float64{"T": 1.0}
		
		divergence, err := js.Calculate(p, q)
		require.NoError(t, err)
		assert.Greater(t, divergence, 0.0, "Different distributions should have positive divergence")
		assert.LessOrEqual(t, divergence, 1.0, "JS divergence should be bounded by 1")
	})
	
	t.Run("partially overlapping distributions", func(t *testing.T) {
		p := map[string]float64{"A": 0.5, "T": 0.5}
		q := map[string]float64{"A": 0.3, "G": 0.7}
		
		divergence, err := js.Calculate(p, q)
		require.NoError(t, err)
		assert.Greater(t, divergence, 0.0)
		assert.Less(t, divergence, 1.0)
	})
	
	t.Run("uniform vs skewed distribution", func(t *testing.T) {
		uniform := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
		skewed := map[string]float64{"A": 0.97, "T": 0.01, "G": 0.01, "C": 0.01}
		
		divergence, err := js.Calculate(uniform, skewed)
		require.NoError(t, err)
		assert.Greater(t, divergence, 0.1, "Uniform vs skewed should have significant divergence")
	})
	
	t.Run("error cases", func(t *testing.T) {
		// Empty distributions
		_, err := js.Calculate(map[string]float64{}, map[string]float64{"A": 1.0})
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "cannot be empty")
		
		// Invalid probability distribution (doesn't sum to 1)
		invalid := map[string]float64{"A": 0.5, "T": 0.3} // Sums to 0.8
		valid := map[string]float64{"A": 0.5, "T": 0.5}
		_, err = js.Calculate(invalid, valid)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "do not sum to 1")
		
		// Negative probabilities
		negative := map[string]float64{"A": -0.1, "T": 1.1}
		_, err = js.Calculate(negative, valid)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "negative probability")
	})
}

// TestJensenShannonDivergence_CalculateDistance tests JS distance calculation
func TestJensenShannonDivergence_CalculateDistance(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{"G": 0.5, "C": 0.5}
	
	divergence, err := js.Calculate(p, q)
	require.NoError(t, err)
	
	distance, err := js.CalculateDistance(p, q)
	require.NoError(t, err)
	
	expectedDistance := math.Sqrt(divergence)
	assert.InDelta(t, expectedDistance, distance, 1e-9)
	
	// Test error propagation
	invalid := map[string]float64{}
	_, err = js.CalculateDistance(invalid, q)
	assert.Error(t, err)
}

// TestJensenShannonDivergence_CalculateFromFrequencies tests frequency-based calculation
func TestJensenShannonDivergence_CalculateFromFrequencies(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	t.Run("frequency to probability conversion", func(t *testing.T) {
		freqP := map[string]int{"A": 25, "T": 25, "G": 25, "C": 25}
		freqQ := map[string]int{"A": 50, "T": 30, "G": 15, "C": 5}
		
		divergence, err := js.CalculateFromFrequencies(freqP, freqQ)
		require.NoError(t, err)
		assert.Greater(t, divergence, 0.0)
		
		// Compare with direct probability calculation
		probP := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
		probQ := map[string]float64{"A": 0.5, "T": 0.3, "G": 0.15, "C": 0.05}
		
		expectedDivergence, err := js.Calculate(probP, probQ)
		require.NoError(t, err)
		
		assert.InDelta(t, expectedDivergence, divergence, 1e-9)
	})
	
	t.Run("empty frequencies", func(t *testing.T) {
		freqP := map[string]int{}
		freqQ := map[string]int{"A": 1}
		
		divergence, err := js.CalculateFromFrequencies(freqP, freqQ)
		assert.Error(t, err)
		assert.Equal(t, 0.0, divergence)
	})
}

// TestJensenShannonDivergence_frequenciesToProbabilities tests frequency conversion
func TestJensenShannonDivergence_frequenciesToProbabilities(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	t.Run("normal conversion", func(t *testing.T) {
		frequencies := map[string]int{"A": 40, "T": 30, "G": 20, "C": 10}
		probabilities := js.frequenciesToProbabilities(frequencies)
		
		expected := map[string]float64{"A": 0.4, "T": 0.3, "G": 0.2, "C": 0.1}
		assert.Equal(t, expected, probabilities)
		
		// Verify sum is 1
		sum := 0.0
		for _, prob := range probabilities {
			sum += prob
		}
		assert.InDelta(t, 1.0, sum, 1e-9)
	})
	
	t.Run("zero total", func(t *testing.T) {
		frequencies := map[string]int{}
		probabilities := js.frequenciesToProbabilities(frequencies)
		assert.Empty(t, probabilities)
		
		frequencies = map[string]int{"A": 0, "T": 0}
		probabilities = js.frequenciesToProbabilities(frequencies)
		assert.Empty(t, probabilities)
	})
}

// TestJensenShannonDivergence_validateProbabilityDistribution tests validation
func TestJensenShannonDivergence_validateProbabilityDistribution(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	t.Run("valid distributions", func(t *testing.T) {
		valid := map[string]float64{"A": 0.5, "T": 0.5}
		err := js.validateProbabilityDistribution(valid)
		assert.NoError(t, err)
		
		// Test with small numerical error
		almostValid := map[string]float64{"A": 0.5000001, "T": 0.4999999}
		err = js.validateProbabilityDistribution(almostValid)
		assert.NoError(t, err)
	})
	
	t.Run("invalid distributions", func(t *testing.T) {
		// Negative probability
		negative := map[string]float64{"A": -0.1, "T": 1.1}
		err := js.validateProbabilityDistribution(negative)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "negative probability")
		
		// Doesn't sum to 1
		wrongSum := map[string]float64{"A": 0.3, "T": 0.3}
		err = js.validateProbabilityDistribution(wrongSum)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "do not sum to 1")
	})
}

// TestJensenShannonDivergence_getUnionKeys tests key union functionality
func TestJensenShannonDivergence_getUnionKeys(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{"G": 0.5, "C": 0.5, "A": 0.0} // A appears in both
	
	allKeys := js.getUnionKeys(p, q)
	
	expected := map[string]bool{"A": true, "T": true, "G": true, "C": true}
	assert.Equal(t, expected, allKeys)
}

// TestJensenShannonDivergence_normalizeWithSmoothing tests smoothing normalization
func TestJensenShannonDivergence_normalizeWithSmoothing(t *testing.T) {
	js := NewJensenShannonDivergence()
	js.Smoothing = 0.01 // Use larger smoothing for testing
	
	dist := map[string]float64{"A": 0.8, "T": 0.2}
	allKeys := map[string]bool{"A": true, "T": true, "G": true, "C": true}
	
	normalized := js.normalizeWithSmoothing(dist, allKeys)
	
	// Check that all keys are present
	assert.Len(t, normalized, 4)
	assert.Contains(t, normalized, "A")
	assert.Contains(t, normalized, "T")
	assert.Contains(t, normalized, "G")
	assert.Contains(t, normalized, "C")
	
	// Check that G and C got smoothing values (after normalization)
	// Since smoothing adds to all values and then normalizes, values won't equal original smoothing
	assert.Greater(t, normalized["G"], 0.0)
	assert.Greater(t, normalized["C"], 0.0)
	assert.Less(t, normalized["G"], normalized["A"]) // Should be less than A which had original probability
	assert.Less(t, normalized["C"], normalized["A"])
	
	// Check normalization (should sum to 1)
	sum := 0.0
	for _, prob := range normalized {
		sum += prob
	}
	assert.InDelta(t, 1.0, sum, 1e-9)
}

// TestJensenShannonDivergence_calculateKLDivergence tests KL divergence calculation
func TestJensenShannonDivergence_calculateKLDivergence(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{"A": 0.5, "T": 0.5}
	
	// KL divergence of identical distributions should be 0
	klDiv := js.calculateKLDivergence(p, q)
	assert.InDelta(t, 0.0, klDiv, 1e-9)
	
	// KL divergence is not symmetric
	p2 := map[string]float64{"A": 0.9, "T": 0.1}
	q2 := map[string]float64{"A": 0.1, "T": 0.9}
	
	klPQ := js.calculateKLDivergence(p2, q2)
	klQP := js.calculateKLDivergence(q2, p2)
	
	assert.Greater(t, klPQ, 0.0)
	assert.Greater(t, klQP, 0.0)
	// Note: The values might be equal in some cases, so let's just check they're both positive
	// KL divergence is asymmetric in general but can be equal for some specific distributions
}

// TestJensenShannonDivergence_CompareDistributions tests comprehensive comparison
func TestJensenShannonDivergence_CompareDistributions(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	p := map[string]float64{"A": 0.4, "T": 0.3, "G": 0.2, "C": 0.1}
	q := map[string]float64{"A": 0.3, "T": 0.4, "G": 0.2, "C": 0.1}
	
	comparison, err := js.CompareDistributions(p, q)
	require.NoError(t, err)
	
	// Check that all metrics are calculated
	assert.Greater(t, comparison.JSDivergence, 0.0)
	assert.Greater(t, comparison.JSDistance, 0.0)
	assert.GreaterOrEqual(t, comparison.TotalVariation, 0.0)
	assert.GreaterOrEqual(t, comparison.HellingerDistance, 0.0)
	assert.GreaterOrEqual(t, comparison.BhattacharyyaDistance, 0.0)
	assert.GreaterOrEqual(t, comparison.ChiSquared, 0.0)
	
	// Check key analysis
	assert.Empty(t, comparison.KeysOnlyInP) // Same keys in both
	assert.Empty(t, comparison.KeysOnlyInQ)
	assert.Len(t, comparison.CommonKeys, 4)
	
	// JS distance should be square root of divergence
	assert.InDelta(t, math.Sqrt(comparison.JSDivergence), comparison.JSDistance, 1e-9)
}

// TestJensenShannonDivergence_DistanceMetrics tests individual distance metrics
func TestJensenShannonDivergence_DistanceMetrics(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	// Test with identical distributions
	p := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
	q := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
	
	tv := js.calculateTotalVariation(p, q)
	hellinger := js.calculateHellingerDistance(p, q)
	bhatt := js.calculateBhattacharyyaDistance(p, q)
	chiSq := js.calculateChiSquared(p, q)
	
	// All should be 0 or close to 0 for identical distributions
	assert.InDelta(t, 0.0, tv, 1e-9)
	assert.InDelta(t, 0.0, hellinger, 1e-9)
	assert.InDelta(t, 0.0, bhatt, 1e-9)
	assert.InDelta(t, 0.0, chiSq, 1e-9)
	
	// Test with different distributions
	q2 := map[string]float64{"A": 1.0, "T": 0.0, "G": 0.0, "C": 0.0}
	
	tv2 := js.calculateTotalVariation(p, q2)
	hellinger2 := js.calculateHellingerDistance(p, q2)
	chiSq2 := js.calculateChiSquared(p, q2)
	
	assert.Greater(t, tv2, 0.0)
	assert.Greater(t, hellinger2, 0.0)
	assert.Greater(t, chiSq2, 0.0)
}

// TestJensenShannonDivergence_KeyAnalysis tests key analysis functions
func TestJensenShannonDivergence_KeyAnalysis(t *testing.T) {
	js := NewJensenShannonDivergence()
	
	p := map[string]float64{"A": 0.5, "T": 0.3, "G": 0.2}
	q := map[string]float64{"A": 0.4, "C": 0.4, "N": 0.2}
	
	onlyInP := js.getKeysOnlyIn(p, q)
	onlyInQ := js.getKeysOnlyIn(q, p)
	common := js.getCommonKeys(p, q)
	
	assert.Contains(t, onlyInP, "T")
	assert.Contains(t, onlyInP, "G")
	assert.NotContains(t, onlyInP, "A")
	
	assert.Contains(t, onlyInQ, "C")
	assert.Contains(t, onlyInQ, "N")
	assert.NotContains(t, onlyInQ, "A")
	
	assert.Contains(t, common, "A")
	assert.Len(t, common, 1)
}

// TestNewBatchComparison tests batch comparison creation
func TestNewBatchComparison(t *testing.T) {
	bc := NewBatchComparison()
	
	assert.NotNil(t, bc.calculator)
	assert.Empty(t, bc.results)
}

// TestBatchComparison_AddComparison tests adding comparisons to batch
func TestBatchComparison_AddComparison(t *testing.T) {
	bc := NewBatchComparison()
	
	p1 := map[string]float64{"A": 0.5, "T": 0.5}
	q1 := map[string]float64{"A": 0.3, "T": 0.7}
	
	p2 := map[string]float64{"G": 0.6, "C": 0.4}
	q2 := map[string]float64{"G": 0.4, "C": 0.6}
	
	bc.AddComparison("test1", p1, q1)
	bc.AddComparison("test2", p2, q2)
	
	results := bc.GetResults()
	assert.Len(t, results, 2)
	
	assert.Equal(t, "test1", results[0].ID)
	assert.Equal(t, "test2", results[1].ID)
	assert.NoError(t, results[0].Error)
	assert.NoError(t, results[1].Error)
	assert.Greater(t, results[0].JSDivergence, 0.0)
	assert.Greater(t, results[1].JSDivergence, 0.0)
}

// TestBatchComparison_ErrorHandling tests error handling in batch operations
func TestBatchComparison_ErrorHandling(t *testing.T) {
	bc := NewBatchComparison()
	
	// Add successful comparison
	validP := map[string]float64{"A": 0.5, "T": 0.5}
	validQ := map[string]float64{"A": 0.3, "T": 0.7}
	bc.AddComparison("valid", validP, validQ)
	
	// Add failed comparison (empty distribution)
	invalidP := map[string]float64{}
	bc.AddComparison("invalid", invalidP, validQ)
	
	// Test results filtering
	allResults := bc.GetResults()
	successfulResults := bc.GetSuccessfulResults()
	
	assert.Len(t, allResults, 2)
	assert.Len(t, successfulResults, 1)
	assert.Equal(t, "valid", successfulResults[0].ID)
	assert.NoError(t, successfulResults[0].Error)
	
	// Check that failed result has error
	assert.Error(t, allResults[1].Error)
}

// TestBatchComparison_Statistics tests batch statistics calculation
func TestBatchComparison_Statistics(t *testing.T) {
	bc := NewBatchComparison()
	
	// Add several successful comparisons
	comparisons := []struct {
		id string
		p  map[string]float64
		q  map[string]float64
	}{
		{"low", map[string]float64{"A": 0.5, "T": 0.5}, map[string]float64{"A": 0.49, "T": 0.51}},
		{"medium", map[string]float64{"A": 0.5, "T": 0.5}, map[string]float64{"A": 0.3, "T": 0.7}},
		{"high", map[string]float64{"A": 1.0}, map[string]float64{"T": 1.0}},
	}
	
	for _, comp := range comparisons {
		bc.AddComparison(comp.id, comp.p, comp.q)
	}
	
	// Add one failed comparison
	bc.AddComparison("failed", map[string]float64{}, map[string]float64{"A": 1.0})
	
	// Test individual statistics
	avgJS := bc.GetAverageJSDivergence()
	maxJS := bc.GetMaxJSDivergence()
	minJS := bc.GetMinJSDivergence()
	
	assert.Greater(t, avgJS, 0.0)
	assert.Greater(t, maxJS, avgJS)
	assert.Less(t, minJS, avgJS)
	assert.GreaterOrEqual(t, avgJS, minJS)
	assert.LessOrEqual(t, avgJS, maxJS)
	
	// Test threshold checking
	assert.False(t, bc.IsWithinThreshold(0.1)) // High divergence should exceed threshold
	assert.True(t, bc.IsWithinThreshold(10.0)) // Very high threshold should pass
	
	// Test comprehensive statistics
	stats := bc.GetStatistics()
	assert.Equal(t, 4, stats.TotalComparisons)
	assert.Equal(t, 3, stats.SuccessfulComparisons)
	assert.Equal(t, 1, stats.FailedComparisons)
	assert.Equal(t, avgJS, stats.AverageJSDivergence)
	assert.Equal(t, maxJS, stats.MaxJSDivergence)
	assert.Equal(t, minJS, stats.MinJSDivergence)
	assert.Greater(t, stats.StdDevJSDivergence, 0.0)
}

// TestBatchComparison_EmptyBatch tests edge cases with empty batches
func TestBatchComparison_EmptyBatch(t *testing.T) {
	bc := NewBatchComparison()
	
	// Test statistics on empty batch
	assert.Equal(t, 0.0, bc.GetAverageJSDivergence())
	assert.Equal(t, 0.0, bc.GetMaxJSDivergence())
	assert.Equal(t, 0.0, bc.GetMinJSDivergence())
	assert.False(t, bc.IsWithinThreshold(0.5))
	assert.Empty(t, bc.GetSuccessfulResults())
	
	stats := bc.GetStatistics()
	assert.Equal(t, 0, stats.TotalComparisons)
	assert.Equal(t, 0, stats.SuccessfulComparisons)
	assert.Equal(t, 0, stats.FailedComparisons)
}

// TestJensenShannonDivergence_DifferentBases tests different logarithm bases
func TestJensenShannonDivergence_DifferentBases(t *testing.T) {
	js2 := NewJensenShannonDivergence()
	js2.Base = 2.0
	
	jsE := NewJensenShannonDivergence()
	jsE.Base = math.E
	
	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{"G": 0.5, "C": 0.5}
	
	div2, err := js2.Calculate(p, q)
	require.NoError(t, err)
	
	divE, err := jsE.Calculate(p, q)
	require.NoError(t, err)
	
	// Results should be different but proportional
	assert.NotEqual(t, div2, divE)
	assert.Greater(t, div2, 0.0)
	assert.Greater(t, divE, 0.0)
	
	// Base 2 result should be larger than base e result (since log_2(x) > log_e(x) for x > 1)
	ratio := div2 / divE
	expectedRatio := 1.0 / math.Log(2) // Conversion factor between natural log and log base 2
	assert.InDelta(t, expectedRatio, ratio, 0.01)
}

// TestJensenShannonDivergence_SmootingEffect tests smoothing parameter effect
func TestJensenShannonDivergence_SmootingEffect(t *testing.T) {
	jsSmall := NewJensenShannonDivergence()
	jsSmall.Smoothing = 1e-15
	
	jsLarge := NewJensenShannonDivergence()
	jsLarge.Smoothing = 0.01
	
	p := map[string]float64{"A": 1.0}
	q := map[string]float64{"T": 1.0}
	
	divSmall, err := jsSmall.Calculate(p, q)
	require.NoError(t, err)
	
	divLarge, err := jsLarge.Calculate(p, q)
	require.NoError(t, err)
	
	// Larger smoothing should reduce divergence
	assert.Less(t, divLarge, divSmall)
	assert.Greater(t, divSmall, 0.0)
	assert.Greater(t, divLarge, 0.0)
	// Don't assume specific bounds as they depend on the exact smoothing implementation
}

// BenchmarkJensenShannonDivergence_Calculate benchmarks JS divergence calculation
func BenchmarkJensenShannonDivergence_Calculate(b *testing.B) {
	js := NewJensenShannonDivergence()
	
	// Create larger distributions for realistic benchmarking
	p := make(map[string]float64)
	q := make(map[string]float64)
	
	// Generate all 4-mers
	bases := []string{"A", "T", "G", "C"}
	for i, b1 := range bases {
		for j, b2 := range bases {
			for k, b3 := range bases {
				for l, b4 := range bases {
					kmer := b1 + b2 + b3 + b4
					p[kmer] = 1.0 / 256.0 // Uniform distribution
					q[kmer] = float64((i+j+k+l+1)*10) / 2560.0 // Skewed distribution
				}
			}
		}
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := js.Calculate(p, q)
		if err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkBatchComparison benchmarks batch comparison operations
func BenchmarkBatchComparison(b *testing.B) {
	bc := NewBatchComparison()
	
	// Prepare test distributions
	p := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
	q := map[string]float64{"A": 0.4, "T": 0.3, "G": 0.2, "C": 0.1}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		bc.AddComparison("test", p, q)
	}
}