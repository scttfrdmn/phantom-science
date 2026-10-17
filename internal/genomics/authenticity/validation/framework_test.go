package validation

import (
	"encoding/json"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// Helper function to create a mock frequency table for testing
func createMockFrequencyTable(frequencies map[string]int, k int) *kmer.FrequencyTable {
	total := 0
	for _, count := range frequencies {
		total += count
	}
	
	return &kmer.FrequencyTable{
		K:           k,
		Frequencies: frequencies,
		TotalCount:  total,
		TotalKmers:  total,
	}
}

// TestNewValidationFramework tests framework creation
func TestNewValidationFramework(t *testing.T) {
	vf := NewValidationFramework()
	
	assert.Equal(t, 0.01, vf.JSThreshold)
	assert.Equal(t, 0.02, vf.WassersteinThreshold)
	assert.Equal(t, 1000, vf.MinSampleSize)
	assert.Equal(t, 100000, vf.MaxSampleSize)
	assert.NotNil(t, vf.jsCalculator)
	assert.NotNil(t, vf.batchComparison)
}

// TestValidationTestType_String tests string representation of validation test types
func TestValidationTestType_String(t *testing.T) {
	testCases := []struct {
		testType ValidationTestType
		expected string
	}{
		{JSdivergenceTest, "Jensen-Shannon Divergence"},
		{WassersteinTest, "Wasserstein Distance"},
		{ComplexityTest, "K-mer Complexity"},
		{EntropyTest, "Shannon Entropy"},
		{GCContentTest, "GC Content"},
		{KmerDistributionTest, "K-mer Distribution"},
		{StatisticalTest, "Statistical Test"},
		{BiologicalTest, "Biological Authenticity"},
		{ValidationTestType(999), "Unknown"},
	}
	
	for _, tc := range testCases {
		assert.Equal(t, tc.expected, tc.testType.String())
	}
}

// TestValidationFramework_ValidateKmerFrequencies tests k-mer frequency validation
func TestValidationFramework_ValidateKmerFrequencies(t *testing.T) {
	vf := NewValidationFramework()
	
	t.Run("identical distributions", func(t *testing.T) {
		// Create identical frequency tables
		frequencies := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25}
		reference := createMockFrequencyTable(frequencies, 2)
		synthetic := createMockFrequencyTable(frequencies, 2)
		
		result := vf.ValidateKmerFrequencies(reference, synthetic)
		
		assert.True(t, result.Passed)
		assert.Equal(t, "K-mer Frequency Validation", result.TestName)
		assert.Equal(t, JSdivergenceTest, result.TestType)
		assert.Less(t, result.Score, vf.JSThreshold)
		assert.Contains(t, result.Message, "PASSED")
		assert.NotNil(t, result.Details)
		assert.Greater(t, result.ExecutionTime, time.Duration(0))
	})
	
	t.Run("different distributions", func(t *testing.T) {
		// Create different frequency tables
		refFreq := map[string]int{"AT": 50, "GC": 30, "TA": 15, "CG": 5}
		synFreq := map[string]int{"AT": 10, "GC": 10, "TA": 40, "CG": 40}
		
		reference := createMockFrequencyTable(refFreq, 2)
		synthetic := createMockFrequencyTable(synFreq, 2)
		
		result := vf.ValidateKmerFrequencies(reference, synthetic)
		
		assert.False(t, result.Passed)
		assert.Greater(t, result.Score, vf.JSThreshold)
		assert.Contains(t, result.Message, "FAILED")
		
		details := result.Details.(map[string]interface{})
		assert.Contains(t, details, "js_divergence")
		assert.Contains(t, details, "reference_kmers")
		assert.Contains(t, details, "synthetic_kmers")
	})
	
	t.Run("empty distributions", func(t *testing.T) {
		// Create empty frequency tables to trigger error
		reference := createMockFrequencyTable(map[string]int{}, 2)
		synthetic := createMockFrequencyTable(map[string]int{"AT": 1}, 2)
		
		result := vf.ValidateKmerFrequencies(reference, synthetic)
		
		assert.False(t, result.Passed)
		assert.Equal(t, 1.0, result.Score) // Worst possible score
		assert.Contains(t, result.Message, "Failed to calculate JS divergence")
	})
}

// TestValidationFramework_ValidateKmerComplexity tests k-mer complexity validation
func TestValidationFramework_ValidateKmerComplexity(t *testing.T) {
	vf := NewValidationFramework()
	
	t.Run("similar complexity", func(t *testing.T) {
		// Create frequency tables with similar complexity
		refFreq := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25}
		synFreq := map[string]int{"AT": 24, "GC": 26, "TA": 25, "CG": 25}
		
		reference := createMockFrequencyTable(refFreq, 2)
		synthetic := createMockFrequencyTable(synFreq, 2)
		
		result := vf.ValidateKmerComplexity(reference, synthetic)
		
		assert.True(t, result.Passed)
		assert.Equal(t, "K-mer Complexity Validation", result.TestName)
		assert.Equal(t, ComplexityTest, result.TestType)
		assert.Less(t, result.Score, 0.1) // Below threshold
		assert.Contains(t, result.Message, "PASSED")
		
		details := result.Details.(map[string]interface{})
		assert.Contains(t, details, "reference_complexity")
		assert.Contains(t, details, "synthetic_complexity")
		assert.Contains(t, details, "difference")
	})
	
	t.Run("different complexity", func(t *testing.T) {
		// Create frequency tables with very different complexity
		refFreq := map[string]int{"AT": 100} // Low complexity (1 unique k-mer)
		synFreq := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25} // High complexity
		
		reference := createMockFrequencyTable(refFreq, 2)
		synthetic := createMockFrequencyTable(synFreq, 2)
		
		result := vf.ValidateKmerComplexity(reference, synthetic)
		
		assert.False(t, result.Passed)
		assert.Greater(t, result.Score, 0.1) // Above threshold
		assert.Contains(t, result.Message, "FAILED")
	})
}

// TestValidationFramework_ValidateEntropy tests entropy validation
func TestValidationFramework_ValidateEntropy(t *testing.T) {
	vf := NewValidationFramework()
	
	t.Run("similar entropy", func(t *testing.T) {
		// Create frequency tables with similar entropy
		refFreq := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25}
		synFreq := map[string]int{"AT": 24, "GC": 26, "TA": 25, "CG": 25}
		
		reference := createMockFrequencyTable(refFreq, 2)
		synthetic := createMockFrequencyTable(synFreq, 2)
		
		result := vf.ValidateEntropy(reference, synthetic)
		
		assert.True(t, result.Passed)
		assert.Equal(t, "Shannon Entropy Validation", result.TestName)
		assert.Equal(t, EntropyTest, result.TestType)
		assert.Less(t, result.Score, 0.05) // Below threshold
		assert.Contains(t, result.Message, "PASSED")
	})
	
	t.Run("different entropy", func(t *testing.T) {
		// Create frequency tables with very different entropy
		refFreq := map[string]int{"AT": 100} // Low entropy
		synFreq := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25} // High entropy
		
		reference := createMockFrequencyTable(refFreq, 2)
		synthetic := createMockFrequencyTable(synFreq, 2)
		
		result := vf.ValidateEntropy(reference, synthetic)
		
		assert.False(t, result.Passed)
		assert.Greater(t, result.Score, 0.05) // Above threshold
		assert.Contains(t, result.Message, "FAILED")
		
		details := result.Details.(map[string]interface{})
		assert.Contains(t, details, "reference_entropy")
		assert.Contains(t, details, "synthetic_entropy")
		assert.Contains(t, details, "difference")
	})
}

// TestValidationFramework_ValidateGCContent tests GC content validation
func TestValidationFramework_ValidateGCContent(t *testing.T) {
	vf := NewValidationFramework()
	
	t.Run("similar GC content", func(t *testing.T) {
		// Create frequency tables with similar GC content
		refFreq := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25} // 50% GC
		synFreq := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25} // 50% GC (identical)
		
		reference := createMockFrequencyTable(refFreq, 2)
		synthetic := createMockFrequencyTable(synFreq, 2)
		
		result := vf.ValidateGCContent(reference, synthetic)
		
		assert.True(t, result.Passed)
		assert.Equal(t, "GC Content Validation", result.TestName)
		assert.Equal(t, GCContentTest, result.TestType)
		assert.Less(t, result.Score, 0.02) // Below threshold
		assert.Contains(t, result.Message, "PASSED")
	})
	
	t.Run("different GC content", func(t *testing.T) {
		// Create frequency tables with very different GC content
		refFreq := map[string]int{"AT": 90, "TA": 10} // 0% GC
		synFreq := map[string]int{"GC": 50, "CG": 50} // 100% GC
		
		reference := createMockFrequencyTable(refFreq, 2)
		synthetic := createMockFrequencyTable(synFreq, 2)
		
		result := vf.ValidateGCContent(reference, synthetic)
		
		assert.False(t, result.Passed)
		assert.Greater(t, result.Score, 0.02) // Above threshold
		assert.Contains(t, result.Message, "FAILED")
		
		details := result.Details.(map[string]interface{})
		assert.Contains(t, details, "reference_gc")
		assert.Contains(t, details, "synthetic_gc")
		assert.Contains(t, details, "difference")
	})
}

// TestValidationFramework_ValidateMultipleKmers tests multiple k-mer validation
func TestValidationFramework_ValidateMultipleKmers(t *testing.T) {
	vf := NewValidationFramework()
	
	t.Run("successful validation", func(t *testing.T) {
		// Create similar frequency tables
		refFreq1 := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25}
		synFreq1 := map[string]int{"AT": 24, "GC": 26, "TA": 25, "CG": 25}
		
		refFreq2 := map[string]int{"ATC": 20, "GCA": 20, "TAT": 20, "CGC": 20, "GGG": 20}
		synFreq2 := map[string]int{"ATC": 19, "GCA": 21, "TAT": 20, "CGC": 20, "GGG": 20}
		
		referenceTables := []*kmer.FrequencyTable{
			createMockFrequencyTable(refFreq1, 2),
			createMockFrequencyTable(refFreq2, 3),
		}
		syntheticTables := []*kmer.FrequencyTable{
			createMockFrequencyTable(synFreq1, 2),
			createMockFrequencyTable(synFreq2, 3),
		}
		
		report := vf.ValidateMultipleKmers(referenceTables, syntheticTables)
		
		assert.True(t, report.OverallPassed)
		assert.Equal(t, 8, report.TotalTests) // 4 tests per k-mer size * 2 sizes
		assert.Equal(t, 8, report.PassedTests)
		assert.Equal(t, 0, report.FailedTests)
		assert.Less(t, report.AverageScore, 0.1) // Should be low for passing tests
		assert.Len(t, report.Results, 8)
		assert.Greater(t, report.ExecutionTime, time.Duration(0))
		assert.Contains(t, report.RecommendedActions, "All validation tests passed successfully!")
	})
	
	t.Run("mismatched table counts", func(t *testing.T) {
		referenceTables := []*kmer.FrequencyTable{
			createMockFrequencyTable(map[string]int{"AT": 50}, 2),
		}
		syntheticTables := []*kmer.FrequencyTable{
			createMockFrequencyTable(map[string]int{"AT": 50}, 2),
			createMockFrequencyTable(map[string]int{"ATC": 50}, 3), // Extra table
		}
		
		report := vf.ValidateMultipleKmers(referenceTables, syntheticTables)
		
		assert.False(t, report.OverallPassed)
		assert.Equal(t, 0, report.TotalTests)
		assert.Equal(t, 1, report.FailedTests)
		assert.Equal(t, 1.0, report.AverageScore)
		assert.Contains(t, report.RecommendedActions[0], "Ensure equal number")
	})
	
	t.Run("mixed results", func(t *testing.T) {
		// Create one passing and one failing comparison
		refFreq1 := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25} // Will pass
		synFreq1 := map[string]int{"AT": 24, "GC": 26, "TA": 25, "CG": 25}
		
		refFreq2 := map[string]int{"AT": 100} // Will fail (very different from synthetic)
		synFreq2 := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25}
		
		referenceTables := []*kmer.FrequencyTable{
			createMockFrequencyTable(refFreq1, 2),
			createMockFrequencyTable(refFreq2, 2),
		}
		syntheticTables := []*kmer.FrequencyTable{
			createMockFrequencyTable(synFreq1, 2),
			createMockFrequencyTable(synFreq2, 2),
		}
		
		report := vf.ValidateMultipleKmers(referenceTables, syntheticTables)
		
		assert.False(t, report.OverallPassed)
		assert.Equal(t, 8, report.TotalTests)
		assert.Less(t, report.PassedTests, report.TotalTests)
		assert.Greater(t, report.FailedTests, 0)
		assert.Greater(t, len(report.RecommendedActions), 0)
		assert.NotContains(t, report.RecommendedActions, "All validation tests passed")
	})
}

// TestValidationFramework_generateRecommendations tests recommendation generation
func TestValidationFramework_generateRecommendations(t *testing.T) {
	vf := NewValidationFramework()
	
	t.Run("all tests passed", func(t *testing.T) {
		results := []ValidationResult{
			{Passed: true, TestType: JSdivergenceTest},
			{Passed: true, TestType: ComplexityTest},
		}
		
		recommendations := vf.generateRecommendations(results)
		
		assert.Len(t, recommendations, 1)
		assert.Contains(t, recommendations[0], "All validation tests passed successfully!")
	})
	
	t.Run("specific test failures", func(t *testing.T) {
		results := []ValidationResult{
			{Passed: false, TestType: JSdivergenceTest},
			{Passed: false, TestType: ComplexityTest},
			{Passed: false, TestType: EntropyTest},
			{Passed: false, TestType: GCContentTest},
			{Passed: true, TestType: StatisticalTest},
		}
		
		recommendations := vf.generateRecommendations(results)
		
		assert.GreaterOrEqual(t, len(recommendations), 4)
		
		// Check for specific recommendations
		recText := joinStringSlice(recommendations)
		assert.Contains(t, recText, "Jensen-Shannon divergence")
		assert.Contains(t, recText, "complexity mismatch")
		assert.Contains(t, recText, "Shannon entropy")
		assert.Contains(t, recText, "GC content")
	})
}

// TestValidationFramework_ThresholdMethods tests threshold getter/setter methods
func TestValidationFramework_ThresholdMethods(t *testing.T) {
	vf := NewValidationFramework()
	
	t.Run("get thresholds", func(t *testing.T) {
		thresholds := vf.GetValidationThresholds()
		
		assert.Equal(t, vf.JSThreshold, thresholds["js_divergence"])
		assert.Equal(t, vf.WassersteinThreshold, thresholds["wasserstein_distance"])
		assert.Contains(t, thresholds, "complexity_diff")
		assert.Contains(t, thresholds, "entropy_diff")
		assert.Contains(t, thresholds, "gc_content_diff")
	})
	
	t.Run("set thresholds", func(t *testing.T) {
		newThresholds := map[string]float64{
			"js_divergence":        0.005,
			"wasserstein_distance": 0.015,
		}
		
		vf.SetValidationThresholds(newThresholds)
		
		assert.Equal(t, 0.005, vf.JSThreshold)
		assert.Equal(t, 0.015, vf.WassersteinThreshold)
	})
}

// TestNewValidationFrameworkWithConfig tests framework creation with custom config
func TestNewValidationFrameworkWithConfig(t *testing.T) {
	config := ValidationConfig{
		JSThreshold:           0.005,
		WassersteinThreshold:  0.015,
		MinSampleSize:         500,
		MaxSampleSize:         50000,
		EnableDetailedReport:  true,
		EnableRecommendations: true,
	}
	
	vf := NewValidationFrameworkWithConfig(config)
	
	assert.Equal(t, config.JSThreshold, vf.JSThreshold)
	assert.Equal(t, config.WassersteinThreshold, vf.WassersteinThreshold)
	assert.Equal(t, config.MinSampleSize, vf.MinSampleSize)
	assert.Equal(t, config.MaxSampleSize, vf.MaxSampleSize)
}

// TestValidationFramework_ExportValidationReport tests report export functionality
func TestValidationFramework_ExportValidationReport(t *testing.T) {
	vf := NewValidationFramework()
	report := &ValidationReport{
		OverallPassed: true,
		TotalTests:    5,
		PassedTests:   5,
		FailedTests:   0,
	}
	
	t.Run("unsupported format", func(t *testing.T) {
		_, err := vf.ExportValidationReport(report, "unsupported")
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "unsupported export format")
	})
	
	t.Run("json export", func(t *testing.T) {
		data, err := vf.ExportValidationReport(report, "json")
		require.NoError(t, err)
		assert.True(t, json.Valid(data))
		assert.Contains(t, string(data), "overall_passed")
		assert.Contains(t, string(data), "total_tests")
	})

	t.Run("csv export", func(t *testing.T) {
		data, err := vf.ExportValidationReport(report, "csv")
		require.NoError(t, err)
		assert.Contains(t, string(data), "test_name")
		assert.Contains(t, string(data), "passed")
	})

	t.Run("xml export", func(t *testing.T) {
		data, err := vf.ExportValidationReport(report, "xml")
		require.NoError(t, err)
		assert.Contains(t, string(data), "<validation_report>")
		assert.Contains(t, string(data), "<overall_passed>")
	})
}

// TestValidationFramework_CalculateComprehensiveMetrics tests comprehensive metrics calculation
func TestValidationFramework_CalculateComprehensiveMetrics(t *testing.T) {
	vf := NewValidationFramework()
	
	t.Run("successful calculation", func(t *testing.T) {
		refFreq := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25}
		synFreq := map[string]int{"AT": 24, "GC": 26, "TA": 25, "CG": 25}
		
		reference := createMockFrequencyTable(refFreq, 2)
		synthetic := createMockFrequencyTable(synFreq, 2)
		
		metrics, err := vf.CalculateComprehensiveMetrics(reference, synthetic)
		require.NoError(t, err)
		
		assert.GreaterOrEqual(t, metrics.JSdivergence, 0.0)
		assert.GreaterOrEqual(t, metrics.WassersteinDistance, 0.0)
		assert.GreaterOrEqual(t, metrics.ComplexityDifference, 0.0)
		assert.GreaterOrEqual(t, metrics.EntropyDifference, 0.0)
		assert.GreaterOrEqual(t, metrics.GCContentDifference, 0.0)
		assert.GreaterOrEqual(t, metrics.TotalVariationDistance, 0.0)
		assert.GreaterOrEqual(t, metrics.HellingerDistance, 0.0)
		assert.GreaterOrEqual(t, metrics.BhattacharyyaDistance, 0.0)
	})
	
	t.Run("error in JS calculation", func(t *testing.T) {
		// Create empty reference table to trigger error
		reference := createMockFrequencyTable(map[string]int{}, 2)
		synthetic := createMockFrequencyTable(map[string]int{"AT": 1}, 2)
		
		_, err := vf.CalculateComprehensiveMetrics(reference, synthetic)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "failed to calculate JS divergence")
	})
}

// TestAbsFunction tests the abs helper function
func TestAbsFunction(t *testing.T) {
	testCases := []struct {
		input    float64
		expected float64
	}{
		{5.0, 5.0},
		{-5.0, 5.0},
		{0.0, 0.0},
		{-0.0, 0.0},
		{3.14159, 3.14159},
		{-3.14159, 3.14159},
	}
	
	for _, tc := range testCases {
		result := abs(tc.input)
		assert.Equal(t, tc.expected, result)
	}
}

// Helper function to join string slice for testing
func joinStringSlice(slice []string) string {
	result := ""
	for _, s := range slice {
		result += s + " "
	}
	return result
}

// BenchmarkValidationFramework_ValidateKmerFrequencies benchmarks k-mer frequency validation
func BenchmarkValidationFramework_ValidateKmerFrequencies(b *testing.B) {
	vf := NewValidationFramework()
	
	// Create realistic frequency tables
	refFreq := map[string]int{
		"AAAA": 100, "AAAT": 95, "AAAC": 90, "AAAG": 85,
		"TTTA": 80, "TTTC": 75, "TTTG": 70, "TTTT": 65,
		"CCCA": 60, "CCCT": 55, "CCCG": 50, "CCCC": 45,
		"GGGA": 40, "GGGT": 35, "GGGC": 30, "GGGG": 25,
	}
	synFreq := map[string]int{
		"AAAA": 98, "AAAT": 93, "AAAC": 88, "AAAG": 83,
		"TTTA": 78, "TTTC": 73, "TTTG": 68, "TTTT": 63,
		"CCCA": 58, "CCCT": 53, "CCCG": 48, "CCCC": 43,
		"GGGA": 38, "GGGT": 33, "GGGC": 28, "GGGG": 23,
	}
	
	reference := createMockFrequencyTable(refFreq, 4)
	synthetic := createMockFrequencyTable(synFreq, 4)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		vf.ValidateKmerFrequencies(reference, synthetic)
	}
}

// BenchmarkValidationFramework_ValidateMultipleKmers benchmarks multiple k-mer validation
func BenchmarkValidationFramework_ValidateMultipleKmers(b *testing.B) {
	vf := NewValidationFramework()
	
	// Create multiple frequency tables
	referenceTables := make([]*kmer.FrequencyTable, 3)
	syntheticTables := make([]*kmer.FrequencyTable, 3)
	
	for i := 0; i < 3; i++ {
		refFreq := map[string]int{"AT": 25, "GC": 25, "TA": 25, "CG": 25}
		synFreq := map[string]int{"AT": 24, "GC": 26, "TA": 25, "CG": 25}
		
		referenceTables[i] = createMockFrequencyTable(refFreq, 2+i)
		syntheticTables[i] = createMockFrequencyTable(synFreq, 2+i)
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		vf.ValidateMultipleKmers(referenceTables, syntheticTables)
	}
}