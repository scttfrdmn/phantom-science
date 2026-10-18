// Package validation compares a synthetic k-mer frequency table with a
// reference table and reports pass/fail per test: Jensen-Shannon divergence,
// k-mer complexity, normalized Shannon entropy and GC content, each against
// a fixed absolute threshold. The thresholds are configuration defaults, not
// derived from data.
package validation

import (
	"bytes"
	"encoding/csv"
	"encoding/json"
	"encoding/xml"
	"fmt"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/errors"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/statistics"
)

// ValidationFramework runs the k-mer table comparisons described in the
// package comment. WassersteinThreshold, MinSampleSize and MaxSampleSize are
// stored but not used by any test.
type ValidationFramework struct {
	JSThreshold          float64
	WassersteinThreshold float64
	MinSampleSize        int
	MaxSampleSize        int
	jsCalculator         *statistics.JensenShannonDivergence
	batchComparison      *statistics.BatchComparison
}

// NewValidationFramework creates a new validation framework
func NewValidationFramework() *ValidationFramework {
	return &ValidationFramework{
		JSThreshold:          0.01,   // Target Jensen-Shannon divergence < 0.01
		WassersteinThreshold: 0.02,   // Target Wasserstein distance < 0.02
		MinSampleSize:        1000,   // Minimum sample size for validation
		MaxSampleSize:        100000, // Maximum sample size for validation
		jsCalculator:         statistics.NewJensenShannonDivergence(),
		batchComparison:      statistics.NewBatchComparison(),
	}
}

// ValidationResult represents the result of a validation test
type ValidationResult struct {
	TestName      string
	TestType      ValidationTestType
	Passed        bool
	Score         float64
	Threshold     float64
	Message       string
	Details       interface{}
	Timestamp     time.Time
	ExecutionTime time.Duration
}

// ValidationTestType defines the type of validation test
type ValidationTestType int

const (
	JSdivergenceTest ValidationTestType = iota
	WassersteinTest
	ComplexityTest
	EntropyTest
	GCContentTest
	KmerDistributionTest
	StatisticalTest
	BiologicalTest
)

// ValidationTestTypeString returns string representation of validation test type
func (vtt ValidationTestType) String() string {
	switch vtt {
	case JSdivergenceTest:
		return "Jensen-Shannon Divergence"
	case WassersteinTest:
		return "Wasserstein Distance"
	case ComplexityTest:
		return "K-mer Complexity"
	case EntropyTest:
		return "Shannon Entropy"
	case GCContentTest:
		return "GC Content"
	case KmerDistributionTest:
		return "K-mer Distribution"
	case StatisticalTest:
		return "Statistical Test"
	case BiologicalTest:
		return "Biological Authenticity"
	default:
		return "Unknown"
	}
}

// ValidationReport contains comprehensive validation results
type ValidationReport struct {
	OverallPassed      bool
	TotalTests         int
	PassedTests        int
	FailedTests        int
	AverageScore       float64
	Results            []ValidationResult
	Timestamp          time.Time
	ExecutionTime      time.Duration
	RecommendedActions []string
}

// ValidateKmerFrequencies validates k-mer frequency distributions
func (vf *ValidationFramework) ValidateKmerFrequencies(reference, synthetic *kmer.FrequencyTable) *ValidationResult {
	start := time.Now()

	// Get probability distributions
	refProbs := reference.GetProbabilityDistribution()
	synProbs := synthetic.GetProbabilityDistribution()

	// Calculate Jensen-Shannon divergence
	jsDivergence, err := vf.jsCalculator.Calculate(refProbs, synProbs)
	if err != nil {
		// Handle calculation error with proper error type
		genomicsErr := errors.WrapError(err, errors.ValidationError, "validation_framework", "validate_kmer_frequencies", "failed to calculate JS divergence")
		errors.HandleError(genomicsErr)
		return &ValidationResult{
			TestName:      "K-mer Frequency Validation",
			TestType:      JSdivergenceTest,
			Passed:        false,
			Score:         1.0, // Worst possible score
			Threshold:     vf.JSThreshold,
			Message:       fmt.Sprintf("Failed to calculate JS divergence: %v", err),
			Timestamp:     time.Now(),
			ExecutionTime: time.Since(start),
		}
	}

	passed := jsDivergence < vf.JSThreshold
	message := fmt.Sprintf("JS divergence: %.6f (threshold: %.6f)", jsDivergence, vf.JSThreshold)

	if passed {
		message += " - PASSED"
	} else {
		message += " - FAILED"
	}

	return &ValidationResult{
		TestName:  "K-mer Frequency Validation",
		TestType:  JSdivergenceTest,
		Passed:    passed,
		Score:     jsDivergence,
		Threshold: vf.JSThreshold,
		Message:   message,
		Details: map[string]interface{}{
			"js_divergence":   jsDivergence,
			"reference_kmers": len(refProbs),
			"synthetic_kmers": len(synProbs),
		},
		Timestamp:     time.Now(),
		ExecutionTime: time.Since(start),
	}
}

// ValidateKmerComplexity validates k-mer complexity metrics
func (vf *ValidationFramework) ValidateKmerComplexity(reference, synthetic *kmer.FrequencyTable) *ValidationResult {
	start := time.Now()

	refComplexity := reference.GetKmerComplexity()
	synComplexity := synthetic.GetKmerComplexity()

	// Absolute difference
	complexityDiff := abs(refComplexity - synComplexity)
	threshold := 0.1

	passed := complexityDiff < threshold
	message := fmt.Sprintf("Complexity difference: %.4f (ref: %.4f, syn: %.4f, threshold: %.4f)",
		complexityDiff, refComplexity, synComplexity, threshold)

	if passed {
		message += " - PASSED"
	} else {
		message += " - FAILED"
	}

	return &ValidationResult{
		TestName:  "K-mer Complexity Validation",
		TestType:  ComplexityTest,
		Passed:    passed,
		Score:     complexityDiff,
		Threshold: threshold,
		Message:   message,
		Details: map[string]interface{}{
			"reference_complexity": refComplexity,
			"synthetic_complexity": synComplexity,
			"difference":           complexityDiff,
		},
		Timestamp:     time.Now(),
		ExecutionTime: time.Since(start),
	}
}

// ValidateEntropy validates Shannon entropy of k-mer distributions
func (vf *ValidationFramework) ValidateEntropy(reference, synthetic *kmer.FrequencyTable) *ValidationResult {
	start := time.Now()

	refEntropy := reference.GetNormalizedEntropy()
	synEntropy := synthetic.GetNormalizedEntropy()

	// Absolute difference
	entropyDiff := abs(refEntropy - synEntropy)
	threshold := 0.05

	passed := entropyDiff < threshold
	message := fmt.Sprintf("Normalized entropy difference: %.4f (ref: %.4f, syn: %.4f, threshold: %.4f)",
		entropyDiff, refEntropy, synEntropy, threshold)

	if passed {
		message += " - PASSED"
	} else {
		message += " - FAILED"
	}

	return &ValidationResult{
		TestName:  "Shannon Entropy Validation",
		TestType:  EntropyTest,
		Passed:    passed,
		Score:     entropyDiff,
		Threshold: threshold,
		Message:   message,
		Details: map[string]interface{}{
			"reference_entropy": refEntropy,
			"synthetic_entropy": synEntropy,
			"difference":        entropyDiff,
		},
		Timestamp:     time.Now(),
		ExecutionTime: time.Since(start),
	}
}

// ValidateGCContent validates GC content similarity
func (vf *ValidationFramework) ValidateGCContent(reference, synthetic *kmer.FrequencyTable) *ValidationResult {
	start := time.Now()

	refGC := reference.GetGCContent()
	synGC := synthetic.GetGCContent()

	// Absolute difference in GC fraction
	gcDiff := abs(refGC - synGC)
	threshold := 0.02

	passed := gcDiff < threshold
	message := fmt.Sprintf("GC content difference: %.4f (ref: %.4f, syn: %.4f, threshold: %.4f)",
		gcDiff, refGC, synGC, threshold)

	if passed {
		message += " - PASSED"
	} else {
		message += " - FAILED"
	}

	return &ValidationResult{
		TestName:  "GC Content Validation",
		TestType:  GCContentTest,
		Passed:    passed,
		Score:     gcDiff,
		Threshold: threshold,
		Message:   message,
		Details: map[string]interface{}{
			"reference_gc": refGC,
			"synthetic_gc": synGC,
			"difference":   gcDiff,
		},
		Timestamp:     time.Now(),
		ExecutionTime: time.Since(start),
	}
}

// ValidateMultipleKmers validates multiple k-mer sizes simultaneously
func (vf *ValidationFramework) ValidateMultipleKmers(referenceTables, syntheticTables []*kmer.FrequencyTable) *ValidationReport {
	start := time.Now()

	if len(referenceTables) != len(syntheticTables) {
		return &ValidationReport{
			OverallPassed:      false,
			TotalTests:         0,
			PassedTests:        0,
			FailedTests:        1,
			AverageScore:       1.0,
			Results:            nil,
			Timestamp:          time.Now(),
			ExecutionTime:      time.Since(start),
			RecommendedActions: []string{"Ensure equal number of reference and synthetic k-mer tables"},
		}
	}

	var results []ValidationResult
	passedCount := 0
	totalScore := 0.0

	for i, refTable := range referenceTables {
		synTable := syntheticTables[i]

		// Validate k-mer frequencies
		freqResult := vf.ValidateKmerFrequencies(refTable, synTable)
		results = append(results, *freqResult)
		if freqResult.Passed {
			passedCount++
		}
		totalScore += freqResult.Score

		// Validate complexity
		complexityResult := vf.ValidateKmerComplexity(refTable, synTable)
		results = append(results, *complexityResult)
		if complexityResult.Passed {
			passedCount++
		}
		totalScore += complexityResult.Score

		// Validate entropy
		entropyResult := vf.ValidateEntropy(refTable, synTable)
		results = append(results, *entropyResult)
		if entropyResult.Passed {
			passedCount++
		}
		totalScore += entropyResult.Score

		// Validate GC content
		gcResult := vf.ValidateGCContent(refTable, synTable)
		results = append(results, *gcResult)
		if gcResult.Passed {
			passedCount++
		}
		totalScore += gcResult.Score
	}

	totalTests := len(results)
	failedTests := totalTests - passedCount
	overallPassed := failedTests == 0
	averageScore := totalScore / float64(totalTests)

	// Generate recommendations
	recommendations := vf.generateRecommendations(results)

	return &ValidationReport{
		OverallPassed:      overallPassed,
		TotalTests:         totalTests,
		PassedTests:        passedCount,
		FailedTests:        failedTests,
		AverageScore:       averageScore,
		Results:            results,
		Timestamp:          time.Now(),
		ExecutionTime:      time.Since(start),
		RecommendedActions: recommendations,
	}
}

// generateRecommendations generates recommendations based on validation results
func (vf *ValidationFramework) generateRecommendations(results []ValidationResult) []string {
	var recommendations []string

	// Count failures by type
	failureCount := make(map[ValidationTestType]int)
	for _, result := range results {
		if !result.Passed {
			failureCount[result.TestType]++
		}
	}

	// Generate specific recommendations
	if failureCount[JSdivergenceTest] > 0 {
		recommendations = append(recommendations,
			"High Jensen-Shannon divergence detected. Consider adjusting k-mer frequency modeling parameters.")
	}

	if failureCount[ComplexityTest] > 0 {
		recommendations = append(recommendations,
			"K-mer complexity mismatch detected. Review sequence generation algorithm for proper diversity.")
	}

	if failureCount[EntropyTest] > 0 {
		recommendations = append(recommendations,
			"Shannon entropy mismatch detected. Check for proper randomness in sequence generation.")
	}

	if failureCount[GCContentTest] > 0 {
		recommendations = append(recommendations,
			"GC content mismatch detected. Adjust base composition parameters in sequence generation.")
	}

	if len(recommendations) == 0 {
		recommendations = append(recommendations, "All validation tests passed successfully!")
	}

	return recommendations
}

// GetValidationThresholds returns the current validation thresholds
func (vf *ValidationFramework) GetValidationThresholds() map[string]float64 {
	return map[string]float64{
		"js_divergence":        vf.JSThreshold,
		"wasserstein_distance": vf.WassersteinThreshold,
		"complexity_diff":      0.1,
		"entropy_diff":         0.05,
		"gc_content_diff":      0.02,
	}
}

// SetValidationThresholds updates the validation thresholds
func (vf *ValidationFramework) SetValidationThresholds(thresholds map[string]float64) {
	if jsThreshold, exists := thresholds["js_divergence"]; exists {
		vf.JSThreshold = jsThreshold
	}
	if wassersteinThreshold, exists := thresholds["wasserstein_distance"]; exists {
		vf.WassersteinThreshold = wassersteinThreshold
	}
}

// ValidationConfig contains configuration for validation framework
type ValidationConfig struct {
	JSThreshold           float64
	WassersteinThreshold  float64
	MinSampleSize         int
	MaxSampleSize         int
	EnableDetailedReport  bool
	EnableRecommendations bool
}

// NewValidationFrameworkWithConfig creates a validation framework with custom configuration
func NewValidationFrameworkWithConfig(config ValidationConfig) *ValidationFramework {
	vf := NewValidationFramework()
	vf.JSThreshold = config.JSThreshold
	vf.WassersteinThreshold = config.WassersteinThreshold
	vf.MinSampleSize = config.MinSampleSize
	vf.MaxSampleSize = config.MaxSampleSize
	return vf
}

// ExportValidationReport exports validation report to different formats
func (vf *ValidationFramework) ExportValidationReport(report *ValidationReport, format string) ([]byte, error) {
	switch format {
	case "json":
		return vf.exportToJSON(report)
	case "csv":
		return vf.exportToCSV(report)
	case "xml":
		return vf.exportToXML(report)
	default:
		return nil, errors.NewValidationError("validation_framework", "export_validation_report", fmt.Sprintf("unsupported export format: %s", format))
	}
}

// exportToJSON exports validation report to JSON format
func (vf *ValidationFramework) exportToJSON(report *ValidationReport) ([]byte, error) {
	type jsonResult struct {
		TestName    string  `json:"test_name"`
		TestType    string  `json:"test_type"`
		Passed      bool    `json:"passed"`
		Score       float64 `json:"score"`
		Threshold   float64 `json:"threshold"`
		Message     string  `json:"message"`
		Timestamp   string  `json:"timestamp"`
		ExecutionMS int64   `json:"execution_time_ms"`
	}
	type jsonReport struct {
		OverallPassed      bool         `json:"overall_passed"`
		TotalTests         int          `json:"total_tests"`
		PassedTests        int          `json:"passed_tests"`
		FailedTests        int          `json:"failed_tests"`
		AverageScore       float64      `json:"average_score"`
		Timestamp          string       `json:"timestamp"`
		ExecutionMS        int64        `json:"execution_time_ms"`
		RecommendedActions []string     `json:"recommended_actions"`
		Results            []jsonResult `json:"results"`
	}

	jr := jsonReport{
		OverallPassed:      report.OverallPassed,
		TotalTests:         report.TotalTests,
		PassedTests:        report.PassedTests,
		FailedTests:        report.FailedTests,
		AverageScore:       report.AverageScore,
		Timestamp:          report.Timestamp.Format(time.RFC3339),
		ExecutionMS:        report.ExecutionTime.Milliseconds(),
		RecommendedActions: report.RecommendedActions,
	}
	for _, r := range report.Results {
		jr.Results = append(jr.Results, jsonResult{
			TestName:    r.TestName,
			TestType:    r.TestType.String(),
			Passed:      r.Passed,
			Score:       r.Score,
			Threshold:   r.Threshold,
			Message:     r.Message,
			Timestamp:   r.Timestamp.Format(time.RFC3339),
			ExecutionMS: r.ExecutionTime.Milliseconds(),
		})
	}
	return json.Marshal(jr)
}

// exportToCSV exports validation report to CSV format
func (vf *ValidationFramework) exportToCSV(report *ValidationReport) ([]byte, error) {
	var buf bytes.Buffer
	w := csv.NewWriter(&buf)

	if err := w.Write([]string{
		"test_name", "test_type", "passed", "score", "threshold", "message", "timestamp", "execution_time_ms",
	}); err != nil {
		return nil, err
	}
	for _, r := range report.Results {
		if err := w.Write([]string{
			r.TestName,
			r.TestType.String(),
			fmt.Sprintf("%v", r.Passed),
			fmt.Sprintf("%.6f", r.Score),
			fmt.Sprintf("%.6f", r.Threshold),
			r.Message,
			r.Timestamp.Format(time.RFC3339),
			fmt.Sprintf("%d", r.ExecutionTime.Milliseconds()),
		}); err != nil {
			return nil, err
		}
	}
	w.Flush()
	if err := w.Error(); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// exportToXML exports validation report to XML format
func (vf *ValidationFramework) exportToXML(report *ValidationReport) ([]byte, error) {
	type xmlResult struct {
		XMLName     xml.Name `xml:"result"`
		TestName    string   `xml:"test_name"`
		TestType    string   `xml:"test_type"`
		Passed      bool     `xml:"passed"`
		Score       float64  `xml:"score"`
		Threshold   float64  `xml:"threshold"`
		Message     string   `xml:"message"`
		Timestamp   string   `xml:"timestamp"`
		ExecutionMS int64    `xml:"execution_time_ms"`
	}
	type xmlReport struct {
		XMLName            xml.Name    `xml:"validation_report"`
		OverallPassed      bool        `xml:"overall_passed"`
		TotalTests         int         `xml:"total_tests"`
		PassedTests        int         `xml:"passed_tests"`
		FailedTests        int         `xml:"failed_tests"`
		AverageScore       float64     `xml:"average_score"`
		Timestamp          string      `xml:"timestamp"`
		ExecutionMS        int64       `xml:"execution_time_ms"`
		RecommendedActions []string    `xml:"recommended_actions>action"`
		Results            []xmlResult `xml:"results>result"`
	}

	xr := xmlReport{
		OverallPassed:      report.OverallPassed,
		TotalTests:         report.TotalTests,
		PassedTests:        report.PassedTests,
		FailedTests:        report.FailedTests,
		AverageScore:       report.AverageScore,
		Timestamp:          report.Timestamp.Format(time.RFC3339),
		ExecutionMS:        report.ExecutionTime.Milliseconds(),
		RecommendedActions: report.RecommendedActions,
	}
	for _, r := range report.Results {
		xr.Results = append(xr.Results, xmlResult{
			TestName:    r.TestName,
			TestType:    r.TestType.String(),
			Passed:      r.Passed,
			Score:       r.Score,
			Threshold:   r.Threshold,
			Message:     r.Message,
			Timestamp:   r.Timestamp.Format(time.RFC3339),
			ExecutionMS: r.ExecutionTime.Milliseconds(),
		})
	}
	return xml.MarshalIndent(xr, "", "  ")
}

// Helper function to calculate absolute value
func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

// ValidationMetrics contains comprehensive metrics for validation
type ValidationMetrics struct {
	JSdivergence           float64
	WassersteinDistance    float64
	ComplexityDifference   float64
	EntropyDifference      float64
	GCContentDifference    float64
	TotalVariationDistance float64
	HellingerDistance      float64
	BhattacharyyaDistance  float64
}

// CalculateComprehensiveMetrics computes every distance metric between the
// two tables without applying thresholds.
func (vf *ValidationFramework) CalculateComprehensiveMetrics(reference, synthetic *kmer.FrequencyTable) (*ValidationMetrics, error) {
	refProbs := reference.GetProbabilityDistribution()
	synProbs := synthetic.GetProbabilityDistribution()

	// Calculate Jensen-Shannon divergence
	jsDivergence, err := vf.jsCalculator.Calculate(refProbs, synProbs)
	if err != nil {
		return nil, errors.WrapError(err, errors.ValidationError, "validation_framework", "calculate_comprehensive_metrics", "failed to calculate JS divergence")
	}

	// Calculate other metrics
	comparison, err := vf.jsCalculator.CompareDistributions(refProbs, synProbs)
	if err != nil {
		return nil, errors.WrapError(err, errors.ValidationError, "validation_framework", "calculate_comprehensive_metrics", "failed to compare distributions")
	}

	metrics := &ValidationMetrics{
		JSdivergence:           jsDivergence,
		WassersteinDistance:    comparison.WassersteinDistance,
		ComplexityDifference:   abs(reference.GetKmerComplexity() - synthetic.GetKmerComplexity()),
		EntropyDifference:      abs(reference.GetNormalizedEntropy() - synthetic.GetNormalizedEntropy()),
		GCContentDifference:    abs(reference.GetGCContent() - synthetic.GetGCContent()),
		TotalVariationDistance: comparison.TotalVariation,
		HellingerDistance:      comparison.HellingerDistance,
		BhattacharyyaDistance:  comparison.BhattacharyyaDistance,
	}

	return metrics, nil
}
