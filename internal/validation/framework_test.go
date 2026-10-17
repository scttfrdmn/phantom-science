package validation

import (
	"fmt"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestNewFramework(t *testing.T) {
	config := types.ValidationConfig{
		Enable:          true,
		StatisticalTests: []string{"hardy_weinberg", "linkage_disequilibrium"},
		PopulationGeneticsTests: []string{"population_structure"},
		PipelineCompatibilityTests: []string{"gatk_validation"},
	}
	
	framework := NewFramework(config)
	if framework == nil {
		t.Fatal("Expected validation framework to be created")
	}
	if !framework.config.Enable {
		t.Error("Expected validation to be enabled")
	}
	if len(framework.config.StatisticalTests) != 2 {
		t.Errorf("Expected 2 statistical tests, got %d", len(framework.config.StatisticalTests))
	}
}

func TestFrameworkRunValidation(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{"hardy_weinberg_equilibrium"},
		PopulationGeneticsTests: []string{"population_structure"},
		PipelineCompatibilityTests: []string{"gatk_best_practices"},
	}
	
	framework := NewFramework(config)
	
	// Mock data - in real implementation this would be actual genomic data
	mockData := []byte("mock genomic data")
	
	results := framework.RunValidation(mockData)
	
	if results.Summary == "" {
		t.Fatal("Expected validation results to be returned")
	}
	
	// Check that tests were run
	expectedTests := len(config.StatisticalTests) + 
					 len(config.PopulationGeneticsTests) + 
					 len(config.PipelineCompatibilityTests)
	
	if len(results.Tests) != expectedTests {
		t.Errorf("Expected %d test results, got %d", expectedTests, len(results.Tests))
	}
	
	// Check overall structure
	if results.ComplianceLevel == "" {
		t.Error("Expected compliance level to be set")
	}
	if results.Summary == "" {
		t.Error("Expected summary to be set")
	}
}

func TestFrameworkDisabled(t *testing.T) {
	config := types.ValidationConfig{
		Enable: false,
		StatisticalTests: []string{"hardy_weinberg"},
	}
	
	framework := NewFramework(config)
	results := framework.RunValidation(nil)
	
	// Should return minimal results when disabled
	if results.OverallPass {
		t.Error("Expected validation to not pass when disabled")
	}
	if len(results.Tests) != 0 {
		t.Errorf("Expected no tests to run when disabled, got %d", len(results.Tests))
	}
}

func TestMockValidationTests(t *testing.T) {
	// Test validation framework with mock implementations
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{"hardy_weinberg_equilibrium"},
	}
	
	framework := NewFramework(config)
	results := framework.RunValidation([]byte("mock data"))
	
	// Basic validation that framework produces results
	if len(results.Tests) == 0 {
		t.Error("Expected at least one test result")
	}
	
	if results.ComplianceLevel == "" {
		t.Error("Expected compliance level to be set")
	}
}

func TestValidationWithRegulatoryCompliance(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{"hardy_weinberg"},
		RegulatoryCompliance: types.RegulatoryConfig{
			FDACompliance: true,
			EMACompliance: true,
			ValidationStudies: []string{"analytical_validation"},
			QualityMetrics: []string{"precision_recall"},
		},
	}
	
	framework := NewFramework(config)
	results := framework.RunValidation([]byte("test data"))
	
	if !results.RegulatoryReady {
		t.Error("Expected regulatory readiness with FDA/EMA compliance")
	}
	
	if results.ComplianceLevel != "regulatory" {
		t.Errorf("Expected compliance level 'regulatory', got %s", results.ComplianceLevel)
	}
}

func TestValidationWithAcademicCompliance(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{"hardy_weinberg", "linkage_disequilibrium"},
		AcademicValidation: types.AcademicConfig{
			PeerReviewReady: true,
			PublicationFormat: "nature_genetics",
			BenchmarkDatasets: []string{"1000_genomes"},
			StatisticalPower: 0.8,
		},
	}
	
	framework := NewFramework(config)
	results := framework.RunValidation([]byte("academic test data"))
	
	if !results.AcademicReady {
		t.Error("Expected academic readiness with peer review settings")
	}
}

func TestValidationResultTimestamp(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{"hardy_weinberg_equilibrium"},
	}
	
	framework := NewFramework(config)
	results := framework.RunValidation([]byte("test"))
	
	if len(results.Tests) > 0 && results.Tests[0].Timestamp.IsZero() {
		t.Error("Expected timestamp to be set in validation result")
	}
}

func TestEmptyTestConfiguration(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{},
		PopulationGeneticsTests: []string{},
		PipelineCompatibilityTests: []string{},
	}
	
	framework := NewFramework(config)
	results := framework.RunValidation([]byte("test data"))
	
	if len(results.Tests) != 0 {
		t.Errorf("Expected no tests with empty configuration, got %d", len(results.Tests))
	}
	
	// Should still have overall structure
	if results.Summary == "" {
		t.Error("Expected summary even with no tests")
	}
}

func TestUnknownTestType(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{"unknown_test_type"},
	}
	
	framework := NewFramework(config)
	results := framework.RunValidation([]byte("test data"))
	
	// Should handle unknown tests gracefully
	if len(results.Tests) > 1 {
		t.Error("Expected unknown tests to be skipped")
	}
}

func TestValidationResultConsistency(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{"hardy_weinberg_equilibrium"},
	}
	
	framework := NewFramework(config)
	
	// Run validation multiple times with same data
	data := []byte("consistent test data")
	
	results1 := framework.RunValidation(data)
	results2 := framework.RunValidation(data)
	
	if len(results1.Tests) != len(results2.Tests) {
		t.Error("Expected consistent number of test results")
	}
	
	if results1.Tests[0].TestName != results2.Tests[0].TestName {
		t.Error("Expected consistent test names")
	}
}

func TestValidationPerformance(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{
			"hardy_weinberg_equilibrium",
			"linkage_disequilibrium",
		},
		PopulationGeneticsTests: []string{
			"population_structure",
		},
		PipelineCompatibilityTests: []string{
			"gatk_best_practices",
		},
	}
	
	framework := NewFramework(config)
	
	// Large mock dataset
	largeData := make([]byte, 10000)
	for i := range largeData {
		largeData[i] = byte('A' + (i % 4)) // Mock DNA sequence
	}
	
	// Validation should complete in reasonable time
	results := framework.RunValidation(largeData)
	
	if results.Summary == "" {
		t.Fatal("Expected validation to complete with large dataset")
	}
	
	if len(results.Tests) == 0 {
		t.Error("Expected tests to run with large dataset")
	}
}

// Test individual validation methods
func TestValidateFASTQFormat(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)
	
	// Valid FASTQ format
	validFASTQ := `@read1
ATCGATCGATCG
+
IIIIIIIIIIII
@read2
GCTAGCTAGCTA
+
HHHHHHHHHHHH`
	
	result := framework.runStatisticalTest("format_compliance", validFASTQ)
	if !result.Passed {
		t.Errorf("Valid FASTQ should pass format validation: %s", result.Message)
	}
	
	// Invalid FASTQ format (wrong line count)
	invalidFASTQ := `@read1
ATCGATCG
+`
	
	result = framework.runStatisticalTest("format_compliance", invalidFASTQ)
	if result.Passed {
		t.Error("Invalid FASTQ should fail format validation")
	}
	
	// Invalid data type
	result = framework.runStatisticalTest("format_compliance", 123)
	if result.Passed {
		t.Error("Non-string data should fail format validation")
	}
}

func TestValidateSequenceComplexity(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)
	
	// Complex sequence (should pass)
	complexFASTQ := `@read1
ATCGATCGATCG
+
IIIIIIIIIIII`
	
	result := framework.runStatisticalTest("sequence_complexity", complexFASTQ)
	if !result.Passed {
		t.Errorf("Complex sequence should pass: %s", result.Message)
	}
	
	// Homopolymer sequence (should fail)
	homopolymerFASTQ := `@read1
AAAAAAAAAA
+
IIIIIIIIII`
	
	result = framework.runStatisticalTest("sequence_complexity", homopolymerFASTQ)
	if result.Passed {
		t.Error("Homopolymer sequence should fail complexity check")
	}
	
	// Invalid data type
	result = framework.runStatisticalTest("sequence_complexity", []int{1, 2, 3})
	if result.Passed {
		t.Error("Non-string data should fail complexity validation")
	}
}

func TestValidateBaseComposition(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)
	
	// Valid DNA bases
	validFASTQ := `@read1
ATCGATCGN
+
IIIIIIIII`
	
	result := framework.runPopulationGeneticsTest("base_composition", validFASTQ)
	if !result.Passed {
		t.Errorf("Valid DNA bases should pass: %s", result.Message)
	}
	
	// Invalid DNA base
	invalidFASTQ := `@read1
ATCGATCGX
+
IIIIIIIII`
	
	result = framework.runPopulationGeneticsTest("base_composition", invalidFASTQ)
	if result.Passed {
		t.Error("Invalid DNA base should fail validation")
	}
	
	// Invalid data type
	result = framework.runPopulationGeneticsTest("base_composition", map[string]int{"test": 1})
	if result.Passed {
		t.Error("Non-string data should fail base composition validation")
	}
}

func TestValidateGCContent(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)
	
	// Normal GC content (50%)
	normalGC := `@read1
ATCGATCG
+
IIIIIIII`
	
	result := framework.runPopulationGeneticsTest("gc_content", normalGC)
	if !result.Passed {
		t.Errorf("Normal GC content should pass: %s", result.Message)
	}
	
	// Very low GC content (0%)
	lowGC := `@read1
ATATATATAT
+
IIIIIIIIII`
	
	result = framework.runPopulationGeneticsTest("gc_content", lowGC)
	if result.Passed {
		t.Error("Very low GC content should fail validation")
	}
	
	// Very high GC content (100%)
	highGC := `@read1
GCGCGCGCGC
+
IIIIIIIIII`
	
	result = framework.runPopulationGeneticsTest("gc_content", highGC)
	if result.Passed {
		t.Error("Very high GC content should fail validation")
	}
	
	// No valid bases
	noBases := `@read1
NNNNNNNNNN
+
IIIIIIIIII`
	
	result = framework.runPopulationGeneticsTest("gc_content", noBases)
	if result.Passed {
		t.Error("No valid bases should fail GC content validation")
	}
	
	// Invalid data type
	result = framework.runPopulationGeneticsTest("gc_content", 42.5)
	if result.Passed {
		t.Error("Non-string data should fail GC content validation")
	}
}

func TestValidateHeaderFormat(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)
	
	// Valid header format
	validFASTQ := `@read1_illumina
ATCGATCG
+
IIIIIIII`
	
	result := framework.runPipelineCompatibilityTest("header_format", validFASTQ)
	if !result.Passed {
		t.Errorf("Valid header should pass: %s", result.Message)
	}
	
	// Invalid header format (missing @)
	invalidFASTQ := `read1_no_prefix
ATCGATCG
+
IIIIIIII`
	
	result = framework.runPipelineCompatibilityTest("header_format", invalidFASTQ)
	if result.Passed {
		t.Error("Invalid header format should fail validation")
	}
	
	// Invalid data type
	result = framework.runPipelineCompatibilityTest("header_format", []string{"test"})
	if result.Passed {
		t.Error("Non-string data should fail header format validation")
	}
}

func TestValidationTestsWithAllMethods(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{
			"format_compliance",
			"quality_distribution", 
			"sequence_complexity",
			"read_length_consistency",
		},
		PopulationGeneticsTests: []string{
			"hardy_weinberg",
			"base_composition",
			"gc_content",
		},
		PipelineCompatibilityTests: []string{
			"fastq_parsing",
			"quality_encoding",
			"header_format",
		},
	}
	
	framework := NewFramework(config)
	
	testFASTQ := `@read1_illumina
ATCGATCGATCG
+read1_illumina
IIIIIIIIIIII
@read2_illumina  
GCTAGCTAGCTA
+read2_illumina
HHHHHHHHHHHH`
	
	results := framework.RunValidation(testFASTQ)
	
	expectedTestCount := len(config.StatisticalTests) + 
						len(config.PopulationGeneticsTests) + 
						len(config.PipelineCompatibilityTests)
	
	if len(results.Tests) != expectedTestCount {
		t.Errorf("Expected %d test results, got %d", expectedTestCount, len(results.Tests))
	}
	
	// Check that all tests have proper structure
	for _, test := range results.Tests {
		if test.TestName == "" {
			t.Error("Test should have a name")
		}
		if test.Confidence < 0 || test.Confidence > 1 {
			t.Errorf("Test %s confidence should be 0-1, got %.2f", test.TestName, test.Confidence)
		}
		if test.Timestamp.IsZero() {
			t.Errorf("Test %s should have timestamp", test.TestName)
		}
	}
}

func TestUnknownTestTypes(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)
	
	// Test unknown statistical test
	result := framework.runStatisticalTest("unknown_stat_test", "data")
	if result.Passed {
		t.Error("Unknown statistical test should fail")
	}
	if result.TestName != "unknown_stat_test" {
		t.Error("Unknown test should preserve test name")
	}
	
	// Test unknown population genetics test
	result = framework.runPopulationGeneticsTest("unknown_popgen_test", "data")
	if result.Passed {
		t.Error("Unknown population genetics test should fail")
	}
	
	// Test unknown pipeline compatibility test
	result = framework.runPipelineCompatibilityTest("unknown_pipeline_test", "data")
	if result.Passed {
		t.Error("Unknown pipeline compatibility test should fail")
	}
}

func TestValidationSummaryGeneration(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{"format_compliance", "sequence_complexity"},
	}
	
	framework := NewFramework(config)
	
	// Test with passing data
	validFASTQ := `@read1
ATCGATCGATCG
+
IIIIIIIIIIII`
	
	results := framework.RunValidation(validFASTQ)
	if !results.OverallPass {
		t.Error("Valid FASTQ should result in overall pass")
	}
	if results.Summary == "" {
		t.Error("Summary should be generated")
	}
	
	// Test with failing data (homopolymer)
	invalidFASTQ := `@read1
AAAAAAAAAAAA
+
IIIIIIIIIIII`
	
	results = framework.RunValidation(invalidFASTQ)
	if results.OverallPass {
		t.Error("Invalid FASTQ should result in overall fail")
	}
	if results.Summary == "" {
		t.Error("Summary should be generated for failures")
	}
}

func TestValidationWithComplexRegulatoryConfig(t *testing.T) {
	config := types.ValidationConfig{
		Enable: true,
		StatisticalTests: []string{"format_compliance"},
		RegulatoryCompliance: types.RegulatoryConfig{
			FDACompliance: true,
			EMACompliance: true,
			ValidationStudies: []string{"analytical_validation", "clinical_validation"},
			QualityMetrics: []string{"precision", "recall", "f1_score"},
		},
	}
	
	framework := NewFramework(config)
	results := framework.RunValidation("@read1\nATCG\n+\nIIII")
	
	if results.ComplianceLevel != "regulatory" {
		t.Errorf("Expected regulatory compliance level, got %s", results.ComplianceLevel)
	}
	
	if !results.RegulatoryReady {
		t.Error("Should be regulatory ready with FDA and EMA compliance")
	}
}

func TestValidationEdgeCases(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)
	
	// Test empty string
	result := framework.runStatisticalTest("format_compliance", "")
	if result.Passed {
		t.Error("Empty string should fail format validation")
	}
	
	// Test single line (not valid FASTQ)
	result = framework.runStatisticalTest("format_compliance", "single line")
	if result.Passed {
		t.Error("Single line should fail FASTQ format validation")
	}
	
	// Test whitespace-only
	result = framework.runStatisticalTest("format_compliance", "   \n  \n  ")
	if result.Passed {
		t.Error("Whitespace-only should fail format validation")
	}
}

// --- Tests for real HWE implementation ---

func TestHardyWeinberg_ExplicitCounts_PassesHWE(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)

	// Perfect HWE at p=0.5: 250 HomRef, 500 Het, 250 HomAlt → χ²=0, p=1.0
	data := NewValidationData("")
	data.GenotypeCounts = []GenotypeCount{{NHomRef: 250, NHet: 500, NHomAlt: 250}}

	result := framework.runPopulationGeneticsTest("hardy_weinberg", data)
	if !result.Passed {
		t.Errorf("Perfect HWE genotypes should pass: %s", result.Message)
	}
}

func TestHardyWeinberg_ExplicitCounts_FailsHWE(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)

	// Extreme excess homozygosity (strong inbreeding): 450 HomRef, 100 Het, 450 HomAlt
	// pHat=0.5, expected: eAA=250, eAa=500, eaa=250 → χ²≈640, p≈0
	data := NewValidationData("")
	data.GenotypeCounts = []GenotypeCount{{NHomRef: 450, NHet: 100, NHomAlt: 450}}

	result := framework.runPopulationGeneticsTest("hardy_weinberg", data)
	if result.Passed {
		t.Errorf("Strongly inbred genotypes should fail HWE: %s", result.Message)
	}
	if result.Confidence > 0.05 {
		t.Errorf("p-value should be < 0.05 for extreme inbreeding, got %.4f", result.Confidence)
	}
}

func TestHardyWeinberg_MultipleVariants_PassesHWE(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)

	// Aggregate across 3 independent loci all in HWE
	data := NewValidationData("")
	data.GenotypeCounts = []GenotypeCount{
		{NHomRef: 250, NHet: 500, NHomAlt: 250}, // locus 1 p=0.5
		{NHomRef: 360, NHet: 480, NHomAlt: 160}, // locus 2 p=0.6
		{NHomRef: 160, NHet: 480, NHomAlt: 360}, // locus 3 p=0.4
	}

	result := framework.runPopulationGeneticsTest("hardy_weinberg", data)
	if !result.Passed {
		t.Errorf("Aggregated HWE-compliant loci should pass: %s", result.Message)
	}
}

func TestHardyWeinberg_NoData_Skipped(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)

	result := framework.runPopulationGeneticsTest("hardy_weinberg", nil)
	if !result.Passed {
		t.Errorf("No-data HWE test should be skipped (pass), got: %s", result.Message)
	}
	if result.Confidence != 0.5 {
		t.Errorf("Skipped HWE should report 0.5 confidence, got %.2f", result.Confidence)
	}
}

func TestHardyWeinberg_FASTQProxy(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)

	// 20 reads with distinct GC content: 5 AT-rich, 10 balanced, 5 GC-rich.
	// GC fractions: 5×0%, 10×50%, 5×100%. mean≈50%.
	// 0% < 35% → HomRef (5), 50% in [35%,65%] → Het (10), 100% > 65% → HomAlt (5)
	var sb strings.Builder
	for i := 0; i < 20; i++ {
		fmt.Fprintf(&sb, "@read%d\n", i)
		switch {
		case i < 5:
			sb.WriteString("AAAAAAAAAAAAAAAAAAAA\n") // 0% GC
		case i < 15:
			sb.WriteString("ATATATATATATAGCGCGCG\n") // 50% GC
		default:
			sb.WriteString("GCGCGCGCGCGCGCGCGCGC\n") // 100% GC
		}
		sb.WriteString("+\nIIIIIIIIIIIIIIIIIIII\n")
	}

	result := framework.runPopulationGeneticsTest("hardy_weinberg", sb.String())
	if result.TestName != "hardy_weinberg" {
		t.Error("Expected hardy_weinberg test result")
	}
	// With 5 HomRef, 10 Het, 5 HomAlt: pHat=0.5, expected: eAA=5, eAa=10, eaa=5 → χ²=0, passes
	if !result.Passed {
		t.Errorf("Balanced FASTQ proxy should pass HWE: %s", result.Message)
	}
}

// --- Tests for real read_length_consistency implementation ---

func TestReadLengthConsistency_SameLength(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)

	fastq := "@read1\nATCGATCG\n+\nIIIIIIII\n@read2\nGCTAGCTA\n+\nHHHHHHHH"
	result := framework.runStatisticalTest("read_length_consistency", fastq)
	if !result.Passed {
		t.Errorf("Consistent read lengths should pass: %s", result.Message)
	}
}

func TestReadLengthConsistency_DifferentLengths(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)

	// Read 1 is 8 bp, read 2 is 12 bp
	fastq := "@read1\nATCGATCG\n+\nIIIIIIII\n@read2\nGCTAGCTAGCTA\n+\nHHHHHHHHHHHH"
	result := framework.runStatisticalTest("read_length_consistency", fastq)
	if result.Passed {
		t.Errorf("Inconsistent read lengths should fail: %s", result.Message)
	}
}

func TestReadLengthConsistency_NilData(t *testing.T) {
	config := types.ValidationConfig{Enable: true}
	framework := NewFramework(config)

	result := framework.runStatisticalTest("read_length_consistency", nil)
	if !result.Passed {
		t.Errorf("Nil data should not fail read length consistency: %s", result.Message)
	}
}

// --- Tests for ValidationData type ---

func TestNewValidationData(t *testing.T) {
	data := NewValidationData("@r1\nATCG\n+\nIIII")
	if data == nil {
		t.Fatal("NewValidationData returned nil")
	}
	if data.FASTQ == "" {
		t.Error("Expected FASTQ to be set")
	}
	if len(data.GenotypeCounts) != 0 {
		t.Error("Expected empty GenotypeCounts by default")
	}
}

// --- extractFASTQString edge cases ---

func TestExtractFASTQString_PlainString(t *testing.T) {
	const input = "@r1\nACGT\n+\nIIII"
	got := extractFASTQString(input)
	if got != input {
		t.Errorf("expected %q, got %q", input, got)
	}
}

func TestExtractFASTQString_ByteSlice(t *testing.T) {
	input := []byte("@r1\nACGT\n+\nIIII")
	got := extractFASTQString(input)
	if got != string(input) {
		t.Errorf("expected %q, got %q", string(input), got)
	}
}

func TestExtractFASTQString_ValidationData(t *testing.T) {
	vd := &ValidationData{FASTQ: "@r1\nACGT\n+\nIIII"}
	got := extractFASTQString(vd)
	if got != vd.FASTQ {
		t.Errorf("expected %q, got %q", vd.FASTQ, got)
	}
}

func TestExtractFASTQString_NilReturnsEmpty(t *testing.T) {
	got := extractFASTQString(nil)
	if got != "" {
		t.Errorf("expected empty string for nil, got %q", got)
	}
}

func TestExtractFASTQString_UnknownTypeReturnsEmpty(t *testing.T) {
	got := extractFASTQString(42)
	if got != "" {
		t.Errorf("expected empty string for unsupported type, got %q", got)
	}
}

// --- computeGenotypeProxy edge cases ---

// makeFASTQRead builds a single FASTQ record with the given sequence.
func makeFASTQRead(seq string) string {
	return "@read\n" + seq + "\n+\n" + strings.Repeat("I", len(seq))
}

func TestComputeGenotypeProxy_TwoReadsTooFew(t *testing.T) {
	// Only 2 reads → should return ok=false
	fastq := makeFASTQRead("AAAAAAAAAA") + "\n" + makeFASTQRead("GGGGGGGGGG")
	_, _, _, ok := computeGenotypeProxy(fastq)
	if ok {
		t.Error("expected ok=false for fewer than 3 reads")
	}
}

func TestComputeGenotypeProxy_EmptyString(t *testing.T) {
	_, _, _, ok := computeGenotypeProxy("")
	if ok {
		t.Error("expected ok=false for empty FASTQ string")
	}
}

func TestComputeGenotypeProxy_AllSameGC_AllHet(t *testing.T) {
	// All reads have identical GC; all within ±0.15 of mean → all Het
	read := makeFASTQRead("ACGTACGTAC") // 40% GC, same for every read
	fastq := read + "\n" + read + "\n" + read
	nHomRef, nHet, nHomAlt, ok := computeGenotypeProxy(fastq)
	if !ok {
		t.Fatal("expected ok=true for 3 reads")
	}
	if nHet != 3 {
		t.Errorf("expected nHet=3, got nHomRef=%d nHet=%d nHomAlt=%d", nHomRef, nHet, nHomAlt)
	}
}

func TestComputeGenotypeProxy_PolarisedReads(t *testing.T) {
	// 3 AT-rich reads (low GC) + 3 balanced + 3 GC-rich → should produce HomRef, Het, HomAlt
	atRich := makeFASTQRead("AAATTTAAAT")    // ~10% GC
	balanced := makeFASTQRead("ACGTACGTAC") // ~40% GC
	gcRich := makeFASTQRead("GCGCGCGCGC")   // ~90% GC
	var sb strings.Builder
	for range 3 {
		sb.WriteString(atRich + "\n")
	}
	for range 3 {
		sb.WriteString(balanced + "\n")
	}
	for range 3 {
		sb.WriteString(gcRich + "\n")
	}
	nHomRef, nHet, nHomAlt, ok := computeGenotypeProxy(sb.String())
	if !ok {
		t.Fatal("expected ok=true")
	}
	if nHomRef == 0 || nHet == 0 || nHomAlt == 0 {
		t.Errorf("expected all three bins non-zero, got nHomRef=%d nHet=%d nHomAlt=%d",
			nHomRef, nHet, nHomAlt)
	}
	if nHomRef+nHet+nHomAlt != 9 {
		t.Errorf("counts should sum to 9, got %d", nHomRef+nHet+nHomAlt)
	}
}

func TestComputeGenotypeProxy_NonDNABases(t *testing.T) {
	// Reads with N (ambiguous) bases — total>0 only from A/T/G/C; Ns are ignored
	// 3 reads with only 'N' produce total=0 → no GC fractions → ok=false
	nRead := makeFASTQRead("NNNNNNNNNN")
	fastq := nRead + "\n" + nRead + "\n" + nRead
	_, _, _, ok := computeGenotypeProxy(fastq)
	if ok {
		t.Error("expected ok=false when all bases are N (no valid GC fractions)")
	}
}

// ---------------------------------------------------------------------------
// validateQualityDistribution
// ---------------------------------------------------------------------------

func makeValidFASTQ(seqs []string, quals []string) string {
	var sb strings.Builder
	for i, seq := range seqs {
		q := quals[i]
		fmt.Fprintf(&sb, "@read%d\n%s\n+\n%s\n", i+1, seq, q)
	}
	return strings.TrimRight(sb.String(), "\n")
}

func TestValidateQualityDistribution_Valid(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	// Mean Phred ≈ 35 ('?' = Phred 30, 'I' = Phred 40)
	fastq := makeValidFASTQ(
		[]string{"ACGT", "TTAG"},
		[]string{"????", "IIII"},
	)
	res := fw.validateQualityDistribution(fastq, time.Now())
	if !res.Passed {
		t.Errorf("expected Passed=true for valid quality distribution, got: %s", res.Message)
	}
}

func TestValidateQualityDistribution_Degenerate(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	// All quality characters identical → degenerate distribution
	fastq := makeValidFASTQ(
		[]string{"ACGT"},
		[]string{"BBBB"}, // all same character
	)
	res := fw.validateQualityDistribution(fastq, time.Now())
	if res.Passed {
		t.Error("expected Passed=false for degenerate quality distribution")
	}
}

func TestValidateQualityDistribution_TooLow(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	// All '!' = Phred 0 — mean well below 10; but also all-same → degenerate check fires first
	// Use a mix of low-quality chars to get mean < 10 without being degenerate
	fastq := makeValidFASTQ(
		[]string{"ACGT", "TTAG"},
		[]string{"!\"#$", "!\"#$"}, // Phred 0,1,2,3 — mean = 1.5
	)
	res := fw.validateQualityDistribution(fastq, time.Now())
	if res.Passed {
		t.Error("expected Passed=false for mean quality below minimum")
	}
}

func TestValidateQualityDistribution_NoData(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	res := fw.validateQualityDistribution("", time.Now())
	if !res.Passed {
		t.Error("expected Passed=true (skip) when no FASTQ data provided")
	}
}

// ---------------------------------------------------------------------------
// validateFASTQParsing
// ---------------------------------------------------------------------------

func TestValidateFASTQParsing_Valid(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	fastq := makeValidFASTQ(
		[]string{"ACGT", "TTAG"},
		[]string{"IIII", "HHHH"},
	)
	res := fw.validateFASTQParsing(fastq, time.Now())
	if !res.Passed {
		t.Errorf("expected Passed=true for valid FASTQ, got: %s", res.Message)
	}
}

func TestValidateFASTQParsing_NotMultipleOf4(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	// 6 lines — not a multiple of 4
	fastq := "@read1\nACGT\n+\nIIII\n@read2\nTTAG"
	res := fw.validateFASTQParsing(fastq, time.Now())
	if res.Passed {
		t.Error("expected Passed=false for non-multiple-of-4 line count")
	}
}

func TestValidateFASTQParsing_MissingAt(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	fastq := "read1\nACGT\n+\nIIII" // header missing '@'
	res := fw.validateFASTQParsing(fastq, time.Now())
	if res.Passed {
		t.Error("expected Passed=false when header line does not start with '@'")
	}
}

func TestValidateFASTQParsing_MissingPlus(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	fastq := "@read1\nACGT\n-\nIIII" // separator is '-' not '+'
	res := fw.validateFASTQParsing(fastq, time.Now())
	if res.Passed {
		t.Error("expected Passed=false when separator line does not start with '+'")
	}
}

func TestValidateFASTQParsing_LengthMismatch(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	fastq := "@read1\nACGT\n+\nIII" // qual 3 chars but seq 4 chars
	res := fw.validateFASTQParsing(fastq, time.Now())
	if res.Passed {
		t.Error("expected Passed=false when sequence and quality lengths differ")
	}
}

func TestValidateFASTQParsing_NoData(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	res := fw.validateFASTQParsing("", time.Now())
	if !res.Passed {
		t.Error("expected Passed=true (skip) when no FASTQ data provided")
	}
}

// ---------------------------------------------------------------------------
// validateQualityEncoding
// ---------------------------------------------------------------------------

func TestValidateQualityEncoding_Valid(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	fastq := makeValidFASTQ(
		[]string{"ACGT"},
		[]string{"IIII"}, // 'I' = ASCII 73, within 33–126
	)
	res := fw.validateQualityEncoding(fastq, time.Now())
	if !res.Passed {
		t.Errorf("expected Passed=true for valid Phred+33 encoding, got: %s", res.Message)
	}
}

func TestValidateQualityEncoding_OutOfRange(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	// Quality line contains a space (ASCII 32, below minimum of 33)
	fastq := "@read1\nACGT\n+\nII I"
	res := fw.validateQualityEncoding(fastq, time.Now())
	if res.Passed {
		t.Error("expected Passed=false for quality character below ASCII 33")
	}
}

func TestValidateQualityEncoding_NoData(t *testing.T) {
	fw := NewFramework(types.ValidationConfig{Enable: true})
	res := fw.validateQualityEncoding("", time.Now())
	if !res.Passed {
		t.Error("expected Passed=true (skip) when no FASTQ data provided")
	}
}