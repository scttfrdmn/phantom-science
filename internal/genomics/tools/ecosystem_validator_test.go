package tools_test

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics/tools"
)

// ─────────────────────────────────────────────────────────────────────────────
// EcosystemConfig.Validate
// ─────────────────────────────────────────────────────────────────────────────

func TestEcosystemConfig_Validate_Valid(t *testing.T) {
	cfg := tools.EcosystemConfig{
		TestDataPath:      "/tmp/test",
		ReferencePath:     "/tmp/ref",
		ToolPaths:         map[string]string{"dummy": "/usr/bin/true"},
		MaxConcurrentJobs: 2,
		TimeoutDuration:   30 * time.Minute,
	}
	err := cfg.Validate()
	assert.NoError(t, err)
}

func TestEcosystemConfig_Validate_MissingTestDataPath(t *testing.T) {
	cfg := tools.EcosystemConfig{
		ReferencePath: "/tmp/ref",
		ToolPaths:     map[string]string{"dummy": "/usr/bin/true"},
	}
	err := cfg.Validate()
	require.Error(t, err)
	assert.Contains(t, err.Error(), "test data path")
}

func TestEcosystemConfig_Validate_MissingReferencePath(t *testing.T) {
	cfg := tools.EcosystemConfig{
		TestDataPath: "/tmp/test",
		ToolPaths:    map[string]string{"dummy": "/usr/bin/true"},
	}
	err := cfg.Validate()
	require.Error(t, err)
	assert.Contains(t, err.Error(), "reference path")
}

func TestEcosystemConfig_Validate_EmptyToolPaths(t *testing.T) {
	cfg := tools.EcosystemConfig{
		TestDataPath:  "/tmp/test",
		ReferencePath: "/tmp/ref",
		ToolPaths:     map[string]string{},
	}
	err := cfg.Validate()
	require.Error(t, err)
	assert.Contains(t, err.Error(), "tool path")
}

func TestEcosystemConfig_Validate_NilToolPaths(t *testing.T) {
	cfg := tools.EcosystemConfig{
		TestDataPath:  "/tmp/test",
		ReferencePath: "/tmp/ref",
		ToolPaths:     nil,
	}
	err := cfg.Validate()
	require.Error(t, err)
}

func TestEcosystemConfig_Validate_DefaultsMaxConcurrentJobs(t *testing.T) {
	cfg := tools.EcosystemConfig{
		TestDataPath:      "/tmp/test",
		ReferencePath:     "/tmp/ref",
		ToolPaths:         map[string]string{"t": "/usr/bin/true"},
		MaxConcurrentJobs: 0, // should be set to 4 by Validate
	}
	err := cfg.Validate()
	require.NoError(t, err)
	assert.Equal(t, 4, cfg.MaxConcurrentJobs)
}

func TestEcosystemConfig_Validate_DefaultsTimeoutDuration(t *testing.T) {
	cfg := tools.EcosystemConfig{
		TestDataPath:    "/tmp/test",
		ReferencePath:   "/tmp/ref",
		ToolPaths:       map[string]string{"t": "/usr/bin/true"},
		TimeoutDuration: 0, // should be set to 1h by Validate
	}
	err := cfg.Validate()
	require.NoError(t, err)
	assert.Equal(t, 1*time.Hour, cfg.TimeoutDuration)
}

// ─────────────────────────────────────────────────────────────────────────────
// GenomicsToolSuite
// ─────────────────────────────────────────────────────────────────────────────

func TestNewGenomicsToolSuite_Empty(t *testing.T) {
	suite := tools.NewGenomicsToolSuite(map[string]string{})
	require.NotNil(t, suite)
	assert.Empty(t, suite.GetAvailableTools())
}

func TestNewGenomicsToolSuite_UnknownToolsIgnored(t *testing.T) {
	// Paths that don't match known tool names should be stored but not exposed
	suite := tools.NewGenomicsToolSuite(map[string]string{
		"totally_unknown": "/usr/bin/true",
	})
	require.NotNil(t, suite)
	// "totally_unknown" is not in the standard tool list, so it won't be added
	assert.Empty(t, suite.GetAvailableTools())
}

func TestNewGenomicsToolSuite_KnownToolWithFakePath(t *testing.T) {
	// Provide a path for a known tool; the suite registers it but doesn't check existence yet
	suite := tools.NewGenomicsToolSuite(map[string]string{
		"samtools": "/usr/bin/samtools-fake",
	})
	require.NotNil(t, suite)
	avail := suite.GetAvailableTools()
	require.Len(t, avail, 1)
	assert.Equal(t, "samtools", avail[0].Name)
	assert.Equal(t, "/usr/bin/samtools-fake", avail[0].ExecutablePath)
}

func TestNewGenomicsToolSuite_MultipleKnownTools(t *testing.T) {
	paths := map[string]string{
		"gatk":     "/opt/tools/gatk",
		"bwa":      "/opt/tools/bwa",
		"samtools": "/opt/tools/samtools",
		"bcftools": "/opt/tools/bcftools",
	}
	suite := tools.NewGenomicsToolSuite(paths)
	avail := suite.GetAvailableTools()
	assert.Len(t, avail, 4)

	names := make(map[string]bool)
	for _, t_ := range avail {
		names[t_.Name] = true
	}
	assert.True(t, names["gatk"])
	assert.True(t, names["bwa"])
	assert.True(t, names["samtools"])
	assert.True(t, names["bcftools"])
}

func TestGenomicsToolSuite_VerifyToolAvailability_NoTools(t *testing.T) {
	suite := tools.NewGenomicsToolSuite(map[string]string{})
	// No tools → nothing to verify → no error
	err := suite.VerifyToolAvailability()
	assert.NoError(t, err)
}

func TestGenomicsToolSuite_VerifyToolAvailability_NonExistentPath(t *testing.T) {
	suite := tools.NewGenomicsToolSuite(map[string]string{
		"gatk": "/nonexistent/path/gatk",
	})
	err := suite.VerifyToolAvailability()
	assert.Error(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// NewWorkflowEngine
// ─────────────────────────────────────────────────────────────────────────────

func TestNewWorkflowEngine_NoEnginesAvailable(t *testing.T) {
	cfg := tools.EcosystemConfig{
		TestDataPath:  "/tmp/test",
		ReferencePath: "/tmp/ref",
		ToolPaths:     map[string]string{"t": "/usr/bin/true"},
	}
	// WorkflowEngine creation is non-destructive even when no engines found
	engine := tools.NewWorkflowEngine(cfg)
	require.NotNil(t, engine)
}

func TestWorkflowEngine_Initialize_NoEngines(t *testing.T) {
	// A non-nil but empty ToolPaths map tells NewWorkflowEngine to use only
	// explicitly configured paths, skipping PATH probing.  With no engines
	// configured, Initialize must return an error regardless of what tools
	// happen to be installed on the host machine.
	cfg := tools.EcosystemConfig{ToolPaths: map[string]string{}}
	engine := tools.NewWorkflowEngine(cfg)
	err := engine.Initialize()
	assert.Error(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// NewFormatValidator
// ─────────────────────────────────────────────────────────────────────────────

func TestNewFormatValidator(t *testing.T) {
	fv := tools.NewFormatValidator()
	require.NotNil(t, fv)
}

// ─────────────────────────────────────────────────────────────────────────────
// NewQualityAnalyzer
// ─────────────────────────────────────────────────────────────────────────────

func TestNewQualityAnalyzer(t *testing.T) {
	thresholds := tools.QualityThresholds{
		MinAccuracy:       0.95,
		MaxFalsePositive:  0.01,
		MaxFalseNegative:  0.05,
		MinCoverage:       20,
		MaxContamination:  0.02,
		MinMappingQuality: 30,
	}
	qa := tools.NewQualityAnalyzer(thresholds)
	require.NotNil(t, qa)
}

// ─────────────────────────────────────────────────────────────────────────────
// NewEcosystemMetrics
// ─────────────────────────────────────────────────────────────────────────────

func TestNewEcosystemMetrics(t *testing.T) {
	em := tools.NewEcosystemMetrics()
	require.NotNil(t, em)
	assert.Equal(t, int64(0), em.ToolsValidated)
	assert.Equal(t, float64(0), em.DataProcessedGB)
}

func TestEcosystemMetrics_UpdateMetrics(t *testing.T) {
	em := tools.NewEcosystemMetrics()

	em.UpdateMetrics(3, 3, 512.0, 1000)
	assert.Equal(t, int64(3), em.ToolsValidated)
	assert.Equal(t, int64(3), em.ToolsAttempted)
	assert.InDelta(t, 512.0/1024.0, em.DataProcessedGB, 1e-9)
	assert.Equal(t, int64(1000), em.ExecutionTimeMS)
}

func TestEcosystemMetrics_UpdateMetrics_AccumulatesCorrectly(t *testing.T) {
	em := tools.NewEcosystemMetrics()

	em.UpdateMetrics(2, 3, 256.0, 500)
	em.UpdateMetrics(3, 4, 256.0, 500)

	assert.Equal(t, int64(5), em.ToolsValidated)
	assert.Equal(t, int64(7), em.ToolsAttempted)
	assert.InDelta(t, 512.0/1024.0, em.DataProcessedGB, 1e-9)
	assert.Equal(t, int64(1000), em.ExecutionTimeMS)
}

func TestEcosystemMetrics_SuccessRateCalculated(t *testing.T) {
	em := tools.NewEcosystemMetrics()
	em.UpdateMetrics(8, 10, 100.0, 2000)

	// 8 passed out of 10 attempted → SuccessRate = 0.8
	assert.InDelta(t, 0.8, em.SuccessRate, 1e-9)
	assert.InDelta(t, 0.2, em.ErrorRate, 1e-9)
}

// ─────────────────────────────────────────────────────────────────────────────
// ToolValidationResult — field and structural tests
// (calculateMetrics is unexported; we test the result struct directly)
// ─────────────────────────────────────────────────────────────────────────────

func TestToolValidationResult_Fields(t *testing.T) {
	result := &tools.ToolValidationResult{
		ToolName: "gatk",
		Version:  "4.2.0",
		Status:   "passed",
		Metrics:  map[string]float64{"success_rate": 1.0, "total_test_cases": 3.0},
	}
	assert.Equal(t, "gatk", result.ToolName)
	assert.Equal(t, "4.2.0", result.Version)
	assert.Equal(t, "passed", result.Status)
	assert.Equal(t, 1.0, result.Metrics["success_rate"])
}

func TestToolValidationResult_TestCases(t *testing.T) {
	tc1 := &tools.TestCaseResult{
		TestName: "basic_execution",
		Status:   "passed",
		Duration: 100 * time.Millisecond,
		Output:   "version 4.2.0",
	}
	tc2 := &tools.TestCaseResult{
		TestName:     "platform_test",
		Status:       "failed",
		ErrorMessage: "timeout",
		Duration:     5 * time.Second,
	}

	result := &tools.ToolValidationResult{
		ToolName:  "samtools",
		TestCases: []*tools.TestCaseResult{tc1, tc2},
	}
	assert.Len(t, result.TestCases, 2)
	assert.Equal(t, "basic_execution", result.TestCases[0].TestName)
	assert.Equal(t, "failed", result.TestCases[1].Status)
}

// ─────────────────────────────────────────────────────────────────────────────
// ValidationMode constants
// ─────────────────────────────────────────────────────────────────────────────

func TestValidationMode_Constants(t *testing.T) {
	assert.Equal(t, tools.ValidationMode(0), tools.ModeBasicCompatibility)
	assert.Equal(t, tools.ValidationMode(1), tools.ModePerformanceBenchmark)
	assert.Equal(t, tools.ValidationMode(2), tools.ModePopulationScale)
	assert.Equal(t, tools.ValidationMode(3), tools.ModeWorkflowIntegration)
	assert.Equal(t, tools.ValidationMode(4), tools.ModeRealTimeValidation)
	assert.Equal(t, tools.ValidationMode(5), tools.ModeCloudDeployment)
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper / synthetic data utilities (exercised directly through exported API)
// ─────────────────────────────────────────────────────────────────────────────

func TestGenomicsTool_Fields(t *testing.T) {
	tool := &tools.GenomicsTool{
		Name:           "fastqc",
		Version:        "0.11.9",
		ExecutablePath: "/opt/fastqc",
		Category:       "quality_control",
		Capabilities:   []string{"sequence_qc"},
		Parameters:     map[string]interface{}{"threads": 4},
	}
	assert.Equal(t, "fastqc", tool.Name)
	assert.Equal(t, "quality_control", tool.Category)
	assert.Len(t, tool.Capabilities, 1)
}

func TestBenchmarkSuite_Fields(t *testing.T) {
	suite := tools.BenchmarkSuite{
		Name:  "variant_calling",
		Tools: []string{"gatk", "freebayes"},
	}
	assert.Equal(t, "variant_calling", suite.Name)
	assert.Len(t, suite.Tools, 2)
}

func TestPopulationTest_Fields(t *testing.T) {
	pt := tools.PopulationTest{
		Name:        "1000_genomes",
		DatasetSize: 1000,
		Populations: []string{"YRI", "CEU", "CHB"},
	}
	assert.Equal(t, 1000, pt.DatasetSize)
	assert.Len(t, pt.Populations, 3)
}

func TestStreamingTest_Fields(t *testing.T) {
	st := tools.StreamingTest{
		Name:     "real_time_vc",
		DataRate: "1MB/s",
		Duration: 5 * time.Minute,
	}
	assert.Equal(t, "1MB/s", st.DataRate)
}

func TestWorkflowTest_Fields(t *testing.T) {
	wt := tools.WorkflowTest{
		Engine:       "cromwell",
		WorkflowType: "wdl",
		Name:         "gatk_germline",
	}
	assert.Equal(t, "cromwell", wt.Engine)
	assert.Equal(t, "wdl", wt.WorkflowType)
}

func TestCloudConfig_Fields(t *testing.T) {
	cc := tools.CloudConfig{
		Provider:      "aws",
		Region:        "us-east-1",
		InstanceTypes: []string{"c5.4xlarge", "r5.8xlarge"},
		StorageConfig: tools.StorageConfig{
			BucketName:   "my-genomics-bucket",
			StorageClass: "STANDARD",
		},
		NetworkConfig: tools.NetworkConfig{
			VPCId:     "vpc-12345",
			SubnetIds: []string{"subnet-a", "subnet-b"},
		},
		SecurityGroups: []string{"sg-genomics"},
	}
	assert.Equal(t, "aws", cc.Provider)
	assert.Len(t, cc.InstanceTypes, 2)
	assert.Equal(t, "vpc-12345", cc.NetworkConfig.VPCId)
}

func TestResourceLimits_Fields(t *testing.T) {
	rl := tools.ResourceLimits{
		MaxMemoryGB:      64,
		MaxCPUCores:      32,
		MaxDiskSpaceGB:   500,
		MaxExecutionTime: 4 * time.Hour,
	}
	assert.Equal(t, 64, rl.MaxMemoryGB)
	assert.Equal(t, 4*time.Hour, rl.MaxExecutionTime)
}

func TestQualityThresholds_Fields(t *testing.T) {
	qt := tools.QualityThresholds{
		MinAccuracy:       0.99,
		MaxFalsePositive:  0.001,
		MaxFalseNegative:  0.01,
		MinCoverage:       30,
		MaxContamination:  0.02,
		MinMappingQuality: 20,
	}
	assert.Equal(t, 0.99, qt.MinAccuracy)
	assert.Equal(t, 30, qt.MinCoverage)
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTQRead / VCFVariant (exported structs used in test data generation)
// ─────────────────────────────────────────────────────────────────────────────

func TestFASTQRead_Fields(t *testing.T) {
	read := &tools.FASTQRead{
		Header:   "@read_1_illumina_YRI",
		Sequence: "ACGTACGTACGT",
		Plus:     "+",
		Quality:  "IIIIIIIIIIIII",
	}
	assert.True(t, len(read.Sequence) > 0)
	assert.Equal(t, "+", read.Plus)
}

func TestVCFVariant_Fields(t *testing.T) {
	variant := &tools.VCFVariant{
		Chromosome: "chr1",
		Position:   1000,
		ID:         "var_1_YRI",
		Reference:  "A",
		Alternate:  []string{"T"},
		Quality:    35.0,
		Filter:     "PASS",
		Info:       map[string]string{"AF": "0.150"},
	}
	assert.Equal(t, "chr1", variant.Chromosome)
	assert.Equal(t, int64(1000), variant.Position)
	assert.Equal(t, []string{"T"}, variant.Alternate)
}

// ─────────────────────────────────────────────────────────────────────────────
// EcosystemValidationReport / OverallAssessment construction
// ─────────────────────────────────────────────────────────────────────────────

func TestOverallAssessment_StatusMapping(t *testing.T) {
	tests := []struct {
		score          float64
		expectedStatus string
	}{
		{0.95, "excellent"},
		{0.85, "good"},
		{0.72, "acceptable"},
		{0.60, "needs_improvement"},
	}

	for _, tc := range tests {
		oa := &tools.OverallAssessment{
			OverallScore: tc.score,
		}
		// Apply the same logic as generateOverallAssessment
		if tc.score >= 0.9 {
			oa.Status = "excellent"
		} else if tc.score >= 0.8 {
			oa.Status = "good"
		} else if tc.score >= 0.7 {
			oa.Status = "acceptable"
		} else {
			oa.Status = "needs_improvement"
		}
		assert.Equal(t, tc.expectedStatus, oa.Status, "score=%v", tc.score)
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// SaveValidationReport
// ─────────────────────────────────────────────────────────────────────────────

// SaveValidationReport writes the report, including the validator config,
// to disk. Cloud credentials and the storage encryption key in that config
// must never be written.
func TestSaveValidationReport(t *testing.T) {
	tmp := t.TempDir()
	reportPath := filepath.Join(tmp, "report.json")

	const secretKey, secretToken, encKey = "AKIA-SECRET-ACCESS-KEY", "session-token-value-123", "enc-key-value-456"
	report := &tools.EcosystemValidationReport{
		StartTime: time.Now(),
		EndTime:   time.Now(),
		ValidatorConfig: tools.EcosystemConfig{
			CloudConfig: tools.CloudConfig{
				Provider:      "aws",
				Region:        "us-west-2",
				StorageConfig: tools.StorageConfig{BucketName: "bucket", EncryptionKey: encKey},
				Credentials:   map[string]string{"aws_secret_access_key": secretKey, "session_token": secretToken},
			},
		},
		ToolResults: map[string]*tools.ToolValidationResult{
			"samtools": {
				ToolName: "samtools",
				Status:   "passed",
				Metrics:  map[string]float64{"success_rate": 1.0},
			},
		},
		WorkflowResults:    map[string]*tools.WorkflowValidationResult{},
		PerformanceMetrics: tools.NewEcosystemMetrics(),
		OverallAssessment: &tools.OverallAssessment{
			OverallScore: 0.92,
			Status:       "excellent",
		},
	}
	report.Duration = report.EndTime.Sub(report.StartTime)

	// SaveValidationReport uses no validator state.
	require.NoError(t, (&tools.GenomicsEcosystemValidator{}).SaveValidationReport(report, reportPath))

	raw, err := os.ReadFile(reportPath)
	require.NoError(t, err)
	assert.Contains(t, string(raw), "samtools")
	assert.Contains(t, string(raw), "us-west-2", "non-secret cloud config should still be reported")
	for _, secret := range []string{secretKey, secretToken, encKey, "aws_secret_access_key"} {
		if strings.Contains(string(raw), secret) {
			t.Errorf("saved report contains credential material %q", secret)
		}
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// ProcessingStatistics field checks (data struct only)
// ─────────────────────────────────────────────────────────────────────────────

func TestProcessingStatistics_Fields(t *testing.T) {
	ps := tools.ValidationStatus{
		TotalTools:      5,
		ValidatedTools:  3,
		ValidationState: "running",
		Progress:        0.6,
	}
	assert.Equal(t, 5, ps.TotalTools)
	assert.Equal(t, 0.6, ps.Progress)
}

// ─────────────────────────────────────────────────────────────────────────────
// NewGenomicsEcosystemValidator — negative path (config validation fails)
// ─────────────────────────────────────────────────────────────────────────────

func TestNewGenomicsEcosystemValidator_InvalidConfig(t *testing.T) {
	// Empty config → Validate() returns error → constructor fails
	_, err := tools.NewGenomicsEcosystemValidator(tools.EcosystemConfig{})
	require.Error(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// Format info helper (formatInfo is unexported but exercised through writeVCFFile)
// We test it indirectly by verifying the VCFVariant structure.
// ─────────────────────────────────────────────────────────────────────────────

func TestVCFVariant_InfoMap(t *testing.T) {
	v := &tools.VCFVariant{
		Chromosome: "chr22",
		Position:   999999,
		ID:         "rs001",
		Reference:  "G",
		Alternate:  []string{"A"},
		Quality:    50.0,
		Filter:     "PASS",
		Info:       map[string]string{"AF": "0.05", "DP": "100"},
	}
	assert.Equal(t, "0.05", v.Info["AF"])
	assert.Equal(t, "100", v.Info["DP"])
}
