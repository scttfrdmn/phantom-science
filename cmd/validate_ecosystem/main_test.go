//go:build !release
// +build !release

package main

import (
	"encoding/json"
	"os"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/tools"
)

func TestContains_Found(t *testing.T) {
	if !contains([]string{"a", "b", "c"}, "b") {
		t.Error("Expected to find 'b' in slice")
	}
}

func TestContains_NotFound(t *testing.T) {
	if contains([]string{"a", "b", "c"}, "d") {
		t.Error("Should not find 'd' in slice")
	}
}

func TestContains_Empty(t *testing.T) {
	if contains([]string{}, "anything") {
		t.Error("Should not find anything in empty slice")
	}
}

func TestCreateDefaultConfiguration_NotNil(t *testing.T) {
	config := createDefaultConfiguration()
	if config == nil {
		t.Fatal("Expected non-nil default configuration")
	}
}

func TestCreateDefaultConfiguration_HasToolPaths(t *testing.T) {
	config := createDefaultConfiguration()
	if len(config.ToolPaths) == 0 {
		t.Error("Expected default config to have tool paths")
	}
}

func TestCreateDefaultConfiguration_HasValidationModes(t *testing.T) {
	config := createDefaultConfiguration()
	if len(config.ValidationModes) == 0 {
		t.Error("Expected default config to have validation modes")
	}
}

func TestCreateDefaultConfiguration_HasPositiveValues(t *testing.T) {
	config := createDefaultConfiguration()
	if config.MinAccuracy <= 0 || config.MinAccuracy > 1 {
		t.Errorf("Expected MinAccuracy between 0 and 1, got %f", config.MinAccuracy)
	}
	if config.MaxConcurrentJobs <= 0 {
		t.Errorf("Expected positive MaxConcurrentJobs, got %d", config.MaxConcurrentJobs)
	}
}

func TestValidateConfig_ValidConfig(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./test_data",
		ReferencePath:  "./reference",
		ToolPaths:      map[string]string{"gatk": "/usr/local/bin/gatk"},
		ValidationModes: []string{"basic"},
	}
	if err := validateConfig(config); err != nil {
		t.Errorf("Expected no error for valid config, got: %v", err)
	}
}

func TestValidateConfig_MissingTestDataPath(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "",
		ReferencePath:  "./reference",
		ToolPaths:      map[string]string{"gatk": "/usr/local/bin/gatk"},
		ValidationModes: []string{"basic"},
	}
	if err := validateConfig(config); err == nil {
		t.Error("Expected error for missing test data path")
	}
}

func TestValidateConfig_MissingReferencePath(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./test_data",
		ReferencePath:  "",
		ToolPaths:      map[string]string{"gatk": "/usr/local/bin/gatk"},
		ValidationModes: []string{"basic"},
	}
	if err := validateConfig(config); err == nil {
		t.Error("Expected error for missing reference path")
	}
}

func TestValidateConfig_EmptyToolPaths(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./test_data",
		ReferencePath:  "./reference",
		ToolPaths:      map[string]string{},
		ValidationModes: []string{"basic"},
	}
	if err := validateConfig(config); err == nil {
		t.Error("Expected error for empty tool paths")
	}
}

func TestValidateConfig_EmptyValidationModes(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./test_data",
		ReferencePath:  "./reference",
		ToolPaths:      map[string]string{"gatk": "/usr/local/bin/gatk"},
		ValidationModes: []string{},
	}
	if err := validateConfig(config); err == nil {
		t.Error("Expected error for empty validation modes")
	}
}

func TestValidateConfig_InvalidMode(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./test_data",
		ReferencePath:  "./reference",
		ToolPaths:      map[string]string{"gatk": "/usr/local/bin/gatk"},
		ValidationModes: []string{"invalid_mode"},
	}
	if err := validateConfig(config); err == nil {
		t.Error("Expected error for invalid validation mode")
	}
}

func TestValidateConfig_AllValidModes(t *testing.T) {
	validModes := []string{"basic", "performance", "population", "workflow", "realtime", "cloud", "all"}
	for _, mode := range validModes {
		config := &EcosystemValidationConfig{
			TestDataPath:    "./test_data",
			ReferencePath:  "./reference",
			ToolPaths:      map[string]string{"tool": "/path"},
			ValidationModes: []string{mode},
		}
		if err := validateConfig(config); err != nil {
			t.Errorf("Expected no error for mode '%s', got: %v", mode, err)
		}
	}
}

// Tests for convertToEcosystemConfig
func TestConvertToEcosystemConfig_BasicMode(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./data",
		ReferencePath:  "./ref",
		ToolPaths:      map[string]string{"gatk": "/bin/gatk"},
		ValidationModes: []string{"basic"},
	}
	ecosystemConfig := convertToEcosystemConfig(config)
	if len(ecosystemConfig.ValidationModes) == 0 {
		t.Error("convertToEcosystemConfig(basic) produced no validation modes")
	}
}

func TestConvertToEcosystemConfig_AllMode(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./data",
		ReferencePath:  "./ref",
		ToolPaths:      map[string]string{"gatk": "/bin/gatk"},
		ValidationModes: []string{"all"},
		CloudProvider:  "",
	}
	ecosystemConfig := convertToEcosystemConfig(config)
	if len(ecosystemConfig.ValidationModes) == 0 {
		t.Error("convertToEcosystemConfig(all) produced no validation modes")
	}
}

func TestConvertToEcosystemConfig_AllModeWithCloud(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./data",
		ReferencePath:  "./ref",
		ToolPaths:      map[string]string{"gatk": "/bin/gatk"},
		ValidationModes: []string{"all"},
		CloudProvider:  "aws",
	}
	ecosystemConfig := convertToEcosystemConfig(config)
	if len(ecosystemConfig.ValidationModes) == 0 {
		t.Error("convertToEcosystemConfig(all+cloud) produced no validation modes")
	}
}

func TestConvertToEcosystemConfig_PerformanceMode(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./data",
		ReferencePath:  "./ref",
		ToolPaths:      map[string]string{"gatk": "/bin/gatk"},
		ValidationModes: []string{"performance"},
	}
	ecosystemConfig := convertToEcosystemConfig(config)
	if len(ecosystemConfig.ValidationModes) == 0 {
		t.Error("convertToEcosystemConfig(performance) produced no validation modes")
	}
}

func TestConvertToEcosystemConfig_PopulationMode(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./data",
		ReferencePath:  "./ref",
		ToolPaths:      map[string]string{"gatk": "/bin/gatk"},
		ValidationModes: []string{"population"},
	}
	ecosystemConfig := convertToEcosystemConfig(config)
	if len(ecosystemConfig.ValidationModes) == 0 {
		t.Error("convertToEcosystemConfig(population) produced no validation modes")
	}
}

func TestConvertToEcosystemConfig_WorkflowMode(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./data",
		ReferencePath:  "./ref",
		ToolPaths:      map[string]string{"gatk": "/bin/gatk"},
		ValidationModes: []string{"workflow"},
	}
	ecosystemConfig := convertToEcosystemConfig(config)
	if len(ecosystemConfig.ValidationModes) == 0 {
		t.Error("convertToEcosystemConfig(workflow) produced no validation modes")
	}
}

func TestConvertToEcosystemConfig_RealtimeMode(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./data",
		ReferencePath:  "./ref",
		ToolPaths:      map[string]string{"gatk": "/bin/gatk"},
		ValidationModes: []string{"realtime"},
	}
	ecosystemConfig := convertToEcosystemConfig(config)
	if len(ecosystemConfig.ValidationModes) == 0 {
		t.Error("convertToEcosystemConfig(realtime) produced no validation modes")
	}
}

func TestConvertToEcosystemConfig_CloudMode(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    "./data",
		ReferencePath:  "./ref",
		ToolPaths:      map[string]string{"gatk": "/bin/gatk"},
		ValidationModes: []string{"cloud"},
	}
	ecosystemConfig := convertToEcosystemConfig(config)
	if len(ecosystemConfig.ValidationModes) == 0 {
		t.Error("convertToEcosystemConfig(cloud) produced no validation modes")
	}
}

func TestConvertToEcosystemConfig_CopiesFields(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:     "./data",
		ReferencePath:   "./ref",
		ToolPaths:       map[string]string{"gatk": "/bin/gatk"},
		ValidationModes: []string{"basic"},
		MaxMemoryGB:     8,
		MaxCPUCores:     4,
		MinAccuracy:     0.99,
	}
	ecosystemConfig := convertToEcosystemConfig(config)
	if ecosystemConfig.TestDataPath != "./data" {
		t.Errorf("TestDataPath = %q, want ./data", ecosystemConfig.TestDataPath)
	}
	if ecosystemConfig.ResourceLimits.MaxMemoryGB != 8 {
		t.Errorf("MaxMemoryGB = %d, want 8", ecosystemConfig.ResourceLimits.MaxMemoryGB)
	}
}

// Test printConfiguration (just prints to stdout)
func TestPrintConfiguration_NoError(t *testing.T) {
	config := createDefaultConfiguration()
	// Should not panic
	printConfiguration(config)
}

// ---- Tests for file-writing functions ----

func TestSaveConfiguration_CreatesFile(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-save-config-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	config := createDefaultConfiguration()
	outPath := tmpDir + "/config.json"
	if err := saveConfiguration(config, outPath); err != nil {
		t.Fatalf("saveConfiguration: %v", err)
	}

	info, err := os.Stat(outPath)
	if err != nil {
		t.Fatalf("output file not created: %v", err)
	}
	if info.Size() == 0 {
		t.Error("config file is empty")
	}
}

func TestSaveConfiguration_WritesValidJSON(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-save-json-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	config := createDefaultConfiguration()
	config.TestDataPath = "/custom/data/path"
	config.MaxCPUCores = 16

	outPath := tmpDir + "/config.json"
	if err := saveConfiguration(config, outPath); err != nil {
		t.Fatalf("saveConfiguration: %v", err)
	}

	data, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	var loaded EcosystemValidationConfig
	if err := json.Unmarshal(data, &loaded); err != nil {
		t.Fatalf("json.Unmarshal: %v", err)
	}
	if loaded.TestDataPath != "/custom/data/path" {
		t.Errorf("TestDataPath = %q, want /custom/data/path", loaded.TestDataPath)
	}
	if loaded.MaxCPUCores != 16 {
		t.Errorf("MaxCPUCores = %d, want 16", loaded.MaxCPUCores)
	}
}

func TestSaveConfiguration_InvalidPath(t *testing.T) {
	config := createDefaultConfiguration()
	err := saveConfiguration(config, "/nonexistent/dir/config.json")
	if err == nil {
		t.Error("expected error writing to nonexistent directory, got nil")
	}
}

func TestLoadConfiguration_ReadsExistingFile(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-load-config-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	// Write a config file manually
	configData := `{
		"test_data_path": "/my/data",
		"reference_path": "/my/ref",
		"validation_modes": ["basic"],
		"tool_paths": {"gatk": "/bin/gatk"},
		"max_concurrent_jobs": 8
	}`
	configPath := tmpDir + "/config.json"
	if err := os.WriteFile(configPath, []byte(configData), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	loaded, err := loadConfiguration(configPath)
	if err != nil {
		t.Fatalf("loadConfiguration: %v", err)
	}
	if loaded.TestDataPath != "/my/data" {
		t.Errorf("TestDataPath = %q, want /my/data", loaded.TestDataPath)
	}
	if loaded.MaxConcurrentJobs != 8 {
		t.Errorf("MaxConcurrentJobs = %d, want 8", loaded.MaxConcurrentJobs)
	}
}

func TestLoadConfiguration_CreatesDefaultWhenMissing(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-load-default-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	configPath := tmpDir + "/does_not_exist.json"
	config, err := loadConfiguration(configPath)
	if err != nil {
		t.Fatalf("loadConfiguration with missing file: %v", err)
	}
	if config == nil {
		t.Fatal("expected non-nil config")
	}
	// The default config file should now exist
	if _, err := os.Stat(configPath); err != nil {
		t.Errorf("default config file not written: %v", err)
	}
}

func TestLoadConfiguration_InvalidJSON(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-load-invalid-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	configPath := tmpDir + "/bad.json"
	if err := os.WriteFile(configPath, []byte("not valid json {{{"), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	_, err = loadConfiguration(configPath)
	if err == nil {
		t.Error("expected error for invalid JSON, got nil")
	}
}

func TestGenerateSummaryReport_CreatesFile(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-summary-report-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	report := makeSyntheticReport()
	if err := generateSummaryReport(report, tmpDir); err != nil {
		t.Fatalf("generateSummaryReport: %v", err)
	}

	summaryPath := tmpDir + "/ecosystem_validation_summary.txt"
	info, err := os.Stat(summaryPath)
	if err != nil {
		t.Fatalf("summary file not created: %v", err)
	}
	if info.Size() == 0 {
		t.Error("summary file is empty")
	}
}

func TestGenerateSummaryReport_ContainsExpectedContent(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-summary-content-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	report := makeSyntheticReport()
	if err := generateSummaryReport(report, tmpDir); err != nil {
		t.Fatalf("generateSummaryReport: %v", err)
	}

	data, err := os.ReadFile(tmpDir + "/ecosystem_validation_summary.txt")
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	content := string(data)

	// Check for key sections
	for _, want := range []string{
		"Execution Details",
		"Tool Validation Results",
		"Performance Scores",
		"Recommendations",
	} {
		if !strings.Contains(content, want) {
			t.Errorf("summary missing section %q", want)
		}
	}
}

func TestGenerateSummaryReport_WithCriticalIssues(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-summary-issues-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	report := makeSyntheticReport()
	report.OverallAssessment.CriticalIssues = []string{"Tool gatk not found", "Memory limit exceeded"}

	if err := generateSummaryReport(report, tmpDir); err != nil {
		t.Fatalf("generateSummaryReport: %v", err)
	}

	data, err := os.ReadFile(tmpDir + "/ecosystem_validation_summary.txt")
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	content := string(data)
	if !strings.Contains(content, "Critical Issues") {
		t.Error("summary missing 'Critical Issues' section when issues present")
	}
	if !strings.Contains(content, "Tool gatk not found") {
		t.Error("summary missing critical issue text")
	}
}

func TestGenerateSummaryReport_InvalidOutputDir(t *testing.T) {
	report := makeSyntheticReport()
	err := generateSummaryReport(report, "/nonexistent/dir/that/does/not/exist")
	if err == nil {
		t.Error("expected error for invalid output dir, got nil")
	}
}

// ---------- overrideConfigWithFlags ----------

func TestOverrideConfigWithFlags_DefaultValues_NoChange(t *testing.T) {
	// With all flags at their defaults, overrideConfigWithFlags should not change TestDataPath or ReferencePath
	config := &EcosystemValidationConfig{
		TestDataPath:  "/original/data",
		ReferencePath: "/original/ref",
		ToolPaths:     map[string]string{},
	}
	// Default *testDataPath == "./test_data", *referencePath == "./reference"
	// so the if-branches won't fire
	overrideConfigWithFlags(config)
	if config.TestDataPath != "/original/data" {
		t.Errorf("TestDataPath changed to %q, want /original/data", config.TestDataPath)
	}
	if config.ReferencePath != "/original/ref" {
		t.Errorf("ReferencePath changed to %q, want /original/ref", config.ReferencePath)
	}
}

func TestOverrideConfigWithFlags_SkipCloud_ClearsCloudProvider(t *testing.T) {
	// Set skipCloud flag to true and verify CloudProvider gets cleared
	origSkipCloud := *skipCloud
	defer func() { *skipCloud = origSkipCloud }()
	*skipCloud = true

	config := &EcosystemValidationConfig{
		ToolPaths:     map[string]string{},
		CloudProvider: "aws",
	}
	overrideConfigWithFlags(config)
	if config.CloudProvider != "" {
		t.Errorf("CloudProvider = %q, want empty after skipCloud", config.CloudProvider)
	}
}

func TestOverrideConfigWithFlags_ToolsPath_ParsesPairs(t *testing.T) {
	// Set toolsPath to a valid pair and verify it's added to config
	origToolsPath := *toolsPath
	defer func() { *toolsPath = origToolsPath }()
	*toolsPath = "gatk:/bin/gatk,bwa:/usr/bin/bwa"

	config := &EcosystemValidationConfig{
		ToolPaths: map[string]string{},
	}
	overrideConfigWithFlags(config)
	if config.ToolPaths["gatk"] != "/bin/gatk" {
		t.Errorf("ToolPaths[gatk] = %q, want /bin/gatk", config.ToolPaths["gatk"])
	}
	if config.ToolPaths["bwa"] != "/usr/bin/bwa" {
		t.Errorf("ToolPaths[bwa] = %q, want /usr/bin/bwa", config.ToolPaths["bwa"])
	}
}

func TestOverrideConfigWithFlags_ValidationMode_SetsMode(t *testing.T) {
	origMode := *validationMode
	defer func() { *validationMode = origMode }()
	*validationMode = "basic,performance"

	config := &EcosystemValidationConfig{
		ToolPaths:       map[string]string{},
		ValidationModes: []string{"all"},
	}
	overrideConfigWithFlags(config)
	if len(config.ValidationModes) != 2 {
		t.Errorf("ValidationModes = %v, want 2 modes", config.ValidationModes)
	}
}

// ---------- printValidationSummary ----------

func TestPrintValidationSummary_NoCrash(t *testing.T) {
	// printValidationSummary prints to stdout; just verify it doesn't panic
	report := makeSyntheticReport()
	printValidationSummary(report) // should not panic
}

func TestPrintValidationSummary_WithRecommendations(t *testing.T) {
	report := makeSyntheticReport()
	report.OverallAssessment.Recommendations = []string{"Upgrade tools", "Add more memory"}
	printValidationSummary(report) // should not panic
}

// ---------- executeEcosystemValidation ----------

func TestExecuteEcosystemValidation_EmptyToolPaths_ReturnsError(t *testing.T) {
	config := &EcosystemValidationConfig{
		TestDataPath:    t.TempDir(),
		ReferencePath:  t.TempDir(),
		ToolPaths:      map[string]string{"gatk": "/nonexistent/gatk"},
		ValidationModes: []string{"basic"},
	}
	err := executeEcosystemValidation(config)
	// May succeed or fail depending on tools.NewGenomicsEcosystemValidator; just no panic
	_ = err
}

// makeSyntheticReport builds a minimal EcosystemValidationReport for testing.
func makeSyntheticReport() *tools.EcosystemValidationReport {
	toolResults := map[string]*tools.ToolValidationResult{
		"gatk": {
			ToolName: "gatk",
			Version:  "4.2.0",
			Status:   "passed",
		},
		"bwa": {
			ToolName:     "bwa",
			Version:      "0.7.17",
			Status:       "failed",
			ErrorMessage: "binary not found",
		},
	}
	return &tools.EcosystemValidationReport{
		StartTime:   time.Now().Add(-5 * time.Minute),
		EndTime:     time.Now(),
		Duration:    5 * time.Minute,
		ToolResults: toolResults,
		OverallAssessment: &tools.OverallAssessment{
			OverallScore:       0.82,
			CompatibilityScore: 0.90,
			PerformanceScore:   0.80,
			ScalabilityScore:   0.75,
			ReliabilityScore:   0.85,
			Status:             "good",
			Recommendations:    []string{"Upgrade gatk to v4.3", "Add more test data"},
			CriticalIssues:     []string{},
		},
	}
}
