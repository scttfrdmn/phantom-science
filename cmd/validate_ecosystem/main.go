package main

import (
	"context"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/tools"
)

// EcosystemValidationConfig contains configuration for ecosystem validation
type EcosystemValidationConfig struct {
	// Tool paths
	ToolPaths map[string]string `json:"tool_paths"`
	
	// Test data configuration
	TestDataPath    string `json:"test_data_path"`
	ReferencePath   string `json:"reference_path"`
	
	// Validation settings
	ValidationModes   []string `json:"validation_modes"`
	PopulationPanels  []string `json:"population_panels"`
	PlatformTargets   []string `json:"platform_targets"`
	
	// Performance settings
	MaxConcurrentJobs int           `json:"max_concurrent_jobs"`
	TimeoutDuration   time.Duration `json:"timeout_duration"`
	
	// Resource limits
	MaxMemoryGB      int           `json:"max_memory_gb"`
	MaxCPUCores      int           `json:"max_cpu_cores"`
	MaxDiskSpaceGB   int           `json:"max_disk_space_gb"`
	MaxExecutionTime time.Duration `json:"max_execution_time"`
	
	// Quality thresholds
	MinAccuracy        float64 `json:"min_accuracy"`
	MaxFalsePositive   float64 `json:"max_false_positive"`
	MaxFalseNegative   float64 `json:"max_false_negative"`
	MinCoverage        int     `json:"min_coverage"`
	MaxContamination   float64 `json:"max_contamination"`
	MinMappingQuality  int     `json:"min_mapping_quality"`
	
	// Cloud configuration
	CloudProvider    string            `json:"cloud_provider"`
	CloudRegion      string            `json:"cloud_region"`
	CloudCredentials map[string]string `json:"cloud_credentials"`
	
	// Output settings
	ReportPath       string `json:"report_path"`
	DetailedReporting bool   `json:"detailed_reporting"`
	JSONOutput       bool   `json:"json_output"`
}

var (
	configFile    = flag.String("config", "ecosystem_validation.json", "Configuration file path")
	outputDir     = flag.String("output", "./ecosystem_validation_output", "Output directory for reports")
	validationMode = flag.String("mode", "all", "Validation mode: basic,performance,population,workflow,realtime,cloud,all")
	toolsPath     = flag.String("tools", "", "Comma-separated tool paths (tool:path)")
	testDataPath  = flag.String("testdata", "./test_data", "Test data directory")
	referencePath = flag.String("reference", "./reference", "Reference genome directory")
	skipCloud     = flag.Bool("skip-cloud", false, "Skip cloud deployment validation")
	verbose       = flag.Bool("verbose", false, "Enable verbose logging")
	dryRun        = flag.Bool("dry-run", false, "Perform dry run without actual execution")
)

func main() {
	flag.Parse()
	
	if *verbose {
		log.Println("🔬 Phantom Science Genomics Ecosystem Validator")
		log.Println("================================================")
	}
	
	// Load configuration
	config, err := loadConfiguration(*configFile)
	if err != nil {
		log.Fatalf("Failed to load configuration: %v", err)
	}
	
	// Override configuration with command-line flags
	overrideConfigWithFlags(config)
	
	// Validate configuration
	if err := validateConfig(config); err != nil {
		log.Fatalf("Configuration validation failed: %v", err)
	}
	
	if *dryRun {
		log.Println("🧪 Dry run mode - configuration validated successfully")
		printConfiguration(config)
		return
	}
	
	// Create output directory
	if err := os.MkdirAll(*outputDir, 0750); err != nil {
		log.Fatalf("Failed to create output directory: %v", err)
	}
	
	// Execute ecosystem validation
	if err := executeEcosystemValidation(config); err != nil {
		log.Fatalf("Ecosystem validation failed: %v", err)
	}
	
	log.Println("✅ Ecosystem validation completed successfully!")
}

// loadConfiguration loads configuration from file or creates default
func loadConfiguration(configPath string) (*EcosystemValidationConfig, error) {
	// Check if config file exists
	if _, err := os.Stat(configPath); os.IsNotExist(err) {
		if *verbose {
			log.Printf("Configuration file not found, creating default: %s", configPath)
		}
		
		config := createDefaultConfiguration()
		
		// Save default configuration
		if err := saveConfiguration(config, configPath); err != nil {
			return nil, fmt.Errorf("failed to save default configuration: %w", err)
		}
		
		return config, nil
	}
	
	// Load existing configuration
	file, err := os.Open(configPath) // #nosec G304
	if err != nil {
		return nil, fmt.Errorf("failed to open configuration file: %w", err)
	}
	defer file.Close()
	
	var config EcosystemValidationConfig
	decoder := json.NewDecoder(file)
	if err := decoder.Decode(&config); err != nil {
		return nil, fmt.Errorf("failed to parse configuration: %w", err)
	}
	
	return &config, nil
}

// createDefaultConfiguration creates a default configuration
func createDefaultConfiguration() *EcosystemValidationConfig {
	return &EcosystemValidationConfig{
		ToolPaths: map[string]string{
			"gatk":      "/usr/local/bin/gatk",
			"bwa":       "/usr/local/bin/bwa",
			"bcftools":  "/usr/local/bin/bcftools",
			"samtools":  "/usr/local/bin/samtools",
			"freebayes": "/usr/local/bin/freebayes",
			"manta":     "/usr/local/bin/manta",
			"fastqc":    "/usr/local/bin/fastqc",
			"multiqc":   "/usr/local/bin/multiqc",
		},
		TestDataPath:      "./test_data",
		ReferencePath:     "./reference",
		ValidationModes:   []string{"basic", "performance", "workflow"},
		PopulationPanels:  []string{"EUR", "AFR", "EAS", "AMR", "SAS"},
		PlatformTargets:   []string{"illumina", "pacbio", "nanopore"},
		MaxConcurrentJobs: 4,
		TimeoutDuration:   1 * time.Hour,
		MaxMemoryGB:       32,
		MaxCPUCores:       8,
		MaxDiskSpaceGB:    100,
		MaxExecutionTime:  2 * time.Hour,
		MinAccuracy:       0.95,
		MaxFalsePositive:  0.05,
		MaxFalseNegative:  0.05,
		MinCoverage:       20,
		MaxContamination:  0.02,
		MinMappingQuality: 20,
		CloudProvider:     "",
		CloudRegion:       "",
		CloudCredentials:  make(map[string]string),
		ReportPath:        "./ecosystem_validation_report.json",
		DetailedReporting: true,
		JSONOutput:        true,
	}
}

// saveConfiguration saves configuration to file
func saveConfiguration(config *EcosystemValidationConfig, configPath string) error {
	file, err := os.Create(configPath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()
	
	encoder := json.NewEncoder(file)
	encoder.SetIndent("", "  ")
	return encoder.Encode(config)
}

// overrideConfigWithFlags overrides configuration with command-line flags
func overrideConfigWithFlags(config *EcosystemValidationConfig) {
	if *testDataPath != "./test_data" {
		config.TestDataPath = *testDataPath
	}
	
	if *referencePath != "./reference" {
		config.ReferencePath = *referencePath
	}
	
	if *toolsPath != "" {
		toolPairs := strings.Split(*toolsPath, ",")
		for _, pair := range toolPairs {
			parts := strings.SplitN(pair, ":", 2)
			if len(parts) == 2 {
				config.ToolPaths[parts[0]] = parts[1]
			}
		}
	}
	
	if *validationMode != "all" {
		modes := strings.Split(*validationMode, ",")
		config.ValidationModes = modes
	}
	
	if *skipCloud {
		config.CloudProvider = ""
	}
	
	// Set output path
	config.ReportPath = filepath.Join(*outputDir, "ecosystem_validation_report.json")
}

// validateConfig validates the configuration
func validateConfig(config *EcosystemValidationConfig) error {
	if config.TestDataPath == "" {
		return fmt.Errorf("test data path is required")
	}
	
	if config.ReferencePath == "" {
		return fmt.Errorf("reference path is required")
	}
	
	if len(config.ToolPaths) == 0 {
		return fmt.Errorf("at least one tool path must be specified")
	}
	
	if len(config.ValidationModes) == 0 {
		return fmt.Errorf("at least one validation mode must be specified")
	}
	
	// Validate validation modes
	validModes := map[string]bool{
		"basic": true, "performance": true, "population": true,
		"workflow": true, "realtime": true, "cloud": true, "all": true,
	}
	
	for _, mode := range config.ValidationModes {
		if !validModes[mode] {
			return fmt.Errorf("invalid validation mode: %s", mode)
		}
	}
	
	return nil
}

// printConfiguration prints the current configuration
func printConfiguration(config *EcosystemValidationConfig) {
	fmt.Println("\n📋 Configuration Summary:")
	fmt.Println("========================")
	fmt.Printf("Test Data Path: %s\n", config.TestDataPath)
	fmt.Printf("Reference Path: %s\n", config.ReferencePath)
	fmt.Printf("Validation Modes: %v\n", config.ValidationModes)
	fmt.Printf("Population Panels: %v\n", config.PopulationPanels)
	fmt.Printf("Platform Targets: %v\n", config.PlatformTargets)
	fmt.Printf("Max Concurrent Jobs: %d\n", config.MaxConcurrentJobs)
	fmt.Printf("Timeout Duration: %v\n", config.TimeoutDuration)
	
	fmt.Println("\n🔧 Tool Paths:")
	for tool, path := range config.ToolPaths {
		fmt.Printf("  %s: %s\n", tool, path)
	}
	
	fmt.Printf("\nReport Path: %s\n", config.ReportPath)
}

// executeEcosystemValidation executes the main ecosystem validation
func executeEcosystemValidation(config *EcosystemValidationConfig) error {
	if *verbose {
		log.Println("🚀 Starting ecosystem validation...")
		printConfiguration(config)
	}
	
	// Convert configuration to ecosystem config
	ecosystemConfig := convertToEcosystemConfig(config)
	
	// Create ecosystem validator
	validator, err := tools.NewGenomicsEcosystemValidator(ecosystemConfig)
	if err != nil {
		return fmt.Errorf("failed to create ecosystem validator: %w", err)
	}
	
	// Execute validation
	ctx, cancel := context.WithTimeout(context.Background(), config.TimeoutDuration)
	defer cancel()
	
	if *verbose {
		log.Println("🔍 Executing comprehensive ecosystem validation...")
	}
	
	report, err := validator.ValidateEcosystem(ctx)
	if err != nil {
		return fmt.Errorf("ecosystem validation failed: %w", err)
	}
	
	// Save detailed report
	if err := validator.SaveValidationReport(report, config.ReportPath); err != nil {
		return fmt.Errorf("failed to save validation report: %w", err)
	}
	
	// Generate summary report
	if err := generateSummaryReport(report, *outputDir); err != nil {
		return fmt.Errorf("failed to generate summary report: %w", err)
	}
	
	// Print results summary
	printValidationSummary(report)
	
	return nil
}

// convertToEcosystemConfig converts the CLI config to ecosystem config
func convertToEcosystemConfig(config *EcosystemValidationConfig) tools.EcosystemConfig {
	// Convert validation modes
	var validationModes []tools.ValidationMode
	for _, mode := range config.ValidationModes {
		switch mode {
		case "basic", "all":
			validationModes = append(validationModes, tools.ModeBasicCompatibility)
		case "performance":
			validationModes = append(validationModes, tools.ModePerformanceBenchmark)
		case "population":
			validationModes = append(validationModes, tools.ModePopulationScale)
		case "workflow":
			validationModes = append(validationModes, tools.ModeWorkflowIntegration)
		case "realtime":
			validationModes = append(validationModes, tools.ModeRealTimeValidation)
		case "cloud":
			validationModes = append(validationModes, tools.ModeCloudDeployment)
		}
	}
	
	// If "all" was specified, include all modes
	if contains(config.ValidationModes, "all") {
		validationModes = []tools.ValidationMode{
			tools.ModeBasicCompatibility,
			tools.ModePerformanceBenchmark,
			tools.ModePopulationScale,
			tools.ModeWorkflowIntegration,
			tools.ModeRealTimeValidation,
		}
		
		// Add cloud mode only if cloud provider is configured
		if config.CloudProvider != "" {
			validationModes = append(validationModes, tools.ModeCloudDeployment)
		}
	}
	
	return tools.EcosystemConfig{
		ToolPaths:         config.ToolPaths,
		TestDataPath:      config.TestDataPath,
		ReferencePath:     config.ReferencePath,
		ValidationModes:   validationModes,
		PopulationPanels:  config.PopulationPanels,
		PlatformTargets:   config.PlatformTargets,
		MaxConcurrentJobs: config.MaxConcurrentJobs,
		TimeoutDuration:   config.TimeoutDuration,
		ResourceLimits: tools.ResourceLimits{
			MaxMemoryGB:      config.MaxMemoryGB,
			MaxCPUCores:      config.MaxCPUCores,
			MaxDiskSpaceGB:   config.MaxDiskSpaceGB,
			MaxExecutionTime: config.MaxExecutionTime,
		},
		QualityThresholds: tools.QualityThresholds{
			MinAccuracy:       config.MinAccuracy,
			MaxFalsePositive:  config.MaxFalsePositive,
			MaxFalseNegative:  config.MaxFalseNegative,
			MinCoverage:       config.MinCoverage,
			MaxContamination:  config.MaxContamination,
			MinMappingQuality: config.MinMappingQuality,
		},
		CloudConfig: tools.CloudConfig{
			Provider:    config.CloudProvider,
			Region:      config.CloudRegion,
			Credentials: config.CloudCredentials,
		},
	}
}

// generateSummaryReport generates a human-readable summary report
func generateSummaryReport(report *tools.EcosystemValidationReport, outputDir string) error {
	summaryPath := filepath.Join(outputDir, "ecosystem_validation_summary.txt")
	
	file, err := os.Create(summaryPath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()
	
	fmt.Fprintf(file, "# Phantom Science Genomics Ecosystem Validation Summary\n\n")
	fmt.Fprintf(file, "## Execution Details\n")
	fmt.Fprintf(file, "- Start Time: %s\n", report.StartTime.Format(time.RFC3339))
	fmt.Fprintf(file, "- End Time: %s\n", report.EndTime.Format(time.RFC3339))
	fmt.Fprintf(file, "- Duration: %v\n", report.Duration)
	fmt.Fprintf(file, "- Overall Score: %.2f/1.0\n", report.OverallAssessment.OverallScore)
	fmt.Fprintf(file, "- Status: %s\n\n", report.OverallAssessment.Status)
	
	fmt.Fprintf(file, "## Tool Validation Results\n")
	for toolName, result := range report.ToolResults {
		fmt.Fprintf(file, "- %s: %s (v%s)\n", toolName, result.Status, result.Version)
		if result.ErrorMessage != "" {
			fmt.Fprintf(file, "  Error: %s\n", result.ErrorMessage)
		}
	}
	
	fmt.Fprintf(file, "\n## Performance Scores\n")
	fmt.Fprintf(file, "- Compatibility: %.2f/1.0\n", report.OverallAssessment.CompatibilityScore)
	fmt.Fprintf(file, "- Performance: %.2f/1.0\n", report.OverallAssessment.PerformanceScore)
	fmt.Fprintf(file, "- Scalability: %.2f/1.0\n", report.OverallAssessment.ScalabilityScore)
	fmt.Fprintf(file, "- Reliability: %.2f/1.0\n", report.OverallAssessment.ReliabilityScore)
	
	fmt.Fprintf(file, "\n## Recommendations\n")
	for _, recommendation := range report.OverallAssessment.Recommendations {
		fmt.Fprintf(file, "- %s\n", recommendation)
	}
	
	if len(report.OverallAssessment.CriticalIssues) > 0 {
		fmt.Fprintf(file, "\n## Critical Issues\n")
		for _, issue := range report.OverallAssessment.CriticalIssues {
			fmt.Fprintf(file, "- %s\n", issue)
		}
	}
	
	return nil
}

// printValidationSummary prints a summary of validation results
func printValidationSummary(report *tools.EcosystemValidationReport) {
	fmt.Println("\n🎯 Validation Summary:")
	fmt.Println("=====================")
	fmt.Printf("Overall Score: %.2f/1.0 (%s)\n", 
		report.OverallAssessment.OverallScore, 
		report.OverallAssessment.Status)
	fmt.Printf("Duration: %v\n", report.Duration)
	fmt.Printf("Tools Validated: %d\n", len(report.ToolResults))
	
	// Count passed/failed tools
	passed := 0
	failed := 0
	for _, result := range report.ToolResults {
		if result.Status == "passed" {
			passed++
		} else {
			failed++
		}
	}
	
	fmt.Printf("Tools Passed: %d\n", passed)
	fmt.Printf("Tools Failed: %d\n", failed)
	
	if len(report.OverallAssessment.Recommendations) > 0 {
		fmt.Println("\n📋 Key Recommendations:")
		for _, rec := range report.OverallAssessment.Recommendations {
			fmt.Printf("- %s\n", rec)
		}
	}
	
	fmt.Printf("\n📊 Full Report: %s\n", *outputDir)
}

// contains checks if a slice contains a string
func contains(slice []string, item string) bool {
	for _, s := range slice {
		if s == item {
			return true
		}
	}
	return false
}