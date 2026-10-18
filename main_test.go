package main

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	phantomcmd "github.com/scttfrdmn/phantom-science/cmd"
)

func TestMainVersionFlag(t *testing.T) {
	// Build the binary first
	tmpDir := t.TempDir()
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err := cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test version flag
	cmd = exec.Command(binaryPath, "--version")
	output, err := cmd.Output()
	if err != nil {
		t.Fatalf("Failed to run version command: %v", err)
	}
	
	outputStr := string(output)
	want := "Phantom Science v" + phantomcmd.Version
	if !strings.Contains(outputStr, want) {
		t.Errorf("Expected version output to contain %q, got: %s", want, outputStr)
	}
	
	if !strings.Contains(outputStr, "Scientific Simulation Platform") {
		t.Errorf("Expected version output to contain platform description, got: %s", outputStr)
	}
	
	if !strings.Contains(outputStr, "Copyright 2025-2026 Scott Friedman") {
		t.Errorf("Expected version output to contain copyright, got: %s", outputStr)
	}
}

func TestMainNoConfig(t *testing.T) {
	tmpDir := t.TempDir()
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err := cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test without config flag - should fail
	cmd = exec.Command(binaryPath)
	output, err := cmd.CombinedOutput()
	if err == nil {
		t.Error("Expected command to fail without config flag")
	}
	
	outputStr := string(output)
	if !strings.Contains(outputStr, "Sequencer configuration file is required") {
		t.Errorf("Expected error message about required config, got: %s", outputStr)
	}
}

func TestMainWithValidConfig(t *testing.T) {
	tmpDir := t.TempDir()
	
	// Create test config
	configContent := `{
  "sequencer_model": "MiSeq",
  "sample_name": "Test_Sample",
  "lane": 1,
  "output_dir": "` + filepath.Join(tmpDir, "output") + `",
  "read_length": 50,
  "reads_per_batch": 10,
  "batch_interval": "1s",
  "total_batches": 1,
  "paired_end": false,
  "base_quality": 30,
  "quality_decay": 0.1,
  "q30_percentage": 85.0
}`
	
	configPath := filepath.Join(tmpDir, "test_config.json")
	err := os.WriteFile(configPath, []byte(configContent), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}
	
	// Build binary
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err = cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test with valid config
	cmd = exec.Command(binaryPath, "--config", configPath, "--command", "generate")
	output, err := cmd.CombinedOutput()
	if err != nil {
		t.Logf("Command output: %s", string(output))
		t.Fatalf("Expected command to succeed with valid config: %v", err)
	}
	
	outputStr := string(output)
	if !strings.Contains(outputStr, "Starting Phantom Science genomic data generation") {
		t.Errorf("Expected generation start message, got: %s", outputStr)
	}
	
	if !strings.Contains(outputStr, "Generation completed successfully") {
		t.Errorf("Expected success message, got: %s", outputStr)
	}
}

func TestMainWithPhantomConfig(t *testing.T) {
	tmpDir := t.TempDir()
	
	// Create test configs
	sequencerConfig := `{
  "sequencer_model": "NovaSeq 6000",
  "sample_name": "Test_Sample",
  "lane": 1,
  "output_dir": "` + filepath.Join(tmpDir, "output") + `",
  "read_length": 50,
  "reads_per_batch": 5,
  "batch_interval": "1s",
  "total_batches": 1,
  "paired_end": false,
  "base_quality": 35,
  "quality_decay": 0.05,
  "q30_percentage": 90.0
}`

	phantomConfig := `{
  "enterprise": {
    "customer_id": "TEST_001",
    "fingerprinting": {
      "enable": false
    }
  },
  "validation": {
    "enable": false
  }
}`
	
	sequencerPath := filepath.Join(tmpDir, "sequencer.json")
	phantomPath := filepath.Join(tmpDir, "phantom.json")
	
	err := os.WriteFile(sequencerPath, []byte(sequencerConfig), 0600)
	if err != nil {
		t.Fatalf("Failed to create sequencer config: %v", err)
	}
	
	err = os.WriteFile(phantomPath, []byte(phantomConfig), 0600)
	if err != nil {
		t.Fatalf("Failed to create phantom config: %v", err)
	}
	
	// Build binary
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err = cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test with both configs
	cmd = exec.Command(binaryPath, "--config", sequencerPath, "--phantom", phantomPath, "--command", "generate")
	output, err := cmd.CombinedOutput()
	if err != nil {
		t.Logf("Command output: %s", string(output))
		t.Fatalf("Expected command to succeed with both configs: %v", err)
	}
	
	outputStr := string(output)
	if !strings.Contains(outputStr, "Generation completed successfully") {
		t.Errorf("Expected success message with phantom config, got: %s", outputStr)
	}
}

func TestMainValidateCommand(t *testing.T) {
	tmpDir := t.TempDir()
	outDir := filepath.Join(tmpDir, "out")

	sequencerConfig := fmt.Sprintf(`{
  "sequencer_model": "MiSeq",
  "sample_name": "Test_Sample",
  "lane": 1,
  "output_dir": %q,
  "read_length": 150,
  "reads_per_batch": 100,
  "total_batches": 1,
  "paired_end": true,
  "base_quality": 30,
  "quality_decay": 0.1
}`, outDir)

	phantomConfig := `{
  "enterprise": {"customer_id": "VALIDATE_TEST", "fingerprinting": {"enable": false}},
  "validation": {
    "enable": true,
    "statistical_tests": ["format_compliance", "read_length_consistency"],
    "pipeline_compatibility_tests": ["fastq_parsing", "quality_encoding"]
  }
}`

	sequencerPath := filepath.Join(tmpDir, "sequencer.json")
	phantomPath := filepath.Join(tmpDir, "phantom.json")
	if err := os.WriteFile(sequencerPath, []byte(sequencerConfig), 0600); err != nil {
		t.Fatalf("Failed to create sequencer config: %v", err)
	}
	if err := os.WriteFile(phantomPath, []byte(phantomConfig), 0600); err != nil {
		t.Fatalf("Failed to create phantom config: %v", err)
	}

	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	if err := exec.Command("go", "build", "-o", binaryPath, ".").Run(); err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}

	if out, err := exec.Command(binaryPath, "generate", "--config", sequencerPath).CombinedOutput(); err != nil {
		t.Fatalf("generate failed: %v\n%s", err, out)
	}
	files, _ := filepath.Glob(filepath.Join(outDir, "*.fastq"))
	if len(files) != 2 {
		t.Fatalf("want 2 generated files, got %v", files)
	}

	// validate checks every file with the phantom config's tests.
	out, err := exec.Command(binaryPath, append([]string{"validate", "--phantom", phantomPath}, files...)...).CombinedOutput()
	if err != nil {
		t.Fatalf("validate failed on generated files: %v\n%s", err, out)
	}
	if got := strings.Count(string(out), "All 4 validation tests passed"); got != 2 {
		t.Errorf("want both files to pass 4 tests, report:\n%s", out)
	}

	// The legacy --command validate cannot take files and must not report success.
	out, err = exec.Command(binaryPath, "--config", sequencerPath, "--phantom", phantomPath, "--command", "validate").CombinedOutput()
	if err == nil {
		t.Fatalf("--command validate without files succeeded:\n%s", out)
	}
	if !strings.Contains(string(out), "usage: phantom-science validate FILE") {
		t.Errorf("missing usage hint:\n%s", out)
	}
}

func TestMainUnknownCommand(t *testing.T) {
	tmpDir := t.TempDir()
	
	configContent := `{
  "sequencer_model": "MiSeq",
  "sample_name": "Test_Sample",
  "lane": 1,
  "output_dir": "./test_output",
  "read_length": 50,
  "reads_per_batch": 10,
  "total_batches": 1,
  "paired_end": false,
  "base_quality": 30,
  "quality_decay": 0.1,
  "q30_percentage": 85.0
}`
	
	configPath := filepath.Join(tmpDir, "test_config.json")
	err := os.WriteFile(configPath, []byte(configContent), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}
	
	// Build binary
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err = cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test with unknown command
	cmd = exec.Command(binaryPath, "--config", configPath, "--command", "unknown")
	output, err := cmd.CombinedOutput()
	if err == nil {
		t.Error("Expected command to fail with unknown command")
	}
	
	outputStr := string(output)
	if !strings.Contains(outputStr, "Unknown command: unknown") {
		t.Errorf("Expected unknown command error, got: %s", outputStr)
	}
}

func TestMainVerboseFlag(t *testing.T) {
	tmpDir := t.TempDir()
	
	configContent := `{
  "sequencer_model": "MiSeq",
  "sample_name": "Test_Sample",
  "lane": 1,
  "output_dir": "` + filepath.Join(tmpDir, "output") + `",
  "read_length": 50,
  "reads_per_batch": 5,
  "batch_interval": "1s",
  "total_batches": 1,
  "paired_end": false,
  "base_quality": 30,
  "quality_decay": 0.1,
  "q30_percentage": 85.0
}`
	
	configPath := filepath.Join(tmpDir, "test_config.json")
	err := os.WriteFile(configPath, []byte(configContent), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}
	
	// Build binary
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err = cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test with verbose flag
	cmd = exec.Command(binaryPath, "--config", configPath, "--verbose", "--command", "generate")
	output, err := cmd.CombinedOutput()
	if err != nil {
		t.Logf("Command output: %s", string(output))
		t.Fatalf("Expected command to succeed with verbose flag: %v", err)
	}
	
	// Verbose mode should produce more detailed output
	outputStr := string(output)
	if !strings.Contains(outputStr, "Phantom Science genomic data generation") {
		t.Errorf("Expected detailed output in verbose mode, got: %s", outputStr)
	}
}

func TestMainHelpFlag(t *testing.T) {
	// Build the binary first
	tmpDir := t.TempDir()
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err := cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test help flag
	cmd = exec.Command(binaryPath, "--help")
	output, err := cmd.Output()
	if err != nil {
		t.Fatalf("Failed to run help command: %v", err)
	}
	
	outputStr := string(output)
	if !strings.Contains(outputStr, "enterprise-grade scientific simulation platform") {
		t.Errorf("Expected help output to contain platform description, got: %s", outputStr)
	}
	
	if !strings.Contains(outputStr, "Usage:") {
		t.Errorf("Expected help output to contain usage information, got: %s", outputStr)
	}
	
	if !strings.Contains(outputStr, "phantom-science") {
		t.Errorf("Expected help output to contain command name, got: %s", outputStr)
	}
}

func TestMainConfigValidation(t *testing.T) {
	// Build the binary first
	tmpDir := t.TempDir()
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err := cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test with missing config - should fail
	cmd = exec.Command(binaryPath)
	output, err := cmd.CombinedOutput()
	if err == nil {
		t.Error("Expected command to fail without config file")
	}
	
	outputStr := string(output)
	if !strings.Contains(outputStr, "configuration file is required") {
		t.Errorf("Expected error about missing config, got: %s", outputStr)
	}
}

func TestMainErrorHandling(t *testing.T) {
	// Build the binary first
	tmpDir := t.TempDir()
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err := cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test with invalid config file
	invalidConfigPath := filepath.Join(tmpDir, "invalid.json")
	err = os.WriteFile(invalidConfigPath, []byte("{invalid json}"), 0600)
	if err != nil {
		t.Fatalf("Failed to create invalid config: %v", err)
	}
	
	cmd = exec.Command(binaryPath, "--config", invalidConfigPath)
	output, err := cmd.CombinedOutput()
	if err == nil {
		t.Error("Expected command to fail with invalid config")
	}
	
	// Should get an error about invalid JSON or config parsing
	outputStr := string(output)
	if !strings.Contains(outputStr, "invalid") && !strings.Contains(outputStr, "error") && !strings.Contains(outputStr, "failed") {
		t.Errorf("Expected error message for invalid config, got: %s", outputStr)
	}
}

func TestMainUnknownFlag(t *testing.T) {
	// Build the binary first
	tmpDir := t.TempDir()
	binaryPath := filepath.Join(tmpDir, "phantom-science-test")
	
	cmd := exec.Command("go", "build", "-o", binaryPath, ".")
	err := cmd.Run()
	if err != nil {
		t.Fatalf("Failed to build binary: %v", err)
	}
	
	// Test with unknown flag
	cmd = exec.Command(binaryPath, "--unknown-flag")
	output, err := cmd.CombinedOutput()
	if err == nil {
		t.Error("Expected command to fail with unknown flag")
	}
	
	outputStr := string(output)
	if !strings.Contains(outputStr, "unknown") && !strings.Contains(outputStr, "flag") {
		t.Errorf("Expected error about unknown flag, got: %s", outputStr)
	}
}