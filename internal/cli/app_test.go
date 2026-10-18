package cli

import (
	"os"
	"path/filepath"
	"testing"
)

func TestNewApp(t *testing.T) {
	// Create temporary config files
	tmpDir := t.TempDir()
	
	sequencerConfig := `{
  "sequencer_model": "MiSeq",
  "sample_name": "Test_Sample",
  "lane": 1,
  "output_dir": "./test_output",
  "read_length": 150,
  "reads_per_batch": 1000,
  "batch_interval": "10m",
  "total_batches": 2,
  "paired_end": true,
  "base_quality": 30,
  "quality_decay": 0.1,
  "q30_percentage": 85.0
}`

	configPath := filepath.Join(tmpDir, "sequencer.json")
	err := os.WriteFile(configPath, []byte(sequencerConfig), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}

	// Test app creation without phantom config
	app, err := NewApp(configPath, "", false)
	if err != nil {
		t.Fatalf("Expected app to be created, got error: %v", err)
	}
	if app == nil {
		t.Fatal("Expected app to be non-nil")
	}
	if app.config == nil {
		t.Fatal("Expected sequencer config to be loaded")
	}
	if app.phantomConfig == nil {
		t.Fatal("Expected phantom config to have defaults")
	}
	if app.phantomConfig.Enterprise.Fingerprinting.Enable {
		t.Error("Expected fingerprinting to be disabled by default")
	}
}

func TestNewAppWithPhantomConfig(t *testing.T) {
	tmpDir := t.TempDir()
	
	sequencerConfig := `{
  "sequencer_model": "NovaSeq 6000",
  "sample_name": "Test_Sample",
  "lane": 1,
  "output_dir": "./test_output",
  "read_length": 150,
  "reads_per_batch": 1000,
  "total_batches": 2,
  "paired_end": true,
  "base_quality": 35,
  "quality_decay": 0.05,
  "q30_percentage": 90.0
}`

	phantomConfig := `{
  "enterprise": {
    "customer_id": "TEST_001",
    "license_key": "test-key",
    "service_tier": "enterprise",
    "compliance_level": "hipaa",
    "fingerprinting": {
      "enable": true,
      "customer_fingerprint": true,
      "dataset_fingerprint": false,
      "steganography_layers": ["quality_scores"],
      "cryptographic_signing": true,
      "tamper_detection": true
    }
  },
  "validation": {
    "enable": true,
    "statistical_tests": ["hardy_weinberg"],
    "population_genetics_tests": [],
    "pipeline_compatibility_tests": []
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

	app, err := NewApp(sequencerPath, phantomPath, true)
	if err != nil {
		t.Fatalf("Expected app to be created, got error: %v", err)
	}
	if app == nil {
		t.Fatal("Expected app to be non-nil")
	}
	if !app.verbose {
		t.Error("Expected verbose mode to be enabled")
	}
	if app.phantomConfig.Enterprise.CustomerID != "TEST_001" {
		t.Errorf("Expected customer ID TEST_001, got %s", app.phantomConfig.Enterprise.CustomerID)
	}
	if !app.phantomConfig.Enterprise.Fingerprinting.Enable {
		t.Error("Expected fingerprinting to be enabled")
	}
}

func TestNewAppInvalidConfig(t *testing.T) {
	tmpDir := t.TempDir()

	invalidConfig := `{invalid json}`
	configPath := filepath.Join(tmpDir, "invalid.json")

	err := os.WriteFile(configPath, []byte(invalidConfig), 0600)
	if err != nil {
		t.Fatalf("Failed to create invalid config: %v", err)
	}

	_, err = NewApp(configPath, "", false)
	if err == nil {
		t.Error("Expected error for invalid config, got nil")
	}
}

func TestNewAppNonExistentConfig(t *testing.T) {
	_, err := NewApp("/nonexistent/config.json", "", false)
	if err == nil {
		t.Error("Expected error for non-existent config, got nil")
	}
}

func TestAppGenerate(t *testing.T) {
	tmpDir := t.TempDir()
	outputDir := filepath.Join(tmpDir, "output")
	
	sequencerConfig := `{
  "sequencer_model": "MiSeq",
  "sample_name": "Test_Sample",
  "lane": 1,
  "output_dir": "` + outputDir + `",
  "read_length": 50,
  "reads_per_batch": 100,
  "batch_interval": "1s",
  "total_batches": 1,
  "paired_end": false,
  "base_quality": 30,
  "quality_decay": 0.1,
  "q30_percentage": 85.0
}`

	configPath := filepath.Join(tmpDir, "sequencer.json")
	err := os.WriteFile(configPath, []byte(sequencerConfig), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}

	app, err := NewApp(configPath, "", false)
	if err != nil {
		t.Fatalf("NewApp failed: %v", err)
	}

	err = app.Generate()
	if err != nil {
		t.Errorf("Expected generation to succeed, got error: %v", err)
	}
	
	// Check if output directory was created
	if _, err := os.Stat(outputDir); os.IsNotExist(err) {
		t.Error("Expected output directory to be created")
	}
}

func TestAppConfigAccess(t *testing.T) {
	tmpDir := t.TempDir()
	
	sequencerConfig := `{
  "sequencer_model": "NextSeq 2000",
  "flow_cell_type": "P3",
  "sample_name": "Test_Sample",
  "lane": 1,
  "output_dir": "./test_output",
  "read_length": 75,
  "reads_per_batch": 5000,
  "total_batches": 3,
  "paired_end": true,
  "base_quality": 33,
  "quality_decay": 0.08,
  "q30_percentage": 90.0
}`

	configPath := filepath.Join(tmpDir, "sequencer.json")
	err := os.WriteFile(configPath, []byte(sequencerConfig), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}

	app, err := NewApp(configPath, "", false)
	if err != nil {
		t.Fatalf("NewApp failed: %v", err)
	}

	// Test that config is properly loaded and accessible
	if app.config.SequencerModel != "NextSeq 2000" {
		t.Errorf("Expected NextSeq 2000, got %s", app.config.SequencerModel)
	}
	if app.config.FlowCellType != "P3" {
		t.Errorf("Expected P3, got %s", app.config.FlowCellType)
	}
	if app.config.ReadLength != 75 {
		t.Errorf("Expected read length 75, got %d", app.config.ReadLength)
	}
	if app.config.ReadsPerBatch != 5000 {
		t.Errorf("Expected 5000 reads per batch, got %d", app.config.ReadsPerBatch)
	}
}