package cli

import (
	"os"
	"path/filepath"
	"testing"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestConfigLoaderYAML(t *testing.T) {
	// Create temporary YAML config
	yamlContent := `
schema_version: "1.0"
last_updated: "2024-07-03"
config_type: "sequencer"

sequencer:
  model: "NovaSeq 6000"
  flow_cell_type: "S2"
  chemistry_version: "v1.5"
  instrument_id: "A00123"
  run_id: "240101_A00123_0001_BH5TNVDRXY"
  flowcell_id: "H5TNVDRXY"

sample:
  name: "Sample_001"
  lane: 1
  output_dir: "./test_output"

reads:
  length: 150
  paired_end: true
  reads_per_batch: 1000000
  total_batches: 5
  batch_interval: "30m"

quality:
  base_quality: 35
  quality_decay: 0.05
  q30_percentage: 90.0

performance:
  run_time_hours: 36
  max_output_tb: 3.0

technology:
  sbs_version: "v1.5"
  four_channel_chemistry: true
`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "test_config.yaml")
	
	err := os.WriteFile(configPath, []byte(yamlContent), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("Failed to load YAML config: %v", err)
	}

	if config.SequencerModel != "NovaSeq 6000" {
		t.Errorf("Expected NovaSeq 6000, got %s", config.SequencerModel)
	}
	if config.ReadLength != 150 {
		t.Errorf("Expected read length 150, got %d", config.ReadLength)
	}
	if config.BaseQuality != 35 {
		t.Errorf("Expected base quality 35, got %d", config.BaseQuality)
	}
}

func TestConfigLoaderJSON(t *testing.T) {
	jsonContent := `{
  "sequencer_model": "MiSeq",
  "flow_cell_type": "v3",
  "sample_name": "Sample_001",
  "lane": 1,
  "output_dir": "./test_output",
  "read_length": 300,
  "reads_per_batch": 100000,
  "batch_interval": "15m",
  "total_batches": 3,
  "paired_end": true,
  "base_quality": 30,
  "quality_decay": 0.12,
  "q30_percentage": 80.0
}`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "test_config.json")
	
	err := os.WriteFile(configPath, []byte(jsonContent), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("Failed to load JSON config: %v", err)
	}

	if config.SequencerModel != "MiSeq" {
		t.Errorf("Expected MiSeq, got %s", config.SequencerModel)
	}
	if config.ReadLength != 300 {
		t.Errorf("Expected read length 300, got %d", config.ReadLength)
	}
	if config.Q30Percentage != 80.0 {
		t.Errorf("Expected Q30 percentage 80.0, got %f", config.Q30Percentage)
	}
}

func TestConfigLoaderPhantomYAML(t *testing.T) {
	yamlContent := `
schema_version: "1.0"
config_type: "phantom_enterprise"

enterprise:
  customer_id: "TEST_001"
  license_key: "test-license"
  service_tier: "enterprise"
  compliance_level: "hipaa_gdpr"

fingerprinting:
  enable: true
  customer_fingerprint: true
  dataset_fingerprint: true
  steganography_layers:
    - "quality_scores"
    - "read_ids"
  key_file: "/etc/phantom/fingerprint.key"
  plaintext: true
  cryptographic_signing: true
  tamper_detection: true

audit_logging:
  enable: true
  log_level: "info"
  retention_days: 2555
  compliance_format: "hipaa"

validation:
  enable: true
  validation_level: "regulatory"
  statistical_tests:
    - "hardy_weinberg_equilibrium"
    - "linkage_disequilibrium"
`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "phantom_config.yaml")
	
	err := os.WriteFile(configPath, []byte(yamlContent), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadPhantomConfig(configPath)
	if err != nil {
		t.Fatalf("Failed to load Phantom YAML config: %v", err)
	}

	if config.Enterprise.CustomerID != "TEST_001" {
		t.Errorf("Expected customer ID TEST_001, got %s", config.Enterprise.CustomerID)
	}
	if !config.Enterprise.Fingerprinting.Enable {
		t.Error("Expected fingerprinting to be enabled")
	}
	if len(config.Enterprise.Fingerprinting.SteganographyLayers) != 2 {
		t.Errorf("Expected 2 steganography layers, got %d", len(config.Enterprise.Fingerprinting.SteganographyLayers))
	}
	if got := config.Enterprise.Fingerprinting.KeyFile; got != "/etc/phantom/fingerprint.key" {
		t.Errorf("key_file not mapped from v1 YAML: got %q", got)
	}
	if !config.Enterprise.Fingerprinting.Plaintext {
		t.Error("plaintext not mapped from v1 YAML")
	}
}

func TestConfigLoaderAutoDetection(t *testing.T) {
	// Test auto-detection with .yaml extension
	yamlContent := `
schema_version: "1.0"
config_type: "sequencer"

sequencer:
  model: "NextSeq 2000"

sample:
  name: "test_sample"
  lane: 1

reads:
  length: 150

quality:
  base_quality: 33
`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "auto_config.yaml")
	
	err := os.WriteFile(configPath, []byte(yamlContent), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("Failed to load auto-detected config: %v", err)
	}

	if config.SequencerModel != "NextSeq 2000" {
		t.Errorf("Expected NextSeq 2000, got %s", config.SequencerModel)
	}
}

func TestConfigLoaderFileNotFound(t *testing.T) {
	loader := NewConfigLoader(false)
	_, err := loader.LoadSequencerConfig("/nonexistent/config.yaml")
	if err == nil {
		t.Error("Expected error for non-existent file")
	}
}

func TestConfigLoaderInvalidYAML(t *testing.T) {
	invalidYAML := `
sequencer:
  model: "NovaSeq"
  invalid_yaml: [unclosed bracket
`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "invalid.yaml")
	
	err := os.WriteFile(configPath, []byte(invalidYAML), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}

	loader := NewConfigLoader(false)
	_, err = loader.LoadSequencerConfig(configPath)
	if err == nil {
		t.Error("Expected error for invalid YAML")
	}
}

func TestConfigLoaderInvalidJSON(t *testing.T) {
	invalidJSON := `{
  "sequencer_model": "NovaSeq",
  "invalid": json
}`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "invalid.json")
	
	err := os.WriteFile(configPath, []byte(invalidJSON), 0600)
	if err != nil {
		t.Fatalf("Failed to create test config: %v", err)
	}

	loader := NewConfigLoader(false)
	_, err = loader.LoadSequencerConfig(configPath)
	if err == nil {
		t.Error("Expected error for invalid JSON")
	}
}

func TestValidateConfigVersion(t *testing.T) {
	loader := NewConfigLoader(false)
	
	// Test valid version
	err := loader.ValidateConfigVersion("1.0")
	if err != nil {
		t.Errorf("Expected version 1.0 to be valid, got error: %v", err)
	}
	
	// Test invalid version
	err = loader.ValidateConfigVersion("2.0")
	if err == nil {
		t.Error("Expected version 2.0 to be invalid")
	}
}

func TestConvertPhantomV1ToLegacy(t *testing.T) {
	loader := NewConfigLoader(false)
	
	v1Config := &types.PhantomConfigV1{
		Enterprise: types.EnterpriseConfigV1{
			CustomerID:      "TEST_123",
			LicenseKey:      "test-key",
			ComplianceLevel: "hipaa",
			ServiceTier:     "enterprise",
		},
		Fingerprinting: types.FingerprintConfigV1{
			Enable:               true,
			CustomerFingerprint:  true,
			DatasetFingerprint:   false,
			SteganographyLayers:  []string{"quality_scores"},
			CryptographicSigning: true,
			TamperDetection:      true,
		},
		Validation: types.ValidationConfigV1{
			Enable:          true,
			ValidationLevel: "academic",
			StatisticalTests: []string{"hardy_weinberg"},
		},
	}
	
	legacy := loader.convertPhantomV1ToLegacy(v1Config)
	
	if legacy.Enterprise.CustomerID != "TEST_123" {
		t.Errorf("Expected customer ID TEST_123, got %s", legacy.Enterprise.CustomerID)
	}
	if !legacy.Enterprise.Fingerprinting.Enable {
		t.Error("Expected fingerprinting to be enabled")
	}
	if legacy.Enterprise.Fingerprinting.DatasetFingerprint {
		t.Error("Expected dataset fingerprinting to be disabled")
	}
	if len(legacy.Enterprise.Fingerprinting.SteganographyLayers) != 1 {
		t.Errorf("Expected 1 steganography layer, got %d", len(legacy.Enterprise.Fingerprinting.SteganographyLayers))
	}
}

// TestConfigLoaderJSON_HMMModelPath verifies that hmm_model_path round-trips through
// JSON config loading.
func TestConfigLoaderJSON_HMMModelPath(t *testing.T) {
	jsonContent := `{
  "sequencer_model": "MiSeq",
  "sample_name": "HMM_Sample",
  "lane": 1,
  "output_dir": "./out",
  "read_length": 150,
  "reads_per_batch": 100,
  "total_batches": 1,
  "base_quality": 30,
  "quality_decay": 0.1,
  "hmm_model_path": "/data/models/yri.hmm.json"
}`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "seq.json")
	if err := os.WriteFile(configPath, []byte(jsonContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("LoadSequencerConfig: %v", err)
	}
	if config.HMMModelPath != "/data/models/yri.hmm.json" {
		t.Errorf("HMMModelPath: want /data/models/yri.hmm.json, got %q", config.HMMModelPath)
	}
}

// TestConfigLoaderLegacyYAML_HMMModelPath verifies that hmm_model_path round-trips
// through the legacy YAML loading path (yaml.Unmarshal directly into SequencerConfig).
func TestConfigLoaderLegacyYAML_HMMModelPath(t *testing.T) {
	yamlContent := `sequencer_model: NextSeq 2000
sample_name: HMM_Sample
lane: 1
output_dir: ./out
read_length: 150
reads_per_batch: 100
total_batches: 1
base_quality: 30
quality_decay: 0.1
hmm_model_path: /data/models/ceu.hmm.json
`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "legacy.yaml")
	if err := os.WriteFile(configPath, []byte(yamlContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("LoadSequencerConfig: %v", err)
	}
	if config.HMMModelPath != "/data/models/ceu.hmm.json" {
		t.Errorf("HMMModelPath: want /data/models/ceu.hmm.json, got %q", config.HMMModelPath)
	}
}

// TestConfigLoaderV1YAML_HMMModelPath verifies that hmm_model_path round-trips through
// the v1 YAML format (SequencerConfigV1 → ToLegacy).
func TestConfigLoaderV1YAML_HMMModelPath(t *testing.T) {
	yamlContent := `schema_version: "1.0"
config_type: sequencer
sequencer:
  model: NovaSeq 6000
sample:
  name: HMM_V1_Sample
  lane: 1
  output_dir: ./out
reads:
  length: 150
  reads_per_batch: 100
  total_batches: 1
quality:
  base_quality: 30
  quality_decay: 0.1
hmm_model_path: /data/models/chb.hmm.json
`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "v1.yaml")
	if err := os.WriteFile(configPath, []byte(yamlContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("LoadSequencerConfig: %v", err)
	}
	if config.HMMModelPath != "/data/models/chb.hmm.json" {
		t.Errorf("HMMModelPath: want /data/models/chb.hmm.json, got %q", config.HMMModelPath)
	}
}

// TestConfigLoader_HMMModelPath_EmptyWhenAbsent verifies that HMMModelPath is empty
// when not present in the config file (no regression on existing configs).
func TestConfigLoader_HMMModelPath_EmptyWhenAbsent(t *testing.T) {
	jsonContent := `{"sequencer_model":"MiSeq","sample_name":"s","lane":1,"read_length":150,"reads_per_batch":10,"total_batches":1,"base_quality":30,"quality_decay":0.1}`
	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "no_hmm.json")
	if err := os.WriteFile(configPath, []byte(jsonContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("LoadSequencerConfig: %v", err)
	}
	if config.HMMModelPath != "" {
		t.Errorf("expected empty HMMModelPath, got %q", config.HMMModelPath)
	}
}

// TestWarnIfHMMModelInvalid_InaccessiblePath verifies that warnIfHMMModelInvalid
// does not return an error — it only prints a warning — when the path is missing.
func TestWarnIfHMMModelInvalid_InaccessiblePath(t *testing.T) {
	loader := NewConfigLoader(false)
	// Must not panic; warning is printed to stdout (not captured here).
	loader.warnIfHMMModelInvalid("/nonexistent/model.hmm.json")
}

// TestWarnIfHMMModelInvalid_EmptyPath verifies no-op for empty path.
func TestWarnIfHMMModelInvalid_EmptyPath(t *testing.T) {
	loader := NewConfigLoader(false)
	loader.warnIfHMMModelInvalid("") // must not panic
}

// TestConfigLoaderJSON_Population verifies that population and model_dir round-trip
// through JSON config loading.
func TestConfigLoaderJSON_Population(t *testing.T) {
	jsonContent := `{
  "sequencer_model": "MiSeq",
  "sample_name": "Pop_Sample",
  "lane": 1,
  "output_dir": "./out",
  "read_length": 150,
  "reads_per_batch": 100,
  "total_batches": 1,
  "base_quality": 30,
  "quality_decay": 0.1,
  "population": "YRI",
  "model_dir": "/data/models"
}`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "seq.json")
	if err := os.WriteFile(configPath, []byte(jsonContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("LoadSequencerConfig: %v", err)
	}
	if config.Population != "YRI" {
		t.Errorf("Population: want YRI, got %q", config.Population)
	}
	if config.ModelDir != "/data/models" {
		t.Errorf("ModelDir: want /data/models, got %q", config.ModelDir)
	}
}

// TestConfigLoaderLegacyYAML_Population verifies that population and model_dir
// round-trip through the legacy YAML loading path.
func TestConfigLoaderLegacyYAML_Population(t *testing.T) {
	yamlContent := `sequencer_model: NextSeq 2000
sample_name: Pop_Sample
lane: 1
output_dir: ./out
read_length: 150
reads_per_batch: 100
total_batches: 1
base_quality: 30
quality_decay: 0.1
population: CEU
model_dir: /data/models
`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "legacy.yaml")
	if err := os.WriteFile(configPath, []byte(yamlContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("LoadSequencerConfig: %v", err)
	}
	if config.Population != "CEU" {
		t.Errorf("Population: want CEU, got %q", config.Population)
	}
	if config.ModelDir != "/data/models" {
		t.Errorf("ModelDir: want /data/models, got %q", config.ModelDir)
	}
}

// TestConfigLoaderV1YAML_Population verifies that population and model_dir
// round-trip through the v1 YAML format (SequencerConfigV1 → ToLegacy).
func TestConfigLoaderV1YAML_Population(t *testing.T) {
	yamlContent := `schema_version: "1.0"
config_type: sequencer
sequencer:
  model: NovaSeq 6000
sample:
  name: Pop_V1_Sample
  lane: 1
  output_dir: ./out
reads:
  length: 150
  reads_per_batch: 100
  total_batches: 1
quality:
  base_quality: 30
  quality_decay: 0.1
population: CHB
model_dir: /data/models
`

	tmpDir := t.TempDir()
	configPath := filepath.Join(tmpDir, "v1.yaml")
	if err := os.WriteFile(configPath, []byte(yamlContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	loader := NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		t.Fatalf("LoadSequencerConfig: %v", err)
	}
	if config.Population != "CHB" {
		t.Errorf("Population: want CHB, got %q", config.Population)
	}
	if config.ModelDir != "/data/models" {
		t.Errorf("ModelDir: want /data/models, got %q", config.ModelDir)
	}
}