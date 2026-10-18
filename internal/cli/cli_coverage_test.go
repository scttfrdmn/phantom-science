//go:build !release
// +build !release

package cli_test

import (
	"encoding/json"
	"os"
	"path/filepath"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestApp_GetConfig(t *testing.T) {
	tmpDir := t.TempDir()
	cfgPath := cliTestSeqJSON(t, tmpDir, "GetConfigModel", 150)
	app, err := cli.NewApp(cfgPath, "", false)
	if err != nil {
		t.Fatalf("NewApp failed: %v", err)
	}
	cfg := app.GetConfig()
	if cfg == nil {
		t.Fatal("GetConfig should not return nil")
	}
	if cfg.SequencerModel != "GetConfigModel" {
		t.Errorf("Expected GetConfigModel, got %s", cfg.SequencerModel)
	}
}

func TestApp_GetPhantomConfig(t *testing.T) {
	tmpDir := t.TempDir()
	cfgPath := cliTestSeqJSON(t, tmpDir, "PhantomModel", 150)
	app, err := cli.NewApp(cfgPath, "", false)
	if err != nil {
		t.Fatalf("NewApp failed: %v", err)
	}
	pcfg := app.GetPhantomConfig()
	if pcfg == nil {
		t.Fatal("GetPhantomConfig should not return nil")
	}
	if pcfg.Enterprise.Fingerprinting.Enable {
		t.Error("Expected fingerprinting disabled")
	}
}

func TestApp_SetConfig(t *testing.T) {
	tmpDir := t.TempDir()
	cfgPath := cliTestSeqJSON(t, tmpDir, "OldModel", 100)
	app, err := cli.NewApp(cfgPath, "", false)
	if err != nil {
		t.Fatalf("NewApp failed: %v", err)
	}
	app.SetConfig(types.SequencerConfig{
		SequencerModel: "UpdatedModel",
		ReadLength:     250,
		BaseQuality:    35,
	})
	cfg := app.GetConfig()
	if cfg == nil {
		t.Fatal("GetConfig nil after SetConfig")
	}
	if cfg.SequencerModel != "UpdatedModel" {
		t.Errorf("Expected UpdatedModel, got %s", cfg.SequencerModel)
	}
	if cfg.ReadLength != 250 {
		t.Errorf("Expected 250, got %d", cfg.ReadLength)
	}
}

func TestApp_SetPhantomConfig(t *testing.T) {
	tmpDir := t.TempDir()
	cfgPath := cliTestSeqJSON(t, tmpDir, "SetModel", 150)
	app, err := cli.NewApp(cfgPath, "", false)
	if err != nil {
		t.Fatalf("NewApp failed: %v", err)
	}
	app.SetPhantomConfig(&types.PhantomConfig{
		Enterprise: types.EnterpriseConfig{
			CustomerID:     "SET_CUST_001",
			Fingerprinting: types.FingerprintConfig{Enable: true},
		},
	})
	pcfg := app.GetPhantomConfig()
	if pcfg == nil {
		t.Fatal("GetPhantomConfig nil after SetPhantomConfig")
	}
	if pcfg.Enterprise.CustomerID != "SET_CUST_001" {
		t.Errorf("Expected SET_CUST_001, got %s", pcfg.Enterprise.CustomerID)
	}
	if !pcfg.Enterprise.Fingerprinting.Enable {
		t.Error("Expected fingerprinting enabled")
	}
}

func TestNewApp_VerboseMode(t *testing.T) {
	tmpDir := t.TempDir()
	cfgPath := cliTestSeqJSON(t, tmpDir, "VerboseModel", 75)
	app, err := cli.NewApp(cfgPath, "", true)
	if err != nil {
		t.Fatalf("NewApp failed: %v", err)
	}
	if app == nil {
		t.Fatal("Expected non-nil app in verbose mode")
	}
}

func TestConfigLoader_Verbose_SequencerJSON(t *testing.T) {
	data, _ := json.Marshal(map[string]interface{}{
		"sequencer_model": "VerboseMi",
		"sample_name":     "s",
		"lane":            1,
		"read_length":     300,
		"base_quality":    30,
	})
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "c.json")
	if err := os.WriteFile(p, data, 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(true)
	cfg, err := loader.LoadSequencerConfig(p)
	if err != nil {
		t.Fatalf("load failed: %v", err)
	}
	if cfg.SequencerModel != "VerboseMi" {
		t.Errorf("Expected VerboseMi, got %s", cfg.SequencerModel)
	}
}

func TestConfigLoader_Verbose_PhantomJSON(t *testing.T) {
	data, _ := json.Marshal(types.PhantomConfig{
		Enterprise: types.EnterpriseConfig{
			CustomerID:     "VJCUST",
			Fingerprinting: types.FingerprintConfig{Enable: true},
		},
		Validation: types.ValidationConfig{Enable: false},
	})
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "ph.json")
	if err := os.WriteFile(p, data, 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(true)
	cfg, err := loader.LoadPhantomConfig(p)
	if err != nil {
		t.Fatalf("load failed: %v", err)
	}
	if cfg.Enterprise.CustomerID != "VJCUST" {
		t.Errorf("Expected VJCUST, got %s", cfg.Enterprise.CustomerID)
	}
}

func TestConfigLoader_PhantomJSON_NonVerbose(t *testing.T) {
	data, _ := json.Marshal(types.PhantomConfig{
		Enterprise: types.EnterpriseConfig{
			CustomerID: "JSON001",
			Fingerprinting: types.FingerprintConfig{
				Enable:             false,
				DatasetFingerprint: true,
			},
		},
		Validation: types.ValidationConfig{Enable: true},
	})
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "ph.json")
	if err := os.WriteFile(p, data, 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(false)
	cfg, err := loader.LoadPhantomConfig(p)
	if err != nil {
		t.Fatalf("load failed: %v", err)
	}
	if cfg.Enterprise.CustomerID != "JSON001" {
		t.Errorf("Expected JSON001, got %s", cfg.Enterprise.CustomerID)
	}
	if cfg.Enterprise.Fingerprinting.Enable {
		t.Error("Expected fingerprinting disabled")
	}
	if !cfg.Enterprise.Fingerprinting.DatasetFingerprint {
		t.Error("Expected DatasetFingerprint enabled")
	}
	if !cfg.Validation.Enable {
		t.Error("Expected validation enabled")
	}
}


func TestConfigLoader_Verbose_SequencerYAML(t *testing.T) {
	yml := "schema_version: 1.0\n" +
		"config_type: sequencer\n" +
		"sequencer:\n" +
		"  model: VerboseNova\n" +
		"sample:\n" +
		"  name: s\n" +
		"  lane: 1\n" +
		"reads:\n" +
		"  length: 150\n" +
		"quality:\n" +
		"  base_quality: 35\n"
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "c.yaml")
	if err := os.WriteFile(p, []byte(yml), 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(true)
	cfg, err := loader.LoadSequencerConfig(p)
	if err != nil {
		t.Fatalf("load failed: %v", err)
	}
	if cfg.SequencerModel != "VerboseNova" {
		t.Errorf("Expected VerboseNova, got %s", cfg.SequencerModel)
	}
}

func TestConfigLoader_Verbose_PhantomYAML(t *testing.T) {
	yml := "schema_version: 1.0\n" +
		"config_type: phantom_enterprise\n" +
		"enterprise:\n" +
		"  customer_id: VCUST\n" +
		"  service_tier: enterprise\n" +
		"fingerprinting:\n" +
		"  enable: false\n" +
		"validation:\n" +
		"  enable: false\n"
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "ph.yaml")
	if err := os.WriteFile(p, []byte(yml), 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(true)
	cfg, err := loader.LoadPhantomConfig(p)
	if err != nil {
		t.Fatalf("load failed: %v", err)
	}
	if cfg.Enterprise.CustomerID != "VCUST" {
		t.Errorf("Expected VCUST, got %s", cfg.Enterprise.CustomerID)
	}
}

func TestConfigLoader_PhantomAutoDetect_YAML(t *testing.T) {
	yml := "schema_version: 1.0\n" +
		"config_type: phantom_enterprise\n" +
		"enterprise:\n" +
		"  customer_id: AUTO_Y\n" +
		"fingerprinting:\n" +
		"  enable: false\n" +
		"validation:\n" +
		"  enable: false\n"
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "ph.cfg")
	if err := os.WriteFile(p, []byte(yml), 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(false)
	cfg, err := loader.LoadPhantomConfig(p)
	if err != nil {
		t.Fatalf("auto-detect YAML failed: %v", err)
	}
	if cfg.Enterprise.CustomerID != "AUTO_Y" {
		t.Errorf("Expected AUTO_Y, got %s", cfg.Enterprise.CustomerID)
	}
}

func TestConfigLoader_PhantomLegacyYAML(t *testing.T) {
	// Legacy YAML path - yaml.v3 maps struct fields without yaml tags using lowercase
	// enterprise.Customer_id does not match CustomerID without yaml tags
	// Test that the legacy YAML path loads without error
	yml := "enterprise:\n" +
		"  fingerprinting:\n" +
		"    enable: true\n" +
		"validation:\n" +
		"  enable: false\n"
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "legacy.yaml")
	if err := os.WriteFile(p, []byte(yml), 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(false)
	cfg, err := loader.LoadPhantomConfig(p)
	if err != nil {
		t.Fatalf("legacy YAML failed: %v", err)
	}
	if cfg == nil {
		t.Fatal("Expected non-nil config from legacy YAML")
	}
}
func TestConfigLoader_PhantomAutoDetect_JSON(t *testing.T) {
	// Use .json extension to test loadPhantomJSON code path directly
	data, _ := json.Marshal(types.PhantomConfig{
		Enterprise: types.EnterpriseConfig{CustomerID: "AUTO_J"},
	})
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "ph.json")
	if err := os.WriteFile(p, data, 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(false)
	cfg, err := loader.LoadPhantomConfig(p)
	if err != nil {
		t.Fatalf("JSON phantom config load failed: %v", err)
	}
	if cfg.Enterprise.CustomerID != "AUTO_J" {
		t.Errorf("Expected AUTO_J, got %s", cfg.Enterprise.CustomerID)
	}
}
func TestConfigLoader_PhantomAutoDetect_Invalid(t *testing.T) {
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "ph.dat")
	if err := os.WriteFile(p, []byte("!!BAD@@##"), 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(false)
	_, err := loader.LoadPhantomConfig(p)
	if err == nil {
		t.Error("Expected error for invalid phantom config")
	}
}

func TestConfigLoader_SequencerAutoDetect_Invalid(t *testing.T) {
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "seq.dat")
	if err := os.WriteFile(p, []byte("!!BAD@@##"), 0600); err != nil {
		t.Fatal(err)
	}
	loader := cli.NewConfigLoader(false)
	_, err := loader.LoadSequencerConfig(p)
	if err == nil {
		t.Error("Expected error for invalid sequencer config")
	}
}

func TestAppGenerate_WithFingerprintingEnabled(t *testing.T) {
	tmpDir := t.TempDir()
	outputDir := filepath.Join(tmpDir, "out_fp")
	configPath := cliTestSeqJSONWithOutput(t, tmpDir, "MiSeq", 50, outputDir)
	phantomPath := cliTestPhantomConfig(t, tmpDir, "FP_CUST", true, false)
	app, err := cli.NewApp(configPath, phantomPath, false)
	if err != nil {
		t.Fatalf("NewApp failed: %v", err)
	}
	if err := app.Generate(); err != nil {
		t.Errorf("Generate with fingerprinting failed: %v", err)
	}
}

func TestAppGenerate_WithValidationEnabled(t *testing.T) {
	tmpDir := t.TempDir()
	outputDir := filepath.Join(tmpDir, "out_val")
	configPath := cliTestSeqJSONWithOutput(t, tmpDir, "NovaSeq 6000", 50, outputDir)
	phantomPath := cliTestPhantomConfig(t, tmpDir, "VAL_CUST", false, true)
	app, err := cli.NewApp(configPath, phantomPath, false)
	if err != nil {
		t.Fatalf("NewApp failed: %v", err)
	}
	if err := app.Generate(); err != nil {
		t.Errorf("Generate with validation failed: %v", err)
	}
}

func TestValidateConfigVersion_AllCases(t *testing.T) {
	loader := cli.NewConfigLoader(false)
	for _, v := range []string{"1.0"} {
		if err := loader.ValidateConfigVersion(v); err != nil {
			t.Errorf("Version %q should be valid: %v", v, err)
		}
	}
	for _, v := range []string{"0.1", "2.0", "1.1", "", "latest"} {
		if err := loader.ValidateConfigVersion(v); err == nil {
			t.Errorf("Version %q should be invalid", v)
		}
	}
}

func cliTestSeqJSON(t *testing.T, tmpDir, model string, readLength int) string {
	t.Helper()
	return cliTestSeqJSONWithOutput(t, tmpDir, model, readLength, "./test_output")
}

func cliTestSeqJSONWithOutput(t *testing.T, tmpDir, model string, readLength int, outputDir string) string {
	t.Helper()
	seqCfg := types.SequencerConfig{
		SequencerModel: model,
		SampleName:     "test",
		Lane:           1,
		OutputDir:      outputDir,
		ReadLength:     readLength,
		ReadsPerBatch:  100,
		BatchInterval:  "1s",
		TotalBatches:   1,
		PairedEnd:      false,
		BaseQuality:    30,
		QualityDecay:   0.1,
		Q30Percentage:  85.0,
	}
	data, err := json.Marshal(seqCfg)
	if err != nil {
		t.Fatalf("marshal seq config: %v", err)
	}
	p := filepath.Join(tmpDir, "sequencer.json")
	if err := os.WriteFile(p, data, 0600); err != nil {
		t.Fatalf("write seq config: %v", err)
	}
	return p
}

func cliTestPhantomConfig(t *testing.T, tmpDir, customerID string, enableFP, enableVal bool) string {
	t.Helper()
	phCfg := types.PhantomConfig{
		Enterprise: types.EnterpriseConfig{
			CustomerID: customerID,
			Fingerprinting: types.FingerprintConfig{
				Enable: enableFP,
			},
		},
		Validation: types.ValidationConfig{
			Enable: enableVal,
		},
	}
	data, err := json.Marshal(phCfg)
	if err != nil {
		t.Fatalf("marshal phantom config: %v", err)
	}
	p := filepath.Join(tmpDir, "phantom.json")
	if err := os.WriteFile(p, data, 0600); err != nil {
		t.Fatalf("write phantom config: %v", err)
	}
	return p
}
