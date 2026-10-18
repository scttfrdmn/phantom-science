package cli

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"gopkg.in/yaml.v3"
	"github.com/scttfrdmn/phantom-science/internal/security"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// ConfigLoader handles loading and parsing configuration files
type ConfigLoader struct {
	verbose bool
}

// NewConfigLoader creates a new configuration loader
func NewConfigLoader(verbose bool) *ConfigLoader {
	return &ConfigLoader{verbose: verbose}
}

// LoadSequencerConfig loads sequencer configuration from JSON or YAML
func (cl *ConfigLoader) LoadSequencerConfig(configPath string) (*types.SequencerConfig, error) {
	if cl.verbose {
		fmt.Printf("Loading sequencer config from: %s\n", configPath)
	}

	// Validate file path for security
	if err := security.ValidateFilePath(configPath, ""); err != nil {
		return nil, fmt.Errorf("invalid config file path: %w", err)
	}

	// Check if file exists
	if _, err := os.Stat(configPath); os.IsNotExist(err) {
		return nil, fmt.Errorf("config file not found: %s", configPath)
	}

	// Read file content
	content, err := os.ReadFile(configPath) // #nosec G304 - path validated by security.ValidateFilePath
	if err != nil {
		return nil, fmt.Errorf("failed to read config file: %w", err)
	}

	// Determine file format based on extension
	ext := strings.ToLower(filepath.Ext(configPath))
	
	var config *types.SequencerConfig
	var loadErr error
	switch ext {
	case ".yaml", ".yml":
		config, loadErr = cl.loadSequencerYAML(content)
	case ".json":
		config, loadErr = cl.loadSequencerJSON(content)
	default:
		config, loadErr = cl.loadSequencerAuto(content)
	}
	if loadErr != nil {
		return nil, loadErr
	}
	cl.warnIfHMMModelInvalid(config.HMMModelPath)
	// Also warn for the population-derived path when no explicit path overrides it.
	if config.HMMModelPath == "" && config.Population != "" && config.ModelDir != "" {
		cl.warnIfHMMModelInvalid(filepath.Join(config.ModelDir, config.Population+".hmm.json"))
	}
	return config, nil
}

// LoadPhantomConfig loads Phantom Science configuration from JSON or YAML
func (cl *ConfigLoader) LoadPhantomConfig(configPath string) (*types.PhantomConfig, error) {
	if cl.verbose {
		fmt.Printf("Loading phantom config from: %s\n", configPath)
	}

	// Validate file path for security
	if err := security.ValidateFilePath(configPath, ""); err != nil {
		return nil, fmt.Errorf("invalid config file path: %w", err)
	}

	// Check if file exists
	if _, err := os.Stat(configPath); os.IsNotExist(err) {
		return nil, fmt.Errorf("config file not found: %s", configPath)
	}

	// Read file content
	content, err := os.ReadFile(configPath) // #nosec G304 - path validated by security.ValidateFilePath
	if err != nil {
		return nil, fmt.Errorf("failed to read config file: %w", err)
	}

	// Determine file format based on extension
	ext := strings.ToLower(filepath.Ext(configPath))
	
	switch ext {
	case ".yaml", ".yml":
		return cl.loadPhantomYAML(content)
	case ".json":
		return cl.loadPhantomJSON(content)
	default:
		// Try to auto-detect format
		return cl.loadPhantomAuto(content)
	}
}

// loadSequencerYAML loads YAML sequencer configuration
func (cl *ConfigLoader) loadSequencerYAML(content []byte) (*types.SequencerConfig, error) {
	// Try to load as v1.0 YAML first
	var configV1 types.SequencerConfigV1
	if err := yaml.Unmarshal(content, &configV1); err == nil {
		// Check if it's a v1.0 config by looking for schema_version
		if configV1.SchemaVersion != "" {
			if cl.verbose {
				fmt.Printf("Loaded v%s YAML sequencer config\n", configV1.SchemaVersion)
			}
			// Convert v1 to legacy format for compatibility
			return configV1.ToLegacy(), nil
		}
	}

	// Fall back to legacy YAML format
	var config types.SequencerConfig
	if err := yaml.Unmarshal(content, &config); err != nil {
		return nil, fmt.Errorf("failed to parse YAML config: %w", err)
	}

	if cl.verbose {
		fmt.Println("Loaded legacy YAML sequencer config")
	}
	return &config, nil
}

// loadSequencerJSON loads JSON sequencer configuration
func (cl *ConfigLoader) loadSequencerJSON(content []byte) (*types.SequencerConfig, error) {
	var config types.SequencerConfig
	if err := json.Unmarshal(content, &config); err != nil {
		return nil, fmt.Errorf("failed to parse JSON config: %w", err)
	}

	if cl.verbose {
		fmt.Println("Loaded JSON sequencer config")
	}
	return &config, nil
}

// loadSequencerAuto attempts to auto-detect format
func (cl *ConfigLoader) loadSequencerAuto(content []byte) (*types.SequencerConfig, error) {
	// Try YAML first (more common in bioinformatics)
	if config, err := cl.loadSequencerYAML(content); err == nil {
		return config, nil
	}

	// Fall back to JSON
	if config, err := cl.loadSequencerJSON(content); err == nil {
		return config, nil
	}

	return nil, fmt.Errorf("unable to parse config file as YAML or JSON")
}

// loadPhantomYAML loads YAML Phantom Science configuration
func (cl *ConfigLoader) loadPhantomYAML(content []byte) (*types.PhantomConfig, error) {
	// Try to load as v1.0 YAML first
	var configV1 types.PhantomConfigV1
	if err := yaml.Unmarshal(content, &configV1); err == nil {
		// Check if it's a v1.0 config
		if configV1.SchemaVersion != "" {
			if cl.verbose {
				fmt.Printf("Loaded v%s YAML phantom config\n", configV1.SchemaVersion)
			}
			// Convert v1 to legacy format for compatibility
			return cl.convertPhantomV1ToLegacy(&configV1), nil
		}
	}

	// Fall back to legacy YAML format
	var config types.PhantomConfig
	if err := yaml.Unmarshal(content, &config); err != nil {
		return nil, fmt.Errorf("failed to parse YAML phantom config: %w", err)
	}

	if cl.verbose {
		fmt.Println("Loaded legacy YAML phantom config")
	}
	return &config, nil
}

// loadPhantomJSON loads JSON Phantom Science configuration
func (cl *ConfigLoader) loadPhantomJSON(content []byte) (*types.PhantomConfig, error) {
	var config types.PhantomConfig
	if err := json.Unmarshal(content, &config); err != nil {
		return nil, fmt.Errorf("failed to parse JSON phantom config: %w", err)
	}

	if cl.verbose {
		fmt.Println("Loaded JSON phantom config")
	}
	return &config, nil
}

// loadPhantomAuto attempts to auto-detect format
func (cl *ConfigLoader) loadPhantomAuto(content []byte) (*types.PhantomConfig, error) {
	// Try YAML first
	if config, err := cl.loadPhantomYAML(content); err == nil {
		return config, nil
	}

	// Fall back to JSON
	if config, err := cl.loadPhantomJSON(content); err == nil {
		return config, nil
	}

	return nil, fmt.Errorf("unable to parse phantom config file as YAML or JSON")
}

// convertPhantomV1ToLegacy converts v1 phantom config to legacy format
func (cl *ConfigLoader) convertPhantomV1ToLegacy(v1 *types.PhantomConfigV1) *types.PhantomConfig {
	return &types.PhantomConfig{
		Enterprise: types.EnterpriseConfig{
			CustomerID:      v1.Enterprise.CustomerID,
			LicenseKey:      v1.Enterprise.LicenseKey,
			ComplianceLevel: v1.Enterprise.ComplianceLevel,
			ServiceTier:     v1.Enterprise.ServiceTier,
			Fingerprinting: types.FingerprintConfig{
				Enable:              v1.Fingerprinting.Enable,
				CustomerFingerprint: v1.Fingerprinting.CustomerFingerprint,
				DatasetFingerprint:  v1.Fingerprinting.DatasetFingerprint,
				SteganographyLayers: v1.Fingerprinting.SteganographyLayers,
				KeyFile:             v1.Fingerprinting.KeyFile,
				Plaintext:           v1.Fingerprinting.Plaintext,
				CryptographicSigning: v1.Fingerprinting.CryptographicSigning,
				TamperDetection:     v1.Fingerprinting.TamperDetection,
			},
			AuditLogging: types.AuditConfig{
				Enable:           v1.AuditLogging.Enable,
				LogLevel:         v1.AuditLogging.LogLevel,
				RetentionDays:    v1.AuditLogging.RetentionDays,
				ComplianceFormat: v1.AuditLogging.ComplianceFormat,
				ExportEndpoint:   v1.AuditLogging.ExportEndpoint,
			},
		},
		Validation: types.ValidationConfig{
			Enable:                     v1.Validation.Enable,
			StatisticalTests:           v1.Validation.StatisticalTests,
			PopulationGeneticsTests:    v1.Validation.PopulationGeneticsTests,
			PipelineCompatibilityTests: v1.Validation.PipelineCompatibilityTests,
			RegulatoryCompliance: types.RegulatoryConfig{
				FDACompliance:     v1.RegulatoryCompliance.FDACompliance,
				EMACompliance:     v1.RegulatoryCompliance.EMACompliance,
				ValidationStudies: v1.RegulatoryCompliance.ValidationStudies,
				QualityMetrics:    v1.RegulatoryCompliance.QualityMetrics,
			},
			AcademicValidation: types.AcademicConfig{
				PeerReviewReady:   v1.AcademicValidation.PeerReviewReady,
				PublicationFormat: v1.AcademicValidation.PublicationFormat,
				BenchmarkDatasets: v1.AcademicValidation.BenchmarkDatasets,
				StatisticalPower:  v1.AcademicValidation.StatisticalPower,
			},
		},
	}
}

// warnIfHMMModelInvalid prints a warning when HMMModelPath is set but the file
// is not accessible. The engine will fall back to crypto/rand regardless; this
// is an early pre-flight notice for the operator.
func (cl *ConfigLoader) warnIfHMMModelInvalid(path string) {
	if path == "" {
		return
	}
	if _, err := os.Stat(path); err != nil {
		fmt.Printf("   [warn] hmm_model_path %q is not accessible: %v\n", path, err)
	}
}

// ValidateConfigVersion checks if the configuration version is supported
func (cl *ConfigLoader) ValidateConfigVersion(version string) error {
	supportedVersions := []string{"1.0"}
	
	for _, supported := range supportedVersions {
		if version == supported {
			return nil
		}
	}
	
	return fmt.Errorf("unsupported config version: %s, supported versions: %v", version, supportedVersions)
}