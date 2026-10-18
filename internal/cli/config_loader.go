package cli

import (
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"reflect"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/security"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"gopkg.in/yaml.v3"
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
		return nil, fmt.Errorf("%s: %w", configPath, loadErr)
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

// errNoSequencerSettings reports a config file that parsed but set nothing.
var errNoSequencerSettings = errors.New("no recognised sequencer settings: keys must match the JSON field names " +
	"(sample_name, read_length, ...); v1 YAML needs schema_version, v2 platform templates need version: \"2.0\" and specifications")

// checkNonZero rejects a config with no settings at all, which is what an
// unrecognised layout decodes to.
func checkNonZero(config *types.SequencerConfig) (*types.SequencerConfig, error) {
	if reflect.ValueOf(*config).IsZero() {
		return nil, errNoSequencerSettings
	}
	return config, nil
}

// loadSequencerYAML loads YAML sequencer configuration. Three layouts are
// recognised: v2 platform templates (top-level version "2.x" and a
// specifications section, see loadSequencerTemplateV2), v1 configs
// (schema_version), and flat configs whose keys are the JSON field names
// of types.SequencerConfig. A file that sets no recognised setting is an
// error.
func (cl *ConfigLoader) loadSequencerYAML(content []byte) (*types.SequencerConfig, error) {
	var peek struct {
		Version        string    `yaml:"version"`
		SchemaVersion  string    `yaml:"schema_version"`
		Specifications yaml.Node `yaml:"specifications"`
	}
	if err := yaml.Unmarshal(content, &peek); err != nil {
		return nil, fmt.Errorf("failed to parse YAML config: %w", err)
	}
	if peek.SchemaVersion == "" && peek.Version != "" && peek.Specifications.Kind != 0 {
		return cl.loadSequencerTemplateV2(content, peek.Version)
	}

	if peek.SchemaVersion != "" {
		var configV1 types.SequencerConfigV1
		if err := yaml.Unmarshal(content, &configV1); err != nil {
			return nil, fmt.Errorf("failed to parse v%s YAML config: %w", peek.SchemaVersion, err)
		}
		if cl.verbose {
			fmt.Printf("Loaded v%s YAML sequencer config\n", configV1.SchemaVersion)
		}
		// Convert v1 to legacy format for compatibility
		return checkNonZero(configV1.ToLegacy())
	}

	// Flat YAML uses the same keys as JSON configs. (SequencerConfig has
	// JSON tags; decoding the YAML directly would silently ignore
	// snake_case keys such as sample_name.)
	data, err := yamlToJSON(content)
	if err != nil {
		return nil, fmt.Errorf("failed to parse YAML config: %w", err)
	}
	var config types.SequencerConfig
	if err := json.Unmarshal(data, &config); err != nil {
		return nil, fmt.Errorf("failed to parse YAML config: %w", err)
	}

	if cl.verbose {
		fmt.Println("Loaded legacy YAML sequencer config")
	}
	return checkNonZero(&config)
}

// yamlToJSON converts a YAML document to JSON, so it can be decoded with
// JSON struct tags.
func yamlToJSON(content []byte) ([]byte, error) {
	var doc any
	if err := yaml.Unmarshal(content, &doc); err != nil {
		return nil, err
	}
	v, err := jsonCompatible(doc)
	if err != nil {
		return nil, err
	}
	if v == nil {
		v = map[string]any{}
	}
	return json.Marshal(v)
}

// jsonCompatible converts YAML-decoded values to values encoding/json can
// marshal (string-keyed maps).
func jsonCompatible(v any) (any, error) {
	switch t := v.(type) {
	case map[string]any:
		for k, e := range t {
			c, err := jsonCompatible(e)
			if err != nil {
				return nil, err
			}
			t[k] = c
		}
		return t, nil
	case map[any]any:
		m := make(map[string]any, len(t))
		for k, e := range t {
			ks, ok := k.(string)
			if !ok {
				return nil, fmt.Errorf("non-string key %v", k)
			}
			c, err := jsonCompatible(e)
			if err != nil {
				return nil, err
			}
			m[ks] = c
		}
		return m, nil
	case []any:
		for i, e := range t {
			c, err := jsonCompatible(e)
			if err != nil {
				return nil, err
			}
			t[i] = c
		}
		return t, nil
	default:
		return v, nil
	}
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
	return checkNonZero(&config)
}

// loadSequencerAuto attempts to auto-detect format
func (cl *ConfigLoader) loadSequencerAuto(content []byte) (*types.SequencerConfig, error) {
	// Try YAML first (more common in bioinformatics)
	config, yamlErr := cl.loadSequencerYAML(content)
	if yamlErr == nil {
		return config, nil
	}

	// Fall back to JSON
	if config, err := cl.loadSequencerJSON(content); err == nil {
		return config, nil
	}

	return nil, fmt.Errorf("unable to parse config file as YAML or JSON: %w", yamlErr)
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

// convertPhantomV1ToLegacy converts a v1 phantom config to the legacy
// format. Settings with no legacy equivalent, such as the security,
// volume_testing, output_customization, performance_targets,
// testing_scenarios and simplified_quality sections, are dropped.
func (cl *ConfigLoader) convertPhantomV1ToLegacy(v1 *types.PhantomConfigV1) *types.PhantomConfig {
	return &types.PhantomConfig{
		Enterprise: types.EnterpriseConfig{
			CustomerID:      v1.Enterprise.CustomerID,
			LicenseKey:      v1.Enterprise.LicenseKey,
			ComplianceLevel: v1.Enterprise.ComplianceLevel,
			ServiceTier:     v1.Enterprise.ServiceTier,
			Fingerprinting: types.FingerprintConfig{
				Enable:               v1.Fingerprinting.Enable,
				CustomerFingerprint:  v1.Fingerprinting.CustomerFingerprint,
				DatasetFingerprint:   v1.Fingerprinting.DatasetFingerprint,
				SteganographyLayers:  v1.Fingerprinting.SteganographyLayers,
				KeyFile:              v1.Fingerprinting.KeyFile,
				Plaintext:            v1.Fingerprinting.Plaintext,
				CryptographicSigning: v1.Fingerprinting.CryptographicSigning,
				TamperDetection:      v1.Fingerprinting.TamperDetection,
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
// is not accessible. The engine falls back to uniform random bases regardless;
// this is an early pre-flight notice for the operator.
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

// Run settings applied to a v2 platform template when its run section does
// not set them. v2 files (configs/sequencers/*.yaml) describe an instrument,
// not a run, so they are completed with these documented defaults; a run
// section, or command-line flags, override them.
const (
	TemplateDefaultSampleName    = "sample"
	TemplateDefaultOutputDir     = "output"
	TemplateDefaultReadsPerBatch = 4000 // MinKNOW's default reads per FASTQ file
	TemplateDefaultTotalBatches  = 1
	TemplateDefaultTier          = "sup"
	TemplateDefaultBaseQuality   = 20
)

// sequencerTemplateV2 is the v2 platform template layout
// (configs/sequencers/*.yaml). Only the keys that drive generation are
// decoded; the rest is documentation.
type sequencerTemplateV2 struct {
	Version        string `yaml:"version"`
	Platform       string `yaml:"platform"`
	Vendor         string `yaml:"vendor"`
	Model          string `yaml:"model"`
	Chemistry      string `yaml:"chemistry"`
	Specifications struct {
		ReadType   string `yaml:"read_type"`
		ReadLength struct {
			Min          int    `yaml:"min"`
			Max          int    `yaml:"max"`
			N50Target    int    `yaml:"n50_target"`
			Median       int    `yaml:"median"`
			Distribution string `yaml:"distribution"`
		} `yaml:"read_length"`
		Accuracy struct {
			DefaultBasecaller string `yaml:"default_basecaller"`
			BasecallerModes   map[string]struct {
				Accuracy float64 `yaml:"accuracy"`
				Quality  int     `yaml:"quality"`
			} `yaml:"basecaller_modes"`
		} `yaml:"accuracy"`
		Throughput struct {
			MaxOutputTB           float64 `yaml:"max_output_tb"`
			ActiveChannelsPerFC   int     `yaml:"active_channels_per_fc"`
			TranslocationSpeedBPS float64 `yaml:"translocation_speed_bps"`
		} `yaml:"throughput"`
	} `yaml:"specifications"`
	Output struct {
		Compression string `yaml:"compression"`
		FilePattern string `yaml:"file_pattern"`
	} `yaml:"output"`
	// Run holds run settings with the same keys as a JSON sequencer config
	// (sample_name, output_dir, reads_per_batch, seed, ont, output, ...).
	Run yaml.Node `yaml:"run"`
}

// isONTTemplate reports whether a v2 template describes an Oxford Nanopore
// instrument.
func (t *sequencerTemplateV2) isONTTemplate() bool {
	return strings.HasPrefix(strings.ToLower(t.Platform), "ont") ||
		strings.Contains(strings.ToLower(t.Vendor), "oxford nanopore")
}

// loadSequencerTemplateV2 converts a v2 platform template to a
// SequencerConfig. Only Oxford Nanopore templates can be generated; other
// vendors' templates fail with an explicit error rather than loading as an
// empty config.
//
// Mapping: platform "ont"; specifications.read_length → the read-length
// distribution (log_normal → lognormal with n50_target, median, min, max;
// ReadLength is set to the median as the nominal length);
// specifications.accuracy.default_basecaller (default "sup") → error-model
// tier, and that mode's quality → BaseQuality (default 20, no decay);
// throughput.active_channels_per_fc and translocation_speed_bps → ONT
// channels and speed (package nanopore defaults when unset);
// output.compression and file_pattern → Output. Then the run section is
// applied, and run settings it leaves unset get the TemplateDefault*
// values.
func (cl *ConfigLoader) loadSequencerTemplateV2(content []byte, version string) (*types.SequencerConfig, error) {
	if !strings.HasPrefix(version, "2.") && version != "2" {
		return nil, fmt.Errorf("unsupported sequencer template version %q: want 2.x", version)
	}
	var t sequencerTemplateV2
	if err := yaml.Unmarshal(content, &t); err != nil {
		return nil, fmt.Errorf("failed to parse v%s sequencer template: %w", version, err)
	}
	if !t.isONTTemplate() {
		return nil, fmt.Errorf("v%s sequencer template for %s %s (platform %q) cannot be generated: only Oxford Nanopore templates are supported; use a v1 or JSON config for short-read runs",
			version, t.Vendor, t.Model, t.Platform)
	}
	spec := t.Specifications
	rl := spec.ReadLength
	c := types.SequencerConfig{
		Platform:         "ont",
		SequencerModel:   strings.TrimSpace(t.Vendor + " " + t.Model),
		ChemistryVersion: t.Chemistry,
		MaxOutputTB:      spec.Throughput.MaxOutputTB,
		ReadLength:       rl.Median,
		ONT: types.ONTSpec{
			Channels: spec.Throughput.ActiveChannelsPerFC,
			SpeedBPS: spec.Throughput.TranslocationSpeedBPS,
		},
		Output: types.OutputSpec{
			Compression: t.Output.Compression,
			FilePattern: t.Output.FilePattern,
		},
	}
	switch d := strings.ToLower(rl.Distribution); d {
	case "log_normal", "lognormal":
		if rl.N50Target <= 0 || rl.Median <= 0 {
			return nil, fmt.Errorf("v%s template read_length: log_normal needs n50_target and median", version)
		}
		c.ReadLengthDistribution = types.ReadLengthSpec{
			Dist: types.ReadLengthLognormal, N50: rl.N50Target, Median: rl.Median, Min: rl.Min, Max: rl.Max,
		}
	case "", "fixed":
		c.ReadLength = max(rl.Median, rl.N50Target)
	default:
		return nil, fmt.Errorf("v%s template read_length: unsupported distribution %q (want log_normal or fixed)", version, rl.Distribution)
	}
	tier := strings.ToLower(spec.Accuracy.DefaultBasecaller)
	if tier == "" {
		tier = TemplateDefaultTier
	}
	c.ErrorModel.Tier = tier
	c.BaseQuality = TemplateDefaultBaseQuality
	if m, ok := spec.Accuracy.BasecallerModes[tier]; ok && m.Quality > 0 {
		c.BaseQuality = m.Quality
	}

	if t.Run.Kind != 0 {
		var run any
		if err := t.Run.Decode(&run); err != nil {
			return nil, fmt.Errorf("v%s template run section: %w", version, err)
		}
		run, err := jsonCompatible(run)
		if err != nil {
			return nil, fmt.Errorf("v%s template run section: %w", version, err)
		}
		if _, ok := run.(map[string]any); !ok {
			return nil, fmt.Errorf("v%s template run section must be a mapping", version)
		}
		data, err := json.Marshal(run)
		if err != nil {
			return nil, fmt.Errorf("v%s template run section: %w", version, err)
		}
		dec := json.NewDecoder(strings.NewReader(string(data)))
		dec.DisallowUnknownFields() // a misspelt run key must not be ignored
		if err := dec.Decode(&c); err != nil {
			return nil, fmt.Errorf("v%s template run section: %w", version, err)
		}
	}
	var defaulted []string
	setDefault := func(name string, unset bool, apply func()) {
		if unset {
			apply()
			defaulted = append(defaulted, name)
		}
	}
	setDefault("sample_name="+TemplateDefaultSampleName, c.SampleName == "", func() { c.SampleName = TemplateDefaultSampleName })
	setDefault("output_dir="+TemplateDefaultOutputDir, c.OutputDir == "", func() { c.OutputDir = TemplateDefaultOutputDir })
	setDefault(fmt.Sprintf("reads_per_batch=%d", TemplateDefaultReadsPerBatch), c.ReadsPerBatch == 0, func() { c.ReadsPerBatch = TemplateDefaultReadsPerBatch })
	setDefault(fmt.Sprintf("total_batches=%d", TemplateDefaultTotalBatches), c.TotalBatches == 0, func() { c.TotalBatches = TemplateDefaultTotalBatches })
	if c.ReadLength <= 0 {
		return nil, fmt.Errorf("v%s template sets no read length: set specifications.read_length or run.read_length", version)
	}
	if cl.verbose {
		fmt.Printf("Loaded v%s sequencer template %s %s\n", version, t.Vendor, t.Model)
		if len(defaulted) > 0 {
			fmt.Printf("   run defaults applied: %s\n", strings.Join(defaulted, ", "))
		}
	}
	return &c, nil
}
