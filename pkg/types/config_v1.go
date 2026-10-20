package types

import "time"

// ConfigVersion represents the configuration schema version
type ConfigVersion struct {
	SchemaVersion    string `yaml:"schema_version" json:"schema_version"`
	LastUpdated      string `yaml:"last_updated" json:"last_updated"`
	ValidatedAgainst string `yaml:"validated_against,omitempty" json:"validated_against,omitempty"`
	ConfigType       string `yaml:"config_type" json:"config_type"`
}

// SequencerConfigV1 represents versioned sequencer configuration
type SequencerConfigV1 struct {
	ConfigVersion `yaml:",inline"`

	// Sequencer Hardware
	Sequencer SequencerHardware `yaml:"sequencer" json:"sequencer"`

	// Sample Configuration
	Sample SampleConfig `yaml:"sample" json:"sample"`

	// Read Structure
	Reads ReadConfig `yaml:"reads" json:"reads"`

	// Quality Parameters
	Quality QualityConfig `yaml:"quality" json:"quality"`

	// Performance Specifications
	Performance PerformanceConfig `yaml:"performance" json:"performance"`

	// Technology Features
	Technology TechnologyConfig `yaml:"technology" json:"technology"`

	// Application Information (optional)
	Applications ApplicationConfig `yaml:"applications,omitempty" json:"applications,omitempty"`

	// HMMModelPath is an optional path to a serialized HMM model (.hmm.json).
	// When set, the engine uses it for biologically authentic sequence generation.
	HMMModelPath string `yaml:"hmm_model_path,omitempty" json:"hmm_model_path,omitempty"`

	// Population is a 1000 Genomes population code (e.g. "YRI", "CEU", "CHB").
	// Used with ModelDir to auto-resolve a model file.
	Population string `yaml:"population,omitempty" json:"population,omitempty"`

	// ModelDir is a directory containing <Population>.hmm.json model files.
	ModelDir string `yaml:"model_dir,omitempty" json:"model_dir,omitempty"`

	// Illumina is the Illumina read model (qualities, errors, insert size),
	// with the same keys as the illumina section of JSON configs.
	Illumina IlluminaSpec `yaml:"illumina,omitempty" json:"illumina,omitzero"`
}

type SequencerHardware struct {
	Model            string `yaml:"model" json:"model"`
	FlowCellType     string `yaml:"flow_cell_type" json:"flow_cell_type"`
	ChemistryVersion string `yaml:"chemistry_version,omitempty" json:"chemistry_version,omitempty"`
	InstrumentID     string `yaml:"instrument_id" json:"instrument_id"`
	RunID            string `yaml:"run_id" json:"run_id"`
	FlowCellID       string `yaml:"flowcell_id" json:"flowcell_id"`
}

type SampleConfig struct {
	Name      string `yaml:"name" json:"name"`
	Lane      int    `yaml:"lane" json:"lane"`
	OutputDir string `yaml:"output_dir" json:"output_dir"`
}

type ReadConfig struct {
	Length        int    `yaml:"length" json:"length"`
	PairedEnd     bool   `yaml:"paired_end" json:"paired_end"`
	ReadsPerBatch int    `yaml:"reads_per_batch" json:"reads_per_batch"`
	TotalBatches  int    `yaml:"total_batches" json:"total_batches"`
	BatchInterval string `yaml:"batch_interval" json:"batch_interval"`
}

type QualityConfig struct {
	BaseQuality           int     `yaml:"base_quality" json:"base_quality"`
	QualityDecay          float64 `yaml:"quality_decay" json:"quality_decay"`
	Q30Percentage         float64 `yaml:"q30_percentage" json:"q30_percentage"`
	ClusterDensityKPerMM2 int     `yaml:"cluster_density_k_per_mm2,omitempty" json:"cluster_density_k_per_mm2,omitempty"`
}

type PerformanceConfig struct {
	RunTimeHours           int     `yaml:"run_time_hours,omitempty" json:"run_time_hours,omitempty"`
	MaxOutputTB            float64 `yaml:"max_output_tb,omitempty" json:"max_output_tb,omitempty"`
	MaxOutputGB            float64 `yaml:"max_output_gb,omitempty" json:"max_output_gb,omitempty"`
	MaxCycles              int     `yaml:"max_cycles,omitempty" json:"max_cycles,omitempty"`
	LoadingConcentrationPM float64 `yaml:"loading_concentration_pm,omitempty" json:"loading_concentration_pm,omitempty"`
	PhixSpikeInPercent     float64 `yaml:"phix_spike_in_percent,omitempty" json:"phix_spike_in_percent,omitempty"`
}

type TechnologyConfig struct {
	SBSVersion            string   `yaml:"sbs_version,omitempty" json:"sbs_version,omitempty"`
	PatternedFlowCell     bool     `yaml:"patterned_flow_cell,omitempty" json:"patterned_flow_cell,omitempty"`
	DualFlowCellCapable   bool     `yaml:"dual_flow_cell_capable,omitempty" json:"dual_flow_cell_capable,omitempty"`
	TwoChannelChemistry   bool     `yaml:"two_channel_chemistry,omitempty" json:"two_channel_chemistry,omitempty"`
	FourChannelChemistry  bool     `yaml:"four_channel_chemistry,omitempty" json:"four_channel_chemistry,omitempty"`
	LowDiversitySensitive bool     `yaml:"low_diversity_sensitive,omitempty" json:"low_diversity_sensitive,omitempty"`
	BenchtopFormat        bool     `yaml:"benchtop_format,omitempty" json:"benchtop_format,omitempty"`
	XLEAPImprovements     []string `yaml:"xleap_improvements,omitempty" json:"xleap_improvements,omitempty"`
	ReagentImprovements   []string `yaml:"reagent_improvements,omitempty" json:"reagent_improvements,omitempty"`
	ClusterGeneration     string   `yaml:"cluster_generation,omitempty" json:"cluster_generation,omitempty"`
	FlowCellDesign        string   `yaml:"flow_cell_design,omitempty" json:"flow_cell_design,omitempty"`
}

type ApplicationConfig struct {
	OptimalFor          []string `yaml:"optimal_for,omitempty" json:"optimal_for,omitempty"`
	ReadLengthAdvantage string   `yaml:"read_length_advantage,omitempty" json:"read_length_advantage,omitempty"`
}

// PhantomConfigV1 represents versioned Phantom Science configuration
type PhantomConfigV1 struct {
	ConfigVersion `yaml:",inline"`

	// Enterprise Features
	Enterprise EnterpriseConfigV1 `yaml:"enterprise" json:"enterprise"`

	// Fingerprinting
	Fingerprinting FingerprintConfigV1 `yaml:"fingerprinting" json:"fingerprinting"`

	// Audit Logging
	AuditLogging AuditConfigV1 `yaml:"audit_logging" json:"audit_logging"`

	// Validation Framework
	Validation ValidationConfigV1 `yaml:"validation" json:"validation"`

	// Regulatory Compliance
	RegulatoryCompliance RegulatoryConfigV1 `yaml:"regulatory_compliance,omitempty" json:"regulatory_compliance,omitempty"`

	// Academic Validation
	AcademicValidation AcademicConfigV1 `yaml:"academic_validation,omitempty" json:"academic_validation,omitempty"`

	// Security Features
	Security SecurityConfigV1 `yaml:"security,omitempty" json:"security,omitempty"`

	// Volume Testing (optional)
	VolumeTesting VolumeTestingConfig `yaml:"volume_testing,omitempty" json:"volume_testing,omitempty"`

	// Output Customization
	OutputCustomization OutputCustomizationConfig `yaml:"output_customization,omitempty" json:"output_customization,omitempty"`

	// Performance Targets
	PerformanceTargets PerformanceTargetsConfig `yaml:"performance_targets,omitempty" json:"performance_targets,omitempty"`

	// Testing Scenarios
	TestingScenarios map[string]TestingScenario `yaml:"testing_scenarios,omitempty" json:"testing_scenarios,omitempty"`

	// Simplified Quality (for volume testing)
	SimplifiedQuality SimplifiedQualityConfig `yaml:"simplified_quality,omitempty" json:"simplified_quality,omitempty"`
}

type EnterpriseConfigV1 struct {
	CustomerID      string `yaml:"customer_id" json:"customer_id"`
	LicenseKey      string `json:"-" yaml:"-"`
	Organization    string `yaml:"organization,omitempty" json:"organization,omitempty"`
	ServiceTier     string `yaml:"service_tier" json:"service_tier"`
	ComplianceLevel string `yaml:"compliance_level" json:"compliance_level"`
}

type FingerprintConfigV1 struct {
	Enable               bool     `yaml:"enable" json:"enable"`
	CustomerFingerprint  bool     `yaml:"customer_fingerprint" json:"customer_fingerprint"`
	DatasetFingerprint   bool     `yaml:"dataset_fingerprint" json:"dataset_fingerprint"`
	SteganographyLayers  []string `yaml:"steganography_layers" json:"steganography_layers"`
	KeyFile              string   `yaml:"key_file,omitempty" json:"key_file,omitempty"`
	Plaintext            bool     `yaml:"plaintext,omitempty" json:"plaintext,omitempty"`
	CryptographicSigning bool     `yaml:"cryptographic_signing" json:"cryptographic_signing"`
	TamperDetection      bool     `yaml:"tamper_detection" json:"tamper_detection"`
}

type AuditConfigV1 struct {
	Enable           bool   `yaml:"enable" json:"enable"`
	LogLevel         string `yaml:"log_level" json:"log_level"`
	RetentionDays    int    `yaml:"retention_days" json:"retention_days"`
	ComplianceFormat string `yaml:"compliance_format" json:"compliance_format"`
	ExportEndpoint   string `yaml:"export_endpoint,omitempty" json:"export_endpoint,omitempty"`
	Encryption       string `yaml:"encryption,omitempty" json:"encryption,omitempty"`
}

type ValidationConfigV1 struct {
	Enable                     bool     `yaml:"enable" json:"enable"`
	ValidationLevel            string   `yaml:"validation_level" json:"validation_level"`
	StatisticalTests           []string `yaml:"statistical_tests" json:"statistical_tests"`
	PopulationGeneticsTests    []string `yaml:"population_genetics_tests" json:"population_genetics_tests"`
	PipelineCompatibilityTests []string `yaml:"pipeline_compatibility_tests" json:"pipeline_compatibility_tests"`
}

type RegulatoryConfigV1 struct {
	FDACompliance     bool     `yaml:"fda_compliance" json:"fda_compliance"`
	EMACompliance     bool     `yaml:"ema_compliance" json:"ema_compliance"`
	GCPCompliance     bool     `yaml:"gcp_compliance,omitempty" json:"gcp_compliance,omitempty"`
	ValidationStudies []string `yaml:"validation_studies" json:"validation_studies"`
	QualityMetrics    []string `yaml:"quality_metrics" json:"quality_metrics"`
}

type AcademicConfigV1 struct {
	PeerReviewReady           bool     `yaml:"peer_review_ready" json:"peer_review_ready"`
	PublicationFormat         string   `yaml:"publication_format" json:"publication_format"`
	BenchmarkDatasets         []string `yaml:"benchmark_datasets" json:"benchmark_datasets"`
	StatisticalPower          float64  `yaml:"statistical_power" json:"statistical_power"`
	EffectSizeReporting       bool     `yaml:"effect_size_reporting,omitempty" json:"effect_size_reporting,omitempty"`
	MultipleTestingCorrection string   `yaml:"multiple_testing_correction,omitempty" json:"multiple_testing_correction,omitempty"`
}

type SecurityConfigV1 struct {
	DataEncryption        string `yaml:"data_encryption" json:"data_encryption"`
	KeyManagement         string `yaml:"key_management" json:"key_management"`
	AccessControl         string `yaml:"access_control" json:"access_control"`
	NetworkSecurity       string `yaml:"network_security" json:"network_security"`
	DataMasking           bool   `yaml:"data_masking" json:"data_masking"`
	ZeroTrustArchitecture bool   `yaml:"zero_trust_architecture" json:"zero_trust_architecture"`
}

type VolumeTestingConfig struct {
	FastGeneration           bool `yaml:"fast_generation" json:"fast_generation"`
	MinimalValidation        bool `yaml:"minimal_validation" json:"minimal_validation"`
	SimpleSequences          bool `yaml:"simple_sequences" json:"simple_sequences"`
	ReducedQualityComplexity bool `yaml:"reduced_quality_complexity" json:"reduced_quality_complexity"`
	SkipCompression          bool `yaml:"skip_compression" json:"skip_compression"`
	BatchProcessing          bool `yaml:"batch_processing" json:"batch_processing"`
	MemoryOptimization       bool `yaml:"memory_optimization" json:"memory_optimization"`
	ParallelGeneration       bool `yaml:"parallel_generation" json:"parallel_generation"`
}

type OutputCustomizationConfig struct {
	Compression             string                 `yaml:"compression" json:"compression"`
	FileFormat              string                 `yaml:"file_format" json:"file_format"`
	BatchProcessing         bool                   `yaml:"batch_processing" json:"batch_processing"`
	RealTimeGeneration      bool                   `yaml:"real_time_generation" json:"real_time_generation"`
	QualityProfiles         string                 `yaml:"quality_profiles" json:"quality_profiles"`
	ContaminationSimulation ContaminationSimConfig `yaml:"contamination_simulation" json:"contamination_simulation"`
	AdapterSimulation       AdapterSimConfig       `yaml:"adapter_simulation" json:"adapter_simulation"`
	ErrorProfiling          ErrorProfilingConfig   `yaml:"error_profiling" json:"error_profiling"`
}

type ContaminationSimConfig struct {
	Enable bool `yaml:"enable" json:"enable"`
}

type AdapterSimConfig struct {
	Enable bool `yaml:"enable" json:"enable"`
}

type ErrorProfilingConfig struct {
	Enable bool `yaml:"enable" json:"enable"`
}

type PerformanceTargetsConfig struct {
	TargetThroughput  string `yaml:"target_throughput" json:"target_throughput"`
	MaxMemoryUsage    string `yaml:"max_memory_usage" json:"max_memory_usage"`
	TargetFileCount   int    `yaml:"target_file_count" json:"target_file_count"`
	ConcurrentBatches int    `yaml:"concurrent_batches" json:"concurrent_batches"`
}

type TestingScenario struct {
	Files        int    `yaml:"files" json:"files"`
	ReadsPerFile int    `yaml:"reads_per_file" json:"reads_per_file"`
	Description  string `yaml:"description" json:"description"`
}

type SimplifiedQualityConfig struct {
	FixedQualityScore      int  `yaml:"fixed_quality_score" json:"fixed_quality_score"`
	SimpleReadIDs          bool `yaml:"simple_read_ids" json:"simple_read_ids"`
	NoDuplicateChecking    bool `yaml:"no_duplicate_checking" json:"no_duplicate_checking"`
	DeterministicSequences bool `yaml:"deterministic_sequences" json:"deterministic_sequences"`
}

// Legacy support - convert old config to new format
func (old *SequencerConfig) ToV1() *SequencerConfigV1 {
	return &SequencerConfigV1{
		ConfigVersion: ConfigVersion{
			SchemaVersion: "1.0",
			LastUpdated:   time.Now().Format("2006-01-02"),
			ConfigType:    "sequencer",
		},
		Sequencer: SequencerHardware{
			Model:            old.SequencerModel,
			FlowCellType:     old.FlowCellType,
			ChemistryVersion: old.ChemistryVersion,
			InstrumentID:     old.InstrumentID,
			RunID:            old.RunID,
			FlowCellID:       old.FlowCellID,
		},
		Sample: SampleConfig{
			Name:      old.SampleName,
			Lane:      old.Lane,
			OutputDir: old.OutputDir,
		},
		Reads: ReadConfig{
			Length:        old.ReadLength,
			PairedEnd:     old.PairedEnd,
			ReadsPerBatch: old.ReadsPerBatch,
			TotalBatches:  old.TotalBatches,
			BatchInterval: old.BatchInterval,
		},
		Quality: QualityConfig{
			BaseQuality:           old.BaseQuality,
			QualityDecay:          old.QualityDecay,
			Q30Percentage:         old.Q30Percentage,
			ClusterDensityKPerMM2: old.ClusterDensityKPerMM2,
		},
		Performance: PerformanceConfig{
			RunTimeHours:           old.RunTimeHours,
			MaxOutputTB:            old.MaxOutputTB,
			MaxOutputGB:            old.MaxOutputGB,
			MaxCycles:              old.MaxCycles,
			LoadingConcentrationPM: old.LoadingConcentrationPM,
			PhixSpikeInPercent:     old.PhixSpikeInPercent,
		},
		Technology: TechnologyConfig{
			TwoChannelChemistry:  old.TwoChannelChemistry,
			FourChannelChemistry: old.FourChannelChemistry,
		},
		HMMModelPath: old.HMMModelPath,
		Population:   old.Population,
		ModelDir:     old.ModelDir,
	}
}

// Convert V1 config back to legacy format for compatibility
func (v1 *SequencerConfigV1) ToLegacy() *SequencerConfig {
	return &SequencerConfig{
		SequencerModel:         v1.Sequencer.Model,
		FlowCellType:           v1.Sequencer.FlowCellType,
		ChemistryVersion:       v1.Sequencer.ChemistryVersion,
		InstrumentID:           v1.Sequencer.InstrumentID,
		RunID:                  v1.Sequencer.RunID,
		FlowCellID:             v1.Sequencer.FlowCellID,
		SampleName:             v1.Sample.Name,
		Lane:                   v1.Sample.Lane,
		OutputDir:              v1.Sample.OutputDir,
		ReadLength:             v1.Reads.Length,
		ReadsPerBatch:          v1.Reads.ReadsPerBatch,
		BatchInterval:          v1.Reads.BatchInterval,
		TotalBatches:           v1.Reads.TotalBatches,
		PairedEnd:              v1.Reads.PairedEnd,
		BaseQuality:            v1.Quality.BaseQuality,
		QualityDecay:           v1.Quality.QualityDecay,
		Q30Percentage:          v1.Quality.Q30Percentage,
		ClusterDensityKPerMM2:  v1.Quality.ClusterDensityKPerMM2,
		LoadingConcentrationPM: v1.Performance.LoadingConcentrationPM,
		PhixSpikeInPercent:     v1.Performance.PhixSpikeInPercent,
		RunTimeHours:           v1.Performance.RunTimeHours,
		MaxOutputTB:            v1.Performance.MaxOutputTB,
		MaxOutputGB:            v1.Performance.MaxOutputGB,
		MaxCycles:              v1.Performance.MaxCycles,
		TwoChannelChemistry:    v1.Technology.TwoChannelChemistry,
		FourChannelChemistry:   v1.Technology.FourChannelChemistry,
		HMMModelPath:           v1.HMMModelPath,
		Population:             v1.Population,
		ModelDir:               v1.ModelDir,
		Illumina:               v1.Illumina,
	}
}
