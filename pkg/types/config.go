// Package types defines shared configuration structures used across the
// phantom-science application, including SequencerConfig and PhantomConfig.
package types

import "time"

// SequencerConfig represents sequencer-realistic configuration
type SequencerConfig struct {
    // Sequencer Identity
    SequencerModel    string  `json:"sequencer_model"`
    FlowCellType      string  `json:"flow_cell_type"`
    ChemistryVersion  string  `json:"chemistry_version,omitempty"`
    InstrumentID      string  `json:"instrument_id"`
    RunID             string  `json:"run_id"`
    FlowCellID        string  `json:"flowcell_id"`
    
    // Basic Sequencer Parameters
    SampleName    string `json:"sample_name"`
    Lane          int    `json:"lane"`
    OutputDir     string `json:"output_dir"`
    ReadLength    int    `json:"read_length"`
    ReadsPerBatch int    `json:"reads_per_batch"`
    BatchInterval string `json:"batch_interval"`
    TotalBatches  int    `json:"total_batches"`
    PairedEnd     bool   `json:"paired_end"`
    BaseQuality   int    `json:"base_quality"`
    QualityDecay  float64 `json:"quality_decay"`
    
    // Dataset Information
    DatasetID     string `json:"dataset_id,omitempty"`
    
    // Enterprise Features
    Enterprise    EnterpriseConfig `json:"enterprise,omitempty"`
    Fingerprinting FingerprintConfig `json:"fingerprinting,omitempty"`
    
    // Quality and Performance Metrics
    Q30Percentage          float64 `json:"q30_percentage"`
    ClusterDensityKPerMM2  int     `json:"cluster_density_k_per_mm2,omitempty"`
    LoadingConcentrationPM float64 `json:"loading_concentration_pm,omitempty"`
    PhixSpikeInPercent     float64 `json:"phix_spike_in_percent,omitempty"`
    RunTimeHours           int     `json:"run_time_hours,omitempty"`
    MaxOutputTB            float64 `json:"max_output_tb,omitempty"`
    MaxOutputGB            float64 `json:"max_output_gb,omitempty"`
    MaxCycles              int     `json:"max_cycles,omitempty"`
    
    // Technology Flags
    TwoChannelChemistry  bool `json:"two_channel_chemistry,omitempty"`
    FourChannelChemistry bool `json:"four_channel_chemistry,omitempty"`

    // HMMModelPath is an optional path to a serialized HMM model (.hmm.json).
    // When set, the engine loads the model and uses it for biologically authentic
    // sequence generation instead of crypto/rand. If empty or the file cannot be
    // loaded, the engine falls back to crypto/rand transparently.
    // Takes priority over Population+ModelDir when both are set.
    HMMModelPath string `json:"hmm_model_path,omitempty" yaml:"hmm_model_path,omitempty"`

    // Population is a 1000 Genomes population code (e.g. "YRI", "CEU", "CHB").
    // When ModelDir is also set, the engine resolves ModelDir/<Population>.hmm.json
    // automatically. Ignored when HMMModelPath is set explicitly.
    Population string `json:"population,omitempty" yaml:"population,omitempty"`

    // ModelDir is a directory containing pre-trained .hmm.json model files named
    // <Population>.hmm.json (e.g. YRI.hmm.json). Used together with Population.
    // Ignored when HMMModelPath is set explicitly.
    ModelDir string `json:"model_dir,omitempty" yaml:"model_dir,omitempty"`
}

// PhantomConfig represents Phantom Science enterprise features
type PhantomConfig struct {
    // Enterprise Features
    Enterprise   EnterpriseConfig `json:"enterprise"`
    Validation   ValidationConfig `json:"validation"`
    
    // Scientific Configuration
    ReferenceGenome     GenomeConfig    `json:"reference_genome"`
    Contamination       ContamConfig    `json:"contamination"`
    AdapterSequences    AdapterConfig   `json:"adapter_sequences"`
    ErrorProfile        ErrorConfig     `json:"error_profile"`
    InsertSizeProfile   InsertConfig    `json:"insert_size_profile"`
    VariantSimulation   VariantConfig   `json:"variant_simulation"`
    AmpliconBias        AmpliconConfig  `json:"amplicon_bias"`
    CoverageProfile     CoverageConfig  `json:"coverage_profile"`
    AWSDataSources      AWSConfig       `json:"aws_data_sources"`
    DiseaseModels       DiseaseConfig   `json:"disease_models"`
    PopulationGenetics  PopGenConfig    `json:"population_genetics"`
    GenomicsIntegration GenomicsDBConfig `json:"genomics_integration"`
}

// EnterpriseConfig contains enterprise-specific settings
type EnterpriseConfig struct {
    CustomerID       string           `json:"customer_id"`
    LicenseKey string `json:"-" yaml:"-"`
    ComplianceLevel  string           `json:"compliance_level"`
    AuditLogging     AuditConfig      `json:"audit_logging"`
    Fingerprinting   FingerprintConfig `json:"fingerprinting"`
    ValidationLevel  string           `json:"validation_level"`
    ServiceTier      string           `json:"service_tier"`
}

type AuditConfig struct {
    Enable           bool   `json:"enable"`
    LogLevel         string `json:"log_level"`
    RetentionDays    int    `json:"retention_days"`
    ComplianceFormat string `json:"compliance_format"`
    ExportEndpoint   string `json:"export_endpoint"`
}

type FingerprintConfig struct {
    Enable              bool     `json:"enable"`
    CustomerFingerprint bool     `json:"customer_fingerprint"`
    DatasetFingerprint  bool     `json:"dataset_fingerprint"`
    SteganographyLayers []string `json:"steganography_layers"`
    CryptographicSigning bool    `json:"cryptographic_signing"`
    TamperDetection     bool     `json:"tamper_detection"`
}

type ValidationConfig struct {
    Enable                    bool                   `json:"enable"`
    StatisticalTests          []string              `json:"statistical_tests"`
    PopulationGeneticsTests   []string              `json:"population_genetics_tests"`
    PipelineCompatibilityTests []string             `json:"pipeline_compatibility_tests"`
    RegulatoryCompliance      RegulatoryConfig      `json:"regulatory_compliance"`
    AcademicValidation        AcademicConfig        `json:"academic_validation"`
}

type RegulatoryConfig struct {
    FDACompliance    bool     `json:"fda_compliance"`
    EMACompliance    bool     `json:"ema_compliance"`
    ValidationStudies []string `json:"validation_studies"`
    QualityMetrics   []string `json:"quality_metrics"`
}

type AcademicConfig struct {
    PeerReviewReady   bool     `json:"peer_review_ready"`
    PublicationFormat string   `json:"publication_format"`
    BenchmarkDatasets []string `json:"benchmark_datasets"`
    StatisticalPower  float64  `json:"statistical_power"`
}

// Additional config types (simplified for MVP)
type GenomeConfig struct {
    UseReference      bool               `json:"use_reference"`
    ReferencePath     string             `json:"reference_path"`
    ChromosomeWeights map[string]float64 `json:"chromosome_weights"`
    GCBias            float64            `json:"gc_bias"`
    DownloadReference bool               `json:"download_reference"`
    ReferenceSource   string             `json:"reference_source"`
}

// Validation result types
type ValidationResult struct {
    TestName     string    `json:"test_name"`
    Passed       bool      `json:"passed"`
    PValue       float64   `json:"p_value,omitempty"`
    EffectSize   float64   `json:"effect_size,omitempty"`
    Confidence   float64   `json:"confidence"`
    Message      string    `json:"message"`
    Timestamp    time.Time `json:"timestamp"`
    ComplianceOK bool      `json:"compliance_ok"`
}

type ValidationResults struct {
    OverallPass      bool               `json:"overall_pass"`
    ComplianceLevel  string             `json:"compliance_level"`
    Tests            []ValidationResult `json:"tests"`
    Summary          string             `json:"summary"`
    RegulatoryReady  bool               `json:"regulatory_ready"`
    AcademicReady    bool               `json:"academic_ready"`
}

// Placeholder types for remaining config sections
type ContamConfig struct{}
type AdapterConfig struct{}
type ErrorConfig struct{}
type InsertConfig struct{}
type VariantConfig struct{}
type AmpliconConfig struct{}
type CoverageConfig struct{}
type AWSConfig struct{}
type DiseaseConfig struct{}
type PopGenConfig struct{}
type GenomicsDBConfig struct{}