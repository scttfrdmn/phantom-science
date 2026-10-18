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

    // Platform selects the read model and output format: "illumina" (the
    // default when empty) or "ont" (Oxford Nanopore; not implemented yet, #12).
    Platform string `json:"platform,omitempty" yaml:"platform,omitempty"`

    // Seed makes generation reproducible: the same seed and configuration
    // produce byte-identical output. Zero means "draw a random seed"; the
    // drawn seed is logged to stderr so the run can be repeated.
    Seed uint64 `json:"seed,omitempty" yaml:"seed,omitempty"`

    // ReadLengthDistribution selects how read lengths are sampled. The zero
    // value (or Dist "fixed") gives every read exactly ReadLength bases.
    ReadLengthDistribution ReadLengthSpec `json:"read_length_distribution,omitzero" yaml:"read_length_distribution,omitempty"`

    // ErrorModel configures sequencing errors applied to generated reads.
    // The zero value applies no errors (#11).
    ErrorModel ErrorModelSpec `json:"error_model,omitzero" yaml:"error_model,omitempty"`

    // ONT holds Oxford Nanopore run metadata used when Platform is "ont".
    ONT ONTSpec `json:"ont,omitzero" yaml:"ont,omitempty"`

    // Output controls file naming, compression and run summary files.
    Output OutputSpec `json:"output,omitzero" yaml:"output,omitempty"`

    // Output-size targets and parallelism (#16). When TotalBases or
    // TotalReads is set, batch counts are derived from the read-length
    // distribution instead of TotalBatches. Zero means unset.
    TotalBases int64 `json:"total_bases,omitempty" yaml:"total_bases,omitempty"`
    TotalReads int64 `json:"total_reads,omitempty" yaml:"total_reads,omitempty"`
    FlowCells  int   `json:"flow_cells,omitempty" yaml:"flow_cells,omitempty"`
    Workers    int   `json:"workers,omitempty" yaml:"workers,omitempty"`
}

// Read-length distribution names for ReadLengthSpec.Dist.
const (
    ReadLengthFixed     = "fixed"
    ReadLengthLognormal = "lognormal"
    ReadLengthGamma     = "gamma"
    ReadLengthEmpirical = "empirical"
)

// ReadLengthSpec describes the read-length distribution (#10). Which fields
// apply depends on Dist: lognormal uses N50 and Median (or Mean and SD),
// gamma uses Shape and Scale (or Mean and SD), empirical loads a length
// histogram from EmpiricalPath. Min and Max clamp every sampler when set.
type ReadLengthSpec struct {
    Dist          string  `json:"dist,omitempty" yaml:"dist,omitempty"` // fixed | lognormal | gamma | empirical
    Mean          float64 `json:"mean,omitempty" yaml:"mean,omitempty"`
    SD            float64 `json:"sd,omitempty" yaml:"sd,omitempty"`
    Min           int     `json:"min,omitempty" yaml:"min,omitempty"`
    Max           int     `json:"max,omitempty" yaml:"max,omitempty"`
    N50           int     `json:"n50,omitempty" yaml:"n50,omitempty"`
    Median        int     `json:"median,omitempty" yaml:"median,omitempty"`
    Shape         float64 `json:"shape,omitempty" yaml:"shape,omitempty"`
    Scale         float64 `json:"scale,omitempty" yaml:"scale,omitempty"`
    EmpiricalPath string  `json:"empirical_path,omitempty" yaml:"empirical_path,omitempty"`
}

// ErrorModelSpec configures the sequencing error model (#11).
type ErrorModelSpec struct {
    // Tier is the ONT basecaller tier: "sup", "hac" or "fast". Empty means
    // no error model.
    Tier string `json:"tier,omitempty" yaml:"tier,omitempty"`
}

// ONTSpec holds Oxford Nanopore run metadata written to read headers and the
// sequencing summary (#12). Empty values are generated.
type ONTSpec struct {
    RunID           string  `json:"run_id,omitempty" yaml:"run_id,omitempty"`
    FlowCellID      string  `json:"flow_cell_id,omitempty" yaml:"flow_cell_id,omitempty"`
    ProtocolGroupID string  `json:"protocol_group_id,omitempty" yaml:"protocol_group_id,omitempty"`
    Channels        int     `json:"channels,omitempty" yaml:"channels,omitempty"`
    SpeedBPS        float64 `json:"speed_bps,omitempty" yaml:"speed_bps,omitempty"` // translocation speed, bases per second
}

// OutputSpec controls output files (#12, #16).
type OutputSpec struct {
    // Compression is "none" (default when empty) or "gzip".
    Compression string `json:"compression,omitempty" yaml:"compression,omitempty"`
    // FilePattern overrides the output file naming; empty keeps the
    // Illumina-style <sample>_S1_L<lane>_[R<mate>_]<batch>.fastq names.
    FilePattern string `json:"file_pattern,omitempty" yaml:"file_pattern,omitempty"`
    // SummaryFile, when set, names a per-run read summary (ONT
    // sequencing_summary.txt) written next to the FASTQ files.
    SummaryFile string `json:"summary_file,omitempty" yaml:"summary_file,omitempty"`
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

// ResolveFingerprinting sets c.Fingerprinting to the fingerprinting settings
// generation should use. A phantom config that enables fingerprinting takes
// precedence (with its customer ID); otherwise settings nested under the
// sequencer config's "enterprise" section are used when the top-level
// "fingerprinting" section does not enable them. phantom may be nil.
func (c *SequencerConfig) ResolveFingerprinting(phantom *PhantomConfig) {
    if phantom != nil && phantom.Enterprise.Fingerprinting.Enable {
        c.Fingerprinting = phantom.Enterprise.Fingerprinting
        if phantom.Enterprise.CustomerID != "" {
            c.Enterprise.CustomerID = phantom.Enterprise.CustomerID
        }
        return
    }
    if !c.Fingerprinting.Enable && c.Enterprise.Fingerprinting.Enable {
        c.Fingerprinting = c.Enterprise.Fingerprinting
    }
}

type FingerprintConfig struct {
    Enable              bool     `json:"enable"`
    CustomerFingerprint bool     `json:"customer_fingerprint"`
    DatasetFingerprint  bool     `json:"dataset_fingerprint"`
    SteganographyLayers []string `json:"steganography_layers"`
    // KeyFile holds the key that encrypts quality-score fingerprints; when
    // empty the PHANTOM_FINGERPRINT_KEY environment variable is used.
    KeyFile             string   `json:"key_file,omitempty"`
    // Plaintext allows unencrypted quality-score fingerprints when no key is set.
    Plaintext           bool     `json:"plaintext,omitempty"`
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