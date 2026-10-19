// Package types defines the configuration and result structures shared by
// the phantom-science packages: SequencerConfig (what to generate),
// PhantomConfig (fingerprinting and validation settings) and the validation
// result types.
//
// Several fields are parsed and carried through config conversion but not
// used by generation; their comments say so. They are kept so that existing
// configuration files keep loading.
package types

import "time"

// SequencerConfig describes one sequencing run to generate: instrument and
// run identity, read structure, quality model and output settings.
type SequencerConfig struct {
	// Sequencer Identity. SequencerModel and FlowCellType are shown in the
	// TUI and used by its timing simulator; ChemistryVersion and
	// InstrumentID are parsed but not used by generation.
	SequencerModel   string `json:"sequencer_model"`
	FlowCellType     string `json:"flow_cell_type"`
	ChemistryVersion string `json:"chemistry_version,omitempty"`
	InstrumentID     string `json:"instrument_id"`
	RunID            string `json:"run_id"`
	FlowCellID       string `json:"flowcell_id"`

	// Basic Sequencer Parameters. BatchInterval is parsed but not used by
	// generation, which writes files as fast as it can. BaseQuality and
	// QualityDecay set the Illumina quality model.
	SampleName    string  `json:"sample_name"`
	Lane          int     `json:"lane"`
	OutputDir     string  `json:"output_dir"`
	ReadLength    int     `json:"read_length"`
	ReadsPerBatch int     `json:"reads_per_batch"`
	BatchInterval string  `json:"batch_interval"`
	TotalBatches  int     `json:"total_batches"`
	PairedEnd     bool    `json:"paired_end"`
	BaseQuality   int     `json:"base_quality"`
	QualityDecay  float64 `json:"quality_decay"`

	// Dataset Information
	DatasetID string `json:"dataset_id,omitempty"`

	// Enterprise holds a customer ID and, for older configs, nested
	// fingerprinting settings (see ResolveFingerprinting). Fingerprinting
	// holds the settings generation uses.
	Enterprise     EnterpriseConfig  `json:"enterprise,omitempty"`
	Fingerprinting FingerprintConfig `json:"fingerprinting,omitempty"`

	// Quality and Performance Metrics. Only MaxOutputTB and MaxOutputGB
	// affect generation (see MaxOutputBases); the other fields here are
	// descriptive instrument specifications, parsed but not used.
	Q30Percentage          float64 `json:"q30_percentage"`
	ClusterDensityKPerMM2  int     `json:"cluster_density_k_per_mm2,omitempty"`
	LoadingConcentrationPM float64 `json:"loading_concentration_pm,omitempty"`
	PhixSpikeInPercent     float64 `json:"phix_spike_in_percent,omitempty"`
	RunTimeHours           int     `json:"run_time_hours,omitempty"`
	MaxOutputTB            float64 `json:"max_output_tb,omitempty"`
	MaxOutputGB            float64 `json:"max_output_gb,omitempty"`
	MaxCycles              int     `json:"max_cycles,omitempty"`

	// Technology Flags: parsed but not used by generation.
	TwoChannelChemistry  bool `json:"two_channel_chemistry,omitempty"`
	FourChannelChemistry bool `json:"four_channel_chemistry,omitempty"`

	// HMMModelPath is an optional path to a serialized HMM model (.hmm.json).
	// When set, the engine loads the model and samples bases from it instead
	// of from the seeded uniform generator. If empty or the file cannot be
	// loaded, the engine falls back to uniform bases (logging the load error).
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
	// default when empty) or "ont" (Oxford Nanopore: MinKNOW-style headers and
	// a sequencing_summary.txt).
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

	// Output-size targets and parallelism (#16). Zero means unset.
	//
	// TotalBases and TotalReads (at most one may be set) replace
	// TotalBatches: output files of ReadsPerBatch reads are written until
	// the target is met. TotalReads is exact (read pairs for paired-end
	// output; the last file is partial). TotalBases stops after the read
	// (pair) that reaches the target, so the output holds at least
	// TotalBases and less than TotalBases plus one read (pair) length.
	// Both apply per flow cell.
	TotalBases int64 `json:"total_bases,omitempty" yaml:"total_bases,omitempty"`
	TotalReads int64 `json:"total_reads,omitempty" yaml:"total_reads,omitempty"`
	// FlowCells > 1 generates that many independent runs, each in its own
	// subdirectory with a run seed, run ID and flow-cell ID derived from
	// Seed, and each sized by the targets above.
	FlowCells int `json:"flow_cells,omitempty" yaml:"flow_cells,omitempty"`
	// Workers is the number of output files generated in parallel (0 or 1:
	// one at a time). Output does not depend on it.
	Workers int `json:"workers,omitempty" yaml:"workers,omitempty"`
}

// MaxOutputBases returns the hard cap on output bases per flow cell set by
// MaxOutputGB and MaxOutputTB, which are read as gigabases and terabases
// (sequencer output is quoted in bases); the smaller applies when both are
// set. Zero means no cap.
func (c SequencerConfig) MaxOutputBases() int64 {
	var limit int64
	for _, v := range []float64{c.MaxOutputGB * 1e9, c.MaxOutputTB * 1e12} {
		if !(v > 0) {
			continue
		}
		n := int64(9e18)
		if v < 9e18 {
			n = max(int64(v), 1)
		}
		if limit == 0 || n < limit {
			limit = n
		}
	}
	return limit
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

// OutputSpec controls output files.
type OutputSpec struct {
	// Compression is "none" (default when empty) or "gzip".
	Compression string `json:"compression,omitempty" yaml:"compression,omitempty"`
	// FilePattern overrides the output file naming with placeholders
	// {sample_name}, {flowcell_id}, {run_id}, {file_number}, {barcode}, {lane}
	// and {mate}. Empty keeps the platform default: Illumina
	// <sample>_S1_L<lane>_[R<mate>_]<batch>.fastq, or ONT
	// <sample>_<flowcell>_<run>_<n>.fastq.
	FilePattern string `json:"file_pattern,omitempty" yaml:"file_pattern,omitempty"`
	// SummaryFile, when set, names a per-run read summary (ONT
	// sequencing_summary.txt) written next to the FASTQ files.
	SummaryFile string `json:"summary_file,omitempty" yaml:"summary_file,omitempty"`
}

// PhantomConfig holds the fingerprinting and validation settings loaded with
// -p. Only Enterprise (customer ID and fingerprinting) and Validation affect
// a run.
type PhantomConfig struct {
	Enterprise EnterpriseConfig `json:"enterprise"`
	Validation ValidationConfig `json:"validation"`

	// Scientific Configuration: placeholders with no fields, not used.
	ReferenceGenome     GenomeConfig     `json:"reference_genome"`
	Contamination       ContamConfig     `json:"contamination"`
	AdapterSequences    AdapterConfig    `json:"adapter_sequences"`
	ErrorProfile        ErrorConfig      `json:"error_profile"`
	InsertSizeProfile   InsertConfig     `json:"insert_size_profile"`
	VariantSimulation   VariantConfig    `json:"variant_simulation"`
	AmpliconBias        AmpliconConfig   `json:"amplicon_bias"`
	CoverageProfile     CoverageConfig   `json:"coverage_profile"`
	AWSDataSources      AWSConfig        `json:"aws_data_sources"`
	DiseaseModels       DiseaseConfig    `json:"disease_models"`
	PopulationGenetics  PopGenConfig     `json:"population_genetics"`
	GenomicsIntegration GenomicsDBConfig `json:"genomics_integration"`
}

// EnterpriseConfig holds the customer ID and fingerprinting settings.
// LicenseKey is never serialized and is not checked; ComplianceLevel,
// AuditLogging, ValidationLevel and ServiceTier are parsed but not used:
// setting them does not enable any compliance or audit behaviour.
type EnterpriseConfig struct {
	CustomerID      string            `json:"customer_id"`
	LicenseKey      string            `json:"-" yaml:"-"`
	ComplianceLevel string            `json:"compliance_level"`
	AuditLogging    AuditConfig       `json:"audit_logging"`
	Fingerprinting  FingerprintConfig `json:"fingerprinting"`
	ValidationLevel string            `json:"validation_level"`
	ServiceTier     string            `json:"service_tier"`
}

// AuditConfig is parsed from configuration but not used; no audit log is
// written.
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

// FingerprintConfig controls fingerprinting of generated data.
type FingerprintConfig struct {
	// Enable turns fingerprinting on.
	Enable bool `json:"enable"`
	// CustomerFingerprint and DatasetFingerprint tag each read ID with the
	// customer and dataset IDs.
	CustomerFingerprint bool `json:"customer_fingerprint"`
	DatasetFingerprint  bool `json:"dataset_fingerprint"`
	// SteganographyLayers: "quality_scores" embeds one fingerprint per FASTQ
	// file in its quality scores (requires the stegano library). Other
	// entries are recorded in the run manifest but have no effect.
	SteganographyLayers []string `json:"steganography_layers"`
	// KeyFile holds the key that encrypts quality-score fingerprints; when
	// empty the PHANTOM_FINGERPRINT_KEY environment variable is used.
	KeyFile string `json:"key_file,omitempty"`
	// Plaintext allows unencrypted quality-score fingerprints when no key is set.
	Plaintext bool `json:"plaintext,omitempty"`
	// CryptographicSigning and TamperDetection are parsed but not
	// implemented: no signature is written and no tampering is detected.
	CryptographicSigning bool `json:"cryptographic_signing"`
	TamperDetection      bool `json:"tamper_detection"`
}

// ValidationConfig lists the validation tests to run (see package
// internal/validation). AcademicValidation is parsed but not used;
// RegulatoryCompliance only sets the labels described on ValidationResults.
type ValidationConfig struct {
	Enable                     bool             `json:"enable"`
	StatisticalTests           []string         `json:"statistical_tests"`
	PopulationGeneticsTests    []string         `json:"population_genetics_tests"`
	PipelineCompatibilityTests []string         `json:"pipeline_compatibility_tests"`
	RegulatoryCompliance       RegulatoryConfig `json:"regulatory_compliance"`
	AcademicValidation         AcademicConfig   `json:"academic_validation"`
}

// RegulatoryConfig records regulatory intentions. The validation framework
// does not implement an FDA/EMA procedure, so these settings never make
// ValidationResults.RegulatoryReady true.
type RegulatoryConfig struct {
	FDACompliance     bool     `json:"fda_compliance"`
	EMACompliance     bool     `json:"ema_compliance"`
	ValidationStudies []string `json:"validation_studies"`
	QualityMetrics    []string `json:"quality_metrics"`
}

// AcademicConfig is parsed from configuration but not used.
type AcademicConfig struct {
	PeerReviewReady   bool     `json:"peer_review_ready"`
	PublicationFormat string   `json:"publication_format"`
	BenchmarkDatasets []string `json:"benchmark_datasets"`
	StatisticalPower  float64  `json:"statistical_power"`
}

// GenomeConfig is parsed from configuration but not used by generation.
type GenomeConfig struct {
	UseReference      bool               `json:"use_reference"`
	ReferencePath     string             `json:"reference_path"`
	ChromosomeWeights map[string]float64 `json:"chromosome_weights"`
	GCBias            float64            `json:"gc_bias"`
	DownloadReference bool               `json:"download_reference"`
	ReferenceSource   string             `json:"reference_source"`
}

// ValidationResult is the outcome of one validation test. Confidence is a
// fixed value chosen per test (the Hardy-Weinberg test reports its p-value),
// not a computed confidence. ComplianceOK equals Passed. EffectSize and
// PValue are not set by any current test.
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

// ValidationResults is the outcome of one validation-framework run.
type ValidationResults struct {
	// OverallPass is true only if at least one test ran and every test passed.
	OverallPass bool `json:"overall_pass"`
	// ComplianceLevel reflects only what ran: "none" (disabled, no tests, or
	// a test failed), "tests_passed" (all tests passed), or "academic" (all
	// tests passed and academic validation was requested). Config flags alone
	// never raise it, and no "regulatory" level is produced.
	ComplianceLevel string             `json:"compliance_level"`
	Tests           []ValidationResult `json:"tests"`
	Summary         string             `json:"summary"`
	// RegulatoryReady is always false: no FDA/EMA validation procedure is
	// implemented, so RegulatoryConfig flags cannot establish readiness. The
	// field is kept for JSON compatibility.
	RegulatoryReady bool `json:"regulatory_ready"`
	// AcademicReady is OverallPass && AcademicValidation.PeerReviewReady.
	AcademicReady bool `json:"academic_ready"`
}

// Placeholder types for PhantomConfig sections that have no settings yet.
// They accept any JSON object and ignore its contents.
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
