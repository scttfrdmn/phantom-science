package types

import (
	"encoding/json"
	"reflect"
	"strings"
	"testing"
	"time"

	"gopkg.in/yaml.v3"
)

// The generation fields added for seeds, read lengths, error models, ONT
// metadata, output and scale must round-trip, and stay out of serialized
// configs that do not set them.
func TestSequencerConfig_GenerationFieldsRoundTrip(t *testing.T) {
	in := SequencerConfig{
		SampleName: "S",
		Platform:   "ont",
		Seed:       42,
		ReadLengthDistribution: ReadLengthSpec{
			Dist: ReadLengthLognormal, N50: 50000, Median: 25000, Min: 200, Max: 2000000,
			Mean: 30000, SD: 20000, Shape: 2, Scale: 750, EmpiricalPath: "h.tsv",
		},
		ErrorModel: ErrorModelSpec{Tier: "sup"},
		ONT:        ONTSpec{RunID: "r", FlowCellID: "PAM12345", ProtocolGroupID: "p", Channels: 3000, SpeedBPS: 400},
		Output:     OutputSpec{Compression: "gzip", FilePattern: "x", SummaryFile: "sequencing_summary.txt"},
		Illumina: IlluminaSpec{QualityModel: "linear", Errors: "quality",
			InsertSize: InsertSizeSpec{Dist: "lognormal", Mean: 400, SD: 60, Min: 150, Max: 1000}},
		TotalBases: 140e9, TotalReads: 1e6, FlowCells: 2, Workers: 8,
	}
	data, err := json.Marshal(in)
	if err != nil {
		t.Fatal(err)
	}
	var back SequencerConfig
	if err := json.Unmarshal(data, &back); err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(back, in) {
		t.Errorf("JSON round trip:\n got %+v\nwant %+v", back, in)
	}
	for _, key := range []string{`"seed":42`, `"platform":"ont"`, `"read_length_distribution":{`, `"n50":50000`, `"error_model":{"tier":"sup"}`, `"total_bases":140000000000`, `"illumina":{"quality_model":"linear"`, `"insert_size":{"dist":"lognormal"`} {
		if !strings.Contains(string(data), key) {
			t.Errorf("JSON lacks %s: %s", key, data)
		}
	}

	ydata, err := yaml.Marshal(in)
	if err != nil {
		t.Fatal(err)
	}
	var yback SequencerConfig
	if err := yaml.Unmarshal(ydata, &yback); err != nil {
		t.Fatal(err)
	}
	if yback.Seed != 42 || yback.Platform != "ont" || yback.ReadLengthDistribution != in.ReadLengthDistribution ||
		yback.ErrorModel != in.ErrorModel || yback.ONT != in.ONT || yback.Output != in.Output || yback.Illumina != in.Illumina ||
		yback.TotalBases != in.TotalBases || yback.TotalReads != in.TotalReads || yback.FlowCells != 2 || yback.Workers != 8 {
		t.Errorf("YAML round trip lost fields: %+v", yback)
	}
}

func TestSequencerConfig_GenerationFieldsOmittedWhenUnset(t *testing.T) {
	data, err := json.Marshal(SequencerConfig{SampleName: "S"})
	if err != nil {
		t.Fatal(err)
	}
	for _, key := range []string{"seed", "platform", "read_length_distribution", "error_model", `"ont"`, `"output"`, "illumina", "total_bases", "total_reads", "flow_cells", "workers"} {
		if strings.Contains(string(data), key) {
			t.Errorf("unset %s serialized: %s", key, data)
		}
	}
}

func TestSequencerConfigValidation(t *testing.T) {
	tests := []struct {
		name   string
		config SequencerConfig
		valid  bool
	}{
		{
			name: "valid NovaSeq config",
			config: SequencerConfig{
				SequencerModel:   "NovaSeq 6000",
				FlowCellType:     "S2",
				ChemistryVersion: "v1.5",
				InstrumentID:     "A00123",
				RunID:            "240101_A00123_0001_BH5TNVDRXY",
				FlowCellID:       "H5TNVDRXY",
				SampleName:       "Sample_001",
				Lane:             1,
				OutputDir:        "./output",
				ReadLength:       150,
				ReadsPerBatch:    1000000,
				BatchInterval:    "30m",
				TotalBatches:     10,
				PairedEnd:        true,
				BaseQuality:      35,
				QualityDecay:     0.05,
				Q30Percentage:    90.0,
			},
			valid: true,
		},
		{
			name: "invalid read length",
			config: SequencerConfig{
				ReadLength: -1,
			},
			valid: false,
		},
		{
			name: "invalid lane number",
			config: SequencerConfig{
				Lane:       0,
				ReadLength: 150,
			},
			valid: false,
		},
		{
			name: "invalid base quality",
			config: SequencerConfig{
				Lane:        1,
				ReadLength:  150,
				BaseQuality: -5,
			},
			valid: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			err := tt.config.Validate()
			if tt.valid && err != nil {
				t.Errorf("Expected valid config, got error: %v", err)
			}
			if !tt.valid && err == nil {
				t.Errorf("Expected invalid config, got no error")
			}
		})
	}
}

func TestSequencerConfigDefaults(t *testing.T) {
	config := &SequencerConfig{}
	config.SetDefaults()

	if config.BaseQuality == 0 {
		t.Error("Expected default base quality to be set")
	}
	if config.QualityDecay == 0 {
		t.Error("Expected default quality decay to be set")
	}
	if config.Q30Percentage == 0 {
		t.Error("Expected default Q30 percentage to be set")
	}
	if config.Lane == 0 {
		t.Error("Expected default lane to be set")
	}
}

func TestSequencerConfigV1Conversion(t *testing.T) {
	original := &SequencerConfig{
		SequencerModel:   "NovaSeq 6000",
		FlowCellType:     "S2",
		ChemistryVersion: "v1.5",
		InstrumentID:     "A00123",
		RunID:            "240101_A00123_0001_BH5TNVDRXY",
		FlowCellID:       "H5TNVDRXY",
		SampleName:       "Sample_001",
		Lane:             1,
		OutputDir:        "./output",
		ReadLength:       150,
		ReadsPerBatch:    1000000,
		BatchInterval:    "30m",
		TotalBatches:     10,
		PairedEnd:        true,
		BaseQuality:      35,
		QualityDecay:     0.05,
		Q30Percentage:    90.0,
		Illumina:         IlluminaSpec{Errors: IlluminaErrorsNone},
		ReferencePath:    "/ref/hg38.fa",
		ReferenceRegion:  "chr22",
		StrictRef:        true,
	}
	defer func() {
		back := original.ToV1().ToLegacy()
		if back.Illumina != original.Illumina || back.ReferencePath != original.ReferencePath ||
			back.ReferenceRegion != original.ReferenceRegion || back.StrictRef != original.StrictRef {
			t.Errorf("v1 round trip lost the illumina or reference settings: %+v", back)
		}
	}()

	// Convert to V1
	v1 := original.ToV1()
	if v1.SchemaVersion != "1.0" {
		t.Errorf("Expected schema version 1.0, got %s", v1.SchemaVersion)
	}
	if v1.ConfigType != "sequencer" {
		t.Errorf("Expected config type sequencer, got %s", v1.ConfigType)
	}
	if v1.Sequencer.Model != original.SequencerModel {
		t.Errorf("Expected model %s, got %s", original.SequencerModel, v1.Sequencer.Model)
	}

	// Convert back to legacy
	converted := v1.ToLegacy()
	if converted.SequencerModel != original.SequencerModel {
		t.Errorf("Expected model %s, got %s", original.SequencerModel, converted.SequencerModel)
	}
	if converted.ReadLength != original.ReadLength {
		t.Errorf("Expected read length %d, got %d", original.ReadLength, converted.ReadLength)
	}
	if converted.BaseQuality != original.BaseQuality {
		t.Errorf("Expected base quality %d, got %d", original.BaseQuality, converted.BaseQuality)
	}
}

func TestValidationResult(t *testing.T) {
	result := ValidationResult{
		TestName:     "test_validation",
		Passed:       true,
		PValue:       0.05,
		EffectSize:   0.2,
		Confidence:   0.95,
		Message:      "Test passed",
		Timestamp:    time.Now(),
		ComplianceOK: true,
	}

	if !result.Passed {
		t.Error("Expected test to pass")
	}
	if result.TestName == "" {
		t.Error("Expected test name to be set")
	}
	if result.Confidence != 0.95 {
		t.Errorf("Expected confidence 0.95, got %f", result.Confidence)
	}
}

func TestValidationResults(t *testing.T) {
	results := ValidationResults{
		OverallPass:     true,
		ComplianceLevel: "regulatory",
		Tests: []ValidationResult{
			{
				TestName:   "hardy_weinberg",
				Passed:     true,
				Confidence: 0.95,
			},
			{
				TestName:   "linkage_disequilibrium",
				Passed:     true,
				Confidence: 0.90,
			},
		},
		Summary:         "All tests passed",
		RegulatoryReady: true,
		AcademicReady:   true,
	}

	if !results.OverallPass {
		t.Error("Expected overall pass to be true")
	}
	if len(results.Tests) != 2 {
		t.Errorf("Expected 2 tests, got %d", len(results.Tests))
	}
	if !results.RegulatoryReady {
		t.Error("Expected regulatory ready to be true")
	}
	if !results.AcademicReady {
		t.Error("Expected academic ready to be true")
	}
}

func TestPhantomConfigV1Validation(t *testing.T) {
	config := PhantomConfigV1{
		ConfigVersion: ConfigVersion{
			SchemaVersion: "1.0",
			ConfigType:    "phantom_enterprise",
		},
		Enterprise: EnterpriseConfigV1{
			CustomerID:      "TEST_001",
			LicenseKey:      "test-license",
			ServiceTier:     "enterprise",
			ComplianceLevel: "hipaa_gdpr",
		},
		Fingerprinting: FingerprintConfigV1{
			Enable:               true,
			CustomerFingerprint:  true,
			DatasetFingerprint:   true,
			SteganographyLayers:  []string{"quality_scores", "read_ids"},
			CryptographicSigning: true,
			TamperDetection:      true,
		},
		Validation: ValidationConfigV1{
			Enable:          true,
			ValidationLevel: "regulatory",
			StatisticalTests: []string{
				"hardy_weinberg_equilibrium",
				"linkage_disequilibrium",
			},
		},
	}

	if config.SchemaVersion != "1.0" {
		t.Errorf("Expected schema version 1.0, got %s", config.SchemaVersion)
	}
	if !config.Fingerprinting.Enable {
		t.Error("Expected fingerprinting to be enabled")
	}
	if len(config.Fingerprinting.SteganographyLayers) != 2 {
		t.Errorf("Expected 2 steganography layers, got %d", len(config.Fingerprinting.SteganographyLayers))
	}
}

// Add validation methods to SequencerConfig for testing
func (c *SequencerConfig) Validate() error {
	if c.ReadLength <= 0 {
		return &ValidationError{Field: "ReadLength", Message: "must be positive"}
	}
	if c.Lane <= 0 {
		return &ValidationError{Field: "Lane", Message: "must be positive"}
	}
	if c.BaseQuality < 0 {
		return &ValidationError{Field: "BaseQuality", Message: "must be non-negative"}
	}
	return nil
}

func (c *SequencerConfig) SetDefaults() {
	if c.BaseQuality == 0 {
		c.BaseQuality = 30
	}
	if c.QualityDecay == 0 {
		c.QualityDecay = 0.1
	}
	if c.Q30Percentage == 0 {
		c.Q30Percentage = 85.0
	}
	if c.Lane == 0 {
		c.Lane = 1
	}
}

// ValidationError for testing
type ValidationError struct {
	Field   string
	Message string
}

func (e *ValidationError) Error() string {
	return e.Field + ": " + e.Message
}
func TestResolveFingerprinting(t *testing.T) {
	quality := FingerprintConfig{Enable: true, SteganographyLayers: []string{"quality_scores"}}
	readIDs := FingerprintConfig{Enable: true, SteganographyLayers: []string{"read_ids"}}

	tests := []struct {
		name         string
		seq          SequencerConfig
		phantom      *PhantomConfig
		wantLayers   []string
		wantEnable   bool
		wantCustomer string
	}{
		{
			name:         "top-level settings kept without phantom config",
			seq:          SequencerConfig{Fingerprinting: readIDs, Enterprise: EnterpriseConfig{CustomerID: "SEQ"}},
			wantEnable:   true,
			wantLayers:   []string{"read_ids"},
			wantCustomer: "SEQ",
		},
		{
			name:         "nested enterprise settings used when top-level is disabled",
			seq:          SequencerConfig{Enterprise: EnterpriseConfig{CustomerID: "SEQ", Fingerprinting: quality}},
			wantEnable:   true,
			wantLayers:   []string{"quality_scores"},
			wantCustomer: "SEQ",
		},
		{
			name:         "enabled phantom config overrides settings and customer",
			seq:          SequencerConfig{Fingerprinting: readIDs, Enterprise: EnterpriseConfig{CustomerID: "SEQ"}},
			phantom:      &PhantomConfig{Enterprise: EnterpriseConfig{CustomerID: "PH", Fingerprinting: quality}},
			wantEnable:   true,
			wantLayers:   []string{"quality_scores"},
			wantCustomer: "PH",
		},
		{
			name:         "disabled phantom config (CLI/TUI default) changes nothing",
			seq:          SequencerConfig{Fingerprinting: readIDs, Enterprise: EnterpriseConfig{CustomerID: "SEQ"}},
			phantom:      &PhantomConfig{Enterprise: EnterpriseConfig{CustomerID: "DEMO"}},
			wantEnable:   true,
			wantLayers:   []string{"read_ids"},
			wantCustomer: "SEQ",
		},
		{
			name:         "nothing enabled anywhere",
			seq:          SequencerConfig{Enterprise: EnterpriseConfig{CustomerID: "SEQ"}},
			phantom:      &PhantomConfig{},
			wantCustomer: "SEQ",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			c := tt.seq
			c.ResolveFingerprinting(tt.phantom)
			if c.Fingerprinting.Enable != tt.wantEnable {
				t.Errorf("Enable = %v, want %v", c.Fingerprinting.Enable, tt.wantEnable)
			}
			if len(c.Fingerprinting.SteganographyLayers) != len(tt.wantLayers) ||
				(len(tt.wantLayers) > 0 && c.Fingerprinting.SteganographyLayers[0] != tt.wantLayers[0]) {
				t.Errorf("layers = %v, want %v", c.Fingerprinting.SteganographyLayers, tt.wantLayers)
			}
			if c.Enterprise.CustomerID != tt.wantCustomer {
				t.Errorf("CustomerID = %q, want %q", c.Enterprise.CustomerID, tt.wantCustomer)
			}
		})
	}
}

func TestSequencerConfig_MaxOutputBases(t *testing.T) {
	for _, tc := range []struct {
		gb, tb float64
		want   int64
	}{
		{0, 0, 0},
		{15, 0, 15e9},
		{0, 3, 3e12},
		{5000, 3, 3e12}, // the smaller cap applies
		{1e-6, 0, 1000},
		{1e-12, 0, 1}, // a set cap never rounds down to "no cap"
		{-1, 0, 0},
	} {
		c := SequencerConfig{MaxOutputGB: tc.gb, MaxOutputTB: tc.tb}
		if got := c.MaxOutputBases(); got != tc.want {
			t.Errorf("gb=%g tb=%g: MaxOutputBases = %d, want %d", tc.gb, tc.tb, got, tc.want)
		}
	}
}

// IlluminaSpec defaults: binned qualities with quality errors, linear with
// none, and inserts of 350 ± 50; Validate rejects unknown names and
// impossible insert ranges.
func TestIlluminaSpec_DefaultsAndValidate(t *testing.T) {
	d := IlluminaSpec{}.WithDefaults()
	if d.QualityModel != IlluminaQualityBinned || d.Errors != IlluminaErrorsQuality ||
		d.InsertSize != (InsertSizeSpec{Dist: InsertSizeNormal, Mean: 350, SD: 50}) {
		t.Errorf("defaults %+v", d)
	}
	if l := (IlluminaSpec{QualityModel: IlluminaQualityLinear}).WithDefaults(); l.Errors != IlluminaErrorsNone {
		t.Errorf("linear defaults to errors %q, want none", l.Errors)
	}
	if l := (IlluminaSpec{QualityModel: IlluminaQualityLinear, Errors: IlluminaErrorsQuality}).WithDefaults(); l.Errors != IlluminaErrorsQuality {
		t.Errorf("explicit errors overridden: %q", l.Errors)
	}
	for _, ok := range []IlluminaSpec{{}, {QualityModel: "linear"}, {InsertSize: InsertSizeSpec{Dist: "lognormal", Mean: 500, SD: 100, Min: 100, Max: 2000}}} {
		if err := ok.Validate(); err != nil {
			t.Errorf("%+v: %v", ok, err)
		}
	}
	for _, bad := range []IlluminaSpec{
		{QualityModel: "rta2"}, {Errors: "indel"}, {InsertSize: InsertSizeSpec{Dist: "uniform"}},
		{InsertSize: InsertSizeSpec{Mean: -1}}, {InsertSize: InsertSizeSpec{SD: -1}},
		{InsertSize: InsertSizeSpec{Min: -1}}, {InsertSize: InsertSizeSpec{Min: 500, Max: 400}},
		{InsertSize: InsertSizeSpec{Max: 100}}, {InsertSize: InsertSizeSpec{Min: 1000}},
	} {
		if err := bad.Validate(); err == nil {
			t.Errorf("%+v accepted", bad)
		}
	}
}
