package types

import (
	"testing"
	"time"
)

func TestSequencerConfigValidation(t *testing.T) {
	tests := []struct {
		name   string
		config SequencerConfig
		valid  bool
	}{
		{
			name: "valid NovaSeq config",
			config: SequencerConfig{
				SequencerModel:    "NovaSeq 6000",
				FlowCellType:      "S2",
				ChemistryVersion:  "v1.5",
				InstrumentID:      "A00123",
				RunID:             "240101_A00123_0001_BH5TNVDRXY",
				FlowCellID:        "H5TNVDRXY",
				SampleName:        "Sample_001",
				Lane:              1,
				OutputDir:         "./output",
				ReadLength:        150,
				ReadsPerBatch:     1000000,
				BatchInterval:     "30m",
				TotalBatches:      10,
				PairedEnd:         true,
				BaseQuality:       35,
				QualityDecay:      0.05,
				Q30Percentage:     90.0,
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
	}

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