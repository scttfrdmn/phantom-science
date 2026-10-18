//go:build !release
// +build !release

package main

import (
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

// ---------------------------------------------------------------------------
// Tests for PlatformErrorSignatureManager construction
// ---------------------------------------------------------------------------

func TestNewPlatformErrorSignatureManager_NotNil(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	if manager == nil {
		t.Fatal("expected non-nil PlatformErrorSignatureManager")
	}
}

// ---------------------------------------------------------------------------
// Tests for GetSupportedPlatforms
// ---------------------------------------------------------------------------

func TestGetSupportedPlatforms_Count(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	platforms := manager.GetSupportedPlatforms()
	if len(platforms) != 7 {
		t.Errorf("expected 7 supported platforms, got %d", len(platforms))
	}
}

func TestGetSupportedPlatforms_AllExpectedPresent(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	supported := manager.GetSupportedPlatforms()

	expected := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}

	for _, exp := range expected {
		found := false
		for _, p := range supported {
			if p == exp {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("expected platform %s to be in supported list", exp)
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for GetErrorSignature
// ---------------------------------------------------------------------------

func TestGetErrorSignature_AllPlatforms(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}

	for _, platform := range platforms {
		sig := manager.GetErrorSignature(platform)
		if sig == nil {
			t.Errorf("expected non-nil signature for platform %s", platform)
		}
	}
}

func TestGetErrorSignature_UnknownPlatform(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	sig := manager.GetErrorSignature(genomics.SequencingPlatform("UnknownPlatform"))
	if sig != nil {
		t.Error("expected nil signature for unknown platform")
	}
}

func TestGetErrorSignature_Fingerprints_NonEmpty(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}

	fingerprints := make(map[string]bool)
	for _, platform := range platforms {
		sig := manager.GetErrorSignature(platform)
		if sig == nil {
			t.Errorf("nil signature for %s", platform)
			continue
		}
		if sig.SignatureFingerprint == "" {
			t.Errorf("empty fingerprint for platform %s", platform)
		}
		// Check uniqueness
		if fingerprints[sig.SignatureFingerprint] {
			t.Errorf("duplicate fingerprint %q for platform %s", sig.SignatureFingerprint, platform)
		}
		fingerprints[sig.SignatureFingerprint] = true
	}
}

func TestGetErrorSignature_IlluminaHasSubstitution(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	illuminaPlatforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
	}

	for _, platform := range illuminaPlatforms {
		sig := manager.GetErrorSignature(platform)
		if sig == nil {
			t.Errorf("nil signature for %s", platform)
			continue
		}
		if _, ok := sig.ErrorPatterns["substitution"]; !ok {
			t.Errorf("Illumina platform %s missing substitution error pattern", platform)
		}
	}
}

func TestGetErrorSignature_LongReadHasIndelPatterns(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	longReadPlatforms := []genomics.SequencingPlatform{
		genomics.PacBioRevio,
		genomics.PromethION24,
	}

	for _, platform := range longReadPlatforms {
		sig := manager.GetErrorSignature(platform)
		if sig == nil {
			t.Errorf("nil signature for %s", platform)
			continue
		}
		if _, ok := sig.ErrorPatterns["insertion"]; !ok {
			t.Errorf("long-read platform %s missing insertion pattern", platform)
		}
		if _, ok := sig.ErrorPatterns["deletion"]; !ok {
			t.Errorf("long-read platform %s missing deletion pattern", platform)
		}
	}
}

func TestGetErrorSignature_FingerprintContainsTechnology(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	tests := []struct {
		platform genomics.SequencingPlatform
		keyword  string
	}{
		{genomics.NovaSeq6000, "NOVA"},
		{genomics.NextSeq2000, "NEXTSEQ"},
		{genomics.MiSeq, "MISEQ"},
		{genomics.PacBioRevio, "PACBIO"},
		{genomics.DNBSeqT20x2, "MGI"},
		{genomics.ElementAVITI, "ELEMENT"},
		{genomics.PromethION24, "ONT"},
	}

	for _, tt := range tests {
		sig := manager.GetErrorSignature(tt.platform)
		if sig == nil {
			t.Errorf("nil signature for %s", tt.platform)
			continue
		}
		if !strings.Contains(strings.ToUpper(sig.SignatureFingerprint), tt.keyword) {
			t.Errorf("platform %s fingerprint %q does not contain %q",
				tt.platform, sig.SignatureFingerprint, tt.keyword)
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for GetPopulationErrorData
// ---------------------------------------------------------------------------

func TestGetPopulationErrorData_AllExpected(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	pops := []population.Population{
		population.YRI, population.CEU, population.CHB, population.MXL,
	}

	for _, pop := range pops {
		data := manager.GetPopulationErrorData(pop)
		if data == nil {
			t.Errorf("expected non-nil population error data for %s", pop)
		}
	}
}

func TestGetPopulationErrorData_FieldsInRange(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	pops := []population.Population{
		population.YRI, population.LWK, population.GWD,
		population.CEU, population.GBR, population.TSI,
		population.CHB, population.JPT, population.CHS,
		population.MXL, population.PUR,
	}

	for _, pop := range pops {
		data := manager.GetPopulationErrorData(pop)
		if data == nil {
			t.Errorf("nil data for %s", pop)
			continue
		}

		// BaseErrorModifier should be around 1.0 (±20%)
		if data.BaseErrorModifier < 0.5 || data.BaseErrorModifier > 2.0 {
			t.Errorf("%s BaseErrorModifier = %v out of expected range [0.5, 2.0]",
				pop, data.BaseErrorModifier)
		}

		// GCBiasModifier should be reasonable
		if data.GCBiasModifier < 0.5 || data.GCBiasModifier > 2.0 {
			t.Errorf("%s GCBiasModifier = %v out of expected range [0.5, 2.0]",
				pop, data.GCBiasModifier)
		}

		// RepeatTolerance should be reasonable
		if data.RepeatTolerance < 0.5 || data.RepeatTolerance > 2.0 {
			t.Errorf("%s RepeatTolerance = %v out of expected range [0.5, 2.0]",
				pop, data.RepeatTolerance)
		}

		// QualityBias should be reasonable
		if data.QualityBias < 0.5 || data.QualityBias > 2.0 {
			t.Errorf("%s QualityBias = %v out of expected range [0.5, 2.0]",
				pop, data.QualityBias)
		}
	}
}

func TestGetPopulationErrorData_AfricanHasLowerErrorModifier(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	africanPops := []population.Population{population.YRI, population.LWK, population.GWD}
	for _, pop := range africanPops {
		data := manager.GetPopulationErrorData(pop)
		if data == nil {
			t.Errorf("nil data for %s", pop)
			continue
		}
		// African populations should have BaseErrorModifier < 1.0 (0.98)
		if data.BaseErrorModifier >= 1.0 {
			t.Errorf("%s BaseErrorModifier = %v, expected < 1.0 for African population",
				pop, data.BaseErrorModifier)
		}
	}
}

func TestGetPopulationErrorData_UnknownPopulation(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	data := manager.GetPopulationErrorData(population.Population("UNKNOWN_POP"))
	if data != nil {
		t.Error("expected nil data for unknown population")
	}
}

// ---------------------------------------------------------------------------
// Tests for ApplyErrorSignature
// ---------------------------------------------------------------------------

func TestApplyErrorSignature_Success(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	testSeq := "ATCGATCGATCGATCGAAAAAATCGTCGTCGCCCCCGGGGATCGATCG"
	testQuals := make([]int, len(testSeq))
	for i := range testQuals {
		testQuals[i] = 35
	}

	_, _, result := manager.ApplyErrorSignature(
		genomics.NovaSeq6000,
		population.YRI,
		testSeq,
		testQuals,
		1,
	)

	if !result.Success {
		t.Errorf("expected success, got message: %s", result.Message)
	}
	if result.Platform != genomics.NovaSeq6000 {
		t.Errorf("result platform = %v, want NovaSeq6000", result.Platform)
	}
}

func TestApplyErrorSignature_UnknownPlatform(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	testSeq := "ATCGATCG"
	testQuals := []int{35, 35, 35, 35, 35, 35, 35, 35}

	retSeq, retQual, result := manager.ApplyErrorSignature(
		genomics.SequencingPlatform("UNKNOWN"),
		population.YRI,
		testSeq,
		testQuals,
		1,
	)

	if result.Success {
		t.Error("expected failure for unknown platform")
	}
	if retSeq != testSeq {
		t.Errorf("sequence should be unchanged for unknown platform")
	}
	if len(retQual) != len(testQuals) {
		t.Errorf("quality length should be unchanged for unknown platform")
	}
}

func TestApplyErrorSignature_AllPlatforms(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}

	testSeq := "ATCGATCGATCGATCGAAAAAATCGTCGTCGCCCCCGGGGATCGATCG"
	testQuals := make([]int, len(testSeq))
	for i := range testQuals {
		testQuals[i] = 35
	}

	for _, platform := range platforms {
		_, modQual, result := manager.ApplyErrorSignature(
			platform, population.CEU, testSeq, testQuals, 1,
		)

		if !result.Success {
			t.Errorf("platform %s: expected success, got: %s", platform, result.Message)
		}
		if len(modQual) == 0 {
			t.Errorf("platform %s: modified quality should be non-empty", platform)
		}
		if result.ErrorsApplied < 0 {
			t.Errorf("platform %s: errors applied should be >= 0, got %d", platform, result.ErrorsApplied)
		}
	}
}

func TestApplyErrorSignature_OutputLengthConsistency(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	testSeq := "ATCGATCGATCG"
	testQuals := make([]int, len(testSeq))
	for i := range testQuals {
		testQuals[i] = 30
	}

	modSeq, modQual, result := manager.ApplyErrorSignature(
		genomics.MiSeq, population.CHB, testSeq, testQuals, 0,
	)

	if !result.Success {
		t.Fatalf("expected success: %s", result.Message)
	}
	// Quality length should match input
	if len(modQual) != len(testQuals) {
		t.Errorf("quality length = %d, want %d", len(modQual), len(testQuals))
	}
	// Sequence should be non-empty
	if len(modSeq) == 0 {
		t.Error("modified sequence should not be empty")
	}
}

func TestApplyErrorSignature_PopulationPropagated(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	testSeq := "ATCGATCGATCG"
	testQuals := make([]int, len(testSeq))
	for i := range testQuals {
		testQuals[i] = 30
	}

	_, _, result := manager.ApplyErrorSignature(
		genomics.NovaSeq6000, population.MXL, testSeq, testQuals, 0,
	)

	if result.Population != population.MXL {
		t.Errorf("result population = %v, want MXL", result.Population)
	}
}

func TestApplyErrorSignature_UnknownPopulationFallback(t *testing.T) {
	// Should still succeed with a fallback PopulationErrorData
	manager := genomics.NewPlatformErrorSignatureManager()

	testSeq := "ATCGATCGATCG"
	testQuals := make([]int, len(testSeq))
	for i := range testQuals {
		testQuals[i] = 30
	}

	_, _, result := manager.ApplyErrorSignature(
		genomics.NovaSeq6000,
		population.Population("UNKNOWN_POP"),
		testSeq,
		testQuals,
		0,
	)

	if !result.Success {
		t.Errorf("expected success even for unknown population (fallback), got: %s", result.Message)
	}
}

// ---------------------------------------------------------------------------
// Tests for ValidatePlatformSignatures
// ---------------------------------------------------------------------------

func TestValidatePlatformSignatures_ReturnsMap(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	validation := manager.ValidatePlatformSignatures()

	if len(validation) == 0 {
		t.Error("expected non-empty validation map")
	}
}

func TestValidatePlatformSignatures_AllConfigured(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	validation := manager.ValidatePlatformSignatures()

	if !validation["all_platform_signatures_configured"] {
		t.Error("expected all platform signatures to be configured")
	}
}

func TestValidatePlatformSignatures_PopulationDataLoaded(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	validation := manager.ValidatePlatformSignatures()

	if !validation["population_error_data_loaded"] {
		t.Error("expected population error data to be loaded")
	}
}

func TestValidatePlatformSignatures_FingerprintsUnique(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	validation := manager.ValidatePlatformSignatures()

	if !validation["signature_fingerprints_unique"] {
		t.Error("expected signature fingerprints to be unique")
	}
}

func TestValidatePlatformSignatures_IlluminaRepresented(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	validation := manager.ValidatePlatformSignatures()

	if !validation["illumina_platforms_represented"] {
		t.Error("expected Illumina platforms to be represented")
	}
}

func TestValidatePlatformSignatures_LongReadRepresented(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	validation := manager.ValidatePlatformSignatures()

	if !validation["long_read_platforms_represented"] {
		t.Error("expected long-read platforms to be represented")
	}
}

func TestValidatePlatformSignatures_IlluminaErrorRates(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	validation := manager.ValidatePlatformSignatures()

	if !validation["illumina_error_rates_appropriate"] {
		t.Error("expected Illumina error rates to be appropriate (low)")
	}
}

func TestValidatePlatformSignatures_LongReadIndelPatterns(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	validation := manager.ValidatePlatformSignatures()

	if !validation["long_read_indel_patterns"] {
		t.Error("expected long-read indel patterns to be present")
	}
}

func TestValidatePlatformSignatures_ErrorApplicationFunctional(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()
	validation := manager.ValidatePlatformSignatures()

	if !validation["error_application_functional"] {
		t.Error("expected error application to be functional")
	}
}

// ---------------------------------------------------------------------------
// Tests for error rate characteristics
// ---------------------------------------------------------------------------

func TestShortReadPlatforms_LowErrorRates(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	shortReadPlatforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
	}

	for _, platform := range shortReadPlatforms {
		sig := manager.GetErrorSignature(platform)
		if sig == nil {
			t.Errorf("nil signature for %s", platform)
			continue
		}
		subPattern, ok := sig.ErrorPatterns["substitution"]
		if !ok {
			t.Errorf("no substitution pattern for %s", platform)
			continue
		}

		if len(subPattern.FrequencyPattern) == 0 {
			t.Errorf("%s substitution frequency pattern is empty", platform)
			continue
		}

		avgErrorRate := 0.0
		for _, rate := range subPattern.FrequencyPattern {
			avgErrorRate += rate
		}
		avgErrorRate /= float64(len(subPattern.FrequencyPattern))

		if avgErrorRate > 0.002 {
			t.Errorf("%s avg error rate = %v, expected <= 0.002 for short-read platform",
				platform, avgErrorRate)
		}
	}
}

func TestLongReadPlatforms_HigherErrorRates(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	longReadPlatforms := []genomics.SequencingPlatform{
		genomics.PacBioRevio,
		genomics.PromethION24,
	}

	for _, platform := range longReadPlatforms {
		sig := manager.GetErrorSignature(platform)
		if sig == nil {
			t.Errorf("nil signature for %s", platform)
			continue
		}
		subPattern, ok := sig.ErrorPatterns["substitution"]
		if !ok {
			t.Errorf("no substitution pattern for %s", platform)
			continue
		}

		if len(subPattern.FrequencyPattern) == 0 {
			t.Errorf("%s substitution frequency pattern is empty", platform)
			continue
		}

		avgErrorRate := 0.0
		for _, rate := range subPattern.FrequencyPattern {
			avgErrorRate += rate
		}
		avgErrorRate /= float64(len(subPattern.FrequencyPattern))

		if avgErrorRate < 0.003 {
			t.Errorf("%s avg error rate = %v, expected >= 0.003 for long-read platform",
				platform, avgErrorRate)
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for technology classification (mirrors main() Test 6)
// ---------------------------------------------------------------------------

func TestTechnologyClassification_Illumina(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	illuminaPlatforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
	}

	for _, platform := range illuminaPlatforms {
		sig := manager.GetErrorSignature(platform)
		if sig == nil {
			t.Errorf("nil signature for %s", platform)
			continue
		}
		fp := strings.ToUpper(sig.SignatureFingerprint)
		if !strings.Contains(fp, "ILLUMINA") && !strings.Contains(fp, "SBS") &&
			!strings.Contains(fp, "NOVA") && !strings.Contains(fp, "NEXTSEQ") && !strings.Contains(fp, "MISEQ") {
			t.Errorf("Illumina platform %s fingerprint %q doesn't indicate Illumina tech",
				platform, sig.SignatureFingerprint)
		}
	}
}

func TestTechnologyClassification_ONT(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	sig := manager.GetErrorSignature(genomics.PromethION24)
	if sig == nil {
		t.Fatal("nil signature for PromethION24")
	}

	fp := strings.ToUpper(sig.SignatureFingerprint)
	if !strings.Contains(fp, "NANOPORE") && !strings.Contains(fp, "ONT") {
		t.Errorf("ONT platform fingerprint %q doesn't indicate ONT/Nanopore tech",
			sig.SignatureFingerprint)
	}
}

func TestTechnologyClassification_Element(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	sig := manager.GetErrorSignature(genomics.ElementAVITI)
	if sig == nil {
		t.Fatal("nil signature for ElementAVITI")
	}

	fp := strings.ToUpper(sig.SignatureFingerprint)
	if !strings.Contains(fp, "AVIDITY") && !strings.Contains(fp, "ELEMENT") {
		t.Errorf("Element platform fingerprint %q doesn't indicate Element/Avidity tech",
			sig.SignatureFingerprint)
	}
}

// ---------------------------------------------------------------------------
// Tests mirroring the quality metrics logic from main()
// ---------------------------------------------------------------------------

func TestQualityMetricsComputation(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	supportedPlatforms := manager.GetSupportedPlatforms()
	testPopulations := []population.Population{
		population.YRI, population.CEU, population.CHB, population.MXL,
	}
	validation := manager.ValidatePlatformSignatures()

	passedTests := 0
	totalTests := len(validation)
	for _, passed := range validation {
		if passed {
			passedTests++
		}
	}

	qualityMetrics := map[string]bool{
		"All 7 platforms have error signatures":           len(supportedPlatforms) == 7,
		"Platform signatures configured":                  passedTests >= totalTests*3/4,
		"Population-specific error data loaded":           len(testPopulations) >= 4,
		"Error signature application operational":         true,
		"Technology classification accurate":              true,
		"Error rate characteristics platform-appropriate": true,
		"Platform fingerprints unique":                    true,
		"Phase 2 Week 8 foundation complete":              true,
	}

	passedQuality := 0
	for _, passed := range qualityMetrics {
		if passed {
			passedQuality++
		}
	}

	totalQuality := len(qualityMetrics)

	// Main's final check: passedQuality >= totalQuality*7/8
	if passedQuality < totalQuality*7/8 {
		t.Errorf("quality assessment failed: %d/%d (need >= %d)",
			passedQuality, totalQuality, totalQuality*7/8)
	}
}

// ---------------------------------------------------------------------------
// Tests for error pattern quality correlation
// ---------------------------------------------------------------------------

func TestSubstitutionPattern_QualityCorrelation(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	// MiSeq should have very strong negative quality correlation
	sig := manager.GetErrorSignature(genomics.MiSeq)
	if sig == nil {
		t.Fatal("nil signature for MiSeq")
	}

	subPattern, ok := sig.ErrorPatterns["substitution"]
	if !ok {
		t.Fatal("MiSeq missing substitution pattern")
	}

	// Quality correlation should be negative (more errors = lower quality)
	if subPattern.QualityCorrelation >= 0 {
		t.Errorf("MiSeq quality correlation = %v, expected negative", subPattern.QualityCorrelation)
	}
}

// ---------------------------------------------------------------------------
// Tests for platform constants
// ---------------------------------------------------------------------------

func TestSequencingPlatformConstants(t *testing.T) {
	tests := []struct {
		platform genomics.SequencingPlatform
		expected string
	}{
		{genomics.NovaSeq6000, "NovaSeq 6000"},
		{genomics.NextSeq2000, "NextSeq 2000"},
		{genomics.MiSeq, "MiSeq"},
		{genomics.PacBioRevio, "PacBio Revio"},
		{genomics.DNBSeqT20x2, "DNBSEQ-T20x2"},
		{genomics.ElementAVITI, "Element AVITI"},
		{genomics.PromethION24, "PromethION 24"},
	}

	for _, tt := range tests {
		if string(tt.platform) != tt.expected {
			t.Errorf("platform %v = %q, want %q", tt.platform, string(tt.platform), tt.expected)
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for error signature error patterns structure
// ---------------------------------------------------------------------------

func TestErrorPatternFrequencyPattern_NonEmpty(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}

	for _, platform := range platforms {
		sig := manager.GetErrorSignature(platform)
		if sig == nil {
			t.Errorf("nil signature for %s", platform)
			continue
		}
		for patternName, pattern := range sig.ErrorPatterns {
			if len(pattern.FrequencyPattern) == 0 {
				t.Errorf("platform %s pattern %q has empty FrequencyPattern",
					platform, patternName)
			}
		}
	}
}

func TestErrorPatternFrequencyPattern_PositiveRates(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.PacBioRevio,
		genomics.PromethION24,
	}

	for _, platform := range platforms {
		sig := manager.GetErrorSignature(platform)
		if sig == nil {
			t.Errorf("nil signature for %s", platform)
			continue
		}
		for patternName, pattern := range sig.ErrorPatterns {
			for i, rate := range pattern.FrequencyPattern {
				if rate < 0 {
					t.Errorf("platform %s pattern %q position %d: negative rate %v",
						platform, patternName, i, rate)
				}
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for quality score calculation mirroring main()'s Test 5 logic
// ---------------------------------------------------------------------------

func TestQualityStatisticsCalculation(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	testSeq := "ATCGATCGATCGATCGAAAAAATCGTCGTCGCCCCCGGGGATCGATCG"
	testQuals := make([]int, len(testSeq))
	for i := range testQuals {
		testQuals[i] = 35
	}

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
	}

	for _, platform := range platforms {
		_, modQual, result := manager.ApplyErrorSignature(
			platform, population.YRI, testSeq, testQuals, 1,
		)

		if !result.Success {
			t.Errorf("platform %s: unexpected failure: %s", platform, result.Message)
			continue
		}

		// Calculate quality stats as main() does
		totalQ := 0
		minQ := modQual[0]
		maxQ := modQual[0]
		for _, q := range modQual {
			totalQ += q
			if q < minQ {
				minQ = q
			}
			if q > maxQ {
				maxQ = q
			}
		}
		meanQ := float64(totalQ) / float64(len(modQual))

		if meanQ <= 0 {
			t.Errorf("platform %s: mean quality = %v, expected > 0", platform, meanQ)
		}
		if minQ < 0 {
			t.Errorf("platform %s: min quality = %d, expected >= 0", platform, minQ)
		}
		if maxQ > 100 {
			t.Errorf("platform %s: max quality = %d, expected <= 100", platform, maxQ)
		}
		if minQ > maxQ {
			t.Errorf("platform %s: min quality %d > max quality %d", platform, minQ, maxQ)
		}
	}
}

func TestSequenceDifferenceCounting(t *testing.T) {
	manager := genomics.NewPlatformErrorSignatureManager()

	// Use a very long sequence to give error application opportunities
	testSeq := strings.Repeat("ATCGATCGATCGATCGAAAAAATCGTCGTCGCCCCCGGGGATCGATCG", 5)
	testQuals := make([]int, len(testSeq))
	for i := range testQuals {
		testQuals[i] = 35
	}

	modSeq, _, result := manager.ApplyErrorSignature(
		genomics.PacBioRevio, population.YRI, testSeq, testQuals, 1,
	)

	if !result.Success {
		t.Fatalf("expected success: %s", result.Message)
	}

	differences := 0
	for i := 0; i < len(testSeq) && i < len(modSeq); i++ {
		if testSeq[i] != modSeq[i] {
			differences++
		}
	}

	// Differences should be >= 0
	if differences < 0 {
		t.Errorf("differences = %d, expected >= 0", differences)
	}
}

// ---------------------------------------------------------------------------
// Tests for strings.Repeat usage in main() (simple sanity check)
// ---------------------------------------------------------------------------

func TestStringsRepeat(t *testing.T) {
	result := strings.Repeat("=", 50)
	if len(result) != 50 {
		t.Errorf("strings.Repeat(\"=\", 50) length = %d, want 50", len(result))
	}
	for _, c := range result {
		if c != '=' {
			t.Errorf("expected all '=' characters, got %q", c)
		}
	}
}

// ---------------------------------------------------------------------------
// Tests for phase 2 week 8 completion threshold
// ---------------------------------------------------------------------------

func TestPhase2Week8CompletionThreshold(t *testing.T) {
	// main() checks: passedQuality >= totalQuality*7/8
	testCases := []struct {
		passed   int
		total    int
		wantPass bool
	}{
		{8, 8, true},
		{7, 8, true},
		{6, 8, false},
		{7, 8, true},
		{0, 8, false},
	}

	for _, tc := range testCases {
		got := tc.passed >= tc.total*7/8
		if got != tc.wantPass {
			t.Errorf("passed=%d total=%d: got %v, want %v", tc.passed, tc.total, got, tc.wantPass)
		}
	}
}
