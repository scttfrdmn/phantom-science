//go:build !release
// +build !release

package populations

import (
	"testing"
)

// White-box tests for private helper methods of AfricanPopulationEnhancer.

func TestGetVariantChromosome_KnownVariant(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	tests := []struct {
		variant string
		want    string
	}{
		{"rs429358", "19"},
		{"rs7412", "19"},
		{"rs1815739", "11"},
		{"rs4988235", "2"},
		{"rs4244285", "10"},
		{"rs776746", "7"},
		{"rs_unknown", "1"}, // default
	}
	for _, tt := range tests {
		got := ape.getVariantChromosome(tt.variant)
		if got != tt.want {
			t.Errorf("getVariantChromosome(%q) = %q, want %q", tt.variant, got, tt.want)
		}
	}
}

func TestGetVariantPosition_KnownVariant(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	tests := []struct {
		variant string
		want    int64
	}{
		{"rs429358", 45411941},
		{"rs7412", 45412079},
		{"rs1815739", 66560624},
		{"rs4988235", 136608646},
		{"rs4244285", 94781859},
		{"rs776746", 99270539},
		{"rs_unknown", 1000000}, // default
	}
	for _, tt := range tests {
		got := ape.getVariantPosition(tt.variant)
		if got != tt.want {
			t.Errorf("getVariantPosition(%q) = %d, want %d", tt.variant, got, tt.want)
		}
	}
}

func TestGenerateZygosity_ReturnsValidValue(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	// randomFloat always returns 0.5; test coverage of each return branch
	validZygosities := map[string]bool{"0/0": true, "0/1": true, "1/1": true}
	for _, freq := range []float64{0.01, 0.5, 0.99} {
		got := ape.generateZygosity(freq)
		if !validZygosities[got] {
			t.Errorf("generateZygosity(%f) = %q, want one of 0/0, 0/1, 1/1", freq, got)
		}
	}
}

func TestAbsPopulations(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	_ = ape // suppress unused warning; abs is package-level
	tests := []struct{ in, want float64 }{
		{1.0, 1.0},
		{-1.0, 1.0},
		{0.0, 0.0},
	}
	for _, tt := range tests {
		got := abs(tt.in)
		if got != tt.want {
			t.Errorf("abs(%f) = %f, want %f", tt.in, got, tt.want)
		}
	}
}

// --- validateEnhancedModel regression tests ---

// TestValidateEnhancedModel_NoDemographicData checks that metrics that require extended
// population data (DemographicParams, LinkagePatterns, StructuralVariants,
// PharmacogenomicVariants) are set to 0.0, not to hardcoded constants, when those
// fields are absent from the population struct.
func TestValidateEnhancedModel_NoDemographicData(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	basePop := &H3AfricaPopulation{
		PopulationID: "TEST",
		AlleleFrequencies: map[string]float64{
			"rs429358": 0.205,
			"rs7412":   0.123,
		},
		// DemographicParams, LinkagePatterns, StructuralVariants,
		// PharmacogenomicVariants intentionally absent.
	}
	params := &SyntheticGenParams{
		StructuralVarRate: 0.0008,
		PopulationStructure: &PopulationStructure{
			Fst: map[string]float64{"EUR": 0.15},
		},
		// AlleleFreqs must match basePop to enable accuracy calculation.
		AlleleFreqs: map[string]float64{
			"rs429358": 0.205,
			"rs7412":   0.123,
		},
	}

	metrics := ape.validateEnhancedModel(basePop, params)

	// AlleleFreqAccuracy must be computed (both variants have non-zero freq).
	if metrics.AlleleFreqAccuracy <= 0 {
		t.Errorf("AlleleFreqAccuracy = %v, want > 0", metrics.AlleleFreqAccuracy)
	}
	// Metrics requiring absent data must be 0.0, not hardcoded constants like 0.88.
	if metrics.DemographicFit != 0.0 {
		t.Errorf("DemographicFit = %v, want 0.0 (no DemographicParams)", metrics.DemographicFit)
	}
	if metrics.LinkageAccuracy != 0.0 {
		t.Errorf("LinkageAccuracy = %v, want 0.0 (no LinkagePatterns)", metrics.LinkageAccuracy)
	}
	if metrics.StructuralVarAccuracy != 0.0 {
		t.Errorf("StructuralVarAccuracy = %v, want 0.0 (no StructuralVariants)", metrics.StructuralVarAccuracy)
	}
	if metrics.PharmacogenomicFit != 0.0 {
		t.Errorf("PharmacogenomicFit = %v, want 0.0 (no PharmacogenomicVariants)", metrics.PharmacogenomicFit)
	}
	// OverallValidationScore must be > 0 because AlleleFreqAccuracy contributes.
	if metrics.OverallValidationScore <= 0 {
		t.Errorf("OverallValidationScore = %v, want > 0", metrics.OverallValidationScore)
	}
}

// TestValidateEnhancedModel_WithDemographicParams verifies that DemographicFit is
// computed (non-zero) when DemographicParams with epochs is provided.
func TestValidateEnhancedModel_WithDemographicParams(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	basePop := &H3AfricaPopulation{
		PopulationID: "TEST",
		AlleleFrequencies: map[string]float64{
			"rs429358": 0.20,
		},
		DemographicParams: &DemographicHistory{
			EffectivePopulationSize: []PopulationEpoch{
				{TimeAgo: 100, PopulationSize: 10000},
				{TimeAgo: 500, PopulationSize: 7000},
				{TimeAgo: 1000, PopulationSize: 5000},
			},
		},
	}
	params := &SyntheticGenParams{StructuralVarRate: 0.0008}

	metrics := ape.validateEnhancedModel(basePop, params)

	if metrics.DemographicFit <= 0 {
		t.Errorf("DemographicFit = %v, want > 0 when DemographicParams available", metrics.DemographicFit)
	}
	if metrics.DemographicFit > 1.0 {
		t.Errorf("DemographicFit = %v, want <= 1.0", metrics.DemographicFit)
	}
}

// TestValidateEnhancedModel_WithPGxVariants verifies PharmacogenomicFit computation.
// When a PGx variant ID is present in the allele frequency map, coverage = 1.0.
func TestValidateEnhancedModel_WithPGxVariants(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	basePop := &H3AfricaPopulation{
		PopulationID: "TEST",
		AlleleFrequencies: map[string]float64{
			"rs4244285": 0.155, // PGx variant present in allele freq table
		},
		PharmacogenomicVariants: []*PGxVariant{
			{VariantID: "rs4244285", GeneSymbol: "CYP2C19"},
		},
	}
	params := &SyntheticGenParams{StructuralVarRate: 0.0008}

	metrics := ape.validateEnhancedModel(basePop, params)

	// All PGx variants are covered in the allele frequency table → fit = 1.0.
	if metrics.PharmacogenomicFit != 1.0 {
		t.Errorf("PharmacogenomicFit = %v, want 1.0 (all PGx variants covered)", metrics.PharmacogenomicFit)
	}
}

// TestValidateEnhancedModel_ZeroFreqVariantsSkipped ensures zero-frequency variants
// do not cause division-by-zero (NaN) in AlleleFreqAccuracy.
func TestValidateEnhancedModel_ZeroFreqVariantsSkipped(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	basePop := &H3AfricaPopulation{
		PopulationID: "TEST",
		AlleleFrequencies: map[string]float64{
			"rs1815739": 0.000, // zero — must be skipped
			"rs429358":  0.205, // non-zero — must be included
		},
	}
	params := &SyntheticGenParams{
		StructuralVarRate: 0.0008,
		// AlleleFreqs matches the non-zero variant so accuracy can be computed.
		AlleleFreqs: map[string]float64{
			"rs429358": 0.205,
		},
	}

	metrics := ape.validateEnhancedModel(basePop, params)

	if metrics.AlleleFreqAccuracy != metrics.AlleleFreqAccuracy {
		t.Error("AlleleFreqAccuracy is NaN — zero-frequency variant not skipped")
	}
	if metrics.AlleleFreqAccuracy <= 0 {
		t.Errorf("AlleleFreqAccuracy = %v, want > 0 (non-zero variant should contribute)", metrics.AlleleFreqAccuracy)
	}
	// Only 1 non-zero variant counted.
	if metrics.ComparisonMetrics["allele_freq_variants"] != 1.0 {
		t.Errorf("allele_freq_variants = %v, want 1.0", metrics.ComparisonMetrics["allele_freq_variants"])
	}
}

// TestValidateEnhancedModel_WithStructuralVariants verifies StructuralVarAccuracy
// is computed when SVs are present.
func TestValidateEnhancedModel_WithStructuralVariants(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	basePop := &H3AfricaPopulation{
		PopulationID: "TEST",
		AlleleFrequencies: map[string]float64{
			"rs429358": 0.20,
		},
		StructuralVariants: []*StructuralVariant{
			{SVID: "SV001", Type: "DEL"},
		},
	}
	// 1 SV over 3,000,000 kb → observed rate = 1/3,000,000 ≈ 3.33e-7 per kb.
	// model rate = 0.0008 per kb → ratio ≈ 4.2e-4 → log ratio ≈ -7.8 → score ≈ 0.
	// Regardless of score, it must be in [0, 1] and must not be 0.0 from the
	// "no data" branch.
	params := &SyntheticGenParams{StructuralVarRate: 0.0008}

	metrics := ape.validateEnhancedModel(basePop, params)

	if metrics.StructuralVarAccuracy < 0 || metrics.StructuralVarAccuracy > 1.0 {
		t.Errorf("StructuralVarAccuracy = %v, want in [0, 1]", metrics.StructuralVarAccuracy)
	}
	if metrics.ComparisonMetrics["structural_variants"] != 1.0 {
		t.Errorf("structural_variants metric = %v, want 1.0", metrics.ComparisonMetrics["structural_variants"])
	}
}

// --- extractGenomicsParameters confidence regression ---

// TestExtractGenomicsParameters_ConfidenceFromData verifies that Confidence is
// derived from entry quality scores, not a hardcoded constant.
func TestExtractGenomicsParameters_ConfidenceFromData(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	basePop := &H3AfricaPopulation{
		PopulationID:      "YRI",
		AlleleFrequencies: map[string]float64{"rs429358": 0.2},
	}

	// Single entry with known quality and sample size.
	entries := []*LiteratureEntry{
		{
			QualityScore: 0.90,
			SampleSize:   500,
			ExtractedData: &ExtractedGenomicsData{
				AlleleFrequencies: map[string]map[string]float64{
					"YRI": {"rs429358": 0.20},
				},
				StatisticalSummaries: map[string]float64{},
			},
		},
	}

	combined, err := ape.extractGenomicsParameters(entries, basePop)
	if err != nil {
		t.Fatalf("extractGenomicsParameters returned error: %v", err)
	}

	// With sizeW = 500/1000 = 0.5 and quality = 0.90:
	// confidence = (0.90 * 0.5) / 0.5 = 0.90
	if combined.Confidence != 0.90 {
		t.Errorf("Confidence = %v, want 0.90 (weighted quality of single entry)", combined.Confidence)
	}
}

// TestExtractGenomicsParameters_NoEntries verifies Confidence is 0.0 when no
// literature entries are provided.
func TestExtractGenomicsParameters_NoEntries(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	basePop := &H3AfricaPopulation{
		PopulationID:      "YRI",
		AlleleFrequencies: map[string]float64{},
	}

	combined, err := ape.extractGenomicsParameters([]*LiteratureEntry{}, basePop)
	if err != nil {
		t.Fatalf("extractGenomicsParameters returned error: %v", err)
	}

	if combined.Confidence != 0.0 {
		t.Errorf("Confidence = %v, want 0.0 (no entries)", combined.Confidence)
	}
}

// TestGenerateAlleles_KnownVariants verifies real dbSNP alleles are returned.
func TestGenerateAlleles_KnownVariants(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	known := map[string][2]string{
		"rs429358":  {"T", "C"},
		"rs7412":    {"C", "T"},
		"rs1815739": {"C", "T"},
		"rs4988235": {"C", "T"},
		"rs4244285": {"G", "A"},
		"rs776746":  {"A", "G"},
	}
	for rsID, want := range known {
		ref, alt := ape.generateAlleles(rsID)
		if ref != want[0] || alt != want[1] {
			t.Errorf("%s: got ref=%s alt=%s, want ref=%s alt=%s", rsID, ref, alt, want[0], want[1])
		}
	}
}

// TestGenerateAlleles_UnknownVariant returns a valid, different REF/ALT pair.
func TestGenerateAlleles_UnknownVariant(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	valid := map[byte]bool{'A': true, 'T': true, 'G': true, 'C': true}
	for i := 0; i < 100; i++ {
		ref, alt := ape.generateAlleles("rs99999999")
		if len(ref) != 1 || len(alt) != 1 {
			t.Fatalf("expected single-base alleles, got ref=%q alt=%q", ref, alt)
		}
		if !valid[ref[0]] || !valid[alt[0]] {
			t.Errorf("non-DNA base: ref=%s alt=%s", ref, alt)
		}
		if ref == alt {
			t.Errorf("ref == alt = %s", ref)
		}
	}
}

// ── africanFSTForPop ─────────────────────────────────────────────────────────

// TestAfricanFSTForPop_KnownPopulations verifies that each catalogued sub-population
// returns its own distinct FST values (not all identical).
func TestAfricanFSTForPop_KnownPopulations(t *testing.T) {
	cases := []struct {
		popID      string
		wantEURLo  float64 // lower bound
		wantEURHi  float64 // upper bound
	}{
		{"YRI", 0.150, 0.156},
		{"LWK", 0.140, 0.146},
		{"GWD", 0.145, 0.151},
		{"MSL", 0.148, 0.153},
		{"ESN", 0.149, 0.154},
		{"AWI", 0.104, 0.111}, // admixed — lower FST_EUR
		{"ASW", 0.104, 0.111},
	}
	for _, tc := range cases {
		fst := africanFSTForPop(tc.popID)
		eurFST := fst["EUR"]
		if eurFST < tc.wantEURLo || eurFST > tc.wantEURHi {
			t.Errorf("africanFSTForPop(%q)[EUR] = %f, want in [%f, %f]",
				tc.popID, eurFST, tc.wantEURLo, tc.wantEURHi)
		}
	}
}

// TestAfricanFSTForPop_UnknownFallsBackToYRI verifies the YRI fallback.
func TestAfricanFSTForPop_UnknownFallsBackToYRI(t *testing.T) {
	yri := africanFSTForPop("YRI")
	unknown := africanFSTForPop("UNKNOWN_POP")
	if yri["EUR"] != unknown["EUR"] {
		t.Errorf("unknown pop FST_EUR = %f, want YRI value %f", unknown["EUR"], yri["EUR"])
	}
}

// TestAfricanFSTForPop_PopulationsAreDistinct checks that YRI and LWK differ.
func TestAfricanFSTForPop_PopulationsAreDistinct(t *testing.T) {
	yri := africanFSTForPop("YRI")
	lwk := africanFSTForPop("LWK")
	if yri["EUR"] == lwk["EUR"] {
		t.Errorf("YRI and LWK should have distinct FST_EUR values, both got %f", yri["EUR"])
	}
}

// TestAfricanFSTForPop_DiasporaLowerThanContinental verifies admixed (AWI) has
// lower FST_EUR than continental (YRI) due to European gene flow.
func TestAfricanFSTForPop_DiasporaLowerThanContinental(t *testing.T) {
	awi := africanFSTForPop("AWI")
	yri := africanFSTForPop("YRI")
	if awi["EUR"] >= yri["EUR"] {
		t.Errorf("AWI FST_EUR (%f) should be less than YRI FST_EUR (%f) due to admixture",
			awi["EUR"], yri["EUR"])
	}
}

// TestPickAlternateAllele_TsTvBias verifies transitions are more common than transversions.
func TestPickAlternateAllele_TsTvBias(t *testing.T) {
	ape := &AfricanPopulationEnhancer{}
	const N = 10000
	transitions := 0
	for i := 0; i < N; i++ {
		alt := ape.pickAlternateAllele('A')
		if alt == 'G' { // A→G is the only transition for ref=A
			transitions++
		}
	}
	// Expected ~51% transitions; require at least 40% (lenient statistical guard).
	ratio := float64(transitions) / N
	if ratio < 0.40 {
		t.Errorf("transition ratio for ref=A too low: %.3f (expected ~0.51)", ratio)
	}
}

// ── Tests for validateEnhancedModel (issue #81) ──────────────────────────────

func TestValidateEnhancedModel_AlleleFreqAccuracy_FromParams(t *testing.T) {
	// When params.AlleleFreqs is an exact copy of basePop.AlleleFrequencies,
	// the accuracy must be 1.0 (no discrepancy between model target and reference).
	ape := &AfricanPopulationEnhancer{}

	basePop := &H3AfricaPopulation{
		PopulationID: "YRI",
		AlleleFrequencies: map[string]float64{
			"rs429358":  0.42,
			"rs1815739": 0.52,
			"rs334":     0.15,
		},
	}
	params := &SyntheticGenParams{
		AlleleFreqs: map[string]float64{
			"rs429358":  0.42,
			"rs1815739": 0.52,
			"rs334":     0.15,
		},
		MutationRates:  make(map[string]float64),
		SelectionCoeff: make(map[string]float64),
	}

	metrics := ape.validateEnhancedModel(basePop, params)
	if metrics.AlleleFreqAccuracy != 1.0 {
		t.Errorf("expected AlleleFreqAccuracy=1.0 when params match reference, got %f",
			metrics.AlleleFreqAccuracy)
	}
	if metrics.ComparisonMetrics["allele_freq_variants"] != 3.0 {
		t.Errorf("expected 3 variants counted, got %f",
			metrics.ComparisonMetrics["allele_freq_variants"])
	}
}

func TestValidateEnhancedModel_AlleleFreqAccuracy_Discrepancy(t *testing.T) {
	// When params.AlleleFreqs differs from the reference, accuracy should be < 1.0.
	ape := &AfricanPopulationEnhancer{}

	basePop := &H3AfricaPopulation{
		PopulationID: "YRI",
		AlleleFrequencies: map[string]float64{
			"rs429358": 0.42,
		},
	}
	params := &SyntheticGenParams{
		AlleleFreqs: map[string]float64{
			"rs429358": 0.20, // 52% relative error → accuracy = 0.48
		},
		MutationRates:  make(map[string]float64),
		SelectionCoeff: make(map[string]float64),
	}

	metrics := ape.validateEnhancedModel(basePop, params)
	if metrics.AlleleFreqAccuracy >= 1.0 {
		t.Errorf("expected AlleleFreqAccuracy < 1.0 for discrepant params, got %f",
			metrics.AlleleFreqAccuracy)
	}
}

func TestValidateEnhancedModel_AlleleFreqAccuracy_MissingParams(t *testing.T) {
	// If params.AlleleFreqs is nil/empty, no variants are matched and
	// AlleleFreqAccuracy stays 0.0 (the zero value — not the circular ~1.0 stub).
	ape := &AfricanPopulationEnhancer{}

	basePop := &H3AfricaPopulation{
		PopulationID: "YRI",
		AlleleFrequencies: map[string]float64{
			"rs429358": 0.42,
		},
	}
	params := &SyntheticGenParams{
		AlleleFreqs:    nil, // missing
		MutationRates:  make(map[string]float64),
		SelectionCoeff: make(map[string]float64),
	}

	metrics := ape.validateEnhancedModel(basePop, params)
	if metrics.AlleleFreqAccuracy != 0.0 {
		t.Errorf("expected AlleleFreqAccuracy=0.0 when params.AlleleFreqs is nil, got %f",
			metrics.AlleleFreqAccuracy)
	}
	if metrics.ComparisonMetrics["allele_freq_variants"] != 0.0 {
		t.Errorf("expected 0 variants counted, got %f",
			metrics.ComparisonMetrics["allele_freq_variants"])
	}
}
