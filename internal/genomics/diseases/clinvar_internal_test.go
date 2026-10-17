//go:build !release

package diseases

import (
	"context"
	"math"
	"net/http"
	"net/http/httptest"
	"testing"
	"time"
)

// TestValidateCohortRealism_EmptyPatients_NoNaN verifies that an empty patients
// slice returns zero-valued metrics without NaN (#113).
func TestValidateCohortRealism_EmptyPatients_NoNaN(t *testing.T) {
	cv := NewClinVarIntegrator("")
	config := &DiseaseSimulationConfig{
		PopulationFrequencies: map[string]float64{"EUR": 0.3, "AFR": 0.2},
		CarrierFrequency:      0.05,
	}
	stats := &CohortStatistics{
		PopulationBreakdown: map[string]int{"EUR": 0},
	}

	metrics := cv.validateCohortRealism(config, nil, stats)

	if math.IsNaN(metrics.CarrierFrequencyAccuracy) {
		t.Error("CarrierFrequencyAccuracy is NaN for empty patients")
	}
	for pop, acc := range metrics.PopulationFrequencyAccuracy {
		if math.IsNaN(acc) {
			t.Errorf("PopulationFrequencyAccuracy[%q] is NaN for empty patients", pop)
		}
	}
	if math.IsNaN(metrics.OverallQualityScore) {
		t.Error("OverallQualityScore is NaN for empty patients")
	}
}

// TestGenerateAgeForDisease_NilModel falls back to the [18,80] uniform range.
func TestGenerateAgeForDisease_NilModel(t *testing.T) {
	cv := NewClinVarIntegrator("")
	config := &DiseaseSimulationConfig{PenetranceModel: nil}
	for i := 0; i < 100; i++ {
		age := cv.generateAgeForDisease(config)
		if age < 18 || age > 80 {
			t.Errorf("age %d out of [18,80] for nil model", age)
		}
	}
}

// TestGenerateAgeForDisease_PenetranceModel checks that ages fall within the
// declared AgeGroup ranges and that groups with higher incremental penetrance
// are more often selected (statistical smoke test over many draws).
func TestGenerateAgeForDisease_PenetranceModel(t *testing.T) {
	cv := NewClinVarIntegrator("")
	config := &DiseaseSimulationConfig{
		PenetranceModel: &PenetranceModel{
			AgeGroups: []AgeGroup{
				{MinAge: 30, MaxAge: 40, Penetrance: 0.1},  // incremental weight 0.1
				{MinAge: 40, MaxAge: 50, Penetrance: 0.5},  // incremental weight 0.4
				{MinAge: 50, MaxAge: 60, Penetrance: 0.8},  // incremental weight 0.3
			},
		},
	}

	counts := make([]int, 3)
	const N = 1000
	for i := 0; i < N; i++ {
		age := cv.generateAgeForDisease(config)
		switch {
		case age >= 30 && age <= 40:
			counts[0]++
		case age > 40 && age <= 50:
			counts[1]++
		case age > 50 && age <= 60:
			counts[2]++
		default:
			t.Errorf("age %d outside any declared group", age)
		}
	}

	// Group 1 (40–50, weight 0.4) should be most common; group 0 least.
	if counts[0] >= counts[1] {
		t.Errorf("expected group[1] > group[0], got %v", counts)
	}
}

// TestGenerateAgeForDisease_EmptyAgeGroups falls back to uniform [18,80].
func TestGenerateAgeForDisease_EmptyAgeGroups(t *testing.T) {
	cv := NewClinVarIntegrator("")
	config := &DiseaseSimulationConfig{
		PenetranceModel: &PenetranceModel{AgeGroups: []AgeGroup{}},
	}
	age := cv.generateAgeForDisease(config)
	if age < 18 || age > 80 {
		t.Errorf("age %d out of [18,80] for empty groups", age)
	}
}

// TestFilterPathogenicVariants tests the internal filterPathogenicVariants function
func TestFilterPathogenicVariants(t *testing.T) {
	cv := NewClinVarIntegrator("")

	variants := []*ClinVarRecord{
		{VariationID: "1", ClinicalSignificance: "Pathogenic"},
		{VariationID: "2", ClinicalSignificance: "Likely pathogenic"},
		{VariationID: "3", ClinicalSignificance: "Benign"},
		{VariationID: "4", ClinicalSignificance: "Uncertain significance"},
		{VariationID: "5", ClinicalSignificance: "Pathogenic/Likely pathogenic"},
		{VariationID: "6", ClinicalSignificance: "Likely benign"},
	}

	result := cv.filterPathogenicVariants(variants)

	if len(result) != 3 {
		t.Errorf("Expected 3 pathogenic variants, got %d", len(result))
	}

	// Verify the right ones were selected
	ids := make(map[string]bool)
	for _, v := range result {
		ids[v.VariationID] = true
	}
	if !ids["1"] {
		t.Error("Expected Pathogenic variant to be included")
	}
	if !ids["2"] {
		t.Error("Expected Likely pathogenic variant to be included")
	}
	if !ids["5"] {
		t.Error("Expected Pathogenic/Likely pathogenic variant to be included")
	}
	if ids["3"] {
		t.Error("Benign variant should not be included")
	}
	if ids["4"] {
		t.Error("Uncertain significance variant should not be included")
	}
}

// TestFilterPathogenicVariants_Empty tests filtering empty list
func TestFilterPathogenicVariants_Empty(t *testing.T) {
	cv := NewClinVarIntegrator("")
	result := cv.filterPathogenicVariants(nil)
	if len(result) != 0 {
		t.Errorf("Expected 0 variants for nil input, got %d", len(result))
	}
}

// TestParseClinVarResponse_CFTR tests parsing for CFTR gene
func TestParseClinVarResponse_CFTR(t *testing.T) {
	cv := NewClinVarIntegrator("")
	variants, err := cv.parseClinVarResponse([]byte(`{}`), "CFTR")
	if err != nil {
		t.Fatalf("Expected no error, got: %v", err)
	}
	if len(variants) == 0 {
		t.Error("Expected CFTR variants to be returned")
	}
	if variants[0].GeneSymbol != "CFTR" {
		t.Errorf("Expected GeneSymbol CFTR, got %s", variants[0].GeneSymbol)
	}
	if variants[0].ClinicalSignificance != "Pathogenic" {
		t.Errorf("Expected Pathogenic significance, got %s", variants[0].ClinicalSignificance)
	}
}

// TestParseClinVarResponse_HTT tests parsing for HTT gene
func TestParseClinVarResponse_HTT(t *testing.T) {
	cv := NewClinVarIntegrator("")
	variants, err := cv.parseClinVarResponse([]byte(`{}`), "HTT")
	if err != nil {
		t.Fatalf("Expected no error, got: %v", err)
	}
	if len(variants) == 0 {
		t.Error("Expected HTT variants to be returned")
	}
	if variants[0].GeneSymbol != "HTT" {
		t.Errorf("Expected GeneSymbol HTT, got %s", variants[0].GeneSymbol)
	}
}

// TestParseClinVarResponse_BRCA1 tests parsing for BRCA1 gene
func TestParseClinVarResponse_BRCA1(t *testing.T) {
	cv := NewClinVarIntegrator("")
	variants, err := cv.parseClinVarResponse([]byte(`{}`), "BRCA1")
	if err != nil {
		t.Fatalf("Expected no error, got: %v", err)
	}
	if len(variants) == 0 {
		t.Error("Expected BRCA1 variants to be returned")
	}
	if variants[0].GeneSymbol != "BRCA1" {
		t.Errorf("Expected GeneSymbol BRCA1, got %s", variants[0].GeneSymbol)
	}
}

// TestParseClinVarResponse_Unknown tests parsing for unknown gene
func TestParseClinVarResponse_Unknown(t *testing.T) {
	cv := NewClinVarIntegrator("")
	variants, err := cv.parseClinVarResponse([]byte(`{}`), "UNKNOWN_GENE")
	if err != nil {
		t.Fatalf("Expected no error for unknown gene, got: %v", err)
	}
	if len(variants) != 0 {
		t.Errorf("Expected 0 variants for unknown gene, got %d", len(variants))
	}
}

// TestFetchVariantsForGene_Cache tests the cache mechanism
func TestFetchVariantsForGene_Cache(t *testing.T) {
	cv := NewClinVarIntegrator("")

	// Pre-populate cache
	cachedVariant := &ClinVarRecord{
		VariationID:          "cached-1",
		GeneSymbol:           "CACHED_GENE",
		ClinicalSignificance: "Pathogenic",
	}
	cv.cache["gene_CACHED_GENE"] = cachedVariant

	// Fetching the cached gene should return from cache without network call
	variants, err := cv.fetchVariantsForGene(context.Background(), "CACHED_GENE")
	if err != nil {
		t.Fatalf("Expected no error, got: %v", err)
	}
	if len(variants) != 1 {
		t.Errorf("Expected 1 cached variant, got %d", len(variants))
	}
	if variants[0].VariationID != "cached-1" {
		t.Errorf("Expected cached variant ID 'cached-1', got %s", variants[0].VariationID)
	}
}

// TestFetchVariantsForGene_NetworkError tests handling of connection refused
func TestFetchVariantsForGene_NetworkError(t *testing.T) {
	cv := &ClinVarIntegrator{
		// Use a port that nothing is listening on
		baseURL:   "http://127.0.0.1:1/",
		client:    &http.Client{Timeout: 100 * time.Millisecond},
		apiKey:    "",
		rateLimit: 0,
		cache:     make(map[string]*ClinVarRecord),
	}

	_, err := cv.fetchVariantsForGene(context.Background(), "CFTR")
	if err == nil {
		t.Error("Expected error for connection refused")
	}
}

// TestFetchPathogenicVariants_NetworkError tests FetchPathogenicVariants with network failure
func TestFetchPathogenicVariants_NetworkError(t *testing.T) {
	cv := &ClinVarIntegrator{
		baseURL:   "http://127.0.0.1:1/",
		client:    &http.Client{Timeout: 100 * time.Millisecond},
		apiKey:    "",
		rateLimit: 0,
		cache:     make(map[string]*ClinVarRecord),
	}

	_, err := cv.FetchPathogenicVariants(context.Background(), []string{"CFTR"})
	if err == nil {
		t.Error("Expected error for network failure")
	}
}

// TestFetchPathogenicVariants_EmptyGenes tests with no gene symbols
func TestFetchPathogenicVariants_EmptyGenes(t *testing.T) {
	cv := NewClinVarIntegrator("")
	variants, err := cv.FetchPathogenicVariants(context.Background(), []string{})
	if err != nil {
		t.Fatalf("Expected no error for empty genes, got: %v", err)
	}
	if len(variants) != 0 {
		t.Errorf("Expected 0 variants for empty genes, got %d", len(variants))
	}
}

// TestGenerateDiseaseGenotype_EmptyVariants covers the error path when no variants are available
func TestGenerateDiseaseGenotype_EmptyVariants(t *testing.T) {
	cv := NewClinVarIntegrator("")
	patient := &SyntheticPatient{
		Genotype: make(map[string]*PatientVariant),
	}
	config := &DiseaseSimulationConfig{
		PathogenicVariants: []*ClinVarRecord{}, // empty
	}
	err := cv.generateDiseaseGenotype(patient, config, true)
	if err == nil {
		t.Error("expected error when PathogenicVariants is empty")
	}
}

// TestSelectVariantByFrequency_SingleVariant returns the only variant directly.
func TestSelectVariantByFrequency_SingleVariant(t *testing.T) {
	cv := NewClinVarIntegrator("")
	v := &ClinVarRecord{VariationID: "v1", GeneSymbol: "CFTR"}
	got := cv.selectVariantByFrequency([]*ClinVarRecord{v}, "EUR")
	if got.VariationID != "v1" {
		t.Errorf("single variant: got %q, want %q", got.VariationID, "v1")
	}
}

// TestSelectVariantByFrequency_FrequencyWeighted verifies that a variant with a much
// higher allele frequency in the patient's population is selected most often.
func TestSelectVariantByFrequency_FrequencyWeighted(t *testing.T) {
	cv := NewClinVarIntegrator("")
	// variant "rare" has frequency 0.01 and "common" has 0.99 in EUR.
	rare := &ClinVarRecord{
		VariationID: "rare",
		Frequencies: map[string]float64{"EUR": 0.01},
	}
	common := &ClinVarRecord{
		VariationID: "common",
		Frequencies: map[string]float64{"EUR": 0.99},
	}

	counts := map[string]int{}
	const N = 1000
	for i := 0; i < N; i++ {
		v := cv.selectVariantByFrequency([]*ClinVarRecord{rare, common}, "EUR")
		counts[v.VariationID]++
	}
	// "common" should be selected far more often than "rare".
	if counts["common"] <= counts["rare"] {
		t.Errorf("expected common > rare, got common=%d rare=%d", counts["common"], counts["rare"])
	}
}

// TestSelectVariantByFrequency_NoFrequencyData falls back to equal weighting.
func TestSelectVariantByFrequency_NoFrequencyData(t *testing.T) {
	cv := NewClinVarIntegrator("")
	v1 := &ClinVarRecord{VariationID: "v1"} // no Frequencies map
	v2 := &ClinVarRecord{VariationID: "v2", Frequencies: map[string]float64{}}
	v3 := &ClinVarRecord{VariationID: "v3", Frequencies: map[string]float64{"AFR": 0.1}}

	// All three should be selectable; no panic.
	counts := map[string]int{}
	for i := 0; i < 300; i++ {
		v := cv.selectVariantByFrequency([]*ClinVarRecord{v1, v2, v3}, "EUR")
		counts[v.VariationID]++
	}
	// All variants should appear at least once over 300 draws.
	for _, id := range []string{"v1", "v2", "v3"} {
		if counts[id] == 0 {
			t.Errorf("variant %q never selected with equal weights (counts=%v)", id, counts)
		}
	}
}

// TestCalculateCohortStatistics_WithCarriers covers the "carrier" case in the switch
func TestCalculateCohortStatistics_WithCarriers(t *testing.T) {
	cv := NewClinVarIntegrator("")

	patients := []*SyntheticPatient{
		{AffectionStatus: "affected", Population: "EUR", Gender: "female", Age: 35},
		{AffectionStatus: "carrier", Population: "AFR", Gender: "female", Age: 42},
		{AffectionStatus: "unaffected", Population: "EUR", Gender: "female", Age: 28},
	}
	stats := cv.calculateCohortStatistics(patients)

	if stats.AffectedCount != 1 {
		t.Errorf("expected 1 affected, got %d", stats.AffectedCount)
	}
	if stats.CarrierCount != 1 {
		t.Errorf("expected 1 carrier, got %d", stats.CarrierCount)
	}
	if stats.UnaffectedCount != 1 {
		t.Errorf("expected 1 unaffected, got %d", stats.UnaffectedCount)
	}
}

// TestFetchVariantsForGene_SuccessPath covers the HTTP success path in fetchVariantsForGene
func TestFetchVariantsForGene_SuccessPath(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Type", "application/json")
		w.WriteHeader(http.StatusOK)
		_, _ = w.Write([]byte(`{}`))
	}))
	defer srv.Close()

	cv := &ClinVarIntegrator{
		baseURL:   srv.URL + "/",
		client:    srv.Client(),
		apiKey:    "",
		rateLimit: 0,
		cache:     make(map[string]*ClinVarRecord),
	}

	variants, err := cv.fetchVariantsForGene(context.Background(), "CFTR")
	if err != nil {
		t.Errorf("expected success for mock server, got error: %v", err)
	}
	// CFTR is a known gene so parseClinVarResponse returns variants even for empty JSON
	_ = variants
}

// TestFetchVariantsForGene_NonOKStatus covers the non-200 HTTP status path
func TestFetchVariantsForGene_NonOKStatus(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.WriteHeader(http.StatusBadRequest)
	}))
	defer srv.Close()

	cv := &ClinVarIntegrator{
		baseURL:   srv.URL + "/",
		client:    srv.Client(),
		apiKey:    "",
		rateLimit: 0,
		cache:     make(map[string]*ClinVarRecord),
	}

	_, err := cv.fetchVariantsForGene(context.Background(), "CFTR")
	if err == nil {
		t.Error("expected error for non-OK HTTP status")
	}
}
