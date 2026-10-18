//go:build !release
// +build !release

package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/diseases"
)

func TestContains_Found(t *testing.T) {
	slice := []string{"apple", "banana", "cherry"}
	if !contains(slice, "banana") {
		t.Error("contains should find 'banana' in slice")
	}
}

func TestContains_NotFound(t *testing.T) {
	slice := []string{"apple", "banana", "cherry"}
	if contains(slice, "grape") {
		t.Error("contains should not find 'grape' in slice")
	}
}

func TestContains_Empty(t *testing.T) {
	if contains([]string{}, "anything") {
		t.Error("contains should return false for empty slice")
	}
}

func TestGetDiseaseConfiguration_CysticFibrosis(t *testing.T) {
	config, err := getDiseaseConfiguration("cystic_fibrosis")
	if err != nil {
		t.Fatalf("Expected no error for cystic_fibrosis, got: %v", err)
	}
	if config == nil {
		t.Fatal("Expected non-nil config")
	}
	if config.DiseaseName != "Cystic Fibrosis" {
		t.Errorf("Expected 'Cystic Fibrosis', got '%s'", config.DiseaseName)
	}
	if len(config.GeneSymbols) == 0 {
		t.Error("Expected at least one gene symbol")
	}
	if config.InheritancePattern != "autosomal_recessive" {
		t.Errorf("Expected 'autosomal_recessive', got '%s'", config.InheritancePattern)
	}
}

func TestGetDiseaseConfiguration_HuntingtonDisease(t *testing.T) {
	config, err := getDiseaseConfiguration("huntington_disease")
	if err != nil {
		t.Fatalf("Expected no error for huntington_disease, got: %v", err)
	}
	if config.InheritancePattern != "autosomal_dominant" {
		t.Errorf("Expected 'autosomal_dominant', got '%s'", config.InheritancePattern)
	}
}

func TestGetDiseaseConfiguration_BRCAHereditary(t *testing.T) {
	config, err := getDiseaseConfiguration("brca_hereditary_cancer")
	if err != nil {
		t.Fatalf("Expected no error for brca_hereditary_cancer, got: %v", err)
	}
	if config.DiseaseName == "" {
		t.Error("Expected non-empty disease name")
	}
}

func TestGetDiseaseConfiguration_Unknown(t *testing.T) {
	_, err := getDiseaseConfiguration("unknown_disease")
	if err == nil {
		t.Error("Expected error for unknown disease type")
	}
}

func TestGetDiseaseConfiguration_Isolation(t *testing.T) {
	// Verify returned configs are independent copies
	config1, err := getDiseaseConfiguration("cystic_fibrosis")
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}
	config2, err := getDiseaseConfiguration("cystic_fibrosis")
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}

	// Modify one, should not affect the other
	config1.SyntheticCohortSize = 999
	if config2.SyntheticCohortSize == 999 {
		t.Error("Expected configurations to be independent copies")
	}
}

func TestValidateInputs_Valid(t *testing.T) {
	// Set global flag values for a valid configuration
	*diseaseType = "cystic_fibrosis"
	*cohortSize = 100
	*affectedRatio = 0.1
	*population = "EUR"
	*format = "json"
	if err := validateInputs(); err != nil {
		t.Errorf("validateInputs() unexpected error: %v", err)
	}
}

func TestValidateInputs_InvalidDisease(t *testing.T) {
	*diseaseType = "nonexistent_disease"
	*cohortSize = 100
	*affectedRatio = 0.1
	*population = "EUR"
	*format = "json"
	if err := validateInputs(); err == nil {
		t.Error("validateInputs() should fail for invalid disease type")
	}
}

func TestValidateInputs_InvalidCohortSize(t *testing.T) {
	*diseaseType = "cystic_fibrosis"
	*cohortSize = -1
	*affectedRatio = 0.1
	*population = "EUR"
	*format = "json"
	if err := validateInputs(); err == nil {
		t.Error("validateInputs() should fail for negative cohort size")
	}
}

func TestValidateInputs_InvalidAffectedRatio(t *testing.T) {
	*diseaseType = "cystic_fibrosis"
	*cohortSize = 100
	*affectedRatio = 1.5
	*population = "EUR"
	*format = "json"
	if err := validateInputs(); err == nil {
		t.Error("validateInputs() should fail for ratio > 1")
	}
}

func TestValidateInputs_InvalidPopulation(t *testing.T) {
	*diseaseType = "cystic_fibrosis"
	*cohortSize = 100
	*affectedRatio = 0.1
	*population = "INVALID_POP"
	*format = "json"
	if err := validateInputs(); err == nil {
		t.Error("validateInputs() should fail for invalid population")
	}
}

func TestValidateInputs_InvalidFormat(t *testing.T) {
	*diseaseType = "cystic_fibrosis"
	*cohortSize = 100
	*affectedRatio = 0.1
	*population = "EUR"
	*format = "csv"
	if err := validateInputs(); err == nil {
		t.Error("validateInputs() should fail for invalid format")
	}
}

// ---------- helper to build a minimal SyntheticDiseaseResult ----------

func makeSyntheticResult() *diseases.SyntheticDiseaseResult {
	return &diseases.SyntheticDiseaseResult{
		DiseaseConfig: &diseases.DiseaseSimulationConfig{
			DiseaseName:        "Test Disease",
			GeneSymbols:       []string{"GENE1"},
			InheritancePattern: "autosomal_recessive",
			PathogenicVariants: []*diseases.ClinVarRecord{},
		},
		GeneratedSamples: []*diseases.SyntheticPatient{
			{
				PatientID:       "P001",
				Gender:          "male",
				AffectionStatus: "affected",
				Genotype:        map[string]*diseases.PatientVariant{},
			},
			{
				PatientID:       "P002",
				Gender:          "female",
				AffectionStatus: "unaffected",
				Genotype:        map[string]*diseases.PatientVariant{},
			},
		},
		Statistics: &diseases.CohortStatistics{
			TotalSamples:        2,
			AffectedCount:       1,
			CarrierCount:        0,
			UnaffectedCount:     1,
			PopulationBreakdown: map[string]int{"EUR": 2},
		},
		ValidationMetrics: &diseases.ValidationMetrics{
			OverallQualityScore:         0.85,
			CarrierFrequencyAccuracy:    0.90,
			PopulationFrequencyAccuracy: map[string]float64{"EUR": 0.88},
		},
		GenerationTime: time.Second * 5,
	}
}

// ---------- saveJSONReport ----------

func TestSaveJSONReport_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "test_result.json")
	result := makeSyntheticResult()
	if err := saveJSONReport(result, path); err != nil {
		t.Fatalf("saveJSONReport returned error: %v", err)
	}
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty JSON file")
	}
}

func TestSaveJSONReport_BadPath(t *testing.T) {
	err := saveJSONReport(makeSyntheticResult(), "/nonexistent/dir/result.json")
	if err == nil {
		t.Error("expected error for bad path")
	}
}

// ---------- saveSummaryReport ----------

func TestSaveSummaryReport_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "summary.txt")
	if err := saveSummaryReport(makeSyntheticResult(), path); err != nil {
		t.Fatalf("saveSummaryReport returned error: %v", err)
	}
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty summary file")
	}
}

func TestSaveSummaryReport_BadPath(t *testing.T) {
	err := saveSummaryReport(makeSyntheticResult(), "/nonexistent/dir/summary.txt")
	if err == nil {
		t.Error("expected error for bad path")
	}
}

// ---------- saveVCFFormat ----------

func TestSaveVCFFormat_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "output.vcf")
	if err := saveVCFFormat(makeSyntheticResult(), path); err != nil {
		t.Fatalf("saveVCFFormat returned error: %v", err)
	}
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("expected VCF file to exist: %v", err)
	}
}

func TestSaveVCFFormat_BadPath(t *testing.T) {
	err := saveVCFFormat(makeSyntheticResult(), "/nonexistent/dir/output.vcf")
	if err == nil {
		t.Error("expected error for bad path")
	}
}

// ---------- savePEDFormat ----------

func TestSavePEDFormat_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "output.ped")
	if err := savePEDFormat(makeSyntheticResult(), path); err != nil {
		t.Fatalf("savePEDFormat returned error: %v", err)
	}
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("expected PED file to exist: %v", err)
	}
}

func TestSavePEDFormat_BadPath(t *testing.T) {
	err := savePEDFormat(makeSyntheticResult(), "/nonexistent/dir/output.ped")
	if err == nil {
		t.Error("expected error for bad path")
	}
}

// ---------- printSimulationSummary ----------

func TestPrintSimulationSummary_NoCrash(t *testing.T) {
	printSimulationSummary(makeSyntheticResult())
}

// ---------- customizeConfiguration ----------

func TestCustomizeConfiguration_MixedPopulation(t *testing.T) {
	config, err := getDiseaseConfiguration("cystic_fibrosis")
	if err != nil {
		t.Fatalf("getDiseaseConfiguration: %v", err)
	}
	*cohortSize = 500
	*affectedRatio = 0.2
	*population = "mixed"
	customizeConfiguration(config)
	if config.SyntheticCohortSize != 500 {
		t.Errorf("SyntheticCohortSize = %d, want 500", config.SyntheticCohortSize)
	}
	if config.AffectedRatio != 0.2 {
		t.Errorf("AffectedRatio = %f, want 0.2", config.AffectedRatio)
	}
}

func TestCustomizeConfiguration_KnownPopulation(t *testing.T) {
	config, err := getDiseaseConfiguration("cystic_fibrosis")
	if err != nil {
		t.Fatalf("getDiseaseConfiguration: %v", err)
	}
	*cohortSize = 200
	*affectedRatio = 0.1
	*population = "EUR"
	customizeConfiguration(config)
	if config.SyntheticCohortSize != 200 {
		t.Errorf("SyntheticCohortSize = %d, want 200", config.SyntheticCohortSize)
	}
	// EUR should appear in the frequencies map
	if _, ok := config.PopulationFrequencies["EUR"]; !ok {
		t.Error("expected EUR to be in PopulationFrequencies")
	}
}

func TestCustomizeConfiguration_UnknownPopulation(t *testing.T) {
	config, err := getDiseaseConfiguration("cystic_fibrosis")
	if err != nil {
		t.Fatalf("getDiseaseConfiguration: %v", err)
	}
	*cohortSize = 100
	*affectedRatio = 0.1
	*population = "UNKNOWN_POP"
	customizeConfiguration(config)
	if config.SyntheticCohortSize != 100 {
		t.Errorf("SyntheticCohortSize = %d, want 100", config.SyntheticCohortSize)
	}
	// Unknown population uses conservative carrier frequency
	if config.CarrierFrequency != 0.01 {
		t.Errorf("CarrierFrequency = %f, want 0.01 (conservative estimate)", config.CarrierFrequency)
	}
}

// ---------- saveResults ----------

func TestSaveResults_JSONFormat(t *testing.T) {
	tmpDir := t.TempDir()
	*outputDir = tmpDir
	*format = "json"
	result := makeSyntheticResult()
	if err := saveResults(result); err != nil {
		t.Fatalf("saveResults(json) returned error: %v", err)
	}
}

func TestSaveResults_VCFFormat(t *testing.T) {
	tmpDir := t.TempDir()
	*outputDir = tmpDir
	*format = "vcf"
	result := makeSyntheticResult()
	if err := saveResults(result); err != nil {
		t.Fatalf("saveResults(vcf) returned error: %v", err)
	}
}

func TestSaveResults_PEDFormat(t *testing.T) {
	tmpDir := t.TempDir()
	*outputDir = tmpDir
	*format = "ped"
	result := makeSyntheticResult()
	if err := saveResults(result); err != nil {
		t.Fatalf("saveResults(ped) returned error: %v", err)
	}
}

func TestSaveResults_BadOutputDir(t *testing.T) {
	*outputDir = "/nonexistent/directory/that/does/not/exist"
	*format = "json"
	err := saveResults(makeSyntheticResult())
	if err == nil {
		t.Error("expected error for non-existent output directory")
	}
}

// ---------- saveVCFFormat with variants and zygosity ----------

func TestSaveVCFFormat_WithVariantsAndZygosity(t *testing.T) {
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "zygosity.vcf")

	result := &diseases.SyntheticDiseaseResult{
		DiseaseConfig: &diseases.DiseaseSimulationConfig{
			DiseaseName: "Test Disease",
			PathogenicVariants: []*diseases.ClinVarRecord{
				{
					VariationID:          "VAR001",
					GeneSymbol:           "GENE1",
					Chromosome:           "7",
					Position:             117548628,
					ReferenceAllele:      "C",
					AlternateAllele:      "T",
					ClinicalSignificance: "Pathogenic",
				},
			},
		},
		GeneratedSamples: []*diseases.SyntheticPatient{
			{
				PatientID:       "P001",
				Gender:          "male",
				AffectionStatus: "affected",
				Genotype: map[string]*diseases.PatientVariant{
					"GENE1": {Zygosity: "heterozygous", ClinicalSignificance: "Pathogenic"},
				},
			},
			{
				PatientID:       "P002",
				Gender:          "female",
				AffectionStatus: "affected",
				Genotype: map[string]*diseases.PatientVariant{
					"GENE1": {Zygosity: "homozygous", ClinicalSignificance: "Pathogenic"},
				},
			},
			{
				PatientID:       "P003",
				Gender:          "male",
				AffectionStatus: "unaffected",
				Genotype:        map[string]*diseases.PatientVariant{}, // no GENE1 → 0/0
			},
		},
		Statistics:        &diseases.CohortStatistics{TotalSamples: 3},
		ValidationMetrics: &diseases.ValidationMetrics{},
	}

	if err := saveVCFFormat(result, path); err != nil {
		t.Fatalf("saveVCFFormat returned error: %v", err)
	}

	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("expected VCF file: %v", err)
	}
	content := string(data)
	if !strings.Contains(content, "0/1") {
		t.Error("expected heterozygous genotype 0/1 in VCF output")
	}
	if !strings.Contains(content, "1/1") {
		t.Error("expected homozygous genotype 1/1 in VCF output")
	}
}
