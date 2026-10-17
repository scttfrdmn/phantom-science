//go:build !release

package diseases_test

import (
	"context"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/diseases"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// --- Constructor tests ---

func TestNewClinVarIntegrator(t *testing.T) {
	cv := diseases.NewClinVarIntegrator("")
	require.NotNil(t, cv)
}

func TestNewClinVarIntegrator_WithAPIKey(t *testing.T) {
	cv := diseases.NewClinVarIntegrator("test-api-key")
	require.NotNil(t, cv)
}

// --- WellCharacterizedDiseases tests ---

func TestWellCharacterizedDiseases_NotEmpty(t *testing.T) {
	assert.NotEmpty(t, diseases.WellCharacterizedDiseases)
}

func TestWellCharacterizedDiseases_CysticFibrosis(t *testing.T) {
	cf, ok := diseases.WellCharacterizedDiseases["cystic_fibrosis"]
	require.True(t, ok)
	assert.Equal(t, "Cystic Fibrosis", cf.DiseaseName)
	assert.Contains(t, cf.GeneSymbols, "CFTR")
	assert.Equal(t, "autosomal_recessive", cf.InheritancePattern)
	assert.Greater(t, cf.CarrierFrequency, 0.0)
	assert.NotEmpty(t, cf.PopulationFrequencies)
}

func TestWellCharacterizedDiseases_HuntingtonDisease(t *testing.T) {
	hd, ok := diseases.WellCharacterizedDiseases["huntington_disease"]
	require.True(t, ok)
	assert.Equal(t, "Huntington Disease", hd.DiseaseName)
	assert.Contains(t, hd.GeneSymbols, "HTT")
	assert.Equal(t, "autosomal_dominant", hd.InheritancePattern)
	require.NotNil(t, hd.PenetranceModel)
	assert.NotEmpty(t, hd.PenetranceModel.AgeGroups)
	assert.Equal(t, "both", hd.PenetranceModel.Gender)
}

func TestWellCharacterizedDiseases_BRCAHereditary(t *testing.T) {
	brca, ok := diseases.WellCharacterizedDiseases["brca_hereditary_cancer"]
	require.True(t, ok)
	assert.Contains(t, brca.GeneSymbols, "BRCA1")
	assert.Contains(t, brca.GeneSymbols, "BRCA2")
	assert.Equal(t, "autosomal_dominant", brca.InheritancePattern)
	require.NotNil(t, brca.PenetranceModel)
	// BRCA has Ashkenazi Jewish founder variant
	ajFreq, exists := brca.PopulationFrequencies["AJ"]
	assert.True(t, exists)
	assert.Greater(t, ajFreq, 0.01)
}

func TestWellCharacterizedDiseases_PopulationFrequencies(t *testing.T) {
	for name, config := range diseases.WellCharacterizedDiseases {
		t.Run(name, func(t *testing.T) {
			assert.NotEmpty(t, config.PopulationFrequencies, "disease %s should have population frequencies", name)
			for pop, freq := range config.PopulationFrequencies {
				assert.GreaterOrEqual(t, freq, 0.0, "population %s frequency should be non-negative", pop)
				assert.LessOrEqual(t, freq, 1.0, "population %s frequency should be <= 1.0", pop)
			}
		})
	}
}

func TestWellCharacterizedDiseases_AffectedRatio(t *testing.T) {
	for name, config := range diseases.WellCharacterizedDiseases {
		t.Run(name, func(t *testing.T) {
			assert.Greater(t, config.AffectedRatio, 0.0, "disease %s should have positive affected ratio", name)
			assert.LessOrEqual(t, config.AffectedRatio, 1.0, "disease %s affected ratio should be <= 1.0", name)
		})
	}
}

// --- Data structure tests ---

func TestClinVarRecord_Construction(t *testing.T) {
	record := &diseases.ClinVarRecord{
		VariationID:          "5",
		AlleleID:             "7",
		Type:                 "single nucleotide variant",
		Name:                 "CFTR:c.1521_1523delCTT",
		GeneSymbol:           "CFTR",
		Chromosome:           "7",
		Position:             117199644,
		ReferenceAllele:      "CTT",
		AlternateAllele:      "-",
		ClinicalSignificance: "Pathogenic",
		ReviewStatus:         "criteria provided, multiple submitters, no conflicts",
	}
	assert.Equal(t, "5", record.VariationID)
	assert.Equal(t, "CFTR", record.GeneSymbol)
	assert.Equal(t, "Pathogenic", record.ClinicalSignificance)
}

func TestClinVarCondition_Construction(t *testing.T) {
	cond := diseases.ClinVarCondition{
		Name:              "Cystic fibrosis",
		ConceptID:         "C0010674",
		DatabaseName:      "MedGen",
		ModeOfInheritance: "Autosomal recessive inheritance",
		Phenotypes:        []string{"lung disease", "pancreatic insufficiency"},
	}
	assert.Equal(t, "Cystic fibrosis", cond.Name)
	assert.Equal(t, 2, len(cond.Phenotypes))
}

func TestDiseaseSimulationConfig_Construction(t *testing.T) {
	config := &diseases.DiseaseSimulationConfig{
		DiseaseName:       "Test Disease",
		GeneSymbols:      []string{"GENE1"},
		CarrierFrequency: 0.01,
		InheritancePattern: "autosomal_recessive",
		SyntheticCohortSize: 100,
		AffectedRatio:   0.25,
		PopulationFrequencies: map[string]float64{
			"EUR": 0.01,
			"AFR": 0.005,
		},
	}
	assert.Equal(t, "Test Disease", config.DiseaseName)
	assert.Equal(t, 0.01, config.CarrierFrequency)
	assert.Equal(t, 100, config.SyntheticCohortSize)
}

func TestPenetranceModel_Construction(t *testing.T) {
	model := &diseases.PenetranceModel{
		AgeGroups: []diseases.AgeGroup{
			{MinAge: 30, MaxAge: 40, Penetrance: 0.1},
			{MinAge: 40, MaxAge: 50, Penetrance: 0.5},
		},
		Gender:           "female",
		ModifyingFactors: map[string]float64{"family_history": 1.5},
	}
	assert.Equal(t, 2, len(model.AgeGroups))
	assert.Equal(t, "female", model.Gender)
	assert.Equal(t, 1.5, model.ModifyingFactors["family_history"])
}

func TestAgeGroup_Construction(t *testing.T) {
	ag := diseases.AgeGroup{
		MinAge:     30,
		MaxAge:     40,
		Penetrance: 0.1,
	}
	assert.Equal(t, 30, ag.MinAge)
	assert.Equal(t, 40, ag.MaxAge)
	assert.Equal(t, 0.1, ag.Penetrance)
}

func TestSyntheticPatient_Construction(t *testing.T) {
	patient := &diseases.SyntheticPatient{
		PatientID:       "SYN_AFFECTED_001",
		Population:      "EUR",
		Age:             45,
		Gender:          "female",
		AffectionStatus: "affected",
		Genotype:        make(map[string]*diseases.PatientVariant),
		PhenotypeProfile: []string{"symptom1", "symptom2"},
	}
	assert.Equal(t, "SYN_AFFECTED_001", patient.PatientID)
	assert.Equal(t, "affected", patient.AffectionStatus)
	assert.Equal(t, 2, len(patient.PhenotypeProfile))
}

func TestPatientVariant_Construction(t *testing.T) {
	pv := &diseases.PatientVariant{
		VariantID:            "5",
		Chromosome:           "7",
		Position:             117199644,
		ReferenceAllele:      "CTT",
		AlternateAllele:      "-",
		Zygosity:             "homozygous",
		ClinicalSignificance: "Pathogenic",
		InheritedFrom:        "maternal",
		PenetranceScore:      0.95,
	}
	assert.Equal(t, "homozygous", pv.Zygosity)
	assert.Equal(t, "maternal", pv.InheritedFrom)
}

func TestFamilyHistory_Construction(t *testing.T) {
	fh := &diseases.FamilyHistory{
		ParentalAffection: map[string]bool{
			"mother": true,
			"father": false,
		},
		SiblingAffection: []bool{true, false, false},
		ExtendedFamily:   map[string]bool{"uncle": true},
	}
	assert.True(t, fh.ParentalAffection["mother"])
	assert.False(t, fh.ParentalAffection["father"])
	assert.Equal(t, 3, len(fh.SiblingAffection))
}

func TestCohortStatistics_Construction(t *testing.T) {
	stats := &diseases.CohortStatistics{
		TotalSamples:        100,
		AffectedCount:       25,
		CarrierCount:        10,
		UnaffectedCount:     65,
		VariantFrequencies:  map[string]float64{"VAR1": 0.1},
		PopulationBreakdown: map[string]int{"EUR": 60, "AFR": 40},
		AgeDistribution:     map[string]int{"30-39": 30, "40-49": 40},
		GenderDistribution:  map[string]int{"male": 50, "female": 50},
	}
	assert.Equal(t, 100, stats.TotalSamples)
	assert.Equal(t, 25, stats.AffectedCount)
	assert.Equal(t, 25+10+65, stats.AffectedCount+stats.CarrierCount+stats.UnaffectedCount)
}

func TestValidationMetrics_Construction(t *testing.T) {
	vm := &diseases.ValidationMetrics{
		CarrierFrequencyAccuracy:   0.95,
		PenetranceModelAccuracy:    0.88,
		InheritancePatternAccuracy: 0.92,
		ClinicalRealism:            0.85,
		BiologicalConsistency:      0.90,
		OverallQualityScore:        0.90,
		PopulationFrequencyAccuracy: map[string]float64{
			"EUR": 0.95,
			"AFR": 0.88,
		},
	}
	assert.Equal(t, 0.95, vm.CarrierFrequencyAccuracy)
	assert.Equal(t, 0.90, vm.OverallQualityScore)
}

// --- GenerateDiseaseSimulation tests (no network needed when variants provided) ---

func TestGenerateDiseaseSimulation_WithPreloadedVariants(t *testing.T) {
	cv := diseases.NewClinVarIntegrator("")

	// Pre-populate with a variant so no network call is made
	cftrVariant := &diseases.ClinVarRecord{
		VariationID:          "5",
		AlleleID:             "7",
		Type:                 "deletion",
		GeneSymbol:           "CFTR",
		Chromosome:           "7",
		Position:             117199644,
		ReferenceAllele:      "CTT",
		AlternateAllele:      "-",
		ClinicalSignificance: "Pathogenic",
		Conditions: []diseases.ClinVarCondition{
			{Name: "Cystic fibrosis"},
		},
	}

	config := &diseases.DiseaseSimulationConfig{
		DiseaseName:       "Cystic Fibrosis",
		GeneSymbols:      []string{"CFTR"},
		PathogenicVariants: []*diseases.ClinVarRecord{cftrVariant},
		CarrierFrequency:  0.04,
		PopulationFrequencies: map[string]float64{
			"EUR": 0.04,
			"AFR": 0.01,
		},
		InheritancePattern:  "autosomal_recessive",
		SyntheticCohortSize: 10,
		AffectedRatio:       0.25,
	}

	result, err := cv.GenerateDiseaseSimulation(context.Background(), config)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.NotNil(t, result.Statistics)
	assert.NotNil(t, result.ValidationMetrics)
	assert.Equal(t, config, result.DiseaseConfig)
	assert.Greater(t, result.GenerationTime.Nanoseconds(), int64(0))
}

func TestGenerateDiseaseSimulation_AutosomeDominant(t *testing.T) {
	cv := diseases.NewClinVarIntegrator("")

	httVariant := &diseases.ClinVarRecord{
		VariationID:          "10",
		GeneSymbol:           "HTT",
		Chromosome:           "4",
		Position:             3074877,
		ReferenceAllele:      "(CAG)35",
		AlternateAllele:      "(CAG)40",
		ClinicalSignificance: "Pathogenic",
	}

	config := &diseases.DiseaseSimulationConfig{
		DiseaseName:       "Huntington Disease",
		GeneSymbols:      []string{"HTT"},
		PathogenicVariants: []*diseases.ClinVarRecord{httVariant},
		CarrierFrequency:  0.0001,
		PopulationFrequencies: map[string]float64{
			"EUR": 0.0001,
		},
		InheritancePattern:  "autosomal_dominant",
		SyntheticCohortSize: 8,
		AffectedRatio:       0.5,
		PenetranceModel: &diseases.PenetranceModel{
			AgeGroups: []diseases.AgeGroup{
				{MinAge: 40, MaxAge: 60, Penetrance: 0.5},
			},
			Gender: "both",
		},
	}

	result, err := cv.GenerateDiseaseSimulation(context.Background(), config)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Len(t, result.GeneratedSamples, 8)
}

func TestGenerateDiseaseSimulation_CohortStatistics(t *testing.T) {
	cv := diseases.NewClinVarIntegrator("")

	variant := &diseases.ClinVarRecord{
		VariationID:          "1",
		GeneSymbol:           "GENE1",
		Chromosome:           "1",
		Position:             100000,
		ReferenceAllele:      "A",
		AlternateAllele:      "T",
		ClinicalSignificance: "Pathogenic",
	}

	cohortSize := 20
	affectedRatio := 0.5
	config := &diseases.DiseaseSimulationConfig{
		DiseaseName:       "Test Disease",
		GeneSymbols:      []string{"GENE1"},
		PathogenicVariants: []*diseases.ClinVarRecord{variant},
		CarrierFrequency:  0.1,
		PopulationFrequencies: map[string]float64{
			"EUR": 0.5,
			"AFR": 0.5,
		},
		InheritancePattern:  "autosomal_dominant",
		SyntheticCohortSize: cohortSize,
		AffectedRatio:       affectedRatio,
	}

	result, err := cv.GenerateDiseaseSimulation(context.Background(), config)
	require.NoError(t, err)
	require.NotNil(t, result.Statistics)

	// Total should equal cohort size
	stats := result.Statistics
	assert.Equal(t, cohortSize, stats.TotalSamples)
	assert.Equal(t, cohortSize, stats.AffectedCount+stats.CarrierCount+stats.UnaffectedCount)
}

func TestGenerateDiseaseSimulation_ValidationMetrics(t *testing.T) {
	cv := diseases.NewClinVarIntegrator("")

	variant := &diseases.ClinVarRecord{
		VariationID:          "1",
		GeneSymbol:           "GENE1",
		Chromosome:           "1",
		Position:             100000,
		ReferenceAllele:      "A",
		AlternateAllele:      "T",
		ClinicalSignificance: "Pathogenic",
	}

	config := &diseases.DiseaseSimulationConfig{
		DiseaseName:       "Test Disease",
		GeneSymbols:      []string{"GENE1"},
		PathogenicVariants: []*diseases.ClinVarRecord{variant},
		CarrierFrequency:  0.1,
		PopulationFrequencies: map[string]float64{"EUR": 1.0},
		InheritancePattern:  "autosomal_dominant",
		SyntheticCohortSize: 10,
		AffectedRatio:       0.5,
	}

	result, err := cv.GenerateDiseaseSimulation(context.Background(), config)
	require.NoError(t, err)
	require.NotNil(t, result.ValidationMetrics)
	assert.NotNil(t, result.ValidationMetrics.PopulationFrequencyAccuracy)
}

func TestGenerateDiseaseSimulation_SamplesHaveRequiredFields(t *testing.T) {
	cv := diseases.NewClinVarIntegrator("")

	variant := &diseases.ClinVarRecord{
		VariationID:          "42",
		GeneSymbol:           "BRCA1",
		Chromosome:           "17",
		Position:             43124016,
		ReferenceAllele:      "AG",
		AlternateAllele:      "-",
		ClinicalSignificance: "Pathogenic",
	}

	config := &diseases.DiseaseSimulationConfig{
		DiseaseName:       "BRCA Cancer",
		GeneSymbols:      []string{"BRCA1"},
		PathogenicVariants: []*diseases.ClinVarRecord{variant},
		CarrierFrequency:  0.002,
		PopulationFrequencies: map[string]float64{"EUR": 1.0},
		InheritancePattern:  "autosomal_dominant",
		SyntheticCohortSize: 5,
		AffectedRatio:       0.4,
	}

	result, err := cv.GenerateDiseaseSimulation(context.Background(), config)
	require.NoError(t, err)

	for _, patient := range result.GeneratedSamples {
		assert.NotEmpty(t, patient.PatientID)
		assert.NotEmpty(t, patient.Population)
		assert.NotEmpty(t, patient.Gender)
		assert.Greater(t, patient.Age, 0)
		assert.Contains(t, []string{"affected", "carrier", "unaffected"}, patient.AffectionStatus)
		assert.NotNil(t, patient.Genotype)
		assert.NotNil(t, patient.FamilyHistory)
	}
}

// --- ClinVarSubmission tests ---

func TestClinVarSubmission_Construction(t *testing.T) {
	sub := diseases.ClinVarSubmission{
		Submitter:            "Test Lab",
		ClinicalSignificance: "Pathogenic",
		ReviewStatus:         "criteria provided",
	}
	assert.Equal(t, "Test Lab", sub.Submitter)
	assert.Equal(t, "Pathogenic", sub.ClinicalSignificance)
}
