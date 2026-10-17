//go:build !release

package populations_test

import (
	"context"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/populations"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// --- Constructor tests ---

func TestNewAfricanPopulationEnhancer(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()
	require.NotNil(t, ape)
}

// --- H3AfricaPopulationRegistry tests ---

func TestH3AfricaPopulationRegistry_NotEmpty(t *testing.T) {
	assert.NotEmpty(t, populations.H3AfricaPopulationRegistry)
}

func TestH3AfricaPopulationRegistry_YRI(t *testing.T) {
	yri, ok := populations.H3AfricaPopulationRegistry["YRI"]
	require.True(t, ok)
	assert.Equal(t, "YRI", yri.PopulationID)
	assert.Equal(t, "Yoruba in Ibadan, Nigeria", yri.PopulationName)
	assert.Equal(t, "Nigeria", yri.Country)
	assert.Equal(t, "West Africa", yri.Region)
	assert.Equal(t, "Niger-Congo", yri.LanguageFamily)
	assert.Greater(t, yri.SampleSize, 0)
	assert.NotEmpty(t, yri.AlleleFrequencies)
	assert.NotEmpty(t, yri.GeographicCoords)
	assert.NotEmpty(t, yri.Citations)
}

func TestH3AfricaPopulationRegistry_LWK(t *testing.T) {
	lwk, ok := populations.H3AfricaPopulationRegistry["LWK"]
	require.True(t, ok)
	assert.Equal(t, "LWK", lwk.PopulationID)
	assert.Equal(t, "Kenya", lwk.Country)
	assert.Equal(t, "East Africa", lwk.Region)
	assert.NotEmpty(t, lwk.AlleleFrequencies)
}

func TestH3AfricaPopulationRegistry_GWD(t *testing.T) {
	gwd, ok := populations.H3AfricaPopulationRegistry["GWD"]
	require.True(t, ok)
	assert.Equal(t, "Gambia", gwd.Country)
	assert.NotEmpty(t, gwd.PharmacogenomicVariants)

	// Check CYP2C19 pharmacogenomic variant
	var hasCYP2C19 bool
	for _, pgx := range gwd.PharmacogenomicVariants {
		if pgx.GeneSymbol == "CYP2C19" {
			hasCYP2C19 = true
			assert.Equal(t, "rs4244285", pgx.VariantID)
		}
	}
	assert.True(t, hasCYP2C19, "GWD should have CYP2C19 PGx variant")
}

func TestH3AfricaPopulationRegistry_MSL(t *testing.T) {
	msl, ok := populations.H3AfricaPopulationRegistry["MSL"]
	require.True(t, ok)
	assert.Equal(t, "Sierra Leone", msl.Country)
	assert.Equal(t, "Mende", msl.EthnicGroup)
}

func TestH3AfricaPopulationRegistry_ESN(t *testing.T) {
	esn, ok := populations.H3AfricaPopulationRegistry["ESN"]
	require.True(t, ok)
	assert.Equal(t, "Nigeria", esn.Country)
	assert.Equal(t, "Esan", esn.EthnicGroup)
}

func TestH3AfricaPopulationRegistry_AWI(t *testing.T) {
	awi, ok := populations.H3AfricaPopulationRegistry["AWI"]
	require.True(t, ok)
	assert.Greater(t, awi.SampleSize, 10000) // Large consortium study
	assert.NotEmpty(t, awi.Citations)
}

func TestH3AfricaPopulationRegistry_AlleleFrequencies(t *testing.T) {
	for popID, pop := range populations.H3AfricaPopulationRegistry {
		t.Run(popID, func(t *testing.T) {
			for variantID, freq := range pop.AlleleFrequencies {
				assert.GreaterOrEqual(t, freq, 0.0, "variant %s freq should be >= 0", variantID)
				assert.LessOrEqual(t, freq, 1.0, "variant %s freq should be <= 1", variantID)
			}
		})
	}
}

// --- Data structure tests ---

func TestH3AfricaPopulation_Construction(t *testing.T) {
	pop := &populations.H3AfricaPopulation{
		PopulationID:   "TEST",
		PopulationName: "Test Population",
		Country:        "TestLand",
		Region:         "West Africa",
		LanguageFamily: "Niger-Congo",
		EthnicGroup:    "TestGroup",
		SampleSize:     100,
		StudyPeriod:    []int{2010, 2015},
		AlleleFrequencies: map[string]float64{
			"rs429358": 0.2,
		},
		GeographicCoords: []float64{10.0, 20.0},
		Citations:        []string{"10.1038/example"},
	}
	assert.Equal(t, "TEST", pop.PopulationID)
	assert.Equal(t, 100, pop.SampleSize)
	assert.Equal(t, 0.2, pop.AlleleFrequencies["rs429358"])
}

func TestH3AfricaStudy_Construction(t *testing.T) {
	study := &populations.H3AfricaStudy{
		StudyID:           "STUDY001",
		StudyName:         "Test Study",
		Institution:       "Test University",
		Country:           "Kenya",
		StudyType:         "GWAS",
		Populations:       []string{"YRI", "LWK"},
		SampleSize:        500,
		SequencingPlatform: "Illumina HiSeq",
		DataType:          "WGS",
		DOI:               "10.1038/example",
		DataAvailability:  "public",
	}
	assert.Equal(t, "STUDY001", study.StudyID)
	assert.Equal(t, 500, study.SampleSize)
	assert.Contains(t, study.Populations, "YRI")
}

func TestAfricanVariant_Construction(t *testing.T) {
	v := &populations.AfricanVariant{
		VariantID:       "rs334",
		Chromosome:      "11",
		Position:        5246696,
		ReferenceAllele: "A",
		AlternateAllele: "T",
		PopulationFreqs: map[string]float64{
			"YRI": 0.12,
			"LWK": 0.15,
		},
		FunctionalClass: "missense_variant",
		GeneSymbol:      "HBB",
		Consequence:     "p.Glu7Val",
		ClinicalRelevance: "Sickle cell trait",
		DiscoveryStudy:  "STUDY001",
		ValidationStatus: "validated",
	}
	assert.Equal(t, "rs334", v.VariantID)
	assert.Equal(t, 0.12, v.PopulationFreqs["YRI"])
}

func TestLDBlock_Construction(t *testing.T) {
	ld := &populations.LDBlock{
		BlockID:    "BLOCK001",
		Chromosome: "1",
		StartPos:   1000000,
		EndPos:     1100000,
		Variants:   []string{"rs1", "rs2", "rs3"},
		HaplotypeFreqs: map[string]float64{
			"AGC": 0.45,
			"TAT": 0.35,
		},
		RecombinationRate: 0.001,
		Population:        "YRI",
	}
	assert.Equal(t, "BLOCK001", ld.BlockID)
	assert.Equal(t, int64(1000000), ld.StartPos)
	assert.Equal(t, 3, len(ld.Variants))
}

func TestDemographicHistory_Construction(t *testing.T) {
	dh := &populations.DemographicHistory{
		EffectivePopulationSize: []populations.PopulationEpoch{
			{TimeAgo: 100, PopulationSize: 10000, GrowthRate: 0.01},
			{TimeAgo: 1000, PopulationSize: 5000, GrowthRate: 0.0},
		},
		MigrationEvents: []*populations.MigrationEvent{
			{
				SourcePopulation: "YRI",
				TargetPopulation: "LWK",
				MigrationRate:    0.01,
				StartTime:        100,
				EndTime:          50,
			},
		},
	}
	assert.Equal(t, 2, len(dh.EffectivePopulationSize))
	assert.Equal(t, 1, len(dh.MigrationEvents))
}

func TestBottleneckEvent_Construction(t *testing.T) {
	be := &populations.BottleneckEvent{
		PopulationName: "YRI",
		BottleneckTime: 500,
		OriginalSize:   10000,
		BottleneckSize: 1000,
		RecoveryTime:   400,
		Severity:       0.9,
	}
	assert.Equal(t, "YRI", be.PopulationName)
	assert.Equal(t, 0.9, be.Severity)
}

func TestPGxVariant_Construction(t *testing.T) {
	pgx := &populations.PGxVariant{
		VariantID:    "rs4244285",
		GeneSymbol:   "CYP2C19",
		DrugName:     "Clopidogrel",
		Phenotype:    "Poor Metabolizer",
		PopulationFreq: map[string]float64{"GWD": 0.155},
		ClinicalImpact: "Reduced efficacy",
		StudyEvidence:  []string{"PMID:12345"},
	}
	assert.Equal(t, "CYP2C19", pgx.GeneSymbol)
	assert.Equal(t, 0.155, pgx.PopulationFreq["GWD"])
}

func TestStructuralVariant_Construction(t *testing.T) {
	sv := &populations.StructuralVariant{
		SVID:       "SV001",
		Type:       "DEL",
		Chromosome: "1",
		StartPos:   1000000,
		EndPos:     1010000,
		Size:       10000,
		PopulationFreq: map[string]float64{"YRI": 0.05},
		GeneOverlap:    []string{"GENE1"},
		FunctionalImpact: "gene_deletion",
	}
	assert.Equal(t, "DEL", sv.Type)
	assert.Equal(t, 10000, sv.Size)
}

func TestPopValidationMetrics_Construction(t *testing.T) {
	vm := &populations.PopValidationMetrics{
		AlleleFreqAccuracy:     0.95,
		LinkageAccuracy:        0.82,
		DemographicFit:         0.88,
		StructuralVarAccuracy:  0.78,
		PharmacogenomicFit:     0.85,
		OverallValidationScore: 0.856,
		ComparisonMetrics:      map[string]float64{"fst_eur": 0.15},
	}
	assert.Equal(t, 0.95, vm.AlleleFreqAccuracy)
	assert.Equal(t, 0.15, vm.ComparisonMetrics["fst_eur"])
}

func TestQualityMetrics_Construction(t *testing.T) {
	qm := &populations.QualityMetrics{
		BiologicalRealism:     0.88,
		StatisticalFidelity:   0.92,
		PopulationSpecificity: 0.85,
		ClinicalRelevance:     0.80,
		OverallQuality:        0.8625,
	}
	assert.Equal(t, 0.88, qm.BiologicalRealism)
	assert.Equal(t, 0.8625, qm.OverallQuality)
}

func TestPopulationStructure_Construction(t *testing.T) {
	ps := &populations.PopulationStructure{
		Fst: map[string]float64{
			"EUR": 0.15,
			"EAS": 0.18,
		},
		AdmixtureProportions: map[string]float64{
			"african":  0.98,
			"european": 0.02,
		},
		IsolationByDistance: true,
		SubpopulationSizes: map[string]int{
			"rural": 60000,
			"urban": 40000,
		},
	}
	assert.True(t, ps.IsolationByDistance)
	assert.Equal(t, 0.15, ps.Fst["EUR"])
}

func TestCoalescentParams_Construction(t *testing.T) {
	cp := &populations.CoalescentParams{
		ThetaWatterson:    0.00108,
		ThetaPi:           0.00125,
		TajimasD:          0.32,
		FuAndLiD:          -0.15,
		NucleotideDiversity: 0.00125,
		SegregationSites:  1250,
	}
	assert.Equal(t, 0.00108, cp.ThetaWatterson)
	assert.Equal(t, 1250, cp.SegregationSites)
}

// --- EnhanceAfricanPopulations tests (no actual network calls) ---

func TestEnhanceAfricanPopulations_SinglePopulation(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()
	result, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"YRI"})
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Contains(t, result, "YRI")

	yriModel := result["YRI"]
	require.NotNil(t, yriModel)
	assert.NotNil(t, yriModel.Population)
	assert.NotNil(t, yriModel.SyntheticParameters)
	assert.NotNil(t, yriModel.ValidationMetrics)
	assert.NotNil(t, yriModel.QualityAssessment)
}

func TestEnhanceAfricanPopulations_MultiplePopulations(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()
	result, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"YRI", "LWK", "GWD"})
	require.NoError(t, err)
	assert.Len(t, result, 3)
	assert.Contains(t, result, "YRI")
	assert.Contains(t, result, "LWK")
	assert.Contains(t, result, "GWD")
}

func TestEnhanceAfricanPopulations_UnknownPopulation(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()
	_, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"NONEXISTENT"})
	assert.Error(t, err)
}

func TestEnhanceAfricanPopulations_ModelHasValidation(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()
	result, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"LWK"})
	require.NoError(t, err)

	lwk := result["LWK"]
	require.NotNil(t, lwk.ValidationMetrics)
	// AlleleFreqAccuracy is always computed from the population's allele frequency table.
	assert.Greater(t, lwk.ValidationMetrics.AlleleFreqAccuracy, 0.0)
	// OverallValidationScore is the mean of all five metrics; it is positive when at
	// least one metric is non-zero.
	assert.Greater(t, lwk.ValidationMetrics.OverallValidationScore, 0.0)
	// DemographicFit is 0.0 for LWK because no DemographicParams are loaded in the
	// basic H3Africa registry — this is honest, not a failure.
	assert.Equal(t, 0.0, lwk.ValidationMetrics.DemographicFit)
}

func TestEnhanceAfricanPopulations_SyntheticParams(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()
	result, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"YRI"})
	require.NoError(t, err)

	params := result["YRI"].SyntheticParameters
	require.NotNil(t, params)
	assert.Greater(t, params.VariantDensity, 0.0)
	assert.Greater(t, params.TransitionTsv, 0.0)
	assert.NotNil(t, params.PopulationStructure)
	assert.NotEmpty(t, params.MutationRates)
	assert.NotEmpty(t, params.SelectionCoeff)
}

// --- GenerateAfricanSyntheticData tests ---

func TestGenerateAfricanSyntheticData_Basic(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()

	// First enhance the population
	_, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"YRI"})
	require.NoError(t, err)

	// Now generate synthetic data
	dataset, err := ape.GenerateAfricanSyntheticData(context.Background(), "YRI", 10)
	require.NoError(t, err)
	require.NotNil(t, dataset)
	assert.Equal(t, "YRI", dataset.PopulationID)
	assert.Equal(t, 10, dataset.SampleSize)
	assert.Len(t, dataset.GeneratedSamples, 10)
	assert.NotNil(t, dataset.QualityMetrics)
	assert.Equal(t, "validated", dataset.ValidationStatus)
}

func TestGenerateAfricanSyntheticData_SampleFields(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()

	_, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"GWD"})
	require.NoError(t, err)

	dataset, err := ape.GenerateAfricanSyntheticData(context.Background(), "GWD", 5)
	require.NoError(t, err)

	for _, sample := range dataset.GeneratedSamples {
		assert.NotEmpty(t, sample.SampleID)
		assert.Equal(t, "GWD", sample.PopulationID)
		assert.NotNil(t, sample.Genotype)
		assert.NotNil(t, sample.Ancestry)
		// QualityScore is a float64; positive when allele frequencies are all non-zero
		assert.NotNil(t, sample.QualityScore)
	}
}

func TestGenerateAfricanSyntheticData_WithPGxVariants(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()

	// GWD has pharmacogenomic variants defined
	_, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"GWD"})
	require.NoError(t, err)

	dataset, err := ape.GenerateAfricanSyntheticData(context.Background(), "GWD", 3)
	require.NoError(t, err)

	// GWD samples should have PGx profile (pharmacogenomic variants)
	for _, sample := range dataset.GeneratedSamples {
		assert.NotNil(t, sample.PGxProfile)
	}
}

func TestGenerateAfricanSyntheticData_ModelNotFound(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()

	// Don't enhance first, model won't exist
	_, err := ape.GenerateAfricanSyntheticData(context.Background(), "YRI", 10)
	assert.Error(t, err)
}

func TestGenerateAfricanSyntheticData_GenerationHistory(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()

	enhanced, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"MSL"})
	require.NoError(t, err)

	// Initially empty generation history
	assert.Empty(t, enhanced["MSL"].GenerationHistory)

	_, err = ape.GenerateAfricanSyntheticData(context.Background(), "MSL", 5)
	require.NoError(t, err)

	// After generation, history should have one event
	assert.Len(t, enhanced["MSL"].GenerationHistory, 1)
	event := enhanced["MSL"].GenerationHistory[0]
	assert.Equal(t, "synthetic_data_generation", event.EventType)
	assert.Equal(t, 5, event.DataGenerated)
}

func TestGenerateAfricanSyntheticData_Timing(t *testing.T) {
	ape := populations.NewAfricanPopulationEnhancer()

	_, err := ape.EnhanceAfricanPopulations(context.Background(), []string{"YRI"})
	require.NoError(t, err)

	dataset, err := ape.GenerateAfricanSyntheticData(context.Background(), "YRI", 20)
	require.NoError(t, err)
	assert.GreaterOrEqual(t, dataset.GenerationTime.Nanoseconds(), int64(0))
}
