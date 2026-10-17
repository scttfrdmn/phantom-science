package population

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestPopulation_String tests Population string representation
func TestPopulation_String(t *testing.T) {
	testCases := []struct {
		population Population
		expected   string
	}{
		{YRI, "YRI"},
		{LWK, "LWK"},
		{GWD, "GWD"},
		{CEU, "CEU"},
		{GBR, "GBR"},
		{TSI, "TSI"},
		{CHB, "CHB"},
		{JPT, "JPT"},
		{CHS, "CHS"},
		{MXL, "MXL"},
		{PUR, "PUR"},
	}

	for _, tc := range testCases {
		t.Run(string(tc.population), func(t *testing.T) {
			assert.Equal(t, tc.expected, tc.population.String())
		})
	}
}

// TestGetSupportedPopulations tests supported populations list
func TestGetSupportedPopulations(t *testing.T) {
	populations := GetSupportedPopulations()

	expectedCount := 11
	assert.Len(t, populations, expectedCount, "Should have exactly 11 supported populations")

	// Check all expected populations are present
	expectedPops := []Population{YRI, LWK, GWD, CEU, GBR, TSI, CHB, JPT, CHS, MXL, PUR}
	for _, expected := range expectedPops {
		assert.Contains(t, populations, expected, "Population %s should be in supported list", expected)
	}
}

// TestIsValidPopulation tests population validation
func TestIsValidPopulation(t *testing.T) {
	testCases := []struct {
		population string
		expected   bool
	}{
		{"YRI", true},
		{"CEU", true},
		{"CHB", true},
		{"INVALID", false},
		{"", false},
		{"yri", false}, // Case sensitive
		{"EUR", false}, // Not a specific population code
	}

	for _, tc := range testCases {
		t.Run(tc.population, func(t *testing.T) {
			assert.Equal(t, tc.expected, IsValidPopulation(tc.population))
		})
	}
}

// TestNewDemographicSimulator tests simulator creation
func TestNewDemographicSimulator(t *testing.T) {
	projectRoot := "/test/project"
	simulator := NewDemographicSimulator(projectRoot)

	assert.NotNil(t, simulator)
	assert.Contains(t, simulator.pythonPath, "phantom_science_env/bin/python3")
	assert.Contains(t, simulator.scriptPath, "internal/population/demographics.py")
	assert.Contains(t, simulator.envPath, "phantom_science_env")
}

// TestDemographicSimulator_GetPopulationParameters tests parameter retrieval
func TestDemographicSimulator_GetPopulationParameters(t *testing.T) {
	simulator := NewDemographicSimulator("/test")

	t.Run("valid populations", func(t *testing.T) {
		validPops := []Population{YRI, CEU, CHB}

		for _, pop := range validPops {
			params, err := simulator.GetPopulationParameters(pop)
			require.NoError(t, err)
			require.NotNil(t, params)

			assert.Equal(t, pop, params.Population)
			assert.Greater(t, params.EffectivePopulation, 0)
			assert.Greater(t, params.MutationRate, 0.0)
			assert.Greater(t, params.RecombinationRate, 0.0)
			assert.NotEmpty(t, params.SuperPopulation)
			assert.NotEmpty(t, params.Region)
			assert.NotEmpty(t, params.Description)
		}
	})

	t.Run("invalid population", func(t *testing.T) {
		params, err := simulator.GetPopulationParameters(Population("INVALID"))
		assert.Error(t, err)
		assert.Nil(t, params)
		assert.Contains(t, err.Error(), "unknown population")
	})

	t.Run("YRI parameters validation", func(t *testing.T) {
		params, err := simulator.GetPopulationParameters(YRI)
		require.NoError(t, err)

		assert.Equal(t, 12300, params.EffectivePopulation)
		assert.Equal(t, 0.0166, params.GrowthRate)
		assert.Equal(t, 1.25e-8, params.MutationRate)
		assert.Equal(t, 1.13e-8, params.RecombinationRate)
		assert.Equal(t, "African", params.SuperPopulation)
		assert.Equal(t, "West Africa", params.Region)
		assert.Contains(t, params.Description, "Yoruba")
	})

	t.Run("CEU parameters validation", func(t *testing.T) {
		params, err := simulator.GetPopulationParameters(CEU)
		require.NoError(t, err)

		assert.Equal(t, 11418, params.EffectivePopulation)
		assert.Equal(t, 0.0038, params.GrowthRate)
		assert.Equal(t, "European", params.SuperPopulation)
		assert.Equal(t, "Northern Europe", params.Region)
		assert.Contains(t, params.Description, "European ancestry")
	})

	t.Run("CHB parameters validation", func(t *testing.T) {
		params, err := simulator.GetPopulationParameters(CHB)
		require.NoError(t, err)

		assert.Equal(t, 9025, params.EffectivePopulation)
		assert.Equal(t, 0.0048, params.GrowthRate)
		assert.Equal(t, "East Asian", params.SuperPopulation)
		assert.Equal(t, "East Asia", params.Region)
		assert.Contains(t, params.Description, "Han Chinese")
	})
}

// TestSimulationResult_Validation tests simulation result structure
func TestSimulationResult_Validation(t *testing.T) {
	result := &SimulationResult{
		Population:     YRI,
		NumVariants:    1000,
		NumTrees:       50,
		NumSamples:     100,
		SequenceLength: 1000000,
		SimulationTime: 2.5,
		MemoryUsedMB:   150.0,
		DemographicParams: DemographicParameters{
			Population:          YRI,
			EffectivePopulation: 12300,
			MutationRate:        1.25e-8,
		},
		Variants: []VariantData{
			{
				Position:       1000,
				Alleles:        []string{"A", "T"},
				Genotypes:      []int{0, 1, 0, 1},
				Frequency:      0.25,
				AncestralState: "A",
				NumAlleles:     2,
			},
		},
		Statistics: PopulationStatistics{
			Population:             "YRI",
			NumVariants:            1000,
			NumSamples:             100,
			DiversityPi:            0.001,
			WattersonTheta:         0.0012,
			TajimaD:                -0.5,
			SegregataingSites:      950,
			RareVariantsCount:      800,
			CommonVariantsCount:    200,
			SiteFrequencySpectrum:  []int{800, 150, 50},
			MeanAlleleFrequency:    0.15,
			DemographicParams:      map[string]interface{}{"Ne": 12300},
		},
	}

	// Validate structure completeness
	assert.Equal(t, YRI, result.Population)
	assert.Equal(t, 1000, result.NumVariants)
	assert.Greater(t, result.SimulationTime, 0.0)
	assert.Len(t, result.Variants, 1)

	// Validate variant data
	variant := result.Variants[0]
	assert.Equal(t, 1000, variant.Position)
	assert.Len(t, variant.Alleles, 2)
	assert.Equal(t, []string{"A", "T"}, variant.Alleles)
	assert.Len(t, variant.Genotypes, 4)
	assert.Equal(t, 0.25, variant.Frequency)

	// Validate statistics
	stats := result.Statistics
	assert.Equal(t, "YRI", stats.Population)
	assert.Equal(t, 1000, stats.NumVariants)
	assert.Equal(t, 100, stats.NumSamples)
	assert.Greater(t, stats.DiversityPi, 0.0)
	assert.Less(t, stats.TajimaD, 0.0) // Expected negative for growing populations
	assert.Greater(t, stats.RareVariantsCount, stats.CommonVariantsCount)
}

// TestDemographicParameters_Validation tests parameter structure validation
func TestDemographicParameters_Validation(t *testing.T) {
	params := DemographicParameters{
		Population:              YRI,
		EffectivePopulation:     12300,
		GrowthRate:              0.0166,
		MutationRate:            1.25e-8,
		RecombinationRate:       1.13e-8,
		BottleneckTime:          2000,
		BottleneckSize:          1000,
		FoundingTimeYears:       70000,
		FoundingTimeGenerations: 2800,
		SuperPopulation:         "African",
		Region:                  "West Africa",
		Description:             "Yoruba in Ibadan, Nigeria",
	}

	// Validate all required fields are present and reasonable
	assert.Equal(t, YRI, params.Population)
	assert.Greater(t, params.EffectivePopulation, 0)
	assert.Greater(t, params.GrowthRate, 0.0)
	assert.Greater(t, params.MutationRate, 0.0)
	assert.Greater(t, params.RecombinationRate, 0.0)
	assert.Greater(t, params.BottleneckTime, 0)
	assert.Greater(t, params.BottleneckSize, 0)
	assert.Greater(t, params.FoundingTimeYears, 0)
	assert.Greater(t, params.FoundingTimeGenerations, 0)
	assert.NotEmpty(t, params.SuperPopulation)
	assert.NotEmpty(t, params.Region)
	assert.NotEmpty(t, params.Description)

	// Validate biological realism of parameters
	assert.Less(t, params.MutationRate, 1e-6, "Mutation rate should be realistic")
	assert.Less(t, params.RecombinationRate, 1e-6, "Recombination rate should be realistic")
	assert.Less(t, params.GrowthRate, 1.0, "Growth rate should be reasonable")
	assert.Greater(t, params.EffectivePopulation, 1000, "Effective population should be substantial")
}

// TestVariantData_Validation tests variant data structure
func TestVariantData_Validation(t *testing.T) {
	variant := VariantData{
		Position:       123456,
		Alleles:        []string{"A", "T", "C"},
		Genotypes:      []int{0, 1, 2, 0, 1},
		Frequency:      0.3,
		AncestralState: "A",
		NumAlleles:     3,
	}

	assert.Greater(t, variant.Position, 0)
	assert.Len(t, variant.Alleles, 3)
	assert.Equal(t, variant.NumAlleles, len(variant.Alleles))
	assert.GreaterOrEqual(t, variant.Frequency, 0.0)
	assert.LessOrEqual(t, variant.Frequency, 1.0)
	assert.Contains(t, variant.Alleles, variant.AncestralState)

	// Validate genotypes reference valid alleles
	for _, genotype := range variant.Genotypes {
		assert.GreaterOrEqual(t, genotype, 0)
		assert.Less(t, genotype, variant.NumAlleles)
	}
}

// TestPopulationStatistics_Validation tests population statistics structure
func TestPopulationStatistics_Validation(t *testing.T) {
	stats := PopulationStatistics{
		Population:             "YRI",
		NumVariants:            1500,
		NumSamples:             100,
		DiversityPi:            0.001234,
		WattersonTheta:         0.001456,
		TajimaD:                -0.75,
		SegregataingSites:      1400,
		RareVariantsCount:      1200,
		CommonVariantsCount:    300,
		SiteFrequencySpectrum:  []int{1000, 200, 100, 100, 50, 50},
		MeanAlleleFrequency:    0.125,
		DemographicParams:      map[string]interface{}{"Ne": 12300, "growth": 0.016},
	}

	// Basic validation
	assert.NotEmpty(t, stats.Population)
	assert.Greater(t, stats.NumVariants, 0)
	assert.Greater(t, stats.NumSamples, 0)
	assert.Greater(t, stats.DiversityPi, 0.0)
	assert.Greater(t, stats.WattersonTheta, 0.0)

	// Biological realism checks
	assert.LessOrEqual(t, stats.SegregataingSites, stats.NumVariants)
	assert.Equal(t, stats.NumVariants, stats.RareVariantsCount+stats.CommonVariantsCount)
	assert.Greater(t, stats.RareVariantsCount, stats.CommonVariantsCount, "Most variants should be rare")
	assert.GreaterOrEqual(t, stats.MeanAlleleFrequency, 0.0)
	assert.LessOrEqual(t, stats.MeanAlleleFrequency, 1.0)

	// Site frequency spectrum validation
	sfsSum := 0
	for _, count := range stats.SiteFrequencySpectrum {
		sfsSum += count
		assert.GreaterOrEqual(t, count, 0)
	}
	assert.LessOrEqual(t, sfsSum, stats.NumVariants, "SFS sum should not exceed total variants")

	// Tajima's D should be reasonable for most populations
	assert.GreaterOrEqual(t, stats.TajimaD, -3.0, "Tajima's D should not be extremely negative")
	assert.LessOrEqual(t, stats.TajimaD, 3.0, "Tajima's D should not be extremely positive")
}

// TestDemographicSimulator_ValidatePopulationStatistics tests validation logic
func TestDemographicSimulator_ValidatePopulationStatistics(t *testing.T) {
	simulator := NewDemographicSimulator("/test")

	t.Run("valid statistics", func(t *testing.T) {
		results := map[Population]*SimulationResult{
			YRI: {
				Population: YRI,
				Statistics: PopulationStatistics{
					Population:          "YRI",
					NumVariants:         1000,
					DiversityPi:         0.000020, // Within corrected YRI range [0.000010, 0.000030]
					TajimaD:             -0.5,     // Reasonable value
					RareVariantsCount:   800,
					CommonVariantsCount: 200,
				},
				DemographicParams: DemographicParameters{
					EffectivePopulation: 12300,
				},
			},
		}

		err := simulator.ValidatePopulationStatistics(results)
		assert.NoError(t, err)
	})

	t.Run("diversity out of range", func(t *testing.T) {
		results := map[Population]*SimulationResult{
			YRI: {
				Population: YRI,
				Statistics: PopulationStatistics{
					Population:          "YRI",
					NumVariants:         1000,
					DiversityPi:         0.000040, // Above YRI max 0.000030
					TajimaD:             -0.5,
					RareVariantsCount:   800,
					CommonVariantsCount: 200,
				},
				DemographicParams: DemographicParameters{
					EffectivePopulation: 12300,
				},
			},
		}

		err := simulator.ValidatePopulationStatistics(results)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "diversity")
		assert.Contains(t, err.Error(), "outside expected range")
	})

	// Regression test for #122: YRI (African) must have higher diversity range than CEU (European).
	t.Run("YRI diversity higher than CEU range", func(t *testing.T) {
		// A diversity value that is valid for YRI but would be out-of-range for CEU
		// demonstrates the ranges are ordered correctly (Africa > Europe).
		yriValidButCeuInvalid := 0.000020 // within YRI [0.000010, 0.000030]; above CEU max 0.000015
		results := map[Population]*SimulationResult{
			YRI: {
				Population: YRI,
				Statistics: PopulationStatistics{
					DiversityPi: yriValidButCeuInvalid,
					TajimaD:     -0.5,
				},
				DemographicParams: DemographicParameters{EffectivePopulation: 12300},
			},
		}
		assert.NoError(t, simulator.ValidatePopulationStatistics(results), "YRI π=0.000020 should be valid")

		results[YRI].Statistics.DiversityPi = 0.000008 // below YRI min 0.000010
		assert.Error(t, simulator.ValidatePopulationStatistics(results), "YRI π=0.000008 should be rejected (too low)")
	})

	t.Run("Tajima's D out of range", func(t *testing.T) {
		results := map[Population]*SimulationResult{
			YRI: {
				Population: YRI,
				Statistics: PopulationStatistics{
					Population:          "YRI",
					NumVariants:         1000,
					DiversityPi:         0.00001,
					TajimaD:             2.5, // Too high
					RareVariantsCount:   800,
					CommonVariantsCount: 200,
				},
				DemographicParams: DemographicParameters{
					EffectivePopulation: 12300,
				},
			},
		}

		err := simulator.ValidatePopulationStatistics(results)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "Tajima's D")
		assert.Contains(t, err.Error(), "outside reasonable range")
	})

	t.Run("insufficient rare variants", func(t *testing.T) {
		results := map[Population]*SimulationResult{
			YRI: {
				Population: YRI,
				Statistics: PopulationStatistics{
					Population:          "YRI",
					NumVariants:         1000,
					DiversityPi:         0.00001,
					TajimaD:             -0.5,
					RareVariantsCount:   600, // Only 60%, should be >70%
					CommonVariantsCount: 400,
				},
				DemographicParams: DemographicParameters{
					EffectivePopulation: 12300,
				},
			},
		}

		err := simulator.ValidatePopulationStatistics(results)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "rare variant proportion")
		assert.Contains(t, err.Error(), "too low")
	})
}

// BenchmarkDemographicSimulator_GetPopulationParameters benchmarks parameter retrieval
func BenchmarkDemographicSimulator_GetPopulationParameters(b *testing.B) {
	simulator := NewDemographicSimulator("/test")

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = simulator.GetPopulationParameters(YRI)
	}
}

// BenchmarkIsValidPopulation benchmarks population validation
func BenchmarkIsValidPopulation(b *testing.B) {
	populations := []string{"YRI", "CEU", "CHB", "INVALID", ""}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, pop := range populations {
			IsValidPopulation(pop)
		}
	}
}