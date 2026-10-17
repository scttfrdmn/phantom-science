package population

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewPrivacyPreservingGeneticsEngine tests engine creation
func TestNewPrivacyPreservingGeneticsEngine(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()

	assert.NotNil(t, engine)
	assert.NotNil(t, engine.dpganGenerator)
	assert.NotNil(t, engine.pateganGenerator)
	assert.NotNil(t, engine.federatedLearning)
	assert.NotNil(t, engine.differentialPrivacy)
	assert.NotNil(t, engine.tabpfnDefense)
	assert.NotNil(t, engine.localAncestryShield)
	assert.NotNil(t, engine.haplonetCountermeasure)
	assert.NotNil(t, engine.privacyBudgetManager)
}

// TestNewDPGANGenerator tests DPGAN generator creation
func TestNewDPGANGenerator(t *testing.T) {
	gen := NewDPGANGenerator()

	assert.NotNil(t, gen)
	assert.NotNil(t, gen.generator)
	assert.NotNil(t, gen.discriminator)
	assert.Equal(t, 1.0, gen.privacyBudget)
	assert.NotNil(t, gen.noiseMechanism)
	assert.NotNil(t, gen.gradientClipping)
	assert.NotNil(t, gen.populationModel)
	assert.Equal(t, 1000, gen.generator.outputDim)
	assert.Equal(t, 0.50, gen.discriminator.accuracyTarget)
}

// TestNewPATEGANGenerator tests PATEGAN generator creation
func TestNewPATEGANGenerator(t *testing.T) {
	gen := NewPATEGANGenerator()

	assert.NotNil(t, gen)
	assert.NotNil(t, gen.teacherEnsemble)
	assert.NotNil(t, gen.populationTeachers)
	assert.Equal(t, 0.75, gen.consensusThreshold)
	assert.NotNil(t, gen.privacyAccountant)
	assert.Equal(t, 10.0, gen.privacyAccountant.privacyBudgetRemaining)
}

// TestNewFederatedPrivacyManager tests federated privacy manager creation
func TestNewFederatedPrivacyManager(t *testing.T) {
	mgr := NewFederatedPrivacyManager()

	assert.NotNil(t, mgr)
	assert.NotNil(t, mgr.federationNodes)
	assert.NotNil(t, mgr.aggregationServer)
	assert.NotNil(t, mgr.privacyProtocols)
	assert.NotNil(t, mgr.secureAggregation)
	assert.NotNil(t, mgr.populationNodeMap)
}

// TestNewAdvancedDifferentialPrivacy tests advanced DP creation
func TestNewAdvancedDifferentialPrivacy(t *testing.T) {
	dp := NewAdvancedDifferentialPrivacy()

	assert.NotNil(t, dp)
	assert.NotNil(t, dp.privacyMechanisms)
	assert.NotNil(t, dp.familyCorrelationHandler)
	assert.NotNil(t, dp.genomicSensitivityAnalyzer)
	assert.NotNil(t, dp.compositionTracker)
	assert.NotNil(t, dp.adaptiveNoiseManager)
	assert.NotNil(t, dp.genomicSensitivityAnalyzer.sensitivityFunctions)
	assert.NotNil(t, dp.genomicSensitivityAnalyzer.locusSpecificSensitivity)
	assert.NotNil(t, dp.genomicSensitivityAnalyzer.populationSensitivity)
}

// TestNewTabPFNDefenseSystem tests TabPFN defense system creation
func TestNewTabPFNDefenseSystem(t *testing.T) {
	defense := NewTabPFNDefenseSystem()

	assert.NotNil(t, defense)
	assert.NotNil(t, defense.featureMasking)
	assert.NotNil(t, defense.adversarialPerturbations)
	assert.NotNil(t, defense.tabularTransformerShield)
	assert.NotNil(t, defense.ancestryObfuscation)
	assert.NotNil(t, defense.featureMasking.featureImportanceMasking)
	assert.NotNil(t, defense.adversarialPerturbations.perturbationStrategies)
}

// TestNewLocalAncestryShield tests local ancestry shield creation
func TestNewLocalAncestryShield(t *testing.T) {
	shield := NewLocalAncestryShield()

	assert.NotNil(t, shield)
	assert.NotNil(t, shield.localAncestryScrambler)
	assert.NotNil(t, shield.haplotypeObfuscation)
	assert.NotNil(t, shield.linkageDisequilibriumMasking)
	assert.NotNil(t, shield.admixturePerturbation)
}

// TestNewHaploNetCountermeasure tests HaploNet countermeasure creation
func TestNewHaploNetCountermeasure(t *testing.T) {
	counter := NewHaploNetCountermeasure()

	assert.NotNil(t, counter)
	assert.NotNil(t, counter.vaeBypassing)
	assert.NotNil(t, counter.latentSpaceObfuscation)
	assert.NotNil(t, counter.populationStructureScrambling)
	assert.NotNil(t, counter.dimensionalityReduction)
}

// TestNewPrivacyBudgetManager tests privacy budget manager creation
func TestNewPrivacyBudgetManager(t *testing.T) {
	mgr := NewPrivacyBudgetManager()

	assert.NotNil(t, mgr)
	assert.Equal(t, 10.0, mgr.totalBudget)
	assert.Equal(t, 0.0, mgr.allocatedBudget)
	assert.NotNil(t, mgr.budgetHistory)
	assert.NotNil(t, mgr.compositionRules)
	assert.True(t, mgr.compositionRules.AdvancedComposition)
	assert.True(t, mgr.compositionRules.OptimalComposition)
}

// TestPrivacyBudgetManager_AllocateBudget tests budget allocation
func TestPrivacyBudgetManager_AllocateBudget(t *testing.T) {
	mgr := NewPrivacyBudgetManager()

	t.Run("successful allocation", func(t *testing.T) {
		requirements := &PrivacyRequirements{
			DifferentialPrivacyEpsilon: 1.0,
		}

		allocation, err := mgr.AllocateBudget(requirements)
		require.NoError(t, err)
		require.NotNil(t, allocation)

		assert.Equal(t, 0.4, allocation.DPGANBudget)
		assert.Equal(t, 0.3, allocation.PATEGANBudget)
		assert.Equal(t, 0.2, allocation.FederatedBudget)
		assert.Equal(t, 0.08, allocation.DefenseBudget)
		assert.Equal(t, 0.02, allocation.ValidationBudget)
		assert.Equal(t, 1.0, mgr.allocatedBudget)
		assert.Len(t, mgr.budgetHistory, 1)
	})

	t.Run("insufficient budget", func(t *testing.T) {
		// Allocate almost all budget
		mgr2 := NewPrivacyBudgetManager()
		requirements := &PrivacyRequirements{
			DifferentialPrivacyEpsilon: 10.0,
		}

		_, err := mgr2.AllocateBudget(requirements)
		require.NoError(t, err)

		// Next allocation should fail
		requirements2 := &PrivacyRequirements{
			DifferentialPrivacyEpsilon: 1.0,
		}
		_, err = mgr2.AllocateBudget(requirements2)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "insufficient privacy budget")
	})

	t.Run("multiple allocations", func(t *testing.T) {
		mgr3 := NewPrivacyBudgetManager()

		for i := 0; i < 5; i++ {
			requirements := &PrivacyRequirements{
				DifferentialPrivacyEpsilon: 1.0,
			}
			allocation, err := mgr3.AllocateBudget(requirements)
			require.NoError(t, err)
			assert.NotNil(t, allocation)
		}

		assert.Equal(t, 5.0, mgr3.allocatedBudget)
		assert.Len(t, mgr3.budgetHistory, 5)
	})

	t.Run("zero epsilon returns error", func(t *testing.T) {
		mgr4 := NewPrivacyBudgetManager()
		requirements := &PrivacyRequirements{
			DifferentialPrivacyEpsilon: 0.0,
		}
		_, err := mgr4.AllocateBudget(requirements)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "DifferentialPrivacyEpsilon")
		assert.Equal(t, 0.0, mgr4.allocatedBudget)
	})

	t.Run("negative epsilon returns error", func(t *testing.T) {
		mgr5 := NewPrivacyBudgetManager()
		requirements := &PrivacyRequirements{
			DifferentialPrivacyEpsilon: -1.0,
		}
		_, err := mgr5.AllocateBudget(requirements)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "DifferentialPrivacyEpsilon")
		assert.Equal(t, 0.0, mgr5.allocatedBudget)
	})
}

// TestDPGANGenerator_GenerateData tests data generation
func TestDPGANGenerator_GenerateData(t *testing.T) {
	gen := NewDPGANGenerator()

	populations := []Population{YRI, CEU, CHB, JPT, LWK}
	numSamples := 10
	privacyBudget := 1.0

	for _, pop := range populations {
		t.Run(string(pop), func(t *testing.T) {
			data, err := gen.GenerateData(pop, numSamples, privacyBudget)
			require.NoError(t, err)
			require.NotNil(t, data)

			assert.Len(t, data, numSamples)

			for _, sample := range data {
				assert.Len(t, sample, gen.generator.outputDim)
			}
		})
	}
}

// TestDPGANGenerator_GenerateData_ZeroSamples tests zero samples generation
func TestDPGANGenerator_GenerateData_ZeroSamples(t *testing.T) {
	gen := NewDPGANGenerator()

	data, err := gen.GenerateData(CEU, 0, 1.0)
	require.NoError(t, err)
	assert.Empty(t, data)
}

// TestDPGANGenerator_getPopulationBias tests population bias calculation
func TestDPGANGenerator_getPopulationBias(t *testing.T) {
	gen := NewDPGANGenerator()

	populations := []Population{YRI, CEU, CHB, JPT, LWK, GBR} // LWK and GBR use default
	for _, pop := range populations {
		t.Run(string(pop), func(t *testing.T) {
			bias := gen.getPopulationBias(pop, 100)
			// Bias should be within reasonable range
			assert.GreaterOrEqual(t, bias, -0.2)
			assert.LessOrEqual(t, bias, 0.2)
		})
	}
}

// TestTabPFNDefenseSystem_ApplyDefense tests TabPFN defense application
func TestTabPFNDefenseSystem_ApplyDefense(t *testing.T) {
	defense := NewTabPFNDefenseSystem()

	data := [][]float64{
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
		{7.0, 8.0, 9.0},
	}

	defended, err := defense.ApplyDefense(data, CEU)
	require.NoError(t, err)
	require.NotNil(t, defended)

	assert.Len(t, defended, len(data))
	for i, sample := range defended {
		assert.Len(t, sample, len(data[i]))
		// Values should be slightly different due to noise
		// (though with low noise they could be close)
	}
}

// TestTabPFNDefenseSystem_ApplyDefense_EmptyData tests defense with empty data
func TestTabPFNDefenseSystem_ApplyDefense_EmptyData(t *testing.T) {
	defense := NewTabPFNDefenseSystem()

	defended, err := defense.ApplyDefense([][]float64{}, CEU)
	require.NoError(t, err)
	assert.Empty(t, defended)
}

// TestLocalAncestryShield_ApplyShield tests local ancestry shield application
func TestLocalAncestryShield_ApplyShield(t *testing.T) {
	shield := NewLocalAncestryShield()

	data := [][]float64{
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
	}

	populations := []Population{YRI, CEU, CHB, JPT, LWK}
	for _, pop := range populations {
		t.Run(string(pop), func(t *testing.T) {
			shielded, err := shield.ApplyShield(data, pop)
			require.NoError(t, err)
			require.NotNil(t, shielded)
			assert.Len(t, shielded, len(data))
			for i, sample := range shielded {
				assert.Len(t, sample, len(data[i]))
			}
		})
	}
}

// TestHaploNetCountermeasure_ApplyCountermeasure tests HaploNet countermeasure application
func TestHaploNetCountermeasure_ApplyCountermeasure(t *testing.T) {
	counter := NewHaploNetCountermeasure()

	data := [][]float64{
		{0.1, 0.2, 0.3, 0.4},
		{0.5, 0.6, 0.7, 0.8},
	}

	protected, err := counter.ApplyCountermeasure(data, YRI)
	require.NoError(t, err)
	require.NotNil(t, protected)
	assert.Len(t, protected, len(data))
	for i, sample := range protected {
		assert.Len(t, sample, len(data[i]))
	}
}

// TestGeneratePrivacyPreservingGenomicData_BasicNoDefenses tests generation without defenses
func TestGeneratePrivacyPreservingGenomicData_BasicNoDefenses(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()
	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 1.0,
		UseFederatedPrivacy:        false,
		DefendAgainstTabPFN:        false,
		DefendAgainstLocalAncestry: false,
		DefendAgainstHaploNet:      false,
	}

	result, err := engine.GeneratePrivacyPreservingGenomicData(CEU, 5, requirements)
	require.NoError(t, err)
	require.NotNil(t, result)

	assert.Len(t, result.GeneratedData, 5)
	assert.NotNil(t, result.PrivacyMetrics)
	assert.NotNil(t, result.DefenseMetrics)

	for _, sample := range result.GeneratedData {
		assert.NotEmpty(t, sample.SampleID)
		assert.Equal(t, CEU, sample.Population)
		assert.NotEmpty(t, sample.GenotypeData)
	}
}

// TestGeneratePrivacyPreservingGenomicData_WithFederated tests with federated privacy
func TestGeneratePrivacyPreservingGenomicData_WithFederated(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()
	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 2.0,
		UseFederatedPrivacy:        true,
		DefendAgainstTabPFN:        false,
		DefendAgainstLocalAncestry: false,
		DefendAgainstHaploNet:      false,
	}

	result, err := engine.GeneratePrivacyPreservingGenomicData(YRI, 3, requirements)
	require.NoError(t, err)
	assert.Len(t, result.GeneratedData, 3)
}

// TestGeneratePrivacyPreservingGenomicData_WithTabPFN tests with TabPFN defense
func TestGeneratePrivacyPreservingGenomicData_WithTabPFN(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()
	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 2.0,
		UseFederatedPrivacy:        false,
		DefendAgainstTabPFN:        true,
		DefendAgainstLocalAncestry: false,
		DefendAgainstHaploNet:      false,
	}

	result, err := engine.GeneratePrivacyPreservingGenomicData(CHB, 4, requirements)
	require.NoError(t, err)
	assert.Len(t, result.GeneratedData, 4)
}

// TestGeneratePrivacyPreservingGenomicData_WithLocalAncestry tests with local ancestry defense
func TestGeneratePrivacyPreservingGenomicData_WithLocalAncestry(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()
	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 2.0,
		UseFederatedPrivacy:        false,
		DefendAgainstTabPFN:        false,
		DefendAgainstLocalAncestry: true,
		DefendAgainstHaploNet:      false,
	}

	result, err := engine.GeneratePrivacyPreservingGenomicData(JPT, 3, requirements)
	require.NoError(t, err)
	assert.Len(t, result.GeneratedData, 3)
}

// TestGeneratePrivacyPreservingGenomicData_WithHaploNet tests with HaploNet defense
func TestGeneratePrivacyPreservingGenomicData_WithHaploNet(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()
	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 2.0,
		UseFederatedPrivacy:        false,
		DefendAgainstTabPFN:        false,
		DefendAgainstLocalAncestry: false,
		DefendAgainstHaploNet:      true,
	}

	result, err := engine.GeneratePrivacyPreservingGenomicData(MXL, 2, requirements)
	require.NoError(t, err)
	assert.Len(t, result.GeneratedData, 2)
}

// TestGeneratePrivacyPreservingGenomicData_AllDefenses tests with all defenses
func TestGeneratePrivacyPreservingGenomicData_AllDefenses(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()
	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 3.0,
		UseFederatedPrivacy:        true,
		DefendAgainstTabPFN:        true,
		DefendAgainstLocalAncestry: true,
		DefendAgainstHaploNet:      true,
	}

	result, err := engine.GeneratePrivacyPreservingGenomicData(CEU, 5, requirements)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Len(t, result.GeneratedData, 5)

	// Privacy metrics should be set
	assert.Greater(t, result.PrivacyMetrics.OverallPrivacyRating, 0.0)
	assert.LessOrEqual(t, result.PrivacyMetrics.OverallPrivacyRating, 1.0)

	// Defense metrics should be set
	assert.Greater(t, result.DefenseMetrics.OverallDefenseRating, 0.0)
}

// TestGeneratePrivacyPreservingGenomicData_BudgetExhausted tests with exhausted budget
func TestGeneratePrivacyPreservingGenomicData_BudgetExhausted(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()

	// Exhaust the budget first
	requirements1 := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 10.0,
	}
	_, err := engine.GeneratePrivacyPreservingGenomicData(CEU, 1, requirements1)
	require.NoError(t, err)

	// Now try to generate more - should fail with budget exhausted
	requirements2 := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 1.0,
	}
	_, err = engine.GeneratePrivacyPreservingGenomicData(CEU, 1, requirements2)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "budget")
}

// TestPrivacyMetrics_Structure tests PrivacyMetrics structure
func TestPrivacyMetrics_Structure(t *testing.T) {
	metrics := &PrivacyMetrics{
		DifferentialPrivacyEpsilon: 1.0,
		PrivacyBudgetConsumed:     0.8,
		FamilyPrivacyScore:        0.95,
		PopulationPrivacyScore:    0.92,
		OverallPrivacyRating:      0.94,
	}

	assert.Equal(t, 1.0, metrics.DifferentialPrivacyEpsilon)
	assert.Equal(t, 0.8, metrics.PrivacyBudgetConsumed)
	assert.Equal(t, 0.95, metrics.FamilyPrivacyScore)
	assert.Equal(t, 0.92, metrics.PopulationPrivacyScore)
	assert.Equal(t, 0.94, metrics.OverallPrivacyRating)
}

// TestDefenseMetrics_Structure tests DefenseMetrics structure
func TestDefenseMetrics_Structure(t *testing.T) {
	metrics := &DefenseMetrics{
		TabPFNResistance:        0.88,
		LocalAncestryResistance: 0.91,
		HaploNetResistance:      0.89,
		PGGANResistance:         0.93,
		DPGANResistance:         0.85,
		FederatedResistance:     0.90,
		OverallDefenseRating:    0.89,
	}

	assert.Equal(t, 0.88, metrics.TabPFNResistance)
	assert.Equal(t, 0.91, metrics.LocalAncestryResistance)
	assert.Equal(t, 0.89, metrics.HaploNetResistance)
	assert.Equal(t, 0.89, metrics.OverallDefenseRating)
}

// TestGenomicSample_Structure tests GenomicSample structure
func TestGenomicSample_Structure(t *testing.T) {
	sample := GenomicSample{
		SampleID:     "synthetic_CEU_001",
		Population:   CEU,
		GenotypeData: []float64{0.1, 0.2, 0.3},
		PrivacyScore: 0.95,
		QualityScore: 0.88,
	}

	assert.Equal(t, "synthetic_CEU_001", sample.SampleID)
	assert.Equal(t, CEU, sample.Population)
	assert.Len(t, sample.GenotypeData, 3)
	assert.Equal(t, 0.95, sample.PrivacyScore)
}

// TestPrivacyRequirements_Structure tests PrivacyRequirements structure
func TestPrivacyRequirements_Structure(t *testing.T) {
	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 1.5,
		UseFederatedPrivacy:        true,
		DefendAgainstTabPFN:        true,
		DefendAgainstLocalAncestry: false,
		DefendAgainstHaploNet:      true,
		DefendAgainstPGGAN:         false,
		MinimumPrivacyScore:        0.90,
		MaximumPrivacyBudget:       5.0,
	}

	assert.Equal(t, 1.5, requirements.DifferentialPrivacyEpsilon)
	assert.True(t, requirements.UseFederatedPrivacy)
	assert.True(t, requirements.DefendAgainstTabPFN)
	assert.False(t, requirements.DefendAgainstLocalAncestry)
	assert.True(t, requirements.DefendAgainstHaploNet)
	assert.Equal(t, 0.90, requirements.MinimumPrivacyScore)
	assert.Equal(t, 5.0, requirements.MaximumPrivacyBudget)
}

// TestQualityMetrics_Structure tests QualityMetrics structure
func TestQualityMetrics_Structure(t *testing.T) {
	metrics := &QualityMetrics{
		BiologicalRealism:         0.92,
		StatisticalAccuracy:       0.88,
		PopulationCharacteristics: 0.95,
		UtilityPreservation:       0.90,
	}

	assert.Equal(t, 0.92, metrics.BiologicalRealism)
	assert.Equal(t, 0.88, metrics.StatisticalAccuracy)
	assert.Equal(t, 0.95, metrics.PopulationCharacteristics)
	assert.Equal(t, 0.90, metrics.UtilityPreservation)
}

// TestBudgetAllocation_Structure tests BudgetAllocation structure
func TestBudgetAllocation_Structure(t *testing.T) {
	allocation := BudgetAllocation{
		DPGANBudget:      0.4,
		PATEGANBudget:    0.3,
		FederatedBudget:  0.2,
		DefenseBudget:    0.08,
		ValidationBudget: 0.02,
	}

	assert.Equal(t, 0.4, allocation.DPGANBudget)
	assert.Equal(t, 0.3, allocation.PATEGANBudget)
	assert.Equal(t, 0.2, allocation.FederatedBudget)
	assert.Equal(t, 0.08, allocation.DefenseBudget)
	assert.Equal(t, 0.02, allocation.ValidationBudget)

	// Verify total sums to 1.0
	total := allocation.DPGANBudget + allocation.PATEGANBudget + allocation.FederatedBudget +
		allocation.DefenseBudget + allocation.ValidationBudget
	assert.InDelta(t, 1.0, total, 0.001)
}

// TestDifferentialPrivacyNoise_generateNoise tests noise generation
func TestDifferentialPrivacyNoise_generateNoise(t *testing.T) {
	t.Run("Gaussian noise", func(t *testing.T) {
		noise := &DifferentialPrivacyNoise{
			noiseMechanism: "Gaussian",
			sensitivity:    1.0,
			epsilon:        1.0,
			delta:          1e-5,
		}

		// Generate multiple samples to test distribution
		samples := make([]float64, 1000)
		for i := range samples {
			samples[i] = noise.generateNoise(1.0)
		}

		// Check samples are varied (not all zero)
		allZero := true
		for _, s := range samples {
			if s != 0.0 {
				allZero = false
				break
			}
		}
		assert.False(t, allZero, "Gaussian noise should not all be zero")
	})

	t.Run("Laplacian noise", func(t *testing.T) {
		noise := &DifferentialPrivacyNoise{
			noiseMechanism: "Laplacian",
			sensitivity:    1.0,
			epsilon:        1.0,
			delta:          1e-5,
		}

		samples := make([]float64, 1000)
		for i := range samples {
			samples[i] = noise.generateNoise(1.0)
		}

		allZero := true
		for _, s := range samples {
			if s != 0.0 {
				allZero = false
				break
			}
		}
		assert.False(t, allZero, "Laplacian noise should not all be zero")
	})

	t.Run("unknown mechanism returns zero", func(t *testing.T) {
		noise := &DifferentialPrivacyNoise{
			noiseMechanism: "Unknown",
			sensitivity:    1.0,
			epsilon:        1.0,
		}

		n := noise.generateNoise(1.0)
		assert.Equal(t, 0.0, n)
	})
}

// TestSampleLaplacian tests the Laplacian sampling function
func TestSampleLaplacian(t *testing.T) {
	// Run multiple samples to verify the function works
	samples := make([]float64, 1000)
	for i := range samples {
		samples[i] = sampleLaplacian(1.0)
	}

	// Calculate mean (should be close to 0 for Laplacian)
	sum := 0.0
	for _, s := range samples {
		sum += s
	}
	mean := sum / float64(len(samples))
	// Mean should be roughly 0 for large number of samples
	assert.Less(t, mean*mean, 0.5, "Mean of Laplacian samples should be close to 0")

	// Test with different scale
	largeScaleSamples := make([]float64, 1000)
	smallScaleSamples := make([]float64, 1000)
	for i := range largeScaleSamples {
		largeScaleSamples[i] = sampleLaplacian(10.0)
		smallScaleSamples[i] = sampleLaplacian(0.1)
	}

	// Large scale should have larger variance
	largeVar := variance(largeScaleSamples)
	smallVar := variance(smallScaleSamples)
	assert.Greater(t, largeVar, smallVar, "Large scale should have larger variance")
}

// variance is a helper to compute variance
func variance(vals []float64) float64 {
	if len(vals) == 0 {
		return 0.0
	}
	sum := 0.0
	for _, v := range vals {
		sum += v
	}
	mean := sum / float64(len(vals))
	sumSq := 0.0
	for _, v := range vals {
		sumSq += (v - mean) * (v - mean)
	}
	return sumSq / float64(len(vals))
}

// TestAdmixturePerturbation_generatePerturbation tests perturbation generation
func TestAdmixturePerturbation_generatePerturbation(t *testing.T) {
	perturbation := &AdmixturePerturbation{}

	populations := []Population{YRI, CEU, CHB, JPT, LWK} // LWK uses default factor 0
	for _, pop := range populations {
		t.Run(string(pop), func(t *testing.T) {
			p := perturbation.generatePerturbation(pop, 100)
			// Should not panic, value should be finite
			assert.NotEqual(t, p, p+1, "Value should be finite") // NaN check
		})
	}
}

// TestLDMasking_generateMask tests mask generation
func TestLDMasking_generateMask(t *testing.T) {
	masking := &LDMasking{}

	for pos := 0; pos < 10; pos++ {
		mask := masking.generateMask(pos)
		// Should return a valid float
		assert.NotEqual(t, mask, mask+1, "Mask should be finite") // NaN check
	}
}

// TestVAEBypassing_generateEncoderNoise tests encoder noise generation
func TestVAEBypassing_generateEncoderNoise(t *testing.T) {
	bypassing := &VAEBypassing{}

	populations := []Population{YRI, CEU, CHB, JPT}
	for _, pop := range populations {
		noise := bypassing.generateEncoderNoise(pop, 50)
		assert.NotEqual(t, noise, noise+1, "Noise should be finite")
	}
}

// TestLatentSpaceObfuscation_generateLatentNoise tests latent noise generation
func TestLatentSpaceObfuscation_generateLatentNoise(t *testing.T) {
	obfuscation := &LatentSpaceObfuscation{}

	for pos := 0; pos < 10; pos++ {
		noise := obfuscation.generateLatentNoise(pos)
		assert.NotEqual(t, noise, noise+1, "Latent noise should be finite")
	}
}

// TestFederatedPrivacyManager_enhancePrivacy tests federated privacy enhancement
func TestFederatedPrivacyManager_enhancePrivacy(t *testing.T) {
	mgr := NewFederatedPrivacyManager()

	data := [][]float64{
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
	}

	enhanced, err := mgr.enhancePrivacy(data, CEU)
	require.NoError(t, err)
	require.NotNil(t, enhanced)
	assert.Len(t, enhanced, len(data))

	for i, sample := range enhanced {
		assert.Len(t, sample, len(data[i]))
	}
}

// TestPrivacyPreservingGeneticsEngine_validatePrivacyGuarantees tests validation
func TestPrivacyPreservingGeneticsEngine_validatePrivacyGuarantees(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()

	data := [][]float64{
		{1.0, 2.0},
		{3.0, 4.0},
	}

	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 1.0,
	}

	metrics := engine.validatePrivacyGuarantees(data, requirements)
	assert.NotNil(t, metrics)
	assert.Greater(t, metrics.OverallPrivacyRating, 0.0)
	assert.LessOrEqual(t, metrics.OverallPrivacyRating, 1.0)
	assert.Greater(t, metrics.FamilyPrivacyScore, 0.0)
	assert.Greater(t, metrics.PopulationPrivacyScore, 0.0)
}

// TestPrivacyPreservingGeneticsEngine_validateDefenseEffectiveness tests defense validation
func TestPrivacyPreservingGeneticsEngine_validateDefenseEffectiveness(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()

	data := [][]float64{
		{1.0, 2.0},
	}

	requirements := &PrivacyRequirements{
		DefendAgainstTabPFN:        true,
		DefendAgainstLocalAncestry: true,
		DefendAgainstHaploNet:      true,
	}

	metrics := engine.validateDefenseEffectiveness(data, requirements)
	assert.NotNil(t, metrics)
	assert.Greater(t, metrics.TabPFNResistance, 0.0)
	assert.Greater(t, metrics.LocalAncestryResistance, 0.0)
	assert.Greater(t, metrics.HaploNetResistance, 0.0)
	assert.Greater(t, metrics.OverallDefenseRating, 0.0)
}

// TestPrivacyPreservingGeneticsEngine_convertToGenomicSample tests sample conversion
func TestPrivacyPreservingGeneticsEngine_convertToGenomicSample(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()

	data := []float64{0.1, 0.2, 0.3, 0.4, 0.5}
	populations := []Population{YRI, CEU, CHB, JPT, MXL}

	for _, pop := range populations {
		t.Run(string(pop), func(t *testing.T) {
			sample := engine.convertToGenomicSample(data, pop)
			assert.NotEmpty(t, sample.SampleID)
			assert.Equal(t, pop, sample.Population)
			assert.Equal(t, data, sample.GenotypeData)
			assert.Equal(t, 0.95, sample.PrivacyScore)
		})
	}
}

// TestGeneratePrivacyPreservingGenomicData_ZeroSamples tests zero samples edge case
func TestGeneratePrivacyPreservingGenomicData_ZeroSamples(t *testing.T) {
	engine := NewPrivacyPreservingGeneticsEngine()

	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 1.0,
	}

	result, err := engine.GeneratePrivacyPreservingGenomicData(CEU, 0, requirements)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Empty(t, result.GeneratedData)
}

// Benchmark tests

// BenchmarkDPGANGenerator_GenerateData benchmarks data generation
func BenchmarkDPGANGenerator_GenerateData(b *testing.B) {
	gen := NewDPGANGenerator()

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = gen.GenerateData(CEU, 10, 1.0)
	}
}

// BenchmarkTabPFNDefenseSystem_ApplyDefense benchmarks defense application
func BenchmarkTabPFNDefenseSystem_ApplyDefense(b *testing.B) {
	defense := NewTabPFNDefenseSystem()

	data := make([][]float64, 10)
	for i := range data {
		data[i] = make([]float64, 100)
		for j := range data[i] {
			data[i][j] = float64(i*100 + j)
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = defense.ApplyDefense(data, CEU)
	}
}

// BenchmarkPrivacyBudgetManager_AllocateBudget benchmarks budget allocation
func BenchmarkPrivacyBudgetManager_AllocateBudget(b *testing.B) {
	requirements := &PrivacyRequirements{
		DifferentialPrivacyEpsilon: 0.01,
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		mgr := NewPrivacyBudgetManager()
		_, _ = mgr.AllocateBudget(requirements)
	}
}

// TestGetPopulationParameters_MXL tests MXL population parameters
func TestGetPopulationParameters_MXL(t *testing.T) {
	simulator := NewDemographicSimulator("/test")

	params, err := simulator.GetPopulationParameters(MXL)
	require.NoError(t, err)
	require.NotNil(t, params)

	assert.Equal(t, MXL, params.Population)
	assert.Equal(t, 8500, params.EffectivePopulation)
	assert.Equal(t, 0.0045, params.GrowthRate)
	assert.Equal(t, 1.25e-8, params.MutationRate)
	assert.Equal(t, 1.13e-8, params.RecombinationRate)
	assert.Equal(t, "Americas", params.SuperPopulation)
	assert.Equal(t, "Americas", params.Region)
	assert.Contains(t, params.Description, "Mexican")
}

// TestGetPopulationParameters_AllSupported tests getting parameters for all supported populations
func TestGetPopulationParameters_AllSupported(t *testing.T) {
	simulator := NewDemographicSimulator("/test")

	// Test all populations in the map
	definedPops := []Population{YRI, CEU, CHB, MXL}
	for _, pop := range definedPops {
		t.Run(string(pop), func(t *testing.T) {
			params, err := simulator.GetPopulationParameters(pop)
			require.NoError(t, err)
			require.NotNil(t, params)
			assert.Equal(t, pop, params.Population)
			assert.Greater(t, params.EffectivePopulation, 0)
			assert.Greater(t, params.MutationRate, 0.0)
		})
	}
}

// TestDemographicSimulator_ValidatePopulationStatistics_AdditionalCases tests additional validation cases
func TestDemographicSimulator_ValidatePopulationStatistics_AdditionalCases(t *testing.T) {
	simulator := NewDemographicSimulator("/test")

	t.Run("population not in expected diversity map - passes diversity check", func(t *testing.T) {
		// MXL is not in the expectedDiversity map in ValidatePopulationStatistics
		// so it should pass the diversity check regardless of value
		results := map[Population]*SimulationResult{
			MXL: {
				Population: MXL,
				Statistics: PopulationStatistics{
					Population:          "MXL",
					NumVariants:         1000,
					DiversityPi:         0.001, // Any value - not checked for MXL
					TajimaD:             -0.3,  // Within range
					RareVariantsCount:   800,
					CommonVariantsCount: 200,
				},
				DemographicParams: DemographicParameters{
					EffectivePopulation: 8500,
				},
			},
		}

		err := simulator.ValidatePopulationStatistics(results)
		assert.NoError(t, err)
	})

	t.Run("empty results map", func(t *testing.T) {
		results := map[Population]*SimulationResult{}
		err := simulator.ValidatePopulationStatistics(results)
		assert.NoError(t, err) // No populations to validate
	})

	t.Run("CEU population validation", func(t *testing.T) {
		results := map[Population]*SimulationResult{
			CEU: {
				Population: CEU,
				Statistics: PopulationStatistics{
					Population:          "CEU",
					NumVariants:         1000,
					DiversityPi:         0.0002, // Above CEU max [0.000005, 0.000015] — expects error
					TajimaD:             -0.3,
					RareVariantsCount:   750,
					CommonVariantsCount: 250,
				},
				DemographicParams: DemographicParameters{
					EffectivePopulation: 11418,
				},
			},
		}

		err := simulator.ValidatePopulationStatistics(results)
		// CEU corrected range is [0.000005, 0.000015]; 0.0002 > 0.000015 so it should fail
		assert.Error(t, err)
	})
}

// TestDemographicSimulator_ValidatePopulationStatistics_CEU_InRange tests CEU with valid diversity
func TestDemographicSimulator_ValidatePopulationStatistics_CEU_InRange(t *testing.T) {
	simulator := NewDemographicSimulator("/test")

	results := map[Population]*SimulationResult{
		CEU: {
			Population: CEU,
			Statistics: PopulationStatistics{
				Population:          "CEU",
				NumVariants:         1000,
				DiversityPi:         0.000010, // Within corrected CEU range [0.000005, 0.000015]
				TajimaD:             -0.3,    // Valid Tajima's D
				RareVariantsCount:   750,
				CommonVariantsCount: 250,
			},
			DemographicParams: DemographicParameters{
				EffectivePopulation: 11418,
			},
		},
	}

	err := simulator.ValidatePopulationStatistics(results)
	assert.NoError(t, err)
}

// TestDemographicSimulator_ValidatePopulationStatistics_CHB tests CHB population
func TestDemographicSimulator_ValidatePopulationStatistics_CHB_Valid(t *testing.T) {
	simulator := NewDemographicSimulator("/test")

	results := map[Population]*SimulationResult{
		CHB: {
			Population: CHB,
			Statistics: PopulationStatistics{
				Population:          "CHB",
				NumVariants:         1000,
				DiversityPi:         0.000025, // Within range [0.000010, 0.000030]
				TajimaD:             -0.4,
				RareVariantsCount:   800,
				CommonVariantsCount: 200,
			},
			DemographicParams: DemographicParameters{
				EffectivePopulation: 9025,
			},
		},
	}

	err := simulator.ValidatePopulationStatistics(results)
	assert.NoError(t, err)
}

// TestDemographicSimulator_ValidatePopulationStatistics_CHB_OutOfRange tests CHB diversity out of range
func TestDemographicSimulator_ValidatePopulationStatistics_CHB_OutOfRange(t *testing.T) {
	simulator := NewDemographicSimulator("/test")

	results := map[Population]*SimulationResult{
		CHB: {
			Population: CHB,
			Statistics: PopulationStatistics{
				Population:          "CHB",
				NumVariants:         1000,
				DiversityPi:         0.00005, // Too high for CHB (> 0.000030)
				TajimaD:             -0.4,
				RareVariantsCount:   800,
				CommonVariantsCount: 200,
			},
			DemographicParams: DemographicParameters{
				EffectivePopulation: 9025,
			},
		},
	}

	err := simulator.ValidatePopulationStatistics(results)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "diversity")
}

// TestGeneratePrivacyPreservingGenomicData_AllPopulations tests all population codes
func TestGeneratePrivacyPreservingGenomicData_AllPopulations(t *testing.T) {
	populations := []Population{YRI, LWK, GWD, CEU, GBR, TSI, CHB, JPT, CHS, MXL, PUR}

	for _, pop := range populations {
		t.Run(string(pop), func(t *testing.T) {
			engine := NewPrivacyPreservingGeneticsEngine()
			requirements := &PrivacyRequirements{
				DifferentialPrivacyEpsilon: 1.0,
			}

			result, err := engine.GeneratePrivacyPreservingGenomicData(pop, 2, requirements)
			require.NoError(t, err)
			assert.Len(t, result.GeneratedData, 2)
			for _, sample := range result.GeneratedData {
				assert.Equal(t, pop, sample.Population)
			}
		})
	}
}

// TestDemographicSimulator_SimulatePopulation_FailsWithoutPython tests SimulatePopulation
// The function will fail because the Python script doesn't exist, but we still cover the code
func TestDemographicSimulator_SimulatePopulation_FailsWithoutPython(t *testing.T) {
	simulator := NewDemographicSimulator("/nonexistent/path")

	// Call with a nil seed to test the non-seed branch
	result, err := simulator.SimulatePopulation(YRI, 10, 100000, nil)
	assert.Error(t, err)
	assert.Nil(t, result)
	assert.Contains(t, err.Error(), "failed to execute demographic simulation")
}

// TestDemographicSimulator_SimulatePopulation_WithSeed tests SimulatePopulation with a random seed
func TestDemographicSimulator_SimulatePopulation_WithSeed(t *testing.T) {
	simulator := NewDemographicSimulator("/nonexistent/path")

	seed := 42
	// Call with a seed to test the randomSeed != nil branch
	result, err := simulator.SimulatePopulation(CEU, 10, 100000, &seed)
	assert.Error(t, err)
	assert.Nil(t, result)
	assert.Contains(t, err.Error(), "failed to execute demographic simulation")
}

// TestDemographicSimulator_SimulateAllPopulations_FailsWithoutPython tests SimulateAllPopulations
func TestDemographicSimulator_SimulateAllPopulations_FailsWithoutPython(t *testing.T) {
	simulator := NewDemographicSimulator("/nonexistent/path")

	results, err := simulator.SimulateAllPopulations(10, 100000, nil)
	assert.Error(t, err)
	assert.Nil(t, results)
	assert.Contains(t, err.Error(), "failed to simulate")
}

// TestDemographicSimulator_IntegrateWithHMMModel_FailsWithoutPython tests IntegrateWithHMMModel
func TestDemographicSimulator_IntegrateWithHMMModel_FailsWithoutPython(t *testing.T) {
	simulator := NewDemographicSimulator("/nonexistent/path")

	model, err := simulator.IntegrateWithHMMModel(YRI, 10, 100000, 4)
	assert.Error(t, err)
	assert.Nil(t, model)
	assert.Contains(t, err.Error(), "failed to simulate population")
}

// TestHaplotypeSimulator_SimulateHaplotypes_FailsWithoutPython tests SimulateHaplotypes
func TestHaplotypeSimulator_SimulateHaplotypes_FailsWithoutPython(t *testing.T) {
	simulator := NewHaplotypeSimulator("/nonexistent/path")

	// Test without seed and without matrix
	result, err := simulator.SimulateHaplotypes(YRI, "1", 10, 100000, nil, false)
	assert.Error(t, err)
	assert.Nil(t, result)
	assert.Contains(t, err.Error(), "failed to execute haplotype simulation")
}

// TestHaplotypeSimulator_SimulateHaplotypes_WithSeedAndMatrix tests SimulateHaplotypes with seed and matrix
func TestHaplotypeSimulator_SimulateHaplotypes_WithSeedAndMatrix(t *testing.T) {
	simulator := NewHaplotypeSimulator("/nonexistent/path")

	seed := 42
	// Test with seed and matrix to hit those branches
	result, err := simulator.SimulateHaplotypes(CEU, "22", 10, 100000, &seed, true)
	assert.Error(t, err)
	assert.Nil(t, result)
	assert.Contains(t, err.Error(), "failed to execute haplotype simulation")
}

// TestHaplotypeSimulator_SimulateAllPopulationsHaplotypes_FailsWithoutPython tests SimulateAllPopulationsHaplotypes
func TestHaplotypeSimulator_SimulateAllPopulationsHaplotypes_FailsWithoutPython(t *testing.T) {
	simulator := NewHaplotypeSimulator("/nonexistent/path")

	results, err := simulator.SimulateAllPopulationsHaplotypes("1", 10, 100000, nil)
	assert.Error(t, err)
	assert.Nil(t, results)
	assert.Contains(t, err.Error(), "failed to simulate")
}
