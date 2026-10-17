package authenticity

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewDemographicGenome tests the DemographicGenome constructor.
func TestNewDemographicGenome(t *testing.T) {
	populations := []string{"CEU", "YRI", "CHB", "JPT"}
	for _, pop := range populations {
		t.Run(pop, func(t *testing.T) {
			genome := NewDemographicGenome(pop)
			require.NotNil(t, genome)
			assert.Empty(t, genome.GetChromosomeNames())
		})
	}
}

// TestDemographicGenome_AddChromosome tests adding chromosomes with variants.
func TestDemographicGenome_AddChromosome(t *testing.T) {
	genome := NewDemographicGenome("CEU")

	variants := []DemographicVariant{
		{Position: 10, ReferenceAllele: 'A', AlternateAllele: 'T', Frequency: 0.5},
		{Position: 50, ReferenceAllele: 'G', AlternateAllele: 'C', Frequency: 0.3},
	}
	genome.AddChromosome("chr1", 200, variants)

	names := genome.GetChromosomeNames()
	assert.Len(t, names, 1)
	assert.Contains(t, names, "chr1")

	chr, err := genome.LoadChromosome("chr1")
	require.NoError(t, err)
	assert.NotNil(t, chr)
	assert.Equal(t, "chr1", chr.GetName())
	assert.Len(t, chr.GetSequence(), 200)
	assert.Equal(t, 200, chr.GetLength())
}

// TestDemographicGenome_LoadChromosome_NotFound tests loading a nonexistent chromosome.
func TestDemographicGenome_LoadChromosome_NotFound(t *testing.T) {
	genome := NewDemographicGenome("CEU")
	_, err := genome.LoadChromosome("nonexistent")
	assert.Error(t, err)
}

// TestDemographicGenome_MultipleChromosomes tests adding and loading multiple chromosomes.
func TestDemographicGenome_MultipleChromosomes(t *testing.T) {
	genome := NewDemographicGenome("EUR")

	for i := 1; i <= 5; i++ {
		name := "chr" + string(rune('0'+i))
		genome.AddChromosome(name, 1000, nil)
	}

	names := genome.GetChromosomeNames()
	assert.Len(t, names, 5)
}

// TestDemographicGenome_Sequence tests that sequence bases are valid DNA.
func TestDemographicGenome_Sequence(t *testing.T) {
	genome := NewDemographicGenome("CEU")
	genome.AddChromosome("chr1", 500, nil)

	chr, err := genome.LoadChromosome("chr1")
	require.NoError(t, err)

	seq := chr.GetSequence()
	assert.Len(t, seq, 500)
	for _, base := range seq {
		assert.Contains(t, "ATGC", string(base),
			"invalid base in sequence: %c", base)
	}
}

// TestDemographicGenome_VariantApplication tests that variants alter the sequence.
func TestDemographicGenome_VariantApplication(t *testing.T) {
	// Use high-frequency variant to ensure it gets applied most of the time.
	// We run multiple iterations to verify that at least sometimes it's applied.
	applied := false
	for attempt := 0; attempt < 50; attempt++ {
		genome := NewDemographicGenome("CEU")
		variants := []DemographicVariant{
			{Position: 0, ReferenceAllele: 'A', AlternateAllele: 'C', Frequency: 0.9},
		}
		genome.AddChromosome("chr1", 100, variants)

		chr, err := genome.LoadChromosome("chr1")
		require.NoError(t, err)

		if chr.GetSequence()[0] == 'C' {
			applied = true
			break
		}
	}
	assert.True(t, applied, "high-frequency variant should be applied at least once in 50 attempts")
}

// TestNewEnhancedHMMKmerModel tests the enhanced model constructor.
func TestNewEnhancedHMMKmerModel(t *testing.T) {
	model := NewEnhancedHMMKmerModel(4, "CEU")
	require.NotNil(t, model)

	assert.Equal(t, "CEU", model.GetPopulation())
	assert.False(t, model.IsEnhanced())
	assert.Equal(t, 4, model.K())
}

// TestEnhancedHMMKmerModel_IntegrateDemographicVariants tests demographic integration.
func TestEnhancedHMMKmerModel_IntegrateDemographicVariants(t *testing.T) {
	model := NewEnhancedHMMKmerModel(4, "CEU")
	require.NotNil(t, model)

	variants := []DemographicVariant{
		{Position: 10, ReferenceAllele: 'A', AlternateAllele: 'T', Frequency: 0.1},
		{Position: 50, ReferenceAllele: 'G', AlternateAllele: 'C', Frequency: 0.05},
	}
	stats := map[string]float64{
		"tajima_d":     0.2,
		"nucleotide_d": 0.001,
	}

	err := model.IntegrateDemographicVariants("chr1", 2000, variants, stats)
	require.NoError(t, err)

	assert.True(t, model.IsEnhanced())
	assert.True(t, model.IsTrained())
}

// TestEnhancedHMMKmerModel_GetEnhancedModelStats tests enhanced model statistics.
func TestEnhancedHMMKmerModel_GetEnhancedModelStats(t *testing.T) {
	model := NewEnhancedHMMKmerModel(4, "CEU")
	require.NotNil(t, model)

	err := model.IntegrateDemographicVariants("chr1", 2000, nil, nil)
	require.NoError(t, err)

	stats := model.GetEnhancedModelStats()
	require.NotNil(t, stats)

	assert.Equal(t, "CEU", stats["population"])
	assert.True(t, stats["demographic_enhanced"].(bool))
	assert.NotNil(t, stats["base_model"])
	assert.Equal(t, 1, stats["chromosomes_loaded"])
}

// TestEnhancedHMMKmerModel_GetPopulationStatistics tests population statistics storage.
func TestEnhancedHMMKmerModel_GetPopulationStatistics(t *testing.T) {
	model := NewEnhancedHMMKmerModel(4, "JPT")
	require.NotNil(t, model)

	stats := map[string]float64{
		"tajima_d":  -0.5,
		"fst":       0.12,
		"pi":        0.0009,
	}
	err := model.IntegrateDemographicVariants("chr1", 1000, nil, stats)
	require.NoError(t, err)

	popStats := model.GetPopulationStatistics()
	assert.InDelta(t, -0.5, popStats["tajima_d"], 1e-9)
	assert.InDelta(t, 0.12, popStats["fst"], 1e-9)
	assert.InDelta(t, 0.0009, popStats["pi"], 1e-9)
}

// TestEnhancedHMMKmerModel_GeneratePopulationSequence_NotEnhanced tests generation guard.
func TestEnhancedHMMKmerModel_GeneratePopulationSequence_NotEnhanced(t *testing.T) {
	model := NewEnhancedHMMKmerModel(4, "CEU")
	require.NotNil(t, model)

	_, err := model.GeneratePopulationSequence(100)
	assert.Error(t, err)
}

// TestEnhancedHMMKmerModel_GeneratePopulationSequence tests successful generation.
func TestEnhancedHMMKmerModel_GeneratePopulationSequence(t *testing.T) {
	model := NewEnhancedHMMKmerModel(4, "CEU")
	require.NotNil(t, model)

	err := model.IntegrateDemographicVariants("chr1", 5000, nil, nil)
	require.NoError(t, err)

	result, err := model.GeneratePopulationSequence(100)
	require.NoError(t, err)
	require.NotNil(t, result)

	assert.Len(t, result.Sequence, 100)
	for _, base := range result.Sequence {
		assert.Contains(t, "ATGC", string(base))
	}
}

// TestEnhancedHMMKmerModel_MultiplePopulations tests models for different populations.
func TestEnhancedHMMKmerModel_MultiplePopulations(t *testing.T) {
	populations := []string{"CEU", "YRI", "CHB", "MXL"}
	for _, pop := range populations {
		t.Run(pop, func(t *testing.T) {
			model := NewEnhancedHMMKmerModel(3, pop)
			require.NotNil(t, model)
			assert.Equal(t, pop, model.GetPopulation())

			err := model.IntegrateDemographicVariants("chr1", 2000, nil, nil)
			require.NoError(t, err)
			assert.True(t, model.IsEnhanced())
		})
	}
}

// TestGenerateReferenceSequence tests the internal sequence generator.
func TestGenerateReferenceSequence(t *testing.T) {
	lengths := []int{100, 500, 1000}
	for _, length := range lengths {
		seq := generateReferenceSequence(length)
		assert.Len(t, seq, length)
		for _, base := range seq {
			assert.Contains(t, "ATGC", string(base),
				"unexpected base in reference sequence: %c", base)
		}
	}
}

// TestCalculateGCContent tests the GC content helper.
func TestCalculateGCContent(t *testing.T) {
	tests := []struct {
		name     string
		sequence []byte
		expected float64
		delta    float64
	}{
		{"all_A", []byte("AAAA"), 0.0, 1e-9},
		{"all_G", []byte("GGGG"), 1.0, 1e-9},
		{"half_GC", []byte("ATGC"), 0.5, 1e-9},
		{"all_C", []byte("CCCC"), 1.0, 1e-9},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := calculateGCContent(tt.sequence)
			assert.InDelta(t, tt.expected, result, tt.delta)
		})
	}
}

// TestAbs tests the internal absolute value helper.
func TestAbs(t *testing.T) {
	assert.InDelta(t, 1.5, abs(-1.5), 1e-9)
	assert.InDelta(t, 1.5, abs(1.5), 1e-9)
	assert.InDelta(t, 0.0, abs(0.0), 1e-9)
}

// TestEnhancedHMMKmerModel_EnhancedStats_MultipleChromosomes tests stats with multiple chromosomes.
func TestEnhancedHMMKmerModel_EnhancedStats_MultipleChromosomes(t *testing.T) {
	model := NewEnhancedHMMKmerModel(4, "GBR")
	require.NotNil(t, model)

	require.NoError(t, model.IntegrateDemographicVariants("chr1", 1000, nil, nil))
	require.NoError(t, model.IntegrateDemographicVariants("chr2", 1000, nil, nil))

	stats := model.GetEnhancedModelStats()
	require.NotNil(t, stats)
	// After two integrations chromosomes_loaded should be 2.
	assert.Equal(t, 2, stats["chromosomes_loaded"])
	assert.NotNil(t, stats["chromosome_names"])
	assert.NotNil(t, stats["total_variants_integrated"])
}

func TestEnhancedHMMKmerModel_AdjustQualityScores_NoPanic(t *testing.T) {
    // adjustQualityScores has an empty body - just verify it doesn't panic
    model := NewEnhancedHMMKmerModel(4, "GBR")
    require.NotNil(t, model)
    qualities := []byte{30, 35, 32, 28, 40}
    model.adjustQualityScores(qualities)
    // If we reach here, it didn't panic
}

func TestEnhancedHMMKmerModel_AdjustQualityScores_Empty(t *testing.T) {
    model := NewEnhancedHMMKmerModel(4, "GBR")
    model.adjustQualityScores([]byte{})
    // Empty input should not panic
}
