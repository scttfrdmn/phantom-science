package population

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewHaplotypeSimulator tests haplotype simulator creation
func TestNewHaplotypeSimulator(t *testing.T) {
	projectRoot := "/test/project"
	simulator := NewHaplotypeSimulator(projectRoot)

	assert.NotNil(t, simulator)
	assert.Contains(t, simulator.pythonPath, "phantom_science_env/bin/python3")
	assert.Contains(t, simulator.scriptPath, "internal/population/haplotype_cli.py")
	assert.Contains(t, simulator.envPath, "phantom_science_env")
}

// TestRecombinationMap_Validation tests recombination map structure
func TestRecombinationMap_Validation(t *testing.T) {
	recombMap := RecombinationMap{
		Chromosome:     "1",
		Positions:      []int{1000, 2000, 3000, 4000, 5000},
		Rates:          []float64{1.2e-8, 1.5e-8, 1.1e-8, 1.3e-8, 1.4e-8},
		MapLength:      2.5, // centiMorgans
		PhysicalLength: 5000,
	}

	// Basic validation
	assert.Equal(t, "1", recombMap.Chromosome)
	assert.Len(t, recombMap.Positions, 5)
	assert.Len(t, recombMap.Rates, 5)
	assert.Equal(t, len(recombMap.Positions), len(recombMap.Rates))
	assert.Greater(t, recombMap.MapLength, 0.0)
	assert.Greater(t, recombMap.PhysicalLength, 0)

	// Validate positions are sorted
	for i := 1; i < len(recombMap.Positions); i++ {
		assert.Greater(t, recombMap.Positions[i], recombMap.Positions[i-1], 
			"Positions should be in ascending order")
	}

	// Validate recombination rates are realistic
	for _, rate := range recombMap.Rates {
		assert.Greater(t, rate, 0.0, "Recombination rates should be positive")
		assert.Less(t, rate, 1e-6, "Recombination rates should be realistic")
	}

	// Validate map length vs physical length ratio
	expectedRatio := recombMap.MapLength / (float64(recombMap.PhysicalLength) / 1e6) // cM/Mb
	assert.Greater(t, expectedRatio, 0.5, "Map length ratio should be reasonable")
	assert.Less(t, expectedRatio, 1000.0, "Map length ratio should not be excessive")
}

// TestHaplotypeBlock_Validation tests haplotype block structure
func TestHaplotypeBlock_Validation(t *testing.T) {
	block := HaplotypeBlock{
		StartPosition: 1000000,
		EndPosition:   1050000,
		Length:        50000,
		VariantCount:  25,
		RSquaredMean:  0.75,
		Population:    YRI,
	}

	// Basic validation
	assert.Greater(t, block.StartPosition, 0)
	assert.Greater(t, block.EndPosition, block.StartPosition)
	assert.Equal(t, block.Length, block.EndPosition-block.StartPosition)
	assert.Greater(t, block.VariantCount, 0)
	assert.GreaterOrEqual(t, block.RSquaredMean, 0.0)
	assert.LessOrEqual(t, block.RSquaredMean, 1.0)
	assert.Equal(t, YRI, block.Population)

	// Validate LD structure expectations
	assert.Greater(t, block.RSquaredMean, 0.5, "High LD blocks should have r² > 0.5")
	assert.Less(t, block.Length, 200000, "Typical LD blocks should be < 200kb")
	
	// Validate variant density
	variantDensity := float64(block.VariantCount) / float64(block.Length) * 1000
	assert.Greater(t, variantDensity, 0.1, "Should have reasonable variant density per kb")
	assert.Less(t, variantDensity, 10.0, "Variant density should not be excessive")
}

// TestHaplotypeVariant_Validation tests haplotype variant structure
func TestHaplotypeVariant_Validation(t *testing.T) {
	variant := HaplotypeVariant{
		Position:       1234567,
		Alleles:        []string{"A", "T"},
		AncestralState: "A",
		Frequency:      0.25,
	}

	// Basic validation
	assert.Greater(t, variant.Position, 0)
	assert.Len(t, variant.Alleles, 2)
	assert.Contains(t, variant.Alleles, variant.AncestralState)
	assert.GreaterOrEqual(t, variant.Frequency, 0.0)
	assert.LessOrEqual(t, variant.Frequency, 1.0)

	// Validate alleles are valid DNA bases
	validBases := map[string]bool{"A": true, "T": true, "G": true, "C": true}
	for _, allele := range variant.Alleles {
		assert.True(t, validBases[allele], "Allele %s should be a valid DNA base", allele)
	}
}

// TestHaplotypeStatistics_Validation tests haplotype statistics structure
func TestHaplotypeStatistics_Validation(t *testing.T) {
	stats := HaplotypeStatistics{
		DiversityPi:         0.001234,
		SegregatingSites:    1500,
		NumVariants:         1600,
		NumHaplotypes:       200,
		MeanRSquared:        0.65,
		MedianRSquared:      0.70,
		HighLDFraction:      0.45,
		EffectiveRecombRate: 1.2e-8,
		GeneticMapLengthCM:  2.5,
		HaplotypeDiversity:  0.95,
		LDPairsAnalyzed:     50000,
	}

	// Basic validation
	assert.Greater(t, stats.DiversityPi, 0.0)
	assert.Greater(t, stats.SegregatingSites, 0)
	assert.Greater(t, stats.NumVariants, 0)
	assert.Greater(t, stats.NumHaplotypes, 0)
	assert.GreaterOrEqual(t, stats.MeanRSquared, 0.0)
	assert.LessOrEqual(t, stats.MeanRSquared, 1.0)
	assert.GreaterOrEqual(t, stats.MedianRSquared, 0.0)
	assert.LessOrEqual(t, stats.MedianRSquared, 1.0)
	assert.GreaterOrEqual(t, stats.HighLDFraction, 0.0)
	assert.LessOrEqual(t, stats.HighLDFraction, 1.0)
	assert.Greater(t, stats.EffectiveRecombRate, 0.0)
	assert.Greater(t, stats.GeneticMapLengthCM, 0.0)
	assert.GreaterOrEqual(t, stats.HaplotypeDiversity, 0.0)
	assert.LessOrEqual(t, stats.HaplotypeDiversity, 1.0)
	assert.Greater(t, stats.LDPairsAnalyzed, 0)

	// Biological realism checks
	assert.LessOrEqual(t, stats.SegregatingSites, stats.NumVariants, 
		"Segregating sites should not exceed total variants")
	assert.Less(t, stats.EffectiveRecombRate, 1e-6, 
		"Recombination rate should be biologically realistic")
	assert.Greater(t, stats.HaplotypeDiversity, 0.8, 
		"Haplotype diversity should be high for most populations")
	assert.Less(t, stats.GeneticMapLengthCM, 10.0, 
		"Genetic map length should be reasonable for typical regions")
}

// TestHaplotypeSimulationResult_Validation tests simulation result structure
func TestHaplotypeSimulationResult_Validation(t *testing.T) {
	result := HaplotypeSimulationResult{
		Population:     CEU,
		Chromosome:     "22",
		NumVariants:    1200,
		NumTrees:       45,
		NumSamples:     100,
		SequenceLength: 1000000,
		SimulationTime: 3.5,
		MemoryUsedMB:   250.0,
		Variants: []HaplotypeVariant{
			{
				Position:       100000,
				Alleles:        []string{"G", "A"},
				AncestralState: "G",
				Frequency:      0.15,
			},
			{
				Position:       200000,
				Alleles:        []string{"C", "T"},
				AncestralState: "C",
				Frequency:      0.35,
			},
		},
		RecombinationMap: RecombinationMap{
			Chromosome:     "22",
			Positions:      []int{0, 250000, 500000, 750000, 1000000},
			Rates:          []float64{1.1e-8, 1.3e-8, 1.2e-8, 1.4e-8, 1.1e-8},
			MapLength:      1.8,
			PhysicalLength: 1000000,
		},
		LDBlocks: []HaplotypeBlock{
			{
				StartPosition: 50000,
				EndPosition:   150000,
				Length:        100000,
				VariantCount:  15,
				RSquaredMean:  0.8,
				Population:    CEU,
			},
		},
		Statistics: HaplotypeStatistics{
			DiversityPi:         0.0008,
			SegregatingSites:    1150,
			NumVariants:         1200,
			NumHaplotypes:       180,
			MeanRSquared:        0.72,
			MedianRSquared:      0.75,
			HighLDFraction:      0.55,
			EffectiveRecombRate: 1.2e-8,
			GeneticMapLengthCM:  1.8,
			HaplotypeDiversity:  0.92,
			LDPairsAnalyzed:     75000,
		},
		HaplotypeMatrix: [][]int{
			{0, 1}, // 2 variants to match len(result.Variants)
			{1, 0},
			{0, 1},
		},
	}

	// Basic structure validation
	assert.Equal(t, CEU, result.Population)
	assert.Equal(t, "22", result.Chromosome)
	assert.Greater(t, result.NumVariants, 0)
	assert.Greater(t, result.NumTrees, 0)
	assert.Greater(t, result.NumSamples, 0)
	assert.Greater(t, result.SequenceLength, 0)
	assert.Greater(t, result.SimulationTime, 0.0)
	assert.Greater(t, result.MemoryUsedMB, 0.0)

	// Validate consistency between components
	assert.Equal(t, len(result.Variants), 2) // We created 2 variants in the test data
	assert.Equal(t, result.Statistics.NumVariants, result.NumVariants)
	assert.Equal(t, result.RecombinationMap.Chromosome, result.Chromosome)
	assert.Equal(t, result.RecombinationMap.PhysicalLength, result.SequenceLength)

	// Validate LD blocks are within sequence
	for _, block := range result.LDBlocks {
		assert.GreaterOrEqual(t, block.StartPosition, 0)
		assert.LessOrEqual(t, block.EndPosition, result.SequenceLength)
		assert.Equal(t, result.Population, block.Population)
	}

	// Validate variants are within sequence and sorted
	for i, variant := range result.Variants {
		assert.GreaterOrEqual(t, variant.Position, 0)
		assert.LessOrEqual(t, variant.Position, result.SequenceLength)
		
		if i > 0 {
			assert.Greater(t, variant.Position, result.Variants[i-1].Position,
				"Variants should be sorted by position")
		}
	}

	// Validate haplotype matrix dimensions
	if len(result.HaplotypeMatrix) > 0 {
		expectedCols := len(result.Variants) // Use actual number of variants in test data
		for i, haplotype := range result.HaplotypeMatrix {
			assert.Len(t, haplotype, expectedCols, 
				"Haplotype %d should have %d variants", i, expectedCols)
			
			// Validate allele indices are valid
			for j, allele := range haplotype {
				assert.GreaterOrEqual(t, allele, 0, 
					"Allele index should be non-negative")
				if j < len(result.Variants) {
					assert.Less(t, allele, len(result.Variants[j].Alleles),
						"Allele index should be valid for variant %d", j)
				}
			}
		}
	}
}

// TestHaplotypeSimulator_ValidateLDPatterns tests LD pattern validation
func TestHaplotypeSimulator_ValidateLDPatterns(t *testing.T) {
	simulator := NewHaplotypeSimulator("/test")

	// Create mock simulation results for testing
	results := map[Population]*HaplotypeSimulationResult{
		YRI: {
			Population: YRI,
			Statistics: HaplotypeStatistics{
				DiversityPi:         0.0012, // Higher diversity (African)
				MeanRSquared:        0.05,
				EffectiveRecombRate: 1.2e-8,
				GeneticMapLengthCM:  2.0,
			},
			SequenceLength: 1000000,
			LDBlocks: []HaplotypeBlock{
				{Length: 50000, VariantCount: 10},
				{Length: 75000, VariantCount: 15},
			},
		},
		CEU: {
			Population: CEU,
			Statistics: HaplotypeStatistics{
				DiversityPi:         0.0008, // Lower diversity (European)
				MeanRSquared:        0.06,
				EffectiveRecombRate: 1.1e-8,
				GeneticMapLengthCM:  1.8,
			},
			SequenceLength: 1000000,
			LDBlocks: []HaplotypeBlock{
				{Length: 80000, VariantCount: 12},
			},
		},
		CHB: {
			Population: CHB,
			Statistics: HaplotypeStatistics{
				DiversityPi:         0.0007, // Lower diversity (East Asian)
				MeanRSquared:        0.07,
				EffectiveRecombRate: 1.0e-8,
				GeneticMapLengthCM:  1.5,
			},
			SequenceLength: 1000000,
			LDBlocks: []HaplotypeBlock{
				{Length: 90000, VariantCount: 8},
				{Length: 60000, VariantCount: 6},
			},
		},
	}

	validation, err := simulator.ValidateLDPatterns(results)
	require.NoError(t, err)
	require.NotNil(t, validation)

	// Check that validation tests are present
	expectedTests := []string{
		"african_highest_diversity",
		"reasonable_ld_structure", 
		"haplotype_blocks_present",
		"recombination_variation",
		"realistic_map_lengths",
	}

	for _, test := range expectedTests {
		_, exists := validation[test]
		assert.True(t, exists, "Validation test %s should be present", test)
	}

	// African populations should have highest diversity
	assert.True(t, validation["african_highest_diversity"], 
		"African populations should have highest diversity")

	// LD structure should be reasonable
	assert.True(t, validation["reasonable_ld_structure"],
		"LD structure should be reasonable")

	// Haplotype blocks should be present
	assert.True(t, validation["haplotype_blocks_present"],
		"Haplotype blocks should be present")

	// There should be variation in recombination rates
	assert.True(t, validation["recombination_variation"],
		"There should be variation in recombination rates")

	// Map lengths should be realistic
	assert.True(t, validation["realistic_map_lengths"],
		"Genetic map lengths should be realistic")
}

// TestHaplotypeSimulator_AnalyzeCrossPopulationLD tests cross-population LD analysis
func TestHaplotypeSimulator_AnalyzeCrossPopulationLD(t *testing.T) {
	simulator := NewHaplotypeSimulator("/test")

	results := map[Population]*HaplotypeSimulationResult{
		YRI: {
			Population: YRI,
			Statistics: HaplotypeStatistics{
				DiversityPi:  0.0012,
				MeanRSquared: 0.05,
			},
			LDBlocks: []HaplotypeBlock{{Length: 50000}, {Length: 60000}},
		},
		CEU: {
			Population: CEU,
			Statistics: HaplotypeStatistics{
				DiversityPi:  0.0008,
				MeanRSquared: 0.07,
			},
			LDBlocks: []HaplotypeBlock{{Length: 80000}},
		},
	}

	analysis := simulator.AnalyzeCrossPopulationLD(results)
	require.NotNil(t, analysis)

	// Check required fields are present
	assert.Contains(t, analysis, "diversity_by_population")
	assert.Contains(t, analysis, "ld_metrics_by_population")
	assert.Contains(t, analysis, "block_counts_by_population")
	assert.Contains(t, analysis, "diversity_mean")
	assert.Contains(t, analysis, "diversity_cv")
	assert.Contains(t, analysis, "ld_mean")
	assert.Contains(t, analysis, "ld_cv")

	// Validate diversity metrics
	diversities := analysis["diversity_by_population"].(map[string]float64)
	assert.Equal(t, 0.0012, diversities["YRI"])
	assert.Equal(t, 0.0008, diversities["CEU"])

	// Validate LD metrics
	ldMetrics := analysis["ld_metrics_by_population"].(map[string]float64)
	assert.Equal(t, 0.05, ldMetrics["YRI"])
	assert.Equal(t, 0.07, ldMetrics["CEU"])

	// Validate block counts
	blockCounts := analysis["block_counts_by_population"].(map[string]int)
	assert.Equal(t, 2, blockCounts["YRI"])
	assert.Equal(t, 1, blockCounts["CEU"])

	// Validate summary statistics
	diversityMean := analysis["diversity_mean"].(float64)
	assert.Equal(t, 0.001, diversityMean) // (0.0012 + 0.0008) / 2

	ldMean := analysis["ld_mean"].(float64)
	assert.InDelta(t, 0.06, ldMean, 1e-10) // (0.05 + 0.07) / 2
}

// TestHaplotypeSimulator_GetRecombinationMap tests recombination map extraction
func TestHaplotypeSimulator_GetRecombinationMap(t *testing.T) {
	simulator := NewHaplotypeSimulator("/test")

	expectedMap := RecombinationMap{
		Chromosome:     "1",
		Positions:      []int{0, 500000, 1000000},
		Rates:          []float64{1.1e-8, 1.3e-8, 1.2e-8},
		MapLength:      2.0,
		PhysicalLength: 1000000,
	}

	result := &HaplotypeSimulationResult{
		RecombinationMap: expectedMap,
	}

	retrievedMap := simulator.GetRecombinationMap(result)

	assert.Equal(t, expectedMap.Chromosome, retrievedMap.Chromosome)
	assert.Equal(t, expectedMap.Positions, retrievedMap.Positions)
	assert.Equal(t, expectedMap.Rates, retrievedMap.Rates)
	assert.Equal(t, expectedMap.MapLength, retrievedMap.MapLength)
	assert.Equal(t, expectedMap.PhysicalLength, retrievedMap.PhysicalLength)
}

// TestHaplotypeSimulator_GetLDStatistics tests LD statistics extraction
func TestHaplotypeSimulator_GetLDStatistics(t *testing.T) {
	simulator := NewHaplotypeSimulator("/test")

	expectedStats := HaplotypeStatistics{
		DiversityPi:         0.001,
		SegregatingSites:    1000,
		NumVariants:         1100,
		NumHaplotypes:       200,
		MeanRSquared:        0.65,
		MedianRSquared:      0.68,
		HighLDFraction:      0.45,
		EffectiveRecombRate: 1.2e-8,
		GeneticMapLengthCM:  2.5,
		HaplotypeDiversity:  0.95,
		LDPairsAnalyzed:     50000,
	}

	result := &HaplotypeSimulationResult{
		Statistics: expectedStats,
	}

	retrievedStats := simulator.GetLDStatistics(result)

	assert.Equal(t, expectedStats.DiversityPi, retrievedStats.DiversityPi)
	assert.Equal(t, expectedStats.SegregatingSites, retrievedStats.SegregatingSites)
	assert.Equal(t, expectedStats.NumVariants, retrievedStats.NumVariants)
	assert.Equal(t, expectedStats.NumHaplotypes, retrievedStats.NumHaplotypes)
	assert.Equal(t, expectedStats.MeanRSquared, retrievedStats.MeanRSquared)
	assert.Equal(t, expectedStats.MedianRSquared, retrievedStats.MedianRSquared)
	assert.Equal(t, expectedStats.HighLDFraction, retrievedStats.HighLDFraction)
	assert.Equal(t, expectedStats.EffectiveRecombRate, retrievedStats.EffectiveRecombRate)
	assert.Equal(t, expectedStats.GeneticMapLengthCM, retrievedStats.GeneticMapLengthCM)
	assert.Equal(t, expectedStats.HaplotypeDiversity, retrievedStats.HaplotypeDiversity)
	assert.Equal(t, expectedStats.LDPairsAnalyzed, retrievedStats.LDPairsAnalyzed)
}

// Test helper functions

// TestCalculateGroupMeanDiversity tests diversity calculation helper
func TestCalculateGroupMeanDiversity(t *testing.T) {
	results := map[Population]*HaplotypeSimulationResult{
		YRI: {Statistics: HaplotypeStatistics{DiversityPi: 0.0012}},
		LWK: {Statistics: HaplotypeStatistics{DiversityPi: 0.0011}},
		GWD: {Statistics: HaplotypeStatistics{DiversityPi: 0.0013}},
		CEU: {Statistics: HaplotypeStatistics{DiversityPi: 0.0008}},
	}

	africanPops := []Population{YRI, LWK, GWD}
	africanDiversity := calculateGroupMeanDiversity(results, africanPops)

	expectedMean := (0.0012 + 0.0011 + 0.0013) / 3.0
	assert.InDelta(t, expectedMean, africanDiversity, 1e-6)

	// Test with missing population
	missingPops := []Population{YRI, Population("MISSING")}
	diversityWithMissing := calculateGroupMeanDiversity(results, missingPops)
	assert.Equal(t, 0.0012, diversityWithMissing) // Should only count YRI

	// Test with no populations
	emptyDiversity := calculateGroupMeanDiversity(results, []Population{})
	assert.Equal(t, 0.0, emptyDiversity)
}

// TestCalculateMean tests mean calculation helper
func TestCalculateMean(t *testing.T) {
	values := []float64{1.0, 2.0, 3.0, 4.0, 5.0}
	mean := calculateMean(values)
	assert.Equal(t, 3.0, mean)

	// Test empty slice
	emptyMean := calculateMean([]float64{})
	assert.Equal(t, 0.0, emptyMean)

	// Test single value
	singleMean := calculateMean([]float64{42.0})
	assert.Equal(t, 42.0, singleMean)
}

// TestCalculateCV tests coefficient of variation calculation
func TestCalculateCV(t *testing.T) {
	// Test normal case
	values := []float64{1.0, 2.0, 3.0, 4.0, 5.0}
	cv := calculateCV(values)
	assert.Greater(t, cv, 0.0)
	assert.Less(t, cv, 1.0) // CV should be reasonable

	// Test with identical values (should be 0)
	identicalValues := []float64{5.0, 5.0, 5.0, 5.0}
	identicalCV := calculateCV(identicalValues)
	assert.Equal(t, 0.0, identicalCV)

	// Test empty slice
	emptyCV := calculateCV([]float64{})
	assert.Equal(t, 0.0, emptyCV)

	// Test single value
	singleCV := calculateCV([]float64{42.0})
	assert.Equal(t, 0.0, singleCV)

	// Test with zero mean
	zeroMeanValues := []float64{-1.0, 0.0, 1.0}
	zeroMeanCV := calculateCV(zeroMeanValues)
	assert.Equal(t, 0.0, zeroMeanCV)
}

// TestSqrt tests square root implementation
func TestSqrt(t *testing.T) {
	testCases := []struct {
		input    float64
		expected float64
		tolerance float64
	}{
		{0.0, 0.0, 0.0},
		{1.0, 1.0, 1e-10},
		{4.0, 2.0, 1e-10},
		{9.0, 3.0, 1e-10},
		{16.0, 4.0, 1e-10},
		{2.0, 1.414213562, 1e-6},
	}

	for _, tc := range testCases {
		t.Run(fmt.Sprintf("sqrt(%.1f)", tc.input), func(t *testing.T) {
			result := sqrt(tc.input)
			assert.InDelta(t, tc.expected, result, tc.tolerance)
		})
	}

	// Test negative input
	negativeResult := sqrt(-1.0)
	assert.Equal(t, 0.0, negativeResult)
}

// BenchmarkCalculateMean benchmarks mean calculation
func BenchmarkCalculateMean(b *testing.B) {
	values := make([]float64, 1000)
	for i := range values {
		values[i] = float64(i)
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		calculateMean(values)
	}
}

// BenchmarkCalculateCV benchmarks coefficient of variation calculation
func BenchmarkCalculateCV(b *testing.B) {
	values := make([]float64, 1000)
	for i := range values {
		values[i] = float64(i % 100)
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		calculateCV(values)
	}
}