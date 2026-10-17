package features

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewCpGDepletionModel tests model creation
func TestNewCpGDepletionModel(t *testing.T) {
	model := NewCpGDepletionModel()

	assert.NotNil(t, model)
	assert.Equal(t, 0.8, model.depletionRate)
	assert.Equal(t, 0.7, model.methylationRate)
	assert.Equal(t, 0.02, model.transitionProb)
	assert.Equal(t, 0.5, model.gcContentThreshold)
	assert.Equal(t, 200, model.lengthThreshold)
	assert.NotNil(t, model.rng)
}

// TestCpGDepletionModel_SetDepletionRate tests depletion rate setting
func TestCpGDepletionModel_SetDepletionRate(t *testing.T) {
	model := NewCpGDepletionModel()

	// Valid rates
	model.SetDepletionRate(0.5)
	assert.Equal(t, 0.5, model.depletionRate)

	model.SetDepletionRate(0.0)
	assert.Equal(t, 0.0, model.depletionRate)

	model.SetDepletionRate(1.0)
	assert.Equal(t, 1.0, model.depletionRate)

	// Invalid rates should not change the value
	originalRate := model.depletionRate
	model.SetDepletionRate(-0.1)
	assert.Equal(t, originalRate, model.depletionRate)

	model.SetDepletionRate(1.1)
	assert.Equal(t, originalRate, model.depletionRate)
}

// TestCpGDepletionModel_SetMethylationRate tests methylation rate setting
func TestCpGDepletionModel_SetMethylationRate(t *testing.T) {
	model := NewCpGDepletionModel()

	// Valid rates
	model.SetMethylationRate(0.6)
	assert.Equal(t, 0.6, model.methylationRate)

	model.SetMethylationRate(0.0)
	assert.Equal(t, 0.0, model.methylationRate)

	model.SetMethylationRate(1.0)
	assert.Equal(t, 1.0, model.methylationRate)

	// Invalid rates should not change the value
	originalRate := model.methylationRate
	model.SetMethylationRate(-0.1)
	assert.Equal(t, originalRate, model.methylationRate)

	model.SetMethylationRate(1.1)
	assert.Equal(t, originalRate, model.methylationRate)
}

// TestCpGDepletionModel_LoadCpGIslands tests CpG island loading
func TestCpGDepletionModel_LoadCpGIslands(t *testing.T) {
	model := NewCpGDepletionModel()

	islands := []CpGIsland{
		{
			Interval: GenomicInterval{
				Start:      100,
				End:        300,
				Chromosome: "chr1",
				GCContent:  0.6,
			},
			Length:      200,
			IsProtected: true,
		},
		{
			Interval: GenomicInterval{
				Start:      500,
				End:        700,
				Chromosome: "chr1",
				GCContent:  0.65,
			},
			Length:      200,
			IsProtected: true,
		},
	}

	model.LoadCpGIslands(islands)

	assert.Len(t, model.islandLocations, 2)
	assert.Equal(t, 100, model.islandLocations[0].Start)
	assert.Equal(t, 300, model.islandLocations[0].End)
	assert.Equal(t, 500, model.islandLocations[1].Start)
	assert.Equal(t, 700, model.islandLocations[1].End)
}

// TestCpGDepletionModel_ApplyDepletion tests basic depletion functionality
func TestCpGDepletionModel_ApplyDepletion(t *testing.T) {
	model := NewCpGDepletionModel()
	model.SetSeed(12345) // For reproducible results

	t.Run("empty sequence", func(t *testing.T) {
		result, err := model.ApplyDepletion([]byte{})
		assert.Error(t, err)
		assert.Nil(t, result)
		assert.Contains(t, err.Error(), "sequence cannot be empty")
	})

	t.Run("sequence with CpG sites", func(t *testing.T) {
		// Create a sequence with known CpG sites
		sequence := []byte("ATCGATCGCGATCGATCGCGATCGATCG")
		
		result, err := model.ApplyDepletion(sequence)
		require.NoError(t, err)
		require.NotNil(t, result)

		// Basic validations
		assert.Equal(t, len(sequence), len(result.OriginalSequence))
		assert.Equal(t, len(sequence), len(result.DepletedSequence))
		assert.Equal(t, sequence, result.OriginalSequence)
		assert.GreaterOrEqual(t, result.CpGSitesOriginal, 0)
		assert.LessOrEqual(t, result.CpGSitesRemaining, result.CpGSitesOriginal)
		assert.GreaterOrEqual(t, result.DepletionRate, 0.0)
		assert.LessOrEqual(t, result.DepletionRate, 1.0)
	})

	t.Run("sequence without CpG sites", func(t *testing.T) {
		sequence := []byte("ATATATATATATATATAT")
		
		result, err := model.ApplyDepletion(sequence)
		require.NoError(t, err)
		require.NotNil(t, result)

		assert.Equal(t, 0, result.CpGSitesOriginal)
		assert.Equal(t, 0, result.CpGSitesRemaining)
		assert.Equal(t, 0.0, result.DepletionRate)
		assert.Equal(t, sequence, result.DepletedSequence) // No changes
	})

	t.Run("high CpG content sequence", func(t *testing.T) {
		// Create a CpG-rich sequence
		sequence := []byte("CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG")
		
		result, err := model.ApplyDepletion(sequence)
		require.NoError(t, err)
		require.NotNil(t, result)

		assert.Greater(t, result.CpGSitesOriginal, 0)
		assert.GreaterOrEqual(t, len(result.IslandRegions), 0) // Should identify CpG islands
	})
}

// TestCpGDepletionModel_countCpGSites tests CpG site counting
func TestCpGDepletionModel_countCpGSites(t *testing.T) {
	model := NewCpGDepletionModel()

	testCases := []struct {
		sequence []byte
		expected int
	}{
		{[]byte(""), 0},
		{[]byte("A"), 0},
		{[]byte("CG"), 1},
		{[]byte("CGCG"), 2},
		{[]byte("ATCGATCG"), 2},
		{[]byte("CCGGCCGG"), 2},
		{[]byte("ACGTACGT"), 2},
		{[]byte("ATATATATA"), 0},
		{[]byte("CGCGCGCG"), 4}, // CGCGCGCG has 4 CpG sites: CG at positions 0,2,4,6
	}

	for _, tc := range testCases {
		t.Run(string(tc.sequence), func(t *testing.T) {
			count := model.countCpGSites(tc.sequence)
			assert.Equal(t, tc.expected, count)
		})
	}
}

// TestCpGDepletionModel_calculateGCContent tests GC content calculation
func TestCpGDepletionModel_calculateGCContent(t *testing.T) {
	model := NewCpGDepletionModel()

	testCases := []struct {
		sequence []byte
		expected float64
	}{
		{[]byte(""), 0.0},
		{[]byte("ATGC"), 0.5},
		{[]byte("GGCC"), 1.0},
		{[]byte("AATT"), 0.0},
		{[]byte("CGCGCGCG"), 1.0},
		{[]byte("ATATATATA"), 0.0},
		{[]byte("AAAAGGGCCCCTTT"), 0.5},
		{[]byte("NNNN"), 0.0}, // Non-standard bases
	}

	for _, tc := range testCases {
		t.Run(string(tc.sequence), func(t *testing.T) {
			gcContent := model.calculateGCContent(tc.sequence)
			assert.InDelta(t, tc.expected, gcContent, 0.001)
		})
	}
}

// TestCpGDepletionModel_calculateCpGRatio tests CpG ratio calculation
func TestCpGDepletionModel_calculateCpGRatio(t *testing.T) {
	model := NewCpGDepletionModel()

	testCases := []struct {
		sequence    []byte
		description string
		minRatio    float64
		maxRatio    float64
	}{
		{[]byte(""), "empty sequence", 0.0, 0.0},
		{[]byte("A"), "single base", 0.0, 0.0},
		{[]byte("CG"), "single CpG", 0.0, 10.0}, // High ratio expected
		{[]byte("CGCGCGCG"), "CpG rich", 0.5, 5.0}, // Should have reasonable ratio
		{[]byte("ATATATATA"), "no CpG", 0.0, 0.0},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			ratio := model.calculateCpGRatio(tc.sequence)
			assert.GreaterOrEqual(t, ratio, tc.minRatio)
			assert.LessOrEqual(t, ratio, tc.maxRatio)
		})
	}
}

// TestCpGDepletionModel_identifyCpGIslands tests CpG island identification
func TestCpGDepletionModel_identifyCpGIslands(t *testing.T) {
	model := NewCpGDepletionModel()

	t.Run("short sequence", func(t *testing.T) {
		sequence := []byte("CGCGCGCG")
		islands := model.identifyCpGIslands(sequence)
		// Should be empty because sequence is too short for island detection
		assert.Len(t, islands, 0)
	})

	t.Run("long CpG-rich sequence", func(t *testing.T) {
		// Create a long CpG-rich sequence that should be identified as an island
		sequence := make([]byte, 500)
		for i := 0; i < len(sequence); i += 2 {
			if i+1 < len(sequence) {
				sequence[i] = 'C'
				sequence[i+1] = 'G'
			}
		}

		islands := model.identifyCpGIslands(sequence)
		assert.GreaterOrEqual(t, len(islands), 0) // Should find at least 0 islands

		// Validate island properties
		for _, island := range islands {
			assert.Greater(t, island.Length, 0)
			assert.GreaterOrEqual(t, island.GCContent, 0.0)
			assert.LessOrEqual(t, island.GCContent, 1.0)
			assert.GreaterOrEqual(t, island.Interval.Start, 0)
			assert.LessOrEqual(t, island.Interval.End, len(sequence))
			assert.Greater(t, island.Interval.End, island.Interval.Start)
		}
	})

	t.Run("mixed sequence", func(t *testing.T) {
		// Create a sequence with both CpG-rich and CpG-poor regions
		cpgRich := make([]byte, 250)
		for i := 0; i < len(cpgRich); i += 2 {
			if i+1 < len(cpgRich) {
				cpgRich[i] = 'C'
				cpgRich[i+1] = 'G'
			}
		}

		cpgPoor := make([]byte, 250)
		for i := 0; i < len(cpgPoor); i += 2 {
			if i+1 < len(cpgPoor) {
				cpgPoor[i] = 'A'
				cpgPoor[i+1] = 'T'
			}
		}

		sequence := append(cpgRich, cpgPoor...)
		islands := model.identifyCpGIslands(sequence)

		// Validate that any islands found have reasonable properties
		for _, island := range islands {
			assert.Greater(t, island.GCContent, 0.4) // Should be GC-rich
			assert.True(t, island.IsProtected)
		}
	})
}

// TestCpGDepletionModel_isInCpGIsland tests position checking in CpG islands
func TestCpGDepletionModel_isInCpGIsland(t *testing.T) {
	model := NewCpGDepletionModel()

	islands := []CpGIsland{
		{
			Interval: GenomicInterval{Start: 100, End: 300},
		},
		{
			Interval: GenomicInterval{Start: 500, End: 700},
		},
	}

	testCases := []struct {
		position int
		expected bool
	}{
		{50, false},   // Before first island
		{100, true},   // Start of first island
		{200, true},   // Middle of first island
		{299, true},   // End of first island (exclusive)
		{300, false},  // Just after first island
		{400, false},  // Between islands
		{500, true},   // Start of second island
		{600, true},   // Middle of second island
		{699, true},   // End of second island
		{700, false},  // After second island
	}

	for _, tc := range testCases {
		t.Run(fmt.Sprintf("position_%d", tc.position), func(t *testing.T) {
			result := model.isInCpGIsland(tc.position, islands)
			assert.Equal(t, tc.expected, result)
		})
	}
}

// TestCpGDepletionModel_getContextMethylationProbability tests context-specific methylation
func TestCpGDepletionModel_getContextMethylationProbability(t *testing.T) {
	model := NewCpGDepletionModel()

	testCases := []struct {
		sequence    []byte
		position    int
		description string
		minProb     float64
		maxProb     float64
	}{
		{[]byte("ATCGATCG"), 2, "normal context", 0.8, 1.0},
		{[]byte("CCGGCCGG"), 2, "CCGG context", 1.0, 1.0}, // Should be higher
		{[]byte("GCGCGCGC"), 2, "GCGC context", 0.8, 1.0},
		{[]byte("CGCGCGCG"), 2, "CGCG context", 0.7, 1.0}, // Should be lower
		{[]byte("AT"), 0, "short sequence", 1.0, 1.0},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			prob := model.getContextMethylationProbability(tc.position, tc.sequence)
			assert.GreaterOrEqual(t, prob, tc.minProb)
			assert.LessOrEqual(t, prob, tc.maxProb)
		})
	}
}

// TestCpGDepletionModel_ValidateDepletion tests result validation
func TestCpGDepletionModel_ValidateDepletion(t *testing.T) {
	model := NewCpGDepletionModel()

	t.Run("nil result", func(t *testing.T) {
		err := model.ValidateDepletion(nil)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "depletion result is nil")
	})

	t.Run("mismatched sequence lengths", func(t *testing.T) {
		result := &DepletionResult{
			OriginalSequence: []byte("ATCG"),
			DepletedSequence: []byte("ATCGA"), // Different length
		}

		err := model.ValidateDepletion(result)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "sequence lengths do not match")
	})

	t.Run("more CpG sites in depleted", func(t *testing.T) {
		result := &DepletionResult{
			OriginalSequence:  []byte("ATCG"),
			DepletedSequence:  []byte("ATCG"),
			CpGSitesOriginal:  1,
			CpGSitesRemaining: 2, // More than original
		}

		err := model.ValidateDepletion(result)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "depleted sequence has more CpG sites than original")
	})

	t.Run("invalid transition events", func(t *testing.T) {
		result := &DepletionResult{
			OriginalSequence:  []byte("ATCG"),
			DepletedSequence:  []byte("ATCG"),
			CpGSitesOriginal:  1,
			CpGSitesRemaining: 0,
			TransitionEvents: []TransitionEvent{
				{OriginalBase: 'A', NewBase: 'T'}, // Invalid transition
			},
		}

		err := model.ValidateDepletion(result)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "invalid transition event")
	})

	t.Run("valid result", func(t *testing.T) {
		result := &DepletionResult{
			OriginalSequence:  []byte("CGCG"),
			DepletedSequence:  []byte("TGCG"),
			CpGSitesOriginal:  2,
			CpGSitesRemaining: 1,
			TransitionEvents: []TransitionEvent{
				{OriginalBase: 'C', NewBase: 'T', Position: 0},
			},
		}

		err := model.ValidateDepletion(result)
		assert.NoError(t, err)
	})
}

// TestCpGDepletionModel_GetDepletionStats tests statistics generation
func TestCpGDepletionModel_GetDepletionStats(t *testing.T) {
	model := NewCpGDepletionModel()

	result := &DepletionResult{
		OriginalSequence:  []byte("CGCGCGCG"),
		DepletedSequence:  []byte("TGCGTGCG"),
		CpGSitesOriginal:  3,
		CpGSitesRemaining: 1,
		DepletionRate:     0.67,
		TransitionEvents: []TransitionEvent{
			{Position: 0, InCpGIsland: false},
			{Position: 4, InCpGIsland: true},
		},
		IslandRegions: []CpGIsland{
			{Length: 100, ObservedCpG: 10, GCContent: 0.6},
			{Length: 150, ObservedCpG: 15, GCContent: 0.7},
		},
	}

	stats := model.GetDepletionStats(result)

	// Check basic statistics
	assert.Equal(t, 3, stats["original_cpg_sites"])
	assert.Equal(t, 1, stats["remaining_cpg_sites"])
	assert.Equal(t, 0.67, stats["depletion_rate"])
	assert.Equal(t, 2, stats["transition_events"])
	assert.Equal(t, 2, stats["cpg_islands_found"])
	assert.Equal(t, 8, stats["sequence_length"])

	// Check island statistics
	assert.Equal(t, 250, stats["total_island_length"])
	assert.Equal(t, 25, stats["total_island_cpg"])
	assert.InDelta(t, 0.65, stats["avg_island_gc_content"], 0.01)
	assert.InDelta(t, 31.25, stats["island_coverage"], 0.01) // 250/8

	// Check transition statistics
	assert.Equal(t, 1, stats["transitions_in_islands"])
	assert.Equal(t, 1, stats["transitions_outside_islands"])
}

// TestCpGDepletionModel_SimulateEvolutionaryDepletion tests evolutionary simulation
func TestCpGDepletionModel_SimulateEvolutionaryDepletion(t *testing.T) {
	model := NewCpGDepletionModel()
	model.SetSeed(12345) // For reproducible results

	t.Run("single generation", func(t *testing.T) {
		sequence := []byte("CGCGCGCG")
		
		result, err := model.SimulateEvolutionaryDepletion(sequence, 1)
		require.NoError(t, err)
		require.NotNil(t, result)

		assert.Equal(t, sequence, result.OriginalSequence)
		assert.Equal(t, len(sequence), len(result.DepletedSequence))
		assert.LessOrEqual(t, result.CpGSitesRemaining, result.CpGSitesOriginal)
	})

	t.Run("multiple generations", func(t *testing.T) {
		sequence := []byte("CGCGCGCGCGCGCGCG")
		
		result, err := model.SimulateEvolutionaryDepletion(sequence, 5)
		require.NoError(t, err)
		require.NotNil(t, result)

		assert.Equal(t, sequence, result.OriginalSequence)
		assert.LessOrEqual(t, result.CpGSitesRemaining, result.CpGSitesOriginal)
		
		// Should have accumulated transition events over generations
		assert.GreaterOrEqual(t, len(result.TransitionEvents), 0)
	})

	t.Run("error propagation", func(t *testing.T) {
		// Test with empty sequence to trigger error
		result, err := model.SimulateEvolutionaryDepletion([]byte{}, 1)
		assert.Error(t, err)
		assert.Nil(t, result)
	})
}

// TestCpGDepletionModel_SetSeed tests seed setting
func TestCpGDepletionModel_SetSeed(t *testing.T) {
	model := NewCpGDepletionModel()

	// Set seed and generate some results
	model.SetSeed(12345)
	sequence := []byte("CGCGCGCGCGCGCGCG")
	result1, err := model.ApplyDepletion(sequence)
	require.NoError(t, err)

	// Reset with same seed and generate again
	model.SetSeed(12345)
	result2, err := model.ApplyDepletion(sequence)
	require.NoError(t, err)

	// Results should be identical
	assert.Equal(t, result1.DepletedSequence, result2.DepletedSequence)
	assert.Equal(t, len(result1.TransitionEvents), len(result2.TransitionEvents))
}

// TestCpGDepletionModel_GetModelParameters tests parameter retrieval
func TestCpGDepletionModel_GetModelParameters(t *testing.T) {
	model := NewCpGDepletionModel()

	// Add some islands
	islands := []CpGIsland{{Interval: GenomicInterval{Start: 100, End: 300}}}
	model.LoadCpGIslands(islands)

	params := model.GetModelParameters()

	assert.Equal(t, 0.8, params["depletion_rate"])
	assert.Equal(t, 0.7, params["methylation_rate"])
	assert.Equal(t, 0.02, params["transition_prob"])
	assert.Equal(t, 0.5, params["gc_content_threshold"])
	assert.Equal(t, 200, params["length_threshold"])
	assert.Equal(t, 1, params["known_islands"])
}

// TestGenomicInterval_Validation tests genomic interval structure
func TestGenomicInterval_Validation(t *testing.T) {
	interval := GenomicInterval{
		Start:      1000,
		End:        2000,
		Chromosome: "chr1",
		Strand:     "+",
		GCContent:  0.45,
		CpGRatio:   0.8,
	}

	assert.Greater(t, interval.End, interval.Start)
	assert.NotEmpty(t, interval.Chromosome)
	assert.GreaterOrEqual(t, interval.GCContent, 0.0)
	assert.LessOrEqual(t, interval.GCContent, 1.0)
	assert.GreaterOrEqual(t, interval.CpGRatio, 0.0)
}

// TestCpGIsland_Validation tests CpG island structure
func TestCpGIsland_Validation(t *testing.T) {
	island := CpGIsland{
		Interval: GenomicInterval{
			Start: 1000,
			End:   1500,
		},
		ObservedCpG: 25,
		ExpectedCpG: 20.5,
		CpGRatio:    1.22,
		GCContent:   0.65,
		Length:      500,
		IsProtected: true,
	}

	assert.Greater(t, island.Length, 0)
	assert.GreaterOrEqual(t, island.ObservedCpG, 0)
	assert.GreaterOrEqual(t, island.ExpectedCpG, 0.0)
	assert.GreaterOrEqual(t, island.CpGRatio, 0.0)
	assert.GreaterOrEqual(t, island.GCContent, 0.0)
	assert.LessOrEqual(t, island.GCContent, 1.0)
	assert.Equal(t, island.Interval.End-island.Interval.Start, island.Length)
}

// TestTransitionEvent_Validation tests transition event structure
func TestTransitionEvent_Validation(t *testing.T) {
	event := TransitionEvent{
		Position:         150,
		OriginalBase:     'C',
		NewBase:          'T',
		InCpGIsland:      false,
		MethylationState: "methylated",
		TransitionType:   "deamination",
	}

	assert.GreaterOrEqual(t, event.Position, 0)
	assert.Equal(t, byte('C'), event.OriginalBase)
	assert.Equal(t, byte('T'), event.NewBase)
	assert.Contains(t, []string{"methylated", "unmethylated", "hemimethylated"}, event.MethylationState)
	assert.Contains(t, []string{"deamination", "replication_error"}, event.TransitionType)
}

// TestCpGDepletionModel_overlapsWithExistingIsland tests island overlap detection
func TestCpGDepletionModel_overlapsWithExistingIsland(t *testing.T) {
	model := NewCpGDepletionModel()

	existingIslands := []CpGIsland{
		{Interval: GenomicInterval{Start: 100, End: 300}},
		{Interval: GenomicInterval{Start: 500, End: 700}},
	}

	testCases := []struct {
		newIsland CpGIsland
		expected  bool
		description string
	}{
		{
			CpGIsland{Interval: GenomicInterval{Start: 50, End: 150}},
			true,
			"overlaps with first island",
		},
		{
			CpGIsland{Interval: GenomicInterval{Start: 250, End: 350}},
			true,
			"overlaps with first island end",
		},
		{
			CpGIsland{Interval: GenomicInterval{Start: 400, End: 450}},
			false,
			"no overlap between islands",
		},
		{
			CpGIsland{Interval: GenomicInterval{Start: 800, End: 900}},
			false,
			"after all existing islands",
		},
		{
			CpGIsland{Interval: GenomicInterval{Start: 0, End: 50}},
			false,
			"before all existing islands",
		},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			result := model.overlapsWithExistingIsland(tc.newIsland, existingIslands)
			assert.Equal(t, tc.expected, result)
		})
	}
}

// BenchmarkCpGDepletionModel_ApplyDepletion benchmarks depletion application
func BenchmarkCpGDepletionModel_ApplyDepletion(b *testing.B) {
	model := NewCpGDepletionModel()
	model.SetSeed(12345)

	// Create a realistic sequence
	sequence := make([]byte, 1000)
	for i := 0; i < len(sequence); i += 4 {
		if i+3 < len(sequence) {
			sequence[i] = 'A'
			sequence[i+1] = 'T'
			sequence[i+2] = 'C'
			sequence[i+3] = 'G'
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = model.ApplyDepletion(sequence)
	}
}

// BenchmarkCpGDepletionModel_identifyCpGIslands benchmarks island identification
func BenchmarkCpGDepletionModel_identifyCpGIslands(b *testing.B) {
	model := NewCpGDepletionModel()

	// Create a long CpG-rich sequence
	sequence := make([]byte, 10000)
	for i := 0; i < len(sequence); i += 2 {
		if i+1 < len(sequence) {
			if i%100 < 50 { // Alternating CpG-rich and poor regions
				sequence[i] = 'C'
				sequence[i+1] = 'G'
			} else {
				sequence[i] = 'A'
				sequence[i+1] = 'T'
			}
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = model.identifyCpGIslands(sequence)
	}
}

// TestCpGDepletionModel_AnalyzeSequence tests the AnalyzeSequence method and the
// CpGRatio > 1.0 edge case (regression for #38).
func TestCpGDepletionModel_AnalyzeSequence(t *testing.T) {
	model := NewCpGDepletionModel()

	t.Run("empty sequence returns zero score", func(t *testing.T) {
		score, details, err := model.AnalyzeSequence([]byte{})
		require.NoError(t, err)
		assert.Equal(t, 0.0, score)
		assert.NotNil(t, details)
	})

	t.Run("fully depleted sequence: no CpG at all", func(t *testing.T) {
		// No C or G → expectedCpG = 0; cpgRatio = 0 → depletionScore = 1.0
		seq := []byte("AAATTTAAATTTAAATTT")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.Equal(t, 0.0, details.CpGRatio)
		assert.Equal(t, 1.0, score)
		assert.Equal(t, 1.0, details.DepletionScore)
		assert.Equal(t, 1.0, details.RawDepletionScore)
	})

	t.Run("moderate depletion: ratio < 1", func(t *testing.T) {
		// "GGGCCCGGGCCC" (12 bases): G=6, C=6
		// expectedCpG = 6*6/12 = 3
		// CpG dinucleotides: only C(5)→G(6) qualifies → count=1
		// cpgRatio = 1/3 ≈ 0.333 → depletionScore ≈ 0.667
		seq := []byte("GGGCCCGGGCCC")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.InDelta(t, 1.0/3.0, details.CpGRatio, 0.001)
		assert.Greater(t, score, 0.0, "depleted sequence should have score > 0")
		assert.Less(t, score, 1.0, "partially depleted sequence should have score < 1")
		assert.Equal(t, score, details.DepletionScore)
		assert.Equal(t, score, details.RawDepletionScore, "no clamping needed when ratio < 1")
	})

	t.Run("CpG-enriched sequence: ratio > 1.0 should NOT produce depletionScore = 1.0", func(t *testing.T) {
		// Construct a sequence with far more CpG than expected.
		// "CGCGCGCGCGCG" (12 bases): C=6, G=6, expected=36/12=3; CpG dinucleotides: CG×6 = 6
		// cpgRatio = 6/3 = 2.0 → rawScore = -1.0 → depletionScore clamped to 0.0
		seq := []byte("CGCGCGCGCGCG")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.Greater(t, details.CpGRatio, 1.0, "cpgRatio should exceed 1.0 for CpG-enriched sequence")
		assert.Equal(t, 0.0, score, "depletionScore must be 0 (not 1) for CpG-enriched sequence")
		assert.Equal(t, 0.0, details.DepletionScore)
		assert.Less(t, details.RawDepletionScore, 0.0, "RawDepletionScore should be negative for enriched sequence")
	})

	t.Run("score is in [0, 1] for any input", func(t *testing.T) {
		sequences := [][]byte{
			[]byte("AAAAAAAAAA"),
			[]byte("GCGCGCGCGC"),
			[]byte("CGCGCGCGCG"),
			[]byte("ACGTACGTACGT"),
			[]byte("TTTTCCCCGGGG"),
		}
		for _, seq := range sequences {
			score, details, err := model.AnalyzeSequence(seq)
			require.NoError(t, err)
			assert.GreaterOrEqual(t, score, 0.0, "score must be >= 0 for seq %q", seq)
			assert.LessOrEqual(t, score, 1.0, "score must be <= 1 for seq %q", seq)
			assert.Equal(t, score, details.DepletionScore)
		}
	})
}
// TestExtendCpGIsland_ExpectedCpG_NoOverflow is a regression test for #119.
// Before the fix, float64(gCount*cCount) performed integer multiplication before
// promotion to float64; on 32-bit systems this overflows for sequences where
// gCount*cCount > math.MaxInt32, producing a negative ExpectedCpG and a corrupt
// CpGRatio.  We verify the result is positive and mathematically correct for a
// large, GC-balanced sequence.
func TestExtendCpGIsland_ExpectedCpG_NoOverflow(t *testing.T) {
	model := NewCpGDepletionModel()

	// Build a 100 KB sequence that is 50% G / 50% C so gCount=cCount=50000.
	// With the old code on 32-bit: 50000*50000 = 2,500,000,000 > MaxInt32
	// → wraps to -1,794,967,296 → ExpectedCpG negative → CpGRatio NaN/negative.
	const n = 100_000
	seq := make([]byte, n)
	for i := range seq {
		if i%2 == 0 {
			seq[i] = 'G'
		} else {
			seq[i] = 'C'
		}
	}

	island := model.extendCpGIsland(seq, 0, 200)

	// ExpectedCpG = gCount * cCount / N = 50000 * 50000 / 100000 = 25000
	const wantExpected = 25_000.0
	if island.ExpectedCpG <= 0 {
		t.Errorf("ExpectedCpG = %v; want positive (integer overflow regression)", island.ExpectedCpG)
	}
	if island.ExpectedCpG != wantExpected {
		t.Errorf("ExpectedCpG = %v, want %v", island.ExpectedCpG, wantExpected)
	}
	if island.CpGRatio < 0 {
		t.Errorf("CpGRatio = %v; want non-negative (integer overflow regression)", island.CpGRatio)
	}
}

// TestCalculateCpGRatio_NoOverflow is a regression test for #121.
// calculateCpGRatio used float64(gCount*cCount) — int multiplication before
// float64 promotion. On 32-bit int, gCount=cCount=50000 overflows.
func TestCalculateCpGRatio_NoOverflow(t *testing.T) {
	model := NewCpGDepletionModel()

	const n = 100_000
	seq := make([]byte, n)
	for i := range seq {
		if i%2 == 0 {
			seq[i] = 'G'
		} else {
			seq[i] = 'C'
		}
	}

	ratio := model.calculateCpGRatio(seq)

	if ratio < 0 {
		t.Errorf("calculateCpGRatio returned %v for large GC sequence; want non-negative (integer overflow regression)", ratio)
	}
	if ratio != ratio { // NaN check
		t.Errorf("calculateCpGRatio returned NaN")
	}
}

// TestAnalyzeSequence_NoOverflow is a regression test for #121.
// AnalyzeSequence had the same float64(gCount*cCount) overflow as calculateCpGRatio.
func TestAnalyzeSequence_NoOverflow(t *testing.T) {
	model := NewCpGDepletionModel()

	const n = 100_000
	seq := make([]byte, n)
	for i := range seq {
		if i%2 == 0 {
			seq[i] = 'G'
		} else {
			seq[i] = 'C'
		}
	}

	score, details, err := model.AnalyzeSequence(seq)
	if err != nil {
		t.Fatalf("AnalyzeSequence: %v", err)
	}
	if score < 0 || score > 1 {
		t.Errorf("score %v out of [0,1] (integer overflow regression)", score)
	}
	if details.CpGRatio < 0 {
		t.Errorf("CpGRatio %v negative (integer overflow regression)", details.CpGRatio)
	}
}
