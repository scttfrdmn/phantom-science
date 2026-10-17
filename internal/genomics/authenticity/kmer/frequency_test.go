package kmer

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewFrequencyExtractor tests the creation of new frequency extractors
func TestNewFrequencyExtractor(t *testing.T) {
	extractor := NewFrequencyExtractor(4)
	
	assert.Equal(t, 4, extractor.K)
	assert.True(t, extractor.IgnoreN)
	assert.False(t, extractor.CaseSensitive)
	assert.False(t, extractor.CircularSequence)
	assert.Equal(t, 1, extractor.MinFrequency)
}

// TestFrequencyExtractor_ExtractFromSequence tests k-mer extraction from single sequences
func TestFrequencyExtractor_ExtractFromSequence(t *testing.T) {
	extractor := NewFrequencyExtractor(3)
	
	t.Run("basic extraction", func(t *testing.T) {
		sequence := []byte("ATCGATCG")
		freqTable, err := extractor.ExtractFromSequence(sequence)
		
		require.NoError(t, err)
		assert.Equal(t, 3, freqTable.K)
		assert.Equal(t, 6, freqTable.TotalKmers) // 8 - 3 + 1 = 6
		
		// Expected k-mers: ATC, TCG, CGA, GAT, ATC, TCG
		assert.Equal(t, 2, freqTable.Frequencies["ATC"])
		assert.Equal(t, 2, freqTable.Frequencies["TCG"])
		assert.Equal(t, 1, freqTable.Frequencies["CGA"])
		assert.Equal(t, 1, freqTable.Frequencies["GAT"])
	})
	
	t.Run("case sensitivity", func(t *testing.T) {
		extractor.CaseSensitive = true
		sequence := []byte("ATCatc")
		freqTable, err := extractor.ExtractFromSequence(sequence)
		
		require.NoError(t, err)
		
		// With case sensitivity, only valid (uppercase) k-mers are kept
		// "ATCatc" produces k-mers: "ATC", "TCa", "Cat", "atc"
		// But isValidKmer only accepts A,T,G,C so lowercase are filtered out
		assert.Equal(t, 1, freqTable.Frequencies["ATC"])
		assert.Equal(t, 1, len(freqTable.Frequencies)) // Only "ATC" is valid
		
		extractor.CaseSensitive = false // Reset
		
		// Without case sensitivity, all converted to uppercase
		freqTable2, err := extractor.ExtractFromSequence(sequence)
		require.NoError(t, err)
		assert.Equal(t, 2, freqTable2.Frequencies["ATC"]) // "ATC" appears twice after conversion
		assert.Equal(t, 1, freqTable2.Frequencies["TCA"]) // "TCa" becomes "TCA"  
		assert.Equal(t, 1, freqTable2.Frequencies["CAT"]) // "Cat" becomes "CAT"
	})
	
	t.Run("ignore N bases", func(t *testing.T) {
		sequence := []byte("ATCNGATCG")
		freqTable, err := extractor.ExtractFromSequence(sequence)
		
		require.NoError(t, err)
		// K-mers containing N should be skipped
		assert.NotContains(t, freqTable.Frequencies, "TCN")
		assert.NotContains(t, freqTable.Frequencies, "CNG")
		assert.NotContains(t, freqTable.Frequencies, "NGA")
	})
	
	t.Run("circular sequence", func(t *testing.T) {
		extractor.CircularSequence = true
		sequence := []byte("ATCG")
		freqTable, err := extractor.ExtractFromSequence(sequence)
		
		require.NoError(t, err)
		assert.Equal(t, 4, freqTable.TotalKmers) // Should be same as sequence length
		assert.Contains(t, freqTable.Frequencies, "GAT") // Wraps around: G + AT
		
		extractor.CircularSequence = false // Reset
	})
	
	t.Run("minimum frequency filter", func(t *testing.T) {
		extractor.MinFrequency = 2
		sequence := []byte("ATCGATCG") // ATC and TCG appear twice, others once
		freqTable, err := extractor.ExtractFromSequence(sequence)
		
		require.NoError(t, err)
		assert.Equal(t, 2, len(freqTable.Frequencies)) // Only ATC and TCG
		assert.Equal(t, 2, freqTable.Frequencies["ATC"])
		assert.Equal(t, 2, freqTable.Frequencies["TCG"])
		
		extractor.MinFrequency = 1 // Reset
	})
	
	t.Run("error cases", func(t *testing.T) {
		// Sequence too short
		_, err := extractor.ExtractFromSequence([]byte("AT"))
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "shorter than k-mer size")
	})
}

// TestFrequencyExtractor_ExtractFromMultipleSequences tests extraction from multiple sequences
func TestFrequencyExtractor_ExtractFromMultipleSequences(t *testing.T) {
	extractor := NewFrequencyExtractor(3)
	
	t.Run("multiple sequences", func(t *testing.T) {
		sequences := [][]byte{
			[]byte("ATCG"),
			[]byte("GATC"),
		}
		
		freqTable, err := extractor.ExtractFromMultipleSequences(sequences)
		require.NoError(t, err)
		
		// First sequence: ATC, TCG
		// Second sequence: GAT, ATC  
		assert.Equal(t, 2, freqTable.Frequencies["ATC"]) // Appears in both
		assert.Equal(t, 1, freqTable.Frequencies["TCG"])
		assert.Equal(t, 1, freqTable.Frequencies["GAT"])
		assert.Equal(t, 4, freqTable.TotalCount)
		assert.Equal(t, 4, freqTable.TotalKmers) // 2 + 2
	})
	
	t.Run("error cases", func(t *testing.T) {
		// No sequences
		_, err := extractor.ExtractFromMultipleSequences([][]byte{})
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "no sequences provided")
		
		// One sequence too short
		sequences := [][]byte{
			[]byte("ATCG"),
			[]byte("AT"), // Too short
		}
		_, err = extractor.ExtractFromMultipleSequences(sequences)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "shorter than k-mer size")
	})
}

// TestFrequencyExtractor_isValidKmer tests k-mer validation
func TestFrequencyExtractor_isValidKmer(t *testing.T) {
	extractor := NewFrequencyExtractor(3)
	
	assert.True(t, extractor.isValidKmer("ATC"))
	assert.True(t, extractor.isValidKmer("GGGG"))
	assert.False(t, extractor.isValidKmer("ATCX")) // Invalid base
	assert.True(t, extractor.isValidKmer(""))      // Empty is considered valid
	
	// Test N handling
	extractor.IgnoreN = false
	assert.True(t, extractor.isValidKmer("ATCN"))
	
	extractor.IgnoreN = true
	assert.False(t, extractor.isValidKmer("ATCN"))
}

// TestFrequencyTable_GetProbabilityDistribution tests probability calculation
func TestFrequencyTable_GetProbabilityDistribution(t *testing.T) {
	freqTable := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 4,
			"TC": 2,
			"CG": 2,
		},
		TotalCount: 8,
		TotalKmers: 8,
	}
	
	probs := freqTable.GetProbabilityDistribution()
	
	assert.InDelta(t, 0.5, probs["AT"], 1e-9)  // 4/8
	assert.InDelta(t, 0.25, probs["TC"], 1e-9) // 2/8
	assert.InDelta(t, 0.25, probs["CG"], 1e-9) // 2/8
	
	// Test empty table
	emptyTable := &FrequencyTable{TotalCount: 0}
	emptyProbs := emptyTable.GetProbabilityDistribution()
	assert.Empty(t, emptyProbs)
}

// TestFrequencyTable_GetTopKmers tests top k-mer retrieval
func TestFrequencyTable_GetTopKmers(t *testing.T) {
	freqTable := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 10,
			"TC": 5,
			"CG": 3,
			"GG": 1,
		},
		TotalCount: 19,
		TotalKmers: 19,
	}
	
	top3 := freqTable.GetTopKmers(3)
	
	require.Len(t, top3, 3)
	assert.Equal(t, "AT", top3[0].Kmer)
	assert.Equal(t, 10, top3[0].Count)
	assert.InDelta(t, 10.0/19.0, top3[0].Frequency, 1e-9)
	
	assert.Equal(t, "TC", top3[1].Kmer)
	assert.Equal(t, 5, top3[1].Count)
	
	assert.Equal(t, "CG", top3[2].Kmer)
	assert.Equal(t, 3, top3[2].Count)
	
	// Test requesting more than available
	all := freqTable.GetTopKmers(10)
	assert.Len(t, all, 4) // Only 4 k-mers available
}

// TestGetExpectedKmers tests expected k-mer generation
func TestGetExpectedKmers(t *testing.T) {
	t.Run("k=1", func(t *testing.T) {
		kmers := GetExpectedKmers(1)
		expected := []string{"A", "T", "G", "C"}
		assert.ElementsMatch(t, expected, kmers)
	})
	
	t.Run("k=2", func(t *testing.T) {
		kmers := GetExpectedKmers(2)
		assert.Len(t, kmers, 16) // 4^2
		assert.Contains(t, kmers, "AA")
		assert.Contains(t, kmers, "AT")
		assert.Contains(t, kmers, "GC")
		assert.Contains(t, kmers, "TT")
	})
	
	t.Run("k=0", func(t *testing.T) {
		kmers := GetExpectedKmers(0)
		assert.Nil(t, kmers)
	})
	
	t.Run("k=3", func(t *testing.T) {
		kmers := GetExpectedKmers(3)
		assert.Len(t, kmers, 64) // 4^3
	})
}

// TestFrequencyTable_GetKmerComplexity tests complexity calculation
func TestFrequencyTable_GetKmerComplexity(t *testing.T) {
	t.Run("full complexity", func(t *testing.T) {
		// All possible 2-mers present
		freqTable := &FrequencyTable{
			K: 2,
			Frequencies: map[string]int{
				"AA": 1, "AT": 1, "AG": 1, "AC": 1,
				"TA": 1, "TT": 1, "TG": 1, "TC": 1,
				"GA": 1, "GT": 1, "GG": 1, "GC": 1,
				"CA": 1, "CT": 1, "CG": 1, "CC": 1,
			},
			TotalCount: 16,
		}
		
		complexity := freqTable.GetKmerComplexity()
		assert.InDelta(t, 1.0, complexity, 1e-9) // 16/16 = 1.0
	})
	
	t.Run("partial complexity", func(t *testing.T) {
		freqTable := &FrequencyTable{
			K: 2,
			Frequencies: map[string]int{
				"AT": 2,
				"GC": 2,
			},
			TotalCount: 4,
		}
		
		complexity := freqTable.GetKmerComplexity()
		assert.InDelta(t, 0.125, complexity, 1e-9) // 2/16 = 0.125
	})
	
	t.Run("empty table", func(t *testing.T) {
		freqTable := &FrequencyTable{TotalCount: 0}
		complexity := freqTable.GetKmerComplexity()
		assert.Equal(t, 0.0, complexity)
	})
}

// TestFrequencyTable_GetEntropyShannon tests Shannon entropy calculation
func TestFrequencyTable_GetEntropyShannon(t *testing.T) {
	t.Run("uniform distribution", func(t *testing.T) {
		freqTable := &FrequencyTable{
			K: 2,
			Frequencies: map[string]int{
				"AT": 2,
				"GC": 2,
				"TA": 2,
				"CG": 2,
			},
			TotalCount: 8,
		}
		
		entropy := freqTable.GetEntropyShannon()
		expectedEntropy := 2.0 // log2(4) for uniform distribution of 4 items
		assert.InDelta(t, expectedEntropy, entropy, 1e-9)
	})
	
	t.Run("skewed distribution", func(t *testing.T) {
		freqTable := &FrequencyTable{
			K: 2,
			Frequencies: map[string]int{
				"AT": 7,
				"GC": 1,
			},
			TotalCount: 8,
		}
		
		entropy := freqTable.GetEntropyShannon()
		// Calculate expected: -(7/8 * log2(7/8) + 1/8 * log2(1/8))
		p1 := 7.0 / 8.0
		p2 := 1.0 / 8.0
		expected := -(p1*math.Log2(p1) + p2*math.Log2(p2))
		assert.InDelta(t, expected, entropy, 1e-9)
	})
	
	t.Run("empty table", func(t *testing.T) {
		freqTable := &FrequencyTable{TotalCount: 0}
		entropy := freqTable.GetEntropyShannon()
		assert.Equal(t, 0.0, entropy)
	})
}

// TestGetMaxEntropy tests maximum entropy calculation
func TestGetMaxEntropy(t *testing.T) {
	assert.Equal(t, 2.0, GetMaxEntropy(1))  // 1 * log2(4) = 1 * 2
	assert.Equal(t, 4.0, GetMaxEntropy(2))  // 2 * log2(4) = 2 * 2
	assert.Equal(t, 6.0, GetMaxEntropy(3))  // 3 * log2(4) = 3 * 2
}

// TestFrequencyTable_GetNormalizedEntropy tests normalized entropy
func TestFrequencyTable_GetNormalizedEntropy(t *testing.T) {
	freqTable := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 2,
			"GC": 2,
			"TA": 2,
			"CG": 2,
		},
		TotalCount: 8,
	}
	
	normalizedEntropy := freqTable.GetNormalizedEntropy()
	// Shannon entropy is 2.0, max entropy for k=2 is 4.0
	assert.InDelta(t, 0.5, normalizedEntropy, 1e-9)
}

// TestMergeFrequencyTables tests table merging
func TestMergeFrequencyTables(t *testing.T) {
	table1 := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 3,
			"GC": 2,
		},
		TotalCount: 5,
		TotalKmers: 5,
	}
	
	table2 := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 1,
			"TA": 3,
		},
		TotalCount: 4,
		TotalKmers: 4,
	}
	
	t.Run("successful merge", func(t *testing.T) {
		merged, err := MergeFrequencyTables(table1, table2)
		require.NoError(t, err)
		
		assert.Equal(t, 2, merged.K)
		assert.Equal(t, 4, merged.Frequencies["AT"]) // 3 + 1
		assert.Equal(t, 2, merged.Frequencies["GC"]) // 2 + 0
		assert.Equal(t, 3, merged.Frequencies["TA"]) // 0 + 3
		assert.Equal(t, 9, merged.TotalCount)        // 5 + 4
		assert.Equal(t, 9, merged.TotalKmers)        // 5 + 4
	})
	
	t.Run("error cases", func(t *testing.T) {
		// No tables
		_, err := MergeFrequencyTables()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "no frequency tables provided")
		
		// Mismatched k values
		table3 := &FrequencyTable{K: 3}
		_, err = MergeFrequencyTables(table1, table3)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "mismatched k values")
	})
}

// TestFrequencyTable_GetGCContent tests GC content calculation
func TestFrequencyTable_GetGCContent(t *testing.T) {
	freqTable := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"GC": 2, // 4 G/C bases
			"AT": 2, // 0 G/C bases
			"GA": 1, // 1 G/C base
		},
		TotalCount: 5,
	}
	
	gcContent := freqTable.GetGCContent()
	// Total bases: 5 * 2 = 10
	// G/C bases: 2*2 + 2*0 + 1*1 = 5
	// GC content: 5/10 = 0.5
	assert.InDelta(t, 0.5, gcContent, 1e-9)
	
	// Test empty table
	emptyTable := &FrequencyTable{TotalCount: 0}
	assert.Equal(t, 0.0, emptyTable.GetGCContent())
}

// TestFrequencyTable_GetKmerStats tests comprehensive statistics
func TestFrequencyTable_GetKmerStats(t *testing.T) {
	freqTable := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 10,
			"GC": 5,
			"TA": 3,
			"CG": 2,
		},
		TotalCount: 20,
		TotalKmers: 25,
	}
	
	stats := freqTable.GetKmerStats()
	
	assert.Equal(t, 2, stats.K)
	assert.Equal(t, 25, stats.TotalKmers)
	assert.Equal(t, 4, stats.UniqueKmers)
	assert.Equal(t, 20, stats.TotalCount)
	assert.Equal(t, 2, stats.MinCount)
	assert.Equal(t, 10, stats.MaxCount)
	assert.InDelta(t, 4.0, stats.MedianCount, 1e-9) // (3+5)/2
	assert.InDelta(t, 5.0, stats.MeanCount, 1e-9)  // (10+5+3+2)/4
	assert.Greater(t, stats.Complexity, 0.0)
	assert.Greater(t, stats.ShannonEntropy, 0.0)
	assert.Greater(t, stats.NormalizedEntropy, 0.0)
	assert.Greater(t, stats.GCContent, 0.0)
}

// TestFrequencyTable_FilterByFrequency tests frequency filtering
func TestFrequencyTable_FilterByFrequency(t *testing.T) {
	freqTable := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 10,
			"GC": 5,
			"TA": 3,
			"CG": 1,
		},
		TotalCount: 19,
		TotalKmers: 20,
	}
	
	filtered := freqTable.FilterByFrequency(3, 6)
	
	assert.Equal(t, 2, len(filtered.Frequencies))
	assert.Equal(t, 5, filtered.Frequencies["GC"])
	assert.Equal(t, 3, filtered.Frequencies["TA"])
	assert.Equal(t, 8, filtered.TotalCount) // 5 + 3
	assert.Equal(t, 20, filtered.TotalKmers) // Preserved
	assert.NotContains(t, filtered.Frequencies, "AT") // Too high
	assert.NotContains(t, filtered.Frequencies, "CG") // Too low
}

// TestFrequencyTable_GetRareKmers tests rare k-mer identification
func TestFrequencyTable_GetRareKmers(t *testing.T) {
	freqTable := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 10,
			"GC": 5,
			"TA": 2,
			"CG": 1,
		},
	}
	
	rare := freqTable.GetRareKmers(3)
	
	assert.Len(t, rare, 2)
	assert.Contains(t, rare, "CG") // count = 1 < 3
	assert.Contains(t, rare, "TA") // count = 2 < 3
	
	// Should be sorted
	assert.True(t, rare[0] <= rare[1])
}

// TestFrequencyTable_GetAbundantKmers tests abundant k-mer identification
func TestFrequencyTable_GetAbundantKmers(t *testing.T) {
	freqTable := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 10,
			"GC": 5,
			"TA": 2,
			"CG": 1,
		},
	}
	
	abundant := freqTable.GetAbundantKmers(4)
	
	assert.Len(t, abundant, 2)
	assert.Contains(t, abundant, "AT") // count = 10 > 4
	assert.Contains(t, abundant, "GC") // count = 5 > 4
	
	// Should be sorted
	assert.True(t, abundant[0] <= abundant[1])
}

// TestFrequencyTable_EdgeCases tests edge cases and boundary conditions
func TestFrequencyTable_EdgeCases(t *testing.T) {
	t.Run("single base sequences", func(t *testing.T) {
		extractor := NewFrequencyExtractor(1)
		sequence := []byte("AAAA")
		
		freqTable, err := extractor.ExtractFromSequence(sequence)
		require.NoError(t, err)
		
		assert.Equal(t, 1, len(freqTable.Frequencies))
		assert.Equal(t, 4, freqTable.Frequencies["A"])
		assert.Equal(t, 4, freqTable.TotalCount)
	})
	
	t.Run("exact k-mer size sequence", func(t *testing.T) {
		extractor := NewFrequencyExtractor(4)
		sequence := []byte("ATCG") // Exactly k=4 length
		
		freqTable, err := extractor.ExtractFromSequence(sequence)
		require.NoError(t, err)
		
		assert.Equal(t, 1, len(freqTable.Frequencies))
		assert.Equal(t, 1, freqTable.Frequencies["ATCG"])
		assert.Equal(t, 1, freqTable.TotalKmers)
	})
	
	t.Run("empty frequency table operations", func(t *testing.T) {
		emptyTable := &FrequencyTable{
			K: 3,
			Frequencies: make(map[string]int),
			TotalCount: 0,
		}
		
		// All operations should handle empty tables gracefully
		assert.Empty(t, emptyTable.GetProbabilityDistribution())
		assert.Empty(t, emptyTable.GetTopKmers(5))
		assert.Equal(t, 0.0, emptyTable.GetKmerComplexity())
		assert.Equal(t, 0.0, emptyTable.GetEntropyShannon())
		assert.Equal(t, 0.0, emptyTable.GetNormalizedEntropy())
		assert.Equal(t, 0.0, emptyTable.GetGCContent())
		assert.Empty(t, emptyTable.GetRareKmers(1))
		assert.Empty(t, emptyTable.GetAbundantKmers(1))
	})
}

// BenchmarkFrequencyExtractor benchmarks performance
func BenchmarkFrequencyExtractor(b *testing.B) {
	extractor := NewFrequencyExtractor(6)
	
	// Generate a realistic sequence (1MB)
	sequence := make([]byte, 1000000)
	bases := []byte{'A', 'T', 'G', 'C'}
	for i := range sequence {
		sequence[i] = bases[i%4]
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := extractor.ExtractFromSequence(sequence)
		if err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkKmerStats benchmarks statistics calculation
func BenchmarkKmerStats(b *testing.B) {
	// Create a large frequency table
	freqTable := &FrequencyTable{
		K: 6,
		Frequencies: make(map[string]int),
		TotalCount: 100000,
		TotalKmers: 100000,
	}
	
	// Add many k-mers
	kmers := GetExpectedKmers(6)
	for i, kmer := range kmers {
		freqTable.Frequencies[kmer] = i + 1
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = freqTable.GetKmerStats()
	}
}