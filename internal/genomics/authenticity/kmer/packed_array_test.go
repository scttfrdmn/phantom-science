package kmer

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewPackedFrequencyTable tests the constructor with valid and invalid inputs.
func TestNewPackedFrequencyTable(t *testing.T) {
	t.Run("valid_k=1", func(t *testing.T) {
		pft, err := NewPackedFrequencyTable(1)
		require.NoError(t, err)
		assert.Equal(t, 1, pft.K)
		assert.Equal(t, 4, pft.MaxKmers) // 4^1
		assert.Equal(t, 0, pft.TotalCount)
	})

	t.Run("valid_k=3", func(t *testing.T) {
		pft, err := NewPackedFrequencyTable(3)
		require.NoError(t, err)
		assert.Equal(t, 3, pft.K)
		assert.Equal(t, 64, pft.MaxKmers) // 4^3
	})

	t.Run("valid_k=4", func(t *testing.T) {
		pft, err := NewPackedFrequencyTable(4)
		require.NoError(t, err)
		assert.Equal(t, 256, pft.MaxKmers) // 4^4
	})

	t.Run("k=0_invalid", func(t *testing.T) {
		_, err := NewPackedFrequencyTable(0)
		assert.Error(t, err)
	})

	t.Run("k=17_invalid", func(t *testing.T) {
		_, err := NewPackedFrequencyTable(17)
		assert.Error(t, err)
	})
}

// TestNewPackedArray tests packed array creation.
func TestNewPackedArray(t *testing.T) {
	t.Run("valid", func(t *testing.T) {
		pa, err := NewPackedArray(64, 32)
		require.NoError(t, err)
		assert.NotNil(t, pa)
	})

	t.Run("invalid_itemCount", func(t *testing.T) {
		_, err := NewPackedArray(0, 8)
		assert.Error(t, err)
	})

	t.Run("invalid_bitsPerItem_zero", func(t *testing.T) {
		_, err := NewPackedArray(64, 0)
		assert.Error(t, err)
	})

	t.Run("invalid_bitsPerItem_too_large", func(t *testing.T) {
		_, err := NewPackedArray(64, 33)
		assert.Error(t, err)
	})
}

// TestPackedArray_GetSet tests basic get and set operations.
func TestPackedArray_GetSet(t *testing.T) {
	pa, err := NewPackedArray(16, 8)
	require.NoError(t, err)

	t.Run("set_and_get_zero", func(t *testing.T) {
		assert.Equal(t, uint32(0), pa.Get(0))
	})

	t.Run("set_and_get_value", func(t *testing.T) {
		require.NoError(t, pa.Set(0, 42))
		assert.Equal(t, uint32(42), pa.Get(0))
	})

	t.Run("set_multiple_indices", func(t *testing.T) {
		require.NoError(t, pa.Set(3, 10))
		require.NoError(t, pa.Set(7, 200))
		assert.Equal(t, uint32(10), pa.Get(3))
		assert.Equal(t, uint32(200), pa.Get(7))
	})

	t.Run("out_of_range_get", func(t *testing.T) {
		assert.Equal(t, uint32(0), pa.Get(-1))
		assert.Equal(t, uint32(0), pa.Get(1000))
	})

	t.Run("out_of_range_set", func(t *testing.T) {
		err := pa.Set(-1, 5)
		assert.Error(t, err)
		err = pa.Set(1000, 5)
		assert.Error(t, err)
	})
}

// TestPackedArray_Add tests the increment operation.
func TestPackedArray_Add(t *testing.T) {
	pa, err := NewPackedArray(16, 16)
	require.NoError(t, err)

	require.NoError(t, pa.Add(0, 10))
	assert.Equal(t, uint32(10), pa.Get(0))

	require.NoError(t, pa.Add(0, 5))
	assert.Equal(t, uint32(15), pa.Get(0))

	require.NoError(t, pa.Add(1, 100))
	assert.Equal(t, uint32(100), pa.Get(1))
}

// TestPackedArray_ExpandBits tests automatic bit expansion.
func TestPackedArray_ExpandBits(t *testing.T) {
	pa, err := NewPackedArray(8, 8)
	require.NoError(t, err)

	// Set a small value first.
	require.NoError(t, pa.Set(0, 1))

	// Now set a large value that forces bit expansion.
	require.NoError(t, pa.Set(1, 1000))
	assert.Equal(t, uint32(1000), pa.Get(1))
	// Previous value should be preserved.
	assert.Equal(t, uint32(1), pa.Get(0))
}

// TestPackedArray_GetMemoryUsage tests memory usage reporting.
func TestPackedArray_GetMemoryUsage(t *testing.T) {
	pa, err := NewPackedArray(256, 8)
	require.NoError(t, err)

	usage := pa.GetMemoryUsage()
	assert.Greater(t, usage, int64(0))
}

// TestKmerStringToIndex tests k-mer string to index conversion.
func TestKmerStringToIndex(t *testing.T) {
	tests := []struct {
		kmer     string
		expected int
		hasError bool
	}{
		{"A", 0, false},
		{"T", 1, false},
		{"G", 2, false},
		{"C", 3, false},
		{"AA", 0, false},
		{"AT", 1, false},
		{"AG", 2, false},
		{"AC", 3, false},
		{"TA", 4, false},
		{"TT", 5, false},
		{"ATGC", 0*64 + 1*16 + 2*4 + 3, false}, // A=0, T=1, G=2, C=3
		{"N", 0, true},   // invalid base
		{"X", 0, true},   // invalid base
		{"ATX", 0, true}, // invalid base in middle
	}

	for _, tt := range tests {
		t.Run(tt.kmer, func(t *testing.T) {
			idx, err := KmerStringToIndex(tt.kmer)
			if tt.hasError {
				assert.Error(t, err)
			} else {
				require.NoError(t, err)
				assert.Equal(t, tt.expected, idx)
			}
		})
	}
}

// TestKmerStringToIndex_LowerCase tests lowercase base handling.
func TestKmerStringToIndex_LowerCase(t *testing.T) {
	// Lowercase should map to the same index as uppercase.
	idxUpper, err := KmerStringToIndex("ATGC")
	require.NoError(t, err)

	idxLower, err := KmerStringToIndex("atgc")
	require.NoError(t, err)

	assert.Equal(t, idxUpper, idxLower)
}

// TestIndexToKmerString tests index to k-mer string conversion.
func TestIndexToKmerString(t *testing.T) {
	// Round-trip tests.
	tests := []string{"A", "T", "G", "C", "AT", "GC", "ATGC", "AAAA"}
	for _, kmer := range tests {
		t.Run(kmer, func(t *testing.T) {
			idx, err := KmerStringToIndex(kmer)
			require.NoError(t, err)
			result := IndexToKmerString(idx, len(kmer))
			assert.Equal(t, kmer, result)
		})
	}
}

// TestPackedFrequencyTable_IncrementKmer tests k-mer frequency incrementing.
func TestPackedFrequencyTable_IncrementKmer(t *testing.T) {
	pft, err := NewPackedFrequencyTable(3)
	require.NoError(t, err)

	require.NoError(t, pft.IncrementKmer("ATG"))
	require.NoError(t, pft.IncrementKmer("ATG"))
	require.NoError(t, pft.IncrementKmer("GCA"))

	assert.Equal(t, 3, pft.TotalCount)
	assert.Equal(t, 2, pft.NonZeroCount)

	freq, err := pft.GetFrequency("ATG")
	require.NoError(t, err)
	assert.Equal(t, 2, freq)

	freq, err = pft.GetFrequency("GCA")
	require.NoError(t, err)
	assert.Equal(t, 1, freq)
}

// TestPackedFrequencyTable_IncrementKmer_WrongLength tests wrong k-mer length error.
func TestPackedFrequencyTable_IncrementKmer_WrongLength(t *testing.T) {
	pft, err := NewPackedFrequencyTable(3)
	require.NoError(t, err)

	err = pft.IncrementKmer("AT") // Wrong length
	assert.Error(t, err)

	err = pft.IncrementKmer("ATGC") // Wrong length
	assert.Error(t, err)
}

// TestPackedFrequencyTable_GetFrequency_WrongLength tests wrong length in GetFrequency.
func TestPackedFrequencyTable_GetFrequency_WrongLength(t *testing.T) {
	pft, err := NewPackedFrequencyTable(3)
	require.NoError(t, err)

	_, err = pft.GetFrequency("AT")
	assert.Error(t, err)
}

// TestPackedFrequencyTable_GetFrequency_InvalidBase tests invalid base in GetFrequency.
func TestPackedFrequencyTable_GetFrequency_InvalidBase(t *testing.T) {
	pft, err := NewPackedFrequencyTable(3)
	require.NoError(t, err)

	_, err = pft.GetFrequency("ATX")
	assert.Error(t, err)
}

// TestPackedFrequencyTable_ToFrequencyTable tests conversion to regular FrequencyTable.
func TestPackedFrequencyTable_ToFrequencyTable(t *testing.T) {
	pft, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)

	kmers := []string{"AT", "GC", "AT", "TA", "GC", "GC"}
	for _, kmer := range kmers {
		require.NoError(t, pft.IncrementKmer(kmer))
	}

	ft := pft.ToFrequencyTable()
	require.NotNil(t, ft)

	assert.Equal(t, 2, ft.K)
	assert.Equal(t, 2, ft.Frequencies["AT"])
	assert.Equal(t, 3, ft.Frequencies["GC"])
	assert.Equal(t, 1, ft.Frequencies["TA"])
	assert.Equal(t, len(kmers), ft.TotalCount)
}

// TestFromFrequencyTable tests creating a PackedFrequencyTable from a FrequencyTable.
func TestFromFrequencyTable(t *testing.T) {
	ft := &FrequencyTable{
		K: 2,
		Frequencies: map[string]int{
			"AT": 5,
			"GC": 3,
			"TA": 1,
		},
		TotalCount: 9,
		TotalKmers: 10,
	}

	pft, err := FromFrequencyTable(ft)
	require.NoError(t, err)
	require.NotNil(t, pft)

	assert.Equal(t, 2, pft.K)
	assert.Equal(t, 9, pft.TotalCount)
	assert.Equal(t, 10, pft.TotalKmers)

	freq, err := pft.GetFrequency("AT")
	require.NoError(t, err)
	assert.Equal(t, 5, freq)

	freq, err = pft.GetFrequency("GC")
	require.NoError(t, err)
	assert.Equal(t, 3, freq)
}

// TestPackedFrequencyTable_GetTopKmers tests top k-mer retrieval.
func TestPackedFrequencyTable_GetTopKmers(t *testing.T) {
	pft, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)

	// Add k-mers with different frequencies.
	for i := 0; i < 10; i++ {
		require.NoError(t, pft.IncrementKmer("AT"))
	}
	for i := 0; i < 5; i++ {
		require.NoError(t, pft.IncrementKmer("GC"))
	}
	require.NoError(t, pft.IncrementKmer("TA"))

	top2 := pft.GetTopKmers(2)
	require.Len(t, top2, 2)
	assert.Equal(t, "AT", top2[0].Kmer)
	assert.Equal(t, 10, top2[0].Count)
	assert.Equal(t, "GC", top2[1].Kmer)

	// Request more than available.
	all := pft.GetTopKmers(100)
	assert.Equal(t, 3, len(all))
}

// TestPackedFrequencyTable_GetMemoryUsage tests memory usage reporting.
func TestPackedFrequencyTable_GetMemoryUsage(t *testing.T) {
	pft, err := NewPackedFrequencyTable(4)
	require.NoError(t, err)

	usage := pft.GetMemoryUsage()
	assert.Greater(t, usage, int64(0))
}

// TestPackedFrequencyTable_GetCompressionRatio tests compression ratio reporting.
func TestPackedFrequencyTable_GetCompressionRatio(t *testing.T) {
	pft, err := NewPackedFrequencyTable(3)
	require.NoError(t, err)

	// Add some data so compression ratio is meaningful.
	for _, kmer := range []string{"ATG", "GCA", "TGC", "CAT"} {
		require.NoError(t, pft.IncrementKmer(kmer))
	}

	ratio := pft.GetCompressionRatio()
	assert.Greater(t, ratio, 0.0)
}

// TestPackedFrequencyTable_GetStatistics tests the statistics summary.
func TestPackedFrequencyTable_GetStatistics(t *testing.T) {
	pft, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)

	for _, kmer := range []string{"AT", "GC", "AT", "TA"} {
		require.NoError(t, pft.IncrementKmer(kmer))
	}

	stats := pft.GetStatistics()
	assert.Equal(t, 2, stats.K)
	assert.Equal(t, 16, stats.MaxKmers) // 4^2
	assert.Equal(t, 3, stats.NonZeroKmers)
	assert.Equal(t, 4, stats.TotalCount)
	assert.Greater(t, stats.MemoryUsage, int64(0))
	assert.GreaterOrEqual(t, stats.Sparsity, 0.0)
	assert.LessOrEqual(t, stats.Sparsity, 1.0)
}

// TestCountLeadingZeros tests the internal leading-zeros counter.
func TestCountLeadingZeros(t *testing.T) {
	tests := []struct {
		x        uint32
		expected int
	}{
		{0, 32},
		{1, 31},
		{2, 30},
		{4, 29},
		{0xFFFF, 16},
		{0xFFFFFFFF, 0},
		{0x80000000, 0},
	}

	for _, tt := range tests {
		t.Run("", func(t *testing.T) {
			result := countLeadingZeros(tt.x)
			assert.Equal(t, tt.expected, result)
		})
	}
}

// TestPackedArray_RoundTrip tests that values survive a set/get round trip for all indices.
func TestPackedArray_RoundTrip(t *testing.T) {
	pft, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)

	// Populate all 16 possible 2-mers.
	kmers2 := GetExpectedKmers(2)
	for i, kmer := range kmers2 {
		for j := 0; j <= i; j++ {
			require.NoError(t, pft.IncrementKmer(kmer))
		}
	}

	// Verify all counts.
	for i, kmer := range kmers2 {
		freq, err := pft.GetFrequency(kmer)
		require.NoError(t, err)
		assert.Equal(t, i+1, freq, "kmer %s should have count %d", kmer, i+1)
	}
}

// TestCalculatePackedArraySize tests the size calculation helper.
func TestCalculatePackedArraySize(t *testing.T) {
	tests := []struct {
		itemCount   int
		bitsPerItem int
		expected    int
	}{
		{32, 1, 1},   // 32 bits = 1 uint32
		{33, 1, 2},   // 33 bits needs 2 uint32s
		{16, 2, 1},   // 32 bits = 1 uint32
		{17, 2, 2},   // 34 bits needs 2 uint32s
		{8, 4, 1},    // 32 bits = 1 uint32
		{4, 8, 1},    // 32 bits = 1 uint32
		{5, 8, 2},    // 40 bits needs 2 uint32s
	}

	for _, tt := range tests {
		result := calculatePackedArraySize(tt.itemCount, tt.bitsPerItem)
		assert.Equal(t, tt.expected, result, "itemCount=%d bitsPerItem=%d", tt.itemCount, tt.bitsPerItem)
	}
}
