package kmer

import (
	"bytes"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewPackedFrequencyExtractor tests the constructor defaults.
func TestNewPackedFrequencyExtractor(t *testing.T) {
	ext := NewPackedFrequencyExtractor(4)
	require.NotNil(t, ext)
	assert.Equal(t, 4, ext.K)
	assert.True(t, ext.IgnoreN)
	assert.False(t, ext.CaseSensitive)
	assert.False(t, ext.CircularSequence)
	assert.Equal(t, 1, ext.MinFrequency)
	assert.True(t, ext.UsePackedStorage)
}

// TestNewPackedFrequencyExtractorWithConfig tests config-based construction.
func TestNewPackedFrequencyExtractorWithConfig(t *testing.T) {
	cfg := PackedFrequencyExtractorConfig{
		K:              3,
		IgnoreN:        false,
		CaseSensitive:  true,
		MinFrequency:   2,
		UsePackedStorage: true,
	}
	ext := NewPackedFrequencyExtractorWithConfig(cfg)
	require.NotNil(t, ext)
	assert.Equal(t, 3, ext.K)
	assert.False(t, ext.IgnoreN)
	assert.True(t, ext.CaseSensitive)
	assert.Equal(t, 2, ext.MinFrequency)
}

// TestNewPackedFrequencyExtractorWithConfig_DefaultMinFreq tests that zero MinFrequency defaults to 1.
func TestNewPackedFrequencyExtractorWithConfig_DefaultMinFreq(t *testing.T) {
	cfg := PackedFrequencyExtractorConfig{
		K:            3,
		MinFrequency: 0, // should default to 1
	}
	ext := NewPackedFrequencyExtractorWithConfig(cfg)
	assert.Equal(t, 1, ext.MinFrequency)
}

// TestPackedFrequencyExtractor_ExtractFromSequence tests basic extraction.
func TestPackedFrequencyExtractor_ExtractFromSequence(t *testing.T) {
	ext := NewPackedFrequencyExtractor(3)

	t.Run("basic_extraction", func(t *testing.T) {
		seq := []byte("ATCGATCG") // 6 3-mers: ATC, TCG, CGA, GAT, ATC, TCG
		pft, err := ext.ExtractFromSequence(seq)
		require.NoError(t, err)
		require.NotNil(t, pft)

		freq, err := pft.GetFrequency("ATC")
		require.NoError(t, err)
		assert.Equal(t, 2, freq)

		freq, err = pft.GetFrequency("TCG")
		require.NoError(t, err)
		assert.Equal(t, 2, freq)

		freq, err = pft.GetFrequency("CGA")
		require.NoError(t, err)
		assert.Equal(t, 1, freq)
	})

	t.Run("sequence_too_short", func(t *testing.T) {
		_, err := ext.ExtractFromSequence([]byte("AT"))
		assert.Error(t, err)
	})
}

// TestPackedFrequencyExtractor_ExtractFromSequence_CaseInsensitive tests case folding.
func TestPackedFrequencyExtractor_ExtractFromSequence_CaseInsensitive(t *testing.T) {
	ext := NewPackedFrequencyExtractor(3)
	ext.CaseSensitive = false

	seq := []byte("atcgatcg")
	pft, err := ext.ExtractFromSequence(seq)
	require.NoError(t, err)

	freq, err := pft.GetFrequency("ATC")
	require.NoError(t, err)
	assert.Equal(t, 2, freq)
}

// TestPackedFrequencyExtractor_ExtractFromSequence_IgnoreN tests N-base skipping.
func TestPackedFrequencyExtractor_ExtractFromSequence_IgnoreN(t *testing.T) {
	ext := NewPackedFrequencyExtractor(3)
	ext.IgnoreN = true

	seq := []byte("ATCNGATCG")
	pft, err := ext.ExtractFromSequence(seq)
	require.NoError(t, err)

	// None of ATC should appear 0 times because they contain N.
	_, err = pft.GetFrequency("TCN")
	assert.Error(t, err) // "N" not a valid base for GetFrequency

	// "ATCNGATCG" with k=3 and N-skip gives: ATC, GAT, ATC, TCG
	// So ATC appears twice.
	freq, err := pft.GetFrequency("ATC")
	require.NoError(t, err)
	assert.Equal(t, 2, freq)
}

// TestPackedFrequencyExtractor_ExtractFromSequence_Circular tests circular sequence.
func TestPackedFrequencyExtractor_ExtractFromSequence_Circular(t *testing.T) {
	ext := NewPackedFrequencyExtractor(3)
	ext.CircularSequence = true

	seq := []byte("ATCG")
	pft, err := ext.ExtractFromSequence(seq)
	require.NoError(t, err)

	// For circular seq of length 4, should have 4 k-mers.
	assert.Equal(t, 4, pft.TotalKmers)
}

// TestPackedFrequencyExtractor_ExtractFromSequence_MinFrequency tests frequency filtering.
func TestPackedFrequencyExtractor_ExtractFromSequence_MinFrequency(t *testing.T) {
	ext := NewPackedFrequencyExtractor(3)
	ext.MinFrequency = 2

	seq := []byte("ATCGATCG") // ATC and TCG appear twice, others once
	pft, err := ext.ExtractFromSequence(seq)
	require.NoError(t, err)

	// Only k-mers with count >= 2 should be retained.
	assert.Equal(t, 2, pft.NonZeroCount)

	freq, err := pft.GetFrequency("ATC")
	require.NoError(t, err)
	assert.Equal(t, 2, freq)

	// CGA appeared only once, should be filtered out.
	freq, err = pft.GetFrequency("CGA")
	require.NoError(t, err)
	assert.Equal(t, 0, freq)
}

// TestPackedFrequencyExtractor_ExtractFromMultipleSequences tests multi-sequence extraction.
func TestPackedFrequencyExtractor_ExtractFromMultipleSequences(t *testing.T) {
	ext := NewPackedFrequencyExtractor(3)

	t.Run("multiple_sequences", func(t *testing.T) {
		seqs := [][]byte{
			[]byte("ATCG"),
			[]byte("GATC"),
		}

		pft, err := ext.ExtractFromMultipleSequences(seqs)
		require.NoError(t, err)
		require.NotNil(t, pft)

		// ATC appears in both sequences.
		freq, err := pft.GetFrequency("ATC")
		require.NoError(t, err)
		assert.Equal(t, 2, freq)
	})

	t.Run("empty_sequences", func(t *testing.T) {
		_, err := ext.ExtractFromMultipleSequences([][]byte{})
		assert.Error(t, err)
	})

	t.Run("one_sequence_too_short", func(t *testing.T) {
		seqs := [][]byte{
			[]byte("ATCG"),
			[]byte("AT"), // Too short
		}
		_, err := ext.ExtractFromMultipleSequences(seqs)
		assert.Error(t, err)
	})
}

// TestPackedFrequencyExtractor_ExtractFromSequenceBatch tests batched extraction.
func TestPackedFrequencyExtractor_ExtractFromSequenceBatch(t *testing.T) {
	ext := NewPackedFrequencyExtractor(3)

	t.Run("batch_processing", func(t *testing.T) {
		seqs := make([][]byte, 10)
		for i := range seqs {
			seqs[i] = []byte("ATCGATCG")
		}

		pft, err := ext.ExtractFromSequenceBatch(seqs, 3)
		require.NoError(t, err)
		require.NotNil(t, pft)

		// ATC should appear 2*10=20 times.
		freq, err := pft.GetFrequency("ATC")
		require.NoError(t, err)
		assert.Equal(t, 20, freq)
	})

	t.Run("empty_sequences", func(t *testing.T) {
		_, err := ext.ExtractFromSequenceBatch([][]byte{}, 10)
		assert.Error(t, err)
	})

	t.Run("default_batch_size", func(t *testing.T) {
		seqs := [][]byte{[]byte("ATCGATCG")}
		pft, err := ext.ExtractFromSequenceBatch(seqs, 0) // 0 triggers default
		require.NoError(t, err)
		require.NotNil(t, pft)
	})
}

// TestPackedFrequencyExtractor_ExtractFromSequenceStreaming tests streaming extraction.
func TestPackedFrequencyExtractor_ExtractFromSequenceStreaming(t *testing.T) {
	ext := NewPackedFrequencyExtractor(3)

	seq := []byte("ATCGATCGATCGATCG")
	reader := bytes.NewReader(seq)

	pft, err := ext.ExtractFromSequenceStreaming(reader, 8)
	require.NoError(t, err)
	require.NotNil(t, pft)

	// Should extract k-mers from the streaming data.
	freq, err := pft.GetFrequency("ATC")
	require.NoError(t, err)
	assert.Greater(t, freq, 0)
}

// TestPackedFrequencyExtractor_GetMemoryEstimate tests memory estimate calculation.
func TestPackedFrequencyExtractor_GetMemoryEstimate(t *testing.T) {
	ext := NewPackedFrequencyExtractor(4)

	estimate := ext.GetMemoryEstimate(100, 1000)
	assert.Greater(t, estimate.MaxKmers, 0)
	assert.Greater(t, estimate.PackedMemoryUsage, int64(0))
	assert.GreaterOrEqual(t, estimate.EstimatedUniqueKmers, 0)
}

// TestPackedFrequencyExtractor_GetOptimalBatchSize tests batch size calculation.
func TestPackedFrequencyExtractor_GetOptimalBatchSize(t *testing.T) {
	ext := NewPackedFrequencyExtractor(4)

	batchSize := ext.GetOptimalBatchSize(100, 1000)
	assert.GreaterOrEqual(t, batchSize, 100)
	assert.LessOrEqual(t, batchSize, 10000)
}

// TestPackedFrequencyExtractor_IsValidKmer tests the validity check.
func TestPackedFrequencyExtractor_IsValidKmer(t *testing.T) {
	ext := NewPackedFrequencyExtractor(3)

	assert.True(t, ext.isValidKmer("ATC"))
	assert.True(t, ext.isValidKmer("GGG"))
	assert.False(t, ext.isValidKmer("ATX")) // invalid base
	assert.True(t, ext.isValidKmer(""))    // empty is valid

	// N handling.
	ext.IgnoreN = true
	assert.False(t, ext.isValidKmer("ATN"))

	ext.IgnoreN = false
	assert.True(t, ext.isValidKmer("ATN"))
}

// TestPackedFrequencyExtractor_LargeSequence tests extraction on a larger sequence.
func TestPackedFrequencyExtractor_LargeSequence(t *testing.T) {
	ext := NewPackedFrequencyExtractor(4)

	// Generate a 10KB deterministic sequence.
	seq := make([]byte, 10000)
	bases := []byte{'A', 'T', 'G', 'C'}
	for i := range seq {
		seq[i] = bases[i%4]
	}

	pft, err := ext.ExtractFromSequence(seq)
	require.NoError(t, err)
	require.NotNil(t, pft)

	assert.Greater(t, pft.TotalCount, 0)
	assert.Greater(t, pft.NonZeroCount, 0)
}

// TestPackedFrequencyExtractor_MergePackedTables tests the internal merge helper.
func TestPackedFrequencyExtractor_MergePackedTables(t *testing.T) {
	ext := NewPackedFrequencyExtractor(2)

	dst, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)
	require.NoError(t, dst.IncrementKmer("AT"))
	require.NoError(t, dst.IncrementKmer("GC"))

	src, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)
	require.NoError(t, src.IncrementKmer("AT"))
	require.NoError(t, src.IncrementKmer("TA"))

	err = ext.mergePackedTables(dst, src)
	require.NoError(t, err)

	freq, err := dst.GetFrequency("AT")
	require.NoError(t, err)
	assert.Equal(t, 2, freq) // 1 + 1

	freq, err = dst.GetFrequency("GC")
	require.NoError(t, err)
	assert.Equal(t, 1, freq) // unchanged

	freq, err = dst.GetFrequency("TA")
	require.NoError(t, err)
	assert.Equal(t, 1, freq) // from src
}

// TestPackedFrequencyExtractor_MergePackedTables_MismatchedK tests k mismatch error.
func TestPackedFrequencyExtractor_MergePackedTables_MismatchedK(t *testing.T) {
	ext := NewPackedFrequencyExtractor(2)

	dst, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)

	src, err := NewPackedFrequencyTable(3)
	require.NoError(t, err)

	err = ext.mergePackedTables(dst, src)
	assert.Error(t, err)
}

// BenchmarkPackedFrequencyExtractor benchmarks extraction performance.
func BenchmarkPackedFrequencyExtractor(b *testing.B) {
	ext := NewPackedFrequencyExtractor(6)

	seq := make([]byte, 100000)
	bases := []byte{'A', 'T', 'G', 'C'}
	for i := range seq {
		seq[i] = bases[i%4]
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := ext.ExtractFromSequence(seq)
		if err != nil {
			b.Fatal(err)
		}
	}
}
