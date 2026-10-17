//go:build !release
// +build !release

package main

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestGenerateTestSequences verifies that generateTestSequences produces the
// requested number of sequences, each with the correct length and only valid
// DNA bases.
func TestGenerateTestSequences(t *testing.T) {
	count, length := 8, 50
	seqs := generateTestSequences(count, length)

	require.Len(t, seqs, count, "must return exactly count sequences")

	validBases := map[byte]bool{'A': true, 'T': true, 'G': true, 'C': true}
	for i, seq := range seqs {
		assert.Len(t, seq, length, "sequence %d must have length %d", i, length)
		for j, b := range seq {
			assert.True(t, validBases[b],
				"sequence %d, position %d: unexpected base %q", i, j, b)
		}
	}
}

// TestGenerateTestSequences_Zero verifies edge-case: zero count returns empty slice.
func TestGenerateTestSequences_Zero(t *testing.T) {
	seqs := generateTestSequences(0, 100)
	assert.Empty(t, seqs, "generateTestSequences(0, 100) must return empty slice")
}

// TestGenerateTestSequences_Deterministic verifies that the sequence generation
// is deterministic (same inputs always produce identical output).
func TestGenerateTestSequences_Deterministic(t *testing.T) {
	seqs1 := generateTestSequences(5, 30)
	seqs2 := generateTestSequences(5, 30)
	require.Len(t, seqs1, len(seqs2))
	for i := range seqs1 {
		assert.Equal(t, seqs1[i], seqs2[i], "sequence %d must be identical across calls", i)
	}
}

// TestRunCPUBenchmark_DoesNotPanic verifies that runCPUBenchmark completes without
// panicking on valid input. CPU k-mer extraction is always available.
func TestRunCPUBenchmark_DoesNotPanic(t *testing.T) {
	seqs := generateTestSequences(10, 100)
	assert.NotPanics(t, func() {
		runCPUBenchmark(seqs, 4, false)
	}, "runCPUBenchmark must not panic")
}

// TestRunCPUBenchmark_Verbose verifies that verbose mode does not panic.
func TestRunCPUBenchmark_Verbose(t *testing.T) {
	seqs := generateTestSequences(5, 80)
	assert.NotPanics(t, func() {
		runCPUBenchmark(seqs, 4, true)
	}, "runCPUBenchmark(verbose=true) must not panic")
}

// TestRunCPUBenchmark_EmptySequences verifies runCPUBenchmark handles an empty
// input slice gracefully (error path inside the function).
func TestRunCPUBenchmark_EmptySequences(t *testing.T) {
	assert.NotPanics(t, func() {
		runCPUBenchmark([][]byte{}, 4, false)
	}, "runCPUBenchmark with empty sequences must not panic")
}

// TestRunGPUBenchmark_DoesNotPanic verifies that runGPUBenchmark handles the case
// where Metal is unavailable (shader library not compiled) without panicking.
func TestRunGPUBenchmark_DoesNotPanic(t *testing.T) {
	seqs := generateTestSequences(5, 80)
	assert.NotPanics(t, func() {
		runGPUBenchmark(seqs, 4, false)
	}, "runGPUBenchmark must not panic when Metal is unavailable")
}

// TestRunGPUBenchmark_Verbose verifies that verbose mode does not panic even
// when Metal initialisation fails (prints error and returns early).
func TestRunGPUBenchmark_Verbose(t *testing.T) {
	seqs := generateTestSequences(5, 80)
	assert.NotPanics(t, func() {
		runGPUBenchmark(seqs, 4, true)
	}, "runGPUBenchmark(verbose=true) must not panic")
}

// TestRunPerformanceComparison_DoesNotPanic verifies that the comparison function
// does not panic when Metal is unavailable in the test environment.
func TestRunPerformanceComparison_DoesNotPanic(t *testing.T) {
	seqs := generateTestSequences(5, 80)
	assert.NotPanics(t, func() {
		runPerformanceComparison(seqs, 4, false)
	}, "runPerformanceComparison must not panic when Metal is unavailable")
}
