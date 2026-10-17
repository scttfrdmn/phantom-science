//go:build !release
// +build !release

package main

import (
	"testing"
)

func TestGenerateTestSequences_Count(t *testing.T) {
	seqs := generateTestSequences(10, 100)
	if len(seqs) != 10 {
		t.Errorf("generateTestSequences(10, 100) returned %d sequences, want 10", len(seqs))
	}
}

func TestGenerateTestSequences_Length(t *testing.T) {
	seqs := generateTestSequences(5, 50)
	for i, seq := range seqs {
		if len(seq) != 50 {
			t.Errorf("sequence[%d] has length %d, want 50", i, len(seq))
		}
	}
}

func TestGenerateTestSequences_ValidBases(t *testing.T) {
	seqs := generateTestSequences(3, 100)
	validBases := map[byte]bool{'A': true, 'T': true, 'G': true, 'C': true}
	for i, seq := range seqs {
		for j, b := range seq {
			if !validBases[b] {
				t.Errorf("sequence[%d][%d] = %c, not a valid base", i, j, b)
			}
		}
	}
}

func TestGenerateTestSequences_Zero(t *testing.T) {
	seqs := generateTestSequences(0, 100)
	if len(seqs) != 0 {
		t.Errorf("generateTestSequences(0, 100) returned %d sequences, want 0", len(seqs))
	}
}

// Tests for benchmark functions that write to stdout but are testable

func TestRunMapBenchmark_NoError(t *testing.T) {
	seqs := generateTestSequences(5, 200)
	// Should complete without panic
	runMapBenchmark(seqs, 4, false)
}

func TestRunMapBenchmark_Verbose(t *testing.T) {
	seqs := generateTestSequences(3, 150)
	runMapBenchmark(seqs, 4, true)
}

func TestRunPackedBenchmark_NoError(t *testing.T) {
	seqs := generateTestSequences(5, 200)
	runPackedBenchmark(seqs, 4, false)
}

func TestRunPackedBenchmark_Verbose(t *testing.T) {
	seqs := generateTestSequences(3, 150)
	runPackedBenchmark(seqs, 4, true)
}

func TestRunMemoryBenchmark_NoError(t *testing.T) {
	seqs := generateTestSequences(5, 200)
	runMemoryBenchmark(seqs, 4, false)
}

func TestRunPerformanceComparison_NoError(t *testing.T) {
	seqs := generateTestSequences(5, 200)
	runPerformanceComparison(seqs, 4, false)
}

func TestRunAccuracyComparison_NoError(t *testing.T) {
	seqs := generateTestSequences(5, 200)
	runAccuracyComparison(seqs, 4, false)
}

func TestRunAccuracyComparison_Verbose(t *testing.T) {
	seqs := generateTestSequences(5, 200)
	runAccuracyComparison(seqs, 4, true)
}
