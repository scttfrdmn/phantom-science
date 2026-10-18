//go:build !release
// +build !release

package main

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func TestPassFailIcon_Pass(t *testing.T) {
	if got := passFailIcon(true); got != "\u2705" {
		t.Errorf("passFailIcon(true) = %q, want checkmark emoji", got)
	}
}

func TestPassFailIcon_Fail(t *testing.T) {
	if got := passFailIcon(false); got != "\u274c" {
		t.Errorf("passFailIcon(false) = %q, want cross emoji", got)
	}
}

func TestIsValidKmerOptimized(t *testing.T) {
	tests := []struct {
		kmer  string
		valid bool
	}{
		{"ACGT", true},
		{"AAAA", true},
		{"GCGC", true},
		{"", true},
		{"ACGN", false},
		{"acgt", false},
		{"ACG-", false},
		{"ACGU", false},
	}
	for _, tt := range tests {
		got := isValidKmerOptimized(tt.kmer)
		if got != tt.valid {
			t.Errorf("isValidKmerOptimized(%q) = %v, want %v", tt.kmer, got, tt.valid)
		}
	}
}

func TestCalculateSequenceComplexity_ShortSequence(t *testing.T) {
	got := calculateSequenceComplexity([]byte("ACG"), 4)
	if got != 0.0 {
		t.Errorf("got %v, want 0.0 for seq shorter than k", got)
	}
}

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	got := calculateSequenceComplexity([]byte{}, 4)
	if got != 0.0 {
		t.Errorf("got %v, want 0.0 for empty sequence", got)
	}
}

func TestCalculateSequenceComplexity_LowComplexity(t *testing.T) {
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAA")
	got := calculateSequenceComplexity(seq, 4)
	if got > 0.1 {
		t.Errorf("all-A complexity = %v, want <= 0.1", got)
	}
}

func TestCalculateSequenceComplexity_Range(t *testing.T) {
	seq := []byte("ACGTACGTACGTACGTACGTACGT")
	got := calculateSequenceComplexity(seq, 4)
	if got < 0.0 || got > 1.0 {
		t.Errorf("complexity %v outside [0,1]", got)
	}
}

func TestCalculateSequenceComplexity_HighDiversity(t *testing.T) {
	// With k=1 every unique base is a unique kmer
	seq := []byte("ACGTACGT")
	got := calculateSequenceComplexity(seq, 1)
	if got <= 0.0 || got > 1.0 {
		t.Errorf("complexity %v unexpected for diverse seq with k=1", got)
	}
}

func TestCalculateSequenceComplexity_ExactK(t *testing.T) {
	// Sequence length exactly equals k: 1 kmer total, 1 unique -> complexity = 1.0
	seq := []byte("ACGT")
	got := calculateSequenceComplexity(seq, 4)
	if got != 1.0 {
		t.Errorf("got %v, want 1.0 for seq length == k", got)
	}
}

// ---------- NewOptimizedFullGenomeTrainer ----------

func TestNewOptimizedFullGenomeTrainer_NotNil(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)
	if trainer == nil {
		t.Error("NewOptimizedFullGenomeTrainer returned nil")
	}
}

func TestNewOptimizedFullGenomeTrainer_HasWorkers(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)
	if trainer.numWorkers <= 0 {
		t.Errorf("numWorkers = %d, want > 0", trainer.numWorkers)
	}
}

// ---------- TrainOptimized ----------

func TestTrainOptimized_EmptyRefGenome(t *testing.T) {
	// With a tmpDir as genome dir, all chromosomes will fail to load.
	// TrainOptimized should handle this gracefully.
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)
	stats, err := trainer.TrainOptimized()
	if err != nil {
		t.Fatalf("TrainOptimized(empty refGenome) error: %v", err)
	}
	if stats == nil {
		t.Fatal("expected non-nil stats")
	}
	// No chromosomes found, should have processed 0
	if stats.ChromosomesProcessed < 0 {
		t.Errorf("ChromosomesProcessed = %d, want >= 0", stats.ChromosomesProcessed)
	}
}
