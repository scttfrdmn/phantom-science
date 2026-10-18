//go:build !release
// +build !release

package main

import (
	"math"
	"testing"
)

func TestFindMaxHomopolymerRun_Empty(t *testing.T) {
	if got := findMaxHomopolymerRun([]byte{}); got != 0 {
		t.Errorf("findMaxHomopolymerRun([]) = %d, want 0", got)
	}
}

func TestFindMaxHomopolymerRun_Single(t *testing.T) {
	if got := findMaxHomopolymerRun([]byte("A")); got != 1 {
		t.Errorf("findMaxHomopolymerRun(A) = %d, want 1", got)
	}
}

func TestFindMaxHomopolymerRun_WithRun(t *testing.T) {
	// ATAAAACG has a run of 4 A's
	if got := findMaxHomopolymerRun([]byte("ATAAAACG")); got != 4 {
		t.Errorf("findMaxHomopolymerRun(ATAAAACG) = %d, want 4", got)
	}
}

func TestFindMaxHomopolymerRun_AllSame(t *testing.T) {
	if got := findMaxHomopolymerRun([]byte("AAAA")); got != 4 {
		t.Errorf("findMaxHomopolymerRun(AAAA) = %d, want 4", got)
	}
}

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	if got := calculateSequenceComplexity([]byte{}); got != 0.0 {
		t.Errorf("calculateSequenceComplexity([]) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_TooShort(t *testing.T) {
	// length < 3 returns 0
	if got := calculateSequenceComplexity([]byte("AT")); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(AT) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_Valid(t *testing.T) {
	// ATCGATCG has multiple distinct 3-mers, complexity > 0
	got := calculateSequenceComplexity([]byte("ATCGATCG"))
	if got <= 0.0 || got > 1.0 {
		t.Errorf("calculateSequenceComplexity(ATCGATCG) = %f, want in (0,1]", got)
	}
}

func TestCalculateCpGDepletion_Empty(t *testing.T) {
	if got := calculateCpGDepletion([][]byte{}); got != 0.0 {
		t.Errorf("calculateCpGDepletion([]) = %f, want 0.0", got)
	}
}

func TestCalculateCpGDepletion_WithCpG(t *testing.T) {
	// ACGT has one CG dinucleotide
	seqs := [][]byte{[]byte("ACGT")}
	got := calculateCpGDepletion(seqs)
	if math.IsNaN(got) || got < 0.0 {
		t.Errorf("calculateCpGDepletion(ACGT) = %f, want non-negative", got)
	}
}

func TestCalculateJSDivergence_Identical(t *testing.T) {
	p := map[string]float64{"AAA": 0.5, "TTT": 0.5}
	q := map[string]float64{"AAA": 0.5, "TTT": 0.5}
	got := calculateJSDivergence(p, q)
	if math.Abs(got) > 1e-9 {
		t.Errorf("calculateJSDivergence(identical) = %f, want ~0", got)
	}
}

func TestCalculateBiologicalMetrics_Empty(t *testing.T) {
	metrics := calculateBiologicalMetrics([][]byte{})
	if metrics.GCContent != 0 && !math.IsNaN(metrics.GCContent) {
		t.Errorf("expected 0 GC for empty seqs, got %f", metrics.GCContent)
	}
}

func TestCalculateBiologicalMetrics_WithSeqs(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCATCGATCGATCG")}
	metrics := calculateBiologicalMetrics(seqs)
	if metrics.GCContent <= 0 || metrics.GCContent > 1 {
		t.Errorf("GCContent %f out of (0,1]", metrics.GCContent)
	}
	if metrics.ATContent <= 0 || metrics.ATContent > 1 {
		t.Errorf("ATContent %f out of (0,1]", metrics.ATContent)
	}
	if metrics.MaxHomopolymer <= 0 {
		t.Errorf("expected positive MaxHomopolymer, got %d", metrics.MaxHomopolymer)
	}
}

func TestCalculateBiologicalMetrics_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	metrics := calculateBiologicalMetrics(seqs)
	if metrics.GCContent != 1.0 {
		t.Errorf("GCContent = %f, want 1.0", metrics.GCContent)
	}
	if metrics.ATContent != 0.0 {
		t.Errorf("ATContent = %f, want 0.0", metrics.ATContent)
	}
}

func TestPrintMetrics_NoCrash(t *testing.T) {
	metrics := BiologicalMetrics{
		GCContent:          0.42,
		ATContent:          0.58,
		MaxHomopolymer:     5,
		SequenceComplexity: 0.7,
		CpGDepletion:       0.25,
	}
	printMetrics(metrics)
}

func TestCalculateQuickJSDivergence_EmptyRefKmers(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCG")}
	got := calculateQuickJSDivergence(seqs, []string{}, 4)
	if got < 0 {
		t.Errorf("JSD should be non-negative, got %f", got)
	}
}

func TestCalculateQuickJSDivergence_WithRefKmers(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCGATCG")}
	refKmers := []string{"ATCG", "TCGA", "CGAT", "GATC", "ATCG"}
	got := calculateQuickJSDivergence(seqs, refKmers, 4)
	if got < 0 || got > 1.0+1e-9 {
		t.Errorf("JSD = %f, should be in [0,1]", got)
	}
}

func TestCalculateQuickJSDivergence_EmptySeqs(t *testing.T) {
	refKmers := []string{"ATCG", "GCTA"}
	got := calculateQuickJSDivergence([][]byte{}, refKmers, 4)
	if got < 0 {
		t.Errorf("JSD should be non-negative for empty seqs, got %f", got)
	}
}
