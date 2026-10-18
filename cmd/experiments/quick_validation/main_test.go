//go:build !release
// +build !release

package main

import (
	"testing"
)

func TestQuickTestChromosome_GetName(t *testing.T) {
	c := &QuickTestChromosome{name: "chr1", sequence: []byte("ATCG")}
	if got := c.GetName(); got != "chr1" {
		t.Errorf("GetName() = %q, want \"chr1\"", got)
	}
}

func TestQuickTestChromosome_GetLength(t *testing.T) {
	c := &QuickTestChromosome{name: "chr1", sequence: []byte("ATCG")}
	if got := c.GetLength(); got != 4 {
		t.Errorf("GetLength() = %d, want 4", got)
	}
}

func TestQuickTestChromosome_GetSequence(t *testing.T) {
	seq := []byte("ATCG")
	c := &QuickTestChromosome{name: "chr1", sequence: seq}
	if string(c.GetSequence()) != "ATCG" {
		t.Errorf("GetSequence() = %q, want \"ATCG\"", c.GetSequence())
	}
}

func TestQuickTestGenome_GetChromosomeNames(t *testing.T) {
	g := &QuickTestGenome{
		chromosomes: map[string][]byte{
			"chr1": []byte("ATCG"),
			"chr2": []byte("GCTA"),
		},
	}
	names := g.GetChromosomeNames()
	if len(names) != 2 {
		t.Errorf("GetChromosomeNames() returned %d names, want 2", len(names))
	}
}

func TestQuickTestGenome_LoadChromosome_Found(t *testing.T) {
	g := &QuickTestGenome{
		chromosomes: map[string][]byte{
			"chr1": []byte("ATCG"),
		},
	}
	chr, err := g.LoadChromosome("chr1")
	if err != nil {
		t.Fatalf("LoadChromosome(chr1) error: %v", err)
	}
	if chr.GetName() != "chr1" {
		t.Errorf("expected name chr1, got %q", chr.GetName())
	}
}

func TestQuickTestGenome_LoadChromosome_NotFound(t *testing.T) {
	g := &QuickTestGenome{
		chromosomes: map[string][]byte{},
	}
	_, err := g.LoadChromosome("chrX")
	if err == nil {
		t.Error("expected error for missing chromosome")
	}
}

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	if got := calculateSequenceComplexity([]byte{}, 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(empty) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_TooShort(t *testing.T) {
	if got := calculateSequenceComplexity([]byte("ATG"), 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(too short) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	// Note: this implementation may return NaN for empty input
	got := calculateGCContent(seqs)
	if got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestCalculateSequenceComplexity_AllSame(t *testing.T) {
	seq := []byte("AAAAAAAAAAAAA")
	got := calculateSequenceComplexity(seq, 4)
	// All same bases: 1 unique kmer out of 10 total -> 0.1
	want := 1.0 / float64(len(seq)-4+1)
	if got-want > 1e-9 || want-got > 1e-9 {
		t.Errorf("calculateSequenceComplexity(all-A) = %f, want %f", got, want)
	}
}

func TestCalculateSequenceComplexity_Range(t *testing.T) {
	seq := []byte("ACGTACGTACGTACGT")
	got := calculateSequenceComplexity(seq, 4)
	if got < 0.0 || got > 1.0 {
		t.Errorf("calculateSequenceComplexity = %f, not in [0,1]", got)
	}
}

func TestCalculateSequenceComplexity_ExactK(t *testing.T) {
	// Length == k: 1 kmer, 1 unique -> 1.0
	seq := []byte("ACGT")
	got := calculateSequenceComplexity(seq, 4)
	if got != 1.0 {
		t.Errorf("calculateSequenceComplexity(ACGT, k=4) = %f, want 1.0", got)
	}
}
