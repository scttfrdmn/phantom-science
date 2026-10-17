//go:build !release
// +build !release

package main

import (
	"testing"
)

// ---------- calculateGCContent - missing branch coverage ----------

func TestCalculateGCContent_Empty(t *testing.T) {
	// Empty input: totalBases == 0 -> should return 0.0
	got := calculateGCContent([][]byte{})
	if got != 0.0 {
		t.Errorf("calculateGCContent(empty) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_EmptySequence(t *testing.T) {
	// Slice containing an empty byte slice
	got := calculateGCContent([][]byte{{}})
	if got != 0.0 {
		t.Errorf("calculateGCContent([empty]) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_AllAT(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATATAT")}
	got := calculateGCContent(seqs)
	if got != 0.0 {
		t.Errorf("calculateGCContent(all AT) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_Mixed(t *testing.T) {
	// 2 G/C out of 4 bases -> 0.5
	seqs := [][]byte{[]byte("ATGC")}
	got := calculateGCContent(seqs)
	if got != 0.5 {
		t.Errorf("calculateGCContent(ATGC) = %f, want 0.5", got)
	}
}

func TestCalculateGCContent_MultipleSequences(t *testing.T) {
	// seq1: GGGG (4 GC), seq2: AAAA (0 GC) -> 4/8 = 0.5
	seqs := [][]byte{[]byte("GGGG"), []byte("AAAA")}
	got := calculateGCContent(seqs)
	if got != 0.5 {
		t.Errorf("calculateGCContent(GGGG+AAAA) = %f, want 0.5", got)
	}
}

func TestCalculateGCContent_LowercaseNotCounted(t *testing.T) {
	// Lowercase bases are not G or C (uppercase comparison)
	seqs := [][]byte{[]byte("gcgc")}
	got := calculateGCContent(seqs)
	if got != 0.0 {
		t.Errorf("calculateGCContent(lowercase gcgc) = %f, want 0.0 (only uppercase counted)", got)
	}
}

// ---------- calculateSequenceComplexity - additional branches ----------

func TestCalculateSequenceComplexity_KEqualsOne(t *testing.T) {
	seq := []byte("ACGTACGT")
	got := calculateSequenceComplexity(seq, 1)
	// 4 unique 1-mers out of 8 total
	want := 4.0 / 8.0
	if got-want > 1e-9 || want-got > 1e-9 {
		t.Errorf("calculateSequenceComplexity(k=1) = %f, want %f", got, want)
	}
}

func TestCalculateSequenceComplexity_LargeK(t *testing.T) {
	// k == len(seq): 1 kmer total, 1 unique -> 1.0
	seq := []byte("ACGTACGT")
	got := calculateSequenceComplexity(seq, len(seq))
	if got != 1.0 {
		t.Errorf("calculateSequenceComplexity(k=len) = %f, want 1.0", got)
	}
}

// ---------- QuickTestGenome - additional coverage ----------

func TestQuickTestGenome_Empty(t *testing.T) {
	g := &QuickTestGenome{chromosomes: map[string][]byte{}}
	names := g.GetChromosomeNames()
	if len(names) != 0 {
		t.Errorf("GetChromosomeNames() = %d, want 0", len(names))
	}
}

func TestQuickTestChromosome_EmptySequence(t *testing.T) {
	c := &QuickTestChromosome{name: "chr0", sequence: []byte{}}
	if c.GetLength() != 0 {
		t.Errorf("GetLength() = %d, want 0", c.GetLength())
	}
	if len(c.GetSequence()) != 0 {
		t.Errorf("GetSequence() length = %d, want 0", len(c.GetSequence()))
	}
}

func TestQuickTestGenome_LoadChromosome_ReturnsCorrectSequence(t *testing.T) {
	seq := []byte("GCTAGCTA")
	g := &QuickTestGenome{
		chromosomes: map[string][]byte{
			"chrY": seq,
		},
	}
	chr, err := g.LoadChromosome("chrY")
	if err != nil {
		t.Fatalf("LoadChromosome error: %v", err)
	}
	if string(chr.GetSequence()) != string(seq) {
		t.Errorf("GetSequence() = %q, want %q", chr.GetSequence(), seq)
	}
	if chr.GetLength() != len(seq) {
		t.Errorf("GetLength() = %d, want %d", chr.GetLength(), len(seq))
	}
}
