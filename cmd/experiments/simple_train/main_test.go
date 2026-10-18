//go:build !release
// +build !release

package main

import (
	"os"
	"testing"
)

// simple_train/main.go contains only file/network functions and main().
// Tests verify the package compiles correctly under the build tag.

func TestPackageCompiles(t *testing.T) {
	// This test verifies the package compiles; file-dependent functions are excluded.
	t.Log("simple_train package compiled successfully")
}

// GetSequence, GetName, GetLength are pure accessor methods on SimpleChromosome.

func TestSimpleChromosome_GetName(t *testing.T) {
	c := &SimpleChromosome{name: "chr1", sequence: []byte("ATGC"), length: 4}
	if got := c.GetName(); got != "chr1" {
		t.Errorf("GetName() = %q, want %q", got, "chr1")
	}
}

func TestSimpleChromosome_GetLength(t *testing.T) {
	c := &SimpleChromosome{name: "chr1", sequence: []byte("ATGC"), length: 4}
	if got := c.GetLength(); got != 4 {
		t.Errorf("GetLength() = %d, want 4", got)
	}
}

func TestSimpleChromosome_GetSequence(t *testing.T) {
	seq := []byte("ATGC")
	c := &SimpleChromosome{name: "chr1", sequence: seq, length: 4}
	got := c.GetSequence()
	if string(got) != "ATGC" {
		t.Errorf("GetSequence() = %q, want %q", got, "ATGC")
	}
}

func TestSimpleChromosome_GetSequence_Empty(t *testing.T) {
	c := &SimpleChromosome{name: "chrX", sequence: nil, length: 0}
	if got := c.GetSequence(); got != nil {
		t.Errorf("GetSequence() on empty chromosome = %v, want nil", got)
	}
}

func TestNewSimpleReferenceGenome_NotNil(t *testing.T) {
	rg := NewSimpleReferenceGenome("/tmp/test")
	if rg == nil {
		t.Fatal("NewSimpleReferenceGenome returned nil")
	}
}

func TestSimpleReferenceGenome_GetChromosomeNames(t *testing.T) {
	rg := NewSimpleReferenceGenome("/tmp/test")
	names := rg.GetChromosomeNames()
	if len(names) == 0 {
		t.Error("GetChromosomeNames() returned empty slice")
	}
	expected := map[string]bool{"chr21": true, "chr22": true}
	for _, n := range names {
		if !expected[n] {
			t.Errorf("unexpected chromosome name: %q", n)
		}
	}
}

// ---------- LoadChromosome ----------

func TestLoadChromosome_CachedChromosome(t *testing.T) {
	// Pre-populate the cache; LoadChromosome should return it without file I/O
	chr := &SimpleChromosome{name: "chr21", sequence: []byte("ATCG"), length: 4}
	rg := &SimpleReferenceGenome{
		chromosomes: map[string]*SimpleChromosome{"chr21": chr},
		refDir:      t.TempDir(),
	}
	got, err := rg.LoadChromosome("chr21")
	if err != nil {
		t.Fatalf("LoadChromosome(cached) error: %v", err)
	}
	if got.GetName() != "chr21" {
		t.Errorf("GetName() = %q, want chr21", got.GetName())
	}
	if string(got.GetSequence()) != "ATCG" {
		t.Errorf("GetSequence() = %q, want ATCG", got.GetSequence())
	}
}

func TestLoadChromosome_MissingFile_ReturnsError(t *testing.T) {
	// Empty tmpDir: no GRCh38.fasta, so LoadChromosome must return an error
	rg := NewSimpleReferenceGenome(t.TempDir())
	_, err := rg.LoadChromosome("chr21")
	if err == nil {
		t.Error("LoadChromosome(missing fasta) should return error, got nil")
	}
}

func TestLoadChromosome_WithFastaFile(t *testing.T) {
	// Create a minimal FASTA file containing chr21
	tmpDir := t.TempDir()
	fastaContent := ">chr21 Homo sapiens chromosome 21\nATCGATCGATCGATCGATCGATCG\n>chr22 end\nGCTA\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	rg := NewSimpleReferenceGenome(tmpDir)
	chr, err := rg.LoadChromosome("chr21")
	if err != nil {
		t.Fatalf("LoadChromosome(with fasta) error: %v", err)
	}
	if chr.GetLength() == 0 {
		t.Error("loaded chromosome has zero length")
	}
}

func TestLoadChromosome_ChromosomeNotInFasta_ReturnsError(t *testing.T) {
	// FASTA file exists but doesn't contain the requested chromosome
	tmpDir := t.TempDir()
	fastaContent := ">chr1 Homo sapiens chromosome 1\nATCGATCGATCG\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	rg := NewSimpleReferenceGenome(tmpDir)
	_, err := rg.LoadChromosome("chr21")
	if err == nil {
		t.Error("LoadChromosome(missing chromosome) should return error")
	}
}
