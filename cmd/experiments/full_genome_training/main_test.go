//go:build !release
// +build !release

package main

import (
	"math"
	"os"
	"path/filepath"
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

func TestPassFailPadded_Pass(t *testing.T) {
	if got := passFailPadded(true); got != "\u2705 PASS" {
		t.Errorf("passFailPadded(true) = %q", got)
	}
}

func TestPassFailPadded_Fail(t *testing.T) {
	if got := passFailPadded(false); got != "\u274c FAIL" {
		t.Errorf("passFailPadded(false) = %q", got)
	}
}

func TestIsValidKmer(t *testing.T) {
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
		got := isValidKmer(tt.kmer)
		if got != tt.valid {
			t.Errorf("isValidKmer(%q) = %v, want %v", tt.kmer, got, tt.valid)
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

func TestCalculateGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	got := calculateGCContent(seqs)
	if got != 1.0 {
		t.Errorf("got %v, want 1.0", got)
	}
}

func TestCalculateGCContent_AllAT(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATAT")}
	got := calculateGCContent(seqs)
	if got != 0.0 {
		t.Errorf("got %v, want 0.0", got)
	}
}

func TestCalculateGCContent_HalfGC(t *testing.T) {
	seqs := [][]byte{[]byte("ACGT")}
	got := calculateGCContent(seqs)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("got %v, want 0.5", got)
	}
}

func TestCalculateGCContent_MultipleSequences(t *testing.T) {
	seqs := [][]byte{[]byte("AAAA"), []byte("CCCC")}
	got := calculateGCContent(seqs)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("got %v, want 0.5", got)
	}
}

func TestFindMaxHomopolymerRun_Empty(t *testing.T) {
	got := findMaxHomopolymerRun([][]byte{})
	if got != 0 {
		t.Errorf("got %v, want 0", got)
	}
}

func TestFindMaxHomopolymerRun_EmptySeq(t *testing.T) {
	got := findMaxHomopolymerRun([][]byte{{}})
	if got != 0 {
		t.Errorf("got %v, want 0", got)
	}
}

func TestFindMaxHomopolymerRun_SingleBase(t *testing.T) {
	got := findMaxHomopolymerRun([][]byte{[]byte("A")})
	if got != 1 {
		t.Errorf("got %v, want 1", got)
	}
}

func TestFindMaxHomopolymerRun_NoRun(t *testing.T) {
	got := findMaxHomopolymerRun([][]byte{[]byte("ACGT")})
	if got != 1 {
		t.Errorf("got %v, want 1", got)
	}
}

func TestFindMaxHomopolymerRun_RunOf4(t *testing.T) {
	got := findMaxHomopolymerRun([][]byte{[]byte("ACGTAAAACGT")})
	if got != 4 {
		t.Errorf("got %v, want 4", got)
	}
}

func TestFindMaxHomopolymerRun_MaxAcrossSequences(t *testing.T) {
	seqs := [][]byte{
		[]byte("AAACGT"),
		[]byte("ACGTTTTTTCGT"),
	}
	got := findMaxHomopolymerRun(seqs)
	if got != 6 {
		t.Errorf("got %v, want 6", got)
	}
}

func TestCalculateCpGDepletion_NoCpG(t *testing.T) {
	seqs := [][]byte{[]byte("AAAAAAA")}
	got := calculateCpGDepletion(seqs)
	if got != 0.0 {
		t.Errorf("got %v, want 0.0", got)
	}
}

func TestCalculateCpGDepletion_Empty(t *testing.T) {
	got := calculateCpGDepletion([][]byte{})
	if got != 0.0 {
		t.Errorf("got %v, want 0.0", got)
	}
}

func TestCalculateCpGDepletion_NonNegative(t *testing.T) {
	seqs := [][]byte{[]byte("ACGACGACGACG")}
	got := calculateCpGDepletion(seqs)
	if got < 0.0 {
		t.Errorf("CpG depletion %v is negative", got)
	}
}

func TestCalculateCpGDepletion_WithCpG(t *testing.T) {
	seqs := [][]byte{[]byte("CGCGCGCG")}
	got := calculateCpGDepletion(seqs)
	if got <= 0.0 {
		t.Errorf("expected positive ratio, got %v", got)
	}
}

func TestCalculateJSDivergenceFromProbs_Identical(t *testing.T) {
	p := map[string]float64{"ACGT": 0.5, "GCTA": 0.5}
	q := map[string]float64{"ACGT": 0.5, "GCTA": 0.5}
	got := calculateJSDivergenceFromProbs(p, q)
	if got > 1e-10 {
		t.Errorf("identical distributions: JS = %v, want ~0", got)
	}
}

func TestCalculateJSDivergenceFromProbs_Disjoint(t *testing.T) {
	p := map[string]float64{"AAAA": 1.0}
	q := map[string]float64{"CCCC": 1.0}
	got := calculateJSDivergenceFromProbs(p, q)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("disjoint distributions: JS = %v, want 1.0", got)
	}
}

func TestCalculateJSDivergenceFromProbs_Empty(t *testing.T) {
	got := calculateJSDivergenceFromProbs(map[string]float64{}, map[string]float64{})
	if got != 0.0 {
		t.Errorf("empty maps: JS = %v, want 0.0", got)
	}
}

func TestCalculateJSDivergenceFromProbs_Range(t *testing.T) {
	p := map[string]float64{"ACGT": 0.3, "GCTA": 0.4, "TTTT": 0.3}
	q := map[string]float64{"ACGT": 0.1, "GCTA": 0.6, "AAAA": 0.3}
	got := calculateJSDivergenceFromProbs(p, q)
	if got < 0.0 || got > 1.0+1e-9 {
		t.Errorf("JS divergence %v outside [0,1]", got)
	}
}

// --- Tests for NewFullGenomeTrainer ---

func TestNewFullGenomeTrainer_NotNil(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewFullGenomeTrainer(model, refGenome)
	if trainer == nil {
		t.Fatal("NewFullGenomeTrainer returned nil")
	}
}

func TestNewFullGenomeTrainer_FieldsSet(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewFullGenomeTrainer(model, refGenome)
	if trainer.model != model {
		t.Error("trainer.model not set correctly")
	}
	if trainer.refGenome != refGenome {
		t.Error("trainer.refGenome not set correctly")
	}
}

// --- Tests for TrainOnFullGenome ---

// trainableRefGenome creates a temp directory with one minimal FASTA file
// so that TrainFromReferenceGenome can succeed (it calls GetChromosomeNames
// and LoadChromosome, which rely on .fasta files in the dir).
func makeFastaDir(t *testing.T, chrName string, seq string) string {
	t.Helper()
	dir := t.TempDir()
	// Write a minimal FASTA file: ">chrN\nSEQUENCE\n"
	content := ">" + chrName + "\n" + seq + "\n"
	path := filepath.Join(dir, chrName+".fasta")
	if err := os.WriteFile(path, []byte(content), 0600); err != nil {
		t.Fatalf("failed to write fasta file: %v", err)
	}
	return dir
}

func TestTrainOnFullGenome_NoChromosomesSucceeds(t *testing.T) {
	// Use an empty temp dir: no chromosome files can be loaded, so the
	// loop skips all chromosomes. TrainFromReferenceGenome will be called
	// with a genome that has no chromosome names -> it still runs and returns nil.
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewFullGenomeTrainer(model, refGenome)
	stats, err := trainer.TrainOnFullGenome()
	if err != nil {
		t.Fatalf("TrainOnFullGenome returned unexpected error: %v", err)
	}
	if stats == nil {
		t.Fatal("TrainOnFullGenome returned nil stats")
	}
	// With no loadable chromosomes the counts stay at zero.
	if stats.ChromosomesProcessed != 0 {
		t.Errorf("expected 0 chromosomes processed, got %d", stats.ChromosomesProcessed)
	}
}

func TestTrainOnFullGenome_WithOneChromosome(t *testing.T) {
	// Create a valid FASTA file for chr1 so at least one chromosome is processed.
	// A short 20-base sequence is sufficient.
	seq := "ACGTACGTACGTACGTACGT"
	dir := makeFastaDir(t, "chr1", seq)
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(dir)
	trainer := NewFullGenomeTrainer(model, refGenome)
	stats, err := trainer.TrainOnFullGenome()
	if err != nil {
		t.Fatalf("TrainOnFullGenome returned unexpected error: %v", err)
	}
	if stats == nil {
		t.Fatal("nil stats returned")
	}
	// The FASTA file has chr1 which is also in the hardcoded chromosome list,
	// so at least one chromosome should be processed.
	if stats.ChromosomesProcessed < 1 {
		t.Errorf("expected at least 1 chromosome processed, got %d", stats.ChromosomesProcessed)
	}
	if stats.TotalKmers == 0 {
		t.Error("expected non-zero TotalKmers")
	}
}

func TestTrainOnFullGenome_StatsFields(t *testing.T) {
	seq := "ACGTACGTACGTACGTACGT"
	dir := makeFastaDir(t, "chr1", seq)
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(dir)
	trainer := NewFullGenomeTrainer(model, refGenome)
	stats, err := trainer.TrainOnFullGenome()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// GenomeCoverage is in [0, 1]
	if stats.GenomeCoverage < 0 || stats.GenomeCoverage > 1 {
		t.Errorf("GenomeCoverage %v outside [0,1]", stats.GenomeCoverage)
	}
	if stats.TrainingTime <= 0 {
		t.Error("TrainingTime should be positive")
	}
}

// --- Tests for extractFullGenomeKmers ---

func TestExtractFullGenomeKmers_EmptyDir(t *testing.T) {
	// An empty directory means no chromosome files exist, so all LoadChromosome
	// calls fail silently and the result is an empty slice.
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	kmers := extractFullGenomeKmers(refGenome, 4, 100)
	// Should be empty (or at most 0) because no fasta files are available.
	if len(kmers) != 0 {
		t.Errorf("expected 0 kmers from empty dir, got %d", len(kmers))
	}
}

func TestExtractFullGenomeKmers_WithChr1(t *testing.T) {
	// Provide a chr1 file with a long enough sequence to yield k-mers.
	seq := "ACGTACGTACGTACGTACGTACGTACGTACGT"
	dir := makeFastaDir(t, "chr1", seq)
	refGenome := reference.NewSimpleReferenceGenome(dir)
	kmers := extractFullGenomeKmers(refGenome, 4, 10)
	if len(kmers) == 0 {
		t.Error("expected at least one k-mer from chr1 sequence")
	}
	// All returned k-mers should be valid (only A/T/G/C).
	for _, kmer := range kmers {
		if !isValidKmer(kmer) {
			t.Errorf("invalid k-mer returned: %q", kmer)
		}
	}
}

func TestExtractFullGenomeKmers_RespectsMaxKmers(t *testing.T) {
	// extractFullGenomeKmers samples from up to 10 chromosomes and computes
	// sampleStep = len(sequence) / (maxKmers / len(chromosomes)).
	// If maxKmers < 10, integer division yields 0 and causes a divide-by-zero.
	// Use maxKmers >= 10 to keep (maxKmers / 10) >= 1 and avoid that path.
	dna := "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	dir := makeFastaDir(t, "chr1", dna)
	refGenome := reference.NewSimpleReferenceGenome(dir)
	maxKmers := 20 // large enough that maxKmers/10 >= 1
	kmers := extractFullGenomeKmers(refGenome, 4, maxKmers)
	if len(kmers) > maxKmers {
		t.Errorf("got %d kmers, want at most %d", len(kmers), maxKmers)
	}
}

// --- Tests for calculateJSDivergence ---

func TestCalculateJSDivergence_IdenticalDistributions(t *testing.T) {
	// When sequences and refKmers share the same k-mer profile, divergence is ~0.
	dna := []byte("ACGTACGTACGTACGT")
	k := 4
	// Build refKmers from the same sequence.
	var refKmers []string
	for i := 0; i <= len(dna)-k; i++ {
		refKmers = append(refKmers, string(dna[i:i+k]))
	}
	got := calculateJSDivergence([][]byte{dna}, refKmers, k)
	if got < 0 {
		t.Errorf("JS divergence must be non-negative, got %v", got)
	}
	// Identical input should yield a small divergence (near zero).
	if got > 0.01 {
		t.Errorf("expected near-zero divergence for identical data, got %v", got)
	}
}

func TestCalculateJSDivergence_EmptySequences(t *testing.T) {
	// No generated sequences -> genTotal is 0, refKmers provides a distribution.
	refKmers := []string{"ACGT", "ACGT", "GCTA"}
	got := calculateJSDivergence([][]byte{}, refKmers, 4)
	if got < 0 {
		t.Errorf("JS divergence must be non-negative, got %v", got)
	}
}

func TestCalculateJSDivergence_EmptyRefKmers(t *testing.T) {
	// No reference k-mers -> refTotal is 0.
	seqs := [][]byte{[]byte("ACGTACGT")}
	got := calculateJSDivergence(seqs, []string{}, 4)
	if got < 0 {
		t.Errorf("JS divergence must be non-negative, got %v", got)
	}
}

func TestCalculateJSDivergence_NonNegativeRange(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACGTACGTACGT"),
		[]byte("TTTTAAAACCCC"),
	}
	refKmers := []string{"ACGT", "TTTT", "AAAA", "CCCC", "GCTA"}
	got := calculateJSDivergence(seqs, refKmers, 4)
	if got < 0 || got > 1.0+1e-9 {
		t.Errorf("JS divergence %v outside [0, 1]", got)
	}
}

