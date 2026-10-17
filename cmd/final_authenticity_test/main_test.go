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

func TestIsValidKmer_Valid(t *testing.T) {
	if !isValidKmer("ATCG") {
		t.Error("ATCG should be a valid kmer")
	}
	if !isValidKmer("GCGCGCGC") {
		t.Error("GCGCGCGC should be a valid kmer")
	}
}

func TestIsValidKmer_Invalid(t *testing.T) {
	if isValidKmer("ATCGN") {
		t.Error("ATCGN should not be a valid kmer (N is invalid)")
	}
	if isValidKmer("atcg") {
		t.Error("lowercase kmer should not be valid")
	}
	if isValidKmer("ATCG1") {
		t.Error("kmer with digit should not be valid")
	}
}

func TestCalculateSequenceComplexity_High(t *testing.T) {
	// Diverse sequence
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	complexity := calculateSequenceComplexity(seq, 4)
	if complexity <= 0 {
		t.Errorf("Expected positive complexity, got %f", complexity)
	}
}

func TestCalculateSequenceComplexity_Low(t *testing.T) {
	// Repetitive sequence
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
	complexity := calculateSequenceComplexity(seq, 4)
	if complexity > 0.1 {
		t.Errorf("Expected very low complexity for repetitive sequence, got %f", complexity)
	}
}

func TestCalculateSequenceComplexity_TooShort(t *testing.T) {
	seq := []byte("AT")
	complexity := calculateSequenceComplexity(seq, 8)
	if complexity != 0.0 {
		t.Errorf("Expected 0.0 for sequence shorter than k, got %f", complexity)
	}
}

func TestCalculateCpGDepletion_WithCpG(t *testing.T) {
	// Sequence with CpG dinucleotides
	seqs := [][]byte{[]byte("ACGTATCG")} // Has 2 CG dinucleotides
	depletion := calculateCpGDepletion(seqs)
	if depletion <= 0 {
		t.Errorf("Expected positive CpG depletion value, got %f", depletion)
	}
}

func TestCalculateCpGDepletion_NoCpG(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATAT")}
	depletion := calculateCpGDepletion(seqs)
	if depletion != 0.0 {
		t.Errorf("Expected 0.0 for sequence with no CpG, got %f", depletion)
	}
}

func TestCalculateCpGDepletion_Empty(t *testing.T) {
	seqs := [][]byte{}
	depletion := calculateCpGDepletion(seqs)
	if depletion != 0.0 {
		t.Errorf("Expected 0.0 for empty sequences, got %f", depletion)
	}
}

func TestFindMaxHomopolymerRun_Single(t *testing.T) {
	max := findMaxHomopolymerRun([]byte("ATCG"))
	if max != 1 {
		t.Errorf("Expected max 1 for no homopolymers, got %d", max)
	}
}

func TestFindMaxHomopolymerRun_WithRun(t *testing.T) {
	max := findMaxHomopolymerRun([]byte("ATAAAATCG"))
	if max != 4 {
		t.Errorf("Expected max 4 for AAAA run, got %d", max)
	}
}

func TestFindMaxHomopolymerRun_AtEnd(t *testing.T) {
	max := findMaxHomopolymerRun([]byte("ATCGGGGG"))
	if max != 5 {
		t.Errorf("Expected max 5 for GGGGG at end, got %d", max)
	}
}

func TestFindMaxHomopolymerRun_Empty(t *testing.T) {
	max := findMaxHomopolymerRun([]byte{})
	if max != 0 {
		t.Errorf("Expected 0 for empty sequence, got %d", max)
	}
}

func TestCalculateJSDivergenceFromProbs_Identical(t *testing.T) {
	p := map[string]float64{"ATCG": 0.5, "GCTA": 0.5}
	q := map[string]float64{"ATCG": 0.5, "GCTA": 0.5}
	js := calculateJSDivergenceFromProbs(p, q)
	if js > 0.001 {
		t.Errorf("Expected near-zero JS divergence for identical distributions, got %f", js)
	}
}

func TestCalculateJSDivergenceFromProbs_Different(t *testing.T) {
	p := map[string]float64{"ATCG": 1.0}
	q := map[string]float64{"GCTA": 1.0}
	js := calculateJSDivergenceFromProbs(p, q)
	if js <= 0 {
		t.Errorf("Expected positive JS divergence for completely different distributions, got %f", js)
	}
}

func TestCalculateJSDivergence_Basic(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCGATCG")}
	refKmers := []string{"ATCG", "TCGA", "CGAT", "GATC"}
	k := 4
	js := calculateJSDivergence(seqs, refKmers, k)
	if js < 0 {
		t.Errorf("JS divergence should be non-negative, got %f", js)
	}
}

func TestCalculateKLDivergence_Basic(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCGATCG")}
	refKmers := []string{"ATCG", "TCGA", "CGAT", "GATC"}
	k := 4
	kl := calculateKLDivergence(seqs, refKmers, k)
	if math.IsNaN(kl) {
		t.Error("KL divergence should not be NaN")
	}
}

func TestCalculateBiologicalRealismScore_Perfect(t *testing.T) {
	metrics := ModelMetrics{
		GCContentError:      0.0,  // Perfect GC content
		SequenceComplexity: 1.0,  // Max complexity
		MaxHomopolymer:      0,    // No homopolymers
		CpGDepletion:        0.25, // Target CpG
	}
	score := calculateBiologicalRealismScore(metrics)
	if score <= 0 {
		t.Errorf("Expected positive score for perfect metrics, got %f", score)
	}
}

func TestCalculateBiologicalRealismScore_Poor(t *testing.T) {
	metrics := ModelMetrics{
		GCContentError:      1.0, // Terrible GC content
		SequenceComplexity: 0.0, // No complexity
		MaxHomopolymer:      100, // Long homopolymers
		CpGDepletion:        5.0, // Very high CpG
	}
	score := calculateBiologicalRealismScore(metrics)
	// Score should be clamped at 0 for terrible metrics
	if score < 0 {
		t.Errorf("Score should not be negative, got %f", score)
	}
}

// ---------- printDetailedComparison, assessProductionReadiness, generateFinalReport ----------
// These functions just log output; tests verify they don't panic

func TestPrintDetailedComparison_NocrashReadyModels(t *testing.T) {
	m := ModelMetrics{
		ModelName:          "Test",
		JSdivergence:       0.005,
		KLdivergence:       0.05,
		GCContent:          0.42,
		GCContentError:     0.02,
		SequenceComplexity: 0.7,
		MaxHomopolymer:     5,
		CpGDepletion:       0.25,
		BiologicalRealism:  85.0,
		ProductionReady:    true,
	}
	// Should not panic
	printDetailedComparison(m, m, m)
}

func TestPrintDetailedComparison_NocrasNotReadyModels(t *testing.T) {
	poor := ModelMetrics{
		ModelName:          "Poor",
		JSdivergence:       0.5,
		GCContentError:     0.2,
		SequenceComplexity: 0.2,
		MaxHomopolymer:     25,
		ProductionReady:    false,
	}
	printDetailedComparison(poor, poor, poor)
}

func TestAssessProductionReadiness_AllReady(t *testing.T) {
	ready := ModelMetrics{
		ModelName:          "Ready",
		JSdivergence:       0.005,
		GCContent:          0.42,
		GCContentError:     0.01,
		SequenceComplexity: 0.8,
		MaxHomopolymer:     4,
		BiologicalRealism:  90.0,
		ProductionReady:    true,
	}
	// Should not panic
	assessProductionReadiness(ready, ready, ready)
}

func TestAssessProductionReadiness_NoneReady(t *testing.T) {
	notReady := ModelMetrics{
		ModelName:          "NotReady",
		JSdivergence:       0.5,
		GCContentError:     0.3,
		SequenceComplexity: 0.1,
		MaxHomopolymer:     30,
		ProductionReady:    false,
	}
	assessProductionReadiness(notReady, notReady, notReady)
}

func TestGenerateFinalReport_ProductionReady(t *testing.T) {
	ready := ModelMetrics{
		ModelName:          "Best",
		JSdivergence:       0.005,
		BiologicalRealism:  92.0,
		ProductionReady:    true,
	}
	worse := ModelMetrics{
		ModelName:          "Worse",
		JSdivergence:       0.05,
		BiologicalRealism:  70.0,
		ProductionReady:    false,
	}
	// Should not panic
	generateFinalReport(worse, ready, worse)
}

func TestGenerateFinalReport_NoneReady(t *testing.T) {
	notReady := ModelMetrics{
		ModelName:          "NotReady",
		JSdivergence:       0.5,
		BiologicalRealism:  40.0,
		ProductionReady:    false,
	}
	generateFinalReport(notReady, notReady, notReady)
}

// ---------- extractReferenceKmers ----------

func TestExtractReferenceKmers_EmptyDir_ReturnsEmpty(t *testing.T) {
	// With empty tmpDir, chr21 and chr22 cannot be loaded; returns empty
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	kmers := extractReferenceKmers(refGenome, 4, 100)
	if kmers == nil {
		kmers = []string{}
	}
	if len(kmers) != 0 {
		t.Errorf("extractReferenceKmers(empty dir) = %d kmers, want 0", len(kmers))
	}
}

// ---------- calculateComprehensiveMetrics ----------

func TestCalculateComprehensiveMetrics_BasicSequences(t *testing.T) {
	seqs := [][]byte{
		[]byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
		[]byte("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
	}
	refKmers := []string{"ATCG", "TCGA", "CGAT", "GATC", "ATCG"}
	metrics := calculateComprehensiveMetrics(seqs, refKmers, 4, "TestModel")
	if metrics.ModelName != "TestModel" {
		t.Errorf("ModelName = %q, want TestModel", metrics.ModelName)
	}
	if metrics.GCContent < 0 || metrics.GCContent > 1 {
		t.Errorf("GCContent = %f, want in [0,1]", metrics.GCContent)
	}
	if metrics.SequenceComplexity < 0 {
		t.Errorf("SequenceComplexity = %f, should be >= 0", metrics.SequenceComplexity)
	}
}

func TestCalculateComprehensiveMetrics_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")}
	metrics := calculateComprehensiveMetrics(seqs, []string{}, 4, "GCModel")
	if math.Abs(metrics.GCContent-1.0) > 0.01 {
		t.Errorf("GCContent = %f, want ~1.0 for all-GC", metrics.GCContent)
	}
}

// ---------- generateSequences ----------

func TestGenerateSequences_TrainedAdaptiveHMM(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	if err := model.TrainFromSequence(seq); err != nil {
		t.Skipf("model training failed: %v", err)
	}
	seqs := generateSequences(model, 50, 2, "adaptive")
	if len(seqs) != 2 {
		t.Errorf("generateSequences returned %d sequences, want 2", len(seqs))
	}
	for i, s := range seqs {
		if len(s) == 0 {
			t.Errorf("sequence[%d] is empty", i)
		}
	}
}

// TestGenerateSequences_HMMKmerModel covers the HMMKmerModel case in the switch statement
func TestGenerateSequences_HMMKmerModel(t *testing.T) {
	// Build a DemographicGenome to use as a reference for training
	dg := authenticity.NewDemographicGenome("EUR")
	dg.AddChromosome("chr1", 500, nil)

	model := authenticity.NewHMMKmerModel(4)
	if err := model.TrainFromReferenceGenome(dg); err != nil {
		t.Skipf("failed to train HMMKmerModel: %v", err)
	}

	seqs := generateSequences(model, 50, 2, "hmm")
	if len(seqs) != 2 {
		t.Errorf("generateSequences returned %d sequences, want 2", len(seqs))
	}
}

// TestGenerateSequences_EnhancedHMMKmerModel covers the EnhancedHMMKmerModel case in the switch
func TestGenerateSequences_EnhancedHMMKmerModel(t *testing.T) {
	model := authenticity.NewEnhancedHMMKmerModel(4, "EUR")
	err := model.IntegrateDemographicVariants("chr1", 500, nil, nil)
	if err != nil {
		t.Skipf("failed to integrate demographic variants: %v", err)
	}

	seqs := generateSequences(model, 50, 2, "enhanced")
	if len(seqs) != 2 {
		t.Errorf("generateSequences returned %d sequences, want 2", len(seqs))
	}
}

// TestExtractReferenceKmers_WithRealChromosome covers the success path (lines 156-166)
func TestExtractReferenceKmers_WithRealChromosome(t *testing.T) {
	tmpDir := t.TempDir()

	// Create chr21.fasta with valid ACGT sequence
	seq := "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
	content := ">chr21\n" + seq + "\n"
	if err := os.WriteFile(filepath.Join(tmpDir, "chr21.fasta"), []byte(content), 0600); err != nil {
		t.Fatalf("failed to create chr21.fasta: %v", err)
	}

	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	kmers := extractReferenceKmers(refGenome, 4, 10)

	if len(kmers) == 0 {
		t.Error("expected k-mers to be extracted from chr21")
	}
	if len(kmers) > 10 {
		t.Errorf("expected at most 10 k-mers, got %d", len(kmers))
	}
}

// TestCalculateBiologicalRealismScore_HighComplexity covers the complexityScore > 25 clamp branch
func TestCalculateBiologicalRealismScore_HighComplexity(t *testing.T) {
	metrics := ModelMetrics{
		GCContentError:     0.0,
		SequenceComplexity: 2.0, // > 1.0 means complexityScore = 50 > 25, triggers clamp
		MaxHomopolymer:     0,
		CpGDepletion:       0.25,
	}
	score := calculateBiologicalRealismScore(metrics)
	// complexityScore is clamped to 25, so total score = 25 + 25 + 25 + 25 = 100 max
	if score <= 0 {
		t.Errorf("Expected positive score, got %f", score)
	}
}

// TestGenerateFinalReport_OriginalBestModel covers the bestModel = original branch (line 544.61-546.3)
// where enhanced.BiologicalRealism < original.BiologicalRealism
func TestGenerateFinalReport_OriginalBestModel(t *testing.T) {
	original := ModelMetrics{
		ModelName:         "Original",
		JSdivergence:      0.005,
		BiologicalRealism: 95.0,
		ProductionReady:   true,
	}
	enhanced := ModelMetrics{
		ModelName:         "Enhanced",
		JSdivergence:      0.05,
		BiologicalRealism: 70.0,
		ProductionReady:   false,
	}
	adaptive := ModelMetrics{
		ModelName:         "Adaptive",
		JSdivergence:      0.05,
		BiologicalRealism: 80.0,
		ProductionReady:   false,
	}
	// enhanced(70) < original(95) -> bestModel = original; also adaptive(80) < enhanced ... wait
	// adaptive(80) < enhanced(70)? No. So bestModel starts as adaptive,
	// then adaptive(80) < enhanced(70)? No. Then enhanced(70) < original(95)? Yes -> bestModel = original
	generateFinalReport(original, enhanced, adaptive)
}

// Keep unused import used
var _ = math.Pi
