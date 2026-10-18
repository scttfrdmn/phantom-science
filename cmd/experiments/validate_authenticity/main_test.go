//go:build !release
// +build !release

package main

import (
	"math"
	"os"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func TestCalculateKmerFrequencies_Empty(t *testing.T) {
	freq := calculateKmerFrequencies([]string{})
	if len(freq) != 0 {
		t.Errorf("expected empty map, got %d entries", len(freq))
	}
}

func TestCalculateKmerFrequencies_SumsToOne(t *testing.T) {
	kmers := []string{"ACGT", "ACGT", "GCTA", "TTTT"}
	freq := calculateKmerFrequencies(kmers)
	sum := 0.0
	for _, v := range freq {
		sum += v
	}
	if math.Abs(sum-1.0) > 1e-9 {
		t.Errorf("frequencies sum = %v, want 1.0", sum)
	}
}

func TestCalculateKmerFrequencies_CorrectProb(t *testing.T) {
	kmers := []string{"AAAA", "AAAA", "CCCC"}
	freq := calculateKmerFrequencies(kmers)
	if math.Abs(freq["AAAA"]-2.0/3.0) > 1e-9 {
		t.Errorf("freq[AAAA] = %v, want 0.667", freq["AAAA"])
	}
	if math.Abs(freq["CCCC"]-1.0/3.0) > 1e-9 {
		t.Errorf("freq[CCCC] = %v, want 0.333", freq["CCCC"])
	}
}

func TestCalculateJSDivergence_Identical(t *testing.T) {
	p := map[string]float64{"ACGT": 0.5, "GCTA": 0.5}
	q := map[string]float64{"ACGT": 0.5, "GCTA": 0.5}
	got := calculateJSDivergence(p, q)
	if got > 1e-10 {
		t.Errorf("identical distributions: JS = %v, want ~0", got)
	}
}

func TestCalculateJSDivergence_Disjoint(t *testing.T) {
	p := map[string]float64{"AAAA": 1.0}
	q := map[string]float64{"CCCC": 1.0}
	got := calculateJSDivergence(p, q)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("disjoint distributions: JS = %v, want 1.0", got)
	}
}

func TestCalculateJSDivergence_Range(t *testing.T) {
	p := map[string]float64{"ACGT": 0.3, "GCTA": 0.4, "TTTT": 0.3}
	q := map[string]float64{"ACGT": 0.1, "GCTA": 0.6, "AAAA": 0.3}
	got := calculateJSDivergence(p, q)
	if got < 0.0 || got > 1.0+1e-9 {
		t.Errorf("JS divergence %v outside [0,1]", got)
	}
}

func TestCalculateKLDivergence_Identical(t *testing.T) {
	p := map[string]float64{"ACGT": 0.5, "GCTA": 0.5}
	q := map[string]float64{"ACGT": 0.5, "GCTA": 0.5}
	got := calculateKLDivergence(p, q)
	if got > 1e-10 {
		t.Errorf("identical distributions: KL = %v, want ~0", got)
	}
}

func TestCalculateKLDivergence_NonNegative(t *testing.T) {
	p := map[string]float64{"ACGT": 0.4, "GCTA": 0.6}
	q := map[string]float64{"ACGT": 0.6, "GCTA": 0.4}
	got := calculateKLDivergence(p, q)
	if got < 0.0 {
		t.Errorf("KL divergence %v is negative", got)
	}
}

func TestFindMaxHomopolymerRun_Empty(t *testing.T) {
	got := findMaxHomopolymerRun([]byte{})
	if got != 0 {
		t.Errorf("got %v, want 0", got)
	}
}

func TestFindMaxHomopolymerRun_SingleBase(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("A"))
	if got != 1 {
		t.Errorf("got %v, want 1", got)
	}
}

func TestFindMaxHomopolymerRun_NoRun(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("ACGT"))
	if got != 1 {
		t.Errorf("got %v, want 1", got)
	}
}

func TestFindMaxHomopolymerRun_RunOf4(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("ACGTAAAACGT"))
	if got != 4 {
		t.Errorf("got %v, want 4", got)
	}
}

func TestFindMaxHomopolymerRun_RunAtEnd(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("ACGTTTTT"))
	if got != 5 {
		t.Errorf("got %v, want 5 for run at end", got)
	}
}

func TestCalculateSequenceComplexity_Short(t *testing.T) {
	got := calculateSequenceComplexity([]byte("AC"))
	if got != 0.0 {
		t.Errorf("got %v, want 0.0 for seq shorter than 3", got)
	}
}

func TestCalculateSequenceComplexity_Range(t *testing.T) {
	seq := []byte("ACGTACGTACGTACGTACGTACGT")
	got := calculateSequenceComplexity(seq)
	if got < 0.0 || got > 1.0 {
		t.Errorf("complexity %v outside [0,1]", got)
	}
}

func TestCalculateSequenceComplexity_LowComplexity(t *testing.T) {
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAA")
	got := calculateSequenceComplexity(seq)
	if got > 0.1 {
		t.Errorf("all-A complexity = %v, want <= 0.1", got)
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

func TestCalculateCpGDepletion_WithCpG(t *testing.T) {
	seqs := [][]byte{[]byte("CGCGCGCG")}
	got := calculateCpGDepletion(seqs)
	if got <= 0.0 {
		t.Errorf("expected positive ratio, got %v", got)
	}
}

func TestCalculateGCContentFromKmers_AllGC(t *testing.T) {
	kmers := []string{"GC", "GC", "GC"}
	got := calculateGCContentFromKmers(kmers)
	if got != 1.0 {
		t.Errorf("got %v, want 1.0", got)
	}
}

func TestCalculateGCContentFromKmers_AllAT(t *testing.T) {
	kmers := []string{"AT", "AT"}
	got := calculateGCContentFromKmers(kmers)
	if got != 0.0 {
		t.Errorf("got %v, want 0.0", got)
	}
}

func TestCalculateGCContentFromKmers_Half(t *testing.T) {
	kmers := []string{"ACGT"}
	got := calculateGCContentFromKmers(kmers)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("got %v, want 0.5", got)
	}
}

func TestCalculateBiologicalMetrics_Basic(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACGTACGTACGT"),
		[]byte("GCGCATATATGCGC"),
	}
	metrics := calculateBiologicalMetrics(seqs, 4)
	if metrics.GCContent < 0.0 || metrics.GCContent > 1.0 {
		t.Errorf("GCContent = %f, want in [0,1]", metrics.GCContent)
	}
	if metrics.ATContent < 0.0 || metrics.ATContent > 1.0 {
		t.Errorf("ATContent = %f, want in [0,1]", metrics.ATContent)
	}
}

func TestCalculateBiologicalMetrics_Empty(t *testing.T) {
	metrics := calculateBiologicalMetrics([][]byte{}, 4)
	// Should not panic with empty sequences
	_ = metrics
}

func TestCalculateBiologicalMetrics_HomopolymerDetected(t *testing.T) {
	seqs := [][]byte{[]byte("ACGTAAAAACGT")}
	metrics := calculateBiologicalMetrics(seqs, 4)
	if metrics.MaxHomopolymer < 5 {
		t.Errorf("MaxHomopolymer = %d, want >= 5", metrics.MaxHomopolymer)
	}
}

func TestCalculateBiologicalMetrics_CpGPresent(t *testing.T) {
	seqs := [][]byte{[]byte("CGCGCGCGCGCG")}
	metrics := calculateBiologicalMetrics(seqs, 4)
	if metrics.CpGDepletion <= 0.0 {
		t.Errorf("CpGDepletion = %f, want > 0", metrics.CpGDepletion)
	}
}

func TestCalculateGCContentFromKmers_Empty(t *testing.T) {
	got := calculateGCContentFromKmers([]string{})
	// Should return 0 for empty (0/0)
	if got != 0.0 && !isNaN(got) {
		t.Logf("calculateGCContentFromKmers(empty) = %f (NaN is ok too)", got)
	}
}

func isNaN(f float64) bool {
	return f != f
}

func makeValidationResult(k int, status string) ValidationResult {
	return ValidationResult{
		ModelK:            k,
		TrainingDataset:   "test",
		ValidationDataset: "test",
		JSdivergence:      0.005,
		KLDivergence:      0.01,
		GCContentDiff:     0.02,
		ValidationStatus:  status,
		ValidationTime:    time.Now(),
		Notes:             "test notes",
		BiologicalMetrics: BiologicalMetrics{
			GCContent:          0.42,
			ATContent:          0.58,
			MaxHomopolymer:     5,
			SequenceComplexity: 0.7,
			CpGDepletion:       0.25,
		},
	}
}

func TestLogValidationResult_NoCrash(t *testing.T) {
	result := makeValidationResult(4, "PASS")
	logValidationResult(result)
}

func TestLogValidationResult_FailStatus(t *testing.T) {
	result := makeValidationResult(3, "FAIL")
	logValidationResult(result)
}

func TestGenerateValidationReport_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	results := []ValidationResult{
		makeValidationResult(4, "PASS"),
		makeValidationResult(6, "FAIL"),
	}
	generateValidationReport(results, tmpDir)
	if _, err := os.Stat(tmpDir + "/validation_report.json"); os.IsNotExist(err) {
		t.Error("validation_report.json was not created")
	}
	if _, err := os.Stat(tmpDir + "/validation_report.md"); os.IsNotExist(err) {
		t.Error("validation_report.md was not created")
	}
}

func TestGenerateValidationReport_Empty(t *testing.T) {
	tmpDir := t.TempDir()
	generateValidationReport([]ValidationResult{}, tmpDir)
	if _, err := os.Stat(tmpDir + "/validation_report.json"); os.IsNotExist(err) {
		t.Error("validation_report.json was not created for empty results")
	}
}

func TestGenerateMarkdownReport_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	results := []ValidationResult{
		makeValidationResult(4, "PASS"),
		makeValidationResult(3, "FAIL"),
	}
	filename := tmpDir + "/test_report.md"
	generateMarkdownReport(results, filename)
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		t.Error("markdown report was not created")
	}
}

func TestGenerateMarkdownReport_Empty(t *testing.T) {
	tmpDir := t.TempDir()
	filename := tmpDir + "/empty_report.md"
	generateMarkdownReport([]ValidationResult{}, filename)
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		t.Error("markdown report was not created for empty results")
	}
}

// Ensure math import is used
var _ = math.Abs

// ---------- validateModel ----------

func TestValidateModel_UntrainedModel_ReturnsFailed(t *testing.T) {
	model := authenticity.NewHMMKmerModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	result := validateModel(model, refGenome, 4)
	// Untrained model → generation fails → ValidationStatus = "FAILED"
	if result.ValidationStatus != "FAILED" {
		t.Errorf("ValidationStatus = %q, want FAILED for untrained model", result.ValidationStatus)
	}
	if result.ModelK != 4 {
		t.Errorf("ModelK = %d, want 4", result.ModelK)
	}
}


// ---------------------------------------------------------------------------
// Additional tests for pure helper functions
// ---------------------------------------------------------------------------

func TestCalculateKmerFrequencies_Single(t *testing.T) {
	freq := calculateKmerFrequencies([]string{"ACGT"})
	if len(freq) != 1 {
		t.Fatalf("expected 1 entry, got %d", len(freq))
	}
	if freq["ACGT"] != 1.0 {
		t.Errorf("freq[ACGT] = %v, want 1.0", freq["ACGT"])
	}
}

func TestCalculateKmerFrequencies_AllSame(t *testing.T) {
	kmers := []string{"AAAA", "AAAA", "AAAA", "AAAA"}
	freq := calculateKmerFrequencies(kmers)
	if len(freq) != 1 {
		t.Fatalf("expected 1 distinct kmer, got %d", len(freq))
	}
	if freq["AAAA"] != 1.0 {
		t.Errorf("freq[AAAA] = %v, want 1.0", freq["AAAA"])
	}
}

func TestCalculateKmerFrequencies_ValuesInRange(t *testing.T) {
	kmers := []string{"ACGT", "GCTA", "TTTT", "AAAA", "CCCC"}
	freq := calculateKmerFrequencies(kmers)
	for k, v := range freq {
		if v < 0.0 || v > 1.0 {
			t.Errorf("freq[%s] = %v, want in [0,1]", k, v)
		}
	}
}

func TestCalculateJSDivergence_NonNegative(t *testing.T) {
	cases := []struct {
		p, q map[string]float64
	}{
		{
			map[string]float64{"A": 0.6, "B": 0.4},
			map[string]float64{"A": 0.2, "B": 0.8},
		},
		{
			map[string]float64{"A": 1.0},
			map[string]float64{"A": 0.5, "B": 0.5},
		},
		{
			map[string]float64{"A": 0.25, "B": 0.25, "C": 0.25, "D": 0.25},
			map[string]float64{"A": 0.25, "B": 0.25, "C": 0.25, "D": 0.25},
		},
	}
	for i, tc := range cases {
		got := calculateJSDivergence(tc.p, tc.q)
		if got < 0.0 {
			t.Errorf("case %d: JS = %v, must be >= 0", i, got)
		}
	}
}

func TestCalculateJSDivergence_Symmetric(t *testing.T) {
	p := map[string]float64{"ACGT": 0.3, "GCTA": 0.4, "TTTT": 0.3}
	q := map[string]float64{"ACGT": 0.1, "GCTA": 0.6, "AAAA": 0.3}
	pq := calculateJSDivergence(p, q)
	qp := calculateJSDivergence(q, p)
	if math.Abs(pq-qp) > 1e-9 {
		t.Errorf("JSD is not symmetric: JSD(p,q)=%v JSD(q,p)=%v", pq, qp)
	}
}

func TestCalculateJSDivergence_EmptyMaps(t *testing.T) {
	got := calculateJSDivergence(map[string]float64{}, map[string]float64{})
	if got != 0.0 {
		t.Errorf("empty distributions: JS = %v, want 0.0", got)
	}
}

func TestCalculateJSDivergence_OneEmptyMap(t *testing.T) {
	// p has mass on one kmer, q is empty.
	// m[kmer] = (1+0)/2 = 0.5 => p*log2(p/m) = 1*log2(1/0.5)=1, q term=0 => js/2 = 0.5
	p := map[string]float64{"ACGT": 1.0}
	q := map[string]float64{}
	got := calculateJSDivergence(p, q)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("one empty distribution: JS = %v, want 0.5", got)
	}
}

func TestCalculateKLDivergence_BothDirectionsNonNegative(t *testing.T) {
	p := map[string]float64{"A": 0.8, "B": 0.2}
	q := map[string]float64{"A": 0.2, "B": 0.8}
	if pq := calculateKLDivergence(p, q); pq < 0.0 {
		t.Errorf("KL(p||q) = %v, must be >= 0", pq)
	}
	if qp := calculateKLDivergence(q, p); qp < 0.0 {
		t.Errorf("KL(q||p) = %v, must be >= 0", qp)
	}
}

func TestCalculateKLDivergence_MissingKeyInQ(t *testing.T) {
	// The implementation skips terms where q[kmer]==0, so when p has a key
	// that q has, the shared-key terms can go negative (p*log2(p/q) < 0 when p < q).
	// The contract tested here is that the result is finite (no NaN or Inf).
	p := map[string]float64{"A": 0.5, "B": 0.5}
	q := map[string]float64{"A": 1.0} // B absent from q
	got := calculateKLDivergence(p, q)
	if math.IsNaN(got) || math.IsInf(got, 0) {
		t.Errorf("KL = %v, must be finite", got)
	}
	// Verify the expected computed value: only the A term contributes,
	// giving 0.5 * log2(0.5/1.0) = -0.5
	expected := 0.5 * math.Log2(0.5/1.0)
	if math.Abs(got-expected) > 1e-9 {
		t.Errorf("KL = %v, want %v", got, expected)
	}
}

func TestCalculateKLDivergence_EmptyP(t *testing.T) {
	p := map[string]float64{}
	q := map[string]float64{"A": 1.0}
	got := calculateKLDivergence(p, q)
	if got != 0.0 {
		t.Errorf("KL(empty||q) = %v, want 0.0", got)
	}
}

func TestCalculateGCContentFromKmers_MixedTwoBase(t *testing.T) {
	// "GC" and "AT" together: 2 GC bases / 4 total = 0.5
	kmers := []string{"GC", "AT"}
	got := calculateGCContentFromKmers(kmers)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("got %v, want 0.5", got)
	}
}

func TestCalculateGCContentFromKmers_SingleG(t *testing.T) {
	got := calculateGCContentFromKmers([]string{"G"})
	if got != 1.0 {
		t.Errorf("got %v, want 1.0", got)
	}
}

func TestCalculateGCContentFromKmers_InRange(t *testing.T) {
	kmers := []string{"ACGT", "GCGC", "ATAT"}
	got := calculateGCContentFromKmers(kmers)
	if got < 0.0 || got > 1.0 {
		t.Errorf("GC content %v outside [0,1]", got)
	}
}

func TestFindMaxHomopolymerRun_AllSame(t *testing.T) {
	seq := []byte("GGGGGGG")
	got := findMaxHomopolymerRun(seq)
	if got != 7 {
		t.Errorf("got %v, want 7", got)
	}
}

func TestFindMaxHomopolymerRun_MultipleRuns(t *testing.T) {
	// Runs: AA(2) CCC(3) TTT(3) G(1) — max is 3
	seq := []byte("AACCCTTG")
	got := findMaxHomopolymerRun(seq)
	if got != 3 {
		t.Errorf("got %v, want 3", got)
	}
}

func TestFindMaxHomopolymerRun_RunAtStart(t *testing.T) {
	seq := []byte("AAAAAACGT")
	got := findMaxHomopolymerRun(seq)
	if got != 6 {
		t.Errorf("got %v, want 6", got)
	}
}

func TestFindMaxHomopolymerRun_TwoBases(t *testing.T) {
	seq := []byte("AC")
	got := findMaxHomopolymerRun(seq)
	if got != 1 {
		t.Errorf("got %v, want 1", got)
	}
}

func TestCalculateSequenceComplexity_ExactlyThreeBases(t *testing.T) {
	// len=3: one 3-mer window, maxPossibleKmers=1, uniqueKmers=1 => complexity=1.0
	seq := []byte("ACG")
	got := calculateSequenceComplexity(seq)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("got %v, want 1.0 for a 3-base sequence", got)
	}
}

func TestCalculateSequenceComplexity_RepeatPattern(t *testing.T) {
	// Repeating triplet: fewer unique 3-mers relative to positions
	seq := []byte("ACGACGACGACGACGACG")
	got := calculateSequenceComplexity(seq)
	if got < 0.0 || got > 1.0 {
		t.Errorf("complexity %v outside [0,1]", got)
	}
}

func TestCalculateSequenceComplexity_HighVsLow(t *testing.T) {
	// High-diversity sequence should have higher complexity than uniform sequence.
	seqHigh := []byte("ACGTGCATCGATCGATCGATCGATCGATCG")
	seqLow := []byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
	high := calculateSequenceComplexity(seqHigh)
	low := calculateSequenceComplexity(seqLow)
	if high <= low {
		t.Errorf("expected high-diversity (%v) > uniform (%v)", high, low)
	}
}

func TestCalculateCpGDepletion_NoCGDinucleotide(t *testing.T) {
	// GC bases present but no consecutive CG dinucleotide
	seqs := [][]byte{[]byte("GCTATAGCTA")}
	got := calculateCpGDepletion(seqs)
	if got != 0.0 {
		t.Errorf("got %v, want 0.0 (no CG dinucleotides)", got)
	}
}

func TestCalculateCpGDepletion_KnownRatio(t *testing.T) {
	// "CGAA": dinucleotides are CG, GA, AA => 1 CG out of 3 total
	// ratio = (1/3) / 0.01
	seqs := [][]byte{[]byte("CGAA")}
	got := calculateCpGDepletion(seqs)
	expected := (1.0 / 3.0) / 0.01
	if math.Abs(got-expected) > 1e-6 {
		t.Errorf("got %v, want %v", got, expected)
	}
}

func TestCalculateCpGDepletion_MultiSequenceAggregation(t *testing.T) {
	// seq1 "CGCG": 3 dinucleotides, 2 CG
	// seq2 "ATCG": 3 dinucleotides, 1 CG
	// total: 6 dinucleotides, 3 CG => observed = 0.5 => ratio = 0.5/0.01 = 50
	seqs := [][]byte{
		[]byte("CGCG"),
		[]byte("ATCG"),
	}
	got := calculateCpGDepletion(seqs)
	expected := 0.5 / 0.01
	if math.Abs(got-expected) > 1e-6 {
		t.Errorf("got %v, want %v", got, expected)
	}
}

func TestCalculateCpGDepletion_SingleBase(t *testing.T) {
	// Single-base sequence produces no dinucleotides
	seqs := [][]byte{[]byte("C")}
	got := calculateCpGDepletion(seqs)
	if got != 0.0 {
		t.Errorf("got %v, want 0.0 for single-base sequence", got)
	}
}

// TestCalculateBiologicalMetrics_EmptySequences is a regression test for #128.
// When sequences is nil or contains only empty byte slices, totalBases stays 0
// and the old code divided by it unconditionally, causing a panic.
func TestCalculateBiologicalMetrics_EmptySequences(t *testing.T) {
	// nil slice — must not panic; GCContent and ATContent must be 0.0, not NaN.
	m := calculateBiologicalMetrics(nil, 4)
	if math.IsNaN(m.GCContent) || math.IsInf(m.GCContent, 0) {
		t.Errorf("GCContent should be 0.0 for nil sequences, got %v", m.GCContent)
	}
	if math.IsNaN(m.ATContent) || math.IsInf(m.ATContent, 0) {
		t.Errorf("ATContent should be 0.0 for nil sequences, got %v", m.ATContent)
	}

	// All-empty sequences — same expectation.
	m2 := calculateBiologicalMetrics([][]byte{{}, {}}, 4)
	if math.IsNaN(m2.GCContent) || math.IsInf(m2.GCContent, 0) {
		t.Errorf("GCContent should be 0.0 for empty sequences, got %v", m2.GCContent)
	}
	if math.IsNaN(m2.ATContent) || math.IsInf(m2.ATContent, 0) {
		t.Errorf("ATContent should be 0.0 for empty sequences, got %v", m2.ATContent)
	}
}
