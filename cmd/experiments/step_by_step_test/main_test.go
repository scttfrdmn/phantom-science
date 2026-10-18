//go:build !release
// +build !release

package main

import (
	"math"
	"os"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func TestPassFailIcon_Pass(t *testing.T) {
	if got := passFailIcon(true); got != "✅" {
		t.Errorf("passFailIcon(true) = %q, want ✅", got)
	}
}

func TestPassFailIcon_Fail(t *testing.T) {
	if got := passFailIcon(false); got != "❌" {
		t.Errorf("passFailIcon(false) = %q, want ❌", got)
	}
}

func TestPassFailPadded_Pass(t *testing.T) {
	if got := passFailPadded(true); got != "✅ PASS" {
		t.Errorf("passFailPadded(true) = %q, want \"✅ PASS\"", got)
	}
}

func TestPassFailPadded_Fail(t *testing.T) {
	if got := passFailPadded(false); got != "❌ FAIL" {
		t.Errorf("passFailPadded(false) = %q, want \"❌ FAIL\"", got)
	}
}

func TestIsValidKmer_Valid(t *testing.T) {
	if !isValidKmer("ATCG") {
		t.Error("ATCG should be valid kmer")
	}
}

func TestIsValidKmer_Invalid(t *testing.T) {
	if isValidKmer("ATCGN") {
		t.Error("ATCGN should not be valid kmer")
	}
}

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	if got := calculateSequenceComplexity([]byte{}, 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(empty) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	if got := calculateGCContent(seqs); got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestFindMaxHomopolymerRun_Empty(t *testing.T) {
	if got := findMaxHomopolymerRun([][]byte{}); got != 0 {
		t.Errorf("findMaxHomopolymerRun([]) = %d, want 0", got)
	}
}

func TestCalculateJSDivergenceFromProbs_Identical(t *testing.T) {
	p := map[string]float64{"AAAA": 0.5, "TTTT": 0.5}
	q := map[string]float64{"AAAA": 0.5, "TTTT": 0.5}
	got := calculateJSDivergenceFromProbs(p, q)
	if math.Abs(got) > 1e-9 {
		t.Errorf("JSD(identical) = %f, want ~0", got)
	}
}

func TestCalculateSequenceComplexity_Diverse(t *testing.T) {
	seq := []byte("ACGTACGTACGTACGT")
	got := calculateSequenceComplexity(seq, 4)
	if got <= 0.0 || got > 1.0 {
		t.Errorf("calculateSequenceComplexity(diverse) = %f, want in (0,1]", got)
	}
}

func TestCalculateSequenceComplexity_AllSame(t *testing.T) {
	seq := make([]byte, 20)
	for i := range seq { seq[i] = 'A' }
	got := calculateSequenceComplexity(seq, 4)
	if got > 0.1 {
		t.Errorf("calculateSequenceComplexity(uniform) = %f, want < 0.1", got)
	}
}

func TestCalculateSequenceComplexity_TooShort(t *testing.T) {
	seq := []byte("ACG")
	if got := calculateSequenceComplexity(seq, 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(too short) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_AllAT(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATAT")}
	if got := calculateGCContent(seqs); got != 0.0 {
		t.Errorf("calculateGCContent(all AT) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_Mixed(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCAT")}
	got := calculateGCContent(seqs)
	if got < 0.0 || got > 1.0 {
		t.Errorf("calculateGCContent(mixed) = %f, want in [0,1]", got)
	}
}

func TestFindMaxHomopolymerRun_SingleBase(t *testing.T) {
	seqs := [][]byte{[]byte("A")}
	if got := findMaxHomopolymerRun(seqs); got != 1 {
		t.Errorf("findMaxHomopolymerRun([A]) = %d, want 1", got)
	}
}

func TestFindMaxHomopolymerRun_LongRun(t *testing.T) {
	seqs := [][]byte{[]byte("ACGTAAAAACGT")}
	if got := findMaxHomopolymerRun(seqs); got != 5 {
		t.Errorf("findMaxHomopolymerRun(5xA) = %d, want 5", got)
	}
}

func TestFindMaxHomopolymerRun_RunAtEnd(t *testing.T) {
	seqs := [][]byte{[]byte("ACGTTTT")}
	if got := findMaxHomopolymerRun(seqs); got != 4 {
		t.Errorf("findMaxHomopolymerRun(run at end) = %d, want 4", got)
	}
}

func TestFindMaxHomopolymerRun_MultipleSeqs(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACGT"),
		[]byte("AAAAACGT"),
	}
	if got := findMaxHomopolymerRun(seqs); got != 5 {
		t.Errorf("findMaxHomopolymerRun(multi seqs) = %d, want 5", got)
	}
}

func TestFindMaxHomopolymerRun_EmptySeqInSlice(t *testing.T) {
	seqs := [][]byte{[]byte{}, []byte("AAAA")}
	if got := findMaxHomopolymerRun(seqs); got != 4 {
		t.Errorf("findMaxHomopolymerRun(with empty) = %d, want 4", got)
	}
}

func TestCalculateCpGDepletion_Empty(t *testing.T) {
	if got := calculateCpGDepletion([][]byte{}); got != 0.0 {
		t.Errorf("calculateCpGDepletion([]) = %f, want 0.0", got)
	}
}

func TestCalculateCpGDepletion_WithCpG(t *testing.T) {
	seqs := [][]byte{[]byte("CGCGCGCG")}
	got := calculateCpGDepletion(seqs)
	if got <= 0.0 {
		t.Errorf("calculateCpGDepletion(CpG-rich) = %f, want > 0", got)
	}
}

func TestCalculateCpGDepletion_NoCpG(t *testing.T) {
	seqs := [][]byte{[]byte("AAATTTAAATTT")}
	if got := calculateCpGDepletion(seqs); got != 0.0 {
		t.Errorf("calculateCpGDepletion(no CpG) = %f, want 0.0", got)
	}
}

func TestCalculateJSDivergence_WithData(t *testing.T) {
	seqs := [][]byte{[]byte("ACGTACGTACGT")}
	refKmers := []string{"ACGT", "CGTA", "GTAC", "TACG", "ACGT"}
	got := calculateJSDivergence(seqs, refKmers, 4)
	if got < 0.0 || got > 1.0+1e-9 {
		t.Errorf("calculateJSDivergence = %f, want in [0,1]", got)
	}
}

func TestCalculateJSDivergence_EmptySeqs(t *testing.T) {
	got := calculateJSDivergence([][]byte{}, []string{"ACGT"}, 4)
	if got < 0.0 || got > 1.0+1e-9 {
		t.Errorf("calculateJSDivergence(empty seqs) = %f", got)
	}
}

func TestIsValidKmer_Empty(t *testing.T) {
	if !isValidKmer("") {
		t.Error("empty kmer should be valid (vacuously true)")
	}
}

func TestIsValidKmer_SingleBase(t *testing.T) {
	for _, base := range []string{"A", "T", "G", "C"} {
		if !isValidKmer(base) {
			t.Errorf("single base %q should be valid kmer", base)
		}
	}
}

// ---------- extractLimitedReferenceKmers ----------

func TestExtractLimitedReferenceKmers_EmptyDir_ReturnsEmpty(t *testing.T) {
	// With empty tmpDir, chr22 cannot be loaded; the function should return empty slice
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	kmers := extractLimitedReferenceKmers(refGenome, 4, 100)
	// When chr22 fails to load, returns empty kmers (no error propagated)
	if kmers == nil {
		kmers = []string{}
	}
	if len(kmers) != 0 {
		t.Errorf("extractLimitedReferenceKmers(empty dir) = %d kmers, want 0", len(kmers))
	}
}

func TestExtractLimitedReferenceKmers_WithChr22FASTA(t *testing.T) {
	// Create a FASTA file for chr22 so the chromosome can be loaded
	tmpDir := t.TempDir()
	fastaContent := ">chr22 Homo sapiens chromosome 22\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	kmers := extractLimitedReferenceKmers(refGenome, 4, 10)
	if len(kmers) == 0 {
		t.Error("extractLimitedReferenceKmers should return kmers when chr22 loads successfully")
	}
	if len(kmers) > 10 {
		t.Errorf("extractLimitedReferenceKmers returned %d kmers, want at most 10", len(kmers))
	}
}

// ---------- computeMetrics ----------

func TestComputeMetrics_Basic(t *testing.T) {
	// Build a simple deterministic sequence to drive all sub-functions.
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCG") // 32 bases
	sequences := [][]byte{seq}
	refKmers := []string{"ATCG", "TCGA", "CGAT", "GATC"}
	m := computeMetrics(sequences, refKmers, 4)
	if m.Complexity < 0 || m.Complexity > 1 {
		t.Errorf("computeMetrics: Complexity %v out of [0,1]", m.Complexity)
	}
	if m.GCContent < 0 || m.GCContent > 1 {
		t.Errorf("computeMetrics: GCContent %v out of [0,1]", m.GCContent)
	}
	if m.GCError < 0 {
		t.Errorf("computeMetrics: GCError %v < 0", m.GCError)
	}
	if m.MaxHomopolymer < 0 {
		t.Errorf("computeMetrics: MaxHomopolymer %v < 0", m.MaxHomopolymer)
	}
	if m.CpGDepletion < 0 {
		t.Errorf("computeMetrics: CpGDepletion %v < 0", m.CpGDepletion)
	}
	if m.JSDivergence < 0 {
		t.Errorf("computeMetrics: JSDivergence %v < 0", m.JSDivergence)
	}
}

func TestComputeMetrics_GCError(t *testing.T) {
	// 50% GC sequence: GCError should be |0.5 - 0.42| = 0.08
	seq := []byte("GGCCAATT") // 4 GC, 4 AT -> 50% GC
	sequences := [][]byte{seq}
	m := computeMetrics(sequences, []string{}, 4)
	want := math.Abs(0.5 - 0.42)
	if math.Abs(m.GCError-want) > 1e-9 {
		t.Errorf("computeMetrics GCError = %v, want ~%v", m.GCError, want)
	}
}

// ---------- evaluatePassing ----------

func TestEvaluatePassing_AllPass(t *testing.T) {
	m := AuthenticityMetrics{
		Complexity:     0.8,  // > 0.5 -> pass
		GCError:        0.02, // < 0.05 -> pass
		MaxHomopolymer: 5,    // < 15 -> pass
		CpGDepletion:   0.3,  // |0.3-0.25|=0.05 < 0.5 -> pass
		JSDivergence:   0.005, // < 0.01 -> pass
	}
	r := evaluatePassing(m)
	if !r.ComplexityPass {
		t.Error("expected ComplexityPass=true")
	}
	if !r.GCPass {
		t.Error("expected GCPass=true")
	}
	if !r.HomopolymerPass {
		t.Error("expected HomopolymerPass=true")
	}
	if !r.CpGPass {
		t.Error("expected CpGPass=true")
	}
	if !r.JSPass {
		t.Error("expected JSPass=true")
	}
	if r.TotalPassed != 4 {
		t.Errorf("TotalPassed = %d, want 4", r.TotalPassed)
	}
}

func TestEvaluatePassing_AllFail(t *testing.T) {
	m := AuthenticityMetrics{
		Complexity:     0.1,  // < 0.5 -> fail
		GCError:        0.1,  // > 0.05 -> fail
		MaxHomopolymer: 20,   // > 15 -> fail
		CpGDepletion:   1.0,  // |1.0-0.25|=0.75 > 0.5 -> fail
		JSDivergence:   0.05, // > 0.01 -> fail
	}
	r := evaluatePassing(m)
	if r.ComplexityPass {
		t.Error("expected ComplexityPass=false")
	}
	if r.GCPass {
		t.Error("expected GCPass=false")
	}
	if r.HomopolymerPass {
		t.Error("expected HomopolymerPass=false")
	}
	if r.CpGPass {
		t.Error("expected CpGPass=false")
	}
	if r.JSPass {
		t.Error("expected JSPass=false")
	}
	if r.TotalPassed != 0 {
		t.Errorf("TotalPassed = %d, want 0", r.TotalPassed)
	}
}

func TestEvaluatePassing_PartialPass(t *testing.T) {
	m := AuthenticityMetrics{
		Complexity:     0.8,  // pass
		GCError:        0.1,  // fail
		MaxHomopolymer: 5,    // pass
		CpGDepletion:   0.3,  // pass (CpGPass not counted in TotalPassed)
		JSDivergence:   0.05, // fail
	}
	r := evaluatePassing(m)
	// TotalPassed counts: complexityPass + gcPass + homopolymerPass + jsPass
	// = 1 + 0 + 1 + 0 = 2
	if r.TotalPassed != 2 {
		t.Errorf("TotalPassed = %d, want 2", r.TotalPassed)
	}
}

// ---------- determineFinalStatus ----------

func TestDetermineFinalStatus_ProductionReady(t *testing.T) {
	if got := determineFinalStatus(4); got != "PRODUCTION READY" {
		t.Errorf("determineFinalStatus(4) = %q, want PRODUCTION READY", got)
	}
}

func TestDetermineFinalStatus_SignificantProgress(t *testing.T) {
	for _, n := range []int{2, 3} {
		if got := determineFinalStatus(n); got != "SIGNIFICANT PROGRESS" {
			t.Errorf("determineFinalStatus(%d) = %q, want SIGNIFICANT PROGRESS", n, got)
		}
	}
}

func TestDetermineFinalStatus_NeedsMajorImprovement(t *testing.T) {
	for _, n := range []int{0, 1} {
		if got := determineFinalStatus(n); got != "NEEDS MAJOR IMPROVEMENT" {
			t.Errorf("determineFinalStatus(%d) = %q, want NEEDS MAJOR IMPROVEMENT", n, got)
		}
	}
}

// ---------- improvementAreas ----------

func TestImprovementAreas_AllPass_NoAreas(t *testing.T) {
	m := AuthenticityMetrics{}
	r := PassResults{ComplexityPass: true, GCPass: true, HomopolymerPass: true, JSPass: true}
	areas := improvementAreas(m, r)
	if len(areas) != 0 {
		t.Errorf("improvementAreas(all pass) = %v, want empty", areas)
	}
}

func TestImprovementAreas_AllFail_FourAreas(t *testing.T) {
	m := AuthenticityMetrics{}
	r := PassResults{ComplexityPass: false, GCPass: false, HomopolymerPass: false, JSPass: false}
	areas := improvementAreas(m, r)
	if len(areas) != 4 {
		t.Errorf("improvementAreas(all fail) = %d areas, want 4", len(areas))
	}
}

func TestImprovementAreas_PartialFail(t *testing.T) {
	m := AuthenticityMetrics{}
	r := PassResults{ComplexityPass: false, GCPass: true, HomopolymerPass: false, JSPass: true}
	areas := improvementAreas(m, r)
	if len(areas) != 2 {
		t.Errorf("improvementAreas(partial fail) = %d areas, want 2", len(areas))
	}
}

func TestImprovementAreas_ContainsExpectedStrings(t *testing.T) {
	m := AuthenticityMetrics{}
	r := PassResults{ComplexityPass: false, GCPass: false, HomopolymerPass: true, JSPass: true}
	areas := improvementAreas(m, r)
	found := map[string]bool{}
	for _, a := range areas {
		found[a] = true
	}
	if !found["Sequence complexity needs improvement"] {
		t.Error("expected 'Sequence complexity needs improvement' in areas")
	}
	if !found["GC content needs adjustment"] {
		t.Error("expected 'GC content needs adjustment' in areas")
	}
}

// ---------- formatPerformanceSummary ----------

func TestFormatPerformanceSummary_Normal(t *testing.T) {
	result := formatPerformanceSummary(10.0, 1000, 50000, 100*1_000_000)
	if result == "" {
		t.Error("formatPerformanceSummary returned empty string")
	}
	// processingRate = 1000/10 = 100
	// ms = 100ms, genRate = 50000/100 = 500
}

func TestFormatPerformanceSummary_ZeroTrainingTime(t *testing.T) {
	// processingRate = 0 (avoid div by zero)
	result := formatPerformanceSummary(0.0, 1000, 50000, 100*1_000_000)
	if result == "" {
		t.Error("formatPerformanceSummary(zero training time) returned empty string")
	}
}

func TestFormatPerformanceSummary_ZeroGenTime(t *testing.T) {
	// avgGenNs = 0 -> msPerSeq = 0 -> genRate = 0 (avoid div by zero)
	result := formatPerformanceSummary(5.0, 500, 50000, 0)
	if result == "" {
		t.Error("formatPerformanceSummary(zero gen time) returned empty string")
	}
}

// ---------- logResultsBanner ----------

func TestLogResultsBanner_DoesNotPanic(t *testing.T) {
	m := AuthenticityMetrics{
		Complexity:     0.7,
		GCError:        0.02,
		MaxHomopolymer: 5,
		JSDivergence:   0.005,
	}
	r := PassResults{ComplexityPass: true, GCPass: true, HomopolymerPass: true, JSPass: true, TotalPassed: 4}
	// Should not panic or error
	logResultsBanner(m, r)
}

func TestLogResultsBanner_FailCase(t *testing.T) {
	m := AuthenticityMetrics{
		Complexity:     0.1,
		GCError:        0.1,
		MaxHomopolymer: 20,
		JSDivergence:   0.05,
	}
	r := PassResults{ComplexityPass: false, GCPass: false, HomopolymerPass: false, JSPass: false, TotalPassed: 0}
	logResultsBanner(m, r)
}

// ---------- logFinalDetermination ----------

func TestLogFinalDetermination_ProductionReady(t *testing.T) {
	m := AuthenticityMetrics{Complexity: 0.8, GCError: 0.02, MaxHomopolymer: 5, JSDivergence: 0.005}
	r := PassResults{ComplexityPass: true, GCPass: true, HomopolymerPass: true, JSPass: true, TotalPassed: 4}
	// Should not panic
	logFinalDetermination(4, m, r)
}

func TestLogFinalDetermination_SignificantProgress_AllFailing(t *testing.T) {
	m := AuthenticityMetrics{Complexity: 0.1, GCError: 0.1, MaxHomopolymer: 20, JSDivergence: 0.05}
	r := PassResults{ComplexityPass: false, GCPass: false, HomopolymerPass: false, JSPass: false, TotalPassed: 2}
	// totalPassed=2 triggers "SIGNIFICANT PROGRESS" with all improvement areas
	logFinalDetermination(2, m, r)
}

func TestLogFinalDetermination_SignificantProgress_SomePassing(t *testing.T) {
	m := AuthenticityMetrics{Complexity: 0.8, GCError: 0.02, MaxHomopolymer: 5, JSDivergence: 0.05}
	r := PassResults{ComplexityPass: true, GCPass: true, HomopolymerPass: true, JSPass: false, TotalPassed: 3}
	logFinalDetermination(3, m, r)
}

func TestLogFinalDetermination_NeedsMajorImprovement(t *testing.T) {
	m := AuthenticityMetrics{Complexity: 0.1, GCError: 0.2, MaxHomopolymer: 25, JSDivergence: 0.1}
	r := PassResults{ComplexityPass: false, GCPass: false, HomopolymerPass: false, JSPass: false, TotalPassed: 1}
	// totalPassed=1 -> default case
	logFinalDetermination(1, m, r)
}

func TestLogFinalDetermination_ZeroPassed(t *testing.T) {
	m := AuthenticityMetrics{}
	r := PassResults{}
	logFinalDetermination(0, m, r)
}

// ---------- logIntroduction ----------

func TestLogIntroduction_DoesNotPanic(t *testing.T) {
	logIntroduction()
}

// ---------- logModelConfig ----------

func TestLogModelConfig_DoesNotPanic(t *testing.T) {
	logModelConfig(4)
}

// ---------- logJSDivergenceHeader ----------

func TestLogJSDivergenceHeader_DoesNotPanic(t *testing.T) {
	logJSDivergenceHeader(12345)
}

func TestLogJSDivergenceHeader_ZeroKmers(t *testing.T) {
	logJSDivergenceHeader(0)
}

// ---------- logModelDetails ----------

func TestLogModelDetails_Trained(t *testing.T) {
	// Should not panic
	logModelDetails(4, 256, true)
}

func TestLogModelDetails_NotTrained(t *testing.T) {
	logModelDetails(4, 0, false)
}

// ---------- logStep5Header ----------

func TestLogStep5Header_DoesNotPanic(t *testing.T) {
	logStep5Header()
}

// ---------- logCompletionFooter ----------

func TestLogCompletionFooter_DoesNotPanic(t *testing.T) {
	logCompletionFooter()
}

// ---------- logMetricsResults ----------

func TestLogMetricsResults_AllPass(t *testing.T) {
	m := AuthenticityMetrics{
		Complexity:     0.8,
		GCContent:      0.42,
		GCError:        0.0,
		MaxHomopolymer: 5,
		CpGDepletion:   0.25,
		JSDivergence:   0.005,
	}
	r := PassResults{ComplexityPass: true, GCPass: true, HomopolymerPass: true, CpGPass: true, JSPass: true}
	// Should not panic
	logMetricsResults(m, r)
}

func TestLogMetricsResults_AllFail(t *testing.T) {
	m := AuthenticityMetrics{
		Complexity:     0.1,
		GCContent:      0.5,
		GCError:        0.08,
		MaxHomopolymer: 20,
		CpGDepletion:   1.0,
		JSDivergence:   0.05,
	}
	r := PassResults{ComplexityPass: false, GCPass: false, HomopolymerPass: false, CpGPass: false, JSPass: false}
	logMetricsResults(m, r)
}

// ---------- logPerformanceSummary9 ----------

func TestLogPerformanceSummary9_Normal(t *testing.T) {
	// Should not panic
	logPerformanceSummary9(5*time.Second, 10000, 50000, 3, 100*time.Millisecond)
}

func TestLogPerformanceSummary9_ZeroTimes(t *testing.T) {
	// Both zero times -> avoid div by zero
	logPerformanceSummary9(0, 0, 50000, 3, 0)
}

// Ensure math and os imports are used
var _ = math.Abs
var _ = os.Stat
