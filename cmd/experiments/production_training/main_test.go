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

// ---------- passFailIcon ----------

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

// ---------- calculateSequenceComplexity ----------

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	if got := calculateSequenceComplexity([]byte{}, 4); got != 0.0 {
		t.Errorf("complexity(empty) = %v, want 0", got)
	}
}

func TestCalculateSequenceComplexity_ShorterThanK(t *testing.T) {
	if got := calculateSequenceComplexity([]byte("AT"), 4); got != 0.0 {
		t.Errorf("complexity(shorter than k) = %v, want 0", got)
	}
}

func TestCalculateSequenceComplexity_SingleKmer(t *testing.T) {
	// seq length == k -> exactly 1 kmer, 1 unique -> 1.0
	got := calculateSequenceComplexity([]byte("ATGC"), 4)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("complexity(single kmer) = %v, want 1.0", got)
	}
}

func TestCalculateSequenceComplexity_RepetitiveSeq(t *testing.T) {
	// All same base, k=4: 1 unique kmer / many total
	seq := []byte("AAAAAAAA") // 8 bases -> 5 kmers, 1 unique
	got := calculateSequenceComplexity(seq, 4)
	want := 1.0 / 5.0
	if math.Abs(got-want) > 1e-9 {
		t.Errorf("complexity(repetitive) = %v, want %v", got, want)
	}
}

func TestCalculateSequenceComplexity_HighComplexity(t *testing.T) {
	seq := []byte("ATCGATCGATCGATCGATCG")
	got := calculateSequenceComplexity(seq, 4)
	if got <= 0 {
		t.Errorf("complexity should be > 0, got %v", got)
	}
	if got > 1.0 {
		t.Errorf("complexity should be <= 1.0, got %v", got)
	}
}

// ---------- calculateGCContent ([][]byte) ----------

func TestCalculateGCContent_Empty(t *testing.T) {
	if got := calculateGCContent([][]byte{}); got != 0.0 {
		t.Errorf("gcContent(empty) = %v, want 0", got)
	}
}

func TestCalculateGCContent_AllAT(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATATAT")}
	if got := calculateGCContent(seqs); got != 0.0 {
		t.Errorf("gcContent(all AT) = %v, want 0", got)
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGCGC")}
	if got := calculateGCContent(seqs); got != 1.0 {
		t.Errorf("gcContent(all GC) = %v, want 1.0", got)
	}
}

func TestCalculateGCContent_Mixed(t *testing.T) {
	seqs := [][]byte{[]byte("AATTGGCC")} // 4 GC out of 8 = 0.5
	got := calculateGCContent(seqs)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("gcContent(AATTGGCC) = %v, want 0.5", got)
	}
}

func TestCalculateGCContent_MultipleSeqs(t *testing.T) {
	seqs := [][]byte{
		[]byte("GGGG"), // 4 GC
		[]byte("AAAA"), // 0 GC
	} // 4 GC out of 8 total = 0.5
	got := calculateGCContent(seqs)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("gcContent(mixed seqs) = %v, want 0.5", got)
	}
}

// ---------- findMaxHomopolymerRun ([][]byte) ----------

func TestFindMaxHomopolymerRun_Empty(t *testing.T) {
	if got := findMaxHomopolymerRun([][]byte{}); got != 0 {
		t.Errorf("maxHomopolymer(empty) = %v, want 0", got)
	}
}

func TestFindMaxHomopolymerRun_NoRun(t *testing.T) {
	seqs := [][]byte{[]byte("ATGC")}
	if got := findMaxHomopolymerRun(seqs); got != 1 {
		t.Errorf("maxHomopolymer(no run) = %v, want 1", got)
	}
}

func TestFindMaxHomopolymerRun_Run(t *testing.T) {
	seqs := [][]byte{[]byte("AAATTTGC")}
	if got := findMaxHomopolymerRun(seqs); got != 3 {
		t.Errorf("maxHomopolymer(AAATTTGC) = %v, want 3", got)
	}
}

func TestFindMaxHomopolymerRun_RunAtEnd(t *testing.T) {
	seqs := [][]byte{[]byte("ATGCCCCC")}
	if got := findMaxHomopolymerRun(seqs); got != 5 {
		t.Errorf("maxHomopolymer(ATGCCCCC) = %v, want 5", got)
	}
}

func TestFindMaxHomopolymerRun_AcrossMultipleSeqs(t *testing.T) {
	seqs := [][]byte{
		[]byte("ATGC"),       // max run = 1
		[]byte("AAAAAATGC"),  // max run = 6
	}
	if got := findMaxHomopolymerRun(seqs); got != 6 {
		t.Errorf("maxHomopolymer(across seqs) = %v, want 6", got)
	}
}

// ---------- calculateCpGDepletion ([][]byte) ----------

func TestCalculateCpGDepletion_Empty(t *testing.T) {
	if got := calculateCpGDepletion([][]byte{}); got != 0.0 {
		t.Errorf("cpgDepletion(empty) = %v, want 0", got)
	}
}

func TestCalculateCpGDepletion_NoCpG(t *testing.T) {
	seqs := [][]byte{[]byte("AAAATTTTAAAATTTT")}
	if got := calculateCpGDepletion(seqs); got != 0.0 {
		t.Errorf("cpgDepletion(no CpG) = %v, want 0", got)
	}
}

func TestCalculateCpGDepletion_WithCpG(t *testing.T) {
	// Build a sequence with known C and G content
	seqs := [][]byte{[]byte("ATCGATCGATCG")}
	got := calculateCpGDepletion(seqs)
	if got <= 0 {
		t.Errorf("cpgDepletion(with CpG) should be > 0, got %v", got)
	}
}

// ---------- extractKmerCounts ----------

func TestExtractKmerCounts_Empty(t *testing.T) {
	counts := extractKmerCounts([][]byte{}, 4)
	if len(counts) != 0 {
		t.Errorf("expected empty counts, got %v", counts)
	}
}

func TestExtractKmerCounts_SingleSeq(t *testing.T) {
	seqs := [][]byte{[]byte("ATGCAT")}
	counts := extractKmerCounts(seqs, 4)
	// Kmers: ATGC, TGCA, GCAT -> 3 unique
	if len(counts) != 3 {
		t.Errorf("expected 3 kmers, got %d: %v", len(counts), counts)
	}
	if counts["ATGC"] != 1 {
		t.Errorf("expected ATGC count=1, got %d", counts["ATGC"])
	}
}

func TestExtractKmerCounts_RepetitiveSeq(t *testing.T) {
	// AAAAAAA k=4: kmer AAAA appears 4 times
	seqs := [][]byte{[]byte("AAAAAAA")}
	counts := extractKmerCounts(seqs, 4)
	if counts["AAAA"] != 4 {
		t.Errorf("expected AAAA count=4, got %d", counts["AAAA"])
	}
}

// ---------- normalizeCounts ----------

func TestNormalizeCounts_Empty(t *testing.T) {
	probs := normalizeCounts(map[string]int{})
	if len(probs) != 0 {
		t.Errorf("expected empty probs, got %v", probs)
	}
}

func TestNormalizeCounts_SingleEntry(t *testing.T) {
	probs := normalizeCounts(map[string]int{"AAAA": 10})
	if math.Abs(probs["AAAA"]-1.0) > 1e-9 {
		t.Errorf("single entry should normalize to 1.0, got %v", probs["AAAA"])
	}
}

func TestNormalizeCounts_MultipleEntries(t *testing.T) {
	counts := map[string]int{"AAAA": 3, "TTTT": 1}
	probs := normalizeCounts(counts)
	if math.Abs(probs["AAAA"]-0.75) > 1e-9 {
		t.Errorf("AAAA prob = %v, want 0.75", probs["AAAA"])
	}
	if math.Abs(probs["TTTT"]-0.25) > 1e-9 {
		t.Errorf("TTTT prob = %v, want 0.25", probs["TTTT"])
	}
}

func TestNormalizeCounts_SumToOne(t *testing.T) {
	counts := map[string]int{"AAAA": 5, "TTTT": 3, "GCGC": 2}
	probs := normalizeCounts(counts)
	total := 0.0
	for _, p := range probs {
		total += p
	}
	if math.Abs(total-1.0) > 1e-9 {
		t.Errorf("probs should sum to 1.0, got %v", total)
	}
}

// ---------- jensenShannonDivergence ----------

func TestJensenShannonDivergence_IdenticalDistributions(t *testing.T) {
	p := map[string]float64{"AAAA": 0.5, "TTTT": 0.5}
	q := map[string]float64{"AAAA": 0.5, "TTTT": 0.5}
	got := jensenShannonDivergence(p, q)
	if math.Abs(got) > 1e-9 {
		t.Errorf("JSD(identical) = %v, want ~0", got)
	}
}

func TestJensenShannonDivergence_CompletelyDifferent(t *testing.T) {
	p := map[string]float64{"AAAA": 1.0}
	q := map[string]float64{"TTTT": 1.0}
	got := jensenShannonDivergence(p, q)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("JSD(completely different) = %v, want 1.0", got)
	}
}

func TestJensenShannonDivergence_NonNegative(t *testing.T) {
	p := map[string]float64{"AAAA": 0.6, "TTTT": 0.4}
	q := map[string]float64{"AAAA": 0.3, "GCGC": 0.7}
	got := jensenShannonDivergence(p, q)
	if got < 0 {
		t.Errorf("JSD should be non-negative, got %v", got)
	}
}

func TestJensenShannonDivergence_BoundedByOne(t *testing.T) {
	p := map[string]float64{"AAAA": 0.6, "TTTT": 0.4}
	q := map[string]float64{"CCCC": 0.5, "GGGG": 0.5}
	got := jensenShannonDivergence(p, q)
	if got > 1.0+1e-9 {
		t.Errorf("JSD should be <= 1.0, got %v", got)
	}
}

// ---------- calculateJensenShannonDivergence ----------

func TestCalculateJensenShannonDivergence_IdenticalSeqs(t *testing.T) {
	seqs := [][]byte{[]byte("ATCGATCGATCGATCG")}
	got := calculateJensenShannonDivergence(seqs, seqs, 4)
	if got < 0 {
		t.Errorf("JSD should be non-negative, got %v", got)
	}
}

func TestCalculateJensenShannonDivergence_EmptySeqs(t *testing.T) {
	got := calculateJensenShannonDivergence([][]byte{}, [][]byte{}, 4)
	if got != 0.0 {
		t.Errorf("JSD(empty,empty) = %v, want 0", got)
	}
}

func TestCalculateJensenShannonDivergence_DifferentSeqs(t *testing.T) {
	p := [][]byte{[]byte("AAAAAAAAAAAAAAAA")}
	q := [][]byte{[]byte("CCCCCCCCCCCCCCCC")}
	got := calculateJensenShannonDivergence(p, q, 4)
	if got < 0 || got > 1.0+1e-9 {
		t.Errorf("JSD should be in [0,1], got %v", got)
	}
}

// ---------- writeProgressReport ----------

func TestWriteProgressReport_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	report := ProgressReport{
		Phase:   "test",
		Status:  "running",
		Message: "test message",
	}
	// Should not panic
	writeProgressReport(report, tmpDir)

	// Verify latest_progress.json was created
	latestPath := tmpDir + "/latest_progress.json"
	if _, err := os.Stat(latestPath); os.IsNotExist(err) {
		t.Error("latest_progress.json was not created")
	}
}

func TestWriteProgressReport_BadDir(t *testing.T) {
	// Writing to a non-existent directory should log a warning but not panic
	writeProgressReport(ProgressReport{Phase: "test"}, "/nonexistent/path/12345")
}

// ---------- loadReferenceSequences ----------

func TestLoadReferenceSequences_EmptyDir_ReturnsError(t *testing.T) {
	// With empty tmpDir, LoadChromosome("chr1") will fail
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	_, err := loadReferenceSequences(refGenome, 5, 1000)
	if err == nil {
		t.Error("loadReferenceSequences(empty refGenome) expected error, got nil")
	}
}

func TestLoadReferenceSequences_WithChr1FASTA(t *testing.T) {
	// Create a GRCh38.fasta file containing chr1 data
	tmpDir := t.TempDir()
	// Build a sequence of 500 bases
	seq := make([]byte, 500)
	for i := range seq {
		bases := []byte("ATCG")
		seq[i] = bases[i%4]
	}
	fastaContent := ">chr1 Homo sapiens chromosome 1\n" + string(seq) + "\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	seqs, err := loadReferenceSequences(refGenome, 3, 100)
	if err != nil {
		t.Fatalf("loadReferenceSequences error: %v", err)
	}
	if len(seqs) == 0 {
		t.Error("expected at least 1 sequence")
	}
	for i, s := range seqs {
		if len(s) == 0 {
			t.Errorf("sequence[%d] is empty", i)
		}
	}
}

// ---------- computeProductionMetrics ----------

func TestComputeProductionMetrics_Basic(t *testing.T) {
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCG") // 32 bases
	seqs := [][]byte{seq}
	m := computeProductionMetrics(seqs, 4, 100*time.Millisecond, 50000)
	if m.Complexity < 0 || m.Complexity > 1 {
		t.Errorf("Complexity %v out of [0,1]", m.Complexity)
	}
	if m.GCContent < 0 || m.GCContent > 1 {
		t.Errorf("GCContent %v out of [0,1]", m.GCContent)
	}
	if m.GCError < 0 {
		t.Errorf("GCError %v < 0", m.GCError)
	}
	if m.MaxHomopolymer < 0 {
		t.Errorf("MaxHomopolymer %v < 0", m.MaxHomopolymer)
	}
	if m.GenerationSpeed <= 0 {
		t.Errorf("GenerationSpeed %v should be > 0 when avgGenTime > 0", m.GenerationSpeed)
	}
	// ComplexityImprovement = complexity / 0.005
	wantImprovement := m.Complexity / 0.005
	if math.Abs(m.ComplexityImprovement-wantImprovement) > 1e-9 {
		t.Errorf("ComplexityImprovement %v != %v", m.ComplexityImprovement, wantImprovement)
	}
}

func TestComputeProductionMetrics_ZeroAvgGenTime(t *testing.T) {
	seq := []byte("ATCGATCG")
	seqs := [][]byte{seq}
	m := computeProductionMetrics(seqs, 4, 0, 50000)
	// GenerationSpeed should be 0 when avgGenTime == 0
	if m.GenerationSpeed != 0.0 {
		t.Errorf("GenerationSpeed should be 0 when avgGenTime=0, got %v", m.GenerationSpeed)
	}
}

// ---------- evaluateProductionPassing ----------

func TestEvaluateProductionPassing_AllPass(t *testing.T) {
	m := ProductionMetrics{
		Complexity:     0.5,  // > 0.02 -> pass
		GCError:        0.01, // < 0.03 -> pass
		MaxHomopolymer: 5,    // < 15 -> pass
		CpGDepletion:   0.3,  // |0.3-0.25|=0.05 < 0.4 -> pass
	}
	r := evaluateProductionPassing(m, 0.01) // jsDivergence=0.01 < 0.02 -> pass
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
	if !r.OverallPass {
		t.Error("expected OverallPass=true")
	}
	if r.OverallStatus != "SUCCESS" {
		t.Errorf("OverallStatus = %q, want SUCCESS", r.OverallStatus)
	}
}

func TestEvaluateProductionPassing_AllFail(t *testing.T) {
	m := ProductionMetrics{
		Complexity:     0.01, // < 0.02 -> fail
		GCError:        0.1,  // > 0.03 -> fail
		MaxHomopolymer: 20,   // > 15 -> fail
		CpGDepletion:   1.0,  // |1.0-0.25|=0.75 > 0.4 -> fail
	}
	r := evaluateProductionPassing(m, 0.05) // jsDivergence > 0.02 -> fail
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
	if r.OverallPass {
		t.Error("expected OverallPass=false")
	}
	if r.OverallStatus != "PARTIAL_SUCCESS" {
		t.Errorf("OverallStatus = %q, want PARTIAL_SUCCESS", r.OverallStatus)
	}
}

func TestEvaluateProductionPassing_NegativeJSDivergence(t *testing.T) {
	// jsDivergence < 0 means "not computed" -> jsPass=false
	m := ProductionMetrics{
		Complexity:     0.5,
		GCError:        0.01,
		MaxHomopolymer: 5,
		CpGDepletion:   0.3,
	}
	r := evaluateProductionPassing(m, -1)
	if r.JSPass {
		t.Error("expected JSPass=false for jsDivergence=-1")
	}
}

// ---------- buildFinalSummary ----------

func TestBuildFinalSummary_ContainsExpectedContent(t *testing.T) {
	m := ProductionMetrics{
		Complexity:            0.1234,
		ComplexityImprovement: 24.68,
		GCContent:             0.42,
		GCError:               0.0,
		MaxHomopolymer:        8,
		CpGDepletion:          0.25,
		GenerationSpeed:       100.0,
	}
	summary := buildFinalSummary("2026-01-01 12:00:00", 5*time.Minute, 1000000, m, 0.015, "SUCCESS")
	if summary == "" {
		t.Error("buildFinalSummary returned empty string")
	}
	if !containsStr(summary, "PRODUCTION FULL GENOME TRAINING SUMMARY") {
		t.Error("summary missing header")
	}
	if !containsStr(summary, "SUCCESS") {
		t.Error("summary missing status")
	}
	if !containsStr(summary, "1000000") {
		t.Error("summary missing kmer count")
	}
}

// containsStr checks if s contains substr (used in tests).
func containsStr(s, substr string) bool {
	return len(s) >= len(substr) && (s == substr ||
		func() bool {
			for i := 0; i <= len(s)-len(substr); i++ {
				if s[i:i+len(substr)] == substr {
					return true
				}
			}
			return false
		}())
}

// ---------- writeFinalSummary ----------

func TestWriteFinalSummary_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	writeFinalSummary(tmpDir, "test summary content")
	path := tmpDir + "/final_summary.txt"
	if _, err := os.Stat(path); os.IsNotExist(err) {
		t.Error("writeFinalSummary did not create final_summary.txt")
	}
}

func TestWriteFinalSummary_BadDir(t *testing.T) {
	// Should log warning but not panic
	writeFinalSummary("/nonexistent/path/xyz123", "test summary")
}

// ---------- buildFinalProgressReport ----------

func TestBuildFinalProgressReport_Fields(t *testing.T) {
	trainingStart := time.Now().Add(-5 * time.Minute)
	trainingTime := 5 * time.Minute
	totalKmers := uint64(1000000)
	m := ProductionMetrics{
		Complexity:            0.1,
		GCContent:             0.42,
		GCError:               0.0,
		MaxHomopolymer:        8,
		CpGDepletion:          0.25,
		GenerationSpeed:       100.0,
	}
	report := buildFinalProgressReport(trainingStart, trainingTime, totalKmers, 24, m, 0.015, "SUCCESS", 20.0)
	if report.Phase != "completed" {
		t.Errorf("Phase = %q, want completed", report.Phase)
	}
	if report.Status != "finished" {
		t.Errorf("Status = %q, want finished", report.Status)
	}
	if report.KmersProcessed != int64(totalKmers) {
		t.Errorf("KmersProcessed = %d, want %d", report.KmersProcessed, totalKmers)
	}
	if report.ChromosomesLoaded != 24 {
		t.Errorf("ChromosomesLoaded = %d, want 24", report.ChromosomesLoaded)
	}
	if report.OverallStatus != "SUCCESS" {
		t.Errorf("OverallStatus = %q, want SUCCESS", report.OverallStatus)
	}
	if report.Message == "" {
		t.Error("Message should not be empty")
	}
}

// ---------- determineOverallStatus ----------

func TestDetermineOverallStatus_Success(t *testing.T) {
	if got := determineOverallStatus(true); got != "SUCCESS" {
		t.Errorf("determineOverallStatus(true) = %q, want SUCCESS", got)
	}
}

func TestDetermineOverallStatus_PartialSuccess(t *testing.T) {
	if got := determineOverallStatus(false); got != "PARTIAL_SUCCESS" {
		t.Errorf("determineOverallStatus(false) = %q, want PARTIAL_SUCCESS", got)
	}
}

// ---------- logProductionResultsSummary ----------

func TestLogProductionResultsSummary_DoesNotPanic(t *testing.T) {
	m := ProductionMetrics{
		Complexity:            0.1,
		ComplexityImprovement: 20.0,
		GCContent:             0.42,
		GCError:               0.0,
		MaxHomopolymer:        8,
		CpGDepletion:          0.25,
		GenerationSpeed:       100.0,
	}
	// Should not panic
	logProductionResultsSummary(5*time.Minute, 1000000, m, 0.015, "SUCCESS")
}

// ---------- logProductionFinalOutcome ----------

func TestLogProductionFinalOutcome_Success(t *testing.T) {
	// Should not panic
	logProductionFinalOutcome(true)
}

func TestLogProductionFinalOutcome_Progress(t *testing.T) {
	// Should not panic
	logProductionFinalOutcome(false)
}

// ---------- logProductionIntroduction ----------

func TestLogProductionIntroduction_DoesNotPanic(t *testing.T) {
	logProductionIntroduction("/tmp/test_training")
}

// ---------- logProductionModelSetup ----------

func TestLogProductionModelSetup_DoesNotPanic(t *testing.T) {
	logProductionModelSetup(4, "data/human", 24)
}

// ---------- logTrainingStepHeader ----------

func TestLogTrainingStepHeader_DoesNotPanic(t *testing.T) {
	logTrainingStepHeader("/tmp/test_training")
}

// ---------- logTrainingComplete ----------

func TestLogTrainingComplete_DoesNotPanic(t *testing.T) {
	logTrainingComplete(5*time.Minute, 1000000)
}

// ---------- logGenerationComplete ----------

func TestLogGenerationComplete_DoesNotPanic(t *testing.T) {
	logGenerationComplete(100*time.Millisecond, 500.0)
}

// ---------- writeStartupProgressReport ----------

func TestWriteStartupProgressReport_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	writeStartupProgressReport(tmpDir)
	if _, err := os.Stat(tmpDir + "/latest_progress.json"); os.IsNotExist(err) {
		t.Error("writeStartupProgressReport did not create latest_progress.json")
	}
}

// ---------- writeTrainingStartReport ----------

func TestWriteTrainingStartReport_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	writeTrainingStartReport(tmpDir)
	if _, err := os.Stat(tmpDir + "/latest_progress.json"); os.IsNotExist(err) {
		t.Error("writeTrainingStartReport did not create latest_progress.json")
	}
}

// ---------- writeValidationStartReport ----------

func TestWriteValidationStartReport_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	writeValidationStartReport(tmpDir, 5, 50000)
	if _, err := os.Stat(tmpDir + "/latest_progress.json"); os.IsNotExist(err) {
		t.Error("writeValidationStartReport did not create latest_progress.json")
	}
}

// ---------- logProductionCompletion ----------

func TestLogProductionCompletion_DoesNotPanic(t *testing.T) {
	logProductionCompletion("/tmp/test_output", 5*time.Minute)
}

// ---------- writeInitProgressReport ----------

func TestWriteInitProgressReport_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	writeInitProgressReport(tmpDir, 24, 4)
	// Should create latest_progress.json for each call
	if _, err := os.Stat(tmpDir + "/latest_progress.json"); os.IsNotExist(err) {
		t.Error("writeInitProgressReport did not create latest_progress.json")
	}
}

// ---------- writeTrainingCompleteReport ----------

func TestWriteTrainingCompleteReport_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	writeTrainingCompleteReport(tmpDir, 5*time.Minute, 1000000)
	if _, err := os.Stat(tmpDir + "/latest_progress.json"); os.IsNotExist(err) {
		t.Error("writeTrainingCompleteReport did not create latest_progress.json")
	}
}

// ---------- logProductionMetricsResults ----------

func TestLogProductionMetricsResults_AllPass(t *testing.T) {
	m := ProductionMetrics{
		Complexity:            0.1,
		ComplexityImprovement: 20.0,
		GCContent:             0.42,
		GCError:               0.0,
		MaxHomopolymer:        5,
		CpGDepletion:          0.25,
	}
	r := ProductionPassResults{ComplexityPass: true, GCPass: true, HomopolymerPass: true, CpGPass: true, JSPass: true}
	// Should not panic
	logProductionMetricsResults(m, r, 0.01)
}

func TestLogProductionMetricsResults_NegativeJS(t *testing.T) {
	m := ProductionMetrics{
		Complexity:            0.01,
		ComplexityImprovement: 2.0,
		GCContent:             0.5,
		GCError:               0.08,
		MaxHomopolymer:        20,
		CpGDepletion:          1.0,
	}
	r := ProductionPassResults{ComplexityPass: false, GCPass: false, HomopolymerPass: false, CpGPass: false, JSPass: false}
	// jsDivergence=-1 -> skip JS log line
	logProductionMetricsResults(m, r, -1)
}

// ---------- calculateCpGDepletion edge cases ----------

func TestCalculateCpGDepletion_OnlyC(t *testing.T) {
	// No G bases -> totalG = 0 -> returns 0.0
	seqs := [][]byte{[]byte("CCCCCCCC")}
	got := calculateCpGDepletion(seqs)
	if got != 0.0 {
		t.Errorf("cpgDepletion(only C) = %v, want 0", got)
	}
}

func TestCalculateCpGDepletion_OnlyG(t *testing.T) {
	// No C bases -> totalC = 0 -> returns 0.0
	seqs := [][]byte{[]byte("GGGGGGGG")}
	got := calculateCpGDepletion(seqs)
	if got != 0.0 {
		t.Errorf("cpgDepletion(only G) = %v, want 0", got)
	}
}

// ---------- loadReferenceSequences edge cases ----------

func TestLoadReferenceSequences_ShortChromosome(t *testing.T) {
	// Create a very short chr1 - only 50 bases, asking for sequences of length 100
	tmpDir := t.TempDir()
	seq := make([]byte, 50)
	for i := range seq {
		seq[i] = []byte("ATCG")[i%4]
	}
	fastaContent := ">chr1 test\n" + string(seq) + "\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	seqs, err := loadReferenceSequences(refGenome, 1, 100)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// Should return at most 1 truncated sequence
	if len(seqs) == 0 {
		t.Error("expected at least 1 sequence from short chromosome")
	}
}

func TestLoadReferenceSequences_ManySequences(t *testing.T) {
	// Create a chr1 sequence long enough to provide multiple sub-sequences
	tmpDir := t.TempDir()
	seq := make([]byte, 1000)
	for i := range seq {
		seq[i] = []byte("ATCG")[i%4]
	}
	fastaContent := ">chr1 test\n" + string(seq) + "\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	seqs, err := loadReferenceSequences(refGenome, 3, 100)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(seqs) == 0 {
		t.Error("expected sequences from long chromosome")
	}
}

func TestLoadReferenceSequences_EmptyChromosome(t *testing.T) {
	// Create a chr1 FASTA with only a header (empty sequence) to trigger
	// the start<0 and end<=start branches.
	tmpDir := t.TempDir()
	// Write a FASTA with an empty sequence body - just the header line
	fastaContent := ">chr1 test\n\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	seqs, err := loadReferenceSequences(refGenome, 1, 100)
	// Either an error or empty sequences is acceptable
	if err != nil {
		// Error path is fine
		return
	}
	// If no error, sequences should be empty (no valid data to return)
	if len(seqs) != 0 {
		t.Errorf("expected 0 sequences from empty chromosome, got %d", len(seqs))
	}
}

func TestLoadReferenceSequences_TinyChromosome(t *testing.T) {
	// Create a chr1 FASTA with 1 base - requesting length=100 should trigger
	// the negative-start branch (chrLength=1 < length=100).
	tmpDir := t.TempDir()
	fastaContent := ">chr1 test\nA\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	// i=0: start=0, start >= chrLength(1) is false, so start stays 0
	// end = 0 + 100 = 100, end > chrLength(1) -> end = 1
	// end(1) > start(0), so we get 1 sequence
	seqs, err := loadReferenceSequences(refGenome, 2, 100)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// With a 1-base chromosome requesting length=100 and 2 sequences:
	// i=0: start=0 (not >=1), end=min(100,1)=1, so 1 byte sequence
	// i=1: start=1*100*10=1000, start>=1, start=1-100=-99 -> start=0, end=100, end>1 -> end=1, end(1)>start(0), 1 byte
	// Actually i=1: start=1000 >= chrLength=1, start = 1-100 = -99, then start<0 -> start=0,
	// end=0+100=100, end>1 -> end=1, end<=start? 1<=0? No. So we get a sequence.
	if len(seqs) == 0 {
		t.Error("expected at least 1 sequence from tiny chromosome")
	}
}

// ---------- jensenShannonDivergence edge case ----------

func TestJensenShannonDivergence_BothZeroProb(t *testing.T) {
	// When both p and q have entries that are 0 the continue branch is taken
	p := map[string]float64{"AAAA": 1.0, "TTTT": 0.0}
	q := map[string]float64{"AAAA": 0.0, "TTTT": 1.0}
	got := jensenShannonDivergence(p, q)
	if got < 0 || got > 1.0+1e-9 {
		t.Errorf("JSD should be in [0,1], got %v", got)
	}
}

// ---------- runValidationAndReport ----------

func TestRunValidationAndReport_WithEmptyRefGenome(t *testing.T) {
	// Use an empty ref genome so loadReferenceSequences returns an error (JS=-1 branch)
	tmpDir := t.TempDir()
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())

	seq := make([]byte, 100)
	for i := range seq {
		seq[i] = []byte("ATCG")[i%4]
	}
	sequences := [][]byte{seq}

	// Should not panic
	runValidationAndReport(
		sequences,
		50*time.Millisecond,
		100,
		4,
		refGenome,
		time.Now().Add(-time.Minute),
		time.Minute,
		uint64(1000),
		24,
		time.Now().Add(-2*time.Minute),
		tmpDir,
	)

	// Verify output files were created
	if _, err := os.Stat(tmpDir + "/latest_progress.json"); os.IsNotExist(err) {
		t.Error("runValidationAndReport did not create latest_progress.json")
	}
	if _, err := os.Stat(tmpDir + "/final_summary.txt"); os.IsNotExist(err) {
		t.Error("runValidationAndReport did not create final_summary.txt")
	}
}

func TestRunValidationAndReport_WithFASTARefGenome(t *testing.T) {
	// Create a FASTA with chr1 so JS divergence is computed (non-error branch)
	tmpDir := t.TempDir()
	refDir := t.TempDir()
	seq := make([]byte, 600)
	for i := range seq {
		seq[i] = []byte("ATCG")[i%4]
	}
	fastaContent := ">chr1 test\n" + string(seq) + "\n"
	if err := os.WriteFile(refDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	refGenome := reference.NewSimpleReferenceGenome(refDir)

	genSeq := make([]byte, 100)
	for i := range genSeq {
		genSeq[i] = []byte("ATCG")[i%4]
	}
	sequences := [][]byte{genSeq}

	// Should not panic
	runValidationAndReport(
		sequences,
		30*time.Millisecond,
		100,
		4,
		refGenome,
		time.Now().Add(-time.Minute),
		time.Minute,
		uint64(5000),
		24,
		time.Now().Add(-2*time.Minute),
		tmpDir,
	)

	// Verify both output files were created
	if _, err := os.Stat(tmpDir + "/latest_progress.json"); os.IsNotExist(err) {
		t.Error("runValidationAndReport (with FASTA) did not create latest_progress.json")
	}
}

// ---------- writeTrainingFailReport ----------

func TestWriteTrainingFailReport_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	trainingStart := time.Now().Add(-10 * time.Second)
	// Should not panic and should create progress files
	writeTrainingFailReport(tmpDir, trainingStart, "Training failed: some error")
	if _, err := os.Stat(tmpDir + "/latest_progress.json"); os.IsNotExist(err) {
		t.Error("writeTrainingFailReport did not create latest_progress.json")
	}
}

func TestWriteTrainingFailReport_BadDir(t *testing.T) {
	// Should log warning but not panic even with bad directory
	writeTrainingFailReport("/nonexistent/path/xyz999", time.Now(), "some error")
}

// ---------- writeGenerationFailReport ----------

func TestWriteGenerationFailReport_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	writeGenerationFailReport(tmpDir, "Generation failed: bad model")
	if _, err := os.Stat(tmpDir + "/latest_progress.json"); os.IsNotExist(err) {
		t.Error("writeGenerationFailReport did not create latest_progress.json")
	}
}

func TestWriteGenerationFailReport_BadDir(t *testing.T) {
	// Should log warning but not panic even with bad directory
	writeGenerationFailReport("/nonexistent/path/xyz888", "some generation error")
}

// ---------- computeJSDivergenceWithLogging ----------

func TestComputeJSDivergenceWithLogging_EmptyRefGenome(t *testing.T) {
	// With an empty temp dir, loadReferenceSequences will fail -> returns -1
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	seqs := [][]byte{[]byte("ATCGATCGATCG")}
	result := computeJSDivergenceWithLogging(seqs, refGenome, 3, 100, 4)
	if result != -1 {
		t.Errorf("expected -1 when reference genome cannot be loaded, got %v", result)
	}
}

func TestComputeJSDivergenceWithLogging_WithFASTA(t *testing.T) {
	// Create a valid FASTA with chr1 so loadReferenceSequences succeeds
	tmpDir := t.TempDir()
	seq := make([]byte, 600)
	for i := range seq {
		seq[i] = []byte("ATCG")[i%4]
	}
	fastaContent := ">chr1 test\n" + string(seq) + "\n"
	if err := os.WriteFile(tmpDir+"/GRCh38.fasta", []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	seqs := [][]byte{[]byte("ATCGATCGATCGATCGATCGATCGATCGATCG")}
	result := computeJSDivergenceWithLogging(seqs, refGenome, 2, 100, 4)
	// Should return a non-negative value when successful
	if result < 0 {
		t.Errorf("expected non-negative JS divergence, got %v", result)
	}
}

// Ensure math import is used
var _ = math.Abs
