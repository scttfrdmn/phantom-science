//go:build !release
// +build !release

package main

import (
	"os"
	"testing"
	"time"
)

func TestParseTrainingTime_Valid(t *testing.T) {
	got := parseTrainingTime("1h30m")
	expected := 90 * time.Minute
	if got != expected {
		t.Errorf("parseTrainingTime(1h30m) = %v, want %v", got, expected)
	}
}

func TestParseTrainingTime_Invalid(t *testing.T) {
	got := parseTrainingTime("not-a-time")
	if got != 0 {
		t.Errorf("parseTrainingTime(invalid) = %v, want 0", got)
	}
}

func TestParseTrainingTime_Zero(t *testing.T) {
	got := parseTrainingTime("0s")
	if got != 0 {
		t.Errorf("parseTrainingTime(0s) = %v, want 0", got)
	}
}

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	if got := calculateSequenceComplexity([]byte{}, 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(empty) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_TooShort(t *testing.T) {
	if got := calculateSequenceComplexity([]byte("ATG"), 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(3-char with k=4) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_Repetitive(t *testing.T) {
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
	got := calculateSequenceComplexity(seq, 4)
	if got >= 0.5 {
		t.Errorf("calculateSequenceComplexity(all-A) = %f, expected < 0.5", got)
	}
}

func TestCalculateSequenceComplexity_Diverse(t *testing.T) {
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCG")
	got := calculateSequenceComplexity(seq, 4)
	if got <= 0 {
		t.Errorf("calculateSequenceComplexity(diverse) = %f, expected > 0", got)
	}
}

func TestCalculateGCContent_Empty(t *testing.T) {
	if got := calculateGCContent([]byte{}); got != 0.0 {
		t.Errorf("calculateGCContent(empty) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	if got := calculateGCContent([]byte("GCGCGCGC")); got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestCalculateGCContent_AllAT(t *testing.T) {
	if got := calculateGCContent([]byte("ATATATAT")); got != 0.0 {
		t.Errorf("calculateGCContent(all AT) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_HalfGC(t *testing.T) {
	got := calculateGCContent([]byte("GCATGCAT"))
	if got < 0.49 || got > 0.51 {
		t.Errorf("calculateGCContent(half GC) = %f, want ~0.5", got)
	}
}

func TestGetAuthenticVCFPopulations_NotEmpty(t *testing.T) {
	pops := getAuthenticVCFPopulations()
	if len(pops) == 0 {
		t.Error("getAuthenticVCFPopulations() returned empty slice")
	}
}

func TestGetAuthenticVCFPopulations_ValidGCContent(t *testing.T) {
	pops := getAuthenticVCFPopulations()
	for _, pop := range pops {
		if pop.ActualGCContent < 0 || pop.ActualGCContent > 1 {
			t.Errorf("population %s has invalid GC content: %f", pop.Population, pop.ActualGCContent)
		}
	}
}

// TestWriteAuthenticVCFFASTA_CreatesFile exercises writeAuthenticVCFFASTA with
// a temporary file path, covering the full happy-path.
func TestWriteAuthenticVCFFASTA_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	outPath := tmpDir + "/test.fasta"

	pop := AuthenticVCFPopulation{
		Population:       "YRI",
		Region:           "Africa",
		ProcessedVariants: 1000,
		TrainingTime:     2 * time.Hour,
		ActualGCContent:  0.42,
		ActualComplexity: 0.94,
	}
	seq := []byte("ACGTACGTACGTACGT")

	err := writeAuthenticVCFFASTA(outPath, pop, len(seq), 1, seq)
	if err != nil {
		t.Fatalf("writeAuthenticVCFFASTA: unexpected error: %v", err)
	}

	// Verify file was created and is non-empty.
	info, err := os.Stat(outPath)
	if err != nil {
		t.Fatalf("stat output file: %v", err)
	}
	if info.Size() == 0 {
		t.Error("writeAuthenticVCFFASTA wrote an empty file")
	}
}

// TestWriteAuthenticVCFFASTA_ErrorOnBadPath verifies an error is returned when
// the destination directory does not exist.
func TestWriteAuthenticVCFFASTA_ErrorOnBadPath(t *testing.T) {
	pop := AuthenticVCFPopulation{Population: "CEU"}
	err := writeAuthenticVCFFASTA("/nonexistent-dir/out.fasta", pop, 10, 1, []byte("ACGT"))
	if err == nil {
		t.Error("writeAuthenticVCFFASTA with bad path: expected error, got nil")
	}
}

// TestGenerateAuthenticSummaryReport_WritesReport exercises
// generateAuthenticSummaryReport end-to-end using a temp directory.
func TestGenerateAuthenticSummaryReport_WritesReport(t *testing.T) {
	tmpDir := t.TempDir()

	pops := []AuthenticVCFPopulation{
		{
			Population:       "YRI",
			Region:           "Africa",
			ProcessedVariants: 500,
			TrainingTime:     time.Hour,
			ActualGCContent:  0.42,
			ActualComplexity: 0.93,
		},
	}

	// Should not panic; writes a Markdown report file.
	generateAuthenticSummaryReport(tmpDir, pops, 15)

	reportPath := tmpDir + "/AUTHENTIC_VCF_TRAINING_REPORT.md"
	info, err := os.Stat(reportPath)
	if err != nil {
		t.Fatalf("report file not created: %v", err)
	}
	if info.Size() == 0 {
		t.Error("generateAuthenticSummaryReport wrote an empty report")
	}
}

// TestGenerateAuthenticSummaryReport_EmptyPopulations verifies no panic occurs
// when an empty population slice is supplied.
func TestGenerateAuthenticSummaryReport_EmptyPopulations(t *testing.T) {
	tmpDir := t.TempDir()
	generateAuthenticSummaryReport(tmpDir, []AuthenticVCFPopulation{}, 0)
}

// TestWriteAuthenticVCFFASTA_LongSequence exercises the line-wrapping code path
// where the sequence is longer than 80 characters. This ensures the branch where
// i+80 <= len(sequence) (a full 80-char line written without clamping) is hit,
// as well as the final partial-line branch.
func TestWriteAuthenticVCFFASTA_LongSequence(t *testing.T) {
	tmpDir := t.TempDir()
	outPath := tmpDir + "/long.fasta"

	pop := AuthenticVCFPopulation{
		Population:        "GBR",
		Region:            "Europe",
		ProcessedVariants: 2000,
		TrainingTime:      time.Hour,
		ActualGCContent:   0.50,
		ActualComplexity:  0.95,
	}

	// Build a 200-base sequence so the loop produces two full 80-char lines and
	// one partial line, exercising both branches of the end-clamping check.
	seq := make([]byte, 200)
	bases := []byte{'A', 'T', 'C', 'G'}
	for i := range seq {
		seq[i] = bases[i%4]
	}

	err := writeAuthenticVCFFASTA(outPath, pop, len(seq), 2, seq)
	if err != nil {
		t.Fatalf("writeAuthenticVCFFASTA (long seq): unexpected error: %v", err)
	}

	info, err := os.Stat(outPath)
	if err != nil {
		t.Fatalf("stat output file: %v", err)
	}
	if info.Size() == 0 {
		t.Error("writeAuthenticVCFFASTA wrote an empty file for long sequence")
	}
}

// TestWriteAuthenticVCFFASTA_ExactMultipleOf80 writes a sequence of exactly 80
// bases so that a single full chunk is written with no clamping.
func TestWriteAuthenticVCFFASTA_ExactMultipleOf80(t *testing.T) {
	tmpDir := t.TempDir()
	outPath := tmpDir + "/exact80.fasta"

	pop := AuthenticVCFPopulation{
		Population:        "TSI",
		Region:            "Europe",
		ProcessedVariants: 500,
		TrainingTime:      30 * time.Minute,
		ActualGCContent:   0.48,
		ActualComplexity:  0.92,
	}

	seq := make([]byte, 80)
	for i := range seq {
		seq[i] = 'A'
	}

	err := writeAuthenticVCFFASTA(outPath, pop, len(seq), 1, seq)
	if err != nil {
		t.Fatalf("writeAuthenticVCFFASTA (80-char seq): unexpected error: %v", err)
	}
}

// TestWriteAuthenticVCFFASTA_EmptySequence writes a zero-length sequence to
// verify no panic or error occurs (the loop body is never entered).
func TestWriteAuthenticVCFFASTA_EmptySequence(t *testing.T) {
	tmpDir := t.TempDir()
	outPath := tmpDir + "/empty_seq.fasta"

	pop := AuthenticVCFPopulation{
		Population: "LWK",
		Region:     "Africa",
	}

	err := writeAuthenticVCFFASTA(outPath, pop, 0, 1, []byte{})
	if err != nil {
		t.Fatalf("writeAuthenticVCFFASTA (empty seq): unexpected error: %v", err)
	}
}

// TestGenerateAuthenticSummaryReport_BadDir exercises the error-logging branch
// inside generateAuthenticSummaryReport by supplying a non-existent output
// directory. The function logs the error and returns; no panic should occur.
func TestGenerateAuthenticSummaryReport_BadDir(t *testing.T) {
	generateAuthenticSummaryReport("/nonexistent-dir-for-test-abc123", []AuthenticVCFPopulation{
		{Population: "YRI", Region: "Africa"},
	}, 5)
	// Success criterion: no panic.
}

// TestGetAuthenticVCFPopulations_Fields verifies that every population entry
// has non-empty Population and Region strings and a positive variant count.
func TestGetAuthenticVCFPopulations_Fields(t *testing.T) {
	pops := getAuthenticVCFPopulations()
	for _, pop := range pops {
		if pop.Population == "" {
			t.Error("population entry has empty Population field")
		}
		if pop.Region == "" {
			t.Errorf("population %s has empty Region field", pop.Population)
		}
		if pop.ProcessedVariants <= 0 {
			t.Errorf("population %s has non-positive ProcessedVariants: %d",
				pop.Population, pop.ProcessedVariants)
		}
	}
}

// TestCalculateSequenceComplexity_ExactK tests a sequence whose length equals k,
// producing exactly one k-mer and a complexity of 1.0.
func TestCalculateSequenceComplexity_ExactK(t *testing.T) {
	seq := []byte("ACGT") // length == k == 4
	got := calculateSequenceComplexity(seq, 4)
	if got != 1.0 {
		t.Errorf("calculateSequenceComplexity(len==k) = %f, want 1.0", got)
	}
}

// TestCalculateGCContent_SingleBase verifies correct results for single-base sequences.
func TestCalculateGCContent_SingleBase(t *testing.T) {
	tests := []struct {
		base byte
		want float64
	}{
		{'G', 1.0},
		{'C', 1.0},
		{'A', 0.0},
		{'T', 0.0},
	}
	for _, tt := range tests {
		got := calculateGCContent([]byte{tt.base})
		if got != tt.want {
			t.Errorf("calculateGCContent([%c]) = %f, want %f", tt.base, got, tt.want)
		}
	}
}
