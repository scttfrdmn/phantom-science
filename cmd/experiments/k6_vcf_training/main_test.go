//go:build !release
// +build !release

package main

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"
	"time"

	authenticity "github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestCalculateStateSpaceSize_K1(t *testing.T) {
	// 4^(1-1) = 4^0 = 1
	if got := calculateStateSpaceSize(1); got != 1 {
		t.Errorf("calculateStateSpaceSize(1) = %d, want 1", got)
	}
}

func TestCalculateStateSpaceSize_K2(t *testing.T) {
	// 4^(2-1) = 4^1 = 4
	if got := calculateStateSpaceSize(2); got != 4 {
		t.Errorf("calculateStateSpaceSize(2) = %d, want 4", got)
	}
}

func TestCalculateStateSpaceSize_K4(t *testing.T) {
	// 4^(4-1) = 4^3 = 64
	if got := calculateStateSpaceSize(4); got != 64 {
		t.Errorf("calculateStateSpaceSize(4) = %d, want 64", got)
	}
}

func TestCalculateStateSpaceSize_K6(t *testing.T) {
	// 4^(6-1) = 4^5 = 1024
	if got := calculateStateSpaceSize(6); got != 1024 {
		t.Errorf("calculateStateSpaceSize(6) = %d, want 1024", got)
	}
}

func TestCalculateComplexity_Empty(t *testing.T) {
	if got := calculateComplexity([][]byte{}); got != 0 {
		t.Errorf("calculateComplexity([]) = %f, want 0", got)
	}
}

func TestCalculateComplexity_ShortSeq(t *testing.T) {
	// Sequence shorter than k=6 should contribute 0
	seqs := [][]byte{[]byte("ATCG")}
	if got := calculateComplexity(seqs); got != 0 {
		t.Errorf("calculateComplexity(short seq) = %f, want 0", got)
	}
}

func TestCalculateComplexity_RepetitiveSeq(t *testing.T) {
	// All same base, low complexity
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
	seqs := [][]byte{seq}
	got := calculateComplexity(seqs)
	if got >= 0.5 {
		t.Errorf("calculateComplexity(all-A) = %f, expected < 0.5", got)
	}
}

func TestCalculateComplexity_DiverseSeq(t *testing.T) {
	// Diverse sequence, higher complexity
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	seqs := [][]byte{seq}
	got := calculateComplexity(seqs)
	if got <= 0 {
		t.Errorf("calculateComplexity(diverse) = %f, expected > 0", got)
	}
}

func TestCalculateGCContent_Empty(t *testing.T) {
	if got := calculateGCContent([][]byte{}); got != 0 {
		t.Errorf("calculateGCContent([]) = %f, want 0", got)
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	if got := calculateGCContent(seqs); got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestCalculateGCContent_AllAT(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATAT")}
	if got := calculateGCContent(seqs); got != 0.0 {
		t.Errorf("calculateGCContent(all AT) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_HalfGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCATGCAT")}
	got := calculateGCContent(seqs)
	if got < 0.49 || got > 0.51 {
		t.Errorf("calculateGCContent(half GC) = %f, want ~0.5", got)
	}
}

func TestFormatBytes_Bytes(t *testing.T) {
	if got := formatBytes(500); got != "500 B" {
		t.Errorf("formatBytes(500) = %q, want \"500 B\"", got)
	}
}

func TestFormatBytes_KB(t *testing.T) {
	got := formatBytes(1024)
	if got != "1.0 KB" {
		t.Errorf("formatBytes(1024) = %q, want \"1.0 KB\"", got)
	}
}

func TestFormatBytes_MB(t *testing.T) {
	got := formatBytes(1024 * 1024)
	if got != "1.0 MB" {
		t.Errorf("formatBytes(1MB) = %q, want \"1.0 MB\"", got)
	}
}

func TestFormatBytes_GB(t *testing.T) {
	got := formatBytes(1024 * 1024 * 1024)
	if got != "1.0 GB" {
		t.Errorf("formatBytes(1GB) = %q, want \"1.0 GB\"", got)
	}
}

func TestAbs_Positive(t *testing.T) {
	if got := abs(5.0); got != 5.0 {
		t.Errorf("abs(5.0) = %f, want 5.0", got)
	}
}

func TestAbs_Negative(t *testing.T) {
	if got := abs(-5.0); got != 5.0 {
		t.Errorf("abs(-5.0) = %f, want 5.0", got)
	}
}

func TestAbs_Zero(t *testing.T) {
	if got := abs(0.0); got != 0.0 {
		t.Errorf("abs(0.0) = %f, want 0.0", got)
	}
}

func TestDetectAppleSiliconCores_Positive(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

// ---------- getNextSteps ----------

func TestGetNextSteps_HighComplexity(t *testing.T) {
	// complexity > 0.015 → enhanced control path
	got := getNextSteps(0.02)
	expected := "1. Implement enhanced complexity control\n2. Train on full 11-population dataset\n3. Prepare for production validation"
	if got != expected {
		t.Errorf("getNextSteps(0.02) = %q, want %q", got, expected)
	}
}

func TestGetNextSteps_MediumComplexity(t *testing.T) {
	// complexity in (0.01, 0.015] → retest path
	got := getNextSteps(0.012)
	expected := "1. Implement enhanced complexity control\n2. Retest with complexity improvements\n3. Consider Phase 2 if insufficient"
	if got != expected {
		t.Errorf("getNextSteps(0.012) = %q, want %q", got, expected)
	}
}

func TestGetNextSteps_LowComplexity(t *testing.T) {
	// complexity <= 0.01 → document limitations path
	got := getNextSteps(0.005)
	expected := "1. Document k=6 limitations\n2. Implement Phase 2 structural improvements\n3. Consider context-aware VCF processing"
	if got != expected {
		t.Errorf("getNextSteps(0.005) = %q, want %q", got, expected)
	}
}

func TestGetNextSteps_ExactBoundary015(t *testing.T) {
	// exactly 0.015 is NOT > 0.015, falls to medium branch
	got := getNextSteps(0.015)
	expected := "1. Implement enhanced complexity control\n2. Retest with complexity improvements\n3. Consider Phase 2 if insufficient"
	if got != expected {
		t.Errorf("getNextSteps(0.015) = %q, want %q", got, expected)
	}
}

func TestGetNextSteps_ExactBoundary01(t *testing.T) {
	// exactly 0.01 is NOT > 0.01, falls to low branch
	got := getNextSteps(0.01)
	expected := "1. Document k=6 limitations\n2. Implement Phase 2 structural improvements\n3. Consider context-aware VCF processing"
	if got != expected {
		t.Errorf("getNextSteps(0.01) = %q, want %q", got, expected)
	}
}

// ---------- getSuccessEmoji ----------

func TestGetSuccessEmoji_HighComplexity(t *testing.T) {
	got := getSuccessEmoji(0.02)
	if got != "✅" {
		t.Errorf("getSuccessEmoji(0.02) = %q, want ✅", got)
	}
}

func TestGetSuccessEmoji_MediumComplexity(t *testing.T) {
	got := getSuccessEmoji(0.012)
	if got != "🟡" {
		t.Errorf("getSuccessEmoji(0.012) = %q, want 🟡", got)
	}
}

func TestGetSuccessEmoji_LowComplexity(t *testing.T) {
	got := getSuccessEmoji(0.005)
	if got != "❌" {
		t.Errorf("getSuccessEmoji(0.005) = %q, want ❌", got)
	}
}

func TestGetSuccessEmoji_ExactBoundary015(t *testing.T) {
	// 0.015 is not > 0.015, should be medium
	got := getSuccessEmoji(0.015)
	if got != "🟡" {
		t.Errorf("getSuccessEmoji(0.015) = %q, want 🟡", got)
	}
}

func TestGetSuccessEmoji_ExactBoundary01(t *testing.T) {
	// 0.01 is not > 0.01, should be low
	got := getSuccessEmoji(0.01)
	if got != "❌" {
		t.Errorf("getSuccessEmoji(0.01) = %q, want ❌", got)
	}
}

// ---------- calculateGCContent additional coverage (k6_vcf_training) ----------

func TestCalculateGCContent_K6_EmptySlice(t *testing.T) {
	got := calculateGCContent([][]byte{})
	if got != 0 {
		t.Errorf("calculateGCContent(empty) = %v, want 0", got)
	}
}

func TestCalculateGCContent_K6_EmptySequenceInSlice(t *testing.T) {
	// Empty inner sequence should be skipped; result averaged over all elements
	seqs := [][]byte{{}, []byte("GGCC")}
	got := calculateGCContent(seqs)
	if got < 0 || got > 1 {
		t.Errorf("calculateGCContent with empty inner = %v, want in [0,1]", got)
	}
}

// ---------- savePhase1Results ----------

func TestSavePhase1Results_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	savePhase1Results(tmpDir, 30*time.Second, 0.082, 0.42, 1000000, 5000, 6)
	path := filepath.Join(tmpDir, "phase1_k6_results.txt")
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty results file")
	}
}

// ---------- VCF processing helper ----------

func createTestGzippedVCFk6(t *testing.T) string {
	t.Helper()
	tmpDir := t.TempDir()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Write([]byte("chr1\t1000\t.\tA\tGCT\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t2000\t.\tCT\tC\t.\tPASS\t.\n"))
	gw.Close()
	vcfPath := filepath.Join(tmpDir, "test.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)
	return vcfPath
}

// ---------- processRealVCFFile tests ----------

func TestProcessRealVCFFile_Basic(t *testing.T) {
	vcfPath := createTestGzippedVCFk6(t)
	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processRealVCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 2 {
		t.Errorf("expected 2 variants, got %d", variants)
	}
}

func TestProcessRealVCFFile_NonExistentFile(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processRealVCFFile(model, "/nonexistent/path/test.vcf.gz", "TEST")
	if fileSize != 0 {
		t.Errorf("expected 0 file size for missing file, got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for missing file, got %d", variants)
	}
}

func TestProcessRealVCFFile_NotGzipped(t *testing.T) {
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "notgzip.vcf.gz")
	os.WriteFile(path, []byte("this is not gzip data"), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processRealVCFFile(model, path, "TEST")
	// Should handle error gracefully and return 0,0
	if fileSize != 0 || variants != 0 {
		t.Errorf("expected 0,0 for invalid gzip, got %d,%d", fileSize, variants)
	}
}

// ---------- processPopulationVCFs tests ----------

func TestProcessPopulationVCFs_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "POP")
	os.MkdirAll(popDir, 0750)

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processPopulationVCFs(model, "POP", tmpDir, 1)
	if fileSize != 0 {
		t.Errorf("expected 0 bytes for empty dir, got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for empty dir, got %d", variants)
	}
}

func TestProcessPopulationVCFs_WithVCFs(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	os.MkdirAll(popDir, 0750)

	// Create a test VCF in the population directory
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Write([]byte("chr1\t100\t.\tA\tT\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t200\t.\tG\tC\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t300\t.\tAT\tGC\t.\tPASS\t.\n"))
	gw.Close()
	vcfPath := filepath.Join(popDir, "chr1.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processPopulationVCFs(model, "YRI", tmpDir, 1)
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 3 {
		t.Errorf("expected 3 variants, got %d", variants)
	}
}

// ---------- generateValidationSequences tests ----------

func TestGenerateValidationSequences_Basic(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	trainingSeq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	if err := model.TrainFromSequence(trainingSeq); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}
	seqs, err := generateValidationSequences(model, 1, 100)
	if err != nil {
		t.Fatalf("generateValidationSequences failed: %v", err)
	}
	if len(seqs) != 1 {
		t.Errorf("expected 1 sequence, got %d", len(seqs))
	}
	if len(seqs[0]) != 100 {
		t.Errorf("expected sequence length 100, got %d", len(seqs[0]))
	}
}

func TestGenerateValidationSequences_Multiple(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	trainingSeq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	if err := model.TrainFromSequence(trainingSeq); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}
	seqs, err := generateValidationSequences(model, 3, 50)
	if err != nil {
		t.Fatalf("generateValidationSequences failed: %v", err)
	}
	if len(seqs) != 3 {
		t.Errorf("expected 3 sequences, got %d", len(seqs))
	}
	for i, s := range seqs {
		if len(s) != 50 {
			t.Errorf("sequence[%d] length = %d, want 50", i, len(s))
		}
	}
}

// ---------- processTestPopulations ----------

func TestProcessTestPopulations_EmptyDir_SkipsAll(t *testing.T) {
	// With an empty tmpDir and populations that don't exist as subdirs,
	// processTestPopulations should skip all and return (0, 0).
	model := authenticity.NewAdaptiveHMMModel(6)
	tmpDir := t.TempDir()
	populations := []string{"YRI", "CEU", "CHB"}
	bytes, variants := processTestPopulations(model, tmpDir, populations, 1)
	if bytes != 0 {
		t.Errorf("expected 0 bytes, got %d", bytes)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants, got %d", variants)
	}
}

func TestProcessTestPopulations_EmptyPopulations_ReturnsZero(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(6)
	bytes, variants := processTestPopulations(model, t.TempDir(), []string{}, 1)
	if bytes != 0 || variants != 0 {
		t.Errorf("expected (0, 0) for empty populations, got (%d, %d)", bytes, variants)
	}
}
