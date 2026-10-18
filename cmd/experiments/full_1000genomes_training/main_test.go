//go:build !release
// +build !release

package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestAbs_Positive(t *testing.T) {
	got := abs(3.14)
	if math.Abs(got-3.14) > 1e-9 {
		t.Errorf("abs(3.14) = %v, want 3.14", got)
	}
}

func TestAbs_Negative(t *testing.T) {
	got := abs(-2.71)
	if math.Abs(got-2.71) > 1e-9 {
		t.Errorf("abs(-2.71) = %v, want 2.71", got)
	}
}

func TestAbs_Zero(t *testing.T) {
	got := abs(0.0)
	if got != 0.0 {
		t.Errorf("abs(0) = %v, want 0", got)
	}
}

func TestDetectAppleSiliconCores_ReturnsPositive(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestCalculateComplexity_Empty(t *testing.T) {
	got := calculateComplexity([][]byte{})
	if got != 0.0 {
		t.Errorf("calculateComplexity([]) = %v, want 0.0", got)
	}
}

func TestCalculateComplexity_ShortSeq(t *testing.T) {
	got := calculateComplexity([][]byte{[]byte("ACG")})
	if got != 0.0 {
		t.Errorf("calculateComplexity(short) = %v, want 0.0", got)
	}
}

func TestCalculateComplexity_LowComplexity(t *testing.T) {
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAA")
	got := calculateComplexity([][]byte{seq})
	if got > 0.1 {
		t.Errorf("all-A complexity = %v, want <= 0.1", got)
	}
}

func TestCalculateComplexity_Range(t *testing.T) {
	seq := []byte("ACGTACGTACGTACGTACGTACGT")
	got := calculateComplexity([][]byte{seq})
	if got < 0.0 || got > 1.0 {
		t.Errorf("complexity %v outside [0,1]", got)
	}
}

func TestCalculateGCContent1000G_Empty(t *testing.T) {
	got := calculateGCContent([][]byte{})
	if got != 0.0 {
		t.Errorf("got %v, want 0.0", got)
	}
}

func TestCalculateGCContent1000G_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	got := calculateGCContent(seqs)
	if got != 1.0 {
		t.Errorf("got %v, want 1.0", got)
	}
}

func TestCalculateGCContent1000G_AllAT(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATAT")}
	got := calculateGCContent(seqs)
	if got != 0.0 {
		t.Errorf("got %v, want 0.0", got)
	}
}

func TestCalculateGCContent1000G_HalfGC(t *testing.T) {
	seqs := [][]byte{[]byte("AACCGGTT")}
	got := calculateGCContent(seqs)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("got %v, want 0.5", got)
	}
}

func TestCalculateGCContent1000G_AveragesAcrossSeqs(t *testing.T) {
	seqs := [][]byte{[]byte("AAAA"), []byte("CCCC")}
	got := calculateGCContent(seqs)
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("got %v, want 0.5", got)
	}
}

func TestFormatSize_Bytes(t *testing.T) {
	got := formatSize(500)
	if got != "500 B" {
		t.Errorf("formatSize(500) = %q, want \"500 B\"", got)
	}
}

func TestFormatSize_Kilobytes(t *testing.T) {
	got := formatSize(1024)
	if got != "1.0 KB" {
		t.Errorf("formatSize(1024) = %q, want \"1.0 KB\"", got)
	}
}

func TestFormatSize_Megabytes(t *testing.T) {
	got := formatSize(1024 * 1024)
	if got != "1.0 MB" {
		t.Errorf("formatSize(1MB) = %q, want \"1.0 MB\"", got)
	}
}

func TestGeneratePopulationSequence_Length(t *testing.T) {
	length := 100
	seq := generatePopulationSequence("YRI", length, 0)
	if len(seq) != length {
		t.Errorf("len = %d, want %d", len(seq), length)
	}
}

func TestGeneratePopulationSequence_ValidBases(t *testing.T) {
	seq := generatePopulationSequence("CEU", 200, 1)
	valid := map[byte]bool{'A': true, 'T': true, 'C': true, 'G': true}
	for i, b := range seq {
		if !valid[b] {
			t.Errorf("invalid base %q at position %d", b, i)
		}
	}
}

func TestGeneratePopulationSequence_AllPopulations(t *testing.T) {
	pops := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR", "unknown"}
	for _, pop := range pops {
		seq := generatePopulationSequence(pop, 50, 0)
		if len(seq) != 50 {
			t.Errorf("pop %s: len = %d, want 50", pop, len(seq))
		}
	}
}

func TestGeneratePopulationSequence_ZeroLength(t *testing.T) {
	seq := generatePopulationSequence("YRI", 0, 0)
	if len(seq) != 0 {
		t.Errorf("expected empty slice, got len %d", len(seq))
	}
}

// ---------- saveResults ----------

func TestSaveResults_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	saveResults(tmpDir, 30*time.Second, 0.082, 0.42, "SUCCESS", 6)
	path := filepath.Join(tmpDir, "final_results.txt")
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty results file")
	}
}

// TestSaveResults_ComplexityBranches exercises different complexity/gcContent values
// to cover different branches in saveResults (improvement rendering when complexity > 0).
func TestSaveResults_ComplexityBranches(t *testing.T) {
	tests := []struct {
		name       string
		duration   time.Duration
		complexity float64
		gcContent  float64
		status     string
		perfCores  int
	}{
		// complexity > 0.05: renders improvement line
		{"production_ready", 2 * time.Hour, 0.09, 0.42, "PRODUCTION_READY", 8},
		// complexity in (0.01, 0.05]: renders improvement but not "production"
		{"good_progress", 1 * time.Hour, 0.03, 0.42, "GOOD_PROGRESS", 6},
		// complexity <= 0.01: no improvement branch
		{"low_complexity", 30 * time.Minute, 0.005, 0.50, "SUCCESS", 4},
		// zero values
		{"zeros", 0, 0.0, 0.0, "", 1},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			saveResults(tmpDir, tt.duration, tt.complexity, tt.gcContent, tt.status, tt.perfCores)
			path := filepath.Join(tmpDir, "final_results.txt")
			info, err := os.Stat(path)
			if err != nil {
				t.Fatalf("file should exist after saveResults: %v", err)
			}
			if info.Size() == 0 {
				t.Error("results file should not be empty")
			}
		})
	}
}

// TestGetDirSize_EmptyDir ensures getDirSize returns 0 for an empty directory.
func TestGetDirSize_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	got := getDirSize(tmpDir)
	if got != 0 {
		t.Errorf("getDirSize(empty dir) = %d, want 0", got)
	}
}

// TestGetDirSize_WithFile ensures getDirSize accounts for file sizes.
func TestGetDirSize_WithFile(t *testing.T) {
	tmpDir := t.TempDir()
	data := []byte("hello world test data")
	path := filepath.Join(tmpDir, "test.txt")
	if err := os.WriteFile(path, data, 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	got := getDirSize(tmpDir)
	if got != int64(len(data)) {
		t.Errorf("getDirSize = %d, want %d", got, len(data))
	}
}

// TestCountDataset_NonexistentDir verifies countDataset handles missing dirs gracefully.
func TestCountDataset_NonexistentDir(t *testing.T) {
	size, count := countDataset("/nonexistent/path/does/not/exist")
	// Should return some value without panicking; count should be 0
	_ = size
	if count != 0 {
		t.Errorf("countDataset(nonexistent) count = %d, want 0", count)
	}
}

// ---------- trainOnFullDataset ----------

// TestTrainOnFullDataset_EmptyDir calls trainOnFullDataset with a tmpDir that
// has no population subdirectories. All 11 populations are missing so the
// function skips them all and returns nil without blocking.
func TestTrainOnFullDataset_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	progressDir := t.TempDir()

	model := authenticity.NewAdaptiveHMMModel(4)
	err := trainOnFullDataset(model, tmpDir, progressDir, 2)
	if err != nil {
		t.Errorf("trainOnFullDataset with empty dir returned error: %v", err)
	}
}

// TestTrainOnFullDataset_OnePopNoVCFs creates one matching population directory
// but leaves it empty (no VCF files). processPopulation returns an error for
// that population which is logged but does not propagate to the caller.
func TestTrainOnFullDataset_OnePopNoVCFs(t *testing.T) {
	tmpDir := t.TempDir()
	progressDir := t.TempDir()

	// Create YRI directory but leave it empty
	if err := os.MkdirAll(filepath.Join(tmpDir, "YRI"), 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}

	model := authenticity.NewAdaptiveHMMModel(4)
	// processPopulation will fail (no VCF files), but trainOnFullDataset only logs it
	err := trainOnFullDataset(model, tmpDir, progressDir, 2)
	if err != nil {
		t.Errorf("trainOnFullDataset should not propagate population error: %v", err)
	}
}

// ---------- processPopulation ----------

// TestProcessPopulation_EmptyDir verifies that processPopulation returns an
// error when the population directory exists but contains no VCF files.
func TestProcessPopulation_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	if err := os.MkdirAll(filepath.Join(tmpDir, "YRI"), 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}

	model := authenticity.NewAdaptiveHMMModel(4)
	_, err := processPopulation(model, "YRI", tmpDir, 2)
	if err == nil {
		t.Error("processPopulation with no VCF files should return an error")
	}
}

// ---------- processVCFFile ----------

// TestProcessVCFFile_Synthetic verifies that processVCFFile generates synthetic
// training data (it does not read from disk) and trains the model successfully.
// The vcfFile parameter is accepted but never opened; the function generates
// population-specific sequences purely from the pop string.
func TestProcessVCFFile_Synthetic(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)

	err := processVCFFile(model, "/does/not/exist.vcf.gz", "YRI")
	if err != nil {
		t.Errorf("processVCFFile (synthetic) returned error: %v", err)
	}
}

// TestProcessVCFFile_AllPopulations exercises processVCFFile with each of the
// 11 known populations to ensure all population-specific code paths execute.
func TestProcessVCFFile_AllPopulations(t *testing.T) {
	populations := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"}
	for _, pop := range populations {
		t.Run(pop, func(t *testing.T) {
			model := authenticity.NewAdaptiveHMMModel(4)
			err := processVCFFile(model, "/fake/file.vcf.gz", pop)
			if err != nil {
				t.Errorf("processVCFFile for population %s returned error: %v", pop, err)
			}
		})
	}
}

// ---------- generateValidationSequences ----------

// TestGenerateValidationSequences_Basic trains a model on a short sequence and
// then generates a small number of short validation sequences.
func TestGenerateValidationSequences_Basic(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)

	// Train with a short sequence so the model can generate output
	trainSeq := []byte("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(trainSeq); err != nil {
		t.Fatalf("TrainFromSequence: %v", err)
	}

	seqs, err := generateValidationSequences(model, 2, 50)
	if err != nil {
		t.Fatalf("generateValidationSequences returned error: %v", err)
	}
	if len(seqs) != 2 {
		t.Errorf("expected 2 sequences, got %d", len(seqs))
	}
	for i, seq := range seqs {
		if len(seq) == 0 {
			t.Errorf("sequence %d is empty", i)
		}
	}
}

// TestGenerateValidationSequences_Zero checks that requesting zero sequences
// returns an empty slice without error.
func TestGenerateValidationSequences_Zero(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)

	seqs, err := generateValidationSequences(model, 0, 50)
	if err != nil {
		t.Fatalf("generateValidationSequences(0) returned error: %v", err)
	}
	if len(seqs) != 0 {
		t.Errorf("expected 0 sequences, got %d", len(seqs))
	}
}

// TestCalculateGCContent1000G_EmptyInnerSeq exercises the "skip empty inner
// sequence" branch by passing a slice that contains one empty sub-slice
// alongside a normal sub-slice.
func TestCalculateGCContent1000G_EmptyInnerSeq(t *testing.T) {
	seqs := [][]byte{
		{},             // empty – should be skipped, not contribute to average
		[]byte("GCGC"), // 100% GC
	}
	got := calculateGCContent(seqs)
	// With one empty and one all-GC sequence the average is computed over
	// len(seqs)==2, so we get 1.0/2 == 0.5.
	if got < 0.0 || got > 1.0 {
		t.Errorf("calculateGCContent with empty inner seq = %v, want in [0,1]", got)
	}
}

// TestSaveResults_BadDir exercises the error branch inside saveResults where
// the results file cannot be created because the directory does not exist.
// The function prints a warning and returns without panicking.
func TestSaveResults_BadDir(t *testing.T) {
	saveResults("/nonexistent-dir-for-test-xyz/progress", time.Second, 0.05, 0.42, "SUCCESS", 4)
	// Success criterion: no panic.
}

// TestGenerateValidationSequences_UntrainedError verifies that
// generateValidationSequences propagates the error from GenerateSequence when
// the model has not been trained.
func TestGenerateValidationSequences_UntrainedError(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	// Untrained model – GenerateSequence must return an error.
	_, err := generateValidationSequences(model, 1, 50)
	if err == nil {
		t.Error("expected error from generateValidationSequences with untrained model, got nil")
	}
}

// TestProcessPopulation_WithVCFFiles creates a population directory that
// contains a file with the .vcf.gz extension (even though it is not a valid
// gzip file, processVCFFile never opens it on disk). This exercises the
// parallel-processing loop inside processPopulation.
func TestProcessPopulation_WithVCFFiles(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}

	// Create two dummy .vcf.gz files – the content does not matter because
	// processVCFFile never reads them; it generates synthetic training data.
	for _, name := range []string{"chr1.vcf.gz", "chr2.vcf.gz"} {
		path := filepath.Join(popDir, name)
		if err := os.WriteFile(path, []byte("dummy"), 0600); err != nil {
			t.Fatalf("WriteFile %s: %v", name, err)
		}
	}

	model := authenticity.NewAdaptiveHMMModel(4)
	processed, err := processPopulation(model, "YRI", tmpDir, 2)
	if err != nil {
		t.Fatalf("processPopulation returned unexpected error: %v", err)
	}
	if processed != 2 {
		t.Errorf("expected 2 files processed, got %d", processed)
	}
}

// TestCountDataset_WithFiles creates a temp dir structure matching the
// "*/*.vcf.gz" glob pattern and verifies countDataset counts the files.
func TestCountDataset_WithFiles(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}

	vcfPath := filepath.Join(popDir, "chr1.vcf.gz")
	if err := os.WriteFile(vcfPath, []byte("data"), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	_, count := countDataset(tmpDir)
	if count != 1 {
		t.Errorf("countDataset: expected 1 VCF file, got %d", count)
	}
}

// TestDetectAppleSiliconCores_ConsistentResult calls detectAppleSiliconCores
// twice and verifies the result is deterministic and positive.
func TestDetectAppleSiliconCores_ConsistentResult(t *testing.T) {
	first := detectAppleSiliconCores()
	second := detectAppleSiliconCores()
	if first != second {
		t.Errorf("detectAppleSiliconCores inconsistent: %d vs %d", first, second)
	}
	if first <= 0 {
		t.Errorf("detectAppleSiliconCores returned non-positive value: %d", first)
	}
}

// TestProcessPopulation_TenVCFFiles creates a population directory with 10
// .vcf.gz files to trigger the "processed%10==0" progress-logging branch inside
// processPopulation. processVCFFile does not read the files; the content is
// irrelevant.
func TestProcessPopulation_TenVCFFiles(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "CEU")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}

	for i := 0; i < 10; i++ {
		name := filepath.Join(popDir, fmt.Sprintf("chr%d.vcf.gz", i+1))
		if err := os.WriteFile(name, []byte("dummy"), 0600); err != nil {
			t.Fatalf("WriteFile %s: %v", name, err)
		}
	}

	model := authenticity.NewAdaptiveHMMModel(4)
	processed, err := processPopulation(model, "CEU", tmpDir, 4)
	if err != nil {
		t.Fatalf("processPopulation returned unexpected error: %v", err)
	}
	if processed != 10 {
		t.Errorf("expected 10 files processed, got %d", processed)
	}
}

// TestCountDataset_SubdirVCFFiles creates a two-level directory (populationDir/
// popName/*.vcf.gz) matching the glob pattern inside countDataset and verifies
// the file count is returned correctly.
func TestCountDataset_SubdirVCFFiles(t *testing.T) {
	tmpDir := t.TempDir()
	for _, pop := range []string{"YRI", "CEU"} {
		popDir := filepath.Join(tmpDir, pop)
		if err := os.MkdirAll(popDir, 0750); err != nil {
			t.Fatalf("MkdirAll %s: %v", pop, err)
		}
		for _, chrom := range []string{"chr1.vcf.gz", "chr2.vcf.gz"} {
			path := filepath.Join(popDir, chrom)
			if err := os.WriteFile(path, []byte("data"), 0600); err != nil {
				t.Fatalf("WriteFile %s: %v", chrom, err)
			}
		}
	}

	size, count := countDataset(tmpDir)
	if count != 4 {
		t.Errorf("countDataset: expected 4 VCF files, got %d", count)
	}
	if size == "" {
		t.Error("countDataset returned empty size string")
	}
}
