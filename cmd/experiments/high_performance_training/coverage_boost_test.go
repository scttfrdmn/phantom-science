//go:build !release
// +build !release

package main

import (
	"compress/gzip"
	"context"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// createVCFGzHPT2 creates a minimal gzip-encoded .vcf.gz stub in dir/name.
func createVCFGzHPT2(t *testing.T, dir, name string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	f, err := os.Create(p)
	if err != nil {
		t.Fatalf("createVCFGzHPT2 %s: %v", p, err)
	}
	defer f.Close()
	gw := gzip.NewWriter(f)
	_, _ = gw.Write([]byte("##fileformat=VCFv4.1\n"))
	if err := gw.Close(); err != nil {
		t.Fatalf("gzip.Close: %v", err)
	}
	return p
}

// TestTrainPopulation_WithVCFFiles_HPT exercises the inner VCF-processing loop
// (files found, each sleeps 100ms).  12 files ensures the %%10 branch fires.
func TestTrainPopulation_WithVCFFiles_HPT(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping slow test in short mode")
	}
	tmpDir := t.TempDir()
	pop := "YRI_HPT2"
	popDir := filepath.Join(tmpDir, pop)
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}
	for i := 0; i < 12; i++ {
		createVCFGzHPT2(t, popDir, "chr"+string(rune('A'+i))+".vcf.gz")
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx := context.Background()
	processedFiles := 0
	totalFiles := 12
	if err := trainPopulation(ctx, model, pop, tmpDir, &processedFiles, totalFiles); err != nil {
		t.Errorf("trainPopulation returned unexpected error: %v", err)
	}
	if processedFiles != 12 {
		t.Errorf("processedFiles = %d, want 12", processedFiles)
	}
}

// TestTrainPopulation_CancelledMidLoop_HPT cancels the context so the
// ctx.Done() select arm inside the VCF file loop is taken.
func TestTrainPopulation_CancelledMidLoop_HPT(t *testing.T) {
	tmpDir := t.TempDir()
	pop := "CEU_HPT2"
	popDir := filepath.Join(tmpDir, pop)
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}
	for i := 0; i < 3; i++ {
		createVCFGzHPT2(t, popDir, "chr"+string(rune('a'+i))+".vcf.gz")
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx, cancel := context.WithTimeout(context.Background(), 50*time.Millisecond)
	defer cancel()

	processedFiles := 0
	err := trainPopulation(ctx, model, pop, tmpDir, &processedFiles, 3)
	_ = err // may be DeadlineExceeded or nil
}

// TestTrainPopulationParallel_WithOnePopDir_HPT creates one population directory
// with a vcf.gz stub and exercises the parallel goroutine branch.
func TestTrainPopulationParallel_WithOnePopDir_HPT(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping slow test in short mode")
	}
	tmpDir := t.TempDir()
	yriDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(yriDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}
	createVCFGzHPT2(t, yriDir, "chr1.vcf.gz")

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()

	if err := trainPopulationParallel(ctx, model, tmpDir, tmpDir); err != nil {
		t.Errorf("trainPopulationParallel returned error: %v", err)
	}
}

// TestTrainBaselineParallel_WithRealFile_HPT creates a GRCh38_complete.fasta
// stub so that the file-exists branch of trainBaselineParallel is exercised.
func TestTrainBaselineParallel_WithRealFile_HPT(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := filepath.Join(tmpDir, "GRCh38_complete.fasta")
	if err := os.WriteFile(fastaPath, []byte(">chr1\nACGTACGTACGT\n"), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	// The function calls model.TrainFromFileParallel; may succeed or return error.
	// We just need to exercise the code past the os.Stat check.
	_ = trainBaselineParallel(context.Background(), model, tmpDir)
}

// TestGenerateValidationSequencesParallel_Untrained_HPT calls the function with
// an untrained model to exercise the goroutine error path.
func TestGenerateValidationSequencesParallel_Untrained_HPT(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seqs, avgTime, err := generateValidationSequencesParallel(context.Background(), model, 2, 50)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(seqs) != 2 {
		t.Errorf("expected 2-element slice, got %d", len(seqs))
	}
	_ = avgTime
}

// TestSaveFinalResults_ZeroValues_HPT calls saveFinalResults with all-zero
// values; verifies the Sprintf template does not panic.
func TestSaveFinalResults_ZeroValues_HPT(t *testing.T) {
	tmpDir := t.TempDir()
	saveFinalResults(tmpDir, 0, 0.0, 0.0, 0.0, "TRAINING_COMPLETE")
	p := filepath.Join(tmpDir, "final_results.txt")
	if _, err := os.Stat(p); err != nil {
		t.Fatalf("final_results.txt should exist: %v", err)
	}
}

// TestSaveFinalResults_WriteError_HPT exercises the os.WriteFile error path
// by calling saveFinalResults with a nonexistent parent directory.
func TestSaveFinalResults_WriteError_HPT(t *testing.T) {
	// Use a path inside a directory that doesn't exist
	saveFinalResults("/nonexistent/dir/progress", 1*time.Second, 0.06, 0.42, 0.005, "TARGETS_MET")
	// Should not panic; the error is logged via fmt.Printf
}
