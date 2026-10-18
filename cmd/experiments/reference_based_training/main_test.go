//go:build !release
// +build !release

package main

import (
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestExtractGenomicWindow(t *testing.T) {
	refSeq := []byte("AAAAACCCCCGGGGGTTTTTNNNNN")

	t.Run("normal window", func(t *testing.T) {
		win := extractGenomicWindow(refSeq, 12, 10)
		if len(win) == 0 {
			t.Error("expected non-empty window")
		}
	})

	t.Run("pos out of range returns nil", func(t *testing.T) {
		win := extractGenomicWindow(refSeq, -1, 10)
		if win != nil {
			t.Error("expected nil for negative pos")
		}
	})

	t.Run("pos at end of sequence returns nil", func(t *testing.T) {
		win := extractGenomicWindow(refSeq, len(refSeq), 10)
		if win != nil {
			t.Error("expected nil for pos >= len(refSeq)")
		}
	})

	t.Run("window clipped at start", func(t *testing.T) {
		win := extractGenomicWindow(refSeq, 1, 10)
		if win == nil {
			t.Error("expected non-nil window near start")
		}
	})

	t.Run("window clipped at end", func(t *testing.T) {
		win := extractGenomicWindow(refSeq, len(refSeq)-2, 10)
		if win == nil {
			t.Error("expected non-nil window near end")
		}
	})
}

func TestApplyVariantToWindow(t *testing.T) {
	window := []byte("AAAAAAAAAA")

	t.Run("SNP substitution", func(t *testing.T) {
		result := applyVariantToWindow(window, 5, "A", "G")
		if result == nil {
			t.Fatal("expected non-nil result for SNP")
		}
		if len(result) != len(window) {
			t.Errorf("SNP should preserve length: got %d, want %d", len(result), len(window))
		}
		if result[5] != 'G' {
			t.Errorf("expected G at pos 5, got %c", result[5])
		}
	})

	t.Run("SNP does not modify original", func(t *testing.T) {
		orig := make([]byte, len(window))
		copy(orig, window)
		applyVariantToWindow(window, 3, "A", "C")
		for i, b := range window {
			if b != orig[i] {
				t.Errorf("original window was modified at pos %d", i)
			}
		}
	})

	t.Run("insertion", func(t *testing.T) {
		result := applyVariantToWindow(window, 3, "A", "ATG")
		if result == nil {
			t.Fatal("expected non-nil result for insertion")
		}
		// Insertion makes sequence longer
		if len(result) <= len(window) {
			t.Errorf("insertion should lengthen sequence, got %d <= %d", len(result), len(window))
		}
	})

	t.Run("deletion", func(t *testing.T) {
		result := applyVariantToWindow(window, 2, "AAA", "A")
		if result == nil {
			t.Fatal("expected non-nil result for deletion")
		}
		// Deletion makes sequence shorter
		if len(result) >= len(window) {
			t.Errorf("deletion should shorten sequence, got %d >= %d", len(result), len(window))
		}
	})

	t.Run("negative relative pos returns nil", func(t *testing.T) {
		result := applyVariantToWindow(window, -1, "A", "G")
		if result != nil {
			t.Error("expected nil for negative relativePos")
		}
	})

	t.Run("relative pos >= window len returns nil", func(t *testing.T) {
		result := applyVariantToWindow(window, len(window), "A", "G")
		if result != nil {
			t.Error("expected nil for relativePos >= len(window)")
		}
	})
}

func TestCalculateComplexityEmpty(t *testing.T) {
	if got := calculateComplexity(nil, 6); got != 0 {
		t.Errorf("expected 0 for nil, got %f", got)
	}
	if got := calculateComplexity([][]byte{}, 6); got != 0 {
		t.Errorf("expected 0 for empty slice, got %f", got)
	}
}

func TestCalculateComplexityBounded(t *testing.T) {
	seq := []byte("ACGTACGTGCATGCAT")
	got := calculateComplexity([][]byte{seq}, 4)
	if got < 0 || got > 1 {
		t.Errorf("complexity should be in [0,1], got %f", got)
	}
}

func TestCalculateGCContentEmpty(t *testing.T) {
	if got := calculateGCContent(nil); got != 0 {
		t.Errorf("expected 0 for nil, got %f", got)
	}
	if got := calculateGCContent([][]byte{}); got != 0 {
		t.Errorf("expected 0 for empty slice, got %f", got)
	}
}

func TestCalculateGCContentAllGC(t *testing.T) {
	seq := []byte("GCGCGCGCGC")
	if got := calculateGCContent([][]byte{seq}); got != 1.0 {
		t.Errorf("expected 1.0, got %f", got)
	}
}

func TestCalculateGCContentHalf(t *testing.T) {
	seq := []byte("AACCGGTT")
	got := calculateGCContent([][]byte{seq})
	if got != 0.5 {
		t.Errorf("expected 0.5, got %f", got)
	}
}

func TestFormatBytes(t *testing.T) {
	tests := []struct {
		input    int64
		expected string
	}{
		{0, "0 B"},
		{512, "512 B"},
		{1023, "1023 B"},
		{1024, "1.0 KB"},
		{1536, "1.5 KB"},
		{1024 * 1024, "1.0 MB"},
		{1024 * 1024 * 1024, "1.0 GB"},
	}
	for _, tt := range tests {
		got := formatBytes(tt.input)
		if got != tt.expected {
			t.Errorf("formatBytes(%d) = %q, want %q", tt.input, got, tt.expected)
		}
	}
}

func TestAbs(t *testing.T) {
	tests := []struct{ in, out float64 }{
		{0, 0},
		{1.5, 1.5},
		{-1.5, 1.5},
		{-0.001, 0.001},
	}
	for _, tt := range tests {
		if got := abs(tt.in); got != tt.out {
			t.Errorf("abs(%f) = %f, want %f", tt.in, got, tt.out)
		}
	}
}

func TestDetectAppleSiliconCoresPositive(t *testing.T) {
	if cores := detectAppleSiliconCores(); cores <= 0 {
		t.Errorf("expected positive core count, got %d", cores)
	}
}

func TestGetSuccessEmoji(t *testing.T) {
	// complexity >= 0.08 AND gc in [0.40,0.44]
	if got := getSuccessEmoji(0.09, 0.42); got != "✅" {
		t.Errorf("expected ✅, got %q", got)
	}
	// gc out of range but better than Phase 2.2
	if got := getSuccessEmoji(0.09, 0.46); got == "✅" {
		t.Errorf("expected non-✅ for gc out of optimal range, got %q", got)
	}
}

func TestGetNextSteps(t *testing.T) {
	// Both thresholds met
	result := getNextSteps(0.09, 0.42)
	if result == "" {
		t.Error("expected non-empty steps string")
	}
	// Partial success
	result2 := getNextSteps(0.09, 0.48)
	if result2 == "" {
		t.Error("expected non-empty steps string for partial success")
	}
	// Nothing meets threshold
	result3 := getNextSteps(0.01, 0.60)
	if result3 == "" {
		t.Error("expected non-empty steps string for failure case")
	}
}

// ---------- saveReferenceBasedResults ----------

func TestSaveReferenceBasedResults_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	saveReferenceBasedResults(tmpDir, 30*time.Second, 0.082, 0.42, 1000000, 5000, 200, 6, "SUCCESS", 1.5, 1000)
	path := filepath.Join(tmpDir, "reference_based_results.txt")
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty results file")
	}
}

// ---------- loadReferenceGenome ----------

// TestLoadReferenceGenome_NonexistentFile verifies that loadReferenceGenome
// returns an error when the given file does not exist.
func TestLoadReferenceGenome_NonexistentFile(t *testing.T) {
	_, err := loadReferenceGenome("/nonexistent/genome.fasta", []string{"chr1"})
	if err == nil {
		t.Error("expected error for nonexistent reference file, got nil")
	}
}

// TestLoadReferenceGenome_EmptyDir verifies that loadReferenceGenome returns
// an error (file not found) when called with a path in an empty tmpDir.
func TestLoadReferenceGenome_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	_, err := loadReferenceGenome(filepath.Join(tmpDir, "genome.fasta"), []string{"chr22"})
	if err == nil {
		t.Error("expected error for missing FASTA file, got nil")
	}
}

// TestLoadReferenceGenome_ValidFASTA verifies that loadReferenceGenome can
// parse a minimal FASTA file and load the requested chromosome.
func TestLoadReferenceGenome_ValidFASTA(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := filepath.Join(tmpDir, "genome.fasta")
	fastaContent := ">chr22\nACGTACGTACGT\nGCATGCATGCAT\n>chr1\nAAAAAAAA\n"
	if err := os.WriteFile(fastaPath, []byte(fastaContent), 0600); err != nil {
		t.Fatalf("failed to write FASTA file: %v", err)
	}

	refGenome, err := loadReferenceGenome(fastaPath, []string{"chr22"})
	if err != nil {
		t.Fatalf("loadReferenceGenome failed: %v", err)
	}
	if _, ok := refGenome.chromosomes["chr22"]; !ok {
		t.Error("expected chr22 to be loaded, but it is missing")
	}
	if _, ok := refGenome.chromosomes["chr1"]; ok {
		t.Error("expected chr1 to be absent (not in requested list), but it was loaded")
	}
}

// ---------- processReferenceBasedPopulations ----------

// TestProcessReferenceBasedPopulations_EmptyDir verifies that
// processReferenceBasedPopulations returns (0, 0, 0) when no population
// subdirectories exist in the given directory.
func TestProcessReferenceBasedPopulations_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := &ReferenceGenome{chromosomes: make(map[string][]byte)}
	populations := []string{"YRI", "CEU"}
	b, v, w := processReferenceBasedPopulations(model, refGenome, tmpDir, populations, 1000, 1)
	if b != 0 || v != 0 || w != 0 {
		t.Errorf("processReferenceBasedPopulations(empty dir) = (%d,%d,%d), want (0,0,0)", b, v, w)
	}
}

// ---------- processReferenceBasedPopulationVCFs ----------

// TestProcessReferenceBasedPopulationVCFs_EmptyDir verifies that
// processReferenceBasedPopulationVCFs returns (0, 0, 0) when the population
// directory exists but contains no VCF files.
func TestProcessReferenceBasedPopulationVCFs_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	// Create the population subdirectory without any VCF files.
	popDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("failed to create pop dir: %v", err)
	}

	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := &ReferenceGenome{chromosomes: make(map[string][]byte)}
	b, v, w := processReferenceBasedPopulationVCFs(model, refGenome, "YRI", tmpDir, 1000, 1)
	if b != 0 || v != 0 || w != 0 {
		t.Errorf("processReferenceBasedPopulationVCFs(no vcf files) = (%d,%d,%d), want (0,0,0)", b, v, w)
	}
}

// ---------- processReferenceBasedVCFFile ----------

// TestProcessReferenceBasedVCFFile_NonexistentFile verifies that
// processReferenceBasedVCFFile returns (0, 0, 0) when the file does not exist.
func TestProcessReferenceBasedVCFFile_NonexistentFile(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := &ReferenceGenome{chromosomes: make(map[string][]byte)}
	b, v, w := processReferenceBasedVCFFile(model, refGenome, "/nonexistent/chr22.vcf.gz", "YRI", 1000)
	if b != 0 || v != 0 || w != 0 {
		t.Errorf("processReferenceBasedVCFFile(nonexistent) = (%d,%d,%d), want (0,0,0)", b, v, w)
	}
}

func TestProcessReferenceBasedVCFFile_NotGzip(t *testing.T) {
	// A plain-text file is not a valid gzip; gzip.NewReader fails → returns (0, 0, 0)
	tmpDir := t.TempDir()
	vcfPath := filepath.Join(tmpDir, "chr21.vcf.gz")
	if err := os.WriteFile(vcfPath, []byte("not a gzip file"), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := &ReferenceGenome{chromosomes: make(map[string][]byte)}
	b, v, w := processReferenceBasedVCFFile(model, refGenome, vcfPath, "YRI", 1000)
	// The function returns (0, 0, 0) when gzip.NewReader fails, despite the file being accessible
	if v != 0 || w != 0 {
		t.Errorf("expected 0 variants/windows for non-gzip file, got v=%d w=%d", v, w)
	}
	// b may be 0 or fileSize depending on implementation; just verify no panic
	_ = b
}

// ---------- generateValidationSequences (reference_based_training) ----------

// TestRefBasedGenerateValidationSequences verifies that
// generateValidationSequences returns the requested number of sequences after
// the model has been trained.
func TestRefBasedGenerateValidationSequences(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	trainingSeq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	if err := model.TrainFromSequence(trainingSeq); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	seqs, err := generateValidationSequences(model, 2, 50)
	if err != nil {
		t.Fatalf("generateValidationSequences returned error: %v", err)
	}
	if len(seqs) != 2 {
		t.Errorf("expected 2 sequences, got %d", len(seqs))
	}
	for i, s := range seqs {
		if len(s) == 0 {
			t.Errorf("sequence %d is empty", i)
		}
	}
}

// TestRefBasedGenerateValidationSequences_UntrainedModel verifies that
// generateValidationSequences returns an error when the model is untrained.
func TestRefBasedGenerateValidationSequences_UntrainedModel(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	_, err := generateValidationSequences(model, 1, 50)
	if err == nil {
		t.Error("expected error for untrained model, got nil")
	}
}
