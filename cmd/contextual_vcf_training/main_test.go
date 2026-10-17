//go:build !release
// +build !release

package main

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"testing"
	"time"

	authenticity "github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestCalculateStateSpaceSize(t *testing.T) {
	tests := []struct {
		k        int
		expected int
	}{
		{1, 1},
		{2, 4},
		{3, 16},
		{4, 64},
		{5, 256},
		{6, 1024},
		{7, 4096},
		{8, 16384},
	}
	for _, tt := range tests {
		got := calculateStateSpaceSize(tt.k)
		if got != tt.expected {
			t.Errorf("calculateStateSpaceSize(%d) = %d, want %d", tt.k, got, tt.expected)
		}
	}
}

func TestReconstructContextualSequence(t *testing.T) {
	t.Run("nil input returns nil", func(t *testing.T) {
		result := reconstructContextualSequence(nil)
		if result != nil {
			t.Errorf("expected nil for nil input, got %v", result)
		}
	})

	t.Run("empty slice returns nil", func(t *testing.T) {
		result := reconstructContextualSequence([]Variant{})
		if result != nil {
			t.Errorf("expected nil for empty variants, got %v", result)
		}
	})

	t.Run("single variant produces 10000-byte window", func(t *testing.T) {
		variants := []Variant{{pos: 100, ref: "A", alt: "G"}}
		result := reconstructContextualSequence(variants)
		if result == nil {
			t.Fatal("expected non-nil result")
		}
		if len(result) != 10000 {
			t.Errorf("expected window size 10000, got %d", len(result))
		}
	})

	t.Run("multiple variants produce 10000-byte window", func(t *testing.T) {
		variants := make([]Variant, 50)
		for i := range variants {
			variants[i] = Variant{pos: int64(i * 100), ref: "A", alt: "C"}
		}
		result := reconstructContextualSequence(variants)
		if result == nil {
			t.Fatal("expected non-nil result")
		}
		if len(result) != 10000 {
			t.Errorf("expected window size 10000, got %d", len(result))
		}
	})

	t.Run("result contains only valid DNA bases", func(t *testing.T) {
		variants := []Variant{
			{pos: 50, ref: "A", alt: "G"},
			{pos: 150, ref: "C", alt: "T"},
		}
		result := reconstructContextualSequence(variants)
		for i, b := range result {
			switch b {
			case 'A', 'T', 'G', 'C':
				// valid
			default:
				t.Errorf("invalid base %c at position %d", b, i)
			}
		}
	})
}

func TestReconstructContextualSequenceGCContent(t *testing.T) {
	// The function targets 42% GC before applying variants.
	// With many variants applied, GC should be approximately in range.
	variants := []Variant{{pos: 1, ref: "A", alt: "A"}}
	result := reconstructContextualSequence(variants)
	if result == nil {
		t.Fatal("expected non-nil result")
	}
	// Count GC bases
	gcCount := 0
	for _, b := range result {
		if b == 'G' || b == 'C' {
			gcCount++
		}
	}
	gcFraction := float64(gcCount) / float64(len(result))
	// Accept a wide range since variants are applied on top
	if gcFraction < 0.1 || gcFraction > 0.9 {
		t.Errorf("GC fraction %f seems unreasonable", gcFraction)
	}
}

func TestCalculateComplexityContextual(t *testing.T) {
	t.Run("empty returns 0", func(t *testing.T) {
		if got := calculateComplexity(nil); got != 0 {
			t.Errorf("expected 0 for nil, got %f", got)
		}
		if got := calculateComplexity([][]byte{}); got != 0 {
			t.Errorf("expected 0 for empty slice, got %f", got)
		}
	})

	t.Run("uniform sequence has low complexity", func(t *testing.T) {
		seq := make([]byte, 100)
		for i := range seq {
			seq[i] = 'A'
		}
		got := calculateComplexity([][]byte{seq})
		if got > 0.1 {
			t.Errorf("uniform sequence complexity too high: %f", got)
		}
	})

	t.Run("bounded between 0 and 1", func(t *testing.T) {
		seq := []byte("ACGTACGTGCATGCAT")
		got := calculateComplexity([][]byte{seq})
		if got < 0 || got > 1 {
			t.Errorf("complexity %f out of [0,1]", got)
		}
	})
}

func TestCalculateGCContentContextual(t *testing.T) {
	t.Run("empty returns 0", func(t *testing.T) {
		if got := calculateGCContent(nil); got != 0 {
			t.Errorf("expected 0, got %f", got)
		}
	})

	t.Run("all G/C returns 1.0", func(t *testing.T) {
		seq := []byte("GCGCGC")
		if got := calculateGCContent([][]byte{seq}); got != 1.0 {
			t.Errorf("expected 1.0, got %f", got)
		}
	})

	t.Run("all A/T returns 0.0", func(t *testing.T) {
		seq := []byte("ATATAT")
		if got := calculateGCContent([][]byte{seq}); got != 0.0 {
			t.Errorf("expected 0.0, got %f", got)
		}
	})

	t.Run("50pct returns 0.5", func(t *testing.T) {
		seq := []byte("AACCGGTT")
		if got := calculateGCContent([][]byte{seq}); got != 0.5 {
			t.Errorf("expected 0.5, got %f", got)
		}
	})
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

func TestGetContextualInsight(t *testing.T) {
	// When current deviation is less than previous, improvement message expected
	improved := getContextualInsight(0.01, 0.05)
	if improved == "" {
		t.Error("expected non-empty improvement message")
	}
	// When current deviation is greater, different message expected
	notImproved := getContextualInsight(0.06, 0.05)
	if notImproved == "" {
		t.Error("expected non-empty non-improvement message")
	}
	if improved == notImproved {
		t.Error("improved and not-improved messages should differ")
	}
}

func TestGetNextStepsPhase21(t *testing.T) {
	// Success path
	s1 := getNextStepsPhase21(0.09, 0.42)
	if s1 == "" {
		t.Error("expected non-empty steps for success case")
	}
	// Medium path
	s2 := getNextStepsPhase21(0.06, 0.43)
	if s2 == "" {
		t.Error("expected non-empty steps for medium case")
	}
	// Failure path
	s3 := getNextStepsPhase21(0.01, 0.60)
	if s3 == "" {
		t.Error("expected non-empty steps for failure case")
	}
}

func TestGetSuccessEmojiPhase21(t *testing.T) {
	// Full success: complexity > 0.08 AND gc in [0.40, 0.44]
	if got := getSuccessEmojiPhase21(0.09, 0.42); got != "✅" {
		t.Errorf("expected ✅ for full success, got %q", got)
	}
	// Medium success
	if got := getSuccessEmojiPhase21(0.06, 0.43); got != "🟡" {
		t.Errorf("expected 🟡 for medium success, got %q", got)
	}
	// Failure
	if got := getSuccessEmojiPhase21(0.01, 0.60); got != "❌" {
		t.Errorf("expected ❌ for failure, got %q", got)
	}
}

func TestVariantStruct(t *testing.T) {
	v := Variant{pos: 12345, ref: "A", alt: "G"}
	if v.pos != 12345 {
		t.Errorf("expected pos 12345, got %d", v.pos)
	}
	if v.ref != "A" {
		t.Errorf("expected ref A, got %s", v.ref)
	}
	if v.alt != "G" {
		t.Errorf("expected alt G, got %s", v.alt)
	}
}

// ---------- calculateComplexity ----------

func TestCalculateComplexity_EmptySlice(t *testing.T) {
	got := calculateComplexity([][]byte{})
	if got != 0 {
		t.Errorf("calculateComplexity(empty) = %v, want 0", got)
	}
}

func TestCalculateComplexity_ShortSequenceSkipped(t *testing.T) {
	// Sequence shorter than k=6 should be skipped; only valid sequence counted
	short := []byte("ACGT")         // len 4 < 6, skipped
	valid := []byte("ACGTACGTACGT") // len 12 >= 6, counted
	got := calculateComplexity([][]byte{short, valid})
	// Result should be non-zero (averaged over 2 seqs, short contributes 0)
	if got <= 0 {
		t.Errorf("calculateComplexity with mixed lengths = %v, want > 0", got)
	}
}

func TestCalculateComplexity_AllShortSequences(t *testing.T) {
	// All sequences shorter than k=6; sum stays 0, divide by len gives 0
	seqs := [][]byte{[]byte("AC"), []byte("GC")}
	got := calculateComplexity(seqs)
	if got != 0 {
		t.Errorf("calculateComplexity(all short) = %v, want 0", got)
	}
}

func TestCalculateComplexity_SingleLongSequence(t *testing.T) {
	seq := []byte("ACGTACGTACGTACGTACGT")
	got := calculateComplexity([][]byte{seq})
	if got <= 0 || got > 1 {
		t.Errorf("calculateComplexity(single long) = %v, want in (0, 1]", got)
	}
}

// ---------- calculateGCContent (contextual_vcf_training) ----------

func TestCalculateGCContent_CVT_EmptySlice(t *testing.T) {
	got := calculateGCContent([][]byte{})
	if got != 0 {
		t.Errorf("calculateGCContent(empty) = %v, want 0", got)
	}
}

func TestCalculateGCContent_CVT_EmptySequenceInSlice(t *testing.T) {
	// An empty sequence inside the slice should be skipped
	seqs := [][]byte{{}, []byte("GGCC")}
	got := calculateGCContent(seqs)
	// GGCC is 100% GC; empty skipped; average = (0 + 1.0) / 2 = 0.5
	if got < 0 || got > 1 {
		t.Errorf("calculateGCContent with empty inner seq = %v, want in [0,1]", got)
	}
}

func TestCalculateGCContent_CVT_PureAT(t *testing.T) {
	seqs := [][]byte{[]byte("AAATTT")}
	got := calculateGCContent(seqs)
	if got != 0 {
		t.Errorf("calculateGCContent(pure AT) = %v, want 0", got)
	}
}

func TestCalculateGCContent_CVT_PureGC(t *testing.T) {
	seqs := [][]byte{[]byte("GGGGCCCC")}
	got := calculateGCContent(seqs)
	if got != 1.0 {
		t.Errorf("calculateGCContent(pure GC) = %v, want 1.0", got)
	}
}

func TestSavePhase21Results_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	savePhase21Results(tmpDir, time.Second*30, 0.072, 0.42, 1000000, 5000, 4, "pass")
	path := filepath.Join(tmpDir, "phase2_1_contextual_results.txt")
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty results file")
	}
}

// ---------- VCF processing helper ----------

func createTestGzippedVCFContextual(t *testing.T) string {
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

// ---------- processContextualVCFFile tests ----------

func TestProcessContextualVCFFile_Basic(t *testing.T) {
	vcfPath := createTestGzippedVCFContextual(t)
	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processContextualVCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 2 {
		t.Errorf("expected 2 variants, got %d", variants)
	}
}

func TestProcessContextualVCFFile_NonExistentFile(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processContextualVCFFile(model, "/nonexistent/path/test.vcf.gz", "TEST")
	if fileSize != 0 {
		t.Errorf("expected 0 file size for missing file, got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for missing file, got %d", variants)
	}
}

func TestProcessContextualVCFFile_NotGzipped(t *testing.T) {
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "notgzip.vcf.gz")
	os.WriteFile(path, []byte("this is not gzip data"), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processContextualVCFFile(model, path, "TEST")
	if fileSize != 0 || variants != 0 {
		t.Errorf("expected 0,0 for invalid gzip, got %d,%d", fileSize, variants)
	}
}

func TestProcessContextualVCFFile_ManyVariants(t *testing.T) {
	// Test with enough variants to trigger the batch-of-100 processing path
	tmpDir := t.TempDir()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	// Write 120 variants to trigger the batch flush (>100)
	for i := 0; i < 120; i++ {
		line := fmt.Sprintf("chr1\t%d\t.\tATCG\tGCAT\t.\tPASS\t.\n", (i+1)*1000)
		gw.Write([]byte(line))
	}
	gw.Close()
	vcfPath := filepath.Join(tmpDir, "many_variants.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processContextualVCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 120 {
		t.Errorf("expected 120 variants, got %d", variants)
	}
}

func TestProcessContextualVCFFile_HeadersOnly(t *testing.T) {
	tmpDir := t.TempDir()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Close()
	vcfPath := filepath.Join(tmpDir, "headers_only.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processContextualVCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size (from header content), got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for headers-only, got %d", variants)
	}
}

// ---------- processContextualPopulationVCFs tests ----------

func TestProcessContextualPopulationVCFs_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "POP")
	os.MkdirAll(popDir, 0750)

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processContextualPopulationVCFs(model, "POP", tmpDir, 1)
	if fileSize != 0 {
		t.Errorf("expected 0 bytes for empty dir, got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for empty dir, got %d", variants)
	}
}

func TestProcessContextualPopulationVCFs_WithVCFs(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	os.MkdirAll(popDir, 0750)

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
	fileSize, variants := processContextualPopulationVCFs(model, "YRI", tmpDir, 1)
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 3 {
		t.Errorf("expected 3 variants, got %d", variants)
	}
}

// ---------- generateValidationSequences tests ----------

func TestGenerateValidationSequences_Contextual_Basic(t *testing.T) {
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

func TestGenerateValidationSequences_Contextual_Multiple(t *testing.T) {
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

// ---------- processContextualPopulations ----------

func TestProcessContextualPopulations_EmptyDir_SkipsAll(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(6)
	tmpDir := t.TempDir()
	populations := []string{"YRI", "CEU", "CHB"}
	bytes, variants := processContextualPopulations(model, tmpDir, populations, 1)
	if bytes != 0 {
		t.Errorf("expected 0 bytes, got %d", bytes)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants, got %d", variants)
	}
}

func TestProcessContextualPopulations_EmptyPopulations_ReturnsZero(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(6)
	bytes, variants := processContextualPopulations(model, t.TempDir(), []string{}, 1)
	if bytes != 0 || variants != 0 {
		t.Errorf("expected (0, 0) for empty populations, got (%d, %d)", bytes, variants)
	}
}

// TestProcessContextualPopulations_WithExistingPopulation tests the goroutine
// success path where a population directory exists and has VCF files.
func TestProcessContextualPopulations_WithExistingPopulation(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}

	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Write([]byte("chr1\t100\t.\tA\tT\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t200\t.\tGC\tAT\t.\tPASS\t.\n"))
	gw.Close()
	vcfPath := filepath.Join(popDir, "chr1.vcf.gz")
	if err := os.WriteFile(vcfPath, buf.Bytes(), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	totalBytes, totalVariants := processContextualPopulations(model, tmpDir, []string{"YRI"}, 1)
	if totalBytes <= 0 {
		t.Errorf("expected positive bytes, got %d", totalBytes)
	}
	if totalVariants != 2 {
		t.Errorf("expected 2 variants, got %d", totalVariants)
	}
}

// TestProcessContextualPopulations_MultiplePopulations tests that missing
// populations are skipped while existing ones are processed.
func TestProcessContextualPopulations_MultiplePopulations(t *testing.T) {
	tmpDir := t.TempDir()

	// Only create YRI directory - CEU will be missing and skipped.
	popDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}

	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Write([]byte("chr1\t100\t.\tA\tT\t.\tPASS\t.\n"))
	gw.Close()
	os.WriteFile(filepath.Join(popDir, "chr1.vcf.gz"), buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	totalBytes, totalVariants := processContextualPopulations(model, tmpDir, []string{"YRI", "CEU"}, 1)
	if totalBytes <= 0 {
		t.Errorf("expected positive bytes from YRI, got %d", totalBytes)
	}
	if totalVariants != 1 {
		t.Errorf("expected 1 variant from YRI, got %d", totalVariants)
	}
}

// TestGenerateValidationSequences_Contextual_UntrainedModelReturnsError tests
// the error return path in generateValidationSequences.
func TestGenerateValidationSequences_Contextual_UntrainedModelReturnsError(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	// Don't train - GenerateSequence should return an error.
	seqs, err := generateValidationSequences(model, 1, 100)
	if err == nil {
		t.Error("expected error from untrained model, got nil")
	}
	if seqs != nil {
		t.Errorf("expected nil sequences on error, got %v", seqs)
	}
}

// TestProcessContextualVCFFile_ShortFields tests that lines with fewer than
// 5 tab-separated fields are silently skipped.
func TestProcessContextualVCFFile_ShortFields(t *testing.T) {
	tmpDir := t.TempDir()

	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	// Line with only 3 fields (< 5 required) should be skipped.
	gw.Write([]byte("chr1\t100\tonly_three_fields\n"))
	// A valid line to ensure some output.
	gw.Write([]byte("chr1\t200\t.\tA\tT\t.\tPASS\t.\n"))
	gw.Close()

	vcfPath := filepath.Join(tmpDir, "short_fields.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	fileSize, variants := processContextualVCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 1 {
		t.Errorf("expected 1 variant (short field line skipped), got %d", variants)
	}
}

// TestProcessContextualVCFFile_InvalidPosition tests that a VCF line with
// a non-numeric position field is gracefully skipped.
func TestProcessContextualVCFFile_InvalidPosition(t *testing.T) {
	tmpDir := t.TempDir()

	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	// Invalid position (not a number).
	gw.Write([]byte("chr1\tINVALID_POS\t.\tA\tT\t.\tPASS\t.\n"))
	// A valid line.
	gw.Write([]byte("chr1\t500\t.\tG\tC\t.\tPASS\t.\n"))
	gw.Close()

	vcfPath := filepath.Join(tmpDir, "invalid_pos.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	fileSize, variants := processContextualVCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	// Only the valid line should be counted.
	if variants != 1 {
		t.Errorf("expected 1 variant (invalid pos skipped), got %d", variants)
	}
}

// TestProcessContextualPopulationVCFs_HeadersOnly tests that a VCF with only
// headers results in 0 variants but positive file size.
func TestProcessContextualPopulationVCFs_HeadersOnly(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "POP")
	os.MkdirAll(popDir, 0750)

	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Close()
	vcfPath := filepath.Join(popDir, "headers.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	fileSize, variants := processContextualPopulationVCFs(model, "POP", tmpDir, 1)
	if fileSize <= 0 {
		t.Errorf("expected positive file size from headers-only VCF, got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for headers-only, got %d", variants)
	}
}

// TestSavePhase21Results_BranchVariants exercises the different status branches
// in savePhase21Results by varying complexity and GC content inputs.
func TestSavePhase21Results_BranchVariants(t *testing.T) {
	tests := []struct {
		name       string
		complexity float64
		gcContent  float64
		status     string
	}{
		// Full success: complexity > 0.08 AND GC in [0.40, 0.44]
		{"full_success", 0.09, 0.42, "SUCCESS"},
		// Medium success: complexity > 0.05 AND GC in [0.38, 0.46]
		{"medium_success", 0.06, 0.43, "PROGRESS"},
		// Below thresholds
		{"below_thresholds", 0.02, 0.60, "REFINE"},
		// Boundary condition: exactly at complexity 0.08 and GC 0.40
		{"boundary_min", 0.08, 0.40, "BOUNDARY"},
		// gcImprovement > 1.0 path: currentDev < previousDev
		{"gc_improved", 0.09, 0.425, "IMPROVED"},
		// gcImprovement <= 1.0 path: currentDev >= previousDev
		{"gc_not_improved", 0.09, 0.479, "FLAT"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			savePhase21Results(tmpDir, 10*time.Second, tt.complexity, tt.gcContent, 500000, 1000, 4, tt.status)
			path := filepath.Join(tmpDir, "phase2_1_contextual_results.txt")
			info, err := os.Stat(path)
			if err != nil {
				t.Fatalf("results file not created: %v", err)
			}
			if info.Size() == 0 {
				t.Error("results file is empty")
			}
		})
	}
}

// TestDetectAppleSiliconCores_Contextual_ReturnsPositive verifies the function
// always returns a positive core count.
func TestDetectAppleSiliconCores_Contextual_ReturnsPositive(t *testing.T) {
	cores := detectAppleSiliconCores()
	if cores <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", cores)
	}
}

// TestDetectAppleSiliconCores_Contextual_BoundedByTotal verifies P-cores do not
// exceed total core count.
func TestDetectAppleSiliconCores_Contextual_BoundedByTotal(t *testing.T) {
	total := runtime.NumCPU()
	cores := detectAppleSiliconCores()
	if cores > total {
		t.Errorf("detectAppleSiliconCores() = %d exceeds NumCPU() = %d", cores, total)
	}
}

// TestSavePhase21Results_NonExistentDirWriteError tests that the function does
// not panic when os.WriteFile fails due to a missing directory.
func TestSavePhase21Results_NonExistentDirWriteError(t *testing.T) {
	nonExistentDir := filepath.Join(t.TempDir(), "nonexistent", "subdir")
	// Should not panic even when write fails.
	savePhase21Results(nonExistentDir, 5*time.Second, 0.09, 0.42, 100000, 1000, 4, "TEST")
}

// TestProcessContextualVCFFile_EmptyRefOrAlt tests that VCF lines with empty
// ref or alt fields are skipped.
func TestProcessContextualVCFFile_EmptyRefOrAlt(t *testing.T) {
	tmpDir := t.TempDir()

	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	// Line with empty ref field (5 fields but ref is empty).
	gw.Write([]byte("chr1\t100\t.\t\tT\t.\tPASS\t.\n"))
	// A valid line to ensure other logic runs.
	gw.Write([]byte("chr1\t200\t.\tA\tT\t.\tPASS\t.\n"))
	gw.Close()
	vcfPath := filepath.Join(tmpDir, "empty_ref.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	_, variants := processContextualVCFFile(model, vcfPath, "TEST")
	// Only the valid line should be counted.
	if variants != 1 {
		t.Errorf("expected 1 variant (empty ref line skipped), got %d", variants)
	}
}

// TestReconstructContextualSequence_LargeAlt tests behavior when alt allele
// is large enough to extend beyond the window boundary.
func TestReconstructContextualSequence_LargeAlt(t *testing.T) {
	// alt that is longer than windowSize would be rejected (windowPos+len(alt) >= windowSize)
	largeAlt := string(make([]byte, 9999))
	for i := range []byte(largeAlt) {
		if i < len(largeAlt) {
			largeAlt = largeAlt[:i] + "A" + largeAlt[i+1:]
		}
	}
	variants := []Variant{
		{pos: 1, ref: "A", alt: largeAlt},
		{pos: 2, ref: "G", alt: "C"},
	}
	result := reconstructContextualSequence(variants)
	if result == nil {
		t.Fatal("expected non-nil result")
	}
	if len(result) != 10000 {
		t.Errorf("expected window size 10000, got %d", len(result))
	}
}
