//go:build !release
// +build !release

package main

import (
	"compress/gzip"
	"fmt"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// createGzipVCFPWT2 creates a minimal gzipped VCF file for testing.
func createGzipVCFPWT2(t *testing.T, dir, name string, vcfLines []string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	f, err := os.Create(p)
	if err != nil {
		t.Fatalf("createGzipVCFPWT2: %v", err)
	}
	defer f.Close()
	gw := gzip.NewWriter(f)
	for _, line := range vcfLines {
		_, _ = gw.Write(append([]byte(line), '\n'))
	}
	if err := gw.Close(); err != nil {
		t.Fatalf("gzip.Close: %v", err)
	}
	return p
}

// TestProcessWeightedPopulationVCFs_WithVCFFile_PWT2 creates a population
// directory with a real vcf.gz file and exercises the gzip-reading code path.
func TestProcessWeightedPopulationVCFs_WithVCFFile_PWT2(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}
	vcfLines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t1000\t.\tACGT\tGCAT\t.\tPASS\t.",
		"chr1\t2000\t.\tA\tG\t.\tPASS\t.",
	}
	createGzipVCFPWT2(t, popDir, "chr1.vcf.gz", vcfLines)

	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedPopulationVCFs(model, "YRI", tmpDir, 2.0, 2)
	if bytes == 0 {
		t.Error("expected non-zero bytes from a real vcf.gz file")
	}
	if variants == 0 {
		t.Error("expected non-zero variants from VCF with data lines")
	}
}

// TestProcessWeightedVCFFile_LargeSequenceChunk_PWT2 exercises the 100KB flush path.
func TestProcessWeightedVCFFile_LargeSequenceChunk_PWT2(t *testing.T) {
	tmpDir := t.TempDir()
	ref := "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	alt := "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
	vcfLines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}
	for i := 0; i < 2200; i++ {
		vcfLines = append(vcfLines, fmt.Sprintf("chr1\t%d\t.\t%s\t%s\t.\tPASS\t.", 1000+i, ref, alt))
	}
	p := createGzipVCFPWT2(t, tmpDir, "bigchunk.vcf.gz", vcfLines)

	model := authenticity.NewAdaptiveHMMModel(4)
	if err := model.TrainFromSequence([]byte("ACGTACGTACGTACGT")); err != nil {
		t.Fatalf("initial train: %v", err)
	}

	bytes, variants := processWeightedVCFFile(model, p, "YRI", 2.0)
	if bytes == 0 {
		t.Error("expected non-zero bytes")
	}
	if variants == 0 {
		t.Error("expected non-zero variants")
	}
}

// TestProcessWeightedVCFFile_FractionalWeight_PWT2 exercises the fractional-weight branch.
func TestProcessWeightedVCFFile_FractionalWeight_PWT2(t *testing.T) {
	tmpDir := t.TempDir()
	vcfLines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t1000\t.\tACGT\tGCAT\t.\tPASS\t.",
	}
	p := createGzipVCFPWT2(t, tmpDir, "frac.vcf.gz", vcfLines)
	model := authenticity.NewAdaptiveHMMModel(4)
	_, variants := processWeightedVCFFile(model, p, "CEU", 1.8)
	if variants == 0 {
		t.Error("expected non-zero variants with weight=1.8")
	}
}

// TestProcessWeightedVCFFile_InvalidGzip_PWT2 exercises the gzip-reader error path.
func TestProcessWeightedVCFFile_InvalidGzip_PWT2(t *testing.T) {
	tmpDir := t.TempDir()
	p := filepath.Join(tmpDir, "bad.vcf.gz")
	if err := os.WriteFile(p, []byte("not gzip data"), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedVCFFile(model, p, "YRI", 1.0)
	if bytes != 0 {
		t.Errorf("expected 0 bytes for invalid gzip, got %d", bytes)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for invalid gzip, got %d", variants)
	}
}

// TestCalculateGCContent_EmptySequenceInSlice_PWT2 verifies empty entries don't panic.
func TestCalculateGCContent_EmptySequenceInSlice_PWT2(t *testing.T) {
	seqs := [][]byte{
		{},
		[]byte("GCGC"),
	}
	got := calculateGCContent(seqs)
	if got < 0 || got > 1 {
		t.Errorf("calculateGCContent = %v, want in [0,1]", got)
	}
}

// TestGenerateValidationSequences_PopWeighted_Error_PWT2 tests the untrained path.
func TestGenerateValidationSequences_PopWeighted_Error_PWT2(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	_, err := generateValidationSequences(model, 1, 50)
	_ = err
}

// TestProcessWeightedPopulations_WithVCFFiles_PWT2 exercises goroutine + semaphore.
func TestProcessWeightedPopulations_WithVCFFiles_PWT2(t *testing.T) {
	tmpDir := t.TempDir()
	weights := initializePopulationWeights()
	vcfLines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t1000\t.\tA\tG\t.\tPASS\t.",
	}
	for _, pop := range []string{"YRI", "CEU"} {
		popDir := filepath.Join(tmpDir, pop)
		if err := os.MkdirAll(popDir, 0750); err != nil {
			t.Fatalf("MkdirAll(%s): %v", pop, err)
		}
		createGzipVCFPWT2(t, popDir, "chr1.vcf.gz", vcfLines)
	}

	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedPopulations(model, tmpDir, []string{"YRI", "CEU"}, weights, 2)
	if bytes == 0 {
		t.Error("expected non-zero bytes from populations with VCF files")
	}
	if variants == 0 {
		t.Error("expected non-zero variants from populations with VCF files")
	}
}

// TestCalculateComplexity_ShortSeq_PWT2 exercises the len(seq) < k path.
func TestCalculateComplexity_ShortSeq_PWT2(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACGTA"),
		[]byte("ACGTACGTACGT"),
	}
	got := calculateComplexity(seqs)
	if got < 0 || got > 1 {
		t.Errorf("calculateComplexity = %v, want in [0,1]", got)
	}
}

// TestSavePhase22Results_WriteError_PWT exercises the os.WriteFile error path
// by using a nonexistent parent directory.
func TestSavePhase22Results_WriteError_PWT(t *testing.T) {
	weights := initializePopulationWeights()
	// Save to a nonexistent directory - WriteFile should fail and print error
	savePhase22Results("/nonexistent/dir/progress", 10*time.Second, 0.08, 0.42, 1000, 100, 4, "test", weights)
	// Should not panic
}

// TestProcessWeightedVCFFile_ShortFields_PWT exercises the len(fields) < 5 branch.
func TestProcessWeightedVCFFile_ShortFields_PWT(t *testing.T) {
	tmpDir := t.TempDir()
	vcfLines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t1000\t.",           // only 3 fields < 5 → skipped
		"chr1\t2000\t.\tA\tG",    // exactly 5 fields
	}
	p := createGzipVCFPWT2(t, tmpDir, "shortfields.vcf.gz", vcfLines)
	model := authenticity.NewAdaptiveHMMModel(4)
	_, variants := processWeightedVCFFile(model, p, "YRI", 1.0)
	// The 5-field line has valid REF and ALT, so 1 variant should be counted
	if variants == 0 {
		t.Error("expected non-zero variants from 5-field VCF line")
	}
}

// TestProcessWeightedVCFFile_WeightExactly2_PWT tests weight=2.0 (integer, no fractional branch).
func TestProcessWeightedVCFFile_WeightExactly2_PWT(t *testing.T) {
	tmpDir := t.TempDir()
	vcfLines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t1000\t.\tA\tG\t.\tPASS\t.",
	}
	p := createGzipVCFPWT2(t, tmpDir, "weight2.vcf.gz", vcfLines)
	model := authenticity.NewAdaptiveHMMModel(4)
	// weight=2.0: int(2.0)=2, 2.0 == float64(2) so no fractional branch
	_, variants := processWeightedVCFFile(model, p, "YRI", 2.0)
	if variants == 0 {
		t.Error("expected non-zero variants with weight=2.0")
	}
}

// TestProcessWeightedVCFFile_EmptyRefAlt_PWT tests lines where ref or alt is empty.
func TestProcessWeightedVCFFile_EmptyRefAlt_PWT(t *testing.T) {
	tmpDir := t.TempDir()
	vcfLines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t1000\t.\t\tG\t.\tPASS\t.",  // empty REF
		"chr1\t2000\t.\tA\t\t.\tPASS\t.",  // empty ALT
	}
	p := createGzipVCFPWT2(t, tmpDir, "emptyallele.vcf.gz", vcfLines)
	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedVCFFile(model, p, "YRI", 1.0)
	// Both lines have empty allele; neither should be counted as variant
	if variants != 0 {
		t.Errorf("expected 0 variants for empty alleles, got %d", variants)
	}
	_ = bytes
}

