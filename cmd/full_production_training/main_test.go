//go:build !release
// +build !release

package main

import (
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
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

func TestCalculateComplexityEmpty(t *testing.T) {
	if got := calculateComplexity(nil, 6); got != 0 {
		t.Errorf("expected 0 for nil, got %f", got)
	}
	if got := calculateComplexity([][]byte{}, 6); got != 0 {
		t.Errorf("expected 0 for empty slice, got %f", got)
	}
}

func TestCalculateComplexityShortSequence(t *testing.T) {
	seqs := [][]byte{[]byte("ACG")}
	if got := calculateComplexity(seqs, 6); got != 0 {
		t.Errorf("expected 0 for sequence shorter than k, got %f", got)
	}
}

func TestCalculateComplexityUniform(t *testing.T) {
	seq := make([]byte, 100)
	for i := range seq {
		seq[i] = 'A'
	}
	got := calculateComplexity([][]byte{seq}, 6)
	if got > 0.1 {
		t.Errorf("uniform sequence should have very low complexity, got %f", got)
	}
}

func TestCalculateComplexityBounded(t *testing.T) {
	seq := []byte("ACGTACGTGCATGCAT")
	got := calculateComplexity([][]byte{seq}, 4)
	if got < 0 || got > 1 {
		t.Errorf("complexity should be in [0,1], got %f", got)
	}
}

func TestCalculateComplexityDeterministic(t *testing.T) {
	seq := []byte("ACGTACGTGCATGCATACGT")
	c1 := calculateComplexity([][]byte{seq}, 4)
	c2 := calculateComplexity([][]byte{seq}, 4)
	if c1 != c2 {
		t.Errorf("complexity not deterministic: %f != %f", c1, c2)
	}
}

func TestCalculateGCContentEmpty(t *testing.T) {
	if got := calculateGCContent(nil); got != 0 {
		t.Errorf("expected 0 for nil, got %f", got)
	}
	if got := calculateGCContent([][]byte{}); got != 0 {
		t.Errorf("expected 0 for empty slice, got %f", got)
	}
	if got := calculateGCContent([][]byte{{}}); got != 0 {
		t.Errorf("expected 0 for empty sequence, got %f", got)
	}
}

func TestCalculateGCContentAllGC(t *testing.T) {
	seq := []byte("GCGCGCGCGC")
	if got := calculateGCContent([][]byte{seq}); got != 1.0 {
		t.Errorf("expected 1.0, got %f", got)
	}
}

func TestCalculateGCContentAllAT(t *testing.T) {
	seq := []byte("ATATATATAT")
	if got := calculateGCContent([][]byte{seq}); got != 0.0 {
		t.Errorf("expected 0.0, got %f", got)
	}
}

func TestCalculateGCContentHalf(t *testing.T) {
	seq := []byte("AACCGGTT")
	got := calculateGCContent([][]byte{seq})
	if got != 0.5 {
		t.Errorf("expected 0.5, got %f", got)
	}
}

func TestCalculateGCContentMultiple(t *testing.T) {
	seqGC := []byte("GCGC")
	seqAT := []byte("ATAT")
	got := calculateGCContent([][]byte{seqGC, seqAT})
	if abs(got-0.5) > 1e-9 {
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

func TestGetOptimalProductionConfig(t *testing.T) {
	cfg := getOptimalProductionConfig()
	if cfg.KmerSize <= 0 {
		t.Errorf("expected positive KmerSize, got %d", cfg.KmerSize)
	}
	if cfg.TargetComplexity <= 0 {
		t.Errorf("expected positive TargetComplexity, got %f", cfg.TargetComplexity)
	}
	if cfg.TargetGCContent <= 0 || cfg.TargetGCContent >= 1 {
		t.Errorf("expected TargetGCContent in (0,1), got %f", cfg.TargetGCContent)
	}
	if len(cfg.PopulationWeights) == 0 {
		t.Error("expected non-empty PopulationWeights")
	}
	if cfg.ProcessingApproach == "" {
		t.Error("expected non-empty ProcessingApproach")
	}
}

func TestPopulationWeightsAllPositive(t *testing.T) {
	cfg := getOptimalProductionConfig()
	for pop, w := range cfg.PopulationWeights {
		if w <= 0 {
			t.Errorf("population %s has non-positive weight %f", pop, w)
		}
	}
}

func TestGetThresholdStatus(t *testing.T) {
	if got := getThresholdStatus(0.09, 0.08); got != "✅" {
		t.Errorf("above threshold: want ✅, got %q", got)
	}
	if got := getThresholdStatus(0.08, 0.08); got != "✅" {
		t.Errorf("at threshold: want ✅, got %q", got)
	}
	if got := getThresholdStatus(0.07, 0.08); got != "❌" {
		t.Errorf("below threshold: want ❌, got %q", got)
	}
}

func TestGetGCStatus(t *testing.T) {
	if got := getGCStatus(0.01); got != "✅" {
		t.Errorf("deviation 0.01: want ✅, got %q", got)
	}
	if got := getGCStatus(0.02); got != "✅" {
		t.Errorf("deviation 0.02: want ✅, got %q", got)
	}
	if got := getGCStatus(0.03); got != "❌" {
		t.Errorf("deviation 0.03: want ❌, got %q", got)
	}
}

func TestGetValidationReadinessDiffers(t *testing.T) {
	r := getValidationReadiness(true)
	nr := getValidationReadiness(false)
	if r == "" || nr == "" {
		t.Error("expected non-empty strings")
	}
	if r == nr {
		t.Error("ready and not-ready should differ")
	}
}

func TestGetDeploymentReadinessDiffers(t *testing.T) {
	r := getDeploymentReadiness(true)
	nr := getDeploymentReadiness(false)
	if r == "" || nr == "" {
		t.Error("expected non-empty strings")
	}
	if r == nr {
		t.Error("ready and not-ready should differ")
	}
}

func TestGetProductionEmoji(t *testing.T) {
	if got := getProductionEmoji(true); got != "✅" {
		t.Errorf("ready=true: want ✅, got %q", got)
	}
	if got := getProductionEmoji(false); got != "🔄" {
		t.Errorf("ready=false: want 🔄, got %q", got)
	}
}

func TestGetNextStepsProductionDiffers(t *testing.T) {
	r := getNextStepsProduction(true)
	nr := getNextStepsProduction(false)
	if r == "" || nr == "" {
		t.Error("expected non-empty strings")
	}
	if r == nr {
		t.Error("ready and not-ready should differ")
	}
}

func TestSaveProductionResults_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	config := ProductionConfig{
		KmerSize:           6,
		PopulationWeights:  map[string]float64{"YRI": 2.0},
		TargetComplexity:   0.08,
		TargetGCContent:    0.42,
		ProcessingApproach: "weighted",
		OptimizationLevel:  "standard",
	}
	saveProductionResults(tmpDir, time.Second*30, 0.082, 0.42, 1000000, 5000, 4, config, true)
	path := filepath.Join(tmpDir, "production_final_results.txt")
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty results file")
	}
}

// --- Tests for generateValidationSequences ---

func TestGenerateValidationSequences_UntrainedModelErrors(t *testing.T) {
	// An untrained AdaptiveHMMModel returns an error from GenerateSequence.
	model := authenticity.NewAdaptiveHMMModel(4)
	_, err := generateValidationSequences(model, 1, 100)
	if err == nil {
		t.Error("expected error from untrained model, got nil")
	}
}

func TestGenerateValidationSequences_ZeroSequences(t *testing.T) {
	// Requesting 0 sequences should succeed immediately and return an empty slice.
	model := authenticity.NewAdaptiveHMMModel(4)
	seqs, err := generateValidationSequences(model, 0, 100)
	if err != nil {
		t.Fatalf("unexpected error for 0 sequences: %v", err)
	}
	if len(seqs) != 0 {
		t.Errorf("expected 0 sequences, got %d", len(seqs))
	}
}

func TestGenerateValidationSequences_TrainedModelSucceeds(t *testing.T) {
	// Train a model on a short sequence so it can generate output.
	model := authenticity.NewAdaptiveHMMModel(4)
	dna := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(dna); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}
	seqs, err := generateValidationSequences(model, 2, 20)
	if err != nil {
		t.Fatalf("expected success with trained model, got: %v", err)
	}
	if len(seqs) != 2 {
		t.Errorf("expected 2 sequences, got %d", len(seqs))
	}
	for i, s := range seqs {
		if len(s) != 20 {
			t.Errorf("sequence %d: expected length 20, got %d", i, len(s))
		}
	}
}

// --- Tests for processProductionPopulations ---

func TestProcessProductionPopulations_NoPopulationDirs(t *testing.T) {
	// populationDir exists but no subdirectories for the populations are present,
	// so all populations are skipped and totals remain 0.
	populationDir := t.TempDir()
	model := authenticity.NewAdaptiveHMMModel(4)
	config := getOptimalProductionConfig()
	populations := []string{"YRI", "CEU"}
	totalBytes, totalVariants := processProductionPopulations(model, populationDir, populations, config, 1)
	if totalBytes != 0 {
		t.Errorf("expected 0 bytes, got %d", totalBytes)
	}
	if totalVariants != 0 {
		t.Errorf("expected 0 variants, got %d", totalVariants)
	}
}

func TestProcessProductionPopulations_EmptyPopulationList(t *testing.T) {
	populationDir := t.TempDir()
	model := authenticity.NewAdaptiveHMMModel(4)
	config := getOptimalProductionConfig()
	totalBytes, totalVariants := processProductionPopulations(model, populationDir, []string{}, config, 1)
	if totalBytes != 0 || totalVariants != 0 {
		t.Errorf("expected (0,0), got (%d,%d)", totalBytes, totalVariants)
	}
}

func TestProcessProductionPopulations_PopDirExistsButEmpty(t *testing.T) {
	// Create a real population subdirectory but with no VCF files.
	populationDir := t.TempDir()
	if err := os.MkdirAll(filepath.Join(populationDir, "YRI"), 0750); err != nil {
		t.Fatal(err)
	}
	model := authenticity.NewAdaptiveHMMModel(4)
	config := getOptimalProductionConfig()
	totalBytes, totalVariants := processProductionPopulations(model, populationDir, []string{"YRI"}, config, 1)
	// No VCF files -> no bytes, no variants processed.
	if totalBytes != 0 {
		t.Errorf("expected 0 bytes, got %d", totalBytes)
	}
	if totalVariants != 0 {
		t.Errorf("expected 0 variants, got %d", totalVariants)
	}
}

// --- Tests for processWeightedPopulationVCFs ---

func TestProcessWeightedPopulationVCFs_NonexistentDir(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	// Directory does not exist: filepath.Glob returns empty, so 0 files processed.
	bytes, variants := processWeightedPopulationVCFs(model, "YRI", "/nonexistent/path/that/does/not/exist", 1.0, 1)
	if bytes != 0 || variants != 0 {
		t.Errorf("expected (0,0) for nonexistent dir, got (%d,%d)", bytes, variants)
	}
}

func TestProcessWeightedPopulationVCFs_EmptyPopDir(t *testing.T) {
	// populationDir/pop exists but is empty (no .vcf.gz files).
	populationDir := t.TempDir()
	if err := os.MkdirAll(filepath.Join(populationDir, "CEU"), 0750); err != nil {
		t.Fatal(err)
	}
	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedPopulationVCFs(model, "CEU", populationDir, 1.0, 1)
	if bytes != 0 || variants != 0 {
		t.Errorf("expected (0,0), got (%d,%d)", bytes, variants)
	}
}

// makeGzippedVCF creates a minimal gzipped VCF file in dir with the given content.
func makeGzippedVCF(t *testing.T, dir, filename string, lines []string) string {
	t.Helper()
	path := filepath.Join(dir, filename)
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("failed to create vcf.gz: %v", err)
	}
	gw := gzip.NewWriter(f)
	for _, line := range lines {
		if _, werr := gw.Write([]byte(line + "\n")); werr != nil {
			t.Fatalf("failed to write gzip data: %v", werr)
		}
	}
	if err := gw.Close(); err != nil {
		t.Fatalf("failed to close gzip writer: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("failed to close file: %v", err)
	}
	return path
}

func TestProcessWeightedPopulationVCFs_WithVCFFile(t *testing.T) {
	populationDir := t.TempDir()
	popDir := filepath.Join(populationDir, "YRI")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatal(err)
	}
	// A minimal VCF: header + two variant lines.
	lines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t100\t.\tACGT\tTGCA\t.\t.\t.",
		"chr1\t200\t.\tAAAA\tCCCC\t.\t.\t.",
	}
	makeGzippedVCF(t, popDir, "test.vcf.gz", lines)
	model := authenticity.NewAdaptiveHMMModel(4)
	if err := model.TrainFromSequence([]byte("ACGTACGTACGT")); err != nil {
		t.Fatal(err)
	}
	bytes, variants := processWeightedPopulationVCFs(model, "YRI", populationDir, 1.0, 1)
	if bytes == 0 {
		t.Error("expected non-zero bytes from VCF file")
	}
	if variants == 0 {
		t.Error("expected non-zero variants from VCF file")
	}
}

// --- Tests for processWeightedVCFFile ---

func TestProcessWeightedVCFFile_NonexistentFile(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedVCFFile(model, "/nonexistent/file.vcf.gz", "YRI", 1.0)
	if bytes != 0 || variants != 0 {
		t.Errorf("expected (0,0) for nonexistent file, got (%d,%d)", bytes, variants)
	}
}

func TestProcessWeightedVCFFile_NotGzipped(t *testing.T) {
	// A plain text file that is not valid gzip causes gzip.NewReader to fail -> returns (0,0).
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "plain.vcf.gz")
	if err := os.WriteFile(path, []byte("this is not gzip data"), 0600); err != nil {
		t.Fatal(err)
	}
	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedVCFFile(model, path, "YRI", 1.0)
	// Invalid gzip returns (0,0) per the error handling in the function.
	if bytes != 0 || variants != 0 {
		t.Errorf("expected (0,0) for invalid gzip, got (%d,%d)", bytes, variants)
	}
}

func TestProcessWeightedVCFFile_ValidVCF(t *testing.T) {
	tmpDir := t.TempDir()
	lines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t100\t.\tACGT\tTGCA\t.\t.\t.",
		"chr1\t150\t.\tA\tT\t.\t.\t.",
	}
	path := makeGzippedVCF(t, tmpDir, "test.vcf.gz", lines)
	model := authenticity.NewAdaptiveHMMModel(4)
	if err := model.TrainFromSequence([]byte("ACGTACGTACGT")); err != nil {
		t.Fatal(err)
	}
	bytes, variants := processWeightedVCFFile(model, path, "YRI", 1.0)
	if bytes == 0 {
		t.Error("expected non-zero file size")
	}
	if variants != 2 {
		t.Errorf("expected 2 variants, got %d", variants)
	}
}

func TestProcessWeightedVCFFile_HeaderOnlyVCF(t *testing.T) {
	// VCF with only header lines -> 0 variants.
	tmpDir := t.TempDir()
	lines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}
	path := makeGzippedVCF(t, tmpDir, "header_only.vcf.gz", lines)
	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedVCFFile(model, path, "YRI", 1.0)
	if bytes == 0 {
		t.Error("expected non-zero file size for a valid gzip file")
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for header-only VCF, got %d", variants)
	}
}

func TestProcessWeightedVCFFile_InsufficientFieldsLine(t *testing.T) {
	// Lines with fewer than 5 tab-separated fields are skipped.
	tmpDir := t.TempDir()
	lines := []string{
		"##fileformat=VCFv4.1",
		"#CHROM\tPOS",
		"chr1\t100\t.\tA", // only 4 fields - should be skipped
	}
	path := makeGzippedVCF(t, tmpDir, "short.vcf.gz", lines)
	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedVCFFile(model, path, "YRI", 1.0)
	if bytes == 0 {
		t.Error("expected non-zero bytes even for skipped-line VCF")
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for short-field lines, got %d", variants)
	}
}

func TestProcessWeightedVCFFile_WeightGreaterThanOne(t *testing.T) {
	// weight=2.5 -> repetitions should be 2 (int(2.5)=2, but check 0.5 > 0.5 is false -> stays 2).
	tmpDir := t.TempDir()
	lines := []string{
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t100\t.\tACGT\tTGCA\t.\t.\t.",
	}
	path := makeGzippedVCF(t, tmpDir, "w2.vcf.gz", lines)
	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedVCFFile(model, path, "YRI", 2.5)
	if bytes == 0 {
		t.Error("expected non-zero bytes")
	}
	if variants != 1 {
		t.Errorf("expected 1 variant, got %d", variants)
	}
}
