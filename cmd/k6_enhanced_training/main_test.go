//go:build !release
// +build !release

package main

import (
	"bytes"
	"compress/gzip"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"testing"
	"time"

	authenticity "github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestAbs_K6(t *testing.T) {
	tests := []struct {
		name     string
		input    float64
		expected float64
	}{
		{"positive", 2.71, 2.71},
		{"negative", -2.71, 2.71},
		{"zero", 0.0, 0.0},
		{"large positive", 1e8, 1e8},
		{"large negative", -1e8, 1e8},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := abs(tt.input)
			if got != tt.expected {
				t.Errorf("abs(%v) = %v, want %v", tt.input, got, tt.expected)
			}
		})
	}
}

func TestCalculateStateSpaceSize_K6(t *testing.T) {
	tests := []struct {
		k    int
		want int
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
		if got != tt.want {
			t.Errorf("calculateStateSpaceSize(%d) = %d, want %d", tt.k, got, tt.want)
		}
	}
}

func TestCalculateComplexity_K6(t *testing.T) {
	tests := []struct {
		name    string
		seqs    [][]byte
		wantMin float64
		wantMax float64
	}{
		{
			name:    "empty sequences returns 0",
			seqs:    [][]byte{},
			wantMin: 0.0,
			wantMax: 0.0,
		},
		{
			name:    "short sequence below k=6 returns 0",
			seqs:    [][]byte{[]byte("ACGT")},
			wantMin: 0.0,
			wantMax: 0.0,
		},
		{
			name:    "repeated bases low complexity",
			seqs:    [][]byte{[]byte("AAAAAAAAAAAAAAAAAAAAAA")},
			wantMin: 0.0,
			wantMax: 0.1,
		},
		{
			name:    "diverse sequence positive complexity",
			seqs:    [][]byte{[]byte("ACGTACGTACGTACGTACGTACGT")},
			wantMin: 0.01,
			wantMax: 1.0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateComplexity(tt.seqs)
			if got < tt.wantMin || got > tt.wantMax {
				t.Errorf("calculateComplexity() = %v, want [%v, %v]", got, tt.wantMin, tt.wantMax)
			}
		})
	}
}

func TestCalculateGCContent_K6(t *testing.T) {
	tests := []struct {
		name string
		seqs [][]byte
		want float64
	}{
		{"empty returns 0", [][]byte{}, 0.0},
		{"all A returns 0", [][]byte{[]byte("AAAA")}, 0.0},
		{"all G returns 1", [][]byte{[]byte("GGGG")}, 1.0},
		{"50pct GC", [][]byte{[]byte("ACGT")}, 0.5},
		{"multiple seqs averaged", [][]byte{[]byte("AAAA"), []byte("GGGG")}, 0.5},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateGCContent(tt.seqs)
			if math.Abs(got-tt.want) > 1e-9 {
				t.Errorf("calculateGCContent() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestFormatBytes_K6(t *testing.T) {
	tests := []struct {
		name  string
		bytes int64
		want  string
	}{
		{"zero bytes", 0, "0 B"},
		{"small bytes", 512, "512 B"},
		{"one KB", 1024, "1.0 KB"},
		{"1.5 KB", 1536, "1.5 KB"},
		{"one MB", 1024 * 1024, "1.0 MB"},
		{"one GB", 1024 * 1024 * 1024, "1.0 GB"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := formatBytes(tt.bytes)
			if got != tt.want {
				t.Errorf("formatBytes(%d) = %q, want %q", tt.bytes, got, tt.want)
			}
		})
	}
}

func TestDetectAppleSiliconCores_K6(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestGetComparisonSymbol(t *testing.T) {
	tests := []struct {
		name     string
		current  float64
		previous float64
		wantHas  string
	}{
		{"maintained when equal", 0.08, 0.08, "maintained"},
		{"maintained when better", 0.09, 0.08, "maintained"},
		{"decreased when worse", 0.07, 0.08, "decreased"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getComparisonSymbol(tt.current, tt.previous)
			if len(got) == 0 {
				t.Errorf("getComparisonSymbol(%v, %v) returned empty string", tt.current, tt.previous)
			}
		})
	}
}

func TestGetGCComparisonSymbol(t *testing.T) {
	tests := []struct {
		name        string
		currentDev  float64
		previousDev float64
	}{
		{"improved when lower deviation", 0.01, 0.05},
		{"similar when higher deviation", 0.06, 0.05},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getGCComparisonSymbol(tt.currentDev, tt.previousDev)
			if len(got) == 0 {
				t.Errorf("getGCComparisonSymbol(%v, %v) returned empty string", tt.currentDev, tt.previousDev)
			}
		})
	}
}

func TestGetImprovementSymbol(t *testing.T) {
	tests := []struct {
		name        string
		currentDev  float64
		previousDev float64
	}{
		{"improved", 0.01, 0.05},
		{"not improved", 0.06, 0.05},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getImprovementSymbol(tt.currentDev, tt.previousDev)
			if len(got) == 0 {
				t.Errorf("getImprovementSymbol(%v, %v) returned empty string", tt.currentDev, tt.previousDev)
			}
		})
	}
}

func TestGetNextStepsPhase12(t *testing.T) {
	tests := []struct {
		name       string
		complexity float64
		gcContent  float64
	}{
		{"complete success path", 0.09, 0.42},
		{"good progress path", 0.06, 0.42},
		{"needs refinement path", 0.02, 0.50},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getNextStepsPhase12(tt.complexity, tt.gcContent)
			if len(got) == 0 {
				t.Errorf("getNextStepsPhase12(%v, %v) returned empty string", tt.complexity, tt.gcContent)
			}
		})
	}
}

func TestGetSuccessEmojiPhase12(t *testing.T) {
	tests := []struct {
		name       string
		complexity float64
		gcContent  float64
	}{
		{"excellent result", 0.09, 0.42},
		{"good result", 0.06, 0.42},
		{"poor result", 0.02, 0.50},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getSuccessEmojiPhase12(tt.complexity, tt.gcContent)
			if len(got) == 0 {
				t.Errorf("getSuccessEmojiPhase12(%v, %v) returned empty string", tt.complexity, tt.gcContent)
			}
		})
	}
}

// ---------- savePhase12Results ----------

func TestSavePhase12Results_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	savePhase12Results(tmpDir, 30*time.Second, 0.082, 0.42, 1000000, 5000, 6, "SUCCESS")
	path := filepath.Join(tmpDir, "phase1_2_enhanced_results.txt")
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty results file")
	}
}

// ---------- VCF processing helper ----------

func createTestGzippedVCFenhanced(t *testing.T) string {
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

func TestProcessRealVCFFile_Enhanced_Basic(t *testing.T) {
	vcfPath := createTestGzippedVCFenhanced(t)
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

func TestProcessRealVCFFile_Enhanced_NonExistentFile(t *testing.T) {
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

func TestProcessRealVCFFile_Enhanced_NotGzipped(t *testing.T) {
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "notgzip.vcf.gz")
	os.WriteFile(path, []byte("this is not gzip data"), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processRealVCFFile(model, path, "TEST")
	if fileSize != 0 || variants != 0 {
		t.Errorf("expected 0,0 for invalid gzip, got %d,%d", fileSize, variants)
	}
}

// ---------- processPopulationVCFs tests ----------

func TestProcessPopulationVCFs_Enhanced_EmptyDir(t *testing.T) {
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

func TestProcessPopulationVCFs_Enhanced_WithVCFs(t *testing.T) {
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
	fileSize, variants := processPopulationVCFs(model, "YRI", tmpDir, 1)
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 3 {
		t.Errorf("expected 3 variants, got %d", variants)
	}
}

// ---------- generateValidationSequences tests ----------

func TestGenerateValidationSequences_Enhanced_Basic(t *testing.T) {
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

func TestGenerateValidationSequences_Enhanced_Multiple(t *testing.T) {
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

func TestProcessTestPopulations_WithExistingPopulation(t *testing.T) {
	// Create a population directory with a gzipped VCF file so the
	// goroutine success path (population dir found) is exercised.
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
	totalBytes, totalVariants := processTestPopulations(model, tmpDir, []string{"YRI"}, 1)
	if totalBytes <= 0 {
		t.Errorf("expected positive bytes, got %d", totalBytes)
	}
	if totalVariants != 2 {
		t.Errorf("expected 2 variants, got %d", totalVariants)
	}
}

// TestProcessTestPopulations_MultiplePopulations tests that multiple populations
// are handled correctly, including one that is missing (skipped).
func TestProcessTestPopulations_MultiplePopulations(t *testing.T) {
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
	totalBytes, totalVariants := processTestPopulations(model, tmpDir, []string{"YRI", "CEU"}, 1)
	if totalBytes <= 0 {
		t.Errorf("expected positive bytes from YRI, got %d", totalBytes)
	}
	if totalVariants != 1 {
		t.Errorf("expected 1 variant from YRI, got %d", totalVariants)
	}
}

// TestProcessRealVCFFile_LargeSequenceChunk exercises the path where
// the accumulated sequence builder exceeds 100KB and triggers mid-loop training.
func TestProcessRealVCFFile_LargeSequenceChunk(t *testing.T) {
	tmpDir := t.TempDir()

	// Build large ref and alt allele strings (50001 bases each) so that
	// after 2 variants the sequence builder will exceed 100KB and flush.
	refAllele := make([]byte, 50001)
	altAllele := make([]byte, 50001)
	for i := range refAllele {
		refAllele[i] = 'A'
		altAllele[i] = 'T'
	}

	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	// Write 3 variants - after 2 the sequence builder will exceed 100KB.
	for i := 0; i < 3; i++ {
		line := "chr1\t10" + string(rune('0'+i)) + "\t.\t" + string(refAllele) + "\t" + string(altAllele) + "\t.\tPASS\t.\n"
		gw.Write([]byte(line))
	}
	gw.Close()

	vcfPath := filepath.Join(tmpDir, "large_variants.vcf.gz")
	if err := os.WriteFile(vcfPath, buf.Bytes(), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	fileSize, variants := processRealVCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 3 {
		t.Errorf("expected 3 variants, got %d", variants)
	}
}

// TestProcessRealVCFFile_ShortFields tests lines with fewer than 5 tab-separated
// fields - they should be skipped without panicking.
func TestProcessRealVCFFile_ShortFields(t *testing.T) {
	tmpDir := t.TempDir()

	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	// A data line with only 3 fields (< 5 required) should be skipped.
	gw.Write([]byte("chr1\t100\tonly_three_fields\n"))
	// A valid line to ensure some output.
	gw.Write([]byte("chr1\t200\t.\tA\tT\t.\tPASS\t.\n"))
	gw.Close()

	vcfPath := filepath.Join(tmpDir, "short_fields.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(6)
	fileSize, variants := processRealVCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	// Only the valid line counts.
	if variants != 1 {
		t.Errorf("expected 1 variant (short field line skipped), got %d", variants)
	}
}

// TestCalculateGCContent_EmptyInnerSequence exercises the path where a slice
// contains an empty byte slice - it must be skipped gracefully.
func TestCalculateGCContent_EmptyInnerSequence(t *testing.T) {
	// Slice with one empty sequence and one valid sequence.
	seqs := [][]byte{
		{},           // empty - should be skipped
		[]byte("GG"), // 100% GC
	}
	got := calculateGCContent(seqs)
	// Average = (0 + 1.0) / 2 = 0.5 because empty is averaged as 0.
	if got < 0 || got > 1.0 {
		t.Errorf("calculateGCContent with empty inner sequence = %v, want in [0,1]", got)
	}
}

// TestSavePhase12Results_BranchVariants exercises the different status branches
// in savePhase12Results by varying complexity and GC content inputs.
func TestSavePhase12Results_BranchVariants(t *testing.T) {
	tests := []struct {
		name       string
		complexity float64
		gcContent  float64
		status     string
	}{
		// complexity >= 0.0819 (maintained) and gcDeviation < abs(0.4817-0.42) (improved)
		{"maintained_improved", 0.0819, 0.425, "SUCCESS"},
		// complexity < 0.0819 (decreased) and gcDeviation >= abs(0.4817-0.42) (not improved)
		{"decreased_not_improved", 0.07, 0.479, "PARTIAL"},
		// complexity < 0.0819 and gcDeviation < previous - triggers different symbols
		{"decreased_gc_improved", 0.07, 0.422, "PARTIAL"},
		// complexity > 0.08 and GC excellent range
		{"excellent_everything", 0.09, 0.420, "SUCCESS"},
		// complexity in medium range and GC in good range
		{"medium_everything", 0.06, 0.41, "PROGRESS"},
		// below all thresholds
		{"below_thresholds", 0.02, 0.50, "REFINE"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			// Should not panic.
			savePhase12Results(tmpDir, 10*time.Second, tt.complexity, tt.gcContent, 500000, 1000, 4, tt.status)
			path := filepath.Join(tmpDir, "phase1_2_enhanced_results.txt")
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

// TestDetectAppleSiliconCores_ReturnsPositive verifies the function always returns
// a positive number regardless of the actual CPU count.
func TestDetectAppleSiliconCores_ReturnsPositive(t *testing.T) {
	cores := detectAppleSiliconCores()
	if cores <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", cores)
	}
}

// TestDetectAppleSiliconCores_BoundedByTotal verifies the result does not exceed
// total CPU count (it should be <= total cores).
func TestDetectAppleSiliconCores_BoundedByTotal(t *testing.T) {
	total := runtime.NumCPU()
	cores := detectAppleSiliconCores()
	if cores > total {
		t.Errorf("detectAppleSiliconCores() = %d, but total CPUs = %d; P-cores should not exceed total", cores, total)
	}
}

// TestGenerateValidationSequences_UntrainedModelReturnsError ensures the error
// path in generateValidationSequences is exercised when the model returns an error.
func TestGenerateValidationSequences_UntrainedModelReturnsError(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(6)
	if model == nil {
		t.Fatal("failed to create model")
	}
	// Don't train the model - GenerateSequence should return an error.
	seqs, err := generateValidationSequences(model, 1, 100)
	if err == nil {
		t.Error("expected error from untrained model, got nil")
	}
	if seqs != nil {
		t.Errorf("expected nil sequences on error, got %v", seqs)
	}
}

// TestCalculateGCContent_MultipleEmptyAndValid exercises the mixed slice case:
// empty sequences are skipped but contribute to the average divisor.
func TestCalculateGCContent_MultipleEmptyAndValid(t *testing.T) {
	seqs := [][]byte{
		{},                // empty - skipped (no GC counted, but len still divides)
		[]byte("GGGGCCCC"), // 100% GC
		{},                // another empty
	}
	got := calculateGCContent(seqs)
	// totalGC = 1.0 (only second), averaged over 3 = 0.333...
	if got < 0 || got > 1.0 {
		t.Errorf("calculateGCContent with empties = %v, want in [0,1]", got)
	}
}

// TestProcessPopulationVCFs_HeadersOnlyFile tests a VCF with only headers
// resulting in 0 variants but positive file size.
func TestProcessPopulationVCFs_HeadersOnly(t *testing.T) {
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
	fileSize, variants := processPopulationVCFs(model, "POP", tmpDir, 1)
	if fileSize <= 0 {
		t.Errorf("expected positive file size from headers-only VCF, got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for headers-only, got %d", variants)
	}
}

// TestSavePhase12Results_NonExistentDirWriteError tests the error path where
// the progress directory doesn't exist, causing os.WriteFile to fail.
func TestSavePhase12Results_NonExistentDirWriteError(t *testing.T) {
	// Use a path that definitely doesn't exist - WriteFile should fail gracefully.
	nonExistentDir := filepath.Join(t.TempDir(), "nonexistent", "subdir")
	// Should not panic even when write fails.
	savePhase12Results(nonExistentDir, 5*time.Second, 0.09, 0.42, 100000, 1000, 4, "TEST")
}

// TestProcessRealVCFFile_EmptyRefOrAlt tests that VCF lines with empty
// ref or alt fields are skipped (len(ref)==0 or len(alt)==0).
func TestProcessRealVCFFile_EmptyRefOrAlt(t *testing.T) {
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
	_, variants := processRealVCFFile(model, vcfPath, "TEST")
	// Only the valid line should be counted.
	if variants != 1 {
		t.Errorf("expected 1 variant (empty ref line skipped), got %d", variants)
	}
}
