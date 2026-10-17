//go:build !release
// +build !release

package main

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"testing"
	"time"

	authenticity "github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestAbs_K8(t *testing.T) {
	tests := []struct {
		name     string
		input    float64
		expected float64
	}{
		{"positive", 1.23, 1.23},
		{"negative", -1.23, 1.23},
		{"zero", 0.0, 0.0},
		{"large positive", 1e9, 1e9},
		{"large negative", -1e9, 1e9},
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

func TestCalculateStateSpaceSize_K8(t *testing.T) {
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

func TestStateSpaceDoublesWithEachKStep(t *testing.T) {
	// Each increment in k multiplies state space by 4
	prev := calculateStateSpaceSize(2)
	for k := 3; k <= 8; k++ {
		curr := calculateStateSpaceSize(k)
		if curr != prev*4 {
			t.Errorf("calculateStateSpaceSize(%d)=%d should be 4x calculateStateSpaceSize(%d)=%d",
				k, curr, k-1, prev)
		}
		prev = curr
	}
}

func TestCalculateMemoryFootprint(t *testing.T) {
	// k=7 and k=8 should have positive memory in MB
	// k=6 and below round to 0 MB due to integer division
	mem7 := calculateMemoryFootprint(7)
	mem8 := calculateMemoryFootprint(8)
	if mem7 <= 0 {
		t.Errorf("calculateMemoryFootprint(7) = %d MB, want > 0", mem7)
	}
	if mem8 <= 0 {
		t.Errorf("calculateMemoryFootprint(8) = %d MB, want > 0", mem8)
	}
	if mem8 <= mem7 {
		t.Errorf("expected memory for k=8 (%d MB) > k=7 (%d MB)", mem8, mem7)
	}
}

func TestCalculateComplexity_K8(t *testing.T) {
	// Note: k=8 means sequences shorter than 8 contribute 0 complexity.
	// For all-same-base sequences, complexity = 1/(len-k+1).
	// Need len >= 100 for ratio to be clearly small.
	allA100 := make([]byte, 100)
	for i := range allA100 {
		allA100[i] = 'A'
	}

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
			name:    "short sequence below k=8 returns 0",
			seqs:    [][]byte{[]byte("ACGT")},
			wantMin: 0.0,
			wantMax: 0.0,
		},
		{
			name:    "all same bases long sequence low complexity",
			seqs:    [][]byte{allA100},
			wantMin: 0.0,
			wantMax: 0.02,
		},
		{
			name:    "diverse sequence positive complexity",
			seqs:    [][]byte{[]byte("ACGTACGTACGTACGTACGTACGTACGT")},
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

func TestCalculateGCContent_K8(t *testing.T) {
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

func TestFormatBytes_K8(t *testing.T) {
	tests := []struct {
		name  string
		bytes int64
		want  string
	}{
		{"zero bytes", 0, "0 B"},
		{"small bytes", 512, "512 B"},
		{"one KB", 1024, "1.0 KB"},
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

func TestDetectAppleSiliconCores_K8(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestSafeMemoryRatioK8vsK7(t *testing.T) {
	// Both k=7 and k=8 have positive memory footprint
	ratio := safeMemoryRatio(8, 7)
	if ratio <= 1.0 {
		t.Errorf("safeMemoryRatio(8, 7) = %v, want > 1.0 (k=8 uses more memory)", ratio)
	}
}

func TestSafeMemoryRatioNoPanic(t *testing.T) {
	// Should not panic even with small k values that produce 0 MB footprint
	// safeMemoryRatio returns absolute value when denominator is 0
	got := safeMemoryRatio(8, 1)
	if got < 0 {
		t.Errorf("safeMemoryRatio(8, 1) = %v, want >= 0", got)
	}
}

func TestSafeMemoryRatioSameK(t *testing.T) {
	// Same k should give a stable result (either 1.0 or the absolute value if 0 denominator)
	got := safeMemoryRatio(8, 8)
	if got < 0 {
		t.Errorf("safeMemoryRatio(8, 8) = %v, want >= 0", got)
	}
}

func TestGetK8OptimizationInsight(t *testing.T) {
	tests := []struct {
		name           string
		complexityVsK6 float64
	}{
		{"significant improvement", 1.5},
		{"modest improvement", 1.1},
		{"no improvement", 0.9},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getK8OptimizationInsight(tt.complexityVsK6)
			if len(got) == 0 {
				t.Errorf("getK8OptimizationInsight(%v) returned empty string", tt.complexityVsK6)
			}
		})
	}
}

func TestGetNextStepsPhase3(t *testing.T) {
	tests := []struct {
		name           string
		complexity     float64
		gcContent      float64
		complexityVsK6 float64
	}{
		{"stretch goal success", 0.11, 0.42, 1.5},
		{"significant improvement", 0.06, 0.42, 1.3},
		{"marginal benefit", 0.06, 0.42, 1.0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getNextStepsPhase3(tt.complexity, tt.gcContent, tt.complexityVsK6)
			if len(got) == 0 {
				t.Errorf("getNextStepsPhase3(%v, %v, %v) returned empty string",
					tt.complexity, tt.gcContent, tt.complexityVsK6)
			}
		})
	}
}

func TestGetComputationalRecommendation(t *testing.T) {
	tests := []struct {
		name           string
		complexityVsK6 float64
	}{
		{"k8 recommended", 1.5},
		{"k6 recommended", 1.0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getComputationalRecommendation(tt.complexityVsK6, 0)
			if len(got) == 0 {
				t.Errorf("getComputationalRecommendation(%v, 0) returned empty string", tt.complexityVsK6)
			}
		})
	}
}

func TestGetSuccessEmojiPhase3(t *testing.T) {
	tests := []struct {
		name           string
		complexity     float64
		complexityVsK6 float64
	}{
		{"stretch goal achieved", 0.11, 1.5},
		{"significant improvement", 0.09, 1.3},
		{"moderate improvement", 0.09, 1.1},
		{"insufficient", 0.04, 0.9},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getSuccessEmojiPhase3(tt.complexity, tt.complexityVsK6)
			if len(got) == 0 {
				t.Errorf("getSuccessEmojiPhase3(%v, %v) returned empty string",
					tt.complexity, tt.complexityVsK6)
			}
		})
	}
}

// TestCalculateGCContent_EmptySequenceInSlice tests that empty sequences within
// a slice are skipped, covering the len(seq)==0 branch.
func TestCalculateGCContent_EmptySequenceInSlice(t *testing.T) {
	// A slice that contains one empty sequence alongside a real one.
	// The empty sequence should be skipped, so the result equals the GC of
	// the non-empty sequence only.
	seqs := [][]byte{
		{},
		[]byte("GCGC"), // 100% GC
	}
	got := calculateGCContent(seqs)
	// Only 1 non-empty sequence contributes, but the divisor is len(sequences)=2.
	// totalGC = 1.0 (from GCGC), then / 2 = 0.5
	want := 0.5
	if math.Abs(got-want) > 1e-9 {
		t.Errorf("calculateGCContent with empty-sequence slot = %v, want %v", got, want)
	}
}

// TestDetectAppleSiliconCores_ValidRange verifies the function returns a
// positive value consistent with the machine's actual core count.
func TestDetectAppleSiliconCores_KnownBounds(t *testing.T) {
	got := detectAppleSiliconCores()
	total := runtime.NumCPU()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
	if got > total {
		t.Errorf("detectAppleSiliconCores() = %d, want <= total CPUs %d", got, total)
	}
}

func TestSavePhase3Results_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	savePhase3Results(tmpDir, time.Second*30, 0.085, 0.42, 1000000, 5000, 4, "pass", 1.5)
	path := filepath.Join(tmpDir, "phase3_k8_advanced_results.txt")
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty results file")
	}
}

// ---------- VCF processing helper ----------

func createTestGzippedVCFK8(t *testing.T) string {
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

// ---------- processK8VCFFile tests ----------

func TestProcessK8VCFFile_Basic(t *testing.T) {
	vcfPath := createTestGzippedVCFK8(t)
	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}
	fileSize, variants := processK8VCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 2 {
		t.Errorf("expected 2 variants, got %d", variants)
	}
}

func TestProcessK8VCFFile_NonExistentFile(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}
	fileSize, variants := processK8VCFFile(model, "/nonexistent/path/test.vcf.gz", "TEST")
	if fileSize != 0 {
		t.Errorf("expected 0 file size for missing file, got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for missing file, got %d", variants)
	}
}

func TestProcessK8VCFFile_NotGzipped(t *testing.T) {
	tmpDir := t.TempDir()
	path := filepath.Join(tmpDir, "notgzip.vcf.gz")
	os.WriteFile(path, []byte("this is not gzip data"), 0600)

	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}
	fileSize, variants := processK8VCFFile(model, path, "TEST")
	if fileSize != 0 || variants != 0 {
		t.Errorf("expected 0,0 for invalid gzip, got %d,%d", fileSize, variants)
	}
}

func TestProcessK8VCFFile_HeadersOnly(t *testing.T) {
	tmpDir := t.TempDir()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Close()
	vcfPath := filepath.Join(tmpDir, "headers_only.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}
	fileSize, variants := processK8VCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size (from header), got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for headers-only file, got %d", variants)
	}
}

// ---------- processK8PopulationVCFs tests ----------

func TestProcessK8PopulationVCFs_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "POP")
	os.MkdirAll(popDir, 0750)

	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}
	fileSize, variants := processK8PopulationVCFs(model, "POP", tmpDir, 1)
	if fileSize != 0 {
		t.Errorf("expected 0 bytes for empty dir, got %d", fileSize)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for empty dir, got %d", variants)
	}
}

func TestProcessK8PopulationVCFs_WithVCFs(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	os.MkdirAll(popDir, 0750)

	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Write([]byte("chr1\t100\t.\tA\tT\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t200\t.\tG\tC\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t300\t.\tATGC\tGCAT\t.\tPASS\t.\n"))
	gw.Close()
	vcfPath := filepath.Join(popDir, "chr1.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}
	fileSize, variants := processK8PopulationVCFs(model, "YRI", tmpDir, 1)
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 3 {
		t.Errorf("expected 3 variants, got %d", variants)
	}
}

// ---------- generateValidationSequences tests (k=8) ----------

func TestGenerateValidationSequences_K8_Basic(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
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

func TestGenerateValidationSequences_K8_Multiple(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
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

// ---------- processK8Populations ----------

func TestProcessK8Populations_EmptyDir_SkipsAll(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(8)
	tmpDir := t.TempDir()
	populations := []string{"YRI", "CEU", "CHB"}
	bytes, variants := processK8Populations(model, tmpDir, populations, 1)
	if bytes != 0 {
		t.Errorf("expected 0 bytes, got %d", bytes)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants, got %d", variants)
	}
}

func TestProcessK8Populations_EmptyPopulations_ReturnsZero(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(8)
	bytes, variants := processK8Populations(model, t.TempDir(), []string{}, 1)
	if bytes != 0 || variants != 0 {
		t.Errorf("expected (0, 0) for empty populations, got (%d, %d)", bytes, variants)
	}
}

// TestProcessK8Populations_WithExistingPop exercises the goroutine path where
// the population directory exists, covering atomic adds and the completion printf.
func TestProcessK8Populations_WithExistingPop(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("failed to create pop dir: %v", err)
	}

	// Create a gzipped VCF file inside the population directory
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Write([]byte("chr1\t100\t.\tACGT\tGCAT\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t200\t.\tT\tA\t.\tPASS\t.\n"))
	gw.Close()
	vcfPath := filepath.Join(popDir, "chr1.vcf.gz")
	if err := os.WriteFile(vcfPath, buf.Bytes(), 0600); err != nil {
		t.Fatalf("failed to write vcf: %v", err)
	}

	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}

	totalBytes, totalVariants := processK8Populations(model, tmpDir, []string{"YRI"}, 1)
	if totalBytes <= 0 {
		t.Errorf("expected positive bytes, got %d", totalBytes)
	}
	if totalVariants != 2 {
		t.Errorf("expected 2 variants, got %d", totalVariants)
	}
}

// TestProcessK8Populations_MultiplePops exercises the case with multiple
// populations where some exist and some do not (both branches of the os.Stat check).
func TestProcessK8Populations_MultiplePops_SomeExist(t *testing.T) {
	tmpDir := t.TempDir()

	// Create YRI population with a VCF file
	yriDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(yriDir, 0750); err != nil {
		t.Fatalf("failed to create YRI dir: %v", err)
	}
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Write([]byte("chr1\t100\t.\tA\tT\t.\tPASS\t.\n"))
	gw.Close()
	vcfPath := filepath.Join(yriDir, "chr1.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	// CEU does not exist - should be skipped with warning print
	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}

	totalBytes, totalVariants := processK8Populations(model, tmpDir, []string{"YRI", "CEU"}, 2)
	if totalBytes <= 0 {
		t.Errorf("expected positive bytes from YRI, got %d", totalBytes)
	}
	if totalVariants != 1 {
		t.Errorf("expected 1 variant from YRI, got %d", totalVariants)
	}
}

// TestSavePhase3Results_AllComplexityBranches covers different complexity and
// gcContent combinations to exercise all decision branches in savePhase3Results.
func TestSavePhase3Results_AllComplexityBranches(t *testing.T) {
	tests := []struct {
		name           string
		complexity     float64
		gcContent      float64
		complexityVsK6 float64
		status         string
	}{
		// Branch: complexity > 0.10 && gcContent in [0.40,0.44] -> STRETCH GOAL SUCCESS
		{"stretch_goal", 0.11, 0.42, 1.5, "PHASE 3 STRETCH GOAL SUCCESS"},
		// Branch: complexity > 0.08 && gcContent in [0.40,0.44] -> EXCELLENT SUCCESS
		{"excellent_success", 0.09, 0.42, 1.2, "PHASE 3 EXCELLENT SUCCESS"},
		// Branch: complexityVsK6 > 1.2 -> SIGNIFICANT IMPROVEMENT
		{"significant_improvement", 0.06, 0.42, 1.3, "PHASE 3 SIGNIFICANT IMPROVEMENT"},
		// Branch: else -> MARGINAL BENEFIT
		{"marginal_benefit", 0.04, 0.42, 1.0, "PHASE 3 MARGINAL BENEFIT"},
		// Non-zero but very small values
		{"zero_values", 0.001, 0.10, 0.5, "PHASE 3 NO BENEFIT"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			// duration=1s to avoid division by zero in throughput calc
			savePhase3Results(tmpDir, time.Second, tt.complexity, tt.gcContent,
				500000, 1000, 4, tt.status, tt.complexityVsK6)
			path := filepath.Join(tmpDir, "phase3_k8_advanced_results.txt")
			info, err := os.Stat(path)
			if err != nil {
				t.Fatalf("expected results file to exist: %v", err)
			}
			if info.Size() == 0 {
				t.Error("results file should not be empty")
			}
		})
	}
}

// TestSavePhase3Results_WriteError triggers the os.WriteFile error branch by
// using a directory path as the results file destination.
func TestSavePhase3Results_WriteError(t *testing.T) {
	tmpDir := t.TempDir()
	// Create a subdirectory where the results file path would collide with a directory
	collisionDir := filepath.Join(tmpDir, "phase3_k8_advanced_results.txt")
	if err := os.MkdirAll(collisionDir, 0750); err != nil {
		t.Fatalf("failed to create collision dir: %v", err)
	}
	// This should attempt to write to the path but fail because it's a directory
	// The function prints a warning but does not panic
	savePhase3Results(tmpDir, time.Second, 0.05, 0.42, 100, 10, 4, "test", 1.0)
	// If we reach here without panic, the error branch was exercised
}

// TestProcessK8VCFFile_LargeVCF tests processing a VCF file with many variants
// to trigger the "process accumulated sequence" path (> 50000 bytes).
func TestProcessK8VCFFile_LargeVCF(t *testing.T) {
	tmpDir := t.TempDir()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))

	// Write enough variants so sequenceBuilder.Len() > 50000
	// Each variant adds ~20 bytes; need > 2500 variants
	longRef := string(make([]byte, 50))
	longAlt := string(make([]byte, 50))
	for i := range []byte(longRef) {
		longRef = longRef[:i] + "A" + longRef[i+1:]
	}
	for i := range []byte(longAlt) {
		longAlt = longAlt[:i] + "T" + longAlt[i+1:]
	}
	// Each line contributes 100 bytes to sequenceBuilder; 600 lines = 60000 bytes
	for i := 0; i < 600; i++ {
		line := "chr1\t" + fmt.Sprintf("%d", i*100) + "\t.\t" + longRef + "\t" + longAlt + "\t.\tPASS\t.\n"
		gw.Write([]byte(line))
	}
	gw.Close()
	vcfPath := filepath.Join(tmpDir, "large.vcf.gz")
	if err := os.WriteFile(vcfPath, buf.Bytes(), 0600); err != nil {
		t.Fatalf("failed to write vcf: %v", err)
	}

	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}
	fileSize, variants := processK8VCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	if variants != 600 {
		t.Errorf("expected 600 variants, got %d", variants)
	}
}

// TestProcessK8VCFFile_ShortFieldsSkipped tests that VCF lines with fewer
// than 5 tab-separated fields are skipped (the len(fields) < 5 branch).
func TestProcessK8VCFFile_ShortFieldsSkipped(t *testing.T) {
	tmpDir := t.TempDir()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	// Lines with < 5 fields should be skipped
	gw.Write([]byte("chr1\t100\t.\n"))          // 3 fields
	gw.Write([]byte("chr1\t200\t.\tA\n"))        // 4 fields
	gw.Write([]byte("chr1\t300\t.\tA\tT\t.\n")) // 6 fields - valid
	gw.Close()
	vcfPath := filepath.Join(tmpDir, "short_fields.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	model := authenticity.NewAdaptiveHMMModel(8)
	fileSize, variants := processK8VCFFile(model, vcfPath, "TEST")
	if fileSize <= 0 {
		t.Errorf("expected positive file size, got %d", fileSize)
	}
	// Only the 6-field line is valid
	if variants != 1 {
		t.Errorf("expected 1 variant (only valid line), got %d", variants)
	}
}

// TestGenerateValidationSequences_K8_UntrainedModel tests that generateValidationSequences
// returns an error when the model has not been trained (covering the error path at line 393).
func TestGenerateValidationSequences_K8_UntrainedModel(t *testing.T) {
	// Create a fresh, untrained k=8 model
	model := authenticity.NewAdaptiveHMMModel(8)
	if model == nil {
		t.Fatal("failed to create k=8 model")
	}
	// GenerateSequence on an untrained model returns an error
	_, err := generateValidationSequences(model, 1, 100)
	if err == nil {
		t.Error("expected error from generateValidationSequences with untrained model, got nil")
	}
}

// TestDetectAppleSiliconCores_K8_OutputConsistency verifies that repeated calls
// return the same value and the result is sensible for the current machine.
func TestDetectAppleSiliconCores_K8_OutputConsistency(t *testing.T) {
	first := detectAppleSiliconCores()
	second := detectAppleSiliconCores()
	if first != second {
		t.Errorf("detectAppleSiliconCores() not deterministic: %d vs %d", first, second)
	}
	total := runtime.NumCPU()
	if first > total {
		t.Errorf("detectAppleSiliconCores() = %d > runtime.NumCPU() = %d", first, total)
	}
	if first <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", first)
	}
}
