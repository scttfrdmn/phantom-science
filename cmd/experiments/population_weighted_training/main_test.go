//go:build !release
// +build !release

package main

import (
	"compress/gzip"
	"math"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestAbs_PopWeight(t *testing.T) {
	tests := []struct {
		name     string
		input    float64
		expected float64
	}{
		{"positive", 0.5, 0.5},
		{"negative", -0.5, 0.5},
		{"zero", 0.0, 0.0},
		{"large positive", 99.9, 99.9},
		{"large negative", -99.9, 99.9},
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

func TestCalculateStateSpaceSize_Pop(t *testing.T) {
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
	}
	for _, tt := range tests {
		got := calculateStateSpaceSize(tt.k)
		if got != tt.want {
			t.Errorf("calculateStateSpaceSize(%d) = %d, want %d", tt.k, got, tt.want)
		}
	}
}

func TestInitializePopulationWeights(t *testing.T) {
	weights := initializePopulationWeights()

	// African populations should have highest weights
	if weights.African["YRI"] <= 1.0 {
		t.Errorf("YRI weight = %v, want > 1.0 (African high diversity)", weights.African["YRI"])
	}
	if weights.African["LWK"] <= 1.0 {
		t.Errorf("LWK weight = %v, want > 1.0", weights.African["LWK"])
	}
	if weights.African["GWD"] <= 1.0 {
		t.Errorf("GWD weight = %v, want > 1.0", weights.African["GWD"])
	}

	// European baseline should be 1.0
	if weights.European["CEU"] != 1.0 {
		t.Errorf("CEU weight = %v, want 1.0 (European baseline)", weights.European["CEU"])
	}

	// Asian populations should be defined
	if _, ok := weights.Asian["CHB"]; !ok {
		t.Error("CHB not found in Asian weights")
	}
	if _, ok := weights.Asian["JPT"]; !ok {
		t.Error("JPT not found in Asian weights")
	}

	// Admixed populations should have high weights
	if weights.Admixed["MXL"] <= 1.0 {
		t.Errorf("MXL weight = %v, want > 1.0 (high admixture)", weights.Admixed["MXL"])
	}
}

func TestGetPopulationWeight(t *testing.T) {
	weights := initializePopulationWeights()

	tests := []struct {
		pop     string
		wantMin float64
		wantMax float64
	}{
		{"YRI", 1.5, 3.0},   // African highest
		{"CEU", 0.9, 1.1},   // European baseline
		{"CHB", 1.0, 2.0},   // Asian moderate
		{"MXL", 1.5, 3.0},   // Admixed high
		{"UNKNOWN", 0.9, 1.1}, // Unknown gets default 1.0
	}
	for _, tt := range tests {
		t.Run(tt.pop, func(t *testing.T) {
			got := getPopulationWeight(tt.pop, weights)
			if got < tt.wantMin || got > tt.wantMax {
				t.Errorf("getPopulationWeight(%q) = %v, want [%v, %v]", tt.pop, got, tt.wantMin, tt.wantMax)
			}
		})
	}
}

func TestGetPopulationWeightDefaultIsOne(t *testing.T) {
	weights := initializePopulationWeights()
	got := getPopulationWeight("NONEXISTENT_POP", weights)
	if got != 1.0 {
		t.Errorf("getPopulationWeight(unknown) = %v, want 1.0 (default)", got)
	}
}

func TestCalculateComplexity_Pop(t *testing.T) {
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
			name:    "short sequences below k=6",
			seqs:    [][]byte{[]byte("ACGT")},
			wantMin: 0.0,
			wantMax: 0.0,
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

func TestCalculateGCContent_Pop(t *testing.T) {
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

func TestFormatBytes_Pop(t *testing.T) {
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

func TestDetectAppleSiliconCores_Pop(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestGetWeightingEffectiveness(t *testing.T) {
	tests := []struct {
		name     string
		current  float64
		previous float64
	}{
		{"improved", 0.09, 0.08},
		{"maintained", 0.080, 0.082},
		{"decreased", 0.07, 0.08},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getWeightingEffectiveness(tt.current, tt.previous)
			if len(got) == 0 {
				t.Errorf("getWeightingEffectiveness(%v, %v) returned empty string", tt.current, tt.previous)
			}
		})
	}
}

func TestGetNextStepsPhase22(t *testing.T) {
	tests := []struct {
		name      string
		complexity float64
		gcContent  float64
	}{
		{"complete success", 0.09, 0.42},
		{"good progress", 0.06, 0.42},
		{"needs advanced methods", 0.03, 0.50},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getNextStepsPhase22(tt.complexity, tt.gcContent)
			if len(got) == 0 {
				t.Errorf("getNextStepsPhase22(%v, %v) returned empty string", tt.complexity, tt.gcContent)
			}
		})
	}
}

func TestGetSuccessEmojiPhase22(t *testing.T) {
	tests := []struct {
		name      string
		complexity float64
		gcContent  float64
	}{
		{"excellent result", 0.09, 0.42},
		{"good result", 0.06, 0.42},
		{"poor result", 0.02, 0.55},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := getSuccessEmojiPhase22(tt.complexity, tt.gcContent)
			if len(got) == 0 {
				t.Errorf("getSuccessEmojiPhase22(%v, %v) returned empty string", tt.complexity, tt.gcContent)
			}
		})
	}
}

func TestAfricanWeightsHighestDiversity(t *testing.T) {
	weights := initializePopulationWeights()
	europeanBaseline := weights.European["CEU"]
	for pop, weight := range weights.African {
		if weight <= europeanBaseline {
			t.Errorf("African population %s weight %v should be > European baseline %v", pop, weight, europeanBaseline)
		}
	}
}

func TestSavePhase22Results_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	weights := PopulationWeights{
		African:  map[string]float64{"YRI": 2.0},
		European: map[string]float64{"GBR": 1.0},
		Asian:    map[string]float64{"CHS": 1.0},
		Admixed:  map[string]float64{"MXL": 1.5},
	}
	savePhase22Results(tmpDir, time.Second*30, 0.078, 0.42, 1000000, 5000, 4, "pass", weights)
	path := filepath.Join(tmpDir, "phase2_2_population_weighted_results.txt")
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty results file")
	}
}

// ---------- processWeightedPopulations ----------

// TestProcessWeightedPopulations_EmptyDir calls processWeightedPopulations with a
// tmpDir that has no population subdirectories. All populations are skipped and
// the function returns (0, 0) without error.
func TestProcessWeightedPopulations_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	model := authenticity.NewAdaptiveHMMModel(4)
	weights := initializePopulationWeights()

	populations := []string{"YRI", "CEU", "CHB"}
	bytes, variants := processWeightedPopulations(model, tmpDir, populations, weights, 2)
	if bytes != 0 {
		t.Errorf("expected 0 bytes, got %d", bytes)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants, got %d", variants)
	}
}

// TestProcessWeightedPopulations_EmptyPopDirs creates population directories
// that exist but have no VCF files. The function processes them but returns zeros.
func TestProcessWeightedPopulations_EmptyPopDirs(t *testing.T) {
	tmpDir := t.TempDir()
	model := authenticity.NewAdaptiveHMMModel(4)
	weights := initializePopulationWeights()

	for _, pop := range []string{"YRI", "CEU"} {
		if err := os.MkdirAll(filepath.Join(tmpDir, pop), 0750); err != nil {
			t.Fatalf("MkdirAll: %v", err)
		}
	}

	bytes, variants := processWeightedPopulations(model, tmpDir, []string{"YRI", "CEU"}, weights, 2)
	if bytes != 0 {
		t.Errorf("expected 0 bytes from empty pop dirs, got %d", bytes)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants from empty pop dirs, got %d", variants)
	}
}

// ---------- processWeightedPopulationVCFs ----------

// TestProcessWeightedPopulationVCFs_EmptyDir verifies that processWeightedPopulationVCFs
// returns (0, 0) when the population directory exists but has no VCF files.
func TestProcessWeightedPopulationVCFs_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	if err := os.MkdirAll(filepath.Join(tmpDir, "YRI"), 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}
	model := authenticity.NewAdaptiveHMMModel(4)

	bytes, variants := processWeightedPopulationVCFs(model, "YRI", tmpDir, 2.0, 2)
	if bytes != 0 {
		t.Errorf("expected 0 bytes, got %d", bytes)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants, got %d", variants)
	}
}

// ---------- processWeightedVCFFile ----------

// TestProcessWeightedVCFFile_WithGzipVCF creates a minimal gzipped VCF file and
// calls processWeightedVCFFile to exercise the gzip reading, VCF parsing, and
// model training code paths.
func TestProcessWeightedVCFFile_WithGzipVCF(t *testing.T) {
	tmpDir := t.TempDir()
	vcfPath := createWeightedTestVCFGz(t, tmpDir)

	model := authenticity.NewAdaptiveHMMModel(4)
	// Train the model with a seed sequence first so it can handle training calls
	if err := model.TrainFromSequence([]byte("ACGTACGTACGTACGT")); err != nil {
		t.Fatalf("initial TrainFromSequence: %v", err)
	}

	bytes, variants := processWeightedVCFFile(model, vcfPath, "YRI", 2.0)
	if bytes == 0 {
		t.Error("expected non-zero bytes from a real gzipped VCF file")
	}
	if variants == 0 {
		t.Error("expected non-zero variant count from VCF with data lines")
	}
}

// TestProcessWeightedVCFFile_NonExistentFile verifies graceful error handling
// when the VCF file does not exist.
func TestProcessWeightedVCFFile_NonExistentFile(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	bytes, variants := processWeightedVCFFile(model, "/no/such/file.vcf.gz", "CEU", 1.0)
	if bytes != 0 {
		t.Errorf("expected 0 bytes for missing file, got %d", bytes)
	}
	if variants != 0 {
		t.Errorf("expected 0 variants for missing file, got %d", variants)
	}
}

// createWeightedTestVCFGz is a helper that creates a minimal gzipped VCF file
// for use in processWeightedVCFFile tests.
func createWeightedTestVCFGz(t *testing.T, dir string) string {
	t.Helper()
	vcfPath := filepath.Join(dir, "test.vcf.gz")
	f, err := os.Create(vcfPath)
	if err != nil {
		t.Fatalf("create vcf.gz: %v", err)
	}
	defer f.Close()

	gw := gzip.NewWriter(f)
	_, _ = gw.Write([]byte("##fileformat=VCFv4.1\n"))
	_, _ = gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	_, _ = gw.Write([]byte("chr1\t1000\t.\tACGT\tGCAT\t.\tPASS\t.\n"))
	_, _ = gw.Write([]byte("chr1\t2000\t.\tA\tG\t.\tPASS\t.\n"))
	if err := gw.Close(); err != nil {
		t.Fatalf("gzip.Close: %v", err)
	}
	return vcfPath
}

// ---------- generateValidationSequences ----------

// TestGenerateValidationSequences_PopWeighted trains a model and generates a
// small set of validation sequences.
func TestGenerateValidationSequences_PopWeighted(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)

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

// TestGenerateValidationSequences_PopWeighted_Zero checks that requesting zero
// sequences returns an empty slice without error.
func TestGenerateValidationSequences_PopWeighted_Zero(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)

	seqs, err := generateValidationSequences(model, 0, 50)
	if err != nil {
		t.Fatalf("generateValidationSequences(0) returned error: %v", err)
	}
	if len(seqs) != 0 {
		t.Errorf("expected 0 sequences, got %d", len(seqs))
	}
}

// TestSavePhase22Results_ComplexityBranches calls savePhase22Results with values
// that trigger different branches in the helper functions called inside it:
// getWeightingEffectiveness (improved/maintained/decreased) and
// getNextStepsPhase22 / getSuccessEmojiPhase22 (3 branches each).
func TestSavePhase22Results_ComplexityBranches(t *testing.T) {
	baseWeights := PopulationWeights{
		African:  map[string]float64{"YRI": 2.0},
		European: map[string]float64{"CEU": 1.0},
		Asian:    map[string]float64{"CHB": 1.2},
		Admixed:  map[string]float64{"MXL": 1.8},
	}

	tests := []struct {
		name          string
		duration      time.Duration
		complexity    float64
		gcContent     float64
		totalBytes    int64
		totalVariants int64
		perfCores     int
		status        string
	}{
		// complexity > 0.08 + gcContent in [0.40, 0.44]: PHASE 2.2 COMPLETE SUCCESS branch
		// getWeightingEffectiveness: complexity(0.09) > previousComplexity(0.0819) → "improved" branch
		// getNextStepsPhase22: complexity>0.08 + gc in range → branch 1
		// getSuccessEmojiPhase22: complexity>0.08 + gc in range → "✅"
		{"complete_success", 30 * time.Second, 0.09, 0.42, 5000000, 10000, 8, "PHASE 2.2 COMPLETE SUCCESS"},

		// complexity > 0.05 + gcContent in [0.38, 0.46]: GOOD PROGRESS branch
		// getWeightingEffectiveness: complexity(0.06) < 0.0819*0.95? No: 0.06 < 0.0779 → "decreased"
		// getNextStepsPhase22: complexity>0.05 + gc in range → branch 2
		// getSuccessEmojiPhase22: complexity>0.05 + gc in range → "🟡"
		{"good_progress", 60 * time.Second, 0.06, 0.43, 2000000, 4000, 6, "PHASE 2.2 GOOD PROGRESS"},

		// complexity <= 0.05 or gcContent out of range: NEEDS REFINEMENT branch
		// getWeightingEffectiveness: complexity(0.04) < 0.0819*0.95=0.0778 → "decreased" branch
		// getNextStepsPhase22: below thresholds → branch 3
		// getSuccessEmojiPhase22: below thresholds → "❌"
		{"needs_refinement", 120 * time.Second, 0.04, 0.55, 1000000, 2000, 4, "PHASE 2.2 NEEDS REFINEMENT"},

		// getWeightingEffectiveness "maintained" branch: currentComplexity in [previous*0.95, previous]
		// 0.0819 * 0.95 = 0.077805 → use 0.079 which is < 0.0819 but >= 0.077805
		{"maintained_complexity", 45 * time.Second, 0.079, 0.42, 3000000, 6000, 8, "PHASE 2.2 GOOD PROGRESS"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			savePhase22Results(tmpDir, tt.duration, tt.complexity, tt.gcContent,
				tt.totalBytes, tt.totalVariants, tt.perfCores, tt.status, baseWeights)
			path := filepath.Join(tmpDir, "phase2_2_population_weighted_results.txt")
			info, err := os.Stat(path)
			if err != nil {
				t.Fatalf("file should exist after savePhase22Results: %v", err)
			}
			if info.Size() == 0 {
				t.Error("results file should not be empty")
			}
		})
	}
}
