//go:build !release
// +build !release

package main

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestIsValidDNASequence(t *testing.T) {
	tests := []struct {
		input    string
		expected bool
	}{
		{"ACGT", true},
		{"acgt", true},    // lowercase accepted via ToUpper
		{"AAAA", true},
		{"GCGCGC", true},
		{"ACGN", false},   // N is not valid
		{"ACGT.", false},  // dot is not valid
		{"", false},        // empty is invalid
		{"A", true},
		{"G", true},
		{"C", true},
		{"T", true},
		{"X", false},
		{"ACGTACGT", true},
	}
	for _, tt := range tests {
		got := isValidDNASequence(tt.input)
		if got != tt.expected {
			t.Errorf("isValidDNASequence(%q) = %v, want %v", tt.input, got, tt.expected)
		}
	}
}

func TestCalculateAverageGCContent(t *testing.T) {
	t.Run("empty returns 0", func(t *testing.T) {
		if got := calculateAverageGCContent(nil); got != 0.0 {
			t.Errorf("expected 0.0 for nil, got %f", got)
		}
		if got := calculateAverageGCContent([][]byte{}); got != 0.0 {
			t.Errorf("expected 0.0 for empty, got %f", got)
		}
	})

	t.Run("all G/C returns 1.0", func(t *testing.T) {
		seq := []byte("GCGCGCGCGC")
		got := calculateAverageGCContent([][]byte{seq})
		if got != 1.0 {
			t.Errorf("expected 1.0, got %f", got)
		}
	})

	t.Run("all A/T returns 0.0", func(t *testing.T) {
		seq := []byte("ATATATATAT")
		got := calculateAverageGCContent([][]byte{seq})
		if got != 0.0 {
			t.Errorf("expected 0.0, got %f", got)
		}
	})

	t.Run("50pct GC returns 0.5", func(t *testing.T) {
		seq := []byte("AACCGGTT")
		got := calculateAverageGCContent([][]byte{seq})
		if got != 0.5 {
			t.Errorf("expected 0.5, got %f", got)
		}
	})

	t.Run("multiple sequences averaged correctly", func(t *testing.T) {
		seqGC := []byte("GCGC") // 100%
		seqAT := []byte("ATAT") // 0%
		got := calculateAverageGCContent([][]byte{seqGC, seqAT})
		if got != 0.5 {
			t.Errorf("expected 0.5, got %f", got)
		}
	})
}

func TestCalculateSequenceComplexity(t *testing.T) {
	t.Run("short sequence returns 0", func(t *testing.T) {
		seq := []byte("ACG")
		if got := calculateSequenceComplexity(seq, 6); got != 0.0 {
			t.Errorf("expected 0.0 for short sequence, got %f", got)
		}
	})

	t.Run("uniform sequence has near-zero complexity", func(t *testing.T) {
		seq := make([]byte, 100)
		for i := range seq {
			seq[i] = 'A'
		}
		got := calculateSequenceComplexity(seq, 4)
		if got > 0.1 {
			t.Errorf("uniform sequence complexity too high: %f", got)
		}
	})

	t.Run("bounded between 0 and 1", func(t *testing.T) {
		seq := []byte("ACGTACGTGCATGCAT")
		got := calculateSequenceComplexity(seq, 4)
		if got < 0 || got > 1 {
			t.Errorf("complexity %f out of [0,1]", got)
		}
	})

	t.Run("deterministic results", func(t *testing.T) {
		seq := []byte("ACGTACGTGCATGCATACGT")
		c1 := calculateSequenceComplexity(seq, 4)
		c2 := calculateSequenceComplexity(seq, 4)
		if c1 != c2 {
			t.Errorf("complexity not deterministic: %f != %f", c1, c2)
		}
	})
}

func TestCalculateAverageComplexity(t *testing.T) {
	t.Run("empty returns 0", func(t *testing.T) {
		if got := calculateAverageComplexity(nil, 4); got != 0.0 {
			t.Errorf("expected 0.0 for nil, got %f", got)
		}
		if got := calculateAverageComplexity([][]byte{}, 4); got != 0.0 {
			t.Errorf("expected 0.0 for empty, got %f", got)
		}
	})

	t.Run("sequences shorter than k excluded", func(t *testing.T) {
		short := []byte("AC")
		got := calculateAverageComplexity([][]byte{short}, 6)
		if got != 0.0 {
			t.Errorf("expected 0.0 when all sequences too short, got %f", got)
		}
	})

	t.Run("positive complexity for valid sequences", func(t *testing.T) {
		seq := []byte("ACGTACGTGCATGCAT")
		got := calculateAverageComplexity([][]byte{seq}, 4)
		if got <= 0 {
			t.Errorf("expected positive complexity, got %f", got)
		}
	})
}

func TestDetectAppleSiliconCoresPositive(t *testing.T) {
	if cores := detectAppleSiliconCores(); cores <= 0 {
		t.Errorf("expected positive core count, got %d", cores)
	}
}

func TestGetParallelPopulationTargets(t *testing.T) {
	pops := getParallelPopulationTargets()
	if len(pops) == 0 {
		t.Error("expected non-empty population targets")
	}
	for _, p := range pops {
		if p.Population == "" {
			t.Error("expected non-empty population name")
		}
		if p.Region == "" {
			t.Error("expected non-empty region")
		}
		if p.Individuals <= 0 {
			t.Errorf("population %s has non-positive individual count %d", p.Population, p.Individuals)
		}
	}
}

func TestParallelVCFTrainingStruct(t *testing.T) {
	training := ParallelVCFTraining{
		Population:  "YRI",
		Region:      "Africa",
		Individuals: 108,
	}
	if training.Population != "YRI" {
		t.Errorf("expected YRI, got %s", training.Population)
	}
	if training.Region != "Africa" {
		t.Errorf("expected Africa, got %s", training.Region)
	}
	if training.Individuals != 108 {
		t.Errorf("expected 108, got %d", training.Individuals)
	}
	if training.Trained {
		t.Error("expected Trained to default to false")
	}
}

// --- Additional tests added to increase coverage ---

func TestIsValidDNASequence_AdditionalCases(t *testing.T) {
	tests := []struct {
		input    string
		expected bool
	}{
		// Boundary: single valid bases
		{"A", true},
		{"T", true},
		{"G", true},
		{"C", true},
		// Mixed case accepted via ToUpper
		{"AcGt", true},
		{"aAtTgGcC", true},
		// Invalid characters
		{"ACGTU", false}, // U (uracil) is not valid DNA
		{"ACG-", false},  // dash not valid
		{"AC GT", false}, // space not valid
		{"ACGT1", false}, // number not valid
		// Empty string is invalid
		{"", false},
		// Very long valid sequence
		{"ACGTACGTACGTACGT", true},
		// Lowercase valid only bases
		{"acgt", true},
	}
	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			got := isValidDNASequence(tt.input)
			if got != tt.expected {
				t.Errorf("isValidDNASequence(%q) = %v, want %v", tt.input, got, tt.expected)
			}
		})
	}
}

func TestCalculateAverageGCContent_EmptyBytesSlice(t *testing.T) {
	// Slice with one empty []byte — totalBases becomes 0, hits the totalBases==0 branch
	got := calculateAverageGCContent([][]byte{{}})
	if got != 0.0 {
		t.Errorf("calculateAverageGCContent([empty]) = %f, want 0.0", got)
	}
}

func TestCalculateAverageGCContent_MultipleEmptySequences(t *testing.T) {
	got := calculateAverageGCContent([][]byte{{}, {}, {}})
	if got != 0.0 {
		t.Errorf("calculateAverageGCContent(multiple empty) = %f, want 0.0", got)
	}
}

func TestCalculateAverageGCContent_MixedEmptyAndValid(t *testing.T) {
	// One empty, one full-GC — result should still reflect only real bases
	got := calculateAverageGCContent([][]byte{{}, []byte("GGGG")})
	if got != 1.0 {
		t.Errorf("calculateAverageGCContent([empty, GGGG]) = %f, want 1.0", got)
	}
}

func TestCalculateAverageGCContent_SingleBase(t *testing.T) {
	tests := []struct {
		name string
		seq  []byte
		want float64
	}{
		{"single G", []byte("G"), 1.0},
		{"single C", []byte("C"), 1.0},
		{"single A", []byte("A"), 0.0},
		{"single T", []byte("T"), 0.0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateAverageGCContent([][]byte{tt.seq})
			if got != tt.want {
				t.Errorf("calculateAverageGCContent([%q]) = %f, want %f", tt.seq, got, tt.want)
			}
		})
	}
}

func TestCalculateSequenceComplexity_PVT_Additional(t *testing.T) {
	tests := []struct {
		name    string
		seq     []byte
		k       int
		wantMin float64
		wantMax float64
	}{
		// k equals length: single unique kmer, totalKmers = 1, ratio = 1.0
		{"k equals length single kmer", []byte("ACGT"), 4, 1.0, 1.0},
		// k = 1, all same base: 1 unique out of n → very low
		{"k=1 all A low complexity", []byte("AAAAAAAAAA"), 1, 0.0, 0.15},
		// k = 1, all four bases: 4 unique out of 4 → 1.0
		{"k=1 all four bases", []byte("ACGT"), 1, 1.0, 1.0},
		// Longer sequence with repeated pattern
		{"repeated ACGT pattern", []byte("ACGTACGTACGT"), 4, 0.0, 1.0},
		// Larger k than sequence: returns 0
		{"k > len returns 0", []byte("ACG"), 10, 0.0, 0.0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateSequenceComplexity(tt.seq, tt.k)
			if got < tt.wantMin || got > tt.wantMax {
				t.Errorf("calculateSequenceComplexity(%q, %d) = %v, want [%v, %v]",
					tt.seq, tt.k, got, tt.wantMin, tt.wantMax)
			}
		})
	}
}

func TestCalculateAverageComplexity_AdditionalCases(t *testing.T) {
	tests := []struct {
		name    string
		seqs    [][]byte
		k       int
		wantMin float64
		wantMax float64
	}{
		{
			name:    "multiple valid sequences average is bounded",
			seqs:    [][]byte{[]byte("ACGTACGT"), []byte("GCTAGCTA"), []byte("TTTAAACCC")},
			k:       4,
			wantMin: 0.0,
			wantMax: 1.0,
		},
		{
			name:    "mix of short and valid sequences only counts valid",
			seqs:    [][]byte{[]byte("AC"), []byte("ACGTACGT")},
			k:       4,
			wantMin: 0.0,
			wantMax: 1.0,
		},
		{
			name:    "all sequences shorter than k returns 0",
			seqs:    [][]byte{[]byte("A"), []byte("AC"), []byte("ACG")},
			k:       4,
			wantMin: 0.0,
			wantMax: 0.0,
		},
		{
			name:    "uniform high-complexity sequences all equal",
			seqs:    [][]byte{[]byte("ACGTACGTACGT"), []byte("ACGTACGTACGT")},
			k:       4,
			wantMin: 0.0,
			wantMax: 1.0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateAverageComplexity(tt.seqs, tt.k)
			if got < tt.wantMin || got > tt.wantMax {
				t.Errorf("calculateAverageComplexity(%v, %d) = %v, want [%v, %v]",
					tt.seqs, tt.k, got, tt.wantMin, tt.wantMax)
			}
		})
	}
}

func TestGetParallelPopulationTargets_Details(t *testing.T) {
	pops := getParallelPopulationTargets()
	if len(pops) == 0 {
		t.Fatal("expected non-empty population targets")
	}

	// Check all required regions are represented
	regions := make(map[string]bool)
	for _, p := range pops {
		regions[p.Region] = true
	}
	for _, region := range []string{"Africa", "Europe", "East Asia"} {
		if !regions[region] {
			t.Errorf("expected region %q to be present in population targets", region)
		}
	}

	// Verify total individual count is reasonable (> 0 overall)
	total := 0
	for _, p := range pops {
		total += p.Individuals
	}
	if total <= 0 {
		t.Errorf("total individuals across all populations = %d, want > 0", total)
	}
}

func TestDetectAppleSiliconCores_PVT_ResultBounded(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
	// Result must never exceed total CPU count (it equals or is a fraction of it)
	if got > 1024 {
		t.Errorf("detectAppleSiliconCores() = %d, unreasonably large", got)
	}
}

func TestParallelVCFTraining_ZeroValueStruct(t *testing.T) {
	var training ParallelVCFTraining
	if training.Population != "" {
		t.Errorf("expected empty Population, got %q", training.Population)
	}
	if training.Individuals != 0 {
		t.Errorf("expected 0 Individuals, got %d", training.Individuals)
	}
	if training.Trained {
		t.Error("expected Trained to be false for zero value")
	}
	if training.ProcessedVariants != 0 {
		t.Errorf("expected 0 ProcessedVariants, got %d", training.ProcessedVariants)
	}
}

func TestCalculateAverageGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{
		[]byte("GCGCGC"),
		[]byte("CCCCGGGG"),
	}
	got := calculateAverageGCContent(seqs)
	if got != 1.0 {
		t.Errorf("calculateAverageGCContent(all-GC seqs) = %f, want 1.0", got)
	}
}

func TestCalculateAverageGCContent_AllAT(t *testing.T) {
	seqs := [][]byte{
		[]byte("AAAAT"),
		[]byte("TTTTT"),
	}
	got := calculateAverageGCContent(seqs)
	if got != 0.0 {
		t.Errorf("calculateAverageGCContent(all-AT seqs) = %f, want 0.0", got)
	}
}

// ---------- processVCFFile ----------

func createTestGzippedVCF(t *testing.T) string {
	t.Helper()
	tmpDir := t.TempDir()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	// VCF header lines
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("##reference=GRCh38\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	// Variant lines with valid DNA ref/alt
	gw.Write([]byte("chr1\t1000\t.\tA\tG\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t2000\t.\tCT\tC\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t3000\t.\tGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC\tA\t.\tPASS\t.\n"))
	// Invalid line (too few fields)
	gw.Write([]byte("chr1\t4000\t.\n"))
	// Line with non-DNA characters (should be skipped)
	gw.Write([]byte("chr1\t5000\t.\tN\tG\t.\tPASS\t.\n"))
	gw.Close()
	vcfPath := filepath.Join(tmpDir, "chr1.vcf.gz")
	if err := os.WriteFile(vcfPath, buf.Bytes(), 0600); err != nil {
		t.Fatalf("failed to write test VCF: %v", err)
	}
	return vcfPath
}

func TestProcessVCFFile_WithValidData(t *testing.T) {
	vcfPath := createTestGzippedVCF(t)
	variants, seqs, err := processVCFFile(vcfPath, "YRI")
	if err != nil {
		t.Fatalf("processVCFFile returned error: %v", err)
	}
	if variants <= 0 {
		t.Errorf("expected > 0 variants, got %d", variants)
	}
	_ = seqs
}

func TestProcessVCFFile_NotFound(t *testing.T) {
	variants, seqs, err := processVCFFile("/nonexistent/file.vcf.gz", "YRI")
	if err == nil {
		t.Error("expected error for nonexistent file")
	}
	if variants != 0 || seqs != nil {
		t.Errorf("expected zero variants and nil seqs for error case")
	}
}

func TestProcessVCFFile_NotGzip(t *testing.T) {
	tmpDir := t.TempDir()
	plainFile := filepath.Join(tmpDir, "test.vcf.gz")
	os.WriteFile(plainFile, []byte("not gzip content"), 0600)
	_, _, err := processVCFFile(plainFile, "YRI")
	if err == nil {
		t.Error("expected error for non-gzip file")
	}
}

// ---------- findPopulationVCFFiles ----------

func TestFindPopulationVCFFiles_Empty(t *testing.T) {
	tmpDir := t.TempDir()
	files, err := findPopulationVCFFiles(tmpDir, "YRI")
	if err != nil {
		t.Fatalf("findPopulationVCFFiles returned error: %v", err)
	}
	if len(files) != 0 {
		t.Errorf("expected 0 files in empty dir, got %d", len(files))
	}
}

func TestFindPopulationVCFFiles_WithPopDir(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	os.MkdirAll(popDir, 0755)
	// Create a VCF file matching the chr*.vcf.gz pattern
	os.WriteFile(filepath.Join(popDir, "chr1.vcf.gz"), []byte("gzip"), 0600)
	files, err := findPopulationVCFFiles(tmpDir, "YRI")
	if err != nil {
		t.Fatalf("findPopulationVCFFiles returned error: %v", err)
	}
	if len(files) == 0 {
		t.Error("expected at least 1 file in populated dir")
	}
}

func TestFindPopulationVCFFiles_WithPrefixFile(t *testing.T) {
	tmpDir := t.TempDir()
	// Create a VCF file matching basePath/YRI_*.vcf.gz pattern
	os.WriteFile(filepath.Join(tmpDir, "YRI_chr1.vcf.gz"), []byte("gzip"), 0600)
	files, err := findPopulationVCFFiles(tmpDir, "YRI")
	if err != nil {
		t.Fatalf("findPopulationVCFFiles returned error: %v", err)
	}
	if len(files) == 0 {
		t.Error("expected at least 1 file from prefix pattern")
	}
}

// ---------- checkVCFDataAvailability ----------

func TestCheckVCFDataAvailability_NoPaths(t *testing.T) {
	// None of the hardcoded paths exist in test environment
	got := checkVCFDataAvailability()
	// Should return "" since no real genomics data is available
	if got != "" {
		t.Logf("checkVCFDataAvailability returned %q (data may be available)", got)
	}
}

// ---------- writeParallelVCFFASTA ----------

func TestWriteParallelVCFFASTA_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	training := &ParallelVCFTraining{
		Population:  "YRI",
		Region:      "Africa",
		Individuals: 108,
	}
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	filePath := filepath.Join(tmpDir, "YRI_PARALLEL_VCF_L82_S01.fasta")
	if err := writeParallelVCFFASTA(filePath, "YRI", training, len(seq), 1, seq); err != nil {
		t.Fatalf("writeParallelVCFFASTA returned error: %v", err)
	}
	info, err := os.Stat(filePath)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty FASTA file")
	}
}

func TestWriteParallelVCFFASTA_BadPath(t *testing.T) {
	training := &ParallelVCFTraining{Population: "YRI", Region: "Africa"}
	err := writeParallelVCFFASTA("/nonexistent/dir/file.fasta", "YRI", training, 10, 1, []byte("ATCG"))
	if err == nil {
		t.Error("expected error for bad path")
	}
}

// ---------- analyzeParallelVCFTrainingResults ----------

func TestAnalyzeParallelVCFTrainingResults_Empty(t *testing.T) {
	// Should not panic with empty map
	analyzeParallelVCFTrainingResults(map[string]*ParallelVCFTraining{}, 0)
}

func TestAnalyzeParallelVCFTrainingResults_AllSucceeded(t *testing.T) {
	models := map[string]*ParallelVCFTraining{
		"YRI": {Population: "YRI", Region: "Africa", Trained: true, ProcessedVariants: 5000, ActualComplexity: 0.08, ActualGCContent: 0.47},
		"CEU": {Population: "CEU", Region: "Europe", Trained: true, ProcessedVariants: 4000, ActualComplexity: 0.07, ActualGCContent: 0.41},
	}
	analyzeParallelVCFTrainingResults(models, 30*time.Second)
}

func TestAnalyzeParallelVCFTrainingResults_SomeFailed(t *testing.T) {
	models := map[string]*ParallelVCFTraining{
		"YRI": {Population: "YRI", Region: "Africa", Trained: true, ProcessedVariants: 5000, ActualComplexity: 0.08, ActualGCContent: 0.47},
		"CEU": {Population: "CEU", Region: "Europe", Trained: false, Error: fmt.Errorf("test error")},
	}
	analyzeParallelVCFTrainingResults(models, 30*time.Second)
}

func TestAnalyzeParallelVCFTrainingResults_MostlyFailed(t *testing.T) {
	models := map[string]*ParallelVCFTraining{
		"YRI": {Population: "YRI", Region: "Africa", Trained: false, Error: fmt.Errorf("failed")},
		"CEU": {Population: "CEU", Region: "Europe", Trained: false, Error: fmt.Errorf("failed")},
		"CHB": {Population: "CHB", Region: "East Asia", Trained: false, Error: fmt.Errorf("failed")},
		"MXL": {Population: "MXL", Region: "Americas", Trained: true, ProcessedVariants: 100},
	}
	analyzeParallelVCFTrainingResults(models, 10*time.Second)
}

// ---------- helpers ----------

// newTestAdaptiveModel creates an AdaptiveHMMModel for tests, skipping the test
// if the model cannot be created (e.g., invalid k).
func newTestAdaptiveModel(t *testing.T, k int) *authenticity.AdaptiveHMMModel {
	t.Helper()
	m := authenticity.NewAdaptiveHMMModel(k)
	if m == nil {
		t.Fatalf("NewAdaptiveHMMModel(%d) returned nil", k)
	}
	return m
}

// ---------- trainPopulationsInParallel ----------

// TestTrainPopulationsInParallel_EmptyPopulationsSlice verifies that calling with
// an empty populations slice returns an empty (non-nil) map without panicking.
func TestTrainPopulationsInParallel_EmptyPopulationsSlice(t *testing.T) {
	tmpDir := t.TempDir()
	result := trainPopulationsInParallel([]ParallelVCFTraining{}, tmpDir, 2)
	if result == nil {
		t.Error("expected non-nil map from trainPopulationsInParallel with empty input")
	}
	if len(result) != 0 {
		t.Errorf("expected empty map, got %d entries", len(result))
	}
}

// TestTrainPopulationsInParallel_NoVCFFiles exercises the path where population
// directories don't exist (no VCF files), so trainPopulationOnVCF returns an
// error and each entry is stored with Trained=false.
func TestTrainPopulationsInParallel_NoVCFFiles(t *testing.T) {
	tmpDir := t.TempDir()
	populations := []ParallelVCFTraining{
		{Population: "YRI", Region: "Africa", Individuals: 108},
		{Population: "CEU", Region: "Europe", Individuals: 99},
	}
	result := trainPopulationsInParallel(populations, tmpDir, 2)
	if len(result) != 2 {
		t.Fatalf("expected 2 result entries, got %d", len(result))
	}
	for pop, training := range result {
		if training.Trained {
			t.Errorf("expected population %s to be untrained (no VCF files in tmpDir)", pop)
		}
		if training.Error == nil {
			t.Errorf("expected error for population %s with no VCF files", pop)
		}
	}
}

// TestTrainPopulationsInParallel_SingleWorker ensures parallelism with maxParallel=1
// still processes all populations.
func TestTrainPopulationsInParallel_SingleWorker(t *testing.T) {
	tmpDir := t.TempDir()
	populations := []ParallelVCFTraining{
		{Population: "GBR", Region: "Europe", Individuals: 91},
	}
	result := trainPopulationsInParallel(populations, tmpDir, 1)
	if _, ok := result["GBR"]; !ok {
		t.Error("expected GBR entry in result map")
	}
}

// ---------- trainPopulationOnVCF ----------

// TestTrainPopulationOnVCF_EmptyDir verifies that when the vcfDataPath contains
// no VCF files for the given population, an error is returned.
func TestTrainPopulationOnVCF_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	training := &ParallelVCFTraining{
		Population:  "YRI",
		Region:      "Africa",
		Individuals: 108,
	}
	err := trainPopulationOnVCF(training, tmpDir)
	if err == nil {
		t.Error("expected error when no VCF files are present")
	}
}

// TestTrainPopulationOnVCF_SetsModel verifies that the function always sets a
// model on the ParallelVCFTraining struct, even when it returns an error due to
// missing files (model is initialized before the file search).
func TestTrainPopulationOnVCF_SetsModel(t *testing.T) {
	tmpDir := t.TempDir()
	training := &ParallelVCFTraining{
		Population:  "LWK",
		Region:      "Africa",
		Individuals: 99,
	}
	// Call will fail due to missing files, but model should still be set.
	_ = trainPopulationOnVCF(training, tmpDir)
	if training.Model == nil {
		t.Error("expected Model to be initialized even on failure")
	}
}

// TestTrainPopulationOnVCF_WithValidVCFFile verifies that when a real VCF file
// is present, training succeeds and populates the training fields.
func TestTrainPopulationOnVCF_WithValidVCFFile(t *testing.T) {
	vcfPath := createTestGzippedVCF(t)
	// The VCF file is in a tmpDir; we need to place it in a sub-directory
	// matching the population name (YRI) OR use the YRI_ prefix pattern.
	// createTestGzippedVCF puts it in t.TempDir()/chr1.vcf.gz but we need
	// to place it where findPopulationVCFFiles can find it.

	// Build a proper directory layout: basePath/YRI/chr1.vcf.gz
	basePath := t.TempDir()
	popDir := filepath.Join(basePath, "YRI")
	if err := os.MkdirAll(popDir, 0755); err != nil {
		t.Fatalf("failed to create pop dir: %v", err)
	}
	destPath := filepath.Join(popDir, "chr1.vcf.gz")
	data, err := os.ReadFile(vcfPath)
	if err != nil {
		t.Fatalf("failed to read test VCF: %v", err)
	}
	if err := os.WriteFile(destPath, data, 0600); err != nil {
		t.Fatalf("failed to write VCF to pop dir: %v", err)
	}

	training := &ParallelVCFTraining{
		Population:  "YRI",
		Region:      "Africa",
		Individuals: 108,
	}
	err = trainPopulationOnVCF(training, basePath)
	if err != nil {
		t.Fatalf("trainPopulationOnVCF returned unexpected error: %v", err)
	}
	if !training.Trained {
		// trainPopulationOnVCF does not set training.Trained; the caller does.
		// So just check that no error was returned and model is set.
	}
	if training.Model == nil {
		t.Error("expected Model to be set after successful training")
	}
	if training.ProcessedVariants == 0 {
		t.Error("expected ProcessedVariants > 0 after processing VCF file")
	}
}

// ---------- generateVCFTrainedSequences ----------

// TestGenerateVCFTrainedSequences_EmptyMap verifies that calling with an empty
// trainedModels map does not panic and produces no output.
func TestGenerateVCFTrainedSequences_EmptyMap(t *testing.T) {
	// Should not panic.
	generateVCFTrainedSequences(map[string]*ParallelVCFTraining{})
}

// TestGenerateVCFTrainedSequences_UntrainedModels verifies that entries with
// Trained=false are skipped gracefully.
func TestGenerateVCFTrainedSequences_UntrainedModels(t *testing.T) {
	models := map[string]*ParallelVCFTraining{
		"YRI": {Population: "YRI", Region: "Africa", Trained: false},
	}
	// Should not panic even though model is nil for untrained entry.
	generateVCFTrainedSequences(models)
}

// TestGenerateVCFTrainedSequences_TrainedModel verifies that a properly trained
// model results in sequence files being generated. We train via TrainFromSequence
// to avoid requiring real VCF files.
func TestGenerateVCFTrainedSequences_TrainedModel(t *testing.T) {
	// We need the output directory to be writable. generateVCFTrainedSequences
	// hardcodes "sequences/parallel_vcf_trained" relative to CWD, so we change
	// the working directory to a temp dir for this test.
	origWD, err := os.Getwd()
	if err != nil {
		t.Fatalf("failed to get working dir: %v", err)
	}
	tmpDir := t.TempDir()
	if err := os.Chdir(tmpDir); err != nil {
		t.Fatalf("failed to chdir to temp dir: %v", err)
	}
	defer func() {
		if err := os.Chdir(origWD); err != nil {
			t.Logf("warning: failed to restore working dir: %v", err)
		}
	}()

	// Build a trained ParallelVCFTraining entry using TrainFromSequence.
	// The Model field type is *authenticity.AdaptiveHMMModel (same package
	// already imported by main.go in this package).
	m := newTestAdaptiveModel(t, 4)
	seq := []byte("ACGTACGTGCATGCATACGTACGTGCATGCATATGCATGCACGTACGTACGT")
	if err := m.TrainFromSequence(seq); err != nil {
		t.Fatalf("TrainFromSequence returned error: %v", err)
	}

	training := &ParallelVCFTraining{
		Population:        "CEU",
		Region:            "Europe",
		Individuals:       99,
		Trained:           true,
		Model:             m,
		ProcessedVariants: 10,
		TrainingSequences: [][]byte{seq},
		ActualGCContent:   0.45,
		ActualComplexity:  0.70,
	}

	models := map[string]*ParallelVCFTraining{
		"CEU": training,
	}

	// Should not panic; may produce files in tmpDir/sequences/parallel_vcf_trained/CEU/
	generateVCFTrainedSequences(models)
}

// ============================================================
// Additional pure-helper function tests
// ============================================================

// --- isValidDNASequence edge cases ---

func TestIsValidDNASequence_PVT_OnlySpaces(t *testing.T) {
	if isValidDNASequence("   ") {
		t.Error("string of spaces should not be valid DNA")
	}
}

func TestIsValidDNASequence_PVT_OnlyNumbers(t *testing.T) {
	if isValidDNASequence("1234") {
		t.Error("string of digits should not be valid DNA")
	}
}

func TestIsValidDNASequence_PVT_MixedValidInvalid(t *testing.T) {
	if isValidDNASequence("ACGTx") {
		t.Error("sequence with trailing invalid char should fail")
	}
}

func TestIsValidDNASequence_PVT_LowercaseOnly(t *testing.T) {
	if !isValidDNASequence("atcg") {
		t.Error("all-lowercase atcg should be valid DNA (accepted via ToUpper)")
	}
}

func TestIsValidDNASequence_PVT_Uracil(t *testing.T) {
	if isValidDNASequence("ACGU") {
		t.Error("U (uracil) should not be valid DNA base")
	}
}

// --- calculateAverageGCContent additional coverage ---

func TestCalculateAverageGCContent_PVT_NilInput(t *testing.T) {
	if got := calculateAverageGCContent(nil); got != 0.0 {
		t.Errorf("calculateAverageGCContent(nil) = %f, want 0.0", got)
	}
}

func TestCalculateAverageGCContent_PVT_KnownFraction(t *testing.T) {
	// "ATCATCATCG": A=3,T=3,C=3,G=1 -> GC=4, total=10 -> 0.4
	seq := []byte("ATCATCATCG")
	got := calculateAverageGCContent([][]byte{seq})
	const want = 0.4
	if got < want-0.0001 || got > want+0.0001 {
		t.Errorf("calculateAverageGCContent = %f, want %f", got, want)
	}
}

func TestCalculateAverageGCContent_PVT_ThreeSequences(t *testing.T) {
	// combined: (2+4+2) GC / (4+8+4) total = 8/16 = 0.5
	seqs := [][]byte{
		[]byte("GCTA"),     // 2 GC / 4
		[]byte("GCGCATAT"), // 4 GC / 8
		[]byte("GCAA"),     // 2 GC / 4
	}
	got := calculateAverageGCContent(seqs)
	const want = 0.5
	if got < want-0.0001 || got > want+0.0001 {
		t.Errorf("calculateAverageGCContent(3 seqs) = %f, want %f", got, want)
	}
}

// --- calculateSequenceComplexity exact-value cases ---

func TestCalculateSequenceComplexity_PVT_KEqualsLength(t *testing.T) {
	// k == len(seq): 1 unique kmer, totalKmers=1 -> ratio = 1.0
	seq := []byte("ACGT")
	got := calculateSequenceComplexity(seq, 4)
	if got != 1.0 {
		t.Errorf("calculateSequenceComplexity(k==len) = %f, want 1.0", got)
	}
}

func TestCalculateSequenceComplexity_PVT_K1AllSameBase(t *testing.T) {
	// k=1, 10 A's: 1 unique kmer / 10 total = 0.1
	seq := []byte("AAAAAAAAAA")
	got := calculateSequenceComplexity(seq, 1)
	const want = 0.1
	if got < want-0.0001 || got > want+0.0001 {
		t.Errorf("calculateSequenceComplexity(k=1, all A) = %f, want %f", got, want)
	}
}

func TestCalculateSequenceComplexity_PVT_K1AllFourBases(t *testing.T) {
	// k=1, ACGT: 4 unique / 4 total = 1.0
	seq := []byte("ACGT")
	got := calculateSequenceComplexity(seq, 1)
	if got != 1.0 {
		t.Errorf("calculateSequenceComplexity(k=1, ACGT) = %f, want 1.0", got)
	}
}

func TestCalculateSequenceComplexity_PVT_KLargerThanSeq(t *testing.T) {
	got := calculateSequenceComplexity([]byte("AC"), 10)
	if got != 0.0 {
		t.Errorf("calculateSequenceComplexity(k>len) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_PVT_RepeatedMotif(t *testing.T) {
	// "ACGTACGT" k=4: kmers ACGT,CGTA,GTAC,TACG,ACGT
	// unique={ACGT,CGTA,GTAC,TACG}=4, total=5 -> 0.8
	seq := []byte("ACGTACGT")
	got := calculateSequenceComplexity(seq, 4)
	const want = 0.8
	if got < want-0.0001 || got > want+0.0001 {
		t.Errorf("calculateSequenceComplexity(ACGTACGT,k=4) = %f, want %f", got, want)
	}
}

// --- calculateAverageComplexity additional coverage ---

func TestCalculateAverageComplexity_PVT_NilInput(t *testing.T) {
	if got := calculateAverageComplexity(nil, 4); got != 0.0 {
		t.Errorf("calculateAverageComplexity(nil) = %f, want 0.0", got)
	}
}

func TestCalculateAverageComplexity_PVT_IdenticalSequencesAverage(t *testing.T) {
	// Two identical sequences: average == single complexity
	seq := []byte("ACGTACGT")
	single := calculateSequenceComplexity(seq, 4)
	got := calculateAverageComplexity([][]byte{seq, seq}, 4)
	if got < single-0.0001 || got > single+0.0001 {
		t.Errorf("average of identical seqs = %f, want %f", got, single)
	}
}

func TestCalculateAverageComplexity_PVT_OnlyShortExcluded(t *testing.T) {
	// Short seq excluded; only long one counts -> result equals single long complexity
	short := []byte("AC")
	long := []byte("ACGTACGT")
	wantSingle := calculateSequenceComplexity(long, 4)
	got := calculateAverageComplexity([][]byte{short, long}, 4)
	if got < wantSingle-0.0001 || got > wantSingle+0.0001 {
		t.Errorf("calculateAverageComplexity(short+long) = %f, want %f", got, wantSingle)
	}
}

func TestCalculateAverageComplexity_PVT_KEqualsSeqLen(t *testing.T) {
	// When k == len of every sequence, each has complexity 1.0 -> average = 1.0
	s1 := []byte("ACGT")
	s2 := []byte("GCTA")
	got := calculateAverageComplexity([][]byte{s1, s2}, 4)
	if got != 1.0 {
		t.Errorf("calculateAverageComplexity(k==len each) = %f, want 1.0", got)
	}
}

// --- detectAppleSiliconCores property tests ---

func TestDetectAppleSiliconCores_PVT_AlwaysPositive(t *testing.T) {
	if got := detectAppleSiliconCores(); got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestDetectAppleSiliconCores_PVT_ReasonableUpperBound(t *testing.T) {
	if got := detectAppleSiliconCores(); got > 256 {
		t.Errorf("detectAppleSiliconCores() = %d, unexpectedly large", got)
	}
}

func TestDetectAppleSiliconCores_PVT_Deterministic(t *testing.T) {
	first := detectAppleSiliconCores()
	second := detectAppleSiliconCores()
	if first != second {
		t.Errorf("detectAppleSiliconCores() not deterministic: %d != %d", first, second)
	}
}
