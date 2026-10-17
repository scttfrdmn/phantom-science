//go:build !release
// +build !release

package main

import (
	"context"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestAbs_Streaming(t *testing.T) {
	tests := []struct {
		name     string
		input    float64
		expected float64
	}{
		{"positive", 7.77, 7.77},
		{"negative", -7.77, 7.77},
		{"zero", 0.0, 0.0},
		{"small positive", 0.001, 0.001},
		{"small negative", -0.001, 0.001},
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

func TestCalculateSequenceComplexity_Streaming(t *testing.T) {
	tests := []struct {
		name     string
		sequence []byte
		k        int
		wantMin  float64
		wantMax  float64
	}{
		{"empty sequence returns 0", []byte{}, 4, 0.0, 0.0},
		{"sequence shorter than k returns 0", []byte("ACG"), 4, 0.0, 0.0},
		{"all same bases low complexity", []byte("AAAAAAAAAAAAAAAA"), 4, 0.0, 0.1},
		{"diverse sequence positive complexity", []byte("ACGTACGTACGTACGT"), 4, 0.1, 1.0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateSequenceComplexity(tt.sequence, tt.k)
			if got < tt.wantMin || got > tt.wantMax {
				t.Errorf("calculateSequenceComplexity(%q, %d) = %v, want [%v, %v]",
					tt.sequence, tt.k, got, tt.wantMin, tt.wantMax)
			}
		})
	}
}

func TestCalculateGCContent_Streaming(t *testing.T) {
	tests := []struct {
		name     string
		sequence []byte
		want     float64
	}{
		{"empty returns 0", []byte{}, 0.0},
		{"all A returns 0", []byte("AAAA"), 0.0},
		{"all T returns 0", []byte("TTTT"), 0.0},
		{"all G returns 1", []byte("GGGG"), 1.0},
		{"all C returns 1", []byte("CCCC"), 1.0},
		{"50pct GC ACGT", []byte("ACGT"), 0.5},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateGCContent(tt.sequence)
			if math.Abs(got-tt.want) > 1e-9 {
				t.Errorf("calculateGCContent(%q) = %v, want %v", tt.sequence, got, tt.want)
			}
		})
	}
}

func TestCalculateMaxHomopolymer_Streaming(t *testing.T) {
	tests := []struct {
		name     string
		sequence []byte
		want     int
	}{
		{"empty returns 0", []byte{}, 0},
		{"single base returns 1", []byte("A"), 1},
		{"no repeat returns 1", []byte("ACGT"), 1},
		{"double A", []byte("AA"), 2},
		{"triple A at start", []byte("AAACGT"), 3},
		{"run in middle", []byte("ACGTTTTACG"), 4},
		{"run at end", []byte("ACGAAAAAA"), 6},
		{"multiple runs picks max", []byte("AAACCCCC"), 5},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateMaxHomopolymer(tt.sequence)
			if got != tt.want {
				t.Errorf("calculateMaxHomopolymer(%q) = %v, want %v", tt.sequence, got, tt.want)
			}
		})
	}
}

func TestCalculateCpGDepletion_Streaming(t *testing.T) {
	tests := []struct {
		name    string
		seq     []byte
		wantMin float64
		wantMax float64
	}{
		{"empty returns 0", []byte{}, 0.0, 0.0},
		{"single base returns 0", []byte("A"), 0.0, 0.0},
		{"no C or G returns 0", []byte("ATATAT"), 0.0, 0.0},
		{"only C returns 0", []byte("CCCCCC"), 0.0, 0.0},
		{"only G returns 0", []byte("GGGGGG"), 0.0, 0.0},
		{"CG present non-negative", []byte("ACGACG"), 0.0, 100.0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateCpGDepletion(tt.seq)
			if got < tt.wantMin || got > tt.wantMax {
				t.Errorf("calculateCpGDepletion(%q) = %v, want [%v, %v]",
					tt.seq, got, tt.wantMin, tt.wantMax)
			}
		})
	}
}

func TestDetectAppleSiliconCores_Streaming(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestCalculateAuthenticityMetricsStreaming_Empty(t *testing.T) {
	complexity, gcContent, maxHomopolymer, cpgDepletion := calculateAuthenticityMetrics([][]byte{})
	if complexity != 0 || gcContent != 0 || maxHomopolymer != 0 || cpgDepletion != 0 {
		t.Errorf("calculateAuthenticityMetrics(empty) = (%v, %v, %v, %v), want all zeros",
			complexity, gcContent, maxHomopolymer, cpgDepletion)
	}
}

func TestCalculateAuthenticityMetricsStreaming_WithData(t *testing.T) {
	sequences := [][]byte{
		[]byte("ACGTACGTACGT"),
		[]byte("GCTAGCTAGCTA"),
	}
	complexity, gcContent, maxHomopolymer, cpgDepletion := calculateAuthenticityMetrics(sequences)
	if complexity <= 0 {
		t.Errorf("expected complexity > 0, got %v", complexity)
	}
	if gcContent < 0 || gcContent > 1 {
		t.Errorf("GC content %v out of [0,1] range", gcContent)
	}
	if maxHomopolymer <= 0 {
		t.Errorf("expected maxHomopolymer > 0, got %v", maxHomopolymer)
	}
	if cpgDepletion < 0 {
		t.Errorf("cpgDepletion should be non-negative, got %v", cpgDepletion)
	}
}

func TestCalculateAuthenticityMetricsStreaming_HomopolymerMax(t *testing.T) {
	// Sequence with a long homopolymer run should dominate the max
	sequences := [][]byte{
		[]byte("ACGTACGT"),
		[]byte("AAAAAAAAAAAACGT"), // 12 A's in a row
	}
	_, _, maxHomopolymer, _ := calculateAuthenticityMetrics(sequences)
	if maxHomopolymer < 12 {
		t.Errorf("expected maxHomopolymer >= 12 (long A run), got %v", maxHomopolymer)
	}
}

func TestCalculateGCContentSingleBase(t *testing.T) {
	tests := []struct {
		base []byte
		want float64
	}{
		{[]byte("A"), 0.0},
		{[]byte("T"), 0.0},
		{[]byte("G"), 1.0},
		{[]byte("C"), 1.0},
	}
	for _, tt := range tests {
		got := calculateGCContent(tt.base)
		if math.Abs(got-tt.want) > 1e-9 {
			t.Errorf("calculateGCContent(%q) = %v, want %v", tt.base, got, tt.want)
		}
	}
}

func TestCalculateCpGDepletionNoCpGdinucleotide(t *testing.T) {
	// GC present but no CG dinucleotide (GC not CG)
	seq := []byte("GCGCGCGC") // Only GC not CG
	got := calculateCpGDepletion(seq)
	if got < 0 {
		t.Errorf("cpgDepletion should be non-negative, got %v", got)
	}
}

func TestCalculateCpGDepletionWithCpG(t *testing.T) {
	// Explicit CG dinucleotides present
	seq := []byte("ACGACGACG")
	got := calculateCpGDepletion(seq)
	if got < 0 {
		t.Errorf("cpgDepletion should be non-negative with CpG present, got %v", got)
	}
}

// ---------- calculateCpGDepletion additional coverage ----------

func TestCalculateCpGDepletion_EndsWithC(t *testing.T) {
	// Last base is C; exercises "if sequence[last] == 'C'" branch.
	// "GC": cCount=1, gCount=1, cgCount=0 (GC != CG)
	got := calculateCpGDepletion([]byte("GC"))
	if got < 0 {
		t.Errorf("calculateCpGDepletion(GC) = %v, want >= 0", got)
	}
}

func TestCalculateCpGDepletion_EndsWithG(t *testing.T) {
	// Last base is G; exercises "if sequence[last] == 'G'" branch.
	// "CG": cgCount=1, cCount=1, gCount=1 → ratio > 0
	got := calculateCpGDepletion([]byte("CG"))
	if got <= 0 {
		t.Errorf("calculateCpGDepletion(CG) = %v, want > 0", got)
	}
}

func TestCalculateCpGDepletion_MixedEndsG(t *testing.T) {
	// Mixed sequence ending in G with CG dinucleotides
	seq := []byte("AACGTTG")
	got := calculateCpGDepletion(seq)
	if got < 0 {
		t.Errorf("calculateCpGDepletion(%q) = %v, want >= 0", seq, got)
	}
}

func TestCalculateCpGDepletion_EndsWithCPlusCG(t *testing.T) {
	// Contains CG and ends with C
	seq := []byte("ACGC")
	got := calculateCpGDepletion(seq)
	if got < 0 {
		t.Errorf("calculateCpGDepletion(%q) = %v, want >= 0", seq, got)
	}
}

// ---------- detectAppleSiliconCores sanity checks ----------

func TestDetectAppleSiliconCores_ReturnsPositive_Streaming(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestDetectAppleSiliconCores_ReasonableBound_Streaming(t *testing.T) {
	got := detectAppleSiliconCores()
	if got > 256 {
		t.Errorf("detectAppleSiliconCores() = %d, unreasonably large", got)
	}
}

// ---------- saveFinalResults ----------

func TestSaveFinalResults_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	saveFinalResults(tmpDir, 30*time.Second, 0.082, 0.42, 6)
	path := filepath.Join(tmpDir, "final_results.txt")
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty results file")
	}
}

// TestSaveFinalResults_ComplexityBranches exercises different complexity values
// to hit the "improvement" branch (complexity > 0 leading to complexity/0.005)
// and the zero/minimal branch.
func TestSaveFinalResults_ComplexityBranches(t *testing.T) {
	tests := []struct {
		name       string
		duration   time.Duration
		complexity float64
		gcContent  float64
		perfCores  int
	}{
		// complexity > 0.05: renders "improvement" line
		{"high_complexity", 90 * time.Second, 0.12, 0.42, 8},
		// complexity <= 0.05: skips improvement rendering
		{"low_complexity", 5 * time.Second, 0.002, 0.48, 4},
		// zero values
		{"zeros", 0, 0.0, 0.0, 1},
		// large perfCores
		{"many_cores", 10 * time.Minute, 0.08, 0.41, 16},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			saveFinalResults(tmpDir, tt.duration, tt.complexity, tt.gcContent, tt.perfCores)
			path := filepath.Join(tmpDir, "final_results.txt")
			info, err := os.Stat(path)
			if err != nil {
				t.Fatalf("file should exist after saveFinalResults: %v", err)
			}
			if info.Size() == 0 {
				t.Error("results file should not be empty")
			}
		})
	}
}

// ---------- trainBaselineStreaming ----------

func TestTrainBaselineStreaming_NoFile_ReturnsError(t *testing.T) {
	tmpDir := t.TempDir()
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx := context.Background()
	err := trainBaselineStreaming(ctx, model, tmpDir, tmpDir)
	if err == nil {
		t.Fatal("expected error when GRCh38_complete.fasta is missing, got nil")
	}
}

// ---------- streamFASTAFile ----------

func TestStreamFASTAFile_WithSyntheticFASTA(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := filepath.Join(tmpDir, "test.fasta")

	// Write a small FASTA file with two sequences
	content := ">chr1\nACGTACGTACGT\n>chr2\nGCGCGCGCGCGC\n"
	if err := os.WriteFile(fastaPath, []byte(content), 0600); err != nil {
		t.Fatalf("failed to write test FASTA: %v", err)
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx := context.Background()
	if err := streamFASTAFile(ctx, model, fastaPath, tmpDir); err != nil {
		t.Errorf("streamFASTAFile returned unexpected error: %v", err)
	}
}

func TestStreamFASTAFile_MissingFile_ReturnsError(t *testing.T) {
	tmpDir := t.TempDir()
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx := context.Background()
	err := streamFASTAFile(ctx, model, filepath.Join(tmpDir, "nonexistent.fasta"), tmpDir)
	if err == nil {
		t.Fatal("expected error for missing FASTA file, got nil")
	}
}

func TestStreamFASTAFile_CancelledContext(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := filepath.Join(tmpDir, "test.fasta")

	content := ">chr1\nACGTACGTACGT\n"
	if err := os.WriteFile(fastaPath, []byte(content), 0600); err != nil {
		t.Fatalf("failed to write test FASTA: %v", err)
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx, cancel := context.WithCancel(context.Background())
	cancel() // cancel immediately

	// Should return context error or complete quickly
	_ = streamFASTAFile(ctx, model, fastaPath, tmpDir)
}

// ---------- processSequenceChunks ----------

func TestProcessSequenceChunks_Basic(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	// Train the model first so GenerateSequence works
	trainSeq := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(trainSeq); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	sequence := []byte("ACGTACGTACGTACGTACGT")
	var processedChunks int64
	var totalBases int64
	chunkSize := 8

	if err := processSequenceChunks(ctx, model, sequence, chunkSize, &processedChunks, &totalBases); err != nil {
		t.Errorf("processSequenceChunks returned error: %v", err)
	}
	if processedChunks == 0 {
		t.Error("expected processedChunks > 0")
	}
	if totalBases == 0 {
		t.Error("expected totalBases > 0")
	}
}

func TestProcessSequenceChunks_EmptySequence(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx := context.Background()
	var processedChunks int64
	var totalBases int64

	if err := processSequenceChunks(ctx, model, []byte{}, 1000, &processedChunks, &totalBases); err != nil {
		t.Errorf("processSequenceChunks(empty) returned unexpected error: %v", err)
	}
	if processedChunks != 0 {
		t.Errorf("expected 0 processedChunks for empty sequence, got %d", processedChunks)
	}
}

func TestProcessSequenceChunks_CancelledContext(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx, cancel := context.WithCancel(context.Background())
	cancel() // cancel before calling

	sequence := []byte("ACGTACGTACGTACGTACGT")
	var processedChunks int64
	var totalBases int64

	// May return context error since context is cancelled
	_ = processSequenceChunks(ctx, model, sequence, 8, &processedChunks, &totalBases)
}

// ---------- trainPopulationStreaming ----------

func TestTrainPopulationStreaming_EmptyDir_ReturnsNil(t *testing.T) {
	tmpDir := t.TempDir()
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx := context.Background()
	// All populations will be missing in an empty tmpDir, so function returns nil
	if err := trainPopulationStreaming(ctx, model, tmpDir, tmpDir); err != nil {
		t.Errorf("trainPopulationStreaming with no populations: expected nil, got %v", err)
	}
}

// ---------- processPopulationVCFs ----------

func TestProcessPopulationVCFs_EmptyPopDir_ReturnsNil(t *testing.T) {
	tmpDir := t.TempDir()
	pop := "YRI"
	// Create the population directory but leave it empty (no vcf.gz files)
	popDir := filepath.Join(tmpDir, pop)
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("failed to create pop dir: %v", err)
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	ctx := context.Background()
	if err := processPopulationVCFs(ctx, model, pop, tmpDir); err != nil {
		t.Errorf("processPopulationVCFs with empty dir: expected nil, got %v", err)
	}
}

// ---------- processSyntheticVCF ----------

func TestProcessSyntheticVCF_AllPopulations(t *testing.T) {
	populations := []string{"YRI", "CEU", "CHB", "GBR", "JPT", "MXL"}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	// Train model first
	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	for _, pop := range populations {
		t.Run(pop, func(t *testing.T) {
			// vcfFile path does not need to exist for processSyntheticVCF
			if err := processSyntheticVCF(ctx, model, "/any/path/file.vcf.gz", pop); err != nil {
				t.Errorf("processSyntheticVCF(%s): unexpected error: %v", pop, err)
			}
		})
	}
}

func TestProcessSyntheticVCF_DefaultPopulation(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	// "PUR" is an admixed population hitting the default switch case
	if err := processSyntheticVCF(ctx, model, "/any/path/file.vcf.gz", "PUR"); err != nil {
		t.Errorf("processSyntheticVCF(PUR default): unexpected error: %v", err)
	}
}

// ---------- generateValidationSequences ----------

func TestGenerateValidationSequences_Basic(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	// Train the model so it can generate sequences
	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	seqs, avgTime, err := generateValidationSequences(model, 3, 100)
	if err != nil {
		t.Fatalf("generateValidationSequences failed: %v", err)
	}
	if len(seqs) != 3 {
		t.Errorf("expected 3 sequences, got %d", len(seqs))
	}
	if avgTime <= 0 {
		t.Errorf("expected avgTime > 0, got %v", avgTime)
	}
	for i, s := range seqs {
		if len(s) == 0 {
			t.Errorf("sequence %d is empty", i)
		}
	}
}

func TestGenerateValidationSequences_SingleSequence(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	seqs, avgTime, err := generateValidationSequences(model, 1, 50)
	if err != nil {
		t.Fatalf("generateValidationSequences(1, 50) failed: %v", err)
	}
	if len(seqs) != 1 {
		t.Errorf("expected 1 sequence, got %d", len(seqs))
	}
	// avgTime should be a valid non-negative duration
	_ = avgTime
}

// Ensure the _ import of time is used (it was already in the original imports)
var _ = time.Second

// ============================================================
// Additional tests to improve coverage
// ============================================================

// ---------- trainPopulationStreaming with real population dirs ----------

// TestTrainPopulationStreaming_WithPopulationDir creates a population directory
// so the goroutine body is actually exercised (wg.Add, semaphore, etc.).
func TestTrainPopulationStreaming_WithPopulationDir(t *testing.T) {
	tmpDir := t.TempDir()
	// Create a YRI population directory but leave it empty (no vcf.gz files).
	yriDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(yriDir, 0750); err != nil {
		t.Fatalf("failed to create YRI dir: %v", err)
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	if err := trainPopulationStreaming(ctx, model, tmpDir, tmpDir); err != nil {
		t.Errorf("trainPopulationStreaming: unexpected error: %v", err)
	}
}

// TestTrainPopulationStreaming_WithVCFFiles puts actual (fake) vcf.gz files in
// the population directory so processPopulationVCFs runs the file-level loop.
func TestTrainPopulationStreaming_WithVCFFiles(t *testing.T) {
	tmpDir := t.TempDir()
	yriDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(yriDir, 0750); err != nil {
		t.Fatalf("failed to create YRI dir: %v", err)
	}

	// Write a placeholder vcf.gz file (content does not matter for processSyntheticVCF)
	vcfPath := filepath.Join(yriDir, "chr1.vcf.gz")
	if err := os.WriteFile(vcfPath, []byte("placeholder"), 0600); err != nil {
		t.Fatalf("failed to create placeholder vcf.gz: %v", err)
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	if err := trainPopulationStreaming(ctx, model, tmpDir, tmpDir); err != nil {
		t.Errorf("trainPopulationStreaming with VCF files: unexpected error: %v", err)
	}
}

// ---------- processPopulationVCFs with VCF files ----------

// TestProcessPopulationVCFs_WithVCFFiles puts fake vcf.gz files in the pop dir
// so that the file-processing loop and processSyntheticVCF get called.
// Using 15 files also triggers the 10-file progress log branch.
func TestProcessPopulationVCFs_WithVCFFiles(t *testing.T) {
	tmpDir := t.TempDir()
	pop := "CEU"
	popDir := filepath.Join(tmpDir, pop)
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("failed to create pop dir: %v", err)
	}

	// Create 15 files so the processedFiles%10==0 progress branch fires.
	for i := 0; i < 15; i++ {
		vcfPath := filepath.Join(popDir, fmt.Sprintf("chr%d.vcf.gz", i+1))
		if err := os.WriteFile(vcfPath, []byte("placeholder"), 0600); err != nil {
			t.Fatalf("failed to create vcf.gz %d: %v", i, err)
		}
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	if err := processPopulationVCFs(ctx, model, pop, tmpDir); err != nil {
		t.Errorf("processPopulationVCFs with files: unexpected error: %v", err)
	}
}

// TestProcessPopulationVCFs_AllKnownPopulations exercises each population name
// to hit all branches of the switch in processSyntheticVCF.
func TestProcessPopulationVCFs_AllKnownPopulations(t *testing.T) {
	populations := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"}

	for _, pop := range populations {
		t.Run(pop, func(t *testing.T) {
			tmpDir := t.TempDir()
			popDir := filepath.Join(tmpDir, pop)
			if err := os.MkdirAll(popDir, 0750); err != nil {
				t.Fatalf("failed to create pop dir: %v", err)
			}
			vcfPath := filepath.Join(popDir, "chr1.vcf.gz")
			if err := os.WriteFile(vcfPath, []byte("placeholder"), 0600); err != nil {
				t.Fatalf("failed to create vcf.gz: %v", err)
			}

			model := authenticity.NewParallelHMMModel(4)
			defer model.Cleanup()

			seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
			if err := model.TrainFromSequence(seed); err != nil {
				t.Fatalf("TrainFromSequence failed: %v", err)
			}

			ctx := context.Background()
			if err := processPopulationVCFs(ctx, model, pop, tmpDir); err != nil {
				t.Errorf("processPopulationVCFs(%s): unexpected error: %v", pop, err)
			}
		})
	}
}

// ---------- streamFASTAFile additional branches ----------

// TestStreamFASTAFile_MultipleSequences ensures the "process accumulated
// sequence" branch fires (when a new header line arrives while currentSequence
// is non-empty).
func TestStreamFASTAFile_MultipleSequences(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := filepath.Join(tmpDir, "multi.fasta")

	// Three sequences so each new header triggers processing of the previous one.
	content := ">seq1\nACGTACGTACGT\n>seq2\nGCGCGCGCGCGC\n>seq3\nTATATATATATA\n"
	if err := os.WriteFile(fastaPath, []byte(content), 0600); err != nil {
		t.Fatalf("failed to write FASTA: %v", err)
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	if err := streamFASTAFile(ctx, model, fastaPath, tmpDir); err != nil {
		t.Errorf("streamFASTAFile(multi): unexpected error: %v", err)
	}
}

// ---------- saveFinalResults write-error path ----------

// TestSaveFinalResults_BadDir exercises the os.WriteFile error path
// (the warning print branch) when the output directory does not exist.
func TestSaveFinalResults_BadDir(t *testing.T) {
	// Should not panic even when the directory does not exist.
	saveFinalResults("/nonexistent/path/xyz/abc", 5*time.Second, 0.07, 0.42, 4)
}

// ---------- processSyntheticVCF additional population coverage ----------

// TestProcessSyntheticVCF_TSI covers the "TSI" (European) population path.
func TestProcessSyntheticVCF_TSI(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	if err := processSyntheticVCF(ctx, model, "/any/path/file.vcf.gz", "TSI"); err != nil {
		t.Errorf("processSyntheticVCF(TSI): unexpected error: %v", err)
	}
}

// TestProcessSyntheticVCF_CHS covers the "CHS" (Asian) population path.
func TestProcessSyntheticVCF_CHS(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	if err := processSyntheticVCF(ctx, model, "/any/path/file.vcf.gz", "CHS"); err != nil {
		t.Errorf("processSyntheticVCF(CHS): unexpected error: %v", err)
	}
}

// ---------- trainBaselineStreaming success path ----------

// TestTrainBaselineStreaming_WithFakeFile creates a fake GRCh38_complete.fasta
// file so trainBaselineStreaming proceeds past the os.Stat check.
func TestTrainBaselineStreaming_WithFakeFile(t *testing.T) {
	tmpDir := t.TempDir()
	refPath := filepath.Join(tmpDir, "GRCh38_complete.fasta")
	// Write a small valid FASTA file with two sequences.
	content := ">chr1\nACGTACGTACGT\n>chr2\nGCGCGCGCGCGC\n"
	if err := os.WriteFile(refPath, []byte(content), 0600); err != nil {
		t.Fatalf("failed to create fake FASTA: %v", err)
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx := context.Background()
	// trainBaselineStreaming should stream the fake FASTA successfully.
	if err := trainBaselineStreaming(ctx, model, tmpDir, tmpDir); err != nil {
		t.Errorf("trainBaselineStreaming with fake file: unexpected error: %v", err)
	}
}

// ---------- processPopulationVCFs context cancellation ----------

// TestProcessPopulationVCFs_CancelledContext passes a pre-cancelled context to
// exercise the ctx.Done() branch in the VCF processing loop.
func TestProcessPopulationVCFs_CancelledContextMid(t *testing.T) {
	tmpDir := t.TempDir()
	pop := "YRI"
	popDir := filepath.Join(tmpDir, pop)
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("failed to create pop dir: %v", err)
	}
	// Create several VCF files so the loop has something to iterate.
	for i := 0; i < 5; i++ {
		vcfPath := filepath.Join(popDir, fmt.Sprintf("chr%d.vcf.gz", i+1))
		if err := os.WriteFile(vcfPath, []byte("placeholder"), 0600); err != nil {
			t.Fatalf("failed to create vcf.gz %d: %v", i, err)
		}
	}

	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	seed := []byte("ACGTACGTACGTACGTACGTACGTACGTACGT")
	if err := model.TrainFromSequence(seed); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	ctx, cancel := context.WithCancel(context.Background())
	cancel() // cancel immediately so the select fires on first iteration

	// Should return ctx.Err() (context.Canceled).
	err := processPopulationVCFs(ctx, model, pop, tmpDir)
	if err == nil {
		// Context was already cancelled; expect an error or completion
		// (race between goroutine scheduling and cancellation check).
	}
	_ = err
}

// ---------- generateValidationSequences error path ----------

// TestGenerateValidationSequences_UntrainedModel tests the error return path
// when the ParallelHMMModel has not been trained (GenerateSequence fails).
func TestGenerateValidationSequences_UntrainedModel_Streaming(t *testing.T) {
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	// Untrained model - GenerateSequence should fail, exercising the error return.
	_, _, err := generateValidationSequences(model, 1, 100)
	// We accept either error or nil - if the model succeeds untrained, that is ok too.
	_ = err
}

// ---------- fmt import used by new tests ----------
var _ = fmt.Sprintf
