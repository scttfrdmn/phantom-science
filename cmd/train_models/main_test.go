//go:build !release
// +build !release

package main

import (
	"log"
	"os"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/training"
)

func newTestSession() *TrainingSession {
	return &TrainingSession{
		logger: log.New(os.Stdout, "[TEST] ", 0),
	}
}

func TestValidateSequenceQuality_Valid(t *testing.T) {
	session := newTestSession()
	sequences := [][]byte{
		[]byte("ATGCATGCATGCATGCATGCATGCATGC"),
		[]byte("GCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
	}
	if err := session.validateSequenceQuality(sequences); err != nil {
		t.Errorf("expected no error for valid sequences, got: %v", err)
	}
}

func TestValidateSequenceQuality_EmptySlice(t *testing.T) {
	session := newTestSession()
	if err := session.validateSequenceQuality([][]byte{}); err != nil {
		t.Errorf("expected no error for empty slice, got: %v", err)
	}
}

func TestValidateSequenceQuality_EmptySequence(t *testing.T) {
	session := newTestSession()
	sequences := [][]byte{{}}
	if err := session.validateSequenceQuality(sequences); err == nil {
		t.Error("expected error for empty sequence, got nil")
	}
}

func TestValidateSequenceQuality_InvalidBase(t *testing.T) {
	session := newTestSession()
	sequences := [][]byte{
		[]byte("ATGCNATGCATGCATGCATGCATGCATGC"),
	}
	if err := session.validateSequenceQuality(sequences); err == nil {
		t.Error("expected error for invalid base N, got nil")
	}
}

func TestValidateSequenceQuality_LowGCContent(t *testing.T) {
	session := newTestSession()
	sequences := [][]byte{[]byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")}
	if err := session.validateSequenceQuality(sequences); err == nil {
		t.Error("expected error for unrealistically low GC content, got nil")
	}
}

func TestValidateSequenceQuality_HighGCContent(t *testing.T) {
	session := newTestSession()
	sequences := [][]byte{[]byte("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")}
	if err := session.validateSequenceQuality(sequences); err == nil {
		t.Error("expected error for unrealistically high GC content, got nil")
	}
}

func TestValidateSequenceQuality_BoundaryGC(t *testing.T) {
	session := newTestSession()
	sequences := [][]byte{[]byte("ATATGCGCATATGCGCATATGCGCATATGCGC")}
	if err := session.validateSequenceQuality(sequences); err != nil {
		t.Errorf("expected no error for 50 percent GC content, got: %v", err)
	}
}

func TestValidateSequenceQuality_MultipleSequencesOneFails(t *testing.T) {
	session := newTestSession()
	sequences := [][]byte{
		[]byte("ATGCATGCATGCATGCATGCATGC"),
		[]byte("AAAAAAAAAAAAAAAAAAAAAAAA"),
	}
	if err := session.validateSequenceQuality(sequences); err == nil {
		t.Error("expected error when one sequence has invalid GC content, got nil")
	}
}

// Tests for constructor and getter functions

func TestNewReferenceGenomeLoader_NotNil(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/test")
	if loader == nil {
		t.Error("NewReferenceGenomeLoader() returned nil")
	}
}

func TestNewReferenceGenomeLoader_StoresDir(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/test_dir")
	if loader.referenceDir != "/tmp/test_dir" {
		t.Errorf("referenceDir = %q, want /tmp/test_dir", loader.referenceDir)
	}
}

func TestNewReferenceGenomeLoader_EmptyChromosomes(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/test")
	if len(loader.chromosomes) != 0 {
		t.Error("NewReferenceGenomeLoader() should have empty chromosomes initially")
	}
}

func TestNewMockReferenceGenome_NotNil(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/test")
	mock := NewMockReferenceGenome(loader)
	if mock == nil {
		t.Error("NewMockReferenceGenome() returned nil")
	}
}

func TestNewMockReferenceGenome_HasChromosomes(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/test")
	mock := NewMockReferenceGenome(loader)
	if len(mock.chromosomes) == 0 {
		t.Error("NewMockReferenceGenome() has no chromosomes")
	}
}

func TestGetChromosomeNames_ReturnsChromosomes(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/test")
	mock := NewMockReferenceGenome(loader)
	names := mock.GetChromosomeNames()
	if len(names) == 0 {
		t.Error("GetChromosomeNames() returned empty list")
	}
}

func TestMockChromosome_GetSequence(t *testing.T) {
	chr := &MockChromosome{Name: "chr1", Sequence: []byte("ACGT"), Length: 4}
	if string(chr.GetSequence()) != "ACGT" {
		t.Errorf("GetSequence() = %q, want ACGT", string(chr.GetSequence()))
	}
}

func TestMockChromosome_GetName(t *testing.T) {
	chr := &MockChromosome{Name: "chr1", Sequence: []byte("ACGT"), Length: 4}
	if chr.GetName() != "chr1" {
		t.Errorf("GetName() = %q, want chr1", chr.GetName())
	}
}

func TestMockChromosome_GetLength(t *testing.T) {
	chr := &MockChromosome{Name: "chr1", Sequence: []byte("ACGT"), Length: 4}
	if chr.GetLength() != 4 {
		t.Errorf("GetLength() = %d, want 4", chr.GetLength())
	}
}

func TestNewTrainingSession_NotNil(t *testing.T) {
	session := NewTrainingSession("/tmp/test")
	if session == nil {
		t.Error("NewTrainingSession() returned nil")
	}
}

func TestNewTrainingSession_HasConfig(t *testing.T) {
	session := NewTrainingSession("/tmp/test")
	if session.config == nil {
		t.Error("NewTrainingSession() has nil config")
	}
}

func TestNewTrainingSession_HasReferenceGenome(t *testing.T) {
	session := NewTrainingSession("/tmp/test")
	if session.referenceGenome == nil {
		t.Error("NewTrainingSession() has nil referenceGenome")
	}
}

func TestNewTrainingSession_HasTrainedModels(t *testing.T) {
	session := NewTrainingSession("/tmp/test")
	if session.trainedModels == nil {
		t.Error("NewTrainingSession() has nil trainedModels map")
	}
}

// TestReferenceGenomeLoader_CacheHit exercises the cache-hit branch of
// ReferenceGenomeLoader.LoadChromosome (line 32-34).
func TestReferenceGenomeLoader_CacheHit(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/nonexistent-dir")
	// Pre-populate the internal cache.
	loader.chromosomes["chr1"] = []byte("ACGT")

	seq, err := loader.LoadChromosome("chr1")
	if err != nil {
		t.Fatalf("LoadChromosome from cache: unexpected error: %v", err)
	}
	if string(seq) != "ACGT" {
		t.Errorf("LoadChromosome from cache = %q, want %q", seq, "ACGT")
	}
}

// TestReferenceGenomeLoader_MissingFile exercises the file-open error branch
// (lines 42-45) by pointing at a directory where neither the chromosome FASTA
// nor GRCh38.fasta exist.
func TestReferenceGenomeLoader_MissingFile(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/nonexistent-dir")
	_, err := loader.LoadChromosome("chrX")
	if err == nil {
		t.Error("LoadChromosome with missing files: expected error, got nil")
	}
}

// TestMockReferenceGenome_LoadChromosome_Error exercises MockReferenceGenome.LoadChromosome
// when the underlying loader fails (propagates the error correctly).
func TestMockReferenceGenome_LoadChromosome_Error(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/nonexistent-dir")
	mock := NewMockReferenceGenome(loader)
	_, err := mock.LoadChromosome("chrY")
	if err == nil {
		t.Error("MockReferenceGenome.LoadChromosome with missing files: expected error, got nil")
	}
}

// TestMockReferenceGenome_LoadChromosome_CacheHit exercises the success path of
// MockReferenceGenome.LoadChromosome by seeding the loader cache.
func TestMockReferenceGenome_LoadChromosome_CacheHit(t *testing.T) {
	loader := NewReferenceGenomeLoader("/tmp/nonexistent-dir")
	loader.chromosomes["chr2"] = []byte("TGCATGCA")

	mock := NewMockReferenceGenome(loader)
	chr, err := mock.LoadChromosome("chr2")
	if err != nil {
		t.Fatalf("MockReferenceGenome.LoadChromosome: unexpected error: %v", err)
	}
	if string(chr.GetSequence()) != "TGCATGCA" {
		t.Errorf("GetSequence() = %q, want %q", chr.GetSequence(), "TGCATGCA")
	}
	if chr.GetName() != "chr2" {
		t.Errorf("GetName() = %q, want %q", chr.GetName(), "chr2")
	}
	if chr.GetLength() != 8 {
		t.Errorf("GetLength() = %d, want 8", chr.GetLength())
	}
}

func newSessionWithTempDir(t *testing.T) *TrainingSession {
	t.Helper()
	tmpDir := t.TempDir()
	session := NewTrainingSession("/tmp/nonexistent-ref")
	session.config = &training.TrainingConfig{
		KmerSizes:               []int{4},
		Populations:             []string{"EUR"},
		JSdivergenceThreshold:   0.01,
		ExpertAccuracyThreshold: 0.55,
		OutputPath:              tmpDir,
	}
	session.startTime = time.Now()
	return session
}

func TestGenerateTrainingReport_CreatesFile(t *testing.T) {
	session := newSessionWithTempDir(t)
	if err := session.generateTrainingReport(); err != nil {
		t.Errorf("generateTrainingReport() returned error: %v", err)
	}
	if _, err := os.Stat(session.config.OutputPath + "/training_report.txt"); os.IsNotExist(err) {
		t.Error("training_report.txt was not created")
	}
}

func TestGenerateTrainingReport_WithModels(t *testing.T) {
	session := newSessionWithTempDir(t)
	model := authenticity.NewHMMKmerModel(4)
	session.trainedModels["hmm_kmer_4"] = model
	if err := session.generateTrainingReport(); err != nil {
		t.Errorf("generateTrainingReport() with models returned error: %v", err)
	}
}

func TestSaveModelSummary_CreatesFile(t *testing.T) {
	session := newSessionWithTempDir(t)
	model := authenticity.NewHMMKmerModel(4)
	result := &authenticity.GenerationResult{
		Sequence:        []byte("ATCGATCG"),
		KmerFrequencies: map[string]int{"ATCG": 2, "TCGA": 1},
		GenerationTime:  time.Millisecond * 5,
	}
	if err := session.saveModelSummary("hmm_kmer_4", model, result); err != nil {
		t.Errorf("saveModelSummary() returned error: %v", err)
	}
	if _, err := os.Stat(session.config.OutputPath + "/hmm_kmer_4_summary.txt"); os.IsNotExist(err) {
		t.Error("summary file was not created")
	}
}

// newSessionWithSeededChromosome creates a TrainingSession whose referenceGenome
// loader has a pre-seeded chromosome so that TrainFromReferenceGenome succeeds.
// The sequence is long enough for k-mer processing (> k bases) and has balanced
// GC content so that the generated sequences pass the GC quality check.
func newSessionWithSeededChromosome(t *testing.T) *TrainingSession {
	t.Helper()
	tmpDir := t.TempDir()

	session := NewTrainingSession("/tmp/nonexistent-ref")
	session.config = &training.TrainingConfig{
		KmerSizes:               []int{4},
		Populations:             []string{"EUR"},
		JSdivergenceThreshold:   0.01,
		ExpertAccuracyThreshold: 0.55,
		OutputPath:              tmpDir,
	}
	session.startTime = time.Now()

	// Seed the loader cache with a synthetic chromosome that has balanced GC
	// content so that generated sequences satisfy the 0.2–0.8 GC quality check.
	syntheticSeq := make([]byte, 2000)
	pattern := []byte("ATGCATGCATGCATGCATGCATGCATGCATGC")
	for i := range syntheticSeq {
		syntheticSeq[i] = pattern[i%len(pattern)]
	}
	session.referenceGenome.loader.chromosomes["chr21"] = syntheticSeq
	session.referenceGenome.loader.chromosomes["chr22"] = syntheticSeq

	return session
}

// TestTrainHMMModels_WithNoRealData exercises TrainHMMModels when the reference
// genome loader has no files to load. All chromosomes fail to load, the model
// trains with empty data (uniform probs) and is still marked trained=true.
func TestTrainHMMModels_WithNoRealData(t *testing.T) {
	session := newSessionWithTempDir(t)
	err := session.TrainHMMModels()
	if err != nil {
		t.Errorf("TrainHMMModels() returned unexpected error: %v", err)
	}
	// At least one model should have been stored.
	if len(session.trainedModels) == 0 {
		t.Error("expected at least one trained model after TrainHMMModels()")
	}
	for key, model := range session.trainedModels {
		if !model.IsTrained() {
			t.Errorf("model %q should be marked trained", key)
		}
	}
}

// TestTrainHMMModels_CreatesOutputDir verifies that TrainHMMModels creates the
// output directory specified in the config.
func TestTrainHMMModels_CreatesOutputDir(t *testing.T) {
	session := newSessionWithTempDir(t)
	// Remove the output dir so we verify it gets created.
	os.RemoveAll(session.config.OutputPath)

	err := session.TrainHMMModels()
	if err != nil {
		t.Errorf("TrainHMMModels() returned unexpected error: %v", err)
	}
	if _, statErr := os.Stat(session.config.OutputPath); os.IsNotExist(statErr) {
		t.Error("expected output directory to be created by TrainHMMModels()")
	}
}

// TestValidateModels_AfterTraining verifies that ValidateModels succeeds when
// called on a session that has been trained via TrainHMMModels.
// We use a seeded chromosome so that the model learns non-trivial probabilities
// producing balanced GC-content sequences.
func TestValidateModels_AfterTraining(t *testing.T) {
	session := newSessionWithSeededChromosome(t)

	if err := session.TrainHMMModels(); err != nil {
		t.Fatalf("TrainHMMModels() returned unexpected error: %v", err)
	}

	err := session.ValidateModels()
	if err != nil {
		t.Errorf("ValidateModels() returned unexpected error: %v", err)
	}
}

// TestValidateModels_EmptyModels verifies that ValidateModels returns nil (no
// error) when the trainedModels map is empty (nothing to validate).
func TestValidateModels_EmptyModels(t *testing.T) {
	session := newSessionWithTempDir(t)
	// trainedModels is already empty after construction.
	err := session.ValidateModels()
	if err != nil {
		t.Errorf("ValidateModels() with empty models returned unexpected error: %v", err)
	}
}

// TestRunTraining_WithSeededChromosome exercises the full pipeline: TrainHMMModels
// -> ValidateModels -> generateTrainingReport.
func TestRunTraining_WithSeededChromosome(t *testing.T) {
	session := newSessionWithSeededChromosome(t)

	err := session.RunTraining()
	if err != nil {
		t.Errorf("RunTraining() returned unexpected error: %v", err)
	}

	// Verify that training report was created.
	reportPath := session.config.OutputPath + "/training_report.txt"
	if _, statErr := os.Stat(reportPath); os.IsNotExist(statErr) {
		t.Error("expected training_report.txt to be created by RunTraining()")
	}
}

// TestRunTraining_BadOutputPath verifies that RunTraining fails gracefully when
// the output path cannot be created (e.g., nested under a file that exists).
func TestRunTraining_BadOutputPath(t *testing.T) {
	// Create a file at the path where the output directory would be, so MkdirAll fails.
	tmpDir := t.TempDir()
	blockingFile := tmpDir + "/blocked"
	if err := os.WriteFile(blockingFile, []byte("x"), 0600); err != nil {
		t.Fatalf("failed to create blocking file: %v", err)
	}

	session := NewTrainingSession("/tmp/nonexistent-ref")
	session.config = &training.TrainingConfig{
		KmerSizes:               []int{4},
		Populations:             []string{"EUR"},
		JSdivergenceThreshold:   0.01,
		ExpertAccuracyThreshold: 0.55,
		// Use a path like "/blockedfile/subdir" which can't be created.
		OutputPath: blockingFile + "/subdir",
	}
	session.startTime = time.Now()

	err := session.RunTraining()
	if err == nil {
		t.Error("expected RunTraining() to fail with uncreateable output path, got nil")
	}
}
