package genomics

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestNewEngine(t *testing.T) {
	config := types.SequencerConfig{
		SampleName:    "Test_Sample",
		Lane:          1,
		OutputDir:     "./test_output",
		ReadLength:    150,
		ReadsPerBatch: 1000,
		BatchInterval: "10m",
		TotalBatches:  2,
		PairedEnd:     true,
		BaseQuality:   30,
		QualityDecay:  0.1,
	}

	engine := NewEngine(config)
	if engine == nil {
		t.Fatal("Expected engine to be created")
	}
	if engine.config.SampleName != "Test_Sample" {
		t.Errorf("Expected sample name Test_Sample, got %s", engine.config.SampleName)
	}
	if engine.config.ReadLength != 150 {
		t.Errorf("Expected read length 150, got %d", engine.config.ReadLength)
	}
}

func TestEngineGenerate(t *testing.T) {
	tmpDir := t.TempDir()
	
	config := types.SequencerConfig{
		SampleName:    "Test_Sample",
		Lane:          1,
		OutputDir:     tmpDir,
		ReadLength:    50,
		ReadsPerBatch: 100,
		BatchInterval: "1s",
		TotalBatches:  2,
		PairedEnd:     true,
		BaseQuality:   30,
		QualityDecay:  0.1,
	}

	engine := NewEngine(config)
	err := engine.Generate()
	if err != nil {
		t.Fatalf("Expected generation to succeed, got error: %v", err)
	}

	// Check that output files were created
	files, err := os.ReadDir(tmpDir)
	if err != nil {
		t.Fatalf("Failed to read output directory: %v", err)
	}

	if len(files) == 0 {
		t.Error("Expected output files to be created")
	}

	// Check for FASTQ files
	fastqCount := 0
	for _, file := range files {
		if strings.HasSuffix(file.Name(), ".fastq") {
			fastqCount++
		}
	}

	expectedFiles := config.TotalBatches * 2 // R1 and R2 for paired-end
	if fastqCount != expectedFiles {
		t.Errorf("Expected %d FASTQ files, got %d", expectedFiles, fastqCount)
	}
}

func TestEngineGenerateSingleEnd(t *testing.T) {
	tmpDir := t.TempDir()
	
	config := types.SequencerConfig{
		SampleName:    "SingleEnd_Sample",
		Lane:          1,
		OutputDir:     tmpDir,
		ReadLength:    75,
		ReadsPerBatch: 50,
		BatchInterval: "1s",
		TotalBatches:  1,
		PairedEnd:     false, // Single-end
		BaseQuality:   25,
		QualityDecay:  0.15,
	}

	engine := NewEngine(config)
	err := engine.Generate()
	if err != nil {
		t.Fatalf("Expected generation to succeed, got error: %v", err)
	}

	// Check that only R1 files were created for single-end
	files, err := os.ReadDir(tmpDir)
	if err != nil {
		t.Fatalf("Failed to read output directory: %v", err)
	}

	fastqCount := 0
	r2Count := 0
	for _, file := range files {
		if strings.HasSuffix(file.Name(), ".fastq") {
			fastqCount++
			if strings.Contains(file.Name(), "_R2_") {
				r2Count++
			}
		}
	}

	if fastqCount != config.TotalBatches {
		t.Errorf("Expected %d FASTQ files for single-end, got %d", config.TotalBatches, fastqCount)
	}
	if r2Count != 0 {
		t.Errorf("Expected no R2 files for single-end, got %d", r2Count)
	}
}

func TestGenerateFASTQFile(t *testing.T) {
	tmpDir := t.TempDir()
	filename := filepath.Join(tmpDir, "test_R1_001.fastq")
	
	config := types.SequencerConfig{
		SampleName:   "Test_Sample",
		ReadLength:   100,
		ReadsPerBatch: 10,
		BaseQuality:  30,
		QualityDecay: 0.1,
	}

	engine := NewEngine(config)
	err := engine.generateFASTQFile(filename, 1, 1)
	if err != nil {
		t.Fatalf("Expected FASTQ generation to succeed, got error: %v", err)
	}

	// Check that file was created
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		t.Error("Expected FASTQ file to be created")
	}

	// Read and validate FASTQ content
	content, err := os.ReadFile(filename)
	if err != nil {
		t.Fatalf("Failed to read FASTQ file: %v", err)
	}

	lines := strings.Split(string(content), "\n")
	if len(lines) < 4 {
		t.Error("Expected at least 4 lines in FASTQ file")
	}

	// Check FASTQ format
	if !strings.HasPrefix(lines[0], "@") {
		t.Error("Expected first line to start with '@' (read header)")
	}
	if !strings.HasPrefix(lines[2], "+") {
		t.Error("Expected third line to start with '+' (quality header)")
	}

	// Check read length
	if len(lines[1]) != config.ReadLength {
		t.Errorf("Expected read length %d, got %d", config.ReadLength, len(lines[1]))
	}

	// Check quality string length matches sequence length
	if len(lines[3]) != len(lines[1]) {
		t.Error("Expected quality string length to match sequence length")
	}

	// Validate DNA sequence (should only contain A, T, C, G)
	for _, char := range lines[1] {
		if char != 'A' && char != 'T' && char != 'C' && char != 'G' {
			t.Errorf("Expected only A, T, C, G in sequence, found %c", char)
		}
	}
}

func TestFASTQFormat(t *testing.T) {
	tmpDir := t.TempDir()
	filename := filepath.Join(tmpDir, "test_format.fastq")
	
	config := types.SequencerConfig{
		SampleName:   "Test_Sample",
		ReadLength:   50,
		ReadsPerBatch: 5,
		BaseQuality:  30,
		QualityDecay: 0.1,
	}

	engine := NewEngine(config)
	err := engine.generateFASTQFile(filename, 1, 1)
	if err != nil {
		t.Fatalf("Expected FASTQ generation to succeed, got error: %v", err)
	}

	// Read and validate basic FASTQ format
	content, err := os.ReadFile(filename)
	if err != nil {
		t.Fatalf("Failed to read FASTQ file: %v", err)
	}

	lines := strings.Split(string(content), "\n")
	
	// Should have at least 4 lines per read (header, sequence, +, quality)
	if len(lines) < 4 {
		t.Error("Expected at least 4 lines in FASTQ file")
	}

	// Check basic FASTQ format
	if len(lines) >= 4 {
		if !strings.HasPrefix(lines[0], "@") {
			t.Error("Expected first line to start with '@' (read header)")
		}
		if !strings.HasPrefix(lines[2], "+") {
			t.Error("Expected third line to start with '+' (quality header)")
		}
		
		// Check read length
		if len(lines[1]) != config.ReadLength {
			t.Errorf("Expected read length %d, got %d", config.ReadLength, len(lines[1]))
		}
	}
}

func TestEngineErrorHandling(t *testing.T) {
	// Test with invalid output directory
	config := types.SequencerConfig{
		SampleName:    "Test_Sample",
		Lane:          1,
		OutputDir:     "/invalid/path/that/does/not/exist",
		ReadLength:    50,
		ReadsPerBatch: 10,
		TotalBatches:  1,
		PairedEnd:     false,
	}

	engine := NewEngine(config)
	err := engine.Generate()
	if err == nil {
		t.Error("Expected error for invalid output directory")
	}
}

func TestEngineConfigValidation(t *testing.T) {
	// Test with zero read length
	config := types.SequencerConfig{
		SampleName:    "Test_Sample",
		Lane:          1,
		OutputDir:     "./test_output",
		ReadLength:    0, // Invalid
		ReadsPerBatch: 10,
		TotalBatches:  1,
	}

	engine := NewEngine(config)
	err := engine.Generate()
	if err == nil {
		t.Error("Expected error for zero read length")
	}
}

func TestFASTQFileNaming(t *testing.T) {
	tmpDir := t.TempDir()
	
	config := types.SequencerConfig{
		SampleName:    "Sample_123",
		Lane:          2,
		OutputDir:     tmpDir,
		ReadLength:    50,
		ReadsPerBatch: 10,
		TotalBatches:  2,
		PairedEnd:     true,
	}

	engine := NewEngine(config)
	err := engine.Generate()
	if err != nil {
		t.Fatalf("Expected generation to succeed, got error: %v", err)
	}

	files, err := os.ReadDir(tmpDir)
	if err != nil {
		t.Fatalf("Failed to read output directory: %v", err)
	}

	// Check naming convention: SampleName_S#_L00#_R#_00#.fastq
	expectedNames := []string{
		"Sample_123_S1_L002_R1_001.fastq",
		"Sample_123_S1_L002_R2_001.fastq",
		"Sample_123_S1_L002_R1_002.fastq",
		"Sample_123_S1_L002_R2_002.fastq",
	}

	fileNames := make([]string, len(files))
	for i, file := range files {
		fileNames[i] = file.Name()
	}

	for _, expected := range expectedNames {
		found := false
		for _, actual := range fileNames {
			if actual == expected {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("Expected file %s not found. Available files: %v", expected, fileNames)
		}
	}
}

func TestLongReadGeneration(t *testing.T) {
	tmpDir := t.TempDir()
	
	config := types.SequencerConfig{
		SampleName:    "LongRead_Sample",
		Lane:          1,
		OutputDir:     tmpDir,
		ReadLength:    300, // Long read like MiSeq v3
		ReadsPerBatch: 5,
		TotalBatches:  1,
		PairedEnd:     true,
		BaseQuality:   30,
		QualityDecay:  0.12, // Higher decay for longer reads
	}

	engine := NewEngine(config)
	err := engine.Generate()
	if err != nil {
		t.Fatalf("Expected generation to succeed, got error: %v", err)
	}

	// Read the generated file and check read length
	files, err := os.ReadDir(tmpDir)
	if err != nil {
		t.Fatalf("Failed to read output directory: %v", err)
	}

	for _, file := range files {
		if strings.HasSuffix(file.Name(), ".fastq") {
			content, err := os.ReadFile(filepath.Join(tmpDir, file.Name()))
			if err != nil {
				t.Fatalf("Failed to read FASTQ file: %v", err)
			}

			lines := strings.Split(string(content), "\n")
			if len(lines) >= 2 {
				// Check first sequence length
				if len(lines[1]) != 300 {
					t.Errorf("Expected read length 300, got %d", len(lines[1]))
				}
			}
			break // Only check first file
		}
	}
}

// makeEngineTestSequence returns a deterministic ATGC byte slice for HMM training.
func makeEngineTestSequence(n int) []byte {
	bases := []byte{'A', 'T', 'G', 'C'}
	s := make([]byte, n)
	for i := range s {
		s[i] = bases[(i+i/7)%4]
	}
	return s
}

// TestEngine_HMMModel_GeneratesValidSequences verifies that when HMMModelPath is set
// to a valid model, the engine uses it and produces valid FASTQ output.
func TestEngine_HMMModel_GeneratesValidSequences(t *testing.T) {
	// Train a small AdaptiveHMMModel and save it.
	model := authenticity.NewAdaptiveHMMModel(3)
	if model == nil {
		t.Fatal("failed to create AdaptiveHMMModel")
	}
	if err := model.TrainFromSequence(makeEngineTestSequence(5000)); err != nil {
		t.Fatalf("TrainFromSequence failed: %v", err)
	}

	modelPath := filepath.Join(t.TempDir(), "test.hmm.json")
	if err := model.Save(modelPath); err != nil {
		t.Fatalf("model.Save failed: %v", err)
	}

	tmpDir := t.TempDir()
	config := types.SequencerConfig{
		SampleName:    "HMMTest_Sample",
		Lane:          1,
		OutputDir:     tmpDir,
		ReadLength:    50,
		ReadsPerBatch: 5,
		TotalBatches:  1,
		PairedEnd:     false,
		BaseQuality:   30,
		QualityDecay:  0.1,
		HMMModelPath:  modelPath,
	}

	engine := NewEngine(config)
	if engine.hmmModel == nil {
		t.Fatal("expected hmmModel to be loaded after NewEngine with valid HMMModelPath")
	}
	if err := engine.Generate(); err != nil {
		t.Fatalf("Generate() failed: %v", err)
	}

	files, err := os.ReadDir(tmpDir)
	if err != nil {
		t.Fatalf("ReadDir failed: %v", err)
	}
	var fastqPath string
	for _, f := range files {
		if strings.HasSuffix(f.Name(), ".fastq") {
			fastqPath = filepath.Join(tmpDir, f.Name())
			break
		}
	}
	if fastqPath == "" {
		t.Fatal("no FASTQ file generated")
	}

	content, err := os.ReadFile(fastqPath)
	if err != nil {
		t.Fatalf("ReadFile failed: %v", err)
	}
	lines := strings.Split(string(content), "\n")
	if len(lines) < 4 {
		t.Fatalf("too few lines in FASTQ: %d", len(lines))
	}

	seq1 := lines[1]
	qual1 := lines[3]
	if len(seq1) != config.ReadLength {
		t.Errorf("sequence length: want %d, got %d", config.ReadLength, len(seq1))
	}
	if len(qual1) != config.ReadLength {
		t.Errorf("quality length: want %d, got %d", config.ReadLength, len(qual1))
	}
	for _, c := range seq1 {
		if c != 'A' && c != 'T' && c != 'G' && c != 'C' {
			t.Errorf("invalid base %c in HMM-generated sequence", c)
		}
	}
}

// TestEngine_HMMModel_FallsBackOnMissingPath verifies silent fallback when
// HMMModelPath points to a non-existent file: hmmModel is nil and Generate() succeeds.
func TestEngine_HMMModel_FallsBackOnMissingPath(t *testing.T) {
	tmpDir := t.TempDir()
	config := types.SequencerConfig{
		SampleName:    "Fallback_Sample",
		Lane:          1,
		OutputDir:     tmpDir,
		ReadLength:    50,
		ReadsPerBatch: 5,
		TotalBatches:  1,
		PairedEnd:     false,
		BaseQuality:   30,
		QualityDecay:  0.1,
		HMMModelPath:  "/nonexistent/path/model.hmm.json",
	}

	engine := NewEngine(config)
	if engine.hmmModel != nil {
		t.Error("expected hmmModel to be nil after missing-path fallback")
	}
	if err := engine.Generate(); err != nil {
		t.Errorf("Generate() should succeed after fallback, got: %v", err)
	}
}

// TestEngine_HMMModel_FallsBackOnCorruptFile verifies silent fallback when
// HMMModelPath points to a corrupt JSON file: hmmModel is nil and Generate() succeeds.
func TestEngine_HMMModel_FallsBackOnCorruptFile(t *testing.T) {
	corruptPath := filepath.Join(t.TempDir(), "corrupt.hmm.json")
	if err := os.WriteFile(corruptPath, []byte("not valid json {{{"), 0600); err != nil {
		t.Fatalf("WriteFile failed: %v", err)
	}

	tmpDir := t.TempDir()
	config := types.SequencerConfig{
		SampleName:    "Corrupt_Sample",
		Lane:          1,
		OutputDir:     tmpDir,
		ReadLength:    50,
		ReadsPerBatch: 5,
		TotalBatches:  1,
		PairedEnd:     false,
		BaseQuality:   30,
		QualityDecay:  0.1,
		HMMModelPath:  corruptPath,
	}

	engine := NewEngine(config)
	if engine.hmmModel != nil {
		t.Error("expected hmmModel to be nil after corrupt-file fallback")
	}
	if err := engine.Generate(); err != nil {
		t.Errorf("Generate() should succeed after fallback, got: %v", err)
	}
}

// --- resolveModelPath unit tests ---

func TestResolveModelPath_ExplicitPath(t *testing.T) {
	cfg := types.SequencerConfig{HMMModelPath: "/models/explicit.hmm.json"}
	if got := resolveModelPath(cfg); got != "/models/explicit.hmm.json" {
		t.Errorf("want /models/explicit.hmm.json, got %q", got)
	}
}

func TestResolveModelPath_PopulationAndDir(t *testing.T) {
	cfg := types.SequencerConfig{Population: "YRI", ModelDir: "/models"}
	want := filepath.Join("/models", "YRI.hmm.json")
	if got := resolveModelPath(cfg); got != want {
		t.Errorf("want %q, got %q", want, got)
	}
}

func TestResolveModelPath_ExplicitTakesPriority(t *testing.T) {
	cfg := types.SequencerConfig{
		HMMModelPath: "/models/explicit.hmm.json",
		Population:   "YRI",
		ModelDir:     "/models",
	}
	if got := resolveModelPath(cfg); got != "/models/explicit.hmm.json" {
		t.Errorf("HMMModelPath should take priority; got %q", got)
	}
}

func TestResolveModelPath_PopulationOnlyNoDir(t *testing.T) {
	cfg := types.SequencerConfig{Population: "YRI"} // ModelDir empty
	if got := resolveModelPath(cfg); got != "" {
		t.Errorf("expected empty path without ModelDir, got %q", got)
	}
}

func TestResolveModelPath_Empty(t *testing.T) {
	cfg := types.SequencerConfig{}
	if got := resolveModelPath(cfg); got != "" {
		t.Errorf("expected empty path for zero config, got %q", got)
	}
}

// --- engine integration tests for population-to-model resolution ---

// TestEngine_PopulationModel_ResolvesFromModelDir verifies that when Population and
// ModelDir are set, the engine loads <ModelDir>/<Population>.hmm.json.
func TestEngine_PopulationModel_ResolvesFromModelDir(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(3)
	if model == nil {
		t.Fatal("failed to create AdaptiveHMMModel")
	}
	if err := model.TrainFromSequence(makeEngineTestSequence(5000)); err != nil {
		t.Fatalf("TrainFromSequence: %v", err)
	}

	modelDir := t.TempDir()
	if err := model.Save(filepath.Join(modelDir, "YRI.hmm.json")); err != nil {
		t.Fatalf("model.Save: %v", err)
	}

	tmpDir := t.TempDir()
	config := types.SequencerConfig{
		SampleName:    "PopTest_Sample",
		Lane:          1,
		OutputDir:     tmpDir,
		ReadLength:    50,
		ReadsPerBatch: 5,
		TotalBatches:  1,
		PairedEnd:     false,
		BaseQuality:   30,
		QualityDecay:  0.1,
		Population:    "YRI",
		ModelDir:      modelDir,
	}

	engine := NewEngine(config)
	if engine.hmmModel == nil {
		t.Fatal("expected hmmModel to be loaded via Population+ModelDir")
	}
	if err := engine.Generate(); err != nil {
		t.Fatalf("Generate() failed: %v", err)
	}
}

// TestEngine_PopulationModel_HMMModelPathTakesPriority verifies that an explicit
// HMMModelPath is used even when Population+ModelDir are also set.
func TestEngine_PopulationModel_HMMModelPathTakesPriority(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(3)
	if model == nil {
		t.Fatal("failed to create AdaptiveHMMModel")
	}
	if err := model.TrainFromSequence(makeEngineTestSequence(5000)); err != nil {
		t.Fatalf("TrainFromSequence: %v", err)
	}

	modelDir := t.TempDir()
	explicitPath := filepath.Join(modelDir, "explicit.hmm.json")
	if err := model.Save(explicitPath); err != nil {
		t.Fatalf("model.Save explicit: %v", err)
	}
	if err := model.Save(filepath.Join(modelDir, "YRI.hmm.json")); err != nil {
		t.Fatalf("model.Save pop: %v", err)
	}

	config := types.SequencerConfig{
		SampleName:    "PriorityTest",
		Lane:          1,
		OutputDir:     t.TempDir(),
		ReadLength:    50,
		ReadsPerBatch: 5,
		TotalBatches:  1,
		BaseQuality:   30,
		QualityDecay:  0.1,
		HMMModelPath:  explicitPath,
		Population:    "YRI",
		ModelDir:      modelDir,
	}

	engine := NewEngine(config)
	if engine.hmmModel == nil {
		t.Fatal("expected hmmModel to be loaded")
	}
	if err := engine.Generate(); err != nil {
		t.Fatalf("Generate() failed: %v", err)
	}
}

// TestEngine_PopulationModel_FallsBackWhenDirMissing verifies that a non-existent
// ModelDir causes a silent fallback: hmmModel is nil and Generate() succeeds.
func TestEngine_PopulationModel_FallsBackWhenDirMissing(t *testing.T) {
	config := types.SequencerConfig{
		SampleName:    "FallbackPop_Sample",
		Lane:          1,
		OutputDir:     t.TempDir(),
		ReadLength:    50,
		ReadsPerBatch: 5,
		TotalBatches:  1,
		BaseQuality:   30,
		QualityDecay:  0.1,
		Population:    "YRI",
		ModelDir:      "/nonexistent/model/dir",
	}

	engine := NewEngine(config)
	if engine.hmmModel != nil {
		t.Error("expected hmmModel to be nil when ModelDir does not exist")
	}
	if err := engine.Generate(); err != nil {
		t.Errorf("Generate() should succeed after fallback, got: %v", err)
	}
}