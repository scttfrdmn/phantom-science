package timing

import (
	"context"
	"testing"
	"time"
	
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestGetSequencerTiming(t *testing.T) {
	testCases := []struct {
		model        string
		flowCellType string
	}{
		{"NovaSeq 6000", "S1"},
		{"NovaSeq 6000", "S2"},
		{"NextSeq 2000", "P1"},
		{"MiSeq v3", "v3"},
		{"Unknown Platform", "Unknown"},
	}
	
	for _, tc := range testCases {
		timing := GetSequencerTiming(tc.model, tc.flowCellType)
		
		if timing.TotalRunTimeHours <= 0 {
			t.Errorf("Platform %s/%s should have positive run time, got %d", tc.model, tc.flowCellType, timing.TotalRunTimeHours)
		}
		
		if timing.ClusterGenMinutes < 0 {
			t.Errorf("Platform %s/%s should have non-negative cluster generation time, got %d", tc.model, tc.flowCellType, timing.ClusterGenMinutes)
		}
		
		if timing.SequencingMinutes <= 0 {
			t.Errorf("Platform %s/%s should have positive sequencing time, got %d", tc.model, tc.flowCellType, timing.SequencingMinutes)
		}
	}
}

func TestNewRunSimulator(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		TotalBatches:   5,
		PairedEnd:      true,
	}
	
	compressions := []TimeCompression{RealTime, Compressed10x, Compressed100x, Compressed1000x}
	
	for _, compression := range compressions {
		simulator := NewRunSimulator(config, compression)
		if simulator == nil {
			t.Errorf("NewRunSimulator should not return nil for compression %v", compression)
		}
	}
}

func TestRunSimulatorEstimateFileSize(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		TotalBatches:   2,
		PairedEnd:      true,
	}
	
	simulator := NewRunSimulator(config, Compressed100x)
	fileSize := simulator.EstimateFileSize(1000) // 1000 reads per file
	
	if fileSize <= 0 {
		t.Errorf("EstimateFileSize should return positive value, got %d", fileSize)
	}
	
	// File size should be reasonable for the given parameters
	// With 1000 reads * 150bp, expect at least a few KB
	expectedMinSize := int64(1000) // 1KB minimum
	if fileSize < expectedMinSize {
		t.Errorf("File size %d seems too small for given parameters", fileSize)
	}
}

func TestRunSimulatorEstimateTotalOutput(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq v3",
		ReadLength:     150,
		ReadsPerBatch:  500,
		TotalBatches:   3,
	}
	
	simulator := NewRunSimulator(config, Compressed10x)
	fileCount, totalBytes := simulator.EstimateTotalOutput()
	
	if fileCount <= 0 {
		t.Error("EstimateTotalOutput should return positive file count")
	}
	
	if totalBytes <= 0 {
		t.Error("EstimateTotalOutput should return positive byte count")
	}
}

func TestRunSimulatorWithDifferentPlatforms(t *testing.T) {
	platforms := []string{"NovaSeq 6000", "NextSeq 2000", "MiSeq v3"}
	
	for _, platform := range platforms {
		config := types.SequencerConfig{
			SequencerModel: platform,
			ReadLength:     100,
			ReadsPerBatch:  100,
			TotalBatches:   1,
		}
		
		simulator := NewRunSimulator(config, Compressed10x)
		if simulator == nil {
			t.Errorf("Should create simulator for platform %s", platform)
			continue
		}
		
		fileSize := simulator.EstimateFileSize(100)
		if fileSize <= 0 {
			t.Errorf("Platform %s should estimate positive file size", platform)
		}
		
		fileCount, totalBytes := simulator.EstimateTotalOutput()
		if fileCount <= 0 || totalBytes <= 0 {
			t.Errorf("Platform %s should estimate positive output", platform)
		}
	}
}

func TestTimeCompressionConstants(t *testing.T) {
	// Test that compression constants are properly defined
	if RealTime != 1 {
		t.Errorf("RealTime should be 1, got %d", RealTime)
	}
	
	if Compressed10x != 10 {
		t.Errorf("Compressed10x should be 10, got %d", Compressed10x)
	}
	
	if Compressed100x != 100 {
		t.Errorf("Compressed100x should be 100, got %d", Compressed100x)
	}
	
	if Compressed1000x != 1000 {
		t.Errorf("Compressed1000x should be 1000, got %d", Compressed1000x)
	}
}

func TestGetSequencerTimingKnownPlatforms(t *testing.T) {
	testCases := []struct {
		model        string
		flowCellType string
		expectedHours int
	}{
		{"NovaSeq 6000", "S1", 13},
		{"NovaSeq 6000", "S2", 36},
		{"NovaSeq 6000", "S4", 44},
		{"NextSeq 2000", "P3", 18},
		{"NextSeq 2000", "P4", 24},
		{"MiSeq", "v3", 65},
		{"MiSeq", "v2", 27},
		{"PacBio Revio", "SMRT", 24},
		{"DNBSEQ-T20x2", "DNB", 48},
		{"Element AVITI", "Dual", 38},
		{"PromethION 24", "R10.4.1", 72},
		{"PromethION 24", "R9.4.1", 72},
	}
	
	for _, tc := range testCases {
		timing := GetSequencerTiming(tc.model, tc.flowCellType)
		
		if timing.Model != tc.model {
			t.Errorf("Expected model %s, got %s", tc.model, timing.Model)
		}
		
		if timing.FlowCellType != tc.flowCellType {
			t.Errorf("Expected flow cell %s, got %s", tc.flowCellType, timing.FlowCellType)
		}
		
		if timing.TotalRunTimeHours != tc.expectedHours {
			t.Errorf("Platform %s/%s: expected %d hours, got %d", 
				tc.model, tc.flowCellType, tc.expectedHours, timing.TotalRunTimeHours)
		}
		
		if timing.BytesPerSecond <= 0 {
			t.Errorf("Platform %s/%s should have positive throughput", tc.model, tc.flowCellType)
		}
	}
}

func TestGetSequencerTimingDefaultFallback(t *testing.T) {
	// Test unknown platform falls back to defaults
	timing := GetSequencerTiming("Unknown Platform", "Unknown FlowCell")
	
	if timing.Model != "Unknown Platform" {
		t.Errorf("Expected model to be preserved, got %s", timing.Model)
	}
	
	if timing.FlowCellType != "Unknown FlowCell" {
		t.Errorf("Expected flow cell to be preserved, got %s", timing.FlowCellType)
	}
	
	if timing.TotalRunTimeHours != 24 {
		t.Errorf("Expected default 24 hours, got %d", timing.TotalRunTimeHours)
	}
	
	if timing.BytesPerSecond != 10*1024*1024 {
		t.Errorf("Expected default 10MB/s throughput, got %d", timing.BytesPerSecond)
	}
}

func TestRunSimulatorCallbacks(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		FlowCellType:   "S1",
		ReadLength:     100,
		ReadsPerBatch:  100,
		TotalBatches:   2,
		PairedEnd:      true,
		SampleName:     "TestSample",
		Lane:           1,
	}
	
	simulator := NewRunSimulator(config, Compressed1000x)
	
	// Test callback setup
	var fileCallbacks []string
	var progressCallbacks []string
	var errorCallbacks []error
	
	simulator.FileCallback = func(filename string, batchNum int, readNum int) {
		fileCallbacks = append(fileCallbacks, filename)
	}
	
	simulator.ProgressCallback = func(stage string, percent float64, eta time.Duration) {
		progressCallbacks = append(progressCallbacks, stage)
	}
	
	simulator.ErrorCallback = func(err error) {
		errorCallbacks = append(errorCallbacks, err)
	}
	
	// Run a very short simulation
	ctx, cancel := context.WithTimeout(context.Background(), 2*time.Second)
	defer cancel()
	
	err := simulator.Run(ctx)
	if err != nil && err != context.DeadlineExceeded {
		t.Errorf("Run should complete or timeout gracefully, got: %v", err)
	}
	
	// Verify callbacks were called
	if len(progressCallbacks) == 0 {
		t.Error("Progress callback should have been called")
	}
}

func TestRunSimulatorFileNaming(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		ReadsPerBatch:  500,
		TotalBatches:   2,
		PairedEnd:      true,
		SampleName:     "TestSample",
		Lane:           3,
	}
	
	simulator := NewRunSimulator(config, Compressed1000x)
	
	var generatedFiles []string
	simulator.FileCallback = func(filename string, batchNum int, readNum int) {
		generatedFiles = append(generatedFiles, filename)
	}
	
	// Run with timeout to generate some files
	ctx, cancel := context.WithTimeout(context.Background(), 500*time.Millisecond)
	defer cancel()
	
	simulator.Run(ctx)
	
	// Check file naming pattern
	for _, filename := range generatedFiles {
		if filename == "" {
			t.Error("Generated filename should not be empty")
		}
		
		// Should contain sample name
		if len(filename) < len("TestSample") {
			t.Errorf("Filename %s seems too short", filename)
		}
	}
}

func TestRunStages(t *testing.T) {
	stages := []RunStage{
		StageInitialization,
		StageClusterGen,
		StageSequencing,
		StageBaseCalling,
		StageComplete,
	}
	
	expectedValues := []string{
		"Initialization",
		"Cluster Generation", 
		"Sequencing",
		"Base Calling",
		"Complete",
	}
	
	for i, stage := range stages {
		if string(stage) != expectedValues[i] {
			t.Errorf("Stage %d: expected %s, got %s", i, expectedValues[i], string(stage))
		}
	}
}

func TestRunProgressStructure(t *testing.T) {
	progress := RunProgress{
		Stage:          StageSequencing,
		Percent:        75.5,
		ETA:            30 * time.Minute,
		FilesGenerated: 15,
		TotalFiles:     20,
		BytesGenerated: 1024 * 1024 * 500, // 500MB
		TotalBytes:     1024 * 1024 * 1000, // 1GB
		CurrentFile:    "sample_S1_L001_R1_001.fastq",
	}
	
	if progress.Stage != StageSequencing {
		t.Error("RunProgress fields should be settable")
	}
	
	if progress.Percent != 75.5 {
		t.Error("RunProgress numeric fields should be settable")
	}
	
	if progress.CurrentFile != "sample_S1_L001_R1_001.fastq" {
		t.Error("RunProgress string fields should be settable")
	}
}

func TestGetCompressionDescription(t *testing.T) {
	testCases := []struct {
		compression TimeCompression
		realHours   int
		expectMinutes bool
	}{
		{Compressed1000x, 24, true},  // 24 hours -> 1.44 minutes
		{Compressed100x, 24, true},   // 24 hours -> 14.4 minutes
		{Compressed10x, 24, false},   // 24 hours -> 2.4 hours
		{RealTime, 24, false},        // 24 hours -> 24 hours
	}
	
	for _, tc := range testCases {
		desc := GetCompressionDescription(tc.compression, tc.realHours)
		
		if desc == "" {
			t.Errorf("Description should not be empty for compression %v", tc.compression)
		}
		
		if tc.expectMinutes && len(desc) < 10 {
			t.Errorf("Description seems too short: %s", desc)
		}
		
		// Should contain compression ratio
		found := false
		for _, char := range desc {
			if char >= '0' && char <= '9' {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("Description should contain numeric compression info: %s", desc)
		}
	}
}

func TestRunSimulatorContextCancellation(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		FlowCellType:   "S1",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		TotalBatches:   10,
		PairedEnd:      true,
		SampleName:     "CancelTest",
		Lane:           1,
	}
	
	simulator := NewRunSimulator(config, Compressed100x)
	
	// Cancel context quickly
	ctx, cancel := context.WithCancel(context.Background())
	cancel() // Cancel immediately
	
	err := simulator.Run(ctx)
	if err != context.Canceled {
		t.Errorf("Expected context.Canceled, got %v", err)
	}
}

func TestEstimateFileSizeVariousConfigs(t *testing.T) {
	testCases := []struct {
		readLength    int
		readsPerFile  int
		expectedMinSize int64
	}{
		{50, 100, 1000},     // Short reads, few reads
		{150, 1000, 150000}, // Standard reads
		{300, 500, 100000},  // Long reads
	}
	
	for _, tc := range testCases {
		config := types.SequencerConfig{
			SequencerModel: "NovaSeq 6000",
			ReadLength:     tc.readLength,
			ReadsPerBatch:  1000,
			TotalBatches:   1,
		}
		
		simulator := NewRunSimulator(config, RealTime)
		size := simulator.EstimateFileSize(tc.readsPerFile)
		
		if size < tc.expectedMinSize {
			t.Errorf("Read length %d, reads %d: expected at least %d bytes, got %d",
				tc.readLength, tc.readsPerFile, tc.expectedMinSize, size)
		}
	}
}

func TestEstimateTotalOutputPairedVsSingle(t *testing.T) {
	baseConfig := types.SequencerConfig{
		SequencerModel: "NextSeq 2000",
		FlowCellType:   "P3",
		ReadLength:     100,
		ReadsPerBatch:  500,
		TotalBatches:   3,
	}
	
	// Test single-end
	singleConfig := baseConfig
	singleConfig.PairedEnd = false
	
	singleSim := NewRunSimulator(singleConfig, RealTime)
	singleFiles, singleBytes := singleSim.EstimateTotalOutput()
	
	// Test paired-end
	pairedConfig := baseConfig
	pairedConfig.PairedEnd = true
	
	pairedSim := NewRunSimulator(pairedConfig, RealTime)
	pairedFiles, pairedBytes := pairedSim.EstimateTotalOutput()
	
	// Paired-end should generate twice as many files
	if pairedFiles != singleFiles*2 {
		t.Errorf("Paired-end should generate 2x files: single=%d, paired=%d", singleFiles, pairedFiles)
	}
	
	// Paired-end should generate twice as much data
	if pairedBytes != singleBytes*2 {
		t.Errorf("Paired-end should generate 2x data: single=%d, paired=%d", singleBytes, pairedBytes)
	}
}

// TestRunSimulator_PairedEnd_ContextCancellation verifies that a cancelled context
// is honoured during the R2 ticker wait and does not hang (#114).
func TestRunSimulator_PairedEnd_ContextCancellation(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		FlowCellType:   "S4",
		TotalBatches:   20,
		ReadsPerBatch:  1000,
		ReadLength:     150,
		PairedEnd:      true,
		SampleName:     "cancel_test",
		Lane:           1,
	}

	sim := NewRunSimulator(config, RealTime*100) // fast but non-zero compression

	// Cancel the context almost immediately — well before the simulation finishes.
	ctx, cancel := context.WithTimeout(context.Background(), 50*time.Millisecond)
	defer cancel()

	done := make(chan error, 1)
	go func() {
		done <- sim.Run(ctx)
	}()

	select {
	case err := <-done:
		// Should have returned a context error, not hung.
		if err == nil {
			t.Error("expected context cancellation error, got nil")
		}
	case <-time.After(5 * time.Second):
		t.Error("Run did not return after context cancellation (goroutine blocked on bare ticker receive)")
	}
}