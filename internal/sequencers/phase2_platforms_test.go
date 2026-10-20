package sequencers

import (
	"os"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/timing"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestElementAVITIPlatform(t *testing.T) {
	// Create temporary output directory
	tempDir, err := os.MkdirTemp("", "element_aviti_test")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tempDir)

	// Test Element AVITI configuration
	config := types.SequencerConfig{
		SequencerModel:   "Element AVITI",
		FlowCellType:     "Dual",
		ChemistryVersion: "Avidity_SBB_v1",
		InstrumentID:     "AVITI001",
		RunID:            "240712_AVITI001_0001",
		FlowCellID:       "AVITI24071201",
		SampleName:       "test_aviti_sample",
		Lane:             1,
		OutputDir:        tempDir,
		ReadLength:       150,
		ReadsPerBatch:    1000,
		BatchInterval:    "30s",
		TotalBatches:     2,
		PairedEnd:        true,
		BaseQuality:      44,
		QualityDecay:     0.02,
		Q30Percentage:    96.2,
	}

	// Test timing profile retrieval
	timingProfile := timing.GetSequencerTiming(config.SequencerModel, config.FlowCellType)
	if timingProfile.Model != "Element AVITI" {
		t.Errorf("Expected Element AVITI timing model, got %s", timingProfile.Model)
	}
	if timingProfile.TotalRunTimeHours != 38 {
		t.Errorf("Expected 38 hour run time, got %d", timingProfile.TotalRunTimeHours)
	}
	if timingProfile.FilesPerHour != 8 {
		t.Errorf("Expected 8 files per hour, got %d", timingProfile.FilesPerHour)
	}

	// Test simulation configuration (without running full simulation)
	simulator := timing.NewRunSimulator(config, timing.Compressed1000x)
	if simulator.Timing.Model != "Element AVITI" {
		t.Errorf("Expected Element AVITI in simulator, got %s", simulator.Timing.Model)
	}

	// Verify file estimates
	totalFiles, totalBytes := simulator.EstimateTotalOutput()
	expectedFiles := 4 // 2 batches × 2 reads (paired-end)
	if totalFiles != expectedFiles {
		t.Errorf("Expected %d files estimate, got %d", expectedFiles, totalFiles)
	}

	if totalBytes <= 0 {
		t.Errorf("Expected positive byte estimate, got %d", totalBytes)
	}

	t.Logf("Element AVITI simulation configured: %d files, %d bytes estimated", totalFiles, totalBytes)
}

func TestPromethION24Platform(t *testing.T) {
	// Create temporary output directory
	tempDir, err := os.MkdirTemp("", "promethion24_test")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tempDir)

	// Test PromethION 24 configuration with R10.4.1 chemistry
	config := types.SequencerConfig{
		SequencerModel:   "PromethION 24",
		FlowCellType:     "R10.4.1",
		ChemistryVersion: "R10.4.1",
		InstrumentID:     "PROM001",
		RunID:            "240712_PROM001_0001",
		FlowCellID:       "PAK12345",
		SampleName:       "test_nanopore_sample",
		Lane:             1,
		OutputDir:        tempDir,
		ReadLength:       25000, // 25kb average for nanopore
		ReadsPerBatch:    100,   // Fewer reads due to ultra-long lengths
		BatchInterval:    "2m",  // Real-time streaming
		TotalBatches:     2,
		PairedEnd:        false, // Nanopore is single-end
		BaseQuality:      20,    // Q20 typical for nanopore SUP
		QualityDecay:     0.0,   // No positional decay for nanopore
	}

	// Test timing profile retrieval
	timingProfile := timing.GetSequencerTiming(config.SequencerModel, config.FlowCellType)
	if timingProfile.Model != "PromethION 24" {
		t.Errorf("Expected PromethION 24 timing model, got %s", timingProfile.Model)
	}
	if timingProfile.TotalRunTimeHours != 72 {
		t.Errorf("Expected 72 hour run time, got %d", timingProfile.TotalRunTimeHours)
	}
	if timingProfile.BaseCallingMinutes != 0 {
		t.Errorf("Expected 0 base calling minutes (real-time), got %d", timingProfile.BaseCallingMinutes)
	}
	if timingProfile.FilesPerHour != 24 {
		t.Errorf("Expected 24 files per hour (24 flow cells), got %d", timingProfile.FilesPerHour)
	}

	// Test simulation configuration (without running full simulation)
	simulator := timing.NewRunSimulator(config, timing.Compressed1000x)
	if simulator.Timing.Model != "PromethION 24" {
		t.Errorf("Expected PromethION 24 in simulator, got %s", simulator.Timing.Model)
	}

	// Verify file estimates for nanopore platform
	totalFiles, totalBytes := simulator.EstimateTotalOutput()
	expectedFiles := 2 // 2 batches × 1 read (single-end nanopore)
	if totalFiles != expectedFiles {
		t.Errorf("Expected %d files estimate, got %d", expectedFiles, totalFiles)
	}

	if totalBytes <= 0 {
		t.Errorf("Expected positive byte estimate, got %d", totalBytes)
	}

	t.Logf("PromethION 24 simulation configured: %d files, %d bytes estimated", totalFiles, totalBytes)
}

func TestPromethION24LegacyChemistry(t *testing.T) {
	// Test with R9.4.1 legacy chemistry
	config := types.SequencerConfig{
		SequencerModel:   "PromethION 24",
		FlowCellType:     "R9.4.1",
		ChemistryVersion: "R9.4.1",
		InstrumentID:     "PROM001",
		RunID:            "240712_PROM001_R9_0001",
		FlowCellID:       "FLO94123",
		SampleName:       "test_nanopore_r9_sample",
		Lane:             1,
		OutputDir:        "/tmp",
		ReadLength:       15000, // Shorter average for R9.4.1
		ReadsPerBatch:    200,
		BatchInterval:    "3m",
		TotalBatches:     1,
		PairedEnd:        false,
		BaseQuality:      12, // Lower quality for R9.4.1
		QualityDecay:     0.0,
	}

	// Test timing profile retrieval for legacy chemistry
	timingProfile := timing.GetSequencerTiming(config.SequencerModel, config.FlowCellType)
	if timingProfile.Model != "PromethION 24" {
		t.Errorf("Expected PromethION 24 timing model, got %s", timingProfile.Model)
	}
	if timingProfile.FlowCellType != "R9.4.1" {
		t.Errorf("Expected R9.4.1 flow cell type, got %s", timingProfile.FlowCellType)
	}
	if timingProfile.FilesPerHour != 20 {
		t.Errorf("Expected 20 files per hour for R9.4.1 (lower than R10.4.1), got %d", timingProfile.FilesPerHour)
	}
	if timingProfile.BytesPerSecond != 600*1024*1024 {
		t.Errorf("Expected 600 MB/s for R9.4.1, got %d", timingProfile.BytesPerSecond)
	}
}

func TestPhase2PlatformCoverage(t *testing.T) {
	// Verify all Phase 2 platforms are supported
	phase2Platforms := []struct {
		model    string
		flowCell string
	}{
		{"Element AVITI", "Dual"},
		{"PromethION 24", "R10.4.1"},
		{"PromethION 24", "R9.4.1"},
	}

	for _, platform := range phase2Platforms {
		timingProfile := timing.GetSequencerTiming(platform.model, platform.flowCell)
		if timingProfile.Model == "" {
			t.Errorf("Platform %s with flow cell %s not supported", platform.model, platform.flowCell)
		}
		if timingProfile.Model != platform.model {
			t.Errorf("Expected model %s, got %s", platform.model, timingProfile.Model)
		}
		if timingProfile.FlowCellType != platform.flowCell {
			t.Errorf("Expected flow cell %s, got %s", platform.flowCell, timingProfile.FlowCellType)
		}
	}
}
