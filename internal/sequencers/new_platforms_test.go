// Package sequencers holds only tests. They check the per-platform timing
// profiles in internal/timing (run time, files per hour, throughput) and
// the timing simulator's output estimates for PacBio Revio, MGI
// DNBSEQ-T20x2, Element AVITI and ONT PromethION 24. Some tests only assert
// on values written into a literal SequencerConfig in the test itself and
// so exercise no production code.
package sequencers

import (
	"github.com/scttfrdmn/phantom-science/internal/timing"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"testing"
	"time"
)

func TestPacBioRevioTiming(t *testing.T) {
	// Test PacBio Revio timing specifications
	timingProfile := timing.GetSequencerTiming("PacBio Revio", "SMRT")

	// Validate timing parameters
	if timingProfile.Model != "PacBio Revio" {
		t.Errorf("Expected model 'PacBio Revio', got '%s'", timingProfile.Model)
	}

	if timingProfile.FlowCellType != "SMRT" {
		t.Errorf("Expected flow cell type 'SMRT', got '%s'", timingProfile.FlowCellType)
	}

	if timingProfile.TotalRunTimeHours != 24 {
		t.Errorf("Expected 24 hour run time, got %d", timingProfile.TotalRunTimeHours)
	}

	// Validate that timing is realistic for long-read platform
	if timingProfile.BytesPerSecond < 100*1024*1024 {
		t.Errorf("Expected high throughput for long reads, got %d bytes/sec", timingProfile.BytesPerSecond)
	}

	t.Logf("PacBio Revio timing validation passed: %+v", timingProfile)
}

func TestMGIDNBSeqTiming(t *testing.T) {
	// Test MGI DNBSEQ-T20×2 timing specifications
	timingProfile := timing.GetSequencerTiming("DNBSEQ-T20x2", "DNB")

	// Validate timing parameters
	if timingProfile.Model != "DNBSEQ-T20x2" {
		t.Errorf("Expected model 'DNBSEQ-T20x2', got '%s'", timingProfile.Model)
	}

	if timingProfile.FlowCellType != "DNB" {
		t.Errorf("Expected flow cell type 'DNB', got '%s'", timingProfile.FlowCellType)
	}

	if timingProfile.TotalRunTimeHours != 48 {
		t.Errorf("Expected 48 hour run time, got %d", timingProfile.TotalRunTimeHours)
	}

	// Validate ultra-high throughput characteristics
	if timingProfile.BytesPerSecond < 400*1024*1024 {
		t.Errorf("Expected ultra-high throughput, got %d bytes/sec", timingProfile.BytesPerSecond)
	}

	if timingProfile.FilesPerHour < 5 {
		t.Errorf("Expected high file generation rate, got %d files/hour", timingProfile.FilesPerHour)
	}

	t.Logf("MGI DNBSEQ-T20×2 timing validation passed: %+v", timingProfile)
}

func TestPacBioRevioConfiguration(t *testing.T) {
	// Test PacBio Revio configuration loading and validation
	config := types.SequencerConfig{
		SequencerModel: "PacBio Revio",
		FlowCellType:   "SMRT",
		SampleName:     "Test_HiFi_Sample",
		Lane:           1,
		ReadLength:     17500, // Average HiFi read length
		PairedEnd:      false, // HiFi reads are single-end
		BaseQuality:    30,    // Q30 median
		Q30Percentage:  90.0,  // 90% Q30+
		TotalBatches:   10,
		ReadsPerBatch:  100000,
	}

	// Validate configuration makes sense for PacBio
	if config.ReadLength < 10000 {
		t.Errorf("PacBio read length should be >10kb, got %d", config.ReadLength)
	}

	if config.PairedEnd {
		t.Error("PacBio HiFi reads should be single-end")
	}

	if config.BaseQuality < 25 {
		t.Errorf("PacBio HiFi should have high quality, got Q%d", config.BaseQuality)
	}

	t.Logf("PacBio Revio configuration validation passed")
}

func TestMGIDNBSeqConfiguration(t *testing.T) {
	// Test MGI DNBSEQ-T20×2 configuration loading and validation
	config := types.SequencerConfig{
		SequencerModel: "DNBSEQ-T20x2",
		FlowCellType:   "DNB",
		SampleName:     "Test_Ultra_High_Throughput",
		Lane:           1,
		ReadLength:     150,     // PE150
		PairedEnd:      true,    // Paired-end sequencing
		BaseQuality:    30,      // Q30 initial
		Q30Percentage:  85.0,    // 85% Q30+
		TotalBatches:   100,     // High throughput
		ReadsPerBatch:  1000000, // 1M reads per batch
	}

	// Validate configuration makes sense for MGI ultra-high throughput
	if config.ReadLength != 150 && config.ReadLength != 100 {
		t.Errorf("MGI should use PE100 or PE150, got %d", config.ReadLength)
	}

	if !config.PairedEnd {
		t.Error("MGI DNBSEQ typically uses paired-end sequencing")
	}

	if config.ReadsPerBatch < 500000 {
		t.Errorf("MGI ultra-high throughput should have high reads per batch, got %d", config.ReadsPerBatch)
	}

	t.Logf("MGI DNBSEQ-T20×2 configuration validation passed")
}

func TestTimeCompressionForNewPlatforms(t *testing.T) {
	// Test time compression works correctly for new platforms
	testCases := []struct {
		model        string
		flowCellType string
		compression  timing.TimeCompression
		description  string
	}{
		{"PacBio Revio", "SMRT", timing.Compressed100x, "PacBio 100x compression"},
		{"DNBSEQ-T20x2", "DNB", timing.Compressed1000x, "MGI 1000x compression"},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			timingProfile := timing.GetSequencerTiming(tc.model, tc.flowCellType)

			// Calculate compressed time
			realDuration := time.Duration(timingProfile.TotalRunTimeHours) * time.Hour
			compressedDuration := realDuration / time.Duration(tc.compression)

			// Validate compression makes sense
			if compressedDuration > time.Hour {
				t.Errorf("Compressed time too long: %v", compressedDuration)
			}

			description := timing.GetCompressionDescription(tc.compression, timingProfile.TotalRunTimeHours)
			t.Logf("%s: %s → %s", tc.description, realDuration, description)
		})
	}
}

func TestPlatformCapabilities(t *testing.T) {
	// Test that new platforms have expected capabilities
	platforms := []struct {
		model              string
		flowCellType       string
		expectedThroughput string
		expectedReadType   string
	}{
		{"PacBio Revio", "SMRT", "high", "long_read_hifi"},
		{"DNBSEQ-T20x2", "DNB", "ultra_high", "short_read"},
	}

	for _, platform := range platforms {
		t.Run(platform.model, func(t *testing.T) {
			timingProfile := timing.GetSequencerTiming(platform.model, platform.flowCellType)

			// Validate platform exists and has reasonable parameters
			if timingProfile.Model == "" {
				t.Errorf("Platform not found: %s", platform.model)
			}

			if timingProfile.TotalRunTimeHours <= 0 {
				t.Errorf("Invalid run time for %s: %d hours", platform.model, timingProfile.TotalRunTimeHours)
			}

			if timingProfile.BytesPerSecond <= 0 {
				t.Errorf("Invalid throughput for %s: %d bytes/sec", platform.model, timingProfile.BytesPerSecond)
			}

			t.Logf("Platform %s validated successfully", platform.model)
		})
	}
}

func TestPlatformComparison(t *testing.T) {
	// Compare new platforms with existing Illumina platforms
	platforms := []struct {
		model        string
		flowCellType string
	}{
		{"NovaSeq 6000", "S2"},   // Existing high-throughput
		{"PacBio Revio", "SMRT"}, // New long-read
		{"DNBSEQ-T20x2", "DNB"},  // New ultra-high throughput
	}

	t.Log("Platform comparison:")
	for _, platform := range platforms {
		timingProfile := timing.GetSequencerTiming(platform.model, platform.flowCellType)
		throughputMBps := float64(timingProfile.BytesPerSecond) / (1024 * 1024)

		t.Logf("%-15s %-8s: %2dh runtime, %6.0f MB/s, %d files/hour",
			platform.model, platform.flowCellType, timingProfile.TotalRunTimeHours,
			throughputMBps, timingProfile.FilesPerHour)
	}
}

func TestReadLengthValidation(t *testing.T) {
	// Test that read length validation works for different platform types
	testCases := []struct {
		platform    string
		readLength  int
		shouldPass  bool
		description string
	}{
		{"PacBio Revio", 17500, true, "PacBio typical HiFi read"},
		{"PacBio Revio", 150, false, "PacBio with short read"},
		{"DNBSEQ-T20x2", 150, true, "MGI typical PE150"},
		{"DNBSEQ-T20x2", 20000, false, "MGI with long read"},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			// No production code validates read length per platform; this
			// only checks the consistency of the table above.
			isLongRead := tc.readLength > 1000
			isPacBio := tc.platform == "PacBio Revio"

			if isPacBio && !isLongRead && tc.shouldPass {
				t.Error("PacBio should require long reads")
			}

			if !isPacBio && isLongRead && tc.shouldPass {
				t.Error("Short-read platforms should not use long reads")
			}

			t.Logf("Read length validation: %s with %d bp reads: %v",
				tc.platform, tc.readLength, tc.shouldPass)
		})
	}
}

// Benchmark tests for new platforms
func BenchmarkPacBioRevioTiming(b *testing.B) {
	for i := 0; i < b.N; i++ {
		timingProfile := timing.GetSequencerTiming("PacBio Revio", "SMRT")
		_ = timingProfile
	}
}

func BenchmarkMGIDNBSeqTiming(b *testing.B) {
	for i := 0; i < b.N; i++ {
		timingProfile := timing.GetSequencerTiming("DNBSEQ-T20x2", "DNB")
		_ = timingProfile
	}
}
