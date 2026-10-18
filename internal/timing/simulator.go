// Package timing paces a simulated sequencing run (initialization, cluster
// generation, sequencing, base calling) and calls back as each output file
// is due, optionally compressed in time. The per-instrument run times and
// data rates in GetSequencerTiming are rough illustrative figures, not
// measurements or vendor specifications.
package timing

import (
	"context"
	"fmt"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"time"
)

// TimeCompression represents different time compression ratios
type TimeCompression int

const (
	RealTime        TimeCompression = 1
	Compressed10x   TimeCompression = 10
	Compressed100x  TimeCompression = 100
	Compressed1000x TimeCompression = 1000
)

// SequencerTiming holds timing information for different sequencer models
type SequencerTiming struct {
	Model              string
	FlowCellType       string
	TotalRunTimeHours  int
	ClusterGenMinutes  int
	SequencingMinutes  int
	BaseCallingMinutes int
	FilesPerHour       int
	BytesPerSecond     int64
}

// GetSequencerTiming returns timing specifications for different sequencers
func GetSequencerTiming(model, flowCellType string) SequencerTiming {
	timings := map[string]map[string]SequencerTiming{
		"NovaSeq 6000": {
			"S1": {
				Model:              "NovaSeq 6000",
				FlowCellType:       "S1",
				TotalRunTimeHours:  13,
				ClusterGenMinutes:  120,
				SequencingMinutes:  600,
				BaseCallingMinutes: 60,
				FilesPerHour:       2,
				BytesPerSecond:     50 * 1024 * 1024, // 50 MB/s
			},
			"S2": {
				Model:              "NovaSeq 6000",
				FlowCellType:       "S2",
				TotalRunTimeHours:  36,
				ClusterGenMinutes:  180,
				SequencingMinutes:  1980,
				BaseCallingMinutes: 240,
				FilesPerHour:       1,
				BytesPerSecond:     80 * 1024 * 1024, // 80 MB/s
			},
			"S4": {
				Model:              "NovaSeq 6000",
				FlowCellType:       "S4",
				TotalRunTimeHours:  44,
				ClusterGenMinutes:  240,
				SequencingMinutes:  2400,
				BaseCallingMinutes: 360,
				FilesPerHour:       1,
				BytesPerSecond:     120 * 1024 * 1024, // 120 MB/s
			},
		},
		"NextSeq 2000": {
			"P3": {
				Model:              "NextSeq 2000",
				FlowCellType:       "P3",
				TotalRunTimeHours:  18,
				ClusterGenMinutes:  60,
				SequencingMinutes:  900,
				BaseCallingMinutes: 120,
				FilesPerHour:       3,
				BytesPerSecond:     25 * 1024 * 1024, // 25 MB/s
			},
			"P4": {
				Model:              "NextSeq 2000",
				FlowCellType:       "P4",
				TotalRunTimeHours:  24,
				ClusterGenMinutes:  90,
				SequencingMinutes:  1200,
				BaseCallingMinutes: 150,
				FilesPerHour:       2,
				BytesPerSecond:     40 * 1024 * 1024, // 40 MB/s
			},
		},
		"MiSeq": {
			"v3": {
				Model:              "MiSeq",
				FlowCellType:       "v3",
				TotalRunTimeHours:  65,
				ClusterGenMinutes:  240,
				SequencingMinutes:  3600,
				BaseCallingMinutes: 60,
				FilesPerHour:       1,
				BytesPerSecond:     2 * 1024 * 1024, // 2 MB/s
			},
			"v2": {
				Model:              "MiSeq",
				FlowCellType:       "v2",
				TotalRunTimeHours:  27,
				ClusterGenMinutes:  120,
				SequencingMinutes:  1440,
				BaseCallingMinutes: 60,
				FilesPerHour:       2,
				BytesPerSecond:     5 * 1024 * 1024, // 5 MB/s
			},
		},
		"PacBio Revio": {
			"SMRT": {
				Model:              "PacBio Revio",
				FlowCellType:       "SMRT",
				TotalRunTimeHours:  24,
				ClusterGenMinutes:  60,                // Sample loading and preparation
				SequencingMinutes:  1200,              // 20 hours active sequencing
				BaseCallingMinutes: 60,                // HiFi consensus calling
				FilesPerHour:       4,                 // Higher throughput for long reads
				BytesPerSecond:     200 * 1024 * 1024, // 200 MB/s (higher due to long reads)
			},
		},
		"DNBSEQ-T20x2": {
			"DNB": {
				Model:              "DNBSEQ-T20x2",
				FlowCellType:       "DNB",
				TotalRunTimeHours:  48,
				ClusterGenMinutes:  240,               // DNB formation takes longer
				SequencingMinutes:  2400,              // 40 hours active sequencing (ultra-high throughput)
				BaseCallingMinutes: 120,               // Data processing
				FilesPerHour:       6,                 // Very high throughput
				BytesPerSecond:     500 * 1024 * 1024, // 500 MB/s (ultra-high throughput)
			},
		},
		"Element AVITI": {
			"Dual": {
				Model:              "Element AVITI",
				FlowCellType:       "Dual",
				TotalRunTimeHours:  38,
				ClusterGenMinutes:  60,                // Rolling circle amplification
				SequencingMinutes:  2040,              // 34 hours active sequencing (dual flow cell)
				BaseCallingMinutes: 120,               // Real-time basecalling with UltraQ
				FilesPerHour:       8,                 // High throughput dual flow cell system
				BytesPerSecond:     300 * 1024 * 1024, // 300 MB/s (high accuracy, high throughput)
			},
		},
		"PromethION 24": {
			"R10.4.1": {
				Model:              "PromethION 24",
				FlowCellType:       "R10.4.1",
				TotalRunTimeHours:  72,
				ClusterGenMinutes:  30,                // Minimal sample preparation for nanopore
				SequencingMinutes:  4290,              // 71.5 hours active sequencing (real-time)
				BaseCallingMinutes: 0,                 // Real-time GPU basecalling during sequencing
				FilesPerHour:       24,                // 24 flow cells, continuous output
				BytesPerSecond:     800 * 1024 * 1024, // 800 MB/s (24 flow cells × ultra-long reads)
			},
			"R9.4.1": {
				Model:              "PromethION 24",
				FlowCellType:       "R9.4.1",
				TotalRunTimeHours:  72,
				ClusterGenMinutes:  30,                // Minimal sample preparation for nanopore
				SequencingMinutes:  4290,              // 71.5 hours active sequencing (real-time)
				BaseCallingMinutes: 0,                 // Real-time GPU basecalling during sequencing
				FilesPerHour:       20,                // Slightly lower throughput with older chemistry
				BytesPerSecond:     600 * 1024 * 1024, // 600 MB/s (older chemistry, lower accuracy)
			},
		},
	}

	if modelTimings, exists := timings[model]; exists {
		if timing, exists := modelTimings[flowCellType]; exists {
			return timing
		}
	}

	// Default timing for unknown models
	return SequencerTiming{
		Model:              model,
		FlowCellType:       flowCellType,
		TotalRunTimeHours:  24,
		ClusterGenMinutes:  120,
		SequencingMinutes:  1200,
		BaseCallingMinutes: 120,
		FilesPerHour:       2,
		BytesPerSecond:     10 * 1024 * 1024,
	}
}

// RunSimulator paces a simulated run using the stage durations in Timing.
type RunSimulator struct {
	Config           types.SequencerConfig
	Timing           SequencerTiming
	Compression      TimeCompression
	FileCallback     func(filename string, batchNum int, readNum int)
	ProgressCallback func(stage string, percent float64, eta time.Duration)
	ErrorCallback    func(err error)
}

// NewRunSimulator creates a new run simulator
func NewRunSimulator(config types.SequencerConfig, compression TimeCompression) *RunSimulator {
	timing := GetSequencerTiming(config.SequencerModel, config.FlowCellType)

	return &RunSimulator{
		Config:      config,
		Timing:      timing,
		Compression: compression,
	}
}

// RunStage represents different stages of a sequencing run
type RunStage string

const (
	StageInitialization RunStage = "Initialization"
	StageClusterGen     RunStage = "Cluster Generation"
	StageSequencing     RunStage = "Sequencing"
	StageBaseCalling    RunStage = "Base Calling"
	StageComplete       RunStage = "Complete"
)

// RunProgress represents the current progress of a sequencing run
type RunProgress struct {
	Stage          RunStage
	Percent        float64
	ETA            time.Duration
	FilesGenerated int
	TotalFiles     int
	BytesGenerated int64
	TotalBytes     int64
	CurrentFile    string
}

// EstimateFileSize calculates expected file size for a given configuration
func (rs *RunSimulator) EstimateFileSize(readsPerFile int) int64 {
	// FASTQ file size estimation
	// Each read: 4 lines (header + sequence + "+" + quality)
	// Header: ~50 chars average
	// Sequence: read length
	// "+": 1 char
	// Quality: read length

	headerSize := int64(50)
	sequenceSize := int64(rs.Config.ReadLength)
	qualitySize := int64(rs.Config.ReadLength)
	lineBreaks := int64(4) // 4 newlines per read

	bytesPerRead := headerSize + sequenceSize + 1 + qualitySize + lineBreaks
	totalBytes := bytesPerRead * int64(readsPerFile)

	// Uncompressed size; gzip output is not estimated.
	return totalBytes
}

// EstimateTotalOutput calculates total expected output for the run
func (rs *RunSimulator) EstimateTotalOutput() (int, int64) {
	totalFiles := rs.Config.TotalBatches
	if rs.Config.PairedEnd {
		totalFiles *= 2 // R1 and R2 files
	}

	bytesPerFile := rs.EstimateFileSize(rs.Config.ReadsPerBatch)
	totalBytes := bytesPerFile * int64(totalFiles)

	return totalFiles, totalBytes
}

// Run steps through the run stages, waiting for each stage's duration
// divided by Compression, and calls FileCallback as files are due.
func (rs *RunSimulator) Run(ctx context.Context) error {
	totalFiles, totalBytes := rs.EstimateTotalOutput()

	// Calculate compressed timing
	realDuration := time.Duration(rs.Timing.TotalRunTimeHours) * time.Hour
	_ = realDuration / time.Duration(rs.Compression) // compressedDuration (not used currently)

	startTime := time.Now()

	// Stage 1: Initialization
	if err := rs.runStage(ctx, StageInitialization, 2*time.Minute/time.Duration(rs.Compression), 0, totalFiles); err != nil {
		return err
	}

	// Stage 2: Cluster Generation
	clusterTime := time.Duration(rs.Timing.ClusterGenMinutes) * time.Minute / time.Duration(rs.Compression)
	if err := rs.runStage(ctx, StageClusterGen, clusterTime, 0, totalFiles); err != nil {
		return err
	}

	// Stage 3: Sequencing (with file generation)
	sequencingTime := time.Duration(rs.Timing.SequencingMinutes) * time.Minute / time.Duration(rs.Compression)
	if err := rs.runSequencingStage(ctx, sequencingTime, totalFiles, totalBytes); err != nil {
		return err
	}

	// Stage 4: Base Calling (final processing)
	baseCallTime := time.Duration(rs.Timing.BaseCallingMinutes) * time.Minute / time.Duration(rs.Compression)
	if err := rs.runStage(ctx, StageBaseCalling, baseCallTime, totalFiles, totalFiles); err != nil {
		return err
	}

	// Stage 5: Complete
	if rs.ProgressCallback != nil {
		rs.ProgressCallback(string(StageComplete), 100.0, 0)
	}

	actualDuration := time.Since(startTime)
	fmt.Printf("Run completed in %v (compressed from %v at %dx speed)\n",
		actualDuration, realDuration, rs.Compression)

	return nil
}

// runStage executes a non-file-generating stage
func (rs *RunSimulator) runStage(ctx context.Context, stage RunStage, duration time.Duration, currentFiles, totalFiles int) error {
	if rs.ProgressCallback != nil {
		rs.ProgressCallback(string(stage), 0.0, duration)
	}

	// Update progress periodically during stage
	ticker := time.NewTicker(duration / 20) // 20 updates per stage
	defer ticker.Stop()

	startTime := time.Now()

	for {
		select {
		case <-ctx.Done():
			return ctx.Err()
		case <-ticker.C:
			elapsed := time.Since(startTime)
			percent := float64(elapsed) / float64(duration) * 100
			remaining := duration - elapsed

			if percent >= 100 {
				if rs.ProgressCallback != nil {
					rs.ProgressCallback(string(stage), 100.0, 0)
				}
				return nil
			}

			if rs.ProgressCallback != nil {
				rs.ProgressCallback(string(stage), percent, remaining)
			}
		}
	}
}

// runSequencingStage executes the sequencing stage with file generation
func (rs *RunSimulator) runSequencingStage(ctx context.Context, duration time.Duration, totalFiles int, totalBytes int64) error {
	if rs.ProgressCallback != nil {
		rs.ProgressCallback(string(StageSequencing), 0.0, duration)
	}

	if totalFiles <= 0 {
		return fmt.Errorf("runSequencingStage: totalFiles must be > 0 (got %d)", totalFiles)
	}

	// Calculate file generation timing
	fileInterval := duration / time.Duration(totalFiles)

	startTime := time.Now()
	filesGenerated := 0

	// Generate files at realistic intervals
	ticker := time.NewTicker(fileInterval)
	defer ticker.Stop()

	for batch := 1; batch <= rs.Config.TotalBatches; batch++ {
		select {
		case <-ctx.Done():
			return ctx.Err()
		case <-ticker.C:
			// Generate R1 file
			if rs.FileCallback != nil {
				filename := fmt.Sprintf("%s_S1_L%03d_R1_%03d.fastq",
					rs.Config.SampleName, rs.Config.Lane, batch)
				rs.FileCallback(filename, batch, 1)
				filesGenerated++
			}

			// Generate R2 file if paired-end
			if rs.Config.PairedEnd {
				select {
				case <-ctx.Done():
					return ctx.Err()
				case <-ticker.C:
				}
				if rs.FileCallback != nil {
					filename := fmt.Sprintf("%s_S1_L%03d_R2_%03d.fastq",
						rs.Config.SampleName, rs.Config.Lane, batch)
					rs.FileCallback(filename, batch, 2)
					filesGenerated++
				}
			}

			// Update progress
			elapsed := time.Since(startTime)
			percent := float64(filesGenerated) / float64(totalFiles) * 100
			remaining := duration - elapsed

			if rs.ProgressCallback != nil {
				rs.ProgressCallback(string(StageSequencing), percent, remaining)
			}

			if filesGenerated >= totalFiles {
				return nil
			}
		}
	}

	return nil
}

// GetCompressionDescription returns a human-readable description of time compression
func GetCompressionDescription(compression TimeCompression, realHours int) string {
	compressedMinutes := float64(realHours*60) / float64(compression)

	if compressedMinutes < 60 {
		return fmt.Sprintf("%.1f minutes (%.0fx compressed from %d hours)",
			compressedMinutes, float64(compression), realHours)
	}

	compressedHours := compressedMinutes / 60
	return fmt.Sprintf("%.1f hours (%.0fx compressed from %d hours)",
		compressedHours, float64(compression), realHours)
}
