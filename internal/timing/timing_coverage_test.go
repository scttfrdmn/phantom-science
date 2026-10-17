//go:build !release

package timing_test

import (
	"context"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/timing"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ---------------------------------------------------------------------------
// TimeCompression constants
// ---------------------------------------------------------------------------

func TestTimeCompression_Constants(t *testing.T) {
	assert.Equal(t, timing.TimeCompression(1), timing.RealTime)
	assert.Equal(t, timing.TimeCompression(10), timing.Compressed10x)
	assert.Equal(t, timing.TimeCompression(100), timing.Compressed100x)
	assert.Equal(t, timing.TimeCompression(1000), timing.Compressed1000x)

	assert.True(t, timing.RealTime < timing.Compressed10x)
	assert.True(t, timing.Compressed10x < timing.Compressed100x)
	assert.True(t, timing.Compressed100x < timing.Compressed1000x)
}

// ---------------------------------------------------------------------------
// SequencerTiming struct
// ---------------------------------------------------------------------------

func TestSequencerTiming_StructFields(t *testing.T) {
	st := timing.SequencerTiming{
		Model:              "TestSeq",
		FlowCellType:       "TestFlow",
		TotalRunTimeHours:  48,
		ClusterGenMinutes:  180,
		SequencingMinutes:  2700,
		BaseCallingMinutes: 120,
		FilesPerHour:       3,
		BytesPerSecond:     100 * 1024 * 1024,
	}

	assert.Equal(t, "TestSeq", st.Model)
	assert.Equal(t, "TestFlow", st.FlowCellType)
	assert.Equal(t, 48, st.TotalRunTimeHours)
	assert.Equal(t, 180, st.ClusterGenMinutes)
	assert.Equal(t, 2700, st.SequencingMinutes)
	assert.Equal(t, 120, st.BaseCallingMinutes)
	assert.Equal(t, 3, st.FilesPerHour)
	assert.Equal(t, int64(100*1024*1024), st.BytesPerSecond)
}

// ---------------------------------------------------------------------------
// GetSequencerTiming — all known platforms
// ---------------------------------------------------------------------------

func TestGetSequencerTiming_NovaSeq6000_S1(t *testing.T) {
	st := timing.GetSequencerTiming("NovaSeq 6000", "S1")
	assert.Equal(t, "NovaSeq 6000", st.Model)
	assert.Equal(t, "S1", st.FlowCellType)
	assert.Equal(t, 13, st.TotalRunTimeHours)
	assert.True(t, st.BytesPerSecond > 0)
}

func TestGetSequencerTiming_NovaSeq6000_S2(t *testing.T) {
	st := timing.GetSequencerTiming("NovaSeq 6000", "S2")
	assert.Equal(t, 36, st.TotalRunTimeHours)
	assert.True(t, st.SequencingMinutes > st.ClusterGenMinutes)
}

func TestGetSequencerTiming_NovaSeq6000_S4(t *testing.T) {
	st := timing.GetSequencerTiming("NovaSeq 6000", "S4")
	assert.Equal(t, 44, st.TotalRunTimeHours)
	assert.Equal(t, int64(120*1024*1024), st.BytesPerSecond)
}

func TestGetSequencerTiming_NextSeq2000_P3(t *testing.T) {
	st := timing.GetSequencerTiming("NextSeq 2000", "P3")
	assert.Equal(t, "NextSeq 2000", st.Model)
	assert.Equal(t, "P3", st.FlowCellType)
	assert.Equal(t, 18, st.TotalRunTimeHours)
}

func TestGetSequencerTiming_NextSeq2000_P4(t *testing.T) {
	st := timing.GetSequencerTiming("NextSeq 2000", "P4")
	assert.Equal(t, 24, st.TotalRunTimeHours)
}

func TestGetSequencerTiming_MiSeq_v3(t *testing.T) {
	st := timing.GetSequencerTiming("MiSeq", "v3")
	assert.Equal(t, "MiSeq", st.Model)
	assert.Equal(t, "v3", st.FlowCellType)
	assert.Equal(t, 65, st.TotalRunTimeHours)
}

func TestGetSequencerTiming_MiSeq_v2(t *testing.T) {
	st := timing.GetSequencerTiming("MiSeq", "v2")
	assert.Equal(t, 27, st.TotalRunTimeHours)
}

func TestGetSequencerTiming_PacBioRevio_SMRT(t *testing.T) {
	st := timing.GetSequencerTiming("PacBio Revio", "SMRT")
	assert.Equal(t, "PacBio Revio", st.Model)
	assert.Equal(t, "SMRT", st.FlowCellType)
	assert.Equal(t, 24, st.TotalRunTimeHours)
	assert.True(t, st.BytesPerSecond > 0)
}

func TestGetSequencerTiming_DNBSEQ_T20x2(t *testing.T) {
	st := timing.GetSequencerTiming("DNBSEQ-T20x2", "DNB")
	assert.Equal(t, "DNBSEQ-T20x2", st.Model)
	assert.Equal(t, 48, st.TotalRunTimeHours)
	assert.Equal(t, int64(500*1024*1024), st.BytesPerSecond)
}

func TestGetSequencerTiming_ElementAVITI(t *testing.T) {
	st := timing.GetSequencerTiming("Element AVITI", "Dual")
	assert.Equal(t, "Element AVITI", st.Model)
	assert.Equal(t, 38, st.TotalRunTimeHours)
}

func TestGetSequencerTiming_PromethION24_R10(t *testing.T) {
	st := timing.GetSequencerTiming("PromethION 24", "R10.4.1")
	assert.Equal(t, "PromethION 24", st.Model)
	assert.Equal(t, 72, st.TotalRunTimeHours)
	assert.Equal(t, 0, st.BaseCallingMinutes) // real-time basecalling
}

func TestGetSequencerTiming_PromethION24_R9(t *testing.T) {
	st := timing.GetSequencerTiming("PromethION 24", "R9.4.1")
	assert.Equal(t, 72, st.TotalRunTimeHours)
}

func TestGetSequencerTiming_UnknownModel_DefaultFallback(t *testing.T) {
	st := timing.GetSequencerTiming("UnknownSeq", "UnknownFlow")
	assert.Equal(t, "UnknownSeq", st.Model)
	assert.Equal(t, "UnknownFlow", st.FlowCellType)
	assert.Equal(t, 24, st.TotalRunTimeHours)
	assert.Equal(t, 120, st.ClusterGenMinutes)
	assert.Equal(t, 1200, st.SequencingMinutes)
	assert.Equal(t, 120, st.BaseCallingMinutes)
	assert.Equal(t, 2, st.FilesPerHour)
	assert.Equal(t, int64(10*1024*1024), st.BytesPerSecond)
}

func TestGetSequencerTiming_KnownModelUnknownFlowCell(t *testing.T) {
	// Known model, but unknown flow cell type — must fall back to defaults.
	st := timing.GetSequencerTiming("NovaSeq 6000", "X99")
	assert.Equal(t, "NovaSeq 6000", st.Model)
	assert.Equal(t, "X99", st.FlowCellType)
	assert.Equal(t, 24, st.TotalRunTimeHours) // default
}

func TestGetSequencerTiming_AllPositiveFields(t *testing.T) {
	platforms := []struct{ model, flow string }{
		{"NovaSeq 6000", "S1"},
		{"NovaSeq 6000", "S2"},
		{"NovaSeq 6000", "S4"},
		{"NextSeq 2000", "P3"},
		{"NextSeq 2000", "P4"},
		{"MiSeq", "v3"},
		{"MiSeq", "v2"},
		{"PacBio Revio", "SMRT"},
		{"DNBSEQ-T20x2", "DNB"},
		{"Element AVITI", "Dual"},
		{"PromethION 24", "R10.4.1"},
		{"PromethION 24", "R9.4.1"},
	}

	for _, p := range platforms {
		st := timing.GetSequencerTiming(p.model, p.flow)
		assert.True(t, st.TotalRunTimeHours > 0, "%s/%s: TotalRunTimeHours must be positive", p.model, p.flow)
		assert.True(t, st.SequencingMinutes > 0, "%s/%s: SequencingMinutes must be positive", p.model, p.flow)
		assert.True(t, st.BytesPerSecond > 0, "%s/%s: BytesPerSecond must be positive", p.model, p.flow)
		assert.True(t, st.FilesPerHour > 0, "%s/%s: FilesPerHour must be positive", p.model, p.flow)
		assert.True(t, st.ClusterGenMinutes >= 0, "%s/%s: ClusterGenMinutes must be non-negative", p.model, p.flow)
		assert.True(t, st.BaseCallingMinutes >= 0, "%s/%s: BaseCallingMinutes must be non-negative", p.model, p.flow)
	}
}

// ---------------------------------------------------------------------------
// NewRunSimulator
// ---------------------------------------------------------------------------

func TestNewRunSimulator_BasicCreation(t *testing.T) {
	cfg := basicConfig()
	sim := timing.NewRunSimulator(cfg, timing.Compressed1000x)
	require.NotNil(t, sim)
}

func TestNewRunSimulator_Fields(t *testing.T) {
	cfg := basicConfig()
	sim := timing.NewRunSimulator(cfg, timing.Compressed100x)

	assert.Equal(t, timing.Compressed100x, sim.Compression)
	assert.Equal(t, cfg.SequencerModel, sim.Config.SequencerModel)
	assert.Equal(t, cfg.ReadLength, sim.Config.ReadLength)
	assert.Equal(t, "NovaSeq 6000", sim.Timing.Model)
	assert.True(t, sim.Timing.TotalRunTimeHours > 0)
}

func TestNewRunSimulator_AllCompressions(t *testing.T) {
	cfg := basicConfig()
	compressions := []timing.TimeCompression{
		timing.RealTime,
		timing.Compressed10x,
		timing.Compressed100x,
		timing.Compressed1000x,
	}

	for _, c := range compressions {
		sim := timing.NewRunSimulator(cfg, c)
		require.NotNil(t, sim, "simulator must not be nil for compression %d", c)
		assert.Equal(t, c, sim.Compression)
	}
}

func TestNewRunSimulator_UnknownPlatform(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "FuturePlatform",
		FlowCellType:   "NextGen",
		ReadLength:     200,
		ReadsPerBatch:  500,
		TotalBatches:   3,
	}
	sim := timing.NewRunSimulator(cfg, timing.Compressed1000x)
	require.NotNil(t, sim)
	// Should have used fallback timing.
	assert.Equal(t, 24, sim.Timing.TotalRunTimeHours)
}

// ---------------------------------------------------------------------------
// EstimateFileSize
// ---------------------------------------------------------------------------

func TestEstimateFileSize_ShortReads(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "MiSeq",
		ReadLength:     50,
		ReadsPerBatch:  1000,
		TotalBatches:   1,
	}
	sim := timing.NewRunSimulator(cfg, timing.RealTime)
	size := sim.EstimateFileSize(100)
	assert.True(t, size > 0, "file size must be positive")
}

func TestEstimateFileSize_StandardReads(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		TotalBatches:   1,
	}
	sim := timing.NewRunSimulator(cfg, timing.RealTime)
	size1000 := sim.EstimateFileSize(1000)
	size100 := sim.EstimateFileSize(100)

	assert.True(t, size1000 > 0)
	assert.True(t, size100 > 0)
	// More reads = larger file.
	assert.True(t, size1000 > size100, "1000 reads must produce a larger file than 100 reads")
}

func TestEstimateFileSize_LongReads(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "PromethION 24",
		FlowCellType:   "R10.4.1",
		ReadLength:     300,
		ReadsPerBatch:  500,
		TotalBatches:   1,
	}
	sim := timing.NewRunSimulator(cfg, timing.RealTime)
	size := sim.EstimateFileSize(500)
	assert.True(t, size > 0)
}

func TestEstimateFileSize_ScalesWithReadLength(t *testing.T) {
	makeSimWithReadLen := func(readLen int) *timing.RunSimulator {
		cfg := types.SequencerConfig{
			SequencerModel: "NovaSeq 6000",
			ReadLength:     readLen,
			ReadsPerBatch:  1000,
			TotalBatches:   1,
		}
		return timing.NewRunSimulator(cfg, timing.RealTime)
	}

	size50 := makeSimWithReadLen(50).EstimateFileSize(1000)
	size150 := makeSimWithReadLen(150).EstimateFileSize(1000)
	size300 := makeSimWithReadLen(300).EstimateFileSize(1000)

	assert.True(t, size150 > size50, "longer reads must produce larger files (150 vs 50)")
	assert.True(t, size300 > size150, "longer reads must produce larger files (300 vs 150)")
}

// ---------------------------------------------------------------------------
// EstimateTotalOutput
// ---------------------------------------------------------------------------

func TestEstimateTotalOutput_SingleEnd(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		TotalBatches:   5,
		PairedEnd:      false,
	}
	sim := timing.NewRunSimulator(cfg, timing.RealTime)
	files, bytes := sim.EstimateTotalOutput()

	assert.Equal(t, 5, files)
	assert.True(t, bytes > 0)
}

func TestEstimateTotalOutput_PairedEnd(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		TotalBatches:   5,
		PairedEnd:      true,
	}
	sim := timing.NewRunSimulator(cfg, timing.RealTime)
	files, bytes := sim.EstimateTotalOutput()

	assert.Equal(t, 10, files) // 5 batches * 2 (R1 + R2)
	assert.True(t, bytes > 0)
}

func TestEstimateTotalOutput_PairedVsSingle_TwiceFiles(t *testing.T) {
	base := types.SequencerConfig{
		SequencerModel: "NextSeq 2000",
		FlowCellType:   "P3",
		ReadLength:     100,
		ReadsPerBatch:  500,
		TotalBatches:   4,
	}

	singleCfg := base
	singleCfg.PairedEnd = false
	singleSim := timing.NewRunSimulator(singleCfg, timing.RealTime)
	singleFiles, singleBytes := singleSim.EstimateTotalOutput()

	pairedCfg := base
	pairedCfg.PairedEnd = true
	pairedSim := timing.NewRunSimulator(pairedCfg, timing.RealTime)
	pairedFiles, pairedBytes := pairedSim.EstimateTotalOutput()

	assert.Equal(t, singleFiles*2, pairedFiles, "paired must produce 2x files")
	assert.Equal(t, singleBytes*2, pairedBytes, "paired must produce 2x bytes")
}

func TestEstimateTotalOutput_OneBatch(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		ReadsPerBatch:  500,
		TotalBatches:   1,
		PairedEnd:      false,
	}
	sim := timing.NewRunSimulator(cfg, timing.RealTime)
	files, bytes := sim.EstimateTotalOutput()

	assert.Equal(t, 1, files)
	assert.True(t, bytes > 0)
}

// ---------------------------------------------------------------------------
// RunStage constants
// ---------------------------------------------------------------------------

func TestRunStage_Constants(t *testing.T) {
	assert.Equal(t, timing.RunStage("Initialization"), timing.StageInitialization)
	assert.Equal(t, timing.RunStage("Cluster Generation"), timing.StageClusterGen)
	assert.Equal(t, timing.RunStage("Sequencing"), timing.StageSequencing)
	assert.Equal(t, timing.RunStage("Base Calling"), timing.StageBaseCalling)
	assert.Equal(t, timing.RunStage("Complete"), timing.StageComplete)
}

func TestRunStage_AsString(t *testing.T) {
	assert.Equal(t, "Initialization", string(timing.StageInitialization))
	assert.Equal(t, "Cluster Generation", string(timing.StageClusterGen))
	assert.Equal(t, "Sequencing", string(timing.StageSequencing))
	assert.Equal(t, "Base Calling", string(timing.StageBaseCalling))
	assert.Equal(t, "Complete", string(timing.StageComplete))
}

// ---------------------------------------------------------------------------
// RunProgress struct
// ---------------------------------------------------------------------------

func TestRunProgress_StructFields(t *testing.T) {
	rp := timing.RunProgress{
		Stage:          timing.StageSequencing,
		Percent:        42.5,
		ETA:            15 * time.Minute,
		FilesGenerated: 8,
		TotalFiles:     20,
		BytesGenerated: 1024 * 1024 * 200,
		TotalBytes:     1024 * 1024 * 500,
		CurrentFile:    "sample_S1_L001_R1_003.fastq",
	}

	assert.Equal(t, timing.StageSequencing, rp.Stage)
	assert.Equal(t, 42.5, rp.Percent)
	assert.Equal(t, 15*time.Minute, rp.ETA)
	assert.Equal(t, 8, rp.FilesGenerated)
	assert.Equal(t, 20, rp.TotalFiles)
	assert.Equal(t, "sample_S1_L001_R1_003.fastq", rp.CurrentFile)
}

// ---------------------------------------------------------------------------
// GetCompressionDescription
// ---------------------------------------------------------------------------

func TestGetCompressionDescription_1000x(t *testing.T) {
	desc := timing.GetCompressionDescription(timing.Compressed1000x, 24)
	assert.NotEmpty(t, desc)
	// 24 hours / 1000 = 1.44 minutes — should mention "minutes".
	assert.True(t, strings.Contains(desc, "minute"), "1000x compression of 24h should be expressed in minutes, got: %s", desc)
}

func TestGetCompressionDescription_100x(t *testing.T) {
	desc := timing.GetCompressionDescription(timing.Compressed100x, 24)
	assert.NotEmpty(t, desc)
	// 24 hours / 100 = 14.4 minutes — still under 60 minutes.
	assert.True(t, strings.Contains(desc, "minute"), "100x compression of 24h should be expressed in minutes, got: %s", desc)
}

func TestGetCompressionDescription_10x(t *testing.T) {
	desc := timing.GetCompressionDescription(timing.Compressed10x, 24)
	assert.NotEmpty(t, desc)
	// 24 hours / 10 = 2.4 hours — should mention "hours".
	assert.True(t, strings.Contains(desc, "hour"), "10x compression of 24h should be expressed in hours, got: %s", desc)
}

func TestGetCompressionDescription_RealTime(t *testing.T) {
	desc := timing.GetCompressionDescription(timing.RealTime, 24)
	assert.NotEmpty(t, desc)
	// 24 hours / 1 = 24 hours — should mention "hours".
	assert.True(t, strings.Contains(desc, "hour"), "realtime for 24h should be expressed in hours, got: %s", desc)
}

func TestGetCompressionDescription_ContainsRealHours(t *testing.T) {
	for _, hours := range []int{13, 24, 36, 48, 65, 72} {
		desc := timing.GetCompressionDescription(timing.Compressed100x, hours)
		assert.True(t, strings.Contains(desc, "from"), "description should contain 'from', got: %s", desc)
		assert.NotEmpty(t, desc)
	}
}

func TestGetCompressionDescription_AllPlatformHours(t *testing.T) {
	platforms := []struct {
		model string
		flow  string
	}{
		{"NovaSeq 6000", "S1"},
		{"MiSeq", "v3"},
		{"PromethION 24", "R10.4.1"},
	}

	for _, p := range platforms {
		st := timing.GetSequencerTiming(p.model, p.flow)
		desc := timing.GetCompressionDescription(timing.Compressed100x, st.TotalRunTimeHours)
		assert.NotEmpty(t, desc, "compression description for %s/%s must not be empty", p.model, p.flow)
	}
}

// ---------------------------------------------------------------------------
// Callbacks are wired correctly (set but not necessarily triggered in fast tests)
// ---------------------------------------------------------------------------

func TestRunSimulator_CallbacksCanBeSet(t *testing.T) {
	sim := timing.NewRunSimulator(basicConfig(), timing.Compressed1000x)

	var fileEvents []string
	var progressEvents []string
	var errorEvents []error

	sim.FileCallback = func(filename string, batch, read int) {
		fileEvents = append(fileEvents, filename)
	}
	sim.ProgressCallback = func(stage string, pct float64, eta time.Duration) {
		progressEvents = append(progressEvents, stage)
	}
	sim.ErrorCallback = func(err error) {
		errorEvents = append(errorEvents, err)
	}

	// Run with immediate cancellation.
	ctx, cancel := context.WithCancel(context.Background())
	cancel()

	err := sim.Run(ctx)
	assert.Error(t, err, "immediately cancelled context must return an error")
	assert.Equal(t, context.Canceled, err)
}

func TestRunSimulator_NilCallbacks_DoNotPanic(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		FlowCellType:   "S1",
		ReadLength:     150,
		ReadsPerBatch:  100,
		TotalBatches:   1,
		PairedEnd:      false,
		SampleName:     "NilTest",
		Lane:           1,
	}
	sim := timing.NewRunSimulator(cfg, timing.Compressed1000x)
	// All callbacks remain nil.

	ctx, cancel := context.WithTimeout(context.Background(), 100*time.Millisecond)
	defer cancel()

	// Should not panic even with nil callbacks.
	err := sim.Run(ctx)
	// May complete or time out; either is fine.
	_ = err
}

// ---------------------------------------------------------------------------
// Context cancellation
// ---------------------------------------------------------------------------

func TestRunSimulator_ImmediateCancel(t *testing.T) {
	sim := timing.NewRunSimulator(basicConfig(), timing.Compressed1000x)

	ctx, cancel := context.WithCancel(context.Background())
	cancel() // cancel before Run

	err := sim.Run(ctx)
	assert.Equal(t, context.Canceled, err)
}

func TestRunSimulator_TimeoutCancel(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		FlowCellType:   "S4",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		TotalBatches:   100,
		PairedEnd:      true,
		SampleName:     "TimeoutTest",
		Lane:           1,
	}
	sim := timing.NewRunSimulator(cfg, timing.Compressed10x)

	ctx, cancel := context.WithTimeout(context.Background(), 50*time.Millisecond)
	defer cancel()

	err := sim.Run(ctx)
	assert.Error(t, err, "a very short timeout must cause an error")
}

// ---------------------------------------------------------------------------
// File naming verification via FileCallback
// ---------------------------------------------------------------------------

func TestRunSimulator_FileNaming_SingleEnd(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v2",
		ReadLength:     150,
		ReadsPerBatch:  200,
		TotalBatches:   2,
		PairedEnd:      false,
		SampleName:     "MySample",
		Lane:           2,
	}
	sim := timing.NewRunSimulator(cfg, timing.Compressed1000x)

	var files []string
	sim.FileCallback = func(name string, batch, read int) {
		files = append(files, name)
	}

	ctx, cancel := context.WithTimeout(context.Background(), 500*time.Millisecond)
	defer cancel()

	_ = sim.Run(ctx)

	for _, f := range files {
		assert.True(t, strings.Contains(f, "MySample"), "filename must contain sample name, got: %s", f)
		assert.True(t, strings.HasSuffix(f, ".fastq"), "filename must end in .fastq, got: %s", f)
		// Should contain "R1" for single-end.
		assert.True(t, strings.Contains(f, "R1"), "single-end filename should contain R1, got: %s", f)
	}
}

func TestRunSimulator_FileNaming_PairedEnd_ContainsR1R2(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "NextSeq 2000",
		FlowCellType:   "P3",
		ReadLength:     100,
		ReadsPerBatch:  100,
		TotalBatches:   2,
		PairedEnd:      true,
		SampleName:     "PairedSample",
		Lane:           1,
	}
	sim := timing.NewRunSimulator(cfg, timing.Compressed1000x)

	var files []string
	sim.FileCallback = func(name string, batch, read int) {
		files = append(files, name)
	}

	ctx, cancel := context.WithTimeout(context.Background(), 1*time.Second)
	defer cancel()

	_ = sim.Run(ctx)

	// If any files were generated, check they are well-formed.
	for _, f := range files {
		assert.True(t, strings.Contains(f, "PairedSample"), "filename must contain sample name, got: %s", f)
		assert.True(t, strings.HasSuffix(f, ".fastq"), "filename must end in .fastq, got: %s", f)
	}
}

// ---------------------------------------------------------------------------
// Progress callback receives expected stage names
// ---------------------------------------------------------------------------

func TestRunSimulator_ProgressCallback_StageNames(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		FlowCellType:   "S1",
		ReadLength:     100,
		ReadsPerBatch:  100,
		TotalBatches:   1,
		PairedEnd:      false,
		SampleName:     "ProgressTest",
		Lane:           1,
	}
	sim := timing.NewRunSimulator(cfg, timing.Compressed1000x)

	stagesObserved := make(map[string]bool)
	sim.ProgressCallback = func(stage string, pct float64, eta time.Duration) {
		stagesObserved[stage] = true
	}

	ctx, cancel := context.WithTimeout(context.Background(), 2*time.Second)
	defer cancel()

	err := sim.Run(ctx)
	// Either completes or times out; we care only that callbacks fired.
	_ = err

	assert.True(t, len(stagesObserved) > 0, "at least one progress stage must have been observed")
}

// ---------------------------------------------------------------------------
// Config fields exposed via RunSimulator
// ---------------------------------------------------------------------------

func TestRunSimulator_ConfigPreservation(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "PacBio Revio",
		FlowCellType:   "SMRT",
		ReadLength:     15000,
		ReadsPerBatch:  50,
		TotalBatches:   3,
		PairedEnd:      false,
		SampleName:     "LongRead",
		Lane:           1,
	}
	sim := timing.NewRunSimulator(cfg, timing.Compressed100x)

	assert.Equal(t, "PacBio Revio", sim.Config.SequencerModel)
	assert.Equal(t, "SMRT", sim.Config.FlowCellType)
	assert.Equal(t, 15000, sim.Config.ReadLength)
	assert.Equal(t, 50, sim.Config.ReadsPerBatch)
	assert.Equal(t, 3, sim.Config.TotalBatches)
	assert.Equal(t, "LongRead", sim.Config.SampleName)
	assert.Equal(t, timing.Compressed100x, sim.Compression)
	assert.Equal(t, "PacBio Revio", sim.Timing.Model)
}

// ---------------------------------------------------------------------------
// EstimateFileSize vs EstimateTotalOutput consistency
// ---------------------------------------------------------------------------

func TestEstimateTotalOutput_ConsistentWithFileSize(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		TotalBatches:   5,
		PairedEnd:      false,
	}
	sim := timing.NewRunSimulator(cfg, timing.RealTime)

	files, totalBytes := sim.EstimateTotalOutput()
	perFileBytes := sim.EstimateFileSize(cfg.ReadsPerBatch)

	assert.Equal(t, cfg.TotalBatches, files)
	assert.Equal(t, perFileBytes*int64(files), totalBytes)
}

// ---------------------------------------------------------------------------
// Deeper simulation coverage — force tickers through small timescales
// ---------------------------------------------------------------------------

// TestRunSimulator_CompleteFullRun runs a complete simulation at maximum
// compression with a tiny dataset so the entire run completes within the test
// timeout while exercising the runStage loop completion branch.
func TestRunSimulator_CompleteFullRun(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v2",
		ReadLength:     50,
		ReadsPerBatch:  10,
		TotalBatches:   1,
		PairedEnd:      false,
		SampleName:     "FullRun",
		Lane:           1,
	}
	// MiSeq v2: 27h total. At 1000x that is 27*60/1000 = 1.62 minutes = ~97s.
	// That is too slow for a unit test.
	// Use a custom very-high compression to get sub-second execution.
	// We override the timing struct after construction.
	sim := timing.NewRunSimulator(cfg, timing.Compressed1000x)

	// Override timing so each stage is extremely short (~milliseconds).
	sim.Timing.TotalRunTimeHours = 1
	sim.Timing.ClusterGenMinutes = 1
	sim.Timing.SequencingMinutes = 1
	sim.Timing.BaseCallingMinutes = 1

	var completedStages []string
	sim.ProgressCallback = func(stage string, pct float64, eta time.Duration) {
		if pct >= 100.0 {
			completedStages = append(completedStages, stage)
		}
	}

	var generatedFiles []string
	sim.FileCallback = func(name string, batch, read int) {
		generatedFiles = append(generatedFiles, name)
	}

	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()

	err := sim.Run(ctx)
	// May complete normally or time out; either is acceptable.
	if err != nil {
		assert.True(t, err == context.DeadlineExceeded || err == context.Canceled,
			"unexpected error: %v", err)
	}
}

// TestRunSimulator_PairedEndFullRun exercises the paired-end sequencing branch
// inside runSequencingStage with a minimal dataset.
func TestRunSimulator_PairedEndFullRun(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v2",
		ReadLength:     50,
		ReadsPerBatch:  10,
		TotalBatches:   1,
		PairedEnd:      true,
		SampleName:     "PairedFullRun",
		Lane:           1,
	}
	sim := timing.NewRunSimulator(cfg, timing.Compressed1000x)
	sim.Timing.TotalRunTimeHours = 1
	sim.Timing.ClusterGenMinutes = 1
	sim.Timing.SequencingMinutes = 1
	sim.Timing.BaseCallingMinutes = 1

	var files []string
	sim.FileCallback = func(name string, batch, read int) {
		files = append(files, name)
	}

	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()

	err := sim.Run(ctx)
	if err != nil {
		assert.True(t, err == context.DeadlineExceeded || err == context.Canceled,
			"unexpected error: %v", err)
	}

	// If the sequencing stage ran, we expect R1 and R2 files.
	for _, f := range files {
		assert.True(t, strings.Contains(f, "R1") || strings.Contains(f, "R2"),
			"paired-end file must contain R1 or R2, got: %s", f)
	}
}

// TestRunSimulator_PromethION_TimingValues verifies that PromethION 24 timing
// is loaded correctly, including the zero BaseCallingMinutes value that
// reflects real-time GPU basecalling during sequencing.
func TestRunSimulator_PromethION_TimingValues(t *testing.T) {
	cfg := types.SequencerConfig{
		SequencerModel: "PromethION 24",
		FlowCellType:   "R10.4.1",
		ReadLength:     50,
		ReadsPerBatch:  5,
		TotalBatches:   1,
		PairedEnd:      false,
		SampleName:     "NanoporeTest",
		Lane:           1,
	}
	sim := timing.NewRunSimulator(cfg, timing.Compressed1000x)

	// Verify PromethION has zero BaseCallingMinutes (real-time basecalling).
	assert.Equal(t, 0, sim.Timing.BaseCallingMinutes,
		"PromethION R10.4.1 should have 0 BaseCallingMinutes (real-time)")
	assert.Equal(t, 72, sim.Timing.TotalRunTimeHours)
	assert.True(t, sim.Timing.SequencingMinutes > 0)
	// BytesPerSecond should be very high for PromethION.
	assert.True(t, sim.Timing.BytesPerSecond >= int64(800*1024*1024))
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

func basicConfig() types.SequencerConfig {
	return types.SequencerConfig{
		SequencerModel: "NovaSeq 6000",
		FlowCellType:   "S1",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		TotalBatches:   5,
		PairedEnd:      true,
		SampleName:     "TestSample",
		Lane:           1,
	}
}
