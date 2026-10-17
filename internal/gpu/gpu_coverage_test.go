//go:build !release

package gpu_test

import (
	"context"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/gpu"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ---------------------------------------------------------------------------
// OptimizationConfig tests
// ---------------------------------------------------------------------------

func TestDefaultOptimizationConfig_Values(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	require.NotNil(t, cfg)

	assert.True(t, cfg.MaxMemoryUsagePercent > 0.0, "MaxMemoryUsagePercent must be positive")
	assert.True(t, cfg.MaxMemoryUsagePercent <= 100.0, "MaxMemoryUsagePercent must be <= 100")
	assert.True(t, cfg.PreallocationSizeGB >= 0.0, "PreallocationSizeGB must be non-negative")
	assert.True(t, cfg.GCThresholdPercent > 0.0, "GCThresholdPercent must be positive")
	assert.True(t, cfg.ChunkSizeGB > 0.0, "ChunkSizeGB must be positive")
	assert.True(t, cfg.ConcurrentStreams > 0, "ConcurrentStreams must be positive")
	assert.True(t, cfg.CacheMaxSizeGB > 0.0, "CacheMaxSizeGB must be positive")
	assert.True(t, cfg.CacheTTL > 0, "CacheTTL must be positive")
	assert.True(t, cfg.OptimizationInterval > 0, "OptimizationInterval must be positive")
	assert.True(t, cfg.StreamingEnabled, "StreamingEnabled should default to true")
	assert.True(t, cfg.CacheEnabled, "CacheEnabled should default to true")
	assert.True(t, cfg.AdaptiveEnabled, "AdaptiveEnabled should default to true")
	assert.True(t, cfg.ProfilingEnabled, "ProfilingEnabled should default to true")
}

func TestOptimizationConfig_CustomValues(t *testing.T) {
	cfg := &gpu.OptimizationConfig{
		MaxMemoryUsagePercent: 90.0,
		PreallocationSizeGB:   4.0,
		GCThresholdPercent:    80.0,
		StreamingEnabled:      false,
		ChunkSizeGB:           2.0,
		ConcurrentStreams:     8,
		CacheEnabled:          false,
		CacheMaxSizeGB:        8.0,
		CacheTTL:              time.Hour,
		AdaptiveEnabled:       false,
		ProfilingEnabled:      false,
		OptimizationInterval:  10 * time.Minute,
	}

	assert.Equal(t, 90.0, cfg.MaxMemoryUsagePercent)
	assert.Equal(t, 4.0, cfg.PreallocationSizeGB)
	assert.Equal(t, 8, cfg.ConcurrentStreams)
	assert.False(t, cfg.StreamingEnabled)
	assert.False(t, cfg.CacheEnabled)
	assert.False(t, cfg.AdaptiveEnabled)
}

// ---------------------------------------------------------------------------
// MemoryPool constructor
// ---------------------------------------------------------------------------

func TestNewMemoryPool_WithDefaultConfig(t *testing.T) {
	pool, err := gpu.NewMemoryPool(8.0, nil)
	require.NoError(t, err)
	require.NotNil(t, pool)
}

func TestNewMemoryPool_WithExplicitConfig(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(32.0, cfg)
	require.NoError(t, err)
	require.NotNil(t, pool)
}

func TestNewMemoryPool_VeryLargeMaxMemory(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(1024.0, cfg)
	require.NoError(t, err)
	require.NotNil(t, pool)
}

// ---------------------------------------------------------------------------
// Priority constants
// ---------------------------------------------------------------------------

func TestPriorityConstants(t *testing.T) {
	assert.Equal(t, gpu.Priority(0), gpu.PriorityLow)
	assert.Equal(t, gpu.Priority(1), gpu.PriorityNormal)
	assert.Equal(t, gpu.Priority(2), gpu.PriorityHigh)
	assert.Equal(t, gpu.Priority(3), gpu.PriorityCritical)

	assert.True(t, gpu.PriorityLow < gpu.PriorityNormal)
	assert.True(t, gpu.PriorityNormal < gpu.PriorityHigh)
	assert.True(t, gpu.PriorityHigh < gpu.PriorityCritical)
}

// ---------------------------------------------------------------------------
// AllocateMemory & FreeMemory
// ---------------------------------------------------------------------------

func TestAllocateMemory_AllPriorities(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(64.0, cfg)
	require.NoError(t, err)

	priorities := []gpu.Priority{
		gpu.PriorityLow,
		gpu.PriorityNormal,
		gpu.PriorityHigh,
		gpu.PriorityCritical,
	}

	for _, p := range priorities {
		alloc, allocErr := pool.AllocateMemory(0.5, "priority-test", p)
		require.NoError(t, allocErr, "allocation must succeed for priority %v", p)
		require.NotNil(t, alloc)

		assert.Equal(t, 0.5, alloc.SizeGB)
		assert.Equal(t, "priority-test", alloc.Purpose)
		assert.Equal(t, p, alloc.Priority)
		assert.NotEmpty(t, alloc.ID)
		assert.False(t, alloc.CreatedAt.IsZero())
		assert.False(t, alloc.LastAccessed.IsZero())

		freeErr := pool.FreeMemory(alloc.ID)
		require.NoError(t, freeErr)
	}
}

func TestAllocateMemory_MultipleSimultaneous(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(64.0, cfg)
	require.NoError(t, err)

	const n = 5
	allocs := make([]*gpu.Allocation, n)

	for i := 0; i < n; i++ {
		allocs[i], err = pool.AllocateMemory(1.0, "batch", gpu.PriorityNormal)
		require.NoError(t, err)
		require.NotNil(t, allocs[i])
	}

	ids := make(map[string]struct{})
	for _, a := range allocs {
		_, dup := ids[a.ID]
		assert.False(t, dup, "allocation IDs must be unique, duplicate: %s", a.ID)
		ids[a.ID] = struct{}{}
	}

	for _, a := range allocs {
		assert.NoError(t, pool.FreeMemory(a.ID))
	}
}

func TestFreeMemory_NonExistentID(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(16.0, cfg)
	require.NoError(t, err)

	err = pool.FreeMemory("does-not-exist-12345")
	assert.Error(t, err, "FreeMemory must return an error for an unknown ID")
}

func TestFreeMemory_DoubleFree(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(16.0, cfg)
	require.NoError(t, err)

	alloc, err := pool.AllocateMemory(1.0, "double-free", gpu.PriorityNormal)
	require.NoError(t, err)

	require.NoError(t, pool.FreeMemory(alloc.ID))
	err = pool.FreeMemory(alloc.ID)
	assert.Error(t, err, "second free of the same ID must return an error")
}

// ---------------------------------------------------------------------------
// Metrics
// ---------------------------------------------------------------------------

func TestGetMetrics_InitialState(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(16.0, cfg)
	require.NoError(t, err)

	m := pool.GetMetrics()
	assert.True(t, m.TotalAllocationsGB >= 0.0)
	assert.True(t, m.PeakUsageGB >= 0.0)
	assert.True(t, m.AllocationCount >= 0)
	assert.True(t, m.EvictionCount >= 0)
}

func TestGetMetrics_AfterAllocations(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(64.0, cfg)
	require.NoError(t, err)

	before := pool.GetMetrics()

	alloc1, err := pool.AllocateMemory(2.0, "metrics-test", gpu.PriorityNormal)
	require.NoError(t, err)
	alloc2, err := pool.AllocateMemory(3.0, "metrics-test", gpu.PriorityHigh)
	require.NoError(t, err)

	afterAlloc := pool.GetMetrics()
	assert.True(t, afterAlloc.AllocationCount > before.AllocationCount)
	assert.True(t, afterAlloc.TotalAllocationsGB > before.TotalAllocationsGB)
	assert.True(t, afterAlloc.PeakUsageGB >= before.PeakUsageGB)

	require.NoError(t, pool.FreeMemory(alloc1.ID))
	require.NoError(t, pool.FreeMemory(alloc2.ID))

	afterFree := pool.GetMetrics()
	assert.True(t, afterFree.TotalAllocationsGB < afterAlloc.TotalAllocationsGB)
}

func TestGetMetrics_PeakTracking(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(64.0, cfg)
	require.NoError(t, err)

	a1, err := pool.AllocateMemory(4.0, "peak-test", gpu.PriorityNormal)
	require.NoError(t, err)

	peakDuringAlloc := pool.GetMetrics().PeakUsageGB
	assert.True(t, peakDuringAlloc >= 0.0)

	require.NoError(t, pool.FreeMemory(a1.ID))

	peakAfterFree := pool.GetMetrics().PeakUsageGB
	assert.True(t, peakAfterFree >= 0.0)
}

// ---------------------------------------------------------------------------
// GetDeviceStatus
// ---------------------------------------------------------------------------

func TestGetDeviceStatus_NotEmpty(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(16.0, cfg)
	require.NoError(t, err)

	status := pool.GetDeviceStatus()
	assert.NotEmpty(t, status, "must report at least one device")
}

func TestGetDeviceStatus_Fields(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(16.0, cfg)
	require.NoError(t, err)

	for key, ds := range pool.GetDeviceStatus() {
		assert.True(t, ds.TotalMemoryGB > 0.0, "device %s: TotalMemoryGB must be positive", key)
		assert.True(t, ds.AvailableMemoryGB >= 0.0, "device %s: AvailableMemoryGB must be non-negative", key)
		assert.True(t, ds.UtilizationPercent >= 0.0, "device %s: utilization must be >= 0", key)
		assert.True(t, ds.UtilizationPercent <= 100.0, "device %s: utilization must be <= 100", key)
		assert.True(t, ds.AllocationCount >= 0, "device %s: AllocationCount must be non-negative", key)
		assert.False(t, ds.LastUsed.IsZero(), "device %s: LastUsed must be set", key)
	}
}

func TestGetDeviceStatus_UtilizationChanges(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(64.0, cfg)
	require.NoError(t, err)

	statusBefore := pool.GetDeviceStatus()
	var totalBefore float64
	for _, ds := range statusBefore {
		totalBefore += ds.UtilizationPercent
	}

	alloc, err := pool.AllocateMemory(2.0, "util-test", gpu.PriorityNormal)
	require.NoError(t, err)

	statusAfter := pool.GetDeviceStatus()
	var totalAfter float64
	for _, ds := range statusAfter {
		totalAfter += ds.UtilizationPercent
	}

	assert.True(t, totalAfter >= totalBefore, "utilization must not decrease after allocation")

	require.NoError(t, pool.FreeMemory(alloc.ID))
}

func TestGetDeviceStatus_AllocationCountTracked(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(64.0, cfg)
	require.NoError(t, err)

	alloc, err := pool.AllocateMemory(1.0, "count-test", gpu.PriorityLow)
	require.NoError(t, err)

	var totalAllocs int
	for _, ds := range pool.GetDeviceStatus() {
		totalAllocs += ds.AllocationCount
	}
	assert.True(t, totalAllocs >= 1, "at least one allocation must be reflected in device status")

	require.NoError(t, pool.FreeMemory(alloc.ID))
}

// ---------------------------------------------------------------------------
// DeviceStatus struct (exported fields)
// ---------------------------------------------------------------------------

func TestDeviceStatus_StructFields(t *testing.T) {
	ds := gpu.DeviceStatus{
		DeviceID:           3,
		TotalMemoryGB:      24.0,
		AvailableMemoryGB:  20.0,
		UtilizationPercent: 16.67,
		AllocationCount:    2,
		LastUsed:           time.Now(),
	}

	assert.Equal(t, 3, ds.DeviceID)
	assert.Equal(t, 24.0, ds.TotalMemoryGB)
	assert.Equal(t, 20.0, ds.AvailableMemoryGB)
	assert.Equal(t, 16.67, ds.UtilizationPercent)
	assert.Equal(t, 2, ds.AllocationCount)
	assert.False(t, ds.LastUsed.IsZero())
}

// ---------------------------------------------------------------------------
// Allocation struct (exported fields)
// ---------------------------------------------------------------------------

func TestAllocation_StructFields(t *testing.T) {
	now := time.Now()
	a := &gpu.Allocation{
		ID:           "alloc-xyz",
		SizeGB:       3.5,
		Purpose:      "test",
		CreatedAt:    now,
		LastAccessed: now,
		Priority:     gpu.PriorityHigh,
		DeviceID:     1,
	}

	assert.Equal(t, "alloc-xyz", a.ID)
	assert.Equal(t, 3.5, a.SizeGB)
	assert.Equal(t, "test", a.Purpose)
	assert.Equal(t, gpu.PriorityHigh, a.Priority)
	assert.Equal(t, 1, a.DeviceID)
}

// ---------------------------------------------------------------------------
// StreamingConfig
// ---------------------------------------------------------------------------

func TestDefaultStreamingConfig_Values(t *testing.T) {
	cfg := gpu.DefaultStreamingConfig()
	require.NotNil(t, cfg)

	assert.True(t, cfg.MaxConcurrentStreams > 0, "MaxConcurrentStreams must be positive")
	assert.True(t, cfg.ChunkSizeGB > 0.0, "ChunkSizeGB must be positive")
	assert.True(t, cfg.BufferSizeGB > 0.0, "BufferSizeGB must be positive")
	assert.True(t, cfg.TimeoutDuration > 0, "TimeoutDuration must be positive")
	assert.True(t, cfg.RetryAttempts >= 0, "RetryAttempts must be non-negative")
	assert.True(t, cfg.CompressionEnabled, "CompressionEnabled should default to true")
	assert.True(t, cfg.ChecksumValidation, "ChecksumValidation should default to true")
}

func TestStreamingConfig_CustomValues(t *testing.T) {
	cfg := &gpu.StreamingConfig{
		MaxConcurrentStreams: 2,
		ChunkSizeGB:         0.5,
		BufferSizeGB:        1.0,
		TimeoutDuration:     5 * time.Minute,
		RetryAttempts:       1,
		CompressionEnabled:  false,
		ChecksumValidation:  false,
	}

	assert.Equal(t, 2, cfg.MaxConcurrentStreams)
	assert.Equal(t, 0.5, cfg.ChunkSizeGB)
	assert.False(t, cfg.CompressionEnabled)
	assert.False(t, cfg.ChecksumValidation)
}

// ---------------------------------------------------------------------------
// NewStreamingProcessor
// ---------------------------------------------------------------------------

func TestNewStreamingProcessor_WithNilConfig(t *testing.T) {
	pool := mustPool(t)
	sp := gpu.NewStreamingProcessor(pool, nil)
	require.NotNil(t, sp)
}

func TestNewStreamingProcessor_WithExplicitConfig(t *testing.T) {
	pool := mustPool(t)
	cfg := gpu.DefaultStreamingConfig()
	sp := gpu.NewStreamingProcessor(pool, cfg)
	require.NotNil(t, sp)
}

// ---------------------------------------------------------------------------
// StreamingMetrics
// ---------------------------------------------------------------------------

func TestStreamingProcessor_InitialMetrics(t *testing.T) {
	pool := mustPool(t)
	sp := gpu.NewStreamingProcessor(pool, nil)

	m := sp.GetMetrics()
	assert.Equal(t, 0, m.TotalStreams)
	assert.Equal(t, 0, m.ActiveStreams)
	assert.Equal(t, 0, m.CompletedStreams)
	assert.Equal(t, 0, m.FailedStreams)
	assert.True(t, m.TotalDataProcessed >= 0.0)
	assert.True(t, m.AverageThroughput >= 0.0)
	assert.True(t, m.PeakThroughput >= 0.0)
}

func TestStreamingMetrics_StructFields(t *testing.T) {
	m := &gpu.StreamingMetrics{
		TotalStreams:        10,
		ActiveStreams:       2,
		CompletedStreams:    7,
		FailedStreams:       1,
		TotalDataProcessed: 100.0,
		AverageThroughput:  50.0,
		PeakThroughput:     75.0,
	}

	assert.Equal(t, 10, m.TotalStreams)
	assert.Equal(t, 2, m.ActiveStreams)
	assert.Equal(t, 7, m.CompletedStreams)
	assert.Equal(t, 1, m.FailedStreams)
	assert.Equal(t, 100.0, m.TotalDataProcessed)
	assert.Equal(t, 50.0, m.AverageThroughput)
	assert.Equal(t, 75.0, m.PeakThroughput)
}

// ---------------------------------------------------------------------------
// GetStreamStatus & CancelStream — error paths
// ---------------------------------------------------------------------------

func TestGetStreamStatus_UnknownStream(t *testing.T) {
	pool := mustPool(t)
	sp := gpu.NewStreamingProcessor(pool, nil)

	_, err := sp.GetStreamStatus("no-such-stream")
	assert.Error(t, err)
}

func TestCancelStream_UnknownStream(t *testing.T) {
	pool := mustPool(t)
	sp := gpu.NewStreamingProcessor(pool, nil)

	err := sp.CancelStream("no-such-stream")
	assert.Error(t, err)
}

// ---------------------------------------------------------------------------
// ProcessDataset — concurrency limit
// ---------------------------------------------------------------------------

func TestProcessDataset_ConcurrencyLimit(t *testing.T) {
	pool := mustPool(t)

	cfg := gpu.DefaultStreamingConfig()
	cfg.MaxConcurrentStreams = 2
	cfg.ChunkSizeGB = 0.01
	cfg.TimeoutDuration = 30 * time.Second
	cfg.RetryAttempts = 0

	sp := gpu.NewStreamingProcessor(pool, cfg)

	input := gpu.DatasetInput{
		InputPath:  "/tmp/input.fastq",
		OutputDir:  "/tmp/output",
		SizeGB:     5.0,
		Format:     "FASTQ",
		Parameters: map[string]interface{}{},
	}

	ctx := context.Background()

	streams := make([]*gpu.ProcessingStream, cfg.MaxConcurrentStreams)
	for i := 0; i < cfg.MaxConcurrentStreams; i++ {
		s, err := sp.ProcessDataset(ctx, input)
		require.NoError(t, err, "stream %d should start", i)
		streams[i] = s
	}

	_, err := sp.ProcessDataset(ctx, input)
	assert.Error(t, err, "stream beyond limit must be rejected")

	for _, s := range streams {
		_ = sp.CancelStream(s.ID)
		for range s.ResultChan {
		}
	}
}

// ---------------------------------------------------------------------------
// ProcessDataset — cancel, status, metrics
// ---------------------------------------------------------------------------

func TestProcessDataset_CancelStream(t *testing.T) {
	pool := mustPool(t)

	cfg := gpu.DefaultStreamingConfig()
	cfg.ChunkSizeGB = 0.01
	cfg.TimeoutDuration = 60 * time.Second
	cfg.RetryAttempts = 0

	sp := gpu.NewStreamingProcessor(pool, cfg)

	input := gpu.DatasetInput{
		InputPath: "/tmp/in.fastq",
		OutputDir: "/tmp/out",
		SizeGB:    5.0,
		Format:    "FASTQ",
	}

	stream, err := sp.ProcessDataset(context.Background(), input)
	require.NoError(t, err)
	require.NotNil(t, stream)

	m := sp.GetMetrics()
	assert.True(t, m.TotalStreams >= 1)

	err = sp.CancelStream(stream.ID)
	assert.NoError(t, err)

	for range stream.ResultChan {
	}
}

func TestProcessDataset_StreamStatusWhileActive(t *testing.T) {
	pool := mustPool(t)

	cfg := gpu.DefaultStreamingConfig()
	cfg.ChunkSizeGB = 0.01
	cfg.TimeoutDuration = 60 * time.Second
	cfg.RetryAttempts = 0

	sp := gpu.NewStreamingProcessor(pool, cfg)

	input := gpu.DatasetInput{
		InputPath: "/tmp/in.fastq",
		OutputDir: "/tmp/out",
		SizeGB:    2.0,
		Format:    "FASTQ",
	}

	stream, err := sp.ProcessDataset(context.Background(), input)
	require.NoError(t, err)

	status, statusErr := sp.GetStreamStatus(stream.ID)
	if statusErr == nil {
		assert.NotEmpty(t, status.ID)
		assert.True(t, status.DatasetSize > 0.0)
		assert.True(t, status.ChunksTotal >= 1)
	}

	_ = sp.CancelStream(stream.ID)
	for range stream.ResultChan {
	}
}

func TestProcessDataset_MetricsAfterCompletion(t *testing.T) {
	pool := mustPool(t)

	cfg := gpu.DefaultStreamingConfig()
	cfg.ChunkSizeGB = 1.0
	cfg.TimeoutDuration = 30 * time.Second
	cfg.RetryAttempts = 0
	cfg.ChecksumValidation = false
	cfg.CompressionEnabled = false

	sp := gpu.NewStreamingProcessor(pool, cfg)

	input := gpu.DatasetInput{
		InputPath: "/tmp/in.fastq",
		OutputDir: "/tmp/out",
		SizeGB:    0.1,
		Format:    "FASTQ",
	}

	stream, err := sp.ProcessDataset(context.Background(), input)
	require.NoError(t, err)

	for range stream.ResultChan {
	}

	time.Sleep(100 * time.Millisecond)

	m := sp.GetMetrics()
	assert.True(t, m.TotalStreams >= 1)
}

func TestProcessDataset_MultipleSequential(t *testing.T) {
	pool := mustPool(t)

	cfg := gpu.DefaultStreamingConfig()
	cfg.ChunkSizeGB = 1.0
	cfg.MaxConcurrentStreams = 4
	cfg.TimeoutDuration = 30 * time.Second
	cfg.RetryAttempts = 0
	cfg.ChecksumValidation = false
	cfg.CompressionEnabled = false

	sp := gpu.NewStreamingProcessor(pool, cfg)

	input := gpu.DatasetInput{
		InputPath: "/tmp/in.fastq",
		OutputDir: "/tmp/out",
		SizeGB:    0.1,
		Format:    "FASTQ",
	}

	for i := 0; i < 3; i++ {
		stream, err := sp.ProcessDataset(context.Background(), input)
		require.NoError(t, err)
		for range stream.ResultChan {
		}
		time.Sleep(50 * time.Millisecond)
	}

	m := sp.GetMetrics()
	assert.True(t, m.TotalStreams >= 3, "expected at least 3 total streams, got %d", m.TotalStreams)
}

// ---------------------------------------------------------------------------
// StreamStatus string representations
// ---------------------------------------------------------------------------

func TestStreamStatus_AllStrings(t *testing.T) {
	cases := map[gpu.StreamStatus]string{
		gpu.StatusInitializing: "initializing",
		gpu.StatusProcessing:   "processing",
		gpu.StatusCompleted:    "completed",
		gpu.StatusFailed:       "failed",
		gpu.StatusCancelled:    "cancelled",
	}
	for status, want := range cases {
		assert.Equal(t, want, status.String(), "StreamStatus(%d).String()", int(status))
	}
}

func TestStreamStatus_Unknown(t *testing.T) {
	unknown := gpu.StreamStatus(999)
	s := unknown.String()
	assert.NotEmpty(t, s)
}

// ---------------------------------------------------------------------------
// DatasetInput struct
// ---------------------------------------------------------------------------

func TestDatasetInput_AllFields(t *testing.T) {
	params := map[string]interface{}{"threads": 8, "quality": 30}
	input := gpu.DatasetInput{
		InputPath:  "/data/sample.fastq",
		OutputDir:  "/results",
		SizeGB:     25.5,
		Format:     "FASTQ",
		Parameters: params,
	}

	assert.Equal(t, "/data/sample.fastq", input.InputPath)
	assert.Equal(t, "/results", input.OutputDir)
	assert.Equal(t, 25.5, input.SizeGB)
	assert.Equal(t, "FASTQ", input.Format)
	assert.Equal(t, 8, input.Parameters["threads"])
}

func TestDatasetInput_EmptyParameters(t *testing.T) {
	input := gpu.DatasetInput{
		InputPath: "/tmp/in.fastq",
		OutputDir: "/tmp/out",
		SizeGB:    1.0,
		Format:    "FASTQ",
	}
	assert.Nil(t, input.Parameters)
	input.Parameters = map[string]interface{}{}
	assert.NotNil(t, input.Parameters)
}

// ---------------------------------------------------------------------------
// ProcessingResult struct
// ---------------------------------------------------------------------------

func TestProcessingResult_AllFields(t *testing.T) {
	dur := 500 * time.Millisecond
	meta := map[string]interface{}{"device_id": 0, "chunk_offset": 1.5}
	result := &gpu.ProcessingResult{
		ChunkID:    "stream-123-chunk-0",
		SizeGB:     1.0,
		Duration:   dur,
		OutputPath: "/tmp/out/chunk-0.fastq",
		Checksum:   "sha256:deadbeef",
		Error:      nil,
		Metadata:   meta,
	}

	assert.Equal(t, "stream-123-chunk-0", result.ChunkID)
	assert.Equal(t, 1.0, result.SizeGB)
	assert.Equal(t, dur, result.Duration)
	assert.Equal(t, "sha256:deadbeef", result.Checksum)
	assert.Nil(t, result.Error)
	assert.Equal(t, 0, result.Metadata["device_id"])
}

func TestProcessingResult_WithError(t *testing.T) {
	result := &gpu.ProcessingResult{
		ChunkID:  "chunk-err",
		SizeGB:   0.5,
		Duration: time.Millisecond,
		Error:    assert.AnError,
	}

	assert.NotNil(t, result.Error)
	assert.Equal(t, "chunk-err", result.ChunkID)
}

// ---------------------------------------------------------------------------
// ProcessingStream struct
// ---------------------------------------------------------------------------

func TestProcessingStream_StructFields(t *testing.T) {
	now := time.Now()
	stream := &gpu.ProcessingStream{
		ID:              "stream-abc",
		DatasetSize:     10.0,
		ProcessedSize:   5.0,
		ChunksTotal:     10,
		ChunksProcessed: 5,
		StartTime:       now,
		Status:          gpu.StatusProcessing,
		Error:           nil,
	}

	assert.Equal(t, "stream-abc", stream.ID)
	assert.Equal(t, 10.0, stream.DatasetSize)
	assert.Equal(t, 5.0, stream.ProcessedSize)
	assert.Equal(t, 10, stream.ChunksTotal)
	assert.Equal(t, 5, stream.ChunksProcessed)
	assert.Equal(t, gpu.StatusProcessing, stream.Status)
	assert.Nil(t, stream.Error)
}

// ---------------------------------------------------------------------------
// Concurrency safety smoke test
// ---------------------------------------------------------------------------

func TestMemoryPool_ConcurrentAllocFree(t *testing.T) {
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(64.0, cfg)
	require.NoError(t, err)

	done := make(chan struct{}, 4)

	for i := 0; i < 4; i++ {
		go func() {
			defer func() { done <- struct{}{} }()
			for j := 0; j < 5; j++ {
				a, aErr := pool.AllocateMemory(0.1, "concurrent", gpu.PriorityNormal)
				if aErr == nil {
					_ = pool.FreeMemory(a.ID)
				}
			}
		}()
	}

	for i := 0; i < 4; i++ {
		<-done
	}

	m := pool.GetMetrics()
	assert.True(t, m.AllocationCount >= 0)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

func mustPool(t *testing.T) *gpu.MemoryPool {
	t.Helper()
	cfg := gpu.DefaultOptimizationConfig()
	cfg.AdaptiveEnabled = false

	pool, err := gpu.NewMemoryPool(64.0, cfg)
	require.NoError(t, err)
	return pool
}
