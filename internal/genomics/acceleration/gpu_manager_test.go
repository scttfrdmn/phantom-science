//go:build !release
// +build !release

package acceleration_test

import (
	"runtime"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/acceleration"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewGPUManager_Succeeds verifies that NewGPUManager never returns an error.
// GPUManager is designed to always succeed by falling back to CPU when GPU is
// unavailable.
func TestNewGPUManager_Succeeds(t *testing.T) {
	mgr, err := acceleration.NewGPUManager()
	require.NoError(t, err, "NewGPUManager() must always succeed (falls back to CPU)")
	require.NotNil(t, mgr, "returned manager must not be nil")
	defer mgr.Close() //nolint:errcheck
}

// TestNewGPUManagerWithBackend_None verifies that explicitly requesting CPU-only
// mode succeeds and yields a usable manager.
func TestNewGPUManagerWithBackend_None(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	require.NotNil(t, mgr)
	defer mgr.Close() //nolint:errcheck
}

// TestNewGPUManagerWithBackend_Auto verifies that auto-detect mode always yields
// a usable manager regardless of the GPU environment.
func TestNewGPUManagerWithBackend_Auto(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendAuto)
	require.NoError(t, err)
	require.NotNil(t, mgr)
	defer mgr.Close() //nolint:errcheck
}

// TestGPUBackend_String verifies that the String() method on every defined
// GPUBackend constant returns a non-empty, human-readable label.
func TestGPUBackend_String(t *testing.T) {
	cases := []struct {
		backend  acceleration.GPUBackend
		expected string
	}{
		{acceleration.GPUBackendNone, "CPU"},
		{acceleration.GPUBackendMetal, "Apple Metal"},
		{acceleration.GPUBackendCUDA, "NVIDIA CUDA"},
		{acceleration.GPUBackendAuto, "Auto-select"},
	}

	for _, tc := range cases {
		label := tc.backend.String()
		assert.Equal(t, tc.expected, label,
			"GPUBackend(%d).String() = %q, want %q", tc.backend, label, tc.expected)
	}
}

// TestGPUBackend_Unknown verifies that an unrecognised backend constant does not
// panic and returns the fallback label.
func TestGPUBackend_Unknown(t *testing.T) {
	unknown := acceleration.GPUBackend(999)
	label := unknown.String()
	assert.Equal(t, "Unknown", label)
}

// TestGetCapabilities_Structure verifies that GetCapabilities returns a struct
// with the expected shape regardless of GPU availability.
func TestGetCapabilities_Structure(t *testing.T) {
	mgr, err := acceleration.NewGPUManager()
	require.NoError(t, err)
	defer mgr.Close() //nolint:errcheck

	caps := mgr.GetCapabilities()

	// CPU backend must always be in AvailableBackends.
	assert.Contains(t, caps.AvailableBackends, acceleration.GPUBackendNone,
		"CPU must always be listed as an available backend")

	// DeviceInfo must always include at least the CPU entry.
	require.NotEmpty(t, caps.DeviceInfo, "DeviceInfo must contain at least a CPU entry")

	// Find the CPU device entry and validate it.
	foundCPU := false
	for _, d := range caps.DeviceInfo {
		if d.Backend == acceleration.GPUBackendNone {
			foundCPU = true
			assert.Equal(t, "CPU", d.Type, "CPU DeviceInfo.Type must be 'CPU'")
			assert.NotEmpty(t, d.Name, "CPU DeviceInfo.Name must not be empty")
			cores, ok := d.Capabilities["cores"]
			assert.True(t, ok, "CPU capabilities must include 'cores'")
			assert.Greater(t, cores.(int), 0, "CPU core count must be positive")
		}
	}
	assert.True(t, foundCPU, "DeviceInfo must include a CPU entry")
}

// TestGetCapabilities_PreferredBackend verifies that when running on macOS
// without CUDA and without CGO the preferred backend is reported as CPU or Metal.
func TestGetCapabilities_PreferredBackend(t *testing.T) {
	mgr, err := acceleration.NewGPUManager()
	require.NoError(t, err)
	defer mgr.Close() //nolint:errcheck

	caps := mgr.GetCapabilities()

	// The preferred backend must be one of the valid constants.
	validBackends := map[acceleration.GPUBackend]bool{
		acceleration.GPUBackendNone:  true,
		acceleration.GPUBackendMetal: true,
		acceleration.GPUBackendCUDA:  true,
	}
	assert.True(t, validBackends[caps.PreferredBackend],
		"PreferredBackend must be a known GPUBackend, got %v", caps.PreferredBackend)
}

// TestGetCapabilities_MacOS_MetalOrCPU verifies that on macOS the system
// correctly identifies Metal availability. With CGO disabled Metal is unavailable,
// so the preferred backend must be CPU.
func TestGetCapabilities_MacOS_MetalOrCPU(t *testing.T) {
	if runtime.GOOS != "darwin" {
		t.Skip("macOS-specific test")
	}

	mgr, err := acceleration.NewGPUManager()
	require.NoError(t, err)
	defer mgr.Close() //nolint:errcheck

	caps := mgr.GetCapabilities()

	// On macOS, CUDA is never available.
	assert.False(t, caps.CudaAvailable, "CUDA must not be available on macOS")

	// The preferred backend is Metal (CGO+darwin) or CPU (no CGO). Both are valid.
	assert.True(t,
		caps.PreferredBackend == acceleration.GPUBackendMetal ||
			caps.PreferredBackend == acceleration.GPUBackendNone,
		"on macOS the preferred backend must be Metal or CPU, got %v", caps.PreferredBackend)
}

// TestExtractKmerFrequencies_CPUFallback verifies that k-mer extraction works
// end-to-end using CPU fallback on a system without GPU acceleration.
func TestExtractKmerFrequencies_CPUFallback(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	defer mgr.Close() //nolint:errcheck

	sequences := [][]byte{
		[]byte("ACGTACGTACGT"),
		[]byte("TTTTAAAA"),
	}

	table, err := mgr.ExtractKmerFrequencies(sequences, 4)
	require.NoError(t, err, "ExtractKmerFrequencies must succeed with CPU fallback")
	require.NotNil(t, table)
	assert.Equal(t, 4, table.K, "frequency table K must match the requested k")
	assert.Greater(t, table.TotalCount, 0, "frequency table must contain k-mer counts")
}

// TestExtractKmerFrequencies_EmptySlice verifies that passing an empty sequence
// slice to ExtractKmerFrequencies is handled gracefully (no panic).
func TestExtractKmerFrequencies_EmptySlice(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	defer mgr.Close() //nolint:errcheck

	// An empty slice is a degenerate case; the function must not panic.
	// It may return an error or an empty table; both are acceptable.
	assert.NotPanics(t, func() {
		_, _ = mgr.ExtractKmerFrequencies([][]byte{}, 4)
	})
}

// TestExtractKmerFrequenciesBatch_CPUFallback verifies batched extraction with
// CPU fallback produces correct results.
func TestExtractKmerFrequenciesBatch_CPUFallback(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	defer mgr.Close() //nolint:errcheck

	sequences := [][]byte{
		[]byte("ACGTACGTACGT"),
		[]byte("GCGCGCGC"),
		[]byte("AAAAGGGGTTTTCCCC"),
	}

	table, err := mgr.ExtractKmerFrequenciesBatch(sequences, 4, 2)
	require.NoError(t, err)
	require.NotNil(t, table)
	assert.Equal(t, 4, table.K)
	assert.Greater(t, table.TotalCount, 0)
}

// TestClose_Idempotent verifies that calling Close twice does not panic or
// return an error.
func TestClose_Idempotent(t *testing.T) {
	mgr, err := acceleration.NewGPUManager()
	require.NoError(t, err)

	err1 := mgr.Close()
	assert.NoError(t, err1)

	// Second close should also be safe.
	err2 := mgr.Close()
	assert.NoError(t, err2)
}

// TestSwitchBackend_ToNone verifies that switching to the CPU backend is always
// safe and succeeds.
func TestSwitchBackend_ToNone(t *testing.T) {
	mgr, err := acceleration.NewGPUManager()
	require.NoError(t, err)
	defer mgr.Close() //nolint:errcheck

	err = mgr.SwitchBackend(acceleration.GPUBackendNone)
	assert.NoError(t, err, "switching to CPU backend must always succeed")
}

// TestSwitchBackend_SameBackend verifies that switching to the current backend
// is a no-op and returns no error.
func TestSwitchBackend_SameBackend(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	defer mgr.Close() //nolint:errcheck

	// Switch to the same backend should be a no-op.
	err = mgr.SwitchBackend(acceleration.GPUBackendNone)
	assert.NoError(t, err)
}

// TestBenchmarkPerformance_CPUOnly verifies that BenchmarkPerformance produces
// a valid result struct when only CPU is available.
func TestBenchmarkPerformance_CPUOnly(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	defer mgr.Close() //nolint:errcheck

	sequences := [][]byte{
		[]byte("ACGTACGTACGTACGT"),
		[]byte("TTTTAAAACCCCGGGG"),
	}

	benchmark, err := mgr.BenchmarkPerformance(sequences, 4)
	require.NoError(t, err, "BenchmarkPerformance must succeed with CPU fallback")
	require.NotNil(t, benchmark)

	assert.Equal(t, len(sequences), benchmark.TestSequences)
	assert.Equal(t, 4, benchmark.KmerSize)
	assert.Contains(t, benchmark.Results, acceleration.GPUBackendNone,
		"CPU benchmark result must be present")

	cpuResult := benchmark.Results[acceleration.GPUBackendNone]
	require.NotNil(t, cpuResult)
	assert.True(t, cpuResult.Success, "CPU benchmark must succeed")
	assert.Greater(t, cpuResult.TotalKmers, 0, "CPU benchmark must count some k-mers")
}
