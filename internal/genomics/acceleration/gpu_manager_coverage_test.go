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

// TestExtractKmerFrequencies_NotInitialized verifies that calling
// ExtractKmerFrequencies on a zero-value GPUManager returns an "not initialized"
// error rather than panicking.
func TestExtractKmerFrequencies_NotInitialized(t *testing.T) {
	var mgr acceleration.GPUManager
	result, err := mgr.ExtractKmerFrequencies([][]byte{[]byte("ACGT")}, 4)
	assert.Error(t, err, "ExtractKmerFrequencies on zero-value GPUManager must return an error")
	assert.Nil(t, result, "result must be nil when not initialized")
	assert.Contains(t, err.Error(), "not initialized",
		"error must mention 'not initialized', got: %v", err)
}

// TestExtractKmerFrequenciesBatch_NotInitialized verifies that calling
// ExtractKmerFrequenciesBatch on a zero-value GPUManager returns an error
// instead of panicking.
func TestExtractKmerFrequenciesBatch_NotInitialized(t *testing.T) {
	var mgr acceleration.GPUManager
	result, err := mgr.ExtractKmerFrequenciesBatch([][]byte{[]byte("ACGT")}, 4, 10)
	assert.Error(t, err, "ExtractKmerFrequenciesBatch on zero-value GPUManager must return an error")
	assert.Nil(t, result, "result must be nil when not initialized")
	assert.Contains(t, err.Error(), "not initialized",
		"error must mention 'not initialized', got: %v", err)
}

// TestBenchmarkPerformance_NotInitialized verifies that calling BenchmarkPerformance
// on a zero-value GPUManager returns an error instead of panicking.
func TestBenchmarkPerformance_NotInitialized(t *testing.T) {
	var mgr acceleration.GPUManager
	result, err := mgr.BenchmarkPerformance([][]byte{[]byte("ACGT")}, 4)
	assert.Error(t, err, "BenchmarkPerformance on zero-value GPUManager must return an error")
	assert.Nil(t, result, "result must be nil when not initialized")
	assert.Contains(t, err.Error(), "not initialized",
		"error must mention 'not initialized', got: %v", err)
}

// TestNewGPUManagerWithBackend_CUDA_MacOS verifies that requesting the CUDA backend
// on macOS gracefully falls back to CPU (GPU init will fail) and still returns a
// usable, non-nil manager.
func TestNewGPUManagerWithBackend_CUDA_MacOS(t *testing.T) {
	if runtime.GOOS != "darwin" {
		t.Skip("macOS-specific test – CUDA is never available on macOS")
	}

	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendCUDA)
	// NewGPUManagerWithBackend is designed to always succeed, falling back to
	// CPU when GPU initialization fails.
	require.NoError(t, err, "NewGPUManagerWithBackend(CUDA) must not return an error on macOS")
	require.NotNil(t, mgr, "manager must not be nil")
	defer mgr.Close() //nolint:errcheck
}

// TestNewGPUManagerWithBackend_UnknownBackend verifies that passing an unrecognised
// backend constant gracefully falls back to CPU and returns a usable manager.
func TestNewGPUManagerWithBackend_UnknownBackend(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackend(999))
	// The GPU processor init will fail for the unknown backend, but the manager
	// falls back to CPU and should not surface an error to the caller.
	require.NoError(t, err, "NewGPUManagerWithBackend(unknown) must not return an error")
	require.NotNil(t, mgr, "manager must not be nil")
	defer mgr.Close() //nolint:errcheck
}

// TestSwitchBackend_FromCPU_ToMetal verifies that switching from the CPU backend
// to Metal on macOS (where the shader library is unavailable in the test
// environment) does not panic. The error returned by SwitchBackend is acceptable.
func TestSwitchBackend_FromCPU_ToMetal(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	require.NotNil(t, mgr)
	defer mgr.Close() //nolint:errcheck

	assert.NotPanics(t, func() {
		// On macOS without a compiled Metal shader library the switch will return
		// an error; that is the expected and acceptable outcome.
		_ = mgr.SwitchBackend(acceleration.GPUBackendMetal)
	}, "SwitchBackend(Metal) must not panic")
}

// TestSwitchBackend_FromCPU_ToCUDA verifies that switching from the CPU backend
// to CUDA on a machine without CUDA returns a non-nil error.
func TestSwitchBackend_FromCPU_ToCUDA(t *testing.T) {
	if acceleration.GPUBackendCUDA.String() == "" {
		t.Skip("CUDA backend constant unavailable")
	}

	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	require.NotNil(t, mgr)
	defer mgr.Close() //nolint:errcheck

	switchErr := mgr.SwitchBackend(acceleration.GPUBackendCUDA)
	// On macOS (or any machine without CUDA) this must return an error because
	// CUDA is not available.
	if runtime.GOOS == "darwin" {
		assert.Error(t, switchErr,
			"SwitchBackend(CUDA) must return an error on macOS where CUDA is absent")
	}
	// On Linux with a real CUDA GPU the switch may succeed; we do not assert
	// error in that case so the test remains portable.
}

// TestExtractKmerFrequencies_WithAutoBackend verifies that k-mer extraction works
// end-to-end when the manager is created with GPUBackendAuto. The system
// automatically selects the best available backend; on CI this will be CPU.
func TestExtractKmerFrequencies_WithAutoBackend(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendAuto)
	require.NoError(t, err)
	require.NotNil(t, mgr)
	defer mgr.Close() //nolint:errcheck

	sequences := [][]byte{
		[]byte("ACGTACGTACGT"),
		[]byte("TTTTGGGGAAAA"),
		[]byte("CCCCGGGGTTTT"),
	}

	table, err := mgr.ExtractKmerFrequencies(sequences, 4)
	require.NoError(t, err, "ExtractKmerFrequencies must succeed via CPU fallback")
	require.NotNil(t, table)
	assert.Equal(t, 4, table.K, "frequency table K must match the requested k")
	assert.Greater(t, table.TotalCount, 0, "frequency table must contain k-mer counts")
}

// TestSwitchBackend_SameBackend_NoOp verifies that switching to the same backend
// returns nil immediately without reinitialising anything.
func TestSwitchBackend_SameBackend_NoOp(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	require.NotNil(t, mgr)
	defer mgr.Close() //nolint:errcheck

	// Switching to the current backend (None/CPU) must be a no-op returning nil.
	switchErr := mgr.SwitchBackend(acceleration.GPUBackendNone)
	assert.NoError(t, switchErr, "SwitchBackend to the same backend must return no error")
}

// TestExtractKmerFrequencies_MetalBackend_NilProcessor verifies that when the Metal
// backend is preferred but the processor is nil (shader library unavailable),
// ExtractKmerFrequencies falls back to CPU and returns valid results.
//
// Scenario: SwitchBackend(Metal) sets preferredGPU=Metal but the Metal processor
// init fails (no shader library in test env), leaving metalProcessor=nil.
// ExtractKmerFrequencies must enter the Metal case of its switch, detect the nil
// processor, and fall back to CPU without panicking or returning an error.
func TestExtractKmerFrequencies_MetalBackend_NilProcessor(t *testing.T) {
	if runtime.GOOS != "darwin" {
		t.Skip("Metal backend test is macOS-only")
	}

	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	require.NotNil(t, mgr)
	defer mgr.Close() //nolint:errcheck

	// Switch to Metal backend. On this machine Metal GPU is present but the shader
	// library is not compiled, so initializeGPUProcessors returns an error and
	// metalProcessor stays nil, but preferredGPU is set to Metal.
	_ = mgr.SwitchBackend(acceleration.GPUBackendMetal)

	sequences := [][]byte{
		[]byte("ACGTACGTACGT"),
		[]byte("TTTTGGGGAAAA"),
	}

	// Must succeed via CPU fallback (Metal case with nil processor).
	table, extractErr := mgr.ExtractKmerFrequencies(sequences, 4)
	require.NoError(t, extractErr, "ExtractKmerFrequencies must succeed via CPU fallback when Metal processor is nil")
	require.NotNil(t, table)
	assert.Equal(t, 4, table.K)
	assert.Greater(t, table.TotalCount, 0)
}

// TestExtractKmerFrequenciesBatch_ZeroBatchSize verifies that passing batchSize=0
// uses the default batch size (1000) and completes successfully via CPU processing.
func TestExtractKmerFrequenciesBatch_ZeroBatchSize(t *testing.T) {
	mgr, err := acceleration.NewGPUManagerWithBackend(acceleration.GPUBackendNone)
	require.NoError(t, err)
	require.NotNil(t, mgr)
	defer mgr.Close() //nolint:errcheck

	sequences := [][]byte{
		[]byte("ACGTACGTACGT"),
		[]byte("GCGCGCGCGCGC"),
		[]byte("TTTTAAAACCCC"),
		[]byte("GGGGCCCCAAAA"),
	}

	// batchSize=0 must trigger the default (1000) and not cause a panic or error.
	table, err := mgr.ExtractKmerFrequenciesBatch(sequences, 4, 0)
	require.NoError(t, err, "ExtractKmerFrequenciesBatch with batchSize=0 must succeed")
	require.NotNil(t, table)
	assert.Equal(t, 4, table.K, "frequency table K must match the requested k")
	assert.Greater(t, table.TotalCount, 0, "frequency table must contain k-mer counts")
}
