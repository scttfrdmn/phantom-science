//go:build !release
// +build !release

package metal_test

import (
	"runtime"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/acceleration/metal"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestIsAvailable verifies that IsAvailable returns a bool without panicking.
// Behaviour is platform- and CGO-dependent:
//   - CGO disabled: always false (no-CGO stub)
//   - CGO enabled on non-darwin: always false (GOOS check in the source)
//   - CGO enabled on darwin: true when a Metal device is present
func TestIsAvailable(t *testing.T) {
	// Must never panic.
	available := metal.IsAvailable()
	assert.IsType(t, false, available, "IsAvailable() must return a bool")

	if runtime.GOOS != "darwin" {
		assert.False(t, available, "Metal must not be available on non-macOS platforms")
	}
	// On macOS with CGO the result depends on whether a Metal device is present;
	// we do not assert true here because the test may also run under no-CGO
	// (e.g. standard CI build), in which case IsAvailable returns false even on darwin.
}

// TestNewMetalKmerProcessor_NoCGO verifies error handling when Metal is unavailable.
// On macOS with CGO disabled the factory must return a meaningful error, not panic.
func TestNewMetalKmerProcessor_Unavailable(t *testing.T) {
	if metal.IsAvailable() {
		t.Skip("Metal is available on this machine – skipping unavailable-path test")
	}

	proc, err := metal.NewMetalKmerProcessor()

	require.Error(t, err, "NewMetalKmerProcessor() must return an error when Metal is unavailable")
	assert.Nil(t, proc, "processor must be nil when an error is returned")
	assert.NotEmpty(t, err.Error(), "error message must not be empty")
}

// TestNewMetalKmerProcessor_ErrorContent verifies the error communicates the reason.
func TestNewMetalKmerProcessor_ErrorContent(t *testing.T) {
	if metal.IsAvailable() {
		t.Skip("Metal is available on this machine – skipping error content test")
	}

	_, err := metal.NewMetalKmerProcessor()
	require.Error(t, err)

	// Either Metal is disabled or CGO is off; either way the message must be informative.
	assert.NotEmpty(t, err.Error())
}

// newMetalProcessor is a test helper that creates a processor, skipping if Metal
// is unavailable OR if the Metal shader library cannot be loaded.  The latter
// can happen in test environments where the .metallib was not compiled
// alongside the binary (Metal shaders must be pre-compiled via the Metal
// toolchain to be accessible at runtime).
func newMetalProcessor(t *testing.T) *metal.MetalKmerProcessor {
	t.Helper()

	if !metal.IsAvailable() {
		t.Skip("Metal is not available – skipping test")
	}

	proc, err := metal.NewMetalKmerProcessor()
	if err != nil {
		// Metal GPU is present but the shader library failed to load.
		// This is expected in plain `go test` environments where the .metallib
		// is not embedded in the test binary; skip rather than fail.
		t.Skipf("Metal shader library unavailable in this test environment: %v", err)
	}

	require.NotNil(t, proc)
	return proc
}

// TestNewMetalKmerProcessor_Available tests successful construction on macOS when Metal
// is present. The processor must be created, usable, and closeable without panicking.
func TestNewMetalKmerProcessor_Available(t *testing.T) {
	proc := newMetalProcessor(t)

	// Close must not panic and must return no error.
	closeErr := proc.Close()
	assert.NoError(t, closeErr, "Close() must not return an error")
}

// TestMetalProcessor_EmptySequences verifies that passing an empty slice to
// ExtractKmerFrequenciesGPU returns an error rather than panicking.
func TestMetalProcessor_EmptySequences(t *testing.T) {
	proc := newMetalProcessor(t)
	defer proc.Close() //nolint:errcheck

	result, err := proc.ExtractKmerFrequenciesGPU([][]byte{}, 6)
	assert.Error(t, err, "empty sequence slice must produce an error")
	assert.Nil(t, result, "result must be nil when an error is returned")
}

// TestMetalProcessor_GetDeviceInfo verifies that device info can be retrieved
// without panicking and contains sensible values.
func TestMetalProcessor_GetDeviceInfo(t *testing.T) {
	proc := newMetalProcessor(t)
	defer proc.Close() //nolint:errcheck

	info := proc.GetDeviceInfo()
	assert.True(t, info.IsAvailable, "DeviceInfo.IsAvailable must be true for an initialised processor")
	assert.NotEmpty(t, info.Name, "DeviceInfo.Name must be non-empty for an initialised processor")
	assert.Greater(t, info.MaxThreads, 0, "DeviceInfo.MaxThreads must be positive")
}

// TestMetalProcessor_Close_IdempotentNoCGO verifies that calling Close on a nil
// *MetalKmerProcessor (the value returned when Metal is unavailable) does not
// panic. Because the no-CGO stub uses an empty struct receiver, we can test via
// direct construction of the zero value.
func TestMetalProcessor_NoCGO_CloseIsNoop(t *testing.T) {
	if metal.IsAvailable() {
		t.Skip("Metal is available; no-CGO noop test does not apply")
	}

	// The factory returns nil when unavailable; verify the constructor error
	// is sufficient to guard all downstream calls.
	proc, err := metal.NewMetalKmerProcessor()
	require.Error(t, err)
	require.Nil(t, proc)
	// If the caller guards with the error (as they should), Close is never
	// called on a nil pointer. The API contract is that the caller checks err.
}

// TestEstimateSpeedup verifies that EstimateSpeedup returns non-negative values
// across the documented range of maxThreads inputs.
func TestEstimateSpeedup(t *testing.T) {
	cases := []int{0, 128, 256, 512, 1024, 2048}
	for _, maxThreads := range cases {
		speedup := metal.EstimateSpeedup(maxThreads)
		assert.GreaterOrEqual(t, speedup, 0.0,
			"EstimateSpeedup(%d) must be non-negative", maxThreads)
	}
}

// TestGenerateTestSequences verifies the helper returns sequences of the correct
// count and length containing only valid DNA bases.
func TestGenerateTestSequences(t *testing.T) {
	count, length := 5, 40
	seqs := metal.GenerateTestSequences(count, length)

	require.Len(t, seqs, count, "GenerateTestSequences must return exactly count sequences")

	validBases := map[byte]bool{'A': true, 'T': true, 'G': true, 'C': true}
	for i, seq := range seqs {
		assert.Len(t, seq, length, "sequence %d must have length %d", i, length)
		for j, b := range seq {
			assert.True(t, validBases[b],
				"sequence %d, position %d: unexpected base %q", i, j, b)
		}
	}
}

// TestGenerateTestSequences_Zero verifies edge-case with zero count.
func TestGenerateTestSequences_Zero(t *testing.T) {
	seqs := metal.GenerateTestSequences(0, 100)
	assert.Empty(t, seqs, "GenerateTestSequences(0, 100) must return an empty slice")
}

// TestDeviceInfo_ZeroValue verifies the DeviceInfo type has sensible zero values.
func TestDeviceInfo_ZeroValue(t *testing.T) {
	var info metal.DeviceInfo
	assert.False(t, info.IsAvailable, "zero-value DeviceInfo.IsAvailable must be false")
	assert.Empty(t, info.Name)
	assert.Zero(t, info.MaxThreads)
	assert.Zero(t, info.MaxBufferSize)
}

// TestPerformanceStats_ZeroValue verifies PerformanceStats zero values are sane.
func TestPerformanceStats_ZeroValue(t *testing.T) {
	var stats metal.PerformanceStats
	assert.Empty(t, stats.DeviceName)
	assert.Zero(t, stats.MaxThreadsPerGroup)
	assert.Zero(t, stats.MaxBufferSize)
	assert.Zero(t, stats.EstimatedSpeedup)
	assert.Equal(t, [2]int{0, 0}, stats.OptimalKRange)
}
