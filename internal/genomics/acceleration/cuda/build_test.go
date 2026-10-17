//go:build !release

package cuda_test

import (
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/acceleration/cuda"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestCheckCudaEnvironment_NoNVCC verifies that on systems without nvcc installed,
// CheckCudaEnvironment returns an error containing "nvcc not found".
func TestCheckCudaEnvironment_NoNVCC(t *testing.T) {
	// Save and clear PATH so nvcc cannot be found.
	origPath := os.Getenv("PATH")
	t.Cleanup(func() { os.Setenv("PATH", origPath) })

	// Set PATH to an empty directory so exec.LookPath("nvcc") fails.
	emptyDir := t.TempDir()
	require.NoError(t, os.Setenv("PATH", emptyDir))

	err := cuda.CheckCudaEnvironment()
	require.Error(t, err, "CheckCudaEnvironment must return an error when nvcc is not on PATH")
	assert.True(t,
		strings.Contains(err.Error(), "nvcc not found"),
		"error message %q must contain \"nvcc not found\"", err.Error(),
	)
}

// TestGetCudaLibraryPaths_NoPaths verifies that without CUDA installed,
// GetCudaLibraryPaths returns an error because no standard CUDA library
// directories exist on the filesystem.
func TestGetCudaLibraryPaths_NoPaths(t *testing.T) {
	// Verify none of the common CUDA paths exist before running.
	commonPaths := []string{
		"/usr/local/cuda/lib64",
		"/usr/local/cuda/lib",
		"/opt/cuda/lib64",
		"/opt/cuda/lib",
	}

	cudaPathEnv := os.Getenv("CUDA_PATH")
	anyCudaPathExists := false
	for _, p := range commonPaths {
		if _, err := os.Stat(p); err == nil {
			anyCudaPathExists = true
			break
		}
	}
	if cudaPathEnv != "" {
		for _, suffix := range []string{"lib64", "lib"} {
			p := cudaPathEnv + string(os.PathSeparator) + suffix
			if _, err := os.Stat(p); err == nil {
				anyCudaPathExists = true
				break
			}
		}
	}

	if anyCudaPathExists {
		t.Skip("At least one CUDA library path exists on this machine; skipping no-paths test")
	}

	// Clear CUDA_PATH so the env-variable branch does not inject paths either.
	origCudaPath := os.Getenv("CUDA_PATH")
	t.Cleanup(func() { os.Setenv("CUDA_PATH", origCudaPath) })
	require.NoError(t, os.Unsetenv("CUDA_PATH"))

	paths, err := cuda.GetCudaLibraryPaths()
	require.Error(t, err, "GetCudaLibraryPaths must return an error when no CUDA paths exist")
	assert.Nil(t, paths, "paths slice must be nil on error")
}

// TestGetCudaLibraryPaths_WithCUDA_PATH verifies that setting CUDA_PATH to a
// temporary directory still results in an error on machines without a real
// CUDA installation, because the expected sub-directories do not exist inside
// the temp dir.
func TestGetCudaLibraryPaths_WithCUDA_PATH(t *testing.T) {
	// Clear any real CUDA paths; if standard paths exist, the function could
	// succeed through them rather than through the env-var branch.
	standardPaths := []string{
		"/usr/local/cuda/lib64",
		"/usr/local/cuda/lib",
		"/opt/cuda/lib64",
		"/opt/cuda/lib",
	}
	for _, p := range standardPaths {
		if _, err := os.Stat(p); err == nil {
			t.Skip("Standard CUDA library path exists; cannot isolate CUDA_PATH-only behaviour")
		}
	}

	tmpDir := t.TempDir() // contains no lib64 / lib sub-directories

	origCudaPath := os.Getenv("CUDA_PATH")
	t.Cleanup(func() { os.Setenv("CUDA_PATH", origCudaPath) })
	require.NoError(t, os.Setenv("CUDA_PATH", tmpDir))

	paths, err := cuda.GetCudaLibraryPaths()
	require.Error(t, err, "GetCudaLibraryPaths must return an error when CUDA_PATH points to a dir without CUDA sub-directories")
	assert.Nil(t, paths, "paths slice must be nil on error")
}

// TestBuildCudaKernels_NoSourceFile verifies that BuildCudaKernels returns an
// error when kmer_kernels.cu does not exist in the current working directory.
func TestBuildCudaKernels_NoSourceFile(t *testing.T) {
	// Change working directory to a temporary directory that has no .cu file.
	origDir, err := os.Getwd()
	require.NoError(t, err)
	t.Cleanup(func() { os.Chdir(origDir) }) //nolint:errcheck

	tmpDir := t.TempDir()
	require.NoError(t, os.Chdir(tmpDir))

	err = cuda.BuildCudaKernels()
	require.Error(t, err, "BuildCudaKernels must return an error when kmer_kernels.cu does not exist in cwd")
	assert.True(t,
		strings.Contains(err.Error(), "CUDA source file not found") ||
			strings.Contains(err.Error(), "kmer_kernels.cu"),
		"error %q must mention the missing source file", err.Error(),
	)
}

// TestGenerateCudaBuildFlags_NoEnvironment verifies that GenerateCudaBuildFlags
// returns an error when no CUDA environment is available.
func TestGenerateCudaBuildFlags_NoEnvironment(t *testing.T) {
	// Restrict PATH so nvcc cannot be found.
	origPath := os.Getenv("PATH")
	t.Cleanup(func() { os.Setenv("PATH", origPath) })

	emptyDir := t.TempDir()
	require.NoError(t, os.Setenv("PATH", emptyDir))

	cflags, ldflags, err := cuda.GenerateCudaBuildFlags()
	require.Error(t, err, "GenerateCudaBuildFlags must return an error when CUDA is unavailable")
	assert.Empty(t, cflags, "cflags must be empty on error")
	assert.Empty(t, ldflags, "ldflags must be empty on error")
}

// TestValidateCudaBuild_NoEnvironment verifies that ValidateCudaBuild returns
// an error when no CUDA environment is present.
func TestValidateCudaBuild_NoEnvironment(t *testing.T) {
	// Restrict PATH so nvcc cannot be found.
	origPath := os.Getenv("PATH")
	t.Cleanup(func() { os.Setenv("PATH", origPath) })

	emptyDir := t.TempDir()
	require.NoError(t, os.Setenv("PATH", emptyDir))

	err := cuda.ValidateCudaBuild()
	require.Error(t, err, "ValidateCudaBuild must return an error when CUDA is unavailable")
}

// TestWriteCudaBuildConfig_NoEnvironment verifies that WriteCudaBuildConfig
// returns an error when GenerateCudaBuildFlags fails due to a missing CUDA
// environment.
func TestWriteCudaBuildConfig_NoEnvironment(t *testing.T) {
	// Restrict PATH so nvcc cannot be found, causing GenerateCudaBuildFlags to fail.
	origPath := os.Getenv("PATH")
	t.Cleanup(func() { os.Setenv("PATH", origPath) })

	emptyDir := t.TempDir()
	require.NoError(t, os.Setenv("PATH", emptyDir))

	configPath := emptyDir + "/cuda_build.mk"
	err := cuda.WriteCudaBuildConfig(configPath)
	require.Error(t, err, "WriteCudaBuildConfig must return an error when GenerateCudaBuildFlags fails")

	// Confirm no config file was written.
	_, statErr := os.Stat(configPath)
	assert.True(t, os.IsNotExist(statErr),
		"config file must not be created when WriteCudaBuildConfig returns an error")
}

// TestDetectCudaCapabilities_NoNvidiaSMI verifies that DetectCudaCapabilities
// returns an error on macOS (where nvidia-smi is not present) or on any
// system where nvidia-smi is absent from PATH.
func TestDetectCudaCapabilities_NoNvidiaSMI(t *testing.T) {
	if runtime.GOOS != "darwin" {
		// On non-macOS, only run when nvidia-smi is genuinely absent.
		origPath := os.Getenv("PATH")
		t.Cleanup(func() { os.Setenv("PATH", origPath) })

		emptyDir := t.TempDir()
		require.NoError(t, os.Setenv("PATH", emptyDir))
	}
	// On macOS, nvidia-smi is inherently absent; no PATH manipulation needed.

	err := cuda.DetectCudaCapabilities()
	require.Error(t, err, "DetectCudaCapabilities must return an error when nvidia-smi is not available")
}

// TestCudaKmerProcessor_ZeroValue_StubMethods verifies that every stub method
// on a zero-value CudaKmerProcessor behaves correctly without panicking.
// This is relevant for the !cgo and cgo&&!cuda builds that run on macOS and
// in standard CI without a CUDA SDK.
func TestCudaKmerProcessor_ZeroValue_StubMethods(t *testing.T) {
	if cuda.IsAvailable() {
		t.Skip("CUDA is available on this machine; stub behaviour test only applies to non-CUDA builds")
	}

	var proc cuda.CudaKmerProcessor

	t.Run("GetDeviceInfo", func(t *testing.T) {
		info := proc.GetDeviceInfo()
		// Stub returns zero-value DeviceInfo; must not panic.
		assert.IsType(t, cuda.DeviceInfo{}, info, "GetDeviceInfo must return a DeviceInfo value")
		assert.Empty(t, info.Name, "stub DeviceInfo.Name must be empty")
		assert.False(t, info.IsAvailable, "stub DeviceInfo.IsAvailable must be false")
	})

	t.Run("GetPerformanceStats", func(t *testing.T) {
		stats := proc.GetPerformanceStats()
		// Stub returns zero-value PerformanceStats; must not panic.
		assert.IsType(t, cuda.PerformanceStats{}, stats, "GetPerformanceStats must return a PerformanceStats value")
		assert.Empty(t, stats.DeviceName, "stub PerformanceStats.DeviceName must be empty")
		assert.Zero(t, stats.EstimatedSpeedup, "stub PerformanceStats.EstimatedSpeedup must be zero")
	})

	t.Run("Close", func(t *testing.T) {
		err := proc.Close()
		assert.NoError(t, err, "stub Close must return nil")
	})

	t.Run("ExtractKmerFrequenciesGPU", func(t *testing.T) {
		result, err := proc.ExtractKmerFrequenciesGPU([][]byte{}, 4)
		require.Error(t, err, "stub ExtractKmerFrequenciesGPU must return an error")
		assert.Nil(t, result, "stub ExtractKmerFrequenciesGPU must return nil result on error")
		assert.True(t,
			strings.Contains(err.Error(), "CUDA"),
			"error %q must mention CUDA", err.Error(),
		)
	})

	t.Run("ExtractKmerFrequenciesGPUBatch", func(t *testing.T) {
		result, err := proc.ExtractKmerFrequenciesGPUBatch([][]byte{}, 4, 100)
		require.Error(t, err, "stub ExtractKmerFrequenciesGPUBatch must return an error")
		assert.Nil(t, result, "stub ExtractKmerFrequenciesGPUBatch must return nil result on error")
		assert.True(t,
			strings.Contains(err.Error(), "CUDA"),
			"error %q must mention CUDA", err.Error(),
		)
	})

	t.Run("BenchmarkGPUPerformance", func(t *testing.T) {
		result, err := proc.BenchmarkGPUPerformance([][]byte{}, 4)
		require.Error(t, err, "stub BenchmarkGPUPerformance must return an error")
		assert.Nil(t, result, "stub BenchmarkGPUPerformance must return nil result on error")
		assert.True(t,
			strings.Contains(err.Error(), "CUDA"),
			"error %q must mention CUDA", err.Error(),
		)
	})
}

// TestGetCudaLibraryPaths_WithFakeLib64 verifies that GetCudaLibraryPaths succeeds
// when CUDA_PATH points to a directory that contains a "lib64" subdirectory.
// This exercises the success path (paths found, return paths) without a real CUDA SDK.
func TestGetCudaLibraryPaths_WithFakeLib64(t *testing.T) {
	if runtime.GOOS == "windows" {
		t.Skip("Windows path handling differs; skipping fake-lib64 test")
	}

	// Skip if any standard CUDA path already exists (would reach success through another route).
	standardPaths := []string{
		"/usr/local/cuda/lib64",
		"/usr/local/cuda/lib",
		"/opt/cuda/lib64",
		"/opt/cuda/lib",
	}
	for _, p := range standardPaths {
		if _, err := os.Stat(p); err == nil {
			t.Skipf("standard CUDA path %s exists; cannot isolate CUDA_PATH-only success path", p)
		}
	}

	// Create a fake CUDA installation directory that contains a lib64 sub-directory.
	tmpDir := t.TempDir()
	lib64Dir := filepath.Join(tmpDir, "lib64")
	require.NoError(t, os.MkdirAll(lib64Dir, 0750))

	origCudaPath := os.Getenv("CUDA_PATH")
	t.Cleanup(func() { os.Setenv("CUDA_PATH", origCudaPath) }) //nolint:errcheck
	require.NoError(t, os.Setenv("CUDA_PATH", tmpDir))

	paths, err := cuda.GetCudaLibraryPaths()
	require.NoError(t, err, "GetCudaLibraryPaths must succeed when CUDA_PATH/lib64 exists")
	require.NotEmpty(t, paths, "at least one path must be returned")
	assert.True(t,
		strings.Contains(paths[0], "lib64") || strings.Contains(paths[0], "lib"),
		"returned path %q must be a CUDA library directory", paths[0],
	)
}
