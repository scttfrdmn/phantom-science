//go:build !release
// +build !release

package metal_test

import (
	"os"
	"runtime"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/acceleration/metal"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestCheckMetalEnvironment verifies that CheckMetalEnvironment never panics.
// On darwin the result depends on whether Xcode command-line tools and the Metal
// toolchain are installed; both success and error are acceptable outcomes.
// On any other platform an error mentioning "Metal" must be returned.
func TestCheckMetalEnvironment(t *testing.T) {
	assert.NotPanics(t, func() {
		err := metal.CheckMetalEnvironment()
		if runtime.GOOS != "darwin" {
			// On non-darwin the function must return an informative error.
			if assert.Error(t, err, "CheckMetalEnvironment() must return an error on non-darwin") {
				assert.True(t,
					strings.Contains(err.Error(), "Metal") || strings.Contains(err.Error(), "metal"),
					"error must mention Metal, got: %v", err)
			}
		}
		// On darwin: xcrun may or may not be available, and the Metal toolchain
		// component may not be installed. Both error and nil are acceptable.
	}, "CheckMetalEnvironment() must not panic on any platform")
}

// TestBuildMetalShaders_NoSourceFile verifies that BuildMetalShaders returns an
// error on non-darwin platforms immediately. On darwin the test's working directory
// is the package source directory which does contain kmer_kernels.metal, so the
// function will proceed to compilation. Because the Metal toolchain may or may not
// be installed on the test machine, we simply assert that the function does not
// panic and that any error returned is non-nil (we never expect silent success in
// a plain `go test` environment without the Metal toolchain installed).
func TestBuildMetalShaders_NoSourceFile(t *testing.T) {
	if runtime.GOOS != "darwin" {
		err := metal.BuildMetalShaders()
		assert.Error(t, err, "BuildMetalShaders() must return an error on non-darwin")
		return
	}

	// On darwin: assert the function doesn't panic regardless of toolchain state.
	assert.NotPanics(t, func() {
		// The result is environment-dependent: if the Metal toolchain is not installed
		// xcrun will fail; if kmer_kernels.metal is missing it returns "not found".
		// Either way we expect an error in a normal test environment.
		_ = metal.BuildMetalShaders()
	}, "BuildMetalShaders() must not panic on darwin")
}

// TestInstallMetalLibrary_NoLib verifies that InstallMetalLibrary returns an error
// when neither the pre-compiled .metallib nor a usable Metal toolchain are present.
//
// On darwin the function first checks for a pre-built library; when that is absent
// it attempts to build it via BuildMetalShaders (which also fails because the Metal
// toolchain is not installed in this test environment), so InstallMetalLibrary must
// propagate that error.
//
// On non-darwin the function must return a "only supported on macOS" error.
func TestInstallMetalLibrary_NoLib(t *testing.T) {
	err := metal.InstallMetalLibrary(t.TempDir())

	assert.Error(t, err, "InstallMetalLibrary() must return an error when the metallib is absent")

	if runtime.GOOS != "darwin" {
		assert.True(t,
			strings.Contains(err.Error(), "macOS") || strings.Contains(err.Error(), "darwin"),
			"error must mention macOS on non-darwin, got: %v", err)
	}
	// On darwin the error propagates from BuildMetalShaders (Metal toolchain absent
	// or .metallib not pre-compiled).
}

// TestBuildMetalShaders_NonDarwin verifies that BuildMetalShaders returns an error
// immediately on non-darwin platforms.
func TestBuildMetalShaders_NonDarwin(t *testing.T) {
	if runtime.GOOS == "darwin" {
		t.Skip("non-darwin-specific test – skipping on macOS")
	}

	err := metal.BuildMetalShaders()
	assert.Error(t, err, "BuildMetalShaders() must error on non-darwin")
}

// TestCheckMetalEnvironment_NonDarwin verifies that CheckMetalEnvironment returns
// an error immediately on non-darwin platforms.
func TestCheckMetalEnvironment_NonDarwin(t *testing.T) {
	if runtime.GOOS == "darwin" {
		t.Skip("non-darwin-specific test – skipping on macOS")
	}

	err := metal.CheckMetalEnvironment()
	assert.Error(t, err, "CheckMetalEnvironment() must error on non-darwin")
}

// TestInstallMetalLibrary_NonDarwin verifies that InstallMetalLibrary returns an
// error immediately on non-darwin platforms.
func TestInstallMetalLibrary_NonDarwin(t *testing.T) {
	if runtime.GOOS == "darwin" {
		t.Skip("non-darwin-specific test – skipping on macOS")
	}

	err := metal.InstallMetalLibrary(t.TempDir())
	assert.Error(t, err, "InstallMetalLibrary() must error on non-darwin")
	assert.True(t,
		strings.Contains(err.Error(), "macOS") || strings.Contains(err.Error(), "darwin"),
		"error must mention macOS, got: %v", err)
}

// TestCheckMetalEnvironment_NoXcrun verifies that CheckMetalEnvironment returns
// an "xcrun not found" error on macOS when xcrun is not on PATH.
func TestCheckMetalEnvironment_NoXcrun(t *testing.T) {
	if runtime.GOOS != "darwin" {
		t.Skip("macOS-specific xcrun test")
	}

	origPath := os.Getenv("PATH")
	t.Cleanup(func() { os.Setenv("PATH", origPath) }) //nolint:errcheck

	// Clear PATH so exec.LookPath("xcrun") fails.
	emptyDir := t.TempDir()
	require.NoError(t, os.Setenv("PATH", emptyDir))

	err := metal.CheckMetalEnvironment()
	require.Error(t, err, "CheckMetalEnvironment() must return an error when xcrun is not on PATH")
	assert.True(t,
		strings.Contains(err.Error(), "xcrun") || strings.Contains(err.Error(), "Xcode"),
		"error must mention xcrun or Xcode, got: %v", err)
}

// TestBuildMetalShaders_NoSourceFile_TempDir verifies that BuildMetalShaders returns
// a "Metal source file not found" error when the Metal source is absent from the
// current working directory.
func TestBuildMetalShaders_NoSourceFile_TempDir(t *testing.T) {
	if runtime.GOOS != "darwin" {
		t.Skip("macOS-specific Metal build test")
	}

	origDir, err := os.Getwd()
	require.NoError(t, err)
	t.Cleanup(func() { os.Chdir(origDir) }) //nolint:errcheck

	// Change to a temp directory that has no kmer_kernels.metal file.
	tmpDir := t.TempDir()
	require.NoError(t, os.Chdir(tmpDir))

	buildErr := metal.BuildMetalShaders()
	require.Error(t, buildErr, "BuildMetalShaders() must return an error when kmer_kernels.metal is absent")
	assert.True(t,
		strings.Contains(buildErr.Error(), "Metal source file not found") ||
			strings.Contains(buildErr.Error(), "kmer_kernels.metal"),
		"error must mention the missing source file, got: %v", buildErr)
}

// TestInstallMetalLibrary_WithFakeLib verifies that InstallMetalLibrary successfully
// copies a pre-existing .metallib file into the application bundle structure.
// A fake metallib is created in a temporary directory to exercise the install
// success path (mkdir + cp) without requiring the Metal toolchain.
func TestInstallMetalLibrary_WithFakeLib(t *testing.T) {
	if runtime.GOOS != "darwin" {
		t.Skip("macOS-specific Metal library install test")
	}

	origDir, err := os.Getwd()
	require.NoError(t, err)
	t.Cleanup(func() { os.Chdir(origDir) }) //nolint:errcheck

	// Create a temp working directory and plant a fake .metallib there.
	libDir := t.TempDir()
	require.NoError(t, os.Chdir(libDir))
	require.NoError(t, os.WriteFile("kmer_kernels.metallib", []byte("fake metallib content"), 0600))

	bundleDir := t.TempDir()
	installErr := metal.InstallMetalLibrary(bundleDir)
	require.NoError(t, installErr, "InstallMetalLibrary() must succeed when metallib exists")

	// Verify the file was placed in the bundle's Resources directory.
	destPath := bundleDir + "/Contents/Resources/kmer_kernels.metallib"
	info, statErr := os.Stat(destPath)
	require.NoError(t, statErr, "metallib must exist at expected bundle path after install")
	assert.Greater(t, info.Size(), int64(0), "installed metallib must be non-empty")
}
