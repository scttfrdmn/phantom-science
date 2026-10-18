// Package metal counts k-mers on Apple GPUs through a Metal compute kernel
// (kmer_kernels.metal) called via an Objective-C bridge. The GPU
// implementation builds only on macOS with CGO; elsewhere the package
// compiles to stubs whose constructors return an error. This file holds
// helpers that compile and install the Metal shader library.
package metal

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
)

// BuildMetalShaders compiles kmer_kernels.metal in the working directory into a
// .metallib with xcrun.
func BuildMetalShaders() error {
	if runtime.GOOS != "darwin" {
		return fmt.Errorf("Metal shaders can only be built on macOS")
	}

	// Get the current package directory
	currentDir, err := os.Getwd()
	if err != nil {
		return fmt.Errorf("failed to get current directory: %w", err)
	}

	metalSourcePath := filepath.Join(currentDir, "kmer_kernels.metal")
	metallibPath := filepath.Join(currentDir, "kmer_kernels.metallib")

	// Check if Metal source file exists
	if _, err := os.Stat(metalSourcePath); os.IsNotExist(err) {
		return fmt.Errorf("Metal source file not found: %s", metalSourcePath)
	}

	// Compile Metal shaders using xcrun metal
	fmt.Printf("Compiling Metal shaders: %s\n", metalSourcePath)

	cmd := exec.Command("xcrun", "metal",
		"-o", metallibPath,
		metalSourcePath) // #nosec G204

	output, err := cmd.CombinedOutput()
	if err != nil {
		return fmt.Errorf("failed to compile Metal shaders: %w\nOutput: %s", err, string(output))
	}

	fmt.Printf("Metal shaders compiled successfully: %s\n", metallibPath)
	return nil
}

// CheckMetalEnvironment verifies that the Metal development environment is available
func CheckMetalEnvironment() error {
	if runtime.GOOS != "darwin" {
		return fmt.Errorf("Metal is only available on macOS")
	}

	// Check for xcrun metal compiler
	if _, err := exec.LookPath("xcrun"); err != nil {
		return fmt.Errorf("xcrun not found - please install Xcode command line tools")
	}

	// Test xcrun metal availability
	cmd := exec.Command("xcrun", "metal", "--version")
	if err := cmd.Run(); err != nil {
		return fmt.Errorf("Metal compiler not available - please install Xcode")
	}

	return nil
}

// InstallMetalLibrary copies the compiled Metal library to the application bundle
func InstallMetalLibrary(bundlePath string) error {
	if runtime.GOOS != "darwin" {
		return fmt.Errorf("Metal library installation only supported on macOS")
	}

	currentDir, err := os.Getwd()
	if err != nil {
		return fmt.Errorf("failed to get current directory: %w", err)
	}

	sourcePath := filepath.Join(currentDir, "kmer_kernels.metallib")
	if _, err := os.Stat(sourcePath); os.IsNotExist(err) {
		// Try to build if not found
		if err := BuildMetalShaders(); err != nil {
			return fmt.Errorf("Metal library not found and build failed: %w", err)
		}
	}

	// Create bundle resources directory if it doesn't exist
	resourcesDir := filepath.Join(bundlePath, "Contents", "Resources")
	if err := os.MkdirAll(resourcesDir, 0750); err != nil {
		return fmt.Errorf("failed to create resources directory: %w", err)
	}

	// Copy Metal library to bundle
	destPath := filepath.Join(resourcesDir, "kmer_kernels.metallib")
	cmd := exec.Command("cp", sourcePath, destPath) // #nosec G204
	if err := cmd.Run(); err != nil {
		return fmt.Errorf("failed to copy Metal library to bundle: %w", err)
	}

	fmt.Printf("Metal library installed to bundle: %s\n", destPath)
	return nil
}
