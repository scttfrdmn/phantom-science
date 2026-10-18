// Package cuda counts k-mers on NVIDIA GPUs through a CUDA kernel
// (kmer_kernels.cu) called via a C++ bridge. The GPU implementation builds
// only with CGO, the cuda build tag and the CUDA Toolkit; other builds
// compile to stubs whose constructors return an error. This file holds
// helpers that check the CUDA environment and compile the kernels with nvcc.
package cuda

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
)

// BuildCudaKernels compiles kmer_kernels.cu in the working directory with nvcc.
func BuildCudaKernels() error {
	// Get the current package directory
	currentDir, err := os.Getwd()
	if err != nil {
		return fmt.Errorf("failed to get current directory: %w", err)
	}

	cudaSourcePath := filepath.Join(currentDir, "kmer_kernels.cu")
	ptxOutputPath := filepath.Join(currentDir, "kmer_kernels.ptx")

	// Check if CUDA source file exists
	if _, err := os.Stat(cudaSourcePath); os.IsNotExist(err) {
		return fmt.Errorf("CUDA source file not found: %s", cudaSourcePath)
	}

	// Compile CUDA kernels using nvcc
	fmt.Printf("Compiling CUDA kernels: %s\n", cudaSourcePath)

	args := []string{
		"--ptx",                                  // Generate PTX intermediate code
		"-O3",                                    // Optimize for performance
		"-use_fast_math",                         // Use fast math operations
		"-arch=sm_50",                            // Target compute capability 5.0+ (Maxwell and later)
		"-gencode", "arch=compute_50,code=sm_50", // Maxwell
		"-gencode", "arch=compute_60,code=sm_60", // Pascal
		"-gencode", "arch=compute_70,code=sm_70", // Volta
		"-gencode", "arch=compute_75,code=sm_75", // Turing
		"-gencode", "arch=compute_80,code=sm_80", // Ampere
		"-gencode", "arch=compute_86,code=sm_86", // Ampere (RTX 30xx)
		"-gencode", "arch=compute_89,code=sm_89", // Ada Lovelace
		"-gencode", "arch=compute_90,code=sm_90", // Hopper
		"-o", ptxOutputPath,
		cudaSourcePath,
	}

	cmd := exec.Command("nvcc", args...) // #nosec G204
	output, err := cmd.CombinedOutput()
	if err != nil {
		return fmt.Errorf("failed to compile CUDA kernels: %w\nOutput: %s", err, string(output))
	}

	fmt.Printf("CUDA kernels compiled successfully: %s\n", ptxOutputPath)
	return nil
}

// CheckCudaEnvironment verifies that the CUDA development environment is available
func CheckCudaEnvironment() error {
	// Check for nvcc compiler
	if _, err := exec.LookPath("nvcc"); err != nil {
		return fmt.Errorf("nvcc not found - please install CUDA Toolkit")
	}

	// Check CUDA version
	cmd := exec.Command("nvcc", "--version")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return fmt.Errorf("failed to get CUDA version: %w", err)
	}

	fmt.Printf("CUDA environment check passed:\n%s\n", string(output))

	// Check for CUDA runtime
	cmd = exec.Command("nvidia-smi")
	if err := cmd.Run(); err != nil {
		return fmt.Errorf("nvidia-smi not found - please install NVIDIA drivers")
	}

	return nil
}

// DetectCudaCapabilities detects available CUDA devices and their capabilities
func DetectCudaCapabilities() error {
	fmt.Printf("Detecting CUDA capabilities...\n")

	// Use nvidia-ml-py equivalent or deviceQuery if available
	cmd := exec.Command("nvidia-smi", "--query-gpu=name,memory.total,compute_cap", "--format=csv,noheader,nounits")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return fmt.Errorf("failed to query GPU information: %w", err)
	}

	lines := strings.Split(strings.TrimSpace(string(output)), "\n")
	for i, line := range lines {
		parts := strings.Split(line, ", ")
		if len(parts) >= 3 {
			fmt.Printf("GPU %d: %s, %s MB, Compute %s\n", i, parts[0], parts[1], parts[2])
		}
	}

	return nil
}

// GetCudaLibraryPaths returns the paths to CUDA libraries for linking
func GetCudaLibraryPaths() ([]string, error) {
	var paths []string

	// Common CUDA installation paths
	commonPaths := []string{
		"/usr/local/cuda/lib64",
		"/usr/local/cuda/lib",
		"/opt/cuda/lib64",
		"/opt/cuda/lib",
	}

	// Check environment variable
	if cudaPath := os.Getenv("CUDA_PATH"); cudaPath != "" {
		commonPaths = append([]string{
			filepath.Join(cudaPath, "lib64"),
			filepath.Join(cudaPath, "lib"),
		}, commonPaths...)
	}

	// Windows paths
	if runtime.GOOS == "windows" {
		if cudaPath := os.Getenv("CUDA_PATH"); cudaPath != "" {
			commonPaths = append(commonPaths, filepath.Join(cudaPath, "lib", "x64"))
		}
		commonPaths = append(commonPaths, "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v11.8\\lib\\x64")
		commonPaths = append(commonPaths, "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.0\\lib\\x64")
	}

	// Find existing paths
	for _, path := range commonPaths {
		if _, err := os.Stat(path); err == nil { // #nosec G703
			paths = append(paths, path)
		}
	}

	if len(paths) == 0 {
		return nil, fmt.Errorf("no CUDA library paths found")
	}

	return paths, nil
}

// GenerateCudaBuildFlags generates the appropriate build flags for CGO
func GenerateCudaBuildFlags() (string, string, error) {
	// Check CUDA availability
	if err := CheckCudaEnvironment(); err != nil {
		return "", "", fmt.Errorf("CUDA environment check failed: %w", err)
	}

	// Get library paths
	libPaths, err := GetCudaLibraryPaths()
	if err != nil {
		return "", "", fmt.Errorf("failed to find CUDA libraries: %w", err)
	}

	// Build CGO flags
	var cflags []string
	var ldflags []string

	// Include paths
	cudaIncludePaths := []string{
		"/usr/local/cuda/include",
		"/opt/cuda/include",
	}

	if cudaPath := os.Getenv("CUDA_PATH"); cudaPath != "" {
		cudaIncludePaths = append([]string{filepath.Join(cudaPath, "include")}, cudaIncludePaths...)
	}

	for _, path := range cudaIncludePaths {
		if _, err := os.Stat(path); err == nil { // #nosec G703
			cflags = append(cflags, "-I"+path)
			break // Use first found path
		}
	}

	// Library paths
	for _, path := range libPaths {
		ldflags = append(ldflags, "-L"+path)
	}

	// Required libraries
	requiredLibs := []string{"cudart", "cuda"}
	for _, lib := range requiredLibs {
		ldflags = append(ldflags, "-l"+lib)
	}

	// Platform-specific flags
	switch runtime.GOOS {
	case "linux":
		ldflags = append(ldflags, "-Wl,-rpath="+libPaths[0])
	case "darwin":
		// macOS specific flags if needed
	case "windows":
		// Windows specific flags if needed
	}

	return strings.Join(cflags, " "), strings.Join(ldflags, " "), nil
}

// WriteCudaBuildConfig writes a build configuration file for easier compilation
func WriteCudaBuildConfig(configPath string) error {
	cflags, ldflags, err := GenerateCudaBuildFlags()
	if err != nil {
		return fmt.Errorf("failed to generate build flags: %w", err)
	}

	config := fmt.Sprintf(`# CUDA Build Configuration
# Generated automatically - do not edit manually

# CGO flags for CUDA compilation
export CGO_CFLAGS="%s"
export CGO_LDFLAGS="%s"

# Build commands
build-cuda:
	CGO_ENABLED=1 go build -tags cuda ./...

test-cuda:
	CGO_ENABLED=1 go test -tags cuda ./...

benchmark-cuda:
	CGO_ENABLED=1 go test -tags cuda -bench=. ./...

# Kernel compilation
compile-kernels:
	nvcc --ptx -O3 -use_fast_math -arch=sm_50 -o kmer_kernels.ptx kmer_kernels.cu

clean-cuda:
	rm -f *.ptx *.cubin
`, cflags, ldflags)

	return os.WriteFile(configPath, []byte(config), 0600)
}

// ValidateCudaBuild performs a basic build validation
func ValidateCudaBuild() error {
	fmt.Printf("Validating CUDA build environment...\n")

	// Check CUDA environment
	if err := CheckCudaEnvironment(); err != nil {
		return err
	}

	// Check library paths
	_, err := GetCudaLibraryPaths()
	if err != nil {
		return err
	}

	// Try to detect capabilities
	if err := DetectCudaCapabilities(); err != nil {
		fmt.Printf("Warning: Could not detect GPU capabilities: %v\n", err)
	}

	fmt.Printf("✅ CUDA build environment validation passed\n")
	return nil
}
