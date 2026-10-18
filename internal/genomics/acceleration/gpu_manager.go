// Package acceleration selects a k-mer counting backend: Apple Metal
// (macOS with CGO), NVIDIA CUDA (CGO and the cuda build tag) or the CPU
// extractor in authenticity/kmer. GPUManager falls back to the CPU when no
// GPU backend initializes or a GPU call fails, and BenchmarkPerformance times
// each available backend against the CPU on caller-supplied sequences.
package acceleration

import (
	"fmt"
	"runtime"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/acceleration/cuda"
	"github.com/scttfrdmn/phantom-science/internal/genomics/acceleration/metal"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// GPUManager provides unified GPU acceleration across different platforms
type GPUManager struct {
	metalProcessor *metal.MetalKmerProcessor
	cudaProcessor  *cuda.CudaKmerProcessor
	cpuExtractor   *kmer.FrequencyExtractor
	preferredGPU   GPUBackend
	isInitialized  bool
}

// GPUBackend represents the available GPU acceleration backends
type GPUBackend int

const (
	GPUBackendNone GPUBackend = iota
	GPUBackendMetal
	GPUBackendCUDA
	GPUBackendAuto // Automatically select best available
)

func (gb GPUBackend) String() string {
	switch gb {
	case GPUBackendNone:
		return "CPU"
	case GPUBackendMetal:
		return "Apple Metal"
	case GPUBackendCUDA:
		return "NVIDIA CUDA"
	case GPUBackendAuto:
		return "Auto-select"
	default:
		return "Unknown"
	}
}

// GPUCapabilities contains information about available GPU acceleration
type GPUCapabilities struct {
	MetalAvailable    bool
	CudaAvailable     bool
	PreferredBackend  GPUBackend
	AvailableBackends []GPUBackend
	DeviceInfo        []DeviceInfo
}

// DeviceInfo contains information about a GPU device
type DeviceInfo struct {
	Backend      GPUBackend
	Name         string
	Type         string
	Capabilities map[string]interface{}
	IsPreferred  bool
}

// NewGPUManager creates a new GPU manager with automatic backend detection
func NewGPUManager() (*GPUManager, error) {
	return NewGPUManagerWithBackend(GPUBackendAuto)
}

// NewGPUManagerWithBackend creates a new GPU manager with specified backend preference
func NewGPUManagerWithBackend(preferredBackend GPUBackend) (*GPUManager, error) {
	manager := &GPUManager{
		preferredGPU: preferredBackend,
	}

	// Initialize CPU fallback
	manager.cpuExtractor = kmer.NewFrequencyExtractor(6) // Default k=6

	// Auto-detect best available GPU backend
	if preferredBackend == GPUBackendAuto {
		manager.preferredGPU = manager.detectBestGPU()
	}

	// Initialize GPU processors based on availability and preference
	if err := manager.initializeGPUProcessors(); err != nil {
		// GPU initialization failed, but we can still use CPU
		fmt.Printf("Warning: GPU initialization failed: %v\n", err)
		fmt.Printf("Falling back to CPU-only processing\n")
		manager.preferredGPU = GPUBackendNone
	}

	manager.isInitialized = true
	return manager, nil
}

// detectBestGPU automatically detects the best available GPU backend
func (gm *GPUManager) detectBestGPU() GPUBackend {
	// Check platform-specific GPU availability
	switch runtime.GOOS {
	case "darwin":
		// macOS: prefer Metal, fallback to CUDA if available
		if metal.IsAvailable() {
			return GPUBackendMetal
		}
		if cuda.IsAvailable() {
			return GPUBackendCUDA
		}
	default:
		// Linux/Windows: prefer CUDA, Metal not available
		if cuda.IsAvailable() {
			return GPUBackendCUDA
		}
	}

	return GPUBackendNone
}

// initializeGPUProcessors initializes the appropriate GPU processors
func (gm *GPUManager) initializeGPUProcessors() error {
	switch gm.preferredGPU {
	case GPUBackendMetal:
		if !metal.IsAvailable() {
			return fmt.Errorf("Metal GPU acceleration not available")
		}
		processor, err := metal.NewMetalKmerProcessor()
		if err != nil {
			return fmt.Errorf("failed to initialize Metal processor: %w", err)
		}
		gm.metalProcessor = processor

	case GPUBackendCUDA:
		if !cuda.IsAvailable() {
			return fmt.Errorf("CUDA GPU acceleration not available")
		}
		processor, err := cuda.NewCudaKmerProcessor()
		if err != nil {
			return fmt.Errorf("failed to initialize CUDA processor: %w", err)
		}
		gm.cudaProcessor = processor

	case GPUBackendNone:
		// CPU-only mode, no GPU initialization needed
		return nil

	default:
		return fmt.Errorf("unsupported GPU backend: %v", gm.preferredGPU)
	}

	return nil
}

// GetCapabilities returns information about available GPU acceleration
func (gm *GPUManager) GetCapabilities() GPUCapabilities {
	capabilities := GPUCapabilities{
		MetalAvailable:    metal.IsAvailable(),
		CudaAvailable:     cuda.IsAvailable(),
		PreferredBackend:  gm.preferredGPU,
		AvailableBackends: []GPUBackend{},
		DeviceInfo:        []DeviceInfo{},
	}

	// Add available backends
	if capabilities.MetalAvailable {
		capabilities.AvailableBackends = append(capabilities.AvailableBackends, GPUBackendMetal)
	}
	if capabilities.CudaAvailable {
		capabilities.AvailableBackends = append(capabilities.AvailableBackends, GPUBackendCUDA)
	}
	capabilities.AvailableBackends = append(capabilities.AvailableBackends, GPUBackendNone)

	// Add device information
	if gm.metalProcessor != nil {
		metalInfo := gm.metalProcessor.GetDeviceInfo()
		capabilities.DeviceInfo = append(capabilities.DeviceInfo, DeviceInfo{
			Backend: GPUBackendMetal,
			Name:    metalInfo.Name,
			Type:    metalInfo.Type,
			Capabilities: map[string]interface{}{
				"maxThreads":    metalInfo.MaxThreads,
				"maxBufferSize": metalInfo.MaxBufferSize,
			},
			IsPreferred: gm.preferredGPU == GPUBackendMetal,
		})
	}

	if gm.cudaProcessor != nil {
		cudaInfo := gm.cudaProcessor.GetDeviceInfo()
		capabilities.DeviceInfo = append(capabilities.DeviceInfo, DeviceInfo{
			Backend: GPUBackendCUDA,
			Name:    cudaInfo.Name,
			Type:    cudaInfo.Type,
			Capabilities: map[string]interface{}{
				"maxThreadsPerBlock":  cudaInfo.MaxThreadsPerBlock,
				"globalMemorySize":    cudaInfo.GlobalMemorySize,
				"multiprocessorCount": cudaInfo.MultiprocessorCount,
				"computeCapability":   cudaInfo.ComputeCapability,
			},
			IsPreferred: gm.preferredGPU == GPUBackendCUDA,
		})
	}

	// Always add CPU info
	capabilities.DeviceInfo = append(capabilities.DeviceInfo, DeviceInfo{
		Backend: GPUBackendNone,
		Name:    fmt.Sprintf("CPU (%s %s)", runtime.GOOS, runtime.GOARCH),
		Type:    "CPU",
		Capabilities: map[string]interface{}{
			"cores": runtime.NumCPU(),
		},
		IsPreferred: gm.preferredGPU == GPUBackendNone,
	})

	return capabilities
}

// ExtractKmerFrequencies extracts k-mer frequencies using the best available processor
func (gm *GPUManager) ExtractKmerFrequencies(sequences [][]byte, k int) (*kmer.FrequencyTable, error) {
	if !gm.isInitialized {
		return nil, fmt.Errorf("GPU manager not initialized")
	}

	// Try GPU acceleration first
	switch gm.preferredGPU {
	case GPUBackendMetal:
		if gm.metalProcessor != nil {
			result, err := gm.metalProcessor.ExtractKmerFrequenciesGPU(sequences, k)
			if err == nil {
				return result, nil
			}
			fmt.Printf("Metal GPU processing failed: %v, falling back to CPU\n", err)
		}

	case GPUBackendCUDA:
		if gm.cudaProcessor != nil {
			result, err := gm.cudaProcessor.ExtractKmerFrequenciesGPU(sequences, k)
			if err == nil {
				return result, nil
			}
			fmt.Printf("CUDA GPU processing failed: %v, falling back to CPU\n", err)
		}
	}

	// Fallback to CPU processing
	gm.cpuExtractor.K = k
	return gm.cpuExtractor.ExtractFromMultipleSequences(sequences)
}

// ExtractKmerFrequenciesBatch processes sequences in batches using GPU acceleration
func (gm *GPUManager) ExtractKmerFrequenciesBatch(sequences [][]byte, k int, batchSize int) (*kmer.FrequencyTable, error) {
	if !gm.isInitialized {
		return nil, fmt.Errorf("GPU manager not initialized")
	}

	// Try GPU batch processing first
	switch gm.preferredGPU {
	case GPUBackendMetal:
		if gm.metalProcessor != nil {
			result, err := gm.metalProcessor.ExtractKmerFrequenciesGPUBatch(sequences, k, batchSize)
			if err == nil {
				return result, nil
			}
			fmt.Printf("Metal GPU batch processing failed: %v, falling back to CPU\n", err)
		}

	case GPUBackendCUDA:
		if gm.cudaProcessor != nil {
			result, err := gm.cudaProcessor.ExtractKmerFrequenciesGPUBatch(sequences, k, batchSize)
			if err == nil {
				return result, nil
			}
			fmt.Printf("CUDA GPU batch processing failed: %v, falling back to CPU\n", err)
		}
	}

	// Fallback to CPU processing with manual batching
	if batchSize <= 0 {
		batchSize = 1000 // Default CPU batch size
	}

	var allTables []*kmer.FrequencyTable
	gm.cpuExtractor.K = k

	for i := 0; i < len(sequences); i += batchSize {
		end := i + batchSize
		if end > len(sequences) {
			end = len(sequences)
		}

		batch := sequences[i:end]
		table, err := gm.cpuExtractor.ExtractFromMultipleSequences(batch)
		if err != nil {
			return nil, fmt.Errorf("CPU batch processing failed: %w", err)
		}

		allTables = append(allTables, table)
	}

	if len(allTables) == 1 {
		return allTables[0], nil
	}

	return kmer.MergeFrequencyTables(allTables...)
}

// BenchmarkPerformance times k-mer counting on each available backend
// (CPU always, Metal and CUDA when initialized) and reports measured
// throughput and speedup over the CPU.
func (gm *GPUManager) BenchmarkPerformance(testSequences [][]byte, k int) (*PerformanceBenchmark, error) {
	if !gm.isInitialized {
		return nil, fmt.Errorf("GPU manager not initialized")
	}

	benchmark := &PerformanceBenchmark{
		TestSequences: len(testSequences),
		KmerSize:      k,
		Results:       make(map[GPUBackend]*BackendResult),
	}

	// Benchmark CPU (always available)
	if result, err := gm.benchmarkCPU(testSequences, k); err == nil {
		benchmark.Results[GPUBackendNone] = result
	}

	// Benchmark Metal if available
	if gm.metalProcessor != nil {
		if result, err := gm.benchmarkMetal(testSequences, k); err == nil {
			benchmark.Results[GPUBackendMetal] = result
		}
	}

	// Benchmark CUDA if available
	if gm.cudaProcessor != nil {
		if result, err := gm.benchmarkCUDA(testSequences, k); err == nil {
			benchmark.Results[GPUBackendCUDA] = result
		}
	}

	// Calculate relative performance
	gm.calculateRelativePerformance(benchmark)

	return benchmark, nil
}

// benchmarkCPU benchmarks CPU k-mer processing
func (gm *GPUManager) benchmarkCPU(sequences [][]byte, k int) (*BackendResult, error) {
	extractor := kmer.NewFrequencyExtractor(k)
	start := time.Now()
	result, err := extractor.ExtractFromMultipleSequences(sequences)
	duration := time.Since(start)

	if err != nil {
		return nil, err
	}

	return &BackendResult{
		Backend:    GPUBackendNone,
		Duration:   duration,
		TotalKmers: result.TotalCount,
		Throughput: gm.calculateThroughput(sequences, duration),
		Success:    true,
	}, nil
}

// benchmarkMetal benchmarks Apple Metal GPU processing
func (gm *GPUManager) benchmarkMetal(sequences [][]byte, k int) (*BackendResult, error) {
	metalResult, err := gm.metalProcessor.BenchmarkGPUPerformance(sequences, k)
	if err != nil {
		return nil, err
	}

	return &BackendResult{
		Backend:      GPUBackendMetal,
		Duration:     metalResult.GPUDuration,
		TotalKmers:   metalResult.GPUKmers,
		Throughput:   metalResult.Throughput,
		Success:      metalResult.ResultsMatch,
		SpeedupVsCPU: metalResult.Speedup,
		DeviceName:   metalResult.DeviceName,
	}, nil
}

// benchmarkCUDA benchmarks NVIDIA CUDA GPU processing
func (gm *GPUManager) benchmarkCUDA(sequences [][]byte, k int) (*BackendResult, error) {
	cudaResult, err := gm.cudaProcessor.BenchmarkGPUPerformance(sequences, k)
	if err != nil {
		return nil, err
	}

	return &BackendResult{
		Backend:      GPUBackendCUDA,
		Duration:     cudaResult.GPUDuration,
		TotalKmers:   cudaResult.GPUKmers,
		Throughput:   cudaResult.Throughput,
		Success:      cudaResult.ResultsMatch,
		SpeedupVsCPU: cudaResult.Speedup,
		DeviceName:   cudaResult.DeviceName,
	}, nil
}

// calculateThroughput calculates throughput in MB/s
func (gm *GPUManager) calculateThroughput(sequences [][]byte, duration time.Duration) float64 {
	totalBytes := 0
	for _, seq := range sequences {
		totalBytes += len(seq)
	}
	seconds := float64(duration.Nanoseconds()) / 1e9
	return float64(totalBytes) / (1024 * 1024) / seconds
}

// calculateRelativePerformance calculates relative performance metrics
func (gm *GPUManager) calculateRelativePerformance(benchmark *PerformanceBenchmark) {
	cpuResult := benchmark.Results[GPUBackendNone]
	if cpuResult == nil {
		return
	}

	// Calculate speedup for GPU backends
	for backend, result := range benchmark.Results {
		if backend != GPUBackendNone && result.Success {
			result.SpeedupVsCPU = float64(cpuResult.Duration.Nanoseconds()) / float64(result.Duration.Nanoseconds())
		}
	}

	// Find best performing backend
	var bestBackend GPUBackend
	var bestThroughput float64

	for backend, result := range benchmark.Results {
		if result.Success && result.Throughput > bestThroughput {
			bestThroughput = result.Throughput
			bestBackend = backend
		}
	}

	benchmark.BestBackend = bestBackend
	benchmark.BestThroughput = bestThroughput
}

// PerformanceBenchmark contains the per-backend results of BenchmarkPerformance.
type PerformanceBenchmark struct {
	TestSequences  int
	KmerSize       int
	Results        map[GPUBackend]*BackendResult
	BestBackend    GPUBackend
	BestThroughput float64
}

// BackendResult contains benchmark results for a specific backend
type BackendResult struct {
	Backend      GPUBackend
	Duration     time.Duration
	TotalKmers   int
	Throughput   float64 // MB/s
	Success      bool
	SpeedupVsCPU float64
	DeviceName   string
}

// Close releases all GPU resources
func (gm *GPUManager) Close() error {
	if gm.metalProcessor != nil {
		gm.metalProcessor.Close()
		gm.metalProcessor = nil
	}

	if gm.cudaProcessor != nil {
		gm.cudaProcessor.Close()
		gm.cudaProcessor = nil
	}

	gm.isInitialized = false
	return nil
}

// SwitchBackend switches to a different GPU backend
func (gm *GPUManager) SwitchBackend(backend GPUBackend) error {
	if backend == gm.preferredGPU {
		return nil // Already using this backend
	}

	// Clean up current GPU processor
	switch gm.preferredGPU {
	case GPUBackendMetal:
		if gm.metalProcessor != nil {
			gm.metalProcessor.Close()
			gm.metalProcessor = nil
		}
	case GPUBackendCUDA:
		if gm.cudaProcessor != nil {
			gm.cudaProcessor.Close()
			gm.cudaProcessor = nil
		}
	}

	// Switch to new backend
	gm.preferredGPU = backend
	return gm.initializeGPUProcessors()
}
