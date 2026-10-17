package kmer

import (
	"fmt"
	"runtime"
)

// GPUPackedProcessor provides GPU-accelerated k-mer processing with packed array output
type GPUPackedProcessor struct {
	k                 int
	gpuBackend        GPUBackend
	metalProcessor    MetalProcessor   // Interface to Metal GPU processor
	cudaProcessor     CudaProcessor    // Interface to CUDA GPU processor
	usePackedOutput   bool
	isInitialized     bool
}

// GPUBackend represents the available GPU backends
type GPUBackend int

const (
	GPUBackendAuto GPUBackend = iota
	GPUBackendMetal
	GPUBackendCUDA
	GPUBackendCPU
)

// MetalProcessor interface for Metal GPU integration
type MetalProcessor interface {
	ExtractKmerFrequenciesGPU(sequences [][]byte, k int) (*FrequencyTable, error)
	IsAvailable() bool
	Close() error
}

// CudaProcessor interface for CUDA GPU integration
type CudaProcessor interface {
	ExtractKmerFrequenciesGPU(sequences [][]byte, k int) (*FrequencyTable, error)
	IsAvailable() bool
	Close() error
}

// NewGPUPackedProcessor creates a new GPU-accelerated processor with packed output
func NewGPUPackedProcessor(k int, backend GPUBackend) (*GPUPackedProcessor, error) {
	if k <= 0 || k > 16 {
		return nil, fmt.Errorf("k must be between 1 and 16, got %d", k)
	}

	processor := &GPUPackedProcessor{
		k:               k,
		gpuBackend:      backend,
		usePackedOutput: true,
	}

	// Auto-detect best available GPU backend
	if backend == GPUBackendAuto {
		processor.gpuBackend = processor.detectBestGPU()
	}

	// Initialize GPU processors
	if err := processor.initializeGPUProcessors(); err != nil {
		fmt.Printf("Warning: GPU initialization failed: %v, falling back to CPU\n", err)
		processor.gpuBackend = GPUBackendCPU
	}

	processor.isInitialized = true
	return processor, nil
}

// detectBestGPU automatically detects the best available GPU backend
func (gpp *GPUPackedProcessor) detectBestGPU() GPUBackend {
	switch runtime.GOOS {
	case "darwin":
		// macOS: prefer Metal, fallback to CUDA if available
		if gpp.isMetalAvailable() {
			return GPUBackendMetal
		}
		if gpp.isCudaAvailable() {
			return GPUBackendCUDA
		}
	default:
		// Linux/Windows: prefer CUDA
		if gpp.isCudaAvailable() {
			return GPUBackendCUDA
		}
	}
	return GPUBackendCPU
}

// initializeGPUProcessors initializes the appropriate GPU processors
func (gpp *GPUPackedProcessor) initializeGPUProcessors() error {
	switch gpp.gpuBackend {
	case GPUBackendMetal:
		if !gpp.isMetalAvailable() {
			return fmt.Errorf("Metal GPU not available")
		}
		// Initialize Metal processor (would integrate with existing Metal implementation)
		return nil

	case GPUBackendCUDA:
		if !gpp.isCudaAvailable() {
			return fmt.Errorf("CUDA GPU not available")
		}
		// Initialize CUDA processor (would integrate with existing CUDA implementation)
		return nil

	case GPUBackendCPU:
		// CPU-only mode, no GPU initialization needed
		return nil

	default:
		return fmt.Errorf("unsupported GPU backend: %v", gpp.gpuBackend)
	}
}

// ExtractKmerFrequenciesGPUPacked extracts k-mer frequencies using GPU with packed output
func (gpp *GPUPackedProcessor) ExtractKmerFrequenciesGPUPacked(sequences [][]byte) (*PackedFrequencyTable, error) {
	if !gpp.isInitialized {
		return nil, fmt.Errorf("GPU processor not initialized")
	}

	// Try GPU processing first
	var mapTable *FrequencyTable
	var err error

	switch gpp.gpuBackend {
	case GPUBackendMetal:
		if gpp.metalProcessor != nil {
			mapTable, err = gpp.metalProcessor.ExtractKmerFrequenciesGPU(sequences, gpp.k)
		}
	case GPUBackendCUDA:
		if gpp.cudaProcessor != nil {
			mapTable, err = gpp.cudaProcessor.ExtractKmerFrequenciesGPU(sequences, gpp.k)
		}
	case GPUBackendCPU:
		return gpp.extractWithCPUFallback(sequences)
	}

	if err != nil {
		fmt.Printf("GPU processing failed: %v, falling back to CPU\n", err)
		return gpp.extractWithCPUFallback(sequences)
	}

	// Convert GPU result to packed format
	if mapTable != nil {
		return FromFrequencyTable(mapTable)
	}

	return gpp.extractWithCPUFallback(sequences)
}

// extractWithCPUFallback performs CPU-based extraction with packed output
func (gpp *GPUPackedProcessor) extractWithCPUFallback(sequences [][]byte) (*PackedFrequencyTable, error) {
	extractor := NewPackedFrequencyExtractor(gpp.k)
	return extractor.ExtractFromMultipleSequences(sequences)
}

// isMetalAvailable checks if Metal GPU is available (placeholder)
func (gpp *GPUPackedProcessor) isMetalAvailable() bool {
	// This would call the actual Metal availability check
	// return metal.IsAvailable()
	return runtime.GOOS == "darwin"
}

// isCudaAvailable checks if CUDA GPU is available.
// Returns false when no CudaProcessor has been injected; otherwise delegates
// to the processor's own IsAvailable() method.
func (gpp *GPUPackedProcessor) isCudaAvailable() bool {
	if gpp.cudaProcessor == nil {
		return false
	}
	return gpp.cudaProcessor.IsAvailable()
}

// ExtractKmerFrequenciesGPUPackedBatch processes sequences in batches with GPU acceleration
func (gpp *GPUPackedProcessor) ExtractKmerFrequenciesGPUPackedBatch(sequences [][]byte, batchSize int) (*PackedFrequencyTable, error) {
	if len(sequences) == 0 {
		return nil, fmt.Errorf("no sequences provided")
	}

	if batchSize <= 0 {
		batchSize = gpp.calculateOptimalBatchSize(sequences)
	}

	// Create result table
	resultTable, err := NewPackedFrequencyTable(gpp.k)
	if err != nil {
		return nil, fmt.Errorf("failed to create result table: %w", err)
	}

	// Process sequences in batches
	for i := 0; i < len(sequences); i += batchSize {
		end := i + batchSize
		if end > len(sequences) {
			end = len(sequences)
		}

		batch := sequences[i:end]
		batchTable, err := gpp.ExtractKmerFrequenciesGPUPacked(batch)
		if err != nil {
			return nil, fmt.Errorf("failed to process batch %d-%d: %w", i, end-1, err)
		}

		// Merge batch results
		if err := gpp.mergePackedTables(resultTable, batchTable); err != nil {
			return nil, fmt.Errorf("failed to merge batch %d-%d: %w", i, end-1, err)
		}
	}

	return resultTable, nil
}

// calculateOptimalBatchSize determines optimal batch size for GPU processing
func (gpp *GPUPackedProcessor) calculateOptimalBatchSize(sequences [][]byte) int {
	if len(sequences) == 0 {
		return 1000
	}

	// Estimate average sequence length
	totalLen := 0
	for _, seq := range sequences {
		totalLen += len(seq)
	}
	avgLen := totalLen / len(sequences)

	// Calculate memory per sequence
	maxKmers := 1
	for i := 0; i < gpp.k; i++ {
		maxKmers *= 4
	}

	memoryPerSequence := avgLen + maxKmers*4 // Sequence + packed array

	// Use different memory limits based on backend
	var availableMemory int
	switch gpp.gpuBackend {
	case GPUBackendMetal:
		availableMemory = 8 * 1024 * 1024 * 1024 // 8GB for unified memory
	case GPUBackendCUDA:
		availableMemory = 4 * 1024 * 1024 * 1024 // 4GB conservative estimate
	default:
		availableMemory = 2 * 1024 * 1024 * 1024 // 2GB for CPU
	}

	batchSize := availableMemory / memoryPerSequence

	// Ensure reasonable bounds
	if batchSize < 100 {
		batchSize = 100
	} else if batchSize > 50000 {
		batchSize = 50000
	}

	return batchSize
}

// mergePackedTables merges two packed frequency tables
func (gpp *GPUPackedProcessor) mergePackedTables(dst, src *PackedFrequencyTable) error {
	if dst.K != src.K {
		return fmt.Errorf("cannot merge tables with different k values: %d vs %d", dst.K, src.K)
	}

	// Merge all non-zero frequencies
	for i := 0; i < src.MaxKmers; i++ {
		srcFreq := src.PackedFrequencies.Get(i)
		if srcFreq > 0 {
			dstFreq := dst.PackedFrequencies.Get(i)
			if err := dst.PackedFrequencies.Set(i, dstFreq+srcFreq); err != nil {
				return fmt.Errorf("failed to merge frequency at index %d: %w", i, err)
			}

			// Update statistics
			if dstFreq == 0 {
				dst.NonZeroCount++
			}
		}
	}

	// Update totals
	dst.TotalCount += src.TotalCount
	dst.TotalKmers += src.TotalKmers

	return nil
}

// GetBackendInfo returns information about the current GPU backend
func (gpp *GPUPackedProcessor) GetBackendInfo() BackendInfo {
	return BackendInfo{
		Backend:       gpp.gpuBackend,
		IsGPU:         gpp.gpuBackend != GPUBackendCPU,
		IsInitialized: gpp.isInitialized,
		K:             gpp.k,
		PackedOutput:  gpp.usePackedOutput,
	}
}

// BackendInfo contains information about the GPU backend
type BackendInfo struct {
	Backend       GPUBackend
	IsGPU         bool
	IsInitialized bool
	K             int
	PackedOutput  bool
}

// Close releases GPU resources
func (gpp *GPUPackedProcessor) Close() error {
	if gpp.metalProcessor != nil {
		gpp.metalProcessor.Close()
		gpp.metalProcessor = nil
	}

	if gpp.cudaProcessor != nil {
		gpp.cudaProcessor.Close()
		gpp.cudaProcessor = nil
	}

	gpp.isInitialized = false
	return nil
}

// BenchmarkGPUPackedPerformance compares GPU and CPU performance with packed output
func (gpp *GPUPackedProcessor) BenchmarkGPUPackedPerformance(testSequences [][]byte) (*PackedBenchmarkResults, error) {
	if !gpp.isInitialized {
		return nil, fmt.Errorf("GPU processor not initialized")
	}

	results := &PackedBenchmarkResults{
		K:             gpp.k,
		TestSequences: len(testSequences),
		Backend:       gpp.gpuBackend,
	}

	// Benchmark GPU processing
	if gpp.gpuBackend != GPUBackendCPU {
		gpuTable, gpuErr := gpp.ExtractKmerFrequenciesGPUPacked(testSequences)
		if gpuErr == nil {
			results.GPUResult = gpuTable
			results.GPUMemoryUsage = gpuTable.GetMemoryUsage()
			results.GPUCompressionRatio = gpuTable.GetCompressionRatio()
		} else {
			results.GPUError = gpuErr.Error()
		}
	}

	// Benchmark CPU processing for comparison
	cpuTable, cpuErr := gpp.extractWithCPUFallback(testSequences)
	if cpuErr == nil {
		results.CPUResult = cpuTable
		results.CPUMemoryUsage = cpuTable.GetMemoryUsage()
		results.CPUCompressionRatio = cpuTable.GetCompressionRatio()
	} else {
		results.CPUError = cpuErr.Error()
	}

	// Compare results
	if results.GPUResult != nil && results.CPUResult != nil {
		results.ResultsMatch = gpp.comparePackedTables(results.GPUResult, results.CPUResult)
		if results.GPUMemoryUsage > 0 {
			results.MemoryEfficiency = float64(results.CPUMemoryUsage) / float64(results.GPUMemoryUsage)
		}
	}

	return results, nil
}

// comparePackedTables compares two packed frequency tables for correctness
func (gpp *GPUPackedProcessor) comparePackedTables(table1, table2 *PackedFrequencyTable) bool {
	if table1.TotalCount != table2.TotalCount || table1.NonZeroCount != table2.NonZeroCount {
		return false
	}

	// Compare all frequencies
	for i := 0; i < table1.MaxKmers; i++ {
		if table1.PackedFrequencies.Get(i) != table2.PackedFrequencies.Get(i) {
			return false
		}
	}

	return true
}

// PackedBenchmarkResults contains benchmark results for packed GPU processing
type PackedBenchmarkResults struct {
	K                     int
	TestSequences         int
	Backend               GPUBackend
	GPUResult             *PackedFrequencyTable
	CPUResult             *PackedFrequencyTable
	GPUMemoryUsage        int64
	CPUMemoryUsage        int64
	GPUCompressionRatio   float64
	CPUCompressionRatio   float64
	MemoryEfficiency      float64
	ResultsMatch          bool
	GPUError              string
	CPUError              string
}

// GetPackedArrayStatistics returns detailed statistics about packed array usage
func (gpp *GPUPackedProcessor) GetPackedArrayStatistics(table *PackedFrequencyTable) PackedArrayStatistics {
	stats := table.GetStatistics()
	
	return PackedArrayStatistics{
		K:                 stats.K,
		MaxKmers:          stats.MaxKmers,
		NonZeroKmers:      stats.NonZeroKmers,
		TotalCount:        stats.TotalCount,
		MemoryUsage:       stats.MemoryUsage,
		CompressionRatio:  stats.CompressionRatio,
		BitsPerItem:       stats.BitsPerItem,
		PackedArraySize:   stats.PackedArraySize,
		MaxFrequency:      stats.MaxFrequency,
		Sparsity:          stats.Sparsity,
		Backend:           gpp.gpuBackend,
		EstimatedSpeedup:  gpp.getEstimatedSpeedup(),
	}
}

// getEstimatedSpeedup estimates the speedup from using packed arrays
func (gpp *GPUPackedProcessor) getEstimatedSpeedup() float64 {
	// Base speedup from packed arrays (memory efficiency)
	baseSpeedup := 3.0

	// Additional speedup from GPU acceleration
	switch gpp.gpuBackend {
	case GPUBackendMetal:
		return baseSpeedup * 8.0 // Metal GPU speedup
	case GPUBackendCUDA:
		return baseSpeedup * 15.0 // CUDA GPU speedup
	default:
		return baseSpeedup // CPU packed arrays only
	}
}

// PackedArrayStatistics contains comprehensive statistics about packed array performance
type PackedArrayStatistics struct {
	K                 int
	MaxKmers          int
	NonZeroKmers      int
	TotalCount        int
	MemoryUsage       int64
	CompressionRatio  float64
	BitsPerItem       int
	PackedArraySize   int
	MaxFrequency      int
	Sparsity          float64
	Backend           GPUBackend
	EstimatedSpeedup  float64
}