//go:build cgo && cuda

package cuda

import (
	"fmt"
	"runtime"
	"time"
	"unsafe"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

/*
#cgo CFLAGS: -I.
#cgo LDFLAGS: -lcudart -lcuda -L/usr/local/cuda/lib64

#include "cuda_bridge.h"
*/
import "C"

// CudaKmerProcessor provides GPU-accelerated k-mer processing using NVIDIA CUDA
type CudaKmerProcessor struct {
	device                *C.CUDADevice
	maxThreadsPerBlock    int
	globalMemorySize      int64
	deviceName            string
	computeCapability     int
	multiprocessorCount   int
	isInitialized         bool
}

// NewCudaKmerProcessor creates a new CUDA-accelerated k-mer processor
func NewCudaKmerProcessor() (*CudaKmerProcessor, error) {
	// Check if CUDA is available
	if !bool(C.isCudaAvailable()) {
		return nil, fmt.Errorf("CUDA is not available on this system")
	}

	// Initialize CUDA device
	device := C.initializeCudaDevice()
	if device == nil {
		return nil, fmt.Errorf("failed to initialize CUDA device")
	}

	processor := &CudaKmerProcessor{
		device:                device,
		maxThreadsPerBlock:    int(C.getCudaMaxThreadsPerBlock(device)),
		globalMemorySize:      int64(C.getCudaGlobalMemorySize(device)),
		deviceName:            C.GoString(C.getCudaDeviceName(device)),
		computeCapability:     int(C.getCudaComputeCapability(device)),
		multiprocessorCount:   int(C.getCudaMultiprocessorCount(device)),
		isInitialized:         true,
	}

	// Set finalizer to clean up CUDA resources
	runtime.SetFinalizer(processor, (*CudaKmerProcessor).cleanup)

	return processor, nil
}

// IsAvailable checks if CUDA GPU acceleration is available
func IsAvailable() bool {
	return bool(C.isCudaAvailable())
}

// GetDeviceInfo returns information about the CUDA device
func (cp *CudaKmerProcessor) GetDeviceInfo() DeviceInfo {
	return DeviceInfo{
		Name:                cp.deviceName,
		Type:                fmt.Sprintf("NVIDIA CUDA GPU (CC %d.%d)", cp.computeCapability/10, cp.computeCapability%10),
		MaxThreadsPerBlock:  cp.maxThreadsPerBlock,
		GlobalMemorySize:    cp.globalMemorySize,
		MultiprocessorCount: cp.multiprocessorCount,
		ComputeCapability:   cp.computeCapability,
		IsAvailable:         cp.isInitialized,
	}
}

// DeviceInfo contains information about the CUDA GPU device
type DeviceInfo struct {
	Name                string
	Type                string
	MaxThreadsPerBlock  int
	GlobalMemorySize    int64
	MultiprocessorCount int
	ComputeCapability   int
	IsAvailable         bool
}

// ExtractKmerFrequenciesGPU extracts k-mer frequencies using CUDA GPU acceleration
func (cp *CudaKmerProcessor) ExtractKmerFrequenciesGPU(sequences [][]byte, k int) (*kmer.FrequencyTable, error) {
	if !cp.isInitialized {
		return nil, fmt.Errorf("CUDA device not initialized")
	}

	if len(sequences) == 0 {
		return nil, fmt.Errorf("no sequences provided")
	}

	if k <= 0 || k > 16 {
		return nil, fmt.Errorf("k must be between 1 and 16 for GPU processing")
	}

	// Convert sequences to GPU-friendly format
	gpuData, err := cp.prepareSequenceData(sequences)
	if err != nil {
		return nil, fmt.Errorf("failed to prepare sequence data: %w", err)
	}
	defer gpuData.cleanup()

	// Calculate maximum possible k-mers (4^k)
	maxKmers := 1
	for i := 0; i < k; i++ {
		maxKmers *= 4
	}

	// Allocate result buffer for k-mer counts
	results := make([]C.int, maxKmers)

	// Process k-mers on GPU
	processedKmers := C.processKmersOnCuda(
		cp.device,
		gpuData.sequences,
		C.int(gpuData.numSequences),
		gpuData.sequenceLengths,
		C.int(k),
		&results[0],
		C.int(maxKmers),
	)

	if processedKmers < 0 {
		return nil, fmt.Errorf("CUDA k-mer processing failed with error code %d", processedKmers)
	}

	// Convert GPU results back to FrequencyTable
	return cp.convertCudaResults(results, k, int(processedKmers)), nil
}

// GPUSequenceData holds sequence data prepared for GPU processing
type GPUSequenceData struct {
	sequences       *C.char
	sequenceLengths *C.int
	numSequences    int
	totalSize       int
}

// prepareSequenceData converts Go byte slices to C arrays for GPU processing
func (cp *CudaKmerProcessor) prepareSequenceData(sequences [][]byte) (*GPUSequenceData, error) {
	totalSize := 0
	for _, seq := range sequences {
		totalSize += len(seq) + 1 // +1 for null terminator
	}

	// Check memory limits (use 80% of available GPU memory)
	if int64(totalSize) > cp.globalMemorySize/5*4 {
		return nil, fmt.Errorf("sequence data size (%d bytes) exceeds GPU memory limit (%d bytes)",
			totalSize, cp.globalMemorySize/5*4)
	}

	// Allocate continuous memory for all sequences
	cSequences := C.malloc(C.size_t(totalSize))
	if cSequences == nil {
		return nil, fmt.Errorf("failed to allocate memory for GPU sequence data")
	}

	// Allocate array for sequence lengths
	cLengths := (*C.int)(C.malloc(C.size_t(len(sequences) * int(unsafe.Sizeof(C.int(0))))))
	if cLengths == nil {
		C.free(cSequences)
		return nil, fmt.Errorf("failed to allocate memory for sequence lengths")
	}

	// Copy sequences to continuous memory buffer
	offset := 0
	lengthArray := (*[1 << 20]C.int)(unsafe.Pointer(cLengths))[:len(sequences):len(sequences)]

	for i, seq := range sequences {
		// Convert sequence to uppercase and validate
		upperSeq := make([]byte, len(seq))
		for j, base := range seq {
			switch base {
			case 'a', 'A':
				upperSeq[j] = 'A'
			case 't', 'T':
				upperSeq[j] = 'T'
			case 'g', 'G':
				upperSeq[j] = 'G'
			case 'c', 'C':
				upperSeq[j] = 'C'
			case 'n', 'N':
				upperSeq[j] = 'N'
			default:
				upperSeq[j] = 'N' // Replace invalid bases with N
			}
		}

		// Copy to C buffer
		C.memcpy(unsafe.Pointer(uintptr(cSequences)+uintptr(offset)),
			unsafe.Pointer(&upperSeq[0]), C.size_t(len(seq)))

		// Add null terminator
		*(*C.char)(unsafe.Pointer(uintptr(cSequences) + uintptr(offset) + uintptr(len(seq)))) = 0

		lengthArray[i] = C.int(len(seq))
		offset += len(seq) + 1
	}

	return &GPUSequenceData{
		sequences:       (*C.char)(cSequences),
		sequenceLengths: cLengths,
		numSequences:    len(sequences),
		totalSize:       totalSize,
	}, nil
}

// cleanup releases GPU sequence data memory
func (gsd *GPUSequenceData) cleanup() {
	if gsd.sequences != nil {
		C.free(unsafe.Pointer(gsd.sequences))
		gsd.sequences = nil
	}
	if gsd.sequenceLengths != nil {
		C.free(unsafe.Pointer(gsd.sequenceLengths))
		gsd.sequenceLengths = nil
	}
}

// convertCudaResults converts CUDA integer results to FrequencyTable
func (cp *CudaKmerProcessor) convertCudaResults(results []C.int, k int, totalKmers int) *kmer.FrequencyTable {
	frequencies := make(map[string]int)

	// Generate all possible k-mers and map CUDA results
	bases := []byte{'A', 'T', 'G', 'C'}
	kmerIndex := 0

	var generateKmers func(current []byte, remaining int)
	generateKmers = func(current []byte, remaining int) {
		if remaining == 0 {
			if int(results[kmerIndex]) > 0 {
				frequencies[string(current)] = int(results[kmerIndex])
			}
			kmerIndex++
			return
		}

		for _, base := range bases {
			newCurrent := make([]byte, len(current)+1)
			copy(newCurrent, current)
			newCurrent[len(current)] = base
			generateKmers(newCurrent, remaining-1)
		}
	}

	generateKmers([]byte{}, k)

	// Calculate total count
	totalCount := 0
	for _, count := range frequencies {
		totalCount += count
	}

	return &kmer.FrequencyTable{
		K:           k,
		Frequencies: frequencies,
		TotalCount:  totalCount,
		TotalKmers:  totalKmers,
	}
}

// ExtractKmerFrequenciesGPUBatch processes multiple sequence batches on GPU
func (cp *CudaKmerProcessor) ExtractKmerFrequenciesGPUBatch(sequences [][]byte, k int, batchSize int) (*kmer.FrequencyTable, error) {
	if !cp.isInitialized {
		return nil, fmt.Errorf("CUDA device not initialized")
	}

	if batchSize <= 0 {
		batchSize = cp.calculateOptimalBatchSize(sequences, k)
	}

	var allTables []*kmer.FrequencyTable

	// Process sequences in batches
	for i := 0; i < len(sequences); i += batchSize {
		end := i + batchSize
		if end > len(sequences) {
			end = len(sequences)
		}

		batch := sequences[i:end]
		table, err := cp.ExtractKmerFrequenciesGPU(batch, k)
		if err != nil {
			return nil, fmt.Errorf("failed to process batch %d-%d: %w", i, end-1, err)
		}

		allTables = append(allTables, table)
	}

	// Merge all frequency tables
	if len(allTables) == 1 {
		return allTables[0], nil
	}

	return kmer.MergeFrequencyTables(allTables...)
}

// calculateOptimalBatchSize determines the optimal batch size based on available GPU memory
func (cp *CudaKmerProcessor) calculateOptimalBatchSize(sequences [][]byte, k int) int {
	if len(sequences) == 0 {
		return 1000 // Default batch size
	}

	// Estimate average sequence length
	totalLength := 0
	for _, seq := range sequences {
		totalLength += len(seq)
	}
	avgLength := totalLength / len(sequences)

	// Calculate memory per sequence (sequence + results)
	memoryPerSequence := avgLength + (1 << (2 * k)) * 4 // 4 bytes per int result

	// Use 70% of available GPU memory for batching (more conservative than Metal)
	availableMemory := int(float64(cp.globalMemorySize) * 0.7)
	batchSize := availableMemory / memoryPerSequence

	// Ensure reasonable bounds
	if batchSize < 100 {
		batchSize = 100
	} else if batchSize > 50000 {
		batchSize = 50000
	}

	return batchSize
}

// cleanup releases CUDA resources
func (cp *CudaKmerProcessor) cleanup() {
	if cp.isInitialized && cp.device != nil {
		C.releaseCudaDevice(cp.device)
		cp.device = nil
		cp.isInitialized = false
	}
}

// Close manually releases CUDA resources
func (cp *CudaKmerProcessor) Close() error {
	cp.cleanup()
	return nil
}

// GetPerformanceStats returns performance statistics for the GPU processor
func (cp *CudaKmerProcessor) GetPerformanceStats() PerformanceStats {
	return PerformanceStats{
		DeviceName:           cp.deviceName,
		MaxThreadsPerBlock:   cp.maxThreadsPerBlock,
		GlobalMemorySize:     cp.globalMemorySize,
		MultiprocessorCount:  cp.multiprocessorCount,
		ComputeCapability:    cp.computeCapability,
		EstimatedSpeedup:     EstimateSpeedup(cp.computeCapability, cp.multiprocessorCount),
		OptimalKRange:        [2]int{6, 14}, // Optimal k-mer sizes for CUDA processing
		RecommendedBatchSize: 10000,
	}
}

// PerformanceStats contains performance information about CUDA GPU processing
type PerformanceStats struct {
	DeviceName           string
	MaxThreadsPerBlock   int
	GlobalMemorySize     int64
	MultiprocessorCount  int
	ComputeCapability    int
	EstimatedSpeedup     float64
	OptimalKRange        [2]int
	RecommendedBatchSize int
}

// EstimateSpeedup estimates the performance improvement over CPU processing
func EstimateSpeedup(computeCapability, smCount int) float64 {
	// Conservative estimate based on compute capability and SM count
	baseSpeedup := 1.0

	// Compute capability factor
	switch {
	case computeCapability >= 80: // Ampere (RTX 30xx, A100)
		baseSpeedup = 30.0
	case computeCapability >= 75: // Turing (RTX 20xx, T4)
		baseSpeedup = 25.0
	case computeCapability >= 70: // Volta (V100)
		baseSpeedup = 20.0
	case computeCapability >= 60: // Pascal (GTX 10xx, P100)
		baseSpeedup = 15.0
	case computeCapability >= 50: // Maxwell (GTX 9xx)
		baseSpeedup = 10.0
	default: // Older architectures
		baseSpeedup = 5.0
	}

	// SM count factor (more SMs = more parallelism)
	smFactor := 1.0
	switch {
	case smCount >= 108: // High-end cards (RTX 4090, A100)
		smFactor = 1.5
	case smCount >= 80: // High-end cards (RTX 3080, V100)
		smFactor = 1.3
	case smCount >= 60: // Mid-high cards (RTX 3070)
		smFactor = 1.2
	case smCount >= 40: // Mid-range cards
		smFactor = 1.1
	case smCount >= 20: // Entry cards
		smFactor = 1.0
	default: // Very low-end or old cards
		smFactor = 0.8
	}

	return baseSpeedup * smFactor
}

// BenchmarkGPUPerformance runs a performance benchmark of GPU k-mer processing
func (cp *CudaKmerProcessor) BenchmarkGPUPerformance(testSequences [][]byte, k int) (*BenchmarkResults, error) {
	if !cp.isInitialized {
		return nil, fmt.Errorf("CUDA device not initialized")
	}

	// Generate test data if none provided
	if len(testSequences) == 0 {
		testSequences = GenerateTestSequences(1000, 1000) // 1000 sequences of 1000bp each
	}

	// Benchmark GPU processing
	startTime := time.Now()
	gpuResult, err := cp.ExtractKmerFrequenciesGPU(testSequences, k)
	gpuDuration := time.Since(startTime)

	if err != nil {
		return nil, fmt.Errorf("GPU benchmark failed: %w", err)
	}

	// Benchmark CPU processing for comparison
	extractor := kmer.NewFrequencyExtractor(k)
	startTime = time.Now()
	cpuResult, err := extractor.ExtractFromMultipleSequences(testSequences)
	cpuDuration := time.Since(startTime)

	if err != nil {
		return nil, fmt.Errorf("CPU benchmark failed: %w", err)
	}

	// Calculate speedup
	speedup := float64(cpuDuration.Nanoseconds()) / float64(gpuDuration.Nanoseconds())

	return &BenchmarkResults{
		DeviceName:     cp.deviceName,
		TestSequences:  len(testSequences),
		KmerSize:       k,
		GPUDuration:    gpuDuration,
		CPUDuration:    cpuDuration,
		Speedup:        speedup,
		GPUKmers:       gpuResult.TotalCount,
		CPUKmers:       cpuResult.TotalCount,
		ResultsMatch:   compareTables(gpuResult, cpuResult),
		Throughput:     calculateThroughput(testSequences, gpuDuration),
	}, nil
}

// BenchmarkResults contains the results of a GPU performance benchmark
type BenchmarkResults struct {
	DeviceName    string
	TestSequences int
	KmerSize      int
	GPUDuration   time.Duration
	CPUDuration   time.Duration
	Speedup       float64
	GPUKmers      int
	CPUKmers      int
	ResultsMatch  bool
	Throughput    float64 // MB/s
}

// Helper functions for benchmarking
func GenerateTestSequences(count, length int) [][]byte {
	sequences := make([][]byte, count)
	bases := []byte{'A', 'T', 'G', 'C'}

	for i := 0; i < count; i++ {
		seq := make([]byte, length)
		for j := 0; j < length; j++ {
			seq[j] = bases[(i*length+j)%4] // Deterministic but varied pattern
		}
		sequences[i] = seq
	}

	return sequences
}

func compareTables(gpu, cpu *kmer.FrequencyTable) bool {
	if gpu.TotalCount != cpu.TotalCount {
		return false
	}

	// Compare frequency maps
	for kmer, count := range gpu.Frequencies {
		if cpu.Frequencies[kmer] != count {
			return false
		}
	}

	return true
}

func calculateThroughput(sequences [][]byte, duration time.Duration) float64 {
	totalBytes := 0
	for _, seq := range sequences {
		totalBytes += len(seq)
	}

	seconds := float64(duration.Nanoseconds()) / 1e9
	return float64(totalBytes) / (1024 * 1024) / seconds // MB/s
}