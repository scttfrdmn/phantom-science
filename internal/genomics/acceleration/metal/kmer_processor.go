//go:build cgo && darwin

package metal

import (
	"fmt"
	"runtime"
	"time"
	"unsafe"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

/*
#cgo CFLAGS: -x objective-c
#cgo LDFLAGS: -framework Metal -framework Foundation

#include "metal_bridge.h"
*/
import "C"

// MetalKmerProcessor provides GPU-accelerated k-mer processing using Apple Metal
type MetalKmerProcessor struct {
	device         *C.MetalDevice
	maxThreads     int
	maxBufferSize  int64
	deviceName     string
	isInitialized  bool
}

// NewMetalKmerProcessor creates a new Metal-accelerated k-mer processor
func NewMetalKmerProcessor() (*MetalKmerProcessor, error) {
	// Check if Metal is available (macOS with Metal support)
	if runtime.GOOS != "darwin" {
		return nil, fmt.Errorf("Metal GPU acceleration is only available on macOS")
	}

	if !bool(C.isMetalAvailable()) {
		return nil, fmt.Errorf("Metal is not available on this device")
	}

	// Initialize Metal device
	device := C.initializeMetalDevice()
	if device == nil {
		return nil, fmt.Errorf("failed to initialize Metal device")
	}

	processor := &MetalKmerProcessor{
		device:        device,
		maxThreads:    int(C.getMetalMaxThreadsPerGroup(device)),
		maxBufferSize: int64(C.getMetalMaxBufferLength(device)),
		deviceName:    C.GoString(C.getMetalDeviceName(device)),
		isInitialized: true,
	}

	// Set finalizer to clean up Metal resources
	runtime.SetFinalizer(processor, (*MetalKmerProcessor).cleanup)

	return processor, nil
}

// IsAvailable checks if Metal GPU acceleration is available
func IsAvailable() bool {
	if runtime.GOOS != "darwin" {
		return false
	}
	return bool(C.isMetalAvailable())
}

// GetDeviceInfo returns information about the Metal device
func (mp *MetalKmerProcessor) GetDeviceInfo() DeviceInfo {
	return DeviceInfo{
		Name:          mp.deviceName,
		Type:          "Apple Metal GPU",
		MaxThreads:    mp.maxThreads,
		MaxBufferSize: mp.maxBufferSize,
		IsAvailable:   mp.isInitialized,
	}
}

// DeviceInfo contains information about the GPU device
type DeviceInfo struct {
	Name          string
	Type          string
	MaxThreads    int
	MaxBufferSize int64
	IsAvailable   bool
}

// ExtractKmerFrequenciesGPU extracts k-mer frequencies using Metal GPU acceleration
func (mp *MetalKmerProcessor) ExtractKmerFrequenciesGPU(sequences [][]byte, k int) (*kmer.FrequencyTable, error) {
	if !mp.isInitialized {
		return nil, fmt.Errorf("Metal device not initialized")
	}

	if len(sequences) == 0 {
		return nil, fmt.Errorf("no sequences provided")
	}

	if k <= 0 || k > 16 {
		return nil, fmt.Errorf("k must be between 1 and 16 for GPU processing")
	}

	// Convert sequences to GPU-friendly format
	gpuData, err := mp.prepareSequenceData(sequences)
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
	processedKmers := C.processKmersOnGPU(
		mp.device,
		gpuData.sequences,
		C.int(gpuData.numSequences),
		gpuData.sequenceLengths,
		C.int(k),
		&results[0],
		C.int(maxKmers),
	)

	if processedKmers < 0 {
		return nil, fmt.Errorf("GPU k-mer processing failed with error code %d", processedKmers)
	}

	// Convert GPU results back to FrequencyTable
	return mp.convertGPUResults(results, k, int(processedKmers)), nil
}

// GPUSequenceData holds sequence data prepared for GPU processing
type GPUSequenceData struct {
	sequences       *C.char
	sequenceLengths *C.int
	numSequences    int
	totalSize       int
}

// prepareSequenceData converts Go byte slices to C arrays for GPU processing
func (mp *MetalKmerProcessor) prepareSequenceData(sequences [][]byte) (*GPUSequenceData, error) {
	totalSize := 0
	for _, seq := range sequences {
		totalSize += len(seq) + 1 // +1 for null terminator
	}

	// Check buffer size limits
	if int64(totalSize) > mp.maxBufferSize {
		return nil, fmt.Errorf("sequence data size (%d bytes) exceeds GPU buffer limit (%d bytes)", 
			totalSize, mp.maxBufferSize)
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
		// Convert sequence to uppercase and copy
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

// convertGPUResults converts GPU integer results to FrequencyTable
func (mp *MetalKmerProcessor) convertGPUResults(results []C.int, k int, totalKmers int) *kmer.FrequencyTable {
	frequencies := make(map[string]int)
	
	// Generate all possible k-mers and map GPU results
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
func (mp *MetalKmerProcessor) ExtractKmerFrequenciesGPUBatch(sequences [][]byte, k int, batchSize int) (*kmer.FrequencyTable, error) {
	if !mp.isInitialized {
		return nil, fmt.Errorf("Metal device not initialized")
	}

	if batchSize <= 0 {
		batchSize = mp.calculateOptimalBatchSize(sequences, k)
	}

	var allTables []*kmer.FrequencyTable
	
	// Process sequences in batches
	for i := 0; i < len(sequences); i += batchSize {
		end := i + batchSize
		if end > len(sequences) {
			end = len(sequences)
		}
		
		batch := sequences[i:end]
		table, err := mp.ExtractKmerFrequenciesGPU(batch, k)
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
func (mp *MetalKmerProcessor) calculateOptimalBatchSize(sequences [][]byte, k int) int {
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
	
	// Use 80% of available GPU memory for batching
	availableMemory := int(float64(mp.maxBufferSize) * 0.8)
	batchSize := availableMemory / memoryPerSequence
	
	// Ensure reasonable bounds
	if batchSize < 100 {
		batchSize = 100
	} else if batchSize > 10000 {
		batchSize = 10000
	}

	return batchSize
}

// cleanup releases Metal resources
func (mp *MetalKmerProcessor) cleanup() {
	if mp.isInitialized && mp.device != nil {
		C.releaseMetalDevice(mp.device)
		mp.device = nil
		mp.isInitialized = false
	}
}

// Close manually releases Metal resources
func (mp *MetalKmerProcessor) Close() error {
	mp.cleanup()
	return nil
}

// GetPerformanceStats returns performance statistics for the GPU processor
func (mp *MetalKmerProcessor) GetPerformanceStats() PerformanceStats {
	return PerformanceStats{
		DeviceName:           mp.deviceName,
		MaxThreadsPerGroup:   mp.maxThreads,
		MaxBufferSize:        mp.maxBufferSize,
		EstimatedSpeedup:     EstimateSpeedup(mp.maxThreads),
		OptimalKRange:        [2]int{4, 12}, // Optimal k-mer sizes for GPU processing
		RecommendedBatchSize: 5000,
	}
}

// PerformanceStats contains performance information about GPU processing
type PerformanceStats struct {
	DeviceName           string
	MaxThreadsPerGroup   int
	MaxBufferSize        int64
	EstimatedSpeedup     float64
	OptimalKRange        [2]int
	RecommendedBatchSize int
}

// EstimateSpeedup estimates the performance improvement over CPU processing
func EstimateSpeedup(maxThreads int) float64 {
	// Conservative estimate based on parallel processing capability
	// Apple Silicon GPUs typically have 256-1024+ cores
	if maxThreads >= 1024 {
		return 8.0 // High-end M2/M3 chips
	} else if maxThreads >= 512 {
		return 6.0 // Mid-range M2 chips
	} else if maxThreads >= 256 {
		return 4.0 // M1 chips
	} else {
		return 2.0 // Older or limited Metal devices
	}
}

// BenchmarkGPUPerformance runs a performance benchmark of GPU k-mer processing
func (mp *MetalKmerProcessor) BenchmarkGPUPerformance(testSequences [][]byte, k int) (*BenchmarkResults, error) {
	if !mp.isInitialized {
		return nil, fmt.Errorf("Metal device not initialized")
	}

	// Generate test data if none provided
	if len(testSequences) == 0 {
		testSequences = GenerateTestSequences(1000, 1000) // 1000 sequences of 1000bp each
	}

	// Benchmark GPU processing
	startTime := time.Now()
	gpuResult, err := mp.ExtractKmerFrequenciesGPU(testSequences, k)
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
		DeviceName:     mp.deviceName,
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

// Additional helper functions would be implemented...

func GenerateTestSequences(count, length int) [][]byte {
	// Implementation would generate random DNA sequences for testing
	// This is a placeholder - full implementation would create realistic sequences
	sequences := make([][]byte, count)
	bases := []byte{'A', 'T', 'G', 'C'}
	
	for i := 0; i < count; i++ {
		seq := make([]byte, length)
		for j := 0; j < length; j++ {
			seq[j] = bases[j%4] // Simple pattern for testing
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