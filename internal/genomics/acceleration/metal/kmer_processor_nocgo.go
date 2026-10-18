//go:build !cgo || !darwin

package metal

import (
	"fmt"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// MetalKmerProcessor is a stub in this build (no CGO, or not macOS); its
// constructor always returns an error.
type MetalKmerProcessor struct{}

// DeviceInfo contains information about the GPU device
type DeviceInfo struct {
	Name          string
	Type          string
	MaxThreads    int
	MaxBufferSize int64
	IsAvailable   bool
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
	Throughput    float64
}

// IsAvailable returns false in this build
func IsAvailable() bool {
	return false
}

// NewMetalKmerProcessor returns an error in this build
func NewMetalKmerProcessor() (*MetalKmerProcessor, error) {
	return nil, fmt.Errorf("Metal not available: requires macOS and CGO")
}

// GetDeviceInfo returns empty device info
func (mp *MetalKmerProcessor) GetDeviceInfo() DeviceInfo {
	return DeviceInfo{}
}

// ExtractKmerFrequenciesGPU returns an error in this build
func (mp *MetalKmerProcessor) ExtractKmerFrequenciesGPU(sequences [][]byte, k int) (*kmer.FrequencyTable, error) {
	return nil, fmt.Errorf("Metal not available: requires macOS and CGO")
}

// ExtractKmerFrequenciesGPUBatch returns an error in this build
func (mp *MetalKmerProcessor) ExtractKmerFrequenciesGPUBatch(sequences [][]byte, k int, batchSize int) (*kmer.FrequencyTable, error) {
	return nil, fmt.Errorf("Metal not available: requires macOS and CGO")
}

// BenchmarkGPUPerformance returns an error in this build
func (mp *MetalKmerProcessor) BenchmarkGPUPerformance(testSequences [][]byte, k int) (*BenchmarkResults, error) {
	return nil, fmt.Errorf("Metal not available: requires macOS and CGO")
}

// Close is a no-op in this build
func (mp *MetalKmerProcessor) Close() error {
	return nil
}

// GetPerformanceStats returns empty stats in this build
func (mp *MetalKmerProcessor) GetPerformanceStats() PerformanceStats {
	return PerformanceStats{}
}

// EstimateSpeedup returns 0 in this build
func EstimateSpeedup(maxThreads int) float64 {
	return 0.0
}

// GenerateTestSequences returns count sequences of the given length, each the
// repeating pattern ATGC (benchmark input only).
func GenerateTestSequences(count, length int) [][]byte {
	sequences := make([][]byte, count)
	bases := []byte{'A', 'T', 'G', 'C'}
	for i := 0; i < count; i++ {
		seq := make([]byte, length)
		for j := 0; j < length; j++ {
			seq[j] = bases[j%4]
		}
		sequences[i] = seq
	}
	return sequences
}
