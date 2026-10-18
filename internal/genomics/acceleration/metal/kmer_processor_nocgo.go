//go:build !cgo || !darwin

package metal

import (
	"fmt"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// MetalKmerProcessor provides GPU-accelerated k-mer processing using Apple Metal
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

// IsAvailable returns false when CGO is disabled
func IsAvailable() bool {
	return false
}

// NewMetalKmerProcessor returns an error when CGO is disabled
func NewMetalKmerProcessor() (*MetalKmerProcessor, error) {
	return nil, fmt.Errorf("Metal not available: CGO disabled")
}

// GetDeviceInfo returns empty device info
func (mp *MetalKmerProcessor) GetDeviceInfo() DeviceInfo {
	return DeviceInfo{}
}

// ExtractKmerFrequenciesGPU returns an error when CGO is disabled
func (mp *MetalKmerProcessor) ExtractKmerFrequenciesGPU(sequences [][]byte, k int) (*kmer.FrequencyTable, error) {
	return nil, fmt.Errorf("Metal not available: CGO disabled")
}

// ExtractKmerFrequenciesGPUBatch returns an error when CGO is disabled
func (mp *MetalKmerProcessor) ExtractKmerFrequenciesGPUBatch(sequences [][]byte, k int, batchSize int) (*kmer.FrequencyTable, error) {
	return nil, fmt.Errorf("Metal not available: CGO disabled")
}

// BenchmarkGPUPerformance returns an error when CGO is disabled
func (mp *MetalKmerProcessor) BenchmarkGPUPerformance(testSequences [][]byte, k int) (*BenchmarkResults, error) {
	return nil, fmt.Errorf("Metal not available: CGO disabled")
}

// Close is a no-op when CGO is disabled
func (mp *MetalKmerProcessor) Close() error {
	return nil
}

// GetPerformanceStats returns empty stats when CGO is disabled
func (mp *MetalKmerProcessor) GetPerformanceStats() PerformanceStats {
	return PerformanceStats{}
}

// EstimateSpeedup returns 0 when CGO is disabled
func EstimateSpeedup(maxThreads int) float64 {
	return 0.0
}

// GenerateTestSequences generates test DNA sequences
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
