//go:build !cgo

package cuda

import (
	"fmt"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// CudaKmerProcessor is a stub in builds without CGO; its constructor always
// returns an error.
type CudaKmerProcessor struct{}

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

// NewCudaKmerProcessor returns an error when CGO is disabled
func NewCudaKmerProcessor() (*CudaKmerProcessor, error) {
	return nil, fmt.Errorf("CUDA not available: CGO disabled")
}

// GetDeviceInfo returns empty device info
func (cp *CudaKmerProcessor) GetDeviceInfo() DeviceInfo {
	return DeviceInfo{}
}

// ExtractKmerFrequenciesGPU returns an error when CGO is disabled
func (cp *CudaKmerProcessor) ExtractKmerFrequenciesGPU(sequences [][]byte, k int) (*kmer.FrequencyTable, error) {
	return nil, fmt.Errorf("CUDA not available: CGO disabled")
}

// ExtractKmerFrequenciesGPUBatch returns an error when CGO is disabled
func (cp *CudaKmerProcessor) ExtractKmerFrequenciesGPUBatch(sequences [][]byte, k int, batchSize int) (*kmer.FrequencyTable, error) {
	return nil, fmt.Errorf("CUDA not available: CGO disabled")
}

// BenchmarkGPUPerformance returns an error when CGO is disabled
func (cp *CudaKmerProcessor) BenchmarkGPUPerformance(testSequences [][]byte, k int) (*BenchmarkResults, error) {
	return nil, fmt.Errorf("CUDA not available: CGO disabled")
}

// Close is a no-op when CGO is disabled
func (cp *CudaKmerProcessor) Close() error {
	return nil
}

// GetPerformanceStats returns empty stats when CGO is disabled
func (cp *CudaKmerProcessor) GetPerformanceStats() PerformanceStats {
	return PerformanceStats{}
}

// EstimateSpeedup returns 0 when CGO is disabled
func EstimateSpeedup(computeCapability, smCount int) float64 {
	return 0.0
}

// GenerateTestSequences returns count sequences of the given length cycling
// through ATGC (benchmark input only).
func GenerateTestSequences(count, length int) [][]byte {
	sequences := make([][]byte, count)
	bases := []byte{'A', 'T', 'G', 'C'}
	for i := 0; i < count; i++ {
		seq := make([]byte, length)
		for j := 0; j < length; j++ {
			seq[j] = bases[(i*length+j)%4]
		}
		sequences[i] = seq
	}
	return sequences
}
