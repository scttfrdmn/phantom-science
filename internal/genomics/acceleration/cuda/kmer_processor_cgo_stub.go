//go:build cgo && !cuda

// kmer_processor_cgo_stub.go is compiled when CGO is enabled but the 'cuda'
// build tag is absent (i.e., the CUDA Toolkit is not installed or not
// requested).  It has two jobs:
//
//  1. Import "C" with no LDFLAGS — this makes the package a CGO package so
//     that cuda_bridge.cpp (which compiles to safe stubs via the
//     #ifdef __has_include(<cuda_runtime.h>) guard) is accepted by the build
//     system instead of producing "C++ source files not allowed" errors.
//
//  2. Provide the same API surface as kmer_processor_nocgo.go so that callers
//     compile regardless of whether CGO is enabled.
//
// On machines with the CUDA SDK installed, build with -tags cuda to get the
// real GPU-accelerated implementation in kmer_processor.go.

package cuda

import (
	"fmt"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

/*
// Empty CGO preamble — no LDFLAGS so the linker does not try to find libcuda.
#include "cuda_bridge.h"
*/
import "C"

// Silence the "C imported and not used" warning by referencing a type.
var _ = C.isCudaAvailable

// CudaKmerProcessor provides GPU-accelerated k-mer processing using NVIDIA CUDA.
// This build is a stub: CUDA is not available without the cuda build tag.
type CudaKmerProcessor struct{}

// DeviceInfo contains information about the CUDA GPU device.
type DeviceInfo struct {
	Name                string
	Type                string
	MaxThreadsPerBlock  int
	GlobalMemorySize    int64
	MultiprocessorCount int
	ComputeCapability   int
	IsAvailable         bool
}

// PerformanceStats contains performance information about CUDA GPU processing.
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

// BenchmarkResults contains the results of a GPU performance benchmark.
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

// IsAvailable returns false when the CUDA Toolkit is absent.
func IsAvailable() bool {
	return bool(C.isCudaAvailable())
}

// NewCudaKmerProcessor returns an error when CUDA is not available.
func NewCudaKmerProcessor() (*CudaKmerProcessor, error) {
	if !IsAvailable() {
		return nil, fmt.Errorf("CUDA not available: no CUDA device found (build with -tags cuda for full CUDA support)")
	}
	return nil, fmt.Errorf("CUDA not available: compiled without cuda build tag")
}

// GetDeviceInfo returns empty device info when CUDA is unavailable.
func (cp *CudaKmerProcessor) GetDeviceInfo() DeviceInfo {
	return DeviceInfo{}
}

// ExtractKmerFrequenciesGPU returns an error when CUDA is unavailable.
func (cp *CudaKmerProcessor) ExtractKmerFrequenciesGPU(sequences [][]byte, k int) (*kmer.FrequencyTable, error) {
	return nil, fmt.Errorf("CUDA not available: build with -tags cuda")
}

// ExtractKmerFrequenciesGPUBatch returns an error when CUDA is unavailable.
func (cp *CudaKmerProcessor) ExtractKmerFrequenciesGPUBatch(sequences [][]byte, k int, batchSize int) (*kmer.FrequencyTable, error) {
	return nil, fmt.Errorf("CUDA not available: build with -tags cuda")
}

// BenchmarkGPUPerformance returns an error when CUDA is unavailable.
func (cp *CudaKmerProcessor) BenchmarkGPUPerformance(testSequences [][]byte, k int) (*BenchmarkResults, error) {
	return nil, fmt.Errorf("CUDA not available: build with -tags cuda")
}

// Close is a no-op when CUDA is unavailable.
func (cp *CudaKmerProcessor) Close() error {
	return nil
}

// GetPerformanceStats returns empty stats when CUDA is unavailable.
func (cp *CudaKmerProcessor) GetPerformanceStats() PerformanceStats {
	return PerformanceStats{}
}

// EstimateSpeedup returns 0 when CUDA is unavailable.
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
