package kmer

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewGPUPackedProcessor tests the constructor with all backends.
func TestNewGPUPackedProcessor(t *testing.T) {
	backends := []GPUBackend{GPUBackendAuto, GPUBackendCPU, GPUBackendMetal}

	for _, backend := range backends {
		t.Run("", func(t *testing.T) {
			proc, err := NewGPUPackedProcessor(4, backend)
			require.NoError(t, err)
			require.NotNil(t, proc)

			info := proc.GetBackendInfo()
			assert.Equal(t, 4, info.K)
			assert.True(t, info.IsInitialized)
			assert.True(t, info.PackedOutput)

			require.NoError(t, proc.Close())
		})
	}
}

// TestNewGPUPackedProcessor_InvalidK tests constructor validation.
func TestNewGPUPackedProcessor_InvalidK(t *testing.T) {
	_, err := NewGPUPackedProcessor(0, GPUBackendCPU)
	assert.Error(t, err)

	_, err = NewGPUPackedProcessor(17, GPUBackendCPU)
	assert.Error(t, err)
}

// TestGPUPackedProcessor_ExtractKmerFrequenciesGPUPacked_CPU tests CPU-mode extraction.
func TestGPUPackedProcessor_ExtractKmerFrequenciesGPUPacked_CPU(t *testing.T) {
	proc, err := NewGPUPackedProcessor(3, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	seqs := [][]byte{
		[]byte("ATCGATCG"),
		[]byte("GCATGCAT"),
	}

	pft, err := proc.ExtractKmerFrequenciesGPUPacked(seqs)
	require.NoError(t, err)
	require.NotNil(t, pft)

	assert.Equal(t, 3, pft.K)
	assert.Greater(t, pft.TotalCount, 0)
}

// TestGPUPackedProcessor_ExtractKmerFrequenciesGPUPacked_Auto tests auto-detect extraction.
func TestGPUPackedProcessor_ExtractKmerFrequenciesGPUPacked_Auto(t *testing.T) {
	proc, err := NewGPUPackedProcessor(4, GPUBackendAuto)
	require.NoError(t, err)
	defer proc.Close()

	seqs := [][]byte{
		[]byte("ATCGATCGATCG"),
	}

	pft, err := proc.ExtractKmerFrequenciesGPUPacked(seqs)
	require.NoError(t, err)
	require.NotNil(t, pft)
	assert.Greater(t, pft.TotalCount, 0)
}

// TestGPUPackedProcessor_ExtractKmerFrequenciesGPUPacked_Uninitialized tests not-init guard.
func TestGPUPackedProcessor_ExtractKmerFrequenciesGPUPacked_Uninitialized(t *testing.T) {
	proc := &GPUPackedProcessor{
		k:             3,
		isInitialized: false,
	}

	_, err := proc.ExtractKmerFrequenciesGPUPacked([][]byte{[]byte("ATCG")})
	assert.Error(t, err)
}

// TestGPUPackedProcessor_ExtractKmerFrequenciesGPUPackedBatch tests batch extraction.
func TestGPUPackedProcessor_ExtractKmerFrequenciesGPUPackedBatch(t *testing.T) {
	proc, err := NewGPUPackedProcessor(3, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	seqs := make([][]byte, 10)
	for i := range seqs {
		seqs[i] = []byte("ATCGATCG")
	}

	t.Run("with_batch_size", func(t *testing.T) {
		pft, err := proc.ExtractKmerFrequenciesGPUPackedBatch(seqs, 3)
		require.NoError(t, err)
		require.NotNil(t, pft)
		assert.Greater(t, pft.TotalCount, 0)
	})

	t.Run("auto_batch_size", func(t *testing.T) {
		pft, err := proc.ExtractKmerFrequenciesGPUPackedBatch(seqs, 0)
		require.NoError(t, err)
		require.NotNil(t, pft)
	})

	t.Run("empty_sequences", func(t *testing.T) {
		_, err := proc.ExtractKmerFrequenciesGPUPackedBatch([][]byte{}, 10)
		assert.Error(t, err)
	})
}

// TestGPUPackedProcessor_GetBackendInfo tests the backend info struct.
func TestGPUPackedProcessor_GetBackendInfo(t *testing.T) {
	proc, err := NewGPUPackedProcessor(4, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	info := proc.GetBackendInfo()
	assert.Equal(t, GPUBackendCPU, info.Backend)
	assert.False(t, info.IsGPU)
	assert.True(t, info.IsInitialized)
	assert.Equal(t, 4, info.K)
	assert.True(t, info.PackedOutput)
}

// TestGPUPackedProcessor_GetBackendInfo_GPU tests GPU backend info.
func TestGPUPackedProcessor_GetBackendInfo_GPU(t *testing.T) {
	proc, err := NewGPUPackedProcessor(4, GPUBackendMetal)
	require.NoError(t, err)
	defer proc.Close()

	info := proc.GetBackendInfo()
	// On macOS this will be Metal (IsGPU=true); on other platforms Metal
	// initialization fails and the processor falls back to CPU.
	if proc.isMetalAvailable() {
		assert.Equal(t, GPUBackendMetal, info.Backend)
		assert.True(t, info.IsGPU)
	} else {
		assert.Equal(t, GPUBackendCPU, info.Backend)
		assert.False(t, info.IsGPU)
	}
}

// TestGPUPackedProcessor_BenchmarkPerformance tests the benchmark helper.
func TestGPUPackedProcessor_BenchmarkPerformance(t *testing.T) {
	proc, err := NewGPUPackedProcessor(3, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	seqs := [][]byte{
		[]byte("ATCGATCG"),
		[]byte("GCATGCAT"),
	}

	results, err := proc.BenchmarkGPUPackedPerformance(seqs)
	require.NoError(t, err)
	require.NotNil(t, results)

	assert.Equal(t, 3, results.K)
	assert.Equal(t, 2, results.TestSequences)
	// CPU result should always be populated.
	assert.NotNil(t, results.CPUResult)
	assert.Greater(t, results.CPUMemoryUsage, int64(0))
}

// TestGPUPackedProcessor_BenchmarkPerformance_Uninitialized tests not-init guard.
func TestGPUPackedProcessor_BenchmarkPerformance_Uninitialized(t *testing.T) {
	proc := &GPUPackedProcessor{
		k:             3,
		isInitialized: false,
	}

	_, err := proc.BenchmarkGPUPackedPerformance(nil)
	assert.Error(t, err)
}

// TestGPUPackedProcessor_GetPackedArrayStatistics tests packed array stats.
func TestGPUPackedProcessor_GetPackedArrayStatistics(t *testing.T) {
	proc, err := NewGPUPackedProcessor(2, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	pft, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)

	for _, kmer := range []string{"AT", "GC", "AT", "TA"} {
		require.NoError(t, pft.IncrementKmer(kmer))
	}

	stats := proc.GetPackedArrayStatistics(pft)
	assert.Equal(t, 2, stats.K)
	assert.Greater(t, stats.MemoryUsage, int64(0))
	assert.Greater(t, stats.EstimatedSpeedup, 0.0)
}

// TestGPUPackedProcessor_Close tests clean shutdown.
func TestGPUPackedProcessor_Close(t *testing.T) {
	proc, err := NewGPUPackedProcessor(3, GPUBackendCPU)
	require.NoError(t, err)

	err = proc.Close()
	assert.NoError(t, err)
	assert.False(t, proc.isInitialized)
}

// TestGPUPackedProcessor_MergePackedTables tests the internal merge.
func TestGPUPackedProcessor_MergePackedTables(t *testing.T) {
	proc, err := NewGPUPackedProcessor(2, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	dst, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)
	require.NoError(t, dst.IncrementKmer("AT"))

	src, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)
	require.NoError(t, src.IncrementKmer("AT"))
	require.NoError(t, src.IncrementKmer("GC"))

	err = proc.mergePackedTables(dst, src)
	require.NoError(t, err)

	freq, err := dst.GetFrequency("AT")
	require.NoError(t, err)
	assert.Equal(t, 2, freq)

	freq, err = dst.GetFrequency("GC")
	require.NoError(t, err)
	assert.Equal(t, 1, freq)
}

// TestGPUPackedProcessor_MergePackedTables_Mismatch tests k-value mismatch error.
func TestGPUPackedProcessor_MergePackedTables_Mismatch(t *testing.T) {
	proc, err := NewGPUPackedProcessor(2, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	dst, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)

	src, err := NewPackedFrequencyTable(3)
	require.NoError(t, err)

	err = proc.mergePackedTables(dst, src)
	assert.Error(t, err)
}

// TestGPUPackedProcessor_ComparePackedTables tests the comparison helper.
func TestGPUPackedProcessor_ComparePackedTables(t *testing.T) {
	proc, err := NewGPUPackedProcessor(2, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	table1, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)

	table2, err := NewPackedFrequencyTable(2)
	require.NoError(t, err)

	// Empty tables should match.
	assert.True(t, proc.comparePackedTables(table1, table2))

	// Add data to table1 only.
	require.NoError(t, table1.IncrementKmer("AT"))
	assert.False(t, proc.comparePackedTables(table1, table2))

	// Add same data to table2.
	require.NoError(t, table2.IncrementKmer("AT"))
	assert.True(t, proc.comparePackedTables(table1, table2))
}

// TestGPUPackedProcessor_GetEstimatedSpeedup tests speedup estimation.
func TestGPUPackedProcessor_GetEstimatedSpeedup(t *testing.T) {
	cpuProc, err := NewGPUPackedProcessor(3, GPUBackendCPU)
	require.NoError(t, err)
	defer cpuProc.Close()
	cpuSpeedup := cpuProc.getEstimatedSpeedup()
	assert.Greater(t, cpuSpeedup, 0.0)

	metalProc, err := NewGPUPackedProcessor(3, GPUBackendMetal)
	require.NoError(t, err)
	defer metalProc.Close()
	metalSpeedup := metalProc.getEstimatedSpeedup()
	assert.Greater(t, metalSpeedup, 0.0)

	// Metal should be faster than pure CPU.
	if metalProc.gpuBackend == GPUBackendMetal {
		assert.Greater(t, metalSpeedup, cpuSpeedup)
	}
}

// TestGPUPackedProcessor_CalculateOptimalBatchSize tests batch size calculation.
func TestGPUPackedProcessor_CalculateOptimalBatchSize(t *testing.T) {
	proc, err := NewGPUPackedProcessor(4, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	seqs := [][]byte{
		[]byte("ATCGATCGATCG"),
		[]byte("GCATGCATGCAT"),
	}

	batchSize := proc.calculateOptimalBatchSize(seqs)
	assert.GreaterOrEqual(t, batchSize, 100)
	assert.LessOrEqual(t, batchSize, 50000)
}

// TestGPUPackedProcessor_CalculateOptimalBatchSize_Empty tests empty input.
func TestGPUPackedProcessor_CalculateOptimalBatchSize_Empty(t *testing.T) {
	proc, err := NewGPUPackedProcessor(4, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	batchSize := proc.calculateOptimalBatchSize(nil)
	assert.Equal(t, 1000, batchSize)
}

// mockCudaProcessor is a minimal CudaProcessor for testing isCudaAvailable.
type mockCudaProcessor struct{ available bool }

func (m *mockCudaProcessor) ExtractKmerFrequenciesGPU(_ [][]byte, _ int) (*FrequencyTable, error) {
	return nil, nil
}
func (m *mockCudaProcessor) IsAvailable() bool { return m.available }
func (m *mockCudaProcessor) Close() error       { return nil }

// TestIsCudaAvailable_NilProcessor returns false when no processor is set.
func TestIsCudaAvailable_NilProcessor(t *testing.T) {
	proc, err := NewGPUPackedProcessor(4, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	proc.cudaProcessor = nil
	assert.False(t, proc.isCudaAvailable())
}

// TestIsCudaAvailable_WithProcessor delegates to processor.IsAvailable().
func TestIsCudaAvailable_WithProcessor(t *testing.T) {
	proc, err := NewGPUPackedProcessor(4, GPUBackendCPU)
	require.NoError(t, err)
	defer proc.Close()

	proc.cudaProcessor = &mockCudaProcessor{available: true}
	assert.True(t, proc.isCudaAvailable())

	proc.cudaProcessor = &mockCudaProcessor{available: false}
	assert.False(t, proc.isCudaAvailable())
}
