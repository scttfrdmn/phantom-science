package authenticity

import (
	"context"
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewParallelHMMModel tests model creation.
func TestNewParallelHMMModel(t *testing.T) {
	model := NewParallelHMMModel(4)
	require.NotNil(t, model)

	stats := model.GetPerformanceStats()
	assert.NotNil(t, stats)
	assert.Greater(t, stats["cpu_cores"].(int), 0)
	assert.False(t, stats["trained"].(bool))

	// Cleanup the worker goroutines to avoid leak.
	model.Cleanup()
}

// TestParallelHMMModel_TrainFromSequence tests training on a sequence directly.
func TestParallelHMMModel_TrainFromSequence(t *testing.T) {
	model := NewParallelHMMModel(3)
	require.NotNil(t, model)
	defer model.Cleanup()

	seq := makeParallelTestSequence(500)
	err := model.TrainFromSequence(seq)
	require.NoError(t, err)

	stats := model.GetPerformanceStats()
	assert.Greater(t, stats["total_kmers"].(uint64), uint64(0))
}

// TestParallelHMMModel_TrainFromSequence_TooShort tests that short sequences are skipped.
func TestParallelHMMModel_TrainFromSequence_TooShort(t *testing.T) {
	model := NewParallelHMMModel(4)
	require.NotNil(t, model)
	defer model.Cleanup()

	err := model.TrainFromSequence([]byte("AT"))
	require.NoError(t, err)
}

// TestParallelHMMModel_GenerateSequence_Untrained tests generation from untrained model.
func TestParallelHMMModel_GenerateSequence_Untrained(t *testing.T) {
	model := NewParallelHMMModel(4)
	require.NotNil(t, model)
	defer model.Cleanup()

	// Untrained model falls back to generateSequenceUntrained which still works.
	seq, err := model.GenerateSequence(100)
	require.NoError(t, err)
	assert.Len(t, seq, 100)
	for _, base := range seq {
		assert.Contains(t, "ATGC", string(base))
	}
}

// TestParallelHMMModel_GenerateSequence_Trained tests generation after training.
func TestParallelHMMModel_GenerateSequence_Trained(t *testing.T) {
	model := NewParallelHMMModel(4)
	require.NotNil(t, model)
	defer model.Cleanup()

	seq := makeParallelTestSequence(1000)
	require.NoError(t, model.TrainFromSequence(seq))

	// Mark as trained so the trained path is exercised.
	model.normalizeMatrices()
	// Access via atomic store to mark trained.
	// Use internal state directly since we're in the same package.
	model.trained = 1

	generated, err := model.GenerateSequence(100)
	require.NoError(t, err)
	assert.Len(t, generated, 100)
	for _, base := range generated {
		assert.Contains(t, "ATGC", string(base))
	}
}

// TestParallelHMMModel_GenerateSequenceParallel tests the parallel generation path.
func TestParallelHMMModel_GenerateSequenceParallel(t *testing.T) {
	model := NewParallelHMMModel(4)
	require.NotNil(t, model)
	defer model.Cleanup()

	require.NoError(t, model.TrainFromSequence(makeParallelTestSequence(1000)))
	model.normalizeMatrices()
	model.trained = 1

	ctx := context.Background()
	generated, err := model.GenerateSequenceParallel(100, ctx)
	require.NoError(t, err)
	assert.Len(t, generated, 100)
}

// TestParallelHMMModel_GenerateSequenceParallel_Untrained tests that untrained returns an error.
func TestParallelHMMModel_GenerateSequenceParallel_Untrained(t *testing.T) {
	model := NewParallelHMMModel(4)
	require.NotNil(t, model)
	defer model.Cleanup()

	ctx := context.Background()
	_, err := model.GenerateSequenceParallel(100, ctx)
	assert.Error(t, err)
}

// TestParallelHMMModel_GetPerformanceStats tests the stats getter.
func TestParallelHMMModel_GetPerformanceStats(t *testing.T) {
	model := NewParallelHMMModel(4)
	require.NotNil(t, model)
	defer model.Cleanup()

	stats := model.GetPerformanceStats()
	// Required keys present.
	assert.Contains(t, stats, "cpu_cores")
	assert.Contains(t, stats, "gpu_enabled")
	assert.Contains(t, stats, "total_kmers")
	assert.Contains(t, stats, "current_gc")
	assert.Contains(t, stats, "current_complexity")
	assert.Contains(t, stats, "trained")
	assert.Contains(t, stats, "chunk_size")
}

// TestParallelHMMModel_MatrixOperations tests CPU matrix multiplication.
func TestParallelHMMModel_MatrixOperations(t *testing.T) {
	model := NewParallelHMMModel(2)
	require.NotNil(t, model)
	defer model.Cleanup()

	a := [][]float64{{1.0, 2.0}, {3.0, 4.0}}
	b := [][]float64{{5.0, 6.0}, {7.0, 8.0}}

	// Test CPU multiplication.
	resultCPU := model.multiplyMatricesCPU(a, b)
	assert.InDelta(t, 19.0, resultCPU[0][0], 1e-9) // 1*5 + 2*7
	assert.InDelta(t, 22.0, resultCPU[0][1], 1e-9) // 1*6 + 2*8
	assert.InDelta(t, 43.0, resultCPU[1][0], 1e-9) // 3*5 + 4*7
	assert.InDelta(t, 50.0, resultCPU[1][1], 1e-9) // 3*6 + 4*8

	// Test GPU path (delegates to CPU or optimized on this platform).
	resultGPU := model.multiplyMatricesGPU(a, b)
	assert.InDelta(t, 19.0, resultGPU[0][0], 1e-9)
}

// TestParallelHMMModel_VectorAddition tests CPU vector addition.
func TestParallelHMMModel_VectorAddition(t *testing.T) {
	model := NewParallelHMMModel(2)
	require.NotNil(t, model)
	defer model.Cleanup()

	a := []float64{1.0, 2.0, 3.0}
	b := []float64{4.0, 5.0, 6.0}

	resultCPU := model.vectorAddCPU(a, b)
	assert.InDelta(t, 5.0, resultCPU[0], 1e-9)
	assert.InDelta(t, 7.0, resultCPU[1], 1e-9)
	assert.InDelta(t, 9.0, resultCPU[2], 1e-9)

	resultGPU := model.vectorAddGPU(a, b)
	assert.InDelta(t, 5.0, resultGPU[0], 1e-9)
}

// TestParallelHMMModel_VectorAddition_Large tests optimized (parallel) vector add.
func TestParallelHMMModel_VectorAddition_Large(t *testing.T) {
	model := NewParallelHMMModel(2)
	require.NotNil(t, model)
	defer model.Cleanup()

	size := 2000
	a := make([]float64, size)
	b := make([]float64, size)
	for i := range a {
		a[i] = float64(i)
		b[i] = float64(i) * 2
	}

	result := model.vectorAddOptimized(a, b)
	assert.Len(t, result, size)
	for i, v := range result {
		assert.InDelta(t, float64(i)*3, v, 1e-9)
	}
}

// TestReadFASTAFile tests reading a FASTA file.
func TestReadFASTAFile(t *testing.T) {
	// Create a temporary FASTA file.
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "test.fa")

	content := ">seq1 description\nATCGATCGATCG\nATCG\n>seq2\nGCGCGCGC\n"
	require.NoError(t, os.WriteFile(fastaPath, []byte(content), 0644))

	sequences, err := ReadFASTAFile(fastaPath)
	require.NoError(t, err)
	require.Len(t, sequences, 2)

	// First sequence should be the concatenation of its two lines.
	assert.Equal(t, []byte("ATCGATCGATCGATCG"), sequences[0])
	assert.Equal(t, []byte("GCGCGCGC"), sequences[1])
}

// TestReadFASTAFile_NotFound tests reading a missing file.
func TestReadFASTAFile_NotFound(t *testing.T) {
	_, err := ReadFASTAFile("/nonexistent/path/file.fa")
	assert.Error(t, err)
}

// TestReadFASTAFile_EmptyFile tests reading an empty FASTA file.
func TestReadFASTAFile_EmptyFile(t *testing.T) {
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "empty.fa")
	require.NoError(t, os.WriteFile(fastaPath, []byte(""), 0644))

	sequences, err := ReadFASTAFile(fastaPath)
	require.NoError(t, err)
	assert.Empty(t, sequences)
}

// TestReadFASTAFile_SingleSequence tests a FASTA file with one sequence.
func TestReadFASTAFile_SingleSequence(t *testing.T) {
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "single.fa")
	content := ">seq1\nATCGATCG\n"
	require.NoError(t, os.WriteFile(fastaPath, []byte(content), 0644))

	sequences, err := ReadFASTAFile(fastaPath)
	require.NoError(t, err)
	require.Len(t, sequences, 1)
	assert.Equal(t, []byte("ATCGATCG"), sequences[0])
}

// TestDetectPerformanceCores tests the core detection helper.
func TestDetectPerformanceCores(t *testing.T) {
	cores := detectPerformanceCores()
	assert.Greater(t, cores, 0)
}

// TestParallelHMMModel_TrainFileParallel tests file-based parallel training.
func TestParallelHMMModel_TrainFileParallel(t *testing.T) {
	model := NewParallelHMMModel(3)
	require.NotNil(t, model)
	defer model.Cleanup()

	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "train.fa")

	// Write a FASTA file with enough sequence to train on.
	seq := makeParallelTestSequence(5000)
	content := ">chr1\n" + string(seq) + "\n"
	require.NoError(t, os.WriteFile(fastaPath, []byte(content), 0644))

	ctx := context.Background()
	err := model.TrainFromFileParallel(fastaPath, ctx)
	require.NoError(t, err)

	stats := model.GetPerformanceStats()
	assert.True(t, stats["trained"].(bool))
}

// TestParallelHMMModel_TrainFileParallel_NotFound tests missing FASTA file error.
func TestParallelHMMModel_TrainFileParallel_NotFound(t *testing.T) {
	model := NewParallelHMMModel(3)
	require.NotNil(t, model)
	defer model.Cleanup()

	ctx := context.Background()
	err := model.TrainFromFileParallel("/nonexistent/train.fa", ctx)
	assert.Error(t, err)
}

// makeParallelTestSequence generates a deterministic test sequence.
func makeParallelTestSequence(length int) []byte {
	bases := []byte{'A', 'T', 'G', 'C'}
	seq := make([]byte, length)
	for i := 0; i < length; i++ {
		seq[i] = bases[(i+i/5+i/17)%4]
	}
	return seq
}
