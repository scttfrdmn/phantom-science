package authenticity

import (
	"fmt"
	"sync"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewAdaptiveHMMModel tests the constructor with valid and invalid inputs.
func TestNewAdaptiveHMMModel(t *testing.T) {
	validKValues := []int{2, 3, 4, 5, 6, 7, 8}
	for _, k := range validKValues {
		t.Run(fmt.Sprintf("valid_k=%d", k), func(t *testing.T) {
			model := NewAdaptiveHMMModel(k)
			require.NotNil(t, model)
			assert.Equal(t, k, model.K())
			assert.False(t, model.IsTrained())
			assert.Equal(t, uint64(0), model.TotalKmers())
			assert.NotEmpty(t, model.StateSpace())
		})
	}

	invalidKValues := []int{0, 1, 9, 10, -1}
	for _, k := range invalidKValues {
		t.Run(fmt.Sprintf("invalid_k=%d", k), func(t *testing.T) {
			model := NewAdaptiveHMMModel(k)
			assert.Nil(t, model)
		})
	}
}

// TestAdaptiveHMMModel_StateSpace tests the state space has the correct size.
func TestAdaptiveHMMModel_StateSpace(t *testing.T) {
	tests := []struct {
		k              int
		expectedStates int
	}{
		{2, 4},   // 4^1
		{3, 16},  // 4^2
		{4, 64},  // 4^3
		{5, 256}, // 4^4
	}

	for _, tt := range tests {
		t.Run(fmt.Sprintf("k=%d", tt.k), func(t *testing.T) {
			model := NewAdaptiveHMMModel(tt.k)
			require.NotNil(t, model)
			states := model.StateSpace()
			assert.Len(t, states, tt.expectedStates)
			// Each state should be a (k-1)-mer of valid bases.
			for _, state := range states {
				assert.Len(t, state, tt.k-1)
				for _, base := range state {
					assert.Contains(t, "ATGC", string(base))
				}
			}
		})
	}
}

// TestAdaptiveHMMModel_TrainFromReferenceGenome tests training with a mock genome.
func TestAdaptiveHMMModel_TrainFromReferenceGenome(t *testing.T) {
	genome := newAdaptiveTestGenome(5000)

	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)

	err := model.TrainFromReferenceGenome(genome)
	require.NoError(t, err)

	assert.True(t, model.IsTrained())
	assert.Greater(t, model.TotalKmers(), uint64(0))
}

// TestAdaptiveHMMModel_TrainFromSequence tests the direct sequence training path.
func TestAdaptiveHMMModel_TrainFromSequence(t *testing.T) {
	model := NewAdaptiveHMMModel(3)
	require.NotNil(t, model)

	seq := makeAdaptiveTestSequence(500)
	err := model.TrainFromSequence(seq)
	require.NoError(t, err)
	assert.True(t, model.IsTrained())
	assert.Greater(t, model.TotalKmers(), uint64(0))
}

// TestAdaptiveHMMModel_TrainFromSequence_TooShort tests short sequence skip.
func TestAdaptiveHMMModel_TrainFromSequence_TooShort(t *testing.T) {
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)

	// Sequence shorter than k should not error, just do nothing.
	err := model.TrainFromSequence([]byte("AT"))
	require.NoError(t, err)
}

// TestAdaptiveHMMModel_ValidateModel_Untrained tests that an untrained model fails validation.
func TestAdaptiveHMMModel_ValidateModel_Untrained(t *testing.T) {
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)

	err := model.ValidateModel()
	assert.Error(t, err)
}

// TestAdaptiveHMMModel_ValidateModel_Trained tests that a trained model passes validation.
func TestAdaptiveHMMModel_ValidateModel_Trained(t *testing.T) {
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)

	err := model.TrainFromSequence(makeAdaptiveTestSequence(500))
	require.NoError(t, err)

	err = model.ValidateModel()
	assert.NoError(t, err)
}

// TestAdaptiveHMMModel_GenerateSequence_Untrained tests that generation fails when untrained.
func TestAdaptiveHMMModel_GenerateSequence_Untrained(t *testing.T) {
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)

	_, err := model.GenerateSequence(100)
	assert.Error(t, err)
}

// TestAdaptiveHMMModel_GenerateSequence_TooShort tests length < k error.
func TestAdaptiveHMMModel_GenerateSequence_TooShort(t *testing.T) {
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)

	require.NoError(t, model.TrainFromSequence(makeAdaptiveTestSequence(500)))

	_, err := model.GenerateSequence(2)
	assert.Error(t, err)
}

// TestAdaptiveHMMModel_GenerateSequence tests successful sequence generation.
func TestAdaptiveHMMModel_GenerateSequence(t *testing.T) {
	genome := newAdaptiveTestGenome(5000)
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)
	require.NoError(t, model.TrainFromReferenceGenome(genome))

	lengths := []int{50, 100, 500}
	for _, length := range lengths {
		t.Run(fmt.Sprintf("length=%d", length), func(t *testing.T) {
			result, err := model.GenerateSequence(length)
			require.NoError(t, err)
			require.NotNil(t, result)

			assert.Len(t, result.Sequence, length)
			assert.Len(t, result.Qualities, length)
			assert.Greater(t, result.GenerationTime.Nanoseconds(), int64(0))

			// All bases should be valid DNA.
			for _, base := range result.Sequence {
				assert.Contains(t, "ATGC", string(base),
					"generated sequence contains invalid base: %c", base)
			}
			// Quality scores should be in valid Phred range.
			for _, q := range result.Qualities {
				assert.LessOrEqual(t, q, byte(41))
			}
		})
	}
}

// TestAdaptiveHMMModel_GenerationStats tests that generation updates statistics.
func TestAdaptiveHMMModel_GenerationStats(t *testing.T) {
	genome := newAdaptiveTestGenome(5000)
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)
	require.NoError(t, model.TrainFromReferenceGenome(genome))

	_, err := model.GenerateSequence(200)
	require.NoError(t, err)

	stats := model.GetGenerationStats()
	assert.Equal(t, 1, stats.TotalSequences)
	assert.Equal(t, int64(200), stats.TotalBases)
	assert.Greater(t, stats.AverageComplexity, 0.0)
	assert.Greater(t, stats.AverageGCContent, 0.0)
	assert.GreaterOrEqual(t, stats.MaxHomopolymer, 1)
}

// TestAdaptiveHMMModel_SetSeed tests seed reproducibility.
func TestAdaptiveHMMModel_SetSeed(t *testing.T) {
	genome := newAdaptiveTestGenome(5000)

	model1 := NewAdaptiveHMMModel(3)
	model2 := NewAdaptiveHMMModel(3)
	require.NotNil(t, model1)
	require.NotNil(t, model2)

	require.NoError(t, model1.TrainFromReferenceGenome(genome))
	require.NoError(t, model2.TrainFromReferenceGenome(genome))

	model1.SetSeed(42)
	model2.SetSeed(42)

	res1, err := model1.GenerateSequence(100)
	require.NoError(t, err)
	res2, err := model2.GenerateSequence(100)
	require.NoError(t, err)

	assert.Equal(t, res1.Sequence, res2.Sequence)
	assert.Equal(t, res1.Qualities, res2.Qualities)
}

// TestAdaptiveHMMModel_MultipleGenerations tests multiple generation calls accumulate stats.
func TestAdaptiveHMMModel_MultipleGenerations(t *testing.T) {
	genome := newAdaptiveTestGenome(3000)
	model := NewAdaptiveHMMModel(3)
	require.NotNil(t, model)
	require.NoError(t, model.TrainFromReferenceGenome(genome))

	const n = 5
	for i := 0; i < n; i++ {
		result, err := model.GenerateSequence(64)
		require.NoError(t, err)
		assert.Len(t, result.Sequence, 64)
	}

	stats := model.GetGenerationStats()
	assert.Equal(t, n, stats.TotalSequences)
}

// TestAdaptiveHMMModel_ConcurrentGeneration verifies that GenerateSequence is
// safe to call from multiple goroutines simultaneously (regression for the
// data race where generateAdaptiveSequence held only RLock while writing
// m.complexityMemory, and GenerateSequence read m.trained without any lock).
func TestAdaptiveHMMModel_ConcurrentGeneration(t *testing.T) {
	seq := makeAdaptiveTestSequence(2000)
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)
	require.NoError(t, model.TrainFromSequence(seq))

	const goroutines = 8
	const callsEach = 10
	errs := make(chan error, goroutines*callsEach)
	var wg sync.WaitGroup
	for g := 0; g < goroutines; g++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := 0; i < callsEach; i++ {
				result, err := model.GenerateSequence(100)
				if err != nil {
					errs <- err
					return
				}
				if len(result.Sequence) != 100 {
					errs <- fmt.Errorf("expected length 100, got %d", len(result.Sequence))
					return
				}
			}
		}()
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		require.NoError(t, err)
	}

	stats := model.GetGenerationStats()
	assert.Equal(t, goroutines*callsEach, stats.TotalSequences)
}

// TestAdaptiveHMMModel_TrainWithNBases tests that N-containing k-mers are skipped.
func TestAdaptiveHMMModel_TrainWithNBases(t *testing.T) {
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)

	// Sequence with N bases mixed in.
	seq := []byte("ATCGATCGNNATCGATCGATCGATCGNNNATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	err := model.TrainFromSequence(seq)
	require.NoError(t, err)
	assert.True(t, model.IsTrained())
}

// TestAdaptiveHMMModel_KmerFrequenciesTracked tests that k-mer frequencies are tracked.
func TestAdaptiveHMMModel_KmerFrequenciesTracked(t *testing.T) {
	genome := newAdaptiveTestGenome(5000)
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)
	require.NoError(t, model.TrainFromReferenceGenome(genome))

	result, err := model.GenerateSequence(200)
	require.NoError(t, err)
	assert.NotEmpty(t, result.KmerFrequencies)
}

// TestAdaptiveHMMModel_DifferentKValues tests generation for different k sizes.
func TestAdaptiveHMMModel_DifferentKValues(t *testing.T) {
	for _, k := range []int{2, 3, 4, 5, 6} {
		t.Run(fmt.Sprintf("k=%d", k), func(t *testing.T) {
			model := NewAdaptiveHMMModel(k)
			require.NotNil(t, model)
			require.NoError(t, model.TrainFromSequence(makeAdaptiveTestSequence(1000)))

			result, err := model.GenerateSequence(50)
			require.NoError(t, err)
			assert.Len(t, result.Sequence, 50)
		})
	}
}

// TestAdaptiveHMMModel_TrainMultipleChromosomes tests training over multiple chromosomes.
func TestAdaptiveHMMModel_TrainMultipleChromosomes(t *testing.T) {
	genome := newAdaptiveMultiChromGenome(3)
	model := NewAdaptiveHMMModel(4)
	require.NotNil(t, model)

	err := model.TrainFromReferenceGenome(genome)
	require.NoError(t, err)
	assert.True(t, model.IsTrained())
	assert.Greater(t, model.TotalKmers(), uint64(0))
}

// --- local test helpers ---

func newAdaptiveTestGenome(length int) *MockReferenceGenome {
	seq := makeAdaptiveTestSequence(length)
	return &MockReferenceGenome{
		chromosomes: map[string]*MockChromosome{
			"chr1": {Name: "chr1", Sequence: seq, Length: len(seq)},
		},
	}
}

func newAdaptiveMultiChromGenome(numChrom int) *MockReferenceGenome {
	chroms := make(map[string]*MockChromosome, numChrom)
	for i := 1; i <= numChrom; i++ {
		name := fmt.Sprintf("chr%d", i)
		seq := makeAdaptiveTestSequence(2000 + i*100)
		chroms[name] = &MockChromosome{Name: name, Sequence: seq, Length: len(seq)}
	}
	return &MockReferenceGenome{chromosomes: chroms}
}

// makeAdaptiveTestSequence generates a deterministic DNA test sequence.
func makeAdaptiveTestSequence(length int) []byte {
	bases := []byte{'A', 'T', 'G', 'C'}
	seq := make([]byte, length)
	for i := 0; i < length; i++ {
		seq[i] = bases[(i+i/7+i/31)%4]
	}
	return seq
}
