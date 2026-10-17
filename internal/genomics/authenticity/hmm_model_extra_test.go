package authenticity

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestHMMKmerModel_GetModelStats tests model statistics both before and after training.
func TestHMMKmerModel_GetModelStats(t *testing.T) {
	t.Run("untrained", func(t *testing.T) {
		model := NewHMMKmerModel(4)
		require.NotNil(t, model)

		stats := model.GetModelStats()
		require.NotNil(t, stats)
		assert.Equal(t, 4, stats["k"])
		assert.Equal(t, false, stats["trained"])
		assert.Equal(t, uint64(0), stats["total_kmers"])
	})

	t.Run("trained", func(t *testing.T) {
		refGenome := createSyntheticReferenceGenome(nil)
		model := NewHMMKmerModel(4)
		require.NotNil(t, model)

		require.NoError(t, model.TrainFromReferenceGenome(refGenome))

		stats := model.GetModelStats()
		assert.Equal(t, true, stats["trained"])
		assert.Greater(t, stats["total_kmers"].(uint64), uint64(0))
		assert.Contains(t, stats, "transition_entropy")
		assert.Contains(t, stats, "most_probable_start_state")
		assert.Contains(t, stats, "max_start_probability")
	})
}

// TestHMMKmerModel_K tests the K() accessor.
func TestHMMKmerModel_K(t *testing.T) {
	for _, k := range []int{2, 3, 4, 5, 6, 7, 8} {
		model := NewHMMKmerModel(k)
		require.NotNil(t, model)
		assert.Equal(t, k, model.K())
	}
}

// TestHMMKmerModel_IsTrained tests the IsTrained() accessor.
func TestHMMKmerModel_IsTrained(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)
	assert.False(t, model.IsTrained())

	refGenome := createSyntheticReferenceGenome(nil)
	require.NoError(t, model.TrainFromReferenceGenome(refGenome))
	assert.True(t, model.IsTrained())
}

// TestHMMKmerModel_TotalKmers tests the TotalKmers() accessor.
func TestHMMKmerModel_TotalKmers(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)
	assert.Equal(t, uint64(0), model.TotalKmers())

	refGenome := createSyntheticReferenceGenome(nil)
	require.NoError(t, model.TrainFromReferenceGenome(refGenome))
	assert.Greater(t, model.TotalKmers(), uint64(0))
}

// TestHMMKmerModel_StateSpace tests the StateSpace() accessor.
func TestHMMKmerModel_StateSpace(t *testing.T) {
	model := NewHMMKmerModel(3)
	require.NotNil(t, model)

	ss := model.StateSpace()
	assert.Len(t, ss, 16) // 4^(3-1) = 4^2 = 16
}

// TestHMMKmerModel_SetSeed_BeforeTraining tests SetSeed on an untrained model.
func TestHMMKmerModel_SetSeed_BeforeTraining(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)
	// Should not panic.
	model.SetSeed(12345)
}

// TestHMMKmerModel_TrainFromReferenceGenome_SkipsShortChromosome tests short chromosome skip.
func TestHMMKmerModel_TrainFromReferenceGenome_SkipsShortChromosome(t *testing.T) {
	// A mock genome with one chromosome shorter than k.
	genome := &MockReferenceGenome{
		chromosomes: map[string]*MockChromosome{
			"chrShort": {Name: "chrShort", Sequence: []byte("AT"), Length: 2},
			"chrLong":  {Name: "chrLong", Sequence: makeAdaptiveTestSequence(500), Length: 500},
		},
	}

	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	err := model.TrainFromReferenceGenome(genome)
	require.NoError(t, err)
	assert.True(t, model.IsTrained())
	// Only chrLong contributes k-mers.
	assert.Greater(t, model.TotalKmers(), uint64(0))
}

// TestHMMKmerModel_TrainFromReferenceGenome_WithNBases tests chromosome with N bases.
func TestHMMKmerModel_TrainFromReferenceGenome_WithNBases(t *testing.T) {
	seq := append(makeAdaptiveTestSequence(200), []byte("NNNNNNN")...)
	seq = append(seq, makeAdaptiveTestSequence(200)...)

	genome := &MockReferenceGenome{
		chromosomes: map[string]*MockChromosome{
			"chr1": {Name: "chr1", Sequence: seq, Length: len(seq)},
		},
	}

	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	err := model.TrainFromReferenceGenome(genome)
	require.NoError(t, err)
	assert.True(t, model.IsTrained())
}

// TestHMMKmerModel_ValidateModel_MatrixNormalization tests that probabilities sum to 1.
func TestHMMKmerModel_ValidateModel_MatrixNormalization(t *testing.T) {
	refGenome := createSyntheticReferenceGenome(nil)
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)
	require.NoError(t, model.TrainFromReferenceGenome(refGenome))

	err := model.ValidateModel()
	assert.NoError(t, err, "trained model should pass validation (all probs sum to 1)")
}

// TestHMMKmerModel_GenerateSequence_k2 tests sequence generation with k=2.
func TestHMMKmerModel_GenerateSequence_k2(t *testing.T) {
	refGenome := createSyntheticReferenceGenome(nil)
	model := NewHMMKmerModel(2)
	require.NotNil(t, model)
	require.NoError(t, model.TrainFromReferenceGenome(refGenome))

	result, err := model.GenerateSequence(50)
	require.NoError(t, err)
	assert.Len(t, result.Sequence, 50)
	for _, base := range result.Sequence {
		assert.Contains(t, "ATGC", string(base))
	}
}

// TestContainsN tests the containsN helper function.
func TestContainsN(t *testing.T) {
	assert.False(t, containsN("ATCG"))
	assert.True(t, containsN("ATNCG"))
	assert.True(t, containsN("NNNN"))
	assert.False(t, containsN(""))
}

// TestHMMKmerModel_NormalDistribution tests the normal distribution calculation.
func TestHMMKmerModel_NormalDistribution(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	// Normal dist at mean should be highest.
	atMean := model.normalDistribution(20.0, 20.0, 5.0)
	above := model.normalDistribution(25.0, 20.0, 5.0)
	below := model.normalDistribution(15.0, 20.0, 5.0)

	assert.Greater(t, atMean, above)
	assert.Greater(t, atMean, below)
	assert.Greater(t, atMean, 0.0)
}

// TestHMMKmerModel_SimulateQuality tests quality simulation for various k-mers.
func TestHMMKmerModel_SimulateQuality(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	tests := []struct {
		kmer     string
		position int
	}{
		{"ATCG", 0},
		{"GGGG", 10},  // homopolymer
		{"AAAA", 100}, // high position, homopolymer
		{"GCGC", 50},  // high GC
	}

	for _, tt := range tests {
		quality := model.simulateQuality(tt.position, tt.kmer)
		assert.GreaterOrEqual(t, quality, byte(0))
		assert.LessOrEqual(t, quality, byte(41))
	}
}

// TestHMMKmerModel_CalculateGCContent tests GC content calculation.
func TestHMMKmerModel_CalculateGCContent(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	assert.InDelta(t, 0.0, model.calculateGCContent("AAAA"), 1e-9)
	assert.InDelta(t, 1.0, model.calculateGCContent("GCGC"), 1e-9)
	assert.InDelta(t, 0.5, model.calculateGCContent("ATGC"), 1e-9)
}

// TestHMMKmerModel_CalculateHomopolymerPenalty tests penalty computation.
func TestHMMKmerModel_CalculateHomopolymerPenalty(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	// No homopolymer should give 0 penalty.
	assert.InDelta(t, 0.0, model.calculateHomopolymerPenalty("ATCG"), 1e-9)
	assert.InDelta(t, 0.0, model.calculateHomopolymerPenalty("AATG"), 1e-9)

	// Three-in-a-row starts penalizing: penalty = (maxRun-2)*2.0
	assert.InDelta(t, 2.0, model.calculateHomopolymerPenalty("AAA"), 1e-9) // maxRun=3, penalty=(3-2)*2=2
	assert.Greater(t, model.calculateHomopolymerPenalty("AAAA"), 0.0)
	assert.Greater(t, model.calculateHomopolymerPenalty("AAAAA"), model.calculateHomopolymerPenalty("AAAA"))
}

// TestHMMKmerModel_GetNextState_k2 tests the k=2 special case in getNextState.
func TestHMMKmerModel_GetNextState_k2(t *testing.T) {
	model := NewHMMKmerModel(2)
	require.NotNil(t, model)

	// For k=2, state is just the next base index.
	assert.Equal(t, 0, model.getNextState(0, 'A'))
	assert.Equal(t, 1, model.getNextState(0, 'T'))
	assert.Equal(t, 2, model.getNextState(0, 'G'))
	assert.Equal(t, 3, model.getNextState(0, 'C'))
}

// TestHMMKmerModel_GetNextState_kGt2 tests the k>2 sliding window in getNextState.
func TestHMMKmerModel_GetNextState_kGt2(t *testing.T) {
	model := NewHMMKmerModel(3)
	require.NotNil(t, model)

	// State "AT" + base 'G' = new state "TG"
	stateIdx := model.getStateIndex("AT")
	require.NotEqual(t, -1, stateIdx)

	nextIdx := model.getNextState(stateIdx, 'G')
	assert.NotEqual(t, -1, nextIdx)
	assert.Equal(t, "TG", model.stateSpace[nextIdx])
}

// TestHMMKmerModel_GetBaseIndex tests all base cases.
func TestHMMKmerModel_GetBaseIndex(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	assert.Equal(t, 0, model.getBaseIndex('A'))
	assert.Equal(t, 1, model.getBaseIndex('T'))
	assert.Equal(t, 2, model.getBaseIndex('G'))
	assert.Equal(t, 3, model.getBaseIndex('C'))
	assert.Equal(t, -1, model.getBaseIndex('N'))
	assert.Equal(t, -1, model.getBaseIndex('X'))
}

// TestHMMKmerModel_GetStateIndex tests state lookup.
func TestHMMKmerModel_GetStateIndex(t *testing.T) {
	model := NewHMMKmerModel(3)
	require.NotNil(t, model)

	// Known states for k=3 are all (k-1)=2-mer combos.
	assert.NotEqual(t, -1, model.getStateIndex("AT"))
	assert.NotEqual(t, -1, model.getStateIndex("GC"))
	assert.Equal(t, -1, model.getStateIndex("XY")) // Invalid state
}

// TestHMMKmerModel_TrainOnChromosome tests the internal chromosome trainer.
func TestHMMKmerModel_TrainOnChromosome(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	seq := makeAdaptiveTestSequence(500)
	chr := &MockChromosome{Name: "chr1", Sequence: seq, Length: len(seq)}

	transitionCounts := make([][]int, len(model.stateSpace))
	emissionCounts := make([][]int, len(model.stateSpace))
	startCounts := make([]int, len(model.stateSpace))

	for i := range transitionCounts {
		transitionCounts[i] = make([]int, 4)
		emissionCounts[i] = make([]int, 42)
	}

	err := model.trainOnChromosome(chr, transitionCounts, emissionCounts, startCounts)
	require.NoError(t, err)

	// Some transitions should have been counted.
	total := 0
	for _, row := range transitionCounts {
		for _, count := range row {
			total += count
		}
	}
	assert.Greater(t, total, 0)
}

// TestHMMKmerModel_TrainOnChromosome_TooShort tests short chromosome skip in trainOnChromosome.
func TestHMMKmerModel_TrainOnChromosome_TooShort(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	chr := &MockChromosome{Name: "chrX", Sequence: []byte("AT"), Length: 2}
	transitionCounts := make([][]int, len(model.stateSpace))
	emissionCounts := make([][]int, len(model.stateSpace))
	startCounts := make([]int, len(model.stateSpace))
	for i := range transitionCounts {
		transitionCounts[i] = make([]int, 4)
		emissionCounts[i] = make([]int, 42)
	}

	err := model.trainOnChromosome(chr, transitionCounts, emissionCounts, startCounts)
	require.NoError(t, err)

	// No counts should have been added.
	for _, row := range transitionCounts {
		for _, count := range row {
			assert.Equal(t, 0, count)
		}
	}
}
