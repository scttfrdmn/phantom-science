//go:build !release

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

// TestHMMKmerModel_SimulateQuality checks the heuristic quality exactly:
// 30 - 10*|GC - 0.5| - 2*(longest run - 2 if >= 3). It does not depend on the
// offset into the training sequence (an earlier version reached 0 after
// about 3,000 bases of a chromosome).
func TestHMMKmerModel_SimulateQuality(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	for kmer, want := range map[string]byte{"ATCG": 30, "GGGG": 21, "AAAA": 21, "GCGC": 25, "AAAT": 23} {
		assert.Equal(t, want, model.simulateQuality(kmer), kmer)
	}
}

// Training qualities must not collapse on long chromosomes: on 20 kb of
// balanced sequence the mean generated quality stays near the heuristic's
// values (>= 20), not near 0.
func TestHMMKmerModel_QualitiesNotDegenerateOnLongChromosome(t *testing.T) {
	seq := makeAdaptiveTestSequence(20000)
	genome := &MockReferenceGenome{chromosomes: map[string]*MockChromosome{
		"chr1": {Name: "chr1", Sequence: seq, Length: len(seq)},
	}}
	model := NewHMMKmerModel(4)
	require.NoError(t, model.TrainFromReferenceGenome(genome))
	model.SetSeed(1)
	res, err := model.GenerateSequence(5000)
	require.NoError(t, err)
	sum := 0
	for _, q := range res.Qualities {
		sum += int(q)
	}
	mean := float64(sum) / float64(len(res.Qualities))
	assert.Greater(t, mean, 20.0)
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

// TestHMMKmerModel_TrainOnChromosome checks the chromosome trainer counts
// exactly the k-mers of A/C/G/T: "ACGTNACGTA" with k=3 has the valid 3-mers
// ACG, CGT, ACG, CGT, GTA (those spanning N are skipped).
func TestHMMKmerModel_TrainOnChromosome(t *testing.T) {
	model := NewHMMKmerModel(3)
	require.NotNil(t, model)

	chr := &MockChromosome{Name: "chr1", Sequence: []byte("ACGTNACGTA"), Length: 10}
	transitionCounts := newCountMatrix(len(model.stateSpace), 4)
	emissionCounts := make([][]int, len(model.stateSpace))
	for i := range emissionCounts {
		emissionCounts[i] = make([]int, 42)
	}

	n := model.trainOnChromosome(chr, transitionCounts, emissionCounts)
	assert.Equal(t, uint64(5), n)
	ac, cg, gt := model.getStateIndex("AC"), model.getStateIndex("CG"), model.getStateIndex("GT")
	assert.Equal(t, 2.0, transitionCounts[ac][baseIndex('G')])
	assert.Equal(t, 2.0, transitionCounts[cg][baseIndex('T')])
	assert.Equal(t, 1.0, transitionCounts[gt][baseIndex('A')])
	total := 0.0
	for _, row := range transitionCounts {
		for _, c := range row {
			total += c
		}
	}
	assert.Equal(t, 5.0, total)

	short := &MockChromosome{Name: "chrX", Sequence: []byte("AT"), Length: 2}
	assert.Equal(t, uint64(0), model.trainOnChromosome(short, transitionCounts, emissionCounts))
}
