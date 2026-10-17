//go:build !release

package authenticity

import (
	"encoding/json"
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestHMMKmerModel_SaveLoad_RoundTrip trains a k=4 model, saves, loads, and verifies
// all matrices are bit-for-bit identical and that generation still works.
func TestHMMKmerModel_SaveLoad_RoundTrip(t *testing.T) {
	refGenome := createSyntheticReferenceGenome(t)

	original := NewHMMKmerModel(4)
	require.NotNil(t, original)
	require.NoError(t, original.TrainFromReferenceGenome(refGenome))

	tmp := t.TempDir()
	path := filepath.Join(tmp, "model.hmm.json")

	require.NoError(t, original.Save(path))

	// Confirm the file exists and is non-empty.
	info, err := os.Stat(path)
	require.NoError(t, err)
	assert.Greater(t, info.Size(), int64(0))

	loaded, err := LoadHMMKmerModel(path)
	require.NoError(t, err)
	require.NotNil(t, loaded)

	// Verify structural fields.
	assert.Equal(t, original.k, loaded.k)
	assert.Equal(t, original.trained, loaded.trained)
	assert.Equal(t, original.totalKmers, loaded.totalKmers)
	assert.Equal(t, original.stateSpace, loaded.stateSpace)

	// Verify all transition probabilities match within float64 tolerance.
	require.Equal(t, len(original.transitionMatrix), len(loaded.transitionMatrix))
	for i := range original.transitionMatrix {
		require.Equal(t, len(original.transitionMatrix[i]), len(loaded.transitionMatrix[i]))
		for j := range original.transitionMatrix[i] {
			assert.InDelta(t, original.transitionMatrix[i][j], loaded.transitionMatrix[i][j], 1e-12,
				"transitionMatrix[%d][%d] mismatch", i, j)
		}
	}

	// Verify all emission probabilities match.
	require.Equal(t, len(original.emissionMatrix), len(loaded.emissionMatrix))
	for i := range original.emissionMatrix {
		require.Equal(t, len(original.emissionMatrix[i]), len(loaded.emissionMatrix[i]))
		for j := range original.emissionMatrix[i] {
			assert.InDelta(t, original.emissionMatrix[i][j], loaded.emissionMatrix[i][j], 1e-12,
				"emissionMatrix[%d][%d] mismatch", i, j)
		}
	}

	// Verify start probabilities match.
	require.Equal(t, len(original.startProbabilities), len(loaded.startProbabilities))
	for i := range original.startProbabilities {
		assert.InDelta(t, original.startProbabilities[i], loaded.startProbabilities[i], 1e-12,
			"startProbabilities[%d] mismatch", i)
	}

	// Verify the loaded model can generate sequences.
	result, err := loaded.GenerateSequence(500)
	require.NoError(t, err)
	assert.Len(t, result.Sequence, 500)
	assert.Len(t, result.Qualities, 500)
	assert.Greater(t, len(result.KmerFrequencies), 0)
}

// TestAdaptiveHMMModel_SaveLoad_RoundTrip trains an AdaptiveHMMModel, saves, loads, and
// verifies biological constraint fields and matrices survive the round-trip.
func TestAdaptiveHMMModel_SaveLoad_RoundTrip(t *testing.T) {
	original := NewAdaptiveHMMModel(4)
	require.NotNil(t, original)
	require.NoError(t, original.TrainFromSequence(makeAdaptiveTestSequence(10000)))

	tmp := t.TempDir()
	path := filepath.Join(tmp, "adaptive.hmm.json")

	require.NoError(t, original.Save(path))

	loaded, err := LoadAdaptiveHMMModel(path)
	require.NoError(t, err)
	require.NotNil(t, loaded)

	// Verify structural fields.
	assert.Equal(t, original.k, loaded.k)
	assert.Equal(t, original.trained, loaded.trained)
	assert.Equal(t, original.totalKmers, loaded.totalKmers)
	assert.Equal(t, original.stateSpace, loaded.stateSpace)

	// Verify biological constraint fields survive serialization.
	assert.InDelta(t, original.gcContentTarget, loaded.gcContentTarget, 1e-12)
	assert.InDelta(t, original.gcContentTolerance, loaded.gcContentTolerance, 1e-12)
	assert.Equal(t, original.homopolymerMaxLength, loaded.homopolymerMaxLength)
	assert.InDelta(t, original.cpgDepletionRatio, loaded.cpgDepletionRatio, 1e-12)
	assert.InDelta(t, original.targetComplexity, loaded.targetComplexity, 1e-12)
	assert.Equal(t, original.complexityWindowSize, loaded.complexityWindowSize)

	// Verify adaptive weights map is preserved.
	assert.Equal(t, len(original.adaptiveWeights), len(loaded.adaptiveWeights))
	for state, w := range original.adaptiveWeights {
		assert.InDelta(t, w, loaded.adaptiveWeights[state], 1e-12, "adaptiveWeights[%q] mismatch", state)
	}

	// Verify complexity memory was reset to empty slice (not serialized).
	assert.Empty(t, loaded.complexityMemory)

	// Verify transition and emission matrices match.
	require.Equal(t, len(original.transitionMatrix), len(loaded.transitionMatrix))
	for i := range original.transitionMatrix {
		for j := range original.transitionMatrix[i] {
			assert.InDelta(t, original.transitionMatrix[i][j], loaded.transitionMatrix[i][j], 1e-12,
				"transitionMatrix[%d][%d] mismatch", i, j)
		}
	}

	// Verify loaded model can generate sequences.
	result, err := loaded.GenerateSequence(500)
	require.NoError(t, err)
	assert.Len(t, result.Sequence, 500)
}

// TestHMMKmerModel_Save_UntrainedReturnsError verifies that saving an untrained model
// returns an error rather than writing a corrupt file.
func TestHMMKmerModel_Save_UntrainedReturnsError(t *testing.T) {
	model := NewHMMKmerModel(4)
	require.NotNil(t, model)

	tmp := t.TempDir()
	err := model.Save(filepath.Join(tmp, "untrained.hmm.json"))
	assert.Error(t, err)

	// Confirm no file was written.
	_, statErr := os.Stat(filepath.Join(tmp, "untrained.hmm.json"))
	assert.True(t, os.IsNotExist(statErr), "file should not exist after failed save")
}

// TestLoadHMMKmerModel_InvalidPath verifies that loading from a non-existent path
// returns an error.
func TestLoadHMMKmerModel_InvalidPath(t *testing.T) {
	_, err := LoadHMMKmerModel("/nonexistent/path/model.hmm.json")
	assert.Error(t, err)
}

// TestLoadHMMKmerModel_CorruptData verifies that loading a file containing garbage JSON
// returns an error.
func TestLoadHMMKmerModel_CorruptData(t *testing.T) {
	tmp := t.TempDir()
	path := filepath.Join(tmp, "corrupt.hmm.json")
	require.NoError(t, os.WriteFile(path, []byte("this is not valid json {{{"), 0600))

	_, err := LoadHMMKmerModel(path)
	assert.Error(t, err)
}

// TestAdaptiveHMMModel_SaveLoad_GeneratesConsistentSequences verifies that sequences
// generated by the original and the loaded model have statistically similar GC content
// and complexity.  The rng is reseeded on load so sequences are not bit-for-bit equal,
// but biological properties should be equivalent.
func TestAdaptiveHMMModel_SaveLoad_GeneratesConsistentSequences(t *testing.T) {
	original := NewAdaptiveHMMModel(4)
	require.NotNil(t, original)
	require.NoError(t, original.TrainFromSequence(makeAdaptiveTestSequence(10000)))

	tmp := t.TempDir()
	path := filepath.Join(tmp, "adaptive.hmm.json")
	require.NoError(t, original.Save(path))

	loaded, err := LoadAdaptiveHMMModel(path)
	require.NoError(t, err)

	const seqLen = 5000
	origResult, err := original.GenerateSequence(seqLen)
	require.NoError(t, err)
	loadResult, err := loaded.GenerateSequence(seqLen)
	require.NoError(t, err)

	origGC := computeGCContent(origResult.Sequence)
	loadGC := computeGCContent(loadResult.Sequence)

	origComplexity := computeKmerComplexity(origResult.Sequence, original.k)
	loadComplexity := computeKmerComplexity(loadResult.Sequence, loaded.k)

	// GC content should be in the same ballpark (within 0.10 of each other).
	assert.InDelta(t, origGC, loadGC, 0.10,
		"GC content of original (%.3f) and loaded (%.3f) models should be similar", origGC, loadGC)

	// Complexity should be similar (within 0.15).
	assert.InDelta(t, origComplexity, loadComplexity, 0.15,
		"complexity of original (%.3f) and loaded (%.3f) models should be similar", origComplexity, loadComplexity)
}

// computeGCContent calculates the fraction of G/C bases in a sequence.
func computeGCContent(seq []byte) float64 {
	if len(seq) == 0 {
		return 0
	}
	gc := 0
	for _, b := range seq {
		if b == 'G' || b == 'C' {
			gc++
		}
	}
	return float64(gc) / float64(len(seq))
}

// computeKmerComplexity calculates the ratio of unique k-mers to total k-mers.
func computeKmerComplexity(seq []byte, k int) float64 {
	if len(seq) < k {
		return 0
	}
	unique := make(map[string]struct{})
	total := 0
	for i := 0; i <= len(seq)-k; i++ {
		unique[string(seq[i:i+k])] = struct{}{}
		total++
	}
	return float64(len(unique)) / float64(total)
}

// TestLoadModel_HMMKmerModel verifies that LoadModel dispatches to LoadHMMKmerModel.
func TestLoadModel_HMMKmerModel(t *testing.T) {
	refGenome := createSyntheticReferenceGenome(t)
	m := NewHMMKmerModel(4)
	require.NotNil(t, m)
	require.NoError(t, m.TrainFromReferenceGenome(refGenome))

	path := filepath.Join(t.TempDir(), "kmer.hmm.json")
	require.NoError(t, m.Save(path))

	gen, err := LoadModel(path)
	require.NoError(t, err)
	require.NotNil(t, gen)
	assert.True(t, gen.IsTrained())

	result, err := gen.GenerateSequence(200)
	require.NoError(t, err)
	assert.Len(t, result.Sequence, 200)
}

// TestLoadModel_AdaptiveHMMModel verifies that LoadModel dispatches to LoadAdaptiveHMMModel.
func TestLoadModel_AdaptiveHMMModel(t *testing.T) {
	m := NewAdaptiveHMMModel(4)
	require.NotNil(t, m)
	require.NoError(t, m.TrainFromSequence(makeAdaptiveTestSequence(5000)))

	path := filepath.Join(t.TempDir(), "adaptive.hmm.json")
	require.NoError(t, m.Save(path))

	gen, err := LoadModel(path)
	require.NoError(t, err)
	require.NotNil(t, gen)
	assert.True(t, gen.IsTrained())

	result, err := gen.GenerateSequence(200)
	require.NoError(t, err)
	assert.Len(t, result.Sequence, 200)
}

// TestLoadModel_MissingPath verifies that LoadModel returns an error for a missing file.
func TestLoadModel_MissingPath(t *testing.T) {
	_, err := LoadModel("/nonexistent/path/model.hmm.json")
	assert.Error(t, err)
}

// TestLoadModel_CorruptJSON verifies that LoadModel returns an error for corrupt JSON.
func TestLoadModel_CorruptJSON(t *testing.T) {
	path := filepath.Join(t.TempDir(), "corrupt.hmm.json")
	require.NoError(t, os.WriteFile(path, []byte("not valid json {{{"), 0600))
	_, err := LoadModel(path)
	assert.Error(t, err)
}

// TestLoadModel_UnknownModelType verifies that LoadModel returns an error for an
// unrecognised ModelType in the envelope.
func TestLoadModel_UnknownModelType(t *testing.T) {
	env := ModelEnvelope{
		Format:    "phantom-science-hmm",
		Version:   "1",
		ModelType: "UnknownModel",
	}
	raw, err := json.Marshal(env)
	require.NoError(t, err)

	path := filepath.Join(t.TempDir(), "unknown.hmm.json")
	require.NoError(t, os.WriteFile(path, raw, 0600))

	_, err = LoadModel(path)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "unknown model type")
}
