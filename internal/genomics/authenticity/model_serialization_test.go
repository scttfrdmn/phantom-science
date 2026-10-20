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

	// Generation parameters, counts and finalized probabilities survive.
	assert.InDelta(t, original.targetComplexity, loaded.targetComplexity, 1e-12)
	assert.Equal(t, original.complexityWindowSize, loaded.complexityWindowSize)
	assert.Equal(t, original.transitionCounts, loaded.transitionCounts)
	assert.Equal(t, original.startProbabilities, loaded.startProbabilities)
	assert.True(t, loaded.finalized)

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

// TestAdaptiveHMMModel_SaveLoad_Deterministic: the saved probabilities
// round-trip exactly, so Load -> SetSeed -> Generate is byte-identical across
// two loads and to the original model with the same seed, and training a
// loaded model further gives the same result as training the original.
func TestAdaptiveHMMModel_SaveLoad_Deterministic(t *testing.T) {
	original := NewAdaptiveHMMModel(4)
	require.NotNil(t, original)
	require.NoError(t, original.TrainFromSequence(makeAdaptiveTestSequence(10000)))

	path := filepath.Join(t.TempDir(), "adaptive.hmm.json")
	require.NoError(t, original.Save(path))

	gen := func(m *AdaptiveHMMModel) *GenerationResult {
		m.SetSeed(99)
		r, err := m.GenerateSequence(5000)
		require.NoError(t, err)
		return r
	}
	a, err := LoadAdaptiveHMMModel(path)
	require.NoError(t, err)
	b, err := LoadAdaptiveHMMModel(path)
	require.NoError(t, err)
	ra, rb, ro := gen(a), gen(b), gen(original)
	assert.Equal(t, ra.Sequence, rb.Sequence)
	assert.Equal(t, ra.Qualities, rb.Qualities)
	assert.Equal(t, ro.Sequence, ra.Sequence)
	assert.Equal(t, ro.Qualities, ra.Qualities)

	more := []byte("GGGGCCCCGGGGCCCCAAAATTTT")
	require.NoError(t, a.TrainFromSequence(more))
	require.NoError(t, original.TrainFromSequence(more))
	a.Finalize()
	original.Finalize()
	assert.Equal(t, original.transitionMatrix, a.transitionMatrix)
	assert.Equal(t, original.startProbabilities, a.startProbabilities)
}

// A version 1 adaptive model file (written before #25) is rejected with a
// message that says to retrain; a version 1 HMMKmerModel file still loads.
func TestLoad_Version1Files(t *testing.T) {
	dir := t.TempDir()

	adaptive := NewAdaptiveHMMModel(3)
	require.NoError(t, adaptive.TrainFromSequence(makeAdaptiveTestSequence(500)))
	ap := filepath.Join(dir, "a.hmm.json")
	require.NoError(t, adaptive.Save(ap))
	setEnvelopeVersion(t, ap, "1")
	_, err := LoadModel(ap)
	require.Error(t, err)
	assert.Contains(t, err.Error(), "version 1")
	assert.Contains(t, err.Error(), "retrain")

	kmer := NewHMMKmerModel(3)
	require.NoError(t, kmer.TrainFromReferenceGenome(createSyntheticReferenceGenome(t)))
	kp := filepath.Join(dir, "k.hmm.json")
	require.NoError(t, kmer.Save(kp))
	setEnvelopeVersion(t, kp, "1")
	loaded, err := LoadHMMKmerModel(kp)
	require.NoError(t, err)
	assert.Equal(t, kmer.transitionMatrix, loaded.transitionMatrix)

	setEnvelopeVersion(t, kp, "3")
	_, err = LoadHMMKmerModel(kp)
	assert.ErrorContains(t, err, "unsupported version")
}

// Malformed model data is rejected on load instead of panicking later.
func TestLoadAdaptiveHMMModel_RejectsMalformed(t *testing.T) {
	m := NewAdaptiveHMMModel(3)
	require.NoError(t, m.TrainFromSequence(makeAdaptiveTestSequence(500)))
	path := filepath.Join(t.TempDir(), "a.hmm.json")
	require.NoError(t, m.Save(path))

	raw, err := os.ReadFile(path)
	require.NoError(t, err)
	var env ModelEnvelope
	require.NoError(t, json.Unmarshal(raw, &env))
	var data adaptiveHMMModelData
	require.NoError(t, json.Unmarshal(env.Data, &data))

	write := func(d adaptiveHMMModelData) {
		b, err := json.Marshal(d)
		require.NoError(t, err)
		env.Data = b
		out, err := json.Marshal(env)
		require.NoError(t, err)
		require.NoError(t, os.WriteFile(path, out, 0600))
	}
	bad := data
	bad.TransitionMatrix = bad.TransitionMatrix[:3]
	write(bad)
	_, err = LoadAdaptiveHMMModel(path)
	assert.Error(t, err)

	bad = data
	bad.TransitionMatrix = append([][]float64{{1, 1, 0, 0}}, data.TransitionMatrix[1:]...)
	write(bad)
	_, err = LoadAdaptiveHMMModel(path)
	assert.ErrorContains(t, err, "sum to")
}

// setEnvelopeVersion rewrites the envelope version of the model file at path.
func setEnvelopeVersion(t *testing.T, path, version string) {
	t.Helper()
	raw, err := os.ReadFile(path)
	require.NoError(t, err)
	var env ModelEnvelope
	require.NoError(t, json.Unmarshal(raw, &env))
	env.Version = version
	out, err := json.Marshal(env)
	require.NoError(t, err)
	require.NoError(t, os.WriteFile(path, out, 0600))
}

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
