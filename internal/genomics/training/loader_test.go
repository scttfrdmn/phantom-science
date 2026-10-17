//go:build !release

package training

import (
	"encoding/json"
	"os"
	"path/filepath"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/quality"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// --- HMMTrainer.loadReferenceGenome production path ---

func TestHMMTrainer_LoadReferenceGenome_NilHook_EmptyPath(t *testing.T) {
	trainer := NewHMMTrainer(&TrainingConfig{})
	_, err := trainer.loadReferenceGenome("YRI")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "ReferenceGenomePath not configured")
}

func TestHMMTrainer_LoadReferenceGenome_NilHook_BaseDir(t *testing.T) {
	dir := makeSyntheticFASTADir(t)
	trainer := NewHMMTrainer(&TrainingConfig{ReferenceGenomePath: dir})
	genome, err := trainer.loadReferenceGenome("YRI")
	require.NoError(t, err)
	assert.NotNil(t, genome)
}

func TestHMMTrainer_LoadReferenceGenome_NilHook_PopSubdir(t *testing.T) {
	// When a population-specific subdirectory exists, it should be preferred.
	baseDir := t.TempDir()
	popDir := filepath.Join(baseDir, "YRI")
	require.NoError(t, os.MkdirAll(popDir, 0750))
	seq := "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
	require.NoError(t, os.WriteFile(filepath.Join(popDir, "chr1.fasta"),
		[]byte(">chr1\n"+seq+"\n"), 0600))

	trainer := NewHMMTrainer(&TrainingConfig{ReferenceGenomePath: baseDir})
	genome, err := trainer.loadReferenceGenome("YRI")
	require.NoError(t, err)
	assert.NotNil(t, genome)
}

// --- QualityTrainer.loadQualityTrainingData production path ---

func TestQualityTrainer_LoadQualityTrainingData_NilHook_EmptyPath(t *testing.T) {
	trainer := NewQualityTrainer(&TrainingConfig{})
	_, err := trainer.loadQualityTrainingData("illumina")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "QualityDataPath not configured")
}

func TestQualityTrainer_LoadQualityTrainingData_NilHook_MissingFile(t *testing.T) {
	dir := t.TempDir()
	trainer := NewQualityTrainer(&TrainingConfig{QualityDataPath: dir})
	_, err := trainer.loadQualityTrainingData("illumina")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "read")
}

func TestQualityTrainer_LoadQualityTrainingData_NilHook_CorruptJSON(t *testing.T) {
	dir := t.TempDir()
	require.NoError(t, os.WriteFile(filepath.Join(dir, "illumina.json"), []byte("{not json}"), 0600))
	trainer := NewQualityTrainer(&TrainingConfig{QualityDataPath: dir})
	_, err := trainer.loadQualityTrainingData("illumina")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "parse")
}

func TestQualityTrainer_LoadQualityTrainingData_NilHook_EmptyArray(t *testing.T) {
	dir := t.TempDir()
	require.NoError(t, os.WriteFile(filepath.Join(dir, "illumina.json"), []byte("[]"), 0600))
	trainer := NewQualityTrainer(&TrainingConfig{QualityDataPath: dir})
	_, err := trainer.loadQualityTrainingData("illumina")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "no data points")
}

func TestQualityTrainer_LoadQualityTrainingData_NilHook_ValidFile(t *testing.T) {
	dir := t.TempDir()
	points := syntheticQualityData("illumina")
	raw, err := json.Marshal(points)
	require.NoError(t, err)
	require.NoError(t, os.WriteFile(filepath.Join(dir, "illumina.json"), raw, 0600))

	trainer := NewQualityTrainer(&TrainingConfig{QualityDataPath: dir})
	loaded, err := trainer.loadQualityTrainingData("illumina")
	require.NoError(t, err)
	assert.Len(t, loaded, len(points))
	assert.Equal(t, points[0].GCContent, loaded[0].GCContent)
}

func TestQualityTrainer_LoadQualityTrainingData_NilHook_PlatformRouting(t *testing.T) {
	// Different platforms read their own files.
	dir := t.TempDir()
	for _, platform := range []string{"illumina", "pacbio", "nanopore"} {
		raw, err := json.Marshal(syntheticQualityData(platform))
		require.NoError(t, err)
		require.NoError(t, os.WriteFile(filepath.Join(dir, platform+".json"), raw, 0600))
	}

	trainer := NewQualityTrainer(&TrainingConfig{QualityDataPath: dir})
	for _, platform := range []string{"illumina", "pacbio", "nanopore"} {
		pts, err := trainer.loadQualityTrainingData(platform)
		require.NoError(t, err, "platform %s", platform)
		assert.NotEmpty(t, pts, "platform %s", platform)
		assert.Equal(t, platform, pts[0].Platform)
	}
}

// --- PopulationTrainer.getPopulationParameters production path ---

func TestPopulationTrainer_GetPopulationParameters_NilHook_EmptyPath(t *testing.T) {
	trainer := NewPopulationTrainer(&TrainingConfig{})
	_, err := trainer.getPopulationParameters("YRI")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "PopulationDataPath not configured")
}

func TestPopulationTrainer_GetPopulationParameters_NilHook_MissingFile(t *testing.T) {
	dir := t.TempDir()
	trainer := NewPopulationTrainer(&TrainingConfig{PopulationDataPath: dir})
	_, err := trainer.getPopulationParameters("YRI")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "read")
}

func TestPopulationTrainer_GetPopulationParameters_NilHook_CorruptJSON(t *testing.T) {
	dir := t.TempDir()
	require.NoError(t, os.WriteFile(filepath.Join(dir, "YRI.json"), []byte("{bad}"), 0600))
	trainer := NewPopulationTrainer(&TrainingConfig{PopulationDataPath: dir})
	_, err := trainer.getPopulationParameters("YRI")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "parse")
}

func TestPopulationTrainer_GetPopulationParameters_NilHook_ValidFile(t *testing.T) {
	dir := t.TempDir()
	params := syntheticPopulationParams()
	raw, err := json.Marshal(params)
	require.NoError(t, err)
	require.NoError(t, os.WriteFile(filepath.Join(dir, "YRI.json"), raw, 0600))

	trainer := NewPopulationTrainer(&TrainingConfig{PopulationDataPath: dir})
	loaded, err := trainer.getPopulationParameters("YRI")
	require.NoError(t, err)
	assert.Equal(t, params.EffectivePopSize, loaded.EffectivePopSize)
	assert.InDelta(t, params.MutationRate, loaded.MutationRate, 1e-12)
	assert.InDelta(t, params.RecombinationRate, loaded.RecombinationRate, 1e-12)
}

func TestPopulationTrainer_GetPopulationParameters_NilHook_PopulationRouting(t *testing.T) {
	dir := t.TempDir()
	for _, pop := range []string{"YRI", "CEU", "CHB"} {
		p := &PopulationParameters{
			EffectivePopSize:  10_000,
			MutationRate:      1.2e-8,
			RecombinationRate: 1.0e-8,
		}
		raw, _ := json.Marshal(p)
		require.NoError(t, os.WriteFile(filepath.Join(dir, pop+".json"), raw, 0600))
	}

	trainer := NewPopulationTrainer(&TrainingConfig{PopulationDataPath: dir})
	for _, pop := range []string{"YRI", "CEU", "CHB"} {
		loaded, err := trainer.getPopulationParameters(pop)
		require.NoError(t, err, "pop %s", pop)
		assert.Equal(t, 10_000, loaded.EffectivePopSize)
	}
}

// jsonRoundTrip verifies GCQualityDataPoint marshals and unmarshals cleanly.
func TestGCQualityDataPoint_JSONRoundTrip(t *testing.T) {
	original := quality.GCQualityDataPoint{
		GCContent:      0.45,
		QualityScore:   37.5,
		SequenceLength: 150,
		Platform:       "illumina",
		Cycle:          3,
		ReadPosition:   75,
	}
	raw, err := json.Marshal(original)
	require.NoError(t, err)
	var decoded quality.GCQualityDataPoint
	require.NoError(t, json.Unmarshal(raw, &decoded))
	assert.Equal(t, original, decoded)
}
