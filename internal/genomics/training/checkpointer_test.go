//go:build !release

package training

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/quality"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// makeTrainedKmerModel creates a trained k=4 HMMKmerModel using a synthetic FASTA genome.
func makeTrainedKmerModel(t *testing.T) *authenticity.HMMKmerModel {
	t.Helper()
	dir := makeSyntheticFASTADir(t)
	genome := reference.NewReferenceGenome(dir)
	model := authenticity.NewHMMKmerModel(4)
	require.NotNil(t, model)
	err := model.TrainFromReferenceGenome(genome)
	require.NoError(t, err)
	return model
}

// makeTrainedAdaptiveModel creates a trained AdaptiveHMMModel using a short synthetic sequence.
func makeTrainedAdaptiveModel(t *testing.T) *authenticity.AdaptiveHMMModel {
	t.Helper()
	model := authenticity.NewAdaptiveHMMModel(4)
	require.NotNil(t, model)
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	err := model.TrainFromSequence(seq)
	require.NoError(t, err)
	return model
}

// --- ModelCheckpointer.saveModel tests ---

func TestModelCheckpointer_SaveModel_EmptyPath(t *testing.T) {
	c := NewModelCheckpointer("")
	model := authenticity.NewHMMKmerModel(4)

	err := c.saveModel(model, "test_key")
	assert.NoError(t, err, "empty output path should be a no-op")
}

func TestModelCheckpointer_SaveModel_HMMKmerModel(t *testing.T) {
	outDir := t.TempDir()
	c := NewModelCheckpointer(outDir)
	model := makeTrainedKmerModel(t)

	err := c.saveModel(model, "4mer_YRI")
	require.NoError(t, err)

	_, statErr := os.Stat(filepath.Join(outDir, "4mer_YRI.hmm.json"))
	assert.NoError(t, statErr, "saved file should exist")
}

func TestModelCheckpointer_SaveModel_AdaptiveHMMModel(t *testing.T) {
	outDir := t.TempDir()
	c := NewModelCheckpointer(outDir)
	model := makeTrainedAdaptiveModel(t)

	err := c.saveModel(model, "adaptive_YRI")
	require.NoError(t, err)

	_, statErr := os.Stat(filepath.Join(outDir, "adaptive_YRI.hmm.json"))
	assert.NoError(t, statErr, "saved file should exist")
}

func TestModelCheckpointer_SaveModel_UnsupportedType(t *testing.T) {
	outDir := t.TempDir()
	c := NewModelCheckpointer(outDir)

	err := c.saveModel("not_a_model", "key")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "unsupported model type")
}

func TestModelCheckpointer_SaveModel_NestedKey(t *testing.T) {
	outDir := t.TempDir()
	c := NewModelCheckpointer(outDir)
	model := makeTrainedKmerModel(t)

	// Key with a subdirectory component
	err := c.saveModel(model, "subdir/4mer_CEU")
	require.NoError(t, err)

	_, statErr := os.Stat(filepath.Join(outDir, "subdir", "4mer_CEU.hmm.json"))
	assert.NoError(t, statErr, "saved file in subdirectory should exist")
}

// --- ModelCheckpointer.loadModel tests ---

func TestModelCheckpointer_LoadModel_MissingFile(t *testing.T) {
	outDir := t.TempDir()
	c := NewModelCheckpointer(outDir)

	_, err := c.loadModel("nonexistent_key")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "checkpointer: read")
}

func TestModelCheckpointer_LoadModel_CorruptJSON(t *testing.T) {
	outDir := t.TempDir()
	c := NewModelCheckpointer(outDir)

	err := os.WriteFile(filepath.Join(outDir, "bad_model.hmm.json"), []byte("{not valid json}"), 0600)
	require.NoError(t, err)

	_, loadErr := c.loadModel("bad_model")
	require.Error(t, loadErr)
	assert.Contains(t, loadErr.Error(), "checkpointer: parse envelope")
}

func TestModelCheckpointer_LoadModel_UnknownModelType(t *testing.T) {
	outDir := t.TempDir()
	c := NewModelCheckpointer(outDir)

	env := authenticity.ModelEnvelope{
		Format:    "phantom-hmm",
		Version:   "1.0",
		ModelType: "UnknownModelType",
	}
	raw, err := json.Marshal(env)
	require.NoError(t, err)
	err = os.WriteFile(filepath.Join(outDir, "unknown.hmm.json"), raw, 0600)
	require.NoError(t, err)

	_, loadErr := c.loadModel("unknown")
	require.Error(t, loadErr)
	assert.Contains(t, loadErr.Error(), "unknown model type")
}

func TestModelCheckpointer_LoadModel_HMMKmerModel(t *testing.T) {
	outDir := t.TempDir()
	c := NewModelCheckpointer(outDir)
	original := makeTrainedKmerModel(t)

	err := c.saveModel(original, "kmer_load_test")
	require.NoError(t, err)

	loaded, err := c.loadModel("kmer_load_test")
	require.NoError(t, err)

	km, ok := loaded.(*authenticity.HMMKmerModel)
	require.True(t, ok, "loaded model should be *HMMKmerModel")
	assert.True(t, km.IsTrained())
}

func TestModelCheckpointer_LoadModel_AdaptiveHMMModel(t *testing.T) {
	outDir := t.TempDir()
	c := NewModelCheckpointer(outDir)
	original := makeTrainedAdaptiveModel(t)

	err := c.saveModel(original, "adaptive_load_test")
	require.NoError(t, err)

	loaded, err := c.loadModel("adaptive_load_test")
	require.NoError(t, err)

	am, ok := loaded.(*authenticity.AdaptiveHMMModel)
	require.True(t, ok, "loaded model should be *AdaptiveHMMModel")
	assert.True(t, am.IsTrained())
}

// --- trainModels success path tests ---

// TestTrainingPipeline_TrainModels_Success injects working loaders into the
// pipeline's sub-trainers and calls trainModels directly.
func TestTrainingPipeline_TrainModels_Success(t *testing.T) {
	fastaDir := makeSyntheticFASTADir(t)

	config := &TrainingConfig{
		KmerSizes:   []int{4},
		Populations: []string{"YRI"},
		Platforms:   []string{"illumina"},
	}
	pipeline := NewTrainingPipeline(config)

	pipeline.modelTrainer.hmmTrainer.genomeLoader = func(pop string) (*reference.ReferenceGenome, error) {
		return reference.NewReferenceGenome(fastaDir), nil
	}
	pipeline.modelTrainer.qualityTrainer.dataLoader = func(platform string) ([]quality.GCQualityDataPoint, error) {
		return syntheticQualityData(platform), nil
	}
	pipeline.modelTrainer.popTrainer.paramsLoader = func(pop string) (*PopulationParameters, error) {
		return syntheticPopulationParams(), nil
	}

	err := pipeline.trainModels()
	require.NoError(t, err)
	// 1 HMM (k=4, YRI) + 1 quality (illumina) + 1 population (YRI) = 3
	assert.Equal(t, 3, pipeline.metrics.ModelsTrained)
}

// TestTrainingPipeline_TrainModels_HMMError verifies trainModels returns an
// error when HMM training fails.
func TestTrainingPipeline_TrainModels_HMMError(t *testing.T) {
	config := &TrainingConfig{
		KmerSizes:   []int{4},
		Populations: []string{"YRI"},
		Platforms:   []string{},
	}
	pipeline := NewTrainingPipeline(config)
	pipeline.modelTrainer.hmmTrainer.genomeLoader = func(pop string) (*reference.ReferenceGenome, error) {
		return nil, fmt.Errorf("synthetic genome error")
	}

	err := pipeline.trainModels()
	require.Error(t, err)
	assert.Contains(t, err.Error(), "failed to train HMM model")
}

// TestTrainingPipeline_TrainModels_QualityError verifies trainModels returns
// an error when quality training fails.
func TestTrainingPipeline_TrainModels_QualityError(t *testing.T) {
	config := &TrainingConfig{
		KmerSizes:   []int{},
		Populations: []string{},
		Platforms:   []string{"illumina"},
	}
	pipeline := NewTrainingPipeline(config)
	pipeline.modelTrainer.qualityTrainer.dataLoader = func(platform string) ([]quality.GCQualityDataPoint, error) {
		return nil, fmt.Errorf("synthetic quality error")
	}

	err := pipeline.trainModels()
	require.Error(t, err)
	assert.Contains(t, err.Error(), "failed to train quality model")
}

// TestTrainingPipeline_TrainModels_PopulationError verifies trainModels returns
// an error when population training fails.
func TestTrainingPipeline_TrainModels_PopulationError(t *testing.T) {
	config := &TrainingConfig{
		KmerSizes:   []int{},
		Populations: []string{"YRI"},
		Platforms:   []string{},
	}
	pipeline := NewTrainingPipeline(config)
	pipeline.modelTrainer.popTrainer.paramsLoader = func(pop string) (*PopulationParameters, error) {
		return nil, fmt.Errorf("synthetic population error")
	}

	err := pipeline.trainModels()
	require.Error(t, err)
	assert.Contains(t, err.Error(), "failed to train population model")
}

// TestTrainingPipeline_Execute_TrainModelsFails verifies Execute returns
// "model training failed" when trainModels encounters an error.
func TestTrainingPipeline_Execute_TrainModelsFails(t *testing.T) {
	config := &TrainingConfig{
		KmerSizes:   []int{4},
		Populations: []string{"YRI"},
		Platforms:   []string{},
	}
	pipeline := NewTrainingPipeline(config)

	// Pre-populate reference cache so prepareData succeeds
	pipeline.dataManager.referenceCache["GRCh38"] = nil
	pipeline.dataManager.referenceCache["GRCh39"] = nil

	// Inject a failing genome loader so trainModels errors
	pipeline.modelTrainer.hmmTrainer.genomeLoader = func(pop string) (*reference.ReferenceGenome, error) {
		return nil, fmt.Errorf("injected genome failure")
	}

	err := pipeline.Execute()
	require.Error(t, err)
	assert.Contains(t, err.Error(), "model training failed")
}

// TestTrainingPipeline_Execute_PrepareDataFails_BadReference exercises Execute
// when prepareData fails because the GRCh38 reference genome directory is missing.
func TestTrainingPipeline_Execute_PrepareDataFails_BadReference(t *testing.T) {
	config := &TrainingConfig{
		ReferenceGenomePath: "/nonexistent/path/that/cannot/exist_xyz",
		KmerSizes:           []int{},
		Populations:         []string{},
		Platforms:           []string{},
	}
	pipeline := NewTrainingPipeline(config)

	err := pipeline.Execute()
	require.Error(t, err)
	assert.Contains(t, err.Error(), "data preparation failed")
}
