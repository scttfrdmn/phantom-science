//go:build !release

package training

import (
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/population"
	"github.com/scttfrdmn/phantom-science/internal/genomics/quality"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// makeSyntheticFASTADir creates a temp directory with a single chr1.fasta file
// containing a 100-base ATGC sequence — enough to train a k=4 HMM.
func makeSyntheticFASTADir(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	seq := "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
	content := fmt.Sprintf(">chr1\n%s\n", seq)
	err := os.WriteFile(filepath.Join(dir, "chr1.fasta"), []byte(content), 0600)
	require.NoError(t, err)
	return dir
}

// syntheticQualityData returns a small slice of GCQualityDataPoints that span
// multiple GC bins so fitBiasParameters produces a valid (≥ 0) biasSlope.
func syntheticQualityData(platform string) []quality.GCQualityDataPoint {
	return []quality.GCQualityDataPoint{
		{GCContent: 0.30, QualityScore: 35.0, SequenceLength: 100, Platform: platform, Cycle: 1, ReadPosition: 1},
		{GCContent: 0.50, QualityScore: 38.0, SequenceLength: 100, Platform: platform, Cycle: 1, ReadPosition: 1},
		{GCContent: 0.70, QualityScore: 32.0, SequenceLength: 100, Platform: platform, Cycle: 1, ReadPosition: 1},
		{GCContent: 0.30, QualityScore: 34.0, SequenceLength: 150, Platform: platform, Cycle: 2, ReadPosition: 50},
		{GCContent: 0.50, QualityScore: 37.0, SequenceLength: 150, Platform: platform, Cycle: 2, ReadPosition: 50},
	}
}

// syntheticPopulationParams returns a *PopulationParameters that passes ValidateParameters.
func syntheticPopulationParams() *PopulationParameters {
	return &PopulationParameters{
		EffectivePopSize:  10_000,
		MutationRate:      1.2e-8,
		RecombinationRate: 1.0e-8,
		DemographicEvents: nil,
	}
}

// --- HMMTrainer success paths ---

// TestHMMTrainer_TrainHMMModel_Success exercises the full success path of
// TrainHMMModel by injecting a genomeLoader that returns a temp-dir genome.
func TestHMMTrainer_TrainHMMModel_Success(t *testing.T) {
	dir := makeSyntheticFASTADir(t)

	config := &TrainingConfig{}
	trainer := NewHMMTrainer(config)
	trainer.genomeLoader = func(population string) (*reference.ReferenceGenome, error) {
		return reference.NewReferenceGenome(dir), nil
	}

	err := trainer.TrainHMMModel(4, "YRI")
	require.NoError(t, err)

	// Model must be stored and marked trained
	key := "4mer_YRI"
	model, ok := trainer.models[key]
	require.True(t, ok, "model should be stored under key %q", key)
	assert.True(t, model.IsTrained(), "HMMKmerModel should report IsTrained() == true")
	assert.True(t, trainer.trainedModels[key])
}

// TestHMMTrainer_TrainHMMModel_Success_WithSave exercises saveModel by providing
// an OutputPath, verifying the .hmm.json file is written to disk.
func TestHMMTrainer_TrainHMMModel_Success_WithSave(t *testing.T) {
	fastaDir := makeSyntheticFASTADir(t)
	outDir := t.TempDir()

	config := &TrainingConfig{OutputPath: outDir}
	trainer := NewHMMTrainer(config)
	trainer.genomeLoader = func(pop string) (*reference.ReferenceGenome, error) {
		return reference.NewReferenceGenome(fastaDir), nil
	}

	err := trainer.TrainHMMModel(4, "CEU")
	require.NoError(t, err)

	savedPath := filepath.Join(outDir, "4mer_CEU.hmm.json")
	_, statErr := os.Stat(savedPath)
	assert.NoError(t, statErr, "saved model file should exist at %s", savedPath)
}

// TestHMMTrainer_TrainHMMModel_GenomeLoaderError verifies that a genomeLoader
// returning an error causes TrainHMMModel to return that error.
func TestHMMTrainer_TrainHMMModel_GenomeLoaderError(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewHMMTrainer(config)
	trainer.genomeLoader = func(pop string) (*reference.ReferenceGenome, error) {
		return nil, fmt.Errorf("injected genome load failure")
	}

	err := trainer.TrainHMMModel(4, "YRI")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "injected genome load failure")
}

// TestHMMTrainer_TrainHMMModel_MultipleK exercises training with different k values.
func TestHMMTrainer_TrainHMMModel_MultipleK(t *testing.T) {
	dir := makeSyntheticFASTADir(t)

	config := &TrainingConfig{}
	trainer := NewHMMTrainer(config)
	trainer.genomeLoader = func(pop string) (*reference.ReferenceGenome, error) {
		return reference.NewReferenceGenome(dir), nil
	}

	for _, k := range []int{2, 3, 4} {
		t.Run(fmt.Sprintf("k=%d", k), func(t *testing.T) {
			err := trainer.TrainHMMModel(k, "EUR")
			require.NoError(t, err)
			key := fmt.Sprintf("%dmer_EUR", k)
			assert.True(t, trainer.trainedModels[key])
		})
	}
}

// --- QualityTrainer success paths ---

// TestQualityTrainer_TrainQualityModel_Success exercises the full success path
// of TrainQualityModel by injecting a dataLoader with synthetic data points.
func TestQualityTrainer_TrainQualityModel_Success(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewQualityTrainer(config)
	trainer.dataLoader = func(platform string) ([]quality.GCQualityDataPoint, error) {
		return syntheticQualityData(platform), nil
	}

	err := trainer.TrainQualityModel("illumina")
	require.NoError(t, err)

	model, ok := trainer.models["illumina"]
	require.True(t, ok, "model should be stored under key 'illumina'")
	assert.NotNil(t, model)
}

// TestQualityTrainer_TrainQualityModel_MultiplePlatforms_WithData exercises
// multiple platforms through the success path.
func TestQualityTrainer_TrainQualityModel_MultiplePlatforms_WithData(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewQualityTrainer(config)
	trainer.dataLoader = func(platform string) ([]quality.GCQualityDataPoint, error) {
		return syntheticQualityData(platform), nil
	}

	for _, platform := range []string{"illumina", "pacbio", "nanopore"} {
		t.Run(platform, func(t *testing.T) {
			err := trainer.TrainQualityModel(platform)
			require.NoError(t, err)
			assert.NotNil(t, trainer.models[platform])
		})
	}
}

// TestQualityTrainer_TrainQualityModel_DataLoaderError verifies error propagation.
func TestQualityTrainer_TrainQualityModel_DataLoaderError(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewQualityTrainer(config)
	trainer.dataLoader = func(platform string) ([]quality.GCQualityDataPoint, error) {
		return nil, fmt.Errorf("injected quality data failure")
	}

	err := trainer.TrainQualityModel("illumina")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "injected quality data failure")
}

// --- PopulationTrainer success paths ---

// TestPopulationTrainer_TrainPopulationModel_Success exercises the full success
// path by injecting a paramsLoader with valid PopulationParameters.
func TestPopulationTrainer_TrainPopulationModel_Success(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewPopulationTrainer(config)
	trainer.paramsLoader = func(pop string) (*PopulationParameters, error) {
		return syntheticPopulationParams(), nil
	}

	err := trainer.TrainPopulationModel("YRI")
	require.NoError(t, err)

	simulator, ok := trainer.models["YRI"]
	require.True(t, ok, "simulator should be stored under key 'YRI'")
	assert.NotNil(t, simulator)
}

// TestPopulationTrainer_TrainPopulationModel_WithDemographicEvents tests the
// DemographicEvents loop inside TrainPopulationModel.
func TestPopulationTrainer_TrainPopulationModel_WithDemographicEvents(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewPopulationTrainer(config)
	trainer.paramsLoader = func(pop string) (*PopulationParameters, error) {
		return &PopulationParameters{
			EffectivePopSize:  10_000,
			MutationRate:      1.2e-8,
			RecombinationRate: 1.0e-8,
			DemographicEvents: []population.DemographicEvent{
				{
					Type:       "bottleneck",
					Time:       1000,
					Parameters: map[string]float64{"reduction": 0.1},
					Populations: []string{"YRI"},
				},
			},
		}, nil
	}

	err := trainer.TrainPopulationModel("YRI")
	require.NoError(t, err)

	sim, ok := trainer.models["YRI"]
	require.True(t, ok)
	events := sim.GetDemographicEvents()
	assert.Len(t, events, 1)
}

// TestPopulationTrainer_TrainPopulationModel_MultiplePopulations_WithData
// exercises multiple populations through the success path.
func TestPopulationTrainer_TrainPopulationModel_MultiplePopulations_WithData(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewPopulationTrainer(config)
	trainer.paramsLoader = func(pop string) (*PopulationParameters, error) {
		return syntheticPopulationParams(), nil
	}

	for _, pop := range []string{"YRI", "CEU", "CHB", "JPT"} {
		t.Run(pop, func(t *testing.T) {
			err := trainer.TrainPopulationModel(pop)
			require.NoError(t, err)
			assert.NotNil(t, trainer.models[pop])
		})
	}
}

// TestPopulationTrainer_TrainPopulationModel_ParamsLoaderError tests error propagation.
func TestPopulationTrainer_TrainPopulationModel_ParamsLoaderError(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewPopulationTrainer(config)
	trainer.paramsLoader = func(pop string) (*PopulationParameters, error) {
		return nil, fmt.Errorf("injected params failure")
	}

	err := trainer.TrainPopulationModel("YRI")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "injected params failure")
}

// TestQualityTrainer_TrainQualityModel_EmptyDataError verifies that when the
// dataLoader succeeds but returns an empty slice, TrainFromRealData returns an
// error — covering the "if err := model.TrainFromRealData" error path.
func TestQualityTrainer_TrainQualityModel_EmptyDataError(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewQualityTrainer(config)
	trainer.dataLoader = func(platform string) ([]quality.GCQualityDataPoint, error) {
		return []quality.GCQualityDataPoint{}, nil // success but empty
	}

	err := trainer.TrainQualityModel("illumina")
	require.Error(t, err)
	// TrainFromRealData returns "no training data provided" for empty slice
}

// TestHMMTrainer_TrainHMMModel_SaveModelError verifies that when OutputPath is
// unwritable, saveModel returns an error — covering the MkdirAll failure path.
func TestHMMTrainer_TrainHMMModel_SaveModelError(t *testing.T) {
	dir := makeSyntheticFASTADir(t)
	// /dev/null is a device file; MkdirAll cannot create a subdirectory inside it
	config := &TrainingConfig{OutputPath: "/dev/null/hmm_models_cannot_exist"}
	trainer := NewHMMTrainer(config)
	trainer.genomeLoader = func(pop string) (*reference.ReferenceGenome, error) {
		return reference.NewReferenceGenome(dir), nil
	}

	err := trainer.TrainHMMModel(4, "YRI")
	require.Error(t, err)
	assert.Contains(t, err.Error(), "create output directory")
}

// TestPopulationTrainer_TrainPopulationModel_InvalidParams verifies that
// invalid PopulationParameters (zero pop size) cause ValidateParameters to fail.
func TestPopulationTrainer_TrainPopulationModel_InvalidParams(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewPopulationTrainer(config)
	trainer.paramsLoader = func(pop string) (*PopulationParameters, error) {
		return &PopulationParameters{
			EffectivePopSize:  0, // invalid
			MutationRate:      1.2e-8,
			RecombinationRate: 1.0e-8,
		}, nil
	}

	err := trainer.TrainPopulationModel("YRI")
	require.Error(t, err) // ValidateParameters rejects EffectivePopSize == 0
}
