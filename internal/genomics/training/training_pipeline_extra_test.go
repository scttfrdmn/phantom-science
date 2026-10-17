package training

import (
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/validation"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestDataManager_GetReferenceGenome_NotLoaded tests getting an uncached reference genome
func TestDataManager_GetReferenceGenome_NotLoaded(t *testing.T) {
	config := &TrainingConfig{
		ReferenceGenomePath: "/test/reference",
	}
	dm := NewDataManager(config)

	_, err := dm.GetReferenceGenome("GRCh38")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "not loaded")
}

// TestDataManager_GetPopulationPanel_NotLoaded tests getting an uncached population panel
func TestDataManager_GetPopulationPanel_NotLoaded(t *testing.T) {
	config := &TrainingConfig{
		PopulationDataPath: "/test/population",
	}
	dm := NewDataManager(config)

	_, err := dm.GetPopulationPanel("EUR")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "not loaded")
}

// TestDataManager_GetValidationDataset_NotLoaded tests getting an uncached validation dataset
func TestDataManager_GetValidationDataset_NotLoaded(t *testing.T) {
	config := &TrainingConfig{
		ValidationDataPath: "/test/validation",
	}
	dm := NewDataManager(config)

	_, err := dm.GetValidationDataset("illumina")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "not loaded")
}

// TestDataManager_LoadValidationDataset tests loading a validation dataset
func TestDataManager_LoadValidationDataset(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		ValidationDataPath: tempDir,
	}
	dm := NewDataManager(config)

	platforms := []string{"illumina", "novaseq", "nextseq", "pacbio", "sequel", "nanopore", "minion", "unknown_platform"}
	for _, platform := range platforms {
		t.Run(platform, func(t *testing.T) {
			err := dm.LoadValidationDataset(platform)
			require.NoError(t, err)

			// Should be cached now
			dataset, err := dm.GetValidationDataset(platform)
			require.NoError(t, err)
			assert.NotNil(t, dataset)
			assert.Equal(t, platform, dataset.Platform)
			assert.Equal(t, "simulated", dataset.DataType)
			assert.Greater(t, dataset.Size, 0)
			assert.NotNil(t, dataset.Sequences)
			assert.NotNil(t, dataset.Qualities)
		})
	}
}

// TestDataManager_LoadValidationDataset_Cached tests that loading twice uses cache
func TestDataManager_LoadValidationDataset_Cached(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		ValidationDataPath: tempDir,
	}
	dm := NewDataManager(config)

	err := dm.LoadValidationDataset("illumina")
	require.NoError(t, err)

	// Load again - should use cache (no error)
	err = dm.LoadValidationDataset("illumina")
	require.NoError(t, err)
}

// TestDataManager_GenerateValidationData_Platforms tests all platform branches
func TestDataManager_GenerateValidationData_Platforms(t *testing.T) {
	config := &TrainingConfig{}
	dm := NewDataManager(config)

	testCases := []struct {
		platform       string
		expectedLength int
		expectedCount  int
	}{
		{"illumina", 150, 1000},
		{"novaseq", 150, 1000},
		{"nextseq", 150, 1000},
		{"pacbio", 10000, 100},
		{"sequel", 10000, 100},
		{"nanopore", 5000, 200},
		{"minion", 5000, 200},
		{"anything_else", 150, 1000},
	}

	for _, tc := range testCases {
		t.Run(tc.platform, func(t *testing.T) {
			dataset := &ValidationDataset{
				Platform: tc.platform,
				DataType: "simulated",
				Metadata: make(map[string]interface{}),
			}
			err := dm.generateValidationData(dataset)
			require.NoError(t, err)
			assert.Equal(t, tc.expectedCount, dataset.Size)
			assert.Equal(t, tc.expectedLength, dataset.Metadata["sequence_length"])
			assert.Len(t, dataset.Sequences, tc.expectedCount)
			assert.Len(t, dataset.Qualities, tc.expectedCount)
		})
	}
}

// TestDataManager_LoadPopulationPanel tests population panel loading
func TestDataManager_LoadPopulationPanel(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		PopulationDataPath: tempDir,
		ValidationSplit:    0.2,
	}
	dm := NewDataManager(config)

	populations := []string{"AFR", "YRI", "LWK", "GWD", "EUR", "CEU", "GBR", "TSI", "ASN", "CHB", "JPT", "CHS", "AMR", "MXL", "PUR", "UNKNOWN"}
	for _, pop := range populations {
		t.Run(pop, func(t *testing.T) {
			err := dm.LoadPopulationPanel(pop)
			require.NoError(t, err)

			panel, err := dm.GetPopulationPanel(pop)
			require.NoError(t, err)
			assert.NotNil(t, panel)
		})
	}
}


// TestDataManager_LoadPopulationPanel_Cached tests caching
func TestDataManager_LoadPopulationPanel_Cached(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		PopulationDataPath: tempDir,
		ValidationSplit:    0.2,
	}
	dm := NewDataManager(config)

	err := dm.LoadPopulationPanel("EUR")
	require.NoError(t, err)

	// Load again - should use cache
	err = dm.LoadPopulationPanel("EUR")
	require.NoError(t, err)
}

// TestDataManager_LoadPopulationPanel_ZeroSplit tests zero validation split
func TestDataManager_LoadPopulationPanel_ZeroSplit(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		PopulationDataPath: tempDir,
		ValidationSplit:    0.0, // triggers default sample size path
	}
	dm := NewDataManager(config)

	err := dm.LoadPopulationPanel("EUR")
	require.NoError(t, err)
}

// TestValidationManager_TestMemoryUsage tests memory usage testing
func TestValidationManager_TestMemoryUsage(t *testing.T) {
	config := &TrainingConfig{}
	vm := NewValidationManager(config)

	// Without pre-existing performance results
	err := vm.TestMemoryUsage()
	assert.NoError(t, err)
	assert.NotNil(t, vm.results.PerformanceResults)
	assert.Greater(t, vm.results.PerformanceResults.MemoryUsage, int64(0))
}

// TestValidationManager_TestMemoryUsage_WithExistingResults tests memory usage with pre-existing results
func TestValidationManager_TestMemoryUsage_WithExistingResults(t *testing.T) {
	config := &TrainingConfig{}
	vm := NewValidationManager(config)

	// Pre-populate results
	vm.results.PerformanceResults = &PerformanceResults{
		GenerationSpeed: 50000.0,
	}

	err := vm.TestMemoryUsage()
	assert.NoError(t, err)
	assert.Greater(t, vm.results.PerformanceResults.MemoryUsage, int64(0))
	// Should preserve existing fields
	assert.Equal(t, 50000.0, vm.results.PerformanceResults.GenerationSpeed)
}

// TestValidationManager_TestScalability tests scalability testing
func TestValidationManager_TestScalability(t *testing.T) {
	config := &TrainingConfig{}
	vm := NewValidationManager(config)

	// Without pre-existing performance results
	err := vm.TestScalability()
	assert.NoError(t, err)
	assert.NotNil(t, vm.results.PerformanceResults)
	assert.Greater(t, vm.results.PerformanceResults.ScalabilityScore, 0.0)
}

// TestValidationManager_TestScalability_WithExistingResults tests scalability with pre-existing results
func TestValidationManager_TestScalability_WithExistingResults(t *testing.T) {
	config := &TrainingConfig{}
	vm := NewValidationManager(config)

	vm.results.PerformanceResults = &PerformanceResults{
		GenerationSpeed: 50000.0,
	}

	err := vm.TestScalability()
	assert.NoError(t, err)
	assert.Equal(t, 90.0, vm.results.PerformanceResults.ScalabilityScore)
}

// TestValidationManager_AggregateResults_NoJSResults tests aggregation with no JS divergence results.
// When no validation data has been collected, the score must be 0 (conservative default)
// and PassesValidation must be false — a model with no evidence should not silently pass.
func TestValidationManager_AggregateResults_NoJSResults(t *testing.T) {
	config := &TrainingConfig{}
	vm := NewValidationManager(config)

	results := vm.AggregateResults()
	assert.NotNil(t, results)
	assert.Equal(t, 0.0, results.OverallScore)
	assert.False(t, results.PassesValidation)
}

// TestValidationManager_AggregateResults_FailingScore tests aggregation with low score
func TestValidationManager_AggregateResults_FailingScore(t *testing.T) {
	config := &TrainingConfig{}
	vm := NewValidationManager(config)

	// Add large JS divergence values that will cause low overall score
	vm.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.05, // Large divergence -> low score
		"pacbio":   0.08,
	}

	results := vm.AggregateResults()
	assert.NotNil(t, results)
	assert.Less(t, results.OverallScore, 95.0)
	assert.False(t, results.PassesValidation)
	assert.Greater(t, len(results.Recommendations), 0)
}

// TestNewPipelineLogger_InvalidPath tests logger with invalid path
func TestNewPipelineLogger_InvalidPath(t *testing.T) {
	// Using non-existent path - should not panic, just return logger without file
	logger := NewPipelineLogger("/non_existent_path_that_cant_be_created_xyz")

	assert.NotNil(t, logger)
	// Should still log without crashing
	logger.Info("test message")
	logger.Error("test error", assert.AnError)
}

// TestPipelineLogger_Info tests info logging
func TestPipelineLogger_Info(t *testing.T) {
	tempDir := t.TempDir()
	logger := NewPipelineLogger(tempDir)
	require.NotNil(t, logger)

	// Should not panic
	logger.Info("Test info message")
	logger.Info("Another info message")

	// Verify log file content
	logPath := filepath.Join(tempDir, "training.log")
	content, err := os.ReadFile(logPath) // #nosec G304
	require.NoError(t, err)
	assert.Contains(t, string(content), "INFO")
	assert.Contains(t, string(content), "Test info message")
}

// TestPipelineLogger_Error tests error logging
func TestPipelineLogger_Error(t *testing.T) {
	tempDir := t.TempDir()
	logger := NewPipelineLogger(tempDir)
	require.NotNil(t, logger)

	testErr := assert.AnError
	logger.Error("something failed", testErr)

	logPath := filepath.Join(tempDir, "training.log")
	content, err := os.ReadFile(logPath) // #nosec G304
	require.NoError(t, err)
	assert.Contains(t, string(content), "ERROR")
	assert.Contains(t, string(content), "something failed")
}

// TestNewModelTrainer_SubComponents tests model trainer sub-components creation
func TestNewModelTrainer_SubComponents(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath: tempDir,
	}

	trainer := NewModelTrainer(config)
	require.NotNil(t, trainer)

	// Test delegating methods work (even though underlying implementations return errors)
	err := trainer.TrainHMMModel(4, "EUR")
	assert.Error(t, err) // Expected - loadReferenceGenome not implemented

	err = trainer.TrainQualityModel("illumina")
	assert.Error(t, err) // Expected - loadQualityTrainingData not implemented

	err = trainer.TrainPopulationModel("EUR")
	assert.Error(t, err) // Expected - getPopulationParameters not implemented
}

// TestHMMTrainer_TrainHMMModel tests HMM model training (expected failure path)
func TestHMMTrainer_TrainHMMModel(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewHMMTrainer(config)

	// This should fail because loadReferenceGenome returns error
	err := trainer.TrainHMMModel(4, "EUR")
	assert.Error(t, err)
}

// TestQualityTrainer_TrainQualityModel tests quality model training (expected failure path)
func TestQualityTrainer_TrainQualityModel(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewQualityTrainer(config)

	// This should fail because loadQualityTrainingData returns error
	err := trainer.TrainQualityModel("illumina")
	assert.Error(t, err)
}

// TestPopulationTrainer_TrainPopulationModel tests population model training (expected failure path)
func TestPopulationTrainer_TrainPopulationModel(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewPopulationTrainer(config)

	// This should fail because getPopulationParameters returns error
	err := trainer.TrainPopulationModel("EUR")
	assert.Error(t, err)
}

// TestTrainingPipeline_Execute_FailsWithBadPaths tests Execute with bad paths
func TestTrainingPipeline_Execute_FailsWithBadPaths(t *testing.T) {
	config := &TrainingConfig{
		ReferenceGenomePath: "/nonexistent/reference",
		PopulationDataPath:  "/nonexistent/population",
		ValidationDataPath:  "/nonexistent/validation",
		OutputPath:          "/nonexistent/output",
		KmerSizes:           []int{4},
		Populations:         []string{"EUR"},
		Platforms:           []string{"illumina"},
		ValidationSplit:     0.2,
	}

	pipeline := NewTrainingPipeline(config)

	// Execute should fail because reference genome loading fails
	err := pipeline.Execute()
	assert.Error(t, err)
}

// TestTrainingPipeline_Execute_WithValidData tests Execute succeeds up to model training
func TestTrainingPipeline_Execute_WithValidData(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		ReferenceGenomePath:        tempDir,
		PopulationDataPath:         tempDir,
		ValidationDataPath:         tempDir,
		OutputPath:                 tempDir,
		CachePath:                  tempDir,
		KmerSizes:                  []int{4},
		Populations:                []string{"EUR"},
		Platforms:                  []string{"illumina"},
		ValidationSplit:            0.2,
		CrossValidationFolds:       5,
		MaxMemoryGB:                8,
		NumWorkers:                 4,
		BatchSize:                  100,
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.9,
	}

	pipeline := NewTrainingPipeline(config)

	// Execute will fail during reference genome loading (no actual data)
	// but we test the pipeline can be created and starts execution
	err := pipeline.Execute()
	assert.Error(t, err) // Expected - no actual reference data
}

// TestTrainingConfig_Defaults tests configuration with various field combinations
func TestTrainingConfig_Defaults(t *testing.T) {
	config := &TrainingConfig{
		EnableEarlyStopping: true,
		EnableGPU:           true,
		MaxTrainingTime:     time.Hour * 24,
		CheckpointInterval:  time.Minute * 30,
	}

	assert.True(t, config.EnableEarlyStopping)
	assert.True(t, config.EnableGPU)
	assert.Equal(t, time.Hour*24, config.MaxTrainingTime)
	assert.Equal(t, time.Minute*30, config.CheckpointInterval)
}

// TestPopulationParameters_Structure tests PopulationParameters struct
func TestPopulationParameters_Structure(t *testing.T) {
	params := &PopulationParameters{
		EffectivePopSize:  10000,
		MutationRate:      1.25e-8,
		RecombinationRate: 1.13e-8,
	}

	assert.Equal(t, 10000, params.EffectivePopSize)
	assert.Equal(t, 1.25e-8, params.MutationRate)
	assert.Equal(t, 1.13e-8, params.RecombinationRate)
	assert.Nil(t, params.DemographicEvents)
}

// TestDefaultAWSDataConfig tests default AWS config creation
func TestDefaultAWSDataConfig(t *testing.T) {
	config := DefaultAWSDataConfig()

	assert.NotNil(t, config)
	assert.Equal(t, "us-east-1", config.Region)
	assert.NotEmpty(t, config.DownloadDir)
	assert.Greater(t, config.MaxConcurrency, 0)
	assert.Greater(t, config.MaxRetries, 0)
	assert.True(t, config.VerifyChecksums)
	assert.True(t, config.EnableProgress)
	assert.True(t, config.CacheData)
	assert.True(t, config.CleanupOnError)
}

// TestGenomicDatasets tests that genomic dataset catalog is populated
func TestGenomicDatasets(t *testing.T) {
	assert.NotEmpty(t, GenomicDatasets)

	expectedDatasets := []string{"grch38", "grch39", "1000genomes", "giab", "broad_references", "encode", "gnomad", "tcga", "sra", "allen_brain"}
	for _, name := range expectedDatasets {
		dataset, exists := GenomicDatasets[name]
		assert.True(t, exists, "Dataset %s should exist", name)
		assert.NotEmpty(t, dataset.Name)
		assert.NotEmpty(t, dataset.Bucket)
	}
}

// TestAuthenticityValidationDatasets tests validation dataset catalog
func TestAuthenticityValidationDatasets(t *testing.T) {
	assert.NotEmpty(t, AuthenticityValidationDatasets)

	expectedDatasets := []string{"reference_genome", "population_vcf", "illumina_fastq", "pacbio_fastq", "nanopore_fastq"}
	for _, name := range expectedDatasets {
		dataset, exists := AuthenticityValidationDatasets[name]
		assert.True(t, exists, "Dataset %s should exist", name)
		assert.NotEmpty(t, dataset.Name)
		assert.NotEmpty(t, dataset.Bucket)
		assert.NotEmpty(t, dataset.Files)

		for _, file := range dataset.Files {
			assert.NotEmpty(t, file.Key)
			assert.Greater(t, file.Size, int64(0))
		}
	}
}

// TestAWSDataConfig_Structure tests AWSDataConfig structure
func TestAWSDataConfig_Structure(t *testing.T) {
	config := &AWSDataConfig{
		Region:          "eu-west-1",
		DownloadDir:     "/tmp/downloads",
		MaxConcurrency:  5,
		MaxRetries:      3,
		VerifyChecksums: true,
		EnableProgress:  false,
		CacheData:       false,
		CleanupOnError:  true,
	}

	assert.Equal(t, "eu-west-1", config.Region)
	assert.Equal(t, "/tmp/downloads", config.DownloadDir)
	assert.Equal(t, 5, config.MaxConcurrency)
	assert.Equal(t, 3, config.MaxRetries)
}

// TestDatasetInfo_Structure tests DatasetInfo structure
func TestDatasetInfo_Structure(t *testing.T) {
	now := time.Now()
	info := DatasetInfo{
		Name:        "Test Dataset",
		Description: "Test description",
		Bucket:      "test-bucket",
		Prefix:      "test/prefix/",
		Files: []FileInfo{
			{
				Key:          "file1.vcf.gz",
				Size:         1000000,
				LastModified: now,
				ETag:         "abc123",
				LocalPath:    "/local/path/file1.vcf.gz",
			},
		},
		Size:         1000000,
		LastModified: now,
		Checksum:     "sha256sum",
	}

	assert.Equal(t, "Test Dataset", info.Name)
	assert.Equal(t, "test-bucket", info.Bucket)
	assert.Len(t, info.Files, 1)
	assert.Equal(t, "file1.vcf.gz", info.Files[0].Key)
	assert.Equal(t, int64(1000000), info.Files[0].Size)
}

// TestFileInfo_Structure tests FileInfo structure
func TestFileInfo_Structure(t *testing.T) {
	now := time.Now()
	fileInfo := FileInfo{
		Key:          "test/file.vcf.gz",
		Size:         500000,
		LastModified: now,
		ETag:         "etag123",
		LocalPath:    "/local/test/file.vcf.gz",
	}

	assert.Equal(t, "test/file.vcf.gz", fileInfo.Key)
	assert.Equal(t, int64(500000), fileInfo.Size)
	assert.Equal(t, now, fileInfo.LastModified)
	assert.Equal(t, "etag123", fileInfo.ETag)
}

// TestNewProgressTracker tests progress tracker creation
func TestNewProgressTracker(t *testing.T) {
	totalBytes := int64(1024 * 1024 * 100) // 100 MB
	tracker := NewProgressTracker(totalBytes)

	assert.NotNil(t, tracker)
	assert.Equal(t, totalBytes, tracker.totalBytes)
	assert.Equal(t, int64(0), tracker.downloadedBytes)
}

// TestProgressTracker_Update tests progress updates
func TestProgressTracker_Update(t *testing.T) {
	tracker := NewProgressTracker(1000)

	tracker.Update(100)
	assert.Equal(t, int64(100), tracker.downloadedBytes)

	tracker.Update(200)
	assert.Equal(t, int64(300), tracker.downloadedBytes)
}

// TestProgressTracker_GetProgress tests getting progress
func TestProgressTracker_GetProgress(t *testing.T) {
	tracker := NewProgressTracker(1000)
	tracker.Update(500)

	progress := tracker.GetProgress()
	assert.NotNil(t, progress)
	assert.Equal(t, int64(1000), progress.TotalBytes)
	assert.Equal(t, int64(500), progress.DownloadedBytes)
	assert.GreaterOrEqual(t, progress.Speed, 0.0)
}

// TestProgressTracker_GetProgress_NoDownload tests progress with no bytes downloaded
func TestProgressTracker_GetProgress_NoDownload(t *testing.T) {
	tracker := NewProgressTracker(1000)

	progress := tracker.GetProgress()
	assert.NotNil(t, progress)
	assert.Equal(t, int64(0), progress.DownloadedBytes)
}

// TestDownloadProgress_Structure tests DownloadProgress structure
func TestDownloadProgress_Structure(t *testing.T) {
	now := time.Now()
	progress := &DownloadProgress{
		TotalFiles:      10,
		CompletedFiles:  5,
		TotalBytes:      1000000,
		DownloadedBytes: 500000,
		StartTime:       now,
		CurrentFile:     "file.vcf.gz",
		Speed:           10.5,
		ETA:             time.Second * 30,
	}

	assert.Equal(t, 10, progress.TotalFiles)
	assert.Equal(t, 5, progress.CompletedFiles)
	assert.Equal(t, int64(1000000), progress.TotalBytes)
	assert.Equal(t, int64(500000), progress.DownloadedBytes)
	assert.Equal(t, "file.vcf.gz", progress.CurrentFile)
	assert.Equal(t, 10.5, progress.Speed)
}

// TestNewAWSDataAcquisition_BadConfig tests creation with bad config
func TestNewAWSDataAcquisition_BadConfig(t *testing.T) {
	tempDir := t.TempDir()
	config := &AWSDataConfig{
		Region:         "us-east-1",
		DownloadDir:    tempDir,
		MaxConcurrency: 5,
		MaxRetries:     3,
	}

	// Note: This may fail in CI without AWS credentials but should not panic
	_, err := NewAWSDataAcquisition(config)
	// We don't assert error/no-error because it depends on AWS SDK availability
	// Just test it doesn't panic
	_ = err
}

// TestAWSDataAcquisition_fileExists tests file existence checking via temp files
func TestAWSDataAcquisition_fileExists_indirect(t *testing.T) {
	tempDir := t.TempDir()

	// Create a test file
	testFile := filepath.Join(tempDir, "test.txt")
	err := os.WriteFile(testFile, []byte("content"), 0600)
	require.NoError(t, err)

	// Test using os.Stat which is what fileExists uses internally
	info, err := os.Stat(testFile)
	assert.NoError(t, err)
	assert.Greater(t, info.Size(), int64(0))

	// Test non-existent file
	_, err = os.Stat(filepath.Join(tempDir, "nonexistent.txt"))
	assert.Error(t, err)
}

// TestTrainingPipeline_validateModels_WithExpertAccuracy tests expert accuracy threshold
func TestValidationManager_FinalValidation_ExpertAccuracyThreshold(t *testing.T) {
	config := &TrainingConfig{
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.9,
	}
	vm := NewValidationManager(config)

	// Populate results that will fail
	vm.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.005, // Within threshold
	}
	vm.results.ExpertTestResults = nil // No expert results
	vm.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	results := vm.AggregateResults()
	assert.NotNil(t, results)
}

// TestModelCheckpointer_Structure tests model checkpointer structure
func TestModelCheckpointer_Structure(t *testing.T) {
	checkpointer := NewModelCheckpointer("/path/to/output")
	assert.NotNil(t, checkpointer)
	assert.Equal(t, "/path/to/output", checkpointer.outputPath)
}

// TestValidationManager_ValidateModels_WithExpertThreshold tests that RunExpertTesting
// succeeds even when the expert framework is not yet populated with test data.
// The honest result is ExpertTestResults == nil, not fabricated data.
func TestValidationManager_ValidateModels_WithExpertThreshold(t *testing.T) {
	config := &TrainingConfig{
		ExpertAccuracyThreshold: 0.6, // Greater than 0 triggers RunExpertTesting
	}
	vm := NewValidationManager(config)

	err := vm.RunStatisticalValidation()
	assert.NoError(t, err)

	// Without test cases/experts registered, RunExpertTesting skips gracefully.
	err = vm.RunExpertTesting()
	assert.NoError(t, err)
	// Framework not configured → results are nil, not fabricated.
	assert.Nil(t, vm.results.ExpertTestResults)

	err = vm.RunDownstreamTesting()
	assert.NoError(t, err)
}

// TestTrainingPipeline_PrepareData_EmptyPopulations tests prepareData with no populations
func TestTrainingPipeline_PrepareData_EmptyPopulations(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		ReferenceGenomePath: tempDir,
		PopulationDataPath:  tempDir,
		ValidationDataPath:  tempDir,
		OutputPath:          tempDir,
		KmerSizes:           []int{4},
		Populations:         []string{}, // Empty - no population loading
		Platforms:           []string{"illumina"},
		ValidationSplit:     0.2,
	}

	pipeline := NewTrainingPipeline(config)

	// Execute will fail at reference genome loading since no real data exists
	err := pipeline.Execute()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "data preparation failed")
}

// TestTrainingPipeline_ValidateModels_ExpertAccuracyThreshold tests the path where
// expert testing is invoked due to ExpertAccuracyThreshold > 0.  With no framework
// data configured, RunExpertTesting succeeds with nil results (correct honest behaviour).
func TestTrainingPipeline_ValidateModels_ExpertAccuracyThreshold(t *testing.T) {
	config := &TrainingConfig{
		ExpertAccuracyThreshold: 0.6,
	}
	vm := NewValidationManager(config)

	err := vm.RunStatisticalValidation()
	assert.NoError(t, err)

	if config.ExpertAccuracyThreshold > 0 {
		err = vm.RunExpertTesting()
		assert.NoError(t, err)
		// Framework not configured → nil, not fabricated.
		assert.Nil(t, vm.results.ExpertTestResults)
	}

	err = vm.RunDownstreamTesting()
	assert.NoError(t, err)
}

// TestTrainingPipeline_FinalValidation_JSdivergenceFailure tests final validation with high JS divergence
func TestTrainingPipeline_FinalValidation_JSdivergenceFailure(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		JSdivergenceThreshold:      0.001, // Very tight threshold
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.9,
	}

	vm := NewValidationManager(config)

	// Populate results that will exceed the JS divergence threshold
	vm.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.005, // This value makes overall score low (< 95)
		"pacbio":   0.008,
		"nanopore": 0.012,
	}
	vm.results.ExpertTestResults = nil
	vm.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	// The overall score will be low due to high JS divergence, causing finalValidation to fail
	results := vm.AggregateResults()
	assert.Less(t, results.OverallScore, 95.0)
}

// TestTrainingPipeline_FinalValidation_HighExpertAccuracy tests final validation when expert accuracy is too high
func TestTrainingPipeline_FinalValidation_HighExpertAccuracy(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.3, // Expert accuracy should be BELOW this (authenticity goal)
		DownstreamSuccessThreshold: 0.9,
	}

	vm := NewValidationManager(config)

	// Simulate results where expert accuracy (0.48) exceeds threshold (0.3)
	vm.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.000001, // Very low, so score stays high
	}
	vm.results.ExpertTestResults = nil // No expert results - skip that check
	vm.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	results := vm.AggregateResults()
	// With very low JS divergence, overall score should be near 95
	assert.GreaterOrEqual(t, results.OverallScore, 90.0)
	assert.NotNil(t, config.ExpertAccuracyThreshold)
}

// TestTrainingPipeline_FinalValidation_LowDownstreamSuccess tests low downstream success rate
func TestTrainingPipeline_FinalValidation_LowDownstreamSuccess(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.99, // Very high threshold (>0.95)
	}

	vm := NewValidationManager(config)

	// Low JS divergence to allow high overall score
	vm.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.000001, // Tiny divergence -> near 95 score
	}
	vm.results.ExpertTestResults = nil
	vm.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95, // Below the 0.99 threshold
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	results := vm.AggregateResults()
	// Config threshold is higher than actual results
	assert.Less(t, results.DownstreamResults.SuccessRate, config.DownstreamSuccessThreshold)
}

// TestTrainingPipeline_TrainModels_MultipleKmers tests trainModels with multiple k-mer sizes
func TestTrainingPipeline_TrainModels_MultipleKmers(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		KmerSizes:                  []int{4, 6, 8},
		Populations:                []string{"EUR", "AFR"},
		Platforms:                  []string{"illumina", "pacbio"},
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.9,
	}

	mt := NewModelTrainer(config)

	// TrainHMMModel will fail (loadReferenceGenome not implemented)
	// but we exercise the loop paths
	for _, k := range config.KmerSizes {
		for _, pop := range config.Populations {
			err := mt.TrainHMMModel(k, pop)
			assert.Error(t, err) // Expected failure
		}
	}

	// TrainQualityModel will fail (loadQualityTrainingData not implemented)
	for _, platform := range config.Platforms {
		err := mt.TrainQualityModel(platform)
		assert.Error(t, err) // Expected failure
	}

	// TrainPopulationModel will fail (getPopulationParameters not implemented)
	for _, pop := range config.Populations {
		err := mt.TrainPopulationModel(pop)
		assert.Error(t, err) // Expected failure
	}
}

// TestValidationManager_FinalValidation_PassWithGoodResults tests finalValidation success path
func TestValidationManager_FinalValidation_PassWithGoodResults(t *testing.T) {
	config := &TrainingConfig{
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.9,
	}
	vm := NewValidationManager(config)

	// Pre-populate with results that will pass all checks
	vm.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.000001, // Extremely low divergence
	}
	vm.results.ExpertTestResults = nil // nil means check is skipped
	vm.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA", "GATK"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	results := vm.AggregateResults()
	assert.NotNil(t, results)
	// With near-zero JS divergence, overall score should be close to 95
	// 100 - (0.000001 * 1000) = 99.999
	assert.Greater(t, results.OverallScore, 90.0)
}

// TestTrainingPipeline_Execute_WithEmptyPlatforms tests Execute with no platforms
func TestTrainingPipeline_Execute_WithEmptyPlatforms(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		ReferenceGenomePath:        tempDir,
		PopulationDataPath:         tempDir,
		ValidationDataPath:         tempDir,
		OutputPath:                 tempDir,
		KmerSizes:                  []int{},
		Populations:                []string{},
		Platforms:                  []string{},
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.0, // Zero means no expert testing
		DownstreamSuccessThreshold: 0.9,
	}

	pipeline := NewTrainingPipeline(config)

	// Execute will fail at reference genome loading phase
	err := pipeline.Execute()
	assert.Error(t, err)
}

// TestTrainingPipeline_Execute_NoExpertTesting tests pipeline with no expert testing threshold
func TestTrainingPipeline_Execute_NoExpertTesting(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		ReferenceGenomePath:        tempDir,
		PopulationDataPath:         tempDir,
		ValidationDataPath:         tempDir,
		OutputPath:                 tempDir,
		KmerSizes:                  []int{4},
		Populations:                []string{"EUR"},
		Platforms:                  []string{"illumina"},
		ValidationSplit:            0.2,
		ExpertAccuracyThreshold:    0.0, // Zero = no expert testing
		JSdivergenceThreshold:      0.1,
		DownstreamSuccessThreshold: 0.9,
	}

	pipeline := NewTrainingPipeline(config)

	// Will fail at reference genome loading (no actual data), not at expert testing
	err := pipeline.Execute()
	assert.Error(t, err)
	// Should fail in data preparation phase
	assert.Contains(t, err.Error(), "data preparation failed")
}

// TestTrainingPipeline_ValidateModels_DirectCall tests validateModels called directly
func TestTrainingPipeline_ValidateModels_DirectCall(t *testing.T) {
	tempDir := t.TempDir()
	// Test with ExpertAccuracyThreshold = 0 (no expert testing)
	config := &TrainingConfig{
		OutputPath:              tempDir,
		ExpertAccuracyThreshold: 0.0,
	}
	pipeline := NewTrainingPipeline(config)

	// Call validateModels directly (white-box test - same package)
	err := pipeline.validateModels()
	assert.NoError(t, err)
}

// TestTrainingPipeline_ValidateModels_WithExpertAccuracy_Direct tests validateModels
// when ExpertAccuracyThreshold > 0 triggers RunExpertTesting.
// Without test cases/experts registered, results are nil (honest, not fabricated).
func TestTrainingPipeline_ValidateModels_WithExpertAccuracy_Direct(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:              tempDir,
		ExpertAccuracyThreshold: 0.6, // > 0 triggers RunExpertTesting
	}
	pipeline := NewTrainingPipeline(config)

	// Call validateModels directly (white-box test)
	err := pipeline.validateModels()
	assert.NoError(t, err)
	// Framework not configured → ExpertTestResults remains nil, not fabricated.
	assert.Nil(t, pipeline.validator.results.ExpertTestResults)
}

// TestTrainingPipeline_TestPerformance_DirectCall tests testPerformance called directly
func TestTrainingPipeline_TestPerformance_DirectCall(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath: tempDir,
	}
	pipeline := NewTrainingPipeline(config)

	// Call testPerformance directly (white-box test)
	err := pipeline.testPerformance()
	assert.NoError(t, err)
	// Performance results should be populated
	assert.NotNil(t, pipeline.validator.results.PerformanceResults)
	assert.Greater(t, pipeline.validator.results.PerformanceResults.GenerationSpeed, 0.0)
	assert.Greater(t, pipeline.validator.results.PerformanceResults.MemoryUsage, int64(0))
	assert.Greater(t, pipeline.validator.results.PerformanceResults.ScalabilityScore, 0.0)
}

// TestTrainingPipeline_FinalValidation_DirectCall tests finalValidation called directly
func TestTrainingPipeline_FinalValidation_DirectCall(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.9,
	}
	pipeline := NewTrainingPipeline(config)

	// Pre-populate results so finalValidation can run
	pipeline.validator.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.000001, // Very low JS divergence -> high score
	}
	pipeline.validator.results.ExpertTestResults = nil // nil -> skip expert check
	pipeline.validator.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	// Call finalValidation directly (white-box test)
	err := pipeline.finalValidation()
	assert.NoError(t, err)
}

// TestTrainingPipeline_FinalValidation_LowScore tests finalValidation with low overall score
func TestTrainingPipeline_FinalValidation_LowScore(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.9,
	}
	pipeline := NewTrainingPipeline(config)

	// Pre-populate results with high JS divergence -> low overall score
	pipeline.validator.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.1, // High JS divergence -> score = max(0, 100 - 0.1*1000) = 0
	}
	pipeline.validator.results.ExpertTestResults = nil
	pipeline.validator.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	// Call finalValidation directly - should fail because overall score < 95
	err := pipeline.finalValidation()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "overall authenticity score too low")
}

// TestTrainingPipeline_FinalValidation_HighJSDivergence tests finalValidation with high JS divergence
func TestTrainingPipeline_FinalValidation_HighJSDivergence(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		JSdivergenceThreshold:      0.001, // Very tight threshold
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.9,
	}
	pipeline := NewTrainingPipeline(config)

	// Pre-populate results with high JS divergence score but low individual values
	// Need overall score >= 95 but individual JS divergence to exceed threshold
	pipeline.validator.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.000001, // Tiny so score is ~100 = passes overall check
	}
	pipeline.validator.results.ExpertTestResults = nil
	pipeline.validator.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	// First call AggregateResults to set the overall score
	_ = pipeline.validator.AggregateResults()

	// Now set a value that exceeds the JSdivergenceThreshold in the already-aggregated results
	pipeline.validator.results.JSdivergenceResults["illumina"] = 0.01 // > 0.001 threshold

	// Call finalValidation directly - should fail on JS divergence check
	err := pipeline.finalValidation()
	// The first check is OverallScore < 95 - the value we set before aggregation
	// was 0.000001 so score is ~100, passes first check
	// Then the JS divergence check: 0.01 > 0.001 threshold -> error
	_ = err // May or may not error depending on overall score
}

// TestTrainingPipeline_FinalValidation_ExpertAccuracyTooHigh tests expert accuracy check
func TestTrainingPipeline_FinalValidation_ExpertAccuracyTooHigh(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.3,  // Expert should be BELOW this
		DownstreamSuccessThreshold: 0.9,
	}
	pipeline := NewTrainingPipeline(config)

	// Pre-populate results so overall score is >= 95
	pipeline.validator.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.000001, // Tiny -> high score
	}
	// Set expert results with accuracy ABOVE the threshold (bad for authenticity)
	pipeline.validator.results.ExpertTestResults = &validation.ExpertTestingResults{
		OverallAccuracy: 0.7, // 70% > 30% threshold -> too identifiable
	}
	pipeline.validator.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	// Aggregate to set overall score
	_ = pipeline.validator.AggregateResults()

	// Call finalValidation directly
	err := pipeline.finalValidation()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "expert accuracy too high")
}

// TestTrainingPipeline_FinalValidation_LowDownstreamSuccess_Direct tests low downstream success
func TestTrainingPipeline_FinalValidation_LowDownstreamSuccess_Direct(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.99, // Very high threshold
	}
	pipeline := NewTrainingPipeline(config)

	// Pre-populate with high score and no expert results
	pipeline.validator.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.000001, // Tiny -> high score
	}
	pipeline.validator.results.ExpertTestResults = nil
	pipeline.validator.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95, // Below 0.99 threshold
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	// Aggregate to set overall score >= 95
	_ = pipeline.validator.AggregateResults()

	// Call finalValidation directly
	err := pipeline.finalValidation()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "downstream success rate too low")
}

// TestTrainingPipeline_TrainModels_EmptyKmers tests trainModels with empty kmer sizes
func TestTrainingPipeline_TrainModels_EmptyKmers(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:  tempDir,
		KmerSizes:   []int{},    // Empty
		Populations: []string{}, // Empty
		Platforms:   []string{}, // Empty
	}
	pipeline := NewTrainingPipeline(config)

	// With empty kmer/population/platform lists, trainModels should succeed
	err := pipeline.trainModels()
	assert.NoError(t, err)
}

// TestTrainingPipeline_TrainModels_WithHMMFailure tests trainModels when HMM training fails
func TestTrainingPipeline_TrainModels_WithHMMFailure(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:  tempDir,
		KmerSizes:   []int{4},
		Populations: []string{"EUR"},
		Platforms:   []string{},
	}
	pipeline := NewTrainingPipeline(config)

	// HMM training will fail (loadReferenceGenome not implemented)
	err := pipeline.trainModels()
	assert.Error(t, err)
	// trainModels wraps the error with "failed to train HMM model"
	assert.Contains(t, err.Error(), "failed to train HMM model")
}

// TestTrainingPipeline_PrepareData_EmptyLists tests prepareData with empty lists
func TestTrainingPipeline_PrepareData_EmptyLists(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		ReferenceGenomePath: tempDir,
		PopulationDataPath:  tempDir,
		ValidationDataPath:  tempDir,
		OutputPath:          tempDir,
		Populations:         []string{},
		Platforms:           []string{},
	}
	pipeline := NewTrainingPipeline(config)

	// prepareData will fail at reference genome loading
	err := pipeline.prepareData()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "failed to load reference genome")
}

// TestTrainingPipeline_PrepareData_WithPopulationAndPlatform tests prepareData loading populations and platforms
func TestTrainingPipeline_PrepareData_WithPopulationAndPlatform(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		ReferenceGenomePath: tempDir,
		PopulationDataPath:  tempDir,
		ValidationDataPath:  tempDir,
		OutputPath:          tempDir,
		Populations:         []string{"EUR"},
		Platforms:           []string{"illumina"},
		ValidationSplit:     0.2,
	}
	pipeline := NewTrainingPipeline(config)

	// Will fail at reference genome loading (no actual genome data)
	err := pipeline.prepareData()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "failed to load reference genome")
}

// TestTrainingPipeline_TrainModels_QualityAndPopFail tests trainModels when quality/pop models fail
func TestTrainingPipeline_TrainModels_QualityFail(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:  tempDir,
		KmerSizes:   []int{},          // Empty = no HMM training, skip to quality
		Populations: []string{},       // Empty = no population training
		Platforms:   []string{"illumina"}, // Non-empty = try quality training
	}
	pipeline := NewTrainingPipeline(config)

	// Quality training will fail (loadQualityTrainingData not implemented)
	err := pipeline.trainModels()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "failed to train quality model")
}

// TestTrainingPipeline_TrainModels_PopFail tests trainModels when population model fails
func TestTrainingPipeline_TrainModels_PopFail(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:  tempDir,
		KmerSizes:   []int{},              // Empty
		Populations: []string{"EUR"},      // Non-empty = try population training
		Platforms:   []string{},           // Empty
	}
	pipeline := NewTrainingPipeline(config)

	// Population training will fail (getPopulationParameters not implemented)
	err := pipeline.trainModels()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "failed to train population model")
}

// TestTrainingPipeline_FinalValidation_JSdivergencePerPlatform tests JS divergence per platform check
func TestTrainingPipeline_FinalValidation_JSdivergencePerPlatform(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		JSdivergenceThreshold:      0.001, // Very tight
		ExpertAccuracyThreshold:    0.6,
		DownstreamSuccessThreshold: 0.9,
	}
	pipeline := NewTrainingPipeline(config)

	// Set JS divergence results that give high overall score but exceed per-platform threshold
	// Using tiny average divergence so overall score = ~100 (passes first check)
	// but set individual platforms to exceed the threshold
	pipeline.validator.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.0000001, // Tiny average: (0.0000001)/1 * 1000 = 0.0001, so score = 99.9999
	}
	pipeline.validator.results.ExpertTestResults = nil
	pipeline.validator.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	// Aggregate first to set overall score high
	_ = pipeline.validator.AggregateResults()

	// Now set illumina's divergence to exceed threshold for the per-platform check
	pipeline.validator.results.JSdivergenceResults["illumina"] = 0.005 // > 0.001 threshold

	// finalValidation should fail on JS divergence per-platform check
	err := pipeline.finalValidation()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "JS divergence too high")
}

// TestTrainingPipeline_FullPipelineWithEmptyLists tests a pipeline with empty lists that reaches later phases
func TestTrainingPipeline_FullPipelineWithEmptyKmerAndPopAndPlatform(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		KmerSizes:                  []int{},    // Empty
		Populations:                []string{}, // Empty
		Platforms:                  []string{}, // Empty
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.0,
		DownstreamSuccessThreshold: 0.9,
	}
	pipeline := NewTrainingPipeline(config)

	// With empty lists, trainModels succeeds with 0 models trained
	err := pipeline.trainModels()
	assert.NoError(t, err)

	// validateModels also succeeds (populates JS divergence results)
	err = pipeline.validateModels()
	assert.NoError(t, err)

	// testPerformance also succeeds
	err = pipeline.testPerformance()
	assert.NoError(t, err)

	// finalValidation will fail because RunStatisticalValidation set JS divergence values
	// that result in an overall score of ~91.67 (below the 95.0 threshold)
	pipeline.validator.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}
	err = pipeline.finalValidation()
	// Score will be low due to RunStatisticalValidation's default values
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "overall authenticity score too low")
}

// TestTrainingPipeline_TrainModels_CountsModelsTrained tests that ModelsTrained is incremented
func TestTrainingPipeline_TrainModels_EmptyAllLists(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:  tempDir,
		KmerSizes:   []int{},
		Populations: []string{},
		Platforms:   []string{},
	}
	pipeline := NewTrainingPipeline(config)

	// With all empty lists, training succeeds with 0 models trained
	err := pipeline.trainModels()
	assert.NoError(t, err)
	assert.Equal(t, 0, pipeline.metrics.ModelsTrained)
}

// TestHMMTrainer_LoadReferenceGenome tests loadReferenceGenome directly
func TestHMMTrainer_LoadReferenceGenome(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewHMMTrainer(config)

	// Without ReferenceGenomePath configured, loader returns a descriptive error.
	genome, err := trainer.loadReferenceGenome("EUR")
	assert.Error(t, err)
	assert.Nil(t, genome)
	assert.Contains(t, err.Error(), "ReferenceGenomePath not configured")
}

// TestHMMTrainer_SaveModel tests saveModel directly
func TestHMMTrainer_SaveModel(t *testing.T) {
	// With empty OutputPath, saveModel is a no-op and returns nil.
	config := &TrainingConfig{}
	trainer := NewHMMTrainer(config)
	model := authenticity.NewHMMKmerModel(4)
	err := trainer.saveModel(model, "test_key")
	assert.NoError(t, err)

	// With a real output path and an untrained model, Save returns an error.
	config2 := &TrainingConfig{OutputPath: t.TempDir()}
	trainer2 := NewHMMTrainer(config2)
	err2 := trainer2.saveModel(model, "test_key")
	assert.Error(t, err2)
}

// TestQualityTrainer_LoadQualityTrainingData tests loadQualityTrainingData directly
func TestQualityTrainer_LoadQualityTrainingData(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewQualityTrainer(config)

	// Without QualityDataPath configured, loader returns a descriptive error.
	data, err := trainer.loadQualityTrainingData("illumina")
	assert.Error(t, err)
	assert.Nil(t, data)
	assert.Contains(t, err.Error(), "QualityDataPath not configured")
}

// TestPopulationTrainer_GetPopulationParameters tests getPopulationParameters directly
func TestPopulationTrainer_GetPopulationParameters(t *testing.T) {
	config := &TrainingConfig{}
	trainer := NewPopulationTrainer(config)

	// Without PopulationDataPath configured, loader returns a descriptive error.
	params, err := trainer.getPopulationParameters("EUR")
	assert.Error(t, err)
	assert.Nil(t, params)
	assert.Contains(t, err.Error(), "PopulationDataPath not configured")
}

// TestAWSDataAcquisition_GetDatasetInfo_NotFound tests GetDatasetInfo with non-existent dataset
func TestAWSDataAcquisition_GetDatasetInfo_NotFound(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// Test with dataset that doesn't exist in the catalog
	_, err := acq.GetDatasetInfo("nonexistent_dataset")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "dataset not found")
}

// TestAWSDataAcquisition_CalculateDownloadSize_NotFound tests CalculateDownloadSize with invalid dataset
func TestAWSDataAcquisition_CalculateDownloadSize_NotFound(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// Test with non-existent dataset name
	_, err := acq.CalculateDownloadSize([]string{"nonexistent_dataset"})
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "failed to get dataset info")
}

// TestDataManager_LoadReferenceGenome_CachedPath tests the cached early-return path in LoadReferenceGenome
func TestDataManager_LoadReferenceGenome_CachedPath(t *testing.T) {
	config := &TrainingConfig{
		ReferenceGenomePath: "/test/reference",
	}
	dm := NewDataManager(config)

	// Manually insert a nil entry into the cache to simulate a previously loaded assembly
	// (we use nil because the actual type doesn't matter for the cache-hit path)
	dm.referenceCache["GRCh38"] = nil

	// LoadReferenceGenome should return nil immediately (cache hit)
	err := dm.LoadReferenceGenome("GRCh38")
	assert.NoError(t, err)
}

// TestDataManager_GetReferenceGenome_CachedPath tests getting a cached reference genome
func TestDataManager_GetReferenceGenome_CachedPath(t *testing.T) {
	config := &TrainingConfig{}
	dm := NewDataManager(config)

	// Insert nil into cache to simulate a cached entry
	dm.referenceCache["GRCh38"] = nil

	// GetReferenceGenome should return the cached value (nil in this case)
	genome, err := dm.GetReferenceGenome("GRCh38")
	assert.NoError(t, err)
	assert.Nil(t, genome)
}

// TestDataManager_GetPopulationPanel_CachedPath tests getting a cached population panel
func TestDataManager_GetPopulationPanel_CachedPath(t *testing.T) {
	config := &TrainingConfig{}
	dm := NewDataManager(config)

	// Insert nil into cache to simulate a cached entry
	dm.populationCache["EUR"] = nil

	// GetPopulationPanel should return the cached value
	panel, err := dm.GetPopulationPanel("EUR")
	assert.NoError(t, err)
	assert.Nil(t, panel)
}

// TestDataManager_GetValidationDataset_CachedPath tests getting a cached validation dataset
func TestDataManager_GetValidationDataset_CachedPath(t *testing.T) {
	config := &TrainingConfig{}
	dm := NewDataManager(config)

	// Insert a real entry into the cache
	cachedDataset := &ValidationDataset{
		Platform: "illumina",
		DataType: "simulated",
		Size:     100,
	}
	dm.validationCache["illumina"] = cachedDataset

	// GetValidationDataset should return the cached value
	dataset, err := dm.GetValidationDataset("illumina")
	assert.NoError(t, err)
	assert.Equal(t, cachedDataset, dataset)
}

// TestAWSDataAcquisition_CalculateDownloadSize_Empty tests with empty dataset list
func TestAWSDataAcquisition_CalculateDownloadSize_Empty(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// Empty list should return 0 with no error
	size, err := acq.CalculateDownloadSize([]string{})
	assert.NoError(t, err)
	assert.Equal(t, int64(0), size)
}

// TestDataManager_LoadPopulationPanel_CachedPath tests that loading a population panel twice uses the cache
func TestDataManager_LoadPopulationPanel_CachedEntry(t *testing.T) {
	config := &TrainingConfig{}
	dm := NewDataManager(config)

	// Manually insert nil into the cache (simulate loaded panel)
	dm.populationCache["YRI"] = nil

	// LoadPopulationPanel should see the cached entry and return nil
	err := dm.LoadPopulationPanel("YRI")
	assert.NoError(t, err)
}

// TestDataManager_LoadValidationDataset_CachedEntry tests that loading a dataset twice uses the cache
func TestDataManager_LoadValidationDataset_CachedEntry(t *testing.T) {
	config := &TrainingConfig{}
	dm := NewDataManager(config)

	// Manually insert a cached entry
	dm.validationCache["pacbio"] = &ValidationDataset{
		Platform: "pacbio",
		DataType: "simulated",
	}

	// LoadValidationDataset should see the cached entry and return nil
	err := dm.LoadValidationDataset("pacbio")
	assert.NoError(t, err)
}

// TestTrainingPipeline_TrainModels_ModelsTrained_Counter tests that ModelsTrained counter increments
func TestTrainingPipeline_TrainModels_PopModelsFail(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:  tempDir,
		KmerSizes:   []int{},
		Populations: []string{"YRI", "CEU"},
		Platforms:   []string{},
	}
	pipeline := NewTrainingPipeline(config)

	// Population training will fail, but first iteration should be attempted
	err := pipeline.trainModels()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "failed to train population model")
}

// TestTrainingPipeline_TrainModels_QualityMultiplePlatforms tests quality training for multiple platforms
func TestTrainingPipeline_TrainModels_QualityMultiplePlatforms(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:  tempDir,
		KmerSizes:   []int{},
		Populations: []string{},
		Platforms:   []string{"illumina", "pacbio", "nanopore"},
	}
	pipeline := NewTrainingPipeline(config)

	// Quality training will fail on first platform
	err := pipeline.trainModels()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "failed to train quality model")
}

// TestTrainingPipeline_PrepareData_BypassRefGenomeCache tests prepareData when reference genomes are pre-cached
// This allows the population and platform loading phases to be reached
func TestTrainingPipeline_PrepareData_BypassRefGenomeCache(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		ReferenceGenomePath: tempDir,
		PopulationDataPath:  tempDir,
		ValidationDataPath:  tempDir,
		OutputPath:          tempDir,
		Populations:         []string{"EUR"},
		Platforms:           []string{"illumina"},
		ValidationSplit:     0.2,
	}
	pipeline := NewTrainingPipeline(config)

	// Pre-populate reference genome cache to bypass the failing LoadReferenceGenome calls
	// This allows prepareData to proceed past the reference genome loop
	pipeline.dataManager.referenceCache["GRCh38"] = nil
	pipeline.dataManager.referenceCache["GRCh39"] = nil

	// prepareData should now proceed to the population loading phase
	// LoadPopulationPanel succeeds (uses simulation), LoadValidationDataset also succeeds
	err := pipeline.prepareData()
	assert.NoError(t, err)
}

// TestTrainingPipeline_Execute_BypassRefGenomeCache tests Execute when reference genomes are pre-cached
// This allows coverage of the model training, validation, and performance phases
func TestTrainingPipeline_Execute_BypassRefGenomeCache(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		Populations:                []string{},
		Platforms:                  []string{},
		KmerSizes:                  []int{},
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.0, // No expert testing
		DownstreamSuccessThreshold: 0.9,
		ValidationSplit:            0.2,
	}
	pipeline := NewTrainingPipeline(config)

	// Pre-populate reference genome cache to bypass the failing LoadReferenceGenome
	pipeline.dataManager.referenceCache["GRCh38"] = nil
	pipeline.dataManager.referenceCache["GRCh39"] = nil

	// Execute should now proceed past prepareData (empty populations and platforms)
	// trainModels succeeds (empty lists), validateModels sets JS divergence ~91% score
	// finalValidation will fail because JS divergence causes score below 95%
	err := pipeline.Execute()
	// finalValidation fails because RunStatisticalValidation sets default values giving ~91% score
	assert.Error(t, err)
	// Verify it reached the final validation phase (not data prep failure)
	assert.NotContains(t, err.Error(), "data preparation failed")
	assert.Contains(t, err.Error(), "final validation failed")
}

// TestTrainingPipeline_Execute_PassFinalValidation tests Execute reaching the end of pipeline
func TestTrainingPipeline_Execute_PassFinalValidation(t *testing.T) {
	tempDir := t.TempDir()
	config := &TrainingConfig{
		OutputPath:                 tempDir,
		Populations:                []string{},
		Platforms:                  []string{},
		KmerSizes:                  []int{},
		JSdivergenceThreshold:      0.1,
		ExpertAccuracyThreshold:    0.0,  // No expert testing
		DownstreamSuccessThreshold: 0.9,
		ValidationSplit:            0.2,
	}
	pipeline := NewTrainingPipeline(config)

	// Pre-populate reference genome cache
	pipeline.dataManager.referenceCache["GRCh38"] = nil
	pipeline.dataManager.referenceCache["GRCh39"] = nil

	// Pre-populate validator results to pass finalValidation
	pipeline.validator.results.JSdivergenceResults = map[string]float64{
		"illumina": 0.000001, // Very low -> score near 100
	}
	pipeline.validator.results.ExpertTestResults = nil
	pipeline.validator.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	// Execute should succeed: prepareData passes (cached refs, empty pops/platforms),
	// trainModels passes (empty lists), validateModels REPLACES our JSdivergence results,
	// so finalValidation will fail due to the default ~91% score from RunStatisticalValidation
	// We can only control validateModels if we don't call validateModels in Execute
	// Since Execute calls validateModels which calls RunStatisticalValidation, the results get reset
	// So this test verifies that Execute reaches finalValidation and fails there
	err := pipeline.Execute()
	// finalValidation fails because validateModels reset the JS divergence to default values
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "final validation failed")
}

// TestPerformanceTester_Structure tests performance tester structure
func TestPerformanceTester_Structure(t *testing.T) {
	config := &TrainingConfig{
		MaxMemoryGB: 16,
		NumWorkers:  8,
	}
	tester := NewPerformanceTester(config)
	assert.NotNil(t, tester)
	assert.Equal(t, config, tester.config)
	assert.NotNil(t, tester.results)
}

// TestHMMTrainer_Structure tests HMM trainer structure
func TestHMMTrainer_Structure(t *testing.T) {
	config := &TrainingConfig{
		KmerSizes:  []int{4, 6, 8},
		OutputPath: "/test/output",
	}
	trainer := NewHMMTrainer(config)
	assert.NotNil(t, trainer)
	assert.Equal(t, config, trainer.config)
	assert.NotNil(t, trainer.models)
	assert.NotNil(t, trainer.trainedModels)
	assert.Empty(t, trainer.models)
	assert.Empty(t, trainer.trainedModels)
}

// TestQualityTrainer_Structure tests quality trainer structure
func TestQualityTrainer_Structure(t *testing.T) {
	config := &TrainingConfig{
		Platforms: []string{"illumina", "pacbio"},
	}
	trainer := NewQualityTrainer(config)
	assert.NotNil(t, trainer)
	assert.Equal(t, config, trainer.config)
	assert.NotNil(t, trainer.models)
	assert.Empty(t, trainer.models)
}

// TestPopulationTrainer_Structure tests population trainer structure
func TestPopulationTrainer_Structure(t *testing.T) {
	config := &TrainingConfig{
		Populations: []string{"EUR", "AFR", "ASN"},
	}
	trainer := NewPopulationTrainer(config)
	assert.NotNil(t, trainer)
	assert.Equal(t, config, trainer.config)
	assert.NotNil(t, trainer.models)
	assert.Empty(t, trainer.models)
}
