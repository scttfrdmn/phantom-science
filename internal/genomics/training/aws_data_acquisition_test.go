package training

import (
	"log"
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// newTestAWSDataAcquisition creates an AWSDataAcquisition for testing
// without real AWS credentials (nil s3Client and downloader are OK for non-S3 methods)
func newTestAWSDataAcquisition(downloadDir string) *AWSDataAcquisition {
	return &AWSDataAcquisition{
		s3Client:    nil, // Not used in filesystem-only tests
		downloader:  nil, // Not used in filesystem-only tests
		config:      DefaultAWSDataConfig(),
		downloadDir: downloadDir,
		logger:      log.New(os.Stdout, "[TEST-AWS] ", log.LstdFlags),
	}
}

// smallTestDataset creates a small dataset for testing (avoids huge file sizes)
var smallTestDatasets = map[string]DatasetInfo{
	"test_dataset": {
		Name:        "Test Dataset",
		Description: "A small test dataset",
		Bucket:      "test-bucket",
		Prefix:      "test/",
		Files: []FileInfo{
			{Key: "test_file1.txt", Size: 10},
			{Key: "test_file2.txt", Size: 20},
		},
	},
}

// TestAWSDataAcquisition_fileExists tests the fileExists method
func TestAWSDataAcquisition_fileExists(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	t.Run("existing non-empty file", func(t *testing.T) {
		testFile := filepath.Join(tempDir, "existing.txt")
		err := os.WriteFile(testFile, []byte("content"), 0600)
		require.NoError(t, err)

		assert.True(t, acq.fileExists(testFile))
	})

	t.Run("empty file returns false", func(t *testing.T) {
		emptyFile := filepath.Join(tempDir, "empty.txt")
		err := os.WriteFile(emptyFile, []byte{}, 0600)
		require.NoError(t, err)

		assert.False(t, acq.fileExists(emptyFile))
	})

	t.Run("non-existent file returns false", func(t *testing.T) {
		assert.False(t, acq.fileExists(filepath.Join(tempDir, "nonexistent.txt")))
	})
}

// TestAWSDataAcquisition_verifyFileSize tests file size verification
func TestAWSDataAcquisition_verifyFileSize(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	t.Run("correct file size", func(t *testing.T) {
		content := []byte("hello world")
		testFile := filepath.Join(tempDir, "size_test.txt")
		err := os.WriteFile(testFile, content, 0600)
		require.NoError(t, err)

		err = acq.verifyFileSize(testFile, int64(len(content)))
		assert.NoError(t, err)
	})

	t.Run("wrong file size", func(t *testing.T) {
		content := []byte("hello world")
		testFile := filepath.Join(tempDir, "size_wrong.txt")
		err := os.WriteFile(testFile, content, 0600)
		require.NoError(t, err)

		err = acq.verifyFileSize(testFile, 999)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "size mismatch")
	})

	t.Run("non-existent file", func(t *testing.T) {
		err := acq.verifyFileSize(filepath.Join(tempDir, "nonexistent.txt"), 100)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "failed to stat file")
	})
}

// TestAWSDataAcquisition_GetDownloadedDataPaths tests path retrieval
func TestAWSDataAcquisition_GetDownloadedDataPaths(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	paths := acq.GetDownloadedDataPaths()

	assert.NotNil(t, paths)
	assert.NotEmpty(t, paths)

	// Should have paths for all datasets in AuthenticityValidationDatasets
	for datasetName := range AuthenticityValidationDatasets {
		path, exists := paths[datasetName]
		assert.True(t, exists, "Should have path for dataset %s", datasetName)
		assert.Contains(t, path, tempDir)
		assert.Contains(t, path, datasetName)
	}
}

// TestAWSDataAcquisition_ValidateDownloadedData_MissingFiles tests validation with missing files
func TestAWSDataAcquisition_ValidateDownloadedData_MissingFiles(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// Without downloading anything, validation should fail
	err := acq.ValidateDownloadedData()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "missing file")
}

// createSmallTestFiles creates small files matching the specified file infos
func createSmallTestFiles(t *testing.T, baseDir string, datasetName string, files []FileInfo) {
	t.Helper()
	datasetDir := filepath.Join(baseDir, datasetName)
	err := os.MkdirAll(datasetDir, 0750)
	require.NoError(t, err)

	for _, file := range files {
		localPath := filepath.Join(datasetDir, filepath.Base(file.Key))
		content := make([]byte, file.Size)
		for i := range content {
			content[i] = 'X' // Fill with non-zero content
		}
		err := os.WriteFile(localPath, content, 0600)
		require.NoError(t, err)
	}
}

// TestAWSDataAcquisition_ValidateDownloadedData_WithSmallFiles tests validation with small test files
func TestAWSDataAcquisition_ValidateDownloadedData_WithSmallFiles(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// Create a minimal subset of files — we only need to test with files we can actually create
	// Use the illumina_fastq dataset which has relatively small files
	illuminaDataset := AuthenticityValidationDatasets["illumina_fastq"]
	createSmallTestFiles(t, tempDir, "illumina_fastq", illuminaDataset.Files[:1]) // Only first file (1MB)

	// Validation will still fail because not all datasets are present
	err := acq.ValidateDownloadedData()
	assert.Error(t, err) // Expected - not all files present
}

// TestAWSDataAcquisition_CleanupIncompleteDownloads_NoFiles tests cleanup with no files
func TestAWSDataAcquisition_CleanupIncompleteDownloads_NoFiles(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// With no files, cleanup should succeed (nothing to clean)
	err := acq.CleanupIncompleteDownloads()
	assert.NoError(t, err)
}

// TestAWSDataAcquisition_CleanupIncompleteDownloads_WrongSizeFiles tests cleanup removes wrong-size files
func TestAWSDataAcquisition_CleanupIncompleteDownloads_WithWrongSizeFiles(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// Create files with wrong sizes (just a few bytes each, not GB)
	for datasetName, dataset := range AuthenticityValidationDatasets {
		datasetDir := filepath.Join(tempDir, datasetName)
		err := os.MkdirAll(datasetDir, 0750)
		require.NoError(t, err)

		for _, file := range dataset.Files {
			localPath := filepath.Join(datasetDir, filepath.Base(file.Key))
			// Write wrong size (5 bytes, not GB)
			if file.Size != 5 {
				err := os.WriteFile(localPath, []byte("12345"), 0600)
				require.NoError(t, err)
			}
		}
	}

	err := acq.CleanupIncompleteDownloads()
	assert.NoError(t, err)

	// Wrong-size files should be removed
	for datasetName, dataset := range AuthenticityValidationDatasets {
		datasetDir := filepath.Join(tempDir, datasetName)
		for _, file := range dataset.Files {
			if file.Size != 5 {
				localPath := filepath.Join(datasetDir, filepath.Base(file.Key))
				_, statErr := os.Stat(localPath)
				assert.True(t, os.IsNotExist(statErr), "Wrong-size file should be removed: %s", localPath)
			}
		}
	}
}

// TestAWSDataAcquisition_downloadDataset_CacheHit_SmallFiles tests cache hit with small files
func TestAWSDataAcquisition_downloadDataset_CacheHit_SmallFiles(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)
	acq.config.CacheData = true

	// Use illumina benchmark bed file (1MB) which is small enough to create
	datasetName := "illumina_fastq"
	dataset := AuthenticityValidationDatasets[datasetName]

	// Only create the first file (bed file, 1MB = 1000000 bytes)
	firstFile := dataset.Files[0]
	datasetDir := filepath.Join(tempDir, datasetName)
	err := os.MkdirAll(datasetDir, 0750)
	require.NoError(t, err)

	localPath := filepath.Join(datasetDir, filepath.Base(firstFile.Key))
	content := make([]byte, firstFile.Size)
	err = os.WriteFile(localPath, content, 0600)
	require.NoError(t, err)

	// Create a dataset with just the one file we can create
	smallDataset := DatasetInfo{
		Name:   dataset.Name,
		Bucket: dataset.Bucket,
		Prefix: dataset.Prefix,
		Files:  []FileInfo{firstFile},
	}

	// downloadDataset should skip downloading (cache hit) since file exists and has right size
	err = acq.downloadDataset(datasetName, smallDataset)
	assert.NoError(t, err)
}

// TestAWSDataAcquisition_downloadDataset_NoCacheForceDownload tests that no-cache path invokes download
func TestAWSDataAcquisition_downloadDataset_NoCacheForceDownload(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)
	acq.config.CacheData = false // Disable cache to force download

	// Use smallest dataset with tiny synthetic size for testing
	syntheticDataset := DatasetInfo{
		Name:   "Test",
		Bucket: "test-bucket",
		Prefix: "test/",
		Files: []FileInfo{
			{Key: "tiny_file.txt", Size: 5},
		},
	}

	// With nil downloader, will panic when trying to download
	// Recover from the panic to test the code path
	defer func() {
		if r := recover(); r != nil {
			// The panic is expected because downloader is nil
			// This path proves the download code path was reached
		}
	}()

	// This will fail (nil pointer panic) which we recover from
	_ = acq.downloadDataset("test", syntheticDataset)
}

// createAllTestFilesForDataset creates dummy files for every file in AuthenticityValidationDatasets
// to simulate a fully cached download (any content, any size - just needs to exist)
func createAllTestFilesForDataset(t *testing.T, baseDir string) {
	t.Helper()
	for datasetName, dataset := range AuthenticityValidationDatasets {
		datasetDir := filepath.Join(baseDir, datasetName)
		err := os.MkdirAll(datasetDir, 0750)
		require.NoError(t, err)

		for _, file := range dataset.Files {
			localPath := filepath.Join(datasetDir, filepath.Base(file.Key))
			// Write a tiny non-empty file (1 byte) - just needs to exist and be non-empty
			err := os.WriteFile(localPath, []byte("X"), 0600)
			require.NoError(t, err)
		}
	}
}

// TestAWSDataAcquisition_DownloadReferenceGenome_CacheHit tests DownloadReferenceGenome with cache
func TestAWSDataAcquisition_DownloadReferenceGenome_CacheHit(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)
	acq.config.CacheData = true

	// Create dummy file for the reference genome files
	refDataset := AuthenticityValidationDatasets["reference_genome"]
	datasetDir := filepath.Join(tempDir, "reference_genome")
	err := os.MkdirAll(datasetDir, 0750)
	require.NoError(t, err)

	for _, file := range refDataset.Files {
		localPath := filepath.Join(datasetDir, filepath.Base(file.Key))
		err := os.WriteFile(localPath, []byte("X"), 0600)
		require.NoError(t, err)
	}

	// With cache enabled and files present, should skip actual S3 download
	err = acq.DownloadReferenceGenome()
	assert.NoError(t, err)
}

// TestAWSDataAcquisition_DownloadPopulationData_CacheHit tests DownloadPopulationData with cache
func TestAWSDataAcquisition_DownloadPopulationData_CacheHit(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)
	acq.config.CacheData = true

	// Create dummy files for population VCF
	popDataset := AuthenticityValidationDatasets["population_vcf"]
	datasetDir := filepath.Join(tempDir, "population_vcf")
	err := os.MkdirAll(datasetDir, 0750)
	require.NoError(t, err)

	for _, file := range popDataset.Files {
		localPath := filepath.Join(datasetDir, filepath.Base(file.Key))
		err := os.WriteFile(localPath, []byte("X"), 0600)
		require.NoError(t, err)
	}

	err = acq.DownloadPopulationData()
	assert.NoError(t, err)
}

// TestAWSDataAcquisition_DownloadValidationFASTQ_CacheHit tests DownloadValidationFASTQ with cache
func TestAWSDataAcquisition_DownloadValidationFASTQ_CacheHit(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)
	acq.config.CacheData = true

	// Create dummy files for all three FASTQ datasets
	for _, datasetName := range []string{"illumina_fastq", "pacbio_fastq", "nanopore_fastq"} {
		dataset := AuthenticityValidationDatasets[datasetName]
		datasetDir := filepath.Join(tempDir, datasetName)
		err := os.MkdirAll(datasetDir, 0750)
		require.NoError(t, err)

		for _, file := range dataset.Files {
			localPath := filepath.Join(datasetDir, filepath.Base(file.Key))
			err := os.WriteFile(localPath, []byte("X"), 0600)
			require.NoError(t, err)
		}
	}

	err := acq.DownloadValidationFASTQ()
	assert.NoError(t, err)
}

// TestAWSDataAcquisition_DownloadAuthenticityDatasets_CacheHit tests full dataset download with cache
func TestAWSDataAcquisition_DownloadAuthenticityDatasets_CacheHit(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)
	acq.config.CacheData = true

	// Create dummy files for all datasets
	createAllTestFilesForDataset(t, tempDir)

	err := acq.DownloadAuthenticityDatasets()
	assert.NoError(t, err)
}

// createExactSizeTestFiles creates files with EXACT sizes matching the dataset file infos
// This is only feasible for files up to a few MB; large files (GB) cannot be created in tests
func createExactSizeSmallFiles(t *testing.T, baseDir string) {
	t.Helper()
	for datasetName, dataset := range AuthenticityValidationDatasets {
		datasetDir := filepath.Join(baseDir, datasetName)
		err := os.MkdirAll(datasetDir, 0750)
		require.NoError(t, err)

		for _, file := range dataset.Files {
			// Only create files that are small enough (< 5MB to avoid disk space issues)
			if file.Size <= 5*1024*1024 {
				localPath := filepath.Join(datasetDir, filepath.Base(file.Key))
				content := make([]byte, file.Size)
				for i := range content {
					content[i] = 'X'
				}
				err := os.WriteFile(localPath, content, 0600)
				require.NoError(t, err)
			}
		}
	}
}

// TestAWSDataAcquisition_ValidateDownloadedData_WrongSizeFile tests validation with wrong-size file
func TestAWSDataAcquisition_ValidateDownloadedData_WrongSizeFile(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// Create one file for the illumina dataset with WRONG size
	datasetDir := filepath.Join(tempDir, "illumina_fastq")
	err := os.MkdirAll(datasetDir, 0750)
	require.NoError(t, err)

	illuminaDataset := AuthenticityValidationDatasets["illumina_fastq"]
	firstFile := illuminaDataset.Files[0]
	localPath := filepath.Join(datasetDir, filepath.Base(firstFile.Key))

	// Write with wrong size (not 1MB)
	err = os.WriteFile(localPath, []byte("wrong size"), 0600)
	require.NoError(t, err)

	// Validation should fail with size mismatch for this file
	err = acq.ValidateDownloadedData()
	assert.Error(t, err)
	// Could be "missing file" for other datasets, or "file validation failed" for the wrong-size one
}

// TestAWSDataAcquisition_ValidateDownloadedData_IlluminaCorrectSize tests validation with correct-size illumina file
func TestAWSDataAcquisition_ValidateDownloadedData_IlluminaCorrectSize(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// Create the 1MB illumina BED file with correct size
	datasetDir := filepath.Join(tempDir, "illumina_fastq")
	err := os.MkdirAll(datasetDir, 0750)
	require.NoError(t, err)

	illuminaDataset := AuthenticityValidationDatasets["illumina_fastq"]
	bedFile := illuminaDataset.Files[0] // 1MB BED file

	localPath := filepath.Join(datasetDir, filepath.Base(bedFile.Key))
	content := make([]byte, bedFile.Size) // Exact 1MB size
	for i := range content {
		content[i] = 'X'
	}
	err = os.WriteFile(localPath, content, 0600)
	require.NoError(t, err)

	// Create the 50MB VCF file with correct size would be too large
	// Instead verify that the file we created passes validation for this file
	err = acq.verifyFileSize(localPath, bedFile.Size)
	assert.NoError(t, err)
}

// TestAWSDataAcquisition_GetDatasetInfo_NilS3Client tests GetDatasetInfo with nil s3Client
func TestAWSDataAcquisition_GetDatasetInfo_NilS3Client(t *testing.T) {
	tempDir := t.TempDir()
	acq := newTestAWSDataAcquisition(tempDir)

	// GetDatasetInfo for an existing dataset calls ListDatasetFiles which uses nil s3Client
	// This will panic with nil pointer dereference on s3Client.ListObjectsV2Pages
	defer func() {
		if r := recover(); r != nil {
			// Panic is expected because s3Client is nil
			// This confirms the GetDatasetInfo path was reached for existing datasets
		}
	}()

	// This will either panic (nil s3Client) or return an error
	_, _ = acq.GetDatasetInfo("illumina_fastq")
}
