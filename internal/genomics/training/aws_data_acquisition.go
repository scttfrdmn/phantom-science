package training

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"sync"
	"time"

	"github.com/aws/aws-sdk-go/aws"
	"github.com/aws/aws-sdk-go/aws/session"
	"github.com/aws/aws-sdk-go/service/s3"
	"github.com/aws/aws-sdk-go/service/s3/s3manager"
)

// AWSDataAcquisition downloads files listed in AuthenticityValidationDatasets
// from public S3 buckets (AWS Open Data).
type AWSDataAcquisition struct {
	s3Client    *s3.S3
	downloader  *s3manager.Downloader
	config      *AWSDataConfig
	downloadDir string
	logger      *log.Logger
	mu          sync.RWMutex
}

// AWSDataConfig contains configuration for AWS data acquisition
type AWSDataConfig struct {
	Region          string
	DownloadDir     string
	MaxConcurrency  int
	MaxRetries      int  // Not used yet
	VerifyChecksums bool // Not used yet: only file sizes are checked
	EnableProgress  bool // Not used yet
	CacheData       bool // Skip files that already exist and are non-empty
	CleanupOnError  bool // Not used yet
}

// DatasetInfo represents information about a genomic dataset
type DatasetInfo struct {
	Name         string
	Description  string
	Bucket       string
	Prefix       string
	Files        []FileInfo
	Size         int64
	LastModified time.Time
	Checksum     string
}

// FileInfo represents information about a single file
type FileInfo struct {
	Key          string
	Size         int64
	LastModified time.Time
	ETag         string
	LocalPath    string
}

// DownloadProgress tracks download progress
type DownloadProgress struct {
	TotalFiles      int
	CompletedFiles  int
	TotalBytes      int64
	DownloadedBytes int64
	StartTime       time.Time
	CurrentFile     string
	Speed           float64 // MB/s
	ETA             time.Duration
}

// GenomicDatasets lists AWS Open Data buckets and prefixes. The entries are
// not checked against S3 by any test; see the notes on individual entries.
var GenomicDatasets = map[string]DatasetInfo{
	"grch38": {
		Name:        "GRCh38 Reference Genome",
		Description: "Human reference genome (the prefix points at hs37d5.fa, which is GRCh37 with decoys, not GRCh38)",
		Bucket:      "aws-roda-hcls-datalake",
		Prefix:      "1000_genomes_dragen_3.7.6/reference/hs37d5.fa",
	},
	"grch39": {
		Name:        "GRCh39 Reference Genome",
		Description: "Reference directory of the 1000 Genomes DRAGEN 3.7.6 release (GRCh39 has not been released)",
		Bucket:      "aws-roda-hcls-datalake",
		Prefix:      "1000_genomes_dragen_3.7.6/reference/",
	},
	"1000genomes": {
		Name:        "1000 Genomes Project",
		Description: "1000 Genomes Project Phase 3 data",
		Bucket:      "aws-roda-hcls-datalake",
		Prefix:      "1000_genomes_dragen_3.7.6/",
	},
	"giab": {
		Name:        "Genome in a Bottle",
		Description: "NIST Genome in a Bottle reference materials",
		Bucket:      "giab",
		Prefix:      "release/",
	},
	"broad_references": {
		Name:        "Broad Institute References",
		Description: "Broad Institute reference genomes and resources",
		Bucket:      "broad-references",
		Prefix:      "hg38/",
	},
	"encode": {
		Name:        "ENCODE Project",
		Description: "ENCODE Project functional genomics data",
		Bucket:      "encode-public",
		Prefix:      "2020/",
	},
	"gnomad": {
		Name:        "gnomAD",
		Description: "Genome Aggregation Database variants",
		Bucket:      "gnomad-public-us-east-1",
		Prefix:      "release/",
	},
	"tcga": {
		Name:        "The Cancer Genome Atlas",
		Description: "TCGA genomic data",
		Bucket:      "tcga-2-open",
		Prefix:      "",
	},
	"sra": {
		Name:        "Sequence Read Archive",
		Description: "NCBI Sequence Read Archive on AWS",
		Bucket:      "sra-pub-run-odp",
		Prefix:      "sra/",
	},
	"allen_brain": {
		Name:        "Allen Brain Observatory",
		Description: "Allen Institute Brain Observatory data",
		Bucket:      "allen-brain-observatory",
		Prefix:      "visual-coding-2p/",
	},
}

// AuthenticityValidationDatasets lists the files that
// DownloadAuthenticityDatasets fetches. File sizes are approximate, but
// downloadDataset and ValidateDownloadedData require an exact size match.
var AuthenticityValidationDatasets = map[string]DatasetInfo{
	"reference_genome": {
		Name:        "GRCh38 Reference Genome",
		Description: "Complete human reference genome for model training",
		Bucket:      "broad-references",
		Prefix:      "hg38/v0/",
		Files: []FileInfo{
			{Key: "Homo_sapiens_assembly38.fasta", Size: 3200000000},
			{Key: "Homo_sapiens_assembly38.fasta.fai", Size: 3000},
			{Key: "Homo_sapiens_assembly38.dict", Size: 7000}, // gitleaks:allow -- S3 object key, not a credential
		},
	},
	"population_vcf": {
		Name:        "1000 Genomes Population VCF",
		Description: "Population genetics data for authenticity training",
		Bucket:      "aws-roda-hcls-datalake",
		Prefix:      "1000_genomes_dragen_3.7.6/vcf/",
		Files: []FileInfo{
			{Key: "ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz", Size: 500000000},
			{Key: "ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz", Size: 400000000},
			{Key: "ALL.chr20.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz", Size: 600000000},
		},
	},
	"illumina_fastq": {
		Name:        "Illumina FASTQ Samples",
		Description: "GIAB HG001 v4.2.1 benchmark BED and VCF (no FASTQ files)",
		Bucket:      "giab",
		Prefix:      "release/NA12878_HG001/NISTv4.2.1/GRCh38/",
		Files: []FileInfo{
			{Key: "HG001_GRCh38_1_22_v4.2.1_benchmark.bed", Size: 1000000},
			{Key: "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz", Size: 50000000},
		},
	},
	"pacbio_fastq": {
		Name:        "PacBio HiFi FASTQ Samples",
		Description: "Real PacBio HiFi sequencing data for validation",
		Bucket:      "human-pangenomics",
		Prefix:      "working/HPRC_PLUS/HG002/raw_data/pacbio_hifi/",
		Files: []FileInfo{
			{Key: "m64011_190830_220126.Q20.fastq.gz", Size: 5000000000},
			{Key: "m64011_190901_095311.Q20.fastq.gz", Size: 5000000000},
		},
	},
	"nanopore_fastq": {
		Name:        "Nanopore FASTQ Samples",
		Description: "Real Oxford Nanopore sequencing data for validation",
		Bucket:      "human-pangenomics",
		Prefix:      "working/HPRC_PLUS/HG002/raw_data/oxford_nanopore/",
		Files: []FileInfo{
			{Key: "HG002_ONT_UL_guppy_4.2.2.fastq.gz", Size: 10000000000},
		},
	},
}

// NewAWSDataAcquisition creates a new AWS data acquisition client
func NewAWSDataAcquisition(config *AWSDataConfig) (*AWSDataAcquisition, error) {
	// Create AWS session
	sess, err := session.NewSession(&aws.Config{
		Region: aws.String(config.Region),
	})
	if err != nil {
		return nil, fmt.Errorf("failed to create AWS session: %w", err)
	}

	// Create S3 client
	s3Client := s3.New(sess)

	// Create downloader with configuration
	downloader := s3manager.NewDownloader(sess, func(d *s3manager.Downloader) {
		d.Concurrency = config.MaxConcurrency
		d.PartSize = 64 * 1024 * 1024 // 64MB parts
	})

	// Create download directory
	if err := os.MkdirAll(config.DownloadDir, 0750); err != nil {
		return nil, fmt.Errorf("failed to create download directory: %w", err)
	}

	return &AWSDataAcquisition{
		s3Client:    s3Client,
		downloader:  downloader,
		config:      config,
		downloadDir: config.DownloadDir,
		logger:      log.New(os.Stdout, "[AWS-DATA] ", log.LstdFlags),
	}, nil
}

// DownloadAuthenticityDatasets downloads all datasets needed for authenticity validation
func (a *AWSDataAcquisition) DownloadAuthenticityDatasets() error {
	a.logger.Println("Starting download of authenticity validation datasets")

	totalSize := int64(0)
	for _, dataset := range AuthenticityValidationDatasets {
		for _, file := range dataset.Files {
			totalSize += file.Size
		}
	}

	a.logger.Printf("Total data to download: %.2f GB", float64(totalSize)/(1024*1024*1024))

	// Download each dataset
	for datasetName, dataset := range AuthenticityValidationDatasets {
		a.logger.Printf("Downloading dataset: %s", dataset.Name)

		if err := a.downloadDataset(datasetName, dataset); err != nil {
			return fmt.Errorf("failed to download dataset %s: %w", datasetName, err)
		}
	}

	a.logger.Println("All authenticity validation datasets downloaded successfully")
	return nil
}

// downloadDataset downloads a specific dataset
func (a *AWSDataAcquisition) downloadDataset(datasetName string, dataset DatasetInfo) error {
	datasetDir := filepath.Join(a.downloadDir, datasetName)
	if err := os.MkdirAll(datasetDir, 0750); err != nil {
		return fmt.Errorf("failed to create dataset directory: %w", err)
	}

	// Download each file in the dataset
	for _, file := range dataset.Files {
		localPath := filepath.Join(datasetDir, filepath.Base(file.Key))

		// Check if file already exists and is valid
		if a.fileExists(localPath) && a.config.CacheData {
			a.logger.Printf("File already exists, skipping: %s", file.Key)
			continue
		}

		if err := a.downloadFile(dataset.Bucket, file.Key, localPath); err != nil {
			return fmt.Errorf("failed to download file %s: %w", file.Key, err)
		}

		// Verify file size
		if err := a.verifyFileSize(localPath, file.Size); err != nil {
			return fmt.Errorf("file size verification failed for %s: %w", file.Key, err)
		}
	}

	return nil
}

// downloadFile downloads a single file from S3
func (a *AWSDataAcquisition) downloadFile(bucket, key, localPath string) error {
	a.logger.Printf("Downloading: s3://%s/%s -> %s", bucket, key, localPath)

	// Create local file
	file, err := os.Create(localPath) // #nosec G304
	if err != nil {
		return fmt.Errorf("failed to create local file: %w", err)
	}
	defer file.Close()

	// Download from S3
	startTime := time.Now()
	bytesDownloaded, err := a.downloader.Download(file, &s3.GetObjectInput{
		Bucket: aws.String(bucket),
		Key:    aws.String(key),
	})

	if err != nil {
		// Clean up partial download
		os.Remove(localPath)
		return fmt.Errorf("S3 download failed: %w", err)
	}

	duration := time.Since(startTime)
	speed := float64(bytesDownloaded) / duration.Seconds() / (1024 * 1024) // MB/s

	a.logger.Printf("Downloaded %s (%.2f MB) in %.2fs at %.2f MB/s",
		key, float64(bytesDownloaded)/(1024*1024), duration.Seconds(), speed)

	return nil
}

// fileExists checks if a file exists and is not empty
func (a *AWSDataAcquisition) fileExists(path string) bool {
	info, err := os.Stat(path)
	return err == nil && info.Size() > 0
}

// verifyFileSize verifies that downloaded file has expected size
func (a *AWSDataAcquisition) verifyFileSize(path string, expectedSize int64) error {
	info, err := os.Stat(path)
	if err != nil {
		return fmt.Errorf("failed to stat file: %w", err)
	}

	if info.Size() != expectedSize {
		return fmt.Errorf("size mismatch: expected %d, got %d", expectedSize, info.Size())
	}

	return nil
}

// ListDatasetFiles lists all files in a dataset
func (a *AWSDataAcquisition) ListDatasetFiles(bucket, prefix string) ([]FileInfo, error) {
	input := &s3.ListObjectsV2Input{
		Bucket: aws.String(bucket),
		Prefix: aws.String(prefix),
	}

	var files []FileInfo

	err := a.s3Client.ListObjectsV2Pages(input, func(page *s3.ListObjectsV2Output, lastPage bool) bool {
		for _, obj := range page.Contents {
			files = append(files, FileInfo{
				Key:          *obj.Key,
				Size:         *obj.Size,
				LastModified: *obj.LastModified,
				ETag:         *obj.ETag,
			})
		}
		return !lastPage
	})

	if err != nil {
		return nil, fmt.Errorf("failed to list objects: %w", err)
	}

	return files, nil
}

// DownloadReferenceGenome downloads reference genome specifically for training
func (a *AWSDataAcquisition) DownloadReferenceGenome() error {
	a.logger.Println("Downloading reference genome for model training")

	// Download GRCh38 reference genome
	dataset := AuthenticityValidationDatasets["reference_genome"]
	return a.downloadDataset("reference_genome", dataset)
}

// DownloadPopulationData downloads population genetics data
func (a *AWSDataAcquisition) DownloadPopulationData() error {
	a.logger.Println("Downloading population genetics data")

	// Download 1000 Genomes population VCF files
	dataset := AuthenticityValidationDatasets["population_vcf"]
	return a.downloadDataset("population_vcf", dataset)
}

// DownloadValidationFASTQ downloads FASTQ files for validation
func (a *AWSDataAcquisition) DownloadValidationFASTQ() error {
	a.logger.Println("Downloading validation FASTQ datasets")

	// Download Illumina FASTQ
	if err := a.downloadDataset("illumina_fastq", AuthenticityValidationDatasets["illumina_fastq"]); err != nil {
		return fmt.Errorf("failed to download Illumina FASTQ: %w", err)
	}

	// Download PacBio FASTQ
	if err := a.downloadDataset("pacbio_fastq", AuthenticityValidationDatasets["pacbio_fastq"]); err != nil {
		return fmt.Errorf("failed to download PacBio FASTQ: %w", err)
	}

	// Download Nanopore FASTQ
	if err := a.downloadDataset("nanopore_fastq", AuthenticityValidationDatasets["nanopore_fastq"]); err != nil {
		return fmt.Errorf("failed to download Nanopore FASTQ: %w", err)
	}

	return nil
}

// GetDatasetInfo returns information about available datasets
func (a *AWSDataAcquisition) GetDatasetInfo(datasetName string) (*DatasetInfo, error) {
	if dataset, exists := AuthenticityValidationDatasets[datasetName]; exists {
		// Update file list from S3
		files, err := a.ListDatasetFiles(dataset.Bucket, dataset.Prefix)
		if err != nil {
			return nil, fmt.Errorf("failed to list dataset files: %w", err)
		}

		dataset.Files = files
		return &dataset, nil
	}

	return nil, fmt.Errorf("dataset not found: %s", datasetName)
}

// CalculateDownloadSize calculates total size of datasets to download
func (a *AWSDataAcquisition) CalculateDownloadSize(datasetNames []string) (int64, error) {
	totalSize := int64(0)

	for _, datasetName := range datasetNames {
		dataset, err := a.GetDatasetInfo(datasetName)
		if err != nil {
			return 0, fmt.Errorf("failed to get dataset info for %s: %w", datasetName, err)
		}

		for _, file := range dataset.Files {
			totalSize += file.Size
		}
	}

	return totalSize, nil
}

// ValidateDownloadedData validates that all downloaded data is intact
func (a *AWSDataAcquisition) ValidateDownloadedData() error {
	a.logger.Println("Validating downloaded data integrity")

	for datasetName, dataset := range AuthenticityValidationDatasets {
		datasetDir := filepath.Join(a.downloadDir, datasetName)

		for _, file := range dataset.Files {
			localPath := filepath.Join(datasetDir, filepath.Base(file.Key))

			// Check file exists
			if !a.fileExists(localPath) {
				return fmt.Errorf("missing file: %s", localPath)
			}

			// Verify file size
			if err := a.verifyFileSize(localPath, file.Size); err != nil {
				return fmt.Errorf("file validation failed for %s: %w", localPath, err)
			}
		}
	}

	a.logger.Println("All downloaded data validated successfully")
	return nil
}

// CleanupIncompleteDownloads removes incomplete or corrupted downloads
func (a *AWSDataAcquisition) CleanupIncompleteDownloads() error {
	a.logger.Println("Cleaning up incomplete downloads")

	for datasetName, dataset := range AuthenticityValidationDatasets {
		datasetDir := filepath.Join(a.downloadDir, datasetName)

		for _, file := range dataset.Files {
			localPath := filepath.Join(datasetDir, filepath.Base(file.Key))

			// Check if file exists but has wrong size
			if info, err := os.Stat(localPath); err == nil {
				if info.Size() != file.Size {
					a.logger.Printf("Removing corrupted file: %s", localPath)
					if err := os.Remove(localPath); err != nil {
						return fmt.Errorf("failed to remove corrupted file: %w", err)
					}
				}
			}
		}
	}

	return nil
}

// GetDownloadedDataPaths returns paths to downloaded datasets
func (a *AWSDataAcquisition) GetDownloadedDataPaths() map[string]string {
	paths := make(map[string]string)

	for datasetName := range AuthenticityValidationDatasets {
		paths[datasetName] = filepath.Join(a.downloadDir, datasetName)
	}

	return paths
}

// DefaultAWSDataConfig returns the default AWS data acquisition settings.
func DefaultAWSDataConfig() *AWSDataConfig {
	return &AWSDataConfig{
		Region:          "us-east-1",
		DownloadDir:     "./data/aws_genomics",
		MaxConcurrency:  10,
		MaxRetries:      3,
		VerifyChecksums: true,
		EnableProgress:  true,
		CacheData:       true,
		CleanupOnError:  true,
	}
}

// ProgressTracker tracks bytes downloaded against a total. Nothing in this
// package updates one yet.
type ProgressTracker struct {
	mu              sync.Mutex
	totalBytes      int64
	downloadedBytes int64
	startTime       time.Time
	lastUpdate      time.Time
}

// NewProgressTracker creates a new progress tracker
func NewProgressTracker(totalBytes int64) *ProgressTracker {
	return &ProgressTracker{
		totalBytes: totalBytes,
		startTime:  time.Now(),
		lastUpdate: time.Now(),
	}
}

// Update updates the progress tracker
func (p *ProgressTracker) Update(bytesDownloaded int64) {
	p.mu.Lock()
	defer p.mu.Unlock()

	p.downloadedBytes += bytesDownloaded
	p.lastUpdate = time.Now()
}

// GetProgress returns current progress information
func (p *ProgressTracker) GetProgress() *DownloadProgress {
	p.mu.Lock()
	defer p.mu.Unlock()

	elapsed := time.Since(p.startTime)
	speed := float64(p.downloadedBytes) / elapsed.Seconds() / (1024 * 1024) // MB/s

	var eta time.Duration
	if speed > 0 {
		remainingBytes := p.totalBytes - p.downloadedBytes
		eta = time.Duration(float64(remainingBytes)/speed/1024/1024) * time.Second
	}

	return &DownloadProgress{
		TotalBytes:      p.totalBytes,
		DownloadedBytes: p.downloadedBytes,
		StartTime:       p.startTime,
		Speed:           speed,
		ETA:             eta,
	}
}
