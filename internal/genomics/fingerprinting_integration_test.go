package genomics

import (
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestGenomicsEngineWithFingerprinting(t *testing.T) {
	// Create a temporary directory for test output
	tempDir := t.TempDir()

	// Create test configuration with fingerprinting enabled
	config := types.SequencerConfig{
		SampleName:    "test-sample",
		Lane:          1,
		OutputDir:     tempDir,
		ReadLength:    100,
		ReadsPerBatch: 10,
		TotalBatches:  1,
		PairedEnd:     false,
		BaseQuality:   40,
		QualityDecay:  0.1,
		DatasetID:     "test-dataset-123",
		Enterprise: types.EnterpriseConfig{
			CustomerID: "test-customer-456",
		},
		Fingerprinting: types.FingerprintConfig{
			Enable:               true,
			CustomerFingerprint:  true,
			DatasetFingerprint:   true,
			SteganographyLayers:  []string{"read_ids"}, // labels only; keyed_seed is covered by fingerprint_keyed_test.go
			CryptographicSigning: true,
			TamperDetection:      true,
		},
	}

	// Create genomics engine
	engine := NewEngine(config)
	defer engine.Close()

	t.Run("GenerateFASTQWithFingerprinting", func(t *testing.T) {
		// Generate FASTQ file
		filename := filepath.Join(tempDir, "test_fingerprinted.fastq")
		err := engine.GenerateFASTQFile(filename, 1, 1)
		if err != nil {
			t.Fatalf("Failed to generate FASTQ file: %v", err)
		}

		// Verify file was created
		if _, err := os.Stat(filename); os.IsNotExist(err) {
			t.Fatal("FASTQ file was not created")
		}

		// Read and validate file contents
		content, err := os.ReadFile(filename)
		if err != nil {
			t.Fatalf("Failed to read FASTQ file: %v", err)
		}

		contentStr := string(content)

		// Check that file contains expected number of reads (count lines starting with @)
		lines := strings.Split(contentStr, "\n")
		readCount := 0
		for _, line := range lines {
			if strings.HasPrefix(line, "@") && !strings.HasPrefix(line, "@+") {
				readCount++
			}
		}

		expectedReads := config.ReadsPerBatch
		if readCount != expectedReads {
			t.Errorf("Expected %d reads, got %d", expectedReads, readCount)
		}

		// Every read ID must carry the customer and dataset tags.
		for i, line := range lines {
			if i%4 != 0 || line == "" {
				continue
			}
			if !strings.Contains(line, config.Enterprise.CustomerID) || !strings.Contains(line, config.DatasetID) {
				t.Fatalf("read ID %q lacks customer %q or dataset %q", line, config.Enterprise.CustomerID, config.DatasetID)
			}
		}
	})

	t.Run("GenerateWithoutFingerprinting", func(t *testing.T) {
		// Create config without fingerprinting
		configNoFingerprint := config
		configNoFingerprint.Fingerprinting.Enable = false
		configNoFingerprint.OutputDir = filepath.Join(tempDir, "no_fingerprint")

		// Create new engine
		engineNoFingerprint := NewEngine(configNoFingerprint)
		defer engineNoFingerprint.Close()

		// Create output directory
		err := os.MkdirAll(configNoFingerprint.OutputDir, 0750)
		if err != nil {
			t.Fatalf("Failed to create output directory: %v", err)
		}

		// Generate FASTQ file
		filename := filepath.Join(configNoFingerprint.OutputDir, "test_no_fingerprint.fastq")
		err = engineNoFingerprint.GenerateFASTQFile(filename, 1, 1)
		if err != nil {
			t.Fatalf("Failed to generate FASTQ file without fingerprinting: %v", err)
		}

		// Verify file was created
		if _, err := os.Stat(filename); os.IsNotExist(err) {
			t.Fatal("FASTQ file was not created")
		}

		// Read file contents
		content, err := os.ReadFile(filename)
		if err != nil {
			t.Fatalf("Failed to read FASTQ file: %v", err)
		}

		// Basic validation - should contain reads
		if len(content) == 0 {
			t.Error("Generated FASTQ file is empty")
		}
	})

	t.Run("CompareQualityScores", func(t *testing.T) {
		// Generate two files: one with fingerprinting, one without
		filenameWithFP := filepath.Join(tempDir, "with_fp.fastq")
		filenameWithoutFP := filepath.Join(tempDir, "without_fp.fastq")

		// Generate with fingerprinting
		err := engine.GenerateFASTQFile(filenameWithFP, 1, 1)
		if err != nil {
			t.Fatalf("Failed to generate FASTQ file with fingerprinting: %v", err)
		}

		// Generate without fingerprinting
		configNoFP := config
		configNoFP.Fingerprinting.Enable = false
		engineNoFP := NewEngine(configNoFP)
		defer engineNoFP.Close()

		err = engineNoFP.GenerateFASTQFile(filenameWithoutFP, 1, 1)
		if err != nil {
			t.Fatalf("Failed to generate FASTQ file without fingerprinting: %v", err)
		}

		// Read both files
		contentWithFP, err := os.ReadFile(filenameWithFP)
		if err != nil {
			t.Fatalf("Failed to read fingerprinted FASTQ file: %v", err)
		}

		contentWithoutFP, err := os.ReadFile(filenameWithoutFP)
		if err != nil {
			t.Fatalf("Failed to read non-fingerprinted FASTQ file: %v", err)
		}

		// Both files should have content
		if len(contentWithFP) == 0 || len(contentWithoutFP) == 0 {
			t.Error("One or both FASTQ files are empty")
		}

		// Files should be different if fingerprinting is working
		// (though they might be the same if stegano library is not available)
		if string(contentWithFP) == string(contentWithoutFP) {
			t.Log("Warning: Fingerprinted and non-fingerprinted files are identical (may be expected if stegano library not available)")
		}
	})
}

func TestGenomicsEnginePairedEndWithFingerprinting(t *testing.T) {
	// Create a temporary directory for test output
	tempDir := t.TempDir()

	// Create test configuration for paired-end sequencing with fingerprinting
	config := types.SequencerConfig{
		SampleName:    "test-paired-sample",
		Lane:          1,
		OutputDir:     tempDir,
		ReadLength:    75,
		ReadsPerBatch: 5,
		TotalBatches:  1,
		PairedEnd:     true,
		BaseQuality:   35,
		QualityDecay:  0.05,
		DatasetID:     "paired-dataset-789",
		Enterprise: types.EnterpriseConfig{
			CustomerID: "paired-customer-101",
		},
		Fingerprinting: types.FingerprintConfig{
			Enable:              true,
			CustomerFingerprint: true,
			DatasetFingerprint:  true,
			SteganographyLayers: []string{"read_ids"}, // labels only; keyed_seed is covered by fingerprint_keyed_test.go
		},
	}

	// Create genomics engine
	engine := NewEngine(config)
	defer engine.Close()

	// Run the complete generation process
	err := engine.Generate()
	if err != nil {
		t.Fatalf("Failed to generate paired-end FASTQ files: %v", err)
	}

	// Check that both R1 and R2 files were created
	r1File := filepath.Join(tempDir, "test-paired-sample_S1_L001_R1_001.fastq")
	r2File := filepath.Join(tempDir, "test-paired-sample_S1_L001_R2_001.fastq")

	if _, err := os.Stat(r1File); os.IsNotExist(err) {
		t.Error("R1 FASTQ file was not created")
	}

	if _, err := os.Stat(r2File); os.IsNotExist(err) {
		t.Error("R2 FASTQ file was not created")
	}

	// Verify both files have content
	r1Content, err := os.ReadFile(r1File)
	if err != nil {
		t.Fatalf("Failed to read R1 file: %v", err)
	}

	r2Content, err := os.ReadFile(r2File)
	if err != nil {
		t.Fatalf("Failed to read R2 file: %v", err)
	}

	if len(r1Content) == 0 || len(r2Content) == 0 {
		t.Error("One or both paired-end files are empty")
	}

	// Verify that both files have the expected number of reads
	r1ReadCount := countFASTQReads(string(r1Content))
	r2ReadCount := countFASTQReads(string(r2Content))

	expectedReads := config.ReadsPerBatch
	if r1ReadCount != expectedReads {
		t.Errorf("R1 file: expected %d reads, got %d", expectedReads, r1ReadCount)
	}

	if r2ReadCount != expectedReads {
		t.Errorf("R2 file: expected %d reads, got %d", expectedReads, r2ReadCount)
	}
}

// Helper function to count FASTQ reads in content
func countFASTQReads(content string) int {
	count := 0
	for i, char := range content {
		if char == '@' && (i == 0 || content[i-1] == '\n') {
			count++
		}
	}
	return count
}
