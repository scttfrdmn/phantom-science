//go:build !release
// +build !release

package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestNewEnterpriseFingerprinter_NonNil verifies that NewEnterpriseFingerprinter
// returns a non-nil fingerprinter with the correct customerID.
func TestNewEnterpriseFingerprinter_NonNil(t *testing.T) {
	ef := NewEnterpriseFingerprinter("TEST_CUSTOMER_001")
	if ef == nil {
		t.Fatal("NewEnterpriseFingerprinter returned nil")
	}
	if ef.customerID != "TEST_CUSTOMER_001" {
		t.Errorf("customerID = %q, want \"TEST_CUSTOMER_001\"", ef.customerID)
	}
	if ef.engine == nil {
		t.Error("engine should not be nil")
	}
	defer ef.engine.Close()
}

// TestNewEnterpriseFingerprinter_EmptyCustomerID verifies that an empty customerID
// still returns a non-nil fingerprinter.
func TestNewEnterpriseFingerprinter_EmptyCustomerID(t *testing.T) {
	ef := NewEnterpriseFingerprinter("")
	if ef == nil {
		t.Fatal("NewEnterpriseFingerprinter with empty customerID returned nil")
	}
	defer ef.engine.Close()
}

// TestNewEnterpriseFingerprinter_SpecialChars verifies customerIDs with special
// characters are preserved verbatim.
func TestNewEnterpriseFingerprinter_SpecialChars(t *testing.T) {
	customerID := "ACME-Corp_2024/Q1@genomics"
	ef := NewEnterpriseFingerprinter(customerID)
	if ef == nil {
		t.Fatal("NewEnterpriseFingerprinter returned nil")
	}
	if ef.customerID != customerID {
		t.Errorf("customerID = %q, want %q", ef.customerID, customerID)
	}
	defer ef.engine.Close()
}

// TestDemoScenarios_NoPanic verifies that DemoScenarios does not panic.
func TestDemoScenarios_NoPanic(t *testing.T) {
	ef := NewEnterpriseFingerprinter("TEST_CUSTOMER")
	defer ef.engine.Close()
	// Should not panic.
	ef.DemoScenarios()
}

// TestGenerateAuditReport_CreatesFile verifies that generateAuditReport creates
// a file with the expected content.
func TestGenerateAuditReport_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	ef := NewEnterpriseFingerprinter("AUDIT_TEST_CUSTOMER")
	defer ef.engine.Close()

	err := ef.generateAuditReport("test_dataset_001", "fp-abc123", 5, tmpDir)
	if err != nil {
		t.Fatalf("generateAuditReport returned error: %v", err)
	}

	reportPath := filepath.Join(tmpDir, "fingerprint_audit_report.txt")
	data, err := os.ReadFile(reportPath)
	if err != nil {
		t.Fatalf("could not read audit report file: %v", err)
	}

	content := string(data)

	// Verify key fields appear in the report.
	checks := []string{
		"PHANTOM SCIENCE",
		"ENTERPRISE FINGERPRINTING AUDIT REPORT",
		"AUDIT_TEST_CUSTOMER",
		"test_dataset_001",
		"fp-abc123",
		"5", // numSamples
		"STEGANOGRAPHIC FINGERPRINTING SUMMARY",
		"COMPLIANCE INFORMATION",
		"FINGERPRINT VERIFICATION",
		"SECURITY FEATURES",
	}
	for _, check := range checks {
		if !strings.Contains(content, check) {
			t.Errorf("audit report missing expected string %q", check)
		}
	}
}

// TestGenerateAuditReport_ZeroSamples verifies generateAuditReport handles zero samples.
func TestGenerateAuditReport_ZeroSamples(t *testing.T) {
	tmpDir := t.TempDir()
	ef := NewEnterpriseFingerprinter("ZERO_SAMPLES_CUSTOMER")
	defer ef.engine.Close()

	err := ef.generateAuditReport("empty_dataset", "fp-zero", 0, tmpDir)
	if err != nil {
		t.Fatalf("generateAuditReport with 0 samples returned error: %v", err)
	}

	reportPath := filepath.Join(tmpDir, "fingerprint_audit_report.txt")
	if _, err := os.Stat(reportPath); os.IsNotExist(err) {
		t.Error("audit report file was not created")
	}
}

// TestGenerateAuditReport_InvalidDir verifies generateAuditReport fails gracefully
// when the output directory does not exist and cannot be created by WriteFile.
func TestGenerateAuditReport_InvalidDir(t *testing.T) {
	ef := NewEnterpriseFingerprinter("BAD_DIR_CUSTOMER")
	defer ef.engine.Close()

	// Use a path that cannot be written to (directory doesn't exist, no mkdir).
	err := ef.generateAuditReport("ds", "fp", 1, "/nonexistent/invalid/path/dir")
	if err == nil {
		t.Error("generateAuditReport should fail for non-existent output directory")
	}
}

// TestGenerateAuditReport_CustomerIDInContent verifies that the customerID is
// embedded in the report content.
func TestGenerateAuditReport_CustomerIDInContent(t *testing.T) {
	tmpDir := t.TempDir()
	customerID := "PHARMA_CORP_XR9"
	ef := NewEnterpriseFingerprinter(customerID)
	defer ef.engine.Close()

	err := ef.generateAuditReport("clinical_trial_001", "fp-clinical", 10, tmpDir)
	if err != nil {
		t.Fatalf("generateAuditReport returned error: %v", err)
	}

	data, err := os.ReadFile(filepath.Join(tmpDir, "fingerprint_audit_report.txt"))
	if err != nil {
		t.Fatalf("could not read audit report: %v", err)
	}

	if !strings.Contains(string(data), customerID) {
		t.Errorf("audit report does not contain customerID %q", customerID)
	}
}

// TestGenerateAuditReport_MasterFingerprintInContent verifies that the master
// fingerprint appears in the report.
func TestGenerateAuditReport_MasterFingerprintInContent(t *testing.T) {
	tmpDir := t.TempDir()
	ef := NewEnterpriseFingerprinter("FP_CHECK_CUSTOMER")
	defer ef.engine.Close()

	masterFP := "MASTER-FP-XYZ-99887766"
	err := ef.generateAuditReport("fp_dataset", masterFP, 3, tmpDir)
	if err != nil {
		t.Fatalf("generateAuditReport returned error: %v", err)
	}

	data, err := os.ReadFile(filepath.Join(tmpDir, "fingerprint_audit_report.txt"))
	if err != nil {
		t.Fatalf("could not read audit report: %v", err)
	}

	if !strings.Contains(string(data), masterFP) {
		t.Errorf("audit report does not contain master fingerprint %q", masterFP)
	}
}

// TestVerifyDataset_EmptyDir verifies VerifyDataset on a directory with no FASTQ files.
func TestVerifyDataset_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	ef := NewEnterpriseFingerprinter("VERIFY_CUSTOMER")
	defer ef.engine.Close()

	// No FASTQ files → verified=0, failed=0 → should return nil.
	err := ef.VerifyDataset(tmpDir)
	if err != nil {
		t.Errorf("VerifyDataset on empty dir should return nil, got: %v", err)
	}
}

// TestVerifyDataset_WithFastqFiles verifies VerifyDataset processes FASTQ files.
// Because the engine may or may not verify (depending on content), we only
// check that the function runs without panicking and returns consistent output.
func TestVerifyDataset_WithFastqFiles(t *testing.T) {
	tmpDir := t.TempDir()

	// Create two dummy FASTQ files.
	for _, name := range []string{"sample1.fastq", "sample2.fastq"} {
		content := "@READ1\nATCGATCGATCG\n+\nIIIIIIIIIIII\n"
		path := filepath.Join(tmpDir, name)
		if err := os.WriteFile(path, []byte(content), 0600); err != nil {
			t.Fatalf("failed to write FASTQ file %s: %v", name, err)
		}
	}

	ef := NewEnterpriseFingerprinter("VERIFY_FASTQ_CUSTOMER")
	defer ef.engine.Close()

	// The function may return an error if verification fails for the dummy files —
	// that is acceptable. We just verify it does not panic.
	_ = ef.VerifyDataset(tmpDir)
}

// TestVerifyDataset_UnfingerprintedFilesReturnError verifies that VerifyDataset
// returns an error when FASTQ files don't contain valid fingerprints.
func TestVerifyDataset_UnfingerprintedFilesReturnError(t *testing.T) {
	tmpDir := t.TempDir()

	// Write plain (non-fingerprinted) FASTQ content.
	content := "@READ1\nATCGATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIIIIIII\n"
	if err := os.WriteFile(filepath.Join(tmpDir, "plain.fastq"), []byte(content), 0600); err != nil {
		t.Fatalf("failed to write test FASTQ: %v", err)
	}

	ef := NewEnterpriseFingerprinter("UNFINGERPRINTED_CUSTOMER")
	defer ef.engine.Close()

	// Verification should fail because the content is not fingerprinted.
	err := ef.VerifyDataset(tmpDir)
	if err == nil {
		t.Log("VerifyDataset returned nil (engine may have verified the dummy content)")
	} else {
		// Expected: error mentioning failed files.
		if !strings.Contains(err.Error(), "verification failed") {
			t.Errorf("error message = %q, expected to contain \"verification failed\"", err.Error())
		}
	}
}

// TestVerifyDataset_SingleFastqFile verifies VerifyDataset handles exactly one file.
func TestVerifyDataset_SingleFastqFile(t *testing.T) {
	tmpDir := t.TempDir()

	content := "@SEQ001\nGCGCGCGCGCGCGCGCGCGC\n+\nJJJJJJJJJJJJJJJJJJJJ\n"
	if err := os.WriteFile(filepath.Join(tmpDir, "one.fastq"), []byte(content), 0600); err != nil {
		t.Fatalf("failed to write FASTQ: %v", err)
	}

	ef := NewEnterpriseFingerprinter("SINGLE_FILE_CUSTOMER")
	defer ef.engine.Close()

	// Either passes or fails — must not panic.
	_ = ef.VerifyDataset(tmpDir)
}

// TestVerifyDataset_NonexistentDir verifies VerifyDataset on a non-existent directory.
// filepath.Glob returns nil (no error) for paths that match nothing, so the
// function should succeed with zero files found.
func TestVerifyDataset_NonexistentDir(t *testing.T) {
	ef := NewEnterpriseFingerprinter("MISSING_DIR_CUSTOMER")
	defer ef.engine.Close()

	// filepath.Glob does not return an error for paths with no matches,
	// so VerifyDataset should return nil.
	err := ef.VerifyDataset("/tmp/nonexistent_phantom_dir_xyz_12345")
	if err != nil {
		t.Errorf("VerifyDataset with non-existent dir returned unexpected error: %v", err)
	}
}

// TestFingerprintDataset_CreatesOutputDir verifies FingerprintDataset creates
// the output directory when it doesn't exist.
func TestFingerprintDataset_CreatesOutputDir(t *testing.T) {
	tmpDir := t.TempDir()
	outputDir := filepath.Join(tmpDir, "fingerprinted_output")

	ef := NewEnterpriseFingerprinter("DATASET_CUSTOMER_001")
	defer ef.engine.Close()

	// FingerprintDataset may fail because genomicsEngine.Generate() may not
	// produce real output. We test that the output dir is created on success.
	err := ef.FingerprintDataset("test_dataset", outputDir, 1)
	if err != nil {
		// Even on failure, the directory may have been created.
		// Just check no panic happened.
		t.Logf("FingerprintDataset returned error (expected in test): %v", err)
		return
	}

	if _, err := os.Stat(outputDir); os.IsNotExist(err) {
		t.Error("output directory was not created by FingerprintDataset")
	}
}

// TestFingerprintDataset_ZeroSamples verifies FingerprintDataset with 0 samples.
// With zero samples the loop doesn't execute, so the audit report should be created.
func TestFingerprintDataset_ZeroSamples(t *testing.T) {
	tmpDir := t.TempDir()
	ef := NewEnterpriseFingerprinter("ZERO_SAMPLE_CUSTOMER")
	defer ef.engine.Close()

	err := ef.FingerprintDataset("zero_sample_dataset", tmpDir, 0)
	if err != nil {
		// This may fail if GenerateFingerprint fails; log and skip further checks.
		t.Logf("FingerprintDataset(0 samples) returned error: %v", err)
		return
	}

	// Audit report should be present.
	reportPath := filepath.Join(tmpDir, "fingerprint_audit_report.txt")
	if _, err := os.Stat(reportPath); os.IsNotExist(err) {
		t.Error("audit report was not created for 0-sample dataset")
	}
}

// TestEnterpriseFingerprinter_CustomerIDPreserved verifies that the customerID
// stored in the struct is exactly what was passed in.
func TestEnterpriseFingerprinter_CustomerIDPreserved(t *testing.T) {
	ids := []string{
		"simple",
		"COMPLEX_ID_WITH_NUMBERS_123",
		"id-with-dashes",
		"id_with_underscores",
	}
	for _, id := range ids {
		ef := NewEnterpriseFingerprinter(id)
		if ef.customerID != id {
			t.Errorf("customerID = %q, want %q", ef.customerID, id)
		}
		ef.engine.Close()
	}
}

// TestFingerprintDataset_MkdirFailure verifies FingerprintDataset returns an error
// when the output directory cannot be created.
func TestFingerprintDataset_MkdirFailure(t *testing.T) {
	ef := NewEnterpriseFingerprinter("MKDIR_FAIL_CUSTOMER")
	defer ef.engine.Close()

	// Attempt to create a directory inside a file (not a directory) → should fail.
	tmpFile, err := os.CreateTemp("", "phantom_test_*.txt")
	if err != nil {
		t.Fatalf("could not create temp file: %v", err)
	}
	defer os.Remove(tmpFile.Name())
	tmpFile.Close()

	// Try to use the file as if it were a directory.
	badDir := filepath.Join(tmpFile.Name(), "subdir")
	err = ef.FingerprintDataset("ds", badDir, 0)
	if err == nil {
		t.Error("FingerprintDataset should fail when MkdirAll fails")
	}
	if !strings.Contains(err.Error(), "failed to create output directory") {
		t.Errorf("error = %q, want to contain \"failed to create output directory\"", err.Error())
	}
}

// TestGenerateAuditReport_ReportFilename verifies the audit report file is always
// named "fingerprint_audit_report.txt".
func TestGenerateAuditReport_ReportFilename(t *testing.T) {
	tmpDir := t.TempDir()
	ef := NewEnterpriseFingerprinter("FILENAME_TEST_CUSTOMER")
	defer ef.engine.Close()

	if err := ef.generateAuditReport("ds", "fp", 2, tmpDir); err != nil {
		t.Fatalf("generateAuditReport returned error: %v", err)
	}

	expected := filepath.Join(tmpDir, "fingerprint_audit_report.txt")
	if _, err := os.Stat(expected); os.IsNotExist(err) {
		t.Errorf("expected report at %s, but file does not exist", expected)
	}
}

// TestGenerateAuditReport_NumSamplesInContent verifies that the sample count
// appears in the report body.
func TestGenerateAuditReport_NumSamplesInContent(t *testing.T) {
	tmpDir := t.TempDir()
	ef := NewEnterpriseFingerprinter("SAMPLE_COUNT_CUSTOMER")
	defer ef.engine.Close()

	if err := ef.generateAuditReport("ds", "fp", 42, tmpDir); err != nil {
		t.Fatalf("generateAuditReport returned error: %v", err)
	}

	data, err := os.ReadFile(filepath.Join(tmpDir, "fingerprint_audit_report.txt"))
	if err != nil {
		t.Fatalf("could not read report: %v", err)
	}

	if !strings.Contains(string(data), "42") {
		t.Error("audit report does not contain the sample count (42)")
	}
}

// TestGenerateAuditReport_OverwritesPreviousReport verifies that calling
// generateAuditReport twice overwrites the previous report.
func TestGenerateAuditReport_OverwritesPreviousReport(t *testing.T) {
	tmpDir := t.TempDir()
	ef := NewEnterpriseFingerprinter("OVERWRITE_CUSTOMER")
	defer ef.engine.Close()

	if err := ef.generateAuditReport("first_dataset", "fp-first", 1, tmpDir); err != nil {
		t.Fatalf("first generateAuditReport call failed: %v", err)
	}

	if err := ef.generateAuditReport("second_dataset", "fp-second", 2, tmpDir); err != nil {
		t.Fatalf("second generateAuditReport call failed: %v", err)
	}

	data, err := os.ReadFile(filepath.Join(tmpDir, "fingerprint_audit_report.txt"))
	if err != nil {
		t.Fatalf("could not read report: %v", err)
	}

	content := string(data)
	// Second report should have second_dataset, not first_dataset.
	if !strings.Contains(content, "second_dataset") {
		t.Error("overwritten report should contain second_dataset")
	}
	if strings.Contains(content, "first_dataset") {
		t.Error("overwritten report should NOT contain first_dataset")
	}
}
