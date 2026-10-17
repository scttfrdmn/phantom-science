//go:build !release
// +build !release

package main

import (
	"os"
	"path/filepath"
	"testing"
)

func TestAnonymizeSequence_WithHeader(t *testing.T) {
	content := []byte(">SEQ001 population=YRI length=100 gc=0.42\nATCGATCGATCG")
	result := anonymizeSequence(content)
	resultStr := string(result)
	// Should not contain population info
	if contains(resultStr, "population=YRI") {
		t.Error("anonymizeSequence should remove population info")
	}
	// Should start with ANONYMOUS
	if !contains(resultStr, "ANONYMOUS") {
		t.Error("anonymizeSequence should add ANONYMOUS header")
	}
}

func TestAnonymizeSequence_NoHeader(t *testing.T) {
	content := []byte("ATCGATCGATCG")
	result := anonymizeSequence(content)
	// Should be unchanged (no header to anonymize)
	if string(result) != "ATCGATCGATCG" {
		t.Errorf("anonymizeSequence without header modified content: %q", string(result))
	}
}

func TestExtractLength_WithLength(t *testing.T) {
	header := ">SEQ001 population=YRI length=500 gc=0.42"
	got := extractLength(header)
	if got != "500" {
		t.Errorf("extractLength = %q, want \"500\"", got)
	}
}

func TestExtractLength_NoLength(t *testing.T) {
	header := ">SEQ001 population=YRI gc=0.42"
	got := extractLength(header)
	if got != "unknown" {
		t.Errorf("extractLength (no length field) = %q, want \"unknown\"", got)
	}
}

func contains(s, substr string) bool {
	return len(s) >= len(substr) && (s == substr || len(s) > 0 && containsHelper(s, substr))
}

func containsHelper(s, substr string) bool {
	for i := 0; i <= len(s)-len(substr); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}

// ---------- createPackageStructure ----------

func TestCreatePackageStructure_CreatesSubdirs(t *testing.T) {
	tmpDir := t.TempDir()
	if err := createPackageStructure(tmpDir); err != nil {
		t.Fatalf("createPackageStructure error: %v", err)
	}
	for _, subdir := range []string{"blind_test_sequences", "reference_sequences", "evaluation_materials", "results"} {
		info, err := os.Stat(filepath.Join(tmpDir, subdir))
		if err != nil {
			t.Errorf("expected subdir %s to exist: %v", subdir, err)
		} else if !info.IsDir() {
			t.Errorf("expected %s to be a directory", subdir)
		}
	}
}

func TestCreatePackageStructure_BadPath(t *testing.T) {
	// Passing a path under a non-existent parent should fail
	err := createPackageStructure("/nonexistent/path/xyz/output")
	if err == nil {
		t.Error("expected error for bad path")
	}
}

// ---------- createAnonymizedSequences ----------

func TestCreateAnonymizedSequences_EmptyPopulations(t *testing.T) {
	tmpDir := t.TempDir()
	_ = createPackageStructure(tmpDir)
	count, err := createAnonymizedSequences("/nonexistent/src", tmpDir, []string{})
	if err != nil {
		t.Fatalf("createAnonymizedSequences(no populations) error: %v", err)
	}
	if count != 0 {
		t.Errorf("expected 0 sequences, got %d", count)
	}
}

func TestCreateAnonymizedSequences_MissingPopDirs(t *testing.T) {
	tmpDir := t.TempDir()
	_ = createPackageStructure(tmpDir)
	// Source dir doesn't exist, populations skip with warning
	count, err := createAnonymizedSequences(tmpDir+"/nosrc", tmpDir, []string{"YRI", "CHB"})
	if err != nil {
		t.Fatalf("createAnonymizedSequences(missing pop dirs) error: %v", err)
	}
	if count != 0 {
		t.Errorf("expected 0 sequences (no pop dirs), got %d", count)
	}
}

func TestCreateAnonymizedSequences_WithFastaFiles(t *testing.T) {
	tmpDir := t.TempDir()
	_ = createPackageStructure(tmpDir)

	// Create fake source with population dir + FASTA file
	srcDir := filepath.Join(tmpDir, "src")
	popDir := filepath.Join(srcDir, "YRI")
	_ = os.MkdirAll(popDir, 0750)
	fastaContent := ">YRI_S01 population=YRI length=20\nATCGATCGATCGATCGATCG\n"
	_ = os.WriteFile(filepath.Join(popDir, "sample.fasta"), []byte(fastaContent), 0600)

	count, err := createAnonymizedSequences(srcDir, tmpDir, []string{"YRI"})
	if err != nil {
		t.Fatalf("createAnonymizedSequences(with files) error: %v", err)
	}
	if count != 1 {
		t.Errorf("expected 1 sequence, got %d", count)
	}
}

// ---------- copyAndAnonymizeFile ----------

func TestCopyAndAnonymizeFile_Success(t *testing.T) {
	tmpDir := t.TempDir()
	srcPath := filepath.Join(tmpDir, "original.fasta")
	dstPath := filepath.Join(tmpDir, "anon.fasta")
	content := ">SEQ001 population=YRI length=20\nATCGATCGATCGATCGATCG\n"
	_ = os.WriteFile(srcPath, []byte(content), 0600)

	if err := copyAndAnonymizeFile(srcPath, dstPath); err != nil {
		t.Fatalf("copyAndAnonymizeFile error: %v", err)
	}
	if _, err := os.Stat(dstPath); err != nil {
		t.Fatal("expected output file to exist")
	}
}

func TestCopyAndAnonymizeFile_BadSrc(t *testing.T) {
	err := copyAndAnonymizeFile("/nonexistent/src.fasta", "/tmp/dst.fasta")
	if err == nil {
		t.Error("expected error for bad source path")
	}
}

// ---------- copyFile ----------

func TestCopyFile_Success(t *testing.T) {
	tmpDir := t.TempDir()
	srcPath := filepath.Join(tmpDir, "src.fasta")
	dstPath := filepath.Join(tmpDir, "dst.fasta")
	_ = os.WriteFile(srcPath, []byte("ATCG"), 0600)

	if err := copyFile(srcPath, dstPath); err != nil {
		t.Fatalf("copyFile error: %v", err)
	}
	data, _ := os.ReadFile(dstPath)
	if string(data) != "ATCG" {
		t.Errorf("copyFile content = %q, want \"ATCG\"", string(data))
	}
}

func TestCopyFile_BadSrc(t *testing.T) {
	err := copyFile("/nonexistent/src.fasta", "/tmp/dst_test.fasta")
	if err == nil {
		t.Error("expected error for bad source path")
	}
}

// ---------- createExpertMaterials ----------

func TestCreateExpertMaterials_Success(t *testing.T) {
	tmpDir := t.TempDir()
	_ = createPackageStructure(tmpDir)
	if err := createExpertMaterials(tmpDir, 42); err != nil {
		t.Fatalf("createExpertMaterials error: %v", err)
	}
}

// ---------- createEvaluationProtocol ----------

func TestCreateEvaluationProtocol_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	if err := createEvaluationProtocol(tmpDir, 10); err != nil {
		t.Fatalf("createEvaluationProtocol error: %v", err)
	}
	if _, err := os.Stat(filepath.Join(tmpDir, "EVALUATION_PROTOCOL.md")); err != nil {
		t.Error("EVALUATION_PROTOCOL.md not created")
	}
}

// ---------- createExpertREADME ----------

func TestCreateExpertREADME_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	if err := createExpertREADME(tmpDir); err != nil {
		t.Fatalf("createExpertREADME error: %v", err)
	}
	if _, err := os.Stat(filepath.Join(tmpDir, "README_FOR_EXPERTS.md")); err != nil {
		t.Error("README_FOR_EXPERTS.md not created")
	}
}
