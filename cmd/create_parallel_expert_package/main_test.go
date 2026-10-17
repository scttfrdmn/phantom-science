//go:build !release
// +build !release

package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestAnonymizeParallelSequence_WithHeader(t *testing.T) {
	content := []byte(">SEQ001 population=YRI length=100 complexity=0.85 gc_content=0.42 vcf_trained=true\nATCGATCGATCG")
	result := anonymizeParallelSequence(content)
	resultStr := string(result)
	// Should not contain population info
	if strings.Contains(resultStr, "population=YRI") {
		t.Error("anonymizeParallelSequence should remove population info")
	}
	// Should start with ANONYMOUS_PARALLEL_VCF_SEQUENCE
	if !strings.Contains(resultStr, "ANONYMOUS_PARALLEL_VCF_SEQUENCE") {
		t.Error("anonymizeParallelSequence should add ANONYMOUS_PARALLEL_VCF_SEQUENCE header")
	}
}

func TestAnonymizeParallelSequence_NoHeader(t *testing.T) {
	content := []byte("ATCGATCGATCG")
	result := anonymizeParallelSequence(content)
	if string(result) != "ATCGATCGATCG" {
		t.Errorf("anonymizeParallelSequence without header modified content: %q", string(result))
	}
}

func TestExtractMetadata_WithKey(t *testing.T) {
	header := ">SEQ001 population=YRI length=500 gc_content=0.42"
	got := extractMetadata(header, "length")
	if got != "500" {
		t.Errorf("extractMetadata(length) = %q, want \"500\"", got)
	}
}

func TestExtractMetadata_WithGCContent(t *testing.T) {
	header := ">SEQ001 population=YRI length=500 gc_content=0.42"
	got := extractMetadata(header, "gc_content")
	if got != "0.42" {
		t.Errorf("extractMetadata(gc_content) = %q, want \"0.42\"", got)
	}
}

func TestExtractMetadata_MissingKey(t *testing.T) {
	header := ">SEQ001 population=YRI"
	got := extractMetadata(header, "length")
	if got != "unknown" {
		t.Errorf("extractMetadata(missing key) = %q, want \"unknown\"", got)
	}
}

// ---------- createPackageStructure ----------

func TestCreatePackageStructure_CreatesSubdirs(t *testing.T) {
	tmpDir := t.TempDir()
	if err := createPackageStructure(tmpDir); err != nil {
		t.Fatalf("createPackageStructure error: %v", err)
	}
	for _, subdir := range []string{
		"blind_test_sequences", "reference_sequences",
		"evaluation_materials", "results", "training_documentation",
	} {
		info, err := os.Stat(filepath.Join(tmpDir, subdir))
		if err != nil {
			t.Errorf("expected subdir %s to exist: %v", subdir, err)
		} else if !info.IsDir() {
			t.Errorf("expected %s to be a directory", subdir)
		}
	}
}

func TestCreatePackageStructure_BadPath(t *testing.T) {
	err := createPackageStructure("/nonexistent/path/xyz/output")
	if err == nil {
		t.Error("expected error for bad path")
	}
}

// ---------- createParallelAnonymizedSequences ----------

func TestCreateParallelAnonymizedSequences_EmptyPopulations(t *testing.T) {
	tmpDir := t.TempDir()
	_ = createPackageStructure(tmpDir)
	count, err := createParallelAnonymizedSequences("/nonexistent/src", tmpDir, []string{})
	if err != nil {
		t.Fatalf("error: %v", err)
	}
	if count != 0 {
		t.Errorf("expected 0 sequences, got %d", count)
	}
}

func TestCreateParallelAnonymizedSequences_MissingPopDirs(t *testing.T) {
	tmpDir := t.TempDir()
	_ = createPackageStructure(tmpDir)
	count, err := createParallelAnonymizedSequences(tmpDir+"/nosrc", tmpDir, []string{"YRI", "CHB"})
	if err != nil {
		t.Fatalf("error: %v", err)
	}
	if count != 0 {
		t.Errorf("expected 0 sequences (no pop dirs), got %d", count)
	}
}

func TestCreateParallelAnonymizedSequences_WithFastaFiles(t *testing.T) {
	tmpDir := t.TempDir()
	_ = createPackageStructure(tmpDir)

	srcDir := filepath.Join(tmpDir, "src")
	popDir := filepath.Join(srcDir, "YRI")
	_ = os.MkdirAll(popDir, 0750)
	fastaContent := ">YRI_S01 population=YRI length=20 complexity=0.80 gc_content=0.42 vcf_trained=true\nATCGATCGATCGATCGATCG\n"
	_ = os.WriteFile(filepath.Join(popDir, "sample.fasta"), []byte(fastaContent), 0600)

	count, err := createParallelAnonymizedSequences(srcDir, tmpDir, []string{"YRI"})
	if err != nil {
		t.Fatalf("error: %v", err)
	}
	if count != 1 {
		t.Errorf("expected 1 sequence, got %d", count)
	}
}

// ---------- copyAndAnonymizeParallelFile ----------

func TestCopyAndAnonymizeParallelFile_Success(t *testing.T) {
	tmpDir := t.TempDir()
	srcPath := filepath.Join(tmpDir, "original.fasta")
	dstPath := filepath.Join(tmpDir, "anon.fasta")
	content := ">SEQ001 population=YRI length=20\nATCGATCGATCGATCGATCG\n"
	_ = os.WriteFile(srcPath, []byte(content), 0600)

	if err := copyAndAnonymizeParallelFile(srcPath, dstPath); err != nil {
		t.Fatalf("error: %v", err)
	}
	if _, err := os.Stat(dstPath); err != nil {
		t.Fatal("expected output file to exist")
	}
}

func TestCopyAndAnonymizeParallelFile_BadSrc(t *testing.T) {
	err := copyAndAnonymizeParallelFile("/nonexistent/src.fasta", "/tmp/dst.fasta")
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
		t.Fatalf("error: %v", err)
	}
	data, _ := os.ReadFile(dstPath)
	if string(data) != "ATCG" {
		t.Errorf("content = %q, want \"ATCG\"", string(data))
	}
}

func TestCopyFile_BadSrc(t *testing.T) {
	err := copyFile("/nonexistent/src.fasta", "/tmp/dst_parallel_test.fasta")
	if err == nil {
		t.Error("expected error for bad source path")
	}
}

// ---------- createParallelEvaluationProtocol ----------

func TestCreateParallelEvaluationProtocol_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	if err := createParallelEvaluationProtocol(tmpDir, 10); err != nil {
		t.Fatalf("error: %v", err)
	}
	if _, err := os.Stat(filepath.Join(tmpDir, "PARALLEL_VCF_EVALUATION_PROTOCOL.md")); err != nil {
		t.Error("PARALLEL_VCF_EVALUATION_PROTOCOL.md not created")
	}
}

// ---------- createParallelExpertREADME ----------

func TestCreateParallelExpertREADME_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	if err := createParallelExpertREADME(tmpDir); err != nil {
		t.Fatalf("error: %v", err)
	}
	if _, err := os.Stat(filepath.Join(tmpDir, "README_PARALLEL_VCF_EXPERT_VALIDATION.md")); err != nil {
		t.Error("README_PARALLEL_VCF_EXPERT_VALIDATION.md not created")
	}
}

// ---------- createTrainingDocumentation ----------

func TestCreateTrainingDocumentation_MissingSourceFile(t *testing.T) {
	tmpDir := t.TempDir()
	_ = createPackageStructure(tmpDir)
	// Source file doesn't exist in test directory → expected error
	err := createTrainingDocumentation(tmpDir)
	if err == nil {
		t.Error("expected error when training documentation source file is missing")
	}
}

// ---------- createParallelExpertMaterials ----------

func TestCreateParallelExpertMaterials_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	// Pre-create the evaluation_materials dir since createParallelExpertMaterials
	// delegates to createParallelEvaluationProtocol which writes directly into it.
	materialsDir := filepath.Join(tmpDir, "evaluation_materials")
	if err := os.MkdirAll(materialsDir, 0750); err != nil {
		t.Fatalf("failed to create evaluation_materials dir: %v", err)
	}
	// createTrainingDocumentation may fail if main.go isn't in the package dir,
	// but the first two (protocol + README) should succeed.
	_ = createParallelExpertMaterials(tmpDir, 100)
	// Verify that the evaluation protocol was created inside the materials dir
	if _, err2 := os.Stat(filepath.Join(materialsDir, "PARALLEL_VCF_EVALUATION_PROTOCOL.md")); err2 != nil {
		t.Errorf("PARALLEL_VCF_EVALUATION_PROTOCOL.md not created: %v", err2)
	}
}

func TestCreateParallelExpertMaterials_ZeroSequences_NoCrash(t *testing.T) {
	tmpDir := t.TempDir()
	// Should not panic with zero sequences
	_ = createParallelExpertMaterials(tmpDir, 0)
}
