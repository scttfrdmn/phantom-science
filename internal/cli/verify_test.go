//go:build !release
// +build !release

package cli

import (
	"bytes"
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// testSequencerConfig returns a small, fingerprint-free config writing to dir.
func testSequencerConfig(dir string) types.SequencerConfig {
	return types.SequencerConfig{
		SampleName:    "verify",
		Lane:          1,
		OutputDir:     dir,
		ReadLength:    150,
		ReadsPerBatch: 200,
		TotalBatches:  1,
		BaseQuality:   38,
		QualityDecay:  0.1,
		DatasetID:     "ds-1",
		Enterprise:    types.EnterpriseConfig{CustomerID: "CUST_1"},
	}
}

// generateFASTQ runs the engine with cfg and returns the files it wrote.
func generateFASTQ(t *testing.T, cfg types.SequencerConfig) []string {
	t.Helper()
	e := genomics.NewEngine(cfg)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatalf("Generate: %v", err)
	}
	files := e.GeneratedFiles()
	if len(files) == 0 {
		t.Fatal("no files generated")
	}
	return files
}

func TestValidateFiles_GeneratedFASTQPasses(t *testing.T) {
	files := generateFASTQ(t, testSequencerConfig(t.TempDir()))
	var out bytes.Buffer
	if err := ValidateFiles(&out, files, DefaultValidationTests); err != nil {
		t.Fatalf("ValidateFiles: %v\n%s", err, out.String())
	}
	if !strings.Contains(out.String(), "All 9 validation tests passed") {
		t.Errorf("unexpected report:\n%s", out.String())
	}
}

func TestValidateFiles_MalformedFails(t *testing.T) {
	path := filepath.Join(t.TempDir(), "bad.fastq")
	if err := os.WriteFile(path, []byte("@r1\nACGT\n+\nII\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if err := ValidateFiles(&out, []string{path}, DefaultValidationTests); err == nil {
		t.Fatalf("malformed FASTQ passed validation:\n%s", out.String())
	}
	if !strings.Contains(out.String(), "❌ fastq_parsing") {
		t.Errorf("report does not flag the length mismatch:\n%s", out.String())
	}
}

// Running zero checks used to report "All validation tests passed".
func TestValidateFiles_NoTestsConfiguredFails(t *testing.T) {
	files := generateFASTQ(t, testSequencerConfig(t.TempDir()))
	var out bytes.Buffer
	err := ValidateFiles(&out, files, types.ValidationConfig{Enable: true})
	if err == nil || !strings.Contains(out.String(), "No validation tests configured") {
		t.Fatalf("err = %v, report:\n%s", err, out.String())
	}
}

func TestValidateAndVerifyRequireFiles(t *testing.T) {
	if err := ValidateFiles(&bytes.Buffer{}, nil, DefaultValidationTests); err == nil {
		t.Error("ValidateFiles with no files succeeded")
	}
	if err := VerifyFiles(&bytes.Buffer{}, nil, VerifyOptions{}); err == nil {
		t.Error("VerifyFiles with no files succeeded")
	}
}

// Regression: Generate validated nil instead of the generated files, so any
// configured test failed with "Invalid data type".
func TestAppGenerate_ValidatesGeneratedFiles(t *testing.T) {
	dir := t.TempDir()
	seqPath := filepath.Join(dir, "seq.json")
	phantomPath := filepath.Join(dir, "phantom.json")
	writeJSON(t, seqPath, testSequencerConfig(filepath.Join(dir, "out")))
	writeJSON(t, phantomPath, types.PhantomConfig{Validation: types.ValidationConfig{
		Enable:                     true,
		StatisticalTests:           []string{"format_compliance", "read_length_consistency"},
		PipelineCompatibilityTests: []string{"fastq_parsing", "quality_encoding"},
	}})

	app, err := NewApp(seqPath, phantomPath, false)
	if err != nil {
		t.Fatalf("NewApp: %v", err)
	}
	if err := app.Generate(); err != nil {
		t.Fatalf("Generate with validation enabled: %v", err)
	}
}

func writeJSON(t *testing.T, path string, v any) {
	t.Helper()
	data, err := json.Marshal(v)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, data, 0o600); err != nil {
		t.Fatal(err)
	}
}

func TestEffectiveValidation(t *testing.T) {
	if _, ok := EffectiveValidation(types.ValidationConfig{}); ok {
		t.Error("disabled validation reported as enabled")
	}
	if got, ok := EffectiveValidation(types.ValidationConfig{Enable: true}); !ok || len(got.StatisticalTests) == 0 {
		t.Errorf("enabled with no tests = %+v, %v; want DefaultValidationTests", got, ok)
	}
	custom := types.ValidationConfig{Enable: true, PipelineCompatibilityTests: []string{"fastq_parsing"}}
	if got, ok := EffectiveValidation(custom); !ok || len(got.StatisticalTests) != 0 || got.PipelineCompatibilityTests[0] != "fastq_parsing" {
		t.Errorf("listed tests not kept: %+v, %v", got, ok)
	}
}

// Gzipped output (output.compression: gzip) must validate like plain FASTQ.
func TestValidateFiles_Gzip(t *testing.T) {
	cfg := testSequencerConfig(t.TempDir())
	cfg.Output.Compression = "gzip"
	files := generateFASTQ(t, cfg)
	if !strings.HasSuffix(files[0], ".fastq.gz") {
		t.Fatalf("expected .fastq.gz output, got %v", files)
	}
	var out bytes.Buffer
	if err := ValidateFiles(&out, files, DefaultValidationTests); err != nil {
		t.Fatalf("ValidateFiles on gzip output: %v\n%s", err, out.String())
	}

	bad := filepath.Join(t.TempDir(), "bad.fastq.gz")
	if err := os.WriteFile(bad, []byte("not gzip"), 0o600); err != nil {
		t.Fatal(err)
	}
	if err := ValidateFiles(&bytes.Buffer{}, []string{bad}, DefaultValidationTests); err == nil {
		t.Error("corrupt .gz file passed validation")
	}
}
