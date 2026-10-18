//go:build !release && cgo
// +build !release,cgo

package genomics

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func qualityFingerprintConfig(dir string, reads, batches int, paired bool) types.SequencerConfig {
	return types.SequencerConfig{
		SampleName:    "fp-sample",
		Lane:          1,
		OutputDir:     dir,
		ReadLength:    100,
		ReadsPerBatch: reads,
		TotalBatches:  batches,
		PairedEnd:     paired,
		BaseQuality:   38,
		QualityDecay:  0.1,
		DatasetID:     "dataset-42",
		Enterprise:    types.EnterpriseConfig{CustomerID: "customer-7"},
		Fingerprinting: types.FingerprintConfig{
			Enable:              true,
			SteganographyLayers: []string{"quality_scores"},
		},
	}
}

func extractFileFingerprint(t *testing.T, e *Engine, path string) (string, error) {
	t.Helper()
	ff, ok := e.fingerprintingEngine.(fingerprinting.FileFingerprinter)
	if !ok {
		t.Fatalf("engine %T cannot extract file fingerprints", e.fingerprintingEngine)
	}
	return ff.ExtractFromFASTQFile(path)
}

// Regression for the per-read embedding bug: every generated file must carry
// one extractable fingerprint, the same for the whole run.
func TestGenerate_EmbedsOneFingerprintPerRun(t *testing.T) {
	dir := t.TempDir()
	cfg := qualityFingerprintConfig(dir, 50, 2, true)
	e := NewEngine(cfg)
	defer e.Close()

	if err := e.Generate(); err != nil {
		t.Fatalf("Generate: %v", err)
	}

	files, err := filepath.Glob(filepath.Join(dir, "*.fastq"))
	if err != nil || len(files) != 4 {
		t.Fatalf("want 4 FASTQ files (2 batches × R1/R2), got %v (err %v)", files, err)
	}

	var want string
	for _, f := range files {
		got, err := extractFileFingerprint(t, e, f)
		if err != nil {
			t.Fatalf("%s carries no fingerprint: %v", filepath.Base(f), err)
		}
		if want == "" {
			want = got
			if !strings.HasPrefix(want, "customer-7:dataset-42:stegano-") {
				t.Errorf("fingerprint %q does not identify customer and dataset", want)
			}
		} else if got != want {
			t.Errorf("%s fingerprint %q differs from %q", filepath.Base(f), got, want)
		}

		info, err := os.Stat(f)
		if err != nil {
			t.Fatal(err)
		}
		if perm := info.Mode().Perm(); perm != 0o600 {
			t.Errorf("%s mode %o, want 0600 preserved through embedding", filepath.Base(f), perm)
		}

		data, err := os.ReadFile(f)
		if err != nil {
			t.Fatal(err)
		}
		lines := strings.Split(strings.TrimSuffix(string(data), "\n"), "\n")
		if len(lines) != 4*cfg.ReadsPerBatch {
			t.Fatalf("%s has %d lines, want %d", filepath.Base(f), len(lines), 4*cfg.ReadsPerBatch)
		}
		for i := 0; i < len(lines); i += 4 {
			if len(lines[i+1]) != cfg.ReadLength || len(lines[i+3]) != cfg.ReadLength {
				t.Fatalf("%s record %d: sequence/quality length %d/%d, want %d",
					filepath.Base(f), i/4, len(lines[i+1]), len(lines[i+3]), cfg.ReadLength)
			}
		}
	}
}

// GenerateFASTQFile (the TUI path) fingerprints each file it writes.
func TestGenerateFASTQFile_Fingerprinted(t *testing.T) {
	dir := t.TempDir()
	e := NewEngine(qualityFingerprintConfig(dir, 50, 1, false))
	defer e.Close()

	path := filepath.Join(dir, "tui.fastq")
	if err := e.GenerateFASTQFile(path, 1, 1); err != nil {
		t.Fatalf("GenerateFASTQFile: %v", err)
	}
	if _, err := extractFileFingerprint(t, e, path); err != nil {
		t.Errorf("file carries no fingerprint: %v", err)
	}
}

// Without the quality_scores layer no fingerprint is embedded, even with
// fingerprinting enabled (read-ID tagging only).
func TestGenerate_ReadIDLayerOnly_NoQualityFingerprint(t *testing.T) {
	dir := t.TempDir()
	cfg := qualityFingerprintConfig(dir, 50, 1, false)
	cfg.Fingerprinting.SteganographyLayers = []string{"read_ids"}
	e := NewEngine(cfg)
	defer e.Close()

	if err := e.Generate(); err != nil {
		t.Fatalf("Generate: %v", err)
	}
	files, _ := filepath.Glob(filepath.Join(dir, "*.fastq"))
	if len(files) != 1 {
		t.Fatalf("want 1 file, got %v", files)
	}
	if fp, err := extractFileFingerprint(t, e, files[0]); err == nil {
		t.Errorf("unexpected fingerprint %q without the quality_scores layer", fp)
	}
}

// A file too small to hold the fingerprint fails generation and is removed,
// instead of being written unmarked.
func TestGenerate_TooSmallForFingerprint(t *testing.T) {
	dir := t.TempDir()
	e := NewEngine(qualityFingerprintConfig(dir, 1, 1, false))
	defer e.Close()

	err := e.Generate()
	if !errors.Is(err, fingerprinting.ErrCapacityExceeded) {
		t.Fatalf("Generate error = %v, want ErrCapacityExceeded", err)
	}
	if files, _ := filepath.Glob(filepath.Join(dir, "*.fastq")); len(files) != 0 {
		t.Errorf("unmarked output left behind: %v", files)
	}
}
