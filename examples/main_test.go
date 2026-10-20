//go:build !release

package main

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// writeKey writes a key file with mode 0600 and returns its path.
func writeKey(t *testing.T, dir, key string) string {
	t.Helper()
	path := filepath.Join(dir, "fingerprint.key")
	if err := os.WriteFile(path, []byte(key+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	return path
}

func fastqFiles(t *testing.T, dir string) []string {
	t.Helper()
	files, err := filepath.Glob(filepath.Join(dir, "*.fastq"))
	if err != nil {
		t.Fatal(err)
	}
	return files
}

// The example's round trip: every generated file is attributed to the
// customer and dataset by regeneration with the key, and only with it.
// No stegano library is needed.
func TestFingerprintDataset_RoundTrip(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "out")
	keyFile := writeKey(t, dir, "example-key-0123456789abcdef")
	t.Setenv(fingerprinting.KeyEnvVar, "")

	files, man, err := FingerprintDataset(Options{OutputDir: out, CustomerID: "ACME", DatasetID: "study_1", KeyFile: keyFile, Reads: 500, Seed: 1})
	if err != nil {
		t.Fatalf("FingerprintDataset: %v", err)
	}
	if len(files) != 2 {
		t.Fatalf("got %d files %v, want an R1/R2 pair", len(files), files)
	}

	var report bytes.Buffer
	if err := VerifyDataset(&report, files, keyFile, man); err != nil {
		t.Fatalf("VerifyDataset: %v\n%s", err, report.String())
	}
	if n := strings.Count(report.String(), "✅"); n != 2 || !strings.Contains(report.String(), "customer=ACME dataset=study_1 salt=1") {
		t.Errorf("report has %d matches, want 2 attributed to ACME/study_1:\n%s", n, report.String())
	}

	// A wrong key and no key both fail.
	wrongKey := writeKey(t, t.TempDir(), "another-key-0123456789abcdef")
	if err := VerifyDataset(&bytes.Buffer{}, files, wrongKey, man); err == nil {
		t.Error("VerifyDataset succeeded with the wrong key")
	}
	if err := VerifyDataset(&bytes.Buffer{}, files, "", man); err == nil {
		t.Error("VerifyDataset succeeded without a key")
	}
}

// The fingerprinted files are exactly the files of an unfingerprinted run
// with the derived seed: fingerprinting changes nothing in the data.
func TestFingerprintDataset_DataUnchanged(t *testing.T) {
	dir := t.TempDir()
	keyFile := writeKey(t, dir, "example-key-0123456789abcdef")
	files, _, err := FingerprintDataset(Options{OutputDir: filepath.Join(dir, "fp"), CustomerID: "ACME", DatasetID: "study_1", KeyFile: keyFile, Reads: 200, Seed: 9})
	if err != nil {
		t.Fatal(err)
	}
	seed, err := fingerprinting.DeriveSeed([]byte("example-key-0123456789abcdef"), "ACME", "study_1", 9)
	if err != nil {
		t.Fatal(err)
	}
	plainDir := filepath.Join(dir, "plain")
	e := genomics.NewEngineLog(types.SequencerConfig{
		SampleName: "fingerprint-example", Lane: 1, OutputDir: plainDir, ReadLength: 150,
		ReadsPerBatch: 200, TotalBatches: 1, PairedEnd: true, BaseQuality: 35, QualityDecay: 0.1, Seed: seed,
	}, nil)
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	for _, f := range files {
		got, _ := os.ReadFile(f)
		want, err := os.ReadFile(filepath.Join(plainDir, filepath.Base(f)))
		if err != nil || !bytes.Equal(got, want) {
			t.Errorf("%s differs from the unfingerprinted run with the derived seed (%v)", filepath.Base(f), err)
		}
	}
}

// Without a key, generation refuses to start and writes no FASTQ files.
func TestFingerprintDataset_RequiresKey(t *testing.T) {
	out := t.TempDir()
	t.Setenv(fingerprinting.KeyEnvVar, "")
	if _, _, err := FingerprintDataset(Options{OutputDir: out, CustomerID: "ACME", DatasetID: "d", Reads: 100, Seed: 1}); err == nil {
		t.Fatal("FingerprintDataset succeeded without a key")
	}
	if files := fastqFiles(t, out); len(files) != 0 {
		t.Errorf("files written without a key: %v", files)
	}
}

// A key file readable by other users is rejected.
func TestFingerprintDataset_RejectsExposedKeyFile(t *testing.T) {
	dir := t.TempDir()
	keyFile := writeKey(t, dir, "example-key-0123456789abcdef")
	if err := os.Chmod(keyFile, 0o644); err != nil {
		t.Fatal(err)
	}
	_, _, err := FingerprintDataset(Options{OutputDir: filepath.Join(dir, "out"), CustomerID: "ACME", DatasetID: "d", KeyFile: keyFile, Reads: 100, Seed: 1})
	if err == nil || !strings.Contains(err.Error(), "chmod 600") {
		t.Fatalf("err = %v, want advice to chmod 600 the key file", err)
	}
}

// Reads that the run did not produce are not attributed to it.
func TestVerifyDataset_UnrelatedFile(t *testing.T) {
	dir := t.TempDir()
	keyFile := writeKey(t, dir, "example-key-0123456789abcdef")
	_, man, err := FingerprintDataset(Options{OutputDir: filepath.Join(dir, "out"), CustomerID: "ACME", DatasetID: "d", KeyFile: keyFile, Reads: 100, Seed: 1})
	if err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(dir, "other.fastq")
	if err := os.WriteFile(path, []byte("@r1\nACGTTGCAAGGCTTAACCGGTTAAGCTAGCTAGGATCCAGT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	var report bytes.Buffer
	if err := VerifyDataset(&report, []string{path}, keyFile, man); err == nil {
		t.Fatalf("VerifyDataset attributed an unrelated file:\n%s", report.String())
	}
}
