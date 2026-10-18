//go:build !release

package main

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
)

// steganoAvailable reports whether this build links the stegano library.
func steganoAvailable() bool {
	_, err := fingerprinting.NewVerifier()
	return err == nil
}

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

// The example's round trip: every generated file carries an encrypted
// fingerprint that names the customer and dataset and is readable only with
// the key it was embedded with.
func TestFingerprintDataset_EncryptedRoundTrip(t *testing.T) {
	if !steganoAvailable() {
		t.Skip("needs the stegano library (CGO build)")
	}
	dir := t.TempDir()
	out := filepath.Join(dir, "out")
	keyFile := writeKey(t, dir, "example-key-0123456789abcdef")
	t.Setenv(fingerprinting.KeyEnvVar, "")

	files, err := FingerprintDataset(Options{OutputDir: out, CustomerID: "ACME", DatasetID: "study_1", KeyFile: keyFile, Reads: 500, Seed: 1})
	if err != nil {
		t.Fatalf("FingerprintDataset: %v", err)
	}
	if len(files) != 2 {
		t.Fatalf("got %d files %v, want an R1/R2 pair", len(files), files)
	}

	var report bytes.Buffer
	if err := VerifyDataset(&report, files, keyFile, "ACME", "study_1"); err != nil {
		t.Fatalf("VerifyDataset: %v\n%s", err, report.String())
	}
	if n := strings.Count(report.String(), "(encrypted)"); n != 2 {
		t.Errorf("report has %d encrypted fingerprints, want 2:\n%s", n, report.String())
	}

	// The wrong customer, a wrong key and no key all fail.
	if err := VerifyDataset(&bytes.Buffer{}, files, keyFile, "OTHER", "study_1"); err == nil {
		t.Error("VerifyDataset accepted a fingerprint for another customer")
	}
	wrongKey := writeKey(t, t.TempDir(), "another-key-0123456789abcdef")
	if err := VerifyDataset(&bytes.Buffer{}, files, wrongKey, "ACME", "study_1"); err == nil {
		t.Error("VerifyDataset succeeded with the wrong key")
	}
	if err := VerifyDataset(&bytes.Buffer{}, files, "", "ACME", "study_1"); err == nil {
		t.Error("VerifyDataset succeeded without a key")
	}
}

// Without a key, generation refuses to embed and writes no FASTQ files.
func TestFingerprintDataset_RequiresKey(t *testing.T) {
	out := t.TempDir()
	t.Setenv(fingerprinting.KeyEnvVar, "")
	if _, err := FingerprintDataset(Options{OutputDir: out, CustomerID: "ACME", DatasetID: "d", Reads: 100, Seed: 1}); err == nil {
		t.Fatal("FingerprintDataset succeeded without a key")
	}
	if files := fastqFiles(t, out); len(files) != 0 {
		t.Errorf("unmarked files left behind: %v", files)
	}
}

// In builds without stegano, generation fails before writing anything
// instead of writing files without the requested fingerprint.
func TestFingerprintDataset_NoSteganoFailsClosed(t *testing.T) {
	if steganoAvailable() {
		t.Skip("stegano library is linked")
	}
	dir := t.TempDir()
	out := filepath.Join(dir, "out")
	keyFile := writeKey(t, dir, "example-key-0123456789abcdef")
	_, err := FingerprintDataset(Options{OutputDir: out, CustomerID: "ACME", DatasetID: "d", KeyFile: keyFile, Reads: 100, Seed: 1})
	if err == nil || !strings.Contains(err.Error(), "stegano") {
		t.Fatalf("err = %v, want an error naming the stegano library", err)
	}
	if files := fastqFiles(t, out); len(files) != 0 {
		t.Errorf("unmarked files written: %v", files)
	}
}

// A key file readable by other users is rejected.
func TestFingerprintDataset_RejectsExposedKeyFile(t *testing.T) {
	dir := t.TempDir()
	keyFile := writeKey(t, dir, "example-key-0123456789abcdef")
	if err := os.Chmod(keyFile, 0o644); err != nil {
		t.Fatal(err)
	}
	_, err := FingerprintDataset(Options{OutputDir: filepath.Join(dir, "out"), CustomerID: "ACME", DatasetID: "d", KeyFile: keyFile, Reads: 100, Seed: 1})
	if err == nil {
		t.Fatal("FingerprintDataset accepted a key file with mode 0644")
	}
	if steganoAvailable() && !strings.Contains(err.Error(), "chmod 600") {
		t.Errorf("err = %v, want advice to chmod 600 the key file", err)
	}
}

func TestVerifyDataset_UnfingerprintedFile(t *testing.T) {
	if !steganoAvailable() {
		t.Skip("needs the stegano library (CGO build)")
	}
	dir := t.TempDir()
	path := filepath.Join(dir, "plain.fastq")
	if err := os.WriteFile(path, []byte("@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	var report bytes.Buffer
	if err := VerifyDataset(&report, []string{path}, writeKey(t, dir, "example-key-0123456789abcdef"), "ACME", "d"); err == nil {
		t.Fatalf("VerifyDataset accepted a file without a fingerprint:\n%s", report.String())
	}
}
