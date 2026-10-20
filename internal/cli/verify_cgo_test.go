//go:build !release && cgo
// +build !release,cgo

package cli

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
)

// Generation no longer writes quality-score fingerprints (they changed
// quality values; docs/CHARTER.md, principle 2), but verify still reads them
// from files made by phantom-science 0.9.0 and earlier. These tests make
// such files the way 0.9.0 did: generate, then embed one fingerprint per
// file with stegano's file-level embedding.

func writeTestKey(t *testing.T, dir, key string) string {
	t.Helper()
	path := filepath.Join(dir, filepath.Base(dir)+".key")
	if err := os.WriteFile(path, []byte(key+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	return path
}

// legacyFingerprint embeds a 0.9.0-style fingerprint for CUST_1/ds-1 in the
// FASTQ file at path, encrypted with key (plaintext when key is empty).
func legacyFingerprint(t *testing.T, path, key string) {
	t.Helper()
	var opts []fingerprinting.Option
	if key != "" {
		opts = append(opts, fingerprinting.WithKey([]byte(key)))
	}
	e, err := fingerprinting.NewSteganoEngine("CUST_1", opts...)
	if err != nil {
		t.Fatal(err)
	}
	fp := fingerprinting.FormatFingerprint("CUST_1", "ds-1", fingerprinting.LibraryVersion(), time.Unix(1700000000, 0))
	if err := e.EmbedInFASTQFile(path, path, fp); err != nil {
		t.Fatalf("embedding a legacy fingerprint: %v", err)
	}
}

func TestVerifyFiles(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	dir := t.TempDir()
	keyFile := writeTestKey(t, dir, "verify-test-key-0123456789")
	wrongKeyFile := writeTestKey(t, t.TempDir(), "a-different-key-0123456789")

	encrypted := generateFASTQ(t, testSequencerConfig(filepath.Join(dir, "enc")))[0]
	legacyFingerprint(t, encrypted, "verify-test-key-0123456789")
	plaintext := generateFASTQ(t, testSequencerConfig(filepath.Join(dir, "plaintext")))[0]
	legacyFingerprint(t, plaintext, "")

	unmarked := generateFASTQ(t, testSequencerConfig(filepath.Join(dir, "unmarked")))[0]

	var out bytes.Buffer
	if err := VerifyFiles(&out, []string{encrypted, plaintext}, VerifyOptions{Customer: "CUST_1", Dataset: "ds-1", KeyFile: keyFile}); err != nil {
		t.Fatalf("VerifyFiles: %v\n%s", err, out.String())
	}
	for _, want := range []string{"legacy quality-score fingerprint", "customer=CUST_1 dataset=ds-1", "(encrypted)", "(plaintext)"} {
		if !strings.Contains(out.String(), want) {
			t.Errorf("report lacks %q:\n%s", want, out.String())
		}
	}

	for name, tc := range map[string]struct {
		paths []string
		opts  VerifyOptions
		want  string
	}{
		"no key":            {[]string{encrypted}, VerifyOptions{}, "fingerprint is encrypted"},
		"wrong key":         {[]string{encrypted}, VerifyOptions{KeyFile: wrongKeyFile}, "cannot be decrypted"},
		"customer mismatch": {[]string{encrypted}, VerifyOptions{KeyFile: keyFile, Customer: "OTHER"}, `expected "OTHER"`},
		"unmarked file":     {[]string{unmarked}, VerifyOptions{KeyFile: keyFile}, "no quality-score fingerprint"},
		"missing file":      {[]string{filepath.Join(dir, "missing.fastq")}, VerifyOptions{}, "cannot verify"},
	} {
		t.Run(name, func(t *testing.T) {
			var out bytes.Buffer
			if err := VerifyFiles(&out, tc.paths, tc.opts); err == nil {
				t.Fatalf("VerifyFiles succeeded:\n%s", out.String())
			}
			if !strings.Contains(out.String(), tc.want) {
				t.Errorf("report lacks %q:\n%s", tc.want, out.String())
			}
		})
	}
}

// verify must read legacy fingerprints from gzipped files, which 0.9.0
// fingerprinted before compression.
func TestVerifyFiles_Gzip(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	dir := t.TempDir()
	keyFile := writeTestKey(t, dir, "verify-gzip-key-0123456789")
	cfg := testSequencerConfig(filepath.Join(dir, "ont"))
	cfg.Platform = "ont"
	cfg.ReadLength = 2000
	plain := generateFASTQ(t, cfg)[0]
	legacyFingerprint(t, plain, "verify-gzip-key-0123456789")
	data, err := os.ReadFile(plain)
	if err != nil {
		t.Fatal(err)
	}
	var gz bytes.Buffer
	zw := gzip.NewWriter(&gz)
	if _, err := zw.Write(data); err != nil {
		t.Fatal(err)
	}
	if err := zw.Close(); err != nil {
		t.Fatal(err)
	}
	gzPath := plain + ".gz"
	if err := os.WriteFile(gzPath, gz.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if err := VerifyFiles(&out, []string{gzPath}, VerifyOptions{Customer: "CUST_1", Dataset: "ds-1", KeyFile: keyFile}); err != nil {
		t.Fatalf("VerifyFiles on gzip output: %v\n%s", err, out.String())
	}
	if !strings.Contains(out.String(), "(encrypted)") {
		t.Errorf("unexpected report:\n%s", out.String())
	}
}
