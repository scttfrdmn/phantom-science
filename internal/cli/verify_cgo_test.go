//go:build !release && cgo
// +build !release,cgo

package cli

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func writeTestKey(t *testing.T, dir, key string) string {
	t.Helper()
	path := filepath.Join(dir, filepath.Base(dir)+".key")
	if err := os.WriteFile(path, []byte(key+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestVerifyFiles(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	dir := t.TempDir()
	keyFile := writeTestKey(t, dir, "verify-test-key-0123456789")
	wrongKeyFile := writeTestKey(t, t.TempDir(), "a-different-key-0123456789")

	cfg := testSequencerConfig(filepath.Join(dir, "enc"))
	cfg.Fingerprinting = types.FingerprintConfig{Enable: true, SteganographyLayers: []string{"quality_scores"}, KeyFile: keyFile}
	encrypted := generateFASTQ(t, cfg)[0]

	plainCfg := testSequencerConfig(filepath.Join(dir, "plaintext"))
	plainCfg.Fingerprinting = types.FingerprintConfig{Enable: true, SteganographyLayers: []string{"quality_scores"}, Plaintext: true}
	plaintext := generateFASTQ(t, plainCfg)[0]

	unmarked := generateFASTQ(t, testSequencerConfig(filepath.Join(dir, "unmarked")))[0]

	var out bytes.Buffer
	if err := VerifyFiles(&out, []string{encrypted, plaintext}, VerifyOptions{Customer: "CUST_1", Dataset: "ds-1", KeyFile: keyFile}); err != nil {
		t.Fatalf("VerifyFiles: %v\n%s", err, out.String())
	}
	for _, want := range []string{"customer=CUST_1 dataset=ds-1", "(encrypted)", "(plaintext)"} {
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
		"unmarked file":     {[]string{unmarked}, VerifyOptions{KeyFile: keyFile}, "no fingerprint"},
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

// verify must read fingerprints from gzipped output, which is fingerprinted
// before compression.
func TestVerifyFiles_Gzip(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	dir := t.TempDir()
	keyFile := writeTestKey(t, dir, "verify-gzip-key-0123456789")
	cfg := testSequencerConfig(filepath.Join(dir, "ont"))
	cfg.Platform = "ont"
	cfg.ReadLength = 2000
	cfg.Output.Compression = "gzip"
	cfg.Fingerprinting = types.FingerprintConfig{Enable: true, SteganographyLayers: []string{"quality_scores"}, KeyFile: keyFile}
	files := generateFASTQ(t, cfg)
	if !strings.HasSuffix(files[0], ".fastq.gz") {
		t.Fatalf("expected .fastq.gz output, got %v", files)
	}
	var out bytes.Buffer
	if err := VerifyFiles(&out, files, VerifyOptions{Customer: "CUST_1", Dataset: "ds-1", KeyFile: keyFile}); err != nil {
		t.Fatalf("VerifyFiles on gzip output: %v\n%s", err, out.String())
	}
	if !strings.Contains(out.String(), "(encrypted)") {
		t.Errorf("unexpected report:\n%s", out.String())
	}
}

