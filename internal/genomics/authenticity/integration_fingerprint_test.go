//go:build !release

package authenticity

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// fakeFingerprinter is a fingerprinting engine without file-level embedding;
// fingerprint generation and per-buffer embedding can be made to fail.
type fakeFingerprinter struct {
	genErr   error
	embedErr error
}

func (f *fakeFingerprinter) GenerateFingerprint(string) (string, error) { return "FP", f.genErr }
func (f *fakeFingerprinter) EmbedInQualities(q []byte, _ string) ([]byte, error) {
	return q, f.embedErr
}
func (f *fakeFingerprinter) EmbedCustomerID(id string) string            { return id }
func (f *fakeFingerprinter) EmbedDatasetFingerprint(id, _ string) string { return id }
func (f *fakeFingerprinter) VerifyFingerprint([]byte) bool               { return false }
func (f *fakeFingerprinter) ExtractFingerprint([]byte) (string, error)   { return "", nil }
func (f *fakeFingerprinter) Close()                                      {}

// fakeFileFingerprinter adds whole-file embedding, appending a marker line
// to the file on success.
type fakeFileFingerprinter struct{ fakeFingerprinter }

func (f *fakeFileFingerprinter) EmbedInFASTQFile(in, out, fp string) error {
	if f.embedErr != nil {
		return f.embedErr
	}
	data, err := os.ReadFile(in)
	if err != nil {
		return err
	}
	return os.WriteFile(out, append(data, []byte("#marked:"+fp+"\n")...), 0o600)
}
func (f *fakeFileFingerprinter) ExtractFromFASTQFile(string) (string, error) { return "FP", nil }

func fingerprintedIntegrationManager(t *testing.T, fp fingerprinting.Engine) (*IntegrationManager, string) {
	t.Helper()
	engine := createInitializedEngine(t)
	t.Cleanup(func() { engine.Close() })
	dir := t.TempDir()
	cfg := types.SequencerConfig{
		SampleName: "S", Lane: 1, ReadLength: 40, BaseQuality: 30, QualityDecay: 0.1,
		ReadsPerBatch: 5, OutputDir: dir, DatasetID: "DS1",
		Fingerprinting: types.FingerprintConfig{
			Enable:              true,
			SteganographyLayers: []string{"quality_scores"},
		},
	}
	im := &IntegrationManager{authenticityEngine: engine, fingerprintingEngine: fp, config: cfg}
	return im, filepath.Join(dir, "out.fastq")
}

// With quality-score fingerprinting enabled, a fingerprint that cannot be
// generated or embedded fails the call and leaves no unmarked file behind.
func TestGenerateAuthenticFASTQFile_FingerprintFailureFailsClosed(t *testing.T) {
	embedErr := errors.New("embed failed")
	cases := map[string]fingerprinting.Engine{
		"embed error":            &fakeFileFingerprinter{fakeFingerprinter{embedErr: embedErr}},
		"generate error":         &fakeFileFingerprinter{fakeFingerprinter{genErr: embedErr}},
		"no file-level embedder": &fakeFingerprinter{},
	}
	for name, fp := range cases {
		t.Run(name, func(t *testing.T) {
			im, path := fingerprintedIntegrationManager(t, fp)
			if err := im.GenerateAuthenticFASTQFile(path, 1, 1); err == nil {
				t.Fatal("GenerateAuthenticFASTQFile = nil, want fingerprinting error")
			}
			if _, err := os.Stat(path); !os.IsNotExist(err) {
				t.Errorf("unmarked file left at %s (stat err %v)", path, err)
			}
		})
	}
}

func TestGenerateAuthenticFASTQFile_FingerprintEmbedded(t *testing.T) {
	im, path := fingerprintedIntegrationManager(t, &fakeFileFingerprinter{})
	if err := im.GenerateAuthenticFASTQFile(path, 1, 1); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.HasSuffix(string(data), "#marked:FP\n") {
		t.Error("file was not passed through the file-level fingerprint embedder")
	}
	if got := strings.Count(string(data), "\n+\n"); got != 5 {
		t.Errorf("%d FASTQ records, want 5", got)
	}
}
