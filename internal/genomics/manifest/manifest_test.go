//go:build !release

package manifest

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestSanitizeConfig(t *testing.T) {
	c := types.SequencerConfig{
		SampleName: "S", OutputDir: "/home/me/out", Workers: 8, Seed: 3,
		Fingerprinting: types.FingerprintConfig{Enable: true, KeyFile: "/secret/k", SteganographyLayers: []string{"quality_scores"}},
		Enterprise: types.EnterpriseConfig{CustomerID: "c", LicenseKey: "LIC-123",
			Fingerprinting: types.FingerprintConfig{KeyFile: "/secret/k2"}},
	}
	s := SanitizeConfig(c)
	if s.OutputDir != "" || s.Workers != 0 || s.Fingerprinting.KeyFile != "" || s.Enterprise.Fingerprinting.KeyFile != "" || s.Enterprise.LicenseKey != "" {
		t.Errorf("not sanitized: %+v", s)
	}
	if s.SampleName != "S" || s.Seed != 3 || !s.Fingerprinting.Enable || s.Enterprise.CustomerID != "c" {
		t.Errorf("settings lost: %+v", s)
	}
	s.Fingerprinting.SteganographyLayers[0] = "x"
	if c.Fingerprinting.SteganographyLayers[0] != "quality_scores" {
		t.Error("SanitizeConfig aliases the caller's slices")
	}
}

func TestManifest_RoundTripAndVerify(t *testing.T) {
	dir := t.TempDir()
	if err := os.MkdirAll(filepath.Join(dir, "FC1"), 0o750); err != nil {
		t.Fatal(err)
	}
	files := map[string]string{"FC1/a.fastq": "@r\nAC\n+\nII\n", "FC1/sequencing_summary.txt": "h\n"}
	m := New("", 7, time.Date(2024, 5, 1, 2, 3, 4, 0, time.FixedZone("x", 3600)), types.SequencerConfig{OutputDir: dir})
	for rel, content := range files {
		p := filepath.Join(dir, filepath.FromSlash(rel))
		if err := os.WriteFile(p, []byte(content), 0o600); err != nil {
			t.Fatal(err)
		}
		n, sum, err := HashFile(p)
		if err != nil {
			t.Fatal(err)
		}
		role := RoleFASTQ
		if strings.HasSuffix(rel, ".txt") {
			role = RoleSummary
		}
		m.Add(Output{Path: rel, Role: role, FlowCell: 1, Bytes: n, SHA256: sum, Reads: 1, Bases: 2})
	}
	if m.Totals.Files != 1 || m.Totals.Reads != 1 || m.Totals.Bases != 2 || m.Totals.Bytes != int64(len(files["FC1/a.fastq"])+2) {
		t.Errorf("totals %+v", m.Totals)
	}
	if m.Tool.Version != "dev" || m.RunStart.Location() != time.UTC {
		t.Errorf("header %+v", m)
	}
	data, err := m.Marshal()
	if err != nil {
		t.Fatal(err)
	}
	if bytes.Contains(data, []byte(dir)) {
		t.Error("manifest records the output directory")
	}
	again, _ := m.Marshal()
	if !bytes.Equal(data, again) {
		t.Error("Marshal is not deterministic")
	}
	path := filepath.Join(dir, FileName)
	if err := os.WriteFile(path, data, 0o600); err != nil {
		t.Fatal(err)
	}
	got, err := Read(path)
	if err != nil {
		t.Fatal(err)
	}
	if err := got.VerifyDir(dir); err != nil {
		t.Errorf("VerifyDir on intact output: %v", err)
	}
	if err := os.WriteFile(filepath.Join(dir, "FC1", "a.fastq"), []byte("@r\nAG\n+\nII\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	if err := got.VerifyDir(dir); err == nil || !strings.Contains(err.Error(), "FC1/a.fastq") {
		t.Errorf("VerifyDir missed a modified file: %v", err)
	}
	if err := os.WriteFile(path, []byte(`{"format":"other"}`), 0o600); err != nil {
		t.Fatal(err)
	}
	if _, err := Read(path); err == nil {
		t.Error("Read accepted a foreign format")
	}
}
