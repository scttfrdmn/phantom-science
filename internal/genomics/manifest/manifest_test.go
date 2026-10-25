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
		VariantsPath: "/home/me/calls.vcf.gz", VariantsRegion: "chr5:1-900",
		ReferencePath: "/home/me/hg38.fa", ReferenceRegion: "chr22",
		HMMModelPath: "/home/me/m.hmm.json", ModelDir: "/home/me/models", Population: "YRI",
		ReadLengthDistribution: types.ReadLengthSpec{Dist: types.ReadLengthEmpirical, EmpiricalPath: "/home/me/lengths.tsv"},
		Fingerprinting:         types.FingerprintConfig{Enable: true, KeyFile: "/secret/k", SteganographyLayers: []string{"quality_scores"}},
		Enterprise: types.EnterpriseConfig{CustomerID: "c", LicenseKey: "LIC-123",
			Fingerprinting: types.FingerprintConfig{KeyFile: "/secret/k2"}},
	}
	s := SanitizeConfig(c)
	if s.OutputDir != "" || s.Workers != 0 || s.Fingerprinting.KeyFile != "" || s.Enterprise.Fingerprinting.KeyFile != "" || s.Enterprise.LicenseKey != "" || s.VariantsPath != "" {
		t.Errorf("not sanitized: %+v", s)
	}
	// No local path is left (#37): the manifest records input files by hash.
	m := New("", 1, time.Unix(0, 0), c)
	data, err := m.Marshal()
	if err != nil {
		t.Fatal(err)
	}
	for _, leak := range []string{"/home/me", "/secret", "LIC-123"} {
		if bytes.Contains(data, []byte(leak)) {
			t.Errorf("manifest holds %s:\n%s", leak, data)
		}
	}
	if s.SampleName != "S" || s.Seed != 3 || !s.Fingerprinting.Enable || s.Enterprise.CustomerID != "c" || s.VariantsRegion != "chr5:1-900" ||
		s.ReferenceRegion != "chr22" || s.Population != "YRI" || s.ReadLengthDistribution.Dist != types.ReadLengthEmpirical {
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

// A signal run has no SequencerConfig: the config is omitted, the signal
// section is recorded, and BLOW5 outputs count as read outputs.
func TestManifest_SignalRun(t *testing.T) {
	m := New("v", 3, time.Date(2024, 1, 1, 0, 0, 0, 0, time.UTC), types.SequencerConfig{})
	m.Signal = &Signal{Profile: "p", Inputs: []Input{{Name: "a.fastq", Bytes: 1, SHA256: "x"}}}
	m.Add(Output{Path: "a_0.blow5", Role: RoleSignal, Bytes: 10, Reads: 4, Bases: 400})
	m.Add(Output{Path: "truth.fa", Role: RoleTruth, Bytes: 5})
	if m.Totals != (Totals{Files: 1, Reads: 4, Bases: 400, Bytes: 15}) {
		t.Errorf("totals %+v", m.Totals)
	}
	data, err := m.Marshal()
	if err != nil {
		t.Fatal(err)
	}
	if bytes.Contains(data, []byte(`"config"`)) || !bytes.Contains(data, []byte(`"signal": {`)) {
		t.Errorf("signal manifest JSON:\n%s", data)
	}
	g := New("v", 3, time.Time{}, types.SequencerConfig{Platform: "ont"})
	if data, _ := g.Marshal(); !bytes.Contains(data, []byte(`"config"`)) || bytes.Contains(data, []byte(`"signal"`)) {
		t.Errorf("generation manifest JSON:\n%s", data)
	}
}
