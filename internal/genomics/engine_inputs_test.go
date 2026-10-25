//go:build !release

package genomics

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
)

// Regenerate checks the reference a manifest records (#37): it refuses an
// engine without a reference or with another FASTA, rather than reproducing
// different reads, and reproduces the run with the recorded one.
func TestRegenerate_ChecksReference(t *testing.T) {
	ref := writeTestReference(t)
	c := refConfig(t, ref)
	c.TotalBatches = 1
	_, e := runTree(t, c)
	m := e.Manifest()
	if m.Reference == nil || m.Reference.Name != "ref.fa" || m.Config.ReferencePath != "" {
		t.Fatalf("manifest reference %+v, config reference_path %q", m.Reference, m.Config.ReferencePath)
	}
	other := filepath.Join(t.TempDir(), "other.fa")
	data, err := os.ReadFile(ref.path)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(other, append(data, []byte(">chr3\nACGTACGTAC\n")...), 0o600); err != nil {
		t.Fatal(err)
	}
	files := []RunFile{{FlowCell: 1, Batch: 1, Mate: 1, Reads: int(m.Outputs[0].Reads)}}
	regen := func(path string) error {
		cfg := m.Config
		cfg.Seed, cfg.ReferencePath = m.Seed, path
		e := NewEngineLog(cfg, nil)
		defer e.Close()
		_, err := e.Regenerate(RegenOptions{Files: files, Reference: m.Reference}, func(RunFile, []byte) error { return nil })
		return err
	}
	for path, want := range map[string]string{"": "configure that FASTA", other: "not other.fa"} {
		if err := regen(path); err == nil || !strings.Contains(err.Error(), want) || !strings.Contains(err.Error(), m.Reference.SHA256) {
			t.Errorf("reference %q: err = %v, want %q", path, err, want)
		}
	}
	if err := regen(ref.path); err != nil {
		t.Errorf("with the recorded reference: %v", err)
	}
	if got, err := manifest.HashInput(ref.path); err != nil || !got.Same(m.Reference) {
		t.Errorf("HashInput = %+v, %v", got, err)
	}
}
