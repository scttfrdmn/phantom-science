//go:build !release

package cmd

import (
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/validation"
)

// repoRoot returns the absolute path of the repository root (the config
// loader rejects relative paths containing "..").
func repoRoot(t *testing.T) string {
	t.Helper()
	root, err := filepath.Abs("..")
	if err != nil {
		t.Fatal(err)
	}
	return root
}

// Every sequencer config shipped in configs/ and examples/ loads to a usable
// config, except the v2 templates documented as description only, which fail
// with an explicit error (configs/README.md).
func TestShippedSequencerConfigsLoad(t *testing.T) {
	root := repoRoot(t)
	var files []string
	for _, pattern := range []string{
		"configs/dev/*.json", "configs/sequencer/*.json", "configs/sequencer/v1.0/*.yaml",
		"configs/sequencers/*.yaml", "examples/*.yaml",
	} {
		m, err := filepath.Glob(filepath.Join(root, pattern))
		if err != nil {
			t.Fatal(err)
		}
		files = append(files, m...)
	}
	if len(files) < 15 {
		t.Fatalf("found only %d configs: %v", len(files), files)
	}
	descriptionOnly := map[string]bool{"pacbio_revio.yaml": true, "element_aviti.yaml": true, "mgi_dnbseq_t20x2.yaml": true}
	for _, f := range files {
		t.Run(filepath.Base(f), func(t *testing.T) {
			c, err := cli.NewConfigLoader(false).LoadSequencerConfig(f)
			if descriptionOnly[filepath.Base(f)] && strings.Contains(f, "configs/sequencers") {
				if err == nil || !strings.Contains(err.Error(), "only Oxford Nanopore templates") {
					t.Fatalf("err = %v, want the description-only error", err)
				}
				return
			}
			if err != nil {
				t.Fatalf("load: %v", err)
			}
			if c.SampleName == "" || c.ReadLength <= 0 || c.ReadsPerBatch <= 0 || c.TotalBatches <= 0 || c.OutputDir == "" {
				t.Errorf("incomplete config: sample %q, read length %d, %d x %d reads, output %q",
					c.SampleName, c.ReadLength, c.TotalBatches, c.ReadsPerBatch, c.OutputDir)
			}
		})
	}
}

// Every shipped phantom config loads and names only validation tests the
// framework knows, so enabling it does not fail with "Unknown ... test".
func TestShippedPhantomConfigsUseKnownTests(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join(repoRoot(t), "configs/phantom/v1.0/*.yaml"))
	if len(files) == 0 {
		t.Fatal("no phantom configs found")
	}
	const fastq = "@r1\nACGTACGTAC\n+\nIIIIIIIIII\n@r2\nTTGCAACGTA\n+\nIIIIIIIIII\n"
	for _, f := range files {
		t.Run(filepath.Base(f), func(t *testing.T) {
			p, err := cli.NewConfigLoader(false).LoadPhantomConfig(f)
			if err != nil {
				t.Fatalf("load: %v", err)
			}
			cfg, ok := cli.EffectiveValidation(p.Validation)
			if !ok {
				return
			}
			for _, r := range validation.NewFramework(cfg).RunValidation(fastq).Tests {
				if strings.HasPrefix(r.Message, "Unknown") {
					t.Errorf("%s: %s", r.TestName, r.Message)
				}
			}
		})
	}
}
