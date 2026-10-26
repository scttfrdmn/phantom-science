//go:build !release

package clinical

import (
	"bufio"
	"context"
	"os"
	"path/filepath"
	"testing"

	runmanifest "github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// fastqSeqs returns the sequence lines of a plain FASTQ file.
func fastqSeqs(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	var out []string
	sc := bufio.NewScanner(f)
	for i := 0; sc.Scan(); i++ {
		if i%4 == 1 {
			out = append(out, sc.Text())
		}
	}
	return out
}

// Regenerate reproduces every read of a scenario, in order, from the run
// manifest's records alone (#31): here a proband-only synthetic scenario
// with interleaved pairs, whose run manifest records the spec, the members
// and the seed.
func TestRegenerate_FromRunManifest(t *testing.T) {
	dir := t.TempDir()
	sc := &ClinicalScenario{
		Phenotype: PhenotypeSpec{PlantedVariant: &PlantedVariantSpec{Chrom: "chr7", Pos: 30001, Ref: "A", Alt: "G", Zygosity: "heterozygous"}},
		Sequencer: SequencerSpec{Platform: "NovaSeq6000", Coverage: 3, ReadLength: 100, PairedEnd: true},
		Output:    OutputSpec{ProbandFASTQ: filepath.Join(dir, "p.fastq"), GroundTruth: filepath.Join(dir, "truth.json")},
	}
	cfg := &types.SequencerConfig{SampleName: "P", Lane: 1, ReadLength: 100, Seed: 12}
	g := NewScenarioGenerator(cfg)
	gt, err := g.Run(context.Background(), sc)
	if err != nil {
		t.Fatal(err)
	}
	if err := WriteManifest(sc.Output.GroundTruth, gt); err != nil {
		t.Fatal(err)
	}
	manPath := filepath.Join(dir, runmanifest.FileName)
	if err := g.WriteRunManifest(manPath, sc.Output.GroundTruth, "test"); err != nil {
		t.Fatal(err)
	}
	m, err := runmanifest.Read(manPath)
	if err != nil {
		t.Fatal(err)
	}
	if err := m.VerifyDir(dir); err != nil {
		t.Fatal(err)
	}
	written := fastqSeqs(t, sc.Output.ProbandFASTQ)
	if m.Seed != 12 || m.Fingerprint != nil || len(m.Outputs) != 2 || m.Outputs[0].Member != "proband" ||
		m.Outputs[0].Reads != int64(len(written)) || m.Outputs[0].Bases != int64(100*len(written)) {
		t.Fatalf("manifest %+v", m)
	}

	c := m.Config
	c.Seed = m.Seed
	i := 0
	n, err := Regenerate(context.Background(), c, m.Scenario, func(member string, seq []byte) error {
		if member != "proband" || i >= len(written) || string(seq) != written[i] {
			t.Fatalf("regenerated read %d (%s) differs", i, member)
		}
		i++
		return nil
	})
	if err != nil || n != len(written) || i != n {
		t.Fatalf("regenerated %d of %d reads: %v", n, len(written), err)
	}
	c.Seed = 0
	if _, err := Regenerate(context.Background(), c, m.Scenario, func(string, []byte) error { return nil }); err == nil {
		t.Error("regenerated without a seed")
	}
}
