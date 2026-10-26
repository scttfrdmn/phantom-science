//go:build !release

package transcriptome

import (
	"io"
	"path/filepath"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
)

// A run written from a RunSpec is regenerated read for read from the spec
// its manifest records (#31), and the manifest describes every output.
func TestRunSpec_RegenerateAndManifest(t *testing.T) {
	spec := RunSpec{
		Reads: 250, Tier: "hac", Kit: DefaultKit, NovelRate: 0.3, FullLengthFraction: DefaultFullLengthFraction,
		Zipf:       ZipfModel{Alpha: DefaultZipfAlpha, Sigma: DefaultExpressionSigma},
		GenomePath: filepath.Join("testdata", "tiny_genome.fna"), GTFPath: filepath.Join("testdata", "tiny.gtf"),
	}
	const seed = 21
	sim, _, err := spec.Simulator(seed, io.Discard)
	if err != nil {
		t.Fatal(err)
	}
	dir := t.TempDir()
	res, err := WriteRun(dir, sim, OutputOptions{})
	if err != nil {
		t.Fatal(err)
	}
	rec, err := spec.Record()
	if err != nil {
		t.Fatal(err)
	}
	if err := WriteManifest("test", seed, nil, rec, sim, res); err != nil {
		t.Fatal(err)
	}
	m, err := manifest.Read(filepath.Join(dir, manifest.FileName))
	if err != nil {
		t.Fatal(err)
	}
	if err := m.VerifyDir(dir); err != nil {
		t.Fatal(err)
	}
	if m.Seed != seed || m.Fingerprint != nil || m.Transcriptome.Genome.Name != "tiny_genome.fna" || m.Totals.Reads != 250 || m.Totals.Bases != res.Bases {
		t.Fatalf("manifest %+v %+v", m, m.Transcriptome)
	}

	back := SpecFromRecord(m.Transcriptome)
	back.GenomePath, back.GTFPath = spec.GenomePath, spec.GTFPath
	written := readFASTQFile(t, res.FASTQ)
	i := 0
	n, err := Regenerate(back, seed, func(seq []byte) error {
		if i >= len(written) || string(seq) != written[i].seq {
			t.Fatalf("regenerated read %d differs", i)
		}
		i++
		return nil
	})
	if err != nil || n != len(written) || i != n {
		t.Fatalf("regenerated %d of %d reads: %v", n, len(written), err)
	}
}
