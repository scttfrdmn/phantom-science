//go:build !release

package clinical

import (
	"bytes"
	"context"
	mathrand "math/rand"
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Haplotype-level Mendelian transmission on simulated haplotypes: each
// proband allele comes from the corresponding parent, at every site.
func TestTrioHaplotypes_MendelianFromSimulation(t *testing.T) {
	sites := 0
	for seed := uint64(1); seed <= 20; seed++ {
		res, err := population.SimulateHaplotypes("YRI", "chr1", 1000, 50000, 6, seed)
		require.NoError(t, err)
		// #nosec G404 G115 -- test RNG
		trio, err := GenerateTrioHaplotypes(res, nil, mathrand.New(mathrand.NewSource(int64(seed))))
		require.NoError(t, err)
		for i := range trio.Variants {
			sites++
			assert.Contains(t, []int{trio.MotherH1[i], trio.MotherH2[i]}, trio.ProbandH1[i], "seed %d site %d: maternal allele", seed, i)
			assert.Contains(t, []int{trio.FatherH1[i], trio.FatherH2[i]}, trio.ProbandH2[i], "seed %d site %d: paternal allele", seed, i)
		}
	}
	require.Greater(t, sites, 100)
}

// Without crossovers each parent transmits either haplotype with probability
// 1/2 (the proband is not always given the parent's first haplotype).
func TestSimulateCrossover_UnbiasedTransmission(t *testing.T) {
	h1, h2 := []int{1, 1, 1}, []int{0, 0, 0}
	pos := []int{10, 20, 30}
	// #nosec G404 -- test RNG
	rng := mathrand.New(mathrand.NewSource(1))
	const n = 2000
	fromH1 := 0
	for i := 0; i < n; i++ {
		got := simulateCrossover(h1, h2, pos, 0, rng)
		require.True(t, got[0] == got[1] && got[1] == got[2], "no crossover: one whole haplotype")
		fromH1 += got[0]
	}
	assert.InDelta(t, 0.5, float64(fromH1)/n, 0.05)
}

// A simulated background variant at the de novo position is dropped, so it
// cannot replace the planted allele.
func TestGenerateTrioHaplotypes_DropsOverlapWithDeNovo(t *testing.T) {
	res := &population.HaplotypeSimulationResult{
		Variants: []population.HaplotypeVariant{
			{Position: 99, Alleles: []string{"A", "G"}},
			{Position: 100, Alleles: []string{"C", "A"}}, // the de novo base (1-based 101)
			{Position: 103, Alleles: []string{"G", "T"}}, // inside a 4-base REF span 100-103
			{Position: 104, Alleles: []string{"T", "C"}},
		},
		HaplotypeMatrix: [][]int{{1, 1, 1, 1}, {0, 0, 0, 0}, {1, 1, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
	}
	dn := &genomics.VariantRecord{Chromosome: "chr1", Position: 101, Ref: "CTTG", Alt: "C"}
	// #nosec G404 -- test RNG
	trio, err := GenerateTrioHaplotypes(res, dn, mathrand.New(mathrand.NewSource(1)))
	require.NoError(t, err)
	var kept []int
	for _, v := range trio.Variants {
		kept = append(kept, v.Position)
	}
	assert.Equal(t, []int{99, 104}, kept)
	assert.Len(t, trio.MotherH1, 2)
	assert.Len(t, trio.ProbandH2, 2)
}

// runTrioScenario runs a de novo trio scenario around the SLC2A1 catalogue
// variant and returns the manifest and the three FASTQs.
func runTrioScenario(t *testing.T, seed uint64, pop string) (*GroundTruthManifest, map[string][]byte) {
	t.Helper()
	dir := t.TempDir()
	sc := &ClinicalScenario{
		Phenotype: PhenotypeSpec{HPOTerms: []string{"HP:0001250", "HP:0011972"}},
		Trio:      &TrioSpec{Generate: true, DeNovo: true},
		Sequencer: SequencerSpec{Coverage: 12, ReadLength: 100},
		Output: OutputSpec{
			ProbandFASTQ: filepath.Join(dir, "proband.fastq"),
			MotherFASTQ:  filepath.Join(dir, "mother.fastq"),
			FatherFASTQ:  filepath.Join(dir, "father.fastq"),
			GroundTruth:  filepath.Join(dir, "truth.json"),
		},
	}
	cfg := &types.SequencerConfig{SampleName: "trio", Lane: 1, ReadLength: 100, BaseQuality: 30, Seed: seed, Population: pop}
	m, err := NewScenarioGenerator(cfg).Run(context.Background(), sc)
	require.NoError(t, err)
	out := map[string][]byte{}
	for member, path := range map[string]string{
		"proband": sc.Output.ProbandFASTQ, "mother": sc.Output.MotherFASTQ, "father": sc.Output.FatherFASTQ,
	} {
		data, err := os.ReadFile(path)
		require.NoError(t, err)
		require.NotEmpty(t, data, member)
		out[member] = data
	}
	return m, out
}

// The scenario-level trio from the catalogue: deterministic for a seed, the
// de novo ALT only in the proband, binned qualities for every member. Runs
// with nothing on PATH: no Python or other external process is involved.
// (The haplotype content, Mendelian consistency and read origin are checked
// in TestScenarioReference_Trio, with and without a reference.)
func TestScenario_TrioPureGo(t *testing.T) {
	t.Setenv("PATH", t.TempDir())
	const seed, pop = 2026, "CEU"

	m, fq := runTrioScenario(t, seed, pop)
	require.NotNil(t, m.Trio)
	assert.True(t, m.Trio.DeNovo)
	assert.Equal(t, "CEU", m.Trio.Population)
	assert.Equal(t, "CEU", m.Trio.DemographicModel)
	assert.False(t, m.Trio.MotherCarries)
	assert.False(t, m.Trio.FatherCarries)
	assert.Equal(t, basesSynthetic, m.Sequencer.Bases)

	m2, fq2 := runTrioScenario(t, seed, pop)
	for member := range fq {
		assert.True(t, bytes.Equal(fq[member], fq2[member]), "%s FASTQ must be byte-identical for the same seed", member)
	}
	assert.Equal(t, m.Trio.BackgroundVariants, m2.Trio.BackgroundVariants)
	_, fq3 := runTrioScenario(t, seed+1, pop)
	assert.False(t, bytes.Equal(fq["proband"], fq3["proband"]), "a different seed gives different reads")

	// De novo SLC2A1 chr1:42929242 C>T: rebuild the window to get its
	// context on the proband's haplotypes.
	cfg := types.SequencerConfig{SampleName: "trio", Lane: 1, ReadLength: 100, Seed: seed, Population: pop}
	dn0 := int64(42929242 - 1)
	region := &genomics.GenomicRegion{Chromosome: "chr1", Start: dn0 - windowHalfSize, End: dn0 + windowHalfSize}
	tp, err := NewScenarioGenerator(&cfg).planTrio(cfg, region, "chr1",
		&genomics.VariantRecord{Chromosome: "chr1", Position: 42929242, Ref: "C", Alt: "T", Freq: 0.5})
	require.NoError(t, err)
	haps, err := tp.memberHaplotypes("proband")
	require.NoError(t, err)
	off := windowHalfSize
	altCtx, refCtx := string(haps[0][off-12:off+13]), string(haps[1][off-12:off+13])
	require.Equal(t, byte('T'), altCtx[12])
	require.Equal(t, byte('C'), refCtx[12])
	for member, data := range fq {
		recs := parseRecords(t, data)
		requireBinnedQualities(t, recs)
		n := countCovering(recs, altCtx, refCtx)
		if member == "proband" {
			// Sequencing errors are on here, so only a loose check.
			assert.Greater(t, n[0], 0, "proband carries the de novo ALT")
			assert.Greater(t, n[1], 0, "proband is heterozygous for the de novo")
		} else {
			assert.Zero(t, n[0], "%s must not carry the de novo ALT", member)
		}
	}
}

func TestScenario_TrioUnsupportedPopulation(t *testing.T) {
	dir := t.TempDir()
	sc := &ClinicalScenario{
		Phenotype: PhenotypeSpec{HPOTerms: []string{"HP:0004923"}},
		Trio:      &TrioSpec{Generate: true, DeNovo: true},
		Sequencer: SequencerSpec{Coverage: 1, ReadLength: 100},
		Output: OutputSpec{
			ProbandFASTQ: filepath.Join(dir, "p.fastq"),
			GroundTruth:  filepath.Join(dir, "m.json"),
		},
	}
	cfg := &types.SequencerConfig{SampleName: "t", Lane: 1, ReadLength: 100, BaseQuality: 30, Seed: 1, Population: "FIN"}
	_, err := NewScenarioGenerator(cfg).Run(context.Background(), sc)
	require.Error(t, err)
	assert.Contains(t, err.Error(), `"FIN"`)
	assert.Contains(t, err.Error(), "supported:")
	assert.NoFileExists(t, sc.Output.ProbandFASTQ)
}

func TestLoadScenario_R2NeedsPairedEnd(t *testing.T) {
	path := filepath.Join(t.TempDir(), "s.yaml")
	require.NoError(t, os.WriteFile(path, []byte(`
phenotype: {hpo_terms: [HP:0004923]}
output: {proband_fastq: p1.fastq, proband_fastq_r2: p2.fastq, ground_truth: m.json}
`), 0o600))
	_, err := LoadScenario(path)
	require.Error(t, err)
	assert.Contains(t, err.Error(), "paired_end")
}

func TestLoadScenario_TrioRequiresDeNovo(t *testing.T) {
	path := filepath.Join(t.TempDir(), "s.yaml")
	require.NoError(t, os.WriteFile(path, []byte(`
phenotype: {hpo_terms: [HP:0004923]}
trio: {generate: true, de_novo: false}
output: {proband_fastq: p.fastq, ground_truth: m.json}
`), 0o600))
	_, err := LoadScenario(path)
	require.Error(t, err)
	assert.Contains(t, err.Error(), "de_novo")
}
