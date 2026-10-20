//go:build !release

package clinical

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	mathrand "math/rand"
	"os"
	"path/filepath"
	"strconv"
	"strings"
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

// readsByStart parses haplotype-FASTQ records into (read start, sequence).
func readsByStart(t *testing.T, data []byte) [][2]string {
	t.Helper()
	var out [][2]string
	sc := bufio.NewScanner(bytes.NewReader(data))
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	var lines []string
	for sc.Scan() {
		lines = append(lines, sc.Text())
	}
	require.Zero(t, len(lines)%4)
	for i := 0; i < len(lines); i += 4 {
		// @name:lane:HXXXXXX:1:1101:X:Y 1:N:0:1 with start = Y*100000 + X
		f := strings.Split(strings.SplitN(lines[i][1:], " ", 2)[0], ":")
		require.GreaterOrEqual(t, len(f), 7)
		x, err := strconv.ParseInt(f[5], 10, 64)
		require.NoError(t, err)
		y, err := strconv.ParseInt(f[6], 10, 64)
		require.NoError(t, err)
		out = append(out, [2]string{fmt.Sprint(y*100000 + x), lines[i+1]})
	}
	return out
}

// allelesAt returns the bases observed at 0-based position pos.
func allelesAt(t *testing.T, reads [][2]string, pos int64) map[byte]int {
	t.Helper()
	seen := map[byte]int{}
	for _, r := range reads {
		start, err := strconv.ParseInt(r[0], 10, 64)
		require.NoError(t, err)
		off := pos - start
		if off >= 0 && off < int64(len(r[1])) {
			seen[r[1][off]]++
		}
	}
	return seen
}

// The scenario-level trio: deterministic for a seed, the de novo ALT appears
// only in the proband, and at every simulated background site the bases each
// member shows are exactly the alleles on that member's haplotypes, with the
// proband's alleles drawn from its parents (Mendelian consistency). Runs with
// nothing on PATH: no Python or other external process is involved.
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

	m2, fq2 := runTrioScenario(t, seed, pop)
	for member := range fq {
		assert.True(t, bytes.Equal(fq[member], fq2[member]), "%s FASTQ must be byte-identical for the same seed", member)
	}
	assert.Equal(t, m.Trio.BackgroundVariants, m2.Trio.BackgroundVariants)
	_, fq3 := runTrioScenario(t, seed+1, pop)
	assert.False(t, bytes.Equal(fq["proband"], fq3["proband"]), "a different seed gives different reads")

	reads := map[string][][2]string{}
	for member, data := range fq {
		reads[member] = readsByStart(t, data)
	}

	// De novo SLC2A1 chr1:42929242 C>T: ALT in the proband only.
	dn := int64(42929242 - 1)
	for member, rs := range reads {
		obs := allelesAt(t, rs, dn)
		if member == "proband" {
			assert.Greater(t, obs['T'], 0, "proband carries the de novo ALT")
			assert.Greater(t, obs['C'], 0, "proband is heterozygous for the de novo")
		} else {
			assert.Zero(t, obs['T'], "%s must not carry the de novo ALT", member)
			assert.Greater(t, obs['C'], 0, "%s shows REF at the de novo site", member)
		}
	}

	// Recreate the haplotypes the run used and check every background site.
	var cfgSeed uint64 = seed
	region0 := int64(42929242-1) - windowHalfSize
	res, err := population.SimulateHaplotypes(pop, "chr1", int(region0), 2*windowHalfSize, 6, cfgSeed^trioSimSeedSalt)
	require.NoError(t, err)
	require.Equal(t, m.Trio.BackgroundVariants, res.NumVariants)
	// #nosec G404 G115 -- same derivation as the generator
	trio, err := GenerateTrioHaplotypes(res, &genomics.VariantRecord{Chromosome: "chr1", Position: 42929242, Ref: "C", Alt: "T"},
		mathrand.New(mathrand.NewSource(int64(cfgSeed))))
	require.NoError(t, err)
	require.NotEmpty(t, trio.Variants)

	haps := map[string][2][]int{
		"mother":  {trio.MotherH1, trio.MotherH2},
		"father":  {trio.FatherH1, trio.FatherH2},
		"proband": {trio.ProbandH1, trio.ProbandH2},
	}
	checked := 0
	for i, v := range trio.Variants {
		ref, alt := v.Alleles[0][0], v.Alleles[1][0]
		allele := func(x int) byte {
			if x == 1 {
				return alt
			}
			return ref
		}
		parental := map[byte]bool{}
		for _, member := range []string{"mother", "father"} {
			parental[allele(haps[member][0][i])] = true
			parental[allele(haps[member][1][i])] = true
		}
		for member, h := range haps {
			obs := allelesAt(t, reads[member], int64(v.Position))
			want := map[byte]bool{allele(h[0][i]): true, allele(h[1][i]): true}
			for b := range obs {
				require.True(t, want[b], "%s site %d: base %q not on its haplotypes", member, v.Position, b)
				if member == "proband" {
					require.True(t, parental[b], "proband site %d: allele %q absent from both parents", v.Position, b)
				}
			}
		}
		checked++
	}
	t.Logf("checked %d background sites; manifest reports %d", checked, m.Trio.BackgroundVariants)
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
