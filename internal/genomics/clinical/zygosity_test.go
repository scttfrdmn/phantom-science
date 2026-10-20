//go:build !release

package clinical_test

// Tests for issue #17 in clinical scenarios: exact homozygous VAF, indel-safe
// haplotype injection without 'N' padding, and seeded reproducibility.

import (
	"bytes"
	"context"
	mathrand "math/rand"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/clinical"
	"github.com/scttfrdmn/phantom-science/internal/population"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// runScenario runs a simple (non-trio) scenario planting one variant and
// returns the proband reads.
func runScenario(t *testing.T, zygosity string, seed uint64) ([]byte, [][2]string) {
	t.Helper()
	dir := t.TempDir()
	sc := &clinical.ClinicalScenario{
		Phenotype: clinical.PhenotypeSpec{PlantedVariant: &clinical.PlantedVariantSpec{
			Chrom: testChrom, Pos: testLocus1Based, Ref: testRef, Alt: testAlt, Zygosity: zygosity,
		}},
		Sequencer: clinical.SequencerSpec{Coverage: 200, ReadLength: 150},
		Output: clinical.OutputSpec{
			ProbandFASTQ: filepath.Join(dir, "proband.fastq"),
			GroundTruth:  filepath.Join(dir, "truth.json"),
		},
	}
	// Sequencing errors off: the test checks variant injection, and a
	// quality-matched substitution at the locus would change the VAF.
	cfg := &types.SequencerConfig{SampleName: "hom", Lane: 1, ReadLength: 150, BaseQuality: 30, Seed: seed,
		Illumina: types.IlluminaSpec{Errors: "none"}}
	m, err := clinical.NewScenarioGenerator(cfg).Run(context.Background(), sc)
	require.NoError(t, err)
	require.Len(t, m.PlantedVariants, 1)
	data, err := os.ReadFile(sc.Output.ProbandFASTQ)
	require.NoError(t, err)
	return data, parseFASTQ(data)
}

// locusAlleles counts ALT and REF bases at the planted locus.
func locusAlleles(t *testing.T, reads [][2]string) (alt, ref int) {
	t.Helper()
	locus0 := testLocus1Based - 1
	for _, rd := range reads {
		rs, err := parseReadStart(rd[0])
		require.NoError(t, err)
		off := locus0 - rs
		if off < 0 || int(off) >= len(rd[1]) {
			continue
		}
		switch rd[1][off] {
		case testAlt[0]:
			alt++
		case testRef[0]:
			ref++
		}
	}
	return alt, ref
}

func TestScenario_HomozygousVAFExactlyOne(t *testing.T) {
	_, reads := runScenario(t, "homozygous", 0)
	alt, ref := locusAlleles(t, reads)
	require.Greater(t, alt, 50, "need locus coverage")
	assert.Zero(t, ref, "homozygous scenario: %d REF reads at the locus, want 0 (VAF exactly 1.0)", ref)
}

func TestScenario_SeedReproducible(t *testing.T) {
	a, _ := runScenario(t, "heterozygous", 1234)
	b, _ := runScenario(t, "heterozygous", 1234)
	c, _ := runScenario(t, "heterozygous", 1235)
	assert.True(t, bytes.Equal(a, b), "same seed must give byte-identical FASTQ")
	assert.False(t, bytes.Equal(a, c), "different seeds should give different FASTQ")
}

// A deletion and an insertion on the haplotype next to the de novo SNV: the
// haplotypes carry exactly those alleles, and the reads sampled from them
// keep their length, contain no N, are substrings of a haplotype, and show
// the de novo allele at VAF 1/2.
func TestHaplotypeFASTQ_IndelsNoNPadding(t *testing.T) {
	const readLength = 100
	locus0 := int(testLocus1Based - 1)
	variants := []population.HaplotypeVariant{
		{Position: locus0 - 20, Alleles: []string{"ACGTA", "A"}}, // deletion -4, left of the locus
		{Position: locus0 + 30, Alleles: []string{"G", "GCCCC"}}, // insertion +4, right of the locus
	}
	matrix := make([][]int, 6)
	for i := range matrix {
		matrix[i] = []int{1, 1} // every haplotype carries both indels
	}
	result := &population.HaplotypeSimulationResult{Variants: variants, HaplotypeMatrix: matrix, NumSamples: 6}
	deNovo := &genomics.VariantRecord{Chromosome: testChrom, Position: testLocus1Based, Ref: testRef, Alt: testAlt, Freq: 1}
	// #nosec G404 -- test RNG
	trio, err := clinical.GenerateTrioHaplotypes(result, deNovo, mathrand.New(mathrand.NewSource(7)))
	require.NoError(t, err)

	// A 400-base window with the locus at offset 200.
	start := int64(locus0) - 200
	base := randomBases(3, 400)
	copy(base[180:], "ACGTA")
	base[200] = testRef[0]
	base[230] = 'G'
	haps, err := clinical.ExportMemberHaplotypes(base, start, trio, "proband")
	require.NoError(t, err)
	want := func(dn byte) string {
		return string(base[:181]) + string(base[185:200]) + string(dn) + string(base[201:231]) + "CCCC" + string(base[231:])
	}
	require.Equal(t, want(testAlt[0]), string(haps[0]), "proband H1: deletion, de novo ALT, insertion")
	require.Equal(t, want(testRef[0]), string(haps[1]), "proband H2: deletion, REF, insertion")

	reads := memberReads(t, base, start, trio, "proband", readLength, 2000, 99)
	require.Len(t, reads, 2000)
	for _, rd := range reads {
		seq := rd[1]
		require.Len(t, seq, readLength)
		require.NotContains(t, seq, "N")
		found := false
		for _, h := range haps {
			if strings.Contains(string(h), seq) || strings.Contains(string(h), revcomp(seq)) {
				found = true
			}
		}
		require.True(t, found, "read is not a substring of either haplotype: %s", seq)
	}
	// The locus is at offset 196 of the haplotypes, after the deletion.
	c := countContexts(reads, haps[0][186:207], haps[1][186:207])
	requireBinomial(t, c[0], c[0]+c[1], 0.5, "proband de novo")
}
