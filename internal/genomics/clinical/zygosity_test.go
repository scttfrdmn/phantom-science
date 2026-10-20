//go:build !release

package clinical_test

// Tests for issue #17 in clinical scenarios: exact homozygous VAF, indel-safe
// haplotype injection without 'N' padding, and seeded reproducibility.

import (
	"bufio"
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
	cfg := &types.SequencerConfig{SampleName: "hom", Lane: 1, ReadLength: 150, BaseQuality: 30, Seed: seed}
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
// reads keep their length, contain no N, and the de novo allele lands at its
// offset after the deletion to its left.
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

	// Reads start exactly 50 bases before the locus.
	region := &genomics.GenomicRegion{Chromosome: testChrom, Start: int64(locus0) - 50, End: int64(locus0) - 50 + readLength + 1}
	cfg := &types.SequencerConfig{SampleName: "indel", Lane: 1, ReadLength: readLength, BaseQuality: 30, Seed: 99}
	engine := genomics.NewEngine(*cfg)
	defer engine.Close()

	var buf bytes.Buffer
	bw := bufio.NewWriter(&buf)
	_, err = clinical.ExportGenerateHaplotypeFASTQ(bw, engine, cfg, trio, "proband", region, readLength, 400)
	require.NoError(t, err)
	require.NoError(t, bw.Flush())

	lines := strings.Split(strings.TrimSpace(buf.String()), "\n")
	require.Equal(t, 400*4, len(lines))
	var alt, ref int
	for i := 0; i < len(lines); i += 4 {
		seq, qual := lines[i+1], lines[i+3]
		require.Len(t, seq, readLength)
		require.Len(t, qual, readLength)
		require.NotContains(t, seq, "N")
		// Deletion at read offset 30 keeps its anchor base 'A'; the locus
		// (offset 50) moves 4 bases left; the insertion follows at 80-4.
		require.Equal(t, byte('A'), seq[30], "deletion anchor")
		require.Equal(t, "GCCCC", seq[76:81], "insertion allele")
		switch seq[46] {
		case testAlt[0]:
			alt++
		case testRef[0]:
			ref++
		default:
			t.Fatalf("locus base %q is neither REF nor ALT: %s", seq[46], seq)
		}
	}
	vaf := float64(alt) / float64(alt+ref)
	assert.InDelta(t, 0.5, vaf, 0.1, "proband de novo VAF")
}
