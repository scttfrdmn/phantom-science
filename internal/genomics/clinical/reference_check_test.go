//go:build !release

package clinical

import (
	"compress/gzip"
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// tinyFASTA is a two-contig reference with 10-base lines. chr1 positions
// (1-based): 1-10 ACGTACGTAC, 11-20 GGGGCCCCAA, 21-30 ttttgggggc (soft-masked).
const tinyFASTA = ">chr1 test contig\n" +
	"ACGTACGTAC\n" +
	"GGGGCCCCAA\n" +
	"ttttgggggc\n" +
	">7\n" +
	"AATCTTTGGC\n"

func writeTinyFASTA(t *testing.T, withIndex bool) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "tiny.fa")
	require.NoError(t, os.WriteFile(path, []byte(tinyFASTA), 0o600))
	if withIndex {
		fp, err := reference.NewFastaParser(path)
		require.NoError(t, err)
		idx, err := fp.BuildIndex()
		require.NoError(t, err)
		require.NoError(t, fp.SaveIndex(idx))
	}
	return path
}

func TestCheckPlantedRef_Match(t *testing.T) {
	for _, withIndex := range []bool{false, true} {
		ref := writeTinyFASTA(t, withIndex)
		ok := []*PlantedVariantSpec{
			{Chrom: "chr1", Pos: 1, Ref: "A", Alt: "G"},
			{Chrom: "chr1", Pos: 11, Ref: "G", Alt: "T"},   // first base of a line
			{Chrom: "chr1", Pos: 9, Ref: "ACGG", Alt: "A"}, // spans a line break
			{Chrom: "chr1", Pos: 21, Ref: "T", Alt: "C"},   // soft-masked base
			{Chrom: "chr7", Pos: 2, Ref: "ATCT", Alt: "A"}, // "chr7" found as "7"
			{Chrom: "1", Pos: 30, Ref: "C", Alt: "A"},      // "1" found as "chr1"; last base
		}
		assert.NoError(t, CheckPlantedRef(ref, ok...), "withIndex=%v", withIndex)
	}
}

func TestCheckPlantedRef_Mismatch(t *testing.T) {
	ref := writeTinyFASTA(t, true)
	cases := []struct {
		v    *PlantedVariantSpec
		want string
	}{
		{&PlantedVariantSpec{Chrom: "chr1", Pos: 2, Ref: "A", Alt: "G"}, `has "C" at chr1:2-2`},
		{&PlantedVariantSpec{Chrom: "chr7", Pos: 2, Ref: "CTT", Alt: "C"}, `has "ATC"`},
		{&PlantedVariantSpec{Chrom: "chr2", Pos: 5, Ref: "A", Alt: "G"}, `chromosome "chr2" is not in reference`},
		{&PlantedVariantSpec{Chrom: "chr1", Pos: 30, Ref: "CA", Alt: "C"}, "invalid coordinates"},
		{&PlantedVariantSpec{Chrom: "chr1", Pos: 0, Ref: "A", Alt: "G"}, "no REF allele or position"},
	}
	for _, c := range cases {
		err := CheckPlantedRef(ref, c.v)
		require.Error(t, err, "%+v", c.v)
		assert.Contains(t, err.Error(), c.want)
	}
}

func TestCheckPlantedRef_Gzip(t *testing.T) {
	path := filepath.Join(t.TempDir(), "tiny.fa.gz")
	f, err := os.Create(path)
	require.NoError(t, err)
	zw := gzip.NewWriter(f)
	_, err = zw.Write([]byte(tinyFASTA))
	require.NoError(t, err)
	require.NoError(t, zw.Close())
	require.NoError(t, f.Close())

	assert.NoError(t, CheckPlantedRef(path, &PlantedVariantSpec{Chrom: "chr1", Pos: 12, Ref: "GGG", Alt: "G"}))
	assert.Error(t, CheckPlantedRef(path, &PlantedVariantSpec{Chrom: "chr1", Pos: 12, Ref: "C", Alt: "G"}))
}

func TestCheckPlantedRef_MissingFile(t *testing.T) {
	err := CheckPlantedRef(filepath.Join(t.TempDir(), "absent.fa"),
		&PlantedVariantSpec{Chrom: "chr1", Pos: 1, Ref: "A", Alt: "G"})
	assert.Error(t, err)
}

// Run with a reference: a matching REF runs and is recorded in the manifest;
// a mismatching REF fails before any reads are written.
func TestRun_ReferenceCheck(t *testing.T) {
	ref := writeTinyFASTA(t, false)
	run := func(refAllele string) (*GroundTruthManifest, string, error) {
		dir := t.TempDir()
		sc := &ClinicalScenario{
			Phenotype: PhenotypeSpec{PlantedVariant: &PlantedVariantSpec{
				Chrom: "chr1", Pos: 15, Ref: refAllele, Alt: "T", Zygosity: "heterozygous",
			}},
			// The reads come from the 30-base chr1, so they are shorter.
			Sequencer: SequencerSpec{Coverage: 2, ReadLength: 20},
			Output: OutputSpec{
				ProbandFASTQ: filepath.Join(dir, "p.fastq"),
				GroundTruth:  filepath.Join(dir, "m.json"),
			},
		}
		cfg := &types.SequencerConfig{SampleName: "ref", Lane: 1, ReadLength: 20, BaseQuality: 30, Seed: 3}
		m, err := NewScenarioGenerator(cfg).WithReference(ref).Run(context.Background(), sc)
		return m, sc.Output.ProbandFASTQ, err
	}

	m, fq, err := run("C")
	require.NoError(t, err)
	assert.Equal(t, ref, m.Reference)
	assert.FileExists(t, fq)

	_, fq, err = run("A")
	require.Error(t, err)
	assert.True(t, strings.Contains(err.Error(), `REF "A" does not match`), err.Error())
	assert.NoFileExists(t, fq, "no reads are written when the REF check fails")
}
