//go:build !release

package clinical

import (
	"context"
	"encoding/json"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestLoadScenario_YAML(t *testing.T) {
	yaml := `
phenotype:
  hpo_terms:
    - HP:0001250
    - HP:0011972
  resolved_gene: SLC2A1
  planted_variant:
    chrom: chr1
    pos: 42929242
    ref: C
    alt: T
    zygosity: heterozygous
    consequence: missense_variant
    acmg_class: pathogenic
trio:
  generate: true
  de_novo: true
sequencer:
  platform: NovaSeq6000
  coverage: 50
  read_length: 150
output:
  proband_fastq: proband.fastq.gz
  mother_fastq: mother.fastq.gz
  father_fastq: father.fastq.gz
  ground_truth: manifest.json
`
	dir := t.TempDir()
	path := filepath.Join(dir, "scenario.yaml")
	require.NoError(t, os.WriteFile(path, []byte(yaml), 0o600))

	sc, err := LoadScenario(path)
	require.NoError(t, err)

	assert.Equal(t, []string{"HP:0001250", "HP:0011972"}, sc.Phenotype.HPOTerms)
	assert.Equal(t, "SLC2A1", sc.Phenotype.ResolvedGene)
	require.NotNil(t, sc.Phenotype.PlantedVariant)
	pv := sc.Phenotype.PlantedVariant
	assert.Equal(t, "chr1", pv.Chrom)
	assert.Equal(t, int64(42929242), pv.Pos)
	assert.Equal(t, "C", pv.Ref)
	assert.Equal(t, "T", pv.Alt)
	assert.Equal(t, "heterozygous", pv.Zygosity)
	assert.Equal(t, "missense_variant", pv.Consequence)
	assert.Equal(t, "pathogenic", pv.ACMGClass)

	require.NotNil(t, sc.Trio)
	assert.True(t, sc.Trio.Generate)
	assert.True(t, sc.Trio.DeNovo)

	assert.Equal(t, "NovaSeq6000", sc.Sequencer.Platform)
	assert.Equal(t, 50, sc.Sequencer.Coverage)
	assert.Equal(t, 150, sc.Sequencer.ReadLength)

	assert.Equal(t, "proband.fastq.gz", sc.Output.ProbandFASTQ)
	assert.Equal(t, "mother.fastq.gz", sc.Output.MotherFASTQ)
	assert.Equal(t, "father.fastq.gz", sc.Output.FatherFASTQ)
	assert.Equal(t, "manifest.json", sc.Output.GroundTruth)
}

func TestLoadScenario_MissingRequired(t *testing.T) {
	dir := t.TempDir()

	// Missing hpo_terms and planted_variant.
	bad := `
sequencer:
  platform: NovaSeq6000
output:
  proband_fastq: p.fastq
  ground_truth: gt.json
`
	p := filepath.Join(dir, "bad.yaml")
	require.NoError(t, os.WriteFile(p, []byte(bad), 0o600))
	_, err := LoadScenario(p)
	assert.Error(t, err)
}

func TestManifest_RoundTrip(t *testing.T) {
	m := &GroundTruthManifest{
		GeneratedAt: time.Date(2026, 3, 11, 12, 0, 0, 0, time.UTC),
		ProbandID:   "proband-001",
		PlantedVariants: []PlantedVariant{{
			Gene:        "SLC2A1",
			Chrom:       "chr1",
			Pos:         42929242,
			Ref:         "C",
			Alt:         "T",
			Zygosity:    "heterozygous",
			ExpectedVAF: 0.5,
			Consequence: "missense_variant",
			ACMGClass:   "pathogenic",
			HPOTerms:    []string{"HP:0001250", "HP:0011972"},
		}},
		Trio: &TrioManifest{
			MotherCarries: false,
			FatherCarries: false,
			DeNovo:        true,
		},
		Sequencer: SequencerManifest{
			Platform:      "NovaSeq6000",
			Coverage:      50,
			ReadLength:    150,
			RegionFocused: "chr1:42904241-42954241",
		},
		ComputeMetrics: ComputeMetrics{
			ReadsGenerated: 10000,
			WallClockNS:    1234567890,
		},
	}

	dir := t.TempDir()
	path := filepath.Join(dir, "manifest.json")
	require.NoError(t, WriteManifest(path, m))

	data, err := os.ReadFile(path)
	require.NoError(t, err)

	var got GroundTruthManifest
	require.NoError(t, json.Unmarshal(data, &got))

	assert.Equal(t, m.ProbandID, got.ProbandID)
	require.Len(t, got.PlantedVariants, 1)
	assert.Equal(t, m.PlantedVariants[0].Gene, got.PlantedVariants[0].Gene)
	assert.Equal(t, m.PlantedVariants[0].Pos, got.PlantedVariants[0].Pos)
	assert.Equal(t, m.PlantedVariants[0].ExpectedVAF, got.PlantedVariants[0].ExpectedVAF)
	assert.Equal(t, m.PlantedVariants[0].HPOTerms, got.PlantedVariants[0].HPOTerms)
	require.NotNil(t, got.Trio)
	assert.True(t, got.Trio.DeNovo)
	assert.Equal(t, m.Sequencer.RegionFocused, got.Sequencer.RegionFocused)
	assert.Equal(t, m.ComputeMetrics.ReadsGenerated, got.ComputeMetrics.ReadsGenerated)
}

// The catalogue entries as verified against HPO (hp/releases/2026-09-01)
// and ClinVar on 2026-10-19; see the Source comment on hpoCatalog.
func TestResolveHPO_VerifiedEntries(t *testing.T) {
	cases := []struct {
		terms                []string
		gene, clinvar, chrom string
		pos                  int64
		ref, alt, csq        string
	}{
		{[]string{"HP:0001250", "HP:0011972"}, "SLC2A1", "VCV000016113", "chr1", 42929242, "C", "T", "missense_variant"},
		{[]string{"HP:0004923"}, "PAH", "VCV000000577", "chr12", 102840493, "G", "A", "missense_variant"},
		{[]string{"HP:0012236"}, "CFTR", "VCV000007105", "chr7", 117559590, "ATCT", "A", "inframe_deletion"},
	}
	for _, c := range cases {
		spec, entry, err := ResolveHPO(c.terms)
		require.NoError(t, err, "%v", c.terms)
		require.NotNil(t, spec)
		require.NotNil(t, entry)
		assert.Equal(t, c.gene, entry.Gene)
		assert.Equal(t, c.clinvar, entry.ClinVarID)
		assert.Equal(t, "GRCh38", entry.Build)
		assert.Equal(t, c.chrom, spec.Chrom)
		assert.Equal(t, c.pos, spec.Pos)
		assert.Equal(t, c.ref, spec.Ref)
		assert.Equal(t, c.alt, spec.Alt)
		assert.Equal(t, c.csq, spec.Consequence)
		assert.Equal(t, "pathogenic", spec.ACMGClass)
		assert.Equal(t, "heterozygous", spec.Zygosity)
	}
}

// Terms the old catalogue mis-mapped (HP:0001253 is an alt_id of Ataxia,
// HP:0001360 is Holoprosencephaly, HP:0002837 is Recurrent bronchitis) no
// longer resolve.
func TestResolveHPO_OldWrongTermsDoNotResolve(t *testing.T) {
	for _, terms := range [][]string{
		{"HP:0001250", "HP:0001253"},
		{"HP:0001360"},
		{"HP:0002837"},
	} {
		_, _, err := ResolveHPO(terms)
		assert.Error(t, err, "%v", terms)
	}
}

// Structural checks on every entry: labels parallel the terms, REF/ALT are
// VCF-normalised (an indel shares its anchor base and no other base), and the
// consequence agrees with the allele lengths.
func TestCatalog_EntriesWellFormed(t *testing.T) {
	cat := Catalog()
	require.Len(t, cat, 3)
	for _, e := range cat {
		assert.Len(t, e.HPOLabels, len(e.Terms), e.Gene)
		assert.NotEmpty(t, e.ClinVarID, e.Gene)
		assert.Regexp(t, `^VCV\d{9}$`, e.ClinVarID)
		assert.Equal(t, "GRCh38", e.Build)
		assert.Regexp(t, `^chr[0-9XY]+$`, e.Chrom)
		assert.Positive(t, e.Pos)
		assert.Regexp(t, `^[ACGT]+$`, e.Ref)
		assert.Regexp(t, `^[ACGT]+$`, e.Alt)
		switch {
		case len(e.Ref) == 1 && len(e.Alt) == 1:
			assert.NotEqual(t, e.Ref, e.Alt)
			assert.Equal(t, "missense_variant", e.Consequence, e.Gene)
		default:
			short, long := e.Alt, e.Ref
			if len(short) > len(long) {
				short, long = long, short
			}
			assert.Len(t, short, 1, "%s: indel must be anchor + inserted/deleted bases", e.Gene)
			assert.Equal(t, short[0], long[0], "%s: indel anchor base", e.Gene)
			if (len(long)-1)%3 == 0 {
				assert.Contains(t, e.Consequence, "inframe", e.Gene)
			} else {
				assert.Equal(t, "frameshift_variant", e.Consequence, e.Gene)
			}
		}
	}
}

func TestResolveHPO_NoMatch(t *testing.T) {
	_, _, err := ResolveHPO([]string{"HP:9999999"})
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "planted_variant")
}

func TestResolveHPO_Superset(t *testing.T) {
	// Extra HPO terms beyond catalog keys should still resolve.
	spec, entry, err := ResolveHPO([]string{
		"HP:0001250", "HP:0011972",
		"HP:0012345", // extra term not in catalog key
	})
	require.NoError(t, err)
	assert.Equal(t, "SLC2A1", entry.Gene)
	require.NotNil(t, spec)
}

// Seizure alone is not enough for the GLUT1 entry.
func TestResolveHPO_PartialTermsDoNotMatch(t *testing.T) {
	_, _, err := ResolveHPO([]string{"HP:0001250"})
	assert.Error(t, err)
}

func TestPlantedVariantSpec_ExpectedVAF(t *testing.T) {
	het := &PlantedVariantSpec{Zygosity: "heterozygous"}
	assert.Equal(t, 0.5, het.ExpectedVAF())

	hom := &PlantedVariantSpec{Zygosity: "homozygous"}
	assert.Equal(t, 1.0, hom.ExpectedVAF())

	// Default (empty) should be het-like.
	empty := &PlantedVariantSpec{}
	assert.Equal(t, 0.5, empty.ExpectedVAF())
}

// An HPO-resolved scenario records the catalogue provenance in the manifest.
func TestRun_HPOResolvedManifest(t *testing.T) {
	dir := t.TempDir()
	sc := &ClinicalScenario{
		Phenotype: PhenotypeSpec{HPOTerms: []string{"HP:0012236"}},
		Sequencer: SequencerSpec{Coverage: 2, ReadLength: 100},
		Output: OutputSpec{
			ProbandFASTQ: filepath.Join(dir, "p.fastq"),
			GroundTruth:  filepath.Join(dir, "m.json"),
		},
	}
	cfg := &types.SequencerConfig{SampleName: "cf", Lane: 1, ReadLength: 100, BaseQuality: 30, Seed: 5}
	m, err := NewScenarioGenerator(cfg, dir).Run(context.Background(), sc)
	require.NoError(t, err)
	require.Len(t, m.PlantedVariants, 1)
	v := m.PlantedVariants[0]
	assert.Equal(t, "CFTR", v.Gene)
	assert.Equal(t, "VCV000007105", v.ClinVarID)
	assert.Equal(t, "GRCh38", v.Build)
	assert.Equal(t, "ATCT", v.Ref)
	assert.Equal(t, "A", v.Alt)
	assert.Equal(t, "inframe_deletion", v.Consequence)
	assert.Contains(t, v.HGVS, "p.Phe508del")
}
