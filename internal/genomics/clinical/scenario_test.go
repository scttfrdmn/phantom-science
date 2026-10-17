//go:build !release

package clinical

import (
    "encoding/json"
    "os"
    "path/filepath"
    "testing"
    "time"

    "github.com/stretchr/testify/assert"
    "github.com/stretchr/testify/require"
)

func TestLoadScenario_YAML(t *testing.T) {
    yaml := `
phenotype:
  hpo_terms:
    - HP:0001250
    - HP:0001253
  resolved_gene: SLC2A1
  planted_variant:
    chrom: chr1
    pos: 43393429
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

    assert.Equal(t, []string{"HP:0001250", "HP:0001253"}, sc.Phenotype.HPOTerms)
    assert.Equal(t, "SLC2A1", sc.Phenotype.ResolvedGene)
    require.NotNil(t, sc.Phenotype.PlantedVariant)
    pv := sc.Phenotype.PlantedVariant
    assert.Equal(t, "chr1", pv.Chrom)
    assert.Equal(t, int64(43393429), pv.Pos)
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
            Pos:         43393429,
            Ref:         "C",
            Alt:         "T",
            Zygosity:    "heterozygous",
            ExpectedVAF: 0.5,
            Consequence: "missense_variant",
            ACMGClass:   "pathogenic",
            HPOTerms:    []string{"HP:0001250", "HP:0001253"},
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
            RegionFocused: "chr1:43368429-43418429",
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

func TestResolveHPO_GLUT1(t *testing.T) {
    spec, gene, err := ResolveHPO(nil, []string{"HP:0001250", "HP:0001253"}, nil)
    require.NoError(t, err)
    assert.Equal(t, "SLC2A1", gene)
    require.NotNil(t, spec)
    assert.Equal(t, "chr1", spec.Chrom)
    assert.Equal(t, int64(43393429), spec.Pos)
    assert.Equal(t, "C", spec.Ref)
    assert.Equal(t, "T", spec.Alt)
    assert.Equal(t, "pathogenic", spec.ACMGClass)
}

func TestResolveHPO_PKU(t *testing.T) {
    spec, gene, err := ResolveHPO(nil, []string{"HP:0001360"}, nil)
    require.NoError(t, err)
    assert.Equal(t, "PAH", gene)
    require.NotNil(t, spec)
    assert.Equal(t, "chr12", spec.Chrom)
}

func TestResolveHPO_CF(t *testing.T) {
    spec, gene, err := ResolveHPO(nil, []string{"HP:0002837"}, nil)
    require.NoError(t, err)
    assert.Equal(t, "CFTR", gene)
    require.NotNil(t, spec)
    assert.Equal(t, "chr7", spec.Chrom)
}

func TestResolveHPO_NoMatch(t *testing.T) {
    _, _, err := ResolveHPO(nil, []string{"HP:9999999"}, nil)
    assert.Error(t, err)
    assert.Contains(t, err.Error(), "planted_variant")
}

func TestResolveHPO_Superset(t *testing.T) {
    // Extra HPO terms beyond catalog keys should still resolve.
    spec, gene, err := ResolveHPO(nil, []string{
        "HP:0001250", "HP:0001253",
        "HP:0012345", // extra term not in catalog key
    }, nil)
    require.NoError(t, err)
    assert.Equal(t, "SLC2A1", gene)
    require.NotNil(t, spec)
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
