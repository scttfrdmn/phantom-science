//go:build !release

// Package clinical implements clinical scenario-driven synthetic FASTQ generation.
// Scenarios are YAML-specified phenotype+sequencer configurations that produce
// positionally-accurate reads with planted variants and machine-readable ground truth.
package clinical

import (
	"fmt"
	"os"

	"gopkg.in/yaml.v3"
)

// ClinicalScenario is the top-level YAML type for a clinical demonstration scenario.
type ClinicalScenario struct {
	Phenotype PhenotypeSpec `yaml:"phenotype" json:"phenotype"`
	Trio      *TrioSpec     `yaml:"trio,omitempty" json:"trio,omitempty"`
	Sequencer SequencerSpec `yaml:"sequencer" json:"sequencer"`
	Output    OutputSpec    `yaml:"output" json:"-"`
}

// PhenotypeSpec describes the phenotype and optional pre-specified variant.
type PhenotypeSpec struct {
	HPOTerms       []string            `yaml:"hpo_terms" json:"hpo_terms,omitempty"`
	ResolvedGene   string              `yaml:"resolved_gene,omitempty" json:"resolved_gene,omitempty"`
	PlantedVariant *PlantedVariantSpec `yaml:"planted_variant,omitempty" json:"planted_variant,omitempty"`
}

// PlantedVariantSpec is an explicit variant to plant (overrides HPO resolution).
type PlantedVariantSpec struct {
	Chrom       string `yaml:"chrom" json:"chrom"`
	Pos         int64  `yaml:"pos" json:"pos"`
	Ref         string `yaml:"ref" json:"ref"`
	Alt         string `yaml:"alt" json:"alt"`
	Zygosity    string `yaml:"zygosity" json:"zygosity"` // "heterozygous" | "homozygous"
	Consequence string `yaml:"consequence,omitempty" json:"consequence,omitempty"`
	ACMGClass   string `yaml:"acmg_class,omitempty" json:"acmg_class,omitempty"`
}

// ExpectedVAF returns 0.5 for heterozygous and 1.0 for homozygous.
func (p *PlantedVariantSpec) ExpectedVAF() float64 {
	if p.Zygosity == "homozygous" {
		return 1.0
	}
	return 0.5
}

// TrioSpec controls optional trio generation.
type TrioSpec struct {
	Generate bool `yaml:"generate" json:"generate"`
	DeNovo   bool `yaml:"de_novo" json:"de_novo"`
}

// SequencerSpec describes the sequencing parameters. Platform is a label
// recorded in the manifest; the read model is the sequencer config's
// (Illumina by default: binned qualities, quality-matched errors).
type SequencerSpec struct {
	Platform   string `yaml:"platform" json:"platform"`
	Coverage   int    `yaml:"coverage" json:"coverage"`
	ReadLength int    `yaml:"read_length" json:"read_length"`
	// PairedEnd writes read pairs, the two ends of one fragment: to the
	// *_fastq and *_fastq_r2 files, or interleaved in *_fastq when no R2
	// file is given.
	PairedEnd bool `yaml:"paired_end,omitempty" json:"paired_end,omitempty"`
}

// OutputSpec specifies output file paths. The *R2 paths are the read-2
// files of paired-end scenarios.
type OutputSpec struct {
	ProbandFASTQ   string `yaml:"proband_fastq"`
	ProbandFASTQR2 string `yaml:"proband_fastq_r2,omitempty"`
	MotherFASTQ    string `yaml:"mother_fastq,omitempty"`
	MotherFASTQR2  string `yaml:"mother_fastq_r2,omitempty"`
	FatherFASTQ    string `yaml:"father_fastq,omitempty"`
	FatherFASTQR2  string `yaml:"father_fastq_r2,omitempty"`
	GroundTruth    string `yaml:"ground_truth"`
}

// LoadScenario reads and unmarshals a YAML scenario specification.
func LoadScenario(path string) (*ClinicalScenario, error) {
	data, err := os.ReadFile(path) // #nosec G304 -- path from CLI flag
	if err != nil {
		return nil, fmt.Errorf("reading scenario file: %w", err)
	}
	var sc ClinicalScenario
	if err := yaml.Unmarshal(data, &sc); err != nil {
		return nil, fmt.Errorf("parsing scenario YAML: %w", err)
	}
	if len(sc.Phenotype.HPOTerms) == 0 && sc.Phenotype.PlantedVariant == nil {
		return nil, fmt.Errorf("scenario must specify hpo_terms or planted_variant")
	}
	if sc.Trio != nil && sc.Trio.Generate && !sc.Trio.DeNovo {
		return nil, fmt.Errorf("trio.generate requires trio.de_novo: true (only de novo trios are generated)")
	}
	if sc.Output.ProbandFASTQ == "" {
		return nil, fmt.Errorf("scenario must specify output.proband_fastq")
	}
	if sc.Output.GroundTruth == "" {
		return nil, fmt.Errorf("scenario must specify output.ground_truth")
	}
	if !sc.Sequencer.PairedEnd {
		for _, f := range []struct{ name, path string }{
			{"proband_fastq_r2", sc.Output.ProbandFASTQR2},
			{"mother_fastq_r2", sc.Output.MotherFASTQR2},
			{"father_fastq_r2", sc.Output.FatherFASTQR2},
		} {
			if f.path != "" {
				return nil, fmt.Errorf("output.%s needs sequencer.paired_end: true", f.name)
			}
		}
	}
	return &sc, nil
}
