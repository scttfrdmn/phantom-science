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
	Phenotype PhenotypeSpec `yaml:"phenotype"`
	Trio      *TrioSpec     `yaml:"trio,omitempty"`
	Sequencer SequencerSpec `yaml:"sequencer"`
	Output    OutputSpec    `yaml:"output"`
}

// PhenotypeSpec describes the phenotype and optional pre-specified variant.
type PhenotypeSpec struct {
	HPOTerms       []string            `yaml:"hpo_terms"`
	ResolvedGene   string              `yaml:"resolved_gene,omitempty"`
	PlantedVariant *PlantedVariantSpec `yaml:"planted_variant,omitempty"`
}

// PlantedVariantSpec is an explicit variant to plant (overrides HPO resolution).
type PlantedVariantSpec struct {
	Chrom       string `yaml:"chrom"`
	Pos         int64  `yaml:"pos"`
	Ref         string `yaml:"ref"`
	Alt         string `yaml:"alt"`
	Zygosity    string `yaml:"zygosity"` // "heterozygous" | "homozygous"
	Consequence string `yaml:"consequence,omitempty"`
	ACMGClass   string `yaml:"acmg_class,omitempty"`
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
	Generate bool `yaml:"generate"`
	DeNovo   bool `yaml:"de_novo"`
}

// SequencerSpec describes the sequencing parameters.
type SequencerSpec struct {
	Platform   string `yaml:"platform"`
	Coverage   int    `yaml:"coverage"`
	ReadLength int    `yaml:"read_length"`
}

// OutputSpec specifies output file paths.
type OutputSpec struct {
	ProbandFASTQ string `yaml:"proband_fastq"`
	MotherFASTQ  string `yaml:"mother_fastq,omitempty"`
	FatherFASTQ  string `yaml:"father_fastq,omitempty"`
	GroundTruth  string `yaml:"ground_truth"`
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
	return &sc, nil
}
