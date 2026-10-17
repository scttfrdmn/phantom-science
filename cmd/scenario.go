//go:build !release

package cmd

import (
	"fmt"
	"os"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/genomics/clinical"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"github.com/spf13/cobra"
)

var scenarioCmd = &cobra.Command{
	Use:   "scenario",
	Short: "Generate clinical scenario FASTQ with planted variants and ground truth manifest",
	Long: `Generate synthetic FASTQ reads for a clinical demonstration scenario.

The scenario is specified in a YAML file that maps HPO terms to a gene and exact
variant coordinates (or an explicit planted_variant). Reads covering the target
locus carry the ALT allele at the expected VAF. A machine-readable ground truth
manifest records what was planted.

Example:
  phantom-science scenario --spec glut1_scenario.yaml --population CEU

The YAML spec format:
  phenotype:
    hpo_terms: [HP:0001250, HP:0001253]
  sequencer:
    platform: NovaSeq6000
    coverage: 30
    read_length: 150
  output:
    proband_fastq: proband.fastq
    ground_truth: manifest.json`,
	RunE: func(cmd *cobra.Command, args []string) error {
		specPath, _ := cmd.Flags().GetString("spec")
		if specPath == "" {
			return fmt.Errorf("--spec is required")
		}
		population, _ := cmd.Flags().GetString("population")
		configPath, _ := cmd.InheritedFlags().GetString("config")

		// Load sequencer config; use minimal defaults when --config not provided.
		var sequencerCfg *types.SequencerConfig
		if configPath != "" {
			loader := cli.NewConfigLoader(false)
			loaded, err := loader.LoadSequencerConfig(configPath)
			if err != nil {
				return fmt.Errorf("loading config: %w", err)
			}
			sequencerCfg = loaded
		} else {
			sequencerCfg = &types.SequencerConfig{
				SampleName:    "proband",
				Lane:          1,
				ReadLength:    150,
				ReadsPerBatch: 1,
				TotalBatches:  1,
				BaseQuality:   30,
			}
		}

		if population != "" {
			sequencerCfg.Population = population
		} else if sequencerCfg.Population == "" {
			sequencerCfg.Population = "CEU"
		}

		scenario, err := clinical.LoadScenario(specPath)
		if err != nil {
			return fmt.Errorf("loading scenario: %w", err)
		}

		projectRoot, _ := os.Getwd()
		gen := clinical.NewScenarioGenerator(sequencerCfg, projectRoot)
		manifest, err := gen.Run(cmd.Context(), scenario)
		if err != nil {
			return fmt.Errorf("running scenario: %w", err)
		}

		// Write manifest.
		if err := clinical.WriteManifest(scenario.Output.GroundTruth, manifest); err != nil {
			return fmt.Errorf("writing manifest: %w", err)
		}

		// Summary to stderr.
		if len(manifest.PlantedVariants) > 0 {
			v := manifest.PlantedVariants[0]
			fmt.Fprintf(os.Stderr, "   [scenario] proband FASTQ: %s\n", scenario.Output.ProbandFASTQ)
			fmt.Fprintf(os.Stderr, "   [scenario] planted: %s %s:%d %s>%s (VAF=%.2f)\n",
				v.Gene, v.Chrom, v.Pos, v.Ref, v.Alt, v.ExpectedVAF)
		}
		fmt.Fprintf(os.Stderr, "   [scenario] reads generated: %d  wall: %.3fs\n",
			manifest.ComputeMetrics.ReadsGenerated,
			float64(manifest.ComputeMetrics.WallClockNS)/1e9)
		fmt.Fprintf(os.Stderr, "   [scenario] ground truth: %s\n", scenario.Output.GroundTruth)

		return nil
	},
}

func init() {
	rootCmd.AddCommand(scenarioCmd)
	scenarioCmd.Flags().String("spec", "", "Path to YAML scenario specification (required)")
	scenarioCmd.Flags().String("population", "", "Override population code (default CEU)")
}
