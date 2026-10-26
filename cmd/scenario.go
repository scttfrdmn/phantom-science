//go:build !release

package cmd

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/genomics/clinical"
	runmanifest "github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"github.com/spf13/cobra"
)

var scenarioCmd = &cobra.Command{
	Use:   "scenario",
	Short: "Generate clinical scenario FASTQ with planted variants and ground truth manifest",
	Long: `Generate synthetic FASTQ reads for a clinical demonstration scenario.

The scenario is a YAML file listing HPO terms, which a built-in catalogue maps
to a pathogenic ClinVar variant (GRCh38, VCF-normalised; no network lookup):

  HP:0001250 Seizure + HP:0011972 Hypoglycorrhachia
      SLC2A1 chr1:42929242 C>T p.Gly314Ser missense   (VCV000016113)
  HP:0004923 Hyperphenylalaninemia
      PAH    chr12:102840493 G>A p.Arg408Trp missense (VCV000000577)
  HP:0012236 Elevated sweat chloride
      CFTR   chr7:117559590 ATCT>A p.Phe508del inframe deletion (VCV000007105)

Other terms are an error unless the spec gives phenotype.planted_variant
explicitly. Reads are placed in a 50 kb window around the variant (±25 kb)
at the requested coverage, and reads covering it carry the ALT allele at the
expected VAF (0.5 heterozygous, 1.0 homozygous). A JSON ground-truth
manifest records what was planted.

Read bases. With --reference genome.fa (uncompressed FASTA; a .fai index is
used when present; chromosome names with or without "chr") the planted REF
allele is checked against the FASTA first, and a mismatch is an error. The
reads are then sampled from the reference over the window (clipped to the
contig), and the variant is planted at its real coordinate. Without
--reference the bases are synthetic (uniform, or the configured HMM model):
in a proband-only scenario each read's bases are drawn independently, so
reads do not overlap consistently and cannot be aligned to any reference;
trio scenarios draw one synthetic window sequence.

Read model. Every scenario FASTQ is written by the same read pipeline as
'generate': Illumina reads have binned qualities (Q2/Q12/Q23/Q37) and
quality-matched substitution and indel errors by default (illumina.errors
in --config; "substitutions" leaves out the indels, "none" gives error-free
bases). With sequencer.paired_end: true
the reads are pairs, the two ends of one fragment (illumina.insert_size),
written to output.<member>_fastq and output.<member>_fastq_r2, or
interleaved in output.<member>_fastq when no R2 file is given. sequencer.platform is
recorded in the manifest; the read model is the --config platform
(Illumina by default).

Trio scenarios (trio: {generate: true, de_novo: true} with output.mother_fastq
and output.father_fastq) simulate six phased haplotypes for the window with the
Go structured coalescent (no Python): the parents get two each, the proband
inherits one from each parent with Poisson crossovers, and the planted variant
is a de novo on the proband only. Each member's reads are sampled from that
member's two haplotype sequences, so every read comes from one haplotype.
With --reference the haplotypes are the reference window carrying the
simulated sites: the coalescent gives each site's position and carriers,
the reference base is the ancestral (REF) allele, and the derived (ALT)
allele is the simulated derived base, or the simulated ancestral base when
the derived one equals the reference. Without --reference the window is one
synthetic sequence carrying the simulated alleles. --population selects the
demographic model: YRI, CEU, CHB or MXL, or a 1000 Genomes code modelled by
the closest of these (e.g. GBR, TSI, JPT, LWK); other codes are an error. It
is used only by trio scenarios.

Example, with glut1_scenario.yaml containing the spec below:
  phantom-science scenario --spec glut1_scenario.yaml --reference hg38.fa

The YAML spec format:
  phenotype:
    hpo_terms: [HP:0001250, HP:0011972]
  sequencer:
    platform: NovaSeq6000
    coverage: 30
    read_length: 150
    paired_end: true            # optional
  output:
    proband_fastq: proband_R1.fastq
    proband_fastq_r2: proband_R2.fastq   # with paired_end
    ground_truth: manifest.json

Run manifest and fingerprints. Every scenario also writes phantom-manifest.json
in the directory of output.ground_truth: the sequencer configuration (without
paths or key), the seed, the scenario without its output paths, the input
files by name, size and sha256, and each member FASTQ by size and sha256.
With fingerprinting enabled in --config (fingerprinting.enable, a key from
fingerprinting.key_file or PHANTOM_FINGERPRINT_KEY, enterprise.customer_id
and dataset_id or --dataset-id) or in a phantom config (-p), as for generate,
the scenario's one run seed -- from which the haplotype simulation, the
synthetic window and every member's reads are all derived -- is derived from
the key, customer, dataset and salt (the configured seed, or a drawn value
logged to stderr). The reads are those of an ordinary scenario with that seed;
the manifest records the scheme, customer, dataset and salt, never the seed
or the key. "phantom-science verify --manifest phantom-manifest.json" (with
--reference for a --reference scenario) regenerates the scenario to attribute
suspect reads.`,
	// Run errors (e.g. a REF mismatch) are not usage errors.
	SilenceUsage: true,
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

		if err := applyFingerprintFlags(cmd, sequencerCfg); err != nil {
			return err
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

		gen := clinical.NewScenarioGenerator(sequencerCfg)
		if refPath, _ := cmd.Flags().GetString("reference"); refPath != "" {
			gen.WithReference(refPath)
		}
		manifest, err := gen.Run(cmd.Context(), scenario)
		if err != nil {
			return fmt.Errorf("running scenario: %w", err)
		}

		// Write manifest.
		if err := clinical.WriteManifest(scenario.Output.GroundTruth, manifest); err != nil {
			return fmt.Errorf("writing manifest: %w", err)
		}
		runManifest := filepath.Join(filepath.Dir(scenario.Output.GroundTruth), runmanifest.FileName)
		if err := gen.WriteRunManifest(runManifest, scenario.Output.GroundTruth, Version); err != nil {
			return err
		}

		// Summary to stderr.
		if len(manifest.PlantedVariants) > 0 {
			v := manifest.PlantedVariants[0]
			fmt.Fprintf(os.Stderr, "   [scenario] proband FASTQ: %s\n", scenario.Output.ProbandFASTQ)
			if scenario.Output.ProbandFASTQR2 != "" {
				fmt.Fprintf(os.Stderr, "   [scenario] proband R2 FASTQ: %s\n", scenario.Output.ProbandFASTQR2)
			}
			fmt.Fprintf(os.Stderr, "   [scenario] read bases: %s\n", manifest.Sequencer.Bases)
			fmt.Fprintf(os.Stderr, "   [scenario] planted: %s %s:%d %s>%s (VAF=%.2f)\n",
				v.Gene, v.Chrom, v.Pos, v.Ref, v.Alt, v.ExpectedVAF)
		}
		fmt.Fprintf(os.Stderr, "   [scenario] reads generated: %d  wall: %.3fs\n",
			manifest.ComputeMetrics.ReadsGenerated,
			float64(manifest.ComputeMetrics.WallClockNS)/1e9)
		fmt.Fprintf(os.Stderr, "   [scenario] ground truth: %s\n", scenario.Output.GroundTruth)
		fmt.Fprintf(os.Stderr, "   [scenario] run manifest: %s\n", runManifest)

		return nil
	},
}

func init() {
	rootCmd.AddCommand(scenarioCmd)
	scenarioCmd.Flags().String("spec", "", "Path to YAML scenario specification (required)")
	scenarioCmd.Flags().String("population", "", "Trio population: YRI, CEU, CHB, MXL or a 1000 Genomes alias of them (default CEU)")
	scenarioCmd.Flags().String("dataset-id", "", "Dataset ID (dataset_id); with fingerprinting it names the dataset the keyed seed is derived for")
	scenarioCmd.Flags().String("reference", "", "Reference FASTA (GRCh38 for catalogue variants, uncompressed): reads are sampled from it, and the planted REF allele must match it or the run fails")
}
