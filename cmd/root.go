package cmd

import (
	"fmt"
	"log"
	"os"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/spf13/cobra"
)

const (
	Version   = "0.11.0"
	BuildDate = "2026-10-20"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "phantom-science",
	Short: "Synthetic sequencing data generator",
	Long: `Phantom Science generates synthetic sequencing data (FASTQ files, ONT cDNA
reads and raw nanopore signal) for testing sequencer pipelines and tools.

Features:
- Illumina-style FASTQ (single or paired end) and Oxford Nanopore FASTQ with
  MinKNOW headers, a sequencing summary, basecaller error tiers and read-length
  distributions (generate)
- Bases from an HMM model trained on VCF data (train), or uniform random bases
  when no model is given
- Byte-reproducible output from --seed, and a phantom-manifest.json for each
  run written to files
- Keyed-seed fingerprints that leave the generated data unchanged, verified
  by regenerating a run from its manifest and the key (verify)
- Format, quality and composition checks on FASTQ files (validate)
- ONT cDNA reads with expression and isoform ground truth (transcriptome) and
  raw signal as BLOW5 (signal)
- Clinical scenarios with planted variants (scenario), throughput benchmarks
  (bench), and an interactive terminal interface for timed sequencer runs (tui)

Running phantom-science with --config and no subcommand runs the legacy
--command (generate, validate or verify).

Copyright 2025-2026 Scott Friedman`,
	Run: func(cmd *cobra.Command, args []string) {
		// If no subcommand is provided, run the legacy generate command
		config, _ := cmd.Flags().GetString("config")
		phantom, _ := cmd.Flags().GetString("phantom")
		verbose, _ := cmd.Flags().GetBool("verbose")
		command, _ := cmd.Flags().GetString("command")
		hmmModelPath, _ := cmd.Flags().GetString("hmm-model-path")
		population, _ := cmd.Flags().GetString("population")
		modelDir, _ := cmd.Flags().GetString("model-dir")

		if config == "" {
			log.Fatal("Sequencer configuration file is required. Use --config flag or try 'phantom-science tui' for interactive mode.")
		}

		if verbose {
			log.SetFlags(log.LstdFlags | log.Lshortfile)
		}

		app, err := cli.NewApp(config, phantom, verbose)
		if err != nil {
			log.Fatalf("Configuration error: %v", err)
		}

		// Flag overrides config-file values; empty flags leave config-file values intact.
		if hmmModelPath != "" || population != "" || modelDir != "" {
			cfg := *app.GetConfig()
			if hmmModelPath != "" {
				cfg.HMMModelPath = hmmModelPath
			}
			if population != "" {
				cfg.Population = population
			}
			if modelDir != "" {
				cfg.ModelDir = modelDir
			}
			app.SetConfig(cfg)
		}

		switch command {
		case "generate":
			if err := app.Generate(); err != nil {
				log.Fatalf("Generation failed: %v", err)
			}
		case "validate":
			if err := app.Validate(args); err != nil {
				log.Fatalf("Validation failed: %v", err)
			}
		case "verify":
			if err := app.Verify(args, cli.VerifyOptions{}); err != nil {
				log.Fatalf("Verification failed: %v", err)
			}
		default:
			log.Fatalf("Unknown command: %s", command)
		}
	},
}

// Execute adds all child commands to the root command and sets flags appropriately.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}

func init() {
	// Execute prints the returned error; without this cobra prints it as well.
	rootCmd.SilenceErrors = true

	// Global flags
	rootCmd.PersistentFlags().StringP("config", "c", "", "Sequencer configuration file path")
	rootCmd.PersistentFlags().StringP("phantom", "p", "", "Phantom Science configuration file path (optional)")
	rootCmd.PersistentFlags().BoolP("verbose", "v", false, "Enable verbose logging")
	rootCmd.PersistentFlags().String("hmm-model-path", "", "Path to serialized HMM model (.hmm.json) — overrides hmm_model_path in config file")
	rootCmd.PersistentFlags().String("population", "", "1000 Genomes population code (e.g. YRI, CEU, CHB) — used with --model-dir")
	rootCmd.PersistentFlags().String("model-dir", "", "Directory containing <population>.hmm.json model files — used with --population")

	// Legacy flags for backward compatibility
	rootCmd.Flags().String("command", "generate", "Command to execute (generate, validate, verify)")
	rootCmd.Flags().Bool("version", false, "Print version information")

	// Handle version flag
	rootCmd.PreRun = func(cmd *cobra.Command, args []string) {
		version, _ := cmd.Flags().GetBool("version")
		if version {
			fmt.Printf("Phantom Science v%s (built: %s)\n", Version, BuildDate)
			fmt.Println("Synthetic sequencing data generator")
			fmt.Println("Copyright 2025-2026 Scott Friedman")
			os.Exit(0)
		}
	}
}
