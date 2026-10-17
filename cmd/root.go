package cmd

import (
	"fmt"
	"log"
	"os"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/spf13/cobra"
)

const (
	Version   = "0.7.0"
	BuildDate = "2026-10-17"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "phantom-science",
	Short: "Scientific Simulation Platform",
	Long: `Phantom Science is an enterprise-grade scientific simulation platform that generates
realistic synthetic FASTQ files for testing, validation, and development purposes.

Features:
- Realistic sequencer simulation (NovaSeq 6000, NextSeq 2000, MiSeq)
- Enterprise fingerprinting and steganographic marking
- Statistical validation frameworks
- Time compression for accelerated testing
- Interactive Terminal User Interface (TUI)

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
			if err := app.Validate(); err != nil {
				log.Fatalf("Validation failed: %v", err)
			}
		case "verify":
			if err := app.Verify(); err != nil {
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
			fmt.Println("Scientific Simulation Platform")
			fmt.Println("Copyright 2025-2026 Scott Friedman")
			os.Exit(0)
		}
	}
}