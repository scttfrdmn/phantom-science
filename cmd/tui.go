package cmd

import (
	"fmt"
	"os"

	"github.com/charmbracelet/bubbletea"
	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/tui"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"github.com/spf13/cobra"
)

// tuiCmd represents the tui command
var tuiCmd = &cobra.Command{
	Use:   "tui",
	Short: "Launch Terminal User Interface for sequencer simulation",
	Long: `Launch an interactive Terminal User Interface (TUI) that runs a timed
sequencer simulation and writes FASTQ files as the run progresses.

What it does:
- Selecting sequencer models (NovaSeq 6000, NextSeq 2000, MiSeq) with the arrow
  keys; the choice replaces sequencer_model and flow_cell_type from --config
- Starting and stopping a run from the Actions panel (Tab moves between panels).
  Files are written on the instrument's simulated schedule at 100x time
  compression, one batch per sample-sheet row, into output_dir
- Fingerprinting settings from --config and --phantom-config (the TUI does
  not run validation)

Not implemented yet: editing the sample sheets (they hold two demo samples),
changing the time compression (the 1x, 10x, 100x and 1000x options are shown
but cannot be selected), the parameter toggles, and live progress (the
progress panel shows the run as complete only when it ends). --enterprise has
no effect, and the global -p flag is not read; use --phantom-config.

Example usage:
  phantom-science tui                                                # MiSeq defaults, ./sequencer_output
  phantom-science tui --config configs/sequencer/v1.0/miseq-v3.yaml  # load a sequencer config
  phantom-science tui --phantom-config configs/phantom/v1.0/enterprise-settings.yaml`,
	Run: runTUI,
}

func init() {
	rootCmd.AddCommand(tuiCmd)
	
	tuiCmd.Flags().String("config", "", "Path to sequencer configuration file (YAML or JSON)")
	tuiCmd.Flags().String("phantom-config", "", "Path to Phantom Science configuration file (YAML or JSON)")
	tuiCmd.Flags().Bool("enterprise", false, "No effect (kept for compatibility)")
}

func runTUI(cmd *cobra.Command, args []string) {
	configPath, _ := cmd.Flags().GetString("config")
	phantomConfigPath, _ := cmd.Flags().GetString("phantom-config")
	_, _ = cmd.Flags().GetBool("enterprise") // enterprise (not used currently)

	// Load configurations
	var sequencerConfig types.SequencerConfig
	var phantomConfig *types.PhantomConfig
	
	if configPath != "" {
		app, err := cli.NewApp(configPath, phantomConfigPath, false)
		if err != nil {
			fmt.Fprintf(os.Stderr, "error: %v\n", err)
			os.Exit(1)
		}
		sequencerConfig = *app.GetConfig()
		phantomConfig = app.GetPhantomConfig()
	} else {
		// Use default configuration for TUI
		sequencerConfig = getDefaultSequencerConfig()
		
		// Set default phantom config
		phantomConfig = &types.PhantomConfig{
			Enterprise: types.EnterpriseConfig{
				CustomerID:   "DEMO",
				ServiceTier:  "standard",
				Fingerprinting: types.FingerprintConfig{
					Enable: false,
				},
			},
			Validation: types.ValidationConfig{
				Enable: false,
			},
		}
		
		if phantomConfigPath != "" {
			loader := cli.NewConfigLoader(false)
			loaded, err := loader.LoadPhantomConfig(phantomConfigPath)
			if err != nil {
				fmt.Fprintf(os.Stderr, "warning: could not load phantom config %s: %v\n",
					phantomConfigPath, err)
			} else {
				phantomConfig = loaded
			}
		}
	}

	// Create TUI model
	model := tui.New(sequencerConfig, phantomConfig)
	
	// Create bubbletea program
	program := tea.NewProgram(
		model,
		tea.WithAltScreen(),       // Use alternate screen buffer
		tea.WithMouseCellMotion(), // Enable mouse support
	)

	// Run the TUI
	if _, err := program.Run(); err != nil {
		fmt.Printf("Error running TUI: %v\n", err)
		os.Exit(1)
	}
}

// getDefaultSequencerConfig returns a default configuration for the TUI
func getDefaultSequencerConfig() types.SequencerConfig {
	return types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		SampleName:     "Demo_Sample",
		Lane:           1,
		OutputDir:      "./sequencer_output",
		ReadLength:     150,
		ReadsPerBatch:  1000,
		BatchInterval:  "5m",
		TotalBatches:   10,
		PairedEnd:      true,
		BaseQuality:    30,
		QualityDecay:   0.1,
		Q30Percentage:  85.0,
	}
}