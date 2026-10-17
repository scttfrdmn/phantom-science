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
	Long: `Launch an interactive Terminal User Interface (TUI) for configuring and running
realistic sequencer simulations. The TUI provides an intuitive interface for:

- Selecting sequencer models (NovaSeq 6000, NextSeq 2000, MiSeq)
- Managing sample sheets with real-time editing
- Configuring run parameters and output settings
- Monitoring simulation progress with realistic timing
- Time compression for faster testing (1x, 10x, 100x, 1000x)

Example usage:
  phantom-science tui                                    # Start with default configuration
  phantom-science tui --config sequencer.yaml          # Load specific sequencer config
  phantom-science tui --phantom-config phantom.yaml    # Load enterprise features`,
	Run: runTUI,
}

func init() {
	rootCmd.AddCommand(tuiCmd)
	
	tuiCmd.Flags().String("config", "", "Path to sequencer configuration file (YAML or JSON)")
	tuiCmd.Flags().String("phantom-config", "", "Path to Phantom Science configuration file (YAML or JSON)")
	tuiCmd.Flags().Bool("enterprise", false, "Enable enterprise features")
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