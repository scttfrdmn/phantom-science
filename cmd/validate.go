package cmd

import (
	"fmt"
	"os"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/spf13/cobra"
)

var validateCmd = &cobra.Command{
	Use:   "validate FILE...",
	Short: "Run the validation framework on FASTQ files",
	Long: `Run the validation framework on each FASTQ file and report every check.

The checks come from the validation section of the phantom config (-p) when it
is enabled, and otherwise from the built-in suite: format, quality, complexity,
read length, base composition, GC content, parsing, encoding and header checks.
The command exits non-zero if any file fails.`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		cfg := cli.DefaultValidationTests
		if phantomPath, _ := cmd.InheritedFlags().GetString("phantom"); phantomPath != "" {
			phantom, err := cli.NewConfigLoader(false).LoadPhantomConfig(phantomPath)
			if err != nil {
				return fmt.Errorf("loading phantom config: %w", err)
			}
			if effective, ok := cli.EffectiveValidation(phantom.Validation); ok {
				cfg = effective
			}
		}
		cmd.SilenceUsage = true
		return cli.ValidateFiles(os.Stdout, args, cfg)
	},
}

func init() {
	rootCmd.AddCommand(validateCmd)
}
