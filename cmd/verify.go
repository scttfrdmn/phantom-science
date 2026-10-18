package cmd

import (
	"os"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/spf13/cobra"
)

var verifyCmd = &cobra.Command{
	Use:   "verify FILE...",
	Short: "Check FASTQ files for a steganographic fingerprint",
	Long: `Extract and report the steganographic fingerprint embedded in each FASTQ file.

The command exits non-zero if any file has no fingerprint, cannot be read, or
does not match --customer / --dataset when they are given. It needs a build
with the stegano library (see "make stegano").`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		customer, _ := cmd.Flags().GetString("customer")
		dataset, _ := cmd.Flags().GetString("dataset")
		cmd.SilenceUsage = true
		return cli.VerifyFiles(os.Stdout, args, cli.VerifyExpectation{Customer: customer, Dataset: dataset})
	},
}

func init() {
	rootCmd.AddCommand(verifyCmd)
	verifyCmd.Flags().String("customer", "", "Fail unless every fingerprint names this customer ID")
	verifyCmd.Flags().String("dataset", "", "Fail unless every fingerprint names this dataset ID")
}
