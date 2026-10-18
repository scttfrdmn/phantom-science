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

Encrypted fingerprints need the key they were embedded with, from --key-file
or the PHANTOM_FINGERPRINT_KEY environment variable. Each line reports whether
the fingerprint was encrypted.

The command exits non-zero if any file has no fingerprint, cannot be read or
decrypted, or does not match --customer / --dataset when they are given. It
needs a build with the stegano library (see "make stegano").`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		customer, _ := cmd.Flags().GetString("customer")
		dataset, _ := cmd.Flags().GetString("dataset")
		keyFile, _ := cmd.Flags().GetString("key-file")
		cmd.SilenceUsage = true
		return cli.VerifyFiles(os.Stdout, args, cli.VerifyOptions{Customer: customer, Dataset: dataset, KeyFile: keyFile})
	},
}

func init() {
	rootCmd.AddCommand(verifyCmd)
	verifyCmd.Flags().String("customer", "", "Fail unless every fingerprint names this customer ID")
	verifyCmd.Flags().String("dataset", "", "Fail unless every fingerprint names this dataset ID")
	verifyCmd.Flags().String("key-file", "", "Key for encrypted fingerprints (default: $PHANTOM_FINGERPRINT_KEY)")
}
