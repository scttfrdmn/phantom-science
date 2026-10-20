package cmd

import (
	"os"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/spf13/cobra"
)

var verifyCmd = &cobra.Command{
	Use:   "verify FILE...",
	Short: "Check whether FASTQ files come from a fingerprinted run",
	Long: `Check whether FASTQ files come from a fingerprinted run.

Keyed-seed fingerprints (the default) do not change the data: a fingerprinted
run's seed is derived from the key, customer, dataset and salt. verify
regenerates each candidate run from the key and its records, without writing
it, and counts the suspect reads it reproduces: exactly, or by 31-mer
containment, so trimmed, re-oriented or quality-rewritten reads still match.
Candidates come from:

  --manifest phantom-manifest.json   (repeatable; each is one candidate)
  --config CFG --customer C --dataset D --salt S
  --config CFG --candidates candidates.tsv   (customer<TAB>dataset<TAB>salt lines)

--config is the global sequencer-config flag; it must describe the run as
generated (generate flags such as --platform or --reference change it, which
is why the manifest, holding the effective configuration, is preferred).
The key comes from --key-file or PHANTOM_FINGERPRINT_KEY. Each file gets one
verdict: the candidate it matches, with the reads checked and matched and the
matches expected by chance. --decoy also regenerates the first candidate with
a random seed to measure the background (useful for HMM-model or
reference-sampled reads, whose k-mers can recur between runs).

Without --manifest, --config or --candidates, verify extracts the legacy
quality-score fingerprint written by phantom-science 0.9.0 and earlier, which
needs a build with the stegano library (see "make stegano"); --customer and
--dataset then require specific values.

The command exits non-zero if any file matches no candidate, has no legacy
fingerprint, or cannot be read.`,
	Example: `  phantom-science verify --key-file fingerprint.key --manifest out/phantom-manifest.json suspect.fastq.gz
  phantom-science verify --key-file fingerprint.key --config run.json --candidates customers.tsv suspect.fastq`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		f := cmd.Flags()
		opts := cli.VerifyOptions{}
		opts.Customer, _ = f.GetString("customer")
		opts.Dataset, _ = f.GetString("dataset")
		opts.KeyFile, _ = f.GetString("key-file")
		opts.Manifests, _ = f.GetStringArray("manifest")
		opts.ConfigPath, _ = f.GetString("config")
		opts.Candidates, _ = f.GetString("candidates")
		opts.Salt, _ = f.GetUint64("salt")
		opts.SaltSet = f.Changed("salt")
		opts.MaxReads, _ = f.GetInt("max-reads")
		opts.MaxFiles, _ = f.GetInt("max-files")
		opts.MaxRegenReads, _ = f.GetInt64("max-regen-reads")
		opts.Decoy, _ = f.GetBool("decoy")
		cmd.SilenceUsage = true
		return cli.VerifyFiles(os.Stdout, args, opts)
	},
}

func init() {
	rootCmd.AddCommand(verifyCmd)
	f := verifyCmd.Flags()
	f.String("key-file", "", "Fingerprint key (default: $PHANTOM_FINGERPRINT_KEY)")
	f.StringArray("manifest", nil, "Run manifest of a candidate run (repeatable)")
	f.String("customer", "", "Candidate customer ID (with --config); for legacy fingerprints, the customer they must name")
	f.String("dataset", "", "Candidate dataset ID (with --config); for legacy fingerprints, the dataset they must name")
	f.Uint64("salt", 0, "Candidate salt: the run's --seed, or the salt in its manifest or log (with --config)")
	f.String("candidates", "", "File of candidate customer<TAB>dataset<TAB>salt lines (with --config or --manifest)")
	f.Int("max-reads", cli.DefaultVerifyMaxReads, "Reads examined per suspect file")
	f.Int("max-files", 0, "Regenerate at most this many files per candidate (0: all)")
	f.Int64("max-regen-reads", 0, "Regenerate at most this many reads per candidate (0: all)")
	f.Bool("decoy", false, "Also match against the first candidate regenerated with a random seed (measured background)")
}
