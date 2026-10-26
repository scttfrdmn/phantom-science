package cmd

import (
	"os"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/spf13/cobra"
)

var verifyCmd = &cobra.Command{
	Use:   "verify FILE...",
	Short: "Check whether FASTQ or BLOW5 files come from a fingerprinted run",
	Long: `Check whether FASTQ or BLOW5 files come from a fingerprinted run.

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
A run is regenerated with the files it read, which the manifest records by
name, size and sha256, never by path: pass the VCF of a --variants run with
--variants, the FASTA of a --reference run with --reference, the model of an
hmm_model_path run with --hmm-model and the histogram of an empirical
read-length run with --read-length-histogram. A file is used only when its
sha256 matches the record; without it, verify refuses.

Runs of the scenario, signal and transcriptome commands are verified the
same way from the phantom-manifest.json they write (#31). A scenario
manifest (in the directory of its ground_truth file) is regenerated with
--reference for a --reference scenario. A signal manifest needs the run's
input reads (--signal-input), pore model (--pore-model) and, if it used one,
--truth file (--signal-truth); a transcriptome manifest its --transcripts,
--genome, --gtf and --expression-profile files, as recorded. Suspect SLOW5 or
BLOW5 files are matched by raw signal, not read IDs: a read matches when its
whole signal, or half of up to 32 windows of 16 consecutive samples sampled
along it, occur in the regenerated signal, so renamed, re-packed or trimmed
reads still match; a signal verdict always includes a decoy run, because
the chance of a window matching depends on the profile's noise. (Manifests of
phantom-science 0.14.0 and earlier record these paths instead: the file
given with the same name, or the recorded path, is used unchecked, with a
warning.) The key
comes from --key-file or PHANTOM_FINGERPRINT_KEY. Each file gets one
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
		opts.Variants, _ = f.GetStringArray("variants")
		opts.References, _ = f.GetStringArray("reference")
		opts.HMMModels, _ = f.GetStringArray("hmm-model")
		if p, _ := f.GetString("hmm-model-path"); p != "" { // the global flag names one more
			opts.HMMModels = append(opts.HMMModels, p)
		}
		opts.ReadLengths, _ = f.GetStringArray("read-length-histogram")
		opts.SignalInputs, _ = f.GetStringArray("signal-input")
		opts.PoreModels, _ = f.GetStringArray("pore-model")
		opts.SignalTruths, _ = f.GetStringArray("signal-truth")
		opts.Transcripts, _ = f.GetStringArray("transcripts")
		opts.Genomes, _ = f.GetStringArray("genome")
		opts.GTFs, _ = f.GetStringArray("gtf")
		opts.ExpressionProfiles, _ = f.GetStringArray("expression-profile")
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
	f.StringArray("variants", nil, "VCF a candidate run planted variants from (generate --variants); matched to a manifest by sha256 (repeatable)")
	f.StringArray("reference", nil, "Reference FASTA a candidate run sampled reads from (generate or scenario --reference); matched to a manifest by sha256 (repeatable)")
	f.StringArray("hmm-model", nil, "HMM model a candidate run sampled bases from (hmm_model_path); matched to a manifest by sha256 (repeatable)")
	f.StringArray("read-length-histogram", nil, "Read-length histogram a candidate run sampled lengths from (read_length_distribution.empirical_path); matched to a manifest by sha256 (repeatable)")
	f.StringArray("signal-input", nil, "Input reads of a candidate signal run (signal --input); matched to a manifest by sha256 (repeatable)")
	f.StringArray("pore-model", nil, "Pore model of a candidate signal run (signal --pore-model); matched to a manifest by sha256 (repeatable)")
	f.StringArray("signal-truth", nil, "Truth file of a candidate signal run (signal --truth); matched to a manifest by sha256 (repeatable)")
	f.StringArray("transcripts", nil, "Transcript FASTA of a candidate transcriptome run (transcriptome --transcripts); matched to a manifest by sha256 (repeatable)")
	f.StringArray("genome", nil, "Genome FASTA of a candidate transcriptome run (transcriptome --genome); matched to a manifest by sha256 (repeatable)")
	f.StringArray("gtf", nil, "GTF of a candidate transcriptome run (transcriptome --gtf); matched to a manifest by sha256 (repeatable)")
	f.StringArray("expression-profile", nil, "Expression table of a candidate transcriptome run (transcriptome --expression-profile); matched to a manifest by sha256 (repeatable)")
}
