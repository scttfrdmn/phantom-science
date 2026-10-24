package cmd

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/illumina"
	"github.com/scttfrdmn/phantom-science/internal/storage/s3sink"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

var generateCmd = &cobra.Command{
	Use:   "generate",
	Short: "Generate synthetic FASTQ reads, optionally streaming to stdout",
	Long: `Generate synthetic FASTQ reads from a sequencer config (--config; see configs/).

Bases come from an HMM model when one is configured (--hmm-model-path, hmm_model_path,
or --model-dir with --population; train one with the train subcommand) and are
otherwise uniform random. Without --stdout, reads are written to output_dir from the
config (or --output). With --stdout, reads are streamed to stdout (paired-end mates
interleaved) for piping into another tool:

  phantom-science generate --config configs/sequencer/simple-run.json --stdout \
    --total-reads 1000 --variants controlled.vcf --sample-id tumor > tumor.fastq

(controlled.vcf is your own VCF.)

--reference genome.fa samples every read (paired-end: every fragment) from a reference
genome, uniformly over its A/C/G/T positions and on either strand, so the reads align
back to it. The FASTA must be uncompressed (gunzip or bgzip -d it); a samtools faidx
.fai next to it is used when present. --region chr:start-end (1-based, inclusive) or
--region chr22 restricts sampling to that part, and only that part is loaded: a whole
GRCh38 takes about 3 GB of memory, a 1 Mb region about 1 MB. File output then also
writes a <name>.truth.tsv per batch (read name, contig, 0-based fragment start, strand,
fragment length), listed in the manifest with role "truth". --reference cannot be
combined with an HMM model.

--variants (variants_path) plants known variants at their allele frequencies (AF= from
INFO; default 0.5), in file and --stdout output alike: each read (paired-end: each
fragment, so both mates) draws ALT with probability AF. With --reference, variants with a
CHROM/POS are applied in reference coordinates to every fragment that covers them; a
variant whose REF does not match the reference is skipped and counted in one warning
(--strict-ref makes it an error, before anything is written). Without --reference,
variants with a CHROM/POS are written at their exact coordinates: reads are placed in
--region chr:start-end (1-based, inclusive; variants_region), or, without --region, in the
span of the variants plus one read length on each side (all variants must then be on one
chromosome); read bases are then generated per read, not copied from a reference genome,
so an aligner will not place the reads at these coordinates in a real reference.
File output also writes a <name>.variants.tsv per batch (read name, CHROM, POS, ID, REF,
ALT, the allele the read or fragment carries, its offset), listed in the manifest with
role "truth"; the manifest records the VCF by name, size and sha256, and "verify
--manifest" needs the same VCF (--variants) to regenerate the run.
--sample-id appends RG:Z:<id> to every read header; bwa mem -C copies this into the SAM RG tag.

--seed makes a run reproducible: the same seed and config give byte-identical output.
Without --seed a random seed is drawn and logged to stderr. Diagnostics always go to
stderr, so --stdout output is clean FASTQ.

With fingerprinting enabled (a phantom config, -p, or the config's fingerprinting
section) the seed is derived from the fingerprint key, the customer ID and the
--dataset-id, with --seed (or a drawn, logged value) as the salt; the output is an
ordinary run with that seed, unchanged by the fingerprint. "phantom-science verify
--manifest" attributes files to the run with the key.

Output size and scale (file output):
  --total-bases 140G   write files until the target is met (within one read of it)
  --total-reads N      exactly N reads (read pairs for paired-end); last file partial
  --flow-cells 2       independent runs in per-flow-cell subdirectories, seed-derived IDs
  --workers 8          generate 8 files in parallel; output is identical for any count
  --output s3://bucket/prefix/   spool, finalize and multipart-upload each file
max_output_gb / max_output_tb in the config cap the bases of each flow cell. Every run
writes phantom-manifest.json (seed, sanitized config, size and sha256 of every output).`,
	RunE: func(cmd *cobra.Command, args []string) error {
		return runGenerate(cmd, os.Stdout)
	},
}

func init() {
	rootCmd.AddCommand(generateCmd)
	generateCmd.SilenceUsage = true
	f := generateCmd.Flags()
	f.Bool("stdout", false, "Stream FASTQ reads to stdout instead of writing files")
	f.String("variants", "", "VCF file of variants to plant at their AF (plain or .gz; variants_path)")
	f.String("sample-id", "", "Sample ID; appends RG:Z:<id> to each read header (use with BWA -C)")
	f.String("region", "", "With --reference: sample reads only from chr:start-end (1-based, inclusive) or a whole contig. Without: place reads in chr:start-end for positioned variants (variants_region)")
	f.String("reference", "", "Sample reads from this uncompressed reference FASTA (reference_path)")
	f.Bool("strict-ref", false, "With --reference, fail when a --variants REF does not match the reference (default: skip and warn)")
	// --config here shadows the persistent root flag so it shows in the subcommand help.
	f.String("config", "", "Sequencer configuration file path")

	f.Uint64("seed", 0, "Random seed for reproducible output (0 or unset: draw one and log it to stderr)")
	f.String("dataset-id", "", "Dataset ID (dataset_id); with fingerprinting it names the dataset the keyed seed is derived for")
	f.String("platform", "", "Sequencing platform: illumina (default) or ont")
	f.String("read-length-dist", "", "Read-length distribution: fixed (default), lognormal, gamma or empirical")
	f.Int("read-length-n50", 0, "Read-length N50 for the lognormal distribution")
	f.Int("read-length-median", 0, "Median read length for the lognormal distribution")
	f.Float64("read-length-mean", 0, "Mean read length")
	f.Float64("read-length-sd", 0, "Read-length standard deviation")
	f.String("basecaller-tier", "", "ONT basecaller tier for the error model: sup, hac or fast (platform ont only)")
	f.String("quality-model", "", "Illumina quality model: binned (default; NovaSeq 6000 bins Q2/Q11/Q25/Q37, fitted to a real run, see --illumina-profile) or linear (base_quality - cycle*quality_decay)")
	f.String("illumina-errors", "", "Illumina sequencing errors: quality (default for binned; substitutions at 10^(-Q/10), Q2 as N) or none")
	f.String("illumina-profile", "", "Illumina run profile of the binned model: "+strings.Join(illumina.ProfileNames(), ", ")+" (default "+illumina.DefaultProfile+"; set explicitly, it also brings its run's insert sizes)")
	f.Float64("insert-size-mean", 0, "Mean paired-end insert (fragment) size in bases (default: 455, lognormal; with --illumina-profile, the profile's measured insert sizes)")
	f.Float64("insert-size-sd", 0, "Standard deviation of the paired-end insert size (default: 100; with --illumina-profile, the profile's)")
	f.String("compress", "", "Output compression: none (default) or gzip")
	addScaleFlags(f)
	f.Int64("total-reads", 0, "Target number of reads (read pairs for paired-end) per flow cell, replacing total_batches")
	f.Int("workers", 0, "Output files generated in parallel (default 1); output is identical for any value")
	f.String("output", "", "Output directory or s3://bucket/prefix/ (overrides output_dir in the config)")
}

// addScaleFlags adds the output-size and parallelism flags shared by
// generate and bench.
func addScaleFlags(f *pflag.FlagSet) {
	f.String("total-bases", "", "Target output bases per flow cell, e.g. 140G (K/M/G/T suffixes); output stops within one read of it")
	f.Int("flow-cells", 0, "Number of independent flow-cell runs, each in its own subdirectory with seed-derived IDs")
}

// runGenerate implements the generate command, streaming to out with --stdout.
func runGenerate(cmd *cobra.Command, out io.Writer) error {
	flags := cmd.Flags()
	configPath, _ := flags.GetString("config")
	if configPath == "" {
		return fmt.Errorf("--config is required")
	}
	stdoutMode, _ := flags.GetBool("stdout")
	variantsFile, _ := flags.GetString("variants")
	sampleID, _ := flags.GetString("sample-id")
	regionFlag, _ := flags.GetString("region")

	// Persistent flags inherited from rootCmd.
	hmmModelPath, _ := cmd.InheritedFlags().GetString("hmm-model-path")
	population, _ := cmd.InheritedFlags().GetString("population")
	modelDir, _ := cmd.InheritedFlags().GetString("model-dir")

	loader := cli.NewConfigLoader(false)
	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		return fmt.Errorf("loading config: %w", err)
	}

	if hmmModelPath != "" {
		config.HMMModelPath = hmmModelPath
	}
	if population != "" {
		config.Population = population
	}
	if modelDir != "" {
		config.ModelDir = modelDir
	}
	if err := applyGenerateFlags(cmd, config); err != nil {
		return err
	}
	if variantsFile != "" {
		config.VariantsPath = variantsFile
	}
	// With a reference, --region restricts sampling; without one it places
	// reads for positioned variants.
	if regionFlag != "" {
		switch {
		case config.ReferencePath != "":
			config.ReferenceRegion = regionFlag
		case !stdoutMode && config.VariantsPath == "":
			return fmt.Errorf("--region without --reference places reads for positioned --variants: pass --variants, or --reference to sample from a region")
		default:
			config.VariantsRegion = regionFlag
		}
	}

	var phantomConfig *types.PhantomConfig
	if phantomPath, _ := cmd.InheritedFlags().GetString("phantom"); phantomPath != "" {
		if phantomConfig, err = loader.LoadPhantomConfig(phantomPath); err != nil {
			return fmt.Errorf("loading phantom config: %w", err)
		}
	}
	config.ResolveFingerprinting(phantomConfig)

	// The engine loads the variants (variants_path) and resolves where reads
	// are placed; a bad VCF, region or --strict-ref mismatch fails before any
	// output is written.
	engine := genomics.NewEngine(*config)
	defer engine.Close()
	engine.SetToolVersion(Version)

	if !stdoutMode {
		if err := engine.Generate(); err != nil {
			return err
		}
		if s3sink.IsURI(config.OutputDir) {
			if phantomConfig != nil {
				if _, ok := cli.EffectiveValidation(phantomConfig.Validation); ok {
					fmt.Fprintf(os.Stderr, "   [warn] validation skipped: output was uploaded to %s\n", config.OutputDir)
				}
			}
			return nil
		}
		// Validate the generated files when the phantom config enables validation.
		if phantomConfig != nil {
			if vcfg, ok := cli.EffectiveValidation(phantomConfig.Validation); ok {
				lengths, err := cli.ReadLengthValidation(*config)
				if err != nil {
					return err
				}
				return cli.ValidateFiles(out, engine.GeneratedFiles(), vcfg, lengths)
			}
		}
		return nil
	}

	opts := genomics.GenerateOpts{SampleID: sampleID}
	w := bufio.NewWriterSize(out, 1<<20) // 1 MB output buffer
	if err := engine.GenerateToWriter(w, opts); err != nil {
		return err
	}
	return w.Flush()
}

// applyGenerateFlags copies explicitly set generation flags into config.
func applyGenerateFlags(cmd *cobra.Command, config *types.SequencerConfig) error {
	flags := cmd.Flags()
	set := func(name string) bool { return flags.Changed(name) }
	if set("seed") {
		config.Seed, _ = flags.GetUint64("seed")
	}
	if set("dataset-id") {
		config.DatasetID, _ = flags.GetString("dataset-id")
	}
	if set("platform") {
		config.Platform, _ = flags.GetString("platform")
	}
	if set("read-length-dist") {
		config.ReadLengthDistribution.Dist, _ = flags.GetString("read-length-dist")
	}
	if set("read-length-n50") {
		config.ReadLengthDistribution.N50, _ = flags.GetInt("read-length-n50")
	}
	if set("read-length-median") {
		config.ReadLengthDistribution.Median, _ = flags.GetInt("read-length-median")
	}
	if set("read-length-mean") {
		config.ReadLengthDistribution.Mean, _ = flags.GetFloat64("read-length-mean")
	}
	if set("read-length-sd") {
		config.ReadLengthDistribution.SD, _ = flags.GetFloat64("read-length-sd")
	}
	if set("basecaller-tier") {
		config.ErrorModel.Tier, _ = flags.GetString("basecaller-tier")
	}
	if set("quality-model") {
		config.Illumina.QualityModel, _ = flags.GetString("quality-model")
	}
	if set("illumina-profile") {
		config.Illumina.Profile, _ = flags.GetString("illumina-profile")
	}
	if set("illumina-errors") {
		config.Illumina.Errors, _ = flags.GetString("illumina-errors")
	}
	if set("insert-size-mean") {
		config.Illumina.InsertSize.Mean, _ = flags.GetFloat64("insert-size-mean")
	}
	if set("insert-size-sd") {
		config.Illumina.InsertSize.SD, _ = flags.GetFloat64("insert-size-sd")
	}
	if set("compress") {
		config.Output.Compression, _ = flags.GetString("compress")
	}
	if set("total-bases") {
		s, _ := flags.GetString("total-bases")
		n, err := parseBaseCount(s)
		if err != nil {
			return fmt.Errorf("--total-bases: %w", err)
		}
		config.TotalBases = n
	}
	if set("total-reads") {
		config.TotalReads, _ = flags.GetInt64("total-reads")
	}
	if set("flow-cells") {
		config.FlowCells, _ = flags.GetInt("flow-cells")
	}
	if set("workers") {
		if n, err := flags.GetInt("workers"); err == nil {
			config.Workers = n
		}
	}
	if set("output") {
		config.OutputDir, _ = flags.GetString("output")
	}
	if set("reference") {
		config.ReferencePath, _ = flags.GetString("reference")
	}
	if set("strict-ref") {
		config.StrictRef, _ = flags.GetBool("strict-ref")
	}
	return nil
}

// parseBaseCount parses a base count with an optional K, M, G or T suffix
// (powers of 1000), e.g. "140G".
func parseBaseCount(s string) (int64, error) {
	t := strings.ToUpper(strings.TrimSpace(s))
	mult := int64(1)
	if n := len(t); n > 0 {
		switch t[n-1] {
		case 'K':
			mult = 1e3
		case 'M':
			mult = 1e6
		case 'G':
			mult = 1e9
		case 'T':
			mult = 1e12
		}
		if mult > 1 {
			t = t[:n-1]
		}
	}
	v, err := strconv.ParseFloat(t, 64)
	if err != nil || v <= 0 || v*float64(mult) > 9e18 {
		return 0, fmt.Errorf("invalid base count %q: want a positive number with an optional K, M, G or T suffix", s)
	}
	return int64(v * float64(mult)), nil
}
