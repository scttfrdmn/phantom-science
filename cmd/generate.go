package cmd

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
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

--variants injects known variants at their allele frequencies (AF= from INFO; default 0.5).
Variants with a CHROM/POS are written at their exact coordinates: reads are placed in
--region chr:start-end (1-based, inclusive), or, without --region, in the span of the
variants plus one read length on each side (all variants must then be on one chromosome).
Read bases are generated per read, not copied from a reference genome, so an aligner
will not place the reads at these coordinates in a real reference.
--sample-id appends RG:Z:<id> to every read header; bwa mem -C copies this into the SAM RG tag.

--seed makes a run reproducible: the same seed and config give byte-identical output.
Without --seed a random seed is drawn and logged to stderr. Diagnostics always go to
stderr, so --stdout output is clean FASTQ.

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
	f := generateCmd.Flags()
	f.Bool("stdout", false, "Stream FASTQ reads to stdout instead of writing files")
	f.String("variants", "", "VCF file of variants to inject (plain or .gz); requires --stdout")
	f.String("sample-id", "", "Sample ID; appends RG:Z:<id> to each read header (use with BWA -C)")
	f.String("region", "", "Place reads in chr:start-end (1-based, inclusive) for positioned variants; requires --stdout")
	// --config here shadows the persistent root flag so it shows in the subcommand help.
	f.String("config", "", "Sequencer configuration file path")

	f.Uint64("seed", 0, "Random seed for reproducible output (0 or unset: draw one and log it to stderr)")
	f.String("platform", "", "Sequencing platform: illumina (default) or ont")
	f.String("read-length-dist", "", "Read-length distribution: fixed (default), lognormal, gamma or empirical")
	f.Int("read-length-n50", 0, "Read-length N50 for the lognormal distribution")
	f.Int("read-length-median", 0, "Median read length for the lognormal distribution")
	f.Float64("read-length-mean", 0, "Mean read length")
	f.Float64("read-length-sd", 0, "Read-length standard deviation")
	f.String("basecaller-tier", "", "ONT basecaller tier for the error model: sup, hac or fast (platform ont only)")
	f.String("quality-model", "", "Illumina quality model: binned (default; RTA3 bins Q2/Q12/Q23/Q37) or linear (base_quality - cycle*quality_decay)")
	f.String("illumina-errors", "", "Illumina sequencing errors: quality (default for binned; substitutions at 10^(-Q/10), Q2 as N) or none")
	f.Float64("insert-size-mean", 0, "Mean paired-end insert (fragment) size in bases (default 350)")
	f.Float64("insert-size-sd", 0, "Standard deviation of the paired-end insert size (default 50)")
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
	if !stdoutMode {
		if variantsFile != "" {
			return fmt.Errorf("--variants requires --stdout")
		}
		if regionFlag != "" {
			return fmt.Errorf("--region requires --stdout")
		}
	}

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

	var phantomConfig *types.PhantomConfig
	if phantomPath, _ := cmd.InheritedFlags().GetString("phantom"); phantomPath != "" {
		if phantomConfig, err = loader.LoadPhantomConfig(phantomPath); err != nil {
			return fmt.Errorf("loading phantom config: %w", err)
		}
	}
	config.ResolveFingerprinting(phantomConfig)

	var variants []genomics.VariantRecord
	if variantsFile != "" {
		variants, err = loadVariantsVCF(variantsFile)
		if err != nil {
			return fmt.Errorf("loading variants from %s: %w", variantsFile, err)
		}
		fmt.Fprintf(os.Stderr, "   [generate] loaded %d variants from %s\n", len(variants), variantsFile)
	}
	var region *genomics.GenomicRegion
	if stdoutMode {
		if region, err = resolveStreamRegion(variants, regionFlag, config.ReadLength); err != nil {
			return err
		}
		if region != nil {
			reportUnplacedVariants(variants, region)
		}
	}

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

	opts := genomics.GenerateOpts{Variants: variants, SampleID: sampleID, RefRegion: region}
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

// parseRegion parses chr:start-end (1-based, inclusive; commas allowed) into
// a 0-based half-open region.
func parseRegion(s string) (*genomics.GenomicRegion, error) {
	bad := fmt.Errorf("invalid --region %q: want chr:start-end (1-based, inclusive)", s)
	i := strings.LastIndexByte(s, ':')
	if i <= 0 {
		return nil, bad
	}
	startStr, endStr, ok := strings.Cut(strings.ReplaceAll(s[i+1:], ",", ""), "-")
	if !ok {
		return nil, bad
	}
	start, err1 := strconv.ParseInt(startStr, 10, 64)
	end, err2 := strconv.ParseInt(endStr, 10, 64)
	if err1 != nil || err2 != nil || start < 1 || end < start {
		return nil, bad
	}
	return &genomics.GenomicRegion{Chromosome: s[:i], Start: start - 1, End: end}, nil
}

// resolveStreamRegion returns the region reads are placed in for --stdout:
// the --region flag when given; otherwise, when variants carry positions, the
// span of the positioned variants extended by one read length on each side
// (clamped at 0); otherwise nil (position-blind injection).
func resolveStreamRegion(variants []genomics.VariantRecord, regionFlag string, readLength int) (*genomics.GenomicRegion, error) {
	if regionFlag != "" {
		region, err := parseRegion(regionFlag)
		if err != nil {
			return nil, err
		}
		if region.End-region.Start <= int64(readLength) {
			return nil, fmt.Errorf("--region %s is %d bases; it must be longer than the read length %d", regionFlag, region.End-region.Start, readLength)
		}
		return region, nil
	}
	var region *genomics.GenomicRegion
	chroms := map[string]bool{}
	for _, v := range variants {
		if v.Chromosome == "" || v.Position <= 0 {
			continue
		}
		chroms[v.Chromosome] = true
		start, end := v.Position-1, v.Position-1+int64(len(v.Ref))
		if region == nil {
			region = &genomics.GenomicRegion{Chromosome: v.Chromosome, Start: start, End: end}
			continue
		}
		region.Start = min(region.Start, start)
		region.End = max(region.End, end)
	}
	if len(chroms) > 1 {
		names := make([]string, 0, len(chroms))
		for c := range chroms {
			names = append(names, c)
		}
		sort.Strings(names)
		return nil, fmt.Errorf("variants are on %d chromosomes (%s): pass --region chr:start-end to choose where reads are placed", len(names), strings.Join(names, ", "))
	}
	if region != nil {
		region.Start = max(0, region.Start-int64(readLength))
		region.End += int64(readLength)
	}
	return region, nil
}

// reportUnplacedVariants warns on stderr about variants that region-mode
// injection cannot place: those without a position or on another chromosome.
func reportUnplacedVariants(variants []genomics.VariantRecord, region *genomics.GenomicRegion) {
	unpositioned, elsewhere := 0, 0
	for _, v := range variants {
		switch {
		case v.Chromosome == "" || v.Position <= 0:
			unpositioned++
		case strings.TrimPrefix(v.Chromosome, "chr") != strings.TrimPrefix(region.Chromosome, "chr"):
			elsewhere++
		}
	}
	fmt.Fprintf(os.Stderr, "   [generate] placing reads in %s:%d-%d\n", region.Chromosome, region.Start+1, region.End)
	if unpositioned > 0 {
		fmt.Fprintf(os.Stderr, "   [warn] %d variants without CHROM/POS are ignored in region mode\n", unpositioned)
	}
	if elsewhere > 0 {
		fmt.Fprintf(os.Stderr, "   [warn] %d variants on chromosomes other than %s are ignored\n", elsewhere, region.Chromosome)
	}
}

// loadVariantsVCF reads REF/ALT from a VCF file (plain or .gz) and returns
// VariantRecords. If AF= is present in the INFO field it is used as the
// injection frequency; otherwise 0.5 (heterozygous default) is used.
func loadVariantsVCF(path string) ([]genomics.VariantRecord, error) {
	f, err := os.Open(path) // #nosec G304 -- path from CLI flag
	if err != nil {
		return nil, err
	}
	defer f.Close()

	sc, gzClose, err := openVCFScanner(f, path)
	if err != nil {
		return nil, err
	}
	if gzClose != nil {
		defer gzClose()
	}

	var records []genomics.VariantRecord
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}
		chrom := fields[0]
		var pos int64
		if _, err := fmt.Sscanf(fields[1], "%d", &pos); err != nil {
			pos = 0 // unparseable position → fall back to position-blind
		}
		ref := fields[3]
		alt := fields[4]
		if ref == "." || alt == "." || ref == "" || alt == "" {
			continue
		}
		// Use only the first ALT for multi-allelic sites.
		if idx := strings.IndexByte(alt, ','); idx >= 0 {
			alt = alt[:idx]
		}
		freq := 0.5 // heterozygous default
		if len(fields) >= 8 {
			freq = parseAFFromInfo(fields[7], freq)
		}
		records = append(records, genomics.VariantRecord{
			Ref: ref, Alt: alt, Freq: freq,
			Chromosome: chrom, Position: pos,
		})
	}
	return records, sc.Err()
}

// parseAFFromInfo extracts the AF=<value> field from a VCF INFO column.
// Returns defaultFreq when AF is absent or unparseable.
func parseAFFromInfo(info string, defaultFreq float64) float64 {
	for _, field := range strings.Split(info, ";") {
		if !strings.HasPrefix(field, "AF=") {
			continue
		}
		val := field[3:]
		if idx := strings.IndexByte(val, ','); idx >= 0 {
			val = val[:idx] // first value for multi-allelic
		}
		var f float64
		if _, err := fmt.Sscanf(val, "%f", &f); err == nil && f >= 0 && f <= 1 {
			return f
		}
	}
	return defaultFreq
}
