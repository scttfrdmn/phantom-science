package cmd

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/spf13/cobra"
)

var generateCmd = &cobra.Command{
	Use:   "generate",
	Short: "Generate synthetic FASTQ reads, optionally streaming to stdout",
	Long: `Generate synthetic FASTQ reads using an HMM model trained with the train subcommand.

Without --stdout, reads are written to the OutputDir specified in --config (default behaviour).
With --stdout, reads are streamed to stdout for direct piping into an aligner:

  phantom-science generate --config seq.yaml --stdout \
    --variants controlled.vcf --sample-id tumor | \
    bwa mem -C ref.fa - | \
    bamq write --stdin s3://bucket/bams/tumor.bam

--variants injects known variants at their allele frequencies (AF= from INFO; default 0.5).
--sample-id appends RG:Z:<id> to every read header; BWA -C copies this into the SAM RG tag.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		configPath, _ := cmd.Flags().GetString("config")
		if configPath == "" {
			return fmt.Errorf("--config is required")
		}
		stdoutMode, _ := cmd.Flags().GetBool("stdout")
		variantsFile, _ := cmd.Flags().GetString("variants")
		sampleID, _ := cmd.Flags().GetString("sample-id")

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

		engine := genomics.NewEngine(*config)
		defer engine.Close()

		if !stdoutMode {
			return engine.Generate()
		}

		var variants []genomics.VariantRecord
		if variantsFile != "" {
			variants, err = loadVariantsVCF(variantsFile)
			if err != nil {
				return fmt.Errorf("loading variants from %s: %w", variantsFile, err)
			}
			fmt.Fprintf(os.Stderr, "   [generate] loaded %d variants from %s\n", len(variants), variantsFile)
		}

		opts := genomics.GenerateOpts{Variants: variants, SampleID: sampleID}
		w := bufio.NewWriterSize(os.Stdout, 1<<20) // 1 MB output buffer
		if err := engine.GenerateToWriter(w, opts); err != nil {
			return err
		}
		return w.Flush()
	},
}

func init() {
	rootCmd.AddCommand(generateCmd)
	generateCmd.Flags().Bool("stdout", false, "Stream FASTQ reads to stdout instead of writing files")
	generateCmd.Flags().String("variants", "", "VCF file of variants to inject (plain or .gz)")
	generateCmd.Flags().String("sample-id", "", "Sample ID; appends RG:Z:<id> to each read header (use with BWA -C)")
	// --config here shadows the persistent root flag so it shows in the subcommand help.
	generateCmd.Flags().String("config", "", "Sequencer configuration file path")
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
