package cmd

import (
	"fmt"
	"io"
	"path/filepath"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore/errmodel"
	"github.com/scttfrdmn/phantom-science/internal/genomics/transcriptome"
	"github.com/spf13/cobra"
)

var transcriptomeCmd = &cobra.Command{
	Use:   "transcriptome",
	Short: "Generate ONT cDNA (PCS kit) reads from transcripts with expression and isoform truth",
	Long: `Generate synthetic Oxford Nanopore cDNA-PCR reads (PCS109/110/111/114) from a
transcript table, with known expression and isoform ground truth (#13).

Inputs (one of):
  --transcripts tx.fa              transcript FASTA (GENCODE transcripts.fa or any
                                   FASTA; plain or gzip). No exon structure, so
                                   novel isoforms are alternative 3' ends only.
  --transcripts tx.fa --gtf a.gtf  as above, with exon structure from the GTF for
                                   transcripts whose exon lengths match.
  --genome g.fa --gtf a.gtf        transcripts spliced from the genome; minus-strand
                                   transcripts are reverse-complemented. Use this
                                   to align reads to the genome (minimap2 -ax splice)
                                   and quantify with Bambu.

Expression: by default a Zipf rank-abundance model with log-normal scatter
(log TPM = c - alpha*log(rank) + sigma*Z; with alpha=1 and sigma=0.5 over
GENCODE's ~250,000 transcripts, a few hundred carry half the reads). --expression-profile loads a table instead:
"transcript_id<TAB>tpm" (header optional), a table with a TPM column (salmon
quant.sf), or Bambu counts_transcript.txt (TXNAME GENEID samples...; samples
summed). Values are rescaled to sum to 1e6; IDs match with or without version.

Each read: a molecule drawn by TPM (counts are multinomial), full length with
--full-length-fraction or else 5'-truncated (3' fragment of 100 + geometric
mean 500 bases), a 20-100 nt poly(A), the kit's strand-switch primer and RT
primer, a random strand, and the --tier error model (none: exact constructs).
With --novel-rate p, a share p of each eligible transcript's molecules comes
from its novel isoform (one internal exon skipped when it has >= 3 exons,
otherwise an alternative 3' end).

Outputs in --out DIR: reads.fastq[.gz] (MinKNOW headers), sequencing_summary.txt,
truth.tsv (one row per read), counts.tsv (per transcript and novel isoform:
true TPM, expected and realised reads), novel_isoforms.tsv, and
novel_isoforms.gtf when exon coordinates are known.

--seed makes a run reproducible (byte-identical output); without it a seed is
drawn and logged to stderr.

Example (GRCh38.fa and gencode.gtf are your own downloads, e.g. from GENCODE):
  phantom-science transcriptome --genome GRCh38.fa --gtf gencode.gtf \
      --reads 1000000 --out sim --seed 1 --compress gzip`,
	RunE: func(cmd *cobra.Command, args []string) error {
		return runTranscriptome(cmd, cmd.OutOrStdout(), cmd.ErrOrStderr())
	},
}

func init() {
	rootCmd.AddCommand(transcriptomeCmd)
	// Errors here are about inputs, not flag syntax: print them without usage.
	transcriptomeCmd.SilenceUsage = true
	f := transcriptomeCmd.Flags()
	f.String("transcripts", "", "Transcript FASTA (plain or gzip)")
	f.String("genome", "", "Genome FASTA (with --gtf; plain or gzip)")
	f.String("gtf", "", "GTF annotation (exon records; plain or gzip)")
	f.String("expression-profile", "", "Expression table (id<TAB>tpm, salmon quant.sf or Bambu counts)")
	f.Int("reads", 0, "Number of reads to generate (required)")
	f.String("out", "", "Output directory (required)")
	f.String("tier", transcriptome.DefaultTier, "Basecaller error tier: sup, hac, fast or none")
	f.Uint64("seed", 0, "Random seed (0 or unset: draw one and log it to stderr)")
	f.Float64("novel-rate", transcriptome.DefaultNovelRate, "Share of each eligible transcript's molecules from its novel isoform")
	f.String("compress", "none", "FASTQ compression: none or gzip")
	f.String("kit", transcriptome.DefaultKit, "cDNA kit primers: "+strings.Join(transcriptome.KitNames(), ", "))
	f.Float64("full-length-fraction", transcriptome.DefaultFullLengthFraction, "Fraction of molecules that are full length")
	f.Float64("zipf-alpha", transcriptome.DefaultZipfAlpha, "Default expression model: rank-abundance exponent")
	f.Float64("expression-sigma", transcriptome.DefaultExpressionSigma, "Default expression model: log-normal scatter (natural log)")
	f.String("sample-id", "", "ONT sample_id written to read headers")
}

func runTranscriptome(cmd *cobra.Command, out, logw io.Writer) error {
	f := cmd.Flags()
	txPath, _ := f.GetString("transcripts")
	genomePath, _ := f.GetString("genome")
	gtfPath, _ := f.GetString("gtf")
	profile, _ := f.GetString("expression-profile")
	reads, _ := f.GetInt("reads")
	outDir, _ := f.GetString("out")
	tier, _ := f.GetString("tier")
	seed, _ := f.GetUint64("seed")
	novel, _ := f.GetFloat64("novel-rate")
	compress, _ := f.GetString("compress")
	kitName, _ := f.GetString("kit")
	fullLen, _ := f.GetFloat64("full-length-fraction")
	alpha, _ := f.GetFloat64("zipf-alpha")
	sigma, _ := f.GetFloat64("expression-sigma")
	sampleID, _ := f.GetString("sample-id")

	switch {
	case outDir == "":
		return fmt.Errorf("--out is required")
	case reads <= 0:
		return fmt.Errorf("--reads must be positive")
	case txPath != "" && genomePath != "":
		return fmt.Errorf("use either --transcripts or --genome with --gtf, not both")
	case txPath == "" && genomePath == "":
		return fmt.Errorf("--transcripts or --genome with --gtf is required")
	case genomePath != "" && gtfPath == "":
		return fmt.Errorf("--genome requires --gtf")
	}
	var gz bool
	switch compress {
	case "", "none":
	case "gzip":
		gz = true
	default:
		return fmt.Errorf("unknown --compress %q: want none or gzip", compress)
	}
	if tier != transcriptome.TierNone {
		if _, err := errmodel.ParseTier(tier); err != nil {
			return fmt.Errorf("--tier: %w (or none)", err)
		}
	}
	kit, err := transcriptome.LookupKit(kitName)
	if err != nil {
		return err
	}
	zipf := transcriptome.ZipfModel{Alpha: alpha, Sigma: sigma}
	if err := zipf.Validate(); err != nil {
		return err
	}

	seed, drawn, err := genomics.ResolveSeed(seed)
	if err != nil {
		return err
	}
	if drawn {
		fmt.Fprintf(logw, "   [info] seed %d (pass --seed %d to reproduce this run)\n", seed, seed)
	}

	cfg := transcriptome.DefaultConfig(reads, seed)
	cfg.Tier, cfg.Kit, cfg.NovelRate, cfg.FullLengthFraction = tier, kit, novel, fullLen
	cfg.Run.SampleID = sampleID
	if err := cfg.Validate(); err != nil {
		return err
	}
	if err := cfg.Run.WithDefaults(seed).Validate(); err != nil {
		return err
	}

	table, err := loadTranscriptTable(txPath, genomePath, gtfPath, logw)
	if err != nil {
		return err
	}
	var tpm []float64
	if profile != "" {
		var st transcriptome.ProfileStats
		tpm, st, err = transcriptome.LoadProfile(profile, table)
		if err != nil {
			return err
		}
		fmt.Fprintf(logw, "   [info] expression profile %s: %d rows, %d matched, %d not in the transcript table (values: %s)\n",
			profile, st.Rows, st.Matched, st.Unmatched, st.Column)
	} else {
		tpm = transcriptome.DefaultTPM(table, zipf, seed)
	}

	sim, err := transcriptome.NewSimulator(table, tpm, cfg)
	if err != nil {
		return err
	}
	novelSpecies := 0
	for _, sp := range sim.Species() {
		if sp.Novel != transcriptome.NovelNone {
			novelSpecies++
		}
	}
	res, err := transcriptome.WriteRun(outDir, sim, transcriptome.OutputOptions{Gzip: gz})
	if err != nil {
		return err
	}
	fmt.Fprintf(out, "Wrote %d cDNA reads (%d from %d novel isoforms) from %d transcripts to %s\n",
		res.Reads, res.NovelReads, novelSpecies, len(table.Transcripts), filepath.Clean(outDir))
	fmt.Fprintf(out, "  %s, %s, %s, %s, %s\n", filepath.Base(res.FASTQ), transcriptome.SummaryName,
		transcriptome.TruthName, transcriptome.CountsName, transcriptome.NovelTSVName)
	return nil
}

// loadTranscriptTable loads the transcript table from the CLI inputs,
// reporting skipped records to logw.
func loadTranscriptTable(txPath, genomePath, gtfPath string, logw io.Writer) (*transcriptome.Table, error) {
	var ann *transcriptome.Annotation
	if gtfPath != "" {
		var err error
		if ann, err = transcriptome.ParseGTF(gtfPath); err != nil {
			return nil, err
		}
		if ann.Skipped > 0 {
			fmt.Fprintf(logw, "   [warn] %s: skipped %d unstranded transcripts\n", gtfPath, ann.Skipped)
		}
	}
	if genomePath != "" {
		table, st, err := transcriptome.BuildFromGenome(genomePath, ann)
		if err != nil {
			return nil, err
		}
		if st.MissingChrom > 0 {
			fmt.Fprintf(logw, "   [warn] %d annotated transcripts are on chromosomes missing from %s and were skipped\n", st.MissingChrom, genomePath)
		}
		return table, nil
	}
	table, err := transcriptome.LoadTranscriptFASTA(txPath)
	if err != nil {
		return nil, err
	}
	if ann != nil {
		att, mis := table.AttachAnnotation(ann)
		fmt.Fprintf(logw, "   [info] exon structure from %s for %d of %d transcripts", gtfPath, att, len(table.Transcripts))
		if mis > 0 {
			fmt.Fprintf(logw, " (%d matched by ID but differ in length)", mis)
		}
		fmt.Fprintln(logw)
	}
	return table, nil
}
