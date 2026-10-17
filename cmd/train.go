package cmd

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/spf13/cobra"
)

var trainCmd = &cobra.Command{
	Use:   "train",
	Short: "Train an HMM model from VCF data and export to .hmm.json",
	Long: `Train a Hidden Markov Model (HMM) from VCF variant call data and save it as a
.hmm.json file for use with --hmm-model-path.

Model types:
  adaptive  Stream VCF sequences through an AdaptiveHMMModel in 100 KB chunks.
            Default — well-suited for large files and population-specific training.
  kmer      Accumulate all VCF sequences before training an HMMKmerModel.
            Best for small to medium VCF datasets.

Example usage:
  phantom-science train --vcf variants.vcf.gz --out YRI.hmm.json
  phantom-science train --vcf YRI.vcf.gz --vcf YRI_chr2.vcf.gz --out YRI.hmm.json
  phantom-science train --vcf variants.vcf --out model.hmm.json --model-type kmer --kmer-size 6
  phantom-science train --vcf variants.vcf.gz --out models/CEU.hmm.json --population CEU`,
	RunE: func(cmd *cobra.Command, args []string) error {
		vcfFiles, _ := cmd.Flags().GetStringArray("vcf")
		outPath, _ := cmd.Flags().GetString("out")
		modelType, _ := cmd.Flags().GetString("model-type")
		kmerSize, _ := cmd.Flags().GetInt("kmer-size")
		population, _ := cmd.Flags().GetString("population")
		return runTrain(vcfFiles, outPath, modelType, population, kmerSize)
	},
}

func init() {
	rootCmd.AddCommand(trainCmd)
	trainCmd.Flags().StringArray("vcf", nil, "VCF input file (.vcf or .vcf.gz); may be repeated")
	trainCmd.Flags().String("out", "", "Output path for the trained model (.hmm.json)")
	trainCmd.Flags().String("model-type", "adaptive", `Model type to train: "adaptive" or "kmer"`)
	trainCmd.Flags().Int("kmer-size", 4, "K-mer size for the model (2–8)")
	_ = trainCmd.MarkFlagRequired("out")
}

// runTrain trains an HMM model from one or more VCF files and writes it to outPath.
func runTrain(vcfFiles []string, outPath, modelType, population string, kmerSize int) error {
	if len(vcfFiles) == 0 {
		return fmt.Errorf("at least one --vcf file is required")
	}
	if kmerSize < 2 || kmerSize > 8 {
		return fmt.Errorf("--kmer-size must be between 2 and 8, got %d", kmerSize)
	}

	fmt.Printf("   [train] model-type=%s  k=%d  population=%q\n", modelType, kmerSize, population)
	fmt.Printf("   [train] VCF inputs: %v\n", vcfFiles)
	fmt.Printf("   [train] output: %s\n", outPath)

	switch modelType {
	case "adaptive":
		return trainAdaptive(vcfFiles, outPath, kmerSize)
	case "kmer":
		return trainKmer(vcfFiles, outPath, kmerSize)
	default:
		return fmt.Errorf("unknown --model-type %q: choose \"adaptive\" or \"kmer\"", modelType)
	}
}

// trainAdaptive streams VCF sequences through an AdaptiveHMMModel in 100 KB chunks.
func trainAdaptive(vcfFiles []string, outPath string, kmerSize int) error {
	model := authenticity.NewAdaptiveHMMModel(kmerSize)
	if model == nil {
		return fmt.Errorf("failed to create adaptive HMM model with k=%d", kmerSize)
	}

	for _, path := range vcfFiles {
		fmt.Printf("   [train] reading %s...\n", path)
		if err := streamVCFToAdaptive(path, model); err != nil {
			return fmt.Errorf("training from %s: %w", path, err)
		}
	}

	if !model.IsTrained() {
		return fmt.Errorf("no sequences were trained; check that VCF files contain variant records")
	}

	if err := model.Save(outPath); err != nil {
		return fmt.Errorf("saving model to %s: %w", outPath, err)
	}
	fmt.Printf("   [train] adaptive HMM model saved to %s\n", outPath)
	return nil
}

// streamVCFToAdaptive reads REF+ALT from a VCF file and trains the model in 100 KB chunks.
func streamVCFToAdaptive(path string, model *authenticity.AdaptiveHMMModel) error {
	f, err := os.Open(path) // #nosec G304 -- path from CLI flag
	if err != nil {
		return err
	}
	defer f.Close()

	scanner, gzClose, err := openVCFScanner(f, path)
	if err != nil {
		return err
	}
	if gzClose != nil {
		defer gzClose()
	}

	var buf strings.Builder
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.SplitN(line, "\t", 6)
		if len(fields) < 5 {
			continue
		}
		appendVCFAlleles(&buf, fields[3], fields[4])
		if buf.Len() >= 100_000 {
			if err := model.TrainFromSequence([]byte(buf.String())); err != nil {
				return err
			}
			buf.Reset()
		}
	}
	if buf.Len() > 0 {
		if err := model.TrainFromSequence([]byte(buf.String())); err != nil {
			return err
		}
	}
	return scanner.Err()
}

// trainKmer accumulates all VCF sequences and trains an HMMKmerModel via TrainFromReferenceGenome.
func trainKmer(vcfFiles []string, outPath string, kmerSize int) error {
	model := authenticity.NewHMMKmerModel(kmerSize)
	if model == nil {
		return fmt.Errorf("failed to create kmer HMM model with k=%d", kmerSize)
	}

	var allSeq strings.Builder
	for _, path := range vcfFiles {
		fmt.Printf("   [train] reading %s...\n", path)
		if err := accumulateVCFSequences(path, &allSeq); err != nil {
			return fmt.Errorf("reading %s: %w", path, err)
		}
	}

	seq := allSeq.String()
	if len(seq) < kmerSize {
		return fmt.Errorf("accumulated sequence too short (%d bases) to train a k=%d model", len(seq), kmerSize)
	}

	genome := &vcfGenome{
		chromosomes: map[string]*vcfChromosome{
			"chr1": {name: "chr1", sequence: []byte(seq)},
		},
		names: []string{"chr1"},
	}
	if err := model.TrainFromReferenceGenome(genome); err != nil {
		return fmt.Errorf("training kmer model: %w", err)
	}

	if err := model.Save(outPath); err != nil {
		return fmt.Errorf("saving model to %s: %w", outPath, err)
	}
	fmt.Printf("   [train] kmer HMM model saved to %s\n", outPath)
	return nil
}

// accumulateVCFSequences appends REF+ALT sequences from a VCF file to buf.
func accumulateVCFSequences(path string, buf *strings.Builder) error {
	f, err := os.Open(path) // #nosec G304
	if err != nil {
		return err
	}
	defer f.Close()

	scanner, gzClose, err := openVCFScanner(f, path)
	if err != nil {
		return err
	}
	if gzClose != nil {
		defer gzClose()
	}

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.SplitN(line, "\t", 6)
		if len(fields) < 5 {
			continue
		}
		appendVCFAlleles(buf, fields[3], fields[4])
	}
	return scanner.Err()
}

// appendVCFAlleles appends REF and the first ALT allele to buf.
// Multi-allelic ALT (comma-separated) uses only the first allele.
func appendVCFAlleles(buf *strings.Builder, ref, alt string) {
	if len(ref) > 0 {
		buf.WriteString(ref)
	}
	if len(alt) > 0 && alt != "." {
		if idx := strings.IndexByte(alt, ','); idx >= 0 {
			alt = alt[:idx]
		}
		buf.WriteString(alt)
	}
}

// openVCFScanner returns a *bufio.Scanner for a VCF file (plain or .gz).
// If the file is gzip-compressed, the returned cleanup func closes the gzip reader
// and must be deferred by the caller.
func openVCFScanner(f *os.File, path string) (*bufio.Scanner, func() error, error) {
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, nil, fmt.Errorf("opening gzip reader for %s: %w", path, err)
		}
		sc := bufio.NewScanner(gz)
		sc.Buffer(make([]byte, 0, 64*1024), 1024*1024)
		return sc, gz.Close, nil
	}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	return sc, nil, nil
}

// vcfChromosome implements authenticity.ChromosomeInterface for in-memory VCF sequences.
type vcfChromosome struct {
	name     string
	sequence []byte
}

func (c *vcfChromosome) GetName() string     { return c.name }
func (c *vcfChromosome) GetSequence() []byte { return c.sequence }
func (c *vcfChromosome) GetLength() int      { return len(c.sequence) }

// vcfGenome implements authenticity.GenomeInterface backed by in-memory sequences.
type vcfGenome struct {
	chromosomes map[string]*vcfChromosome
	names       []string
}

func (g *vcfGenome) GetChromosomeNames() []string { return g.names }
func (g *vcfGenome) LoadChromosome(name string) (authenticity.ChromosomeInterface, error) {
	chr, ok := g.chromosomes[name]
	if !ok {
		return nil, fmt.Errorf("chromosome %q not found in VCF genome", name)
	}
	return chr, nil
}
