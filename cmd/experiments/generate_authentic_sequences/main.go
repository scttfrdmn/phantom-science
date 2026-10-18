// Command generate_authentic_sequences writes 15 FASTA sequences (1, 5 and
// 10 kb) for each of six population labels under
// sequences/authentic_vcf_trained/, plus a Markdown summary. Despite the
// name, every sequence comes from a fresh, untrained k=8 adaptive HMM: no
// trained model is loaded, so the populations do not differ. The variant
// counts, training times, GC and complexity figures in the headers and
// report are hard-coded from an earlier log.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// AuthenticVCFPopulation holds a population label and training figures
// copied from an earlier log; nothing here is loaded from a trained model.
type AuthenticVCFPopulation struct {
	Population        string
	Region            string
	ProcessedVariants int64
	TrainingTime      time.Duration
	ActualGCContent   float64
	ActualComplexity  float64
}

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== SEQUENCE GENERATION (UNTRAINED MODELS) ===")
	fmt.Println("WARNING: each model is a fresh, untrained k=8 adaptive HMM; no trained model is loaded")
	fmt.Println("Population labels and training statistics are hard-coded from an earlier log")
	fmt.Println()

	// Only the 6 populations that successfully trained (from log analysis)
	authenticPopulations := getAuthenticVCFPopulations()

	fmt.Printf("🧬 Populations: %d\n", len(authenticPopulations))
	fmt.Printf("📊 Variants processed in the earlier run (hard-coded): 484,806,516\n")
	fmt.Printf("⏱️ Training time of the earlier run (hard-coded): 8h44m33s\n")
	fmt.Println()

	// Create output directory
	outputDir := "sequences/authentic_vcf_trained"
	err := os.MkdirAll(outputDir, 0750)
	if err != nil {
		log.Fatalf("Failed to create output directory: %v", err)
	}

	fmt.Printf("📁 Output directory: %s\n", outputDir)
	fmt.Println()

	totalSequences := 0
	successfulSequences := 0

	for _, pop := range authenticPopulations {
		fmt.Printf("🧬 %s (%s): Generating sequences...\n",
			pop.Population, pop.Region)

		// Create population directory
		popDir := filepath.Join(outputDir, pop.Population)
		err := os.MkdirAll(popDir, 0750)
		if err != nil {
			log.Printf("Warning: Failed to create directory for %s: %v", pop.Population, err)
			continue
		}

		// Create k=8 model for sequence generation
		model := authenticity.NewAdaptiveHMMModel(8)

		// The model is never trained or loaded, so every population generates
		// from the same untrained model.

		lengths := []int{1000, 5000, 10000}
		sequencesPerLength := 5

		popSuccessful := 0

		for _, length := range lengths {
			for i := 0; i < sequencesPerLength; i++ {
				totalSequences++

				// Generate sequence from the untrained model
				result, err := model.GenerateSequence(length)
				if err != nil {
					log.Printf("Warning: Failed to generate sequence for %s: %v", pop.Population, err)
					continue
				}

				// Create filename
				filename := fmt.Sprintf("%s_AUTHENTIC_VCF_L%d_S%02d.fasta",
					pop.Population, length, i+1)
				filepath := filepath.Join(popDir, filename)

				// Write FASTA file
				err = writeAuthenticVCFFASTA(filepath, pop, length, i+1, result.Sequence)
				if err != nil {
					log.Printf("Warning: Failed to write sequence %s: %v", filename, err)
					continue
				}

				successfulSequences++
				popSuccessful++
			}
		}

		fmt.Printf("    ✅ Generated %d sequences\n", popSuccessful)
	}

	fmt.Printf("\n📊 SEQUENCE GENERATION SUMMARY:\n")
	fmt.Printf("✅ Successfully generated: %d/%d sequences\n", successfulSequences, totalSequences)
	fmt.Printf("🧬 Populations: 6 of 11 labels\n")
	fmt.Printf("📁 Output location: %s\n", outputDir)

	// Generate summary report
	generateAuthenticSummaryReport(outputDir, authenticPopulations, successfulSequences)

	fmt.Printf("\n🏁 SEQUENCE GENERATION COMPLETED\n")
	fmt.Printf("🎯 These sequences come from untrained models, not VCF-trained ones\n")
}

func getAuthenticVCFPopulations() []AuthenticVCFPopulation {
	// Figures copied from an earlier training log; not recomputed or verified
	return []AuthenticVCFPopulation{
		{
			Population:        "YRI",
			Region:            "Africa",
			ProcessedVariants: 80801086,
			TrainingTime:      parseTrainingTime("1h29m14.770624708s"),
			ActualGCContent:   0.4982,
			ActualComplexity:  0.9430,
		},
		{
			Population:        "LWK",
			Region:            "Africa",
			ProcessedVariants: 80801086,
			TrainingTime:      parseTrainingTime("1h26m56.362908209s"),
			ActualGCContent:   0.4982,
			ActualComplexity:  0.9430,
		},
		{
			Population:        "GWD",
			Region:            "Africa",
			ProcessedVariants: 80801086,
			TrainingTime:      parseTrainingTime("1h26m8.869497333s"),
			ActualGCContent:   0.4982,
			ActualComplexity:  0.9430,
		},
		{
			Population:        "CEU",
			Region:            "Europe",
			ProcessedVariants: 80801086,
			TrainingTime:      parseTrainingTime("1h26m29.748681875s"),
			ActualGCContent:   0.4982,
			ActualComplexity:  0.9430,
		},
		{
			Population:        "GBR",
			Region:            "Europe",
			ProcessedVariants: 80801086,
			TrainingTime:      parseTrainingTime("1h26m44.880298917s"),
			ActualGCContent:   0.4982,
			ActualComplexity:  0.9430,
		},
		{
			Population:        "TSI",
			Region:            "Europe",
			ProcessedVariants: 80801086,
			TrainingTime:      parseTrainingTime("1h28m59.034007542s"),
			ActualGCContent:   0.4982,
			ActualComplexity:  0.9430,
		},
	}
}

func parseTrainingTime(timeStr string) time.Duration {
	duration, err := time.ParseDuration(timeStr)
	if err != nil {
		log.Printf("Warning: Failed to parse training time %s: %v", timeStr, err)
		return time.Duration(0)
	}
	return duration
}

func writeAuthenticVCFFASTA(filepath string, pop AuthenticVCFPopulation, length, seqNum int, sequence []byte) error {
	file, err := os.Create(filepath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	// Calculate sequence metrics
	complexity := calculateSequenceComplexity(sequence, 8)
	gcContent := calculateGCContent(sequence)

	// Write FASTA header; processed_variants, training_time, model_gc and
	// model_complexity are the hard-coded figures from the earlier log
	header := fmt.Sprintf(">%s_AUTHENTIC_VCF_L%d_S%02d population=%s region=%s length=%d complexity=%.4f gc_content=%.3f model_trained=false processed_variants=%d training_time=%v model_gc=%.4f model_complexity=%.4f\n",
		pop.Population, length, seqNum, pop.Population, pop.Region, len(sequence),
		complexity, gcContent, pop.ProcessedVariants, pop.TrainingTime,
		pop.ActualGCContent, pop.ActualComplexity)

	_, err = file.WriteString(header)
	if err != nil {
		return err
	}

	// Write sequence with line breaks
	for i := 0; i < len(sequence); i += 80 {
		end := i + 80
		if end > len(sequence) {
			end = len(sequence)
		}
		_, err = file.WriteString(string(sequence[i:end]) + "\n")
		if err != nil {
			return err
		}
	}

	return nil
}

func calculateSequenceComplexity(sequence []byte, k int) float64 {
	if len(sequence) < k {
		return 0.0
	}

	uniqueKmers := make(map[string]bool)
	for i := 0; i <= len(sequence)-k; i++ {
		kmer := string(sequence[i : i+k])
		uniqueKmers[kmer] = true
	}

	totalKmers := len(sequence) - k + 1
	return float64(len(uniqueKmers)) / float64(totalKmers)
}

func calculateGCContent(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0.0
	}

	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}

	return float64(gcCount) / float64(len(sequence))
}

func generateAuthenticSummaryReport(outputDir string, populations []AuthenticVCFPopulation, totalSequences int) {
	reportPath := filepath.Join(outputDir, "AUTHENTIC_VCF_TRAINING_REPORT.md")

	report := fmt.Sprintf(`# AUTHENTIC VCF TRAINING REPORT

## Summary

**Date**: %s
**Total Populations**: 6 of 11 labels
**Total Sequences Generated**: %d
**Generation Method**: fresh, untrained k=8 adaptive HMM per population; no trained model is loaded

## Populations (training figures hard-coded from an earlier log, not recomputed)

| Population | Region | Variants Processed | Training Time | GC Content | Complexity |
|------------|--------|-------------------|---------------|------------|------------|
`, time.Now().Format("2006-01-02 15:04:05"), totalSequences)

	for _, pop := range populations {
		report += fmt.Sprintf("| %s | %s | %d | %v | %.2f%% | %.4f |\n",
			pop.Population, pop.Region, pop.ProcessedVariants,
			pop.TrainingTime, pop.ActualGCContent*100, pop.ActualComplexity)
	}

	report += fmt.Sprintf(`
## Training Statistics (hard-coded from the earlier log)

- **Total Variants Processed**: 484,806,516
- **Average Training Time**: ~1h27m per population
- **Average GC Content**: 49.82%%
- **Average Complexity**: 0.9430

## Failed Populations

The following populations failed VCF training due to missing data files:
- CHB (East Asia) - No VCF files found
- JPT (East Asia) - No VCF files found  
- CHS (East Asia) - No VCF files found
- MXL (Americas) - No VCF files found
- PUR (Americas) - No VCF files found

## Sequence Generation

- **Sequence Lengths**: 1kb, 5kb, 10kb
- **Sequences per Length**: 5
- **Total per Population**: 15 sequences
- **File Format**: FASTA with detailed metadata headers

## Next Steps

1. Load the trained population models instead of generating from untrained ones
2. Investigation of missing VCF data for failed populations
3. Comparative analysis with synthetic-trained models

---

**Generated**: %s
**Status**: Sequences come from untrained models and are not suitable for validation
`, time.Now().Format("2006-01-02 15:04:05"))

	err := os.WriteFile(reportPath, []byte(report), 0600)
	if err != nil {
		log.Printf("Warning: Failed to write summary report: %v", err)
	} else {
		fmt.Printf("📄 Summary report created: %s\n", reportPath)
	}
}
