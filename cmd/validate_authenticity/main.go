package main

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"sort"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// ValidationResult represents the results of authenticity validation
type ValidationResult struct {
	ModelK              int                        `json:"model_k"`
	TrainingDataset     string                     `json:"training_dataset"`
	ValidationDataset   string                     `json:"validation_dataset"`
	JSdivergence        float64                    `json:"js_divergence"`
	KLDivergence        float64                    `json:"kl_divergence"`
	GCContentDiff       float64                    `json:"gc_content_diff"`
	KmerDistribution    map[string]float64         `json:"kmer_distribution"`
	BiologicalMetrics   BiologicalMetrics          `json:"biological_metrics"`
	ValidationStatus    string                     `json:"validation_status"`
	ValidationTime      time.Time                  `json:"validation_time"`
	Notes               string                     `json:"notes"`
}

// BiologicalMetrics represents biological realism metrics
type BiologicalMetrics struct {
	GCContent           float64 `json:"gc_content"`
	ATContent           float64 `json:"at_content"`
	MaxHomopolymer      int     `json:"max_homopolymer"`
	SequenceComplexity  float64 `json:"sequence_complexity"`
	CpGDepletion        float64 `json:"cpg_depletion"`
	CodonUsageBias      float64 `json:"codon_usage_bias"`
}

// JSThreshold is the target threshold for Jensen-Shannon divergence
const JSThreshold = 0.01

func main() {
	log.Printf("=== Phantom Science Authenticity Validation Pipeline ===")
	log.Printf("Target JS Divergence: < %.3f", JSThreshold)
	
	// Create validation results directory
	resultsDir := "validation_results"
	if err := os.MkdirAll(resultsDir, 0750); err != nil {
		log.Fatalf("Failed to create results directory: %v", err)
	}
	
	// Load reference genome for validation
	refGenomeDir := "data/validation/reference_genome"
	if _, err := os.Stat(refGenomeDir); os.IsNotExist(err) {
		log.Fatalf("Reference genome directory not found: %s", refGenomeDir)
	}
	
	refGenome := reference.NewSimpleReferenceGenome(refGenomeDir)
	
	// Test all trained models
	modelSizes := []int{3, 4, 5, 6}
	var allResults []ValidationResult
	
	for _, k := range modelSizes {
		log.Printf("\n=== Validating %d-mer model ===", k)
		
		// Load trained model
		model := authenticity.NewHMMKmerModel(k)
		err := model.TrainFromReferenceGenome(refGenome)
		if err != nil {
			log.Printf("Failed to train %d-mer model: %v", k, err)
			continue
		}
		
		// Run validation
		result := validateModel(model, refGenome, k)
		allResults = append(allResults, result)
		
		// Log results
		logValidationResult(result)
	}
	
	// Generate comprehensive report
	generateValidationReport(allResults, resultsDir)
	
	// Check if any models meet the authenticity threshold
	passedModels := 0
	for _, result := range allResults {
		if result.JSdivergence < JSThreshold {
			passedModels++
		}
	}
	
	log.Printf("\n=== Final Validation Summary ===")
	log.Printf("Models tested: %d", len(allResults))
	log.Printf("Models passing JS threshold (< %.3f): %d", JSThreshold, passedModels)
	
	if passedModels > 0 {
		log.Printf("✓ AUTHENTICITY VALIDATION PASSED")
		log.Printf("✓ Generated sequences are statistically indistinguishable from real genomic data")
	} else {
		log.Printf("⚠ AUTHENTICITY VALIDATION REQUIRES IMPROVEMENT")
		log.Printf("⚠ Generated sequences show detectable statistical differences")
	}
	
	log.Printf("\nValidation results saved to: %s", resultsDir)
}

// validateModel performs comprehensive validation of a trained model
func validateModel(model *authenticity.HMMKmerModel, refGenome *reference.SimpleReferenceGenome, k int) ValidationResult {
	result := ValidationResult{
		ModelK:            k,
		TrainingDataset:   "GRCh38_chr21_chr22",
		ValidationDataset: "GRCh38_chr21_chr22",
		ValidationTime:    time.Now(),
	}
	
	// Generate test sequences
	log.Printf("Generating validation sequences...")
	sequenceLength := 100000
	numSequences := 10
	
	var generatedSequences [][]byte
	var generatedKmers []string
	
	for i := 0; i < numSequences; i++ {
		genResult, err := model.GenerateSequence(sequenceLength)
		if err != nil {
			log.Printf("Generation failed: %v", err)
			result.ValidationStatus = "FAILED"
			result.Notes = fmt.Sprintf("Generation failed: %v", err)
			return result
		}
		
		generatedSequences = append(generatedSequences, genResult.Sequence)
		
		// Extract k-mers from generated sequence
		for j := 0; j <= len(genResult.Sequence)-k; j++ {
			kmer := string(genResult.Sequence[j : j+k])
			generatedKmers = append(generatedKmers, kmer)
		}
	}
	
	// Extract k-mers from reference genome
	log.Printf("Extracting k-mers from reference genome...")
	var referenceKmers []string
	
	// Process chromosomes 21 and 22
	for _, chrName := range []string{"chr21", "chr22"} {
		chr, err := refGenome.LoadChromosome(chrName)
		if err != nil {
			log.Printf("Failed to load chromosome %s: %v", chrName, err)
			continue
		}
		
		sequence := chr.GetSequence()
		for i := 0; i <= len(sequence)-k; i++ {
			kmer := string(sequence[i : i+k])
			referenceKmers = append(referenceKmers, kmer)
		}
	}
	
	// Calculate k-mer frequency distributions
	log.Printf("Calculating k-mer frequency distributions...")
	genFreq := calculateKmerFrequencies(generatedKmers)
	refFreq := calculateKmerFrequencies(referenceKmers)
	
	// Calculate Jensen-Shannon divergence
	result.JSdivergence = calculateJSDivergence(genFreq, refFreq)
	result.KLDivergence = calculateKLDivergence(genFreq, refFreq)
	result.KmerDistribution = genFreq
	
	// Calculate biological metrics
	result.BiologicalMetrics = calculateBiologicalMetrics(generatedSequences, k)
	
	// Calculate GC content difference
	genGC := result.BiologicalMetrics.GCContent
	refGC := calculateGCContentFromKmers(referenceKmers)
	result.GCContentDiff = math.Abs(genGC - refGC)
	
	// Determine validation status
	if result.JSdivergence < JSThreshold {
		result.ValidationStatus = "PASSED"
		result.Notes = fmt.Sprintf("JS divergence %.6f < threshold %.3f", result.JSdivergence, JSThreshold)
	} else {
		result.ValidationStatus = "FAILED"
		result.Notes = fmt.Sprintf("JS divergence %.6f > threshold %.3f", result.JSdivergence, JSThreshold)
	}
	
	return result
}

// calculateKmerFrequencies calculates k-mer frequency distribution
func calculateKmerFrequencies(kmers []string) map[string]float64 {
	counts := make(map[string]int)
	total := 0
	
	for _, kmer := range kmers {
		counts[kmer]++
		total++
	}
	
	frequencies := make(map[string]float64)
	for kmer, count := range counts {
		frequencies[kmer] = float64(count) / float64(total)
	}
	
	return frequencies
}

// calculateJSDivergence calculates Jensen-Shannon divergence between two distributions
func calculateJSDivergence(p, q map[string]float64) float64 {
	// Get all unique k-mers
	allKmers := make(map[string]bool)
	for kmer := range p {
		allKmers[kmer] = true
	}
	for kmer := range q {
		allKmers[kmer] = true
	}
	
	// Calculate M = (P + Q) / 2
	m := make(map[string]float64)
	for kmer := range allKmers {
		pVal := p[kmer]
		qVal := q[kmer]
		m[kmer] = (pVal + qVal) / 2.0
	}
	
	// Calculate JS divergence
	js := 0.0
	for kmer := range allKmers {
		pVal := p[kmer]
		qVal := q[kmer]
		mVal := m[kmer]
		
		if pVal > 0 && mVal > 0 {
			js += pVal * math.Log2(pVal/mVal)
		}
		if qVal > 0 && mVal > 0 {
			js += qVal * math.Log2(qVal/mVal)
		}
	}
	
	return js / 2.0
}

// calculateKLDivergence calculates Kullback-Leibler divergence
func calculateKLDivergence(p, q map[string]float64) float64 {
	kl := 0.0
	for kmer, pVal := range p {
		qVal := q[kmer]
		if pVal > 0 && qVal > 0 {
			kl += pVal * math.Log2(pVal/qVal)
		}
	}
	return kl
}

// calculateBiologicalMetrics calculates biological realism metrics
func calculateBiologicalMetrics(sequences [][]byte, k int) BiologicalMetrics {
	var metrics BiologicalMetrics
	
	totalBases := 0
	gcCount := 0
	atCount := 0
	maxHomopolymer := 0
	
	for _, seq := range sequences {
		// GC/AT content
		for _, base := range seq {
			totalBases++
			if base == 'G' || base == 'C' {
				gcCount++
			} else if base == 'A' || base == 'T' {
				atCount++
			}
		}
		
		// Homopolymer runs
		maxRun := findMaxHomopolymerRun(seq)
		if maxRun > maxHomopolymer {
			maxHomopolymer = maxRun
		}
	}
	
	if totalBases > 0 {
		metrics.GCContent = float64(gcCount) / float64(totalBases)
		metrics.ATContent = float64(atCount) / float64(totalBases)
	}
	metrics.MaxHomopolymer = maxHomopolymer
	
	// Calculate sequence complexity (using first sequence as representative)
	if len(sequences) > 0 {
		metrics.SequenceComplexity = calculateSequenceComplexity(sequences[0])
	}
	
	// Calculate CpG depletion
	metrics.CpGDepletion = calculateCpGDepletion(sequences)
	
	// Codon usage bias (simplified)
	metrics.CodonUsageBias = 0.0 // Placeholder
	
	return metrics
}

// calculateGCContentFromKmers calculates GC content from k-mer list
func calculateGCContentFromKmers(kmers []string) float64 {
	gcCount := 0
	totalBases := 0
	
	for _, kmer := range kmers {
		for _, base := range kmer {
			totalBases++
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}
	}
	
	return float64(gcCount) / float64(totalBases)
}

// findMaxHomopolymerRun finds the maximum homopolymer run in a sequence
func findMaxHomopolymerRun(sequence []byte) int {
	if len(sequence) == 0 {
		return 0
	}
	
	maxRun := 1
	currentRun := 1
	
	for i := 1; i < len(sequence); i++ {
		if sequence[i] == sequence[i-1] {
			currentRun++
		} else {
			if currentRun > maxRun {
				maxRun = currentRun
			}
			currentRun = 1
		}
	}
	
	if currentRun > maxRun {
		maxRun = currentRun
	}
	
	return maxRun
}

// calculateSequenceComplexity calculates sequence complexity using k-mer diversity
func calculateSequenceComplexity(sequence []byte) float64 {
	if len(sequence) < 3 {
		return 0.0
	}
	
	kmerCount := make(map[string]int)
	for i := 0; i <= len(sequence)-3; i++ {
		kmer := string(sequence[i : i+3])
		kmerCount[kmer]++
	}
	
	maxPossibleKmers := len(sequence) - 2
	return float64(len(kmerCount)) / float64(maxPossibleKmers)
}

// calculateCpGDepletion calculates CpG depletion ratio
func calculateCpGDepletion(sequences [][]byte) float64 {
	cgCount := 0
	totalDinucleotides := 0
	
	for _, seq := range sequences {
		for i := 0; i < len(seq)-1; i++ {
			totalDinucleotides++
			if seq[i] == 'C' && seq[i+1] == 'G' {
				cgCount++
			}
		}
	}
	
	if totalDinucleotides == 0 {
		return 0.0
	}
	
	// Expected CpG frequency is ~1% in human genome
	expectedCpG := 0.01
	observedCpG := float64(cgCount) / float64(totalDinucleotides)
	
	return observedCpG / expectedCpG
}

// logValidationResult logs the validation results
func logValidationResult(result ValidationResult) {
	log.Printf("=== %d-mer Model Validation Results ===", result.ModelK)
	log.Printf("JS Divergence: %.6f (threshold: %.3f)", result.JSdivergence, JSThreshold)
	log.Printf("KL Divergence: %.6f", result.KLDivergence)
	log.Printf("GC Content Diff: %.4f", result.GCContentDiff)
	log.Printf("Status: %s", result.ValidationStatus)
	log.Printf("Biological Metrics:")
	log.Printf("  GC Content: %.2f%%", result.BiologicalMetrics.GCContent*100)
	log.Printf("  AT Content: %.2f%%", result.BiologicalMetrics.ATContent*100)
	log.Printf("  Max Homopolymer: %d", result.BiologicalMetrics.MaxHomopolymer)
	log.Printf("  Sequence Complexity: %.3f", result.BiologicalMetrics.SequenceComplexity)
	log.Printf("  CpG Depletion: %.3f", result.BiologicalMetrics.CpGDepletion)
	log.Printf("Notes: %s", result.Notes)
	log.Printf("")
}

// generateValidationReport generates a comprehensive validation report
func generateValidationReport(results []ValidationResult, outputDir string) {
	// Save JSON report
	jsonFile := filepath.Join(outputDir, "validation_report.json")
	jsonData, err := json.MarshalIndent(results, "", "  ")
	if err != nil {
		log.Printf("Failed to marshal JSON: %v", err)
		return
	}
	
	if err := os.WriteFile(jsonFile, jsonData, 0600); err != nil {
		log.Printf("Failed to write JSON report: %v", err)
		return
	}
	
	// Generate markdown report
	mdFile := filepath.Join(outputDir, "validation_report.md")
	generateMarkdownReport(results, mdFile)
	
	log.Printf("Validation report saved to: %s", outputDir)
}

// generateMarkdownReport generates a markdown validation report
func generateMarkdownReport(results []ValidationResult, filename string) {
	report := `# Phantom Science Authenticity Validation Report

## Summary

This report contains the results of authenticity validation for the Phantom Science genomics engine. The validation measures how statistically similar generated sequences are to real genomic data using Jensen-Shannon divergence.

**Target Threshold:** JS Divergence < 0.01

## Results

| Model | JS Divergence | KL Divergence | GC Content Diff | Status | Notes |
|-------|---------------|---------------|-----------------|--------|-------|
`

	// Sort results by model k
	sort.Slice(results, func(i, j int) bool {
		return results[i].ModelK < results[j].ModelK
	})
	
	for _, result := range results {
		report += fmt.Sprintf("| %d-mer | %.6f | %.6f | %.4f | %s | %s |\n",
			result.ModelK,
			result.JSdivergence,
			result.KLDivergence,
			result.GCContentDiff,
			result.ValidationStatus,
			result.Notes)
	}
	
	report += "\n## Biological Realism Metrics\n\n"
	
	for _, result := range results {
		report += fmt.Sprintf("### %d-mer Model\n\n", result.ModelK)
		report += fmt.Sprintf("- **GC Content:** %.2f%%\n", result.BiologicalMetrics.GCContent*100)
		report += fmt.Sprintf("- **AT Content:** %.2f%%\n", result.BiologicalMetrics.ATContent*100)
		report += fmt.Sprintf("- **Max Homopolymer:** %d\n", result.BiologicalMetrics.MaxHomopolymer)
		report += fmt.Sprintf("- **Sequence Complexity:** %.3f\n", result.BiologicalMetrics.SequenceComplexity)
		report += fmt.Sprintf("- **CpG Depletion:** %.3f\n", result.BiologicalMetrics.CpGDepletion)
		report += "\n"
	}
	
	report += "## Validation Methodology\n\n"
	report += "1. **Training:** HMM k-mer models trained on GRCh38 chromosomes 21 and 22\n"
	report += "2. **Generation:** 10 sequences of 100,000 bases each generated per model\n"
	report += "3. **Comparison:** k-mer frequency distributions compared between generated and reference sequences\n"
	report += "4. **Metrics:** Jensen-Shannon divergence, KL divergence, and biological realism metrics calculated\n"
	report += "5. **Threshold:** JS divergence < 0.01 indicates statistically indistinguishable sequences\n"
	
	if err := os.WriteFile(filename, []byte(report), 0600); err != nil {
		log.Printf("Failed to write markdown report: %v", err)
	}
}