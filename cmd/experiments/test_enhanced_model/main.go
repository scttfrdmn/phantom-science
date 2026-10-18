// Command test_enhanced_model trains the original and enhanced 4-mer HMM
// models on the reference genome in data/validation/reference_genome
// (relative to the current directory), generates sequences from each and
// logs GC content, homopolymer length, complexity, CpG ratio and a
// Jensen-Shannon divergence against k-mers sampled from chr21/chr22.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"log"
	"math"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func main() {
	log.Printf("=== Enhanced HMM Model Comparison ===")

	// Load reference genome
	refGenomeDir := "data/validation/reference_genome"
	refGenome := reference.NewSimpleReferenceGenome(refGenomeDir)

	// Test 4-mer models (original vs enhanced)
	k := 4

	log.Printf("\n=== Training Original Model ===")
	originalModel := authenticity.NewHMMKmerModel(k)
	start := time.Now()
	err := originalModel.TrainFromReferenceGenome(refGenome)
	if err != nil {
		log.Fatalf("Original model training failed: %v", err)
	}
	originalTrainingTime := time.Since(start)
	log.Printf("Original model trained in: %v", originalTrainingTime)

	log.Printf("\n=== Training Enhanced Model ===")
	enhancedModel := authenticity.NewEnhancedHMMKmerModel(k, "CEU")
	start = time.Now()
	err = enhancedModel.TrainFromReferenceGenome(refGenome)
	if err != nil {
		log.Fatalf("Enhanced model training failed: %v", err)
	}
	enhancedTrainingTime := time.Since(start)
	log.Printf("Enhanced model trained in: %v", enhancedTrainingTime)

	// Generate sequences from both models
	sequenceLength := 50000
	numSequences := 5

	log.Printf("\n=== Generating Sequences ===")

	// Original model sequences
	var originalSequences [][]byte
	var originalGenerationTime time.Duration

	for i := 0; i < numSequences; i++ {
		result, err := originalModel.GenerateSequence(sequenceLength)
		if err != nil {
			log.Fatalf("Original generation failed: %v", err)
		}
		originalSequences = append(originalSequences, result.Sequence)
		originalGenerationTime += result.GenerationTime
	}

	// Enhanced model sequences
	var enhancedSequences [][]byte
	var enhancedGenerationTime time.Duration

	for i := 0; i < numSequences; i++ {
		result, err := enhancedModel.GenerateSequence(sequenceLength)
		if err != nil {
			log.Fatalf("Enhanced generation failed: %v", err)
		}
		enhancedSequences = append(enhancedSequences, result.Sequence)
		enhancedGenerationTime += result.GenerationTime
	}

	// Compare biological realism
	log.Printf("\n=== Biological Realism Comparison ===")

	originalMetrics := calculateBiologicalMetrics(originalSequences)
	enhancedMetrics := calculateBiologicalMetrics(enhancedSequences)

	log.Printf("Original Model Metrics:")
	printMetrics(originalMetrics)
	log.Printf("Average generation time: %v", originalGenerationTime/time.Duration(numSequences))

	log.Printf("\nEnhanced Model Metrics:")
	printMetrics(enhancedMetrics)
	log.Printf("Average generation time: %v", enhancedGenerationTime/time.Duration(numSequences))

	// Calculate improvements
	log.Printf("\n=== Improvements ===")

	gcImprovement := math.Abs(enhancedMetrics.GCContent-0.42) < math.Abs(originalMetrics.GCContent-0.42)
	complexityImprovement := enhancedMetrics.SequenceComplexity > originalMetrics.SequenceComplexity
	homopolymerImprovement := enhancedMetrics.MaxHomopolymer < originalMetrics.MaxHomopolymer

	log.Printf("GC Content closer to 42%%: %v", gcImprovement)
	log.Printf("Sequence complexity improved: %v (%.3f -> %.3f)", complexityImprovement, originalMetrics.SequenceComplexity, enhancedMetrics.SequenceComplexity)
	log.Printf("Max homopolymer reduced: %v (%d -> %d)", homopolymerImprovement, originalMetrics.MaxHomopolymer, enhancedMetrics.MaxHomopolymer)

	// Quick JS divergence estimate
	log.Printf("\n=== Quick JS Divergence Estimate ===")

	// Extract k-mers from reference (sample)
	refKmers := extractReferenceKmers(refGenome, k, 100000)

	// Calculate JS divergence for both models
	originalJSD := calculateQuickJSDivergence(originalSequences, refKmers, k)
	enhancedJSD := calculateQuickJSDivergence(enhancedSequences, refKmers, k)

	log.Printf("Original JS Divergence: %.6f", originalJSD)
	log.Printf("Enhanced JS Divergence: %.6f", enhancedJSD)

	improvement := originalJSD - enhancedJSD
	log.Printf("JS Divergence improvement: %.6f", improvement)

	if enhancedJSD < 0.01 {
		log.Printf("✓ Enhanced model JS divergence is below this program's 0.01 threshold")
	} else if improvement > 0 {
		log.Printf("✓ Enhanced model shows improvement but needs further refinement")
	} else {
		log.Printf("⚠ Enhanced model needs additional improvements")
	}

	log.Printf("\n=== Summary ===")
	log.Printf("Enhanced model shows:")
	if gcImprovement {
		log.Printf("✓ Better GC content accuracy")
	}
	if complexityImprovement {
		log.Printf("✓ Improved sequence complexity")
	}
	if homopolymerImprovement {
		log.Printf("✓ Reduced homopolymer runs")
	}
	if improvement > 0 {
		log.Printf("✓ Lower Jensen-Shannon divergence")
	}
}

type BiologicalMetrics struct {
	GCContent          float64
	ATContent          float64
	MaxHomopolymer     int
	SequenceComplexity float64
	CpGDepletion       float64
}

func calculateBiologicalMetrics(sequences [][]byte) BiologicalMetrics {
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

	metrics.GCContent = float64(gcCount) / float64(totalBases)
	metrics.ATContent = float64(atCount) / float64(totalBases)
	metrics.MaxHomopolymer = maxHomopolymer

	// Calculate sequence complexity (using first sequence as representative)
	if len(sequences) > 0 {
		metrics.SequenceComplexity = calculateSequenceComplexity(sequences[0])
	}

	// Calculate CpG depletion
	metrics.CpGDepletion = calculateCpGDepletion(sequences)

	return metrics
}

func printMetrics(metrics BiologicalMetrics) {
	log.Printf("  GC Content: %.2f%%", metrics.GCContent*100)
	log.Printf("  AT Content: %.2f%%", metrics.ATContent*100)
	log.Printf("  Max Homopolymer: %d", metrics.MaxHomopolymer)
	log.Printf("  Sequence Complexity: %.3f", metrics.SequenceComplexity)
	log.Printf("  CpG Depletion: %.3f", metrics.CpGDepletion)
}

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

	expectedCpG := 0.01
	observedCpG := float64(cgCount) / float64(totalDinucleotides)

	return observedCpG / expectedCpG
}

func extractReferenceKmers(refGenome *reference.SimpleReferenceGenome, k int, maxKmers int) []string {
	var kmers []string
	count := 0

	// Process chromosomes 21 and 22
	for _, chrName := range []string{"chr21", "chr22"} {
		if count >= maxKmers {
			break
		}

		chr, err := refGenome.LoadChromosome(chrName)
		if err != nil {
			continue
		}

		sequence := chr.GetSequence()
		for i := 0; i <= len(sequence)-k && count < maxKmers; i++ {
			kmer := string(sequence[i : i+k])
			kmers = append(kmers, kmer)
			count++
		}
	}

	return kmers
}

func calculateQuickJSDivergence(sequences [][]byte, refKmers []string, k int) float64 {
	// Calculate k-mer frequencies for generated sequences
	genKmers := make(map[string]int)
	genTotal := 0

	for _, seq := range sequences {
		for i := 0; i <= len(seq)-k; i++ {
			kmer := string(seq[i : i+k])
			genKmers[kmer]++
			genTotal++
		}
	}

	// Calculate k-mer frequencies for reference
	refKmerMap := make(map[string]int)
	for _, kmer := range refKmers {
		refKmerMap[kmer]++
	}
	refTotal := len(refKmers)

	// Convert to probabilities
	genProb := make(map[string]float64)
	for kmer, count := range genKmers {
		genProb[kmer] = float64(count) / float64(genTotal)
	}

	refProb := make(map[string]float64)
	for kmer, count := range refKmerMap {
		refProb[kmer] = float64(count) / float64(refTotal)
	}

	// Calculate JS divergence
	return calculateJSDivergence(genProb, refProb)
}

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
