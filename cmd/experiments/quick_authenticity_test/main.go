// Command quick_authenticity_test trains a k=4 adaptive HMM on the reference
// genome in data/validation/reference_genome (path hard-coded, relative to
// the working directory), generates three 50 kb sequences and logs their
// complexity, GC content, longest homopolymer, CpG ratio and a 4-mer
// Jensen-Shannon divergence against a sample of reference k-mers, each
// against a fixed threshold.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
// Passing its four thresholds does not show the sequences are
// indistinguishable from real data.
package main

import (
	"log"
	"math"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func main() {
	log.Printf("=== Quick Authenticity Assessment ===")
	log.Printf("Threshold used here: JS Divergence < 0.01")

	// Load reference genome
	refGenomeDir := "data/validation/reference_genome"
	refGenome := reference.NewSimpleReferenceGenome(refGenomeDir)

	// Test 4-mer adaptive model (most promising)
	k := 4

	log.Printf("\n=== Testing Adaptive Model ===")
	adaptiveModel := authenticity.NewAdaptiveHMMModel(k)

	start := time.Now()
	err := adaptiveModel.TrainFromReferenceGenome(refGenome)
	if err != nil {
		log.Fatalf("Adaptive model training failed: %v", err)
	}
	trainingTime := time.Since(start)
	log.Printf("Adaptive model trained in: %v", trainingTime)

	// Generate test sequences
	log.Printf("\n=== Generating Test Sequences ===")
	sequenceLength := 50000
	numSequences := 3

	var sequences [][]byte
	var totalGenTime time.Duration

	for i := 0; i < numSequences; i++ {
		result, err := adaptiveModel.GenerateSequence(sequenceLength)
		if err != nil {
			log.Fatalf("Generation failed: %v", err)
		}
		sequences = append(sequences, result.Sequence)
		totalGenTime += result.GenerationTime
	}

	avgGenTime := totalGenTime / time.Duration(numSequences)
	log.Printf("Generated %d sequences, avg time: %v", numSequences, avgGenTime)

	// Calculate key metrics
	log.Printf("\n=== Key Metrics ===")

	// 1. Sequence Complexity
	complexity := calculateSequenceComplexity(sequences[0], k)
	log.Printf("Sequence Complexity: %.3f (target: > 0.5)", complexity)

	// 2. GC Content
	gcContent := calculateGCContent(sequences)
	gcError := math.Abs(gcContent - 0.42)
	log.Printf("GC Content: %.2f%% (target: 42%%, error: %.4f)", gcContent*100, gcError)

	// 3. Homopolymer Analysis
	maxHomopolymer := findMaxHomopolymerRun(sequences)
	log.Printf("Max Homopolymer: %d (target: < 15)", maxHomopolymer)

	// 4. CpG Depletion
	cpgDepletion := calculateCpGDepletion(sequences)
	log.Printf("CpG Depletion: %.3f (target: ~0.25)", cpgDepletion)

	// 5. Quick JS Divergence Estimate
	log.Printf("\n=== Quick JS Divergence Estimate ===")
	refKmers := extractLimitedReferenceKmers(refGenome, k, 100000)
	log.Printf("Using %d reference k-mers for estimation", len(refKmers))

	jsDivergence := calculateJSDivergence(sequences, refKmers, k)
	log.Printf("JS Divergence: %.6f (target: < 0.01)", jsDivergence)

	// Overall Assessment
	log.Printf("\n=== Overall Assessment ===")

	complexityPass := complexity > 0.5
	gcPass := gcError < 0.05
	homopolymerPass := maxHomopolymer < 15
	jsPass := jsDivergence < 0.01

	log.Printf("Complexity Check:   %v (%.3f > 0.5)", complexityPass, complexity)
	log.Printf("GC Content Check:   %v (%.4f < 0.05)", gcPass, gcError)
	log.Printf("Homopolymer Check:  %v (%d < 15)", homopolymerPass, maxHomopolymer)
	log.Printf("JS Divergence Check: %v (%.6f < 0.01)", jsPass, jsDivergence)

	totalPassed := 0
	if complexityPass {
		totalPassed++
	}
	if gcPass {
		totalPassed++
	}
	if homopolymerPass {
		totalPassed++
	}
	if jsPass {
		totalPassed++
	}

	log.Printf("\nPASSED: %d/4 criteria", totalPassed)

	if totalPassed == 4 {
		log.Printf("✅ All four threshold checks passed")
		log.Printf("   These are summary statistics only; they do not show the sequences are indistinguishable from real data")
	} else if totalPassed >= 2 {
		log.Printf("⚠️  Some threshold checks passed")
	} else {
		log.Printf("❌ Fewer than two threshold checks passed")
	}

	// Performance Summary
	log.Printf("\n=== Performance Summary ===")
	log.Printf("Training Speed: %.2f seconds", trainingTime.Seconds())
	log.Printf("Generation Speed: %.2f ms per 1K bases", float64(avgGenTime.Nanoseconds())/1000000.0/float64(sequenceLength/1000))
	log.Printf("Memory Usage: not measured (k=%d model)", k)

	// Next Steps
	log.Printf("\n=== Next Steps ===")
	if jsPass {
		log.Printf("1. Compare against held-out real data with a classifier")
		log.Printf("2. Performance benchmarking")
		log.Printf("3. Integration with the sequencer simulation")
	} else {
		log.Printf("1. 🔄 Continue model refinement for JS divergence")
		log.Printf("2. 🔄 Expand training data if needed")
		log.Printf("3. 🔄 Fine-tune biological constraints")
	}

	// Generation Stats (if available)
	if stats := adaptiveModel.GetGenerationStats(); stats.TotalSequences > 0 {
		log.Printf("\n=== Generation Statistics ===")
		log.Printf("Total Sequences: %d", stats.TotalSequences)
		log.Printf("Average Complexity: %.3f", stats.AverageComplexity)
		log.Printf("Average GC Content: %.2f%%", stats.AverageGCContent*100)
		log.Printf("Max Homopolymer: %d", stats.MaxHomopolymer)
	}
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

func calculateGCContent(sequences [][]byte) float64 {
	totalBases := 0
	gcCount := 0

	for _, seq := range sequences {
		for _, base := range seq {
			totalBases++
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}
	}

	return float64(gcCount) / float64(totalBases)
}

func findMaxHomopolymerRun(sequences [][]byte) int {
	maxHomopolymer := 0

	for _, seq := range sequences {
		if len(seq) == 0 {
			continue
		}

		maxRun := 1
		currentRun := 1

		for i := 1; i < len(seq); i++ {
			if seq[i] == seq[i-1] {
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

		if maxRun > maxHomopolymer {
			maxHomopolymer = maxRun
		}
	}

	return maxHomopolymer
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

func extractLimitedReferenceKmers(refGenome *reference.SimpleReferenceGenome, k int, maxKmers int) []string {
	var kmers []string
	count := 0

	// Process chromosome 22 (smaller)
	chr, err := refGenome.LoadChromosome("chr22")
	if err != nil {
		log.Printf("Warning: Could not load chromosome chr22: %v", err)
		return kmers
	}

	sequence := chr.GetSequence()
	log.Printf("Processing chromosome chr22 (%d bases)", len(sequence))

	for i := 0; i <= len(sequence)-k && count < maxKmers; i++ {
		kmer := string(sequence[i : i+k])
		// Only include valid k-mers
		if isValidKmer(kmer) {
			kmers = append(kmers, kmer)
			count++
		}
	}

	return kmers
}

func isValidKmer(kmer string) bool {
	for _, base := range kmer {
		if base != 'A' && base != 'T' && base != 'G' && base != 'C' {
			return false
		}
	}
	return true
}

func calculateJSDivergence(sequences [][]byte, refKmers []string, k int) float64 {
	// Calculate k-mer frequencies for generated sequences
	genKmers := make(map[string]int)
	genTotal := 0

	for _, seq := range sequences {
		for i := 0; i <= len(seq)-k; i++ {
			kmer := string(seq[i : i+k])
			if isValidKmer(kmer) {
				genKmers[kmer]++
				genTotal++
			}
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
	return calculateJSDivergenceFromProbs(genProb, refProb)
}

func calculateJSDivergenceFromProbs(p, q map[string]float64) float64 {
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
