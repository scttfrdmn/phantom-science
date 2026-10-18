// Command astrapi_full_genome_training trains a k=4 adaptive HMM on the
// reference genome in the hard-coded directory data/astrapi_training/human,
// generates five 50 kb sequences and prints their complexity, GC content,
// longest homopolymer, CpG ratio and 4-mer Jensen-Shannon divergence against
// chr1. The ratios it prints compare with baselines (complexity 0.005, JS
// 0.055) hard-coded from earlier runs, and complexity is computed on the
// first sequence only.
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
	log.Printf("=== Astrapi High-Capacity Full Genome Training ===")
	log.Printf("Training on the reference genome in data/astrapi_training/human")
	log.Printf("Baseline complexity for comparison: 0.005 (hard-coded from an earlier run)")
	log.Printf("Target: JS divergence < 0.01")

	// Step 1: Initialize with astrapi genome data
	log.Printf("\n[STEP 1] Initializing full genome access from astrapi...")
	refGenomeDir := "data/astrapi_training/human"
	refGenome := reference.NewSimpleReferenceGenome(refGenomeDir)
	log.Printf("✓ Astrapi genome access initialized")
	log.Printf("✓ Data location: %s", refGenomeDir)

	// Step 2: Create adaptive model for full genome training
	log.Printf("\n[STEP 2] Creating high-capacity adaptive HMM model...")
	k := 4
	adaptiveModel := authenticity.NewAdaptiveHMMModel(k)
	log.Printf("✓ Adaptive model created with k=%d", k)
	log.Printf("✓ Model configured for full genome training:")
	log.Printf("  - Target complexity: > 0.1 (baseline 0.005)")
	log.Printf("  - Target GC content: 42%%")
	log.Printf("  - Max homopolymer: 8 bases")
	log.Printf("  - CpG depletion: 25%%")
	log.Printf("  - Training data: every chromosome the reference genome lists")

	// Step 3: Full genome training (this will take significant time)
	log.Printf("\n[STEP 3] Training adaptive model on complete human genome...")
	log.Printf("Training on the reference genome - this can take a long time and a lot of memory...")

	start := time.Now()
	err := adaptiveModel.TrainFromReferenceGenome(refGenome)
	if err != nil {
		log.Fatalf("❌ Training failed: %v", err)
	}
	trainingTime := time.Since(start)

	log.Printf("✅ TRAINING COMPLETED SUCCESSFULLY!")
	log.Printf("✅ Training time: %v", trainingTime)
	log.Printf("✅ Total k-mers processed: %d", adaptiveModel.TotalKmers())

	// Step 4: Generate test sequences for validation
	log.Printf("\n[STEP 4] Generating test sequences from full genome model...")
	sequenceLength := 50000
	numSequences := 5 // More sequences for better validation

	var sequences [][]byte
	var totalGenTime time.Duration

	for i := 0; i < numSequences; i++ {
		log.Printf("  Generating sequence %d/%d (%d bases)...", i+1, numSequences, sequenceLength)

		result, err := adaptiveModel.GenerateSequence(sequenceLength)
		if err != nil {
			log.Fatalf("❌ Generation failed: %v", err)
		}

		sequences = append(sequences, result.Sequence)
		totalGenTime += result.GenerationTime

		log.Printf("    ✓ Generated in %v", result.GenerationTime)
	}

	avgGenTime := totalGenTime / time.Duration(numSequences)
	generationSpeed := float64(sequenceLength) / avgGenTime.Seconds() / 1000.0
	log.Printf("✓ All sequences generated, average time: %v", avgGenTime)
	log.Printf("✓ Generation speed: %.1f kb/s", generationSpeed)

	// Step 5: Comprehensive biological metrics validation
	log.Printf("\n[STEP 5] Comprehensive biological realism validation...")

	// 5.1: Sequence Complexity (CRITICAL METRIC)
	log.Printf("  Calculating sequence complexity...")
	complexity := calculateSequenceComplexity(sequences[0], k)
	complexityImprovement := complexity / 0.005 // vs current 0.005
	complexityPass := complexity > 0.1          // Much higher threshold for full genome
	log.Printf("    Complexity: %.3f (target: > 0.1, %.1fx the 0.005 baseline) %s",
		complexity, complexityImprovement, passFailIcon(complexityPass))

	// 5.2: GC Content
	log.Printf("  Calculating GC content...")
	gcContent := calculateGCContent(sequences)
	gcError := math.Abs(gcContent - 0.42)
	gcPass := gcError < 0.02 // Tighter tolerance for full genome
	log.Printf("    GC Content: %.2f%% (target: 42%%, error: %.4f) %s",
		gcContent*100, gcError, passFailIcon(gcPass))

	// 5.3: Homopolymer Analysis
	log.Printf("  Analyzing homopolymer runs...")
	maxHomopolymer := findMaxHomopolymerRun(sequences)
	homopolymerPass := maxHomopolymer < 12 // Stricter for full genome
	log.Printf("    Max Homopolymer: %d (target: < 12) %s",
		maxHomopolymer, passFailIcon(homopolymerPass))

	// 5.4: CpG Depletion
	log.Printf("  Calculating CpG depletion...")
	cpgDepletion := calculateCpGDepletion(sequences)
	cpgPass := math.Abs(cpgDepletion-0.25) < 0.3
	log.Printf("    CpG Depletion: %.3f (target: ~0.25) %s",
		cpgDepletion, passFailIcon(cpgPass))

	// Step 6: Jensen-Shannon Divergence (CRITICAL AUTHENTICITY METRIC)
	log.Printf("\n[STEP 6] Calculating Jensen-Shannon divergence vs real genome...")
	log.Printf("  Loading reference sequences for comparison...")

	// Load a sample of reference genome for comparison
	referenceSequences, err := loadReferenceSequences(refGenome, 5, sequenceLength)
	if err != nil {
		log.Printf("⚠ Warning: Could not load reference sequences: %v", err)
	} else {
		jsDivergence := calculateJensenShannonDivergence(sequences, referenceSequences, k)
		jsPass := jsDivergence < 0.01
		jsImprovement := 0.055 / jsDivergence // vs current 0.055
		log.Printf("    JS Divergence: %.6f (target: < 0.01, %.1fx better than the 0.055 baseline) %s",
			jsDivergence, jsImprovement, passFailIcon(jsPass))
	}

	// Step 7: Overall assessment
	log.Printf("\n[STEP 7] Full Genome Training Results Summary...")
	overallPass := complexityPass && gcPass && homopolymerPass && cpgPass

	log.Printf("\n🧬 FULL GENOME TRAINING RESULTS:")
	log.Printf("====================================")
	log.Printf("Training Data:       %s", refGenomeDir)
	log.Printf("Training Time:       %v", trainingTime)
	log.Printf("K-mers Processed:    %d", adaptiveModel.TotalKmers())
	log.Printf("Complexity:          %.3f (%.1fx baseline)", complexity, complexityImprovement)
	log.Printf("Generation Speed:    %.1f kb/s", generationSpeed)
	log.Printf("Overall Status:      %s", passFailIcon(overallPass))

	if overallPass {
		log.Printf("\n✅ All four sequence metrics passed their thresholds")
	} else {
		log.Printf("\n⚠ One or more sequence metrics failed their thresholds")
	}

	log.Printf("\n📊 Next Steps:")
	log.Printf("1. Download 1000 Genomes population data")
	log.Printf("2. Train a population-aware model")
	log.Printf("3. Add mammalian genomes for universal patterns")
	log.Printf("4. Target complexity: 0.5+")
}

// Helper functions

func passFailIcon(pass bool) string {
	if pass {
		return "✅"
	}
	return "❌"
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

	if totalBases == 0 {
		return 0.0
	}
	return float64(gcCount) / float64(totalBases)
}

func findMaxHomopolymerRun(sequences [][]byte) int {
	maxRun := 0

	for _, seq := range sequences {
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
	}

	return maxRun
}

func calculateCpGDepletion(sequences [][]byte) float64 {
	totalCpG := 0
	totalC := 0
	totalG := 0

	for _, seq := range sequences {
		for i := 0; i < len(seq)-1; i++ {
			if seq[i] == 'C' && seq[i+1] == 'G' {
				totalCpG++
			}
			if seq[i] == 'C' {
				totalC++
			}
			if seq[i] == 'G' {
				totalG++
			}
		}
	}

	if totalC == 0 || totalG == 0 {
		return 0.0
	}

	expectedCpG := float64(totalC*totalG) / float64(len(sequences)*50000) // approximate
	if expectedCpG == 0 {
		return 0.0
	}

	return float64(totalCpG) / expectedCpG
}

func loadReferenceSequences(refGenome *reference.SimpleReferenceGenome, numSequences, length int) ([][]byte, error) {
	var sequences [][]byte

	// Load chr1 for reference sequences
	chr1, err := refGenome.LoadChromosome("chr1")
	if err != nil {
		return nil, err
	}

	chrSequence := chr1.GetSequence()
	chrLength := chr1.GetLength()

	for i := 0; i < numSequences; i++ {
		// Sample from different regions of chr1
		start := i * length * 10
		if start >= chrLength {
			start = chrLength - length
		}
		if start < 0 {
			start = 0
		}

		end := start + length
		if end > chrLength {
			end = chrLength
		}

		if end <= start {
			break
		}

		sequence := make([]byte, end-start)
		copy(sequence, chrSequence[start:end])
		sequences = append(sequences, sequence)
	}

	return sequences, nil
}

func calculateJensenShannonDivergence(syntheticSeqs, realSeqs [][]byte, k int) float64 {
	syntheticKmers := extractKmerCounts(syntheticSeqs, k)
	realKmers := extractKmerCounts(realSeqs, k)

	// Normalize to probabilities
	syntheticProbs := normalizeCounts(syntheticKmers)
	realProbs := normalizeCounts(realKmers)

	// Calculate Jensen-Shannon divergence
	return jensenShannonDivergence(syntheticProbs, realProbs)
}

func extractKmerCounts(sequences [][]byte, k int) map[string]int {
	counts := make(map[string]int)

	for _, seq := range sequences {
		for i := 0; i <= len(seq)-k; i++ {
			kmer := string(seq[i : i+k])
			counts[kmer]++
		}
	}

	return counts
}

func normalizeCounts(counts map[string]int) map[string]float64 {
	total := 0
	for _, count := range counts {
		total += count
	}

	probs := make(map[string]float64)
	for kmer, count := range counts {
		probs[kmer] = float64(count) / float64(total)
	}

	return probs
}

func jensenShannonDivergence(p, q map[string]float64) float64 {
	// Get all unique kmers
	allKmers := make(map[string]bool)
	for kmer := range p {
		allKmers[kmer] = true
	}
	for kmer := range q {
		allKmers[kmer] = true
	}

	jsd := 0.0
	for kmer := range allKmers {
		pVal := p[kmer]
		qVal := q[kmer]

		if pVal == 0 && qVal == 0 {
			continue
		}

		m := (pVal + qVal) / 2.0

		if pVal > 0 && m > 0 {
			jsd += pVal * math.Log2(pVal/m)
		}
		if qVal > 0 && m > 0 {
			jsd += qVal * math.Log2(qVal/m)
		}
	}

	return jsd / 2.0
}
