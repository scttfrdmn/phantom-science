// Command validate_models trains a 4-mer HMMKmerModel on the reference
// genome in data/validation/reference_genome (relative to the current
// directory), generates ten 10 kb sequences and logs their GC content, base
// composition and longest homopolymer against fixed ranges (GC 35-50%,
// homopolymer < 20).
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"log"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func main() {
	log.Printf("=== Phantom Science Model Validation ===")

	// Create reference genome
	refGenome := reference.NewSimpleReferenceGenome("data/validation/reference_genome")

	// Create and train a 4-mer model
	model := authenticity.NewHMMKmerModel(4)

	log.Printf("Training 4-mer model for validation...")
	start := time.Now()
	err := model.TrainFromReferenceGenome(refGenome)
	if err != nil {
		log.Fatalf("Training failed: %v", err)
	}
	trainingTime := time.Since(start)

	log.Printf("Model trained in %v", trainingTime)
	log.Printf("Total k-mers processed: %d", model.TotalKmers())
	log.Printf("State space size: %d", len(model.StateSpace()))

	// Generate multiple sequences and analyze them
	log.Printf("Generating sequences for validation...")

	sequences := make([][]byte, 10)
	var totalGenTime time.Duration

	for i := 0; i < 10; i++ {
		result, err := model.GenerateSequence(10000)
		if err != nil {
			log.Fatalf("Generation failed: %v", err)
		}
		sequences[i] = result.Sequence
		totalGenTime += result.GenerationTime
	}

	log.Printf("Generated 10 sequences of 10,000 bases each")
	log.Printf("Average generation time: %v", totalGenTime/10)

	// Analyze sequence properties
	log.Printf("\n=== Sequence Analysis ===")

	totalGC := 0
	totalBases := 0

	for i, seq := range sequences {
		gcCount := 0
		for _, base := range seq {
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}
		gcContent := float64(gcCount) / float64(len(seq)) * 100
		totalGC += gcCount
		totalBases += len(seq)

		log.Printf("Sequence %d: GC content = %.1f%%", i+1, gcContent)
	}

	avgGC := float64(totalGC) / float64(totalBases) * 100
	log.Printf("Average GC content: %.1f%%", avgGC)

	// Check for biological realism
	log.Printf("\n=== Biological Realism Check ===")

	// Human genome has approximately 40-42% GC content
	if avgGC >= 35.0 && avgGC <= 50.0 {
		log.Print("✓ GC content is within realistic range (35-50%)")
	} else {
		log.Printf("⚠ GC content outside realistic range: %.1f%%", avgGC)
	}

	// Check for base distribution
	baseCounts := make(map[byte]int)
	for _, seq := range sequences {
		for _, base := range seq {
			baseCounts[base]++
		}
	}

	log.Printf("Base distribution:")
	for base, count := range baseCounts {
		percentage := float64(count) / float64(totalBases) * 100
		log.Printf("  %c: %.1f%%", base, percentage)
	}

	// Check for homopolymer runs
	log.Printf("\n=== Homopolymer Analysis ===")

	maxHomopolymer := 0
	for _, seq := range sequences {
		maxRun := findMaxHomopolymerRun(seq)
		if maxRun > maxHomopolymer {
			maxHomopolymer = maxRun
		}
	}

	log.Printf("Maximum homopolymer run: %d", maxHomopolymer)

	if maxHomopolymer < 20 {
		log.Printf("✓ Homopolymer runs are within realistic range")
	} else {
		log.Printf("⚠ Homopolymer runs are too long: %d", maxHomopolymer)
	}

	log.Printf("\n=== Validation Complete ===")
	log.Printf("Checks above cover GC content and homopolymer length only")
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
