// Command quick_validation trains a k=4 adaptive HMM on chr22 from
// data/astrapi_training/human (path hard-coded, relative to the working
// directory), generates 10 kb, and logs training time, k-mer count,
// complexity and GC content, plus a linear extrapolation of training time to
// a 3.2 Gb genome.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"fmt"
	"log"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func main() {
	log.Printf("=== Quick Astrapi Training Validation ===")
	log.Printf("Quick test to validate astrapi server training capability")

	// Test with single chromosome first
	log.Printf("\n[TEST] Single chromosome training validation...")
	refGenomeDir := "data/astrapi_training/human"
	refGenome := reference.NewSimpleReferenceGenome(refGenomeDir)

	// Check chromosomes available
	chromosomes := refGenome.GetChromosomeNames()
	log.Printf("✓ Chromosomes available: %d", len(chromosomes))
	log.Printf("✓ First few: %v", chromosomes[:5])

	// Load chr22 (smallest autosome for quick test)
	log.Printf("\n[TEST] Loading chromosome 22 for quick validation...")
	chr22, err := refGenome.LoadChromosome("chr22")
	if err != nil {
		log.Fatalf("❌ Failed to load chr22: %v", err)
	}

	sequence := chr22.GetSequence()
	log.Printf("✓ chr22 loaded: %d bases", len(sequence))

	// Quick training test with simple model
	log.Printf("\n[TEST] Quick training test...")
	k := 4
	model := authenticity.NewAdaptiveHMMModel(k)

	start := time.Now()

	// Create a simple genome interface for training
	testGenome := &QuickTestGenome{
		chromosomes: map[string][]byte{
			"chr22": sequence,
		},
	}

	err = model.TrainFromReferenceGenome(testGenome)
	if err != nil {
		log.Fatalf("❌ Quick training failed: %v", err)
	}

	trainingTime := time.Since(start)
	log.Printf("✅ Quick training completed in %v", trainingTime)
	log.Printf("✅ K-mers processed: %d", model.TotalKmers())

	// Generate test sequence
	log.Printf("\n[TEST] Generating test sequence...")
	result, err := model.GenerateSequence(10000)
	if err != nil {
		log.Fatalf("❌ Generation failed: %v", err)
	}

	log.Printf("✓ Generated %d bases in %v", len(result.Sequence), result.GenerationTime)

	// Quick metrics
	complexity := calculateSequenceComplexity(result.Sequence, k)
	gcContent := calculateGCContent([][]byte{result.Sequence})

	log.Printf("\n📊 QUICK VALIDATION RESULTS:")
	log.Printf("============================")
	log.Printf("Training Data:    chr22 only (%d bases)", len(sequence))
	log.Printf("Training Time:    %v", trainingTime)
	log.Printf("K-mers:          %d", model.TotalKmers())
	log.Printf("Complexity:      %.4f", complexity)
	log.Printf("GC Content:      %.2f%%", gcContent*100)

	// Performance projections
	log.Printf("\n🚀 ASTRAPI SERVER PERFORMANCE VALIDATION:")
	log.Printf("=========================================")

	basesPerSecond := float64(len(sequence)) / trainingTime.Seconds()
	fullGenomeEstimate := time.Duration(float64(3.2e9)/basesPerSecond) * time.Second

	log.Printf("chr22 Training Rate:   %.1f Mbases/sec", basesPerSecond/1e6)
	log.Printf("chr22 Training Rate:   %.1f bases/sec", basesPerSecond)
	log.Printf("Full Genome Estimate:  %v (linear extrapolation)", fullGenomeEstimate)

	if complexity > 0.01 {
		log.Printf("✅ Training and generation completed (complexity %.4f > 0.01)", complexity)
		log.Printf("✅ Full genome training estimated at %v", fullGenomeEstimate)
	} else {
		log.Printf("⚠ Complexity %.4f is at or below 0.01", complexity)
	}

	log.Printf("\n📈 Linear estimate for a complete genome:")
	log.Printf("- Training time: %v", fullGenomeEstimate)
}

// QuickTestGenome implements GenomeInterface for testing
type QuickTestGenome struct {
	chromosomes map[string][]byte
}

func (g *QuickTestGenome) GetChromosomeNames() []string {
	var names []string
	for name := range g.chromosomes {
		names = append(names, name)
	}
	return names
}

func (g *QuickTestGenome) LoadChromosome(name string) (authenticity.ChromosomeInterface, error) {
	sequence, exists := g.chromosomes[name]
	if !exists {
		return nil, fmt.Errorf("chromosome %s not found", name)
	}

	return &QuickTestChromosome{
		name:     name,
		sequence: sequence,
	}, nil
}

// QuickTestChromosome implements ChromosomeInterface
type QuickTestChromosome struct {
	name     string
	sequence []byte
}

func (c *QuickTestChromosome) GetSequence() []byte {
	return c.sequence
}

func (c *QuickTestChromosome) GetName() string {
	return c.name
}

func (c *QuickTestChromosome) GetLength() int {
	return len(c.sequence)
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
