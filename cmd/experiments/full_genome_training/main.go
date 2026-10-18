// Command full_genome_training counts 4-mers in chr1-22, X and Y of the
// reference genome in the hard-coded directory data/validation/reference_genome,
// trains a k=4 adaptive HMM on it with TrainFromReferenceGenome, generates
// three 50 kb sequences and prints complexity, GC content, longest
// homopolymer, CpG ratio and the 4-mer JS divergence against the reference.
// Comparisons with an earlier chr21+chr22 run use hard-coded figures.
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
	log.Printf("=== Full Genome Training for Maximum Complexity ===")
	log.Printf("Training on every chromosome (chr1-22, X, Y) found in the reference directory")
	log.Printf("Target: Sequence complexity > 0.5, JS divergence < 0.01")

	// Step 1: Initialize full genome access
	log.Printf("\n[STEP 1] Initializing full genome access...")
	refGenomeDir := "data/validation/reference_genome"
	refGenome := reference.NewSimpleReferenceGenome(refGenomeDir)

	// List all available chromosomes
	chrNames := refGenome.GetChromosomeNames()
	log.Printf("✓ Available chromosomes: %v", chrNames)

	// Step 2: Create enhanced adaptive model
	log.Printf("\n[STEP 2] Creating enhanced adaptive model...")
	k := 4
	model := authenticity.NewAdaptiveHMMModel(k)
	log.Printf("✓ Model created with k=%d", k)
	log.Printf("✓ State space: %d states", int(math.Pow(4, float64(k-1))))

	// Step 3: Create full-genome trainer
	log.Printf("\n[STEP 3] Creating full-genome trainer...")
	trainer := NewFullGenomeTrainer(model, refGenome)

	// Step 4: Train on full genome
	log.Printf("\n[STEP 4] Training on full human genome...")
	log.Printf("⚠️  This can take a long time and a lot of memory")

	start := time.Now()
	stats, err := trainer.TrainOnFullGenome()
	if err != nil {
		log.Fatalf("❌ Full genome training failed: %v", err)
	}
	trainingTime := time.Since(start)

	log.Printf("✅ FULL GENOME TRAINING COMPLETE!")
	log.Printf("✓ Training time: %v", trainingTime)
	log.Printf("✓ Chromosomes processed: %d", stats.ChromosomesProcessed)
	log.Printf("✓ Total k-mers processed: %d", stats.TotalKmers)
	log.Printf("✓ Unique k-mers observed: %d", stats.UniqueKmers)
	log.Printf("✓ Processing rate: %.0f k-mers/second", float64(stats.TotalKmers)/trainingTime.Seconds())
	log.Printf("✓ K-mer space coverage: %.2f%% of possible %d-mers observed", stats.GenomeCoverage*100, k)

	// Step 5: Test sequence complexity
	log.Printf("\n[STEP 5] Testing sequence complexity with full genome training...")

	sequenceLength := 50000
	numSequences := 3

	var sequences [][]byte
	var totalGenTime time.Duration

	for i := 0; i < numSequences; i++ {
		log.Printf("  Generating sequence %d/%d...", i+1, numSequences)

		result, err := model.GenerateSequence(sequenceLength)
		if err != nil {
			log.Fatalf("❌ Generation failed: %v", err)
		}

		sequences = append(sequences, result.Sequence)
		totalGenTime += result.GenerationTime

		// Calculate complexity for this sequence
		complexity := calculateSequenceComplexity(result.Sequence, k)
		log.Printf("    Sequence %d complexity: %.3f", i+1, complexity)
	}

	avgGenTime := totalGenTime / time.Duration(numSequences)
	log.Printf("✓ Average generation time: %v", avgGenTime)

	// Step 6: Comprehensive metrics
	log.Printf("\n[STEP 6] Comprehensive metrics with full genome training...")

	// 6.1: Sequence Complexity
	avgComplexity := 0.0
	for _, seq := range sequences {
		complexity := calculateSequenceComplexity(seq, k)
		avgComplexity += complexity
	}
	avgComplexity /= float64(len(sequences))

	complexityPass := avgComplexity > 0.5
	log.Printf("Average Sequence Complexity: %.3f (target: > 0.5) %s", avgComplexity, passFailIcon(complexityPass))

	// 6.2: GC Content
	gcContent := calculateGCContent(sequences)
	gcError := math.Abs(gcContent - 0.42)
	gcPass := gcError < 0.05
	log.Printf("GC Content: %.2f%% (target: 42%%, error: %.4f) %s", gcContent*100, gcError, passFailIcon(gcPass))

	// 6.3: Homopolymer Analysis
	maxHomopolymer := findMaxHomopolymerRun(sequences)
	homopolymerPass := maxHomopolymer < 15
	log.Printf("Max Homopolymer: %d (target: < 15) %s", maxHomopolymer, passFailIcon(homopolymerPass))

	// 6.4: CpG Depletion
	cpgDepletion := calculateCpGDepletion(sequences)
	cpgPass := math.Abs(cpgDepletion-0.25) < 0.5
	log.Printf("CpG Depletion: %.3f (target: ~0.25) %s", cpgDepletion, passFailIcon(cpgPass))

	// Step 7: JS Divergence with full genome reference
	log.Printf("\n[STEP 7] JS Divergence with full genome reference...")
	log.Printf("  Extracting reference k-mers from multiple chromosomes...")

	refKmers := extractFullGenomeKmers(refGenome, k, 500000)
	log.Printf("  ✓ Extracted %d reference k-mers from full genome", len(refKmers))

	jsDivergence := calculateJSDivergence(sequences, refKmers, k)
	jsPass := jsDivergence < 0.01
	log.Printf("JS Divergence: %.6f (target: < 0.01) %s", jsDivergence, passFailIcon(jsPass))

	// Step 8: Final Assessment
	log.Printf("\n[STEP 8] Final Assessment with Full Genome Training")

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

	log.Printf("╔════════════════════════════════════════════════════════════════════════════════════════╗")
	log.Printf("║                        FULL GENOME TRAINING RESULTS                                   ║")
	log.Printf("╠════════════════════════════════════════════════════════════════════════════════════════╣")
	log.Printf("║ Chromosomes Loaded:    %-2d                                                            ║", stats.ChromosomesProcessed)
	log.Printf("║ Sequence Complexity:   %.3f > 0.5     %s                                    ║", avgComplexity, passFailPadded(complexityPass))
	log.Printf("║ GC Content Error:      %.4f < 0.05    %s                                    ║", gcError, passFailPadded(gcPass))
	log.Printf("║ Max Homopolymer:       %d < 15       %s                                    ║", maxHomopolymer, passFailPadded(homopolymerPass))
	log.Printf("║ JS Divergence:         %.6f < 0.01 %s                                    ║", jsDivergence, passFailPadded(jsPass))
	log.Printf("╠════════════════════════════════════════════════════════════════════════════════════════╣")
	log.Printf("║ CRITERIA PASSED:       %d/4                                                          ║", totalPassed)
	log.Printf("║ (an earlier chr21+chr22 run passed 1/4; hard-coded, not re-run)                       ║")
	log.Printf("╚════════════════════════════════════════════════════════════════════════════════════════╝")

	// Step 9: Comparison with previous results
	log.Printf("\n[STEP 9] Comparison with Limited Training")
	log.Printf("📊 Training Data Comparison (earlier figures hard-coded):")
	log.Printf("   • Previous: chr21 + chr22 = ~97.5MB")
	log.Printf("   • Current: %d chromosomes loaded", stats.ChromosomesProcessed)

	log.Printf("📈 Complexity Improvement:")
	log.Printf("   • Previous complexity (hard-coded): 0.005")
	log.Printf("   • Current complexity: %.3f", avgComplexity)
	log.Printf("   • Ratio to previous: %.0fx", avgComplexity/0.005)

	log.Printf("🎯 Threshold Summary:")
	if totalPassed == 4 {
		log.Printf("   ✅ All four thresholds met")
	} else if totalPassed >= 2 {
		log.Printf("   ⚠️  %d/4 thresholds met", totalPassed)
	} else {
		log.Printf("   ❌ Fewer than 2/4 thresholds met")
	}

	log.Printf("\n=== Full Genome Training Complete ===")
}

// FullGenomeTrainer handles training on the complete human genome
type FullGenomeTrainer struct {
	model     *authenticity.AdaptiveHMMModel
	refGenome *reference.SimpleReferenceGenome
}

// TrainingStats contains statistics from full genome training
type TrainingStats struct {
	ChromosomesProcessed int
	TotalKmers           uint64
	UniqueKmers          int
	GenomeCoverage       float64
	TrainingTime         time.Duration
}

// NewFullGenomeTrainer creates a new full genome trainer
func NewFullGenomeTrainer(model *authenticity.AdaptiveHMMModel, refGenome *reference.SimpleReferenceGenome) *FullGenomeTrainer {
	return &FullGenomeTrainer{
		model:     model,
		refGenome: refGenome,
	}
}

// TrainOnFullGenome trains the model on the complete human genome
func (t *FullGenomeTrainer) TrainOnFullGenome() (*TrainingStats, error) {
	log.Printf("Starting full genome training...")

	// Define chromosomes to process (all autosomes + X, Y)
	chromosomes := []string{
		"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
		"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
		"chr21", "chr22", "chrX", "chrY",
	}

	stats := &TrainingStats{}
	uniqueKmers := make(map[string]bool)

	start := time.Now()

	// Process each chromosome
	for i, chrName := range chromosomes {
		log.Printf("  Processing %s (%d/%d)...", chrName, i+1, len(chromosomes))

		chr, err := t.refGenome.LoadChromosome(chrName)
		if err != nil {
			log.Printf("    Warning: Could not load %s: %v", chrName, err)
			continue
		}

		sequence := chr.GetSequence()
		log.Printf("    Chromosome %s: %d bases", chrName, len(sequence))

		// Count k-mers for this chromosome
		kmerCount := 0
		k := t.model.K()

		for j := 0; j <= len(sequence)-k; j++ {
			kmer := string(sequence[j : j+k])
			if isValidKmer(kmer) {
				uniqueKmers[kmer] = true
				kmerCount++
			}
		}

		stats.TotalKmers += uint64(kmerCount)
		stats.ChromosomesProcessed++

		log.Printf("    ✓ %s: %d k-mers processed", chrName, kmerCount)
	}

	// Train the model using the standard interface
	err := t.model.TrainFromReferenceGenome(t.refGenome)
	if err != nil {
		return nil, err
	}

	stats.TrainingTime = time.Since(start)
	stats.UniqueKmers = len(uniqueKmers)

	// Calculate theoretical maximum k-mers (4^k)
	maxPossibleKmers := int(math.Pow(4, float64(t.model.K())))
	stats.GenomeCoverage = float64(stats.UniqueKmers) / float64(maxPossibleKmers)

	return stats, nil
}

// Helper functions

func isValidKmer(kmer string) bool {
	for _, base := range kmer {
		if base != 'A' && base != 'T' && base != 'G' && base != 'C' {
			return false
		}
	}
	return true
}

func passFailIcon(pass bool) string {
	if pass {
		return "✅"
	}
	return "❌"
}

func passFailPadded(pass bool) string {
	if pass {
		return "✅ PASS"
	}
	return "❌ FAIL"
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

func extractFullGenomeKmers(refGenome *reference.SimpleReferenceGenome, k int, maxKmers int) []string {
	var kmers []string
	count := 0

	// Sample from multiple chromosomes for better representation
	chromosomes := []string{"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"}

	for _, chrName := range chromosomes {
		if count >= maxKmers {
			break
		}

		chr, err := refGenome.LoadChromosome(chrName)
		if err != nil {
			continue
		}

		sequence := chr.GetSequence()

		// Sample evenly across the chromosome
		sampleStep := len(sequence) / (maxKmers / len(chromosomes))
		if sampleStep < 1 {
			sampleStep = 1
		}

		for i := 0; i <= len(sequence)-k && count < maxKmers; i += sampleStep {
			kmer := string(sequence[i : i+k])
			if isValidKmer(kmer) {
				kmers = append(kmers, kmer)
				count++
			}
		}
	}

	return kmers
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
