// Command final_authenticity_test trains the original, enhanced and
// adaptive k=4 HMMs on the reference genome in the hard-coded directory
// data/validation/reference_genome, generates five 100 kb sequences from
// each, and compares them with the first 500,000 4-mers of chr21/chr22 by
// JS and KL divergence, GC content, complexity, longest homopolymer and CpG
// ratio. The "biological realism" score is an ad hoc 0-100 sum of those
// metrics, and "production ready" (ModelMetrics.ProductionReady) only means
// all four pass thresholds were met.
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
	log.Printf("=== Final Authenticity Model Comparison ===")
	log.Printf("Target: JS divergence < 0.01")

	// Load reference genome
	refGenomeDir := "data/validation/reference_genome"
	refGenome := reference.NewSimpleReferenceGenome(refGenomeDir)

	// Test 4-mer models (all three versions)
	k := 4

	log.Printf("\n=== Training All Models ===")

	// Original model
	log.Printf("Training Original HMM Model...")
	originalModel := authenticity.NewHMMKmerModel(k)
	start := time.Now()
	err := originalModel.TrainFromReferenceGenome(refGenome)
	if err != nil {
		log.Fatalf("Original model training failed: %v", err)
	}
	log.Printf("Original model trained in: %v", time.Since(start))

	// Enhanced model
	log.Printf("Training Enhanced HMM Model...")
	enhancedModel := authenticity.NewEnhancedHMMKmerModel(k, "CEU")
	start = time.Now()
	err = enhancedModel.TrainFromReferenceGenome(refGenome)
	if err != nil {
		log.Fatalf("Enhanced model training failed: %v", err)
	}
	log.Printf("Enhanced model trained in: %v", time.Since(start))

	// Adaptive model
	log.Printf("Training Adaptive HMM Model...")
	adaptiveModel := authenticity.NewAdaptiveHMMModel(k)
	start = time.Now()
	err = adaptiveModel.TrainFromReferenceGenome(refGenome)
	if err != nil {
		log.Fatalf("Adaptive model training failed: %v", err)
	}
	log.Printf("Adaptive model trained in: %v", time.Since(start))

	// Generate sequences from all models
	sequenceLength := 100000
	numSequences := 5

	log.Printf("\n=== Generating Sequences ===")

	// Original model sequences
	originalSequences := generateSequences(originalModel, sequenceLength, numSequences, "Original")

	// Enhanced model sequences
	enhancedSequences := generateSequences(enhancedModel, sequenceLength, numSequences, "Enhanced")

	// Adaptive model sequences
	adaptiveSequences := generateSequences(adaptiveModel, sequenceLength, numSequences, "Adaptive")

	// Extract reference k-mers for JS divergence calculation
	log.Printf("\n=== Extracting Reference K-mers ===")
	refKmers := extractReferenceKmers(refGenome, k, 500000)
	log.Printf("Extracted %d reference k-mers", len(refKmers))

	// Calculate comprehensive metrics for all models
	log.Printf("\n=== Comprehensive Metrics Comparison ===")

	originalMetrics := calculateComprehensiveMetrics(originalSequences, refKmers, k, "Original")
	enhancedMetrics := calculateComprehensiveMetrics(enhancedSequences, refKmers, k, "Enhanced")
	adaptiveMetrics := calculateComprehensiveMetrics(adaptiveSequences, refKmers, k, "Adaptive")

	// Print detailed comparison
	log.Printf("\n=== Detailed Results ===")
	printDetailedComparison(originalMetrics, enhancedMetrics, adaptiveMetrics)

	// Determine best model and whether it meets every threshold
	log.Printf("\n=== Threshold Assessment ===")
	assessProductionReadiness(originalMetrics, enhancedMetrics, adaptiveMetrics)

	// Generate final authenticity report
	log.Printf("\n=== Final Authenticity Report ===")
	generateFinalReport(originalMetrics, enhancedMetrics, adaptiveMetrics)
}

type ModelMetrics struct {
	ModelName          string
	JSdivergence       float64
	KLdivergence       float64
	GCContent          float64
	GCContentTarget    float64
	GCContentError     float64
	SequenceComplexity float64
	MaxHomopolymer     int
	CpGDepletion       float64
	GenerationTime     time.Duration
	BiologicalRealism  float64 // Overall biological realism score
	ProductionReady    bool
}

func generateSequences(model interface{}, length, count int, modelName string) [][]byte {
	var sequences [][]byte

	log.Printf("Generating %d sequences of %d bases each for %s model...", count, length, modelName)

	for i := 0; i < count; i++ {
		var result *authenticity.GenerationResult
		var err error

		switch m := model.(type) {
		case *authenticity.HMMKmerModel:
			result, err = m.GenerateSequence(length)
		case *authenticity.EnhancedHMMKmerModel:
			result, err = m.GenerateSequence(length)
		case *authenticity.AdaptiveHMMModel:
			result, err = m.GenerateSequence(length)
		default:
			log.Fatalf("Unknown model type")
		}

		if err != nil {
			log.Fatalf("%s generation failed: %v", modelName, err)
		}

		sequences = append(sequences, result.Sequence)
	}

	return sequences
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
			log.Printf("Warning: Could not load chromosome %s: %v", chrName, err)
			continue
		}

		sequence := chr.GetSequence()
		log.Printf("Processing chromosome %s (%d bases)", chrName, len(sequence))

		for i := 0; i <= len(sequence)-k && count < maxKmers; i++ {
			kmer := string(sequence[i : i+k])
			// Only include valid k-mers
			if isValidKmer(kmer) {
				kmers = append(kmers, kmer)
				count++
			}
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

func calculateComprehensiveMetrics(sequences [][]byte, refKmers []string, k int, modelName string) ModelMetrics {
	log.Printf("Calculating comprehensive metrics for %s model...", modelName)

	var metrics ModelMetrics
	metrics.ModelName = modelName

	// Calculate basic biological metrics
	totalBases := 0
	gcCount := 0
	maxHomopolymer := 0

	for _, seq := range sequences {
		for _, base := range seq {
			totalBases++
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}

		// Max homopolymer
		maxRun := findMaxHomopolymerRun(seq)
		if maxRun > maxHomopolymer {
			maxHomopolymer = maxRun
		}
	}

	metrics.GCContent = float64(gcCount) / float64(totalBases)
	metrics.GCContentTarget = 0.42
	metrics.GCContentError = math.Abs(metrics.GCContent - metrics.GCContentTarget)
	metrics.MaxHomopolymer = maxHomopolymer

	// Calculate sequence complexity
	metrics.SequenceComplexity = calculateSequenceComplexity(sequences[0], k)

	// Calculate CpG depletion
	metrics.CpGDepletion = calculateCpGDepletion(sequences)

	// Calculate JS divergence
	metrics.JSdivergence = calculateJSDivergence(sequences, refKmers, k)
	metrics.KLdivergence = calculateKLDivergence(sequences, refKmers, k)

	// Calculate overall biological realism score
	metrics.BiologicalRealism = calculateBiologicalRealismScore(metrics)

	// ProductionReady means all four thresholds are met; nothing more
	metrics.ProductionReady = metrics.JSdivergence < 0.01 &&
		metrics.GCContentError < 0.05 &&
		metrics.SequenceComplexity > 0.5 &&
		metrics.MaxHomopolymer < 15

	return metrics
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

func calculateKLDivergence(sequences [][]byte, refKmers []string, k int) float64 {
	// Similar to JS divergence but calculate KL divergence
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

	refKmerMap := make(map[string]int)
	for _, kmer := range refKmers {
		refKmerMap[kmer]++
	}
	refTotal := len(refKmers)

	// Calculate KL divergence
	kl := 0.0
	for kmer, count := range genKmers {
		pVal := float64(count) / float64(genTotal)
		qVal := float64(refKmerMap[kmer]) / float64(refTotal)

		if pVal > 0 && qVal > 0 {
			kl += pVal * math.Log2(pVal/qVal)
		}
	}

	return kl
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

func calculateBiologicalRealismScore(metrics ModelMetrics) float64 {
	score := 0.0

	// GC content score (0-25 points)
	gcScore := 25.0 * (1.0 - metrics.GCContentError/0.1)
	if gcScore < 0 {
		gcScore = 0
	}
	score += gcScore

	// Sequence complexity score (0-25 points)
	complexityScore := 25.0 * metrics.SequenceComplexity
	if complexityScore > 25 {
		complexityScore = 25
	}
	score += complexityScore

	// Homopolymer score (0-25 points)
	homopolymerScore := 25.0 * (1.0 - float64(metrics.MaxHomopolymer)/20.0)
	if homopolymerScore < 0 {
		homopolymerScore = 0
	}
	score += homopolymerScore

	// CpG depletion score (0-25 points)
	cpgTarget := 0.25
	cpgError := math.Abs(metrics.CpGDepletion - cpgTarget)
	cpgScore := 25.0 * (1.0 - cpgError/1.0)
	if cpgScore < 0 {
		cpgScore = 0
	}
	score += cpgScore

	return score
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

func printDetailedComparison(original, enhanced, adaptive ModelMetrics) {
	log.Printf("┌─────────────────────────────────────────────────────────────────────────────────────┐")
	log.Printf("│                               DETAILED COMPARISON                                   │")
	log.Printf("├─────────────────────────────────────────────────────────────────────────────────────┤")
	log.Printf("│ Metric                  │ Original    │ Enhanced    │ Adaptive    │ Target      │")
	log.Printf("├─────────────────────────────────────────────────────────────────────────────────────┤")
	log.Printf("│ JS Divergence           │ %.6f  │ %.6f  │ %.6f  │ < 0.010000  │",
		original.JSdivergence, enhanced.JSdivergence, adaptive.JSdivergence)
	log.Printf("│ KL Divergence           │ %.6f  │ %.6f  │ %.6f  │ < 0.100000  │",
		original.KLdivergence, enhanced.KLdivergence, adaptive.KLdivergence)
	log.Printf("│ GC Content              │ %.2f%%     │ %.2f%%     │ %.2f%%     │ 42.00%%     │",
		original.GCContent*100, enhanced.GCContent*100, adaptive.GCContent*100)
	log.Printf("│ GC Content Error        │ %.4f    │ %.4f    │ %.4f    │ < 0.0500    │",
		original.GCContentError, enhanced.GCContentError, adaptive.GCContentError)
	log.Printf("│ Sequence Complexity     │ %.3f     │ %.3f     │ %.3f     │ > 0.500     │",
		original.SequenceComplexity, enhanced.SequenceComplexity, adaptive.SequenceComplexity)
	log.Printf("│ Max Homopolymer         │ %d         │ %d         │ %d         │ < 15        │",
		original.MaxHomopolymer, enhanced.MaxHomopolymer, adaptive.MaxHomopolymer)
	log.Printf("│ CpG Depletion           │ %.3f     │ %.3f     │ %.3f     │ ~0.250      │",
		original.CpGDepletion, enhanced.CpGDepletion, adaptive.CpGDepletion)
	log.Printf("│ Biological Realism      │ %.1f/100   │ %.1f/100   │ %.1f/100   │ > 80.0      │",
		original.BiologicalRealism, enhanced.BiologicalRealism, adaptive.BiologicalRealism)
	log.Printf("│ All Thresholds Met      │ %v       │ %v       │ %v       │ true        │",
		original.ProductionReady, enhanced.ProductionReady, adaptive.ProductionReady)
	log.Printf("└─────────────────────────────────────────────────────────────────────────────────────┘")
}

func assessProductionReadiness(original, enhanced, adaptive ModelMetrics) {
	models := []ModelMetrics{original, enhanced, adaptive}

	var bestModel *ModelMetrics
	bestScore := -1.0

	for i := range models {
		if models[i].BiologicalRealism > bestScore {
			bestScore = models[i].BiologicalRealism
			bestModel = &models[i]
		}
	}

	log.Printf("Best Model: %s", bestModel.ModelName)
	log.Printf("Best Score: %.1f/100", bestModel.BiologicalRealism)

	if bestModel.ProductionReady {
		log.Printf("✅ ALL THRESHOLDS MET - %s model meets every criterion below", bestModel.ModelName)
		log.Printf("   • JS Divergence: %.6f < 0.01 ✓", bestModel.JSdivergence)
		log.Printf("   • GC Content: %.2f%% (error: %.4f) ✓", bestModel.GCContent*100, bestModel.GCContentError)
		log.Printf("   • Sequence Complexity: %.3f > 0.5 ✓", bestModel.SequenceComplexity)
		log.Printf("   • Max Homopolymer: %d < 15 ✓", bestModel.MaxHomopolymer)
	} else {
		log.Printf("⚠️  NEEDS IMPROVEMENT - %s model does not meet all criteria", bestModel.ModelName)

		if bestModel.JSdivergence >= 0.01 {
			log.Printf("   • JS Divergence: %.6f ≥ 0.01 ✗", bestModel.JSdivergence)
		}
		if bestModel.GCContentError >= 0.05 {
			log.Printf("   • GC Content Error: %.4f ≥ 0.05 ✗", bestModel.GCContentError)
		}
		if bestModel.SequenceComplexity <= 0.5 {
			log.Printf("   • Sequence Complexity: %.3f ≤ 0.5 ✗", bestModel.SequenceComplexity)
		}
		if bestModel.MaxHomopolymer >= 15 {
			log.Printf("   • Max Homopolymer: %d ≥ 15 ✗", bestModel.MaxHomopolymer)
		}
	}
}

func generateFinalReport(original, enhanced, adaptive ModelMetrics) {
	log.Printf("=== PHANTOM SCIENCE AUTHENTICITY ENGINE FINAL REPORT ===")
	log.Printf("")
	log.Printf("OBJECTIVE: Compare the 4-mer distribution and sequence statistics of generated sequences with the reference genome")
	log.Printf("SUCCESS CRITERIA: Jensen-Shannon divergence < 0.01")
	log.Printf("")
	log.Printf("MODELS TESTED:")
	log.Printf("1. Original HMM Model - Basic k-mer HMM implementation")
	log.Printf("2. Enhanced HMM Model - Biological realism constraints")
	log.Printf("3. Adaptive HMM Model - Dynamic complexity control")
	log.Printf("")
	log.Printf("RESULTS SUMMARY:")

	models := []ModelMetrics{original, enhanced, adaptive}
	for _, model := range models {
		log.Printf("")
		log.Printf("%s Model:", model.ModelName)
		log.Printf("  JS Divergence: %.6f", model.JSdivergence)
		log.Printf("  Biological Realism: %.1f/100", model.BiologicalRealism)
		log.Printf("  All Thresholds Met: %v", model.ProductionReady)

		if model.ProductionReady {
			log.Printf("  ✅ Meets all thresholds")
		} else {
			log.Printf("  ⚠️  Needs further improvement")
		}
	}

	log.Printf("")
	log.Printf("RECOMMENDATIONS:")

	bestModel := adaptive
	if adaptive.BiologicalRealism < enhanced.BiologicalRealism {
		bestModel = enhanced
	}
	if enhanced.BiologicalRealism < original.BiologicalRealism {
		bestModel = original
	}

	if bestModel.ProductionReady {
		log.Printf("• Use the %s model for further testing", bestModel.ModelName)
		log.Printf("• Begin expert blind testing phase")
		log.Printf("• Integrate with existing sequencer simulation")
	} else {
		log.Printf("• Continue development on %s model", bestModel.ModelName)
		log.Printf("• Focus on sequence complexity improvement")
		log.Printf("• Consider additional training data")
	}

	log.Printf("")
	log.Printf("NEXT STEPS:")
	log.Printf("1. Expert blind testing with real bioinformaticians")
	log.Printf("2. Performance benchmarking and optimization")
	log.Printf("3. Integration with existing pipeline")
}
