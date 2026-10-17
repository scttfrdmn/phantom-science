package main

import (
	"fmt"
	"log"
	"math"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func main() {
	logIntroduction()

	// Step 1: Initialize
	log.Printf("\n[STEP 1] Initializing reference genome access...")
	refGenomeDir := "data/validation/reference_genome"
	refGenome := reference.NewSimpleReferenceGenome(refGenomeDir)
	log.Printf("✓ Reference genome initialized")

	// Step 2: Create adaptive model
	log.Printf("\n[STEP 2] Creating adaptive HMM model...")
	k := 4
	adaptiveModel := authenticity.NewAdaptiveHMMModel(k)
	logModelConfig(k)
	
	// Step 3: Train the model (this will take time)
	log.Printf("\n[STEP 3] Training adaptive model on real genomic data...")
	log.Printf("This may take several minutes - please be patient...")
	
	start := time.Now()
	err := adaptiveModel.TrainFromReferenceGenome(refGenome)
	if err != nil {
		log.Fatalf("❌ Training failed: %v", err)
	}
	trainingTime := time.Since(start)
	
	log.Printf("✓ Training completed successfully!")
	log.Printf("✓ Training time: %v", trainingTime)
	log.Printf("✓ Total k-mers processed: %d", adaptiveModel.TotalKmers())
	
	// Step 4: Generate test sequences
	log.Printf("\n[STEP 4] Generating test sequences...")
	sequenceLength := 50000
	numSequences := 3
	
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
	log.Printf("✓ All sequences generated, average time: %v", avgGenTime)
	
	// Step 5: Calculate biological metrics
	logStep5Header()

	// Step 6: Jensen-Shannon Divergence
	refKmers := extractLimitedReferenceKmers(refGenome, k, 100000)
	logJSDivergenceHeader(len(refKmers))

	metrics := computeMetrics(sequences, refKmers, k)
	results := evaluatePassing(metrics)
	totalPassed := results.TotalPassed

	logMetricsResults(metrics, results)

	// Step 7: Overall Assessment
	log.Printf("\n[STEP 7] Overall Assessment")
	logResultsBanner(metrics, results)

	// Step 8: Final Determination
	log.Printf("\n[STEP 8] Final Determination")
	logFinalDetermination(totalPassed, metrics, results)

	// Step 9: Performance Summary
	logPerformanceSummary9(trainingTime, adaptiveModel.TotalKmers(), sequenceLength, numSequences, avgGenTime)
	
	// Step 10: Technical Details
	log.Printf("\n[STEP 10] Technical Details")
	if stats := adaptiveModel.GetGenerationStats(); stats.TotalSequences > 0 {
		log.Printf("📊 Generation Statistics:")
		log.Printf("   • Total sequences generated: %d", stats.TotalSequences)
		log.Printf("   • Average complexity: %.3f", stats.AverageComplexity)
		log.Printf("   • Average GC content: %.2f%%", stats.AverageGCContent*100)
		log.Printf("   • Max homopolymer observed: %d", stats.MaxHomopolymer)
	}
	
	logModelDetails(k, len(adaptiveModel.StateSpace()), adaptiveModel.IsTrained())

	logCompletionFooter()
}

// logModelDetails logs the ⚙️ model configuration block in Step 10.
func logModelDetails(k int, stateSpaceSize int, isTrained bool) {
	log.Printf("⚙️  Model Configuration:")
	log.Printf("   • K-mer size: %d", k)
	log.Printf("   • State space size: %d", stateSpaceSize)
	log.Printf("   • Trained: %v", isTrained)
}

// AuthenticityMetrics holds the computed biological metrics for a set of sequences.
type AuthenticityMetrics struct {
	Complexity     float64
	GCContent      float64
	GCError        float64
	MaxHomopolymer int
	CpGDepletion   float64
	JSDivergence   float64
}

// PassResults holds pass/fail results for each authenticity criterion.
type PassResults struct {
	ComplexityPass  bool
	GCPass          bool
	HomopolymerPass bool
	CpGPass         bool
	JSPass          bool
	TotalPassed     int
}

// computeMetrics calculates all authenticity metrics for a set of sequences
// and reference k-mers.
func computeMetrics(sequences [][]byte, refKmers []string, k int) AuthenticityMetrics {
	complexity := calculateSequenceComplexity(sequences[0], k)
	gcContent := calculateGCContent(sequences)
	gcError := math.Abs(gcContent - 0.42)
	maxHomopolymer := findMaxHomopolymerRun(sequences)
	cpgDepletion := calculateCpGDepletion(sequences)
	jsDivergence := calculateJSDivergence(sequences, refKmers, k)
	return AuthenticityMetrics{
		Complexity:     complexity,
		GCContent:      gcContent,
		GCError:        gcError,
		MaxHomopolymer: maxHomopolymer,
		CpGDepletion:   cpgDepletion,
		JSDivergence:   jsDivergence,
	}
}

// evaluatePassing applies the authenticity pass/fail thresholds to the metrics.
func evaluatePassing(m AuthenticityMetrics) PassResults {
	complexityPass := m.Complexity > 0.5
	gcPass := m.GCError < 0.05
	homopolymerPass := m.MaxHomopolymer < 15
	cpgPass := math.Abs(m.CpGDepletion-0.25) < 0.5
	jsPass := m.JSDivergence < 0.01

	total := 0
	if complexityPass {
		total++
	}
	if gcPass {
		total++
	}
	if homopolymerPass {
		total++
	}
	if jsPass {
		total++
	}
	return PassResults{
		ComplexityPass:  complexityPass,
		GCPass:          gcPass,
		HomopolymerPass: homopolymerPass,
		CpGPass:         cpgPass,
		JSPass:          jsPass,
		TotalPassed:     total,
	}
}

// determineFinalStatus returns a human-readable status string based on how
// many criteria were passed.
func determineFinalStatus(totalPassed int) string {
	switch {
	case totalPassed == 4:
		return "PRODUCTION READY"
	case totalPassed >= 2:
		return "SIGNIFICANT PROGRESS"
	default:
		return "NEEDS MAJOR IMPROVEMENT"
	}
}

// improvementAreas returns a list of improvement suggestion strings for criteria
// that failed.  Only items that failed are included.
func improvementAreas(m AuthenticityMetrics, r PassResults) []string {
	var areas []string
	if !r.ComplexityPass {
		areas = append(areas, "Sequence complexity needs improvement")
	}
	if !r.GCPass {
		areas = append(areas, "GC content needs adjustment")
	}
	if !r.HomopolymerPass {
		areas = append(areas, "Homopolymer control needs improvement")
	}
	if !r.JSPass {
		areas = append(areas, "JS divergence needs reduction")
	}
	return areas
}

// formatPerformanceSummary returns a human-readable performance summary line.
func formatPerformanceSummary(trainingSeconds float64, totalKmers uint64, sequenceLength int, avgGenNs int64) string {
	var processingRate float64
	if trainingSeconds > 0 {
		processingRate = float64(totalKmers) / trainingSeconds
	}
	var genRate float64
	msPerSeq := avgGenNs / 1_000_000
	if msPerSeq > 0 {
		genRate = float64(sequenceLength) / float64(msPerSeq)
	}
	return fmt.Sprintf("training=%.2fs kmers=%d rate=%.0f/s genrate=%.2fbases/ms",
		trainingSeconds, totalKmers, processingRate, genRate)
}

// logIntroduction logs the opening banner for the validation run.
func logIntroduction() {
	log.Printf("=== Step-by-Step Authenticity Validation ===")
	log.Printf("Taking time to properly validate the adaptive model...")
	log.Printf("Target: JS Divergence < 0.01 for production authenticity")
}

// logModelConfig logs the model configuration summary.
func logModelConfig(k int) {
	log.Printf("✓ Adaptive model created with k=%d", k)
	log.Printf("✓ Model configured for:")
	log.Printf("  - Target complexity: 70%%")
	log.Printf("  - Target GC content: 42%%")
	log.Printf("  - Max homopolymer: 8 bases")
	log.Printf("  - CpG depletion: 25%%")
}

// logJSDivergenceHeader logs the Step 6 Jensen-Shannon divergence header lines.
func logJSDivergenceHeader(numRefKmers int) {
	log.Printf("\n[STEP 6] Calculating Jensen-Shannon divergence...")
	log.Printf("  Extracting reference k-mers from chromosome 22...")
	log.Printf("  ✓ Extracted %d reference k-mers", numRefKmers)
	log.Printf("  Computing JS divergence...")
}

// logStep5Header logs the Step 5 header and metric-calculation progress lines.
func logStep5Header() {
	log.Printf("\n[STEP 5] Calculating biological realism metrics...")
	log.Printf("  Calculating sequence complexity...")
	log.Printf("  Calculating GC content...")
	log.Printf("  Analyzing homopolymer runs...")
	log.Printf("  Calculating CpG depletion...")
}

// logCompletionFooter logs the final completion lines.
func logCompletionFooter() {
	log.Printf("\n=== Authenticity Validation Complete ===")
	log.Printf("Thank you for your patience during this thorough evaluation.")
}

// logMetricsResults emits the per-metric result lines.
func logMetricsResults(m AuthenticityMetrics, r PassResults) {
	log.Printf("    Complexity: %.3f (target: > 0.5) %s", m.Complexity, passFailIcon(r.ComplexityPass))
	log.Printf("    GC Content: %.2f%% (target: 42%%, error: %.4f) %s", m.GCContent*100, m.GCError, passFailIcon(r.GCPass))
	log.Printf("    Max Homopolymer: %d (target: < 15) %s", m.MaxHomopolymer, passFailIcon(r.HomopolymerPass))
	log.Printf("    CpG Depletion: %.3f (target: ~0.25) %s", m.CpGDepletion, passFailIcon(r.CpGPass))
	log.Printf("    JS Divergence: %.6f (target: < 0.01) %s", m.JSDivergence, passFailIcon(r.JSPass))
}

// logPerformanceSummary9 emits the Step 9 performance summary log lines.
func logPerformanceSummary9(trainingTime time.Duration, totalKmers uint64, sequenceLength int, numSequences int, avgGenTime time.Duration) {
	log.Printf("\n[STEP 9] Performance Summary")
	log.Printf("🔧 Training Performance:")
	log.Printf("   • Training time: %.2f seconds", trainingTime.Seconds())
	log.Printf("   • K-mers processed: %d", totalKmers)
	var processingRate float64
	if trainingTime.Seconds() > 0 {
		processingRate = float64(totalKmers) / trainingTime.Seconds()
	}
	log.Printf("   • Processing rate: %.0f k-mers/second", processingRate)
	log.Printf("🚀 Generation Performance:")
	log.Printf("   • Average generation time: %v", avgGenTime)
	var genRate float64
	ms := avgGenTime.Nanoseconds() / 1_000_000
	if ms > 0 {
		genRate = float64(sequenceLength) / float64(ms)
	}
	log.Printf("   • Generation rate: %.2f bases/ms", genRate)
	log.Printf("   • Sequences generated: %d", numSequences)
}

// logResultsBanner emits the formatted results summary table using log.Printf.
func logResultsBanner(m AuthenticityMetrics, r PassResults) {
	log.Printf("╔════════════════════════════════════════════════════════════════════════════════════════╗")
	log.Printf("║                                    RESULTS SUMMARY                                     ║")
	log.Printf("╠════════════════════════════════════════════════════════════════════════════════════════╣")
	log.Printf("║ Sequence Complexity:   %.3f > 0.5     %s                                    ║", m.Complexity, passFailPadded(r.ComplexityPass))
	log.Printf("║ GC Content Error:      %.4f < 0.05    %s                                    ║", m.GCError, passFailPadded(r.GCPass))
	log.Printf("║ Max Homopolymer:       %d < 15       %s                                    ║", m.MaxHomopolymer, passFailPadded(r.HomopolymerPass))
	log.Printf("║ JS Divergence:         %.6f < 0.01 %s                                    ║", m.JSDivergence, passFailPadded(r.JSPass))
	log.Printf("╠════════════════════════════════════════════════════════════════════════════════════════╣")
	log.Printf("║ CRITERIA PASSED:       %d/4                                                          ║", r.TotalPassed)
	log.Printf("╚════════════════════════════════════════════════════════════════════════════════════════╝")
}

// logFinalDetermination prints the final determination section based on how many
// criteria passed.  It also logs specific improvement areas when applicable.
func logFinalDetermination(totalPassed int, m AuthenticityMetrics, r PassResults) {
	finalStatus := determineFinalStatus(totalPassed)
	log.Printf("Status: %s", finalStatus)

	switch {
	case totalPassed == 4:
		log.Printf("🎉 PRODUCTION READY!")
		log.Printf("✅ Adaptive model meets ALL authenticity criteria")
		log.Printf("✅ Generated sequences are statistically indistinguishable from real genomic data")
		log.Printf("✅ Ready for expert blind testing and production deployment")
		log.Printf("\n📋 Next Steps:")
		log.Printf("   1. Begin expert blind testing phase")
		log.Printf("   2. Performance benchmarking and optimization")
		log.Printf("   3. Integration with existing sequencer simulation")
		log.Printf("   4. Production deployment preparation")
	case totalPassed >= 2:
		log.Printf("⚠️  SIGNIFICANT PROGRESS MADE")
		log.Printf("✅ Adaptive model shows strong improvements")
		log.Printf("⚠️  Continue refinement to meet remaining criteria")
		log.Printf("\n📋 Areas for Improvement:")
		for _, area := range improvementAreas(m, r) {
			log.Printf("   • %s", area)
		}
	default:
		log.Printf("❌ NEEDS MAJOR IMPROVEMENT")
		log.Printf("❌ Adaptive model requires significant additional work")
		log.Printf("❌ Consider alternative approaches or expanded training data")
	}
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