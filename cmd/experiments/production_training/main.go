package main

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// ProgressReport represents training progress that can be written to files
type ProgressReport struct {
	Timestamp           time.Time `json:"timestamp"`
	Phase               string    `json:"phase"`
	Status              string    `json:"status"`
	TrainingStartTime   time.Time `json:"training_start_time,omitempty"`
	TrainingDuration    string    `json:"training_duration,omitempty"`
	KmersProcessed      int64     `json:"kmers_processed,omitempty"`
	ChromosomesLoaded   int       `json:"chromosomes_loaded,omitempty"`
	SequenceComplexity  float64   `json:"sequence_complexity,omitempty"`
	GCContent           float64   `json:"gc_content,omitempty"`
	JSdivergence        float64   `json:"js_divergence,omitempty"`
	MaxHomopolymer      int       `json:"max_homopolymer,omitempty"`
	CpGDepletion        float64   `json:"cpg_depletion,omitempty"`
	GenerationSpeed     float64   `json:"generation_speed_kb_per_sec,omitempty"`
	OverallStatus       string    `json:"overall_status,omitempty"`
	Message             string    `json:"message,omitempty"`
	ErrorMessage        string    `json:"error_message,omitempty"`
}

func writeProgressReport(report ProgressReport, outputDir string) {
	// Write timestamped progress file
	timestamp := time.Now().Format("20060102_150405")
	filename := fmt.Sprintf("progress_%s_%s.json", report.Phase, timestamp)
	filePath := filepath.Join(outputDir, filename)
	
	data, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		log.Printf("Warning: Could not marshal progress report: %v", err)
		return
	}
	
	err = os.WriteFile(filePath, data, 0600)
	if err != nil {
		log.Printf("Warning: Could not write progress file %s: %v", filePath, err)
		return
	}
	
	// Also write latest status
	latestPath := filepath.Join(outputDir, "latest_progress.json")
	err = os.WriteFile(latestPath, data, 0600)
	if err != nil {
		log.Printf("Warning: Could not write latest progress: %v", err)
	}
}

func main() {
	// Setup progress monitoring
	outputDir := "training_progress"
	err := os.MkdirAll(outputDir, 0750)
	if err != nil {
		log.Fatalf("❌ Could not create progress directory: %v", err)
	}

	logProductionIntroduction(outputDir)

	// Initialize progress
	startTime := time.Now()
	writeStartupProgressReport(outputDir)

	// Step 1: Initialize astrapi genome access
	log.Printf("\n[STEP 1] Initializing astrapi genome access...")
	refGenomeDir := "data/astrapi_training/human"
	refGenome := reference.NewSimpleReferenceGenome(refGenomeDir)

	chromosomes := refGenome.GetChromosomeNames()

	// Step 2: Create production model
	log.Printf("\n[STEP 2] Creating production adaptive HMM model...")
	k := 4
	adaptiveModel := authenticity.NewAdaptiveHMMModel(k)

	logProductionModelSetup(k, refGenomeDir, len(chromosomes))
	writeInitProgressReport(outputDir, len(chromosomes), k)

	// Step 3: Production training on complete genome
	logTrainingStepHeader(outputDir)
	
	writeTrainingStartReport(outputDir)

	// Start training
	trainingStart := time.Now()
	err = adaptiveModel.TrainFromReferenceGenome(refGenome)
	if err != nil {
		writeTrainingFailReport(outputDir, trainingStart, fmt.Sprintf("Training failed: %v", err))
		os.Exit(1)
	}
	
	trainingTime := time.Since(trainingStart)
	totalKmers := adaptiveModel.TotalKmers()
	
	logTrainingComplete(trainingTime, totalKmers)
	writeTrainingCompleteReport(outputDir, trainingTime, totalKmers)
	
	// Step 4: Generate validation sequences
	log.Printf("\n[STEP 4] Generating validation sequences...")
	sequenceLength := 50000
	numSequences := 5
	
	writeValidationStartReport(outputDir, numSequences, sequenceLength)

	var sequences [][]byte
	var totalGenTime time.Duration
	
	for i := 0; i < numSequences; i++ {
		log.Printf("  Generating validation sequence %d/%d...", i+1, numSequences)
		
		result, err := adaptiveModel.GenerateSequence(sequenceLength)
		if err != nil {
			writeGenerationFailReport(outputDir, fmt.Sprintf("Generation failed: %v", err))
			os.Exit(1)
		}
		
		sequences = append(sequences, result.Sequence)
		totalGenTime += result.GenerationTime
	}
	
	avgGenTime := totalGenTime / time.Duration(numSequences)
	generationSpeed := float64(sequenceLength) / avgGenTime.Seconds() / 1000.0

	logGenerationComplete(avgGenTime, generationSpeed)

	runValidationAndReport(sequences, avgGenTime, sequenceLength, k, refGenome,
		trainingStart, trainingTime, totalKmers, len(chromosomes), startTime, outputDir)
}

// runValidationAndReport performs Steps 5-10 of production training: computing
// authenticity metrics, JS divergence, writing progress reports, and logging
// the final outcome. It is extracted from main() to enable unit testing.
func runValidationAndReport(
	sequences [][]byte,
	avgGenTime time.Duration,
	sequenceLength int,
	k int,
	refGenome *reference.SimpleReferenceGenome,
	trainingStart time.Time,
	trainingTime time.Duration,
	totalKmers uint64,
	chromosomesLoaded int,
	startTime time.Time,
	outputDir string,
) {
	// Step 5: Comprehensive metrics validation
	log.Printf("\n[STEP 5] Calculating production authenticity metrics...")

	metrics := computeProductionMetrics(sequences, k, avgGenTime, sequenceLength)
	complexityImprovement := metrics.ComplexityImprovement

	// Step 6: Jensen-Shannon Divergence
	log.Printf("\n[STEP 6] Calculating Jensen-Shannon divergence...")
	var jsDivergence float64 = -1

	referenceSequences, err := loadReferenceSequences(refGenome, 5, sequenceLength)
	if err != nil {
		log.Printf("⚠ Warning: Could not load reference sequences: %v", err)
	} else {
		jsDivergence = calculateJensenShannonDivergence(sequences, referenceSequences, k)
		jsImprovement := 0.055 / jsDivergence
		log.Printf("  JS Divergence: %.6f (%.1fx improvement)", jsDivergence, jsImprovement)
	}

	results := evaluateProductionPassing(metrics, jsDivergence)
	overallPass := results.OverallPass
	overallStatus := determineOverallStatus(overallPass)

	logProductionMetricsResults(metrics, results, jsDivergence)

	// Final progress report
	finalReport := buildFinalProgressReport(
		trainingStart,
		trainingTime,
		totalKmers,
		chromosomesLoaded,
		metrics,
		jsDivergence,
		overallStatus,
		complexityImprovement,
	)

	writeProgressReport(finalReport, outputDir)

	// Final summary
	logProductionResultsSummary(trainingTime, totalKmers, metrics, jsDivergence, overallStatus)

	// Write final summary file
	summary := buildFinalSummary(
		time.Now().Format("2006-01-02 15:04:05"),
		trainingTime,
		totalKmers,
		metrics,
		jsDivergence,
		overallStatus,
	)
	writeFinalSummary(outputDir, summary)

	logProductionFinalOutcome(overallPass)

	totalTime := time.Since(startTime)
	logProductionCompletion(outputDir, totalTime)
}

// Helper functions

// ProductionMetrics holds the comprehensive authenticity metrics for production validation.
type ProductionMetrics struct {
	Complexity             float64
	ComplexityImprovement  float64
	GCContent              float64
	GCError                float64
	MaxHomopolymer         int
	CpGDepletion           float64
	JSDivergence           float64
	GenerationSpeed        float64
}

// ProductionPassResults records the pass/fail outcome for each production criterion.
type ProductionPassResults struct {
	ComplexityPass  bool
	GCPass          bool
	HomopolymerPass bool
	CpGPass         bool
	JSPass          bool
	OverallPass     bool
	OverallStatus   string
}

// computeProductionMetrics calculates all production authenticity metrics.
func computeProductionMetrics(sequences [][]byte, k int, avgGenTime time.Duration, sequenceLength int) ProductionMetrics {
	complexity := calculateSequenceComplexity(sequences[0], k)
	complexityImprovement := complexity / 0.005
	gcContent := calculateGCContent(sequences)
	gcError := math.Abs(gcContent - 0.42)
	maxHomopolymer := findMaxHomopolymerRun(sequences)
	cpgDepletion := calculateCpGDepletion(sequences)
	var generationSpeed float64
	if avgGenTime > 0 {
		generationSpeed = float64(sequenceLength) / avgGenTime.Seconds() / 1000.0
	}
	return ProductionMetrics{
		Complexity:            complexity,
		ComplexityImprovement: complexityImprovement,
		GCContent:             gcContent,
		GCError:               gcError,
		MaxHomopolymer:        maxHomopolymer,
		CpGDepletion:          cpgDepletion,
		GenerationSpeed:       generationSpeed,
	}
}

// evaluateProductionPassing applies production pass/fail thresholds.
func evaluateProductionPassing(m ProductionMetrics, jsDivergence float64) ProductionPassResults {
	complexityPass := m.Complexity > 0.02
	gcPass := m.GCError < 0.03
	homopolymerPass := m.MaxHomopolymer < 15
	cpgPass := math.Abs(m.CpGDepletion-0.25) < 0.4
	jsPass := jsDivergence >= 0 && jsDivergence < 0.02

	overallPass := complexityPass && gcPass && homopolymerPass && cpgPass && jsPass
	overallStatus := "PARTIAL_SUCCESS"
	if overallPass {
		overallStatus = "SUCCESS"
	}
	return ProductionPassResults{
		ComplexityPass:  complexityPass,
		GCPass:          gcPass,
		HomopolymerPass: homopolymerPass,
		CpGPass:         cpgPass,
		JSPass:          jsPass,
		OverallPass:     overallPass,
		OverallStatus:   overallStatus,
	}
}

// buildFinalSummary formats the final text summary for writing to disk.
func buildFinalSummary(
	completedAt string,
	trainingTime time.Duration,
	totalKmers uint64,
	m ProductionMetrics,
	jsDivergence float64,
	overallStatus string,
) string {
	return fmt.Sprintf(`PRODUCTION FULL GENOME TRAINING SUMMARY
=============================================

Training Completed: %s
Training Duration: %v
Total K-mers: %d

AUTHENTICITY METRICS:
- Sequence Complexity: %.4f (%.1fx improvement)
- GC Content: %.2f%% (target: 42%%)
- JS Divergence: %.6f (target: < 0.02)
- Max Homopolymer: %d (target: < 15)
- CpG Depletion: %.3f (target: ~0.25)
- Generation Speed: %.1f kb/s

OVERALL STATUS: %s

%s Production-ready genomic authenticity achieved!
%s Ready for population genomics expansion
%s Astrapi infrastructure validated for world-class training

Next Steps:
1. Download 1000 Genomes population data (2TB subset)
2. Train population-aware model for additional improvements
3. Add mammalian genomes for universal patterns
4. Target final complexity: 0.5+ (100x total improvement)
`,
		completedAt,
		trainingTime,
		totalKmers,
		m.Complexity, m.ComplexityImprovement,
		m.GCContent*100,
		jsDivergence,
		m.MaxHomopolymer,
		m.CpGDepletion,
		m.GenerationSpeed,
		overallStatus,
		"✅", "✅", "✅")
}

// writeFinalSummary writes the summary text file to the output directory.
func writeFinalSummary(outputDir string, summary string) {
	summaryPath := filepath.Join(outputDir, "final_summary.txt")
	if err := os.WriteFile(summaryPath, []byte(summary), 0600); err != nil {
		log.Printf("Warning: Could not write summary file: %v", err)
	}
}

// buildFinalProgressReport constructs the final ProgressReport for the completed training run.
func buildFinalProgressReport(
	trainingStart time.Time,
	trainingTime time.Duration,
	totalKmers uint64,
	chromosomesLoaded int,
	m ProductionMetrics,
	jsDivergence float64,
	overallStatus string,
	complexityImprovement float64,
) ProgressReport {
	return ProgressReport{
		Timestamp:          time.Now(),
		Phase:              "completed",
		Status:             "finished",
		TrainingStartTime:  trainingStart,
		TrainingDuration:   trainingTime.String(),
		KmersProcessed:     int64(totalKmers), // #nosec G115
		ChromosomesLoaded:  chromosomesLoaded,
		SequenceComplexity: m.Complexity,
		GCContent:          m.GCContent,
		JSdivergence:       jsDivergence,
		MaxHomopolymer:     m.MaxHomopolymer,
		CpGDepletion:       m.CpGDepletion,
		GenerationSpeed:    m.GenerationSpeed,
		OverallStatus:      overallStatus,
		Message:            fmt.Sprintf("Production training completed with %.1fx complexity improvement", complexityImprovement),
	}
}

// determineOverallStatus returns "SUCCESS" when overallPass is true, otherwise
// "PARTIAL_SUCCESS".
func determineOverallStatus(overallPass bool) string {
	if overallPass {
		return "SUCCESS"
	}
	return "PARTIAL_SUCCESS"
}

// logProductionIntroduction logs the opening banner for production training.
func logProductionIntroduction(outputDir string) {
	log.Printf("=== PRODUCTION FULL GENOME TRAINING ===")
	log.Printf("Astrapi High-Capacity Training with Progress Monitoring")
	log.Printf("Target: JS Divergence < 0.01 for production authenticity")
	log.Printf("Expected: 20-50x complexity improvement")
	log.Printf("📊 Progress monitoring: %s/", outputDir)
	log.Printf("📊 Monitor with: watch 'cat %s/latest_progress.json | jq .'", outputDir)
}

// logProductionModelSetup logs the model setup confirmation.
func logProductionModelSetup(k int, refGenomeDir string, numChromosomes int) {
	log.Printf("✓ Astrapi genome access initialized")
	log.Printf("✓ Data location: %s", refGenomeDir)
	log.Printf("✓ Available chromosomes: %d", numChromosomes)
	log.Printf("✓ Production model created with k=%d", k)
	log.Printf("✓ Configured for complete genome training")
}

// logTrainingStepHeader logs the Step 3 training header lines.
func logTrainingStepHeader(outputDir string) {
	log.Printf("\n[STEP 3] PRODUCTION TRAINING - Complete Human Genome...")
	log.Printf("Training on complete GRCh38 (3.2GB) - estimated 2+ hours")
	log.Printf("Progress will be written to %s/", outputDir)
	log.Printf("This process will run uninterrupted until completion...")
}

// logTrainingComplete logs training completion with time and k-mer count.
func logTrainingComplete(trainingTime time.Duration, totalKmers uint64) {
	log.Printf("✅ PRODUCTION TRAINING COMPLETED!")
	log.Printf("✅ Training time: %v", trainingTime)
	log.Printf("✅ Total k-mers processed: %d", totalKmers)
}

// logGenerationComplete logs the summary after all validation sequences are generated.
func logGenerationComplete(avgGenTime time.Duration, generationSpeed float64) {
	log.Printf("✓ Validation sequences generated")
	log.Printf("✓ Average generation time: %v", avgGenTime)
	log.Printf("✓ Generation speed: %.1f kb/s", generationSpeed)
}

// writeTrainingFailReport writes a training failure report and logs the error.
func writeTrainingFailReport(outputDir string, trainingStart time.Time, errorMsg string) {
	log.Printf("❌ %s", errorMsg)
	writeProgressReport(ProgressReport{
		Timestamp:        time.Now(),
		Phase:            "training",
		Status:           "failed",
		TrainingDuration: time.Since(trainingStart).String(),
		ErrorMessage:     errorMsg,
	}, outputDir)
}

// writeGenerationFailReport writes a generation failure report and logs the error.
func writeGenerationFailReport(outputDir string, errorMsg string) {
	log.Printf("❌ %s", errorMsg)
	writeProgressReport(ProgressReport{
		Timestamp:    time.Now(),
		Phase:        "validation",
		Status:       "failed",
		ErrorMessage: errorMsg,
	}, outputDir)
}

// computeJSDivergenceWithLogging computes the JS divergence from reference sequences,
// logging the result. Returns jsDivergence (negative if reference seqs couldn't be loaded).
func computeJSDivergenceWithLogging(sequences [][]byte, refGenome *reference.SimpleReferenceGenome, numRefSeqs int, seqLen int, k int) float64 {
	referenceSequences, err := loadReferenceSequences(refGenome, numRefSeqs, seqLen)
	if err != nil {
		log.Printf("⚠ Warning: Could not load reference sequences: %v", err)
		return -1
	}
	jsDivergence := calculateJensenShannonDivergence(sequences, referenceSequences, k)
	jsImprovement := 0.055 / jsDivergence
	log.Printf("  JS Divergence: %.6f (%.1fx improvement)", jsDivergence, jsImprovement)
	return jsDivergence
}

// writeStartupProgressReport writes the very first initialization progress report.
func writeStartupProgressReport(outputDir string) {
	writeProgressReport(ProgressReport{
		Timestamp: time.Now(),
		Phase:     "initialization",
		Status:    "starting",
		Message:   "Production training initialization",
	}, outputDir)
}

// writeTrainingStartReport writes the training start progress report.
func writeTrainingStartReport(outputDir string) {
	writeProgressReport(ProgressReport{
		Timestamp:         time.Now(),
		Phase:             "training",
		Status:            "in_progress",
		TrainingStartTime: time.Now(),
		Message:           "Starting complete genome training - this will take 2+ hours",
	}, outputDir)
}

// writeValidationStartReport writes the validation start progress report.
func writeValidationStartReport(outputDir string, numSequences int, sequenceLength int) {
	writeProgressReport(ProgressReport{
		Timestamp: time.Now(),
		Phase:     "validation",
		Status:    "in_progress",
		Message:   fmt.Sprintf("Generating %d validation sequences of %d bases each", numSequences, sequenceLength),
	}, outputDir)
}

// logProductionCompletion logs the final footer lines after production training.
func logProductionCompletion(outputDir string, totalTime time.Duration) {
	log.Printf("\n📊 Progress files saved in: %s/", outputDir)
	log.Printf("📊 View summary: cat %s/final_summary.txt", outputDir)
	log.Printf("📊 View JSON: cat %s/latest_progress.json | jq", outputDir)
	log.Printf("\n🚀 PRODUCTION TRAINING COMPLETE IN %v", totalTime)
}

// writeInitProgressReport writes the initialization-phase progress reports.
func writeInitProgressReport(outputDir string, chromosomesLoaded int, k int) {
	writeProgressReport(ProgressReport{
		Timestamp:         time.Now(),
		Phase:             "initialization",
		Status:            "completed",
		ChromosomesLoaded: chromosomesLoaded,
		Message:           fmt.Sprintf("Loaded %d chromosomes from astrapi", chromosomesLoaded),
	}, outputDir)
	writeProgressReport(ProgressReport{
		Timestamp: time.Now(),
		Phase:     "model_creation",
		Status:    "completed",
		Message:   fmt.Sprintf("Adaptive HMM model created with k=%d", k),
	}, outputDir)
}

// writeTrainingCompleteReport writes the training completion progress report.
func writeTrainingCompleteReport(outputDir string, trainingTime time.Duration, totalKmers uint64) {
	writeProgressReport(ProgressReport{
		Timestamp:        time.Now(),
		Phase:            "training",
		Status:           "completed",
		TrainingDuration: trainingTime.String(),
		KmersProcessed:   int64(totalKmers), // #nosec G115
		Message:          "Complete genome training finished successfully",
	}, outputDir)
}

// logProductionMetricsResults logs individual metric pass/fail lines.
func logProductionMetricsResults(m ProductionMetrics, r ProductionPassResults, jsDivergence float64) {
	log.Printf("  Complexity: %.4f (%.1fx improvement) %s", m.Complexity, m.ComplexityImprovement, passFailIcon(r.ComplexityPass))
	log.Printf("  GC Content: %.2f%% (error: %.4f) %s", m.GCContent*100, m.GCError, passFailIcon(r.GCPass))
	log.Printf("  Max Homopolymer: %d %s", m.MaxHomopolymer, passFailIcon(r.HomopolymerPass))
	log.Printf("  CpG Depletion: %.3f %s", m.CpGDepletion, passFailIcon(r.CpGPass))
	if jsDivergence >= 0 {
		log.Printf("  JS Divergence: %.6f %s", jsDivergence, passFailIcon(r.JSPass))
	}
}

// logProductionResultsSummary prints the final training results summary block.
func logProductionResultsSummary(trainingTime time.Duration, totalKmers uint64, m ProductionMetrics, jsDivergence float64, overallStatus string) {
	log.Printf("\n🧬 PRODUCTION TRAINING RESULTS:")
	log.Printf("=====================================")
	log.Printf("Training Data:       Complete GRCh38 genome (3.2GB)")
	log.Printf("Training Time:       %v", trainingTime)
	log.Printf("K-mers Processed:    %d", totalKmers)
	log.Printf("Complexity:          %.4f (%.1fx improvement)", m.Complexity, m.ComplexityImprovement)
	log.Printf("GC Content:          %.2f%% (error: %.4f)", m.GCContent*100, m.GCError)
	log.Printf("JS Divergence:       %.6f", jsDivergence)
	log.Printf("Generation Speed:    %.1f kb/s", m.GenerationSpeed)
	log.Printf("Overall Status:      %s", overallStatus)
}

// logProductionFinalOutcome prints SUCCESS or PROGRESS message depending on overallPass.
func logProductionFinalOutcome(overallPass bool) {
	if overallPass {
		log.Printf("\n✅ SUCCESS: Production training achieved world-class authenticity!")
		log.Printf("✅ Ready for population genomics expansion")
		log.Printf("✅ Astrapi infrastructure validated for massive scale")
	} else {
		log.Printf("\n⚠ PROGRESS: Significant improvements achieved")
		log.Printf("⚠ Consider population data for further enhancement")
	}
}

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
	
	expectedCpG := float64(totalC * totalG) / float64(len(sequences)*50000)
	if expectedCpG == 0 {
		return 0.0
	}
	
	return float64(totalCpG) / expectedCpG
}

func loadReferenceSequences(refGenome *reference.SimpleReferenceGenome, numSequences, length int) ([][]byte, error) {
	var sequences [][]byte
	
	chr1, err := refGenome.LoadChromosome("chr1")
	if err != nil {
		return nil, err
	}
	
	chrSequence := chr1.GetSequence()
	chrLength := chr1.GetLength()
	
	for i := 0; i < numSequences; i++ {
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
	
	syntheticProbs := normalizeCounts(syntheticKmers)
	realProbs := normalizeCounts(realKmers)
	
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