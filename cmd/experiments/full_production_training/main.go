// Command full_production_training trains one k=6 adaptive HMM on the
// gzipped VCFs of 11 populations under the hard-coded directory
// /Volumes/Public/genomics_training/populations, generates sequences, checks
// complexity and GC content against fixed thresholds, and writes
// full_production_progress/production_final_results.txt. The training input
// is each record's REF and ALT strings concatenated, repeated by a
// hard-coded per-population weight (rounded to 1 or 2), so it is not
// genomic sequence. Training errors are ignored. "Production" in the
// identifiers means only that both thresholds were met.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"sync/atomic"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// ProductionConfig represents the final optimized configuration
type ProductionConfig struct {
	KmerSize           int
	PopulationWeights  map[string]float64
	TargetComplexity   float64
	TargetGCContent    float64
	ProcessingApproach string
	OptimizationLevel  string
}

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== FULL POPULATION-WEIGHTED TRAINING ===")
	fmt.Println("11-population 1000 Genomes VCF training (concatenated REF+ALT alleles, population-weighted)")
	fmt.Println("Pass thresholds: complexity and GC content (see below)")

	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 CPU: %d cores (%d assumed performance cores + %d others)\n",
		totalCores, perfCores, totalCores-perfCores)
	fmt.Printf("⚡ GOMAXPROCS set to %d\n", perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)

	// Configuration - Will be optimized based on phase results
	config := getOptimalProductionConfig()

	// Configuration
	populationDir := "/Volumes/Public/genomics_training/populations"
	progressDir := "full_production_progress"

	// Create progress directory
	if err := os.MkdirAll(progressDir, 0750); err != nil {
		log.Fatalf("Failed to create progress directory: %v", err)
	}

	// Check dataset
	if _, err := os.Stat(populationDir); os.IsNotExist(err) {
		log.Fatalf("1000 Genomes dataset not found: %s", populationDir)
	}

	// Full 11-population training set
	allPopulations := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"}

	fmt.Printf("📊 Processing %d populations\n", len(allPopulations))
	fmt.Printf("⚖️  Population weighting: %s\n", config.ProcessingApproach)
	fmt.Printf("🎯 Target complexity: %.3f\n", config.TargetComplexity)
	fmt.Printf("🎯 Target GC content: %.1f%% ±2%%\n", config.TargetGCContent*100)
	fmt.Println()

	// Initialize model with optimal configuration
	fmt.Println("[STEP 1] Initializing model...")
	model := authenticity.NewAdaptiveHMMModel(config.KmerSize)
	fmt.Printf("✅ Model initialized\n")
	fmt.Printf("📊 K-mer size: %d\n", config.KmerSize)
	fmt.Printf("📊 State space size: %d states\n", calculateStateSpaceSize(config.KmerSize))
	fmt.Printf("⚖️  Weighting strategy: %s\n", config.ProcessingApproach)
	fmt.Printf("🎯 Optimization level: %s\n", config.OptimizationLevel)
	fmt.Println()

	// VCF processing
	fmt.Println("[STEP 2] VCF PROCESSING...")
	fmt.Printf("🧬 Training input: REF+ALT allele strings, repeated by population weight\n")
	fmt.Printf("⚡ Workers: %d\n", perfCores)
	fmt.Printf("💾 I/O intensive: Reading compressed VCF files\n")
	fmt.Printf("🎯 Target complexity: %.3f\n", config.TargetComplexity)
	fmt.Printf("🎯 Target GC content: %.1f%% ±2%%\n", config.TargetGCContent*100)
	fmt.Println()

	start := time.Now()
	totalBytes, totalVariants := processProductionPopulations(model, populationDir, allPopulations, config, perfCores)
	duration := time.Since(start)

	fmt.Printf("✅ VCF PROCESSING COMPLETED\n")
	fmt.Printf("✅ Processing time: %v\n", duration)
	fmt.Printf("✅ Data processed: %s\n", formatBytes(totalBytes))
	fmt.Printf("✅ Variants processed: %d\n", totalVariants)
	fmt.Printf("✅ Throughput: %.2f MB/s\n", float64(totalBytes)/duration.Seconds()/1024/1024)
	fmt.Println()

	// Generate final validation sequences
	fmt.Println("[STEP 3] Generating validation sequences...")
	sequences, err := generateValidationSequences(model, 10, 100000) // Larger validation set
	if err != nil {
		log.Fatalf("Failed to generate validation sequences: %v", err)
	}
	fmt.Printf("✅ Generated %d validation sequences\n", len(sequences))
	fmt.Println()

	// Calculate final authenticity metrics
	fmt.Println("[STEP 4] Calculating sequence metrics...")
	complexity := calculateComplexity(sequences, config.KmerSize)
	gcContent := calculateGCContent(sequences)

	// Advanced metrics for expert validation
	gcDeviation := abs(gcContent - config.TargetGCContent)
	productionReady := complexity >= config.TargetComplexity && gcDeviation <= 0.02

	fmt.Printf("📊 RESULTS:\n")
	fmt.Printf("  Complexity: %.4f", complexity)
	if complexity >= config.TargetComplexity {
		fmt.Printf(" ✅ ABOVE THRESHOLD")
	} else {
		fmt.Printf(" ❌ BELOW THRESHOLD")
	}
	fmt.Println()

	fmt.Printf("  GC Content: %.2f%%", gcContent*100)
	if gcDeviation <= 0.02 {
		fmt.Printf(" (deviation: %.3f%%) ✅ WITHIN 2%%", gcDeviation*100)
	} else if gcDeviation <= 0.04 {
		fmt.Printf(" (deviation: %.3f%%) 🟡 ACCEPTABLE", gcDeviation*100)
	} else {
		fmt.Printf(" (deviation: %.3f%%) ❌ NEEDS REFINEMENT", gcDeviation*100)
	}
	fmt.Println()

	fmt.Printf("  Overall Status: ")
	if productionReady {
		fmt.Printf("✅ BOTH THRESHOLDS MET\n")
	} else {
		fmt.Printf("🔄 THRESHOLDS NOT MET\n")
	}
	fmt.Println()

	// Threshold assessment
	fmt.Println("🎉 TRAINING RESULTS:")
	fmt.Println("===================================")
	fmt.Printf("Training Method:     %s\n", config.ProcessingApproach)
	fmt.Printf("Dataset:             %s\n", populationDir)
	fmt.Printf("Populations:         %d\n", len(allPopulations))
	fmt.Printf("Data Processed:      %s of compressed VCF\n", formatBytes(totalBytes))
	fmt.Printf("Processing Time:     %v\n", duration)
	fmt.Printf("Variants Processed:  %d\n", totalVariants)
	fmt.Printf("Workers:             %d\n", perfCores)
	fmt.Printf("Model Configuration: k=%d, %d states\n", config.KmerSize, calculateStateSpaceSize(config.KmerSize))
	fmt.Printf("Optimization:        %s\n", config.OptimizationLevel)
	fmt.Printf("Complexity:          %.4f (target: %.3f)\n", complexity, config.TargetComplexity)
	fmt.Printf("GC Content:          %.2f%% (target: %.1f%% ±2%%)\n", gcContent*100, config.TargetGCContent*100)

	// Expert validation readiness
	fmt.Println()
	fmt.Println("🔬 THRESHOLD CHECK:")
	fmt.Println("===============================")
	if productionReady {
		fmt.Printf("✅ Complexity threshold achieved: %.4f >= %.3f\n", complexity, config.TargetComplexity)
		fmt.Printf("✅ GC content within tolerance: %.3f%% deviation <= 2%%\n", gcDeviation*100)

		fmt.Println()
		fmt.Printf("🎯 NEXT STEPS:\n")
		fmt.Printf("1. Expert blind testing with real bioinformaticians\n")
		fmt.Printf("2. Integration with existing sequencer simulation\n")
		fmt.Printf("3. Performance benchmarking and optimization\n")
	} else {
		fmt.Printf("🔄 Additional optimization required:\n")
		if complexity < config.TargetComplexity {
			fmt.Printf("   - Complexity below threshold (%.4f < %.3f)\n", complexity, config.TargetComplexity)
		}
		if gcDeviation > 0.02 {
			fmt.Printf("   - GC content deviation too high (%.3f%% > 2%%)\n", gcDeviation*100)
		}
		fmt.Printf("   - Consider advanced techniques or parameter tuning\n")
	}

	fmt.Println()
	fmt.Printf("📊 Results saved in: %s/\n", progressDir)
	fmt.Printf("🏁 TRAINING COMPLETE IN %v\n", duration)

	// Save results
	saveProductionResults(progressDir, duration, complexity, gcContent, totalBytes, totalVariants, perfCores, config, productionReady)
}

// getOptimalProductionConfig returns the best configuration based on phase testing results
func getOptimalProductionConfig() ProductionConfig {
	// This will be updated based on the results from our phase testing
	// k=6 is the setting that passed the complexity threshold in earlier runs

	return ProductionConfig{
		KmerSize: 6,
		PopulationWeights: map[string]float64{
			// African populations (highest diversity)
			"YRI": 2.0, "LWK": 1.9, "GWD": 1.8,
			// European populations (moderate diversity)
			"CEU": 1.0, "GBR": 1.0, "TSI": 1.1,
			// Asian populations (moderate diversity)
			"CHB": 1.2, "JPT": 1.1, "CHS": 1.2,
			// Admixed populations (high diversity)
			"MXL": 1.8, "PUR": 1.7,
		},
		TargetComplexity:   0.08, // pass threshold
		TargetGCContent:    0.42, // Human genomic average
		ProcessingApproach: "Population-Weighted k=6",
		OptimizationLevel:  "Experimental",
	}
}

func calculateStateSpaceSize(k int) int {
	result := 1
	for i := 0; i < k-1; i++ {
		result *= 4
	}
	return result
}

func processProductionPopulations(model *authenticity.AdaptiveHMMModel, populationDir string, populations []string, config ProductionConfig, numCores int) (int64, int64) {
	var totalBytes int64
	var totalVariants int64
	var wg sync.WaitGroup

	fmt.Printf("🔥 Starting VCF processing with %d workers\n", numCores)
	fmt.Printf("🌍 Processing %d populations\n", len(populations))

	// Process up to three populations at a time
	semaphore := make(chan struct{}, 3) // 3 concurrent populations

	startTime := time.Now()

	for i, pop := range populations {
		popDir := filepath.Join(populationDir, pop)
		if _, err := os.Stat(popDir); os.IsNotExist(err) {
			fmt.Printf("⚠ Population %s not found, skipping\n", pop)
			continue
		}

		weight := config.PopulationWeights[pop]
		if weight == 0 {
			weight = 1.0 // Default weight
		}

		wg.Add(1)
		go func(population string, popWeight float64, index int) {
			defer wg.Done()

			semaphore <- struct{}{}        // Acquire
			defer func() { <-semaphore }() // Release

			popStart := time.Now()
			bytes, variants := processWeightedPopulationVCFs(model, population, populationDir, popWeight, numCores)
			popDuration := time.Since(popStart)

			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)

			// Progress reporting
			elapsed := time.Since(startTime)
			remaining := time.Duration(float64(elapsed) * float64(len(populations)-index-1) / float64(index+1))

			fmt.Printf("✅ Population %s complete (%d/%d, %.1fx weight): %s processed, %d variants in %v\n",
				population, index+1, len(populations), popWeight, formatBytes(bytes), variants, popDuration)
			fmt.Printf("📊 Progress: %.1f%% complete, ~%v remaining\n",
				float64(index+1)*100.0/float64(len(populations)), remaining)
		}(pop, weight, i)
	}

	wg.Wait()

	return totalBytes, totalVariants
}

func processWeightedPopulationVCFs(model *authenticity.AdaptiveHMMModel, pop, populationDir string, weight float64, numCores int) (int64, int64) {
	popDir := filepath.Join(populationDir, pop)
	fmt.Printf("🧬 Processing population %s with %.1fx weight...\n", pop, weight)

	vcfFiles, err := filepath.Glob(filepath.Join(popDir, "*.vcf.gz"))
	if err != nil {
		fmt.Printf("⚠ Error finding VCF files for %s: %v\n", pop, err)
		return 0, 0
	}

	var totalBytes int64
	var totalVariants int64
	var wg sync.WaitGroup

	// Up to six VCF files at a time
	semaphore := make(chan struct{}, 6) // 6 concurrent files

	for _, vcfFile := range vcfFiles {
		wg.Add(1)
		go func(file string) {
			defer wg.Done()

			semaphore <- struct{}{}        // Acquire
			defer func() { <-semaphore }() // Release

			bytes, variants := processWeightedVCFFile(model, file, pop, weight)
			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)
		}(vcfFile)
	}

	wg.Wait()

	return totalBytes, totalVariants
}

func processWeightedVCFFile(model *authenticity.AdaptiveHMMModel, vcfFile, pop string, weight float64) (int64, int64) {
	file, err := os.Open(vcfFile) // #nosec G304
	if err != nil {
		return 0, 0
	}
	defer file.Close()

	stat, err := file.Stat()
	if err != nil {
		return 0, 0
	}
	fileSize := stat.Size()

	gzipReader, err := gzip.NewReader(file)
	if err != nil {
		return 0, 0
	}
	defer gzipReader.Close()

	scanner := bufio.NewScanner(gzipReader)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)

	var variantCount int64
	var sequenceBuilder strings.Builder

	for scanner.Scan() {
		line := scanner.Text()

		if strings.HasPrefix(line, "#") {
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}

		ref := fields[3]
		alt := fields[4]

		if len(ref) > 0 && len(alt) > 0 {
			// Apply population weight
			repetitions := int(weight)
			if weight > float64(repetitions) && (int(weight*10)%10) > 5 {
				repetitions++
			}

			for i := 0; i < repetitions; i++ {
				sequenceBuilder.WriteString(ref)
				sequenceBuilder.WriteString(alt)
			}
			variantCount++
		}

		if sequenceBuilder.Len() > 100000 {
			sequence := sequenceBuilder.String()
			if err := model.TrainFromSequence([]byte(sequence)); err != nil {
				// errors are ignored
			}
			sequenceBuilder.Reset()
		}
	}

	if sequenceBuilder.Len() > 0 {
		sequence := sequenceBuilder.String()
		model.TrainFromSequence([]byte(sequence))
	}

	return fileSize, variantCount
}

func generateValidationSequences(model *authenticity.AdaptiveHMMModel, numSequences, length int) ([][]byte, error) {
	sequences := make([][]byte, numSequences)

	for i := 0; i < numSequences; i++ {
		result, err := model.GenerateSequence(length)
		if err != nil {
			return nil, fmt.Errorf("failed to generate sequence %d: %w", i+1, err)
		}
		sequences[i] = result.Sequence
	}

	return sequences, nil
}

func calculateComplexity(sequences [][]byte, k int) float64 {
	if len(sequences) == 0 {
		return 0
	}

	totalComplexity := 0.0
	for _, seq := range sequences {
		if len(seq) < k {
			continue
		}

		uniqueKmers := make(map[string]bool)
		for i := 0; i <= len(seq)-k; i++ {
			kmer := string(seq[i : i+k])
			uniqueKmers[kmer] = true
		}

		complexity := float64(len(uniqueKmers)) / float64(len(seq)-k+1)
		totalComplexity += complexity
	}

	return totalComplexity / float64(len(sequences))
}

func calculateGCContent(sequences [][]byte) float64 {
	if len(sequences) == 0 {
		return 0
	}

	totalGC := 0.0
	for _, seq := range sequences {
		if len(seq) == 0 {
			continue
		}

		gcCount := 0
		for _, base := range seq {
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}

		gc := float64(gcCount) / float64(len(seq))
		totalGC += gc
	}

	return totalGC / float64(len(sequences))
}

func formatBytes(bytes int64) string {
	const unit = 1024
	if bytes < unit {
		return fmt.Sprintf("%d B", bytes)
	}
	div, exp := int64(unit), 0
	for n := bytes / unit; n >= unit; n /= unit {
		div *= unit
		exp++
	}
	return fmt.Sprintf("%.1f %cB", float64(bytes)/float64(div), "KMGTPE"[exp])
}

func detectAppleSiliconCores() int {
	totalCores := runtime.NumCPU()

	if totalCores >= 8 {
		switch totalCores {
		case 8:
			return 4
		case 10:
			return 6
		case 12:
			return 8
		case 16:
			return 12
		case 20:
			return 16
		case 24:
			return 16
		default:
			return (totalCores * 2) / 3
		}
	}

	if totalCores > 4 {
		return totalCores - 2
	}

	return totalCores
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func saveProductionResults(progressDir string, duration time.Duration, complexity, gcContent float64, totalBytes, totalVariants int64, perfCores int, config ProductionConfig, productionReady bool) {
	resultsFile := filepath.Join(progressDir, "production_final_results.txt")

	complexityImprovement := complexity / 0.0051
	gcDeviation := abs(gcContent - config.TargetGCContent)

	status := "THRESHOLDS NOT MET"
	if productionReady {
		status = "THRESHOLDS MET"
	}

	results := fmt.Sprintf(`FULL POPULATION-WEIGHTED TRAINING RESULTS (experiment)
======================================================

Training Completed: %s
Training Data: 1000 Genomes VCF REF+ALT allele strings, 11 populations
Total Training Time: %v
Workers (GOMAXPROCS): %d

CONFIGURATION:
- K-mer Size: %d
- State Space: %d states
- Processing Approach: %s
- Optimization Level: %s
- Population Weighting: each record's alleles are repeated round(weight) times
  (weights are hard-coded; the rounding turns them into 1 or 2)

METRICS:
- Sequence Complexity: %.4f (%.1fx the hard-coded 0.0051 baseline)
- Target Complexity: %.3f %s
- GC Content: %.2f%% (deviation: %.3f%% from %.1f%% target)
- GC Tolerance: ±2%% %s
- Threshold Status: %s

DATASET PROCESSING:
- Data Processed: %s of compressed VCF
- Variants Processed: %d
- Processing Time: %v
- Throughput: %.2f MB/s
- Workers: %d

POPULATION WEIGHTS (configured):
- African: YRI (2.0x), LWK (1.9x), GWD (1.8x)
- European: CEU (1.0x), GBR (1.0x), TSI (1.1x)
- Asian: CHB (1.2x), JPT (1.1x), CHS (1.2x)
- Admixed: MXL (1.8x), PUR (1.7x)

THRESHOLD CHECK:
%s

LIMITATIONS:
%s

NOTES:
- Training input is concatenated allele strings, not reference sequence
- Per-population weights are hard-coded, not derived from the data
%s Threshold check (complexity and GC content)

NEXT STEPS:
%s

Final Status: %s
Processing Complete: %s
`, time.Now().Format("2006-01-02 15:04:05"), duration, perfCores,
		config.KmerSize, calculateStateSpaceSize(config.KmerSize), config.ProcessingApproach, config.OptimizationLevel,
		complexity, complexityImprovement, config.TargetComplexity, getThresholdStatus(complexity, config.TargetComplexity),
		gcContent*100, gcDeviation*100, config.TargetGCContent*100, getGCStatus(gcDeviation), status,
		formatBytes(totalBytes), totalVariants, duration, float64(totalBytes)/duration.Seconds()/1024/1024, perfCores,
		getValidationReadiness(productionReady), getDeploymentReadiness(productionReady),
		getProductionEmoji(productionReady), getNextStepsProduction(productionReady), status, time.Now().Format("2006-01-02 15:04:05"))

	if err := os.WriteFile(resultsFile, []byte(results), 0600); err != nil {
		fmt.Printf("⚠ Could not save results file: %v\n", err)
	}
}

func getThresholdStatus(complexity, target float64) string {
	if complexity >= target {
		return "✅"
	}
	return "❌"
}

func getGCStatus(deviation float64) string {
	if deviation <= 0.02 {
		return "✅"
	}
	return "❌"
}

func getValidationReadiness(ready bool) string {
	if ready {
		return "✅ Complexity and GC content within the configured thresholds"
	}
	return "🔄 Complexity or GC content outside the configured thresholds"
}

func getDeploymentReadiness(ready bool) string {
	if ready {
		return "- Thresholds met, but the model is trained on allele strings, not sequence\n- No comparison with real reads or reference sequence was made"
	}
	return "- Thresholds not met\n- The model is trained on allele strings, not sequence"
}

func getProductionEmoji(ready bool) string {
	if ready {
		return "✅"
	}
	return "🔄"
}

func getNextStepsProduction(ready bool) string {
	if ready {
		return "1. Expert blind testing\n2. Integration with existing sequencer simulation\n3. Performance benchmarking"
	}
	return "1. Analyze suboptimal metrics and identify improvements\n2. Consider advanced techniques (k=8, neural networks)\n3. Fine-tune configuration parameters\n4. Rerun training with the new settings"
}
