// Command population_weighted_training trains a k=6
// authenticity.AdaptiveHMMModel on the YRI, CEU and CHB VCFs under the
// hard-coded path /Volumes/Public/genomics_training/populations, repeating
// each population's training sequences according to a fixed diversity weight
// (YRI 2.0, CEU 1.0, CHB 1.2). It prints complexity and GC content of
// generated sequences and writes a summary to population_weighted_progress/.
// The weights are chosen by hand, and earlier-phase figures (complexity
// 0.0051 and 0.0819) are constants from earlier runs.
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

// PopulationWeights represents diversity-based training weights
type PopulationWeights struct {
	// African populations (highest diversity)
	African map[string]float64
	// European populations (moderate diversity)
	European map[string]float64
	// Asian populations (moderate diversity)
	Asian map[string]float64
	// Admixed populations (high diversity)
	Admixed map[string]float64
}

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== PHASE 2.2 POPULATION-WEIGHTED TRAINING ===")
	fmt.Println("Population-specific complexity weighting based on genetic diversity")
	fmt.Println("Target: Complexity >0.08 + GC Content 42% ±2%")

	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n",
		totalCores, perfCores, totalCores-perfCores)
	fmt.Printf("⚡ Using %d P-cores for weighted processing\n", perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)

	// Configuration
	populationDir := "/Volumes/Public/genomics_training/populations"
	progressDir := "population_weighted_progress"

	// Create progress directory
	if err := os.MkdirAll(progressDir, 0750); err != nil {
		log.Fatalf("Failed to create progress directory: %v", err)
	}

	// Check dataset
	if _, err := os.Stat(populationDir); os.IsNotExist(err) {
		log.Fatalf("1000 Genomes dataset not found: %s", populationDir)
	}

	// Initialize population weights based on genetic diversity
	weights := initializePopulationWeights()

	// Test on multiple populations with diversity weighting
	testPopulations := []string{"YRI", "CEU", "CHB"} // High, medium, medium diversity
	fmt.Printf("📊 PHASE 2.2 TEST: Processing %v populations with diversity weighting\n", testPopulations)
	fmt.Printf("🧬 Population weighting based on known genetic diversity patterns\n")
	fmt.Printf("⚖️  African populations: %.1fx weight (highest diversity)\n", weights.African["YRI"])
	fmt.Printf("⚖️  European populations: %.1fx weight (moderate diversity)\n", weights.European["CEU"])
	fmt.Printf("⚖️  Asian populations: %.1fx weight (moderate diversity)\n", weights.Asian["CHB"])
	fmt.Println()

	// Initialize model with k=6 (proven to work)
	fmt.Println("[STEP 1] Initializing k=6 + population weighting model...")
	model := authenticity.NewAdaptiveHMMModel(6)
	fmt.Printf("✅ Model initialized with k=6 + POPULATION WEIGHTING\n")
	fmt.Printf("📊 State space size: %d states\n", calculateStateSpaceSize(6))
	fmt.Printf("⚖️  Diversity weighting: African 2.0x, European 1.0x, Asian 1.2x, Admixed 1.8x\n")
	fmt.Printf("🎯 Hypothesis: high-diversity populations raise complexity\n")
	fmt.Println()

	// Process population-weighted VCF files
	fmt.Println("[STEP 2] PROCESSING VCF FILES WITH POPULATION WEIGHTING...")
	fmt.Printf("🧬 Applying diversity-based training weights\n")
	fmt.Printf("⚡ Apple Silicon optimization: %d P-cores\n", perfCores)
	fmt.Printf("💾 I/O intensive: Reading compressed VCF files\n")
	fmt.Printf("🎯 Target complexity: >0.08 (leverage high-diversity populations)\n")
	fmt.Printf("🎯 Target GC content: 42%% ±2%% (population-balanced)\n")
	fmt.Println()

	start := time.Now()
	totalBytes, totalVariants := processWeightedPopulations(model, populationDir, testPopulations, weights, perfCores)
	duration := time.Since(start)

	fmt.Printf("✅ POPULATION-WEIGHTED PROCESSING COMPLETED!\n")
	fmt.Printf("✅ Processing time: %v\n", duration)
	fmt.Printf("✅ Data processed: %s\n", formatBytes(totalBytes))
	fmt.Printf("✅ Variants processed: %d\n", totalVariants)
	fmt.Printf("✅ Throughput: %.2f MB/s\n", float64(totalBytes)/duration.Seconds()/1024/1024)
	fmt.Println()

	// Generate validation sequences
	fmt.Println("[STEP 3] Generating population-weighted validation sequences...")
	sequences, err := generateValidationSequences(model, 5, 50000)
	if err != nil {
		log.Fatalf("Failed to generate validation sequences: %v", err)
	}
	fmt.Printf("✅ Generated %d validation sequences\n", len(sequences))
	fmt.Println()

	// Calculate authenticity metrics
	fmt.Println("[STEP 4] Calculating population-weighted authenticity metrics...")
	complexity := calculateComplexity(sequences)
	gcContent := calculateGCContent(sequences)

	// Calculate improvements
	complexityImprovement := complexity / 0.0051
	gcDeviation := abs(gcContent - 0.42)

	fmt.Printf("📊 PHASE 2.2 POPULATION-WEIGHTED RESULTS:\n")
	fmt.Printf("  Complexity: %.4f", complexity)
	if complexity > 0.08 {
		fmt.Printf(" (%.1fx improvement) ✅ EXCELLENT", complexityImprovement)
	} else if complexity > 0.05 {
		fmt.Printf(" (%.1fx improvement) 🟡 GOOD", complexityImprovement)
	} else {
		fmt.Printf(" (%.1fx improvement) ❌ INSUFFICIENT", complexityImprovement)
	}
	fmt.Println()

	fmt.Printf("  GC Content: %.2f%%", gcContent*100)
	if gcContent >= 0.40 && gcContent <= 0.44 {
		fmt.Printf(" (deviation: %.3f%%) ✅ EXCELLENT GC CONTROL", gcDeviation*100)
	} else if gcContent >= 0.38 && gcContent <= 0.46 {
		fmt.Printf(" (deviation: %.3f%%) 🟡 GOOD GC CONTROL", gcDeviation*100)
	} else {
		fmt.Printf(" (deviation: %.3f%%) ❌ NEEDS MORE WORK", gcDeviation*100)
	}
	fmt.Println()
	fmt.Println()

	// Phase 2.2 Assessment
	fmt.Println("🎉 PHASE 2.2 POPULATION-WEIGHTED TRAINING RESULTS:")
	fmt.Println("=================================================")
	fmt.Printf("Training Method:     k=6 + Population Diversity Weighting + Apple Silicon\n")
	fmt.Printf("Test Populations:    YRI (2.0x), CEU (1.0x), CHB (1.2x)\n")
	fmt.Printf("Data Processed:      %s real genomic data\n", formatBytes(totalBytes))
	fmt.Printf("Processing Time:     %v\n", duration)
	fmt.Printf("Variants Processed:  %d real genomic variants\n", totalVariants)
	fmt.Printf("Apple Silicon:       %d P-cores utilized\n", perfCores)
	fmt.Printf("K-mer Size:          6\n")
	fmt.Printf("State Space:         %d states\n", calculateStateSpaceSize(6))
	fmt.Printf("Weighting Strategy:  Diversity-based population scaling\n")
	fmt.Printf("Complexity:          %.4f (%.1fx vs baseline)\n", complexity, complexityImprovement)
	fmt.Printf("GC Content:          %.2f%% (deviation: %.3f%%)\n", gcContent*100, gcDeviation*100)

	// Decision Logic
	nextStep := ""
	status := ""
	if complexity > 0.08 && gcContent >= 0.40 && gcContent <= 0.44 {
		nextStep = "PROCEED TO FULL 11-POPULATION WEIGHTED TRAINING"
		status = "PHASE 2.2 COMPLETE SUCCESS"
		fmt.Printf("🎯 STATUS: PHASE 2.2 TARGETS MET - next step is full 11-population training\n")
	} else if complexity > 0.05 && gcContent >= 0.38 && gcContent <= 0.46 {
		nextStep = "PROCEED TO FULL TRAINING WITH WEIGHT TUNING"
		status = "PHASE 2.2 GOOD PROGRESS"
		fmt.Printf("🔄 STATUS: GOOD PROGRESS - Fine-tune weights and scale up\n")
	} else {
		nextStep = "CONSIDER PHASE 3 ADVANCED APPROACHES (K=8 OR NEURAL)"
		status = "PHASE 2.2 NEEDS REFINEMENT"
		fmt.Printf("⚠  STATUS: NEEDS REFINEMENT - May need advanced methods\n")
	}

	fmt.Printf("Next Action:         %s\n", nextStep)
	fmt.Println()
	fmt.Printf("📊 Results saved in: %s/\n", progressDir)
	fmt.Printf("🏁 PHASE 2.2 POPULATION-WEIGHTED TRAINING COMPLETE IN %v\n", duration)

	// Population weighting insights
	fmt.Println()
	fmt.Println("📈 POPULATION WEIGHTING INSIGHTS:")
	fmt.Println("=================================")
	fmt.Printf("African populations (YRI): Highest genetic diversity → 2.0x training weight\n")
	fmt.Printf("European populations (CEU): Moderate diversity → 1.0x baseline weight\n")
	fmt.Printf("Asian populations (CHB): Moderate diversity → 1.2x weight\n")
	fmt.Printf("Admixed populations: High diversity → 1.8x weight\n")
	fmt.Println()
	fmt.Printf("Hypothesis (not tested here): weighting diverse populations raises complexity\n")

	// Save detailed results for documentation
	savePhase22Results(progressDir, duration, complexity, gcContent, totalBytes, totalVariants, perfCores, status, weights)
}

func initializePopulationWeights() PopulationWeights {
	return PopulationWeights{
		// African populations have highest genetic diversity
		African: map[string]float64{
			"YRI": 2.0, // Yoruba - highest diversity
			"LWK": 1.9, // Luhya - very high diversity
			"GWD": 1.8, // Gambian - high diversity
		},
		// European populations have moderate diversity
		European: map[string]float64{
			"CEU": 1.0, // Utah Europeans - baseline
			"GBR": 1.0, // British - baseline
			"TSI": 1.1, // Tuscan - slightly higher
		},
		// Asian populations have moderate diversity
		Asian: map[string]float64{
			"CHB": 1.2, // Han Chinese Beijing
			"JPT": 1.1, // Japanese Tokyo
			"CHS": 1.2, // Han Chinese South
		},
		// Admixed populations have high diversity
		Admixed: map[string]float64{
			"MXL": 1.8, // Mexican - high admixture
			"PUR": 1.7, // Puerto Rican - high admixture
		},
	}
}

func getPopulationWeight(pop string, weights PopulationWeights) float64 {
	// Check each population category
	if weight, exists := weights.African[pop]; exists {
		return weight
	}
	if weight, exists := weights.European[pop]; exists {
		return weight
	}
	if weight, exists := weights.Asian[pop]; exists {
		return weight
	}
	if weight, exists := weights.Admixed[pop]; exists {
		return weight
	}
	// Default weight for unknown populations
	return 1.0
}

func calculateStateSpaceSize(k int) int {
	// State space is 4^(k-1) for (k-1)-mers
	result := 1
	for i := 0; i < k-1; i++ {
		result *= 4
	}
	return result
}

func processWeightedPopulations(model *authenticity.AdaptiveHMMModel, populationDir string, populations []string, weights PopulationWeights, numCores int) (int64, int64) {
	var totalBytes int64
	var totalVariants int64
	var wg sync.WaitGroup

	fmt.Printf("🔥 Starting population-weighted VCF processing with %d P-cores\n", numCores)

	// Process populations with limited concurrency to avoid overwhelming I/O
	semaphore := make(chan struct{}, 2) // Limit to 2 concurrent populations

	for _, pop := range populations {
		popDir := filepath.Join(populationDir, pop)
		if _, err := os.Stat(popDir); os.IsNotExist(err) {
			fmt.Printf("⚠ Population %s not found, skipping\n", pop)
			continue
		}

		weight := getPopulationWeight(pop, weights)

		wg.Add(1)
		go func(population string, popWeight float64) {
			defer wg.Done()

			semaphore <- struct{}{}        // Acquire
			defer func() { <-semaphore }() // Release

			bytes, variants := processWeightedPopulationVCFs(model, population, populationDir, popWeight, numCores)
			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)

			fmt.Printf("✅ Population %s complete (%.1fx weight): %s processed, %d variants\n",
				population, popWeight, formatBytes(bytes), variants)
		}(pop, weight)
	}

	wg.Wait()

	return totalBytes, totalVariants
}

func processWeightedPopulationVCFs(model *authenticity.AdaptiveHMMModel, pop, populationDir string, weight float64, numCores int) (int64, int64) {
	popDir := filepath.Join(populationDir, pop)
	fmt.Printf("🧬 Processing population %s with %.1fx diversity weight (real VCF files)...\n", pop, weight)

	// Get all VCF files for this population
	vcfFiles, err := filepath.Glob(filepath.Join(popDir, "*.vcf.gz"))
	if err != nil {
		fmt.Printf("⚠ Error finding VCF files for %s: %v\n", pop, err)
		return 0, 0
	}

	var totalBytes int64
	var totalVariants int64
	var wg sync.WaitGroup

	// Process VCF files with limited concurrency for I/O
	semaphore := make(chan struct{}, 4) // Limit concurrent file processing

	for _, vcfFile := range vcfFiles {
		wg.Add(1)
		go func(file string) {
			defer wg.Done()

			semaphore <- struct{}{}        // Acquire
			defer func() { <-semaphore }() // Release

			bytes, variants := processWeightedVCFFile(model, file, pop, weight)
			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)

			if variants > 0 {
				fmt.Printf("  📄 %s: %s, %d variants (%.1fx weighted)\n",
					filepath.Base(file), formatBytes(bytes), variants, weight)
			}
		}(vcfFile)
	}

	wg.Wait()

	return totalBytes, totalVariants
}

func processWeightedVCFFile(model *authenticity.AdaptiveHMMModel, vcfFile, pop string, weight float64) (int64, int64) {
	// Open the gzipped VCF file
	file, err := os.Open(vcfFile) // #nosec G304
	if err != nil {
		fmt.Printf("⚠ Error opening %s: %v\n", vcfFile, err)
		return 0, 0
	}
	defer file.Close()

	// Get file size
	stat, err := file.Stat()
	if err != nil {
		fmt.Printf("⚠ Error getting file size for %s: %v\n", vcfFile, err)
		return 0, 0
	}
	fileSize := stat.Size()

	// Create gzip reader
	gzipReader, err := gzip.NewReader(file)
	if err != nil {
		fmt.Printf("⚠ Error creating gzip reader for %s: %v\n", vcfFile, err)
		return 0, 0
	}
	defer gzipReader.Close()

	// Read and process the VCF file with population weighting
	scanner := bufio.NewScanner(gzipReader)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024) // 1MB buffer

	var variantCount int64
	var sequenceBuilder strings.Builder

	for scanner.Scan() {
		line := scanner.Text()

		// Skip header lines
		if strings.HasPrefix(line, "#") {
			continue
		}

		// Parse VCF line
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}

		// Extract reference and alternate alleles
		ref := fields[3]
		alt := fields[4]

		// Add sequences to our training data with population weighting
		if len(ref) > 0 && len(alt) > 0 {
			// Apply population weight by repeating sequences
			repetitions := int(weight)
			if weight > float64(repetitions) {
				// Handle fractional weights probabilistically
				if float64(int(weight*10)%10) > 5 {
					repetitions++
				}
			}

			for i := 0; i < repetitions; i++ {
				sequenceBuilder.WriteString(ref)
				sequenceBuilder.WriteString(alt)
			}
			variantCount++
		}

		// Process accumulated sequence when we have enough data
		if sequenceBuilder.Len() > 100000 { // 100KB chunks
			sequence := sequenceBuilder.String()
			if err := model.TrainFromSequence([]byte(sequence)); err != nil {
				fmt.Printf("⚠ Error training on sequence from %s: %v\n", vcfFile, err)
			}
			sequenceBuilder.Reset()
		}
	}

	// Process any remaining sequence
	if sequenceBuilder.Len() > 0 {
		sequence := sequenceBuilder.String()
		if err := model.TrainFromSequence([]byte(sequence)); err != nil {
			fmt.Printf("⚠ Error training on final sequence from %s: %v\n", vcfFile, err)
		}
	}

	if err := scanner.Err(); err != nil {
		fmt.Printf("⚠ Error reading %s: %v\n", vcfFile, err)
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

func calculateComplexity(sequences [][]byte) float64 {
	if len(sequences) == 0 {
		return 0
	}

	totalComplexity := 0.0
	for _, seq := range sequences {
		k := 6 // Using k=6 for complexity calculation
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
		case 8: // M1: 4P + 4E
			return 4
		case 10: // M1 Pro: 6P + 4E
			return 6
		case 12: // M1 Max: 8P + 4E
			return 8
		case 16: // M2 Ultra: 12P + 4E
			return 12
		case 20: // M1 Ultra: 16P + 4E
			return 16
		case 24: // M2 Ultra: 16P + 8E
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

func savePhase22Results(progressDir string, duration time.Duration, complexity, gcContent float64, totalBytes, totalVariants int64, perfCores int, status string, weights PopulationWeights) {
	resultsFile := filepath.Join(progressDir, "phase2_2_population_weighted_results.txt")

	complexityImprovement := complexity / 0.0051
	gcDeviation := abs(gcContent - 0.42)

	results := fmt.Sprintf(`PHASE 2.2 POPULATION-WEIGHTED TRAINING RESULTS
==============================================

Training Completed: %s
Phase: 2.2 (Population-Specific Diversity Weighting)
Test Populations: YRI (2.0x), CEU (1.0x), CHB (1.2x)
Total Training Time: %v
Apple Silicon P-Cores: %d

MODEL CONFIGURATION:
- K-mer Size: 6
- State Space: %d states
- Weighting Strategy: Diversity-based population scaling
- African Populations: %.1fx weight (highest diversity)
- European Populations: %.1fx weight (baseline)
- Asian Populations: %.1fx weight (moderate diversity)
- Admixed Populations: %.1fx weight (high diversity)

AUTHENTICITY METRICS:
- Sequence Complexity: %.4f (%.1fx improvement vs baseline)
- GC Content: %.2f%% (deviation: %.3f%% from 42%% target)
- Phase 2.2 Status: %s

PERFORMANCE SUMMARY:
- Data Processed: %s
- Variants Processed: %d
- Processing Time: %v
- Throughput: %.2f MB/s
- Apple Silicon: %d P-cores utilized

POPULATION WEIGHTING RATIONALE:
- YRI (Yoruba): 2.0x weight - Highest genetic diversity globally
- CEU (Utah Europeans): 1.0x weight - Baseline diversity reference
- CHB (Han Chinese): 1.2x weight - Moderate diversity, large population
- Strategy: Leverage high-diversity populations for complexity enhancement

PHASE PROGRESSION SUMMARY (earlier phases are constants from earlier runs):
- k=4 Baseline: Complexity 0.0051, GC 48.37%%
- Phase 1 (k=6): Complexity 0.0819, GC 48.17%%
- Phase 1.2 (Enhanced GC): Complexity 0.0819, GC 47.90%%
- Phase 2.2 (Pop Weighted): Complexity %.4f, GC %.2f%%

WEIGHTING EFFECTIVENESS:
%s

DECISION MATRIX:
- Complexity >0.08 + GC 40-44%%: Proceed to full 11-population weighted training
- Complexity >0.05 + GC 38-46%%: Proceed with weight tuning
- Below thresholds: Consider Phase 3 advanced approaches (k=8 or neural)

CURRENT RESULT: %.4f complexity, %.2f%% GC → %s

NEXT STEPS:
%s

TARGET CHECK:
%s Complexity >0.08 and GC 40-44%%

Training Status: PHASE 2.2 %s
`, time.Now().Format("2006-01-02 15:04:05"), duration, perfCores,
		calculateStateSpaceSize(6), weights.African["YRI"], weights.European["CEU"],
		weights.Asian["CHB"], weights.Admixed["MXL"],
		complexity, complexityImprovement, gcContent*100, gcDeviation*100, status,
		formatBytes(totalBytes), totalVariants, duration,
		float64(totalBytes)/duration.Seconds()/1024/1024, perfCores,
		complexity, gcContent*100,
		getWeightingEffectiveness(complexity, 0.0819),
		complexity, gcContent*100, status,
		getNextStepsPhase22(complexity, gcContent),
		getSuccessEmojiPhase22(complexity, gcContent), status)

	if err := os.WriteFile(resultsFile, []byte(results), 0600); err != nil {
		fmt.Printf("⚠ Could not save results file: %v\n", err)
	}
}

func getWeightingEffectiveness(currentComplexity, previousComplexity float64) string {
	if currentComplexity > previousComplexity {
		improvement := (currentComplexity - previousComplexity) / previousComplexity * 100
		return fmt.Sprintf("✅ Population weighting increased complexity by %.1f%%", improvement)
	} else if currentComplexity >= previousComplexity*0.95 {
		return "🟡 Population weighting maintained complexity levels"
	} else {
		decrease := (previousComplexity - currentComplexity) / previousComplexity * 100
		return fmt.Sprintf("⚠ Population weighting decreased complexity by %.1f%%\n   May need weight tuning or different approach", decrease)
	}
}

func getNextStepsPhase22(complexity, gcContent float64) string {
	if complexity > 0.08 && gcContent >= 0.40 && gcContent <= 0.44 {
		return "1. Scale to full 11-population weighted training\n2. Prepare for expert validation"
	} else if complexity > 0.05 && gcContent >= 0.38 && gcContent <= 0.46 {
		return "1. Fine-tune population weights\n2. Test on full dataset with adjusted weights\n3. Consider combining with contextual processing"
	} else {
		return "1. Implement Phase 3 advanced approaches\n2. Test k=8 for higher complexity\n3. Evaluate neural network methods"
	}
}

func getSuccessEmojiPhase22(complexity, gcContent float64) string {
	if complexity > 0.08 && gcContent >= 0.40 && gcContent <= 0.44 {
		return "✅"
	} else if complexity > 0.05 && gcContent >= 0.38 && gcContent <= 0.46 {
		return "🟡"
	} else {
		return "❌"
	}
}
