// Command k6_enhanced_training trains a k=6 authenticity.AdaptiveHMMModel
// on the YRI VCFs under the hard-coded path
// /Volumes/Public/genomics_training/populations, generates test sequences,
// prints complexity and GC content against a 42% ±1.5% target, and writes a
// summary to k6_enhanced_progress/. The "Phase 1" figures it compares
// against (complexity 0.0819, GC 48.17%) are constants from an earlier run,
// not recomputed.
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

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== PHASE 1.2 K=6 + ENHANCED GC CONTROL ===")
	fmt.Println("Testing k=6 with enhanced GC content balancing")
	fmt.Println("Target: Complexity >0.08 + GC Content 42% ±1.5%")

	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n",
		totalCores, perfCores, totalCores-perfCores)
	fmt.Printf("⚡ Using %d P-cores for I/O intensive VCF processing\n", perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)

	// Configuration
	populationDir := "/Volumes/Public/genomics_training/populations"
	progressDir := "k6_enhanced_progress"

	// Create progress directory
	if err := os.MkdirAll(progressDir, 0750); err != nil {
		log.Fatalf("Failed to create progress directory: %v", err)
	}

	// Check dataset
	if _, err := os.Stat(populationDir); os.IsNotExist(err) {
		log.Fatalf("1000 Genomes dataset not found: %s", populationDir)
	}

	// Test on single population first (YRI - highest diversity)
	testPopulations := []string{"YRI"}
	fmt.Printf("📊 PHASE 1.2 TEST: Processing %v population for enhanced GC control\n", testPopulations)
	fmt.Printf("🧬 YRI population chosen for highest genetic diversity\n")
	fmt.Printf("🎯 Enhanced GC target: 42%% ±1.5%% (40.5-43.5%% range)\n")
	fmt.Println()

	// Initialize model with k=6 and enhanced GC control
	fmt.Println("[STEP 1] Initializing k=6 + enhanced GC control model...")
	model := authenticity.NewAdaptiveHMMModel(6)
	fmt.Printf("✅ Model initialized with k=6 + ENHANCED GC CONTROL\n")
	fmt.Printf("📊 State space size: %d states\n", calculateStateSpaceSize(6))
	fmt.Printf("🎯 GC tolerance: ±1.5%% (enhanced from ±2%%)\n")
	fmt.Printf("🔧 Progressive GC correction: 3-tier filtering system\n")
	fmt.Printf("⚖️  Adaptive biological weighting: GC-responsive penalties/bonuses\n")
	fmt.Println()

	// Process real VCF files
	fmt.Println("[STEP 2] PROCESSING VCF FILES WITH ENHANCED GC CONTROL...")
	fmt.Printf("🧬 Reading actual genomic variants from disk\n")
	fmt.Printf("⚡ Apple Silicon optimization: %d P-cores\n", perfCores)
	fmt.Printf("💾 I/O intensive: Reading compressed VCF files\n")
	fmt.Printf("🎯 Target complexity: >0.08 (maintain Phase 1 success)\n")
	fmt.Printf("🎯 Target GC content: 40.5-43.5%% (enhanced control)\n")
	fmt.Println()

	start := time.Now()
	totalBytes, totalVariants := processTestPopulations(model, populationDir, testPopulations, perfCores)
	duration := time.Since(start)

	fmt.Printf("✅ K=6 + ENHANCED GC CONTROL COMPLETED!\n")
	fmt.Printf("✅ Processing time: %v\n", duration)
	fmt.Printf("✅ Data processed: %s\n", formatBytes(totalBytes))
	fmt.Printf("✅ Variants processed: %d\n", totalVariants)
	fmt.Printf("✅ Throughput: %.2f MB/s\n", float64(totalBytes)/duration.Seconds()/1024/1024)
	fmt.Println()

	// Generate validation sequences
	fmt.Println("[STEP 3] Generating enhanced validation sequences...")
	sequences, err := generateValidationSequences(model, 5, 50000)
	if err != nil {
		log.Fatalf("Failed to generate validation sequences: %v", err)
	}
	fmt.Printf("✅ Generated %d validation sequences\n", len(sequences))
	fmt.Println()

	// Calculate authenticity metrics
	fmt.Println("[STEP 4] Calculating enhanced authenticity metrics...")
	complexity := calculateComplexity(sequences)
	gcContent := calculateGCContent(sequences)

	// Calculate improvements
	complexityImprovement := complexity / 0.0051
	gcDeviation := abs(gcContent - 0.42)
	gcImprovementFromPhase1 := abs(0.4817-0.42) / gcDeviation // Improvement over Phase 1

	fmt.Printf("📊 PHASE 1.2 ENHANCED RESULTS:\n")
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
	if gcContent >= 0.405 && gcContent <= 0.435 {
		fmt.Printf(" (deviation: %.3f%%) ✅ EXCELLENT GC CONTROL", gcDeviation*100)
	} else if gcContent >= 0.40 && gcContent <= 0.44 {
		fmt.Printf(" (deviation: %.3f%%) 🟡 GOOD GC CONTROL", gcDeviation*100)
	} else {
		fmt.Printf(" (deviation: %.3f%%) ❌ NEEDS MORE WORK", gcDeviation*100)
	}
	fmt.Println()

	if gcImprovementFromPhase1 > 1.0 {
		fmt.Printf("  GC Improvement: %.1fx better than Phase 1 🎯\n", gcImprovementFromPhase1)
	}
	fmt.Println()

	// Phase 1.2 Assessment
	fmt.Println("🎉 PHASE 1.2 ENHANCED TRAINING RESULTS:")
	fmt.Println("=======================================")
	fmt.Printf("Training Method:     k=6 + Enhanced GC Control + Apple Silicon\n")
	fmt.Printf("Test Population:     YRI (highest diversity)\n")
	fmt.Printf("Data Processed:      %s real genomic data\n", formatBytes(totalBytes))
	fmt.Printf("Processing Time:     %v\n", duration)
	fmt.Printf("Variants Processed:  %d real genomic variants\n", totalVariants)
	fmt.Printf("Apple Silicon:       %d P-cores utilized\n", perfCores)
	fmt.Printf("K-mer Size:          6\n")
	fmt.Printf("State Space:         %d states\n", calculateStateSpaceSize(6))
	fmt.Printf("GC Control:          Enhanced ±1.5%% tolerance + progressive correction\n")
	fmt.Printf("Complexity:          %.4f (%.1fx vs baseline)\n", complexity, complexityImprovement)
	fmt.Printf("GC Content:          %.2f%% (deviation: %.3f%%)\n", gcContent*100, gcDeviation*100)

	// Decision Logic
	nextStep := ""
	status := ""
	if complexity > 0.08 && gcContent >= 0.405 && gcContent <= 0.435 {
		nextStep = "PROCEED TO FULL 11-POPULATION TRAINING"
		status = "PHASE 1.2 COMPLETE SUCCESS"
		fmt.Printf("🎯 STATUS: PHASE 1.2 TARGETS MET - next step is full 11-population training\n")
	} else if complexity > 0.05 && gcContent >= 0.40 && gcContent <= 0.44 {
		nextStep = "PROCEED TO FULL TRAINING WITH MONITORING"
		status = "PHASE 1.2 GOOD PROGRESS"
		fmt.Printf("🔄 STATUS: GOOD PROGRESS - Continue with careful monitoring\n")
	} else {
		nextStep = "CONSIDER PHASE 2 STRUCTURAL IMPROVEMENTS"
		status = "PHASE 1.2 NEEDS REFINEMENT"
		fmt.Printf("⚠  STATUS: NEEDS REFINEMENT - May need additional improvements\n")
	}

	fmt.Printf("Next Action:         %s\n", nextStep)
	fmt.Println()
	fmt.Printf("📊 Results saved in: %s/\n", progressDir)
	fmt.Printf("🏁 PHASE 1.2 ENHANCED TRAINING COMPLETE IN %v\n", duration)

	// Comparison with Phase 1
	fmt.Println()
	fmt.Println("📈 PHASE 1 vs PHASE 1.2 COMPARISON:")
	fmt.Println("(Phase 1 values are constants from an earlier run)")
	fmt.Println("====================================")
	fmt.Printf("Complexity:    0.0819 → %.4f", complexity)
	if complexity >= 0.0819 {
		fmt.Printf(" (maintained ✅)\n")
	} else {
		fmt.Printf(" (decreased %.1f%% ⚠)\n", (0.0819-complexity)/0.0819*100)
	}
	fmt.Printf("GC Content:    48.17%% → %.2f%%", gcContent*100)
	if gcDeviation < abs(0.4817-0.42) {
		fmt.Printf(" (improved %.1fx ✅)\n", gcImprovementFromPhase1)
	} else {
		fmt.Printf(" (similar deviation ⚠)\n")
	}

	// Save detailed results for documentation
	savePhase12Results(progressDir, duration, complexity, gcContent, totalBytes, totalVariants, perfCores, status)
}

func calculateStateSpaceSize(k int) int {
	// State space is 4^(k-1) for (k-1)-mers
	result := 1
	for i := 0; i < k-1; i++ {
		result *= 4
	}
	return result
}

func processTestPopulations(model *authenticity.AdaptiveHMMModel, populationDir string, populations []string, numCores int) (int64, int64) {
	var totalBytes int64
	var totalVariants int64
	var wg sync.WaitGroup

	fmt.Printf("🔥 Starting k=6 + enhanced GC control with %d P-cores\n", numCores)

	// Process populations with limited concurrency to avoid overwhelming I/O
	semaphore := make(chan struct{}, 2) // Limit to 2 concurrent populations

	for _, pop := range populations {
		popDir := filepath.Join(populationDir, pop)
		if _, err := os.Stat(popDir); os.IsNotExist(err) {
			fmt.Printf("⚠ Population %s not found, skipping\n", pop)
			continue
		}

		wg.Add(1)
		go func(population string) {
			defer wg.Done()

			semaphore <- struct{}{}        // Acquire
			defer func() { <-semaphore }() // Release

			bytes, variants := processPopulationVCFs(model, population, populationDir, numCores)
			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)

			fmt.Printf("✅ Population %s complete: %s processed, %d variants\n",
				population, formatBytes(bytes), variants)
		}(pop)
	}

	wg.Wait()

	return totalBytes, totalVariants
}

func processPopulationVCFs(model *authenticity.AdaptiveHMMModel, pop, populationDir string, numCores int) (int64, int64) {
	popDir := filepath.Join(populationDir, pop)
	fmt.Printf("🧬 Processing population %s with enhanced GC control (real VCF files)...\n", pop)

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

			bytes, variants := processRealVCFFile(model, file, pop)
			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)

			if variants > 0 {
				fmt.Printf("  📄 %s: %s, %d variants\n",
					filepath.Base(file), formatBytes(bytes), variants)
			}
		}(vcfFile)
	}

	wg.Wait()

	return totalBytes, totalVariants
}

func processRealVCFFile(model *authenticity.AdaptiveHMMModel, vcfFile, pop string) (int64, int64) {
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

	// Read and process the VCF file
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

		// Add sequences to our training data
		if len(ref) > 0 && len(alt) > 0 {
			sequenceBuilder.WriteString(ref)
			sequenceBuilder.WriteString(alt)
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

func savePhase12Results(progressDir string, duration time.Duration, complexity, gcContent float64, totalBytes, totalVariants int64, perfCores int, status string) {
	resultsFile := filepath.Join(progressDir, "phase1_2_enhanced_results.txt")

	complexityImprovement := complexity / 0.0051
	gcDeviation := abs(gcContent - 0.42)
	gcImprovementFromPhase1 := abs(0.4817-0.42) / gcDeviation

	results := fmt.Sprintf(`PHASE 1.2 ENHANCED GC CONTROL RESULTS
====================================

Training Completed: %s
Phase: 1.2 (Enhanced GC Content Control)
Test Population: YRI (highest diversity)
Total Training Time: %v
Apple Silicon P-Cores: %d

MODEL CONFIGURATION:
- K-mer Size: 6
- State Space: %d states
- GC Tolerance: ±1.5%% (enhanced from ±2%%)
- Progressive GC Correction: 3-tier filtering system
- Adaptive Biological Weighting: GC-responsive penalties/bonuses

AUTHENTICITY METRICS:
- Sequence Complexity: %.4f (%.1fx improvement vs baseline)
- GC Content: %.2f%% (deviation: %.3f%% from 42%% target)
- GC Improvement vs Phase 1: %.1fx better
- Phase 1.2 Status: %s

PERFORMANCE SUMMARY:
- Data Processed: %s
- Variants Processed: %d
- Processing Time: %v
- Throughput: %.2f MB/s
- Apple Silicon: %d P-cores utilized

COMPARISON WITH PHASE 1 (Phase 1 values are constants from an earlier run):
- Complexity: 0.0819 → %.4f %s
- GC Content: 48.17%% → %.2f%% %s
- GC Deviation: 6.17%% → %.3f%% %s

ENHANCED GC CONTROL FEATURES:
✅ Tighter tolerance: ±1.5%% vs ±2%%
✅ Progressive correction: 3-tier filtering (slight/moderate/severe)
✅ Adaptive weighting: GC-responsive penalties and bonuses
✅ Local GC monitoring: State-based GC content tracking
✅ Multi-level response: Different strategies by deviation severity

DECISION MATRIX:
- Complexity >0.08 + GC 40.5-43.5%%: Proceed to full 11-population training
- Complexity >0.05 + GC 40-44%%: Proceed with monitoring
- Below thresholds: Consider Phase 2 improvements

CURRENT RESULT: %.4f complexity, %.2f%% GC → %s

NEXT STEPS:
%s

TARGET CHECK:
%s Complexity >0.08 and GC 40.5-43.5%%

Training Status: PHASE 1.2 %s
`, time.Now().Format("2006-01-02 15:04:05"), duration, perfCores,
		calculateStateSpaceSize(6), complexity, complexityImprovement, gcContent*100, gcDeviation*100,
		gcImprovementFromPhase1, status, formatBytes(totalBytes), totalVariants, duration,
		float64(totalBytes)/duration.Seconds()/1024/1024, perfCores,
		complexity, getComparisonSymbol(complexity, 0.0819), gcContent*100,
		getGCComparisonSymbol(gcDeviation, abs(0.4817-0.42)), gcDeviation*100,
		getImprovementSymbol(gcDeviation, abs(0.4817-0.42)),
		complexity, gcContent*100, status,
		getNextStepsPhase12(complexity, gcContent),
		getSuccessEmojiPhase12(complexity, gcContent), status)

	if err := os.WriteFile(resultsFile, []byte(results), 0600); err != nil {
		fmt.Printf("⚠ Could not save results file: %v\n", err)
	}
}

func getComparisonSymbol(current, previous float64) string {
	if current >= previous {
		return "(maintained ✅)"
	} else {
		decrease := (previous - current) / previous * 100
		return fmt.Sprintf("(decreased %.1f%% ⚠)", decrease)
	}
}

func getGCComparisonSymbol(currentDev, previousDev float64) string {
	if currentDev < previousDev {
		improvement := previousDev / currentDev
		return fmt.Sprintf("(improved %.1fx ✅)", improvement)
	} else {
		return "(similar deviation ⚠)"
	}
}

func getImprovementSymbol(currentDev, previousDev float64) string {
	if currentDev < previousDev {
		return "✅"
	} else {
		return "⚠"
	}
}

func getNextStepsPhase12(complexity, gcContent float64) string {
	if complexity > 0.08 && gcContent >= 0.405 && gcContent <= 0.435 {
		return "1. Proceed to full 11-population training\n2. Prepare for expert validation"
	} else if complexity > 0.05 && gcContent >= 0.40 && gcContent <= 0.44 {
		return "1. Proceed to full training with careful GC monitoring\n2. Additional GC control tuning if needed\n3. Validate on full dataset"
	} else {
		return "1. Consider Phase 2 structural improvements\n2. Alternative VCF processing approaches\n3. Context-aware training methods"
	}
}

func getSuccessEmojiPhase12(complexity, gcContent float64) string {
	if complexity > 0.08 && gcContent >= 0.405 && gcContent <= 0.435 {
		return "✅"
	} else if complexity > 0.05 && gcContent >= 0.40 && gcContent <= 0.44 {
		return "🟡"
	} else {
		return "❌"
	}
}
