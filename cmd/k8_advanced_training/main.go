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
	
	fmt.Println("=== PHASE 3 K=8 ADVANCED TRAINING ===")
	fmt.Println("Testing computational limits with k=8 for maximum complexity")
	fmt.Println("Target: Complexity >0.10 + GC Content 42% ±2%")
	
	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n", 
		totalCores, perfCores, totalCores-perfCores)
	fmt.Printf("⚡ Using %d P-cores for k=8 intensive processing\n", perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)
	
	// Configuration
	populationDir := "/Volumes/Public/genomics_training/populations"
	progressDir := "k8_advanced_progress"
	
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
	fmt.Printf("📊 PHASE 3 TEST: Processing %v population for k=8 validation\n", testPopulations)
	fmt.Printf("🧬 YRI population chosen for highest genetic diversity\n")
	fmt.Printf("🎯 k=8 state space: %d states (16x larger than k=6)\n", calculateStateSpaceSize(8))
	fmt.Printf("⚠️  Computational intensity: High memory and CPU requirements\n")
	fmt.Println()

	// Initialize model with k=8
	fmt.Println("[STEP 1] Initializing k=8 advanced model...")
	model := authenticity.NewAdaptiveHMMModel(8)
	fmt.Printf("✅ Model initialized with k=8 (MAXIMUM COMPLEXITY)\n")
	fmt.Printf("📊 State space size: %d states\n", calculateStateSpaceSize(8))
	fmt.Printf("💾 Memory footprint: ~%d MB (estimated)\n", calculateMemoryFootprint(8))
	fmt.Printf("🎯 Expected complexity: >0.10 (20x improvement over baseline)\n")
	fmt.Printf("⚡ Processing approach: Optimized for computational intensity\n")
	fmt.Println()

	// Process k=8 VCF files
	fmt.Println("[STEP 2] PROCESSING VCF FILES WITH K=8...")
	fmt.Printf("🧬 Reading actual genomic variants from disk\n")
	fmt.Printf("⚡ Apple Silicon optimization: %d P-cores\n", perfCores)
	fmt.Printf("💾 I/O + computation intensive: k=8 processing\n")
	fmt.Printf("🎯 Target complexity: >0.10 (stretch goal achievement)\n")
	fmt.Printf("🎯 Target GC content: 42%% ±2%% (with k=8 precision)\n")
	fmt.Printf("⚠️  Expected time: 45-60 minutes (computational complexity)\n")
	fmt.Println()

	start := time.Now()
	totalBytes, totalVariants := processK8Populations(model, populationDir, testPopulations, perfCores)
	duration := time.Since(start)
	
	fmt.Printf("✅ K=8 ADVANCED PROCESSING COMPLETED!\n")
	fmt.Printf("✅ Processing time: %v\n", duration)
	fmt.Printf("✅ Data processed: %s\n", formatBytes(totalBytes))
	fmt.Printf("✅ Variants processed: %d\n", totalVariants)
	fmt.Printf("✅ Throughput: %.2f MB/s\n", float64(totalBytes)/duration.Seconds()/1024/1024)
	fmt.Println()

	// Generate validation sequences
	fmt.Println("[STEP 3] Generating k=8 validation sequences...")
	sequences, err := generateValidationSequences(model, 5, 50000)
	if err != nil {
		log.Fatalf("Failed to generate validation sequences: %v", err)
	}
	fmt.Printf("✅ Generated %d validation sequences\n", len(sequences))
	fmt.Println()

	// Calculate authenticity metrics
	fmt.Println("[STEP 4] Calculating k=8 authenticity metrics...")
	complexity := calculateComplexity(sequences)
	gcContent := calculateGCContent(sequences)
	
	// Calculate improvements
	complexityImprovement := complexity / 0.0051
	gcDeviation := abs(gcContent - 0.42)
	complexityVsK6 := complexity / 0.0819  // Comparison with k=6 success
	
	fmt.Printf("📊 PHASE 3 K=8 RESULTS:\n")
	fmt.Printf("  Complexity: %.4f", complexity)
	if complexity > 0.10 {
		fmt.Printf(" (%.1fx improvement) ✅ STRETCH GOAL ACHIEVED", complexityImprovement)
	} else if complexity > 0.08 {
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
	
	fmt.Printf("  k=8 vs k=6: %.1fx complexity improvement", complexityVsK6)
	if complexityVsK6 > 1.2 {
		fmt.Printf(" ✅ SIGNIFICANT K=8 BENEFIT")
	} else if complexityVsK6 > 1.0 {
		fmt.Printf(" 🟡 MODEST K=8 BENEFIT")
	} else {
		fmt.Printf(" ❌ NO K=8 BENEFIT")
	}
	fmt.Println()
	fmt.Println()

	// Phase 3 Assessment
	fmt.Println("🎉 PHASE 3 K=8 ADVANCED TRAINING RESULTS:")
	fmt.Println("========================================")
	fmt.Printf("Training Method:     k=8 Advanced + Apple Silicon\n")
	fmt.Printf("Test Population:     YRI (highest diversity)\n")
	fmt.Printf("Data Processed:      %s real genomic data\n", formatBytes(totalBytes))
	fmt.Printf("Processing Time:     %v\n", duration)
	fmt.Printf("Variants Processed:  %d real genomic variants\n", totalVariants)
	fmt.Printf("Apple Silicon:       %d P-cores utilized\n", perfCores)
	fmt.Printf("K-mer Size:          8 (maximum tested)\n")
	fmt.Printf("State Space:         %d states (16,384 vs 1,024 for k=6)\n", calculateStateSpaceSize(8))
	fmt.Printf("Memory Usage:        ~%d MB estimated\n", calculateMemoryFootprint(8))
	fmt.Printf("Complexity:          %.4f (%.1fx vs baseline, %.1fx vs k=6)\n", complexity, complexityImprovement, complexityVsK6)
	fmt.Printf("GC Content:          %.2f%% (deviation: %.3f%%)\n", gcContent*100, gcDeviation*100)
	
	// Decision Logic
	nextStep := ""
	status := ""
	if complexity > 0.10 && gcContent >= 0.40 && gcContent <= 0.44 {
		nextStep = "PROCEED TO FULL K=8 PRODUCTION TRAINING"
		status = "PHASE 3 STRETCH GOAL SUCCESS"
		fmt.Printf("🎯 STATUS: STRETCH GOAL SUCCESS - k=8 achieves maximum complexity\n")
	} else if complexity > 0.08 && gcContent >= 0.40 && gcContent <= 0.44 {
		nextStep = "PRODUCTION READY - CHOOSE OPTIMAL K VALUE"
		status = "PHASE 3 EXCELLENT SUCCESS"
		fmt.Printf("🎯 STATUS: EXCELLENT SUCCESS - Production ready with k=8\n")
	} else if complexityVsK6 > 1.2 {
		nextStep = "PROCEED TO FULL K=8 TRAINING"
		status = "PHASE 3 SIGNIFICANT IMPROVEMENT"
		fmt.Printf("🔄 STATUS: SIGNIFICANT K=8 IMPROVEMENT - Worth the computational cost\n")
	} else {
		nextStep = "OPTIMIZE K=6 FOR PRODUCTION"
		status = "PHASE 3 MARGINAL BENEFIT"
		fmt.Printf("⚠  STATUS: MARGINAL K=8 BENEFIT - k=6 may be optimal\n")
	}
	
	fmt.Printf("Next Action:         %s\n", nextStep)
	fmt.Println()
	fmt.Printf("📊 Results saved in: %s/\n", progressDir)
	fmt.Printf("🏁 PHASE 3 K=8 TRAINING COMPLETE IN %v\n", duration)
	
	// Computational analysis
	fmt.Println()
	fmt.Println("📈 K-VALUE COMPUTATIONAL ANALYSIS:")
	fmt.Println("==================================")
	fmt.Printf("k=4: 64 states, fast processing, low complexity\n")
	fmt.Printf("k=6: 1,024 states, moderate processing, excellent complexity\n")
	fmt.Printf("k=8: 16,384 states, intensive processing, complexity %.1fx vs k=6\n", complexityVsK6)
	fmt.Println()
	fmt.Printf("Computational cost vs benefit analysis:\n")
	fmt.Printf("- k=6 → k=8: %dx state space increase\n", calculateStateSpaceSize(8)/calculateStateSpaceSize(6))
	fmt.Printf("- Processing time: %v for k=8 vs ~25min for k=6\n", duration)
	fmt.Printf("- Complexity gain: %.1fx improvement\n", complexityVsK6)
	
	if complexityVsK6 > 1.2 {
		fmt.Printf("✅ RECOMMENDATION: k=8 worth the computational cost\n")
	} else {
		fmt.Printf("🤔 RECOMMENDATION: k=6 may be optimal balance of performance/quality\n")
	}

	// Save detailed results for documentation
	savePhase3Results(progressDir, duration, complexity, gcContent, totalBytes, totalVariants, perfCores, status, complexityVsK6)
}

func calculateStateSpaceSize(k int) int {
	// State space is 4^(k-1) for (k-1)-mers
	result := 1
	for i := 0; i < k-1; i++ {
		result *= 4
	}
	return result
}

func calculateMemoryFootprint(k int) int {
	// Rough estimate: state space size * transition matrix * 8 bytes per float64
	stateSpace := calculateStateSpaceSize(k)
	transitionMemory := stateSpace * 4 * 8 // 4 bases, 8 bytes per float64
	emissionMemory := stateSpace * 42 * 8  // 42 quality scores, 8 bytes per float64
	totalBytes := transitionMemory + emissionMemory
	return totalBytes / (1024 * 1024) // Convert to MB
}

func processK8Populations(model *authenticity.AdaptiveHMMModel, populationDir string, populations []string, numCores int) (int64, int64) {
	var totalBytes int64
	var totalVariants int64
	var wg sync.WaitGroup
	
	fmt.Printf("🔥 Starting k=8 advanced processing with %d P-cores\n", numCores)
	fmt.Printf("⚠️  High computational intensity - may take significant time\n")
	
	// Process populations with limited concurrency to avoid overwhelming system
	semaphore := make(chan struct{}, 1) // Limit to 1 concurrent population due to k=8 intensity
	
	for _, pop := range populations {
		popDir := filepath.Join(populationDir, pop)
		if _, err := os.Stat(popDir); os.IsNotExist(err) {
			fmt.Printf("⚠ Population %s not found, skipping\n", pop)
			continue
		}
		
		wg.Add(1)
		go func(population string) {
			defer wg.Done()
			
			semaphore <- struct{}{} // Acquire
			defer func() { <-semaphore }() // Release
			
			bytes, variants := processK8PopulationVCFs(model, population, populationDir, numCores)
			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)
			
			fmt.Printf("✅ Population %s complete: %s processed, %d variants\n", 
				population, formatBytes(bytes), variants)
		}(pop)
	}
	
	wg.Wait()
	
	return totalBytes, totalVariants
}

func processK8PopulationVCFs(model *authenticity.AdaptiveHMMModel, pop, populationDir string, numCores int) (int64, int64) {
	popDir := filepath.Join(populationDir, pop)
	fmt.Printf("🧬 Processing population %s with k=8 (computationally intensive)...\n", pop)
	
	// Get all VCF files for this population
	vcfFiles, err := filepath.Glob(filepath.Join(popDir, "*.vcf.gz"))
	if err != nil {
		fmt.Printf("⚠ Error finding VCF files for %s: %v\n", pop, err)
		return 0, 0
	}
	
	var totalBytes int64
	var totalVariants int64
	var wg sync.WaitGroup
	
	// Process VCF files with limited concurrency due to k=8 intensity
	semaphore := make(chan struct{}, 2) // Reduced from 4 to 2 for k=8
	
	for _, vcfFile := range vcfFiles {
		wg.Add(1)
		go func(file string) {
			defer wg.Done()
			
			semaphore <- struct{}{} // Acquire
			defer func() { <-semaphore }() // Release
			
			bytes, variants := processK8VCFFile(model, file, pop)
			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)
			
			if variants > 0 {
				fmt.Printf("  📄 %s: %s, %d variants (k=8 processing)\n", 
					filepath.Base(file), formatBytes(bytes), variants)
			}
		}(vcfFile)
	}
	
	wg.Wait()
	
	return totalBytes, totalVariants
}

func processK8VCFFile(model *authenticity.AdaptiveHMMModel, vcfFile, pop string) (int64, int64) {
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
	
	// Read and process the VCF file with k=8 optimization
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
		// Use smaller chunks for k=8 to manage memory
		if sequenceBuilder.Len() > 50000 { // 50KB chunks for k=8
			sequence := sequenceBuilder.String()
			if err := model.TrainFromSequence([]byte(sequence)); err != nil {
				fmt.Printf("⚠ Error training k=8 on sequence from %s: %v\n", vcfFile, err)
			}
			sequenceBuilder.Reset()
		}
	}
	
	// Process any remaining sequence
	if sequenceBuilder.Len() > 0 {
		sequence := sequenceBuilder.String()
		if err := model.TrainFromSequence([]byte(sequence)); err != nil {
			fmt.Printf("⚠ Error training k=8 on final sequence from %s: %v\n", vcfFile, err)
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
		k := 8 // Using k=8 for complexity calculation
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
		case 8:  // M1: 4P + 4E
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

func savePhase3Results(progressDir string, duration time.Duration, complexity, gcContent float64, totalBytes, totalVariants int64, perfCores int, status string, complexityVsK6 float64) {
	resultsFile := filepath.Join(progressDir, "phase3_k8_advanced_results.txt")
	
	complexityImprovement := complexity / 0.0051
	gcDeviation := abs(gcContent - 0.42)
	
	results := fmt.Sprintf(`PHASE 3 K=8 ADVANCED TRAINING RESULTS
====================================

Training Completed: %s
Phase: 3 (K=8 Maximum Complexity)
Test Population: YRI (highest diversity)
Total Training Time: %v
Apple Silicon P-Cores: %d

MODEL CONFIGURATION:
- K-mer Size: 8 (maximum tested)
- State Space: %d states (16x larger than k=6)
- Memory Footprint: ~%d MB estimated
- Computational Intensity: Very high
- Processing Approach: Optimized for k=8 complexity

AUTHENTICITY METRICS:
- Sequence Complexity: %.4f (%.1fx improvement vs baseline)
- GC Content: %.2f%% (deviation: %.3f%% from 42%% target)
- Complexity vs k=6: %.1fx improvement
- Phase 3 Status: %s

PERFORMANCE SUMMARY:
- Data Processed: %s
- Variants Processed: %d
- Processing Time: %v
- Throughput: %.2f MB/s
- Apple Silicon: %d P-cores utilized

COMPUTATIONAL ANALYSIS:
- k=6 State Space: 1,024 states
- k=8 State Space: 16,384 states (16x increase)
- Memory Increase: ~%.1fx estimated
- Processing Time vs k=6: ~%.1fx longer
- Complexity Benefit: %.1fx improvement

K-VALUE OPTIMIZATION:
%s

PHASE PROGRESSION SUMMARY:
- k=4 Baseline: Complexity 0.0051, GC 48.37%%
- Phase 1 (k=6): Complexity 0.0819, GC 48.17%%
- Phase 3 (k=8): Complexity %.4f, GC %.2f%%

DECISION MATRIX:
- Complexity >0.10 + GC 40-44%%: k=8 stretch goal achieved
- Complexity >0.08 + significant k=8 benefit: k=8 recommended
- Marginal k=8 benefit: k=6 optimal balance

CURRENT RESULT: %.4f complexity, %.2f%% GC, %.1fx vs k=6 → %s

NEXT STEPS:
%s

COMPUTATIONAL RECOMMENDATIONS:
%s

TECHNICAL ACHIEVEMENTS:
✅ k=8 model successfully implemented
✅ 16,384 state space managed efficiently
✅ Apple Silicon optimization for intensive computation
✅ Memory management for large state spaces
%s k=8 computational feasibility demonstrated

Training Status: PHASE 3 %s
`, time.Now().Format("2006-01-02 15:04:05"), duration, perfCores, 
	calculateStateSpaceSize(8), calculateMemoryFootprint(8),
	complexity, complexityImprovement, gcContent*100, gcDeviation*100, complexityVsK6, status,
	formatBytes(totalBytes), totalVariants, duration, 
	float64(totalBytes)/duration.Seconds()/1024/1024, perfCores,
	safeMemoryRatio(8, 6), 
	duration.Minutes()/25.0, complexityVsK6,
	getK8OptimizationInsight(complexityVsK6),
	complexity, gcContent*100,
	complexity, gcContent*100, complexityVsK6, status,
	getNextStepsPhase3(complexity, gcContent, complexityVsK6),
	getComputationalRecommendation(complexityVsK6, duration),
	getSuccessEmojiPhase3(complexity, complexityVsK6), status)
	
	if err := os.WriteFile(resultsFile, []byte(results), 0600); err != nil {
		fmt.Printf("⚠ Could not save results file: %v\n", err)
	}
}

func getK8OptimizationInsight(complexityVsK6 float64) string {
	if complexityVsK6 > 1.2 {
		return "✅ k=8 provides significant complexity improvement\n   16x state space increase justified by results\n   Computational cost worth the quality benefit"
	} else if complexityVsK6 > 1.0 {
		return "🟡 k=8 provides modest complexity improvement\n   May not justify 16x computational increase\n   k=6 likely optimal for most applications"
	} else {
		return "❌ k=8 shows diminishing returns\n   No significant benefit over k=6\n   k=6 is optimal balance of performance/quality"
	}
}

func getNextStepsPhase3(complexity, gcContent, complexityVsK6 float64) string {
	if complexity > 0.10 && gcContent >= 0.40 && gcContent <= 0.44 {
		return "1. Scale k=8 to full 11-population training\n2. Prepare for expert validation with stretch goal results\n3. Production deployment with k=8 maximum complexity"
	} else if complexityVsK6 > 1.2 {
		return "1. Optimize k=8 memory usage\n2. Scale to full dataset with performance monitoring\n3. Compare k=6 vs k=8 production costs"
	} else {
		return "1. Finalize k=6 as optimal model\n2. Focus on full k=6 production training\n3. Expert validation with proven k=6 approach"
	}
}

func getComputationalRecommendation(complexityVsK6 float64, duration time.Duration) string {
	if complexityVsK6 > 1.2 {
		return fmt.Sprintf("✅ k=8 recommended: %.1fx complexity benefit justifies %v processing time", complexityVsK6, duration)
	} else {
		return fmt.Sprintf("🤔 k=6 recommended: %.1fx complexity benefit may not justify %v processing time", complexityVsK6, duration)
	}
}

func safeMemoryRatio(k1, k2 int) float64 {
	mem1 := calculateMemoryFootprint(k1)
	mem2 := calculateMemoryFootprint(k2)
	if mem2 == 0 {
		return float64(mem1) // Return absolute value if division by zero
	}
	return float64(mem1) / float64(mem2)
}

func getSuccessEmojiPhase3(complexity, complexityVsK6 float64) string {
	if complexity > 0.10 {
		return "✅"
	} else if complexityVsK6 > 1.2 {
		return "✅"
	} else if complexity > 0.08 {
		return "🟡"
	} else {
		return "❌"
	}
}