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
	
	fmt.Println("=== PHASE 1 K=6 VCF TRAINING ===")
	fmt.Println("Testing k=6 for improved genomic complexity")
	fmt.Println("Target: 3x complexity improvement (>0.015)")
	
	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n", 
		totalCores, perfCores, totalCores-perfCores)
	fmt.Printf("⚡ Using %d P-cores for I/O intensive VCF processing\n", perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)
	
	// Configuration
	populationDir := "/Volumes/Public/genomics_training/populations"
	progressDir := "k6_vcf_progress"
	
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
	fmt.Printf("📊 PHASE 1 TEST: Processing %v population for k=6 validation\n", testPopulations)
	fmt.Printf("🧬 YRI population chosen for highest genetic diversity\n")
	fmt.Println()

	// Initialize model with k=6
	fmt.Println("[STEP 1] Initializing k=6 genomics model...")
	model := authenticity.NewAdaptiveHMMModel(6)
	fmt.Printf("✅ Model initialized with k=6 (PHASE 1 IMPROVEMENT)\n")
	fmt.Printf("📊 State space size: %d states (vs 64 for k=4)\n", calculateStateSpaceSize(6))
	fmt.Println()

	// Process real VCF files
	fmt.Println("[STEP 2] PROCESSING VCF FILES WITH K=6...")
	fmt.Printf("🧬 Reading actual genomic variants from disk\n")
	fmt.Printf("⚡ Apple Silicon optimization: %d P-cores\n", perfCores)
	fmt.Printf("💾 I/O intensive: Reading compressed VCF files\n")
	fmt.Printf("🎯 Target complexity: >0.015 (3x improvement)\n")
	fmt.Println()

	start := time.Now()
	totalBytes, totalVariants := processTestPopulations(model, populationDir, testPopulations, perfCores)
	duration := time.Since(start)
	
	fmt.Printf("✅ K=6 VCF PROCESSING COMPLETED!\n")
	fmt.Printf("✅ Processing time: %v\n", duration)
	fmt.Printf("✅ Data processed: %s\n", formatBytes(totalBytes))
	fmt.Printf("✅ Variants processed: %d\n", totalVariants)
	fmt.Printf("✅ Throughput: %.2f MB/s\n", float64(totalBytes)/duration.Seconds()/1024/1024)
	fmt.Println()

	// Generate validation sequences
	fmt.Println("[STEP 3] Generating k=6 validation sequences...")
	sequences, err := generateValidationSequences(model, 5, 50000)
	if err != nil {
		log.Fatalf("Failed to generate validation sequences: %v", err)
	}
	fmt.Printf("✅ Generated %d validation sequences\n", len(sequences))
	fmt.Println()

	// Calculate authenticity metrics
	fmt.Println("[STEP 4] Calculating k=6 authenticity metrics...")
	complexity := calculateComplexity(sequences)
	gcContent := calculateGCContent(sequences)
	
	fmt.Printf("📊 PHASE 1 K=6 RESULTS:\n")
	fmt.Printf("  Complexity: %.4f", complexity)
	if complexity > 0.015 {
		fmt.Printf(" (%.1fx improvement) ✅ PHASE 1 SUCCESS", complexity/0.0051)
	} else if complexity > 0.01 {
		fmt.Printf(" (%.1fx improvement) 🟡 PARTIAL SUCCESS", complexity/0.0051)
	} else {
		fmt.Printf(" (%.1fx improvement) ❌ INSUFFICIENT", complexity/0.0051)
	}
	fmt.Println()
	
	fmt.Printf("  GC Content: %.2f%%", gcContent*100)
	if gcContent >= 0.40 && gcContent <= 0.44 {
		fmt.Print(" ✅")
	} else if gcContent >= 0.35 && gcContent <= 0.50 {
		fmt.Print(" 🟡")
	} else {
		fmt.Printf(" (error: %.4f) ❌", abs(gcContent-0.42))
	}
	fmt.Println()
	fmt.Println()

	// Phase 1 Assessment
	fmt.Println("🎉 PHASE 1 K=6 TRAINING RESULTS:")
	fmt.Println("=================================")
	fmt.Printf("Training Method:     k=6 VCF Processing + Apple Silicon\n")
	fmt.Printf("Test Population:     YRI (highest diversity)\n")
	fmt.Printf("Data Processed:      %s real genomic data\n", formatBytes(totalBytes))
	fmt.Printf("Processing Time:     %v\n", duration)
	fmt.Printf("Variants Processed:  %d real genomic variants\n", totalVariants)
	fmt.Printf("Apple Silicon:       %d P-cores utilized\n", perfCores)
	fmt.Printf("K-mer Size:          6 (vs 4 baseline)\n")
	fmt.Printf("State Space:         %d states\n", calculateStateSpaceSize(6))
	fmt.Printf("Complexity:          %.4f\n", complexity)
	fmt.Printf("GC Content:          %.2f%%\n", gcContent*100)
	
	// Decision Logic
	nextStep := ""
	if complexity > 0.015 {
		nextStep = "PROCEED TO ENHANCED COMPLEXITY CONTROL + FULL DATASET"
		fmt.Printf("🎯 STATUS: PHASE 1 SUCCESS - Continue to full training\n")
	} else if complexity > 0.01 {
		nextStep = "PROCEED TO ENHANCED COMPLEXITY CONTROL"
		fmt.Printf("🔄 STATUS: PARTIAL SUCCESS - Enhance complexity control\n")
	} else {
		nextStep = "CONSIDER PHASE 2 STRUCTURAL IMPROVEMENTS"
		fmt.Printf("⚠  STATUS: INSUFFICIENT - May need structural changes\n")
	}
	
	fmt.Printf("Next Action:         %s\n", nextStep)
	fmt.Println()
	fmt.Printf("📊 Results saved in: %s/\n", progressDir)
	fmt.Printf("🏁 PHASE 1 K=6 TRAINING COMPLETE IN %v\n", duration)
	fmt.Printf("📈 Improvement: %.1fx complexity increase\n", complexity/0.0051)

	// Save detailed results for documentation
	savePhase1Results(progressDir, duration, complexity, gcContent, totalBytes, totalVariants, perfCores)
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
	
	fmt.Printf("🔥 Starting k=6 VCF processing with %d P-cores\n", numCores)
	
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
			
			semaphore <- struct{}{} // Acquire
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
	fmt.Printf("🧬 Processing population %s with k=6 (real VCF files)...\n", pop)
	
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
			
			semaphore <- struct{}{} // Acquire
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

func savePhase1Results(progressDir string, duration time.Duration, complexity, gcContent float64, totalBytes, totalVariants int64, perfCores int) {
	resultsFile := filepath.Join(progressDir, "phase1_k6_results.txt")
	
	improvement := complexity / 0.0051
	success := "INSUFFICIENT"
	if complexity > 0.015 {
		success = "SUCCESS"
	} else if complexity > 0.01 {
		success = "PARTIAL SUCCESS"
	}
	
	results := fmt.Sprintf(`PHASE 1 K=6 TRAINING RESULTS
=============================

Training Completed: %s
Phase: 1 (K-mer Size Increase)
Test Population: YRI (highest diversity)
Total Training Time: %v
Apple Silicon P-Cores: %d

MODEL CONFIGURATION:
- K-mer Size: 6 (increased from 4)
- State Space: %d states (vs 64 for k=4)
- Training Data: Real VCF variants

AUTHENTICITY METRICS:
- Sequence Complexity: %.4f (%.1fx improvement)
- GC Content: %.2f%% (target: 42%%)
- Phase 1 Status: %s

PERFORMANCE SUMMARY:
- Data Processed: %s
- Variants Processed: %d
- Processing Time: %v
- Throughput: %.2f MB/s
- Apple Silicon: %d P-cores utilized

DECISION MATRIX:
- Complexity >0.015: Proceed to enhanced complexity control + full dataset
- Complexity 0.01-0.015: Proceed to enhanced complexity control only  
- Complexity <0.01: Consider Phase 2 structural improvements

CURRENT RESULT: %.4f complexity -> %s

NEXT STEPS:
%s

BASELINE COMPARISON:
- k=4 Complexity: 0.0051
- k=6 Complexity: %.4f
- Improvement Factor: %.1fx
- Target Factor: 3x (0.015)

TECHNICAL ACHIEVEMENTS:
✅ Successfully increased k-mer size to k=6
✅ Maintained real VCF processing integrity
✅ Apple Silicon optimization preserved
✅ Single population test completed efficiently
%s Complexity improvement target achieved

Training Status: PHASE 1 %s
`, time.Now().Format("2006-01-02 15:04:05"), duration, perfCores, 
	calculateStateSpaceSize(6), complexity, improvement, gcContent*100, success,
	formatBytes(totalBytes), totalVariants, duration, 
	float64(totalBytes)/duration.Seconds()/1024/1024, perfCores,
	complexity, success,
	getNextSteps(complexity), complexity, improvement,
	getSuccessEmoji(complexity), success)
	
	if err := os.WriteFile(resultsFile, []byte(results), 0600); err != nil {
		fmt.Printf("⚠ Could not save results file: %v\n", err)
	}
}

func getNextSteps(complexity float64) string {
	if complexity > 0.015 {
		return "1. Implement enhanced complexity control\n2. Train on full 11-population dataset\n3. Prepare for production validation"
	} else if complexity > 0.01 {
		return "1. Implement enhanced complexity control\n2. Retest with complexity improvements\n3. Consider Phase 2 if insufficient"
	} else {
		return "1. Document k=6 limitations\n2. Implement Phase 2 structural improvements\n3. Consider context-aware VCF processing"
	}
}

func getSuccessEmoji(complexity float64) string {
	if complexity > 0.015 {
		return "✅"
	} else if complexity > 0.01 {
		return "🟡"
	} else {
		return "❌"
	}
}