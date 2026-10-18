// Command reference_based_training loads chr22 from a hard-coded GRCh38 FASTA
// (/Volumes/Public/genomics_training/human/GRCh38_complete.fasta), and for
// each YRI VCF record under /Volumes/Public/genomics_training/populations
// trains a k=6 adaptive HMM on the 10 kb reference window around the variant
// and on the same window with the ALT allele applied. It then generates five
// 50 kb sequences, prints complexity and GC content against fixed targets,
// and writes a summary to reference_based_progress/.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
// The "REF+ALT" comparison figures (47.05% GC, complexity baseline 0.0051)
// are constants from an earlier run that is not reproduced here.
package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// ReferenceGenome holds chromosome sequences
type ReferenceGenome struct {
	chromosomes map[string][]byte
	mu          sync.RWMutex
}

// Variant represents a genomic variant
type Variant struct {
	Chromosome string
	Position   int
	Ref        string
	Alt        string
}

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== REFERENCE-BASED VCF TRAINING ===")
	fmt.Println("Training on reference windows around VCF variants")
	fmt.Println("Target: Complexity >0.08 + GC Content 42% ±2%")

	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n",
		totalCores, perfCores, totalCores-perfCores)
	fmt.Printf("⚡ Using %d P-cores for reference-based processing\n", perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)

	// Configuration
	referenceFile := "/Volumes/Public/genomics_training/human/GRCh38_complete.fasta"
	populationDir := "/Volumes/Public/genomics_training/populations"
	progressDir := "reference_based_progress"
	windowSize := 10000 // 10KB genomic windows

	// Create progress directory
	if err := os.MkdirAll(progressDir, 0750); err != nil {
		log.Fatalf("Failed to create progress directory: %v", err)
	}

	// Check files exist
	if _, err := os.Stat(referenceFile); os.IsNotExist(err) {
		log.Fatalf("Reference genome not found: %s", referenceFile)
	}
	if _, err := os.Stat(populationDir); os.IsNotExist(err) {
		log.Fatalf("1000 Genomes dataset not found: %s", populationDir)
	}

	// Test on single population first (YRI - highest diversity)
	testPopulations := []string{"YRI"}
	fmt.Printf("📊 REFERENCE-BASED TEST: Processing %v population\n", testPopulations)
	fmt.Printf("🧬 YRI population chosen for highest genetic diversity\n")
	fmt.Printf("🏗️  Window size: %d bp\n", windowSize)
	fmt.Printf("📚 Reference: GRCh38 (chr22 only)\n")
	fmt.Printf("🎯 GC target: 42%%\n")
	fmt.Println()

	// Load reference genome (subset for testing)
	fmt.Println("[STEP 1] Loading reference genome...")
	refGenome, err := loadReferenceGenome(referenceFile, []string{"chr22"}) // Start with chr22 for speed
	if err != nil {
		log.Fatalf("Failed to load reference genome: %v", err)
	}
	fmt.Printf("✅ Reference genome loaded for testing chromosomes\n")
	fmt.Printf("📚 Chromosomes available: %d\n", len(refGenome.chromosomes))
	fmt.Println()

	// Initialize model with k=6
	fmt.Println("[STEP 2] Initializing k=6 + reference-based model...")
	model := authenticity.NewAdaptiveHMMModel(6)
	fmt.Printf("✅ Model initialized with k=6 + REFERENCE-BASED PROCESSING\n")
	fmt.Printf("📊 Possible 6-mers: 4096\n")
	fmt.Printf("🏗️  Processing method: Variants applied to reference windows\n")
	fmt.Printf("🎯 Complexity target: >0.08\n")
	fmt.Printf("🎯 GC content target: 42%% ±2%%\n")
	fmt.Println()

	// Process reference-based VCF files
	fmt.Println("[STEP 3] PROCESSING VCF FILES WITH REFERENCE-BASED APPROACH...")
	fmt.Printf("🧬 Training on reference and ALT-applied windows\n")
	fmt.Printf("⚡ Apple Silicon optimization: %d P-cores\n", perfCores)
	fmt.Printf("💾 I/O + genomic processing: Reference + variants\n")
	fmt.Printf("🏗️  Window approach: %d bp genomic contexts\n", windowSize)
	fmt.Printf("🎯 Goal: reduce the GC deviation of the earlier REF+ALT run\n")
	fmt.Println()

	start := time.Now()
	totalBytes, totalVariants, totalWindows := processReferenceBasedPopulations(
		model, refGenome, populationDir, testPopulations, windowSize, perfCores)
	duration := time.Since(start)

	fmt.Printf("✅ REFERENCE-BASED PROCESSING COMPLETED!\n")
	fmt.Printf("✅ Processing time: %v\n", duration)
	fmt.Printf("✅ Data processed: %s\n", formatBytes(totalBytes))
	fmt.Printf("✅ Variants processed: %d\n", totalVariants)
	fmt.Printf("✅ Genomic windows generated: %d\n", totalWindows)
	fmt.Printf("✅ Throughput: %.2f MB/s\n", float64(totalBytes)/duration.Seconds()/1024/1024)
	fmt.Println()

	// Generate validation sequences
	fmt.Println("[STEP 4] Generating reference-based validation sequences...")
	sequences, err := generateValidationSequences(model, 5, 50000)
	if err != nil {
		log.Fatalf("Failed to generate validation sequences: %v", err)
	}
	fmt.Printf("✅ Generated %d validation sequences\n", len(sequences))
	fmt.Println()

	// Calculate authenticity metrics
	fmt.Println("[STEP 5] Calculating reference-based authenticity metrics...")
	complexity := calculateComplexity(sequences, 6)
	gcContent := calculateGCContent(sequences)

	// Ratios against constants recorded from an earlier REF+ALT run
	// (complexity 0.0051, GC 47.05%), which is not reproduced here.
	complexityImprovement := complexity / 0.0051
	gcDeviation := abs(gcContent - 0.42)
	gcImprovement := abs(0.4705-0.42) - gcDeviation

	fmt.Printf("📊 REFERENCE-BASED RESULTS:\n")
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
		fmt.Printf(" (deviation: %.3f%%) ✅ WITHIN TARGET", gcDeviation*100)
	} else if gcDeviation < 0.0505 { // Better than Phase 2.2
		fmt.Printf(" (deviation: %.3f%%) 🟡 IMPROVED", gcDeviation*100)
	} else {
		fmt.Printf(" (deviation: %.3f%%) ❌ NOT IMPROVED", gcDeviation*100)
	}
	fmt.Println()

	if gcImprovement > 0 {
		fmt.Printf("  GC Improvement: %.3f%% smaller deviation than the earlier REF+ALT run ✅\n", gcImprovement*100)
	} else {
		fmt.Printf("  GC Improvement: %.3f%% vs the earlier REF+ALT run ❌\n", gcImprovement*100)
	}
	fmt.Println()

	// Assessment against the fixed targets
	productionReady := complexity >= 0.08 && gcDeviation <= 0.02
	gcSolved := gcDeviation < 0.0505 // smaller deviation than the earlier REF+ALT run

	fmt.Println("🎉 REFERENCE-BASED TRAINING RESULTS:")
	fmt.Println("=====================================")
	fmt.Printf("Training Method:     k=6 + reference windows around variants\n")
	fmt.Printf("Test Population:     YRI\n")
	fmt.Printf("Data Processed:      %s of compressed VCF\n", formatBytes(totalBytes))
	fmt.Printf("Processing Time:     %v\n", duration)
	fmt.Printf("Variants Processed:  %d VCF records\n", totalVariants)
	fmt.Printf("Genomic Windows:     %d training windows (REF and ALT-applied)\n", totalWindows)
	fmt.Printf("Window Size:         %d bp\n", windowSize)
	fmt.Printf("Apple Silicon:       %d P-cores utilized\n", perfCores)
	fmt.Printf("Reference Genome:    GRCh38 (chr22 only)\n")
	fmt.Printf("Complexity:          %.4f (%.1fx the earlier run's 0.0051)\n", complexity, complexityImprovement)
	fmt.Printf("GC Content:          %.2f%% (deviation: %.3f%%)\n", gcContent*100, gcDeviation*100)

	// Decision Logic
	status := ""
	nextStep := ""
	if productionReady {
		status = "COMPLEXITY AND GC TARGETS MET"
		nextStep = "SCALE TO FULL 11-POPULATION TRAINING"
		fmt.Printf("🎯 STATUS: Both complexity and GC targets met\n")
	} else if gcSolved && complexity > 0.05 {
		status = "GC DEVIATION REDUCED - COMPLEXITY ABOVE 0.05"
		nextStep = "OPTIMIZE AND SCALE TO FULL DATASET"
		fmt.Printf("🔥 STATUS: GC deviation reduced, complexity above 0.05\n")
	} else if gcSolved {
		status = "GC DEVIATION REDUCED - COMPLEXITY BELOW TARGET"
		nextStep = "COMBINE WITH K=8 OR POPULATION WEIGHTING"
		fmt.Printf("🟡 STATUS: GC deviation reduced, complexity below target\n")
	} else {
		status = "TARGETS NOT MET"
		nextStep = "ANALYZE REFERENCE PROCESSING PARAMETERS"
		fmt.Printf("⚠  STATUS: Targets not met\n")
	}

	fmt.Printf("Next Action:         %s\n", nextStep)
	fmt.Println()
	fmt.Printf("📊 Results saved in: %s/\n", progressDir)
	fmt.Printf("🏁 REFERENCE-BASED TRAINING COMPLETE IN %v\n", duration)

	// Save detailed results
	saveReferenceBasedResults(progressDir, duration, complexity, gcContent, totalBytes,
		totalVariants, totalWindows, perfCores, status, gcImprovement, windowSize)
}

func loadReferenceGenome(referenceFile string, chromosomes []string) (*ReferenceGenome, error) {
	refGenome := &ReferenceGenome{
		chromosomes: make(map[string][]byte),
	}

	file, err := os.Open(referenceFile) // #nosec G304
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024) // 10MB buffer for large chromosomes

	var currentChrom string
	var currentSeq strings.Builder
	chromMap := make(map[string]bool)
	for _, chr := range chromosomes {
		chromMap[chr] = true
	}

	for scanner.Scan() {
		line := scanner.Text()

		if strings.HasPrefix(line, ">") {
			// Save previous chromosome if it was one we wanted
			if currentChrom != "" && chromMap[currentChrom] {
				refGenome.chromosomes[currentChrom] = []byte(strings.ToUpper(currentSeq.String()))
				fmt.Printf("  📚 Loaded %s: %d bp\n", currentChrom, currentSeq.Len())
			}

			// Start new chromosome
			parts := strings.Fields(line)
			currentChrom = strings.TrimPrefix(parts[0], ">")
			currentSeq.Reset()

			// Only process chromosomes we need
			if !chromMap[currentChrom] {
				currentChrom = ""
			}
		} else if currentChrom != "" {
			currentSeq.WriteString(line)
		}
	}

	// Save last chromosome
	if currentChrom != "" && chromMap[currentChrom] {
		refGenome.chromosomes[currentChrom] = []byte(strings.ToUpper(currentSeq.String()))
		fmt.Printf("  📚 Loaded %s: %d bp\n", currentChrom, currentSeq.Len())
	}

	return refGenome, scanner.Err()
}

func processReferenceBasedPopulations(model *authenticity.AdaptiveHMMModel, refGenome *ReferenceGenome,
	populationDir string, populations []string, windowSize, numCores int) (int64, int64, int64) {

	var totalBytes int64
	var totalVariants int64
	var totalWindows int64
	var wg sync.WaitGroup

	fmt.Printf("🔥 Starting reference-based processing with %d P-cores\n", numCores)
	fmt.Printf("🏗️  Creating reference windows around variants\n")

	// Process populations with controlled concurrency
	semaphore := make(chan struct{}, 2) // Allow 2 concurrent populations

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

			bytes, variants, windows := processReferenceBasedPopulationVCFs(
				model, refGenome, population, populationDir, windowSize, numCores)
			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)
			atomic.AddInt64(&totalWindows, windows)

			fmt.Printf("✅ Population %s complete: %s processed, %d variants, %d windows\n",
				population, formatBytes(bytes), variants, windows)
		}(pop)
	}

	wg.Wait()

	return totalBytes, totalVariants, totalWindows
}

func processReferenceBasedPopulationVCFs(model *authenticity.AdaptiveHMMModel, refGenome *ReferenceGenome,
	pop, populationDir string, windowSize, numCores int) (int64, int64, int64) {

	popDir := filepath.Join(populationDir, pop)
	fmt.Printf("🧬 Processing population %s with reference-based approach...\n", pop)

	// Only process chr22 for this test
	vcfFiles := []string{filepath.Join(popDir, "chr22.vcf.gz")}

	var totalBytes int64
	var totalVariants int64
	var totalWindows int64
	var wg sync.WaitGroup

	// Process VCF files with limited concurrency
	semaphore := make(chan struct{}, 3)

	for _, vcfFile := range vcfFiles {
		if _, err := os.Stat(vcfFile); os.IsNotExist(err) {
			continue
		}

		wg.Add(1)
		go func(file string) {
			defer wg.Done()

			semaphore <- struct{}{}        // Acquire
			defer func() { <-semaphore }() // Release

			bytes, variants, windows := processReferenceBasedVCFFile(
				model, refGenome, file, pop, windowSize)
			atomic.AddInt64(&totalBytes, bytes)
			atomic.AddInt64(&totalVariants, variants)
			atomic.AddInt64(&totalWindows, windows)

			if variants > 0 {
				fmt.Printf("  📄 %s: %s, %d variants, %d windows (reference-based)\n",
					filepath.Base(file), formatBytes(bytes), variants, windows)
			}
		}(vcfFile)
	}

	wg.Wait()

	return totalBytes, totalVariants, totalWindows
}

func processReferenceBasedVCFFile(model *authenticity.AdaptiveHMMModel, refGenome *ReferenceGenome,
	vcfFile, pop string, windowSize int) (int64, int64, int64) {

	file, err := os.Open(vcfFile) // #nosec G304
	if err != nil {
		return 0, 0, 0
	}
	defer file.Close()

	stat, err := file.Stat()
	if err != nil {
		return 0, 0, 0
	}
	fileSize := stat.Size()

	gzipReader, err := gzip.NewReader(file)
	if err != nil {
		return 0, 0, 0
	}
	defer gzipReader.Close()

	scanner := bufio.NewScanner(gzipReader)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)

	var variantCount int64
	var windowCount int64

	// Get chromosome from filename
	chromosome := strings.TrimSuffix(filepath.Base(vcfFile), ".vcf.gz")

	// Get reference sequence for this chromosome
	refGenome.mu.RLock()
	refSeq, exists := refGenome.chromosomes[chromosome]
	refGenome.mu.RUnlock()

	if !exists {
		fmt.Printf("⚠ Reference sequence not found for %s\n", chromosome)
		return fileSize, 0, 0
	}

	for scanner.Scan() {
		line := scanner.Text()

		if strings.HasPrefix(line, "#") {
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}

		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			continue
		}

		ref := fields[3]
		alt := fields[4]

		if len(ref) > 0 && len(alt) > 0 {
			// Create genomic window around variant
			windowStart := pos - 1 - windowSize/2 // Convert to 0-based
			window := extractGenomicWindow(refSeq, pos-1, windowSize)
			if len(window) > 0 {
				// Calculate relative position within window
				relativePos := windowSize / 2
				if windowStart < 0 {
					relativePos = pos - 1 // Adjust for start-of-chromosome windows
				}

				// Apply variant to create alternate haplotype
				altWindow := applyVariantToWindow(window, relativePos, ref, alt)

				// Train on both reference and alternate windows
				if err := model.TrainFromSequence(window); err == nil {
					windowCount++
				}
				if len(altWindow) > 0 {
					if err := model.TrainFromSequence(altWindow); err == nil {
						windowCount++
					}
				}
			}
			variantCount++

			// Progress reporting every 10K variants
			if variantCount%10000 == 0 {
				fmt.Printf("  📊 Processed %d variants, generated %d windows\n", variantCount, windowCount)
			}
		}
	}

	return fileSize, variantCount, windowCount
}

func extractGenomicWindow(refSeq []byte, pos, windowSize int) []byte {
	if pos < 0 || pos >= len(refSeq) {
		return nil
	}

	start := pos - windowSize/2
	end := pos + windowSize/2

	if start < 0 {
		start = 0
	}
	if end > len(refSeq) {
		end = len(refSeq)
	}

	if end-start < windowSize/4 { // Minimum window size
		return nil
	}

	return refSeq[start:end]
}

func applyVariantToWindow(window []byte, relativePos int, ref, alt string) []byte {
	if relativePos < 0 || relativePos >= len(window) {
		return nil
	}

	// Handle SNPs (most common)
	if len(ref) == 1 && len(alt) == 1 {
		altWindow := make([]byte, len(window))
		copy(altWindow, window)
		altWindow[relativePos] = alt[0]
		return altWindow
	}

	// Handle simple indels (basic implementation)
	if relativePos+len(ref) <= len(window) {
		// Create new window with variant applied
		altWindow := make([]byte, 0, len(window)+len(alt)-len(ref))

		// Copy sequence before variant
		altWindow = append(altWindow, window[:relativePos]...)

		// Insert alternate allele
		altWindow = append(altWindow, []byte(alt)...)

		// Copy sequence after variant
		altWindow = append(altWindow, window[relativePos+len(ref):]...)

		return altWindow
	}

	return nil
}

// ... (rest of helper functions: generateValidationSequences, calculateComplexity, etc.)
// These are similar to existing implementations

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

func saveReferenceBasedResults(progressDir string, duration time.Duration, complexity, gcContent float64,
	totalBytes, totalVariants, totalWindows int64, perfCores int, status string, gcImprovement float64, windowSize int) {

	resultsFile := filepath.Join(progressDir, "reference_based_results.txt")

	complexityImprovement := complexity / 0.0051
	gcDeviation := abs(gcContent - 0.42)

	results := fmt.Sprintf(`REFERENCE-BASED VCF PROCESSING RESULTS
====================================

Training Completed: %s
Approach: k=6 HMM trained on reference windows, with and without the ALT allele
Test Population: YRI
Total Training Time: %v
Apple Silicon P-Cores: %d

REFERENCE-BASED CONFIGURATION:
- Reference Genome: GRCh38 (chr22 only)
- Window Size: %d bp
- Processing Method: Variants applied to reference sequences
- K-mer Size: 6
- Possible 6-mers: 4096

SEQUENCE STATISTICS:
- Sequence Complexity: %.4f (%.1fx the earlier run's 0.0051)
- GC Content: %.2f%% (deviation: %.3f%% from 42%% target)
- GC Improvement: %.3f%% vs the earlier REF+ALT run
- Status: %s

PROCESSING SUMMARY:
- Data Processed: %s
- Variants Processed: %d
- Genomic Windows: %d training windows
- Processing Time: %v
- Throughput: %.2f MB/s
- Apple Silicon: %d P-cores utilized

COMPARISON WITH THE EARLIER REF+ALT RUN (constant, not reproduced here):
- REF+ALT Concatenation: 47.05%% GC (5.05%% deviation)
- Reference-Based: %.2f%% GC (%.3f%% deviation)
- GC Improvement: %.3f%%

Targets met (complexity >= 0.08 and GC 40-44%%): %s

Next Steps: %s

Processing Status: %s
Completion Time: %s
`, time.Now().Format("2006-01-02 15:04:05"), duration, perfCores,
		windowSize, complexity, complexityImprovement, gcContent*100, gcDeviation*100, gcImprovement*100, status,
		formatBytes(totalBytes), totalVariants, totalWindows, duration,
		float64(totalBytes)/duration.Seconds()/1024/1024, perfCores,
		gcContent*100, gcDeviation*100, gcImprovement*100,
		getSuccessEmoji(complexity, gcContent), getNextSteps(complexity, gcContent),
		status, time.Now().Format("2006-01-02 15:04:05"))

	if err := os.WriteFile(resultsFile, []byte(results), 0600); err != nil {
		fmt.Printf("⚠ Could not save results file: %v\n", err)
	}
}

func getSuccessEmoji(complexity, gcContent float64) string {
	if complexity >= 0.08 && gcContent >= 0.40 && gcContent <= 0.44 {
		return "✅"
	} else if abs(gcContent-0.42) < 0.0505 {
		return "🟡"
	}
	return "🔄"
}

func getNextSteps(complexity, gcContent float64) string {
	gcDeviation := abs(gcContent - 0.42)

	if complexity >= 0.08 && gcDeviation <= 0.02 {
		return "1. Scale to full 11-population training\n2. Compare against held-out real data"
	} else if gcDeviation < 0.0505 {
		return "1. Optimize reference processing parameters\n2. Scale to multiple chromosomes"
	} else {
		return "1. Debug reference genome integration\n2. Validate variant application logic"
	}
}
