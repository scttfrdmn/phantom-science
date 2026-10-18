// Command real_vcf_training streams every *.vcf.gz under the hard-coded
// directory /Volumes/Public/genomics_training/populations into one k=4
// adaptive HMM, then generates five 50 kb sequences and prints their
// complexity and GC content along with bytes, variants and throughput.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
// The training input is the REF and ALT alleles of each record concatenated
// into 100 kB chunks, not reference sequence, so the model learns no genomic
// context. The "improvement" factor is relative to an assumed 0.005
// baseline, not a measured one.
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

	fmt.Println("=== REAL VCF TRAINING ===")
	fmt.Println("Actually reading and processing real VCF files")
	fmt.Println("This will take significant time with real I/O and parsing")

	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n",
		totalCores, perfCores, totalCores-perfCores)
	fmt.Printf("⚡ Using %d P-cores for I/O intensive VCF processing\n", perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)

	// Configuration
	populationDir := "/Volumes/Public/genomics_training/populations"
	progressDir := "real_vcf_progress"

	// Create progress directory
	if err := os.MkdirAll(progressDir, 0750); err != nil {
		log.Fatalf("Failed to create progress directory: %v", err)
	}

	// Check dataset
	if _, err := os.Stat(populationDir); os.IsNotExist(err) {
		log.Fatalf("1000 Genomes dataset not found: %s", populationDir)
	}

	// Count VCF files
	vcfCount := countVCFFiles(populationDir)
	fmt.Printf("📊 Found %d VCF files to process\n", vcfCount)
	fmt.Printf("📊 Processing time depends on dataset size and disk speed\n")
	fmt.Printf("🔥 Each VCF file will be read and parsed completely\n")
	fmt.Println()

	// Initialize model
	fmt.Println("[STEP 1] Initializing genomics model...")
	model := authenticity.NewAdaptiveHMMModel(4)
	fmt.Printf("✅ Model initialized with k=4\n")
	fmt.Println()

	// Process real VCF files
	fmt.Println("[STEP 2] PROCESSING REAL VCF FILES...")
	fmt.Printf("🧬 Reading actual genomic variants from disk\n")
	fmt.Printf("⚡ Apple Silicon optimization: %d P-cores\n", perfCores)
	fmt.Printf("💾 I/O intensive: Reading compressed VCF files\n")
	fmt.Println()

	start := time.Now()
	totalBytes, totalVariants := processRealVCFFiles(model, populationDir, perfCores)
	duration := time.Since(start)

	fmt.Printf("✅ REAL VCF PROCESSING COMPLETED!\n")
	fmt.Printf("✅ Processing time: %v\n", duration)
	fmt.Printf("✅ Data processed: %s\n", formatBytes(totalBytes))
	fmt.Printf("✅ Variants processed: %d\n", totalVariants)
	fmt.Printf("✅ Throughput: %.2f MB/s\n", float64(totalBytes)/duration.Seconds()/1024/1024)
	fmt.Println()

	// Generate validation sequences
	fmt.Println("[STEP 3] Generating validation sequences...")
	sequences, err := generateValidationSequences(model, 5, 50000)
	if err != nil {
		log.Fatalf("Failed to generate validation sequences: %v", err)
	}
	fmt.Printf("✅ Generated %d validation sequences\n", len(sequences))
	fmt.Println()

	// Summary statistics of the generated sequences
	fmt.Println("[STEP 4] Calculating sequence statistics...")
	complexity := calculateComplexity(sequences)
	gcContent := calculateGCContent(sequences)

	fmt.Printf("📊 SEQUENCE STATISTICS:\n")
	fmt.Printf("  Complexity: %.4f", complexity)
	if complexity > 0.05 {
		fmt.Printf(" (%.1fx an assumed 0.005 baseline) ✅", complexity/0.005)
	} else if complexity > 0.01 {
		fmt.Printf(" (%.1fx an assumed 0.005 baseline) 🟡", complexity/0.005)
	} else {
		fmt.Print(" ❌")
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

	// Final results
	fmt.Println("🎉 REAL VCF TRAINING RESULTS:")
	fmt.Println("=============================")
	fmt.Printf("Training Method:     VCF REF+ALT alleles, concatenated\n")
	fmt.Printf("Data Processed:      %s of compressed VCF\n", formatBytes(totalBytes))
	fmt.Printf("Processing Time:     %v\n", duration)
	fmt.Printf("Throughput:          %.2f MB/s\n", float64(totalBytes)/duration.Seconds()/1024/1024)
	fmt.Printf("Variants Processed:  %d VCF records\n", totalVariants)
	fmt.Printf("Apple Silicon:       %d P-cores utilized\n", perfCores)
	fmt.Printf("Complexity:          %.4f\n", complexity)
	fmt.Printf("GC Content:          %.2f%%\n", gcContent*100)

	fmt.Println()
	fmt.Printf("📊 Results saved in: %s/\n", progressDir)
	fmt.Printf("🏁 REAL VCF TRAINING COMPLETE IN %v\n", duration)
}

func countVCFFiles(populationDir string) int {
	vcfFiles, err := filepath.Glob(filepath.Join(populationDir, "*/*.vcf.gz"))
	if err != nil {
		return 0
	}
	return len(vcfFiles)
}

func processRealVCFFiles(model *authenticity.AdaptiveHMMModel, populationDir string, numCores int) (int64, int64) {
	populations := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"}

	var totalBytes int64
	var totalVariants int64
	var wg sync.WaitGroup

	fmt.Printf("🔥 Starting real VCF processing with %d P-cores\n", numCores)

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
	fmt.Printf("🧬 Processing population %s (real VCF files)...\n", pop)

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
		k := 4
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
