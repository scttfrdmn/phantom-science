// Command streaming_training trains a k=4 parallel HMM by streaming the
// hard-coded reference /Volumes/Public/genomics_training/human/GRCh38_complete.fasta,
// then runs a "population" pass over the *.vcf.gz files under
// /Volumes/Public/genomics_training/populations, generates five 50 kb
// sequences and prints complexity, GC, homopolymer and CpG figures.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
// The population pass does not read the VCF files: for each file it trains on
// five periodic synthetic sequences (processSyntheticVCF), which corrupts the
// reference-trained model rather than adding population signal. The
// "improvement" factor is relative to an assumed 0.005 baseline, and no JS
// divergence is computed.
package main

import (
	"bufio"
	"context"
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

	fmt.Println("=== STREAMING HIGH-PERFORMANCE GENOMICS TRAINING ===")
	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon detected: %d total cores (%d P-cores + %d E-cores)\n",
		totalCores, perfCores, totalCores-perfCores)
	fmt.Printf("🎯 Streaming training experiment (JS divergence is not computed)\n")
	fmt.Printf("⚡ Using %d P-cores with streaming I/O\n", perfCores)
	fmt.Println()

	// Set CPU utilization to performance cores only
	runtime.GOMAXPROCS(perfCores)

	// Configuration
	populationDir := "/Volumes/Public/genomics_training/populations"
	humanGenomeDir := "/Volumes/Public/genomics_training/human"
	progressDir := "streaming_progress"

	// Create progress directory
	if err := os.MkdirAll(progressDir, 0750); err != nil {
		log.Fatalf("Failed to create progress directory: %v", err)
	}

	fmt.Printf("📊 Progress monitoring: %s/\n", progressDir)
	fmt.Printf("🚀 Streaming mode: Processing genome in chunks\n")
	fmt.Println()

	// Create context for cancellation
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	// Initialize streaming parallel model
	fmt.Println("[STEP 1] Initializing streaming parallel model...")
	model := authenticity.NewParallelHMMModel(4)
	defer model.Cleanup()

	fmt.Printf("✓ Streaming model initialized with %d P-cores\n", perfCores)
	fmt.Println()

	// Step 2: Streaming baseline training
	fmt.Println("[STEP 2] Streaming baseline training...")
	baselineStart := time.Now()
	if err := trainBaselineStreaming(ctx, model, humanGenomeDir, progressDir); err != nil {
		log.Fatalf("Failed baseline training: %v", err)
	}
	baselineTime := time.Since(baselineStart)
	fmt.Printf("✅ Baseline training completed in %v\n", baselineTime)
	fmt.Println()

	// Step 3: Streaming population training
	fmt.Println("[STEP 3] STREAMING POPULATION TRAINING...")
	fmt.Println("🧬 For each VCF file: training on synthetic periodic sequences (VCF contents are not read)")
	fmt.Printf("⚡ Apple Silicon acceleration: %d P-cores\n", perfCores)
	fmt.Println()

	populationStart := time.Now()
	if err := trainPopulationStreaming(ctx, model, populationDir, progressDir); err != nil {
		log.Fatalf("Failed population training: %v", err)
	}
	populationTime := time.Since(populationStart)

	totalTime := baselineTime + populationTime
	fmt.Printf("✅ STREAMING TRAINING COMPLETED!\n")
	fmt.Printf("✅ Total training time: %v\n", totalTime)
	fmt.Printf("✅ Baseline time: %v\n", baselineTime)
	fmt.Printf("✅ Population time: %v\n", populationTime)
	fmt.Println()

	// Step 4: Validation
	fmt.Println("[STEP 4] Generating validation sequences...")
	sequences, avgTime, err := generateValidationSequences(model, 5, 50000)
	if err != nil {
		log.Fatalf("Failed validation: %v", err)
	}

	fmt.Printf("✓ Generated %d validation sequences\n", len(sequences))
	fmt.Printf("✓ Average generation time: %v\n", avgTime)
	fmt.Println()

	// Step 5: Authenticity metrics
	fmt.Println("[STEP 5] Calculating sequence statistics...")
	complexity, gcContent, maxHomopolymer, cpgDepletion := calculateAuthenticityMetrics(sequences)

	fmt.Printf("📊 SEQUENCE STATISTICS:\n")
	fmt.Printf("  Complexity: %.4f", complexity)
	if complexity > 0.05 {
		fmt.Printf(" (%.1fx an assumed 0.005 baseline) ✅", complexity/0.005)
	} else {
		fmt.Print(" ❌")
	}
	fmt.Println()

	fmt.Printf("  GC Content: %.2f%%", gcContent*100)
	if gcContent >= 0.40 && gcContent <= 0.44 {
		fmt.Print(" ✅")
	} else {
		fmt.Printf(" (error: %.4f) ❌", abs(gcContent-0.42))
	}
	fmt.Println()

	fmt.Printf("  Max Homopolymer: %d", maxHomopolymer)
	if maxHomopolymer < 15 {
		fmt.Print(" ✅")
	} else {
		fmt.Print(" ❌")
	}
	fmt.Println()

	fmt.Printf("  CpG Depletion: %.3f", cpgDepletion)
	if cpgDepletion >= 0.15 && cpgDepletion <= 0.35 {
		fmt.Print(" ✅")
	} else {
		fmt.Print(" ❌")
	}
	fmt.Println()
	fmt.Println()

	// Final results
	fmt.Println("🚀 STREAMING TRAINING RESULTS:")
	fmt.Println("==============================")
	fmt.Printf("Training Method:     Streaming + Apple Silicon (%d P-cores)\n", perfCores)
	fmt.Printf("Dataset:             GRCh38 reference + synthetic sequences per VCF file\n")
	fmt.Printf("Total Training Time: %v\n", totalTime)
	fmt.Printf("Complexity:          %.4f", complexity)
	if complexity > 0.05 {
		fmt.Printf(" (%.1fx an assumed 0.005 baseline)", complexity/0.005)
	}
	fmt.Println()
	fmt.Printf("GC Content:          %.2f%% (error: %.4f)\n", gcContent*100, abs(gcContent-0.42))

	// Save final results
	saveFinalResults(progressDir, totalTime, complexity, gcContent, perfCores)

	fmt.Println()
	fmt.Printf("📊 Results saved in: %s/\n", progressDir)
	fmt.Printf("🏁 STREAMING TRAINING COMPLETE IN %v\n", totalTime)
}

func trainBaselineStreaming(ctx context.Context, model *authenticity.ParallelHMMModel, humanDir, progressDir string) error {
	referenceFile := filepath.Join(humanDir, "GRCh38_complete.fasta")
	if _, err := os.Stat(referenceFile); os.IsNotExist(err) {
		return fmt.Errorf("reference genome not found: %s", referenceFile)
	}

	fmt.Printf("🧬 Streaming training on reference genome: %s\n", referenceFile)

	// Stream the large FASTA file
	return streamFASTAFile(ctx, model, referenceFile, progressDir)
}

func streamFASTAFile(ctx context.Context, model *authenticity.ParallelHMMModel, filename, progressDir string) error {
	file, err := os.Open(filename) // #nosec G304
	if err != nil {
		return fmt.Errorf("failed to open FASTA file: %w", err)
	}
	defer file.Close()

	fmt.Printf("📖 Streaming large FASTA file: %s\n", filename)

	var processedChunks int64
	var totalBases int64

	scanner := bufio.NewScanner(file)
	buf := make([]byte, 0, 64*1024) // 64KB buffer
	scanner.Buffer(buf, 1024*1024)  // 1MB max token size

	var currentSequence []byte
	chunkSize := 1000000 // 1MB chunks

	for scanner.Scan() {
		select {
		case <-ctx.Done():
			return ctx.Err()
		default:
		}

		line := strings.TrimSpace(scanner.Text())

		if strings.HasPrefix(line, ">") {
			// Process accumulated sequence
			if len(currentSequence) > 0 {
				if err := processSequenceChunks(ctx, model, currentSequence, chunkSize, &processedChunks, &totalBases); err != nil {
					return err
				}
				currentSequence = nil
			}
		} else {
			// Accumulate sequence data
			currentSequence = append(currentSequence, []byte(strings.ToUpper(line))...)
		}

		// Progress reporting
		if processedChunks%100 == 0 && processedChunks > 0 {
			fmt.Printf("📊 Processed %d chunks, %d bases\n", processedChunks, totalBases)
		}
	}

	// Process the last sequence
	if len(currentSequence) > 0 {
		if err := processSequenceChunks(ctx, model, currentSequence, chunkSize, &processedChunks, &totalBases); err != nil {
			return err
		}
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("error reading FASTA file: %w", err)
	}

	fmt.Printf("✅ Streaming complete: %d chunks, %d total bases\n", processedChunks, totalBases)
	return nil
}

func processSequenceChunks(ctx context.Context, model *authenticity.ParallelHMMModel, sequence []byte, chunkSize int, processedChunks, totalBases *int64) error {
	for i := 0; i < len(sequence); i += chunkSize {
		select {
		case <-ctx.Done():
			return ctx.Err()
		default:
		}

		end := i + chunkSize
		if end > len(sequence) {
			end = len(sequence)
		}

		chunk := sequence[i:end]

		// Train on this chunk
		if err := model.TrainFromSequence(chunk); err != nil {
			return fmt.Errorf("failed to train on chunk: %w", err)
		}

		atomic.AddInt64(processedChunks, 1)
		atomic.AddInt64(totalBases, int64(len(chunk)))
	}

	return nil
}

func trainPopulationStreaming(ctx context.Context, model *authenticity.ParallelHMMModel, populationDir, progressDir string) error {
	populations := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"}

	fmt.Printf("📊 Processing %d populations with streaming\n", len(populations))

	var wg sync.WaitGroup
	semaphore := make(chan struct{}, 4) // Limit concurrent populations

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

			if err := processPopulationVCFs(ctx, model, population, populationDir); err != nil {
				fmt.Printf("⚠ Error processing population %s: %v\n", population, err)
			}
		}(pop)
	}

	wg.Wait()
	fmt.Printf("✅ All populations processed\n")
	return nil
}

func processPopulationVCFs(ctx context.Context, model *authenticity.ParallelHMMModel, pop, populationDir string) error {
	popDir := filepath.Join(populationDir, pop)
	fmt.Printf("🧬 Processing population %s...\n", pop)

	vcfFiles, err := filepath.Glob(filepath.Join(popDir, "*.vcf.gz"))
	if err != nil {
		return err
	}

	processedFiles := 0
	for _, vcfFile := range vcfFiles {
		select {
		case <-ctx.Done():
			return ctx.Err()
		default:
		}

		// The VCF file is not read: processSyntheticVCF trains on synthetic
		// periodic sequences instead.
		if err := processSyntheticVCF(ctx, model, vcfFile, pop); err != nil {
			fmt.Printf("⚠ Error processing %s: %v\n", vcfFile, err)
			continue
		}

		processedFiles++
		if processedFiles%10 == 0 {
			fmt.Printf("  📊 %s: %d/%d files processed\n", pop, processedFiles, len(vcfFiles))
		}
	}

	fmt.Printf("✅ Population %s completed: %d files\n", pop, processedFiles)
	return nil
}

func processSyntheticVCF(ctx context.Context, model *authenticity.ParallelHMMModel, vcfFile, pop string) error {
	// Generate synthetic population-specific sequences
	bases := []byte{'A', 'T', 'C', 'G'}

	// Create 5 synthetic sequences per VCF file
	for i := 0; i < 5; i++ {
		sequence := make([]byte, 50000) // 50KB synthetic sequence

		// Generate population-specific patterns
		for j := range sequence {
			// Add slight population bias
			switch pop {
			case "YRI", "LWK", "GWD": // African populations - higher diversity
				sequence[j] = bases[(j*3+i)%4]
			case "CEU", "GBR", "TSI": // European populations
				sequence[j] = bases[(j*2+i)%4]
			case "CHB", "JPT", "CHS": // Asian populations
				sequence[j] = bases[(j*5+i)%4]
			default: // Admixed populations
				sequence[j] = bases[(j*7+i)%4]
			}
		}

		// Train on synthetic sequence
		if err := model.TrainFromSequence(sequence); err != nil {
			return err
		}
	}

	return nil
}

func generateValidationSequences(model *authenticity.ParallelHMMModel, numSequences, length int) ([][]byte, time.Duration, error) {
	sequences := make([][]byte, numSequences)
	totalTime := time.Duration(0)

	for i := 0; i < numSequences; i++ {
		start := time.Now()
		seq, err := model.GenerateSequence(length)
		if err != nil {
			return nil, 0, fmt.Errorf("failed to generate sequence %d: %w", i+1, err)
		}
		duration := time.Since(start)

		sequences[i] = seq
		totalTime += duration
	}

	avgTime := totalTime / time.Duration(numSequences)
	return sequences, avgTime, nil
}

func calculateAuthenticityMetrics(sequences [][]byte) (complexity, gcContent float64, maxHomopolymer int, cpgDepletion float64) {
	if len(sequences) == 0 {
		return 0, 0, 0, 0
	}

	totalComplexity := 0.0
	totalGC := 0.0
	totalMaxHomopolymer := 0
	totalCpGDepletion := 0.0

	for _, seq := range sequences {
		totalComplexity += calculateSequenceComplexity(seq, 4)
		totalGC += calculateGCContent(seq)

		homopolymer := calculateMaxHomopolymer(seq)
		if homopolymer > totalMaxHomopolymer {
			totalMaxHomopolymer = homopolymer
		}

		totalCpGDepletion += calculateCpGDepletion(seq)
	}

	return totalComplexity / float64(len(sequences)),
		totalGC / float64(len(sequences)),
		totalMaxHomopolymer,
		totalCpGDepletion / float64(len(sequences))
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

func calculateGCContent(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0.0
	}

	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}

	return float64(gcCount) / float64(len(sequence))
}

func calculateMaxHomopolymer(sequence []byte) int {
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

func calculateCpGDepletion(sequence []byte) float64 {
	if len(sequence) < 2 {
		return 0.0
	}

	cgCount := 0
	cCount := 0
	gCount := 0

	for i := 0; i < len(sequence)-1; i++ {
		if sequence[i] == 'C' && sequence[i+1] == 'G' {
			cgCount++
		}
		if sequence[i] == 'C' {
			cCount++
		}
		if sequence[i] == 'G' {
			gCount++
		}
	}

	if sequence[len(sequence)-1] == 'C' {
		cCount++
	}
	if sequence[len(sequence)-1] == 'G' {
		gCount++
	}

	if cCount == 0 || gCount == 0 {
		return 0.0
	}

	expected := float64(cCount) * float64(gCount) / float64(len(sequence))
	if expected == 0 {
		return 0.0
	}

	return float64(cgCount) / expected
}

func saveFinalResults(progressDir string, totalTime time.Duration, complexity, gcContent float64, perfCores int) {
	resultsFile := filepath.Join(progressDir, "final_results.txt")

	results := fmt.Sprintf(`STREAMING HIGH-PERFORMANCE GENOMICS TRAINING RESULTS
===================================================

Training Completed: %s
Total Training Time: %v
Apple Silicon P-Cores: %d

SEQUENCE STATISTICS:
- Sequence Complexity: %.4f (%.1fx an assumed 0.005 baseline)
- GC Content: %.2f%% (target: 42%%)

RUN CONFIGURATION:
- Streaming I/O in chunks (memory use not measured)
- Apple Silicon: P-cores only
- Population pass trains on synthetic sequences; VCF contents are not read

Run completed.
`, time.Now().Format("2006-01-02 15:04:05"), totalTime, perfCores,
		complexity, complexity/0.005, gcContent*100)

	if err := os.WriteFile(resultsFile, []byte(results), 0600); err != nil {
		fmt.Printf("⚠ Could not save results file: %v\n", err)
	}
}

func detectAppleSiliconCores() int {
	totalCores := runtime.NumCPU()

	if totalCores >= 8 {
		switch totalCores {
		case 8: // M1: 4P + 4E
			return 4
		case 10: // M1 Pro: 6P + 4E or 8P + 2E
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
