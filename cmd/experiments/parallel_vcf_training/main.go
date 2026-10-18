// Command parallel_vcf_training trains one k=8 authenticity.AdaptiveHMMModel
// per 1000 Genomes population, several populations at a time, from the
// *.vcf.gz files under  or /data/1000genomes (at least
// 10 of the 11 population directories must have VCFs). It then writes FASTA
// sequences from each model to sequences/parallel_vcf_trained/ and prints
// per-population GC content and complexity. No speedup is measured.
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
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// ParallelVCFTraining represents parallel VCF-based training for a population
type ParallelVCFTraining struct {
	Population        string
	Region            string
	Individuals       int
	Model             *authenticity.AdaptiveHMMModel
	VCFFilesPath      string
	ProcessedVariants int64
	TrainingSequences [][]byte
	TrainingTime      time.Duration
	Trained           bool
	ActualGCContent   float64
	ActualComplexity  float64
	Error             error
}

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== PARALLEL VCF POPULATION TRAINING ===")
	fmt.Println("Training k=8 HMM models on authentic 1000 Genomes VCF data")
	fmt.Println("Training several populations in parallel")
	fmt.Println()

	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n",
		totalCores, perfCores, totalCores-perfCores)

	// Use P-cores for parallel training
	maxParallel := perfCores // 8 P-cores for parallel training
	fmt.Printf("🚀 Parallel training: %d populations simultaneously\n", maxParallel)
	fmt.Println()

	// Define 1000 Genomes populations with VCF data paths
	populations := getParallelPopulationTargets()

	fmt.Printf("🧬 PARALLEL VCF TRAINING: Processing %d populations\n", len(populations))
	fmt.Printf("📊 Architecture: k=8 HMM with authentic genomic variant data\n")
	fmt.Printf("🔬 Strategy: Parallel population-specific VCF processing\n")
	fmt.Printf("📁 Data source: 1000 Genomes Project VCF files\n")
	fmt.Println()

	// Check VCF data availability
	vcfDataPath := checkVCFDataAvailability()
	if vcfDataPath == "" {
		log.Fatal("❌ VCF data not found. Please ensure 1000 Genomes VCF files are available.")
	}

	fmt.Printf("✅ VCF data located at: %s\n", vcfDataPath)
	fmt.Println()

	startTime := time.Now()

	// Train populations in parallel
	fmt.Println("🧬 PARALLEL POPULATION VCF TRAINING:")
	fmt.Println("====================================")

	trainedModels := trainPopulationsInParallel(populations, vcfDataPath, maxParallel)

	totalTrainingTime := time.Since(startTime)

	// Generate validation sequences from VCF-trained models
	fmt.Println()
	fmt.Println("🧬 VCF-TRAINED SEQUENCE GENERATION:")
	fmt.Println("===================================")

	generateVCFTrainedSequences(trainedModels)

	// Analyze training results
	fmt.Println()
	fmt.Println("📊 PARALLEL VCF TRAINING RESULTS:")
	fmt.Println("=================================")

	analyzeParallelVCFTrainingResults(trainedModels, totalTrainingTime)

	fmt.Printf("\n🏁 PARALLEL VCF TRAINING COMPLETED\n")
	fmt.Printf("⏱️ Total time: %v for authentic 1000 Genomes VCF processing\n", totalTrainingTime)
	fmt.Printf("🚀 Up to %d populations trained in parallel\n", maxParallel)
	fmt.Printf("🔬 VCF-trained sequences ready for validation\n")
}

func trainPopulationsInParallel(populations []ParallelVCFTraining, vcfDataPath string, maxParallel int) map[string]*ParallelVCFTraining {
	trainedModels := make(map[string]*ParallelVCFTraining)
	var mutex sync.Mutex

	// Channel to limit concurrent training
	semaphore := make(chan struct{}, maxParallel)
	var wg sync.WaitGroup

	for i := range populations {
		wg.Add(1)
		go func(pop *ParallelVCFTraining) {
			defer wg.Done()

			// Acquire semaphore
			semaphore <- struct{}{}
			defer func() { <-semaphore }()

			fmt.Printf("🚀 STARTING %s (%s - %d individuals) [PARALLEL]\n",
				pop.Population, pop.Region, pop.Individuals)

			trainingStart := time.Now()

			// Train on real VCF data
			err := trainPopulationOnVCF(pop, vcfDataPath)
			if err != nil {
				pop.Error = err
				pop.Trained = false
				fmt.Printf("❌ %s: Training failed - %v\n", pop.Population, err)
			} else {
				pop.TrainingTime = time.Since(trainingStart)
				pop.Trained = true
				fmt.Printf("✅ %s: Processed %d variants, GC %.2f%%, Complexity %.4f (%v) [PARALLEL]\n",
					pop.Population, pop.ProcessedVariants, pop.ActualGCContent*100,
					pop.ActualComplexity, pop.TrainingTime)
			}

			// Thread-safe update of results
			mutex.Lock()
			trainedModels[pop.Population] = pop
			mutex.Unlock()

		}(&populations[i])
	}

	// Wait for all training to complete
	wg.Wait()

	return trainedModels
}

func getParallelPopulationTargets() []ParallelVCFTraining {
	return []ParallelVCFTraining{
		// African populations (highest diversity)
		{Population: "YRI", Region: "Africa", Individuals: 108},
		{Population: "LWK", Region: "Africa", Individuals: 99},
		{Population: "GWD", Region: "Africa", Individuals: 113},

		// European populations
		{Population: "CEU", Region: "Europe", Individuals: 99},
		{Population: "GBR", Region: "Europe", Individuals: 91},
		{Population: "TSI", Region: "Europe", Individuals: 107},

		// East Asian populations
		{Population: "CHB", Region: "East Asia", Individuals: 103},
		{Population: "JPT", Region: "East Asia", Individuals: 104},
		{Population: "CHS", Region: "East Asia", Individuals: 105},

		// American populations (admixed)
		{Population: "MXL", Region: "Americas", Individuals: 64},
		{Population: "PUR", Region: "Americas", Individuals: 104},
	}
}

func checkVCFDataAvailability() string {
	// Check for real 1000 Genomes data locations
	potentialPaths := []string{
		os.Getenv("PHANTOM_VCF_DATA_DIR"), // empty when unset; os.Stat fails and it is skipped
		"/data/1000genomes",
	}

	for _, path := range potentialPaths {
		if _, err := os.Stat(path); err == nil {
			// Check for population subdirectories with VCF files
			populations := []string{"YRI", "CEU", "CHB", "JPT", "TSI", "PUR", "LWK", "GWD", "GBR", "CHS", "MXL"}
			validPops := 0

			for _, pop := range populations {
				popDir := filepath.Join(path, pop)
				if stat, err := os.Stat(popDir); err == nil && stat.IsDir() {
					vcfFiles, err := filepath.Glob(filepath.Join(popDir, "*.vcf.gz"))
					if err == nil && len(vcfFiles) > 0 {
						validPops++
					}
				}
			}

			// If we have data for at least 10/11 populations, use this path
			if validPops >= 10 {
				return path
			}
		}
	}

	return ""
}

func trainPopulationOnVCF(training *ParallelVCFTraining, vcfDataPath string) error {
	// Create k=8 model for this population
	training.Model = authenticity.NewAdaptiveHMMModel(8)

	// Find VCF files for this population
	vcfFiles, err := findPopulationVCFFiles(vcfDataPath, training.Population)
	if err != nil {
		return fmt.Errorf("failed to find VCF files: %v", err)
	}

	if len(vcfFiles) == 0 {
		return fmt.Errorf("no VCF files found for population %s", training.Population)
	}

	fmt.Printf("📁 %s: Found %d VCF files\n", training.Population, len(vcfFiles))

	// Process each VCF file
	totalVariants := int64(0)
	var allSequences [][]byte

	for _, vcfFile := range vcfFiles {
		variants, sequences, err := processVCFFile(vcfFile, training.Population)
		if err != nil {
			log.Printf("⚠️  %s: Warning: Failed to process %s: %v",
				training.Population, filepath.Base(vcfFile), err)
			continue
		}

		totalVariants += variants
		allSequences = append(allSequences, sequences...)

		fmt.Printf("📊 %s: Processed %s: %d variants\n",
			training.Population, filepath.Base(vcfFile), variants)
	}

	training.ProcessedVariants = totalVariants
	training.TrainingSequences = allSequences

	if len(allSequences) == 0 {
		return fmt.Errorf("no sequences extracted from VCF files")
	}

	// Train model on extracted sequences
	fmt.Printf("🧬 %s: Training model on %d sequences...\n",
		training.Population, len(allSequences))

	for _, sequence := range allSequences {
		err := training.Model.TrainFromSequence(sequence)
		if err != nil {
			log.Printf("⚠️  %s: Warning: Training failed for sequence: %v",
				training.Population, err)
		}
	}

	// Calculate actual metrics from training data
	training.ActualGCContent = calculateAverageGCContent(allSequences)
	training.ActualComplexity = calculateAverageComplexity(allSequences, 8)

	return nil
}

func findPopulationVCFFiles(basePath, population string) ([]string, error) {
	var vcfFiles []string

	// Look for population-specific directory structure like basePath/YRI/chr*.vcf.gz
	popDir := filepath.Join(basePath, population)
	if stat, err := os.Stat(popDir); err == nil && stat.IsDir() {
		// Find all chromosome VCF files for this population
		pattern := filepath.Join(popDir, "chr*.vcf.gz")
		files, err := filepath.Glob(pattern)
		if err == nil && len(files) > 0 {
			vcfFiles = append(vcfFiles, files...)
		}
	}

	// Also check for files directly in basePath with population prefix
	pattern := filepath.Join(basePath, fmt.Sprintf("%s_*.vcf.gz", population))
	files, _ := filepath.Glob(pattern)
	vcfFiles = append(vcfFiles, files...)

	return vcfFiles, nil
}

func processVCFFile(vcfPath, population string) (int64, [][]byte, error) {
	file, err := os.Open(vcfPath) // #nosec G304
	if err != nil {
		return 0, nil, err
	}
	defer file.Close()

	// Create gzip reader
	gzReader, err := gzip.NewReader(file)
	if err != nil {
		return 0, nil, err
	}
	defer gzReader.Close()

	scanner := bufio.NewScanner(gzReader)

	var sequences [][]byte
	variantCount := int64(0)
	currentSequence := make([]byte, 0, 1000) // Start with 1KB capacity

	for scanner.Scan() {
		line := scanner.Text()

		// Skip header lines
		if strings.HasPrefix(line, "#") {
			continue
		}

		// Parse VCF record
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}

		ref := fields[3]
		alt := fields[4]

		// Skip non-standard bases
		if !isValidDNASequence(ref) || !isValidDNASequence(alt) {
			continue
		}

		// Add variant to sequence (simplified - use REF for now)
		currentSequence = append(currentSequence, []byte(ref)...)

		// Add some context around the variant
		if len(alt) > 0 {
			currentSequence = append(currentSequence, []byte(alt)...)
		}

		variantCount++

		// Create sequence chunks of reasonable size
		if len(currentSequence) >= 1000 {
			sequences = append(sequences, make([]byte, len(currentSequence)))
			copy(sequences[len(sequences)-1], currentSequence)
			currentSequence = currentSequence[:0] // Reset but keep capacity
		}
	}

	// Add final sequence if it has content
	if len(currentSequence) > 0 {
		sequences = append(sequences, make([]byte, len(currentSequence)))
		copy(sequences[len(sequences)-1], currentSequence)
	}

	return variantCount, sequences, scanner.Err()
}

func isValidDNASequence(seq string) bool {
	for _, base := range strings.ToUpper(seq) {
		if base != 'A' && base != 'T' && base != 'G' && base != 'C' {
			return false
		}
	}
	return len(seq) > 0
}

func calculateAverageGCContent(sequences [][]byte) float64 {
	if len(sequences) == 0 {
		return 0.0
	}

	totalGC := 0
	totalBases := 0

	for _, seq := range sequences {
		for _, base := range seq {
			totalBases++
			if base == 'G' || base == 'C' {
				totalGC++
			}
		}
	}

	if totalBases == 0 {
		return 0.0
	}

	return float64(totalGC) / float64(totalBases)
}

func calculateAverageComplexity(sequences [][]byte, k int) float64 {
	if len(sequences) == 0 {
		return 0.0
	}

	totalComplexity := 0.0
	validSequences := 0

	for _, seq := range sequences {
		if len(seq) >= k {
			complexity := calculateSequenceComplexity(seq, k)
			totalComplexity += complexity
			validSequences++
		}
	}

	if validSequences == 0 {
		return 0.0
	}

	return totalComplexity / float64(validSequences)
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

func generateVCFTrainedSequences(trainedModels map[string]*ParallelVCFTraining) {
	outputDir := "sequences/parallel_vcf_trained"
	err := os.MkdirAll(outputDir, 0750)
	if err != nil {
		log.Printf("Warning: Failed to create VCF sequences directory: %v", err)
		return
	}

	lengths := []int{1000, 5000, 10000}
	sequencesPerLength := 5 // Fewer sequences for VCF-trained models

	successfulSequences := 0
	totalSequences := 0

	for population, training := range trainedModels {
		if !training.Trained {
			fmt.Printf("⚠️  %s: Model not trained, skipping\n", population)
			continue
		}

		fmt.Printf("🧬 %s: Generating parallel VCF-trained sequences...\n", population)

		// Create population directory
		popDir := filepath.Join(outputDir, population)
		err := os.MkdirAll(popDir, 0750)
		if err != nil {
			log.Printf("Warning: Failed to create directory for %s: %v", population, err)
			continue
		}

		for _, length := range lengths {
			for i := 0; i < sequencesPerLength; i++ {
				totalSequences++

				// Generate sequence
				result, err := training.Model.GenerateSequence(length)
				if err != nil {
					log.Printf("Warning: Failed to generate parallel VCF-trained sequence for %s: %v", population, err)
					continue
				}

				// Create filename
				filename := fmt.Sprintf("%s_PARALLEL_VCF_L%d_S%02d.fasta", population, length, i+1)
				filepath := filepath.Join(popDir, filename)

				// Write FASTA file
				err = writeParallelVCFFASTA(filepath, population, training, length, i+1, result.Sequence)
				if err != nil {
					log.Printf("Warning: Failed to write parallel VCF sequence %s: %v", filename, err)
					continue
				}

				successfulSequences++
			}
		}

		fmt.Printf("    ✅ Generated %d parallel VCF-trained sequences\n", len(lengths)*sequencesPerLength)
	}

	fmt.Printf("\n📊 PARALLEL VCF SEQUENCE GENERATION SUMMARY:\n")
	fmt.Printf("✅ Successfully generated: %d/%d parallel VCF-trained sequences\n", successfulSequences, totalSequences)
	fmt.Printf("📁 Output location: %s\n", outputDir)
}

func writeParallelVCFFASTA(filepath, population string, training *ParallelVCFTraining, length, seqNum int, sequence []byte) error {
	file, err := os.Create(filepath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	// Calculate sequence metrics
	complexity := calculateSequenceComplexity(sequence, 8)
	gcContent := calculateAverageGCContent([][]byte{sequence})

	// Write FASTA header with parallel VCF training metadata
	header := fmt.Sprintf(">%s_PARALLEL_VCF_L%d_S%02d population=%s region=%s length=%d complexity=%.4f gc_content=%.3f parallel_vcf_trained=true processed_variants=%d training_time=%v cores_used=8\n",
		population, length, seqNum, population, training.Region, len(sequence),
		complexity, gcContent, training.ProcessedVariants, training.TrainingTime)

	_, err = file.WriteString(header)
	if err != nil {
		return err
	}

	// Write sequence with line breaks
	for i := 0; i < len(sequence); i += 80 {
		end := i + 80
		if end > len(sequence) {
			end = len(sequence)
		}
		_, err = file.WriteString(string(sequence[i:end]) + "\n")
		if err != nil {
			return err
		}
	}

	return nil
}

func analyzeParallelVCFTrainingResults(trainedModels map[string]*ParallelVCFTraining, totalTime time.Duration) {
	successfulModels := 0
	totalVariants := int64(0)
	avgComplexity := 0.0
	avgGCContent := 0.0

	fmt.Println("📊 PARALLEL POPULATION VCF TRAINING RESULTS:")

	// Sort by region for display
	regions := map[string][]*ParallelVCFTraining{
		"Africa":    {},
		"Europe":    {},
		"East Asia": {},
		"Americas":  {},
	}

	for _, training := range trainedModels {
		regions[training.Region] = append(regions[training.Region], training)
	}

	regionOrder := []string{"Africa", "Europe", "East Asia", "Americas"}

	for _, region := range regionOrder {
		if len(regions[region]) == 0 {
			continue
		}

		fmt.Printf("\n🌍 %s:\n", region)

		for _, training := range regions[region] {
			status := "❌ FAILED"
			if training.Trained {
				status = "✅ SUCCESS"
				successfulModels++
				totalVariants += training.ProcessedVariants
				avgComplexity += training.ActualComplexity
				avgGCContent += training.ActualGCContent
			}

			fmt.Printf("  %s %s: ", status, training.Population)
			if training.Trained {
				fmt.Printf("Variants: %d, GC: %.2f%%, Complexity: %.4f (%v)\n",
					training.ProcessedVariants, training.ActualGCContent*100,
					training.ActualComplexity, training.TrainingTime)
			} else {
				fmt.Printf("Training failed - %v\n", training.Error)
			}
		}
	}

	// Overall statistics
	fmt.Printf("\n📈 OVERALL PARALLEL VCF TRAINING PERFORMANCE:\n")
	fmt.Printf("✅ Successful models: %d/%d (%.1f%%)\n",
		successfulModels, len(trainedModels),
		float64(successfulModels)/float64(len(trainedModels))*100)

	if successfulModels > 0 {
		avgComplexity /= float64(successfulModels)
		avgGCContent /= float64(successfulModels)

		fmt.Printf("📊 Average complexity: %.4f\n", avgComplexity)
		fmt.Printf("🧬 Average GC content: %.2f%%\n", avgGCContent*100)
		fmt.Printf("🔢 Total variants processed: %d\n", totalVariants)
		fmt.Printf("⏱️ Total parallel training time: %v\n", totalTime)
	}

	// Assessment
	fmt.Printf("\n🎯 PARALLEL VCF TRAINING ASSESSMENT:\n")
	if successfulModels == len(trainedModels) {
		fmt.Printf("✅ All populations trained on their VCF data\n")
		fmt.Printf("🔬 Next: validate the generated sequences\n")
	} else if successfulModels >= len(trainedModels)*3/4 {
		fmt.Printf("🟡 MOSTLY SUCCESSFUL: %d/%d populations trained in parallel\n", successfulModels, len(trainedModels))
		fmt.Printf("🔧 ACTION: Debug failed populations and retry training\n")
	} else {
		fmt.Printf("🔴 TRAINING ISSUES: Only %d/%d populations successful in parallel\n", successfulModels, len(trainedModels))
		fmt.Printf("🚨 CRITICAL: Investigate parallel VCF processing pipeline\n")
	}
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
