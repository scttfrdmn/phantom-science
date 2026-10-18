// Command complete_11_population_training trains one k=8 adaptive HMM per
// 1000 Genomes population label (11 of them) and writes sample FASTA
// sequences under sequences/validation/. No 1000 Genomes data is read: each
// model is trained on a 50 kb random sequence steered towards a GC target
// hard-coded per population, so the models differ only in GC bias. The
// per-population GC targets and individual counts are hard-coded, and the
// performance-core count is guessed from runtime.NumCPU.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"fmt"
	"log"
	"math/rand/v2"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// PopulationModel represents a population-specific genomic model
type PopulationModel struct {
	Population  string
	Model       *authenticity.AdaptiveHMMModel
	TargetGC    float64
	Complexity  float64
	Trained     bool
	TrainTime   time.Duration
	Region      string
	Individuals int
}

// PopulationTarget defines training targets for each 1000 Genomes population
type PopulationTarget struct {
	Population  string
	TargetGC    float64
	Tolerance   float64
	Region      string
	Individuals int
	Rationale   string
}

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== 11-POPULATION SYNTHETIC-DATA TRAINING ===")
	fmt.Println("k=8 adaptive HMM per population, trained on synthetic sequences")
	fmt.Println("11 1000 Genomes population labels; no 1000 Genomes data is read")

	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 CPU: %d cores (%d assumed performance cores + %d others)\n",
		totalCores, perfCores, totalCores-perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)

	// Define all 11 1000 Genomes populations with research-based targets
	targets := getComplete1000GenomesTargets()

	fmt.Printf("📊 Processing %d populations\n", len(targets))
	fmt.Printf("🚀 Model: k=8 adaptive HMM\n")
	fmt.Printf("🎯 Training data: 50 kb random sequence per population, steered to a hard-coded GC target\n")
	fmt.Println()

	// Display population overview
	fmt.Println("🌍 POPULATION OVERVIEW:")
	displayPopulationOverview(targets)
	fmt.Println()

	// Train each population model
	models := make(map[string]*PopulationModel)

	fmt.Println("🧬 POPULATION MODEL TRAINING:")
	fmt.Println("==============================")

	startTime := time.Now()

	for i, target := range targets {
		fmt.Printf("[%d/%d] Training %s (%s - %d individuals)...\n",
			i+1, len(targets), target.Population, target.Region, target.Individuals)

		modelStart := time.Now()
		model := trainPopulationModel(target)
		modelDuration := time.Since(modelStart)

		model.TrainTime = modelDuration
		models[target.Population] = model

		statusEmoji := "✅"
		if !model.Trained {
			statusEmoji = "❌"
		}

		fmt.Printf("  %s %s: Complexity %.4f, GC %.2f%% (%v)\n",
			statusEmoji, target.Population, model.Complexity, model.TargetGC*100, modelDuration)
	}

	totalTrainingTime := time.Since(startTime)

	// Generate comprehensive results
	fmt.Println()
	fmt.Println("🎉 COMPLETE 11-POPULATION RESULTS:")
	fmt.Println("==================================")

	analyzePopulationResults(models, targets)

	// Regional analysis
	fmt.Println()
	fmt.Println("🌍 REGIONAL ANALYSIS:")
	fmt.Println("====================")
	analyzeRegionalResults(models, targets)

	// Overall system performance
	fmt.Println()
	fmt.Println("📈 OVERALL SYSTEM PERFORMANCE:")
	fmt.Println("==============================")

	successCount := 0
	totalComplexity := 0.0
	totalGCAccuracy := 0.0
	totalIndividuals := 0

	for _, target := range targets {
		model := models[target.Population]
		if model.Trained {
			successCount++
			totalComplexity += model.Complexity
			gcAccuracy := 1.0 - abs(model.TargetGC-target.TargetGC)/target.TargetGC
			totalGCAccuracy += gcAccuracy
			totalIndividuals += target.Individuals
		}
	}

	if successCount > 0 {
		avgComplexity := totalComplexity / float64(successCount)
		avgGCAccuracy := totalGCAccuracy / float64(successCount)

		fmt.Printf("✅ Successful Models: %d/%d (%.1f%%)\n",
			successCount, len(targets), float64(successCount)/float64(len(targets))*100)
		fmt.Printf("📊 Average Complexity: %.4f (earlier k=8 run, hard-coded: 0.6677)\n", avgComplexity)
		fmt.Printf("🎯 Average GC Accuracy: %.1f%%\n", avgGCAccuracy*100)
		fmt.Printf("👥 Total Coverage: %d individuals across %d populations\n",
			totalIndividuals, successCount)
		fmt.Printf("⏱️ Total Training Time: %v (~%.1f min per population)\n",
			totalTrainingTime, totalTrainingTime.Minutes()/float64(len(targets)))
	}

	// Strategic assessment
	fmt.Println()
	fmt.Println("🎯 STRATEGIC ASSESSMENT:")
	fmt.Println("========================")

	if successCount == len(targets) {
		fmt.Printf("✅ All %d population models trained\n", len(targets))
	} else if successCount >= len(targets)*3/4 {
		fmt.Printf("🟡 %d/%d population models trained\n", successCount, len(targets))
		fmt.Printf("🔧 ACTION: Debug failed populations and retry training\n")
		fmt.Printf("📋 Train every population before using the output\n")
	} else {
		fmt.Printf("🔴 Only %d/%d population models trained\n", successCount, len(targets))
		fmt.Printf("🚨 Investigate the training failures above\n")
		fmt.Printf("🔍 ACTION: Review k=8 architecture for population scalability\n")
	}

	// Generate sample sequences for validation
	fmt.Println()
	fmt.Println("🧬 SEQUENCE GENERATION FOR VALIDATION:")
	fmt.Println("======================================")

	generateValidationSequences(models, targets)

	fmt.Printf("\n🏁 COMPLETE 11-POPULATION TRAINING FINISHED\n")
	fmt.Printf("⚡ Total time: %v\n", totalTrainingTime)
	fmt.Printf("🧬 Sample sequences written\n")
	fmt.Printf("📁 Sequences saved to: sequences/validation/\n")
}

func getComplete1000GenomesTargets() []PopulationTarget {
	return []PopulationTarget{
		// African populations (highest diversity)
		{
			Population: "YRI", TargetGC: 0.47, Tolerance: 0.02, Region: "Africa",
			Individuals: 108, Rationale: "Yoruba - highest genetic diversity, accept authentic GC pattern",
		},
		{
			Population: "LWK", TargetGC: 0.465, Tolerance: 0.02, Region: "Africa",
			Individuals: 99, Rationale: "Luhya - East African diversity, high GC functional variants",
		},
		{
			Population: "GWD", TargetGC: 0.468, Tolerance: 0.02, Region: "Africa",
			Individuals: 113, Rationale: "Gambian - West African coastal diversity",
		},

		// European populations
		{
			Population: "CEU", TargetGC: 0.45, Tolerance: 0.02, Region: "Europe",
			Individuals: 99, Rationale: "Utah Europeans - European-specific genomic patterns",
		},
		{
			Population: "GBR", TargetGC: 0.448, Tolerance: 0.02, Region: "Europe",
			Individuals: 91, Rationale: "British - Northern European population structure",
		},
		{
			Population: "TSI", TargetGC: 0.452, Tolerance: 0.02, Region: "Europe",
			Individuals: 107, Rationale: "Tuscan - Southern European Mediterranean patterns",
		},

		// East Asian populations
		{
			Population: "CHB", TargetGC: 0.46, Tolerance: 0.02, Region: "East Asia",
			Individuals: 103, Rationale: "Han Chinese Beijing - East Asian population genomics",
		},
		{
			Population: "JPT", TargetGC: 0.458, Tolerance: 0.02, Region: "East Asia",
			Individuals: 104, Rationale: "Japanese - Island population genomic structure",
		},
		{
			Population: "CHS", TargetGC: 0.462, Tolerance: 0.02, Region: "East Asia",
			Individuals: 105, Rationale: "Han Chinese South - Southern Chinese diversity",
		},

		// American populations (admixed)
		{
			Population: "MXL", TargetGC: 0.445, Tolerance: 0.02, Region: "Americas",
			Individuals: 64, Rationale: "Mexican - Native American and European admixture",
		},
		{
			Population: "PUR", TargetGC: 0.448, Tolerance: 0.02, Region: "Americas",
			Individuals: 104, Rationale: "Puerto Rican - Caribbean admixed population",
		},
	}
}

func displayPopulationOverview(targets []PopulationTarget) {
	regions := make(map[string][]PopulationTarget)
	for _, target := range targets {
		regions[target.Region] = append(regions[target.Region], target)
	}

	regionOrder := []string{"Africa", "Europe", "East Asia", "Americas"}

	for _, region := range regionOrder {
		pops := regions[region]
		if len(pops) == 0 {
			continue
		}

		totalIndividuals := 0
		for _, pop := range pops {
			totalIndividuals += pop.Individuals
		}

		fmt.Printf("🌍 %s (%d populations, %d individuals):\n", region, len(pops), totalIndividuals)

		sort.Slice(pops, func(i, j int) bool {
			return pops[i].Population < pops[j].Population
		})

		for _, pop := range pops {
			fmt.Printf("  • %s: %d individuals, %.1f%% GC target\n",
				pop.Population, pop.Individuals, pop.TargetGC*100)
		}
	}
}

func trainPopulationModel(target PopulationTarget) *PopulationModel {
	// Create k=8 model for this population
	model := authenticity.NewAdaptiveHMMModel(8)

	// Population-specific GC content will be handled during training

	// Train with synthetic population-representative data
	err := trainModelWithSyntheticData(model, target)
	if err != nil {
		log.Printf("Warning: Training failed for %s: %v", target.Population, err)
		return &PopulationModel{
			Population:  target.Population,
			Model:       model,
			TargetGC:    target.TargetGC,
			Region:      target.Region,
			Individuals: target.Individuals,
			Trained:     false,
			Complexity:  0.0,
		}
	}

	// Validate training by generating an actual test sequence
	testResult, err := model.GenerateSequence(10000) // 10kb test sequence
	if err != nil {
		log.Printf("Warning: Sequence generation failed for %s: %v", target.Population, err)
		return &PopulationModel{
			Population:  target.Population,
			Model:       model,
			TargetGC:    target.TargetGC,
			Region:      target.Region,
			Individuals: target.Individuals,
			Trained:     false,
			Complexity:  0.0,
		}
	}

	// Calculate actual complexity from generated sequence
	complexity := calculateSequenceComplexity(testResult.Sequence, 8)

	// Calculate actual GC content from generated sequence
	gcContent := calculateSequenceGCContent(testResult.Sequence)

	populationModel := &PopulationModel{
		Population:  target.Population,
		Model:       model,
		TargetGC:    gcContent, // Use actual measured GC content from generation
		Region:      target.Region,
		Individuals: target.Individuals,
		Trained:     true,
		Complexity:  complexity,
	}

	return populationModel
}

func analyzePopulationResults(models map[string]*PopulationModel, targets []PopulationTarget) {
	// Sort populations by complexity for display
	type PopResult struct {
		Target PopulationTarget
		Model  *PopulationModel
	}

	var results []PopResult
	for _, target := range targets {
		if model, exists := models[target.Population]; exists {
			results = append(results, PopResult{Target: target, Model: model})
		}
	}

	sort.Slice(results, func(i, j int) bool {
		return results[i].Model.Complexity > results[j].Model.Complexity
	})

	for i, result := range results {
		target := result.Target
		model := result.Model

		complexitySuccess := model.Complexity >= 0.50
		gcSuccess := abs(model.TargetGC-target.TargetGC) <= target.Tolerance

		status := "✅ SUCCESS"
		if !complexitySuccess || !gcSuccess || !model.Trained {
			status = "❌ FAILED"
		}

		fmt.Printf("[%2d] %s (%s, %d ind.):\n", i+1, target.Population, target.Region, target.Individuals)
		fmt.Printf("     Complexity: %.4f %s | GC: %.2f%% (target: %.1f%%) %s\n",
			model.Complexity, getStatusEmoji(complexitySuccess),
			model.TargetGC*100, target.TargetGC*100, getStatusEmoji(gcSuccess))
		fmt.Printf("     Training: %v | Status: %s\n", model.TrainTime, status)
		fmt.Println()
	}
}

func analyzeRegionalResults(models map[string]*PopulationModel, targets []PopulationTarget) {
	regions := make(map[string][]PopulationTarget)
	for _, target := range targets {
		regions[target.Region] = append(regions[target.Region], target)
	}

	regionOrder := []string{"Africa", "Europe", "East Asia", "Americas"}

	for _, region := range regionOrder {
		pops := regions[region]
		if len(pops) == 0 {
			continue
		}

		totalComplexity := 0.0
		totalIndividuals := 0
		successCount := 0

		for _, pop := range pops {
			if model, exists := models[pop.Population]; exists && model.Trained {
				totalComplexity += model.Complexity
				totalIndividuals += pop.Individuals
				successCount++
			}
		}

		if successCount > 0 {
			avgComplexity := totalComplexity / float64(successCount)
			fmt.Printf("🌍 %s: %.4f avg complexity (%d/%d populations, %d individuals)\n",
				region, avgComplexity, successCount, len(pops), totalIndividuals)
		}
	}
}

func getStatusEmoji(success bool) string {
	if success {
		return "✅"
	}
	return "❌"
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
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

// trainModelWithSyntheticData trains the model with population-representative synthetic data
func trainModelWithSyntheticData(model *authenticity.AdaptiveHMMModel, target PopulationTarget) error {
	// Create synthetic training sequence based on population characteristics
	// This simulates what would be done with real VCF data

	baseSequenceLength := 50000 // 50kb training sequence

	// Generate population-specific synthetic sequence
	syntheticSequence := generatePopulationSyntheticSequence(target, baseSequenceLength)

	// Train the model on the synthetic sequence
	return model.TrainFromSequence(syntheticSequence)
}

// generatePopulationSyntheticSequence generates a synthetic sequence with population characteristics
func generatePopulationSyntheticSequence(target PopulationTarget, length int) []byte {
	sequence := make([]byte, length)
	bases := []byte{'A', 'T', 'G', 'C'}

	// Set random seed for reproducible population-specific sequences
	rng := rand.New(rand.NewPCG(uint64(hashString(target.Population)), 0)) // #nosec G115 G404

	// Generate sequence with target GC content
	gcTarget := target.TargetGC

	for i := 0; i < length; i++ {
		// Calculate current GC content
		if i > 100 { // After initial 100 bases
			currentGC := calculateGCInWindow(sequence, i, 100)

			// Adjust base selection based on GC target
			if currentGC < gcTarget-0.02 {
				// Need more GC content
				if rng.Float64() < 0.7 {
					sequence[i] = []byte{'G', 'C'}[rng.IntN(2)]
				} else {
					sequence[i] = []byte{'A', 'T'}[rng.IntN(2)]
				}
			} else if currentGC > gcTarget+0.02 {
				// Need more AT content
				if rng.Float64() < 0.7 {
					sequence[i] = []byte{'A', 'T'}[rng.IntN(2)]
				} else {
					sequence[i] = []byte{'G', 'C'}[rng.IntN(2)]
				}
			} else {
				// Within target range
				sequence[i] = bases[rng.IntN(4)]
			}
		} else {
			// Initial bases
			sequence[i] = bases[rng.IntN(4)]
		}
	}

	return sequence
}

// calculateGCInWindow calculates GC content in a window
func calculateGCInWindow(sequence []byte, pos, windowSize int) float64 {
	start := pos - windowSize
	if start < 0 {
		start = 0
	}

	gcCount := 0
	totalBases := 0

	for i := start; i < pos && i < len(sequence); i++ {
		totalBases++
		if sequence[i] == 'G' || sequence[i] == 'C' {
			gcCount++
		}
	}

	if totalBases == 0 {
		return 0.0
	}

	return float64(gcCount) / float64(totalBases)
}

// calculateSequenceComplexity calculates k-mer complexity of a sequence
func calculateSequenceComplexity(sequence []byte, k int) float64 {
	if len(sequence) < k {
		return 0.0
	}

	// Count unique k-mers
	uniqueKmers := make(map[string]bool)
	for i := 0; i <= len(sequence)-k; i++ {
		kmer := string(sequence[i : i+k])
		uniqueKmers[kmer] = true
	}

	totalKmers := len(sequence) - k + 1
	return float64(len(uniqueKmers)) / float64(totalKmers)
}

// calculateSequenceGCContent calculates GC content of a sequence
func calculateSequenceGCContent(sequence []byte) float64 {
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

// hashString creates a simple hash of a string for seeding
func hashString(s string) int {
	hash := 0
	for _, c := range s {
		hash = hash*31 + int(c)
	}
	if hash < 0 {
		hash = -hash
	}
	return hash
}

// generateValidationSequences generates sequences for expert validation
func generateValidationSequences(models map[string]*PopulationModel, targets []PopulationTarget) {
	// Create output directory
	outputDir := "sequences/validation"
	err := os.MkdirAll(outputDir, 0750)
	if err != nil {
		log.Printf("Warning: Failed to create output directory: %v", err)
		return
	}

	lengths := []int{1000, 5000, 10000} // 1kb, 5kb, 10kb
	sequencesPerLength := 10

	totalSequences := 0
	successfulSequences := 0

	for _, target := range targets {
		model, exists := models[target.Population]
		if !exists || !model.Trained {
			fmt.Printf("  ⚠️  %s: Model not available, skipping\n", target.Population)
			continue
		}

		fmt.Printf("  🧬 %s: Generating sequences...\n", target.Population)

		// Create population-specific directory
		popDir := filepath.Join(outputDir, target.Population)
		err := os.MkdirAll(popDir, 0750)
		if err != nil {
			log.Printf("Warning: Failed to create population directory: %v", err)
			continue
		}

		for _, length := range lengths {
			for i := 0; i < sequencesPerLength; i++ {
				totalSequences++

				// Generate sequence
				result, err := model.Model.GenerateSequence(length)
				if err != nil {
					log.Printf("Warning: Failed to generate sequence for %s (length %d): %v", target.Population, length, err)
					continue
				}

				// Create filename
				filename := fmt.Sprintf("%s_L%d_S%02d.fasta", target.Population, length, i+1)
				filepath := filepath.Join(popDir, filename)

				// Write FASTA file
				err = writeFASTASequence(filepath, target.Population, length, i+1, result.Sequence, model)
				if err != nil {
					log.Printf("Warning: Failed to write sequence file %s: %v", filename, err)
					continue
				}

				successfulSequences++
			}
		}

		fmt.Printf("    ✅ Generated %d sequences (1kb, 5kb, 10kb × %d each)\n", len(lengths)*sequencesPerLength, sequencesPerLength)
	}

	fmt.Printf("\n📊 SEQUENCE GENERATION SUMMARY:\n")
	fmt.Printf("✅ Successfully generated: %d/%d sequences\n", successfulSequences, totalSequences)
	fmt.Printf("📁 Output location: %s\n", outputDir)
	fmt.Printf("🔍 Sequences come from models trained on synthetic data\n")
}

// writeFASTASequence writes a generated sequence to a FASTA file with metadata
func writeFASTASequence(filepath, population string, length, sequenceNum int, sequence []byte, model *PopulationModel) error {
	file, err := os.Create(filepath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	// Calculate sequence metrics
	complexity := calculateSequenceComplexity(sequence, 8)
	gcContent := calculateSequenceGCContent(sequence)

	// Write FASTA header with metadata
	header := fmt.Sprintf(">%s_L%d_S%02d population=%s region=%s length=%d complexity=%.4f gc_content=%.3f model_target_gc=%.3f\n",
		population, length, sequenceNum, population, model.Region, len(sequence), complexity, gcContent, model.TargetGC)

	_, err = file.WriteString(header)
	if err != nil {
		return err
	}

	// Write sequence with line breaks every 80 characters
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
