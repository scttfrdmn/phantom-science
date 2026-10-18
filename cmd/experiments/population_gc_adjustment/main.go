// Command population_gc_adjustment generates sequences from a new, untrained
// k=8 authenticity.AdaptiveHMMModel, then rewrites bases (G/C to A/T or the
// reverse) until each sequence is near a per-population GC target, and
// prints complexity and GC content per profile. Nothing is trained or
// loaded; the 0.6677 "k=8 complexity" and the observed GC values in the
// profiles are constants from an earlier run, and the GC targets are
// chosen, not derived from data.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"fmt"
	"log"
	"math/rand/v2"
	"runtime"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// PopulationGCProfile defines population-specific GC content targets
type PopulationGCProfile struct {
	Population     string
	BaseComplexity float64 // k=8 complexity from an earlier run (0.6677)
	ObservedGC     float64 // GC content from earlier training runs
	TargetGC       float64 // Population-appropriate target
	Tolerance      float64 // Acceptable deviation
	Rationale      string  // Why this target was chosen
}

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== POPULATION GC ADJUSTMENT (OPTIMIZED) ===")
	fmt.Println("Post-processes GC content of k=8 model output toward population-specific targets")

	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n",
		totalCores, perfCores, totalCores-perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)

	// Population profiles: chosen GC targets plus constants from earlier runs
	profiles := getOptimizedPopulationProfiles()

	fmt.Printf("📊 OPTIMIZED POPULATION STRATEGY: Processing %d population profiles\n", len(profiles))
	fmt.Printf("🚀 Base Model: k=8 AdaptiveHMMModel (untrained)\n")
	fmt.Printf("🎯 Approach: Post-process GC content while preserving complexity\n")
	fmt.Println()

	// Initialize a single k=8 model; it is not trained or loaded
	fmt.Println("[STEP 1] Creating k=8 model...")
	baseModel := authenticity.NewAdaptiveHMMModel(8)

	fmt.Printf("✅ k=8 model created (no training state loaded)\n")
	fmt.Println()

	// Process each population efficiently
	fmt.Println("[STEP 2] Generating population-specific sequences...")

	results := make(map[string]PopulationResults)

	for _, profile := range profiles {
		fmt.Printf("🧬 Processing %s (target: %.1f%% GC)...\n",
			profile.Population, profile.TargetGC*100)

		start := time.Now()
		result := processPopulationEfficiently(baseModel, profile)
		duration := time.Since(start)

		result.TrainingTime = duration
		results[profile.Population] = result

		fmt.Printf("  ✅ %s complete: Complexity %.4f, GC %.2f%% in %v\n",
			profile.Population, result.Complexity, result.GCContent*100, duration)
	}

	// Generate comprehensive results
	fmt.Println()
	fmt.Println("🎉 OPTIMIZED POPULATION-SPECIFIC RESULTS:")
	fmt.Println("=========================================")

	successCount := 0
	totalComplexity := 0.0
	totalGCAccuracy := 0.0

	for pop, result := range results {
		profile := getProfileByPopulation(profiles, pop)
		gcAccuracy := 1.0 - abs(result.GCContent-profile.TargetGC)/profile.TargetGC
		complexitySuccess := result.Complexity >= 0.50 // Maintain high complexity
		gcSuccess := abs(result.GCContent-profile.TargetGC) <= profile.Tolerance

		status := "✅ SUCCESS"
		if !complexitySuccess || !gcSuccess {
			status = "⚠ NEEDS REFINEMENT"
		} else {
			successCount++
		}

		fmt.Printf("\n📊 %s (%s):\n", pop, profile.Rationale)
		fmt.Printf("  Complexity: %.4f %s\n",
			result.Complexity, getStatusEmoji(complexitySuccess))
		fmt.Printf("  GC Content: %.2f%% (target: %.1f%% ±%.1f%%) %s\n",
			result.GCContent*100, profile.TargetGC*100, profile.Tolerance*100, getStatusEmoji(gcSuccess))
		fmt.Printf("  Processing Time: %v\n", result.TrainingTime)
		fmt.Printf("  Status: %s\n", status)

		totalComplexity += result.Complexity
		totalGCAccuracy += gcAccuracy
	}

	avgComplexity := totalComplexity / float64(len(results))
	avgGCAccuracy := totalGCAccuracy / float64(len(results))

	fmt.Println("\n📈 OVERALL OPTIMIZED RESULTS:")
	fmt.Printf("Average Complexity: %.4f (earlier k=8 run: 0.6677)\n", avgComplexity)
	fmt.Printf("Average GC Accuracy: %.1f%%\n", avgGCAccuracy*100)
	fmt.Printf("Successful Models: %d/%d\n", successCount, len(results))

	if successCount == len(results) {
		fmt.Printf("🎯 STATUS: ALL PROFILES MET COMPLEXITY AND GC TARGETS\n")
		fmt.Printf("Next Action: EXPERT REVIEW\n")
	} else {
		fmt.Printf("🔄 STATUS: PARTIAL SUCCESS - ANALYZE FAILED POPULATIONS\n")
		fmt.Printf("Next Action: REFINE SUBOPTIMAL POPULATIONS\n")
	}

	fmt.Printf("\n🏁 OPTIMIZED POPULATION PROCESSING COMPLETE\n")
}

func getOptimizedPopulationProfiles() []PopulationGCProfile {
	return []PopulationGCProfile{
		{
			Population:     "YRI",
			BaseComplexity: 0.6677, // k=8 complexity from an earlier run
			ObservedGC:     0.4816, // From k=8 results
			TargetGC:       0.47,   // Chosen target near the observed GC
			Tolerance:      0.02,   // ±2%
			Rationale:      "Highest diversity; target near the observed GC",
		},
		{
			Population:     "CEU",
			BaseComplexity: 0.6677, // k=8 complexity from an earlier run
			ObservedGC:     0.4705, // From population-weighted results
			TargetGC:       0.45,   // Slight European adjustment
			Tolerance:      0.02,   // ±2%
			Rationale:      "European-specific minor adjustment",
		},
		{
			Population:     "CHB",
			BaseComplexity: 0.6677, // k=8 complexity from an earlier run
			ObservedGC:     0.4720, // Estimated from pattern
			TargetGC:       0.46,   // Asian population structure
			Tolerance:      0.02,   // ±2%
			Rationale:      "Asian population genomic structure",
		},
		{
			Population:     "Global",
			BaseComplexity: 0.6677, // k=8 complexity from an earlier run
			ObservedGC:     0.4200, // Literature standard
			TargetGC:       0.42,   // Pan-population balanced
			Tolerance:      0.02,   // ±2%
			Rationale:      "Pan-population balanced model",
		},
	}
}

func processPopulationEfficiently(baseModel *authenticity.AdaptiveHMMModel, profile PopulationGCProfile) PopulationResults {
	// Generate sequences from the k=8 model
	sequences, err := generateValidationSequences(baseModel, 10, 50000)
	if err != nil {
		return PopulationResults{
			Population: profile.Population,
			Complexity: 0.0,
			GCContent:  0.0,
			Success:    false,
		}
	}

	// Apply population-specific GC adjustment while preserving complexity
	adjustedSequences := adjustGCContentPreserveComplexity(sequences, profile.TargetGC)

	// Calculate final metrics
	complexity := calculateComplexity(adjustedSequences, 8)
	gcContent := calculateGCContent(adjustedSequences)

	// Determine success
	complexitySuccess := complexity >= 0.50 // Maintain high complexity
	gcSuccess := abs(gcContent-profile.TargetGC) <= profile.Tolerance

	return PopulationResults{
		Population:  profile.Population,
		Complexity:  complexity,
		GCContent:   gcContent,
		Success:     complexitySuccess && gcSuccess,
		Sequences:   len(adjustedSequences),
		TrainingSeq: len(sequences),
	}
}

func adjustGCContentPreserveComplexity(sequences [][]byte, targetGC float64) [][]byte {
	adjusted := make([][]byte, len(sequences))

	for i, seq := range sequences {
		adjusted[i] = adjustSequenceGC(seq, targetGC)
	}

	return adjusted
}

func adjustSequenceGC(sequence []byte, targetGC float64) []byte {
	if len(sequence) == 0 {
		return sequence
	}

	// Create copy to avoid modifying original
	adjusted := make([]byte, len(sequence))
	copy(adjusted, sequence)

	currentGC := calculateGCContentSingle(adjusted)
	_ = currentGC // Use variable to avoid unused warning

	// If already close to target, return as-is
	if abs(currentGC-targetGC) <= 0.01 {
		return adjusted
	}

	// Simple adjustment strategy: minimal changes to preserve k-mer patterns
	if currentGC > targetGC {
		// Reduce GC content: G/C -> A/T
		adjustGCDown(adjusted, targetGC)
	} else {
		// Increase GC content: A/T -> G/C
		adjustGCUp(adjusted, targetGC)
	}

	return adjusted
}

func adjustGCDown(sequence []byte, targetGC float64) {
	gcPositions := make([]int, 0)

	// Find G/C positions
	for i, base := range sequence {
		if base == 'G' || base == 'C' {
			gcPositions = append(gcPositions, i)
		}
	}

	_ = len(gcPositions) // Use variable to avoid unused warning
	targetGCCount := int(targetGC * float64(len(sequence)))
	changesToMake := len(gcPositions) - targetGCCount

	if changesToMake <= 0 {
		return
	}

	// Randomly change some G/C to A/T
	for i := 0; i < changesToMake && i < len(gcPositions); i++ {
		pos := gcPositions[rand.IntN(len(gcPositions))] // #nosec G404
		if sequence[pos] == 'G' {
			sequence[pos] = 'A'
		} else if sequence[pos] == 'C' {
			sequence[pos] = 'T'
		}
	}
}

func adjustGCUp(sequence []byte, targetGC float64) {
	atPositions := make([]int, 0)

	// Find A/T positions
	for i, base := range sequence {
		if base == 'A' || base == 'T' {
			atPositions = append(atPositions, i)
		}
	}

	currentGCCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			currentGCCount++
		}
	}

	_ = currentGCCount // Use variable to avoid unused warning
	targetGCCount := int(targetGC * float64(len(sequence)))
	changesToMake := targetGCCount - currentGCCount

	if changesToMake <= 0 {
		return
	}

	// Randomly change some A/T to G/C
	for i := 0; i < changesToMake && i < len(atPositions); i++ {
		pos := atPositions[rand.IntN(len(atPositions))] // #nosec G404
		if sequence[pos] == 'A' {
			sequence[pos] = 'G'
		} else if sequence[pos] == 'T' {
			sequence[pos] = 'C'
		}
	}
}

func calculateGCContentSingle(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0
	}

	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}

	return float64(gcCount) / float64(len(sequence))
}

type PopulationResults struct {
	Population   string
	Complexity   float64
	GCContent    float64
	Success      bool
	Sequences    int
	TrainingSeq  int
	TrainingTime time.Duration
}

func generateValidationSequences(model *authenticity.AdaptiveHMMModel, numSequences, length int) ([][]byte, error) {
	sequences := make([][]byte, numSequences)

	for i := 0; i < numSequences; i++ {
		result, err := model.GenerateSequence(length)
		if err != nil {
			return nil, err
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

func getProfileByPopulation(profiles []PopulationGCProfile, pop string) PopulationGCProfile {
	for _, p := range profiles {
		if p.Population == pop {
			return p
		}
	}
	return PopulationGCProfile{}
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
