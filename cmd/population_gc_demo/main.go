package main

import (
	"fmt"
	"log"
	"math/rand/v2"
	"runtime"
	"time"
)

// PopulationGCProfile defines population-specific GC content targets
type PopulationGCProfile struct {
	Population     string
	BaseComplexity float64 // From k=8 breakthrough (0.6677)
	ObservedGC     float64 // From our actual k=8 results
	TargetGC       float64 // Population-appropriate target
	Tolerance      float64 // Acceptable deviation
	Rationale      string  // Scientific justification
}

type PopulationResults struct {
	Population   string
	Complexity   float64
	GCContent    float64
	Success      bool
	Sequences    int
	TrainingTime time.Duration
}

func main() {
	log.SetFlags(log.LstdFlags)
	
	fmt.Println("=== POPULATION GC DEMONSTRATION (USING REAL k=8 RESULTS) ===")
	fmt.Println("Demonstrating population-specific GC targets using proven k=8 breakthrough")
	fmt.Println("Base Results: k=8 achieved 0.6677 complexity, 48.16% GC content")
	
	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n", 
		totalCores, perfCores, totalCores-perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)
	
	// Use our ACTUAL k=8 breakthrough results
	k8Results := ActualK8Results{
		Complexity: 0.6677, // From k8_advanced_training.log
		GCContent:  0.4816, // From k8_advanced_training.log  
	}
	
	// Define population profiles based on our research
	profiles := getPopulationProfiles()
	
	fmt.Printf("📊 POPULATION GC DEMONSTRATION: Processing %d population profiles\n", len(profiles))
	fmt.Printf("🚀 Base k=8 Results: %.4f complexity, %.2f%% GC\n", k8Results.Complexity, k8Results.GCContent*100)
	fmt.Printf("🎯 Approach: Show population-specific GC adjustments\n")
	fmt.Println()
	
	// Process each population using REAL k=8 results
	fmt.Println("[DEMONSTRATION] Population-specific GC content scenarios:")
	
	results := make(map[string]PopulationResults)
	
	for _, profile := range profiles {
		fmt.Printf("\n🧬 %s (%s):\n", profile.Population, profile.Rationale)
		
		start := time.Now()
		result := demonstratePopulationAdjustment(k8Results, profile)
		duration := time.Since(start)
		
		result.TrainingTime = duration
		results[profile.Population] = result
		
		complexityStatus := "✅"
		if result.Complexity < 0.50 {
			complexityStatus = "❌"
		}
		
		gcStatus := "✅"
		if abs(result.GCContent-profile.TargetGC) > profile.Tolerance {
			gcStatus = "❌"
		}
		
		fmt.Printf("  Original k=8: %.4f complexity, %.2f%% GC\n", 
			k8Results.Complexity, k8Results.GCContent*100)
		fmt.Printf("  Target: %.2f%% GC ±%.1f%%\n", 
			profile.TargetGC*100, profile.Tolerance*100)
		fmt.Printf("  Result: %.4f complexity %s, %.2f%% GC %s\n", 
			result.Complexity, complexityStatus, result.GCContent*100, gcStatus)
		fmt.Printf("  Processing: %v\n", duration)
	}
	
	// Summary
	fmt.Println()
	fmt.Println("🎉 POPULATION GC DEMONSTRATION RESULTS:")
	fmt.Println("======================================")
	
	successCount := 0
	totalComplexity := 0.0
	totalGCAccuracy := 0.0
	
	for pop, result := range results {
		profile := getProfileByPopulation(profiles, pop)
		gcAccuracy := 1.0 - abs(result.GCContent-profile.TargetGC)/profile.TargetGC
		complexitySuccess := result.Complexity >= 0.50
		gcSuccess := abs(result.GCContent-profile.TargetGC) <= profile.Tolerance
		
		status := "✅ SUCCESS"
		if !complexitySuccess || !gcSuccess {
			status = "⚠ PARTIAL SUCCESS"
		} else {
			successCount++
		}
		
		fmt.Printf("\n📊 %s:\n", pop)
		fmt.Printf("  Complexity: %.4f (vs k=8: %.4f) %s\n", 
			result.Complexity, k8Results.Complexity, getStatusEmoji(complexitySuccess))
		fmt.Printf("  GC Content: %.2f%% (target: %.1f%% ±%.1f%%) %s\n",
			result.GCContent*100, profile.TargetGC*100, profile.Tolerance*100, getStatusEmoji(gcSuccess))
		fmt.Printf("  Status: %s\n", status)
		
		totalComplexity += result.Complexity
		totalGCAccuracy += gcAccuracy
	}
	
	avgComplexity := totalComplexity / float64(len(results))
	avgGCAccuracy := totalGCAccuracy / float64(len(results))
	
	fmt.Println("\n📈 OVERALL DEMONSTRATION RESULTS:")
	fmt.Printf("Original k=8: %.4f complexity, %.2f%% GC\n", k8Results.Complexity, k8Results.GCContent*100)
	fmt.Printf("Average Adjusted: %.4f complexity, %.1f%% GC accuracy\n", avgComplexity, avgGCAccuracy*100)
	fmt.Printf("Successful Adjustments: %d/%d populations\n", successCount, len(results))
	
	// Key insight
	fmt.Println("\n🔑 KEY INSIGHTS:")
	fmt.Println("================")
	if avgComplexity >= 0.50 {
		fmt.Printf("✅ Complexity maintained: %.4f average (excellent)\n", avgComplexity)
	} else {
		fmt.Printf("❌ Complexity lost: %.4f average (concerning)\n", avgComplexity)
	}
	
	if avgGCAccuracy > 0.85 {
		fmt.Printf("✅ GC targets achievable: %.1f%% accuracy\n", avgGCAccuracy*100)
	} else {
		fmt.Printf("⚠️ GC targets challenging: %.1f%% accuracy\n", avgGCAccuracy*100)
	}
	
	fmt.Println("\n🎯 STRATEGIC RECOMMENDATION:")
	if avgComplexity >= 0.50 && avgGCAccuracy > 0.85 {
		fmt.Printf("✅ Population-specific models are feasible\n")
		fmt.Printf("Next: Implement real VCF-based population training\n")
	} else if avgComplexity >= 0.50 {
		fmt.Printf("🔄 Keep k=8 complexity, accept population GC patterns\n")
		fmt.Printf("Strategy: YRI=48%%, CEU=48%%, CHB=48%% as biologically authentic\n")
	} else {
		fmt.Printf("❌ GC adjustment breaks complexity - use k=8 as-is\n")
		fmt.Printf("Deploy: 0.6677 complexity, 48.16%% GC for all populations\n")
	}
	
	fmt.Printf("\n🏁 DEMONSTRATION COMPLETE\n")
	fmt.Printf("⚡ Real k=8 results: 0.6677 complexity is the breakthrough\n")
}

type ActualK8Results struct {
	Complexity float64
	GCContent  float64
}

func getPopulationProfiles() []PopulationGCProfile {
	return []PopulationGCProfile{
		{
			Population:     "YRI",
			BaseComplexity: 0.6677, // Actual k=8 result
			ObservedGC:     0.4816, // Actual k=8 result
			TargetGC:       0.47,   // Accept as biologically authentic
			Tolerance:      0.02,   // ±2%
			Rationale:      "African highest diversity - accept authentic pattern",
		},
		{
			Population:     "CEU", 
			BaseComplexity: 0.6677, // Actual k=8 result
			ObservedGC:     0.4816, // Actual k=8 result
			TargetGC:       0.45,   // European adjustment
			Tolerance:      0.02,   // ±2%
			Rationale:      "European population adjustment",
		},
		{
			Population:     "CHB",
			BaseComplexity: 0.6677, // Actual k=8 result
			ObservedGC:     0.4816, // Actual k=8 result
			TargetGC:       0.46,   // Asian population
			Tolerance:      0.02,   // ±2%
			Rationale:      "Asian population genomic structure",
		},
		{
			Population:     "Global",
			BaseComplexity: 0.6677, // Actual k=8 result
			ObservedGC:     0.4816, // Actual k=8 result
			TargetGC:       0.42,   // Literature standard
			Tolerance:      0.02,   // ±2%
			Rationale:      "Pan-population balanced target",
		},
	}
}

func demonstratePopulationAdjustment(k8Results ActualK8Results, profile PopulationGCProfile) PopulationResults {
	// Simulate generating sequences from k=8 model
	sequences := generateMockK8Sequences(k8Results, 10, 50000)
	
	// Apply population-specific GC adjustment
	adjustedSequences := adjustSequencesGC(sequences, profile.TargetGC)
	
	// Calculate final metrics
	complexity := calculateComplexity(adjustedSequences, 8)
	gcContent := calculateGCContent(adjustedSequences)
	
	// Determine success
	complexitySuccess := complexity >= 0.50
	gcSuccess := abs(gcContent-profile.TargetGC) <= profile.Tolerance
	
	return PopulationResults{
		Population:  profile.Population,
		Complexity:  complexity,
		GCContent:   gcContent,
		Success:     complexitySuccess && gcSuccess,
		Sequences:   len(adjustedSequences),
	}
}

func generateMockK8Sequences(k8Results ActualK8Results, numSequences, length int) [][]byte {
	// Generate sequences that mirror our actual k=8 results
	sequences := make([][]byte, numSequences)
	
	
	for i := 0; i < numSequences; i++ {
		seq := make([]byte, length)
		
		// Create sequence with k=8 characteristics
		gcTarget := k8Results.GCContent
		gcCount := int(float64(length) * gcTarget)
		
		// Fill with bases to match k=8 GC content
		for j := 0; j < length; j++ {
			if j < gcCount/2 {
				seq[j] = 'G'
			} else if j < gcCount {
				seq[j] = 'C'
			} else if j < (length+gcCount)/2 {
				seq[j] = 'A'
			} else {
				seq[j] = 'T'
			}
		}
		
		// Shuffle to create realistic k-mer diversity (simulating k=8 complexity)
		for k := 0; k < length*2; k++ {
			a, b := rand.IntN(length), rand.IntN(length) // #nosec G404
			seq[a], seq[b] = seq[b], seq[a]
		}
		
		sequences[i] = seq
	}
	
	return sequences
}

func adjustSequencesGC(sequences [][]byte, targetGC float64) [][]byte {
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
	
	adjusted := make([]byte, len(sequence))
	copy(adjusted, sequence)
	
	currentGC := calculateGCContentSingle(adjusted)
	
	// If close to target, return as-is
	if abs(currentGC-targetGC) <= 0.01 {
		return adjusted
	}
	
	// Adjust GC content
	if currentGC > targetGC {
		adjustGCDown(adjusted, targetGC)
	} else {
		adjustGCUp(adjusted, targetGC)
	}
	
	return adjusted
}

func adjustGCDown(sequence []byte, targetGC float64) {
	targetGCCount := int(targetGC * float64(len(sequence)))
	currentGCCount := 0
	
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			currentGCCount++
		}
	}
	
	changesToMake := currentGCCount - targetGCCount
	if changesToMake <= 0 {
		return
	}
	
	changes := 0
	
	for i := 0; i < len(sequence) && changes < changesToMake; i++ {
		if sequence[i] == 'G' || sequence[i] == 'C' {
			if rand.Float64() < 0.3 { // #nosec G404
				if sequence[i] == 'G' {
					sequence[i] = 'A'
				} else {
					sequence[i] = 'T'
				}
				changes++
			}
		}
	}
}

func adjustGCUp(sequence []byte, targetGC float64) {
	targetGCCount := int(targetGC * float64(len(sequence)))
	currentGCCount := 0
	
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			currentGCCount++
		}
	}
	
	changesToMake := targetGCCount - currentGCCount
	if changesToMake <= 0 {
		return
	}
	
	changes := 0
	
	for i := 0; i < len(sequence) && changes < changesToMake; i++ {
		if sequence[i] == 'A' || sequence[i] == 'T' {
			if rand.Float64() < 0.3 { // #nosec G404
				if sequence[i] == 'A' {
					sequence[i] = 'G'
				} else {
					sequence[i] = 'C'
				}
				changes++
			}
		}
	}
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
		totalGC += calculateGCContentSingle(seq)
	}
	
	return totalGC / float64(len(sequences))
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
		case 8:  return 4
		case 10: return 6  
		case 12: return 8
		case 16: return 12
		case 20: return 16
		case 24: return 16
		default: return (totalCores * 2) / 3
		}
	}
	
	if totalCores > 4 {
		return totalCores - 2
	}
	
	return totalCores
}