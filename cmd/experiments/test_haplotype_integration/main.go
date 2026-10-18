// Command test_haplotype_integration runs the Python haplotype simulator
// (from the current directory, which must be the repository root with
// phantom_science_env) for several populations on chr22 and prints variant
// counts, LD blocks, r², diversity, recombination-map length and block-size
// statistics.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"fmt"
	"log"
	"os"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/population"
)

func main() {
	fmt.Println("🧬 TESTING HAPLOTYPE SIMULATION INTEGRATION")
	fmt.Println("=" + strings.Repeat("=", 47))

	// Get project root
	projectRoot, err := os.Getwd()
	if err != nil {
		log.Fatalf("Failed to get project root: %v", err)
	}

	// Create haplotype simulator
	simulator := population.NewHaplotypeSimulator(projectRoot)

	// Test 1: Single population haplotype simulation
	fmt.Println("\n📍 Testing Single Population Haplotype Simulation (YRI)...")
	randomSeed := 42
	yriResult, err := simulator.SimulateHaplotypes(
		population.YRI,
		"chr22",
		30,         // sample size
		10_000_000, // 10Mb sequence
		&randomSeed,
		false, // don't include matrix for now
	)
	if err != nil {
		log.Fatalf("Failed to simulate YRI haplotypes: %v", err)
	}

	fmt.Printf("✅ YRI Haplotype Results:\n")
	fmt.Printf("   Variants: %d\n", yriResult.NumVariants)
	fmt.Printf("   LD Blocks: %d\n", len(yriResult.LDBlocks))
	fmt.Printf("   Genetic Map: %.1f cM\n", yriResult.Statistics.GeneticMapLengthCM)
	fmt.Printf("   Mean r²: %.6f\n", yriResult.Statistics.MeanRSquared)
	fmt.Printf("   Diversity (π): %.6f\n", yriResult.Statistics.DiversityPi)
	fmt.Printf("   Haplotype Diversity: %.6f\n", yriResult.Statistics.HaplotypeDiversity)
	fmt.Printf("   Simulation Time: %.3fs\n", yriResult.SimulationTime)
	fmt.Printf("   Memory Used: %.1fMB\n", yriResult.MemoryUsedMB)

	// Test 2: Multiple populations haplotype simulation
	fmt.Println("\n🌍 Testing Multiple Population Haplotype Simulations...")
	testPops := []population.Population{
		population.YRI, // African
		population.CEU, // European
		population.CHB, // East Asian
		population.MXL, // Americas
	}

	results := make(map[population.Population]*population.HaplotypeSimulationResult)
	totalStart := time.Now()

	for _, pop := range testPops {
		result, err := simulator.SimulateHaplotypes(
			pop,
			"chr22",
			20,        // smaller sample for speed
			5_000_000, // 5Mb sequence
			&randomSeed,
			false,
		)
		if err != nil {
			log.Printf("Failed to simulate %s haplotypes: %v", pop, err)
			continue
		}

		results[pop] = result
		fmt.Printf("   ✅ %s: %d variants, %d blocks, π=%.6f, r²=%.3f\n",
			pop, result.NumVariants, len(result.LDBlocks),
			result.Statistics.DiversityPi, result.Statistics.MeanRSquared)
	}

	totalTime := time.Since(totalStart)
	fmt.Printf("\n🎉 All population haplotypes simulated in %.2fs\n", totalTime.Seconds())

	// Test 3: LD pattern validation
	fmt.Println("\n📊 Validating LD Patterns Across Populations...")
	validation, err := simulator.ValidateLDPatterns(results)
	if err != nil {
		log.Printf("LD validation failed: %v", err)
	}

	passedTests := 0
	totalTests := len(validation)
	for _, passed := range validation {
		if passed {
			passedTests++
		}
	}

	fmt.Printf("🎯 LD Pattern Validation: %d/%d tests passed\n", passedTests, totalTests)

	// Test 4: Cross-population LD analysis
	fmt.Println("\n🔍 Cross-Population LD Analysis...")
	analysis := simulator.AnalyzeCrossPopulationLD(results)

	if diversityMean, ok := analysis["diversity_mean"].(float64); ok {
		fmt.Printf("   Mean Population Diversity: %.6f\n", diversityMean)
	}
	if diversityCV, ok := analysis["diversity_cv"].(float64); ok {
		fmt.Printf("   Diversity Coefficient of Variation: %.3f\n", diversityCV)
	}
	if ldMean, ok := analysis["ld_mean"].(float64); ok {
		fmt.Printf("   Mean Population LD: %.6f\n", ldMean)
	}
	if ldCV, ok := analysis["ld_cv"].(float64); ok {
		fmt.Printf("   LD Coefficient of Variation: %.3f\n", ldCV)
	}

	// Test 5: Recombination map analysis
	fmt.Println("\n🧭 Recombination Map Analysis...")
	for pop, result := range results {
		recombMap := simulator.GetRecombinationMap(result)
		ldStats := simulator.GetLDStatistics(result)

		fmt.Printf("   %s Recombination:\n", pop)
		fmt.Printf("     Map Length: %.1f cM\n", recombMap.MapLength)
		fmt.Printf("     Effective Rate: %.2f cM/Mb\n", ldStats.EffectiveRecombRate)
		fmt.Printf("     Intervals: %d\n", len(recombMap.Positions))
	}

	// Test 6: Haplotype block structure analysis
	fmt.Println("\n🧱 Haplotype Block Structure Analysis...")
	for pop, result := range results {
		if len(result.LDBlocks) > 0 {
			totalBlockLength := 0
			totalVariants := 0
			maxBlockSize := 0
			minBlockSize := result.LDBlocks[0].Length

			for _, block := range result.LDBlocks {
				totalBlockLength += block.Length
				totalVariants += block.VariantCount
				if block.Length > maxBlockSize {
					maxBlockSize = block.Length
				}
				if block.Length < minBlockSize {
					minBlockSize = block.Length
				}
			}

			meanBlockSize := totalBlockLength / len(result.LDBlocks)
			genomeCoverage := float64(totalBlockLength) / float64(result.SequenceLength) * 100

			fmt.Printf("   %s Haplotype Blocks:\n", pop)
			fmt.Printf("     Total Blocks: %d\n", len(result.LDBlocks))
			fmt.Printf("     Mean Block Size: %d bp\n", meanBlockSize)
			fmt.Printf("     Size Range: %d - %d bp\n", minBlockSize, maxBlockSize)
			fmt.Printf("     Genome Coverage: %.1f%%\n", genomeCoverage)
			fmt.Printf("     Mean Variants per Block: %.1f\n", float64(totalVariants)/float64(len(result.LDBlocks)))
		} else {
			fmt.Printf("   %s: No haplotype blocks identified\n", pop)
		}
	}

	// Summary
	fmt.Println("\n📋 HAPLOTYPE INTEGRATION TEST SUMMARY")
	fmt.Println("=" + strings.Repeat("=", 40))
	fmt.Printf("LD pattern validation: %d/%d tests passed\n", passedTests, totalTests)
	fmt.Printf("Other steps have no pass/fail criterion; per-population failures are logged above.\n")

	// Performance summary
	fmt.Println("\n⚡ Performance Summary:")
	totalVariants := 0
	totalBlocks := 0
	totalMemory := 0.0
	totalSimTime := 0.0

	for _, result := range results {
		totalVariants += result.NumVariants
		totalBlocks += len(result.LDBlocks)
		totalMemory += result.MemoryUsedMB
		totalSimTime += result.SimulationTime
	}

	fmt.Printf("   Total Variants Simulated: %d\n", totalVariants)
	fmt.Printf("   Total LD Blocks Identified: %d\n", totalBlocks)
	fmt.Printf("   Total Memory Used: %.1f MB\n", totalMemory)
	fmt.Printf("   Total Simulation Time: %.2f seconds\n", totalSimTime)
	fmt.Printf("   Average Memory per Population: %.1f MB\n", totalMemory/float64(len(results)))
	fmt.Printf("   Average Time per Population: %.2f seconds\n", totalSimTime/float64(len(results)))
}
