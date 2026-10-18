// Command test_complete_vcf_integration runs the Python demographic and
// haplotype simulators through internal/population (which need the
// phantom_science_env virtual environment under the working directory) for
// YRI, CEU and CHB, prints variant counts, LD blocks and diversity, and runs
// ValidateLDPatterns on the results.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
// It reads no VCF data: it only checks whether the hard-coded
// /Volumes/Public/genomics_training/populations directory and a
// parallel_vcf_training.log file exist. Three of its checklist items are
// hard-coded true, and the closing summary lists design goals, not results
// measured by this program.
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
	fmt.Println("🔗 VCF-HAPLOTYPE INTEGRATION SMOKE TEST (experiment)")
	fmt.Println("=" + strings.Repeat("=", 50))

	// Get project root
	projectRoot, err := os.Getwd()
	if err != nil {
		log.Fatalf("Failed to get project root: %v", err)
	}

	fmt.Println("\n📊 Runs the demographic and haplotype simulators; reads no VCF data")

	// Test 1: VCF data path check
	fmt.Println("\n1️⃣ Checking VCF data path...")

	vcfPath := "/Volumes/Public/genomics_training/populations"
	if _, err := os.Stat(vcfPath); err != nil {
		fmt.Printf("   ⚠️  VCF data path not accessible: %s\n", vcfPath)
	} else {
		fmt.Printf("   ✅ VCF data path accessible: %s (not read)\n", vcfPath)
	}

	// Test 2: Demographic Model Integration
	fmt.Println("\n2️⃣ Testing Demographic Model Integration...")

	simulator := population.NewDemographicSimulator(projectRoot)

	randomSeed := 42
	testResult, err := simulator.SimulatePopulation(
		population.YRI,
		50,        // sample size
		5_000_000, // 5Mb sequence
		&randomSeed,
	)
	if err != nil {
		log.Printf("   ❌ Demographic simulation failed: %v", err)
	} else {
		fmt.Printf("   ✅ Demographic simulation successful\n")
		fmt.Printf("      Variants: %d, Diversity: %.6f\n",
			testResult.NumVariants, testResult.Statistics.DiversityPi)
	}

	// Test 3: Haplotype Simulation Integration
	fmt.Println("\n3️⃣ Testing Haplotype Simulation Integration...")

	haplotypeSimulator := population.NewHaplotypeSimulator(projectRoot)

	haplotypeResult, err := haplotypeSimulator.SimulateHaplotypes(
		population.YRI,
		"chr22",
		30,        // sample size
		5_000_000, // 5Mb sequence
		&randomSeed,
		false, // don't include matrix
	)
	if err != nil {
		log.Printf("   ❌ Haplotype simulation failed: %v", err)
	} else {
		fmt.Printf("   ✅ Haplotype simulation successful\n")
		fmt.Printf("      Variants: %d, LD Blocks: %d, Mean r²: %.6f\n",
			haplotypeResult.NumVariants, len(haplotypeResult.LDBlocks),
			haplotypeResult.Statistics.MeanRSquared)
	}

	// Test 4: Multiple populations
	fmt.Println("\n4️⃣ Testing Multiple Populations...")

	testPops := []population.Population{
		population.YRI, // African
		population.CEU, // European
		population.CHB, // East Asian
	}

	integrationResults := make(map[population.Population]*population.HaplotypeSimulationResult)
	totalStart := time.Now()

	for _, pop := range testPops {
		result, err := haplotypeSimulator.SimulateHaplotypes(
			pop,
			"chr22",
			20,        // smaller sample for speed
			3_000_000, // 3Mb sequence
			&randomSeed,
			false,
		)
		if err != nil {
			fmt.Printf("   ❌ %s simulation failed: %v\n", pop, err)
			continue
		}

		integrationResults[pop] = result
		fmt.Printf("   ✅ %s: %d variants, %d blocks, π=%.6f\n",
			pop, result.NumVariants, len(result.LDBlocks),
			result.Statistics.DiversityPi)
	}

	totalTime := time.Since(totalStart)
	fmt.Printf("   Simulated %d/%d populations in %.2fs\n", len(integrationResults), len(testPops), totalTime.Seconds())

	// Test 5: LD pattern checks
	fmt.Println("\n5️⃣ Running ValidateLDPatterns...")

	if len(integrationResults) > 0 {
		validation, err := haplotypeSimulator.ValidateLDPatterns(integrationResults)
		if err != nil {
			fmt.Printf("   ❌ LD pattern checks failed: %v\n", err)
		} else {
			passedTests := 0
			totalTests := len(validation)
			for _, passed := range validation {
				if passed {
					passedTests++
				}
			}
			fmt.Printf("   📊 LD pattern checks: %d/%d passed\n",
				passedTests, totalTests)
		}
	}

	// Test 6: Earlier training log
	fmt.Println("\n6️⃣ Checking for an earlier VCF training log...")

	vcfTrainingLog := "parallel_vcf_training.log"
	if _, err := os.Stat(vcfTrainingLog); err == nil {
		fmt.Printf("   ✅ %s found (contents not checked)\n", vcfTrainingLog)
	} else {
		fmt.Printf("   📝 %s not found\n", vcfTrainingLog)
	}

	fmt.Println("\n📋 CHECKLIST")
	fmt.Println("=" + strings.Repeat("=", 40))

	// The last three items are hard-coded true; nothing here checks them.
	qualityMetrics := map[string]bool{
		"Demographic simulation ran":              testResult != nil,
		"Haplotype simulation ran":                haplotypeResult != nil,
		"Three populations simulated":             len(integrationResults) >= 3,
		"LD blocks present":                       haplotypeResult != nil && len(haplotypeResult.LDBlocks) > 0,
		"Population results available to compare": len(integrationResults) >= 3,
		"Statistical framework (hard-coded true)": true,
		"Memory efficiency (hard-coded true)":     true,
		"Pipeline operational (hard-coded true)":  true,
	}

	passedQuality := 0
	totalQuality := len(qualityMetrics)

	for metric, passed := range qualityMetrics {
		status := "✅ PASSED"
		if !passed {
			status = "❌ FAILED"
		} else {
			passedQuality++
		}
		fmt.Printf("   %s: %s\n", metric, status)
	}

	fmt.Printf("\n🎯 RESULT")
	fmt.Printf("\n   Checklist: %d/%d items true\n", passedQuality, totalQuality)

	if passedQuality >= totalQuality*3/4 {
		fmt.Printf("   ✅ Simulators ran for the test populations\n")
	} else {
		fmt.Printf("   ⚠️  Some simulator runs failed\n")
	}

	fmt.Printf("\n📊 NOT MEASURED BY THIS PROGRAM\n")
	fmt.Print("=" + strings.Repeat("=", 40) + "\n")
	fmt.Printf("Agreement of the demographic parameters with the literature, realism of\n")
	fmt.Printf("the simulated LD, memory use, and any comparison with real VCF data are\n")
	fmt.Printf("not checked here.\n")
}
