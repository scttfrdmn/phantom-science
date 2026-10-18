// Command test_population_clustering_validation runs
// internal/population/pca_clustering_validator.py on its own (it is not
// given the generated data), simulates haplotypes for several populations
// with the Python haplotype simulator, and prints diversity ordering, LD
// validation counts and a clustering "readiness" tally from variant and
// block counts. Run it from the repository root.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
// No PCA clustering of the generated sequences is performed, and some items
// in its final assessment are hard-coded to true.
package main

import (
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/population"
)

func main() {
	fmt.Println("🔬 TESTING POPULATION CLUSTERING VALIDATION")
	fmt.Println("=" + strings.Repeat("=", 50))

	// Get project root
	projectRoot, err := os.Getwd()
	if err != nil {
		log.Fatalf("Failed to get project root: %v", err)
	}

	fmt.Println("\n🎯 PHASE 2 WEEK 6 VALIDATION SUMMARY")
	fmt.Println("PCA Clustering Pattern Validation with Generated Sequences")

	// Test 1: Check PCA clustering validation dependencies
	fmt.Println("\n1️⃣ Testing PCA Clustering Validation Framework...")

	clusteringScriptPath := fmt.Sprintf("%s/internal/population/pca_clustering_validator.py", projectRoot)
	if _, err := os.Stat(clusteringScriptPath); err != nil {
		fmt.Printf("   ❌ PCA clustering validation script not found: %s\n", clusteringScriptPath)
	} else {
		fmt.Printf("   ✅ PCA clustering validation script available: %s\n", clusteringScriptPath)

		// Test PCA clustering validation script
		cmd := exec.Command("python3", clusteringScriptPath) // #nosec G204
		cmd.Dir = projectRoot

		output, err := cmd.CombinedOutput()
		if err != nil {
			fmt.Printf("   ⚠️  PCA clustering validation test failed: %v\n", err)
			fmt.Printf("   📝 Note: Dependencies may need installation (sklearn, scipy)\n")
		} else {
			fmt.Printf("   ✅ PCA clustering validation test successful\n")
			// Print key results from output
			lines := strings.Split(string(output), "\n")
			for _, line := range lines {
				if strings.Contains(line, "PCA Clustering Validation:") || strings.Contains(line, "tests passed") {
					fmt.Printf("   📊 %s\n", strings.TrimSpace(line))
				}
			}
		}
	}

	// Test 2: Generate sequences for clustering validation
	fmt.Println("\n2️⃣ Testing Sequence Generation for Clustering Analysis...")

	simulator := population.NewHaplotypeSimulator(projectRoot)

	// Generate haplotype sequences for multiple populations
	randomSeed := 42
	testPopulations := []population.Population{
		population.YRI, // African
		population.CEU, // European
		population.CHB, // East Asian
		population.MXL, // Americas
	}

	haplotypeResults := make(map[population.Population]*population.HaplotypeSimulationResult)
	generationStart := time.Now()

	for _, pop := range testPopulations {
		result, err := simulator.SimulateHaplotypes(
			pop,
			"chr22",
			25,        // sample size for clustering
			3_000_000, // 3Mb sequence for analysis
			&randomSeed,
			false, // don't include matrix for clustering analysis
		)
		if err != nil {
			fmt.Printf("   ❌ %s haplotype generation failed: %v\n", pop, err)
			continue
		}

		haplotypeResults[pop] = result
		fmt.Printf("   ✅ %s: %d variants, diversity=%.6f, LD blocks=%d\n",
			pop, result.NumVariants, result.Statistics.DiversityPi, len(result.LDBlocks))
	}

	generationTime := time.Since(generationStart)
	fmt.Printf("   🎉 Sequence generation completed in %.2fs\n", generationTime.Seconds())

	// Test 3: Population Differentiation Assessment
	fmt.Println("\n3️⃣ Testing Population Differentiation Assessment...")

	if len(haplotypeResults) >= 2 {
		// Calculate basic population differentiation metrics
		diversities := make(map[population.Population]float64)
		variantCounts := make(map[population.Population]int)

		for pop, result := range haplotypeResults {
			diversities[pop] = result.Statistics.DiversityPi
			variantCounts[pop] = result.NumVariants
		}

		// Calculate diversity statistics
		var diversitySum float64
		var diversityValues []float64
		for _, diversity := range diversities {
			diversitySum += diversity
			diversityValues = append(diversityValues, diversity)
		}

		diversityMean := diversitySum / float64(len(diversities))

		// Calculate diversity range and coefficient of variation
		minDiversity, maxDiversity := diversityValues[0], diversityValues[0]
		for _, d := range diversityValues {
			if d < minDiversity {
				minDiversity = d
			}
			if d > maxDiversity {
				maxDiversity = d
			}
		}

		diversityRange := maxDiversity - minDiversity
		var diversityCV float64
		if diversityMean > 0 {
			diversityCV = diversityRange / diversityMean
		} else {
			diversityCV = 0
		}

		fmt.Printf("   📊 Population Differentiation Analysis:\n")
		fmt.Printf("      Mean Diversity: %.6f\n", diversityMean)
		fmt.Printf("      Diversity Range: %.6f\n", diversityRange)
		fmt.Printf("      Diversity CV: %.3f\n", diversityCV)

		// Test expected diversity hierarchy (YRI > CEU > CHB)
		yriDiversity := diversities[population.YRI]
		ceuDiversity, hasEuropean := diversities[population.CEU]
		chbDiversity, hasEastAsian := diversities[population.CHB]

		hierarchyTests := 0
		hierarchyPassed := 0

		if hasEuropean {
			hierarchyTests++
			if yriDiversity > ceuDiversity {
				hierarchyPassed++
				fmt.Printf("      ✅ African > European diversity: %.6f > %.6f\n", yriDiversity, ceuDiversity)
			} else {
				fmt.Printf("      ❌ African ≤ European diversity: %.6f ≤ %.6f\n", yriDiversity, ceuDiversity)
			}
		}

		if hasEastAsian {
			hierarchyTests++
			if yriDiversity > chbDiversity {
				hierarchyPassed++
				fmt.Printf("      ✅ African > East Asian diversity: %.6f > %.6f\n", yriDiversity, chbDiversity)
			} else {
				fmt.Printf("      ❌ African ≤ East Asian diversity: %.6f ≤ %.6f\n", yriDiversity, chbDiversity)
			}
		}

		if hasEuropean && hasEastAsian {
			hierarchyTests++
			if ceuDiversity > chbDiversity {
				hierarchyPassed++
				fmt.Printf("      ✅ European > East Asian diversity: %.6f > %.6f\n", ceuDiversity, chbDiversity)
			} else {
				fmt.Printf("      ❌ European ≤ East Asian diversity: %.6f ≤ %.6f\n", ceuDiversity, chbDiversity)
			}
		}

		fmt.Printf("      🎯 Diversity Hierarchy: %d/%d tests passed\n", hierarchyPassed, hierarchyTests)

		if hierarchyPassed >= hierarchyTests*2/3 {
			fmt.Printf("   ✅ Population diversity hierarchy matches expectations\n")
		} else {
			fmt.Printf("   ⚠️  Population diversity hierarchy needs refinement\n")
		}
	}

	// Test 4: Cross-Population LD Structure Analysis
	fmt.Println("\n4️⃣ Testing Cross-Population LD Structure Analysis...")

	if len(haplotypeResults) >= 2 {
		ldValidation, err := simulator.ValidateLDPatterns(haplotypeResults)
		if err != nil {
			fmt.Printf("   ❌ LD pattern validation failed: %v\n", err)
		} else {
			passedLDTests := 0
			totalLDTests := len(ldValidation)
			for _, passed := range ldValidation {
				if passed {
					passedLDTests++
				}
			}

			fmt.Printf("   📊 LD Pattern Validation: %d/%d tests passed\n", passedLDTests, totalLDTests)

			// Analyze LD statistics across populations
			fmt.Printf("   📈 Cross-Population LD Analysis:\n")
			for pop, result := range haplotypeResults {
				fmt.Printf("      %s: Mean r²=%.6f, LD blocks=%d\n",
					pop, result.Statistics.MeanRSquared, len(result.LDBlocks))
			}

			if passedLDTests >= totalLDTests*3/4 {
				fmt.Printf("   ✅ LD patterns show proper population structure\n")
			} else {
				fmt.Printf("   ⚠️  LD patterns may need enhancement for clustering\n")
			}
		}
	}

	// Test 5: Clustering Readiness Assessment
	fmt.Println("\n5️⃣ Testing Clustering Readiness Assessment...")

	clusteringReadyCount := 0
	insufficientVariants := 0

	for pop, result := range haplotypeResults {
		// Assess readiness for PCA clustering analysis
		sufficientVariants := result.NumVariants > 200                 // Need enough variants for PCA
		sufficientLDStructure := len(result.LDBlocks) > 5              // Need LD structure
		sufficientDiversity := result.Statistics.DiversityPi > 0.00001 // Need genetic diversity

		readinessScore := 0
		if sufficientVariants {
			readinessScore++
		} else {
			insufficientVariants++
		}
		if sufficientLDStructure {
			readinessScore++
		}
		if sufficientDiversity {
			readinessScore++
		}

		if readinessScore >= 2 {
			clusteringReadyCount++
			fmt.Printf("   ✅ %s ready for clustering: %d variants, %d blocks, π=%.6f\n",
				pop, result.NumVariants, len(result.LDBlocks), result.Statistics.DiversityPi)
		} else {
			fmt.Printf("   ⚠️  %s may need optimization: %d variants, %d blocks, π=%.6f\n",
				pop, result.NumVariants, len(result.LDBlocks), result.Statistics.DiversityPi)
		}
	}

	fmt.Printf("   📊 Clustering Readiness: %d/%d populations ready\n",
		clusteringReadyCount, len(haplotypeResults))

	if clusteringReadyCount >= len(haplotypeResults)*3/4 {
		fmt.Printf("   ✅ Sufficient populations ready for PCA clustering validation\n")
	} else {
		fmt.Printf("   ⚠️  More populations need optimization for robust clustering analysis\n")
	}

	// Test 6: Continental Grouping Preparation
	fmt.Println("\n6️⃣ Testing Continental Grouping Preparation...")

	continentalGroups := map[string][]population.Population{
		"African":    {population.YRI},
		"European":   {population.CEU},
		"East_Asian": {population.CHB},
		"Americas":   {population.MXL},
	}

	continentsRepresented := 0
	for continent, pops := range continentalGroups {
		continentReady := false
		for _, pop := range pops {
			if _, exists := haplotypeResults[pop]; exists {
				continentReady = true
				break
			}
		}

		if continentReady {
			continentsRepresented++
			fmt.Printf("   ✅ %s populations represented for clustering\n", continent)
		} else {
			fmt.Printf("   ❌ %s populations missing for clustering\n", continent)
		}
	}

	fmt.Printf("   🌍 Continental Representation: %d/%d continents ready\n",
		continentsRepresented, len(continentalGroups))

	if continentsRepresented >= 3 {
		fmt.Printf("   ✅ Sufficient continental diversity for PCA clustering validation\n")
	} else {
		fmt.Printf("   ⚠️  Need more continental groups for comprehensive clustering analysis\n")
	}

	// Validate integration quality
	fmt.Println("\n📋 POPULATION CLUSTERING VALIDATION ASSESSMENT")
	fmt.Println("=" + strings.Repeat("=", 50))

	// Quality metrics
	qualityMetrics := map[string]bool{
		"PCA clustering framework (not checked)":                        true,
		"Sequence generation (>= 3 populations simulated)":              len(haplotypeResults) >= 3,
		"Population differentiation (only checks >= 3 populations)":     len(haplotypeResults) >= 3,
		"LD structure (>= 2 populations clustering-ready)":              clusteringReadyCount >= 2,
		"Continental groups represented":                                continentsRepresented >= 3,
		"Clustering readiness (>= 2/3 of populations)":                  clusteringReadyCount >= len(haplotypeResults)*2/3,
		"Population diversity hierarchy (not checked here; see step 3)": true,
		"Phase 2 Week 6 foundation (not checked)":                       true,
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

	// Final assessment
	fmt.Printf("\n🎯 PHASE 2 WEEK 6 COMPLETION ASSESSMENT\n")
	fmt.Printf("   Clustering Validation Quality: %d/%d checks true (3 of them hard-coded)\n", passedQuality, totalQuality)

	if passedQuality >= totalQuality*3/4 {
		fmt.Printf("   Most checks true; no PCA clustering of generated sequences was performed\n")
	} else {
		fmt.Printf("   ⚠️  Several checks failed; see above\n")
	}
}
