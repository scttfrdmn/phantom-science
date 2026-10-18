// Command test_demographic_hmm_integration simulates YRI, CEU and CHB
// variants with the Python demographic simulator (run from the current
// directory, which must be the repository root with phantom_science_env),
// feeds them into an EnhancedHMMKmerModel, generates a 1 kb sequence per
// population and prints its GC content against a fixed per-population
// target.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"fmt"
	"log"
	"os"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

func main() {
	fmt.Println("🧬 TESTING DEMOGRAPHIC HMM INTEGRATION")
	fmt.Println("=" + strings.Repeat("=", 41))

	// Get project root
	projectRoot, err := os.Getwd()
	if err != nil {
		log.Fatalf("Failed to get project root: %v", err)
	}

	// Create demographic simulator
	simulator := population.NewDemographicSimulator(projectRoot)

	// Test populations
	testPops := []population.Population{
		population.YRI, // African
		population.CEU, // European
		population.CHB, // East Asian
	}

	fmt.Println("\n📍 Testing Demographic HMM Integration for 3 Populations...")

	for _, pop := range testPops {
		fmt.Printf("\n🧬 Processing %s Population...\n", pop)

		// Simulate population with demographic model
		randomSeed := 42
		result, err := simulator.SimulatePopulation(
			pop,
			50,        // sample size
			5_000_000, // 5Mb sequence for speed
			&randomSeed,
		)
		if err != nil {
			log.Printf("Failed to simulate %s: %v", pop, err)
			continue
		}

		fmt.Printf("   📊 Simulated: %d variants, π=%.6f, D=%.3f\n",
			result.NumVariants, result.Statistics.DiversityPi, result.Statistics.TajimaD)

		// Create enhanced HMM model
		enhancedModel := authenticity.NewEnhancedHMMKmerModel(8, string(pop))

		// Convert demographic variants to HMM format
		hmmVariants := convertToHMMVariants(result.Variants)

		// Integrate demographic data
		err = enhancedModel.IntegrateDemographicVariants(
			"chr1",
			result.SequenceLength,
			hmmVariants,
			convertStatistics(result.Statistics),
		)
		if err != nil {
			log.Printf("Failed to integrate %s variants: %v", pop, err)
			continue
		}

		// Test sequence generation
		fmt.Printf("   🔬 Testing population-specific sequence generation...\n")
		generatedSeq, err := enhancedModel.GeneratePopulationSequence(1000)
		if err != nil {
			log.Printf("Failed to generate %s sequence: %v", pop, err)
			continue
		}

		// Analyze generated sequence
		gcContent := calculateGCContent(generatedSeq.Sequence)
		fmt.Printf("   📈 Generated sequence: length=%d, GC=%.3f\n",
			len(generatedSeq.Sequence), gcContent)

		// Validate model statistics
		modelStats := enhancedModel.GetEnhancedModelStats()
		fmt.Printf("   ✅ Model enhanced: %v, k-mers: %d\n",
			modelStats["demographic_enhanced"], enhancedModel.TotalKmers())

		// Validate population-specific characteristics
		expectedGC := getExpectedGCContent(string(pop))
		gcDiff := abs(gcContent - expectedGC)
		if gcDiff < 0.05 {
			fmt.Printf("   ✅ GC content within expected range (target=%.3f, actual=%.3f)\n",
				expectedGC, gcContent)
		} else {
			fmt.Printf("   ⚠️  GC content outside expected range (target=%.3f, actual=%.3f)\n",
				expectedGC, gcContent)
		}
	}

	// Integration validation
	fmt.Println("\n🔍 Validating Integration Quality...")

	// Test 1: Cross-population comparison
	fmt.Println("   📊 Cross-population sequence comparison...")
	yriModel := authenticity.NewEnhancedHMMKmerModel(8, "YRI")
	ceuModel := authenticity.NewEnhancedHMMKmerModel(8, "CEU")

	// Generate test sequences
	yriSeq, _ := generateTestSequence(yriModel, 1000)
	ceuSeq, _ := generateTestSequence(ceuModel, 1000)

	if yriSeq != nil && ceuSeq != nil {
		yriGC := calculateGCContent(yriSeq)
		ceuGC := calculateGCContent(ceuSeq)

		fmt.Printf("   📈 YRI GC: %.3f, CEU GC: %.3f, diff: %.3f\n",
			yriGC, ceuGC, abs(yriGC-ceuGC))

		if abs(yriGC-ceuGC) > 0.005 {
			fmt.Printf("   ✅ Population-specific differences detected\n")
		} else {
			fmt.Printf("   ⚠️  Population differences may be too subtle\n")
		}
	}

	// Test 2: Demographic vs VCF comparison
	fmt.Println("   🔄 Comparing demographic vs VCF-trained models...")

	// Load existing VCF-trained model for comparison
	// vcfModel := authenticity.NewHMMKmerModel(8)
	// Note: In practice, would load actual VCF-trained model

	fmt.Printf("   📋 Not implemented: no VCF-trained model is loaded, so no comparison is made\n")

	// Summary
	fmt.Println("\n📋 DEMOGRAPHIC HMM INTEGRATION SUMMARY")
	fmt.Println("=" + strings.Repeat("=", 40))
	fmt.Printf("Run finished. See the per-population lines above for GC checks;\n")
	fmt.Printf("this summary does not itself validate anything.\n")
}

// Helper functions

func convertToHMMVariants(variants []population.VariantData) []authenticity.DemographicVariant {
	hmmVariants := make([]authenticity.DemographicVariant, len(variants))

	for i, variant := range variants {
		var refAllele, altAllele byte = 'A', 'T'

		if len(variant.Alleles) >= 2 {
			if len(variant.Alleles[0]) > 0 {
				refAllele = variant.Alleles[0][0]
			}
			if len(variant.Alleles[1]) > 0 {
				altAllele = variant.Alleles[1][0]
			}
		}

		hmmVariants[i] = authenticity.DemographicVariant{
			Position:        variant.Position,
			ReferenceAllele: refAllele,
			AlternateAllele: altAllele,
			Frequency:       variant.Frequency,
			Genotypes:       variant.Genotypes,
		}
	}

	return hmmVariants
}

func convertStatistics(stats population.PopulationStatistics) map[string]float64 {
	return map[string]float64{
		"diversity_pi":    stats.DiversityPi,
		"watterson_theta": stats.WattersonTheta,
		"tajima_d":        stats.TajimaD,
		"mean_frequency":  stats.MeanAlleleFrequency,
		"num_variants":    float64(stats.NumVariants),
		"rare_variants":   float64(stats.RareVariantsCount),
	}
}

func calculateGCContent(sequence []byte) float64 {
	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}
	return float64(gcCount) / float64(len(sequence))
}

func getExpectedGCContent(population string) float64 {
	expectedGC := map[string]float64{
		"YRI": 0.41,
		"CEU": 0.42,
		"CHB": 0.40,
	}

	if gc, exists := expectedGC[population]; exists {
		return gc
	}
	return 0.41
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func generateTestSequence(model *authenticity.EnhancedHMMKmerModel, length int) ([]byte, error) {
	if !model.IsEnhanced() {
		return nil, fmt.Errorf("model not enhanced")
	}

	result, err := model.GeneratePopulationSequence(length)
	if err != nil {
		return nil, err
	}

	return result.Sequence, nil
}
