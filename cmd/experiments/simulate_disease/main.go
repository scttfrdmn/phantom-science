// Command simulate_disease generates a synthetic patient cohort for one of
// the built-in disease configurations (cystic fibrosis, Huntington disease,
// BRCA hereditary cancer) with diseases.ClinVarIntegrator, which queries
// NCBI ClinVar over the network for the genes' pathogenic variants, and
// writes the cohort as JSON, VCF or PED plus a summary.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
// The cohorts are synthetic and not suitable for clinical use; the
// "validation" metrics compare the cohort against the configuration's own
// expected values; they are always computed.
package main

import (
	"context"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/diseases"
)

var (
	diseaseType          = flag.String("disease", "cystic_fibrosis", "Disease to simulate: cystic_fibrosis, huntington_disease, brca_hereditary_cancer")
	cohortSize           = flag.Int("cohort-size", 1000, "Number of synthetic patients to generate")
	affectedRatio        = flag.Float64("affected-ratio", 0.1, "Ratio of affected patients in cohort")
	population           = flag.String("population", "mixed", "Target population: EUR, AFR, EAS, AMR, SAS, AJ, mixed")
	outputDir            = flag.String("output", "./disease_simulation_output", "Output directory for generated cohorts")
	format               = flag.String("format", "json", "Output format: json, vcf, ped")
	clinvarAPI           = flag.String("clinvar-api-key", "", "NCBI API key for ClinVar access")
	verbose              = flag.Bool("verbose", false, "Enable verbose logging")
	includePhenotypes    = flag.Bool("phenotypes", true, "Include phenotype profiles in output")
	includeFamilyHistory = flag.Bool("family-history", true, "Include family history in synthetic patients")
)

func main() {
	flag.Parse()

	if *verbose {
		log.Println("🧬 Phantom Science Disease Simulation")
		log.Println("=====================================")
	}

	// Validate input parameters
	if err := validateInputs(); err != nil {
		log.Fatalf("Input validation failed: %v", err)
	}

	// Create output directory
	if err := os.MkdirAll(*outputDir, 0750); err != nil {
		log.Fatalf("Failed to create output directory: %v", err)
	}

	// Initialize ClinVar integrator
	integrator := diseases.NewClinVarIntegrator(*clinvarAPI)

	// Execute disease simulation
	if err := executeSimulation(integrator); err != nil {
		log.Fatalf("Disease simulation failed: %v", err)
	}

	log.Println("✅ Disease simulation completed successfully!")
}

func validateInputs() error {
	validDiseases := []string{"cystic_fibrosis", "huntington_disease", "brca_hereditary_cancer"}
	if !contains(validDiseases, *diseaseType) {
		return fmt.Errorf("invalid disease type: %s. Valid options: %s", *diseaseType, strings.Join(validDiseases, ", "))
	}

	if *cohortSize <= 0 {
		return fmt.Errorf("cohort size must be positive")
	}

	if *affectedRatio < 0 || *affectedRatio > 1 {
		return fmt.Errorf("affected ratio must be between 0 and 1")
	}

	validPopulations := []string{"EUR", "AFR", "EAS", "AMR", "SAS", "AJ", "mixed"}
	if !contains(validPopulations, *population) {
		return fmt.Errorf("invalid population: %s. Valid options: %s", *population, strings.Join(validPopulations, ", "))
	}

	validFormats := []string{"json", "vcf", "ped"}
	if !contains(validFormats, *format) {
		return fmt.Errorf("invalid format: %s. Valid options: %s", *format, strings.Join(validFormats, ", "))
	}

	return nil
}

func executeSimulation(integrator *diseases.ClinVarIntegrator) error {
	if *verbose {
		log.Printf("🎯 Simulating %s disease with %d patients (%.1f%% affected)",
			*diseaseType, *cohortSize, *affectedRatio*100)
	}

	// Get disease configuration
	config, err := getDiseaseConfiguration(*diseaseType)
	if err != nil {
		return fmt.Errorf("failed to get disease configuration: %w", err)
	}

	// Customize configuration based on command-line parameters
	customizeConfiguration(config)

	// Generate synthetic disease cohort
	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Minute)
	defer cancel()

	if *verbose {
		log.Println("🔬 Generating synthetic disease cohort...")
	}

	result, err := integrator.GenerateDiseaseSimulation(ctx, config)
	if err != nil {
		return fmt.Errorf("failed to generate disease simulation: %w", err)
	}

	// Save results in requested format
	if err := saveResults(result); err != nil {
		return fmt.Errorf("failed to save results: %w", err)
	}

	// Print summary
	printSimulationSummary(result)

	return nil
}

func getDiseaseConfiguration(diseaseType string) (*diseases.DiseaseSimulationConfig, error) {
	// Get base configuration from well-characterized diseases
	baseConfig, exists := diseases.WellCharacterizedDiseases[diseaseType]
	if !exists {
		return nil, fmt.Errorf("disease configuration not found for: %s", diseaseType)
	}

	// Create a copy to avoid modifying the original
	config := &diseases.DiseaseSimulationConfig{
		DiseaseName:           baseConfig.DiseaseName,
		GeneSymbols:           append([]string{}, baseConfig.GeneSymbols...),
		CarrierFrequency:      baseConfig.CarrierFrequency,
		PopulationFrequencies: make(map[string]float64),
		InheritancePattern:    baseConfig.InheritancePattern,
		AffectedRatio:         baseConfig.AffectedRatio,
	}

	// Copy population frequencies
	for pop, freq := range baseConfig.PopulationFrequencies {
		config.PopulationFrequencies[pop] = freq
	}

	// Copy penetrance model if present
	if baseConfig.PenetranceModel != nil {
		config.PenetranceModel = &diseases.PenetranceModel{
			AgeGroups:        append([]diseases.AgeGroup{}, baseConfig.PenetranceModel.AgeGroups...),
			Gender:           baseConfig.PenetranceModel.Gender,
			ModifyingFactors: make(map[string]float64),
		}

		for factor, value := range baseConfig.PenetranceModel.ModifyingFactors {
			config.PenetranceModel.ModifyingFactors[factor] = value
		}
	}

	return config, nil
}

func customizeConfiguration(config *diseases.DiseaseSimulationConfig) {
	// Set cohort size and affected ratio
	config.SyntheticCohortSize = *cohortSize
	config.AffectedRatio = *affectedRatio

	// Adjust population frequencies if specific population requested
	if *population != "mixed" {
		// Focus on single population
		newFreqs := make(map[string]float64)
		if freq, exists := config.PopulationFrequencies[*population]; exists {
			newFreqs[*population] = 1.0
			config.CarrierFrequency = freq
		} else {
			// Use default frequency for unknown populations
			newFreqs[*population] = 1.0
			config.CarrierFrequency = 0.01 // Conservative estimate
		}
		config.PopulationFrequencies = newFreqs
	}
}

func saveResults(result *diseases.SyntheticDiseaseResult) error {
	// Save detailed JSON report
	jsonPath := filepath.Join(*outputDir, fmt.Sprintf("%s_cohort_detailed.json",
		strings.ToLower(strings.ReplaceAll(result.DiseaseConfig.DiseaseName, " ", "_"))))

	if err := saveJSONReport(result, jsonPath); err != nil {
		return fmt.Errorf("failed to save JSON report: %w", err)
	}

	// Save summary report
	summaryPath := filepath.Join(*outputDir, fmt.Sprintf("%s_cohort_summary.txt",
		strings.ToLower(strings.ReplaceAll(result.DiseaseConfig.DiseaseName, " ", "_"))))

	if err := saveSummaryReport(result, summaryPath); err != nil {
		return fmt.Errorf("failed to save summary report: %w", err)
	}

	// Save in additional formats if requested
	switch *format {
	case "vcf":
		vcfPath := filepath.Join(*outputDir, fmt.Sprintf("%s_cohort.vcf",
			strings.ToLower(strings.ReplaceAll(result.DiseaseConfig.DiseaseName, " ", "_"))))
		if err := saveVCFFormat(result, vcfPath); err != nil {
			return fmt.Errorf("failed to save VCF format: %w", err)
		}

	case "ped":
		pedPath := filepath.Join(*outputDir, fmt.Sprintf("%s_cohort.ped",
			strings.ToLower(strings.ReplaceAll(result.DiseaseConfig.DiseaseName, " ", "_"))))
		if err := savePEDFormat(result, pedPath); err != nil {
			return fmt.Errorf("failed to save PED format: %w", err)
		}
	}

	return nil
}

func saveJSONReport(result *diseases.SyntheticDiseaseResult, filePath string) error {
	file, err := os.Create(filePath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	encoder := json.NewEncoder(file)
	encoder.SetIndent("", "  ")
	return encoder.Encode(result)
}

func saveSummaryReport(result *diseases.SyntheticDiseaseResult, filePath string) error {
	file, err := os.Create(filePath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	fmt.Fprintf(file, "# %s Synthetic Disease Cohort Summary\n\n", result.DiseaseConfig.DiseaseName)
	fmt.Fprintf(file, "## Generation Details\n")
	fmt.Fprintf(file, "- Generation Time: %v\n", result.GenerationTime)
	fmt.Fprintf(file, "- Inheritance Pattern: %s\n", result.DiseaseConfig.InheritancePattern)
	fmt.Fprintf(file, "- Target Genes: %s\n", strings.Join(result.DiseaseConfig.GeneSymbols, ", "))
	fmt.Fprintf(file, "\n## Cohort Statistics\n")
	fmt.Fprintf(file, "- Total Samples: %d\n", result.Statistics.TotalSamples)
	fmt.Fprintf(file, "- Affected: %d (%.1f%%)\n", result.Statistics.AffectedCount,
		float64(result.Statistics.AffectedCount)/float64(result.Statistics.TotalSamples)*100)
	fmt.Fprintf(file, "- Carriers: %d (%.1f%%)\n", result.Statistics.CarrierCount,
		float64(result.Statistics.CarrierCount)/float64(result.Statistics.TotalSamples)*100)
	fmt.Fprintf(file, "- Unaffected: %d (%.1f%%)\n", result.Statistics.UnaffectedCount,
		float64(result.Statistics.UnaffectedCount)/float64(result.Statistics.TotalSamples)*100)

	fmt.Fprintf(file, "\n## Population Breakdown\n")
	for pop, count := range result.Statistics.PopulationBreakdown {
		fmt.Fprintf(file, "- %s: %d (%.1f%%)\n", pop, count,
			float64(count)/float64(result.Statistics.TotalSamples)*100)
	}

	fmt.Fprintf(file, "\n## Validation Metrics\n")
	fmt.Fprintf(file, "- Overall Quality Score: %.3f\n", result.ValidationMetrics.OverallQualityScore)
	fmt.Fprintf(file, "- Carrier Frequency Accuracy: %.3f\n", result.ValidationMetrics.CarrierFrequencyAccuracy)

	for pop, accuracy := range result.ValidationMetrics.PopulationFrequencyAccuracy {
		fmt.Fprintf(file, "- %s Population Frequency Accuracy: %.3f\n", pop, accuracy)
	}

	return nil
}

func saveVCFFormat(result *diseases.SyntheticDiseaseResult, filePath string) error {
	file, err := os.Create(filePath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	// Write VCF header
	fmt.Fprintf(file, "##fileformat=VCFv4.2\n")
	fmt.Fprintf(file, "##source=PhantomScience_DiseaseSimulation\n")
	fmt.Fprintf(file, "##reference=GRCh38\n")
	fmt.Fprintf(file, "##INFO=<ID=DISEASE,Number=1,Type=String,Description=\"Associated disease\">\n")
	fmt.Fprintf(file, "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene symbol\">\n")
	fmt.Fprintf(file, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	fmt.Fprintf(file, "##FORMAT=<ID=CLNSIG,Number=1,Type=String,Description=\"Clinical significance\">\n")

	// Write sample header
	fmt.Fprintf(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
	for _, sample := range result.GeneratedSamples {
		fmt.Fprintf(file, "\t%s", sample.PatientID)
	}
	fmt.Fprintf(file, "\n")

	// Write variant records (simplified)
	for _, variant := range result.DiseaseConfig.PathogenicVariants {
		fmt.Fprintf(file, "%s\t%d\t%s\t%s\t%s\t999\tPASS\tDISEASE=%s;GENE=%s\tGT:CLNSIG",
			variant.Chromosome, variant.Position, variant.VariationID,
			variant.ReferenceAllele, variant.AlternateAllele,
			result.DiseaseConfig.DiseaseName, variant.GeneSymbol)

		// Write genotypes for each sample
		for _, sample := range result.GeneratedSamples {
			genotype := "0/0" // Default unaffected
			if patientVariant, exists := sample.Genotype[variant.GeneSymbol]; exists {
				if patientVariant.Zygosity == "heterozygous" {
					genotype = "0/1"
				} else if patientVariant.Zygosity == "homozygous" {
					genotype = "1/1"
				}
			}
			fmt.Fprintf(file, "\t%s:%s", genotype, variant.ClinicalSignificance)
		}
		fmt.Fprintf(file, "\n")
	}

	return nil
}

func savePEDFormat(result *diseases.SyntheticDiseaseResult, filePath string) error {
	file, err := os.Create(filePath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	// Write PED header comment
	fmt.Fprintf(file, "# Family_ID Individual_ID Paternal_ID Maternal_ID Sex Phenotype\n")

	// Write sample records
	for i, sample := range result.GeneratedSamples {
		familyID := fmt.Sprintf("FAM_%04d", i+1)
		sex := "1" // Male
		if sample.Gender == "female" {
			sex = "2"
		}

		phenotype := "1" // Unaffected
		if sample.AffectionStatus == "affected" {
			phenotype = "2" // Affected
		}

		fmt.Fprintf(file, "%s\t%s\t0\t0\t%s\t%s\n",
			familyID, sample.PatientID, sex, phenotype)
	}

	return nil
}

func printSimulationSummary(result *diseases.SyntheticDiseaseResult) {
	fmt.Println("\n🎯 Disease Simulation Summary:")
	fmt.Println("=============================")
	fmt.Printf("Disease: %s\n", result.DiseaseConfig.DiseaseName)
	fmt.Printf("Generation Time: %v\n", result.GenerationTime)
	fmt.Printf("Total Samples: %d\n", result.Statistics.TotalSamples)
	fmt.Printf("Affected: %d (%.1f%%)\n", result.Statistics.AffectedCount,
		float64(result.Statistics.AffectedCount)/float64(result.Statistics.TotalSamples)*100)
	fmt.Printf("Carriers: %d (%.1f%%)\n", result.Statistics.CarrierCount,
		float64(result.Statistics.CarrierCount)/float64(result.Statistics.TotalSamples)*100)
	fmt.Printf("Overall Quality Score: %.3f\n", result.ValidationMetrics.OverallQualityScore)

	fmt.Println("\n📊 Population Distribution:")
	for pop, count := range result.Statistics.PopulationBreakdown {
		fmt.Printf("- %s: %d (%.1f%%)\n", pop, count,
			float64(count)/float64(result.Statistics.TotalSamples)*100)
	}

	fmt.Printf("\n📁 Results saved to: %s\n", *outputDir)
}

func contains(slice []string, item string) bool {
	for _, s := range slice {
		if s == item {
			return true
		}
	}
	return false
}
