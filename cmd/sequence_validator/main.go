package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

// ValidationResult represents validation results for a single sequence
type ValidationResult struct {
	Population     string  `json:"population"`
	Region         string  `json:"region"`
	Length         int     `json:"length"`
	SequenceNumber int     `json:"sequence_number"`
	Filename       string  `json:"filename"`
	
	// Sequence Quality Metrics
	ActualLength    int     `json:"actual_length"`
	GCContent       float64 `json:"gc_content"`
	Complexity      float64 `json:"complexity"`
	ValidBases      int     `json:"valid_bases"`
	InvalidBases    int     `json:"invalid_bases"`
	
	// Genomic Validity
	MaxHomopolymer  int     `json:"max_homopolymer"`
	HasStopCodons   bool    `json:"has_stop_codons"`
	
	// Population-Specific Metrics
	TargetGC        float64 `json:"target_gc"`
	GCDeviation     float64 `json:"gc_deviation"`
	ComplexityScore string  `json:"complexity_score"`
	
	// Overall Assessment
	IsValid         bool    `json:"is_valid"`
	ValidationScore float64 `json:"validation_score"`
	Issues          []string `json:"issues"`
}

// PopulationSummary represents summary statistics for a population
type PopulationSummary struct {
	Population      string  `json:"population"`
	Region          string  `json:"region"`
	TotalSequences  int     `json:"total_sequences"`
	ValidSequences  int     `json:"valid_sequences"`
	AvgGCContent    float64 `json:"avg_gc_content"`
	AvgComplexity   float64 `json:"avg_complexity"`
	AvgScore        float64 `json:"avg_score"`
	Issues          map[string]int `json:"issues"`
}

// ValidationReport represents the complete validation report
type ValidationReport struct {
	TotalSequences     int                 `json:"total_sequences"`
	ValidSequences     int                 `json:"valid_sequences"`
	OverallScore       float64             `json:"overall_score"`
	PopulationResults  []ValidationResult  `json:"population_results"`
	PopulationSummary  []PopulationSummary `json:"population_summary"`
	RegionalAnalysis   map[string]float64  `json:"regional_analysis"`
	ValidationCriteria string              `json:"validation_criteria"`
	Timestamp          string              `json:"timestamp"`
}

func main() {
	log.SetFlags(log.LstdFlags)
	
	fmt.Println("=== GENOMIC SEQUENCE VALIDATION ===")
	fmt.Println("Validating generated sequences for biological authenticity")
	fmt.Println()
	
	inputDir := "sequences/validation"
	outputFile := "validation_report.json"
	
	// Check if input directory exists
	if _, err := os.Stat(inputDir); os.IsNotExist(err) {
		log.Fatalf("Input directory does not exist: %s", inputDir)
	}
	
	fmt.Printf("📁 Input directory: %s\n", inputDir)
	fmt.Printf("📄 Output report: %s\n", outputFile)
	fmt.Println()
	
	// Validate all sequences
	results, err := validateAllSequences(inputDir)
	if err != nil {
		log.Fatalf("Validation failed: %v", err)
	}
	
	// Generate summary statistics
	report := generateValidationReport(results)
	
	// Display results
	displayValidationResults(report)
	
	// Save JSON report
	err = saveValidationReport(report, outputFile)
	if err != nil {
		log.Printf("Warning: Failed to save validation report: %v", err)
	} else {
		fmt.Printf("\n✅ Validation report saved to: %s\n", outputFile)
	}
	
	fmt.Println("\n🔬 Sequence validation completed")
}

func validateAllSequences(inputDir string) ([]ValidationResult, error) {
	var results []ValidationResult
	
	// Walk through all population directories
	err := filepath.Walk(inputDir, func(path string, info os.FileInfo, err error) error {
		if err != nil {
			return err
		}
		
		// Skip directories and non-FASTA files
		if info.IsDir() || !strings.HasSuffix(strings.ToLower(info.Name()), ".fasta") {
			return nil
		}
		
		// Parse filename to extract metadata
		population, length, seqNum, err := parseFilename(info.Name())
		if err != nil {
			log.Printf("Warning: Could not parse filename %s: %v", info.Name(), err)
			return nil
		}
		
		// Validate sequence
		result, err := validateSequenceFile(path, population, length, seqNum)
		if err != nil {
			log.Printf("Warning: Failed to validate %s: %v", path, err)
			return nil
		}
		
		results = append(results, *result)
		return nil
	})
	
	return results, err
}

func parseFilename(filename string) (population string, length int, seqNum int, err error) {
	// Expected format: POP_L1000_S01.fasta
	name := strings.TrimSuffix(filename, ".fasta")
	parts := strings.Split(name, "_")
	
	if len(parts) != 3 {
		return "", 0, 0, fmt.Errorf("unexpected filename format: %s", filename)
	}
	
	population = parts[0]
	
	// Parse length (L1000 -> 1000)
	lengthStr := strings.TrimPrefix(parts[1], "L")
	length, err = strconv.Atoi(lengthStr)
	if err != nil {
		return "", 0, 0, fmt.Errorf("invalid length in filename: %s", parts[1])
	}
	
	// Parse sequence number (S01 -> 1)
	seqNumStr := strings.TrimPrefix(parts[2], "S")
	seqNum, err = strconv.Atoi(seqNumStr)
	if err != nil {
		return "", 0, 0, fmt.Errorf("invalid sequence number in filename: %s", parts[2])
	}
	
	return population, length, seqNum, nil
}

func validateSequenceFile(filepath, population string, expectedLength, seqNum int) (*ValidationResult, error) {
	file, err := os.Open(filepath) // #nosec G304
	if err != nil {
		return nil, err
	}
	defer file.Close()
	
	scanner := bufio.NewScanner(file)
	
	// Read header line
	if !scanner.Scan() {
		return nil, fmt.Errorf("empty file or missing header")
	}
	
	header := scanner.Text()
	if !strings.HasPrefix(header, ">") {
		return nil, fmt.Errorf("invalid FASTA format: missing header")
	}
	
	// Parse header metadata
	targetGC, region := parseHeaderMetadata(header)
	
	// Read sequence
	var sequence strings.Builder
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" {
			sequence.WriteString(line)
		}
	}
	
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	
	seqData := sequence.String()
	if len(seqData) == 0 {
		return nil, fmt.Errorf("empty sequence")
	}
	
	// Validate sequence
	return validateSequence(seqData, population, region, expectedLength, seqNum, targetGC, filepath)
}

func parseHeaderMetadata(header string) (targetGC float64, region string) {
	// Extract target_gc and region from header
	// Format: >POP_L1000_S01 population=POP region=Region ... model_target_gc=0.496
	parts := strings.Fields(header)
	
	for _, part := range parts {
		if strings.HasPrefix(part, "region=") {
			region = strings.TrimPrefix(part, "region=")
		} else if strings.HasPrefix(part, "model_target_gc=") {
			gcStr := strings.TrimPrefix(part, "model_target_gc=")
			if gc, err := strconv.ParseFloat(gcStr, 64); err == nil {
				targetGC = gc
			}
		}
	}
	
	return targetGC, region
}

func validateSequence(sequence, population, region string, expectedLength, seqNum int, targetGC float64, filename string) (*ValidationResult, error) {
	result := &ValidationResult{
		Population:     population,
		Region:         region,
		Length:         expectedLength,
		SequenceNumber: seqNum,
		Filename:       filepath.Base(filename),
		TargetGC:       targetGC,
		Issues:         []string{},
	}
	
	// Basic sequence validation
	result.ActualLength = len(sequence)
	result.ValidBases, result.InvalidBases = countValidBases(sequence)
	result.GCContent = calculateGCContent(sequence)
	result.Complexity = calculateComplexity(sequence, 8)
	result.MaxHomopolymer = findMaxHomopolymer(sequence)
	result.HasStopCodons = hasStopCodons(sequence)
	
	// Calculate deviations
	result.GCDeviation = result.GCContent - targetGC
	
	// Assign complexity score
	result.ComplexityScore = classifyComplexity(result.Complexity)
	
	// Check for issues
	if result.ActualLength != expectedLength {
		result.Issues = append(result.Issues, fmt.Sprintf("Length mismatch: expected %d, got %d", expectedLength, result.ActualLength))
	}
	
	if result.InvalidBases > 0 {
		result.Issues = append(result.Issues, fmt.Sprintf("Contains %d invalid bases", result.InvalidBases))
	}
	
	if abs(result.GCDeviation) > 0.05 { // 5% tolerance
		result.Issues = append(result.Issues, fmt.Sprintf("GC content deviation: %.3f (target: %.3f)", result.GCContent, targetGC))
	}
	
	if result.Complexity < 0.3 {
		result.Issues = append(result.Issues, fmt.Sprintf("Low complexity: %.4f", result.Complexity))
	}
	
	if result.MaxHomopolymer > 10 {
		result.Issues = append(result.Issues, fmt.Sprintf("Long homopolymer run: %d bases", result.MaxHomopolymer))
	}
	
	// Calculate overall validation score
	result.ValidationScore = calculateValidationScore(result)
	result.IsValid = result.ValidationScore >= 0.7 && len(result.Issues) == 0
	
	return result, nil
}

func countValidBases(sequence string) (valid, invalid int) {
	for _, base := range sequence {
		switch base {
		case 'A', 'T', 'G', 'C', 'a', 't', 'g', 'c':
			valid++
		default:
			invalid++
		}
	}
	return valid, invalid
}

func calculateGCContent(sequence string) float64 {
	if len(sequence) == 0 {
		return 0.0
	}
	
	gcCount := 0
	for _, base := range strings.ToUpper(sequence) {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}
	
	return float64(gcCount) / float64(len(sequence))
}

func calculateComplexity(sequence string, k int) float64 {
	if len(sequence) < k {
		return 0.0
	}
	
	uniqueKmers := make(map[string]bool)
	sequence = strings.ToUpper(sequence)
	
	for i := 0; i <= len(sequence)-k; i++ {
		kmer := sequence[i : i+k]
		uniqueKmers[kmer] = true
	}
	
	totalKmers := len(sequence) - k + 1
	return float64(len(uniqueKmers)) / float64(totalKmers)
}

func findMaxHomopolymer(sequence string) int {
	if len(sequence) == 0 {
		return 0
	}
	
	maxRun := 1
	currentRun := 1
	sequence = strings.ToUpper(sequence)
	
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

func hasStopCodons(sequence string) bool {
	stopCodons := []string{"TAA", "TAG", "TGA"}
	sequence = strings.ToUpper(sequence)
	
	for _, codon := range stopCodons {
		if strings.Contains(sequence, codon) {
			return true
		}
	}
	
	return false
}

func classifyComplexity(complexity float64) string {
	switch {
	case complexity >= 0.8:
		return "Very High"
	case complexity >= 0.6:
		return "High"
	case complexity >= 0.4:
		return "Medium"
	case complexity >= 0.2:
		return "Low"
	default:
		return "Very Low"
	}
}

func calculateValidationScore(result *ValidationResult) float64 {
	score := 1.0
	
	// Penalize for issues
	score -= float64(len(result.Issues)) * 0.1
	
	// Reward high complexity
	if result.Complexity >= 0.6 {
		score += 0.2
	} else if result.Complexity < 0.3 {
		score -= 0.3
	}
	
	// Penalize GC deviation
	gcPenalty := abs(result.GCDeviation) * 2.0
	score -= gcPenalty
	
	// Penalize invalid bases
	if result.InvalidBases > 0 {
		score -= 0.5
	}
	
	// Clamp between 0 and 1
	if score < 0 {
		score = 0
	}
	if score > 1 {
		score = 1
	}
	
	return score
}

func generateValidationReport(results []ValidationResult) *ValidationReport {
	report := &ValidationReport{
		TotalSequences:     len(results),
		PopulationResults:  results,
		RegionalAnalysis:   make(map[string]float64),
		ValidationCriteria: "Biological authenticity based on GC content, k-mer complexity, homopolymer analysis",
		Timestamp:          fmt.Sprintf("%v", os.Getenv("NOW")),
	}
	
	// Calculate overall statistics
	validCount := 0
	totalScore := 0.0
	
	for _, result := range results {
		if result.IsValid {
			validCount++
		}
		totalScore += result.ValidationScore
	}
	
	report.ValidSequences = validCount
	if len(results) > 0 {
		report.OverallScore = totalScore / float64(len(results))
	}
	
	// Generate population summaries
	popStats := make(map[string]*PopulationSummary)
	
	for _, result := range results {
		if _, exists := popStats[result.Population]; !exists {
			popStats[result.Population] = &PopulationSummary{
				Population: result.Population,
				Region:     result.Region,
				Issues:     make(map[string]int),
			}
		}
		
		summary := popStats[result.Population]
		summary.TotalSequences++
		
		if result.IsValid {
			summary.ValidSequences++
		}
		
		summary.AvgGCContent += result.GCContent
		summary.AvgComplexity += result.Complexity
		summary.AvgScore += result.ValidationScore
		
		for _, issue := range result.Issues {
			summary.Issues[issue]++
		}
	}
	
	// Finalize population summaries
	for _, summary := range popStats {
		if summary.TotalSequences > 0 {
			summary.AvgGCContent /= float64(summary.TotalSequences)
			summary.AvgComplexity /= float64(summary.TotalSequences)
			summary.AvgScore /= float64(summary.TotalSequences)
		}
		report.PopulationSummary = append(report.PopulationSummary, *summary)
	}
	
	// Sort population summaries
	sort.Slice(report.PopulationSummary, func(i, j int) bool {
		return report.PopulationSummary[i].AvgScore > report.PopulationSummary[j].AvgScore
	})
	
	// Regional analysis
	regionStats := make(map[string][]float64)
	for _, result := range results {
		regionStats[result.Region] = append(regionStats[result.Region], result.ValidationScore)
	}
	
	for region, scores := range regionStats {
		total := 0.0
		for _, score := range scores {
			total += score
		}
		report.RegionalAnalysis[region] = total / float64(len(scores))
	}
	
	return report
}

func displayValidationResults(report *ValidationReport) {
	fmt.Println("🔬 VALIDATION RESULTS SUMMARY:")
	fmt.Println("==============================")
	fmt.Printf("✅ Valid sequences: %d/%d (%.1f%%)\n", 
		report.ValidSequences, report.TotalSequences, 
		float64(report.ValidSequences)/float64(report.TotalSequences)*100)
	fmt.Printf("📊 Overall validation score: %.3f\n", report.OverallScore)
	fmt.Println()
	
	fmt.Println("🌍 POPULATION PERFORMANCE:")
	fmt.Println("==========================")
	for i, summary := range report.PopulationSummary {
		fmt.Printf("[%2d] %s (%s):\n", i+1, summary.Population, summary.Region)
		fmt.Printf("     Valid: %d/%d (%.1f%%) | Score: %.3f\n", 
			summary.ValidSequences, summary.TotalSequences,
			float64(summary.ValidSequences)/float64(summary.TotalSequences)*100,
			summary.AvgScore)
		fmt.Printf("     GC: %.3f | Complexity: %.3f\n", 
			summary.AvgGCContent, summary.AvgComplexity)
		
		if len(summary.Issues) > 0 {
			fmt.Printf("     Issues: ")
			issueCount := 0
			for issue, count := range summary.Issues {
				if issueCount > 0 {
					fmt.Printf(", ")
				}
				fmt.Printf("%s (%d)", issue, count)
				issueCount++
				if issueCount >= 2 { // Limit display
					break
				}
			}
			fmt.Println()
		}
		fmt.Println()
	}
	
	fmt.Println("🌍 REGIONAL ANALYSIS:")
	fmt.Println("====================")
	regions := []string{"Africa", "Europe", "East Asia", "Americas"}
	for _, region := range regions {
		if score, exists := report.RegionalAnalysis[region]; exists {
			fmt.Printf("🌍 %s: %.3f average validation score\n", region, score)
		}
	}
}

func saveValidationReport(report *ValidationReport, filename string) error {
	data, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		return err
	}
	
	return os.WriteFile(filename, data, 0600)
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}