package population

import (
	"encoding/json"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"time"
)

// RecombinationMap represents a recombination map for a chromosome
type RecombinationMap struct {
	Chromosome     string    `json:"chromosome"`
	Positions      []int     `json:"positions"`
	Rates          []float64 `json:"rates"`
	MapLength      float64   `json:"map_length"`
	PhysicalLength int       `json:"physical_length"`
}

// HaplotypeBlock represents a haplotype block with LD structure
type HaplotypeBlock struct {
	StartPosition int        `json:"start_position"`
	EndPosition   int        `json:"end_position"`
	Length        int        `json:"length"`
	VariantCount  int        `json:"variant_count"`
	RSquaredMean  float64    `json:"r_squared_mean"`
	Population    Population `json:"population"`
}

// HaplotypeVariant represents a variant in haplotype simulation
type HaplotypeVariant struct {
	Position        int      `json:"position"`
	Alleles         []string `json:"alleles"`
	AncestralState  string   `json:"ancestral_state"`
	Frequency       float64  `json:"frequency"`
}

// HaplotypeStatistics holds comprehensive haplotype statistics
type HaplotypeStatistics struct {
	DiversityPi         float64 `json:"diversity_pi"`
	SegregatingSites    int     `json:"segregating_sites"`
	NumVariants         int     `json:"n_variants"`
	NumHaplotypes       int     `json:"n_haplotypes"`
	MeanRSquared        float64 `json:"mean_r_squared"`
	MedianRSquared      float64 `json:"median_r_squared"`
	HighLDFraction      float64 `json:"high_ld_fraction"`
	EffectiveRecombRate float64 `json:"effective_recomb_rate"`
	GeneticMapLengthCM  float64 `json:"genetic_map_length_cm"`
	HaplotypeDiversity  float64 `json:"haplotype_diversity"`
	LDPairsAnalyzed     int     `json:"ld_pairs_analyzed"`
}

// HaplotypeSimulationResult holds complete haplotype simulation results
type HaplotypeSimulationResult struct {
	Population        Population             `json:"population"`
	Chromosome        string                 `json:"chromosome"`
	NumVariants       int                    `json:"num_variants"`
	NumTrees          int                    `json:"num_trees"`
	NumSamples        int                    `json:"num_samples"`
	SequenceLength    int                    `json:"sequence_length"`
	SimulationTime    float64                `json:"simulation_time"`
	MemoryUsedMB      float64                `json:"memory_used_mb"`
	Variants          []HaplotypeVariant     `json:"variants"`
	RecombinationMap  RecombinationMap       `json:"recombination_map"`
	LDBlocks          []HaplotypeBlock       `json:"ld_blocks"`
	Statistics        HaplotypeStatistics    `json:"statistics"`
	HaplotypeMatrix   [][]int                `json:"haplotype_matrix,omitempty"` // Optional - can be large
}

// HaplotypeSimulator provides interface to Python haplotype simulation
type HaplotypeSimulator struct {
	pythonPath string
	scriptPath string
	envPath    string
}

// NewHaplotypeSimulator creates a new haplotype simulator
func NewHaplotypeSimulator(projectRoot string) *HaplotypeSimulator {
	return &HaplotypeSimulator{
		pythonPath: filepath.Join(projectRoot, "phantom_science_env", "bin", "python3"),
		scriptPath: filepath.Join(projectRoot, "internal", "population", "haplotype_cli.py"),
		envPath:    filepath.Join(projectRoot, "phantom_science_env"),
	}
}

// SimulateHaplotypes simulates haplotypes with LD structure for a population
func (hs *HaplotypeSimulator) SimulateHaplotypes(
	population Population,
	chromosome string,
	sampleSize int,
	sequenceLength int,
	randomSeed *int,
	includeMatrix bool,
) (*HaplotypeSimulationResult, error) {
	
	// Prepare command arguments
	args := []string{
		hs.scriptPath,
		"--population", string(population),
		"--chromosome", chromosome,
		"--sample-size", strconv.Itoa(sampleSize),
		"--sequence-length", strconv.Itoa(sequenceLength),
		"--output-format", "json",
	}
	
	if randomSeed != nil {
		args = append(args, "--random-seed", strconv.Itoa(*randomSeed))
	}
	
	if includeMatrix {
		args = append(args, "--include-haplotype-matrix")
	}
	
	// Set up command with virtual environment
	cmd := exec.Command(hs.pythonPath, args...) // #nosec G204
	cmd.Env = append(os.Environ(),
		fmt.Sprintf("VIRTUAL_ENV=%s", hs.envPath),
		fmt.Sprintf("PATH=%s/bin:%s", hs.envPath, os.Getenv("PATH")),
	)
	
	fmt.Printf("🧬 Simulating %s haplotypes for %s...\n", population, chromosome)
	
	// Execute simulation
	output, err := cmd.Output()
	if err != nil {
		return nil, fmt.Errorf("failed to execute haplotype simulation: %w", err)
	}
	
	// Parse JSON result
	var result HaplotypeSimulationResult
	if err := json.Unmarshal(output, &result); err != nil {
		return nil, fmt.Errorf("failed to parse haplotype simulation result: %w", err)
	}
	
	fmt.Printf("   ✅ %s: %d variants, %d LD blocks, r²=%.3f\n",
		population, result.NumVariants, len(result.LDBlocks), result.Statistics.MeanRSquared)
	
	return &result, nil
}

// SimulateAllPopulationsHaplotypes simulates haplotypes for all 11 populations
func (hs *HaplotypeSimulator) SimulateAllPopulationsHaplotypes(
	chromosome string,
	sampleSize int,
	sequenceLength int,
	randomSeed *int,
) (map[Population]*HaplotypeSimulationResult, error) {
	
	populations := GetSupportedPopulations()
	results := make(map[Population]*HaplotypeSimulationResult)
	
	fmt.Printf("🌍 Simulating haplotypes for all 11 populations (%s)...\n", chromosome)
	startTime := time.Now()
	
	for _, pop := range populations {
		result, err := hs.SimulateHaplotypes(
			pop,
			chromosome,
			sampleSize,
			sequenceLength,
			randomSeed,
			false, // Don't include matrix for all populations (too much memory)
		)
		if err != nil {
			return nil, fmt.Errorf("failed to simulate %s haplotypes: %w", pop, err)
		}
		
		results[pop] = result
	}
	
	totalTime := time.Since(startTime)
	fmt.Printf("🎉 All population haplotypes simulated in %.2fs\n", totalTime.Seconds())
	
	return results, nil
}

// ValidateLDPatterns validates LD patterns across populations
func (hs *HaplotypeSimulator) ValidateLDPatterns(
	results map[Population]*HaplotypeSimulationResult,
) (map[string]bool, error) {
	
	fmt.Println("📊 Validating LD patterns across populations...")
	
	validation := make(map[string]bool)
	
	// Extract statistics for analysis
	africanPops := []Population{YRI, LWK, GWD}
	europeanPops := []Population{CEU, GBR, TSI}
	eastAsianPops := []Population{CHB, JPT, CHS}
	
	// Calculate population group statistics
	africanDiversity := calculateGroupMeanDiversity(results, africanPops)
	europeanDiversity := calculateGroupMeanDiversity(results, europeanPops)
	eastAsianDiversity := calculateGroupMeanDiversity(results, eastAsianPops)
	
	// Test 1: African populations should have highest diversity
	validation["african_highest_diversity"] = africanDiversity > europeanDiversity && 
											  africanDiversity > eastAsianDiversity
	
	// Test 2: LD structure should be reasonable
	allRSquared := []float64{}
	for _, result := range results {
		allRSquared = append(allRSquared, result.Statistics.MeanRSquared)
	}
	meanRSquared := calculateMean(allRSquared)
	validation["reasonable_ld_structure"] = meanRSquared > 0.01 && meanRSquared < 0.1
	
	// Test 3: Haplotype blocks should be present
	totalBlocks := 0
	for _, result := range results {
		totalBlocks += len(result.LDBlocks)
	}
	validation["haplotype_blocks_present"] = totalBlocks > 0
	
	// Test 4: Population differentiation in recombination
	recombRates := []float64{}
	for _, result := range results {
		recombRates = append(recombRates, result.Statistics.EffectiveRecombRate)
	}
	recombCV := calculateCV(recombRates)
	validation["recombination_variation"] = recombCV > 0.05
	
	// Test 5: Realistic genetic map lengths
	mapLengths := []float64{}
	var avgSequenceLength float64
	for _, result := range results {
		mapLengths = append(mapLengths, result.Statistics.GeneticMapLengthCM)
		avgSequenceLength += float64(result.SequenceLength)
	}
	if len(results) > 0 {
		avgSequenceLength /= float64(len(results))
	}
	meanMapLength := calculateMean(mapLengths)
	// For sequence, expect ~1.5 cM/Mb average
	expectedRange := avgSequenceLength / 1_000_000 * 1.5 // ~1.5 cM/Mb average
	validation["realistic_map_lengths"] = meanMapLength > expectedRange*0.5 && 
										  meanMapLength < expectedRange*2.0
	
	// Print validation results
	passed := 0
	total := len(validation)
	
	for test, result := range validation {
		status := "✅ PASSED"
		if !result {
			status = "❌ FAILED"
		} else {
			passed++
		}
		fmt.Printf("   %s: %s\n", test, status)
	}
	
	fmt.Printf("\n🎯 LD Pattern Validation: %d/%d tests passed\n", passed, total)
	
	return validation, nil
}

// GetRecombinationMap extracts recombination map from haplotype result
func (hs *HaplotypeSimulator) GetRecombinationMap(
	result *HaplotypeSimulationResult,
) RecombinationMap {
	return result.RecombinationMap
}

// GetLDStatistics extracts LD statistics from haplotype result
func (hs *HaplotypeSimulator) GetLDStatistics(
	result *HaplotypeSimulationResult,
) HaplotypeStatistics {
	return result.Statistics
}

// AnalyzeCrossPopulationLD analyzes LD patterns across populations
func (hs *HaplotypeSimulator) AnalyzeCrossPopulationLD(
	results map[Population]*HaplotypeSimulationResult,
) map[string]interface{} {
	
	analysis := make(map[string]interface{})
	
	// Collect statistics across populations
	diversities := make(map[string]float64)
	ldMetrics := make(map[string]float64)
	blockCounts := make(map[string]int)
	
	for pop, result := range results {
		diversities[string(pop)] = result.Statistics.DiversityPi
		ldMetrics[string(pop)] = result.Statistics.MeanRSquared
		blockCounts[string(pop)] = len(result.LDBlocks)
	}
	
	analysis["diversity_by_population"] = diversities
	analysis["ld_metrics_by_population"] = ldMetrics
	analysis["block_counts_by_population"] = blockCounts
	
	// Calculate summary statistics
	allDiversities := make([]float64, 0, len(diversities))
	allLDMetrics := make([]float64, 0, len(ldMetrics))
	
	for _, div := range diversities {
		allDiversities = append(allDiversities, div)
	}
	for _, ld := range ldMetrics {
		allLDMetrics = append(allLDMetrics, ld)
	}
	
	analysis["diversity_mean"] = calculateMean(allDiversities)
	analysis["diversity_cv"] = calculateCV(allDiversities)
	analysis["ld_mean"] = calculateMean(allLDMetrics)
	analysis["ld_cv"] = calculateCV(allLDMetrics)
	
	return analysis
}

// Helper functions

func calculateGroupMeanDiversity(
	results map[Population]*HaplotypeSimulationResult,
	populations []Population,
) float64 {
	total := 0.0
	count := 0
	
	for _, pop := range populations {
		if result, exists := results[pop]; exists {
			total += result.Statistics.DiversityPi
			count++
		}
	}
	
	if count == 0 {
		return 0.0
	}
	
	return total / float64(count)
}

func calculateMean(values []float64) float64 {
	if len(values) == 0 {
		return 0.0
	}
	
	sum := 0.0
	for _, v := range values {
		sum += v
	}
	
	return sum / float64(len(values))
}

func calculateCV(values []float64) float64 {
	if len(values) <= 1 {
		return 0.0
	}
	
	mean := calculateMean(values)
	if mean == 0.0 {
		return 0.0
	}
	
	sumSquares := 0.0
	for _, v := range values {
		diff := v - mean
		sumSquares += diff * diff
	}
	
	variance := sumSquares / float64(len(values)-1)
	stddev := sqrt(variance)
	
	return stddev / mean
}

// Simple square root implementation
func sqrt(x float64) float64 {
	if x < 0 {
		return 0
	}
	if x == 0 {
		return 0
	}
	
	// Newton's method for square root
	z := x
	for i := 0; i < 10; i++ {
		z = (z + x/z) / 2
	}
	return z
}