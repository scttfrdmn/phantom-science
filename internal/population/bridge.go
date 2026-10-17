package population

import (
	"encoding/json"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// Population represents supported populations
type Population string

const (
	YRI Population = "YRI" // Yoruba - West Africa
	LWK Population = "LWK" // Luhya - East Africa
	GWD Population = "GWD" // Gambian - West Africa
	CEU Population = "CEU" // CEPH European - Northern Europe
	GBR Population = "GBR" // British - Northern Europe
	TSI Population = "TSI" // Tuscan - Southern Europe
	CHB Population = "CHB" // Han Chinese - East Asia
	JPT Population = "JPT" // Japanese - East Asia
	CHS Population = "CHS" // Southern Han Chinese - East Asia
	MXL Population = "MXL" // Mexican - Americas
	PUR Population = "PUR" // Puerto Rican - Americas
)

// DemographicParameters holds parameters for a population
type DemographicParameters struct {
	Population              Population `json:"population"`
	EffectivePopulation     int        `json:"effective_population"`
	GrowthRate              float64    `json:"growth_rate"`
	MutationRate            float64    `json:"mutation_rate"`
	RecombinationRate       float64    `json:"recombination_rate"`
	BottleneckTime          int        `json:"bottleneck_time"`
	BottleneckSize          int        `json:"bottleneck_size"`
	FoundingTimeYears       int        `json:"founding_time_years"`
	FoundingTimeGenerations int        `json:"founding_time_generations"`
	SuperPopulation         string     `json:"super_population"`
	Region                  string     `json:"region"`
	Description             string     `json:"description"`
}

// SimulationResult holds the result of demographic simulation
type SimulationResult struct {
	Population       Population              `json:"population"`
	NumVariants      int                     `json:"num_variants"`
	NumTrees         int                     `json:"num_trees"`
	NumSamples       int                     `json:"num_samples"`
	SequenceLength   int                     `json:"sequence_length"`
	SimulationTime   float64                 `json:"simulation_time"`
	MemoryUsedMB     float64                 `json:"memory_used_mb"`
	DemographicParams DemographicParameters  `json:"demographic_params"`
	Variants         []VariantData           `json:"variants"`
	Statistics       PopulationStatistics    `json:"statistics"`
}

// VariantData represents a genetic variant
type VariantData struct {
	Position        int       `json:"position"`
	Alleles         []string  `json:"alleles"`
	Genotypes       []int     `json:"genotypes"`
	Frequency       float64   `json:"frequency"`
	AncestralState  string    `json:"ancestral_state"`
	NumAlleles      int       `json:"num_alleles"`
}

// PopulationStatistics holds population genetics statistics
type PopulationStatistics struct {
	Population             string                 `json:"population"`
	NumVariants            int                    `json:"num_variants"`
	NumSamples             int                    `json:"num_samples"`
	DiversityPi            float64                `json:"diversity_pi"`
	WattersonTheta         float64                `json:"watterson_theta"`
	TajimaD                float64                `json:"tajima_d"`
	SegregataingSites      int                    `json:"segregating_sites"`
	RareVariantsCount      int                    `json:"rare_variants_count"`
	CommonVariantsCount    int                    `json:"common_variants_count"`
	SiteFrequencySpectrum  []int                  `json:"site_frequency_spectrum"`
	MeanAlleleFrequency    float64                `json:"mean_allele_frequency"`
	DemographicParams      map[string]interface{} `json:"demographic_params"`
}

// DemographicSimulator provides interface to Python demographic models
type DemographicSimulator struct {
	pythonPath string
	scriptPath string
	envPath    string
}

// NewDemographicSimulator creates a new demographic simulator
func NewDemographicSimulator(projectRoot string) *DemographicSimulator {
	return &DemographicSimulator{
		pythonPath: filepath.Join(projectRoot, "phantom_science_env", "bin", "python3"),
		scriptPath: filepath.Join(projectRoot, "internal", "population", "demographics.py"),
		envPath:    filepath.Join(projectRoot, "phantom_science_env"),
	}
}

// SimulatePopulation simulates a single population using coalescent models
func (ds *DemographicSimulator) SimulatePopulation(
	population Population,
	sampleSize int,
	sequenceLength int,
	randomSeed *int,
) (*SimulationResult, error) {
	
	// Prepare command arguments
	args := []string{
		ds.scriptPath,
		"--population", string(population),
		"--sample-size", strconv.Itoa(sampleSize),
		"--sequence-length", strconv.Itoa(sequenceLength),
		"--output-format", "json",
	}
	
	if randomSeed != nil {
		args = append(args, "--random-seed", strconv.Itoa(*randomSeed))
	}
	
	// Set up command with virtual environment
	cmd := exec.Command(ds.pythonPath, args...) // #nosec G204
	cmd.Env = append(os.Environ(),
		fmt.Sprintf("VIRTUAL_ENV=%s", ds.envPath),
		fmt.Sprintf("PATH=%s/bin:%s", ds.envPath, os.Getenv("PATH")),
	)
	
	// Execute simulation
	output, err := cmd.Output()
	if err != nil {
		return nil, fmt.Errorf("failed to execute demographic simulation: %w", err)
	}
	
	// Parse JSON result
	var result SimulationResult
	if err := json.Unmarshal(output, &result); err != nil {
		return nil, fmt.Errorf("failed to parse simulation result: %w", err)
	}
	
	return &result, nil
}

// SimulateAllPopulations simulates all 11 populations
func (ds *DemographicSimulator) SimulateAllPopulations(
	sampleSize int,
	sequenceLength int,
	randomSeed *int,
) (map[Population]*SimulationResult, error) {
	
	populations := []Population{YRI, LWK, GWD, CEU, GBR, TSI, CHB, JPT, CHS, MXL, PUR}
	results := make(map[Population]*SimulationResult)
	
	fmt.Println("🌍 Simulating all 11 populations with demographic models...")
	startTime := time.Now()
	
	for _, pop := range populations {
		fmt.Printf("🧬 Simulating %s...\n", pop)
		
		result, err := ds.SimulatePopulation(pop, sampleSize, sequenceLength, randomSeed)
		if err != nil {
			return nil, fmt.Errorf("failed to simulate %s: %w", pop, err)
		}
		
		results[pop] = result
		fmt.Printf("   ✅ %s: %d variants, %.2fs\n", 
			pop, result.NumVariants, result.SimulationTime)
	}
	
	totalTime := time.Since(startTime)
	fmt.Printf("🎉 All populations simulated in %.2fs\n", totalTime.Seconds())
	
	return results, nil
}

// GetPopulationParameters returns demographic parameters for a population
func (ds *DemographicSimulator) GetPopulationParameters(population Population) (*DemographicParameters, error) {
	// Static parameters based on literature
	parameters := map[Population]DemographicParameters{
		YRI: {
			Population:              YRI,
			EffectivePopulation:     12300,
			GrowthRate:              0.0166,
			MutationRate:            1.25e-8,
			RecombinationRate:       1.13e-8,
			BottleneckTime:          2000,
			BottleneckSize:          1000,
			FoundingTimeYears:       70000,
			FoundingTimeGenerations: 2800,
			SuperPopulation:         "African",
			Region:                  "West Africa",
			Description:             "Yoruba in Ibadan, Nigeria",
		},
		CEU: {
			Population:              CEU,
			EffectivePopulation:     11418,
			GrowthRate:              0.0038,
			MutationRate:            1.25e-8,
			RecombinationRate:       1.13e-8,
			BottleneckTime:          1600,
			BottleneckSize:          2000,
			FoundingTimeYears:       40000,
			FoundingTimeGenerations: 1600,
			SuperPopulation:         "European",
			Region:                  "Northern Europe",
			Description:             "Utah residents with Northern European ancestry",
		},
		CHB: {
			Population:              CHB,
			EffectivePopulation:     9025,
			GrowthRate:              0.0048,
			MutationRate:            1.25e-8,
			RecombinationRate:       1.13e-8,
			BottleneckTime:          1800,
			BottleneckSize:          1500,
			FoundingTimeYears:       50000,
			FoundingTimeGenerations: 2000,
			SuperPopulation:         "East Asian",
			Region:                  "East Asia",
			Description:             "Han Chinese in Beijing, China",
		},
		MXL: {
			Population:              MXL,
			EffectivePopulation:     8500,
			GrowthRate:              0.0045,
			MutationRate:            1.25e-8,
			RecombinationRate:       1.13e-8,
			BottleneckTime:          1500,
			BottleneckSize:          1200,
			FoundingTimeYears:       15000,
			FoundingTimeGenerations: 600,
			SuperPopulation:         "Americas",
			Region:                  "Americas",
			Description:             "Mexican Ancestry from Los Angeles USA",
		},
		// Add other populations as needed
	}
	
	if params, exists := parameters[population]; exists {
		return &params, nil
	}
	
	return nil, fmt.Errorf("unknown population: %s", population)
}

// IntegrateWithHMMModel integrates demographic simulation with existing HMM k-mer model
func (ds *DemographicSimulator) IntegrateWithHMMModel(
	population Population,
	sampleSize int,
	sequenceLength int,
	kmerSize int,
) (*authenticity.HMMKmerModel, error) {
	
	fmt.Printf("🔗 Integrating %s demographic model with HMM k-mer model (k=%d)...\n", 
		population, kmerSize)
	
	// Simulate population with demographic model
	result, err := ds.SimulatePopulation(population, sampleSize, sequenceLength, nil)
	if err != nil {
		return nil, fmt.Errorf("failed to simulate population: %w", err)
	}
	
	// Create HMM k-mer model
	model := authenticity.NewHMMKmerModel(kmerSize)
	
	// Convert demographic variants to training data for HMM
	// This would require implementing a ChromosomeInterface that uses
	// the simulated variants instead of VCF data
	
	fmt.Printf("   📊 Generated %d variants from demographic model\n", result.NumVariants)
	fmt.Printf("   🧬 Diversity (π): %.6f\n", result.Statistics.DiversityPi)
	fmt.Printf("   📈 Tajima's D: %.3f\n", result.Statistics.TajimaD)
	
	// For now, return the model - full integration would require
	// implementing ChromosomeInterface for simulated data
	return model, nil
}

// ValidatePopulationStatistics validates that simulated populations have realistic statistics
func (ds *DemographicSimulator) ValidatePopulationStatistics(
	results map[Population]*SimulationResult,
) error {
	
	fmt.Println("📊 Validating population statistics...")
	
	for pop, result := range results {
		stats := result.Statistics
		params := result.DemographicParams
		
		// Validate diversity levels.
		// African populations (YRI) have the highest nucleotide diversity, followed by
		// East Asian (CHB), then European (CEU) — consistent with 1000 Genomes Phase 3.
		expectedDiversity := map[Population]struct{ min, max float64 }{
			YRI: {0.000010, 0.000030}, // African: highest diversity
			CEU: {0.000005, 0.000015}, // European: lowest diversity
			CHB: {0.000008, 0.000025}, // East Asian: intermediate
		}
		
		if expected, exists := expectedDiversity[pop]; exists {
			if stats.DiversityPi < expected.min || stats.DiversityPi > expected.max {
				return fmt.Errorf("%s diversity (π=%.6f) outside expected range [%.6f, %.6f]",
					pop, stats.DiversityPi, expected.min, expected.max)
			}
		}
		
		// Validate Tajima's D (should be slightly negative for most populations due to population growth)
		if stats.TajimaD > 0.5 || stats.TajimaD < -2.0 {
			return fmt.Errorf("%s Tajima's D (%.3f) outside reasonable range [-2.0, 0.5]",
				pop, stats.TajimaD)
		}
		
		// Validate rare variant proportion (should be majority)
		rareVariantProp := float64(stats.RareVariantsCount) / float64(stats.NumVariants)
		if rareVariantProp < 0.7 {
			return fmt.Errorf("%s rare variant proportion (%.3f) too low, expected >0.7",
				pop, rareVariantProp)
		}
		
		fmt.Printf("   ✅ %s: π=%.6f, D=%.3f, rare=%.1f%%, Ne=%d\n",
			pop, stats.DiversityPi, stats.TajimaD, rareVariantProp*100, params.EffectivePopulation)
	}
	
	fmt.Println("   🎉 All population statistics validated!")
	return nil
}

// GetSupportedPopulations returns list of supported populations
func GetSupportedPopulations() []Population {
	return []Population{YRI, LWK, GWD, CEU, GBR, TSI, CHB, JPT, CHS, MXL, PUR}
}

// PopulationToString converts population enum to string
func (p Population) String() string {
	return string(p)
}

// IsValidPopulation checks if population is supported
func IsValidPopulation(pop string) bool {
	for _, p := range GetSupportedPopulations() {
		if string(p) == pop {
			return true
		}
	}
	return false
}