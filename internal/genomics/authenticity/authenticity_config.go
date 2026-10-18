package authenticity

import (
	"fmt"
	"path/filepath"
)

// AuthenticityConfig is a configuration schema with validation and presets.
// Nothing outside this file consumes it; the generation code in this package
// is configured through Config instead.
type AuthenticityConfig struct {
	Enabled                  bool    `json:"enabled"`
	ReferenceGenome          string  `json:"reference_genome"`        // "GRCh38", "GRCh39"
	ReferenceGenomePath      string  `json:"reference_genome_path"`   // Path to FASTA file
	Population               string  `json:"population"`              // "EUR", "AFR", "EAS", "SAS", "AMR"
	KmerModelSizes           []int   `json:"kmer_model_sizes"`        // [2, 4, 6, 8]
	JSdivergenceThreshold    float64 `json:"js_divergence_threshold"` // 0.01
	QualityRealism           string  `json:"quality_realism"`         // "basic", "enhanced", "expert"
	ValidationLevel          string  `json:"validation_level"`        // "basic", "comprehensive", "expert"
	EnablePopulationGenetics bool    `json:"enable_population_genetics"`
	EnableQualityPatterns    bool    `json:"enable_quality_patterns"`
	EnableGenomicFeatures    bool    `json:"enable_genomic_features"`
	CacheDirectory           string  `json:"cache_directory"`
	MaxMemoryUsage           int64   `json:"max_memory_usage"` // Bytes
	ParallelProcessing       bool    `json:"parallel_processing"`
	NumWorkers               int     `json:"num_workers"`
}

// PopulationConfig represents population genetics configuration
type PopulationConfig struct {
	PopulationName      string  `json:"population_name"`
	EffectivePopSize    int     `json:"effective_pop_size"`
	MutationRate        float64 `json:"mutation_rate"`
	RecombinationRate   float64 `json:"recombination_rate"`
	TajimaDTarget       float64 `json:"tajima_d_target"`
	NucleotideDiversity float64 `json:"nucleotide_diversity"`
}

// QualityConfig represents quality pattern configuration
type QualityConfig struct {
	Platform                string  `json:"platform"` // "illumina", "pacbio", "nanopore"
	EnableGCBias            bool    `json:"enable_gc_bias"`
	EnableHomopolymerErrors bool    `json:"enable_homopolymer_errors"`
	EnableContextEffects    bool    `json:"enable_context_effects"`
	BaseQuality             float64 `json:"base_quality"`
	QualityDecay            float64 `json:"quality_decay"`
	ErrorRate               float64 `json:"error_rate"`
}

// ValidationConfig represents validation configuration
type ValidationConfig struct {
	EnableJSdivergence    bool    `json:"enable_js_divergence"`
	JSdivergenceThreshold float64 `json:"js_divergence_threshold"`
	EnableWassersteinDist bool    `json:"enable_wasserstein_dist"`
	EnableCramerVonMises  bool    `json:"enable_cramer_von_mises"`
	EnableExpertTesting   bool    `json:"enable_expert_testing"`
	SampleSize            int     `json:"sample_size"`
	ConfidenceLevel       float64 `json:"confidence_level"`
}

// DefaultAuthenticityConfig returns a default configuration
func DefaultAuthenticityConfig() *AuthenticityConfig {
	return &AuthenticityConfig{
		Enabled:                  true,
		ReferenceGenome:          "GRCh38",
		ReferenceGenomePath:      "", // Must be set by user
		Population:               "EUR",
		KmerModelSizes:           []int{2, 4},
		JSdivergenceThreshold:    0.01,
		QualityRealism:           "basic",
		ValidationLevel:          "basic",
		EnablePopulationGenetics: false,
		EnableQualityPatterns:    true,
		EnableGenomicFeatures:    false,
		CacheDirectory:           "./cache/authenticity",
		MaxMemoryUsage:           1024 * 1024 * 1024, // 1GB
		ParallelProcessing:       true,
		NumWorkers:               4,
	}
}

// EnhancedAuthenticityConfig returns DefaultAuthenticityConfig with more
// k-mer sizes, the "enhanced"/"comprehensive" levels and more memory and
// workers.
func EnhancedAuthenticityConfig() *AuthenticityConfig {
	config := DefaultAuthenticityConfig()
	config.KmerModelSizes = []int{2, 4, 6}
	config.QualityRealism = "enhanced"
	config.ValidationLevel = "comprehensive"
	config.EnablePopulationGenetics = true
	config.EnableGenomicFeatures = true
	config.MaxMemoryUsage = 2 * 1024 * 1024 * 1024 // 2GB
	config.NumWorkers = 8
	return config
}

// ExpertAuthenticityConfig returns EnhancedAuthenticityConfig with k up to 8,
// the "expert" levels and a tighter JS threshold.
func ExpertAuthenticityConfig() *AuthenticityConfig {
	config := EnhancedAuthenticityConfig()
	config.KmerModelSizes = []int{2, 4, 6, 8}
	config.QualityRealism = "expert"
	config.ValidationLevel = "expert"
	config.JSdivergenceThreshold = 0.005
	config.MaxMemoryUsage = 4 * 1024 * 1024 * 1024 // 4GB
	config.NumWorkers = 16
	return config
}

// GetPopulationConfig returns hard-coded parameters for the configured
// super-population (EUR if unknown). The values have no cited source and
// are not used by any model.
func (c *AuthenticityConfig) GetPopulationConfig() *PopulationConfig {
	configs := map[string]*PopulationConfig{
		"EUR": {
			PopulationName:      "European",
			EffectivePopSize:    10000,
			MutationRate:        2.3e-8,
			RecombinationRate:   1.2e-8,
			TajimaDTarget:       0.0,
			NucleotideDiversity: 0.001,
		},
		"AFR": {
			PopulationName:      "African",
			EffectivePopSize:    12000,
			MutationRate:        2.3e-8,
			RecombinationRate:   1.2e-8,
			TajimaDTarget:       0.2,
			NucleotideDiversity: 0.0015,
		},
		"EAS": {
			PopulationName:      "East Asian",
			EffectivePopSize:    8000,
			MutationRate:        2.3e-8,
			RecombinationRate:   1.2e-8,
			TajimaDTarget:       -0.3,
			NucleotideDiversity: 0.0008,
		},
		"SAS": {
			PopulationName:      "South Asian",
			EffectivePopSize:    9000,
			MutationRate:        2.3e-8,
			RecombinationRate:   1.2e-8,
			TajimaDTarget:       -0.1,
			NucleotideDiversity: 0.0009,
		},
		"AMR": {
			PopulationName:      "Admixed American",
			EffectivePopSize:    7000,
			MutationRate:        2.3e-8,
			RecombinationRate:   1.2e-8,
			TajimaDTarget:       -0.2,
			NucleotideDiversity: 0.0012,
		},
	}

	if config, exists := configs[c.Population]; exists {
		return config
	}

	// Default to European if population not found
	return configs["EUR"]
}

// GetQualityConfig returns quality configuration based on platform
func (c *AuthenticityConfig) GetQualityConfig(platform string) *QualityConfig {
	configs := map[string]*QualityConfig{
		"illumina": {
			Platform:                "illumina",
			EnableGCBias:            true,
			EnableHomopolymerErrors: false,
			EnableContextEffects:    true,
			BaseQuality:             30.0,
			QualityDecay:            0.1,
			ErrorRate:               0.001,
		},
		"pacbio": {
			Platform:                "pacbio",
			EnableGCBias:            false,
			EnableHomopolymerErrors: true,
			EnableContextEffects:    true,
			BaseQuality:             20.0,
			QualityDecay:            0.05,
			ErrorRate:               0.01,
		},
		"nanopore": {
			Platform:                "nanopore",
			EnableGCBias:            false,
			EnableHomopolymerErrors: true,
			EnableContextEffects:    true,
			BaseQuality:             10.0,
			QualityDecay:            0.2,
			ErrorRate:               0.05,
		},
	}

	if config, exists := configs[platform]; exists {
		return config
	}

	// Default to Illumina if platform not found
	return configs["illumina"]
}

// GetValidationConfig returns validation configuration
func (c *AuthenticityConfig) GetValidationConfig() *ValidationConfig {
	configs := map[string]*ValidationConfig{
		"basic": {
			EnableJSdivergence:    true,
			JSdivergenceThreshold: c.JSdivergenceThreshold,
			EnableWassersteinDist: false,
			EnableCramerVonMises:  false,
			EnableExpertTesting:   false,
			SampleSize:            1000,
			ConfidenceLevel:       0.95,
		},
		"comprehensive": {
			EnableJSdivergence:    true,
			JSdivergenceThreshold: c.JSdivergenceThreshold,
			EnableWassersteinDist: true,
			EnableCramerVonMises:  true,
			EnableExpertTesting:   false,
			SampleSize:            5000,
			ConfidenceLevel:       0.95,
		},
		"expert": {
			EnableJSdivergence:    true,
			JSdivergenceThreshold: c.JSdivergenceThreshold,
			EnableWassersteinDist: true,
			EnableCramerVonMises:  true,
			EnableExpertTesting:   true,
			SampleSize:            10000,
			ConfidenceLevel:       0.99,
		},
	}

	if config, exists := configs[c.ValidationLevel]; exists {
		return config
	}

	// Default to basic if validation level not found
	return configs["basic"]
}

// Validate validates the authenticity configuration
func (c *AuthenticityConfig) Validate() error {
	if !c.Enabled {
		return nil // No validation needed if disabled
	}

	// Validate reference genome
	if c.ReferenceGenome != "GRCh38" && c.ReferenceGenome != "GRCh39" {
		return fmt.Errorf("unsupported reference genome: %s", c.ReferenceGenome)
	}

	// Validate reference genome path
	if c.ReferenceGenomePath == "" {
		return fmt.Errorf("reference genome path must be specified")
	}

	// Validate population
	validPopulations := []string{"EUR", "AFR", "EAS", "SAS", "AMR"}
	validPop := false
	for _, pop := range validPopulations {
		if c.Population == pop {
			validPop = true
			break
		}
	}
	if !validPop {
		return fmt.Errorf("unsupported population: %s", c.Population)
	}

	// Validate k-mer model sizes
	if len(c.KmerModelSizes) == 0 {
		return fmt.Errorf("at least one k-mer model size must be specified")
	}
	for _, k := range c.KmerModelSizes {
		if k < 2 || k > 8 {
			return fmt.Errorf("k-mer size must be between 2 and 8, got: %d", k)
		}
	}

	// Validate JS divergence threshold
	if c.JSdivergenceThreshold <= 0 || c.JSdivergenceThreshold > 1 {
		return fmt.Errorf("JS divergence threshold must be between 0 and 1, got: %f", c.JSdivergenceThreshold)
	}

	// Validate quality realism level
	validQualityLevels := []string{"basic", "enhanced", "expert"}
	validQuality := false
	for _, level := range validQualityLevels {
		if c.QualityRealism == level {
			validQuality = true
			break
		}
	}
	if !validQuality {
		return fmt.Errorf("unsupported quality realism level: %s", c.QualityRealism)
	}

	// Validate validation level
	validValidationLevels := []string{"basic", "comprehensive", "expert"}
	validValidation := false
	for _, level := range validValidationLevels {
		if c.ValidationLevel == level {
			validValidation = true
			break
		}
	}
	if !validValidation {
		return fmt.Errorf("unsupported validation level: %s", c.ValidationLevel)
	}

	// Validate memory usage
	if c.MaxMemoryUsage < 100*1024*1024 { // Minimum 100MB
		return fmt.Errorf("max memory usage too low: %d bytes", c.MaxMemoryUsage)
	}

	// Validate number of workers
	if c.NumWorkers < 1 || c.NumWorkers > 64 {
		return fmt.Errorf("number of workers must be between 1 and 64, got: %d", c.NumWorkers)
	}

	return nil
}

// GetCacheFilePath returns the path for a cache file
func (c *AuthenticityConfig) GetCacheFilePath(filename string) string {
	return filepath.Join(c.CacheDirectory, filename)
}

// GetKmerCachePath returns the path for k-mer cache files
func (c *AuthenticityConfig) GetKmerCachePath(k int) string {
	filename := fmt.Sprintf("kmer_%s_%s_k%d.cache", c.ReferenceGenome, c.Population, k)
	return c.GetCacheFilePath(filename)
}

// GetModelCachePath returns the path for model cache files
func (c *AuthenticityConfig) GetModelCachePath(modelType string) string {
	filename := fmt.Sprintf("model_%s_%s_%s.cache", modelType, c.ReferenceGenome, c.Population)
	return c.GetCacheFilePath(filename)
}

// IsHighMemoryMode returns true if configuration requires high memory usage
func (c *AuthenticityConfig) IsHighMemoryMode() bool {
	return c.MaxMemoryUsage > 2*1024*1024*1024 // > 2GB
}

// IsHighPerformanceMode reports whether parallel processing is on with at
// least 8 workers.
func (c *AuthenticityConfig) IsHighPerformanceMode() bool {
	return c.ParallelProcessing && c.NumWorkers >= 8
}

// GetMaxKmerSize returns the maximum k-mer size configured
func (c *AuthenticityConfig) GetMaxKmerSize() int {
	maxK := 0
	for _, k := range c.KmerModelSizes {
		if k > maxK {
			maxK = k
		}
	}
	return maxK
}

// GetMinKmerSize returns the minimum k-mer size configured
func (c *AuthenticityConfig) GetMinKmerSize() int {
	if len(c.KmerModelSizes) == 0 {
		return 2
	}

	minK := c.KmerModelSizes[0]
	for _, k := range c.KmerModelSizes {
		if k < minK {
			minK = k
		}
	}
	return minK
}
