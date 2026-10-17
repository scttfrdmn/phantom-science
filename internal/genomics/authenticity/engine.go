package authenticity

import (
	"fmt"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/errors"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/reference"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/statistics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/validation"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Engine represents the Genomics Statistical Authenticity Engine
type Engine struct {
	config          *Config
	referenceGenome *reference.GRCh38
	lazyLoader      *reference.LazyLoadingManager
	kmerExtractor   *kmer.FrequencyExtractor
	jsCalculator    *statistics.JensenShannonDivergence
	validator       *validation.ValidationFramework
	
	// Cached frequency tables for different k-mer sizes
	referenceKmerTables map[int]*kmer.FrequencyTable
	
	// Thread safety
	mu sync.RWMutex
	
	// Performance metrics
	initialized   bool
	lastValidated time.Time
}

// Config contains configuration for the authenticity engine
type Config struct {
	// Reference genome configuration
	ReferenceGenomePath string
	Assembly            string // "GRCh38" or "GRCh39"
	
	// K-mer analysis configuration
	KmerSizes           []int   // Default: [2, 3, 4, 5, 6, 7, 8]
	UsePopulationData   bool    // Use population variation data
	PopulationName      string  // e.g., "1000G_EUR"
	
	// Validation configuration
	JSThreshold         float64 // Jensen-Shannon divergence threshold
	WassersteinThreshold float64 // Wasserstein distance threshold
	EnableValidation    bool    // Enable real-time validation
	ValidationInterval  time.Duration // How often to validate
	
	// Performance configuration
	MaxLoadedChromosomes int  // Maximum chromosomes to keep in memory
	EnableLazyLoading    bool // Enable lazy loading of chromosome data
	ParallelProcessing   bool // Enable parallel k-mer processing
	
	// Integration configuration
	IntegrateWithFingerprinting bool // Maintain fingerprinting compatibility
	QualityScoreAuthenticity    bool // Apply authenticity to quality scores
	
	// Debugging and logging
	EnableDebugLogging bool
	LogLevel          string // "debug", "info", "warn", "error"
}

// DefaultConfig returns a default configuration for the authenticity engine
func DefaultConfig() *Config {
	return &Config{
		ReferenceGenomePath:         "/data/reference/GRCh38",
		Assembly:                    "GRCh38",
		KmerSizes:                   []int{2, 3, 4, 5, 6, 7, 8},
		UsePopulationData:           false,
		PopulationName:              "",
		JSThreshold:                 0.01,
		WassersteinThreshold:        0.02,
		EnableValidation:            true,
		ValidationInterval:          5 * time.Minute,
		MaxLoadedChromosomes:        3,
		EnableLazyLoading:           true,
		ParallelProcessing:          true,
		IntegrateWithFingerprinting: true,
		QualityScoreAuthenticity:    true,
		EnableDebugLogging:          false,
		LogLevel:                    "info",
	}
}

// NewEngine creates a new authenticity engine
func NewEngine(config *Config) *Engine {
	if config == nil {
		config = DefaultConfig()
	}
	
	// Initialize reference genome
	refGenome := reference.NewGRCh38(config.ReferenceGenomePath)
	
	// Initialize lazy loading manager
	var lazyLoader *reference.LazyLoadingManager
	if config.EnableLazyLoading {
		lazyLoader = reference.NewLazyLoadingManager(refGenome, config.MaxLoadedChromosomes)
	}
	
	// Initialize k-mer extractor (default to 4-mer)
	kmerExtractor := kmer.NewFrequencyExtractor(4)
	
	// Initialize validation framework
	validationConfig := validation.ValidationConfig{
		JSThreshold:          config.JSThreshold,
		WassersteinThreshold: config.WassersteinThreshold,
		EnableDetailedReport: config.EnableDebugLogging,
		EnableRecommendations: true,
	}
	validator := validation.NewValidationFrameworkWithConfig(validationConfig)
	
	return &Engine{
		config:                config,
		referenceGenome:       refGenome,
		lazyLoader:           lazyLoader,
		kmerExtractor:        kmerExtractor,
		jsCalculator:         statistics.NewJensenShannonDivergence(),
		validator:            validator,
		referenceKmerTables:  make(map[int]*kmer.FrequencyTable),
		initialized:          false,
	}
}

// Initialize prepares the authenticity engine for use
func (e *Engine) Initialize() error {
	e.mu.Lock()
	defer e.mu.Unlock()
	
	if e.initialized {
		return nil
	}
	
	// Load reference k-mer frequencies for all configured k-mer sizes
	if err := e.loadReferenceKmerFrequencies(); err != nil {
		return errors.WrapError(err, errors.InitializationError, "authenticity_engine", "initialize", "failed to load reference k-mer frequencies")
	}
	
	e.initialized = true
	return nil
}

// loadReferenceKmerFrequencies loads k-mer frequencies from the reference genome
func (e *Engine) loadReferenceKmerFrequencies() error {
	// For now, load data from a subset of chromosomes (chr1, chr2, chr3)
	// In a full implementation, this would load from all chromosomes
	testChromosomes := []string{"1", "2", "3"}
	
	for _, k := range e.config.KmerSizes {
		// Configure k-mer extractor for this k-mer size
		e.kmerExtractor.K = k
		
		var allSequences [][]byte
		
		// Load sequences from test chromosomes
		for _, chrName := range testChromosomes {
			var chr *reference.Chromosome
			var err error
			
			if e.config.EnableLazyLoading {
				chr, err = e.lazyLoader.LoadChromosome(chrName)
			} else {
				chr, err = e.referenceGenome.LoadChromosome(chrName)
			}
			
			if err != nil {
				return errors.WrapError(err, errors.ReferenceLoadError, "authenticity_engine", "load_reference_kmers", fmt.Sprintf("failed to load chromosome %s", chrName))
			}
			
			// Sample a subset of the chromosome for performance
			sampleSize := 100000 // Sample 100kb
			if len(chr.Sequence) > sampleSize {
				allSequences = append(allSequences, chr.Sequence[:sampleSize])
			} else {
				allSequences = append(allSequences, chr.Sequence)
			}
		}
		
		// Extract k-mer frequencies
		freqTable, err := e.kmerExtractor.ExtractFromMultipleSequences(allSequences)
		if err != nil {
			return errors.WrapError(err, errors.KmerExtractionError, "authenticity_engine", "load_reference_kmers", fmt.Sprintf("failed to extract %d-mer frequencies", k))
		}
		
		e.referenceKmerTables[k] = freqTable
	}
	
	return nil
}

// GenerateAuthenticSequence generates a biologically authentic DNA sequence
func (e *Engine) GenerateAuthenticSequence(length int, chromosome string) ([]byte, error) {
	if !e.initialized {
		if err := e.Initialize(); err != nil {
			return nil, errors.WrapError(err, errors.InitializationError, "authenticity_engine", "generate_authentic_sequence", "failed to initialize authenticity engine")
		}
	}
	
	e.mu.RLock()
	defer e.mu.RUnlock()
	
	// For Phase 3a.1, we'll use a simple approach based on 4-mer frequencies
	// This will be enhanced in later phases with HMM and population genetics
	
	refTable, exists := e.referenceKmerTables[4]
	if !exists {
		return nil, errors.NewResourceError("authenticity_engine", "generate_authentic_sequence", "reference 4-mer frequency table not available")
	}
	
	// Get probability distribution for 4-mers
	probDist := refTable.GetProbabilityDistribution()
	
	// Generate sequence using weighted random selection
	sequence := make([]byte, length)
	
	// Start with a random 4-mer
	startKmer := e.selectRandomKmer(probDist)
	if len(startKmer) < 4 {
		return nil, errors.NewValidationError("authenticity_engine", "generate_authentic_sequence", fmt.Sprintf("invalid starting k-mer: %s", startKmer))
	}
	
	// Copy the starting k-mer
	copy(sequence[0:4], startKmer)
	
	// Generate remaining bases using Markov chain approach
	for i := 4; i < length; i++ {
		// Get the previous 3 bases (for 4-mer context)
		context := string(sequence[i-3:i])
		
		// Find all 4-mers that start with this context
		candidates := make(map[string]float64)
		for kmer, prob := range probDist {
			if len(kmer) >= 4 && kmer[:3] == context {
				candidates[kmer] = prob
			}
		}
		
		if len(candidates) == 0 {
			// Fallback to random base if no candidates found
			bases := []byte{'A', 'T', 'G', 'C'}
			sequence[i] = bases[i%4]
			continue
		}
		
		// Select next base based on weighted probabilities
		selectedKmer := e.selectRandomKmer(candidates)
		if len(selectedKmer) >= 4 {
			sequence[i] = selectedKmer[3]
		} else {
			sequence[i] = 'A' // Fallback
		}
	}
	
	return sequence, nil
}

// selectRandomKmer selects a k-mer based on probability distribution
func (e *Engine) selectRandomKmer(probDist map[string]float64) string {
	// Simple implementation for Phase 3a.1
	// In later phases, this will use proper weighted random selection
	
	// For now, return the most frequent k-mer
	var maxKmer string
	var maxProb float64
	
	for kmer, prob := range probDist {
		if prob > maxProb {
			maxProb = prob
			maxKmer = kmer
		}
	}
	
	return maxKmer
}

// ValidateSequence validates a generated sequence against reference authenticity metrics
func (e *Engine) ValidateSequence(sequence []byte) (*validation.ValidationReport, error) {
	if !e.initialized {
		if err := e.Initialize(); err != nil {
			return nil, errors.WrapError(err, errors.InitializationError, "authenticity_engine", "validate_sequence", "failed to initialize authenticity engine")
		}
	}
	
	e.mu.RLock()
	defer e.mu.RUnlock()
	
	// Extract k-mer frequencies from the generated sequence
	var syntheticTables []*kmer.FrequencyTable
	var referenceTables []*kmer.FrequencyTable
	
	for _, k := range e.config.KmerSizes {
		// Configure extractor for this k-mer size
		e.kmerExtractor.K = k
		
		// Extract synthetic k-mer frequencies
		synTable, err := e.kmerExtractor.ExtractFromSequence(sequence)
		if err != nil {
			return nil, errors.WrapError(err, errors.KmerExtractionError, "authenticity_engine", "validate_sequence", fmt.Sprintf("failed to extract %d-mer frequencies from synthetic sequence", k))
		}
		
		// Get reference k-mer frequencies
		refTable, exists := e.referenceKmerTables[k]
		if !exists {
			return nil, errors.NewResourceError("authenticity_engine", "validate_sequence", fmt.Sprintf("reference %d-mer frequency table not available", k))
		}
		
		syntheticTables = append(syntheticTables, synTable)
		referenceTables = append(referenceTables, refTable)
	}
	
	// Validate using the validation framework
	report := e.validator.ValidateMultipleKmers(referenceTables, syntheticTables)
	
	// Update last validated timestamp
	e.lastValidated = time.Now()
	
	return report, nil
}

// GetAuthenticityScore calculates an overall authenticity score for a sequence
func (e *Engine) GetAuthenticityScore(sequence []byte) (float64, error) {
	report, err := e.ValidateSequence(sequence)
	if err != nil {
		return 0.0, err
	}
	
	// Calculate authenticity score based on validation results
	if report.TotalTests == 0 {
		return 0.0, errors.NewValidationError("authenticity_engine", "get_authenticity_score", "no validation tests were performed")
	}
	
	// Authenticity score is based on percentage of passed tests
	// and inverse of average score (lower scores are better)
	passRate := float64(report.PassedTests) / float64(report.TotalTests)
	
	// Normalize average score (assuming lower is better)
	normalizedScore := 1.0 - (report.AverageScore / 1.0)
	
	// Combine pass rate and normalized score
	authenticityScore := (passRate + normalizedScore) / 2.0
	
	return authenticityScore, nil
}

// GetReferenceStatistics returns statistics about the reference genome data
func (e *Engine) GetReferenceStatistics() (map[string]interface{}, error) {
	if !e.initialized {
		if err := e.Initialize(); err != nil {
			return nil, errors.WrapError(err, errors.InitializationError, "authenticity_engine", "get_reference_statistics", "failed to initialize authenticity engine")
		}
	}
	
	e.mu.RLock()
	defer e.mu.RUnlock()
	
	stats := make(map[string]interface{})
	
	// Add k-mer statistics for each k-mer size
	for k, table := range e.referenceKmerTables {
		kmerStats := table.GetKmerStats()
		stats[fmt.Sprintf("kmer_%d", k)] = map[string]interface{}{
			"total_kmers":     kmerStats.TotalKmers,
			"unique_kmers":    kmerStats.UniqueKmers,
			"complexity":      kmerStats.Complexity,
			"shannon_entropy": kmerStats.ShannonEntropy,
			"gc_content":      kmerStats.GCContent,
		}
	}
	
	// Add general statistics
	stats["assembly"] = e.config.Assembly
	stats["loaded_chromosomes"] = len(e.referenceGenome.Chromosomes)
	stats["last_validated"] = e.lastValidated
	stats["initialized"] = e.initialized
	
	return stats, nil
}

// UpdateFromSequencerConfig updates the authenticity engine with sequencer configuration
func (e *Engine) UpdateFromSequencerConfig(config types.SequencerConfig) error {
	e.mu.Lock()
	defer e.mu.Unlock()
	
	// Update authenticity-specific parameters based on sequencer config
	// This allows integration with existing sequencer simulation
	
	if config.SequencerModel != "" {
		// Adjust authenticity parameters based on sequencing platform
		switch config.SequencerModel {
		case "illumina":
			e.config.JSThreshold = 0.01 // Strict for Illumina
		case "pacbio":
			e.config.JSThreshold = 0.02 // Slightly relaxed for PacBio
		case "nanopore":
			e.config.JSThreshold = 0.03 // More relaxed for Nanopore
		}
	}
	
	// Update validation interval based on read length
	if config.ReadLength > 0 {
		if config.ReadLength < 100 {
			e.config.ValidationInterval = 1 * time.Minute // More frequent for short reads
		} else {
			e.config.ValidationInterval = 5 * time.Minute // Less frequent for long reads
		}
	}
	
	return nil
}

// Close releases resources used by the authenticity engine
func (e *Engine) Close() error {
	e.mu.Lock()
	defer e.mu.Unlock()
	
	// Close reference genome
	if e.referenceGenome != nil {
		if err := e.referenceGenome.Close(); err != nil {
			return errors.WrapError(err, errors.ResourceError, "authenticity_engine", "close", "failed to close reference genome")
		}
	}
	
	// Clear cached data
	e.referenceKmerTables = make(map[int]*kmer.FrequencyTable)
	e.initialized = false
	
	return nil
}

// GetStatus returns the current status of the authenticity engine
func (e *Engine) GetStatus() map[string]interface{} {
	e.mu.RLock()
	defer e.mu.RUnlock()
	
	status := map[string]interface{}{
		"initialized":           e.initialized,
		"assembly":             e.config.Assembly,
		"kmer_sizes":           e.config.KmerSizes,
		"js_threshold":         e.config.JSThreshold,
		"validation_enabled":   e.config.EnableValidation,
		"lazy_loading_enabled": e.config.EnableLazyLoading,
		"last_validated":       e.lastValidated,
		"reference_tables":     len(e.referenceKmerTables),
	}
	
	return status
}

// IsInitialized returns whether the engine is initialized
func (e *Engine) IsInitialized() bool {
	e.mu.RLock()
	defer e.mu.RUnlock()
	return e.initialized
}

// GetConfig returns the current configuration
func (e *Engine) GetConfig() *Config {
	e.mu.RLock()
	defer e.mu.RUnlock()
	
	// Return a copy to prevent external modification
	configCopy := *e.config
	return &configCopy
}