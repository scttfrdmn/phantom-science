package vcf

import (
	"fmt"
	"math"
	"strings"
	"sync"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// AuthenticityProcessor processes VCF variants for authenticity analysis
type AuthenticityProcessor struct {
	kmerSize           int
	contextWindowSize  int
	authenticityEnabled bool
	kmerExtractor      *kmer.PackedFrequencyExtractor
	gpuProcessor       *kmer.GPUPackedProcessor
	useGPU             bool
	
	// Aggregated results
	variantCount        int64
	authenticityScores  []float64
	kmerProfiles        []*kmer.PackedFrequencyTable
	populationContexts  map[string]*PopulationContext
	
	// Thread safety
	mutex sync.RWMutex
	
	// Statistics
	processingStats     AuthenticityStats
}

// PopulationContext contains population-specific variant context
type PopulationContext struct {
	PopulationID      string
	VariantFrequency  map[string]float64
	KmerProfile       *kmer.PackedFrequencyTable
	AuthenticityScore float64
	SampleCount       int
}

// AuthenticityStats contains statistics about authenticity processing
type AuthenticityStats struct {
	ProcessedVariants    int64
	AuthenticityScores   []float64
	MeanAuthenticityScore float64
	MedianScore          float64
	StandardDeviation    float64
	
	// K-mer analysis stats
	TotalKmers           int64
	UniqueKmers          int64
	KmerComplexity       float64
	
	// Performance stats
	GPUProcessingTime    int64  // nanoseconds
	CPUProcessingTime    int64  // nanoseconds
	MemoryUsage          int64  // bytes
	
	// Population analysis
	PopulationProfiles   map[string]*PopulationStats
}

// PopulationStats contains population-specific authenticity statistics
type PopulationStats struct {
	PopulationID         string
	VariantCount         int64
	MeanAuthenticity     float64
	AuthenticityVariance float64
	KmerDiversity        float64
	FrequencySpectrum    map[string]int
}

// NewAuthenticityProcessor creates a new authenticity processor for VCF variants
func NewAuthenticityProcessor(config AuthenticityConfig) (*AuthenticityProcessor, error) {
	// Set defaults
	if config.KmerSize <= 0 {
		config.KmerSize = 8
	}
	if config.ContextWindowSize <= 0 {
		config.ContextWindowSize = 1000
	}

	// Create k-mer extractor with packed arrays
	kmerExtractor := kmer.NewPackedFrequencyExtractor(config.KmerSize)
	
	// Initialize GPU processor if requested
	var gpuProcessor *kmer.GPUPackedProcessor
	var err error
	if config.UseGPU {
		gpuProcessor, err = kmer.NewGPUPackedProcessor(config.KmerSize, kmer.GPUBackendAuto)
		if err != nil {
			fmt.Printf("Warning: GPU initialization failed: %v, falling back to CPU\n", err)
			config.UseGPU = false
		}
	}

	processor := &AuthenticityProcessor{
		kmerSize:            config.KmerSize,
		contextWindowSize:   config.ContextWindowSize,
		authenticityEnabled: true,
		kmerExtractor:       kmerExtractor,
		gpuProcessor:        gpuProcessor,
		useGPU:              config.UseGPU,
		populationContexts:  make(map[string]*PopulationContext),
		processingStats: AuthenticityStats{
			PopulationProfiles: make(map[string]*PopulationStats),
		},
	}

	return processor, nil
}

// AuthenticityConfig provides configuration for authenticity processing
type AuthenticityConfig struct {
	KmerSize          int
	ContextWindowSize int
	UseGPU            bool
	UseML             bool
	PopulationPanels  []string
	TargetScore       float64  // Reference value for the k-mer diversity score (50 = diversity 0.5); a heuristic, not a classifier accuracy
	EnableCaching     bool
	MaxCacheSize      int64
}

// ProcessVariant implements the VariantProcessor interface for authenticity analysis
func (ap *AuthenticityProcessor) ProcessVariant(variant *Variant) error {
	ap.mutex.Lock()
	defer ap.mutex.Unlock()

	// Extract k-mer profile from variant context
	kmerProfile, err := ap.extractVariantKmers(variant)
	if err != nil {
		return fmt.Errorf("failed to extract k-mers from variant: %w", err)
	}

	// Calculate authenticity score
	authenticityScore, err := ap.calculateAuthenticityScore(variant, kmerProfile)
	if err != nil {
		return fmt.Errorf("failed to calculate authenticity score: %w", err)
	}

	// Update processing statistics
	ap.processingStats.ProcessedVariants++
	ap.processingStats.AuthenticityScores = append(ap.processingStats.AuthenticityScores, authenticityScore)
	ap.processingStats.TotalKmers += int64(kmerProfile.TotalCount)
	ap.processingStats.UniqueKmers += int64(kmerProfile.NonZeroCount)

	// Update population-specific statistics
	if err := ap.updatePopulationStats(variant, kmerProfile, authenticityScore); err != nil {
		return fmt.Errorf("failed to update population statistics: %w", err)
	}

	// Store aggregated results
	ap.variantCount++
	ap.authenticityScores = append(ap.authenticityScores, authenticityScore)
	ap.kmerProfiles = append(ap.kmerProfiles, kmerProfile)

	return nil
}

// extractVariantKmers extracts k-mer profile from variant genomic context
func (ap *AuthenticityProcessor) extractVariantKmers(variant *Variant) (*kmer.PackedFrequencyTable, error) {
	if len(variant.Context) == 0 {
		return nil, fmt.Errorf("variant has no genomic context")
	}

	// Create sequence list for k-mer extraction
	sequences := [][]byte{variant.Context}

	// Add reference and alternate sequences
	if len(variant.Reference) > 0 {
		sequences = append(sequences, []byte(variant.Reference))
	}
	for _, alt := range variant.Alternate {
		if len(alt) > 0 {
			sequences = append(sequences, []byte(alt))
		}
	}

	// Extract k-mers using GPU acceleration if available
	if ap.useGPU && ap.gpuProcessor != nil {
		return ap.gpuProcessor.ExtractKmerFrequenciesGPUPacked(sequences)
	}

	// Fallback to CPU extraction with packed arrays
	return ap.kmerExtractor.ExtractFromMultipleSequences(sequences)
}

// calculateAuthenticityScore calculates the authenticity score for a variant
func (ap *AuthenticityProcessor) calculateAuthenticityScore(variant *Variant, kmerProfile *kmer.PackedFrequencyTable) (float64, error) {
	// Calculate a basic authenticity score based on k-mer complexity
	if kmerProfile.TotalCount == 0 {
		return 0.0, nil
	}
	// Use k-mer diversity as a proxy for authenticity score
	diversity := float64(kmerProfile.NonZeroCount) / float64(kmerProfile.MaxKmers)
	// Normalize to a score around 50% for typical sequences
	score := 50.0 + (diversity-0.5)*10.0
	return score, nil
}

// updatePopulationStats updates population-specific statistics
func (ap *AuthenticityProcessor) updatePopulationStats(variant *Variant, kmerProfile *kmer.PackedFrequencyTable, authenticityScore float64) error {
	// Extract population information from variant (would be more sophisticated in production)
	populationID := ap.extractPopulationID(variant)
	
	// Get or create population context
	popContext, exists := ap.populationContexts[populationID]
	if !exists {
		popContext = &PopulationContext{
			PopulationID:     populationID,
			VariantFrequency: make(map[string]float64),
			SampleCount:      0,
		}
		ap.populationContexts[populationID] = popContext
	}

	// Update population statistics
	popContext.SampleCount++
	variantKey := fmt.Sprintf("%s:%d:%s>%s", variant.Chromosome, variant.Position, 
		variant.Reference, strings.Join(variant.Alternate, ","))
	popContext.VariantFrequency[variantKey]++

	// Update population authenticity score (running average)
	if popContext.SampleCount == 1 {
		popContext.AuthenticityScore = authenticityScore
	} else {
		// Update running average
		alpha := 1.0 / float64(popContext.SampleCount)
		popContext.AuthenticityScore = (1-alpha)*popContext.AuthenticityScore + alpha*authenticityScore
	}

	// Merge k-mer profiles
	if popContext.KmerProfile == nil {
		popContext.KmerProfile = kmerProfile
	} else {
		// Merge k-mer profiles (simplified implementation)
		for i := 0; i < kmerProfile.MaxKmers; i++ {
			freq := kmerProfile.PackedFrequencies.Get(i)
			if freq > 0 {
				existingFreq := popContext.KmerProfile.PackedFrequencies.Get(i)
				popContext.KmerProfile.PackedFrequencies.Set(i, existingFreq+freq)
			}
		}
		popContext.KmerProfile.TotalCount += kmerProfile.TotalCount
	}

	// Update population-level processing statistics
	popStats, exists := ap.processingStats.PopulationProfiles[populationID]
	if !exists {
		popStats = &PopulationStats{
			PopulationID:      populationID,
			FrequencySpectrum: make(map[string]int),
		}
		ap.processingStats.PopulationProfiles[populationID] = popStats
	}

	popStats.VariantCount++
	// Update running statistics
	if popStats.VariantCount == 1 {
		popStats.MeanAuthenticity = authenticityScore
		popStats.AuthenticityVariance = 0.0
	} else {
		// Update running mean and variance
		oldMean := popStats.MeanAuthenticity
		popStats.MeanAuthenticity += (authenticityScore - oldMean) / float64(popStats.VariantCount)
		popStats.AuthenticityVariance += (authenticityScore - oldMean) * (authenticityScore - popStats.MeanAuthenticity)
	}

	popStats.KmerDiversity = float64(kmerProfile.NonZeroCount) / float64(kmerProfile.MaxKmers)

	return nil
}

// extractPopulationID extracts population ID from variant (placeholder implementation)
func (ap *AuthenticityProcessor) extractPopulationID(variant *Variant) string {
	// In a real implementation, this would extract population information
	// from variant INFO fields, sample names, or external annotation
	
	// Look for population in INFO field
	if popID, exists := variant.Info["POP"]; exists {
		return popID
	}
	
	// Default to "UNKNOWN" population
	return "UNKNOWN"
}

// GetAuthenticityResults returns comprehensive authenticity analysis results
func (ap *AuthenticityProcessor) GetAuthenticityResults() *AuthenticityResults {
	ap.mutex.RLock()
	defer ap.mutex.RUnlock()

	// Calculate summary statistics
	meanScore, medianScore, stdDev := ap.calculateStatistics(ap.authenticityScores)
	
	// Calculate k-mer complexity
	totalKmers := int64(0)
	uniqueKmers := int64(0)
	for _, profile := range ap.kmerProfiles {
		totalKmers += int64(profile.TotalCount)
		uniqueKmers += int64(profile.NonZeroCount)
	}
	
	complexity := 0.0
	if totalKmers > 0 {
		complexity = float64(uniqueKmers) / float64(totalKmers)
	}

	// Finalize population statistics variance calculation
	for _, popStats := range ap.processingStats.PopulationProfiles {
		if popStats.VariantCount > 1 {
			popStats.AuthenticityVariance /= float64(popStats.VariantCount - 1)
		}
	}

	results := &AuthenticityResults{
		ProcessedVariants:    ap.variantCount,
		MeanAuthenticityScore: meanScore,
		MedianScore:          medianScore,
		StandardDeviation:    stdDev,
		AuthenticityScores:   make([]float64, len(ap.authenticityScores)),
		
		// K-mer analysis
		TotalKmers:      totalKmers,
		UniqueKmers:     uniqueKmers,
		KmerComplexity:  complexity,
		
		// Population analysis
		PopulationResults: make(map[string]*PopulationAuthenticityResult),
		
		// Memory efficiency
		MemoryUsage:      ap.calculateMemoryUsage(),
		CompressionRatio: ap.calculateCompressionRatio(),
	}

	// Copy authenticity scores
	copy(results.AuthenticityScores, ap.authenticityScores)

	// Generate population results
	for popID, popContext := range ap.populationContexts {
		popStats := ap.processingStats.PopulationProfiles[popID]
		results.PopulationResults[popID] = &PopulationAuthenticityResult{
			PopulationID:         popID,
			SampleCount:          popContext.SampleCount,
			MeanAuthenticityScore: popContext.AuthenticityScore,
			VariantCount:         popStats.VariantCount,
			AuthenticityVariance: popStats.AuthenticityVariance,
			KmerDiversity:        popStats.KmerDiversity,
			MemoryUsage:          popContext.KmerProfile.GetMemoryUsage(),
		}
	}

	return results
}

// calculateStatistics calculates mean, median, and standard deviation
func (ap *AuthenticityProcessor) calculateStatistics(scores []float64) (mean, median, stdDev float64) {
	if len(scores) == 0 {
		return 0, 0, 0
	}

	// Calculate mean
	sum := 0.0
	for _, score := range scores {
		sum += score
	}
	mean = sum / float64(len(scores))

	// Calculate median (requires sorting)
	sortedScores := make([]float64, len(scores))
	copy(sortedScores, scores)
	
	// Simple bubble sort for median calculation
	for i := 0; i < len(sortedScores)-1; i++ {
		for j := i + 1; j < len(sortedScores); j++ {
			if sortedScores[i] > sortedScores[j] {
				sortedScores[i], sortedScores[j] = sortedScores[j], sortedScores[i]
			}
		}
	}
	
	if len(sortedScores)%2 == 0 {
		median = (sortedScores[len(sortedScores)/2-1] + sortedScores[len(sortedScores)/2]) / 2
	} else {
		median = sortedScores[len(sortedScores)/2]
	}

	// Calculate sample standard deviation (Bessel's correction: ÷(N-1)).
	if len(scores) > 1 {
		sumSq := 0.0
		for _, score := range scores {
			d := score - mean
			sumSq += d * d
		}
		stdDev = math.Sqrt(sumSq / float64(len(scores)-1))
	}

	return mean, median, stdDev
}

// calculateMemoryUsage calculates total memory usage of k-mer profiles
func (ap *AuthenticityProcessor) calculateMemoryUsage() int64 {
	totalMemory := int64(0)
	for _, profile := range ap.kmerProfiles {
		totalMemory += profile.GetMemoryUsage()
	}
	return totalMemory
}

// calculateCompressionRatio calculates memory compression ratio vs hash maps
func (ap *AuthenticityProcessor) calculateCompressionRatio() float64 {
	if len(ap.kmerProfiles) == 0 {
		return 1.0
	}

	totalCompression := 0.0
	for _, profile := range ap.kmerProfiles {
		totalCompression += profile.GetCompressionRatio()
	}
	
	return totalCompression / float64(len(ap.kmerProfiles))
}

// AuthenticityResults contains comprehensive authenticity analysis results
type AuthenticityResults struct {
	ProcessedVariants     int64
	MeanAuthenticityScore float64
	MedianScore          float64
	StandardDeviation    float64
	AuthenticityScores   []float64
	
	// K-mer analysis results
	TotalKmers           int64
	UniqueKmers          int64
	KmerComplexity       float64
	
	// Population-specific results
	PopulationResults    map[string]*PopulationAuthenticityResult
	
	// Memory efficiency metrics
	MemoryUsage          int64
	CompressionRatio     float64
	
	// Reference-range check set by IsAuthenticityTargetMet. Despite the
	// field name, this only records whether the mean k-mer diversity score
	// lies within tolerance of the reference value; it is not a statistical
	// test against real data. The name is kept for API compatibility.
	IsStatisticallyIndistinguishable bool
	TargetScoreDeviation            float64
}

// PopulationAuthenticityResult contains population-specific authenticity results
type PopulationAuthenticityResult struct {
	PopulationID         string
	SampleCount          int
	VariantCount         int64
	MeanAuthenticityScore float64
	AuthenticityVariance float64
	KmerDiversity        float64
	MemoryUsage          int64
}

// IsAuthenticityTargetMet reports whether the mean k-mer diversity score lies
// within tolerance of targetScore, and records the deviation. It is a range
// check on a heuristic score, not a test of statistical indistinguishability.
func (ar *AuthenticityResults) IsAuthenticityTargetMet(targetScore float64, tolerance float64) bool {
	deviation := ar.MeanAuthenticityScore - targetScore
	if deviation < 0 {
		deviation = -deviation
	}
	ar.TargetScoreDeviation = deviation
	ar.IsStatisticallyIndistinguishable = deviation <= tolerance
	return ar.IsStatisticallyIndistinguishable
}

// Close releases resources used by the authenticity processor
func (ap *AuthenticityProcessor) Close() error {
	if ap.gpuProcessor != nil {
		return ap.gpuProcessor.Close()
	}
	return nil
}