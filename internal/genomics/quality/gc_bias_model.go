// Package quality has two heuristic per-platform models: GC-content effects
// on base qualities (GCbiasQualityModel) and homopolymer-length-dependent
// indel and substitution errors (HomopolymerErrorModel).
//
// The per-platform parameters are hand-set, not fitted to published error
// profiles, unless a GC model is trained with TrainFromRealData. Both models
// seed their random source from the clock until SetSeed is called. The FASTQ
// generation pipeline does not use this package: ONT errors come from
// internal/genomics/nanopore/errmodel. training.QualityTrainer uses the GC
// model.
package quality

import (
	"fmt"
	"math"
	"math/rand/v2"
	"strings"
	"time"
)

// GCbiasQualityModel adjusts base qualities by local GC content, plus
// hand-set amplification, hexamer-priming and position-periodic terms. Its
// parameters are heuristics per platform family.
type GCbiasQualityModel struct {
	gcContentWindows  []float64           // GC content window values
	qualityAdjustment map[float64]float64 // GC content -> quality adjustment
	platform          string              // Sequencing platform
	amplificationBias bool                // Whether to model amplification bias
	windowSize        int                 // Window size for GC content calculation
	baseQuality       float64             // Base quality score
	biasStrength      float64             // Strength of GC bias effect
	rng               *rand.Rand          // Random number generator

	// Platform-specific parameters
	optimalGCContent float64 // Optimal GC content for platform
	biasSlope        float64 // Slope of bias curve
	saturationType   string  // "linear", "exponential" (actually quadratic), "sigmoid"

	// Training data
	trainingData []GCQualityDataPoint
	modelTrained bool
}

// GCQualityDataPoint represents a training data point
type GCQualityDataPoint struct {
	GCContent      float64
	QualityScore   float64
	SequenceLength int
	Platform       string
	Cycle          int
	ReadPosition   int
}

// QualityPredictionResult represents the result of quality prediction
type QualityPredictionResult struct {
	OriginalQualities []byte
	AdjustedQualities []byte
	GCContentProfile  []float64
	QualityProfile    []float64
	BiasCorrection    []float64
	PlatformEffects   []PlatformEffect
}

// PlatformEffect represents platform-specific quality effects
type PlatformEffect struct {
	Position    int
	EffectType  string // "gc_bias", "amplification", "priming", "cluster_density"
	Magnitude   float64
	Description string
}

// NewGCbiasQualityModel creates a new GC-bias quality model
func NewGCbiasQualityModel(platform string) *GCbiasQualityModel {
	model := &GCbiasQualityModel{
		gcContentWindows:  make([]float64, 0),
		qualityAdjustment: make(map[float64]float64),
		platform:          platform,
		amplificationBias: true,
		windowSize:        50,
		baseQuality:       30.0,
		biasStrength:      1.0,
		rng:               rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
		trainingData:      make([]GCQualityDataPoint, 0),
		modelTrained:      false,
	}

	// Set platform-specific parameters
	model.setPlatformParameters()

	return model
}

// setPlatformParameters sets platform-specific parameters
func (m *GCbiasQualityModel) setPlatformParameters() {
	switch m.platform {
	case "illumina", "novaseq", "nextseq", "miseq":
		m.optimalGCContent = 0.45 // Illumina platforms have GC bias
		m.biasSlope = 2.0
		m.saturationType = "exponential"
		m.amplificationBias = true
		m.windowSize = 50

	case "pacbio", "pacbio_revio":
		m.optimalGCContent = 0.50 // PacBio has less GC bias
		m.biasSlope = 0.5
		m.saturationType = "linear"
		m.amplificationBias = false
		m.windowSize = 100

	case "nanopore", "promethion":
		m.optimalGCContent = 0.50 // Nanopore has minimal GC bias
		m.biasSlope = 0.3
		m.saturationType = "linear"
		m.amplificationBias = false
		m.windowSize = 200

	default:
		// Default to Illumina-like behavior
		m.optimalGCContent = 0.45
		m.biasSlope = 2.0
		m.saturationType = "exponential"
		m.amplificationBias = true
		m.windowSize = 50
	}
}

// PredictQuality adjusts baseQualities (raw Phred values, not ASCII) by the
// GC term plus Gaussian noise (σ = 0.5) and the platform terms, clamped to
// [0, 41]. The per-bin adjustments learned by TrainFromRealData are not used;
// training only changes the optimal GC content and slope.
func (m *GCbiasQualityModel) PredictQuality(sequence []byte, baseQualities []byte) (*QualityPredictionResult, error) {
	if len(sequence) == 0 {
		return nil, fmt.Errorf("sequence cannot be empty")
	}

	if len(baseQualities) != len(sequence) {
		return nil, fmt.Errorf("sequence and quality lengths must match")
	}

	result := &QualityPredictionResult{
		OriginalQualities: make([]byte, len(baseQualities)),
		AdjustedQualities: make([]byte, len(baseQualities)),
		GCContentProfile:  make([]float64, len(sequence)),
		QualityProfile:    make([]float64, len(sequence)),
		BiasCorrection:    make([]float64, len(sequence)),
		PlatformEffects:   make([]PlatformEffect, 0),
	}

	copy(result.OriginalQualities, baseQualities)
	copy(result.AdjustedQualities, baseQualities)

	// Calculate GC content profile
	for i := 0; i < len(sequence); i++ {
		gcContent := m.calculateLocalGCContent(sequence, i)
		result.GCContentProfile[i] = gcContent
	}

	// Apply GC bias corrections
	for i := 0; i < len(sequence); i++ {
		gcContent := result.GCContentProfile[i]

		// Calculate quality adjustment
		qualityAdjustment := m.calculateGCBiasAdjustment(gcContent)
		result.BiasCorrection[i] = qualityAdjustment

		// Apply adjustment
		originalQuality := float64(baseQualities[i])
		adjustedQuality := originalQuality + qualityAdjustment

		// Apply platform-specific effects
		platformAdjustment := m.calculatePlatformEffects(sequence, i)
		adjustedQuality += platformAdjustment

		// Clamp to valid range
		if adjustedQuality < 0 {
			adjustedQuality = 0
		}
		if adjustedQuality > 41 {
			adjustedQuality = 41
		}

		result.AdjustedQualities[i] = byte(adjustedQuality)
		result.QualityProfile[i] = adjustedQuality

		// Record significant platform effects
		if math.Abs(platformAdjustment) > 2.0 {
			effect := PlatformEffect{
				Position:    i,
				EffectType:  m.getEffectType(sequence, i),
				Magnitude:   platformAdjustment,
				Description: m.getEffectDescription(sequence, i),
			}
			result.PlatformEffects = append(result.PlatformEffects, effect)
		}
	}

	return result, nil
}

// calculateLocalGCContent calculates GC content in a local window
func (m *GCbiasQualityModel) calculateLocalGCContent(sequence []byte, position int) float64 {
	// Define window boundaries
	windowStart := position - m.windowSize/2
	windowEnd := position + m.windowSize/2

	// Adjust boundaries to stay within sequence
	if windowStart < 0 {
		windowStart = 0
	}
	if windowEnd > len(sequence) {
		windowEnd = len(sequence)
	}

	// Count GC bases in window
	gcCount := 0
	totalCount := 0

	for i := windowStart; i < windowEnd; i++ {
		base := sequence[i]
		if base == 'G' || base == 'C' {
			gcCount++
		}
		if base == 'A' || base == 'T' || base == 'G' || base == 'C' {
			totalCount++
		}
	}

	if totalCount == 0 {
		return 0.5 // Default to 50% GC
	}

	return float64(gcCount) / float64(totalCount)
}

// calculateGCBiasAdjustment calculates quality adjustment due to GC bias
func (m *GCbiasQualityModel) calculateGCBiasAdjustment(gcContent float64) float64 {
	// Calculate distance from optimal GC content
	gcDeviation := gcContent - m.optimalGCContent

	// Apply bias function based on saturation type
	var biasEffect float64
	switch m.saturationType {
	case "linear":
		biasEffect = -m.biasSlope * math.Abs(gcDeviation)

	case "exponential":
		biasEffect = -m.biasSlope * math.Pow(math.Abs(gcDeviation), 2)

	case "sigmoid":
		// Sigmoid function centered at optimal GC content
		x := (gcContent - m.optimalGCContent) * 10 // Scale factor
		biasEffect = -m.biasSlope * (2.0/(1.0+math.Exp(-x)) - 1.0)

	default:
		biasEffect = -m.biasSlope * math.Abs(gcDeviation)
	}

	// Apply bias strength scaling
	biasEffect *= m.biasStrength

	// Add Gaussian noise (σ = 0.5 quality units)
	noise := m.rng.NormFloat64() * 0.5
	biasEffect += noise

	return biasEffect
}

// calculatePlatformEffects calculates platform-specific quality effects
func (m *GCbiasQualityModel) calculatePlatformEffects(sequence []byte, position int) float64 {
	totalEffect := 0.0

	// Amplification bias (for PCR-based platforms)
	if m.amplificationBias {
		ampEffect := m.calculateAmplificationBias(sequence, position)
		totalEffect += ampEffect
	}

	// Priming bias
	primingEffect := m.calculatePrimingBias(sequence, position)
	totalEffect += primingEffect

	// Cluster density effects (Illumina-specific)
	if m.platform == "illumina" || m.platform == "novaseq" || m.platform == "nextseq" {
		clusterEffect := m.calculateClusterDensityEffect(sequence, position)
		totalEffect += clusterEffect
	}

	return totalEffect
}

// calculateAmplificationBias calculates amplification bias effects
func (m *GCbiasQualityModel) calculateAmplificationBias(sequence []byte, position int) float64 {
	// Amplification bias affects GC-rich regions more
	gcContent := m.calculateLocalGCContent(sequence, position)

	// Model bias as function of GC content
	if gcContent > 0.6 {
		// High GC content reduces amplification efficiency
		return -2.0 * (gcContent - 0.6)
	} else if gcContent < 0.3 {
		// Low GC content can also reduce efficiency
		return -1.0 * (0.3 - gcContent)
	}

	return 0.0
}

// calculatePrimingBias applies fixed quality offsets after six hexamers
// (homopolymers and two others). The list and offsets are illustrative, not
// a measured random-hexamer priming profile.
func (m *GCbiasQualityModel) calculatePrimingBias(sequence []byte, position int) float64 {
	if position >= 6 {
		hexamer := string(sequence[position-6 : position])

		// Hand-picked hexamers and offsets
		favoredHexamers := map[string]float64{
			"AAAAAA": -1.0, // A-rich hexamers are disfavored
			"TTTTTT": -1.0,
			"GGGGGG": -2.0, // G-rich hexamers are strongly disfavored
			"CCCCCC": -2.0,
			"GCGCGC": 1.0, // Balanced hexamers are favored
			"ATCGAT": 1.0,
		}

		if bias, exists := favoredHexamers[hexamer]; exists {
			return bias
		}
	}

	return 0.0
}

// calculateClusterDensityEffect is a stand-in for Illumina cluster-density
// effects: a sinusoid of the position within the read. It does not model
// flow-cell cluster density.
func (m *GCbiasQualityModel) calculateClusterDensityEffect(sequence []byte, position int) float64 {
	// Position-based sinusoidal variation
	densityVariation := 0.5 * math.Sin(float64(position)/100.0)

	// High density regions have reduced quality
	if densityVariation > 0.3 {
		return -1.0 * densityVariation
	}

	return 0.0
}

// getEffectType determines the type of platform effect
func (m *GCbiasQualityModel) getEffectType(sequence []byte, position int) string {
	gcContent := m.calculateLocalGCContent(sequence, position)

	if gcContent > 0.6 || gcContent < 0.3 {
		return "gc_bias"
	}

	if m.amplificationBias && (gcContent > 0.65 || gcContent < 0.25) {
		return "amplification"
	}

	if position >= 6 {
		hexamer := string(sequence[position-6 : position])
		if strings.Contains(hexamer, "AAAA") || strings.Contains(hexamer, "GGGG") {
			return "priming"
		}
	}

	return "cluster_density"
}

// getEffectDescription provides a description of the platform effect
func (m *GCbiasQualityModel) getEffectDescription(sequence []byte, position int) string {
	gcContent := m.calculateLocalGCContent(sequence, position)

	if gcContent > 0.6 {
		return "High GC content reduces polymerase efficiency"
	} else if gcContent < 0.3 {
		return "Low GC content affects secondary structure"
	}

	return "Platform-specific sequencing artifact"
}

// TrainFromRealData bins data points by GC content (0.01 bins), sets the
// optimal GC content to the bin with the highest mean quality, and sets the
// slope to Σ(quality drop) / Σ(GC deviation²) over bins. The slope is a
// ratio estimate, not a least-squares fit.
func (m *GCbiasQualityModel) TrainFromRealData(dataPoints []GCQualityDataPoint) error {
	if len(dataPoints) == 0 {
		return fmt.Errorf("no training data provided")
	}

	m.trainingData = dataPoints

	// Build GC content -> quality adjustment mapping
	gcBins := make(map[float64][]float64)

	for _, point := range dataPoints {
		// Round GC content to nearest 0.01
		gcBin := math.Round(point.GCContent*100) / 100
		gcBins[gcBin] = append(gcBins[gcBin], point.QualityScore)
	}

	// Calculate average quality for each GC bin
	for gcContent, qualities := range gcBins {
		avgQuality := 0.0
		for _, quality := range qualities {
			avgQuality += quality
		}
		avgQuality /= float64(len(qualities))

		// Calculate adjustment relative to base quality
		adjustment := avgQuality - m.baseQuality
		m.qualityAdjustment[gcContent] = adjustment
	}

	// Fit bias parameters
	m.fitBiasParameters()

	m.modelTrained = true
	return nil
}

// fitBiasParameters fits bias parameters to training data
func (m *GCbiasQualityModel) fitBiasParameters() {
	if len(m.trainingData) == 0 {
		return
	}

	// Find optimal GC content (GC content with highest average quality)
	gcQualityMap := make(map[float64][]float64)

	for _, point := range m.trainingData {
		gcBin := math.Round(point.GCContent*100) / 100
		gcQualityMap[gcBin] = append(gcQualityMap[gcBin], point.QualityScore)
	}

	maxAvgQuality := 0.0
	optimalGC := 0.5

	for gcContent, qualities := range gcQualityMap {
		avgQuality := 0.0
		for _, quality := range qualities {
			avgQuality += quality
		}
		avgQuality /= float64(len(qualities))

		if avgQuality > maxAvgQuality {
			maxAvgQuality = avgQuality
			optimalGC = gcContent
		}
	}

	m.optimalGCContent = optimalGC

	// Estimate bias slope
	sumSquaredDeviations := 0.0
	sumQualityDeviations := 0.0
	count := 0

	for gcContent, qualities := range gcQualityMap {
		avgQuality := 0.0
		for _, quality := range qualities {
			avgQuality += quality
		}
		avgQuality /= float64(len(qualities))

		gcDeviation := math.Abs(gcContent - m.optimalGCContent)
		qualityDeviation := maxAvgQuality - avgQuality

		sumSquaredDeviations += gcDeviation * gcDeviation
		sumQualityDeviations += qualityDeviation
		count++
	}

	if sumSquaredDeviations > 0 {
		m.biasSlope = sumQualityDeviations / sumSquaredDeviations
	}
}

// SetBiasStrength sets the strength of GC bias effects
func (m *GCbiasQualityModel) SetBiasStrength(strength float64) {
	if strength >= 0 {
		m.biasStrength = strength
	}
}

// SetWindowSize sets the window size for GC content calculation
func (m *GCbiasQualityModel) SetWindowSize(size int) {
	if size > 0 {
		m.windowSize = size
	}
}

// ValidateModel validates the trained model
func (m *GCbiasQualityModel) ValidateModel() error {
	if !m.modelTrained {
		return fmt.Errorf("model has not been trained")
	}

	if len(m.qualityAdjustment) == 0 {
		return fmt.Errorf("no quality adjustments available")
	}

	if m.optimalGCContent < 0 || m.optimalGCContent > 1 {
		return fmt.Errorf("invalid optimal GC content: %f", m.optimalGCContent)
	}

	if m.biasSlope < 0 {
		return fmt.Errorf("invalid bias slope: %f", m.biasSlope)
	}

	return nil
}

// GetModelStats returns statistics about the trained model
func (m *GCbiasQualityModel) GetModelStats() map[string]interface{} {
	stats := map[string]interface{}{
		"platform":             m.platform,
		"model_trained":        m.modelTrained,
		"optimal_gc_content":   m.optimalGCContent,
		"bias_slope":           m.biasSlope,
		"saturation_type":      m.saturationType,
		"bias_strength":        m.biasStrength,
		"window_size":          m.windowSize,
		"amplification_bias":   m.amplificationBias,
		"training_data_points": len(m.trainingData),
		"gc_bins":              len(m.qualityAdjustment),
	}

	if len(m.qualityAdjustment) > 0 {
		maxAdjustment := 0.0
		minAdjustment := 0.0

		for _, adjustment := range m.qualityAdjustment {
			if adjustment > maxAdjustment {
				maxAdjustment = adjustment
			}
			if adjustment < minAdjustment {
				minAdjustment = adjustment
			}
		}

		stats["max_quality_adjustment"] = maxAdjustment
		stats["min_quality_adjustment"] = minAdjustment
		stats["adjustment_range"] = maxAdjustment - minAdjustment
	}

	return stats
}

// SetSeed sets the random seed for reproducible quality prediction
func (m *GCbiasQualityModel) SetSeed(seed int64) {
	m.rng = rand.New(rand.NewPCG(uint64(seed), 0)) // #nosec G115 G404
}

// GetPlatformParameters returns platform-specific parameters
func (m *GCbiasQualityModel) GetPlatformParameters() map[string]interface{} {
	return map[string]interface{}{
		"platform":           m.platform,
		"optimal_gc_content": m.optimalGCContent,
		"bias_slope":         m.biasSlope,
		"saturation_type":    m.saturationType,
		"amplification_bias": m.amplificationBias,
		"window_size":        m.windowSize,
	}
}

// SimulateQualityProfile simulates a quality profile for a sequence
func (m *GCbiasQualityModel) SimulateQualityProfile(sequence []byte, targetQuality float64) ([]byte, error) {
	if len(sequence) == 0 {
		return nil, fmt.Errorf("sequence cannot be empty")
	}

	// Generate base qualities
	baseQualities := make([]byte, len(sequence))
	for i := range baseQualities {
		baseQualities[i] = byte(targetQuality)
	}

	// Apply GC bias effects
	result, err := m.PredictQuality(sequence, baseQualities)
	if err != nil {
		return nil, err
	}

	return result.AdjustedQualities, nil
}

// Additional utility functions

// calculateSequenceComplexity calculates sequence complexity
func (m *GCbiasQualityModel) calculateSequenceComplexity(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0.0
	}

	// Count unique k-mers (k=3)
	kmerCount := make(map[string]int)
	for i := 0; i <= len(sequence)-3; i++ {
		kmer := string(sequence[i : i+3])
		kmerCount[kmer]++
	}

	// Calculate complexity as fraction of unique k-mers
	maxPossibleKmers := len(sequence) - 2
	if maxPossibleKmers <= 0 {
		return 0.0
	}

	return float64(len(kmerCount)) / float64(maxPossibleKmers)
}

// calculateEntropy calculates sequence entropy
func (m *GCbiasQualityModel) calculateEntropy(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0.0
	}

	// Count base frequencies
	baseCounts := make(map[byte]int)
	for _, base := range sequence {
		baseCounts[base]++
	}

	// Calculate entropy
	entropy := 0.0
	total := float64(len(sequence))

	for _, count := range baseCounts {
		if count > 0 {
			p := float64(count) / total
			entropy -= p * math.Log2(p)
		}
	}

	return entropy
}

// Helper function to calculate string index
func stringIndex(s, substr string) int {
	return strings.Index(s, substr)
}
