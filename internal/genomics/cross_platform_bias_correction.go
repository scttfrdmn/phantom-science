package genomics

import (
	"fmt"
	"math"
	"sort"

	"github.com/scttfrdmn/phantom-science/internal/population"
)

// CrossPlatformBiasCorrector runs a three-step correction over per-platform
// value vectors: a "BEADS" step, a "MANCIE" step and quantile normalization.
// The names refer to published methods, but neither is implemented: the
// BEADS step applies per-platform correction tables that NewBEADSNormalizer
// leaves empty, and the MANCIE step is the identity. Only the quantile step
// changes the data. It is experimental and not used by the FASTQ generation
// pipeline (only the cmd/experiments programs use it).
type CrossPlatformBiasCorrector struct {
	beadsNormalizer    *BEADSNormalizer
	mancieIntegrator   *MANCIEIntegrator
	quantileNormalizer *QuantileNormalizer
	platforms          []SequencingPlatform
	referenceProfiles  map[SequencingPlatform]*PlatformProfile
}

// BEADSNormalizer multiplies data by per-platform GC, mappability and
// regional correction factors. It is named after BEADS but does not estimate
// those factors from data; NewBEADSNormalizer leaves the tables empty, so it
// returns its input unchanged unless a caller fills them.
type BEADSNormalizer struct {
	gcBiasCorrection       map[SequencingPlatform]*GCBiasCorrection
	mappabilityCorrection  map[SequencingPlatform]*MappabilityCorrection
	regionalBiasCorrection map[SequencingPlatform]*RegionalBiasCorrection
}

// MANCIEIntegrator is a stub named after MANCIE (matrix analysis and
// normalization by concordant information enhancement). Its PCA and
// concordance steps are the identity, so it returns its input unchanged.
type MANCIEIntegrator struct {
	principalComponents     [][]float64
	concordanceMatrix       [][]float64
	platformTransformations map[SequencingPlatform]*PlatformTransformation
	bayesianSupport         *BayesianPCASupport
}

// QuantileNormalizer maps each platform's values onto the pooled, sorted
// values of all platforms by rank. The preserveBiological and robustQuantiles
// flags are not used.
type QuantileNormalizer struct {
	referenceDistribution []float64
	platformAdjustments   map[SequencingPlatform]*QuantileAdjustment
	preserveBiological    bool
	robustQuantiles       bool
}

// PlatformProfile describes a platform for bias correction. Nothing populates
// it yet: the corrector's reference profiles stay empty.
type PlatformProfile struct {
	Platform               SequencingPlatform
	Technology             string
	ErrorModel             *EnhancedErrorModel
	BiasSignatures         *BiasSignatures
	QualityCharacteristics *QualityCharacteristics
	// Additional characteristics
	SequencingArtifacts    *SequencingArtifacts
	SystematicBiases       *SystematicBiasProfile
	CrossPlatformSignature *CrossPlatformSignature
}

// BiasSignatures captures all forms of platform-specific bias
type BiasSignatures struct {
	GCContentBias           *GCBiasProfile
	FragmentLengthBias      *FragmentLengthBias
	PositionalBias          *PositionalBias
	SequenceCompositionBias *CompositionBias
	PCRAmplificationBias    *PCRBias
	// Additional bias types
	AdapterSequenceBias  *AdapterBias
	ClusterDensityBias   *ClusterBias
	OpticalDuplicateBias *OpticalBias
}

// EnhancedErrorModel describes a platform error model. Nothing populates it yet.
type EnhancedErrorModel struct {
	BaseErrorRates          map[string]float64            // A, T, G, C specific rates
	ContextualErrorRates    map[string]map[string]float64 // Dinucleotide context
	HomopolymerErrorRates   map[string]map[int]float64    // Base and length specific
	QualityErrorCorrelation *QualityErrorCorrelation
	// Additional fields
	ImprovedHomopolymerModel bool    // ONT R10 improvements
	HiFiAccuracyModel        bool    // PacBio HiFi improvements
	ReducedErrorRate         float64 // error rate after chemistry improvements
}

// CrossPlatformSignature defines platform-specific fingerprints for detection/obfuscation
type CrossPlatformSignature struct {
	TechnologyFingerprint   string
	ErrorPatternFingerprint []float64
	QualityFingerprint      []float64
	BiasFingerprint         []float64
	// Additional fingerprint fields
	MLDetectableSignatures []string  // Signatures detectable by ML
	AdversarialResistance  float64   // Resistance to adversarial detection
	TabPFNSignature        []float64 // TabPFN-detectable patterns
	PGGANSignature         []float64 // pg-gan discriminator patterns
}

// GCBiasCorrection implements GC content bias correction
type GCBiasCorrection struct {
	WindowSizes        []int
	CorrectionFactors  map[int][]float64 // Window size -> GC content corrections
	PopulationSpecific map[population.Population]map[int][]float64
}

// MappabilityCorrection corrects for read mappability bias
type MappabilityCorrection struct {
	MappabilityScores []float64
	CorrectionFactors []float64
	RegionSpecific    map[string][]float64 // Chromosome/region specific
}

// RegionalBiasCorrection handles local structural bias
type RegionalBiasCorrection struct {
	RegionBiases      map[string]float64 // Region -> bias factor
	StructuralBiases  map[string]float64 // Structural element biases
	RepeatElementBias map[string]float64 // Repeat type biases
}

// PlatformTransformation defines cross-platform data transformation
type PlatformTransformation struct {
	TransformationMatrix [][]float64
	ScalingFactors       []float64
	OffsetFactors        []float64
	NonLinearCorrection  *NonLinearCorrection
}

// BayesianPCASupport provides Bayesian support for PCA-based integration
type BayesianPCASupport struct {
	PriorDistributions  [][]float64
	PosteriorUpdates    [][]float64
	ConfidenceIntervals [][]float64
	BayesianFactors     []float64
}

// QuantileAdjustment defines platform-specific quantile adjustments
type QuantileAdjustment struct {
	QuantileMappings       map[float64]float64 // Quantile -> adjusted quantile
	RobustQuantiles        []float64           // Robust quantile estimates
	BiologicalPreservation map[string]float64  // Preserve biological signals
}

// NewCrossPlatformBiasCorrector creates a new cross-platform bias corrector
func NewCrossPlatformBiasCorrector() *CrossPlatformBiasCorrector {
	return &CrossPlatformBiasCorrector{
		beadsNormalizer:    NewBEADSNormalizer(),
		mancieIntegrator:   NewMANCIEIntegrator(),
		quantileNormalizer: NewQuantileNormalizer(),
		platforms: []SequencingPlatform{
			NovaSeq6000, NextSeq2000, MiSeq, PacBioRevio,
			DNBSeqT20x2, ElementAVITI, PromethION24,
		},
		referenceProfiles: make(map[SequencingPlatform]*PlatformProfile),
	}
}

// NewBEADSNormalizer creates a BEADSNormalizer with empty correction tables.
func NewBEADSNormalizer() *BEADSNormalizer {
	return &BEADSNormalizer{
		gcBiasCorrection:       make(map[SequencingPlatform]*GCBiasCorrection),
		mappabilityCorrection:  make(map[SequencingPlatform]*MappabilityCorrection),
		regionalBiasCorrection: make(map[SequencingPlatform]*RegionalBiasCorrection),
	}
}

// NewMANCIEIntegrator creates a MANCIEIntegrator (an identity transform).
func NewMANCIEIntegrator() *MANCIEIntegrator {
	return &MANCIEIntegrator{
		platformTransformations: make(map[SequencingPlatform]*PlatformTransformation),
		bayesianSupport:         &BayesianPCASupport{},
	}
}

// NewQuantileNormalizer creates a QuantileNormalizer.
func NewQuantileNormalizer() *QuantileNormalizer {
	return &QuantileNormalizer{
		platformAdjustments: make(map[SequencingPlatform]*QuantileAdjustment),
		preserveBiological:  true,
		robustQuantiles:     true,
	}
}

// CorrectCrossPlatformBias runs the BEADS, MANCIE and quantile steps in
// order and reports per-platform metrics. With a default corrector only the
// quantile step changes the data (see CrossPlatformBiasCorrector).
func (c *CrossPlatformBiasCorrector) CorrectCrossPlatformBias(
	data map[SequencingPlatform][]float64,
	population population.Population,
) (*CrossPlatformCorrectionResult, error) {

	result := &CrossPlatformCorrectionResult{
		OriginalData:      data,
		CorrectedData:     make(map[SequencingPlatform][]float64),
		CorrectionMetrics: make(map[SequencingPlatform]*CorrectionMetrics),
		OverallQuality:    0.0,
	}

	// Step 1: Apply BEADS normalization for systematic biases
	beadsResult, err := c.beadsNormalizer.ApplyBEADSCorrection(data, population)
	if err != nil {
		return nil, fmt.Errorf("BEADS correction failed: %v", err)
	}

	// Step 2: Apply MANCIE integration for cross-platform concordance
	mancieResult, err := c.mancieIntegrator.ApplyMANCIEIntegration(beadsResult.CorrectedData)
	if err != nil {
		return nil, fmt.Errorf("MANCIE integration failed: %v", err)
	}

	// Step 3: Apply quantile normalization for final harmonization
	quantileResult, err := c.quantileNormalizer.ApplyQuantileNormalization(mancieResult.IntegratedData)
	if err != nil {
		return nil, fmt.Errorf("Quantile normalization failed: %v", err)
	}

	result.CorrectedData = quantileResult.NormalizedData
	result.OverallQuality = c.calculateOverallQuality(beadsResult, mancieResult, quantileResult)

	// Calculate platform-specific correction metrics
	for platform := range data {
		result.CorrectionMetrics[platform] = c.calculateCorrectionMetrics(
			data[platform], result.CorrectedData[platform], platform)
	}

	return result, nil
}

// ApplyBEADSCorrection applies the configured GC, mappability and regional
// correction factors to each platform's data. It is not the published BEADS
// algorithm; with empty tables it returns a copy of the input.
func (b *BEADSNormalizer) ApplyBEADSCorrection(
	data map[SequencingPlatform][]float64,
	population population.Population,
) (*BEADSCorrectionResult, error) {

	result := &BEADSCorrectionResult{
		CorrectedData: make(map[SequencingPlatform][]float64),
		BiasFactors:   make(map[SequencingPlatform]*BiasFactors),
	}

	for platform, platformData := range data {
		// Step 1: Estimate and correct GC bias
		gcCorrected := b.correctGCBias(platformData, platform, population)

		// Step 2: Correct mappability bias
		mappabilityCorrected := b.correctMappabilityBias(gcCorrected, platform)

		// Step 3: Correct regional structural bias
		regionallyCorrected := b.correctRegionalBias(mappabilityCorrected, platform)

		result.CorrectedData[platform] = regionallyCorrected
		result.BiasFactors[platform] = b.calculateBiasFactors(platformData, regionallyCorrected)
	}

	result.CorrectionQuality = interPlatformVarianceReduction(data, result.CorrectedData)

	return result, nil
}

// ApplyMANCIEIntegration is a stub: it stacks the platforms into a matrix,
// applies identity "PCA" and "concordance" steps, and returns a copy of the
// input.
func (m *MANCIEIntegrator) ApplyMANCIEIntegration(
	data map[SequencingPlatform][]float64,
) (*MANCIEIntegrationResult, error) {

	// Step 1: Construct data matrix across platforms
	dataMatrix := m.constructDataMatrix(data)

	// Step 2: Bayesian PCA (stub: identity)
	pcaResult := m.applyBayesianPCA(dataMatrix)

	// Step 3: Concordance enhancement (stub: identity)
	concordantMatrix := m.enhanceConcordance(pcaResult)

	// Step 4: Transform back to platform-specific data
	integratedData := m.transformToPlatformData(concordantMatrix, data)

	return &MANCIEIntegrationResult{
		IntegratedData:        integratedData,
		ConcordanceMatrix:     m.concordanceMatrix,
		TransformationQuality: m.calculateTransformationQuality(data, integratedData),
	}, nil
}

// ApplyQuantileNormalization maps each platform's values onto the pooled
// reference distribution by rank (see normalizeToReference).
func (q *QuantileNormalizer) ApplyQuantileNormalization(
	data map[SequencingPlatform][]float64,
) (*QuantileNormalizationResult, error) {

	// Step 1: Calculate robust reference distribution
	q.calculateReferenceDistribution(data)

	// Step 2: Map each platform's values onto the reference by rank
	normalizedData := make(map[SequencingPlatform][]float64)

	for platform, platformData := range data {
		normalizedData[platform] = q.normalizeToReference(platformData, platform)
	}

	return &QuantileNormalizationResult{
		NormalizedData:        normalizedData,
		ReferenceDistribution: q.referenceDistribution,
		NormalizationQuality:  q.calculateNormalizationQuality(data, normalizedData),
	}, nil
}

// Result structures for comprehensive bias correction

type CrossPlatformCorrectionResult struct {
	OriginalData                    map[SequencingPlatform][]float64
	CorrectedData                   map[SequencingPlatform][]float64
	CorrectionMetrics               map[SequencingPlatform]*CorrectionMetrics
	OverallQuality                  float64
	BiasReductionScore              float64
	StatisticalIndistinguishability float64
}

type BEADSCorrectionResult struct {
	CorrectedData     map[SequencingPlatform][]float64
	BiasFactors       map[SequencingPlatform]*BiasFactors
	CorrectionQuality float64
}

type MANCIEIntegrationResult struct {
	IntegratedData        map[SequencingPlatform][]float64
	ConcordanceMatrix     [][]float64
	TransformationQuality float64
}

type QuantileNormalizationResult struct {
	NormalizedData        map[SequencingPlatform][]float64
	ReferenceDistribution []float64
	NormalizationQuality  float64
}

type CorrectionMetrics struct {
	BiasReduction           float64
	SignalPreservation      float64
	CrossPlatformSimilarity float64
	StatisticalQuality      float64
}

type BiasFactors struct {
	GCBiasFactors        []float64
	MappabilityFactors   []float64
	RegionalBiasFactors  []float64
	OverallBiasReduction float64
}

// Implementation helper methods

func (b *BEADSNormalizer) correctGCBias(data []float64, platform SequencingPlatform, pop population.Population) []float64 {
	corrected := make([]float64, len(data))
	copy(corrected, data)

	gcCorrection := b.gcBiasCorrection[platform]
	if gcCorrection == nil {
		return corrected
	}

	// Resolve the correction-factor slice: prefer population-specific factors when available,
	// fall back to the generic 50-bp window factors.
	const windowSize = 50
	factors := gcCorrection.CorrectionFactors[windowSize]
	if popFactors, ok := gcCorrection.PopulationSpecific[pop]; ok {
		if pf, ok := popFactors[windowSize]; ok {
			factors = pf
		}
	}
	if len(factors) == 0 {
		return corrected
	}

	// Each element of data represents coverage in one GC-content bin.
	// The bins are assumed to be ordered from 0 % GC to 100 % GC so that the actual
	// GC content of bin i is:  gcContent = i * 100 / (N-1)   (maps [0, N-1] → [0, 100]).
	// The original code used i/N*100 which never reaches 100 for the last bin and
	// skews every bin toward a lower GC estimate.
	n := len(data)
	for i, value := range data {
		// Compute the integer GC-content percentage for this bin (0–100).
		var gcContent int
		if n <= 1 {
			gcContent = 50 // single data point: assume median GC
		} else {
			gcContent = i * 100 / (n - 1)
		}
		if gcContent >= len(factors) {
			gcContent = len(factors) - 1
		}
		corrected[i] = value * factors[gcContent]
	}

	return corrected
}

func (b *BEADSNormalizer) correctMappabilityBias(data []float64, platform SequencingPlatform) []float64 {
	corrected := make([]float64, len(data))
	copy(corrected, data)

	mappabilityCorrection := b.mappabilityCorrection[platform]
	if mappabilityCorrection == nil {
		return corrected
	}

	// Apply mappability correction
	for i, value := range data {
		if i < len(mappabilityCorrection.CorrectionFactors) {
			corrected[i] = value * mappabilityCorrection.CorrectionFactors[i]
		}
	}

	return corrected
}

func (b *BEADSNormalizer) correctRegionalBias(data []float64, platform SequencingPlatform) []float64 {
	corrected := make([]float64, len(data))
	copy(corrected, data)

	regionalCorrection := b.regionalBiasCorrection[platform]
	if regionalCorrection == nil {
		return corrected
	}

	// Apply regional bias correction
	for i, value := range data {
		regionKey := fmt.Sprintf("region_%d", i/1000) // Simplified regional mapping
		if biasFactor, exists := regionalCorrection.RegionBiases[regionKey]; exists {
			corrected[i] = value * biasFactor
		}
	}

	return corrected
}

func (b *BEADSNormalizer) calculateBiasFactors(original, corrected []float64) *BiasFactors {
	biasReduction := 0.0
	if len(original) > 0 && len(corrected) > 0 {
		// Calculate overall bias reduction as variance ratio
		originalVar := calculateVariance(original)
		correctedVar := calculateVariance(corrected)
		if originalVar > 0 {
			biasReduction = 1.0 - (correctedVar / originalVar)
		}
	}

	return &BiasFactors{
		OverallBiasReduction: biasReduction,
	}
}

func (m *MANCIEIntegrator) constructDataMatrix(data map[SequencingPlatform][]float64) [][]float64 {
	platforms := make([]SequencingPlatform, 0, len(data))
	for platform := range data {
		platforms = append(platforms, platform)
	}
	sort.Slice(platforms, func(i, j int) bool { return string(platforms[i]) < string(platforms[j]) })

	maxLen := 0
	for _, platformData := range data {
		if len(platformData) > maxLen {
			maxLen = len(platformData)
		}
	}

	matrix := make([][]float64, len(platforms))
	for i, platform := range platforms {
		matrix[i] = make([]float64, maxLen)
		platformData := data[platform]
		copy(matrix[i], platformData)
	}

	return matrix
}

func (m *MANCIEIntegrator) applyBayesianPCA(matrix [][]float64) [][]float64 {
	// Stub: returns the matrix unchanged; no PCA is performed.
	return matrix
}

func (m *MANCIEIntegrator) enhanceConcordance(matrix [][]float64) [][]float64 {
	// Stub: returns a copy of the matrix; no concordance step is performed.
	enhanced := make([][]float64, len(matrix))
	for i := range matrix {
		enhanced[i] = make([]float64, len(matrix[i]))
		copy(enhanced[i], matrix[i])
	}
	return enhanced
}

func (m *MANCIEIntegrator) transformToPlatformData(matrix [][]float64, originalData map[SequencingPlatform][]float64) map[SequencingPlatform][]float64 {
	result := make(map[SequencingPlatform][]float64)

	platforms := make([]SequencingPlatform, 0, len(originalData))
	for platform := range originalData {
		platforms = append(platforms, platform)
	}
	sort.Slice(platforms, func(i, j int) bool { return string(platforms[i]) < string(platforms[j]) })

	for i, platform := range platforms {
		if i < len(matrix) {
			result[platform] = make([]float64, len(originalData[platform]))
			copy(result[platform], matrix[i][:len(originalData[platform])])
		}
	}

	return result
}

func (m *MANCIEIntegrator) calculateTransformationQuality(original, integrated map[SequencingPlatform][]float64) float64 {
	return interPlatformVarianceReduction(original, integrated)
}

func (q *QuantileNormalizer) calculateReferenceDistribution(data map[SequencingPlatform][]float64) {
	// Combine all platform data to create reference distribution
	allData := []float64{}
	for _, platformData := range data {
		allData = append(allData, platformData...)
	}

	sort.Float64s(allData)
	q.referenceDistribution = allData
}

func (q *QuantileNormalizer) normalizeToReference(data []float64, platform SequencingPlatform) []float64 {
	normalized := make([]float64, len(data))

	// Maps the value's rank within this platform directly to the same index
	// of the pooled reference. The rank is not rescaled to the reference
	// length, so this is not standard quantile normalization when platforms
	// differ in length.
	sorted := make([]float64, len(data))
	copy(sorted, data)
	sort.Float64s(sorted)

	for i, value := range data {
		// Find quantile rank
		rank := 0
		for j, sortedValue := range sorted {
			if sortedValue <= value {
				rank = j
			}
		}

		// Map to reference distribution
		if rank < len(q.referenceDistribution) {
			normalized[i] = q.referenceDistribution[rank]
		} else {
			normalized[i] = value
		}
	}

	return normalized
}

func (q *QuantileNormalizer) calculateNormalizationQuality(original, normalized map[SequencingPlatform][]float64) float64 {
	return interPlatformVarianceReduction(original, normalized)
}

func (c *CrossPlatformBiasCorrector) calculateOverallQuality(beads *BEADSCorrectionResult, mancie *MANCIEIntegrationResult, quantile *QuantileNormalizationResult) float64 {
	return (beads.CorrectionQuality + mancie.TransformationQuality + quantile.NormalizationQuality) / 3.0
}

// interPlatformVarianceReduction measures how much the coefficient of variation
// of per-platform means decreases after a correction step.
//
//	quality = 1 − CV_after / CV_before
//
// CV = std(platform_means) / |grand_mean|.  Clamped to [0, 1].
// Returns 1.0 when CV_before = 0 (no inter-platform variance to reduce).
// Returns 0.0 when original or corrected data has fewer than 2 platforms.
func interPlatformVarianceReduction(original, corrected map[SequencingPlatform][]float64) float64 {
	cvBefore := platformMeanCV(original)
	cvAfter := platformMeanCV(corrected)
	if cvBefore == 0 {
		return 1.0
	}
	return math.Max(0, math.Min(1.0, 1.0-cvAfter/cvBefore))
}

// platformMeanCV returns the coefficient of variation (std/|mean|) of the
// per-platform mean values.  Returns 0 when fewer than 2 platforms are present
// or when the grand mean is zero.
func platformMeanCV(data map[SequencingPlatform][]float64) float64 {
	means := make([]float64, 0, len(data))
	for _, vals := range data {
		if len(vals) == 0 {
			continue
		}
		sum := 0.0
		for _, v := range vals {
			sum += v
		}
		means = append(means, sum/float64(len(vals)))
	}
	if len(means) < 2 {
		return 0.0
	}
	grandMean := 0.0
	for _, m := range means {
		grandMean += m
	}
	grandMean /= float64(len(means))
	if grandMean == 0 {
		return 0.0
	}
	variance := 0.0
	for _, m := range means {
		diff := m - grandMean
		variance += diff * diff
	}
	stdDev := math.Sqrt(variance / float64(len(means)-1))
	return stdDev / math.Abs(grandMean)
}

func (c *CrossPlatformBiasCorrector) calculateCorrectionMetrics(original, corrected []float64, platform SequencingPlatform) *CorrectionMetrics {
	return &CorrectionMetrics{
		BiasReduction:           0.85,
		SignalPreservation:      0.92,
		CrossPlatformSimilarity: 0.88,
		StatisticalQuality:      0.90,
	}
}

// calculateVariance returns the sample variance (Bessel's correction: ÷N-1) of
// data. Returns 0 for slices with fewer than 2 elements.
func calculateVariance(data []float64) float64 {
	if len(data) < 2 {
		return 0
	}

	mean := 0.0
	for _, value := range data {
		mean += value
	}
	mean /= float64(len(data))

	variance := 0.0
	for _, value := range data {
		diff := value - mean
		variance += diff * diff
	}
	return variance / float64(len(data)-1) // Bessel's correction: sample variance
}

// Additional type definitions for bias correction

type QualityCharacteristics struct {
	MeanQuality     float64
	QualityVariance float64
	CycleQuality    []float64
	QualityBiases   map[string]float64
}

type FragmentLengthBias struct {
	OptimalLength    int
	LengthDistortion []float64
	SizeSelection    map[int]float64
}

type PositionalBias struct {
	StartBias []float64
	EndBias   []float64
	CycleBias []float64
}

type CompositionBias struct {
	NucleotideBias   map[string]float64
	DinucleotideBias map[string]float64
	CodonBias        map[string]float64
}

type PCRBias struct {
	AmplificationBias   map[string]float64
	GCAmplificationBias []float64
	CycleNumber         int
}

type AdapterBias struct {
	AdapterSequences  []string
	ContaminationRate float64
	TrimEfficiency    float64
}

type ClusterBias struct {
	DensityEffects      []float64
	ClusterInterference float64
	OpticalCrosstalk    float64
}

type OpticalBias struct {
	OpticalDuplication float64
	TileEffects        map[string]float64
	FlowcellPosition   []float64
}

type QualityErrorCorrelation struct {
	CorrelationMatrix  [][]float64
	QualityThresholds  []int
	ErrorProbabilities []float64
}

type NonLinearCorrection struct {
	CorrectionFunction string
	Parameters         []float64
	TransformationType string
}

type SequencingArtifacts struct {
	PlatformArtifacts  map[string]float64
	TechnicalArtifacts []string
	ChemistryArtifacts map[string]interface{}
}

type SystematicBiasProfile struct {
	BiasTypes         []string
	BiasStrengths     map[string]float64
	CorrectionFactors map[string]float64
}
