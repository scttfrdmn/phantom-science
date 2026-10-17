package authenticity

import (
	"fmt"
	"math"
	"math/rand/v2"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

// FingerprintObfuscationSystem implements advanced platform signature elimination
type FingerprintObfuscationSystem struct {
	platformManager     *genomics.PlatformErrorSignatureManager
	qualityModeler      *genomics.MultiPlatformQualityModeler
	adversarialDefense  *AdversarialDefenseSystem
	obfuscationStrategies map[string]*ObfuscationStrategy
	fingerprintAnalyzers  map[string]*FingerprintAnalyzer
	crossPlatformModels   map[string]*CrossPlatformModel
	rng                 *rand.Rand
}

// ObfuscationStrategy represents a specific platform signature obfuscation technique
type ObfuscationStrategy struct {
	StrategyName          string
	ObfuscationType       string // "Error_Pattern_Masking", "Quality_Homogenization", "Signature_Morphing", "Platform_Mimicry"
	TargetPlatforms       []genomics.SequencingPlatform
	ObfuscationStrength   float64 // 0.0-1.0 scale
	PreservationScore     float64 // Data quality preservation
	ComputationalCost     float64
	DetectionResistance   float64 // Resistance to fingerprint detection
	ObfuscationParameters map[string]interface{}
}

// FingerprintAnalyzer detects platform-specific signatures in genomic data
type FingerprintAnalyzer struct {
	AnalyzerName       string
	AnalysisType       string // "Error_Pattern", "Quality_Distribution", "GC_Bias", "Read_Length", "Technology_Signature"
	DetectionMethod    string // "Statistical", "ML_Based", "Pattern_Matching", "Ensemble"
	Sensitivity        float64
	Specificity        float64
	AccuracyScore      float64
	FeatureWeights     map[string]float64
	DetectionThreshold float64
}

// CrossPlatformModel represents a unified model that appears platform-agnostic
type CrossPlatformModel struct {
	ModelName            string
	UnifiedErrorProfile  UnifiedErrorProfile
	UnifiedQualityProfile UnifiedQualityProfile
	PopulationAdaptation  map[population.Population]*PopulationAdaptation
	PlatformCompatibility map[genomics.SequencingPlatform]float64
	IndistinguishabilityScore float64
}

// UnifiedErrorProfile creates platform-agnostic error characteristics
type UnifiedErrorProfile struct {
	BaseErrorRate        float64
	ErrorRateVariance    float64
	SubstitutionPattern  map[string]float64
	IndelPattern         map[string]float64
	ContextDependency    map[string]float64
	QualityCorrelation   float64
	HomopolymerHandling  HomopolymerUnification
	GCBiasNeutralization GCBiasNeutralization
}

// UnifiedQualityProfile creates platform-agnostic quality characteristics
type UnifiedQualityProfile struct {
	TargetMeanQuality    float64
	TargetStdQuality     float64
	QualityRange         [2]int
	CycleDecayPattern    []float64
	UniformDistribution  bool
	PopulationNeutrality bool
	PlatformNormalization map[genomics.SequencingPlatform]QualityNormalization
}

// PopulationAdaptation maintains population genetics while obfuscating platform signatures
type PopulationAdaptation struct {
	Population            population.Population
	GeneticSignatureStrength float64 // Preserve population genetics
	PlatformSignatureStrength float64 // Minimize platform signatures
	AdaptationParameters     map[string]float64
	BiologicalConstraints    BiologicalConstraints
}

// ObfuscationResult stores results from platform fingerprint obfuscation
type ObfuscationResult struct {
	OriginalPlatform         genomics.SequencingPlatform
	ObfuscatedSequences      []string
	ObfuscatedQualities      [][]int
	DetectionProbability     float64
	FingerprintStrength      float64
	QualityPreservation      float64
	PopulationPreservation   float64
	ObfuscationEffectiveness float64
	AnalysisResults          map[string]*FingerprintAnalysisResult
}

// FingerprintAnalysisResult stores detailed fingerprint detection results
type FingerprintAnalysisResult struct {
	AnalyzerName         string
	PlatformDetected     genomics.SequencingPlatform
	ConfidenceScore      float64
	SignatureFeatures    map[string]float64
	DetectionCertainty   float64
	ObfuscationSuccess   bool
}

// HomopolymerUnification handles homopolymer run normalization across platforms
type HomopolymerUnification struct {
	TargetErrorRate      float64
	MaxRunLength         int
	UnifiedPenalties     map[int]int
	BaseNeutralization   bool
	AdaptiveCorrection   bool
}

// GCBiasNeutralization eliminates platform-specific GC bias patterns
type GCBiasNeutralization struct {
	TargetGCContent      float64
	BiasStrength         float64
	WindowSize           int
	NeutralizationMethod string // "Linear", "Sigmoid", "Uniform"
	PopulationPreservation bool
}

// QualityNormalization normalizes quality scores across platforms
type QualityNormalization struct {
	SourcePlatform       genomics.SequencingPlatform
	TargetDistribution   string // "Uniform", "Normal", "Beta"
	MeanAdjustment       float64
	VarianceAdjustment   float64
	RangeMapping         [2][2]int // [source_range][target_range]
}

// BiologicalConstraints ensures obfuscation preserves biological realism
type BiologicalConstraints struct {
	MinDiversityPreservation float64
	MaxErrorRateDeviation    float64
	GeneticDriftTolerance    float64
	PopulationStructureIntegrity float64
	BiologicalPlausibility   float64
}

// NewFingerprintObfuscationSystem creates a new platform fingerprint obfuscation system
func NewFingerprintObfuscationSystem() *FingerprintObfuscationSystem {
	return &FingerprintObfuscationSystem{
		platformManager:       genomics.NewPlatformErrorSignatureManager(),
		qualityModeler:        genomics.NewMultiPlatformQualityModeler(),
		adversarialDefense:    NewAdversarialDefenseSystem(),
		obfuscationStrategies: make(map[string]*ObfuscationStrategy),
		fingerprintAnalyzers:  make(map[string]*FingerprintAnalyzer),
		crossPlatformModels:   make(map[string]*CrossPlatformModel),
		rng:                   rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
	}
}

// InitializeObfuscationStrategies sets up comprehensive platform signature obfuscation techniques
func (fos *FingerprintObfuscationSystem) InitializeObfuscationStrategies() {
	
	// Error Pattern Masking - Neutralize platform-specific error signatures
	fos.obfuscationStrategies["Error_Pattern_Masking"] = &ObfuscationStrategy{
		StrategyName:      "Error Pattern Masking",
		ObfuscationType:   "Error_Pattern_Masking",
		TargetPlatforms:   []genomics.SequencingPlatform{genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq, genomics.PacBioRevio, genomics.DNBSeqT20x2, genomics.ElementAVITI, genomics.PromethION24},
		ObfuscationStrength: 0.85,
		PreservationScore:   0.92,
		ComputationalCost:   0.40,
		DetectionResistance: 0.88,
		ObfuscationParameters: map[string]interface{}{
			"error_rate_normalization": 0.001,     // Target unified error rate
			"pattern_randomization":    0.3,       // Randomization strength
			"context_neutralization":   true,      // Remove context-specific patterns
			"signature_morphing":       0.7,       // Transform signatures
			"biological_preservation":  0.95,      // Preserve biological realism
		},
	}

	// Quality Score Homogenization - Eliminate platform-specific quality distributions
	fos.obfuscationStrategies["Quality_Homogenization"] = &ObfuscationStrategy{
		StrategyName:      "Quality Score Homogenization",
		ObfuscationType:   "Quality_Homogenization",
		TargetPlatforms:   []genomics.SequencingPlatform{genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq, genomics.PacBioRevio, genomics.DNBSeqT20x2, genomics.ElementAVITI, genomics.PromethION24},
		ObfuscationStrength: 0.90,
		PreservationScore:   0.88,
		ComputationalCost:   0.30,
		DetectionResistance: 0.92,
		ObfuscationParameters: map[string]interface{}{
			"target_mean_quality":    35.0,        // Unified mean quality
			"target_std_quality":     4.0,         // Unified standard deviation
			"quality_range":          [2]int{2, 41}, // Unified quality range
			"distribution_type":      "normal",     // Target distribution
			"cycle_normalization":    true,        // Normalize cycle effects
			"platform_mapping":       true,        // Map between platforms
		},
	}

	// Advanced Signature Morphing - Transform platform signatures into neutral patterns
	fos.obfuscationStrategies["Advanced_Signature_Morphing"] = &ObfuscationStrategy{
		StrategyName:      "Advanced Signature Morphing",
		ObfuscationType:   "Signature_Morphing",
		TargetPlatforms:   []genomics.SequencingPlatform{genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq, genomics.PacBioRevio, genomics.DNBSeqT20x2, genomics.ElementAVITI, genomics.PromethION24},
		ObfuscationStrength: 0.93,
		PreservationScore:   0.85,
		ComputationalCost:   0.70,
		DetectionResistance: 0.95,
		ObfuscationParameters: map[string]interface{}{
			"morphing_complexity":     5,          // Complexity of transformation
			"signature_scrambling":    0.8,        // Scrambling intensity
			"temporal_variation":      true,       // Time-dependent morphing
			"adaptive_morphing":       true,       // Adaptive to detection attempts
			"cross_platform_mimicry":  0.6,       // Mimic other platforms
			"steganographic_embedding": false,     // Hide signatures in noise
		},
	}

	// Platform Mimicry - Make sequences appear from different platforms
	fos.obfuscationStrategies["Cross_Platform_Mimicry"] = &ObfuscationStrategy{
		StrategyName:      "Cross-Platform Mimicry",
		ObfuscationType:   "Platform_Mimicry",
		TargetPlatforms:   []genomics.SequencingPlatform{genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq, genomics.PacBioRevio, genomics.DNBSeqT20x2, genomics.ElementAVITI, genomics.PromethION24},
		ObfuscationStrength: 0.88,
		PreservationScore:   0.90,
		ComputationalCost:   0.60,
		DetectionResistance: 0.91,
		ObfuscationParameters: map[string]interface{}{
			"target_platform_rotation": []string{"NovaSeq6000", "PacBioRevio", "PromethION24"},
			"mimicry_fidelity":         0.85,      // Accuracy of mimicry
			"signature_transplantation": true,     // Transplant signatures
			"hybrid_signatures":        0.4,       // Create hybrid patterns
			"detection_evasion":        0.9,       // Evasion strength
		},
	}

	// GC Bias Neutralization - Eliminate GC-content-based platform detection
	fos.obfuscationStrategies["GC_Bias_Neutralization"] = &ObfuscationStrategy{
		StrategyName:      "GC Bias Neutralization",
		ObfuscationType:   "GC_Bias_Neutralization",
		TargetPlatforms:   []genomics.SequencingPlatform{genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq, genomics.PacBioRevio, genomics.DNBSeqT20x2, genomics.ElementAVITI, genomics.PromethION24},
		ObfuscationStrength: 0.82,
		PreservationScore:   0.94,
		ComputationalCost:   0.35,
		DetectionResistance: 0.86,
		ObfuscationParameters: map[string]interface{}{
			"target_gc_content":       0.41,       // Neutral GC content
			"bias_strength_reduction": 0.95,       // Reduce bias strength
			"window_size":             100,        // Analysis window size
			"population_preservation": true,       // Preserve population GC patterns
			"adaptive_correction":     true,       // Adaptive bias correction
		},
	}

	// Homopolymer Normalization - Standardize homopolymer error patterns
	fos.obfuscationStrategies["Homopolymer_Normalization"] = &ObfuscationStrategy{
		StrategyName:      "Homopolymer Pattern Normalization",
		ObfuscationType:   "Homopolymer_Normalization",
		TargetPlatforms:   []genomics.SequencingPlatform{genomics.PacBioRevio, genomics.PromethION24, genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq},
		ObfuscationStrength: 0.87,
		PreservationScore:   0.91,
		ComputationalCost:   0.45,
		DetectionResistance: 0.89,
		ObfuscationParameters: map[string]interface{}{
			"max_run_length":          15,         // Maximum homopolymer length
			"error_rate_normalization": 0.02,     // Normalized error rate
			"base_specific_correction": true,      // Base-specific correction
			"quality_penalty_removal":  0.8,      // Remove platform-specific penalties
			"biological_realism":       0.93,     // Maintain biological realism
		},
	}
}

// InitializeFingerprintAnalyzers sets up advanced platform signature detection systems
func (fos *FingerprintObfuscationSystem) InitializeFingerprintAnalyzers() {
	
	// Error Pattern Analyzer - Detects platform-specific error signatures
	fos.fingerprintAnalyzers["Error_Pattern_Analyzer"] = &FingerprintAnalyzer{
		AnalyzerName:       "Error Pattern Analyzer",
		AnalysisType:       "Error_Pattern",
		DetectionMethod:    "ML_Based",
		Sensitivity:        0.92,
		Specificity:        0.88,
		AccuracyScore:      0.90,
		DetectionThreshold: 0.75,
		FeatureWeights: map[string]float64{
			"substitution_rate":     0.30,
			"indel_rate":           0.25,
			"context_dependency":   0.20,
			"quality_correlation":  0.15,
			"homopolymer_patterns": 0.10,
		},
	}

	// Quality Distribution Analyzer - Detects platform-specific quality patterns
	fos.fingerprintAnalyzers["Quality_Distribution_Analyzer"] = &FingerprintAnalyzer{
		AnalyzerName:       "Quality Distribution Analyzer",
		AnalysisType:       "Quality_Distribution",
		DetectionMethod:    "Statistical",
		Sensitivity:        0.89,
		Specificity:        0.91,
		AccuracyScore:      0.90,
		DetectionThreshold: 0.70,
		FeatureWeights: map[string]float64{
			"mean_quality":       0.25,
			"quality_variance":   0.25,
			"quality_range":      0.20,
			"cycle_decay":        0.20,
			"distribution_shape": 0.10,
		},
	}

	// Technology Signature Analyzer - Detects technology-specific patterns
	fos.fingerprintAnalyzers["Technology_Signature_Analyzer"] = &FingerprintAnalyzer{
		AnalyzerName:       "Technology Signature Analyzer",
		AnalysisType:       "Technology_Signature",
		DetectionMethod:    "Ensemble",
		Sensitivity:        0.94,
		Specificity:        0.86,
		AccuracyScore:      0.90,
		DetectionThreshold: 0.80,
		FeatureWeights: map[string]float64{
			"sequencing_chemistry": 0.35,
			"read_length_pattern":  0.25,
			"error_profile":        0.20,
			"quality_signature":    0.15,
			"platform_artifacts":   0.05,
		},
	}

	// GC Bias Analyzer - Detects GC-content-based platform signatures
	fos.fingerprintAnalyzers["GC_Bias_Analyzer"] = &FingerprintAnalyzer{
		AnalyzerName:       "GC Bias Analyzer",
		AnalysisType:       "GC_Bias",
		DetectionMethod:    "Pattern_Matching",
		Sensitivity:        0.85,
		Specificity:        0.93,
		AccuracyScore:      0.89,
		DetectionThreshold: 0.65,
		FeatureWeights: map[string]float64{
			"gc_content_mean":      0.30,
			"gc_content_variance":  0.25,
			"gc_bias_strength":     0.25,
			"window_patterns":      0.15,
			"population_deviation": 0.05,
		},
	}

	// Advanced Ensemble Analyzer - Combines multiple detection methods
	fos.fingerprintAnalyzers["Advanced_Ensemble_Analyzer"] = &FingerprintAnalyzer{
		AnalyzerName:       "Advanced Ensemble Analyzer",
		AnalysisType:       "Multi_Modal",
		DetectionMethod:    "Ensemble",
		Sensitivity:        0.96,
		Specificity:        0.84,
		AccuracyScore:      0.92,
		DetectionThreshold: 0.85,
		FeatureWeights: map[string]float64{
			"error_patterns":       0.25,
			"quality_distributions": 0.25,
			"technology_signatures": 0.20,
			"gc_bias_patterns":     0.15,
			"read_length_analysis": 0.10,
			"ensemble_consensus":   0.05,
		},
	}
}

// InitializeCrossPlatformModels creates unified platform-agnostic models
func (fos *FingerprintObfuscationSystem) InitializeCrossPlatformModels() {
	
	// Universal Neutral Model - Maximally platform-agnostic
	fos.crossPlatformModels["Universal_Neutral_Model"] = &CrossPlatformModel{
		ModelName: "Universal Neutral Model",
		UnifiedErrorProfile: UnifiedErrorProfile{
			BaseErrorRate:     0.001,  // Conservative neutral error rate
			ErrorRateVariance: 0.0002, // Low variance for consistency
			SubstitutionPattern: map[string]float64{
				"A->T": 0.25, "A->G": 0.25, "A->C": 0.25,
				"T->A": 0.25, "T->G": 0.25, "T->C": 0.25,
				"G->A": 0.25, "G->T": 0.25, "G->C": 0.25,
				"C->A": 0.25, "C->T": 0.25, "C->G": 0.25,
			},
			IndelPattern: map[string]float64{
				"insertion": 0.0001,
				"deletion":  0.0001,
			},
			ContextDependency: map[string]float64{
				"GC": 1.0, "AT": 1.0, "CpG": 1.0, // Neutral context
			},
			QualityCorrelation: -0.75,
			HomopolymerHandling: HomopolymerUnification{
				TargetErrorRate:    0.01,
				MaxRunLength:       10,
				UnifiedPenalties:   map[int]int{3: 2, 5: 4, 10: 8},
				BaseNeutralization: true,
				AdaptiveCorrection: true,
			},
			GCBiasNeutralization: GCBiasNeutralization{
				TargetGCContent:        0.41,
				BiasStrength:           0.05,  // Minimal bias
				WindowSize:             100,
				NeutralizationMethod:   "Linear",
				PopulationPreservation: true,
			},
		},
		UnifiedQualityProfile: UnifiedQualityProfile{
			TargetMeanQuality:    35.0,
			TargetStdQuality:     4.0,
			QualityRange:         [2]int{2, 41},
			CycleDecayPattern:    fos.generateNeutralCycleDecay(150),
			UniformDistribution:  false,
			PopulationNeutrality: true,
			PlatformNormalization: map[genomics.SequencingPlatform]QualityNormalization{
				genomics.NovaSeq6000: {
					SourcePlatform:     genomics.NovaSeq6000,
					TargetDistribution: "Normal",
					MeanAdjustment:     0.0,
					VarianceAdjustment: 1.0,
					RangeMapping:       [2][2]int{{2, 41}, {2, 41}},
				},
				// Additional platform mappings...
			},
		},
		IndistinguishabilityScore: 0.95,
	}

	// Population-Aware Neutral Model - Preserves population genetics while eliminating platform signatures
	fos.crossPlatformModels["Population_Aware_Neutral_Model"] = &CrossPlatformModel{
		ModelName: "Population-Aware Neutral Model",
		UnifiedErrorProfile: UnifiedErrorProfile{
			BaseErrorRate:     0.0012, // Slightly higher for population adaptation
			ErrorRateVariance: 0.0003,
			SubstitutionPattern: map[string]float64{
				"A->T": 0.26, "A->G": 0.24, "A->C": 0.25,
				"T->A": 0.26, "T->G": 0.24, "T->C": 0.25,
				"G->A": 0.24, "G->T": 0.26, "G->C": 0.25,
				"C->A": 0.25, "C->T": 0.26, "C->G": 0.24,
			},
			QualityCorrelation: -0.78,
		},
		PopulationAdaptation: map[population.Population]*PopulationAdaptation{
			population.YRI: {
				Population:                population.YRI,
				GeneticSignatureStrength:  0.95, // Preserve African diversity
				PlatformSignatureStrength: 0.05, // Minimize platform signatures
				AdaptationParameters: map[string]float64{
					"diversity_preservation": 0.98,
					"gc_bias_adjustment":     0.95,
					"error_rate_modifier":    0.98,
				},
				BiologicalConstraints: BiologicalConstraints{
					MinDiversityPreservation:     0.95,
					MaxErrorRateDeviation:        0.05,
					GeneticDriftTolerance:        0.02,
					PopulationStructureIntegrity: 0.98,
					BiologicalPlausibility:       0.96,
				},
			},
			population.CEU: {
				Population:                population.CEU,
				GeneticSignatureStrength:  0.92, // Preserve European patterns
				PlatformSignatureStrength: 0.05,
				AdaptationParameters: map[string]float64{
					"diversity_preservation": 0.92,
					"gc_bias_adjustment":     1.00,
					"error_rate_modifier":    1.00,
				},
				BiologicalConstraints: BiologicalConstraints{
					MinDiversityPreservation:     0.90,
					MaxErrorRateDeviation:        0.05,
					GeneticDriftTolerance:        0.03,
					PopulationStructureIntegrity: 0.95,
					BiologicalPlausibility:       0.94,
				},
			},
			// Additional population adaptations...
		},
		IndistinguishabilityScore: 0.93,
	}

	// Calculate platform compatibility for all models
	for _, model := range fos.crossPlatformModels {
		model.PlatformCompatibility = fos.calculatePlatformCompatibility(model)
	}
}

// ObfuscatePlatformFingerprints applies comprehensive platform signature obfuscation
func (fos *FingerprintObfuscationSystem) ObfuscatePlatformFingerprints(
	sequences map[genomics.SequencingPlatform]map[population.Population][]string,
	targetObfuscationLevel float64,
) (*map[genomics.SequencingPlatform]*ObfuscationResult, error) {

	if len(sequences) == 0 {
		return nil, fmt.Errorf("no sequences provided for obfuscation")
	}

	// Initialize strategies and analyzers if not already done
	if len(fos.obfuscationStrategies) == 0 {
		fos.InitializeObfuscationStrategies()
	}
	if len(fos.fingerprintAnalyzers) == 0 {
		fos.InitializeFingerprintAnalyzers()
	}
	if len(fos.crossPlatformModels) == 0 {
		fos.InitializeCrossPlatformModels()
	}

	results := make(map[genomics.SequencingPlatform]*ObfuscationResult)

	// Process each platform
	for platform, populationSeqs := range sequences {
		fmt.Printf("Obfuscating platform fingerprints for %s...\n", platform)

		// Analyze original fingerprint strength
		originalFingerprint := fos.analyzePlatformFingerprint(platform, populationSeqs)

		// Select optimal obfuscation strategy
		strategy := fos.selectOptimalObfuscationStrategy(platform, originalFingerprint, targetObfuscationLevel)

		// Apply obfuscation across all populations for this platform
		obfuscatedSeqs := []string{}
		obfuscatedQuals := [][]int{}
		totalPreservation := 0.0
		preservationCount := 0

		for pop, seqs := range populationSeqs {
			// Apply obfuscation to population sequences
			popObfuscated, popQualities, preservation, err := fos.applyObfuscationStrategy(
				strategy, platform, pop, seqs)
			if err != nil {
				return nil, fmt.Errorf("obfuscation strategy %q failed for platform %s: %w",
					strategy.ObfuscationType, platform, err)
			}

			obfuscatedSeqs = append(obfuscatedSeqs, popObfuscated...)
			obfuscatedQuals = append(obfuscatedQuals, popQualities...)
			totalPreservation += preservation
			preservationCount++
		}

		// Analyze obfuscated fingerprint strength
		obfuscatedSequences := map[population.Population][]string{}
		seqIndex := 0
		for pop, seqs := range populationSeqs {
			popSeqs := obfuscatedSeqs[seqIndex : seqIndex+len(seqs)]
			obfuscatedSequences[pop] = popSeqs
			seqIndex += len(seqs)
		}
		
		obfuscatedFingerprint := fos.analyzePlatformFingerprint(platform, obfuscatedSequences)

		// Perform comprehensive fingerprint analysis
		analysisResults := fos.performComprehensiveFingerprintAnalysis(
			platform, obfuscatedSequences, obfuscatedQuals)

		// Calculate metrics
		avgPreservation := totalPreservation / float64(preservationCount)
		detectionProb := fos.calculateDetectionProbability(analysisResults)
		obfuscationEffectiveness := (originalFingerprint - obfuscatedFingerprint) / originalFingerprint

		results[platform] = &ObfuscationResult{
			OriginalPlatform:         platform,
			ObfuscatedSequences:      obfuscatedSeqs,
			ObfuscatedQualities:      obfuscatedQuals,
			DetectionProbability:     detectionProb,
			FingerprintStrength:      obfuscatedFingerprint,
			QualityPreservation:      avgPreservation,
			PopulationPreservation:   fos.calculatePopulationPreservation(populationSeqs, obfuscatedSequences),
			ObfuscationEffectiveness: obfuscationEffectiveness,
			AnalysisResults:          analysisResults,
		}
	}

	return &results, nil
}

// analyzePlatformFingerprint analyzes the strength of platform-specific signatures
func (fos *FingerprintObfuscationSystem) analyzePlatformFingerprint(
	platform genomics.SequencingPlatform,
	populationSeqs map[population.Population][]string,
) float64 {

	// Get platform signature
	signature := fos.platformManager.GetErrorSignature(platform)
	if signature == nil {
		return 0.5 // Unknown platform
	}

	// Calculate fingerprint strength based on multiple factors
	errorPatternStrength := fos.calculateErrorPatternStrength(signature)
	qualityPatternStrength := fos.calculateQualityPatternStrength(platform)
	technologySignatureStrength := fos.calculateTechnologySignatureStrength(platform)

	// Weight the components
	overallStrength := (errorPatternStrength*0.4 + 
					   qualityPatternStrength*0.35 + 
					   technologySignatureStrength*0.25)

	return math.Min(1.0, overallStrength)
}

// selectOptimalObfuscationStrategy selects the best obfuscation strategy for a platform
func (fos *FingerprintObfuscationSystem) selectOptimalObfuscationStrategy(
	platform genomics.SequencingPlatform,
	fingerprintStrength float64,
	targetLevel float64,
) *ObfuscationStrategy {

	bestStrategy := (*ObfuscationStrategy)(nil)
	bestScore := 0.0

	for _, strategy := range fos.obfuscationStrategies {
		// Check if strategy is applicable to this platform
		applicable := false
		for _, targetPlatform := range strategy.TargetPlatforms {
			if targetPlatform == platform {
				applicable = true
				break
			}
		}
		if !applicable {
			continue
		}

		// Calculate strategy suitability score
		effectivenessGap := math.Abs(strategy.ObfuscationStrength - targetLevel)
		resistanceScore := strategy.DetectionResistance
		preservationScore := strategy.PreservationScore
		costPenalty := 1.0 - (strategy.ComputationalCost * 0.3)

		suitabilityScore := resistanceScore * preservationScore * costPenalty * (1.0 - effectivenessGap)

		if suitabilityScore > bestScore {
			bestScore = suitabilityScore
			bestStrategy = strategy
		}
	}

	if bestStrategy == nil {
		// Fallback to a default strategy
		return fos.obfuscationStrategies["Error_Pattern_Masking"]
	}

	return bestStrategy
}

// applyObfuscationStrategy applies a specific obfuscation strategy to sequences
func (fos *FingerprintObfuscationSystem) applyObfuscationStrategy(
	strategy *ObfuscationStrategy,
	platform genomics.SequencingPlatform,
	pop population.Population,
	sequences []string,
) ([]string, [][]int, float64, error) {

	switch strategy.ObfuscationType {
	case "Error_Pattern_Masking":
		return fos.applyErrorPatternMasking(strategy, platform, pop, sequences)
	case "Quality_Homogenization":
		return fos.applyQualityHomogenization(strategy, platform, pop, sequences)
	case "Signature_Morphing":
		return fos.applySignatureMorphing(strategy, platform, pop, sequences)
	case "Platform_Mimicry":
		return fos.applyPlatformMimicry(strategy, platform, pop, sequences)
	case "GC_Bias_Neutralization":
		return fos.applyGCBiasNeutralization(strategy, platform, pop, sequences)
	case "Homopolymer_Normalization":
		return fos.applyHomopolymerNormalization(strategy, platform, pop, sequences)
	default:
		return sequences, fos.generateDefaultQualities(sequences), 1.0, nil
	}
}

// Obfuscation strategy implementations

func (fos *FingerprintObfuscationSystem) applyErrorPatternMasking(
	strategy *ObfuscationStrategy,
	platform genomics.SequencingPlatform,
	pop population.Population,
	sequences []string,
) ([]string, [][]int, float64, error) {

	targetErrorRate, ok := strategy.ObfuscationParameters["error_rate_normalization"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter error_rate_normalization")
	}
	randomizationStrength, ok := strategy.ObfuscationParameters["pattern_randomization"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter pattern_randomization")
	}
	morphingStrength, ok := strategy.ObfuscationParameters["signature_morphing"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter signature_morphing")
	}

	obfuscatedSeqs := make([]string, len(sequences))
	obfuscatedQuals := make([][]int, len(sequences))
	preservationSum := 0.0

	for i, seq := range sequences {
		// Generate original qualities
		originalQualities := fos.qualityModeler.GenerateQualityScores(platform, pop, seq, i+1)

		// Apply unified error pattern
		maskedSeq, maskedQual := fos.applyUnifiedErrorPattern(seq, originalQualities, targetErrorRate, randomizationStrength)

		// Apply signature morphing
		if morphingStrength > 0 {
			maskedSeq = fos.morphSignature(maskedSeq, morphingStrength)
		}

		obfuscatedSeqs[i] = maskedSeq
		obfuscatedQuals[i] = maskedQual

		// Calculate preservation score
		preservation := fos.calculateSequencePreservation(seq, maskedSeq)
		preservationSum += preservation
	}

	avgPreservation := preservationSum / float64(len(sequences))
	return obfuscatedSeqs, obfuscatedQuals, avgPreservation, nil
}

func (fos *FingerprintObfuscationSystem) applyQualityHomogenization(
	strategy *ObfuscationStrategy,
	platform genomics.SequencingPlatform,
	pop population.Population,
	sequences []string,
) ([]string, [][]int, float64, error) {

	targetMean, ok := strategy.ObfuscationParameters["target_mean_quality"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter target_mean_quality")
	}
	targetStd, ok := strategy.ObfuscationParameters["target_std_quality"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter target_std_quality")
	}
	qualityRange, ok := strategy.ObfuscationParameters["quality_range"].([2]int)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing [2]int parameter quality_range")
	}

	obfuscatedSeqs := make([]string, len(sequences))
	obfuscatedQuals := make([][]int, len(sequences))

	for i, seq := range sequences {
		// Keep sequence unchanged for quality homogenization
		obfuscatedSeqs[i] = seq

		// Generate homogenized quality scores
		homogenizedQual := fos.generateHomogenizedQualities(len(seq), targetMean, targetStd, qualityRange)
		obfuscatedQuals[i] = homogenizedQual
	}

	// Quality homogenization preserves sequences perfectly
	return obfuscatedSeqs, obfuscatedQuals, 1.0, nil
}

func (fos *FingerprintObfuscationSystem) applySignatureMorphing(
	strategy *ObfuscationStrategy,
	platform genomics.SequencingPlatform,
	pop population.Population,
	sequences []string,
) ([]string, [][]int, float64, error) {

	morphingComplexity, ok := strategy.ObfuscationParameters["morphing_complexity"].(int)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing int parameter morphing_complexity")
	}
	scramblingIntensity, ok := strategy.ObfuscationParameters["signature_scrambling"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter signature_scrambling")
	}
	mimicryStrength, ok := strategy.ObfuscationParameters["cross_platform_mimicry"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter cross_platform_mimicry")
	}

	obfuscatedSeqs := make([]string, len(sequences))
	obfuscatedQuals := make([][]int, len(sequences))
	preservationSum := 0.0

	for i, seq := range sequences {
		originalQualities := fos.qualityModeler.GenerateQualityScores(platform, pop, seq, i+1)

		// Apply complex signature morphing
		morphedSeq := fos.applyComplexMorphing(seq, morphingComplexity, scramblingIntensity)

		// Apply cross-platform mimicry
		if mimicryStrength > 0 {
			morphedSeq = fos.applyCrossPlatformMimicry(morphedSeq, platform, mimicryStrength)
		}

		// Generate corresponding quality scores
		morphedQual := fos.morphQualityScores(originalQualities, morphingComplexity)

		obfuscatedSeqs[i] = morphedSeq
		obfuscatedQuals[i] = morphedQual

		preservation := fos.calculateSequencePreservation(seq, morphedSeq)
		preservationSum += preservation
	}

	avgPreservation := preservationSum / float64(len(sequences))
	return obfuscatedSeqs, obfuscatedQuals, avgPreservation, nil
}

func (fos *FingerprintObfuscationSystem) applyPlatformMimicry(
	strategy *ObfuscationStrategy,
	platform genomics.SequencingPlatform,
	pop population.Population,
	sequences []string,
) ([]string, [][]int, float64, error) {

	targetPlatforms, ok := strategy.ObfuscationParameters["target_platform_rotation"].([]string)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing []string parameter target_platform_rotation")
	}
	mimicryFidelity, ok := strategy.ObfuscationParameters["mimicry_fidelity"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter mimicry_fidelity")
	}

	obfuscatedSeqs := make([]string, len(sequences))
	obfuscatedQuals := make([][]int, len(sequences))
	preservationSum := 0.0

	for i, seq := range sequences {
		// Select target platform for mimicry
		targetPlatformName := targetPlatforms[i%len(targetPlatforms)]
		targetPlatform := fos.stringToPlatform(targetPlatformName)

		// Apply platform mimicry
		mimicSeq := fos.mimicPlatformSignature(seq, platform, targetPlatform, mimicryFidelity)

		// Generate target platform quality scores
		mimicQual := fos.qualityModeler.GenerateQualityScores(targetPlatform, pop, mimicSeq, i+1)

		obfuscatedSeqs[i] = mimicSeq
		obfuscatedQuals[i] = mimicQual

		preservation := fos.calculateSequencePreservation(seq, mimicSeq)
		preservationSum += preservation
	}

	avgPreservation := preservationSum / float64(len(sequences))
	return obfuscatedSeqs, obfuscatedQuals, avgPreservation, nil
}

func (fos *FingerprintObfuscationSystem) applyGCBiasNeutralization(
	strategy *ObfuscationStrategy,
	platform genomics.SequencingPlatform,
	pop population.Population,
	sequences []string,
) ([]string, [][]int, float64, error) {

	targetGC, ok := strategy.ObfuscationParameters["target_gc_content"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter target_gc_content")
	}
	biasReduction, ok := strategy.ObfuscationParameters["bias_strength_reduction"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter bias_strength_reduction")
	}

	obfuscatedSeqs := make([]string, len(sequences))
	obfuscatedQuals := make([][]int, len(sequences))
	preservationSum := 0.0

	for i, seq := range sequences {
		// Apply GC bias neutralization
		neutralizedSeq := fos.neutralizeGCBias(seq, targetGC, biasReduction)

		// Generate quality scores for neutralized sequence
		neutralizedQual := fos.qualityModeler.GenerateQualityScores(platform, pop, neutralizedSeq, i+1)

		obfuscatedSeqs[i] = neutralizedSeq
		obfuscatedQuals[i] = neutralizedQual

		preservation := fos.calculateSequencePreservation(seq, neutralizedSeq)
		preservationSum += preservation
	}

	avgPreservation := preservationSum / float64(len(sequences))
	return obfuscatedSeqs, obfuscatedQuals, avgPreservation, nil
}

func (fos *FingerprintObfuscationSystem) applyHomopolymerNormalization(
	strategy *ObfuscationStrategy,
	platform genomics.SequencingPlatform,
	pop population.Population,
	sequences []string,
) ([]string, [][]int, float64, error) {

	maxRunLength, ok := strategy.ObfuscationParameters["max_run_length"].(int)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing int parameter max_run_length")
	}
	errorRateNormalization, ok := strategy.ObfuscationParameters["error_rate_normalization"].(float64)
	if !ok {
		return nil, nil, 0, fmt.Errorf("obfuscation strategy missing float64 parameter error_rate_normalization")
	}

	obfuscatedSeqs := make([]string, len(sequences))
	obfuscatedQuals := make([][]int, len(sequences))
	preservationSum := 0.0

	for i, seq := range sequences {
		// Apply homopolymer normalization
		normalizedSeq := fos.normalizeHomopolymers(seq, maxRunLength, errorRateNormalization)

		// Generate quality scores for normalized sequence
		normalizedQual := fos.qualityModeler.GenerateQualityScores(platform, pop, normalizedSeq, i+1)

		obfuscatedSeqs[i] = normalizedSeq
		obfuscatedQuals[i] = normalizedQual

		preservation := fos.calculateSequencePreservation(seq, normalizedSeq)
		preservationSum += preservation
	}

	avgPreservation := preservationSum / float64(len(sequences))
	return obfuscatedSeqs, obfuscatedQuals, avgPreservation, nil
}

// Helper methods for obfuscation implementation

func (fos *FingerprintObfuscationSystem) applyUnifiedErrorPattern(
	sequence string, qualities []int, targetErrorRate, randomizationStrength float64,
) (string, []int) {

	modifiedSeq := []byte(sequence)
	modifiedQual := make([]int, len(qualities))
	copy(modifiedQual, qualities)

	for i, base := range modifiedSeq {
		// Apply unified error rate with randomization
		errorProb := targetErrorRate * (1.0 + (fos.rng.Float64()-0.5)*randomizationStrength)
		
		if fos.rng.Float64() < errorProb {
			// Apply error
			bases := []byte{'A', 'T', 'G', 'C'}
			for {
				newBase := bases[fos.rng.IntN(4)]
				if newBase != base {
					modifiedSeq[i] = newBase
					// Adjust quality score
					if i < len(modifiedQual) {
						qualityReduction := 3 + fos.rng.IntN(5)
						modifiedQual[i] = int(math.Max(2, float64(modifiedQual[i])-float64(qualityReduction)))
					}
					break
				}
			}
		}
	}

	return string(modifiedSeq), modifiedQual
}

func (fos *FingerprintObfuscationSystem) morphSignature(sequence string, morphingStrength float64) string {
	if morphingStrength <= 0 {
		return sequence
	}

	modifiedSeq := []byte(sequence)
	morphingRate := morphingStrength * 0.1 // Scale morphing strength

	for i := range modifiedSeq {
		if fos.rng.Float64() < morphingRate {
			// Apply subtle morphing
			bases := []byte{'A', 'T', 'G', 'C'}
			newBase := bases[fos.rng.IntN(4)]
			modifiedSeq[i] = newBase
		}
	}

	return string(modifiedSeq)
}

func (fos *FingerprintObfuscationSystem) generateHomogenizedQualities(
	length int, targetMean, targetStd float64, qualityRange [2]int,
) []int {

	qualities := make([]int, length)
	
	for i := 0; i < length; i++ {
		// Generate quality from normal distribution
		quality := fos.rng.NormFloat64()*targetStd + targetMean
		
		// Clamp to range
		intQuality := int(math.Max(float64(qualityRange[0]), 
			math.Min(float64(qualityRange[1]), quality)))
		
		qualities[i] = intQuality
	}

	return qualities
}

func (fos *FingerprintObfuscationSystem) applyComplexMorphing(
	sequence string, complexity int, scramblingIntensity float64,
) string {

	modifiedSeq := []byte(sequence)
	
	// Apply multiple rounds of morphing based on complexity
	for round := 0; round < complexity; round++ {
		for i := range modifiedSeq {
			if fos.rng.Float64() < scramblingIntensity/float64(complexity) {
				// Complex morphing pattern
				if i > 0 && i < len(modifiedSeq)-1 {
					// Context-aware morphing
					prevBase := modifiedSeq[i-1]
					nextBase := modifiedSeq[i+1]
					
					// Choose morphing based on context
					if prevBase == nextBase {
						// Break homopolymer patterns
						bases := []byte{'A', 'T', 'G', 'C'}
						for {
							newBase := bases[fos.rng.IntN(4)]
							if newBase != prevBase {
								modifiedSeq[i] = newBase
								break
							}
						}
					}
				}
			}
		}
	}

	return string(modifiedSeq)
}

func (fos *FingerprintObfuscationSystem) applyCrossPlatformMimicry(
	sequence string, sourcePlatform genomics.SequencingPlatform, mimicryStrength float64,
) string {

	// Select random target platform for mimicry
	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq,
		genomics.PacBioRevio, genomics.DNBSeqT20x2, genomics.ElementAVITI, genomics.PromethION24,
	}
	
	var targetPlatform genomics.SequencingPlatform
	for {
		targetPlatform = platforms[fos.rng.IntN(len(platforms))]
		if targetPlatform != sourcePlatform {
			break
		}
	}

	return fos.mimicPlatformSignature(sequence, sourcePlatform, targetPlatform, mimicryStrength)
}

func (fos *FingerprintObfuscationSystem) morphQualityScores(qualities []int, complexity int) []int {
	morphedQual := make([]int, len(qualities))
	copy(morphedQual, qualities)

	morphingRate := 0.1 / float64(complexity)

	for i := range morphedQual {
		if fos.rng.Float64() < morphingRate {
			// Apply quality morphing
			adjustment := fos.rng.IntN(6) - 3 // -3 to +3
			morphedQual[i] = int(math.Max(2, 
				math.Min(41, float64(morphedQual[i])+float64(adjustment))))
		}
	}

	return morphedQual
}

func (fos *FingerprintObfuscationSystem) mimicPlatformSignature(
	sequence string, sourcePlatform, targetPlatform genomics.SequencingPlatform, fidelity float64,
) string {

	// Get target platform signature characteristics
	targetSignature := fos.platformManager.GetErrorSignature(targetPlatform)
	if targetSignature == nil {
		return sequence // No signature to mimic
	}

	modifiedSeq := []byte(sequence)

	// Apply target platform characteristics with specified fidelity
	for i := range modifiedSeq {
		if fos.rng.Float64() < fidelity*0.01 { // Scale fidelity
			// Apply target platform error patterns
			if subPattern, exists := targetSignature.ErrorPatterns["substitution"]; exists {
				if i < len(subPattern.FrequencyPattern) {
					errorRate := subPattern.FrequencyPattern[i]
					if fos.rng.Float64() < errorRate {
						// Apply target platform substitution pattern
						bases := []byte{'A', 'T', 'G', 'C'}
						for {
							newBase := bases[fos.rng.IntN(4)]
							if newBase != modifiedSeq[i] {
								modifiedSeq[i] = newBase
								break
							}
						}
					}
				}
			}
		}
	}

	return string(modifiedSeq)
}

func (fos *FingerprintObfuscationSystem) neutralizeGCBias(
	sequence string, targetGC, biasReduction float64,
) string {

	currentGC := fos.calculateGCContent(sequence)
	if math.Abs(currentGC-targetGC) < 0.01 {
		return sequence // Already neutral
	}

	modifiedSeq := []byte(sequence)
	gcDifference := targetGC - currentGC
	adjustmentRate := math.Abs(gcDifference) * biasReduction

	for i := range modifiedSeq {
		if fos.rng.Float64() < adjustmentRate {
			currentBase := modifiedSeq[i]
			
			if gcDifference > 0 { // Need more GC
				if currentBase == 'A' || currentBase == 'T' {
					// Replace with G or C
					modifiedSeq[i] = []byte{'G', 'C'}[fos.rng.IntN(2)]
				}
			} else { // Need less GC
				if currentBase == 'G' || currentBase == 'C' {
					// Replace with A or T
					modifiedSeq[i] = []byte{'A', 'T'}[fos.rng.IntN(2)]
				}
			}
		}
	}

	return string(modifiedSeq)
}

func (fos *FingerprintObfuscationSystem) normalizeHomopolymers(
	sequence string, maxRunLength int, errorRate float64,
) string {

	modifiedSeq := []byte(sequence)
	
	// Find and break long homopolymer runs
	i := 0
	for i < len(modifiedSeq) {
		runLength := 1
		currentBase := modifiedSeq[i]
		
		// Count run length
		for j := i + 1; j < len(modifiedSeq) && modifiedSeq[j] == currentBase; j++ {
			runLength++
		}
		
		// Break run if too long
		if runLength > maxRunLength {
			// Insert errors to break the run
			for j := i + maxRunLength - 1; j < i + runLength; j += maxRunLength {
				if j < len(modifiedSeq) && fos.rng.Float64() < errorRate {
					// Replace with different base
					bases := []byte{'A', 'T', 'G', 'C'}
					for {
						newBase := bases[fos.rng.IntN(4)]
						if newBase != currentBase {
							modifiedSeq[j] = newBase
							break
						}
					}
				}
			}
		}
		
		i += runLength
	}

	return string(modifiedSeq)
}

// Analysis and assessment methods

func (fos *FingerprintObfuscationSystem) performComprehensiveFingerprintAnalysis(
	platform genomics.SequencingPlatform,
	sequences map[population.Population][]string,
	qualities [][]int,
) map[string]*FingerprintAnalysisResult {

	results := make(map[string]*FingerprintAnalysisResult)

	for analyzerName, analyzer := range fos.fingerprintAnalyzers {
		result := fos.runFingerprintAnalyzer(analyzer, platform, sequences, qualities)
		results[analyzerName] = result
	}

	return results
}

func (fos *FingerprintObfuscationSystem) runFingerprintAnalyzer(
	analyzer *FingerprintAnalyzer,
	truePlatform genomics.SequencingPlatform,
	sequences map[population.Population][]string,
	qualities [][]int,
) *FingerprintAnalysisResult {

	// Extract features for analysis
	features := fos.extractAnalysisFeatures(analyzer, sequences, qualities)
	
	// Simulate detection analysis
	detectionScore := fos.simulateDetectionAnalysis(analyzer, features)
	
	// Determine detected platform
	detectedPlatform := fos.simulatePlatformDetection(analyzer, features, detectionScore)
	
	// Calculate confidence
	confidence := detectionScore * analyzer.AccuracyScore

	return &FingerprintAnalysisResult{
		AnalyzerName:       analyzer.AnalyzerName,
		PlatformDetected:   detectedPlatform,
		ConfidenceScore:    confidence,
		SignatureFeatures:  features,
		DetectionCertainty: detectionScore,
		ObfuscationSuccess: detectedPlatform != truePlatform || confidence < analyzer.DetectionThreshold,
	}
}

// Helper calculation methods

func (fos *FingerprintObfuscationSystem) calculateErrorPatternStrength(signature *genomics.ErrorSignature) float64 {
	if signature == nil {
		return 0.0
	}

	strength := 0.0
	patternCount := 0

	for _, pattern := range signature.ErrorPatterns {
		avgRate := 0.0
		for _, rate := range pattern.FrequencyPattern {
			avgRate += rate
		}
		if len(pattern.FrequencyPattern) > 0 {
			avgRate /= float64(len(pattern.FrequencyPattern))
		}
		strength += avgRate
		patternCount++
	}

	if patternCount > 0 {
		strength /= float64(patternCount)
	}

	return math.Min(1.0, strength*1000) // Scale to 0-1 range
}

func (fos *FingerprintObfuscationSystem) calculateQualityPatternStrength(platform genomics.SequencingPlatform) float64 {
	profile := fos.qualityModeler.GetPlatformProfile(platform)
	if profile == nil {
		return 0.0
	}

	// Calculate distinctiveness based on quality characteristics
	meanQuality := profile.QualityScoreDistribution.MeanQuality
	stdQuality := profile.QualityScoreDistribution.StdQuality
	errorRate := profile.ErrorProfile.BaseErrorRate

	// Higher deviation from neutral values = higher fingerprint strength
	meanDeviation := math.Abs(meanQuality-30.0) / 30.0  // Neutral mean ~30
	stdDeviation := math.Abs(stdQuality-5.0) / 5.0      // Neutral std ~5
	errorDeviation := math.Abs(errorRate-0.001) / 0.001 // Neutral error ~0.001

	strength := (meanDeviation + stdDeviation + errorDeviation) / 3.0
	return math.Min(1.0, strength)
}

func (fos *FingerprintObfuscationSystem) calculateTechnologySignatureStrength(platform genomics.SequencingPlatform) float64 {
	// Technology-specific signature strength based on platform type
	technologyStrengths := map[genomics.SequencingPlatform]float64{
		genomics.NovaSeq6000:   0.75, // Strong Illumina signature
		genomics.NextSeq2000:   0.72, // Strong Illumina signature
		genomics.MiSeq:         0.70, // Moderate Illumina signature
		genomics.PacBioRevio:   0.85, // Very strong PacBio signature
		genomics.DNBSeqT20x2:   0.78, // Strong MGI signature
		genomics.ElementAVITI:  0.80, // Strong Element signature
		genomics.PromethION24:  0.88, // Very strong ONT signature
	}

	if strength, exists := technologyStrengths[platform]; exists {
		return strength
	}
	return 0.5 // Unknown platform
}

func (fos *FingerprintObfuscationSystem) generateNeutralCycleDecay(readLength int) []float64 {
	decay := make([]float64, readLength)
	for i := 0; i < readLength; i++ {
		// Minimal linear decay
		decay[i] = 35.0 - (float64(i)/float64(readLength))*2.0
	}
	return decay
}

func (fos *FingerprintObfuscationSystem) calculatePlatformCompatibility(model *CrossPlatformModel) map[genomics.SequencingPlatform]float64 {
	compatibility := make(map[genomics.SequencingPlatform]float64)
	
	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq, genomics.PacBioRevio,
		genomics.DNBSeqT20x2, genomics.ElementAVITI, genomics.PromethION24,
	}

	for _, platform := range platforms {
		// Calculate compatibility based on unified profiles
		errorCompatibility := 1.0 - math.Abs(model.UnifiedErrorProfile.BaseErrorRate-0.001)
		qualityCompatibility := 1.0 - math.Abs(model.UnifiedQualityProfile.TargetMeanQuality-35.0)/35.0
		
		compatibility[platform] = math.Min(1.0, (errorCompatibility+qualityCompatibility)/2.0)
	}

	return compatibility
}

func (fos *FingerprintObfuscationSystem) calculateSequencePreservation(original, modified string) float64 {
	if len(original) != len(modified) {
		return 0.0
	}

	matches := 0
	for i := 0; i < len(original); i++ {
		if original[i] == modified[i] {
			matches++
		}
	}

	return float64(matches) / float64(len(original))
}

func (fos *FingerprintObfuscationSystem) calculatePopulationPreservation(
	original map[population.Population][]string,
	obfuscated map[population.Population][]string,
) float64 {

	totalPreservation := 0.0
	popCount := 0

	for pop, originalSeqs := range original {
		if obfuscatedSeqs, exists := obfuscated[pop]; exists {
			minLen := int(math.Min(float64(len(originalSeqs)), float64(len(obfuscatedSeqs))))
			
			popPreservation := 0.0
			for i := 0; i < minLen; i++ {
				preservation := fos.calculateSequencePreservation(originalSeqs[i], obfuscatedSeqs[i])
				popPreservation += preservation
			}
			
			if minLen > 0 {
				popPreservation /= float64(minLen)
			}
			
			totalPreservation += popPreservation
			popCount++
		}
	}

	if popCount == 0 {
		return 0.0
	}

	return totalPreservation / float64(popCount)
}

func (fos *FingerprintObfuscationSystem) calculateDetectionProbability(
	analysisResults map[string]*FingerprintAnalysisResult,
) float64 {

	successfulObfuscations := 0
	totalAnalyses := len(analysisResults)

	for _, result := range analysisResults {
		if result.ObfuscationSuccess {
			successfulObfuscations++
		}
	}

	if totalAnalyses == 0 {
		return 1.0 // Maximum detection probability if no analysis
	}

	// Detection probability is inverse of obfuscation success rate
	obfuscationSuccessRate := float64(successfulObfuscations) / float64(totalAnalyses)
	return 1.0 - obfuscationSuccessRate
}

func (fos *FingerprintObfuscationSystem) calculateGCContent(sequence string) float64 {
	if len(sequence) == 0 {
		return 0.5
	}

	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}

	return float64(gcCount) / float64(len(sequence))
}

func (fos *FingerprintObfuscationSystem) generateDefaultQualities(sequences []string) [][]int {
	qualities := make([][]int, len(sequences))
	
	for i, seq := range sequences {
		qualities[i] = make([]int, len(seq))
		for j := range qualities[i] {
			qualities[i][j] = 30 + fos.rng.IntN(10) // Q30-Q40
		}
	}

	return qualities
}

func (fos *FingerprintObfuscationSystem) stringToPlatform(platformName string) genomics.SequencingPlatform {
	switch platformName {
	case "NovaSeq6000":
		return genomics.NovaSeq6000
	case "NextSeq2000":
		return genomics.NextSeq2000
	case "MiSeq":
		return genomics.MiSeq
	case "PacBioRevio":
		return genomics.PacBioRevio
	case "DNBSeqT20x2":
		return genomics.DNBSeqT20x2
	case "ElementAVITI":
		return genomics.ElementAVITI
	case "PromethION24":
		return genomics.PromethION24
	default:
		return genomics.NovaSeq6000 // Default fallback
	}
}

func (fos *FingerprintObfuscationSystem) extractAnalysisFeatures(
	analyzer *FingerprintAnalyzer,
	sequences map[population.Population][]string,
	qualities [][]int,
) map[string]float64 {

	features := make(map[string]float64)

	// Extract features based on analyzer type
	switch analyzer.AnalysisType {
	case "Error_Pattern":
		features["substitution_rate"] = fos.calculateAverageSubstitutionRate(sequences)
		features["indel_rate"] = fos.calculateAverageIndelRate(sequences)
		features["context_dependency"] = fos.calculateContextDependency(sequences)
		features["quality_correlation"] = fos.calculateQualityCorrelation(sequences, qualities)
		features["homopolymer_patterns"] = fos.calculateHomopolymerPatterns(sequences)

	case "Quality_Distribution":
		features["mean_quality"] = fos.calculateMeanQuality(qualities)
		features["quality_variance"] = fos.calculateQualityVariance(qualities)
		features["quality_range"] = fos.calculateQualityRange(qualities)
		features["cycle_decay"] = fos.calculateCycleDecay(qualities)
		features["distribution_shape"] = fos.calculateDistributionShape(qualities)

	case "Technology_Signature":
		features["sequencing_chemistry"] = fos.calculateChemistrySignature(sequences)
		features["read_length_pattern"] = fos.calculateReadLengthPattern(sequences)
		features["error_profile"] = fos.calculateErrorProfile(sequences)
		features["quality_signature"] = fos.calculateQualitySignature(qualities)
		features["platform_artifacts"] = fos.calculatePlatformArtifacts(sequences)

	case "GC_Bias":
		features["gc_content_mean"] = fos.calculateMeanGCContent(sequences)
		features["gc_content_variance"] = fos.calculateGCContentVariance(sequences)
		features["gc_bias_strength"] = fos.calculateGCBiasStrength(sequences)
		features["window_patterns"] = fos.calculateWindowPatterns(sequences)
		features["population_deviation"] = fos.calculatePopulationGCDeviation(sequences)

	default:
		// Multi-modal analysis
		features["error_patterns"] = fos.calculateAverageSubstitutionRate(sequences)
		features["quality_distributions"] = fos.calculateMeanQuality(qualities)
		features["technology_signatures"] = fos.calculateChemistrySignature(sequences)
		features["gc_bias_patterns"] = fos.calculateMeanGCContent(sequences)
		features["read_length_analysis"] = fos.calculateReadLengthPattern(sequences)
		features["ensemble_consensus"] = fos.calculateEnsembleConsensus(sequences, qualities)
	}

	return features
}

func (fos *FingerprintObfuscationSystem) simulateDetectionAnalysis(
	analyzer *FingerprintAnalyzer, features map[string]float64,
) float64 {

	// Simulate detection analysis based on analyzer characteristics
	detectionScore := 0.0
	weightSum := 0.0

	for feature, value := range features {
		if weight, exists := analyzer.FeatureWeights[feature]; exists {
			// Normalize feature value and apply weight
			normalizedValue := math.Min(1.0, math.Max(0.0, value))
			detectionScore += normalizedValue * weight
			weightSum += weight
		}
	}

	if weightSum > 0 {
		detectionScore /= weightSum
	}

	// Apply analyzer sensitivity and specificity
	adjustedScore := detectionScore * analyzer.Sensitivity
	
	// Add some noise based on specificity
	noise := (1.0 - analyzer.Specificity) * fos.rng.Float64() * 0.1
	adjustedScore += noise

	return math.Min(1.0, math.Max(0.0, adjustedScore))
}

func (fos *FingerprintObfuscationSystem) simulatePlatformDetection(
	analyzer *FingerprintAnalyzer, features map[string]float64, detectionScore float64,
) genomics.SequencingPlatform {

	// Simple platform detection simulation
	if detectionScore < analyzer.DetectionThreshold {
		// Random platform if detection confidence is low
		platforms := []genomics.SequencingPlatform{
			genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq, genomics.PacBioRevio,
			genomics.DNBSeqT20x2, genomics.ElementAVITI, genomics.PromethION24,
		}
		return platforms[fos.rng.IntN(len(platforms))]
	}

	// Simulate platform detection based on features
	if features["error_patterns"] > 0.005 || features["quality_distributions"] < 20 {
		return genomics.PacBioRevio // Long-read characteristics
	}
	if features["gc_bias_patterns"] > 0.45 || features["technology_signatures"] > 0.8 {
		return genomics.NovaSeq6000 // Illumina characteristics
	}

	return genomics.ElementAVITI // Default fallback
}

// Feature calculation helper methods (simplified implementations)

// calculateAverageSubstitutionRate estimates per-base substitution density as the mean
// absolute deviation of per-sequence GC content from the human somatic average (42%).
// Sequences deviating further from 42% GC have experienced more context-specific
// mutational pressure, serving as a platform-neutral proxy for substitution load.
func (fos *FingerprintObfuscationSystem) calculateAverageSubstitutionRate(sequences map[population.Population][]string) float64 {
	totalDev := 0.0
	count := 0
	const expectedGC = 0.42
	for _, seqs := range sequences {
		for _, seq := range seqs {
			if len(seq) == 0 {
				continue
			}
			gc := 0
			for i := 0; i < len(seq); i++ {
				b := seq[i] | 0x20
				if b == 'g' || b == 'c' {
					gc++
				}
			}
			gcFrac := float64(gc) / float64(len(seq))
			totalDev += math.Abs(gcFrac - expectedGC)
			count++
		}
	}
	if count == 0 {
		return 0.001
	}
	// Scale: mean GC deviation of ~0.04 ≈ substitution rate of 0.001 (1/40 scaling).
	return totalDev / float64(count) / 40.0
}

// calculateAverageIndelRate returns the fraction of bases that fall inside homopolymer
// runs of length ≥ 5. Such runs are the primary source of insertion/deletion errors
// on Illumina (slippage) and are longer on PacBio/Nanopore (STR expansion errors).
func (fos *FingerprintObfuscationSystem) calculateAverageIndelRate(sequences map[population.Population][]string) float64 {
	totalBases := 0
	indelBases := 0
	const minRun = 5
	for _, seqs := range sequences {
		for _, seq := range seqs {
			totalBases += len(seq)
			i := 0
			for i < len(seq) {
				j := i + 1
				for j < len(seq) && (seq[j]|0x20) == (seq[i]|0x20) {
					j++
				}
				if j-i >= minRun {
					indelBases += j - i
				}
				i = j
			}
		}
	}
	if totalBases == 0 {
		return 0.0001
	}
	return float64(indelBases) / float64(totalBases)
}

// calculateContextDependency computes the CpG observed/expected ratio across sequences.
// Human genomic CpG is depleted (~0.25 O/E) due to methylation-driven deamination;
// synthetic sequences without this model approach 1.0. Platform errors can create
// spurious CpG sites.
func (fos *FingerprintObfuscationSystem) calculateContextDependency(sequences map[population.Population][]string) float64 {
	var nC, nG, nCG, nTotal int64
	for _, seqs := range sequences {
		for _, seq := range seqs {
			n := len(seq)
			nTotal += int64(n)
			for i := 0; i < n; i++ {
				b := seq[i] | 0x20
				if b == 'c' {
					nC++
				} else if b == 'g' {
					nG++
				}
				if i > 0 {
					prev := seq[i-1] | 0x20
					if prev == 'c' && b == 'g' {
						nCG++
					}
				}
			}
		}
	}
	if nTotal < 2 || nC == 0 || nG == 0 {
		return 1.0
	}
	// O/E = (n_CG * N) / (n_C * n_G)  where N = total dinucleotide positions = nTotal-seqCount
	// Approximate N as nTotal - number-of-sequences (first position of each read excluded).
	seqCount := int64(0)
	for _, seqs := range sequences {
		seqCount += int64(len(seqs))
	}
	dinucN := nTotal - seqCount
	if dinucN <= 0 {
		return 1.0
	}
	oe := float64(nCG) * float64(dinucN) / (float64(nC) * float64(nG))
	return oe
}

// calculateQualityCorrelation computes the Pearson correlation between read position
// (0-indexed within the read) and Phred quality score. Illumina shows a characteristic
// negative correlation (quality decline toward the 3' end); PacBio/Nanopore differ.
func (fos *FingerprintObfuscationSystem) calculateQualityCorrelation(sequences map[population.Population][]string, qualities [][]int) float64 {
	if len(qualities) == 0 {
		return -0.5 // typical Illumina default
	}
	var sumX, sumY, sumXY, sumX2, sumY2 float64
	var n float64
	for _, quals := range qualities {
		for pos, q := range quals {
			x := float64(pos)
			y := float64(q)
			sumX += x
			sumY += y
			sumXY += x * y
			sumX2 += x * x
			sumY2 += y * y
			n++
		}
	}
	if n < 2 {
		return -0.5
	}
	// Pearson r = (n*ΣXY - ΣX*ΣY) / sqrt((n*ΣX²-(ΣX)²)*(n*ΣY²-(ΣY)²))
	num := n*sumXY - sumX*sumY
	denX := n*sumX2 - sumX*sumX
	denY := n*sumY2 - sumY*sumY
	if denX <= 0 || denY <= 0 {
		return 0.0
	}
	r := num / math.Sqrt(denX*denY)
	if r < -1 {
		r = -1
	} else if r > 1 {
		r = 1
	}
	return r
}

// calculateHomopolymerPatterns returns the fraction of bases in each sequence that
// belong to homopolymer runs of length ≥ 4, averaged across all sequences. Long
// homopolymers are platform-characteristic (PacBio/Nanopore struggle with them).
func (fos *FingerprintObfuscationSystem) calculateHomopolymerPatterns(sequences map[population.Population][]string) float64 {
	totalBases := 0
	hpBases := 0
	const minRun = 4
	for _, seqs := range sequences {
		for _, seq := range seqs {
			totalBases += len(seq)
			i := 0
			for i < len(seq) {
				j := i + 1
				for j < len(seq) && (seq[j]|0x20) == (seq[i]|0x20) {
					j++
				}
				if j-i >= minRun {
					hpBases += j - i
				}
				i = j
			}
		}
	}
	if totalBases == 0 {
		return 0.03
	}
	return float64(hpBases) / float64(totalBases)
}

func (fos *FingerprintObfuscationSystem) calculateMeanQuality(qualities [][]int) float64 {
	if len(qualities) == 0 {
		return 30.0
	}
	
	totalQuality := 0.0
	totalBases := 0
	
	for _, qualArray := range qualities {
		for _, q := range qualArray {
			totalQuality += float64(q)
			totalBases++
		}
	}
	
	if totalBases == 0 {
		return 30.0
	}
	
	return totalQuality / float64(totalBases)
}

func (fos *FingerprintObfuscationSystem) calculateQualityVariance(qualities [][]int) float64 {
	mean := fos.calculateMeanQuality(qualities)
	
	totalVariance := 0.0
	totalBases := 0
	
	for _, qualArray := range qualities {
		for _, q := range qualArray {
			diff := float64(q) - mean
			totalVariance += diff * diff
			totalBases++
		}
	}
	
	if totalBases == 0 {
		return 25.0
	}
	
	return totalVariance / float64(totalBases)
}

func (fos *FingerprintObfuscationSystem) calculateQualityRange(qualities [][]int) float64 {
	if len(qualities) == 0 {
		return 40.0
	}
	
	minQ := 100
	maxQ := 0
	
	for _, qualArray := range qualities {
		for _, q := range qualArray {
			if q < minQ {
				minQ = q
			}
			if q > maxQ {
				maxQ = q
			}
		}
	}
	
	return float64(maxQ - minQ)
}

// calculateCycleDecay measures per-position quality decay from the first 20% of
// cycle positions to the last 20%.  Returns a value in [0, 1]: 0 = no decay,
// 1 = complete quality collapse.
func (fos *FingerprintObfuscationSystem) calculateCycleDecay(qualities [][]int) float64 {
	if len(qualities) == 0 {
		return 0.05
	}
	// Collect all read lengths and find the max position count.
	maxLen := 0
	for _, q := range qualities {
		if len(q) > maxLen {
			maxLen = len(q)
		}
	}
	if maxLen == 0 {
		return 0.05
	}
	earlyEnd := maxLen / 5
	if earlyEnd < 1 {
		earlyEnd = 1
	}
	lateStart := maxLen - maxLen/5
	if lateStart <= earlyEnd {
		lateStart = earlyEnd
	}

	earlySum, earlyN := 0.0, 0
	lateSum, lateN := 0.0, 0
	for _, q := range qualities {
		for i, v := range q {
			if i < earlyEnd {
				earlySum += float64(v)
				earlyN++
			}
			if i >= lateStart {
				lateSum += float64(v)
				lateN++
			}
		}
	}
	if earlyN == 0 || lateN == 0 {
		return 0.05
	}
	earlyMean := earlySum / float64(earlyN)
	lateMean := lateSum / float64(lateN)
	if earlyMean <= 0 {
		return 0.0
	}
	decay := (earlyMean - lateMean) / earlyMean
	if decay < 0 {
		return 0.0
	}
	if decay > 1 {
		return 1.0
	}
	return decay
}

// calculateDistributionShape returns the fraction of quality scores that fall
// outside the interval [mean − σ, mean + σ].  A Gaussian distribution has
// ~32% of values outside this interval; values further from Gaussian indicate
// platform-specific quality score clumping.
func (fos *FingerprintObfuscationSystem) calculateDistributionShape(qualities [][]int) float64 {
	mean := fos.calculateMeanQuality(qualities)
	variance := fos.calculateQualityVariance(qualities)
	std := math.Sqrt(variance)

	total, outside := 0, 0
	for _, q := range qualities {
		for _, v := range q {
			total++
			fv := float64(v)
			if fv < mean-std || fv > mean+std {
				outside++
			}
		}
	}
	if total == 0 {
		return 0.32 // Gaussian baseline
	}
	return float64(outside) / float64(total)
}

// calculateChemistrySignature returns the fraction of reads whose GC content
// falls outside the chemically "normal" range [30%, 70%].  High values
// indicate chemistry artefacts such as PCR duplicates or adapter bias.
func (fos *FingerprintObfuscationSystem) calculateChemistrySignature(sequences map[population.Population][]string) float64 {
	extreme, total := 0, 0
	for _, seqs := range sequences {
		for _, seq := range seqs {
			total++
			gc := fos.calculateGCContent(seq)
			if gc < 0.30 || gc > 0.70 {
				extreme++
			}
		}
	}
	if total == 0 {
		return 0.05
	}
	return float64(extreme) / float64(total)
}

// calculateReadLengthPattern returns the coefficient of variation (std/mean)
// of read lengths across all sequences.  A value near 0 means uniform read
// lengths; higher values indicate variable-length reads (e.g. Nanopore).
func (fos *FingerprintObfuscationSystem) calculateReadLengthPattern(sequences map[population.Population][]string) float64 {
	var lengths []float64
	for _, seqs := range sequences {
		for _, seq := range seqs {
			lengths = append(lengths, float64(len(seq)))
		}
	}
	if len(lengths) < 2 {
		return 0.0
	}
	sum := 0.0
	for _, l := range lengths {
		sum += l
	}
	mean := sum / float64(len(lengths))
	if mean == 0 {
		return 0.0
	}
	varSum := 0.0
	for _, l := range lengths {
		d := l - mean
		varSum += d * d
	}
	std := math.Sqrt(varSum / float64(len(lengths)-1))
	cv := std / mean
	if cv > 1.0 {
		return 1.0
	}
	return cv
}

// calculateErrorProfile returns the fraction of bases in homopolymer runs of
// length ≥ 3.  Short homopolymers (length 3–4) are common in all platforms
// and represent a baseline error-prone context; longer runs increase with
// read length and are more common in long-read technologies.
func (fos *FingerprintObfuscationSystem) calculateErrorProfile(sequences map[population.Population][]string) float64 {
	const minRun = 3
	totalBases, hpBases := 0, 0
	for _, seqs := range sequences {
		for _, seq := range seqs {
			totalBases += len(seq)
			i := 0
			for i < len(seq) {
				j := i + 1
				for j < len(seq) && (seq[j]|0x20) == (seq[i]|0x20) {
					j++
				}
				if j-i >= minRun {
					hpBases += j - i
				}
				i = j
			}
		}
	}
	if totalBases == 0 {
		return 0.05
	}
	return float64(hpBases) / float64(totalBases)
}

func (fos *FingerprintObfuscationSystem) calculateQualitySignature(qualities [][]int) float64 {
	mean := fos.calculateMeanQuality(qualities)
	return mean / 40.0 // Normalize to 0-1
}

// calculatePlatformArtifacts returns the fraction of reads with severely
// aberrant GC content (< 20% or > 80%).  These extreme outliers are
// typically artefactual in human WGS data and serve as a platform quality
// signal distinct from the broader chemistry-signature threshold.
func (fos *FingerprintObfuscationSystem) calculatePlatformArtifacts(sequences map[population.Population][]string) float64 {
	extreme, total := 0, 0
	for _, seqs := range sequences {
		for _, seq := range seqs {
			total++
			gc := fos.calculateGCContent(seq)
			if gc < 0.20 || gc > 0.80 {
				extreme++
			}
		}
	}
	if total == 0 {
		return 0.02
	}
	return float64(extreme) / float64(total)
}

func (fos *FingerprintObfuscationSystem) calculateMeanGCContent(sequences map[population.Population][]string) float64 {
	totalGC := 0.0
	totalSeqs := 0
	
	for _, seqs := range sequences {
		for _, seq := range seqs {
			gc := fos.calculateGCContent(seq)
			totalGC += gc
			totalSeqs++
		}
	}
	
	if totalSeqs == 0 {
		return 0.41
	}
	
	return totalGC / float64(totalSeqs)
}

func (fos *FingerprintObfuscationSystem) calculateGCContentVariance(sequences map[population.Population][]string) float64 {
	mean := fos.calculateMeanGCContent(sequences)
	
	totalVariance := 0.0
	totalSeqs := 0
	
	for _, seqs := range sequences {
		for _, seq := range seqs {
			gc := fos.calculateGCContent(seq)
			diff := gc - mean
			totalVariance += diff * diff
			totalSeqs++
		}
	}
	
	if totalSeqs == 0 {
		return 0.01
	}
	
	return totalVariance / float64(totalSeqs)
}

// calculateGCBiasStrength returns the std dev of per-read GC content across
// all sequences, normalised by 0.20 (a typical upper bound for a sequencer
// with strong GC bias).  Result is clamped to [0, 1].
func (fos *FingerprintObfuscationSystem) calculateGCBiasStrength(sequences map[population.Population][]string) float64 {
	var gcValues []float64
	for _, seqs := range sequences {
		for _, seq := range seqs {
			gcValues = append(gcValues, fos.calculateGCContent(seq))
		}
	}
	if len(gcValues) < 2 {
		return 0.0
	}
	sum := 0.0
	for _, v := range gcValues {
		sum += v
	}
	mean := sum / float64(len(gcValues))
	varSum := 0.0
	for _, v := range gcValues {
		d := v - mean
		varSum += d * d
	}
	std := math.Sqrt(varSum / float64(len(gcValues)-1))
	normalised := std / 0.20
	if normalised > 1.0 {
		return 1.0
	}
	return normalised
}

// calculateWindowPatterns returns the mean local GC-content variance measured
// in 20-base sliding windows (step = 10) across all sequences.  High variance
// indicates non-uniform GC composition — a known fingerprinting signal.
func (fos *FingerprintObfuscationSystem) calculateWindowPatterns(sequences map[population.Population][]string) float64 {
	const windowSize = 20
	const stepSize = 10

	totalVariance := 0.0
	seqCount := 0

	for _, seqs := range sequences {
		for _, seq := range seqs {
			seqCount++
			if len(seq) < windowSize {
				totalVariance += 0.0
				continue
			}
			var windowGCs []float64
			for start := 0; start+windowSize <= len(seq); start += stepSize {
				gc := fos.calculateGCContent(seq[start : start+windowSize])
				windowGCs = append(windowGCs, gc)
			}
			if len(windowGCs) < 2 {
				continue
			}
			wSum := 0.0
			for _, v := range windowGCs {
				wSum += v
			}
			wMean := wSum / float64(len(windowGCs))
			wVarSum := 0.0
			for _, v := range windowGCs {
				d := v - wMean
				wVarSum += d * d
			}
			totalVariance += wVarSum / float64(len(windowGCs)-1)
		}
	}
	if seqCount == 0 {
		return 0.0
	}
	return totalVariance / float64(seqCount)
}

// calculatePopulationGCDeviation returns the std dev of per-population mean
// GC contents.  A value near 0 means all populations have similar GC
// distributions; larger values indicate population-specific biases.
func (fos *FingerprintObfuscationSystem) calculatePopulationGCDeviation(sequences map[population.Population][]string) float64 {
	var popMeans []float64
	for _, seqs := range sequences {
		if len(seqs) == 0 {
			continue
		}
		sum := 0.0
		for _, seq := range seqs {
			sum += fos.calculateGCContent(seq)
		}
		popMeans = append(popMeans, sum/float64(len(seqs)))
	}
	if len(popMeans) < 2 {
		return 0.0
	}
	sum := 0.0
	for _, v := range popMeans {
		sum += v
	}
	mean := sum / float64(len(popMeans))
	varSum := 0.0
	for _, v := range popMeans {
		d := v - mean
		varSum += d * d
	}
	return math.Sqrt(varSum / float64(len(popMeans)-1))
}

// calculateEnsembleConsensus returns the Pearson correlation between per-read
// GC content and mean per-read quality, mapped from [-1, 1] to [0, 1] via
// (r + 1) / 2.  A value near 0.5 indicates no correlation; lower/higher
// values indicate GC-quality anti-/correlation (a platform artefact).
func (fos *FingerprintObfuscationSystem) calculateEnsembleConsensus(sequences map[population.Population][]string, qualities [][]int) float64 {
	// Build flat list of (gc, meanQuality) pairs.
	// Pair sequence index with a flattened quality-read index.
	var allSeqs []string
	for _, seqs := range sequences {
		allSeqs = append(allSeqs, seqs...)
	}
	n := len(allSeqs)
	if n == 0 || len(qualities) == 0 {
		return 0.5
	}
	// Pair by index up to min(n, len(qualities)).
	pairs := n
	if len(qualities) < pairs {
		pairs = len(qualities)
	}

	var sumX, sumY, sumXY, sumX2, sumY2 float64
	count := 0
	for i := 0; i < pairs; i++ {
		gc := fos.calculateGCContent(allSeqs[i])
		if len(qualities[i]) == 0 {
			continue
		}
		qSum := 0.0
		for _, v := range qualities[i] {
			qSum += float64(v)
		}
		meanQ := qSum / float64(len(qualities[i]))
		sumX += gc
		sumY += meanQ
		sumXY += gc * meanQ
		sumX2 += gc * gc
		sumY2 += meanQ * meanQ
		count++
	}
	if count < 2 {
		return 0.5
	}
	fcount := float64(count)
	num := fcount*sumXY - sumX*sumY
	denX := fcount*sumX2 - sumX*sumX
	denY := fcount*sumY2 - sumY*sumY
	if denX <= 0 || denY <= 0 {
		return 0.5
	}
	r := num / math.Sqrt(denX*denY)
	if r < -1 {
		r = -1
	} else if r > 1 {
		r = 1
	}
	return (r + 1.0) / 2.0
}

// Getter methods
func (fos *FingerprintObfuscationSystem) GetObfuscationStrategies() map[string]*ObfuscationStrategy {
	return fos.obfuscationStrategies
}

func (fos *FingerprintObfuscationSystem) GetFingerprintAnalyzers() map[string]*FingerprintAnalyzer {
	return fos.fingerprintAnalyzers
}

func (fos *FingerprintObfuscationSystem) GetCrossPlatformModels() map[string]*CrossPlatformModel {
	return fos.crossPlatformModels
}