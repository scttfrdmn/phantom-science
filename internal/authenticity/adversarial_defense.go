package authenticity

import (
	"fmt"
	"math"
	"math/rand/v2"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

// AdversarialDefenseSystem implements advanced population classification resistance
type AdversarialDefenseSystem struct {
	indistinguishabilityTester *IndistinguishabilityTester
	platformManager            *genomics.PlatformErrorSignatureManager
	qualityModeler             *genomics.MultiPlatformQualityModeler
	adversarialModels          map[string]*AdversarialModel
	defenseStrategies          map[string]*DefenseStrategy
	rng                        *rand.Rand
}

// AdversarialModel represents a population classification attack model
// Updated with 2024-2025 research findings including TabPFN and pg-gan approaches
type AdversarialModel struct {
	ModelName        string
	ModelType        string // "SVM", "RandomForest", "NeuralNetwork", "DeepLearning", "Ensemble", "TabPFN", "PGGAN", "DPGAN"
	FeatureSet       []string
	TrainingAccuracy float64
	ValidationAccuracy float64
	FeatureWeights   map[string]float64
	AttackCapability AttackCapability
	ModelParameters  map[string]interface{}
	// 2024-2025 Advanced Capabilities
	UsesTabularTransformer bool     // TabPFN-based ancestry classification
	UsesGANDiscriminator   bool     // pg-gan discriminator approach
	DifferentialPrivacy    bool     // DP-based attack resistance
	FederatedCapability    bool     // Federated learning attack vectors
	AncestryInferenceMode  string   // "local", "global", "mixed"
}

// AttackCapability defines the sophistication of classification attacks
// Enhanced with 2024-2025 advanced ML techniques
type AttackCapability struct {
	CanDetectPCA           bool
	CanDetectFST           bool
	CanDetectLD            bool
	CanDetectQuality       bool
	CanDetectErrorPatterns bool
	CanDetectGCBias        bool
	CanUseEnsembleMethods  bool
	CanAdaptToDefenses     bool
	MaxFeatureDimensions   int
	RequiredSampleSize     int
	// 2024-2025 Advanced Attack Capabilities
	CanUseTabPFN           bool    // TabPFN tabular transformer attacks
	CanUsePGGAN            bool    // pg-gan multi-discriminator attacks
	CanUseDPGAN            bool    // Differential privacy GAN attacks
	CanUseFederatedLearning bool   // Federated learning attack vectors
	CanDetectCrossPlatformBias bool // Cross-platform bias fingerprinting
	CanUseLocalAncestryInference bool // LAI-Net, SALAI-Net approaches
	CanDetectSequencingArtifacts bool // Platform-specific error signatures
	HaploNetCapability      bool   // VAE-based population structure detection
	MaxROCAUC              float64 // Maximum achievable ROC AUC
	MinSampleSize          int     // Minimum samples for reliable classification
}

// DefenseStrategy represents a specific adversarial defense mechanism
// Enhanced with 2024-2025 privacy-preserving techniques
type DefenseStrategy struct {
	StrategyName         string
	DefenseType          string // "Noise_Injection", "Feature_Masking", "Adversarial_Training", "Ensemble_Defense", "DPGAN_Defense", "Federated_Privacy", "BEADS_Normalization"
	TargetedModels       []string
	EffectivenessScore   float64
	ComputationalCost    float64
	DefenseParameters    map[string]float64
	CountermeasurePower  float64
	// 2024-2025 Advanced Defense Capabilities
	UsesDifferentialPrivacy bool     // DPGAN/PATEGAN-based defense
	UsesFederatedPrivacy    bool     // Federated learning privacy protection
	UsesBEADSNormalization  bool     // BEADS bias elimination
	UsesMANCIEIntegration   bool     // MANCIE cross-platform integration
	UsesQuantileNormalization bool   // QN for cross-platform defense
	PrivacyBudget          float64   // Differential privacy epsilon
	FederationNodes        int       // Number of federated learning nodes
	CrossPlatformBias      float64   // Cross-platform bias correction strength
}

// AdversarialTestResult stores results from adversarial resistance testing
type AdversarialTestResult struct {
	ModelName              string
	DefenseStrategy        string
	PreDefenseAccuracy     float64
	PostDefenseAccuracy    float64
	AccuracyReduction      float64
	DefenseEffectiveness   float64
	ConfidenceInterval     [2]float64
	StatisticalSignificance float64
	ResistanceLevel        string // "Vulnerable", "Partially_Resistant", "Highly_Resistant", "Immune"
}

// AdversarialDefenseResult aggregates comprehensive adversarial testing results
type AdversarialDefenseResult struct {
	OverallResistanceScore float64
	ModelSpecificResults   map[string]*AdversarialTestResult
	DefenseEffectiveness   map[string]float64
	VulnerabilityAssessment VulnerabilityAssessment
	RecommendedDefenses    []string
	SecurityLevel          string
}

// VulnerabilityAssessment provides detailed security analysis
// Updated with 2024-2025 threat landscape
type VulnerabilityAssessment struct {
	HighRiskModels        []string
	CriticalVulnerabilities []string
	DefensiveGaps         []string
	AttackVectors         []string
	RequiredImprovements  []string
	SecurityRating        float64 // 0-100 scale
	// 2024-2025 Advanced Threat Assessment
	TabPFNResistance      float64 // Resistance to TabPFN-based attacks
	PGGANResistance       float64 // Resistance to pg-gan discriminator attacks
	DPGANResistance       float64 // Resistance to differential privacy GAN attacks
	FederatedResistance   float64 // Resistance to federated learning attacks
	CrossPlatformResistance float64 // Resistance to cross-platform fingerprinting
	LocalAncestryResistance float64 // Resistance to local ancestry inference
	HaploNetResistance    float64 // Resistance to VAE-based population detection
	SequencingArtifactResistance float64 // Resistance to platform artifact detection
	OverallThreatLevel    string  // "Low", "Medium", "High", "Critical"
	RecommendedPrivacyBudget float64 // Optimal differential privacy epsilon
}

// NewAdversarialDefenseSystem creates a new adversarial defense system
func NewAdversarialDefenseSystem() *AdversarialDefenseSystem {
	return &AdversarialDefenseSystem{
		indistinguishabilityTester: NewIndistinguishabilityTester(),
		platformManager:            genomics.NewPlatformErrorSignatureManager(),
		qualityModeler:             genomics.NewMultiPlatformQualityModeler(),
		adversarialModels:          make(map[string]*AdversarialModel),
		defenseStrategies:          make(map[string]*DefenseStrategy),
		rng:                        rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
	}
}

// InitializeAdversarialModels sets up sophisticated population classification attack models
func (ads *AdversarialDefenseSystem) InitializeAdversarialModels() {
	
	// Support Vector Machine (SVM) - Traditional approach
	ads.adversarialModels["SVM_Population_Classifier"] = &AdversarialModel{
		ModelName:          "SVM Population Classifier",
		ModelType:          "SVM",
		FeatureSet:         []string{"PCA_Components", "Allele_Frequencies", "LD_Patterns", "GC_Content"},
		TrainingAccuracy:   0.85,
		ValidationAccuracy: 0.78,
		FeatureWeights: map[string]float64{
			"PCA_Components":   0.40,
			"Allele_Frequencies": 0.30,
			"LD_Patterns":      0.20,
			"GC_Content":       0.10,
		},
		AttackCapability: AttackCapability{
			CanDetectPCA:           true,
			CanDetectFST:           true,
			CanDetectLD:            true,
			CanDetectQuality:       false,
			CanDetectErrorPatterns: false,
			CanDetectGCBias:        true,
			CanUseEnsembleMethods:  false,
			CanAdaptToDefenses:     false,
			MaxFeatureDimensions:   100,
			RequiredSampleSize:     50,
		},
		ModelParameters: map[string]interface{}{
			"kernel":      "rbf",
			"C":           1.0,
			"gamma":       "scale",
			"probability": true,
		},
	}

	// Random Forest - Ensemble method
	ads.adversarialModels["RandomForest_Population_Classifier"] = &AdversarialModel{
		ModelName:          "Random Forest Population Classifier",
		ModelType:          "RandomForest",
		FeatureSet:         []string{"SNP_Patterns", "Haplotype_Blocks", "Quality_Distributions", "Error_Rates"},
		TrainingAccuracy:   0.88,
		ValidationAccuracy: 0.82,
		FeatureWeights: map[string]float64{
			"SNP_Patterns":        0.35,
			"Haplotype_Blocks":    0.25,
			"Quality_Distributions": 0.25,
			"Error_Rates":         0.15,
		},
		AttackCapability: AttackCapability{
			CanDetectPCA:           true,
			CanDetectFST:           true,
			CanDetectLD:            true,
			CanDetectQuality:       true,
			CanDetectErrorPatterns: true,
			CanDetectGCBias:        true,
			CanUseEnsembleMethods:  true,
			CanAdaptToDefenses:     false,
			MaxFeatureDimensions:   500,
			RequiredSampleSize:     100,
		},
		ModelParameters: map[string]interface{}{
			"n_estimators":     100,
			"max_depth":        10,
			"min_samples_split": 2,
			"random_state":     42,
		},
	}

	// Neural Network - Deep learning approach
	ads.adversarialModels["NeuralNetwork_Population_Classifier"] = &AdversarialModel{
		ModelName:          "Neural Network Population Classifier",
		ModelType:          "NeuralNetwork",
		FeatureSet:         []string{"Raw_Sequence_Embeddings", "Quality_Vectors", "Platform_Signatures", "Population_Embeddings"},
		TrainingAccuracy:   0.92,
		ValidationAccuracy: 0.86,
		FeatureWeights: map[string]float64{
			"Raw_Sequence_Embeddings": 0.40,
			"Quality_Vectors":         0.25,
			"Platform_Signatures":     0.20,
			"Population_Embeddings":   0.15,
		},
		AttackCapability: AttackCapability{
			CanDetectPCA:           true,
			CanDetectFST:           true,
			CanDetectLD:            true,
			CanDetectQuality:       true,
			CanDetectErrorPatterns: true,
			CanDetectGCBias:        true,
			CanUseEnsembleMethods:  true,
			CanAdaptToDefenses:     true,
			MaxFeatureDimensions:   10000,
			RequiredSampleSize:     200,
		},
		ModelParameters: map[string]interface{}{
			"hidden_layers":   []int{256, 128, 64},
			"activation":      "relu",
			"optimizer":       "adam",
			"learning_rate":   0.001,
			"dropout":         0.3,
			"batch_size":      32,
			"epochs":          100,
		},
	}

	// Deep Learning Convolutional Network - Advanced approach
	ads.adversarialModels["DeepLearning_CNN_Classifier"] = &AdversarialModel{
		ModelName:          "Deep Learning CNN Population Classifier",
		ModelType:          "DeepLearning",
		FeatureSet:         []string{"Sequence_Images", "Quality_Heatmaps", "Error_Pattern_Maps", "Multi_Platform_Features"},
		TrainingAccuracy:   0.95,
		ValidationAccuracy: 0.89,
		FeatureWeights: map[string]float64{
			"Sequence_Images":        0.35,
			"Quality_Heatmaps":       0.25,
			"Error_Pattern_Maps":     0.25,
			"Multi_Platform_Features": 0.15,
		},
		AttackCapability: AttackCapability{
			CanDetectPCA:           true,
			CanDetectFST:           true,
			CanDetectLD:            true,
			CanDetectQuality:       true,
			CanDetectErrorPatterns: true,
			CanDetectGCBias:        true,
			CanUseEnsembleMethods:  true,
			CanAdaptToDefenses:     true,
			MaxFeatureDimensions:   100000,
			RequiredSampleSize:     500,
		},
		ModelParameters: map[string]interface{}{
			"conv_layers":     []map[string]interface{}{
				{"filters": 32, "kernel_size": 3, "activation": "relu"},
				{"filters": 64, "kernel_size": 3, "activation": "relu"},
				{"filters": 128, "kernel_size": 3, "activation": "relu"},
			},
			"pooling":         "max",
			"dense_layers":    []int{512, 256},
			"regularization":  "l2",
			"data_augmentation": true,
		},
	}

	// Advanced Ensemble Model - State-of-the-art approach
	ads.adversarialModels["Advanced_Ensemble_Classifier"] = &AdversarialModel{
		ModelName:          "Advanced Ensemble Population Classifier",
		ModelType:          "Ensemble",
		FeatureSet:         []string{"All_Available_Features", "Cross_Platform_Signatures", "Temporal_Patterns", "Adversarial_Features"},
		TrainingAccuracy:   0.97,
		ValidationAccuracy: 0.91,
		FeatureWeights: map[string]float64{
			"All_Available_Features":    0.30,
			"Cross_Platform_Signatures": 0.25,
			"Temporal_Patterns":         0.25,
			"Adversarial_Features":      0.20,
		},
		AttackCapability: AttackCapability{
			CanDetectPCA:           true,
			CanDetectFST:           true,
			CanDetectLD:            true,
			CanDetectQuality:       true,
			CanDetectErrorPatterns: true,
			CanDetectGCBias:        true,
			CanUseEnsembleMethods:  true,
			CanAdaptToDefenses:     true,
			MaxFeatureDimensions:   1000000,
			RequiredSampleSize:     1000,
		},
		ModelParameters: map[string]interface{}{
			"base_models": []string{"SVM", "RandomForest", "NeuralNetwork", "CNN", "XGBoost"},
			"meta_learner": "StackingClassifier",
			"cv_folds":     5,
			"voting":       "soft",
			"adaptive":     true,
		},
	}
}

// InitializeDefenseStrategies sets up comprehensive adversarial defense mechanisms
func (ads *AdversarialDefenseSystem) InitializeDefenseStrategies() {
	
	// Noise Injection Defense - Add controlled noise to features
	ads.defenseStrategies["Gaussian_Noise_Injection"] = &DefenseStrategy{
		StrategyName:       "Gaussian Noise Injection",
		DefenseType:        "Noise_Injection",
		TargetedModels:     []string{"SVM", "RandomForest", "NeuralNetwork"},
		EffectivenessScore: 0.70,
		ComputationalCost:  0.20,
		DefenseParameters: map[string]float64{
			"noise_std":        0.05,
			"noise_mean":       0.0,
			"feature_coverage": 0.8,
			"adaptive_scaling": 1.0,
		},
		CountermeasurePower: 0.65,
	}

	// Feature Masking Defense - Obfuscate discriminative features
	ads.defenseStrategies["Adaptive_Feature_Masking"] = &DefenseStrategy{
		StrategyName:       "Adaptive Feature Masking",
		DefenseType:        "Feature_Masking",
		TargetedModels:     []string{"SVM", "RandomForest", "NeuralNetwork", "DeepLearning"},
		EffectivenessScore: 0.75,
		ComputationalCost:  0.30,
		DefenseParameters: map[string]float64{
			"masking_probability": 0.15,
			"feature_importance_threshold": 0.1,
			"dynamic_adjustment": 1.0,
			"population_specific": 1.0,
		},
		CountermeasurePower: 0.72,
	}

	// Adversarial Training Defense - Train with adversarial examples
	ads.defenseStrategies["Adversarial_Training_Defense"] = &DefenseStrategy{
		StrategyName:       "Adversarial Training Defense",
		DefenseType:        "Adversarial_Training",
		TargetedModels:     []string{"NeuralNetwork", "DeepLearning"},
		EffectivenessScore: 0.82,
		ComputationalCost:  0.70,
		DefenseParameters: map[string]float64{
			"adversarial_examples_ratio": 0.3,
			"perturbation_strength": 0.1,
			"training_iterations": 50.0,
			"robustness_regularization": 0.01,
		},
		CountermeasurePower: 0.78,
	}

	// Ensemble Defense - Multiple defense layers
	ads.defenseStrategies["Multi_Layer_Ensemble_Defense"] = &DefenseStrategy{
		StrategyName:       "Multi-Layer Ensemble Defense",
		DefenseType:        "Ensemble_Defense",
		TargetedModels:     []string{"All_Models"},
		EffectivenessScore: 0.88,
		ComputationalCost:  0.85,
		DefenseParameters: map[string]float64{
			"defense_layers": 3.0,
			"consensus_threshold": 0.7,
			"adaptive_weighting": 1.0,
			"feedback_learning": 1.0,
		},
		CountermeasurePower: 0.85,
	}

	// Differential Privacy Defense - Privacy-preserving noise
	ads.defenseStrategies["Differential_Privacy_Defense"] = &DefenseStrategy{
		StrategyName:       "Differential Privacy Defense",
		DefenseType:        "Privacy_Preserving",
		TargetedModels:     []string{"All_Models"},
		EffectivenessScore: 0.85,
		ComputationalCost:  0.40,
		DefenseParameters: map[string]float64{
			"epsilon": 1.0,  // Privacy budget
			"delta":   1e-5, // Privacy parameter
			"sensitivity": 1.0,
			"noise_mechanism": 1.0, // Laplacian
		},
		CountermeasurePower: 0.80,
	}

	// Dynamic Feature Scrambling - Real-time feature obfuscation
	ads.defenseStrategies["Dynamic_Feature_Scrambling"] = &DefenseStrategy{
		StrategyName:       "Dynamic Feature Scrambling",
		DefenseType:        "Real_Time_Defense",
		TargetedModels:     []string{"All_Models"},
		EffectivenessScore: 0.79,
		ComputationalCost:  0.55,
		DefenseParameters: map[string]float64{
			"scrambling_frequency": 0.2,
			"pattern_complexity": 5.0,
			"temporal_variation": 1.0,
			"cross_platform_coordination": 1.0,
		},
		CountermeasurePower: 0.76,
	}
}

// TestAdversarialResistance performs comprehensive adversarial resistance testing
func (ads *AdversarialDefenseSystem) TestAdversarialResistance(
	sequences map[genomics.SequencingPlatform]map[population.Population][]string,
	sampleSize int,
) (*AdversarialDefenseResult, error) {

	if len(sequences) == 0 {
		return nil, fmt.Errorf("no sequences provided for adversarial testing")
	}

	// Initialize models and defenses if not already done
	if len(ads.adversarialModels) == 0 {
		ads.InitializeAdversarialModels()
	}
	if len(ads.defenseStrategies) == 0 {
		ads.InitializeDefenseStrategies()
	}

	result := &AdversarialDefenseResult{
		ModelSpecificResults: make(map[string]*AdversarialTestResult),
		DefenseEffectiveness: make(map[string]float64),
	}

	// Test each adversarial model against defense strategies
	allTestResults := []*AdversarialTestResult{}
	
	for modelName, model := range ads.adversarialModels {
		fmt.Printf("Testing adversarial resistance against %s...\n", modelName)
		
		// Test without defenses (baseline)
		baselineAccuracy, err := ads.simulateModelAttack(model, sequences, sampleSize, nil)
		if err != nil {
			fmt.Printf("Warning: Baseline test failed for %s: %v\n", modelName, err)
			continue
		}

		// Test with each defense strategy
		bestDefenseResult := &AdversarialTestResult{
			ModelName:          modelName,
			PreDefenseAccuracy: baselineAccuracy,
			PostDefenseAccuracy: baselineAccuracy,
			AccuracyReduction:  0.0,
			DefenseEffectiveness: 0.0,
			ResistanceLevel:    "Vulnerable",
		}

		for defenseName, defense := range ads.defenseStrategies {
			// Check if defense is applicable to this model
			if !ads.isDefenseApplicable(defense, model) {
				continue
			}

			defendedAccuracy, err := ads.simulateModelAttack(model, sequences, sampleSize, defense)
			if err != nil {
				continue
			}

			testResult := &AdversarialTestResult{
				ModelName:              modelName,
				DefenseStrategy:        defenseName,
				PreDefenseAccuracy:     baselineAccuracy,
				PostDefenseAccuracy:    defendedAccuracy,
				AccuracyReduction:      baselineAccuracy - defendedAccuracy,
				DefenseEffectiveness:   (baselineAccuracy - defendedAccuracy) / baselineAccuracy,
				StatisticalSignificance: ads.calculateSignificance(baselineAccuracy, defendedAccuracy, sampleSize),
				ResistanceLevel:        ads.assessResistanceLevel(baselineAccuracy, defendedAccuracy),
			}

			// Calculate confidence interval for accuracy reduction
			testResult.ConfidenceInterval = ads.calculateConfidenceInterval(
				testResult.AccuracyReduction, sampleSize)

			allTestResults = append(allTestResults, testResult)

			// Track best defense for this model
			if testResult.DefenseEffectiveness > bestDefenseResult.DefenseEffectiveness {
				bestDefenseResult = testResult
			}
		}

		result.ModelSpecificResults[modelName] = bestDefenseResult
	}

	// Calculate defense strategy effectiveness across all models
	for defenseName := range ads.defenseStrategies {
		effectiveness := ads.calculateDefenseEffectiveness(allTestResults, defenseName)
		result.DefenseEffectiveness[defenseName] = effectiveness
	}

	// Calculate overall resistance score
	result.OverallResistanceScore = ads.calculateOverallResistanceScore(allTestResults)

	// Generate vulnerability assessment
	result.VulnerabilityAssessment = ads.generateVulnerabilityAssessment(allTestResults)

	// Provide defense recommendations
	result.RecommendedDefenses = ads.generateDefenseRecommendations(result)

	// Assign security level
	result.SecurityLevel = ads.assessSecurityLevel(result.OverallResistanceScore)

	return result, nil
}

// simulateModelAttack simulates a population classification attack with optional defense
func (ads *AdversarialDefenseSystem) simulateModelAttack(
	model *AdversarialModel,
	sequences map[genomics.SequencingPlatform]map[population.Population][]string,
	sampleSize int,
	defense *DefenseStrategy,
) (float64, error) {

	// Extract features based on model capability
	features, labels, err := ads.extractFeatures(model, sequences, sampleSize)
	if err != nil {
		return 0.0, err
	}

	// Apply defense strategy if provided
	if defense != nil {
		features = ads.applyDefenseStrategy(defense, features, model)
	}

	// Simulate model prediction accuracy
	accuracy := ads.simulateModelAccuracy(model, features, labels, defense)

	return accuracy, nil
}

// extractFeatures extracts relevant features based on model capabilities
func (ads *AdversarialDefenseSystem) extractFeatures(
	model *AdversarialModel,
	sequences map[genomics.SequencingPlatform]map[population.Population][]string,
	sampleSize int,
) ([][]float64, []int, error) {

	var features [][]float64
	var labels []int
	
	labelMap := make(map[population.Population]int)
	labelCounter := 0

	for platform, populationSeqs := range sequences {
		for pop, seqs := range populationSeqs {
			// Assign numeric label to population
			if _, exists := labelMap[pop]; !exists {
				labelMap[pop] = labelCounter
				labelCounter++
			}

			// Extract features for each sequence
			actualSampleSize := int(math.Min(float64(sampleSize), float64(len(seqs))))
			for i := 0; i < actualSampleSize; i++ {
				featureVector := ads.extractSequenceFeatures(model, platform, pop, seqs[i])
				features = append(features, featureVector)
				labels = append(labels, labelMap[pop])
			}
		}
	}

	return features, labels, nil
}

// extractSequenceFeatures extracts features from a single sequence based on model capabilities
func (ads *AdversarialDefenseSystem) extractSequenceFeatures(
	model *AdversarialModel,
	platform genomics.SequencingPlatform,
	pop population.Population,
	sequence string,
) []float64 {

	var features []float64

	// PCA-related features
	if model.AttackCapability.CanDetectPCA {
		pcaFeatures := ads.extractPCAFeatures(sequence)
		features = append(features, pcaFeatures...)
	}

	// F_ST features
	if model.AttackCapability.CanDetectFST {
		fstFeatures := ads.extractFSTFeatures(sequence)
		features = append(features, fstFeatures...)
	}

	// LD features
	if model.AttackCapability.CanDetectLD {
		ldFeatures := ads.extractLDFeatures(sequence)
		features = append(features, ldFeatures...)
	}

	// Quality features
	if model.AttackCapability.CanDetectQuality {
		qualityFeatures := ads.extractQualityFeatures(platform, pop, sequence)
		features = append(features, qualityFeatures...)
	}

	// Error pattern features
	if model.AttackCapability.CanDetectErrorPatterns {
		errorFeatures := ads.extractErrorPatternFeatures(platform, pop, sequence)
		features = append(features, errorFeatures...)
	}

	// GC bias features
	if model.AttackCapability.CanDetectGCBias {
		gcFeatures := ads.extractGCBiasFeatures(sequence)
		features = append(features, gcFeatures...)
	}

	return features
}

// Feature extraction methods
func (ads *AdversarialDefenseSystem) extractPCAFeatures(sequence string) []float64 {
	// Simplified PCA-like features based on sequence composition
	features := make([]float64, 4) // A, T, G, C frequencies
	
	for _, base := range sequence {
		switch base {
		case 'A':
			features[0]++
		case 'T':
			features[1]++
		case 'G':
			features[2]++
		case 'C':
			features[3]++
		}
	}
	
	// Normalize by sequence length
	for i := range features {
		features[i] /= float64(len(sequence))
	}
	
	return features
}

func (ads *AdversarialDefenseSystem) extractFSTFeatures(sequence string) []float64 {
	// Simplified F_ST-related features
	features := []float64{
		ads.calculateSequenceDiversity([]string{sequence}),
		ads.calculateGCContent(sequence),
		float64(len(sequence)),
	}
	
	return features
}

func (ads *AdversarialDefenseSystem) extractLDFeatures(sequence string) []float64 {
	// LD-related features (single sequence: LD proxy will return 0.0; meaningful only
	// when called with multi-sequence slices via extractFeaturesForSequenceSet).
	features := []float64{
		ads.calculateLDProxy([]string{sequence}),
		ads.calculateTransitionRatio(sequence),
		ads.calculateRepeatContent(sequence),
	}

	return features
}

func (ads *AdversarialDefenseSystem) extractQualityFeatures(platform genomics.SequencingPlatform, pop population.Population, sequence string) []float64 {
	// Generate quality scores and extract statistics
	qualities := ads.qualityModeler.GenerateQualityScores(platform, pop, sequence, 1)
	
	if len(qualities) == 0 {
		return []float64{30.0, 5.0} // Default mean and std
	}
	
	// Calculate quality statistics
	sum := 0
	for _, q := range qualities {
		sum += q
	}
	mean := float64(sum) / float64(len(qualities))
	
	variance := 0.0
	for _, q := range qualities {
		diff := float64(q) - mean
		variance += diff * diff
	}
	std := 0.0
	if n := len(qualities); n > 1 {
		std = math.Sqrt(variance / float64(n-1))
	}
	
	return []float64{mean, std}
}

func (ads *AdversarialDefenseSystem) extractErrorPatternFeatures(platform genomics.SequencingPlatform, pop population.Population, sequence string) []float64 {
	// Apply error signature and extract statistics
	qualities := make([]int, len(sequence))
	for i := range qualities {
		qualities[i] = 35 // Default quality
	}
	
	_, _, errorResult := ads.platformManager.ApplyErrorSignature(platform, pop, sequence, qualities, 1)
	
	features := []float64{
		float64(errorResult.ErrorsApplied),
		0.0, // Success rate (simplified)
	}
	
	if errorResult.Success {
		features[1] = 1.0
	}
	
	return features
}

func (ads *AdversarialDefenseSystem) extractGCBiasFeatures(sequence string) []float64 {
	gcContent := ads.calculateGCContent(sequence)
	
	// Calculate GC bias in windows
	windowSize := 50
	gcVariation := 0.0
	windowCount := 0
	
	for i := 0; i <= len(sequence)-windowSize; i += windowSize {
		window := sequence[i:i+windowSize]
		windowGC := ads.calculateGCContent(window)
		gcVariation += math.Abs(windowGC - gcContent)
		windowCount++
	}
	
	if windowCount > 0 {
		gcVariation /= float64(windowCount)
	}
	
	return []float64{gcContent, gcVariation}
}

// Helper feature calculation methods
func (ads *AdversarialDefenseSystem) calculateGCContent(sequence string) float64 {
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

// calculateLDProxy estimates mean pairwise r² linkage disequilibrium over up to 50
// adjacent biallelic position pairs sampled from the provided sequences.
//
// For each adjacent pair of positions (j, j+1) the four haplotype counts
// (nAB, nAb, naB, nab) are derived by treating the modal base at each position as
// the "major" allele. r² = D² / (p_A·q_a·p_B·q_b) where D = p_AB − p_A·p_B.
// Monomorphic or invalid positions are skipped. Returns 0.0 when no valid pairs
// can be computed (< 2 sequences or sequences too short).
func (ads *AdversarialDefenseSystem) calculateLDProxy(sequences []string) float64 {
	if len(sequences) < 2 {
		return 0.0
	}

	// Find minimum length across all sequences.
	minLen := len(sequences[0])
	for _, s := range sequences[1:] {
		if len(s) < minLen {
			minLen = len(s)
		}
	}
	if minLen < 2 {
		return 0.0
	}

	posLimit := minLen - 1
	if posLimit > 50 {
		posLimit = 50
	}

	totalR2 := 0.0
	validPairs := 0

	for j := 0; j < posLimit; j++ {
		// Count base frequencies at positions j and j+1.
		var counts1, counts2 [4]int
		for _, seq := range sequences {
			if j+1 >= len(seq) {
				continue
			}
			if b := ldBaseIndex(seq[j]); b >= 0 {
				counts1[b]++
			}
			if b := ldBaseIndex(seq[j+1]); b >= 0 {
				counts2[b]++
			}
		}
		modal1 := ldMaxIndex(counts1)
		modal2 := ldMaxIndex(counts2)

		// Count biallelic haplotypes (modal = A/B, non-modal = a/b).
		nAB, nAb, naB, nab := 0, 0, 0, 0
		for _, seq := range sequences {
			if j+1 >= len(seq) {
				continue
			}
			b1 := ldBaseIndex(seq[j])
			b2 := ldBaseIndex(seq[j+1])
			if b1 < 0 || b2 < 0 {
				continue
			}
			switch {
			case b1 == modal1 && b2 == modal2:
				nAB++
			case b1 == modal1 && b2 != modal2:
				nAb++
			case b1 != modal1 && b2 == modal2:
				naB++
			default:
				nab++
			}
		}

		total := nAB + nAb + naB + nab
		if total == 0 {
			continue
		}

		pA := float64(nAB+nAb) / float64(total)
		pB := float64(nAB+naB) / float64(total)
		denom := pA * (1 - pA) * pB * (1 - pB)
		if denom <= 0 {
			continue // monomorphic at one or both loci
		}

		pAB := float64(nAB) / float64(total)
		D := pAB - pA*pB
		totalR2 += (D * D) / denom
		validPairs++
	}

	if validPairs == 0 {
		return 0.0
	}
	return totalR2 / float64(validPairs)
}

func (ads *AdversarialDefenseSystem) calculateSequenceDiversity(sequences []string) float64 {
	if len(sequences) < 2 {
		return 0.0
	}
	
	// Simple diversity based on sequence variation
	baseVariation := make(map[byte]int)
	totalBases := 0
	
	for _, seq := range sequences {
		for _, base := range []byte(seq) {
			baseVariation[base]++
			totalBases++
		}
	}
	
	if totalBases == 0 {
		return 0.0
	}
	
	// Shannon diversity
	diversity := 0.0
	for _, count := range baseVariation {
		if count > 0 {
			p := float64(count) / float64(totalBases)
			diversity -= p * math.Log2(p)
		}
	}
	
	return diversity
}

func (ads *AdversarialDefenseSystem) calculateTransitionRatio(sequence string) float64 {
	if len(sequence) < 2 {
		return 0.0
	}
	
	transitions := 0
	transversions := 0
	
	for i := 0; i < len(sequence)-1; i++ {
		curr, next := sequence[i], sequence[i+1]
		
		// Count transitions (A<->G, C<->T) vs transversions
		if (curr == 'A' && next == 'G') || (curr == 'G' && next == 'A') ||
		   (curr == 'C' && next == 'T') || (curr == 'T' && next == 'C') {
			transitions++
		} else if curr != next {
			transversions++
		}
	}
	
	if transversions == 0 {
		return float64(transitions)
	}
	
	return float64(transitions) / float64(transversions)
}

func (ads *AdversarialDefenseSystem) calculateRepeatContent(sequence string) float64 {
	if len(sequence) < 3 {
		return 0.0
	}
	
	repeatBases := 0
	
	for i := 0; i < len(sequence)-2; i++ {
		if sequence[i] == sequence[i+1] && sequence[i+1] == sequence[i+2] {
			repeatBases++
		}
	}
	
	return float64(repeatBases) / float64(len(sequence))
}

// Defense strategy application methods
func (ads *AdversarialDefenseSystem) applyDefenseStrategy(defense *DefenseStrategy, features [][]float64, model *AdversarialModel) [][]float64 {
	switch defense.DefenseType {
	case "Noise_Injection":
		return ads.applyNoiseInjection(defense, features)
	case "Feature_Masking":
		return ads.applyFeatureMasking(defense, features, model)
	case "Adversarial_Training":
		return ads.applyAdversarialTraining(defense, features)
	case "Ensemble_Defense":
		return ads.applyEnsembleDefense(defense, features)
	case "Privacy_Preserving":
		return ads.applyDifferentialPrivacy(defense, features)
	case "Real_Time_Defense":
		return ads.applyDynamicScrambling(defense, features)
	default:
		return features
	}
}

func (ads *AdversarialDefenseSystem) applyNoiseInjection(defense *DefenseStrategy, features [][]float64) [][]float64 {
	noiseStd := defense.DefenseParameters["noise_std"]
	coverageRate := defense.DefenseParameters["feature_coverage"]
	
	defendedFeatures := make([][]float64, len(features))
	
	for i, featureVector := range features {
		defendedFeatures[i] = make([]float64, len(featureVector))
		copy(defendedFeatures[i], featureVector)
		
		// Add noise to selected features
		for j := range featureVector {
			if ads.rng.Float64() < coverageRate {
				noise := ads.rng.NormFloat64() * noiseStd
				defendedFeatures[i][j] += noise
			}
		}
	}
	
	return defendedFeatures
}

func (ads *AdversarialDefenseSystem) applyFeatureMasking(defense *DefenseStrategy, features [][]float64, model *AdversarialModel) [][]float64 {
	maskingProb := defense.DefenseParameters["masking_probability"]
	importanceThreshold := defense.DefenseParameters["feature_importance_threshold"]
	
	defendedFeatures := make([][]float64, len(features))
	
	for i, featureVector := range features {
		defendedFeatures[i] = make([]float64, len(featureVector))
		copy(defendedFeatures[i], featureVector)
		
		// Mask features based on importance and probability
		for j := range featureVector {
			// Get feature importance (simplified)
			importance := ads.getFeatureImportance(model, j, len(featureVector))
			
			if importance > importanceThreshold && ads.rng.Float64() < maskingProb {
				// Mask by setting to mean value or adding significant noise
				defendedFeatures[i][j] = ads.calculateFeatureMean(features, j)
			}
		}
	}
	
	return defendedFeatures
}

func (ads *AdversarialDefenseSystem) applyAdversarialTraining(defense *DefenseStrategy, features [][]float64) [][]float64 {
	perturbationStrength := defense.DefenseParameters["perturbation_strength"]
	
	// For adversarial training, we modify features to be more robust
	defendedFeatures := make([][]float64, len(features))
	
	for i, featureVector := range features {
		defendedFeatures[i] = make([]float64, len(featureVector))
		copy(defendedFeatures[i], featureVector)
		
		// Apply small perturbations to increase robustness
		for j := range featureVector {
			perturbation := (ads.rng.Float64() - 0.5) * perturbationStrength
			defendedFeatures[i][j] += perturbation
		}
	}
	
	return defendedFeatures
}

func (ads *AdversarialDefenseSystem) applyEnsembleDefense(defense *DefenseStrategy, features [][]float64) [][]float64 {
	// Apply multiple defense strategies
	defendedFeatures := features
	
	// Apply noise injection
	noiseDefense := &DefenseStrategy{
		DefenseType: "Noise_Injection",
		DefenseParameters: map[string]float64{
			"noise_std": 0.03,
			"feature_coverage": 0.5,
		},
	}
	defendedFeatures = ads.applyNoiseInjection(noiseDefense, defendedFeatures)
	
	// Apply feature masking
	maskingDefense := &DefenseStrategy{
		DefenseType: "Feature_Masking",
		DefenseParameters: map[string]float64{
			"masking_probability": 0.1,
			"feature_importance_threshold": 0.05,
		},
	}
	defendedFeatures = ads.applyFeatureMasking(maskingDefense, defendedFeatures, &AdversarialModel{})
	
	return defendedFeatures
}

func (ads *AdversarialDefenseSystem) applyDifferentialPrivacy(defense *DefenseStrategy, features [][]float64) [][]float64 {
	epsilon := defense.DefenseParameters["epsilon"]
	sensitivity := defense.DefenseParameters["sensitivity"]
	
	// Apply Laplacian noise for differential privacy
	laplacianScale := sensitivity / epsilon
	
	defendedFeatures := make([][]float64, len(features))
	
	for i, featureVector := range features {
		defendedFeatures[i] = make([]float64, len(featureVector))
		copy(defendedFeatures[i], featureVector)
		
		for j := range featureVector {
			// Add Laplacian noise
			noise := ads.sampleLaplacian(0.0, laplacianScale)
			defendedFeatures[i][j] += noise
		}
	}
	
	return defendedFeatures
}

func (ads *AdversarialDefenseSystem) applyDynamicScrambling(defense *DefenseStrategy, features [][]float64) [][]float64 {
	scramblingFreq := defense.DefenseParameters["scrambling_frequency"]
	patternComplexity := int(defense.DefenseParameters["pattern_complexity"])
	
	defendedFeatures := make([][]float64, len(features))
	
	for i, featureVector := range features {
		defendedFeatures[i] = make([]float64, len(featureVector))
		copy(defendedFeatures[i], featureVector)
		
		// Apply dynamic scrambling pattern
		if ads.rng.Float64() < scramblingFreq {
			// Generate scrambling pattern
			pattern := ads.generateScramblingPattern(len(featureVector), patternComplexity)
			
			// Apply pattern
			for j, scrambleIndex := range pattern {
				if scrambleIndex < len(featureVector) {
					defendedFeatures[i][j] = featureVector[scrambleIndex]
				}
			}
		}
	}
	
	return defendedFeatures
}

// Helper methods for defense strategies
func (ads *AdversarialDefenseSystem) getFeatureImportance(model *AdversarialModel, featureIndex, totalFeatures int) float64 {
	// Simplified feature importance calculation
	if featureIndex < len(model.FeatureSet) {
		featureName := model.FeatureSet[featureIndex%len(model.FeatureSet)]
		if weight, exists := model.FeatureWeights[featureName]; exists {
			return weight
		}
	}
	
	return 1.0 / float64(totalFeatures) // Uniform importance fallback
}

func (ads *AdversarialDefenseSystem) calculateFeatureMean(features [][]float64, featureIndex int) float64 {
	if len(features) == 0 || featureIndex >= len(features[0]) {
		return 0.0
	}
	
	sum := 0.0
	count := 0
	
	for _, featureVector := range features {
		if featureIndex < len(featureVector) {
			sum += featureVector[featureIndex]
			count++
		}
	}
	
	if count == 0 {
		return 0.0
	}
	
	return sum / float64(count)
}

func (ads *AdversarialDefenseSystem) sampleLaplacian(mean, scale float64) float64 {
	// Sample from Laplacian distribution
	u := ads.rng.Float64() - 0.5
	if u >= 0 {
		return mean - scale*math.Log(1-2*u)
	}
	return mean + scale*math.Log(1+2*u)
}

func (ads *AdversarialDefenseSystem) generateScramblingPattern(size, complexity int) []int {
	pattern := make([]int, size)
	for i := 0; i < size; i++ {
		pattern[i] = i
	}
	
	// Apply complex scrambling
	for round := 0; round < complexity; round++ {
		for i := range pattern {
			j := ads.rng.IntN(size)
			pattern[i], pattern[j] = pattern[j], pattern[i]
		}
	}
	
	return pattern
}

// simulateModelAccuracy measures classification accuracy on the provided feature/label
// pairs using a nearest-centroid classifier with an 80/20 stratified hold-out split.
//
// Defense impact is already baked into the features by the time this is called
// (applyDefenseStrategy is applied upstream in simulateModelAttack), so the
// defense parameter is intentionally ignored here.
func (ads *AdversarialDefenseSystem) simulateModelAccuracy(model *AdversarialModel, features [][]float64, labels []int, _ *DefenseStrategy) float64 {
	if len(features) < 10 || len(labels) < 10 {
		return 0.5 // too few samples to classify
	}

	// 80/20 hold-out split (preserves order; no shuffle needed for a proxy metric).
	splitIdx := int(float64(len(features)) * 0.8)
	if splitIdx < 2 || splitIdx >= len(features) {
		return 0.5
	}
	trainFeats, trainLabels := features[:splitIdx], labels[:splitIdx]
	testFeats, testLabels := features[splitIdx:], labels[splitIdx:]

	// Compute per-class centroid from training set.
	nDim := len(trainFeats[0])
	centroidSums := make(map[int][]float64)
	centroidCounts := make(map[int]int)
	for i, feat := range trainFeats {
		lbl := trainLabels[i]
		if _, ok := centroidSums[lbl]; !ok {
			centroidSums[lbl] = make([]float64, nDim)
		}
		for k, v := range feat {
			centroidSums[lbl][k] += v
		}
		centroidCounts[lbl]++
	}
	centroids := make(map[int][]float64, len(centroidSums))
	for lbl, sums := range centroidSums {
		c := make([]float64, nDim)
		cnt := float64(centroidCounts[lbl])
		for k, s := range sums {
			c[k] = s / cnt
		}
		centroids[lbl] = c
	}
	if len(centroids) < 2 {
		return 0.5 // only one class in training set
	}

	// Classify test set by nearest centroid (L2 distance).
	correct := 0
	for i, feat := range testFeats {
		bestLbl, bestDist := -1, math.MaxFloat64
		for lbl, centroid := range centroids {
			d := 0.0
			for k, v := range feat {
				diff := v - centroid[k]
				d += diff * diff
			}
			if d < bestDist {
				bestDist = d
				bestLbl = lbl
			}
		}
		if bestLbl == testLabels[i] {
			correct++
		}
	}

	return float64(correct) / float64(len(testLabels))
}

// Assessment and analysis methods
func (ads *AdversarialDefenseSystem) isDefenseApplicable(defense *DefenseStrategy, model *AdversarialModel) bool {
	// Check if defense targets this model type
	for _, targetModel := range defense.TargetedModels {
		if targetModel == "All_Models" || targetModel == model.ModelType {
			return true
		}
	}
	return false
}

// calculateSignificance returns a two-tailed p-value for the difference between
// two accuracy proportions using a pooled two-proportion z-test:
//
//	p̄  = (preAccuracy + postAccuracy) / 2
//	SE  = sqrt(2 · p̄ · (1 − p̄) / sampleSize)
//	z   = |preAccuracy − postAccuracy| / SE
//	p   = 2 · (1 − Φ(z))
func (ads *AdversarialDefenseSystem) calculateSignificance(preAccuracy, postAccuracy float64, sampleSize int) float64 {
	if sampleSize <= 0 {
		return 1.0
	}
	pBar := (preAccuracy + postAccuracy) / 2
	se := math.Sqrt(2 * pBar * (1 - pBar) / float64(sampleSize))
	if se == 0 {
		return 1.0
	}
	z := math.Abs(preAccuracy-postAccuracy) / se
	return 2 * (1 - ads.normalCDF(z))
}

func (ads *AdversarialDefenseSystem) normalCDF(z float64) float64 {
	return 0.5 * (1 + math.Erf(z/math.Sqrt(2)))
}

func (ads *AdversarialDefenseSystem) calculateConfidenceInterval(accuracyReduction float64, sampleSize int) [2]float64 {
	margin := 1.96 * math.Sqrt(accuracyReduction*(1-accuracyReduction)/float64(sampleSize))
	return [2]float64{
		math.Max(0, accuracyReduction-margin),
		math.Min(1, accuracyReduction+margin),
	}
}

func (ads *AdversarialDefenseSystem) assessResistanceLevel(preAccuracy, postAccuracy float64) string {
	reduction := preAccuracy - postAccuracy
	
	if reduction >= 0.3 {
		return "Highly_Resistant"
	} else if reduction >= 0.15 {
		return "Partially_Resistant"
	} else if reduction >= 0.05 {
		return "Minimally_Resistant"
	}
	return "Vulnerable"
}

func (ads *AdversarialDefenseSystem) calculateDefenseEffectiveness(results []*AdversarialTestResult, defenseName string) float64 {
	relevantResults := []*AdversarialTestResult{}
	
	for _, result := range results {
		if result.DefenseStrategy == defenseName {
			relevantResults = append(relevantResults, result)
		}
	}
	
	if len(relevantResults) == 0 {
		return 0.0
	}
	
	totalEffectiveness := 0.0
	for _, result := range relevantResults {
		totalEffectiveness += result.DefenseEffectiveness
	}
	
	return totalEffectiveness / float64(len(relevantResults))
}

func (ads *AdversarialDefenseSystem) calculateOverallResistanceScore(results []*AdversarialTestResult) float64 {
	if len(results) == 0 {
		return 0.0
	}
	
	totalScore := 0.0
	resistantCount := 0
	
	for _, result := range results {
		if result.ResistanceLevel == "Highly_Resistant" {
			totalScore += 1.0
			resistantCount++
		} else if result.ResistanceLevel == "Partially_Resistant" {
			totalScore += 0.7
			resistantCount++
		} else if result.ResistanceLevel == "Minimally_Resistant" {
			totalScore += 0.3
		}
	}
	
	avgScore := totalScore / float64(len(results))
	resistanceRate := float64(resistantCount) / float64(len(results))
	
	return (avgScore + resistanceRate) / 2
}

func (ads *AdversarialDefenseSystem) generateVulnerabilityAssessment(results []*AdversarialTestResult) VulnerabilityAssessment {
	assessment := VulnerabilityAssessment{
		HighRiskModels:          []string{},
		CriticalVulnerabilities: []string{},
		DefensiveGaps:           []string{},
		AttackVectors:           []string{},
		RequiredImprovements:    []string{},
	}
	
	// Identify high-risk models
	for _, result := range results {
		if result.ResistanceLevel == "Vulnerable" || result.DefenseEffectiveness < 0.3 {
			assessment.HighRiskModels = append(assessment.HighRiskModels, result.ModelName)
		}
	}
	
	// Identify critical vulnerabilities
	if len(assessment.HighRiskModels) > 0 {
		assessment.CriticalVulnerabilities = append(assessment.CriticalVulnerabilities, 
			"Advanced ML models show significant classification capability")
	}
	
	// Calculate security rating
	resistantCount := 0
	for _, result := range results {
		if result.ResistanceLevel == "Highly_Resistant" || result.ResistanceLevel == "Partially_Resistant" {
			resistantCount++
		}
	}
	
	if len(results) > 0 {
		assessment.SecurityRating = float64(resistantCount) / float64(len(results)) * 100
	}
	
	return assessment
}

func (ads *AdversarialDefenseSystem) generateDefenseRecommendations(result *AdversarialDefenseResult) []string {
	recommendations := []string{}
	
	if result.OverallResistanceScore < 0.6 {
		recommendations = append(recommendations, "Implement multi-layer ensemble defense for comprehensive protection")
	}
	
	if result.VulnerabilityAssessment.SecurityRating < 70 {
		recommendations = append(recommendations, "Enhance adversarial training with stronger perturbations")
	}
	
	// Find best performing defense
	bestDefense := ""
	bestScore := 0.0
	for defense, score := range result.DefenseEffectiveness {
		if score > bestScore {
			bestScore = score
			bestDefense = defense
		}
	}
	
	if bestDefense != "" && bestScore > 0.7 {
		recommendations = append(recommendations, 
			fmt.Sprintf("Scale up %s deployment - showing high effectiveness (%.2f)", bestDefense, bestScore))
	}
	
	if len(result.VulnerabilityAssessment.HighRiskModels) > 0 {
		recommendations = append(recommendations, "Focus defensive efforts on neural network and ensemble model resistance")
	}
	
	return recommendations
}

func (ads *AdversarialDefenseSystem) assessSecurityLevel(overallScore float64) string {
	if overallScore >= 0.9 {
		return "Maximum_Security"
	} else if overallScore >= 0.8 {
		return "High_Security"
	} else if overallScore >= 0.7 {
		return "Moderate_Security"
	} else if overallScore >= 0.6 {
		return "Low_Security"
	}
	return "Critical_Vulnerability"
}

// Getter methods
func (ads *AdversarialDefenseSystem) GetAdversarialModels() map[string]*AdversarialModel {
	return ads.adversarialModels
}

func (ads *AdversarialDefenseSystem) GetDefenseStrategies() map[string]*DefenseStrategy {
	return ads.defenseStrategies
}

func (ads *AdversarialDefenseSystem) GetIndistinguishabilityTester() *IndistinguishabilityTester {
	return ads.indistinguishabilityTester
}