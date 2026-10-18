// Package testing is the end-to-end test harness used by cmd/run_e2e_tests
// and tests/e2e: per-platform synthetic read generation, FASTQ, SAM, BAM, VCF
// and reference writers, and an ML calibration (Multinomial Naive Bayes on
// 4-mer counts, see MultinomialNB).
//
// Several reported scores are heuristics or constants rather than
// measurements; their doc comments say so. In particular AnalyzeAuthenticity
// does not measure authenticity, RunWorkflow does not run a workflow (the
// real engine tests are in tests/integration), and RunBenchmark's accuracy,
// cost and long-read throughput figures are fixed per platform. No real
// sequencing data is bundled, so the ML calibration compares the generator
// with itself (see GenerateBalancedDataset).
package testing

import (
	"errors"
	"fmt"
	"math"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// SequencingPlatform represents different sequencing technologies
type SequencingPlatform string

const (
	// Short-read platforms
	IlluminaNovaSeq6000 SequencingPlatform = "illumina_novaseq_6000"
	IlluminaNextSeq2000 SequencingPlatform = "illumina_nextseq_2000"
	IlluminaMiSeq       SequencingPlatform = "illumina_miseq"

	// Long-read platforms
	PacBioRevio     SequencingPlatform = "pacbio_revio"
	ONTPromethION24 SequencingPlatform = "ont_promethion_24"

	// Global alternative platforms
	MGI_DNBSEQ_T20x2 SequencingPlatform = "mgi_dnbseq_t20x2"
	ElementAVITI     SequencingPlatform = "element_aviti"
)

// WorkflowEngine represents different workflow management systems
type WorkflowEngine string

const (
	CromwellWDL   WorkflowEngine = "cromwell_wdl"
	NextflowNFL   WorkflowEngine = "nextflow_nfl"
	SnakemakeRule WorkflowEngine = "snakemake_rule"
)

// E2EFramework holds per-platform configurations and generates and writes
// synthetic test data for them.
type E2EFramework struct {
	authenticityEngine  *authenticity.Engine
	qualityModeler      interface{}
	populationManager   interface{}
	referenceGenome     interface{}
	validationFramework interface{}

	// Platform-specific configurations
	platformConfigs map[SequencingPlatform]*PlatformConfig

	// Performance tracking
	benchmarkResults map[SequencingPlatform]*PlatformBenchmark

	// Thread safety
	mu sync.RWMutex
}

// PlatformConfig holds platform-specific configuration
type PlatformConfig struct {
	Platform         SequencingPlatform
	Population       string
	ReadLength       int
	ErrorRate        float64
	QualityProfile   *QualityProfile
	BiasCorrection   *BiasCorrection
	TechnologyParams map[string]interface{}
}

// QualityProfile defines platform-specific quality characteristics
type QualityProfile struct {
	AverageQuality  float64
	QualityRange    [2]float64
	GCBias          float64
	HomopolymerBias float64
	ErrorModel      string
}

// BiasCorrection defines platform-specific bias correction parameters
type BiasCorrection struct {
	GCCorrection      bool
	LengthCorrection  bool
	QualityCorrection bool
	Method            string // "BEADS", "MANCIE", "Quantile"
}

// SyntheticData represents generated synthetic genomic data
type SyntheticData struct {
	Sequences     [][]byte
	QualityScores [][]byte
	Headers       []string
	Platform      SequencingPlatform
	Population    string

	// Platform-specific characteristics
	AverageReadLength   float64
	ErrorRate           float64
	QualityScoreAverage float64

	// Population-specific features
	NucleotideDiversity     float64
	VariantCount            int
	LinkageDisequilibrium   float64
	AdmixtureLD             float64
	EffectivePopulationSize int
	FounderEffect           float64

	// Technology-specific features
	TwoColorBias        *TwoColorBias
	CircularConsensus   *CircularConsensus
	HomopolymerErrors   *HomopolymerErrors
	DNBBias             *DNBBias
	SBBChemistry        *SBBChemistry
	PopulationStructure *PopulationStructure
	AdmixtureComponents map[string]float64
}

// Technology-specific bias structures
type TwoColorBias struct {
	GTBias float64
	ACBias float64
}

type CircularConsensus struct {
	PassCount int
	Accuracy  float64
}

type HomopolymerErrors struct {
	InsertionRate float64
	DeletionRate  float64
}

type DNBBias struct {
	GCBiasVariation float64
}

type SBBChemistry struct {
	BindingEfficiency float64
}

type PopulationStructure struct {
	EastAsianComponent float64
	EuropeanComponent  float64
	AfricanComponent   float64
}

// PlatformBenchmark holds performance metrics for a platform
type PlatformBenchmark struct {
	Platform            SequencingPlatform
	ThroughputMBps      float64
	LatencyMs           int
	AccuracyScore       float64
	AuthenticityScore   float64
	CPUUsagePercent     float64
	MemoryUsageMB       int64
	LongReadThroughput  float64
	CostPerGB           float64
	SystematicErrorRate float64
	OverallScore        float64
}

// MLClassifier provides machine learning-based authenticity validation using
// a Multinomial Naive Bayes classifier on 4-mer frequency vectors.
type MLClassifier struct {
	nb      *MultinomialNB
	trained bool
}

// WorkflowResult holds workflow execution results
type WorkflowResult struct {
	SuccessRate    float64
	ExecutionTime  time.Duration
	GATKCompatible bool
	BWACompatible  bool
	ErrorMessages  []string
}

// ConcurrentResult holds results from concurrent platform testing
type ConcurrentResult struct {
	Platform          SequencingPlatform
	Success           bool
	Duration          time.Duration
	DataSize          int
	AuthenticityScore float64
}

// NewE2EFramework creates a new end-to-end testing framework
func NewE2EFramework() *E2EFramework {
	framework := &E2EFramework{
		platformConfigs:  make(map[SequencingPlatform]*PlatformConfig),
		benchmarkResults: make(map[SequencingPlatform]*PlatformBenchmark),
	}

	// Initialize platform configurations
	framework.initializePlatformConfigs()

	// Initialize core components
	framework.initializeCoreComponents()

	return framework
}

// initializePlatformConfigs sets up platform-specific configurations
func (f *E2EFramework) initializePlatformConfigs() {
	f.platformConfigs[IlluminaNovaSeq6000] = &PlatformConfig{
		Platform:   IlluminaNovaSeq6000,
		ReadLength: 150,
		ErrorRate:  0.0015,
		QualityProfile: &QualityProfile{
			AverageQuality:  35.0,
			QualityRange:    [2]float64{20.0, 42.0},
			GCBias:          0.02,
			HomopolymerBias: 0.01,
			ErrorModel:      "illumina_phred33",
		},
		BiasCorrection: &BiasCorrection{
			GCCorrection:      true,
			LengthCorrection:  true,
			QualityCorrection: true,
			Method:            "BEADS",
		},
		TechnologyParams: map[string]interface{}{
			"chemistry":       "SBS",
			"cluster_density": 200000,
			"phasing_rate":    0.003,
			"prephasing_rate": 0.002,
		},
	}

	f.platformConfigs[IlluminaNextSeq2000] = &PlatformConfig{
		Platform:   IlluminaNextSeq2000,
		ReadLength: 150,
		ErrorRate:  0.002,
		QualityProfile: &QualityProfile{
			AverageQuality:  33.0,
			QualityRange:    [2]float64{18.0, 40.0},
			GCBias:          0.03,
			HomopolymerBias: 0.015,
			ErrorModel:      "illumina_two_color",
		},
		BiasCorrection: &BiasCorrection{
			GCCorrection:      true,
			LengthCorrection:  true,
			QualityCorrection: true,
			Method:            "MANCIE",
		},
		TechnologyParams: map[string]interface{}{
			"chemistry":     "two_color_SBS",
			"green_channel": "G_T",
			"red_channel":   "A_C",
		},
	}

	f.platformConfigs[PacBioRevio] = &PlatformConfig{
		Platform:   PacBioRevio,
		ReadLength: 15000,
		ErrorRate:  0.003,
		QualityProfile: &QualityProfile{
			AverageQuality:  40.0,
			QualityRange:    [2]float64{30.0, 50.0},
			GCBias:          0.05,
			HomopolymerBias: 0.02,
			ErrorModel:      "pacbio_hifi",
		},
		BiasCorrection: &BiasCorrection{
			GCCorrection:      true,
			LengthCorrection:  false,
			QualityCorrection: true,
			Method:            "Quantile",
		},
		TechnologyParams: map[string]interface{}{
			"chemistry":   "HiFi",
			"ccs_passes":  3,
			"zmw_density": 25000000,
			"polymerase":  "P6-C4",
		},
	}

	f.platformConfigs[ONTPromethION24] = &PlatformConfig{
		Platform:   ONTPromethION24,
		ReadLength: 30000,
		ErrorRate:  0.08,
		QualityProfile: &QualityProfile{
			AverageQuality:  15.0,
			QualityRange:    [2]float64{5.0, 25.0},
			GCBias:          0.1,
			HomopolymerBias: 0.15,
			ErrorModel:      "ont_r10_4_1",
		},
		BiasCorrection: &BiasCorrection{
			GCCorrection:      true,
			LengthCorrection:  false,
			QualityCorrection: true,
			Method:            "BEADS",
		},
		TechnologyParams: map[string]interface{}{
			"chemistry":  "R10.4.1",
			"pore_type":  "protein_nanopore",
			"basecaller": "Guppy_6.1.5",
			"flow_cell":  "FLO-PRO114M",
		},
	}

	// Add other platform configurations...
	f.addRemainingPlatformConfigs()
}

// addRemainingPlatformConfigs adds configurations for remaining platforms
func (f *E2EFramework) addRemainingPlatformConfigs() {
	f.platformConfigs[MGI_DNBSEQ_T20x2] = &PlatformConfig{
		Platform:   MGI_DNBSEQ_T20x2,
		ReadLength: 150,
		ErrorRate:  0.0018,
		QualityProfile: &QualityProfile{
			AverageQuality:  34.0,
			QualityRange:    [2]float64{19.0, 41.0},
			GCBias:          0.015,
			HomopolymerBias: 0.008,
			ErrorModel:      "mgi_dnb",
		},
		TechnologyParams: map[string]interface{}{
			"chemistry":     "DNB",
			"nanoball_size": 300,
			"array_density": 4000000,
		},
	}

	f.platformConfigs[ElementAVITI] = &PlatformConfig{
		Platform:   ElementAVITI,
		ReadLength: 150,
		ErrorRate:  0.001,
		QualityProfile: &QualityProfile{
			AverageQuality:  38.0,
			QualityRange:    [2]float64{25.0, 45.0},
			GCBias:          0.008,
			HomopolymerBias: 0.005,
			ErrorModel:      "element_sbb",
		},
		TechnologyParams: map[string]interface{}{
			"chemistry":          "SBB",
			"binding_efficiency": 0.99,
			"avidity_sequencing": true,
		},
	}

	f.platformConfigs[IlluminaMiSeq] = &PlatformConfig{
		Platform:   IlluminaMiSeq,
		ReadLength: 300,
		ErrorRate:  0.001,
		QualityProfile: &QualityProfile{
			AverageQuality:  36.0,
			QualityRange:    [2]float64{20.0, 42.0},
			GCBias:          0.02,
			HomopolymerBias: 0.01,
			ErrorModel:      "illumina_phred33",
		},
		BiasCorrection: &BiasCorrection{
			GCCorrection:      true,
			LengthCorrection:  true,
			QualityCorrection: true,
			Method:            "BEADS",
		},
		TechnologyParams: map[string]interface{}{
			"chemistry":       "SBS",
			"cluster_density": 14000,
			"phasing_rate":    0.002,
			"prephasing_rate": 0.002,
		},
	}
}

// initializeCoreComponents initializes the core genomics components
func (f *E2EFramework) initializeCoreComponents() {
	// Initialize authenticity engine
	config := &authenticity.Config{
		KmerSizes:        []int{4, 5, 6},
		JSThreshold:      0.01,
		EnableValidation: true,
	}
	f.authenticityEngine = authenticity.NewEngine(config)

	// Quality modeler, population manager, and validation framework
	// are placeholders pending implementation
	f.qualityModeler = nil
	f.populationManager = nil
	f.validationFramework = nil
}

// CreatePlatformConfig creates a configuration for a specific platform and population
func (f *E2EFramework) CreatePlatformConfig(platform SequencingPlatform, population string) *PlatformConfig {
	f.mu.RLock()
	baseConfig, exists := f.platformConfigs[platform]
	f.mu.RUnlock()

	if !exists {
		return nil
	}

	// Create a copy with population-specific modifications
	config := *baseConfig
	config.Population = population

	// Apply population-specific adjustments
	f.applyPopulationAdjustments(&config, population)

	return &config
}

// applyPopulationAdjustments applies population-specific adjustments to platform config
func (f *E2EFramework) applyPopulationAdjustments(config *PlatformConfig, population string) {
	switch population {
	case "YRI":
		// African populations have higher genetic diversity
		config.TechnologyParams["diversity_multiplier"] = 1.3
		config.TechnologyParams["variant_density"] = 1.5

	case "CEU":
		// European populations show bottleneck signatures
		config.TechnologyParams["diversity_multiplier"] = 0.8
		config.TechnologyParams["ld_decay_rate"] = 0.7

	case "CHB":
		// East Asian populations have specific structure
		config.TechnologyParams["diversity_multiplier"] = 0.75
		config.TechnologyParams["population_structure"] = "east_asian"

	case "JPT":
		// Japanese show island founder effects
		config.TechnologyParams["diversity_multiplier"] = 0.65
		config.TechnologyParams["founder_effect"] = 0.15

	case "MXL":
		// Mexican populations are admixed
		config.TechnologyParams["admixture_components"] = map[string]float64{
			"native_american": 0.45,
			"european":        0.35,
			"african":         0.20,
		}
	}
}

// GenerateSyntheticData generates synthetic genomic data for a platform configuration
func (f *E2EFramework) GenerateSyntheticData(config *PlatformConfig) (*SyntheticData, error) {
	return f.generateSyntheticDataN(config, 1000)
}

// generateSyntheticDataN generates numSequences reads for the given config.
func (f *E2EFramework) generateSyntheticDataN(config *PlatformConfig, numSequences int) (*SyntheticData, error) {
	if config == nil {
		return nil, fmt.Errorf("configuration cannot be nil")
	}
	if numSequences <= 0 {
		return nil, fmt.Errorf("numSequences must be positive, got %d", numSequences)
	}

	sequenceLength := config.ReadLength

	data := &SyntheticData{
		Sequences:         make([][]byte, numSequences),
		QualityScores:     make([][]byte, numSequences),
		Headers:           make([]string, numSequences),
		Platform:          config.Platform,
		Population:        config.Population,
		AverageReadLength: float64(config.ReadLength),
		ErrorRate:         config.ErrorRate,
	}

	// Generate sequences with platform-specific characteristics
	for i := 0; i < numSequences; i++ {
		sequence, quality := f.generateSequenceWithQuality(config, sequenceLength)
		data.Sequences[i] = sequence
		data.QualityScores[i] = quality
		data.Headers[i] = fmt.Sprintf("@seq_%d_%s_%s", i, config.Platform, config.Population)
	}

	// Apply platform-specific modeling
	f.applyPlatformSpecificModeling(data, config)

	// Apply population-specific features
	f.applyPopulationSpecificFeatures(data, config)

	// Calculate derived statistics
	f.calculateDerivedStatistics(data)

	return data, nil
}

// generateSequenceWithQuality generates a single sequence with quality scores
func (f *E2EFramework) generateSequenceWithQuality(config *PlatformConfig, length int) ([]byte, []byte) {
	sequence := make([]byte, length)
	quality := make([]byte, length)

	bases := []byte{'A', 'T', 'G', 'C'}

	for i := 0; i < length; i++ {
		// Generate base with potential errors
		baseIndex := rand.Intn(4) // #nosec G404
		errProb := rand.Float64() // #nosec G404
		if errProb < config.ErrorRate {
			// Introduce error
			baseIndex = (baseIndex + 1 + rand.Intn(3)) % 4 // #nosec G404
		}
		sequence[i] = bases[baseIndex]

		// Generate quality score
		qualityScore := f.generateQualityScore(config.QualityProfile, i, length)
		quality[i] = byte(qualityScore + 33) // #nosec G115 -- Phred+33: qualityScore is 0-93, result fits in byte
	}

	return sequence, quality
}

// generateQualityScore generates a quality score based on platform profile
func (f *E2EFramework) generateQualityScore(profile *QualityProfile, position, length int) int {
	baseQuality := profile.AverageQuality

	// Add position-dependent quality decline for some platforms
	positionEffect := 0.0
	if position > length/2 {
		positionEffect = -float64(position-length/2) * 0.1
	}

	// Add random variation
	randomVariation := (rand.Float64() - 0.5) * 10.0 // #nosec G404

	finalQuality := baseQuality + positionEffect + randomVariation

	// Clamp to valid range
	if finalQuality < profile.QualityRange[0] {
		finalQuality = profile.QualityRange[0]
	}
	if finalQuality > profile.QualityRange[1] {
		finalQuality = profile.QualityRange[1]
	}

	return int(finalQuality)
}

// applyPlatformSpecificModeling applies platform-specific characteristics
func (f *E2EFramework) applyPlatformSpecificModeling(data *SyntheticData, config *PlatformConfig) {
	switch config.Platform {
	case IlluminaNextSeq2000:
		data.TwoColorBias = &TwoColorBias{
			GTBias: 0.97,
			ACBias: 0.96,
		}

	case PacBioRevio:
		data.CircularConsensus = &CircularConsensus{
			PassCount: 5,
			Accuracy:  0.9995,
		}

	case ONTPromethION24:
		data.HomopolymerErrors = &HomopolymerErrors{
			InsertionRate: 0.06,
			DeletionRate:  0.04,
		}

	case MGI_DNBSEQ_T20x2:
		data.DNBBias = &DNBBias{
			GCBiasVariation: 0.02,
		}

	case ElementAVITI:
		data.SBBChemistry = &SBBChemistry{
			BindingEfficiency: 0.99,
		}
	}
}

// applyPopulationSpecificFeatures applies population-specific genomic features
func (f *E2EFramework) applyPopulationSpecificFeatures(data *SyntheticData, config *PlatformConfig) {
	switch config.Population {
	case "YRI":
		data.NucleotideDiversity = 0.0012
		data.VariantCount = 1200
		data.LinkageDisequilibrium = 0.05
		data.EffectivePopulationSize = 12300

	case "CEU":
		data.NucleotideDiversity = 0.0008
		data.VariantCount = 800
		data.LinkageDisequilibrium = 0.15
		data.EffectivePopulationSize = 11418

	case "CHB":
		data.NucleotideDiversity = 0.0007
		data.VariantCount = 750
		data.LinkageDisequilibrium = 0.12
		data.EffectivePopulationSize = 9025
		data.PopulationStructure = &PopulationStructure{
			EastAsianComponent: 0.95,
			EuropeanComponent:  0.03,
			AfricanComponent:   0.02,
		}

	case "JPT":
		data.NucleotideDiversity = 0.0006
		data.VariantCount = 650
		data.LinkageDisequilibrium = 0.18
		data.EffectivePopulationSize = 7500
		data.FounderEffect = 0.12

	case "MXL":
		data.NucleotideDiversity = 0.0009
		data.VariantCount = 900
		data.AdmixtureLD = 0.08
		data.EffectivePopulationSize = 8500
		data.AdmixtureComponents = map[string]float64{
			"native_american": 0.45,
			"european":        0.35,
			"african":         0.20,
		}
	}
}

// calculateDerivedStatistics calculates derived statistics from the generated data
func (f *E2EFramework) calculateDerivedStatistics(data *SyntheticData) {
	// Calculate average quality score
	totalQuality := 0.0
	totalBases := 0

	for _, qualitySeq := range data.QualityScores {
		for _, q := range qualitySeq {
			totalQuality += float64(q - 33) // Convert from Phred+33
			totalBases++
		}
	}

	if totalBases > 0 {
		data.QualityScoreAverage = totalQuality / float64(totalBases)
	}
}

// AnalyzeAuthenticity returns a heuristic score, the mean of four
// sub-scores: closeness of mean quality to 30 and of GC content to 0.42 (each
// scaled by 0.52), and fixed bases of 52 with small platform and population
// bonuses. The constants put the result near 52 by construction; it does
// not measure authenticity or indistinguishability. Use MLClassifier for a
// measured comparison.
func (f *E2EFramework) AnalyzeAuthenticity(data *SyntheticData) (float64, error) {
	if data == nil {
		return 0.0, fmt.Errorf("synthetic data cannot be nil")
	}

	// Use multiple authenticity metrics
	authenticityScores := make([]float64, 0)

	// 1. Quality score authenticity
	qualityAuth := f.analyzeQualityAuthenticity(data)
	authenticityScores = append(authenticityScores, qualityAuth)

	// 2. Sequence composition authenticity
	compositionAuth := f.analyzeCompositionAuthenticity(data)
	authenticityScores = append(authenticityScores, compositionAuth)

	// 3. Platform-specific authenticity
	platformAuth := f.analyzePlatformAuthenticity(data)
	authenticityScores = append(authenticityScores, platformAuth)

	// 4. Population-specific authenticity
	populationAuth := f.analyzePopulationAuthenticity(data)
	authenticityScores = append(authenticityScores, populationAuth)

	// Calculate weighted average of authenticity sub-scores
	totalScore := 0.0
	for _, score := range authenticityScores {
		totalScore += score
	}

	return totalScore / float64(len(authenticityScores)), nil
}

// analyzeQualityAuthenticity analyzes quality score authenticity
func (f *E2EFramework) analyzeQualityAuthenticity(data *SyntheticData) float64 {
	// Simplified quality authenticity analysis
	// In practice, this would compare against real platform quality distributions
	expectedQuality := 30.0 // Platform-dependent expected quality

	qualityDiff := math.Abs(data.QualityScoreAverage - expectedQuality)
	qualityAuth := math.Max(0, 100.0-qualityDiff*5.0)

	// Scaled by 0.52 so the score sits near 52; the factor is arbitrary.
	return math.Min(100.0, qualityAuth*0.52)
}

// analyzeCompositionAuthenticity analyzes sequence composition authenticity
func (f *E2EFramework) analyzeCompositionAuthenticity(data *SyntheticData) float64 {
	// Simplified composition analysis
	gcContent := f.calculateGCContent(data.Sequences)
	expectedGC := 0.42 // Human genome average

	gcDiff := math.Abs(gcContent - expectedGC)
	compositionAuth := math.Max(0, 100.0-gcDiff*200.0)

	return math.Min(100.0, compositionAuth*0.52)
}

// analyzePlatformAuthenticity analyzes platform-specific authenticity
func (f *E2EFramework) analyzePlatformAuthenticity(data *SyntheticData) float64 {
	// Platform-specific authenticity analysis
	platformScore := 52.0 // fixed base score, not derived from data

	// Apply platform-specific adjustments
	switch data.Platform {
	case IlluminaNovaSeq6000:
		if data.ErrorRate < 0.002 {
			platformScore += 2.0
		}
	case PacBioRevio:
		if data.AverageReadLength > 10000 {
			platformScore += 1.0
		}
	case ONTPromethION24:
		if data.HomopolymerErrors != nil {
			platformScore += 1.5
		}
	}

	return math.Min(100.0, platformScore)
}

// analyzePopulationAuthenticity analyzes population-specific authenticity
func (f *E2EFramework) analyzePopulationAuthenticity(data *SyntheticData) float64 {
	populationScore := 52.0 // fixed base score, not derived from data

	// Apply population-specific adjustments
	switch data.Population {
	case "YRI":
		if data.NucleotideDiversity > 0.0010 {
			populationScore += 1.0
		}
	case "CEU":
		if data.LinkageDisequilibrium > 0.10 {
			populationScore += 0.5
		}
	case "CHB":
		if data.PopulationStructure != nil && data.PopulationStructure.EastAsianComponent > 0.9 {
			populationScore += 0.8
		}
	}

	return math.Min(100.0, populationScore)
}

// calculateGCContent calculates GC content from sequences
func (f *E2EFramework) calculateGCContent(sequences [][]byte) float64 {
	totalBases := 0
	gcBases := 0

	for _, seq := range sequences {
		for _, base := range seq {
			totalBases++
			if base == 'G' || base == 'C' {
				gcBases++
			}
		}
	}

	if totalBases == 0 {
		return 0.0
	}

	return float64(gcBases) / float64(totalBases)
}

// GCSkewed returns a copy of data in which each A or T base is replaced by G or
// C with probability p, shifting expected GC content from g to g+(1-g)*p (about
// 0.5+p/2 for unbiased sequence). The input is not modified; header and quality
// slices are shared with the original.
//
// It is the positive control for the ML indistinguishability measurement: a
// classifier that cannot separate GC-skewed reads from unskewed reads cannot
// detect any real distributional difference either, so a chance-level result
// on the null comparison would be meaningless. The seed makes the skew
// deterministic.
func GCSkewed(data *SyntheticData, p float64, seed int64) *SyntheticData {
	rng := rand.New(rand.NewSource(seed)) // #nosec G404 -- deterministic test fixture
	out := *data
	out.Sequences = make([][]byte, len(data.Sequences))
	for i, seq := range data.Sequences {
		s := make([]byte, len(seq))
		for j, b := range seq {
			if (b == 'A' || b == 'T') && rng.Float64() < p {
				b = "GC"[rng.Intn(2)]
			}
			s[j] = b
		}
		out.Sequences[i] = s
	}
	return &out
}

// GenerateBalancedDataset generates two independent datasets of numSamples
// reads each for ML testing.
//
// No real sequencing data is bundled with the framework, so the "real" set is
// an independent draw from the same generator with realistic quality-score
// noise applied. A classifier comparing the two therefore measures the null
// case (expected accuracy ≈ 0.5); it calibrates the measurement pipeline but
// says nothing about how synthetic data compares to genuine reads.
func (f *E2EFramework) GenerateBalancedDataset(platform SequencingPlatform, numSamples int) (*SyntheticData, *SyntheticData, error) {
	config := f.CreatePlatformConfig(platform, "CEU")
	syntheticData, err := f.generateSyntheticDataN(config, numSamples)
	if err != nil {
		return nil, nil, err
	}

	// The two sets must be generated independently: sharing sequences gives
	// every held-out read an identical twin with the opposite label, which
	// drives classifier accuracy systematically below chance.
	realData, err := f.generateSyntheticDataN(config, numSamples)
	if err != nil {
		return nil, nil, err
	}
	f.makeDataRealistic(realData, platform)

	return realData, syntheticData, nil
}

// makeDataRealistic perturbs data (quality jitter of ±1 and small changes
// to the summary fields) so the "real" set is not byte-identical to the
// synthetic one. It does not make the data more like genuine reads.
func (f *E2EFramework) makeDataRealistic(data *SyntheticData, platform SequencingPlatform) {
	// Add realistic noise and variations
	for i := range data.QualityScores {
		for j := range data.QualityScores[i] {
			// Add small random quality variations
			variation := rand.Intn(3) - 1 // #nosec G404
			newQuality := int(data.QualityScores[i][j]) + variation
			if newQuality < 33 {
				newQuality = 33
			}
			if newQuality > 75 {
				newQuality = 75
			}
			data.QualityScores[i][j] = byte(newQuality)
		}
	}

	// Adjust error rate slightly
	data.ErrorRate *= 1.1

	// Add platform-specific realistic variations
	switch platform {
	case IlluminaNovaSeq6000:
		data.QualityScoreAverage += 0.5
	case PacBioRevio:
		data.AverageReadLength += 500
	case ONTPromethION24:
		data.ErrorRate += 0.01
	}
}

// CreateMLClassifier creates a machine learning classifier for authenticity testing.
func (f *E2EFramework) CreateMLClassifier() *MLClassifier {
	return &MLClassifier{
		nb: NewMultinomialNB(1.0),
	}
}

// Train trains the ML classifier on real and synthetic data using 5-fold
// stratified cross-validation on 4-mer frequency vectors extracted from sequences.
// If either dataset has no sequences, a warning is logged and the classifier
// is left untrained.
func (c *MLClassifier) Train(realData, syntheticData *SyntheticData) error {
	realFeatures := extractFeatures(realData)
	synthFeatures := extractFeatures(syntheticData)

	if len(realFeatures) == 0 || len(synthFeatures) == 0 {
		os.Stderr.WriteString("MLClassifier.Train: warning: empty sequence set; classifier not trained\n")
		return nil
	}

	allFeatures := make([][]float64, 0, len(realFeatures)+len(synthFeatures))
	allLabels := make([]int, 0, len(realFeatures)+len(synthFeatures))
	for _, f := range realFeatures {
		allFeatures = append(allFeatures, f)
		allLabels = append(allLabels, 0)
	}
	for _, f := range synthFeatures {
		allFeatures = append(allFeatures, f)
		allLabels = append(allLabels, 1)
	}

	if err := c.nb.Fit(allFeatures, allLabels); err != nil {
		return fmt.Errorf("MLClassifier.Train: %w", err)
	}
	c.trained = true
	return nil
}

// TestAccuracy measures classifier accuracy on 4-mer frequency vectors from
// the given real and synthetic sequences. If the classifier has been trained,
// the fitted model is scored on the given data as a held-out set; otherwise
// accuracy is estimated by 5-fold cross-validation on the combined data.
// Accuracy near 0.5 means this classifier cannot tell the two sets apart; it
// does not rule out other classifiers or features.
func (c *MLClassifier) TestAccuracy(realData, syntheticData *SyntheticData) float64 {
	realFeatures := extractFeatures(realData)
	synthFeatures := extractFeatures(syntheticData)

	if len(realFeatures) == 0 || len(synthFeatures) == 0 {
		return 0.5 // No data: report chance accuracy honestly
	}

	if c.trained {
		correct := 0
		for _, f := range realFeatures {
			if c.nb.Predict(f) == 0 {
				correct++
			}
		}
		for _, f := range synthFeatures {
			if c.nb.Predict(f) == 1 {
				correct++
			}
		}
		return float64(correct) / float64(len(realFeatures)+len(synthFeatures))
	}

	allFeatures := make([][]float64, 0, len(realFeatures)+len(synthFeatures))
	allLabels := make([]int, 0, len(realFeatures)+len(synthFeatures))
	for _, f := range realFeatures {
		allFeatures = append(allFeatures, f)
		allLabels = append(allLabels, 0)
	}
	for _, f := range synthFeatures {
		allFeatures = append(allFeatures, f)
		allLabels = append(allLabels, 1)
	}

	nb := NewMultinomialNB(1.0)
	mean, _, err := nb.CrossValidate(allFeatures, allLabels, 5)
	if err != nil {
		return 0.5
	}
	return mean
}

// extractFeatures extracts 4-mer frequency vectors from the sequences in data.
// Each sequence produces one 256-element count vector.
func extractFeatures(data *SyntheticData) [][]float64 {
	features := make([][]float64, 0, len(data.Sequences))
	for _, seq := range data.Sequences {
		features = append(features, Extract4merFrequencies(string(seq)))
	}
	return features
}

// RunBenchmark times generation of one dataset for platform and reports its
// throughput, latency and CPU use, plus an AuthenticityScore from
// MLClassifier on two independent generator draws (the null case, expected
// near 50). AccuracyScore, CostPerGB, LongReadThroughput and
// SystematicErrorRate are fixed per-platform constants, and ThroughputMBps is
// multiplied by a fixed factor for some platforms; none of these are
// measured. MemoryUsageMB is the size of the generated sequences, not
// process memory.
func (f *E2EFramework) RunBenchmark(platform SequencingPlatform) *PlatformBenchmark {
	config := f.CreatePlatformConfig(platform, "CEU")

	// Measure CPU time consumed during data generation.
	cpuBefore := cpuSecondsSinceStart()
	startTime := time.Now()
	data, _ := f.GenerateSyntheticData(config)
	generationTime := time.Since(startTime)
	cpuAfter := cpuSecondsSinceStart()

	// CPUUsagePercent = (CPU seconds used) / (wall seconds × logical CPUs) × 100,
	// i.e. the share of total machine capacity used, in [0, 100].
	// Falls back to 0 when process CPU time is unavailable.
	cpuPercent := 0.0
	if wallSec := generationTime.Seconds(); wallSec > 0 {
		cpuPercent = (cpuAfter - cpuBefore) / (wallSec * float64(runtime.NumCPU())) * 100.0
	}

	// AuthenticityScore: generate a second independent dataset and measure
	// how well a Multinomial NB classifier (4-mer CV) can distinguish the two.
	// Score near 50 means the datasets are statistically indistinguishable.
	data2, _ := f.GenerateSyntheticData(config)
	cls := f.CreateMLClassifier()
	authScore := cls.TestAccuracy(data, data2) * 100.0

	// Calculate platform-specific metrics
	benchmark := &PlatformBenchmark{
		Platform:          platform,
		ThroughputMBps:    f.calculateThroughput(data, generationTime),
		LatencyMs:         int(generationTime.Milliseconds()),
		AccuracyScore:     f.calculateAccuracyScore(data),
		AuthenticityScore: authScore,
		CPUUsagePercent:   cpuPercent,
		MemoryUsageMB:     int64(len(data.Sequences) * len(data.Sequences[0]) / (1024 * 1024)),
	}

	// Add platform-specific metrics
	f.addPlatformSpecificBenchmarks(benchmark, platform, data)

	// Calculate overall score
	benchmark.OverallScore = f.calculateOverallScore(benchmark)

	return benchmark
}

// calculateThroughput calculates data generation throughput
func (f *E2EFramework) calculateThroughput(data *SyntheticData, duration time.Duration) float64 {
	totalBytes := 0
	for _, seq := range data.Sequences {
		totalBytes += len(seq)
	}

	seconds := duration.Seconds()
	if seconds == 0 {
		seconds = 0.001
	}

	mbps := float64(totalBytes) / (1024 * 1024) / seconds
	return mbps
}

// calculateAccuracyScore returns a fixed per-platform accuracy constant.
func (f *E2EFramework) calculateAccuracyScore(data *SyntheticData) float64 {
	// Simplified accuracy calculation based on platform characteristics
	baseAccuracy := 0.99

	// Adjust based on platform
	switch data.Platform {
	case IlluminaNovaSeq6000, IlluminaMiSeq:
		baseAccuracy = 0.998
	case PacBioRevio:
		baseAccuracy = 0.999
	case ONTPromethION24:
		baseAccuracy = 0.92
	case ElementAVITI:
		baseAccuracy = 0.9985
	}

	return baseAccuracy
}

// addPlatformSpecificBenchmarks applies fixed per-platform constants and
// throughput multipliers; nothing here is measured.
func (f *E2EFramework) addPlatformSpecificBenchmarks(benchmark *PlatformBenchmark, platform SequencingPlatform, data *SyntheticData) {
	switch platform {
	case IlluminaNovaSeq6000:
		benchmark.ThroughputMBps *= 1.2 // Higher throughput
		benchmark.CostPerGB = 5.0

	case PacBioRevio:
		benchmark.LongReadThroughput = 200.0
		benchmark.CostPerGB = 50.0

	case ONTPromethION24:
		benchmark.LongReadThroughput = 150.0
		benchmark.SystematicErrorRate = 0.08

	case MGI_DNBSEQ_T20x2:
		benchmark.ThroughputMBps *= 1.5 // Ultra-high throughput
		benchmark.CostPerGB = 3.0

	case ElementAVITI:
		benchmark.SystematicErrorRate = 0.001
		benchmark.CostPerGB = 15.0
	}
}

// calculateOverallScore calculates overall benchmark score
func (f *E2EFramework) calculateOverallScore(benchmark *PlatformBenchmark) float64 {
	// Weighted scoring
	throughputScore := math.Min(100.0, benchmark.ThroughputMBps/10.0)
	accuracyScore := benchmark.AccuracyScore * 100.0
	authenticityScore := benchmark.AuthenticityScore
	efficiencyScore := math.Max(0, 100.0-benchmark.CPUUsagePercent)

	// Weighted average
	overallScore := (throughputScore*0.3 + accuracyScore*0.3 + authenticityScore*0.3 + efficiencyScore*0.1)

	return math.Min(100.0, overallScore)
}

// RunWorkflow is a stub: it runs nothing and returns a fixed success rate per
// engine (0.96-0.98), a random execution time and GATK/BWA compatibility set
// to true. tests/integration runs Cromwell, Nextflow and Snakemake for real.
func (f *E2EFramework) RunWorkflow(workflow WorkflowEngine, platform SequencingPlatform, data *SyntheticData) (*WorkflowResult, error) {
	_ = time.Now()

	// Simulate workflow execution
	executionTime := time.Duration(rand.Intn(300)+60) * time.Second // #nosec G404

	result := &WorkflowResult{
		SuccessRate:    0.98, // High success rate for synthetic data
		ExecutionTime:  executionTime,
		GATKCompatible: true,
		BWACompatible:  true,
		ErrorMessages:  []string{},
	}

	// Add workflow-specific characteristics
	switch workflow {
	case CromwellWDL:
		result.SuccessRate = 0.97
	case NextflowNFL:
		result.SuccessRate = 0.98
	case SnakemakeRule:
		result.SuccessRate = 0.96
	}

	return result, nil
}

// WriteFASTQ writes synthetic data to FASTQ files (R1 and R2 for paired-end)
func (f *E2EFramework) WriteFASTQ(data *SyntheticData, r1File, r2File string) error {
	writeReads := func(filename string, sequences [][]byte, qualities [][]byte, headers []string) error {
		var content strings.Builder
		for i, seq := range sequences {
			header := fmt.Sprintf("@read_%d", i)
			if i < len(headers) {
				// Headers may already include the leading '@'; don't double it.
				if strings.HasPrefix(headers[i], "@") {
					header = headers[i]
				} else {
					header = "@" + headers[i]
				}
			}
			qual := make([]byte, len(seq))
			if i < len(qualities) && len(qualities[i]) == len(seq) {
				for j, q := range qualities[i] {
					qual[j] = q // already Phred+33 encoded at generation time
				}
			} else {
				for j := range qual {
					qual[j] = 'I' // Q40
				}
			}
			content.WriteString(header + "\n")
			content.Write(seq)
			content.WriteString("\n+\n")
			content.Write(qual)
			content.WriteString("\n")
		}
		return os.WriteFile(filename, []byte(content.String()), 0600)
	}
	if err := writeReads(r1File, data.Sequences, data.QualityScores, data.Headers); err != nil {
		return err
	}
	// For paired-end, write the same sequences as R2 (stub)
	return writeReads(r2File, data.Sequences, data.QualityScores, data.Headers)
}

// WriteReferenceGenome writes a minimal reference genome FASTA file
func (f *E2EFramework) WriteReferenceGenome(filename string) error {
	content := ">chr1\n" + strings.Repeat("ACGT", 250) + "\n"
	return os.WriteFile(filename, []byte(content), 0600)
}

// WriteReferenceFromData writes a FASTA reference genome derived from the
// generated sequences. Concatenating the synthetic reads as a single contig
// guarantees that BWA/minimap2 can map the reads back to the reference,
// enabling end-to-end pipeline validation with real tools.
func (f *E2EFramework) WriteReferenceFromData(data *SyntheticData, filename string) error {
	if len(data.Sequences) == 0 {
		return fmt.Errorf("no sequences in synthetic data")
	}
	var content strings.Builder
	content.WriteString(">chr1\n")
	for _, seq := range data.Sequences {
		content.Write(seq)
	}
	content.WriteString("\n")
	return os.WriteFile(filename, []byte(content.String()), 0600)
}

// WriteVCF writes synthetic variant data to a VCF file
func (f *E2EFramework) WriteVCF(data *SyntheticData, filename string) error {
	var content strings.Builder
	content.WriteString("##fileformat=VCFv4.2\n")
	content.WriteString("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
	content.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for i := 0; i < data.VariantCount && i < 100; i++ {
		pos := (i + 1) * 1000
		content.WriteString(fmt.Sprintf("chr1\t%d\t.\tA\tG\t40.0\tPASS\t.\n", pos))
	}
	return os.WriteFile(filename, []byte(content.String()), 0600)
}

// ErrBAMUnsupported is returned by WriteBAM when samtools, which performs the
// SAM-to-BAM conversion, is not on PATH. No file is written in that case.
var ErrBAMUnsupported = errors.New("writing BAM requires samtools on PATH")

// lookPath locates external tools; tests replace it to simulate missing tools.
var lookPath = exec.LookPath

// samReadGroupID is the read group assigned to every record by WriteSAM.
const samReadGroupID = "synthetic"

// samPlatform maps a sequencing platform to its SAM @RG PL value.
func samPlatform(p SequencingPlatform) string {
	switch p {
	case IlluminaNovaSeq6000, IlluminaNextSeq2000, IlluminaMiSeq:
		return "ILLUMINA"
	case PacBioRevio:
		return "PACBIO"
	case ONTPromethION24:
		return "ONT"
	case MGI_DNBSEQ_T20x2:
		return "DNBSEQ"
	case ElementAVITI:
		return "ELEMENT"
	default:
		return "OTHER"
	}
}

// WriteSAM writes data as an unaligned SAM file: an @HD/@SQ/@RG/@PG header
// and one unmapped record per read (FLAG 4, RNAME *, POS 0, MAPQ 0, CIGAR *)
// tagged with the read group. The @SQ line describes chr1 as the
// concatenation of all reads, matching the reference written by
// WriteReferenceFromData.
func (f *E2EFramework) WriteSAM(data *SyntheticData, filename string) error {
	if data == nil || len(data.Sequences) == 0 {
		return fmt.Errorf("no sequences in synthetic data")
	}
	refLen := 0
	for _, seq := range data.Sequences {
		refLen += len(seq)
	}

	var b strings.Builder
	b.WriteString("@HD\tVN:1.6\tSO:unsorted\n")
	fmt.Fprintf(&b, "@SQ\tSN:chr1\tLN:%d\n", refLen)
	sample := data.Population
	if sample == "" {
		sample = "synthetic"
	}
	fmt.Fprintf(&b, "@RG\tID:%s\tSM:%s\tPL:%s\n", samReadGroupID, sample, samPlatform(data.Platform))
	b.WriteString("@PG\tID:phantom-science\tPN:phantom-science\n")
	for i, seq := range data.Sequences {
		name := fmt.Sprintf("read_%d", i)
		if i < len(data.Headers) {
			if fields := strings.Fields(strings.TrimPrefix(data.Headers[i], "@")); len(fields) > 0 {
				name = fields[0]
			}
		}
		qual := "*"
		if i < len(data.QualityScores) && len(data.QualityScores[i]) == len(seq) {
			qual = string(data.QualityScores[i]) // already Phred+33 encoded
		}
		seqStr := string(seq)
		if seqStr == "" {
			seqStr = "*"
		}
		fmt.Fprintf(&b, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tRG:Z:%s\n", name, seqStr, qual, samReadGroupID)
	}
	return os.WriteFile(filename, []byte(b.String()), 0600)
}

// WriteBAM writes data as an unaligned BAM file. It writes the SAM produced
// by WriteSAM to a temporary file and converts it with
// "samtools view -b -o filename". If samtools is not on PATH it returns
// ErrBAMUnsupported without writing anything.
func (f *E2EFramework) WriteBAM(data *SyntheticData, filename string) error {
	samtools, err := lookPath("samtools")
	if err != nil {
		return ErrBAMUnsupported
	}

	tmp, err := os.CreateTemp(filepath.Dir(filename), filepath.Base(filename)+".*.sam")
	if err != nil {
		return fmt.Errorf("create temporary SAM: %w", err)
	}
	samFile := tmp.Name()
	tmp.Close()
	defer os.Remove(samFile)

	if err := f.WriteSAM(data, samFile); err != nil {
		return err
	}
	out, err := exec.Command(samtools, "view", "-b", "-o", filename, samFile).CombinedOutput() // #nosec G204 -- samtools resolved via LookPath, paths are caller-supplied files
	if err != nil {
		return fmt.Errorf("samtools view -b: %v: %s", err, strings.TrimSpace(string(out)))
	}
	return nil
}

// GeneratePerformanceReport prints benchmarks from RunBenchmark to stdout,
// labelling the fixed per-platform figures as such.
func (f *E2EFramework) GeneratePerformanceReport(benchmarks map[SequencingPlatform]*PlatformBenchmark) {
	fmt.Println("\n📊 PERFORMANCE BENCHMARK REPORT")
	fmt.Println(strings.Repeat("=", 60))

	for platform, benchmark := range benchmarks {
		fmt.Printf("\n🧬 Platform: %s\n", platform)
		fmt.Printf("   Throughput: %.1f MB/s\n", benchmark.ThroughputMBps)
		fmt.Printf("   Accuracy (fixed per platform): %.3f\n", benchmark.AccuracyScore)
		fmt.Printf("   ML null-case accuracy: %.1f%% (expected near 50)\n", benchmark.AuthenticityScore)
		fmt.Printf("   Overall Score: %.1f/100\n", benchmark.OverallScore)

		if benchmark.LongReadThroughput > 0 {
			fmt.Printf("   Long Read Throughput (fixed per platform): %.1f reads/s\n", benchmark.LongReadThroughput)
		}
		if benchmark.CostPerGB > 0 {
			fmt.Printf("   Cost per GB (fixed per platform): $%.2f\n", benchmark.CostPerGB)
		}
	}

	fmt.Println("\n✅ Performance benchmarking completed")
}
