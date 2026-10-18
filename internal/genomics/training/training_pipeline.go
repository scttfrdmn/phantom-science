// Package training wires the HMM k-mer, GC-quality and population-genetics
// trainers into a five-phase pipeline (data preparation, training,
// validation, performance testing, final threshold checks), and downloads
// reference and validation files from AWS Open Data.
//
// Several validation stages are stubs: RunDownstreamTesting records a fixed
// result without running any tool, the platform validation datasets are
// uniform random sequences, and without trained models the statistical and
// performance stages store fixed placeholder values.
package training

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"runtime"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/population"
	"github.com/scttfrdmn/phantom-science/internal/genomics/quality"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/scttfrdmn/phantom-science/internal/genomics/validation"
	"github.com/scttfrdmn/phantom-science/internal/security"
)

// TrainingPipeline runs model training and the validation stages in order.
type TrainingPipeline struct {
	config       *TrainingConfig
	dataManager  *DataManager
	modelTrainer *ModelTrainer
	validator    *ValidationManager
	logger       *PipelineLogger
	metrics      *TrainingMetrics
	mu           sync.RWMutex
}

// TrainingConfig contains configuration for the training pipeline
type TrainingConfig struct {
	// Data paths
	ReferenceGenomePath string
	QualityDataPath     string // directory containing <platform>.json quality data files
	PopulationDataPath  string
	ValidationDataPath  string
	OutputPath          string
	CachePath           string

	// Training parameters
	KmerSizes            []int
	Populations          []string
	Platforms            []string
	ValidationSplit      float64
	CrossValidationFolds int

	// Performance parameters
	MaxMemoryGB int
	NumWorkers  int
	BatchSize   int
	EnableGPU   bool

	// Validation thresholds
	JSdivergenceThreshold      float64
	ExpertAccuracyThreshold    float64
	DownstreamSuccessThreshold float64

	// Training controls
	EnableEarlyStopping bool
	MaxTrainingTime     time.Duration
	CheckpointInterval  time.Duration
}

// DataManager handles data loading and preprocessing
type DataManager struct {
	config          *TrainingConfig
	referenceCache  map[string]*reference.ReferenceGenome
	populationCache map[string]*population.PopulationPanel
	validationCache map[string]*ValidationDataset
	mu              sync.RWMutex
}

// ModelTrainer handles training of all model components
type ModelTrainer struct {
	config         *TrainingConfig
	hmmTrainer     *HMMTrainer
	qualityTrainer *QualityTrainer
	popTrainer     *PopulationTrainer
	checkpointer   *ModelCheckpointer
}

// ValidationManager runs the statistical, expert, downstream and
// performance validation stages.
type ValidationManager struct {
	config               *TrainingConfig
	statisticalValidator *validation.StatisticalValidator
	expertTesting        *validation.ExpertTestingFramework
	performanceTester    *PerformanceTester
	results              *ValidationResults
	// seqGenerators holds trained models wired in after the training phase.
	// Keys use the same "<k>mer_<pop>" format as HMMTrainer.models.
	seqGenerators map[string]authenticity.SequenceGenerator
}

// SetSequenceGenerators provides the ValidationManager with trained HMM models
// so RunStatisticalValidation and TestGeneration* can use real generation.
func (vm *ValidationManager) SetSequenceGenerators(gens map[string]authenticity.SequenceGenerator) {
	vm.seqGenerators = gens
}

// TrainingMetrics tracks training progress and performance
type TrainingMetrics struct {
	StartTime          time.Time
	TrainingDuration   time.Duration
	ModelsTrained      int
	ValidationScores   map[string]float64
	PerformanceMetrics map[string]float64
	ErrorLog           []TrainingError
	mu                 sync.RWMutex
}

// ValidationDataset represents a dataset for validation
type ValidationDataset struct {
	Platform  string
	DataType  string
	FilePath  string
	Sequences [][]byte
	Qualities [][]byte
	Metadata  map[string]interface{}
	Size      int
	Checksum  string
}

// TrainingError represents a training error
type TrainingError struct {
	Timestamp time.Time
	Component string
	Error     error
	Severity  string
	Context   map[string]interface{}
}

// ValidationResults stores comprehensive validation results
type ValidationResults struct {
	OverallScore        float64
	JSdivergenceResults map[string]float64
	ExpertTestResults   *validation.ExpertTestingResults
	PerformanceResults  *PerformanceResults
	DownstreamResults   *DownstreamResults
	PassesValidation    bool
	Recommendations     []string
}

// NewTrainingPipeline creates a new training pipeline
func NewTrainingPipeline(config *TrainingConfig) *TrainingPipeline {
	return &TrainingPipeline{
		config:       config,
		dataManager:  NewDataManager(config),
		modelTrainer: NewModelTrainer(config),
		validator:    NewValidationManager(config),
		logger:       NewPipelineLogger(config.OutputPath),
		metrics:      NewTrainingMetrics(),
	}
}

// NewDataManager creates a new data manager
func NewDataManager(config *TrainingConfig) *DataManager {
	return &DataManager{
		config:          config,
		referenceCache:  make(map[string]*reference.ReferenceGenome),
		populationCache: make(map[string]*population.PopulationPanel),
		validationCache: make(map[string]*ValidationDataset),
	}
}

// LoadReferenceGenome loads a reference genome assembly
func (dm *DataManager) LoadReferenceGenome(assembly string) error {
	dm.mu.Lock()
	defer dm.mu.Unlock()

	// Check if already cached
	if _, exists := dm.referenceCache[assembly]; exists {
		return nil
	}

	// Determine base directory for the assembly
	baseDir := filepath.Join(dm.config.ReferenceGenomePath, assembly)

	// Create reference genome instance
	refGenome := reference.NewReferenceGenome(baseDir)

	// Validate that the reference genome has data
	chromosomes := refGenome.GetChromosomeNames()
	if len(chromosomes) == 0 {
		return fmt.Errorf("no chromosomes found for assembly %s", assembly)
	}

	// Test loading one chromosome to validate format
	testChr, err := refGenome.LoadChromosome(chromosomes[0])
	if err != nil {
		return fmt.Errorf("failed to load test chromosome for %s: %w", assembly, err)
	}
	if testChr.GetLength() == 0 {
		return fmt.Errorf("test chromosome has zero length for assembly %s", assembly)
	}

	// Cache the reference genome
	dm.referenceCache[assembly] = refGenome

	log.Printf("Loaded reference genome %s with %d chromosomes", assembly, len(chromosomes))
	return nil
}

// LoadPopulationPanel loads a population panel
func (dm *DataManager) LoadPopulationPanel(populationName string) error {
	dm.mu.Lock()
	defer dm.mu.Unlock()

	// Check if already cached
	if _, exists := dm.populationCache[populationName]; exists {
		return nil
	}

	// Create population panel with default sample size
	sampleSize := 100 // Default sample size
	if dm.config.ValidationSplit > 0 {
		sampleSize = int(100.0 / dm.config.ValidationSplit)
	}

	// The panel's demography (Ne, growth, bottleneck, split from its parent
	// population) comes from the presets in internal/genomics/population.
	// Codes without their own preset (LWK, GBR, JPT, ...) use the closest
	// modelled population as a proxy; anything else gets a constant-size
	// default model.
	panel := population.NewPopulationPanel(populationName, sampleSize)
	if !panel.HasPreset() {
		log.Printf("No demographic preset for population %s; using a constant-size default model", populationName)
	} else if panel.PresetName() != populationName {
		log.Printf("Population %s simulated with the %s demographic preset as a proxy", populationName, panel.PresetName())
	}

	// Load simulation data
	regionLength := 1000000 // 1MB default region
	if err := panel.LoadFromSimulation(regionLength); err != nil {
		return fmt.Errorf("failed to load simulation data for %s: %w", populationName, err)
	}

	// Cache the population panel
	dm.populationCache[populationName] = panel

	stats := panel.GetPopulationStats()
	log.Printf("Loaded population panel %s with %v variants", populationName, stats["variants"])
	return nil
}

// LoadValidationDataset builds (and caches) a synthetic validation dataset
// for a platform; see generateValidationData.
func (dm *DataManager) LoadValidationDataset(platform string) error {
	dm.mu.Lock()
	defer dm.mu.Unlock()

	// Check if already cached
	if _, exists := dm.validationCache[platform]; exists {
		return nil
	}

	// Create validation dataset
	dataset := &ValidationDataset{
		Platform: platform,
		DataType: "simulated",
		FilePath: filepath.Join(dm.config.ValidationDataPath, platform),
		Metadata: make(map[string]interface{}),
	}

	// Generate synthetic validation data
	if err := dm.generateValidationData(dataset); err != nil {
		return fmt.Errorf("failed to generate validation data for %s: %w", platform, err)
	}

	// Cache the validation dataset
	dm.validationCache[platform] = dataset

	log.Printf("Loaded validation dataset for %s with %d sequences", platform, dataset.Size)
	return nil
}

// generateValidationData fills dataset with uniform random bases and raw
// quality values (30, or 10 at the platform's error rate). It is not real
// sequencing data.
func (dm *DataManager) generateValidationData(dataset *ValidationDataset) error {
	// Platform-specific sequence generation parameters
	var sequenceLength int
	var numSequences int
	var errorRate float64

	switch dataset.Platform {
	case "illumina", "novaseq", "nextseq":
		sequenceLength = 150
		numSequences = 1000
		errorRate = 0.001
	case "pacbio", "sequel":
		sequenceLength = 10000
		numSequences = 100
		errorRate = 0.1
	case "nanopore", "minion":
		sequenceLength = 5000
		numSequences = 200
		errorRate = 0.05
	default:
		sequenceLength = 150
		numSequences = 1000
		errorRate = 0.001
	}

	// Generate sequences
	dataset.Sequences = make([][]byte, numSequences)
	dataset.Qualities = make([][]byte, numSequences)

	for i := 0; i < numSequences; i++ {
		// Generate random sequence
		sequence := make([]byte, sequenceLength)
		qualities := make([]byte, sequenceLength)

		for j := 0; j < sequenceLength; j++ {
			// Random base
			bases := []byte{'A', 'T', 'G', 'C'}
			sequence[j] = bases[rand.IntN(4)] // #nosec G404

			// Platform-specific quality score
			baseQuality := 30.0             // Phred score 30
			if rand.Float64() < errorRate { // #nosec G404
				baseQuality = 10.0 // Lower quality for errors
			}
			qualities[j] = byte(baseQuality)
		}

		dataset.Sequences[i] = sequence
		dataset.Qualities[i] = qualities
	}

	dataset.Size = numSequences
	dataset.Metadata["sequence_length"] = sequenceLength
	dataset.Metadata["error_rate"] = errorRate
	dataset.Metadata["generated"] = true

	return nil
}

// GetReferenceGenome returns a cached reference genome
func (dm *DataManager) GetReferenceGenome(assembly string) (*reference.ReferenceGenome, error) {
	dm.mu.RLock()
	defer dm.mu.RUnlock()

	if genome, exists := dm.referenceCache[assembly]; exists {
		return genome, nil
	}

	return nil, fmt.Errorf("reference genome %s not loaded", assembly)
}

// GetPopulationPanel returns a cached population panel
func (dm *DataManager) GetPopulationPanel(populationName string) (*population.PopulationPanel, error) {
	dm.mu.RLock()
	defer dm.mu.RUnlock()

	if panel, exists := dm.populationCache[populationName]; exists {
		return panel, nil
	}

	return nil, fmt.Errorf("population panel %s not loaded", populationName)
}

// GetValidationDataset returns a cached validation dataset
func (dm *DataManager) GetValidationDataset(platform string) (*ValidationDataset, error) {
	dm.mu.RLock()
	defer dm.mu.RUnlock()

	if dataset, exists := dm.validationCache[platform]; exists {
		return dataset, nil
	}

	return nil, fmt.Errorf("validation dataset %s not loaded", platform)
}

// NewModelTrainer creates a new model trainer
func NewModelTrainer(config *TrainingConfig) *ModelTrainer {
	return &ModelTrainer{
		config:         config,
		hmmTrainer:     NewHMMTrainer(config),
		qualityTrainer: NewQualityTrainer(config),
		popTrainer:     NewPopulationTrainer(config),
		checkpointer:   NewModelCheckpointer(config.OutputPath),
	}
}

// TrainHMMModel trains an HMM k-mer model
func (mt *ModelTrainer) TrainHMMModel(k int, population string) error {
	return mt.hmmTrainer.TrainHMMModel(k, population)
}

// TrainQualityModel trains a quality model for a platform
func (mt *ModelTrainer) TrainQualityModel(platform string) error {
	return mt.qualityTrainer.TrainQualityModel(platform)
}

// TrainPopulationModel trains a population genetics model
func (mt *ModelTrainer) TrainPopulationModel(population string) error {
	return mt.popTrainer.TrainPopulationModel(population)
}

// NewValidationManager creates a new validation manager
func NewValidationManager(config *TrainingConfig) *ValidationManager {
	return &ValidationManager{
		config:               config,
		statisticalValidator: validation.NewStatisticalValidator(),
		expertTesting:        validation.NewExpertTestingFramework(),
		performanceTester:    NewPerformanceTester(config),
		results:              &ValidationResults{},
	}
}

// RunStatisticalValidation runs statistical validation tests.
//
// Without trained models (see SetSequenceGenerators) it stores fixed
// placeholder divergences for illumina, pacbio and nanopore and returns.
//
// With models, it generates 50 sequences of 200 bases per model, counts
// 4-mers, and computes the Jensen-Shannon divergence against a uniform
// distribution over the observed 4-mers. That measures k-mer uniformity; it
// is not a comparison with the training data or a held-out reference.
func (vm *ValidationManager) RunStatisticalValidation() error {
	log.Printf("Running statistical validation tests...")

	if len(vm.seqGenerators) == 0 {
		// No trained models available yet; store sentinel values.
		vm.results.JSdivergenceResults = map[string]float64{
			"illumina": 0.005,
			"pacbio":   0.008,
			"nanopore": 0.012,
		}
		return nil
	}

	const sampleSeqs = 50
	const seqLen = 200
	jsResults := make(map[string]float64, len(vm.seqGenerators))

	for modelKey, gen := range vm.seqGenerators {
		// Build k-mer frequency map from generated sequences.
		generatedFreqs := make(map[string]float64)
		for i := 0; i < sampleSeqs; i++ {
			result, err := gen.GenerateSequence(seqLen)
			if err != nil {
				return fmt.Errorf("statistical validation: generate seq for %s: %w", modelKey, err)
			}
			k := len(result.Sequence) / seqLen // approximate k from sequence context
			if k < 2 {
				k = 4
			}
			for j := 0; j <= len(result.Sequence)-k; j++ {
				kmer := string(result.Sequence[j : j+k])
				generatedFreqs[kmer]++
			}
		}
		// Normalise to a probability distribution.
		total := 0.0
		for _, v := range generatedFreqs {
			total += v
		}
		if total > 0 {
			for kmer := range generatedFreqs {
				generatedFreqs[kmer] /= total
			}
		}

		// Reference: uniform over the observed k-mers (not the training
		// distribution).
		referenceFreqs := make(map[string]float64, len(generatedFreqs))
		for kmer := range generatedFreqs {
			referenceFreqs[kmer] = 1.0 / float64(len(generatedFreqs))
		}
		result, err := vm.statisticalValidator.ValidateDistributions(generatedFreqs, referenceFreqs, "kmer_uniformity_"+modelKey)
		if err != nil {
			return fmt.Errorf("statistical validation: JS divergence for %s: %w", modelKey, err)
		}
		jsResults[modelKey] = result.JSdivergence
	}

	vm.results.JSdivergenceResults = jsResults
	return nil
}

// RunExpertTesting runs expert testing using the ExpertTestingFramework.
// When the framework has not been populated with test cases and registered experts,
// expert testing is skipped and ExpertTestResults remains nil — this is the honest
// result when no test data is available.  Populate the framework via
// expertTesting.GenerateTestCases / RegisterExpert / CollectExpertFeedback before
// calling RunExpertTesting to obtain real accuracy measurements.
func (vm *ValidationManager) RunExpertTesting() error {
	log.Printf("Running expert testing...")

	if err := vm.expertTesting.ValidateFramework(); err != nil {
		// Framework not configured (no test cases or experts) — skip gracefully.
		log.Printf("Expert testing skipped: framework not configured: %v", err)
		return nil
	}

	if err := vm.expertTesting.AggregateResults(); err != nil {
		return fmt.Errorf("expert testing: aggregate results: %w", err)
	}

	stats := vm.expertTesting.GetFrameworkStats()
	authenticityScore := vm.expertTesting.GetAuthenticityScore()
	overallAccuracy := 1.0 - authenticityScore/100.0

	numExperts, _ := stats["registered_experts"].(int)
	numCases, _ := stats["total_test_cases"].(int)
	numConsensus, _ := stats["consensus_results"].(int)

	vm.results.ExpertTestResults = &validation.ExpertTestingResults{
		OverallAccuracy:      overallAccuracy,
		AuthenticityScore:    authenticityScore,
		ParticipatingExperts: numExperts,
		TestCasesEvaluated:   numCases,
	}

	log.Printf("Expert testing complete: authenticity=%.1f%% (%d experts, %d cases, %d consensus)",
		authenticityScore, numExperts, numCases, numConsensus)

	return nil
}

// RunDownstreamTesting runs downstream compatibility testing.
//
// Stub: no tool is run. It always records a 0.95 success rate with BWA,
// GATK, BCFtools and SAMtools listed as compatible, so finalValidation's
// downstream check does not test anything yet.
func (vm *ValidationManager) RunDownstreamTesting() error {
	log.Printf("Running downstream compatibility testing...")

	vm.results.DownstreamResults = &DownstreamResults{
		SuccessRate:     0.95,
		CompatibleTools: []string{"BWA", "GATK", "BCFtools", "SAMtools"},
		FailedTools:     []string{},
		ErrorAnalysis:   map[string]string{},
	}

	return nil
}

// TestGenerationSpeed benchmarks sequence generation throughput using trained
// HMM models when available (100 sequences of 150 bases from the first
// model), otherwise it stores fixed placeholder values.
func (vm *ValidationManager) TestGenerationSpeed() error {
	log.Printf("Testing generation speed...")

	if vm.results.PerformanceResults == nil {
		vm.results.PerformanceResults = &PerformanceResults{}
	}

	// Pick the first available trained model for benchmarking.
	var gen authenticity.SequenceGenerator
	for _, g := range vm.seqGenerators {
		gen = g
		break
	}
	if gen == nil {
		// No models; keep existing placeholder values.
		vm.results.PerformanceResults.GenerationSpeed = 50000.0
		vm.results.PerformanceResults.Throughput = 45000.0
		vm.results.PerformanceResults.Latency = time.Millisecond * 10
		return nil
	}

	const bench = 100
	const readLen = 150
	start := time.Now()
	for i := 0; i < bench; i++ {
		if _, err := gen.GenerateSequence(readLen); err != nil {
			return fmt.Errorf("TestGenerationSpeed: generate: %w", err)
		}
	}
	elapsed := time.Since(start)
	perSeq := elapsed / bench
	basesPerSec := float64(readLen) / perSeq.Seconds()

	vm.results.PerformanceResults.GenerationSpeed = basesPerSec
	vm.results.PerformanceResults.Throughput = basesPerSec
	vm.results.PerformanceResults.Latency = perSeq
	return nil
}

// TestMemoryUsage measures heap allocation growth during a generation batch.
func (vm *ValidationManager) TestMemoryUsage() error {
	log.Printf("Testing memory usage...")

	if vm.results.PerformanceResults == nil {
		vm.results.PerformanceResults = &PerformanceResults{}
	}

	var gen authenticity.SequenceGenerator
	for _, g := range vm.seqGenerators {
		gen = g
		break
	}
	if gen == nil {
		vm.results.PerformanceResults.MemoryUsage = 1024 * 1024 * 150 // fallback
		return nil
	}

	var before, after runtime.MemStats
	runtime.GC()
	runtime.ReadMemStats(&before)
	for i := 0; i < 200; i++ {
		if _, err := gen.GenerateSequence(150); err != nil {
			return fmt.Errorf("TestMemoryUsage: generate: %w", err)
		}
	}
	runtime.ReadMemStats(&after)
	vm.results.PerformanceResults.MemoryUsage = int64(after.TotalAlloc - before.TotalAlloc)
	return nil
}

// TestScalability measures throughput efficiency at increasing goroutine counts.
// Scalability score = (throughput at 4 workers) / (4 × throughput at 1 worker) × 100.
func (vm *ValidationManager) TestScalability() error {
	log.Printf("Testing scalability...")

	if vm.results.PerformanceResults == nil {
		vm.results.PerformanceResults = &PerformanceResults{}
	}

	var gen authenticity.SequenceGenerator
	for _, g := range vm.seqGenerators {
		gen = g
		break
	}
	if gen == nil {
		vm.results.PerformanceResults.ScalabilityScore = 90.0 // fallback
		return nil
	}

	const iters = 40
	const readLen = 150

	throughput := func(workers int) float64 {
		start := time.Now()
		var wg sync.WaitGroup
		for w := 0; w < workers; w++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				for i := 0; i < iters/workers; i++ {
					gen.GenerateSequence(readLen) //nolint:errcheck
				}
			}()
		}
		wg.Wait()
		return float64(iters) / time.Since(start).Seconds()
	}

	t1 := throughput(1)
	t4 := throughput(4)
	score := 100.0
	if t1 > 0 {
		score = math.Min(100.0, (t4/(4*t1))*100.0)
	}
	vm.results.PerformanceResults.ScalabilityScore = score
	return nil
}

// AggregateResults aggregates all validation results
func (vm *ValidationManager) AggregateResults() *ValidationResults {
	log.Printf("Aggregating validation results...")

	// Calculate overall score; default to 0 so that missing data is conservative
	overallScore := 0.0

	// Adjust based on JS divergence results
	if vm.results.JSdivergenceResults != nil {
		avgJS := 0.0
		count := 0
		for _, js := range vm.results.JSdivergenceResults {
			avgJS += js
			count++
		}
		if count > 0 {
			avgJS /= float64(count)
			// Lower JS divergence = higher score
			overallScore = math.Max(0, 100.0-(avgJS*1000.0))
		}
	}

	vm.results.OverallScore = overallScore
	vm.results.PassesValidation = vm.results.JSdivergenceResults != nil && overallScore >= 95.0

	if !vm.results.PassesValidation {
		vm.results.Recommendations = append(vm.results.Recommendations,
			"Overall score below threshold - review model parameters")
	}

	return vm.results
}

// Execute runs the complete training and validation pipeline
func (p *TrainingPipeline) Execute() error {
	p.logger.Info("Starting Genomics Authenticity Training Pipeline")
	p.metrics.StartTime = time.Now()

	// Phase 1: Data Preparation
	if err := p.prepareData(); err != nil {
		return fmt.Errorf("data preparation failed: %w", err)
	}

	// Phase 2: Model Training
	if err := p.trainModels(); err != nil {
		return fmt.Errorf("model training failed: %w", err)
	}

	// Phase 3: Validation
	if err := p.validateModels(); err != nil {
		return fmt.Errorf("model validation failed: %w", err)
	}

	// Phase 4: Performance Testing
	if err := p.testPerformance(); err != nil {
		return fmt.Errorf("performance testing failed: %w", err)
	}

	// Phase 5: Final Validation
	if err := p.finalValidation(); err != nil {
		return fmt.Errorf("final validation failed: %w", err)
	}

	p.metrics.TrainingDuration = time.Since(p.metrics.StartTime)
	p.logger.Info("Training pipeline completed successfully")

	return nil
}

// prepareData loads and preprocesses all required data
func (p *TrainingPipeline) prepareData() error {
	p.logger.Info("Phase 1: Data Preparation")

	// Load reference genomes
	for _, assembly := range []string{"GRCh38", "GRCh39"} {
		if err := p.dataManager.LoadReferenceGenome(assembly); err != nil {
			return fmt.Errorf("failed to load reference genome %s: %w", assembly, err)
		}
	}

	// Load population panels
	for _, population := range p.config.Populations {
		if err := p.dataManager.LoadPopulationPanel(population); err != nil {
			return fmt.Errorf("failed to load population panel %s: %w", population, err)
		}
	}

	// Load validation datasets
	for _, platform := range p.config.Platforms {
		if err := p.dataManager.LoadValidationDataset(platform); err != nil {
			return fmt.Errorf("failed to load validation dataset %s: %w", platform, err)
		}
	}

	p.logger.Info("Data preparation completed")
	return nil
}

// trainModels trains all model components
func (p *TrainingPipeline) trainModels() error {
	p.logger.Info("Phase 2: Model Training")

	// Train HMM k-mer models
	for _, k := range p.config.KmerSizes {
		for _, population := range p.config.Populations {
			if err := p.modelTrainer.TrainHMMModel(k, population); err != nil {
				return fmt.Errorf("failed to train HMM model k=%d, pop=%s: %w", k, population, err)
			}
			p.metrics.ModelsTrained++
		}
	}

	// Train quality models
	for _, platform := range p.config.Platforms {
		if err := p.modelTrainer.TrainQualityModel(platform); err != nil {
			return fmt.Errorf("failed to train quality model for %s: %w", platform, err)
		}
		p.metrics.ModelsTrained++
	}

	// Train population genetics models
	for _, population := range p.config.Populations {
		if err := p.modelTrainer.TrainPopulationModel(population); err != nil {
			return fmt.Errorf("failed to train population model for %s: %w", population, err)
		}
		p.metrics.ModelsTrained++
	}

	// Wire trained HMM models into the validator so RunStatisticalValidation
	// and the performance benchmarks can use real GenerateSequence calls.
	hmmModels := p.modelTrainer.hmmTrainer.GetModels()
	gens := make(map[string]authenticity.SequenceGenerator, len(hmmModels))
	for key, m := range hmmModels {
		gens[key] = m
	}
	p.validator.SetSequenceGenerators(gens)

	p.logger.Info("Model training completed")
	return nil
}

// validateModels performs comprehensive model validation
func (p *TrainingPipeline) validateModels() error {
	p.logger.Info("Phase 3: Model Validation")

	// Statistical validation
	if err := p.validator.RunStatisticalValidation(); err != nil {
		return fmt.Errorf("statistical validation failed: %w", err)
	}

	// Expert testing (if enabled)
	if p.config.ExpertAccuracyThreshold > 0 {
		if err := p.validator.RunExpertTesting(); err != nil {
			return fmt.Errorf("expert testing failed: %w", err)
		}
	}

	// Downstream compatibility testing
	if err := p.validator.RunDownstreamTesting(); err != nil {
		return fmt.Errorf("downstream testing failed: %w", err)
	}

	p.logger.Info("Model validation completed")
	return nil
}

// testPerformance performs performance testing
func (p *TrainingPipeline) testPerformance() error {
	p.logger.Info("Phase 4: Performance Testing")

	// Generation speed testing
	if err := p.validator.TestGenerationSpeed(); err != nil {
		return fmt.Errorf("generation speed testing failed: %w", err)
	}

	// Memory usage testing
	if err := p.validator.TestMemoryUsage(); err != nil {
		return fmt.Errorf("memory usage testing failed: %w", err)
	}

	// Scalability testing
	if err := p.validator.TestScalability(); err != nil {
		return fmt.Errorf("scalability testing failed: %w", err)
	}

	p.logger.Info("Performance testing completed")
	return nil
}

// finalValidation performs final validation against all criteria
func (p *TrainingPipeline) finalValidation() error {
	p.logger.Info("Phase 5: Final Validation")

	// Aggregate all results
	results := p.validator.AggregateResults()

	// Check against thresholds
	if results.OverallScore < 95.0 {
		return fmt.Errorf("overall authenticity score too low: %.2f", results.OverallScore)
	}

	// Check JS divergence
	for platform, jsDivergence := range results.JSdivergenceResults {
		if jsDivergence > p.config.JSdivergenceThreshold {
			return fmt.Errorf("JS divergence too high for %s: %.6f", platform, jsDivergence)
		}
	}

	// Check expert testing results
	if results.ExpertTestResults != nil {
		if results.ExpertTestResults.OverallAccuracy > p.config.ExpertAccuracyThreshold {
			return fmt.Errorf("expert accuracy too high: %.3f", results.ExpertTestResults.OverallAccuracy)
		}
	}

	// Check downstream compatibility
	if results.DownstreamResults.SuccessRate < p.config.DownstreamSuccessThreshold {
		return fmt.Errorf("downstream success rate too low: %.3f", results.DownstreamResults.SuccessRate)
	}

	p.logger.Info("Final validation thresholds met (downstream compatibility is a stub and was not tested)")
	return nil
}

// Additional training components

// HMMTrainer handles HMM k-mer model training
type HMMTrainer struct {
	config        *TrainingConfig
	models        map[string]*authenticity.HMMKmerModel
	trainedModels map[string]bool
	mu            sync.RWMutex
	// genomeLoader is an optional override injected in tests (nil = stub/not implemented).
	genomeLoader func(population string) (*reference.ReferenceGenome, error)
}

// GetModels returns a snapshot of the trained HMMKmerModel map.
func (t *HMMTrainer) GetModels() map[string]*authenticity.HMMKmerModel {
	t.mu.RLock()
	defer t.mu.RUnlock()
	out := make(map[string]*authenticity.HMMKmerModel, len(t.models))
	for k, v := range t.models {
		out[k] = v
	}
	return out
}

// NewHMMTrainer creates a new HMM trainer
func NewHMMTrainer(config *TrainingConfig) *HMMTrainer {
	return &HMMTrainer{
		config:        config,
		models:        make(map[string]*authenticity.HMMKmerModel),
		trainedModels: make(map[string]bool),
	}
}

// TrainHMMModel trains an HMM k-mer model
func (t *HMMTrainer) TrainHMMModel(k int, population string) error {
	t.mu.Lock()
	defer t.mu.Unlock()

	modelKey := fmt.Sprintf("%dmer_%s", k, population)

	// Create model
	model := authenticity.NewHMMKmerModel(k)

	// Load reference genome for population
	refGenome, err := t.loadReferenceGenome(population)
	if err != nil {
		return err
	}

	// Train model
	if err := model.TrainFromReferenceGenome(refGenome); err != nil {
		return err
	}

	// Validate model
	if err := model.ValidateModel(); err != nil {
		return err
	}

	// Store model
	t.models[modelKey] = model
	t.trainedModels[modelKey] = true

	// Save model to disk
	if err := t.saveModel(model, modelKey); err != nil {
		return err
	}

	return nil
}

// loadReferenceGenome loads the reference genome for the given population tag.
// In production it reads FASTA files from TrainingConfig.ReferenceGenomePath.
// Tests override this via the genomeLoader DI hook (see trainer_di_test.go).
func (t *HMMTrainer) loadReferenceGenome(pop string) (*reference.ReferenceGenome, error) {
	if t.genomeLoader != nil {
		return t.genomeLoader(pop)
	}
	if t.config.ReferenceGenomePath == "" {
		return nil, fmt.Errorf("loadReferenceGenome: ReferenceGenomePath not configured")
	}
	// Resolve population-specific subdirectory if it exists (e.g. GRCh38/YRI/),
	// otherwise fall back to the base path so a shared reference works too.
	popDir := filepath.Join(t.config.ReferenceGenomePath, pop)
	if info, err := os.Stat(popDir); err == nil && info.IsDir() {
		return reference.NewReferenceGenome(popDir), nil
	}
	return reference.NewReferenceGenome(t.config.ReferenceGenomePath), nil
}

// saveModel saves a trained HMMKmerModel to disk under the configured output path.
func (t *HMMTrainer) saveModel(model *authenticity.HMMKmerModel, key string) error {
	if t.config.OutputPath == "" {
		return nil // no output path configured; skip persistence
	}
	if err := security.SecureMkdirAll(t.config.OutputPath, filepath.Dir(t.config.OutputPath)); err != nil {
		return fmt.Errorf("create output directory: %w", err)
	}
	path := filepath.Join(t.config.OutputPath, key+".hmm.json")
	return model.Save(path)
}

// QualityTrainer handles quality model training
type QualityTrainer struct {
	config *TrainingConfig
	models map[string]*quality.GCbiasQualityModel
	mu     sync.RWMutex
	// dataLoader is an optional override injected in tests (nil = stub/not implemented).
	dataLoader func(platform string) ([]quality.GCQualityDataPoint, error)
}

// NewQualityTrainer creates a new quality trainer
func NewQualityTrainer(config *TrainingConfig) *QualityTrainer {
	return &QualityTrainer{
		config: config,
		models: make(map[string]*quality.GCbiasQualityModel),
	}
}

// TrainQualityModel trains a quality model for a platform
func (t *QualityTrainer) TrainQualityModel(platform string) error {
	t.mu.Lock()
	defer t.mu.Unlock()

	// Create model
	model := quality.NewGCbiasQualityModel(platform)

	// Load training data
	trainingData, err := t.loadQualityTrainingData(platform)
	if err != nil {
		return err
	}

	// Train model
	if err := model.TrainFromRealData(trainingData); err != nil {
		return err
	}

	// Validate model
	if err := model.ValidateModel(); err != nil {
		return err
	}

	// Store model
	t.models[platform] = model

	return nil
}

// loadQualityTrainingData loads GC-quality data points for a sequencing platform.
// It reads <QualityDataPath>/<platform>.json, which must contain a JSON array of
// GCQualityDataPoint objects. Tests override this via the dataLoader DI hook.
func (t *QualityTrainer) loadQualityTrainingData(platform string) ([]quality.GCQualityDataPoint, error) {
	if t.dataLoader != nil {
		return t.dataLoader(platform)
	}
	if t.config.QualityDataPath == "" {
		return nil, fmt.Errorf("loadQualityTrainingData: QualityDataPath not configured")
	}
	path := filepath.Join(t.config.QualityDataPath, platform+".json")
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("loadQualityTrainingData: read %s: %w", path, err)
	}
	var points []quality.GCQualityDataPoint
	if err := json.Unmarshal(raw, &points); err != nil {
		return nil, fmt.Errorf("loadQualityTrainingData: parse %s: %w", path, err)
	}
	if len(points) == 0 {
		return nil, fmt.Errorf("loadQualityTrainingData: %s contains no data points", path)
	}
	return points, nil
}

// PopulationTrainer handles population genetics model training
type PopulationTrainer struct {
	config *TrainingConfig
	models map[string]*population.CoalescentSimulator
	mu     sync.RWMutex
	// paramsLoader is an optional override injected in tests (nil = stub/not implemented).
	paramsLoader func(populationName string) (*PopulationParameters, error)
}

// NewPopulationTrainer creates a new population trainer
func NewPopulationTrainer(config *TrainingConfig) *PopulationTrainer {
	return &PopulationTrainer{
		config: config,
		models: make(map[string]*population.CoalescentSimulator),
	}
}

// TrainPopulationModel trains a population genetics model
func (t *PopulationTrainer) TrainPopulationModel(populationName string) error {
	t.mu.Lock()
	defer t.mu.Unlock()

	// Get population parameters
	params, err := t.getPopulationParameters(populationName)
	if err != nil {
		return err
	}

	// Create simulator
	simulator := population.NewCoalescentSimulator(
		params.EffectivePopSize,
		params.MutationRate,
		params.RecombinationRate,
	)

	// Add demographic events
	for _, event := range params.DemographicEvents {
		simulator.AddDemographicEvent(event)
	}

	// Validate parameters
	if err := simulator.ValidateParameters(); err != nil {
		return err
	}

	// Store model
	t.models[populationName] = simulator

	return nil
}

// getPopulationParameters loads genetics parameters for the named population.
// It reads <PopulationDataPath>/<populationName>.json, which must contain a
// JSON object matching PopulationParameters. Tests override via paramsLoader.
func (t *PopulationTrainer) getPopulationParameters(populationName string) (*PopulationParameters, error) {
	if t.paramsLoader != nil {
		return t.paramsLoader(populationName)
	}
	if t.config.PopulationDataPath == "" {
		return nil, fmt.Errorf("getPopulationParameters: PopulationDataPath not configured")
	}
	path := filepath.Join(t.config.PopulationDataPath, populationName+".json")
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("getPopulationParameters: read %s: %w", path, err)
	}
	var params PopulationParameters
	if err := json.Unmarshal(raw, &params); err != nil {
		return nil, fmt.Errorf("getPopulationParameters: parse %s: %w", path, err)
	}
	return &params, nil
}

// PopulationParameters represents population-specific parameters
type PopulationParameters struct {
	EffectivePopSize  int                           `json:"effective_pop_size"`
	MutationRate      float64                       `json:"mutation_rate"`
	RecombinationRate float64                       `json:"recombination_rate"`
	DemographicEvents []population.DemographicEvent `json:"demographic_events,omitempty"`
}

// PipelineLogger handles logging for the training pipeline
type PipelineLogger struct {
	logFile *os.File
	mu      sync.Mutex
}

// NewPipelineLogger creates a new pipeline logger
func NewPipelineLogger(outputPath string) *PipelineLogger {
	logPath := filepath.Join(outputPath, "training.log")
	logFile, err := os.OpenFile(logPath, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0600) // #nosec G304
	if err != nil {
		log.Printf("Failed to open log file: %v", err)
		return &PipelineLogger{}
	}

	return &PipelineLogger{
		logFile: logFile,
	}
}

// Info logs an info message
func (l *PipelineLogger) Info(message string) {
	l.mu.Lock()
	defer l.mu.Unlock()

	timestamp := time.Now().Format("2006-01-02 15:04:05")
	logMessage := fmt.Sprintf("[INFO] %s: %s\n", timestamp, message)

	if l.logFile != nil {
		l.logFile.WriteString(logMessage)
	}
	log.Print(logMessage)
}

// Error logs an error message
func (l *PipelineLogger) Error(message string, err error) {
	l.mu.Lock()
	defer l.mu.Unlock()

	timestamp := time.Now().Format("2006-01-02 15:04:05")
	logMessage := fmt.Sprintf("[ERROR] %s: %s: %v\n", timestamp, message, err)

	if l.logFile != nil {
		l.logFile.WriteString(logMessage)
	}
	log.Print(logMessage)
}

// NewTrainingMetrics creates new training metrics
func NewTrainingMetrics() *TrainingMetrics {
	return &TrainingMetrics{
		ValidationScores:   make(map[string]float64),
		PerformanceMetrics: make(map[string]float64),
		ErrorLog:           make([]TrainingError, 0),
	}
}

// Additional performance testing structures

// PerformanceTester handles performance testing
type PerformanceTester struct {
	config  *TrainingConfig
	results *PerformanceResults
}

// PerformanceResults stores performance test results
type PerformanceResults struct {
	GenerationSpeed  float64 // bases/second
	MemoryUsage      int64   // bytes
	ScalabilityScore float64 // 0-100
	Latency          time.Duration
	Throughput       float64
}

// DownstreamResults stores downstream compatibility test results
type DownstreamResults struct {
	SuccessRate     float64
	CompatibleTools []string
	FailedTools     []string
	ErrorAnalysis   map[string]string
}

// NewPerformanceTester creates a new performance tester
func NewPerformanceTester(config *TrainingConfig) *PerformanceTester {
	return &PerformanceTester{
		config:  config,
		results: &PerformanceResults{},
	}
}

// ModelCheckpointer handles model checkpointing
type ModelCheckpointer struct {
	outputPath string
	mu         sync.RWMutex
}

// NewModelCheckpointer creates a new model checkpointer
func NewModelCheckpointer(outputPath string) *ModelCheckpointer {
	return &ModelCheckpointer{
		outputPath: outputPath,
	}
}

// saveModel persists a trained model to the output directory.
// The key encodes "population/platform" and becomes the file stem.
func (c *ModelCheckpointer) saveModel(model interface{}, key string) error {
	c.mu.Lock()
	defer c.mu.Unlock()

	if c.outputPath == "" {
		return nil
	}
	dir := filepath.Join(c.outputPath, filepath.Dir(key))
	if err := security.SecureMkdirAll(dir, c.outputPath); err != nil {
		return fmt.Errorf("checkpointer: create directory: %w", err)
	}
	path := filepath.Join(c.outputPath, key+".hmm.json")

	switch m := model.(type) {
	case *authenticity.AdaptiveHMMModel:
		return m.Save(path)
	case *authenticity.HMMKmerModel:
		return m.Save(path)
	default:
		return fmt.Errorf("checkpointer: unsupported model type: %T", model)
	}
}

// loadModel loads a previously saved model from the output directory.
// It inspects the ModelEnvelope to determine the concrete type and returns
// either a *authenticity.AdaptiveHMMModel or *authenticity.HMMKmerModel.
func (c *ModelCheckpointer) loadModel(key string) (interface{}, error) {
	c.mu.RLock()
	defer c.mu.RUnlock()

	path := filepath.Join(c.outputPath, key+".hmm.json")

	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("checkpointer: read %s: %w", path, err)
	}

	// Peek at the model_type field without full deserialization.
	var env authenticity.ModelEnvelope
	if err := json.Unmarshal(raw, &env); err != nil {
		return nil, fmt.Errorf("checkpointer: parse envelope: %w", err)
	}

	switch env.ModelType {
	case "AdaptiveHMMModel":
		return authenticity.LoadAdaptiveHMMModel(path)
	case "HMMKmerModel":
		return authenticity.LoadHMMKmerModel(path)
	default:
		return nil, fmt.Errorf("checkpointer: unknown model type %q", env.ModelType)
	}
}
