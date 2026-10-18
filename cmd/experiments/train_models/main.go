// Command train_models trains HMMKmerModel models for k = 3..6 on chr21 and
// chr22 FASTA files from the directory given as its only argument (falling
// back to GRCh38.fasta there), checks the GC content and alphabet of a few
// generated sequences, and writes per-model summaries and a training report
// to ./models/trained. The trained models themselves are not saved.
//
// It is a one-off experiment, not part of the supported phantom-science CLI;
// use "phantom-science train" instead. The JS-divergence and expert-accuracy
// thresholds in its configuration are recorded in the report but not
// checked.
package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/training"
)

// ReferenceGenomeLoader loads reference genome data from FASTA files
type ReferenceGenomeLoader struct {
	referenceDir string
	chromosomes  map[string][]byte
}

// NewReferenceGenomeLoader creates a new loader for reference genome data
func NewReferenceGenomeLoader(referenceDir string) *ReferenceGenomeLoader {
	return &ReferenceGenomeLoader{
		referenceDir: referenceDir,
		chromosomes:  make(map[string][]byte),
	}
}

// LoadChromosome loads a chromosome from FASTA file
func (r *ReferenceGenomeLoader) LoadChromosome(chrName string) ([]byte, error) {
	if sequence, exists := r.chromosomes[chrName]; exists {
		return sequence, nil
	}

	filePath := filepath.Join(r.referenceDir, chrName+".fasta")
	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		// Try alternative naming
		filePath = filepath.Join(r.referenceDir, "GRCh38.fasta")
	}

	file, err := os.Open(filePath) // #nosec G304
	if err != nil {
		return nil, fmt.Errorf("failed to open file %s: %w", filePath, err)
	}
	defer file.Close()

	var sequence []byte
	scanner := bufio.NewScanner(file)
	inTargetChromosome := false
	targetHeader := ">" + chrName + " "

	for scanner.Scan() {
		line := scanner.Text()

		if strings.HasPrefix(line, ">") {
			if strings.HasPrefix(line, targetHeader) {
				inTargetChromosome = true
				log.Printf("Found chromosome %s", chrName)
				continue
			} else if inTargetChromosome {
				// Reached next chromosome
				break
			}
		} else if inTargetChromosome {
			// Process sequence line
			for _, base := range line {
				if base == 'A' || base == 'T' || base == 'G' || base == 'C' {
					sequence = append(sequence, byte(base))
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	if len(sequence) == 0 {
		return nil, fmt.Errorf("no sequence found for chromosome %s", chrName)
	}

	log.Printf("Loaded chromosome %s with %d bases", chrName, len(sequence))
	r.chromosomes[chrName] = sequence
	return sequence, nil
}

// MockReferenceGenome implements a reference genome interface for training
type MockReferenceGenome struct {
	loader      *ReferenceGenomeLoader
	chromosomes []string
}

// NewMockReferenceGenome creates a new mock reference genome
func NewMockReferenceGenome(loader *ReferenceGenomeLoader) *MockReferenceGenome {
	return &MockReferenceGenome{
		loader:      loader,
		chromosomes: []string{"chr21", "chr22"}, // Start with smaller chromosomes
	}
}

// GetChromosomeNames returns available chromosome names
func (m *MockReferenceGenome) GetChromosomeNames() []string {
	return m.chromosomes
}

// LoadChromosome loads a chromosome for training
func (m *MockReferenceGenome) LoadChromosome(name string) (authenticity.ChromosomeInterface, error) {
	sequence, err := m.loader.LoadChromosome(name)
	if err != nil {
		return nil, err
	}

	return &MockChromosome{
		Name:     name,
		Sequence: sequence,
		Length:   len(sequence),
	}, nil
}

// MockChromosome represents a chromosome for training
type MockChromosome struct {
	Name     string
	Sequence []byte
	Length   int
}

// GetSequence returns the chromosome sequence
func (c *MockChromosome) GetSequence() []byte {
	return c.Sequence
}

// GetName returns the chromosome name
func (c *MockChromosome) GetName() string {
	return c.Name
}

// GetLength returns the chromosome length
func (c *MockChromosome) GetLength() int {
	return c.Length
}

// TrainingSession manages the model training process
type TrainingSession struct {
	config          *training.TrainingConfig
	referenceGenome *MockReferenceGenome
	trainedModels   map[string]*authenticity.HMMKmerModel
	logger          *log.Logger
	startTime       time.Time
}

// NewTrainingSession creates a new training session
func NewTrainingSession(referenceDir string) *TrainingSession {
	loader := NewReferenceGenomeLoader(referenceDir)
	refGenome := NewMockReferenceGenome(loader)

	config := &training.TrainingConfig{
		KmerSizes:               []int{3, 4, 5, 6},
		Populations:             []string{"EUR", "AFR", "EAS"},
		JSdivergenceThreshold:   0.01,
		ExpertAccuracyThreshold: 0.55,
		MaxMemoryGB:             8,
		NumWorkers:              4,
		OutputPath:              "./models/trained",
	}

	return &TrainingSession{
		config:          config,
		referenceGenome: refGenome,
		trainedModels:   make(map[string]*authenticity.HMMKmerModel),
		logger:          log.New(os.Stdout, "[TRAINING] ", log.LstdFlags),
		startTime:       time.Now(),
	}
}

// TrainHMMModels trains HMM k-mer models on real genomic data
func (s *TrainingSession) TrainHMMModels() error {
	s.logger.Println("Starting HMM k-mer model training with real genomic data")

	// Create output directory
	if err := os.MkdirAll(s.config.OutputPath, 0750); err != nil {
		return fmt.Errorf("failed to create output directory: %w", err)
	}

	// Train models for each k-mer size
	for _, k := range s.config.KmerSizes {
		modelKey := fmt.Sprintf("hmm_kmer_%d", k)
		s.logger.Printf("Training %d-mer HMM model", k)

		// Create model
		model := authenticity.NewHMMKmerModel(k)

		// Train model on reference genome
		if err := model.TrainFromReferenceGenome(s.referenceGenome); err != nil {
			return fmt.Errorf("failed to train %d-mer model: %w", k, err)
		}

		// Validate model
		if err := model.ValidateModel(); err != nil {
			return fmt.Errorf("model validation failed for %d-mer: %w", k, err)
		}

		// Test generation
		s.logger.Printf("Testing %d-mer model generation", k)
		testResult, err := model.GenerateSequence(1000)
		if err != nil {
			return fmt.Errorf("generation test failed for %d-mer: %w", k, err)
		}

		s.logger.Printf("✓ %d-mer model trained successfully", k)
		s.logger.Printf("  - Generated %d bases", len(testResult.Sequence))
		s.logger.Printf("  - K-mer frequencies: %d unique k-mers", len(testResult.KmerFrequencies))
		s.logger.Printf("  - Generation time: %v", testResult.GenerationTime)

		// Store model
		s.trainedModels[modelKey] = model

		// Save model summary
		if err := s.saveModelSummary(modelKey, model, testResult); err != nil {
			s.logger.Printf("Warning: failed to save model summary: %v", err)
		}
	}

	s.logger.Printf("All HMM models trained successfully in %v", time.Since(s.startTime))
	return nil
}

// saveModelSummary saves a summary of the trained model
func (s *TrainingSession) saveModelSummary(modelKey string, model *authenticity.HMMKmerModel, testResult *authenticity.GenerationResult) error {
	summaryPath := filepath.Join(s.config.OutputPath, modelKey+"_summary.txt")
	file, err := os.Create(summaryPath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	fmt.Fprintf(file, "Model: %s\n", modelKey)
	fmt.Fprintf(file, "K-mer size: %d\n", model.K())
	fmt.Fprintf(file, "Trained: %v\n", model.IsTrained())
	fmt.Fprintf(file, "Total k-mers processed: %d\n", model.TotalKmers())
	fmt.Fprintf(file, "State space size: %d\n", len(model.StateSpace()))
	fmt.Fprintf(file, "\nTest generation results:\n")
	fmt.Fprintf(file, "- Sequence length: %d\n", len(testResult.Sequence))
	fmt.Fprintf(file, "- Unique k-mers: %d\n", len(testResult.KmerFrequencies))
	fmt.Fprintf(file, "- Generation time: %v\n", testResult.GenerationTime)
	fmt.Fprintf(file, "\nTop 10 k-mers:\n")

	// Show top k-mers
	count := 0
	for kmer, freq := range testResult.KmerFrequencies {
		if count >= 10 {
			break
		}
		fmt.Fprintf(file, "  %s: %d\n", kmer, freq)
		count++
	}

	return nil
}

// ValidateModels performs validation on trained models
func (s *TrainingSession) ValidateModels() error {
	s.logger.Println("Starting model validation")

	for modelKey, model := range s.trainedModels {
		s.logger.Printf("Validating %s", modelKey)

		// Generate multiple sequences for validation
		sequences := make([][]byte, 10)
		for i := 0; i < 10; i++ {
			result, err := model.GenerateSequence(5000)
			if err != nil {
				return fmt.Errorf("validation generation failed for %s: %w", modelKey, err)
			}
			sequences[i] = result.Sequence
		}

		// Basic validation checks
		if err := s.validateSequenceQuality(sequences); err != nil {
			return fmt.Errorf("sequence quality validation failed for %s: %w", modelKey, err)
		}

		s.logger.Printf("✓ %s validation passed", modelKey)
	}

	return nil
}

// validateSequenceQuality performs basic quality checks on generated sequences
func (s *TrainingSession) validateSequenceQuality(sequences [][]byte) error {
	for i, seq := range sequences {
		// Check length
		if len(seq) == 0 {
			return fmt.Errorf("sequence %d is empty", i)
		}

		// Check base composition
		baseCounts := make(map[byte]int)
		for _, base := range seq {
			baseCounts[base]++
		}

		// Ensure all bases are valid
		for base := range baseCounts {
			if base != 'A' && base != 'T' && base != 'G' && base != 'C' {
				return fmt.Errorf("sequence %d contains invalid base: %c", i, base)
			}
		}

		// Check for reasonable GC content (20-80%)
		gcCount := baseCounts['G'] + baseCounts['C']
		gcContent := float64(gcCount) / float64(len(seq))
		if gcContent < 0.2 || gcContent > 0.8 {
			return fmt.Errorf("sequence %d has unrealistic GC content: %.2f", i, gcContent)
		}
	}

	return nil
}

// RunTraining executes the complete training pipeline
func (s *TrainingSession) RunTraining() error {
	s.logger.Println("=== Phantom Science HMM Model Training ===")

	// Phase 1: Train HMM models
	if err := s.TrainHMMModels(); err != nil {
		return fmt.Errorf("HMM training failed: %w", err)
	}

	// Phase 2: Validate models
	if err := s.ValidateModels(); err != nil {
		return fmt.Errorf("model validation failed: %w", err)
	}

	// Phase 3: Generate training report
	if err := s.generateTrainingReport(); err != nil {
		return fmt.Errorf("report generation failed: %w", err)
	}

	s.logger.Printf("Training pipeline completed successfully in %v", time.Since(s.startTime))
	return nil
}

// generateTrainingReport generates a comprehensive training report
func (s *TrainingSession) generateTrainingReport() error {
	reportPath := filepath.Join(s.config.OutputPath, "training_report.txt")
	file, err := os.Create(reportPath) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	fmt.Fprintf(file, "Phantom Science HMM Model Training Report\n")
	fmt.Fprintf(file, "=========================================\n\n")
	fmt.Fprintf(file, "Training started: %v\n", s.startTime)
	fmt.Fprintf(file, "Training duration: %v\n", time.Since(s.startTime))
	fmt.Fprintf(file, "Models trained: %d\n\n", len(s.trainedModels))

	fmt.Fprintf(file, "Configuration:\n")
	fmt.Fprintf(file, "- K-mer sizes: %v\n", s.config.KmerSizes)
	fmt.Fprintf(file, "- JS divergence threshold: %.6f\n", s.config.JSdivergenceThreshold)
	fmt.Fprintf(file, "- Expert accuracy threshold: %.3f\n", s.config.ExpertAccuracyThreshold)
	fmt.Fprintf(file, "- Output path: %s\n\n", s.config.OutputPath)

	fmt.Fprintf(file, "Trained Models:\n")
	for modelKey, model := range s.trainedModels {
		fmt.Fprintf(file, "- %s: k=%d, states=%d, trained=%v\n",
			modelKey, model.K(), len(model.StateSpace()), model.IsTrained())
	}

	fmt.Fprintf(file, "\nNext Steps:\n")
	fmt.Fprintf(file, "1. Validate against real FASTQ datasets\n")
	fmt.Fprintf(file, "2. Measure JS divergence against reference data\n")
	fmt.Fprintf(file, "3. Conduct expert blind testing\n")
	fmt.Fprintf(file, "4. Performance benchmarking\n")
	fmt.Fprintf(file, "5. Integration with existing systems\n")

	return nil
}

func main() {
	if len(os.Args) < 2 {
		log.Fatal("Usage: go run main.go <reference_genome_directory>")
	}

	referenceDir := os.Args[1]
	if _, err := os.Stat(referenceDir); os.IsNotExist(err) { // #nosec G703
		log.Fatalf("Reference genome directory does not exist: %s", referenceDir)
	}

	// Create training session
	session := NewTrainingSession(referenceDir)

	// Run training
	if err := session.RunTraining(); err != nil {
		log.Fatalf("Training failed: %v", err)
	}

	log.Println("Training completed successfully!")
}
