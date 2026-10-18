package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// SimpleReferenceGenome implements a simple reference genome
type SimpleReferenceGenome struct {
	chromosomes map[string]*SimpleChromosome
	refDir      string
}

// SimpleChromosome implements a simple chromosome
type SimpleChromosome struct {
	name     string
	sequence []byte
	length   int
}

// GetSequence returns the chromosome sequence
func (c *SimpleChromosome) GetSequence() []byte {
	return c.sequence
}

// GetName returns the chromosome name
func (c *SimpleChromosome) GetName() string {
	return c.name
}

// GetLength returns the chromosome length
func (c *SimpleChromosome) GetLength() int {
	return c.length
}

// NewSimpleReferenceGenome creates a new simple reference genome
func NewSimpleReferenceGenome(refDir string) *SimpleReferenceGenome {
	return &SimpleReferenceGenome{
		chromosomes: make(map[string]*SimpleChromosome),
		refDir:      refDir,
	}
}

// GetChromosomeNames returns available chromosome names
func (r *SimpleReferenceGenome) GetChromosomeNames() []string {
	return []string{"chr21", "chr22"}
}

// LoadChromosome loads a chromosome from the reference directory
func (r *SimpleReferenceGenome) LoadChromosome(name string) (authenticity.ChromosomeInterface, error) {
	if chr, exists := r.chromosomes[name]; exists {
		return chr, nil
	}

	// Load chromosome from FASTA file
	sequence, err := r.loadChromosomeFromFasta(name)
	if err != nil {
		return nil, err
	}

	chr := &SimpleChromosome{
		name:     name,
		sequence: sequence,
		length:   len(sequence),
	}

	r.chromosomes[name] = chr
	return chr, nil
}

// loadChromosomeFromFasta loads a chromosome sequence from FASTA file
func (r *SimpleReferenceGenome) loadChromosomeFromFasta(chrName string) ([]byte, error) {
	// Try to load from GRCh38.fasta (contains all chromosomes)
	fastaPath := r.refDir + "/GRCh38.fasta"
	
	file, err := os.Open(fastaPath) // #nosec G304
	if err != nil {
		return nil, fmt.Errorf("failed to open FASTA file: %w", err)
	}
	defer file.Close()

	var sequence []byte
	scanner := bufio.NewScanner(file)
	inTargetChromosome := false
	targetHeader := ">" + chrName + " "

	log.Printf("Loading chromosome %s from %s", chrName, fastaPath)

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
			// Process sequence line, only keep valid bases
			for _, base := range line {
				if base == 'A' || base == 'T' || base == 'G' || base == 'C' {
					sequence = append(sequence, byte(base))
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading FASTA file: %w", err)
	}

	if len(sequence) == 0 {
		return nil, fmt.Errorf("no sequence found for chromosome %s", chrName)
	}

	log.Printf("Loaded chromosome %s with %d bases", chrName, len(sequence))
	return sequence, nil
}

func main() {
	if len(os.Args) < 2 {
		log.Fatal("Usage: go run main.go <reference_genome_directory>")
	}

	referenceDir := os.Args[1]
	log.Printf("=== Phantom Science HMM Model Training ===")
	log.Printf("Reference genome directory: %s", referenceDir)

	// Create reference genome
	refGenome := NewSimpleReferenceGenome(referenceDir)

	// Create output directory
	os.MkdirAll("./models/trained", 0750)

	// Train models for different k-mer sizes
	kmerSizes := []int{3, 4, 5, 6}
	
	for _, k := range kmerSizes {
		log.Printf("Training %d-mer HMM model...", k)
		
		// Create model
		model := authenticity.NewHMMKmerModel(k)
		
		// Train model
		start := time.Now()
		err := model.TrainFromReferenceGenome(refGenome)
		if err != nil {
			log.Fatalf("Training failed for %d-mer model: %v", k, err)
		}
		trainingTime := time.Since(start)
		
		log.Printf("✓ %d-mer model trained in %v", k, trainingTime)
		
		// Validate model
		err = model.ValidateModel()
		if err != nil {
			log.Fatalf("Validation failed for %d-mer model: %v", k, err)
		}
		
		log.Printf("✓ %d-mer model validated", k)
		
		// Test generation
		log.Printf("Testing %d-mer model generation...", k)
		testResult, err := model.GenerateSequence(10000)
		if err != nil {
			log.Fatalf("Generation test failed for %d-mer model: %v", k, err)
		}
		
		log.Printf("✓ %d-mer model generated %d bases in %v", k, len(testResult.Sequence), testResult.GenerationTime)
		log.Printf("  - Unique k-mers: %d", len(testResult.KmerFrequencies))
		
		// Calculate basic statistics
		gcCount := 0
		for _, base := range testResult.Sequence {
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}
		gcContent := float64(gcCount) / float64(len(testResult.Sequence))
		log.Printf("  - GC content: %.1f%%", gcContent*100)
		
		// Save model summary
		summaryPath := fmt.Sprintf("./models/trained/hmm_%dmer_summary.txt", k)
		file, err := os.Create(summaryPath) // #nosec G304
		if err == nil {
			fmt.Fprintf(file, "HMM %d-mer Model Summary\n", k)
			fmt.Fprintf(file, "======================\n")
			fmt.Fprintf(file, "Training time: %v\n", trainingTime)
			fmt.Fprintf(file, "State space size: %d\n", len(model.StateSpace()))
			fmt.Fprintf(file, "Total k-mers processed: %d\n", model.TotalKmers())
			fmt.Fprintf(file, "Test generation: %d bases\n", len(testResult.Sequence))
			fmt.Fprintf(file, "GC content: %.1f%%\n", gcContent*100)
			fmt.Fprintf(file, "Unique k-mers: %d\n", len(testResult.KmerFrequencies))
			fmt.Fprintf(file, "Generation time: %v\n", testResult.GenerationTime)
			file.Close()
		}
		
		log.Printf("✓ %d-mer model summary saved to %s", k, summaryPath)
	}
	
	log.Printf("=== Training Complete ===")
	log.Printf("All HMM models trained successfully!")
	log.Printf("Next steps:")
	log.Printf("1. Run tests: go test ./internal/genomics/authenticity/...")
	log.Printf("2. Validate against real FASTQ datasets")
	log.Printf("3. Measure JS divergence")
	log.Printf("4. Conduct expert blind testing")
}