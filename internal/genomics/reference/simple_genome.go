package reference

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// SimpleReferenceGenome provides a simple implementation of reference genome access
type SimpleReferenceGenome struct {
	baseDir string
}

// NewSimpleReferenceGenome creates a new SimpleReferenceGenome instance
func NewSimpleReferenceGenome(baseDir string) *SimpleReferenceGenome {
	return &SimpleReferenceGenome{
		baseDir: baseDir,
	}
}

// GetChromosomeNames returns the names of available chromosomes
func (g *SimpleReferenceGenome) GetChromosomeNames() []string {
	files, err := os.ReadDir(g.baseDir)
	if err != nil {
		return nil
	}

	var chromosomes []string
	for _, file := range files {
		if strings.HasSuffix(file.Name(), ".fasta") && strings.HasPrefix(file.Name(), "chr") {
			// Extract chromosome name from filename (e.g., "chr21.fasta" -> "chr21")
			chrName := strings.TrimSuffix(file.Name(), ".fasta")
			chromosomes = append(chromosomes, chrName)
		}
	}

	return chromosomes
}

// LoadChromosome reads <baseDir>/<name>.fasta, or <baseDir>/GRCh38.fasta when
// that file does not exist (see readFastaFile for how records are matched).
func (g *SimpleReferenceGenome) LoadChromosome(name string) (authenticity.ChromosomeInterface, error) {
	fastaFile := filepath.Join(g.baseDir, name+".fasta")

	// Check if file exists
	if _, err := os.Stat(fastaFile); os.IsNotExist(err) {
		// Try the main GRCh38.fasta file
		fastaFile = filepath.Join(g.baseDir, "GRCh38.fasta")
		if _, err := os.Stat(fastaFile); os.IsNotExist(err) {
			return nil, fmt.Errorf("chromosome file not found: %s", name)
		}
	}

	// Read FASTA file
	sequence, err := readFastaFile(fastaFile, name)
	if err != nil {
		return nil, fmt.Errorf("failed to read chromosome %s: %v", name, err)
	}

	return &SimpleChromosome{
		name:     name,
		sequence: sequence,
		length:   len(sequence),
	}, nil
}

// SimpleChromosome represents a simple chromosome implementation
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

// readFastaFile reads the whole FASTA file into memory and returns the
// upper-cased sequence of every record whose header contains
// targetChromosome as a substring, so "chr1" also collects chr10–chr19 when
// they are in the same file.
func readFastaFile(filename, targetChromosome string) ([]byte, error) {
	data, err := os.ReadFile(filename) // #nosec G304
	if err != nil {
		return nil, err
	}

	lines := strings.Split(string(data), "\n")
	var sequence []byte
	inTargetChromosome := false

	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}

		if strings.HasPrefix(line, ">") {
			// Header line
			if strings.Contains(line, targetChromosome) {
				inTargetChromosome = true
			} else {
				inTargetChromosome = false
			}
		} else if inTargetChromosome {
			// Sequence line
			sequence = append(sequence, []byte(strings.ToUpper(line))...)
		}
	}

	if len(sequence) == 0 {
		return nil, fmt.Errorf("chromosome %s not found in FASTA file", targetChromosome)
	}

	return sequence, nil
}

// ChromosomeInterface defines the interface for chromosome access
type ChromosomeInterface interface {
	GetSequence() []byte
	GetName() string
	GetLength() int
}

// ReferenceGenomeInterface defines the interface for reference genome access
type ReferenceGenomeInterface interface {
	GetChromosomeNames() []string
	LoadChromosome(name string) (ChromosomeInterface, error)
}

// ReferenceGenome is an alias for SimpleReferenceGenome to maintain compatibility
type ReferenceGenome = SimpleReferenceGenome

// NewReferenceGenome creates a new reference genome instance
func NewReferenceGenome(baseDir string) *ReferenceGenome {
	return NewSimpleReferenceGenome(baseDir)
}
