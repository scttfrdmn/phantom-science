package reference

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
	"sync"
)

// ReferenceGenomeParser handles parsing and indexing of reference genome FASTA files
type ReferenceGenomeParser struct {
	fastaPath    string
	indexPath    string
	chromSizes   map[string]int
	chromOffsets map[string]int64
	memoryMapped bool
	mu           sync.RWMutex
}

// Chromosome represents a single chromosome with its sequence and metadata
type Chromosome struct {
	Name      string
	Sequence  []byte
	Length    int
	GCContent float64
	NRegions  []GenomicInterval
}

// GenomicInterval represents a genomic interval [start, end)
type GenomicInterval struct {
	Start int
	End   int
}

// NewReferenceGenomeParser creates a new reference genome parser
func NewReferenceGenomeParser(fastaPath string) *ReferenceGenomeParser {
	return &ReferenceGenomeParser{
		fastaPath:    fastaPath,
		indexPath:    fastaPath + ".fai",
		chromSizes:   make(map[string]int),
		chromOffsets: make(map[string]int64),
		memoryMapped: true,
	}
}

// LoadIndex loads the FASTA index file (.fai) for rapid chromosome lookup
func (p *ReferenceGenomeParser) LoadIndex() error {
	p.mu.Lock()
	defer p.mu.Unlock()

	indexFile, err := os.Open(p.indexPath)
	if err != nil {
		return fmt.Errorf("failed to open index file %s: %w", p.indexPath, err)
	}
	defer indexFile.Close()

	scanner := bufio.NewScanner(indexFile)
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}

		chrName := fields[0]
		chrLength, err := strconv.Atoi(fields[1])
		if err != nil {
			continue
		}
		chrOffset, err := strconv.ParseInt(fields[2], 10, 64)
		if err != nil {
			continue
		}

		p.chromSizes[chrName] = chrLength
		p.chromOffsets[chrName] = chrOffset
	}

	return scanner.Err()
}

// LoadChromosome reads one chromosome into memory by scanning the FASTA from
// the start; the .fai offsets are recorded by LoadIndex but not used to seek,
// and nothing is memory-mapped. The chromosome must be listed in the index.
func (p *ReferenceGenomeParser) LoadChromosome(chrName string) (*Chromosome, error) {
	p.mu.RLock()
	defer p.mu.RUnlock()

	size, exists := p.chromSizes[chrName]
	if !exists {
		return nil, fmt.Errorf("chromosome %s not found in reference genome", chrName)
	}

	// Load chromosome sequence from FASTA file
	sequence, err := p.loadChromosomeSequence(chrName, size)
	if err != nil {
		return nil, fmt.Errorf("failed to load chromosome %s: %w", chrName, err)
	}

	// Create chromosome object
	chr := &Chromosome{
		Name:     chrName,
		Sequence: sequence,
		Length:   len(sequence),
	}

	// Calculate GC content
	chr.GCContent = p.calculateGCContent(sequence)

	// Identify N regions (gaps)
	chr.NRegions = p.identifyNRegions(sequence)

	return chr, nil
}

// loadChromosomeSequence loads the actual sequence data for a chromosome
func (p *ReferenceGenomeParser) loadChromosomeSequence(chrName string, expectedSize int) ([]byte, error) {
	file, err := os.Open(p.fastaPath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var sequence []byte
	foundChromosome := false

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			chrHeader := strings.TrimPrefix(line, ">")
			headerFields := strings.Fields(chrHeader)
			if len(headerFields) > 0 && headerFields[0] == chrName {
				foundChromosome = true
				continue
			} else if foundChromosome {
				// Found next chromosome, stop reading
				break
			}
		}

		if foundChromosome {
			// Convert to uppercase and append to sequence
			sequence = append(sequence, []byte(strings.ToUpper(line))...)
		}
	}

	if !foundChromosome {
		return nil, fmt.Errorf("chromosome %s not found in FASTA file", chrName)
	}

	return sequence, scanner.Err()
}

// calculateGCContent calculates the GC content of a sequence
func (p *ReferenceGenomeParser) calculateGCContent(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0.0
	}

	gcCount := 0
	totalCount := 0

	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
		if base == 'A' || base == 'T' || base == 'G' || base == 'C' {
			totalCount++
		}
	}

	if totalCount == 0 {
		return 0.0
	}

	return float64(gcCount) / float64(totalCount)
}

// identifyNRegions identifies regions of N bases (gaps) in the sequence
func (p *ReferenceGenomeParser) identifyNRegions(sequence []byte) []GenomicInterval {
	var nRegions []GenomicInterval
	inNRegion := false
	var start int

	for i, base := range sequence {
		if base == 'N' {
			if !inNRegion {
				start = i
				inNRegion = true
			}
		} else {
			if inNRegion {
				nRegions = append(nRegions, GenomicInterval{Start: start, End: i})
				inNRegion = false
			}
		}
	}

	// Handle N region at end of sequence
	if inNRegion {
		nRegions = append(nRegions, GenomicInterval{Start: start, End: len(sequence)})
	}

	return nRegions
}

// GetChromosomeNames returns a list of available chromosome names
func (p *ReferenceGenomeParser) GetChromosomeNames() []string {
	p.mu.RLock()
	defer p.mu.RUnlock()

	names := make([]string, 0, len(p.chromSizes))
	for name := range p.chromSizes {
		names = append(names, name)
	}
	return names
}

// GetChromosomeSize returns the size of a specific chromosome
func (p *ReferenceGenomeParser) GetChromosomeSize(chrName string) (int, error) {
	p.mu.RLock()
	defer p.mu.RUnlock()

	size, exists := p.chromSizes[chrName]
	if !exists {
		return 0, fmt.Errorf("chromosome %s not found", chrName)
	}
	return size, nil
}

// ValidateGenomeFile validates that the genome file exists and is accessible
func (p *ReferenceGenomeParser) ValidateGenomeFile() error {
	// Check if FASTA file exists
	if _, err := os.Stat(p.fastaPath); os.IsNotExist(err) {
		return fmt.Errorf("reference genome file not found: %s", p.fastaPath)
	}

	// Check if index file exists
	if _, err := os.Stat(p.indexPath); os.IsNotExist(err) {
		return fmt.Errorf("reference genome index file not found: %s", p.indexPath)
	}

	return nil
}

// CreateIndex creates a FASTA index file for rapid chromosome access
func (p *ReferenceGenomeParser) CreateIndex() error {
	fastaFile, err := os.Open(p.fastaPath)
	if err != nil {
		return fmt.Errorf("failed to open FASTA file: %w", err)
	}
	defer fastaFile.Close()

	indexFile, err := os.Create(p.indexPath)
	if err != nil {
		return fmt.Errorf("failed to create index file: %w", err)
	}
	defer indexFile.Close()

	scanner := bufio.NewScanner(fastaFile)
	var currentChr string
	var currentLength int
	var sequenceStartOffset int64
	var lineLength int
	var bytesPerLine int
	var currentOffset int64

	for scanner.Scan() {
		line := scanner.Text()

		if strings.HasPrefix(line, ">") {
			// Write previous chromosome to index
			if currentChr != "" {
				fmt.Fprintf(indexFile, "%s\t%d\t%d\t%d\t%d\n", // #nosec G705
					currentChr, currentLength, sequenceStartOffset,
					lineLength, bytesPerLine)
			}

			// Start new chromosome
			chrHeader := strings.TrimPrefix(line, ">")
			headerFields := strings.Fields(chrHeader)
			if len(headerFields) > 0 {
				currentChr = headerFields[0]
			}
			currentLength = 0
			lineLength = 0
			bytesPerLine = 0
			sequenceStartOffset = currentOffset + int64(len(line)+1) // Next line after header
		} else {
			// Process sequence line
			currentLength += len(line)
			if lineLength == 0 {
				lineLength = len(line)
				bytesPerLine = len(line) + 1 // +1 for newline
			}
		}

		currentOffset += int64(len(line) + 1) // +1 for newline
	}

	// Write last chromosome to index
	if currentChr != "" {
		fmt.Fprintf(indexFile, "%s\t%d\t%d\t%d\t%d\n", // #nosec G705
			currentChr, currentLength, sequenceStartOffset,
			lineLength, bytesPerLine)
	}

	return scanner.Err()
}

// GetGenomeStats returns basic statistics about the reference genome
func (p *ReferenceGenomeParser) GetGenomeStats() map[string]interface{} {
	p.mu.RLock()
	defer p.mu.RUnlock()

	totalSize := 0
	for _, size := range p.chromSizes {
		totalSize += size
	}

	return map[string]interface{}{
		"total_chromosomes": len(p.chromSizes),
		"total_size":        totalSize,
		"chromosome_sizes":  p.chromSizes,
		"fasta_path":        p.fastaPath,
		"index_path":        p.indexPath,
	}
}
