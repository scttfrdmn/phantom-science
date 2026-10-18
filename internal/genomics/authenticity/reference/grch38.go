// Package reference loads GRCh38 chromosomes from a directory of per-chromosome
// FASTA files (chr1.fa ... chrX.fa, chrY.fa, chrM.fa), computes their GC
// content and N gaps, and keeps a bounded number in memory with
// LazyLoadingManager. It does not read repeat or CpG-island annotations; the
// RepeatElement and CpGIsland types are not populated here.
package reference

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

// GRCh38 loads chromosomes from per-chromosome FASTA files under BasePath.
// It does not check that the files are actually GRCh38.
type GRCh38 struct {
	Assembly    string
	Chromosomes map[string]*Chromosome
	BasePath    string
	IsLoaded    bool
}

// Chromosome represents a single chromosome from the reference genome
type Chromosome struct {
	Name      string
	Sequence  []byte
	Length    int
	GCContent float64
	NRegions  []NRegion
	IsLoaded  bool
	FilePath  string
}

// NRegion represents a region of N bases (gaps) in the chromosome
type NRegion struct {
	Start  int
	End    int
	Length int
}

// RepeatElement represents a repeat element annotation
type RepeatElement struct {
	Start        int
	End          int
	RepeatType   string
	RepeatFamily string
	RepeatName   string
}

// CpGIsland represents a CpG island annotation
type CpGIsland struct {
	Start     int
	End       int
	Length    int
	GCContent float64
	CpGRatio  float64
}

// NewGRCh38 creates a new GRCh38 reference genome parser
func NewGRCh38(basePath string) *GRCh38 {
	return &GRCh38{
		Assembly:    "GRCh38",
		Chromosomes: make(map[string]*Chromosome),
		BasePath:    basePath,
		IsLoaded:    false,
	}
}

// LoadChromosome loads a specific chromosome from FASTA file
func (g *GRCh38) LoadChromosome(chrName string) (*Chromosome, error) {
	if chr, exists := g.Chromosomes[chrName]; exists && chr.IsLoaded {
		return chr, nil
	}

	// Construct the chromosome file path
	// Standard naming convention: chr1.fa, chr2.fa, etc.
	chrFile := fmt.Sprintf("chr%s.fa", chrName)
	if chrName == "23" {
		chrFile = "chrX.fa"
	} else if chrName == "24" {
		chrFile = "chrY.fa"
	} else if chrName == "MT" || chrName == "M" {
		chrFile = "chrM.fa"
	}

	filePath := fmt.Sprintf("%s/%s", g.BasePath, chrFile)

	chr := &Chromosome{
		Name:     chrName,
		FilePath: filePath,
		IsLoaded: false,
	}

	// Load the chromosome sequence from FASTA file
	sequence, err := g.loadFASTASequence(filePath)
	if err != nil {
		return nil, fmt.Errorf("failed to load chromosome %s: %w", chrName, err)
	}

	chr.Sequence = sequence
	chr.Length = len(sequence)
	chr.IsLoaded = true

	// Calculate GC content
	chr.GCContent = g.calculateGCContent(sequence)

	// Identify N regions (gaps)
	chr.NRegions = g.findNRegions(sequence)

	// Store in chromosomes map
	g.Chromosomes[chrName] = chr

	return chr, nil
}

// loadFASTASequence loads a DNA sequence from a FASTA file
func (g *GRCh38) loadFASTASequence(filePath string) ([]byte, error) {
	file, err := os.Open(filePath) // #nosec G304
	if err != nil {
		return nil, fmt.Errorf("failed to open FASTA file %s: %w", filePath, err)
	}
	defer file.Close()

	var sequence []byte
	scanner := bufio.NewScanner(file)

	// Skip header line(s) starting with '>'
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			continue
		}

		// Append sequence line, converting to uppercase
		sequence = append(sequence, []byte(strings.ToUpper(line))...)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading FASTA file: %w", err)
	}

	if len(sequence) == 0 {
		return nil, fmt.Errorf("empty sequence in FASTA file: %s", filePath)
	}

	return sequence, nil
}

// calculateGCContent calculates the GC content of a DNA sequence
func (g *GRCh38) calculateGCContent(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0.0
	}

	gcCount := 0
	totalCount := 0

	for _, base := range sequence {
		switch base {
		case 'G', 'C':
			gcCount++
			totalCount++
		case 'A', 'T':
			totalCount++
			// Skip 'N' bases from GC calculation
		}
	}

	if totalCount == 0 {
		return 0.0
	}

	return float64(gcCount) / float64(totalCount)
}

// findNRegions identifies regions of N bases (gaps) in the sequence
func (g *GRCh38) findNRegions(sequence []byte) []NRegion {
	var regions []NRegion
	var currentStart = -1

	for i, base := range sequence {
		if base == 'N' {
			if currentStart == -1 {
				currentStart = i
			}
		} else {
			if currentStart != -1 {
				// End of N region
				regions = append(regions, NRegion{
					Start:  currentStart,
					End:    i - 1,
					Length: i - currentStart,
				})
				currentStart = -1
			}
		}
	}

	// Handle N region at the end of sequence
	if currentStart != -1 {
		regions = append(regions, NRegion{
			Start:  currentStart,
			End:    len(sequence) - 1,
			Length: len(sequence) - currentStart,
		})
	}

	return regions
}

// GetSequence returns a subsequence from the specified chromosome
func (g *GRCh38) GetSequence(chrName string, start, end int) ([]byte, error) {
	chr, err := g.LoadChromosome(chrName)
	if err != nil {
		return nil, err
	}

	if start < 0 || end >= chr.Length || start > end {
		return nil, fmt.Errorf("invalid coordinates: start=%d, end=%d, chromosome length=%d", start, end, chr.Length)
	}

	return chr.Sequence[start : end+1], nil
}

// GetChromosomeNames returns a list of available chromosome names
func (g *GRCh38) GetChromosomeNames() []string {
	// Standard human chromosomes
	chromosomes := []string{
		"1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
		"11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
		"21", "22", "X", "Y", "MT",
	}
	return chromosomes
}

// GetChromosomeInfo returns basic information about a chromosome
func (g *GRCh38) GetChromosomeInfo(chrName string) (*ChromosomeInfo, error) {
	chr, err := g.LoadChromosome(chrName)
	if err != nil {
		return nil, err
	}

	return &ChromosomeInfo{
		Name:      chr.Name,
		Length:    chr.Length,
		GCContent: chr.GCContent,
		NRegions:  len(chr.NRegions),
	}, nil
}

// ChromosomeInfo contains summary information about a chromosome
type ChromosomeInfo struct {
	Name      string
	Length    int
	GCContent float64
	NRegions  int
}

// IsValidBase checks if a base is valid (A, T, G, C, N)
func IsValidBase(base byte) bool {
	switch base {
	case 'A', 'T', 'G', 'C', 'N':
		return true
	default:
		return false
	}
}

// ComplementBase returns the complement of a DNA base
func ComplementBase(base byte) byte {
	switch base {
	case 'A':
		return 'T'
	case 'T':
		return 'A'
	case 'G':
		return 'C'
	case 'C':
		return 'G'
	case 'N':
		return 'N'
	default:
		return 'N'
	}
}

// ReverseComplement returns the reverse complement of a DNA sequence
func ReverseComplement(sequence []byte) []byte {
	result := make([]byte, len(sequence))
	for i, base := range sequence {
		result[len(sequence)-1-i] = ComplementBase(base)
	}
	return result
}

// Close releases resources used by the reference genome
func (g *GRCh38) Close() error {
	// Clear chromosome data to free memory
	for chrName := range g.Chromosomes {
		g.Chromosomes[chrName].Sequence = nil
		g.Chromosomes[chrName].IsLoaded = false
	}
	g.IsLoaded = false
	return nil
}

// LazyLoadingManager keeps at most MaxLoadedChromosomes chromosome sequences
// in memory, unloading the earliest loaded first.
type LazyLoadingManager struct {
	MaxLoadedChromosomes int
	LoadedChromosomes    []string
	Reference            *GRCh38
}

// NewLazyLoadingManager creates a new lazy loading manager
func NewLazyLoadingManager(reference *GRCh38, maxLoaded int) *LazyLoadingManager {
	return &LazyLoadingManager{
		MaxLoadedChromosomes: maxLoaded,
		LoadedChromosomes:    make([]string, 0),
		Reference:            reference,
	}
}

// LoadChromosome loads a chromosome with lazy loading management
func (lm *LazyLoadingManager) LoadChromosome(chrName string) (*Chromosome, error) {
	// Check if chromosome is already loaded
	chr, exists := lm.Reference.Chromosomes[chrName]
	if exists && chr.IsLoaded {
		return chr, nil
	}

	// If we've reached the maximum number of loaded chromosomes, unload the oldest one
	if len(lm.LoadedChromosomes) >= lm.MaxLoadedChromosomes {
		oldestChr := lm.LoadedChromosomes[0]
		lm.UnloadChromosome(oldestChr)
		lm.LoadedChromosomes = lm.LoadedChromosomes[1:]
	}

	// Load the chromosome
	chr, err := lm.Reference.LoadChromosome(chrName)
	if err != nil {
		return nil, err
	}

	// Add to loaded chromosomes list
	lm.LoadedChromosomes = append(lm.LoadedChromosomes, chrName)

	return chr, nil
}

// UnloadChromosome unloads a chromosome from memory
func (lm *LazyLoadingManager) UnloadChromosome(chrName string) {
	if chr, exists := lm.Reference.Chromosomes[chrName]; exists {
		chr.Sequence = nil
		chr.IsLoaded = false
	}
}
