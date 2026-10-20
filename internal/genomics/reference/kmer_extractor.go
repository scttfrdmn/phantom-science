package reference

import (
	"fmt"
	"strings"
	"sync"

	"github.com/scttfrdmn/phantom-science/internal/detmath"
)

// KmerExtractor handles k-mer frequency extraction from genomic sequences
type KmerExtractor struct {
	k              int
	frequencyTable map[string]uint64
	totalKmers     uint64
	mu             sync.RWMutex
	windowSize     int
	skipNRegions   bool
}

// KmerFrequencyTable represents k-mer frequencies for a given k
type KmerFrequencyTable struct {
	K           int
	Frequencies map[string]float64
	TotalKmers  uint64
	Coverage    float64
}

// NewKmerExtractor creates a new k-mer extractor
func NewKmerExtractor(k int) *KmerExtractor {
	return &KmerExtractor{
		k:              k,
		frequencyTable: make(map[string]uint64),
		totalKmers:     0,
		windowSize:     1000,
		skipNRegions:   true,
	}
}

// SetWindowSize stores a window size. Extraction does not use it; every
// k-mer position is counted.
func (e *KmerExtractor) SetWindowSize(size int) {
	e.mu.Lock()
	defer e.mu.Unlock()
	e.windowSize = size
}

// SetSkipNRegions sets whether to skip N regions during extraction
func (e *KmerExtractor) SetSkipNRegions(skip bool) {
	e.mu.Lock()
	defer e.mu.Unlock()
	e.skipNRegions = skip
}

// ExtractKmerFrequencies replaces the table with counts of every
// overlapping k-mer made only of upper-case A, C, G and T. K-mers with N or
// lower-case (soft-masked) bases are always skipped, whatever SetSkipNRegions
// says.
func (e *KmerExtractor) ExtractKmerFrequencies(sequence []byte) error {
	e.mu.Lock()
	defer e.mu.Unlock()

	if len(sequence) < e.k {
		return fmt.Errorf("sequence length %d is less than k=%d", len(sequence), e.k)
	}

	// Reset frequency table
	e.frequencyTable = make(map[string]uint64)
	e.totalKmers = 0

	// Sliding window k-mer enumeration
	for i := 0; i <= len(sequence)-e.k; i++ {
		kmer := string(sequence[i : i+e.k])

		// Skip k-mers containing N bases if configured
		if e.skipNRegions && strings.Contains(kmer, "N") {
			continue
		}

		// Only count valid DNA bases
		if e.isValidKmer(kmer) {
			e.frequencyTable[kmer]++
			e.totalKmers++
		}
	}

	return nil
}

// ExtractKmerFrequenciesFromChromosome extracts k-mer frequencies from a chromosome
func (e *KmerExtractor) ExtractKmerFrequenciesFromChromosome(chr *Chromosome) error {
	return e.ExtractKmerFrequencies(chr.Sequence)
}

// GetFrequencyTable returns the current frequency table as probabilities
func (e *KmerExtractor) GetFrequencyTable() *KmerFrequencyTable {
	e.mu.RLock()
	defer e.mu.RUnlock()

	frequencies := make(map[string]float64)
	for kmer, count := range e.frequencyTable {
		frequencies[kmer] = float64(count) / float64(e.totalKmers)
	}

	coverage := e.calculateCoverage()

	return &KmerFrequencyTable{
		K:           e.k,
		Frequencies: frequencies,
		TotalKmers:  e.totalKmers,
		Coverage:    coverage,
	}
}

// MergeFrequencyTables merges multiple frequency tables
func (e *KmerExtractor) MergeFrequencyTables(tables []*KmerFrequencyTable) (*KmerFrequencyTable, error) {
	if len(tables) == 0 {
		return nil, fmt.Errorf("no frequency tables to merge")
	}

	k := tables[0].K
	mergedFreqs := make(map[string]uint64)
	totalKmers := uint64(0)

	// Merge all tables
	for _, table := range tables {
		if table.K != k {
			return nil, fmt.Errorf("cannot merge tables with different k values: %d vs %d", k, table.K)
		}

		for kmer, freq := range table.Frequencies {
			count := uint64(float64(freq * float64(table.TotalKmers)))
			mergedFreqs[kmer] += count
			totalKmers += count
		}
	}

	// Convert to probabilities; guard against all counts truncating to zero
	if totalKmers == 0 {
		return &KmerFrequencyTable{K: k, Frequencies: make(map[string]float64), TotalKmers: 0}, nil
	}
	frequencies := make(map[string]float64)
	for kmer, count := range mergedFreqs {
		if count > 0 {
			frequencies[kmer] = float64(count) / float64(totalKmers)
		}
	}

	return &KmerFrequencyTable{
		K:           k,
		Frequencies: frequencies,
		TotalKmers:  totalKmers,
		Coverage:    e.calculateCoverageFromFreqs(frequencies, k),
	}, nil
}

// ComputeExpectedKmers computes the expected number of unique k-mers for given k
func (e *KmerExtractor) ComputeExpectedKmers(k int) int {
	// For DNA: 4^k possible k-mers
	return int(detmath.Pow(4, float64(k)))
}

// GetObservedKmers returns the number of observed unique k-mers
func (e *KmerExtractor) GetObservedKmers() int {
	e.mu.RLock()
	defer e.mu.RUnlock()
	return len(e.frequencyTable)
}

// isValidKmer checks if a k-mer contains only valid DNA bases
func (e *KmerExtractor) isValidKmer(kmer string) bool {
	for _, base := range kmer {
		if base != 'A' && base != 'T' && base != 'G' && base != 'C' {
			return false
		}
	}
	return true
}

// calculateCoverage calculates the coverage of possible k-mers
func (e *KmerExtractor) calculateCoverage() float64 {
	expected := e.ComputeExpectedKmers(e.k)
	observed := len(e.frequencyTable)
	return float64(observed) / float64(expected)
}

// calculateCoverageFromFreqs calculates coverage from frequency map
func (e *KmerExtractor) calculateCoverageFromFreqs(frequencies map[string]float64, k int) float64 {
	expected := e.ComputeExpectedKmers(k)
	observed := len(frequencies)
	return float64(observed) / float64(expected)
}

// GenerateAllPossibleKmers generates all possible k-mers for given k
func (e *KmerExtractor) GenerateAllPossibleKmers(k int) []string {
	bases := []string{"A", "T", "G", "C"}
	kmers := []string{""}

	for i := 0; i < k; i++ {
		var newKmers []string
		for _, kmer := range kmers {
			for _, base := range bases {
				newKmers = append(newKmers, kmer+base)
			}
		}
		kmers = newKmers
	}

	return kmers
}

// NormalizeFrequencies ensures frequencies sum to 1.0
func (e *KmerExtractor) NormalizeFrequencies(frequencies map[string]float64) map[string]float64 {
	sum := 0.0
	for _, freq := range frequencies {
		sum += freq
	}

	if sum == 0.0 {
		return frequencies
	}

	normalized := make(map[string]float64)
	for kmer, freq := range frequencies {
		normalized[kmer] = freq / sum
	}

	return normalized
}

// GetKmerComplexity calculates the complexity of k-mer distribution
func (e *KmerExtractor) GetKmerComplexity() float64 {
	e.mu.RLock()
	defer e.mu.RUnlock()

	if e.totalKmers == 0 {
		return 0.0
	}

	// Calculate Shannon entropy
	entropy := 0.0
	for _, count := range e.frequencyTable {
		prob := float64(count) / float64(e.totalKmers)
		if prob > 0 {
			entropy -= float64(prob * detmath.Log2(prob))
		}
	}

	// Normalize by maximum possible entropy
	maxEntropy := detmath.Log2(float64(e.ComputeExpectedKmers(e.k)))
	if maxEntropy == 0 {
		return 0.0
	}

	return entropy / maxEntropy
}

// KmerFrequency represents a k-mer with its frequency
type KmerFrequency struct {
	Kmer      string
	Frequency int
}

// GetTopKmers returns the most frequent k-mers
func (e *KmerExtractor) GetTopKmers(n int) []KmerFrequency {
	e.mu.RLock()
	defer e.mu.RUnlock()

	var kmerFreqs []KmerFrequency
	for kmer, count := range e.frequencyTable {
		kmerFreqs = append(kmerFreqs, KmerFrequency{
			Kmer:      kmer,
			Frequency: int(count), // #nosec G115
		})
	}

	// Sort by frequency (descending)
	for i := 0; i < len(kmerFreqs)-1; i++ {
		for j := i + 1; j < len(kmerFreqs); j++ {
			if kmerFreqs[i].Frequency < kmerFreqs[j].Frequency {
				kmerFreqs[i], kmerFreqs[j] = kmerFreqs[j], kmerFreqs[i]
			}
		}
	}

	// Return top n
	if n > len(kmerFreqs) {
		n = len(kmerFreqs)
	}

	return kmerFreqs[:n]
}

// CompareFrequencies compares two frequency tables and returns differences
func (e *KmerExtractor) CompareFrequencies(table1, table2 *KmerFrequencyTable) map[string]float64 {
	if table1.K != table2.K {
		return nil
	}

	differences := make(map[string]float64)

	// Get all unique k-mers from both tables
	allKmers := make(map[string]bool)
	for kmer := range table1.Frequencies {
		allKmers[kmer] = true
	}
	for kmer := range table2.Frequencies {
		allKmers[kmer] = true
	}

	// Calculate differences
	for kmer := range allKmers {
		freq1 := table1.Frequencies[kmer]
		freq2 := table2.Frequencies[kmer]
		differences[kmer] = freq2 - freq1
	}

	return differences
}

// GetKmerStats returns statistics about the k-mer distribution
func (e *KmerExtractor) GetKmerStats() map[string]interface{} {
	e.mu.RLock()
	defer e.mu.RUnlock()

	stats := make(map[string]interface{})
	stats["k"] = e.k
	stats["total_kmers"] = e.totalKmers
	stats["unique_kmers"] = len(e.frequencyTable)
	stats["expected_kmers"] = e.ComputeExpectedKmers(e.k)
	stats["coverage"] = e.calculateCoverage()
	stats["complexity"] = e.GetKmerComplexity()

	return stats
}
