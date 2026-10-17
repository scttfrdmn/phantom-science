package kmer

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/errors"
)

// FrequencyTable represents k-mer frequency data
type FrequencyTable struct {
	K           int
	Frequencies map[string]int
	TotalCount  int
	TotalKmers  int
}

// FrequencyExtractor extracts k-mer frequencies from DNA sequences
type FrequencyExtractor struct {
	K                  int
	IgnoreN            bool
	CaseSensitive      bool
	CircularSequence   bool
	MinFrequency       int
}

// NewFrequencyExtractor creates a new k-mer frequency extractor
func NewFrequencyExtractor(k int) *FrequencyExtractor {
	return &FrequencyExtractor{
		K:                k,
		IgnoreN:          true,
		CaseSensitive:    false,
		CircularSequence: false,
		MinFrequency:     1,
	}
}

// ExtractFromSequence extracts k-mer frequencies from a single DNA sequence
func (fe *FrequencyExtractor) ExtractFromSequence(sequence []byte) (*FrequencyTable, error) {
	if len(sequence) < fe.K {
		return nil, errors.NewValidationError("kmer_frequency_extractor", "extract_from_sequence", fmt.Sprintf("sequence length (%d) is shorter than k-mer size (%d)", len(sequence), fe.K))
	}

	frequencies := make(map[string]int)
	totalCount := 0

	// Convert to uppercase if not case sensitive
	if !fe.CaseSensitive {
		sequence = []byte(strings.ToUpper(string(sequence)))
	}

	// Calculate the number of k-mers to extract
	numKmers := len(sequence) - fe.K + 1
	if fe.CircularSequence {
		numKmers = len(sequence) // For circular sequences, we can extract len(sequence) k-mers
	}

	// Extract k-mers using sliding window
	for i := 0; i < numKmers; i++ {
		var kmer string
		
		if fe.CircularSequence && i+fe.K > len(sequence) {
			// Handle circular sequence wrapping
			kmer = string(sequence[i:]) + string(sequence[:i+fe.K-len(sequence)])
		} else {
			kmer = string(sequence[i:i+fe.K])
		}

		// Skip k-mers containing N bases if IgnoreN is true
		if fe.IgnoreN && strings.Contains(kmer, "N") {
			continue
		}

		// Validate k-mer bases
		if !fe.isValidKmer(kmer) {
			continue
		}

		frequencies[kmer]++
		totalCount++
	}

	// Filter by minimum frequency
	if fe.MinFrequency > 1 {
		for kmer, count := range frequencies {
			if count < fe.MinFrequency {
				delete(frequencies, kmer)
				totalCount -= count
			}
		}
	}

	return &FrequencyTable{
		K:           fe.K,
		Frequencies: frequencies,
		TotalCount:  totalCount,
		TotalKmers:  numKmers,
	}, nil
}

// ExtractFromMultipleSequences extracts k-mer frequencies from multiple DNA sequences
func (fe *FrequencyExtractor) ExtractFromMultipleSequences(sequences [][]byte) (*FrequencyTable, error) {
	if len(sequences) == 0 {
		return nil, errors.NewValidationError("kmer_frequency_extractor", "extract_from_multiple_sequences", "no sequences provided")
	}

	combinedFrequencies := make(map[string]int)
	totalCount := 0
	totalKmers := 0

	for i, sequence := range sequences {
		freqTable, err := fe.ExtractFromSequence(sequence)
		if err != nil {
			return nil, errors.WrapError(err, errors.KmerExtractionError, "kmer_frequency_extractor", "extract_from_multiple_sequences", fmt.Sprintf("failed to extract k-mers from sequence %d", i))
		}

		// Combine frequencies
		for kmer, count := range freqTable.Frequencies {
			combinedFrequencies[kmer] += count
		}
		
		totalCount += freqTable.TotalCount
		totalKmers += freqTable.TotalKmers
	}

	return &FrequencyTable{
		K:           fe.K,
		Frequencies: combinedFrequencies,
		TotalCount:  totalCount,
		TotalKmers:  totalKmers,
	}, nil
}

// isValidKmer checks if a k-mer contains only valid DNA bases
func (fe *FrequencyExtractor) isValidKmer(kmer string) bool {
	for _, base := range kmer {
		switch base {
		case 'A', 'T', 'G', 'C':
			continue
		case 'N':
			if !fe.IgnoreN {
				continue
			}
			return false
		default:
			return false
		}
	}
	return true
}

// GetProbabilityDistribution returns the probability distribution of k-mers
func (ft *FrequencyTable) GetProbabilityDistribution() map[string]float64 {
	probabilities := make(map[string]float64)
	
	if ft.TotalCount == 0 {
		return probabilities
	}

	for kmer, count := range ft.Frequencies {
		probabilities[kmer] = float64(count) / float64(ft.TotalCount)
	}

	return probabilities
}

// GetTopKmers returns the top N most frequent k-mers
func (ft *FrequencyTable) GetTopKmers(n int) []KmerFrequency {
	var kmerFreqs []KmerFrequency
	
	for kmer, count := range ft.Frequencies {
		kmerFreqs = append(kmerFreqs, KmerFrequency{
			Kmer:      kmer,
			Count:     count,
			Frequency: float64(count) / float64(ft.TotalCount),
		})
	}

	// Sort by count in descending order
	sort.Slice(kmerFreqs, func(i, j int) bool {
		return kmerFreqs[i].Count > kmerFreqs[j].Count
	})

	// Return top N
	if n > len(kmerFreqs) {
		n = len(kmerFreqs)
	}

	return kmerFreqs[:n]
}

// KmerFrequency represents a k-mer with its frequency information
type KmerFrequency struct {
	Kmer      string
	Count     int
	Frequency float64
}

// GetExpectedKmers returns all possible k-mers for the given k
func GetExpectedKmers(k int) []string {
	if k <= 0 {
		return nil
	}

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

// GetKmerComplexity calculates the complexity of a k-mer frequency distribution
func (ft *FrequencyTable) GetKmerComplexity() float64 {
	if ft.TotalCount == 0 {
		return 0.0
	}

	expectedKmers := int(math.Pow(4, float64(ft.K)))
	observedKmers := len(ft.Frequencies)

	return float64(observedKmers) / float64(expectedKmers)
}

// GetEntropyShannon calculates the Shannon entropy of the k-mer distribution
func (ft *FrequencyTable) GetEntropyShannon() float64 {
	if ft.TotalCount == 0 {
		return 0.0
	}

	entropy := 0.0
	for _, count := range ft.Frequencies {
		if count > 0 {
			p := float64(count) / float64(ft.TotalCount)
			entropy -= p * math.Log2(p)
		}
	}

	return entropy
}

// GetMaxEntropy returns the maximum possible entropy for k-mers of size k
func GetMaxEntropy(k int) float64 {
	return float64(k) * 2.0 // log2(4) = 2 for 4 DNA bases
}

// GetNormalizedEntropy returns the normalized Shannon entropy (0-1 scale)
func (ft *FrequencyTable) GetNormalizedEntropy() float64 {
	shannon := ft.GetEntropyShannon()
	maxEntropy := GetMaxEntropy(ft.K)
	
	if maxEntropy == 0 {
		return 0.0
	}
	
	return shannon / maxEntropy
}

// MergeFrequencyTables merges multiple frequency tables
func MergeFrequencyTables(tables ...*FrequencyTable) (*FrequencyTable, error) {
	if len(tables) == 0 {
		return nil, errors.NewValidationError("kmer_frequency_extractor", "merge_frequency_tables", "no frequency tables provided")
	}

	// Check that all tables have the same k
	k := tables[0].K
	for i, table := range tables {
		if table.K != k {
			return nil, errors.NewValidationError("kmer_frequency_extractor", "merge_frequency_tables", fmt.Sprintf("mismatched k values: table 0 has k=%d, table %d has k=%d", k, i, table.K))
		}
	}

	merged := &FrequencyTable{
		K:           k,
		Frequencies: make(map[string]int),
		TotalCount:  0,
		TotalKmers:  0,
	}

	for _, table := range tables {
		for kmer, count := range table.Frequencies {
			merged.Frequencies[kmer] += count
		}
		merged.TotalCount += table.TotalCount
		merged.TotalKmers += table.TotalKmers
	}

	return merged, nil
}

// GetGCContent calculates the GC content of the k-mer frequency table
func (ft *FrequencyTable) GetGCContent() float64 {
	if ft.TotalCount == 0 {
		return 0.0
	}

	gcCount := 0
	for kmer, count := range ft.Frequencies {
		for _, base := range kmer {
			if base == 'G' || base == 'C' {
				gcCount += count
			}
		}
	}

	totalBases := ft.TotalCount * ft.K
	return float64(gcCount) / float64(totalBases)
}

// GetKmerStats returns comprehensive statistics about the k-mer frequency table
func (ft *FrequencyTable) GetKmerStats() KmerStats {
	_ = ft.GetProbabilityDistribution()
	
	var counts []int
	for _, count := range ft.Frequencies {
		counts = append(counts, count)
	}
	
	sort.Ints(counts)
	
	stats := KmerStats{
		K:                 ft.K,
		TotalKmers:        ft.TotalKmers,
		UniqueKmers:       len(ft.Frequencies),
		TotalCount:        ft.TotalCount,
		Complexity:        ft.GetKmerComplexity(),
		ShannonEntropy:    ft.GetEntropyShannon(),
		NormalizedEntropy: ft.GetNormalizedEntropy(),
		GCContent:         ft.GetGCContent(),
	}

	if len(counts) > 0 {
		stats.MinCount = counts[0]
		stats.MaxCount = counts[len(counts)-1]
		
		// Calculate median
		if len(counts)%2 == 0 {
			stats.MedianCount = float64(counts[len(counts)/2-1]+counts[len(counts)/2]) / 2.0
		} else {
			stats.MedianCount = float64(counts[len(counts)/2])
		}
		
		// Calculate mean
		sum := 0
		for _, count := range counts {
			sum += count
		}
		stats.MeanCount = float64(sum) / float64(len(counts))
	}

	return stats
}

// KmerStats contains comprehensive statistics about a k-mer frequency table
type KmerStats struct {
	K                 int
	TotalKmers        int
	UniqueKmers       int
	TotalCount        int
	MinCount          int
	MaxCount          int
	MeanCount         float64
	MedianCount       float64
	Complexity        float64
	ShannonEntropy    float64
	NormalizedEntropy float64
	GCContent         float64
}

// FilterByFrequency filters k-mers by frequency range
func (ft *FrequencyTable) FilterByFrequency(minFreq, maxFreq int) *FrequencyTable {
	filtered := &FrequencyTable{
		K:           ft.K,
		Frequencies: make(map[string]int),
		TotalCount:  0,
		TotalKmers:  ft.TotalKmers,
	}

	for kmer, count := range ft.Frequencies {
		if count >= minFreq && count <= maxFreq {
			filtered.Frequencies[kmer] = count
			filtered.TotalCount += count
		}
	}

	return filtered
}

// GetRareKmers returns k-mers that appear less than the threshold
func (ft *FrequencyTable) GetRareKmers(threshold int) []string {
	var rare []string
	
	for kmer, count := range ft.Frequencies {
		if count < threshold {
			rare = append(rare, kmer)
		}
	}
	
	sort.Strings(rare)
	return rare
}

// GetAbundantKmers returns k-mers that appear more than the threshold
func (ft *FrequencyTable) GetAbundantKmers(threshold int) []string {
	var abundant []string
	
	for kmer, count := range ft.Frequencies {
		if count > threshold {
			abundant = append(abundant, kmer)
		}
	}
	
	sort.Strings(abundant)
	return abundant
}