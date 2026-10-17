package kmer

import (
	"fmt"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/errors"
)

// PackedFrequencyExtractor extracts k-mer frequencies using memory-efficient packed arrays
type PackedFrequencyExtractor struct {
	K                  int
	IgnoreN            bool
	CaseSensitive      bool
	CircularSequence   bool
	MinFrequency       int
	UsePackedStorage   bool // Enable packed array storage
}

// NewPackedFrequencyExtractor creates a new k-mer frequency extractor using packed arrays
func NewPackedFrequencyExtractor(k int) *PackedFrequencyExtractor {
	return &PackedFrequencyExtractor{
		K:                k,
		IgnoreN:          true,
		CaseSensitive:    false,
		CircularSequence: false,
		MinFrequency:     1,
		UsePackedStorage: true,
	}
}

// ExtractFromSequence extracts k-mer frequencies from a single DNA sequence using packed storage
func (pfe *PackedFrequencyExtractor) ExtractFromSequence(sequence []byte) (*PackedFrequencyTable, error) {
	if len(sequence) < pfe.K {
		return nil, errors.NewValidationError("packed_kmer_extractor", "extract_from_sequence", fmt.Sprintf("sequence length (%d) is shorter than k-mer size (%d)", len(sequence), pfe.K))
	}

	// Create packed frequency table
	pft, err := NewPackedFrequencyTable(pfe.K)
	if err != nil {
		return nil, errors.WrapError(err, errors.MemoryError, "packed_kmer_extractor", "extract_from_sequence", "failed to create packed frequency table")
	}

	// Convert to uppercase if not case sensitive
	if !pfe.CaseSensitive {
		sequence = []byte(strings.ToUpper(string(sequence)))
	}

	// Calculate the number of k-mers to extract
	numKmers := len(sequence) - pfe.K + 1
	if pfe.CircularSequence {
		numKmers = len(sequence) // For circular sequences, we can extract len(sequence) k-mers
	}

	// Extract k-mers using sliding window
	for i := 0; i < numKmers; i++ {
		var kmer string

		if pfe.CircularSequence && i+pfe.K > len(sequence) {
			// Handle circular sequence wrapping
			kmer = string(sequence[i:]) + string(sequence[:i+pfe.K-len(sequence)])
		} else {
			kmer = string(sequence[i : i+pfe.K])
		}

		// Skip k-mers containing N bases if IgnoreN is true
		if pfe.IgnoreN && strings.Contains(kmer, "N") {
			continue
		}

		// Validate k-mer bases
		if !pfe.isValidKmer(kmer) {
			continue
		}

		// Increment k-mer frequency in packed table
		if err := pft.IncrementKmer(kmer); err != nil {
			continue // Skip invalid k-mers
		}

		pft.TotalKmers++
	}

	// Filter by minimum frequency if needed
	if pfe.MinFrequency > 1 {
		pft = pfe.applyMinimumFrequencyFilter(pft)
	}

	return pft, nil
}

// ExtractFromMultipleSequences extracts k-mer frequencies from multiple DNA sequences
func (pfe *PackedFrequencyExtractor) ExtractFromMultipleSequences(sequences [][]byte) (*PackedFrequencyTable, error) {
	if len(sequences) == 0 {
		return nil, errors.NewValidationError("packed_kmer_extractor", "extract_from_multiple_sequences", "no sequences provided")
	}

	// Create combined packed frequency table
	combinedTable, err := NewPackedFrequencyTable(pfe.K)
	if err != nil {
		return nil, errors.WrapError(err, errors.MemoryError, "packed_kmer_extractor", "extract_from_multiple_sequences", "failed to create combined frequency table")
	}

	// Process each sequence and merge results
	for i, sequence := range sequences {
		seqTable, err := pfe.ExtractFromSequence(sequence)
		if err != nil {
			return nil, errors.WrapError(err, errors.KmerExtractionError, "packed_kmer_extractor", "extract_from_multiple_sequences", fmt.Sprintf("failed to extract k-mers from sequence %d", i))
		}

		// Merge frequencies into combined table
		if err := pfe.mergePackedTables(combinedTable, seqTable); err != nil {
			return nil, errors.WrapError(err, errors.MemoryError, "packed_kmer_extractor", "extract_from_multiple_sequences", fmt.Sprintf("failed to merge sequence %d frequencies", i))
		}
	}

	return combinedTable, nil
}

// mergePackedTables merges frequencies from source table into destination table
func (pfe *PackedFrequencyExtractor) mergePackedTables(dst, src *PackedFrequencyTable) error {
	if dst.K != src.K {
		return errors.NewValidationError("packed_kmer_extractor", "merge_packed_tables", fmt.Sprintf("cannot merge tables with different k values: %d vs %d", dst.K, src.K))
	}

	// Merge all non-zero frequencies
	for i := 0; i < src.MaxKmers; i++ {
		srcFreq := src.PackedFrequencies.Get(i)
		if srcFreq > 0 {
			dstFreq := dst.PackedFrequencies.Get(i)
			if err := dst.PackedFrequencies.Set(i, dstFreq+srcFreq); err != nil {
				return errors.WrapError(err, errors.MemoryError, "packed_kmer_extractor", "merge_packed_tables", fmt.Sprintf("failed to merge frequency at index %d", i))
			}

			// Update statistics
			if dstFreq == 0 {
				dst.NonZeroCount++
			}
		}
	}

	// Update totals
	dst.TotalCount += src.TotalCount
	dst.TotalKmers += src.TotalKmers

	return nil
}

// applyMinimumFrequencyFilter filters out k-mers below minimum frequency
func (pfe *PackedFrequencyExtractor) applyMinimumFrequencyFilter(pft *PackedFrequencyTable) *PackedFrequencyTable {
	filteredTable, err := NewPackedFrequencyTable(pft.K)
	if err != nil {
		return pft // Return original on error
	}

	filteredCount := 0
	filteredTotalCount := 0

	// Copy only k-mers meeting minimum frequency requirement
	for i := 0; i < pft.MaxKmers; i++ {
		freq := pft.PackedFrequencies.Get(i)
		if freq >= uint32(pfe.MinFrequency) { // #nosec G115
			if err := filteredTable.PackedFrequencies.Set(i, freq); err != nil {
				continue // Skip on error
			}
			filteredCount++
			filteredTotalCount += int(freq)
		}
	}

	filteredTable.NonZeroCount = filteredCount
	filteredTable.TotalCount = filteredTotalCount
	filteredTable.TotalKmers = pft.TotalKmers

	return filteredTable
}

// isValidKmer checks if a k-mer contains only valid DNA bases
func (pfe *PackedFrequencyExtractor) isValidKmer(kmer string) bool {
	for _, base := range kmer {
		switch base {
		case 'A', 'T', 'G', 'C':
			continue
		case 'N':
			if !pfe.IgnoreN {
				continue
			}
			return false
		default:
			return false
		}
	}
	return true
}

// ExtractFromSequenceStreaming extracts k-mers from a sequence in streaming fashion
// to handle very large sequences without loading everything into memory
func (pfe *PackedFrequencyExtractor) ExtractFromSequenceStreaming(sequenceReader SequenceReader, bufferSize int) (*PackedFrequencyTable, error) {
	pft, err := NewPackedFrequencyTable(pfe.K)
	if err != nil {
		return nil, errors.WrapError(err, errors.MemoryError, "packed_kmer_extractor", "extract_from_sequence_streaming", "failed to create packed frequency table")
	}

	if bufferSize < pfe.K {
		bufferSize = pfe.K * 1000 // Ensure buffer is large enough
	}

	buffer := make([]byte, bufferSize)
	overlap := make([]byte, pfe.K-1) // Overlap from previous chunk
	overlapLen := 0

	for {
		n, err := sequenceReader.Read(buffer[overlapLen:])
		if n == 0 {
			break
		}

		// Add overlap from previous chunk
		if overlapLen > 0 {
			copy(buffer[:overlapLen], overlap)
		}

		totalLen := overlapLen + n
		currentBuffer := buffer[:totalLen]

		// Process k-mers in current buffer
		if err := pfe.processBufferKmers(pft, currentBuffer); err != nil {
			return nil, errors.WrapError(err, errors.KmerExtractionError, "packed_kmer_extractor", "extract_from_sequence_streaming", "failed to process buffer k-mers")
		}

		// Save overlap for next iteration
		if totalLen >= pfe.K-1 {
			overlapLen = pfe.K - 1
			copy(overlap, currentBuffer[totalLen-overlapLen:])
		} else {
			overlapLen = totalLen
			copy(overlap, currentBuffer)
		}

		if err != nil {
			break // End of stream
		}
	}

	return pft, nil
}

// processBufferKmers processes k-mers from a buffer
func (pfe *PackedFrequencyExtractor) processBufferKmers(pft *PackedFrequencyTable, buffer []byte) error {
	if len(buffer) < pfe.K {
		return nil
	}

	// Convert to uppercase if needed
	if !pfe.CaseSensitive {
		for i := range buffer {
			if buffer[i] >= 'a' && buffer[i] <= 'z' {
				buffer[i] -= 32
			}
		}
	}

	// Extract k-mers from buffer
	for i := 0; i <= len(buffer)-pfe.K; i++ {
		kmer := string(buffer[i : i+pfe.K])

		// Skip k-mers containing N bases if IgnoreN is true
		if pfe.IgnoreN && strings.Contains(kmer, "N") {
			continue
		}

		// Validate k-mer bases
		if !pfe.isValidKmer(kmer) {
			continue
		}

		// Increment k-mer frequency
		if err := pft.IncrementKmer(kmer); err != nil {
			continue // Skip invalid k-mers
		}

		pft.TotalKmers++
	}

	return nil
}

// SequenceReader interface for streaming sequence processing
type SequenceReader interface {
	Read([]byte) (int, error)
}

// ExtractFromSequenceBatch processes sequences in batches to optimize memory usage
func (pfe *PackedFrequencyExtractor) ExtractFromSequenceBatch(sequences [][]byte, batchSize int) (*PackedFrequencyTable, error) {
	if len(sequences) == 0 {
		return nil, fmt.Errorf("no sequences provided")
	}

	if batchSize <= 0 {
		batchSize = 1000 // Default batch size
	}

	// Create result table
	resultTable, err := NewPackedFrequencyTable(pfe.K)
	if err != nil {
		return nil, fmt.Errorf("failed to create result table: %w", err)
	}

	// Process sequences in batches
	for i := 0; i < len(sequences); i += batchSize {
		end := i + batchSize
		if end > len(sequences) {
			end = len(sequences)
		}

		batch := sequences[i:end]
		batchTable, err := pfe.ExtractFromMultipleSequences(batch)
		if err != nil {
			return nil, fmt.Errorf("failed to process batch %d-%d: %w", i, end-1, err)
		}

		// Merge batch results
		if err := pfe.mergePackedTables(resultTable, batchTable); err != nil {
			return nil, fmt.Errorf("failed to merge batch %d-%d: %w", i, end-1, err)
		}
	}

	return resultTable, nil
}

// PackedFrequencyExtractorConfig provides configuration options for packed extraction
type PackedFrequencyExtractorConfig struct {
	K                       int
	IgnoreN                 bool
	CaseSensitive          bool
	CircularSequence       bool
	MinFrequency           int
	UsePackedStorage       bool
	StreamingBufferSize    int
	BatchSize              int
	EnableMemoryOptimization bool
}

// NewPackedFrequencyExtractorWithConfig creates an extractor with custom configuration
func NewPackedFrequencyExtractorWithConfig(config PackedFrequencyExtractorConfig) *PackedFrequencyExtractor {
	extractor := &PackedFrequencyExtractor{
		K:                config.K,
		IgnoreN:          config.IgnoreN,
		CaseSensitive:    config.CaseSensitive,
		CircularSequence: config.CircularSequence,
		MinFrequency:     config.MinFrequency,
		UsePackedStorage: config.UsePackedStorage,
	}

	// Set defaults for unspecified values
	if extractor.MinFrequency == 0 {
		extractor.MinFrequency = 1
	}

	return extractor
}

// GetMemoryEstimate estimates memory usage for processing given sequences
func (pfe *PackedFrequencyExtractor) GetMemoryEstimate(numSequences, avgSequenceLength int) MemoryEstimate {
	maxKmers := 1
	for i := 0; i < pfe.K; i++ {
		maxKmers *= 4
	}

	// Estimate number of unique k-mers (rough approximation)
	totalKmers := numSequences * (avgSequenceLength - pfe.K + 1)
	estimatedUniqueKmers := int(float64(maxKmers) * (1.0 - float64(maxKmers-totalKmers)/float64(maxKmers)))
	if estimatedUniqueKmers > maxKmers {
		estimatedUniqueKmers = maxKmers
	}

	// Calculate packed array memory usage
	bitsPerItem := 32 // Start with maximum, will be optimized down
	arraySize := calculatePackedArraySize(maxKmers, bitsPerItem)
	packedMemory := int64(arraySize * 4) // 4 bytes per uint32

	// Calculate equivalent map memory usage
	mapMemory := int64(estimatedUniqueKmers * 48) // ~48 bytes per map entry

	return MemoryEstimate{
		MaxKmers:              maxKmers,
		EstimatedUniqueKmers:  estimatedUniqueKmers,
		PackedMemoryUsage:     packedMemory,
		MapMemoryUsage:        mapMemory,
		MemorySavings:         mapMemory - packedMemory,
		CompressionRatio:      float64(mapMemory) / float64(packedMemory),
	}
}

// MemoryEstimate contains memory usage estimates for packed k-mer extraction
type MemoryEstimate struct {
	MaxKmers              int
	EstimatedUniqueKmers  int
	PackedMemoryUsage     int64
	MapMemoryUsage        int64
	MemorySavings         int64
	CompressionRatio      float64
}

// GetOptimalBatchSize calculates optimal batch size based on available memory
func (pfe *PackedFrequencyExtractor) GetOptimalBatchSize(availableMemoryMB int, avgSequenceLength int) int {
	estimate := pfe.GetMemoryEstimate(1, avgSequenceLength)
	memoryPerSequence := estimate.PackedMemoryUsage
	
	availableBytes := int64(availableMemoryMB) * 1024 * 1024
	batchSize := int(availableBytes / memoryPerSequence)
	
	// Ensure reasonable bounds
	if batchSize < 100 {
		batchSize = 100
	} else if batchSize > 10000 {
		batchSize = 10000
	}
	
	return batchSize
}