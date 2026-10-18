package kmer

import (
	"fmt"
	"math"
	"sort"
	"unsafe"
)

// PackedFrequencyTable stores k-mer frequencies in a dense 4^k array whose
// bit width per entry grows to fit the largest count. It is smaller than a
// map when most k-mers occur and larger for sparse, high-k tables;
// GetCompressionRatio estimates the ratio. No benchmark measures speed.
type PackedFrequencyTable struct {
	K                 int
	PackedFrequencies *PackedArray // Packed binary array for k-mer frequencies
	TotalCount        int
	TotalKmers        int
	MaxKmers          int // 4^k maximum possible k-mers
	NonZeroCount      int // Number of k-mers with non-zero frequency
}

// PackedArray is a fixed-length array of unsigned integers packed at a
// shared bit width that widens on demand.
type PackedArray struct {
	data        []uint32 // Packed frequency data
	maxValue    uint32   // Maximum frequency value stored
	bitsPerItem int      // Bits required per frequency value
	itemCount   int      // Total number of items (4^k)
	mask        uint32   // Bit mask for extracting values
}

// NewPackedFrequencyTable creates a new packed frequency table for k-mer size k
func NewPackedFrequencyTable(k int) (*PackedFrequencyTable, error) {
	if k <= 0 || k > 16 {
		return nil, fmt.Errorf("k must be between 1 and 16, got %d", k)
	}

	maxKmers := int(math.Pow(4, float64(k)))
	packedArray, err := NewPackedArray(maxKmers, 32) // Start with 32-bit precision
	if err != nil {
		return nil, fmt.Errorf("failed to create packed array: %w", err)
	}

	return &PackedFrequencyTable{
		K:                 k,
		PackedFrequencies: packedArray,
		TotalCount:        0,
		TotalKmers:        0,
		MaxKmers:          maxKmers,
		NonZeroCount:      0,
	}, nil
}

// NewPackedArray creates a new packed array with specified item count and max bits per item
func NewPackedArray(itemCount int, maxBitsPerItem int) (*PackedArray, error) {
	if itemCount <= 0 {
		return nil, fmt.Errorf("itemCount must be positive, got %d", itemCount)
	}

	if maxBitsPerItem <= 0 || maxBitsPerItem > 32 {
		return nil, fmt.Errorf("maxBitsPerItem must be between 1 and 32, got %d", maxBitsPerItem)
	}

	// Start with minimum bits per item (1 bit) and expand as needed
	bitsPerItem := 1
	dataSize := calculatePackedArraySize(itemCount, bitsPerItem)

	return &PackedArray{
		data:        make([]uint32, dataSize),
		maxValue:    0,
		bitsPerItem: bitsPerItem,
		itemCount:   itemCount,
		mask:        (1 << bitsPerItem) - 1,
	}, nil
}

// calculatePackedArraySize calculates the required array size for packed storage
func calculatePackedArraySize(itemCount, bitsPerItem int) int {
	totalBits := itemCount * bitsPerItem
	return (totalBits + 31) / 32 // Round up to nearest uint32
}

// Get retrieves the frequency for a k-mer index
func (pa *PackedArray) Get(index int) uint32 {
	if index < 0 || index >= pa.itemCount {
		return 0
	}

	// Calculate bit position
	bitPos := index * pa.bitsPerItem
	wordIndex := bitPos / 32
	bitOffset := bitPos % 32

	if wordIndex >= len(pa.data) {
		return 0
	}

	// Handle values that span across uint32 boundaries
	if bitOffset+pa.bitsPerItem <= 32 {
		// Value fits in single uint32
		return (pa.data[wordIndex] >> bitOffset) & pa.mask
	} else {
		// Value spans two uint32s
		bitsInFirstWord := 32 - bitOffset
		bitsInSecondWord := pa.bitsPerItem - bitsInFirstWord

		if wordIndex+1 >= len(pa.data) {
			return (pa.data[wordIndex] >> bitOffset) & ((1 << bitsInFirstWord) - 1)
		}

		firstPart := (pa.data[wordIndex] >> bitOffset) & ((1 << bitsInFirstWord) - 1)
		secondPart := (pa.data[wordIndex+1] & ((1 << bitsInSecondWord) - 1)) << bitsInFirstWord

		return firstPart | secondPart
	}
}

// Set stores a frequency value for a k-mer index
func (pa *PackedArray) Set(index int, value uint32) error {
	if index < 0 || index >= pa.itemCount {
		return fmt.Errorf("index %d out of range [0, %d)", index, pa.itemCount)
	}

	// Check if we need to expand the array for larger values
	if value > pa.mask {
		if err := pa.expandBits(value); err != nil {
			return fmt.Errorf("failed to expand array: %w", err)
		}
	}

	// Update maximum value
	if value > pa.maxValue {
		pa.maxValue = value
	}

	// Calculate bit position
	bitPos := index * pa.bitsPerItem
	wordIndex := bitPos / 32
	bitOffset := bitPos % 32

	// Ensure data array is large enough
	requiredSize := calculatePackedArraySize(pa.itemCount, pa.bitsPerItem)
	if len(pa.data) < requiredSize {
		newData := make([]uint32, requiredSize)
		copy(newData, pa.data)
		pa.data = newData
	}

	// Clear existing value
	if bitOffset+pa.bitsPerItem <= 32 {
		// Value fits in single uint32
		pa.data[wordIndex] &= ^(pa.mask << bitOffset)
		pa.data[wordIndex] |= (value & pa.mask) << bitOffset
	} else {
		// Value spans two uint32s
		bitsInFirstWord := 32 - bitOffset
		bitsInSecondWord := pa.bitsPerItem - bitsInFirstWord

		// Clear and set first part
		firstMask := uint32(((1 << bitsInFirstWord) - 1) << bitOffset)
		pa.data[wordIndex] &= ^firstMask
		pa.data[wordIndex] |= uint32((value & ((1 << bitsInFirstWord) - 1)) << bitOffset)

		// Clear and set second part
		if wordIndex+1 < len(pa.data) {
			secondMask := uint32((1 << bitsInSecondWord) - 1)
			pa.data[wordIndex+1] &= ^secondMask
			pa.data[wordIndex+1] |= uint32((value >> bitsInFirstWord) & uint32(secondMask))
		}
	}

	return nil
}

// expandBits expands the packed array to accommodate larger values
func (pa *PackedArray) expandBits(newMaxValue uint32) error {
	// Calculate required bits
	newBitsPerItem := 32 - countLeadingZeros(newMaxValue)
	if newBitsPerItem <= pa.bitsPerItem {
		return nil // No expansion needed
	}

	// Create new expanded array
	oldData := pa.data
	oldBitsPerItem := pa.bitsPerItem
	oldMask := pa.mask

	pa.bitsPerItem = newBitsPerItem
	pa.mask = (1 << newBitsPerItem) - 1
	newSize := calculatePackedArraySize(pa.itemCount, pa.bitsPerItem)
	pa.data = make([]uint32, newSize)

	// Copy existing values to expanded array
	for i := 0; i < pa.itemCount; i++ {
		// Get value using old format
		bitPos := i * oldBitsPerItem
		wordIndex := bitPos / 32
		bitOffset := bitPos % 32

		if wordIndex >= len(oldData) {
			continue
		}

		var oldValue uint32
		if bitOffset+oldBitsPerItem <= 32 {
			oldValue = (oldData[wordIndex] >> bitOffset) & oldMask
		} else {
			bitsInFirstWord := 32 - bitOffset
			bitsInSecondWord := oldBitsPerItem - bitsInFirstWord

			firstPart := (oldData[wordIndex] >> bitOffset) & ((1 << bitsInFirstWord) - 1)
			var secondPart uint32
			if wordIndex+1 < len(oldData) {
				secondPart = (oldData[wordIndex+1] & ((1 << bitsInSecondWord) - 1)) << bitsInFirstWord
			}
			oldValue = firstPart | secondPart
		}

		// Set value using new format
		if oldValue > 0 {
			pa.Set(i, oldValue)
		}
	}

	return nil
}

// countLeadingZeros counts the number of leading zeros in a uint32
func countLeadingZeros(x uint32) int {
	if x == 0 {
		return 32
	}

	n := 0
	if x <= 0x0000FFFF {
		n += 16
		x <<= 16
	}
	if x <= 0x00FFFFFF {
		n += 8
		x <<= 8
	}
	if x <= 0x0FFFFFFF {
		n += 4
		x <<= 4
	}
	if x <= 0x3FFFFFFF {
		n += 2
		x <<= 2
	}
	if x <= 0x7FFFFFFF {
		n += 1
	}

	return n
}

// Add increments the frequency for a k-mer index
func (pa *PackedArray) Add(index int, delta uint32) error {
	currentValue := pa.Get(index)
	return pa.Set(index, currentValue+delta)
}

// GetMemoryUsage returns the memory usage of the packed array in bytes
func (pa *PackedArray) GetMemoryUsage() int64 {
	return int64(len(pa.data)*4 + int(unsafe.Sizeof(*pa)))
}

// IncrementKmer increments the frequency of a k-mer string
func (pft *PackedFrequencyTable) IncrementKmer(kmerString string) error {
	if len(kmerString) != pft.K {
		return fmt.Errorf("k-mer length %d does not match table k=%d", len(kmerString), pft.K)
	}

	index, err := KmerStringToIndex(kmerString)
	if err != nil {
		return fmt.Errorf("invalid k-mer string: %w", err)
	}

	oldValue := pft.PackedFrequencies.Get(index)
	if err := pft.PackedFrequencies.Add(index, 1); err != nil {
		return fmt.Errorf("failed to increment k-mer: %w", err)
	}

	// Update statistics
	pft.TotalCount++
	if oldValue == 0 {
		pft.NonZeroCount++
	}

	return nil
}

// GetFrequency returns the frequency of a k-mer string
func (pft *PackedFrequencyTable) GetFrequency(kmerString string) (int, error) {
	if len(kmerString) != pft.K {
		return 0, fmt.Errorf("k-mer length %d does not match table k=%d", len(kmerString), pft.K)
	}

	index, err := KmerStringToIndex(kmerString)
	if err != nil {
		return 0, fmt.Errorf("invalid k-mer string: %w", err)
	}

	return int(pft.PackedFrequencies.Get(index)), nil
}

// KmerStringToIndex converts a k-mer string to its packed array index
func KmerStringToIndex(kmer string) (int, error) {
	index := 0
	for i, base := range kmer {
		var baseValue int
		switch base {
		case 'A', 'a':
			baseValue = 0
		case 'T', 't':
			baseValue = 1
		case 'G', 'g':
			baseValue = 2
		case 'C', 'c':
			baseValue = 3
		default:
			return 0, fmt.Errorf("invalid DNA base '%c' at position %d", base, i)
		}
		index = index*4 + baseValue
	}
	return index, nil
}

// IndexToKmerString converts a packed array index to its k-mer string
func IndexToKmerString(index, k int) string {
	bases := []byte{'A', 'T', 'G', 'C'}
	kmer := make([]byte, k)

	for i := k - 1; i >= 0; i-- {
		kmer[i] = bases[index%4]
		index /= 4
	}

	return string(kmer)
}

// ToFrequencyTable converts the packed frequency table to a standard map-based frequency table
func (pft *PackedFrequencyTable) ToFrequencyTable() *FrequencyTable {
	frequencies := make(map[string]int)

	for i := 0; i < pft.MaxKmers; i++ {
		freq := pft.PackedFrequencies.Get(i)
		if freq > 0 {
			kmerString := IndexToKmerString(i, pft.K)
			frequencies[kmerString] = int(freq)
		}
	}

	return &FrequencyTable{
		K:           pft.K,
		Frequencies: frequencies,
		TotalCount:  pft.TotalCount,
		TotalKmers:  pft.TotalKmers,
	}
}

// FromFrequencyTable creates a packed frequency table from a standard frequency table
func FromFrequencyTable(ft *FrequencyTable) (*PackedFrequencyTable, error) {
	pft, err := NewPackedFrequencyTable(ft.K)
	if err != nil {
		return nil, err
	}

	// Copy all frequencies
	for kmerString, freq := range ft.Frequencies {
		index, err := KmerStringToIndex(kmerString)
		if err != nil {
			continue // Skip invalid k-mers
		}

		if err := pft.PackedFrequencies.Set(index, uint32(freq)); err != nil { // #nosec G115
			return nil, fmt.Errorf("failed to set frequency for %s: %w", kmerString, err)
		}

		if freq > 0 {
			pft.NonZeroCount++
		}
	}

	pft.TotalCount = ft.TotalCount
	pft.TotalKmers = ft.TotalKmers

	return pft, nil
}

// GetTopKmers returns the top N most frequent k-mers from the packed table
func (pft *PackedFrequencyTable) GetTopKmers(n int) []KmerFrequency {
	type kmerFreq struct {
		index int
		count uint32
	}

	var allKmers []kmerFreq
	for i := 0; i < pft.MaxKmers; i++ {
		freq := pft.PackedFrequencies.Get(i)
		if freq > 0 {
			allKmers = append(allKmers, kmerFreq{index: i, count: freq})
		}
	}

	// Sort by frequency in descending order
	sort.Slice(allKmers, func(i, j int) bool {
		return allKmers[i].count > allKmers[j].count
	})

	// Convert to result format
	if n > len(allKmers) {
		n = len(allKmers)
	}

	result := make([]KmerFrequency, n)
	for i := 0; i < n; i++ {
		kmerString := IndexToKmerString(allKmers[i].index, pft.K)
		result[i] = KmerFrequency{
			Kmer:      kmerString,
			Count:     int(allKmers[i].count),
			Frequency: float64(allKmers[i].count) / float64(pft.TotalCount),
		}
	}

	return result
}

// GetMemoryUsage returns the total memory usage of the packed frequency table
func (pft *PackedFrequencyTable) GetMemoryUsage() int64 {
	return pft.PackedFrequencies.GetMemoryUsage() + int64(unsafe.Sizeof(*pft))
}

// GetCompressionRatio returns an estimate of map-based memory (48 bytes per
// non-zero k-mer plus 8 bytes per possible k-mer) divided by the packed
// table's actual memory.
func (pft *PackedFrequencyTable) GetCompressionRatio() float64 {
	// Estimate map-based memory usage
	// Each map entry: ~48 bytes (key string + value int + map overhead)
	mapMemory := int64(pft.NonZeroCount * 48)

	// Add base map structure overhead (~40 bytes + bucket array)
	mapMemory += 40 + int64(pft.MaxKmers*8) // Approximate hash map overhead

	packedMemory := pft.GetMemoryUsage()

	if packedMemory == 0 {
		return 1.0
	}

	return float64(mapMemory) / float64(packedMemory)
}

// GetStatistics returns comprehensive statistics about the packed frequency table
func (pft *PackedFrequencyTable) GetStatistics() PackedTableStatistics {
	return PackedTableStatistics{
		K:                pft.K,
		MaxKmers:         pft.MaxKmers,
		NonZeroKmers:     pft.NonZeroCount,
		TotalCount:       pft.TotalCount,
		MemoryUsage:      pft.GetMemoryUsage(),
		CompressionRatio: pft.GetCompressionRatio(),
		BitsPerItem:      pft.PackedFrequencies.bitsPerItem,
		PackedArraySize:  len(pft.PackedFrequencies.data),
		MaxFrequency:     int(pft.PackedFrequencies.maxValue),
		Sparsity:         1.0 - (float64(pft.NonZeroCount) / float64(pft.MaxKmers)),
	}
}

// PackedTableStatistics contains statistics about packed frequency table usage
type PackedTableStatistics struct {
	K                int
	MaxKmers         int
	NonZeroKmers     int
	TotalCount       int
	MemoryUsage      int64
	CompressionRatio float64
	BitsPerItem      int
	PackedArraySize  int
	MaxFrequency     int
	Sparsity         float64 // Fraction of k-mers with zero frequency
}
