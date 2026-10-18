package reference

import (
	"fmt"
	"os"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/errors"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// LazyReferenceLoader is an LRU cache of reference regions. It does not read
// the reference file: it only checks that the file exists, and every region
// is placeholder sequence generated from the coordinates (see
// loadRegionSequence). Use FastaParser or ReferenceGenomeParser for real
// sequence.
type LazyReferenceLoader struct {
	referenceFile      string
	maxCacheSize       int64
	compressionLevel   int
	enablePackedArrays bool

	// Cache management
	regionCache      map[string]*CachedRegion
	cacheMutex       sync.RWMutex
	currentCacheSize int64
	accessCounter    int64

	// Statistics
	cacheHits        int64
	cacheMisses      int64
	loadingTime      int64 // nanoseconds
	compressionRatio float64
}

// CachedRegion represents a cached genomic region with lazy loading capabilities
type CachedRegion struct {
	Chromosome       string
	StartPos         int64
	EndPos           int64
	Sequence         []byte
	KmerProfile      *kmer.PackedFrequencyTable
	LastAccessed     int64
	LoadTime         int64
	Size             int64
	CompressionRatio float64
	IsLoaded         bool
	LoadMutex        sync.RWMutex
}

// LazyLoaderConfig provides configuration for a LazyReferenceLoader.
type LazyLoaderConfig struct {
	ReferenceFile      string   // Must exist; its contents are not read
	MaxCacheSize       int64    // Maximum cache size in bytes
	RegionSize         int64    // Defaulted but not used
	CompressionLevel   int      // Stored but not used (nothing is compressed)
	EnablePackedArrays bool     // Use packed arrays for k-mer storage
	PreloadRegions     []string // Regions to preload on initialization
	KmerSize           int      // K-mer size for profile generation
}

// NewLazyReferenceLoader creates a new lazy reference genome loader
func NewLazyReferenceLoader(config LazyLoaderConfig) (*LazyReferenceLoader, error) {
	// Validate configuration
	if config.ReferenceFile == "" {
		return nil, errors.NewValidationError("lazy_reference_loader", "new_loader", "reference file is required")
	}

	// Check if reference file exists
	if _, err := os.Stat(config.ReferenceFile); os.IsNotExist(err) {
		return nil, errors.NewReferenceError("lazy_reference_loader", "new_loader", fmt.Sprintf("reference file does not exist: %s", config.ReferenceFile))
	}

	// Set defaults
	if config.MaxCacheSize <= 0 {
		config.MaxCacheSize = 2 * 1024 * 1024 * 1024 // 2GB default
	}
	if config.RegionSize <= 0 {
		config.RegionSize = 10 * 1024 * 1024 // 10MB regions
	}
	if config.CompressionLevel < 0 || config.CompressionLevel > 9 {
		config.CompressionLevel = 6 // Balanced compression
	}
	if config.KmerSize <= 0 {
		config.KmerSize = 8
	}

	loader := &LazyReferenceLoader{
		referenceFile:      config.ReferenceFile,
		maxCacheSize:       config.MaxCacheSize,
		compressionLevel:   config.CompressionLevel,
		enablePackedArrays: config.EnablePackedArrays,
		regionCache:        make(map[string]*CachedRegion),
	}

	// Preload specified regions if configured
	if len(config.PreloadRegions) > 0 {
		if err := loader.preloadRegions(config.PreloadRegions, config.KmerSize); err != nil {
			return nil, errors.WrapError(err, errors.ReferenceLoadError, "lazy_reference_loader", "new_loader", "failed to preload regions")
		}
	}

	return loader, nil
}

// GetSequence returns the cached sequence for the region, generating and
// caching placeholder sequence on a miss. The result is not reference
// sequence.
func (lrl *LazyReferenceLoader) GetSequence(chromosome string, start, end int64) ([]byte, error) {
	regionKey := fmt.Sprintf("%s:%d-%d", chromosome, start, end)

	// Check cache first
	lrl.cacheMutex.RLock()
	region, found := lrl.regionCache[regionKey]
	lrl.cacheMutex.RUnlock()

	if found && region.IsLoaded {
		lrl.cacheHits++
		lrl.updateAccessTime(region)
		return region.Sequence, nil
	}

	lrl.cacheMisses++

	// Load region on demand
	sequence, err := lrl.loadRegionSequence(chromosome, start, end)
	if err != nil {
		return nil, errors.WrapError(err, errors.ReferenceLoadError, "lazy_reference_loader", "get_sequence", fmt.Sprintf("failed to load sequence for %s", regionKey))
	}

	// Cache the loaded region
	lrl.cacheRegion(regionKey, chromosome, start, end, sequence)

	return sequence, nil
}

// GetKmerProfile returns the k-mer profile of the region's (placeholder)
// sequence, caching it with the region.
func (lrl *LazyReferenceLoader) GetKmerProfile(chromosome string, start, end int64, kmerSize int) (*kmer.PackedFrequencyTable, error) {
	regionKey := fmt.Sprintf("%s:%d-%d", chromosome, start, end)

	// Check if region is cached with k-mer profile
	lrl.cacheMutex.RLock()
	region, found := lrl.regionCache[regionKey]
	lrl.cacheMutex.RUnlock()

	if found && region.IsLoaded && region.KmerProfile != nil {
		lrl.cacheHits++
		lrl.updateAccessTime(region)
		return region.KmerProfile, nil
	}

	// Load sequence if not already loaded
	sequence, err := lrl.GetSequence(chromosome, start, end)
	if err != nil {
		return nil, errors.WrapError(err, errors.ReferenceLoadError, "lazy_reference_loader", "get_kmer_profile", "failed to get sequence for k-mer profile")
	}

	// Generate k-mer profile
	kmerProfile, err := lrl.generateKmerProfile(sequence, kmerSize)
	if err != nil {
		return nil, errors.WrapError(err, errors.KmerExtractionError, "lazy_reference_loader", "get_kmer_profile", "failed to generate k-mer profile")
	}

	// Update cached region with k-mer profile
	lrl.cacheMutex.Lock()
	if region, exists := lrl.regionCache[regionKey]; exists {
		region.KmerProfile = kmerProfile
	}
	lrl.cacheMutex.Unlock()

	return kmerProfile, nil
}

// loadRegionSequence is a stub: it does not read the reference file. It
// generates i.i.d. bases (A/T 0.30, G/C 0.20) from a linear congruential
// generator seeded by start and len(chromosome), so chromosomes whose names
// have the same length (chr1, chr2, ...) get identical sequence.
func (lrl *LazyReferenceLoader) loadRegionSequence(chromosome string, start, end int64) ([]byte, error) {
	startTime := time.Now()
	defer func() {
		lrl.loadingTime += time.Since(startTime).Nanoseconds()
	}()

	// Placeholder sequence; no file access
	sequenceLength := end - start
	if sequenceLength <= 0 || sequenceLength > 100*1024*1024 { // Max 100MB regions
		return nil, errors.NewValidationError("lazy_reference_loader", "load_region_sequence", fmt.Sprintf("invalid sequence length: %d", sequenceLength))
	}

	sequence := make([]byte, sequenceLength)
	bases := []byte{'A', 'T', 'G', 'C'}

	// Base composition close to the human genome (~41% GC); no other
	// structure (no dinucleotide bias, repeats or CpG depletion)
	baseProbs := []float64{0.30, 0.30, 0.20, 0.20} // A, T, G, C
	cumProbs := make([]float64, 4)
	cumProbs[0] = baseProbs[0]
	for i := 1; i < 4; i++ {
		cumProbs[i] = cumProbs[i-1] + baseProbs[i]
	}

	// Deterministic seed from start and the chromosome name's length
	seed := uint64(start*31 + int64(len(chromosome))*17) // #nosec G115

	for i := int64(0); i < sequenceLength; i++ {
		// Simple linear congruential generator for deterministic random bases
		seed = (seed*1103515245 + 12345) & 0x7fffffff
		rand := float64(seed) / float64(0x7fffffff)

		baseIndex := 0
		for j := 0; j < 4; j++ {
			if rand <= cumProbs[j] {
				baseIndex = j
				break
			}
		}
		sequence[i] = bases[baseIndex]
	}

	return sequence, nil
}

// generateKmerProfile generates k-mer frequency profile from sequence
func (lrl *LazyReferenceLoader) generateKmerProfile(sequence []byte, kmerSize int) (*kmer.PackedFrequencyTable, error) {
	if lrl.enablePackedArrays {
		// Use packed array k-mer extractor for memory efficiency
		extractor := kmer.NewPackedFrequencyExtractor(kmerSize)
		return extractor.ExtractFromSequence(sequence)
	}

	// Same packed extractor; enablePackedArrays has no effect
	extractor := kmer.NewPackedFrequencyExtractor(kmerSize)
	return extractor.ExtractFromSequence(sequence)
}

// cacheRegion adds a region to the cache with eviction management
func (lrl *LazyReferenceLoader) cacheRegion(regionKey, chromosome string, start, end int64, sequence []byte) {
	region := &CachedRegion{
		Chromosome:   chromosome,
		StartPos:     start,
		EndPos:       end,
		Sequence:     sequence,
		LastAccessed: lrl.getAccessCounter(),
		LoadTime:     time.Now().UnixNano(),
		Size:         int64(len(sequence)),
		IsLoaded:     true,
	}

	lrl.cacheMutex.Lock()
	defer lrl.cacheMutex.Unlock()

	// Check cache size and evict if necessary
	lrl.currentCacheSize += region.Size
	if lrl.currentCacheSize > lrl.maxCacheSize {
		lrl.evictLeastRecentlyUsed()
	}

	lrl.regionCache[regionKey] = region
}

// evictLeastRecentlyUsed removes least recently used regions from cache
func (lrl *LazyReferenceLoader) evictLeastRecentlyUsed() {
	type regionEntry struct {
		key    string
		region *CachedRegion
	}

	var regions []regionEntry
	for key, region := range lrl.regionCache {
		regions = append(regions, regionEntry{key, region})
	}

	// Sort by LastAccessed (ascending - oldest first)
	for i := 0; i < len(regions)-1; i++ {
		for j := i + 1; j < len(regions); j++ {
			if regions[i].region.LastAccessed > regions[j].region.LastAccessed {
				regions[i], regions[j] = regions[j], regions[i]
			}
		}
	}

	// Remove oldest regions until under cache limit
	targetSize := lrl.maxCacheSize * 3 / 4 // Evict to 75% of limit
	for _, entry := range regions {
		if lrl.currentCacheSize <= targetSize {
			break
		}

		delete(lrl.regionCache, entry.key)
		lrl.currentCacheSize -= entry.region.Size
	}
}

// updateAccessTime updates the access time for a cached region
func (lrl *LazyReferenceLoader) updateAccessTime(region *CachedRegion) {
	region.LoadMutex.Lock()
	region.LastAccessed = lrl.getAccessCounter()
	region.LoadMutex.Unlock()
}

// getAccessCounter returns an incremented access counter for LRU tracking
func (lrl *LazyReferenceLoader) getAccessCounter() int64 {
	lrl.accessCounter++
	return lrl.accessCounter
}

// preloadRegions preloads specified genomic regions into cache
func (lrl *LazyReferenceLoader) preloadRegions(regions []string, kmerSize int) error {
	for _, regionStr := range regions {
		// Parse region string (format: chr:start-end)
		var chromosome string
		var start, end int64

		if _, err := fmt.Sscanf(regionStr, "%[^:]:%d-%d", &chromosome, &start, &end); err != nil {
			return fmt.Errorf("invalid region format %s: %w", regionStr, err)
		}

		// Preload sequence
		if _, err := lrl.GetSequence(chromosome, start, end); err != nil {
			return fmt.Errorf("failed to preload region %s: %w", regionStr, err)
		}

		// Preload k-mer profile if specified
		if kmerSize > 0 {
			if _, err := lrl.GetKmerProfile(chromosome, start, end, kmerSize); err != nil {
				return fmt.Errorf("failed to preload k-mer profile for %s: %w", regionStr, err)
			}
		}
	}

	return nil
}

// GetCacheStatistics returns comprehensive cache performance statistics
func (lrl *LazyReferenceLoader) GetCacheStatistics() CacheStatistics {
	lrl.cacheMutex.RLock()
	defer lrl.cacheMutex.RUnlock()

	totalRequests := lrl.cacheHits + lrl.cacheMisses
	hitRatio := 0.0
	if totalRequests > 0 {
		hitRatio = float64(lrl.cacheHits) / float64(totalRequests)
	}

	avgLoadTime := 0.0
	if lrl.cacheMisses > 0 {
		avgLoadTime = float64(lrl.loadingTime) / float64(lrl.cacheMisses)
	}

	cacheUtil := 0.0
	if lrl.maxCacheSize > 0 {
		cacheUtil = float64(lrl.currentCacheSize) / float64(lrl.maxCacheSize)
	}

	return CacheStatistics{
		CacheHits:        lrl.cacheHits,
		CacheMisses:      lrl.cacheMisses,
		HitRatio:         hitRatio,
		CachedRegions:    int64(len(lrl.regionCache)),
		CacheSize:        lrl.currentCacheSize,
		MaxCacheSize:     lrl.maxCacheSize,
		CacheUtilization: cacheUtil,
		AvgLoadTime:      avgLoadTime, // nanoseconds per miss; 0 when no misses have occurred
	}
}

// CacheStatistics contains comprehensive cache performance metrics
type CacheStatistics struct {
	CacheHits        int64
	CacheMisses      int64
	HitRatio         float64
	CachedRegions    int64
	CacheSize        int64
	MaxCacheSize     int64
	CacheUtilization float64
	AvgLoadTime      float64 // Average loading time in nanoseconds
}

// WarmupCache preloads frequently accessed genomic regions
func (lrl *LazyReferenceLoader) WarmupCache(warmupRegions []GenomicRegion, kmerSize int) error {
	for _, region := range warmupRegions {
		// Load sequence
		if _, err := lrl.GetSequence(region.Chromosome, region.Start, region.End); err != nil {
			return fmt.Errorf("failed to warmup region %s:%d-%d: %w",
				region.Chromosome, region.Start, region.End, err)
		}

		// Generate k-mer profile if requested
		if kmerSize > 0 {
			if _, err := lrl.GetKmerProfile(region.Chromosome, region.Start, region.End, kmerSize); err != nil {
				return fmt.Errorf("failed to warmup k-mer profile for %s:%d-%d: %w",
					region.Chromosome, region.Start, region.End, err)
			}
		}
	}

	return nil
}

// GenomicRegion represents a genomic coordinate range
type GenomicRegion struct {
	Chromosome string
	Start      int64
	End        int64
}

// ClearCache removes all cached regions and resets statistics
func (lrl *LazyReferenceLoader) ClearCache() {
	lrl.cacheMutex.Lock()
	defer lrl.cacheMutex.Unlock()

	lrl.regionCache = make(map[string]*CachedRegion)
	lrl.currentCacheSize = 0
	lrl.cacheHits = 0
	lrl.cacheMisses = 0
	lrl.loadingTime = 0
	lrl.accessCounter = 0
}

// Close releases resources used by the lazy loader
func (lrl *LazyReferenceLoader) Close() error {
	lrl.ClearCache()
	return nil
}
