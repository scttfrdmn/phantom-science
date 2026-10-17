package reference

import (
	"fmt"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// ReferenceManager provides centralized management of reference genome data with lazy loading
type ReferenceManager struct {
	defaultLoader     *LazyReferenceLoader
	customLoaders     map[string]*LazyReferenceLoader
	loaderMutex       sync.RWMutex
	
	// Performance optimization
	enablePrefetching bool
	prefetchWindow    int64
	
	// Statistics
	totalRequests     int64
	totalCacheHits    int64
	totalLoadTime     int64
	
	// Configuration
	config           ManagerConfig
}

// ManagerConfig provides configuration for the reference manager
type ManagerConfig struct {
	DefaultReferenceFile string
	DefaultCacheSize     int64
	EnablePrefetching    bool
	PrefetchWindow       int64  // Bases to prefetch around requested region
	MaxLoaders           int    // Maximum number of concurrent loaders
	CompressionLevel     int
	EnablePackedArrays   bool
}

// NewReferenceManager creates a new reference genome manager
func NewReferenceManager(config ManagerConfig) (*ReferenceManager, error) {
	// Set defaults
	if config.DefaultCacheSize <= 0 {
		config.DefaultCacheSize = 4 * 1024 * 1024 * 1024 // 4GB default
	}
	if config.PrefetchWindow <= 0 {
		config.PrefetchWindow = 10000 // 10KB prefetch window
	}
	if config.MaxLoaders <= 0 {
		config.MaxLoaders = 10
	}
	if config.CompressionLevel < 0 || config.CompressionLevel > 9 {
		config.CompressionLevel = 6
	}
	
	manager := &ReferenceManager{
		customLoaders:     make(map[string]*LazyReferenceLoader),
		enablePrefetching: config.EnablePrefetching,
		prefetchWindow:    config.PrefetchWindow,
		config:           config,
	}
	
	// Create default loader if reference file specified
	if config.DefaultReferenceFile != "" {
		defaultConfig := LazyLoaderConfig{
			ReferenceFile:      config.DefaultReferenceFile,
			MaxCacheSize:       config.DefaultCacheSize,
			CompressionLevel:   config.CompressionLevel,
			EnablePackedArrays: config.EnablePackedArrays,
		}
		
		var err error
		manager.defaultLoader, err = NewLazyReferenceLoader(defaultConfig)
		if err != nil {
			return nil, fmt.Errorf("failed to create default reference loader: %w", err)
		}
	}
	
	return manager, nil
}

// RegisterReferenceGenome registers a custom reference genome with the manager
func (rm *ReferenceManager) RegisterReferenceGenome(name, filepath string, cacheSize int64) error {
	rm.loaderMutex.Lock()
	defer rm.loaderMutex.Unlock()
	
	// Check loader limit
	if len(rm.customLoaders) >= rm.config.MaxLoaders {
		return fmt.Errorf("maximum number of loaders (%d) reached", rm.config.MaxLoaders)
	}
	
	if cacheSize <= 0 {
		cacheSize = rm.config.DefaultCacheSize / 2 // Smaller cache for custom loaders
	}
	
	config := LazyLoaderConfig{
		ReferenceFile:      filepath,
		MaxCacheSize:       cacheSize,
		CompressionLevel:   rm.config.CompressionLevel,
		EnablePackedArrays: rm.config.EnablePackedArrays,
	}
	
	loader, err := NewLazyReferenceLoader(config)
	if err != nil {
		return fmt.Errorf("failed to create reference loader for %s: %w", name, err)
	}
	
	rm.customLoaders[name] = loader
	return nil
}

// GetSequence retrieves genomic sequence using the specified or default reference
func (rm *ReferenceManager) GetSequence(chromosome string, start, end int64, referenceName ...string) ([]byte, error) {
	startTime := time.Now()
	defer func() {
		rm.totalRequests++
		rm.totalLoadTime += time.Since(startTime).Nanoseconds()
	}()
	
	loader := rm.getLoader(referenceName...)
	if loader == nil {
		return nil, fmt.Errorf("no reference loader available")
	}
	
	// Apply prefetching if enabled
	actualStart, actualEnd := start, end
	if rm.enablePrefetching {
		actualStart = start - rm.prefetchWindow
		actualEnd = end + rm.prefetchWindow
		if actualStart < 0 {
			actualStart = 0
		}
	}
	
	sequence, err := loader.GetSequence(chromosome, actualStart, actualEnd)
	if err != nil {
		return nil, err
	}
	
	// Extract requested region if prefetching was applied
	if rm.enablePrefetching {
		prefetchOffset := start - actualStart
		requestedLength := end - start
		if prefetchOffset >= 0 && prefetchOffset+requestedLength <= int64(len(sequence)) {
			return sequence[prefetchOffset : prefetchOffset+requestedLength], nil
		}
	}
	
	return sequence, nil
}

// GetKmerProfile retrieves k-mer profile using the specified or default reference
func (rm *ReferenceManager) GetKmerProfile(chromosome string, start, end int64, kmerSize int, referenceName ...string) (*kmer.PackedFrequencyTable, error) {
	loader := rm.getLoader(referenceName...)
	if loader == nil {
		return nil, fmt.Errorf("no reference loader available")
	}
	
	return loader.GetKmerProfile(chromosome, start, end, kmerSize)
}

// getLoader returns the specified loader or default loader
func (rm *ReferenceManager) getLoader(referenceName ...string) *LazyReferenceLoader {
	if len(referenceName) > 0 && referenceName[0] != "" {
		rm.loaderMutex.RLock()
		loader, exists := rm.customLoaders[referenceName[0]]
		rm.loaderMutex.RUnlock()
		if exists {
			return loader
		}
	}
	
	return rm.defaultLoader
}

// WarmupReference preloads frequently accessed regions for better performance
func (rm *ReferenceManager) WarmupReference(regions []GenomicRegion, kmerSize int, referenceName ...string) error {
	loader := rm.getLoader(referenceName...)
	if loader == nil {
		return fmt.Errorf("no reference loader available for warmup")
	}
	
	return loader.WarmupCache(regions, kmerSize)
}

// GetCombinedStatistics returns aggregate statistics across all loaders
func (rm *ReferenceManager) GetCombinedStatistics() CombinedCacheStatistics {
	rm.loaderMutex.RLock()
	defer rm.loaderMutex.RUnlock()
	
	stats := CombinedCacheStatistics{
		TotalRequests: rm.totalRequests,
		AvgLoadTime:   0,
		LoaderStats:   make(map[string]CacheStatistics),
	}
	
	if rm.totalRequests > 0 {
		stats.AvgLoadTime = float64(rm.totalLoadTime) / float64(rm.totalRequests)
	}
	
	// Collect statistics from default loader
	if rm.defaultLoader != nil {
		defaultStats := rm.defaultLoader.GetCacheStatistics()
		stats.LoaderStats["default"] = defaultStats
		stats.TotalCacheHits += defaultStats.CacheHits
		stats.TotalCacheMisses += defaultStats.CacheMisses
		stats.TotalCacheSize += defaultStats.CacheSize
		stats.TotalCachedRegions += defaultStats.CachedRegions
	}
	
	// Collect statistics from custom loaders
	for name, loader := range rm.customLoaders {
		loaderStats := loader.GetCacheStatistics()
		stats.LoaderStats[name] = loaderStats
		stats.TotalCacheHits += loaderStats.CacheHits
		stats.TotalCacheMisses += loaderStats.CacheMisses
		stats.TotalCacheSize += loaderStats.CacheSize
		stats.TotalCachedRegions += loaderStats.CachedRegions
	}
	
	// Calculate overall hit ratio
	totalRequests := stats.TotalCacheHits + stats.TotalCacheMisses
	if totalRequests > 0 {
		stats.OverallHitRatio = float64(stats.TotalCacheHits) / float64(totalRequests)
	}
	
	return stats
}

// CombinedCacheStatistics contains aggregate statistics from all reference loaders
type CombinedCacheStatistics struct {
	TotalRequests       int64
	TotalCacheHits      int64
	TotalCacheMisses    int64
	OverallHitRatio     float64
	TotalCachedRegions  int64
	TotalCacheSize      int64
	AvgLoadTime         float64
	LoaderStats         map[string]CacheStatistics
}

// OptimizeCache performs cache optimization across all loaders
func (rm *ReferenceManager) OptimizeCache() {
	rm.loaderMutex.RLock()
	defer rm.loaderMutex.RUnlock()
	
	// Get current statistics
	stats := rm.GetCombinedStatistics()
	
	// If overall hit ratio is low, suggest cache size increases
	if stats.OverallHitRatio < 0.8 {
		fmt.Printf("Cache optimization: Hit ratio %.1f%% is below target. Consider increasing cache sizes.\n", 
			stats.OverallHitRatio*100)
	}
	
	// Check individual loader performance
	for name, loaderStats := range stats.LoaderStats {
		if loaderStats.CacheUtilization > 0.95 {
			fmt.Printf("Cache optimization: Loader '%s' is at %.1f%% utilization. Consider increasing cache size.\n",
				name, loaderStats.CacheUtilization*100)
		}
	}
}

// ClearAllCaches clears caches for all reference loaders
func (rm *ReferenceManager) ClearAllCaches() {
	rm.loaderMutex.RLock()
	defer rm.loaderMutex.RUnlock()
	
	if rm.defaultLoader != nil {
		rm.defaultLoader.ClearCache()
	}
	
	for _, loader := range rm.customLoaders {
		loader.ClearCache()
	}
	
	// Reset manager statistics
	rm.totalRequests = 0
	rm.totalCacheHits = 0
	rm.totalLoadTime = 0
}

// Close releases all resources used by the reference manager
func (rm *ReferenceManager) Close() error {
	rm.loaderMutex.Lock()
	defer rm.loaderMutex.Unlock()
	
	var errors []error
	
	// Close default loader
	if rm.defaultLoader != nil {
		if err := rm.defaultLoader.Close(); err != nil {
			errors = append(errors, fmt.Errorf("failed to close default loader: %w", err))
		}
	}
	
	// Close custom loaders
	for name, loader := range rm.customLoaders {
		if err := loader.Close(); err != nil {
			errors = append(errors, fmt.Errorf("failed to close loader %s: %w", name, err))
		}
	}
	
	// Clear loader maps
	rm.customLoaders = make(map[string]*LazyReferenceLoader)
	rm.defaultLoader = nil
	
	// Return first error if any occurred
	if len(errors) > 0 {
		return errors[0]
	}
	
	return nil
}

// GetReferenceInfo returns information about available reference genomes
func (rm *ReferenceManager) GetReferenceInfo() ReferenceInfo {
	rm.loaderMutex.RLock()
	defer rm.loaderMutex.RUnlock()
	
	info := ReferenceInfo{
		HasDefault:        rm.defaultLoader != nil,
		CustomReferences:  make([]string, 0, len(rm.customLoaders)),
		TotalLoaders:      len(rm.customLoaders),
		MaxLoaders:        rm.config.MaxLoaders,
	}
	
	for name := range rm.customLoaders {
		info.CustomReferences = append(info.CustomReferences, name)
	}
	
	if rm.defaultLoader != nil {
		info.TotalLoaders++
	}
	
	return info
}

// ReferenceInfo contains information about available reference genomes
type ReferenceInfo struct {
	HasDefault        bool
	CustomReferences  []string
	TotalLoaders      int
	MaxLoaders        int
}

// GetOptimalRegionSize calculates optimal region size based on access patterns
func (rm *ReferenceManager) GetOptimalRegionSize() int64 {
	stats := rm.GetCombinedStatistics()
	
	// Base calculation on cache hit ratio and average access patterns
	baseSize := int64(1024 * 1024) // 1MB base
	
	if stats.OverallHitRatio > 0.9 {
		// High hit ratio - can afford smaller regions
		return baseSize / 2
	} else if stats.OverallHitRatio < 0.7 {
		// Low hit ratio - need larger regions
		return baseSize * 2
	}
	
	return baseSize
}

// PrefetchRegions proactively loads regions that are likely to be accessed
func (rm *ReferenceManager) PrefetchRegions(regions []GenomicRegion, kmerSize int, referenceName ...string) error {
	loader := rm.getLoader(referenceName...)
	if loader == nil {
		return fmt.Errorf("no reference loader available for prefetching")
	}
	
	// Load regions asynchronously in background
	go func() {
		for _, region := range regions {
			// Load sequence
			_, _ = loader.GetSequence(region.Chromosome, region.Start, region.End)
			
			// Load k-mer profile if specified
			if kmerSize > 0 {
				_, _ = loader.GetKmerProfile(region.Chromosome, region.Start, region.End, kmerSize)
			}
		}
	}()
	
	return nil
}