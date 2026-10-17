package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

func main() {
	var (
		referenceFile    = flag.String("reference", "", "Reference genome file (required)")
		cacheSize        = flag.Int64("cache", 2147483648, "Cache size in bytes (2GB default)")
		kmerSize         = flag.Int("k", 8, "K-mer size for profile generation")
		numRegions       = flag.Int("regions", 1000, "Number of test regions")
		regionSize       = flag.Int64("region-size", 10000, "Size of test regions in bases")
		benchmarkType    = flag.String("benchmark", "all", "Benchmark type: cache, prefetch, kmer, or all")
		verbose          = flag.Bool("verbose", false, "Enable verbose output")
		enablePrefetch   = flag.Bool("prefetch", true, "Enable prefetching optimization")
		compressionLevel = flag.Int("compression", 6, "Compression level (0-9)")
		warmupRegions    = flag.Int("warmup", 100, "Number of regions to warmup cache")
	)
	flag.Parse()

	fmt.Printf("🧬 Phantom Science Reference Genome Lazy Loading Benchmark\n")
	fmt.Printf("==========================================================\n\n")

	// Validate required parameters
	if *referenceFile == "" {
		// Create a test reference file for benchmarking
		*referenceFile = createTestReferenceFile()
		fmt.Printf("📁 Using generated test reference: %s\n", *referenceFile)
		defer os.Remove(*referenceFile) // Clean up after benchmarking
	}

	// Display configuration
	displayConfiguration(*referenceFile, *cacheSize, *kmerSize, *numRegions, *regionSize,
		*benchmarkType, *enablePrefetch, *compressionLevel, *warmupRegions)

	// Create reference manager
	managerConfig := reference.ManagerConfig{
		DefaultReferenceFile: *referenceFile,
		DefaultCacheSize:     *cacheSize,
		EnablePrefetching:    *enablePrefetch,
		CompressionLevel:     *compressionLevel,
		EnablePackedArrays:   true,
	}

	manager, err := reference.NewReferenceManager(managerConfig)
	if err != nil {
		log.Fatalf("Failed to create reference manager: %v", err)
	}
	defer manager.Close()

	// Run specified benchmarks
	switch *benchmarkType {
	case "cache":
		runCacheBenchmark(manager, *numRegions, *regionSize, *warmupRegions, *verbose)
	case "prefetch":
		runPrefetchBenchmark(manager, *numRegions, *regionSize, *verbose)
	case "kmer":
		runKmerBenchmark(manager, *numRegions, *regionSize, *kmerSize, *verbose)
	case "all":
		runComprehensiveBenchmark(manager, *numRegions, *regionSize, *kmerSize, *warmupRegions, *verbose)
	default:
		log.Fatalf("Unknown benchmark type: %s", *benchmarkType)
	}

	// Display final statistics
	displayFinalStatistics(manager)
}

func displayConfiguration(referenceFile string, cacheSize int64, kmerSize, numRegions int,
	regionSize int64, benchmarkType string, enablePrefetch bool, compressionLevel, warmupRegions int) {

	fmt.Printf("📁 Reference file: %s\n", referenceFile)
	fmt.Printf("💾 Cache size: %.1f MB\n", float64(cacheSize)/(1024*1024))
	fmt.Printf("🧬 K-mer size: %d\n", kmerSize)
	fmt.Printf("📍 Test regions: %d (%.1f KB each)\n", numRegions, float64(regionSize)/1024)
	fmt.Printf("🔧 Benchmark type: %s\n", benchmarkType)
	fmt.Printf("🚀 Prefetching: %t\n", enablePrefetch)
	fmt.Printf("🗜️ Compression level: %d\n", compressionLevel)
	fmt.Printf("🔥 Warmup regions: %d\n", warmupRegions)
	fmt.Printf("\n")
}

func createTestReferenceFile() string {
	// Create a temporary test reference file
	file, err := os.CreateTemp("", "test_reference_*.fasta")
	if err != nil {
		log.Fatalf("Failed to create test reference file: %v", err)
	}

	// Write basic FASTA content for testing
	chromosomes := []string{"chr1", "chr2", "chr3", "chr4", "chr5"}
	bases := []byte{'A', 'T', 'G', 'C'}

	for _, chr := range chromosomes {
		fmt.Fprintf(file, ">%s\n", chr)
		
		// Generate 10MB of sequence per chromosome
		sequenceLength := 10 * 1024 * 1024
		for i := 0; i < sequenceLength; i++ {
			fmt.Fprintf(file, "%c", bases[i%4])
			if (i+1)%80 == 0 {
				fmt.Fprintf(file, "\n")
			}
		}
		fmt.Fprintf(file, "\n")
	}

	file.Close()
	return file.Name()
}

func runCacheBenchmark(manager *reference.ReferenceManager, numRegions int, regionSize int64, warmupRegions int, verbose bool) {
	fmt.Printf("🗄️ Running Cache Performance Benchmark\n")
	fmt.Printf("======================================\n\n")

	chromosomes := []string{"chr1", "chr2", "chr3", "chr4", "chr5"}
	
	// Phase 1: Cold cache performance
	fmt.Printf("🧊 Phase 1: Cold Cache Performance\n")
	coldStart := time.Now()
	
	for i := 0; i < numRegions/2; i++ {
		chr := chromosomes[i%len(chromosomes)]
		start := int64(i * 50000) // Spread regions across genome
		end := start + regionSize
		
		_, err := manager.GetSequence(chr, start, end)
		if err != nil && verbose {
			fmt.Printf("Warning: Failed to get sequence for %s:%d-%d: %v\n", chr, start, end, err)
		}
	}
	
	coldDuration := time.Since(coldStart)
	fmt.Printf("   Cold cache time: %v (%.2f regions/sec)\n", 
		coldDuration, float64(numRegions/2)/coldDuration.Seconds())

	// Get statistics after cold cache
	coldStats := manager.GetCombinedStatistics()
	fmt.Printf("   Cold cache hit ratio: %.1f%%\n", coldStats.OverallHitRatio*100)

	// Phase 2: Warmup cache
	fmt.Printf("\n🔥 Phase 2: Cache Warmup\n")
	warmupRegionsList := make([]reference.GenomicRegion, warmupRegions)
	for i := 0; i < warmupRegions; i++ {
		chr := chromosomes[i%len(chromosomes)]
		start := int64(i * 10000) // Overlapping regions for better hit ratio
		end := start + regionSize
		warmupRegionsList[i] = reference.GenomicRegion{
			Chromosome: chr,
			Start:      start,
			End:        end,
		}
	}
	
	warmupStart := time.Now()
	err := manager.WarmupReference(warmupRegionsList, 0)
	if err != nil {
		fmt.Printf("Warning: Cache warmup failed: %v\n", err)
	}
	warmupDuration := time.Since(warmupStart)
	fmt.Printf("   Warmup time: %v\n", warmupDuration)

	// Phase 3: Hot cache performance
	fmt.Printf("\n🔥 Phase 3: Hot Cache Performance\n")
	hotStart := time.Now()
	
	// Access the same regions that were warmed up
	for _, region := range warmupRegionsList {
		_, err := manager.GetSequence(region.Chromosome, region.Start, region.End)
		if err != nil && verbose {
			fmt.Printf("Warning: Failed to get sequence for %s:%d-%d: %v\n", 
				region.Chromosome, region.Start, region.End, err)
		}
	}
	
	hotDuration := time.Since(hotStart)
	fmt.Printf("   Hot cache time: %v (%.2f regions/sec)\n", 
		hotDuration, float64(warmupRegions)/hotDuration.Seconds())

	// Get final statistics
	hotStats := manager.GetCombinedStatistics()
	fmt.Printf("   Hot cache hit ratio: %.1f%%\n", hotStats.OverallHitRatio*100)
	
	// Performance improvement
	speedupRatio := coldDuration.Seconds() / hotDuration.Seconds()
	fmt.Printf("   Cache speedup: %.2fx\n", speedupRatio)
	
	fmt.Printf("\n")
}

func runPrefetchBenchmark(manager *reference.ReferenceManager, numRegions int, regionSize int64, verbose bool) {
	fmt.Printf("🚀 Running Prefetch Performance Benchmark\n")
	fmt.Printf("==========================================\n\n")

	chromosomes := []string{"chr1", "chr2", "chr3", "chr4", "chr5"}
	
	// Create test regions for prefetching
	testRegions := make([]reference.GenomicRegion, numRegions)
	for i := 0; i < numRegions; i++ {
		chr := chromosomes[i%len(chromosomes)]
		start := int64(i * 20000) // Sequential regions for prefetch testing
		end := start + regionSize
		testRegions[i] = reference.GenomicRegion{
			Chromosome: chr,
			Start:      start,
			End:        end,
		}
	}

	// Phase 1: Sequential access without prefetching
	fmt.Printf("📈 Phase 1: Sequential Access (No Prefetch)\n")
	sequentialStart := time.Now()
	
	for _, region := range testRegions[:numRegions/2] {
		_, err := manager.GetSequence(region.Chromosome, region.Start, region.End)
		if err != nil && verbose {
			fmt.Printf("Warning: Failed to get sequence for %s:%d-%d: %v\n", 
				region.Chromosome, region.Start, region.End, err)
		}
	}
	
	sequentialDuration := time.Since(sequentialStart)
	fmt.Printf("   Sequential time: %v (%.2f regions/sec)\n", 
		sequentialDuration, float64(numRegions/2)/sequentialDuration.Seconds())

	// Phase 2: Access with prefetching
	fmt.Printf("\n🚀 Phase 2: Prefetch-Optimized Access\n")
	
	// Start prefetching
	prefetchRegions := testRegions[numRegions/2:]
	prefetchStart := time.Now()
	
	err := manager.PrefetchRegions(prefetchRegions, 0)
	if err != nil {
		fmt.Printf("Warning: Prefetch failed: %v\n", err)
	}
	
	// Small delay to allow prefetching to work
	time.Sleep(100 * time.Millisecond)
	
	// Access prefetched regions
	for _, region := range prefetchRegions {
		_, err := manager.GetSequence(region.Chromosome, region.Start, region.End)
		if err != nil && verbose {
			fmt.Printf("Warning: Failed to get sequence for %s:%d-%d: %v\n", 
				region.Chromosome, region.Start, region.End, err)
		}
	}
	
	prefetchDuration := time.Since(prefetchStart)
	fmt.Printf("   Prefetch time: %v (%.2f regions/sec)\n", 
		prefetchDuration, float64(len(prefetchRegions))/prefetchDuration.Seconds())

	// Performance comparison
	if sequentialDuration > 0 && prefetchDuration > 0 {
		prefetchSpeedup := sequentialDuration.Seconds() / prefetchDuration.Seconds()
		fmt.Printf("   Prefetch speedup: %.2fx\n", prefetchSpeedup)
	}
	
	fmt.Printf("\n")
}

func runKmerBenchmark(manager *reference.ReferenceManager, numRegions int, regionSize int64, kmerSize int, verbose bool) {
	fmt.Printf("🧬 Running K-mer Profile Benchmark\n")
	fmt.Printf("==================================\n\n")

	chromosomes := []string{"chr1", "chr2", "chr3", "chr4", "chr5"}
	
	// Phase 1: K-mer profile generation
	fmt.Printf("🔬 Phase 1: K-mer Profile Generation\n")
	kmerStart := time.Now()
	
	var totalKmers int64
	var totalMemory int64
	
	for i := 0; i < numRegions; i++ {
		chr := chromosomes[i%len(chromosomes)]
		start := int64(i * 30000)
		end := start + regionSize
		
		profile, err := manager.GetKmerProfile(chr, start, end, kmerSize)
		if err != nil {
			if verbose {
				fmt.Printf("Warning: Failed to get k-mer profile for %s:%d-%d: %v\n", chr, start, end, err)
			}
			continue
		}
		
		if profile != nil {
			totalKmers += int64(profile.TotalCount)
			totalMemory += profile.GetMemoryUsage()
		}
	}
	
	kmerDuration := time.Since(kmerStart)
	fmt.Printf("   K-mer generation time: %v\n", kmerDuration)
	fmt.Printf("   Total k-mers processed: %d\n", totalKmers)
	fmt.Printf("   Memory usage: %.2f MB\n", float64(totalMemory)/(1024*1024))
	fmt.Printf("   K-mer processing rate: %.0f k-mers/sec\n", 
		float64(totalKmers)/kmerDuration.Seconds())

	// Phase 2: K-mer profile caching efficiency
	fmt.Printf("\n🗄️ Phase 2: K-mer Profile Caching\n")
	cacheStart := time.Now()
	
	// Access the same regions again to test caching
	for i := 0; i < numRegions/2; i++ {
		chr := chromosomes[i%len(chromosomes)]
		start := int64(i * 30000)
		end := start + regionSize
		
		_, err := manager.GetKmerProfile(chr, start, end, kmerSize)
		if err != nil && verbose {
			fmt.Printf("Warning: Failed to get cached k-mer profile for %s:%d-%d: %v\n", chr, start, end, err)
		}
	}
	
	cacheDuration := time.Since(cacheStart)
	fmt.Printf("   Cached k-mer access time: %v\n", cacheDuration)
	
	// Performance comparison
	if kmerDuration > 0 && cacheDuration > 0 {
		cacheSpeedup := (kmerDuration.Seconds() / float64(numRegions)) / 
			(cacheDuration.Seconds() / float64(numRegions/2))
		fmt.Printf("   Cache speedup: %.2fx\n", cacheSpeedup)
	}
	
	fmt.Printf("\n")
}

func runComprehensiveBenchmark(manager *reference.ReferenceManager, numRegions int, regionSize int64, 
	kmerSize int, warmupRegions int, verbose bool) {
	
	fmt.Printf("🎯 Running Comprehensive Reference Genome Benchmark\n")
	fmt.Printf("===================================================\n\n")

	// Memory usage before benchmarks
	var memBefore runtime.MemStats
	runtime.GC()
	runtime.ReadMemStats(&memBefore)

	// Run all benchmark phases
	runCacheBenchmark(manager, numRegions, regionSize, warmupRegions, verbose)
	runPrefetchBenchmark(manager, numRegions, regionSize, verbose)
	runKmerBenchmark(manager, numRegions, regionSize, kmerSize, verbose)

	// Memory usage after benchmarks
	var memAfter runtime.MemStats
	runtime.GC()
	runtime.ReadMemStats(&memAfter)

	// Display comprehensive results
	fmt.Printf("📊 Comprehensive Results Summary\n")
	fmt.Printf("================================\n")
	
	stats := manager.GetCombinedStatistics()
	fmt.Printf("Overall cache hit ratio: %.1f%%\n", stats.OverallHitRatio*100)
	fmt.Printf("Total regions cached: %d\n", stats.TotalCachedRegions)
	fmt.Printf("Cache utilization: %.1f MB / %.1f MB (%.1f%%)\n",
		float64(stats.TotalCacheSize)/(1024*1024),
		float64(stats.TotalCacheSize + (2*1024*1024*1024 - stats.TotalCacheSize))/(1024*1024),
		float64(stats.TotalCacheSize)/float64(2*1024*1024*1024)*100)
	fmt.Printf("Average load time: %.2f ms\n", stats.AvgLoadTime/1e6)

	// Memory efficiency analysis
	memUsed := memAfter.Alloc - memBefore.Alloc
	fmt.Printf("Memory overhead: %.2f MB\n", float64(memUsed)/(1024*1024))
	
	// Performance recommendations
	fmt.Printf("\n💡 Performance Recommendations:\n")
	if stats.OverallHitRatio < 0.8 {
		fmt.Printf("   • Consider increasing cache size for better hit ratio\n")
	}
	if stats.AvgLoadTime > 10e6 { // 10ms
		fmt.Printf("   • Average load time is high - consider SSD storage or prefetching\n")
	}
	if float64(stats.TotalCacheSize)/float64(2*1024*1024*1024) > 0.9 {
		fmt.Printf("   • Cache utilization is high - consider increasing cache size\n")
	}
	
	fmt.Printf("\n")
}

func displayFinalStatistics(manager *reference.ReferenceManager) {
	fmt.Printf("📈 Final Reference Manager Statistics\n")
	fmt.Printf("====================================\n")

	stats := manager.GetCombinedStatistics()
	
	fmt.Printf("Total requests: %d\n", stats.TotalRequests)
	fmt.Printf("Cache hits: %d\n", stats.TotalCacheHits)
	fmt.Printf("Cache misses: %d\n", stats.TotalCacheMisses)
	fmt.Printf("Overall hit ratio: %.1f%%\n", stats.OverallHitRatio*100)
	fmt.Printf("Average request time: %.2f ms\n", stats.AvgLoadTime/1e6)

	// Display per-loader statistics
	if len(stats.LoaderStats) > 0 {
		fmt.Printf("\nPer-loader statistics:\n")
		for name, loaderStats := range stats.LoaderStats {
			fmt.Printf("  %s:\n", name)
			fmt.Printf("    Cache hits: %d\n", loaderStats.CacheHits)
			fmt.Printf("    Cache misses: %d\n", loaderStats.CacheMisses)
			fmt.Printf("    Hit ratio: %.1f%%\n", loaderStats.HitRatio*100)
			fmt.Printf("    Cache utilization: %.1f%%\n", loaderStats.CacheUtilization*100)
			fmt.Printf("    Cached regions: %d\n", loaderStats.CachedRegions)
		}
	}

	// Reference genome information
	refInfo := manager.GetReferenceInfo()
	fmt.Printf("\nReference genome info:\n")
	fmt.Printf("  Default reference: %t\n", refInfo.HasDefault)
	fmt.Printf("  Custom references: %d\n", len(refInfo.CustomReferences))
	fmt.Printf("  Total loaders: %d / %d\n", refInfo.TotalLoaders, refInfo.MaxLoaders)

	fmt.Printf("\n🎉 Benchmark completed successfully!\n")
}