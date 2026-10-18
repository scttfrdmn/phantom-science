package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/vcf"
)

func main() {
	var (
		inputFile      = flag.String("input", "", "Input VCF file (required)")
		outputFile     = flag.String("output", "", "Output analysis file (optional)")
		kmerSize       = flag.Int("k", 8, "K-mer size for authenticity analysis")
		bufferSize     = flag.Int("buffer", 65536, "Buffer size for streaming (bytes)")
		maxCache       = flag.Int64("cache", 1073741824, "Maximum cache size in bytes (1GB default)")
		numWorkers     = flag.Int("workers", 4, "Number of worker goroutines")
		useGPU         = flag.Bool("gpu", false, "Enable GPU acceleration")
		useML          = flag.Bool("ml", false, "Enable ML-based authenticity validation")
		targetScore    = flag.Float64("target", 50.0, "Reference value for the mean k-mer diversity score (50 = half of all possible k-mers observed per context window); a heuristic, not a classifier accuracy")
		tolerance      = flag.Float64("tolerance", 2.1, "Allowed absolute deviation of the mean score from -target")
		verbose        = flag.Bool("verbose", false, "Enable verbose output")
		showProgress   = flag.Bool("progress", false, "Show processing progress")
		contextWindow  = flag.Int("context", 1000, "Genomic context window size")
		populationPanels = flag.String("populations", "", "Comma-separated population panel list")
	)
	flag.Parse()

	fmt.Printf("🧬 Phantom Science VCF Streaming Authenticity Analyzer\n")
	fmt.Printf("======================================================\n\n")

	// Validate required parameters
	if *inputFile == "" {
		log.Fatal("Error: Input VCF file is required (-input)")
	}

	// Check if input file exists
	if _, err := os.Stat(*inputFile); os.IsNotExist(err) {
		log.Fatalf("Error: Input file does not exist: %s", *inputFile)
	}

	// Display configuration
	fmt.Printf("📁 Input VCF: %s\n", *inputFile)
	if *outputFile != "" {
		fmt.Printf("📄 Output file: %s\n", *outputFile)
	}
	fmt.Printf("🧬 K-mer size: %d\n", *kmerSize)
	fmt.Printf("💾 Buffer size: %d bytes (%.1f KB)\n", *bufferSize, float64(*bufferSize)/1024)
	fmt.Printf("🗄️  Cache size: %d bytes (%.1f MB)\n", *maxCache, float64(*maxCache)/(1024*1024))
	fmt.Printf("⚙️  Workers: %d\n", *numWorkers)
	fmt.Printf("🎯 Score reference range: %.1f ± %.1f\n", *targetScore, *tolerance)
	
	if *useGPU {
		fmt.Printf("🚀 GPU acceleration: enabled\n")
	}
	if *useML {
		fmt.Printf("🤖 ML validation: enabled\n")
	}
	if *populationPanels != "" {
		fmt.Printf("👥 Population panels: %s\n", *populationPanels)
	}
	fmt.Printf("\n")

	// Create stream processor configuration
	streamConfig := vcf.StreamConfig{
		BufferSize:        *bufferSize,
		MaxCacheSize:      *maxCache,
		EnableCompression: true,
		PackedArrays:      true, // Always use memory-efficient packed arrays
		GPUAcceleration:   *useGPU,
		NumWorkers:        *numWorkers,
		ContextWindowSize: *contextWindow,
		KmerSize:          *kmerSize,
		CacheRegionSize:   1000000, // 1MB regions
	}

	// Create stream processor
	streamProcessor, err := vcf.NewStreamProcessor(streamConfig)
	if err != nil {
		log.Fatalf("Failed to create stream processor: %v", err)
	}

	// Parse population panels
	var populationList []string
	if *populationPanels != "" {
		populationList = strings.Split(*populationPanels, ",")
		for i, pop := range populationList {
			populationList[i] = strings.TrimSpace(pop)
		}
	}

	// Create authenticity processor configuration
	authConfig := vcf.AuthenticityConfig{
		KmerSize:          *kmerSize,
		ContextWindowSize: *contextWindow,
		UseGPU:            *useGPU,
		UseML:             *useML,
		PopulationPanels:  populationList,
		TargetScore:       *targetScore,
		EnableCaching:     true,
		MaxCacheSize:      *maxCache,
	}

	// Create authenticity processor
	authProcessor, err := vcf.NewAuthenticityProcessor(authConfig)
	if err != nil {
		log.Fatalf("Failed to create authenticity processor: %v", err)
	}
	defer authProcessor.Close()

	// Process VCF file
	fmt.Printf("🔄 Processing VCF file...\n")
	startTime := time.Now()

	// Start progress monitoring if requested
	var progressTicker *time.Ticker
	if *showProgress {
		progressTicker = time.NewTicker(5 * time.Second)
		defer progressTicker.Stop()
		
		go func() {
			for range progressTicker.C {
				stats := streamProcessor.GetStatistics()
				fmt.Printf("📊 Progress: %d variants, %.2f MB processed, %.1f%% cache hit ratio\n",
					stats.ProcessedVariants,
					float64(stats.ProcessedBytes)/(1024*1024),
					stats.CacheHitRatio*100)
			}
		}()
	}

	// Process the VCF file
	err = streamProcessor.ProcessVCFFile(*inputFile, authProcessor)
	if err != nil {
		log.Fatalf("Failed to process VCF file: %v", err)
	}

	processingTime := time.Since(startTime)

	// Stop progress monitoring
	if progressTicker != nil {
		progressTicker.Stop()
	}

	fmt.Printf("✅ Processing completed in %v\n\n", processingTime)

	// Get processing statistics
	streamStats := streamProcessor.GetStatistics()
	authResults := authProcessor.GetAuthenticityResults()

	// Compare against the reference range before display so the summary
	// reports the actual outcome.
	withinRange := authResults.IsAuthenticityTargetMet(*targetScore, *tolerance)

	// Display results
	displayResults(streamStats, authResults, processingTime, *targetScore, *tolerance, *verbose)

	// Write output file if specified
	if *outputFile != "" {
		if err := writeOutputFile(*outputFile, streamStats, authResults); err != nil {
			log.Printf("Warning: Failed to write output file: %v", err)
		} else {
			fmt.Printf("📄 Results written to: %s\n", *outputFile)
		}
	}

	// Report the reference-range check. The score is a k-mer diversity
	// heuristic; it does not measure distinguishability from real data.
	if withinRange {
		fmt.Printf("\n✅ Mean k-mer diversity score within reference range\n")
		fmt.Printf("   Reference: %.1f ± %.1f\n", *targetScore, *tolerance)
		fmt.Printf("   Achieved: %.1f (deviation: %.1f)\n", 
			authResults.MeanAuthenticityScore, authResults.TargetScoreDeviation)
		os.Exit(0)
	} else {
		fmt.Printf("\n⚠️  WARNING: Mean k-mer diversity score outside reference range\n")
		fmt.Printf("   Reference: %.1f ± %.1f\n", *targetScore, *tolerance)
		fmt.Printf("   Achieved: %.1f (deviation: %.1f)\n", 
			authResults.MeanAuthenticityScore, authResults.TargetScoreDeviation)
		os.Exit(1)
	}
}

func displayResults(streamStats vcf.ProcessingStatistics, authResults *vcf.AuthenticityResults, 
					processingTime time.Duration, targetScore, tolerance float64, verbose bool) {
	
	fmt.Printf("📊 Processing Results\n")
	fmt.Printf("====================\n")
	
	// Stream processing stats
	fmt.Printf("🔄 Stream Processing:\n")
	fmt.Printf("   Variants processed: %d\n", streamStats.ProcessedVariants)
	fmt.Printf("   Data processed: %.2f MB\n", float64(streamStats.ProcessedBytes)/(1024*1024))
	fmt.Printf("   Processing time: %v\n", processingTime)
	fmt.Printf("   Throughput: %.2f variants/sec\n", 
		float64(streamStats.ProcessedVariants)/processingTime.Seconds())
	fmt.Printf("   Cache hit ratio: %.1f%%\n", streamStats.CacheHitRatio*100)
	fmt.Printf("   Cache utilization: %.1f%%\n", streamStats.CacheUtilization*100)

	fmt.Printf("\n🎯 Authenticity Analysis:\n")
	fmt.Printf("   Mean authenticity score: %.2f\n", authResults.MeanAuthenticityScore)
	fmt.Printf("   Median authenticity score: %.2f\n", authResults.MedianScore)
	fmt.Printf("   Standard deviation: %.2f\n", authResults.StandardDeviation)
	fmt.Printf("   Reference range: %.1f ± %.1f\n", targetScore, tolerance)

	if authResults.IsStatisticallyIndistinguishable {
		fmt.Printf("   ✅ Within reference range\n")
	} else {
		fmt.Printf("   ❌ Outside reference range\n")
	}
	fmt.Printf("   Deviation from reference: %.2f\n", authResults.TargetScoreDeviation)

	fmt.Printf("\n🧬 K-mer Analysis:\n")
	fmt.Printf("   Total k-mers: %d\n", authResults.TotalKmers)
	fmt.Printf("   Unique k-mers: %d\n", authResults.UniqueKmers)
	fmt.Printf("   K-mer complexity: %.4f\n", authResults.KmerComplexity)

	fmt.Printf("\n💾 Memory Efficiency:\n")
	fmt.Printf("   Memory usage: %.2f MB\n", float64(authResults.MemoryUsage)/(1024*1024))
	fmt.Printf("   Compression ratio: %.2fx\n", authResults.CompressionRatio)

	// Population-specific results
	if len(authResults.PopulationResults) > 1 {
		fmt.Printf("\n👥 Population Analysis:\n")
		for popID, popResult := range authResults.PopulationResults {
			fmt.Printf("   %s:\n", popID)
			fmt.Printf("     Samples: %d\n", popResult.SampleCount)
			fmt.Printf("     Variants: %d\n", popResult.VariantCount)
			fmt.Printf("     Mean authenticity: %.2f\n", popResult.MeanAuthenticityScore)
			fmt.Printf("     K-mer diversity: %.4f\n", popResult.KmerDiversity)
			if verbose {
				fmt.Printf("     Authenticity variance: %.4f\n", popResult.AuthenticityVariance)
				fmt.Printf("     Memory usage: %.2f KB\n", float64(popResult.MemoryUsage)/1024)
			}
		}
	}

	if verbose {
		fmt.Printf("\n🔍 Detailed Statistics:\n")
		fmt.Printf("   Cache hits: %d\n", streamStats.CacheHits)
		fmt.Printf("   Cache misses: %d\n", streamStats.CacheMisses)
		fmt.Printf("   Cached regions: %d\n", streamStats.CachedRegions)
		fmt.Printf("   Cache size: %.2f MB\n", float64(streamStats.CacheSize)/(1024*1024))
		
		if len(authResults.AuthenticityScores) > 0 {
			fmt.Printf("   Score range: %.2f - %.2f\n", 
				findMin(authResults.AuthenticityScores),
				findMax(authResults.AuthenticityScores))
		}
	}
}

func writeOutputFile(filename string, streamStats vcf.ProcessingStatistics, authResults *vcf.AuthenticityResults) error {
	file, err := os.Create(filename) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()

	// Write JSON-formatted results
	fmt.Fprintf(file, "{\n")
	fmt.Fprintf(file, "  \"processing_statistics\": {\n")
	fmt.Fprintf(file, "    \"processed_variants\": %d,\n", streamStats.ProcessedVariants)
	fmt.Fprintf(file, "    \"processed_bytes\": %d,\n", streamStats.ProcessedBytes)
	fmt.Fprintf(file, "    \"cache_hit_ratio\": %.4f,\n", streamStats.CacheHitRatio)
	fmt.Fprintf(file, "    \"cache_utilization\": %.4f\n", streamStats.CacheUtilization)
	fmt.Fprintf(file, "  },\n")
	
	fmt.Fprintf(file, "  \"authenticity_results\": {\n")
	fmt.Fprintf(file, "    \"processed_variants\": %d,\n", authResults.ProcessedVariants)
	fmt.Fprintf(file, "    \"mean_authenticity_score\": %.4f,\n", authResults.MeanAuthenticityScore)
	fmt.Fprintf(file, "    \"median_score\": %.4f,\n", authResults.MedianScore)
	fmt.Fprintf(file, "    \"standard_deviation\": %.4f,\n", authResults.StandardDeviation)
	fmt.Fprintf(file, "    \"total_kmers\": %d,\n", authResults.TotalKmers)
	fmt.Fprintf(file, "    \"unique_kmers\": %d,\n", authResults.UniqueKmers)
	fmt.Fprintf(file, "    \"kmer_complexity\": %.4f,\n", authResults.KmerComplexity)
	fmt.Fprintf(file, "    \"memory_usage\": %d,\n", authResults.MemoryUsage)
	fmt.Fprintf(file, "    \"compression_ratio\": %.4f,\n", authResults.CompressionRatio)
	fmt.Fprintf(file, "    \"statistically_indistinguishable\": %t,\n", authResults.IsStatisticallyIndistinguishable)
	fmt.Fprintf(file, "    \"target_score_deviation\": %.4f\n", authResults.TargetScoreDeviation)
	fmt.Fprintf(file, "  }\n")
	fmt.Fprintf(file, "}\n")

	return nil
}

func findMin(scores []float64) float64 {
	if len(scores) == 0 {
		return 0
	}
	min := scores[0]
	for _, score := range scores[1:] {
		if score < min {
			min = score
		}
	}
	return min
}

func findMax(scores []float64) float64 {
	if len(scores) == 0 {
		return 0
	}
	max := scores[0]
	for _, score := range scores[1:] {
		if score > max {
			max = score
		}
	}
	return max
}

