package main

import (
	"flag"
	"fmt"
	"log"
	"runtime"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

func main() {
	var (
		k                 = flag.Int("k", 8, "K-mer size (1-16)")
		numSequences      = flag.Int("sequences", 1000, "Number of test sequences")
		sequenceLength    = flag.Int("length", 1000, "Length of each test sequence")
		verbose           = flag.Bool("verbose", false, "Enable verbose output")
		packedOnly        = flag.Bool("packed-only", false, "Test packed arrays only")
		mapOnly           = flag.Bool("map-only", false, "Test hash maps only")
		benchmarkMemory   = flag.Bool("memory", false, "Run memory usage benchmarks")
		compareAccuracy   = flag.Bool("accuracy", false, "Compare accuracy between methods")
	)
	flag.Parse()

	fmt.Printf("🧬 Phantom Science Memory Optimization Benchmark\n")
	fmt.Printf("================================================\n\n")

	// Validate parameters
	if *k < 1 || *k > 16 {
		log.Fatalf("Error: k must be between 1 and 16, got %d", *k)
	}

	if *numSequences < 1 || *sequenceLength < *k {
		log.Fatalf("Error: invalid sequence parameters (sequences=%d, length=%d, k=%d)",
			*numSequences, *sequenceLength, *k)
	}

	// Generate test data
	fmt.Printf("🔧 Generating test data...\n")
	fmt.Printf("   Sequences: %d\n", *numSequences)
	fmt.Printf("   Length: %d bp each\n", *sequenceLength)
	fmt.Printf("   K-mer size: %d\n", *k)
	fmt.Printf("   Total data: %.2f MB\n\n", float64(*numSequences**sequenceLength)/(1024*1024))

	testSequences := generateTestSequences(*numSequences, *sequenceLength)
	if *verbose {
		fmt.Printf("✅ Test data generated: %d sequences\n\n", len(testSequences))
	}

	// Run benchmarks based on flags
	if !*packedOnly {
		fmt.Printf("📊 Hash Map-Based K-mer Processing\n")
		fmt.Printf("==================================\n")
		runMapBenchmark(testSequences, *k, *verbose)
		fmt.Printf("\n")
	}

	if !*mapOnly {
		fmt.Printf("🗜️  Packed Array K-mer Processing\n")
		fmt.Printf("=================================\n")
		runPackedBenchmark(testSequences, *k, *verbose)
		fmt.Printf("\n")
	}

	if !*packedOnly && !*mapOnly {
		fmt.Printf("⚖️  Performance Comparison\n")
		fmt.Printf("==========================\n")
		runPerformanceComparison(testSequences, *k, *verbose)
		fmt.Printf("\n")
	}

	if *benchmarkMemory {
		fmt.Printf("🧠 Memory Usage Analysis\n")
		fmt.Printf("========================\n")
		runMemoryBenchmark(testSequences, *k, *verbose)
		fmt.Printf("\n")
	}

	if *compareAccuracy {
		fmt.Printf("🎯 Accuracy Verification\n")
		fmt.Printf("========================\n")
		runAccuracyComparison(testSequences, *k, *verbose)
		fmt.Printf("\n")
	}

	// Memory usage analysis
	if *verbose {
		fmt.Printf("💾 Runtime Memory Statistics\n")
		fmt.Printf("============================\n")
		var m runtime.MemStats
		runtime.GC()
		runtime.ReadMemStats(&m)
		fmt.Printf("Alloc: %.2f MB\n", float64(m.Alloc)/(1024*1024))
		fmt.Printf("TotalAlloc: %.2f MB\n", float64(m.TotalAlloc)/(1024*1024))
		fmt.Printf("Sys: %.2f MB\n", float64(m.Sys)/(1024*1024))
		fmt.Printf("NumGC: %d\n", m.NumGC)
	}
}

func generateTestSequences(count, length int) [][]byte {
	sequences := make([][]byte, count)
	bases := []byte{'A', 'T', 'G', 'C'}

	for i := 0; i < count; i++ {
		seq := make([]byte, length)
		for j := 0; j < length; j++ {
			// Create varied but deterministic sequences
			baseIndex := (i*length + j) % 4
			if j%10 == 0 {
				baseIndex = (baseIndex + 1) % 4 // Add some variation
			}
			seq[j] = bases[baseIndex]
		}
		sequences[i] = seq
	}

	return sequences
}

func runMapBenchmark(sequences [][]byte, k int, verbose bool) {
	extractor := kmer.NewFrequencyExtractor(k)

	// Measure memory before
	var m1 runtime.MemStats
	runtime.GC()
	runtime.ReadMemStats(&m1)

	start := time.Now()
	result, err := extractor.ExtractFromMultipleSequences(sequences)
	duration := time.Since(start)

	// Measure memory after
	var m2 runtime.MemStats
	runtime.ReadMemStats(&m2)

	if err != nil {
		fmt.Printf("❌ Hash map benchmark failed: %v\n", err)
		return
	}

	// Calculate throughput
	totalBases := 0
	for _, seq := range sequences {
		totalBases += len(seq)
	}

	throughputMBps := float64(totalBases) / (1024 * 1024) / duration.Seconds()
	kmersPerSecond := float64(result.TotalCount) / duration.Seconds()
	memoryUsed := int64(m2.Alloc - m1.Alloc) // #nosec G115

	fmt.Printf("⏱️  Duration: %v\n", duration)
	fmt.Printf("📈 Throughput: %.2f MB/s\n", throughputMBps)
	fmt.Printf("🧬 K-mers/sec: %.0f\n", kmersPerSecond)
	fmt.Printf("📊 Total k-mers: %d\n", result.TotalCount)
	fmt.Printf("🎯 Unique k-mers: %d\n", len(result.Frequencies))
	fmt.Printf("💾 Memory used: %.2f MB\n", float64(memoryUsed)/(1024*1024))

	if verbose {
		stats := result.GetKmerStats()
		fmt.Printf("   Complexity: %.4f\n", stats.Complexity)
		fmt.Printf("   Shannon entropy: %.4f\n", stats.ShannonEntropy)
		fmt.Printf("   GC content: %.4f\n", stats.GCContent)
		fmt.Printf("   Mean count: %.2f\n", stats.MeanCount)
		fmt.Printf("   Max count: %d\n", stats.MaxCount)
	}
}

func runPackedBenchmark(sequences [][]byte, k int, verbose bool) {
	extractor := kmer.NewPackedFrequencyExtractor(k)

	// Measure memory before
	var m1 runtime.MemStats
	runtime.GC()
	runtime.ReadMemStats(&m1)

	start := time.Now()
	result, err := extractor.ExtractFromMultipleSequences(sequences)
	duration := time.Since(start)

	// Measure memory after
	var m2 runtime.MemStats
	runtime.ReadMemStats(&m2)

	if err != nil {
		fmt.Printf("❌ Packed array benchmark failed: %v\n", err)
		return
	}

	// Calculate throughput
	totalBases := 0
	for _, seq := range sequences {
		totalBases += len(seq)
	}

	throughputMBps := float64(totalBases) / (1024 * 1024) / duration.Seconds()
	kmersPerSecond := float64(result.TotalCount) / duration.Seconds()
	memoryUsed := int64(m2.Alloc - m1.Alloc) // #nosec G115

	fmt.Printf("⏱️  Duration: %v\n", duration)
	fmt.Printf("📈 Throughput: %.2f MB/s\n", throughputMBps)
	fmt.Printf("🧬 K-mers/sec: %.0f\n", kmersPerSecond)
	fmt.Printf("📊 Total k-mers: %d\n", result.TotalCount)
	fmt.Printf("🎯 Unique k-mers: %d\n", result.NonZeroCount)
	fmt.Printf("💾 Memory used: %.2f MB\n", float64(memoryUsed)/(1024*1024))
	fmt.Printf("🗜️  Compression ratio: %.2fx\n", result.GetCompressionRatio())

	if verbose {
		stats := result.GetStatistics()
		fmt.Printf("   Packed memory: %.2f MB\n", float64(stats.MemoryUsage)/(1024*1024))
		fmt.Printf("   Bits per item: %d\n", stats.BitsPerItem)
		fmt.Printf("   Array size: %d uint32s\n", stats.PackedArraySize)
		fmt.Printf("   Max frequency: %d\n", stats.MaxFrequency)
		fmt.Printf("   Sparsity: %.4f\n", stats.Sparsity)
	}
}

func runPerformanceComparison(sequences [][]byte, k int, verbose bool) {
	// Run both benchmarks for comparison
	fmt.Printf("🏁 Running comparison benchmarks...\n")

	// Hash map benchmark
	mapExtractor := kmer.NewFrequencyExtractor(k)
	mapStart := time.Now()
	mapResult, mapErr := mapExtractor.ExtractFromMultipleSequences(sequences)
	mapDuration := time.Since(mapStart)

	// Packed array benchmark
	packedExtractor := kmer.NewPackedFrequencyExtractor(k)
	packedStart := time.Now()
	packedResult, packedErr := packedExtractor.ExtractFromMultipleSequences(sequences)
	packedDuration := time.Since(packedStart)

	if mapErr != nil || packedErr != nil {
		fmt.Printf("❌ Comparison failed: map error=%v, packed error=%v\n", mapErr, packedErr)
		return
	}

	// Calculate performance metrics
	speedRatio := float64(mapDuration.Nanoseconds()) / float64(packedDuration.Nanoseconds())
	memoryRatio := packedResult.GetCompressionRatio()

	// Verify correctness
	correctness := verifyResults(mapResult, packedResult, k)

	fmt.Printf("🏆 Results Summary:\n")
	fmt.Printf("   Hash Map Duration: %v\n", mapDuration)
	fmt.Printf("   Packed Array Duration: %v\n", packedDuration)
	if speedRatio > 1.0 {
		fmt.Printf("   🚀 Packed arrays %.2fx FASTER\n", speedRatio)
	} else {
		fmt.Printf("   ⚠️  Hash maps %.2fx faster\n", 1.0/speedRatio)
	}
	fmt.Printf("   💾 Memory compression: %.2fx\n", memoryRatio)
	fmt.Printf("   ✅ Results correct: %t\n", correctness)

	if verbose && correctness {
		fmt.Printf("\n📋 Detailed Comparison:\n")
		fmt.Printf("   Hash map k-mers: %d\n", mapResult.TotalCount)
		fmt.Printf("   Packed array k-mers: %d\n", packedResult.TotalCount)
		fmt.Printf("   Hash map unique: %d\n", len(mapResult.Frequencies))
		fmt.Printf("   Packed array unique: %d\n", packedResult.NonZeroCount)
		fmt.Printf("   Memory efficiency: %.1fx improvement\n", memoryRatio)
	}
}

func runMemoryBenchmark(sequences [][]byte, k int, verbose bool) {
	// Test different k-mer sizes to show memory scaling
	kSizes := []int{4, 6, 8, 10, 12}
	
	fmt.Printf("K-mer Size | Hash Map Memory | Packed Memory | Compression Ratio\n")
	fmt.Printf("-----------|-----------------|---------------|------------------\n")

	for _, testK := range kSizes {
		if testK > 16 {
			continue
		}

		// Hash map memory usage
		mapExtractor := kmer.NewFrequencyExtractor(testK)
		var m1, m2 runtime.MemStats
		runtime.GC()
		runtime.ReadMemStats(&m1)
		_, mapErr := mapExtractor.ExtractFromMultipleSequences(sequences)
		runtime.ReadMemStats(&m2)
		mapMemory := int64(m2.Alloc - m1.Alloc) // #nosec G115

		// Packed array memory usage
		packedExtractor := kmer.NewPackedFrequencyExtractor(testK)
		runtime.GC()
		runtime.ReadMemStats(&m1)
		_, packedErr := packedExtractor.ExtractFromMultipleSequences(sequences)
		runtime.ReadMemStats(&m2)
		packedMemory := int64(m2.Alloc - m1.Alloc) // #nosec G115

		if mapErr != nil || packedErr != nil {
			continue
		}

		compression := float64(mapMemory) / float64(packedMemory)
		
		fmt.Printf("%-10d | %-15.2f | %-13.2f | %.2fx\n",
			testK,
			float64(mapMemory)/(1024*1024),
			float64(packedMemory)/(1024*1024),
			compression)
	}
}

func runAccuracyComparison(sequences [][]byte, k int, verbose bool) {
	// Extract with both methods
	mapExtractor := kmer.NewFrequencyExtractor(k)
	mapResult, mapErr := mapExtractor.ExtractFromMultipleSequences(sequences)

	packedExtractor := kmer.NewPackedFrequencyExtractor(k)
	packedResult, packedErr := packedExtractor.ExtractFromMultipleSequences(sequences)

	if mapErr != nil || packedErr != nil {
		fmt.Printf("❌ Accuracy test failed: map error=%v, packed error=%v\n", mapErr, packedErr)
		return
	}

	// Detailed accuracy verification
	correctness := verifyResults(mapResult, packedResult, k)
	
	if correctness {
		fmt.Printf("✅ Perfect accuracy: All k-mer frequencies match exactly\n")
		fmt.Printf("   Total k-mers: %d (both methods)\n", mapResult.TotalCount)
		fmt.Printf("   Unique k-mers: %d (both methods)\n", len(mapResult.Frequencies))
		
		if verbose {
			// Show some example k-mers
			topKmers := mapResult.GetTopKmers(5)
			fmt.Printf("\n🔍 Top 5 k-mers (verification):\n")
			for i, kmerFreq := range topKmers {
				packedFreq, _ := packedResult.GetFrequency(kmerFreq.Kmer)
				fmt.Printf("   %d. %s: %d (hash map) = %d (packed)\n",
					i+1, kmerFreq.Kmer, kmerFreq.Count, packedFreq)
			}
		}
	} else {
		fmt.Printf("❌ Accuracy mismatch detected\n")
		// Could add detailed mismatch analysis here
	}
}

func verifyResults(mapResult *kmer.FrequencyTable, packedResult *kmer.PackedFrequencyTable, k int) bool {
	// Check basic statistics
	if mapResult.TotalCount != packedResult.TotalCount {
		return false
	}

	if len(mapResult.Frequencies) != packedResult.NonZeroCount {
		return false
	}

	// Check all frequencies match
	for kmerStr, mapFreq := range mapResult.Frequencies {
		packedFreq, err := packedResult.GetFrequency(kmerStr)
		if err != nil || mapFreq != packedFreq {
			return false
		}
	}

	// Check all packed frequencies exist in map
	for i := 0; i < packedResult.MaxKmers; i++ {
		packedFreq := packedResult.PackedFrequencies.Get(i)
		if packedFreq > 0 {
			kmerStr := kmer.IndexToKmerString(i, k)
			mapFreq, exists := mapResult.Frequencies[kmerStr]
			if !exists || mapFreq != int(packedFreq) {
				return false
			}
		}
	}

	return true
}