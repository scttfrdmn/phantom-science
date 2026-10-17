package main

import (
	"flag"
	"fmt"
	"log"
	"runtime"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/acceleration/metal"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

func main() {
	var (
		k              = flag.Int("k", 6, "K-mer size (1-16)")
		numSequences   = flag.Int("sequences", 1000, "Number of test sequences")
		sequenceLength = flag.Int("length", 1000, "Length of each test sequence")
		verbose        = flag.Bool("verbose", false, "Enable verbose output")
		cpuOnly        = flag.Bool("cpu-only", false, "Run CPU benchmark only")
		gpuOnly        = flag.Bool("gpu-only", false, "Run GPU benchmark only")
	)
	flag.Parse()

	fmt.Printf("🧬 Phantom Science GPU K-mer Processing Benchmark\n")
	fmt.Printf("================================================\n\n")

	// Validate parameters
	if *k < 1 || *k > 16 {
		log.Fatalf("Error: k must be between 1 and 16, got %d", *k)
	}

	if *numSequences < 1 || *sequenceLength < *k {
		log.Fatalf("Error: invalid sequence parameters (sequences=%d, length=%d, k=%d)", 
			*numSequences, *sequenceLength, *k)
	}

	// Check platform
	if runtime.GOOS != "darwin" {
		fmt.Printf("❌ Metal GPU acceleration is only available on macOS\n")
		fmt.Printf("🔄 Running CPU-only benchmark on %s\n\n", runtime.GOOS)
		*cpuOnly = true
	}

	// Check Metal availability
	var metalAvailable bool
	if !*cpuOnly {
		metalAvailable = metal.IsAvailable()
		if !metalAvailable {
			fmt.Printf("❌ Metal GPU not available on this system\n")
			fmt.Printf("🔄 Running CPU-only benchmark\n\n")
			*cpuOnly = true
		} else {
			fmt.Printf("✅ Metal GPU acceleration available\n\n")
		}
	}

	// Generate test data
	fmt.Printf("🔧 Generating test data...\n")
	fmt.Printf("   Sequences: %d\n", *numSequences)
	fmt.Printf("   Length: %d bp each\n", *sequenceLength)
	fmt.Printf("   K-mer size: %d\n", *k)
	fmt.Printf("   Total data: %.2f MB\n\n", float64(*numSequences * *sequenceLength) / (1024 * 1024))

	testSequences := generateTestSequences(*numSequences, *sequenceLength)
	if *verbose {
		fmt.Printf("✅ Test data generated: %d sequences\n\n", len(testSequences))
	}

	// Run CPU benchmark
	if !*gpuOnly {
		fmt.Printf("🖥️  CPU Benchmark\n")
		fmt.Printf("================\n")
		runCPUBenchmark(testSequences, *k, *verbose)
		fmt.Printf("\n")
	}

	// Run GPU benchmark
	if !*cpuOnly && metalAvailable {
		fmt.Printf("🚀 GPU Benchmark (Apple Metal)\n")
		fmt.Printf("==============================\n")
		runGPUBenchmark(testSequences, *k, *verbose)
		fmt.Printf("\n")
	}

	// Run comparison if both are available
	if !*cpuOnly && !*gpuOnly && metalAvailable {
		fmt.Printf("📊 Performance Comparison\n")
		fmt.Printf("=========================\n")
		runPerformanceComparison(testSequences, *k, *verbose)
	}
}

func generateTestSequences(count, length int) [][]byte {
	sequences := make([][]byte, count)
	bases := []byte{'A', 'T', 'G', 'C'}
	
	for i := 0; i < count; i++ {
		seq := make([]byte, length)
		for j := 0; j < length; j++ {
			seq[j] = bases[(i*length+j)%4] // Deterministic but varied pattern
		}
		sequences[i] = seq
	}
	
	return sequences
}

func runCPUBenchmark(sequences [][]byte, k int, verbose bool) {
	extractor := kmer.NewFrequencyExtractor(k)
	
	start := time.Now()
	result, err := extractor.ExtractFromMultipleSequences(sequences)
	duration := time.Since(start)
	
	if err != nil {
		fmt.Printf("❌ CPU benchmark failed: %v\n", err)
		return
	}

	// Calculate throughput
	totalBases := 0
	for _, seq := range sequences {
		totalBases += len(seq)
	}
	
	throughputMBps := float64(totalBases) / (1024 * 1024) / duration.Seconds()
	kmersPerSecond := float64(result.TotalCount) / duration.Seconds()

	fmt.Printf("⏱️  Duration: %v\n", duration)
	fmt.Printf("📈 Throughput: %.2f MB/s\n", throughputMBps)
	fmt.Printf("🧬 K-mers/sec: %.0f\n", kmersPerSecond)
	fmt.Printf("📊 Total k-mers: %d\n", result.TotalCount)
	fmt.Printf("🎯 Unique k-mers: %d\n", len(result.Frequencies))

	if verbose {
		stats := result.GetKmerStats()
		fmt.Printf("   Complexity: %.4f\n", stats.Complexity)
		fmt.Printf("   Shannon entropy: %.4f\n", stats.ShannonEntropy)
		fmt.Printf("   GC content: %.4f\n", stats.GCContent)
	}
}

func runGPUBenchmark(sequences [][]byte, k int, verbose bool) {
	processor, err := metal.NewMetalKmerProcessor()
	if err != nil {
		fmt.Printf("❌ Failed to initialize Metal processor: %v\n", err)
		return
	}
	defer processor.Close()

	if verbose {
		deviceInfo := processor.GetDeviceInfo()
		fmt.Printf("🔧 Device: %s\n", deviceInfo.Name)
		fmt.Printf("   Max threads: %d\n", deviceInfo.MaxThreads)
		fmt.Printf("   Max buffer: %.2f MB\n", float64(deviceInfo.MaxBufferSize)/(1024*1024))
	}

	start := time.Now()
	result, err := processor.ExtractKmerFrequenciesGPU(sequences, k)
	duration := time.Since(start)

	if err != nil {
		fmt.Printf("❌ GPU benchmark failed: %v\n", err)
		return
	}

	// Calculate throughput
	totalBases := 0
	for _, seq := range sequences {
		totalBases += len(seq)
	}
	
	throughputMBps := float64(totalBases) / (1024 * 1024) / duration.Seconds()
	kmersPerSecond := float64(result.TotalCount) / duration.Seconds()

	fmt.Printf("⏱️  Duration: %v\n", duration)
	fmt.Printf("📈 Throughput: %.2f MB/s\n", throughputMBps)
	fmt.Printf("🧬 K-mers/sec: %.0f\n", kmersPerSecond)
	fmt.Printf("📊 Total k-mers: %d\n", result.TotalCount)
	fmt.Printf("🎯 Unique k-mers: %d\n", len(result.Frequencies))

	if verbose {
		stats := result.GetKmerStats()
		fmt.Printf("   Complexity: %.4f\n", stats.Complexity)
		fmt.Printf("   Shannon entropy: %.4f\n", stats.ShannonEntropy)
		fmt.Printf("   GC content: %.4f\n", stats.GCContent)
	}
}

func runPerformanceComparison(sequences [][]byte, k int, verbose bool) {
	// Run comprehensive benchmark comparison
	processor, err := metal.NewMetalKmerProcessor()
	if err != nil {
		fmt.Printf("❌ Failed to initialize Metal processor: %v\n", err)
		return
	}
	defer processor.Close()

	benchmarkResults, err := processor.BenchmarkGPUPerformance(sequences, k)
	if err != nil {
		fmt.Printf("❌ Benchmark comparison failed: %v\n", err)
		return
	}

	fmt.Printf("🏁 Results Summary:\n")
	fmt.Printf("   CPU Duration: %v\n", benchmarkResults.CPUDuration)
	fmt.Printf("   GPU Duration: %v\n", benchmarkResults.GPUDuration)
	fmt.Printf("   🚀 Speedup: %.2fx\n", benchmarkResults.Speedup)
	fmt.Printf("   📈 GPU Throughput: %.2f MB/s\n", benchmarkResults.Throughput)
	fmt.Printf("   ✅ Results Match: %t\n", benchmarkResults.ResultsMatch)

	if benchmarkResults.Speedup > 1.0 {
		fmt.Printf("   🎉 GPU acceleration achieved %.2fx performance improvement!\n", benchmarkResults.Speedup)
	} else {
		fmt.Printf("   ⚠️  GPU slower than CPU - consider larger datasets or different k-mer sizes\n")
	}

	if verbose {
		fmt.Printf("\n📋 Detailed Results:\n")
		fmt.Printf("   Device: %s\n", benchmarkResults.DeviceName)
		fmt.Printf("   Test sequences: %d\n", benchmarkResults.TestSequences)
		fmt.Printf("   K-mer size: %d\n", benchmarkResults.KmerSize)
		fmt.Printf("   GPU k-mers: %d\n", benchmarkResults.GPUKmers)
		fmt.Printf("   CPU k-mers: %d\n", benchmarkResults.CPUKmers)
	}
}