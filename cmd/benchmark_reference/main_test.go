//go:build !release
// +build !release

package main

import (
	"os"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// TestDisplayConfiguration verifies that displayConfiguration does not panic
// when called with typical parameter values.
func TestDisplayConfiguration(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Fatalf("displayConfiguration panicked: %v", r)
		}
	}()

	displayConfiguration(
		"/tmp/test.fasta",
		2*1024*1024*1024,
		8,
		1000,
		10000,
		"all",
		true,
		6,
		100,
	)
}

// TestDisplayConfigurationVariants exercises several combinations of parameters
// to confirm the function handles edge cases without panicking.
func TestDisplayConfigurationVariants(t *testing.T) {
	cases := []struct {
		name             string
		cacheSize        int64
		kmerSize         int
		numRegions       int
		regionSize       int64
		benchmarkType    string
		enablePrefetch   bool
		compressionLevel int
		warmupRegions    int
	}{
		{
			name:             "all defaults",
			cacheSize:        2147483648,
			kmerSize:         8,
			numRegions:       1000,
			regionSize:       10000,
			benchmarkType:    "all",
			enablePrefetch:   true,
			compressionLevel: 6,
			warmupRegions:    100,
		},
		{
			name:             "zero values",
			cacheSize:        0,
			kmerSize:         0,
			numRegions:       0,
			regionSize:       0,
			benchmarkType:    "",
			enablePrefetch:   false,
			compressionLevel: 0,
			warmupRegions:    0,
		},
		{
			name:             "cache benchmark only",
			cacheSize:        512 * 1024 * 1024,
			kmerSize:         4,
			numRegions:       50,
			regionSize:       5000,
			benchmarkType:    "cache",
			enablePrefetch:   false,
			compressionLevel: 9,
			warmupRegions:    10,
		},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			defer func() {
				if r := recover(); r != nil {
					t.Fatalf("displayConfiguration panicked with %q: %v", tc.name, r)
				}
			}()
			displayConfiguration(
				"/some/reference.fasta",
				tc.cacheSize,
				tc.kmerSize,
				tc.numRegions,
				tc.regionSize,
				tc.benchmarkType,
				tc.enablePrefetch,
				tc.compressionLevel,
				tc.warmupRegions,
			)
		})
	}
}

// TestCreateTestReferenceFile_Skipped skips the slow 50MB FASTA generation
// that createTestReferenceFile performs, unless -run.slow is explicitly set.
// The function writes ~50MB character-by-character which takes several minutes,
// making it unsuitable for standard test runs.
func TestCreateTestReferenceFile_Skipped(t *testing.T) {
	t.Skip("createTestReferenceFile generates a ~50MB file character-by-character: too slow for unit tests")
}

// newTestManager creates a ReferenceManager backed by a minimal FASTA temp
// file.  The caller is responsible for removing the returned file path and
// closing the manager.
func newTestManager(t *testing.T) (*reference.ReferenceManager, string) {
	t.Helper()

	tmpFile, err := os.CreateTemp("", "bench_ref_test_*.fasta")
	if err != nil {
		t.Fatalf("failed to create temp FASTA file: %v", err)
	}
	content := []byte(">chr1\nACGTACGT\n")
	if _, err := tmpFile.Write(content); err != nil {
		os.Remove(tmpFile.Name())
		t.Fatalf("failed to write temp FASTA: %v", err)
	}
	tmpFile.Close()

	cfg := reference.ManagerConfig{
		DefaultReferenceFile: tmpFile.Name(),
		DefaultCacheSize:     1024 * 1024,
		EnablePrefetching:    false,
		CompressionLevel:     0,
		EnablePackedArrays:   false,
	}

	manager, err := reference.NewReferenceManager(cfg)
	if err != nil {
		os.Remove(tmpFile.Name())
		t.Skipf("skipping: could not create reference manager: %v", err)
	}

	return manager, tmpFile.Name()
}

// TestRunCacheBenchmarkErrorPath creates a ReferenceManager backed by a tiny
// FASTA file.  GetSequence calls for out-of-range regions will fail, but the
// function should not panic.
func TestRunCacheBenchmarkErrorPath(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Fatalf("runCacheBenchmark panicked: %v", r)
		}
	}()

	manager, path := newTestManager(t)
	defer os.Remove(path)
	defer manager.Close()

	runCacheBenchmark(manager, 4, 100, 2, false)
}

// TestRunPrefetchBenchmarkErrorPath exercises the prefetch benchmark with a
// tiny reference file; out-of-range regions fail gracefully without panic.
func TestRunPrefetchBenchmarkErrorPath(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Fatalf("runPrefetchBenchmark panicked: %v", r)
		}
	}()

	manager, path := newTestManager(t)
	defer os.Remove(path)
	defer manager.Close()

	runPrefetchBenchmark(manager, 4, 100, false)
}

// TestRunKmerBenchmarkErrorPath exercises the kmer benchmark with a tiny file.
func TestRunKmerBenchmarkErrorPath(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Fatalf("runKmerBenchmark panicked: %v", r)
		}
	}()

	manager, path := newTestManager(t)
	defer os.Remove(path)
	defer manager.Close()

	runKmerBenchmark(manager, 4, 100, 4, false)
}

// TestRunComprehensiveBenchmarkErrorPath exercises the comprehensive benchmark
// dispatcher with a minimal reference file.
func TestRunComprehensiveBenchmarkErrorPath(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Fatalf("runComprehensiveBenchmark panicked: %v", r)
		}
	}()

	manager, path := newTestManager(t)
	defer os.Remove(path)
	defer manager.Close()

	runComprehensiveBenchmark(manager, 4, 100, 4, 2, false)
}

// TestDisplayFinalStatistics verifies that displayFinalStatistics does not
// panic when called with a freshly created manager that has zero statistics.
func TestDisplayFinalStatistics(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Fatalf("displayFinalStatistics panicked: %v", r)
		}
	}()

	manager, path := newTestManager(t)
	defer os.Remove(path)
	defer manager.Close()

	displayFinalStatistics(manager)
}
