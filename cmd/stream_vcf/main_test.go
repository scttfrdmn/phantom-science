//go:build !release
// +build !release

package main

import (
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/vcf"
)

func TestFindMin_SingleElement(t *testing.T) {
	scores := []float64{5.0}
	if got := findMin(scores); got != 5.0 {
		t.Errorf("findMin([5.0]) = %f, want 5.0", got)
	}
}

func TestFindMin_MultipleElements(t *testing.T) {
	scores := []float64{3.0, 1.0, 4.0, 1.5, 9.0}
	if got := findMin(scores); got != 1.0 {
		t.Errorf("findMin(multiple) = %f, want 1.0", got)
	}
}

func TestFindMin_Empty(t *testing.T) {
	got := findMin([]float64{})
	if got != 0 {
		t.Errorf("findMin([]) = %f, want 0", got)
	}
}

func TestFindMax_SingleElement(t *testing.T) {
	scores := []float64{5.0}
	if got := findMax(scores); got != 5.0 {
		t.Errorf("findMax([5.0]) = %f, want 5.0", got)
	}
}

func TestFindMax_MultipleElements(t *testing.T) {
	scores := []float64{3.0, 1.0, 4.0, 1.5, 9.0}
	if got := findMax(scores); got != 9.0 {
		t.Errorf("findMax(multiple) = %f, want 9.0", got)
	}
}

func TestFindMax_Empty(t *testing.T) {
	got := findMax([]float64{})
	if got != 0 {
		t.Errorf("findMax([]) = %f, want 0", got)
	}
}

// ---------- displayResults ----------

func TestDisplayResults_NoCrash(t *testing.T) {
	streamStats := vcf.ProcessingStatistics{
		ProcessedVariants: 100,
		ProcessedBytes:    4096,
		CacheHitRatio:     0.8,
		CacheUtilization:  0.5,
	}
	authResults := &vcf.AuthenticityResults{
		MeanAuthenticityScore:            50.0,
		MedianScore:                      49.5,
		StandardDeviation:                2.1,
		IsStatisticallyIndistinguishable: true,
		TargetScoreDeviation:             0.5,
	}
	// Should not panic
	displayResults(streamStats, authResults, time.Second, 50.0, 2.1, false)
}

func TestDisplayResults_Verbose(t *testing.T) {
	streamStats := vcf.ProcessingStatistics{
		ProcessedVariants: 50,
		ProcessedBytes:    2048,
		CacheHitRatio:     0.6,
		CacheUtilization:  0.4,
		CacheHits:         30,
		CacheMisses:       20,
		CachedRegions:     5,
		CacheSize:         1024,
	}
	authResults := &vcf.AuthenticityResults{
		MeanAuthenticityScore: 48.0,
		AuthenticityScores:    []float64{45.0, 48.0, 51.0},
	}
	// Should not panic with verbose=true
	displayResults(streamStats, authResults, 2*time.Second, 50.0, 2.1, true)
}

func TestDisplayResults_IndistinguishableFalse(t *testing.T) {
	streamStats := vcf.ProcessingStatistics{ProcessedVariants: 10, ProcessedBytes: 512}
	authResults := &vcf.AuthenticityResults{
		IsStatisticallyIndistinguishable: false,
	}
	displayResults(streamStats, authResults, 100*time.Millisecond, 50.0, 2.0, false)
}

func TestDisplayResults_WithPopulationResults(t *testing.T) {
	streamStats := vcf.ProcessingStatistics{ProcessedVariants: 200}
	authResults := &vcf.AuthenticityResults{
		MeanAuthenticityScore: 50.0,
		PopulationResults: map[string]*vcf.PopulationAuthenticityResult{
			"YRI": {SampleCount: 5, VariantCount: 100, MeanAuthenticityScore: 49.0, KmerDiversity: 0.7},
			"CEU": {SampleCount: 5, VariantCount: 90, MeanAuthenticityScore: 51.0, KmerDiversity: 0.72},
		},
	}
	// Should not panic with multiple population results
	displayResults(streamStats, authResults, time.Second, 50.0, 2.0, false)
}

func TestDisplayResults_WithPopulationResults_Verbose(t *testing.T) {
	streamStats := vcf.ProcessingStatistics{ProcessedVariants: 200}
	authResults := &vcf.AuthenticityResults{
		MeanAuthenticityScore: 50.0,
		PopulationResults: map[string]*vcf.PopulationAuthenticityResult{
			"YRI": {SampleCount: 3, VariantCount: 80, MeanAuthenticityScore: 50.0, MemoryUsage: 4096},
			"CEU": {SampleCount: 3, VariantCount: 75, MeanAuthenticityScore: 50.5, MemoryUsage: 4096},
		},
		AuthenticityScores: []float64{48.0, 50.0, 52.0},
	}
	// verbose=true exercises the verbose block within population results
	displayResults(streamStats, authResults, time.Second, 50.0, 2.0, true)
}

// ---------- writeOutputFile ----------

func TestWriteOutputFile_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	outPath := filepath.Join(tmpDir, "results.json")
	streamStats := vcf.ProcessingStatistics{
		ProcessedVariants: 200,
		ProcessedBytes:    8192,
		CacheHitRatio:     0.75,
		CacheUtilization:  0.6,
	}
	authResults := &vcf.AuthenticityResults{
		ProcessedVariants:                200,
		MeanAuthenticityScore:            51.0,
		MedianScore:                      50.5,
		StandardDeviation:                1.5,
		TotalKmers:                       1000,
		UniqueKmers:                      800,
		KmerComplexity:                   0.8,
		MemoryUsage:                      1024,
		CompressionRatio:                 1.5,
		IsStatisticallyIndistinguishable: true,
		TargetScoreDeviation:             0.3,
	}
	err := writeOutputFile(outPath, streamStats, authResults)
	if err != nil {
		t.Fatalf("writeOutputFile error: %v", err)
	}
	info, err := os.Stat(outPath)
	if err != nil {
		t.Fatalf("output file not created: %v", err)
	}
	if info.Size() == 0 {
		t.Error("output file is empty")
	}
}

func TestWriteOutputFile_BadPath(t *testing.T) {
	streamStats := vcf.ProcessingStatistics{}
	authResults := &vcf.AuthenticityResults{}
	err := writeOutputFile("/nonexistent/dir/results.json", streamStats, authResults)
	if err == nil {
		t.Error("writeOutputFile(bad path) expected error, got nil")
	}
}

// Ensure time and vcf imports are used
var _ = time.Second
