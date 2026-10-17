//go:build !release

package statistics

import (
	"math"
	"testing"
)

func TestComputeFST_PerfectFixation(t *testing.T) {
	// One population fixed for the derived allele, the other fixed for the ancestral.
	// Expected F_ST = 1.
	pops := []PopulationSample{
		{AlleleFreq: 1.0, SampleSize: 50},
		{AlleleFreq: 0.0, SampleSize: 50},
	}
	res := ComputeFST(pops)
	if !res.Valid {
		t.Fatal("expected Valid=true for perfect fixation")
	}
	if math.Abs(res.FST-1.0) > 1e-6 {
		t.Errorf("expected FST=1.0 for perfect fixation, got %.6f", res.FST)
	}
}

func TestComputeFST_IdenticalFrequencies(t *testing.T) {
	// Both populations have the same allele frequency → F_ST = 0.
	pops := []PopulationSample{
		{AlleleFreq: 0.3, SampleSize: 100},
		{AlleleFreq: 0.3, SampleSize: 100},
	}
	res := ComputeFST(pops)
	// When MSG = 0, Valid is false (no between-pop variance) and FST = 0.
	if res.FST != 0 {
		t.Errorf("expected FST=0 for identical frequencies, got %.6f", res.FST)
	}
}

func TestComputeFST_LiteratureExample(t *testing.T) {
	// Approximate YRI vs CEU differentiation at a well-known SNP.
	// Typical published F_ST for high-differentiation variants is in the range [0.10, 0.35].
	// Use frequencies roughly consistent with APOE4 (rs429358):
	//   YRI ~0.22, CEU ~0.14 (Auton et al. 2015 1000G Phase 3).
	pops := []PopulationSample{
		{AlleleFreq: 0.22, SampleSize: 108}, // YRI (1000G n≈108)
		{AlleleFreq: 0.14, SampleSize: 99},  // CEU (1000G n≈99)
	}
	res := ComputeFST(pops)
	if !res.Valid {
		t.Fatal("expected Valid=true for literature example")
	}
	// F_ST at this locus should be moderate (< 0.10 is reasonable for a ~8% freq diff).
	if res.FST < 0 || res.FST > 0.50 {
		t.Errorf("literature example FST out of expected range [0, 0.50]: got %.4f", res.FST)
	}
}

func TestComputeFST_ThreePopulations(t *testing.T) {
	// Three populations with different frequencies.
	pops := []PopulationSample{
		{AlleleFreq: 0.1, SampleSize: 50},
		{AlleleFreq: 0.5, SampleSize: 50},
		{AlleleFreq: 0.9, SampleSize: 50},
	}
	res := ComputeFST(pops)
	if !res.Valid {
		t.Fatal("expected Valid=true for three populations")
	}
	if res.FST <= 0 || res.FST > 1 {
		t.Errorf("FST out of [0,1]: %.4f", res.FST)
	}
}

func TestComputeFST_SinglePopulation(t *testing.T) {
	// Fewer than 2 populations → invalid.
	pops := []PopulationSample{
		{AlleleFreq: 0.4, SampleSize: 100},
	}
	res := ComputeFST(pops)
	if res.Valid {
		t.Error("expected Valid=false for single population")
	}
}

func TestComputeFST_ZeroSampleSize(t *testing.T) {
	pops := []PopulationSample{
		{AlleleFreq: 0.5, SampleSize: 0},
		{AlleleFreq: 0.5, SampleSize: 0},
	}
	res := ComputeFST(pops)
	if res.Valid {
		t.Error("expected Valid=false when all sample sizes are zero")
	}
}

func TestComputeMultiLocusFST_Average(t *testing.T) {
	loci := [][]PopulationSample{
		// Locus 1: fixation → FST = 1
		{
			{AlleleFreq: 1.0, SampleSize: 50},
			{AlleleFreq: 0.0, SampleSize: 50},
		},
		// Locus 2: identical → FST = 0 (invalid, excluded)
		{
			{AlleleFreq: 0.3, SampleSize: 50},
			{AlleleFreq: 0.3, SampleSize: 50},
		},
		// Locus 3: partial differentiation
		{
			{AlleleFreq: 0.2, SampleSize: 50},
			{AlleleFreq: 0.8, SampleSize: 50},
		},
	}
	res := ComputeMultiLocusFST(loci)
	if !res.Valid {
		t.Fatal("expected Valid=true for multi-locus with at least two valid loci")
	}
	// Average of fixation locus (1.0) and partial locus (some value > 0).
	// Result should be between 0 and 1.
	if res.FST <= 0 || res.FST > 1 {
		t.Errorf("multi-locus FST out of [0,1]: %.4f", res.FST)
	}
}

func TestComputeMultiLocusFST_AllInvalid(t *testing.T) {
	loci := [][]PopulationSample{
		// Both loci have identical frequencies → all invalid.
		{
			{AlleleFreq: 0.3, SampleSize: 50},
			{AlleleFreq: 0.3, SampleSize: 50},
		},
	}
	res := ComputeMultiLocusFST(loci)
	if res.Valid {
		t.Error("expected Valid=false when all loci are invalid")
	}
}

func TestComputeMultiLocusFST_ManualCalculation(t *testing.T) {
	// Two loci both giving fixation → average should be exactly 1.0.
	loci := [][]PopulationSample{
		{
			{AlleleFreq: 1.0, SampleSize: 50},
			{AlleleFreq: 0.0, SampleSize: 50},
		},
		{
			{AlleleFreq: 1.0, SampleSize: 50},
			{AlleleFreq: 0.0, SampleSize: 50},
		},
	}
	res := ComputeMultiLocusFST(loci)
	if !res.Valid {
		t.Fatal("expected Valid=true")
	}
	if math.Abs(res.FST-1.0) > 1e-6 {
		t.Errorf("expected FST=1.0, got %.6f", res.FST)
	}
}
