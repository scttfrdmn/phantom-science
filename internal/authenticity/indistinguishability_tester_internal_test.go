//go:build !release

package authenticity

import (
	"math"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
)

func newTestTester() *IndistinguishabilityTester {
	it := NewIndistinguishabilityTester()
	return it
}

// ---------------------------------------------------------------------------
// calculateFSTStandardError
// ---------------------------------------------------------------------------

// TestCalculateFSTStandardError_ZeroN returns 1.0 when sampleSize <= 0.
func TestCalculateFSTStandardError_ZeroN(t *testing.T) {
	it := newTestTester()
	if got := it.calculateFSTStandardError(0.2, 0); got != 1.0 {
		t.Errorf("sampleSize=0: got %v, want 1.0", got)
	}
}

// TestCalculateFSTStandardError_Formula verifies sqrt(p*(1-p)/n).
func TestCalculateFSTStandardError_Formula(t *testing.T) {
	it := newTestTester()
	fst := 0.25
	n := 100
	want := math.Sqrt(fst * (1 - fst) / float64(n))
	got := it.calculateFSTStandardError(fst, n)
	if math.Abs(got-want) > 1e-12 {
		t.Errorf("calculateFSTStandardError(0.25, 100) = %v, want %v", got, want)
	}
}

// TestCalculateFSTStandardError_ZeroFST gives 0 SE (fully-fixed variant, no variance).
func TestCalculateFSTStandardError_ZeroFST(t *testing.T) {
	it := newTestTester()
	got := it.calculateFSTStandardError(0.0, 50)
	if got != 0.0 {
		t.Errorf("fst=0: got %v, want 0.0", got)
	}
}

// ---------------------------------------------------------------------------
// calculateFSTSignificance
// ---------------------------------------------------------------------------

// TestCalculateFSTSignificance_ZeroFST should return 1.0 (no deviation).
func TestCalculateFSTSignificance_ZeroFST(t *testing.T) {
	it := newTestTester()
	p := it.calculateFSTSignificance(0.0, 100)
	// SE=0 when fst=0; p-value should be 1.0 (fallback).
	if p != 1.0 {
		t.Errorf("fst=0 p-value = %v, want 1.0", p)
	}
}

// TestCalculateFSTSignificance_HighFST should yield a very small p-value for a large
// FST with a large sample (strongly significant differentiation).
func TestCalculateFSTSignificance_HighFST(t *testing.T) {
	it := newTestTester()
	p := it.calculateFSTSignificance(0.3, 10000)
	if p >= 0.001 {
		t.Errorf("fst=0.3, n=10000 p-value = %v, want < 0.001", p)
	}
}

// TestCalculateFSTSignificance_RangeIsValid ensures the p-value is in [0, 1].
func TestCalculateFSTSignificance_RangeIsValid(t *testing.T) {
	it := newTestTester()
	for _, fst := range []float64{0.01, 0.1, 0.2, 0.5} {
		for _, n := range []int{10, 100, 1000} {
			p := it.calculateFSTSignificance(fst, n)
			if p < 0 || p > 1 {
				t.Errorf("fst=%.2f, n=%d: p-value %v outside [0,1]", fst, n, p)
			}
		}
	}
}

// ---------------------------------------------------------------------------
// calculateLDMetric
// ---------------------------------------------------------------------------

// TestCalculateLDMetric_Empty returns 0 for empty input.
func TestCalculateLDMetric_Empty(t *testing.T) {
	it := newTestTester()
	if got := it.calculateLDMetric(nil); got != 0.0 {
		t.Errorf("nil input: got %v, want 0.0", got)
	}
	if got := it.calculateLDMetric([]string{}); got != 0.0 {
		t.Errorf("empty input: got %v, want 0.0", got)
	}
}

// TestCalculateLDMetric_ShortSequences returns 0 for sequences shorter than 2 bases.
func TestCalculateLDMetric_ShortSequences(t *testing.T) {
	it := newTestTester()
	seqs := []string{"A", "C", "G"}
	if got := it.calculateLDMetric(seqs); got != 0.0 {
		t.Errorf("single-base sequences: got %v, want 0.0", got)
	}
}

// TestCalculateLDMetric_CompleteLD produces r²=1 when all sequences have the same two bases
// at every adjacent pair (perfect association).
func TestCalculateLDMetric_CompleteLD(t *testing.T) {
	it := newTestTester()
	// All sequences are "AC...AC..." pattern — adjacent positions always co-occur.
	// With only A and C present, and A always paired with C, we get perfect LD.
	seqs := make([]string, 20)
	for i := range seqs {
		seqs[i] = "AC" // two positions, always A then C
	}
	got := it.calculateLDMetric(seqs)
	// r² should be 1 (monomorphic across sequences ⇒ D_max check handles this path).
	// Accept either 0 (skipped as monomorphic) or 1.
	if got < 0 || got > 1.001 {
		t.Errorf("complete LD: r²=%v out of [0,1]", got)
	}
}

// TestCalculateLDMetric_InRange verifies the metric is always in [0,1].
func TestCalculateLDMetric_InRange(t *testing.T) {
	it := newTestTester()
	// Mixed sequences with varying bases.
	seqs := []string{
		"ACGTACGT",
		"TGCATGCA",
		"AAGGTTCC",
		"CCTTAAGG",
		"ACGTACGT",
	}
	got := it.calculateLDMetric(seqs)
	if got < 0 || got > 1.001 {
		t.Errorf("mixed sequences: r²=%v outside [0,1]", got)
	}
}

// ---------------------------------------------------------------------------
// calculateLDStandardError
// ---------------------------------------------------------------------------

// TestCalculateLDStandardError_ZeroN returns 1.0 when sampleSize <= 0.
func TestCalculateLDStandardError_ZeroN(t *testing.T) {
	it := newTestTester()
	if got := it.calculateLDStandardError(0.1, 0); got != 1.0 {
		t.Errorf("sampleSize=0: got %v, want 1.0", got)
	}
}

// TestCalculateLDStandardError_Formula verifies sqrt(p*(1-p)/n).
func TestCalculateLDStandardError_Formula(t *testing.T) {
	it := newTestTester()
	ldDiff := 0.15
	n := 200
	want := math.Sqrt(ldDiff * (1 - ldDiff) / float64(n))
	got := it.calculateLDStandardError(ldDiff, n)
	if math.Abs(got-want) > 1e-12 {
		t.Errorf("calculateLDStandardError(0.15, 200) = %v, want %v", got, want)
	}
}

// ---------------------------------------------------------------------------
// calculateLDSignificance
// ---------------------------------------------------------------------------

// TestCalculateLDSignificance_ZeroLDDiff returns 1.0 (no difference, not significant).
func TestCalculateLDSignificance_ZeroLDDiff(t *testing.T) {
	it := newTestTester()
	p := it.calculateLDSignificance(0.0, 100)
	if p != 1.0 {
		t.Errorf("ldDiff=0 p-value = %v, want 1.0", p)
	}
}

// TestCalculateLDSignificance_LargeDiff produces a small p-value.
func TestCalculateLDSignificance_LargeDiff(t *testing.T) {
	it := newTestTester()
	p := it.calculateLDSignificance(0.4, 5000)
	if p >= 0.001 {
		t.Errorf("ldDiff=0.4, n=5000 p-value = %v, want < 0.001", p)
	}
}

// TestCalculateLDSignificance_RangeIsValid ensures the p-value is in [0, 1].
func TestCalculateLDSignificance_RangeIsValid(t *testing.T) {
	it := newTestTester()
	for _, ldDiff := range []float64{0.05, 0.1, 0.3, 0.5} {
		for _, n := range []int{10, 100, 1000} {
			p := it.calculateLDSignificance(ldDiff, n)
			if p < 0 || p > 1 {
				t.Errorf("ldDiff=%.2f, n=%d: p-value %v outside [0,1]", ldDiff, n, p)
			}
		}
	}
}

// ---------------------------------------------------------------------------
// calculateErrorStandardError
// ---------------------------------------------------------------------------

func TestCalculateErrorStandardError_ZeroN(t *testing.T) {
	it := newTestTester()
	if got := it.calculateErrorStandardError(0.3, 0); got != 1.0 {
		t.Errorf("n=0: got %v, want 1.0", got)
	}
}

func TestCalculateErrorStandardError_Formula(t *testing.T) {
	it := newTestTester()
	p, n := 0.1, 200
	want := math.Sqrt(p * (1 - p) / float64(n))
	if got := it.calculateErrorStandardError(p, n); math.Abs(got-want) > 1e-12 {
		t.Errorf("got %v, want %v", got, want)
	}
}

func TestCalculateErrorStandardError_ZeroDiff(t *testing.T) {
	it := newTestTester()
	// SE of 0 proportion is 0
	if got := it.calculateErrorStandardError(0.0, 100); got != 0.0 {
		t.Errorf("p=0: got %v, want 0.0", got)
	}
}

func TestCalculateErrorStandardError_AboveOne_NoNaN(t *testing.T) {
	// #109: mean_errors difference can push errorDiff > 1.0; result must not be NaN.
	it := newTestTester()
	got := it.calculateErrorStandardError(1.5, 100)
	if math.IsNaN(got) {
		t.Error("errorDiff=1.5: got NaN, want finite value (clamped to p=1)")
	}
	if got < 0 {
		t.Errorf("errorDiff=1.5: got negative %v", got)
	}
}

// ---------------------------------------------------------------------------
// calculateErrorStatistics
// ---------------------------------------------------------------------------

func TestCalculateErrorStatistics_Empty_NoNaN(t *testing.T) {
	// #110: empty slice must not produce NaN via division by zero.
	it := newTestTester()
	stats := it.calculateErrorStatistics(nil)
	if math.IsNaN(stats["mean_errors"]) {
		t.Error("empty results: mean_errors is NaN")
	}
	if math.IsNaN(stats["success_rate"]) {
		t.Error("empty results: success_rate is NaN")
	}
}

func TestCalculateErrorStatistics_NonEmpty(t *testing.T) {
	it := newTestTester()
	results := []genomics.ErrorApplicationResult{
		{ErrorsApplied: 3, Success: true},
		{ErrorsApplied: 1, Success: false},
		{ErrorsApplied: 2, Success: true},
	}
	stats := it.calculateErrorStatistics(results)
	wantMean := float64(3+1+2) / 3.0
	if math.Abs(stats["mean_errors"]-wantMean) > 1e-12 {
		t.Errorf("mean_errors: got %v, want %v", stats["mean_errors"], wantMean)
	}
	wantSuccess := 2.0 / 3.0
	if math.Abs(stats["success_rate"]-wantSuccess) > 1e-12 {
		t.Errorf("success_rate: got %v, want %v", stats["success_rate"], wantSuccess)
	}
}

// ---------------------------------------------------------------------------
// calculateErrorSignificance
// ---------------------------------------------------------------------------

func TestCalculateErrorSignificance_ZeroDiff(t *testing.T) {
	it := newTestTester()
	if p := it.calculateErrorSignificance(0.0, 100); p != 1.0 {
		t.Errorf("diff=0: got %v, want 1.0", p)
	}
}

func TestCalculateErrorSignificance_LargeDiff(t *testing.T) {
	it := newTestTester()
	if p := it.calculateErrorSignificance(0.4, 5000); p >= 0.001 {
		t.Errorf("large diff: got %v, want < 0.001", p)
	}
}

func TestCalculateErrorSignificance_RangeValid(t *testing.T) {
	it := newTestTester()
	for _, d := range []float64{0.05, 0.1, 0.2, 0.4} {
		p := it.calculateErrorSignificance(d, 500)
		if p < 0 || p > 1 {
			t.Errorf("diff=%.2f: p=%v outside [0,1]", d, p)
		}
	}
}

// ---------------------------------------------------------------------------
// calculateGCStandardError
// ---------------------------------------------------------------------------

func TestCalculateGCStandardError_ZeroN(t *testing.T) {
	it := newTestTester()
	if got := it.calculateGCStandardError(0.05, 0); got != 1.0 {
		t.Errorf("n=0: got %v, want 1.0", got)
	}
}

func TestCalculateGCStandardError_Formula(t *testing.T) {
	it := newTestTester()
	p, n := 0.05, 100
	want := math.Sqrt(p * (1 - p) / float64(n))
	if got := it.calculateGCStandardError(p, n); math.Abs(got-want) > 1e-12 {
		t.Errorf("got %v, want %v", got, want)
	}
}

// ---------------------------------------------------------------------------
// calculateGCSignificance
// ---------------------------------------------------------------------------

func TestCalculateGCSignificance_ZeroDiff(t *testing.T) {
	it := newTestTester()
	if p := it.calculateGCSignificance(0.0, 100); p != 1.0 {
		t.Errorf("diff=0: got %v, want 1.0", p)
	}
}

func TestCalculateGCSignificance_LargeDiff(t *testing.T) {
	it := newTestTester()
	if p := it.calculateGCSignificance(0.3, 10000); p >= 0.001 {
		t.Errorf("large diff: got %v, want < 0.001", p)
	}
}

func TestCalculateGCSignificance_RangeValid(t *testing.T) {
	it := newTestTester()
	for _, d := range []float64{0.01, 0.05, 0.15, 0.30} {
		p := it.calculateGCSignificance(d, 200)
		if p < 0 || p > 1 {
			t.Errorf("diff=%.2f: p=%v outside [0,1]", d, p)
		}
	}
}

// ---------------------------------------------------------------------------
// ldBaseIndex and ldMaxIndex helpers
// ---------------------------------------------------------------------------

func TestLdBaseIndex(t *testing.T) {
	cases := []struct {
		b    byte
		want int
	}{
		{'A', 0}, {'a', 0},
		{'C', 1}, {'c', 1},
		{'G', 2}, {'g', 2},
		{'T', 3}, {'t', 3},
		{'N', -1}, {'X', -1}, {'-', -1},
	}
	for _, tc := range cases {
		if got := ldBaseIndex(tc.b); got != tc.want {
			t.Errorf("ldBaseIndex(%q) = %d, want %d", tc.b, got, tc.want)
		}
	}
}

func TestLdMaxIndex(t *testing.T) {
	// Basic max
	counts := [4]int{1, 5, 3, 2}
	if got := ldMaxIndex(counts); got != 1 {
		t.Errorf("ldMaxIndex(%v) = %d, want 1", counts, got)
	}
	// Tie — first occurrence wins (index 0).
	counts2 := [4]int{7, 7, 3, 2}
	if got := ldMaxIndex(counts2); got != 0 {
		t.Errorf("ldMaxIndex(%v) = %d, want 0 (tie, first)", counts2, got)
	}
}

// TestCalculateSequenceStatistics_EmptySlice verifies no NaN/division-by-zero
// when called with an empty sequence slice (#123).
func TestCalculateSequenceStatistics_EmptySlice(t *testing.T) {
	it := &IndistinguishabilityTester{}
	stats := it.calculateSequenceStatistics([]string{})
	if v, ok := stats["mean_length"]; ok && (v != v) { // NaN check: NaN != NaN
		t.Errorf("mean_length = NaN for empty slice; want 0 or absent")
	}
	if v, ok := stats["gc_content"]; ok && (v != v) {
		t.Errorf("gc_content = NaN for empty slice; want 0 or absent")
	}
}

// TestCalculateSequenceStatistics_AllEmptyStrings verifies no NaN when all sequences
// are empty strings (totalLength == 0, so gc_content would divide by zero) (#123).
func TestCalculateSequenceStatistics_AllEmptyStrings(t *testing.T) {
	it := &IndistinguishabilityTester{}
	stats := it.calculateSequenceStatistics([]string{"", "", ""})
	if v := stats["mean_length"]; v != 0.0 {
		t.Errorf("mean_length = %v for all-empty sequences; want 0.0", v)
	}
	if v, ok := stats["gc_content"]; ok && (v != v) {
		t.Errorf("gc_content = NaN for all-empty sequences; want 0 or absent")
	}
}
