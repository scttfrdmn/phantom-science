//go:build !release
// +build !release

package genomics

import (
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// calculateVariance — Bessel's correction (#94)
// ---------------------------------------------------------------------------

// TestCalculateVariance_UsesSampleVariance verifies ÷(N-1), not ÷N.
//
// Values [2, 4, 6]:  mean = 4
//   sample variance  = ((4+0+4) / 2) = 4.0
//   pop    variance  = ((4+0+4) / 3) ≈ 2.667
func TestCalculateVariance_UsesSampleVariance(t *testing.T) {
	data := []float64{2, 4, 6}
	got := calculateVariance(data)
	wantSample := 4.0
	wantPop := 8.0 / 3.0

	if math.Abs(got-wantSample) > 1e-9 {
		t.Errorf("variance = %v, want %v (sample ÷N-1); population ÷N would give %v",
			got, wantSample, wantPop)
	}
}

// TestCalculateVariance_ZeroForSingleElement verifies len<2 → 0.
func TestCalculateVariance_ZeroForSingleElement(t *testing.T) {
	for _, data := range [][]float64{nil, {}, {42.0}} {
		got := calculateVariance(data)
		if got != 0 {
			t.Errorf("data=%v: got %v, want 0", data, got)
		}
	}
}

// TestCalculateVariance_IdenticalValues verifies all-equal data → variance 0.
func TestCalculateVariance_IdenticalValues(t *testing.T) {
	data := []float64{3.14, 3.14, 3.14, 3.14}
	got := calculateVariance(data)
	if math.Abs(got) > 1e-12 {
		t.Errorf("identical values: got %v, want 0", got)
	}
}

// TestCalculateVariance_KnownValue verifies a hand-computed two-element case.
//
// Values [1, 3]: mean = 2; sample variance = ((1+1)/1) = 2.
func TestCalculateVariance_TwoElements(t *testing.T) {
	got := calculateVariance([]float64{1, 3})
	if math.Abs(got-2.0) > 1e-9 {
		t.Errorf("got %v, want 2.0", got)
	}
}

// TestCalculateVariance_NonNegative verifies variance is always ≥ 0.
func TestCalculateVariance_NonNegative(t *testing.T) {
	cases := [][]float64{
		{0.1, 0.9},
		{1.0, 2.0, 3.0, 4.0, 5.0},
		{-1.0, 0.0, 1.0},
	}
	for _, data := range cases {
		got := calculateVariance(data)
		if got < 0 {
			t.Errorf("data=%v: got negative variance %v", data, got)
		}
	}
}
