//go:build !release
// +build !release

package features

import (
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// calculateSelectionConsistency — Bessel's correction (#87)
// ---------------------------------------------------------------------------

// buildSPMWithValues returns a SelectionPressureModel whose purifyingSelection
// map is pre-populated so that calculateSelectionConsistency(chromosome, start, end)
// will draw exactly the given set of values.
//
// For start=0, end=15000 and chromosome="chr1":
//   numSamples = 15000/5000 = 3, step = 5000
//   i=0: pos=0    → key "chr1:0-10000"
//   i=1: pos=5000  → key "chr1:0-10000"
//   i=2: pos=10000 → key "chr1:10000-20000"
//
// So setting purifyingSelection["chr1:0-10000"]=v1 and
// purifyingSelection["chr1:10000-20000"]=v2 gives values [v1, v1, v2].
func buildSPMWithValues(v1, v2 float64) *SelectionPressureModel {
	return &SelectionPressureModel{
		purifyingSelection: map[string]float64{
			"chr1:0-10000":     v1,
			"chr1:10000-20000": v2,
		},
		positiveSelection:  make(map[string]float64),
		balancingSelection: make(map[string]float64),
		dNdSRatios:         make(map[string]float64),
		conservationScores: make(map[string]float64),
	}
}

// TestCalculateSelectionConsistency_UsesSampleVariance verifies that the
// coefficient of variation is computed with Bessel's correction (÷N-1, not ÷N).
//
// Values: [0.8, 0.8, 0.2]  mean = 0.6
//   sample variance = 0.24/2 = 0.12  → std ≈ 0.34641 → cv ≈ 0.57735
//   pop    variance = 0.24/3 = 0.08  → std ≈ 0.28284 → cv ≈ 0.47140
//   sample consistency = exp(-0.57735) ≈ 0.5614
//   pop    consistency = exp(-0.47140) ≈ 0.6243
func TestCalculateSelectionConsistency_UsesSampleVariance(t *testing.T) {
	spm := buildSPMWithValues(0.8, 0.2)
	got := spm.calculateSelectionConsistency("chr1", 0, 15000)

	// Expected with Bessel's correction
	wantSample := math.Exp(-math.Sqrt(0.12) / 0.6)
	// Wrong value if population variance had been used
	wantPop := math.Exp(-math.Sqrt(0.08) / 0.6)

	if math.Abs(got-wantSample) > 1e-9 {
		t.Errorf("consistency = %v, want %v (sample variance); population-variance would give %v",
			got, wantSample, wantPop)
	}
}

// TestCalculateSelectionConsistency_IdenticalValues checks that perfectly
// uniform selection pressure gives consistency = 1.0 (CV = 0).
func TestCalculateSelectionConsistency_IdenticalValues(t *testing.T) {
	spm := buildSPMWithValues(0.5, 0.5)
	got := spm.calculateSelectionConsistency("chr1", 0, 15000)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("identical values: got %v, want 1.0", got)
	}
}

// TestCalculateSelectionConsistency_SmallRegion checks that regions that
// produce < 2 samples return 1.0 (the existing guard).
func TestCalculateSelectionConsistency_SmallRegion(t *testing.T) {
	spm := buildSPMWithValues(0.5, 0.3)
	got := spm.calculateSelectionConsistency("chr1", 0, 1000) // numSamples = 0 < 2
	if got != 1.0 {
		t.Errorf("small region: got %v, want 1.0", got)
	}
}

// TestCalculateSelectionConsistency_InRange verifies the result is always in (0, 1].
func TestCalculateSelectionConsistency_InRange(t *testing.T) {
	cases := [][2]float64{
		{0.9, 0.1},
		{0.5, 0.5},
		{1.0, 0.0},
		{0.3, 0.7},
	}
	for _, c := range cases {
		spm := buildSPMWithValues(c[0], c[1])
		got := spm.calculateSelectionConsistency("chr1", 0, 15000)
		if got <= 0 || got > 1 {
			t.Errorf("v1=%v v2=%v: consistency %v outside (0,1]", c[0], c[1], got)
		}
	}
}
