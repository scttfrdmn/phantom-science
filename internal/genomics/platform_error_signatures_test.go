//go:build !release
// +build !release

package genomics

import (
	"math"
	"testing"
)

// TestGenerateLinearDecay_LengthOne is a regression test for #118: when
// length==1, the old code divided by float64(1-1)==0, producing +Inf slope
// and NaN via Inf*0 in the first (and only) element.
func TestGenerateLinearDecay_LengthOne(t *testing.T) {
	m := NewPlatformErrorSignatureManager()
	pattern := m.generateLinearDecay(1, 0.01, 0.02)
	if len(pattern) != 1 {
		t.Fatalf("expected len 1, got %d", len(pattern))
	}
	if math.IsNaN(pattern[0]) || math.IsInf(pattern[0], 0) {
		t.Errorf("generateLinearDecay(1, 0.01, 0.02)[0] = %v; want finite (startRate)", pattern[0])
	}
	// Single-element pattern should return startRate.
	if pattern[0] != 0.01 {
		t.Errorf("expected startRate 0.01, got %v", pattern[0])
	}
}

// TestGenerateLinearDecay_LengthZero verifies an empty slice is returned.
func TestGenerateLinearDecay_LengthZero(t *testing.T) {
	m := NewPlatformErrorSignatureManager()
	pattern := m.generateLinearDecay(0, 0.01, 0.02)
	if len(pattern) != 0 {
		t.Errorf("expected empty slice for length 0, got len %d", len(pattern))
	}
}

// TestGenerateLinearDecay_Normal verifies the standard multi-element case.
func TestGenerateLinearDecay_Normal(t *testing.T) {
	m := NewPlatformErrorSignatureManager()
	pattern := m.generateLinearDecay(5, 0.0, 1.0)
	if len(pattern) != 5 {
		t.Fatalf("expected len 5, got %d", len(pattern))
	}
	// Should linearly interpolate: 0.0, 0.25, 0.5, 0.75, 1.0
	expected := []float64{0.0, 0.25, 0.5, 0.75, 1.0}
	for i, want := range expected {
		if math.Abs(pattern[i]-want) > 1e-9 {
			t.Errorf("pattern[%d] = %v, want %v", i, pattern[i], want)
		}
	}
}
