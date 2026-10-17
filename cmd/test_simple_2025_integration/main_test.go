//go:build !release
// +build !release

package main

import (
	"testing"
)

func TestGenerateTestSequence_Length(t *testing.T) {
	seq := generateTestSequence(50, 0.1)
	if len(seq) != 50 {
		t.Errorf("generateTestSequence(50, 0.1) len = %d, want 50", len(seq))
	}
}

func TestGenerateTestSequence_BaseValue(t *testing.T) {
	// With error rate 0.5, every other element might be modified
	seq := generateTestSequence(10, 0.5)
	for i, v := range seq {
		if v < 0.0 || v > 1.0 {
			t.Errorf("seq[%d] = %f, out of expected range [0,1]", i, v)
		}
	}
}

func TestCalculateComprehensiveAssessment_Positive(t *testing.T) {
	score := calculateComprehensiveAssessment()
	if score <= 0.0 {
		t.Errorf("calculateComprehensiveAssessment() = %f, want > 0", score)
	}
}

func TestCalculateIndistinguishabilityScore_NearHalf(t *testing.T) {
	score := calculateIndistinguishabilityScore()
	// Should be close to 0.5 (target is ~50% classification accuracy)
	if score < 0.4 || score > 0.6 {
		t.Errorf("calculateIndistinguishabilityScore() = %f, expected near 0.5", score)
	}
}

func TestAbs_Positive(t *testing.T) {
	if got := abs(5.5); got != 5.5 {
		t.Errorf("abs(5.5) = %f, want 5.5", got)
	}
}

func TestAbs_Negative(t *testing.T) {
	if got := abs(-5.5); got != 5.5 {
		t.Errorf("abs(-5.5) = %f, want 5.5", got)
	}
}
