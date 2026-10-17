//go:build !release
// +build !release

package features

import (
	"math"
	"testing"
)

// TestAnalyzeSubstitutionPatterns_PureTransitions is a regression test for the
// bug fixed in #117: when a sequence contains only transitions (A↔G, C↔T) and
// no transversions, the function previously computed float64(transitions)/0,
// producing +Inf which propagated through math.Exp to NaN.
func TestAnalyzeSubstitutionPatterns_PureTransitions(t *testing.T) {
	model, err := NewPhylogeneticDistanceModel(0.5)
	if err != nil {
		t.Fatalf("NewPhylogeneticDistanceModel: %v", err)
	}

	// Sequence of alternating AG pairs → pure transitions, zero transversions.
	seq := []byte("AGAGAGAGAGAGAGAGAGAG") // 20 bases, len>=10 passes the early guard

	score := model.analyzeSubstitutionPatterns(seq)

	if math.IsNaN(score) || math.IsInf(score, 0) {
		t.Errorf("analyzeSubstitutionPatterns returned %v for pure-transition sequence; want finite [0,1]", score)
	}
	if score < 0 || score > 1 {
		t.Errorf("score %v out of expected range [0,1]", score)
	}
	// When transversions==0 we return the neutral 0.5
	if score != 0.5 {
		t.Errorf("expected 0.5 for undefined ratio, got %v", score)
	}
}

// TestAnalyzeSubstitutionPatterns_ShortSeq verifies the early-return guard.
func TestAnalyzeSubstitutionPatterns_ShortSeq(t *testing.T) {
	model, err := NewPhylogeneticDistanceModel(0.5)
	if err != nil {
		t.Fatalf("NewPhylogeneticDistanceModel: %v", err)
	}

	score := model.analyzeSubstitutionPatterns([]byte("ACGT"))
	if score != 0.0 {
		t.Errorf("expected 0.0 for short sequence, got %v", score)
	}
}

// TestAnalyzeSubstitutionPatterns_BothTypes verifies a mixed sequence returns
// a finite score in [0,1].
func TestAnalyzeSubstitutionPatterns_BothTypes(t *testing.T) {
	model, err := NewPhylogeneticDistanceModel(0.5)
	if err != nil {
		t.Fatalf("NewPhylogeneticDistanceModel: %v", err)
	}

	// Mix: transitions (AG) and transversions (AC, AT, GC, GT)
	seq := []byte("ACAGTATGCGAGATCTATGC")

	score := model.analyzeSubstitutionPatterns(seq)

	if math.IsNaN(score) || math.IsInf(score, 0) {
		t.Errorf("analyzeSubstitutionPatterns returned %v; want finite", score)
	}
	if score < 0 || score > 1 {
		t.Errorf("score %v out of [0,1]", score)
	}
}
