//go:build !release
// +build !release

package features

import (
	"math"
	"testing"
)

// makeTrinucleotideFreq builds a trinucleotideFreq map where every trinucleotide
// listed in 'keys' has the given frequency, and all others are absent (zero).
func makeTrinucleotideFreq(keys []string, freq float64) map[string]float64 {
	m := make(map[string]float64, len(keys))
	for _, k := range keys {
		m[k] = freq
	}
	return m
}

// allTrinucleotides returns all 64 trinucleotide strings over ACGT.
func allTrinucleotides() []string {
	var out []string
	for _, b1 := range []byte("ACGT") {
		for _, b2 := range []byte("ACGT") {
			for _, b3 := range []byte("ACGT") {
				out = append(out, string([]byte{b1, b2, b3}))
			}
		}
	}
	return out
}

// ── Tests for calculateTrinucleotideBias ─────────────────────────────────────

func TestCalculateTrinucleotideBias_EmptyMap(t *testing.T) {
	msm := &MutationalSignatureModel{}
	score := msm.calculateTrinucleotideBias(map[string]float64{})
	// All 64 observed values are 0.0; expected = 1/64. chi-square is non-zero
	// (maximum bias: no trinucleotide observed at all) → score approaches 1.
	if score <= 0.0 || score > 1.0 {
		t.Errorf("empty map: bias score out of [0,1]: %f", score)
	}
}

func TestCalculateTrinucleotideBias_UniformDistribution(t *testing.T) {
	// When all 64 trinucleotides have equal frequency 1/64, chi-square == 0
	// and bias score should be exactly 0.
	freqs := make(map[string]float64, 64)
	for _, tri := range allTrinucleotides() {
		freqs[tri] = 1.0 / 64.0
	}
	msm := &MutationalSignatureModel{}
	score := msm.calculateTrinucleotideBias(freqs)
	if score > 1e-9 {
		t.Errorf("uniform distribution: expected bias score ~0, got %f", score)
	}
}

func TestCalculateTrinucleotideBias_SingleTrinucleotide(t *testing.T) {
	// All probability mass on one trinucleotide (e.g., "AAA" = 1.0, rest = 0).
	// This is the maximally biased distribution. With df=63 normalisation and a
	// /10 dampening factor the score is ~0.095, which must exceed the uniform
	// baseline (score = 0).
	freqs := map[string]float64{"AAA": 1.0}
	msm := &MutationalSignatureModel{}
	score := msm.calculateTrinucleotideBias(freqs)
	if score <= 0.0 {
		t.Errorf("single trinucleotide: expected bias score > 0, got %f", score)
	}
	if score > 1.0 {
		t.Errorf("single trinucleotide: expected bias score ≤ 1.0, got %f", score)
	}
}

func TestCalculateTrinucleotideBias_AbsentTrinucleotidesRaiseScore(t *testing.T) {
	// Distributing the same total probability over fewer trinucleotides should
	// give a higher bias score than a uniform distribution, because the absent
	// ones contribute (0 - 1/64)^2 / (1/64) each to chi-square.
	//
	// Half the trinucleotides present, each at 2/64 (so they sum to 1.0):
	tris := allTrinucleotides()
	half := make(map[string]float64, 32)
	for i := 0; i < 32; i++ {
		half[tris[i]] = 2.0 / 64.0
	}
	msm := &MutationalSignatureModel{}
	halfScore := msm.calculateTrinucleotideBias(half)

	// Uniform (all 64 present at 1/64) should have score == 0.
	fullUniform := make(map[string]float64, 64)
	for _, tri := range tris {
		fullUniform[tri] = 1.0 / 64.0
	}
	uniformScore := msm.calculateTrinucleotideBias(fullUniform)

	if halfScore <= uniformScore {
		t.Errorf("half-support map (score %f) should be more biased than uniform (score %f)",
			halfScore, uniformScore)
	}
}

func TestCalculateTrinucleotideBias_ScoreInRange(t *testing.T) {
	// Bias score must always be in [0, 1].
	cases := []map[string]float64{
		{},
		{"ACG": 0.5, "TGA": 0.5},
		{"AAA": 0.3, "TTT": 0.3, "GGG": 0.2, "CCC": 0.2},
	}
	msm := &MutationalSignatureModel{}
	for i, c := range cases {
		score := msm.calculateTrinucleotideBias(c)
		if score < 0.0 || score > 1.0 {
			t.Errorf("case %d: score %f out of [0,1]", i, score)
		}
	}
}

// ---------------------------------------------------------------------------
// assessSignatureConsistency — Bessel's correction (#93)
// ---------------------------------------------------------------------------

// TestAssessSignatureConsistency_UsesSampleVariance verifies ÷(N-1) is used.
//
// Contributions: {A: 0.8, B: 0.8, C: 0.2}  mean = 0.6
//   sample variance  = ((0.2²+0.2²+0.4²) / 2) = 0.24/2 = 0.12 → cv = sqrt(0.12)/0.6 ≈ 0.5774
//   pop    variance  = 0.24/3 = 0.08                            → cv = sqrt(0.08)/0.6 ≈ 0.4714
//   sample consistency = exp(-0.5774) ≈ 0.5614
//   pop    consistency = exp(-0.4714) ≈ 0.6243
func TestAssessSignatureConsistency_UsesSampleVariance(t *testing.T) {
	msm := &MutationalSignatureModel{}
	contribs := map[string]float64{"A": 0.8, "B": 0.8, "C": 0.2}
	got := msm.assessSignatureConsistency(contribs)

	wantSample := math.Exp(-math.Sqrt(0.12) / 0.6)
	wantPop := math.Exp(-math.Sqrt(0.08) / 0.6)

	if math.Abs(got-wantSample) > 1e-9 {
		t.Errorf("consistency = %v, want %v (sample variance); population variance would give %v",
			got, wantSample, wantPop)
	}
}

// TestAssessSignatureConsistency_IdenticalValues verifies CV=0 → consistency=1.
func TestAssessSignatureConsistency_IdenticalValues(t *testing.T) {
	msm := &MutationalSignatureModel{}
	contribs := map[string]float64{"A": 0.5, "B": 0.5, "C": 0.5}
	got := msm.assessSignatureConsistency(contribs)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("identical values: got %v, want 1.0", got)
	}
}

// TestAssessSignatureConsistency_SingleSignature verifies count<2 → 1.0.
func TestAssessSignatureConsistency_SingleSignature(t *testing.T) {
	msm := &MutationalSignatureModel{}
	got := msm.assessSignatureConsistency(map[string]float64{"SBS1": 0.9})
	if got != 1.0 {
		t.Errorf("single signature: got %v, want 1.0", got)
	}
}

// TestAssessSignatureConsistency_InRange verifies result is always in (0,1].
func TestAssessSignatureConsistency_InRange(t *testing.T) {
	msm := &MutationalSignatureModel{}
	cases := []map[string]float64{
		{"A": 0.9, "B": 0.1},
		{"A": 0.5, "B": 0.5},
		{"A": 0.3, "B": 0.6, "C": 0.1},
		{"A": 1.0, "B": 0.0},
	}
	for _, c := range cases {
		got := msm.assessSignatureConsistency(c)
		if got <= 0 || got > 1 {
			t.Errorf("contributions=%v: consistency %v outside (0,1]", c, got)
		}
	}
}
