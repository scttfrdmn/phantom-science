//go:build !release
// +build !release

package authenticity

import (
	"math"
	"testing"
	"time"

	"math/rand/v2"
)

// ---------------------------------------------------------------------------
// simulateModelAccuracy — nearest-centroid classifier (#90)
// ---------------------------------------------------------------------------

func newTestADS() *AdversarialDefenseSystem {
	src := rand.NewPCG(uint64(time.Now().UnixNano()), 0)
	return &AdversarialDefenseSystem{
		rng: rand.New(src),
	}
}

// buildSeparableFeatures returns len*2 feature vectors: class 0 clusters near (0,0),
// class 1 clusters near (1,1). With 80/20 split the centroid classifier should
// achieve ~100% accuracy.
func buildSeparableFeatures(n int) ([][]float64, []int) {
	features := make([][]float64, 2*n)
	labels := make([]int, 2*n)
	for i := 0; i < n; i++ {
		features[i] = []float64{0.0 + float64(i)*0.001, 0.0 + float64(i)*0.001}
		labels[i] = 0
	}
	for i := 0; i < n; i++ {
		features[n+i] = []float64{1.0 + float64(i)*0.001, 1.0 + float64(i)*0.001}
		labels[n+i] = 1
	}
	return features, labels
}

// buildIdenticalFeatures returns len*2 feature vectors all equal — classifier
// cannot distinguish, accuracy should be ~0.5.
func buildIdenticalFeatures(n int) ([][]float64, []int) {
	features := make([][]float64, 2*n)
	labels := make([]int, 2*n)
	for i := 0; i < n; i++ {
		features[i] = []float64{0.5, 0.5}
		labels[i] = 0
	}
	for i := 0; i < n; i++ {
		features[n+i] = []float64{0.5, 0.5}
		labels[n+i] = 1
	}
	return features, labels
}

// TestSimulateModelAccuracy_SeparableClasses verifies that clearly separable
// classes yield accuracy well above random (>= 0.8).
func TestSimulateModelAccuracy_SeparableClasses(t *testing.T) {
	ads := newTestADS()
	features, labels := buildSeparableFeatures(20) // 40 samples total
	got := ads.simulateModelAccuracy(nil, features, labels, nil)
	if got < 0.8 {
		t.Errorf("separable classes: accuracy = %.4f, want >= 0.80", got)
	}
}

// TestSimulateModelAccuracy_IdenticalFeatures verifies that indistinguishable
// classes yield accuracy near random (the centroid for both classes is the same
// point, so ties are broken by map iteration — result is in [0,1]).
// We just check that it's a valid probability.
func TestSimulateModelAccuracy_IdenticalFeatures(t *testing.T) {
	ads := newTestADS()
	features, labels := buildIdenticalFeatures(20)
	got := ads.simulateModelAccuracy(nil, features, labels, nil)
	if got < 0 || got > 1 {
		t.Errorf("identical features: accuracy = %.4f outside [0,1]", got)
	}
}

// TestSimulateModelAccuracy_TooFewSamples verifies the guard for < 10 samples.
func TestSimulateModelAccuracy_TooFewSamples(t *testing.T) {
	ads := newTestADS()
	features := [][]float64{{1, 0}, {0, 1}}
	labels := []int{0, 1}
	got := ads.simulateModelAccuracy(nil, features, labels, nil)
	if got != 0.5 {
		t.Errorf("too few samples: got %.4f, want 0.5", got)
	}
}

// TestSimulateModelAccuracy_DefenseIgnored verifies that passing a non-nil defense
// does not change the outcome (defense is applied upstream, not here).
func TestSimulateModelAccuracy_DefenseIgnored(t *testing.T) {
	ads := newTestADS()
	features, labels := buildSeparableFeatures(20)
	noDefense := ads.simulateModelAccuracy(nil, features, labels, nil)
	withDefense := ads.simulateModelAccuracy(nil, features, labels, &DefenseStrategy{
		CountermeasurePower: 1.0,
		EffectivenessScore:  1.0,
	})
	if noDefense != withDefense {
		t.Errorf("defense should be ignored: noDefense=%.4f, withDefense=%.4f", noDefense, withDefense)
	}
}

// TestSimulateModelAccuracy_MultiClass verifies three-class classification.
func TestSimulateModelAccuracy_MultiClass(t *testing.T) {
	ads := newTestADS()
	// Three clearly separated clusters of 10 points each = 30 total.
	n := 10
	features := make([][]float64, 3*n)
	labels := make([]int, 3*n)
	for i := 0; i < n; i++ {
		features[i] = []float64{0.0, float64(i) * 0.001}
		labels[i] = 0
		features[n+i] = []float64{10.0, float64(i) * 0.001}
		labels[n+i] = 1
		features[2*n+i] = []float64{20.0, float64(i) * 0.001}
		labels[2*n+i] = 2
	}
	got := ads.simulateModelAccuracy(nil, features, labels, nil)
	if got < 0.8 {
		t.Errorf("3-class separable: accuracy = %.4f, want >= 0.80", got)
	}
}

// ---------------------------------------------------------------------------
// calculateSignificance — two-proportion z-test (#91)
// ---------------------------------------------------------------------------

// TestCalculateSignificance_IdenticalAccuracies verifies p=1.0 when diff=0.
func TestCalculateSignificance_IdenticalAccuracies(t *testing.T) {
	ads := newTestADS()
	got := ads.calculateSignificance(0.7, 0.7, 100)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("identical accuracies: p = %.6f, want 1.0", got)
	}
}

// TestCalculateSignificance_ZeroSampleSize verifies p=1.0 for n<=0.
func TestCalculateSignificance_ZeroSampleSize(t *testing.T) {
	ads := newTestADS()
	for _, n := range []int{0, -1} {
		got := ads.calculateSignificance(0.7, 0.5, n)
		if got != 1.0 {
			t.Errorf("sampleSize=%d: p = %.6f, want 1.0", n, got)
		}
	}
}

// TestCalculateSignificance_LargeDiff verifies a large difference at big n gives p<<0.05.
func TestCalculateSignificance_LargeDiff(t *testing.T) {
	ads := newTestADS()
	// diff=0.3, n=1000: SE = sqrt(2*0.65*0.35/1000)≈0.0213; z≈14 → p≈0
	got := ads.calculateSignificance(0.8, 0.5, 1000)
	if got >= 0.05 {
		t.Errorf("large diff: p = %.6f, want < 0.05", got)
	}
}

// TestCalculateSignificance_SmallDiff verifies a tiny difference at small n gives p>>0.05.
func TestCalculateSignificance_SmallDiff(t *testing.T) {
	ads := newTestADS()
	// diff=0.01, n=100: SE=sqrt(2*0.505*0.495/100)≈0.0707; z=0.141 → p≈0.888
	got := ads.calculateSignificance(0.51, 0.50, 100)
	if got <= 0.5 {
		t.Errorf("small diff: p = %.6f, want > 0.5", got)
	}
}

// TestCalculateSignificance_UsesPooledSE verifies the formula is not /0.05.
// At n=100, diff=0.01: old formula gave z=2.0 (p≈0.046); correct formula gives z≈0.14 (p≈0.888).
func TestCalculateSignificance_UsesPooledSE(t *testing.T) {
	ads := newTestADS()
	p := ads.calculateSignificance(0.51, 0.50, 100)
	// If the hardcoded /0.05 were still used, z=2.0 and p≈0.046 (< 0.5).
	if p <= 0.5 {
		t.Errorf("wrong SE: p = %.4f suggests hardcoded 0.05 divisor still in use (want > 0.5)", p)
	}
}

// ---------------------------------------------------------------------------
// calculateLDProxy — real LD metric (#92)
// ---------------------------------------------------------------------------

// TestCalculateLDProxy_EmptyInput verifies 0.0 for empty/single-sequence inputs.
func TestCalculateLDProxy_EmptyInput(t *testing.T) {
	ads := newTestADS()
	for _, seqs := range [][]string{nil, {}, {"ACGT"}} {
		got := ads.calculateLDProxy(seqs)
		if got != 0.0 {
			t.Errorf("seqs=%v: got %.4f, want 0.0", seqs, got)
		}
	}
}

// TestCalculateLDProxy_CompleteLD verifies r²=1 when all sequences are identical
// — every adjacent pair has the same haplotype, so D=0 and r²=0 (monomorphic → skipped).
// Actually with identical sequences there's no variation, so result = 0.0.
func TestCalculateLDProxy_MonomorphicIsZero(t *testing.T) {
	ads := newTestADS()
	seqs := []string{"AAAA", "AAAA", "AAAA", "AAAA"}
	got := ads.calculateLDProxy(seqs)
	if got != 0.0 {
		t.Errorf("monomorphic: got %.4f, want 0.0", got)
	}
}

// TestCalculateLDProxy_PerfectLD constructs sequences where two adjacent positions
// are in perfect linkage: position 0 is always A or T; position 1 tracks position 0
// perfectly (A→A, T→T). r² should be 1.0.
func TestCalculateLDProxy_PerfectLD(t *testing.T) {
	ads := newTestADS()
	// Half: "AA..." — half "TT..." — perfect co-segregation at positions 0,1.
	n := 10
	seqs := make([]string, 2*n)
	for i := 0; i < n; i++ {
		seqs[i] = "AACGT"
		seqs[n+i] = "TTCGT"
	}
	got := ads.calculateLDProxy(seqs)
	// Position 0,1: nAB=n (A,A), nAb=0, naB=0, nab=n (T,T) → D=1/2-1/2*1/2=1/4; r²=1
	if math.Abs(got-1.0) > 0.01 {
		t.Errorf("perfect LD at pos 0,1: got %.4f, want ≈1.0", got)
	}
}

// TestCalculateLDProxy_UsesAllSequences verifies that the function uses ALL
// sequences (not just sequences[0]). With 5 identical "AA" and 5 identical "TT"
// there IS variation across sequences, so LD proxy > 0.
func TestCalculateLDProxy_UsesAllSequences(t *testing.T) {
	ads := newTestADS()
	seqs := []string{"AACG", "AACG", "AACG", "AACG", "AACG", "TTCG", "TTCG", "TTCG", "TTCG", "TTCG"}
	got := ads.calculateLDProxy(seqs)
	// Old bug: only used sequences[0]="AACG" → no cross-sequence variation → returns 0.
	// New implementation: sees variation at pos 0 and 1 → r²=1.0 for that pair.
	if got <= 0 {
		t.Errorf("multi-sequence variation: got %.4f, want > 0 (old bug returned 0 by only reading sequences[0])", got)
	}
}

// TestCalculateLDProxy_InRange verifies result is always in [0,1].
func TestCalculateLDProxy_InRange(t *testing.T) {
	ads := newTestADS()
	cases := [][]string{
		{"ACGT", "TGCA", "AAAA", "CCCC", "GGGG", "TTTT"},
		{"ACGTACGT", "TGCATGCA", "AAACCCGG", "TTTAGGCC"},
	}
	for _, seqs := range cases {
		got := ads.calculateLDProxy(seqs)
		if got < 0 || got > 1.0+1e-9 {
			t.Errorf("seqs=%v: got %.4f outside [0,1]", seqs, got)
		}
	}
}
