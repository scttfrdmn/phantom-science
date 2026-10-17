//go:build !release
// +build !release

package features

import (
	"math"
	"testing"
)

// ---- calculateCorrelation tests ----

// TestCalculateCorrelation_PerfectPositive verifies r=1.0 for perfectly correlated data.
func TestCalculateCorrelation_PerfectPositive(t *testing.T) {
	x := []float64{1.0, 2.0, 3.0, 4.0, 5.0}
	y := []float64{2.0, 4.0, 6.0, 8.0, 10.0} // y = 2*x

	got := calculateCorrelation(x, y)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("perfect positive correlation: got %v, want 1.0", got)
	}
}

// TestCalculateCorrelation_PerfectNegative verifies r=-1.0 for perfectly inverse data.
func TestCalculateCorrelation_PerfectNegative(t *testing.T) {
	x := []float64{1.0, 2.0, 3.0, 4.0, 5.0}
	y := []float64{10.0, 8.0, 6.0, 4.0, 2.0} // y = 12 - 2*x

	got := calculateCorrelation(x, y)
	if math.Abs(got-(-1.0)) > 1e-9 {
		t.Errorf("perfect negative correlation: got %v, want -1.0", got)
	}
}

// TestCalculateCorrelation_ZeroCorrelation verifies that data constructed to be
// exactly orthogonal returns exactly 0.0.
//
// With x=[1,2,3] and y=[1,-2,1]:
//   meanX=2, meanY=0
//   sum = (1-2)*(1-0) + (2-2)*(-2-0) + (3-2)*(1-0) = -1 + 0 + 1 = 0
//   → correlation = 0
func TestCalculateCorrelation_ZeroCorrelation(t *testing.T) {
	x := []float64{1.0, 2.0, 3.0}
	y := []float64{1.0, -2.0, 1.0}

	got := calculateCorrelation(x, y)
	if got != 0.0 {
		t.Errorf("zero correlation: got %v, want 0.0", got)
	}
}

// TestCalculateCorrelation_ExactZeroConstantX verifies that a constant x series
// (zero variance in x) returns 0.0.
func TestCalculateCorrelation_ExactZeroConstantX(t *testing.T) {
	x := []float64{5.0, 5.0, 5.0, 5.0}
	y := []float64{1.0, 2.0, 3.0, 4.0}

	got := calculateCorrelation(x, y)
	if got != 0.0 {
		t.Errorf("constant x series: got %v, want 0.0", got)
	}
}

// TestCalculateCorrelation_ExactZeroConstantY verifies that a constant y series
// (zero variance in y) returns 0.0.
func TestCalculateCorrelation_ExactZeroConstantY(t *testing.T) {
	x := []float64{1.0, 2.0, 3.0, 4.0}
	y := []float64{7.0, 7.0, 7.0, 7.0}

	got := calculateCorrelation(x, y)
	if got != 0.0 {
		t.Errorf("constant y series: got %v, want 0.0", got)
	}
}

// TestCalculateCorrelation_BothConstant verifies that two constant series return 0.0.
func TestCalculateCorrelation_BothConstant(t *testing.T) {
	x := []float64{3.0, 3.0, 3.0}
	y := []float64{7.0, 7.0, 7.0}

	got := calculateCorrelation(x, y)
	if got != 0.0 {
		t.Errorf("both constant: got %v, want 0.0", got)
	}
}

// TestCalculateCorrelation_LengthMismatch verifies that mismatched slice lengths return 0.0.
func TestCalculateCorrelation_LengthMismatch(t *testing.T) {
	x := []float64{1.0, 2.0, 3.0}
	y := []float64{1.0, 2.0}

	got := calculateCorrelation(x, y)
	if got != 0.0 {
		t.Errorf("length mismatch: got %v, want 0.0", got)
	}
}

// TestCalculateCorrelation_TooShort verifies that a single-element slice returns 0.0.
func TestCalculateCorrelation_TooShort(t *testing.T) {
	x := []float64{1.0}
	y := []float64{2.0}

	got := calculateCorrelation(x, y)
	if got != 0.0 {
		t.Errorf("single element: got %v, want 0.0", got)
	}
}

// TestCalculateCorrelation_EmptySlices verifies that empty slices return 0.0.
func TestCalculateCorrelation_EmptySlices(t *testing.T) {
	got := calculateCorrelation([]float64{}, []float64{})
	if got != 0.0 {
		t.Errorf("empty slices: got %v, want 0.0", got)
	}
}

// TestCalculateCorrelation_TwoPoints verifies exact correlation with two data points.
func TestCalculateCorrelation_TwoPoints(t *testing.T) {
	// Two points always yield perfect correlation (±1) unless constant.
	x := []float64{0.0, 1.0}
	y := []float64{0.0, 2.0}

	got := calculateCorrelation(x, y)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("two point positive: got %v, want 1.0", got)
	}
}

// TestCalculateCorrelation_Range verifies the output is always in [-1, 1].
func TestCalculateCorrelation_Range(t *testing.T) {
	cases := []struct {
		x, y []float64
	}{
		{[]float64{1, 2, 3, 4, 5}, []float64{2, 4, 6, 8, 10}},
		{[]float64{1, 2, 3, 4, 5}, []float64{10, 8, 6, 4, 2}},
		{[]float64{1, 1, 2, 3, 5, 8}, []float64{0, 1, 1, 2, 3, 5}},
		{[]float64{2, 4, 4, 4, 5, 5, 7, 9}, []float64{4, 4, 4, 6, 5, 5, 5, 7}},
	}

	for _, c := range cases {
		got := calculateCorrelation(c.x, c.y)
		if got < -1.0-1e-9 || got > 1.0+1e-9 {
			t.Errorf("correlation out of range [-1,1]: got %v for x=%v y=%v", got, c.x, c.y)
		}
	}
}

// ---- identifyEvolutionaryOrigin tests ----

// newTestPhylogeneticModel creates a minimal PhylogeneticDistanceModel for testing
// identifyEvolutionaryOrigin without requiring the full initialization.
func newTestPhylogeneticModel() *PhylogeneticDistanceModel {
	return &PhylogeneticDistanceModel{
		speciesDistances:  make(map[string]map[string]float64),
		substitutionRates: make(map[string]float64),
		branchLengths:     make(map[string]float64),
	}
}

// TestIdentifyEvolutionaryOrigin_SingleMinimum verifies that the species with the
// strictly smallest distance is returned.
func TestIdentifyEvolutionaryOrigin_SingleMinimum(t *testing.T) {
	pdm := newTestPhylogeneticModel()

	distances := map[string]float64{
		"Human":      0.05,
		"Chimpanzee": 0.12,
		"Gorilla":    0.25,
		"Macaque":    0.80,
	}

	got := pdm.identifyEvolutionaryOrigin(distances)
	if got != "Human" {
		t.Errorf("identifyEvolutionaryOrigin: got %q, want %q", got, "Human")
	}
}

// TestIdentifyEvolutionaryOrigin_AnotherMinimum verifies a different minimum distance.
func TestIdentifyEvolutionaryOrigin_AnotherMinimum(t *testing.T) {
	pdm := newTestPhylogeneticModel()

	distances := map[string]float64{
		"Human":      0.90,
		"Chimpanzee": 0.88,
		"Marmoset":   0.02, // smallest
		"Macaque":    0.50,
	}

	got := pdm.identifyEvolutionaryOrigin(distances)
	if got != "Marmoset" {
		t.Errorf("identifyEvolutionaryOrigin: got %q, want %q", got, "Marmoset")
	}
}

// TestIdentifyEvolutionaryOrigin_ZeroDistance verifies that a zero-distance species
// is correctly identified as the closest.
func TestIdentifyEvolutionaryOrigin_ZeroDistance(t *testing.T) {
	pdm := newTestPhylogeneticModel()

	distances := map[string]float64{
		"Human":      0.0,
		"Chimpanzee": 0.1,
	}

	got := pdm.identifyEvolutionaryOrigin(distances)
	if got != "Human" {
		t.Errorf("zero-distance species: got %q, want %q", got, "Human")
	}
}

// TestIdentifyEvolutionaryOrigin_EmptyMap verifies that an empty map returns an
// empty string (no closest species found).
func TestIdentifyEvolutionaryOrigin_EmptyMap(t *testing.T) {
	pdm := newTestPhylogeneticModel()

	got := pdm.identifyEvolutionaryOrigin(map[string]float64{})
	if got != "" {
		t.Errorf("empty map: got %q, want %q", got, "")
	}
}

// TestIdentifyEvolutionaryOrigin_SingleEntry verifies a map with one species.
func TestIdentifyEvolutionaryOrigin_SingleEntry(t *testing.T) {
	pdm := newTestPhylogeneticModel()

	distances := map[string]float64{
		"Orangutan": 0.42,
	}

	got := pdm.identifyEvolutionaryOrigin(distances)
	if got != "Orangutan" {
		t.Errorf("single entry: got %q, want %q", got, "Orangutan")
	}
}

// TestIdentifyEvolutionaryOrigin_NegativeDistance verifies behaviour when distances
// include a negative value (degenerate but should still find the minimum).
func TestIdentifyEvolutionaryOrigin_NegativeDistance(t *testing.T) {
	pdm := newTestPhylogeneticModel()

	distances := map[string]float64{
		"A": 0.5,
		"B": -0.1, // smallest (most negative)
		"C": 0.3,
	}

	got := pdm.identifyEvolutionaryOrigin(distances)
	if got != "B" {
		t.Errorf("negative distance: got %q, want %q", got, "B")
	}
}

// ---- calculateOverallScore tests ----

// newTestBiologicalModelingEngine returns a zero-value BiologicalModelingEngine
// suitable for calling calculateOverallScore and categorizeRealism directly.
func newTestBiologicalModelingEngine() *BiologicalModelingEngine {
	return &BiologicalModelingEngine{}
}

// TestCalculateOverallScore_EmptyMap verifies that an empty feature map returns 0.0.
func TestCalculateOverallScore_EmptyMap(t *testing.T) {
	bme := newTestBiologicalModelingEngine()
	got := bme.calculateOverallScore(map[string]float64{})
	if got != 0.0 {
		t.Errorf("empty map: got %v, want 0.0", got)
	}
}

// TestCalculateOverallScore_UnknownFeatureIsIgnored verifies that features not in
// the weight table do not contribute to the score.
func TestCalculateOverallScore_UnknownFeatureIsIgnored(t *testing.T) {
	bme := newTestBiologicalModelingEngine()
	got := bme.calculateOverallScore(map[string]float64{
		"unknown_feature": 1.0,
		"another_unknown": 0.5,
	})
	// None of these feature names appear in the weight table, so total weight = 0
	if got != 0.0 {
		t.Errorf("unknown features only: got %v, want 0.0", got)
	}
}

// TestCalculateOverallScore_SingleKnownFeature verifies the correct weighted average
// when only one known feature is present.
//
// With a single feature "gc_bias" (weight=0.15) and score=1.0:
//   result = (1.0 * 0.15) / 0.15 = 1.0
func TestCalculateOverallScore_SingleKnownFeature(t *testing.T) {
	bme := newTestBiologicalModelingEngine()
	got := bme.calculateOverallScore(map[string]float64{
		"gc_bias": 1.0,
	})
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("single feature gc_bias=1.0: got %v, want 1.0", got)
	}
}

// TestCalculateOverallScore_SingleFeatureHalfScore verifies that a half-score
// single feature correctly produces a 0.5 result.
func TestCalculateOverallScore_SingleFeatureHalfScore(t *testing.T) {
	bme := newTestBiologicalModelingEngine()
	got := bme.calculateOverallScore(map[string]float64{
		"cpg_depletion": 0.5,
	})
	if math.Abs(got-0.5) > 1e-9 {
		t.Errorf("single feature cpg_depletion=0.5: got %v, want 0.5", got)
	}
}

// TestCalculateOverallScore_AllFeaturesEqualOne verifies that all-1.0 scores produce 1.0.
func TestCalculateOverallScore_AllFeaturesEqualOne(t *testing.T) {
	bme := newTestBiologicalModelingEngine()
	allOnes := map[string]float64{
		"gc_bias":                 1.0,
		"cpg_depletion":           1.0,
		"repetitive_elements":     1.0,
		"codon_usage":             1.0,
		"phylogenetic_distance":   1.0,
		"chromatin_accessibility": 1.0,
		"mutational_signature":    1.0,
		"selection_pressure":      1.0,
	}
	got := bme.calculateOverallScore(allOnes)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("all features=1.0: got %v, want 1.0", got)
	}
}

// TestCalculateOverallScore_AllFeaturesEqualZero verifies that all-0.0 scores produce 0.0.
func TestCalculateOverallScore_AllFeaturesEqualZero(t *testing.T) {
	bme := newTestBiologicalModelingEngine()
	allZeros := map[string]float64{
		"gc_bias":                 0.0,
		"cpg_depletion":           0.0,
		"repetitive_elements":     0.0,
		"codon_usage":             0.0,
		"phylogenetic_distance":   0.0,
		"chromatin_accessibility": 0.0,
		"mutational_signature":    0.0,
		"selection_pressure":      0.0,
	}
	got := bme.calculateOverallScore(allZeros)
	if math.Abs(got) > 1e-9 {
		t.Errorf("all features=0.0: got %v, want 0.0", got)
	}
}

// TestCalculateOverallScore_WeightedAverage verifies the weighted-average arithmetic
// using two known-weight features.
//
// "gc_bias" weight=0.15, "phylogenetic_distance" weight=0.20
// With scores 1.0 and 0.0:
//   weighted_sum = 1.0*0.15 + 0.0*0.20 = 0.15
//   total_weight = 0.35
//   result = 0.15 / 0.35 ≈ 0.42857...
func TestCalculateOverallScore_WeightedAverage(t *testing.T) {
	bme := newTestBiologicalModelingEngine()
	got := bme.calculateOverallScore(map[string]float64{
		"gc_bias":               1.0, // weight 0.15
		"phylogenetic_distance": 0.0, // weight 0.20
	})
	want := 0.15 / 0.35
	if math.Abs(got-want) > 1e-9 {
		t.Errorf("weighted average: got %v, want %v", got, want)
	}
}

// TestCalculateOverallScore_MixedKnownAndUnknown verifies that unknown features
// are silently ignored in the weighted average computation.
func TestCalculateOverallScore_MixedKnownAndUnknown(t *testing.T) {
	bme := newTestBiologicalModelingEngine()
	// "gc_bias" (weight 0.15, score 0.8) + "totally_unknown" (ignored)
	// result = (0.8 * 0.15) / 0.15 = 0.8
	got := bme.calculateOverallScore(map[string]float64{
		"gc_bias":         0.8,
		"totally_unknown": 0.9,
	})
	want := 0.8
	if math.Abs(got-want) > 1e-9 {
		t.Errorf("mixed known/unknown: got %v, want %v", got, want)
	}
}

// ---- categorizeRealism tests ----

// TestCategorizeRealism_HighlyRealistic verifies scores at and above the 0.85 threshold.
func TestCategorizeRealism_HighlyRealistic(t *testing.T) {
	bme := newTestBiologicalModelingEngine()

	cases := []float64{0.85, 0.90, 0.95, 1.0}
	for _, score := range cases {
		got := bme.categorizeRealism(score)
		if got != HighlyRealistic {
			t.Errorf("score=%.2f: got %v (%q), want HighlyRealistic", score, got, got.String())
		}
	}
}

// TestCategorizeRealism_ModeratelyRealistic verifies the [0.70, 0.85) range.
func TestCategorizeRealism_ModeratelyRealistic(t *testing.T) {
	bme := newTestBiologicalModelingEngine()

	cases := []float64{0.70, 0.75, 0.80, 0.84, 0.8499}
	for _, score := range cases {
		got := bme.categorizeRealism(score)
		if got != ModeratelyRealistic {
			t.Errorf("score=%.4f: got %v (%q), want ModeratelyRealistic", score, got, got.String())
		}
	}
}

// TestCategorizeRealism_SomewhatRealistic verifies the [0.50, 0.70) range.
func TestCategorizeRealism_SomewhatRealistic(t *testing.T) {
	bme := newTestBiologicalModelingEngine()

	cases := []float64{0.50, 0.55, 0.60, 0.65, 0.6999}
	for _, score := range cases {
		got := bme.categorizeRealism(score)
		if got != SomewhatRealistic {
			t.Errorf("score=%.4f: got %v (%q), want SomewhatRealistic", score, got, got.String())
		}
	}
}

// TestCategorizeRealism_PoorlyRealistic verifies scores below 0.50.
func TestCategorizeRealism_PoorlyRealistic(t *testing.T) {
	bme := newTestBiologicalModelingEngine()

	cases := []float64{0.0, 0.10, 0.25, 0.49, 0.4999}
	for _, score := range cases {
		got := bme.categorizeRealism(score)
		if got != PoorlyRealistic {
			t.Errorf("score=%.4f: got %v (%q), want PoorlyRealistic", score, got, got.String())
		}
	}
}

// TestCategorizeRealism_ExactBoundaries verifies the exact threshold values.
// 0.85 → HighlyRealistic (>=0.85)
// 0.70 → ModeratelyRealistic (>=0.70 but <0.85)
// 0.50 → SomewhatRealistic (>=0.50 but <0.70)
// 0.49 → PoorlyRealistic (<0.50)
func TestCategorizeRealism_ExactBoundaries(t *testing.T) {
	bme := newTestBiologicalModelingEngine()

	tests := []struct {
		score float64
		want  RealismCategory
	}{
		{0.85, HighlyRealistic},
		{0.84999, ModeratelyRealistic},
		{0.70, ModeratelyRealistic},
		{0.69999, SomewhatRealistic},
		{0.50, SomewhatRealistic},
		{0.49999, PoorlyRealistic},
	}

	for _, tc := range tests {
		got := bme.categorizeRealism(tc.score)
		if got != tc.want {
			t.Errorf("categorizeRealism(%.5f) = %v (%q), want %v (%q)",
				tc.score, got, got.String(), tc.want, tc.want.String())
		}
	}
}

// TestCategorizeRealism_AllCategoriesReturned verifies that every RealismCategory
// constant can be produced by categorizeRealism.
func TestCategorizeRealism_AllCategoriesReturned(t *testing.T) {
	bme := newTestBiologicalModelingEngine()

	got := map[RealismCategory]bool{}
	for _, score := range []float64{0.0, 0.50, 0.70, 0.85} {
		got[bme.categorizeRealism(score)] = true
	}

	for _, cat := range []RealismCategory{PoorlyRealistic, SomewhatRealistic, ModeratelyRealistic, HighlyRealistic} {
		if !got[cat] {
			t.Errorf("category %v (%q) was never returned by categorizeRealism", cat, cat.String())
		}
	}
}

// ---- calculateAccessibilityScore tests (via ChromatinAccessibilityModel) ----

// TestCalculateAccessibilityScore_AllStates verifies that every defined ChromatinState
// produces the expected accessibility value from the lookup table.
func TestCalculateAccessibilityScore_AllStates(t *testing.T) {
	// Create a minimal model (no need to call the full constructor; the method
	// only reads the hard-coded accessibilityMap inside it).
	cam := &ChromatinAccessibilityModel{}

	tests := []struct {
		state ChromatinState
		want  float64
	}{
		{ActivePromoter, 0.90},
		{WeakPromoter, 0.70},
		{PoisedPromoter, 0.60},
		{StrongEnhancer, 0.85},
		{WeakEnhancer, 0.65},
		{Txn_Transition, 0.50},
		{Txn_Elongation, 0.55},
		{Weak_Txn, 0.40},
		{Polycomb_Repressed, 0.15},
		{Heterochromatin_low, 0.10},
		{Repetitive_CNV, 0.08},
		{Inactive_Low, 0.20},
	}

	for _, tc := range tests {
		t.Run(tc.state.String(), func(t *testing.T) {
			got := cam.calculateAccessibilityScore(tc.state)
			if got != tc.want {
				t.Errorf("state %v (%q): got %v, want %v", tc.state, tc.state.String(), got, tc.want)
			}
		})
	}
}

// TestCalculateAccessibilityScore_UnknownStateReturnsDefault verifies that an
// undefined ChromatinState value returns the default 0.20.
func TestCalculateAccessibilityScore_UnknownStateReturnsDefault(t *testing.T) {
	cam := &ChromatinAccessibilityModel{}

	// Use a value beyond the last defined constant (Inactive_Low = 11).
	unknownState := ChromatinState(99)
	got := cam.calculateAccessibilityScore(unknownState)
	want := 0.20
	if got != want {
		t.Errorf("unknown state %d: got %v, want %v", unknownState, got, want)
	}
}

// TestCalculateAccessibilityScore_ActivePromoterHighestOpen verifies that
// ActivePromoter has strictly higher accessibility than closed chromatin states.
func TestCalculateAccessibilityScore_ActivePromoterHighestOpen(t *testing.T) {
	cam := &ChromatinAccessibilityModel{}

	activeScore := cam.calculateAccessibilityScore(ActivePromoter)
	closedStates := []ChromatinState{Polycomb_Repressed, Heterochromatin_low, Repetitive_CNV, Inactive_Low}

	for _, state := range closedStates {
		if cam.calculateAccessibilityScore(state) >= activeScore {
			t.Errorf("expected ActivePromoter (%.2f) > %v (%.2f)",
				activeScore, state, cam.calculateAccessibilityScore(state))
		}
	}
}

// TestCalculateAccessibilityScore_Range verifies every score is within [0, 1].
func TestCalculateAccessibilityScore_Range(t *testing.T) {
	cam := &ChromatinAccessibilityModel{}

	allStates := []ChromatinState{
		ActivePromoter, WeakPromoter, PoisedPromoter, StrongEnhancer, WeakEnhancer,
		Txn_Transition, Txn_Elongation, Weak_Txn, Polycomb_Repressed,
		Heterochromatin_low, Repetitive_CNV, Inactive_Low,
	}

	for _, state := range allStates {
		got := cam.calculateAccessibilityScore(state)
		if got < 0.0 || got > 1.0 {
			t.Errorf("state %v: score %v out of [0, 1]", state, got)
		}
	}
}

// ---- calculateDnDsConsistency tests ----

func newSPM() *SelectionPressureModel {
	m, _ := NewSelectionPressureModel()
	return m
}

func TestDnDsConsistency_NilInfo(t *testing.T) {
	spm := newSPM()
	details := &SelectionAnalysisDetails{DnDsInformation: nil}
	if got := spm.calculateDnDsConsistency(details); got != 0.0 {
		t.Errorf("nil DnDsInformation: expected 0.0, got %v", got)
	}
}

func TestDnDsConsistency_PerfectMatch(t *testing.T) {
	spm := newSPM()
	// StrongPurifying expected dN/dS = 0.1 — perfect match → score = 1.0
	details := &SelectionAnalysisDetails{
		DominantSelection: StrongPurifying,
		DnDsInformation:   &DnDsInfo{DnDsRatio: 0.1},
	}
	got := spm.calculateDnDsConsistency(details)
	if math.Abs(got-1.0) > 1e-9 {
		t.Errorf("perfect match: expected 1.0, got %v", got)
	}
}

func TestDnDsConsistency_CloseMatch(t *testing.T) {
	spm := newSPM()
	// NeutralEvolution expected 1.0; observed 1.1 → small penalty
	details := &SelectionAnalysisDetails{
		DominantSelection: NeutralEvolution,
		DnDsInformation:   &DnDsInfo{DnDsRatio: 1.1},
	}
	got := spm.calculateDnDsConsistency(details)
	if got <= 0.9 || got >= 1.0 {
		t.Errorf("close match (1.1 vs 1.0): expected in (0.9, 1.0), got %v", got)
	}
}

func TestDnDsConsistency_LargeDeviation(t *testing.T) {
	spm := newSPM()
	// StrongPurifying expected 0.1; observed 2.0 → large log-ratio penalty
	details := &SelectionAnalysisDetails{
		DominantSelection: StrongPurifying,
		DnDsInformation:   &DnDsInfo{DnDsRatio: 2.0},
	}
	got := spm.calculateDnDsConsistency(details)
	if got >= 0.3 {
		t.Errorf("large deviation (2.0 vs 0.1): expected < 0.3, got %v", got)
	}
}

func TestDnDsConsistency_ZeroRatio(t *testing.T) {
	spm := newSPM()
	details := &SelectionAnalysisDetails{
		DominantSelection: NeutralEvolution,
		DnDsInformation:   &DnDsInfo{DnDsRatio: 0.0},
	}
	if got := spm.calculateDnDsConsistency(details); got != 0.0 {
		t.Errorf("zero ratio: expected 0.0, got %v", got)
	}
}

func TestDnDsConsistency_Continuous(t *testing.T) {
	spm := newSPM()
	// Score must decrease monotonically as the observed ratio diverges from expected
	base := 1.0 // NeutralEvolution expected
	prev := 1.1 // start just above 1.0
	for _, obs := range []float64{2.0, 5.0, 10.0, 50.0} {
		d := &SelectionAnalysisDetails{
			DominantSelection: NeutralEvolution,
			DnDsInformation:   &DnDsInfo{DnDsRatio: obs},
		}
		got := spm.calculateDnDsConsistency(d)
		_ = base
		if got >= prev {
			t.Errorf("expected score to decrease as ratio grows: obs=%v score=%v prev=%v", obs, got, prev)
		}
		prev = got
	}
}
