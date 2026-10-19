//go:build !release

package population

import (
	"math"
	"testing"
)

// TestCalculateTajimaD_NegativeVariance_NoNaN verifies that a non-positive variance
// returns 0.0 rather than NaN (#112).
func TestCalculateTajimaD_NegativeVariance_NoNaN(t *testing.T) {
	s := &CoalescentSimulator{}
	// Passing a large negative thetaWatterson is not realistic, but it exercises
	// the variance <= 0 guard: e2*thetaW^2 can dominate and go negative for unusual
	// parameter combinations. Verify the result is always finite.
	got := s.calculateTajimaD(0.1, -0.5, 10)
	if math.IsNaN(got) {
		t.Errorf("negative variance path: got NaN, want 0.0")
	}
	if math.IsInf(got, 0) {
		t.Errorf("negative variance path: got Inf, want 0.0")
	}
}

// TestCalculateTajimaD_ZeroVariance_IsZero ensures variance == 0 still returns 0.
func TestCalculateTajimaD_ZeroVariance_IsZero(t *testing.T) {
	s := &CoalescentSimulator{}
	// thetaWatterson = 0 → variance = 0; function must return 0 not NaN.
	got := s.calculateTajimaD(0.0, 0.0, 10)
	if got != 0.0 {
		t.Errorf("zero thetaWatterson: got %v, want 0.0", got)
	}
}

// TestGetPopulationAtTime_NoHistory falls back to "ancestral" with empty history.
func TestGetPopulationAtTime_NoHistory(t *testing.T) {
	s := &CoalescentSimulator{}
	if got := s.getPopulationAtTime(100); got != "ancestral" {
		t.Errorf("empty history: got %q, want \"ancestral\"", got)
	}
}

// TestGetPopulationAtTime_NoPopulations skips events with empty Populations slice.
func TestGetPopulationAtTime_NoPopulations(t *testing.T) {
	s := &CoalescentSimulator{
		demographicHistory: []DemographicEvent{
			{Type: "bottleneck", Time: 50, Duration: 100, Populations: nil},
		},
	}
	if got := s.getPopulationAtTime(80); got != "ancestral" {
		t.Errorf("no populations in event: got %q, want \"ancestral\"", got)
	}
}

// TestGetPopulationAtTime_PointEvent matches t == event.Time.
func TestGetPopulationAtTime_PointEvent(t *testing.T) {
	s := &CoalescentSimulator{
		demographicHistory: []DemographicEvent{
			{Type: "expansion", Time: 100, Duration: 0, Populations: []string{"EUR"}},
		},
	}
	if got := s.getPopulationAtTime(100); got != "EUR" {
		t.Errorf("point event at 100: got %q, want \"EUR\"", got)
	}
	// Time before event → ancestral
	if got := s.getPopulationAtTime(50); got != "ancestral" {
		t.Errorf("before point event: got %q, want \"ancestral\"", got)
	}
}

// TestGetPopulationAtTime_RangedEvent matches times within [start, start+duration].
func TestGetPopulationAtTime_RangedEvent(t *testing.T) {
	s := &CoalescentSimulator{
		demographicHistory: []DemographicEvent{
			{Type: "migration", Time: 50, Duration: 100, Populations: []string{"AFR"}},
		},
	}
	for _, tc := range []struct {
		t    float64
		want string
	}{
		{50, "AFR"},
		{100, "AFR"},
		{150, "AFR"},
		{200, "ancestral"}, // beyond end
		{10, "ancestral"},  // before start
	} {
		if got := s.getPopulationAtTime(tc.t); got != tc.want {
			t.Errorf("t=%.0f: got %q, want %q", tc.t, got, tc.want)
		}
	}
}

// TestGetPopulationAtTime_MostRecentEventWins picks smallest Time among overlapping events.
func TestGetPopulationAtTime_MostRecentEventWins(t *testing.T) {
	s := &CoalescentSimulator{
		demographicHistory: []DemographicEvent{
			{Type: "migration", Time: 200, Duration: 500, Populations: []string{"ANCESTRAL_EUR"}},
			{Type: "migration", Time: 50, Duration: 100, Populations: []string{"CEU"}},
		},
	}
	// t=120 is within both events ([200,700] and [50,150]).
	// Wait: 200..700 does include 120? No: start=200, end=700; t=120 < start=200 so it doesn't.
	// Let me use an overlap that works: t=250 is in [200,700] but not [50,150].
	if got := s.getPopulationAtTime(250); got != "ANCESTRAL_EUR" {
		t.Errorf("only old event: got %q, want \"ANCESTRAL_EUR\"", got)
	}
	// t=75 is in [50,150] only.
	if got := s.getPopulationAtTime(75); got != "CEU" {
		t.Errorf("only recent event: got %q, want \"CEU\"", got)
	}
	// t=250 is in [200,700] only — test above covers this.
	// For a true overlap test, use [100,400] and [50,200] → both cover t=150.
	s2 := &CoalescentSimulator{
		demographicHistory: []DemographicEvent{
			{Type: "expansion", Time: 100, Duration: 300, Populations: []string{"OLD"}},
			{Type: "migration", Time: 50, Duration: 200, Populations: []string{"RECENT"}},
		},
	}
	// t=150 is in [50,250] (RECENT) and [100,400] (OLD). Most recent wins (Time=50).
	if got := s2.getPopulationAtTime(150); got != "RECENT" {
		t.Errorf("overlapping events: got %q, want \"RECENT\"", got)
	}
}

// TestCalculateThetaPi_MeanPairwiseDifferences checks θ_π against its
// definition (Tajima 1983): the mean number of pairwise differences,
// Σ_i i(n−i)ξ_i / C(n,2). For n = 4 with one singleton and one doubleton the
// singleton differs in 3 of the 6 pairs and the doubleton in 4, so
// θ_π = 7/6. Fay & Wu's (2000) θ_H = Σ_i i²ξ_i / C(n,2) = 5/6, so
// H = θ_π − θ_H = 2/6.
func TestCalculateThetaPi_MeanPairwiseDifferences(t *testing.T) {
	s := NewCoalescentSimulator(1000, 1e-8, 0)
	n := 4
	sfs := &AlleleFrequencySpectrum{
		Frequencies: []float64{0.25, 0.5, 0.75, 1.0},
		Counts:      []int{1, 1, 0, 0},
	}
	pi := s.calculateThetaPi(sfs, n)
	if want := 7.0 / 6.0; math.Abs(pi-want) > 1e-12 {
		t.Fatalf("θ_π = %v, want mean pairwise differences %v", pi, want)
	}
	sfs.ThetaPi = pi
	if h, want := s.calculateFayWuH(sfs, n), 2.0/6.0; math.Abs(h-want) > 1e-12 {
		t.Fatalf("Fay-Wu H = %v, want %v", h, want)
	}
}

// TestComputeSiteFrequencySpectrum_BinsExactCounts guards the SFS binning: a
// mutation carried by k of n leaves must land in bin k−1. For n = 22,
// (15.0/22)*22 truncates to 14 in floating point.
func TestComputeSiteFrequencySpectrum_BinsExactCounts(t *testing.T) {
	s := NewCoalescentSimulator(1000, 1e-8, 0)
	n := 22
	g := &Genealogy{MutationEvents: []MutationEvent{{Mutation: Mutation{Frequency: 15.0 / 22.0}}}}
	sfs := s.ComputeSiteFrequencySpectrum(g, n)
	if sfs.Counts[14] != 1 {
		t.Fatalf("15/22 mutation binned as %v, want bin 14", sfs.Counts)
	}
}

// TestCoalescent_NeutralEstimatorsUnbiased checks the standard neutral
// expectations for a constant-size population: E[θ_π] = E[θ_W] = E[θ_H] =
// θ = 4Nμ·L (Tajima 1983; Watterson 1975; Fay & Wu 2000), so E[H] = 0 and
// Tajima's D averages ≈ 0. Tolerances are 4 standard errors of the replicate
// mean.
func TestCoalescent_NeutralEstimatorsUnbiased(t *testing.T) {
	const (
		N    = 1000
		mu   = 2.5e-6
		L    = 1000
		n    = 20
		reps = 2000 // ~50 ms
	)
	theta := 4.0 * N * mu * L // 10

	s := NewCoalescentSimulator(N, mu, 0)
	s.SetSeed(20001)

	var pi, tw, h, d []float64
	for r := 0; r < reps; r++ {
		g, err := s.SimulateGenealogy("rep", n, L)
		if err != nil {
			t.Fatal(err)
		}
		sfs := s.ComputeSiteFrequencySpectrum(g, n)
		pi = append(pi, sfs.ThetaPi)
		tw = append(tw, sfs.ThetaWatterson)
		h = append(h, sfs.FayWuH)
		d = append(d, sfs.TajimaD)
	}

	check := func(name string, xs []float64, want, slack float64) {
		mean, se := replicateMeanSE(xs)
		tol := 4*se + slack
		t.Logf("%s: mean %.3f, SE %.3f (want %.3f ± %.3f)", name, mean, se, want, tol)
		if math.Abs(mean-want) > tol {
			t.Errorf("%s mean = %.3f, want %.3f ± %.3f", name, mean, want, tol)
		}
	}
	check("θ_π", pi, theta, 0)
	check("θ_W", tw, theta, 0)
	check("Fay-Wu H", h, 0, 0)
	// D is a ratio, so E[D] is only approximately 0: at finite θ its mean
	// is slightly negative (about −0.05 to −0.1 here). The 0.1 slack covers
	// that; halving θ_π, the bug this guards, gives a mean D near −2.
	check("Tajima D", d, 0, 0.1)
}

func replicateMeanSE(xs []float64) (float64, float64) {
	m := 0.0
	for _, x := range xs {
		m += x
	}
	m /= float64(len(xs))
	v := 0.0
	for _, x := range xs {
		v += (x - m) * (x - m)
	}
	v /= float64(len(xs) - 1)
	return m, math.Sqrt(v / float64(len(xs)))
}
