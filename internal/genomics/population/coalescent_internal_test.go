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
