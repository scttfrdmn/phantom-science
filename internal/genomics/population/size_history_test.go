//go:build !release

package population

import (
	"math"
	"math/rand/v2"
	"testing"
)

// adjusted returns N(at), rounded, for a present size ne and events.
func adjusted(t *testing.T, ne int, events []DemographicEvent, at float64) int {
	t.Helper()
	n, _ := sizeAt(float64(ne), sortedEvents(events), at)
	return int(math.Round(n))
}

// An earlier size function kept a bottleneck in force forever; it must end
// after Duration generations.
func TestGetAdjustedPopSize_BottleneckRespectsDuration(t *testing.T) {
	ev := []DemographicEvent{{
		Type: "bottleneck", Time: 100, Duration: 50,
		Parameters: map[string]float64{"population_size": 10},
	}}
	for _, tc := range []struct {
		at   float64
		want int
	}{{0, 1000}, {99, 1000}, {100, 10}, {149, 10}, {150, 1000}, {10000, 1000}} {
		if got := adjusted(t, 1000, ev, tc.at); got != tc.want {
			t.Errorf("N(%v) = %d, want %d", tc.at, got, tc.want)
		}
	}
}

// Duration 0 is a persistent size change from Time into the past.
func TestGetAdjustedPopSize_ZeroDurationBottleneckIsPersistent(t *testing.T) {
	ev := []DemographicEvent{{Type: "bottleneck", Time: 100, Parameters: map[string]float64{"population_size": 10}}}
	if got := adjusted(t, 1000, ev, 1e6); got != 10 {
		t.Errorf("N(1e6) = %d, want 10", got)
	}
}

// An expansion that ended Time generations ago means the population was
// SMALLER before it (the old code multiplied N further back in time).
func TestGetAdjustedPopSize_ExpansionShrinksBackwards(t *testing.T) {
	ev := []DemographicEvent{{Type: "expansion", Time: 100, Parameters: map[string]float64{"expansion_factor": 10}}}
	if got := adjusted(t, 1000, ev, 50); got != 1000 {
		t.Errorf("N(50) = %d, want 1000", got)
	}
	if got := adjusted(t, 1000, ev, 150); got != 100 {
		t.Errorf("N(150) = %d, want 100", got)
	}
}

// Expansion with a Duration is exponential growth over the window.
func TestGetAdjustedPopSize_ExponentialGrowth(t *testing.T) {
	ev := []DemographicEvent{{Type: "expansion", Time: 0, Duration: 1000, Parameters: map[string]float64{"expansion_factor": math.E}}}
	for _, tc := range []struct {
		at   float64
		want float64
	}{{0, 10000}, {500, 10000 * math.Exp(-0.5)}, {1000, 10000 / math.E}, {5000, 10000 / math.E}} {
		if got := adjusted(t, 10000, ev, tc.at); math.Abs(float64(got)-tc.want) > 1 {
			t.Errorf("N(%v) = %d, want %.0f", tc.at, got, tc.want)
		}
	}
	// growth_rate is used when no expansion_factor is given.
	ev = []DemographicEvent{{Type: "expansion", Time: 0, Duration: 100, Parameters: map[string]float64{"growth_rate": 0.01}}}
	if got, want := adjusted(t, 10000, ev, 100), 10000*math.Exp(-1); math.Abs(float64(got)-want) > 1 {
		t.Errorf("growth_rate: N(100) = %d, want %.0f", got, want)
	}
}

func TestGetAdjustedPopSize_IgnoresNonSizeEvents(t *testing.T) {
	ev := []DemographicEvent{
		{Type: "migration", Time: 10, Duration: 100, Parameters: map[string]float64{"rate": 0.1}},
		{Type: "admixture", Time: 10, Parameters: map[string]float64{"population_size": 3}},
	}
	if got := adjusted(t, 1000, ev, 50); got != 1000 {
		t.Errorf("N(50) = %d, want 1000", got)
	}
}

// nextEventTime must invert the cumulative hazard exactly: for a constant
// epoch followed by a different constant epoch, P(T > t) = exp(-Λ(t)).
func TestNextEventTime_ExactAcrossEpochs(t *testing.T) {
	epochs := []epoch{{start: 0, n0: 100}, {start: 50, n0: 1000}}
	// Λ(50) = 50/100 = 0.5; an exponential draw of 0.25 lands at 25, 0.75
	// lands at 50 + 0.25·1000 = 300.
	if got := nextEventTime(epochs, 1, 0, math.Inf(1), 0.25); math.Abs(got-25) > 1e-9 {
		t.Errorf("got %v, want 25", got)
	}
	if got := nextEventTime(epochs, 1, 0, math.Inf(1), 0.75); math.Abs(got-300) > 1e-9 {
		t.Errorf("got %v, want 300", got)
	}
	if got := nextEventTime(epochs, 1, 0, 200, 0.75); !math.IsInf(got, 1) {
		t.Errorf("event beyond tmax: got %v, want +Inf", got)
	}
	// Exponential epoch N(t) = 1000·exp(-0.01 t): Λ(u) = (exp(0.01u)-1)/(1000·0.01).
	g := []epoch{{start: 0, n0: 1000, g: 0.01}}
	u := 100.0
	h := (math.Exp(0.01*u) - 1) / (1000 * 0.01)
	if got := nextEventTime(g, 1, 0, math.Inf(1), h); math.Abs(got-u) > 1e-6 {
		t.Errorf("exponential epoch: got %v, want %v", got, u)
	}
}

// meanT2 returns the mean and standard error of the pairwise coalescence time
// drawn with nextEventTime over the epochs of the given size history. Two
// lineages in a diploid population of size N(t) coalesce at rate 1/(2N(t)).
func meanT2(t *testing.T, ne int, events []DemographicEvent, reps int) (float64, float64) {
	t.Helper()
	epochs := buildEpochs(float64(ne), events)
	rng := rand.New(rand.NewPCG(11, 0))
	sum, sumSq := 0.0, 0.0
	for i := 0; i < reps; i++ {
		t2 := nextEventTime(epochs, 0.5, 0, math.Inf(1), rng.ExpFloat64())
		if math.IsInf(t2, 1) {
			t.Fatal("no coalescence")
		}
		sum += t2
		sumSq += t2 * t2
	}
	m := sum / float64(reps)
	sd := math.Sqrt((sumSq - float64(reps)*m*m) / float64(reps-1))
	return m, sd / math.Sqrt(float64(reps))
}

// The pairwise coalescence time must match the analytic
// expectation E[T2] = ∫ exp(-∫ 1/(2N(s)) ds) dt under each size history.
func TestPairwiseCoalescenceTime_MatchesTheory(t *testing.T) {
	const reps = 6000
	cases := []struct {
		name   string
		events []DemographicEvent
		want   float64
	}{
		{"constant", nil, 2000},
		{
			// N=100 for t >= 1000: 2000(1-e^-0.5) + e^-0.5·200.
			"instantaneous expansion",
			[]DemographicEvent{{Type: "expansion", Time: 1000, Parameters: map[string]float64{"expansion_factor": 10}}},
			2000*(1-math.Exp(-0.5)) + math.Exp(-0.5)*200,
		},
		{
			// N=100 during [200,300) only.
			"finite bottleneck",
			[]DemographicEvent{{Type: "bottleneck", Time: 200, Duration: 100, Parameters: map[string]float64{"population_size": 100}}},
			2000*(1-math.Exp(-0.1)) + math.Exp(-0.1)*(200*(1-math.Exp(-0.5))+math.Exp(-0.5)*2000),
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m, se := meanT2(t, 1000, tc.events, reps)
			if math.Abs(m-tc.want) > 4*se {
				t.Errorf("E[T2] = %.1f ± %.1f, want %.1f (4σ)", m, se, tc.want)
			}
		})
	}
}

func TestPoissonExact_MeanAndVariance(t *testing.T) {
	rng := rand.New(rand.NewPCG(1, 2))
	for _, lambda := range []float64{0.5, 7, 45} {
		const n = 20000
		sum, sumSq := 0.0, 0.0
		for i := 0; i < n; i++ {
			k := float64(poissonExact(rng, lambda))
			sum += k
			sumSq += k * k
		}
		m := sum / n
		v := sumSq/n - m*m
		if se := math.Sqrt(lambda / n); math.Abs(m-lambda) > 4*se {
			t.Errorf("λ=%v: mean %.3f, want %.3f ± %.3f", lambda, m, lambda, 4*se)
		}
		if math.Abs(v/lambda-1) > 0.1 {
			t.Errorf("λ=%v: variance %.3f, want ≈ %.3f", lambda, v, lambda)
		}
	}
}
