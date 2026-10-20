package population

import (
	"math"
	"sort"
)

// Population-size histories.
//
// Time runs backwards: t is measured in generations before the present. A size
// history is a piecewise function N(t) made of epochs; inside an epoch that
// starts at time a the size is
//
//	N(t) = n0 · exp(−g·(t − a))
//
// so g > 0 means the population was smaller further back in time (i.e. it has
// been growing forward in time) and g = 0 is a constant-size epoch.
//
// DemographicEvent semantics (all times are generations ago):
//
//   - "bottleneck": N = Parameters["population_size"] for t in
//     [Time, Time+Duration). When Duration <= 0 the event is a persistent size
//     change: N = population_size for every t >= Time.
//   - "size_change": N = Parameters["population_size"] for every t >= Time.
//   - "expansion": a forward-time expansion that ENDED at Time. Going back in
//     time the population shrinks. With Duration > 0 the growth is exponential
//     over [Time, Time+Duration): the rate is ln(expansion_factor)/Duration when
//     expansion_factor is given, otherwise Parameters["growth_rate"]. With
//     Duration <= 0 the expansion is instantaneous: N(t) = N/expansion_factor
//     for t >= Time.
//   - Any other type ("migration", "admixture", ...) does not change the size
//     and is ignored by the size history (migration and admixture are not
//     modelled by these simulators).
//
// Bottleneck windows take precedence over other events they overlap; when two
// growth windows overlap, the more recent one wins. Sizes are clamped to >= 1.

// DemographicEvent is one event in a population's size history; see the
// semantics above.
type DemographicEvent struct {
	Type        string             // "bottleneck", "expansion", "size_change", "migration", "admixture"
	Time        int                // Generations ago
	Parameters  map[string]float64 // Event-specific parameters
	Populations []string           // Affected populations
	Duration    int                // Duration in generations
}

// epoch is one piece of a piecewise size history.
type epoch struct {
	start float64 // generations ago at which this epoch begins (backwards in time)
	n0    float64 // size at start
	g     float64 // backward decay rate: N(t) = n0·exp(−g·(t−start))
}

// sortedEvents returns a copy of events sorted by Time (most recent first).
func sortedEvents(events []DemographicEvent) []DemographicEvent {
	out := append([]DemographicEvent(nil), events...)
	sort.SliceStable(out, func(i, j int) bool { return out[i].Time < out[j].Time })
	return out
}

// sizeAt evaluates the size history defined by a present-day size and events
// (sorted by Time) at time t. It returns the size and the local backward decay
// rate g (see epoch).
func sizeAt(present float64, events []DemographicEvent, t float64) (float64, float64) {
	n, g := present, 0.0
	baselineDone := false
	override := -1.0

	for _, e := range events {
		start := float64(e.Time)
		if t < start {
			break
		}
		dur := float64(e.Duration)

		switch e.Type {
		case "bottleneck":
			size, ok := e.Parameters["population_size"]
			if !ok || size <= 0 {
				continue
			}
			if dur > 0 {
				if t < start+dur && override < 0 {
					override = size
				}
				continue
			}
			if !baselineDone {
				n = size
			}
		case "size_change":
			size, ok := e.Parameters["population_size"]
			if ok && size > 0 && !baselineDone {
				n = size
			}
		case "expansion":
			if baselineDone {
				continue
			}
			factor := e.Parameters["expansion_factor"]
			if dur > 0 {
				r := 0.0
				if factor > 0 {
					r = math.Log(factor) / dur
				} else {
					r = e.Parameters["growth_rate"]
				}
				if r == 0 {
					continue
				}
				if t < start+dur {
					n = n * math.Exp(-r*(t-start))
					g = r
					baselineDone = true
					continue
				}
				n = n * math.Exp(-r*dur)
			} else if factor > 0 {
				n = n / factor
			}
		}
	}

	if override > 0 {
		return math.Max(override, 1), 0
	}
	if n < 1 {
		// Clamped sizes are treated as constant to keep the epoch well defined.
		return 1, 0
	}
	return n, g
}

// buildEpochs converts a present-day size and a list of events into a
// piecewise epoch list covering [0, ∞).
func buildEpochs(present float64, events []DemographicEvent) []epoch {
	evs := sortedEvents(events)

	points := map[float64]bool{0: true}
	for _, e := range evs {
		if e.Time < 0 {
			continue
		}
		points[float64(e.Time)] = true
		if e.Duration > 0 {
			points[float64(e.Time+e.Duration)] = true
		}
	}
	starts := make([]float64, 0, len(points))
	for p := range points {
		starts = append(starts, p)
	}
	sort.Float64s(starts)

	epochs := make([]epoch, 0, len(starts))
	for _, a := range starts {
		n, g := sizeAt(present, evs, a)
		if len(epochs) > 0 {
			prev := epochs[len(epochs)-1]
			if prev.g == 0 && g == 0 && prev.n0 == n {
				continue // merge identical constant epochs
			}
		}
		epochs = append(epochs, epoch{start: a, n0: n, g: g})
	}
	// A growing final epoch would have a finite total hazard; it cannot occur
	// because every growth window has a finite end, but guard anyway.
	if last := &epochs[len(epochs)-1]; last.g != 0 {
		last.g = 0
	}
	return epochs
}

// sizeOf returns N(t) for an epoch list.
func sizeOf(epochs []epoch, t float64) float64 {
	e := epochs[epochIndex(epochs, t)]
	return e.n0 * math.Exp(-e.g*(t-e.start))
}

// epochIndex returns the index of the epoch containing t.
func epochIndex(epochs []epoch, t float64) int {
	i := sort.Search(len(epochs), func(i int) bool { return epochs[i].start > t }) - 1
	if i < 0 {
		return 0
	}
	return i
}

// cumHazard is ∫_t^u 1/N(s) ds inside a single epoch.
func (e epoch) cumHazard(t, u float64) float64 {
	if math.IsInf(u, 1) {
		if e.g < 0 {
			// N grows without bound backwards: finite integral.
			return math.Exp(e.g*(t-e.start)) / (e.n0 * -e.g)
		}
		return math.Inf(1)
	}
	if e.g == 0 {
		return (u - t) / e.n0
	}
	return (math.Exp(e.g*(u-e.start)) - math.Exp(e.g*(t-e.start))) / (e.n0 * e.g)
}

// invertHazard returns u such that ∫_t^u 1/N(s) ds = h inside this epoch, or
// +Inf if the epoch can never accumulate that much hazard.
func (e epoch) invertHazard(t, h float64) float64 {
	if e.g == 0 {
		return t + h*e.n0
	}
	arg := math.Exp(e.g*(t-e.start)) + h*e.n0*e.g
	if arg <= 0 {
		return math.Inf(1)
	}
	return e.start + math.Log(arg)/e.g
}

// nextEventTime returns the absolute time of the next event of an
// inhomogeneous Poisson process with intensity rate/N(s), starting at t, given
// a standard exponential draw expDraw. It returns +Inf when the event would
// happen at or after tmax. Sampling is exact (analytic inversion per epoch).
func nextEventTime(epochs []epoch, rate, t, tmax, expDraw float64) float64 {
	if rate <= 0 {
		return math.Inf(1)
	}
	h := expDraw / rate
	for i := epochIndex(epochs, t); i < len(epochs); i++ {
		e := epochs[i]
		end := math.Inf(1)
		if i+1 < len(epochs) {
			end = epochs[i+1].start
		}
		segEnd := math.Min(end, tmax)
		H := e.cumHazard(t, segEnd)
		if h < H {
			u := e.invertHazard(t, h)
			if u >= tmax {
				return math.Inf(1)
			}
			return u
		}
		if segEnd >= tmax {
			return math.Inf(1)
		}
		h -= H
		t = segEnd
	}
	return math.Inf(1)
}
