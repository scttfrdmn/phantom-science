//go:build !release

package fingerprinting

import (
	"math/rand/v2"
	"testing"
)

// noisySignal returns n samples of a stepped level with Gaussian-ish noise,
// like synthetic nanopore signal.
func noisySignal(r *rand.Rand, n int) []int16 {
	s := make([]int16, n)
	level := 500
	for i := range s {
		if r.IntN(10) == 0 {
			level = 300 + r.IntN(600)
		}
		s[i] = int16(level + r.IntN(13) - 6) // #nosec G115 -- small values
	}
	return s
}

// The rolling window hash equals the direct hash at every position.
func TestSignalWindowRollingHash(t *testing.T) {
	r := rand.New(rand.NewPCG(1, 1))
	raw := noisySignal(r, 500)
	m := NewSignalMatcher()
	h := windowHash(raw[:SignalWindow])
	for i := SignalWindow; i < len(raw); i++ {
		h = h*signalPrime + sample(raw[i]) - sample(raw[i-SignalWindow])*m.pow
		if want := windowHash(raw[i-SignalWindow+1 : i+1]); h != want {
			t.Fatalf("position %d: rolling %x, direct %x", i, h, want)
		}
	}
}

func TestSignalMatcher(t *testing.T) {
	r := rand.New(rand.NewPCG(2, 2))
	var regen [][]int16
	for range 20 {
		regen = append(regen, noisySignal(r, 2000+r.IntN(2000)))
	}
	m := NewSignalMatcher()
	whole := m.AddFile("whole")
	for _, s := range regen[:10] {
		m.AddRead(whole, s)
	}
	trimmed := m.AddFile("trimmed")
	for _, s := range regen[10:] {
		m.AddRead(trimmed, s[150:len(s)-77])
	}
	other := m.AddFile("other")
	for range 10 {
		m.AddRead(other, noisySignal(r, 3000))
	}
	// A suspect read that shares a short stretch (fewer than half its
	// windows) with a regenerated read does not match.
	part := m.AddFile("partial")
	s := append(append([]int16(nil), noisySignal(r, 3000)...), regen[0][:200]...)
	m.AddRead(part, s)

	for _, s := range regen {
		m.Observe(s)
	}
	got := m.Results()
	for _, c := range []struct {
		i, exact, contained int
	}{{whole, 10, 0}, {trimmed, 0, 10}, {other, 0, 0}, {part, 0, 0}} {
		if f := got[c.i]; c.i == trimmed && f.Whole != 10 {
			t.Errorf("trimmed: %d whole, want 10", f.Whole)
		}
		if f := got[c.i]; f.Exact != c.exact || f.Contained != c.contained || f.ExpectedChance != 0 {
			t.Errorf("%s: %+v, want %d exact, %d contained", f.Name, f, c.exact, c.contained)
		}
	}
	if n, samples := m.Regenerated(); n != 20 || samples == 0 {
		t.Errorf("regenerated %d reads, %d samples", n, samples)
	}
	m.Reset()
	if f := m.Results()[whole]; f.Matched() != 0 || f.Checked != 10 {
		t.Errorf("after Reset: %+v", f)
	}
}
