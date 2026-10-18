//go:build !release

package quality

import (
	"bytes"
	"reflect"
	"testing"
)

// TestSimulateErrors_DeterministicGivenSeed checks that two models seeded
// alike produce identical errors. chooseErrorType used to walk the error-type
// map in Go's randomised iteration order, so the same draw picked different
// error types from run to run.
//
// The sequence is short with a homopolymer-free tail: SimulateErrors indexes
// later runs by their original positions after earlier deletions (known
// drift, out of scope here), which can run past the end of a long read.
func TestSimulateErrors_DeterministicGivenSeed(t *testing.T) {
	seq := []byte("AAAAACCCGTTTTTTGGGGACGTACGTACGT")
	qual := bytes.Repeat([]byte{30}, len(seq))
	simulate := func(platform string) []*ErrorSimulationResult {
		m := NewHomopolymerErrorModel(platform)
		m.SetSeed(42)
		var out []*ErrorSimulationResult
		for i := 0; i < 300; i++ {
			res, err := m.SimulateErrors(seq, qual)
			if err != nil {
				t.Fatal(err)
			}
			out = append(out, res)
		}
		return out
	}
	for _, platform := range []string{"nanopore", "pacbio", "illumina"} {
		first := simulate(platform)
		types := map[string]bool{}
		for _, res := range first {
			for _, e := range res.Errors {
				types[e.ErrorType] = true
			}
		}
		if len(types) < 2 {
			t.Fatalf("%s: only error types %v simulated; the test needs several", platform, types)
		}
		for run := 1; run < 10; run++ {
			again := simulate(platform)
			for i := range again {
				a, b := again[i], first[i]
				if !bytes.Equal(a.ErrorSequence, b.ErrorSequence) ||
					!bytes.Equal(a.ErrorQualities, b.ErrorQualities) ||
					!reflect.DeepEqual(a.Errors, b.Errors) {
					t.Fatalf("%s: run %d call %d differs from run 0 with the same seed", platform, run, i)
				}
			}
		}
	}
}

// TestChooseErrorType_FollowsDistribution checks that the deterministic
// order still draws each type with its configured probability.
func TestChooseErrorType_FollowsDistribution(t *testing.T) {
	m := NewHomopolymerErrorModel("nanopore")
	m.SetSeed(7)
	const n = 100_000
	got := map[string]int{}
	for i := 0; i < n; i++ {
		got[m.chooseErrorType()]++
	}
	for typ, p := range m.errorDistribution {
		rate := float64(got[typ]) / n
		if d := rate - p; d > 0.01 || d < -0.01 { // > 6σ for p ≤ 0.7
			t.Errorf("%s: rate %.4f, want %.4f", typ, rate, p)
		}
	}
}
