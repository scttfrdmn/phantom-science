//go:build !release

package illumina

import (
	"bytes"
	"math"
	"math/rand/v2"
	"testing"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func newRand(seed uint64) *rand.Rand { return rand.New(rand.NewPCG(seed, 0x111)) } // #nosec G404

// z99 is the two-sided 99 % normal quantile.
const z99 = 2.5758

func binIndex(q byte) int {
	for i, b := range Bins {
		if q == b {
			return i
		}
	}
	return -1
}

// The class weights sum to 1 and the weighted class factors average to 1,
// so the class mixture has the mean per-cycle fractions.
func TestQualityParams_ClassesNormalised(t *testing.T) {
	for mate, p := range []QualityParams{Read1, Read2} {
		var w, low, mid float64
		for _, k := range p.Classes {
			w += k.Weight
			low += k.Weight * k.Low
			mid += k.Weight * k.Mid
		}
		for name, v := range map[string]float64{"weights": w, "Q11 factors": low, "Q25 factors": mid} {
			if math.Abs(v-1) > 1e-4 {
				t.Errorf("read %d: %s sum to %.6f, want 1", mate+1, name, v)
			}
		}
		if len(p.Low) != len(p.Knots) || len(p.Mid) != len(p.Knots) {
			t.Errorf("read %d: %d knots, %d/%d values", mate+1, len(p.Knots), len(p.Low), len(p.Mid))
		}
	}
}

// Only the four NovaSeq values appear, for any read length; the chain
// reproduces from the seed.
func TestQualities_AlphabetAndDeterminism(t *testing.T) {
	r := newRand(1)
	for _, n := range []int{0, 1, 2, 50, 151, 250, 300} {
		for range 2000 {
			q := Read2.Qualities(r, n)
			if len(q) != n {
				t.Fatalf("n=%d: %d qualities", n, len(q))
			}
			for _, v := range q {
				if binIndex(v) < 0 {
					t.Fatalf("quality %d is not a NovaSeq value", v)
				}
			}
		}
	}
	if !bytes.Equal(Read1.Qualities(newRand(7), 151), Read1.Qualities(newRand(7), 151)) {
		t.Error("same seed gave different qualities")
	}
}

// genStats are statistics of generated quality strings.
type genStats struct {
	reads     int
	perCycle  [][4]int // bin counts per cycle
	lowCounts []int    // reads by number of bases below Q37
	from      [4]int   // transitions out of each bin
	same      [4]int   // transitions into the same bin
	sum       float64  // sum of Phred values
}

func generate(p QualityParams, seed uint64, reads, n int) genStats {
	r := newRand(seed)
	s := genStats{reads: reads, perCycle: make([][4]int, n), lowCounts: make([]int, n+1)}
	for range reads {
		q := p.Qualities(r, n)
		low := 0
		for c, v := range q {
			i := binIndex(v)
			s.perCycle[c][i]++
			s.sum += float64(v)
			if v != QHigh {
				low++
			}
			if c > 0 {
				j := binIndex(q[c-1])
				s.from[j]++
				if j == i {
					s.same[j]++
				}
			}
		}
		s.lowCounts[low]++
	}
	return s
}

// The generated per-cycle bin frequencies match Target(c, n) within the
// binomial sampling error, for the measured read length and for longer
// (extrapolated) reads; the stickiness does not shift the distribution.
func TestQualities_FrequenciesMatchTarget(t *testing.T) {
	const reads = 60000
	for _, tc := range []struct {
		p QualityParams
		n int
	}{{Read1, 151}, {Read2, 151}, {Read2, 251}} {
		s := generate(tc.p, 2, reads, tc.n)
		bad := 0
		for c := range tc.n {
			want := tc.p.Target(c, tc.n)
			for i := 1; i < 4; i++ {
				got := float64(s.perCycle[c][i]) / reads
				// 4σ per cell (about 450-750 cells) plus 0.001 for the lag
				// of the sticky chain behind a moving target.
				if math.Abs(got-want[i]) > 4*math.Sqrt(want[i]*(1-want[i])/reads)+0.001 {
					bad++
					t.Errorf("n=%d cycle %d: P(Q%d) = %.4f, target %.4f", tc.n, c+1, Bins[i], got, want[i])
				}
			}
			if s := want[0] + want[1] + want[2] + want[3]; math.Abs(s-1) > 1e-12 {
				t.Errorf("cycle %d: target sums to %v", c+1, s)
			}
			if bad > 5 {
				t.FailNow()
			}
		}
	}
}

// The substitution rate of each quality and cycle matches Rate(q, c) times
// the base factor within the 99 % binomial CI; substituted bases follow
// Matrix; Q2 bases become N; other qualities use the Phred rate.
func TestApplyErrors_RateAndSpectrum(t *testing.T) {
	e := NovaSeqErrors
	r := newRand(5)
	const n = 150
	for _, q := range []byte{QLow, QMid, QHigh, 30} {
		for x := range 4 {
			reads := 20000
			if q != QLow {
				reads = 40000
			}
			var subs [n]int
			var to [4]int
			seq := make([]byte, n)
			qual := bytes.Repeat([]byte{q}, n)
			for range reads {
				for i := range seq {
					seq[i] = acgt[x]
				}
				ApplyErrors(r, seq, qual)
				for i, b := range seq {
					if b != acgt[x] {
						subs[i]++
						to[baseIndex[b]]++
					}
				}
			}
			// Rates per block of 30 cycles.
			for lo := 0; lo < n; lo += 30 {
				var k int
				var p float64
				for c := lo; c < lo+30; c++ {
					k += subs[c]
					p += e.Rate(q, c) * e.BaseFactor[x]
				}
				m := float64(30 * reads)
				p /= 30
				got := float64(k) / m
				if ci := z99 * math.Sqrt(p*(1-p)/m); math.Abs(got-p) > ci {
					t.Errorf("Q%d %c cycles %d-%d: rate %.6f, want %.6f ± %.6f", q, acgt[x], lo+1, lo+30, got, p, ci)
				}
			}
			total := 0
			for _, v := range to {
				total += v
			}
			for y, v := range to {
				w := e.Matrix[x][y] * float64(total)
				if math.Abs(float64(v)-w) > 4*math.Sqrt(w+1)+1 {
					t.Errorf("Q%d %c→%c: %d substitutions, want ≈ %.0f", q, acgt[x], acgt[y], v, w)
				}
			}
		}
	}
	if got, want := e.Rate(30, 10), ErrorProb(30); got != want {
		t.Errorf("Q30 rate %v, want Phred %v", got, want)
	}
	if got := SubstitutionRate(QHigh, 0); got != e.High {
		t.Errorf("SubstitutionRate(Q37) = %v, want %v", got, e.High)
	}
	seq := []byte("ACGTACGN")
	qual := []byte{2, 2, 37, 37, 2, 37, 37, 37}
	ApplyErrors(newRand(6), seq, qual)
	for i, q := range qual {
		if (q == QNoCall) != (seq[i] == 'N') && i != 7 {
			t.Errorf("base %d (Q%d) = %c", i, q, seq[i])
		}
	}
	if seq[7] != 'N' {
		t.Errorf("N base changed to %c", seq[7])
	}
}

// The base factors average to 1 over the human base composition and each
// row of the substitution matrix is a distribution over the other bases.
func TestErrorParams_Normalised(t *testing.T) {
	e := NovaSeqErrors
	comp := [4]float64{0.295, 0.205, 0.205, 0.295}
	m := 0.0
	for x, f := range e.BaseFactor {
		m += comp[x] * f
		s := 0.0
		for y, w := range e.Matrix[x] {
			if y == x && w != 0 {
				t.Errorf("Matrix[%c][%c] = %v, want 0", acgt[x], acgt[y], w)
			}
			s += w
		}
		if math.Abs(s-1) > 1e-4 {
			t.Errorf("Matrix row %c sums to %v", acgt[x], s)
		}
	}
	if math.Abs(m-1) > 1e-4 {
		t.Errorf("base factors average to %v over the human composition", m)
	}
}

// Insert sizes: mean and SD within the 99 % CI of the configured values for
// normal and lognormal; draws never fall below the read length or above Max.
func TestInsertSampler_Moments(t *testing.T) {
	for _, tc := range []struct {
		spec    types.InsertSizeSpec
		readLen int
	}{
		{types.InsertSizeSpec{}, 150},
		{types.InsertSizeSpec{Dist: types.InsertSizeLognormal, Mean: 500, SD: 120}, 100},
		{types.InsertSizeSpec{Mean: 250, SD: 30}, 150},
	} {
		s, err := NewInsertSampler(tc.spec)
		if err != nil {
			t.Fatal(err)
		}
		r := newRand(8)
		const n = 100000
		var sum, sq float64
		for range n {
			v := s.Sample(r, tc.readLen)
			if v < tc.readLen {
				t.Fatalf("%v: insert %d below read length %d", s, v, tc.readLen)
			}
			sum += float64(v)
			sq += float64(v) * float64(v)
		}
		m := sum / n
		sd := math.Sqrt(sq/n - m*m)
		want := s.Spec()
		// SE of the mean σ/√n; SE of the SD ≈ σ/√(2n) (normal), wider for the
		// lognormal's heavier tail, so allow 2×. Rounding to integers adds
		// variance 1/12, negligible here.
		if math.Abs(m-want.Mean) > z99*want.SD/math.Sqrt(n) {
			t.Errorf("%v: mean %.2f, want %.2f", s, m, want.Mean)
		}
		if math.Abs(sd-want.SD) > 2*z99*want.SD/math.Sqrt(2*n) {
			t.Errorf("%v: sd %.2f, want %.2f", s, sd, want.SD)
		}
		t.Logf("%v: mean %.2f sd %.2f", s, m, sd)
	}

	// Truncation: mean 250, sd 30 with read length 250 keeps only the upper
	// half, and Max caps.
	s, _ := NewInsertSampler(types.InsertSizeSpec{Mean: 250, SD: 30, Max: 280})
	r := newRand(9)
	for range 10000 {
		if v := s.Sample(r, 250); v < 250 || v > 280 {
			t.Fatalf("insert %d outside [250, 280]", v)
		}
	}
}

func TestInsertSampler_Validation(t *testing.T) {
	for _, spec := range []types.InsertSizeSpec{
		{Dist: "uniform"}, {Mean: -1}, {SD: -5}, {Min: 500, Max: 400}, {Max: 50}, {Min: 2000},
	} {
		if _, err := NewInsertSampler(spec); err == nil {
			t.Errorf("%+v accepted", spec)
		}
	}
}
