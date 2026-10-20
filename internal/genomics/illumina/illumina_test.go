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

func isBin(q byte) bool {
	for _, b := range Bins {
		if q == b {
			return true
		}
	}
	return false
}

// Only the four RTA3 bins appear; the cycle-0 bin frequencies match π(0)
// within the 99 % binomial CI; the chain reproduces from the seed.
func TestQualities_BinsAndCycleZero(t *testing.T) {
	const reads, n = 20000, 150
	r := newRand(1)
	var first [4]int
	for range reads {
		q := Read1.Qualities(r, n)
		for _, v := range q {
			if !isBin(v) {
				t.Fatalf("quality %d is not an RTA3 bin", v)
			}
		}
		for i, b := range Bins {
			if q[0] == b {
				first[i]++
			}
		}
	}
	for i, p := range Read1.Target(0) {
		got := float64(first[i]) / reads
		ci := z99*math.Sqrt(p*(1-p)/reads) + 1e-9
		if math.Abs(got-p) > ci {
			t.Errorf("cycle 0 P(Q%d) = %.4f, want %.4f ± %.4f", Bins[i], got, p, ci)
		}
	}
	if !bytes.Equal(Read1.Qualities(newRand(7), n), Read1.Qualities(newRand(7), n)) {
		t.Error("same seed gave different qualities")
	}
}

// meanByCycle returns the mean Phred score per cycle over reads draws.
func meanByCycle(p QualityParams, seed uint64, reads, n int) []float64 {
	r := newRand(seed)
	sum := make([]float64, n)
	for range reads {
		for c, v := range p.Qualities(r, n) {
			sum[c] += float64(v)
		}
	}
	for c := range sum {
		sum[c] /= float64(reads)
	}
	return sum
}

func mean(v []float64) float64 {
	s := 0.0
	for _, x := range v {
		s += x
	}
	return s / float64(len(v))
}

// Mean quality falls from early to late cycles, read 2 is worse than read 1,
// and the per-cycle mean tracks the target distribution's mean. The ≥ Q30
// fraction is logged (NovaSeq specifies ≥ 85 % at 2×150).
func TestQualities_DegradeWithCycleAndRead2Lower(t *testing.T) {
	const reads, n = 20000, 150
	m1 := meanByCycle(Read1, 2, reads, n)
	m2 := meanByCycle(Read2, 3, reads, n)
	early, late := mean(m1[:20]), mean(m1[n-20:])
	if !(early > late+1) {
		t.Errorf("R1 mean quality: cycles 1-20 %.2f, last 20 %.2f; want a decrease", early, late)
	}
	if !(mean(m2) < mean(m1)-0.5) {
		t.Errorf("R2 mean %.2f not below R1 mean %.2f", mean(m2), mean(m1))
	}
	for _, c := range []int{0, 75, 149} {
		tg := Read1.Target(c)
		want := 0.0
		for i, b := range Bins {
			want += tg[i] * float64(b)
		}
		// The chain lags π by about 1/(1−Persistence) cycles; π moves by
		// < 0.1 Phred per cycle, so ±0.5 covers lag and sampling noise.
		if math.Abs(m1[c]-want) > 0.5 {
			t.Errorf("cycle %d mean %.2f, target mean %.2f", c, m1[c], want)
		}
	}
	q30 := func(p QualityParams) float64 {
		s := 0.0
		for c := range n {
			s += p.Target(c)[3]
		}
		return s / n
	}
	t.Logf("R1 mean %.2f (early %.2f, late %.2f), R2 mean %.2f; target >=Q30: R1 %.3f, R2 %.3f",
		mean(m1), early, late, mean(m2), q30(Read1), q30(Read2))
}

// Low-quality calls come in runs: P(next is Q12 | Q12) is well above the
// marginal P(Q12).
func TestQualities_LowQualityRuns(t *testing.T) {
	r := newRand(4)
	var low, lowLow, pairs, lowFirst int
	for range 5000 {
		q := Read2.Qualities(r, 150)
		for c := 1; c < len(q); c++ {
			pairs++
			if q[c] == QLow {
				low++
			}
			if q[c-1] == QLow {
				lowFirst++
				if q[c] == QLow {
					lowLow++
				}
			}
		}
	}
	marg, cond := float64(low)/float64(pairs), float64(lowLow)/float64(lowFirst)
	if !(cond > 5*marg) {
		t.Errorf("P(Q12|Q12) = %.3f, marginal %.3f: no runs", cond, marg)
	}
}

// The empirical substitution rate of each bin lies in the 99 % binomial CI
// of 10^(−Q/10); substitutions never keep the base and pick the three other
// bases uniformly; Q2 bases become N.
func TestApplyErrors_RatePerBin(t *testing.T) {
	r := newRand(5)
	const acgt = "ACGT"
	for _, q := range []byte{QLow, QMid, QHigh, 30, 40} {
		n := 2_000_000
		if q >= QHigh {
			n = 10_000_000
		}
		truth := make([]byte, n)
		for i := range truth {
			truth[i] = acgt[r.IntN(4)]
		}
		seq := bytes.Clone(truth)
		qual := bytes.Repeat([]byte{q}, n)
		ApplyErrors(r, seq, qual)
		var subs int
		var to [4]int
		for i := range seq {
			if seq[i] != truth[i] {
				subs++
				to[bytes.IndexByte([]byte(acgt), seq[i])]++
			}
		}
		p := ErrorProb(q)
		got := float64(subs) / float64(n)
		ci := z99 * math.Sqrt(p*(1-p)/float64(n))
		if math.Abs(got-p) > ci {
			t.Errorf("Q%d: substitution rate %.6f, want %.6f ± %.6f", q, got, p, ci)
		} else {
			t.Logf("Q%d: substitution rate %.6f (expected %.6f ± %.6f, n=%d)", q, got, p, ci, n)
		}
		// Uniform target: each base receives about 1/4 of substitutions.
		for i, k := range to {
			e := float64(subs) / 4
			if math.Abs(float64(k)-e) > 4*math.Sqrt(e)+1 {
				t.Errorf("Q%d: %d substitutions to %c, want ≈ %.0f", q, k, acgt[i], e)
			}
		}
	}
	seq := []byte("ACGTACGT")
	qual := []byte{2, 2, 37, 37, 2, 37, 37, 2}
	ApplyErrors(newRand(6), seq, qual)
	for i, q := range qual {
		if (q == QNoCall) != (seq[i] == 'N') {
			t.Errorf("base %d (Q%d) = %c", i, q, seq[i])
		}
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
