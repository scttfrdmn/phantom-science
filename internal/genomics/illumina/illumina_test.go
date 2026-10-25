//go:build !release

package illumina

import (
	"bytes"
	"math"
	"math/rand/v2"
	"strings"
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

// The default profile is novaseq6000-pooled (#36) and is listed first; the
// package's default models are its; the linear model's substitution model
// (ApplyErrors) stays novaseq6000-drr276626's, so linear output does not
// change with the default.
func TestDefaultProfile(t *testing.T) {
	if DefaultProfile != "novaseq6000-pooled" || ProfileNames()[0] != DefaultProfile {
		t.Fatalf("default %q, names %v", DefaultProfile, ProfileNames())
	}
	p, err := LookupProfile("")
	if err != nil || p != &profilePOOLED {
		t.Fatalf("LookupProfile(\"\") = %v, %v", p, err)
	}
	if NovaSeqErrors != profilePOOLED.Errors || &Read1.Classes[0] != &profilePOOLED.Read1.Classes[0] ||
		&Read2.Classes[0] != &profilePOOLED.Read2.Classes[0] || &NovaSeqPairs.Joint[0] != &profilePOOLED.Pairs.Joint[0] {
		t.Error("the package's default models are not the default profile's")
	}
	if linearErrors != profileDRR276626.Errors {
		t.Error("ApplyErrors is not novaseq6000-drr276626's substitution model")
	}
}

// The pooled profile (#36) resolves by name, case-insensitively; it
// describes the pool as the join of its runs, fitted to equal numbers of
// pairs of each; its qualities, for the read lengths of its runs (2×151 and
// 2×150), take only the NovaSeq values and reproduce from the seed.
func TestPooledProfile(t *testing.T) {
	p, err := LookupProfile("NovaSeq6000-Pooled")
	if err != nil || p != &profilePOOLED {
		t.Fatalf("LookupProfile(novaseq6000-pooled) = %v, %v", p, err)
	}
	if len(p.Runs) != 3 {
		t.Fatalf("%d runs, want 3", len(p.Runs))
	}
	var accs, studies []string
	pairs, maxLen := 0, 0
	for _, r := range p.Runs {
		if r.Pairs != p.Runs[0].Pairs {
			t.Errorf("%s: %d pairs, %s %d: the runs are not pooled equally", r.Accession, r.Pairs, p.Runs[0].Accession, p.Runs[0].Pairs)
		}
		if r.Instrument != "Illumina NovaSeq 6000" {
			t.Errorf("%s: instrument %q", r.Accession, r.Instrument)
		}
		accs, studies = append(accs, r.Accession), append(studies, r.Study)
		pairs += r.Pairs
		maxLen = max(maxLen, r.ReadLength)
	}
	if got := strings.Join(accs, ","); p.Run.Accession != got || p.Run.Study != strings.Join(studies, ",") {
		t.Errorf("pool %q / %q, runs %q / %q", p.Run.Accession, p.Run.Study, got, strings.Join(studies, ","))
	}
	if p.Run.Pairs != pairs || p.Run.ReadLength != maxLen {
		t.Errorf("pool %d pairs 2x%d, runs %d pairs, longest %d", p.Run.Pairs, p.Run.ReadLength, pairs, maxLen)
	}
	for _, m := range []int{1, 2} {
		q := p.Params(m)
		for _, n := range []int{150, 151} {
			r := newRand(3)
			for range 2000 {
				for _, v := range q.Qualities(r, n) {
					if binIndex(v) < 0 {
						t.Fatalf("read %d, %d cycles: quality %d is not a NovaSeq value", m, n, v)
					}
				}
			}
			if !bytes.Equal(q.Qualities(newRand(9), n), q.Qualities(newRand(9), n)) {
				t.Errorf("read %d, %d cycles: same seed gave different qualities", m, n)
			}
		}
	}
}

// The class weights sum to 1 and the weighted class factors average to 1,
// so the class mixture has the mean per-cycle fractions (every profile).
func TestQualityParams_ClassesNormalised(t *testing.T) {
	for _, pr := range profiles {
		testClassesNormalised(t, pr)
	}
}

func testClassesNormalised(t *testing.T, pr *Profile) {
	for mate, p := range []QualityParams{pr.Read1, pr.Read2} {
		var w, low, mid float64
		for _, k := range p.Classes {
			w += k.Weight
			low += k.Weight * k.Low
			mid += k.Weight * k.Mid
		}
		for name, v := range map[string]float64{"weights": w, "Q11 factors": low, "Q25 factors": mid} {
			if math.Abs(v-1) > 1e-4 {
				t.Errorf("%s read %d: %s sum to %.6f, want 1", pr.Name, mate+1, name, v)
			}
		}
		if len(p.Low) != len(p.Knots) || len(p.Mid) != len(p.Knots) {
			t.Errorf("%s read %d: %d knots, %d/%d values", pr.Name, mate+1, len(p.Knots), len(p.Low), len(p.Mid))
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

// generate draws reads quality strings of n cycles: without a sequence
// (Qualities) when bases is nil, else for reads drawn by bases
// (QualitiesFor).
func generate(p QualityParams, seed uint64, reads, n int, bases func(r *rand.Rand, n int) []byte) genStats {
	r := newRand(seed)
	sr := newRand(seed + 1000)
	s := genStats{reads: reads, perCycle: make([][4]int, n), lowCounts: make([]int, n+1)}
	for range reads {
		var q []byte
		if bases == nil {
			q = p.Qualities(r, n)
		} else {
			q = p.QualitiesFor(r, bases(sr, n))
		}
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
		s := generate(tc.p, 2, reads, tc.n, nil)
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
// Matrix; Q2 bases become N; other qualities use the Phred rate. ApplyErrors
// is the linear model's substitution model, novaseq6000-drr276626's.
func TestApplyErrors_RateAndSpectrum(t *testing.T) {
	e := profileDRR276626.Errors
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
// row of the substitution matrix is a distribution over the other bases
// (every profile).
func TestErrorParams_Normalised(t *testing.T) {
	for _, pr := range profiles {
		testErrorsNormalised(t, pr.Name, pr.Errors)
	}
}

func testErrorsNormalised(t *testing.T, name string, e ErrorParams) {
	comp := [4]float64{0.295, 0.205, 0.205, 0.295}
	m := 0.0
	for x, f := range e.BaseFactor {
		m += comp[x] * f
		s := 0.0
		for y, w := range e.Matrix[x] {
			if y == x && w != 0 {
				t.Errorf("%s: Matrix[%c][%c] = %v, want 0", name, acgt[x], acgt[y], w)
			}
			s += w
		}
		if math.Abs(s-1) > 1e-4 {
			t.Errorf("%s: Matrix row %c sums to %v", name, acgt[x], s)
		}
	}
	if math.Abs(m-1) > 1e-4 {
		t.Errorf("%s: base factors average to %v over the human composition", name, m)
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
