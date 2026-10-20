// Package illumina models Illumina short-read output: per-cycle binned base
// qualities, quality-matched substitution errors and paired-end insert
// (fragment) sizes. Every function draws only from the *rand.Rand it is
// given, so output is a pure function of the caller's seed.
//
// The default parameters (calibration.go) are fitted to real data: the
// first 200,000 read pairs (2×151) of ENA run DRR276626, Illumina NovaSeq
// 6000, human whole-genome sequencing, aligned to GRCh38 with minimap2.
// scripts/acceptance/illumina_calibrate.sh downloads the reads, measures
// them and regenerates calibration.go; the measured tables are committed in
// scripts/acceptance/calibrate/illumina/results and the tests compare the
// model with them.
//
// # Qualities
//
// NovaSeq 6000 (RTA3) writes four quality values: Q2, Q11, Q25 and Q37,
// '#', ',', ':' and 'F' in Phred+33. These are the values in DRR276626 and
// in current NovaSeq 6000 output. Illumina's technical note "NovaSeq 6000
// System Quality Scores and RTA3 Software" (770-2017-010, 2017) lists
// 12, 23 and 37 (plus Q2 for no-calls), the values of the first instrument
// software; Illumina told users that NovaSeq Control Software v1.1 changed
// the bins to 2, 11, 25 and 37 without updating that note (reported in
// github.com/benjjneb/dada2 issue #791, 2019). A Q2 base is a no-call,
// written as N.
//
// A read's quality string is generated as follows.
//
//   - The read draws a quality class k (weight w_k). Reads differ a lot:
//     in DRR276626 read 1, 11 % of reads have no base below Q37, the median
//     read has 3, and 5 % have more than 30. The classes multiply the
//     per-cycle fractions of Q11 and Q25 by m_k,11 and m_k,25.
//   - The mean fractions of Q11 and Q25 at cycle c, l_11(c) and l_25(c),
//     are piecewise linear in c between knots (fitted by least squares to
//     cycles 1-150), extended linearly past the last knot, and multiplied
//     by FinalLow and FinalMid at a read's final cycle (the last cycle of a
//     NovaSeq read is about 25 % worse than the one before). So the class
//     k read has P(Q11) = m_k,11·l_11(c) and P(Q25) = m_k,25·l_25(c) at
//     cycle c (their sum capped at 0.98), P(Q2) = NoCall and Q37 otherwise.
//     The weights and multipliers were fitted by EM to the per-read counts
//     of Q11 and Q25 bases; they average to 1, so the mixture reproduces
//     l(c).
//   - Given its class a read's cycles are nearly independent: in the real
//     reads P(low | low k cycles earlier) is the same for every k from 1 to
//     80, so the clustering of low-quality calls comes from the read, not
//     from neighbouring cycles. A small stickiness (a Q11 cycle repeats
//     with probability StickLow, a Q25 cycle with StickMid, with the other
//     draws adjusted so the cycle's distribution is unchanged) matches the
//     measured P(Q11 | Q11) and P(Q25 | Q25) of adjacent cycles.
//
// Read 2 has its own parameters (it has more Q11 and Q25 than read 1, as on
// every Illumina instrument).
//
// # Errors
//
// Real Illumina bins are not exactly their Phred values. The substitution
// rate of each bin was measured as the mismatch rate against GRCh38
// (primary alignments, MAPQ ≥ 30, soft clips, indels and N excluded) less
// the rate of the sample's own variants (0.092 % per base). The variant
// rate is the Q37 mismatch rate less the Q37 error rate where mates overlap
// (a mismatch in one mate where the other mate has the reference base). Q37
// bases have one rate, 3.2·10⁻⁴ (Phred 35, not 37). Q25 and Q11 bases get
// worse along the read: e(c) = late − (late − early)·exp(−(c−1)/τ), Q11
// from 2.9 % at cycle 1 to 8.3 % at cycle 150 (Phred Q11 is 7.9 %), Q25
// from 0.07 % to 0.25 % (Phred Q25 is 0.32 %). The rate also depends on the
// base: C and T are miscalled more often than A and G (BaseFactor), and the
// substituted base follows the measured spectrum (Matrix: C→A, G→T, T→A,
// T→G and A→C dominate; C→T and T→C are rare), both measured from Q11
// mismatches, of which about 99 % are sequencing errors, and applied to
// every bin. Insertions and deletions are not modelled: the measured indel
// error rate is 7·10⁻⁵ per base (25 events where mates overlap, against
// about 1.7·10⁻³ substitutions per base).
//
// # Limitations
//
// The reads come from one tile of one lane of one run and one library;
// other runs, instruments (NovaSeq X uses other bins) and read lengths are
// not measured. Beyond 151 cycles the curves are extrapolated. Qualities do
// not depend on the sequence (in real reads low bins are commoner at T).
// No-calls are placed independently at the measured rate (real no-calls
// come in whole-cycle events: the 11 N of DRR276626's sample are all at
// read 1 cycle 22).
//
// # Insert sizes
//
// InsertSampler draws fragment lengths from a normal or log-normal
// distribution, redrawing values outside [max(Min, read length), Max].
// Fragments shorter than the read (adapter read-through) are not modelled.
package illumina

import (
	"fmt"
	"math"
	"math/rand/v2"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// NovaSeq 6000 (RTA3) quality values (Phred scores).
const (
	QNoCall = 2  // '#': no-call, written as N
	QLow    = 11 // ','
	QMid    = 25 // ':'
	QHigh   = 37 // 'F'
)

// Bins lists the quality values in increasing order.
var Bins = [4]byte{QNoCall, QLow, QMid, QHigh}

// maxLow caps P(Q11) + P(Q25) of a class at a cycle.
const maxLow = 0.98

// ReadClass is a read quality class: a share of reads and the factors by
// which their Q11 and Q25 fractions differ from the mean.
type ReadClass struct {
	Weight   float64
	Low, Mid float64
}

// QualityParams parameterises the binned quality model of one mate (see
// the package documentation).
type QualityParams struct {
	// Knots are 1-based cycles in increasing order; Low and Mid are the
	// mean fractions of Q11 and Q25 at those cycles.
	Knots    []int
	Low, Mid []float64
	// FinalLow and FinalMid multiply the fractions at a read's last cycle.
	FinalLow, FinalMid float64
	// Classes are the read quality classes; their weights sum to 1 and
	// their weighted factors average to 1.
	Classes []ReadClass
	// StickLow and StickMid are the probabilities that a Q11 or Q25 cycle
	// is repeated at the next cycle beyond its share of that cycle.
	StickLow, StickMid float64
	// NoCall is P(Q2) per base.
	NoCall float64
}

// ParamsFor returns the default quality parameters of mate (1 or 2).
func ParamsFor(mate int) QualityParams {
	if mate == 2 {
		return Read2
	}
	return Read1
}

// curve evaluates the piecewise-linear curve vals over p.Knots at 1-based
// cycle x, starting the knot search at *j (cycles must not decrease between
// calls sharing j).
func (p QualityParams) curve(vals []float64, x int, j *int) float64 {
	k := p.Knots
	if x <= k[0] || len(k) == 1 {
		return vals[0]
	}
	for *j < len(k)-2 && x > k[*j+1] {
		*j++
	}
	t := float64(x-k[*j]) / float64(k[*j+1]-k[*j])
	return max(0, vals[*j]+t*(vals[*j+1]-vals[*j]))
}

// means returns the mean Q11 and Q25 fractions at 0-based cycle c of a read
// of n cycles.
func (p QualityParams) means(c, n int, j *int) (low, mid float64) {
	low, mid = p.curve(p.Low, c+1, j), p.curve(p.Mid, c+1, j)
	if c == n-1 {
		low, mid = low*p.FinalLow, mid*p.FinalMid
	}
	return low, mid
}

// classTarget returns class k's distribution over Bins given the mean
// fractions.
func (p QualityParams) classTarget(k ReadClass, low, mid float64) [4]float64 {
	a, b := k.Low*low, k.Mid*mid
	if s := a + b; s > maxLow {
		a, b = a*maxLow/s, b*maxLow/s
	}
	nc := p.NoCall
	return [4]float64{nc, a, b, max(0, 1-nc-a-b)}
}

// Target returns the distribution over Bins at 0-based cycle c of a read of
// n cycles, averaged over the read classes.
func (p QualityParams) Target(c, n int) [4]float64 {
	j := 0
	low, mid := p.means(c, n, &j)
	var t [4]float64
	w := 0.0
	for _, k := range p.Classes {
		ct := p.classTarget(k, low, mid)
		for i := range t {
			t[i] += k.Weight * ct[i]
		}
		w += k.Weight
	}
	for i := range t {
		t[i] /= w
	}
	return t
}

// class returns the class of a uniform u.
func (p QualityParams) class(u float64) ReadClass {
	acc := 0.0
	for _, k := range p.Classes {
		acc += k.Weight
		if u < acc {
			return k
		}
	}
	return p.Classes[len(p.Classes)-1]
}

// pick returns the bin of a uniform u under weights w (which sum to 1).
func pick(w [4]float64, u float64) byte {
	acc := 0.0
	for i, x := range w[:3] {
		acc += x
		if u < acc {
			return Bins[i]
		}
	}
	return QHigh
}

// Qualities returns n binned Phred scores (not ASCII) for one read, drawing
// one value from r for the read's class and one per cycle.
func (p QualityParams) Qualities(r *rand.Rand, n int) []byte {
	q := make([]byte, n)
	if n == 0 {
		return q
	}
	k := p.class(r.Float64())
	j := 0
	for c := range q {
		low, mid := p.means(c, n, &j)
		t := p.classTarget(k, low, mid)
		u := r.Float64()
		if c == 0 {
			q[c] = pick(t, u)
			continue
		}
		// Keep a sticky previous bin with probability rho; otherwise draw
		// from nu(b) = t(b)(1 − rho_b) / (1 − Σ rho t), which leaves the
		// cycle's distribution at t.
		rho := 0.0
		switch q[c-1] {
		case QLow:
			rho = p.StickLow
		case QMid:
			rho = p.StickMid
		}
		if u < rho {
			q[c] = q[c-1]
			continue
		}
		u = (u - rho) / (1 - rho)
		z := 1 - p.StickLow*t[1] - p.StickMid*t[2]
		t[1] *= (1 - p.StickLow) / z
		t[2] *= (1 - p.StickMid) / z
		t[0] /= z
		t[3] /= z
		q[c] = pick(t, u)
	}
	return q
}

// ErrorParams is a substitution model for binned qualities (see the package
// documentation).
type ErrorParams struct {
	// High is the rate of Q37 bases.
	High float64
	// Q25 and Q11 rates at cycle c (1-based):
	// late − (late − early)·exp(−(c−1)/tau).
	MidEarly, MidLate, MidTau float64
	LowEarly, LowLate, LowTau float64
	// BaseFactor multiplies the rate of a base (A, C, G, T); it averages to
	// 1 over the human genome's base composition.
	BaseFactor [4]float64
	// Matrix[x][y] is P(read base y | true base x, substitution), bases in
	// the order A, C, G, T.
	Matrix [4][4]float64
}

// Rate returns the substitution rate of a base with quality q at 0-based
// cycle c, averaged over bases. Qualities other than the three called bins
// get their Phred rate 10^(−q/10).
func (e ErrorParams) Rate(q byte, c int) float64 {
	x := float64(c)
	switch q {
	case QHigh:
		return e.High
	case QMid:
		return e.MidLate - (e.MidLate-e.MidEarly)*math.Exp(-x/e.MidTau)
	case QLow:
		return e.LowLate - (e.LowLate-e.LowEarly)*math.Exp(-x/e.LowTau)
	}
	return ErrorProb(q)
}

// errorProb[q] is the Phred error probability 10^(−q/10).
var errorProb = func() (t [94]float64) {
	for q := range t {
		t[q] = math.Pow(10, -float64(q)/10)
	}
	return t
}()

// ErrorProb returns the Phred error probability 10^(−q/10) of q. The model
// does not use it for the binned qualities, whose measured rates differ
// (see SubstitutionRate).
func ErrorProb(q byte) float64 { return errorProb[min(int(q), len(errorProb)-1)] }

// SubstitutionRate returns the substitution rate that ApplyErrors applies,
// on average over bases, to a base of quality q at 0-based cycle c.
func SubstitutionRate(q byte, c int) float64 { return NovaSeqErrors.Rate(q, c) }

// baseIndex maps A, C, G, T to 0-3 and every other byte to -1.
var baseIndex = func() (t [256]int8) {
	for i := range t {
		t[i] = -1
	}
	t['A'], t['C'], t['G'], t['T'] = 0, 1, 2, 3
	return t
}()

const acgt = "ACGT"

// Apply applies substitution errors to seq in place and returns seq and
// qual: a base with quality QNoCall (or lower) becomes N, and any other
// A/C/G/T base at index (cycle) i is replaced, with probability
// Rate(qual[i], i)·BaseFactor[base], by a base drawn from Matrix. Other
// bytes (N) are kept. Each non-no-call A/C/G/T base draws one value from r,
// a substitution one more.
func (e ErrorParams) Apply(r *rand.Rand, seq, qual []byte) ([]byte, []byte) {
	for i, b := range seq {
		if qual[i] <= QNoCall {
			seq[i] = 'N'
			continue
		}
		x := baseIndex[b]
		if x < 0 {
			continue
		}
		if r.Float64() >= e.Rate(qual[i], i)*e.BaseFactor[x] {
			continue
		}
		u, acc := r.Float64(), 0.0
		to := byte(0)
		for y, w := range e.Matrix[x] {
			if w == 0 {
				continue
			}
			to = acgt[y]
			acc += w
			if u < acc {
				break
			}
		}
		seq[i] = to
	}
	return seq, qual
}

// ApplyErrors applies the NovaSeq substitution model (NovaSeqErrors.Apply)
// to seq in place, treating index i as sequencing cycle i (the engine passes
// reads in sequencing orientation).
func ApplyErrors(r *rand.Rand, seq, qual []byte) ([]byte, []byte) {
	return NovaSeqErrors.Apply(r, seq, qual)
}

// maxRedraws bounds the redraws of an out-of-range insert size.
const maxRedraws = 64

// InsertSampler draws paired-end insert (fragment) sizes.
type InsertSampler struct {
	spec      types.InsertSizeSpec
	mu, sigma float64 // lognormal parameters
}

// NewInsertSampler returns the sampler for spec (defaults applied: normal,
// mean 350, SD 50).
func NewInsertSampler(spec types.InsertSizeSpec) (*InsertSampler, error) {
	full := types.IlluminaSpec{InsertSize: spec}.WithDefaults()
	if err := full.Validate(); err != nil {
		return nil, err
	}
	s := &InsertSampler{spec: full.InsertSize}
	if s.spec.Dist == types.InsertSizeLognormal {
		m, sd := s.spec.Mean, s.spec.SD
		s.sigma = math.Sqrt(math.Log(1 + sd*sd/(m*m)))
		s.mu = math.Log(m) - s.sigma*s.sigma/2
	}
	return s, nil
}

// Spec returns the sampler's distribution with defaults applied.
func (s *InsertSampler) Spec() types.InsertSizeSpec { return s.spec }

// Sample returns an insert size of at least max(Min, readLen) and at most Max
// (when set). Out-of-range draws are redrawn; after maxRedraws draws the
// nearest bound of the mean is returned.
func (s *InsertSampler) Sample(r *rand.Rand, readLen int) int {
	lo := max(s.spec.Min, readLen, 1)
	hi := math.MaxInt32
	if s.spec.Max > 0 {
		hi = max(s.spec.Max, lo)
	}
	for range maxRedraws {
		var v float64
		if s.spec.Dist == types.InsertSizeLognormal {
			v = math.Exp(s.mu + s.sigma*r.NormFloat64())
		} else {
			v = s.spec.Mean + s.spec.SD*r.NormFloat64()
		}
		if n := math.Round(v); n >= float64(lo) && n <= float64(hi) {
			return int(n)
		}
	}
	return min(max(int(math.Round(s.spec.Mean)), lo), hi)
}

// String describes the sampler, e.g. "normal(350, 50)".
func (s *InsertSampler) String() string {
	return fmt.Sprintf("%s(%g, %g)", s.spec.Dist, s.spec.Mean, s.spec.SD)
}
