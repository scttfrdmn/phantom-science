// Package illumina models Illumina short-read output: per-cycle binned base
// qualities, quality-matched substitution errors and paired-end insert
// (fragment) sizes. Every function draws only from the *rand.Rand it is
// given, so output is a pure function of the caller's seed.
//
// # Qualities
//
// Current Illumina instruments (RTA3: NovaSeq 6000/X, NextSeq 1000/2000,
// iSeq) report binned quality scores. The bins used here are the RTA3 bins
// Q2, Q12, Q23 and Q37, written as '#', '-', '8' and 'F' in Phred+33 (see
// Illumina, "Quality Score Encoding" (support.illumina.com, Sequencing
// Analysis Viewer / BaseSpace help) and the technical note "Reducing Whole-
// Genome Data Storage Footprint" (Illumina, 2014), which introduced quality
// binning). A Q2 base is a no-call: Illumina software writes no-calls as N
// with quality 2 ('#'), and ApplyErrors does the same.
//
// The per-read quality string is a first-order Markov chain over the four
// bins. At cycle c the chain keeps its previous bin with probability
// Persistence (so low-quality calls come in short runs, as in real reads),
// and otherwise draws a bin from the cycle's target distribution π(c):
//
//	P(Q37) = max(HighFloor, HighStart − (HighStart − HighEnd)·(c/RefCycles)^Shape)
//
// with the remaining mass split between Q23, Q12 and Q2 in the fixed
// proportions MidShare, LowShare and NoCallShare. Quality therefore falls
// with cycle (faster beyond RefCycles, so 2×300 reads end worse than 2×150
// reads), and read 2, which has the lower HighStart and HighEnd, is worse
// than read 1, as on every Illumina instrument. The shape of this curve and
// its parameters are a model, not a measurement: they are to be calibrated
// against real NovaSeq data (issue #27 follow-up), and the defaults were
// chosen to give a NovaSeq-like ≥ Q30 fraction of about 90 % for 2×150.
//
// # Errors
//
// ApplyErrors substitutes each base with probability 10^(−Q/10) (the Phred
// definition) by one of the three other bases chosen uniformly, and writes
// Q2 bases as N. Insertions and deletions are not modelled; they are rare in
// Illumina data (roughly 10⁻⁵ to 10⁻⁴ per base, much lower than substitutions).
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

// RTA3 quality bins (Phred scores).
const (
	QNoCall = 2  // '#': no-call, written as N
	QLow    = 12 // '-'
	QMid    = 23 // '8'
	QHigh   = 37 // 'F'
)

// Bins lists the quality bins in increasing order.
var Bins = [4]byte{QNoCall, QLow, QMid, QHigh}

// QualityParams parameterises the binned quality chain (see the package
// documentation).
type QualityParams struct {
	// HighStart and HighEnd are P(Q37) at cycle 0 and at cycle RefCycles;
	// HighFloor bounds it from below for long reads.
	HighStart, HighEnd, HighFloor float64
	// RefCycles is the cycle at which P(Q37) reaches HighEnd; Shape is the
	// exponent of the decay (2: slow early, faster late).
	RefCycles, Shape float64
	// MidShare, LowShare and NoCallShare split 1 − P(Q37) between Q23, Q12
	// and Q2; they sum to 1.
	MidShare, LowShare, NoCallShare float64
	// Persistence is the probability that a cycle keeps the previous
	// cycle's bin.
	Persistence float64
}

// Read1 and Read2 are the default parameters for the two mates of a pair
// (Read1 also for single-end reads).
var (
	Read1 = QualityParams{
		HighStart: 0.95, HighEnd: 0.85, HighFloor: 0.30, RefCycles: 150, Shape: 2,
		MidShare: 0.60, LowShare: 0.39, NoCallShare: 0.01, Persistence: 0.6,
	}
	Read2 = QualityParams{
		HighStart: 0.93, HighEnd: 0.76, HighFloor: 0.25, RefCycles: 150, Shape: 2,
		MidShare: 0.60, LowShare: 0.39, NoCallShare: 0.01, Persistence: 0.6,
	}
)

// ParamsFor returns the default quality parameters of mate (1 or 2).
func ParamsFor(mate int) QualityParams {
	if mate == 2 {
		return Read2
	}
	return Read1
}

// high returns P(Q37) at cycle c.
func (p QualityParams) high(c int) float64 {
	x := float64(c) / p.RefCycles
	return max(p.HighFloor, p.HighStart-(p.HighStart-p.HighEnd)*math.Pow(x, p.Shape))
}

// Target returns the cycle's target bin distribution π(c), in the order of
// Bins.
func (p QualityParams) Target(c int) [4]float64 {
	h := p.high(c)
	rest := 1 - h
	return [4]float64{rest * p.NoCallShare, rest * p.LowShare, rest * p.MidShare, h}
}

// draw returns the bin of a uniform u in [0, 1) under π(c).
func (p QualityParams) draw(c int, u float64) byte {
	t := p.Target(c)
	acc := 0.0
	for i, w := range t[:3] {
		acc += w
		if u < acc {
			return Bins[i]
		}
	}
	return QHigh
}

// Qualities returns n binned Phred scores (not ASCII) for one read, drawing
// one value from r per cycle.
func (p QualityParams) Qualities(r *rand.Rand, n int) []byte {
	q := make([]byte, n)
	for c := range q {
		u := r.Float64()
		if c > 0 && u < p.Persistence {
			q[c] = q[c-1]
			continue
		}
		if c > 0 {
			u = (u - p.Persistence) / (1 - p.Persistence)
		}
		q[c] = p.draw(c, u)
	}
	return q
}

// errorProb[q] is the Phred error probability 10^(−q/10).
var errorProb = func() (t [94]float64) {
	for q := range t {
		t[q] = math.Pow(10, -float64(q)/10)
	}
	return t
}()

// ErrorProb returns the substitution probability 10^(−q/10) of Phred q.
func ErrorProb(q byte) float64 { return errorProb[min(int(q), len(errorProb)-1)] }

// otherBases[b] lists the three bases other than b (zero for non-ACGT b).
var otherBases = func() (t [256][3]byte) {
	t['A'], t['C'], t['G'], t['T'] = [3]byte{'C', 'G', 'T'}, [3]byte{'A', 'G', 'T'}, [3]byte{'A', 'C', 'T'}, [3]byte{'A', 'C', 'G'}
	return t
}()

// ApplyErrors applies quality-matched sequencing errors to seq in place and
// returns seq and qual: a base with quality QNoCall (or lower) becomes N,
// and any other A/C/G/T base is replaced, with probability 10^(−Q/10), by
// one of the three other bases chosen uniformly. Other bytes (N) are kept.
// Each non-no-call base draws one value from r, a substitution one more.
func ApplyErrors(r *rand.Rand, seq, qual []byte) ([]byte, []byte) {
	for i, b := range seq {
		if qual[i] <= QNoCall {
			seq[i] = 'N'
			continue
		}
		others := otherBases[b]
		if others[0] == 0 {
			continue
		}
		if r.Float64() < ErrorProb(qual[i]) {
			seq[i] = others[r.IntN(3)]
		}
	}
	return seq, qual
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
