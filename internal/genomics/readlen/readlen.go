// Package readlen samples synthetic read lengths from a configured
// distribution (#10): fixed, log-normal, gamma or an empirical histogram.
//
// A Distribution is built once from a types.ReadLengthSpec with New and then
// sampled once per read with a caller-supplied generator, so output is a
// pure function of the generator's seed. Fixed distributions never draw from
// the generator.
//
// Parametrisations:
//
//   - fixed: every read is the fixedLength passed to New
//     (SequencerConfig.ReadLength). No other spec field may be set.
//   - lognormal from N50 and Median (ONT style). The length-weighted
//     distribution of LN(μ, σ²) is LN(μ+σ², σ²), so median = e^μ and
//     N50 = e^(μ+σ²); hence μ = ln(median), σ² = ln(N50/median), which
//     requires N50 > median.
//   - lognormal from Mean and SD: σ² = ln(1 + SD²/Mean²), μ = ln(Mean) − σ²/2.
//   - gamma from Shape and Scale, or from Mean and SD
//     (shape = Mean²/SD², scale = SD²/Mean), sampled with Marsaglia–Tsang.
//   - empirical: a length histogram loaded from EmpiricalPath (see
//     LoadHistogram for the file format), sampled by inverse CDF.
//
// Continuous draws are rounded to the nearest integer. Lengths are truncated
// to [Min, Max] (Min defaults to 1, Max 0 means unbounded): a draw outside
// the bounds is rejected and redrawn, up to maxTries times, after which the
// last draw is clamped. New rejects bounds that keep less than
// minAcceptance of the distribution, so the clamp fallback has probability
// below (1 − minAcceptance)^maxTries ≈ 4e-5 and the sampled distribution is
// the truncated one. Empirical histograms are truncated exactly by dropping
// bins outside the bounds.
package readlen

import (
	"errors"
	"fmt"
	"math"
	"math/rand/v2"
	"sort"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Sampler draws one read length (>= 1) per call, using only r.
type Sampler interface {
	Sample(r *rand.Rand) int
}

const (
	// maxTries bounds the rejection loop that enforces Min/Max.
	maxTries = 1000
	// minAcceptance is the smallest fraction of a continuous distribution
	// the Min/Max bounds may keep.
	minAcceptance = 0.01
	// maxLength caps continuous draws before conversion to int.
	maxLength = 1 << 31
)

// ErrInvalidSpec marks an invalid read-length specification. Test with
// errors.Is.
var ErrInvalidSpec = errors.New("invalid read-length distribution")

func invalid(format string, args ...any) error {
	return fmt.Errorf("%w: %s", ErrInvalidSpec, fmt.Sprintf(format, args...))
}

// Distribution is a read-length distribution ready for sampling. It is
// immutable and safe for concurrent use (each caller supplies its own
// generator).
type Distribution struct {
	dist     string
	fixed    int
	min, max int // max 0 = unbounded

	// continuous distributions (lognormal, gamma)
	mu, sigma    float64 // lognormal
	shape, scale float64 // gamma

	// empirical histogram restricted to [min, max]: lengths ascending and
	// the cumulative count up to and including each length.
	lengths []int
	cum     []float64

	mean, sd float64 // of the sampled (truncated, integer) lengths
	accept   float64 // probability one draw falls within the bounds
}

// New validates spec and returns its distribution. fixedLength is the read
// length used when spec.Dist is "" or "fixed" (SequencerConfig.ReadLength);
// it is ignored otherwise.
func New(spec types.ReadLengthSpec, fixedLength int) (*Distribution, error) {
	switch spec.Dist {
	case "", types.ReadLengthFixed:
		return newFixed(spec, fixedLength)
	case types.ReadLengthLognormal, types.ReadLengthGamma, types.ReadLengthEmpirical:
	default:
		return nil, fmt.Errorf("unknown read-length distribution %q: want fixed, lognormal, gamma or empirical", spec.Dist)
	}

	d := &Distribution{dist: spec.Dist, min: spec.Min, max: spec.Max}
	switch {
	case spec.Min < 0:
		return nil, invalid("min %d is negative", spec.Min)
	case spec.Max < 0:
		return nil, invalid("max %d is negative", spec.Max)
	case spec.Min == 0:
		d.min = 1
	}
	if d.max != 0 && d.max < d.min {
		return nil, invalid("max %d is less than min %d", d.max, d.min)
	}

	var err error
	switch spec.Dist {
	case types.ReadLengthLognormal:
		err = d.initLognormal(spec)
	case types.ReadLengthGamma:
		err = d.initGamma(spec)
	case types.ReadLengthEmpirical:
		err = d.initEmpirical(spec)
	}
	if err != nil {
		return nil, err
	}
	return d, nil
}

func newFixed(spec types.ReadLengthSpec, n int) (*Distribution, error) {
	if spec.Mean != 0 || spec.SD != 0 || spec.N50 != 0 || spec.Median != 0 || spec.Shape != 0 ||
		spec.Scale != 0 || spec.EmpiricalPath != "" || spec.Min != 0 || spec.Max != 0 {
		return nil, invalid("distribution parameters are set but the distribution is fixed; set dist to lognormal, gamma or empirical")
	}
	if n < 1 {
		return nil, invalid("fixed read length %d must be at least 1", n)
	}
	return &Distribution{dist: types.ReadLengthFixed, fixed: n, min: n, max: n, mean: float64(n), accept: 1}, nil
}

func (d *Distribution) initLognormal(spec types.ReadLengthSpec) error {
	byN50 := spec.N50 != 0 || spec.Median != 0
	byMoments := spec.Mean != 0 || spec.SD != 0
	switch {
	case spec.Shape != 0 || spec.Scale != 0 || spec.EmpiricalPath != "":
		return invalid("lognormal takes n50+median or mean+sd, not shape, scale or empirical_path")
	case byN50 && byMoments:
		return invalid("lognormal takes either n50+median or mean+sd, not both")
	case byN50:
		if spec.Median <= 0 || spec.N50 <= 0 {
			return invalid("lognormal needs both n50 and median > 0 (got n50 %d, median %d)", spec.N50, spec.Median)
		}
		if spec.N50 <= spec.Median {
			return invalid("lognormal n50 %d must exceed median %d (N50 > median for every log-normal)", spec.N50, spec.Median)
		}
		s2 := math.Log(float64(spec.N50) / float64(spec.Median))
		d.mu, d.sigma = math.Log(float64(spec.Median)), math.Sqrt(s2)
	case byMoments:
		if !(spec.Mean > 0) || !(spec.SD > 0) || math.IsInf(spec.Mean, 0) || math.IsInf(spec.SD, 0) {
			return invalid("lognormal needs mean > 0 and sd > 0 (got mean %g, sd %g)", spec.Mean, spec.SD)
		}
		s2 := math.Log1p(spec.SD * spec.SD / (spec.Mean * spec.Mean))
		d.mu, d.sigma = math.Log(spec.Mean)-s2/2, math.Sqrt(s2)
	default:
		return invalid("lognormal needs n50+median or mean+sd")
	}
	return d.initContinuous()
}

func (d *Distribution) initGamma(spec types.ReadLengthSpec) error {
	byShape := spec.Shape != 0 || spec.Scale != 0
	byMoments := spec.Mean != 0 || spec.SD != 0
	switch {
	case spec.N50 != 0 || spec.Median != 0 || spec.EmpiricalPath != "":
		return invalid("gamma takes shape+scale or mean+sd, not n50, median or empirical_path")
	case byShape && byMoments:
		return invalid("gamma takes either shape+scale or mean+sd, not both")
	case byShape:
		if !(spec.Shape > 0) || !(spec.Scale > 0) || math.IsInf(spec.Shape, 0) || math.IsInf(spec.Scale, 0) {
			return invalid("gamma needs shape > 0 and scale > 0 (got shape %g, scale %g)", spec.Shape, spec.Scale)
		}
		d.shape, d.scale = spec.Shape, spec.Scale
	case byMoments:
		if !(spec.Mean > 0) || !(spec.SD > 0) || math.IsInf(spec.Mean, 0) || math.IsInf(spec.SD, 0) {
			return invalid("gamma needs mean > 0 and sd > 0 (got mean %g, sd %g)", spec.Mean, spec.SD)
		}
		d.shape, d.scale = (spec.Mean*spec.Mean)/(spec.SD*spec.SD), spec.SD*spec.SD/spec.Mean
	default:
		return invalid("gamma needs shape+scale or mean+sd")
	}
	return d.initContinuous()
}

// initContinuous computes the acceptance probability and the moments of the
// truncated, rounded distribution. A rounded draw k is accepted when
// min <= k <= max, i.e. the continuous draw lies in [min−½, max+½).
func (d *Distribution) initContinuous() error {
	a, b := float64(d.min)-0.5, math.Inf(1)
	if d.max != 0 {
		b = float64(d.max) + 0.5
	}
	d.accept = d.cdf(b) - d.cdf(a)
	if !(d.accept >= minAcceptance) {
		return invalid("min %d / max %d keep only %.3g%% of the %s distribution (need at least %g%%); widen the bounds",
			d.min, d.max, float64(100*d.accept), d.dist, 100*minAcceptance)
	}
	m1, m2 := d.partialMoments(a, b)
	d.mean = m1 / d.accept
	// Sheppard's correction: rounding to integers adds about 1/12 to the
	// variance of a smooth distribution.
	d.sd = math.Sqrt(math.Max(0, m2/d.accept-float64(d.mean*d.mean)) + 1.0/12)
	return nil
}

// partialMoments returns E[X·1{a<X<b}] and E[X²·1{a<X<b}] of the
// untruncated continuous distribution.
func (d *Distribution) partialMoments(a, b float64) (m1, m2 float64) {
	switch d.dist {
	case types.ReadLengthLognormal:
		mu, s := d.mu, d.sigma
		// E[X^k 1{a<X<b}] = e^(kμ + k²σ²/2) [Φ((ln b − μ − kσ²)/σ) − Φ((ln a − μ − kσ²)/σ)]
		part := func(k float64) float64 {
			lo := normCDF((math.Log(a) - mu - float64(k*s*s)) / s)
			hi := 1.0
			if !math.IsInf(b, 1) {
				hi = normCDF((math.Log(b) - mu - float64(k*s*s)) / s)
			}
			return float64(math.Exp(float64(k*mu)+k*k*s*s/2) * (hi - lo))
		}
		return part(1), part(2)
	default: // gamma
		k, th := d.shape, d.scale
		// E[X 1{..}] = kθ [P(k+1, b/θ) − P(k+1, a/θ)]; E[X² 1{..}] = k(k+1)θ² [P(k+2, ·)]
		part := func(j float64) float64 {
			hi := 1.0
			if !math.IsInf(b, 1) {
				hi = gammaP(k+j, b/th)
			}
			return hi - gammaP(k+j, a/th)
		}
		return float64(k * th * part(1)), float64(k * (k + 1) * th * th * part(2))
	}
}

// cdf is the CDF of the untruncated continuous distribution.
func (d *Distribution) cdf(x float64) float64 {
	if x <= 0 {
		return 0
	}
	if math.IsInf(x, 1) {
		return 1
	}
	switch d.dist {
	case types.ReadLengthLognormal:
		return normCDF((math.Log(x) - d.mu) / d.sigma)
	case types.ReadLengthGamma:
		return gammaP(d.shape, x/d.scale)
	}
	return math.NaN()
}

func (d *Distribution) initEmpirical(spec types.ReadLengthSpec) error {
	if spec.Mean != 0 || spec.SD != 0 || spec.N50 != 0 || spec.Median != 0 || spec.Shape != 0 || spec.Scale != 0 {
		return invalid("empirical takes only empirical_path (and min/max)")
	}
	if spec.EmpiricalPath == "" {
		return invalid("empirical needs empirical_path, a length histogram file")
	}
	h, err := LoadHistogram(spec.EmpiricalPath)
	if err != nil {
		return err
	}
	var total, s1, s2 float64
	for _, l := range h.Lengths() {
		c := h[l]
		if l < d.min || (d.max != 0 && l > d.max) || c == 0 {
			continue
		}
		total += c
		d.lengths = append(d.lengths, l)
		d.cum = append(d.cum, total)
		s1 += float64(c * float64(l))
		s2 += float64(c * float64(l) * float64(l))
	}
	if total == 0 {
		return invalid("histogram %s has no counts within min %d / max %d", spec.EmpiricalPath, d.min, d.max)
	}
	d.accept = 1
	d.mean = s1 / total
	d.sd = math.Sqrt(math.Max(0, s2/total-float64(d.mean*d.mean)))
	return nil
}

// Sample draws one read length using only r. Fixed distributions do not
// draw from r.
func (d *Distribution) Sample(r *rand.Rand) int {
	switch d.dist {
	case types.ReadLengthFixed:
		return d.fixed
	case types.ReadLengthEmpirical:
		u := r.Float64() * d.cum[len(d.cum)-1]
		i := sort.Search(len(d.cum), func(i int) bool { return d.cum[i] > u })
		return d.lengths[min(i, len(d.lengths)-1)]
	}
	var n int
	for range maxTries {
		n = d.draw(r)
		if n >= d.min && (d.max == 0 || n <= d.max) {
			return n
		}
	}
	n = max(n, d.min)
	if d.max != 0 {
		n = min(n, d.max)
	}
	return n
}

// draw returns one untruncated continuous draw rounded to an integer.
func (d *Distribution) draw(r *rand.Rand) int {
	var x float64
	if d.dist == types.ReadLengthLognormal {
		x = math.Exp(d.mu + float64(d.sigma*r.NormFloat64()))
	} else {
		x = d.scale * sampleGamma(r, d.shape)
	}
	return int(math.Min(math.Round(x), maxLength))
}

// sampleGamma draws from Gamma(shape, 1) by Marsaglia & Tsang (2000), "A
// simple method for generating gamma variables". For shape < 1 it draws
// Gamma(shape+1) and multiplies by U^(1/shape).
func sampleGamma(r *rand.Rand, shape float64) float64 {
	if shape < 1 {
		u := r.Float64()
		return float64(sampleGamma(r, shape+1) * math.Pow(u, 1/shape))
	}
	d := shape - 1.0/3
	c := 1 / math.Sqrt(9*d)
	for {
		x := r.NormFloat64()
		v := 1 + float64(c*x)
		if v <= 0 {
			continue
		}
		v = float64(v * v * v)
		u := r.Float64()
		x2 := x * x
		if u < 1-float64(0.0331*x2*x2) {
			return float64(d * v)
		}
		if math.Log(u) < float64(0.5*x2)+float64(d*(1-v+math.Log(v))) {
			return float64(d * v)
		}
	}
}

// Dist returns the distribution name (fixed, lognormal, gamma or empirical).
func (d *Distribution) Dist() string { return d.dist }

// Fixed reports whether every sample has the same length.
func (d *Distribution) Fixed() bool { return d.dist == types.ReadLengthFixed }

// Bounds returns the inclusive length bounds; max 0 means unbounded.
func (d *Distribution) Bounds() (minLen, maxLen int) { return d.min, d.max }

// Mean returns the expected sampled length (after truncation and rounding).
func (d *Distribution) Mean() float64 { return d.mean }

// SD returns the standard deviation of the sampled length (after
// truncation; for continuous distributions including the rounding
// variance of 1/12).
func (d *Distribution) SD() float64 { return d.sd }

// CDF returns P(L <= k) for a sampled length L. For continuous
// distributions the rounded, truncated draw gives
// P(L <= k) = (F(k+½) − F(min−½)) / (F(max+½) − F(min−½)); for empirical
// histograms it is the cumulative bin weight.
func (d *Distribution) CDF(k int) float64 {
	if k < d.min {
		return 0
	}
	if d.max != 0 && k >= d.max {
		return 1
	}
	switch d.dist {
	case types.ReadLengthFixed:
		return 1 // k >= min = fixed
	case types.ReadLengthEmpirical:
		i := sort.SearchInts(d.lengths, k+1) // bins with length <= k
		if i == 0 {
			return 0
		}
		return d.cum[i-1] / d.cum[len(d.cum)-1]
	}
	lo := d.cdf(float64(d.min) - 0.5)
	return math.Min(1, (d.cdf(float64(k)+0.5)-lo)/d.accept)
}

// normCDF is the standard normal CDF.
func normCDF(z float64) float64 { return float64(0.5 * math.Erfc(-z/math.Sqrt2)) }

// gammaP is the regularized lower incomplete gamma function P(a, x)
// (Numerical Recipes §6.2: series for x < a+1, continued fraction otherwise).
func gammaP(a, x float64) float64 {
	if x <= 0 {
		return 0
	}
	if math.IsInf(x, 1) {
		return 1
	}
	lg, _ := math.Lgamma(a)
	front := math.Exp(float64(a*math.Log(x)) - x - lg)
	if x < a+1 {
		sum, term := 1/a, 1/a
		for n := 1; n < 10000; n++ {
			term = float64(term * (x / (a + float64(n))))
			sum += term
			if math.Abs(term) < math.Abs(sum)*1e-15 {
				break
			}
		}
		return math.Min(1, front*sum)
	}
	// Lentz's method for the continued fraction of Q(a, x).
	const tiny = 1e-300
	b := x + 1 - a
	c := 1 / tiny
	dd := 1 / b
	h := dd
	for i := 1; i < 10000; i++ {
		an := -float64(i) * (float64(i) - a)
		b += 2
		dd = float64(an*dd) + b
		if math.Abs(dd) < tiny {
			dd = tiny
		}
		c = b + an/c
		if math.Abs(c) < tiny {
			c = tiny
		}
		dd = 1 / dd
		del := float64(dd * c)
		h *= del
		if math.Abs(del-1) < 1e-15 {
			break
		}
	}
	return math.Max(0, 1-float64(front*h))
}
