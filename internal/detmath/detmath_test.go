//go:build !release

package detmath

import (
	"math"
	"math/rand/v2"
	"runtime"
	"runtime/debug"
	"testing"
)

// ulps returns the distance between a and b in units in the last place, or
// 0 when both are the same NaN class, or MaxUint64 when exactly one is NaN.
func ulps(a, b float64) uint64 {
	if math.IsNaN(a) || math.IsNaN(b) {
		if math.IsNaN(a) && math.IsNaN(b) {
			return 0
		}
		return math.MaxUint64
	}
	if a == b {
		return 0 // includes +0 == -0
	}
	// Map the floats onto a monotonic integer line.
	ord := func(x float64) int64 {
		b := int64(math.Float64bits(x))
		if b < 0 {
			b = math.MinInt64 - b
		}
		return b
	}
	d := ord(a) - ord(b)
	if d < 0 {
		d = -d
	}
	return uint64(d)
}

// sameBits reports whether a and b are bit-identical, treating all NaNs as
// equal.
func sameBits(a, b float64) bool {
	if math.IsNaN(a) && math.IsNaN(b) {
		return true
	}
	return math.Float64bits(a) == math.Float64bits(b)
}

// anyFinite returns a float64 with uniformly random sign, exponent and
// mantissa (every finite magnitude, including subnormals, is equally likely
// per binade).
func anyFinite(r *rand.Rand) float64 {
	for {
		x := math.Float64frombits(r.Uint64())
		if !math.IsNaN(x) && !math.IsInf(x, 0) {
			return x
		}
	}
}

func uniform(r *rand.Rand, lo, hi float64) float64 {
	return lo + float64(r.Float64()*(hi-lo))
}

// logUniform returns a positive value with log10 uniform in [lo, hi].
func logUniform(r *rand.Rand, lo, hi float64) float64 {
	return math.Pow(10, uniform(r, lo, hi))
}

var specials = []float64{
	0, math.Copysign(0, -1), 1, -1, 0.5, -0.5, 2, -2, 3, -3, 10, -10,
	math.Inf(1), math.Inf(-1), math.NaN(),
	math.SmallestNonzeroFloat64, -math.SmallestNonzeroFloat64,
	0x1p-1022, -0x1p-1022, 0x1p-1023, 0x1p-1074 * 3,
	math.MaxFloat64, -math.MaxFloat64,
	0x1p-28, -0x1p-28, 0x1p-29, 0x1p-54, -0x1p-54,
	709.782712893383973096, 709.79, -745.13321910194110842, -745.14, -708.39, -708.4,
	1023.99, 1024, -1074, -1075, 1e-300, 1e300, -1e300,
	1 + 0x1p-52, 1 - 0x1p-53, math.Sqrt2 / 2, math.Ln2, math.Ln2 / 2, math.Pi, math.Pi / 2, math.Pi / 4,
	0x1p29, 0x1p29 + 1, 1e9, 1e22, -1e22,
	0.84375, -0.84375, 1.25, -1.25, 6, -6, 28, -28, 1 / 0.35,
	0.9, 1.1, 1.9, 2.1, 7.9, 8, 0x1p58, 0x1p59, -0.5 + 0x1p-40, -2.5, -3.5, -100.5, -1e15 - 0.5,
}

type fn1 struct {
	name string
	got  func(float64) float64
	std  func(float64) float64
	gen  func(*rand.Rand) float64
}

func lgammaV(f func(float64) (float64, int)) func(float64) float64 {
	return func(x float64) float64 { v, _ := f(x); return v }
}

var unary = []fn1{
	{"Exp", Exp, math.Exp, func(r *rand.Rand) float64 { return uniform(r, -746, 710) }},
	{"Exp2", Exp2, math.Exp2, func(r *rand.Rand) float64 { return uniform(r, -1076, 1025) }},
	{"Expm1", Expm1, math.Expm1, func(r *rand.Rand) float64 { return uniform(r, -40, 710) }},
	{"Log", Log, math.Log, func(r *rand.Rand) float64 { return logUniform(r, -320, 308) }},
	{"Log2", Log2, math.Log2, func(r *rand.Rand) float64 { return logUniform(r, -320, 308) }},
	{"Log10", Log10, math.Log10, func(r *rand.Rand) float64 { return logUniform(r, -320, 308) }},
	{"Log1p", Log1p, math.Log1p, func(r *rand.Rand) float64 {
		if r.IntN(2) == 0 {
			return -logUniform(r, -320, 0) // (-1, 0)
		}
		return logUniform(r, -320, 308)
	}},
	{"Erf", Erf, math.Erf, func(r *rand.Rand) float64 { return uniform(r, -7, 7) }},
	{"Erfc", Erfc, math.Erfc, func(r *rand.Rand) float64 { return uniform(r, -7, 28) }},
	{"Lgamma", lgammaV(Lgamma), lgammaV(math.Lgamma), func(r *rand.Rand) float64 {
		if r.IntN(4) == 0 {
			return uniform(r, -60, 0)
		}
		return logUniform(r, -10, 10)
	}},
	{"Sin", Sin, math.Sin, func(r *rand.Rand) float64 { return uniform(r, -1e4, 1e4) }},
	{"Cos", Cos, math.Cos, func(r *rand.Rand) float64 { return uniform(r, -1e4, 1e4) }},
}

// bounds are the largest differences from the math package accepted for
// each function, in ulps of the result. The algorithms are the same; the
// math package differs only by fusing multiply-adds (or, for Exp and Log on
// some CPUs, by using assembly), which moves the core functions by at most
// one ulp. Functions built from several of those steps (Log10 = Log·(1/ln10),
// Erfc = Exp·Exp/x, Lgamma) accumulate a little more; Pow is checked
// separately (powClose). Where a function
// approaches a zero by cancellation (Lgamma near 1 and 2 and between its
// negative poles, Sin and Cos near multiples of pi/2, Log2 near 1) the
// result is tiny and the difference is bounded instead in ulps of 1 (cancel).
var bounds = map[string]struct {
	ulp    uint64
	cancel bool
}{
	"Exp": {2, false}, "Exp2": {1, false}, "Expm1": {1, false},
	"Log": {1, false}, "Log1p": {1, false}, "Log2": {1, true}, "Log10": {2, false},
	"Erf": {1, false}, "Erfc": {4, false},
	"Lgamma": {4, true}, "Sin": {1, true}, "Cos": {1, true},
}

// ulpOf returns the spacing of float64 values at |x|.
func ulpOf(x float64) float64 {
	x = math.Abs(x)
	return math.Nextafter(x, math.Inf(1)) - x
}

// closeness classifies got against want under the function's bound: 0 when
// within the ulp bound, 1 when within the absolute bound near a zero, and -1
// when outside both.
func closeness(name string, got, want float64) int {
	b := bounds[name]
	if ulps(got, want) <= b.ulp {
		return 0
	}
	finite := !math.IsNaN(got) && !math.IsNaN(want) && !math.IsInf(got, 0) && !math.IsInf(want, 0)
	if b.cancel && finite && math.Abs(got-want) <= float64(b.ulp)*ulpOf(math.Max(math.Abs(want), 1)) {
		return 1
	}
	return -1
}

// powClose reports whether Pow(x, y) = got is close enough to math.Pow's
// want. Pow computes x^y as Exp(yf·Log(x))·x^yi; a one-ulp difference in
// Log(x) is multiplied by |yf·ln x| ≤ |y·ln x| in the exponent, so the
// accepted difference grows with that magnitude (twice it, since an ulp is
// between 2^-53 and 2^-52 of the value).
func powClose(x, y, got, want float64) bool {
	d := ulps(got, want)
	return d == 0 || float64(d) <= 4+2*math.Abs(y*math.Log(math.Abs(x)))
}

// mathDefect reports whether x is an input on which the math package itself
// is wrong on this CPU, so that comparing with it says nothing. Both are in
// the amd64 assembly: math.Exp returns +Inf for x above about 709.44, where
// the result is finite but the scaled exponent reaches 1024 (exp_amd64.s
// treats that as overflow), and math.Log ignores the mantissa of subnormal
// inputs (log_amd64.s does not normalise them), which Log10, Pow and Lgamma
// inherit.
func mathDefect(name string, x float64) bool {
	if runtime.GOARCH != "amd64" {
		return false
	}
	subnormal := x != 0 && math.Abs(x) < 0x1p-1022
	switch name {
	case "Exp":
		return x > 709.43 && x <= 709.782712893383973096
	case "Log", "Log10", "Lgamma", "Pow":
		return subnormal
	}
	return false
}

// identityRef returns F(x) computed from F at a scaled argument, for the
// inputs where math is wrong: e^x = (e^(x/2))² and
// ln x = ln(x·2^54) − 54·ln 2.
func identityRef(name string, x float64) (float64, bool) {
	switch name {
	case "Exp":
		h := Exp(x / 2)
		return h * h, true
	case "Log":
		return Log(x*0x1p54) - float64(54*math.Ln2), true
	case "Log10":
		return float64(Log(x*0x1p54)-float64(54*math.Ln2)) * (1 / math.Ln10), true
	}
	return 0, false
}

// isSpecialInput reports whether x is one of the inputs the functions treat
// as special cases (NaN, ±Inf, ±0), whose results must match exactly.
func isSpecialInput(x float64) bool {
	return math.IsNaN(x) || math.IsInf(x, 0) || x == 0
}

// pureGoAMD64 reports whether this is amd64 at the default GOAMD64 level,
// where the compiler does not fuse and the math package's pure-Go functions
// are therefore the same unfused evaluation as this package.
func pureGoAMD64() bool {
	if runtime.GOARCH != "amd64" {
		return false
	}
	if bi, ok := debug.ReadBuildInfo(); ok {
		for _, s := range bi.Settings {
			if s.Key == "GOAMD64" && s.Value != "v1" && s.Value != "v2" {
				return false
			}
		}
	}
	return true
}

// pureOnAMD64 are the functions whose math package implementation on amd64
// is the pure-Go algorithm with no call to the assembly Exp or Log.
var pureOnAMD64 = map[string]bool{"Expm1": true, "Log1p": true, "Sin": true, "Cos": true}

// TestUnaryAgreesWithMath compares each function with the math package on
// edge cases and on random inputs: special inputs must give exactly the same
// result, and every result must be within the function's bound. On amd64 at
// GOAMD64=v1, where the compiler does not fuse, the functions whose math
// implementation is pure Go must match bit for bit.
func TestUnaryAgreesWithMath(t *testing.T) {
	n := 1_000_000
	if testing.Short() {
		n = 50_000
	}
	exactRef := pureGoAMD64()
	for _, f := range unary {
		t.Run(f.name, func(t *testing.T) {
			r := rand.New(rand.NewPCG(1, uint64(len(f.name))))
			diff, nearZero, bad, defects := 0, 0, 0, 0
			var maxd uint64
			check := func(x float64) {
				got, want := f.got(x), f.std(x)
				if sameBits(got, want) {
					return
				}
				diff++
				if mathDefect(f.name, x) {
					defects++
					if ref, ok := identityRef(f.name, x); ok && ulps(got, ref) > 4 {
						t.Errorf("%s(%v) = %v; math is wrong here on %s (%v) and the identity gives %v",
							f.name, x, got, runtime.GOARCH, want, ref)
					}
					return
				}
				if isSpecialInput(x) {
					t.Errorf("%s(%v) = %v, math gives %v (special case)", f.name, x, got, want)
				}
				if exactRef && pureOnAMD64[f.name] {
					t.Errorf("%s(%v) = %#x, pure-Go math gives %#x on amd64", f.name, x,
						math.Float64bits(got), math.Float64bits(want))
				}
				switch closeness(f.name, got, want) {
				case 0:
					maxd = max(maxd, ulps(got, want))
				case 1:
					nearZero++
				default:
					bad++
					if bad <= 20 {
						t.Errorf("%s(%v) = %v (%#x), math gives %v (%#x): %d ulp", f.name, x, got,
							math.Float64bits(got), want, math.Float64bits(want), ulps(got, want))
					}
				}
			}
			for _, x := range specials {
				check(x)
			}
			for i := 0; i < n; i++ {
				if i%4 == 0 {
					check(anyFinite(r))
				} else {
					check(f.gen(r))
				}
			}
			t.Logf("%s: %d of %d inputs differ from math.%s on %s: max %d ulp, %d more near a zero within %d ulp of 1, %d where math is wrong",
				f.name, diff, n+len(specials), f.name, runtime.GOARCH, maxd, nearZero, bounds[f.name].ulp, defects)
		})
	}
}

// TestPowAgreesWithMath is the two-argument version of
// TestUnaryAgreesWithMath for Pow.
func TestPowAgreesWithMath(t *testing.T) {
	n := 1_000_000
	if testing.Short() {
		n = 50_000
	}
	vals := append([]float64{0.25, 4, 1e-10, 1e10, 7, -7, 0.1}, specials...)
	for _, x := range vals {
		for _, y := range vals {
			if mathDefect("Pow", x) {
				continue
			}
			got, want := Pow(x, y), math.Pow(x, y)
			if math.IsNaN(want) || math.IsInf(want, 0) || want == 0 {
				if !sameBits(got, want) {
					t.Errorf("Pow(%v, %v) = %v, math gives %v (special case)", x, y, got, want)
				}
			} else if !powClose(x, y, got, want) {
				t.Errorf("Pow(%v, %v) = %v, math gives %v: %d ulp", x, y, got, want, ulps(got, want))
			}
		}
	}
	r := rand.New(rand.NewPCG(2, 3))
	diff, maxd := 0, uint64(0)
	for i := 0; i < n; i++ {
		var x, y float64
		switch i % 4 {
		case 0:
			x, y = anyFinite(r), anyFinite(r)
		case 1:
			x, y = -logUniform(r, -3, 3), float64(r.IntN(200)-100) // negative base, integer power
		default:
			x, y = logUniform(r, -30, 30), uniform(r, -10, 10)
		}
		if mathDefect("Pow", x) {
			continue
		}
		got, want := Pow(x, y), math.Pow(x, y)
		if !sameBits(got, want) {
			diff++
			maxd = max(maxd, ulps(got, want))
		}
		if !powClose(x, y, got, want) {
			t.Fatalf("Pow(%v, %v) = %v (%#x), math gives %v (%#x): %d ulp", x, y, got,
				math.Float64bits(got), want, math.Float64bits(want), ulps(got, want))
		}
	}
	t.Logf("Pow: %d of %d random inputs differ from math.Pow on %s: max %d ulp", diff, n, runtime.GOARCH, maxd)
	// Integer powers of two and of small integers are exact, as the output
	// path relies on (4^k, 10^-q/10 at q = 0).
	for k := 0; k <= 26; k++ {
		if got := Pow(4, float64(k)); got != float64(int64(1)<<(2*k)) {
			t.Errorf("Pow(4, %d) = %v", k, got)
		}
	}
}

// mul, add and sub are kept out of line so that the reference evaluation
// below rounds after every operation on every architecture, independently
// of the explicit conversions in the code under test.
//
//go:noinline
func mul(a, b float64) float64 { return a * b }

//go:noinline
func add(a, b float64) float64 { return a + b }

//go:noinline
func sub(a, b float64) float64 { return a - b }

// refExp is the math package's pure-Go exp (src/math/exp.go) with every
// operation an out-of-line call.
func refExp(x float64) float64 {
	const (
		Ln2Hi     = 6.93147180369123816490e-01
		Ln2Lo     = 1.90821492927058770002e-10
		Log2e     = 1.44269504088896338700e+00
		Overflow  = 7.09782712893383973096e+02
		Underflow = -7.45133219101941108420e+02
		NearZero  = 1.0 / (1 << 28)
		P1        = 1.66666666666666657415e-01
		P2        = -2.77777777770155933842e-03
		P3        = 6.61375632143793436117e-05
		P4        = -1.65339022054652515390e-06
		P5        = 4.13813679705723846039e-08
	)
	switch {
	case math.IsNaN(x):
		return x
	case x > Overflow:
		return math.Inf(1)
	case x < Underflow:
		return 0
	case -NearZero < x && x < NearZero:
		return add(1, x)
	}
	var k int
	switch {
	case x < 0:
		k = int(sub(mul(Log2e, x), 0.5))
	case x > 0:
		k = int(add(mul(Log2e, x), 0.5))
	}
	hi := sub(x, mul(float64(k), Ln2Hi))
	lo := mul(float64(k), Ln2Lo)
	r := sub(hi, lo)
	t := mul(r, r)
	c := sub(r, mul(t, add(P1, mul(t, add(P2, mul(t, add(P3, mul(t, add(P4, mul(t, P5))))))))))
	y := sub(1, sub(sub(lo, mul(r, c)/sub(2, c)), hi))
	return math.Ldexp(y, k)
}

// refLog is the math package's pure-Go log (src/math/log.go) with every
// operation an out-of-line call.
func refLog(x float64) float64 {
	const (
		Ln2Hi = 6.93147180369123816490e-01
		Ln2Lo = 1.90821492927058770002e-10
		L1    = 6.666666666666735130e-01
		L2    = 3.999999999940941908e-01
		L3    = 2.857142874366239149e-01
		L4    = 2.222219843214978396e-01
		L5    = 1.818357216161805012e-01
		L6    = 1.531383769920937332e-01
		L7    = 1.479819860511658591e-01
	)
	switch {
	case math.IsNaN(x) || math.IsInf(x, 1):
		return x
	case x < 0:
		return math.NaN()
	case x == 0:
		return math.Inf(-1)
	}
	f1, ki := math.Frexp(x)
	if f1 < math.Sqrt2/2 {
		f1 = mul(f1, 2)
		ki--
	}
	f := sub(f1, 1)
	k := float64(ki)
	s := f / add(2, f)
	s2 := mul(s, s)
	s4 := mul(s2, s2)
	t1 := mul(s2, add(L1, mul(s4, add(L3, mul(s4, add(L5, mul(s4, L7)))))))
	t2 := mul(s4, add(L2, mul(s4, add(L4, mul(s4, L6)))))
	R := add(t1, t2)
	hfsq := mul(mul(0.5, f), f)
	return sub(mul(k, Ln2Hi), sub(sub(hfsq, add(mul(s, add(hfsq, R)), mul(k, Ln2Lo))), f))
}

// TestExpLogAreUnfusedEvaluation checks Exp and Log bit for bit against an
// independent step-by-step IEEE evaluation of the same algorithm, on every
// architecture. It also reports how often the math package differs, which is
// the CPU dependence this package removes.
func TestExpLogAreUnfusedEvaluation(t *testing.T) {
	n := 2_000_000
	if testing.Short() {
		n = 100_000
	}
	r := rand.New(rand.NewPCG(4, 5))
	for _, tc := range []struct {
		name     string
		got, ref func(float64) float64
		std      func(float64) float64
		gen      func() float64
	}{
		{"Exp", Exp, refExp, math.Exp, func() float64 { return uniform(r, -746, 710) }},
		{"Log", Log, refLog, math.Log, func() float64 { return logUniform(r, -320, 308) }},
	} {
		stdDiff := 0
		for i := 0; i < n+len(specials); i++ {
			var x float64
			switch {
			case i < len(specials):
				x = specials[i]
			case i%8 == 0:
				x = anyFinite(r)
			default:
				x = tc.gen()
			}
			got, ref := tc.got(x), tc.ref(x)
			if !sameBits(got, ref) {
				t.Fatalf("%s(%v) = %#x, unfused reference gives %#x", tc.name, x,
					math.Float64bits(got), math.Float64bits(ref))
			}
			if !sameBits(got, tc.std(x)) {
				stdDiff++
			}
		}
		t.Logf("%s: identical to the unfused reference on %d inputs; math.%s differs on %d (%.3f%%) on %s",
			tc.name, n+len(specials), tc.name, stdDiff, float64(100*stdDiff)/float64(n), runtime.GOARCH)
	}
}

// TestVariatesMatchRand checks that NormFloat64 and ExpFloat64 return what
// rand.(*Rand).NormFloat64 and ExpFloat64 return from the same state and
// leave the generator in the same state, except on the rare draws where the
// math package's Exp or Log (or its fused float32 arithmetic) decides a
// rejection test or a tail value differently.
func TestVariatesMatchRand(t *testing.T) {
	n := 2_000_000
	if testing.Short() {
		n = 100_000
	}
	for _, tc := range []struct {
		name string
		got  func(*rand.Rand) float64
		std  func(*rand.Rand) float64
	}{
		{"NormFloat64", NormFloat64, (*rand.Rand).NormFloat64},
		{"ExpFloat64", ExpFloat64, (*rand.Rand).ExpFloat64},
	} {
		t.Run(tc.name, func(t *testing.T) {
			// One long stream, compared draw by draw. After each pair of
			// draws the next Uint64 of both generators must agree, which
			// shows the same number of values was consumed. A mismatch can
			// change that number, so both are then reseeded together.
			pa, pb := rand.NewPCG(7, 0), rand.NewPCG(7, 0)
			ra, rb := rand.New(pa), rand.New(pb)
			valueDiff, stateDiff := 0, 0
			for i := 0; i < n; i++ {
				a, b := tc.got(ra), tc.std(rb)
				same := sameBits(a, b)
				if !same {
					valueDiff++
				}
				if ra.Uint64() != rb.Uint64() {
					stateDiff++
					if same {
						t.Fatalf("draw %d: equal values %v but different generator state", i, a)
					}
				}
				if !same {
					pa.Seed(7, uint64(i+1))
					pb.Seed(7, uint64(i+1))
				}
			}
			frac := float64(valueDiff) / float64(n)
			t.Logf("%s: %d of %d draws (%.5f%%) differ from rand.%s on %s (%d of them consumed a different number of values)",
				tc.name, valueDiff, n, 100*frac, tc.name, runtime.GOARCH, stateDiff)
			if frac > 1e-4 {
				t.Errorf("%s differs from rand on %.5f%% of draws, want at most 0.01%%", tc.name, 100*frac)
			}
		})
	}
}

// TestVariateDistributions checks the first two moments and the tail mass
// beyond the ziggurat's base strip, which is drawn by the code that calls
// Exp and Log.
func TestVariateDistributions(t *testing.T) {
	const n = 2_000_000
	r := rand.New(rand.NewPCG(11, 12))
	check := func(name string, got, want, sd float64) {
		t.Helper()
		if math.Abs(got-want) > 4*sd {
			t.Errorf("%s = %v, want %v ± %v (4σ)", name, got, want, 4*sd)
		}
	}
	var s, s2 float64
	tail := 0
	for i := 0; i < n; i++ {
		x := NormFloat64(r)
		s += x
		s2 += float64(x * x)
		if math.Abs(x) > rn {
			tail++
		}
	}
	mean := s / n
	check("normal mean", mean, 0, 1/math.Sqrt(n))
	check("normal variance", s2/n-float64(mean*mean), 1, math.Sqrt(2.0/n)) // (μ4 − σ⁴)/n with μ4 = 3
	pTail := math.Erfc(rn / math.Sqrt2)
	check("normal tail beyond 3.4426", float64(tail), n*pTail, math.Sqrt(n*pTail))

	s, s2, tail = 0, 0, 0
	for i := 0; i < n; i++ {
		x := ExpFloat64(r)
		if x <= 0 {
			t.Fatalf("ExpFloat64 = %v", x)
		}
		s += x
		s2 += float64(x * x)
		if x > re {
			tail++
		}
	}
	mean = s / n
	check("exponential mean", mean, 1, 1/math.Sqrt(n))
	check("exponential variance", s2/n-float64(mean*mean), 1, math.Sqrt(8.0/n)) // (μ4 − σ⁴)/n with μ4 = 9
	pTail = math.Exp(-re)
	check("exponential tail beyond 7.697", float64(tail), n*pTail, math.Sqrt(n*pTail))
}

func BenchmarkExp(b *testing.B) {
	x, s := 0.0, 0.0
	for i := 0; i < b.N; i++ {
		s += Exp(x)
		x += 0x1p-10
		if x > 5 {
			x = -5
		}
	}
	sink = s
}

func BenchmarkMathExp(b *testing.B) {
	x, s := 0.0, 0.0
	for i := 0; i < b.N; i++ {
		s += math.Exp(x)
		x += 0x1p-10
		if x > 5 {
			x = -5
		}
	}
	sink = s
}

func BenchmarkLog(b *testing.B) {
	x, s := 1e-3, 0.0
	for i := 0; i < b.N; i++ {
		s += Log(x)
		x += 0x1p-10
		if x > 100 {
			x = 1e-3
		}
	}
	sink = s
}

func BenchmarkMathLog(b *testing.B) {
	x, s := 1e-3, 0.0
	for i := 0; i < b.N; i++ {
		s += math.Log(x)
		x += 0x1p-10
		if x > 100 {
			x = 1e-3
		}
	}
	sink = s
}

func BenchmarkNormFloat64(b *testing.B) {
	r := rand.New(rand.NewPCG(1, 2))
	s := 0.0
	for i := 0; i < b.N; i++ {
		s += NormFloat64(r)
	}
	sink = s
}

func BenchmarkRandNormFloat64(b *testing.B) {
	r := rand.New(rand.NewPCG(1, 2))
	s := 0.0
	for i := 0; i < b.N; i++ {
		s += r.NormFloat64()
	}
	sink = s
}

var sink float64
