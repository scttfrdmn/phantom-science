//go:build !release

package readlen

import (
	"errors"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

const nSamples = 100_000

// ksCritical is the one-sample Kolmogorov–Smirnov critical value at α = 0.01
// (conservative for discrete distributions).
func ksCritical(n int) float64 { return 1.63 / math.Sqrt(float64(n)) }

func mustNew(t *testing.T, spec types.ReadLengthSpec) *Distribution {
	t.Helper()
	d, err := New(spec, 0)
	if err != nil {
		t.Fatalf("New(%+v): %v", spec, err)
	}
	return d
}

func draw(d *Distribution, seed uint64, n int) []int {
	r := rand.New(rand.NewPCG(seed, 0xda7a))
	xs := make([]int, n)
	for i := range xs {
		xs[i] = d.Sample(r)
	}
	return xs
}

// ksDistance is the KS statistic of sorted integer samples against a
// discrete CDF, evaluated on both sides of every jump.
func ksDistance(xs []int, cdf func(int) float64) float64 {
	xs = slices.Clone(xs)
	slices.Sort(xs)
	n := float64(len(xs))
	var dmax float64
	for i := 0; i < len(xs); {
		j := i
		for j < len(xs) && xs[j] == xs[i] {
			j++
		}
		below, at := float64(i)/n, float64(j)/n
		dmax = math.Max(dmax, math.Abs(at-cdf(xs[i])))
		dmax = math.Max(dmax, math.Abs(below-cdf(xs[i]-1)))
		i = j
	}
	return dmax
}

func meanSD(xs []int) (mean, sd float64) {
	var s, s2 float64
	for _, x := range xs {
		s += float64(x)
	}
	mean = s / float64(len(xs))
	for _, x := range xs {
		s2 += (float64(x) - mean) * (float64(x) - mean)
	}
	return mean, math.Sqrt(s2 / float64(len(xs)-1))
}

// n50 is the length L such that reads of length >= L hold at least half of
// all bases.
func n50(xs []int) int {
	xs = slices.Clone(xs)
	slices.Sort(xs)
	var total, acc int64
	for _, x := range xs {
		total += int64(x)
	}
	for i := len(xs) - 1; i >= 0; i-- {
		acc += int64(xs[i])
		if 2*acc >= total {
			return xs[i]
		}
	}
	return 0
}

func median(xs []int) float64 {
	xs = slices.Clone(xs)
	slices.Sort(xs)
	return float64(xs[len(xs)/2])
}

// TestContinuous_KS checks sampled lengths against the analytic CDF of each
// parametrisation, and the sample mean against the analytic mean.
func TestContinuous_KS(t *testing.T) {
	cases := []struct {
		name     string
		spec     types.ReadLengthSpec
		wantMean float64 // analytic untruncated mean
	}{
		{"lognormal-n50-median", types.ReadLengthSpec{Dist: "lognormal", N50: 50000, Median: 25000}, 25000 * math.Sqrt2},
		{"lognormal-mean-sd", types.ReadLengthSpec{Dist: "lognormal", Mean: 1500, SD: 1000}, 1500},
		{"gamma-promethion", types.ReadLengthSpec{Dist: "gamma", Shape: 1.5, Scale: 10000}, 15000},
		// Shape 0.5 puts 1.3% of the mass below ½, which Min 1 truncates:
		// mean = kθ (1 − P(1.5, x)) / (1 − P(0.5, x)) with x = 0.5/θ,
		// P(0.5, x) = erf(√x), P(1.5, x) = erf(√x) − 2√(x/π) e^−x.
		{"gamma-shape-below-1", types.ReadLengthSpec{Dist: "gamma", Shape: 0.5, Scale: 4000}, func() float64 {
			x := 0.5 / 4000
			p05 := math.Erf(math.Sqrt(x))
			p15 := p05 - 2*math.Sqrt(x/math.Pi)*math.Exp(-x)
			return 2000 * (1 - p15) / (1 - p05)
		}()},
		{"gamma-mean-sd", types.ReadLengthSpec{Dist: "gamma", Mean: 1500, SD: 1200}, 1500},
		{"lognormal-truncated", types.ReadLengthSpec{Dist: "lognormal", N50: 50000, Median: 25000, Min: 5000, Max: 60000}, 0},
		{"gamma-truncated", types.ReadLengthSpec{Dist: "gamma", Shape: 2, Scale: 12500, Min: 1000, Max: 40000}, 0},
	}
	for i, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			d := mustNew(t, tc.spec)
			if tc.wantMean != 0 && math.Abs(d.Mean()-tc.wantMean) > 1e-3*tc.wantMean {
				t.Errorf("Mean() = %.2f, analytic %.2f", d.Mean(), tc.wantMean)
			}
			xs := draw(d, uint64(100+i), nSamples)
			lo, hi := d.Bounds()
			for _, x := range xs {
				if x < lo || (hi != 0 && x > hi) {
					t.Fatalf("sample %d outside [%d, %d]", x, lo, hi)
				}
			}
			ks := ksDistance(xs, func(k int) float64 { return d.CDF(k) })
			if crit := ksCritical(len(xs)); ks >= crit {
				t.Errorf("KS D = %.5f >= %.5f", ks, crit)
			}
			m, _ := meanSD(xs)
			if tol := 4 * d.SD() / math.Sqrt(float64(len(xs))); math.Abs(m-d.Mean()) > tol {
				t.Errorf("sample mean %.1f, want %.1f ± %.1f", m, d.Mean(), tol)
			}
			t.Logf("KS D = %.5f (crit %.5f), mean %.1f (want %.1f), sd %.1f (want %.1f)",
				ks, ksCritical(len(xs)), m, d.Mean(), func() float64 { _, s := meanSD(xs); return s }(), d.SD())
		})
	}
}

// TestKS_DetectsWrongDistribution guards the KS harness itself: samples from
// one distribution must be rejected against a different one.
func TestKS_DetectsWrongDistribution(t *testing.T) {
	a := mustNew(t, types.ReadLengthSpec{Dist: "lognormal", N50: 50000, Median: 25000})
	b := mustNew(t, types.ReadLengthSpec{Dist: "lognormal", N50: 52000, Median: 25000})
	xs := draw(a, 1, nSamples)
	if ks := ksDistance(xs, func(k int) float64 { return b.CDF(k) }); ks < ksCritical(len(xs)) {
		t.Errorf("KS D = %.5f did not reject a 4%% N50 error", ks)
	}
}

// TestLognormal_N50Median checks the sample N50 and median of 10⁵ reads
// against the analytic values, with a tolerance of 4σ of the spread of 50
// replicate samples.
func TestLognormal_N50Median(t *testing.T) {
	const wantN50, wantMedian, reps = 50000, 25000, 50
	d := mustNew(t, types.ReadLengthSpec{Dist: "lognormal", N50: wantN50, Median: wantMedian})
	n50s := make([]float64, reps)
	meds := make([]float64, reps)
	for i := range reps {
		xs := draw(d, uint64(1000+i), nSamples)
		n50s[i], meds[i] = float64(n50(xs)), median(xs)
	}
	sdOf := func(v []float64) float64 {
		var m, s float64
		for _, x := range v {
			m += x
		}
		m /= float64(len(v))
		for _, x := range v {
			s += (x - m) * (x - m)
		}
		return math.Sqrt(s / float64(len(v)-1))
	}
	sdN50, sdMed := sdOf(n50s), sdOf(meds)
	xs := draw(d, 7, nSamples)
	gotN50, gotMed := float64(n50(xs)), median(xs)
	t.Logf("N50 %.0f (want %d, replicate sd %.1f); median %.0f (want %d, replicate sd %.1f)",
		gotN50, wantN50, sdN50, gotMed, wantMedian, sdMed)
	if math.Abs(gotN50-wantN50) > 4*sdN50 {
		t.Errorf("N50 = %.0f, want %d ± %.0f", gotN50, wantN50, 4*sdN50)
	}
	if math.Abs(gotMed-wantMedian) > 4*sdMed {
		t.Errorf("median = %.0f, want %d ± %.0f", gotMed, wantMedian, 4*sdMed)
	}
	// Sanity: the replicate spread is a sampling spread, not a huge number
	// that would make the check vacuous (~0.5% of N50 at n = 10⁵).
	if sdN50 > 0.02*wantN50 || sdMed > 0.02*wantMedian {
		t.Errorf("replicate spread too large: sd(N50) %.0f, sd(median) %.0f", sdN50, sdMed)
	}
}

func TestEmpirical_MatchesHistogram(t *testing.T) {
	path := filepath.Join("testdata", "cdna_hist.tsv")
	h, err := LoadHistogram(path)
	if err != nil {
		t.Fatal(err)
	}
	for _, tc := range []struct {
		name     string
		min, max int
	}{{"full", 0, 0}, {"truncated", 500, 5000}} {
		t.Run(tc.name, func(t *testing.T) {
			d := mustNew(t, types.ReadLengthSpec{Dist: "empirical", EmpiricalPath: path, Min: tc.min, Max: tc.max})
			in := func(l int) bool { return l >= max(tc.min, 1) && (tc.max == 0 || l <= tc.max) }
			var total float64
			for l, c := range h {
				if in(l) {
					total += c
				}
			}
			xs := draw(d, 11, nSamples)
			obs := map[int]int{}
			for _, x := range xs {
				if _, ok := h[x]; !ok || !in(x) {
					t.Fatalf("sampled length %d is not an allowed histogram bin", x)
				}
				obs[x]++
			}
			n := float64(len(xs))
			chi2, bins := 0.0, 0
			for _, l := range h.Lengths() {
				p := 0.0
				if in(l) {
					p = h[l] / total
				}
				if p == 0 {
					if obs[l] != 0 {
						t.Errorf("length %d has probability 0 but was sampled %d times", l, obs[l])
					}
					continue
				}
				exp := n * p
				if tol := 4 * math.Sqrt(n*p*(1-p)); math.Abs(float64(obs[l])-exp) > tol {
					t.Errorf("length %d: %d samples, want %.0f ± %.0f", l, obs[l], exp, tol)
				}
				chi2 += (float64(obs[l]) - exp) * (float64(obs[l]) - exp) / exp
				bins++
			}
			// χ² with k = bins−1 degrees of freedom: mean k, sd √(2k).
			k := float64(bins - 1)
			if chi2 > k+4*math.Sqrt(2*k) {
				t.Errorf("χ² = %.1f with %v dof", chi2, k)
			}
			if ks, crit := ksDistance(xs, d.CDF), ksCritical(len(xs)); ks >= crit {
				t.Errorf("KS D = %.5f >= %.5f", ks, crit)
			}
			m, _ := meanSD(xs)
			if tol := 4 * d.SD() / math.Sqrt(n); math.Abs(m-d.Mean()) > tol {
				t.Errorf("sample mean %.1f, want %.1f ± %.1f", m, d.Mean(), tol)
			}
			t.Logf("χ² = %.1f (%v dof), mean %.1f (want %.1f)", chi2, k, m, d.Mean())
		})
	}
}

func TestLoadHistogram_Formats(t *testing.T) {
	dir := t.TempDir()
	write := func(name, body string) string {
		p := filepath.Join(dir, name)
		if err := os.WriteFile(p, []byte(body), 0o600); err != nil {
			t.Fatal(err)
		}
		return p
	}
	h, err := LoadHistogram(write("ok.csv", "length,count\n# comment\n\n100,2\n200, 3.5\n100 1\n"))
	if err != nil {
		t.Fatal(err)
	}
	if len(h) != 2 || h[100] != 3 || h[200] != 3.5 {
		t.Errorf("histogram = %v, want map[100:3 200:3.5]", h)
	}
	for name, body := range map[string]string{
		"empty":        "# nothing\n",
		"zero-total":   "100\t0\n",
		"neg-count":    "100\t-1\n",
		"zero-length":  "0\t5\n",
		"float-length": "1.5\t5\n",
		"three-fields": "100\t5\t7\n",
		"late-header":  "100\t5\nlength\tcount\n",
	} {
		if _, err := LoadHistogram(write(name, body)); !errors.Is(err, ErrInvalidSpec) {
			t.Errorf("%s: err = %v, want ErrInvalidSpec", name, err)
		}
	}
}

func TestNew_Invalid(t *testing.T) {
	missing := filepath.Join(t.TempDir(), "missing.tsv")
	cases := map[string]types.ReadLengthSpec{
		"n50-equals-median":   {Dist: "lognormal", N50: 25000, Median: 25000},
		"n50-below-median":    {Dist: "lognormal", N50: 20000, Median: 25000},
		"n50-without-median":  {Dist: "lognormal", N50: 50000},
		"lognormal-sd-zero":   {Dist: "lognormal", Mean: 1000},
		"lognormal-sd-neg":    {Dist: "lognormal", Mean: 1000, SD: -1},
		"lognormal-mean-neg":  {Dist: "lognormal", Mean: -1000, SD: 10},
		"lognormal-both":      {Dist: "lognormal", N50: 50000, Median: 25000, Mean: 1000, SD: 10},
		"lognormal-no-params": {Dist: "lognormal"},
		"lognormal-shape":     {Dist: "lognormal", Shape: 2, Scale: 10},
		"gamma-shape-zero":    {Dist: "gamma", Scale: 1000},
		"gamma-shape-neg":     {Dist: "gamma", Shape: -1, Scale: 1000},
		"gamma-scale-zero":    {Dist: "gamma", Shape: 2},
		"gamma-sd-zero":       {Dist: "gamma", Mean: 1000},
		"gamma-both":          {Dist: "gamma", Shape: 2, Scale: 10, Mean: 20, SD: 14},
		"gamma-n50":           {Dist: "gamma", N50: 50000, Median: 25000},
		"gamma-nan":           {Dist: "gamma", Shape: math.NaN(), Scale: 10},
		"empirical-no-path":   {Dist: "empirical"},
		"empirical-missing":   {Dist: "empirical", EmpiricalPath: missing},
		"empirical-mean":      {Dist: "empirical", EmpiricalPath: filepath.Join("testdata", "cdna_hist.tsv"), Mean: 3},
		"empirical-no-mass":   {Dist: "empirical", EmpiricalPath: filepath.Join("testdata", "cdna_hist.tsv"), Min: 21000},
		"min-neg":             {Dist: "gamma", Shape: 2, Scale: 10, Min: -1},
		"max-neg":             {Dist: "gamma", Shape: 2, Scale: 10, Max: -1},
		"max-below-min":       {Dist: "gamma", Shape: 2, Scale: 10, Min: 100, Max: 50},
		"bounds-exclude-mass": {Dist: "lognormal", N50: 50000, Median: 25000, Min: 2_000_000},
		"fixed-with-params":   {Dist: "fixed", Mean: 1000},
		"default-with-params": {N50: 50000, Median: 25000},
	}
	for name, spec := range cases {
		if _, err := New(spec, 150); !errors.Is(err, ErrInvalidSpec) {
			t.Errorf("%s: err = %v, want ErrInvalidSpec", name, err)
		}
	}
	if _, err := New(types.ReadLengthSpec{Dist: "fixed"}, 0); !errors.Is(err, ErrInvalidSpec) {
		t.Errorf("fixed length 0: err = %v, want ErrInvalidSpec", err)
	}
	if _, err := New(types.ReadLengthSpec{Dist: "weibull"}, 150); err == nil || !strings.Contains(err.Error(), "unknown") {
		t.Errorf("unknown dist: err = %v", err)
	}
	if _, err := New(types.ReadLengthSpec{Dist: "empirical", EmpiricalPath: missing}, 150); err == nil || !strings.Contains(err.Error(), "missing.tsv") {
		t.Errorf("missing file error does not name the file: %v", err)
	}
}

func TestFixed_NoRandomness(t *testing.T) {
	for _, dist := range []string{"", "fixed"} {
		d, err := New(types.ReadLengthSpec{Dist: dist}, 150)
		if err != nil {
			t.Fatal(err)
		}
		r := rand.New(rand.NewPCG(1, 2))
		ref := rand.New(rand.NewPCG(1, 2))
		for range 100 {
			if n := d.Sample(r); n != 150 {
				t.Fatalf("fixed sample %d, want 150", n)
			}
		}
		if r.Uint64() != ref.Uint64() {
			t.Error("fixed sampler consumed randomness")
		}
		if !d.Fixed() || d.Mean() != 150 || d.SD() != 0 {
			t.Errorf("fixed: Fixed %v Mean %v SD %v", d.Fixed(), d.Mean(), d.SD())
		}
	}
}

func TestSample_Reproducible(t *testing.T) {
	for _, spec := range []types.ReadLengthSpec{
		{Dist: "lognormal", N50: 50000, Median: 25000},
		{Dist: "gamma", Shape: 0.7, Scale: 3000},
		{Dist: "empirical", EmpiricalPath: filepath.Join("testdata", "cdna_hist.tsv")},
	} {
		d := mustNew(t, spec)
		if a, b := draw(d, 5, 1000), draw(d, 5, 1000); !slices.Equal(a, b) {
			t.Errorf("%s: same seed gave different lengths", spec.Dist)
		}
		if a, b := draw(d, 5, 1000), draw(d, 6, 1000); slices.Equal(a, b) {
			t.Errorf("%s: different seeds gave identical lengths", spec.Dist)
		}
	}
}

// TestTruncation_TightBounds uses bounds that reject most draws: every
// sample must still be in bounds and the conditional distribution correct.
func TestTruncation_TightBounds(t *testing.T) {
	d := mustNew(t, types.ReadLengthSpec{Dist: "gamma", Shape: 1.5, Scale: 10000, Min: 40000, Max: 45000})
	if d.accept > 0.05 {
		t.Fatalf("test needs a low acceptance rate, got %.3f", d.accept)
	}
	xs := draw(d, 3, 20000)
	for _, x := range xs {
		if x < 40000 || x > 45000 {
			t.Fatalf("sample %d outside [40000, 45000]", x)
		}
	}
	if ks := ksDistance(xs, func(k int) float64 { return d.CDF(k) }); ks >= ksCritical(len(xs)) {
		t.Errorf("KS D = %.5f >= %.5f", ks, ksCritical(len(xs)))
	}
}

// TestGammaP checks the incomplete gamma function the gamma KS test and
// moments depend on against closed forms.
func TestGammaP(t *testing.T) {
	for _, x := range []float64{0.01, 0.3, 1, 2.5, 7, 30, 200} {
		if got, want := gammaP(1, x), -math.Expm1(-x); math.Abs(got-want) > 1e-12 {
			t.Errorf("P(1, %g) = %.15g, want %.15g", x, got, want)
		}
		if got, want := gammaP(0.5, x), math.Erf(math.Sqrt(x)); math.Abs(got-want) > 1e-12 {
			t.Errorf("P(0.5, %g) = %.15g, want %.15g", x, got, want)
		}
		// P(2, x) = 1 − e^−x (1 + x)
		if got, want := gammaP(2, x), 1-math.Exp(-x)*(1+x); math.Abs(got-want) > 1e-12 {
			t.Errorf("P(2, %g) = %.15g, want %.15g", x, got, want)
		}
	}
}
