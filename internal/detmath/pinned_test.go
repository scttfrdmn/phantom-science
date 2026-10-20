//go:build !release

package detmath

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"hash"
	"math"
	"math/rand/v2"
	"runtime"
	"sort"
	"testing"
)

// portable returns n inputs built from random bits only (no floating-point
// arithmetic, so they are the same on every CPU): sign chosen by signs
// (+1, -1 or 0 for both), binary exponent uniform in [minExp, maxExp],
// random mantissa; every eighth input is any finite float64.
func portable(r *rand.Rand, n, minExp, maxExp, signs int) []float64 {
	out := make([]float64, n)
	for i := range out {
		if i%8 == 7 {
			out[i] = anyFinite(r)
			continue
		}
		e := minExp + r.IntN(maxExp-minExp+1)
		bits := uint64(e+1023)<<52 | r.Uint64()>>12
		if signs < 0 || signs == 0 && r.IntN(2) == 0 {
			bits |= 1 << 63
		}
		out[i] = math.Float64frombits(bits)
	}
	return out
}

// hashBits feeds the bits of v to h, with every NaN written as one value:
// the payload of a NaN produced by arithmetic differs between CPUs (amd64
// sets the sign bit, arm64 does not), and no output depends on it.
func hashBits(h hash.Hash, v float64) {
	if math.IsNaN(v) {
		v = math.NaN()
	}
	var b [8]byte
	binary.LittleEndian.PutUint64(b[:], math.Float64bits(v))
	h.Write(b[:])
}

// pinnedBits are SHA-256 digests of every function's results on fixed
// inputs, and of 1e6 draws of each variate. CI computes them on amd64 and
// arm64; they must never differ. They were recorded on darwin/arm64 and
// change only if an algorithm in this package changes.
var pinnedBits = map[string]string{
	"Cos":         "4a7730079b07ba7dae20b999e4829d337d4f4c4d5ed0d7432b2b9ee2ce650e4c",
	"Erf":         "eb203ab7d30fcf30ee986dd9d41bce4457013e89929bb1119a3324bb4977f10f",
	"Erfc":        "d55d5d97cff5df120fb3d07ab49e1efc10480658f11666d4917e19e3bba89e9e",
	"Exp":         "d8d49507226b452643a00e164c26b6f1fd7b38def441b6681f57d6c35a47516e",
	"Exp2":        "3bca5c4de4c15d2d6b46cb10f4f8197ff2675aa586208863f681b5c1ad467871",
	"ExpFloat64":  "39c72a34e464ae71885dd985fa17556607791343a223272b1bc19a64e6a1804d",
	"Expm1":       "70e74c57dea2f334d0a2c5e2827991dce77ccdb4d3a4f3605780380050379e1d",
	"Lgamma":      "68fa04b42ba2bd524dc7fb226705b65915189969a5c021d30a604ccc1602b687",
	"Log":         "b1e94d0b951871e7871252c5622d948c350375a4f9f38e707544d81299a5d6b9",
	"Log10":       "c2ec39bbacc325db32df4b1badfe5e50e6af323eab20c773b6860bdecce8498e",
	"Log1p":       "34018959837bf31776a2dcfac539796453dfbad42cac3724e91047671f2d0275",
	"Log2":        "7b098bff94ba291d5e24698f5f5e41ca48447f0d6e264d3f03031fe47be65989",
	"NormFloat64": "9ac967ac5292ae4f07d5d83de44a2b72f5b9756301421650c3d09524517cf372",
	"Pow":         "a43e6565690d38beab36e853f06b366257bfa051a443648cd1d0f12b64791a4e",
	"Sin":         "f5d3322136c75f9c50f586b6566c92bb1dce066b5f41b3adfb4a3bfdeb370e82",
}

// TestPinnedBits checks that each function returns the same bits on every
// CPU: the digests over 200,000 inputs per function are pinned.
func TestPinnedBits(t *testing.T) {
	r := rand.New(rand.NewPCG(2026, 33))
	const n = 200_000
	cases := []struct {
		name string
		f    func(float64) float64
		in   []float64
	}{
		{"Exp", Exp, portable(r, n, -30, 10, 0)},
		{"Exp2", Exp2, portable(r, n, -30, 11, 0)},
		{"Expm1", Expm1, portable(r, n, -60, 10, 0)},
		{"Log", Log, portable(r, n, -1022, 1023, 1)},
		{"Log2", Log2, portable(r, n, -1022, 1023, 1)},
		{"Log10", Log10, portable(r, n, -1022, 1023, 1)},
		{"Log1p", Log1p, append(portable(r, n/2, -60, 1023, 1), portable(r, n/2, -60, -1, -1)...)},
		{"Erf", Erf, portable(r, n, -30, 5, 0)},
		{"Erfc", Erfc, portable(r, n, -30, 5, 0)},
		{"Lgamma", lgammaV(Lgamma), portable(r, n, -30, 20, 0)},
		{"Sin", Sin, portable(r, n, -30, 40, 0)},
		{"Cos", Cos, portable(r, n, -30, 40, 0)},
	}
	got := map[string]string{}
	for _, c := range cases {
		h := sha256.New()
		for _, x := range append(c.in, specials...) {
			hashBits(h, c.f(x))
		}
		got[c.name] = hex.EncodeToString(h.Sum(nil))
	}
	h := sha256.New()
	xs, ys := portable(r, n, -60, 60, 1), portable(r, n, -10, 5, 0)
	for i := range xs {
		hashBits(h, Pow(xs[i], ys[i]))
	}
	for _, x := range specials {
		for _, y := range specials {
			hashBits(h, Pow(x, y))
		}
	}
	got["Pow"] = hex.EncodeToString(h.Sum(nil))
	for _, v := range []struct {
		name string
		f    func(*rand.Rand) float64
	}{{"NormFloat64", NormFloat64}, {"ExpFloat64", ExpFloat64}} {
		rr := rand.New(rand.NewPCG(1, 2))
		h := sha256.New()
		for i := 0; i < 1_000_000; i++ {
			hashBits(h, v.f(rr))
		}
		// The generator's state after the draws.
		var b [8]byte
		binary.LittleEndian.PutUint64(b[:], rr.Uint64())
		h.Write(b[:])
		got[v.name] = hex.EncodeToString(h.Sum(nil))
	}
	names := make([]string, 0, len(got))
	for k := range got {
		names = append(names, k)
	}
	sort.Strings(names)
	for _, k := range names {
		if want, ok := pinnedBits[k]; !ok {
			t.Errorf("%s: no pinned digest; got %q", k, got[k])
		} else if got[k] != want {
			t.Errorf("%s: digest %s on %s, want %s", k, got[k], runtime.GOARCH, want)
		}
	}
}
