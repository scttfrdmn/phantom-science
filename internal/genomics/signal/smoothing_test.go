//go:build !release

package signal

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"math"
	"strings"
	"testing"
)

// rawHash is the SHA-256 of a read's raw samples (little-endian int16) and
// its calibration, so a pinned hash detects any change in synthesis.
func rawHash(rd *Read) string {
	h := sha256.New()
	b := make([]byte, 2)
	for _, v := range rd.Raw {
		binary.LittleEndian.PutUint16(b, uint16(v))
		h.Write(b)
	}
	f := make([]byte, 8)
	for _, x := range []float64{rd.Offset, rd.MedianBefore} {
		binary.LittleEndian.PutUint64(f, math.Float64bits(x))
		h.Write(f)
	}
	return hex.EncodeToString(h.Sum(nil))
}

// TestLevelTauZeroIsStepFunction pins the output of both presets, with
// every feature on (adapter, poly tail, drift, stalls, skips, read jitter,
// white noise), to the hashes of the step-function synthesiser from before
// level_tau existed: level_tau = 0 must stay byte-identical. The pins are
// the unfused (no multiply-add) values, the same on every architecture; the
// RNA004 pin was first taken on arm64 with fused arithmetic and re-pinned to
// the amd64 value when fusion was removed (internal/fmacheck).
func TestLevelTauZeroIsStepFunction(t *testing.T) {
	want := map[string]string{
		PresetDNAR10: "ba6bccdcf7cbfe9406908b422ba5debcd4fd2bb488799c99d40fccc6d0ba56b4",
		PresetRNA004: "5b2e9da2ab68b441b7a3e1e847cf5d0fff38b186404565e3f709039d3db0dc94",
	}
	for _, name := range PresetNames() {
		p, _ := Preset(name)
		// Every modelling parameter is fixed here so that recalibrating a
		// preset does not move the pin; only the device geometry and
		// direction come from the preset.
		p.ModelK = 5
		p.SpeedBPS = p.SamplingRate / 11
		p.DwellCV, p.LevelNoiseScale, p.WhiteNoisePA = 0.4, 1.2, 0.7
		p.DriftPA, p.DriftTau, p.ReadShiftSD, p.ReadScaleSD = 0.8, 1.5, 2, 0.02
		p.StallProb, p.StallMin, p.StallMax, p.SkipProb = 0.005, 4, 12, 0.003
		p.AdapterSeconds, p.AdapterLevelPA, p.AdapterNoisePA = 0.1, 105, 6
		p.AdapterSeq, p.PolyANoisePA = kit14Adapter, 1.1
		p.LevelTau = 0
		seq := randomSeq(testRand(3), 3000)
		rd, err := SynthesizeRead(42, "7e1c3a52-0f7a-4c1e-9d55-3b0b8d9a2c11", seq, p, testModel(t), WithPolyA(30))
		if err != nil {
			t.Fatal(err)
		}
		if got := rawHash(rd); got != want[name] {
			t.Errorf("%s: level_tau=0 hash %s, want %s", name, got, want[name])
		}
	}
}

// TestLevelTauExponentialApproach checks, with noise off, that each sample
// of a k-mer segment sits at target + (prev − target)·exp(−(i+1)/τ), where
// prev is the level at the end of the previous segment, and that the
// events still record the target level.
func TestLevelTauExponentialApproach(t *testing.T) {
	for _, tau := range []float64{0.5, 2, 6} {
		p := cleanProfile(t)
		p.LevelNoiseScale, p.WhiteNoisePA = 0, 0
		p.DwellCV = 0.3
		if err := p.Set("level_tau", tau); err != nil {
			t.Fatal(err)
		}
		m := testModel(t)
		seq := randomSeq(testRand(9), 400)
		rd, err := Synthesize(testRand(5), seq, p, m, WithEvents())
		if err != nil {
			t.Fatal(err)
		}
		pa := rd.PA()
		tol := p.Range/(2*p.Digitisation) + 1e-9
		ev := kmerEvents(rd)
		prev := ev[0].LevelPA
		var steps, checked int
		for _, e := range ev {
			if want := float64(m.Mean[e.Kmer]); math.Abs(e.LevelPA-want) > 1e-4 {
				t.Fatalf("tau %g: event level %g, want target %g", tau, e.LevelPA, want)
			}
			if math.Abs(prev-e.LevelPA) > 5 {
				steps++
			}
			for i := 0; i < e.Len; i++ {
				want := e.LevelPA + (prev-e.LevelPA)*math.Exp(-float64(i+1)/tau)
				if got := pa[e.Start+i]; math.Abs(got-want) > tol {
					t.Fatalf("tau %g: sample %d of k-mer at %d: %g pA, want %g", tau, i, e.Pos, got, want)
				}
				checked++
			}
			if e.Len > 0 {
				prev = e.LevelPA + (prev-e.LevelPA)*math.Exp(-float64(e.Len)/tau)
			}
		}
		if steps < 100 || checked < 3000 {
			t.Fatalf("tau %g: only %d level steps, %d samples checked", tau, steps, checked)
		}
	}
}

// TestLevelTauWithNoise checks that with noise on the smoothed level is the
// per-sample mean: averaged over many identical steps, the first samples
// after a step follow the exponential.
func TestLevelTauWithNoise(t *testing.T) {
	p := cleanProfile(t)
	p.DwellCV = 0
	p.SpeedBPS = p.SamplingRate / 20 // 20 samples per k-mer
	p.LevelTau = 3
	m := testModel(t)
	// Alternate homopolymer runs of A and C: each run's k-mers end in the
	// homopolymer k-mer, so every A→C step has the same size.
	var seq []byte
	for i := 0; i < 400; i++ {
		b := "AC"[i%2]
		for j := 0; j < 12; j++ {
			seq = append(seq, b)
		}
	}
	rd, err := Synthesize(testRand(8), seq, p, m, WithEvents())
	if err != nil {
		t.Fatal(err)
	}
	pa := rd.PA()
	aIdx := KmerIndex([]byte("AAAAA"))
	ev := kmerEvents(rd)
	// After a run of AAAAA k-mers the next k-mer starts from the AAAAA
	// level (20 samples ≫ τ).
	sum := make([]float64, 6)
	n := 0
	var prevLevel, target float64
	for j := 1; j < len(ev); j++ {
		if ev[j-1].Kmer != aIdx || ev[j].Kmer == aIdx {
			continue
		}
		prevLevel, target = float64(m.Mean[aIdx]), ev[j].LevelPA
		// Use only steps into the first k-mer after the run, normalised
		// to a unit step.
		for i := range sum {
			sum[i] += (pa[ev[j].Start+i] - target) / (prevLevel - target)
		}
		n++
	}
	if n < 150 {
		t.Fatalf("only %d steps out of an A run", n)
	}
	// Per-sample noise SD ≤ 2.2 pA (testModel); the normalised mean has SE
	// sd/(|step|·sqrt(n)).
	step := math.Abs(prevLevel - target)
	if step < 5 {
		t.Fatalf("A-run step %.2f pA too small to resolve the approach", step)
	}
	for i := range sum {
		got := sum[i] / float64(n)
		want := math.Exp(-float64(i+1) / p.LevelTau)
		if se := 2.2 / (step * math.Sqrt(float64(n))); math.Abs(got-want) > 5*se+0.01 {
			t.Errorf("sample %d after the step: normalised mean %.4f, want %.4f (SE %.4f)", i, got, want, se)
		}
	}
}

func TestLevelTauValidate(t *testing.T) {
	for _, v := range []float64{-1, math.NaN(), math.Inf(1)} {
		p := cleanProfile(t)
		if err := p.Set("level_tau", v); err != nil {
			t.Fatal(err)
		}
		if err := p.Validate(); err == nil {
			t.Errorf("level_tau %g validated", v)
		}
		if _, err := Synthesize(testRand(1), []byte("ACGTACGTACGT"), p, testModel(t)); err == nil {
			t.Errorf("level_tau %g synthesised", v)
		}
	}
	for _, n := range PresetNames() {
		if p, _ := Preset(n); !(p.LevelTau >= 0) {
			t.Errorf("%s: level_tau %g", n, p.LevelTau)
		}
	}
}

// TestModelStdv checks that a profile's model_stdv is the stdv of a
// two-column table's k-mers (in table units, scaled like the levels), that
// a table's own stdv column wins, that the DNA preset keeps
// DefaultModelStdv (so its output is unchanged) and that Validate rejects
// bad values.
func TestModelStdv(t *testing.T) {
	dna, _ := Preset(PresetDNAR10)
	if dna.ModelStdv != DefaultModelStdv {
		t.Errorf("DNA preset model_stdv %g, want DefaultModelStdv %g", dna.ModelStdv, DefaultModelStdv)
	}
	rna, _ := Preset(PresetRNA004)
	// 0.176 standardised units x 17 pA: the ~3.0 pA per-sample noise of
	// real RNA004 reads.
	if got := rna.ModelStdv * rna.ModelScale; math.Abs(got-2.992) > 1e-9 {
		t.Errorf("RNA004 preset stdv %g pA, want 2.992", got)
	}
	two := testModelText(9, 2, false)
	rna.ModelK = 9
	m, err := ReadPoreModel(strings.NewReader(two), rna.ModelOptions()...)
	if err != nil {
		t.Fatal(err)
	}
	if got := float64(m.Stdv[12345]); math.Abs(got-0.176*17) > 1e-4 {
		t.Errorf("two-column table stdv %g pA, want %g", got, 0.176*17)
	}
	if err := rna.Set("model_stdv", 0.05); err != nil {
		t.Fatal(err)
	}
	if m, _ = ReadPoreModel(strings.NewReader(two), rna.ModelOptions()...); math.Abs(float64(m.Stdv[7])-0.05*17) > 1e-4 {
		t.Errorf("model_stdv=0.05: stdv %g pA", m.Stdv[7])
	}
	three := testModelText(9, 3, false)
	m3, err := ReadPoreModel(strings.NewReader(three), WithDefaultStdv(rna.ModelStdv))
	if err != nil {
		t.Fatal(err)
	}
	ref, _ := ReadPoreModel(strings.NewReader(three))
	if m3.Stdv[99] != ref.Stdv[99] {
		t.Errorf("table stdv column overridden: %g, want %g", m3.Stdv[99], ref.Stdv[99])
	}
	for _, v := range []float64{-0.1, math.NaN(), math.Inf(1)} {
		p := rna
		p.ModelStdv = v
		if err := p.Validate(); err == nil {
			t.Errorf("model_stdv %g validated", v)
		}
	}
}
