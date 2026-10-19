//go:build !release

package signal

import (
	"math"
	"math/rand/v2"
	"slices"
	"sync"
	"testing"
)

func testRand(seed uint64) *rand.Rand { return rand.New(rand.NewPCG(seed, 11)) }

func randomSeq(r *rand.Rand, n int) []byte {
	b := make([]byte, n)
	for i := range b {
		b[i] = "ACGT"[r.IntN(4)]
	}
	return b
}

// cleanProfile is the R10 preset without adapter, drift, read-level
// jitter, stalls or skips, for isolating one effect at a time.
func cleanProfile(t testing.TB) Profile {
	t.Helper()
	p, err := Preset(PresetDNAR10)
	if err != nil {
		t.Fatal(err)
	}
	p.AdapterSeconds = 0
	p.DriftPA = 0
	p.ReadShiftSD, p.ReadScaleSD = 0, 0
	p.StallProb, p.SkipProb = 0, 0
	p.ModelK = 5
	return p
}

func meanVar(xs []float64) (mean, variance float64) {
	for _, x := range xs {
		mean += x
	}
	mean /= float64(len(xs))
	for _, x := range xs {
		variance += (x - mean) * (x - mean)
	}
	return mean, variance / float64(len(xs)-1)
}

func kmerEvents(rd *Read) []Event {
	var out []Event
	for _, e := range rd.Events {
		if e.Kind == EventKmer {
			out = append(out, e)
		}
	}
	return out
}

// TestDwellStatistics checks that samples per base average
// SamplingRate/SpeedBPS and have the configured gamma CV.
func TestDwellStatistics(t *testing.T) {
	model := testModel(t)
	for _, cv := range []float64{0.5, 0.7} {
		p := cleanProfile(t)
		p.DwellCV = cv
		r := testRand(1)
		rd, err := Synthesize(r, randomSeq(r, 40000), p, model, WithEvents())
		if err != nil {
			t.Fatal(err)
		}
		ev := kmerEvents(rd)
		if len(ev) != rd.Kmers || rd.Kmers != 40000-5+1 {
			t.Fatalf("events %d, kmers %d", len(ev), rd.Kmers)
		}
		lens := make([]float64, len(ev))
		next := rd.SeqStart
		for i, e := range ev {
			if e.Start != next {
				t.Fatalf("event %d starts at %d, want contiguous %d", i, e.Start, next)
			}
			next += e.Len
			lens[i] = float64(e.Len)
		}
		if next != rd.SeqEnd || rd.SeqEnd != rd.NumSamples {
			t.Fatalf("events end at %d, SeqEnd %d, NumSamples %d", next, rd.SeqEnd, rd.NumSamples)
		}
		mu := p.MeanDwell()
		// Samples per k-mer = gamma dwell + discretisation (variance ≈ 1/6).
		wantVar := mu*mu*cv*cv + 1.0/6
		mean, v := meanVar(lens)
		n := float64(len(lens))
		if d := math.Abs(mean - mu); d > 4*math.Sqrt(wantVar/n) {
			t.Errorf("cv=%g: mean dwell %.4f, want %.4f ± %.4f", cv, mean, mu, 4*math.Sqrt(wantVar/n))
		}
		// SE of the sample variance: var·sqrt((2 + excess kurtosis)/n);
		// gamma excess kurtosis is 6/shape = 6cv².
		se := wantVar * math.Sqrt((2+6*cv*cv)/n)
		if d := math.Abs(v - wantVar); d > 4*se {
			t.Errorf("cv=%g: dwell variance %.3f, want %.3f ± %.3f (measured CV %.3f)", cv, v, wantVar, 4*se, math.Sqrt(v)/mean)
		}
		t.Logf("cv=%g: mean dwell %.3f (want %.3f), measured CV %.4f", cv, mean, mu, math.Sqrt(v)/mean)
	}
}

// TestLevelStatistics checks segment means against the model levels and
// segment spread against the configured noise, using the ground truth.
func TestLevelStatistics(t *testing.T) {
	model := testModel(t)
	p := cleanProfile(t)
	p.WhiteNoisePA = 1.5
	r := testRand(2)
	seq := randomSeq(r, 8000)
	rd, err := Synthesize(r, seq, p, model, WithEvents())
	if err != nil {
		t.Fatal(err)
	}
	pa := rd.PA()
	q := rd.Range / rd.Digitisation // quantisation step
	var zs []float64
	var ss, df float64
	maxZ := 0.0
	for _, e := range kmerEvents(rd) {
		if want := KmerIndex(seq[e.Pos : e.Pos+5]); e.Kmer != want {
			t.Fatalf("event at %d: kmer %d, want %d", e.Pos, e.Kmer, want)
		}
		wantLevel := float64(model.Mean[e.Kmer])
		lsd := float64(model.Stdv[e.Kmer])
		wantSD := math.Sqrt(lsd*lsd + p.WhiteNoisePA*p.WhiteNoisePA)
		if math.Abs(e.LevelPA-wantLevel) > 1e-9 || math.Abs(e.NoisePA-wantSD) > 1e-9 {
			t.Fatalf("event truth level %g sd %g, want %g %g", e.LevelPA, e.NoisePA, wantLevel, wantSD)
		}
		if e.Len < 2 {
			continue
		}
		seg := pa[e.Start : e.Start+e.Len]
		m, v := meanVar(seg)
		sd := math.Sqrt(wantSD*wantSD + q*q/12)
		z := (m - wantLevel) / (sd / math.Sqrt(float64(e.Len)))
		zs = append(zs, z)
		maxZ = max(maxZ, math.Abs(z))
		ss += v / (sd * sd) * float64(e.Len-1)
		df += float64(e.Len - 1)
	}
	zm, zv := meanVar(zs)
	n := float64(len(zs))
	if math.Abs(zm) > 4/math.Sqrt(n) {
		t.Errorf("mean standardised segment residual %.4f, want 0 ± %.4f", zm, 4/math.Sqrt(n))
	}
	if math.Abs(zv-1) > 4*math.Sqrt(2/(n-1)) {
		t.Errorf("variance of standardised residuals %.4f, want 1 ± %.4f", zv, 4*math.Sqrt(2/(n-1)))
	}
	// Max |z| over n normal draws: P(|z| > 5.5) ≈ 4e-8 each.
	if maxZ > 5.5 {
		t.Errorf("max |z| = %.2f over %d segments", maxZ, len(zs))
	}
	// Pooled within-segment variance / configured variance ≈ 1.
	ratio := ss / df
	if math.Abs(ratio-1) > 4*math.Sqrt(2/df) {
		t.Errorf("pooled variance ratio %.4f, want 1 ± %.4f", ratio, 4*math.Sqrt(2/df))
	}
	t.Logf("%d segments: z mean %.4f var %.4f max %.2f; pooled var ratio %.4f", len(zs), zm, zv, maxZ, ratio)
}

func TestReadLevelScaleShiftInTruth(t *testing.T) {
	model := testModel(t)
	p := cleanProfile(t)
	p.ReadShiftSD, p.ReadScaleSD = 5, 0.05
	r := testRand(3)
	rd, err := Synthesize(r, randomSeq(r, 2000), p, model, WithEvents())
	if err != nil {
		t.Fatal(err)
	}
	ev := kmerEvents(rd)
	// LevelPA = mean·scale + shift for one (scale, shift) per read: fit it.
	a, b := ev[0], ev[len(ev)-1]
	for _, e := range ev {
		if math.Abs(float64(model.Mean[e.Kmer]-model.Mean[a.Kmer])) > 10 {
			b = e
			break
		}
	}
	scale := (b.LevelPA - a.LevelPA) / float64(model.Mean[b.Kmer]-model.Mean[a.Kmer])
	shift := a.LevelPA - float64(model.Mean[a.Kmer])*scale
	if scale == 1 && shift == 0 {
		t.Error("read-level scale/shift not applied")
	}
	for _, e := range ev {
		if want := float64(model.Mean[e.Kmer])*scale + shift; math.Abs(e.LevelPA-want) > 1e-3 {
			t.Fatalf("level %g, want %g (scale %g shift %g)", e.LevelPA, want, scale, shift)
		}
	}
}

func TestDigitiserRoundTripAndClamp(t *testing.T) {
	p, _ := Preset(PresetDNAR10)
	d := p.Digitiser(-127.5)
	tol := d.Range / (2 * d.Digitisation)
	r := testRand(4)
	for i := 0; i < 100000; i++ {
		x := -50 + 400*r.Float64()
		raw, clipped := d.Raw(x)
		if clipped {
			t.Fatalf("%g pA clipped", x)
		}
		if e := math.Abs(d.PA(raw) - x); e > tol+1e-9 {
			t.Fatalf("round trip %g pA → %d → %g: error %g > %g", x, raw, d.PA(raw), e, tol)
		}
	}
	// The formula itself, at a known point: raw = pA·dig/range − offset.
	if raw, _ := d.Raw(100); raw != int16(math.Round(100*2048/281.345551+127.5)) {
		t.Errorf("Raw(100) = %d", raw)
	}
	for _, c := range []struct {
		pa   float64
		want int16
	}{{1e7, math.MaxInt16}, {-1e7, math.MinInt16}} {
		if raw, clipped := d.Raw(c.pa); raw != c.want || !clipped {
			t.Errorf("Raw(%g) = %d clipped=%v, want %d clipped", c.pa, raw, clipped, c.want)
		}
	}

	// A profile whose range is far too small clamps in synthesis too.
	model := testModel(t)
	cp := cleanProfile(t)
	cp.Range = 2 // 1024 raw units per pA: levels of ~85 pA overflow
	rd, err := Synthesize(testRand(5), randomSeq(testRand(6), 200), cp, model)
	if err != nil {
		t.Fatal(err)
	}
	high := 0
	for _, v := range rd.Raw {
		if v == math.MaxInt16 {
			high++
		}
	}
	if rd.Clipped == 0 || high != rd.Clipped {
		t.Errorf("Clipped = %d, samples at MaxInt16 = %d", rd.Clipped, high)
	}
}

func TestStallSkipRates(t *testing.T) {
	model := testModel(t)
	p := cleanProfile(t)
	p.StallProb, p.SkipProb = 0.03, 0.05
	r := testRand(7)
	rd, err := Synthesize(r, randomSeq(r, 30000), p, model, WithEvents())
	if err != nil {
		t.Fatal(err)
	}
	ev := kmerEvents(rd)
	n := float64(len(ev))
	var stalls, skips int
	var stallLens []float64
	for _, e := range ev {
		switch {
		case e.Skip:
			skips++
			if e.Len != 0 {
				t.Fatalf("skipped k-mer at %d has %d samples", e.Pos, e.Len)
			}
		case e.Stall:
			stalls++
			stallLens = append(stallLens, float64(e.Len))
		}
	}
	if stalls != rd.Stalls || skips != rd.Skips {
		t.Errorf("Read counts stalls %d skips %d, events %d %d", rd.Stalls, rd.Skips, stalls, skips)
	}
	for _, c := range []struct {
		name string
		got  int
		p    float64
	}{{"stall", stalls, p.StallProb}, {"skip", skips, p.SkipProb}} {
		want := n * c.p
		if d := math.Abs(float64(c.got) - want); d > 4*math.Sqrt(n*c.p*(1-c.p)) {
			t.Errorf("%s count %d, want %.0f ± %.0f", c.name, c.got, want, 4*math.Sqrt(n*c.p*(1-c.p)))
		}
	}
	// Stall dwell = gamma × U(5, 20): mean 12.5·μ.
	mu := p.MeanDwell()
	mult := (p.StallMin + p.StallMax) / 2
	// Var(XY) for independent X (gamma) and Y (uniform).
	ex2 := mu * mu * (1 + p.DwellCV*p.DwellCV)
	ey2 := mult*mult + (p.StallMax-p.StallMin)*(p.StallMax-p.StallMin)/12
	sd := math.Sqrt(ex2*ey2 - mu*mu*mult*mult)
	m, _ := meanVar(stallLens)
	if d := math.Abs(m - mu*mult); d > 4*sd/math.Sqrt(float64(len(stallLens))) {
		t.Errorf("mean stall length %.1f, want %.1f ± %.1f", m, mu*mult, 4*sd/math.Sqrt(float64(len(stallLens))))
	}
	t.Logf("%d k-mers: %d stalls (%.4f), %d skips (%.4f), mean stall %.1f samples", len(ev), stalls, float64(stalls)/n, skips, float64(skips)/n, m)
}

func TestDeterminismPerReadID(t *testing.T) {
	model := testModel(t)
	p, _ := Preset(PresetDNAR10)
	p.ModelK = 5
	seqA := randomSeq(testRand(8), 3000)
	seqB := randomSeq(testRand(9), 3000)
	const seed = 42
	idA, idB := "0b9d7e5c-3c55-4f43-9d3b-6f4c6c7f3a10", "8a2b1f44-1e0d-4c43-a1bb-92d52e8e7c01"

	gen := func(seed uint64, id string, seq []byte) *Read {
		rd, err := SynthesizeRead(seed, id, seq, p, model)
		if err != nil {
			t.Fatal(err)
		}
		return rd
	}
	a1 := gen(seed, idA, seqA)
	_ = gen(seed, idB, seqB) // interleave another read
	a2 := gen(seed, idA, seqA)
	if !slices.Equal(a1.Raw, a2.Raw) || a1.Offset != a2.Offset || a1.ID != idA {
		t.Fatal("same (seed, read ID) gave different signal")
	}
	b1 := gen(seed, idB, seqB)
	b2 := gen(seed, idB, seqB)
	if !slices.Equal(b1.Raw, b2.Raw) {
		t.Fatal("order dependence: read B differs between runs")
	}
	if c := gen(seed, idB, seqA); slices.Equal(c.Raw, a1.Raw) {
		t.Error("different read ID gave identical signal for the same sequence")
	}
	if c := gen(seed+1, idA, seqA); slices.Equal(c.Raw, a1.Raw) {
		t.Error("different run seed gave identical signal")
	}
	// Recording events must not change the signal.
	e, err := SynthesizeRead(seed, idA, seqA, p, model, WithEvents())
	if err != nil {
		t.Fatal(err)
	}
	if !slices.Equal(e.Raw, a1.Raw) {
		t.Error("WithEvents changed the signal")
	}
}

func TestConcurrentSynthesisRace(t *testing.T) {
	model := testModel(t)
	p, _ := Preset(PresetDNAR10)
	p.ModelK = 5
	const n = 16
	seqs := make([][]byte, n)
	ids := make([]string, n)
	r := testRand(10)
	for i := range seqs {
		seqs[i] = randomSeq(r, 1500)
		ids[i] = string(rune('a'+i)) + "-read"
	}
	serial := make([][]int16, n)
	for i := range seqs {
		rd, err := SynthesizeRead(7, ids[i], seqs[i], p, model)
		if err != nil {
			t.Fatal(err)
		}
		serial[i] = rd.Raw
	}
	got := make([][]int16, n)
	var wg sync.WaitGroup
	for i := n - 1; i >= 0; i-- {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			rd, err := SynthesizeRead(7, ids[i], seqs[i], p, model, WithEvents())
			if err != nil {
				t.Error(err)
				return
			}
			got[i] = rd.Raw
		}(i)
	}
	wg.Wait()
	for i := range got {
		if !slices.Equal(got[i], serial[i]) {
			t.Errorf("read %d differs when synthesised concurrently", i)
		}
	}
}

func TestSegmentsAdapterPolyADuration(t *testing.T) {
	model := testModel(t)
	p, _ := Preset(PresetRNA004)
	p.ModelK = 5
	p.DriftPA = 0
	r := testRand(11)
	seq := randomSeq(r, 1200)
	const polyA = 60
	rd, err := Synthesize(r, seq, p, model, WithPolyA(polyA), WithEvents())
	if err != nil {
		t.Fatal(err)
	}
	if rd.Duration != float64(rd.NumSamples)/rd.SamplingRate || rd.NumSamples != len(rd.Raw) {
		t.Errorf("Duration %g, NumSamples %d, len %d, rate %g", rd.Duration, rd.NumSamples, len(rd.Raw), rd.SamplingRate)
	}
	if want := int(math.Round(p.AdapterSeconds * p.SamplingRate)); rd.AdapterSamples != want {
		t.Errorf("adapter %d samples, want %d", rd.AdapterSamples, want)
	}
	wantPoly := int(math.Round(polyA * p.MeanDwell()))
	if rd.PolyStart != rd.AdapterSamples || rd.PolyEnd-rd.PolyStart != wantPoly || rd.SeqStart != rd.PolyEnd {
		t.Errorf("poly-A [%d,%d), want start %d length %d; SeqStart %d", rd.PolyStart, rd.PolyEnd, rd.AdapterSamples, wantPoly, rd.SeqStart)
	}
	pa := rd.PA()
	_, polyVar := meanVar(pa[rd.PolyStart:rd.PolyEnd])
	q := rd.Range / rd.Digitisation
	wantVar := p.PolyANoisePA*p.PolyANoisePA + q*q/12
	nPoly := float64(rd.PolyEnd - rd.PolyStart)
	if math.Abs(polyVar-wantVar) > 4*wantVar*math.Sqrt(2/nPoly) {
		t.Errorf("poly-A variance %.3f, want %.3f", polyVar, wantVar)
	}
	_, seqVar := meanVar(pa[rd.SeqStart:rd.SeqEnd])
	if polyVar*10 > seqVar {
		t.Errorf("poly-A plateau not low-variance: %.2f vs sequence %.2f", polyVar, seqVar)
	}
	// RNA translocates 3'→5': the first sequence k-mer is the last one.
	ev := kmerEvents(rd)
	if ev[0].Pos != len(seq)-5 || ev[len(ev)-1].Pos != 0 {
		t.Errorf("RNA walk order: first pos %d last %d", ev[0].Pos, ev[len(ev)-1].Pos)
	}
	if rd.Events[0].Kind != EventLeader || rd.Events[1].Kind != EventPolyTail {
		t.Errorf("event order %v %v", rd.Events[0].Kind, rd.Events[1].Kind)
	}

	// DNA: adapter k-mers follow the leader; poly-T at the end; override adapter time.
	pd, _ := Preset(PresetDNAR10)
	pd.ModelK = 5
	rd2, err := Synthesize(testRand(12), seq, pd, model, WithAdapterSeconds(0.5), WithPolyTail(PolyTail{Length: 20, Base: 'T', AtEnd: true}), WithEvents())
	if err != nil {
		t.Fatal(err)
	}
	if rd2.AdapterSamples != 2500 {
		t.Errorf("DNA adapter %d samples, want 2500", rd2.AdapterSamples)
	}
	if rd2.PolyStart != rd2.SeqEnd || rd2.PolyEnd != rd2.NumSamples || rd2.PolyEnd-rd2.PolyStart != int(math.Round(20*pd.MeanDwell())) {
		t.Errorf("poly-T at end: [%d,%d) SeqEnd %d N %d", rd2.PolyStart, rd2.PolyEnd, rd2.SeqEnd, rd2.NumSamples)
	}
	last := rd2.Events[len(rd2.Events)-1]
	if want := float64(model.Mean[KmerIndex([]byte("TTTTT"))]); last.Kind != EventPolyTail || math.Abs(last.LevelPA-want) > 20 {
		t.Errorf("poly-T level %g, want near %g", last.LevelPA, want)
	}
	nAdapter := 0
	for _, e := range rd2.Events {
		if e.Kind == EventAdapter {
			nAdapter++
		}
	}
	if nAdapter != len(kit14Adapter)-5+1 {
		t.Errorf("%d adapter k-mers, want %d", nAdapter, len(kit14Adapter)-4)
	}
	if rd2.Events[0].Kind != EventLeader || rd2.Events[0].Len+sumLen(rd2.Events, EventAdapter) != 2500 {
		t.Errorf("leader + adapter k-mers != adapter samples")
	}
	if rd2.SeqStart != rd2.AdapterSamples || rd2.PolyStart == rd2.AdapterSamples {
		t.Errorf("DNA layout: SeqStart %d Adapter %d", rd2.SeqStart, rd2.AdapterSamples)
	}
}

func sumLen(ev []Event, k EventKind) int {
	n := 0
	for _, e := range ev {
		if e.Kind == k {
			n += e.Len
		}
	}
	return n
}

func TestNonACGTBases(t *testing.T) {
	model := testModel(t)
	p := cleanProfile(t)
	seq := []byte("ACGTNNNNNacgtRYKMacguACGTACGT")
	orig := append([]byte(nil), seq...)
	a, err := Synthesize(testRand(13), seq, p, model, WithEvents())
	if err != nil {
		t.Fatal(err)
	}
	b, _ := Synthesize(testRand(13), seq, p, model)
	if !slices.Equal(a.Raw, b.Raw) {
		t.Error("N replacement not deterministic")
	}
	if !slices.Equal(seq, orig) {
		t.Error("input sequence modified")
	}
	for _, e := range kmerEvents(a) {
		if e.Kmer < 0 || e.Kmer >= model.NumKmers() {
			t.Fatalf("invalid k-mer index %d", e.Kmer)
		}
	}
	if _, err := Synthesize(testRand(1), []byte("ACG"), p, model); err == nil {
		t.Error("sequence shorter than k: want error")
	}
}

func TestDriftAmplitude(t *testing.T) {
	model := testModel(t)
	p := cleanProfile(t)
	p.DriftPA, p.DriftTau = 3, 0.002 // 10-sample correlation at 5 kHz
	p.AdapterSeconds, p.AdapterNoisePA, p.AdapterLevelPA, p.AdapterSeq = 40, 0, 100, ""
	p.Range = 2048 / 32.0 // fine quantisation (1/32 pA) for this check
	rd, err := Synthesize(testRand(14), randomSeq(testRand(15), 10), p, model)
	if err != nil {
		t.Fatal(err)
	}
	pa := rd.PA()[:rd.AdapterSamples]
	_, v := meanVar(pa)
	a := math.Exp(-1 / (p.DriftTau * p.SamplingRate))
	n := float64(len(pa))
	want := p.DriftPA * p.DriftPA
	// Var of the sample variance of an AR(1) process ≈ 2σ⁴(1+a²)/((1−a²)n).
	se := want * math.Sqrt(2*(1+a*a)/((1-a*a)*n))
	if math.Abs(v-want) > 4*se {
		t.Errorf("drift variance %.3f, want %.3f ± %.3f", v, want, 4*se)
	}
}

func TestGammaSampler(t *testing.T) {
	r := testRand(16)
	for _, c := range []struct{ mean, cv float64 }{{12.5, 0.5}, {30, 0.7}, {5, 1.5}} {
		const n = 100000
		xs := make([]float64, n)
		for i := range xs {
			xs[i] = gammaMeanCV(r, c.mean, c.cv)
			if xs[i] < 0 {
				t.Fatal("negative gamma draw")
			}
		}
		m, v := meanVar(xs)
		sd := c.mean * c.cv
		if math.Abs(m-c.mean) > 4*sd/math.Sqrt(n) {
			t.Errorf("gamma(mean %g cv %g): mean %.4f", c.mean, c.cv, m)
		}
		if se := sd * sd * math.Sqrt((2+6*c.cv*c.cv)/n); math.Abs(v-sd*sd) > 4*se {
			t.Errorf("gamma(mean %g cv %g): var %.4f want %.4f ± %.4f", c.mean, c.cv, v, sd*sd, 4*se)
		}
	}
	if gammaMeanCV(r, 7, 0) != 7 {
		t.Error("cv=0 should be deterministic")
	}
}

func TestPresetsAndValidate(t *testing.T) {
	if names := PresetNames(); len(names) != 2 || names[0] != PresetDNAR10 {
		t.Errorf("PresetNames = %v", names)
	}
	for _, n := range PresetNames() {
		p, _ := Preset(n)
		if err := p.Validate(); err != nil {
			t.Errorf("%s: %v", n, err)
		}
	}
	dna, _ := Preset(PresetDNAR10)
	rna, _ := Preset(PresetRNA004)
	if dna.MeanDwell() != 12.5 || dna.Digitisation != 2048 || dna.Reverse {
		t.Errorf("DNA preset: dwell %g dig %g reverse %v", dna.MeanDwell(), dna.Digitisation, dna.Reverse)
	}
	if math.Abs(rna.MeanDwell()-4000.0/130) > 1e-12 || !rna.Reverse {
		t.Errorf("RNA preset: dwell %g reverse %v", rna.MeanDwell(), rna.Reverse)
	}
	if _, err := Preset("r9"); err == nil {
		t.Error("unknown preset: want error")
	}
	bad := []func(*Profile){
		func(p *Profile) { p.SamplingRate = 0 },
		func(p *Profile) { p.SpeedBPS = math.Inf(1) },
		func(p *Profile) { p.DwellCV = -1 },
		func(p *Profile) { p.SkipProb = 1.5 },
		func(p *Profile) { p.StallProb, p.SkipProb = 0.6, 0.6 },
		func(p *Profile) { p.StallMin, p.StallMax = 20, 5 },
		func(p *Profile) { p.DriftPA, p.DriftTau = 1, 0 },
		func(p *Profile) { p.SpeedBPS = 1e6 },
		func(p *Profile) { p.WhiteNoisePA = math.NaN() },
	}
	for i, f := range bad {
		p := dna
		f(&p)
		if err := p.Validate(); err == nil {
			t.Errorf("bad profile %d validated", i)
		}
		if _, err := Synthesize(testRand(1), []byte("ACGTACGTACGT"), p, testModel(t)); err == nil {
			t.Errorf("bad profile %d synthesised", i)
		}
	}
	if _, err := Synthesize(testRand(1), []byte("ACGTACGTACGT"), dna, nil); err == nil {
		t.Error("nil model: want error")
	}
	if _, err := Synthesize(testRand(1), []byte("ACGTACGTACGT"), dna, testModel(t), WithPolyA(-1)); err == nil {
		t.Error("negative poly-A: want error")
	}
}
