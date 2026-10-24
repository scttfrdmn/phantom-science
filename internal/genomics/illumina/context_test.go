//go:build !release

package illumina

import (
	"math"
	"math/rand/v2"
	"strconv"
	"strings"
	"testing"
)

// Read cells of class_reads.tsv: GC bin floor(20·GC) (GC among the
// read's A, C, G and T bases) and longest homopolymer, clamped to 3-20.
const (
	gcBins          = 20
	minHomopolymer  = 3
	maxHomopolymer  = 20
	homopolymerBins = maxHomopolymer - minHomopolymer + 1
)

// readBins returns the cell of seq (GC bin, homopolymer bin), or ok false
// for a read without A, C, G or T bases.
func readBins(seq []byte) (g, h int, ok bool) {
	var gc, n, run, longest int
	last := uint8(4)
	for _, b := range seq {
		s := contextSymbol[b]
		switch {
		case s == 4:
			run = 0
		case s == last:
			run++
		default:
			run = 1
		}
		last = s
		longest = max(longest, run)
		if s == 1 || s == 2 {
			gc++
		}
		if s < 4 {
			n++
		}
	}
	if n == 0 {
		return 0, 0, false
	}
	return gc * gcBins / n, min(max(longest, minHomopolymer), maxHomopolymer) - minHomopolymer, true
}

// readCells is the measured distribution of reads over the
// ReadClassWeights cells (GC bin, homopolymer bin) of each mate, from
// class_reads.tsv, indexed g·homopolymerBins + h.
var readCells [2][(gcBins + 1) * homopolymerBins]float64

func loadReadCells(t *testing.T) {
	t.Helper()
	readCells = [2][(gcBins + 1) * homopolymerBins]float64{}
	for _, r := range readTSV(t, "class_reads.tsv") {
		g := int(num(t, r["gc_bin"]))
		h := min(max(int(num(t, r["homopolymer"])), minHomopolymer), maxHomopolymer) - minHomopolymer
		readCells[int(num(t, r["mate"]))-1][g*homopolymerBins+h] += num(t, r["reads"])
	}
}

// humanBasesMate returns n bases with the 4-base context frequencies
// (cur.freq) of the measured reads and the measured spread of read GC
// and longest homopolymer of mate (0 or 1): a cell (GC bin, homopolymer
// bin) is drawn from the measured distribution, and the read is drawn from
// a third-order Markov chain with the context frequencies, its G and C
// steps weighted to that GC, with a homopolymer of the cell's length (A or
// T with probability 0.93, as measured for runs of 10 or more) placed at a
// random position, until its cell is the drawn one (at most 30 tries).
func humanBasesMate(r *rand.Rand, n, mate int) []byte {
	cell := pickIndex(r, readCells[mate][:])
	g, h := cell/homopolymerBins, cell%homopolymerBins
	gc := min(max((float64(g)+0.5)/gcBins, 0.01), 0.99)
	const g0 = 0.41
	tilt := gc / (1 - gc) * (1 - g0) / g0
	var seq []byte
	for range 30 {
		seq = markovBases(r, n, tilt)
		if run := h + minHomopolymer; h > 0 && run < n {
			b := "AT"[r.IntN(2)]
			if r.Float64() >= 0.93 {
				b = "CG"[r.IntN(2)]
			}
			pos := r.IntN(n - run + 1)
			for i := pos; i < pos+run; i++ {
				seq[i] = b
			}
		}
		if gb, hb, _ := readBins(seq); gb == g && hb == h {
			break
		}
	}
	return seq
}

// markovBases returns n bases of the third-order Markov chain of
// cur.freq with the weights of G and C steps multiplied by tilt.
func markovBases(r *rand.Rand, n int, tilt float64) []byte {
	seq := make([]byte, n)
	w := func(v [4]float64) []float64 {
		v[1] *= tilt
		v[2] *= tilt
		return v[:]
	}
	// The first three bases from the 3-mer marginal.
	var tri [64]float64
	freq := cur.freq
	for x, p := range freq {
		tri[x>>2] += p
	}
	x := pickIndex(r, tri[:])
	for i := 0; i < min(3, n); i++ {
		seq[i] = acgt[x>>(2*(2-i))&3]
	}
	for i := 3; i < n; i++ {
		prev := x & 63
		x = prev<<2 | pickIndex(r, w([4]float64(freq[prev<<2:prev<<2+4])))
		seq[i] = acgt[x&3]
	}
	return seq
}

// pickIndex draws an index of w with probability proportional to w.
func pickIndex(r *rand.Rand, w []float64) int {
	s := 0.0
	for _, v := range w {
		s += v
	}
	u := r.Float64() * s
	for i, v := range w {
		if u < v {
			return i
		}
		u -= v
	}
	return len(w) - 1
}

// contextCounts are quality-bin counts by context: per mate, by 4-base
// context index (interior cycles) and by called base.
type contextCounts struct {
	ctx  [2][256][4]float64
	base [2][4][4]float64
	// same and from count adjacent-cycle repeats of each bin.
	same, from [2][4]float64
}

// realContext loads the measured counts of context_quals.tsv.
func realContext(t *testing.T) *contextCounts {
	t.Helper()
	cc := &contextCounts{}
	for _, r := range readTSV(t, "context_quals.tsv") {
		m := int(num(t, r["mate"])) - 1
		var v [4]float64
		for b, q := range Bins {
			v[b] = count(t, r["Q"+strconv.Itoa(int(q))])
		}
		if x := ContextIndex(r["prev2"][0], r["prev"][0], r["base"][0], r["next"][0]); x >= 0 {
			addBins(&cc.ctx[m][x], v)
		}
		if i := baseIndex[r["base"][0]]; i >= 0 {
			addBins(&cc.base[m][i], v)
		}
	}
	return cc
}

func addBins(a *[4]float64, v [4]float64) {
	for i := range a {
		a[i] += v[i]
	}
}

// genContext tallies the qualities QualitiesFor gives humanBases reads.
func genContext(p QualityParams, m int, seed uint64, reads, n int, cc *contextCounts) {
	r, sr := newRand(seed), newRand(seed+1)
	for range reads {
		seq := humanBasesMate(sr, n, m)
		q := p.QualitiesFor(r, seq)
		for c, v := range q {
			var one [4]float64
			one[binIndex(v)] = 1
			if c >= 2 && c+1 < n {
				addBins(&cc.ctx[m][ContextIndex(seq[c-2], seq[c-1], seq[c], seq[c+1])], one)
			}
			addBins(&cc.base[m][baseIndex[seq[c]]], one)
			if c > 0 {
				j := binIndex(q[c-1])
				cc.from[m][j]++
				if q[c-1] == v {
					cc.same[m][j]++
				}
			}
		}
	}
}

func meanQ(v [4]float64) float64 {
	s, n := 0.0, 0.0
	for i, x := range v {
		s += float64(Bins[i]) * x
		n += x
	}
	return s / n
}

func binFrac(v [4]float64, i int) float64 { return v[i] / (v[0] + v[1] + v[2] + v[3]) }

// Generated qualities depend on the base and its context as measured in
// the profile's run (#35; the tolerances below are DRR276626's, see
// contextTol for the other runs). Per base: mean Q within 0.05 (read 1) and 0.09 (read 2)
// of the measured value, P(Q11) within 10 % and 12 %, P(Q25) within 4 %,
// and the spread of mean Q over the four bases within 25 % of the measured
// spread (0.27 and 0.35; the model without contexts has none). Over the
// 256 4-base contexts the generated P(Q25) follows the measured one
// (r > 0.95 and 0.88, mean relative difference < 8 %). P(Q11) follows it
// less closely (r > 0.6 and 0.25, mean relative difference < 16 % and
// 26 %): most Q11 calls are in poor reads, and in the measured reads poor
// reads are commoner where the sequence is GC-rich or has a long
// homopolymer (and, on read 2, poly-G), which the model, whose read
// classes do not depend on the sequence, does not reproduce
// (class_reads.tsv). The test reads (humanBasesMate) have the measured
// 4-mer composition and spread of read GC and longest homopolymer, and no
// sequencing errors, so their called bases are the true bases.
func TestCalibration_ContextMatchesNovaSeq(t *testing.T) {
	forProfiles(t, testContextMatches)
}

func testContextMatches(t *testing.T, pr *Profile) {
	real := realContext(t)
	gen := &contextCounts{}
	const reads = 100000
	for m, p := range []QualityParams{pr.Read1, pr.Read2} {
		genContext(p, m, uint64(20+m), reads, pr.Run.ReadLength, gen)
	}
	tol := contextTol(pr)
	for m := range 2 {
		lo, hi := [2]float64{99, 99}, [2]float64{}
		for b := range 4 {
			rv, gv := real.base[m][b], gen.base[m][b]
			t.Logf("read %d %c: mean Q %.3f (measured %.3f), P(Q11) %.5f (%.5f), P(Q25) %.5f (%.5f)", m+1, acgt[b],
				meanQ(gv), meanQ(rv), binFrac(gv, 1), binFrac(rv, 1), binFrac(gv, 2), binFrac(rv, 2))
			if math.Abs(meanQ(gv)-meanQ(rv)) > tol[m].meanQ {
				t.Errorf("read %d %c: mean Q %.3f, measured %.3f", m+1, acgt[b], meanQ(gv), meanQ(rv))
			}
			for i, rel := range [3]float64{0, tol[m].low, tol[m].mid} {
				if i == 0 {
					continue
				}
				if g, w := binFrac(gv, i), binFrac(rv, i); math.Abs(g-w) > rel*w {
					t.Errorf("read %d %c: P(Q%d) %.5f, measured %.5f", m+1, acgt[b], Bins[i], g, w)
				}
			}
			for i, v := range [2]float64{meanQ(rv), meanQ(gv)} {
				lo[i], hi[i] = min(lo[i], v), max(hi[i], v)
			}
		}
		spreadReal, spreadGen := hi[0]-lo[0], hi[1]-lo[1]
		t.Logf("read %d: spread of mean Q over bases %.3f (measured %.3f)", m+1, spreadGen, spreadReal)
		if math.Abs(spreadGen-spreadReal) > 0.25*spreadReal {
			t.Errorf("read %d: spread of mean Q over bases %.3f, measured %.3f", m+1, spreadGen, spreadReal)
		}
		for i := 1; i <= 2; i++ {
			var sx, sy, sxx, syy, sxy, sw, dev float64
			for x := range 256 {
				rv, gv := real.ctx[m][x], gen.ctx[m][x]
				w := rv[0] + rv[1] + rv[2] + rv[3]
				a, b := binFrac(rv, i), binFrac(gv, i)
				sx += w * a
				sy += w * b
				sxx += w * a * a
				syy += w * b * b
				sxy += w * a * b
				sw += w
				dev += w * math.Abs(b-a) / a
			}
			mx, my := sx/sw, sy/sw
			corr := (sxy/sw - mx*my) / math.Sqrt((sxx/sw-mx*mx)*(syy/sw-my*my))
			dev /= sw
			t.Logf("read %d Q%d by 4-base context: r = %.3f, mean |relative difference| %.3f", m+1, Bins[i], corr, dev)
			minR, maxDev := tol[m].midR, 0.08
			if i == 1 {
				minR, maxDev = tol[m].lowR, tol[m].lowDev
			}
			if corr < minR || dev > maxDev {
				t.Errorf("read %d Q%d by context: r = %.3f, mean relative difference %.3f", m+1, Bins[i], corr, dev)
			}
		}
	}
}

// ctxTol are the tolerances of TestCalibration_ContextMatchesNovaSeq for
// one mate: of mean Q per base, P(Q11) and P(Q25) per base (relative), and
// of P(Q11) (r, mean relative difference) and P(Q25) (r) over contexts.
type ctxTol struct{ meanQ, low, mid, lowR, lowDev, midR float64 }

// contextTol returns the tolerances of profile p, per mate. The other
// runs' read 2 has poorer G calls than the model gives it (mean Q of G
// 0.20 lower in DRR394916, 0.33 in ERR15600747, P(Q11) 11 % and 15 %
// higher): like poly-G and the GC-rich reads, a read-level effect that the
// model's classes, which do not depend on the sequence, do not carry; its
// P(Q25) by context follows the measured one less closely too.
func contextTol(p *Profile) [2]ctxTol {
	if p == &profileDRR276626 {
		return [2]ctxTol{{0.05, 0.10, 0.04, 0.6, 0.16, 0.95}, {0.09, 0.12, 0.04, 0.25, 0.26, 0.88}}
	}
	return [2]ctxTol{{0.07, 0.10, 0.04, 0.6, 0.16, 0.95}, {0.4, 0.18, 0.05, 0.25, 0.26, 0.78}}
}

// QualitiesFor draws each cycle's bin with the probability TargetFor gives
// for the read's bases (within 4σ of the binomial error plus 0.002 for the
// sticky chain's lag behind a target that now changes with the context),
// for reads of every composition: the measured mix, poly-A, GC-rich and
// with N.
func TestQualitiesFor_FrequenciesMatchTarget(t *testing.T) {
	loadReadCells(t)
	const reads = 40000
	sr := newRand(5)
	gc := []byte("GCGGCCGCGGGCCCGCGCGGCGCCGGGCGCGCCCGG")
	for _, tc := range []struct {
		name string
		p    QualityParams
		seq  []byte
	}{
		{"read 1 sample", Read1, humanBasesMate(sr, 151, 0)},
		{"read 2 sample", Read2, humanBasesMate(sr, 151, 1)},
		{"read 1 poly-A", Read1, []byte(strings.Repeat("A", 60) + string(humanBasesMate(sr, 91, 0)))},
		{"read 2 GC-rich", Read2, []byte(strings.Repeat(string(gc), 5)[:151])},
		{"read 1 with N, 40 bases", Read1, []byte("ACGTNACGTTTGGGNNCCAATTGGCCAAGGTTCCAGTCGA")},
	} {
		n := len(tc.seq)
		counts := make([][4]float64, n)
		r := newRand(6)
		for range reads {
			for c, v := range tc.p.QualitiesFor(r, tc.seq) {
				counts[c][binIndex(v)]++
			}
		}
		bad := 0
		for c := range n {
			want := tc.p.TargetFor(tc.seq, c)
			if s := want[0] + want[1] + want[2] + want[3]; math.Abs(s-1) > 1e-12 {
				t.Errorf("%s cycle %d: target sums to %v", tc.name, c+1, s)
			}
			for i := 1; i < 4; i++ {
				got := counts[c][i] / reads
				if math.Abs(got-want[i]) > 4*math.Sqrt(want[i]*(1-want[i])/reads)+0.002 {
					bad++
					t.Errorf("%s cycle %d (%q): P(Q%d) = %.4f, target %.4f", tc.name, c+1, tc.seq[max(c-2, 0):min(c+2, n)], Bins[i], got, want[i])
				}
			}
			if bad > 5 {
				t.FailNow()
			}
		}
	}
}

// The context factors average to 1 over the measured context frequencies
// (so a read of the sample's composition keeps the per-cycle fractions);
// a context with unknown bases gets the frequency-weighted mean of the
// contexts it matches, and a fully unknown one the overall mean, 1.
func TestContextFactors_Normalised(t *testing.T) {
	for _, pr := range profiles {
		t.Run(pr.Name, func(t *testing.T) { testContextNormalised(t, pr) })
	}
	if ContextIndex('A', 'C', 'G', 'T') != 0<<6|1<<4|2<<2|3 || ContextIndex('A', 'C', 'N', 'T') != -1 {
		t.Error("ContextIndex")
	}
}

func testContextNormalised(t *testing.T, pr *Profile) {
	contextFreq := pr.freq
	sum := 0.0
	for _, v := range contextFreq {
		sum += v
	}
	if math.Abs(sum-1) > 1e-3 {
		t.Errorf("context frequencies sum to %v", sum)
	}
	for mate, p := range []QualityParams{pr.Read1, pr.Read2} {
		for k, c := range p.Classes {
			f := c.Context
			if f == nil {
				t.Fatalf("read %d class %d has no context factors", mate+1, k)
			}
			var l, m float64
			for x := range 256 {
				l += contextFreq[x] * f.Low[x]
				m += contextFreq[x] * f.Mid[x]
				if f.Low[x] <= 0 || f.Mid[x] <= 0 {
					t.Errorf("read %d class %d context %d: factors %v %v", mate+1, k, x, f.Low[x], f.Mid[x])
				}
			}
			if math.Abs(l/sum-1) > 5e-3 || math.Abs(m/sum-1) > 5e-3 {
				t.Errorf("read %d class %d: factors average %.5f, %.5f; want 1", mate+1, k, l/sum, m/sum)
			}
			// A known context maps to the table; unknown positions average.
			x := ContextIndex('G', 'G', 'T', 'A')
			if got := f.ext[contextAt([]byte("GGTA"), 2)]; math.Abs(got[0]-f.Low[x]) > 1e-12*f.Low[x] || math.Abs(got[1]-f.Mid[x]) > 1e-12*f.Mid[x] {
				t.Errorf("read %d class %d: GGTA factors %v, table %v %v", mate+1, k, got, f.Low[x], f.Mid[x])
			}
			if got := f.ext[contextAt([]byte("N"), 0)]; math.Abs(got[0]-l/sum) > 1e-9 || math.Abs(got[1]-m/sum) > 1e-9 {
				t.Errorf("read %d class %d: unknown context factors %v, want %v %v", mate+1, k, got, l/sum, m/sum)
			}
			var w, acc float64
			for b := range 4 {
				y := ContextIndex('G', 'G', 'T', acgt[b])
				w += contextFreq[y]
				acc += contextFreq[y] * f.Low[y]
			}
			if got := f.ext[contextAt([]byte("GGT"), 2)][0]; math.Abs(got-acc/w) > 1e-12 {
				t.Errorf("read %d class %d: GGT- factor %v, want %v", mate+1, k, got, acc/w)
			}
		}
	}
}
