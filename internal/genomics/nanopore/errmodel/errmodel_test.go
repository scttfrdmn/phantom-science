//go:build !release

package errmodel

import (
	"bytes"
	"math"
	"math/rand/v2"
	"sort"
	"sync"
	"testing"
)

func newRand(seed uint64) *rand.Rand { return rand.New(rand.NewPCG(seed, 0x5eed)) } // #nosec G404

// randomTruth returns n uniform ACGT bases.
func randomTruth(r *rand.Rand, n int) []byte {
	const acgt = "ACGT"
	b := make([]byte, n)
	for i := range b {
		b[i] = acgt[r.IntN(4)]
	}
	return b
}

// runsOf returns runs of exactly runLen identical bases, each run's base
// differing from its neighbours, totalling about n bases.
func runsOf(r *rand.Rand, runLen, n int) []byte {
	const acgt = "ACGT"
	b := make([]byte, 0, n+runLen)
	prev := -1
	for len(b) < n {
		k := r.IntN(4)
		for k == prev {
			k = r.IntN(4)
		}
		prev = k
		for j := 0; j < runLen; j++ {
			b = append(b, acgt[k])
		}
	}
	return b
}

// runLengths returns, per truth position, the length of its homopolymer run.
func runLengths(truth []byte) []int {
	rl := make([]int, len(truth))
	for i := 0; i < len(truth); {
		j := i + 1
		for j < len(truth) && truth[j] == truth[i] {
			j++
		}
		for k := i; k < j; k++ {
			rl[k] = j - i
		}
		i = j
	}
	return rl
}

// counts tallies ops per truth position in [lo, hi): subs and dels of truth
// bases in the window and insertions after them; n is the number of truth
// bases in the window and delP the summed deletion probability.
type counts struct {
	n, sub, ins, del int
	delP, delV       float64 // Σ p and Σ p(1−p) of deletion probabilities
}

func tally(m *Model, truth []byte, ops []Op, lo, hi int, c *counts) {
	rl := runLengths(truth)
	t := 0
	for _, op := range ops {
		if op.Kind == Ins {
			if t-1 >= lo && t-1 < hi { // insertion after truth base t-1
				c.ins++
			}
			continue
		}
		if t >= lo && t < hi {
			c.n++
			p := m.DeletionRate(rl[t])
			c.delP += p
			c.delV += p * (1 - p)
			switch op.Kind {
			case Sub:
				c.sub++
			case Del:
				c.del++
			}
		}
		t++
	}
}

// within4σ checks a binomial count against probability p over n trials.
func within4σ(t *testing.T, what string, got, n int, p float64) {
	t.Helper()
	tol := 4 * math.Sqrt(p*(1-p)/float64(n))
	rate := float64(got) / float64(n)
	if math.Abs(rate-p) > tol {
		t.Errorf("%s rate %.5f, want %.5f ± %.5f (4σ, n=%d)", what, rate, p, tol, n)
	} else {
		t.Logf("%s rate %.5f (want %.5f ± %.5f, n=%d)", what, rate, p, tol, n)
	}
}

// checkCounts verifies sub/ins/del rates of c against m's table.
func checkCounts(t *testing.T, m *Model, c counts) {
	t.Helper()
	p := m.Params()
	within4σ(t, "substitution", c.sub, c.n, p.Sub)
	within4σ(t, "insertion", c.ins, c.n, p.Ins)
	// Deletions are Poisson-binomial: per-base probability depends on the
	// homopolymer run length.
	tol := 4 * math.Sqrt(c.delV)
	if d := math.Abs(float64(c.del) - c.delP); d > tol {
		t.Errorf("deletions %d, want %.0f ± %.0f (4σ, n=%d)", c.del, c.delP, tol, c.n)
	} else {
		t.Logf("deletion rate %.5f (want %.5f, n=%d)", float64(c.del)/float64(c.n), c.delP/float64(c.n), c.n)
	}
}

func bigN(short, long int) int {
	if testing.Short() {
		return short
	}
	return long
}

func TestRatesPerTier(t *testing.T) {
	total := bigN(200_000, 1_000_000)
	for i, tier := range Tiers {
		t.Run(string(tier), func(t *testing.T) {
			m, err := New(tier)
			if err != nil {
				t.Fatal(err)
			}
			r := newRand(uint64(100 + i))
			var c counts
			var st Stats
			for got := 0; got < total; {
				truth := randomTruth(r, 1000+r.IntN(19000))
				seq, qual, ops, s := m.Apply(r, truth)
				if len(qual) != len(seq) {
					t.Fatalf("len(qual) %d != len(seq) %d", len(qual), len(seq))
				}
				tally(m, truth, ops, 0, len(truth), &c)
				st.Subs += s.Subs
				st.Ins += s.Ins
				st.Dels += s.Dels
				got += len(truth)
			}
			checkCounts(t, m, c)
			if st.Subs != c.sub || st.Ins != c.ins || st.Dels != c.del {
				t.Errorf("Stats counts %d/%d/%d disagree with ops %d/%d/%d", st.Subs, st.Ins, st.Dels, c.sub, c.ins, c.del)
			}
		})
	}
}

// TestRatesAtLongReadTail measures rates over positions 40–50 kb of 50 kb
// reads: a model with a length cap or positional drift fails here.
func TestRatesAtLongReadTail(t *testing.T) {
	reads := bigN(10, 40)
	for i, tier := range Tiers {
		t.Run(string(tier), func(t *testing.T) {
			m, _ := New(tier)
			r := newRand(uint64(200 + i))
			var c counts
			for k := 0; k < reads; k++ {
				truth := randomTruth(r, 50_000)
				_, _, ops, _ := m.Apply(r, truth)
				tally(m, truth, ops, 40_000, 50_000, &c)
			}
			checkCounts(t, m, c)
		})
	}
}

func TestReplayReproducesRead(t *testing.T) {
	r := newRand(3)
	for _, tier := range Tiers {
		m, _ := New(tier)
		for k := 0; k < 50; k++ {
			truth := randomTruth(r, r.IntN(5000))
			if k%5 == 0 {
				truth = runsOf(r, 1+r.IntN(10), 3000)
			}
			seq, qual, ops, st := m.Apply(r, truth)
			got, err := Replay(truth, ops)
			if err != nil {
				t.Fatalf("%s: Replay: %v", tier, err)
			}
			if !bytes.Equal(got, seq) {
				t.Fatalf("%s: replaying ops does not reproduce the read", tier)
			}
			if len(qual) != len(seq) || st.ReadLen != len(seq) || st.TruthLen != len(truth) {
				t.Fatalf("%s: lengths seq %d qual %d stats %d/%d truth %d", tier, len(seq), len(qual), st.ReadLen, st.TruthLen, len(truth))
			}
			if st.Matches+st.Subs+st.Dels != len(truth) || st.Matches+st.Subs+st.Ins != len(seq) {
				t.Fatalf("%s: op counts %+v inconsistent with lengths", tier, st)
			}
			for i, op := range ops {
				if op.Kind == Sub && op.Base == truthAt(truth, ops, i) {
					t.Fatalf("%s: substitution to the same base at op %d", tier, i)
				}
			}
		}
	}
}

// truthAt returns the truth base consumed by ops[i].
func truthAt(truth []byte, ops []Op, i int) byte {
	t := 0
	for _, op := range ops[:i] {
		if op.Kind != Ins {
			t++
		}
	}
	return truth[t]
}

func TestReplayRejectsBadOps(t *testing.T) {
	truth := []byte("ACGT")
	for name, ops := range map[string][]Op{
		"short": {{Kind: Match}},
		"long":  {{Kind: Match}, {Kind: Match}, {Kind: Match}, {Kind: Match}, {Kind: Del}},
		"kind":  {{Kind: 9}},
	} {
		if _, err := Replay(truth, ops); err == nil {
			t.Errorf("%s: Replay accepted invalid ops", name)
		}
	}
}

// TestHomopolymerDeletions is the positive control: deletion rates in runs of
// length 2, 4, 6 and 8 follow the table and strictly increase, while a
// sequence without homopolymers stays at the base deletion rate.
func TestHomopolymerDeletions(t *testing.T) {
	n := bigN(200_000, 1_000_000)
	for i, tier := range Tiers {
		t.Run(string(tier), func(t *testing.T) {
			m, _ := New(tier)
			r := newRand(uint64(300 + i))
			prevRate, prevTable := -1.0, -1.0
			for _, L := range []int{1, 2, 4, 6, 8} {
				var c counts
				for c.n < n {
					truth := runsOf(r, L, 10_000)
					_, _, ops, _ := m.Apply(r, truth)
					tally(m, truth, ops, 0, len(truth), &c)
				}
				want := m.DeletionRate(L)
				if L == 1 && want != m.Params().Del {
					t.Fatalf("DeletionRate(1) = %v, want base rate %v", want, m.Params().Del)
				}
				within4σ(t, "deletion in runs of "+string(rune('0'+L)), c.del, c.n, want)
				// The other error types are not homopolymer dependent.
				within4σ(t, "substitution in runs of "+string(rune('0'+L)), c.sub, c.n, m.Params().Sub)
				rate := float64(c.del) / float64(c.n)
				if L > 1 && (rate <= prevRate || want <= prevTable) {
					t.Errorf("deletion rate does not increase at run length %d: %.5f (table %.5f) after %.5f (table %.5f)", L, rate, want, prevRate, prevTable)
				}
				prevRate, prevTable = rate, want
			}
		})
	}
}

func TestHomopolymerExcessSaturates(t *testing.T) {
	for _, tier := range Tiers {
		m, _ := New(tier)
		hpCap := m.Params().HPCap
		prev := 0.0
		for L := 2; L <= 60; L++ {
			e := m.HomopolymerExcess(L)
			if e <= prev || e > hpCap {
				t.Fatalf("%s: excess(%d) = %v, previous %v, cap %v", tier, L, e, prev, hpCap)
			}
			prev = e
		}
		if prev < 0.99*hpCap {
			t.Errorf("%s: excess(60) = %v does not approach the cap %v", tier, prev, hpCap)
		}
	}
}

// TestQualities checks the ONT mean quality per tier: per-read MeanQ has its
// median near the tier mode, tiers are ordered sup > hac > fast, values are
// clamped, and error bases get lower qualities than matches.
func TestQualities(t *testing.T) {
	var medians []float64
	for i, tier := range Tiers {
		m, _ := New(tier)
		r := newRand(uint64(400 + i))
		var meanQs []float64
		var errSum, matchSum float64
		var errN, matchN int
		for k := 0; k < 400; k++ {
			truth := randomTruth(r, 2000)
			_, qual, ops, st := m.Apply(r, truth)
			if math.Abs(st.MeanQ-MeanQ(qual)) > 1e-9 {
				t.Fatalf("Stats.MeanQ %v != MeanQ(qual) %v", st.MeanQ, MeanQ(qual))
			}
			meanQs = append(meanQs, st.MeanQ)
			j := 0
			for _, op := range ops {
				if op.Kind == Del {
					continue
				}
				q := qual[j]
				if q < MinPhred || q > MaxPhred {
					t.Fatalf("quality %d outside [%d, %d]", q, MinPhred, MaxPhred)
				}
				if op.Kind == Match {
					matchSum += float64(q)
					matchN++
				} else {
					errSum += float64(q)
					errN++
				}
				j++
			}
		}
		sort.Float64s(meanQs)
		med := meanQs[len(meanQs)/2]
		mode := m.Params().ModalQ
		if math.Abs(med-mode) > 0.75 {
			t.Errorf("%s: median read MeanQ %.2f, want the mode Q%.0f ± 0.75", tier, med, mode)
		}
		if em, mm := errSum/float64(errN), matchSum/float64(matchN); em >= mm-3 {
			t.Errorf("%s: mean quality at errors %.1f not well below matches %.1f", tier, em, mm)
		}
		t.Logf("%s: median read MeanQ %.2f (mode %.0f)", tier, med, mode)
		medians = append(medians, med)
	}
	if !(medians[0] > medians[1] && medians[1] > medians[2]) {
		t.Errorf("median MeanQ not ordered sup > hac > fast: %v", medians)
	}
}

func TestQualitiesOnlyPath(t *testing.T) {
	for i, tier := range Tiers {
		m, _ := New(tier)
		r := newRand(uint64(500 + i))
		var meds []float64
		for k := 0; k < 300; k++ {
			q := m.Qualities(r, 1000)
			if len(q) != 1000 {
				t.Fatalf("len %d", len(q))
			}
			meds = append(meds, MeanQ(q))
		}
		// Baselines are calibrated for the error penalties, so without
		// errors the median read sits a little above the mode.
		sort.Float64s(meds)
		if med, mode := meds[len(meds)/2], m.Params().ModalQ; med < mode || med > mode+2.5 {
			t.Errorf("%s: median error-free read MeanQ %.2f, want in [Q%.0f, Q%.0f+2.5]", tier, med, mode, mode)
		} else {
			t.Logf("%s: median error-free read MeanQ %.2f", tier, med)
		}
	}
	m, _ := New(Sup)
	if q := m.Qualities(newRand(1), 0); len(q) != 0 {
		t.Errorf("Qualities(0) returned %d values", len(q))
	}
}

// TestApplyQualitiesMatchesApply pins the draw contract: Apply equals
// Qualities followed by ApplyQualities on the same generator, and
// ApplyQualities does not modify its input.
func TestApplyQualitiesMatchesApply(t *testing.T) {
	m, _ := New(Hac)
	truth := randomTruth(newRand(9), 30_000)
	seq1, qual1, _, _ := m.Apply(newRand(10), truth)
	r := newRand(10)
	base := m.Qualities(r, len(truth))
	saved := append([]byte(nil), base...)
	seq2, qual2 := m.ApplyQualities(r, truth, base)
	if !bytes.Equal(seq1, seq2) || !bytes.Equal(qual1, qual2) {
		t.Error("Apply != Qualities + ApplyQualities")
	}
	if !bytes.Equal(base, saved) {
		t.Error("ApplyQualities modified its input qualities")
	}
	// A mismatched quality slice falls back to drawn baselines.
	seq3, qual3 := m.ApplyQualities(newRand(10), truth, nil)
	if !bytes.Equal(seq1, seq3) || !bytes.Equal(qual1, qual3) {
		t.Error("ApplyQualities(nil) != Apply")
	}
}

func TestDeterministicAndConcurrent(t *testing.T) {
	m, _ := New(Fast)
	truth := randomTruth(newRand(11), 20_000)
	type out struct{ seq, qual []byte }
	want := make([]out, 8)
	for g := range want {
		s, q, _, _ := m.Apply(newRand(uint64(g)), truth)
		want[g] = out{s, q}
	}
	if bytes.Equal(want[0].seq, want[1].seq) {
		t.Error("different seeds gave identical reads")
	}
	got := make([]out, len(want))
	var wg sync.WaitGroup
	for g := range got {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for k := 0; k < 5; k++ { // repeat to give the race detector work
				s, q, _, _ := m.Apply(newRand(uint64(g)), truth)
				got[g] = out{s, q}
			}
		}()
	}
	wg.Wait()
	for g := range got {
		if !bytes.Equal(got[g].seq, want[g].seq) || !bytes.Equal(got[g].qual, want[g].qual) {
			t.Errorf("seed %d: concurrent output differs from sequential", g)
		}
	}
}

func TestNonACGTTruth(t *testing.T) {
	m, _ := New(Fast)
	truth := bytes.Repeat([]byte("NNacgtN"), 2000)
	seq, qual, ops, _ := m.Apply(newRand(12), truth)
	if got, err := Replay(truth, ops); err != nil || !bytes.Equal(got, seq) || len(qual) != len(seq) {
		t.Fatalf("Replay mismatch on non-ACGT truth: %v", err)
	}
	for _, op := range ops {
		if (op.Kind == Sub || op.Kind == Ins) && bytes.IndexByte([]byte("ACGTNacgt"), op.Base) < 0 {
			t.Fatalf("emitted invalid base %q", op.Base)
		}
	}
}

func TestParseTierAndParams(t *testing.T) {
	for _, tier := range Tiers {
		if got, err := ParseTier(string(tier)); err != nil || got != tier {
			t.Errorf("ParseTier(%q) = %q, %v", tier, got, err)
		}
	}
	for _, bad := range []string{"", "SUP", "ultra"} {
		if _, err := ParseTier(bad); err == nil {
			t.Errorf("ParseTier(%q) accepted", bad)
		}
	}
	if _, err := New("ultra"); err == nil {
		t.Error("New accepted an unknown tier")
	}
	ok := Sup.Params()
	for name, edit := range map[string]func(*Params){
		"negative": func(p *Params) { p.Sub = -0.1 },
		"nan":      func(p *Params) { p.Del = math.NaN() },
		"sum":      func(p *Params) { p.Sub, p.Del = 0.6, 0.6 },
		"ins":      func(p *Params) { p.Ins = 1 },
		"scale":    func(p *Params) { p.HPScale = 0 },
		"modal":    func(p *Params) { p.ModalQ = 90 },
	} {
		p := ok
		edit(&p)
		if _, err := NewParams(p); err == nil {
			t.Errorf("%s: NewParams accepted invalid parameters", name)
		}
	}
	// Tier rates are ordered sup < hac < fast.
	for i := 1; i < len(Tiers); i++ {
		a, b := Tiers[i-1].Params(), Tiers[i].Params()
		if !(a.Sub < b.Sub && a.Ins < b.Ins && a.Del < b.Del && a.HPCap < b.HPCap && a.ModalQ > b.ModalQ) {
			t.Errorf("tier %s is not uniformly more accurate than %s", Tiers[i-1], Tiers[i])
		}
	}
}

func TestMeanQ(t *testing.T) {
	if MeanQ(nil) != 0 {
		t.Error("MeanQ(nil) != 0")
	}
	if got := MeanQ([]byte{20, 20}); math.Abs(got-20) > 1e-9 {
		t.Errorf("MeanQ(20,20) = %v", got)
	}
	// One Q10 base among Q30s dominates: mean error (0.1+0.001)/2.
	if got, want := MeanQ([]byte{10, 30}), -10*math.Log10(0.101/2); math.Abs(got-want) > 1e-9 {
		t.Errorf("MeanQ(10,30) = %v, want %v", got, want)
	}
}

func BenchmarkApply50kb(b *testing.B) {
	m, _ := New(Sup)
	r := newRand(1)
	truth := randomTruth(r, 50_000)
	b.SetBytes(int64(len(truth)))
	for b.Loop() {
		m.Apply(r, truth)
	}
}
