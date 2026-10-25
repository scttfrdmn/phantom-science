package illumina

// Indel errors (#30).
//
// The rates were measured where the two mates of a proper pair overlap
// (scripts/acceptance/illumina_indels.sh; illumina_stats.py indels): both
// mates read the same molecule there, so an indel of the sample is in
// both, and an indel in one mate only, after left-normalisation, is a
// sequencing error. Left out as alignment artefacts are mates that end
// inside a repeat (they cannot show its indel), indels the other mate's
// alignment places within 10 bases or as mismatches, an insertion and a
// deletion of one mate within 10 bases (mismatches aligned as two gaps),
// and places where errors of different pairs cluster (misaligned repeats);
// reads generated without indel errors, measured the same way, then show
// 3·10⁻⁷ per base (the method's floor). The three runs the profiles are
// fitted to were measured on 1,000,000 pairs (ERR15600747) and 4,000,000
// (DRR394916 and DRR276626, whose long inserts overlap rarely), their pool
// on the first 1,000,000 pairs of each, and three runs no profile is
// fitted to on 1,000,000 each; the results, with exact Poisson 95 %
// intervals, are in scripts/acceptance/calibrate/illumina/results/indels.
//
// Indel errors are rare: 1.21·10⁻⁵ per base in the pool (95 % 1.08-1.35;
// insertions 4.6·10⁻⁶, deletions 7.5·10⁻⁶, from 312 errors), 0.4 % of the
// sequencing errors where mates overlap (substitutions are 3.4·10⁻³ per
// base). The runs differ, 0.9·10⁻⁵ to 1.6·10⁻⁵ per base (summary.tsv).
// About half are of one base, a fifth of two. They are concentrated in
// repeats: a homopolymer run of 20 bases or more loses a base in about 4 %
// of reads, a single base in 0.0005 %. So the model is a slip per repeat
// plus a background:
//
//   - homopolymer runs (one base repeated h times, h ≥ 1): an insertion of
//     the run's base or a deletion from the run with probability HPIns(h)
//     and HPDel(h) per run, log-linear in h (×1.48 per base) up to a cap;
//     the sizes (mostly one base, longer in runs of LongRun or more) from
//     HPInsSize and HPDelSize;
//   - short tandem repeats (a unit of 2-6 bases repeated c ≥ 2 times): an
//     insertion or deletion of whole units with probability STRIns(c) and
//     STRDel(c) per tract, log-linear in c (×1.7 per copy);
//   - background: an insertion of bases unlike their neighbours, or a
//     deletion, at BgIns and BgDel per base, with sizes from BgInsSize and
//     BgDelSize (background deletions are of 2 or more bases: a deletion of
//     one base is a homopolymer slip).
//
// The shape (the log-linear fits, the background rates and the sizes) is
// fitted to the pool; each profile scales its insertion and deletion rates
// to its own run's, so it has its run's rates on the run's sequence. The
// measurement does not find every error: reads generated with
// novaseq6000-err15600747 and measured the same way had 0.85 times its
// run's insertion rate and 0.94 times its deletion rate, so the rates are
// raised by that much (×1.18 and ×1.06; fit_indels.py). With the final
// model, generated reads measure insertions and deletions at 4.6·10⁻⁶ and
// 7.4·10⁻⁶ (novaseq6000-pooled; its runs 4.6·10⁻⁶ and 7.5·10⁻⁶) and
// 5.2·10⁻⁶ and 7.8·10⁻⁶ (novaseq6000-err15600747; its run 4.5·10⁻⁶ and
// 7.8·10⁻⁶), each inside the other's 95 % interval or overlapping it (the
// closed loop, TestCalibration_IndelClosedLoop).
//
// The read's qualities are drawn for the bases it reads (an inserted base
// gets a quality like any other); then the quality of each inserted base is
// redrawn from InsQual and those of the two bases either side of a
// deletion from DelQual: in the real reads 25 % of inserted bases and 47 %
// of the bases next to a deletion are Q11, against 3 % of all bases.
//
// Illumina reads are a fixed number of cycles, so an insertion makes the
// read take one template base fewer and a deletion one more: the read
// reads on into the fragment. A deletion that would leave the read short of
// bases (a fragment barely longer than the read) is not made.
//
// Not modelled, because it was not measured to matter: a dependence on
// cycle (chi-square of five cycle blocks against one rate p = 0.32 for
// insertions, 0.33 for deletions; the first 30 cycles are measured only in
// short inserts) or on mate (read 1 and read 2: 61 and 57 insertions, 95
// and 99 deletions). Not modelled either: a dependence on the run's base;
// how the rates vary between runs and tiles beyond the profiles' own (of
// three runs no profile is fitted to, two have 0.66 and 0.81 times the
// deletions the pooled model expects, outside their 95 % intervals, and
// all three 1.1-1.3 times its insertions, inside them); and indel errors at
// the same place in both mates (counted as the sample's indels by the
// measurement, and drawn independently here).

import (
	"encoding/binary"
	"math"
	"math/bits"
	"math/rand/v2"

	"github.com/scttfrdmn/phantom-science/internal/detmath"
)

// IndelParams is the indel-error model of a profile (see above). The rate
// tables are indexed by run length h − 1 (HP) or copies c − 2 (STR), the
// last entry applying to longer repeats; the size tables give P(size =
// i + 1), in bases (homopolymer and background) or units (tandem repeats).
type IndelParams struct {
	HPIns, HPDel   []float64
	STRIns, STRDel []float64
	BgIns, BgDel   float64
	// HPInsSize and HPDelSize are the sizes of slips in runs shorter than
	// LongRun ([0]) and in the others ([1]).
	HPInsSize, HPDelSize   [2][]float64
	STRInsSize, STRDelSize []float64
	BgInsSize, BgDelSize   []float64
	// InsQual and DelQual are the quality distributions (over Bins) of an
	// inserted base and of the bases either side of a deletion.
	InsQual, DelQual [4]float64
}

// LongRun is the run length from which homopolymer slips have their own
// sizes (IndelParams.HPInsSize[1], HPDelSize[1]).
const LongRun = 10

// Tandem repeats have units of strMinUnit to strMaxUnit bases.
const (
	strMinUnit = 2
	strMaxUnit = 6
)

// Runs of hpDirect bases or more, and tracts of strDirect copies or more,
// draw one value each; shorter ones, with rates below 10⁻³, are thinned
// (candidate) so that a read draws a few values, not one per run.
const (
	hpDirect  = 10
	strDirect = 5
)

// maxGap caps a geometric gap (opportunities to the next candidate).
const maxGap = 1 << 40

// Enabled reports whether p makes any indel errors.
func (p *IndelParams) Enabled() bool {
	if p.BgIns > 0 || p.BgDel > 0 {
		return true
	}
	for _, t := range [][]float64{p.HPIns, p.HPDel, p.STRIns, p.STRDel} {
		for _, x := range t {
			if x > 0 {
				return true
			}
		}
	}
	return false
}

// IndelClass is the kind of place an indel error is made.
type IndelClass uint8

// Indel classes.
const (
	IndelHomopolymer IndelClass = iota + 1
	IndelTandem
	IndelBackground
)

// Indel is one indel error of a read.
type Indel struct {
	// Cycle is the 0-based cycle of the first inserted base, or of the
	// base read after a deletion.
	Cycle int
	// Pos is the 0-based template index of the error: an insertion comes
	// before template base Pos, a deletion removes template bases
	// [Pos, Pos+Len).
	Pos int
	// Len is the inserted or deleted bases.
	Len       int
	Insertion bool
	Class     IndelClass
	// Repeat is the run length (homopolymer) or copies (tandem repeat) of
	// the template where the error is made; 0 for background errors.
	Repeat int
}

// IndelModel is IndelParams prepared for drawing.
type IndelModel struct {
	p IndelParams
	// q is the candidate rate of the thinned short runs (hp), short tracts
	// (str) and background bases (bg); logq is log(1 − q).
	hpQ, strQ, bgQ          float64
	hpLogQ, strLogQ, bgLogQ float64
}

// NewIndelModel prepares p for drawing, or returns nil when p makes no
// errors.
func NewIndelModel(p IndelParams) *IndelModel {
	if !p.Enabled() {
		return nil
	}
	m := &IndelModel{p: p}
	for h := 1; h < hpDirect; h++ {
		m.hpQ = max(m.hpQ, rateAt(p.HPIns, h-1)+rateAt(p.HPDel, h-1))
	}
	for c := 2; c < strDirect; c++ {
		m.strQ = max(m.strQ, rateAt(p.STRIns, c-2)+rateAt(p.STRDel, c-2))
	}
	m.bgQ = p.BgIns + p.BgDel
	m.hpLogQ, m.strLogQ, m.bgLogQ = logComplement(m.hpQ), logComplement(m.strQ), logComplement(m.bgQ)
	return m
}

// Params returns the model's parameters.
func (m *IndelModel) Params() IndelParams { return m.p }

// logComplement returns log(1 − q), or 0 for q = 0.
func logComplement(q float64) float64 {
	if q <= 0 {
		return 0
	}
	return detmath.Log(1 - min(q, 0.5))
}

// rateAt returns t[i], or t's last entry for i past its end (0 for an
// empty t).
func rateAt(t []float64, i int) float64 {
	if len(t) == 0 {
		return 0
	}
	return t[min(i, len(t)-1)]
}

// gap draws the number of opportunities before the next candidate of a
// process with candidate rate q per opportunity (logq = log(1 − q)).
func gap(r *rand.Rand, logq float64) int {
	// r.Float64 is inlined as an exact product by 2^-53; the conversion
	// keeps it from being fused into the subtraction (internal/fmacheck).
	u := 1 - float64(r.Float64()) // (0, 1]
	if logq >= 0 {
		return maxGap
	}
	g := detmath.Log(u) / logq
	if !(g < maxGap) {
		return maxGap
	}
	return int(math.Floor(g))
}

// draw returns the 1-based index drawn from the distribution w (which sums
// to 1) for a uniform u.
func draw(w []float64, u float64) int {
	acc := 0.0
	for i, x := range w {
		acc += x
		if u < acc {
			return i + 1
		}
	}
	return max(len(w), 1)
}

// isACGT reports whether b is A, C, G or T.
func isACGT(b byte) bool { return baseIndex[b] >= 0 }

// tractUnit reports whether the k bases at t[i:] are a primitive unit (not
// a repeat of a shorter one) of A, C, G and T that the next k bases repeat.
func tractUnit(t []byte, i, k int) bool {
	if i+2*k > len(t) {
		return false
	}
	for j := range k {
		if t[i+j] != t[i+k+j] || !isACGT(t[i+j]) {
			return false
		}
	}
	for d := 1; d < k; d++ {
		if k%d != 0 {
			continue
		}
		periodic := true
		for j := d; j < k && periodic; j++ {
			periodic = t[i+j] == t[i+j-d]
		}
		if periodic {
			return false
		}
	}
	return true
}

// copiesAt returns the copies of the k-base unit at t[i:] from i.
func copiesAt(t []byte, i, k int) int {
	c := 1
	for i+(c+1)*k <= len(t) {
		for j := range k {
			if t[i+c*k+j] != t[i+j] {
				return c
			}
		}
		c++
	}
	return c
}

// tractStart reports whether a tandem repeat tract of a primitive k-base
// unit starts at t[i]: the unit is repeated at least twice from i, and the
// period-k repetition does not reach back to i − 1 (so a tract starts only
// at its left end).
func tractStart(t []byte, i, k int) bool {
	if i+2*k > len(t) || t[i] != t[i+k] || (i > 0 && t[i-1] == t[i-1+k]) {
		return false
	}
	return tractUnit(t, i, k)
}

// longTractUnits returns a bit set of the unit lengths k (bit k) for which
// a tract of strDirect or more copies may start at t[i]: t[i] is repeated
// at i + k, i + 2k, i + 3k and i + 4k (strDirect is 5). It is a quick
// filter (the bits are necessary, not sufficient) with constant offsets
// the compiler checks once.
func longTractUnits(t []byte, i int) uint8 {
	if i+5*strMaxUnit > len(t) {
		// Near the end of the template: the general test.
		var m uint8
		for k := strMinUnit; k <= strMaxUnit; k++ {
			if i+strDirect*k <= len(t) && t[i] == t[i+k] && t[i] == t[i+2*k] && t[i] == t[i+3*k] && t[i] == t[i+4*k] {
				m |= 1 << k
			}
		}
		return m
	}
	u := t[i : i+5*strMaxUnit : i+5*strMaxUnit]
	c := u[0]
	var m uint8
	if c == u[2] && c == u[4] && c == u[6] && c == u[8] {
		m |= 1 << 2
	}
	if c == u[3] && c == u[6] && c == u[9] && c == u[12] {
		m |= 1 << 3
	}
	if c == u[4] && c == u[8] && c == u[12] && c == u[16] {
		m |= 1 << 4
	}
	if c == u[5] && c == u[10] && c == u[15] && c == u[20] {
		m |= 1 << 5
	}
	if c == u[6] && c == u[12] && c == u[18] && c == u[24] {
		m |= 1 << 6
	}
	return m
}

// nextLong returns the first index in [i, end) where a homopolymer run of
// hpDirect or more bases or a tandem repeat of strDirect or more copies may
// start (a quick test: necessary, not sufficient), or end. It tests eight
// positions at a time: a long run at x has t[x] = t[x+1] = t[x+5] =
// t[x+hpDirect-1], a long tract of a k-base unit t[x] = t[x+k] = t[x+2k] =
// t[x+3k] = t[x+4k] (strDirect is 5).
func nextLong(t []byte, i, end int) int {
	x := i
	const span = 8 + 4*strMaxUnit // bytes read past x by one step
	le := binary.LittleEndian
	for ; x < end && x+span <= len(t); x += 8 {
		w := le.Uint64(t[x:])
		hit := zeroBytes((w ^ le.Uint64(t[x+1:])) | (w ^ le.Uint64(t[x+5:])) | (w ^ le.Uint64(t[x+hpDirect-1:])))
		for k := strMinUnit; k <= strMaxUnit; k++ {
			hit |= zeroBytes((w ^ le.Uint64(t[x+k:])) | (w ^ le.Uint64(t[x+2*k:])) | (w ^ le.Uint64(t[x+3*k:])) | (w ^ le.Uint64(t[x+4*k:])))
		}
		for hit != 0 {
			j := x + bits.TrailingZeros64(hit)/8
			if j >= end {
				return end
			}
			if longAt(t, j) {
				return j
			}
			hit &= hit - 1 // the flag bit of byte j (one per byte)
		}
	}
	for ; x < end; x++ {
		if longAt(t, x) {
			return x
		}
	}
	return min(x, end)
}

// zeroBytes returns a mask with the high bit of each zero byte of v set
// (and possibly of a byte above a zero byte: a necessary test).
func zeroBytes(v uint64) uint64 {
	const lo, hi = 0x0101010101010101, 0x8080808080808080
	return (v - lo) &^ v & hi
}

// longAt is nextLong's test at one position.
func longAt(t []byte, x int) bool {
	c := t[x]
	if x+hpDirect <= len(t) && t[x+1] == c && t[x+hpDirect-1] == c && (x == 0 || t[x-1] != c) {
		return true
	}
	return longTractUnits(t, x) != 0
}

// walk is the state of one read's indel draw.
type walk struct {
	m     *IndelModel
	r     *rand.Rand
	t     []byte
	more  func(k int) []byte
	owned bool // t is a copy that may be appended to
	n     int  // read bases; -1: the whole template
	out   []byte
	edits []Indel
}

// need returns how many more read bases are wanted (-1: the whole
// template).
func (w *walk) need() int {
	if w.n < 0 {
		return -1
	}
	return w.n - len(w.out)
}

// extend asks more for at least k template bases beyond len(t) and reports
// whether t now holds them.
func (w *walk) extend(k int) bool {
	if w.more == nil || k <= 0 {
		return k <= 0
	}
	ext := w.more(k)
	if len(ext) == 0 {
		return false
	}
	if !w.owned {
		w.t = append(make([]byte, 0, len(w.t)+len(ext)+16), w.t...)
		w.owned = true
	}
	w.t = append(w.t, ext...)
	return len(ext) >= k
}

// canDelete reports whether k template bases from i can be deleted with
// enough left (extending the template if needed) to fill the read.
func (w *walk) canDelete(i, k int) bool {
	if w.n < 0 {
		return i+k <= len(w.t)
	}
	short := w.need() - (len(w.t) - i - k)
	return short <= 0 || w.extend(short)
}

// insert writes the inserted bases b (at cycle len(out); as many as the
// read has room for) and records the error.
func (w *walk) insert(i int, b []byte, c IndelClass, rep int) {
	if n := w.need(); n >= 0 && len(b) > n {
		b = b[:n]
	}
	if len(b) == 0 {
		return
	}
	w.edits = append(w.edits, Indel{Cycle: len(w.out), Pos: i, Len: len(b), Insertion: true, Class: c, Repeat: rep})
	w.out = append(w.out, b...)
}

// remove records a deletion of k template bases at i (the caller skips
// them).
func (w *walk) remove(i, k int, c IndelClass, rep int) {
	w.edits = append(w.edits, Indel{Cycle: len(w.out), Pos: i, Len: k, Class: c, Repeat: rep})
}

// decide returns +1 (an insertion), -1 (a deletion) or 0 for a uniform u
// scaled to the opportunity's probabilities pi and pd.
func decide(u, pi, pd float64) int {
	switch {
	case u < pi:
		return 1
	case u < pi+pd:
		return -1
	}
	return 0
}

// candidate is a thinned process over template positions: each position is
// a candidate with probability q, independently (geometric gaps), and a
// candidate that holds an opportunity of probabilities pi, pd is accepted
// as an insertion with probability pi/q and a deletion with pd/q, so each
// opportunity has its own probabilities. next is the next candidate
// position.
type candidate struct {
	q, logq float64
	next    int
}

// refresh replaces a candidate passed over (its position deleted, or a
// position the process was not asked about) by a draw from i.
func (c *candidate) refresh(r *rand.Rand, i int) {
	if c.next < i {
		c.next = i + gap(r, c.logq)
	}
}

// at reports whether position i is a candidate, drawing the next one.
func (c *candidate) at(r *rand.Rand, i int) bool {
	c.refresh(r, i)
	if c.next != i {
		return false
	}
	c.next = i + 1 + gap(r, c.logq)
	return true
}

// Read returns the bases of a read of template t (the bases in sequencing
// orientation, t[0] first) with indel errors, and its errors in cycle
// order. With n ≥ 0 the read is n bases: it reads on into t past n bases
// after deletions, and asks more (nil: no more) for at least k bases beyond
// t when t runs out; a deletion that would leave the read short is not
// made. With n < 0 the read is the whole of t with the errors (its length
// changes by them). t is not modified.
//
// Homopolymer runs of hpDirect bases or more and tandem repeats of
// strDirect copies or more draw one value each; the rarer opportunities
// (shorter runs and tracts, background bases) are thinned (candidate), so a
// read draws three values and a few more per error.
func (m *IndelModel) Read(r *rand.Rand, t []byte, n int, more func(k int) []byte) ([]byte, []Indel) {
	if m == nil {
		if n < 0 || n > len(t) {
			n = len(t)
		}
		return append([]byte(nil), t[:n]...), nil
	}
	size := n
	if size < 0 {
		size = len(t)
	}
	w := walk{m: m, r: r, t: t, more: more, n: n, out: make([]byte, 0, size)}
	hp := candidate{q: m.hpQ, logq: m.hpLogQ, next: gap(r, m.hpLogQ)}
	str := candidate{q: m.strQ, logq: m.strLogQ, next: gap(r, m.strLogQ)}
	bg := candidate{q: m.bgQ, logq: m.bgLogQ, next: gap(r, m.bgLogQ)}
	p := &m.p
	i := 0
loop:
	for w.need() != 0 {
		if i >= len(w.t) {
			if w.n < 0 || !w.extend(w.need()) {
				break
			}
		}
		// Copy the template up to the next position where an error can
		// be made: a candidate of a thinned process, or where a long run
		// or tract may start.
		hp.refresh(w.r, i)
		str.refresh(w.r, i)
		bg.refresh(w.r, i)
		end := len(w.t)
		if w.n >= 0 {
			end = min(end, i+w.need())
		}
		if j := nextLong(w.t, i, min(end, hp.next, str.next, bg.next)); j > i {
			w.out = append(w.out, w.t[i:j]...)
			i = j
			continue
		}
		b := w.t[i]
		// A homopolymer run starts at i: long runs draw, short ones are
		// thinned.
		if isACGT(b) && (i == 0 || w.t[i-1] != b) {
			h := 1
			for i+h < len(w.t) && w.t[i+h] == b {
				h++
			}
			ev := 0
			pi, pd := rateAt(p.HPIns, h-1), rateAt(p.HPDel, h-1)
			if h >= hpDirect {
				ev = decide(w.r.Float64(), pi, pd)
			} else if hp.next <= i && hp.at(w.r, i) {
				ev = decide(float64(w.r.Float64()*hp.q), pi, pd)
			}
			long := 0
			if h >= LongRun {
				long = 1
			}
			switch ev {
			case 1:
				k := draw(p.HPInsSize[long], w.r.Float64())
				ins := make([]byte, k)
				for j := range ins {
					ins[j] = b
				}
				w.insert(i, ins, IndelHomopolymer, h)
			case -1:
				k := min(draw(p.HPDelSize[long], w.r.Float64()), h)
				if w.canDelete(i, k) {
					w.remove(i, k, IndelHomopolymer, h)
					i += k
					continue loop
				}
			}
		} else if hp.next <= i {
			hp.at(w.r, i) // a candidate position without a run start
		}
		if w.need() == 0 {
			break
		}
		// Tandem repeat tracts that start at i: long ones draw, short
		// ones are thinned (one candidate draw per position). Outside a
		// candidate only tracts of strDirect copies or more count, so a
		// base that the next copies do not repeat is passed over at once.
		isCand := str.next <= i && str.at(w.r, i)
		units := uint8(0xff) // unit lengths that may start a tract here
		if !isCand {
			units = longTractUnits(w.t, i)
		}
		for k := strMinUnit; units != 0 && k <= strMaxUnit; k++ {
			t := w.t
			if units&(1<<k) == 0 || !tractStart(t, i, k) {
				continue
			}
			c := copiesAt(t, i, k)
			if c < strDirect && !isCand {
				continue
			}
			ev := 0
			pi, pd := rateAt(p.STRIns, c-2), rateAt(p.STRDel, c-2)
			switch {
			case c >= strDirect:
				ev = decide(w.r.Float64(), pi, pd)
			case isCand:
				isCand = false // one short tract per candidate
				ev = decide(float64(w.r.Float64()*str.q), pi, pd)
			}
			switch ev {
			case 1:
				units := draw(p.STRInsSize, w.r.Float64())
				ins := make([]byte, 0, units*k)
				for range units {
					ins = append(ins, w.t[i:i+k]...)
				}
				w.insert(i, ins, IndelTandem, c)
			case -1:
				units := min(draw(p.STRDelSize, w.r.Float64()), c)
				if w.canDelete(i, units*k) {
					w.remove(i, units*k, IndelTandem, c)
					i += units * k
					continue loop
				}
			}
			if w.need() == 0 {
				break loop
			}
		}
		// Background.
		if bg.next <= i && bg.at(w.r, i) {
			switch decide(float64(w.r.Float64()*bg.q), p.BgIns, p.BgDel) {
			case 1:
				k := draw(p.BgInsSize, w.r.Float64())
				ins := make([]byte, k)
				prev := byte(0)
				if len(w.out) > 0 {
					prev = w.out[len(w.out)-1]
				}
				for j := range ins {
					// Unlike the base before; the last also unlike the next.
					var opts [4]byte
					no := 0
					for x := range 4 {
						c := acgt[x]
						if c != prev && (j < k-1 || c != b) {
							opts[no] = c
							no++
						}
					}
					ins[j] = opts[w.r.IntN(no)]
					prev = ins[j]
				}
				w.insert(i, ins, IndelBackground, 0)
			case -1:
				k := draw(p.BgDelSize, w.r.Float64())
				if w.canDelete(i, k) {
					w.remove(i, k, IndelBackground, 0)
					i += k
					continue loop
				}
			}
			if w.need() == 0 {
				break
			}
		}
		w.out = append(w.out, b)
		i++
	}
	return w.out, w.edits
}

// Qualities redraws the qualities qual (one read's, by cycle) at its indel
// errors edits: each inserted base's from InsQual and the bases either side
// of a deletion from DelQual, one value from r each.
func (m *IndelModel) Qualities(r *rand.Rand, qual []byte, edits []Indel) {
	for _, e := range edits {
		if e.Insertion {
			for c := e.Cycle; c < e.Cycle+e.Len && c < len(qual); c++ {
				qual[c] = pick(m.p.InsQual, r.Float64())
			}
			continue
		}
		for _, c := range [2]int{e.Cycle - 1, e.Cycle} {
			if c >= 0 && c < len(qual) {
				qual[c] = pick(m.p.DelQual, r.Float64())
			}
		}
	}
}

// ReadOffset returns the 0-based offset in the read of template index x,
// given the read's errors edits: x plus the bases inserted before it, less
// those deleted before it (a deleted base maps to where it would have
// been).
func ReadOffset(x int, edits []Indel) int {
	off := x
	for _, e := range edits {
		switch {
		case e.Insertion && e.Pos <= x:
			off += e.Len
		case !e.Insertion && e.Pos < x:
			off -= min(e.Len, x-e.Pos)
		}
	}
	return off
}

// TemplateUsed returns how many template bases a read of n bases with the
// errors edits reads: n less the inserted bases, plus the deleted ones.
func TemplateUsed(n int, edits []Indel) int {
	for _, e := range edits {
		if e.Insertion {
			n -= e.Len
		} else {
			n += e.Len
		}
	}
	return n
}
