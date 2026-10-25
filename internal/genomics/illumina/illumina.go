// Package illumina models Illumina short-read output: per-cycle binned base
// qualities, quality-matched substitution errors, indel errors and
// paired-end insert (fragment) sizes. Every function draws only from the *rand.Rand it is
// given, so output is a pure function of the caller's seed.
//
// The parameters are fitted to real data, one Profile per run (#35): the
// first 200,000 read pairs of an ENA run, Illumina NovaSeq 6000, human
// whole-genome sequencing, aligned to GRCh38 with minimap2.
// scripts/acceptance/illumina_calibrate.sh downloads the reads, measures
// them and writes the run's profile_<accession>.go; the measured tables are
// committed in scripts/acceptance/calibrate/illumina/results (DRR276626)
// and results/profiles/<accession>, and the tests compare each profile
// with its run.
//
// # Run profiles
//
// NovaSeq 6000 runs differ in their qualities more than the model's
// sequence-level effects: the quality level, the share of Q11, how often
// the quality changes between cycles and the Q37 error rate are properties
// of the run (results/heldout/calibration_stats.tsv), and the insert size
// is the library's. So the whole model (per-cycle curves, read classes and
// their mate coupling, stickiness, context factors, substitution model and
// insert sizes) is fitted to each run by the same pipeline, and selected by
// name (illumina.profile):
//
//   - novaseq6000-drr394916 (the default in 0.13.0–0.14.0): DRR394916, NCCHD, TruSeq DNA
//     PCR-Free, 2×151; mean Q 35.29 / 34.60 (read 1 / read 2), 90.3 % /
//     86.9 % at Q30 or more, Q11 4.0 % / 6.0 %, insert median 559;
//   - novaseq6000-drr276626 (the default of earlier releases): DRR276626,
//     SAGA_U, library not recorded, 2×151; mean Q 36.22 / 36.08, 95.4 % /
//     94.7 %, Q11 1.6 % / 2.0 %, insert median 446;
//   - novaseq6000-err15600747: ERR15600747, EBI / Southampton, Agilent
//     SureSelect QXT, 2×150; mean Q 34.57 / 34.31, 86.5 % / 85.5 %, Q11
//     5.8 % / 6.8 %, insert median 361;
//   - novaseq6000-pooled (the default, #36): the three runs above pooled, the same
//     200,000 pairs of each fitted together as one run, so that the read
//     classes span the runs (2×151; the 2×150 run's reads keep their
//     length); mean Q 35.36 / 34.99, 90.7 % / 89.0 %, Q11 3.8 % / 4.9 %,
//     inserts the pooled pairs' (median 462, lognormal 465 ± 140). Against
//     the NovaSeq 6000 runs no profile was fitted to it is the closest
//     profile (scripts/acceptance/calibrate/illumina/results/heldout,
//     profile_spread.tsv); against each fitted run it is further than that
//     run's own profile.
//
// The insert size belongs to the profile because a profile reproduces one
// run, and the run's library is part of it; illumina.insert_size overrides
// it. With no profile set, though, the insert size is lognormal 455 ± 100
// (DRR276626's), not the default profile's 465 ± 140: see "The default" in
// profiles.go for the exact rule. The figures in the rest of this
// documentation are DRR276626's, the run the model was first fitted to.
//
// # Qualities
//
// NovaSeq 6000 (RTA3) writes four quality values: Q2, Q11, Q25 and Q37,
// '#', ',', ':' and 'F' in Phred+33. These are the values in DRR276626 and
// in current NovaSeq 6000 output. Illumina's technical note "NovaSeq 6000
// System Quality Scores and RTA3 Software" (770-2017-010, 2017) lists
// 12, 23 and 37 (plus Q2 for no-calls), the values of the first instrument
// software; Illumina told users that NovaSeq Control Software v1.1 changed
// the bins to 2, 11, 25 and 37 without updating that note (reported in
// github.com/benjjneb/dada2 issue #791, 2019). A Q2 base is a no-call,
// written as N.
//
// A read's quality string is generated as follows.
//
//   - The read draws a quality class k (weight w_k). Reads differ a lot:
//     in DRR276626 read 1, 11 % of reads have no base below Q37, the median
//     read has 3, and 5 % have more than 30. The classes multiply the
//     per-cycle fractions of Q11 and Q25 by m_k,11 and m_k,25.
//   - The mean fractions of Q11 and Q25 at cycle c, l_11(c) and l_25(c),
//     are piecewise linear in c between knots (fitted by least squares to
//     cycles 1-150), extended linearly past the last knot, and multiplied
//     by FinalLow and FinalMid at a read's final cycle (the last cycle of a
//     NovaSeq read is about 25 % worse than the one before). So the class
//     k read has P(Q11) = m_k,11·l_11(c) and P(Q25) = m_k,25·l_25(c) at
//     cycle c (their sum capped at 0.98), P(Q2) = NoCall and Q37 otherwise.
//     The weights and multipliers were fitted by EM to the per-read counts
//     of Q11 and Q25 bases; they average to 1, so the mixture reproduces
//     l(c).
//   - Given its class a read's cycles are nearly independent: in the real
//     reads P(low | low k cycles earlier) is the same for every k from 1 to
//     80, so the clustering of low-quality calls comes from the read, not
//     from neighbouring cycles. A small stickiness (a Q11 cycle repeats
//     with probability StickLow, a Q25 cycle with StickMid, with the other
//     draws adjusted so the cycle's distribution is unchanged) matches the
//     measured P(Q11 | Q11) and P(Q25 | Q25) of adjacent cycles.
//   - Qualities depend on the read's bases (QualitiesFor, #35). In
//     DRR276626 read 1 the mean Q of T is 36.09, of A 36.21, of G 36.28 and
//     of C 36.36, and the 4-base context (two bases before, the base, one
//     after) matters more: P(Q25) is 0.047 at T followed by T and 0.019 at
//     C followed by A; after GG, read 2 has P(Q11) 0.060 at G and 0.018 at
//     C. The effect differs between read classes: in good reads low calls
//     fall at A/T-rich contexts (TTT about twice the average), in the
//     poorest at G/C-rich ones. So each class k multiplies its P(Q11) and
//     P(Q25) at a cycle by factors F_k(x) of the cycle's context x
//     (ReadClass.Context, 256 contexts), measured within the reads of each
//     class against the calls expected from the read's own counts, and
//     normalised to average 1 over the sample's context frequencies (the
//     per-cycle fractions and class levels are unchanged for reads of that
//     composition). Contexts past the read's ends or with N average over
//     the bases they could be.
//
// The context factors were measured on called bases; the engine applies
// them to the true bases of the fragment and then draws substitutions from
// the qualities, which moves some low calls to the substituted base (by up
// to 0.04 in a base's mean Q).
//
// Read 2 has its own parameters (it has more Q11 and Q25 than read 1, as on
// every Illumina instrument).
//
// The two reads of a pair come from one cluster, and their qualities go
// together: in DRR276626 the correlation of the mates' counts of bases
// below Q30 is 0.37 (below Q20, 0.28), and pairs with both mates poor (5 or
// more bases below Q20) are 3.95 times as common as if the mates were
// independent; five NovaSeq 6000 runs from other centres show 0.42-0.57
// below Q20 (results/heldout/calibration_stats.tsv). So a pair draws its
// two read classes together (PairParams): Joint[k1][k2] was fitted by EM to
// the pairs' Q11 and Q25 counts with each mate's class levels fixed, and
// raked so its margins are the per-mate class weights, so each mate's
// qualities are distributed as before. The generated pairs of
// illumina_heldout.sh have the measured coupling: 0.37 and 0.29, a
// both-poor ratio of 3.87, and a mean correlation of the mates' 10-cycle
// block mean qualities of 0.062 (measured 0.064).
// Within their classes the mates' cycles are drawn independently: a
// low-quality event at the same cluster position in both reads is not
// modelled beyond the shared class.
//
// # Errors
//
// Real Illumina bins are not exactly their Phred values. The substitution
// rate of each bin was measured as the mismatch rate against GRCh38
// (primary alignments, MAPQ ≥ 30, soft clips, indels and N excluded) less
// the rate of the sample's own variants (0.092 % per base). The variant
// rate is the Q37 mismatch rate less the Q37 error rate where mates overlap
// (a mismatch in one mate where the other mate has the reference base). Q37
// bases have one rate, 3.2·10⁻⁴ (Phred 35, not 37). Q25 and Q11 bases get
// worse along the read: e(c) = late − (late − early)·exp(−(c−1)/τ), Q11
// from 2.9 % at cycle 1 to 8.3 % at cycle 150 (Phred Q11 is 7.9 %), Q25
// from 0.07 % to 0.25 % (Phred Q25 is 0.32 %); in DRR394916 and
// ERR15600747 Q25 rises ever faster (τ < 0), from 0.13 % and 0.07 % to
// 0.72 % and 0.85 %. The rate also depends on the base: C and T are
// miscalled more often than A and G (BaseFactor), and the
// substituted base follows the measured spectrum (Matrix: C→A, G→T, T→A,
// T→G and A→C dominate; C→T and T→C are rare), both measured from Q11
// mismatches, of which about 99 % are sequencing errors, and applied to
// every bin.
//
// Indel errors (#30) are modelled too (indels.go, IndelParams): 1.2·10⁻⁵
// per base in the pooled runs (insertions 4.6·10⁻⁶, deletions 7.5·10⁻⁶;
// 0.4 % of the sequencing errors where mates overlap, which are 3.4·10⁻³
// per base), mostly slips in long homopolymers and tandem repeats. An
// earlier estimate of 7·10⁻⁵ (25 events, DRR276626) counted mates that end
// inside a long homopolymer, misaligned repeats, indels the two mates'
// alignments place differently and mismatches aligned as gaps, which the
// measurement now leaves out (scripts/acceptance/illumina_indels.sh).
//
// # Limitations
//
// Each profile's reads come from the first tile of one lane of one run and
// one library; other instruments (NovaSeq X uses other bins) and read
// lengths are not measured, and a profile is a copy of its run, not a
// model of how runs vary. Beyond 151 cycles (150 for ERR15600747) the
// curves are extrapolated. In the lower-quality runs read 2's G calls are
// poorer than the model makes them (mean Q of G 0.20 lower in DRR394916,
// 0.33 in ERR15600747), a read-level effect the classes do not carry. The context
// factors do not change with cycle (some do in the real reads: after GG on
// read 2 the excess of Q11 is largest in the first 30 cycles). A read's
// class does not depend on its sequence, but real reads that are GC-rich
// (above about 0.65), GC-poor (below 0.2) or have a long homopolymer (12
// or more) are poorer, as are artefacts such as read 2 poly-G
// (results/class_reads.tsv); so P(Q11) by context follows the measured one
// less closely than P(Q25). Making the class depend on GC and homopolymer
// length was tried: it did not measurably reduce the real-versus-generated
// accuracy (#35) and, with reads sampled from GRCh38, which lack the artefacts,
// raised read 2's mean Q by 0.07. No-calls are placed independently at the
// measured rate (real no-calls come in whole-cycle events: the 11 N of
// DRR276626's sample are all at read 1 cycle 22).
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

	"github.com/scttfrdmn/phantom-science/internal/detmath"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// NovaSeq 6000 (RTA3) quality values (Phred scores).
const (
	QNoCall = 2  // '#': no-call, written as N
	QLow    = 11 // ','
	QMid    = 25 // ':'
	QHigh   = 37 // 'F'
)

// Bins lists the quality values in increasing order.
var Bins = [4]byte{QNoCall, QLow, QMid, QHigh}

// maxLow caps P(Q11) + P(Q25) of a class at a cycle.
const maxLow = 0.98

// ReadClass is a read quality class: a share of reads and the factors by
// which their Q11 and Q25 fractions differ from the mean.
type ReadClass struct {
	Weight   float64
	Low, Mid float64
	// Context, when set, multiplies the class's Q11 and Q25 probabilities
	// at a cycle by factors of the cycle's sequence context
	// (QualitiesFor); nil means no context effect.
	Context *ContextFactors
}

// ContextFactors are the factors by which a read class's P(Q11) and P(Q25)
// at a cycle depend on the cycle's 4-base context: the bases two before,
// one before, at and one after the cycle, in sequencing orientation.
// Low[x] and Mid[x] are indexed by ContextIndex.
type ContextFactors struct {
	Low, Mid [256]float64
	// ext[e] holds the Low and Mid factors of a context with unknown
	// positions (past the read's ends, or not A/C/G/T), averaged over the
	// contexts that match its known bases, weighted by the context
	// frequencies of the measured reads (NewContextFactors); e is
	// 125·s(−2) + 25·s(−1) + 5·s(0) + s(+1) with s = 4 for unknown. For
	// fully known contexts it holds Low and Mid.
	ext [625][2]float64
}

// NewContextFactors returns the factors low and mid (indexed by
// ContextIndex) with the averages for contexts with unknown bases
// precomputed, weighted by the context frequencies freq of the reads they
// were measured on.
func NewContextFactors(freq *[256]float64, low, mid [256]float64) *ContextFactors {
	f := &ContextFactors{Low: low, Mid: mid}
	for e := range f.ext {
		s := [4]int{e / 125, e / 25 % 5, e / 5 % 5, e % 5}
		var l, m, w float64
		for x := range 256 {
			d := [4]int{x >> 6, x >> 4 & 3, x >> 2 & 3, x & 3}
			if (s[0] != 4 && s[0] != d[0]) || (s[1] != 4 && s[1] != d[1]) ||
				(s[2] != 4 && s[2] != d[2]) || (s[3] != 4 && s[3] != d[3]) {
				continue
			}
			pi := freq[x]
			l += float64(pi * low[x])
			m += float64(pi * mid[x])
			w += pi
		}
		f.ext[e] = [2]float64{1, 1}
		if w > 0 {
			f.ext[e] = [2]float64{l / w, m / w}
		}
	}
	return f
}

// ContextIndex returns the index 64·b(−2) + 16·b(−1) + 4·b(0) + b(+1) of
// the context of A/C/G/T bases b (A, C, G, T = 0-3) in ContextFactors, or
// -1 if a byte is not A, C, G or T.
func ContextIndex(prev2, prev, base, next byte) int {
	x := 0
	for _, b := range [4]byte{prev2, prev, base, next} {
		i := baseIndex[b]
		if i < 0 {
			return -1
		}
		x = x<<2 | int(i)
	}
	return x
}

// contextSymbol maps A, C, G, T (either case) to 0-3 and every other byte
// to 4 (unknown).
var contextSymbol = func() (t [256]uint8) {
	for i := range t {
		t[i] = 4
	}
	for i := range len(acgt) {
		b := acgt[i]
		t[b], t[b+'a'-'A'] = uint8(i), uint8(i) // #nosec G115 -- i < 4
	}
	return t
}()

// contextAt returns the extended context index (see ContextFactors.ext) of
// 0-based cycle c of seq.
func contextAt(seq []byte, c int) int {
	e := 0
	for i := c - 2; i <= c+1; i++ {
		s := 4
		if i >= 0 && i < len(seq) {
			s = int(contextSymbol[seq[i]])
		}
		e = e*5 + s
	}
	return e
}

// QualityParams parameterises the binned quality model of one mate (see
// the package documentation).
type QualityParams struct {
	// Knots are 1-based cycles in increasing order; Low and Mid are the
	// mean fractions of Q11 and Q25 at those cycles.
	Knots    []int
	Low, Mid []float64
	// FinalLow and FinalMid multiply the fractions at a read's last cycle.
	FinalLow, FinalMid float64
	// Classes are the read quality classes; their weights sum to 1 and
	// their weighted factors average to 1.
	Classes []ReadClass
	// StickLow and StickMid are the probabilities that a Q11 or Q25 cycle
	// is repeated at the next cycle beyond its share of that cycle.
	StickLow, StickMid float64
	// NoCall is P(Q2) per base.
	NoCall float64
}

// curve evaluates the piecewise-linear curve vals over p.Knots at 1-based
// cycle x, starting the knot search at *j (cycles must not decrease between
// calls sharing j).
func (p QualityParams) curve(vals []float64, x int, j *int) float64 {
	k := p.Knots
	if x <= k[0] || len(k) == 1 {
		return vals[0]
	}
	for *j < len(k)-2 && x > k[*j+1] {
		*j++
	}
	t := float64(x-k[*j]) / float64(k[*j+1]-k[*j])
	return max(0, vals[*j]+float64(t*(vals[*j+1]-vals[*j])))
}

// means returns the mean Q11 and Q25 fractions at 0-based cycle c of a read
// of n cycles.
func (p QualityParams) means(c, n int, j *int) (low, mid float64) {
	low, mid = p.curve(p.Low, c+1, j), p.curve(p.Mid, c+1, j)
	if c == n-1 {
		low, mid = float64(low*p.FinalLow), float64(mid*p.FinalMid)
	}
	return low, mid
}

// classTarget returns class k's distribution over Bins given the mean
// fractions.
func (p QualityParams) classTarget(k ReadClass, low, mid float64) [4]float64 {
	return p.contextTarget(k, low, mid, [2]float64{1, 1})
}

// contextTarget returns class k's distribution over Bins given the mean
// fractions and the context factors f of Q11 and Q25 (1, 1 for none).
func (p QualityParams) contextTarget(k ReadClass, low, mid float64, f [2]float64) [4]float64 {
	a, b := float64(float64(k.Low*low)*f[0]), float64(float64(k.Mid*mid)*f[1])
	if s := a + b; s > maxLow {
		a, b = a*maxLow/s, b*maxLow/s
	}
	nc := p.NoCall
	return [4]float64{nc, a, b, max(0, 1-nc-a-b)}
}

// Target returns the distribution over Bins at 0-based cycle c of a read of
// n cycles, averaged over the read classes.
func (p QualityParams) Target(c, n int) [4]float64 {
	j := 0
	low, mid := p.means(c, n, &j)
	var t [4]float64
	w := 0.0
	for _, k := range p.Classes {
		ct := p.classTarget(k, low, mid)
		for i := range t {
			t[i] += float64(k.Weight * ct[i])
		}
		w += k.Weight
	}
	for i := range t {
		t[i] /= w
	}
	return t
}

// class returns the index of the class of a uniform u.
func (p QualityParams) class(u float64) int {
	acc := 0.0
	for i, k := range p.Classes {
		acc += k.Weight
		if u < acc {
			return i
		}
	}
	return len(p.Classes) - 1
}

// PairParams couples the read classes of the two reads of a pair (#35):
// Joint[k1][k2] is the share of pairs whose read 1 is in class k1 of its
// QualityParams and read 2 in class k2 of its own. The shares sum to 1 and
// the margins are the classes' weights, so each mate's reads keep their
// class mixture; a nil Joint leaves the mates independent.
type PairParams struct {
	Joint [][]float64
}

// Coupled reports whether p couples the mates' classes.
func (p PairParams) Coupled() bool { return len(p.Joint) > 0 }

// Classes returns the read-1 and read-2 classes of a pair for a uniform u
// (cells in row-major order).
func (p PairParams) Classes(u float64) (k1, k2 int) {
	acc := 0.0
	for i, row := range p.Joint {
		for j, w := range row {
			acc += w
			if u < acc {
				return i, j
			}
		}
	}
	last := len(p.Joint) - 1
	return last, len(p.Joint[last]) - 1
}

// pick returns the bin of a uniform u under weights w (which sum to 1).
func pick(w [4]float64, u float64) byte {
	acc := 0.0
	for i, x := range w[:3] {
		acc += x
		if u < acc {
			return Bins[i]
		}
	}
	return QHigh
}

// TargetFor returns the distribution over Bins at 0-based cycle c of the
// read seq (sequencing orientation), averaged over the read classes, with
// each class's factors of the cycle's context in seq.
func (p QualityParams) TargetFor(seq []byte, c int) [4]float64 {
	j := 0
	low, mid := p.means(c, len(seq), &j)
	e := contextAt(seq, c)
	var t [4]float64
	w := 0.0
	for _, k := range p.Classes {
		f := [2]float64{1, 1}
		if k.Context != nil {
			f = k.Context.ext[e]
		}
		ct := p.contextTarget(k, low, mid, f)
		for i := range t {
			t[i] += float64(k.Weight * ct[i])
		}
		w += k.Weight
	}
	for i := range t {
		t[i] /= w
	}
	return t
}

// Qualities returns n binned Phred scores (not ASCII) for one read without
// a sequence (no context factors), drawing one value from r for the read's
// class and one per cycle.
func (p QualityParams) Qualities(r *rand.Rand, n int) []byte {
	return p.qualities(r, nil, n, -1)
}

// QualitiesFor returns binned Phred scores (not ASCII) for the bases seq of
// one read in sequencing orientation (seq[i] is read at cycle i): the
// probabilities of Q11 and Q25 at each cycle are those of the read's class
// times the class's factors of the cycle's context (ReadClass.Context). It
// draws from r as Qualities does: one value for the read's class and one
// per cycle.
func (p QualityParams) QualitiesFor(r *rand.Rand, seq []byte) []byte {
	return p.qualities(r, seq, len(seq), -1)
}

// QualitiesForClass is QualitiesFor for a read of class k (an index into
// Classes, e.g. from PairParams.Classes): it draws only the per-cycle
// values from r. A k outside Classes draws the class from r, as
// QualitiesFor does.
func (p QualityParams) QualitiesForClass(r *rand.Rand, seq []byte, k int) []byte {
	return p.qualities(r, seq, len(seq), k)
}

// qualities draws n qualities of class ki (drawn from r if ki is not a
// class index), with the context factors of seq unless seq is nil.
func (p QualityParams) qualities(r *rand.Rand, seq []byte, n, ki int) []byte {
	q := make([]byte, n)
	if n == 0 {
		return q
	}
	if ki < 0 || ki >= len(p.Classes) {
		ki = p.class(r.Float64())
	}
	k := p.Classes[ki]
	j := 0
	for c := range q {
		low, mid := p.means(c, n, &j)
		f := [2]float64{1, 1}
		if seq != nil && k.Context != nil {
			f = k.Context.ext[contextAt(seq, c)]
		}
		t := p.contextTarget(k, low, mid, f)
		// r.Float64 is inlined as an exact product by 2^-53; the conversion
		// keeps it from being fused into u - rho below (internal/fmacheck).
		u := float64(r.Float64())
		if c == 0 {
			q[c] = pick(t, u)
			continue
		}
		// Keep a sticky previous bin with probability rho; otherwise draw
		// from nu(b) = t(b)(1 − rho_b) / (1 − Σ rho t), which leaves the
		// cycle's distribution at t.
		rho := 0.0
		switch q[c-1] {
		case QLow:
			rho = p.StickLow
		case QMid:
			rho = p.StickMid
		}
		if u < rho {
			q[c] = q[c-1]
			continue
		}
		u = (u - rho) / (1 - rho)
		z := 1 - float64(p.StickLow*t[1]) - float64(p.StickMid*t[2])
		t[1] = float64(t[1] * ((1 - p.StickLow) / z))
		t[2] = float64(t[2] * ((1 - p.StickMid) / z))
		t[0] /= z
		t[3] /= z
		q[c] = pick(t, u)
	}
	return q
}

// ErrorParams is a substitution model for binned qualities (see the package
// documentation).
type ErrorParams struct {
	// High is the rate of Q37 bases.
	High float64
	// Q25 and Q11 rates at cycle c (1-based):
	// late − (late − early)·exp(−(c−1)/tau), at least 0. A negative tau
	// is a rate that grows ever faster along the read (convex); late is
	// then not a late-cycle rate but the asymptote towards early cycles.
	MidEarly, MidLate, MidTau float64
	LowEarly, LowLate, LowTau float64
	// BaseFactor multiplies the rate of a base (A, C, G, T); it averages to
	// 1 over the human genome's base composition.
	BaseFactor [4]float64
	// Matrix[x][y] is P(read base y | true base x, substitution), bases in
	// the order A, C, G, T.
	Matrix [4][4]float64
}

// Rate returns the substitution rate of a base with quality q at 0-based
// cycle c, averaged over bases. Qualities other than the three called bins
// get their Phred rate 10^(−q/10).
func (e ErrorParams) Rate(q byte, c int) float64 {
	x := float64(c)
	switch q {
	case QHigh:
		return e.High
	case QMid:
		return max(0, e.MidLate-float64((e.MidLate-e.MidEarly)*detmath.Exp(-x/e.MidTau)))
	case QLow:
		return max(0, e.LowLate-float64((e.LowLate-e.LowEarly)*detmath.Exp(-x/e.LowTau)))
	}
	return ErrorProb(q)
}

// errorProb[q] is the Phred error probability 10^(−q/10).
var errorProb = func() (t [94]float64) {
	for q := range t {
		t[q] = detmath.Pow(10, -float64(q)/10)
	}
	return t
}()

// ErrorProb returns the Phred error probability 10^(−q/10) of q. The model
// does not use it for the binned qualities, whose measured rates differ
// (see SubstitutionRate).
func ErrorProb(q byte) float64 { return errorProb[min(int(q), len(errorProb)-1)] }

// SubstitutionRate returns the substitution rate that ApplyErrors (the
// linear quality model's errors) applies, on average over bases, to a base
// of quality q at 0-based cycle c. The binned model's rates are its
// profile's, Profile.Errors.Rate.
func SubstitutionRate(q byte, c int) float64 { return linearErrors.Rate(q, c) }

// baseIndex maps A, C, G, T to 0-3 and every other byte to -1.
var baseIndex = func() (t [256]int8) {
	for i := range t {
		t[i] = -1
	}
	t['A'], t['C'], t['G'], t['T'] = 0, 1, 2, 3
	return t
}()

const acgt = "ACGT"

// Apply applies substitution errors to seq in place and returns seq and
// qual: a base with quality QNoCall (or lower) becomes N, and any other
// A/C/G/T base at index (cycle) i is replaced, with probability
// Rate(qual[i], i)·BaseFactor[base], by a base drawn from Matrix. Other
// bytes (N) are kept. Each non-no-call A/C/G/T base draws one value from r,
// a substitution one more.
func (e ErrorParams) Apply(r *rand.Rand, seq, qual []byte) ([]byte, []byte) {
	for i, b := range seq {
		if qual[i] <= QNoCall {
			seq[i] = 'N'
			continue
		}
		x := baseIndex[b]
		if x < 0 {
			continue
		}
		if r.Float64() >= e.Rate(qual[i], i)*e.BaseFactor[x] {
			continue
		}
		u, acc := r.Float64(), 0.0
		to := byte(0)
		for y, w := range e.Matrix[x] {
			if w == 0 {
				continue
			}
			to = acgt[y]
			acc += w
			if u < acc {
				break
			}
		}
		seq[i] = to
	}
	return seq, qual
}

// ApplyErrors applies the substitution model of the linear quality model
// (illumina.errors quality with quality_model linear) to seq in place,
// treating index i as sequencing cycle i (the engine passes reads in
// sequencing orientation). It is novaseq6000-drr276626's model
// (Profile.Errors.Apply), whatever the default profile, so linear output
// does not change with it; the binned model applies its profile's.
func ApplyErrors(r *rand.Rand, seq, qual []byte) ([]byte, []byte) {
	return linearErrors.Apply(r, seq, qual)
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
		s.sigma = math.Sqrt(detmath.Log(1 + sd*sd/(m*m)))
		s.mu = detmath.Log(m) - float64(s.sigma*s.sigma/2)
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
			v = detmath.Exp(s.mu + float64(s.sigma*detmath.NormFloat64(r)))
		} else {
			v = s.spec.Mean + float64(s.spec.SD*detmath.NormFloat64(r))
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
