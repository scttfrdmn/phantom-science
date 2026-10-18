// Package errmodel is a seedable Oxford Nanopore (R10.4.1) sequencing-error
// and base-quality model for synthetic reads (#11).
//
// A Model turns a truth sequence into an observed read in one left-to-right
// pass. For every truth base it draws one of
//
//	deletion      with probability Del + HomopolymerExcess(run length)
//	substitution  with probability Sub
//	match         otherwise
//
// and then, independently, an insertion after that base with probability
// Ins. Every rate is therefore a per-truth-base probability, there is no
// positional drift after indels (the pass walks the truth, not the read) and
// no length cap (50 kb+ reads behave like 1 kb reads).
//
// Qualities are Phred values (not ASCII). Each read draws a read-level
// quality around the tier's mode; every base then gets Gaussian jitter.
// Bases emitted by a substitution or insertion get roughly half the baseline
// quality, and the bases on either side of any error lose NeighbourDrop, the
// way basecallers report low confidence around errors. The per-read target is
// the ONT read quality −10·log10(mean(10^(−q/10))), not the arithmetic mean
// of q: baselines are calibrated so that a read's ONT mean Q after the error
// penalties sits at its drawn read quality, and the per-read MeanQ
// distribution peaks at the tier mode.
//
// A Model is immutable after construction and safe for concurrent use; all
// randomness comes from the *rand.Rand passed to each call.
package errmodel

import (
	"errors"
	"fmt"
	"math"
	"math/rand/v2"
)

// Tier is an ONT basecaller model tier.
type Tier string

// Basecaller tiers.
const (
	Sup  Tier = "sup"
	Hac  Tier = "hac"
	Fast Tier = "fast"
)

// Tiers lists the supported tiers from most to least accurate.
var Tiers = []Tier{Sup, Hac, Fast}

// Phred bounds of generated qualities. Dorado reports up to about Q50.
const (
	MinPhred = 1
	MaxPhred = 50
)

// Params are the parameters of a Model. Rates are per truth base.
type Params struct {
	Sub float64 // substitution probability
	Ins float64 // insertion probability (one inserted base after the truth base)
	Del float64 // deletion probability outside homopolymers

	// Homopolymer deletion excess for a base in a run of length L >= 2:
	// HPCap · (1 − exp(−(L−1)/HPScale)). It is zero for L = 1, strictly
	// increasing in L and saturates at HPCap, the observed ONT behaviour
	// (run lengths are under-called increasingly with length, Delahaye &
	// Nicolas 2021).
	HPCap   float64
	HPScale float64

	ModalQ        float64 // mode of the per-read ONT mean quality
	ReadSD        float64 // SD of the per-read quality around ModalQ
	BaseSD        float64 // SD of the per-base jitter around the read quality
	ErrFactor     float64 // quality of substituted/inserted bases = ErrFactor × baseline
	NeighbourDrop float64 // quality lost by the bases next to an error
}

// The tier table.
//
// The SUP rates are the project roadmap's (0.8% substitution, 0.8%
// insertion, 1.0% deletion: deletions dominate ONT errors, see Delahaye &
// Nicolas, "Sequencing DNA with nanopores: Troubles and biases", PLoS One
// 16:e0257521, 2021; Sereika et al., Nat Methods 19:823–826, 2022, for
// R10.4 Q20+ chemistry). HAC and FAST scale these up, keeping the
// deletion-dominant split, so the total base error is 2.6% / 5.0% / 15.0%,
// i.e. per-base accuracy ≈ 0.974 / 0.95 / 0.85 before homopolymer excess.
// HAC and FAST match the accuracies of configs/sequencers/ont_promethion24.yaml;
// SUP is the roadmap value (the YAML's 0.99 is the modal read Q20).
//
// Modal read qualities Q20 / Q15 / Q10 are the PromethION YAML's
// basecaller_modes and the commonly quoted R10.4.1 "Q20+" SUP, HAC and FAST
// levels. Current Dorado v5 models report higher modal accuracies (ONT model
// release notes quote ≈Q25 SUP); these values are deliberately conservative.
// ONT qualities are only loosely calibrated, so the reported mode and the
// simulated error rates need not satisfy Q = −10·log10(error).
//
// Homopolymer excess caps at 4% / 6% / 10% per base: R10.4.1's dual reader
// head reduced but did not remove length-dependent homopolymer deletions,
// which grow with run length and remain a leading residual error type. The
// cap and scale are modelling choices, not fitted values.
var tierParams = map[Tier]Params{
	Sup:  {Sub: 0.008, Ins: 0.008, Del: 0.010, HPCap: 0.04, HPScale: 4, ModalQ: 20, ReadSD: 2.5, BaseSD: 4, ErrFactor: 0.5, NeighbourDrop: 3},
	Hac:  {Sub: 0.016, Ins: 0.014, Del: 0.020, HPCap: 0.06, HPScale: 4, ModalQ: 15, ReadSD: 2.5, BaseSD: 4, ErrFactor: 0.5, NeighbourDrop: 3},
	Fast: {Sub: 0.050, Ins: 0.035, Del: 0.065, HPCap: 0.10, HPScale: 4, ModalQ: 10, ReadSD: 2, BaseSD: 3, ErrFactor: 0.5, NeighbourDrop: 2},
}

// ParseTier parses a tier name.
func ParseTier(s string) (Tier, error) {
	t := Tier(s)
	if _, ok := tierParams[t]; !ok {
		return "", fmt.Errorf("unknown basecaller tier %q: want sup, hac or fast", s)
	}
	return t, nil
}

// Params returns the tier's parameters (zero for an unknown tier).
func (t Tier) Params() Params { return tierParams[t] }

// Model applies ONT errors and qualities. Create one with New or NewParams.
type Model struct {
	p Params
	// meanHP is the mean homopolymer deletion excess per base of a uniform
	// random sequence, used to calibrate qualities (see centre).
	meanHP float64
}

// New returns the model of tier t.
func New(t Tier) (*Model, error) {
	p, ok := tierParams[t]
	if !ok {
		return nil, fmt.Errorf("unknown basecaller tier %q: want sup, hac or fast", string(t))
	}
	return NewParams(p)
}

// NewParams returns a model with custom parameters.
func NewParams(p Params) (*Model, error) {
	if err := p.validate(); err != nil {
		return nil, err
	}
	m := &Model{p: p}
	// In i.i.d. uniform sequence a base lies in a run of length L with
	// probability L·(3/4)²·(1/4)^(L−1).
	for L := 2; L <= 64; L++ {
		m.meanHP += float64(L) * 9 / 16 * math.Pow(0.25, float64(L-1)) * m.HomopolymerExcess(L)
	}
	return m, nil
}

func (p Params) validate() error {
	for _, v := range []float64{p.Sub, p.Ins, p.Del, p.HPCap, p.ReadSD, p.BaseSD, p.ErrFactor, p.NeighbourDrop} {
		if !(v >= 0) || math.IsInf(v, 0) {
			return errors.New("errmodel: parameters must be finite and non-negative")
		}
	}
	switch {
	case p.Ins >= 1:
		return errors.New("errmodel: insertion rate must be < 1")
	case p.Sub+p.Del+p.HPCap > 1:
		return errors.New("errmodel: substitution + deletion + homopolymer cap must be <= 1")
	case p.HPCap > 0 && !(p.HPScale > 0):
		return errors.New("errmodel: homopolymer scale must be positive")
	case !(p.ModalQ >= MinPhred && p.ModalQ <= MaxPhred): // also rejects NaN
		return fmt.Errorf("errmodel: modal quality must be in [%d, %d]", MinPhred, MaxPhred)
	case p.ErrFactor > 1:
		return errors.New("errmodel: error quality factor must be <= 1")
	}
	return nil
}

// Params returns the model's parameters.
func (m *Model) Params() Params { return m.p }

// HomopolymerExcess is the extra deletion probability of a base in a
// homopolymer run of length runLen (0 for runLen <= 1).
func (m *Model) HomopolymerExcess(runLen int) float64 {
	if runLen <= 1 || m.p.HPCap == 0 {
		return 0
	}
	return m.p.HPCap * -math.Expm1(-float64(runLen-1)/m.p.HPScale)
}

// DeletionRate is the deletion probability of a base in a run of length
// runLen (1 outside homopolymers).
func (m *Model) DeletionRate(runLen int) float64 {
	return m.p.Del + m.HomopolymerExcess(runLen)
}

// OpKind is the kind of an edit operation.
type OpKind uint8

// Edit operations.
const (
	Match OpKind = iota // truth base emitted unchanged
	Sub                 // truth base replaced by Op.Base
	Ins                 // Op.Base emitted; consumes no truth base
	Del                 // truth base dropped
)

func (k OpKind) String() string {
	switch k {
	case Match:
		return "match"
	case Sub:
		return "sub"
	case Ins:
		return "ins"
	case Del:
		return "del"
	}
	return fmt.Sprintf("OpKind(%d)", uint8(k))
}

// Op is one edit operation. Base is the emitted base for Sub and Ins.
type Op struct {
	Kind OpKind
	Base byte
}

// Stats summarises one Apply call.
type Stats struct {
	TruthLen, ReadLen        int
	Matches, Subs, Ins, Dels int
	// MeanQ is the ONT read quality −10·log10(mean(10^(−q/10))) of the
	// output qualities (0 for an empty read).
	MeanQ float64
}

// Replay applies ops to truth and returns the observed read.
func Replay(truth []byte, ops []Op) ([]byte, error) {
	out := make([]byte, 0, len(truth))
	t := 0
	for i, op := range ops {
		switch op.Kind {
		case Ins:
			out = append(out, op.Base)
			continue
		case Match, Sub, Del:
		default:
			return nil, fmt.Errorf("op %d: invalid kind %v", i, op.Kind)
		}
		if t >= len(truth) {
			return nil, fmt.Errorf("op %d (%v) past end of truth (%d bases)", i, op.Kind, len(truth))
		}
		switch op.Kind {
		case Match:
			out = append(out, truth[t])
		case Sub:
			out = append(out, op.Base)
		}
		t++
	}
	if t != len(truth) {
		return nil, fmt.Errorf("ops consume %d of %d truth bases", t, len(truth))
	}
	return out, nil
}

// MeanQ returns the ONT mean quality −10·log10(mean(10^(−q/10))) of Phred
// values q (0 for none).
func MeanQ(qual []byte) float64 {
	if len(qual) == 0 {
		return 0
	}
	var s float64
	for _, q := range qual {
		s += math.Pow(10, -float64(q)/10)
	}
	return -10 * math.Log10(s/float64(len(qual)))
}

// ReadQuality draws one read's target quality (ONT mean) around the mode.
func (m *Model) ReadQuality(r *rand.Rand) float64 {
	q := m.p.ModalQ + m.p.ReadSD*r.NormFloat64()
	return math.Min(MaxPhred, math.Max(MinPhred, q))
}

// Qualities returns n baseline Phred qualities for one read: a read quality
// from ReadQuality, then per-base Gaussian jitter around the calibrated
// centre (see centre). It draws one normal for the read and one per base.
// This is the quality-only path; ApplyQualities penalises error positions.
func (m *Model) Qualities(r *rand.Rand, n int) []byte {
	q := make([]byte, n)
	centre := m.centre(m.ReadQuality(r))
	for i := range q {
		q[i] = clampQ(centre + m.p.BaseSD*r.NormFloat64())
	}
	return q
}

// centre returns the arithmetic centre c of the per-base jitter for a read
// of ONT mean quality q. For baseline b ~ N(c, σ²) and a = ln10/10,
// E[10^(−b/10)] = exp(−a·c + a²σ²/2). Apply then lowers error bases to
// ErrFactor·b and neighbouring bases by NeighbourDrop, so c is solved (by
// bisection) for the expected mean error probability of the read after
// those penalties, at the tier's rates on random sequence, to equal
// 10^(−q/10). The read's ONT mean quality after Apply is then ≈ q; without
// errors it is a little above q. Clamping is ignored.
func (m *Model) centre(q float64) float64 {
	p := m.p
	a := math.Ln10 / 10
	v := a * a * p.BaseSD * p.BaseSD / 2
	del := p.Del + m.meanHP
	errEmit := p.Sub + p.Ins                     // emitted error bases per truth base
	nb := math.Min(2*(del+errEmit), 1-del-p.Sub) // bases bordering an error
	match := 1 - del - p.Sub - nb                // unpenalised bases
	emitted := 1 - del + p.Ins
	target := math.Exp(-a * q)
	mean := func(c float64) float64 {
		b := math.Exp(-a*c + v)
		e := math.Exp(-a*p.ErrFactor*c + v*p.ErrFactor*p.ErrFactor)
		return (match*b + nb*b*math.Exp(a*p.NeighbourDrop) + errEmit*e) / emitted
	}
	lo, hi := q-10, q+40 // mean is decreasing in c
	for range 50 {
		mid := (lo + hi) / 2
		if mean(mid) > target {
			lo = mid
		} else {
			hi = mid
		}
	}
	return (lo + hi) / 2
}

func clampQ(v float64) byte {
	v = math.Round(v)
	if !(v >= MinPhred) { // also catches NaN
		return MinPhred
	}
	if v > MaxPhred {
		return MaxPhred
	}
	return byte(v)
}

// Apply returns the observed read for truth: bases, Phred qualities (same
// length), the edit operations (Replay(truth, ops) == seq) and statistics.
// Baseline qualities come from Qualities; error positions are then
// penalised.
func (m *Model) Apply(r *rand.Rand, truth []byte) (seq, qual []byte, ops []Op, st Stats) {
	base := m.Qualities(r, len(truth))
	ops = make([]Op, 0, len(truth)+len(truth)/16)
	seq, qual = m.apply(r, truth, base, &ops, &st)
	return seq, qual, ops, st
}

// ApplyQualities applies errors to truth keeping the given baseline
// qualities (one per truth base, e.g. from Qualities or an HMM) and
// penalising error positions. It draws exactly what Apply draws after its
// baseline qualities. qual is not modified; if its length differs from
// truth's, baselines are drawn from Qualities instead.
func (m *Model) ApplyQualities(r *rand.Rand, truth, qual []byte) ([]byte, []byte) {
	if len(qual) != len(truth) {
		qual = m.Qualities(r, len(truth))
	}
	return m.apply(r, truth, qual, nil, nil)
}

// apply is the single pass. ops and st may be nil.
func (m *Model) apply(r *rand.Rand, truth, base []byte, ops *[]Op, st *Stats) ([]byte, []byte) {
	const acgt = "ACGT"
	p := m.p
	n := len(truth)
	seq := make([]byte, 0, n+n/16+1)
	qual := make([]byte, 0, n+n/16+1)
	errQ := func(b byte) byte { return clampQ(p.ErrFactor * float64(b)) }
	// Each emitted base is penalised at most once: error bases get errQ,
	// bases bordering an error lose NeighbourDrop. marked is the index of
	// the last emitted base already penalised.
	marked := -1
	// mark penalises the last emitted base, which precedes an error.
	mark := func() {
		if k := len(qual) - 1; k >= 0 && k != marked {
			qual[k] = clampQ(float64(qual[k]) - p.NeighbourDrop)
			marked = k
		}
	}
	afterErr := false // the next emitted base follows an error
	emit := func(b, q byte, isErr bool) {
		if isErr {
			mark()
		}
		seq = append(seq, b)
		qual = append(qual, q)
		if isErr {
			marked = len(qual) - 1
		} else if afterErr {
			qual[len(qual)-1] = clampQ(float64(q) - p.NeighbourDrop)
			marked = len(qual) - 1
		}
		afterErr = isErr
	}
	record := func(k OpKind, b byte) {
		if ops != nil {
			*ops = append(*ops, Op{Kind: k, Base: b})
		}
		if st != nil {
			switch k {
			case Match:
				st.Matches++
			case Sub:
				st.Subs++
			case Ins:
				st.Ins++
			case Del:
				st.Dels++
			}
		}
	}

	runEnd, runLen := 0, 0 // current homopolymer run is truth[runEnd-runLen:runEnd]
	for i := 0; i < n; i++ {
		t := truth[i]
		if i >= runEnd {
			j := i + 1
			for j < n && truth[j] == t {
				j++
			}
			runEnd, runLen = j, j-i
		}
		pDel := p.Del + m.HomopolymerExcess(runLen)
		u := r.Float64()
		switch {
		case u < pDel:
			mark()
			afterErr = true
			record(Del, 0)
		case u < pDel+p.Sub:
			b := otherBase(r, t)
			emit(b, errQ(base[i]), true)
			record(Sub, b)
		default:
			emit(t, base[i], false)
			record(Match, 0)
		}
		if r.Float64() < p.Ins {
			// Half of the insertions repeat the truth base (ONT insertions
			// favour extending the local base, Delahaye & Nicolas 2021); the
			// rest are uniform over ACGT.
			b := t
			if c := r.IntN(8); c < 4 {
				b = acgt[c]
			}
			emit(b, errQ(base[i]), true)
			record(Ins, b)
		}
	}
	if st != nil {
		st.TruthLen, st.ReadLen = n, len(seq)
		st.MeanQ = MeanQ(qual)
	}
	return seq, qual
}

// otherBase draws a base from ACGT different from t, uniformly (any of the
// four when t is not an upper-case ACGT base). It draws once.
func otherBase(r *rand.Rand, t byte) byte {
	const acgt = "ACGT"
	k := -1
	switch t {
	case 'A':
		k = 0
	case 'C':
		k = 1
	case 'G':
		k = 2
	case 'T':
		k = 3
	}
	if k < 0 {
		return acgt[r.IntN(4)]
	}
	return acgt[(k+1+r.IntN(3))&3]
}
