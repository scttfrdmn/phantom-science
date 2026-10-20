package signal

import (
	"errors"
	"fmt"
	"hash/fnv"
	"math"
	"math/rand/v2"
)

// Digitiser converts between pA and raw ADC values:
//
//	pA  = (raw + Offset) × Range / Digitisation
//	raw = round(pA × Digitisation / Range − Offset), clamped to int16
//
// These are the slow5/pod5 calibration fields of the same names.
type Digitiser struct {
	Digitisation, Offset, Range float64
}

// Raw digitises pA, reporting whether the value was clamped to the int16
// range.
func (d Digitiser) Raw(pA float64) (raw int16, clipped bool) {
	v := math.Round(pA*d.Digitisation/d.Range - d.Offset)
	switch {
	case v > math.MaxInt16:
		return math.MaxInt16, true
	case v < math.MinInt16:
		return math.MinInt16, true
	case math.IsNaN(v):
		return 0, true
	}
	return int16(v), false
}

// PA converts a raw value back to pA. For unclamped values
// |PA(Raw(x)) − x| ≤ Range/(2·Digitisation).
func (d Digitiser) PA(raw int16) float64 {
	return (float64(raw) + d.Offset) * d.Range / d.Digitisation
}

// EventKind labels a ground-truth segment.
type EventKind uint8

// Segment kinds, in the order they appear in a read.
const (
	EventLeader   EventKind = iota // adapter leader plateau
	EventAdapter                   // an adapter k-mer
	EventPolyTail                  // poly-A/poly-T plateau
	EventKmer                      // a k-mer of the read sequence
)

func (k EventKind) String() string {
	switch k {
	case EventLeader:
		return "leader"
	case EventAdapter:
		return "adapter"
	case EventPolyTail:
		return "polytail"
	case EventKmer:
		return "kmer"
	}
	return fmt.Sprintf("EventKind(%d)", uint8(k))
}

// Event is one ground-truth segment of a synthetic read: samples
// Raw[Start:Start+Len] were drawn around LevelPA with per-sample SD NoisePA
// (plus the slow drift).
type Event struct {
	Kind    EventKind
	Pos     int // k-mer start in the read (EventKmer) or adapter (EventAdapter); -1 for plateaus
	Kmer    int // k-mer index (KmerIndex); -1 for plateaus
	Start   int // first sample
	Len     int // number of samples (0 for a skip)
	LevelPA float64
	NoisePA float64
	Stall   bool
	Skip    bool
}

// Read is a synthetic raw-signal read. The calibration fields map directly
// to slow5/pod5 read fields (read_id, digitisation, offset, range,
// sampling_rate, len_raw_signal = NumSamples, raw_signal = Raw,
// median_before).
type Read struct {
	ID           string
	Raw          []int16
	Digitisation float64
	Offset       float64
	Range        float64
	SamplingRate float64
	MedianBefore float64 // open-pore current before the read, pA
	NumSamples   int     // len(Raw)
	Duration     float64 // seconds, NumSamples / SamplingRate

	// Segment layout, as sample indices into Raw.
	AdapterSamples int // [0, AdapterSamples) is leader + adapter
	PolyStart      int // [PolyStart, PolyEnd) is the poly-A/T plateau (empty if none)
	PolyEnd        int
	SeqStart       int // [SeqStart, SeqEnd) is the read sequence
	SeqEnd         int

	Kmers   int // read k-mers walked (len(seq) − K + 1)
	Stalls  int // read k-mers that stalled
	Skips   int // read k-mers skipped
	Clipped int // samples clamped to the int16 range

	// Events is the per-segment ground truth, recorded only with
	// WithEvents.
	Events []Event
}

// Digitiser returns the read's ADC conversion.
func (r *Read) Digitiser() Digitiser {
	return Digitiser{Digitisation: r.Digitisation, Offset: r.Offset, Range: r.Range}
}

// PA returns the read's signal in pA.
func (r *Read) PA() []float64 {
	d := r.Digitiser()
	out := make([]float64, len(r.Raw))
	for i, v := range r.Raw {
		out[i] = d.PA(v)
	}
	return out
}

// PolyTail is a homopolymer plateau (poly-A tail of direct RNA or cDNA).
type PolyTail struct {
	Length int  // tail length in bases; the plateau lasts Length × mean dwell samples
	Base   byte // 'A' (default) or 'T': the plateau level is the model level of the K-mer of Base
	AtEnd  bool // place after the read sequence instead of before it
}

// Option configures one Synthesize call.
type Option func(*synthConfig)

type synthConfig struct {
	id             string
	events         bool
	poly           PolyTail
	adapterSeconds float64
	adapterSet     bool
}

// WithID sets Read.ID.
func WithID(id string) Option { return func(c *synthConfig) { c.id = id } }

// WithEvents records per-segment ground truth in Read.Events.
func WithEvents() Option { return func(c *synthConfig) { c.events = true } }

// WithPolyA adds a poly-A plateau of n bases between the adapter and the
// read sequence: where the 3' poly-A tail of a direct-RNA read appears,
// since RNA translocates 3'→5'.
func WithPolyA(n int) Option { return WithPolyTail(PolyTail{Length: n, Base: 'A'}) }

// WithPolyTail adds a homopolymer plateau (see PolyTail). For cDNA, a
// forward-strand read ends in poly-A (AtEnd), a reverse-strand read starts
// with poly-T.
func WithPolyTail(t PolyTail) Option { return func(c *synthConfig) { c.poly = t } }

// WithAdapterSeconds overrides Profile.AdapterSeconds for this read, e.g. to
// match an adapter time already chosen by nanopore.ChannelModel.
func WithAdapterSeconds(s float64) Option {
	return func(c *synthConfig) { c.adapterSeconds, c.adapterSet = s, true }
}

// readStream is the PCG stream of per-read generators.
const readStream = 0x7369676e616c // "signal"

// ReadRand returns the generator for one read: PCG seeded with
// runSeed XOR FNV-1a-64(readID). It depends only on (runSeed, readID).
func ReadRand(runSeed uint64, readID string) *rand.Rand {
	h := fnv.New64a()
	_, _ = h.Write([]byte(readID))
	// #nosec G404 -- simulation, reproducible by design
	return rand.New(rand.NewPCG(runSeed^h.Sum64(), readStream))
}

// SynthesizeRead synthesises the signal of read readID with a generator
// derived from (runSeed, readID) (see ReadRand), so a read's signal can be
// regenerated from its ID without storing it, in any order or concurrently.
// Read.ID is set to readID.
func SynthesizeRead(runSeed uint64, readID string, seq []byte, p Profile, model *PoreModel, opts ...Option) (*Read, error) {
	return Synthesize(ReadRand(runSeed, readID), seq, p, model, append([]Option{WithID(readID)}, opts...)...)
}

// ErrShortRead is returned for sequences shorter than the model's k.
var ErrShortRead = errors.New("sequence shorter than the pore model k")

// Synthesize turns seq into a raw-signal read, drawing only from r.
//
// Bases other than A/C/G/T/U (either case) — N or IUPAC codes — are replaced
// by a uniformly random base drawn from r (in sequence order, before any
// other per-base draw), so a masked base still produces plausible signal.
// seq is not modified.
//
// Per read it draws, in order: the ADC offset ~ N(OffsetMean, OffsetStd),
// median_before ~ N(MedianBeforeMean, MedianBeforeStd), a level scale
// 1 + N(0, ReadScaleSD) and shift N(0, ReadShiftSD), and the initial
// drift. Every k-mer level μ and stdv σ from the model becomes
// μ·scale + shift and σ·LevelNoiseScale·scale; each sample is
// N(level, sqrt(σ'² + WhiteNoisePA²)) plus the drift. With LevelTau > 0 the
// sample means approach each new level exponentially (see Profile.LevelTau);
// this draws nothing extra from r.
//
// Dwell times are gamma(mean = SamplingRate/SpeedBPS, CV = DwellCV); a stall
// multiplies the dwell by U(StallMin, StallMax), a skip makes it zero.
// Dwells are accumulated as continuous time and each k-mer gets the samples
// whose index falls in its interval, so the mean samples per base is exact
// in expectation (no rounding bias) and a skip gives exactly zero samples.
func Synthesize(r *rand.Rand, seq []byte, p Profile, model *PoreModel, opts ...Option) (*Read, error) {
	if err := p.Validate(); err != nil {
		return nil, err
	}
	if model == nil || model.K < 1 || len(model.Mean) != model.NumKmers() || len(model.Stdv) != model.NumKmers() {
		return nil, fmt.Errorf("signal: invalid pore model")
	}
	cfg := synthConfig{}
	for _, o := range opts {
		o(&cfg)
	}
	adapterSeconds := p.AdapterSeconds
	if cfg.adapterSet {
		adapterSeconds = cfg.adapterSeconds
	}
	if !(adapterSeconds >= 0) || math.IsInf(adapterSeconds, 0) {
		return nil, fmt.Errorf("signal: adapter seconds must be non-negative, got %g", adapterSeconds)
	}
	k := model.K
	if len(seq) < k {
		return nil, fmt.Errorf("signal: %w (%d < %d)", ErrShortRead, len(seq), k)
	}
	if cfg.poly.Length < 0 {
		return nil, fmt.Errorf("signal: negative poly tail length %d", cfg.poly.Length)
	}
	polyBase := cfg.poly.Base
	if polyBase == 0 {
		polyBase = 'A'
	}
	if cfg.poly.Length > 0 && baseCode[polyBase] < 0 {
		return nil, fmt.Errorf("signal: poly tail base %q is not a nucleotide", polyBase)
	}

	s := synth{r: r, p: p, model: model, events: cfg.events}
	s.seq = sanitize(r, seq)

	offset := p.OffsetMean + p.OffsetStd*r.NormFloat64()
	medianBefore := p.MedianBeforeMean + p.MedianBeforeStd*r.NormFloat64()
	s.scale = 1 + p.ReadScaleSD*r.NormFloat64()
	if s.scale < 0.5 {
		s.scale = 0.5
	}
	s.shift = p.ReadShiftSD * r.NormFloat64()
	if p.DriftPA > 0 {
		s.driftA = math.Exp(-1 / (p.DriftTau * p.SamplingRate))
		s.driftB = p.DriftPA * math.Sqrt(1-s.driftA*s.driftA)
		s.drift = p.DriftPA * r.NormFloat64()
	}
	if p.LevelTau > 0 {
		s.levelA = math.Exp(-1 / p.LevelTau)
	}
	s.dig = p.Digitiser(offset)
	nk := len(s.seq) - k + 1
	est := int(float64(nk)*p.MeanDwell()*1.1) + int(adapterSeconds*p.SamplingRate) + int(float64(cfg.poly.Length)*p.MeanDwell()) + 16
	s.raw = make([]int16, 0, est)
	if cfg.events {
		s.ev = make([]Event, 0, nk+64)
	}

	out := &Read{
		ID: cfg.id, Digitisation: p.Digitisation, Offset: offset, Range: p.Range,
		SamplingRate: p.SamplingRate, MedianBefore: medianBefore, Kmers: nk,
	}

	// Adapter: leader plateau, then the adapter k-mers, totalling
	// round(adapterSeconds × SamplingRate) samples (the adapter k-mers alone
	// if they take longer).
	if target := int(math.Round(adapterSeconds * p.SamplingRate)); target > 0 {
		ad := sanitize(r, []byte(p.AdapterSeq))
		var dw []int
		if len(ad) >= k {
			dw = s.dwells(len(ad)-k+1, false)
		}
		used := 0
		for _, d := range dw {
			used += d
		}
		if lead := target - used; lead > 0 {
			s.plateau(EventLeader, lead, p.AdapterLevelPA*s.scale+s.shift, p.AdapterNoisePA)
		}
		for i, d := range dw {
			s.kmer(EventAdapter, ad, i, d, false, false)
		}
	}
	out.AdapterSamples = len(s.raw)

	polyLen := int(math.Round(float64(cfg.poly.Length) * p.MeanDwell()))
	polyLevel := 0.0
	if polyLen > 0 {
		hp := make([]byte, k)
		for i := range hp {
			hp[i] = polyBase
		}
		polyLevel = float64(model.Mean[KmerIndex(hp)])*s.scale + s.shift
	}
	out.PolyStart, out.PolyEnd = len(s.raw), len(s.raw)
	if polyLen > 0 && !cfg.poly.AtEnd {
		s.plateau(EventPolyTail, polyLen, polyLevel, p.PolyANoisePA)
		out.PolyEnd = len(s.raw)
	}

	out.SeqStart = len(s.raw)
	dw := s.dwells(nk, true)
	for j := 0; j < nk; j++ {
		pos := j
		if p.Reverse {
			pos = nk - 1 - j
		}
		s.kmer(EventKmer, s.seq, pos, dw[j], s.stall[j], s.skip[j])
		if s.stall[j] {
			out.Stalls++
		}
		if s.skip[j] {
			out.Skips++
		}
	}
	out.SeqEnd = len(s.raw)

	if polyLen > 0 && cfg.poly.AtEnd {
		out.PolyStart = len(s.raw)
		s.plateau(EventPolyTail, polyLen, polyLevel, p.PolyANoisePA)
		out.PolyEnd = len(s.raw)
	}

	out.Raw = s.raw
	out.NumSamples = len(s.raw)
	out.Duration = float64(out.NumSamples) / p.SamplingRate
	out.Clipped = s.clipped
	out.Events = s.ev
	return out, nil
}

// sanitize returns seq with every non-ACGTU byte replaced by a random base.
// It returns seq itself when no replacement is needed.
func sanitize(r *rand.Rand, seq []byte) []byte {
	var out []byte
	for i, b := range seq {
		if baseCode[b] >= 0 {
			continue
		}
		if out == nil {
			out = append([]byte(nil), seq...)
		}
		out[i] = "ACGT"[r.IntN(4)]
	}
	if out == nil {
		return seq
	}
	return out
}

type synth struct {
	r            *rand.Rand
	p            Profile
	model        *PoreModel
	seq          []byte
	scale, shift float64
	drift        float64
	driftA       float64
	driftB       float64
	// level is the smoothed mean level of the last sample emitted and
	// levelA = exp(−1/LevelTau) its per-sample decay (0 when LevelTau is 0);
	// levelSet is false until the first segment.
	level       float64
	levelA      float64
	levelSet    bool
	dig         Digitiser
	raw         []int16
	clipped     int
	events      bool
	ev          []Event
	stall, skip []bool
}

// dwells draws n k-mer dwells (in samples). With events=true it also draws
// stalls and skips and records them in s.stall / s.skip.
func (s *synth) dwells(n int, withEvents bool) []int {
	mean := s.p.MeanDwell()
	out := make([]int, n)
	if withEvents {
		s.stall, s.skip = make([]bool, n), make([]bool, n)
	}
	t := 0.0
	for i := range out {
		d := gammaMeanCV(s.r, mean, s.p.DwellCV)
		if withEvents && (s.p.StallProb > 0 || s.p.SkipProb > 0) {
			u := s.r.Float64()
			switch {
			case u < s.p.SkipProb:
				d = 0
				s.skip[i] = true
			case u < s.p.SkipProb+s.p.StallProb:
				d *= s.p.StallMin + (s.p.StallMax-s.p.StallMin)*s.r.Float64()
				s.stall[i] = true
			}
		}
		next := t + d
		out[i] = int(math.Floor(next)) - int(math.Floor(t))
		t = next
	}
	return out
}

func (s *synth) emit(level, sd float64) {
	v := level + sd*s.r.NormFloat64()
	if s.driftB > 0 {
		s.drift = s.driftA*s.drift + s.driftB*s.r.NormFloat64()
		v += s.drift
	}
	raw, clipped := s.dig.Raw(v)
	if clipped {
		s.clipped++
	}
	s.raw = append(s.raw, raw)
}

// segment emits n samples around target with SD sd. With LevelTau > 0 the
// mean approaches target from the previous segment's last level (see
// Profile.LevelTau); the first segment of a read starts at its target.
func (s *synth) segment(n int, target, sd float64) {
	if s.levelA == 0 {
		for i := 0; i < n; i++ {
			s.emit(target, sd)
		}
		return
	}
	if !s.levelSet {
		s.level, s.levelSet = target, true
	}
	for i := 0; i < n; i++ {
		s.level = target + (s.level-target)*s.levelA
		s.emit(s.level, sd)
	}
}

func (s *synth) plateau(kind EventKind, n int, level, sd float64) {
	start := len(s.raw)
	s.segment(n, level, sd)
	if s.events {
		s.ev = append(s.ev, Event{Kind: kind, Pos: -1, Kmer: -1, Start: start, Len: n, LevelPA: level, NoisePA: sd})
	}
}

func (s *synth) kmer(kind EventKind, seq []byte, pos, n int, stall, skip bool) {
	idx := KmerIndex(seq[pos : pos+s.model.K])
	level := float64(s.model.Mean[idx])*s.scale + s.shift
	lsd := float64(s.model.Stdv[idx]) * s.p.LevelNoiseScale * s.scale
	sd := math.Sqrt(lsd*lsd + s.p.WhiteNoisePA*s.p.WhiteNoisePA)
	start := len(s.raw)
	s.segment(n, level, sd)
	if s.events {
		s.ev = append(s.ev, Event{Kind: kind, Pos: pos, Kmer: idx, Start: start, Len: n, LevelPA: level, NoisePA: sd, Stall: stall, Skip: skip})
	}
}

// gammaMeanCV draws from a gamma distribution with the given mean and
// coefficient of variation (shape 1/cv², scale mean·cv²). cv = 0 returns
// mean.
func gammaMeanCV(r *rand.Rand, mean, cv float64) float64 {
	if cv <= 0 {
		return mean
	}
	shape := 1 / (cv * cv)
	return gammaSample(r, shape) * mean / shape
}

// gammaSample draws Gamma(shape, 1) by Marsaglia & Tsang (2000), boosting
// shape < 1 with U^(1/shape).
func gammaSample(r *rand.Rand, shape float64) float64 {
	if shape < 1 {
		u := r.Float64()
		return gammaSample(r, shape+1) * math.Pow(u, 1/shape)
	}
	d := shape - 1.0/3
	c := 1 / math.Sqrt(9*d)
	for {
		x := r.NormFloat64()
		v := 1 + c*x
		if v <= 0 {
			continue
		}
		v = v * v * v
		u := r.Float64()
		if u < 1-0.0331*x*x*x*x || math.Log(u) < 0.5*x*x+d*(1-v+math.Log(v)) {
			return d * v
		}
	}
}
