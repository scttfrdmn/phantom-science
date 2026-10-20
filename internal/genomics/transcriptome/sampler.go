package transcriptome

import (
	"bytes"
	"fmt"
	"math"
	"math/rand/v2"
	"sort"
	"strconv"

	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore"
	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore/errmodel"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
)

// TierNone disables sequencing errors: reads are the exact constructs, with
// SUP-tier qualities.
const TierNone = "none"

// Default molecule-model parameters. They are modelling choices informed by
// published ONT cDNA data, not fitted values:
//
//   - Full-length fraction 0.7. Reverse transcription and strand switching
//     often stop early, so a large share of cDNA reads lack the 5′ end; with
//     PCS kits the 3′ end is captured by the oligo-dT primer and coverage
//     decays toward the 5′ end (Sessegolo et al., Sci Rep 9:14908, 2019).
//     The realised fraction varies with sample quality and kit; 0.7 is a
//     modelling choice, not a fitted value.
//   - Truncated molecules keep a 3′ fragment of MinFragment + K bases, K
//     geometric with mean ≈ TruncMean (a constant per-base stop hazard),
//     truncated below the transcript length.
//   - Poly(A) tails of 20–100 nt: PolyAMin plus a geometric excess with
//     mean PolyAMean = 60 before truncation at PolyAMax, giving a mean of
//     ≈51 nt. Steady-state mammalian poly(A) tails are typically tens of
//     nucleotides to ≈100 nt, with medians of roughly 50–100 nt
//     (Chang et al., Mol Cell 53:1044, 2014; Subtelny et al., Nature
//     508:66, 2014).
//   - Novel-isoform rate 0.05 of reads.
const (
	DefaultFullLengthFraction = 0.7
	DefaultTruncMean          = 500
	DefaultMinFragment        = 100
	DefaultPolyAMin           = 20
	DefaultPolyAMax           = 100
	DefaultPolyAMean          = 60
	DefaultNovelRate          = 0.05
	DefaultTier               = "sup"
)

// alt3MinExon is the shortest last exon (or transcript, without structure)
// that gets an alternative 3′ end; the new end keeps 20–80% of it.
const alt3MinExon = 50

// Config configures a cDNA run. Start from DefaultConfig.
type Config struct {
	Reads int
	Seed  uint64 // must be non-zero; resolve a random seed first
	Tier  string // sup, hac, fast or none

	Kit Kit

	FullLengthFraction float64 // probability a molecule is full length
	TruncMean          float64 // mean of the geometric excess of truncated fragments
	MinFragment        int     // shortest truncated fragment, bases

	PolyAMin, PolyAMax int     // poly(A) bounds, bases
	PolyAMean          float64 // mean of the geometric excess over PolyAMin before truncation

	NovelRate float64 // probability a molecule comes from its transcript's novel isoform

	Run nanopore.RunInfo // ONT run metadata; unset fields default from Seed
}

// DefaultConfig returns the default configuration for reads reads.
func DefaultConfig(reads int, seed uint64) Config {
	return Config{
		Reads: reads, Seed: seed, Tier: DefaultTier, Kit: Kits[DefaultKit],
		FullLengthFraction: DefaultFullLengthFraction, TruncMean: DefaultTruncMean, MinFragment: DefaultMinFragment,
		PolyAMin: DefaultPolyAMin, PolyAMax: DefaultPolyAMax, PolyAMean: DefaultPolyAMean,
		NovelRate: DefaultNovelRate,
	}
}

// Validate checks the configuration (the tier is checked by NewSimulator).
func (c Config) Validate() error {
	switch {
	case c.Reads <= 0:
		return fmt.Errorf("reads must be positive, got %d", c.Reads)
	case c.Seed == 0:
		return fmt.Errorf("seed must be non-zero")
	case !(c.FullLengthFraction >= 0 && c.FullLengthFraction <= 1):
		return fmt.Errorf("full-length fraction must be in [0, 1], got %g", c.FullLengthFraction)
	case !(c.TruncMean > 0) || math.IsInf(c.TruncMean, 0):
		return fmt.Errorf("truncation mean must be positive, got %g", c.TruncMean)
	case c.MinFragment < 1:
		return fmt.Errorf("minimum fragment must be >= 1, got %d", c.MinFragment)
	case c.PolyAMin < 0 || c.PolyAMax < c.PolyAMin:
		return fmt.Errorf("poly(A) bounds must satisfy 0 <= min <= max, got %d..%d", c.PolyAMin, c.PolyAMax)
	case !(c.PolyAMean >= 0) || math.IsInf(c.PolyAMean, 0):
		return fmt.Errorf("poly(A) mean must be >= 0, got %g", c.PolyAMean)
	case !(c.NovelRate >= 0 && c.NovelRate <= 1):
		return fmt.Errorf("novel-isoform rate must be in [0, 1], got %g", c.NovelRate)
	}
	return c.Kit.validate()
}

// RNG streams derived from the run seed.
const (
	streamExpression = 1
	streamIsoforms   = 2
	streamReads      = 3
	streamHeaders    = 4
)

// Rand returns the deterministic generator for stream of the run seeded
// with seed (PCG with a splitmix64-scrambled stream ID).
func Rand(seed, stream uint64) *rand.Rand {
	// #nosec G404 -- simulation randomness, reproducible by design
	return rand.New(rand.NewPCG(seed, splitmix64(stream)))
}

func splitmix64(x uint64) uint64 {
	x += 0x9e3779b97f4a7c15
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9
	x = (x ^ (x >> 27)) * 0x94d049bb133111eb
	return x ^ (x >> 31)
}

// DefaultTPM draws the default Zipf/log-normal expression for t from the
// run seed's expression stream.
func DefaultTPM(t *Table, m ZipfModel, seed uint64) []float64 {
	tpm, _ := m.TPM(Rand(seed, streamExpression), len(t.Transcripts))
	return tpm
}

// Novel-isoform types.
const (
	NovelNone      = ""
	NovelSkipExon  = "skip_exon"
	NovelAlt3Prime = "alt_3prime"
)

// Species is one kind of molecule: an annotated transcript or the novel
// isoform derived from one.
type Species struct {
	ID     string
	Parent *Transcript
	Novel  string // NovelNone, NovelSkipExon or NovelAlt3Prime
	// Detail describes a novel isoform: "skipped_exon=<k>" (1-based, in
	// transcript order) or "end=<length>" (the new transcript length).
	Detail string
	Seq    []byte
	Exons  []Exon // genomic exons in transcript order; nil without structure
	TPM    float64
}

// IntronChain returns the introns as "start-end" (1-based, inclusive,
// ascending genomic order) joined by commas; "" for single-exon or
// unstructured species.
func (s *Species) IntronChain() string {
	if len(s.Exons) < 2 {
		return ""
	}
	ex := append([]Exon(nil), s.Exons...)
	sort.Slice(ex, func(i, j int) bool { return ex[i].Start < ex[j].Start })
	var b []byte
	for i := 1; i < len(ex); i++ {
		if i > 1 {
			b = append(b, ',')
		}
		b = strconv.AppendInt(b, ex[i-1].End+1, 10)
		b = append(b, '-')
		b = strconv.AppendInt(b, ex[i].Start-1, 10)
	}
	return string(b)
}

// Truth is the ground truth of one read.
type Truth struct {
	Species    *Species
	Reverse    bool // read is the reverse complement of the sense construct
	FullLength bool
	Trunc5p    int // 5′ bases missing from the molecule (= Start)
	PolyA      int
	Start, End int // 0-based half-open interval of Species.Seq in the read
}

// Read is one generated read.
type Read struct {
	Record seqread.Record
	Truth  Truth
	// Construct is the error-free molecule as sequenced (after strand
	// orientation), before the error model.
	Construct []byte
}

// Simulator generates cDNA reads.
type Simulator struct {
	cfg     Config
	table   *Table
	species []*Species
	cum     []float64
	errs    *errmodel.Model // nil for TierNone
	quals   *errmodel.Model
	vnpRC   []byte
	info    nanopore.RunInfo
	maxLen  int
}

// NewSimulator prepares a run over table with per-transcript TPMs (one per
// table transcript, any positive scale).
func NewSimulator(table *Table, tpm []float64, cfg Config) (*Simulator, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	if len(tpm) != len(table.Transcripts) {
		return nil, fmt.Errorf("got %d TPMs for %d transcripts", len(tpm), len(table.Transcripts))
	}
	s := &Simulator{cfg: cfg, table: table, vnpRC: ReverseComplement([]byte(cfg.Kit.VNP))}
	var err error
	if cfg.Tier == TierNone {
		s.quals, err = errmodel.New(errmodel.Sup)
	} else {
		var t errmodel.Tier
		if t, err = errmodel.ParseTier(cfg.Tier); err == nil {
			s.errs, err = errmodel.New(t)
			s.quals = s.errs
		}
	}
	if err != nil {
		return nil, err
	}
	s.info = cfg.Run.WithDefaults(cfg.Seed)
	if err := s.info.Validate(); err != nil {
		return nil, err
	}

	total := 0.0
	for i, v := range tpm {
		if !(v >= 0) || math.IsInf(v, 0) {
			return nil, fmt.Errorf("transcript %s: bad TPM %g", table.Transcripts[i].ID, v)
		}
		total += v
	}
	if !(total > 0) {
		return nil, fmt.Errorf("total expression is zero")
	}
	var novel []*Species
	if cfg.NovelRate > 0 {
		novel = designNovel(Rand(cfg.Seed, streamIsoforms), table)
	}
	sum := 0.0
	for i, tx := range table.Transcripts {
		w := float64(tpm[i] / total * 1e6)
		sp := &Species{ID: tx.ID, Parent: tx, Seq: tx.Seq, Exons: tx.Exons, TPM: w}
		s.species = append(s.species, sp)
		if novel != nil && novel[i] != nil {
			n := novel[i]
			n.TPM = float64(w * cfg.NovelRate)
			sp.TPM = w - n.TPM
			s.species = append(s.species, n)
		}
	}
	s.cum = make([]float64, len(s.species))
	for i, sp := range s.species {
		sum += sp.TPM
		s.cum[i] = sum
		s.maxLen = max(s.maxLen, len(sp.Seq))
	}
	return s, nil
}

// Species returns every molecule species, annotated transcripts in table
// order each followed by its novel isoform (if any).
func (s *Simulator) Species() []*Species { return s.species }

// RunInfo returns the run metadata (defaults filled in).
func (s *Simulator) RunInfo() nanopore.RunInfo { return s.info }

// Config returns the configuration.
func (s *Simulator) Config() Config { return s.cfg }

// slotBases sizes the ONT channel model's time windows.
func (s *Simulator) slotBases() int {
	k := s.cfg.Kit
	return s.maxLen + s.cfg.PolyAMax + len(k.SSP) + len(k.VNP) + len(k.Adapter5) + len(k.Adapter3)
}

// designNovel returns, per table transcript, its novel isoform or nil. A
// transcript with ≥ 3 exons skips one internal exon (chosen at random among
// those whose skipped sequence differs from every annotated transcript and
// every novel isoform already designed for the gene); otherwise, or when no
// skip is new, its 3′ end moves to a random point keeping 20–80% of the last
// exon (of the whole transcript without structure), if that exon has at
// least alt3MinExon bases and the result is new.
func designNovel(r *rand.Rand, t *Table) []*Species {
	byGene := map[string][]*Transcript{}
	for _, tx := range t.Transcripts {
		byGene[tx.GeneID] = append(byGene[tx.GeneID], tx)
	}
	novelByGene := map[string][][]byte{}
	isNew := func(tx *Transcript, seq []byte) bool {
		for _, o := range byGene[tx.GeneID] {
			if bytes.Equal(o.Seq, seq) {
				return false
			}
		}
		for _, o := range novelByGene[tx.GeneID] {
			if bytes.Equal(o, seq) {
				return false
			}
		}
		return true
	}
	out := make([]*Species, len(t.Transcripts))
	for i, tx := range t.Transcripts {
		if n := len(tx.Exons); n >= 3 {
			off := tx.exonOffsets()
			for _, j := range r.Perm(n - 2) {
				k := j + 1
				seq := make([]byte, 0, len(tx.Seq)-tx.Exons[k].Len())
				seq = append(seq, tx.Seq[:off[k]]...)
				seq = append(seq, tx.Seq[off[k]+tx.Exons[k].Len():]...)
				if !isNew(tx, seq) {
					continue
				}
				exons := append(append([]Exon(nil), tx.Exons[:k]...), tx.Exons[k+1:]...)
				out[i] = &Species{
					ID: tx.ID + "_skip" + strconv.Itoa(k+1), Parent: tx, Novel: NovelSkipExon,
					Detail: "skipped_exon=" + strconv.Itoa(k+1), Seq: seq, Exons: exons,
				}
				break
			}
			if out[i] != nil {
				novelByGene[tx.GeneID] = append(novelByGene[tx.GeneID], out[i].Seq)
				continue
			}
		}
		e := len(tx.Seq)
		if len(tx.Exons) > 0 {
			e = tx.Exons[len(tx.Exons)-1].Len()
		}
		if e < alt3MinExon {
			continue
		}
		lo, hi := (e*2+9)/10, e*8/10
		keep := lo + r.IntN(hi-lo+1)
		n := len(tx.Seq) - e + keep
		seq := tx.Seq[:n:n]
		if !isNew(tx, seq) {
			continue
		}
		var exons []Exon
		if len(tx.Exons) > 0 {
			exons = append([]Exon(nil), tx.Exons...)
			last := &exons[len(exons)-1]
			if tx.Strand == '-' {
				last.Start = last.End - int64(keep) + 1
			} else {
				last.End = last.Start + int64(keep) - 1
			}
		}
		out[i] = &Species{
			ID: tx.ID + "_alt3", Parent: tx, Novel: NovelAlt3Prime,
			Detail: "end=" + strconv.Itoa(n), Seq: seq, Exons: exons,
		}
		novelByGene[tx.GeneID] = append(novelByGene[tx.GeneID], seq)
	}
	return out
}

// truncGeom draws k in {0..K} with P(k) ∝ q^k (q in [0, 1)) by inversion.
func truncGeom(r *rand.Rand, q float64, K int) int {
	if K <= 0 || q <= 0 {
		return 0
	}
	u := r.Float64()
	lq := math.Log(q)
	// P(k <= x) = (1 − q^{x+1}) / (1 − q^{K+1}).
	tail := -math.Expm1(float64(K+1) * lq) // 1 − q^{K+1}
	k := int(math.Floor(math.Log1p(-u*tail) / lq))
	return min(max(k, 0), K)
}

// geomQ returns q for a geometric distribution on {0, 1, …} with mean m.
func geomQ(m float64) float64 { return m / (1 + m) }

// fragment draws the molecule's 5′ truncation for a species of length L.
func (s *Simulator) fragment(r *rand.Rand, L int) (start int, full bool) {
	u := r.Float64()
	if u < s.cfg.FullLengthFraction || L-1 < s.cfg.MinFragment {
		return 0, true
	}
	// Keep MinFragment + k bases, k in 0..L−1−MinFragment.
	k := truncGeom(r, math.Exp(-1/s.cfg.TruncMean), L-1-s.cfg.MinFragment)
	return L - (s.cfg.MinFragment + k), false
}

// polyA draws a poly(A) length.
func (s *Simulator) polyA(r *rand.Rand) int {
	return s.cfg.PolyAMin + truncGeom(r, geomQ(s.cfg.PolyAMean), s.cfg.PolyAMax-s.cfg.PolyAMin)
}

// pick returns a species drawn by abundance.
func (s *Simulator) pick(r *rand.Rand) *Species {
	u := r.Float64() * s.cum[len(s.cum)-1]
	i := sort.Search(len(s.cum), func(i int) bool { return s.cum[i] > u })
	return s.species[min(i, len(s.species)-1)]
}

// construct returns the sense molecule for a fragment.
func (s *Simulator) construct(sp *Species, start, polyA int) []byte {
	k := s.cfg.Kit
	frag := sp.Seq[start:]
	b := make([]byte, 0, len(k.Adapter5)+len(k.SSP)+len(frag)+polyA+len(s.vnpRC)+len(k.Adapter3))
	b = append(b, k.Adapter5...)
	b = append(b, k.SSP...)
	b = append(b, frag...)
	b = append(b, bytes.Repeat([]byte{'A'}, polyA)...)
	b = append(b, s.vnpRC...)
	b = append(b, k.Adapter3...)
	return b
}

// Generate produces cfg.Reads reads in order, calling fn for each. The Read
// (and its slices) may be retained by fn. Output depends only on the table,
// TPMs and configuration.
func (s *Simulator) Generate(fn func(*Read) error) error {
	r := Rand(s.cfg.Seed, streamReads)
	hr := Rand(s.cfg.Seed, streamHeaders)
	header := nanopore.NewHeaderFormatter(s.info, s.cfg.Reads, s.slotBases())
	for i := range s.cfg.Reads {
		sp := s.pick(r)
		start, full := s.fragment(r, len(sp.Seq))
		pa := s.polyA(r)
		rev := r.IntN(2) == 1
		mol := s.construct(sp, start, pa)
		if rev {
			mol = ReverseComplement(mol)
		}
		var seq, qual []byte
		if s.errs == nil {
			seq, qual = mol, s.quals.Qualities(r, len(mol))
		} else {
			seq, qual, _, _ = s.errs.Apply(r, mol)
		}
		rd := &Read{
			Record: seqread.Record{Seq: seq, Qual: qual, Batch: 1, Mate: 1, Index: i + 1, Serial: int64(i)},
			Truth: Truth{
				Species: sp, Reverse: rev, FullLength: full, Trunc5p: start, PolyA: pa,
				Start: start, End: len(sp.Seq),
			},
			Construct: mol,
		}
		rd.Record.Name, rd.Record.Comment = header(hr, &rd.Record)
		if err := fn(rd); err != nil {
			return err
		}
	}
	return nil
}
