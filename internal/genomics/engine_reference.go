package genomics

// Reads sampled from a reference genome (#27).
//
// With SequencerConfig.ReferencePath set, the engine loads one
// reference.Sampler (the whole FASTA, or only ReferenceRegion) when it is
// created and plugs Sampler.Fragment into the fragment hook
// (FragmentSampler), so every single-end read and every paired-end fragment
// is a reference fragment:
//
//   - single-end: a fragment of the sampled read length L;
//   - paired-end: a fragment of the sampled insert size F, R1 = frag[:L],
//     R2 = revcomp(frag[F−L:]);
//   - a minus-strand fragment is the reverse complement of the reference.
//
// The sampler draws from the pipeline's shared stream (contig, position,
// strand; see reference.Sampler.Fragment), so output stays a pure function
// of the seed and the same for any number of workers. A reference cannot be
// combined with an HMM model: bases come from one or the other.
//
// Positioned variants (GenerateOpts.Variants with CHROM/POS) are applied in
// reference coordinates: a variant is applied to a fragment when its
// reference span lies inside the fragment (plus the slack its deletions
// pull in), on the forward strand before a minus-strand fragment is
// reverse-complemented, right to left, skipping overlaps, as in region mode.
// Each variant's REF is checked against the reference once per call: a
// mismatch is skipped and counted in one warning, or with StrictRef fails
// the call before anything is written. Position-blind variants are not used.
//
// File output also writes a truth table per output unit, next to the R1
// (or single-end) file: one line per read (pair) with the read name and the
// fragment's origin, recorded in the manifest with role "truth".

import (
	"bufio"
	"errors"
	"fmt"
	"math/rand/v2"
	"sort"
	"strconv"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
)

// errNoReference is returned by ReferenceBases without a reference.
var errNoReference = errors.New("no reference genome configured (reference_path / --reference)")

// loadReference validates the reference settings and loads the sampler.
func (e *Engine) loadReference() {
	c := e.config
	switch {
	case c.ReferencePath == "" && c.ReferenceRegion != "":
		e.refErr = fmt.Errorf("reference_region %q needs reference_path", c.ReferenceRegion)
		return
	case c.ReferencePath == "":
		return
	case resolveModelPath(c) != "":
		e.refErr = fmt.Errorf("reference_path %s cannot be combined with an HMM model (%s): bases come from the reference or from the model, not both",
			c.ReferencePath, resolveModelPath(c))
		return
	}
	s, err := reference.NewSampler(c.ReferencePath, c.ReferenceRegion)
	if err != nil {
		e.refErr = fmt.Errorf("reference: %w", err)
		return
	}
	e.ref = s
	e.fragments = s.Fragment
	where := c.ReferencePath
	if c.ReferenceRegion != "" {
		where += " region " + c.ReferenceRegion
	}
	fmt.Fprintf(e.logw, "   [info] Sampling reads from reference %s\n", where)
	if !e.fixedLengthReads() && !c.PairedEnd {
		if e.refMaxRead, err = e.longestCleanStretch(); err != nil {
			e.refErr = fmt.Errorf("reference: %w", err)
			return
		}
		fmt.Fprintf(e.logw, "   [info] Read lengths above %d bases (the longest stretch free of N) are redrawn\n", e.refMaxRead)
	}
}

// maxLengthRedraws bounds the redraws of a read length that no reference
// stretch can supply; the length is then capped.
const maxLengthRedraws = 1000

// fitLengths wraps a read-length sampler so that, in reference mode with a
// length distribution, lengths no N-free reference stretch can supply are
// redrawn: read lengths follow the distribution truncated at the longest
// stretch (relevant for long reads from a small region).
func (e *Engine) fitLengths(sample func(r *rand.Rand) int) func(r *rand.Rand) int {
	limit := e.refMaxRead
	if limit <= 0 {
		return sample
	}
	return func(r *rand.Rand) int {
		for range maxLengthRedraws {
			if n := sample(r); int64(n) <= limit {
				return n
			}
		}
		return int(limit)
	}
}

// longestCleanStretch returns the longest run of A/C/G/T in the loaded
// reference (the region, when one is set), scanning it in windows.
func (e *Engine) longestCleanStretch() (int64, error) {
	var best int64
	for _, c := range e.ref.Contigs() {
		lo, hi, ok := e.loadedRange(c)
		if !ok {
			continue
		}
		var run int64
		const window = 4 << 20
		for a := lo; a < hi; a += window {
			b, err := e.ref.Bases(c.Name, a, min(a+window, hi))
			if err != nil {
				return 0, err
			}
			for _, x := range b {
				switch x {
				case 'A', 'C', 'G', 'T':
					run++
					best = max(best, run)
				default:
					run = 0
				}
			}
		}
	}
	return best, nil
}

// loadedRange returns the 0-based half-open range of contig c the sampler
// holds: all of it, or the configured region's part of it.
func (e *Engine) loadedRange(c reference.Contig) (int64, int64, bool) {
	region := e.config.ReferenceRegion
	if region == "" || c.Length == 0 {
		return 0, c.Length, c.Length > 0
	}
	if _, err := e.ref.Bases(c.Name, 0, c.Length); err == nil {
		return 0, c.Length, true // the region is this whole contig
	}
	colon := strings.LastIndexByte(region, ':')
	if colon < 0 || chromKey(region[:colon]) != chromKey(c.Name) {
		return 0, 0, false
	}
	from, to, _ := strings.Cut(strings.ReplaceAll(region[colon+1:], ",", ""), "-")
	first, err := strconv.ParseInt(from, 10, 64)
	if err != nil {
		return 0, 0, false
	}
	last := c.Length
	if to != "" {
		if last, err = strconv.ParseInt(to, 10, 64); err != nil {
			return 0, 0, false
		}
	}
	return first - 1, last, true
}

// ReferenceBases returns the upper-cased forward-strand bases [start, end)
// (0-based, half-open) of chrom in the configured reference (with or without
// a "chr" prefix). With reference_region set, the range must lie inside the
// region. It fails when no reference is configured.
func (e *Engine) ReferenceBases(chrom string, start, end int64) ([]byte, error) {
	if e.refErr != nil {
		return nil, e.refErr
	}
	if e.ref == nil {
		return nil, errNoReference
	}
	return e.ref.Bases(chrom, start, end)
}

// chromKey normalises a contig name for matching variants to fragments:
// without a "chr" prefix, and M for MT.
func chromKey(name string) string {
	k := strings.TrimPrefix(name, "chr")
	if k == "MT" {
		return "M"
	}
	return k
}

// posVariants are positioned variants on one contig, sorted by descending
// position, with their reference spans.
type posVariants struct {
	vars     []VariantRecord
	spans    []seqread.Span
	maxSlack int // total deletion length of all the variants
}

// newPosVariants sorts vs (positioned variants on one contig) by descending
// position and indexes them.
func newPosVariants(vs []VariantRecord) *posVariants {
	pv := &posVariants{vars: vs}
	sort.SliceStable(pv.vars, func(i, j int) bool { return pv.vars[i].Position > pv.vars[j].Position })
	for _, v := range pv.vars {
		pv.spans = append(pv.spans, seqread.Span{Pos: v.Position - 1, RefLen: len(v.Ref), AltLen: len(v.Alt)})
		pv.maxSlack += max(0, len(v.Ref)-len(v.Alt))
	}
	return pv
}

// window returns the variants (and their spans) whose 0-based start lies in
// [lo, hi), still in descending order.
func (pv *posVariants) window(lo, hi int64) ([]VariantRecord, []seqread.Span) {
	if pv == nil {
		return nil, nil
	}
	a := sort.Search(len(pv.vars), func(i int) bool { return pv.vars[i].Position-1 < hi })
	b := sort.Search(len(pv.vars), func(i int) bool { return pv.vars[i].Position-1 < lo })
	return pv.vars[a:b], pv.spans[a:b]
}

// refVariants are the positioned variants of one call in reference mode.
type refVariants struct {
	byChrom map[string]*posVariants // chromKey -> variants
	// mismatched counts variants whose REF differs from the reference
	// (skipped); outside those outside the loaded reference (never
	// reached); blind those without CHROM/POS (not used).
	mismatched, outside, blind int
	firstMismatch              string
}

// referenceVariants checks variants against the reference and indexes the
// usable ones by contig. With strict_ref, a REF mismatch is an error.
func (e *Engine) referenceVariants(variants []VariantRecord) (*refVariants, error) {
	rv := &refVariants{byChrom: map[string]*posVariants{}}
	if e.ref == nil {
		return rv, nil
	}
	groups := map[string][]VariantRecord{}
	for _, v := range variants {
		if len(v.Ref) == 0 || len(v.Alt) == 0 {
			continue
		}
		if v.Chromosome == "" || v.Position <= 0 {
			rv.blind++
			continue
		}
		got, err := e.ref.Bases(v.Chromosome, v.Position-1, v.Position-1+int64(len(v.Ref)))
		if err != nil {
			rv.outside++
			continue
		}
		if !strings.EqualFold(string(got), v.Ref) {
			if rv.mismatched == 0 {
				rv.firstMismatch = fmt.Sprintf("%s:%d REF %s, reference %s", v.Chromosome, v.Position, v.Ref, got)
			}
			rv.mismatched++
			continue
		}
		v.Ref, v.Alt = strings.ToUpper(v.Ref), strings.ToUpper(v.Alt)
		k := chromKey(v.Chromosome)
		groups[k] = append(groups[k], v)
	}
	if rv.mismatched > 0 && e.config.StrictRef {
		return nil, fmt.Errorf("strict_ref: %d variants have a REF that does not match the reference %s (first: %s)",
			rv.mismatched, e.config.ReferencePath, rv.firstMismatch)
	}
	for k, vs := range groups {
		rv.byChrom[k] = newPosVariants(vs)
	}
	return rv, nil
}

// report writes the one-line warnings for variants that were not used.
func (rv *refVariants) report(e *Engine) {
	if rv.mismatched > 0 {
		fmt.Fprintf(e.logw, "   [warn] %d variants skipped: their REF does not match the reference (first: %s); set strict_ref / --strict-ref to make this an error\n",
			rv.mismatched, rv.firstMismatch)
	}
	if rv.outside > 0 {
		fmt.Fprintf(e.logw, "   [warn] %d variants lie outside the sampled reference and are never applied\n", rv.outside)
	}
	if rv.blind > 0 {
		fmt.Fprintf(e.logw, "   [warn] %d variants without CHROM/POS are ignored with a reference\n", rv.blind)
	}
}

// origin is where the last sampled fragment came from.
type origin struct {
	chrom  string
	start  int64
	minus  bool
	length int // fragment (or single-end read) length drawn from the sampler
}

// applyReferenceVariants applies the positioned variants that fall in a
// sampled fragment. seq is the fragment as R1 reads it (minus: the reverse
// complement of the reference at [start, start+n)). With trim, the
// fragment is extended to the right by the slack its deletions need
// (from the reference) and trimmed back to n afterwards. Allele draws come
// from r, the stream the fragment was drawn from.
func (p *readPipeline) applyReferenceVariants(r *rand.Rand, pv *posVariants, o origin, seq []byte, trim bool) []byte {
	n := len(seq)
	hi := o.start + int64(n)
	if trim {
		hi += int64(pv.maxSlack)
	}
	local, spans := pv.window(o.start, hi)
	if len(local) == 0 {
		return seq
	}
	fwd := seq
	if o.minus {
		fwd = seqread.ReverseComplement(seq)
	}
	plain := fwd
	if trim {
		if slack := seqread.WindowSlack(o.start, n, spans); slack > 0 {
			// At the end of the loaded sequence, or into N, the fragment
			// cannot be extended; a deletion then leaves it short and the
			// fragment is used without variants (below).
			if ext, err := p.e.ref.Bases(o.chrom, o.start+int64(n), o.start+int64(n+slack)); err == nil && allACGT(ext) {
				fwd = append(append(make([]byte, 0, n+slack), fwd...), ext...)
			}
		}
	}
	out, _, skipped := injectAtPosition(r, fwd, nil, o.start, local)
	if skipped > 0 {
		p.e.overlapSkips.Add(int64(skipped))
	}
	if trim {
		if len(out) < n {
			out = plain
		}
		out = out[:n]
	}
	if o.minus {
		return seqread.ReverseComplement(out)
	}
	return out
}

// allACGT reports whether b holds only A, C, G and T.
func allACGT(b []byte) bool {
	for _, c := range b {
		switch c {
		case 'A', 'C', 'G', 'T':
		default:
			return false
		}
	}
	return true
}

// truthHeader is the first line of a truth table.
const truthHeader = "#read_name\tchrom\tstart\tstrand\tfragment_length\n"

// writeTruth writes one truth-table line for the read (pair) rec: its name,
// the fragment's contig, 0-based leftmost reference position, strand (as R1
// reads it) and length before variants.
func writeTruth(w *bufio.Writer, rec *seqread.Record, o origin) {
	strand := byte('+')
	if o.minus {
		strand = '-'
	}
	w.WriteString(rec.Name) // #nosec G104 -- bufio errors are sticky and surface at Flush
	w.WriteByte('\t')
	w.WriteString(o.chrom)
	w.WriteByte('\t')
	w.WriteString(strconv.FormatInt(o.start, 10))
	w.WriteByte('\t')
	w.WriteByte(strand)
	w.WriteByte('\t')
	w.WriteString(strconv.Itoa(o.length))
	w.WriteByte('\n')
}

// truthPath returns the truth table written next to the FASTQ file fastq.
func truthPath(fastq string) string {
	name := strings.TrimSuffix(fastq, ".gz")
	for _, ext := range []string{".fastq", ".fq"} {
		if strings.HasSuffix(name, ext) {
			name = strings.TrimSuffix(name, ext)
			break
		}
	}
	return name + ".truth.tsv"
}
