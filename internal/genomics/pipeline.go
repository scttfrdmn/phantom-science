package genomics

// The shared per-read pipeline used by file output (Generate,
// GenerateFASTQFile) and streamed output (GenerateToWriter).
//
// Per output stream the hooks are instantiated once (lengthSampler,
// qualityModel, errorModel, headerFormatter; see engine_readlen.go,
// engine_errors.go and engine_output.go). Then, per read (pair), in this
// draw order:
//
// Single-end, legacy qualities (ONT, and Illumina quality_model linear):
//
//	shared: L → (region) read start
//	r:      L+slack bases → qualities (unless HMM gave them) → variants →
//	        trim to L → errors → header
//
// Single-end, binned Illumina qualities:
//
//	shared: L → (region) read start
//	r:      L+slack bases → variants → trim to L → qualities → errors → header
//
// Paired-end (Illumina): each pair is the two ends of one fragment.
//
//	shared: L → insert size F (≥ L) → (region) fragment start →
//	        F+slack fragment bases → variants, applied once → trim to F
//	per mate, from r: R1 = frag[:L], R2 = revcomp(frag[F−L:]) →
//	        qualities (of that mate) → errors → header
//
// Each stream draws only from its own generator (see rng.go), so output is a
// pure function of (seed, configuration, stream ID). The "shared" generator
// is the stream itself when one stream holds both mates (interleaved
// output) or the read is single-end, else the batch's pair stream, so the R1
// and R2 files of a batch regenerate the same fragments independently. For
// paired file output the HMM model (when configured) is seeded from the pair
// stream for the same reason. Fragment HMM qualities are not used: each
// mate's qualities come from the quality model.
//
// Output-size targets, flow cells and parallel file generation (#16) are at
// the end of this file (see generateRun).

import (
	"bufio"
	"context"
	"errors"
	"fmt"
	"io"
	"math/rand/v2"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync/atomic"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/illumina"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
	"github.com/scttfrdmn/phantom-science/internal/security"
	"github.com/scttfrdmn/phantom-science/internal/storage/s3sink"
)

// ErrNotImplemented marks configuration values for planned features that the
// engine does not support yet. Test with errors.Is.
var ErrNotImplemented = errors.New("not implemented yet")

// errAborted stops speculative generation that the run no longer needs.
var errAborted = errors.New("generation aborted")

// checkScale validates the output-size targets, flow cells and workers.
func (e *Engine) checkScale() error {
	c := e.config
	switch {
	case c.TotalBases < 0:
		return fmt.Errorf("invalid total_bases %d", c.TotalBases)
	case c.TotalReads < 0:
		return fmt.Errorf("invalid total_reads %d", c.TotalReads)
	case c.TotalBases > 0 && c.TotalReads > 0:
		return fmt.Errorf("set total_bases or total_reads, not both")
	case c.FlowCells < 0:
		return fmt.Errorf("invalid flow_cells %d", c.FlowCells)
	case c.Workers < 0:
		return fmt.Errorf("invalid workers %d", c.Workers)
	}
	if limit := c.MaxOutputBases(); limit > 0 && c.TotalBases > limit {
		return fmt.Errorf("total_bases %d exceeds the max output of %d bases per flow cell (max_output_gb / max_output_tb)", c.TotalBases, limit)
	}
	if c.FlowCells > 1 && e.isONT() {
		for name, v := range map[string]string{"ont.run_id": c.ONT.RunID, "ont.flow_cell_id": c.ONT.FlowCellID, "flowcell_id": c.FlowCellID} {
			if v != "" {
				return fmt.Errorf("%s %q cannot be shared by %d flow cells: leave it unset to derive one per flow cell from the seed", name, v, c.FlowCells)
			}
		}
	}
	return nil
}

// checkExtensions validates every configured extension before any output is
// written, so unsupported settings fail fast with a clear error.
func (e *Engine) checkExtensions() error {
	if e.seedErr != nil {
		return e.seedErr
	}
	if err := e.checkScale(); err != nil {
		return err
	}
	if err := e.checkOutput(); err != nil {
		return err
	}
	_, err := e.newPipeline(GenerateOpts{}, nil)
	return err
}

// readPipeline holds one output stream's hooks and prepared variants.
type readPipeline struct {
	e            *Engine
	sampleLength func(r *rand.Rand) int
	fixedLength  bool
	qualities    func(r *rand.Rand, n, mate int) []byte
	applyErrors  func(r *rand.Rand, seq, qual []byte) ([]byte, []byte)
	header       func(r *rand.Rand, rec *seqread.Record) (string, string)
	observer     ReadObserver
	maxPhred     byte
	// legacy selects the single-end flow that draws qualities before variant
	// injection (ONT, linear Illumina qualities); binned Illumina qualities
	// are drawn per read after it.
	legacy bool
	// inserts samples paired-end fragment lengths (nil for single-end).
	inserts *illumina.InsertSampler
	// fragments, when set, supplies fragment bases (and their reference
	// position) in place of the HMM model or uniform bases.
	fragments FragmentSampler
	// tags is appended to every header comment: the read-ID fingerprint
	// tags ("\tXC:Z:<customer>\tXD:Z:<dataset>"), or "".
	tags string

	// model generates bases (nil: uniform random bases). It is the
	// stream's own instance for file output.
	model authenticity.SequenceGenerator
	// admit, when set, is asked before each read (pair) is written with its
	// total bases; false stops the batch before that read.
	admit func(bases int64) bool
	// tally, when set, is told the total bases of each written read (pair).
	tally func(bases int)
	// abort, when set and true, stops the batch with errAborted.
	abort *atomic.Bool

	sampleID string
	region   *GenomicRegion
	// positioned holds the variants applied by coordinate (region mode),
	// sorted by descending position; spans are their reference spans.
	positioned []VariantRecord
	spans      []seqread.Span
	maxSlack   int // total deletion length of all positioned variants
	// blind holds position-blind variants; blindSlack is their total
	// deletion length.
	blind      []VariantRecord
	blindSlack int

	pending []seqread.Record // the mates of the read being written
}

// FragmentSampler supplies the bases of one fragment (or single-end read)
// of n bases: seq is the fragment as read 1 reads it, chrom and start (0-based)
// its reference position, and minus reports that seq is the reverse
// complement of the reference there. It must draw only from r. It is the
// hook for sampling reads from a reference genome; nil (the default) uses
// the HMM model or uniform bases.
type FragmentSampler func(r *rand.Rand, n int) (chrom string, start int64, minus bool, seq []byte, err error)

// newPipeline instantiates the hooks for one output stream.
func (e *Engine) newPipeline(opts GenerateOpts, obs ReadObserver) (*readPipeline, error) {
	p := &readPipeline{e: e, observer: obs, sampleID: opts.SampleID, region: opts.RefRegion, maxPhred: e.maxPhred(),
		fixedLength: e.fixedLengthReads(), legacy: !e.binnedQualities(), fragments: e.fragments, tags: e.readTags()}
	var err error
	if p.sampleLength, err = e.lengthSampler(); err != nil {
		return nil, err
	}
	if p.qualities, err = e.qualityModel(); err != nil {
		return nil, err
	}
	if p.applyErrors, err = e.errorModel(); err != nil {
		return nil, err
	}
	if p.header, err = e.headerFormatter(); err != nil {
		return nil, err
	}
	if e.config.PairedEnd {
		if p.inserts, err = illumina.NewInsertSampler(e.illuminaSpec().InsertSize); err != nil {
			return nil, err
		}
	}
	for _, v := range opts.Variants {
		if len(v.Ref) == 0 || len(v.Alt) == 0 {
			continue
		}
		if p.region == nil {
			p.blind = append(p.blind, v)
			p.blindSlack += max(0, len(v.Ref)-len(v.Alt))
			continue
		}
		if v.Chromosome == "" || v.Position <= 0 || !sameChromosome(v.Chromosome, p.region.Chromosome) {
			continue // region mode applies only positioned variants on the region's chromosome
		}
		p.positioned = append(p.positioned, v)
	}
	sort.SliceStable(p.positioned, func(i, j int) bool { return p.positioned[i].Position > p.positioned[j].Position })
	for _, v := range p.positioned {
		p.spans = append(p.spans, seqread.Span{Pos: v.Position - 1, RefLen: len(v.Ref), AltLen: len(v.Alt)})
		p.maxSlack += max(0, len(v.Ref)-len(v.Alt))
	}
	return p, nil
}

// window returns the positioned variants (and their spans) whose 0-based
// start lies in [lo, hi), still in descending order.
func (p *readPipeline) window(lo, hi int64) ([]VariantRecord, []seqread.Span) {
	a := sort.Search(len(p.positioned), func(i int) bool { return p.positioned[i].Position-1 < hi })
	b := sort.Search(len(p.positioned), func(i int) bool { return p.positioned[i].Position-1 < lo })
	return p.positioned[a:b], p.spans[a:b]
}

// sameChromosome compares chromosome names, ignoring a "chr" prefix. An
// empty region chromosome matches every variant.
func sameChromosome(a, region string) bool {
	if region == "" {
		return true
	}
	return strings.TrimPrefix(a, "chr") == strings.TrimPrefix(region, "chr")
}

// batchSpec says what one writeBatch call writes.
type batchSpec struct {
	batch int
	reads int   // reads (pairs) to write
	mates []int // {1} or {2} for one file of a batch, {1, 2} for interleaved pairs
	file  string
	// shared draws the per-read values both mates share: the read length,
	// and for paired-end reads the fragment (insert size, start, bases and
	// variant alleles); for single-end reads, in region mode, the read start.
	shared *rand.Rand
}

// writeBatch writes one batch of reads to bw, drawing from r (and s.shared).
// file is the output file name recorded on each read ("" for streams). It
// reports whether p.admit stopped the batch early.
func (p *readPipeline) writeBatch(bw *bufio.Writer, r *rand.Rand, s batchSpec, serial *int64) (stopped bool, err error) {
	rec := seqread.Record{Batch: s.batch, ReadGroup: p.sampleID, File: s.file}
	if p.region != nil {
		rec.Positioned = true
		rec.Chromosome = p.region.Chromosome
	}
	for readNum := 1; readNum <= s.reads; readNum++ {
		if p.abort != nil && p.abort.Load() {
			return false, errAborted
		}
		length := p.sampleLength(s.shared)
		// frag is the pair's fragment, or for single-end reads from a
		// FragmentSampler the read itself; nil for other single-end reads.
		var frag []byte
		var readStart int64
		switch {
		case p.inserts != nil:
			frag, readStart, err = p.fragment(s.shared, length, &rec)
		case p.fragments != nil:
			frag, readStart, err = p.sampled(s.shared, length, &rec)
		case p.region != nil:
			readStart, err = p.start(s.shared, length)
		}
		if err != nil {
			return false, err
		}
		rec.RefStart = readStart
		p.pending = p.pending[:0]
		var bases int
		for i, mate := range s.mates {
			var seq, qual []byte
			if frag != nil {
				seq, qual = p.mate(r, frag, length, mate)
			} else if seq, qual, err = p.read(r, length, readStart); err != nil {
				return false, err
			}
			rec.Index, rec.Mate, rec.Serial = readNum, mate, *serial+int64(i)
			rec.Seq, rec.Qual = seq, qual
			rec.Name, rec.Comment = p.header(r, &rec)
			if p.sampleID != "" {
				rec.Comment = joinComment(rec.Comment, "RG:Z:"+p.sampleID)
			}
			if p.tags != "" {
				rec.Comment = joinComment(rec.Comment, p.tags)
			}
			p.pending = append(p.pending, rec)
			bases += len(seq)
		}
		if p.admit != nil && !p.admit(int64(bases)) {
			return true, nil
		}
		for i := range p.pending {
			if err := p.writeRecord(bw, &p.pending[i]); err != nil {
				return false, err
			}
			if p.observer != nil {
				if err := p.observer.ObserveRead(&p.pending[i]); err != nil {
					return false, fmt.Errorf("observing read: %w", err)
				}
			}
			*serial++
		}
		if p.tally != nil {
			p.tally(bases)
		}
	}
	return false, nil
}

// joinComment appends tab-separated fields to a header comment.
func joinComment(comment, fields string) string {
	if comment == "" {
		return fields
	}
	return comment + "\t" + fields
}

// start draws the 0-based start of a read (or fragment) of n bases in the
// region.
func (p *readPipeline) start(r *rand.Rand, n int) (int64, error) {
	span := p.region.End - p.region.Start - int64(n)
	if span <= 0 {
		return 0, fmt.Errorf("region [%s:%d-%d) too small for a read or fragment of %d bases",
			p.region.Chromosome, p.region.Start, p.region.End, n)
	}
	return p.region.Start + r.Int64N(span), nil
}

// fragment draws one paired-end fragment from r (the pair's shared stream):
// its insert size F >= length, its start in region mode, F bases with the
// variants applied once, trimmed back to F. With a FragmentSampler the bases
// and position come from it instead, and rec's position is set.
func (p *readPipeline) fragment(r *rand.Rand, length int, rec *seqread.Record) ([]byte, int64, error) {
	n := p.inserts.Sample(r, length)
	if p.fragments != nil {
		return p.sampled(r, n, rec)
	}
	var start int64
	if p.region != nil {
		var err error
		if start, err = p.start(r, n); err != nil {
			return nil, 0, err
		}
	}
	seq, _, err := p.bases(r, n, start, true, nil)
	return seq, start, err
}

// sampled returns n bases from the FragmentSampler and records their
// position on rec.
func (p *readPipeline) sampled(r *rand.Rand, n int, rec *seqread.Record) ([]byte, int64, error) {
	chrom, start, _, seq, err := p.fragments(r, n)
	if err != nil {
		return nil, 0, fmt.Errorf("sampling a fragment: %w", err)
	}
	if len(seq) != n {
		return nil, 0, fmt.Errorf("fragment sampler returned %d bases, want %d", len(seq), n)
	}
	rec.Positioned, rec.Chromosome, rec.RefStart = true, chrom, start
	return seq, start, nil
}

// mate returns the bases and qualities of one mate of frag: R1 reads the
// first length bases, R2 the reverse complement of the last length bases.
// Qualities and errors come from r, the mate's stream.
func (p *readPipeline) mate(r *rand.Rand, frag []byte, length, mate int) ([]byte, []byte) {
	var seq []byte
	if mate == 2 {
		seq = seqread.ReverseComplement(frag[len(frag)-length:])
	} else {
		seq = append([]byte(nil), frag[:length]...)
	}
	qual := p.qualities(r, length, mate)
	if p.applyErrors != nil {
		seq, qual = p.applyErrors(r, seq, qual)
	}
	return seq, qual
}

// bases generates the bases of a read or fragment of length bases starting
// at the 0-based reference position start, from r: length plus slack bases
// (enough to absorb the deletions that can apply), then the variants, then,
// when trim is set, trimmed back to length. qualFn, when set, gives the
// qualities of bases without HMM qualities; they are generated before the
// variants and kept aligned with the bases (the legacy order). Without
// qualFn the returned qualities are nil.
func (p *readPipeline) bases(r *rand.Rand, length int, start int64, trim bool, qualFn func(n int) []byte) ([]byte, []byte, error) {
	e := p.e
	n := length
	var local []VariantRecord // positioned variants that can touch this read
	if len(p.positioned) > 0 {
		hi := start + int64(length)
		if trim {
			hi += int64(p.maxSlack)
		}
		var spans []seqread.Span
		local, spans = p.window(start, hi)
		if trim {
			n += seqread.WindowSlack(start, length, spans)
		}
	} else if trim && p.region == nil {
		n += p.blindSlack
	}
	gen := func(n int) ([]byte, []byte, error) {
		seq, qual, err := generateFrom(p.model, r, n)
		if err != nil {
			return nil, nil, fmt.Errorf("failed to generate sequence: %w", err)
		}
		switch {
		case qualFn == nil:
			qual = nil
		case len(qual) != len(seq):
			qual = qualFn(len(seq))
		}
		return seq, qual, nil
	}
	seq, qual, err := gen(n)
	if err != nil {
		return nil, nil, err
	}
	switch {
	case len(local) > 0:
		var skipped int
		seq, qual, skipped = injectAtPosition(r, seq, qual, start, local)
		if skipped > 0 {
			e.overlapSkips.Add(int64(skipped))
		}
	case len(p.blind) > 0:
		seq, qual = injectVariants(r, seq, qual, p.blind)
	}
	if trim {
		if len(seq) < length { // cannot happen with correct slack; never pad with N
			extra, extraQ, err := gen(length - len(seq))
			if err != nil {
				return nil, nil, err
			}
			seq = append(seq, extra...)
			if qual != nil {
				qual = append(qual, extraQ...)
			}
		}
		seq = seq[:length]
		if qual != nil {
			qual = qual[:length]
		}
	}
	return seq, qual, nil
}

// read produces the final bases and Phred qualities of one single-end read,
// drawing everything from r.
func (p *readPipeline) read(r *rand.Rand, length int, readStart int64) ([]byte, []byte, error) {
	var seq, qual []byte
	var err error
	if p.legacy {
		seq, qual, err = p.bases(r, length, readStart, p.fixedLength, func(n int) []byte { return p.qualities(r, n, 1) })
	} else if seq, _, err = p.bases(r, length, readStart, p.fixedLength, nil); err == nil {
		qual = p.qualities(r, len(seq), 1)
	}
	if err != nil {
		return nil, nil, err
	}
	if p.applyErrors != nil {
		seq, qual = p.applyErrors(r, seq, qual)
	}
	return seq, qual, nil
}

// writeRecord writes rec as a FASTQ record.
func (p *readPipeline) writeRecord(bw *bufio.Writer, rec *seqread.Record) error {
	bw.WriteString(rec.Header()) // #nosec G104 -- bufio errors are sticky and surface below
	bw.WriteByte('\n')
	bw.Write(rec.Seq)
	bw.WriteString("\n+\n")
	for _, q := range rec.Qual {
		bw.WriteByte('!' + min(q, p.maxPhred))
	}
	if err := bw.WriteByte('\n'); err != nil {
		return fmt.Errorf("writing FASTQ record: %w", err)
	}
	return nil
}

// injectVariants applies position-blind variants in order: each variant is
// drawn ALT with probability Freq (r.Float64() < Freq, so Freq >= 1 always
// and Freq <= 0 never), then its first REF occurrence is replaced by ALT.
// qual (may be nil) stays aligned with seq; the read length changes by the
// indel size, and the pipeline trims fixed-length reads afterwards.
func injectVariants(r *rand.Rand, seq, qual []byte, variants []VariantRecord) ([]byte, []byte) {
	for _, v := range variants {
		if len(v.Ref) == 0 || len(v.Alt) == 0 {
			continue
		}
		if !(r.Float64() < v.Freq) {
			continue
		}
		idx := strings.Index(string(seq), v.Ref)
		if idx < 0 {
			continue
		}
		seq, qual = seqread.Splice(seq, qual, idx, len(v.Ref), []byte(v.Alt))
	}
	return seq, qual
}

// injectAtPosition writes each positioned variant's allele (ALT with
// probability Freq, else REF) at its offset in a read whose first base is the
// 0-based reference position readStart. variants must be sorted by
// descending position: applying right to left keeps the offsets of variants
// further left valid after an indel. Variants that do not fit entirely in
// the read are ignored; a variant whose reference span overlaps one already
// applied is skipped and counted. qual (may be nil) stays aligned with seq.
func injectAtPosition(r *rand.Rand, seq, qual []byte, readStart int64, variants []VariantRecord) ([]byte, []byte, int) {
	window := int64(len(seq))
	limit := window // reference offset where the applied variant furthest left starts
	skipped := 0
	for _, v := range variants {
		if v.Position <= 0 || len(v.Ref) == 0 || len(v.Alt) == 0 {
			continue
		}
		off := (v.Position - 1) - readStart
		end := off + int64(len(v.Ref))
		if off < 0 || end > window {
			continue
		}
		if end > limit {
			skipped++
			continue
		}
		allele := v.Ref
		if r.Float64() < v.Freq {
			allele = v.Alt
		}
		seq, qual = seqread.Splice(seq, qual, int(off), len(v.Ref), []byte(allele))
		limit = off
	}
	return seq, qual, skipped
}

// generateFrom returns n bases and their Phred qualities from model, or
// uniform random bases drawn from r (with nil qualities) when model is nil
// or fails.
func generateFrom(model authenticity.SequenceGenerator, r *rand.Rand, n int) ([]byte, []byte, error) {
	if model != nil {
		result, err := model.GenerateSequence(n)
		if err == nil {
			return result.Sequence, result.Qualities, nil
		}
		// fall through to uniform bases on per-call HMM failure
	}
	return randomBases(r, n), nil, nil
}

// runStart is the run's seed-derived acquisition start (the ONT start_time
// origin), recorded in the manifest instead of a wall-clock time.
func (e *Engine) runStart() time.Time {
	return nanopore.RunInfo{}.WithDefaults(e.seed).Start
}

// Output plan, flow cells and parallel file generation (#16).
//
// File plan. Generate writes units: unit n is batch n, one FASTQ file (two
// for paired-end, R1 and R2) of up to ReadsPerBatch reads (read pairs). A
// unit's bytes are a pure function of (flow-cell seed, configuration, n,
// number of reads), and a unit cut to k reads is exactly the first k reads
// of the full unit. Per flow cell the plan is:
//
//   - default: TotalBatches full units;
//   - TotalReads: ceil(TotalReads / ReadsPerBatch) units, the last partial;
//     fixed up front, exact;
//   - TotalBases: units until the target is met. This is streaming: the
//     number of files depends on the generated reads and is not known up
//     front, but it is still a pure function of the seed.
//
// The budget (the TotalBases target and the MaxOutputBases cap) is applied
// to the actual bases of each read (pair), after error-model indels, in
// output order: a read is written while the total is below the target and
// writing it does not exceed the cap. Output therefore stops within one read
// (pair) length of the target, and never exceeds the cap.
//
// Parallelism. Up to Workers units are generated at once, each by its own
// goroutine into a spool directory, with its own RNG streams and its own HMM
// model instance (loaded from the model path and seeded from the unit's
// stream), so workers share no mutable state. One coordinator commits units
// strictly in unit order (moving them into place, or uploading them) and
// applies the budget: when the budget stops inside unit n, unit n is
// regenerated with the admitted number of reads and the speculative units
// after it are discarded. ONT sequencing-summary rows are written per unit
// to a part file and concatenated in unit order. Files, summaries and the
// manifest are therefore the same for any number of workers.

// plan is the output plan of one flow cell.
type plan struct {
	units     int   // units to write; 0: until the budget stops the output
	perUnit   int   // reads (pairs) per unit
	lastReads int   // reads in unit `units` (TotalReads); 0: perUnit
	bases     int64 // TotalBases target; 0: none
	cap       int64 // max output bases; 0: none
}

// outputPlan returns the configured output plan (per flow cell).
func (e *Engine) outputPlan() plan {
	c := e.config
	p := plan{perUnit: c.ReadsPerBatch, bases: c.TotalBases, cap: c.MaxOutputBases()}
	switch per := int64(max(c.ReadsPerBatch, 1)); {
	case c.TotalBases > 0:
	case c.TotalReads > 0:
		p.units = int((c.TotalReads + per - 1) / per)
		p.lastReads = int(c.TotalReads - int64(p.units-1)*per)
	default:
		p.units = c.TotalBatches
	}
	return p
}

// has reports whether unit is part of the plan.
func (p plan) has(unit int) bool { return p.units == 0 || unit <= p.units }

// reads returns the number of reads (pairs) planned for unit.
func (p plan) reads(unit int) int {
	if p.lastReads > 0 && unit == p.units {
		return p.lastReads
	}
	return p.perUnit
}

// multiFile reports whether the run may write more than one file per mate.
func (e *Engine) multiFile() bool {
	c := e.config
	switch {
	case c.TotalBases > 0:
		return true
	case c.TotalReads > 0:
		return c.TotalReads > int64(c.ReadsPerBatch)
	}
	return c.TotalBatches > 1
}

// budget applies a plan's bases target and cap, read by read in output order.
type budget struct {
	p            plan
	reads, bases int64 // admitted reads (pairs) and bases
	capped       bool  // the cap stopped the output
}

// admit reports whether the next read (pair), of n bases, is written; false
// ends the output before it.
func (b *budget) admit(n int64) bool {
	if b.p.bases > 0 && b.bases >= b.p.bases {
		return false
	}
	if b.p.cap > 0 && b.bases+n > b.p.cap {
		b.capped = true
		return false
	}
	b.reads++
	b.bases += n
	return true
}

// mates returns the mate numbers of one batch.
func (e *Engine) mates() []int {
	if e.config.PairedEnd {
		return []int{1, 2}
	}
	return []int{1}
}

// destination receives a run's finished files.
type destination interface {
	// spool is the local directory files are written and finalized in.
	spool() string
	// commit moves the finished local file to rel (slash-separated, relative
	// to the output root) and returns its final location and S3 key.
	commit(local, rel string, tags map[string]string) (final, key string, err error)
	uri() string // S3 URI of the output root; "" for local output
	close() error
}

// localDest writes to an output directory through a spool inside it, so
// commits are renames on one file system.
type localDest struct{ root, dir string }

func (d *localDest) spool() string { return d.dir }
func (d *localDest) uri() string   { return "" }
func (d *localDest) close() error  { return os.RemoveAll(d.dir) }

func (d *localDest) commit(local, rel string, _ map[string]string) (string, string, error) {
	dst := filepath.Join(d.root, filepath.FromSlash(rel))
	if err := security.SecureMkdirAll(filepath.Dir(dst), d.root); err != nil {
		return "", "", err
	}
	if err := os.Rename(local, dst); err != nil {
		return "", "", fmt.Errorf("moving %s into place: %w", rel, err)
	}
	return dst, "", nil
}

// s3Dest uploads files from a temporary spool directory.
type s3Dest struct {
	sink *s3sink.Sink
	dir  string
}

func (d *s3Dest) spool() string { return d.dir }
func (d *s3Dest) uri() string   { return d.sink.Location().URI() }
func (d *s3Dest) close() error  { return os.RemoveAll(d.dir) }

func (d *s3Dest) commit(local, rel string, tags map[string]string) (string, string, error) {
	key, err := d.sink.Put(context.Background(), local, rel, tags)
	if err != nil {
		return "", "", err
	}
	return s3sink.Scheme + d.sink.Location().Bucket + "/" + key, key, nil
}

// openDestination prepares the configured output: a local directory, or an
// S3 prefix when OutputDir is an s3:// URI.
func (e *Engine) openDestination() (destination, error) {
	out := e.config.OutputDir
	if s3sink.IsURI(out) {
		loc, err := s3sink.ParseURI(out)
		if err != nil {
			return nil, err
		}
		up := e.uploader
		if up == nil {
			if up, err = s3sink.NewDefaultUploader(); err != nil {
				return nil, err
			}
		}
		dir, err := os.MkdirTemp("", "phantom-spool-")
		if err != nil {
			return nil, fmt.Errorf("creating spool directory: %w", err)
		}
		return &s3Dest{sink: s3sink.New(loc, up), dir: dir}, nil
	}
	if err := security.SecureMkdirAll(out, ""); err != nil {
		return nil, fmt.Errorf("failed to create output directory: %w", err)
	}
	dir, err := os.MkdirTemp(out, ".phantom-spool-")
	if err != nil {
		return nil, fmt.Errorf("creating spool directory: %w", err)
	}
	return &localDest{root: out, dir: dir}, nil
}

// generateRun writes every flow cell of the run and the manifest, which is
// committed last.
func (e *Engine) generateRun() (err error) {
	dest, err := e.openDestination()
	if err != nil {
		return err
	}
	defer func() { err = errors.Join(err, dest.close()) }()

	m := manifest.New(e.toolVersion, e.seed, e.runStart(), e.config)
	m.Config.Seed = e.seed
	m.Destination = dest.uri()
	m.Encrypted = e.encrypted && e.qualityFingerprintingEnabled()
	n := max(e.config.FlowCells, 1)
	seen := map[string]bool{}
	for k := 1; k <= n; k++ {
		fc, rel, err := e.flowCellEngine(k, n, dest.spool(), seen)
		if err != nil {
			return err
		}
		if err := e.runFlowCell(fc, k, rel, dest, m); err != nil {
			return err
		}
	}
	data, err := m.Marshal()
	if err != nil {
		return err
	}
	local := filepath.Join(dest.spool(), manifest.FileName)
	if err := security.SecureWriteFile(local, data, dest.spool()); err != nil {
		return fmt.Errorf("writing manifest: %w", err)
	}
	if _, _, err := dest.commit(local, manifest.FileName, nil); err != nil {
		return err
	}
	e.manifest = m
	return nil
}

// flowCellSeedAttempts bounds the search for distinct flow-cell IDs.
const flowCellSeedAttempts = 64

// flowCellEngine returns the engine generating flow cell k of n into its
// spool subdirectory, and that subdirectory relative to the output root
// ("" for a single flow cell). With n > 1 each flow cell is an independent
// run with a seed derived from the run seed (and, for ONT, distinct derived
// run and flow-cell IDs; explicit Illumina run_id / flowcell_id values get a
// "-<k>" suffix).
func (e *Engine) flowCellEngine(k, n int, spool string, seen map[string]bool) (*Engine, string, error) {
	c := e.config
	c.FlowCells, c.Workers = 0, 0
	if n > 1 {
		if c.RunID != "" {
			c.RunID += "-" + strconv.Itoa(k)
		}
		if c.FlowCellID != "" {
			c.FlowCellID += "-" + strconv.Itoa(k)
		}
	}
	fc := &Engine{
		config: c, parent: e, logw: e.logw, seed: e.seed,
		fingerprintingEngine: e.fingerprintingEngine, keyErr: e.keyErr, encrypted: e.encrypted,
		hmmModel: e.hmmModel, modelPath: e.modelPath, toolVersion: e.toolVersion,
		customerID: e.customerID, fragments: e.fragments,
	}
	rel := ""
	if n > 1 {
		for attempt := 0; ; attempt++ {
			if attempt == flowCellSeedAttempts {
				return nil, "", fmt.Errorf("flow cell %d: no distinct flow-cell ID derived from seed %d", k, e.seed)
			}
			fc.seed = flowCellSeed(e.seed, k, attempt)
			rel = fc.flowCellDir(k)
			id := rel
			if fc.isONT() {
				id += "/" + fc.ontRunInfo().RunID
			}
			if !seen[rel] && !seen[id] {
				seen[rel], seen[id] = true, true
				break
			}
		}
		if rel == "" || rel == "." || rel == ".." || filepath.Base(rel) != rel {
			return nil, "", fmt.Errorf("flow cell %d: %q is not a valid directory name", k, rel)
		}
	}
	fc.config.Seed = fc.seed
	fc.config.OutputDir = filepath.Join(spool, rel)
	if err := security.SecureMkdirAll(fc.config.OutputDir, spool); err != nil {
		return nil, "", err
	}
	return fc, rel, nil
}

// flowCellDir names flow cell k's subdirectory: the ONT flow-cell ID, the
// Illumina flowcell_id, or flowcell_<k>.
func (e *Engine) flowCellDir(k int) string {
	switch {
	case e.isONT():
		return e.ontRunInfo().FlowCellID
	case e.config.FlowCellID != "":
		return e.config.FlowCellID
	}
	return "flowcell_" + strconv.Itoa(k)
}

// unitFile is one finished file of a unit, in the spool.
type unitFile struct {
	path         string
	bytes        int64
	sha256       string
	reads, bases int64
}

// unitResult is the outcome of generating one unit.
type unitResult struct {
	unit    int
	files   []unitFile
	pairs   []int64 // total bases of each read (pair)
	part    string  // sequencing-summary part file; "" for none
	elapsed time.Duration
	err     error
}

// generateUnit writes unit (batch) with reads reads (pairs) into the spool
// and hashes the finished files.
func (e *Engine) generateUnit(unit, reads int, abort *atomic.Bool) unitResult {
	start := time.Now()
	res := unitResult{unit: unit}
	var obs ReadObserver
	if e.isONT() {
		res.part = filepath.Join(e.config.OutputDir, fmt.Sprintf(".summary-%06d.part", unit))
		obs = e.summaryObserver(res.part)
	}
	for i, mate := range e.mates() {
		lens := make([]int64, 0, reads)
		final, err := e.generateFASTQFile(fileJob{
			path: e.outputPath(unit, mate), batch: unit, mate: mate, reads: reads,
			obs: obs, abort: abort, tally: func(n int) { lens = append(lens, int64(n)) },
		})
		if err == nil {
			f := unitFile{path: final, reads: int64(len(lens))}
			for _, n := range lens {
				f.bases += n
			}
			f.bytes, f.sha256, err = manifest.HashFile(final)
			res.files = append(res.files, f)
		}
		if err != nil {
			closeObserver(obs) // #nosec G104 -- the generation error is returned
			if e.config.PairedEnd {
				err = fmt.Errorf("failed to generate R%d file: %w", mate, err)
			} else {
				err = fmt.Errorf("failed to generate FASTQ file: %w", err)
			}
			res.err = err
			return res
		}
		if i == 0 {
			res.pairs = lens
		} else {
			for j := range min(len(lens), len(res.pairs)) {
				res.pairs[j] += lens[j]
			}
		}
	}
	res.err = closeObserver(obs)
	res.elapsed = time.Since(start)
	return res
}

// discardUnit removes a unit's spooled files.
func discardUnit(r unitResult) {
	for _, f := range r.files {
		os.Remove(f.path) // #nosec G104 -- best effort; the spool is removed at the end
	}
	if r.part != "" {
		os.Remove(r.part) // #nosec G104
	}
}

// runFlowCell generates and commits the units of flow cell k.
func (e *Engine) runFlowCell(fc *Engine, k int, rel string, dest destination, m *manifest.Manifest) error {
	p := fc.outputPlan()
	b := &budget{p: p}
	workers := max(e.config.Workers, 1)
	if fc.hmmModel != nil && fc.modelPath == "" {
		workers = 1 // a model without a path cannot be instantiated per worker
	}
	entry := manifest.FlowCell{Index: k, Seed: fc.seed, Dir: rel, RunID: fc.config.RunID, FlowCellID: fc.config.FlowCellID}
	if fc.isONT() {
		info := fc.ontRunInfo()
		entry.RunID, entry.FlowCellID = info.RunID, info.FlowCellID
	}
	summary := fc.newSummaryMerger()
	defer summary.discard()

	var abort atomic.Bool
	results := make(chan unitResult, workers)
	pending := map[int]unitResult{}
	inflight, next, commit := 0, 1, 1
	done := !p.has(1)
	defer func() {
		abort.Store(true)
		for ; inflight > 0; inflight-- {
			discardUnit(<-results)
		}
		for _, r := range pending {
			discardUnit(r)
		}
	}()
	for {
		for !done && inflight < workers && p.has(next) && next-commit < 2*workers {
			unit, reads := next, p.reads(next)
			go func() { results <- fc.generateUnit(unit, reads, &abort) }()
			next++
			inflight++
		}
		if inflight == 0 {
			break
		}
		r := <-results
		inflight--
		pending[r.unit] = r
		for !done {
			r, ok := pending[commit]
			if !ok {
				break
			}
			delete(pending, commit)
			if r.err != nil {
				discardUnit(r)
				return r.err
			}
			cut := 0
			for cut < len(r.pairs) && b.admit(r.pairs[cut]) {
				cut++
			}
			if cut < len(r.pairs) {
				done = true
				abort.Store(true)
				discardUnit(r)
				if cut == 0 {
					break
				}
				if r = fc.generateUnit(commit, cut, nil); r.err != nil {
					discardUnit(r)
					return r.err
				}
				if len(r.pairs) != cut {
					discardUnit(r)
					return fmt.Errorf("batch %d regenerated with %d reads, want %d", commit, len(r.pairs), cut)
				}
			}
			if err := e.commitUnit(fc, r, k, rel, dest, summary, m, &entry); err != nil {
				return err
			}
			if p.units > 0 {
				fmt.Printf("   Wrote batch %d/%d\n", commit, p.units)
			} else {
				fmt.Printf("   Wrote batch %d (%d bases so far)\n", commit, b.bases)
			}
			commit++
			if !p.has(commit) || (p.bases > 0 && b.bases >= p.bases) {
				done = true
			}
		}
		if done {
			abort.Store(true)
		}
	}
	if entry.Reads == 0 {
		if b.capped {
			return fmt.Errorf("flow cell %d: no reads written: the max output of %d bases is smaller than the first read", k, p.cap)
		}
		return fmt.Errorf("flow cell %d: no reads written", k)
	}
	if b.capped {
		entry.Capped = true
		fmt.Fprintf(e.logw, "   [warn] flow cell %d stopped at the max output of %d bases (max_output_gb / max_output_tb)\n", k, p.cap)
	}
	path, err := summary.finish()
	if err != nil {
		return err
	}
	if path != "" {
		n, sum, err := manifest.HashFile(path)
		if err != nil {
			return err
		}
		relPath := joinRel(rel, filepath.Base(path))
		_, key, err := dest.commit(path, relPath, fc.objectTags(n, 0))
		if err != nil {
			return err
		}
		m.Add(manifest.Output{Path: relPath, Key: key, Role: manifest.RoleSummary, FlowCell: k, Bytes: n, SHA256: sum})
	}
	m.FlowCells = append(m.FlowCells, entry)
	return nil
}

// commitUnit commits a unit's files and summary rows and records them.
func (e *Engine) commitUnit(fc *Engine, r unitResult, k int, rel string, dest destination, summary *summaryMerger, m *manifest.Manifest, entry *manifest.FlowCell) error {
	for _, f := range r.files {
		relPath := joinRel(rel, filepath.Base(f.path))
		final, key, err := dest.commit(f.path, relPath, fc.objectTags(f.bytes, r.elapsed))
		if err != nil {
			return err
		}
		m.Add(manifest.Output{Path: relPath, Key: key, Role: manifest.RoleFASTQ, FlowCell: k, Bytes: f.bytes, SHA256: f.sha256, Reads: f.reads, Bases: f.bases})
		entry.Files++
		entry.Reads += f.reads
		entry.Bases += f.bases
		e.filesMu.Lock()
		e.files = append(e.files, final)
		e.filesMu.Unlock()
	}
	return summary.add(r.part)
}

// joinRel joins an optional subdirectory and a file name with "/".
func joinRel(dir, name string) string {
	if dir == "" {
		return name
	}
	return dir + "/" + name
}

// objectTags are the S3 object tags of an uploaded output: it is
// regenerable from the seed, its size, and the resource time spent making
// it (wall-clock hours of the worker that generated its batch).
func (e *Engine) objectTags(bytes int64, elapsed time.Duration) map[string]string {
	return map[string]string{
		"phantom:regenerable":    "seed",
		"phantom:seed":           strconv.FormatUint(e.seed, 10),
		"phantom:output-bytes":   strconv.FormatInt(bytes, 10),
		"phantom:resource-hours": strconv.FormatFloat(elapsed.Hours(), 'f', 6, 64),
	}
}

// summaryMerger concatenates per-unit sequencing-summary parts, in unit
// order, into the flow cell's summary (one header line).
type summaryMerger struct {
	path string
	f    *os.File
	bw   *bufio.Writer
	dir  string
}

// newSummaryMerger returns the flow cell's merger, or nil when the platform
// writes no summary.
func (e *Engine) newSummaryMerger() *summaryMerger {
	if !e.isONT() {
		return nil
	}
	name := e.config.Output.SummaryFile
	if name == "" {
		name = nanopore.SummaryFileName
	}
	return &summaryMerger{path: filepath.Join(e.config.OutputDir, name), dir: e.config.OutputDir}
}

// add appends the part file (if it exists) and removes it.
func (s *summaryMerger) add(part string) error {
	if s == nil || part == "" {
		return nil
	}
	in, err := os.Open(part) // #nosec G304 -- a spooled part file of this run
	if errors.Is(err, os.ErrNotExist) {
		return nil // the unit wrote no reads
	}
	if err != nil {
		return err
	}
	defer os.Remove(part) // #nosec G104 -- spool cleanup
	defer in.Close()
	br := bufio.NewReader(in)
	if s.f == nil {
		if s.f, err = security.SecureFileCreate(s.path, s.dir); err != nil {
			return fmt.Errorf("creating sequencing summary: %w", err)
		}
		s.bw = bufio.NewWriterSize(s.f, 256*1024)
	} else if _, err := br.ReadSlice('\n'); err != nil && !errors.Is(err, io.EOF) {
		return fmt.Errorf("reading summary part: %w", err) // skip the repeated header
	}
	if _, err := io.Copy(s.bw, br); err != nil {
		return fmt.Errorf("writing sequencing summary: %w", err)
	}
	return nil
}

// finish closes the merged summary and returns its path ("" when no unit
// wrote a summary).
func (s *summaryMerger) finish() (string, error) {
	if s == nil || s.f == nil {
		return "", nil
	}
	err := errors.Join(s.bw.Flush(), s.f.Close())
	s.f = nil
	if err != nil {
		return "", fmt.Errorf("writing sequencing summary: %w", err)
	}
	return s.path, nil
}

// discard closes an unfinished summary.
func (s *summaryMerger) discard() {
	if s != nil && s.f != nil {
		s.f.Close() // #nosec G104 -- abandoned with the spool
		s.f = nil
	}
}
