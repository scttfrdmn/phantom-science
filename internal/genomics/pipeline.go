package genomics

// The shared per-read pipeline used by file output (Generate,
// GenerateFASTQFile) and streamed output (GenerateToWriter).
//
// Per output stream the hooks are instantiated once (lengthSampler,
// qualityModel, errorModel, headerFormatter; see engine_readlen.go,
// engine_errors.go and engine_output.go). Then, per read:
//
//	sample length L → (region) draw read start → per mate:
//	generate L+slack bases → inject variants → trim to L →
//	apply error model → format header → write → notify observer
//
// Each stream draws only from its own generator (see rng.go), in that order,
// so output is a pure function of (seed, configuration, stream ID).

import (
	"bufio"
	"errors"
	"fmt"
	"math/rand/v2"
	"sort"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
)

// ErrNotImplemented marks configuration values for planned features that the
// engine does not support yet. Test with errors.Is.
var ErrNotImplemented = errors.New("not implemented yet")

// checkScale rejects output-size targets and parallelism settings (#16).
func (e *Engine) checkScale() error {
	c := e.config
	switch {
	case c.TotalBases != 0:
		return fmt.Errorf("total_bases output target is %w (#16)", ErrNotImplemented)
	case c.TotalReads != 0:
		return fmt.Errorf("total_reads output target is %w (#16)", ErrNotImplemented)
	case c.FlowCells != 0:
		return fmt.Errorf("flow_cells is %w (#16)", ErrNotImplemented)
	case c.Workers > 1:
		return fmt.Errorf("parallel generation (workers=%d) is %w (#16)", c.Workers, ErrNotImplemented)
	case c.Workers < 0:
		return fmt.Errorf("invalid workers %d", c.Workers)
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
	qualities    func(r *rand.Rand, n int) []byte
	applyErrors  func(r *rand.Rand, seq, qual []byte) ([]byte, []byte)
	header       func(r *rand.Rand, rec *seqread.Record) (string, string)
	observer     ReadObserver
	maxPhred     byte

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
}

// newPipeline instantiates the hooks for one output stream.
func (e *Engine) newPipeline(opts GenerateOpts, obs ReadObserver) (*readPipeline, error) {
	p := &readPipeline{e: e, observer: obs, sampleID: opts.SampleID, region: opts.RefRegion, maxPhred: e.maxPhred(), fixedLength: e.fixedLengthReads()}
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

// writeBatch writes one batch of reads to bw. mates lists the mates written
// per read: {1} or {2} for one file of a batch, {1, 2} for interleaved pairs.
// file is the output file name recorded on each read ("" for streams).
func (p *readPipeline) writeBatch(bw *bufio.Writer, r *rand.Rand, batch int, mates []int, file string, serial *int64) error {
	e := p.e
	rec := seqread.Record{Batch: batch, ReadGroup: p.sampleID, File: file}
	if p.region != nil {
		rec.Positioned = true
		rec.Chromosome = p.region.Chromosome
	}
	for readNum := 1; readNum <= e.config.ReadsPerBatch; readNum++ {
		length := p.sampleLength(r)
		var readStart int64
		if p.region != nil {
			span := p.region.End - p.region.Start - int64(length)
			if span <= 0 {
				return fmt.Errorf("region [%s:%d-%d) too small for read length %d",
					p.region.Chromosome, p.region.Start, p.region.End, length)
			}
			readStart = p.region.Start + r.Int64N(span)
		}
		for _, mate := range mates {
			seq, qual, err := p.read(r, length, readStart)
			if err != nil {
				return err
			}
			rec.Index, rec.Mate, rec.Serial, rec.RefStart = readNum, mate, *serial, readStart
			rec.Seq, rec.Qual = seq, qual
			rec.Name, rec.Comment = p.header(r, &rec)
			if p.sampleID != "" {
				if rec.Comment != "" {
					rec.Comment += "\t"
				}
				rec.Comment += "RG:Z:" + p.sampleID
			}
			if err := p.writeRecord(bw, &rec); err != nil {
				return err
			}
			if p.observer != nil {
				if err := p.observer.ObserveRead(&rec); err != nil {
					return fmt.Errorf("observing read: %w", err)
				}
			}
			*serial++
		}
	}
	return nil
}

// read produces the final bases and Phred qualities of one read.
func (p *readPipeline) read(r *rand.Rand, length int, readStart int64) ([]byte, []byte, error) {
	e := p.e
	n := length
	var local []VariantRecord // positioned variants that can touch this read
	if len(p.positioned) > 0 {
		hi := readStart + int64(length)
		if p.fixedLength {
			hi += int64(p.maxSlack)
		}
		var spans []seqread.Span
		local, spans = p.window(readStart, hi)
		if p.fixedLength {
			n += seqread.WindowSlack(readStart, length, spans)
		}
	} else if p.fixedLength && p.region == nil {
		n += p.blindSlack
	}
	seq, qual, err := e.generateWith(r, n)
	if err != nil {
		return nil, nil, fmt.Errorf("failed to generate sequence: %w", err)
	}
	if len(qual) != len(seq) {
		qual = p.qualities(r, len(seq))
	}
	switch {
	case len(local) > 0:
		var skipped int
		seq, qual, skipped = injectAtPosition(r, seq, qual, readStart, local)
		if skipped > 0 {
			e.overlapSkips.Add(int64(skipped))
		}
	case len(p.blind) > 0:
		seq, qual = injectVariants(r, seq, qual, p.blind)
	}
	if p.fixedLength {
		if len(seq) < length { // cannot happen with correct slack; never pad with N
			extra, extraQ, err := e.generateWith(r, length-len(seq))
			if err != nil {
				return nil, nil, fmt.Errorf("failed to generate sequence: %w", err)
			}
			if len(extraQ) != len(extra) {
				extraQ = p.qualities(r, len(extra))
			}
			seq, qual = append(seq, extra...), append(qual, extraQ...)
		}
		seq, qual = seq[:length], qual[:length]
	}
	if p.applyErrors != nil {
		seq, qual = p.applyErrors(r, seq, qual)
	}
	return seq, qual, nil
}

// writeRecord writes rec as a FASTQ record, applying read-ID fingerprinting.
func (p *readPipeline) writeRecord(bw *bufio.Writer, rec *seqread.Record) error {
	e := p.e
	header := rec.Header()
	if e.config.Fingerprinting.Enable && e.fingerprintingEngine != nil {
		if e.config.Fingerprinting.CustomerFingerprint {
			header = e.fingerprintingEngine.EmbedCustomerID(header)
		}
		if e.config.Fingerprinting.DatasetFingerprint && e.config.DatasetID != "" {
			header = e.fingerprintingEngine.EmbedDatasetFingerprint(header, e.config.DatasetID)
		}
	}
	bw.WriteString(header) // #nosec G104 -- bufio errors are sticky and surface below
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
