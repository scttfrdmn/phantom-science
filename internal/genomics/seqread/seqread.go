// Package seqread defines the read record that flows through the FASTQ
// generation pipeline, plus allele-splicing helpers that keep base qualities
// aligned with bases. It is a leaf package: it imports only the standard
// library, so the engine, clinical scenarios, platform models and output
// writers can all share it.
package seqread

import "time"

// Record is one synthetic read as it moves through the generation pipeline.
//
// The pipeline fills the provenance fields (Batch .. ReadGroup) before the
// header hook runs, sets Seq and Qual after generation, variant injection,
// trimming and the error model, sets Name and Comment from the header hook,
// and finally hands the record to the read observer. The FASTQ header line is
// "@" + Name, followed by " " + Comment when Comment is not empty.
type Record struct {
	// Name is the read ID (without '@'); Comment is the rest of the header
	// line after the first space (e.g. "1:N:0:1" or ONT key=value fields).
	Name    string
	Comment string

	// Seq holds the bases; Qual holds Phred scores (0..93), NOT ASCII. They
	// always have the same length.
	Seq  []byte
	Qual []byte

	// Batch is the 1-based batch (output file) number, Mate is 1 or 2, and
	// Index is the 1-based read number within the batch (shared by both
	// mates of a pair).
	Batch int
	Mate  int
	Index int
	// Serial is the 0-based number of this read within its output stream
	// (counting each mate separately).
	Serial int64

	// Positioned reports whether the read was placed on a reference region;
	// RefStart is then its 0-based start and Chromosome the region's.
	Positioned bool
	RefStart   int64
	Chromosome string

	// ReadGroup is the sample ID written as RG:Z:<id>; empty for none.
	ReadGroup string

	// File is the output file the read is written to; empty for streamed
	// output.
	File string

	// Platform-specific fields set by header hooks (ONT); zero for Illumina.
	Channel   int
	StartTime time.Time
	Duration  time.Duration
}

// Header returns the FASTQ header line for r, without the trailing newline.
func (r *Record) Header() string {
	if r.Comment == "" {
		return "@" + r.Name
	}
	return "@" + r.Name + " " + r.Comment
}

// MeanQuality returns the mean Phred score of r, or 0 for an empty read.
// (This is the arithmetic mean of Phred values, not the error-probability
// mean that some basecallers report.)
func (r *Record) MeanQuality() float64 {
	if len(r.Qual) == 0 {
		return 0
	}
	sum := 0
	for _, q := range r.Qual {
		sum += int(q)
	}
	return float64(sum) / float64(len(r.Qual))
}

// Splice replaces the refLen bases of seq starting at off with allele and
// returns the new sequence and qualities. qual may be nil; when it is not, it
// must have the same length as seq and stays aligned: bases that replace
// reference bases keep those positions' qualities, inserted bases (allele
// longer than refLen) take the quality of the last replaced base, and the
// qualities of deleted bases are dropped. The inputs are not modified.
//
// The caller must ensure 0 <= off and off+refLen <= len(seq).
func Splice(seq, qual []byte, off, refLen int, allele []byte) ([]byte, []byte) {
	n := len(seq) - refLen + len(allele)
	out := make([]byte, 0, n)
	out = append(out, seq[:off]...)
	out = append(out, allele...)
	out = append(out, seq[off+refLen:]...)
	if qual == nil {
		return out, nil
	}
	q := make([]byte, 0, n)
	q = append(q, qual[:off]...)
	common := min(refLen, len(allele))
	q = append(q, qual[off:off+common]...)
	if len(allele) > refLen {
		fill := byte(0)
		switch {
		case refLen > 0:
			fill = qual[off+refLen-1]
		case off > 0:
			fill = qual[off-1]
		case off < len(qual):
			fill = qual[off]
		}
		for i := refLen; i < len(allele); i++ {
			q = append(q, fill)
		}
	}
	q = append(q, qual[off+refLen:]...)
	return out, q
}

// complement maps a base to its complement; bytes other than ACGT (either
// case) map to N.
var complement = func() (t [256]byte) {
	for i := range t {
		t[i] = 'N'
	}
	for _, p := range []string{"AT", "CG", "GC", "TA", "at", "cg", "gc", "ta"} {
		t[p[0]] = p[1]
	}
	return t
}()

// ReverseComplement returns the reverse complement of seq in a new slice.
// N (and any byte that is not a base) becomes N.
func ReverseComplement(seq []byte) []byte {
	out := make([]byte, len(seq))
	for i, b := range seq {
		out[len(seq)-1-i] = complement[b]
	}
	return out
}

// Span is a reference interval [Pos, Pos+RefLen) (0-based) that a variant
// replaces with an allele of AltLen bases.
type Span struct {
	Pos    int64
	RefLen int
	AltLen int
}

// WindowSlack returns how many bases beyond length a read starting at the
// 0-based position start must be generated so that, after every deletion in
// its window is applied, at least length bases remain. Only spans starting at
// or after start can shorten the read; the window grows as deletions pull
// more reference into it, so the computation iterates to a fixed point.
func WindowSlack(start int64, length int, spans []Span) int {
	slack := 0
	for {
		next := 0
		end := start + int64(length+slack)
		for _, s := range spans {
			if s.Pos >= start && s.Pos < end && s.RefLen > s.AltLen {
				next += s.RefLen - s.AltLen
			}
		}
		if next <= slack {
			return slack
		}
		slack = next
	}
}
