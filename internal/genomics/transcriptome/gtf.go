package transcriptome

import (
	"bytes"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// AnnotatedTranscript is one transcript's exon structure from a GTF.
type AnnotatedTranscript struct {
	ID, GeneID string
	Chrom      string
	Strand     byte
	Exons      []Exon // transcript order (5′→3′)
}

// Annotation is the exon structure of every transcript in a GTF, in order of
// first appearance.
type Annotation struct {
	Transcripts []*AnnotatedTranscript
	byID        map[string]int
	// Skipped counts transcripts dropped for an unknown strand ('.').
	Skipped int
}

// Lookup returns the annotated transcript with id, or nil.
func (a *Annotation) Lookup(id string) *AnnotatedTranscript {
	if i, ok := a.byID[id]; ok {
		return a.Transcripts[i]
	}
	return nil
}

// gtfAttr returns the value of key in a GTF attribute column
// (`key "value"; key2 "value2";`), or "".
func gtfAttr(attrs, key string) string {
	for _, f := range strings.Split(attrs, ";") {
		f = strings.TrimSpace(f)
		k, v, ok := strings.Cut(f, " ")
		if !ok || k != key {
			continue
		}
		return strings.Trim(strings.TrimSpace(v), `"`)
	}
	return ""
}

// ParseGTF reads the exon records of a GTF (GENCODE/Ensembl GTF 2.2; plain
// or gzip). Only feature "exon" is used; each exon must carry transcript_id
// (gene_id defaults to the transcript ID). All exons of a transcript must
// share chromosome and strand and must not overlap. Transcripts with strand
// '.' are skipped and counted in Skipped.
func ParseGTF(path string) (*Annotation, error) {
	rc, err := openInput(path)
	if err != nil {
		return nil, err
	}
	defer rc.Close() // #nosec G307 -- read-only
	a, err := parseGTF(rc)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return a, nil
}

func parseGTF(r io.Reader) (*Annotation, error) {
	a := &Annotation{byID: make(map[string]int)}
	unstranded := map[string]bool{}
	lr := newLineReader(r)
	for {
		line, err := lr.next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		f := bytes.Split(line, []byte{'\t'})
		if len(f) < 9 {
			return nil, fmt.Errorf("line %d: GTF needs 9 tab-separated columns, got %d", lr.line, len(f))
		}
		if string(f[2]) != "exon" {
			continue
		}
		start, err1 := strconv.ParseInt(string(f[3]), 10, 64)
		end, err2 := strconv.ParseInt(string(f[4]), 10, 64)
		if err1 != nil || err2 != nil || start < 1 || end < start {
			return nil, fmt.Errorf("line %d: bad exon coordinates %q-%q", lr.line, f[3], f[4])
		}
		attrs := string(f[8])
		id := gtfAttr(attrs, "transcript_id")
		if id == "" {
			return nil, fmt.Errorf("line %d: exon without transcript_id", lr.line)
		}
		gene := gtfAttr(attrs, "gene_id")
		chrom := string(f[0])
		var strand byte
		switch string(f[6]) {
		case "+":
			strand = '+'
		case "-":
			strand = '-'
		default:
			unstranded[id] = true
			continue
		}
		i, ok := a.byID[id]
		if !ok {
			i = len(a.Transcripts)
			a.byID[id] = i
			a.Transcripts = append(a.Transcripts, &AnnotatedTranscript{ID: id, GeneID: gene, Chrom: chrom, Strand: strand})
		}
		t := a.Transcripts[i]
		if t.Chrom != chrom || t.Strand != strand {
			return nil, fmt.Errorf("line %d: transcript %s has exons on %s%c and %s%c", lr.line, id, t.Chrom, t.Strand, chrom, strand)
		}
		if t.GeneID == "" {
			t.GeneID = gene
		}
		t.Exons = append(t.Exons, Exon{Start: start, End: end})
	}
	for _, t := range a.Transcripts {
		if t.GeneID == "" {
			t.GeneID = t.ID
		}
		sortExonsTranscriptOrder(t.Exons, t.Strand)
		for j := 1; j < len(t.Exons); j++ {
			lo, hi := t.Exons[j-1], t.Exons[j]
			if t.Strand == '-' {
				lo, hi = hi, lo
			}
			if hi.Start <= lo.End {
				return nil, fmt.Errorf("transcript %s has overlapping exons %d-%d and %d-%d", t.ID, lo.Start, lo.End, hi.Start, hi.End)
			}
		}
	}
	for id := range unstranded {
		if _, ok := a.byID[id]; !ok {
			a.Skipped++
		}
	}
	if len(a.Transcripts) == 0 {
		return nil, fmt.Errorf("no stranded exon records")
	}
	return a, nil
}

// BuildStats reports what BuildFromGenome could not build.
type BuildStats struct {
	MissingChrom int // transcripts on chromosomes absent from the genome
}

// BuildFromGenome extracts every annotated transcript from a genome FASTA
// (plain or gzip): exons are concatenated in transcript order and minus-strand
// transcripts are reverse-complemented. The genome is streamed one
// chromosome at a time. Transcripts on chromosomes the genome lacks are
// skipped and counted; an exon past the end of its chromosome is an error.
func BuildFromGenome(genomePath string, ann *Annotation) (*Table, BuildStats, error) {
	var st BuildStats
	byChrom := map[string][]int{}
	for i, t := range ann.Transcripts {
		byChrom[t.Chrom] = append(byChrom[t.Chrom], i)
	}
	seqs := make([][]byte, len(ann.Transcripts))
	rc, err := openInput(genomePath)
	if err != nil {
		return nil, st, err
	}
	defer rc.Close() // #nosec G307 -- read-only
	err = readFASTA(rc, func(h string, chrom []byte) error {
		name := h
		if i := strings.IndexAny(h, " \t"); i >= 0 {
			name = h[:i]
		}
		for _, i := range byChrom[name] {
			t := ann.Transcripts[i]
			s, err := spliceTranscript(chrom, t.Exons, t.Strand)
			if err != nil {
				return fmt.Errorf("transcript %s: %w", t.ID, err)
			}
			seqs[i] = s
		}
		delete(byChrom, name)
		return nil
	})
	if err != nil {
		return nil, st, fmt.Errorf("%s: %w", genomePath, err)
	}
	var ts []*Transcript
	for i, a := range ann.Transcripts {
		if seqs[i] == nil {
			st.MissingChrom++
			continue
		}
		ts = append(ts, &Transcript{
			ID: a.ID, GeneID: a.GeneID, Seq: seqs[i],
			Chrom: a.Chrom, Strand: a.Strand, Exons: append([]Exon(nil), a.Exons...),
		})
	}
	if len(ts) == 0 {
		return nil, st, fmt.Errorf("no annotated transcript lies on a chromosome of %s", genomePath)
	}
	t, err := NewTable(ts)
	return t, st, err
}

// spliceTranscript returns the mRNA-sense sequence of exons (transcript
// order) on chrom.
func spliceTranscript(chrom []byte, exons []Exon, strand byte) ([]byte, error) {
	out := make([]byte, 0, exonTotal(exons))
	if strand == '-' {
		// Transcript order is descending; build the genomic (ascending)
		// concatenation, then reverse-complement it.
		for i := len(exons) - 1; i >= 0; i-- {
			e := exons[i]
			if e.End > int64(len(chrom)) {
				return nil, fmt.Errorf("exon %d-%d past chromosome end %d", e.Start, e.End, len(chrom))
			}
			out = append(out, chrom[e.Start-1:e.End]...)
		}
		return ReverseComplement(out), nil
	}
	for _, e := range exons {
		if e.End > int64(len(chrom)) {
			return nil, fmt.Errorf("exon %d-%d past chromosome end %d", e.Start, e.End, len(chrom))
		}
		out = append(out, chrom[e.Start-1:e.End]...)
	}
	return out, nil
}
