package transcriptome

import (
	"fmt"
	"sort"
	"strings"
)

// Exon is one exon in genomic coordinates: 1-based, inclusive, as in GTF.
type Exon struct {
	Start, End int64
}

// Len returns the exon length in bases.
func (e Exon) Len() int { return int(e.End - e.Start + 1) }

// Transcript is one annotated transcript.
type Transcript struct {
	ID     string
	GeneID string
	Seq    []byte // mRNA sense, 5′→3′, upper case

	// Exon structure, when known (from a GTF). Exons are in transcript
	// order, 5′→3′: ascending genomic position on '+', descending on '−'.
	// Their lengths sum to len(Seq). Chrom and Strand are empty / 0 when the
	// structure is unknown.
	Chrom  string
	Strand byte // '+' or '-'
	Exons  []Exon
}

// Len returns the transcript length.
func (t *Transcript) Len() int { return len(t.Seq) }

// HasStructure reports whether the exon structure is known.
func (t *Transcript) HasStructure() bool { return len(t.Exons) > 0 }

// exonOffsets returns the transcript coordinate (0-based) at which each
// exon starts.
func (t *Transcript) exonOffsets() []int {
	off := make([]int, len(t.Exons))
	n := 0
	for i, e := range t.Exons {
		off[i] = n
		n += e.Len()
	}
	return off
}

// Table is the in-memory transcript table, in input order.
type Table struct {
	Transcripts []*Transcript
	byID        map[string]int
}

// NewTable builds a table from transcripts, rejecting empty, duplicate or
// inconsistent entries.
func NewTable(ts []*Transcript) (*Table, error) {
	t := &Table{byID: make(map[string]int, len(ts))}
	for _, tx := range ts {
		if err := t.add(tx); err != nil {
			return nil, err
		}
	}
	if len(t.Transcripts) == 0 {
		return nil, fmt.Errorf("no transcripts")
	}
	return t, nil
}

func (t *Table) add(tx *Transcript) error {
	if tx.ID == "" {
		return fmt.Errorf("transcript with an empty ID")
	}
	if len(tx.Seq) == 0 {
		return fmt.Errorf("transcript %s has an empty sequence", tx.ID)
	}
	if _, dup := t.byID[tx.ID]; dup {
		return fmt.Errorf("duplicate transcript ID %s", tx.ID)
	}
	if tx.GeneID == "" {
		tx.GeneID = tx.ID
	}
	if tx.HasStructure() {
		if n := exonTotal(tx.Exons); n != len(tx.Seq) {
			return fmt.Errorf("transcript %s: exon lengths sum to %d but the sequence has %d bases", tx.ID, n, len(tx.Seq))
		}
	}
	t.byID[tx.ID] = len(t.Transcripts)
	t.Transcripts = append(t.Transcripts, tx)
	return nil
}

// Lookup returns the transcript with id, or nil.
func (t *Table) Lookup(id string) *Transcript {
	if i, ok := t.byID[id]; ok {
		return t.Transcripts[i]
	}
	return nil
}

// Index returns the position of id in the table, or -1.
func (t *Table) Index(id string) int {
	if i, ok := t.byID[id]; ok {
		return i
	}
	return -1
}

func exonTotal(exons []Exon) int {
	n := 0
	for _, e := range exons {
		n += e.Len()
	}
	return n
}

// ParseFASTAHeader extracts the transcript and gene IDs from a transcript
// FASTA header line (without '>'):
//
//   - GENCODE transcripts.fa: "ENST…|ENSG…|OTTHUMG…|OTTHUMT…|name-201|name|len|biotype|"
//     gives the first two '|' fields;
//   - Ensembl cDNA: "ENST… cdna chromosome:… gene:ENSG… …" gives the first
//     token and the gene: tag;
//   - anything else: the first whitespace-delimited token is the transcript
//     ID and the gene ID equals it.
func ParseFASTAHeader(h string) (id, gene string) {
	tok := h
	rest := ""
	if i := strings.IndexAny(h, " \t"); i >= 0 {
		tok, rest = h[:i], h[i+1:]
	}
	if strings.Contains(tok, "|") {
		f := strings.Split(tok, "|")
		id = f[0]
		if len(f) > 1 && f[1] != "" && f[1] != "-" {
			gene = f[1]
		}
	} else {
		id = tok
	}
	if gene == "" {
		for _, w := range strings.Fields(rest) {
			if g, ok := strings.CutPrefix(w, "gene:"); ok && g != "" {
				gene = g
				break
			}
		}
	}
	if gene == "" {
		gene = id
	}
	return id, gene
}

// LoadTranscriptFASTA reads a transcript FASTA (plain or gzip) into a
// table. Sequences are upper-cased; U is read as T.
func LoadTranscriptFASTA(path string) (*Table, error) {
	rc, err := openInput(path)
	if err != nil {
		return nil, err
	}
	defer rc.Close() // #nosec G307 -- read-only
	var ts []*Transcript
	err = readFASTA(rc, func(h string, seq []byte) error {
		id, gene := ParseFASTAHeader(h)
		for i, c := range seq {
			if c == 'U' {
				seq[i] = 'T'
			}
		}
		ts = append(ts, &Transcript{ID: id, GeneID: gene, Seq: seq})
		return nil
	})
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	t, err := NewTable(ts)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return t, nil
}

// AttachAnnotation copies exon structures (and gene IDs) from ann onto the
// table's transcripts with the same ID whose sequence length equals the
// annotated exon total. It returns how many transcripts were annotated and
// how many matched by ID but had a different length (left without
// structure).
func (t *Table) AttachAnnotation(ann *Annotation) (attached, mismatched int) {
	for _, tx := range t.Transcripts {
		a := ann.Lookup(tx.ID)
		if a == nil {
			continue
		}
		if exonTotal(a.Exons) != len(tx.Seq) {
			mismatched++
			continue
		}
		tx.Chrom, tx.Strand = a.Chrom, a.Strand
		tx.Exons = append([]Exon(nil), a.Exons...)
		if a.GeneID != "" {
			tx.GeneID = a.GeneID
		}
		attached++
	}
	return attached, mismatched
}

// sortExonsTranscriptOrder sorts exons ascending by start, then reverses
// them for minus-strand transcripts, giving 5′→3′ transcript order.
func sortExonsTranscriptOrder(exons []Exon, strand byte) {
	sort.Slice(exons, func(i, j int) bool { return exons[i].Start < exons[j].Start })
	if strand == '-' {
		for i, j := 0, len(exons)-1; i < j; i, j = i+1, j-1 {
			exons[i], exons[j] = exons[j], exons[i]
		}
	}
}
