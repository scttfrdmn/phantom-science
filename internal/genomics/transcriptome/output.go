package transcriptome

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"

	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore"
	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore/errmodel"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
)

// Output file names inside the run directory.
const (
	FASTQName      = "reads.fastq"
	TruthName      = "truth.tsv"
	CountsName     = "counts.tsv"
	NovelTSVName   = "novel_isoforms.tsv"
	NovelGTFName   = "novel_isoforms.gtf"
	SummaryName    = nanopore.SummaryFileName
	gtfSourceField = "phantom_science"
)

// TruthColumns are the truth.tsv columns. Coordinates are 0-based,
// half-open, on the isoform's sequence (the annotated transcript, or the
// novel isoform for novel reads): the read covers [start, end) and end is
// always isoform_length (3′ capture).
var TruthColumns = []string{
	"read_id", "transcript_id", "gene_id", "isoform_id", "novel_isoform",
	"strand", "full_length", "trunc5p", "polya_len", "start", "end",
	"isoform_length", "read_length",
}

// CountsColumns are the counts.tsv columns: one row per molecule species
// (every annotated transcript, then novel isoforms after their parent).
// tpm is the true molecule share × 10^6 (sums to 10^6); expected_reads is
// reads × tpm / 10^6; reads is the realised count (sums to the run's reads).
var CountsColumns = []string{
	"transcript_id", "gene_id", "parent_transcript_id", "novel_isoform",
	"length", "tpm", "expected_reads", "reads",
}

// NovelColumns are the novel_isoforms.tsv columns.
var NovelColumns = []string{
	"isoform_id", "parent_transcript_id", "gene_id", "novel_isoform", "detail",
	"length", "chrom", "strand", "intron_chain", "tpm", "reads",
}

// OutputOptions configures WriteRun.
type OutputOptions struct {
	Gzip bool // gzip the FASTQ (reads.fastq.gz)
}

// RunResult summarises a written run.
type RunResult struct {
	Dir        string
	FASTQ      string
	Reads      int
	Bases      int64
	NovelReads int
	Counts     map[*Species]int
}

// WriteRun generates the run into dir (created if needed): the FASTQ,
// sequencing_summary.txt, truth.tsv, counts.tsv and novel_isoforms.tsv
// (plus novel_isoforms.gtf when isoforms have genomic coordinates).
func WriteRun(dir string, s *Simulator, opts OutputOptions) (res *RunResult, err error) {
	if err := os.MkdirAll(dir, 0o750); err != nil {
		return nil, err
	}
	res = &RunResult{Dir: dir, Counts: map[*Species]int{}}
	fqPath := filepath.Join(dir, FASTQName)
	if opts.Gzip {
		fqPath += ".gz"
	}
	res.FASTQ = fqPath

	var closers []io.Closer
	defer func() {
		for i := len(closers) - 1; i >= 0; i-- {
			err = errors.Join(err, closers[i].Close())
		}
	}()
	create := func(path string) (*bufio.Writer, error) {
		f, err := os.Create(path) // #nosec G304 -- output path chosen by the user
		if err != nil {
			return nil, err
		}
		bw := bufio.NewWriterSize(f, 1<<20)
		closers = append(closers, flushCloser{bw, f})
		return bw, nil
	}

	fqFile, err := os.Create(fqPath) // #nosec G304 -- output path chosen by the user
	if err != nil {
		return nil, err
	}
	var fq *bufio.Writer
	if opts.Gzip {
		bw := bufio.NewWriterSize(fqFile, 1<<20)
		gz := gzip.NewWriter(bw)
		closers = append(closers, flushCloser{bw, fqFile})
		fq = bufio.NewWriterSize(gz, 1<<20)
		closers = append(closers, flushCloser{fq, gz})
	} else {
		fq = bufio.NewWriterSize(fqFile, 1<<20)
		closers = append(closers, flushCloser{fq, fqFile})
	}

	summary := nanopore.NewSummaryWriter(s.RunInfo(), errmodel.MaxPhred, func(*seqread.Record) (io.WriteCloser, error) {
		return os.Create(filepath.Join(dir, SummaryName)) // #nosec G304 -- output path chosen by the user
	}, nil)
	closers = append(closers, summary)

	truth, err := create(filepath.Join(dir, TruthName))
	if err != nil {
		return nil, err
	}
	writeRow(truth, TruthColumns...)

	err = s.Generate(func(rd *Read) error {
		rd.Record.File = fqPath
		writeFASTQ(fq, &rd.Record)
		if err := summary.ObserveRead(&rd.Record); err != nil {
			return err
		}
		t := rd.Truth
		sp := t.Species
		strand := "+"
		if t.Reverse {
			strand = "-"
		}
		writeRow(truth, rd.Record.Name, sp.Parent.ID, sp.Parent.GeneID, sp.ID, sp.Novel, strand,
			strconv.FormatBool(t.FullLength), strconv.Itoa(t.Trunc5p), strconv.Itoa(t.PolyA),
			strconv.Itoa(t.Start), strconv.Itoa(t.End), strconv.Itoa(len(sp.Seq)), strconv.Itoa(len(rd.Record.Seq)))
		res.Counts[sp]++
		res.Reads++
		res.Bases += int64(len(rd.Record.Seq))
		if sp.Novel != NovelNone {
			res.NovelReads++
		}
		return nil
	})
	if err != nil {
		return nil, err
	}
	if err := writeCounts(dir, s, res, create); err != nil {
		return nil, err
	}
	return res, nil
}

func writeCounts(dir string, s *Simulator, res *RunResult, create func(string) (*bufio.Writer, error)) error {
	c, err := create(filepath.Join(dir, CountsName))
	if err != nil {
		return err
	}
	n, err := create(filepath.Join(dir, NovelTSVName))
	if err != nil {
		return err
	}
	writeRow(c, CountsColumns...)
	writeRow(n, NovelColumns...)
	var gtf *bufio.Writer
	N := float64(s.cfg.Reads)
	for _, sp := range s.species {
		parent := ""
		if sp.Novel != NovelNone {
			parent = sp.Parent.ID
		}
		writeRow(c, sp.ID, sp.Parent.GeneID, parent, sp.Novel, strconv.Itoa(len(sp.Seq)),
			fmtFloat(sp.TPM), fmtFloat(N*sp.TPM/1e6), strconv.Itoa(res.Counts[sp]))
		if sp.Novel == NovelNone {
			continue
		}
		strand := ""
		if sp.Parent.Strand != 0 {
			strand = string(sp.Parent.Strand)
		}
		writeRow(n, sp.ID, sp.Parent.ID, sp.Parent.GeneID, sp.Novel, sp.Detail, strconv.Itoa(len(sp.Seq)),
			sp.Parent.Chrom, strand, sp.IntronChain(), fmtFloat(sp.TPM), strconv.Itoa(res.Counts[sp]))
		if len(sp.Exons) == 0 {
			continue
		}
		if gtf == nil {
			if gtf, err = create(filepath.Join(dir, NovelGTFName)); err != nil {
				return err
			}
			gtf.WriteString("##description: novel isoforms injected by phantom-science transcriptome\n") // #nosec G104 -- sticky
		}
		writeGTF(gtf, sp)
	}
	return nil
}

// writeGTF writes a transcript line and its exon lines.
func writeGTF(w *bufio.Writer, sp *Species) {
	lo, hi := sp.Exons[0].Start, sp.Exons[0].End
	for _, e := range sp.Exons {
		lo, hi = min(lo, e.Start), max(hi, e.End)
	}
	attrs := fmt.Sprintf(`gene_id "%s"; transcript_id "%s"; parent_transcript_id "%s"; novel_isoform "%s";`,
		sp.Parent.GeneID, sp.ID, sp.Parent.ID, sp.Novel)
	strand := string(sp.Parent.Strand)
	writeRow(w, sp.Parent.Chrom, gtfSourceField, "transcript", strconv.FormatInt(lo, 10), strconv.FormatInt(hi, 10), ".", strand, ".", attrs)
	for i, e := range sp.Exons {
		writeRow(w, sp.Parent.Chrom, gtfSourceField, "exon", strconv.FormatInt(e.Start, 10), strconv.FormatInt(e.End, 10), ".", strand, ".",
			attrs+fmt.Sprintf(` exon_number "%d";`, i+1))
	}
}

func fmtFloat(v float64) string { return strconv.FormatFloat(v, 'g', 10, 64) }

// writeRow writes tab-separated fields and a newline. bufio errors are
// sticky and reported when the writer is flushed.
func writeRow(w *bufio.Writer, fields ...string) {
	for i, f := range fields {
		if i > 0 {
			w.WriteByte('\t') // #nosec G104 -- sticky
		}
		w.WriteString(f) // #nosec G104 -- sticky
	}
	w.WriteByte('\n') // #nosec G104 -- sticky
}

// writeFASTQ writes rec with Phred+33 qualities.
func writeFASTQ(w *bufio.Writer, rec *seqread.Record) {
	w.WriteString(rec.Header()) // #nosec G104 -- sticky
	w.WriteByte('\n')           // #nosec G104
	w.Write(rec.Seq)            // #nosec G104
	w.WriteString("\n+\n")      // #nosec G104
	for _, q := range rec.Qual {
		w.WriteByte(min(q, errmodel.MaxPhred) + 33) // #nosec G104
	}
	w.WriteByte('\n') // #nosec G104
}

// flushCloser flushes a buffer and closes what it wraps.
type flushCloser struct {
	bw *bufio.Writer
	c  io.Closer
}

func (f flushCloser) Close() error { return errors.Join(f.bw.Flush(), f.c.Close()) }
