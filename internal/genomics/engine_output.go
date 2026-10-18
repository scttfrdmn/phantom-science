package genomics

// Output extension points: read headers, file names, stream wrapping, file
// finalisation and read observers (#12, #16).
//
// Contract:
//
//   - checkOutput validates the platform and output settings before any
//     output is written. Planned-but-missing features wrap ErrNotImplemented.
//   - headerFormatter is called once per output stream and returns a function
//     called once per read (each mate separately) with the stream's generator
//     and the read's record (provenance, Seq and Qual already final). It
//     returns the read name and header comment; the header line is
//     "@name comment". It may also set platform fields on rec (Channel,
//     StartTime, Duration) for observers. It must draw only from r. The
//     pipeline appends "\tRG:Z:<id>" to the comment when a read group is set,
//     then applies read-ID fingerprinting to the whole line. The default
//     writes Illumina headers "<sample>:<lane>:HXXXXXX:<batch>:1101:<x>:<y>"
//     "<mate>:N:0:1"; x/y are read*1000/read*2000, or, for reads placed on a
//     region, start%1e5 and (start/1e5)%1e5.
//   - maxPhred caps Phred scores when they are encoded as ASCII (Phred+33).
//   - outputPath returns the file written for one batch and mate (mate is 1,
//     or 2 for the second file of a paired-end batch).
//   - wrapOutput wraps every output stream (each file, or the GenerateToWriter
//     writer) before records are written, e.g. for streaming compression. The
//     pipeline closes the returned writer (not the underlying one) when the
//     stream is complete. The default passes writes through.
//   - finalize runs on each closed output file and returns the path of the
//     finished file (which may differ, e.g. after compressing). The default
//     embeds the quality-score fingerprint when that layer is enabled; any
//     compression of fingerprinted files must happen after embedding, here.
//   - observer returns the run's ReadObserver (nil for none, the default). It
//     is opened once per Generate / GenerateToWriter / GenerateFASTQFile call,
//     sees every written read, and is closed at the end of that call (e.g. to
//     write an ONT sequencing_summary.txt).

import (
	"fmt"
	"io"
	"path/filepath"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
	"math/rand/v2"
)

// Platform names for SequencerConfig.Platform.
const (
	platformIllumina = "illumina"
	platformONT      = "ont"
)

// maxIlluminaPhred is the highest Phred score written for Illumina reads ('J').
const maxIlluminaPhred = 41

// ReadObserver sees every read the pipeline writes. Implementations must be
// safe for concurrent use: output streams may be generated in parallel.
type ReadObserver interface {
	// ObserveRead is called after rec has been written. rec and its slices
	// are only valid during the call.
	ObserveRead(rec *seqread.Record) error
	// Close is called once when the run's output is complete.
	Close() error
}

// checkOutput validates the platform and output configuration.
func (e *Engine) checkOutput() error {
	switch p := strings.ToLower(e.config.Platform); p {
	case "", platformIllumina:
	case platformONT:
		return fmt.Errorf("platform %q output is %w (#12)", e.config.Platform, ErrNotImplemented)
	default:
		return fmt.Errorf("unknown platform %q: want illumina or ont", e.config.Platform)
	}
	out := e.config.Output
	switch c := strings.ToLower(out.Compression); c {
	case "", "none":
	case "gzip":
		return fmt.Errorf("output compression %q is %w (#12)", out.Compression, ErrNotImplemented)
	default:
		return fmt.Errorf("unknown output compression %q: want none or gzip", out.Compression)
	}
	if out.FilePattern != "" {
		return fmt.Errorf("output file_pattern is %w (#16)", ErrNotImplemented)
	}
	if out.SummaryFile != "" {
		return fmt.Errorf("output summary_file (sequencing_summary.txt) is %w (#12)", ErrNotImplemented)
	}
	return nil
}

// headerFormatter returns the per-read header function for the platform.
func (e *Engine) headerFormatter() (func(r *rand.Rand, rec *seqread.Record) (name, comment string), error) {
	sample, lane := e.config.SampleName, e.config.Lane
	return func(_ *rand.Rand, rec *seqread.Record) (string, string) {
		x, y := rec.Index*1000, rec.Index*2000
		if rec.Positioned {
			x, y = int(rec.RefStart%100000), int((rec.RefStart/100000)%100000)
		}
		return fmt.Sprintf("%s:%d:HXXXXXX:%d:1101:%d:%d", sample, lane, rec.Batch, x, y),
			fmt.Sprintf("%d:N:0:1", rec.Mate)
	}, nil
}

// maxPhred returns the highest Phred score encoded in FASTQ output.
func (e *Engine) maxPhred() byte {
	return maxIlluminaPhred
}

// outputPath returns the FASTQ file for batch and mate.
func (e *Engine) outputPath(batch, mate int) string {
	c := e.config
	if c.PairedEnd {
		return filepath.Join(c.OutputDir, fmt.Sprintf("%s_S1_L%03d_R%d_%03d.fastq", c.SampleName, c.Lane, mate, batch))
	}
	return filepath.Join(c.OutputDir, fmt.Sprintf("%s_S1_L%03d_%03d.fastq", c.SampleName, c.Lane, batch))
}

// wrapOutput wraps an output stream before records are written to it.
func (e *Engine) wrapOutput(w io.Writer) (io.WriteCloser, error) {
	return nopWriteCloser{w}, nil
}

// finalize finishes a closed output file and returns its final path.
func (e *Engine) finalize(path string) (string, error) {
	if e.qualityFingerprintingEnabled() {
		if err := e.embedFileFingerprint(path); err != nil {
			return "", err
		}
	}
	return path, nil
}

// observer returns the run's read observer, or nil for none.
func (e *Engine) observer() (ReadObserver, error) {
	return nil, nil
}

type nopWriteCloser struct{ io.Writer }

func (nopWriteCloser) Close() error { return nil }
