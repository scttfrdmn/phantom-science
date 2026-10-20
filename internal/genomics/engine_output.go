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
//     then, with read-ID fingerprinting enabled, the tab-separated SAM-style
//     fields "\tXC:Z:<customer>" and "\tXD:Z:<dataset>" (see
//     fingerprinting.ReadTags); the read name is never changed. The default
//     writes Illumina headers "<sample>:<lane>:HXXXXXX:<batch>:1101:<x>:<y>"
//     "<mate>:N:0:1"; x/y are read*1000/read*2000, or, for reads placed on a
//     region (GenerateOpts.RefRegion), start%1e5 and (start/1e5)%1e5. Reads
//     sampled from a reference keep read*1000/read*2000, so names stay
//     unique.
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
//     is opened once per GenerateToWriter / GenerateFASTQFile call, sees
//     every written read, and is closed at the end of that call (e.g. to
//     write an ONT sequencing_summary.txt). Generate instead opens one
//     summaryObserver per output unit, writing to a part file, and
//     concatenates the parts in unit order (see generateRun), so parallel
//     generation gives the same summary.
//
// Oxford Nanopore output (Platform "ont", #12):
//
//   - Headers are MinKNOW's: "@<uuid4> runid=<hex> read=<n> ch=<channel>
//     start_time=<RFC 3339> flow_cell_id=<id> protocol_group_id=<id>
//     sample_id=<id>" (see package nanopore for the per-channel timing
//     model). Read-ID fingerprint tags, when enabled, follow as tab-separated
//     SAM-style fields (XC:Z:, XD:Z:), as for Illumina.
//   - Reads are single-end; paired_end is rejected. Phred scores are capped
//     at 50 (Dorado's maximum) instead of Illumina's 41.
//   - The observer writes sequencing_summary.txt: for file output to
//     Output.SummaryFile (default "sequencing_summary.txt") in the output
//     directory; for GenerateToWriter only when Output.SummaryFile is set
//     (relative to the output directory when one is configured, else to the
//     working directory). Each call truncates the summary, so callers that
//     write one file per GenerateFASTQFile call get the last call's reads.
//   - Default file names are "<sample>_<flowcell>_<run8>_<n>.fastq", n
//     0-based as in MinKNOW.
//
// Compression "gzip" (any platform): files are named "<name>.gz". Without the
// quality-score fingerprint layer the stream is gzip-compressed while it is
// written (wrapOutput). With it, the plain FASTQ is written, the fingerprint
// is embedded (stegano works on plain FASTQ), then finalize compresses it to
// "<name>.gz" and removes the plain file. GeneratedFiles lists the final
// .gz paths. GenerateToWriter output is gzip-compressed as a stream.
//
// Output.FilePattern (any platform) names output files with placeholders
// {sample_name}, {flowcell_id}, {run_id} (8 characters for ONT, as MinKNOW
// file names), {file_number} (0-based batch), {barcode} ("barcode01"),
// {lane} and {mate}. Multi-batch runs need {file_number}, paired-end runs
// {mate}; ".gz" is appended for gzip output when the pattern lacks it.

import (
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"math/rand/v2"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
	"github.com/scttfrdmn/phantom-science/internal/security"
	"github.com/scttfrdmn/phantom-science/internal/storage/s3sink"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Platform names for SequencerConfig.Platform.
const (
	platformIllumina = "illumina"
	platformONT      = "ont"
)

// maxIlluminaPhred is the highest Phred score written for Illumina reads ('J').
const maxIlluminaPhred = 41

// maxONTPhred is the highest Phred score written for ONT reads (Dorado caps
// its quality scores at 50, 'S').
const maxONTPhred = 50

// Output compression names for OutputSpec.Compression.
const (
	compressionNone = "none"
	compressionGzip = "gzip"
)

// defaultONTFilePattern names ONT output files when no pattern is set.
const defaultONTFilePattern = "{sample_name}_{flowcell_id}_{run_id}_{file_number}.fastq"

// defaultBarcode is the {barcode} file-name placeholder value.
const defaultBarcode = "barcode01"

// ReadObserver sees every read the pipeline writes. Implementations must be
// safe for concurrent use: output streams may be generated in parallel.
type ReadObserver interface {
	// ObserveRead is called after rec has been written. rec and its slices
	// are only valid during the call.
	ObserveRead(rec *seqread.Record) error
	// Close is called once when the run's output is complete.
	Close() error
}

// platform returns the normalised platform name ("" means illumina).
func (e *Engine) platform() string {
	p := strings.ToLower(e.config.Platform)
	if p == "" {
		return platformIllumina
	}
	return p
}

func (e *Engine) isONT() bool { return e.platform() == platformONT }

// gzipOutput reports whether output is gzip-compressed.
func (e *Engine) gzipOutput() bool {
	return strings.EqualFold(e.config.Output.Compression, compressionGzip)
}

// compressAfterEmbed reports whether files are written plain and compressed
// by finalize, after the quality-score fingerprint is embedded.
func (e *Engine) compressAfterEmbed() bool {
	return e.gzipOutput() && e.qualityFingerprintingEnabled()
}

// checkOutput validates the platform and output configuration.
func (e *Engine) checkOutput() error {
	if err := e.checkIllumina(); err != nil {
		return err
	}
	switch e.platform() {
	case platformIllumina:
	case platformONT:
		if e.config.PairedEnd {
			return fmt.Errorf("platform %q reads are single-end: paired_end is not supported", e.config.Platform)
		}
		if err := e.ontRunInfo().Validate(); err != nil {
			return err
		}
	default:
		return fmt.Errorf("unknown platform %q: want illumina or ont", e.config.Platform)
	}
	if s3sink.IsURI(e.config.OutputDir) {
		if _, err := s3sink.ParseURI(e.config.OutputDir); err != nil {
			return fmt.Errorf("output_dir: %w", err)
		}
	}
	out := e.config.Output
	switch c := strings.ToLower(out.Compression); c {
	case "", compressionNone, compressionGzip:
	default:
		return fmt.Errorf("unknown output compression %q: want none or gzip", out.Compression)
	}
	if out.SummaryFile != "" {
		if !e.isONT() {
			return fmt.Errorf("output summary_file is only written for platform ont")
		}
		if filepath.Base(out.SummaryFile) != out.SummaryFile || out.SummaryFile == "." || out.SummaryFile == ".." {
			return fmt.Errorf("output summary_file %q must be a file name (it is written in the output directory)", out.SummaryFile)
		}
	}
	return e.checkFilePattern()
}

// ontRunInfo returns the ONT run metadata, with defaults derived from the
// run seed.
func (e *Engine) ontRunInfo() nanopore.RunInfo {
	o := e.config.ONT
	flowCell := o.FlowCellID
	if flowCell == "" {
		flowCell = e.config.FlowCellID
	}
	return nanopore.RunInfo{
		RunID:           o.RunID,
		FlowCellID:      flowCell,
		ProtocolGroupID: o.ProtocolGroupID,
		SampleID:        e.config.SampleName,
		Channels:        o.Channels,
		SpeedBPS:        o.SpeedBPS,
	}.WithDefaults(e.seed)
}

// slotBases is the read length the ONT channel model sizes each output
// stream's time window for: the fixed read length, or, for distributions,
// twice the largest configured length scale capped at Max. Reads longer
// than this may overlap the next file's reads on the same channel (see
// nanopore.ChannelModel).
func (e *Engine) slotBases() int {
	d := e.config.ReadLengthDistribution
	if d.Dist == "" || d.Dist == types.ReadLengthFixed {
		return max(e.config.ReadLength, 1)
	}
	n := 2 * max(e.config.ReadLength, d.N50, d.Median, int(d.Mean+2*d.SD), 1)
	if d.Max > 0 && n > d.Max {
		n = d.Max
	}
	return n
}

// headerFormatter returns the per-read header function for the platform.
func (e *Engine) headerFormatter() (func(r *rand.Rand, rec *seqread.Record) (name, comment string), error) {
	if e.isONT() {
		return nanopore.NewHeaderFormatter(e.ontRunInfo(), e.config.ReadsPerBatch, e.slotBases()), nil
	}
	sample, lane := e.config.SampleName, e.config.Lane
	// Reads sampled by a FragmentSampler (a reference genome) keep the
	// per-index tile coordinates: their starts are not unique, and read names
	// must be. Their origin is in the truth table.
	byStart := e.fragments == nil
	return func(_ *rand.Rand, rec *seqread.Record) (string, string) {
		x, y := rec.Index*1000, rec.Index*2000
		if rec.Positioned && byStart {
			x, y = int(rec.RefStart%100000), int((rec.RefStart/100000)%100000)
		}
		return fmt.Sprintf("%s:%d:HXXXXXX:%d:1101:%d:%d", sample, lane, rec.Batch, x, y),
			fmt.Sprintf("%d:N:0:1", rec.Mate)
	}, nil
}

// readTags returns the read-ID fingerprint fields appended to every header
// comment ("\tXC:Z:<customer>\tXD:Z:<dataset>" without the leading tab), or
// "" when read-ID fingerprinting is off.
func (e *Engine) readTags() string {
	fp := e.config.Fingerprinting
	if !fp.Enable {
		return ""
	}
	var cust, ds string
	if fp.CustomerFingerprint {
		cust = e.customerID
	}
	if fp.DatasetFingerprint {
		ds = e.config.DatasetID
	}
	return strings.Join(fingerprinting.ReadTags(cust, ds), "\t")
}

// maxPhred returns the highest Phred score encoded in FASTQ output.
func (e *Engine) maxPhred() byte {
	if e.isONT() {
		return maxONTPhred
	}
	return maxIlluminaPhred
}

// placeholderRE matches one {name} file-pattern placeholder.
var placeholderRE = regexp.MustCompile(`\{([^{}]*)\}`)

// filePattern returns the configured or platform default file pattern, or
// "" for the Illumina default names.
func (e *Engine) filePattern() string {
	if p := e.config.Output.FilePattern; p != "" {
		return p
	}
	if e.isONT() {
		return defaultONTFilePattern
	}
	return ""
}

// expandPattern substitutes the file-pattern placeholders for batch and mate.
func (e *Engine) expandPattern(pattern string, batch, mate int) (string, error) {
	c := e.config
	var info nanopore.RunInfo
	if e.isONT() {
		info = e.ontRunInfo()
	}
	var err error
	name := placeholderRE.ReplaceAllStringFunc(pattern, func(m string) string {
		key := m[1 : len(m)-1]
		var v string
		switch key {
		case "sample_name":
			v = c.SampleName
		case "flowcell_id", "flow_cell_id":
			v = info.FlowCellID
			if !e.isONT() {
				v = c.FlowCellID
			}
		case "run_id":
			v = info.ShortRunID()
			if !e.isONT() {
				v = c.RunID
			}
		case "file_number":
			v = strconv.Itoa(batch - 1)
		case "barcode":
			v = defaultBarcode
		case "lane":
			v = strconv.Itoa(c.Lane)
		case "mate":
			v = strconv.Itoa(mate)
		default:
			err = errors.Join(err, fmt.Errorf("unknown placeholder %s (want {sample_name}, {flowcell_id}, {run_id}, {file_number}, {barcode}, {lane} or {mate})", m))
			return m
		}
		if v == "" {
			err = errors.Join(err, fmt.Errorf("placeholder %s has no value for this configuration", m))
		}
		return v
	})
	if err != nil {
		return "", fmt.Errorf("output file_pattern %q: %w", pattern, err)
	}
	if strings.ContainsAny(name, "{}") {
		return "", fmt.Errorf("output file_pattern %q: unbalanced braces", pattern)
	}
	return name, nil
}

// checkFilePattern validates the file pattern: known placeholders, a plain
// file name, and distinct names for every batch and mate.
func (e *Engine) checkFilePattern() error {
	pattern := e.filePattern()
	if pattern == "" {
		return nil
	}
	name, err := e.expandPattern(pattern, 1, 1)
	if err != nil {
		return err
	}
	if name == "" || name == "." || name == ".." || filepath.Base(name) != name {
		return fmt.Errorf("output file_pattern %q must give a file name without directories, got %q", pattern, name)
	}
	if e.multiFile() && !strings.Contains(pattern, "{file_number}") {
		return fmt.Errorf("output file_pattern %q needs {file_number}: the run writes more than one file (total_batches %d, total_reads %d, total_bases %d)", pattern, e.config.TotalBatches, e.config.TotalReads, e.config.TotalBases)
	}
	if e.config.PairedEnd && !strings.Contains(pattern, "{mate}") {
		return fmt.Errorf("output file_pattern %q needs {mate} for paired-end output", pattern)
	}
	if !e.gzipOutput() && strings.HasSuffix(name, ".gz") {
		return fmt.Errorf("output file_pattern %q ends in .gz but compression is not gzip", pattern)
	}
	return nil
}

// outputPath returns the FASTQ file for batch and mate: the file the
// pipeline writes, which for compress-after-embed output is the plain file
// that finalize compresses.
func (e *Engine) outputPath(batch, mate int) string {
	c := e.config
	var name string
	if pattern := e.filePattern(); pattern != "" {
		name, _ = e.expandPattern(pattern, batch, mate) // validated by checkOutput
	} else if c.PairedEnd {
		name = fmt.Sprintf("%s_S1_L%03d_R%d_%03d.fastq", c.SampleName, c.Lane, mate, batch)
	} else {
		name = fmt.Sprintf("%s_S1_L%03d_%03d.fastq", c.SampleName, c.Lane, batch)
	}
	if e.gzipOutput() {
		name = strings.TrimSuffix(name, ".gz")
		if !e.compressAfterEmbed() {
			name += ".gz"
		}
	}
	return filepath.Join(c.OutputDir, name)
}

// wrapOutput wraps an output stream before records are written to it.
func (e *Engine) wrapOutput(w io.Writer) (io.WriteCloser, error) {
	if e.gzipOutput() && !e.compressAfterEmbed() {
		return gzip.NewWriter(w), nil
	}
	return nopWriteCloser{w}, nil
}

// finalize finishes a closed output file and returns its final path.
func (e *Engine) finalize(path string) (string, error) {
	if e.qualityFingerprintingEnabled() {
		if err := e.embedFileFingerprint(path); err != nil {
			return "", err
		}
	}
	if !e.gzipOutput() {
		return path, nil
	}
	final := gzName(path)
	if e.compressAfterEmbed() {
		return e.gzipFile(path, final)
	}
	if final != path { // compressed while written, under a caller-chosen plain name
		if err := os.Rename(path, final); err != nil {
			os.Remove(path) // #nosec G104 -- best effort cleanup; the error is returned
			return "", fmt.Errorf("renaming %s: %w", filepath.Base(path), err)
		}
	}
	return final, nil
}

// gzName returns name with a ".gz" suffix.
func gzName(name string) string {
	if strings.HasSuffix(name, ".gz") {
		return name
	}
	return name + ".gz"
}

// gzipFile compresses the plain file src to dst (which may equal src) and
// removes src. On failure both are removed, so a run never leaves a file
// that silently lacks the requested compression.
func (e *Engine) gzipFile(src, dst string) (string, error) {
	tmp := dst + ".tmp"
	fail := func(err error) (string, error) {
		os.Remove(tmp) // #nosec G104 -- best effort cleanup; the error is returned
		os.Remove(src) // #nosec G104
		return "", fmt.Errorf("compressing %s: %w", filepath.Base(src), err)
	}
	in, err := os.Open(src) // #nosec G304 -- src is an output file this engine wrote
	if err != nil {
		return fail(err)
	}
	out, err := security.SecureFileCreate(tmp, e.config.OutputDir)
	if err != nil {
		in.Close() // #nosec G104 -- read-only
		return fail(err)
	}
	gw := gzip.NewWriter(out)
	_, err = io.Copy(gw, in)
	err = errors.Join(err, gw.Close(), out.Close(), in.Close())
	if err != nil {
		return fail(err)
	}
	if err := os.Rename(tmp, dst); err != nil {
		return fail(err)
	}
	if src != dst {
		if err := os.Remove(src); err != nil {
			return "", fmt.Errorf("removing uncompressed %s: %w", filepath.Base(src), err)
		}
	}
	return dst, nil
}

// observer returns the run's read observer, or nil for none.
func (e *Engine) observer() (ReadObserver, error) {
	if !e.isONT() {
		return nil, nil
	}
	explicit := e.config.Output.SummaryFile
	outDir := e.config.OutputDir
	open := func(rec *seqread.Record) (io.WriteCloser, error) {
		name := explicit
		if name == "" {
			if rec.File == "" {
				return nil, nil // streamed output: no summary unless requested
			}
			name = nanopore.SummaryFileName
		}
		return security.SecureFileCreate(filepath.Join(outDir, name), outDir)
	}
	return nanopore.NewSummaryWriter(e.ontRunInfo(), e.maxPhred(), open, e.summaryFileName), nil
}

// summaryFileName is the filename_fastq value of a read written to file.
func (e *Engine) summaryFileName(file string) string {
	if file == "" {
		return ""
	}
	if e.gzipOutput() {
		return gzName(filepath.Base(file))
	}
	return filepath.Base(file)
}

// summaryObserver returns a summary writer for one output unit of Generate,
// writing its rows (with a header line) to part; Generate concatenates the
// parts in unit order, so the summary does not depend on scheduling.
func (e *Engine) summaryObserver(part string) ReadObserver {
	dir := e.config.OutputDir
	open := func(*seqread.Record) (io.WriteCloser, error) { return security.SecureFileCreate(part, dir) }
	return nanopore.NewSummaryWriter(e.ontRunInfo(), e.maxPhred(), open, e.summaryFileName)
}

type nopWriteCloser struct{ io.Writer }

func (nopWriteCloser) Close() error { return nil }
