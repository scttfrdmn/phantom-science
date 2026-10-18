package slow5

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// DefaultReadsPerFile is the number of reads per output file of a chunked
// run. MinKNOW on a PromethION writes about 4,000 reads per pod5 file, which
// gives the ~400-file shape of a flow cell that downstream task granularity
// is built on.
const DefaultReadsPerFile = 4000

// ChunkOptions configures a ChunkedWriter.
type ChunkOptions struct {
	// ReadsPerFile is the number of records per file; 0 means
	// DefaultReadsPerFile.
	ReadsPerFile int
	// BLOW5 selects the compression of every file.
	BLOW5 BLOW5Options
	// FlowCellID and RunID name the files. Empty values are taken from the
	// flow_cell_id and run_id attributes of read group 0.
	FlowCellID, RunID string
	// Overwrite allows replacing existing files. By default a file that
	// already exists is an error, so a rerun into a used directory cannot
	// leave a mix of old and new chunks.
	Overwrite bool
}

// ChunkedWriter writes a run as a series of BLOW5 files in one directory,
// opening a new file every ReadsPerFile records. Files are named
// <flowcell>_<runid>_<n>.blow5 with n counting from 0; every file carries the
// full header, so each is a valid stand-alone BLOW5 file. Output is
// deterministic: the same header, options and records give byte-identical
// files.
type ChunkedWriter struct {
	dir     string
	h       Header
	opts    ChunkOptions
	prefix  string
	f       *os.File
	w       *BLOW5Writer
	inFile  int
	files   []string
	records int
	err     error
	closed  bool
}

// NewChunkedWriter validates the header and options and creates dir if
// needed. No file is created until the first Write.
func NewChunkedWriter(dir string, h Header, opts ChunkOptions) (*ChunkedWriter, error) {
	h = cloneHeader(h)
	if err := h.validate(); err != nil {
		return nil, err
	}
	if opts.ReadsPerFile < 0 {
		return nil, formatErr("ReadsPerFile %d is negative", opts.ReadsPerFile)
	}
	if opts.ReadsPerFile == 0 {
		opts.ReadsPerFile = DefaultReadsPerFile
	}
	if opts.FlowCellID == "" {
		opts.FlowCellID = h.ReadGroups[0]["flow_cell_id"]
	}
	if opts.RunID == "" {
		opts.RunID = h.ReadGroups[0]["run_id"]
	}
	for name, v := range map[string]string{"flow cell ID": opts.FlowCellID, "run ID": opts.RunID} {
		if !safeNamePart(v) {
			return nil, formatErr("%s %q cannot be used in a file name (set it, and use letters, digits, '-' or '.')", name, v)
		}
	}
	// Check the compression options now rather than at the first file.
	if _, err := NewBLOW5Writer(discard{}, h, opts.BLOW5); err != nil {
		return nil, err
	}
	if err := os.MkdirAll(dir, 0o750); err != nil {
		return nil, err
	}
	return &ChunkedWriter{
		dir:    dir,
		h:      h,
		opts:   opts,
		prefix: opts.FlowCellID + "_" + opts.RunID + "_",
	}, nil
}

type discard struct{}

func (discard) Write(p []byte) (int, error) { return len(p), nil }

func safeNamePart(s string) bool {
	if s == "" || s == "." || s == ".." {
		return false
	}
	for i := 0; i < len(s); i++ {
		c := s[i]
		if !(c == '-' || c == '.' || c >= '0' && c <= '9' || c >= 'a' && c <= 'z' || c >= 'A' && c <= 'Z') {
			return false
		}
	}
	return !strings.HasPrefix(s, ".")
}

// ChunkFileName returns the name of chunk n (0-based) for a flow cell and
// run: <flowcell>_<runid>_<n>.blow5.
func ChunkFileName(flowCellID, runID string, n int) string {
	return fmt.Sprintf("%s_%s_%d.blow5", flowCellID, runID, n)
}

// Write appends a record, starting a new file when the current one is full.
func (cw *ChunkedWriter) Write(r *Record) error {
	if cw.err != nil {
		return cw.err
	}
	if cw.closed {
		return formatErr("write after Close")
	}
	// Validate before opening a file, so a bad first record of a chunk does
	// not leave an empty file behind.
	if err := cw.h.validateRecord(r); err != nil {
		return err
	}
	if cw.w == nil || cw.inFile == cw.opts.ReadsPerFile {
		if err := cw.finishFile(); err != nil {
			return cw.fail(err)
		}
		if err := cw.openFile(); err != nil {
			return cw.fail(err)
		}
	}
	if err := cw.w.Write(r); err != nil {
		if errors.Is(err, ErrFormat) {
			return err // bad record (e.g. aux value); file still consistent
		}
		return cw.fail(err)
	}
	cw.inFile++
	cw.records++
	return nil
}

func (cw *ChunkedWriter) fail(err error) error {
	cw.err = err
	return err
}

func (cw *ChunkedWriter) openFile() error {
	name := filepath.Join(cw.dir, ChunkFileName(cw.opts.FlowCellID, cw.opts.RunID, len(cw.files)))
	flags := os.O_WRONLY | os.O_CREATE | os.O_EXCL
	if cw.opts.Overwrite {
		flags = os.O_WRONLY | os.O_CREATE | os.O_TRUNC
	}
	f, err := os.OpenFile(filepath.Clean(name), flags, 0o600)
	if err != nil {
		return err
	}
	w, err := NewBLOW5Writer(f, cw.h, cw.opts.BLOW5)
	if err != nil {
		_ = f.Close()
		return err
	}
	cw.f, cw.w, cw.inFile = f, w, 0
	cw.files = append(cw.files, name)
	return nil
}

func (cw *ChunkedWriter) finishFile() error {
	if cw.w == nil {
		return nil
	}
	err := cw.w.Close()
	if cerr := cw.f.Close(); err == nil {
		err = cerr
	}
	cw.w, cw.f = nil, nil
	return err
}

// Close finishes the last file. A run with no records produces no files.
func (cw *ChunkedWriter) Close() error {
	if cw.closed {
		return cw.err
	}
	cw.closed = true
	if err := cw.finishFile(); err != nil && cw.err == nil {
		cw.err = err
	}
	return cw.err
}

// Files returns the paths of the files written so far, in order.
func (cw *ChunkedWriter) Files() []string { return append([]string(nil), cw.files...) }

// Records returns the number of records written.
func (cw *ChunkedWriter) Records() int { return cw.records }
