package transcriptome

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"os"
)

// openInput opens path for reading, transparently decompressing gzip
// (detected from the magic bytes, not the file name).
func openInput(path string) (io.ReadCloser, error) {
	f, err := os.Open(path) // #nosec G304 -- user-supplied input file
	if err != nil {
		return nil, err
	}
	br := bufio.NewReaderSize(f, 1<<20)
	magic, _ := br.Peek(2)
	if len(magic) == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
		gz, err := gzip.NewReader(br)
		if err != nil {
			f.Close() // #nosec G104 -- already failing
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		return &gzReadCloser{gz: gz, f: f}, nil
	}
	return &plainReadCloser{Reader: br, f: f}, nil
}

type gzReadCloser struct {
	gz *gzip.Reader
	f  *os.File
}

func (g *gzReadCloser) Read(p []byte) (int, error) { return g.gz.Read(p) }
func (g *gzReadCloser) Close() error               { return errors.Join(g.gz.Close(), g.f.Close()) }

type plainReadCloser struct {
	*bufio.Reader
	f *os.File
}

func (p *plainReadCloser) Close() error { return p.f.Close() }

// lineReader returns lines of any length without the trailing newline (and
// carriage return).
type lineReader struct {
	br   *bufio.Reader
	buf  []byte
	line int
}

func newLineReader(r io.Reader) *lineReader {
	return &lineReader{br: bufio.NewReaderSize(r, 1<<20)}
}

// next returns the next line; the slice is valid until the following call.
// It returns io.EOF after the last line.
func (l *lineReader) next() ([]byte, error) {
	l.buf = l.buf[:0]
	for {
		chunk, err := l.br.ReadSlice('\n')
		l.buf = append(l.buf, chunk...)
		if err == bufio.ErrBufferFull {
			continue
		}
		if err != nil && err != io.EOF {
			return nil, err
		}
		if err == io.EOF && len(l.buf) == 0 {
			return nil, io.EOF
		}
		l.line++
		b := bytes.TrimSuffix(l.buf, []byte{'\n'})
		return bytes.TrimSuffix(b, []byte{'\r'}), nil
	}
}

// readFASTA calls fn for each record of r with its header line (without
// '>') and its sequence in upper case. seq is owned by fn.
func readFASTA(r io.Reader, fn func(header string, seq []byte) error) error {
	lr := newLineReader(r)
	var header string
	var seq []byte
	have := false
	flush := func() error {
		if !have {
			return nil
		}
		return fn(header, seq)
	}
	for {
		line, err := lr.next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' {
			if err := flush(); err != nil {
				return err
			}
			header, seq, have = string(line[1:]), nil, true
			continue
		}
		if line[0] == ';' {
			continue
		}
		if !have {
			return fmt.Errorf("line %d: sequence before the first FASTA header", lr.line)
		}
		for _, c := range line {
			switch {
			case c >= 'a' && c <= 'z':
				seq = append(seq, c-'a'+'A')
			case c == ' ' || c == '\t':
			default:
				seq = append(seq, c)
			}
		}
	}
	return flush()
}

// complement maps a base to its complement (IUPAC-aware; unknown bytes map
// to N).
var complement = func() (t [256]byte) {
	for i := range t {
		t[i] = 'N'
	}
	pairs := "ACGTUMRWSYKVHDBN"
	comps := "TGCAAKYWSRMBDHVN"
	for i := range pairs {
		t[pairs[i]] = comps[i]
		t[pairs[i]+'a'-'A'] = comps[i] + 'a' - 'A'
	}
	return t
}()

// ReverseComplement returns the reverse complement of s as a new slice.
func ReverseComplement(s []byte) []byte {
	out := make([]byte, len(s))
	for i, c := range s {
		out[len(s)-1-i] = complement[c]
	}
	return out
}
