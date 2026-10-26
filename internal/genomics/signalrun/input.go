package signalrun

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"crypto/sha1" // #nosec G505 -- UUIDv5 name hashing (RFC 4122), not security
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"fmt"
	"hash"
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
)

// Namespace is the UUIDv5 namespace of mapped read IDs: the
// RFC 4122 URL namespace applied to the project URL.
var Namespace = UUIDv5(namespaceURL, "https://github.com/scttfrdmn/phantom-science/signal")

// namespaceURL is the RFC 4122 URL namespace, 6ba7b811-9dad-11d1-80b4-00c04fd430c8.
var namespaceURL = [16]byte{0x6b, 0xa7, 0xb8, 0x11, 0x9d, 0xad, 0x11, 0xd1, 0x80, 0xb4, 0x00, 0xc0, 0x4f, 0xd4, 0x30, 0xc8}

// UUIDv5 is the RFC 4122 name-based (SHA-1) UUID of name in namespace ns.
func UUIDv5(ns [16]byte, name string) [16]byte {
	h := sha1.New() // #nosec G401 -- RFC 4122 version 5 is defined over SHA-1
	h.Write(ns[:])
	h.Write([]byte(name))
	var u [16]byte
	copy(u[:], h.Sum(nil))
	u[6] = u[6]&0x0f | 0x50 // version 5
	u[8] = u[8]&0x3f | 0x80 // RFC 4122 variant
	return u
}

func FormatUUID(u [16]byte) string {
	var b [36]byte
	hex.Encode(b[0:8], u[0:4])
	b[8] = '-'
	hex.Encode(b[9:13], u[4:6])
	b[13] = '-'
	hex.Encode(b[14:18], u[6:8])
	b[18] = '-'
	hex.Encode(b[19:23], u[8:10])
	b[23] = '-'
	hex.Encode(b[24:36], u[10:16])
	return string(b[:])
}

// IsUUID reports whether s has the 8-4-4-4-12 hex form of a UUID.
func IsUUID(s string) bool {
	if len(s) != 36 {
		return false
	}
	for i := 0; i < 36; i++ {
		c := s[i]
		if i == 8 || i == 13 || i == 18 || i == 23 {
			if c != '-' {
				return false
			}
			continue
		}
		if !(c >= '0' && c <= '9' || c >= 'a' && c <= 'f' || c >= 'A' && c <= 'F') {
			return false
		}
	}
	return true
}

// ReadID returns the BLOW5 read ID for an input read ID: the ID itself
// (lowercased) if it is a UUID, otherwise its UUIDv5. mapped reports whether
// the result differs from name.
func ReadID(name string) (id string, mapped bool) {
	if IsUUID(name) {
		id = strings.ToLower(name)
		return id, id != name
	}
	return FormatUUID(UUIDv5(Namespace, name)), true
}

// SeqFile reads FASTQ (4-line records) or FASTA (multi-line) records from a
// plain or gzip file, hashing the file's bytes as they are read.
type SeqFile struct {
	path  string
	f     *os.File
	raw   io.Reader // file bytes (through the hash)
	br    *bufio.Reader
	gz    *gzip.Reader
	h     hash.Hash
	n     int64
	fasta bool
	began bool
	line  int
}

type countWriter struct{ n *int64 }

func (c countWriter) Write(p []byte) (int, error) { *c.n += int64(len(p)); return len(p), nil }

func OpenSeqFile(path string) (*SeqFile, error) {
	f, err := os.Open(path) // #nosec G304 -- user-supplied input path
	if err != nil {
		return nil, err
	}
	s := &SeqFile{path: path, f: f, h: sha256.New()}
	rawBuf := bufio.NewReaderSize(io.TeeReader(f, io.MultiWriter(s.h, countWriter{&s.n})), 1<<20)
	s.raw = rawBuf
	if magic, _ := rawBuf.Peek(2); len(magic) == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
		gz, err := gzip.NewReader(rawBuf)
		if err != nil {
			_ = f.Close()
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		s.gz = gz
		s.br = bufio.NewReaderSize(gz, 1<<20)
	} else {
		s.br = rawBuf
	}
	return s, nil
}

// finish reads the rest of the file (so the hash covers all of it) and
// closes it, returning the size and hex SHA-256.
func (s *SeqFile) Finish() (manifest.Input, error) {
	_, err := io.Copy(io.Discard, s.raw)
	if s.gz != nil {
		err = errors.Join(err, s.gz.Close())
	}
	err = errors.Join(err, s.f.Close())
	return manifest.Input{Name: filepath.Base(s.path), Bytes: s.n, SHA256: hex.EncodeToString(s.h.Sum(nil))}, err
}

func (s *SeqFile) Close() { _ = s.f.Close() }

func (s *SeqFile) readLine() ([]byte, error) {
	b, err := s.br.ReadBytes('\n')
	if err == io.EOF && len(b) > 0 {
		err = nil
	}
	if err != nil {
		return nil, err
	}
	s.line++
	return bytes.TrimRight(b, "\r\n"), nil
}

func (s *SeqFile) errorf(format string, args ...any) error {
	return fmt.Errorf("%s:%d: %s", s.path, s.line, fmt.Sprintf(format, args...))
}

// next returns the next record's ID, header comment and sequence, or io.EOF.
func (s *SeqFile) Next() (name, comment string, seq []byte, err error) {
	var head []byte
	for {
		if head, err = s.readLine(); err != nil {
			return "", "", nil, err
		}
		if len(head) > 0 {
			break
		}
	}
	if !s.began {
		s.began = true
		switch head[0] {
		case '@':
		case '>':
			s.fasta = true
		default:
			return "", "", nil, s.errorf("not FASTQ or FASTA (line starts with %q)", head[0])
		}
	}
	marker := byte('@')
	if s.fasta {
		marker = '>'
	}
	if head[0] != marker {
		return "", "", nil, s.errorf("expected a %q header line", marker)
	}
	name, comment, _ = strings.Cut(string(head[1:]), " ")
	if i := strings.IndexByte(name, '\t'); i >= 0 {
		name, comment = name[:i], name[i+1:]+" "+comment
	}
	if name == "" {
		return "", "", nil, s.errorf("empty read ID")
	}
	if s.fasta {
		for {
			p, perr := s.br.Peek(1)
			if perr != nil || p[0] == '>' {
				break
			}
			l, lerr := s.readLine()
			if lerr != nil {
				return "", "", nil, lerr
			}
			seq = append(seq, l...)
		}
		return name, comment, seq, nil
	}
	if seq, err = s.readLine(); err != nil {
		return "", "", nil, s.errorf("truncated FASTQ record %s", name)
	}
	plus, err := s.readLine()
	if err != nil || len(plus) == 0 || plus[0] != '+' {
		return "", "", nil, s.errorf("FASTQ record %s: missing '+' line", name)
	}
	qual, err := s.readLine()
	if err != nil || len(qual) != len(seq) {
		return "", "", nil, s.errorf("FASTQ record %s: quality length differs from sequence length", name)
	}
	return name, comment, append([]byte(nil), seq...), nil
}

// LoadTruth reads a FASTA/FASTQ of true sequences keyed by record ID.
func LoadTruth(path string) (map[string][]byte, manifest.Input, error) {
	s, err := OpenSeqFile(path)
	if err != nil {
		return nil, manifest.Input{}, err
	}
	defer s.Close()
	truth := map[string][]byte{}
	for {
		name, _, seq, err := s.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, manifest.Input{}, err
		}
		if _, dup := truth[name]; dup {
			return nil, manifest.Input{}, fmt.Errorf("%s: duplicate record %s", path, name)
		}
		truth[name] = seq
	}
	in, err := s.Finish()
	return truth, in, err
}

// DetectModelUnits returns "pa" for a squigulator/nanopolish table (header
// row or level_stdv column) and "standardised" for a two-column ONT
// kmer_models table.
func DetectModelUnits(path string) (string, error) {
	f, err := os.Open(path) // #nosec G304 -- user-supplied model path
	if err != nil {
		return "", fmt.Errorf("pore model: %w", err)
	}
	defer f.Close()
	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return "", fmt.Errorf("pore model %s: %w", path, err)
		}
		defer gz.Close()
		r = gz
	}
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		fields := strings.Fields(sc.Text())
		if len(fields) == 0 || strings.HasPrefix(fields[0], "#") {
			continue
		}
		if strings.EqualFold(fields[0], "kmer") || len(fields) >= 3 {
			return "pa", nil
		}
		return "standardised", nil
	}
	if err := sc.Err(); err != nil {
		return "", fmt.Errorf("pore model %s: %w", path, err)
	}
	return "", fmt.Errorf("pore model %s: no k-mer rows", path)
}

// Source resolves an input record to its output read ID, the sequence to
// synthesise and whether that came from --truth.
func Source(name string, seq []byte, truth map[string][]byte) (id string, mapped bool, src []byte, fromTruth bool) {
	id, mapped = ReadID(name)
	if t, ok := truth[name]; ok {
		return id, mapped, t, true
	}
	if t, ok := truth[id]; ok {
		return id, mapped, t, true
	}
	return id, mapped, seq, false
}
