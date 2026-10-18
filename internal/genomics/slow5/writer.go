package slow5

import (
	"bufio"
	"bytes"
	"compress/zlib"
	"encoding/binary"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/klauspost/compress/zstd"
)

// Writer writes SLOW5 or BLOW5 records. Close must be called to finish the
// file (for BLOW5 it writes the end-of-file marker); it does not close the
// underlying io.Writer.
type Writer interface {
	Write(r *Record) error
	Close() error
}

// BLOW5Options selects the BLOW5 compression methods.
type BLOW5Options struct {
	RecordCompression RecordCompression
	SignalCompression SignalCompression
}

// headerText renders the part of the header shared by SLOW5 ASCII and
// BLOW5: the data-header lines (one "@key\tv0\tv1..." line per attribute,
// sorted by key, "." for a group lacking the key), the type line and the
// column-name line (spec, SLOW5 Header; BLOW5 Header: the BLOW5 copy omits
// the #slow5_version and #num_read_groups lines).
func headerText(h *Header) []byte {
	var b bytes.Buffer
	for _, k := range h.attrKeys() {
		b.WriteByte('@')
		b.WriteString(k)
		for _, g := range h.ReadGroups {
			b.WriteByte('\t')
			if v, ok := g[k]; ok {
				b.WriteString(v)
			} else {
				b.WriteString(missingASCII)
			}
		}
		b.WriteByte('\n')
	}
	b.WriteString(primaryTypeLine)
	for _, f := range h.Aux {
		b.WriteByte('\t')
		b.WriteString(f.typeString())
	}
	b.WriteByte('\n')
	b.WriteString(primaryNameLine)
	for _, f := range h.Aux {
		b.WriteByte('\t')
		b.WriteString(f.Name)
	}
	b.WriteByte('\n')
	return b.Bytes()
}

// cloneHeader copies h so later changes by the caller do not affect a
// writer.
func cloneHeader(h Header) Header {
	c := Header{Version: h.Version}
	for _, g := range h.ReadGroups {
		ng := make(ReadGroup, len(g))
		for k, v := range g {
			ng[k] = v
		}
		c.ReadGroups = append(c.ReadGroups, ng)
	}
	for _, f := range h.Aux {
		f.EnumLabels = append([]string(nil), f.EnumLabels...)
		c.Aux = append(c.Aux, f)
	}
	return c
}

// ---------------------------------------------------------------- BLOW5

// BLOW5Writer writes a BLOW5 file.
type BLOW5Writer struct {
	w    *bufio.Writer
	h    Header
	opts BLOW5Options
	zenc *zstd.Encoder
	rec  []byte // reused record buffer
	out  []byte // reused compressed-record buffer
	err  error
	done bool
}

// NewBLOW5Writer validates h and opts and writes the BLOW5 header to w.
func NewBLOW5Writer(w io.Writer, h Header, opts BLOW5Options) (*BLOW5Writer, error) {
	h = cloneHeader(h)
	if err := h.validate(); err != nil {
		return nil, err
	}
	switch opts.RecordCompression {
	case RecordNone, RecordZlib, RecordZstd:
	default:
		return nil, formatErr("unsupported record compression %s", opts.RecordCompression)
	}
	switch opts.SignalCompression {
	case SignalNone, SignalSvbZd:
	default:
		return nil, formatErr("unsupported signal compression %s", opts.SignalCompression)
	}
	// slow5lib slow5_press_version_bump: signal compression and record
	// compression other than none/zlib need v0.2.0.
	if h.version().Less(Version020) &&
		(opts.SignalCompression != SignalNone || opts.RecordCompression == RecordZstd) {
		h.Version = Version020
	}
	bw := &BLOW5Writer{w: bufio.NewWriterSize(w, 1<<20), h: h, opts: opts}
	if opts.RecordCompression == RecordZstd {
		// Single-segment frames always carry Frame_Content_Size, which
		// slow5lib requires (ptr_depress_zstd calls
		// ZSTD_getFrameContentSize and fails on "unknown"). klauspost
		// omits it for non-single-segment frames under 256 bytes. Level
		// SpeedFastest is the closest to slow5lib's
		// SLOW5_ZSTD_COMPRESS_LEVEL 1.
		enc, err := zstd.NewWriter(nil,
			zstd.WithEncoderLevel(zstd.SpeedFastest),
			zstd.WithEncoderConcurrency(1),
			zstd.WithSingleSegment(true))
		if err != nil {
			return nil, err
		}
		bw.zenc = enc
	}
	if _, err := bw.w.Write(blow5HeaderBytes(&bw.h, opts)); err != nil {
		return nil, err
	}
	return bw, nil
}

// blow5HeaderBytes renders the BLOW5 header (spec, BLOW5 Header, Table 6):
//
//	offset 0  char[6]  magic "BLOW5\1"
//	offset 6  uint8    major, 7 minor, 8 patch version
//	offset 9  uint8    record compression
//	offset 10 uint32   number of read groups
//	offset 14 uint8    signal compression (v0.2.0 on)
//	offset 15 ...      reserved, zero up to offset 64
//	offset 64 uint32   size of the header text that follows
//	offset 68 char[]   header text, no NUL
//
// Table 6 lists "14 / 50 / Reserved" after the signal compression byte at
// offset 14, but the text says "49 bytes are reserved"; the reserved bytes
// are 15..63 (49 bytes), as in slow5lib slow5_hdr_to_mem, which zero-fills
// up to SLOW5_BINARY_HDR_SIZE_OFFSET (64).
func blow5HeaderBytes(h *Header, opts BLOW5Options) []byte {
	v := h.version()
	b := make([]byte, headerSizeAt+4, 1024)
	copy(b, blow5Magic)
	b[6], b[7], b[8] = v.Major, v.Minor, v.Patch
	b[9] = byte(opts.RecordCompression)
	binary.LittleEndian.PutUint32(b[10:], uint32(len(h.ReadGroups)))
	b[14] = byte(opts.SignalCompression)
	text := headerText(h)
	binary.LittleEndian.PutUint32(b[headerSizeAt:], uint32(len(text)))
	return append(b, text...)
}

// Write encodes and writes one record. The record is not modified.
func (bw *BLOW5Writer) Write(r *Record) error {
	if bw.err != nil {
		return bw.err
	}
	if bw.done {
		return formatErr("write after Close")
	}
	rec, err := encodeBLOW5Record(bw.rec[:0], &bw.h, r, bw.opts.SignalCompression)
	if err != nil {
		return err // a bad record does not poison the writer
	}
	bw.rec = rec
	payload := rec
	switch bw.opts.RecordCompression {
	case RecordZstd:
		bw.out = bw.zenc.EncodeAll(rec, bw.out[:0])
		payload = bw.out
	case RecordZlib:
		// One complete zlib stream (RFC 1950, deflateInit2 with
		// MAX_WBITS) per record, like slow5lib ptr_compress_zlib with
		// Z_FINISH.
		var zb bytes.Buffer
		zw, _ := zlib.NewWriterLevel(&zb, zlib.DefaultCompression)
		if _, err := zw.Write(rec); err != nil {
			return bw.fail(err)
		}
		if err := zw.Close(); err != nil {
			return bw.fail(err)
		}
		payload = zb.Bytes()
	}
	// "Each BLOW5 record is preceded by the size of the upcoming BLOW5
	// record in bytes (the size of the compressed record if compressed),
	// which is an 8-byte uint64_t" (spec, BLOW5 Data).
	var size [8]byte
	binary.LittleEndian.PutUint64(size[:], uint64(len(payload)))
	if _, err := bw.w.Write(size[:]); err != nil {
		return bw.fail(err)
	}
	if _, err := bw.w.Write(payload); err != nil {
		return bw.fail(err)
	}
	return nil
}

func (bw *BLOW5Writer) fail(err error) error {
	bw.err = err
	return err
}

// Close writes the end-of-file marker "5WOLB" (spec, BLOW5 Footer) and
// flushes. It does not close the underlying writer.
func (bw *BLOW5Writer) Close() error {
	if bw.done {
		return bw.err
	}
	bw.done = true
	if bw.zenc != nil {
		_ = bw.zenc.Close()
	}
	if bw.err != nil {
		return bw.err
	}
	if _, err := bw.w.Write(blow5EOF); err != nil {
		return bw.fail(err)
	}
	if err := bw.w.Flush(); err != nil {
		return bw.fail(err)
	}
	return nil
}

// encodeBLOW5Record appends the uncompressed binary record (spec, BLOW5
// Data, Table 7) to b:
//
//	uint16   read_id length, then read_id bytes (no NUL)
//	uint32   read_group
//	double   digitisation, offset, range, sampling_rate
//	uint64   len_raw_signal (svb-zd: compressed byte count, see svbzd.go)
//	int16[]  raw_signal (or the svb-zd blob)
//	aux fields in header order: primitives as Table 8 (sentinel when
//	missing); arrays as a uint64 element count then the elements (count 0
//	when missing)
func encodeBLOW5Record(b []byte, h *Header, r *Record, sc SignalCompression) ([]byte, error) {
	if err := h.validateRecord(r); err != nil {
		return nil, err
	}
	le := binary.LittleEndian
	b = le.AppendUint16(b, uint16(len(r.ReadID)))
	b = append(b, r.ReadID...)
	b = le.AppendUint32(b, r.ReadGroup)
	for _, x := range [...]float64{r.Digitisation, r.Offset, r.Range, r.SamplingRate} {
		b = le.AppendUint64(b, math.Float64bits(x))
	}
	if sc == SignalSvbZd {
		if uint64(len(r.RawSignal)) > svbZdMaxLen {
			return nil, formatErr("read %s: %d samples exceed the svb-zd limit", r.ReadID, len(r.RawSignal))
		}
		at := len(b)
		b = le.AppendUint64(b, 0)
		b = svbZdEncode(b, r.RawSignal)
		le.PutUint64(b[at:], uint64(len(b)-at-8))
	} else {
		b = le.AppendUint64(b, uint64(len(r.RawSignal)))
		for _, s := range r.RawSignal {
			b = le.AppendUint16(b, uint16(s))
		}
	}
	for _, f := range h.Aux {
		v, present := r.Aux[f.Name]
		arr := f.Type.IsArray()
		if !present {
			if arr {
				b = le.AppendUint64(b, 0)
			} else {
				b = appendSentinel(b, f.Type)
			}
			continue
		}
		elems, err := checkAuxValue(f, v)
		if err != nil {
			return nil, err
		}
		if f.Type == AuxString {
			s := v.(string)
			b = le.AppendUint64(b, uint64(len(s)))
			b = append(b, s...)
			continue
		}
		if arr {
			b = le.AppendUint64(b, uint64(len(elems)))
		}
		for _, e := range elems {
			b = appendPrim(b, e)
		}
	}
	return b, nil
}

// ---------------------------------------------------------------- SLOW5 ASCII

// SLOW5Writer writes a SLOW5 ASCII file.
type SLOW5Writer struct {
	w    *bufio.Writer
	h    Header
	line []byte
	err  error
	done bool
}

// NewSLOW5Writer validates h and writes the SLOW5 ASCII header to w (spec,
// SLOW5 Header: "#slow5_version\t<v>", "#num_read_groups\t<n>", the data
// header, the type line and the column-name line).
func NewSLOW5Writer(w io.Writer, h Header) (*SLOW5Writer, error) {
	h = cloneHeader(h)
	if err := h.validate(); err != nil {
		return nil, err
	}
	sw := &SLOW5Writer{w: bufio.NewWriterSize(w, 1<<20), h: h}
	var b bytes.Buffer
	b.WriteString("#slow5_version\t" + h.version().String() + "\n")
	b.WriteString("#num_read_groups\t" + strconv.Itoa(len(h.ReadGroups)) + "\n")
	b.Write(headerText(&h))
	if _, err := sw.w.Write(b.Bytes()); err != nil {
		return nil, err
	}
	return sw, nil
}

// Write formats and writes one record as a tab-separated line (spec, SLOW5
// Data), with the raw signal comma-separated.
func (sw *SLOW5Writer) Write(r *Record) error {
	if sw.err != nil {
		return sw.err
	}
	if sw.done {
		return formatErr("write after Close")
	}
	line, err := encodeSLOW5Record(sw.line[:0], &sw.h, r)
	if err != nil {
		return err
	}
	sw.line = line
	if _, err := sw.w.Write(line); err != nil {
		sw.err = err
		return err
	}
	return nil
}

// Close flushes the output. SLOW5 ASCII has no end-of-file marker.
func (sw *SLOW5Writer) Close() error {
	if sw.done {
		return sw.err
	}
	sw.done = true
	if sw.err != nil {
		return sw.err
	}
	if err := sw.w.Flush(); err != nil {
		sw.err = err
	}
	return sw.err
}

func encodeSLOW5Record(b []byte, h *Header, r *Record) ([]byte, error) {
	if err := h.validateRecord(r); err != nil {
		return nil, err
	}
	b = append(b, r.ReadID...)
	b = append(b, '\t')
	b = strconv.AppendUint(b, uint64(r.ReadGroup), 10)
	for _, x := range [...]float64{r.Digitisation, r.Offset, r.Range, r.SamplingRate} {
		b = append(b, '\t')
		b = append(b, formatDouble(x, 64)...)
	}
	b = append(b, '\t')
	b = strconv.AppendUint(b, uint64(len(r.RawSignal)), 10)
	b = append(b, '\t')
	// slow5lib writes an empty field for an empty signal.
	for i, s := range r.RawSignal {
		if i > 0 {
			b = append(b, ',')
		}
		b = strconv.AppendInt(b, int64(s), 10)
	}
	for _, f := range h.Aux {
		b = append(b, '\t')
		v, present := r.Aux[f.Name]
		if !present {
			b = append(b, missingASCII...)
			continue
		}
		elems, err := checkAuxValue(f, v)
		if err != nil {
			return nil, err
		}
		if f.Type == AuxString {
			b = append(b, v.(string)...)
			continue
		}
		strs := make([]string, len(elems))
		for i, e := range elems {
			strs[i] = formatPrim(e)
		}
		b = append(b, strings.Join(strs, ",")...)
	}
	return append(b, '\n'), nil
}
