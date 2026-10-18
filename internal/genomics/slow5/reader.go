package slow5

import (
	"bufio"
	"bytes"
	"compress/zlib"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/klauspost/compress/zstd"
)

// Format is the on-disk flavour of a SLOW5 file.
type Format uint8

// File formats.
const (
	FormatSLOW5 Format = iota + 1 // SLOW5 ASCII (.slow5)
	FormatBLOW5                   // binary (.blow5)
)

func (f Format) String() string {
	switch f {
	case FormatSLOW5:
		return "slow5"
	case FormatBLOW5:
		return "blow5"
	}
	return fmt.Sprintf("format(%d)", uint8(f))
}

// maxDecompressedRecord bounds the memory a single zstd record may claim, so
// a corrupt Frame_Content_Size cannot trigger a huge allocation. 2 GiB is
// a read of about one billion samples.
const maxDecompressedRecord = 1 << 31

// maxReadGroups bounds num_read_groups when reading.
const maxReadGroups = 1 << 20

// Reader reads a SLOW5 ASCII or BLOW5 file, detected from its first bytes.
type Reader struct {
	br     *bufio.Reader
	h      Header
	format Format
	rc     RecordCompression
	sc     SignalCompression
	zdec   *zstd.Decoder
	buf    bytes.Buffer
	n      int // records read
	done   bool
	err    error
}

// NewReader reads and validates the header of r. Use Read to iterate the
// records.
func NewReader(r io.Reader) (*Reader, error) {
	rd := &Reader{br: bufio.NewReaderSize(r, 1<<20)}
	magic, err := rd.br.Peek(len(blow5Magic))
	switch {
	case err == nil && bytes.Equal(magic, blow5Magic):
		rd.format = FormatBLOW5
		err = rd.readBLOW5Header()
	case err == nil && bytes.HasPrefix(magic, []byte("#slow5")):
		rd.format = FormatSLOW5
		err = rd.readSLOW5Header()
	case err == nil || errors.Is(err, io.EOF) || errors.Is(err, bufio.ErrBufferFull):
		err = formatErr("not a SLOW5 or BLOW5 file (no BLOW5 magic or #slow5_version line)")
	}
	if err != nil {
		return nil, err
	}
	return rd, nil
}

// Header returns a copy of the file header.
func (rd *Reader) Header() Header { return cloneHeader(rd.h) }

// Format reports whether the file is SLOW5 ASCII or BLOW5.
func (rd *Reader) Format() Format { return rd.format }

// RecordCompression returns the BLOW5 record compression (RecordNone for
// SLOW5 ASCII).
func (rd *Reader) RecordCompression() RecordCompression { return rd.rc }

// SignalCompression returns the BLOW5 signal compression (SignalNone for
// SLOW5 ASCII and for BLOW5 files older than v0.2.0).
func (rd *Reader) SignalCompression() SignalCompression { return rd.sc }

// Close releases decoder resources. It does not close the underlying reader.
func (rd *Reader) Close() error {
	if rd.zdec != nil {
		rd.zdec.Close()
		rd.zdec = nil
	}
	return nil
}

// Read returns the next record, or io.EOF after the last one. For BLOW5,
// io.EOF is returned only once the "5WOLB" end-of-file marker has been
// found as the last five bytes of the file; a file without it is reported
// as truncated (spec, BLOW5 Footer).
func (rd *Reader) Read() (*Record, error) {
	if rd.err != nil {
		return nil, rd.err
	}
	if rd.done {
		return nil, io.EOF
	}
	var (
		rec *Record
		err error
	)
	if rd.format == FormatBLOW5 {
		rec, err = rd.readBLOW5Record()
	} else {
		rec, err = rd.readSLOW5Record()
	}
	if err == io.EOF {
		rd.done = true
		return nil, io.EOF
	}
	if err != nil {
		rd.err = fmt.Errorf("record %d: %w", rd.n, err)
		return nil, rd.err
	}
	if uint64(rec.ReadGroup) >= uint64(len(rd.h.ReadGroups)) {
		// "0 <= read_group < num_read_groups" (spec, Primary fields).
		rd.err = fmt.Errorf("record %d: %w", rd.n, formatErr("read %s: read_group %d >= num_read_groups %d", rec.ReadID, rec.ReadGroup, len(rd.h.ReadGroups)))
		return nil, rd.err
	}
	rd.n++
	return rec, nil
}

// ReadAll reads a whole SLOW5/BLOW5 stream.
func ReadAll(r io.Reader) (Header, []*Record, error) {
	rd, err := NewReader(r)
	if err != nil {
		return Header{}, nil, err
	}
	defer rd.Close()
	var recs []*Record
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			return rd.Header(), recs, nil
		}
		if err != nil {
			return rd.Header(), recs, err
		}
		recs = append(recs, rec)
	}
}

// ---------------------------------------------------------------- headers

func parseVersion(s string) (Version, error) {
	parts := strings.Split(s, ".")
	if len(parts) != 3 {
		return Version{}, formatErr("bad version %q", s)
	}
	var v [3]uint8
	for i, p := range parts {
		n, err := strconv.ParseUint(p, 10, 8)
		if err != nil {
			return Version{}, formatErr("bad version %q", s)
		}
		v[i] = uint8(n)
	}
	return Version{v[0], v[1], v[2]}, nil
}

// checkVersion rejects files newer than this package supports, as slow5lib
// does (slow5_is_version_compatible against SLOW5_VERSION 0.2.0).
func checkVersion(v Version) error {
	if Version020.Less(v) {
		return formatErr("file version %s is newer than the supported %s", v, Version020)
	}
	return nil
}

func (rd *Reader) readBLOW5Header() error {
	var fixed [headerSizeAt + 4]byte
	if _, err := io.ReadFull(rd.br, fixed[:]); err != nil {
		return formatErr("truncated BLOW5 header: %v", err)
	}
	v := Version{fixed[6], fixed[7], fixed[8]}
	if err := checkVersion(v); err != nil {
		return err
	}
	rd.h.Version = v
	rd.rc = RecordCompression(fixed[9])
	nrg := binary.LittleEndian.Uint32(fixed[10:])
	// The signal compression byte exists from v0.2.0 on (spec, Table 6;
	// slow5lib reads it only if slow5_signal_press_version_cmp >= 0).
	if !v.Less(Version020) {
		rd.sc = SignalCompression(fixed[14])
	}
	switch rd.rc {
	case RecordNone, RecordZlib:
	case RecordZstd:
		dec, err := zstd.NewReader(nil, zstd.WithDecoderConcurrency(1),
			zstd.WithDecoderMaxMemory(maxDecompressedRecord))
		if err != nil {
			return err
		}
		rd.zdec = dec
	default:
		return formatErr("unsupported record compression code %d", uint8(rd.rc))
	}
	switch rd.sc {
	case SignalNone, SignalSvbZd:
	default:
		return formatErr("unsupported signal compression code %d (ex-zd is not supported)", uint8(rd.sc))
	}
	if nrg == 0 {
		return formatErr("num_read_groups is 0")
	}
	size := binary.LittleEndian.Uint32(fixed[headerSizeAt:])
	var text bytes.Buffer
	if n, err := io.CopyN(&text, rd.br, int64(size)); err != nil {
		return formatErr("truncated BLOW5 header text: read %d of %d bytes", n, size)
	}
	t := text.String()
	if !strings.HasSuffix(t, "\n") {
		return formatErr("BLOW5 header text does not end with a newline")
	}
	lines := strings.Split(strings.TrimSuffix(t, "\n"), "\n")
	return rd.parseHeaderLines(lines, nrg)
}

func (rd *Reader) readSLOW5Header() error {
	readLine := func() (string, error) {
		s, err := rd.br.ReadString('\n')
		if err != nil {
			return "", formatErr("truncated SLOW5 header")
		}
		return strings.TrimSuffix(s, "\n"), nil
	}
	l, err := readLine()
	if err != nil {
		return err
	}
	key, val, ok := strings.Cut(l, "\t")
	if !ok || key != "#slow5_version" {
		return formatErr("first line must be #slow5_version, got %q", l)
	}
	v, err := parseVersion(val)
	if err != nil {
		return err
	}
	if err := checkVersion(v); err != nil {
		return err
	}
	rd.h.Version = v
	if l, err = readLine(); err != nil {
		return err
	}
	key, val, ok = strings.Cut(l, "\t")
	if !ok || key != "#num_read_groups" {
		return formatErr("second line must be #num_read_groups, got %q", l)
	}
	nrg, err := strconv.ParseUint(val, 10, 32)
	if err != nil || nrg == 0 {
		return formatErr("invalid num_read_groups %q", val)
	}
	// The header ends with the column-name line, which starts "#read_id".
	var lines []string
	for {
		if l, err = readLine(); err != nil {
			return err
		}
		lines = append(lines, l)
		if strings.HasPrefix(l, "#read_id") {
			break
		}
	}
	return rd.parseHeaderLines(lines, uint32(nrg))
}

// parseHeaderLines parses the data header, type line and column-name line.
func (rd *Reader) parseHeaderLines(lines []string, nrg uint32) error {
	// Bound the allocation a corrupt num_read_groups could cause; real
	// files have a handful of runs.
	if nrg > maxReadGroups {
		return formatErr("num_read_groups %d exceeds the supported %d", nrg, maxReadGroups)
	}
	rd.h.ReadGroups = make([]ReadGroup, nrg)
	for i := range rd.h.ReadGroups {
		rd.h.ReadGroups[i] = ReadGroup{}
	}
	i := 0
	seen := map[string]bool{}
	for ; i < len(lines) && strings.HasPrefix(lines[i], "@"); i++ {
		fields := strings.Split(lines[i][1:], "\t")
		key := fields[0]
		if key == "" || seen[key] {
			return formatErr("empty or duplicate header attribute %q", key)
		}
		seen[key] = true
		if uint64(len(fields)-1) != uint64(nrg) {
			return formatErr("header attribute %s has %d values for %d read groups", key, len(fields)-1, nrg)
		}
		for g, v := range fields[1:] {
			if v != missingASCII {
				rd.h.ReadGroups[g][key] = v
			}
		}
	}
	if len(lines)-i != 2 {
		return formatErr("header must end with the type line and the column-name line")
	}
	types, names := lines[i], lines[i+1]
	if !hasColumns(types, primaryTypeLine) {
		return formatErr("bad type line %q", truncate(types))
	}
	if !hasColumns(names, primaryNameLine) {
		return formatErr("bad column-name line %q", truncate(names))
	}
	auxTypesStr := splitRest(types, primaryTypeLine)
	auxNames := splitRest(names, primaryNameLine)
	if len(auxTypesStr) != len(auxNames) {
		return formatErr("%d auxiliary types but %d auxiliary names", len(auxTypesStr), len(auxNames))
	}
	for k, ts := range auxTypesStr {
		f, err := parseAuxType(auxNames[k], ts)
		if err != nil {
			return err
		}
		rd.h.Aux = append(rd.h.Aux, f)
	}
	return nil
}

func hasColumns(line, prefix string) bool {
	return line == prefix || strings.HasPrefix(line, prefix+"\t")
}

func splitRest(line, prefix string) []string {
	if line == prefix {
		return nil
	}
	return strings.Split(line[len(prefix)+1:], "\t")
}

func truncate(s string) string {
	if len(s) > 80 {
		return s[:80] + "..."
	}
	return s
}

// parseAuxType parses a header type string such as "uint8_t", "double*" or
// "enum{a,b,c}" / "enum*{a,b,c}" (the array form as slow5lib writes it).
func parseAuxType(name, ts string) (AuxField, error) {
	f := AuxField{Name: name}
	if name == "" {
		return f, formatErr("empty auxiliary field name")
	}
	if body, ok := strings.CutPrefix(ts, "enum"); ok && body != "" && body != "*" {
		f.Type = AuxEnum
		if b, arr := strings.CutPrefix(body, "*"); arr {
			f.Type, body = AuxEnumArray, b
		}
		labels, ok := strings.CutPrefix(body, "{")
		if ok {
			labels, ok = strings.CutSuffix(labels, "}")
		}
		if !ok || labels == "" || strings.ContainsAny(labels, "{}") {
			return f, formatErr("field %s: bad enum type %q", name, ts)
		}
		f.EnumLabels = strings.Split(labels, ",")
		if len(f.EnumLabels) > maxEnumLabels {
			return f, formatErr("field %s: %d enum labels (max %d)", name, len(f.EnumLabels), maxEnumLabels)
		}
		return f, nil
	}
	for t := AuxInt8; int(t) < len(auxTypes); t++ {
		if t != AuxEnum && t != AuxEnumArray && auxTypes[t].name == ts {
			f.Type = t
			return f, nil
		}
	}
	return f, formatErr("field %s: unknown type %q", name, ts)
}

// ---------------------------------------------------------------- records

func (rd *Reader) readBLOW5Record() (*Record, error) {
	var size [8]byte
	n, err := io.ReadFull(rd.br, size[:])
	switch {
	case err == nil:
	case n == len(blow5EOF) && bytes.Equal(size[:n], blow5EOF):
		// slow5lib slow5_get_next_mem: a short read of exactly the five
		// marker bytes at end of file is the end-of-file marker.
		return nil, io.EOF
	case n == 0:
		return nil, formatErr("missing BLOW5 end-of-file marker (file truncated)")
	default:
		return nil, formatErr("truncated BLOW5 record size (%d bytes)", n)
	}
	recSize := binary.LittleEndian.Uint64(size[:])
	if recSize > maxDecompressedRecord {
		return nil, formatErr("record size %d is implausible", recSize)
	}
	rd.buf.Reset()
	if got, err := io.CopyN(&rd.buf, rd.br, int64(recSize)); err != nil {
		return nil, formatErr("truncated BLOW5 record: %d of %d bytes", got, recSize)
	}
	payload := rd.buf.Bytes()
	switch rd.rc {
	case RecordZstd:
		payload, err = rd.zdec.DecodeAll(payload, nil)
		if err != nil {
			return nil, formatErr("zstd: %v", err)
		}
	case RecordZlib:
		zr, err := zlib.NewReader(bytes.NewReader(payload))
		if err != nil {
			return nil, formatErr("zlib: %v", err)
		}
		payload, err = io.ReadAll(io.LimitReader(zr, maxDecompressedRecord))
		if err != nil {
			return nil, formatErr("zlib: %v", err)
		}
	}
	return decodeBLOW5Record(payload, &rd.h, rd.sc)
}

// decodeBLOW5Record parses one uncompressed binary record (spec, Table 7).
// Like slow5lib slow5_rec_parse, the fields must consume the record exactly.
func decodeBLOW5Record(b []byte, h *Header, sc SignalCompression) (*Record, error) {
	le := binary.LittleEndian
	p := 0
	take := func(n uint64) ([]byte, error) {
		if n > uint64(len(b)-p) {
			return nil, formatErr("record truncated at byte %d (need %d more)", p, n)
		}
		s := b[p : p+int(n)]
		p += int(n)
		return s, nil
	}
	s, err := take(2)
	if err != nil {
		return nil, err
	}
	id, err := take(uint64(le.Uint16(s)))
	if err != nil {
		return nil, err
	}
	r := &Record{ReadID: string(id)}
	s, err = take(4 + 4*8 + 8)
	if err != nil {
		return nil, err
	}
	r.ReadGroup = le.Uint32(s)
	r.Digitisation = math.Float64frombits(le.Uint64(s[4:]))
	r.Offset = math.Float64frombits(le.Uint64(s[12:]))
	r.Range = math.Float64frombits(le.Uint64(s[20:]))
	r.SamplingRate = math.Float64frombits(le.Uint64(s[28:]))
	nsig := le.Uint64(s[36:])
	if sc == SignalSvbZd {
		blob, err := take(nsig)
		if err != nil {
			return nil, err
		}
		if r.RawSignal, err = svbZdDecode(blob); err != nil {
			return nil, err
		}
	} else {
		if nsig > uint64(len(b)) {
			return nil, formatErr("len_raw_signal %d exceeds record size", nsig)
		}
		raw, err := take(2 * nsig)
		if err != nil {
			return nil, err
		}
		r.RawSignal = make([]int16, nsig)
		for i := range r.RawSignal {
			r.RawSignal[i] = int16(le.Uint16(raw[2*i:]))
		}
	}
	for _, f := range h.Aux {
		info := f.Type.info()
		count := uint64(1)
		if info.arr {
			s, err := take(8)
			if err != nil {
				return nil, err
			}
			count = le.Uint64(s)
		}
		if count > uint64(len(b)) {
			return nil, formatErr("auxiliary field %s: length %d exceeds record size", f.Name, count)
		}
		data, err := take(count * uint64(info.size))
		if err != nil {
			return nil, err
		}
		if info.arr && count == 0 {
			continue // missing
		}
		if f.Type == AuxString {
			// The spec allows a stored NUL terminator counted in the
			// length ("In case it is stored, the preceding size should
			// include the NULL character").
			setAux(r, f.Name, string(bytes.TrimSuffix(data, []byte{0})))
			continue
		}
		elems := make([]any, count)
		for i := range elems {
			elems[i] = decodePrim(info.prim, data[i*info.size:(i+1)*info.size])
		}
		if !info.arr {
			if isSentinel(elems[0]) {
				continue // missing (spec, Table 8)
			}
			setAux(r, f.Name, elems[0])
			continue
		}
		setAux(r, f.Name, makeArray(f.Type, elems))
	}
	if p != len(b) {
		return nil, formatErr("%d unparsed bytes at the end of the record", len(b)-p)
	}
	return r, nil
}

func setAux(r *Record, k string, v any) {
	if r.Aux == nil {
		r.Aux = map[string]any{}
	}
	r.Aux[k] = v
}

func (rd *Reader) readSLOW5Record() (*Record, error) {
	line, err := rd.br.ReadString('\n')
	if err == io.EOF && line == "" {
		return nil, io.EOF
	}
	if err != nil {
		return nil, formatErr("truncated SLOW5 record (no trailing newline)")
	}
	return decodeSLOW5Record(strings.TrimSuffix(line, "\n"), &rd.h)
}

func decodeSLOW5Record(line string, h *Header) (*Record, error) {
	fields := strings.Split(line, "\t")
	if len(fields) != numPrimary+len(h.Aux) {
		return nil, formatErr("record has %d fields, header defines %d", len(fields), numPrimary+len(h.Aux))
	}
	r := &Record{ReadID: fields[0]}
	if r.ReadID == "" {
		return nil, formatErr("empty read_id")
	}
	rg, err := strconv.ParseUint(fields[1], 10, 32)
	if err != nil {
		return nil, formatErr("read %s: bad read_group %q", r.ReadID, fields[1])
	}
	r.ReadGroup = uint32(rg)
	dst := [...]*float64{&r.Digitisation, &r.Offset, &r.Range, &r.SamplingRate}
	for i, d := range dst {
		if *d, err = strconv.ParseFloat(fields[2+i], 64); err != nil {
			return nil, formatErr("read %s: bad %s %q", r.ReadID, []string{"digitisation", "offset", "range", "sampling_rate"}[i], fields[2+i])
		}
	}
	nsig, err := strconv.ParseUint(fields[6], 10, 64)
	if err != nil {
		return nil, formatErr("read %s: bad len_raw_signal %q", r.ReadID, fields[6])
	}
	if r.RawSignal, err = parseSignal(fields[7], nsig); err != nil {
		return nil, formatErr("read %s: %v", r.ReadID, err)
	}
	for k, f := range h.Aux {
		s := fields[numPrimary+k]
		if s == missingASCII {
			continue
		}
		if f.Type == AuxString {
			setAux(r, f.Name, s)
			continue
		}
		info := f.Type.info()
		if !info.arr {
			v, err := parsePrim(f.Type, s)
			if err != nil {
				return nil, formatErr("read %s: field %s: %v", r.ReadID, f.Name, err)
			}
			if !isSentinel(v) {
				setAux(r, f.Name, v)
			}
			continue
		}
		parts := strings.Split(s, ",")
		elems := make([]any, len(parts))
		for i, ps := range parts {
			if elems[i], err = parsePrim(info.prim, ps); err != nil {
				return nil, formatErr("read %s: field %s element %d: %v", r.ReadID, f.Name, i, err)
			}
		}
		setAux(r, f.Name, makeArray(f.Type, elems))
	}
	return r, nil
}

// parseSignal parses the comma-separated raw signal, which must hold
// exactly n values. slow5lib writes an empty field for n == 0; "." is also
// accepted.
func parseSignal(s string, n uint64) ([]int16, error) {
	if n == 0 {
		if s == "" || s == missingASCII {
			return []int16{}, nil
		}
		return nil, formatErr("len_raw_signal is 0 but raw_signal is %q", truncate(s))
	}
	if n > uint64(len(s)+1)/2 {
		return nil, formatErr("raw_signal has fewer than len_raw_signal=%d values", n)
	}
	out := make([]int16, 0, n)
	for len(s) > 0 || len(out) == 0 {
		tok := s
		if i := strings.IndexByte(s, ','); i >= 0 {
			tok, s = s[:i], s[i+1:]
			if s == "" {
				return nil, formatErr("raw_signal ends with a comma")
			}
		} else {
			s = ""
		}
		v, err := strconv.ParseInt(tok, 10, 16)
		if err != nil {
			return nil, formatErr("bad raw_signal value %q", tok)
		}
		out = append(out, int16(v))
	}
	if uint64(len(out)) != n {
		return nil, formatErr("raw_signal has %d values, len_raw_signal is %d", len(out), n)
	}
	return out, nil
}
