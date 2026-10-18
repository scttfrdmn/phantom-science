package slow5

import (
	"errors"
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"
)

// Version is a SLOW5 file format version (spec, MISCELLANEOUS, SLOW5
// versioning). It is independent of the slow5lib/slow5tools version.
type Version struct {
	Major, Minor, Patch uint8
}

// Versions known to this package. Version020 is what new files are written
// as unless Header.Version says otherwise.
var (
	Version010 = Version{0, 1, 0}
	Version020 = Version{0, 2, 0}
)

// String formats v as "major.minor.patch", the form of the #slow5_version
// header line.
func (v Version) String() string {
	return fmt.Sprintf("%d.%d.%d", v.Major, v.Minor, v.Patch)
}

// Less reports whether v is an older version than o.
func (v Version) Less(o Version) bool {
	if v.Major != o.Major {
		return v.Major < o.Major
	}
	if v.Minor != o.Minor {
		return v.Minor < o.Minor
	}
	return v.Patch < o.Patch
}

func (v Version) isZero() bool { return v == Version{} }

// RecordCompression is the BLOW5 record compression method, stored as the
// byte at offset 9 of the BLOW5 header (spec, BLOW5 Header, Table 6). The
// constant values are the on-disk codes.
type RecordCompression uint8

// Record compression methods (spec, Table 6 footnote 1).
const (
	RecordNone RecordCompression = 0
	RecordZlib RecordCompression = 1
	RecordZstd RecordCompression = 2
)

func (c RecordCompression) String() string {
	switch c {
	case RecordNone:
		return "none"
	case RecordZlib:
		return "zlib"
	case RecordZstd:
		return "zstd"
	}
	return fmt.Sprintf("record-compression(%d)", uint8(c))
}

// SignalCompression is the BLOW5 raw-signal compression method, stored as
// the byte at offset 14 of the BLOW5 header from v0.2.0 on (spec, BLOW5
// Header, Table 6). The constant values are the on-disk codes.
type SignalCompression uint8

// Signal compression methods (spec, Table 6 footnote 2).
const (
	SignalNone  SignalCompression = 0
	SignalSvbZd SignalCompression = 1
)

func (c SignalCompression) String() string {
	switch c {
	case SignalNone:
		return "none"
	case SignalSvbZd:
		return "svb-zd"
	}
	return fmt.Sprintf("signal-compression(%d)", uint8(c))
}

// AuxType is the data type of an auxiliary field (spec, SLOW5 Data: the
// supported data types; BLOW5 Data, Table 8).
type AuxType uint8

// Auxiliary field types. The Go type of a value in Record.Aux is given in
// brackets; values must have exactly that type.
const (
	AuxInvalid     AuxType = iota
	AuxInt8                // int8_t [int8]
	AuxInt16               // int16_t [int16]
	AuxInt32               // int32_t [int32]
	AuxInt64               // int64_t [int64]
	AuxUint8               // uint8_t [uint8]
	AuxUint16              // uint16_t [uint16]
	AuxUint32              // uint32_t [uint32]
	AuxUint64              // uint64_t [uint64]
	AuxFloat               // float [float32]
	AuxDouble              // double [float64]
	AuxChar                // char [Char]
	AuxEnum                // enum{...} [Enum]
	AuxInt8Array           // int8_t* [[]int8]
	AuxInt16Array          // int16_t* [[]int16]
	AuxInt32Array          // int32_t* [[]int32]
	AuxInt64Array          // int64_t* [[]int64]
	AuxUint8Array          // uint8_t* [[]uint8]
	AuxUint16Array         // uint16_t* [[]uint16]
	AuxUint32Array         // uint32_t* [[]uint32]
	AuxUint64Array         // uint64_t* [[]uint64]
	AuxFloatArray          // float* [[]float32]
	AuxDoubleArray         // double* [[]float64]
	AuxString              // char* [string]
	AuxEnumArray           // enum* [[]Enum]
)

// Char is the value type of a char auxiliary field (one ASCII byte). It is
// distinct from uint8 so that char and uint8_t fields cannot be confused.
type Char byte

// Enum is the value type of an enum auxiliary field: the 0-based index of a
// label in AuxField.EnumLabels (spec, SLOW5 Data: "8-bit enumeration type").
type Enum uint8

type auxTypeInfo struct {
	name string  // type string in the header, e.g. "uint32_t" or "double*"
	size int     // bytes per (element) value, spec Table 8
	prim AuxType // element type for arrays; the type itself for primitives
	arr  bool
}

var auxTypes = [...]auxTypeInfo{
	AuxInt8:        {"int8_t", 1, AuxInt8, false},
	AuxInt16:       {"int16_t", 2, AuxInt16, false},
	AuxInt32:       {"int32_t", 4, AuxInt32, false},
	AuxInt64:       {"int64_t", 8, AuxInt64, false},
	AuxUint8:       {"uint8_t", 1, AuxUint8, false},
	AuxUint16:      {"uint16_t", 2, AuxUint16, false},
	AuxUint32:      {"uint32_t", 4, AuxUint32, false},
	AuxUint64:      {"uint64_t", 8, AuxUint64, false},
	AuxFloat:       {"float", 4, AuxFloat, false},
	AuxDouble:      {"double", 8, AuxDouble, false},
	AuxChar:        {"char", 1, AuxChar, false},
	AuxEnum:        {"enum", 1, AuxEnum, false},
	AuxInt8Array:   {"int8_t*", 1, AuxInt8, true},
	AuxInt16Array:  {"int16_t*", 2, AuxInt16, true},
	AuxInt32Array:  {"int32_t*", 4, AuxInt32, true},
	AuxInt64Array:  {"int64_t*", 8, AuxInt64, true},
	AuxUint8Array:  {"uint8_t*", 1, AuxUint8, true},
	AuxUint16Array: {"uint16_t*", 2, AuxUint16, true},
	AuxUint32Array: {"uint32_t*", 4, AuxUint32, true},
	AuxUint64Array: {"uint64_t*", 8, AuxUint64, true},
	AuxFloatArray:  {"float*", 4, AuxFloat, true},
	AuxDoubleArray: {"double*", 8, AuxDouble, true},
	AuxString:      {"char*", 1, AuxChar, true},
	AuxEnumArray:   {"enum*", 1, AuxEnum, true},
}

func (t AuxType) valid() bool { return t > AuxInvalid && int(t) < len(auxTypes) }

func (t AuxType) info() auxTypeInfo {
	if !t.valid() {
		return auxTypeInfo{}
	}
	return auxTypes[t]
}

// IsArray reports whether t is a 1D array type (including char*). Arrays are
// stored in BLOW5 with a uint64_t element count prefix (spec, BLOW5 Data).
func (t AuxType) IsArray() bool { return t.info().arr }

// String returns the header type string of t ("uint8_t", "char*", ...). For
// enums the labels are not included; see AuxField.typeString.
func (t AuxType) String() string {
	if !t.valid() {
		return fmt.Sprintf("aux-type(%d)", uint8(t))
	}
	return auxTypes[t].name
}

// AuxField is one auxiliary column of the per-read schema: its name (last
// header line) and type (second-last header line).
type AuxField struct {
	Name string
	Type AuxType
	// EnumLabels are the labels of an AuxEnum or AuxEnumArray field, in
	// value order (label i has value i). Required for those types, ignored
	// otherwise.
	EnumLabels []string
}

// typeString renders the field type for the header type line, including the
// enum labels: enum{a,b,c}, and enum*{a,b,c} for an enum array. The spec
// (SLOW5 Data) shows only the scalar form; the array form is slow5lib's
// (slow5_hdr_types_to_str appends "{labels}" to the type name "enum*").
func (f AuxField) typeString() string {
	s := f.Type.String()
	if f.Type == AuxEnum || f.Type == AuxEnumArray {
		return s + "{" + strings.Join(f.EnumLabels, ",") + "}"
	}
	return s
}

// ONTAuxFields is the auxiliary schema of the standard per-read ONT
// attributes (spec, Auxiliary fields, Table 5), in the order slow5tools f2s
// writes them. ONTAux converts between this schema and a typed struct.
var ONTAuxFields = []AuxField{
	{Name: "channel_number", Type: AuxString},
	{Name: "median_before", Type: AuxDouble},
	{Name: "read_number", Type: AuxInt32},
	{Name: "start_mux", Type: AuxUint8},
	{Name: "start_time", Type: AuxUint64},
}

// ONTAux holds the standard ONT auxiliary fields of one read (spec, Table 5).
type ONTAux struct {
	ChannelNumber string  // channel_number, char*; e.g. "504"
	MedianBefore  float64 // median_before, double; open-pore current, pA
	ReadNumber    int32   // read_number, int32_t; per-channel read counter
	StartMux      uint8   // start_mux, uint8_t; well (MUX) setting, 1-4
	StartTime     uint64  // start_time, uint64_t; in samples since run start
}

// Map returns a as a Record.Aux map for the ONTAuxFields schema.
func (a ONTAux) Map() map[string]any {
	return map[string]any{
		"channel_number": a.ChannelNumber,
		"median_before":  a.MedianBefore,
		"read_number":    a.ReadNumber,
		"start_mux":      a.StartMux,
		"start_time":     a.StartTime,
	}
}

// ONTAuxOf extracts the standard ONT auxiliary fields from r. ok is false if
// any of them is missing or has the wrong type.
func ONTAuxOf(r *Record) (a ONTAux, ok bool) {
	var ok1, ok2, ok3, ok4, ok5 bool
	a.ChannelNumber, ok1 = r.Aux["channel_number"].(string)
	a.MedianBefore, ok2 = r.Aux["median_before"].(float64)
	a.ReadNumber, ok3 = r.Aux["read_number"].(int32)
	a.StartMux, ok4 = r.Aux["start_mux"].(uint8)
	a.StartTime, ok5 = r.Aux["start_time"].(uint64)
	return a, ok1 && ok2 && ok3 && ok4 && ok5
}

// ReadGroup holds the data-header attributes of one read group (one
// sequencing run): run_id, flow_cell_id, asic_id, exp_start_time,
// sample_frequency and so on (spec, Data header, Table 3). All values are
// strings, as in the file.
type ReadGroup map[string]string

// Header is a SLOW5/BLOW5 header.
type Header struct {
	// Version is the file format version. The zero value means Version020.
	// Readers return the version found in the file. A BLOW5 writer bumps a
	// 0.1.0 header to 0.2.0 when zstd or signal compression is requested,
	// as slow5lib does (slow5_press_version_bump).
	Version Version
	// ReadGroups holds the data-header attributes per read group; the read
	// group of a record indexes this slice, so its length is
	// num_read_groups and must be at least 1. A key missing from a group
	// is written as "." (spec, Data header); readers omit such keys.
	ReadGroups []ReadGroup
	// Aux is the auxiliary field schema, in column order.
	Aux []AuxField
}

// Record is one read (spec, Primary fields, Table 4; Auxiliary fields).
type Record struct {
	ReadID       string  // read_id, typically a UUID4
	ReadGroup    uint32  // read_group, index into Header.ReadGroups
	Digitisation float64 // digitisation, ADC quantisation levels
	Offset       float64 // offset, ADC offset
	Range        float64 // range, full-scale range in pA
	SamplingRate float64 // sampling_rate, Hz
	// RawSignal is raw_signal; len_raw_signal is its length. Picoamperes
	// are (RawSignal[i] + Offset) * Range / Digitisation.
	RawSignal []int16
	// Aux holds the auxiliary field values by name; the Go type of each
	// value is fixed by its AuxType. A field of the header schema that is
	// absent from the map is missing ("." in SLOW5 ASCII, the Table 8
	// sentinel in BLOW5). Keys not in the schema are an error.
	Aux map[string]any
}

// PicoAmps converts sample i of the raw signal to picoamperes (spec, Primary
// fields): (raw + offset) * range / digitisation.
func (r *Record) PicoAmps(i int) float64 {
	return (float64(r.RawSignal[i]) + r.Offset) * r.Range / r.Digitisation
}

// primaryTypeLine and primaryNameLine are the fixed starts of the last two
// header lines (spec, Primary fields, Table 4; slow5lib SLOW5_COLS).
const (
	primaryTypeLine = "#char*\tuint32_t\tdouble\tdouble\tdouble\tdouble\tuint64_t\tint16_t*"
	primaryNameLine = "#read_id\tread_group\tdigitisation\toffset\trange\tsampling_rate\tlen_raw_signal\traw_signal"
	numPrimary      = 8
	missingASCII    = "." // spec: missing values are "."
	maxReadIDLen    = math.MaxUint16
	headerSizeAt    = 64 // spec, Table 6: offset of the header size field
	maxEnumLabels   = math.MaxUint8
)

var (
	blow5Magic = []byte{'B', 'L', 'O', 'W', '5', 1} // spec, Table 6, offset 0
	blow5EOF   = []byte{'5', 'W', 'O', 'L', 'B'}    // spec, BLOW5 Footer
)

// ErrFormat is wrapped by every error about malformed input or data that
// SLOW5 cannot represent.
var ErrFormat = errors.New("slow5")

func formatErr(format string, args ...any) error {
	return fmt.Errorf("%w: %s", ErrFormat, fmt.Sprintf(format, args...))
}

// badChars reports whether s contains a character the format reserves:
// tab and newline are delimiters (spec, SLOW5 ASCII) and '\r' is not
// allowed ("'\r' or '\r\n' are not allowed").
func badChars(s string) bool { return strings.ContainsAny(s, "\t\n\r") }

// attrKeys returns the sorted union of the attribute keys of all read
// groups. The data header is sorted in ascending native byte order (spec,
// Data header); sort.Strings compares bytes.
func (h *Header) attrKeys() []string {
	set := map[string]struct{}{}
	for _, g := range h.ReadGroups {
		for k := range g {
			set[k] = struct{}{}
		}
	}
	keys := make([]string, 0, len(set))
	for k := range set {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

// validate checks h for things the format cannot represent.
func (h *Header) validate() error {
	if len(h.ReadGroups) == 0 {
		return formatErr("header has no read groups (num_read_groups must be >= 1)")
	}
	if uint64(len(h.ReadGroups)) > math.MaxUint32 {
		return formatErr("too many read groups: %d", len(h.ReadGroups))
	}
	if Version020.Less(h.version()) {
		return formatErr("version %s is newer than the supported %s", h.version(), Version020)
	}
	for i, g := range h.ReadGroups {
		for k, v := range g {
			if k == "" || badChars(k) {
				return formatErr("read group %d: invalid attribute name %q", i, k)
			}
			if v == "" || v == missingASCII || badChars(v) {
				return formatErr("read group %d: attribute %s has value %q, which SLOW5 cannot store (empty, \".\" or containing tab/newline); omit the key for a missing value", i, k, v)
			}
		}
	}
	seen := map[string]bool{}
	for _, name := range strings.Split(primaryNameLine[1:], "\t") {
		seen[name] = true
	}
	for _, f := range h.Aux {
		if f.Name == "" || badChars(f.Name) || strings.ContainsAny(f.Name, ",") {
			return formatErr("invalid auxiliary field name %q", f.Name)
		}
		if seen[f.Name] {
			return formatErr("duplicate field name %q", f.Name)
		}
		seen[f.Name] = true
		if !f.Type.valid() {
			return formatErr("auxiliary field %s: invalid type %d", f.Name, f.Type)
		}
		if f.Type == AuxEnum || f.Type == AuxEnumArray {
			if len(f.EnumLabels) == 0 || len(f.EnumLabels) > maxEnumLabels {
				return formatErr("enum field %s needs 1..%d labels, has %d", f.Name, maxEnumLabels, len(f.EnumLabels))
			}
			if h.version().Less(Version020) {
				return formatErr("enum field %s needs file version 0.2.0 (enums were added in v0.2.0)", f.Name)
			}
			for _, l := range f.EnumLabels {
				if !validEnumLabel(l) {
					return formatErr("enum field %s: label %q must be alphanumeric/underscore (spec, SLOW5 Data)", f.Name, l)
				}
			}
		}
	}
	return nil
}

// validEnumLabel: "the integral-constant names are restricted to
// alphanumeric characters plus underscores" (spec, SLOW5 Data).
func validEnumLabel(l string) bool {
	if l == "" {
		return false
	}
	for i := 0; i < len(l); i++ {
		c := l[i]
		if !(c == '_' || c >= '0' && c <= '9' || c >= 'a' && c <= 'z' || c >= 'A' && c <= 'Z') {
			return false
		}
	}
	return true
}

func (h *Header) version() Version {
	if h.Version.isZero() {
		return Version020
	}
	return h.Version
}

// validateRecord checks r against h before it is encoded.
func (h *Header) validateRecord(r *Record) error {
	if r == nil {
		return formatErr("nil record")
	}
	if r.ReadID == "" || len(r.ReadID) > maxReadIDLen || badChars(r.ReadID) ||
		r.ReadID[0] == '#' || r.ReadID[0] == '@' {
		// '#' and '@' start header lines (spec, rationale: "Both of those
		// characters are not supposed to be used in read identifiers").
		return formatErr("invalid read_id %q", r.ReadID)
	}
	if uint64(r.ReadGroup) >= uint64(len(h.ReadGroups)) {
		return formatErr("read %s: read_group %d out of range (num_read_groups %d)", r.ReadID, r.ReadGroup, len(h.ReadGroups))
	}
	for k := range r.Aux {
		if h.auxIndex(k) < 0 {
			return formatErr("read %s: auxiliary field %q is not in the header schema", r.ReadID, k)
		}
	}
	return nil
}

func (h *Header) auxIndex(name string) int {
	for i, f := range h.Aux {
		if f.Name == name {
			return i
		}
	}
	return -1
}

// formatDouble renders a float for SLOW5 ASCII. slow5lib prints doubles with
// "%f" and strips trailing zeros (slow5_double_to_str), which loses digits
// beyond the sixth decimal; the spec acknowledges ASCII may be lossy
// (rationale: "Do any of the float/double fields in SLOW5 ASCII become
// lossy"). To stay byte-identical to slow5tools when that form is exact and
// lossless otherwise, the slow5lib form is used when it parses back to the
// same value and the shortest exact decimal (never an exponent) otherwise.
// bits is 32 for float, 64 for double.
func formatDouble(x float64, bits int) string {
	switch {
	case math.IsNaN(x):
		return "nan"
	case math.IsInf(x, 1):
		return "inf"
	case math.IsInf(x, -1):
		return "-inf"
	}
	s := strconv.FormatFloat(x, 'f', 6, 64)
	// slow5lib: strip trailing zeros and the point ("%f" always has
	// one); "-0" becomes "0".
	s = strings.TrimSuffix(strings.TrimRight(s, "0"), ".")
	if s == "-0" {
		s = "0"
	}
	if y, err := strconv.ParseFloat(s, bits); err == nil && y == x {
		return s
	}
	return strconv.FormatFloat(x, 'f', -1, bits)
}
