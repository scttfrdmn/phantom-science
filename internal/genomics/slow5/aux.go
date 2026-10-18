package slow5

import (
	"encoding/binary"
	"math"
	"strconv"
	"strings"
)

// Missing values (spec, BLOW5 Data, Table 8): each primitive type has a
// sentinel that means "missing" ("." in SLOW5 ASCII) and "cannot be used to
// represent the real value". Arrays and strings are missing when their
// length is 0. Writers reject real values equal to a sentinel rather than
// silently turning them into missing values; readers drop sentinel values
// from Record.Aux.

// isSentinel reports whether the primitive value v (already type-checked)
// is the Table 8 missing-value sentinel of its type.
func isSentinel(v any) bool {
	switch x := v.(type) {
	case int8:
		return x == math.MaxInt8
	case int16:
		return x == math.MaxInt16
	case int32:
		return x == math.MaxInt32
	case int64:
		return x == math.MaxInt64
	case uint8:
		return x == math.MaxUint8
	case uint16:
		return x == math.MaxUint16
	case uint32:
		return x == math.MaxUint32
	case uint64:
		return x == math.MaxUint64
	case float32:
		return math.IsNaN(float64(x))
	case float64:
		return math.IsNaN(x)
	case Char:
		return x == 0
	case Enum:
		return x == math.MaxUint8
	}
	return false
}

// primValues returns the elements of an aux value as a []any of the
// primitive Go type, checking that v has the Go type t requires. For
// primitive types it returns a single element.
func primValues(f AuxField, v any) ([]any, error) {
	bad := func() ([]any, error) {
		return nil, formatErr("auxiliary field %s (%s): value has Go type %T", f.Name, f.Type, v)
	}
	switch f.Type {
	case AuxInt8:
		if x, ok := v.(int8); ok {
			return []any{x}, nil
		}
	case AuxInt16:
		if x, ok := v.(int16); ok {
			return []any{x}, nil
		}
	case AuxInt32:
		if x, ok := v.(int32); ok {
			return []any{x}, nil
		}
	case AuxInt64:
		if x, ok := v.(int64); ok {
			return []any{x}, nil
		}
	case AuxUint8:
		if x, ok := v.(uint8); ok {
			return []any{x}, nil
		}
	case AuxUint16:
		if x, ok := v.(uint16); ok {
			return []any{x}, nil
		}
	case AuxUint32:
		if x, ok := v.(uint32); ok {
			return []any{x}, nil
		}
	case AuxUint64:
		if x, ok := v.(uint64); ok {
			return []any{x}, nil
		}
	case AuxFloat:
		if x, ok := v.(float32); ok {
			return []any{x}, nil
		}
	case AuxDouble:
		if x, ok := v.(float64); ok {
			return []any{x}, nil
		}
	case AuxChar:
		if x, ok := v.(Char); ok {
			return []any{x}, nil
		}
	case AuxEnum:
		if x, ok := v.(Enum); ok {
			return []any{x}, nil
		}
	case AuxInt8Array:
		return sliceAny[int8](v, bad)
	case AuxInt16Array:
		return sliceAny[int16](v, bad)
	case AuxInt32Array:
		return sliceAny[int32](v, bad)
	case AuxInt64Array:
		return sliceAny[int64](v, bad)
	case AuxUint8Array:
		return sliceAny[uint8](v, bad)
	case AuxUint16Array:
		return sliceAny[uint16](v, bad)
	case AuxUint32Array:
		return sliceAny[uint32](v, bad)
	case AuxUint64Array:
		return sliceAny[uint64](v, bad)
	case AuxFloatArray:
		return sliceAny[float32](v, bad)
	case AuxDoubleArray:
		return sliceAny[float64](v, bad)
	case AuxEnumArray:
		return sliceAny[Enum](v, bad)
	}
	return bad()
}

func sliceAny[T any](v any, bad func() ([]any, error)) ([]any, error) {
	s, ok := v.([]T)
	if !ok {
		return bad()
	}
	out := make([]any, len(s))
	for i, x := range s {
		out[i] = x
	}
	return out, nil
}

// checkAuxValue validates a present aux value: Go type, no sentinels, no
// empty arrays, enum values within the labels, and strings that SLOW5 ASCII
// can hold. For char* fields it returns no elements; callers use the string.
func checkAuxValue(f AuxField, v any) ([]any, error) {
	if f.Type == AuxString {
		s, ok := v.(string)
		if !ok {
			return nil, formatErr("auxiliary field %s (char*): value has Go type %T", f.Name, v)
		}
		if s == "" {
			return nil, formatErr("auxiliary field %s: empty char* is the missing value; omit the field instead", f.Name)
		}
		// NUL is allowed by BLOW5 only as a terminator; "." is ASCII's
		// missing value.
		if s == missingASCII || badChars(s) || strings.IndexByte(s, 0) >= 0 {
			return nil, formatErr("auxiliary field %s: string %q cannot be stored (\".\", tab, newline or NUL)", f.Name, s)
		}
		return nil, nil
	}
	elems, err := primValues(f, v)
	if err != nil {
		return nil, err
	}
	if f.Type.IsArray() && len(elems) == 0 {
		return nil, formatErr("auxiliary field %s: empty %s is the missing value; omit the field instead", f.Name, f.Type)
	}
	for _, e := range elems {
		if isSentinel(e) {
			return nil, formatErr("auxiliary field %s: value %v is the SLOW5 missing-value sentinel for %s (spec Table 8); omit the field instead", f.Name, e, f.Type.info().prim)
		}
		if c, ok := e.(Char); ok && (c == '\t' || c == '\n' || c == '\r' || c == '.' || c >= 0x80) {
			return nil, formatErr("auxiliary field %s: char %q cannot be stored", f.Name, byte(c))
		}
		if en, ok := e.(Enum); ok && int(en) >= len(f.EnumLabels) {
			return nil, formatErr("auxiliary field %s: enum value %d has no label (%d labels)", f.Name, en, len(f.EnumLabels))
		}
	}
	return elems, nil
}

// appendPrim appends the little-endian encoding of one primitive value
// (spec, BLOW5: "all multi-byte numbers are stored in little-endian").
func appendPrim(b []byte, v any) []byte {
	switch x := v.(type) {
	case int8:
		return append(b, byte(x))
	case uint8:
		return append(b, x)
	case Char:
		return append(b, byte(x))
	case Enum:
		return append(b, byte(x))
	case int16:
		return binary.LittleEndian.AppendUint16(b, uint16(x))
	case uint16:
		return binary.LittleEndian.AppendUint16(b, x)
	case int32:
		return binary.LittleEndian.AppendUint32(b, uint32(x))
	case uint32:
		return binary.LittleEndian.AppendUint32(b, x)
	case float32:
		return binary.LittleEndian.AppendUint32(b, math.Float32bits(x))
	case int64:
		return binary.LittleEndian.AppendUint64(b, uint64(x))
	case uint64:
		return binary.LittleEndian.AppendUint64(b, x)
	case float64:
		return binary.LittleEndian.AppendUint64(b, math.Float64bits(x))
	}
	panic("slow5: appendPrim: unexpected type")
}

// appendSentinel appends the Table 8 missing value of primitive type t.
// slow5lib writes NaN as nan(""), whose bits are the quiet NaN
// 0x7ff8000000000000 (float: 0x7fc00000) on the platforms it supports.
func appendSentinel(b []byte, t AuxType) []byte {
	switch t {
	case AuxInt8:
		return append(b, math.MaxInt8)
	case AuxUint8, AuxEnum:
		return append(b, math.MaxUint8)
	case AuxChar:
		return append(b, 0)
	case AuxInt16:
		return binary.LittleEndian.AppendUint16(b, math.MaxInt16)
	case AuxUint16:
		return binary.LittleEndian.AppendUint16(b, math.MaxUint16)
	case AuxInt32:
		return binary.LittleEndian.AppendUint32(b, math.MaxInt32)
	case AuxUint32:
		return binary.LittleEndian.AppendUint32(b, math.MaxUint32)
	case AuxFloat:
		return binary.LittleEndian.AppendUint32(b, 0x7fc00000)
	case AuxInt64:
		return binary.LittleEndian.AppendUint64(b, math.MaxInt64)
	case AuxUint64:
		return binary.LittleEndian.AppendUint64(b, math.MaxUint64)
	case AuxDouble:
		return binary.LittleEndian.AppendUint64(b, 0x7ff8000000000000)
	}
	panic("slow5: appendSentinel: unexpected type")
}

// decodePrim decodes one primitive value of type t from b (len(b) ==
// size of t).
func decodePrim(t AuxType, b []byte) any {
	le := binary.LittleEndian
	switch t {
	case AuxInt8:
		return int8(b[0])
	case AuxUint8:
		return b[0]
	case AuxChar:
		return Char(b[0])
	case AuxEnum:
		return Enum(b[0])
	case AuxInt16:
		return int16(le.Uint16(b))
	case AuxUint16:
		return le.Uint16(b)
	case AuxInt32:
		return int32(le.Uint32(b))
	case AuxUint32:
		return le.Uint32(b)
	case AuxFloat:
		return math.Float32frombits(le.Uint32(b))
	case AuxInt64:
		return int64(le.Uint64(b))
	case AuxUint64:
		return le.Uint64(b)
	case AuxDouble:
		return math.Float64frombits(le.Uint64(b))
	}
	panic("slow5: decodePrim: unexpected type")
}

// makeArray builds the typed Go slice for array type t from decoded
// primitive elements.
func makeArray(t AuxType, elems []any) any {
	switch t {
	case AuxInt8Array:
		return typedSlice[int8](elems)
	case AuxInt16Array:
		return typedSlice[int16](elems)
	case AuxInt32Array:
		return typedSlice[int32](elems)
	case AuxInt64Array:
		return typedSlice[int64](elems)
	case AuxUint8Array:
		return typedSlice[uint8](elems)
	case AuxUint16Array:
		return typedSlice[uint16](elems)
	case AuxUint32Array:
		return typedSlice[uint32](elems)
	case AuxUint64Array:
		return typedSlice[uint64](elems)
	case AuxFloatArray:
		return typedSlice[float32](elems)
	case AuxDoubleArray:
		return typedSlice[float64](elems)
	case AuxEnumArray:
		return typedSlice[Enum](elems)
	}
	panic("slow5: makeArray: unexpected type")
}

func typedSlice[T any](elems []any) []T {
	out := make([]T, len(elems))
	for i, e := range elems {
		out[i] = e.(T)
	}
	return out
}

// formatPrim renders one primitive value for SLOW5 ASCII (slow5lib
// slow5_data_to_str: integers in decimal, enums as their numeric value,
// char as the character, float/double as in formatDouble).
func formatPrim(v any) string {
	switch x := v.(type) {
	case int8:
		return strconv.FormatInt(int64(x), 10)
	case int16:
		return strconv.FormatInt(int64(x), 10)
	case int32:
		return strconv.FormatInt(int64(x), 10)
	case int64:
		return strconv.FormatInt(x, 10)
	case uint8:
		return strconv.FormatUint(uint64(x), 10)
	case uint16:
		return strconv.FormatUint(uint64(x), 10)
	case uint32:
		return strconv.FormatUint(uint64(x), 10)
	case uint64:
		return strconv.FormatUint(x, 10)
	case Enum:
		return strconv.FormatUint(uint64(x), 10)
	case Char:
		return string([]byte{byte(x)})
	case float32:
		return formatDouble(float64(x), 32)
	case float64:
		return formatDouble(x, 64)
	}
	panic("slow5: formatPrim: unexpected type")
}

// parsePrim parses one primitive value of type t from SLOW5 ASCII.
func parsePrim(t AuxType, s string) (any, error) {
	var (
		v   any
		err error
	)
	switch t {
	case AuxInt8, AuxInt16, AuxInt32, AuxInt64:
		bits := auxTypes[t].size * 8
		var n int64
		n, err = strconv.ParseInt(s, 10, bits)
		switch t {
		case AuxInt8:
			v = int8(n)
		case AuxInt16:
			v = int16(n)
		case AuxInt32:
			v = int32(n)
		default:
			v = n
		}
	case AuxUint8, AuxUint16, AuxUint32, AuxUint64, AuxEnum:
		bits := auxTypes[t].size * 8
		var n uint64
		n, err = strconv.ParseUint(s, 10, bits)
		switch t {
		case AuxUint8:
			v = uint8(n)
		case AuxEnum:
			v = Enum(n)
		case AuxUint16:
			v = uint16(n)
		case AuxUint32:
			v = uint32(n)
		default:
			v = n
		}
	case AuxFloat:
		var f float64
		f, err = strconv.ParseFloat(s, 32)
		v = float32(f)
	case AuxDouble:
		v, err = strconv.ParseFloat(s, 64)
	case AuxChar:
		if len(s) != 1 {
			return nil, formatErr("char value %q is not one character", s)
		}
		v = Char(s[0])
	default:
		return nil, formatErr("parsePrim: unexpected type %s", t)
	}
	if err != nil {
		return nil, formatErr("cannot parse %s from %q", t, s)
	}
	return v, nil
}
