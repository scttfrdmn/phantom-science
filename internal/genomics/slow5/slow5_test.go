//go:build !release

package slow5

import (
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

// ------------------------------------------------------------ helpers

func testRand(seed uint64) *rand.Rand { return rand.New(rand.NewPCG(seed, 0x510775)) }

// allAuxFields covers every auxiliary type of the specification.
var allAuxFields = []AuxField{
	{Name: "channel_number", Type: AuxString},
	{Name: "median_before", Type: AuxDouble},
	{Name: "read_number", Type: AuxInt32},
	{Name: "start_mux", Type: AuxUint8},
	{Name: "start_time", Type: AuxUint64},
	{Name: "i8", Type: AuxInt8},
	{Name: "i16", Type: AuxInt16},
	{Name: "i64", Type: AuxInt64},
	{Name: "u16", Type: AuxUint16},
	{Name: "u32", Type: AuxUint32},
	{Name: "f32", Type: AuxFloat},
	{Name: "ch", Type: AuxChar},
	{Name: "end_reason", Type: AuxEnum, EnumLabels: []string{"unknown", "partial", "mux_change", "unblock_mux_change", "signal_positive", "signal_negative"}},
	{Name: "a_i8", Type: AuxInt8Array},
	{Name: "a_i16", Type: AuxInt16Array},
	{Name: "a_i32", Type: AuxInt32Array},
	{Name: "a_i64", Type: AuxInt64Array},
	{Name: "a_u8", Type: AuxUint8Array},
	{Name: "a_u16", Type: AuxUint16Array},
	{Name: "a_u32", Type: AuxUint32Array},
	{Name: "a_u64", Type: AuxUint64Array},
	{Name: "a_f32", Type: AuxFloatArray},
	{Name: "a_f64", Type: AuxDoubleArray},
	{Name: "a_enum", Type: AuxEnumArray, EnumLabels: []string{"a", "b", "c"}},
}

func multiGroupHeader() Header {
	return Header{
		ReadGroups: []ReadGroup{
			{"run_id": "d6e473a6d513ec6bfc150c60fd4556d72f0e6d18", "flow_cell_id": "PAW12345", "asic_id": "0004A30B00232BEC", "exp_start_time": "2024-03-01T10:00:00Z", "sample_frequency": "5000", "sample_id": "NA12878"},
			{"run_id": "1111111111111111111111111111111111111111", "flow_cell_id": "PAW12346", "asic_id": "1004A30B00232BEC", "exp_start_time": "2024-03-02T10:00:00Z", "sample_frequency": "5000"},
			{"run_id": "2222222222222222222222222222222222222222", "flow_cell_id": "PAW12347", "sample_frequency": "4000", "device_type": "promethion"},
		},
		Aux: allAuxFields,
	}
}

// smoothSignal is a random walk around a pore level, like nanopore current.
func smoothSignal(r *rand.Rand, n int) []int16 {
	s := make([]int16, n)
	level := 500.0
	for i := range s {
		if r.IntN(10) == 0 {
			level = 400 + 300*r.Float64()
		}
		v := level + 8*r.NormFloat64()
		s[i] = int16(math.Round(v))
	}
	return s
}

// fullRecord fills every aux field with non-sentinel values.
func fullRecord(r *rand.Rand, i int, group uint32, nsig int) *Record {
	return &Record{
		ReadID:       fmt.Sprintf("%08x-c43d-492a-b6a1-%012x", r.Uint32(), i),
		ReadGroup:    group,
		Digitisation: 2048,
		Offset:       float64(-250 + r.IntN(50)),
		Range:        748.580139160156,
		SamplingRate: 5000,
		RawSignal:    smoothSignal(r, nsig),
		Aux: map[string]any{
			"channel_number": fmt.Sprint(1 + r.IntN(2675)),
			"median_before":  200 + 50*r.Float64(),
			"read_number":    int32(r.IntN(100000)),
			"start_mux":      uint8(1 + r.IntN(4)),
			"start_time":     r.Uint64() >> 20,
			"i8":             int8(-128),
			"i16":            int16(-32768),
			"i64":            int64(math.MinInt64),
			"u16":            uint16(65534),
			"u32":            uint32(4294967294),
			"f32":            float32(0.1),
			"ch":             Char('Z'),
			"end_reason":     Enum(r.IntN(6)),
			"a_i8":           []int8{-128, 0, 126},
			"a_i16":          []int16{-32768, 1, 32766},
			"a_i32":          []int32{math.MinInt32, 7},
			"a_i64":          []int64{math.MinInt64, 0, 42},
			"a_u8":           []uint8{0, 254},
			"a_u16":          []uint16{0, 65534},
			"a_u32":          []uint32{0, 4294967294},
			"a_u64":          []uint64{0, math.MaxUint64 - 1},
			"a_f32":          []float32{1.5, -0.25, 3.4028235e38},
			"a_f64":          []float64{math.Pi, -1e-300, 1e300},
			"a_enum":         []Enum{0, 2, 1},
		},
	}
}

func sameFloat(a, b float64) bool {
	return math.Float64bits(a) == math.Float64bits(b) || (math.IsNaN(a) && math.IsNaN(b))
}

// recordDiff returns "" if a and b are identical (floats bit-for-bit).
func recordDiff(a, b *Record) string {
	switch {
	case a.ReadID != b.ReadID:
		return fmt.Sprintf("read_id %q vs %q", a.ReadID, b.ReadID)
	case a.ReadGroup != b.ReadGroup:
		return fmt.Sprintf("read_group %d vs %d", a.ReadGroup, b.ReadGroup)
	case !sameFloat(a.Digitisation, b.Digitisation), !sameFloat(a.Offset, b.Offset),
		!sameFloat(a.Range, b.Range), !sameFloat(a.SamplingRate, b.SamplingRate):
		return fmt.Sprintf("calibration %v/%v/%v/%v vs %v/%v/%v/%v", a.Digitisation, a.Offset, a.Range, a.SamplingRate, b.Digitisation, b.Offset, b.Range, b.SamplingRate)
	case len(a.RawSignal) != len(b.RawSignal):
		return fmt.Sprintf("signal length %d vs %d", len(a.RawSignal), len(b.RawSignal))
	}
	for i := range a.RawSignal {
		if a.RawSignal[i] != b.RawSignal[i] {
			return fmt.Sprintf("signal[%d] %d vs %d", i, a.RawSignal[i], b.RawSignal[i])
		}
	}
	if len(a.Aux) != len(b.Aux) {
		return fmt.Sprintf("aux %v vs %v", a.Aux, b.Aux)
	}
	for k, v := range a.Aux {
		if !reflect.DeepEqual(v, b.Aux[k]) {
			return fmt.Sprintf("aux %s: %#v vs %#v", k, v, b.Aux[k])
		}
	}
	return ""
}

func headerDiff(a, b Header) string {
	if a.Version != b.Version && !(a.Version.isZero() || b.Version.isZero()) {
		return fmt.Sprintf("version %s vs %s", a.Version, b.Version)
	}
	if len(a.ReadGroups) != len(b.ReadGroups) {
		return fmt.Sprintf("%d vs %d read groups", len(a.ReadGroups), len(b.ReadGroups))
	}
	for i := range a.ReadGroups {
		if !reflect.DeepEqual(a.ReadGroups[i], b.ReadGroups[i]) {
			return fmt.Sprintf("read group %d: %v vs %v", i, a.ReadGroups[i], b.ReadGroups[i])
		}
	}
	if len(a.Aux) != len(b.Aux) {
		return fmt.Sprintf("aux schema %v vs %v", a.Aux, b.Aux)
	}
	for i := range a.Aux {
		if a.Aux[i].Name != b.Aux[i].Name || a.Aux[i].Type != b.Aux[i].Type ||
			strings.Join(a.Aux[i].EnumLabels, ",") != strings.Join(b.Aux[i].EnumLabels, ",") {
			return fmt.Sprintf("aux field %d: %+v vs %+v", i, a.Aux[i], b.Aux[i])
		}
	}
	return ""
}

type codec struct {
	name  string
	ascii bool
	opts  BLOW5Options
}

var allCodecs = []codec{
	{"slow5", true, BLOW5Options{}},
	{"blow5-none-none", false, BLOW5Options{RecordNone, SignalNone}},
	{"blow5-none-svbzd", false, BLOW5Options{RecordNone, SignalSvbZd}},
	{"blow5-zlib-none", false, BLOW5Options{RecordZlib, SignalNone}},
	{"blow5-zlib-svbzd", false, BLOW5Options{RecordZlib, SignalSvbZd}},
	{"blow5-zstd-none", false, BLOW5Options{RecordZstd, SignalNone}},
	{"blow5-zstd-svbzd", false, BLOW5Options{RecordZstd, SignalSvbZd}},
}

func writeAll(t testing.TB, c codec, h Header, recs []*Record) []byte {
	t.Helper()
	var buf bytes.Buffer
	var w Writer
	var err error
	if c.ascii {
		w, err = NewSLOW5Writer(&buf, h)
	} else {
		w, err = NewBLOW5Writer(&buf, h, c.opts)
	}
	if err != nil {
		t.Fatalf("%s: new writer: %v", c.name, err)
	}
	for _, r := range recs {
		if err := w.Write(r); err != nil {
			t.Fatalf("%s: write %s: %v", c.name, r.ReadID, err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("%s: close: %v", c.name, err)
	}
	return buf.Bytes()
}

func readAllBytes(t testing.TB, b []byte) (Header, []*Record) {
	t.Helper()
	h, recs, err := ReadAll(bytes.NewReader(b))
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	return h, recs
}

// ------------------------------------------------------------ round trips

func TestRoundTripAllCodecs(t *testing.T) {
	r := testRand(1)
	h := multiGroupHeader()
	var recs []*Record
	// Short, medium and empty signals across all read groups, and records
	// with every aux field missing or only some present.
	for i, n := range []int{1, 2, 3, 4, 5, 7, 8, 100, 4000, 0, 33333} {
		recs = append(recs, fullRecord(r, i, uint32(i%3), n))
	}
	sparse := fullRecord(r, 99, 1, 50)
	sparse.Aux = map[string]any{"channel_number": "7", "a_u8": []uint8{1}}
	recs = append(recs, sparse)
	none := fullRecord(r, 100, 2, 10)
	none.Aux = nil
	recs = append(recs, none)

	for _, c := range allCodecs {
		t.Run(c.name, func(t *testing.T) {
			data := writeAll(t, c, h, recs)
			gotH, got, err := ReadAll(bytes.NewReader(data))
			if err != nil {
				t.Fatal(err)
			}
			if d := headerDiff(h, gotH); d != "" {
				t.Fatalf("header: %s", d)
			}
			if gotH.Version != Version020 {
				t.Errorf("version %s, want 0.2.0", gotH.Version)
			}
			if len(got) != len(recs) {
				t.Fatalf("read %d records, wrote %d", len(got), len(recs))
			}
			for i := range recs {
				if d := recordDiff(recs[i], got[i]); d != "" {
					t.Errorf("record %d: %s", i, d)
				}
			}
		})
	}
}

func TestASCIIAndBLOW5DecodeIdentically(t *testing.T) {
	r := testRand(2)
	h := multiGroupHeader()
	var recs []*Record
	for i := 0; i < 20; i++ {
		recs = append(recs, fullRecord(r, i, uint32(i%3), 1+r.IntN(5000)))
	}
	_, ascii := readAllBytes(t, writeAll(t, allCodecs[0], h, recs))
	for _, c := range allCodecs[1:] {
		_, bin := readAllBytes(t, writeAll(t, c, h, recs))
		for i := range ascii {
			if d := recordDiff(ascii[i], bin[i]); d != "" {
				t.Errorf("%s record %d: ASCII vs BLOW5: %s", c.name, i, d)
			}
		}
	}
}

func TestLongSignalRoundTrip(t *testing.T) {
	if testing.Short() {
		t.Skip("long signal")
	}
	r := testRand(3)
	h := Header{ReadGroups: []ReadGroup{{"run_id": "r"}}, Aux: ONTAuxFields}
	rec := &Record{ReadID: "long", Digitisation: 8192, Offset: 6, Range: 1467.6, SamplingRate: 4000,
		RawSignal: smoothSignal(r, 1_200_000),
		Aux:       ONTAux{ChannelNumber: "1", MedianBefore: 230.5, ReadNumber: 9, StartMux: 2, StartTime: 1 << 40}.Map()}
	// Extremes at the start and in the middle exercise 3-byte svb codes.
	rec.RawSignal[0], rec.RawSignal[1] = math.MaxInt16, math.MinInt16
	rec.RawSignal[600000], rec.RawSignal[600001] = math.MinInt16, math.MaxInt16
	for _, c := range allCodecs {
		_, got := readAllBytes(t, writeAll(t, c, h, []*Record{rec}))
		if d := recordDiff(rec, got[0]); d != "" {
			t.Errorf("%s: %s", c.name, d)
		}
	}
}

func TestONTAux(t *testing.T) {
	a := ONTAux{ChannelNumber: "504", MedianBefore: 238.78225708007812, ReadNumber: 17981, StartMux: 4, StartTime: 335845487}
	h := Header{ReadGroups: []ReadGroup{{"run_id": "x"}}, Aux: ONTAuxFields}
	rec := &Record{ReadID: "r", Digitisation: 1, Range: 1, SamplingRate: 1, Aux: a.Map()}
	_, got := readAllBytes(t, writeAll(t, allCodecs[6], h, []*Record{rec}))
	b, ok := ONTAuxOf(got[0])
	if !ok || b != a {
		t.Fatalf("ONTAuxOf = %+v, %v; want %+v", b, ok, a)
	}
	if _, ok := ONTAuxOf(&Record{}); ok {
		t.Error("ONTAuxOf of an empty record reported ok")
	}
}

func TestPicoAmps(t *testing.T) {
	// spec, Primary fields: (raw + offset) * range / digitisation.
	r := &Record{Digitisation: 8192, Offset: 23, Range: 1467.61, RawSignal: []int16{1039}}
	want := (1039 + 23) * 1467.61 / 8192
	if got := r.PicoAmps(0); math.Abs(got-want) > 1e-12 {
		t.Fatalf("PicoAmps = %v, want %v", got, want)
	}
}

// ------------------------------------------------------------ svb-zd

func TestSvbZdRoundTripProperties(t *testing.T) {
	r := testRand(4)
	cases := map[string][]int16{
		"empty":    {},
		"one":      {-5},
		"extremes": {math.MaxInt16, math.MinInt16, math.MaxInt16, math.MinInt16, 0, -1, 1, math.MinInt16},
		"constant": make([]int16, 1001),
	}
	ramp := make([]int16, 70000)
	for i := range ramp {
		ramp[i] = int16(i - 35000)
	}
	cases["ramp"] = ramp
	down := make([]int16, 65536)
	for i := range down {
		down[i] = int16(math.MaxInt16 - i)
	}
	cases["down-ramp-full-range"] = down
	for k := 0; k < 200; k++ {
		n := r.IntN(37)
		s := make([]int16, n)
		for i := range s {
			s[i] = int16(r.Uint32())
		}
		cases[fmt.Sprintf("random-%d", k)] = s
	}
	for name, sig := range cases {
		enc := svbZdEncode(nil, sig)
		if got := binary.LittleEndian.Uint32(enc); int(got) != len(sig) {
			t.Errorf("%s: length prefix %d, want %d", name, got, len(sig))
		}
		dec, err := svbZdDecode(enc)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if !reflect.DeepEqual(dec, sig) && !(len(dec) == 0 && len(sig) == 0) {
			t.Fatalf("%s: round trip mismatch", name)
		}
	}
}

func TestSvbZdKnownBytes(t *testing.T) {
	// Hand-computed from slow5lib's algorithm: samples 1, 3, 2, 300, -32768.
	// deltas 1, 2, -1, 298, -33068; zig-zag 2, 4, 1, 596, 66135.
	// 596 = 0x0254 (2 bytes, code 1); 66135 = 0x010257 (3 bytes, code 2).
	// key0 = 0|0<<2|0<<4|1<<6 = 0x40; key1 = 2.
	got := svbZdEncode(nil, []int16{1, 3, 2, 300, -32768})
	want := []byte{5, 0, 0, 0, 0x40, 0x02, 2, 4, 1, 0x54, 0x02, 0x57, 0x02, 0x01}
	if !bytes.Equal(got, want) {
		t.Fatalf("svb-zd bytes\n got % x\nwant % x", got, want)
	}
}

func TestSvbZdSmallerThanRawForSmoothSignal(t *testing.T) {
	sig := smoothSignal(testRand(5), 100000)
	enc := svbZdEncode(nil, sig)
	if raw := 2 * len(sig); len(enc) >= raw*3/4 {
		t.Fatalf("svb-zd %d bytes is not well below raw %d bytes for a smooth signal", len(enc), raw)
	}
}

func TestSvbZdDecodeRejectsCorruption(t *testing.T) {
	enc := svbZdEncode(nil, []int16{10, 20, 30, 40, 50})
	bad := map[string][]byte{
		"short":         enc[:3],
		"truncated":     enc[:len(enc)-1],
		"trailing":      append(append([]byte{}, enc...), 0),
		"huge-count":    {0xff, 0xff, 0xff, 0x7f, 0},
		"out-of-int16":  {2, 0, 0, 0, 0x02, 0xfe, 0xff, 0x01, 0xfe, 0xff, 0x01}, // two deltas of +65535
		"zero-but-data": {0, 0, 0, 0, 1},
	}
	for name, b := range bad {
		if _, err := svbZdDecode(b); err == nil {
			t.Errorf("%s: corrupt svb-zd accepted", name)
		}
	}
}

// ------------------------------------------------------------ byte layout

func TestBLOW5HeaderBytes(t *testing.T) {
	h := Header{
		ReadGroups: []ReadGroup{{"run_id": "abc", "flow_cell_id": "FAX1"}, {"run_id": "def"}},
		Aux:        []AuxField{{Name: "start_mux", Type: AuxUint8}, {Name: "end_reason", Type: AuxEnum, EnumLabels: []string{"unknown", "partial"}}},
	}
	rec := &Record{ReadID: "r0", ReadGroup: 1, Digitisation: 8192, Offset: 6, Range: 1467.6, SamplingRate: 4000, RawSignal: []int16{1, -2, 3},
		Aux: map[string]any{"start_mux": uint8(3), "end_reason": Enum(1)}}
	data := writeAll(t, codec{"x", false, BLOW5Options{RecordNone, SignalSvbZd}}, h, []*Record{rec})

	// spec, Table 6.
	if !bytes.Equal(data[:6], []byte("BLOW5\x01")) {
		t.Fatalf("magic % x", data[:6])
	}
	if data[6] != 0 || data[7] != 2 || data[8] != 0 {
		t.Errorf("version bytes %d.%d.%d, want 0.2.0", data[6], data[7], data[8])
	}
	if data[9] != 0 {
		t.Errorf("record compression byte %d, want 0", data[9])
	}
	if n := binary.LittleEndian.Uint32(data[10:]); n != 2 {
		t.Errorf("num_read_groups %d, want 2", n)
	}
	if data[14] != 1 {
		t.Errorf("signal compression byte %d, want 1 (svb-zd)", data[14])
	}
	for i := 15; i < 64; i++ {
		if data[i] != 0 {
			t.Fatalf("reserved byte %d = %d, must be zero", i, data[i])
		}
	}
	size := binary.LittleEndian.Uint32(data[64:])
	text := string(data[68 : 68+size])
	wantText := "@flow_cell_id\tFAX1\t.\n" +
		"@run_id\tabc\tdef\n" +
		"#char*\tuint32_t\tdouble\tdouble\tdouble\tdouble\tuint64_t\tint16_t*\tuint8_t\tenum{unknown,partial}\n" +
		"#read_id\tread_group\tdigitisation\toffset\trange\tsampling_rate\tlen_raw_signal\traw_signal\tstart_mux\tend_reason\n"
	if text != wantText {
		t.Fatalf("header text\n%q\nwant\n%q", text, wantText)
	}
	// One record: uint64 size, then exactly that many bytes, then 5WOLB.
	off := 68 + int(size)
	recSize := binary.LittleEndian.Uint64(data[off:])
	off += 8
	if off+int(recSize)+5 != len(data) {
		t.Fatalf("record size %d does not reach the EOF marker (file %d bytes)", recSize, len(data))
	}
	if string(data[len(data)-5:]) != "5WOLB" {
		t.Fatalf("EOF marker %q", data[len(data)-5:])
	}
	rb := data[off : off+int(recSize)]
	// spec, Table 7: uint16 id length, id, uint32 read group, 4 doubles,
	// uint64 len_raw_signal (here the svb-zd byte count), signal, aux.
	if binary.LittleEndian.Uint16(rb) != 2 || string(rb[2:4]) != "r0" || binary.LittleEndian.Uint32(rb[4:]) != 1 {
		t.Fatalf("record head % x", rb[:8])
	}
	if math.Float64frombits(binary.LittleEndian.Uint64(rb[8:])) != 8192 || math.Float64frombits(binary.LittleEndian.Uint64(rb[32:])) != 4000 {
		t.Fatal("digitisation/sampling_rate not little-endian doubles at the spec offsets")
	}
	blobLen := binary.LittleEndian.Uint64(rb[40:])
	blob := svbZdEncode(nil, []int16{1, -2, 3})
	if blobLen != uint64(len(blob)) || !bytes.Equal(rb[48:48+blobLen], blob) {
		t.Fatalf("len_raw_signal %d, signal % x; want %d, % x", blobLen, rb[48:48+blobLen], len(blob), blob)
	}
	if tail := rb[48+blobLen:]; !bytes.Equal(tail, []byte{3, 1}) {
		t.Fatalf("aux bytes % x, want 03 01", tail)
	}
}

func TestBLOW5MissingAuxSentinels(t *testing.T) {
	h := Header{ReadGroups: []ReadGroup{{"run_id": "x"}}, Aux: []AuxField{
		{Name: "s", Type: AuxString}, {Name: "d", Type: AuxDouble}, {Name: "i32", Type: AuxInt32},
		{Name: "u8", Type: AuxUint8}, {Name: "u64", Type: AuxUint64}, {Name: "arr", Type: AuxInt16Array},
	}}
	rec := &Record{ReadID: "r", Digitisation: 1, Range: 1, SamplingRate: 1}
	data := writeAll(t, allCodecs[1], h, []*Record{rec})
	// Aux bytes are the last 8+8+4+1+8+8 = 37 bytes before 5WOLB (spec Table 8).
	aux := data[len(data)-5-37 : len(data)-5]
	want := []byte{}
	want = binary.LittleEndian.AppendUint64(want, 0)                  // char* length 0
	want = binary.LittleEndian.AppendUint64(want, 0x7ff8000000000000) // NaN
	want = binary.LittleEndian.AppendUint32(want, math.MaxInt32)
	want = append(want, math.MaxUint8)
	want = binary.LittleEndian.AppendUint64(want, math.MaxUint64)
	want = binary.LittleEndian.AppendUint64(want, 0) // array length 0
	if !bytes.Equal(aux, want) {
		t.Fatalf("missing aux bytes\n got % x\nwant % x", aux, want)
	}
	_, got := readAllBytes(t, data)
	if len(got[0].Aux) != 0 {
		t.Fatalf("missing values read back as %v", got[0].Aux)
	}
	// And in ASCII they are ".".
	ascii := string(writeAll(t, allCodecs[0], h, []*Record{rec}))
	if !strings.HasSuffix(ascii, "\t.\t.\t.\t.\t.\t.\n") {
		t.Fatalf("ASCII missing values: %q", ascii[strings.LastIndex(ascii, "r\t"):])
	}
}

func TestZstdFramesCarryContentSize(t *testing.T) {
	// slow5lib decodes with ZSTD_getFrameContentSize and fails when the
	// frame header lacks Frame_Content_Size. A tiny record is the case
	// klauspost would otherwise encode without it.
	h := Header{ReadGroups: []ReadGroup{{"run_id": "x"}}}
	for _, n := range []int{0, 1, 100, 600, 70000} {
		rec := &Record{ReadID: "r", Digitisation: 1, Range: 1, SamplingRate: 1, RawSignal: smoothSignal(testRand(6), n)}
		data := writeAll(t, codec{"zstd", false, BLOW5Options{RecordZstd, SignalNone}}, h, []*Record{rec})
		size := binary.LittleEndian.Uint32(data[64:])
		frame := data[68+int(size)+8:]
		if binary.LittleEndian.Uint32(frame) != 0xFD2FB528 {
			t.Fatalf("n=%d: not a zstd frame", n)
		}
		fhd := frame[4]
		singleSegment := fhd>>5&1 == 1
		fcsFlag := fhd >> 6
		if !singleSegment && fcsFlag == 0 {
			t.Errorf("n=%d: zstd frame has no Frame_Content_Size", n)
		}
	}
}

func TestVersionBump(t *testing.T) {
	h := Header{Version: Version010, ReadGroups: []ReadGroup{{"run_id": "x"}}}
	for _, c := range allCodecs[1:] {
		data := writeAll(t, c, h, nil)
		want := Version010
		if c.opts.SignalCompression != SignalNone || c.opts.RecordCompression == RecordZstd {
			want = Version020 // slow5lib slow5_press_version_bump
		}
		if got := (Version{data[6], data[7], data[8]}); got != want {
			t.Errorf("%s: version %s, want %s", c.name, got, want)
		}
	}
}

// ------------------------------------------------------------ golden files

func fixture(t *testing.T, name string) []byte {
	t.Helper()
	b, err := os.ReadFile(filepath.Join("testdata", name))
	if err != nil {
		t.Fatal(err)
	}
	if strings.HasSuffix(name, ".gz") {
		zr, err := gzip.NewReader(bytes.NewReader(b))
		if err != nil {
			t.Fatal(err)
		}
		if b, err = io.ReadAll(zr); err != nil {
			t.Fatal(err)
		}
	}
	return b
}

func TestGoldenSlow5libFiles(t *testing.T) {
	ref := fixture(t, "exp_1_lossless.blow5")
	refH, refRecs := readAllBytes(t, ref)
	if refH.Version != Version010 || len(refRecs) != 1 {
		t.Fatalf("reference: version %s, %d records", refH.Version, len(refRecs))
	}
	r := refRecs[0]
	// Values from slow5tools view of the same file.
	if r.ReadID != "a649a4ae-c43d-492a-b6a1-a5b8b8076be4" || r.ReadGroup != 0 ||
		r.Digitisation != 8192 || r.Offset != 23 || r.SamplingRate != 4000 ||
		len(r.RawSignal) != 59676 || r.RawSignal[0] != 1039 || r.RawSignal[1] != 588 {
		t.Fatalf("reference record: %+v (signal %d)", r.ReadID, len(r.RawSignal))
	}
	a, ok := ONTAuxOf(r)
	if !ok || a.ChannelNumber != "115" || a.ReadNumber != 222 || a.StartMux != 1 || a.StartTime != 2817564 ||
		formatDouble(a.MedianBefore, 64) != "225.114" {
		t.Fatalf("reference aux: %+v ok=%v", a, ok)
	}
	if refH.ReadGroups[0]["run_id"] != "d6e473a6d513ec6bfc150c60fd4556d72f0e6d18" ||
		refH.ReadGroups[0]["flow_cell_id"] != "FAB43577" || len(refH.ReadGroups[0]) != 32 {
		t.Fatalf("reference read group: %v", refH.ReadGroups[0])
	}
	if d := headerDiff(Header{Aux: ONTAuxFields, ReadGroups: refH.ReadGroups}, refH); d != "" {
		t.Fatalf("reference aux schema: %s", d)
	}

	for _, tc := range []struct {
		file string
		rc   RecordCompression
		sc   SignalCompression
	}{
		{"lossless_none_svb-zd.blow5", RecordNone, SignalSvbZd},
		{"lossless_zlib_svb-zd.blow5", RecordZlib, SignalSvbZd},
		{"lossless_zstd_svb-zd.blow5", RecordZstd, SignalSvbZd},
	} {
		data := fixture(t, tc.file)
		rd, err := NewReader(bytes.NewReader(data))
		if err != nil {
			t.Fatalf("%s: %v", tc.file, err)
		}
		if rd.Format() != FormatBLOW5 || rd.RecordCompression() != tc.rc || rd.SignalCompression() != tc.sc {
			t.Errorf("%s: format %s %s/%s", tc.file, rd.Format(), rd.RecordCompression(), rd.SignalCompression())
		}
		rd.Close()
		h, recs := readAllBytes(t, data)
		// slow5tools bumps 0.1.0 to 0.2.0 for signal compression.
		if h.Version != Version020 {
			t.Errorf("%s: version %s, want 0.2.0", tc.file, h.Version)
		}
		h.Version = refH.Version
		if d := headerDiff(refH, h); d != "" {
			t.Errorf("%s header: %s", tc.file, d)
		}
		if len(recs) != 1 {
			t.Fatalf("%s: %d records", tc.file, len(recs))
		}
		if d := recordDiff(r, recs[0]); d != "" {
			t.Errorf("%s: %s", tc.file, d)
		}
	}

	// slow5tools' ASCII rendering decodes to the same record. slow5lib
	// prints doubles with %f, so this only holds because these values
	// survive six decimals; compare the doubles to that precision.
	h, recs := readAllBytes(t, fixture(t, "lossless.slow5.gz"))
	if d := headerDiff(refH, h); d != "" {
		t.Errorf("ascii header: %s", d)
	}
	got := recs[0]
	for _, p := range [][2]float64{{got.Range, r.Range}, {got.Aux["median_before"].(float64), r.Aux["median_before"].(float64)}} {
		if math.Abs(p[0]-p[1]) > 5e-7 {
			t.Errorf("ascii double %v vs binary %v", p[0], p[1])
		}
	}
	got.Range, got.Aux["median_before"] = r.Range, r.Aux["median_before"]
	if d := recordDiff(r, got); d != "" {
		t.Errorf("ascii: %s", d)
	}
}

func TestGoldenRewriteReproducesBytes(t *testing.T) {
	// Uncompressed and svb-zd files must be reproduced byte for byte; the
	// header (sorted attributes, padding, sizes), record layout, svb-zd
	// codec and EOF marker all have to match slow5lib for this to pass.
	for _, tc := range []struct {
		file string
		opts BLOW5Options
	}{
		{"exp_1_lossless.blow5", BLOW5Options{RecordNone, SignalNone}},
		{"lossless_none_svb-zd.blow5", BLOW5Options{RecordNone, SignalSvbZd}},
	} {
		orig := fixture(t, tc.file)
		h, recs := readAllBytes(t, orig)
		got := writeAll(t, codec{tc.file, false, tc.opts}, h, recs)
		if !bytes.Equal(got, orig) {
			t.Errorf("%s: rewrite differs (got %d bytes, want %d; first diff at %d)", tc.file, len(got), len(orig), firstDiff(got, orig))
		}
	}
	// SLOW5 ASCII written by slow5tools is reproduced as well: the double
	// formatting matches slow5lib's %f form when that form is exact.
	orig := fixture(t, "lossless.slow5.gz")
	h, recs := readAllBytes(t, orig)
	if got := writeAll(t, allCodecs[0], h, recs); !bytes.Equal(got, orig) {
		t.Errorf("slow5 ASCII rewrite differs at byte %d", firstDiff(got, orig))
	}
	// zlib and zstd output is not byte-identical to slow5lib (different
	// compressor implementations and levels), but decodes identically and
	// the uncompressed records inside are identical.
	for _, tc := range []struct {
		file string
		opts BLOW5Options
	}{
		{"lossless_zlib_svb-zd.blow5", BLOW5Options{RecordZlib, SignalSvbZd}},
		{"lossless_zstd_svb-zd.blow5", BLOW5Options{RecordZstd, SignalSvbZd}},
	} {
		orig := fixture(t, tc.file)
		h, recs := readAllBytes(t, orig)
		_, again := readAllBytes(t, writeAll(t, codec{tc.file, false, tc.opts}, h, recs))
		if d := recordDiff(recs[0], again[0]); d != "" {
			t.Errorf("%s: rewrite decodes differently: %s", tc.file, d)
		}
	}
}

func firstDiff(a, b []byte) int {
	for i := 0; i < len(a) && i < len(b); i++ {
		if a[i] != b[i] {
			return i
		}
	}
	return min(len(a), len(b))
}

// ------------------------------------------------------------ errors

func TestReaderRejectsMalformedInput(t *testing.T) {
	h := Header{ReadGroups: []ReadGroup{{"run_id": "x"}}, Aux: ONTAuxFields}
	rec := &Record{ReadID: "r", Digitisation: 1, Range: 1, SamplingRate: 1, RawSignal: []int16{1, 2, 3}, Aux: ONTAux{ChannelNumber: "1"}.Map()}
	good := writeAll(t, allCodecs[2], h, []*Record{rec})
	hdrEnd := 68 + int(binary.LittleEndian.Uint32(good[64:]))
	mutate := func(f func(b []byte) []byte) []byte { return f(append([]byte(nil), good...)) }
	bad := map[string][]byte{
		"no-eof":           good[:len(good)-5],
		"partial-eof":      good[:len(good)-2],
		"after-eof":        append(append([]byte(nil), good...), 'x'),
		"truncated-record": append(append([]byte(nil), good[:len(good)-9]...), blow5EOF...),
		"bad-magic":        mutate(func(b []byte) []byte { b[5] = 2; return b }),
		"newer-version":    mutate(func(b []byte) []byte { b[7] = 3; return b }),
		"zero-groups":      mutate(func(b []byte) []byte { b[10] = 0; return b }),
		"huge-groups":      mutate(func(b []byte) []byte { b[13] = 0xff; return b }),
		"ex-zd":            mutate(func(b []byte) []byte { b[14] = 2; return b }),
		"rec-codec-250":    mutate(func(b []byte) []byte { b[9] = 250; return b }),
		"header-size":      mutate(func(b []byte) []byte { b[64]--; return b }),
		"record-size-big":  mutate(func(b []byte) []byte { b[hdrEnd]++; return b }),
		"record-size-small": mutate(func(b []byte) []byte {
			b[hdrEnd]--
			return b
		}),
		"read-group-range": mutate(func(b []byte) []byte { b[hdrEnd+8+2+1] = 5; return b }),
		"empty":            {},
		"text":             []byte("hello world\n"),
	}
	for name, b := range bad {
		if _, _, err := ReadAll(bytes.NewReader(b)); err == nil {
			t.Errorf("%s: malformed BLOW5 accepted", name)
		}
	}
	ascii := string(writeAll(t, allCodecs[0], h, []*Record{rec}))
	badASCII := map[string]string{
		"no-newline":   strings.TrimSuffix(ascii, "\n"),
		"short-signal": strings.Replace(ascii, "\t3\t1,2,3\t", "\t4\t1,2,3\t", 1),
		"long-signal":  strings.Replace(ascii, "\t3\t1,2,3\t", "\t2\t1,2,3\t", 1),
		"bad-sample":   strings.Replace(ascii, "1,2,3", "1,2,40000", 1),
		"few-fields":   strings.Replace(ascii, "\t1\t", "\t", 1),
		"bad-version":  strings.Replace(ascii, "0.2.0", "0.3.0", 1),
		"bad-groups":   strings.Replace(ascii, "#num_read_groups\t1", "#num_read_groups\t2", 1),
		"bad-type":     strings.Replace(ascii, "uint8_t", "uint9_t", 1),
	}
	for name, s := range badASCII {
		if _, _, err := ReadAll(strings.NewReader(s)); err == nil {
			t.Errorf("%s: malformed SLOW5 accepted", name)
		}
	}
}

func TestWriterRejectsUnrepresentableData(t *testing.T) {
	h := Header{ReadGroups: []ReadGroup{{"run_id": "x"}}, Aux: allAuxFields}
	base := func() *Record {
		return &Record{ReadID: "r", Digitisation: 1, Range: 1, SamplingRate: 1}
	}
	cases := map[string]func(r *Record){
		"empty-id":         func(r *Record) { r.ReadID = "" },
		"tab-id":           func(r *Record) { r.ReadID = "a\tb" },
		"hash-id":          func(r *Record) { r.ReadID = "#a" },
		"long-id":          func(r *Record) { r.ReadID = strings.Repeat("a", 70000) },
		"group-range":      func(r *Record) { r.ReadGroup = 1 },
		"unknown-aux":      func(r *Record) { r.Aux = map[string]any{"nope": int32(1)} },
		"wrong-go-type":    func(r *Record) { r.Aux = map[string]any{"read_number": 5} },
		"sentinel-int32":   func(r *Record) { r.Aux = map[string]any{"read_number": int32(math.MaxInt32)} },
		"sentinel-uint8":   func(r *Record) { r.Aux = map[string]any{"start_mux": uint8(255)} },
		"sentinel-nan":     func(r *Record) { r.Aux = map[string]any{"median_before": math.NaN()} },
		"sentinel-in-arr":  func(r *Record) { r.Aux = map[string]any{"a_u16": []uint16{1, 65535}} },
		"empty-string":     func(r *Record) { r.Aux = map[string]any{"channel_number": ""} },
		"dot-string":       func(r *Record) { r.Aux = map[string]any{"channel_number": "."} },
		"newline-string":   func(r *Record) { r.Aux = map[string]any{"channel_number": "1\n2"} },
		"empty-array":      func(r *Record) { r.Aux = map[string]any{"a_f64": []float64{}} },
		"enum-no-label":    func(r *Record) { r.Aux = map[string]any{"end_reason": Enum(6)} },
		"char-nul":         func(r *Record) { r.Aux = map[string]any{"ch": Char(0)} },
		"char-dot":         func(r *Record) { r.Aux = map[string]any{"ch": Char('.')} },
		"uint8-as-char":    func(r *Record) { r.Aux = map[string]any{"ch": uint8('A')} },
		"char-as-uint8":    func(r *Record) { r.Aux = map[string]any{"start_mux": Char(1)} },
		"array-wrong-type": func(r *Record) { r.Aux = map[string]any{"a_i32": []int64{1}} },
	}
	for _, c := range []codec{allCodecs[0], allCodecs[6]} {
		var buf bytes.Buffer
		var w Writer
		var err error
		if c.ascii {
			w, err = NewSLOW5Writer(&buf, h)
		} else {
			w, err = NewBLOW5Writer(&buf, h, c.opts)
		}
		if err != nil {
			t.Fatal(err)
		}
		for name, mut := range cases {
			r := base()
			mut(r)
			if err := w.Write(r); !errors.Is(err, ErrFormat) {
				t.Errorf("%s/%s: got %v, want ErrFormat", c.name, name, err)
			}
		}
		// Rejected records leave the writer usable and the file valid.
		if err := w.Write(base()); err != nil {
			t.Fatal(err)
		}
		if err := w.Close(); err != nil {
			t.Fatal(err)
		}
		if _, recs := readAllBytes(t, buf.Bytes()); len(recs) != 1 {
			t.Fatalf("%s: %d records after rejected writes", c.name, len(recs))
		}
	}

	badHeaders := map[string]Header{
		"no-groups":     {},
		"empty-value":   {ReadGroups: []ReadGroup{{"run_id": ""}}},
		"dot-value":     {ReadGroups: []ReadGroup{{"run_id": "."}}},
		"tab-key":       {ReadGroups: []ReadGroup{{"a\tb": "x"}}},
		"dup-aux":       {ReadGroups: []ReadGroup{{}}, Aux: []AuxField{{Name: "a", Type: AuxInt8}, {Name: "a", Type: AuxInt8}}},
		"primary-name":  {ReadGroups: []ReadGroup{{}}, Aux: []AuxField{{Name: "offset", Type: AuxInt8}}},
		"bad-type":      {ReadGroups: []ReadGroup{{}}, Aux: []AuxField{{Name: "a"}}},
		"enum-nolabels": {ReadGroups: []ReadGroup{{}}, Aux: []AuxField{{Name: "a", Type: AuxEnum}}},
		"enum-badlabel": {ReadGroups: []ReadGroup{{}}, Aux: []AuxField{{Name: "a", Type: AuxEnum, EnumLabels: []string{"x y"}}}},
		"enum-v010":     {Version: Version010, ReadGroups: []ReadGroup{{}}, Aux: []AuxField{{Name: "a", Type: AuxEnum, EnumLabels: []string{"x"}}}},
		"version-new":   {Version: Version{0, 3, 0}, ReadGroups: []ReadGroup{{}}},
	}
	for name, bh := range badHeaders {
		if _, err := NewBLOW5Writer(io.Discard, bh, BLOW5Options{}); !errors.Is(err, ErrFormat) {
			t.Errorf("%s: NewBLOW5Writer: %v, want ErrFormat", name, err)
		}
		if _, err := NewSLOW5Writer(io.Discard, bh); !errors.Is(err, ErrFormat) {
			t.Errorf("%s: NewSLOW5Writer: %v, want ErrFormat", name, err)
		}
	}
	if _, err := NewBLOW5Writer(io.Discard, h, BLOW5Options{RecordCompression: 3}); err == nil {
		t.Error("unknown record compression accepted")
	}
	if _, err := NewBLOW5Writer(io.Discard, h, BLOW5Options{SignalCompression: 2}); err == nil {
		t.Error("ex-zd signal compression accepted")
	}
}

func TestWriterDoesNotAliasCallerData(t *testing.T) {
	h := Header{ReadGroups: []ReadGroup{{"run_id": "x"}}}
	var buf bytes.Buffer
	w, err := NewBLOW5Writer(&buf, h, BLOW5Options{RecordZstd, SignalSvbZd})
	if err != nil {
		t.Fatal(err)
	}
	h.ReadGroups[0]["run_id"] = "changed" // after the header was written
	sig := []int16{1, 2, 3}
	rec := &Record{ReadID: "r", Digitisation: 1, Range: 1, SamplingRate: 1, RawSignal: sig}
	if err := w.Write(rec); err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(rec.RawSignal, []int16{1, 2, 3}) {
		t.Fatal("Write modified the record's signal")
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	if err := w.Write(rec); err == nil {
		t.Error("write after Close accepted")
	}
	gotH, _ := readAllBytes(t, buf.Bytes())
	if gotH.ReadGroups[0]["run_id"] != "x" {
		t.Error("header changed after NewBLOW5Writer affected the file")
	}
}

type failWriter struct{ n int }

func (f *failWriter) Write(p []byte) (int, error) {
	if f.n -= len(p); f.n < 0 {
		return 0, errors.New("disk full")
	}
	return len(p), nil
}

func TestWriterReportsIOErrors(t *testing.T) {
	h := Header{ReadGroups: []ReadGroup{{"run_id": "x"}}}
	// The header is buffered, so the failure surfaces when a 2 MB record
	// overflows the 1 MB buffer; it must then stick.
	w, err := NewBLOW5Writer(&failWriter{n: 10}, h, BLOW5Options{})
	if err != nil {
		t.Fatal(err)
	}
	rec := &Record{ReadID: "r", Digitisation: 1, Range: 1, SamplingRate: 1, RawSignal: make([]int16, 1<<20)}
	if err := w.Write(rec); err == nil {
		t.Fatal("Write did not report the I/O error")
	}
	small := &Record{ReadID: "s", Digitisation: 1, Range: 1, SamplingRate: 1}
	if err := w.Write(small); err == nil {
		t.Error("Write after an I/O error succeeded")
	}
	if err := w.Close(); err == nil {
		t.Error("Close after an I/O error succeeded")
	}
	// Small output that fails only at the final flush is reported by Close.
	w, err = NewBLOW5Writer(&failWriter{n: 10}, h, BLOW5Options{})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err == nil {
		t.Error("Close did not report the flush error")
	}
}

// ------------------------------------------------------------ formatting

func TestFormatDouble(t *testing.T) {
	cases := []struct {
		x    float64
		bits int
		want string
	}{
		{8192, 64, "8192"},
		{4000, 64, "4000"},
		{1467.61, 64, "1467.61"},                     // slow5lib %f form, exact
		{225.114, 64, "225.114"},                     //
		{0.1, 64, "0.1"},                             //
		{-23, 64, "-23"},                             //
		{math.Copysign(0, -1), 64, "0"},              // slow5lib strips "-0"
		{1441.389892578125, 64, "1441.389892578125"}, // %f would round; keep exact
		{1e-9, 64, "0.000000001"},                    // %f gives 0; keep exact
		{1e22, 64, "10000000000000000000000"},
		{float64(float32(0.1)), 32, "0.1"},
		{float64(float32(238.78226)), 32, "238.782257"}, // slow5lib %f of the float
		{math.NaN(), 64, "nan"},
		{math.Inf(-1), 64, "-inf"},
	}
	for _, c := range cases {
		if got := formatDouble(c.x, c.bits); got != c.want {
			t.Errorf("formatDouble(%v, %d) = %q, want %q", c.x, c.bits, got, c.want)
		}
	}
	r := testRand(7)
	for i := 0; i < 10000; i++ {
		x := math.Float64frombits(r.Uint64())
		if math.IsNaN(x) || math.IsInf(x, 0) {
			continue
		}
		s := formatDouble(x, 64)
		if strings.ContainsAny(s, "eE") {
			t.Fatalf("formatDouble(%v) = %q has an exponent", x, s)
		}
		var y float64
		if _, err := fmt.Sscan(s, &y); err != nil || y != x {
			t.Fatalf("formatDouble(%v) = %q does not round-trip", x, s)
		}
	}
}

// ------------------------------------------------------------ chunking

func chunkRecords(n int) []*Record {
	r := testRand(8)
	recs := make([]*Record, n)
	for i := range recs {
		// UUID-form IDs: pod5 (blue-crab s2p) parses read IDs as UUIDs.
		recs[i] = &Record{ReadID: fmt.Sprintf("6a1f7c2e-4b9d-4e3a-8c5f-%012x", i), Digitisation: 2048, Offset: -240, Range: 748.5, SamplingRate: 5000,
			RawSignal: smoothSignal(r, 50+r.IntN(500)),
			Aux:       ONTAux{ChannelNumber: fmt.Sprint(1 + i%512), MedianBefore: 210.25, ReadNumber: int32(i), StartMux: uint8(1 + i%4), StartTime: uint64(i) * 5000}.Map()}
	}
	return recs
}

func writeChunks(t *testing.T, dir string, opts ChunkOptions, recs []*Record) *ChunkedWriter {
	t.Helper()
	h := Header{ReadGroups: []ReadGroup{{"run_id": "8f1c2b3d", "flow_cell_id": "PAW12345"}}, Aux: ONTAuxFields}
	cw, err := NewChunkedWriter(dir, h, opts)
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range recs {
		if err := cw.Write(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := cw.Close(); err != nil {
		t.Fatal(err)
	}
	return cw
}

func TestChunkedWriter(t *testing.T) {
	recs := chunkRecords(10)
	dir := t.TempDir()
	cw := writeChunks(t, dir, ChunkOptions{ReadsPerFile: 4, BLOW5: BLOW5Options{RecordZstd, SignalSvbZd}}, recs)
	wantNames := []string{"PAW12345_8f1c2b3d_0.blow5", "PAW12345_8f1c2b3d_1.blow5", "PAW12345_8f1c2b3d_2.blow5"}
	files := cw.Files()
	if len(files) != len(wantNames) || cw.Records() != 10 {
		t.Fatalf("files %v, records %d", files, cw.Records())
	}
	var all []*Record
	for i, f := range files {
		if filepath.Base(f) != wantNames[i] {
			t.Errorf("file %d named %s, want %s", i, filepath.Base(f), wantNames[i])
		}
		b, err := os.ReadFile(f)
		if err != nil {
			t.Fatal(err)
		}
		h, got := readAllBytes(t, b)
		if h.ReadGroups[0]["run_id"] != "8f1c2b3d" || len(h.Aux) != 5 {
			t.Errorf("%s: header not complete: %+v", f, h)
		}
		if want := []int{4, 4, 2}[i]; len(got) != want {
			t.Errorf("%s: %d reads, want %d", f, len(got), want)
		}
		all = append(all, got...)
	}
	for i := range recs {
		if d := recordDiff(recs[i], all[i]); d != "" {
			t.Errorf("read %d: %s", i, d)
		}
	}
	entries, _ := os.ReadDir(dir)
	if len(entries) != 3 {
		t.Errorf("directory has %d entries, want 3", len(entries))
	}

	// Deterministic: a second run gives byte-identical files.
	dir2 := t.TempDir()
	cw2 := writeChunks(t, dir2, ChunkOptions{ReadsPerFile: 4, BLOW5: BLOW5Options{RecordZstd, SignalSvbZd}}, recs)
	for i := range files {
		a, _ := os.ReadFile(files[i])
		b, _ := os.ReadFile(cw2.Files()[i])
		if !bytes.Equal(a, b) {
			t.Errorf("chunk %d differs between identical runs", i)
		}
	}

	// A rerun into the same directory fails instead of mixing chunks...
	h := Header{ReadGroups: []ReadGroup{{"run_id": "8f1c2b3d", "flow_cell_id": "PAW12345"}}, Aux: ONTAuxFields}
	cw3, err := NewChunkedWriter(dir, h, ChunkOptions{ReadsPerFile: 4})
	if err != nil {
		t.Fatal(err)
	}
	if err := cw3.Write(recs[0]); err == nil {
		t.Error("existing chunk overwritten without Overwrite")
	}
	_ = cw3.Close()
	// ...unless Overwrite is set.
	writeChunks(t, dir, ChunkOptions{ReadsPerFile: 4, Overwrite: true}, recs[:2])
}

func TestChunkedWriterEdges(t *testing.T) {
	recs := chunkRecords(8)
	// Exact multiple: no empty trailing file.
	cw := writeChunks(t, t.TempDir(), ChunkOptions{ReadsPerFile: 4}, recs)
	if len(cw.Files()) != 2 {
		t.Errorf("8 reads at 4/file gave %d files", len(cw.Files()))
	}
	// No records: no files.
	cw = writeChunks(t, t.TempDir(), ChunkOptions{}, nil)
	if len(cw.Files()) != 0 {
		t.Errorf("empty run wrote %v", cw.Files())
	}
	// Default size.
	cw = writeChunks(t, t.TempDir(), ChunkOptions{}, recs)
	if len(cw.Files()) != 1 || DefaultReadsPerFile != 4000 {
		t.Errorf("default chunking: %v", cw.Files())
	}
	// Explicit IDs override the header; unsafe names are rejected.
	cw = writeChunks(t, t.TempDir(), ChunkOptions{FlowCellID: "FAX99999", RunID: "abc"}, recs[:1])
	if got := filepath.Base(cw.Files()[0]); got != "FAX99999_abc_0.blow5" {
		t.Errorf("name %s", got)
	}
	h := Header{ReadGroups: []ReadGroup{{"run_id": "r"}}}
	for _, o := range []ChunkOptions{{}, {FlowCellID: "../x", RunID: "r"}, {FlowCellID: "a/b"}, {FlowCellID: "F", ReadsPerFile: -1}} {
		if _, err := NewChunkedWriter(t.TempDir(), h, o); err == nil {
			t.Errorf("options %+v accepted", o)
		}
	}
	if got := ChunkFileName("PAW1", "run", 399); got != "PAW1_run_399.blow5" {
		t.Errorf("ChunkFileName = %s", got)
	}
}
