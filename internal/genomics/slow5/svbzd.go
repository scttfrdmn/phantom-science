package slow5

import (
	"encoding/binary"
	"math"
)

// svb-zd is the BLOW5 signal compression with code 1: "StreamVByte with
// zig-zag delta encoding" (spec, BLOW5 Header, Table 6 footnote 2; Lemire,
// Kurz and Rupp, "Stream VByte: Faster byte-oriented integer compression",
// Information Processing Letters 130 (2018) 1-6).
//
// The specification does not define the byte layout; it is fixed by
// slow5lib (src/slow5_press.c ptr_compress_svb_zd / ptr_compress_svb, and
// thirdparty/streamvbyte), which this follows exactly:
//
//  1. Each int16 sample is widened to int32 and delta-coded against the
//     previous sample, starting from 0: d[i] = s[i] - s[i-1] (s[-1] = 0).
//  2. Each delta is zig-zag coded to uint32: z = (d << 1) ^ (d >> 31).
//  3. The uint32 values are StreamVByte-coded (the standard 1-to-4-byte
//     variant, not streamvbyte_encode_0124): ceil(n/4) control bytes with
//     2 bits per value holding (byte length - 1), value i in bits
//     2*(i%4) of control byte i/4, then each value's low bytes,
//     little-endian.
//  4. The output is prefixed with n, the number of samples, as a
//     little-endian uint32 ("copy original length of ptr (needed for
//     depress)").
//
// In a BLOW5 record with svb-zd the len_raw_signal field holds the size of
// this compressed blob in bytes, not the number of samples (slow5lib
// slow5_rec_to_mem sets len_raw_signal = compressed bytes; slow5_rec_parse
// reads that many bytes and recovers the sample count from the prefix).
// Deltas of int16 samples fit in 17 bits, so no value takes 4 bytes.

// svbZdMaxLen is the largest sample count the uint32 prefix can hold.
const svbZdMaxLen = math.MaxUint32

// svbZdEncode appends the svb-zd encoding of sig to dst.
func svbZdEncode(dst []byte, sig []int16) []byte {
	n := len(sig)
	keyLen := (n + 3) / 4
	dst = binary.LittleEndian.AppendUint32(dst, uint32(n))
	keyAt := len(dst)
	// Control bytes first, then up to 3 data bytes per sample.
	need := keyLen + 3*n
	if cap(dst)-len(dst) < need {
		grown := make([]byte, len(dst), len(dst)+need)
		copy(grown, dst)
		dst = grown
	}
	dst = dst[:keyAt+keyLen]
	for i := keyAt; i < keyAt+keyLen; i++ {
		dst[i] = 0
	}
	var prev int32
	for i, s := range sig {
		d := int32(s) - prev
		prev = int32(s)
		z := uint32(d<<1) ^ uint32(d>>31)
		var code byte
		switch {
		case z < 1<<8:
			dst = append(dst, byte(z))
		case z < 1<<16:
			dst = append(dst, byte(z), byte(z>>8))
			code = 1
		case z < 1<<24:
			dst = append(dst, byte(z), byte(z>>8), byte(z>>16))
			code = 2
		default: // unreachable for int16 input; kept for format completeness
			dst = append(dst, byte(z), byte(z>>8), byte(z>>16), byte(z>>24))
			code = 3
		}
		dst[keyAt+i/4] |= code << (2 * uint(i%4))
	}
	return dst
}

// svbZdDecode decodes a complete svb-zd blob. Like slow5lib it requires the
// StreamVByte data to consume exactly the blob ("Expected
// streamvbyte_decode to read ... bytes"); it also rejects samples outside
// the int16 range, which no encoder of int16 data produces.
func svbZdDecode(b []byte) ([]int16, error) {
	if len(b) < 4 {
		return nil, formatErr("svb-zd: blob of %d bytes has no length prefix", len(b))
	}
	n64 := uint64(binary.LittleEndian.Uint32(b))
	b = b[4:]
	keyLen := (n64 + 3) / 4
	// Every value takes at least one data byte: bound n before allocating.
	if keyLen+n64 > uint64(len(b)) {
		return nil, formatErr("svb-zd: %d samples cannot fit in %d bytes", n64, len(b))
	}
	n := int(n64)
	keys, data := b[:keyLen], b[keyLen:]
	out := make([]int16, n)
	var prev int32
	p := 0
	for i := 0; i < n; i++ {
		code := (keys[i/4] >> (2 * uint(i%4))) & 3
		w := int(code) + 1
		if p+w > len(data) {
			return nil, formatErr("svb-zd: data truncated at sample %d", i)
		}
		var z uint32
		for k := 0; k < w; k++ {
			z |= uint32(data[p+k]) << (8 * uint(k))
		}
		p += w
		d := int32(z>>1) ^ -int32(z&1)
		prev += d
		if prev < math.MinInt16 || prev > math.MaxInt16 {
			return nil, formatErr("svb-zd: sample %d decodes to %d, outside int16", i, prev)
		}
		out[i] = int16(prev)
	}
	if p != len(data) {
		return nil, formatErr("svb-zd: %d trailing bytes after %d samples", len(data)-p, n)
	}
	return out, nil
}
