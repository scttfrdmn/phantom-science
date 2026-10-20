//go:build !release

package seqread

import (
	"bytes"
	"testing"
)

func TestReverseComplement(t *testing.T) {
	for in, want := range map[string]string{"": "", "A": "T", "ACGTN": "NACGT", "AACG": "CGTT", "acgX": "Ncgt"} {
		if got := string(ReverseComplement([]byte(in))); got != want {
			t.Errorf("ReverseComplement(%q) = %q, want %q", in, got, want)
		}
	}
	s := []byte("GATTACA")
	if got := ReverseComplement(ReverseComplement(s)); !bytes.Equal(got, s) {
		t.Errorf("double reverse complement %q, want %q", got, s)
	}
}

func TestSplice(t *testing.T) {
	seq := []byte("ACGTACGT")
	qual := []byte{10, 11, 12, 13, 14, 15, 16, 17}
	cases := []struct {
		name     string
		off, ref int
		allele   string
		wantSeq  string
		wantQual []byte
	}{
		{"snv", 2, 1, "T", "ACTTACGT", []byte{10, 11, 12, 13, 14, 15, 16, 17}},
		{"deletion", 1, 3, "C", "ACACGT", []byte{10, 11, 14, 15, 16, 17}},
		{"insertion", 3, 1, "TGGG", "ACGTGGGACGT", []byte{10, 11, 12, 13, 13, 13, 13, 14, 15, 16, 17}},
		{"at end", 7, 1, "TAA", "ACGTACGTAA", []byte{10, 11, 12, 13, 14, 15, 16, 17, 17, 17}},
		{"mnv", 0, 2, "GG", "GGGTACGT", []byte{10, 11, 12, 13, 14, 15, 16, 17}},
	}
	for _, tc := range cases {
		s, q := Splice(seq, qual, tc.off, tc.ref, []byte(tc.allele))
		if string(s) != tc.wantSeq || !bytes.Equal(q, tc.wantQual) {
			t.Errorf("%s: got %q %v, want %q %v", tc.name, s, q, tc.wantSeq, tc.wantQual)
		}
		if len(s) != len(seq)-tc.ref+len(tc.allele) {
			t.Errorf("%s: length %d", tc.name, len(s))
		}
	}
	if string(seq) != "ACGTACGT" || qual[0] != 10 {
		t.Error("inputs modified")
	}
	if s, q := Splice(seq, nil, 0, 1, []byte("T")); string(s) != "TCGTACGT" || q != nil {
		t.Errorf("nil qualities: %q %v", s, q)
	}
}

func TestWindowSlack(t *testing.T) {
	spans := []Span{
		{Pos: 10, RefLen: 5, AltLen: 1},  // -4, inside [0,20)
		{Pos: 21, RefLen: 3, AltLen: 1},  // -2, pulled in by the first deletion
		{Pos: 30, RefLen: 1, AltLen: 6},  // insertion: no slack
		{Pos: -5, RefLen: 9, AltLen: 1},  // starts before the read: cannot shorten it
		{Pos: 100, RefLen: 9, AltLen: 1}, // far away
	}
	if got := WindowSlack(0, 20, spans); got != 6 {
		t.Errorf("WindowSlack = %d, want 6", got)
	}
	if got := WindowSlack(50, 20, spans); got != 0 {
		t.Errorf("WindowSlack(no variants) = %d, want 0", got)
	}
}

func TestRecordHeaderAndMeanQuality(t *testing.T) {
	r := Record{Name: "r1", Comment: "1:N:0:1", Qual: []byte{10, 20, 30}}
	if r.Header() != "@r1 1:N:0:1" {
		t.Errorf("Header = %q", r.Header())
	}
	r.Comment = ""
	if r.Header() != "@r1" {
		t.Errorf("Header = %q", r.Header())
	}
	if r.MeanQuality() != 20 {
		t.Errorf("MeanQuality = %v", r.MeanQuality())
	}
}
