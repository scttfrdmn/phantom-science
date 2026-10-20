//go:build !release

package fingerprinting

import (
	"slices"
	"strings"
	"testing"
)

func TestReadTags_RoundTrip(t *testing.T) {
	for _, tc := range []struct {
		cust, ds string
		want     []string
	}{
		{"CUST_001", "ds-1", []string{"XC:Z:CUST_001", "XD:Z:ds-1"}},
		{"acme", "", []string{"XC:Z:acme"}},
		{"", "run:2026:7", []string{"XD:Z:run:2026:7"}},
		{"two words", "a\tb", []string{"XC:Z:two_words", "XD:Z:a_b"}},
		{"", "", nil},
	} {
		got := ReadTags(tc.cust, tc.ds)
		if !slices.Equal(got, tc.want) {
			t.Errorf("ReadTags(%q, %q) = %q, want %q", tc.cust, tc.ds, got, tc.want)
		}
		header := "@S:1:HXXXXXX:1:1101:1000:2000 1:N:0:1\tRG:Z:s1"
		if len(got) > 0 {
			header += "\t" + strings.Join(got, "\t")
		}
		c, d := ParseReadTags(header)
		if wc, wd := tagValue(tc.cust), tagValue(tc.ds); c != wc || d != wd {
			t.Errorf("ParseReadTags(%q) = %q, %q; want %q, %q", header, c, d, wc, wd)
		}
	}
}
