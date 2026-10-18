//go:build !release
// +build !release

package fingerprinting

import (
	"testing"
	"time"
)

func TestFingerprintRoundTrip(t *testing.T) {
	created := time.Date(2026, 10, 18, 4, 33, 0, 0, time.UTC)
	for _, tc := range []struct{ customer, dataset string }{
		{"CUST_1", "ds-1"},
		{"acme", "run:2026:batch-7"}, // dataset IDs may contain ':'
	} {
		s := FormatFingerprint(tc.customer, tc.dataset, "0.8.0", created)
		fp, err := ParseFingerprint(s)
		if err != nil {
			t.Fatalf("ParseFingerprint(%q): %v", s, err)
		}
		want := Fingerprint{Customer: tc.customer, Dataset: tc.dataset, SteganoVersion: "0.8.0", Created: created}
		if fp != want {
			t.Errorf("ParseFingerprint(%q) = %+v, want %+v", s, fp, want)
		}
	}
}

func TestParseFingerprint_Rejects(t *testing.T) {
	for _, s := range []string{
		"",
		"cust:ds",
		"cust:ds:v0.8.0:ts-1",              // version field missing "stegano-" prefix
		"cust:ds:stegano-0.8.0:1792288985", // timestamp missing "ts-" prefix
		"cust:ds:stegano-0.8.0:ts-soon",
		"cust:ds:1234567890", // legacy engine format
	} {
		if fp, err := ParseFingerprint(s); err == nil {
			t.Errorf("ParseFingerprint(%q) = %+v, want error", s, fp)
		}
	}
}
