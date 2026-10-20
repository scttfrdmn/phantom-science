//go:build !release

package fingerprinting

import (
	"errors"
	"strings"
	"testing"
)

const seedTestKey = "test-key-0123456789abcdef" // gitleaks:allow -- test fixture, not a real key

// Known answer, computed independently with Python's hmac module:
//
//	hmac.new(b"test-key-0123456789abcdef",
//	         b"phantom-science/fingerprint/v1\x00ACME\x00study_1\x00" + struct.pack(">Q", 42),
//	         hashlib.sha256).digest()[:8]  # 39f899af64bcec57
func TestDeriveSeed_KnownAnswer(t *testing.T) {
	got, err := DeriveSeed([]byte(seedTestKey), "ACME", "study_1", 42)
	if err != nil {
		t.Fatal(err)
	}
	const want = uint64(0x39f899af64bcec57)
	if got != want {
		t.Fatalf("DeriveSeed = %#x, want %#x", got, want)
	}
}

// Every input changes the seed, including moving a byte across the
// customer/dataset boundary (the NUL separators make the encoding
// unambiguous).
func TestDeriveSeed_EveryInputMatters(t *testing.T) {
	base, _ := DeriveSeed([]byte(seedTestKey), "ACME", "study_1", 42)
	cases := []struct {
		name     string
		key      string
		cust, ds string
		salt     uint64
	}{
		{"key", "other-key-0123456789abcdef", "ACME", "study_1", 42}, // gitleaks:allow -- test fixture
		{"customer", seedTestKey, "ACMF", "study_1", 42},
		{"dataset", seedTestKey, "ACME", "study_2", 42},
		{"salt", seedTestKey, "ACME", "study_1", 43},
		{"boundary", seedTestKey, "ACMEs", "tudy_1", 42},
	}
	for _, tc := range cases {
		got, err := DeriveSeed([]byte(tc.key), tc.cust, tc.ds, tc.salt)
		if err != nil {
			t.Fatal(err)
		}
		if got == base {
			t.Errorf("%s: changing it left the seed at %#x", tc.name, got)
		}
	}
	again, _ := DeriveSeed([]byte(seedTestKey), "ACME", "study_1", 42)
	if again != base {
		t.Error("DeriveSeed is not deterministic")
	}
}

func TestDeriveSeed_Rejects(t *testing.T) {
	for name, f := range map[string]func() error{
		"short key":      func() error { _, err := DeriveSeed([]byte("short"), "A", "d", 1); return err },
		"empty customer": func() error { _, err := DeriveSeed([]byte(seedTestKey), "", "d", 1); return err },
		"empty dataset":  func() error { _, err := DeriveSeed([]byte(seedTestKey), "A", "", 1); return err },
		"NUL":            func() error { _, err := DeriveSeed([]byte(seedTestKey), "A\x00B", "d", 1); return err },
	} {
		if f() == nil {
			t.Errorf("%s: no error", name)
		}
	}
}

func TestLayers(t *testing.T) {
	for _, tc := range []struct {
		layers []string
		keyed  bool
		err    bool
	}{
		{nil, true, false},
		{[]string{"keyed_seed"}, true, false},
		{[]string{"read_ids"}, false, false},
		{[]string{"read_ids", "keyed_seed"}, true, false},
		{[]string{"quality_scores"}, false, true},
		{[]string{"keyed_seed", "quality_scores"}, false, true},
		{[]string{"headers"}, false, true},
	} {
		keyed, err := Layers(tc.layers)
		if keyed != tc.keyed || (err != nil) != tc.err {
			t.Errorf("Layers(%v) = %v, %v", tc.layers, keyed, err)
		}
	}
	_, err := Layers([]string{"quality_scores"})
	if !errors.Is(err, ErrQualityScoresLayer) || !strings.Contains(err.Error(), "CHARTER") || !strings.Contains(err.Error(), "keyed_seed") {
		t.Errorf("quality_scores error %q should point to keyed_seed and the charter", err)
	}
}
