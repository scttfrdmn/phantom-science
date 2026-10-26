//go:build !release

package fingerprinting

import (
	"math/rand/v2"
	"strings"
	"testing"
)

func randomSeq(r *rand.Rand, n int) []byte {
	s := make([]byte, n)
	for i := range s {
		s[i] = "ACGT"[r.IntN(4)]
	}
	return s
}

func revComp(s []byte) []byte {
	out := make([]byte, len(s))
	for i, b := range s {
		out[len(s)-1-i] = complement(b)
	}
	return out
}

// Exact, reverse-complement, trimmed, lower-case and lightly edited reads
// match; unrelated reads do not.
func TestMatcher(t *testing.T) {
	r := rand.New(rand.NewPCG(1, 2))
	regen := make([][]byte, 200)
	for i := range regen {
		regen[i] = randomSeq(r, 150)
	}
	m := NewMatcher()
	f := m.AddFile("suspect")
	other := m.AddFile("other")
	m.AddRead(f, regen[0])
	m.AddRead(f, revComp(regen[1]))
	m.AddRead(f, regen[2][10:140])
	m.AddRead(f, []byte(strings.ToLower(string(regen[3]))))
	edited := append([]byte(nil), regen[4]...)
	edited[75] = complement(edited[75])
	m.AddRead(f, edited)
	for range 5 {
		m.AddRead(other, randomSeq(r, 150))
	}
	for _, s := range regen {
		m.Observe(s)
	}
	res := m.Results()
	if g := res[f]; g.Checked != 5 || g.Exact != 3 || g.Contained != 2 || g.Whole != 1 {
		t.Errorf("suspect: %+v, want 5 checked, 3 exact, 2 contained (the trimmed read whole, the edited one not)", g)
	}
	if g := res[other]; g.Checked != 5 || g.Matched() != 0 {
		t.Errorf("other: %+v, want no match", g)
	}
	if g := res[other]; g.ExpectedChance > 1e-20 || g.ExpectedChance < 0 {
		t.Errorf("expected chance matches %g, want ~0 for random 150 bp reads", g.ExpectedChance)
	}
	m.Reset()
	if g := m.Results()[f]; g.Matched() != 0 || g.Checked != 5 {
		t.Errorf("after Reset: %+v", g)
	}
	if reads, _ := m.Regenerated(); reads != 0 {
		t.Errorf("Regenerated after Reset = %d", reads)
	}
}

func TestCanonicalKmer(t *testing.T) {
	r := rand.New(rand.NewPCG(3, 4))
	s := randomSeq(r, MatchK)
	a, ok1 := canonicalKmer(s)
	b, ok2 := canonicalKmer(revComp(s))
	if !ok1 || !ok2 || a != b {
		t.Errorf("canonical k-mer differs between strands: %x %x", a, b)
	}
	s[5] = 'N'
	if _, ok := canonicalKmer(s); ok {
		t.Error("k-mer with N accepted")
	}
}
