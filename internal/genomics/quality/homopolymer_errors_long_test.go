//go:build !release

package quality

import (
	"bytes"
	"strings"
	"testing"
)

// replayErrors applies recorded errors to the original sequence from the last
// to the first, so each error's original Position is still valid.
func replayErrors(t *testing.T, original []byte, errs []HomopolymerError) []byte {
	t.Helper()
	seq := append([]byte(nil), original...)
	for i := len(errs) - 1; i >= 0; i-- {
		e := errs[i]
		switch e.ErrorType {
		case "deletion":
			seq = append(seq[:e.Position], seq[e.Position+1:]...)
		case "insertion":
			seq = append(seq[:e.Position+1], append([]byte{e.Base}, seq[e.Position+1:]...)...)
		case "substitution":
			for _, b := range []byte("ATGC") {
				if b != e.Base {
					seq[e.Position] = b
					break
				}
			}
		default:
			t.Fatalf("unknown error type %q", e.ErrorType)
		}
	}
	return seq
}

// Regression for #21: runs were located on the original sequence but edited in
// the already-modified one, so after an early deletion later edits landed one
// base off and eventually indexed past the end (panic on long inputs).
func TestSimulateErrors_LongHomopolymerRichSequence(t *testing.T) {
	var b strings.Builder
	for i := 0; b.Len() < 6000; i++ {
		b.WriteString(strings.Repeat(string("ACGT"[i%4]), 2+i%9))
		b.WriteString("ACGTTGCA")
	}
	seq := []byte(b.String())
	qual := bytes.Repeat([]byte{30}, len(seq))

	for _, platform := range []string{"nanopore", "pacbio", "illumina"} {
		for seed := int64(1); seed <= 40; seed++ {
			m := NewHomopolymerErrorModel(platform)
			m.SetSeed(seed)
			res, err := m.SimulateErrors(seq, qual)
			if err != nil {
				t.Fatalf("%s seed %d: %v", platform, seed, err)
			}
			if len(res.ErrorQualities) != len(res.ErrorSequence) {
				t.Fatalf("%s seed %d: %d qualities for %d bases", platform, seed, len(res.ErrorQualities), len(res.ErrorSequence))
			}
			if want := replayErrors(t, seq, res.Errors); !bytes.Equal(res.ErrorSequence, want) {
				t.Fatalf("%s seed %d: ErrorSequence differs from replaying its %d recorded errors", platform, seed, len(res.Errors))
			}
		}
	}
}
