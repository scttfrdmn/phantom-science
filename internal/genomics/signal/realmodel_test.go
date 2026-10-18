//go:build !release

package signal

import (
	"os"
	"testing"
)

// TestRealModelSmoke loads ONT tables downloaded with
// scripts/fetch-pore-model.sh, when pointed at them:
//
//	PHANTOM_PORE_MODEL_DNA=~/.cache/phantom-science/pore-models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt
//	PHANTOM_PORE_MODEL_RNA=~/.cache/phantom-science/pore-models/rna004/9mer_levels_v1.txt
//
// It is skipped otherwise (models are not vendored and tests never use the
// network).
func TestRealModelSmoke(t *testing.T) {
	for _, c := range []struct{ env, preset string }{
		{"PHANTOM_PORE_MODEL_DNA", PresetDNAR10},
		{"PHANTOM_PORE_MODEL_RNA", PresetRNA004},
	} {
		path := os.Getenv(c.env)
		if path == "" {
			t.Logf("%s not set; skipping %s", c.env, c.preset)
			continue
		}
		p, _ := Preset(c.preset)
		m, err := LoadPoreModel(path, p.ModelOptions()...)
		if err != nil {
			t.Fatal(err)
		}
		if m.K != 9 || m.Missing != 0 {
			t.Fatalf("%s: K=%d missing=%d", path, m.K, m.Missing)
		}
		lo, hi := m.Mean[0], m.Mean[0]
		for _, v := range m.Mean {
			lo, hi = min(lo, v), max(hi, v)
		}
		// Levels must sit well inside the ADC window and below open pore.
		if lo < 10 || hi > float32(p.MedianBeforeMean) {
			t.Errorf("%s: levels %.1f–%.1f pA look mis-scaled", c.preset, lo, hi)
		}
		r := testRand(1)
		rd, err := Synthesize(r, randomSeq(r, 5000), p, m, WithPolyA(40))
		if err != nil {
			t.Fatal(err)
		}
		if rd.Clipped != 0 {
			t.Errorf("%s: %d clipped samples", c.preset, rd.Clipped)
		}
		t.Logf("%s: levels %.1f–%.1f pA; 5 kb read → %d samples (%.2f s)", c.preset, lo, hi, rd.NumSamples, rd.Duration)
	}
}
