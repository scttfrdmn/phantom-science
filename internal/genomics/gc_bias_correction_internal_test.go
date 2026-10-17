//go:build !release

package genomics

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/population"
)

func newTestBEADS() *BEADSNormalizer {
	return NewBEADSNormalizer()
}

// TestCorrectGCBias_NilCorrection returns data unchanged when no GCBiasCorrection
// is registered for the platform (the common case in production).
func TestCorrectGCBias_NilCorrection(t *testing.T) {
	b := newTestBEADS()
	data := []float64{1.0, 2.0, 3.0}
	got := b.correctGCBias(data, NovaSeq6000, population.CEU)
	for i, v := range got {
		if v != data[i] {
			t.Errorf("no correction: got[%d]=%v, want %v", i, v, data[i])
		}
	}
}

// TestCorrectGCBias_GCContentRange verifies that the GC content computed for each bin
// falls in [0, 100] and that the last bin always maps to exactly 100.
// We do this by injecting a correction whose length is 101 (GC 0–100%) and whose
// factors are f[gc] = 1 + gc/100.  For bin i in a 101-element dataset the expected
// corrected value is data[i] * (1 + i/100).
func TestCorrectGCBias_GCContentRange(t *testing.T) {
	b := newTestBEADS()

	const n = 101
	factors := make([]float64, n) // factors[gc] = 1 + gc/100
	for gc := 0; gc < n; gc++ {
		factors[gc] = 1.0 + float64(gc)/100.0
	}
	b.gcBiasCorrection[NovaSeq6000] = &GCBiasCorrection{
		CorrectionFactors: map[int][]float64{50: factors},
	}

	data := make([]float64, n)
	for i := range data {
		data[i] = 1.0
	}
	got := b.correctGCBias(data, NovaSeq6000, population.CEU)

	// data[0]  → gc=0  → factor=1.0  → 1.0 * 1.0 = 1.0
	// data[50] → gc=50 → factor=1.5  → 1.0 * 1.5 = 1.5
	// data[100]→ gc=100→ factor=2.0  → 1.0 * 2.0 = 2.0
	cases := []struct {
		idx     int
		wantGC  int
		want    float64
	}{
		{0, 0, 1.0},
		{50, 50, 1.5},
		{100, 100, 2.0},
	}
	for _, tc := range cases {
		expected := data[tc.idx] * (1.0 + float64(tc.wantGC)/100.0)
		if got[tc.idx] != expected {
			t.Errorf("bin %d (gc=%d%%): got %v, want %v", tc.idx, tc.wantGC, got[tc.idx], expected)
		}
	}
}

// TestCorrectGCBias_SingleElement maps the only bin to GC=50 (the mid-point default).
func TestCorrectGCBias_SingleElement(t *testing.T) {
	b := newTestBEADS()
	factors := make([]float64, 101)
	for gc := range factors {
		factors[gc] = 1.0
	}
	factors[50] = 3.0 // GC=50% gets a factor of 3
	b.gcBiasCorrection[NovaSeq6000] = &GCBiasCorrection{
		CorrectionFactors: map[int][]float64{50: factors},
	}

	got := b.correctGCBias([]float64{2.0}, NovaSeq6000, population.CEU)
	if got[0] != 6.0 { // 2.0 * 3.0
		t.Errorf("single element: got %v, want 6.0", got[0])
	}
}

// TestCorrectGCBias_PopulationSpecificPreferred verifies that population-specific
// correction factors are preferred over the generic ones when present.
func TestCorrectGCBias_PopulationSpecificPreferred(t *testing.T) {
	b := newTestBEADS()

	// Generic factor for GC=0 is 1.0; population-specific for CEU is 5.0.
	generic := []float64{1.0, 1.0} // index 0 and 1
	specific := []float64{5.0, 1.0}
	b.gcBiasCorrection[NovaSeq6000] = &GCBiasCorrection{
		CorrectionFactors: map[int][]float64{50: generic},
		PopulationSpecific: map[population.Population]map[int][]float64{
			population.CEU: {50: specific},
		},
	}

	// Single data point: maps to GC=50 → index clamped to len-1=1 → factor=1.0 for both.
	// Use two data points so first maps to GC=0 (factor=5.0 for CEU vs 1.0 generic).
	got := b.correctGCBias([]float64{2.0, 2.0}, NovaSeq6000, population.CEU)
	// First bin: gc=0 → population-specific factor=5.0 → 10.0
	if got[0] != 10.0 {
		t.Errorf("population-specific not used: got[0]=%v, want 10.0", got[0])
	}

	// Same data with YRI (no specific) should use generic → factor=1.0 → 2.0
	gotGeneric := b.correctGCBias([]float64{2.0, 2.0}, NovaSeq6000, population.YRI)
	if gotGeneric[0] != 2.0 {
		t.Errorf("generic factors not used for YRI: got[0]=%v, want 2.0", gotGeneric[0])
	}
}

// TestCorrectGCBias_EmptyFactors returns data unchanged when the factor slice is empty.
func TestCorrectGCBias_EmptyFactors(t *testing.T) {
	b := newTestBEADS()
	b.gcBiasCorrection[NovaSeq6000] = &GCBiasCorrection{
		CorrectionFactors: map[int][]float64{}, // no window-size entry
	}
	data := []float64{3.0, 4.0}
	got := b.correctGCBias(data, NovaSeq6000, population.CEU)
	for i, v := range got {
		if v != data[i] {
			t.Errorf("empty factors: got[%d]=%v, want %v", i, v, data[i])
		}
	}
}
