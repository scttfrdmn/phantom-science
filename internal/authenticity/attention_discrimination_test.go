//go:build !release
// +build !release

package authenticity

import (
	"math"
	"testing"
)

// TestAttentionDiscrimination_EmptyMatrix returns 0.0 when the attention
// matrix has no rows.
func TestAttentionDiscrimination_EmptyMatrix(t *testing.T) {
	got := calculateAttentionDiscrimination(nil, nil)
	if got != 0.0 {
		t.Errorf("empty matrix: expected 0.0, got %v", got)
	}
	got = calculateAttentionDiscrimination([][]float64{}, nil)
	if got != 0.0 {
		t.Errorf("zero-length matrix: expected 0.0, got %v", got)
	}
}

// TestAttentionDiscrimination_AllEmptyRows returns 0.0 when every row is nil
// or empty (the current state of CalculateAttention, which is not yet
// implemented).  This is the honest result: no attention data → no
// discrimination signal.
func TestAttentionDiscrimination_AllEmptyRows(t *testing.T) {
	scores := [][]float64{nil, nil, nil}
	got := calculateAttentionDiscrimination(scores, nil)
	if got != 0.0 {
		t.Errorf("all-nil rows: expected 0.0, got %v", got)
	}
}

// TestAttentionDiscrimination_UniformRows returns 0.0 when every row has
// identical weight on each position (maximum entropy, no selectivity, no
// variation across rows → CV = 0).
func TestAttentionDiscrimination_UniformRows(t *testing.T) {
	n := 4
	w := 1.0 / float64(n)
	row := make([]float64, n)
	for i := range row {
		row[i] = w
	}
	scores := [][]float64{row, row, row}
	got := calculateAttentionDiscrimination(scores, nil)
	// All rows identical → std = 0 → CV = 0 → score = 0.0
	if got != 0.0 {
		t.Errorf("uniform rows: expected 0.0, got %v", got)
	}
}

// TestAttentionDiscrimination_SelectiveRows returns a high score when some
// rows concentrate all weight on a single position (entropy ≈ 0) while
// others spread it uniformly (high entropy), producing high CV.
func TestAttentionDiscrimination_SelectiveRows(t *testing.T) {
	n := 8
	// Uniform row: every position gets weight 1/n
	uniform := make([]float64, n)
	for i := range uniform {
		uniform[i] = 1.0 / float64(n)
	}
	// Spike row: all weight on first position — entropy = 0
	spike := make([]float64, n)
	spike[0] = 1.0

	// Mix of spikes and uniforms → large variance in entropy → high CV
	scores := [][]float64{spike, uniform, spike, uniform, spike}
	got := calculateAttentionDiscrimination(scores, nil)
	if got < 0.5 {
		t.Errorf("mixed spike/uniform rows: expected score >= 0.5, got %v", got)
	}
}

// TestAttentionDiscrimination_OutputRange verifies the result is always in [0, 1].
func TestAttentionDiscrimination_OutputRange(t *testing.T) {
	cases := [][][]float64{
		// Single row with all weight on one element
		{{1.0, 0.0, 0.0, 0.0}},
		// Two rows that are very different
		{{1.0, 0.0}, {0.0, 1.0}},
		// Large uniform matrix
		func() [][]float64 {
			n := 16
			w := 1.0 / float64(n)
			row := make([]float64, n)
			for i := range row {
				row[i] = w
			}
			mat := make([][]float64, 10)
			for i := range mat {
				mat[i] = row
			}
			return mat
		}(),
	}
	for i, c := range cases {
		got := calculateAttentionDiscrimination(c, nil)
		if got < 0.0 || got > 1.0 {
			t.Errorf("case %d: score %v out of [0, 1]", i, got)
		}
		if math.IsNaN(got) || math.IsInf(got, 0) {
			t.Errorf("case %d: non-finite score %v", i, got)
		}
	}
}
