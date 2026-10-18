//go:build !release
// +build !release

package genomics

import (
	"bytes"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Streamed output cannot carry a quality-score fingerprint, so requesting one
// must fail before any reads are written rather than stream unmarked data.
func TestGenerateToWriter_RejectsQualityLayer(t *testing.T) {
	cfg := types.SequencerConfig{
		SampleName:    "stream",
		ReadLength:    50,
		ReadsPerBatch: 10,
		TotalBatches:  1,
		BaseQuality:   38,
		Fingerprinting: types.FingerprintConfig{
			Enable:              true,
			SteganographyLayers: []string{"read_ids", "quality_scores"},
		},
	}
	e := NewEngine(cfg)
	defer e.Close()

	var buf bytes.Buffer
	err := e.GenerateToWriter(&buf, GenerateOpts{})
	if err == nil || !strings.Contains(err.Error(), "quality_scores") {
		t.Fatalf("GenerateToWriter error = %v, want a quality_scores layer error", err)
	}
	if buf.Len() != 0 {
		t.Errorf("%d bytes streamed before the error", buf.Len())
	}

	// Read-ID tagging still streams.
	cfg.Fingerprinting.SteganographyLayers = []string{"read_ids"}
	e2 := NewEngine(cfg)
	defer e2.Close()
	if err := e2.GenerateToWriter(&buf, GenerateOpts{}); err != nil {
		t.Fatalf("GenerateToWriter with read_ids only: %v", err)
	}
	if got := strings.Count(buf.String(), "\n+\n"); got != 10 {
		t.Errorf("streamed %d records, want 10", got)
	}
}

func TestQualityFingerprintingEnabled(t *testing.T) {
	for _, tc := range []struct {
		enable bool
		layers []string
		want   bool
	}{
		{true, []string{"quality_scores"}, true},
		{true, []string{"read_ids", "quality_scores"}, true},
		{true, []string{"read_ids"}, false},
		{true, nil, false},
		{false, []string{"quality_scores"}, false},
	} {
		e := &Engine{config: types.SequencerConfig{Fingerprinting: types.FingerprintConfig{
			Enable: tc.enable, SteganographyLayers: tc.layers}}}
		if got := e.qualityFingerprintingEnabled(); got != tc.want {
			t.Errorf("enable=%v layers=%v: got %v, want %v", tc.enable, tc.layers, got, tc.want)
		}
	}
}
