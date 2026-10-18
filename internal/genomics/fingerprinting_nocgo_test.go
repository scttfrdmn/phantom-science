//go:build !release && !cgo
// +build !release,!cgo

package genomics

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Without CGO there is no stegano library, so requesting quality-score
// fingerprints must fail up front instead of writing unmarked files.
func TestGenerate_QualityLayerWithoutStegano(t *testing.T) {
	dir := t.TempDir()
	e := NewEngine(types.SequencerConfig{
		SampleName:    "nocgo",
		OutputDir:     dir,
		ReadLength:    100,
		ReadsPerBatch: 50,
		TotalBatches:  1,
		BaseQuality:   38,
		Fingerprinting: types.FingerprintConfig{
			Enable:              true,
			SteganographyLayers: []string{"quality_scores"},
		},
	})
	defer e.Close()

	err := e.Generate()
	if err == nil || !strings.Contains(err.Error(), "stegano") {
		t.Fatalf("Generate error = %v, want an error naming the stegano requirement", err)
	}
	entries, _ := os.ReadDir(dir)
	if len(entries) != 0 {
		t.Errorf("files written despite the error: %v", entries)
	}

	if err := e.GenerateFASTQFile(filepath.Join(dir, "x.fastq"), 1, 1); err == nil {
		t.Error("GenerateFASTQFile succeeded without stegano")
	}
}
