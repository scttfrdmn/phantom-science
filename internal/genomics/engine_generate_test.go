//go:build !release
// +build !release

package genomics

import (
	"bytes"
	"math/rand/v2"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// testRand returns a fixed-seed generator for unit tests.
func testRand() *rand.Rand { return newStream(1, 0) }

func minimalConfig(t *testing.T) types.SequencerConfig {
	t.Helper()
	return types.SequencerConfig{
		SampleName:    "TestSample",
		Lane:          1,
		OutputDir:     t.TempDir(),
		ReadLength:    50,
		ReadsPerBatch: 10,
		TotalBatches:  1,
		BaseQuality:   30,
		QualityDecay:  0.1,
	}
}

// TestGenerateToWriter_Basic verifies that FASTQ output is written to the
// provided writer in the expected format (4 lines per read).
func TestGenerateToWriter_Basic(t *testing.T) {
	cfg := minimalConfig(t)
	engine := NewEngine(cfg)

	var buf bytes.Buffer
	if err := engine.GenerateToWriter(&buf, GenerateOpts{}); err != nil {
		t.Fatalf("GenerateToWriter: %v", err)
	}

	output := buf.String()
	lines := strings.Split(strings.TrimRight(output, "\n"), "\n")
	// 10 reads × 4 lines each = 40 lines.
	expectedLines := cfg.ReadsPerBatch * cfg.TotalBatches * 4
	if len(lines) != expectedLines {
		t.Errorf("expected %d lines, got %d", expectedLines, len(lines))
	}
	// Every 4th line starting at 0 should begin with '@'.
	for i := 0; i < len(lines); i += 4 {
		if !strings.HasPrefix(lines[i], "@") {
			t.Errorf("line %d should be a FASTQ header (start with @), got: %q", i, lines[i])
		}
		if lines[i+2] != "+" {
			t.Errorf("line %d should be '+', got: %q", i+2, lines[i+2])
		}
	}
}

// TestGenerateToWriter_SampleID verifies that the RG:Z: tag appears in every header.
func TestGenerateToWriter_SampleID(t *testing.T) {
	cfg := minimalConfig(t)
	engine := NewEngine(cfg)

	var buf bytes.Buffer
	opts := GenerateOpts{SampleID: "tumor_42"}
	if err := engine.GenerateToWriter(&buf, opts); err != nil {
		t.Fatalf("GenerateToWriter: %v", err)
	}

	output := buf.String()
	lines := strings.Split(strings.TrimRight(output, "\n"), "\n")
	for i := 0; i < len(lines); i += 4 {
		header := lines[i]
		if !strings.Contains(header, "RG:Z:tumor_42") {
			t.Errorf("header line %d missing RG tag: %q", i, header)
		}
	}
}

// TestGenerateToWriter_PairedEnd verifies paired-end output doubles the read count.
func TestGenerateToWriter_PairedEnd(t *testing.T) {
	cfg := minimalConfig(t)
	cfg.PairedEnd = true
	engine := NewEngine(cfg)

	var buf bytes.Buffer
	if err := engine.GenerateToWriter(&buf, GenerateOpts{}); err != nil {
		t.Fatalf("GenerateToWriter: %v", err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	expected := cfg.ReadsPerBatch * cfg.TotalBatches * 2 * 4 // ×2 for paired
	if len(lines) != expected {
		t.Errorf("expected %d lines for paired-end, got %d", expected, len(lines))
	}
}

// TestGenerateToWriter_Variants verifies that variant injection runs without error
// and does not corrupt read length.
func TestGenerateToWriter_Variants(t *testing.T) {
	cfg := minimalConfig(t)
	engine := NewEngine(cfg)

	variants := []VariantRecord{
		{Ref: "ACGT", Alt: "TTTT", Freq: 1.0}, // always inject
	}
	var buf bytes.Buffer
	opts := GenerateOpts{Variants: variants}
	if err := engine.GenerateToWriter(&buf, opts); err != nil {
		t.Fatalf("GenerateToWriter with variants: %v", err)
	}

	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	for i := 1; i < len(lines); i += 4 {
		seq := lines[i]
		if len(seq) != cfg.ReadLength {
			t.Errorf("read at line %d has length %d, want %d", i, len(seq), cfg.ReadLength)
		}
	}
}

// TestInjectVariants_AlwaysInject verifies that Freq=1.0 always substitutes when REF is found.
func TestInjectVariants_AlwaysInject(t *testing.T) {
	seq := []byte("AAACGTAAA")
	variants := []VariantRecord{{Ref: "ACGT", Alt: "TTTT", Freq: 1.0}}
	result, _ := injectVariants(testRand(), seq, nil, variants)
	if len(result) != len(seq) {
		t.Errorf("result length %d != original %d", len(result), len(seq))
	}
	if !bytes.Contains(result, []byte("TTTT")) {
		t.Errorf("expected TTTT in result, got %q", result)
	}
}

// TestInjectVariants_NeverInject verifies that Freq=0.0 never substitutes.
func TestInjectVariants_NeverInject(t *testing.T) {
	seq := []byte("AAACGTAAA")
	original := make([]byte, len(seq))
	copy(original, seq)
	variants := []VariantRecord{{Ref: "ACGT", Alt: "TTTT", Freq: 0.0}}
	result, _ := injectVariants(testRand(), seq, nil, variants)
	if !bytes.Equal(result, original) {
		t.Errorf("Freq=0 should not modify sequence; got %q, want %q", result, original)
	}
}

// TestInjectVariants_MissingRef verifies that when REF is not found the sequence is unchanged.
func TestInjectVariants_MissingRef(t *testing.T) {
	seq := []byte("AAAAAAAA")
	original := make([]byte, len(seq))
	copy(original, seq)
	variants := []VariantRecord{{Ref: "CCCC", Alt: "TTTT", Freq: 1.0}}
	result, _ := injectVariants(testRand(), seq, nil, variants)
	if !bytes.Equal(result, original) {
		t.Errorf("absent REF should leave sequence unchanged; got %q", result)
	}
}

// TestInjectVariants_IndelChangesLength verifies that injection splices the
// allele without truncating or padding: the read length changes by exactly
// the indel size (the pipeline trims fixed-length reads afterwards), and
// qualities stay aligned with the bases.
func TestInjectVariants_IndelChangesLength(t *testing.T) {
	seq := []byte("AAATGAAA") // len=8
	qual := []byte{1, 2, 3, 4, 5, 6, 7, 8}

	// Longer ALT: replace "ATG" (3) with "GGGGGG" (6) → len 11.
	res, q := injectVariants(testRand(), seq, qual, []VariantRecord{{Ref: "ATG", Alt: "GGGGGG", Freq: 1.0}})
	if string(res) != "AAGGGGGGAAA" {
		t.Errorf("longer ALT: got %q, want AAGGGGGGAAA", res)
	}
	if !bytes.Equal(q, []byte{1, 2, 3, 4, 5, 5, 5, 5, 6, 7, 8}) {
		t.Errorf("longer ALT qualities: got %v", q)
	}

	// Shorter ALT: replace "ATG" (3) with "G" (1) → len 6, no N padding.
	res2, q2 := injectVariants(testRand(), []byte("AAATGAAA"), qual, []VariantRecord{{Ref: "ATG", Alt: "G", Freq: 1.0}})
	if string(res2) != "AAGAAA" {
		t.Errorf("shorter ALT: got %q, want AAGAAA", res2)
	}
	if !bytes.Equal(q2, []byte{1, 2, 3, 6, 7, 8}) {
		t.Errorf("shorter ALT qualities: got %v", q2)
	}
	if string(seq) != "AAATGAAA" {
		t.Errorf("input modified: %q", seq)
	}
}
