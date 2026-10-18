//go:build !release

package cli

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/validation"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func lognormalSequencerConfig(dir string) types.SequencerConfig {
	c := testSequencerConfig(dir)
	c.Seed = 1
	c.ReadsPerBatch = 300
	// The linear-decay quality model drops to Q2 after ~150 bases, which
	// fails quality_distribution on long reads; the long-read quality model
	// is #11. Decay slowly so these tests isolate read lengths.
	c.QualityDecay = 0.002
	c.ReadLengthDistribution = types.ReadLengthSpec{Dist: types.ReadLengthLognormal, N50: 6000, Median: 3000, Min: 500, Max: 50000}
	return c
}

func mustReadLengthValidation(t *testing.T, c types.SequencerConfig) validation.Option {
	t.Helper()
	opt, err := ReadLengthValidation(c)
	if err != nil {
		t.Fatal(err)
	}
	return opt
}

// Variable-length output passes the default suite (lengths reported, not
// judged) and the distribution check derived from its config, but fails the
// strict fixed-length check.
func TestValidateFiles_VariableLengths(t *testing.T) {
	c := lognormalSequencerConfig(t.TempDir())
	files := generateFASTQ(t, c)

	var out bytes.Buffer
	if err := ValidateFiles(&out, files, DefaultValidationTests); err != nil {
		t.Fatalf("default suite rejected variable-length reads: %v\n%s", err, out.String())
	}
	if !strings.Contains(out.String(), "Variable read lengths") {
		t.Errorf("read lengths not reported:\n%s", out.String())
	}

	out.Reset()
	if err := ValidateFiles(&out, files, DefaultValidationTests, mustReadLengthValidation(t, c)); err != nil {
		t.Fatalf("distribution check: %v\n%s", err, out.String())
	}
	if !strings.Contains(out.String(), "match the lognormal distribution") {
		t.Errorf("distribution check not reported:\n%s", out.String())
	}

	out.Reset()
	strict := validation.WithReadLengths(validation.ReadLengthExpectation{Mode: validation.ReadLengthsUniform})
	if err := ValidateFiles(&out, files, DefaultValidationTests, strict); err == nil {
		t.Errorf("strict fixed-length check accepted variable lengths:\n%s", out.String())
	}

	// The expectation of a fixed-length config rejects the same files.
	out.Reset()
	if err := ValidateFiles(&out, files, DefaultValidationTests, mustReadLengthValidation(t, testSequencerConfig(""))); err == nil {
		t.Errorf("fixed-length expectation accepted variable lengths:\n%s", out.String())
	}
}

func TestValidateFiles_FixedLengthsStrict(t *testing.T) {
	c := testSequencerConfig(t.TempDir())
	files := generateFASTQ(t, c)
	var out bytes.Buffer
	if err := ValidateFiles(&out, files, DefaultValidationTests, mustReadLengthValidation(t, c)); err != nil {
		t.Fatalf("fixed-length expectation: %v\n%s", err, out.String())
	}
	c.ReadLength = 100
	out.Reset()
	if err := ValidateFiles(&out, files, DefaultValidationTests, mustReadLengthValidation(t, c)); err == nil {
		t.Errorf("150-base reads accepted for read length 100:\n%s", out.String())
	}
}

func TestValidateFiles_OutOfBoundsLengthFails(t *testing.T) {
	c := lognormalSequencerConfig(t.TempDir())
	files := generateFASTQ(t, c)
	data, err := os.ReadFile(files[0])
	if err != nil {
		t.Fatal(err)
	}
	long := strings.Repeat("A", 50001)
	data = append(data, []byte("@extra\n"+long+"\n+\n"+strings.Repeat("I", len(long))+"\n")...)
	path := filepath.Join(t.TempDir(), "extra.fastq")
	if err := os.WriteFile(path, data, 0o600); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if err := ValidateFiles(&out, []string{path}, types.ValidationConfig{Enable: true, StatisticalTests: []string{"read_length_consistency"}},
		mustReadLengthValidation(t, c)); err == nil || !strings.Contains(out.String(), "longer than max 50000") {
		t.Errorf("read above max accepted: err %v\n%s", err, out.String())
	}
}

func TestReadLengthValidation_Modes(t *testing.T) {
	if _, err := ReadLengthValidation(types.SequencerConfig{ReadLength: 150, ReadLengthDistribution: types.ReadLengthSpec{Dist: "lognormal"}}); err == nil {
		t.Error("invalid distribution accepted")
	}
	// An error model changes lengths: lengths are only reported.
	c := lognormalSequencerConfig("")
	c.ErrorModel.Tier = "sup"
	f := validation.NewFramework(types.ValidationConfig{Enable: true, StatisticalTests: []string{"read_length_consistency"}}, mustReadLengthValidation(t, c))
	if res := f.RunValidation("@a\nACGT\n+\nIIII\n@b\nACGTACGT\n+\nIIIIIIII\n"); !res.OverallPass {
		t.Errorf("error-model expectation rejected variable lengths: %+v", res.Tests)
	}
}
