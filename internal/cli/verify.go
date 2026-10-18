package cli

import (
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics/readlen"
	"github.com/scttfrdmn/phantom-science/internal/validation"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// VerifyOptions configures VerifyFiles. Customer and Dataset, when set, are
// required values. KeyFile holds the key for encrypted fingerprints; when empty
// the fingerprinting.KeyEnvVar environment variable is used.
type VerifyOptions struct {
	Customer string
	Dataset  string
	KeyFile  string
}

// VerifyFiles checks each FASTQ file in paths for a steganographic
// fingerprint and reports one line per file to w. It returns an error if any
// file has no fingerprint, a fingerprint that does not match opts, or cannot
// be checked (including an encrypted fingerprint without the right key).
func VerifyFiles(w io.Writer, paths []string, opts VerifyOptions) error {
	if len(paths) == 0 {
		return fmt.Errorf("no files to verify (usage: phantom-science verify FILE...)")
	}
	key, err := fingerprinting.ResolveKey(opts.KeyFile)
	if err != nil {
		return err
	}
	var vopts []fingerprinting.Option
	if key != nil {
		vopts = append(vopts, fingerprinting.WithKey(key))
	}
	verifier, err := fingerprinting.NewVerifier(vopts...)
	if err != nil {
		return fmt.Errorf("fingerprint verification needs the stegano library (build with CGO after `make stegano`): %w", err)
	}

	failed := 0
	for _, path := range paths {
		result, err := verifyPossiblyGzipped(verifier, path)
		switch {
		case errors.Is(err, fingerprinting.ErrKeyRequired):
			fmt.Fprintf(w, "❌ %s: fingerprint is encrypted; pass --key-file or set %s\n", path, fingerprinting.KeyEnvVar)
			failed++
		case errors.Is(err, fingerprinting.ErrWrongKey):
			fmt.Fprintf(w, "❌ %s: fingerprint cannot be decrypted with this key\n", path)
			failed++
		case err != nil:
			fmt.Fprintf(w, "❌ %s: cannot verify: %v\n", path, err)
			failed++
		case !result.Found:
			fmt.Fprintf(w, "❌ %s: no fingerprint\n", path)
			failed++
		default:
			if msg := checkFingerprint(result.Fingerprint, opts); msg != "" {
				fmt.Fprintf(w, "❌ %s: %s\n", path, msg)
				failed++
			} else {
				protection := "plaintext"
				if result.Encrypted {
					protection = "encrypted"
				}
				fmt.Fprintf(w, "✅ %s: %s (%s)\n", path, describeFingerprint(result.Fingerprint), protection)
			}
		}
	}
	if failed > 0 {
		return fmt.Errorf("%d of %d files failed verification", failed, len(paths))
	}
	return nil
}

// checkFingerprint returns a reason payload fails expect, or "" if it passes.
func checkFingerprint(payload string, expect VerifyOptions) string {
	if expect.Customer == "" && expect.Dataset == "" {
		return ""
	}
	fp, err := fingerprinting.ParseFingerprint(payload)
	if err != nil {
		return fmt.Sprintf("fingerprint %q cannot be checked against expectations: %v", payload, err)
	}
	if expect.Customer != "" && fp.Customer != expect.Customer {
		return fmt.Sprintf("fingerprint customer %q, expected %q", fp.Customer, expect.Customer)
	}
	if expect.Dataset != "" && fp.Dataset != expect.Dataset {
		return fmt.Sprintf("fingerprint dataset %q, expected %q", fp.Dataset, expect.Dataset)
	}
	return ""
}

// describeFingerprint renders payload's fields, or the raw payload if it is
// not in the standard format.
func describeFingerprint(payload string) string {
	fp, err := fingerprinting.ParseFingerprint(payload)
	if err != nil {
		return fmt.Sprintf("fingerprint %q", payload)
	}
	return fmt.Sprintf("customer=%s dataset=%s created=%s stegano=%s",
		fp.Customer, fp.Dataset, fp.Created.Format(time.RFC3339), fp.SteganoVersion)
}

// DefaultValidationTests is the validation suite used when no phantom config
// enables validation: every built-in check that works on a single FASTQ file.
var DefaultValidationTests = types.ValidationConfig{
	Enable:                     true,
	StatisticalTests:           []string{"format_compliance", "quality_distribution", "sequence_complexity", "read_length_consistency"},
	PopulationGeneticsTests:    []string{"base_composition", "gc_content"},
	PipelineCompatibilityTests: []string{"fastq_parsing", "quality_encoding", "header_format"},
}

// EffectiveValidation returns the validation tests to run for a phantom
// config's validation section: the section itself when it lists tests, and
// DefaultValidationTests when it is enabled without listing any. ok is false
// when validation is not enabled.
func EffectiveValidation(cfg types.ValidationConfig) (effective types.ValidationConfig, ok bool) {
	if !cfg.Enable {
		return types.ValidationConfig{}, false
	}
	if len(cfg.StatisticalTests)+len(cfg.PopulationGeneticsTests)+len(cfg.PipelineCompatibilityTests) == 0 {
		return DefaultValidationTests, true
	}
	return cfg, true
}

// ValidateFiles runs the validation framework configured by cfg on each FASTQ
// file in paths and reports the results to w. It returns an error if any
// file cannot be read or fails validation.
//
// Without options the read_length_consistency test accepts variable read
// lengths and reports them informationally (validation.ReadLengthsAny), so
// long-read files whose expected distribution is unknown do not fail just
// because lengths vary. Pass ReadLengthValidation(config) after generating
// with config to check lengths against the configured distribution instead.
func ValidateFiles(w io.Writer, paths []string, cfg types.ValidationConfig, opts ...validation.Option) error {
	if len(paths) == 0 {
		return fmt.Errorf("no files to validate (usage: phantom-science validate FILE...)")
	}
	opts = append([]validation.Option{validation.WithReadLengths(validation.ReadLengthExpectation{Mode: validation.ReadLengthsAny})}, opts...)
	framework := validation.NewFramework(cfg, opts...)

	failed := 0
	for _, path := range paths {
		data, err := readFASTQ(path)
		if err != nil {
			fmt.Fprintf(w, "❌ %s: %v\n", path, err)
			failed++
			continue
		}
		results := framework.RunValidation(string(data))
		status := "✅"
		if !results.OverallPass {
			status = "❌"
			failed++
		}
		fmt.Fprintf(w, "%s %s: %s\n", status, path, results.Summary)
		for _, test := range results.Tests {
			mark := "✅"
			if !test.Passed {
				mark = "❌"
			}
			fmt.Fprintf(w, "    %s %s: %s\n", mark, test.TestName, test.Message)
		}
	}
	if failed > 0 {
		return fmt.Errorf("%d of %d files failed validation", failed, len(paths))
	}
	return nil
}

// ReadLengthValidation returns the read_length_consistency expectation for
// files generated with config, for use as a ValidateFiles option:
//
//   - fixed distribution: every read exactly config.ReadLength bases;
//   - lognormal, gamma or empirical: lengths within the distribution's
//     [Min, Max] (Min defaults to 1) and a sample mean within 4·SD/√n of the
//     distribution's mean (after truncation);
//   - with a sequencing error model (ErrorModel.Tier set), whose indels change
//     read lengths, lengths are only reported (validation.ReadLengthsAny).
//
// It returns an error if config's distribution is invalid.
func ReadLengthValidation(config types.SequencerConfig) (validation.Option, error) {
	d, err := readlen.New(config.ReadLengthDistribution, config.ReadLength)
	if err != nil {
		return nil, fmt.Errorf("read-length distribution: %w", err)
	}
	exp := validation.ReadLengthExpectation{Name: d.Dist()}
	switch {
	case config.ErrorModel.Tier != "":
		exp.Mode = validation.ReadLengthsAny
	case d.Fixed():
		exp.Mode, exp.Length = validation.ReadLengthsUniform, config.ReadLength
	default:
		exp.Mode = validation.ReadLengthsDistributed
		exp.Min, exp.Max = d.Bounds()
		exp.Mean, exp.SD = d.Mean(), d.SD()
	}
	return validation.WithReadLengths(exp), nil
}

// isGzip reports whether path names gzip-compressed output (.gz).
func isGzip(path string) bool {
	return strings.HasSuffix(path, ".gz")
}

// readFASTQ returns the contents of a FASTQ file, decompressing .gz files.
func readFASTQ(path string) ([]byte, error) {
	f, err := os.Open(path) // #nosec G304 -- user-supplied input files are the point of this command
	if err != nil {
		return nil, err
	}
	defer f.Close()
	if !isGzip(path) {
		return io.ReadAll(f)
	}
	zr, err := gzip.NewReader(f)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	defer zr.Close()
	return io.ReadAll(zr)
}

// verifyPossiblyGzipped verifies path, first decompressing a .gz file to a
// private temporary file: fingerprints are embedded in, and extracted from,
// the plain FASTQ.
func verifyPossiblyGzipped(v fingerprinting.Verifier, path string) (fingerprinting.VerifyResult, error) {
	if !isGzip(path) {
		return v.Verify(path)
	}
	in, err := os.Open(path) // #nosec G304 -- user-supplied input files are the point of this command
	if err != nil {
		return fingerprinting.VerifyResult{}, err
	}
	defer in.Close()
	zr, err := gzip.NewReader(in)
	if err != nil {
		return fingerprinting.VerifyResult{}, fmt.Errorf("%s: %w", path, err)
	}
	defer zr.Close()
	tmp, err := os.CreateTemp("", "phantom-verify-*.fastq") // created with mode 0600
	if err != nil {
		return fingerprinting.VerifyResult{}, err
	}
	defer os.Remove(tmp.Name())
	if _, err := io.Copy(tmp, zr); err != nil {
		tmp.Close()
		return fingerprinting.VerifyResult{}, fmt.Errorf("decompressing %s: %w", path, err)
	}
	if err := tmp.Close(); err != nil {
		return fingerprinting.VerifyResult{}, err
	}
	return v.Verify(tmp.Name())
}

