package cli

import (
	"errors"
	"fmt"
	"io"
	"os"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
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
		result, err := verifier.Verify(path)
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
func ValidateFiles(w io.Writer, paths []string, cfg types.ValidationConfig) error {
	if len(paths) == 0 {
		return fmt.Errorf("no files to validate (usage: phantom-science validate FILE...)")
	}
	framework := validation.NewFramework(cfg)

	failed := 0
	for _, path := range paths {
		data, err := os.ReadFile(path) // #nosec G304 -- user-supplied input files are the point of this command
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
