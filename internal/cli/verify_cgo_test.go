//go:build !release && cgo
// +build !release,cgo

package cli

import (
	"bytes"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestVerifyFiles(t *testing.T) {
	dir := t.TempDir()
	cfg := testSequencerConfig(filepath.Join(dir, "fp"))
	cfg.Fingerprinting = types.FingerprintConfig{Enable: true, SteganographyLayers: []string{"quality_scores"}}
	marked := generateFASTQ(t, cfg)[0]
	plain := generateFASTQ(t, testSequencerConfig(filepath.Join(dir, "plain")))[0]

	var out bytes.Buffer
	if err := VerifyFiles(&out, []string{marked}, VerifyExpectation{Customer: "CUST_1", Dataset: "ds-1"}); err != nil {
		t.Fatalf("VerifyFiles: %v\n%s", err, out.String())
	}
	if !strings.Contains(out.String(), "✅") || !strings.Contains(out.String(), "customer=CUST_1 dataset=ds-1") {
		t.Errorf("unexpected report:\n%s", out.String())
	}

	out.Reset()
	if err := VerifyFiles(&out, []string{marked}, VerifyExpectation{Customer: "OTHER"}); err == nil {
		t.Errorf("customer mismatch passed:\n%s", out.String())
	} else if !strings.Contains(out.String(), `expected "OTHER"`) {
		t.Errorf("mismatch not reported:\n%s", out.String())
	}

	out.Reset()
	err := VerifyFiles(&out, []string{marked, plain, filepath.Join(dir, "missing.fastq")}, VerifyExpectation{})
	if err == nil || !strings.Contains(err.Error(), "2 of 3") {
		t.Fatalf("err = %v, want 2 of 3 failures\n%s", err, out.String())
	}
	if !strings.Contains(out.String(), "no fingerprint") || !strings.Contains(out.String(), "cannot verify") {
		t.Errorf("unmarked or missing file not reported:\n%s", out.String())
	}
}
