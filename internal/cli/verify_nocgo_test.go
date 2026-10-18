//go:build !release && !cgo
// +build !release,!cgo

package cli

import (
	"bytes"
	"strings"
	"testing"
)

// Without the stegano library verification must fail, not report success.
func TestVerifyFiles_NeedsStegano(t *testing.T) {
	files := generateFASTQ(t, testSequencerConfig(t.TempDir()))
	err := VerifyFiles(&bytes.Buffer{}, files, VerifyOptions{})
	if err == nil || !strings.Contains(err.Error(), "stegano") {
		t.Fatalf("VerifyFiles error = %v, want a stegano requirement error", err)
	}
}
