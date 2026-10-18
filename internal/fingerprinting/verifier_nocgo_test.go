//go:build !release && !cgo
// +build !release,!cgo

package fingerprinting

import (
	"errors"
	"testing"
)

// Without the stegano library NewVerifier must refuse, not return a verifier
// that reports success.
func TestNewVerifier_UnavailableWithoutCGO(t *testing.T) {
	v, err := NewVerifier()
	if err == nil || v != nil {
		t.Fatalf("NewVerifier() = %v, %v; want nil verifier and an error", v, err)
	}
	if !errors.Is(err, ErrLibraryNotAvailable) {
		t.Errorf("error %v does not wrap ErrLibraryNotAvailable", err)
	}
}
