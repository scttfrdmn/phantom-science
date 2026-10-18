//go:build !cgo
// +build !cgo

package fingerprinting

import "fmt"

// newDefaultEngine reports that the stegano engine is unavailable: it needs
// CGO to link the stegano library.
func newDefaultEngine(customerID string) (Engine, error) {
	return nil, fmt.Errorf("stegano engine not available (CGO disabled)")
}

// newDefaultVerifier reports that verification is unavailable: extracting a
// fingerprint needs the stegano library, which needs CGO.
func newDefaultVerifier() (Verifier, error) {
	return nil, NewFingerprintingError("verify", ErrLibraryNotAvailable, ErrorCodeSystem).
		WithContext("reason", "built without CGO")
}
