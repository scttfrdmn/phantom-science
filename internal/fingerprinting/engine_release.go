//go:build release
// +build release

package fingerprinting

import (
	"log"
)

// For release builds, we only use the legacy implementation to avoid CGO dependencies

// New creates a new fingerprinting engine using legacy implementation
func New(customerID string, opts ...Option) Engine {
	log.Printf("Using legacy fingerprinting engine for release build")
	return NewLegacy(customerID)
}

// NewVerifier reports that verification is unavailable in release builds,
// which do not link the stegano library.
func NewVerifier(opts ...Option) (Verifier, error) {
	return nil, NewFingerprintingError("verify", ErrLibraryNotAvailable, ErrorCodeSystem).
		WithContext("reason", "release build")
}
