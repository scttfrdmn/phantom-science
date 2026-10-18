//go:build release
// +build release

package fingerprinting

import (
	"log"
)

// Release builds do not link the stegano library, so they only have the
// legacy engine.

// New returns a LegacyEngine, which tags read IDs but cannot embed
// quality-score fingerprints.
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
