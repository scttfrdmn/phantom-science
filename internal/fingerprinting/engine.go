//go:build !release
// +build !release

package fingerprinting

import (
	"log"
)

// New creates a fingerprinting engine for customerID, backed by the stegano
// library when built with CGO and falling back to the legacy pure-Go engine
// otherwise.
func New(customerID string) Engine {
	engine, err := newDefaultEngine(customerID)
	if err != nil {
		log.Printf("Stegano library not available, falling back to legacy implementation: %v", err)
		return NewLegacy(customerID)
	}
	return engine
}

// NewVerifier returns a Verifier backed by the stegano library, or an error
// when the build has no stegano support (CGO disabled).
func NewVerifier() (Verifier, error) {
	return newDefaultVerifier()
}
