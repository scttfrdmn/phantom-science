//go:build release
// +build release

package fingerprinting

import (
	"log"
)

// For release builds, we only use the legacy implementation to avoid CGO dependencies

// New creates a new fingerprinting engine using legacy implementation
func New(customerID string) Engine {
	log.Printf("Using legacy fingerprinting engine for release build")
	return NewLegacy(customerID)
}
