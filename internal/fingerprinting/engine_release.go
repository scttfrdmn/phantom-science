//go:build release
// +build release

package fingerprinting

// NewVerifier reports that legacy quality-score fingerprints cannot be read
// in release builds, which do not link the stegano library. Keyed-seed
// fingerprints need no library.
func NewVerifier(opts ...Option) (Verifier, error) {
	return nil, NewFingerprintingError("verify", ErrLibraryNotAvailable, ErrorCodeSystem).
		WithContext("reason", "release build")
}
