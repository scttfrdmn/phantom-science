//go:build !cgo
// +build !cgo

package fingerprinting

// newDefaultVerifier reports that legacy verification is unavailable:
// extracting a quality-score fingerprint needs the stegano library, which
// needs CGO.
func newDefaultVerifier(opts []Option) (Verifier, error) {
	return nil, NewFingerprintingError("verify", ErrLibraryNotAvailable, ErrorCodeSystem).
		WithContext("reason", "built without CGO")
}
