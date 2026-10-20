//go:build !release
// +build !release

package fingerprinting

// NewVerifier returns a Verifier for legacy quality-score fingerprints,
// backed by the stegano library, or an error when the build has no stegano
// support (CGO disabled).
func NewVerifier(opts ...Option) (Verifier, error) {
	return newDefaultVerifier(opts)
}
