package fingerprinting

import (
	"bytes"
	"fmt"
	"os"
	"runtime"
)

// KeyEnvVar names the environment variable that supplies the fingerprint
// encryption key when no key file is configured.
const KeyEnvVar = "PHANTOM_FINGERPRINT_KEY"

// MinKeyLength is the minimum key length in bytes.
const MinKeyLength = 16

// ValidateKey checks that key is usable as a fingerprint encryption key. The
// key is passed to stegano as a C string, so it must not contain NUL bytes;
// text keys such as `openssl rand -base64 32` output are expected.
func ValidateKey(key []byte) error {
	if len(key) < MinKeyLength {
		return fmt.Errorf("fingerprint key is %d bytes; at least %d are required", len(key), MinKeyLength)
	}
	if bytes.IndexByte(key, 0) >= 0 {
		return fmt.Errorf("fingerprint key contains a NUL byte; use a text key (for example `openssl rand -base64 32`)")
	}
	return nil
}

// LoadKeyFile reads a fingerprint key from path, dropping one trailing
// newline. On Unix the file must not be accessible by group or others.
func LoadKeyFile(path string) ([]byte, error) {
	info, err := os.Stat(path)
	if err != nil {
		return nil, fmt.Errorf("fingerprint key file: %w", err)
	}
	if runtime.GOOS != "windows" && info.Mode().Perm()&0o077 != 0 {
		return nil, fmt.Errorf("fingerprint key file %s is accessible by other users (mode %04o); run `chmod 600 %s`",
			path, info.Mode().Perm(), path)
	}
	key, err := os.ReadFile(path) // #nosec G304 -- the key file path is operator configuration
	if err != nil {
		return nil, fmt.Errorf("fingerprint key file: %w", err)
	}
	key = bytes.TrimSuffix(key, []byte("\n"))
	key = bytes.TrimSuffix(key, []byte("\r"))
	if err := ValidateKey(key); err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return key, nil
}

// ResolveKey returns the fingerprint key from keyFile when it is set, and
// otherwise from the KeyEnvVar environment variable. It returns nil and no
// error when neither is set.
func ResolveKey(keyFile string) ([]byte, error) {
	if keyFile != "" {
		return LoadKeyFile(keyFile)
	}
	if env := os.Getenv(KeyEnvVar); env != "" {
		key := []byte(env)
		if err := ValidateKey(key); err != nil {
			return nil, fmt.Errorf("%s: %w", KeyEnvVar, err)
		}
		return key, nil
	}
	return nil, nil
}

// Option configures an Engine or Verifier.
type Option func(*options)

type options struct {
	key []byte
}

// WithKey encrypts embedded fingerprints with key and decrypts them on
// extraction. key should pass ValidateKey. Builds without the stegano library
// ignore it.
func WithKey(key []byte) Option {
	return func(o *options) { o.key = key }
}

func applyOptions(opts []Option) options {
	var o options
	for _, opt := range opts {
		opt(&o)
	}
	return o
}
