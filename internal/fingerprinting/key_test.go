//go:build !release
// +build !release

package fingerprinting

import (
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"
)

func writeKey(t *testing.T, content string, mode os.FileMode) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "fingerprint.key")
	if err := os.WriteFile(path, []byte(content), mode); err != nil {
		t.Fatal(err)
	}
	if err := os.Chmod(path, mode); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestLoadKeyFile(t *testing.T) {
	const key = "0123456789abcdef-key" // gitleaks:allow -- test fixture, not a real key

	got, err := LoadKeyFile(writeKey(t, key+"\n", 0o600))
	if err != nil || string(got) != key {
		t.Errorf("LoadKeyFile = %q, %v; want %q (trailing newline dropped)", got, err, key)
	}
	if got, err := LoadKeyFile(writeKey(t, key+"\r\n", 0o600)); err != nil || string(got) != key {
		t.Errorf("CRLF key = %q, %v; want %q", got, err, key)
	}

	for name, tc := range map[string]struct {
		content string
		mode    os.FileMode
		want    string
	}{
		"too short":  {"short\n", 0o600, "at least"},
		"NUL byte":   {"0123456789abc\x00def", 0o600, "NUL"},
		"group read": {key, 0o640, "chmod 600"},
		"world read": {key, 0o604, "chmod 600"},
	} {
		t.Run(name, func(t *testing.T) {
			if runtime.GOOS == "windows" && strings.Contains(name, "read") {
				t.Skip("permission bits are not checked on Windows")
			}
			if _, err := LoadKeyFile(writeKey(t, tc.content, tc.mode)); err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Errorf("LoadKeyFile error = %v, want one mentioning %q", err, tc.want)
			}
		})
	}

	if _, err := LoadKeyFile(filepath.Join(t.TempDir(), "missing.key")); err == nil {
		t.Error("missing key file accepted")
	}
}

func TestResolveKey(t *testing.T) {
	const fileKey, envKey = "file-key-0123456789", "env-key-0123456789"
	path := writeKey(t, fileKey, 0o600)

	t.Setenv(KeyEnvVar, envKey)
	if got, err := ResolveKey(path); err != nil || string(got) != fileKey {
		t.Errorf("ResolveKey(file) = %q, %v; the key file takes precedence", got, err)
	}
	if got, err := ResolveKey(""); err != nil || string(got) != envKey {
		t.Errorf("ResolveKey(\"\") = %q, %v; want the %s value", got, err, KeyEnvVar)
	}

	t.Setenv(KeyEnvVar, "short")
	if _, err := ResolveKey(""); err == nil || !strings.Contains(err.Error(), KeyEnvVar) {
		t.Errorf("short %s accepted: %v", KeyEnvVar, err)
	}

	t.Setenv(KeyEnvVar, "")
	if got, err := ResolveKey(""); err != nil || got != nil {
		t.Errorf("no key configured = %q, %v; want nil, nil", got, err)
	}
}
