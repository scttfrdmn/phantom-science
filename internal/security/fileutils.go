// Package security has path-checked wrappers for creating files (mode 0600)
// and directories (mode 0750).
//
// A path is rejected if any of its components is "..", so names such as
// "a..b" are allowed. With an allowed base directory, the path must lie
// inside it by path components (base "/data" does not admit "/data2/x"),
// compared after resolving symbolic links in the parts of both paths that
// exist, so a link inside the base cannot point outside it.
package security

import (
	"errors"
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"strings"
)

// SecureFileCreate creates a file with security validation to prevent directory traversal
func SecureFileCreate(filename string, allowedBaseDir string) (*os.File, error) {
	cleanPath, err := checkPath(filename, allowedBaseDir)
	if err != nil {
		return nil, err
	}
	// Create with secure permissions (0600 - owner read/write only)
	return os.OpenFile(cleanPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0600)
}

// SecureMkdirAll creates directories with security validation to prevent directory traversal
func SecureMkdirAll(dir string, allowedBaseDir string) error {
	cleanPath, err := checkPath(dir, allowedBaseDir)
	if err != nil {
		return err
	}
	// Create with secure permissions (0750 - owner read/write/execute, group read/execute)
	return os.MkdirAll(cleanPath, 0750)
}

// SecureWriteFile writes a file with security validation to prevent directory traversal
func SecureWriteFile(filename string, data []byte, allowedBaseDir string) error {
	cleanPath, err := checkPath(filename, allowedBaseDir)
	if err != nil {
		return err
	}
	// Write with secure permissions (0600 - owner read/write only)
	return os.WriteFile(cleanPath, data, 0600)
}

// ValidateFilePath validates a file path for security concerns.
//
// A path is rejected if, after filepath.Clean, it still has a ".." path
// component (a relative path that climbs above its starting directory).
// Names that merely contain ".." (for example "a..b.fastq") are allowed.
//
// If allowedBaseDir is non-empty, the path must also be allowedBaseDir
// itself or lie beneath it. Containment is checked on path components, not
// string prefixes (base "/data" does not admit "/data2/x"), and symlinks in
// the existing part of either path are resolved first, so a symlink inside
// the base that points outside it is rejected.
func ValidateFilePath(path string, allowedBaseDir string) error {
	_, err := checkPath(path, allowedBaseDir)
	return err
}

// checkPath implements ValidateFilePath and returns the cleaned path.
func checkPath(path, allowedBaseDir string) (string, error) {
	cleanPath := filepath.Clean(path)

	if hasDotDotComponent(cleanPath) {
		return "", fmt.Errorf("path traversal not allowed: %s", path)
	}

	if allowedBaseDir == "" {
		return cleanPath, nil
	}

	absPath, err := filepath.Abs(cleanPath)
	if err != nil {
		return "", fmt.Errorf("failed to resolve path %s: %w", path, err)
	}
	absAllowed, err := filepath.Abs(allowedBaseDir)
	if err != nil {
		return "", fmt.Errorf("failed to resolve allowed directory %s: %w", allowedBaseDir, err)
	}

	realPath, err := resolveExisting(absPath)
	if err != nil {
		return "", fmt.Errorf("failed to resolve path %s: %w", path, err)
	}
	realAllowed, err := resolveExisting(absAllowed)
	if err != nil {
		return "", fmt.Errorf("failed to resolve allowed directory %s: %w", allowedBaseDir, err)
	}

	if !within(realAllowed, realPath) {
		return "", fmt.Errorf("path %s outside allowed directory %s", path, allowedBaseDir)
	}
	return cleanPath, nil
}

// hasDotDotComponent reports whether a cleaned path has a ".." element.
func hasDotDotComponent(cleanPath string) bool {
	for _, part := range strings.Split(filepath.ToSlash(cleanPath), "/") {
		if part == ".." {
			return true
		}
	}
	return false
}

// within reports whether target is base or lies beneath it. Both must be
// absolute and clean.
func within(base, target string) bool {
	rel, err := filepath.Rel(base, target)
	if err != nil {
		return false
	}
	return rel != ".." && !strings.HasPrefix(rel, ".."+string(filepath.Separator)) && !filepath.IsAbs(rel)
}

// resolveExisting resolves symlinks in the longest existing prefix of an
// absolute path and re-appends the not-yet-existing remainder.
func resolveExisting(absPath string) (string, error) {
	existing := absPath
	var rest []string
	for {
		resolved, err := filepath.EvalSymlinks(existing)
		if err == nil {
			return filepath.Join(append([]string{resolved}, rest...)...), nil
		}
		if !errors.Is(err, fs.ErrNotExist) {
			return "", err
		}
		parent := filepath.Dir(existing)
		if parent == existing {
			// Nothing exists, not even the root; compare lexically.
			return absPath, nil
		}
		rest = append([]string{filepath.Base(existing)}, rest...)
		existing = parent
	}
}
