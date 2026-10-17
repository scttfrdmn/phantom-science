package security

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// SecureFileCreate creates a file with security validation to prevent directory traversal
func SecureFileCreate(filename string, allowedBaseDir string) (*os.File, error) {
	// Clean and validate the path
	cleanPath := filepath.Clean(filename)
	
	// Prevent directory traversal
	if strings.Contains(cleanPath, "..") {
		return nil, fmt.Errorf("path traversal not allowed: %s", filename)
	}
	
	// Ensure path is within allowed directory if specified
	if allowedBaseDir != "" {
		absPath, err := filepath.Abs(cleanPath)
		if err != nil {
			return nil, fmt.Errorf("failed to resolve path %s: %w", filename, err)
		}
		
		absAllowed, err := filepath.Abs(allowedBaseDir)
		if err != nil {
			return nil, fmt.Errorf("failed to resolve allowed directory %s: %w", allowedBaseDir, err)
		}
		
		if !strings.HasPrefix(absPath, absAllowed) {
			return nil, fmt.Errorf("path %s outside allowed directory %s", filename, allowedBaseDir)
		}
	}
	
	// Create with secure permissions (0600 - owner read/write only)
	return os.OpenFile(cleanPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0600)
}

// SecureMkdirAll creates directories with security validation to prevent directory traversal
func SecureMkdirAll(dir string, allowedBaseDir string) error {
	// Clean and validate the path
	cleanPath := filepath.Clean(dir)
	
	// Prevent directory traversal
	if strings.Contains(cleanPath, "..") {
		return fmt.Errorf("path traversal not allowed: %s", dir)
	}
	
	// Ensure path is within allowed directory if specified
	if allowedBaseDir != "" {
		absPath, err := filepath.Abs(cleanPath)
		if err != nil {
			return fmt.Errorf("failed to resolve path %s: %w", dir, err)
		}
		
		absAllowed, err := filepath.Abs(allowedBaseDir)
		if err != nil {
			return fmt.Errorf("failed to resolve allowed directory %s: %w", allowedBaseDir, err)
		}
		
		if !strings.HasPrefix(absPath, absAllowed) {
			return fmt.Errorf("path %s outside allowed directory %s", dir, allowedBaseDir)
		}
	}
	
	// Create with secure permissions (0750 - owner read/write/execute, group read/execute)
	return os.MkdirAll(cleanPath, 0750)
}

// SecureWriteFile writes a file with security validation to prevent directory traversal
func SecureWriteFile(filename string, data []byte, allowedBaseDir string) error {
	// Clean and validate the path
	cleanPath := filepath.Clean(filename)
	
	// Prevent directory traversal
	if strings.Contains(cleanPath, "..") {
		return fmt.Errorf("path traversal not allowed: %s", filename)
	}
	
	// Ensure path is within allowed directory if specified
	if allowedBaseDir != "" {
		absPath, err := filepath.Abs(cleanPath)
		if err != nil {
			return fmt.Errorf("failed to resolve path %s: %w", filename, err)
		}
		
		absAllowed, err := filepath.Abs(allowedBaseDir)
		if err != nil {
			return fmt.Errorf("failed to resolve allowed directory %s: %w", allowedBaseDir, err)
		}
		
		if !strings.HasPrefix(absPath, absAllowed) {
			return fmt.Errorf("path %s outside allowed directory %s", filename, allowedBaseDir)
		}
	}
	
	// Write with secure permissions (0600 - owner read/write only)
	return os.WriteFile(cleanPath, data, 0600)
}

// ValidateFilePath validates a file path for security concerns
func ValidateFilePath(path string, allowedBaseDir string) error {
	// Clean and validate the path
	cleanPath := filepath.Clean(path)
	
	// Prevent directory traversal
	if strings.Contains(cleanPath, "..") {
		return fmt.Errorf("path traversal not allowed: %s", path)
	}
	
	// Ensure path is within allowed directory if specified
	if allowedBaseDir != "" {
		absPath, err := filepath.Abs(cleanPath)
		if err != nil {
			return fmt.Errorf("failed to resolve path %s: %w", path, err)
		}
		
		absAllowed, err := filepath.Abs(allowedBaseDir)
		if err != nil {
			return fmt.Errorf("failed to resolve allowed directory %s: %w", allowedBaseDir, err)
		}
		
		if !strings.HasPrefix(absPath, absAllowed) {
			return fmt.Errorf("path %s outside allowed directory %s", path, allowedBaseDir)
		}
	}
	
	return nil
}