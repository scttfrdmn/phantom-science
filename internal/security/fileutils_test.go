package security

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestValidateFilePath(t *testing.T) {
	tests := []struct {
		name          string
		path          string
		allowedBase   string
		expectError   bool
		errorContains string
	}{
		{
			name:        "valid simple path",
			path:        "test.txt",
			allowedBase: "",
			expectError: false,
		},
		{
			name:          "path traversal attack",
			path:          "../../../etc/passwd",
			allowedBase:   "",
			expectError:   true,
			errorContains: "path traversal not allowed",
		},
		{
			name:          "path with .. in middle",
			path:          "dir/../../../secret",
			allowedBase:   "",
			expectError:   true,
			errorContains: "path traversal not allowed",
		},
		{
			name:        "valid path within allowed base",
			path:        "/tmp/test/file.txt",
			allowedBase: "/tmp",
			expectError: false,
		},
		{
			name:          "path outside allowed base",
			path:          "/etc/passwd",
			allowedBase:   "/tmp",
			expectError:   true,
			errorContains: "outside allowed directory",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			err := ValidateFilePath(tt.path, tt.allowedBase)
			
			if tt.expectError {
				if err == nil {
					t.Errorf("Expected error for path %s, got none", tt.path)
				} else if tt.errorContains != "" && !strings.Contains(err.Error(), tt.errorContains) {
					t.Errorf("Expected error to contain %s, got: %s", tt.errorContains, err.Error())
				}
			} else {
				if err != nil {
					t.Errorf("Expected no error for path %s, got: %v", tt.path, err)
				}
			}
		})
	}
}

func TestSecureFileCreate(t *testing.T) {
	tmpDir := t.TempDir()
	
	tests := []struct {
		name          string
		filename      string
		allowedBase   string
		expectError   bool
		errorContains string
	}{
		{
			name:        "valid file creation",
			filename:    filepath.Join(tmpDir, "test.txt"),
			allowedBase: tmpDir,
			expectError: false,
		},
		{
			name:          "path traversal in filename",
			filename:      filepath.Join(tmpDir, "../../../etc/malicious"),
			allowedBase:   tmpDir,
			expectError:   true,
			errorContains: "outside allowed directory",
		},
		{
			name:          "file outside allowed directory",
			filename:      "/tmp/outside.txt",
			allowedBase:   tmpDir,
			expectError:   true,
			errorContains: "outside allowed directory",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			file, err := SecureFileCreate(tt.filename, tt.allowedBase)
			
			if tt.expectError {
				if err == nil {
					t.Errorf("Expected error for filename %s, got none", tt.filename)
					if file != nil {
						file.Close()
						os.Remove(tt.filename)
					}
				} else if tt.errorContains != "" && !strings.Contains(err.Error(), tt.errorContains) {
					t.Errorf("Expected error to contain %s, got: %s", tt.errorContains, err.Error())
				}
			} else {
				if err != nil {
					t.Errorf("Expected no error for filename %s, got: %v", tt.filename, err)
				} else {
					if file == nil {
						t.Error("Expected file to be created")
					} else {
						file.Close()
						
						// Check file was created with secure permissions
						info, err := os.Stat(tt.filename)
						if err != nil {
							t.Errorf("File should exist: %v", err)
						} else {
							mode := info.Mode().Perm()
							expectedMode := os.FileMode(0600)
							if mode != expectedMode {
								t.Errorf("Expected file mode %o, got %o", expectedMode, mode)
							}
						}
						
						// Clean up
						os.Remove(tt.filename)
					}
				}
			}
		})
	}
}

func TestSecureMkdirAll(t *testing.T) {
	tmpDir := t.TempDir()
	
	tests := []struct {
		name          string
		dir           string
		allowedBase   string
		expectError   bool
		errorContains string
	}{
		{
			name:        "valid directory creation",
			dir:         filepath.Join(tmpDir, "subdir", "nested"),
			allowedBase: tmpDir,
			expectError: false,
		},
		{
			name:          "path traversal in directory",
			dir:           filepath.Join(tmpDir, "../../../etc"),
			allowedBase:   tmpDir,
			expectError:   true,
			errorContains: "outside allowed directory",
		},
		{
			name:          "directory outside allowed base",
			dir:           "/tmp/outside",
			allowedBase:   tmpDir,
			expectError:   true,
			errorContains: "outside allowed directory",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			err := SecureMkdirAll(tt.dir, tt.allowedBase)
			
			if tt.expectError {
				if err == nil {
					t.Errorf("Expected error for dir %s, got none", tt.dir)
					// Clean up if directory was somehow created
					os.RemoveAll(tt.dir)
				} else if tt.errorContains != "" && !strings.Contains(err.Error(), tt.errorContains) {
					t.Errorf("Expected error to contain %s, got: %s", tt.errorContains, err.Error())
				}
			} else {
				if err != nil {
					t.Errorf("Expected no error for dir %s, got: %v", tt.dir, err)
				} else {
					// Check directory was created with secure permissions
					info, err := os.Stat(tt.dir)
					if err != nil {
						t.Errorf("Directory should exist: %v", err)
					} else {
						mode := info.Mode().Perm()
						expectedMode := os.FileMode(0750)
						if mode != expectedMode {
							t.Errorf("Expected directory mode %o, got %o", expectedMode, mode)
						}
					}
					
					// Clean up
					os.RemoveAll(tt.dir)
				}
			}
		})
	}
}

func TestSecureWriteFile(t *testing.T) {
	tmpDir := t.TempDir()
	testData := []byte("test data content")
	
	tests := []struct {
		name          string
		filename      string
		allowedBase   string
		expectError   bool
		errorContains string
	}{
		{
			name:        "valid file write",
			filename:    filepath.Join(tmpDir, "test.txt"),
			allowedBase: tmpDir,
			expectError: false,
		},
		{
			name:          "path traversal in filename",
			filename:      filepath.Join(tmpDir, "../../../etc/malicious"),
			allowedBase:   tmpDir,
			expectError:   true,
			errorContains: "outside allowed directory",
		},
		{
			name:          "file outside allowed directory",
			filename:      "/tmp/outside.txt",
			allowedBase:   tmpDir,
			expectError:   true,
			errorContains: "outside allowed directory",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			err := SecureWriteFile(tt.filename, testData, tt.allowedBase)
			
			if tt.expectError {
				if err == nil {
					t.Errorf("Expected error for filename %s, got none", tt.filename)
					// Clean up if file was somehow created
					os.Remove(tt.filename)
				} else if tt.errorContains != "" && !strings.Contains(err.Error(), tt.errorContains) {
					t.Errorf("Expected error to contain %s, got: %s", tt.errorContains, err.Error())
				}
			} else {
				if err != nil {
					t.Errorf("Expected no error for filename %s, got: %v", tt.filename, err)
				} else {
					// Check file was created with secure permissions and correct content
					info, err := os.Stat(tt.filename)
					if err != nil {
						t.Errorf("File should exist: %v", err)
					} else {
						mode := info.Mode().Perm()
						expectedMode := os.FileMode(0600)
						if mode != expectedMode {
							t.Errorf("Expected file mode %o, got %o", expectedMode, mode)
						}
						
						// Check content
						content, err := os.ReadFile(tt.filename)
						if err != nil {
							t.Errorf("Should be able to read file: %v", err)
						} else if string(content) != string(testData) {
							t.Errorf("File content mismatch. Expected %s, got %s", testData, content)
						}
					}
					
					// Clean up
					os.Remove(tt.filename)
				}
			}
		})
	}
}

func TestSecurityWithEmptyAllowedBase(t *testing.T) {
	// Test that functions work correctly when no allowedBase is specified
	tmpDir := t.TempDir()

	// Test ValidateFilePath with empty allowedBase
	err := ValidateFilePath("test.txt", "")
	if err != nil {
		t.Errorf("ValidateFilePath should allow any path when allowedBase is empty, got: %v", err)
	}

	// Test with path traversal - should still be blocked
	err = ValidateFilePath("../../../etc/passwd", "")
	if err == nil {
		t.Error("ValidateFilePath should block path traversal even with empty allowedBase")
	}

	// Test SecureFileCreate with empty allowedBase
	testFile := filepath.Join(tmpDir, "test.txt")
	file, err := SecureFileCreate(testFile, "")
	if err != nil {
		t.Errorf("SecureFileCreate should work with empty allowedBase: %v", err)
	} else {
		file.Close()
		os.Remove(testFile)
	}
}

// TestPathTraversalInFilename tests that raw ".." in path is blocked.
// We pass a path that, after filepath.Clean, still contains ".." elements.
func TestSecureFileCreate_PathTraversalDotDot(t *testing.T) {
	// Pass a raw relative path with ".." that Clean will NOT fully resolve to an absolute path.
	// "subdir/../../etc/passwd" after Clean becomes "../etc/passwd" which contains ".."
	_, err := SecureFileCreate("subdir/../../etc/passwd", "")
	if err == nil {
		t.Error("SecureFileCreate should reject path with .. after Clean")
	}
	if err != nil && !strings.Contains(err.Error(), "path traversal") {
		t.Errorf("expected path traversal error, got: %v", err)
	}
}

func TestSecureMkdirAll_PathTraversalDotDot(t *testing.T) {
	err := SecureMkdirAll("subdir/../../etc", "")
	if err == nil {
		t.Error("SecureMkdirAll should reject path with .. after Clean")
	}
	if err != nil && !strings.Contains(err.Error(), "path traversal") {
		t.Errorf("expected path traversal error, got: %v", err)
	}
}

func TestSecureWriteFile_PathTraversalDotDot(t *testing.T) {
	err := SecureWriteFile("subdir/../../etc/passwd", []byte("test"), "")
	if err == nil {
		t.Error("SecureWriteFile should reject path with .. after Clean")
	}
	if err != nil && !strings.Contains(err.Error(), "path traversal") {
		t.Errorf("expected path traversal error, got: %v", err)
	}
}