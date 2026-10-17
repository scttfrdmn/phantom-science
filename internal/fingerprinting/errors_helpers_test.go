//go:build !release && cgo
// +build !release,cgo

package fingerprinting

import (
	"testing"
)

// TestIsValidFASTQData verifies the FASTQ format validator across valid, invalid,
// empty, and malformed inputs.
func TestIsValidFASTQData(t *testing.T) {
	tests := []struct {
		name  string
		data  []byte
		want  bool
	}{
		{
			name: "valid single-record FASTQ",
			data: []byte("@read1\nATCGATCGATCG\n+\nIIIIIIIIIIII\n"),
			want: true,
		},
		{
			name: "valid multi-record FASTQ",
			data: []byte("@read1\nATCGATCGATCG\n+\nIIIIIIIIIIII\n@read2\nGGGGCCCCAAAA\n+\nJJJJJJJJJJJJ\n"),
			want: true,
		},
		{
			name: "valid FASTQ with sequence ID in + line",
			data: []byte("@SRR001666.1\nATCGATCGATCGATCG\n+SRR001666.1\nIIIIIIIIIIIIIIII\n"),
			want: true,
		},
		{
			name:  "empty data",
			data:  []byte{},
			want:  false,
		},
		{
			name:  "nil data treated as empty",
			data:  nil,
			want:  false,
		},
		{
			name:  "data shorter than MinDataLength",
			data:  []byte("@s\nA\n+\nI\n"),
			want:  false,
		},
		{
			name: "starts with @ but missing + separator",
			// Has at least 4 lines but no + separator line
			data: []byte("@read1\nATCGATCGATCGATCGATCG\nsome_line\nanother_line\n"),
			want: false,
		},
		{
			name: "does not start with @",
			// FASTA format instead of FASTQ
			data: []byte(">read1\nATCGATCGATCGATCGATCG\nATCGATCGATCGATCGATCG\nATCGATCGATCGATCGATCG\n"),
			want: false,
		},
		{
			name: "fewer than 4 lines - only header and sequence",
			// Only 2 newline-separated parts: "@read1" and "ATCG" — splits to 3 items,
			// no "+" separator present
			data: []byte("@read1\nATCGATCGATCGATCGATCGATCGATCGATCG"),
			want: false,
		},
		{
			name: "whitespace-prefixed @ still valid",
			// First line has leading whitespace; TrimSpace handles it
			data: []byte("  @read1\nATCGATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIIIIII\n"),
			want: true,
		},
		{
			name: "+ separator on a later line",
			// The + appears as the third line — standard FASTQ layout
			data: []byte("@read1\nATCGATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIIIIII\n"),
			want: true,
		},
		{
			name: "multiple reads, + in middle records",
			data: []byte(
				"@r1\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n" +
					"@r2\nCCCCCCCCCCCCCCCCCCCC\n+\nJJJJJJJJJJJJJJJJJJJJ\n",
			),
			want: true,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := isValidFASTQData(tc.data)
			if got != tc.want {
				t.Errorf("isValidFASTQData(%q) = %v, want %v", tc.data, got, tc.want)
			}
		})
	}
}

// TestContainsInvalidCharacters verifies that the character validator correctly
// identifies printable ASCII, control characters, null bytes, and extended bytes.
func TestContainsInvalidCharacters(t *testing.T) {
	tests := []struct {
		name        string
		fingerprint string
		want        bool // true means contains invalid characters
	}{
		{
			name:        "empty string has no invalid characters",
			fingerprint: "",
			want:        false,
		},
		{
			name:        "printable ASCII only",
			fingerprint: "customer-id:dataset:2024",
			want:        false,
		},
		{
			name:        "printable ASCII full range boundary low (space=32)",
			fingerprint: " ",
			want:        false,
		},
		{
			name:        "printable ASCII full range boundary high (tilde=126)",
			fingerprint: "~",
			want:        false,
		},
		{
			name:        "alphanumeric plus symbols",
			fingerprint: "abcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()",
			want:        false,
		},
		{
			name:        "null byte is invalid",
			fingerprint: "abc\x00def",
			want:        true,
		},
		{
			name:        "newline (control char) is invalid",
			fingerprint: "abc\ndef",
			want:        true,
		},
		{
			name:        "tab (control char) is invalid",
			fingerprint: "abc\tdef",
			want:        true,
		},
		{
			name:        "carriage return is invalid",
			fingerprint: "abc\rdef",
			want:        true,
		},
		{
			name:        "DEL character (127) is invalid",
			fingerprint: "abc\x7fdef",
			want:        true,
		},
		{
			name:        "extended ASCII byte (128) is invalid",
			fingerprint: "abc\x80def",
			want:        true,
		},
		{
			name:        "extended ASCII byte (255) is invalid",
			fingerprint: "abc\xffdef",
			want:        true,
		},
		{
			name:        "SOH control character (1) is invalid",
			fingerprint: "\x01",
			want:        true,
		},
		{
			name:        "escape character (27) is invalid",
			fingerprint: "\x1b[1;31m",
			want:        true,
		},
		{
			name:        "multi-rune Unicode is invalid (non-ASCII runes > 127)",
			fingerprint: "café",
			want:        true,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := containsInvalidCharacters(tc.fingerprint)
			if got != tc.want {
				t.Errorf("containsInvalidCharacters(%q) = %v, want %v", tc.fingerprint, got, tc.want)
			}
		})
	}
}

// TestIsRecoverableError verifies the recoverability mapping for every defined ErrorCode.
func TestIsRecoverableError(t *testing.T) {
	tests := []struct {
		name string
		code ErrorCode
		want bool
	}{
		// Recoverable: can fallback to CPU or retry
		{name: "GPU error is recoverable", code: ErrorCodeGPU, want: true},
		{name: "Timeout is recoverable", code: ErrorCodeTimeout, want: true},
		{name: "Capacity exceeded is recoverable", code: ErrorCodeCapacity, want: true},

		// Not recoverable: data issues or system-level failures
		{name: "Validation error is not recoverable", code: ErrorCodeValidation, want: false},
		{name: "Corruption is not recoverable", code: ErrorCodeCorruption, want: false},
		{name: "System error is not recoverable", code: ErrorCodeSystem, want: false},

		// Default (unknown) falls through to not recoverable
		{name: "Unknown code is not recoverable", code: ErrorCodeUnknown, want: false},
		{name: "Permission error falls to default (not recoverable)", code: ErrorCodePermission, want: false},

		// Verify the iota values are what we expect
		{name: "Zero value (Unknown) is not recoverable", code: 0, want: false},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := isRecoverableError(tc.code)
			if got != tc.want {
				t.Errorf("isRecoverableError(%d) = %v, want %v", tc.code, got, tc.want)
			}
		})
	}
}

// TestIsRecoverableError_AllDefinedCodes ensures coverage of every named constant.
func TestIsRecoverableError_AllDefinedCodes(t *testing.T) {
	// Verify exact recoverability contract for every named ErrorCode constant.
	recoverableCodes := map[ErrorCode]bool{
		ErrorCodeUnknown:    false,
		ErrorCodeValidation: false,
		ErrorCodeCapacity:   true,
		ErrorCodeSystem:     false,
		ErrorCodeGPU:        true,
		ErrorCodeTimeout:    true,
		ErrorCodeCorruption: false,
		ErrorCodePermission: false,
	}

	for code, wantRecoverable := range recoverableCodes {
		got := isRecoverableError(code)
		if got != wantRecoverable {
			t.Errorf("isRecoverableError(ErrorCode=%d) = %v, want %v", code, got, wantRecoverable)
		}
	}
}

