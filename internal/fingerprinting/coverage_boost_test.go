//go:build !release
// +build !release

package fingerprinting_test

import (
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestFingerprintingError_WithContext_NilContext tests the nil Context branch in WithContext.
// When a FingerprintingError is created manually with Context=nil, WithContext should
// initialize the map before setting the value.
func TestFingerprintingError_WithContext_NilContext(t *testing.T) {
	// Create a FingerprintingError with Context explicitly set to nil
	fpErr := &fingerprinting.FingerprintingError{
		Op:      "test_op",
		Err:     fingerprinting.ErrEmptyInput,
		Code:    fingerprinting.ErrorCodeValidation,
		Context: nil,
	}

	// WithContext should initialize the nil map before inserting
	result := fpErr.WithContext("key1", "value1")
	require.NotNil(t, result)
	assert.NotNil(t, result.Context)
	assert.Equal(t, "value1", result.Context["key1"])

	// Chain another WithContext call
	result2 := result.WithContext("key2", 42)
	assert.Equal(t, 42, result2.Context["key2"])
}

// TestIsValidFASTQData_EdgeCases tests edge cases in FASTQ validation.
// These test the isValidFASTQData function's uncovered branches via ValidateInput.
func TestIsValidFASTQData_ViaValidateInput_EdgeCases(t *testing.T) {
	// Valid FASTQ with enough reads to pass capacity check
	var validBuilder strings.Builder
	for i := 0; i < 10; i++ {
		validBuilder.WriteString("@seq1\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
	}
	validFASTQ := []byte(validBuilder.String())
	err := fingerprinting.ValidateInput(validFASTQ, "fp", "test")
	assert.NoError(t, err)

	// FASTQ that starts with @ but has insufficient lines
	tooFewLines := []byte("@seq1\nATCG\n")
	err = fingerprinting.ValidateInput(tooFewLines, "fp", "test")
	assert.Error(t, err)

	// Data that doesn't start with @
	noAtSign := []byte("FASTQ without @ symbol\nATCG\n+\nIIII\n")
	err = fingerprinting.ValidateInput(noAtSign, "fp", "test")
	assert.Error(t, err)
}
