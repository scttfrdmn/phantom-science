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

// TestLegacyEngine_Close_ViaNewLegacy tests Close on the legacy engine.
// While the function body has no statements, calling it exercises the method dispatch.
func TestLegacyEngine_Close_ViaNewLegacy(t *testing.T) {
	engine := fingerprinting.NewLegacy("close-test-cust")
	require.NotNil(t, engine)
	// Close should not panic and should be safe to call multiple times
	engine.Close()
	engine.Close()
}

// TestLegacyEngine_EmbedInQualities_AllQualityBranches tests the quality score
// branching in EmbedInQualities (>= 30 vs < 30, and fingerprint[i] > 'M' vs <= 'M').
func TestLegacyEngine_EmbedInQualities_AllBranches(t *testing.T) {
	engine := fingerprinting.NewLegacy("branch-cust")

	// Qualities mixing below-30 and above-30 values
	// Quality 28 (< 30) should be unchanged
	// Quality 32 (>= 30) should have LSB set/cleared based on fingerprint char
	qualities := []byte{28, 32, 35, 20, 40, 25}

	// Fingerprint with chars both > 'M' (78) and <= 'M' (e.g., 'A'=65)
	fingerprint := "MAN" // 'M'=77<=77, 'A'=65<='M'=77, 'N'=78>'M'

	out, err := engine.EmbedInQualities(qualities, fingerprint)
	require.NoError(t, err)
	assert.Equal(t, len(qualities), len(out))

	// Verify that quality scores below 30 are not modified
	assert.Equal(t, byte(28), out[0]) // unchanged: 28 < 30
	assert.Equal(t, byte(20), out[3]) // unchanged: 20 < 30
	assert.Equal(t, byte(25), out[5]) // unchanged: 25 < 30
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
