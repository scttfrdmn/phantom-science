//go:build !release

package fingerprinting_test

import (
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// --- Error types tests ---

func TestFingerprintingError_Error(t *testing.T) {
	err := fingerprinting.NewFingerprintingError("embed", fingerprinting.ErrEmptyInput, fingerprinting.ErrorCodeValidation)
	msg := err.Error()
	assert.Contains(t, msg, "embed")
	assert.Contains(t, msg, "empty input")
}

func TestFingerprintingError_WithContext(t *testing.T) {
	err := fingerprinting.NewFingerprintingError("op", fingerprinting.ErrEmbeddingFailed, fingerprinting.ErrorCodeSystem)
	err.WithContext("key", "value").WithContext("size", 42)
	msg := err.Error()
	assert.Contains(t, msg, "key=value")
	assert.Contains(t, msg, "size=42")
}

func TestFingerprintingError_WithRecoverable(t *testing.T) {
	err := fingerprinting.NewFingerprintingError("op", fingerprinting.ErrTimeout, fingerprinting.ErrorCodeTimeout)
	assert.True(t, err.Recoverable)
	err.WithRecoverable(false)
	assert.False(t, err.Recoverable)
}

func TestFingerprintingError_Unwrap(t *testing.T) {
	base := fingerprinting.ErrEmptyInput
	wrapped := fingerprinting.NewFingerprintingError("op", base, fingerprinting.ErrorCodeValidation)
	assert.Equal(t, base, wrapped.Unwrap())
}

func TestFingerprintingError_Is(t *testing.T) {
	wrapped := fingerprinting.NewFingerprintingError("op", fingerprinting.ErrEmptyInput, fingerprinting.ErrorCodeValidation)
	assert.True(t, wrapped.Is(fingerprinting.ErrEmptyInput))
	assert.False(t, wrapped.Is(fingerprinting.ErrTimeout))
	assert.False(t, wrapped.Is(nil))
}

func TestNewFingerprintingError_RecoverableCodes(t *testing.T) {
	tests := []struct {
		code        fingerprinting.ErrorCode
		recoverable bool
	}{
		{fingerprinting.ErrorCodeGPU, true},
		{fingerprinting.ErrorCodeTimeout, true},
		{fingerprinting.ErrorCodeCapacity, true},
		{fingerprinting.ErrorCodeValidation, false},
		{fingerprinting.ErrorCodeCorruption, false},
		{fingerprinting.ErrorCodeSystem, false},
		{fingerprinting.ErrorCodeUnknown, false},
	}
	for _, tc := range tests {
		err := fingerprinting.NewFingerprintingError("op", fingerprinting.ErrTimeout, tc.code)
		assert.Equal(t, tc.recoverable, err.Recoverable, "code %d", tc.code)
	}
}

// --- ValidateInput tests ---

func makeFASTQ(seq string) []byte {
	return []byte("@read1\n" + seq + "\n+\n" + strings.Repeat("I", len(seq)) + "\n")
}

func TestValidateInput_Valid(t *testing.T) {
	data := makeFASTQ("ATCGATCGATCG")
	err := fingerprinting.ValidateInput(data, "fp:123", "test")
	assert.NoError(t, err)
}

func TestValidateInput_EmptyData(t *testing.T) {
	err := fingerprinting.ValidateInput([]byte{}, "fp:123", "test")
	assert.Error(t, err)
}

func TestValidateInput_EmptyFingerprint(t *testing.T) {
	data := makeFASTQ("ATCGATCG")
	err := fingerprinting.ValidateInput(data, "", "test")
	assert.Error(t, err)
}

func TestValidateInput_InvalidFormat(t *testing.T) {
	err := fingerprinting.ValidateInput([]byte("not fastq data at all"), "fp:123", "test")
	assert.Error(t, err)
}

func TestValidateInput_FingerprintTooLong(t *testing.T) {
	data := makeFASTQ("ATCGATCG")
	fp := strings.Repeat("x", fingerprinting.MaxFingerprintLength+1)
	err := fingerprinting.ValidateInput(data, fp, "test")
	assert.Error(t, err)
}

func TestValidateInput_InvalidCharacters(t *testing.T) {
	data := makeFASTQ("ATCGATCG")
	err := fingerprinting.ValidateInput(data, "fp\x00null", "test")
	assert.Error(t, err)
}

// --- EstimateCapacity tests ---

func TestEstimateCapacity_Valid(t *testing.T) {
	// 10KB of data, small fingerprint
	data := make([]byte, 10000)
	err := fingerprinting.EstimateCapacity(data, "shortfp", "balanced")
	assert.NoError(t, err)
}

func TestEstimateCapacity_Exceeded(t *testing.T) {
	// Tiny data, large fingerprint
	data := make([]byte, 100)
	fp := strings.Repeat("x", 200)
	err := fingerprinting.EstimateCapacity(data, fp, "balanced")
	assert.Error(t, err)
}

// --- RecoverableOperation tests ---

func TestRecoverableOperation_Success(t *testing.T) {
	calls := 0
	result, err := fingerprinting.RecoverableOperation(
		func() (int, error) { calls++; return 42, nil },
		nil,
		3,
	)
	require.NoError(t, err)
	assert.Equal(t, 42, result)
	assert.Equal(t, 1, calls)
}

func TestRecoverableOperation_FallbackSuccess(t *testing.T) {
	fpErr := fingerprinting.NewFingerprintingError("op", fingerprinting.ErrGPUError, fingerprinting.ErrorCodeGPU)
	result, err := fingerprinting.RecoverableOperation(
		func() (string, error) { return "", fpErr },
		func() (string, error) { return "fallback", nil },
		1,
	)
	require.NoError(t, err)
	assert.Equal(t, "fallback", result)
}

func TestRecoverableOperation_BothFail(t *testing.T) {
	fpErr := fingerprinting.NewFingerprintingError("op", fingerprinting.ErrCorruptedData, fingerprinting.ErrorCodeCorruption)
	_, err := fingerprinting.RecoverableOperation(
		func() (string, error) { return "", fpErr },
		func() (string, error) { return "", fingerprinting.ErrExtractionFailed },
		0,
	)
	assert.Error(t, err)
}
