//go:build !release

package fingerprinting_test

import (
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// --- LegacyEngine tests ---

func TestNewLegacy(t *testing.T) {
	e := fingerprinting.NewLegacy("customer-123")
	require.NotNil(t, e)
}

func TestLegacyEngine_GenerateFingerprint(t *testing.T) {
	e := fingerprinting.NewLegacy("acme")
	fp, err := e.GenerateFingerprint("dataset-1")
	require.NoError(t, err)
	assert.Contains(t, fp, "acme")
	assert.Contains(t, fp, "dataset-1")
}

func TestLegacyEngine_EmbedInQualities(t *testing.T) {
	e := fingerprinting.NewLegacy("cust")
	qualities := []byte{40, 38, 35, 30, 32, 40}
	out, err := e.EmbedInQualities(qualities, "testfp")
	require.NoError(t, err)
	assert.Equal(t, len(qualities), len(out))
}

func TestLegacyEngine_EmbedInQualities_Empty(t *testing.T) {
	e := fingerprinting.NewLegacy("cust")
	out, err := e.EmbedInQualities([]byte{}, "testfp")
	require.NoError(t, err)
	assert.Empty(t, out)
}

func TestLegacyEngine_EmbedCustomerID(t *testing.T) {
	e := fingerprinting.NewLegacy("cust-42")
	result := e.EmbedCustomerID("read001")
	assert.Contains(t, result, "read001")
	assert.Contains(t, result, "cust-42")
}

func TestLegacyEngine_EmbedDatasetFingerprint(t *testing.T) {
	e := fingerprinting.NewLegacy("cust")
	result := e.EmbedDatasetFingerprint("read001", "ds-99")
	assert.Contains(t, result, "read001")
	assert.Contains(t, result, "ds-99")
}

func TestLegacyEngine_VerifyFingerprint(t *testing.T) {
	e := fingerprinting.NewLegacy("cust")
	assert.True(t, e.VerifyFingerprint([]byte{40, 38}))
	assert.False(t, e.VerifyFingerprint([]byte{}))
	assert.False(t, e.VerifyFingerprint([]byte{20, 25})) // below threshold 30
}

func TestLegacyEngine_ExtractFingerprint(t *testing.T) {
	e := fingerprinting.NewLegacy("cust-1")
	fp, err := e.ExtractFingerprint([]byte{40, 38, 35})
	require.NoError(t, err)
	assert.Contains(t, fp, "cust-1")
}

func TestLegacyEngine_ExtractFingerprint_Empty(t *testing.T) {
	e := fingerprinting.NewLegacy("cust")
	_, err := e.ExtractFingerprint([]byte{})
	assert.Error(t, err)
}

func TestLegacyEngine_Close(t *testing.T) {
	e := fingerprinting.NewLegacy("cust")
	e.Close() // should not panic
}

// --- LegacyVerifier tests ---

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
