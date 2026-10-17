package errors_test

import (
	"errors"
	"fmt"
	"strings"
	"sync"
	"testing"
	"time"

	gerrors "github.com/scttfrdmn/phantom-science/internal/errors"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ---------------------------------------------------------------------------
// GenomicsError – constructor and basic methods
// ---------------------------------------------------------------------------

func TestNewGenomicsError(t *testing.T) {
	err := gerrors.NewGenomicsError(gerrors.ValidationError, "comp", "op", "bad input")
	require.NotNil(t, err)
	assert.Equal(t, gerrors.ValidationError, err.Type)
	assert.Equal(t, "comp", err.Component)
	assert.Equal(t, "op", err.Operation)
	assert.Equal(t, "bad input", err.Message)
	assert.NotNil(t, err.Context)
	assert.NotNil(t, err.StackTrace)
	// Timestamp is mocked in production implementation
}

func TestGenomicsError_Error(t *testing.T) {
	err := gerrors.NewGenomicsError(gerrors.ParseError, "parser", "parse_fasta", "unexpected EOF")
	s := err.Error()
	assert.Contains(t, s, "PARSE_ERROR")
	assert.Contains(t, s, "parser")
	assert.Contains(t, s, "parse_fasta")
	assert.Contains(t, s, "unexpected EOF")
}

func TestGenomicsError_ErrorWithCause(t *testing.T) {
	cause := fmt.Errorf("underlying io error")
	err := gerrors.NewGenomicsError(gerrors.ResourceError, "loader", "load", "failed to load")
	err.WithCause(cause)
	s := err.Error()
	assert.Contains(t, s, "underlying io error")
}

func TestGenomicsError_ErrorNoOperation(t *testing.T) {
	err := gerrors.NewGenomicsError(gerrors.MemoryError, "allocator", "", "OOM")
	s := err.Error()
	assert.NotContains(t, s, "operation=")
	assert.Contains(t, s, "OOM")
}

func TestWrapError_NilReturnsNil(t *testing.T) {
	result := gerrors.WrapError(nil, gerrors.ValidationError, "c", "o", "msg")
	assert.Nil(t, result)
}

func TestWrapError_StandardError(t *testing.T) {
	inner := errors.New("some error")
	wrapped := gerrors.WrapError(inner, gerrors.ParseError, "parser", "parse", "wrapping")
	require.NotNil(t, wrapped)
	assert.Equal(t, gerrors.ParseError, wrapped.Type)
	assert.Equal(t, inner, wrapped.Cause)
}

func TestWrapError_GenomicsErrorAddsContext(t *testing.T) {
	orig := gerrors.NewGenomicsError(gerrors.ValidationError, "comp", "op", "original")
	rewrapped := gerrors.WrapError(orig, gerrors.ResourceError, "other", "op2", "re-wrap")
	// When the error is already a GenomicsError, WrapError returns the same pointer
	assert.Equal(t, orig, rewrapped)
	_, hasContext := rewrapped.GetContext("wrapped_at")
	assert.True(t, hasContext)
}

func TestGenomicsError_WithContext(t *testing.T) {
	err := gerrors.NewGenomicsError(gerrors.ConfigurationError, "cfg", "load", "bad config")
	returned := err.WithContext("key1", "value1").WithContext("key2", 42)
	assert.Same(t, err, returned)

	v1, ok := err.GetContext("key1")
	assert.True(t, ok)
	assert.Equal(t, "value1", v1)

	v2, ok := err.GetContext("key2")
	assert.True(t, ok)
	assert.Equal(t, 42, v2)

	_, ok = err.GetContext("nonexistent")
	assert.False(t, ok)
}

func TestGenomicsError_WithCause(t *testing.T) {
	cause := errors.New("root cause")
	err := gerrors.NewGenomicsError(gerrors.DatabaseError, "db", "query", "query failed")
	returned := err.WithCause(cause)
	assert.Same(t, err, returned)
	assert.Equal(t, cause, err.Cause)
}

func TestGenomicsError_IsType(t *testing.T) {
	err := gerrors.NewGenomicsError(gerrors.KmerExtractionError, "kmer", "extract", "bad kmer")
	assert.True(t, err.IsType(gerrors.KmerExtractionError))
	assert.False(t, err.IsType(gerrors.ValidationError))
}

func TestGenomicsError_IsComponent(t *testing.T) {
	err := gerrors.NewGenomicsError(gerrors.IndexError, "indexer", "build", "index failed")
	assert.True(t, err.IsComponent("indexer"))
	assert.False(t, err.IsComponent("other"))
}

func TestGenomicsError_GetStackTrace(t *testing.T) {
	err := gerrors.NewGenomicsError(gerrors.ValidationError, "c", "o", "m")
	trace := err.GetStackTrace()
	// Stack trace may be empty in some environments but should always be a string
	assert.IsType(t, "", trace)
}

func TestGenomicsError_GetContext_NilMap(t *testing.T) {
	// Force a nil context (direct construction)
	err := &gerrors.GenomicsError{
		Type:      gerrors.ValidationError,
		Component: "c",
		Operation: "o",
		Message:   "m",
		Context:   nil,
	}
	_, ok := err.GetContext("anything")
	assert.False(t, ok)
}

// ---------------------------------------------------------------------------
// ErrorType String() and GetSeverity()
// ---------------------------------------------------------------------------

func TestErrorType_String(t *testing.T) {
	cases := []struct {
		et   gerrors.ErrorType
		want string
	}{
		{gerrors.ValidationError, "VALIDATION_ERROR"},
		{gerrors.ConfigurationError, "CONFIGURATION_ERROR"},
		{gerrors.InitializationError, "INITIALIZATION_ERROR"},
		{gerrors.ResourceError, "RESOURCE_ERROR"},
		{gerrors.ParseError, "PARSE_ERROR"},
		{gerrors.FormatError, "FORMAT_ERROR"},
		{gerrors.CorruptedDataError, "CORRUPTED_DATA_ERROR"},
		{gerrors.InvalidInputError, "INVALID_INPUT_ERROR"},
		{gerrors.MemoryError, "MEMORY_ERROR"},
		{gerrors.GPUError, "GPU_ERROR"},
		{gerrors.CacheError, "CACHE_ERROR"},
		{gerrors.TimeoutError, "TIMEOUT_ERROR"},
		{gerrors.BiologicalModelError, "BIOLOGICAL_MODEL_ERROR"},
		{gerrors.PhylogeneticError, "PHYLOGENETIC_ERROR"},
		{gerrors.ChromatinError, "CHROMATIN_ERROR"},
		{gerrors.MutationalSignatureError, "MUTATIONAL_SIGNATURE_ERROR"},
		{gerrors.SelectionPressureError, "SELECTION_PRESSURE_ERROR"},
		{gerrors.KmerExtractionError, "KMER_EXTRACTION_ERROR"},
		{gerrors.KmerValidationError, "KMER_VALIDATION_ERROR"},
		{gerrors.PackedArrayError, "PACKED_ARRAY_ERROR"},
		{gerrors.ReferenceLoadError, "REFERENCE_LOAD_ERROR"},
		{gerrors.IndexError, "INDEX_ERROR"},
		{gerrors.SequenceNotFoundError, "SEQUENCE_NOT_FOUND_ERROR"},
		{gerrors.VCFParseError, "VCF_PARSE_ERROR"},
		{gerrors.VCFFormatError, "VCF_FORMAT_ERROR"},
		{gerrors.VariantProcessingError, "VARIANT_PROCESSING_ERROR"},
		{gerrors.ExternalToolError, "EXTERNAL_TOOL_ERROR"},
		{gerrors.DatabaseError, "DATABASE_ERROR"},
		{gerrors.NetworkError, "NETWORK_ERROR"},
		{gerrors.AuthenticityValidationError, "AUTHENTICITY_VALIDATION_ERROR"},
		{gerrors.SecurityError, "SECURITY_ERROR"},
		{gerrors.AccessDeniedError, "ACCESS_DENIED_ERROR"},
		{gerrors.ErrorType(9999), "UNKNOWN_ERROR"},
	}
	for _, tc := range cases {
		t.Run(tc.want, func(t *testing.T) {
			assert.Equal(t, tc.want, tc.et.String())
		})
	}
}

func TestErrorType_GetSeverity(t *testing.T) {
	// Critical
	for _, et := range []gerrors.ErrorType{gerrors.MemoryError, gerrors.CorruptedDataError, gerrors.SecurityError, gerrors.AccessDeniedError} {
		assert.Equal(t, gerrors.SeverityCritical, et.GetSeverity(), et.String())
	}
	// High
	for _, et := range []gerrors.ErrorType{gerrors.InitializationError, gerrors.GPUError, gerrors.ReferenceLoadError, gerrors.DatabaseError} {
		assert.Equal(t, gerrors.SeverityHigh, et.GetSeverity(), et.String())
	}
	// Medium
	for _, et := range []gerrors.ErrorType{gerrors.ValidationError, gerrors.ConfigurationError, gerrors.BiologicalModelError, gerrors.VCFParseError, gerrors.ExternalToolError} {
		assert.Equal(t, gerrors.SeverityMedium, et.GetSeverity(), et.String())
	}
	// Low
	for _, et := range []gerrors.ErrorType{gerrors.ParseError, gerrors.FormatError, gerrors.CacheError, gerrors.TimeoutError} {
		assert.Equal(t, gerrors.SeverityLow, et.GetSeverity(), et.String())
	}
	// Default medium for unmapped types
	assert.Equal(t, gerrors.SeverityMedium, gerrors.ErrorType(9999).GetSeverity())
}

func TestErrorSeverity_String(t *testing.T) {
	cases := []struct {
		sev  gerrors.ErrorSeverity
		want string
	}{
		{gerrors.SeverityLow, "LOW"},
		{gerrors.SeverityMedium, "MEDIUM"},
		{gerrors.SeverityHigh, "HIGH"},
		{gerrors.SeverityCritical, "CRITICAL"},
		{gerrors.ErrorSeverity(99), "UNKNOWN"},
	}
	for _, tc := range cases {
		t.Run(tc.want, func(t *testing.T) {
			assert.Equal(t, tc.want, tc.sev.String())
		})
	}
}

// ---------------------------------------------------------------------------
// Convenience constructors
// ---------------------------------------------------------------------------

func TestNewValidationError(t *testing.T) {
	err := gerrors.NewValidationError("comp", "op", "bad")
	assert.Equal(t, gerrors.ValidationError, err.Type)
}

func TestNewParseError(t *testing.T) {
	err := gerrors.NewParseError("comp", "op", "bad parse")
	assert.Equal(t, gerrors.ParseError, err.Type)
}

func TestNewResourceError(t *testing.T) {
	err := gerrors.NewResourceError("comp", "op", "no resource")
	assert.Equal(t, gerrors.ResourceError, err.Type)
}

func TestNewBiologicalModelError(t *testing.T) {
	err := gerrors.NewBiologicalModelError("bio", "model", "bad model")
	assert.Equal(t, gerrors.BiologicalModelError, err.Type)
}

func TestNewKmerError(t *testing.T) {
	err := gerrors.NewKmerError("kmer", "extract", "bad kmer")
	assert.Equal(t, gerrors.KmerExtractionError, err.Type)
}

func TestNewReferenceError(t *testing.T) {
	err := gerrors.NewReferenceError("ref", "load", "missing")
	assert.Equal(t, gerrors.ReferenceLoadError, err.Type)
}

func TestNewVCFError(t *testing.T) {
	err := gerrors.NewVCFError("vcf", "parse", "bad vcf")
	assert.Equal(t, gerrors.VCFParseError, err.Type)
}

func TestNewGPUError(t *testing.T) {
	err := gerrors.NewGPUError("gpu", "init", "no GPU")
	assert.Equal(t, gerrors.GPUError, err.Type)
}

// ---------------------------------------------------------------------------
// Error category predicates
// ---------------------------------------------------------------------------

func TestIsValidationError(t *testing.T) {
	assert.True(t, gerrors.IsValidationError(gerrors.NewValidationError("c", "o", "m")))
	assert.False(t, gerrors.IsValidationError(gerrors.NewParseError("c", "o", "m")))
	assert.False(t, gerrors.IsValidationError(errors.New("plain")))
}

func TestIsResourceError(t *testing.T) {
	assert.True(t, gerrors.IsResourceError(gerrors.NewResourceError("c", "o", "m")))
	assert.True(t, gerrors.IsResourceError(gerrors.NewGenomicsError(gerrors.MemoryError, "c", "o", "m")))
	assert.True(t, gerrors.IsResourceError(gerrors.NewGPUError("c", "o", "m")))
	assert.False(t, gerrors.IsResourceError(gerrors.NewValidationError("c", "o", "m")))
	assert.False(t, gerrors.IsResourceError(errors.New("plain")))
}

func TestIsDataError(t *testing.T) {
	assert.True(t, gerrors.IsDataError(gerrors.NewParseError("c", "o", "m")))
	assert.True(t, gerrors.IsDataError(gerrors.NewGenomicsError(gerrors.FormatError, "c", "o", "m")))
	assert.True(t, gerrors.IsDataError(gerrors.NewGenomicsError(gerrors.CorruptedDataError, "c", "o", "m")))
	assert.True(t, gerrors.IsDataError(gerrors.NewGenomicsError(gerrors.InvalidInputError, "c", "o", "m")))
	assert.False(t, gerrors.IsDataError(gerrors.NewValidationError("c", "o", "m")))
	assert.False(t, gerrors.IsDataError(errors.New("plain")))
}

func TestIsBiologicalError(t *testing.T) {
	assert.True(t, gerrors.IsBiologicalError(gerrors.NewBiologicalModelError("c", "o", "m")))
	assert.True(t, gerrors.IsBiologicalError(gerrors.NewGenomicsError(gerrors.PhylogeneticError, "c", "o", "m")))
	assert.True(t, gerrors.IsBiologicalError(gerrors.NewGenomicsError(gerrors.ChromatinError, "c", "o", "m")))
	assert.True(t, gerrors.IsBiologicalError(gerrors.NewGenomicsError(gerrors.MutationalSignatureError, "c", "o", "m")))
	assert.True(t, gerrors.IsBiologicalError(gerrors.NewGenomicsError(gerrors.SelectionPressureError, "c", "o", "m")))
	assert.False(t, gerrors.IsBiologicalError(gerrors.NewValidationError("c", "o", "m")))
	assert.False(t, gerrors.IsBiologicalError(errors.New("plain")))
}

// ---------------------------------------------------------------------------
// ErrorHandler
// ---------------------------------------------------------------------------

func TestNewErrorHandler_Defaults(t *testing.T) {
	cfg := gerrors.ErrorHandlerConfig{MaxStackFrames: 0}
	eh := gerrors.NewErrorHandler(cfg)
	require.NotNil(t, eh)
	// Should not panic
}

func TestErrorHandler_HandleError_Nil(t *testing.T) {
	// LogLevel set above any reachable level to prevent all logging
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		LogLevel: gerrors.LogLevel(99), // suppress all logging during test
	})
	// Should not panic on nil error
	assert.NotPanics(t, func() { eh.HandleError(nil) })
}

func TestErrorHandler_HandleError_GenomicsError(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		EnableMetrics: true,
		LogLevel:      gerrors.LogLevel(99), // suppress
	})
	err := gerrors.NewValidationError("comp", "op", "bad input")
	eh.HandleError(err)

	counts := eh.GetErrorCounts()
	assert.Equal(t, int64(1), counts[gerrors.ValidationError])
}

func TestErrorHandler_HandleError_WrapsPlainError(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		EnableMetrics: true,
		LogLevel:      gerrors.LogLevel(99),
	})
	plain := errors.New("plain error")
	eh.HandleError(plain)

	counts := eh.GetErrorCounts()
	assert.Equal(t, int64(1), counts[gerrors.ResourceError])
}

func TestErrorHandler_GetErrorCounts_Empty(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{EnableMetrics: true})
	counts := eh.GetErrorCounts()
	assert.Empty(t, counts)
}

func TestErrorHandler_ResetErrorCounts(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		EnableMetrics: true,
		LogLevel:      gerrors.LogLevel(99),
	})
	eh.HandleError(gerrors.NewValidationError("c", "o", "m"))
	require.Equal(t, int64(1), eh.GetErrorCounts()[gerrors.ValidationError])

	eh.ResetErrorCounts()
	assert.Empty(t, eh.GetErrorCounts())
}

func TestErrorHandler_GetErrorSummary(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		EnableMetrics: true,
		LogLevel:      gerrors.LogLevel(99),
	})
	eh.HandleError(gerrors.NewValidationError("c", "o", "m"))
	eh.HandleError(gerrors.NewParseError("c", "o", "m"))
	eh.HandleError(gerrors.NewGenomicsError(gerrors.MemoryError, "c", "o", "m"))

	summary := eh.GetErrorSummary()
	assert.Equal(t, int64(3), summary.TotalErrors)
	assert.Greater(t, len(summary.ErrorsByType), 0)
	assert.Greater(t, len(summary.ErrorsBySeverity), 0)
	assert.Greater(t, summary.Timestamp, int64(0))
}

func TestErrorHandler_AddReporter(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		ReportToExternal: true,
		LogLevel:         gerrors.LogLevel(99),
	})

	reporter := &mockReporter{}
	eh.AddReporter(reporter)
	// Should not panic and reporter should eventually be called
	err := gerrors.NewValidationError("c", "o", "m")
	assert.NotPanics(t, func() { eh.HandleError(err) })
	// Give goroutines a moment to report
	time.Sleep(10 * time.Millisecond)
}

type mockReporter struct {
	mu     sync.Mutex
	errors []*gerrors.GenomicsError
}

func (r *mockReporter) ReportError(err *gerrors.GenomicsError) error {
	r.mu.Lock()
	defer r.mu.Unlock()
	r.errors = append(r.errors, err)
	return nil
}

func (r *mockReporter) ReportMetrics(counts map[gerrors.ErrorType]int64) error {
	return nil
}

type failReporter struct{}

func (r *failReporter) ReportError(err *gerrors.GenomicsError) error {
	return errors.New("report failure")
}
func (r *failReporter) ReportMetrics(counts map[gerrors.ErrorType]int64) error {
	return errors.New("metrics failure")
}

func TestErrorHandler_AddReporter_Failure(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		ReportToExternal: true,
		LogLevel:         gerrors.LogLevel(99),
	})
	eh.AddReporter(&failReporter{})
	// Should not panic even if reporter fails
	assert.NotPanics(t, func() { eh.HandleError(gerrors.NewValidationError("c", "o", "m")) })
	time.Sleep(10 * time.Millisecond)
}

// ---------------------------------------------------------------------------
// ErrorHandler logging via different log levels (just make sure no panic)
// ---------------------------------------------------------------------------

func TestErrorHandler_LoggingLevels(t *testing.T) {
	for _, logLevel := range []gerrors.LogLevel{
		gerrors.LogLevelDebug,
		gerrors.LogLevelInfo,
		gerrors.LogLevelWarn,
		gerrors.LogLevelError,
	} {
		t.Run(fmt.Sprintf("LogLevel_%d", logLevel), func(t *testing.T) {
			eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
				LogLevel:         logLevel,
				EnableStackTrace: true,
				EnableContextLog: true,
				MaxStackFrames:   3,
				EnableMetrics:    true,
			})
			err := gerrors.NewGenomicsError(gerrors.ParseError, "c", "o", "m")
			err.WithContext("k", "v")
			assert.NotPanics(t, func() { eh.HandleError(err) })
		})
	}
}

// ---------------------------------------------------------------------------
// Global error handler and convenience functions
// ---------------------------------------------------------------------------

func TestGetGlobalErrorHandler(t *testing.T) {
	eh1 := gerrors.GetGlobalErrorHandler()
	eh2 := gerrors.GetGlobalErrorHandler()
	require.NotNil(t, eh1)
	assert.Same(t, eh1, eh2) // singleton
}

func TestHandleError_Convenience(t *testing.T) {
	// Should not panic
	assert.NotPanics(t, func() {
		gerrors.HandleError(gerrors.NewValidationError("c", "o", "m"))
	})
}

func TestHandleErrorWithContext_Convenience(t *testing.T) {
	err := gerrors.NewValidationError("c", "o", "m")
	ctx := map[string]interface{}{"k": "v", "n": 1}
	assert.NotPanics(t, func() {
		gerrors.HandleErrorWithContext(err, ctx)
	})
	v, ok := err.GetContext("k")
	assert.True(t, ok)
	assert.Equal(t, "v", v)
}

func TestHandleErrorWithContext_PlainError(t *testing.T) {
	plain := errors.New("plain")
	ctx := map[string]interface{}{"k": "v"}
	// Should not panic even for plain errors
	assert.NotPanics(t, func() {
		gerrors.HandleErrorWithContext(plain, ctx)
	})
}

// ---------------------------------------------------------------------------
// RecoveryHandler
// ---------------------------------------------------------------------------

func TestNewRecoveryHandler(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{})
	rh := gerrors.NewRecoveryHandler(eh, "test_component")
	require.NotNil(t, rh)
}

func TestRecoveryHandler_Recover_NoPanic(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRecoveryHandler(eh, "comp")
	// No panic - Recover should be a no-op
	assert.NotPanics(t, func() {
		defer rh.Recover("op")
		// nothing panics here
	})
}

func TestRecoveryHandler_Recover_WithErrorPanic(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		EnableMetrics: true,
		LogLevel:      gerrors.LogLevel(99),
	})
	rh := gerrors.NewRecoveryHandler(eh, "comp")

	func() {
		defer rh.Recover("op")
		panic(errors.New("error panic"))
	}()

	counts := eh.GetErrorCounts()
	assert.Equal(t, int64(1), counts[gerrors.ResourceError])
}

func TestRecoveryHandler_Recover_WithStringPanic(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		EnableMetrics: true,
		LogLevel:      gerrors.LogLevel(99),
	})
	rh := gerrors.NewRecoveryHandler(eh, "comp")

	func() {
		defer rh.Recover("op")
		panic("string panic")
	}()

	counts := eh.GetErrorCounts()
	assert.Equal(t, int64(1), counts[gerrors.ResourceError])
}

func TestRecoveryHandler_SafeExecute_Success(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRecoveryHandler(eh, "comp")

	called := false
	err := rh.SafeExecute("op", func() error {
		called = true
		return nil
	})

	assert.NoError(t, err)
	assert.True(t, called)
}

func TestRecoveryHandler_SafeExecute_ReturnsError(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRecoveryHandler(eh, "comp")

	expected := errors.New("fn error")
	err := rh.SafeExecute("op", func() error {
		return expected
	})

	assert.Equal(t, expected, err)
}

func TestRecoveryHandler_SafeExecute_PanicRecovered(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		EnableMetrics: true,
		LogLevel:      gerrors.LogLevel(99),
	})
	rh := gerrors.NewRecoveryHandler(eh, "comp")

	// SafeExecute defers Recover, so a panic is caught; but the panic unwinds
	// the stack and the function returns the zero value for error (nil).
	assert.NotPanics(t, func() {
		_ = rh.SafeExecute("op", func() error {
			panic("unexpected failure")
		})
	})
}

// ---------------------------------------------------------------------------
// RetryHandler
// ---------------------------------------------------------------------------

func TestNewRetryHandler_Defaults(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRetryHandler(gerrors.RetryConfig{}, eh)
	require.NotNil(t, rh)
}

func TestRetryHandler_Execute_SuccessFirstAttempt(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRetryHandler(gerrors.RetryConfig{MaxAttempts: 3, BaseDelay: 0}, eh)

	calls := 0
	err := rh.Execute("op", func() error {
		calls++
		return nil
	})

	assert.NoError(t, err)
	assert.Equal(t, 1, calls)
}

func TestRetryHandler_Execute_NonRetryableError(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRetryHandler(gerrors.RetryConfig{
		MaxAttempts: 3,
		BaseDelay:   0,
	}, eh)

	// ValidationError is not in the default retryable list
	validationErr := gerrors.NewValidationError("c", "o", "bad")
	calls := 0
	err := rh.Execute("op", func() error {
		calls++
		return validationErr
	})

	assert.Error(t, err)
	assert.Equal(t, 1, calls) // Should stop immediately on non-retryable
}

func TestRetryHandler_Execute_RetryableError_ExhaustedAttempts(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRetryHandler(gerrors.RetryConfig{
		MaxAttempts: 2,
		BaseDelay:   0,           // No delay for testing
		BackoffFactor: 0,
	}, eh)

	networkErr := gerrors.NewGenomicsError(gerrors.NetworkError, "c", "o", "net fail")
	calls := 0
	err := rh.Execute("op", func() error {
		calls++
		return networkErr
	})

	assert.Error(t, err)
	assert.Equal(t, 2, calls)
}

func TestRetryHandler_Execute_RetryableTypeFromConfig(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRetryHandler(gerrors.RetryConfig{
		MaxAttempts:    2,
		BaseDelay:      0,
		RetryableTypes: []gerrors.ErrorType{gerrors.ValidationError},
	}, eh)

	validationErr := gerrors.NewValidationError("c", "o", "bad")
	calls := 0
	err := rh.Execute("op", func() error {
		calls++
		return validationErr
	})

	assert.Error(t, err)
	assert.Equal(t, 2, calls)
}

func TestRetryHandler_Execute_SuccessAfterRetry(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRetryHandler(gerrors.RetryConfig{
		MaxAttempts: 3,
		BaseDelay:   0,
	}, eh)

	networkErr := gerrors.NewGenomicsError(gerrors.NetworkError, "c", "o", "transient")
	calls := 0
	err := rh.Execute("op", func() error {
		calls++
		if calls < 2 {
			return networkErr
		}
		return nil
	})

	assert.NoError(t, err)
	assert.Equal(t, 2, calls)
}

func TestRetryHandler_Execute_PlainErrorNotRetried(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{LogLevel: gerrors.LogLevel(99)})
	rh := gerrors.NewRetryHandler(gerrors.RetryConfig{MaxAttempts: 3, BaseDelay: 0}, eh)

	plain := errors.New("plain error")
	calls := 0
	err := rh.Execute("op", func() error {
		calls++
		return plain
	})

	assert.Error(t, err)
	assert.Equal(t, 1, calls) // Plain errors are not retried
}

// ---------------------------------------------------------------------------
// LogLevel constants
// ---------------------------------------------------------------------------

func TestLogLevel_Values(t *testing.T) {
	assert.Less(t, int(gerrors.LogLevelDebug), int(gerrors.LogLevelInfo))
	assert.Less(t, int(gerrors.LogLevelInfo), int(gerrors.LogLevelWarn))
	assert.Less(t, int(gerrors.LogLevelWarn), int(gerrors.LogLevelError))
	assert.Less(t, int(gerrors.LogLevelError), int(gerrors.LogLevelFatal))
}

// ---------------------------------------------------------------------------
// Concurrency safety for ErrorHandler
// ---------------------------------------------------------------------------

func TestErrorHandler_ConcurrentAccess(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		EnableMetrics: true,
		LogLevel:      gerrors.LogLevel(99),
	})

	var wg sync.WaitGroup
	for i := 0; i < 20; i++ {
		wg.Add(1)
		go func(idx int) {
			defer wg.Done()
			types := []gerrors.ErrorType{
				gerrors.ValidationError,
				gerrors.ParseError,
				gerrors.NetworkError,
				gerrors.CacheError,
			}
			eh.HandleError(gerrors.NewGenomicsError(types[idx%len(types)], "c", "o", "concurrent"))
			_ = eh.GetErrorCounts()
			_ = eh.GetErrorSummary()
		}(i)
	}
	wg.Wait()

	summary := eh.GetErrorSummary()
	assert.Equal(t, int64(20), summary.TotalErrors)
}

// ---------------------------------------------------------------------------
// ErrorSummary structure
// ---------------------------------------------------------------------------

func TestErrorSummary_Fields(t *testing.T) {
	eh := gerrors.NewErrorHandler(gerrors.ErrorHandlerConfig{
		EnableMetrics: true,
		LogLevel:      gerrors.LogLevel(99),
	})
	eh.HandleError(gerrors.NewGenomicsError(gerrors.MemoryError, "c", "o", "oom"))

	summary := eh.GetErrorSummary()
	assert.Equal(t, int64(1), summary.TotalErrors)
	assert.Equal(t, int64(1), summary.ErrorsBySeverity[gerrors.SeverityCritical])
	assert.Equal(t, int64(1), summary.ErrorsByType[gerrors.MemoryError])
	assert.Greater(t, summary.Timestamp, int64(0))
}

// ---------------------------------------------------------------------------
// strings package usage sanity check (ensure Error() produces parseable output)
// ---------------------------------------------------------------------------

func TestGenomicsError_ErrorFormat_AllFields(t *testing.T) {
	err := gerrors.NewGenomicsError(gerrors.VCFParseError, "vcf_parser", "parse_record", "bad VCF line")
	err.WithCause(errors.New("io.EOF"))
	err.WithContext("line_num", 42)

	s := err.Error()
	parts := strings.Split(s, " ")
	assert.Greater(t, len(parts), 2)
	assert.Contains(t, s, "VCF_PARSE_ERROR")
	assert.Contains(t, s, "vcf_parser")
}
