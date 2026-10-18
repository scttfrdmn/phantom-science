package fingerprinting

import (
	"errors"
	"fmt"
	"strings"
)

// Error types for different categories of fingerprinting failures
var (
	// Input validation errors
	ErrEmptyInput          = errors.New("empty input data")
	ErrInvalidFormat       = errors.New("invalid data format")
	ErrInvalidFingerprint  = errors.New("invalid fingerprint")
	ErrInsufficientData    = errors.New("insufficient data for embedding")
	
	// System errors
	ErrBridgeNotInitialized = errors.New("stegano bridge not initialized")
	ErrLibraryNotAvailable  = errors.New("stegano library not available")
	ErrGPUNotAvailable     = errors.New("GPU acceleration not available")
	
	// Operation errors
	ErrEmbeddingFailed     = errors.New("embedding operation failed")
	ErrExtractionFailed    = errors.New("extraction operation failed")
	ErrCapacityExceeded    = errors.New("fingerprint exceeds embedding capacity")
	ErrCorruptedData       = errors.New("data appears to be corrupted")
	ErrKeyRequired         = errors.New("fingerprint is encrypted: a key is required")
	ErrWrongKey            = errors.New("fingerprint cannot be decrypted with this key")
	
	// Performance errors
	ErrTimeout             = errors.New("operation timed out")
	ErrMemoryExhausted     = errors.New("insufficient memory for operation")
	ErrGPUError           = errors.New("GPU operation failed")
)

// FingerprintingError provides detailed error information with context
type FingerprintingError struct {
	Op          string                 // Operation that failed
	Err         error                  // Underlying error
	Code        ErrorCode              // Error code for programmatic handling
	Context     map[string]interface{} // Additional context
	Recoverable bool                   // Whether the error can be recovered from
}

// ErrorCode represents different types of fingerprinting errors
type ErrorCode int

const (
	ErrorCodeUnknown ErrorCode = iota
	ErrorCodeValidation
	ErrorCodeCapacity
	ErrorCodeSystem
	ErrorCodeGPU
	ErrorCodeTimeout
	ErrorCodeCorruption
	ErrorCodePermission
)

func (e *FingerprintingError) Error() string {
	var parts []string
	
	if e.Op != "" {
		parts = append(parts, fmt.Sprintf("operation %s", e.Op))
	}
	
	if e.Err != nil {
		parts = append(parts, e.Err.Error())
	}
	
	if len(e.Context) > 0 {
		var contextParts []string
		for key, value := range e.Context {
			contextParts = append(contextParts, fmt.Sprintf("%s=%v", key, value))
		}
		parts = append(parts, fmt.Sprintf("context: %s", strings.Join(contextParts, ", ")))
	}
	
	return strings.Join(parts, ": ")
}

func (e *FingerprintingError) Unwrap() error {
	return e.Err
}

func (e *FingerprintingError) Is(target error) bool {
	if target == nil {
		return false
	}
	return errors.Is(e.Err, target)
}

// NewFingerprintingError creates a new fingerprinting error with context
func NewFingerprintingError(op string, err error, code ErrorCode) *FingerprintingError {
	return &FingerprintingError{
		Op:          op,
		Err:         err,
		Code:        code,
		Context:     make(map[string]interface{}),
		Recoverable: isRecoverableError(code),
	}
}

// WithContext adds context information to the error
func (e *FingerprintingError) WithContext(key string, value interface{}) *FingerprintingError {
	if e.Context == nil {
		e.Context = make(map[string]interface{})
	}
	e.Context[key] = value
	return e
}

// WithRecoverable sets whether the error is recoverable
func (e *FingerprintingError) WithRecoverable(recoverable bool) *FingerprintingError {
	e.Recoverable = recoverable
	return e
}

// isRecoverableError determines if an error type is typically recoverable
func isRecoverableError(code ErrorCode) bool {
	switch code {
	case ErrorCodeGPU, ErrorCodeTimeout:
		return true // Can fallback to CPU or retry
	case ErrorCodeValidation, ErrorCodeCorruption:
		return false // Data issues are not recoverable
	case ErrorCodeCapacity:
		return true // Can try different strategy or smaller payload
	case ErrorCodeSystem:
		return false // System-level issues usually require intervention
	default:
		return false
	}
}

// ValidateInput performs comprehensive input validation
func ValidateInput(data []byte, fingerprint string, operation string) error {
	if len(data) == 0 {
		return NewFingerprintingError(operation, ErrEmptyInput, ErrorCodeValidation).
			WithContext("data_length", 0)
	}

	if err := validateFingerprint(fingerprint, operation); err != nil {
		return err
	}

	// Validate FASTQ format structure
	if !isValidFASTQData(data) {
		return NewFingerprintingError(operation, ErrInvalidFormat, ErrorCodeValidation).
			WithContext("data_length", len(data)).
			WithContext("format", "FASTQ")
	}

	return nil
}

// validateFingerprint checks that fingerprint is non-empty, at most
// MaxFingerprintLength bytes, and printable ASCII.
func validateFingerprint(fingerprint string, operation string) error {
	if fingerprint == "" {
		return NewFingerprintingError(operation, ErrInvalidFingerprint, ErrorCodeValidation).
			WithContext("fingerprint_length", 0)
	}
	if len(fingerprint) > MaxFingerprintLength {
		return NewFingerprintingError(operation, ErrInvalidFingerprint, ErrorCodeValidation).
			WithContext("fingerprint_length", len(fingerprint)).
			WithContext("max_length", MaxFingerprintLength)
	}
	if containsInvalidCharacters(fingerprint) {
		return NewFingerprintingError(operation, ErrInvalidFingerprint, ErrorCodeValidation).
			WithContext("fingerprint", fingerprint)
	}
	return nil
}

// EstimateCapacity validates that the fingerprint can fit in the available data
func EstimateCapacity(data []byte, fingerprint string, strategy string) error {
	dataLength := len(data)
	fingerprintLength := len(fingerprint)
	
	// Rough capacity estimation based on FASTQ structure
	// Assumes approximately 40% of data is quality scores
	estimatedQualityBytes := int(float64(dataLength) * 0.4)
	
	// Conservative capacity estimation (1 bit per quality score)
	conservativeCapacity := estimatedQualityBytes
	
	if fingerprintLength > conservativeCapacity {
		return NewFingerprintingError("capacity_check", ErrCapacityExceeded, ErrorCodeCapacity).
			WithContext("fingerprint_length", fingerprintLength).
			WithContext("estimated_capacity", conservativeCapacity).
			WithContext("data_length", dataLength).
			WithContext("strategy", strategy)
	}
	
	return nil
}

// Constants for validation
const (
	MaxFingerprintLength = 1024 // Maximum allowed fingerprint length
	MinDataLength        = 20   // Minimum data length for embedding
)

// isValidFASTQData performs basic FASTQ format validation
func isValidFASTQData(data []byte) bool {
	if len(data) < MinDataLength {
		return false
	}
	
	// Check for FASTQ structure (starts with @, has + separator)
	dataStr := string(data)
	lines := strings.Split(dataStr, "\n")
	
	if len(lines) < 4 {
		return false
	}
	
	// Basic FASTQ validation: first line starts with @, third line starts with +
	firstLine := strings.TrimSpace(lines[0])
	if !strings.HasPrefix(firstLine, "@") {
		return false
	}
	
	// Look for quality separator
	hasQualitySeparator := false
	for _, line := range lines {
		if strings.HasPrefix(strings.TrimSpace(line), "+") {
			hasQualitySeparator = true
			break
		}
	}
	
	return hasQualitySeparator
}

// containsInvalidCharacters checks for characters that might cause issues
func containsInvalidCharacters(fingerprint string) bool {
	// Check for null bytes, control characters, or other problematic content
	for _, char := range fingerprint {
		if char < 32 || char > 126 {
			// Non-printable ASCII characters
			return true
		}
	}
	return false
}

// RecoverableOperation wraps operations with automatic recovery for recoverable errors
func RecoverableOperation[T any](operation func() (T, error), fallback func() (T, error), maxRetries int) (T, error) {
	var lastErr error
	var zero T
	
	for attempt := 0; attempt <= maxRetries; attempt++ {
		result, err := operation()
		if err == nil {
			return result, nil
		}
		
		lastErr = err
		
		// Check if error is recoverable
		var fpErr *FingerprintingError
		if errors.As(err, &fpErr) && fpErr.Recoverable && attempt < maxRetries {
			continue // Retry the operation
		}
		
		// If not recoverable or max retries reached, try fallback
		if fallback != nil {
			result, fallbackErr := fallback()
			if fallbackErr == nil {
				return result, nil
			}
			// If fallback also fails, return original error
		}
		
		break
	}
	
	return zero, lastErr
}