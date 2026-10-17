package errors

import (
	"fmt"
	"runtime"
	"strings"
	"time"
)

// GenomicsError represents a standardized error type for the genomics authenticity engine
type GenomicsError struct {
	Type        ErrorType
	Operation   string
	Component   string
	Message     string
	Cause       error
	Context     map[string]interface{}
	StackTrace  []string
	Timestamp   int64
}

// ErrorType represents different categories of errors in the genomics system
type ErrorType int

const (
	// Core system errors
	ValidationError ErrorType = iota
	ConfigurationError
	InitializationError
	ResourceError
	
	// Data processing errors
	ParseError
	FormatError
	CorruptedDataError
	InvalidInputError
	
	// Performance/Resource errors
	MemoryError
	GPUError
	CacheError
	TimeoutError
	
	// Biological analysis errors
	BiologicalModelError
	PhylogeneticError
	ChromatinError
	MutationalSignatureError
	SelectionPressureError
	
	// K-mer analysis errors
	KmerExtractionError
	KmerValidationError
	PackedArrayError
	
	// Reference genome errors
	ReferenceLoadError
	IndexError
	SequenceNotFoundError
	
	// VCF processing errors
	VCFParseError
	VCFFormatError
	VariantProcessingError
	
	// Integration errors
	ExternalToolError
	DatabaseError
	NetworkError
	
	// Authentication/Security errors
	AuthenticityValidationError
	SecurityError
	AccessDeniedError
)

// Error implements the error interface
func (ge *GenomicsError) Error() string {
	var parts []string
	
	// Add error type and component
	parts = append(parts, fmt.Sprintf("[%s:%s]", ge.Type.String(), ge.Component))
	
	// Add operation if specified
	if ge.Operation != "" {
		parts = append(parts, fmt.Sprintf("operation=%s", ge.Operation))
	}
	
	// Add main message
	parts = append(parts, ge.Message)
	
	// Add cause if present
	if ge.Cause != nil {
		parts = append(parts, fmt.Sprintf("cause: %v", ge.Cause))
	}
	
	return strings.Join(parts, " ")
}

// String returns a string representation of the error type
func (et ErrorType) String() string {
	errorTypeNames := map[ErrorType]string{
		ValidationError:            "VALIDATION_ERROR",
		ConfigurationError:         "CONFIGURATION_ERROR", 
		InitializationError:        "INITIALIZATION_ERROR",
		ResourceError:              "RESOURCE_ERROR",
		ParseError:                 "PARSE_ERROR",
		FormatError:                "FORMAT_ERROR",
		CorruptedDataError:         "CORRUPTED_DATA_ERROR",
		InvalidInputError:          "INVALID_INPUT_ERROR",
		MemoryError:                "MEMORY_ERROR",
		GPUError:                   "GPU_ERROR",
		CacheError:                 "CACHE_ERROR",
		TimeoutError:               "TIMEOUT_ERROR",
		BiologicalModelError:       "BIOLOGICAL_MODEL_ERROR",
		PhylogeneticError:          "PHYLOGENETIC_ERROR",
		ChromatinError:             "CHROMATIN_ERROR",
		MutationalSignatureError:   "MUTATIONAL_SIGNATURE_ERROR",
		SelectionPressureError:     "SELECTION_PRESSURE_ERROR",
		KmerExtractionError:        "KMER_EXTRACTION_ERROR",
		KmerValidationError:        "KMER_VALIDATION_ERROR",
		PackedArrayError:           "PACKED_ARRAY_ERROR",
		ReferenceLoadError:         "REFERENCE_LOAD_ERROR",
		IndexError:                 "INDEX_ERROR",
		SequenceNotFoundError:      "SEQUENCE_NOT_FOUND_ERROR",
		VCFParseError:              "VCF_PARSE_ERROR",
		VCFFormatError:             "VCF_FORMAT_ERROR",
		VariantProcessingError:     "VARIANT_PROCESSING_ERROR",
		ExternalToolError:          "EXTERNAL_TOOL_ERROR",
		DatabaseError:              "DATABASE_ERROR",
		NetworkError:               "NETWORK_ERROR",
		AuthenticityValidationError: "AUTHENTICITY_VALIDATION_ERROR",
		SecurityError:              "SECURITY_ERROR",
		AccessDeniedError:          "ACCESS_DENIED_ERROR",
	}
	
	if name, exists := errorTypeNames[et]; exists {
		return name
	}
	return "UNKNOWN_ERROR"
}

// NewGenomicsError creates a new standardized genomics error
func NewGenomicsError(errorType ErrorType, component, operation, message string) *GenomicsError {
	return &GenomicsError{
		Type:       errorType,
		Component:  component,
		Operation:  operation,
		Message:    message,
		Context:    make(map[string]interface{}),
		StackTrace: captureStackTrace(),
		Timestamp:  getCurrentTimestamp(),
	}
}

// WrapError wraps an existing error with genomics-specific context
func WrapError(err error, errorType ErrorType, component, operation, message string) *GenomicsError {
	if err == nil {
		return nil
	}
	
	// If it's already a GenomicsError, add context
	if ge, ok := err.(*GenomicsError); ok {
		ge.Context["wrapped_at"] = fmt.Sprintf("%s:%s", component, operation)
		return ge
	}
	
	return &GenomicsError{
		Type:       errorType,
		Component:  component,
		Operation:  operation,
		Message:    message,
		Cause:      err,
		Context:    make(map[string]interface{}),
		StackTrace: captureStackTrace(),
		Timestamp:  getCurrentTimestamp(),
	}
}

// WithContext adds contextual information to the error
func (ge *GenomicsError) WithContext(key string, value interface{}) *GenomicsError {
	if ge.Context == nil {
		ge.Context = make(map[string]interface{})
	}
	ge.Context[key] = value
	return ge
}

// WithCause sets the underlying cause of the error
func (ge *GenomicsError) WithCause(cause error) *GenomicsError {
	ge.Cause = cause
	return ge
}

// IsType checks if the error is of a specific type
func (ge *GenomicsError) IsType(errorType ErrorType) bool {
	return ge.Type == errorType
}

// IsComponent checks if the error originated from a specific component
func (ge *GenomicsError) IsComponent(component string) bool {
	return ge.Component == component
}

// GetContext retrieves contextual information from the error
func (ge *GenomicsError) GetContext(key string) (interface{}, bool) {
	if ge.Context == nil {
		return nil, false
	}
	value, exists := ge.Context[key]
	return value, exists
}

// GetStackTrace returns the stack trace as a formatted string
func (ge *GenomicsError) GetStackTrace() string {
	return strings.Join(ge.StackTrace, "\n")
}

// captureStackTrace captures the current stack trace
func captureStackTrace() []string {
	var stackTrace []string
	
	// Skip the first few frames (this function and error creation)
	for i := 2; i < 10; i++ {
		pc, file, line, ok := runtime.Caller(i)
		if !ok {
			break
		}
		
		// Get function name
		fn := runtime.FuncForPC(pc)
		var funcName string
		if fn != nil {
			funcName = fn.Name()
		} else {
			funcName = "unknown"
		}
		
		// Simplify file path
		if idx := strings.LastIndex(file, "/"); idx >= 0 {
			file = file[idx+1:]
		}
		
		stackTrace = append(stackTrace, fmt.Sprintf("%s:%d %s", file, line, funcName))
	}
	
	return stackTrace
}

// getCurrentTimestamp returns the current Unix timestamp in nanoseconds
func getCurrentTimestamp() int64 {
	return time.Now().UnixNano()
}

// Convenience functions for common error patterns

// NewValidationError creates a validation error
func NewValidationError(component, operation, message string) *GenomicsError {
	return NewGenomicsError(ValidationError, component, operation, message)
}

// NewParseError creates a parse error
func NewParseError(component, operation, message string) *GenomicsError {
	return NewGenomicsError(ParseError, component, operation, message)
}

// NewResourceError creates a resource error
func NewResourceError(component, operation, message string) *GenomicsError {
	return NewGenomicsError(ResourceError, component, operation, message)
}

// NewBiologicalModelError creates a biological model error
func NewBiologicalModelError(component, operation, message string) *GenomicsError {
	return NewGenomicsError(BiologicalModelError, component, operation, message)
}

// NewKmerError creates a k-mer processing error
func NewKmerError(component, operation, message string) *GenomicsError {
	return NewGenomicsError(KmerExtractionError, component, operation, message)
}

// NewReferenceError creates a reference genome error
func NewReferenceError(component, operation, message string) *GenomicsError {
	return NewGenomicsError(ReferenceLoadError, component, operation, message)
}

// NewVCFError creates a VCF processing error
func NewVCFError(component, operation, message string) *GenomicsError {
	return NewGenomicsError(VCFParseError, component, operation, message)
}

// NewGPUError creates a GPU processing error
func NewGPUError(component, operation, message string) *GenomicsError {
	return NewGenomicsError(GPUError, component, operation, message)
}

// Error checking utilities

// IsValidationError checks if an error is a validation error
func IsValidationError(err error) bool {
	if ge, ok := err.(*GenomicsError); ok {
		return ge.IsType(ValidationError)
	}
	return false
}

// IsResourceError checks if an error is a resource error
func IsResourceError(err error) bool {
	if ge, ok := err.(*GenomicsError); ok {
		return ge.IsType(ResourceError) || ge.IsType(MemoryError) || ge.IsType(GPUError)
	}
	return false
}

// IsDataError checks if an error is related to data processing
func IsDataError(err error) bool {
	if ge, ok := err.(*GenomicsError); ok {
		return ge.IsType(ParseError) || ge.IsType(FormatError) || 
			   ge.IsType(CorruptedDataError) || ge.IsType(InvalidInputError)
	}
	return false
}

// IsBiologicalError checks if an error is related to biological modeling
func IsBiologicalError(err error) bool {
	if ge, ok := err.(*GenomicsError); ok {
		return ge.IsType(BiologicalModelError) || ge.IsType(PhylogeneticError) ||
			   ge.IsType(ChromatinError) || ge.IsType(MutationalSignatureError) ||
			   ge.IsType(SelectionPressureError)
	}
	return false
}

// ErrorSeverity represents the severity level of an error
type ErrorSeverity int

const (
	SeverityLow ErrorSeverity = iota
	SeverityMedium
	SeverityHigh
	SeverityCritical
)

// GetSeverity returns the severity level of an error type
func (et ErrorType) GetSeverity() ErrorSeverity {
	severityMap := map[ErrorType]ErrorSeverity{
		// Critical errors - system cannot continue
		MemoryError:                SeverityCritical,
		CorruptedDataError:         SeverityCritical,
		SecurityError:              SeverityCritical,
		AccessDeniedError:          SeverityCritical,
		
		// High severity - major functionality impacted
		InitializationError:        SeverityHigh,
		GPUError:                   SeverityHigh,
		ReferenceLoadError:         SeverityHigh,
		DatabaseError:              SeverityHigh,
		
		// Medium severity - some functionality impacted
		ValidationError:            SeverityMedium,
		ConfigurationError:         SeverityMedium,
		BiologicalModelError:       SeverityMedium,
		VCFParseError:              SeverityMedium,
		ExternalToolError:          SeverityMedium,
		
		// Low severity - minor issues, system can continue
		ParseError:                 SeverityLow,
		FormatError:                SeverityLow,
		CacheError:                 SeverityLow,
		TimeoutError:               SeverityLow,
	}
	
	if severity, exists := severityMap[et]; exists {
		return severity
	}
	return SeverityMedium // Default to medium
}

// String returns a string representation of the severity
func (es ErrorSeverity) String() string {
	switch es {
	case SeverityLow:
		return "LOW"
	case SeverityMedium:
		return "MEDIUM"
	case SeverityHigh:
		return "HIGH"
	case SeverityCritical:
		return "CRITICAL"
	default:
		return "UNKNOWN"
	}
}