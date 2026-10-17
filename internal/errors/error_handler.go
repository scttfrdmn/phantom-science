package errors

import (
	"fmt"
	"log"
	"os"
	"sync"
	"time"
)

// ErrorHandler provides centralized error handling and logging
type ErrorHandler struct {
	logger          *log.Logger
	errorCounts     map[ErrorType]int64
	errorCountMutex sync.RWMutex
	
	// Configuration
	config ErrorHandlerConfig
	
	// Error reporting
	reporters []ErrorReporter
}

// ErrorHandlerConfig configures the error handler behavior
type ErrorHandlerConfig struct {
	EnableStackTrace   bool
	EnableContextLog   bool
	LogLevel          LogLevel
	MaxStackFrames    int
	EnableMetrics     bool
	ReportToExternal  bool
}

// LogLevel represents different logging levels
type LogLevel int

const (
	LogLevelDebug LogLevel = iota
	LogLevelInfo
	LogLevelWarn
	LogLevelError
	LogLevelFatal
)

// ErrorReporter interface for external error reporting
type ErrorReporter interface {
	ReportError(err *GenomicsError) error
	ReportMetrics(counts map[ErrorType]int64) error
}

// NewErrorHandler creates a new error handler
func NewErrorHandler(config ErrorHandlerConfig) *ErrorHandler {
	// Set defaults
	if config.MaxStackFrames <= 0 {
		config.MaxStackFrames = 10
	}
	
	logger := log.New(os.Stderr, "[GENOMICS_ERROR] ", log.LstdFlags|log.Lshortfile)
	
	return &ErrorHandler{
		logger:      logger,
		errorCounts: make(map[ErrorType]int64),
		config:      config,
		reporters:   make([]ErrorReporter, 0),
	}
}

// HandleError processes and logs an error according to configuration
func (eh *ErrorHandler) HandleError(err error) {
	if err == nil {
		return
	}
	
	// Convert to GenomicsError if not already
	var genomicsErr *GenomicsError
	if ge, ok := err.(*GenomicsError); ok {
		genomicsErr = ge
	} else {
		genomicsErr = WrapError(err, ResourceError, "error_handler", "handle_error", "wrapped external error")
	}
	
	// Update error counts
	if eh.config.EnableMetrics {
		eh.updateErrorCounts(genomicsErr.Type)
	}
	
	// Determine log level based on error severity
	severity := genomicsErr.Type.GetSeverity()
	logLevel := eh.severityToLogLevel(severity)
	
	// Log the error if it meets the minimum log level
	if logLevel >= eh.config.LogLevel {
		eh.logError(genomicsErr, logLevel)
	}
	
	// Report to external systems
	if eh.config.ReportToExternal {
		eh.reportToExternal(genomicsErr)
	}
}

// updateErrorCounts updates the internal error count metrics
func (eh *ErrorHandler) updateErrorCounts(errorType ErrorType) {
	eh.errorCountMutex.Lock()
	defer eh.errorCountMutex.Unlock()
	
	eh.errorCounts[errorType]++
}

// severityToLogLevel converts error severity to log level
func (eh *ErrorHandler) severityToLogLevel(severity ErrorSeverity) LogLevel {
	switch severity {
	case SeverityLow:
		return LogLevelInfo
	case SeverityMedium:
		return LogLevelWarn
	case SeverityHigh:
		return LogLevelError
	case SeverityCritical:
		return LogLevelFatal
	default:
		return LogLevelWarn
	}
}

// logError logs the error with appropriate formatting
func (eh *ErrorHandler) logError(err *GenomicsError, level LogLevel) {
	levelStr := eh.logLevelToString(level)
	
	// Basic error information
	logMsg := fmt.Sprintf("[%s] %s", levelStr, err.Error())
	
	// Add context if enabled
	if eh.config.EnableContextLog && len(err.Context) > 0 {
		logMsg += fmt.Sprintf(" | Context: %+v", err.Context)
	}
	
	// Add stack trace if enabled
	if eh.config.EnableStackTrace && len(err.StackTrace) > 0 {
		maxFrames := eh.config.MaxStackFrames
		if maxFrames > len(err.StackTrace) {
			maxFrames = len(err.StackTrace)
		}
		
		logMsg += "\nStack trace:"
		for i := 0; i < maxFrames; i++ {
			logMsg += fmt.Sprintf("\n  %s", err.StackTrace[i])
		}
	}
	
	// Log based on level
	switch level {
	case LogLevelFatal:
		eh.logger.Fatal(logMsg)
	case LogLevelError:
		eh.logger.Println(logMsg)
	case LogLevelWarn:
		eh.logger.Println(logMsg)
	case LogLevelInfo:
		eh.logger.Println(logMsg)
	case LogLevelDebug:
		eh.logger.Println(logMsg)
	}
}

// logLevelToString converts log level to string
func (eh *ErrorHandler) logLevelToString(level LogLevel) string {
	switch level {
	case LogLevelDebug:
		return "DEBUG"
	case LogLevelInfo:
		return "INFO"
	case LogLevelWarn:
		return "WARN"
	case LogLevelError:
		return "ERROR"
	case LogLevelFatal:
		return "FATAL"
	default:
		return "UNKNOWN"
	}
}

// reportToExternal reports errors to external systems
func (eh *ErrorHandler) reportToExternal(err *GenomicsError) {
	for _, reporter := range eh.reporters {
		go func(r ErrorReporter) {
			if reportErr := r.ReportError(err); reportErr != nil {
				// Log reporting failure (avoid recursion)
				eh.logger.Printf("Failed to report error to external system: %v", reportErr)
			}
		}(reporter)
	}
}

// AddReporter adds an external error reporter
func (eh *ErrorHandler) AddReporter(reporter ErrorReporter) {
	eh.reporters = append(eh.reporters, reporter)
}

// GetErrorCounts returns the current error counts
func (eh *ErrorHandler) GetErrorCounts() map[ErrorType]int64 {
	eh.errorCountMutex.RLock()
	defer eh.errorCountMutex.RUnlock()
	
	// Return a copy to prevent concurrent access issues
	counts := make(map[ErrorType]int64)
	for errorType, count := range eh.errorCounts {
		counts[errorType] = count
	}
	
	return counts
}

// ResetErrorCounts resets all error counts to zero
func (eh *ErrorHandler) ResetErrorCounts() {
	eh.errorCountMutex.Lock()
	defer eh.errorCountMutex.Unlock()
	
	eh.errorCounts = make(map[ErrorType]int64)
}

// GetErrorSummary returns a summary of error statistics
func (eh *ErrorHandler) GetErrorSummary() ErrorSummary {
	eh.errorCountMutex.RLock()
	defer eh.errorCountMutex.RUnlock()
	
	summary := ErrorSummary{
		TotalErrors:      0,
		ErrorsBySeverity: make(map[ErrorSeverity]int64),
		ErrorsByType:     make(map[ErrorType]int64),
		Timestamp:        time.Now().Unix(),
	}
	
	for errorType, count := range eh.errorCounts {
		summary.TotalErrors += count
		summary.ErrorsByType[errorType] = count
		
		severity := errorType.GetSeverity()
		summary.ErrorsBySeverity[severity] += count
	}
	
	return summary
}

// ErrorSummary provides a summary of error statistics
type ErrorSummary struct {
	TotalErrors      int64
	ErrorsBySeverity map[ErrorSeverity]int64
	ErrorsByType     map[ErrorType]int64
	Timestamp        int64
}

// RecoveryHandler provides panic recovery with error reporting
type RecoveryHandler struct {
	errorHandler *ErrorHandler
	component    string
}

// NewRecoveryHandler creates a new recovery handler
func NewRecoveryHandler(errorHandler *ErrorHandler, component string) *RecoveryHandler {
	return &RecoveryHandler{
		errorHandler: errorHandler,
		component:    component,
	}
}

// Recover recovers from panics and reports them as errors
func (rh *RecoveryHandler) Recover(operation string) {
	if r := recover(); r != nil {
		// Convert panic to error
		var err error
		if e, ok := r.(error); ok {
			err = e
		} else {
			err = fmt.Errorf("panic: %v", r)
		}
		
		// Wrap as genomics error
		genomicsErr := WrapError(err, ResourceError, rh.component, operation, "recovered from panic")
		genomicsErr.WithContext("panic_value", r)
		
		// Report the error
		rh.errorHandler.HandleError(genomicsErr)
	}
}

// SafeExecute executes a function with panic recovery
func (rh *RecoveryHandler) SafeExecute(operation string, fn func() error) error {
	defer rh.Recover(operation)
	
	return fn()
}

// Retry mechanism for transient errors
type RetryConfig struct {
	MaxAttempts    int
	BaseDelay      time.Duration
	MaxDelay       time.Duration
	BackoffFactor  float64
	RetryableTypes []ErrorType
}

// RetryHandler handles retryable operations
type RetryHandler struct {
	config       RetryConfig
	errorHandler *ErrorHandler
}

// NewRetryHandler creates a new retry handler
func NewRetryHandler(config RetryConfig, errorHandler *ErrorHandler) *RetryHandler {
	// Set defaults
	if config.MaxAttempts <= 0 {
		config.MaxAttempts = 3
	}
	if config.BaseDelay <= 0 {
		config.BaseDelay = 100 * time.Millisecond
	}
	if config.MaxDelay <= 0 {
		config.MaxDelay = 5 * time.Second
	}
	if config.BackoffFactor <= 0 {
		config.BackoffFactor = 2.0
	}
	
	return &RetryHandler{
		config:       config,
		errorHandler: errorHandler,
	}
}

// Execute executes a function with retry logic
func (rh *RetryHandler) Execute(operation string, fn func() error) error {
	var lastErr error
	
	for attempt := 1; attempt <= rh.config.MaxAttempts; attempt++ {
		err := fn()
		if err == nil {
			return nil // Success
		}
		
		lastErr = err
		
		// Check if error is retryable
		if !rh.isRetryable(err) {
			rh.errorHandler.HandleError(err)
			return err
		}
		
		// Don't delay on the last attempt
		if attempt < rh.config.MaxAttempts {
			delay := rh.calculateDelay(attempt)
			time.Sleep(delay)
		}
		
		// Log retry attempt
		if rh.errorHandler != nil {
			retryErr := NewGenomicsError(ResourceError, "retry_handler", operation, 
				fmt.Sprintf("retry attempt %d/%d", attempt, rh.config.MaxAttempts))
			retryErr.WithCause(err)
			retryErr.WithContext("attempt", attempt)
			rh.errorHandler.HandleError(retryErr)
		}
	}
	
	// All attempts failed
	finalErr := NewGenomicsError(ResourceError, "retry_handler", operation, 
		fmt.Sprintf("operation failed after %d attempts", rh.config.MaxAttempts))
	finalErr.WithCause(lastErr)
	
	rh.errorHandler.HandleError(finalErr)
	return finalErr
}

// isRetryable checks if an error type is retryable
func (rh *RetryHandler) isRetryable(err error) bool {
	if ge, ok := err.(*GenomicsError); ok {
		for _, retryableType := range rh.config.RetryableTypes {
			if ge.Type == retryableType {
				return true
			}
		}
	}
	
	// Default retryable types
	defaultRetryable := []ErrorType{
		NetworkError, TimeoutError, CacheError, ResourceError,
	}
	
	if ge, ok := err.(*GenomicsError); ok {
		for _, retryableType := range defaultRetryable {
			if ge.Type == retryableType {
				return true
			}
		}
	}
	
	return false
}

// calculateDelay calculates the delay for a retry attempt
func (rh *RetryHandler) calculateDelay(attempt int) time.Duration {
	delay := time.Duration(float64(rh.config.BaseDelay) * 
		float64(attempt-1) * rh.config.BackoffFactor)
	
	if delay > rh.config.MaxDelay {
		delay = rh.config.MaxDelay
	}
	
	return delay
}

// Global error handler instance
var globalErrorHandler *ErrorHandler
var globalErrorHandlerOnce sync.Once

// GetGlobalErrorHandler returns the global error handler instance
func GetGlobalErrorHandler() *ErrorHandler {
	globalErrorHandlerOnce.Do(func() {
		config := ErrorHandlerConfig{
			EnableStackTrace: true,
			EnableContextLog: true,
			LogLevel:        LogLevelWarn,
			MaxStackFrames:  5,
			EnableMetrics:   true,
			ReportToExternal: false,
		}
		globalErrorHandler = NewErrorHandler(config)
	})
	
	return globalErrorHandler
}

// Convenience functions using global error handler

// HandleError handles an error using the global error handler
func HandleError(err error) {
	GetGlobalErrorHandler().HandleError(err)
}

// HandleErrorWithContext handles an error with additional context
func HandleErrorWithContext(err error, context map[string]interface{}) {
	if ge, ok := err.(*GenomicsError); ok {
		for key, value := range context {
			ge.WithContext(key, value)
		}
	}
	HandleError(err)
}