//go:build !release && cgo
// +build !release,cgo

package fingerprinting

import (
	"errors"
	"fmt"
	"time"

	stegano "github.com/scttfrdmn/stegano/go"
)

// SteganoEngine implements Engine and SAMFingerprintingEngine on top of the
// stegano library (github.com/scttfrdmn/stegano/go).
//
// Payloads are embedded unencrypted, so anyone with the library can read a
// fingerprint back. Payloads written by stegano 0.5.x and earlier use an
// older wire format and cannot be extracted.
type SteganoEngine struct {
	customerID string
}

// NewSteganoEngine creates an engine backed by the stegano library.
func NewSteganoEngine(customerID string) (*SteganoEngine, error) {
	if customerID == "" {
		return nil, NewFingerprintingError("new_engine", errors.New("customer ID is required"), ErrorCodeValidation)
	}
	return &SteganoEngine{customerID: customerID}, nil
}

// newDefaultEngine returns the stegano-backed engine used by New.
func newDefaultEngine(customerID string) (Engine, error) {
	return NewSteganoEngine(customerID)
}

// LibraryVersion returns the version of the linked stegano library.
func LibraryVersion() string {
	return stegano.Version()
}

// GenerateFingerprint creates a fingerprint payload for a dataset.
func (e *SteganoEngine) GenerateFingerprint(datasetID string) (string, error) {
	return fmt.Sprintf("%s:%s:stegano-%s:ts-%d",
		e.customerID, datasetID, stegano.Version(), time.Now().Unix()), nil
}

// EmbedInQualities embeds fingerprint in the quality scores of fastqData,
// which must be complete FASTQ records, and returns the modified records.
func (e *SteganoEngine) EmbedInQualities(fastqData []byte, fingerprint string) ([]byte, error) {
	if err := ValidateInput(fastqData, fingerprint, "embed"); err != nil {
		return nil, err
	}
	out, err := stegano.Embed(fastqData, []byte(fingerprint), stegano.WithFormat(stegano.FormatFASTQ))
	if err != nil {
		return nil, translateSteganoError("embed", err).
			WithContext("fingerprint_length", len(fingerprint)).
			WithContext("data_length", len(fastqData))
	}
	return out, nil
}

// ExtractFingerprint extracts a fingerprint from FASTQ records.
func (e *SteganoEngine) ExtractFingerprint(fastqData []byte) (string, error) {
	if len(fastqData) == 0 {
		return "", NewFingerprintingError("extract", ErrEmptyInput, ErrorCodeValidation)
	}
	if !isValidFASTQData(fastqData) {
		return "", NewFingerprintingError("extract", ErrInvalidFormat, ErrorCodeValidation).
			WithContext("data_length", len(fastqData))
	}
	payload, err := stegano.Extract(fastqData, stegano.WithFormat(stegano.FormatFASTQ))
	if err != nil {
		return "", translateSteganoError("extract", err).WithContext("data_length", len(fastqData))
	}
	return string(payload), nil
}

// VerifyFingerprint reports whether fastqData carries an extractable fingerprint.
func (e *SteganoEngine) VerifyFingerprint(fastqData []byte) bool {
	_, err := e.ExtractFingerprint(fastqData)
	return err == nil
}

// EmbedCustomerID tags a read ID with the customer identifier.
func (e *SteganoEngine) EmbedCustomerID(readID string) string {
	return fmt.Sprintf("%s:cust:%s", readID, e.customerID)
}

// EmbedDatasetFingerprint tags a read ID with the dataset identifier.
func (e *SteganoEngine) EmbedDatasetFingerprint(readID, datasetID string) string {
	return fmt.Sprintf("%s:ds:%s", readID, datasetID)
}

// EmbedInSAMRecords embeds fingerprint in SAM records using the auxiliary
// strategy (a zs:Z: optional tag), leaving sequences and qualities untouched.
// The returned data is longer than samData.
func (e *SteganoEngine) EmbedInSAMRecords(samData []byte, fingerprint string) ([]byte, error) {
	if len(samData) == 0 {
		return nil, NewFingerprintingError("embed_sam", ErrEmptyInput, ErrorCodeValidation)
	}
	if fingerprint == "" || len(fingerprint) > MaxFingerprintLength || containsInvalidCharacters(fingerprint) {
		return nil, NewFingerprintingError("embed_sam", ErrInvalidFingerprint, ErrorCodeValidation).
			WithContext("fingerprint_length", len(fingerprint))
	}
	out, err := stegano.Embed(samData, []byte(fingerprint),
		stegano.WithFormat(stegano.FormatSAM), stegano.WithStrategy("auxiliary"))
	if err != nil {
		return nil, translateSteganoError("embed_sam", err).WithContext("data_length", len(samData))
	}
	return out, nil
}

// ExtractFromSAMRecords extracts a fingerprint embedded by EmbedInSAMRecords.
func (e *SteganoEngine) ExtractFromSAMRecords(samData []byte) (string, error) {
	if len(samData) == 0 {
		return "", NewFingerprintingError("extract_sam", ErrEmptyInput, ErrorCodeValidation)
	}
	payload, err := stegano.Extract(samData,
		stegano.WithFormat(stegano.FormatSAM), stegano.WithStrategy("auxiliary"))
	if err != nil {
		return "", translateSteganoError("extract_sam", err).WithContext("data_length", len(samData))
	}
	return string(payload), nil
}

// Close releases engine resources. The stegano library holds no per-engine
// state, so this is a no-op.
func (e *SteganoEngine) Close() {}

// translateSteganoError converts a stegano library error into a
// FingerprintingError with the matching error code.
func translateSteganoError(op string, err error) *FingerprintingError {
	var se *stegano.Error
	if !errors.As(err, &se) {
		return NewFingerprintingError(op, err, ErrorCodeUnknown)
	}

	opFailed := ErrEmbeddingFailed
	if op == "extract" || op == "extract_sam" {
		opFailed = ErrExtractionFailed
	}

	var fe *FingerprintingError
	switch se.Code {
	case stegano.ErrInvalidParameter, stegano.ErrFormatNotSupported:
		fe = NewFingerprintingError(op, ErrInvalidFormat, ErrorCodeValidation)
	case stegano.ErrInsufficientCapacity:
		fe = NewFingerprintingError(op, ErrCapacityExceeded, ErrorCodeCapacity)
	case stegano.ErrDataCorrupted:
		fe = NewFingerprintingError(op, ErrCorruptedData, ErrorCodeCorruption)
	case stegano.ErrGPU:
		fe = NewFingerprintingError(op, ErrGPUError, ErrorCodeGPU)
	case stegano.ErrCryptographic, stegano.ErrIO, stegano.ErrOutputSizeChanged, stegano.ErrPanic:
		fe = NewFingerprintingError(op, opFailed, ErrorCodeSystem)
	default:
		fe = NewFingerprintingError(op, opFailed, ErrorCodeUnknown)
	}
	return fe.WithContext("stegano_error", se.Message)
}
