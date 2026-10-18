//go:build !release && cgo
// +build !release,cgo

package fingerprinting

import (
	"errors"
	"fmt"
	"os"
	"strings"
	"time"

	stegano "github.com/scttfrdmn/stegano/go"
)

// SteganoEngine implements Engine and SAMFingerprintingEngine on top of the
// stegano library (github.com/scttfrdmn/stegano/go).
//
// With WithKey, payloads are encrypted (Argon2id key derivation and an AEAD
// cipher, as implemented by stegano) and can only be read back with the same
// key. Without a key they are embedded in plaintext and anyone with the
// library can read them. Payloads written by stegano 0.5.x and earlier use an
// older wire format and cannot be extracted.
type SteganoEngine struct {
	customerID string
	key        []byte // nil: plaintext payloads
}

// NewSteganoEngine creates an engine backed by the stegano library.
func NewSteganoEngine(customerID string, opts ...Option) (*SteganoEngine, error) {
	if customerID == "" {
		return nil, NewFingerprintingError("new_engine", errors.New("customer ID is required"), ErrorCodeValidation)
	}
	return &SteganoEngine{customerID: customerID, key: applyOptions(opts).key}, nil
}

// newDefaultEngine returns the stegano-backed engine used by New.
func newDefaultEngine(customerID string, opts []Option) (Engine, error) {
	return NewSteganoEngine(customerID, opts...)
}

// newDefaultVerifier returns the stegano-backed Verifier used by NewVerifier.
func newDefaultVerifier(opts []Option) (Verifier, error) {
	return steganoVerifier{key: applyOptions(opts).key}, nil
}

// steganoVerifier implements Verifier with ExtractFromFASTQFile.
type steganoVerifier struct {
	key []byte // nil: plaintext fingerprints only
}

// Verify extracts the fingerprint from the FASTQ file at path. Missing or
// unreadable files are errors; a readable file with no fingerprint is not.
//
// stegano reports the same cryptographic error for a wrong key, a missing key
// and a key supplied for a plaintext payload, so a failed keyed extraction is
// retried without the key to tell a plaintext fingerprint from a wrong key.
func (v steganoVerifier) Verify(path string) (VerifyResult, error) {
	if _, err := os.Stat(path); err != nil {
		return VerifyResult{}, err
	}
	keyed := &SteganoEngine{key: v.key}
	fp, err := keyed.ExtractFromFASTQFile(path)
	if err == nil {
		return VerifyResult{Found: true, Fingerprint: fp, Encrypted: v.key != nil}, nil
	}
	if errors.Is(err, ErrWrongKey) {
		if fp, plainErr := (&SteganoEngine{}).ExtractFromFASTQFile(path); plainErr == nil {
			return VerifyResult{Found: true, Fingerprint: fp}, nil
		}
		return VerifyResult{}, err
	}
	var fe *FingerprintingError
	if errors.As(err, &fe) && fe.Code == ErrorCodeCorruption {
		return VerifyResult{}, nil
	}
	return VerifyResult{}, err
}

// LibraryVersion returns the version of the linked stegano library.
func LibraryVersion() string {
	return stegano.Version()
}

// steganoOptions returns opts plus the engine's encryption key, if any.
func (e *SteganoEngine) steganoOptions(opts ...stegano.Option) []stegano.Option {
	if e.key != nil {
		opts = append(opts, stegano.WithPassword(string(e.key)))
	}
	return opts
}

// GenerateFingerprint creates a fingerprint payload for a dataset.
func (e *SteganoEngine) GenerateFingerprint(datasetID string) (string, error) {
	return FormatFingerprint(e.customerID, datasetID, stegano.Version(), time.Now()), nil
}

// EmbedInQualities embeds fingerprint in the quality scores of fastqData,
// which must be complete FASTQ records, and returns the modified records.
func (e *SteganoEngine) EmbedInQualities(fastqData []byte, fingerprint string) ([]byte, error) {
	if err := ValidateInput(fastqData, fingerprint, "embed"); err != nil {
		return nil, err
	}
	out, err := stegano.Embed(fastqData, []byte(fingerprint), e.steganoOptions(stegano.WithFormat(stegano.FormatFASTQ))...)
	if err != nil {
		return nil, e.translateError("embed", err).
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
	payload, err := stegano.Extract(fastqData, e.steganoOptions(stegano.WithFormat(stegano.FormatFASTQ))...)
	if err != nil {
		return "", e.translateError("extract", err).WithContext("data_length", len(fastqData))
	}
	return string(payload), nil
}

// EmbedInFASTQFile embeds fingerprint in the FASTQ file inPath and writes the
// result to outPath, which may be inPath. The file is streamed, so memory use
// does not grow with its size. outPath is replaced atomically and keeps
// inPath's permission bits; on error it is left untouched.
func (e *SteganoEngine) EmbedInFASTQFile(inPath, outPath, fingerprint string) error {
	if err := validateFingerprint(fingerprint, "embed_file"); err != nil {
		return err
	}
	info, err := os.Stat(inPath)
	if err != nil {
		return NewFingerprintingError("embed_file", err, ErrorCodeSystem).WithContext("path", inPath)
	}
	if err := stegano.EmbedFile(inPath, outPath, []byte(fingerprint), e.steganoOptions(stegano.WithFormat(stegano.FormatFASTQ))...); err != nil {
		return e.translateError("embed_file", err).
			WithContext("path", inPath).
			WithContext("file_size", info.Size())
	}
	if err := os.Chmod(outPath, info.Mode().Perm()); err != nil {
		return NewFingerprintingError("embed_file", err, ErrorCodeSystem).WithContext("path", outPath)
	}
	return nil
}

// ExtractFromFASTQFile extracts a fingerprint from the FASTQ file at path.
// Only the start of the file is read.
func (e *SteganoEngine) ExtractFromFASTQFile(path string) (string, error) {
	payload, err := stegano.ExtractFile(path, e.steganoOptions(stegano.WithFormat(stegano.FormatFASTQ))...)
	if err != nil {
		return "", e.translateError("extract_file", err).WithContext("path", path)
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
	if err := validateFingerprint(fingerprint, "embed_sam"); err != nil {
		return nil, err
	}
	out, err := stegano.Embed(samData, []byte(fingerprint),
		e.steganoOptions(stegano.WithFormat(stegano.FormatSAM), stegano.WithStrategy("auxiliary"))...)
	if err != nil {
		return nil, e.translateError("embed_sam", err).WithContext("data_length", len(samData))
	}
	return out, nil
}

// ExtractFromSAMRecords extracts a fingerprint embedded by EmbedInSAMRecords.
func (e *SteganoEngine) ExtractFromSAMRecords(samData []byte) (string, error) {
	if len(samData) == 0 {
		return "", NewFingerprintingError("extract_sam", ErrEmptyInput, ErrorCodeValidation)
	}
	payload, err := stegano.Extract(samData,
		e.steganoOptions(stegano.WithFormat(stegano.FormatSAM), stegano.WithStrategy("auxiliary"))...)
	if err != nil {
		return "", e.translateError("extract_sam", err).WithContext("data_length", len(samData))
	}
	return string(payload), nil
}

// Close releases engine resources. The stegano library holds no per-engine
// state, so this is a no-op.
func (e *SteganoEngine) Close() {}

// translateError converts a stegano library error into a FingerprintingError
// with the matching error code. A cryptographic failure while extracting is a
// key problem: ErrWrongKey when the engine has a key, ErrKeyRequired when not.
func (e *SteganoEngine) translateError(op string, err error) *FingerprintingError {
	var se *stegano.Error
	if !errors.As(err, &se) {
		return NewFingerprintingError(op, err, ErrorCodeUnknown)
	}

	extracting := strings.HasPrefix(op, "extract")
	opFailed := ErrEmbeddingFailed
	if extracting {
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
	case stegano.ErrCryptographic:
		switch {
		case extracting && e.key != nil:
			fe = NewFingerprintingError(op, ErrWrongKey, ErrorCodePermission)
		case extracting:
			fe = NewFingerprintingError(op, ErrKeyRequired, ErrorCodePermission)
		default:
			fe = NewFingerprintingError(op, opFailed, ErrorCodeSystem)
		}
	case stegano.ErrIO, stegano.ErrOutputSizeChanged, stegano.ErrPanic:
		fe = NewFingerprintingError(op, opFailed, ErrorCodeSystem)
	default:
		fe = NewFingerprintingError(op, opFailed, ErrorCodeUnknown)
	}
	return fe.WithContext("stegano_error", se.Message)
}
