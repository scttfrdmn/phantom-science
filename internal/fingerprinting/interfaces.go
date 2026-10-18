// Package fingerprinting marks generated data with a customer and dataset
// fingerprint.
//
// SteganoEngine, built with CGO and without the release tag, uses the stegano
// library (github.com/scttfrdmn/stegano/go, v0.8.0). Generation embeds one
// fingerprint per FASTQ file in its quality scores with EmbedInFASTQFile;
// SAM data carries it in a zs:Z: tag. With WithKey the payload is encrypted
// (Argon2id key derivation, AES-256-GCM); ResolveKey reads the key from a key
// file or the PHANTOM_FINGERPRINT_KEY environment variable. Without a key the
// payload is plaintext and anyone with the library can read it; generation
// writes plaintext fingerprints only when the configuration sets
// plaintext: true. The Verifier returned by NewVerifier extracts
// and, with a key, decrypts the fingerprint from a FASTQ file.
//
// Builds without the stegano library (CGO_ENABLED=0 or -tags release) fall
// back to LegacyEngine, which only tags read IDs; it cannot embed or recover
// a fingerprint, and NewVerifier returns an error.
package fingerprinting

// Engine is the set of fingerprinting operations used by the generator.
type Engine interface {
	GenerateFingerprint(datasetID string) (string, error)
	EmbedInQualities(qualities []byte, fingerprint string) ([]byte, error)
	EmbedCustomerID(readID string) string
	EmbedDatasetFingerprint(readID, datasetID string) string
	VerifyFingerprint(qualities []byte) bool
	ExtractFingerprint(qualities []byte) (string, error)
	Close()
}

// VerifyResult is the outcome of checking one file for a fingerprint.
type VerifyResult struct {
	Found       bool   // a fingerprint was extracted
	Fingerprint string // the extracted payload
	Encrypted   bool   // the payload was encrypted and decrypted with the verifier's key
}

// Verifier checks files for an embedded fingerprint.
type Verifier interface {
	// Verify checks the FASTQ file at path. A readable file without a
	// fingerprint gives a zero VerifyResult and no error. err is non-nil when
	// the file cannot be checked, including ErrKeyRequired for an encrypted
	// fingerprint without a key and ErrWrongKey when the key does not match.
	Verify(path string) (VerifyResult, error)
}

// FileFingerprinter embeds and extracts fingerprints in whole FASTQ files.
// Implemented by SteganoEngine; use a type assertion to check.
type FileFingerprinter interface {
	// EmbedInFASTQFile embeds fingerprint in the quality scores of the FASTQ
	// file inPath and writes the result to outPath, which may be inPath.
	EmbedInFASTQFile(inPath, outPath, fingerprint string) error
	// ExtractFromFASTQFile extracts a fingerprint embedded by EmbedInFASTQFile.
	ExtractFromFASTQFile(path string) (string, error)
}

// SAMFingerprintingEngine extends Engine with SAM-level fingerprinting capabilities.
// Implemented by SteganoEngine; use a type assertion to check.
type SAMFingerprintingEngine interface {
	Engine
	// EmbedInSAMRecords embeds a fingerprint in SAM alignment records via a zs:Z: optional tag.
	EmbedInSAMRecords(samData []byte, fingerprint string) ([]byte, error)
	// ExtractFromSAMRecords extracts a fingerprint from SAM alignment records.
	ExtractFromSAMRecords(samData []byte) (string, error)
}
