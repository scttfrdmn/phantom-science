package fingerprinting

// Engine interface defines the fingerprinting operations
type Engine interface {
	GenerateFingerprint(datasetID string) (string, error)
	EmbedInQualities(qualities []byte, fingerprint string) ([]byte, error)
	EmbedCustomerID(readID string) string
	EmbedDatasetFingerprint(readID, datasetID string) string
	VerifyFingerprint(qualities []byte) bool
	ExtractFingerprint(qualities []byte) (string, error)
	Close()
}

// Verifier checks files for an embedded fingerprint.
type Verifier interface {
	// Verify reports whether the FASTQ file at path carries a fingerprint and
	// returns it. A readable file without one gives (false, "", nil); err is
	// non-nil only when the file cannot be checked.
	Verify(path string) (found bool, fingerprint string, err error)
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
