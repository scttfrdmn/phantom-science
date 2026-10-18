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

// Verifier interface for file verification operations
type Verifier interface {
	Verify(filePath string) (bool, string, error)
	ValidateFile(filePath string) bool
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
