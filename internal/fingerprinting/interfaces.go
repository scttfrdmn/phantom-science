// Package fingerprinting attributes generated data to a customer and
// dataset without changing it (docs/CHARTER.md, principle 2).
//
// Generation uses keyed-seed fingerprints (keyedseed.go): the run seed is
// DeriveSeed(key, customer, dataset, salt), an HMAC-SHA256 under the
// fingerprint key, so the output is an ordinary run whose seed only the key
// holder can compute. It is pure Go and needs no CGO. Verification
// regenerates the run from the key and the run's records and matches suspect
// reads against it (Matcher). ResolveKey reads the key from a key file or the
// PHANTOM_FINGERPRINT_KEY environment variable.
//
// The quality-score fingerprints written by phantom-science 0.9.0 and
// earlier changed quality values, so generation no longer writes them
// (ErrQualityScoresLayer). SteganoEngine, built with CGO and without the
// release tag, wraps the stegano library (github.com/scttfrdmn/stegano/go);
// the Verifier returned by NewVerifier extracts those legacy fingerprints
// from old files, decrypting them with WithKey. Builds without the stegano
// library cannot read them: NewVerifier returns an error.
//
// ReadTags formats the optional read-ID labels (XC:Z:, XD:Z: header
// fields). They are visible text, not a fingerprint.
package fingerprinting

// VerifyResult is the outcome of checking one file for a legacy
// quality-score fingerprint.
type VerifyResult struct {
	Found       bool   // a fingerprint was extracted
	Fingerprint string // the extracted payload
	Encrypted   bool   // the payload was encrypted and decrypted with the verifier's key
}

// Verifier checks files for a legacy quality-score fingerprint.
type Verifier interface {
	// Verify checks the FASTQ file at path. A readable file without a
	// fingerprint gives a zero VerifyResult and no error. err is non-nil when
	// the file cannot be checked, including ErrKeyRequired for an encrypted
	// fingerprint without a key and ErrWrongKey when the key does not match.
	Verify(path string) (VerifyResult, error)
}

// FileFingerprinter embeds and extracts legacy quality-score fingerprints in
// whole FASTQ files. Implemented by SteganoEngine. Generation does not use
// it; tests use it to make files like those of phantom-science 0.9.0.
type FileFingerprinter interface {
	// EmbedInFASTQFile embeds fingerprint in the quality scores of the FASTQ
	// file inPath and writes the result to outPath, which may be inPath.
	EmbedInFASTQFile(inPath, outPath, fingerprint string) error
	// ExtractFromFASTQFile extracts a fingerprint embedded by EmbedInFASTQFile.
	ExtractFromFASTQFile(path string) (string, error)
}

// SAMFingerprintingEngine embeds and extracts fingerprints in SAM records.
// Implemented by SteganoEngine; no command writes SAM.
type SAMFingerprintingEngine interface {
	// EmbedInSAMRecords embeds a fingerprint in SAM alignment records via a zs:Z: optional tag.
	EmbedInSAMRecords(samData []byte, fingerprint string) ([]byte, error)
	// ExtractFromSAMRecords extracts a fingerprint from SAM alignment records.
	ExtractFromSAMRecords(samData []byte) (string, error)
}
