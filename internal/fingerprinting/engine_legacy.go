package fingerprinting

import (
	"fmt"
)

// LegacyEngine is the pure-Go fallback used when the stegano library is not
// linked. Only its read-ID tagging (EmbedCustomerID, EmbedDatasetFingerprint)
// is meaningful: its quality-score methods are placeholders that cannot embed
// a recoverable fingerprint, so generation refuses quality-score
// fingerprinting in builds that use it.
type LegacyEngine struct {
	customerID string
}

// Note: The New() function is defined in engine.go to avoid duplication

// NewLegacy creates a legacy engine. New and the release build use it when
// stegano is unavailable.
func NewLegacy(customerID string) *LegacyEngine {
	return &LegacyEngine{customerID: customerID}
}

// GenerateFingerprint returns "<customer>:<dataset>:1234567890"; the
// timestamp field is a fixed placeholder.
func (e *LegacyEngine) GenerateFingerprint(datasetID string) (string, error) {
	// Pure Go implementation for demo
	fingerprint := fmt.Sprintf("%s:%s:%d", e.customerID, datasetID, 1234567890)
	return fingerprint, nil
}

// EmbedInQualities is a placeholder: it sets or clears the low bit of the
// first len(fingerprint) quality bytes that are at least 30, and nothing can
// recover the fingerprint from the result.
func (e *LegacyEngine) EmbedInQualities(qualities []byte, fingerprint string) ([]byte, error) {
	// Pure Go implementation for demo - simple LSB embedding
	output := make([]byte, len(qualities))
	copy(output, qualities)

	// Simulate embedding fingerprint in quality scores
	for i := 0; i < len(output) && i < len(fingerprint); i++ {
		if output[i] >= 30 {
			// Simple LSB modification
			if fingerprint[i%len(fingerprint)] > 'M' {
				output[i] = output[i] | 1
			} else {
				output[i] = output[i] &^ 1
			}
		}
	}

	return output, nil
}

// EmbedCustomerID appends ":<customer>" to readID.
func (e *LegacyEngine) EmbedCustomerID(readID string) string {
	// Simple demo implementation - append customer ID
	return fmt.Sprintf("%s:%s", readID, e.customerID)
}

// EmbedDatasetFingerprint appends ":DS:<dataset>" to readID.
func (e *LegacyEngine) EmbedDatasetFingerprint(readID, datasetID string) string {
	// Simple demo implementation - append dataset ID
	return fmt.Sprintf("%s:DS:%s", readID, datasetID)
}

// VerifyFingerprint is a placeholder: it reports whether the first quality
// byte is at least 30, which says nothing about a fingerprint.
func (e *LegacyEngine) VerifyFingerprint(qualities []byte) bool {
	// Legacy implementation - basic check on quality scores
	return len(qualities) > 0 && qualities[0] >= 30
}

// ExtractFingerprint is a placeholder: for any non-empty input it returns
// "<customer>:legacy:1234567890" without reading the data.
func (e *LegacyEngine) ExtractFingerprint(qualities []byte) (string, error) {
	// Legacy implementation - return a simple fingerprint
	if len(qualities) == 0 {
		return "", fmt.Errorf("empty quality scores")
	}
	return fmt.Sprintf("%s:legacy:1234567890", e.customerID), nil
}

// Close is a no-op.
func (e *LegacyEngine) Close() {
	// Cleanup if needed
}
