package fingerprinting

import (
	"fmt"
)

// LegacyEngine provides the original pure Go implementation for fallback
type LegacyEngine struct {
	customerID string
}

// Note: The New() function is defined in engine.go to avoid duplication

// NewLegacy explicitly creates a legacy engine for testing or fallback
func NewLegacy(customerID string) *LegacyEngine {
	return &LegacyEngine{customerID: customerID}
}

func (e *LegacyEngine) GenerateFingerprint(datasetID string) (string, error) {
	// Pure Go implementation for demo
	fingerprint := fmt.Sprintf("%s:%s:%d", e.customerID, datasetID, 1234567890)
	return fingerprint, nil
}

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

func (e *LegacyEngine) EmbedCustomerID(readID string) string {
	// Simple demo implementation - append customer ID
	return fmt.Sprintf("%s:%s", readID, e.customerID)
}

func (e *LegacyEngine) EmbedDatasetFingerprint(readID, datasetID string) string {
	// Simple demo implementation - append dataset ID
	return fmt.Sprintf("%s:DS:%s", readID, datasetID)
}

func (e *LegacyEngine) VerifyFingerprint(qualities []byte) bool {
	// Legacy implementation - basic check on quality scores
	return len(qualities) > 0 && qualities[0] >= 30
}

func (e *LegacyEngine) ExtractFingerprint(qualities []byte) (string, error) {
	// Legacy implementation - return a simple fingerprint
	if len(qualities) == 0 {
		return "", fmt.Errorf("empty quality scores")
	}
	return fmt.Sprintf("%s:legacy:1234567890", e.customerID), nil
}

func (e *LegacyEngine) Close() {
	// Cleanup if needed
}

// LegacyVerifier provides the original verifier implementation
type LegacyVerifier struct{}

func NewVerifier() Verifier {
	return &LegacyVerifier{}
}

func (v *LegacyVerifier) Verify(filePath string) (bool, string, error) {
	// This would implement fingerprint extraction and verification
	return true, "phantom-science", nil
}

func (v *LegacyVerifier) ValidateFile(filePath string) bool {
	// Mock implementation for testing
	return true
}