package authenticity

import (
	"fmt"

	"github.com/scttfrdmn/phantom-science/internal/errors"
	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/security"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// IntegrationManager manages integration between authenticity engine and existing systems
type IntegrationManager struct {
	authenticityEngine   *Engine
	fingerprintingEngine fingerprinting.Engine
	config              types.SequencerConfig
}

// NewIntegrationManager creates a new integration manager
func NewIntegrationManager(authenticityEngine *Engine, fingerprintingEngine fingerprinting.Engine, config types.SequencerConfig) *IntegrationManager {
	return &IntegrationManager{
		authenticityEngine:   authenticityEngine,
		fingerprintingEngine: fingerprintingEngine,
		config:              config,
	}
}

// GenerateAuthenticFASTQFile generates a FASTQ file with biologically authentic sequences
func (im *IntegrationManager) GenerateAuthenticFASTQFile(filename string, batch, read int) error {
	// Validate file path against the configured output directory
	if err := security.ValidateFilePath(filename, im.config.OutputDir); err != nil {
		return errors.WrapError(err, errors.SecurityError, "integration_manager", "generate_authentic_fastq", "invalid file path")
	}
	
	file, err := security.SecureFileCreate(filename, im.config.OutputDir)
	if err != nil {
		return errors.WrapError(err, errors.ResourceError, "integration_manager", "generate_authentic_fastq", "failed to create secure file")
	}
	defer file.Close()
	
	for readNum := 1; readNum <= im.config.ReadsPerBatch; readNum++ {
		// Generate read ID with batch and read info
		baseReadID := fmt.Sprintf("@%s:%d:HXXXXXX:%d:1101:%d:%d %d:N:0:1", 
			im.config.SampleName, im.config.Lane, batch, readNum*1000, readNum*2000, read)
		
		// Apply fingerprinting to read ID if enabled
		readID := baseReadID
		if im.config.Fingerprinting.Enable && im.fingerprintingEngine != nil {
			if im.config.Fingerprinting.CustomerFingerprint {
				readID = im.fingerprintingEngine.EmbedCustomerID(baseReadID)
			}
			if im.config.Fingerprinting.DatasetFingerprint && im.config.DatasetID != "" {
				readID = im.fingerprintingEngine.EmbedDatasetFingerprint(readID, im.config.DatasetID)
			}
		}
		
		// Generate biologically authentic DNA sequence
		sequence, err := im.authenticityEngine.GenerateAuthenticSequence(im.config.ReadLength, "1")
		if err != nil {
			return errors.WrapError(err, errors.BiologicalModelError, "integration_manager", "generate_authentic_fastq", "failed to generate authentic sequence")
		}
		
		// Generate quality scores with authenticity-aware patterns
		qualities, err := im.generateAuthenticQualityScores(sequence)
		if err != nil {
			return errors.WrapError(err, errors.BiologicalModelError, "integration_manager", "generate_authentic_fastq", "failed to generate authentic quality scores")
		}
		
		// Apply steganographic fingerprinting to quality scores if enabled
		if im.config.Fingerprinting.Enable && im.fingerprintingEngine != nil {
			// Generate fingerprint for this dataset
			fingerprint, err := im.fingerprintingEngine.GenerateFingerprint(im.config.DatasetID)
			if err == nil {
				// Embed fingerprint in quality scores using stegano
				if modifiedQualities, err := im.fingerprintingEngine.EmbedInQualities(qualities, fingerprint); err == nil {
					qualities = modifiedQualities
				}
			}
		}
		
		// Write FASTQ record
		fmt.Fprintf(file, "%s\n%s\n+\n%s\n", readID, string(sequence), string(qualities))
	}
	
	return nil
}

// generateAuthenticQualityScores generates quality scores that reflect biological authenticity
func (im *IntegrationManager) generateAuthenticQualityScores(sequence []byte) ([]byte, error) {
	qualityChars := "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
	qualities := make([]byte, len(sequence))
	
	for i, base := range sequence {
		// Base quality score with decay
		qualityValue := im.config.BaseQuality - int(float64(i)*im.config.QualityDecay)
		
		// Apply biological authenticity patterns
		qualityValue = im.applyBiologicalQualityPatterns(qualityValue, base, i, sequence)
		
		// Ensure quality is within valid range
		if qualityValue < 2 {
			qualityValue = 2
		}
		if qualityValue > len(qualityChars)-1 {
			qualityValue = len(qualityChars) - 1
		}
		
		qualities[i] = qualityChars[qualityValue]
	}
	
	return qualities, nil
}

// applyBiologicalQualityPatterns applies biological authenticity patterns to quality scores
func (im *IntegrationManager) applyBiologicalQualityPatterns(baseQuality int, base byte, position int, sequence []byte) int {
	adjustedQuality := baseQuality
	
	// GC bias: lower quality for GC-rich regions
	if base == 'G' || base == 'C' {
		gcBias := im.calculateLocalGCBias(position, sequence)
		if gcBias > 0.6 { // High GC content
			adjustedQuality -= 2
		} else if gcBias > 0.4 { // Medium GC content
			adjustedQuality -= 1
		}
	}
	
	// Homopolymer runs: lower quality for long runs
	homopolymerLength := im.calculateHomopolymerLength(position, sequence)
	if homopolymerLength > 3 {
		adjustedQuality -= (homopolymerLength - 3) * 2
	}
	
	// Position-specific effects: quality degradation toward ends
	if position < 10 || position > len(sequence)-10 {
		adjustedQuality -= 1
	}
	
	return adjustedQuality
}

// calculateLocalGCBias calculates local GC content around a position
func (im *IntegrationManager) calculateLocalGCBias(position int, sequence []byte) float64 {
	windowSize := 20
	start := position - windowSize/2
	end := position + windowSize/2
	
	if start < 0 {
		start = 0
	}
	if end >= len(sequence) {
		end = len(sequence) - 1
	}
	
	gcCount := 0
	total := 0
	
	for i := start; i <= end; i++ {
		if sequence[i] == 'G' || sequence[i] == 'C' {
			gcCount++
		}
		total++
	}
	
	if total == 0 {
		return 0.0
	}
	
	return float64(gcCount) / float64(total)
}

// calculateHomopolymerLength calculates the length of homopolymer run at a position
func (im *IntegrationManager) calculateHomopolymerLength(position int, sequence []byte) int {
	if position >= len(sequence) {
		return 0
	}
	
	base := sequence[position]
	length := 1
	
	// Count backwards
	for i := position - 1; i >= 0 && sequence[i] == base; i-- {
		length++
	}
	
	// Count forwards
	for i := position + 1; i < len(sequence) && sequence[i] == base; i++ {
		length++
	}
	
	return length
}

// EnhancedGenomicsEngine wraps the existing genomics engine with authenticity features
type EnhancedGenomicsEngine struct {
	originalEngine       genomicsEngine // Would be the original engine
	authenticityEngine   *Engine
	integrationManager   *IntegrationManager
	useAuthenticity      bool
}

// genomicsEngine represents the original genomics engine interface
type genomicsEngine interface {
	Generate() error
	GenerateFASTQFile(filename string, batch, read int) error
	Close()
}

// NewEnhancedGenomicsEngine creates an enhanced genomics engine with authenticity features
func NewEnhancedGenomicsEngine(originalEngine genomicsEngine, authenticityEngine *Engine, config types.SequencerConfig) *EnhancedGenomicsEngine {
	// Create a dummy fingerprinting engine for integration
	fingerprintingEngine := fingerprinting.New("phantom-science")
	
	integrationManager := NewIntegrationManager(authenticityEngine, fingerprintingEngine, config)
	
	return &EnhancedGenomicsEngine{
		originalEngine:     originalEngine,
		authenticityEngine: authenticityEngine,
		integrationManager: integrationManager,
		useAuthenticity:    true, // Default to using authenticity
	}
}

// Generate generates genomic data with optional authenticity features
func (ege *EnhancedGenomicsEngine) Generate() error {
	// Initialize authenticity engine if not already done
	if ege.useAuthenticity && !ege.authenticityEngine.IsInitialized() {
		if err := ege.authenticityEngine.Initialize(); err != nil {
			// Fall back to original engine if authenticity initialization fails
			genomicsErr := errors.WrapError(err, errors.InitializationError, "enhanced_genomics_engine", "generate", "failed to initialize authenticity engine, falling back to original")
			errors.HandleError(genomicsErr)
			ege.useAuthenticity = false
		}
	}
	
	// Use original engine's Generate method
	// In a full implementation, this would be modified to use authentic sequence generation
	return ege.originalEngine.Generate()
}

// GenerateFASTQFile generates a FASTQ file with authenticity features
func (ege *EnhancedGenomicsEngine) GenerateFASTQFile(filename string, batch, read int) error {
	if ege.useAuthenticity {
		return ege.integrationManager.GenerateAuthenticFASTQFile(filename, batch, read)
	}
	
	// Fall back to original engine
	return ege.originalEngine.GenerateFASTQFile(filename, batch, read)
}

// Close releases resources
func (ege *EnhancedGenomicsEngine) Close() {
	if ege.authenticityEngine != nil {
		ege.authenticityEngine.Close()
	}
	if ege.originalEngine != nil {
		ege.originalEngine.Close()
	}
}

// SetAuthenticityEnabled enables or disables authenticity features
func (ege *EnhancedGenomicsEngine) SetAuthenticityEnabled(enabled bool) {
	ege.useAuthenticity = enabled
}

// IsAuthenticityEnabled returns whether authenticity features are enabled
func (ege *EnhancedGenomicsEngine) IsAuthenticityEnabled() bool {
	return ege.useAuthenticity
}

// GetAuthenticityStatus returns the status of the authenticity engine
func (ege *EnhancedGenomicsEngine) GetAuthenticityStatus() map[string]interface{} {
	if ege.authenticityEngine == nil {
		return map[string]interface{}{
			"available": false,
			"message":   "Authenticity engine not initialized",
		}
	}
	
	status := ege.authenticityEngine.GetStatus()
	status["available"] = true
	status["enabled"] = ege.useAuthenticity
	
	return status
}

// ValidateLastGeneratedSequence validates the last generated sequence
func (ege *EnhancedGenomicsEngine) ValidateLastGeneratedSequence(sequence []byte) (map[string]interface{}, error) {
	if !ege.useAuthenticity || ege.authenticityEngine == nil {
		return nil, errors.NewResourceError("enhanced_genomics_engine", "validate_last_generated_sequence", "authenticity engine not available")
	}
	
	report, err := ege.authenticityEngine.ValidateSequence(sequence)
	if err != nil {
		return nil, err
	}
	
	// Convert validation report to map for easier consumption
	result := map[string]interface{}{
		"overall_passed": report.OverallPassed,
		"total_tests":    report.TotalTests,
		"passed_tests":   report.PassedTests,
		"failed_tests":   report.FailedTests,
		"average_score":  report.AverageScore,
		"timestamp":      report.Timestamp,
		"execution_time": report.ExecutionTime,
		"recommendations": report.RecommendedActions,
	}
	
	// Add individual test results
	var testResults []map[string]interface{}
	for _, testResult := range report.Results {
		testResults = append(testResults, map[string]interface{}{
			"test_name":      testResult.TestName,
			"test_type":      testResult.TestType.String(),
			"passed":         testResult.Passed,
			"score":          testResult.Score,
			"threshold":      testResult.Threshold,
			"message":        testResult.Message,
			"execution_time": testResult.ExecutionTime,
		})
	}
	result["test_results"] = testResults
	
	return result, nil
}

// GetAuthenticityScore returns the authenticity score for a sequence
func (ege *EnhancedGenomicsEngine) GetAuthenticityScore(sequence []byte) (float64, error) {
	if !ege.useAuthenticity || ege.authenticityEngine == nil {
		return 0.0, errors.NewResourceError("enhanced_genomics_engine", "get_authenticity_score", "authenticity engine not available")
	}
	
	return ege.authenticityEngine.GetAuthenticityScore(sequence)
}

// UpdateAuthenticityConfig updates the authenticity configuration
func (ege *EnhancedGenomicsEngine) UpdateAuthenticityConfig(config types.SequencerConfig) error {
	if ege.authenticityEngine == nil {
		return errors.NewResourceError("enhanced_genomics_engine", "update_authenticity_config", "authenticity engine not available")
	}
	
	return ege.authenticityEngine.UpdateFromSequencerConfig(config)
}

// GetReferenceStatistics returns statistics about the reference genome data
func (ege *EnhancedGenomicsEngine) GetReferenceStatistics() (map[string]interface{}, error) {
	if !ege.useAuthenticity || ege.authenticityEngine == nil {
		return nil, errors.NewResourceError("enhanced_genomics_engine", "get_reference_statistics", "authenticity engine not available")
	}
	
	return ege.authenticityEngine.GetReferenceStatistics()
}

// ConfigureForPlatform configures the authenticity engine for a specific sequencing platform
func (ege *EnhancedGenomicsEngine) ConfigureForPlatform(platform string) error {
	if ege.authenticityEngine == nil {
		return errors.NewResourceError("enhanced_genomics_engine", "configure_for_platform", "authenticity engine not available")
	}
	
	// Create a temporary sequencer config for platform-specific settings
	tempConfig := types.SequencerConfig{
		SequencerModel: platform,
	}
	
	return ege.authenticityEngine.UpdateFromSequencerConfig(tempConfig)
}

// GetIntegrationInfo returns information about the integration with existing systems
func (ege *EnhancedGenomicsEngine) GetIntegrationInfo() map[string]interface{} {
	return map[string]interface{}{
		"authenticity_available": ege.authenticityEngine != nil,
		"authenticity_enabled":   ege.useAuthenticity,
		"integration_version":    "Phase 3a.1",
		"supported_platforms":    []string{"illumina", "pacbio", "nanopore"},
		"features": map[string]bool{
			"authentic_sequence_generation": true,
			"kmer_frequency_matching":       true,
			"jensen_shannon_validation":     true,
			"gc_bias_quality_patterns":      true,
			"homopolymer_error_modeling":    true,
			"fingerprinting_integration":    true,
		},
	}
}

// BackwardsCompatibilityWrapper provides backwards compatibility with existing code
type BackwardsCompatibilityWrapper struct {
	enhanced *EnhancedGenomicsEngine
}

// NewBackwardsCompatibilityWrapper creates a wrapper for backwards compatibility
func NewBackwardsCompatibilityWrapper(enhanced *EnhancedGenomicsEngine) *BackwardsCompatibilityWrapper {
	return &BackwardsCompatibilityWrapper{enhanced: enhanced}
}

// Generate implements the original Generate interface
func (bcw *BackwardsCompatibilityWrapper) Generate() error {
	return bcw.enhanced.Generate()
}

// GenerateFASTQFile implements the original GenerateFASTQFile interface
func (bcw *BackwardsCompatibilityWrapper) GenerateFASTQFile(filename string, batch, read int) error {
	return bcw.enhanced.GenerateFASTQFile(filename, batch, read)
}

// Close implements the original Close interface
func (bcw *BackwardsCompatibilityWrapper) Close() {
	bcw.enhanced.Close()
}

// GetEnhancedFeatures provides access to enhanced features
func (bcw *BackwardsCompatibilityWrapper) GetEnhancedFeatures() *EnhancedGenomicsEngine {
	return bcw.enhanced
}