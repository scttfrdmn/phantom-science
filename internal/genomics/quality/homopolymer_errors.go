package quality

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sort"
	"strings"
	"time"
)

// HomopolymerErrorModel models homopolymer-induced sequencing errors
type HomopolymerErrorModel struct {
	errorRates     map[int]float64    // Length -> error rate
	platform       string             // Sequencing platform
	errorTypes     []string           // Error types: insertion, deletion, substitution
	contextEffects bool               // Whether to consider sequence context
	baseErrorRate  float64            // Base error rate for single bases
	maxRunLength   int                // Maximum run length to model
	rng            *rand.Rand         // Random number generator
	
	// Platform-specific parameters
	deletionBias   float64            // Bias toward deletion errors
	insertionBias  float64            // Bias toward insertion errors
	qualityDecay   float64            // Quality decay per homopolymer length
	
	// Error distribution parameters
	errorDistribution map[string]float64 // Error type -> probability
	runLengthWeights  []float64          // Weight for each run length
	
	// Training data
	observedErrors   []HomopolymerError
	modelTrained     bool
}

// HomopolymerError represents an observed homopolymer error
type HomopolymerError struct {
	Position      int
	RunLength     int
	Base          byte
	ErrorType     string  // "insertion", "deletion", "substitution"
	Quality       float64
	Context       string  // Surrounding sequence context
	Platform      string
	TrueLength    int     // True homopolymer length
	ObservedLength int    // Observed homopolymer length
}

// HomopolymerRun represents a homopolymer run in a sequence
type HomopolymerRun struct {
	Start      int
	End        int
	Length     int
	Base       byte
	Context    string
	ErrorProne bool
	ErrorRate  float64
}

// ErrorSimulationResult represents the result of error simulation
type ErrorSimulationResult struct {
	OriginalSequence   []byte
	ErrorSequence      []byte
	OriginalQualities  []byte
	ErrorQualities     []byte
	Errors             []HomopolymerError
	HomopolymerRuns    []HomopolymerRun
	ErrorStatistics    ErrorStatistics
}

// ErrorStatistics contains statistics about introduced errors
type ErrorStatistics struct {
	TotalErrors       int
	InsertionErrors   int
	DeletionErrors    int
	SubstitutionErrors int
	AvgRunLength      float64
	MaxRunLength      int
	ErrorRate         float64
}

// NewHomopolymerErrorModel creates a new homopolymer error model
func NewHomopolymerErrorModel(platform string) *HomopolymerErrorModel {
	model := &HomopolymerErrorModel{
		errorRates:      make(map[int]float64),
		platform:        platform,
		errorTypes:      []string{"insertion", "deletion", "substitution"},
		contextEffects:  true,
		baseErrorRate:   0.001,
		maxRunLength:    20,
		rng:             rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
		observedErrors:  make([]HomopolymerError, 0),
		modelTrained:    false,
	}

	// Set platform-specific parameters
	model.setPlatformParameters()
	
	// Initialize error rates
	model.initializeErrorRates()
	
	return model
}

// setPlatformParameters sets platform-specific parameters
func (m *HomopolymerErrorModel) setPlatformParameters() {
	switch m.platform {
	case "illumina", "novaseq", "nextseq", "miseq":
		// Illumina has relatively low homopolymer error rates
		m.deletionBias = 0.3
		m.insertionBias = 0.2
		m.qualityDecay = 1.0
		m.maxRunLength = 15
		m.baseErrorRate = 0.001
		
	case "pacbio", "pacbio_revio":
		// PacBio has higher homopolymer error rates, mainly deletions
		m.deletionBias = 0.6
		m.insertionBias = 0.3
		m.qualityDecay = 2.0
		m.maxRunLength = 25
		m.baseErrorRate = 0.01
		
	case "nanopore", "promethion":
		// Nanopore has highest homopolymer error rates
		m.deletionBias = 0.7
		m.insertionBias = 0.25
		m.qualityDecay = 3.0
		m.maxRunLength = 30
		m.baseErrorRate = 0.05
		
	default:
		// Default to Illumina-like behavior
		m.deletionBias = 0.3
		m.insertionBias = 0.2
		m.qualityDecay = 1.0
		m.maxRunLength = 15
		m.baseErrorRate = 0.001
	}
	
	// Set error distribution
	m.errorDistribution = map[string]float64{
		"deletion":     m.deletionBias,
		"insertion":    m.insertionBias,
		"substitution": 1.0 - m.deletionBias - m.insertionBias,
	}
}

// initializeErrorRates initializes error rates based on homopolymer length
func (m *HomopolymerErrorModel) initializeErrorRates() {
	// Error rate increases exponentially with homopolymer length
	for length := 1; length <= m.maxRunLength; length++ {
		var errorRate float64
		
		if length == 1 {
			errorRate = m.baseErrorRate
		} else {
			// Exponential increase with run length
			switch m.platform {
			case "illumina", "novaseq", "nextseq", "miseq":
				errorRate = m.baseErrorRate * math.Pow(1.5, float64(length-1))
			case "pacbio", "pacbio_revio":
				errorRate = m.baseErrorRate * math.Pow(2.0, float64(length-1))
			case "nanopore", "promethion":
				errorRate = m.baseErrorRate * math.Pow(2.5, float64(length-1))
			default:
				errorRate = m.baseErrorRate * math.Pow(1.5, float64(length-1))
			}
		}
		
		// Cap at maximum error rate
		if errorRate > 0.5 {
			errorRate = 0.5
		}
		
		m.errorRates[length] = errorRate
	}
}

// FindHomopolymerRuns finds homopolymer runs in a sequence
func (m *HomopolymerErrorModel) FindHomopolymerRuns(sequence []byte) []HomopolymerRun {
	if len(sequence) == 0 {
		return []HomopolymerRun{}
	}
	
	runs := make([]HomopolymerRun, 0)
	currentRun := HomopolymerRun{
		Start: 0,
		Base:  sequence[0],
	}
	
	for i := 1; i < len(sequence); i++ {
		if sequence[i] == currentRun.Base {
			// Continue current run
			continue
		} else {
			// End current run
			currentRun.End = i
			currentRun.Length = i - currentRun.Start
			
			// Add context
			currentRun.Context = m.getContext(sequence, currentRun.Start, currentRun.End)
			
			// Determine if error-prone
			currentRun.ErrorProne = currentRun.Length >= 2
			currentRun.ErrorRate = m.getErrorRate(currentRun.Length)
			
			// Add run if it's a homopolymer (length >= 2)
			if currentRun.Length >= 2 {
				runs = append(runs, currentRun)
			}
			
			// Start new run
			currentRun = HomopolymerRun{
				Start: i,
				Base:  sequence[i],
			}
		}
	}
	
	// Handle last run
	currentRun.End = len(sequence)
	currentRun.Length = len(sequence) - currentRun.Start
	currentRun.Context = m.getContext(sequence, currentRun.Start, currentRun.End)
	currentRun.ErrorProne = currentRun.Length >= 2
	currentRun.ErrorRate = m.getErrorRate(currentRun.Length)
	
	if currentRun.Length >= 2 {
		runs = append(runs, currentRun)
	}
	
	return runs
}

// getContext gets the sequence context around a homopolymer run
func (m *HomopolymerErrorModel) getContext(sequence []byte, start, end int) string {
	contextSize := 3
	contextStart := start - contextSize
	contextEnd := end + contextSize
	
	if contextStart < 0 {
		contextStart = 0
	}
	if contextEnd > len(sequence) {
		contextEnd = len(sequence)
	}
	
	return string(sequence[contextStart:contextEnd])
}

// getErrorRate gets the error rate for a given homopolymer length
func (m *HomopolymerErrorModel) getErrorRate(length int) float64 {
	if length > m.maxRunLength {
		length = m.maxRunLength
	}
	
	if rate, exists := m.errorRates[length]; exists {
		return rate
	}
	
	return m.baseErrorRate
}

// SimulateErrors simulates homopolymer errors in a sequence
func (m *HomopolymerErrorModel) SimulateErrors(sequence []byte, qualities []byte) (*ErrorSimulationResult, error) {
	if len(sequence) == 0 {
		return nil, fmt.Errorf("sequence cannot be empty")
	}
	
	if len(qualities) != len(sequence) {
		return nil, fmt.Errorf("sequence and quality lengths must match")
	}
	
	result := &ErrorSimulationResult{
		OriginalSequence:  make([]byte, len(sequence)),
		ErrorSequence:     make([]byte, len(sequence)),
		OriginalQualities: make([]byte, len(qualities)),
		ErrorQualities:    make([]byte, len(qualities)),
		Errors:            make([]HomopolymerError, 0),
		HomopolymerRuns:   make([]HomopolymerRun, 0),
	}
	
	copy(result.OriginalSequence, sequence)
	copy(result.ErrorSequence, sequence)
	copy(result.OriginalQualities, qualities)
	copy(result.ErrorQualities, qualities)
	
	// Find homopolymer runs
	runs := m.FindHomopolymerRuns(sequence)
	result.HomopolymerRuns = runs
	
	// Apply errors to each run. Runs are found left to right on the original
	// sequence, so an edit lands at its original position plus the net length
	// change of the edits already applied (shift).
	shift := 0
	for _, run := range runs {
		if m.shouldIntroduceError(run) {
			error := m.introduceError(result, run, shift)
			if error != nil {
				result.Errors = append(result.Errors, *error)
				shift += error.ObservedLength - error.TrueLength
			}
		}
	}
	
	// Calculate statistics
	result.ErrorStatistics = m.calculateErrorStatistics(result)
	
	return result, nil
}

// shouldIntroduceError determines if an error should be introduced
func (m *HomopolymerErrorModel) shouldIntroduceError(run HomopolymerRun) bool {
	if !run.ErrorProne {
		return false
	}
	
	// Base probability on error rate
	baseProb := run.ErrorRate
	
	// Adjust for context effects
	if m.contextEffects {
		contextMod := m.getContextModifier(run.Context)
		baseProb *= contextMod
	}
	
	return m.rng.Float64() < baseProb
}

// getContextModifier calculates context-dependent error rate modifier
func (m *HomopolymerErrorModel) getContextModifier(context string) float64 {
	modifier := 1.0
	
	// Check for specific sequence contexts that increase error rates
	if strings.Contains(context, "CG") {
		modifier *= 1.2 // CG context increases error rate
	}
	
	if strings.Contains(context, "AT") {
		modifier *= 0.9 // AT context decreases error rate
	}
	
	// Check for tandem repeats
	if len(context) >= 4 {
		for i := 0; i < len(context)-3; i++ {
			if context[i] == context[i+2] && context[i+1] == context[i+3] {
				modifier *= 1.3 // Tandem repeat increases error rate
				break
			}
		}
	}
	
	return modifier
}

// introduceError introduces an error into the sequence
// introduceError applies one error inside run. error.Position is in original
// sequence coordinates; shift maps it into the current ErrorSequence.
func (m *HomopolymerErrorModel) introduceError(result *ErrorSimulationResult, run HomopolymerRun, shift int) *HomopolymerError {
	// Choose error type
	errorType := m.chooseErrorType()
	
	// Choose position within run
	errorPos := run.Start + m.rng.IntN(run.Length)
	
	// Create error object
	error := &HomopolymerError{
		Position:       errorPos,
		RunLength:      run.Length,
		Base:           run.Base,
		ErrorType:      errorType,
		Quality:        float64(result.OriginalQualities[errorPos]),
		Context:        run.Context,
		Platform:       m.platform,
		TrueLength:     run.Length,
		ObservedLength: run.Length,
	}
	
	// Apply error based on type, at the position in the edited sequence
	at := errorPos + shift
	switch errorType {
	case "deletion":
		m.applyDeletion(result, error, at)
	case "insertion":
		m.applyInsertion(result, error, at)
	case "substitution":
		m.applySubstitution(result, error, at)
	}
	
	// Adjust quality scores
	m.adjustQualityScores(result, at)
	
	return error
}

// chooseErrorType chooses an error type based on platform distribution.
// Types are visited in sorted order (not Go's randomised map order) so the
// choice is a deterministic function of the seed.
func (m *HomopolymerErrorModel) chooseErrorType() string {
	r := m.rng.Float64()
	cumulative := 0.0

	errorTypes := make([]string, 0, len(m.errorDistribution))
	for errorType := range m.errorDistribution {
		errorTypes = append(errorTypes, errorType)
	}
	sort.Strings(errorTypes)
	for _, errorType := range errorTypes {
		cumulative += m.errorDistribution[errorType]
		if r <= cumulative {
			return errorType
		}
	}
	
	return "deletion" // Default
}

// applyDeletion applies a deletion error
func (m *HomopolymerErrorModel) applyDeletion(result *ErrorSimulationResult, error *HomopolymerError, at int) {
	// Remove base at error position
	newSequence := make([]byte, 0, len(result.ErrorSequence)-1)
	newQualities := make([]byte, 0, len(result.ErrorQualities)-1)
	
	for i := 0; i < len(result.ErrorSequence); i++ {
		if i != at {
			newSequence = append(newSequence, result.ErrorSequence[i])
			newQualities = append(newQualities, result.ErrorQualities[i])
		}
	}
	
	result.ErrorSequence = newSequence
	result.ErrorQualities = newQualities
	error.ObservedLength = error.TrueLength - 1
}

// applyInsertion applies an insertion error
func (m *HomopolymerErrorModel) applyInsertion(result *ErrorSimulationResult, error *HomopolymerError, at int) {
	// Insert base at error position
	newSequence := make([]byte, 0, len(result.ErrorSequence)+1)
	newQualities := make([]byte, 0, len(result.ErrorQualities)+1)
	
	for i := 0; i < len(result.ErrorSequence); i++ {
		newSequence = append(newSequence, result.ErrorSequence[i])
		newQualities = append(newQualities, result.ErrorQualities[i])
		
		if i == at {
			// Insert same base (homopolymer extension)
			newSequence = append(newSequence, error.Base)
			// Quality of inserted base is typically lower
			insertQuality := float64(result.ErrorQualities[i]) - 5.0
			if insertQuality < 0 {
				insertQuality = 0
			}
			newQualities = append(newQualities, byte(insertQuality))
		}
	}
	
	result.ErrorSequence = newSequence
	result.ErrorQualities = newQualities
	error.ObservedLength = error.TrueLength + 1
}

// applySubstitution applies a substitution error
func (m *HomopolymerErrorModel) applySubstitution(result *ErrorSimulationResult, error *HomopolymerError, at int) {
	// Choose different base
	bases := []byte{'A', 'T', 'G', 'C'}
	var newBase byte
	
	for _, base := range bases {
		if base != error.Base {
			newBase = base
			break
		}
	}
	
	// Apply substitution
	result.ErrorSequence[at] = newBase
	error.ObservedLength = error.TrueLength // Length doesn't change
}

// adjustQualityScores adjusts quality scores based on error introduction
func (m *HomopolymerErrorModel) adjustQualityScores(result *ErrorSimulationResult, at int) {
	// Reduce quality scores around error position
	reductionAmount := 3.0
	windowSize := 2
	
	start := at - windowSize
	end := at + windowSize
	
	if start < 0 {
		start = 0
	}
	if end >= len(result.ErrorQualities) {
		end = len(result.ErrorQualities) - 1
	}
	
	for i := start; i <= end; i++ {
		originalQuality := float64(result.ErrorQualities[i])
		newQuality := originalQuality - reductionAmount
		
		if newQuality < 0 {
			newQuality = 0
		}
		
		result.ErrorQualities[i] = byte(newQuality)
	}
}

// calculateErrorStatistics calculates error statistics
func (m *HomopolymerErrorModel) calculateErrorStatistics(result *ErrorSimulationResult) ErrorStatistics {
	stats := ErrorStatistics{
		TotalErrors: len(result.Errors),
	}
	
	if len(result.Errors) == 0 {
		return stats
	}
	
	// Count error types
	runLengthSum := 0
	maxRunLength := 0
	
	for _, error := range result.Errors {
		switch error.ErrorType {
		case "insertion":
			stats.InsertionErrors++
		case "deletion":
			stats.DeletionErrors++
		case "substitution":
			stats.SubstitutionErrors++
		}
		
		runLengthSum += error.RunLength
		if error.RunLength > maxRunLength {
			maxRunLength = error.RunLength
		}
	}
	
	if len(result.Errors) > 0 {
		stats.AvgRunLength = float64(runLengthSum) / float64(len(result.Errors))
	}
	stats.MaxRunLength = maxRunLength
	if len(result.OriginalSequence) > 0 {
		stats.ErrorRate = float64(len(result.Errors)) / float64(len(result.OriginalSequence))
	}
	
	return stats
}

// TrainFromObservedErrors trains the error-type distribution from observed error records.
// Because the input contains only errors (not total run observations), it is impossible to
// compute a meaningful error rate (errors / observations). The per-length errorRates set by
// initializeErrorRates() are therefore preserved unchanged. Only errorDistribution (the
// fraction of each error type) is updated, as that can be correctly computed from errors alone.
// To update errorRates when total run counts are available, use TrainFromRunObservations.
func (m *HomopolymerErrorModel) TrainFromObservedErrors(errors []HomopolymerError) error {
	if len(errors) == 0 {
		return fmt.Errorf("no training data provided")
	}

	m.observedErrors = errors

	// Update error type distribution from observed errors.
	typeCounts := make(map[string]int)
	for _, e := range errors {
		typeCounts[e.ErrorType]++
	}

	total := len(errors)
	for errorType, count := range typeCounts {
		m.errorDistribution[errorType] = float64(count) / float64(total)
	}

	m.modelTrained = true
	return nil
}

// TrainFromRunObservations updates per-length error rates when both observed errors and
// total run observation counts are known. totalRunCounts maps homopolymer run length to
// the number of times that run length was observed in the sequenced data (errors + correct).
func (m *HomopolymerErrorModel) TrainFromRunObservations(errors []HomopolymerError, totalRunCounts map[int]int) error {
	// Update error type distribution first.
	if err := m.TrainFromObservedErrors(errors); err != nil {
		return err
	}

	// Tally errors per run length.
	runLengthErrors := make(map[int]int)
	for _, e := range errors {
		runLengthErrors[e.RunLength]++
	}

	// Update error rates: errors / total observations for that run length.
	for length, total := range totalRunCounts {
		if total > 0 {
			rate := float64(runLengthErrors[length]) / float64(total)
			// Cap at maximum error rate to maintain physical plausibility.
			if rate > 0.5 {
				rate = 0.5
			}
			m.errorRates[length] = rate
		}
	}

	return nil
}

// SetSeed sets the random seed for reproducible error simulation
func (m *HomopolymerErrorModel) SetSeed(seed int64) {
	m.rng = rand.New(rand.NewPCG(uint64(seed), 0)) // #nosec G115 G404
}

// GetModelStats returns statistics about the model
func (m *HomopolymerErrorModel) GetModelStats() map[string]interface{} {
	stats := map[string]interface{}{
		"platform":         m.platform,
		"model_trained":    m.modelTrained,
		"base_error_rate":  m.baseErrorRate,
		"max_run_length":   m.maxRunLength,
		"deletion_bias":    m.deletionBias,
		"insertion_bias":   m.insertionBias,
		"quality_decay":    m.qualityDecay,
		"context_effects":  m.contextEffects,
		"observed_errors":  len(m.observedErrors),
	}
	
	// Add error rate statistics
	if len(m.errorRates) > 0 {
		maxErrorRate := 0.0
		minErrorRate := 1.0
		
		for _, rate := range m.errorRates {
			if rate > maxErrorRate {
				maxErrorRate = rate
			}
			if rate < minErrorRate {
				minErrorRate = rate
			}
		}
		
		stats["max_error_rate"] = maxErrorRate
		stats["min_error_rate"] = minErrorRate
		stats["error_rate_entries"] = len(m.errorRates)
	}
	
	return stats
}

// ValidateModel validates the model configuration
func (m *HomopolymerErrorModel) ValidateModel() error {
	if m.baseErrorRate < 0 || m.baseErrorRate > 1 {
		return fmt.Errorf("invalid base error rate: %f", m.baseErrorRate)
	}
	
	if m.maxRunLength <= 0 {
		return fmt.Errorf("invalid max run length: %d", m.maxRunLength)
	}
	
	if m.deletionBias < 0 || m.insertionBias < 0 || m.deletionBias+m.insertionBias > 1 {
		return fmt.Errorf("invalid error bias values: deletion=%f, insertion=%f", m.deletionBias, m.insertionBias)
	}
	
	// Validate error rates
	for length, rate := range m.errorRates {
		if rate < 0 || rate > 1 {
			return fmt.Errorf("invalid error rate for length %d: %f", length, rate)
		}
	}
	
	return nil
}

// GetPlatformParameters returns platform-specific parameters
func (m *HomopolymerErrorModel) GetPlatformParameters() map[string]interface{} {
	return map[string]interface{}{
		"platform":        m.platform,
		"deletion_bias":   m.deletionBias,
		"insertion_bias":  m.insertionBias,
		"quality_decay":   m.qualityDecay,
		"base_error_rate": m.baseErrorRate,
		"max_run_length":  m.maxRunLength,
	}
}

// EstimateErrorImpact estimates the impact of homopolymer errors on sequence
func (m *HomopolymerErrorModel) EstimateErrorImpact(sequence []byte) map[string]interface{} {
	runs := m.FindHomopolymerRuns(sequence)
	
	totalRuns := len(runs)
	errorProneRuns := 0
	totalErrorProbability := 0.0
	avgRunLength := 0.0
	maxRunLength := 0
	
	for _, run := range runs {
		if run.ErrorProne {
			errorProneRuns++
			totalErrorProbability += run.ErrorRate
		}
		
		avgRunLength += float64(run.Length)
		if run.Length > maxRunLength {
			maxRunLength = run.Length
		}
	}
	
	errorProneFraction := 0.0
	if totalRuns > 0 {
		avgRunLength /= float64(totalRuns)
		errorProneFraction = float64(errorProneRuns) / float64(totalRuns)
	}
	homopolymerDensity := 0.0
	if len(sequence) > 0 {
		homopolymerDensity = float64(totalRuns) / float64(len(sequence))
	}

	return map[string]interface{}{
		"total_homopolymer_runs":    totalRuns,
		"error_prone_runs":          errorProneRuns,
		"expected_errors":           totalErrorProbability,
		"avg_run_length":            avgRunLength,
		"max_run_length":            maxRunLength,
		"error_prone_fraction":      errorProneFraction,
		"homopolymer_density":       homopolymerDensity,
	}
}