// Package validation holds frameworks for checking synthetic data against
// real data: Jensen-Shannon divergence between distributions, blind
// expert-review sessions (test-case generation, response collection and
// majority-vote consensus), and record-keeping for academic collaboration and
// publication workflows.
//
// The expert, academic and publication parts manage in-memory records and
// generate text from templates; they do not run a study, contact a journal or
// submit anything. Network clients (PubMed, ORCID, CrossRef) are used only
// where noted on the individual methods.
package validation

import (
	"fmt"
	"math"
	"math/rand/v2"
	"time"
)

// ExpertTestingFramework manages expert validation and blind testing
type ExpertTestingFramework struct {
	testCases          []BlindTestCase
	expertInterface    *ExpertInterface
	resultAggregator   *ResultAggregator
	performanceTracker *PerformanceTracker
	discriminator      *SyntheticDataDiscriminator
	rng                *rand.Rand
}

// BlindTestCase represents a single blind test case
type BlindTestCase struct {
	ID               string
	SyntheticData    []byte
	RealData         []byte
	DataType         string // "fastq", "fasta", "quality_scores"
	Metadata         map[string]string
	ExpertID         string
	ResponseTime     time.Duration
	Classification   string  // "synthetic", "real", "unsure"
	Confidence       float64 // 0.0 to 1.0
	Justification    string
	CorrectAnswer    string
	IsCorrect        bool
	DifficultyLevel  string // "easy", "medium", "hard"
	GenerationMethod string // "random", "hmm", "population_genetics"
}

// ExpertInterface manages expert interactions
type ExpertInterface struct {
	experts            map[string]*Expert
	testSessions       map[string]*TestSession
	calibrationTests   []CalibrationTest
	performanceHistory map[string][]ExpertPerformance
}

// Expert represents an expert evaluator
type Expert struct {
	ID                  string
	Name                string
	Email               string // recipient address for notifications
	Institution         string
	ExpertiseLevel      string // "junior", "senior", "expert"
	SpecializationAreas []string
	CalibrationScore    float64
	PerformanceHistory  []ExpertPerformance
	BiasProfile         BiasProfile
	IsActive            bool
}

// TestSession represents a testing session
type TestSession struct {
	ID             string
	ExpertID       string
	StartTime      time.Time
	EndTime        time.Time
	TestCases      []string // Test case IDs
	CompletedCases []string
	Performance    SessionPerformance
	Interruptions  []Interruption
	Environment    TestEnvironment
}

// ExpertPerformance tracks expert performance over time
type ExpertPerformance struct {
	Date             time.Time
	TestCaseID       string
	Accuracy         float64
	ResponseTime     time.Duration
	Confidence       float64
	DifficultyLevel  string
	CalibrationError float64 // |confidence - accuracy|
}

// BiasProfile captures systematic biases in expert evaluation
type BiasProfile struct {
	SyntheticBias         float64 // Tendency to classify as synthetic
	ConfidenceBias        float64 // Tendency to be over/under confident
	SpeedAccuracyTradeoff float64 // Trade-off between speed and accuracy
	ConsistencyScore      float64 // How consistent are their judgments
	LearningRate          float64 // How quickly they adapt to new patterns
}

// ResultAggregator aggregates results from multiple experts
type ResultAggregator struct {
	aggregationMethod    string // "majority", "weighted", "bayesian"
	expertWeights        map[string]float64
	consensusResults     []ConsensusResult
	disagreementAnalysis []DisagreementAnalysis
}

// ConsensusResult represents aggregated expert consensus
type ConsensusResult struct {
	TestCaseID              string
	ConsensusAnswer         string
	Confidence              float64
	AgreementLevel          float64
	ParticipatingExperts    []string
	DissentingOpinions      []string
	StatisticalSignificance float64
}

// DisagreementAnalysis analyzes cases where experts disagree
type DisagreementAnalysis struct {
	TestCaseID          string
	DisagreementType    string // "classification", "confidence", "justification"
	ConflictingAnswers  []string
	PossibleCauses      []string
	DataCharacteristics map[string]float64
	ResolutionMethod    string
}

// PerformanceTracker tracks framework performance over time
type PerformanceTracker struct {
	overallAccuracy  float64
	accuracyTrends   []AccuracyTrend
	expertMetrics    map[string]ExpertMetrics
	systemMetrics    SystemMetrics
	improvementAreas []ImprovementArea
}

// AccuracyTrend tracks accuracy changes over time
type AccuracyTrend struct {
	Date            time.Time
	OverallAccuracy float64
	ByDifficulty    map[string]float64
	ByDataType      map[string]float64
	ByGenMethod     map[string]float64
	SampleSize      int
}

// SyntheticDataDiscriminator is a placeholder for an ML classifier of
// synthetic versus real data. No model or feature extractor is set and
// nothing in this package calls it.
type SyntheticDataDiscriminator struct {
	model             MachineLearningModel
	features          []FeatureExtractor
	trainingData      []LabeledExample
	performance       DiscriminatorPerformance
	featureImportance map[string]float64
}

// FeatureExtractor turns raw data into a feature vector for a discriminator.
type FeatureExtractor interface {
	ExtractFeatures(data []byte) []float64
	GetFeatureNames() []string
}

// MachineLearningModel is the classifier interface a discriminator would use.
type MachineLearningModel interface {
	Train(features [][]float64, labels []int) error
	Predict(features []float64) (float64, error)
	GetAccuracy() float64
	GetFeatureImportance() map[string]float64
}

// NewExpertTestingFramework creates a new expert testing framework
func NewExpertTestingFramework() *ExpertTestingFramework {
	return &ExpertTestingFramework{
		testCases:          make([]BlindTestCase, 0),
		expertInterface:    NewExpertInterface(),
		resultAggregator:   NewResultAggregator(),
		performanceTracker: NewPerformanceTracker(),
		discriminator:      NewSyntheticDataDiscriminator(),
		rng:                rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
	}
}

// NewExpertInterface creates a new expert interface
func NewExpertInterface() *ExpertInterface {
	return &ExpertInterface{
		experts:            make(map[string]*Expert),
		testSessions:       make(map[string]*TestSession),
		calibrationTests:   make([]CalibrationTest, 0),
		performanceHistory: make(map[string][]ExpertPerformance),
	}
}

// NewResultAggregator creates a new result aggregator
func NewResultAggregator() *ResultAggregator {
	return &ResultAggregator{
		aggregationMethod:    "weighted",
		expertWeights:        make(map[string]float64),
		consensusResults:     make([]ConsensusResult, 0),
		disagreementAnalysis: make([]DisagreementAnalysis, 0),
	}
}

// NewPerformanceTracker creates a new performance tracker
func NewPerformanceTracker() *PerformanceTracker {
	return &PerformanceTracker{
		accuracyTrends:   make([]AccuracyTrend, 0),
		expertMetrics:    make(map[string]ExpertMetrics),
		improvementAreas: make([]ImprovementArea, 0),
	}
}

// NewSyntheticDataDiscriminator creates a new discriminator
func NewSyntheticDataDiscriminator() *SyntheticDataDiscriminator {
	return &SyntheticDataDiscriminator{
		features:          make([]FeatureExtractor, 0),
		trainingData:      make([]LabeledExample, 0),
		featureImportance: make(map[string]float64),
	}
}

// GenerateTestCases generates blind test cases for expert evaluation
func (f *ExpertTestingFramework) GenerateTestCases(numCases int, syntheticData, realData [][]byte) error {
	if len(syntheticData) == 0 || len(realData) == 0 {
		return fmt.Errorf("both synthetic and real data are required")
	}

	if numCases <= 0 {
		return fmt.Errorf("number of test cases must be positive")
	}

	// Generate balanced test cases
	syntheticCases := numCases / 2
	realCases := numCases - syntheticCases

	// Generate synthetic test cases
	for i := 0; i < syntheticCases; i++ {
		if i >= len(syntheticData) {
			break
		}

		testCase := BlindTestCase{
			ID:               fmt.Sprintf("synthetic_%d", i),
			SyntheticData:    syntheticData[i],
			RealData:         nil,
			DataType:         "fastq",
			Metadata:         make(map[string]string),
			CorrectAnswer:    "synthetic",
			DifficultyLevel:  f.assessDifficultyLevel(syntheticData[i]),
			GenerationMethod: "hmm", // Would be determined by actual generation method
		}

		// Add metadata
		testCase.Metadata["sequence_length"] = fmt.Sprintf("%d", len(syntheticData[i]))
		testCase.Metadata["gc_content"] = fmt.Sprintf("%.3f", f.calculateGCContent(syntheticData[i]))

		f.testCases = append(f.testCases, testCase)
	}

	// Generate real test cases
	for i := 0; i < realCases; i++ {
		if i >= len(realData) {
			break
		}

		testCase := BlindTestCase{
			ID:               fmt.Sprintf("real_%d", i),
			SyntheticData:    nil,
			RealData:         realData[i],
			DataType:         "fastq",
			Metadata:         make(map[string]string),
			CorrectAnswer:    "real",
			DifficultyLevel:  f.assessDifficultyLevel(realData[i]),
			GenerationMethod: "real",
		}

		// Add metadata
		testCase.Metadata["sequence_length"] = fmt.Sprintf("%d", len(realData[i]))
		testCase.Metadata["gc_content"] = fmt.Sprintf("%.3f", f.calculateGCContent(realData[i]))

		f.testCases = append(f.testCases, testCase)
	}

	// Shuffle test cases to randomize order
	f.shuffleTestCases()

	return nil
}

// assessDifficultyLevel assesses the difficulty level of distinguishing a sequence
func (f *ExpertTestingFramework) assessDifficultyLevel(data []byte) string {
	// Simple heuristic based on sequence characteristics
	gcContent := f.calculateGCContent(data)
	complexity := f.calculateComplexity(data)

	// Easy: Very high/low GC content or low complexity
	if gcContent < 0.3 || gcContent > 0.7 || complexity < 0.5 {
		return "easy"
	}

	// Hard: Balanced GC content and high complexity
	if gcContent >= 0.45 && gcContent <= 0.55 && complexity > 0.8 {
		return "hard"
	}

	// Medium: Everything else
	return "medium"
}

// calculateGCContent calculates GC content of a sequence
func (f *ExpertTestingFramework) calculateGCContent(data []byte) float64 {
	if len(data) == 0 {
		return 0.0
	}

	gcCount := 0
	for _, base := range data {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}

	return float64(gcCount) / float64(len(data))
}

// calculateComplexity calculates sequence complexity
func (f *ExpertTestingFramework) calculateComplexity(data []byte) float64 {
	if len(data) < 3 {
		return 0.0
	}

	// Count unique k-mers (k=3)
	kmerCount := make(map[string]int)
	for i := 0; i <= len(data)-3; i++ {
		kmer := string(data[i : i+3])
		kmerCount[kmer]++
	}

	// Calculate complexity as fraction of unique k-mers
	maxPossibleKmers := len(data) - 2
	return float64(len(kmerCount)) / float64(maxPossibleKmers)
}

// shuffleTestCases randomly shuffles test cases
func (f *ExpertTestingFramework) shuffleTestCases() {
	for i := len(f.testCases) - 1; i > 0; i-- {
		j := f.rng.IntN(i + 1)
		f.testCases[i], f.testCases[j] = f.testCases[j], f.testCases[i]
	}
}

// RegisterExpert registers a new expert evaluator
func (f *ExpertTestingFramework) RegisterExpert(expert *Expert) {
	f.expertInterface.experts[expert.ID] = expert
	f.expertInterface.performanceHistory[expert.ID] = make([]ExpertPerformance, 0)
}

// CollectExpertFeedback collects feedback from an expert
func (f *ExpertTestingFramework) CollectExpertFeedback(testCaseID, expertID string, classification string, confidence float64, justification string) error {
	// Find test case
	var testCase *BlindTestCase
	for i := range f.testCases {
		if f.testCases[i].ID == testCaseID {
			testCase = &f.testCases[i]
			break
		}
	}

	if testCase == nil {
		return fmt.Errorf("test case not found: %s", testCaseID)
	}

	// Find expert
	expert, exists := f.expertInterface.experts[expertID]
	if !exists {
		return fmt.Errorf("expert not found: %s", expertID)
	}

	// Record response
	testCase.ExpertID = expertID
	testCase.Classification = classification
	testCase.Confidence = confidence
	testCase.Justification = justification
	testCase.IsCorrect = classification == testCase.CorrectAnswer

	// Record performance
	performance := ExpertPerformance{
		Date:             time.Now(),
		TestCaseID:       testCaseID,
		Accuracy:         0.0,
		ResponseTime:     testCase.ResponseTime,
		Confidence:       confidence,
		DifficultyLevel:  testCase.DifficultyLevel,
		CalibrationError: math.Abs(confidence - f.boolToFloat(testCase.IsCorrect)),
	}

	if testCase.IsCorrect {
		performance.Accuracy = 1.0
	}

	expert.PerformanceHistory = append(expert.PerformanceHistory, performance)
	f.expertInterface.performanceHistory[expertID] = append(f.expertInterface.performanceHistory[expertID], performance)

	return nil
}

// boolToFloat converts boolean to float64
func (f *ExpertTestingFramework) boolToFloat(b bool) float64 {
	if b {
		return 1.0
	}
	return 0.0
}

// AggregateResults aggregates results from multiple experts
func (f *ExpertTestingFramework) AggregateResults() error {
	// Group test cases by ID
	testCaseGroups := make(map[string][]BlindTestCase)
	for _, testCase := range f.testCases {
		testCaseGroups[testCase.ID] = append(testCaseGroups[testCase.ID], testCase)
	}

	// Aggregate each test case
	for testCaseID, testCases := range testCaseGroups {
		if len(testCases) == 0 {
			continue
		}

		consensus := f.calculateConsensus(testCases)
		f.resultAggregator.consensusResults = append(f.resultAggregator.consensusResults, consensus)

		// Analyze disagreements
		if consensus.AgreementLevel < 0.8 {
			disagreement := f.analyzeDisagreement(testCaseID, testCases)
			f.resultAggregator.disagreementAnalysis = append(f.resultAggregator.disagreementAnalysis, disagreement)
		}
	}

	return nil
}

// calculateConsensus calculates consensus from expert responses
func (f *ExpertTestingFramework) calculateConsensus(testCases []BlindTestCase) ConsensusResult {
	if len(testCases) == 0 {
		return ConsensusResult{}
	}

	// Count votes
	votes := make(map[string]int)
	confidenceSum := 0.0
	participatingExperts := make([]string, 0)

	for _, testCase := range testCases {
		votes[testCase.Classification]++
		confidenceSum += testCase.Confidence
		participatingExperts = append(participatingExperts, testCase.ExpertID)
	}

	// Find majority answer
	maxVotes := 0
	consensusAnswer := ""
	for answer, count := range votes {
		if count > maxVotes {
			maxVotes = count
			consensusAnswer = answer
		}
	}

	// Calculate agreement level
	agreementLevel := float64(maxVotes) / float64(len(testCases))

	// Calculate average confidence
	avgConfidence := confidenceSum / float64(len(testCases))

	return ConsensusResult{
		TestCaseID:              testCases[0].ID,
		ConsensusAnswer:         consensusAnswer,
		Confidence:              avgConfidence,
		AgreementLevel:          agreementLevel,
		ParticipatingExperts:    participatingExperts,
		StatisticalSignificance: f.calculateStatisticalSignificance(agreementLevel, len(testCases)),
	}
}

// calculateStatisticalSignificance calculates statistical significance of consensus
func (f *ExpertTestingFramework) calculateStatisticalSignificance(agreementLevel float64, sampleSize int) float64 {
	// Two-sided p-value for agreement vs 0.5, normal approximation to the
	// binomial
	if sampleSize == 0 {
		return 0.0
	}

	// Z-score for binomial test
	p := 0.5 // null hypothesis probability
	z := (agreementLevel - p) / math.Sqrt(p*(1-p)/float64(sampleSize))

	return math.Erfc(math.Abs(z) / math.Sqrt(2))
}

// analyzeDisagreement analyzes disagreement patterns
func (f *ExpertTestingFramework) analyzeDisagreement(testCaseID string, testCases []BlindTestCase) DisagreementAnalysis {
	conflictingAnswers := make([]string, 0)
	answerCounts := make(map[string]int)

	for _, testCase := range testCases {
		answerCounts[testCase.Classification]++
	}

	for answer := range answerCounts {
		conflictingAnswers = append(conflictingAnswers, answer)
	}

	// Analyze possible causes
	possibleCauses := f.identifyDisagreementCauses(testCases)

	return DisagreementAnalysis{
		TestCaseID:         testCaseID,
		DisagreementType:   "classification",
		ConflictingAnswers: conflictingAnswers,
		PossibleCauses:     possibleCauses,
		ResolutionMethod:   "majority_vote",
	}
}

// identifyDisagreementCauses identifies possible causes of disagreement
func (f *ExpertTestingFramework) identifyDisagreementCauses(testCases []BlindTestCase) []string {
	causes := make([]string, 0)

	// Check if difficulty level might be a factor
	if len(testCases) > 0 && testCases[0].DifficultyLevel == "hard" {
		causes = append(causes, "high_difficulty")
	}

	// Check for expertise level differences
	expertiseLevels := make(map[string]int)
	for _, testCase := range testCases {
		if expert, exists := f.expertInterface.experts[testCase.ExpertID]; exists {
			expertiseLevels[expert.ExpertiseLevel]++
		}
	}

	if len(expertiseLevels) > 1 {
		causes = append(causes, "expertise_level_difference")
	}

	// Check for confidence spread
	confidences := make([]float64, 0)
	for _, testCase := range testCases {
		confidences = append(confidences, testCase.Confidence)
	}

	if f.calculateStandardDeviation(confidences) > 0.3 {
		causes = append(causes, "confidence_uncertainty")
	}

	return causes
}

// calculateStandardDeviation calculates standard deviation of a slice
func (f *ExpertTestingFramework) calculateStandardDeviation(values []float64) float64 {
	if len(values) <= 1 {
		return 0.0
	}

	// Calculate mean
	sum := 0.0
	for _, value := range values {
		sum += value
	}
	mean := sum / float64(len(values))

	// Calculate variance
	variance := 0.0
	for _, value := range values {
		variance += math.Pow(value-mean, 2)
	}
	variance /= float64(len(values) - 1)

	return math.Sqrt(variance)
}

// GetAuthenticityScore returns (1 - consensus accuracy) * 100 over the
// aggregated test cases. 50 corresponds to chance-level discrimination; values
// above 50 mean the consensus was more often wrong than right, which is not
// evidence of greater authenticity. With no consensus results it returns 0.
func (f *ExpertTestingFramework) GetAuthenticityScore() float64 {
	if len(f.resultAggregator.consensusResults) == 0 {
		return 0.0
	}

	correctCount := 0
	totalCount := 0

	for _, consensus := range f.resultAggregator.consensusResults {
		totalCount++

		// Find original test case to get correct answer
		for _, testCase := range f.testCases {
			if testCase.ID == consensus.TestCaseID {
				if consensus.ConsensusAnswer == testCase.CorrectAnswer {
					correctCount++
				}
				break
			}
		}
	}

	if totalCount == 0 {
		return 0.0
	}

	// Convert to authenticity score (lower expert accuracy = higher authenticity)
	expertAccuracy := float64(correctCount) / float64(totalCount)
	authenticityScore := (1.0 - expertAccuracy) * 100.0

	return authenticityScore
}

// GetFrameworkStats returns comprehensive framework statistics
func (f *ExpertTestingFramework) GetFrameworkStats() map[string]interface{} {
	stats := map[string]interface{}{
		"total_test_cases":           len(f.testCases),
		"registered_experts":         len(f.expertInterface.experts),
		"completed_evaluations":      f.countCompletedEvaluations(),
		"consensus_results":          len(f.resultAggregator.consensusResults),
		"disagreement_cases":         len(f.resultAggregator.disagreementAnalysis),
		"overall_authenticity_score": f.GetAuthenticityScore(),
	}

	// Add difficulty breakdown
	difficultyBreakdown := make(map[string]int)
	for _, testCase := range f.testCases {
		difficultyBreakdown[testCase.DifficultyLevel]++
	}
	stats["difficulty_breakdown"] = difficultyBreakdown

	// Add expert performance summary
	expertPerformance := make(map[string]float64)
	for expertID, history := range f.expertInterface.performanceHistory {
		if len(history) > 0 {
			accuracy := 0.0
			for _, perf := range history {
				accuracy += perf.Accuracy
			}
			expertPerformance[expertID] = accuracy / float64(len(history))
		}
	}
	stats["expert_performance"] = expertPerformance

	return stats
}

// countCompletedEvaluations counts completed evaluations
func (f *ExpertTestingFramework) countCompletedEvaluations() int {
	count := 0
	for _, testCase := range f.testCases {
		if testCase.ExpertID != "" {
			count++
		}
	}
	return count
}

// ValidateFramework validates the framework configuration
func (f *ExpertTestingFramework) ValidateFramework() error {
	if len(f.testCases) == 0 {
		return fmt.Errorf("no test cases available")
	}

	if len(f.expertInterface.experts) == 0 {
		return fmt.Errorf("no experts registered")
	}

	// Check for balanced test cases
	syntheticCount := 0
	realCount := 0

	for _, testCase := range f.testCases {
		if testCase.CorrectAnswer == "synthetic" {
			syntheticCount++
		} else if testCase.CorrectAnswer == "real" {
			realCount++
		}
	}

	if syntheticCount == 0 || realCount == 0 {
		return fmt.Errorf("test cases must include both synthetic and real examples")
	}

	return nil
}

// Additional types and structures for comprehensive expert testing

type CalibrationTest struct {
	ID           string
	TestCaseID   string
	KnownAnswer  string
	ExpertAnswer string
	Confidence   float64
	IsCorrect    bool
}

type TestEnvironment struct {
	Location     string
	Lighting     string
	NoiseLevel   string
	Distractions []string
}

type Interruption struct {
	Time        time.Time
	Type        string
	Duration    time.Duration
	Description string
}

type SessionPerformance struct {
	Accuracy         float64
	AvgResponseTime  time.Duration
	AvgConfidence    float64
	ConsistencyScore float64
	FatigueIndex     float64
}

type ExpertMetrics struct {
	OverallAccuracy   float64
	CalibrationScore  float64
	ConsistencyScore  float64
	ResponseTimeStats ResponseTimeStats
	SpecialtyAccuracy map[string]float64
}

type ResponseTimeStats struct {
	Mean   time.Duration
	Median time.Duration
	StdDev time.Duration
	Min    time.Duration
	Max    time.Duration
}

type SystemMetrics struct {
	AvgAuthenticityScore float64
	ScoreVariability     float64
	ExpertAgreementRate  float64
	SystemReliability    float64
}

type ImprovementArea struct {
	Area           string
	Priority       string
	Description    string
	Recommendation string
	ExpectedImpact float64
}

type LabeledExample struct {
	Features []float64
	Label    int // 0 = real, 1 = synthetic
	Source   string
}

type DiscriminatorPerformance struct {
	Accuracy        float64
	Precision       float64
	Recall          float64
	F1Score         float64
	ROCAUCScore     float64
	ConfusionMatrix [][]int
}

// ExpertTestingResults contains overall results from expert testing
type ExpertTestingResults struct {
	OverallAccuracy         float64
	AuthenticityScore       float64
	ExpertConsensus         []ConsensusResult
	DisagreementCases       []DisagreementAnalysis
	ParticipatingExperts    int
	TestCasesEvaluated      int
	ConfidenceStatistics    ConfidenceStatistics
	PerformanceByDifficulty map[string]float64
	BiasAnalysis            BiasAnalysis
}

// ConfidenceStatistics contains statistics about expert confidence
type ConfidenceStatistics struct {
	MeanConfidence   float64
	MedianConfidence float64
	StdDevConfidence float64
	CalibrationError float64
}

// BiasAnalysis contains analysis of systematic biases
type BiasAnalysis struct {
	SyntheticBias      float64 // Tendency to classify as synthetic
	RealBias           float64 // Tendency to classify as real
	OverconfidenceBias float64 // Tendency to be overconfident
	ConsistencyIndex   float64 // How consistent experts are
}

// StatisticalValidator compares two distributions by Jensen-Shannon
// divergence against ValidationThresholds.JSdivergenceThreshold (0.01) and
// keeps a history of the comparisons. It runs no other test.
type StatisticalValidator struct {
	jsCalculator *JSdivergenceCalculator
	thresholds   ValidationThresholds
	testHistory  []StatisticalTest
	cache        map[string]*ValidationResult
}

// ValidationThresholds contains thresholds for statistical tests
type ValidationThresholds struct {
	JSdivergenceThreshold float64
	PValueThreshold       float64 // not used yet
	SampleSizeMinimum     int     // not used yet
	EffectSizeThreshold   float64 // not used yet
}

// StatisticalTest represents a single statistical test
type StatisticalTest struct {
	TestType   string
	TestDate   time.Time
	Platform   string
	Population string
	Result     *ValidationResult
	DataSize   int
	Passed     bool
}

// NewStatisticalValidator creates a new statistical validator
func NewStatisticalValidator() *StatisticalValidator {
	return &StatisticalValidator{
		jsCalculator: NewJSdivergenceCalculator(),
		thresholds: ValidationThresholds{
			JSdivergenceThreshold: 0.01,
			PValueThreshold:       0.05,
			SampleSizeMinimum:     100,
			EffectSizeThreshold:   0.1,
		},
		testHistory: make([]StatisticalTest, 0),
		cache:       make(map[string]*ValidationResult),
	}
}

// ValidateDistributions computes the JS divergence between the two
// distributions and marks the comparison passed when it is at or below the
// threshold.
func (sv *StatisticalValidator) ValidateDistributions(syntheticData, realData map[string]float64, testType string) (*ValidationResult, error) {
	// Use JS divergence calculator
	result, err := sv.jsCalculator.Calculate(syntheticData, realData)
	if err != nil {
		return nil, fmt.Errorf("failed to calculate JS divergence: %w", err)
	}

	// Check against thresholds
	result.IsSignificant = result.JSdivergence > sv.thresholds.JSdivergenceThreshold

	// Add recommendations
	if result.IsSignificant {
		result.Recommendations = append(result.Recommendations,
			"JS divergence exceeds threshold - consider model tuning")
	}

	// Record test
	test := StatisticalTest{
		TestType: testType,
		TestDate: time.Now(),
		Result:   result,
		DataSize: len(syntheticData) + len(realData),
		Passed:   !result.IsSignificant,
	}
	sv.testHistory = append(sv.testHistory, test)

	return result, nil
}

// GetValidationStats returns validation statistics
func (sv *StatisticalValidator) GetValidationStats() map[string]interface{} {
	stats := map[string]interface{}{
		"total_tests":       len(sv.testHistory),
		"tests_passed":      sv.countPassedTests(),
		"avg_js_divergence": sv.calculateAverageJSDivergence(),
		"thresholds":        sv.thresholds,
	}

	// Add test breakdown by type
	testBreakdown := make(map[string]int)
	for _, test := range sv.testHistory {
		testBreakdown[test.TestType]++
	}
	stats["test_breakdown"] = testBreakdown

	return stats
}

// countPassedTests counts tests that passed
func (sv *StatisticalValidator) countPassedTests() int {
	count := 0
	for _, test := range sv.testHistory {
		if test.Passed {
			count++
		}
	}
	return count
}

// calculateAverageJSDivergence calculates average JS divergence
func (sv *StatisticalValidator) calculateAverageJSDivergence() float64 {
	if len(sv.testHistory) == 0 {
		return 0.0
	}

	sum := 0.0
	for _, test := range sv.testHistory {
		if test.Result != nil {
			sum += test.Result.JSdivergence
		}
	}

	return sum / float64(len(sv.testHistory))
}
