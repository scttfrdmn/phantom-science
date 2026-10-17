package validation

import (
	"fmt"
	"math"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewExpertTestingFramework tests framework creation
func TestNewExpertTestingFramework(t *testing.T) {
	framework := NewExpertTestingFramework()

	assert.NotNil(t, framework)
	assert.NotNil(t, framework.testCases)
	assert.NotNil(t, framework.expertInterface)
	assert.NotNil(t, framework.resultAggregator)
	assert.NotNil(t, framework.performanceTracker)
	assert.NotNil(t, framework.discriminator)
	assert.NotNil(t, framework.rng)
	assert.Len(t, framework.testCases, 0)
}

// TestNewExpertInterface tests expert interface creation
func TestNewExpertInterface(t *testing.T) {
	expertInterface := NewExpertInterface()

	assert.NotNil(t, expertInterface)
	assert.NotNil(t, expertInterface.experts)
	assert.NotNil(t, expertInterface.testSessions)
	assert.NotNil(t, expertInterface.calibrationTests)
	assert.NotNil(t, expertInterface.performanceHistory)
	assert.Len(t, expertInterface.experts, 0)
}

// TestNewResultAggregator tests result aggregator creation
func TestNewResultAggregator(t *testing.T) {
	aggregator := NewResultAggregator()

	assert.NotNil(t, aggregator)
	assert.Equal(t, "weighted", aggregator.aggregationMethod)
	assert.NotNil(t, aggregator.expertWeights)
	assert.NotNil(t, aggregator.consensusResults)
	assert.NotNil(t, aggregator.disagreementAnalysis)
	assert.Len(t, aggregator.consensusResults, 0)
}

// TestNewPerformanceTracker tests performance tracker creation
func TestNewPerformanceTracker(t *testing.T) {
	tracker := NewPerformanceTracker()

	assert.NotNil(t, tracker)
	assert.NotNil(t, tracker.accuracyTrends)
	assert.NotNil(t, tracker.expertMetrics)
	assert.NotNil(t, tracker.improvementAreas)
	assert.Len(t, tracker.accuracyTrends, 0)
}

// TestNewSyntheticDataDiscriminator tests discriminator creation
func TestNewSyntheticDataDiscriminator(t *testing.T) {
	discriminator := NewSyntheticDataDiscriminator()

	assert.NotNil(t, discriminator)
	assert.NotNil(t, discriminator.features)
	assert.NotNil(t, discriminator.trainingData)
	assert.NotNil(t, discriminator.featureImportance)
	assert.Len(t, discriminator.features, 0)
}

// TestGenerateTestCases tests test case generation
func TestGenerateTestCases(t *testing.T) {
	framework := NewExpertTestingFramework()

	// Create mock data
	syntheticData := [][]byte{
		[]byte("ATCGATCGATCGATCGATCGATCG"),
		[]byte("GCGCGCGCGCGCGCGCGCGCGCGC"),
		[]byte("ATATATATATATATATATATATA"),
	}

	realData := [][]byte{
		[]byte("ATCGNNNATCGATCGATCGNNNTC"),
		[]byte("GCATTACGATCGATCGACGTATGC"),
		[]byte("NNNATCGATCGATCGATCGATNNN"),
	}

	numCases := 6

	err := framework.GenerateTestCases(numCases, syntheticData, realData)
	require.NoError(t, err)

	// Check test cases were generated
	assert.Len(t, framework.testCases, numCases)

	// Check balance between synthetic and real
	syntheticCount := 0
	realCount := 0
	for _, testCase := range framework.testCases {
		if testCase.CorrectAnswer == "synthetic" {
			syntheticCount++
		} else if testCase.CorrectAnswer == "real" {
			realCount++
		}
	}

	assert.Equal(t, 3, syntheticCount)
	assert.Equal(t, 3, realCount)

	// Check test case properties
	for _, testCase := range framework.testCases {
		assert.NotEmpty(t, testCase.ID)
		assert.Equal(t, "fastq", testCase.DataType)
		assert.NotNil(t, testCase.Metadata)
		assert.Contains(t, []string{"synthetic", "real"}, testCase.CorrectAnswer)
		assert.Contains(t, []string{"easy", "medium", "hard"}, testCase.DifficultyLevel)
		
		// Check metadata
		assert.Contains(t, testCase.Metadata, "sequence_length")
		assert.Contains(t, testCase.Metadata, "gc_content")
	}
}

// TestGenerateTestCases_ErrorHandling tests error handling in test case generation
func TestGenerateTestCases_ErrorHandling(t *testing.T) {
	framework := NewExpertTestingFramework()

	// Test with empty synthetic data
	err := framework.GenerateTestCases(10, [][]byte{}, [][]byte{[]byte("ATCG")})
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "both synthetic and real data are required")

	// Test with empty real data
	err = framework.GenerateTestCases(10, [][]byte{[]byte("ATCG")}, [][]byte{})
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "both synthetic and real data are required")

	// Test with non-positive number of cases
	err = framework.GenerateTestCases(0, [][]byte{[]byte("ATCG")}, [][]byte{[]byte("ATCG")})
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "number of test cases must be positive")
}

// TestAssessDifficultyLevel tests difficulty level assessment
func TestAssessDifficultyLevel(t *testing.T) {
	framework := NewExpertTestingFramework()

	testCases := []struct {
		sequence         []byte
		expectedLevel    string
		description      string
	}{
		{[]byte("AAAAAAAAAAAAAAAAAAAAAA"), "easy", "low GC content"},
		{[]byte("GGGGGGGGGGGGGGGGGGGGGG"), "easy", "high GC content"},
		{[]byte("AAAAAAAAAAAAAAAAAAAAAA"), "easy", "low complexity"},
		{[]byte("ATCGATCGATCGATCGATCGATCG"), "easy", "balanced GC but low complexity (repetitive)"},
		{[]byte("ATGCGATCGCTAGCTAGCATGCATGCTAGC"), "medium", "balanced GC and high complexity"},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			level := framework.assessDifficultyLevel(tc.sequence)
			assert.Equal(t, tc.expectedLevel, level)
		})
	}
}

// TestCalculateGCContent tests GC content calculation
func TestCalculateGCContent(t *testing.T) {
	framework := NewExpertTestingFramework()

	testCases := []struct {
		sequence []byte
		expected float64
		description string
	}{
		{[]byte(""), 0.0, "empty sequence"},
		{[]byte("ATCG"), 0.5, "balanced sequence"},
		{[]byte("GGCC"), 1.0, "all GC"},
		{[]byte("AATT"), 0.0, "all AT"},
		{[]byte("ATCGATCG"), 0.5, "mixed sequence"},
		{[]byte("NNNNATCG"), 0.25, "sequence with Ns - GC calculated from all positions"},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			gcContent := framework.calculateGCContent(tc.sequence)
			assert.InDelta(t, tc.expected, gcContent, 0.001)
		})
	}
}

// TestCalculateComplexity tests sequence complexity calculation
func TestCalculateComplexity(t *testing.T) {
	framework := NewExpertTestingFramework()

	testCases := []struct {
		sequence []byte
		minComplexity float64
		maxComplexity float64
		description string
	}{
		{[]byte(""), 0.0, 0.0, "empty sequence"},
		{[]byte("AT"), 0.0, 0.0, "too short"},
		{[]byte("AAAAAAA"), 0.0, 0.2, "low complexity"},
		{[]byte("ATCGATCG"), 0.4, 1.0, "high complexity"},
		{[]byte("ATGCATGCATGC"), 0.2, 0.8, "medium complexity"},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			complexity := framework.calculateComplexity(tc.sequence)
			assert.GreaterOrEqual(t, complexity, tc.minComplexity)
			assert.LessOrEqual(t, complexity, tc.maxComplexity)
		})
	}
}

// TestRegisterExpert tests expert registration
func TestRegisterExpert(t *testing.T) {
	framework := NewExpertTestingFramework()

	expert := &Expert{
		ID:                "expert1",
		Name:              "Dr. Smith",
		Institution:       "University",
		ExpertiseLevel:    "senior",
		SpecializationAreas: []string{"genomics", "bioinformatics"},
		CalibrationScore:  0.85,
		PerformanceHistory: make([]ExpertPerformance, 0),
		IsActive:          true,
	}

	framework.RegisterExpert(expert)

	// Check expert was registered
	registeredExpert, exists := framework.expertInterface.experts[expert.ID]
	assert.True(t, exists)
	assert.Equal(t, expert, registeredExpert)

	// Check performance history was initialized
	history, exists := framework.expertInterface.performanceHistory[expert.ID]
	assert.True(t, exists)
	assert.Len(t, history, 0)
}

// TestCollectExpertFeedback tests expert feedback collection
func TestCollectExpertFeedback(t *testing.T) {
	framework := NewExpertTestingFramework()

	// Register an expert
	expert := &Expert{
		ID:             "expert1",
		Name:           "Dr. Smith",
		ExpertiseLevel: "senior",
		PerformanceHistory: make([]ExpertPerformance, 0),
	}
	framework.RegisterExpert(expert)

	// Create a test case
	testCase := BlindTestCase{
		ID:            "test1",
		CorrectAnswer: "synthetic",
		DifficultyLevel: "medium",
	}
	framework.testCases = append(framework.testCases, testCase)

	// Collect feedback
	err := framework.CollectExpertFeedback(
		"test1", "expert1", "synthetic", 0.8, "Pattern suggests artificial generation")
	require.NoError(t, err)

	// Check test case was updated
	updatedTestCase := framework.testCases[0]
	assert.Equal(t, "expert1", updatedTestCase.ExpertID)
	assert.Equal(t, "synthetic", updatedTestCase.Classification)
	assert.Equal(t, 0.8, updatedTestCase.Confidence)
	assert.Equal(t, "Pattern suggests artificial generation", updatedTestCase.Justification)
	assert.True(t, updatedTestCase.IsCorrect)

	// Check expert performance was recorded
	history := framework.expertInterface.performanceHistory["expert1"]
	assert.Len(t, history, 1)
	assert.Equal(t, 1.0, history[0].Accuracy)
	assert.Equal(t, 0.8, history[0].Confidence)
}

// TestCollectExpertFeedback_ErrorHandling tests error handling in feedback collection
func TestCollectExpertFeedback_ErrorHandling(t *testing.T) {
	framework := NewExpertTestingFramework()

	// Test with nonexistent test case
	err := framework.CollectExpertFeedback("nonexistent", "expert1", "synthetic", 0.8, "test")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "test case not found")

	// Test with nonexistent expert
	testCase := BlindTestCase{ID: "test1", CorrectAnswer: "synthetic"}
	framework.testCases = append(framework.testCases, testCase)

	err = framework.CollectExpertFeedback("test1", "nonexistent", "synthetic", 0.8, "test")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "expert not found")
}

// TestCalculateConsensus tests consensus calculation
func TestCalculateConsensus(t *testing.T) {
	framework := NewExpertTestingFramework()

	// Create test cases with expert responses
	testCases := []BlindTestCase{
		{ID: "test1", Classification: "synthetic", Confidence: 0.8, ExpertID: "expert1"},
		{ID: "test1", Classification: "synthetic", Confidence: 0.9, ExpertID: "expert2"},
		{ID: "test1", Classification: "real", Confidence: 0.6, ExpertID: "expert3"},
	}

	consensus := framework.calculateConsensus(testCases)

	assert.Equal(t, "test1", consensus.TestCaseID)
	assert.Equal(t, "synthetic", consensus.ConsensusAnswer)
	assert.InDelta(t, 0.77, consensus.Confidence, 0.01) // (0.8 + 0.9 + 0.6) / 3
	assert.InDelta(t, 0.67, consensus.AgreementLevel, 0.01) // 2/3
	assert.Len(t, consensus.ParticipatingExperts, 3)
	assert.Greater(t, consensus.StatisticalSignificance, 0.0)
}

// TestCalculateStatisticalSignificance tests statistical significance calculation
func TestCalculateStatisticalSignificance(t *testing.T) {
	framework := NewExpertTestingFramework()

	testCases := []struct {
		agreementLevel float64
		sampleSize     int
		description    string
	}{
		{0.5, 10, "random agreement"},
		{0.8, 10, "high agreement"},
		{0.9, 20, "very high agreement large sample"},
		{0.0, 0, "zero sample size"},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			significance := framework.calculateStatisticalSignificance(tc.agreementLevel, tc.sampleSize)
			
			if tc.sampleSize == 0 {
				assert.Equal(t, 0.0, significance)
			} else {
				assert.GreaterOrEqual(t, significance, 0.0)
				assert.LessOrEqual(t, significance, 1.0)
			}
		})
	}
}

// TestIdentifyDisagreementCauses tests disagreement cause identification
func TestIdentifyDisagreementCauses(t *testing.T) {
	framework := NewExpertTestingFramework()

	// Register experts with different expertise levels
	framework.RegisterExpert(&Expert{ID: "expert1", ExpertiseLevel: "junior"})
	framework.RegisterExpert(&Expert{ID: "expert2", ExpertiseLevel: "senior"})
	framework.RegisterExpert(&Expert{ID: "expert3", ExpertiseLevel: "expert"})

	testCases := []BlindTestCase{
		{ExpertID: "expert1", Confidence: 0.3, DifficultyLevel: "hard"},
		{ExpertID: "expert2", Confidence: 0.8, DifficultyLevel: "hard"},
		{ExpertID: "expert3", Confidence: 0.9, DifficultyLevel: "hard"},
	}

	causes := framework.identifyDisagreementCauses(testCases)

	assert.Greater(t, len(causes), 0)
	assert.Contains(t, causes, "high_difficulty")
	assert.Contains(t, causes, "expertise_level_difference")
	assert.Contains(t, causes, "confidence_uncertainty")
}

// TestCalculateStandardDeviation tests standard deviation calculation
func TestCalculateStandardDeviation(t *testing.T) {
	framework := NewExpertTestingFramework()

	testCases := []struct {
		values   []float64
		expected float64
		description string
	}{
		{[]float64{}, 0.0, "empty slice"},
		{[]float64{1.0}, 0.0, "single value"},
		{[]float64{1.0, 1.0, 1.0}, 0.0, "identical values"},
		{[]float64{1.0, 2.0, 3.0}, math.Sqrt(1.0), "simple case"},
		{[]float64{0.3, 0.8, 0.9}, 0.321, "confidence values"},
	}

	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			stdDev := framework.calculateStandardDeviation(tc.values)
			assert.InDelta(t, tc.expected, stdDev, 0.01)
		})
	}
}

// TestGetAuthenticityScore tests authenticity score calculation
func TestGetAuthenticityScore(t *testing.T) {
	framework := NewExpertTestingFramework()

	// Create test cases and consensus results
	framework.testCases = []BlindTestCase{
		{ID: "test1", CorrectAnswer: "synthetic"},
		{ID: "test2", CorrectAnswer: "real"},
		{ID: "test3", CorrectAnswer: "synthetic"},
	}

	framework.resultAggregator.consensusResults = []ConsensusResult{
		{TestCaseID: "test1", ConsensusAnswer: "synthetic"}, // correct
		{TestCaseID: "test2", ConsensusAnswer: "synthetic"}, // incorrect
		{TestCaseID: "test3", ConsensusAnswer: "real"},      // incorrect
	}

	score := framework.GetAuthenticityScore()

	// Expert accuracy: 1/3 = 0.33
	// Authenticity score: (1 - 0.33) * 100 = 67
	assert.InDelta(t, 66.67, score, 0.1)
}

// TestGetFrameworkStats tests framework statistics
func TestGetFrameworkStats(t *testing.T) {
	framework := NewExpertTestingFramework()

	// Add some data
	framework.testCases = []BlindTestCase{
		{ID: "test1", DifficultyLevel: "easy", ExpertID: "expert1"},
		{ID: "test2", DifficultyLevel: "medium", ExpertID: "expert1"},
		{ID: "test3", DifficultyLevel: "hard", ExpertID: ""},
	}

	framework.RegisterExpert(&Expert{ID: "expert1"})
	framework.RegisterExpert(&Expert{ID: "expert2"})

	framework.expertInterface.performanceHistory["expert1"] = []ExpertPerformance{
		{Accuracy: 0.8},
		{Accuracy: 0.9},
	}

	stats := framework.GetFrameworkStats()

	assert.Equal(t, 3, stats["total_test_cases"])
	assert.Equal(t, 2, stats["registered_experts"])
	assert.Equal(t, 2, stats["completed_evaluations"])

	// Check difficulty breakdown
	difficultyBreakdown := stats["difficulty_breakdown"].(map[string]int)
	assert.Equal(t, 1, difficultyBreakdown["easy"])
	assert.Equal(t, 1, difficultyBreakdown["medium"])
	assert.Equal(t, 1, difficultyBreakdown["hard"])

	// Check expert performance
	expertPerformance := stats["expert_performance"].(map[string]float64)
	assert.InDelta(t, 0.85, expertPerformance["expert1"], 0.01) // (0.8 + 0.9) / 2
}

// TestValidateFramework tests framework validation
func TestValidateFramework(t *testing.T) {
	framework := NewExpertTestingFramework()

	// Test empty framework
	err := framework.ValidateFramework()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "no test cases available")

	// Add test cases but no experts
	framework.testCases = []BlindTestCase{
		{CorrectAnswer: "synthetic"},
		{CorrectAnswer: "real"},
	}

	err = framework.ValidateFramework()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "no experts registered")

	// Add experts
	framework.RegisterExpert(&Expert{ID: "expert1"})

	err = framework.ValidateFramework()
	assert.NoError(t, err)

	// Test imbalanced test cases
	framework.testCases = []BlindTestCase{
		{CorrectAnswer: "synthetic"},
		{CorrectAnswer: "synthetic"},
	}

	err = framework.ValidateFramework()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "test cases must include both synthetic and real examples")
}

// TestExpert tests expert structure
func TestExpert(t *testing.T) {
	expert := &Expert{
		ID:              "expert1",
		Name:            "Dr. Smith",
		Institution:     "University",
		ExpertiseLevel:  "senior",
		SpecializationAreas: []string{"genomics", "bioinformatics"},
		CalibrationScore: 0.85,
		PerformanceHistory: make([]ExpertPerformance, 0),
		BiasProfile: BiasProfile{
			SyntheticBias:         0.1,
			ConfidenceBias:        -0.05,
			SpeedAccuracyTradeoff: 0.3,
			ConsistencyScore:      0.9,
			LearningRate:          0.2,
		},
		IsActive: true,
	}

	assert.Equal(t, "expert1", expert.ID)
	assert.Equal(t, "Dr. Smith", expert.Name)
	assert.Equal(t, "senior", expert.ExpertiseLevel)
	assert.Contains(t, expert.SpecializationAreas, "genomics")
	assert.Equal(t, 0.85, expert.CalibrationScore)
	assert.True(t, expert.IsActive)
	assert.Equal(t, 0.1, expert.BiasProfile.SyntheticBias)
}

// TestBlindTestCase tests blind test case structure
func TestBlindTestCase(t *testing.T) {
	testCase := BlindTestCase{
		ID:               "test1",
		SyntheticData:    []byte("ATCGATCG"),
		RealData:         nil,
		DataType:         "fastq",
		Metadata:         map[string]string{"platform": "illumina"},
		ExpertID:         "expert1",
		ResponseTime:     time.Second * 30,
		Classification:   "synthetic",
		Confidence:       0.8,
		Justification:    "Pattern analysis",
		CorrectAnswer:    "synthetic",
		IsCorrect:        true,
		DifficultyLevel:  "medium",
		GenerationMethod: "hmm",
	}

	assert.Equal(t, "test1", testCase.ID)
	assert.Equal(t, "fastq", testCase.DataType)
	assert.Equal(t, "expert1", testCase.ExpertID)
	assert.Equal(t, "synthetic", testCase.Classification)
	assert.Equal(t, 0.8, testCase.Confidence)
	assert.True(t, testCase.IsCorrect)
	assert.Equal(t, "medium", testCase.DifficultyLevel)
}

// TestNewStatisticalValidator tests statistical validator creation
func TestNewStatisticalValidator(t *testing.T) {
	validator := NewStatisticalValidator()

	assert.NotNil(t, validator)
	assert.NotNil(t, validator.jsCalculator)
	assert.NotNil(t, validator.testHistory)
	assert.NotNil(t, validator.cache)

	// Check default thresholds
	assert.Equal(t, 0.01, validator.thresholds.JSdivergenceThreshold)
	assert.Equal(t, 0.05, validator.thresholds.PValueThreshold)
	assert.Equal(t, 100, validator.thresholds.SampleSizeMinimum)
	assert.Equal(t, 0.1, validator.thresholds.EffectSizeThreshold)
}

// TestValidationThresholds tests validation thresholds structure
func TestValidationThresholds(t *testing.T) {
	thresholds := ValidationThresholds{
		JSdivergenceThreshold:   0.01,
		PValueThreshold:        0.05,
		SampleSizeMinimum:      100,
		EffectSizeThreshold:    0.1,
	}

	assert.Equal(t, 0.01, thresholds.JSdivergenceThreshold)
	assert.Equal(t, 0.05, thresholds.PValueThreshold)
	assert.Equal(t, 100, thresholds.SampleSizeMinimum)
	assert.Equal(t, 0.1, thresholds.EffectSizeThreshold)
}

// TestExpertTestingResults tests expert testing results structure
func TestExpertTestingResults(t *testing.T) {
	results := &ExpertTestingResults{
		OverallAccuracy:      0.48,
		AuthenticityScore:    52.0,
		ParticipatingExperts: 5,
		TestCasesEvaluated:   100,
		ConfidenceStatistics: ConfidenceStatistics{
			MeanConfidence:   0.75,
			MedianConfidence: 0.80,
			StdDevConfidence: 0.15,
			CalibrationError: 0.12,
		},
		PerformanceByDifficulty: map[string]float64{
			"easy":   0.65,
			"medium": 0.48,
			"hard":   0.32,
		},
		BiasAnalysis: BiasAnalysis{
			SyntheticBias:      0.05,
			RealBias:          -0.03,
			OverconfidenceBias: 0.08,
			ConsistencyIndex:   0.72,
		},
	}

	assert.Equal(t, 0.48, results.OverallAccuracy)
	assert.Equal(t, 52.0, results.AuthenticityScore)
	assert.Equal(t, 5, results.ParticipatingExperts)
	assert.Equal(t, 100, results.TestCasesEvaluated)
	
	// Check confidence statistics
	assert.Equal(t, 0.75, results.ConfidenceStatistics.MeanConfidence)
	assert.Equal(t, 0.80, results.ConfidenceStatistics.MedianConfidence)
	
	// Check performance by difficulty
	assert.Equal(t, 0.65, results.PerformanceByDifficulty["easy"])
	assert.Equal(t, 0.48, results.PerformanceByDifficulty["medium"])
	assert.Equal(t, 0.32, results.PerformanceByDifficulty["hard"])
	
	// Check bias analysis
	assert.Equal(t, 0.05, results.BiasAnalysis.SyntheticBias)
	assert.Equal(t, -0.03, results.BiasAnalysis.RealBias)
}

// Benchmark tests

// BenchmarkGenerateTestCases benchmarks test case generation
func BenchmarkGenerateTestCases(b *testing.B) {
	framework := NewExpertTestingFramework()

	syntheticData := make([][]byte, 100)
	realData := make([][]byte, 100)

	for i := 0; i < 100; i++ {
		syntheticData[i] = []byte(fmt.Sprintf("ATCGATCGATCGATCGATCG%d", i))
		realData[i] = []byte(fmt.Sprintf("GCATTACGATCGATCGACGT%d", i))
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		err := framework.GenerateTestCases(50, syntheticData, realData)
		if err != nil {
			b.Fatal(err)
		}
		framework.testCases = framework.testCases[:0] // Reset for next iteration
	}
}

// BenchmarkCalculateComplexity benchmarks complexity calculation
func BenchmarkCalculateComplexity(b *testing.B) {
	framework := NewExpertTestingFramework()
	sequence := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		framework.calculateComplexity(sequence)
	}
}

// BenchmarkCollectExpertFeedback benchmarks feedback collection
func BenchmarkCollectExpertFeedback(b *testing.B) {
	framework := NewExpertTestingFramework()

	// Setup
	expert := &Expert{ID: "expert1", PerformanceHistory: make([]ExpertPerformance, 0)}
	framework.RegisterExpert(expert)

	testCase := BlindTestCase{ID: "test1", CorrectAnswer: "synthetic", DifficultyLevel: "medium"}
	framework.testCases = append(framework.testCases, testCase)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		err := framework.CollectExpertFeedback("test1", "expert1", "synthetic", 0.8, "test")
		if err != nil {
			b.Fatal(err)
		}
	}
}

// Mock structures for testing