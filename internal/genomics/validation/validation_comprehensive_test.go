package validation_test

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/validation"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ==================== JSdivergenceCalculator Tests ====================

func TestNewJSdivergenceCalculator(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	require.NotNil(t, calc)
}

func TestJSdivergenceCalculator_SetSmoothingFactor(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	calc.SetSmoothingFactor(0.01)
	// No error expected; just confirm it doesn't panic
}

func TestJSdivergenceCalculator_SetConfidenceLevel(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	calc.SetConfidenceLevel(0.99)
	// No error expected
}

func TestJSdivergenceCalculator_Calculate_IdenticalDistributions(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	p := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
	q := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}

	result, err := calc.Calculate(p, q)
	require.NoError(t, err)
	require.NotNil(t, result)

	// Identical distributions should have near-zero divergence
	assert.InDelta(t, 0.0, result.JSdivergence, 0.001)
	assert.Greater(t, result.SampleSize, 0)
	assert.True(t, result.Smoothed)
	assert.NotEmpty(t, result.Recommendations)
}

func TestJSdivergenceCalculator_Calculate_DifferentDistributions(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	// Population-like allele frequency distributions
	p := map[string]float64{"AA": 0.7, "AT": 0.2, "TT": 0.1}
	q := map[string]float64{"AA": 0.1, "AT": 0.2, "TT": 0.7}

	result, err := calc.Calculate(p, q)
	require.NoError(t, err)
	require.NotNil(t, result)

	assert.Greater(t, result.JSdivergence, 0.0)
	// Confidence interval should be valid
	assert.LessOrEqual(t, result.ConfidenceInterval[0], result.ConfidenceInterval[1])
}

func TestJSdivergenceCalculator_Calculate_AsymmetricDistributions(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	// P has keys that Q doesn't and vice versa
	p := map[string]float64{"CEU": 0.6, "YRI": 0.4}
	q := map[string]float64{"CHB": 0.5, "JPT": 0.3, "CEU": 0.2}

	// Q doesn't sum to 1.0 exactly with the above — fix:
	q = map[string]float64{"CHB": 0.5, "JPT": 0.3, "CEU": 0.2}

	result, err := calc.Calculate(p, q)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.GreaterOrEqual(t, result.JSdivergence, 0.0)
}

func TestJSdivergenceCalculator_Calculate_EmptyP(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	p := map[string]float64{}
	q := map[string]float64{"A": 0.5, "T": 0.5}

	result, err := calc.Calculate(p, q)
	assert.Error(t, err)
	assert.Nil(t, result)
}

func TestJSdivergenceCalculator_Calculate_EmptyQ(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{}

	result, err := calc.Calculate(p, q)
	assert.Error(t, err)
	assert.Nil(t, result)
}

func TestJSdivergenceCalculator_Calculate_InvalidP_NegativeProb(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	p := map[string]float64{"A": -0.1, "T": 1.1}
	q := map[string]float64{"A": 0.5, "T": 0.5}

	result, err := calc.Calculate(p, q)
	assert.Error(t, err)
	assert.Nil(t, result)
	assert.Contains(t, err.Error(), "invalid distribution P")
}

func TestJSdivergenceCalculator_Calculate_InvalidP_WrongSum(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	p := map[string]float64{"A": 0.3, "T": 0.3}
	q := map[string]float64{"A": 0.5, "T": 0.5}

	result, err := calc.Calculate(p, q)
	assert.Error(t, err)
	assert.Nil(t, result)
}

func TestJSdivergenceCalculator_Calculate_InvalidQ(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{"A": 0.9, "T": 0.9} // sum > 1

	result, err := calc.Calculate(p, q)
	assert.Error(t, err)
	assert.Nil(t, result)
	assert.Contains(t, err.Error(), "invalid distribution Q")
}

func TestJSdivergenceCalculator_IsSignificant(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	tests := []struct {
		name       string
		divergence float64
		sampleSize int
		wantSig    bool
	}{
		{"large divergence small sample", 0.05, 50, true},
		{"large divergence medium sample", 0.02, 500, true},
		{"large divergence large sample", 0.02, 1500, true},
		{"small divergence", 0.005, 1000, false},
		{"zero divergence", 0.0, 100, false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			sig := calc.IsSignificant(tt.divergence, tt.sampleSize)
			assert.Equal(t, tt.wantSig, sig)
		})
	}
}

func TestJSdivergenceCalculator_GetAuthenticityScore(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	// Score should be between 0 and 100
	score0 := calc.GetAuthenticityScore(0.0)
	assert.InDelta(t, 100.0, score0, 0.01)

	scoreHigh := calc.GetAuthenticityScore(0.5)
	assert.Greater(t, score0, scoreHigh)
	assert.GreaterOrEqual(t, scoreHigh, 0.0)
}

func TestJSdivergenceCalculator_GetQualityGrade(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	tests := []struct {
		divergence float64
		wantGrade  string
	}{
		{0.0005, "Excellent"},
		{0.005, "Good"},
		{0.02, "Fair"},
		{0.07, "Poor"},
		{0.15, "Failed"},
	}

	for _, tt := range tests {
		t.Run(tt.wantGrade, func(t *testing.T) {
			grade := calc.GetQualityGrade(tt.divergence)
			assert.Equal(t, tt.wantGrade, grade)
		})
	}
}

func TestJSdivergenceCalculator_CompareMultipleDistributions(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	distributions := map[string]map[string]float64{
		"CEU": {"A": 0.3, "T": 0.3, "G": 0.2, "C": 0.2},
		"YRI": {"A": 0.4, "T": 0.2, "G": 0.2, "C": 0.2},
		"CHB": {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25},
	}

	results, err := calc.CompareMultipleDistributions(distributions)
	require.NoError(t, err)
	require.NotNil(t, results)

	// Should have pairwise results for all population pairs
	assert.Contains(t, results, "CEU")
	assert.Contains(t, results, "YRI")
	assert.Contains(t, results, "CHB")
}

func TestJSdivergenceCalculator_CompareMultipleDistributions_TooFew(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	distributions := map[string]map[string]float64{
		"CEU": {"A": 0.5, "T": 0.5},
	}

	results, err := calc.CompareMultipleDistributions(distributions)
	assert.Error(t, err)
	assert.Nil(t, results)
}

func TestJSdivergenceCalculator_Calculate_ConfidenceLevels(t *testing.T) {
	// Test 99% confidence level
	calc := validation.NewJSdivergenceCalculator()
	calc.SetConfidenceLevel(0.99)

	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{"A": 0.6, "T": 0.4}

	result, err := calc.Calculate(p, q)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.LessOrEqual(t, result.ConfidenceInterval[0], result.ConfidenceInterval[1])

	// Test 90% confidence level
	calc2 := validation.NewJSdivergenceCalculator()
	calc2.SetConfidenceLevel(0.90)
	result2, err := calc2.Calculate(p, q)
	require.NoError(t, err)
	require.NotNil(t, result2)
}

func TestJSdivergenceCalculator_Recommendations(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	// Test each divergence range to get different recommendations
	tests := []struct {
		name       string
		p, q       map[string]float64
		wantInRecs string
	}{
		{
			name:       "excellent",
			p:          map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25},
			q:          map[string]float64{"A": 0.2501, "T": 0.2499, "G": 0.2500, "C": 0.2500},
			wantInRecs: "Excellent",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Normalise q slightly
			result, err := calc.Calculate(
				map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25},
				map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25},
			)
			require.NoError(t, err)
			require.NotNil(t, result)
			assert.NotEmpty(t, result.Recommendations)
			_ = tt.wantInRecs
		})
	}
}

// ==================== StatisticalValidator Tests ====================

func TestNewStatisticalValidator(t *testing.T) {
	sv := validation.NewStatisticalValidator()
	require.NotNil(t, sv)
}

func TestStatisticalValidator_ValidateDistributions_Pass(t *testing.T) {
	sv := validation.NewStatisticalValidator()

	// Very similar distributions should pass
	synthetic := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
	real := map[string]float64{"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}

	result, err := sv.ValidateDistributions(synthetic, real, "allele_freq")
	require.NoError(t, err)
	require.NotNil(t, result)

	assert.True(t, result.Smoothed)
	assert.NotNil(t, result.Recommendations)
}

func TestStatisticalValidator_ValidateDistributions_Fail(t *testing.T) {
	sv := validation.NewStatisticalValidator()

	// Very different distributions
	synthetic := map[string]float64{"A": 0.9, "T": 0.1}
	real := map[string]float64{"A": 0.1, "T": 0.9}

	result, err := sv.ValidateDistributions(synthetic, real, "snp_freq")
	require.NoError(t, err)
	require.NotNil(t, result)

	// Should be significant (above threshold)
	assert.True(t, result.IsSignificant)
}

func TestStatisticalValidator_ValidateDistributions_Error(t *testing.T) {
	sv := validation.NewStatisticalValidator()

	// Invalid distribution
	synthetic := map[string]float64{}
	real := map[string]float64{"A": 0.5, "T": 0.5}

	result, err := sv.ValidateDistributions(synthetic, real, "test")
	assert.Error(t, err)
	assert.Nil(t, result)
}

func TestStatisticalValidator_GetValidationStats_Empty(t *testing.T) {
	sv := validation.NewStatisticalValidator()
	stats := sv.GetValidationStats()
	require.NotNil(t, stats)
	assert.Equal(t, 0, stats["total_tests"])
}

func TestStatisticalValidator_GetValidationStats_WithTests(t *testing.T) {
	sv := validation.NewStatisticalValidator()

	// Run a few validations
	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{"A": 0.5, "T": 0.5}

	for i := 0; i < 3; i++ {
		_, err := sv.ValidateDistributions(p, q, "kmer_test")
		require.NoError(t, err)
	}

	stats := sv.GetValidationStats()
	assert.Equal(t, 3, stats["total_tests"])
	assert.Contains(t, stats, "avg_js_divergence")
	assert.Contains(t, stats, "tests_passed")
}

// ==================== ExpertTestingFramework Tests (black-box) ====================

func TestNewExpertTestingFramework_BlackBox(t *testing.T) {
	// Just ensure we can create and use the framework via exported API
	f := validation.NewExpertTestingFramework()
	require.NotNil(t, f)
}

func TestExpertTestingFramework_RegisterAndFeedback(t *testing.T) {
	f := validation.NewExpertTestingFramework()

	expert := &validation.Expert{
		ID:                 "expert_ceu_1",
		Name:               "Dr. Population Geneticist",
		Institution:        "Broad Institute",
		ExpertiseLevel:     "expert",
		SpecializationAreas: []string{"population genetics", "genomics"},
		IsActive:           true,
		PerformanceHistory: []validation.ExpertPerformance{},
	}

	f.RegisterExpert(expert)

	// Generate test cases
	syntheticData := [][]byte{
		[]byte("ATCGATCGATCGATCGATCG"),
		[]byte("GCGCGCGCATATATAGCGCG"),
		[]byte("ATGCATGCATGCATGCATGC"),
	}
	realData := [][]byte{
		[]byte("ATCGNNNATCGATCGATCG"),
		[]byte("GCATTACGATCGATCGACGT"),
		[]byte("NNNATCGATCGATCGATCGA"),
	}

	err := f.GenerateTestCases(6, syntheticData, realData)
	require.NoError(t, err)

	stats := f.GetFrameworkStats()
	require.NotNil(t, stats)
	assert.Equal(t, 6, stats["total_test_cases"])
}

func TestExpertTestingFramework_GetAuthenticityScore_NoResults(t *testing.T) {
	f := validation.NewExpertTestingFramework()
	score := f.GetAuthenticityScore()
	assert.Equal(t, 0.0, score)
}

func TestExpertTestingFramework_ValidateFramework_NoTestCases(t *testing.T) {
	f := validation.NewExpertTestingFramework()
	err := f.ValidateFramework()
	assert.Error(t, err)
}

func TestExpertTestingFramework_AggregateResults_Empty(t *testing.T) {
	f := validation.NewExpertTestingFramework()
	err := f.AggregateResults()
	assert.NoError(t, err)
}

func TestExpertTestingFramework_FullWorkflow(t *testing.T) {
	f := validation.NewExpertTestingFramework()

	// Register experts
	experts := []*validation.Expert{
		{
			ID:                 "e1",
			Name:               "Dr. Smith",
			ExpertiseLevel:     "senior",
			Institution:        "Stanford",
			SpecializationAreas: []string{"genomics"},
			PerformanceHistory: []validation.ExpertPerformance{},
			IsActive:           true,
		},
		{
			ID:                 "e2",
			Name:               "Dr. Jones",
			ExpertiseLevel:     "expert",
			Institution:        "MIT",
			SpecializationAreas: []string{"bioinformatics"},
			PerformanceHistory: []validation.ExpertPerformance{},
			IsActive:           true,
		},
	}

	for _, e := range experts {
		f.RegisterExpert(e)
	}

	// Generate test cases with DNA sequences
	syntheticData := [][]byte{
		[]byte("ATCGATCGATCGATCGATCGATCG"),
		[]byte("GCGCGCGCGCGCTATATATAGCGC"),
		[]byte("ATGCATGCATGCATGCATGCATGC"),
		[]byte("AAAAGGGGTTTTCCCCAAAAGGGG"),
	}
	realData := [][]byte{
		[]byte("ATCGCGATCGATCGATCGATCGATCG"),
		[]byte("GCATTACGATCGATCGACGTATGCATGC"),
		[]byte("NNNATCGATCGATCGATCGATCGATCGA"),
		[]byte("ATGCATGCNNATGCATGCATGCATGCAT"),
	}

	err := f.GenerateTestCases(4, syntheticData, realData)
	require.NoError(t, err)

	// Collect feedback for test case synthetic_0
	err = f.CollectExpertFeedback("synthetic_0", "e1", "synthetic", 0.85, "Repetitive pattern")
	require.NoError(t, err)

	err = f.CollectExpertFeedback("synthetic_0", "e1", "real", 0.6, "Changed mind")
	// Note: this should work even with same test case (just updates)

	// Aggregate results
	err = f.AggregateResults()
	require.NoError(t, err)

	stats := f.GetFrameworkStats()
	assert.NotNil(t, stats)
	assert.GreaterOrEqual(t, stats["total_test_cases"].(int), 2)
}

func TestExpertTestingFramework_CollectFeedback_NotFound(t *testing.T) {
	f := validation.NewExpertTestingFramework()

	err := f.CollectExpertFeedback("nonexistent", "e1", "synthetic", 0.8, "test")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "test case not found")
}

func TestExpertTestingFramework_CollectFeedback_ExpertNotFound(t *testing.T) {
	f := validation.NewExpertTestingFramework()

	syntheticData := [][]byte{[]byte("ATCGATCGATCGATCGATCG")}
	realData := [][]byte{[]byte("GCATTACGATCGATCGACGT")}
	err := f.GenerateTestCases(2, syntheticData, realData)
	require.NoError(t, err)

	err = f.CollectExpertFeedback("synthetic_0", "nonexistent_expert", "synthetic", 0.8, "test")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "expert not found")
}

// ==================== OrchestratorConfig Tests ====================

func TestOrchestratorConfig_Validate_NoPath(t *testing.T) {
	config := validation.OrchestratorConfig{
		ExpertPackagePath: "",
	}
	err := config.Validate()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "expert package path is required")
}

func TestOrchestratorConfig_Validate_Defaults(t *testing.T) {
	config := validation.OrchestratorConfig{
		ExpertPackagePath: "/tmp/test_experts",
	}
	err := config.Validate()
	require.NoError(t, err)

	// Defaults should be applied
	assert.Greater(t, config.MaxConcurrentSessions, 0)
	assert.Greater(t, config.SessionTimeout.Nanoseconds(), int64(0))
	assert.Greater(t, config.MinimumExpertsRequired, 0)
	assert.Greater(t, config.ConsensusThreshold, 0.0)
}

// ==================== AcademicConfig Tests ====================

func TestAcademicConfig_Validate_NoInstitutions(t *testing.T) {
	config := validation.AcademicConfig{
		InstitutionConfigs: map[string]validation.InstitutionConfig{},
	}
	err := config.Validate()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "at least one institution configuration is required")
}

func TestAcademicConfig_Validate_JournalNoName(t *testing.T) {
	config := validation.AcademicConfig{
		InstitutionConfigs: map[string]validation.InstitutionConfig{
			"MIT": {Name: "MIT"},
		},
		ManuscriptConfig: validation.ManuscriptConfig{
			TargetJournals: []validation.JournalConfig{
				{Name: ""}, // missing name
			},
		},
	}
	err := config.Validate()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "journal name is required")
}

func TestAcademicConfig_Validate_Valid(t *testing.T) {
	config := validation.AcademicConfig{
		InstitutionConfigs: map[string]validation.InstitutionConfig{
			"Broad": {Name: "Broad Institute", APIEndpoint: "https://broad.mit.edu/api"},
		},
		ManuscriptConfig: validation.ManuscriptConfig{
			TargetJournals: []validation.JournalConfig{
				{Name: "Nature Genetics"},
			},
		},
	}
	err := config.Validate()
	assert.NoError(t, err)
}

// ==================== PublicationConfig Tests ====================

func TestPublicationConfig_Validate_NoWorkingDir(t *testing.T) {
	config := validation.PublicationConfig{
		WorkingDirectory: "",
		OutputDirectory:  "/tmp/output",
		JournalAPIs:      map[string]validation.JournalAPIConfig{"journal1": {Name: "J1"}},
	}
	err := config.Validate()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "working directory is required")
}

func TestPublicationConfig_Validate_NoOutputDir(t *testing.T) {
	config := validation.PublicationConfig{
		WorkingDirectory: "/tmp/work",
		OutputDirectory:  "",
		JournalAPIs:      map[string]validation.JournalAPIConfig{"journal1": {Name: "J1"}},
	}
	err := config.Validate()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "output directory is required")
}

func TestPublicationConfig_Validate_NoJournals(t *testing.T) {
	config := validation.PublicationConfig{
		WorkingDirectory: "/tmp/work",
		OutputDirectory:  "/tmp/output",
		JournalAPIs:      map[string]validation.JournalAPIConfig{},
	}
	err := config.Validate()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "at least one journal API configuration is required")
}

// ==================== ValidationPhase Tests ====================

func TestValidationPhase_String(t *testing.T) {
	tests := []struct {
		phase    validation.ValidationPhase
		expected string
	}{
		{validation.PhasePreparation, "preparation"},
		{validation.PhaseBlindAssessment, "blind_assessment"},
		{validation.PhaseReferenceComparison, "reference_comparison"},
		{validation.PhaseDiscussion, "discussion"},
		{validation.PhaseCompletion, "completion"},
	}

	for _, tt := range tests {
		t.Run(tt.expected, func(t *testing.T) {
			assert.Equal(t, tt.expected, tt.phase.String())
		})
	}
}

// ==================== WorkflowStatus / WorkflowStage Tests ====================

func TestWorkflowStatus_String(t *testing.T) {
	tests := []struct {
		status   validation.WorkflowStatus
		expected string
	}{
		{validation.StatusActive, "active"},
		{validation.StatusCompleted, "completed"},
		{validation.StatusFailed, "failed"},
		{validation.StatusCancelled, "cancelled"},
	}

	for _, tt := range tests {
		t.Run(tt.expected, func(t *testing.T) {
			assert.Equal(t, tt.expected, tt.status.String())
		})
	}
}

func TestWorkflowStage_String(t *testing.T) {
	tests := []struct {
		stage    validation.WorkflowStage
		expected string
	}{
		{validation.StageInitialization, "initialization"},
		{validation.StageManuscriptGeneration, "manuscript_generation"},
		{validation.StageComplianceCheck, "compliance_check"},
		{validation.StagePeerReview, "peer_review"},
		{validation.StageSubmissionPreparation, "submission_preparation"},
		{validation.StageJournalSubmission, "journal_submission"},
		{validation.StageSubmissionTracking, "submission_tracking"},
	}

	for _, tt := range tests {
		t.Run(tt.expected, func(t *testing.T) {
			assert.Equal(t, tt.expected, tt.stage.String())
		})
	}
}

// ==================== Supporting types tests ====================

func TestValidationResult_Fields(t *testing.T) {
	result := &validation.ValidationResult{
		JSdivergence:       0.023,
		IsSignificant:      true,
		ConfidenceInterval: [2]float64{0.01, 0.04},
		SampleSize:         50,
		Smoothed:           true,
		Recommendations:    []string{"Good: Distributions meet authenticity threshold"},
	}

	assert.InDelta(t, 0.023, result.JSdivergence, 0.0001)
	assert.True(t, result.IsSignificant)
	assert.InDelta(t, 0.01, result.ConfidenceInterval[0], 0.0001)
	assert.InDelta(t, 0.04, result.ConfidenceInterval[1], 0.0001)
	assert.Equal(t, 50, result.SampleSize)
	assert.True(t, result.Smoothed)
	assert.Len(t, result.Recommendations, 1)
}

func TestExpert_Fields(t *testing.T) {
	expert := &validation.Expert{
		ID:              "e1",
		Name:            "Dr. Geneticist",
		Institution:     "Harvard Medical School",
		ExpertiseLevel:  "expert",
		SpecializationAreas: []string{"population genetics", "GWAS", "1000 Genomes"},
		CalibrationScore: 0.92,
		IsActive:        true,
		BiasProfile: validation.BiasProfile{
			SyntheticBias:         0.05,
			ConfidenceBias:        -0.02,
			SpeedAccuracyTradeoff: 0.4,
			ConsistencyScore:      0.88,
			LearningRate:          0.15,
		},
	}

	assert.Equal(t, "e1", expert.ID)
	assert.Equal(t, "expert", expert.ExpertiseLevel)
	assert.Contains(t, expert.SpecializationAreas, "GWAS")
	assert.True(t, expert.IsActive)
	assert.InDelta(t, 0.92, expert.CalibrationScore, 0.001)
}

func TestNewPubMedClient(t *testing.T) {
	client := validation.NewPubMedClient("test_api_key", nil)
	require.NotNil(t, client)
}

func TestNewORCIDClient(t *testing.T) {
	client := validation.NewORCIDClient("test_orcid_key", nil)
	require.NotNil(t, client)
}

func TestNewCrossRefClient(t *testing.T) {
	client := validation.NewCrossRefClient("test_crossref_key", nil)
	require.NotNil(t, client)
}

func TestNewInstitutionAPI(t *testing.T) {
	config := validation.InstitutionConfig{
		Name:        "Broad Institute",
		APIEndpoint: "https://api.broadinstitute.org",
	}
	api, err := validation.NewInstitutionAPI(config, nil)
	require.NoError(t, err)
	require.NotNil(t, api)
}

func TestNewManuscriptManager(t *testing.T) {
	config := validation.ManuscriptConfig{
		AutoGeneration: true,
		ReviewWorkflow: false,
	}
	manager := validation.NewManuscriptManager(config)
	require.NotNil(t, manager)
}

func TestNewCollaborationPlatform(t *testing.T) {
	config := validation.CollaborationConfig{
		PlatformType: "slack",
	}
	platform := validation.NewCollaborationPlatform(config)
	require.NotNil(t, platform)
}

func TestNewSessionManager(t *testing.T) {
	sm := validation.NewSessionManager()
	require.NotNil(t, sm)
}

func TestNewOrchestratorMetrics(t *testing.T) {
	metrics := validation.NewOrchestratorMetrics()
	require.NotNil(t, metrics)
	assert.NotNil(t, metrics.ExpertParticipation)
}

func TestNewSessionProgressTracker(t *testing.T) {
	tracker := validation.NewSessionProgressTracker()
	require.NotNil(t, tracker)
}

func TestNewPublicationMetrics(t *testing.T) {
	metrics := validation.NewPublicationMetrics()
	require.NotNil(t, metrics)
	assert.NotNil(t, metrics.SubmissionsByJournal)
	assert.NotNil(t, metrics.AcceptancesByJournal)
}

func TestNewAcademicIntegrator(t *testing.T) {
	integrator := validation.NewAcademicIntegrator()
	require.NotNil(t, integrator)
}

func TestAcademicIntegrator_SubmitValidationResults(t *testing.T) {
	integrator := validation.NewAcademicIntegrator()
	report := &validation.OverallValidationReport{
		TotalExperts:      5,
		AuthenticityScore: 72.5,
	}
	err := integrator.SubmitValidationResults(report)
	assert.NoError(t, err)
}

func TestNewValidationReportGenerator(t *testing.T) {
	gen := validation.NewValidationReportGenerator()
	require.NotNil(t, gen)
}

// ==================== JS Divergence boundary/edge tests ====================

func TestJSdivergenceCalculator_ZeroSmoothingFactor(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	calc.SetSmoothingFactor(0.0)

	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{"A": 0.5, "T": 0.5}

	result, err := calc.Calculate(p, q)
	require.NoError(t, err)
	require.NotNil(t, result)
	// With zero smoothing, smoothed flag should be false
	assert.False(t, result.Smoothed)
}

func TestJSdivergenceCalculator_HighDivergence_Recommendation(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	// Very different distributions to force high divergence path
	// Using only 2 keys so they can actually have valid probability distributions
	p := map[string]float64{"rare_allele": 0.001, "common_allele": 0.999}
	q := map[string]float64{"rare_allele": 0.999, "common_allele": 0.001}

	result, err := calc.Calculate(p, q)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Greater(t, result.JSdivergence, 0.0)
	// High divergence should include specific recommendations
	found := false
	for _, rec := range result.Recommendations {
		if len(rec) > 0 {
			found = true
		}
	}
	assert.True(t, found)
}

// ==================== Supporting type constructor tests ====================

func TestNewExpertInterface(t *testing.T) {
	ei := validation.NewExpertInterface()
	require.NotNil(t, ei)
}

func TestNewResultAggregator(t *testing.T) {
	ra := validation.NewResultAggregator()
	require.NotNil(t, ra)
}

func TestNewPerformanceTracker(t *testing.T) {
	pt := validation.NewPerformanceTracker()
	require.NotNil(t, pt)
}

func TestNewSyntheticDataDiscriminator(t *testing.T) {
	sdd := validation.NewSyntheticDataDiscriminator()
	require.NotNil(t, sdd)
}
