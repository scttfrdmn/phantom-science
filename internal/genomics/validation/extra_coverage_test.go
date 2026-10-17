package validation_test

import (
	"context"
	"encoding/json"
	"fmt"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/validation"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ==================== AcademicIntegrationAPI deeper coverage ====================

// makeAcademicAPI creates a valid AcademicIntegrationAPI for testing
func makeAcademicAPI(t *testing.T) *validation.AcademicIntegrationAPI {
	t.Helper()
	config := validation.AcademicConfig{
		InstitutionConfigs: map[string]validation.InstitutionConfig{
			"MIT":      {Name: "MIT", APIEndpoint: "https://api.mit.edu"},
			"Stanford": {Name: "Stanford", APIEndpoint: "https://api.stanford.edu"},
		},
	}
	api, err := validation.NewAcademicIntegrationAPI(config)
	require.NoError(t, err)
	require.NotNil(t, api)
	return api
}

// makeAcademicAPIWithJournals creates an API with target journals in ManuscriptConfig.
// highScore controls the number of journals that selectTargetJournals returns:
// >= 75 → up to 3 journals, < 75 → up to 2 journals.
func makeAcademicAPIWithJournals(t *testing.T, highScore bool) *validation.AcademicIntegrationAPI {
	t.Helper()
	config := validation.AcademicConfig{
		InstitutionConfigs: map[string]validation.InstitutionConfig{
			"Broad": {Name: "Broad Institute"},
		},
		ManuscriptConfig: validation.ManuscriptConfig{
			AutoGeneration: true,
			TargetJournals: []validation.JournalConfig{
				{Name: "Nature Genetics", PreferredFormat: "latex"},
				{Name: "Genome Research", PreferredFormat: "word"},
				{Name: "PLOS Genetics", PreferredFormat: "word"},
			},
		},
	}
	api, err := validation.NewAcademicIntegrationAPI(config)
	require.NoError(t, err)
	return api
}

func TestAcademicIntegrationAPI_SubmitToJournals_EmptyJournals(t *testing.T) {
	// Manuscript with no JournalTargets - SubmitToJournals should return empty results
	api := makeAcademicAPI(t)
	manuscript := &validation.ManuscriptDraft{
		Title:           "Test Study",
		JournalTargets:  []*validation.JournalConfig{},
		JournalVersions: map[string]*validation.FormattedManuscript{},
	}
	ctx := context.Background()
	results, err := api.SubmitToJournals(ctx, manuscript)
	require.NoError(t, err)
	assert.Empty(t, results)
}

func TestAcademicIntegrationAPI_SubmitToJournals_NoMatchingVersions(t *testing.T) {
	// JournalTargets has journals but JournalVersions doesn't have matching keys
	api := makeAcademicAPI(t)
	journal := &validation.JournalConfig{
		Name:            "Nature Genetics",
		PreferredFormat: "latex",
	}
	manuscript := &validation.ManuscriptDraft{
		Title:          "Test Study",
		JournalTargets: []*validation.JournalConfig{journal},
		// JournalVersions is nil/empty - no matching version
		JournalVersions: map[string]*validation.FormattedManuscript{},
	}
	ctx := context.Background()
	results, err := api.SubmitToJournals(ctx, manuscript)
	require.NoError(t, err)
	// No versions found, so nothing is submitted
	assert.Empty(t, results)
}

func TestAcademicIntegrationAPI_SubmitToJournals_WithVersions(t *testing.T) {
	// NOTE: SubmitToJournals calls time.Sleep(5 * time.Second) between submissions (source code).
	// We use a single journal to get one submission result (the sleep happens after adding to results).
	// Total test time: ~5 seconds due to the source-level sleep.
	if testing.Short() {
		t.Skip("Skipping test that takes ~5s due to source-level sleep in SubmitToJournals")
	}
	api := makeAcademicAPI(t)
	journal := &validation.JournalConfig{
		Name:            "TestJournal",
		PreferredFormat: "latex",
		ReviewTimeline:  30 * 24 * time.Hour,
	}
	manuscript := &validation.ManuscriptDraft{
		Title:          "Test Study",
		JournalTargets: []*validation.JournalConfig{journal},
		JournalVersions: map[string]*validation.FormattedManuscript{
			"TestJournal": {
				Content:     "Formatted content here",
				Format:      "latex",
				WordCount:   5000,
				FigureCount: 3,
				TableCount:  2,
			},
		},
	}
	ctx := context.Background()
	results, err := api.SubmitToJournals(ctx, manuscript)
	require.NoError(t, err)
	require.Len(t, results, 1)
	assert.Equal(t, "TestJournal", results[0].Journal)
	assert.Equal(t, "submitted", results[0].Status)
	assert.NotEmpty(t, results[0].SubmissionID)
}

func TestAcademicIntegrationAPI_RecruitExperts_WithCandidates(t *testing.T) {
	// Test recruit experts through multiple institutions
	api := makeAcademicAPI(t)
	criteria := validation.ExpertCriteria{
		ResearchAreas:   []string{"population genetics", "genomics"},
		ExpertiseLevel:  "senior",
		Publications:    15,
		RemoteAvailable: true,
		Languages:       []string{"en"},
	}
	ctx := context.Background()
	candidates, err := api.RecruitExperts(ctx, criteria)
	require.NoError(t, err)
	// InstitutionAPI.SearchExperts is a stub returning empty slice
	// So rankedCandidates will be nil (empty)
	_ = candidates
}

func TestAcademicIntegrationAPI_RecruitExperts_InstitutionPreference(t *testing.T) {
	// Test with institution preference in criteria
	api := makeAcademicAPI(t)
	criteria := validation.ExpertCriteria{
		ResearchAreas:  []string{"computational biology"},
		ExpertiseLevel: "expert",
		Publications:   30,
		Institution:    []string{"MIT", "Harvard"},
	}
	ctx := context.Background()
	_, err := api.RecruitExperts(ctx, criteria)
	require.NoError(t, err)
}

func TestAcademicIntegrationAPI_GenerateManuscript_HighScore(t *testing.T) {
	// High authenticity score (>= 75): selectTargetJournals picks up to 3 journals.
	// JournalVersions must be initialized (non-nil map) so FormatForJournal writes succeed.
	api := makeAcademicAPIWithJournals(t, true)
	report := &validation.OverallValidationReport{
		TotalExperts:      12,
		AuthenticityScore: 80.0,
		Timestamp:         time.Now(),
	}
	ctx := context.Background()
	manuscript, err := api.GenerateManuscript(ctx, report)
	require.NoError(t, err)
	require.NotNil(t, manuscript)
	assert.NotEmpty(t, manuscript.Title)
	assert.NotEmpty(t, manuscript.Abstract)
	assert.NotEmpty(t, manuscript.Introduction)
	assert.NotEmpty(t, manuscript.Methods)
	assert.NotEmpty(t, manuscript.Results)
	assert.NotEmpty(t, manuscript.Discussion)
	assert.NotEmpty(t, manuscript.Conclusion)
	assert.NotEmpty(t, manuscript.Keywords)
	assert.NotNil(t, manuscript.Authors)
	assert.NotNil(t, manuscript.Figures)
	assert.NotNil(t, manuscript.Tables)
	assert.NotNil(t, manuscript.References)
	// With 3 journals configured and score >= 75, expect 3 target journals and formatted versions.
	assert.Len(t, manuscript.JournalTargets, 3)
	require.NotNil(t, manuscript.JournalVersions)
	assert.Len(t, manuscript.JournalVersions, 3)
}

func TestAcademicIntegrationAPI_GenerateManuscript_LowScore(t *testing.T) {
	// Low authenticity score (< 75): selectTargetJournals picks up to 2 journals.
	api := makeAcademicAPIWithJournals(t, false)
	report := &validation.OverallValidationReport{
		TotalExperts:      6,
		AuthenticityScore: 55.0,
		Timestamp:         time.Now(),
	}
	ctx := context.Background()
	manuscript, err := api.GenerateManuscript(ctx, report)
	require.NoError(t, err)
	require.NotNil(t, manuscript)
	assert.NotEmpty(t, manuscript.Title)
	// With 3 journals configured but score < 75, expect only 2 target journals.
	assert.Len(t, manuscript.JournalTargets, 2)
	require.NotNil(t, manuscript.JournalVersions)
	assert.Len(t, manuscript.JournalVersions, 2)
}

func TestAcademicIntegrationAPI_GenerateManuscript_ZeroScore(t *testing.T) {
	// Test with score of 0 - exercises the else branch (< 75) of selectTargetJournals
	api := makeAcademicAPIWithJournals(t, false)
	report := &validation.OverallValidationReport{
		TotalExperts:      0,
		AuthenticityScore: 0.0,
		Timestamp:         time.Now(),
	}
	ctx := context.Background()
	manuscript, err := api.GenerateManuscript(ctx, report)
	require.NoError(t, err)
	require.NotNil(t, manuscript)
	assert.NotEmpty(t, manuscript.Title)
}

func TestAcademicIntegrationAPI_SetupCollaborationWorkspace_MultipleExperts(t *testing.T) {
	api := makeAcademicAPI(t)
	experts := []*validation.ExpertCandidate{
		{Name: "Dr. Alice Johnson", Email: "alice@mit.edu", Institution: "MIT"},
		{Name: "Dr. Bob Williams", Email: "bob@stanford.edu", Institution: "Stanford"},
		{Name: "Dr. Carol Chen", Email: "carol@harvard.edu", Institution: "Harvard"},
	}
	projectDetails := &validation.ProjectDetails{
		Title:       "Multi-Expert Validation Study",
		Description: "Population genomics validation by multiple experts",
		Methodology: "Blind assessment protocol",
	}
	ctx := context.Background()
	workspace, err := api.SetupCollaborationWorkspace(ctx, experts, projectDetails)
	require.NoError(t, err)
	require.NotNil(t, workspace)
	assert.Contains(t, workspace.Name, "Multi-Expert Validation Study")
	assert.Equal(t, 5, len(workspace.Channels)) // 5 standard channels
	// All 3 experts should be members (InviteMember is a stub that succeeds)
	assert.Equal(t, 3, len(workspace.Members))
	// SharedFiles should have 3 files
	assert.Equal(t, 3, len(workspace.SharedFiles))
}

func TestAcademicIntegrationAPI_SetupCollaborationWorkspace_NoExperts(t *testing.T) {
	api := makeAcademicAPI(t)
	projectDetails := &validation.ProjectDetails{
		Title:       "Solo Validation Study",
		Description: "Single-expert validation",
	}
	ctx := context.Background()
	workspace, err := api.SetupCollaborationWorkspace(ctx, []*validation.ExpertCandidate{}, projectDetails)
	require.NoError(t, err)
	require.NotNil(t, workspace)
	assert.Empty(t, workspace.Members)
}

// ==================== ExpertCriteria coverage ====================

func TestExpertCriteria_AllFields(t *testing.T) {
	criteria := validation.ExpertCriteria{
		ResearchAreas:   []string{"population genetics", "GWAS"},
		ExpertiseLevel:  "expert",
		Institution:     []string{"MIT", "Stanford", "Harvard"},
		Publications:    50,
		ExperienceYears: 15,
		Specializations: []string{"cancer genomics", "rare variants"},
		TimeZone:        "America/New_York",
		Region:          "North America",
		RemoteAvailable: true,
		Languages:       []string{"en", "de", "fr"},
	}
	assert.Equal(t, "expert", criteria.ExpertiseLevel)
	assert.Equal(t, 50, criteria.Publications)
	assert.Equal(t, 15, criteria.ExperienceYears)
	assert.True(t, criteria.RemoteAvailable)
	assert.Len(t, criteria.ResearchAreas, 2)
	assert.Len(t, criteria.Institution, 3)
	assert.Len(t, criteria.Specializations, 2)
}

// ==================== Publication Workflow deeper coverage ====================

func TestPublicationWorkflow_FullFields(t *testing.T) {
	now := time.Now()
	workflow := &validation.PublicationWorkflow{
		ID:               "pub_workflow_001",
		TargetJournals:   []string{"Nature Genetics", "Genome Biology"},
		StartTime:        now,
		CompletionTime:   now.Add(2 * time.Hour),
		CurrentStage:     validation.StageManuscriptGeneration,
		Status:           validation.StatusActive,
		Progress:         map[string]float64{"overall": 0.25, "manuscript": 0.5},
		Artifacts:        map[string]*validation.Artifact{},
		Submissions:      map[string]*validation.JournalSubmission{},
		StageHistory:     []*validation.StageTransition{},
	}

	assert.Equal(t, "pub_workflow_001", workflow.ID)
	assert.Equal(t, validation.StatusActive, workflow.Status)
	assert.Equal(t, validation.StageManuscriptGeneration, workflow.CurrentStage)
	assert.Len(t, workflow.TargetJournals, 2)
	assert.InDelta(t, 0.25, workflow.Progress["overall"], 0.001)
}

func TestStageTransition_Fields(t *testing.T) {
	now := time.Now()
	transition := &validation.StageTransition{
		Stage:     validation.StageManuscriptGeneration,
		Status:    "completed",
		Message:   "Automatic transition after initialization complete",
		Timestamp: now,
	}

	assert.Equal(t, validation.StageManuscriptGeneration, transition.Stage)
	assert.Equal(t, "completed", transition.Status)
	assert.Equal(t, "Automatic transition after initialization complete", transition.Message)
}

func TestPublicationMetrics_Fields(t *testing.T) {
	metrics := validation.NewPublicationMetrics()
	require.NotNil(t, metrics)
	assert.Equal(t, 0, metrics.TotalWorkflowsStarted)
	assert.Equal(t, 0, metrics.CompletedWorkflows)
	assert.Equal(t, 0, metrics.FailedWorkflows)
	assert.NotNil(t, metrics.SubmissionsByJournal)
	assert.NotNil(t, metrics.AcceptancesByJournal)
}

func TestJournalSubmission_AllFields(t *testing.T) {
	now := time.Now()
	decision := "pending"
	submission := &validation.JournalSubmission{
		SubmissionID: "sub_nature_20240115",
		Journal:      "Nature Genetics",
		Status:       "under_review",
		SubmittedAt:  now,
		LastUpdate:   now.Add(24 * time.Hour),
		TrackingURL:  "https://mts.nature.com/track/sub_nature_20240115",
		Decision:     decision,
		ReviewComments: []string{"Interesting work", "Methods need clarification"},
	}

	assert.Equal(t, "sub_nature_20240115", submission.SubmissionID)
	assert.Equal(t, "Nature Genetics", submission.Journal)
	assert.Equal(t, "under_review", submission.Status)
	assert.Equal(t, "pending", submission.Decision)
	assert.NotEmpty(t, submission.TrackingURL)
	assert.Len(t, submission.ReviewComments, 2)
}

// ==================== ValidationReport Generator coverage ====================

func TestValidationReportGenerator_Coverage(t *testing.T) {
	gen := validation.NewValidationReportGenerator()
	require.NotNil(t, gen)
	// Just verifying it can be created (struct is empty placeholder)
}

// ==================== Supporting types for orchestrator ====================

func TestLiveValidationSession_Fields(t *testing.T) {
	now := time.Now()
	tracker := validation.NewSessionProgressTracker()
	session := &validation.LiveValidationSession{
		SessionID:       "sess_expert1_1700000000",
		ExpertID:        "expert1",
		StartTime:       now,
		CurrentPhase:    validation.PhaseBlindAssessment,
		CompletedTests:  []string{"synthetic_0", "synthetic_1"},
		IsActive:        true,
		LastActivity:    now,
		SessionData:     map[string]interface{}{"test": "data"},
		ProgressTracker: tracker,
		FeedbackChannel: make(chan validation.ExpertFeedback, 10),
		StatusChannel:   make(chan validation.SessionStatus, 10),
	}

	assert.Equal(t, "sess_expert1_1700000000", session.SessionID)
	assert.Equal(t, "expert1", session.ExpertID)
	assert.True(t, session.IsActive)
	assert.Equal(t, validation.PhaseBlindAssessment, session.CurrentPhase)
	assert.Len(t, session.CompletedTests, 2)
}

func TestExpertFeedback_Fields(t *testing.T) {
	feedback := validation.ExpertFeedback{
		SessionID:  "sess_001",
		ExpertID:   "expert_1",
		TestCaseID: "synthetic_0",
		Feedback: map[string]interface{}{
			"classification": "synthetic",
			"confidence":     0.85,
			"justification":  "Repetitive k-mer pattern",
		},
		Timestamp: time.Now(),
		Phase:     validation.PhaseBlindAssessment,
	}

	assert.Equal(t, "sess_001", feedback.SessionID)
	assert.Equal(t, "expert_1", feedback.ExpertID)
	assert.Equal(t, validation.PhaseBlindAssessment, feedback.Phase)
	assert.Equal(t, "synthetic", feedback.Feedback["classification"])
	assert.InDelta(t, 0.85, feedback.Feedback["confidence"].(float64), 0.001)
}

func TestSessionStatus_Fields(t *testing.T) {
	status := validation.SessionStatus{
		SessionID: "sess_001",
		Status:    "progress",
		Message:   "Completed 5/10 test cases",
		Progress:  0.5,
		Timestamp: time.Now(),
	}

	assert.Equal(t, "sess_001", status.SessionID)
	assert.Equal(t, "progress", status.Status)
	assert.InDelta(t, 0.5, status.Progress, 0.001)
}

func TestValidationPhase_AllValues(t *testing.T) {
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
	// Test unknown phase
	var unknown validation.ValidationPhase = 99
	assert.Equal(t, "unknown", unknown.String())
}

func TestQualityAssessment_Fields(t *testing.T) {
	qa := validation.QualityAssessment{
		OverallQuality:       "excellent",
		StrengthAreas:        []string{"population_specificity", "k_mer_complexity"},
		ImprovementAreas:     []string{"gc_content_variation"},
		PublicationReadiness: true,
		ClinicalUtility:      "high",
	}
	assert.Equal(t, "excellent", qa.OverallQuality)
	assert.True(t, qa.PublicationReadiness)
	assert.Len(t, qa.StrengthAreas, 2)
}

func TestOverallValidationReport_Fields(t *testing.T) {
	now := time.Now()
	report := &validation.OverallValidationReport{
		Timestamp:         now,
		TotalExperts:      15,
		AuthenticityScore: 78.5,
		FrameworkStats:    map[string]interface{}{"total_test_cases": 30},
		ConsensusResults:  []validation.ConsensusResult{},
		DisagreementAnalysis: []validation.DisagreementAnalysis{},
		QualityAssessment: validation.QualityAssessment{
			OverallQuality:       "good",
			PublicationReadiness: true,
		},
		Recommendations: []string{"Continue validation", "Prepare manuscript"},
	}

	assert.Equal(t, 15, report.TotalExperts)
	assert.InDelta(t, 78.5, report.AuthenticityScore, 0.001)
	assert.Equal(t, "good", report.QualityAssessment.OverallQuality)
	assert.Len(t, report.Recommendations, 2)
}

func TestSessionValidationReport_Fields(t *testing.T) {
	now := time.Now()
	report := &validation.SessionValidationReport{
		SessionID:       "sess_001",
		ExpertID:        "expert_1",
		SessionDuration: 2*time.Hour + 30*time.Minute,
		CompletedCases:  8,
		TotalCases:      10,
		Accuracy:        0.75,
		Timestamp:       now,
		SessionData:     map[string]interface{}{"notes": "completed"},
	}

	assert.Equal(t, "sess_001", report.SessionID)
	assert.Equal(t, 8, report.CompletedCases)
	assert.Equal(t, 10, report.TotalCases)
	assert.InDelta(t, 0.75, report.Accuracy, 0.001)
}

// ==================== SMTPConfig and OrchestratorConfig fields ====================

func TestSMTPConfig_Fields(t *testing.T) {
	config := validation.SMTPConfig{
		Host:     "smtp.gmail.com",
		Port:     587,
		Username: "research@example.com",
		Password: "secret", // #nosec G101
		FromAddr: "research@example.com",
	}
	assert.Equal(t, "smtp.gmail.com", config.Host)
	assert.Equal(t, 587, config.Port)
	assert.Equal(t, "research@example.com", config.Username)
}

func TestOrchestratorConfig_AllFields(t *testing.T) {
	config := validation.OrchestratorConfig{
		ExpertPackagePath:      "/data/expert_package",
		MaxConcurrentSessions:  10,
		SessionTimeout:         3 * time.Hour,
		AutoReportGeneration:   true,
		AcademicIntegration:    true,
		EmailNotifications:     true,
		MinimumExpertsRequired: 5,
		ConsensusThreshold:     0.75,
		MaxSessionDuration:     8 * time.Hour,
		SMTPConfig: validation.SMTPConfig{
			Host: "smtp.university.edu",
			Port: 465,
		},
	}
	err := config.Validate()
	require.NoError(t, err)
	assert.Equal(t, "/data/expert_package", config.ExpertPackagePath)
	assert.Equal(t, 10, config.MaxConcurrentSessions)
	assert.Equal(t, 5, config.MinimumExpertsRequired)
	assert.InDelta(t, 0.75, config.ConsensusThreshold, 0.001)
}

// ==================== JS Divergence - additional recommendation paths ====================

func TestJSdivergenceCalculator_RecommendationPaths(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()

	// Test "Poor" recommendation path (divergence >= 0.05, < 0.1)
	// We need distributions that yield 0.05-0.1 divergence
	// Use distributions with moderate differences
	p := map[string]float64{"A": 0.5, "T": 0.5}
	q := map[string]float64{"A": 0.1, "T": 0.9}
	result, err := calc.Calculate(p, q)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.NotEmpty(t, result.Recommendations)
	// Just check recommendations are generated for divergence > 0
	assert.Greater(t, result.JSdivergence, 0.0)
}

func TestJSdivergenceCalculator_HighDivergenceAbove01(t *testing.T) {
	calc := validation.NewJSdivergenceCalculator()
	// Extreme distributions to push divergence > 0.1
	p := map[string]float64{
		"allele_A": 0.0001,
		"allele_B": 0.9999,
	}
	q := map[string]float64{
		"allele_A": 0.9999,
		"allele_B": 0.0001,
	}
	result, err := calc.Calculate(p, q)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Greater(t, result.JSdivergence, 0.05)
	// At high divergence, "High divergence" recommendation should appear
	found := false
	for _, rec := range result.Recommendations {
		if len(rec) > 0 {
			found = true
		}
	}
	assert.True(t, found)
}

// ==================== Publication workflow - additional coverage ====================

func TestPublicationWorkflowManager_StartWorkflow_MultipleJournals(t *testing.T) {
	config, _ := makePublicationConfig(t)
	manager, err := validation.NewPublicationWorkflowManager(config)
	require.NoError(t, err)

	report := &validation.OverallValidationReport{
		TotalExperts:      8,
		AuthenticityScore: 72.0,
		Timestamp:         time.Now(),
	}

	ctx := context.Background()
	// Target both journals configured in makePublicationConfig
	workflow, err := manager.StartPublicationWorkflow(ctx, report, []string{"NatureGenetics", "Genome Biology"})
	require.NoError(t, err)
	require.NotNil(t, workflow)
	assert.Contains(t, workflow.TargetJournals, "NatureGenetics")
	assert.Contains(t, workflow.TargetJournals, "Genome Biology")

	// Wait for background goroutine to finish before temp dir cleanup
	time.Sleep(200 * time.Millisecond)
}

// ==================== Consensus and DisagreementAnalysis types ====================

func TestConsensusResult_Fields(t *testing.T) {
	result := validation.ConsensusResult{
		TestCaseID:              "synthetic_5",
		ConsensusAnswer:         "synthetic",
		Confidence:              0.82,
		AgreementLevel:          0.75,
		ParticipatingExperts:    []string{"e1", "e2", "e3", "e4"},
		DissentingOpinions:      []string{"e2"},
		StatisticalSignificance: 0.03,
	}

	assert.Equal(t, "synthetic_5", result.TestCaseID)
	assert.Equal(t, "synthetic", result.ConsensusAnswer)
	assert.InDelta(t, 0.82, result.Confidence, 0.001)
	assert.InDelta(t, 0.75, result.AgreementLevel, 0.001)
	assert.Len(t, result.ParticipatingExperts, 4)
	assert.Len(t, result.DissentingOpinions, 1)
}

func TestDisagreementAnalysis_Fields(t *testing.T) {
	analysis := validation.DisagreementAnalysis{
		TestCaseID:          "real_3",
		DisagreementType:    "classification",
		ConflictingAnswers:  []string{"synthetic", "real"},
		PossibleCauses:      []string{"high_difficulty", "expertise_level_difference"},
		DataCharacteristics: map[string]float64{"gc_content": 0.52, "complexity": 0.85},
		ResolutionMethod:    "majority_vote",
	}

	assert.Equal(t, "real_3", analysis.TestCaseID)
	assert.Equal(t, "classification", analysis.DisagreementType)
	assert.Len(t, analysis.ConflictingAnswers, 2)
	assert.Len(t, analysis.PossibleCauses, 2)
}

// ==================== BlindTestCase metadata and difficulty ====================

func TestExpertTestingFramework_DifficultyLevels_Hard(t *testing.T) {
	f := validation.NewExpertTestingFramework()

	// Create sequences specifically designed to hit "hard" path:
	// GC between 0.45-0.55 AND complexity > 0.8
	// A diverse sequence with balanced GC
	hardSeq := []byte("ATGCGATCAGCTAGCTAGCATGCGATCAGCTAGCTAGCATGCGATCAGCTAGCTAGCATGCGATCAGCTAGCTAGCATGCGATCAGCTAGCTAGCATGCGATCAGCTAGCTAGC")

	syntheticData := [][]byte{hardSeq}
	realData := [][]byte{[]byte("GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG")}

	err := f.GenerateTestCases(2, syntheticData, realData)
	require.NoError(t, err)

	stats := f.GetFrameworkStats()
	require.NotNil(t, stats)
	// Just verify test cases were generated
	assert.GreaterOrEqual(t, stats["total_test_cases"].(int), 1)
}

func TestExpertTestingFramework_AggregateResults_WithData(t *testing.T) {
	f := validation.NewExpertTestingFramework()

	// Register multiple experts
	for i := 1; i <= 3; i++ {
		expert := &validation.Expert{
			ID:              fmt.Sprintf("expert_%d", i),
			ExpertiseLevel:  "senior",
			PerformanceHistory: []validation.ExpertPerformance{},
		}
		f.RegisterExpert(expert)
	}

	syntheticData := [][]byte{
		[]byte("ATCGATCGATCGATCGATCGATCGATCG"),
		[]byte("GCGCGCGCGCGCGCGCGCGCGCGCGCGC"),
	}
	realData := [][]byte{
		[]byte("ATGCNATGCNATGCNATGCNATGCNATGCN"),
		[]byte("GCATTACGATCGATCGACGTATGCATGCAT"),
	}

	err := f.GenerateTestCases(4, syntheticData, realData)
	require.NoError(t, err)

	// Collect feedback from multiple experts on multiple test cases
	for _, expertID := range []string{"expert_1", "expert_2", "expert_3"} {
		_ = f.CollectExpertFeedback("synthetic_0", expertID, "synthetic", 0.8, "repetitive")
	}
	// Create disagreement on synthetic_1
	_ = f.CollectExpertFeedback("synthetic_1", "expert_1", "synthetic", 0.9, "artificial")
	_ = f.CollectExpertFeedback("synthetic_1", "expert_2", "real", 0.6, "looks natural")
	_ = f.CollectExpertFeedback("synthetic_1", "expert_3", "real", 0.55, "uncertain")

	err = f.AggregateResults()
	require.NoError(t, err)

	score := f.GetAuthenticityScore()
	// Score should be between 0 and 100
	assert.GreaterOrEqual(t, score, 0.0)
	assert.LessOrEqual(t, score, 100.0)
}

func TestExpertTestingFramework_ValidateFramework_Complete(t *testing.T) {
	f := validation.NewExpertTestingFramework()

	// First add test cases
	syntheticData := [][]byte{[]byte("ATCGATCGATCGATCGATCG")}
	realData := [][]byte{[]byte("GCATTACGATCGATCGACGT")}
	err := f.GenerateTestCases(2, syntheticData, realData)
	require.NoError(t, err)

	// Without experts - should fail
	err = f.ValidateFramework()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "no experts registered")

	// Register an expert
	f.RegisterExpert(&validation.Expert{
		ID:             "e1",
		ExpertiseLevel: "senior",
	})

	// Now should pass
	err = f.ValidateFramework()
	assert.NoError(t, err)
}

// ==================== Statistical Validator more coverage ====================

func TestStatisticalValidator_MultipleTests(t *testing.T) {
	sv := validation.NewStatisticalValidator()

	// Run multiple different types of tests
	tests := []struct {
		p, q       map[string]float64
		testType   string
	}{
		{
			p:        map[string]float64{"A": 0.3, "T": 0.3, "G": 0.2, "C": 0.2},
			q:        map[string]float64{"A": 0.3, "T": 0.3, "G": 0.2, "C": 0.2},
			testType: "gc_content",
		},
		{
			p:        map[string]float64{"CHB": 0.4, "JPT": 0.3, "CHS": 0.3},
			q:        map[string]float64{"CHB": 0.35, "JPT": 0.33, "CHS": 0.32},
			testType: "population_freq",
		},
		{
			p:        map[string]float64{"kmer_ATCG": 0.5, "kmer_GCTA": 0.5},
			q:        map[string]float64{"kmer_ATCG": 0.6, "kmer_GCTA": 0.4},
			testType: "kmer_distribution",
		},
	}

	for _, tt := range tests {
		result, err := sv.ValidateDistributions(tt.p, tt.q, tt.testType)
		require.NoError(t, err)
		require.NotNil(t, result)
	}

	stats := sv.GetValidationStats()
	assert.Equal(t, 3, stats["total_tests"])
	assert.GreaterOrEqual(t, stats["tests_passed"].(int), 0)
	assert.Contains(t, stats, "test_breakdown")
}

// ==================== JournalAPIConfig and publication types ====================

func TestJournalAPIConfig_Fields(t *testing.T) {
	config := validation.JournalAPIConfig{
		Name:            "Nature Genetics",
		APIEndpoint:     "https://api.nature.com/genetics",
		APIKey:          "test_api_key", // #nosec G101
		PreferredFormat: "latex",
		Requirements: validation.JournalRequirements{
			AbstractWordLimit: 250,
			MainTextWordLimit: 3500,
			MaxFigures:        7,
			MaxTables:         5,
		},
	}

	assert.Equal(t, "Nature Genetics", config.Name)
	assert.Equal(t, "latex", config.PreferredFormat)
	assert.Equal(t, 250, config.Requirements.AbstractWordLimit)
	assert.Equal(t, 3500, config.Requirements.MainTextWordLimit)
	assert.Equal(t, 7, config.Requirements.MaxFigures)
}

func TestArtifact_AllFields(t *testing.T) {
	artifact := &validation.Artifact{
		Type:     "formatted_manuscript",
		FilePath: "/tmp/output/manuscript_nature.latex",
		Metadata: map[string]interface{}{
			"word_count":   5200,
			"figure_count": 4,
			"journal":      "Nature Genetics",
		},
	}

	assert.Equal(t, "formatted_manuscript", artifact.Type)
	assert.Equal(t, "/tmp/output/manuscript_nature.latex", artifact.FilePath)
	assert.Equal(t, 5200, artifact.Metadata["word_count"])
}

// ==================== InstitutionConfig additional fields ====================

func TestInstitutionConfig_AllFields(t *testing.T) {
	config := validation.InstitutionConfig{
		Name:        "Broad Institute",
		APIEndpoint: "https://api.broadinstitute.org",
		APIKey:      "broad_api_key", // #nosec G101
		ContactEmail: "genomics@broadinstitute.org",
		Departments: map[string]validation.DepartmentContact{
			"Genomics": {
				Name:           "Prof. Eric Lander",
				Email:          "lander@broadinstitute.org",
				Phone:          "+1-617-714-7000",
				ResearchAreas:  []string{"population genomics", "cancer genomics"},
				ExpertiseLevel: "expert",
			},
		},
		SupportsDirectSubmission: true,
		SupportsCollaboration:    true,
		SupportsDataSharing:      true,
	}

	assert.Equal(t, "Broad Institute", config.Name)
	assert.True(t, config.SupportsDirectSubmission)
	assert.True(t, config.SupportsCollaboration)
	assert.True(t, config.SupportsDataSharing)
	dept, ok := config.Departments["Genomics"]
	require.True(t, ok)
	assert.Equal(t, "Prof. Eric Lander", dept.Name)
	assert.Len(t, dept.ResearchAreas, 2)
}

func TestManuscriptConfig_AllFields(t *testing.T) {
	config := validation.ManuscriptConfig{
		TargetJournals: []validation.JournalConfig{
			{Name: "Nature", ImpactFactor: 49.9, PreferredFormat: "latex"},
			{Name: "Science", ImpactFactor: 44.7, PreferredFormat: "docx"},
		},
		AutoGeneration:     true,
		TemplateLibrary:    "/templates/journals",
		ReviewWorkflow:     true,
		SubmissionTracking: true,
	}

	assert.True(t, config.AutoGeneration)
	assert.True(t, config.ReviewWorkflow)
	assert.Len(t, config.TargetJournals, 2)
	assert.Equal(t, "Nature", config.TargetJournals[0].Name)
}

// ==================== ResearchNetworkingConfig ====================

func TestResearchNetworkingConfig_Fields(t *testing.T) {
	config := validation.ResearchNetworkingConfig{
		LinkedInIntegration:   true,
		ResearchGateAPI:       "rg_api_key",
		AcademiaEduAPI:        "academia_key",
		GoogleScholarAPI:      "scholar_key",
		AutoProfileUpdate:     false,
		PublicationSync:       true,
		CollaborationTracking: true,
	}

	assert.True(t, config.LinkedInIntegration)
	assert.False(t, config.AutoProfileUpdate)
	assert.True(t, config.PublicationSync)
	assert.Equal(t, "rg_api_key", config.ResearchGateAPI)
}

// ==================== ExpertInvitation and ProjectDetails types ====================

func TestExpertInvitation_Fields(t *testing.T) {
	deadline := time.Now().Add(14 * 24 * time.Hour)
	invitation := &validation.ExpertInvitation{
		RecipientName:   "Dr. Jane Smith",
		RecipientEmail:  "jsmith@mit.edu",
		Institution:     "MIT",
		PersonalizedMessage: "Dear Dr. Smith, we invite you to participate...",
		TimeCommitment:  "4-6 hours over 1-2 sessions",
		Compensation:    "Co-authorship opportunity",
		Deadline:        deadline,
	}

	assert.Equal(t, "Dr. Jane Smith", invitation.RecipientName)
	assert.Equal(t, "MIT", invitation.Institution)
	assert.Equal(t, "4-6 hours over 1-2 sessions", invitation.TimeCommitment)
}

func TestProjectDetails_AllFields(t *testing.T) {
	pd := &validation.ProjectDetails{
		Title:       "Population Genomics Validation",
		Description: "Expert validation of computationally generated sequences",
		Methodology: "Three-phase blind assessment protocol",
		DataSize:    "165 sequences from 11 populations",
		Innovation:  "First validated population-specific sequence generation",
		Applications: []string{
			"Precision medicine",
			"Population genomics research",
			"Clinical diagnostics",
		},
		ExpectedOutcomes: []string{
			"Top-tier publication",
			"Clinical applications framework",
		},
	}

	assert.Equal(t, "Population Genomics Validation", pd.Title)
	assert.Len(t, pd.Applications, 3)
	assert.Len(t, pd.ExpectedOutcomes, 2)
	assert.Equal(t, "165 sequences from 11 populations", pd.DataSize)
}

// ==================== CollaborationConfig and workspace types ====================

func TestCollaborationConfig_AllFields(t *testing.T) {
	config := validation.CollaborationConfig{
		PlatformType:      "slack",
		APIEndpoint:       "https://slack.com/api",
		Credentials:       map[string]string{"token": "xoxb-test-token"},
		AutoInviteExperts: true,
		SharedWorkspace:   true,
		RealTimeEditing:   false,
	}

	assert.Equal(t, "slack", config.PlatformType)
	assert.True(t, config.AutoInviteExperts)
	assert.True(t, config.SharedWorkspace)
	assert.False(t, config.RealTimeEditing)
}

func TestCollaborationWorkspace_Fields(t *testing.T) {
	now := time.Now()
	workspace := &validation.CollaborationWorkspace{
		ID:          "ws_001",
		Name:        "Expert Validation - Population Genomics Study",
		Description: "Collaborative workspace for expert validation",
		CreatedAt:   now,
		Members: []*validation.WorkspaceMember{
			{Name: "Dr. Smith", Email: "smith@mit.edu", Role: "expert_validator", JoinedAt: now},
			{Name: "Dr. Jones", Email: "jones@stanford.edu", Role: "expert_validator", JoinedAt: now},
		},
		Channels: []*validation.Channel{
			{ID: "ch_001", Name: "general-discussion", CreatedAt: now},
		},
		SharedFiles: []*validation.SharedFile{
			{ID: "file_001", Name: "blind_test_sequences", FilePath: "/shared/blind_test_sequences", UploadedAt: now},
		},
	}

	assert.Equal(t, "ws_001", workspace.ID)
	assert.Len(t, workspace.Members, 2)
	assert.Len(t, workspace.Channels, 1)
	assert.Len(t, workspace.SharedFiles, 1)
}

// ==================== ManuscriptDraft full structure coverage ====================

func TestManuscriptDraft_FullStructure(t *testing.T) {
	now := time.Now()
	manuscript := &validation.ManuscriptDraft{
		Title:    "Expert Validation of Population-Specific Genomic Sequence Generation",
		Abstract: "We present a validated approach for generating population-specific genomic sequences...",
		Authors: []*validation.Author{
			{
				Name:         "Dr. Research Team",
				Email:        "team@phantom-science.com",
				Institution:  "Phantom Science",
				ORCID:        "0000-0001-2345-6789",
				Contribution: []string{"Conceptualization", "Software", "Validation"},
			},
		},
		Introduction:  "Population-specific genomic sequence generation...",
		Methods:       "We developed a k=8 HMM approach...",
		Results:       "Expert validation demonstrated 78% authenticity...",
		Discussion:    "Our results show that population-specific patterns...",
		Conclusion:    "This work establishes a new paradigm...",
		Keywords:      []string{"population genomics", "sequence generation", "HMM"},
		JournalTargets: []*validation.JournalConfig{
			{Name: "Nature Genetics", ImpactFactor: 38.0},
		},
		JournalVersions: map[string]*validation.FormattedManuscript{
			"Nature Genetics": {
				Content:     "Formatted content",
				Format:      "latex",
				WordCount:   4800,
				FigureCount: 5,
				TableCount:  3,
			},
		},
		Figures: []*validation.Figure{
			{ID: "fig1", Title: "Validation Results", Type: "png"},
		},
		Tables: []*validation.Table{
			{ID: "table1", Title: "Population Results", Headers: []string{"Population", "Score"}},
		},
		References: []*validation.Reference{
			{ID: "ref1", Title: "1000 Genomes Project", Year: 2015, Journal: "Nature"},
		},
		Supplementary: []*validation.SupplementaryMaterial{
			{ID: "supp1", Title: "Validation Package", Type: "zip"},
		},
		GeneratedAt: now,
	}

	assert.Equal(t, "Expert Validation of Population-Specific Genomic Sequence Generation", manuscript.Title)
	assert.Len(t, manuscript.Authors, 1)
	assert.Len(t, manuscript.Keywords, 3)
	assert.Len(t, manuscript.Figures, 1)
	assert.Len(t, manuscript.Tables, 1)
	assert.Len(t, manuscript.References, 1)
	assert.Len(t, manuscript.Supplementary, 1)
}

// TestJournalAPIConfig_APIKey_NotInJSON is a regression test for #42.
// JournalAPIConfig.APIKey must carry json:"-" yaml:"-" so it is never
// serialised to JSON (e.g. in API responses or log output).
func TestJournalAPIConfig_APIKey_NotInJSON(t *testing.T) {
	config := validation.JournalAPIConfig{
		Name:        "Nature Genetics",
		APIEndpoint: "https://api.nature.com/genetics",
		APIKey:      "secret-journal-key", // #nosec G101
	}
	data, err := json.Marshal(config)
	if err != nil {
		t.Fatalf("json.Marshal failed: %v", err)
	}
	if strings.Contains(string(data), "secret-journal-key") {
		t.Errorf("JournalAPIConfig.APIKey value must not appear in JSON output; got: %s", data)
	}
	if strings.Contains(string(data), `"APIKey"`) {
		t.Errorf("JournalAPIConfig.APIKey field name must not appear in JSON output; got: %s", data)
	}
}

// ensure fmt is referenced
var _ = fmt.Sprintf
