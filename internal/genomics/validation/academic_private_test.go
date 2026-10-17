package validation

// White-box tests for private methods in academic_integration.go and expert_orchestrator.go
// These tests use internal package access to cover methods that cannot be reached via the exported API.

import (
	"context"
	"fmt"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ==================== AcademicIntegrationAPI private method tests ====================

func makeTestAcademicIntegrationAPI(t *testing.T) *AcademicIntegrationAPI {
	t.Helper()
	config := AcademicConfig{
		InstitutionConfigs: map[string]InstitutionConfig{
			"MIT": {Name: "MIT", APIEndpoint: "https://api.mit.edu"},
		},
		ManuscriptConfig: ManuscriptConfig{
			AutoGeneration: true,
		},
	}
	api, err := NewAcademicIntegrationAPI(config)
	require.NoError(t, err)
	return api
}

func TestAcademicIntegrationAPI_rankExpertCandidates(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	candidates := []*ExpertCandidate{
		{Name: "Dr. Low Score", Publications: 5, HIndex: 2, ResearchAreas: []string{"biochemistry"}},
		{Name: "Dr. High Score", Publications: 80, HIndex: 25, ResearchAreas: []string{"population genetics", "genomics"}},
		{Name: "Dr. Mid Score", Publications: 25, HIndex: 10, ResearchAreas: []string{"genomics"}},
	}

	criteria := ExpertCriteria{
		ResearchAreas:  []string{"population genetics", "genomics"},
		ExpertiseLevel: "expert",
		Publications:   20,
	}

	ranked := api.rankExpertCandidates(candidates, criteria)
	require.Len(t, ranked, 3)
	// High-scoring candidate should come first
	assert.Equal(t, "Dr. High Score", ranked[0].Name)
	// All candidates should have match scores
	for _, c := range ranked {
		assert.GreaterOrEqual(t, c.MatchScore, 0.0)
		assert.LessOrEqual(t, c.MatchScore, 1.0)
	}
}

func TestAcademicIntegrationAPI_calculateExpertMatchScore(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	// Test with all matching criteria
	candidate := &ExpertCandidate{
		Name:          "Dr. Perfect Match",
		Institution:   "MIT",
		Publications:  100,
		HIndex:        30,
		ResearchAreas: []string{"population genetics", "genomics"},
	}

	criteria := ExpertCriteria{
		ResearchAreas:  []string{"population genetics", "genomics"},
		ExpertiseLevel: "expert",
		Publications:   50,
		Institution:    []string{"MIT", "Stanford"},
	}

	score := api.calculateExpertMatchScore(candidate, criteria)
	assert.GreaterOrEqual(t, score, 0.0)
	assert.LessOrEqual(t, score, 1.0)
	// Should be a high score (all criteria match)
	assert.Greater(t, score, 0.5)
}

func TestAcademicIntegrationAPI_calculateExpertMatchScore_NoInstitutionPref(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	// Test with no institution preference (should get full institution score)
	candidate := &ExpertCandidate{
		Name:          "Dr. Any Institution",
		Institution:   "Harvard",
		Publications:  30,
		HIndex:        12,
		ResearchAreas: []string{"genomics"},
	}

	criteria := ExpertCriteria{
		ResearchAreas:  []string{"genomics"},
		ExpertiseLevel: "senior",
		Publications:   20,
		Institution:    []string{}, // No preference
	}

	score := api.calculateExpertMatchScore(candidate, criteria)
	assert.Greater(t, score, 0.0)
}

func TestAcademicIntegrationAPI_calculateExpertMatchScore_WrongInstitution(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	// Test when institution doesn't match preferred list
	candidate := &ExpertCandidate{
		Name:          "Dr. Wrong Institution",
		Institution:   "Oxford",
		Publications:  40,
		HIndex:        15,
		ResearchAreas: []string{"population genetics"},
	}

	criteria := ExpertCriteria{
		ResearchAreas:  []string{"population genetics"},
		ExpertiseLevel: "senior",
		Publications:   30,
		Institution:    []string{"MIT", "Stanford"}, // Preferred but not Oxford
	}

	score := api.calculateExpertMatchScore(candidate, criteria)
	// Score should be lower due to institution mismatch
	assert.GreaterOrEqual(t, score, 0.0)
}

func TestAcademicIntegrationAPI_calculateResearchAreaMatch(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	// Perfect match
	score := api.calculateResearchAreaMatch(
		[]string{"population genetics", "genomics"},
		[]string{"population genetics", "genomics"},
	)
	assert.InDelta(t, 1.0, score, 0.001)

	// Partial match (1 of 2)
	score = api.calculateResearchAreaMatch(
		[]string{"population genetics", "bioinformatics"},
		[]string{"population genetics", "genomics"},
	)
	assert.InDelta(t, 0.5, score, 0.001)

	// No match
	score = api.calculateResearchAreaMatch(
		[]string{"biochemistry"},
		[]string{"population genetics", "genomics"},
	)
	assert.InDelta(t, 0.0, score, 0.001)

	// Empty criteria - should return 1.0 (no requirements)
	score = api.calculateResearchAreaMatch(
		[]string{"population genetics"},
		[]string{},
	)
	assert.InDelta(t, 1.0, score, 0.001)
}

func TestAcademicIntegrationAPI_calculateExpertiseLevelMatch(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	// Expert candidate that qualifies as expert level (pubs > 50, h-index > 15)
	expertCandidate := &ExpertCandidate{Publications: 80, HIndex: 20}
	score := api.calculateExpertiseLevelMatch(expertCandidate, "expert")
	assert.InDelta(t, 1.0, score, 0.001)

	// Senior candidate (pubs > 20, h-index > 8)
	seniorCandidate := &ExpertCandidate{Publications: 30, HIndex: 12}
	score = api.calculateExpertiseLevelMatch(seniorCandidate, "senior")
	assert.InDelta(t, 1.0, score, 0.001)

	// Junior candidate (few pubs)
	juniorCandidate := &ExpertCandidate{Publications: 5, HIndex: 2}
	score = api.calculateExpertiseLevelMatch(juniorCandidate, "junior")
	assert.InDelta(t, 1.0, score, 0.001)

	// Junior candidate required but we need expert level - low score
	score = api.calculateExpertiseLevelMatch(juniorCandidate, "expert")
	assert.Less(t, score, 1.0)
	assert.Greater(t, score, 0.0)

	// Expert candidate required but senior available - different fraction
	score = api.calculateExpertiseLevelMatch(seniorCandidate, "expert")
	// Senior (2) / expert (3) = 0.667
	assert.InDelta(t, 0.667, score, 0.01)
}

func TestAcademicIntegrationAPI_generateProjectDetails(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	criteria := ExpertCriteria{
		ResearchAreas: []string{"population genetics"},
	}

	details := api.generateProjectDetails(criteria)
	require.NotNil(t, details)
	assert.NotEmpty(t, details.Title)
	assert.NotEmpty(t, details.Description)
	assert.NotEmpty(t, details.Methodology)
	assert.NotEmpty(t, details.DataSize)
	assert.NotEmpty(t, details.Innovation)
	assert.NotEmpty(t, details.Applications)
	assert.NotEmpty(t, details.ExpectedOutcomes)
}

func TestAcademicIntegrationAPI_generatePersonalizedMessage(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	candidate := &ExpertCandidate{
		Name:          "Dr. Jane Smith",
		ResearchAreas: []string{"population genetics", "genomics"},
	}
	criteria := ExpertCriteria{
		ResearchAreas: []string{"population genetics"},
	}

	msg := api.generatePersonalizedMessage(candidate, criteria)
	assert.NotEmpty(t, msg)
	assert.Contains(t, msg, "Dr. Jane Smith")
}

func TestAcademicIntegrationAPI_generatePersonalizedInvitation(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	candidate := &ExpertCandidate{
		Name:          "Dr. Bob Kumar",
		Email:         "bkumar@stanford.edu",
		Institution:   "Stanford",
		ResearchAreas: []string{"population genomics"},
	}
	criteria := ExpertCriteria{
		ResearchAreas:  []string{"population genetics"},
		ExpertiseLevel: "senior",
	}

	invitation := api.generatePersonalizedInvitation(candidate, criteria)
	require.NotNil(t, invitation)
	assert.Equal(t, "Dr. Bob Kumar", invitation.RecipientName)
	assert.Equal(t, "bkumar@stanford.edu", invitation.RecipientEmail)
	assert.Equal(t, "Stanford", invitation.Institution)
	assert.NotEmpty(t, invitation.PersonalizedMessage)
	assert.NotNil(t, invitation.ProjectDetails)
	assert.NotEmpty(t, invitation.TimeCommitment)
	assert.NotEmpty(t, invitation.Compensation)
	assert.False(t, invitation.Deadline.IsZero())
	assert.NotNil(t, invitation.ContactInfo)
}

func TestAcademicIntegrationAPI_getContactInformation(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	contact := api.getContactInformation()
	require.NotNil(t, contact)
	assert.NotEmpty(t, contact.Name)
	assert.NotEmpty(t, contact.Email)
	assert.NotEmpty(t, contact.Institution)
}

func TestAcademicIntegrationAPI_formatInvitationEmail(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	invitation := &ExpertInvitation{
		RecipientName:       "Dr. Smith",
		PersonalizedMessage: "Dear Dr. Smith, we invite you...",
		ProjectDetails: &ProjectDetails{
			Title:       "Test Project",
			Description: "A genomics study",
		},
		TimeCommitment: "4 hours",
		Deadline:       time.Now().Add(14 * 24 * time.Hour),
	}

	emailHTML := api.formatInvitationEmail(invitation)
	assert.NotEmpty(t, emailHTML)
	assert.Contains(t, emailHTML, "Test Project")
	assert.Contains(t, emailHTML, "4 hours")
}

func TestAcademicIntegrationAPI_sendExpertInvitation(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	candidate := &ExpertCandidate{
		Name:  "Dr. Test Expert",
		Email: "test@university.edu",
	}
	invitation := &ExpertInvitation{
		RecipientName:       "Dr. Test Expert",
		RecipientEmail:      "test@university.edu",
		PersonalizedMessage: "Test message",
		ProjectDetails: &ProjectDetails{
			Title:       "Test Study",
			Description: "Study description",
		},
		TimeCommitment: "4 hours",
		Deadline:       time.Now().Add(14 * 24 * time.Hour),
	}

	ctx := context.Background()
	// sendExpertInvitation calls sendEmail which is a no-op placeholder
	err := api.sendExpertInvitation(ctx, candidate, invitation)
	assert.NoError(t, err)
}

func TestAcademicIntegrationAPI_isResearchAreaMatch(t *testing.T) {
	api := makeTestAcademicIntegrationAPI(t)

	// Exact match
	assert.True(t, api.isResearchAreaMatch("population genetics", "population genetics"))
	// No match
	assert.False(t, api.isResearchAreaMatch("biochemistry", "population genetics"))
	// Case sensitive
	assert.False(t, api.isResearchAreaMatch("Population Genetics", "population genetics"))
}

func TestFormatResearchAreas(t *testing.T) {
	// Test the package-level helper function
	// Empty
	result := formatResearchAreas([]string{})
	assert.Equal(t, "genomics research", result)

	// Single
	result = formatResearchAreas([]string{"population genetics"})
	assert.Equal(t, "population genetics", result)

	// Multiple
	result = formatResearchAreas([]string{"population genetics", "bioinformatics", "cancer genomics"})
	assert.Contains(t, result, "population genetics")
	assert.Contains(t, result, "cancer genomics")
}

func TestMinMaxFunctions(t *testing.T) {
	// Test the min/max helper functions
	assert.Equal(t, 3, min(3, 5))
	assert.Equal(t, 3, min(5, 3))
	assert.Equal(t, 3, min(3, 3))

	assert.Equal(t, 5, max(3, 5))
	assert.Equal(t, 5, max(5, 3))
	assert.Equal(t, 5, max(5, 5))
}

// ==================== Expert Orchestrator private method tests ====================

func makeOrchestratorDirect(t *testing.T) *ExpertValidationOrchestrator {
	t.Helper()
	// Create orchestrator directly without calling initialize()
	// This bypasses the FASTA file loading requirement
	config := OrchestratorConfig{
		ExpertPackagePath:      t.TempDir(),
		MaxConcurrentSessions:  5,
		SessionTimeout:         4 * time.Hour,
		MinimumExpertsRequired: 3,
		ConsensusThreshold:     0.70,
		MaxSessionDuration:     6 * time.Hour,
	}

	framework := NewExpertTestingFramework()
	// Pre-populate with test cases so framework is usable
	syntheticData := [][]byte{
		[]byte("ATCGATCGATCGATCGATCGATCGATCG"),
		[]byte("GCGCGCGCGCGCGCGCGCGCGCGCGCGC"),
	}
	realData := [][]byte{
		[]byte("ATCGCGATCGATCGATCGATCGATCGATCG"),
		[]byte("GCATTACGATCGATCGACGTATGCATGCAT"),
	}
	_ = framework.GenerateTestCases(4, syntheticData, realData)

	orchestrator := &ExpertValidationOrchestrator{
		framework:       framework,
		sessionManager:  NewSessionManager(),
		materialManager: NewValidationMaterialManager(config.ExpertPackagePath),
		reportGenerator: NewValidationReportGenerator(),
		academicIntegrator: NewAcademicIntegrator(),
		config:         config,
		activeSessions: make(map[string]*LiveValidationSession),
		orchestratorMetrics: NewOrchestratorMetrics(),
	}

	return orchestrator
}

func TestExpertValidationOrchestrator_GetActiveSessionStats_Empty(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	stats := orchestrator.GetActiveSessionStats()
	require.NotNil(t, stats)
	assert.Equal(t, 0, stats["active_sessions"])
	assert.Equal(t, 0, stats["completed_sessions"])
	assert.Equal(t, 0, stats["expired_sessions"])
	assert.Equal(t, 0, stats["total_started"])
	assert.NotNil(t, stats["phase_breakdown"])
}

func TestExpertValidationOrchestrator_StartExpertValidationSession(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	expert := &Expert{
		ID:             "expert_1",
		Name:           "Dr. Jane Smith",
		Institution:    "MIT",
		ExpertiseLevel: "expert",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"population genetics"},
		IsActive:       true,
	}

	session, err := orchestrator.StartExpertValidationSession("expert_1", expert)
	require.NoError(t, err)
	require.NotNil(t, session)

	assert.Equal(t, "expert_1", session.ExpertID)
	assert.Equal(t, expert, session.ExpertInfo)
	assert.True(t, session.IsActive)
	assert.Equal(t, PhasePreparation, session.CurrentPhase)
	assert.NotEmpty(t, session.SessionID)
	assert.NotNil(t, session.ProgressTracker)
	assert.NotNil(t, session.FeedbackChannel)
	assert.NotNil(t, session.StatusChannel)

	// Drain the status channel (welcome message)
	select {
	case status := <-session.StatusChannel:
		assert.Equal(t, "started", status.Status)
	case <-time.After(100 * time.Millisecond):
		t.Log("No initial status message received")
	}
}

func TestExpertValidationOrchestrator_StartSession_MaxLimit(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)
	orchestrator.config.MaxConcurrentSessions = 2

	expert := &Expert{
		ID:             "e1",
		ExpertiseLevel: "senior",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}

	// Start session 1
	sess1, err := orchestrator.StartExpertValidationSession("e1", expert)
	require.NoError(t, err)
	require.NotNil(t, sess1)

	// Start session 2
	expert2 := &Expert{
		ID: "e2", ExpertiseLevel: "senior",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}
	sess2, err := orchestrator.StartExpertValidationSession("e2", expert2)
	require.NoError(t, err)
	require.NotNil(t, sess2)

	// Session 3 should fail (limit reached)
	expert3 := &Expert{
		ID: "e3", ExpertiseLevel: "junior",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"bioinformatics"},
	}
	sess3, err := orchestrator.StartExpertValidationSession("e3", expert3)
	assert.Error(t, err)
	assert.Nil(t, sess3)
}

func TestExpertValidationOrchestrator_ProcessFeedback_InvalidSession(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	feedback := ExpertFeedback{
		SessionID:  "nonexistent_session",
		ExpertID:   "expert_1",
		TestCaseID: "test_1",
		Feedback: map[string]interface{}{
			"classification": "synthetic",
			"confidence":     0.8,
		},
		Phase: PhaseBlindAssessment,
	}

	err := orchestrator.ProcessExpertFeedback(feedback)
	assert.Error(t, err)
}

func TestExpertValidationOrchestrator_ProcessFeedback_WrongPhase(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	expert := &Expert{
		ID: "expert_1", ExpertiseLevel: "expert",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}

	session, err := orchestrator.StartExpertValidationSession("expert_1", expert)
	require.NoError(t, err)
	// Drain the status channel
	select {
	case <-session.StatusChannel:
	default:
	}

	// Session starts in PhasePreparation - sending feedback now should fail (wrong phase)
	feedback := ExpertFeedback{
		SessionID:  session.SessionID,
		ExpertID:   "expert_1",
		TestCaseID: "synthetic_0",
		Feedback: map[string]interface{}{
			"classification": "synthetic",
			"confidence":     0.8,
		},
		Phase: PhasePreparation, // Wrong phase for feedback
	}

	err = orchestrator.ProcessExpertFeedback(feedback)
	// PhasePreparation is the default case in switch, so should error
	assert.Error(t, err)
}

func TestExpertValidationOrchestrator_calculatePhaseProgress(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	expert := &Expert{ID: "e1", ExpertiseLevel: "senior", PerformanceHistory: []ExpertPerformance{}}
	session, err := orchestrator.StartExpertValidationSession("e1", expert)
	require.NoError(t, err)

	// Test each phase progress calculation
	session.ProgressTracker.TotalTestCases = 10
	session.ProgressTracker.CompletedCases = 5

	phases := []ValidationPhase{
		PhasePreparation,
		PhaseBlindAssessment,
		PhaseReferenceComparison,
		PhaseDiscussion,
	}

	for _, phase := range phases {
		session.CurrentPhase = phase
		progress := orchestrator.calculatePhaseProgress(session)
		assert.GreaterOrEqual(t, progress, 0.0)
		assert.LessOrEqual(t, progress, 1.0+0.001) // slight margin for float comparison
	}
}

func TestExpertValidationOrchestrator_GetActiveSessionStats_WithSessions(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	// Start multiple sessions
	for i := 1; i <= 3; i++ {
		expert := &Expert{
			ID:             fmt.Sprintf("expert_%d", i),
			ExpertiseLevel: "senior",
			PerformanceHistory: []ExpertPerformance{},
			SpecializationAreas: []string{"genomics"},
		}
		_, err := orchestrator.StartExpertValidationSession(fmt.Sprintf("expert_%d", i), expert)
		require.NoError(t, err)
	}

	stats := orchestrator.GetActiveSessionStats()
	assert.Equal(t, 3, stats["active_sessions"])
	assert.Equal(t, 3, stats["total_started"])
	phaseBreakdown := stats["phase_breakdown"].(map[string]int)
	assert.Equal(t, 3, phaseBreakdown["preparation"])
}

func TestExpertValidationOrchestrator_assessOverallQuality(t *testing.T) {
	// Test all quality branches by controlling authenticity score via framework state.
	// Authenticity = (1 - expert_accuracy) * 100
	// Branches: >= 75 = excellent, >= 60 = good, >= 45 = acceptable, else needs_improvement
	//
	// Test cases:
	//  - 80% authenticity: 2 wrong, 8 right correct->20% accuracy -> out of 10 cases, 2 correct
	//    Actually: 2 of 10 correct -> accuracy = 0.2 -> authenticity = 80
	//  - 65% authenticity: out of 20 cases, 7 correct -> accuracy = 0.35 -> authenticity = 65
	//  - 50% authenticity: out of 10 cases, 5 correct -> accuracy = 0.5 -> authenticity = 50
	//  - 30% authenticity: out of 10 cases, 7 correct -> accuracy = 0.7 -> authenticity = 30

	type qualityTest struct {
		totalCases      int
		correctCases    int
		expectedQuality string
		expectedPubReady bool
	}
	tests := []qualityTest{
		// 2 out of 10 correct -> accuracy 0.2 -> authenticity 80 -> excellent
		{10, 2, "excellent", true},
		// 7 out of 20 correct -> accuracy 0.35 -> authenticity 65 -> good
		{20, 7, "good", true},
		// 5 out of 10 correct -> accuracy 0.5 -> authenticity 50 -> acceptable
		{10, 5, "acceptable", false},
		// 7 out of 10 correct -> accuracy 0.7 -> authenticity 30 -> needs_improvement
		{10, 7, "needs_improvement", false},
	}

	for i, tt := range tests {
		t.Run(fmt.Sprintf("quality_%s", tt.expectedQuality), func(t *testing.T) {
			orchestrator := makeOrchestratorDirect(t)

			// Build test cases (all with CorrectAnswer = "synthetic")
			testCases := make([]BlindTestCase, tt.totalCases)
			for j := 0; j < tt.totalCases; j++ {
				testCases[j] = BlindTestCase{
					ID:            fmt.Sprintf("t%d_%d", i, j+1),
					CorrectAnswer: "synthetic",
				}
			}
			orchestrator.framework.testCases = testCases

			// Build consensus results: first correctCases get "synthetic" (correct),
			// remaining get "real" (incorrect)
			consensusResults := make([]ConsensusResult, tt.totalCases)
			for j := 0; j < tt.totalCases; j++ {
				answer := "real" // incorrect by default
				if j < tt.correctCases {
					answer = "synthetic" // correct
				}
				consensusResults[j] = ConsensusResult{
					TestCaseID:      fmt.Sprintf("t%d_%d", i, j+1),
					ConsensusAnswer: answer,
				}
			}
			orchestrator.framework.resultAggregator.consensusResults = consensusResults

			qa := orchestrator.assessOverallQuality()
			assert.Equal(t, tt.expectedQuality, qa.OverallQuality)
			assert.Equal(t, tt.expectedPubReady, qa.PublicationReadiness)
			assert.NotEmpty(t, qa.StrengthAreas)
			assert.NotEmpty(t, qa.ImprovementAreas)
		})
	}
}

func TestExpertValidationOrchestrator_generateRecommendations(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	recommendations := orchestrator.generateRecommendations()
	assert.NotEmpty(t, recommendations)
	assert.Greater(t, len(recommendations), 0)
}

func TestExpertValidationOrchestrator_ValidationPhaseString_Unknown(t *testing.T) {
	var unknown ValidationPhase = 99
	assert.Equal(t, "unknown", unknown.String())
}

func TestExpertValidationOrchestrator_ProcessFeedback_BlindAssessment(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	expert := &Expert{
		ID:             "expert_ba",
		ExpertiseLevel: "expert",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}

	session, err := orchestrator.StartExpertValidationSession("expert_ba", expert)
	require.NoError(t, err)

	// Move session to blind assessment phase
	session.CurrentPhase = PhaseBlindAssessment

	// Get a test case ID from the framework
	testCases := orchestrator.framework.testCases
	require.Greater(t, len(testCases), 0, "need at least one test case")
	firstCaseID := testCases[0].ID

	// Drain status channel
	select {
	case <-session.StatusChannel:
	default:
	}

	// Send blind assessment feedback
	feedback := ExpertFeedback{
		SessionID:  session.SessionID,
		ExpertID:   "expert_ba",
		TestCaseID: firstCaseID,
		Feedback: map[string]interface{}{
			"classification": "synthetic",
			"confidence":     0.85,
			"justification":  "Repetitive k-mer pattern",
		},
		Phase: PhaseBlindAssessment,
	}

	err = orchestrator.ProcessExpertFeedback(feedback)
	require.NoError(t, err)

	// Check progress was updated
	assert.Equal(t, 1, session.ProgressTracker.CompletedCases)
}

func TestExpertValidationOrchestrator_ProcessFeedback_MissingClassification(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	expert := &Expert{
		ID: "expert_mc", ExpertiseLevel: "senior",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}

	session, err := orchestrator.StartExpertValidationSession("expert_mc", expert)
	require.NoError(t, err)
	session.CurrentPhase = PhaseBlindAssessment

	testCases := orchestrator.framework.testCases
	require.Greater(t, len(testCases), 0)

	feedback := ExpertFeedback{
		SessionID:  session.SessionID,
		ExpertID:   "expert_mc",
		TestCaseID: testCases[0].ID,
		Feedback: map[string]interface{}{
			// Missing "classification"
			"confidence": 0.8,
		},
		Phase: PhaseBlindAssessment,
	}

	err = orchestrator.ProcessExpertFeedback(feedback)
	assert.Error(t, err)
}

func TestExpertValidationOrchestrator_ProcessFeedback_MissingConfidence(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	expert := &Expert{
		ID: "expert_mconf", ExpertiseLevel: "senior",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}

	session, err := orchestrator.StartExpertValidationSession("expert_mconf", expert)
	require.NoError(t, err)
	session.CurrentPhase = PhaseBlindAssessment

	testCases := orchestrator.framework.testCases
	require.Greater(t, len(testCases), 0)

	feedback := ExpertFeedback{
		SessionID:  session.SessionID,
		ExpertID:   "expert_mconf",
		TestCaseID: testCases[0].ID,
		Feedback: map[string]interface{}{
			"classification": "synthetic",
			// Missing "confidence"
		},
		Phase: PhaseBlindAssessment,
	}

	err = orchestrator.ProcessExpertFeedback(feedback)
	assert.Error(t, err)
}

func TestExpertValidationOrchestrator_ProcessFeedback_ReferenceComparison(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	expert := &Expert{
		ID: "expert_rc", ExpertiseLevel: "expert",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}

	session, err := orchestrator.StartExpertValidationSession("expert_rc", expert)
	require.NoError(t, err)
	session.CurrentPhase = PhaseReferenceComparison

	// Drain initial status messages
	for i := 0; i < 5; i++ {
		select {
		case <-session.StatusChannel:
		default:
		}
	}

	feedback := ExpertFeedback{
		SessionID:  session.SessionID,
		ExpertID:   "expert_rc",
		TestCaseID: "ref_test_1",
		Feedback: map[string]interface{}{
			"notes":           "Sequences look authentic",
			"confidence":      0.9,
			"population_bias": "CHB",
		},
		Phase: PhaseReferenceComparison,
	}

	err = orchestrator.ProcessExpertFeedback(feedback)
	require.NoError(t, err)
	// After reference comparison feedback, session should advance to discussion
	assert.Equal(t, PhaseDiscussion, session.CurrentPhase)
}

func TestExpertValidationOrchestrator_ProcessFeedback_Discussion(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	expert := &Expert{
		ID: "expert_disc", ExpertiseLevel: "expert",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}

	session, err := orchestrator.StartExpertValidationSession("expert_disc", expert)
	require.NoError(t, err)
	session.CurrentPhase = PhaseDiscussion

	// Drain initial status messages
	for i := 0; i < 5; i++ {
		select {
		case <-session.StatusChannel:
		default:
		}
	}

	feedback := ExpertFeedback{
		SessionID:  session.SessionID,
		ExpertID:   "expert_disc",
		TestCaseID: "discussion_1",
		Feedback: map[string]interface{}{
			"overall_assessment": "Excellent authenticity",
			"publication_ready":  true,
		},
		Phase: PhaseDiscussion,
	}

	// Discussion feedback completes the session
	err = orchestrator.ProcessExpertFeedback(feedback)
	// The session gets completed - may or may not error based on session completion logic
	_ = err

	// Session should no longer be active (was completed)
	assert.False(t, session.IsActive)
}

func TestExpertValidationOrchestrator_SessionMetrics(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)

	// Verify initial metrics
	assert.Equal(t, 0, orchestrator.orchestratorMetrics.TotalSessionsStarted)
	assert.Equal(t, 0, orchestrator.orchestratorMetrics.ActiveSessions)
	assert.Equal(t, 0, orchestrator.orchestratorMetrics.CompletedSessions)
	assert.Equal(t, 0, orchestrator.orchestratorMetrics.ExpiredSessions)

	// Start a session
	expert := &Expert{ID: "e_metric", ExpertiseLevel: "senior", PerformanceHistory: []ExpertPerformance{}}
	_, err := orchestrator.StartExpertValidationSession("e_metric", expert)
	require.NoError(t, err)

	assert.Equal(t, 1, orchestrator.orchestratorMetrics.TotalSessionsStarted)
	assert.Equal(t, 1, orchestrator.orchestratorMetrics.ActiveSessions)
}

// ---------- email stub functions (sendWelcomeEmail, sendCompletionEmail, sendTimeoutWarningEmail) ----------

func TestExpertValidationOrchestrator_SendWelcomeEmail(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)
	expert := &Expert{ID: "expert-email-1"}
	err := orchestrator.sendWelcomeEmail(expert, "session-123")
	if err != nil {
		t.Errorf("sendWelcomeEmail returned unexpected error: %v", err)
	}
}

func TestExpertValidationOrchestrator_SendCompletionEmail(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)
	expert := &Expert{ID: "expert-email-2"}
	report := &SessionValidationReport{}
	err := orchestrator.sendCompletionEmail(expert, "session-456", report)
	if err != nil {
		t.Errorf("sendCompletionEmail returned unexpected error: %v", err)
	}
}

func TestExpertValidationOrchestrator_SendTimeoutWarningEmail(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)
	expert := &Expert{ID: "expert-email-3"}
	err := orchestrator.sendTimeoutWarningEmail(expert, "session-789")
	if err != nil {
		t.Errorf("sendTimeoutWarningEmail returned unexpected error: %v", err)
	}
}

// ---------- generatePublicationMaterials (stub) ----------

func TestExpertValidationOrchestrator_GeneratePublicationMaterials(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)
	report := &OverallValidationReport{}
	err := orchestrator.generatePublicationMaterials(report)
	if err != nil {
		t.Errorf("generatePublicationMaterials returned unexpected error: %v", err)
	}
}

// ---------- saveOverallReport ----------

func TestExpertValidationOrchestrator_SaveOverallReport_CreatesFile(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)
	report := &OverallValidationReport{
		Timestamp:    time.Now(),
		TotalExperts: 3,
	}
	err := orchestrator.saveOverallReport(report)
	if err != nil {
		t.Errorf("saveOverallReport returned unexpected error: %v", err)
	}
}

// ---------- generateOverallValidationReport ----------

func TestExpertValidationOrchestrator_GenerateOverallValidationReport(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)
	err := orchestrator.generateOverallValidationReport()
	// May succeed or fail depending on framework state; just verify no panic
	_ = err
}

// ---------- handleSessionTimeout and handleSessionExpiry ----------

func TestExpertValidationOrchestrator_HandleSessionTimeout(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)
	expert := &Expert{
		ID:             "expert-timeout",
		ExpertiseLevel: "expert",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}
	session, err := orchestrator.StartExpertValidationSession("expert-timeout", expert)
	if err != nil {
		t.Fatalf("StartExpertValidationSession failed: %v", err)
	}
	// handleSessionTimeout sends a timeout warning; session remains active
	orchestrator.handleSessionTimeout(session.SessionID, session)
	// Drain the status channel to avoid blocking
	for len(session.StatusChannel) > 0 {
		<-session.StatusChannel
	}
	// Session should still be active after timeout warning
	if !session.IsActive {
		t.Error("handleSessionTimeout should not deactivate session (it only sends a warning)")
	}
}

func TestExpertValidationOrchestrator_HandleSessionExpiry(t *testing.T) {
	orchestrator := makeOrchestratorDirect(t)
	expert := &Expert{
		ID:             "expert-expiry",
		ExpertiseLevel: "expert",
		PerformanceHistory: []ExpertPerformance{},
		SpecializationAreas: []string{"genomics"},
	}
	session, err := orchestrator.StartExpertValidationSession("expert-expiry", expert)
	if err != nil {
		t.Fatalf("StartExpertValidationSession failed: %v", err)
	}
	// handleSessionExpiry should mark session as expired
	orchestrator.handleSessionExpiry(session.SessionID, session)
	if session.IsActive {
		t.Error("handleSessionExpiry should deactivate session")
	}
}

// TestSendStatus_NonBlocking_OverflowDoesNotBlock is a regression test for #41.
// It verifies that sending >50 status updates (the old buffer size) to a session
// whose StatusChannel has no consumer never blocks.
func TestSendStatus_NonBlocking_OverflowDoesNotBlock(t *testing.T) {
	ch := make(chan SessionStatus, 50) // same buffer size as production code

	done := make(chan struct{})
	go func() {
		defer close(done)
		// Send 100 updates — 50 beyond the buffer capacity.
		for i := 0; i < 100; i++ {
			sendStatus(ch, SessionStatus{
				SessionID: "test",
				Status:    "progress",
				Message:   fmt.Sprintf("update %d", i),
				Progress:  float64(i) / 100.0,
				Timestamp: time.Now(),
			})
		}
	}()

	select {
	case <-done:
		// All 100 sends completed without blocking — pass.
	case <-time.After(time.Second):
		t.Fatal("sendStatus blocked: goroutine did not finish within 1s (buffer overflow not handled)")
	}
}
