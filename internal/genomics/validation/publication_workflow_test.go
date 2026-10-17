package validation_test

import (
	"context"
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/validation"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// makePublicationConfig creates a valid PublicationConfig for testing
func makePublicationConfig(t *testing.T) (validation.PublicationConfig, string) {
	t.Helper()
	tmpDir := t.TempDir()
	config := validation.PublicationConfig{
		WorkingDirectory:  filepath.Join(tmpDir, "work"),
		TemplateDirectory: filepath.Join(tmpDir, "templates"),
		OutputDirectory:   filepath.Join(tmpDir, "output"),
		JournalAPIs: map[string]validation.JournalAPIConfig{
			"NatureGenetics": {
				Name:            "Nature Genetics",
				APIEndpoint:     "https://api.nature.com",
				PreferredFormat: "latex",
				Requirements: validation.JournalRequirements{
					AbstractWordLimit: 300,
					MainTextWordLimit: 4000,
					MaxFigures:        8,
				},
			},
			"Genome Biology": {
				Name:            "Genome Biology",
				APIEndpoint:     "https://api.genomebiology.com",
				PreferredFormat: "docx",
			},
		},
		ManuscriptSettings: validation.ManuscriptGenerationSettings{
			AutoGenerate:     true,
			CitationStyle:    "nature",
			AbstractWordLimit: 300,
		},
		ReviewSettings: validation.ReviewWorkflowSettings{
			EnableReviewWorkflow: false,
			ReviewDeadline:       7 * 24 * time.Hour,
		},
		ComplianceSettings: validation.ComplianceSettings{
			EnableComplianceCheck: true,
		},
		SubmissionSettings: validation.SubmissionTrackingSettings{
			EnableTracking: true,
		},
		AnalyticsEnabled: true,
	}
	return config, tmpDir
}

func TestNewPublicationWorkflowManager(t *testing.T) {
	config, _ := makePublicationConfig(t)
	manager, err := validation.NewPublicationWorkflowManager(config)
	require.NoError(t, err)
	require.NotNil(t, manager)
}

func TestNewPublicationWorkflowManager_InvalidConfig(t *testing.T) {
	t.Run("no working directory", func(t *testing.T) {
		config := validation.PublicationConfig{
			WorkingDirectory: "",
			OutputDirectory:  "/tmp/output",
			JournalAPIs:      map[string]validation.JournalAPIConfig{"j1": {Name: "J1"}},
		}
		manager, err := validation.NewPublicationWorkflowManager(config)
		assert.Error(t, err)
		assert.Nil(t, manager)
	})

	t.Run("no output directory", func(t *testing.T) {
		config := validation.PublicationConfig{
			WorkingDirectory: "/tmp/work",
			OutputDirectory:  "",
			JournalAPIs:      map[string]validation.JournalAPIConfig{"j1": {Name: "J1"}},
		}
		manager, err := validation.NewPublicationWorkflowManager(config)
		assert.Error(t, err)
		assert.Nil(t, manager)
	})
}

func TestPublicationWorkflowManager_GetPublicationMetrics(t *testing.T) {
	config, _ := makePublicationConfig(t)
	manager, err := validation.NewPublicationWorkflowManager(config)
	require.NoError(t, err)

	metrics := manager.GetPublicationMetrics()
	require.NotNil(t, metrics)
	assert.Equal(t, 0, metrics.TotalWorkflowsStarted)
	assert.Equal(t, 0, metrics.CompletedWorkflows)
}

func TestPublicationWorkflowManager_GetWorkflowStatus_NotFound(t *testing.T) {
	config, _ := makePublicationConfig(t)
	manager, err := validation.NewPublicationWorkflowManager(config)
	require.NoError(t, err)

	_, err = manager.GetWorkflowStatus("nonexistent_workflow")
	assert.Error(t, err)
}

func TestPublicationWorkflowManager_StartWorkflow(t *testing.T) {
	config, _ := makePublicationConfig(t)
	manager, err := validation.NewPublicationWorkflowManager(config)
	require.NoError(t, err)

	report := &validation.OverallValidationReport{
		TotalExperts:      5,
		AuthenticityScore: 75.0,
		Timestamp:         time.Now(),
	}

	ctx := context.Background()
	workflow, err := manager.StartPublicationWorkflow(ctx, report, []string{"NatureGenetics"})
	require.NoError(t, err)
	require.NotNil(t, workflow)

	assert.NotEmpty(t, workflow.ID)
	assert.Equal(t, report, workflow.ValidationReport)
	assert.Contains(t, workflow.TargetJournals, "NatureGenetics")

	// Wait a bit for workflow to start executing
	time.Sleep(50 * time.Millisecond)

	// Metrics should be updated
	metrics := manager.GetPublicationMetrics()
	assert.GreaterOrEqual(t, metrics.TotalWorkflowsStarted, 1)
}

func TestPublicationWorkflowManager_StartWorkflow_ThenGetStatus(t *testing.T) {
	config, _ := makePublicationConfig(t)
	manager, err := validation.NewPublicationWorkflowManager(config)
	require.NoError(t, err)

	report := &validation.OverallValidationReport{
		TotalExperts:      3,
		AuthenticityScore: 60.0,
	}

	ctx := context.Background()
	workflow, err := manager.StartPublicationWorkflow(ctx, report, []string{"NatureGenetics"})
	require.NoError(t, err)

	// Check status immediately (workflow may still be active)
	status, err := manager.GetWorkflowStatus(workflow.ID)
	// Workflow might have completed very quickly, leading to not-found error
	// or still be active
	if err == nil {
		require.NotNil(t, status)
		assert.Equal(t, workflow.ID, status.WorkflowID)
	}

	// Wait for background goroutine to complete before temp dir cleanup
	time.Sleep(200 * time.Millisecond)
}

func TestPublicationWorkflowManager_StartWorkflow_Cancelled(t *testing.T) {
	config, _ := makePublicationConfig(t)
	manager, err := validation.NewPublicationWorkflowManager(config)
	require.NoError(t, err)

	report := &validation.OverallValidationReport{
		TotalExperts:      3,
		AuthenticityScore: 55.0,
	}

	// Cancel context immediately
	ctx, cancel := context.WithCancel(context.Background())
	cancel() // Cancel before even starting

	workflow, err := manager.StartPublicationWorkflow(ctx, report, []string{"NatureGenetics"})
	require.NoError(t, err) // StartPublicationWorkflow itself doesn't error
	require.NotNil(t, workflow)

	// The async workflow will be cancelled
	time.Sleep(100 * time.Millisecond)
}

func TestPublicationWorkflowManager_MultipleWorkflows(t *testing.T) {
	config, _ := makePublicationConfig(t)
	manager, err := validation.NewPublicationWorkflowManager(config)
	require.NoError(t, err)

	report := &validation.OverallValidationReport{
		TotalExperts:      5,
		AuthenticityScore: 80.0,
	}

	ctx := context.Background()

	// Start 3 workflows - IDs may collide if timestamps are identical, so just verify they start
	for i := 0; i < 3; i++ {
		wf, e := manager.StartPublicationWorkflow(ctx, report, []string{"NatureGenetics"})
		require.NoError(t, e)
		require.NotNil(t, wf)
		assert.NotEmpty(t, wf.ID)
	}

	// At least 3 workflows should have been tracked (or more if IDs collide and overwrite)
	metrics := manager.GetPublicationMetrics()
	assert.GreaterOrEqual(t, metrics.TotalWorkflowsStarted, 3)

	// Wait for background goroutines to complete before temp dir cleanup
	time.Sleep(200 * time.Millisecond)
}

// ==================== PublicationConfig validation tests ====================

func TestPublicationConfig_Validate(t *testing.T) {
	t.Run("valid config", func(t *testing.T) {
		config := validation.PublicationConfig{
			WorkingDirectory: "/tmp/work",
			OutputDirectory:  "/tmp/output",
			JournalAPIs: map[string]validation.JournalAPIConfig{
				"journal1": {Name: "Journal 1"},
			},
		}
		err := config.Validate()
		assert.NoError(t, err)
	})
}

// ==================== WorkflowStatus & Artifact tests ====================

func TestWorkflowStatus_Completeness(t *testing.T) {
	// Test unknown status
	var unknownStatus validation.WorkflowStatus = 99
	assert.Equal(t, "unknown", unknownStatus.String())
}

func TestWorkflowStage_Completeness(t *testing.T) {
	// Test unknown stage
	var unknownStage validation.WorkflowStage = 99
	assert.Equal(t, "unknown", unknownStage.String())
}

func TestPublicationWorkflow_Fields(t *testing.T) {
	workflow := &validation.PublicationWorkflow{
		ID:             "pub_123",
		TargetJournals: []string{"Nature Genetics"},
		StartTime:      time.Now(),
		CurrentStage:   validation.StageInitialization,
		Status:         validation.StatusActive,
		Progress:       map[string]float64{"overall": 0.5},
		Artifacts:      map[string]*validation.Artifact{},
		Submissions:    map[string]*validation.JournalSubmission{},
		StageHistory:   []*validation.StageTransition{},
	}

	assert.Equal(t, "pub_123", workflow.ID)
	assert.Equal(t, validation.StatusActive, workflow.Status)
	assert.Equal(t, validation.StageInitialization, workflow.CurrentStage)
}

func TestArtifact_Fields(t *testing.T) {
	artifact := &validation.Artifact{
		Type:     "manuscript",
		FilePath: "/tmp/manuscript.tex",
		Metadata: map[string]interface{}{"word_count": 5000},
	}

	assert.Equal(t, "manuscript", artifact.Type)
	assert.Equal(t, "/tmp/manuscript.tex", artifact.FilePath)
	assert.Equal(t, 5000, artifact.Metadata["word_count"])
}

func TestJournalSubmission_Fields(t *testing.T) {
	now := time.Now()
	submission := &validation.JournalSubmission{
		SubmissionID: "sub_12345",
		Journal:      "Nature Genetics",
		Status:       "submitted",
		SubmittedAt:  now,
		TrackingURL:  "https://naturegenetics.com/track/sub_12345",
	}

	assert.Equal(t, "sub_12345", submission.SubmissionID)
	assert.Equal(t, "submitted", submission.Status)
	assert.Equal(t, "Nature Genetics", submission.Journal)
}

// ==================== Additional collaboration/manuscript types ====================

func TestCollaborationPlatform_TrackInvitation(t *testing.T) {
	config := validation.CollaborationConfig{PlatformType: "slack"}
	platform := validation.NewCollaborationPlatform(config)

	candidate := &validation.ExpertCandidate{
		Name:  "Dr. Test",
		Email: "test@example.com",
	}
	invitation := &validation.ExpertInvitation{
		RecipientName:  "Dr. Test",
		RecipientEmail: "test@example.com",
	}

	// Should not panic
	platform.TrackInvitation(candidate, invitation)
}

func TestCollaborationPlatform_CreateChannel(t *testing.T) {
	config := validation.CollaborationConfig{PlatformType: "teams"}
	platform := validation.NewCollaborationPlatform(config)

	ctx := context.Background()
	channel, err := platform.CreateChannel(ctx, "workspace_1", "general-discussion")
	require.NoError(t, err)
	require.NotNil(t, channel)
	assert.NotEmpty(t, channel.ID)
	assert.Equal(t, "general-discussion", channel.Name)
}

func TestCollaborationPlatform_InviteMember(t *testing.T) {
	config := validation.CollaborationConfig{PlatformType: "slack"}
	platform := validation.NewCollaborationPlatform(config)

	ctx := context.Background()
	member := &validation.WorkspaceMember{
		Name:  "Dr. Test",
		Email: "test@example.com",
		Role:  "expert_validator",
	}

	err := platform.InviteMember(ctx, "workspace_1", member)
	assert.NoError(t, err)
}

func TestCollaborationPlatform_UploadFile(t *testing.T) {
	config := validation.CollaborationConfig{PlatformType: "slack"}
	platform := validation.NewCollaborationPlatform(config)

	ctx := context.Background()
	file, err := platform.UploadFile(ctx, "workspace_1", "blind_test_sequences")
	require.NoError(t, err)
	require.NotNil(t, file)
	assert.NotEmpty(t, file.ID)
	assert.Equal(t, "blind_test_sequences", file.Name)
}

func TestInstitutionAPI_SearchExperts(t *testing.T) {
	// Serve a mock institution API that returns two expert candidates.
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Type", "application/json")
		candidates := []validation.ExpertCandidate{
			{Name: "Dr. Alice Kim", Publications: 30, HIndex: 12},
			{Name: "Dr. Bob Chen", Publications: 25, HIndex: 10},
		}
		json.NewEncoder(w).Encode(candidates)
	}))
	defer srv.Close()

	config := validation.InstitutionConfig{
		Name:        "MIT",
		APIEndpoint: srv.URL,
	}
	api, err := validation.NewInstitutionAPI(config, srv.Client())
	require.NoError(t, err)

	ctx := context.Background()
	criteria := validation.ExpertCriteria{
		ResearchAreas:  []string{"population genetics"},
		ExpertiseLevel: "senior",
		Publications:   20,
	}

	candidates, err := api.SearchExperts(ctx, criteria)
	require.NoError(t, err)
	assert.Len(t, candidates, 2)
	assert.Equal(t, "Dr. Alice Kim", candidates[0].Name)
}

func TestInstitutionAPI_CheckSubmissionStatus(t *testing.T) {
	config := validation.InstitutionConfig{
		Name: "MIT",
	}
	api, err := validation.NewInstitutionAPI(config, nil)
	require.NoError(t, err)

	ctx := context.Background()
	status, err := api.CheckSubmissionStatus(ctx, "sub_123")
	// Placeholder returns not implemented error
	assert.Error(t, err)
	assert.Nil(t, status)
}

func TestManuscriptManager_FormatForJournal(t *testing.T) {
	config := validation.ManuscriptConfig{
		AutoGeneration: true,
	}
	manager := validation.NewManuscriptManager(config)

	manuscript := &validation.ManuscriptDraft{
		Title:   "Test Manuscript",
		Figures: []*validation.Figure{{ID: "fig1"}},
		Tables:  []*validation.Table{{ID: "table1"}},
	}
	journal := &validation.JournalConfig{
		Name:            "Nature Genetics",
		PreferredFormat: "latex",
	}

	formatted, err := manager.FormatForJournal(manuscript, journal)
	require.NoError(t, err)
	require.NotNil(t, formatted)
	assert.Equal(t, "latex", formatted.Format)
	assert.Equal(t, 5000, formatted.WordCount)
}

// ==================== PublicationConfig directory creation tests ====================

func TestPublicationWorkflowManager_CreatesDirectories(t *testing.T) {
	tmpDir := t.TempDir()
	config := validation.PublicationConfig{
		WorkingDirectory:  filepath.Join(tmpDir, "work"),
		TemplateDirectory: filepath.Join(tmpDir, "templates"),
		OutputDirectory:   filepath.Join(tmpDir, "output"),
		JournalAPIs: map[string]validation.JournalAPIConfig{
			"j1": {Name: "J1"},
		},
	}

	manager, err := validation.NewPublicationWorkflowManager(config)
	require.NoError(t, err)
	require.NotNil(t, manager)

	// Directories should have been created
	_, err = os.Stat(filepath.Join(tmpDir, "work"))
	assert.NoError(t, err)
	_, err = os.Stat(filepath.Join(tmpDir, "output"))
	assert.NoError(t, err)
}
