package validation_test

import (
	"context"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/validation"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ==================== AcademicConfig validation tests ====================

func TestAcademicConfig_Validate(t *testing.T) {
	t.Run("valid config with one institution", func(t *testing.T) {
		config := validation.AcademicConfig{
			InstitutionConfigs: map[string]validation.InstitutionConfig{
				"MIT": {Name: "MIT", APIEndpoint: "https://api.mit.edu"},
			},
		}
		err := config.Validate()
		assert.NoError(t, err)
	})

	t.Run("no institutions", func(t *testing.T) {
		config := validation.AcademicConfig{
			InstitutionConfigs: map[string]validation.InstitutionConfig{},
		}
		err := config.Validate()
		assert.Error(t, err)
	})

	t.Run("nil institutions", func(t *testing.T) {
		config := validation.AcademicConfig{}
		err := config.Validate()
		assert.Error(t, err)
	})

	t.Run("valid config with journal targets", func(t *testing.T) {
		config := validation.AcademicConfig{
			InstitutionConfigs: map[string]validation.InstitutionConfig{
				"Harvard": {Name: "Harvard"},
			},
			ManuscriptConfig: validation.ManuscriptConfig{
				TargetJournals: []validation.JournalConfig{
					{Name: "Nature", ImpactFactor: 49.9},
					{Name: "Science", ImpactFactor: 44.7},
				},
			},
		}
		err := config.Validate()
		assert.NoError(t, err)
	})

	t.Run("journal with empty name fails", func(t *testing.T) {
		config := validation.AcademicConfig{
			InstitutionConfigs: map[string]validation.InstitutionConfig{
				"Harvard": {Name: "Harvard"},
			},
			ManuscriptConfig: validation.ManuscriptConfig{
				TargetJournals: []validation.JournalConfig{
					{Name: ""}, // Empty name
				},
			},
		}
		err := config.Validate()
		assert.Error(t, err)
	})
}

func TestNewAcademicIntegrationAPI(t *testing.T) {
	t.Run("valid config", func(t *testing.T) {
		config := validation.AcademicConfig{
			InstitutionConfigs: map[string]validation.InstitutionConfig{
				"MIT": {Name: "MIT"},
			},
		}
		api, err := validation.NewAcademicIntegrationAPI(config)
		require.NoError(t, err)
		require.NotNil(t, api)
	})

	t.Run("invalid config - no institutions", func(t *testing.T) {
		config := validation.AcademicConfig{}
		api, err := validation.NewAcademicIntegrationAPI(config)
		assert.Error(t, err)
		assert.Nil(t, api)
	})

	t.Run("config with multiple institutions", func(t *testing.T) {
		config := validation.AcademicConfig{
			InstitutionConfigs: map[string]validation.InstitutionConfig{
				"MIT":      {Name: "MIT"},
				"Stanford": {Name: "Stanford"},
				"Harvard":  {Name: "Harvard"},
			},
		}
		api, err := validation.NewAcademicIntegrationAPI(config)
		require.NoError(t, err)
		require.NotNil(t, api)
	})
}

func TestAcademicIntegrationAPI_RecruitExperts(t *testing.T) {
	config := validation.AcademicConfig{
		InstitutionConfigs: map[string]validation.InstitutionConfig{
			"MIT": {Name: "MIT"},
		},
	}
	api, err := validation.NewAcademicIntegrationAPI(config)
	require.NoError(t, err)

	criteria := validation.ExpertCriteria{
		ResearchAreas:  []string{"population genetics", "genomics"},
		ExpertiseLevel: "senior",
		Publications:   20,
	}

	ctx := context.Background()
	candidates, err := api.RecruitExperts(ctx, criteria)
	require.NoError(t, err)
	// The placeholder InstitutionAPI returns empty slice; candidates may be nil or empty
	// Just verify no error occurs
	_ = candidates
}

func TestAcademicIntegrationAPI_TrackSubmissionStatus(t *testing.T) {
	config := validation.AcademicConfig{
		InstitutionConfigs: map[string]validation.InstitutionConfig{
			"MIT": {Name: "MIT"},
		},
	}
	api, err := validation.NewAcademicIntegrationAPI(config)
	require.NoError(t, err)

	ctx := context.Background()
	// Should error since submission is not found
	status, err := api.TrackSubmissionStatus(ctx, "nonexistent_submission")
	assert.Error(t, err)
	assert.Nil(t, status)
}

func TestAcademicIntegrationAPI_SetupCollaborationWorkspace(t *testing.T) {
	config := validation.AcademicConfig{
		InstitutionConfigs: map[string]validation.InstitutionConfig{
			"MIT": {Name: "MIT"},
		},
	}
	api, err := validation.NewAcademicIntegrationAPI(config)
	require.NoError(t, err)

	experts := []*validation.ExpertCandidate{
		{
			Name:          "Dr. Jane Smith",
			Email:         "jsmith@mit.edu",
			Institution:   "MIT",
			ResearchAreas: []string{"population genetics"},
			Publications:  45,
			HIndex:        18,
		},
	}

	projectDetails := &validation.ProjectDetails{
		Title:       "Expert Validation Study",
		Description: "Validation of population-specific genomic sequences",
		Methodology: "Blind assessment and expert panel review",
	}

	ctx := context.Background()
	workspace, err := api.SetupCollaborationWorkspace(ctx, experts, projectDetails)
	require.NoError(t, err)
	require.NotNil(t, workspace)
	assert.Contains(t, workspace.Name, projectDetails.Title)
	assert.Greater(t, len(workspace.Channels), 0)
	// Should have at least some members (depends on InviteMember implementation)
	assert.NotNil(t, workspace.Members)
}

func TestAcademicIntegrationAPI_GenerateManuscript(t *testing.T) {
	config := validation.AcademicConfig{
		InstitutionConfigs: map[string]validation.InstitutionConfig{
			"MIT": {Name: "MIT"},
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

	report := &validation.OverallValidationReport{
		TotalExperts:      10,
		AuthenticityScore: 82.5,
		Timestamp:         time.Now(),
	}

	ctx := context.Background()
	manuscript, err := api.GenerateManuscript(ctx, report)
	require.NoError(t, err)
	require.NotNil(t, manuscript)
	assert.NotEmpty(t, manuscript.Title)
	assert.NotEmpty(t, manuscript.Abstract)
	assert.NotNil(t, manuscript.Authors)
	assert.NotEmpty(t, manuscript.Keywords)
	assert.NotNil(t, manuscript.Figures)
	assert.NotNil(t, manuscript.Tables)
	assert.NotNil(t, manuscript.References)
	// With target journals configured and a high authenticity score (82.5 >= 75),
	// GenerateManuscript should populate JournalVersions for up to 3 journals.
	require.NotNil(t, manuscript.JournalVersions, "JournalVersions map must be initialized")
	assert.Len(t, manuscript.JournalTargets, 3, "high score selects 3 journals")
	for _, journal := range manuscript.JournalTargets {
		_, ok := manuscript.JournalVersions[journal.Name]
		assert.True(t, ok, "JournalVersions should have entry for %s", journal.Name)
	}
}

// ==================== OrchestratorConfig tests ====================

func TestOrchestratorConfig_Validate(t *testing.T) {
	t.Run("valid config with package path", func(t *testing.T) {
		config := validation.OrchestratorConfig{
			ExpertPackagePath: "/tmp/expert_package",
		}
		err := config.Validate()
		assert.NoError(t, err)
		// Defaults should be set
		assert.Equal(t, 5, config.MaxConcurrentSessions)
		assert.Equal(t, 3, config.MinimumExpertsRequired)
		assert.Equal(t, 0.70, config.ConsensusThreshold)
	})

	t.Run("missing expert package path", func(t *testing.T) {
		config := validation.OrchestratorConfig{
			ExpertPackagePath: "",
		}
		err := config.Validate()
		assert.Error(t, err)
	})

	t.Run("custom values are preserved", func(t *testing.T) {
		config := validation.OrchestratorConfig{
			ExpertPackagePath:      "/tmp/package",
			MaxConcurrentSessions:  10,
			SessionTimeout:         2 * time.Hour,
			MinimumExpertsRequired: 5,
			ConsensusThreshold:     0.80,
			MaxSessionDuration:     8 * time.Hour,
		}
		err := config.Validate()
		assert.NoError(t, err)
		assert.Equal(t, 10, config.MaxConcurrentSessions)
		assert.Equal(t, 5, config.MinimumExpertsRequired)
		assert.Equal(t, 0.80, config.ConsensusThreshold)
	})
}

func TestNewExpertValidationOrchestrator(t *testing.T) {
	t.Run("missing expert package path fails", func(t *testing.T) {
		config := validation.OrchestratorConfig{
			ExpertPackagePath: "",
		}
		orchestrator, err := validation.NewExpertValidationOrchestrator(config)
		assert.Error(t, err)
		assert.Nil(t, orchestrator)
	})

	t.Run("valid path but no sequence data fails", func(t *testing.T) {
		// The orchestrator initialization requires both synthetic and real sequence data.
		// Without pre-loaded FASTA files in the right structure, initialization fails.
		tmpDir := t.TempDir()
		require.NoError(t, os.MkdirAll(filepath.Join(tmpDir, "blind_test_sequences"), 0750))
		require.NoError(t, os.MkdirAll(filepath.Join(tmpDir, "reference_sequences"), 0750))

		config := validation.OrchestratorConfig{
			ExpertPackagePath:      tmpDir,
			MaxConcurrentSessions:  5,
			MinimumExpertsRequired: 3,
			ConsensusThreshold:     0.70,
		}
		// This will fail because no FASTA data is present to generate test cases
		orchestrator, err := validation.NewExpertValidationOrchestrator(config)
		assert.Error(t, err)
		assert.Nil(t, orchestrator)
	})
}

// ==================== ValidationPhase.String() additional tests ====================

func TestValidationPhase_PhaseCompletion(t *testing.T) {
	// Test specifically that PhaseCompletion returns "completion"
	assert.Equal(t, "completion", validation.PhaseCompletion.String())
}

// ==================== Supporting types tests ====================

func TestSessionProgressTracker_Fields(t *testing.T) {
	tracker := validation.NewSessionProgressTracker()
	require.NotNil(t, tracker)
	assert.Equal(t, 0, tracker.TotalTestCases)
	assert.Equal(t, 0, tracker.CompletedCases)
}

func TestOrchestratorMetrics_New(t *testing.T) {
	metrics := validation.NewOrchestratorMetrics()
	require.NotNil(t, metrics)
	assert.Equal(t, 0, metrics.TotalSessionsStarted)
	assert.Equal(t, 0, metrics.ActiveSessions)
	assert.Equal(t, 0, metrics.CompletedSessions)
	assert.NotNil(t, metrics.ExpertParticipation)
}

func TestSessionManager_New(t *testing.T) {
	sm := validation.NewSessionManager()
	require.NotNil(t, sm)
}

func TestValidationReportGenerator_New(t *testing.T) {
	gen := validation.NewValidationReportGenerator()
	require.NotNil(t, gen)
}

func TestAcademicIntegrator_New(t *testing.T) {
	integrator := validation.NewAcademicIntegrator()
	require.NotNil(t, integrator)
}

func TestAcademicIntegrator_SubmitValidationResultsOK(t *testing.T) {
	integrator := validation.NewAcademicIntegrator()
	report := &validation.OverallValidationReport{
		TotalExperts:      10,
		AuthenticityScore: 80.0,
	}
	err := integrator.SubmitValidationResults(report)
	assert.NoError(t, err)
}

func TestExpertCriteria_Fields(t *testing.T) {
	criteria := validation.ExpertCriteria{
		ResearchAreas:  []string{"population genetics", "computational biology"},
		ExpertiseLevel: "senior",
		Publications:   25,
		ExperienceYears: 10,
		RemoteAvailable: true,
		Languages:      []string{"en", "de"},
	}
	assert.Equal(t, "senior", criteria.ExpertiseLevel)
	assert.Equal(t, 25, criteria.Publications)
	assert.True(t, criteria.RemoteAvailable)
}

func TestExpertCandidate_Fields(t *testing.T) {
	candidate := &validation.ExpertCandidate{
		Name:           "Dr. Bob Kumar",
		Email:          "bkumar@stanford.edu",
		Institution:    "Stanford",
		Department:     "Genetics",
		ResearchAreas:  []string{"population genomics"},
		Publications:   60,
		HIndex:         22,
		ExpertiseScore: 0.92,
		MatchScore:     0.88,
		ORCID:          "0000-0001-2345-6789",
	}
	assert.Equal(t, "Dr. Bob Kumar", candidate.Name)
	assert.Equal(t, 60, candidate.Publications)
	assert.Equal(t, 22, candidate.HIndex)
	assert.InDelta(t, 0.88, candidate.MatchScore, 0.001)
}

func TestPublicationWorkflowManager_NoJournals_Config(t *testing.T) {
	// Test that config with no journal APIs fails validation
	config := validation.PublicationConfig{
		WorkingDirectory: "/tmp/work",
		OutputDirectory:  "/tmp/output",
		JournalAPIs:      map[string]validation.JournalAPIConfig{},
	}
	manager, err := validation.NewPublicationWorkflowManager(config)
	assert.Error(t, err)
	assert.Nil(t, manager)
}

func TestWorkflowStatus_AllValues(t *testing.T) {
	// Test all known WorkflowStatus values
	assert.Equal(t, "active", validation.StatusActive.String())
	assert.Equal(t, "completed", validation.StatusCompleted.String())
	assert.Equal(t, "failed", validation.StatusFailed.String())
	assert.Equal(t, "cancelled", validation.StatusCancelled.String())
}

func TestWorkflowStage_AllValues(t *testing.T) {
	// Test all known WorkflowStage values
	assert.Equal(t, "initialization", validation.StageInitialization.String())
	assert.Equal(t, "manuscript_generation", validation.StageManuscriptGeneration.String())
	assert.Equal(t, "compliance_check", validation.StageComplianceCheck.String())
	assert.Equal(t, "peer_review", validation.StagePeerReview.String())
	assert.Equal(t, "submission_preparation", validation.StageSubmissionPreparation.String())
	assert.Equal(t, "journal_submission", validation.StageJournalSubmission.String())
	assert.Equal(t, "submission_tracking", validation.StageSubmissionTracking.String())
}
