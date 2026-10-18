package validation

import (
	"context"
	"encoding/json"
	"fmt"
	"net/smtp"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/errors"
	"github.com/scttfrdmn/phantom-science/internal/security"
)

// ExpertValidationOrchestrator runs expert blind-assessment sessions in
// memory: it builds test cases from blind/ and reference/ sequence files,
// records feedback through the session phases, enforces timeouts, optionally
// sends SMTP notifications, and writes JSON reports. There is no web portal
// or remote API; feedback is passed in by the caller.
type ExpertValidationOrchestrator struct {
	framework          *ExpertTestingFramework
	sessionManager     *SessionManager
	materialManager    *ValidationMaterialManager
	reportGenerator    *ValidationReportGenerator
	academicIntegrator *AcademicIntegrator

	// Configuration
	config OrchestratorConfig

	// State management
	activeSessions map[string]*LiveValidationSession
	sessionMutex   sync.RWMutex

	// Background goroutine lifecycle
	ctx    context.Context
	cancel context.CancelFunc
	wg     sync.WaitGroup

	// Metrics and monitoring
	orchestratorMetrics *OrchestratorMetrics
}

// OrchestratorConfig provides configuration for the validation orchestrator
type OrchestratorConfig struct {
	ExpertPackagePath     string
	MaxConcurrentSessions int
	SessionTimeout        time.Duration
	AutoReportGeneration  bool
	AcademicIntegration   bool

	// Email and communication settings
	SMTPConfig         SMTPConfig
	EmailNotifications bool

	// Validation parameters
	MinimumExpertsRequired int
	ConsensusThreshold     float64
	MaxSessionDuration     time.Duration
}

// SMTPConfig contains email configuration
type SMTPConfig struct {
	Host     string
	Port     int
	Username string
	Password string `json:"-" yaml:"-"` // #nosec G117
	FromAddr string
}

// LiveValidationSession represents an active expert validation session
type LiveValidationSession struct {
	SessionID      string
	ExpertID       string
	ExpertInfo     *Expert
	StartTime      time.Time
	CurrentPhase   ValidationPhase
	CompletedTests []string

	// Session state
	IsActive     bool
	LastActivity time.Time
	SessionData  map[string]interface{}

	// Progress tracking
	ProgressTracker *SessionProgressTracker

	// Communication
	FeedbackChannel chan ExpertFeedback
	StatusChannel   chan SessionStatus
}

// ValidationPhase represents different phases of expert validation
type ValidationPhase int

const (
	PhasePreparation ValidationPhase = iota
	PhaseBlindAssessment
	PhaseReferenceComparison
	PhaseDiscussion
	PhaseCompletion
)

// SessionProgressTracker tracks session progress and performance
type SessionProgressTracker struct {
	TotalTestCases      int
	CompletedCases      int
	CurrentCaseID       string
	PhaseStartTime      time.Time
	EstimatedCompletion time.Time

	// Performance metrics
	AverageResponseTime time.Duration
	CurrentAccuracy     float64
	ConsistencyScore    float64
}

// ExpertFeedback represents real-time feedback from experts
type ExpertFeedback struct {
	SessionID  string
	ExpertID   string
	TestCaseID string
	Feedback   map[string]interface{}
	Timestamp  time.Time
	Phase      ValidationPhase
}

// SessionStatus represents session status updates
type SessionStatus struct {
	SessionID string
	Status    string
	Message   string
	Progress  float64
	Timestamp time.Time
}

// NewExpertValidationOrchestrator creates a new expert validation orchestrator
func NewExpertValidationOrchestrator(config OrchestratorConfig) (*ExpertValidationOrchestrator, error) {
	// Validate configuration
	if err := config.Validate(); err != nil {
		return nil, errors.WrapError(err, errors.ConfigurationError, "expert_orchestrator", "new_orchestrator", "invalid configuration")
	}

	ctx, cancel := context.WithCancel(context.Background())
	orchestrator := &ExpertValidationOrchestrator{
		framework:           NewExpertTestingFramework(),
		sessionManager:      NewSessionManager(),
		materialManager:     NewValidationMaterialManager(config.ExpertPackagePath),
		reportGenerator:     NewValidationReportGenerator(),
		academicIntegrator:  NewAcademicIntegrator(),
		config:              config,
		activeSessions:      make(map[string]*LiveValidationSession),
		orchestratorMetrics: NewOrchestratorMetrics(),
		ctx:                 ctx,
		cancel:              cancel,
	}

	// Initialize components
	if err := orchestrator.initialize(); err != nil {
		return nil, errors.WrapError(err, errors.InitializationError, "expert_orchestrator", "new_orchestrator", "failed to initialize components")
	}

	return orchestrator, nil
}

// Validate validates the orchestrator configuration
func (config *OrchestratorConfig) Validate() error {
	if config.ExpertPackagePath == "" {
		return fmt.Errorf("expert package path is required")
	}

	if config.MaxConcurrentSessions <= 0 {
		config.MaxConcurrentSessions = 5 // Default
	}

	if config.SessionTimeout <= 0 {
		config.SessionTimeout = 4 * time.Hour // Default 4 hours
	}

	if config.MinimumExpertsRequired <= 0 {
		config.MinimumExpertsRequired = 3 // Default minimum 3 experts
	}

	if config.ConsensusThreshold <= 0 {
		config.ConsensusThreshold = 0.70 // Default 70% consensus
	}

	if config.MaxSessionDuration <= 0 {
		config.MaxSessionDuration = 6 * time.Hour // Default max 6 hours
	}

	return nil
}

// initialize sets up the orchestrator components
func (evo *ExpertValidationOrchestrator) initialize() error {
	// Load existing expert validation materials
	if err := evo.materialManager.LoadValidationMaterials(); err != nil {
		return fmt.Errorf("failed to load validation materials: %w", err)
	}

	// Initialize expert testing framework with materials
	testCases, err := evo.materialManager.GetTestCases()
	if err != nil {
		return fmt.Errorf("failed to get test cases: %w", err)
	}

	// Load synthetic and real data for test case generation
	syntheticData, realData, err := evo.materialManager.LoadSequenceData()
	if err != nil {
		return fmt.Errorf("failed to load sequence data: %w", err)
	}

	// Use the number of pre-defined test cases, or fall back to the number of available
	// data files when no test cases have been explicitly defined.
	numCases := len(testCases)
	if numCases == 0 {
		numCases = len(syntheticData) + len(realData)
	}
	if err := evo.framework.GenerateTestCases(numCases, syntheticData, realData); err != nil {
		return fmt.Errorf("failed to generate test cases: %w", err)
	}

	// Start background monitoring
	evo.wg.Add(1)
	go evo.monitorSessions()

	return nil
}

// StartExpertValidationSession starts a new expert validation session
func (evo *ExpertValidationOrchestrator) StartExpertValidationSession(expertID string, expertInfo *Expert) (*LiveValidationSession, error) {
	evo.sessionMutex.Lock()
	defer evo.sessionMutex.Unlock()

	// Check concurrent session limit
	if len(evo.activeSessions) >= evo.config.MaxConcurrentSessions {
		return nil, errors.NewResourceError("expert_orchestrator", "start_session", "maximum concurrent sessions reached")
	}

	// Generate session ID
	sessionID := fmt.Sprintf("session_%s_%d", expertID, time.Now().Unix())

	// Create session
	session := &LiveValidationSession{
		SessionID:       sessionID,
		ExpertID:        expertID,
		ExpertInfo:      expertInfo,
		StartTime:       time.Now(),
		CurrentPhase:    PhasePreparation,
		CompletedTests:  make([]string, 0),
		IsActive:        true,
		LastActivity:    time.Now(),
		SessionData:     make(map[string]interface{}),
		ProgressTracker: NewSessionProgressTracker(),
		FeedbackChannel: make(chan ExpertFeedback, 100),
		StatusChannel:   make(chan SessionStatus, 50),
	}

	// Initialize progress tracker
	testCases := evo.framework.testCases
	session.ProgressTracker.TotalTestCases = len(testCases)

	// Register expert with framework
	evo.framework.RegisterExpert(expertInfo)

	// Store active session
	evo.activeSessions[sessionID] = session

	// Send welcome notification
	if evo.config.EmailNotifications {
		evo.wg.Add(1)
		go func() {
			defer evo.wg.Done()
			if err := evo.sendWelcomeEmail(expertInfo, sessionID); err != nil {
				fmt.Printf("warning: failed to send welcome email for session %s: %v\n", sessionID, err)
			}
		}()
	}

	// Update metrics
	evo.orchestratorMetrics.ActiveSessions++
	evo.orchestratorMetrics.TotalSessionsStarted++

	// Send initial status
	sendStatus(session.StatusChannel, SessionStatus{
		SessionID: sessionID,
		Status:    "started",
		Message:   "Expert validation session initiated",
		Progress:  0.0,
		Timestamp: time.Now(),
	})

	return session, nil
}

// ProcessExpertFeedback processes real-time feedback from experts
func (evo *ExpertValidationOrchestrator) ProcessExpertFeedback(feedback ExpertFeedback) error {
	sessionID := feedback.SessionID

	evo.sessionMutex.RLock()
	session, exists := evo.activeSessions[sessionID]
	evo.sessionMutex.RUnlock()

	if !exists {
		return errors.NewValidationError("expert_orchestrator", "process_feedback", fmt.Sprintf("session not found: %s", sessionID))
	}

	// Update session activity
	session.LastActivity = time.Now()

	// Process feedback based on current phase
	switch session.CurrentPhase {
	case PhaseBlindAssessment:
		return evo.processBlindAssessmentFeedback(session, feedback)
	case PhaseReferenceComparison:
		return evo.processReferenceComparisonFeedback(session, feedback)
	case PhaseDiscussion:
		return evo.processDiscussionFeedback(session, feedback)
	default:
		return errors.NewValidationError("expert_orchestrator", "process_feedback", fmt.Sprintf("invalid phase for feedback: %v", session.CurrentPhase))
	}
}

// processBlindAssessmentFeedback processes feedback during blind assessment phase
func (evo *ExpertValidationOrchestrator) processBlindAssessmentFeedback(session *LiveValidationSession, feedback ExpertFeedback) error {
	testCaseID := feedback.TestCaseID

	// Extract feedback components
	classification, ok := feedback.Feedback["classification"].(string)
	if !ok {
		return errors.NewValidationError("expert_orchestrator", "process_blind_feedback", "classification missing from feedback")
	}

	confidence, ok := feedback.Feedback["confidence"].(float64)
	if !ok {
		return errors.NewValidationError("expert_orchestrator", "process_blind_feedback", "confidence missing from feedback")
	}

	justification, ok := feedback.Feedback["justification"].(string)
	if !ok {
		justification = "" // Optional field
	}

	// Record feedback with framework
	if err := evo.framework.CollectExpertFeedback(testCaseID, session.ExpertID, classification, confidence, justification); err != nil {
		return errors.WrapError(err, errors.ValidationError, "expert_orchestrator", "process_blind_feedback", "failed to collect expert feedback")
	}

	// Update session progress
	session.CompletedTests = append(session.CompletedTests, testCaseID)
	session.ProgressTracker.CompletedCases++
	session.ProgressTracker.CurrentCaseID = testCaseID

	// Calculate progress
	progress := float64(session.ProgressTracker.CompletedCases) / float64(session.ProgressTracker.TotalTestCases)

	// Send progress update
	sendStatus(session.StatusChannel, SessionStatus{
		SessionID: session.SessionID,
		Status:    "progress",
		Message:   fmt.Sprintf("Completed test case %s (%d/%d)", testCaseID, session.ProgressTracker.CompletedCases, session.ProgressTracker.TotalTestCases),
		Progress:  progress,
		Timestamp: time.Now(),
	})

	// Check if blind assessment phase is complete
	if session.ProgressTracker.CompletedCases >= session.ProgressTracker.TotalTestCases {
		return evo.advanceToNextPhase(session)
	}

	return nil
}

// processReferenceComparisonFeedback processes feedback during reference comparison phase
func (evo *ExpertValidationOrchestrator) processReferenceComparisonFeedback(session *LiveValidationSession, feedback ExpertFeedback) error {
	// Store reference comparison data
	session.SessionData["reference_comparison"] = feedback.Feedback

	// Advance to discussion phase
	return evo.advanceToNextPhase(session)
}

// processDiscussionFeedback processes feedback during discussion phase
func (evo *ExpertValidationOrchestrator) processDiscussionFeedback(session *LiveValidationSession, feedback ExpertFeedback) error {
	// Store discussion feedback
	session.SessionData["discussion_feedback"] = feedback.Feedback

	// Complete the session
	return evo.completeSession(session)
}

// advanceToNextPhase advances the session to the next validation phase
func (evo *ExpertValidationOrchestrator) advanceToNextPhase(session *LiveValidationSession) error {
	switch session.CurrentPhase {
	case PhasePreparation:
		session.CurrentPhase = PhaseBlindAssessment
	case PhaseBlindAssessment:
		session.CurrentPhase = PhaseReferenceComparison
	case PhaseReferenceComparison:
		session.CurrentPhase = PhaseDiscussion
	case PhaseDiscussion:
		session.CurrentPhase = PhaseCompletion
		return evo.completeSession(session)
	default:
		return errors.NewValidationError("expert_orchestrator", "advance_phase", fmt.Sprintf("cannot advance from phase: %v", session.CurrentPhase))
	}

	session.ProgressTracker.PhaseStartTime = time.Now()

	// Send phase change notification
	sendStatus(session.StatusChannel, SessionStatus{
		SessionID: session.SessionID,
		Status:    "phase_change",
		Message:   fmt.Sprintf("Advanced to phase: %v", session.CurrentPhase),
		Progress:  evo.calculatePhaseProgress(session),
		Timestamp: time.Now(),
	})

	return nil
}

// completeSession completes an expert validation session
func (evo *ExpertValidationOrchestrator) completeSession(session *LiveValidationSession) error {
	evo.sessionMutex.Lock()
	defer evo.sessionMutex.Unlock()

	// Mark session as completed
	session.IsActive = false
	session.CurrentPhase = PhaseCompletion

	// Generate individual session report
	sessionReport, err := evo.generateSessionReport(session)
	if err != nil {
		return errors.WrapError(err, errors.ResourceError, "expert_orchestrator", "complete_session", "failed to generate session report")
	}

	// Store session report
	session.SessionData["final_report"] = sessionReport

	// Remove from active sessions
	delete(evo.activeSessions, session.SessionID)

	// Update metrics
	evo.orchestratorMetrics.ActiveSessions--
	evo.orchestratorMetrics.CompletedSessions++

	// Send completion notification
	sendStatus(session.StatusChannel, SessionStatus{
		SessionID: session.SessionID,
		Status:    "completed",
		Message:   "Expert validation session completed successfully",
		Progress:  1.0,
		Timestamp: time.Now(),
	})

	// Send completion email
	if evo.config.EmailNotifications {
		expertInfo := session.ExpertInfo
		sessionID := session.SessionID
		evo.wg.Add(1)
		go func() {
			defer evo.wg.Done()
			if err := evo.sendCompletionEmail(expertInfo, sessionID, sessionReport); err != nil {
				fmt.Printf("warning: failed to send completion email for session %s: %v\n", sessionID, err)
			}
		}()
	}

	// Check if we have enough experts to generate overall validation report
	if evo.orchestratorMetrics.CompletedSessions >= evo.config.MinimumExpertsRequired {
		evo.wg.Add(1)
		go func() {
			defer evo.wg.Done()
			if err := evo.generateOverallValidationReport(); err != nil {
				fmt.Printf("warning: failed to generate overall validation report: %v\n", err)
			}
		}()
	}

	return nil
}

// generateSessionReport generates a report for a completed session
func (evo *ExpertValidationOrchestrator) generateSessionReport(session *LiveValidationSession) (*SessionValidationReport, error) {
	// Calculate session metrics
	sessionDuration := time.Since(session.StartTime)

	// Get expert performance from framework
	expertStats := evo.framework.GetFrameworkStats()
	expertPerformance, exists := expertStats["expert_performance"].(map[string]float64)
	if !exists {
		expertPerformance = make(map[string]float64)
	}

	accuracy := expertPerformance[session.ExpertID]

	report := &SessionValidationReport{
		SessionID:       session.SessionID,
		ExpertID:        session.ExpertID,
		ExpertInfo:      session.ExpertInfo,
		SessionDuration: sessionDuration,
		CompletedCases:  len(session.CompletedTests),
		TotalCases:      session.ProgressTracker.TotalTestCases,
		Accuracy:        accuracy,
		Timestamp:       time.Now(),
		SessionData:     session.SessionData,
	}

	// Add phase-specific analysis
	if blindData, exists := session.SessionData["blind_assessment"]; exists {
		report.BlindAssessmentResults = blindData
	}

	if refData, exists := session.SessionData["reference_comparison"]; exists {
		report.ReferenceComparisonResults = refData
	}

	if discussionData, exists := session.SessionData["discussion_feedback"]; exists {
		report.DiscussionResults = discussionData
	}

	return report, nil
}

// generateOverallValidationReport generates overall validation report from all completed sessions
func (evo *ExpertValidationOrchestrator) generateOverallValidationReport() error {
	// Aggregate results from framework
	if err := evo.framework.AggregateResults(); err != nil {
		return errors.WrapError(err, errors.ValidationError, "expert_orchestrator", "generate_overall_report", "failed to aggregate results")
	}

	// Generate comprehensive report
	overallReport := &OverallValidationReport{
		Timestamp:            time.Now(),
		TotalExperts:         evo.orchestratorMetrics.CompletedSessions,
		FrameworkStats:       evo.framework.GetFrameworkStats(),
		ConsensusResults:     evo.framework.resultAggregator.consensusResults,
		DisagreementAnalysis: evo.framework.resultAggregator.disagreementAnalysis,
		AuthenticityScore:    evo.framework.GetAuthenticityScore(),
		OrchestratorMetrics:  evo.orchestratorMetrics,
	}

	// Add detailed analysis
	overallReport.QualityAssessment = evo.assessOverallQuality()
	overallReport.Recommendations = evo.generateRecommendations()

	// Save report to file system
	if err := evo.saveOverallReport(overallReport); err != nil {
		return errors.WrapError(err, errors.ResourceError, "expert_orchestrator", "generate_overall_report", "failed to save overall report")
	}

	// Generate publication materials if enabled
	if evo.config.AutoReportGeneration {
		evo.wg.Add(1)
		go func() {
			defer evo.wg.Done()
			if err := evo.generatePublicationMaterials(overallReport); err != nil {
				fmt.Printf("warning: failed to generate publication materials: %v\n", err)
			}
		}()
	}

	// Integrate with academic systems if enabled
	if evo.config.AcademicIntegration {
		evo.wg.Add(1)
		go func() {
			defer evo.wg.Done()
			if err := evo.academicIntegrator.SubmitValidationResults(overallReport); err != nil {
				fmt.Printf("warning: failed to submit validation results to academic system: %v\n", err)
			}
		}()
	}

	return nil
}

// calculatePhaseProgress calculates progress within the current phase
func (evo *ExpertValidationOrchestrator) calculatePhaseProgress(session *LiveValidationSession) float64 {
	baseProgress := float64(int(session.CurrentPhase)) / 4.0 // 4 main phases

	switch session.CurrentPhase {
	case PhaseBlindAssessment:
		phaseProgress := float64(session.ProgressTracker.CompletedCases) / float64(session.ProgressTracker.TotalTestCases)
		return baseProgress + (phaseProgress * 0.25) // Blind assessment is 25% of total
	case PhaseReferenceComparison:
		return baseProgress + 0.25 // Constant progress for comparison phase
	case PhaseDiscussion:
		return baseProgress + 0.25 // Constant progress for discussion phase
	default:
		return baseProgress
	}
}

// monitorSessions monitors active sessions for timeouts and issues.
// It runs until the orchestrator's context is cancelled.
func (evo *ExpertValidationOrchestrator) monitorSessions() {
	defer evo.wg.Done()
	ticker := time.NewTicker(1 * time.Minute)
	defer ticker.Stop()

	for {
		select {
		case <-evo.ctx.Done():
			return
		case <-ticker.C:
			evo.sessionMutex.RLock()
			for sessionID, session := range evo.activeSessions {
				// Check for session timeout
				if time.Since(session.LastActivity) > evo.config.SessionTimeout {
					evo.wg.Add(1)
					go func(id string, s *LiveValidationSession) {
						defer evo.wg.Done()
						evo.handleSessionTimeout(id, s)
					}(sessionID, session)
				}

				// Check for maximum session duration
				if time.Since(session.StartTime) > evo.config.MaxSessionDuration {
					evo.wg.Add(1)
					go func(id string, s *LiveValidationSession) {
						defer evo.wg.Done()
						evo.handleSessionExpiry(id, s)
					}(sessionID, session)
				}
			}
			evo.sessionMutex.RUnlock()
		}
	}
}

// Shutdown cancels the orchestrator's context and waits for all background
// goroutines to finish. Call this when the orchestrator is no longer needed.
func (evo *ExpertValidationOrchestrator) Shutdown() {
	evo.cancel()
	evo.wg.Wait()
}

// handleSessionTimeout handles session timeouts
func (evo *ExpertValidationOrchestrator) handleSessionTimeout(sessionID string, session *LiveValidationSession) {
	// Send timeout notification
	sendStatus(session.StatusChannel, SessionStatus{
		SessionID: sessionID,
		Status:    "timeout_warning",
		Message:   "Session approaching timeout due to inactivity",
		Progress:  evo.calculatePhaseProgress(session),
		Timestamp: time.Now(),
	})

	// Send email notification if configured
	if evo.config.EmailNotifications {
		evo.sendTimeoutWarningEmail(session.ExpertInfo, sessionID)
	}
}

// handleSessionExpiry handles session expiry
func (evo *ExpertValidationOrchestrator) handleSessionExpiry(sessionID string, session *LiveValidationSession) {
	evo.sessionMutex.Lock()
	defer evo.sessionMutex.Unlock()

	// Mark session as expired
	session.IsActive = false

	// Generate partial report
	partialReport, _ := evo.generateSessionReport(session)
	session.SessionData["partial_report"] = partialReport
	session.SessionData["expiry_reason"] = "maximum_duration_exceeded"

	// Remove from active sessions
	delete(evo.activeSessions, sessionID)

	// Update metrics
	evo.orchestratorMetrics.ActiveSessions--
	evo.orchestratorMetrics.ExpiredSessions++

	// Send expiry notification
	sendStatus(session.StatusChannel, SessionStatus{
		SessionID: sessionID,
		Status:    "expired",
		Message:   "Session expired due to maximum duration limit",
		Progress:  evo.calculatePhaseProgress(session),
		Timestamp: time.Now(),
	})
}

// sendStatus delivers a status update to a session's StatusChannel.
// The send is non-blocking: if the buffer is full and no consumer is reading,
// the update is dropped rather than stalling the caller goroutine.
func sendStatus(ch chan SessionStatus, status SessionStatus) {
	select {
	case ch <- status:
	default:
	}
}

// GetActiveSessionStats returns statistics about active sessions
func (evo *ExpertValidationOrchestrator) GetActiveSessionStats() map[string]interface{} {
	evo.sessionMutex.RLock()
	defer evo.sessionMutex.RUnlock()

	stats := map[string]interface{}{
		"active_sessions":    len(evo.activeSessions),
		"completed_sessions": evo.orchestratorMetrics.CompletedSessions,
		"expired_sessions":   evo.orchestratorMetrics.ExpiredSessions,
		"total_started":      evo.orchestratorMetrics.TotalSessionsStarted,
	}

	// Add per-phase breakdown
	phaseBreakdown := make(map[string]int)
	for _, session := range evo.activeSessions {
		phaseBreakdown[session.CurrentPhase.String()]++
	}
	stats["phase_breakdown"] = phaseBreakdown

	return stats
}

// String returns string representation of validation phase
func (vp ValidationPhase) String() string {
	switch vp {
	case PhasePreparation:
		return "preparation"
	case PhaseBlindAssessment:
		return "blind_assessment"
	case PhaseReferenceComparison:
		return "reference_comparison"
	case PhaseDiscussion:
		return "discussion"
	case PhaseCompletion:
		return "completion"
	default:
		return "unknown"
	}
}

// Supporting types and structures

type SessionManager struct {
	sessions map[string]*LiveValidationSession
	mutex    sync.RWMutex
}

func NewSessionManager() *SessionManager {
	return &SessionManager{
		sessions: make(map[string]*LiveValidationSession),
	}
}

type ValidationMaterialManager struct {
	basePath     string
	testCases    []TestCase
	sequenceData map[string][]byte
}

func NewValidationMaterialManager(basePath string) *ValidationMaterialManager {
	return &ValidationMaterialManager{
		basePath:     basePath,
		testCases:    make([]TestCase, 0),
		sequenceData: make(map[string][]byte),
	}
}

func (vmm *ValidationMaterialManager) LoadValidationMaterials() error {
	// Load blind test sequences
	blindPath := filepath.Join(vmm.basePath, "blind_test_sequences")
	if err := vmm.loadSequencesFromDirectory(blindPath, "blind"); err != nil {
		return fmt.Errorf("failed to load blind sequences: %w", err)
	}

	// Load reference sequences
	refPath := filepath.Join(vmm.basePath, "reference_sequences")
	if err := vmm.loadSequencesFromDirectory(refPath, "reference"); err != nil {
		return fmt.Errorf("failed to load reference sequences: %w", err)
	}

	return nil
}

func (vmm *ValidationMaterialManager) loadSequencesFromDirectory(dirPath, seqType string) error {
	files, err := os.ReadDir(dirPath)
	if err != nil {
		return fmt.Errorf("failed to read directory %s: %w", dirPath, err)
	}

	for _, file := range files {
		if file.IsDir() || filepath.Ext(file.Name()) != ".fasta" {
			continue
		}

		filePath := filepath.Join(dirPath, file.Name())
		content, err := os.ReadFile(filePath) // #nosec G304
		if err != nil {
			return fmt.Errorf("failed to read file %s: %w", filePath, err)
		}

		vmm.sequenceData[seqType+"/"+file.Name()] = content
	}

	return nil
}

func (vmm *ValidationMaterialManager) GetTestCases() ([]TestCase, error) {
	return vmm.testCases, nil
}

func (vmm *ValidationMaterialManager) LoadSequenceData() ([][]byte, [][]byte, error) {
	var synthetic, real [][]byte

	for filename, content := range vmm.sequenceData {
		// Files are keyed as "seqType/filename", where seqType is "blind" or "reference"
		if strings.HasPrefix(filename, "blind/") {
			synthetic = append(synthetic, content)
		} else {
			real = append(real, content)
		}
	}

	return synthetic, real, nil
}

type TestCase struct {
	ID           string
	FilePath     string
	SequenceType string
	Population   string
}

type ValidationReportGenerator struct{}

func NewValidationReportGenerator() *ValidationReportGenerator {
	return &ValidationReportGenerator{}
}

type AcademicIntegrator struct{}

func NewAcademicIntegrator() *AcademicIntegrator {
	return &AcademicIntegrator{}
}

func (ai *AcademicIntegrator) SubmitValidationResults(report *OverallValidationReport) error {
	// Stub: nothing is submitted.
	return nil
}

type OrchestratorMetrics struct {
	TotalSessionsStarted int
	ActiveSessions       int
	CompletedSessions    int
	ExpiredSessions      int

	// Performance metrics
	AverageSessionDuration time.Duration
	AverageAccuracy        float64
	ExpertParticipation    map[string]int
}

func NewOrchestratorMetrics() *OrchestratorMetrics {
	return &OrchestratorMetrics{
		ExpertParticipation: make(map[string]int),
	}
}

func NewSessionProgressTracker() *SessionProgressTracker {
	return &SessionProgressTracker{
		PhaseStartTime: time.Now(),
	}
}

// Report structures

type SessionValidationReport struct {
	SessionID       string
	ExpertID        string
	ExpertInfo      *Expert
	SessionDuration time.Duration
	CompletedCases  int
	TotalCases      int
	Accuracy        float64
	Timestamp       time.Time
	SessionData     map[string]interface{}

	// Phase-specific results
	BlindAssessmentResults     interface{}
	ReferenceComparisonResults interface{}
	DiscussionResults          interface{}
}

type OverallValidationReport struct {
	Timestamp            time.Time
	TotalExperts         int
	FrameworkStats       map[string]interface{}
	ConsensusResults     []ConsensusResult
	DisagreementAnalysis []DisagreementAnalysis
	AuthenticityScore    float64
	OrchestratorMetrics  *OrchestratorMetrics

	// Analysis results
	QualityAssessment QualityAssessment
	Recommendations   []string
}

type QualityAssessment struct {
	OverallQuality       string
	StrengthAreas        []string
	ImprovementAreas     []string
	PublicationReadiness bool
	ClinicalUtility      string
}

// Email and notification methods

// sendEmail sends a plain-text email via the configured SMTP server.
// Returns nil (no-op) when EmailNotifications is disabled or no recipient address
// is set. The SMTPConfig.Password is never logged or included in error messages.
func (evo *ExpertValidationOrchestrator) sendEmail(to, subject, body string) error {
	cfg := evo.config.SMTPConfig
	if !evo.config.EmailNotifications || to == "" {
		return nil
	}
	if cfg.Host == "" {
		return fmt.Errorf("sendEmail: SMTPConfig.Host not configured")
	}
	addr := fmt.Sprintf("%s:%d", cfg.Host, cfg.Port)
	from := cfg.FromAddr
	if from == "" {
		from = cfg.Username
	}
	msg := fmt.Sprintf("From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n%s", from, to, subject, body)
	auth := smtp.PlainAuth("", cfg.Username, cfg.Password, cfg.Host)
	if err := smtp.SendMail(addr, auth, from, []string{to}, []byte(msg)); err != nil {
		return fmt.Errorf("sendEmail to %s: %w", to, err)
	}
	return nil
}

func (evo *ExpertValidationOrchestrator) sendWelcomeEmail(expert *Expert, sessionID string) error {
	subject := fmt.Sprintf("[PhantomScience] Expert Validation Session %s — Welcome", sessionID)
	body := fmt.Sprintf(
		"Dear %s,\n\nYou have been invited to participate in an expert validation session.\n\n"+
			"Session ID: %s\n"+
			"Institution: %s\n\n"+
			"Please quote the session ID above with your blind assessment.\n\n"+
			"Thank you for your contribution.\n\nPhantomScience Validation Team\n",
		expert.Name, sessionID, expert.Institution,
	)
	return evo.sendEmail(expert.Email, subject, body)
}

func (evo *ExpertValidationOrchestrator) sendCompletionEmail(expert *Expert, sessionID string, report *SessionValidationReport) error {
	subject := fmt.Sprintf("[PhantomScience] Session %s — Completed", sessionID)
	body := fmt.Sprintf(
		"Dear %s,\n\nYour validation session has been completed.\n\n"+
			"Session ID:      %s\n"+
			"Duration:        %s\n"+
			"Cases completed: %d / %d\n"+
			"Accuracy:        %.1f%%\n\n"+
			"Thank you for your participation.\n\nPhantomScience Validation Team\n",
		expert.Name, sessionID,
		report.SessionDuration.Round(time.Second),
		report.CompletedCases, report.TotalCases,
		report.Accuracy*100,
	)
	return evo.sendEmail(expert.Email, subject, body)
}

func (evo *ExpertValidationOrchestrator) sendTimeoutWarningEmail(expert *Expert, sessionID string) error {
	subject := fmt.Sprintf("[PhantomScience] Session %s — Inactivity Warning", sessionID)
	body := fmt.Sprintf(
		"Dear %s,\n\nYour validation session %s is approaching the inactivity timeout.\n\n"+
			"Please continue your assessment before the session expires.\n\n"+
			"PhantomScience Validation Team\n",
		expert.Name, sessionID,
	)
	return evo.sendEmail(expert.Email, subject, body)
}

// Analysis methods

// assessOverallQuality grades the expert authenticity score with fixed
// cut-offs. PublicationReadiness is only "score >= 60"; strength and
// improvement areas and clinical utility are not assessed.
func (evo *ExpertValidationOrchestrator) assessOverallQuality() QualityAssessment {
	authenticityScore := evo.framework.GetAuthenticityScore()
	_ = evo.framework.GetFrameworkStats()

	var quality string
	var publicationReady bool

	if authenticityScore >= 75.0 {
		quality = "excellent"
		publicationReady = true
	} else if authenticityScore >= 60.0 {
		quality = "good"
		publicationReady = true
	} else if authenticityScore >= 45.0 {
		quality = "acceptable"
		publicationReady = false
	} else {
		quality = "needs_improvement"
		publicationReady = false
	}

	return QualityAssessment{
		OverallQuality:       quality,
		StrengthAreas:        []string{"not assessed"},
		ImprovementAreas:     []string{"not assessed"},
		PublicationReadiness: publicationReady,
		ClinicalUtility:      "not assessed",
	}
}

// generateRecommendations returns a fixed list; it does not depend on the
// results.
func (evo *ExpertValidationOrchestrator) generateRecommendations() []string {
	recommendations := []string{
		"Continue validation with additional expert panels",
		"Optimize HMM parameters for improved authenticity scores",
	}

	return recommendations
}

func (evo *ExpertValidationOrchestrator) generatePublicationMaterials(report *OverallValidationReport) error {
	// Stub: no publication material is generated.
	return nil
}

func (evo *ExpertValidationOrchestrator) saveOverallReport(report *OverallValidationReport) error {
	// Create reports directory if it doesn't exist
	reportsDir := filepath.Join(evo.config.ExpertPackagePath, "reports")
	if err := security.SecureMkdirAll(reportsDir, evo.config.ExpertPackagePath); err != nil {
		return fmt.Errorf("failed to create reports directory: %w", err)
	}

	// Generate filename with timestamp
	filename := fmt.Sprintf("overall_validation_report_%d.json", time.Now().Unix())
	filePath := filepath.Join(reportsDir, filename)

	// Marshal report to JSON
	reportData, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		return fmt.Errorf("failed to marshal report: %w", err)
	}

	// Write to file
	if err := os.WriteFile(filePath, reportData, 0600); err != nil {
		return fmt.Errorf("failed to write report file: %w", err)
	}

	return nil
}
