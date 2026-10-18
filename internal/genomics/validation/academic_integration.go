package validation

import (
	"context"
	"encoding/json"
	"fmt"
	"io"
	"math"
	"net/http"
	"net/url"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/errors"
)

// AcademicIntegrationAPI ranks expert candidates and drafts invitation and
// manuscript text. The PubMed, ORCID, CrossRef and institution clients make
// real HTTP requests. Email sending, journal formatting and submission, and
// the collaboration platform are stubs (nothing is sent or uploaded), and the
// manuscript sections are fixed template text.
type AcademicIntegrationAPI struct {
	config         AcademicConfig
	httpClient     *http.Client
	pubmedClient   *PubMedClient
	orcidClient    *ORCIDClient
	crossrefClient *CrossRefClient

	// Institution connections
	institutionAPIs map[string]*InstitutionAPI

	// Manuscript management
	manuscriptManager *ManuscriptManager

	// Collaboration platform
	collaborationPlatform *CollaborationPlatform
}

// AcademicConfig contains configuration for academic system integration
type AcademicConfig struct {
	// API Keys and credentials
	PubMedAPIKey   string `json:"-" yaml:"-"`
	ORCIDAPIKey    string `json:"-" yaml:"-"`
	CrossRefAPIKey string `json:"-" yaml:"-"`

	// Institution-specific configurations
	InstitutionConfigs map[string]InstitutionConfig

	// Manuscript submission settings
	ManuscriptConfig ManuscriptConfig

	// Collaboration settings
	CollaborationConfig CollaborationConfig

	// Research networking
	ResearchNetworking ResearchNetworkingConfig
}

// InstitutionConfig contains configuration for specific academic institutions
type InstitutionConfig struct {
	Name         string
	APIEndpoint  string
	APIKey       string `json:"-" yaml:"-"` // #nosec G117
	ContactEmail string

	// Department-specific contacts
	Departments map[string]DepartmentContact

	// Integration capabilities
	SupportsDirectSubmission bool
	SupportsCollaboration    bool
	SupportsDataSharing      bool
}

// DepartmentContact contains contact information for academic departments
type DepartmentContact struct {
	Name           string
	Email          string
	Phone          string
	ResearchAreas  []string
	ExpertiseLevel string
}

// ManuscriptConfig contains configuration for manuscript management
type ManuscriptConfig struct {
	TargetJournals     []JournalConfig
	AutoGeneration     bool
	TemplateLibrary    string
	ReviewWorkflow     bool
	SubmissionTracking bool
}

// JournalConfig contains configuration for target journals
type JournalConfig struct {
	Name            string
	ImpactFactor    float64
	SubmissionAPI   string
	APIKey          string `json:"-" yaml:"-"` // #nosec G117
	PreferredFormat string
	ReviewTimeline  time.Duration

	// Submission requirements
	Requirements SubmissionRequirements
}

// SubmissionRequirements contains journal submission requirements
type SubmissionRequirements struct {
	WordLimit           int
	FigureLimit         int
	SupplementaryLimit  int
	EthicsApproval      bool
	DataAvailability    bool
	ConflictDeclaration bool
}

// CollaborationConfig contains configuration for research collaboration
type CollaborationConfig struct {
	PlatformType      string // "slack", "teams", "discord", "custom"
	APIEndpoint       string
	Credentials       map[string]string
	AutoInviteExperts bool
	SharedWorkspace   bool
	RealTimeEditing   bool
}

// ResearchNetworkingConfig contains configuration for research networking
type ResearchNetworkingConfig struct {
	LinkedInIntegration bool
	ResearchGateAPI     string
	AcademiaEduAPI      string
	GoogleScholarAPI    string

	// Automatic profile updates
	AutoProfileUpdate     bool
	PublicationSync       bool
	CollaborationTracking bool
}

// NewAcademicIntegrationAPI creates a new academic integration API
func NewAcademicIntegrationAPI(config AcademicConfig) (*AcademicIntegrationAPI, error) {
	// Validate configuration
	if err := config.Validate(); err != nil {
		return nil, errors.WrapError(err, errors.ConfigurationError, "academic_integration", "new_api", "invalid configuration")
	}

	// Create HTTP client with timeout
	httpClient := &http.Client{
		Timeout: 30 * time.Second,
	}

	api := &AcademicIntegrationAPI{
		config:                config,
		httpClient:            httpClient,
		pubmedClient:          NewPubMedClient(config.PubMedAPIKey, httpClient),
		orcidClient:           NewORCIDClient(config.ORCIDAPIKey, httpClient),
		crossrefClient:        NewCrossRefClient(config.CrossRefAPIKey, httpClient),
		institutionAPIs:       make(map[string]*InstitutionAPI),
		manuscriptManager:     NewManuscriptManager(config.ManuscriptConfig),
		collaborationPlatform: NewCollaborationPlatform(config.CollaborationConfig),
	}

	// Initialize institution APIs
	for name, instConfig := range config.InstitutionConfigs {
		institutionAPI, err := NewInstitutionAPI(instConfig, httpClient)
		if err != nil {
			return nil, fmt.Errorf("failed to initialize institution API for %s: %w", name, err)
		}
		api.institutionAPIs[name] = institutionAPI
	}

	return api, nil
}

// Validate validates the academic configuration
func (config *AcademicConfig) Validate() error {
	if len(config.InstitutionConfigs) == 0 {
		return fmt.Errorf("at least one institution configuration is required")
	}

	// Validate journal configurations
	for _, journal := range config.ManuscriptConfig.TargetJournals {
		if journal.Name == "" {
			return fmt.Errorf("journal name is required")
		}
	}

	return nil
}

// RecruitExperts searches the configured institution APIs, ranks the
// candidates, and starts a background invitation of the top 10. Invitations
// are not delivered yet because sendEmail is a stub.
func (api *AcademicIntegrationAPI) RecruitExperts(ctx context.Context, criteria ExpertCriteria) ([]*ExpertCandidate, error) {
	var allCandidates []*ExpertCandidate

	// Search across all configured institutions
	for institutionName, institutionAPI := range api.institutionAPIs {
		candidates, err := institutionAPI.SearchExperts(ctx, criteria)
		if err != nil {
			// Log error but continue with other institutions
			fmt.Printf("Failed to search experts at %s: %v\n", institutionName, err)
			continue
		}

		// Tag candidates with institution
		for _, candidate := range candidates {
			candidate.Institution = institutionName
		}

		allCandidates = append(allCandidates, candidates...)
	}

	// Rank candidates based on expertise match
	rankedCandidates := api.rankExpertCandidates(allCandidates, criteria)

	// Auto-invite top candidates if configured
	if len(rankedCandidates) > 0 {
		go api.autoInviteExperts(ctx, rankedCandidates[:min(10, len(rankedCandidates))], criteria)
	}

	return rankedCandidates, nil
}

// ExpertCriteria defines criteria for expert recruitment
type ExpertCriteria struct {
	ResearchAreas   []string
	ExpertiseLevel  string   // "junior", "senior", "expert"
	Institution     []string // Preferred institutions
	Publications    int      // Minimum publication count
	ExperienceYears int      // Minimum years of experience
	Specializations []string // Specific genomics specializations

	// Geographic preferences
	TimeZone string
	Region   string

	// Collaboration preferences
	RemoteAvailable bool
	Languages       []string
}

// ExpertCandidate represents a potential expert for validation
type ExpertCandidate struct {
	Name           string
	Email          string
	Institution    string
	Department     string
	ResearchAreas  []string
	Publications   int
	HIndex         int
	ExpertiseScore float64

	// Contact information
	Phone           string
	ORCID           string
	LinkedInProfile string

	// Availability
	AvailableFrom  time.Time
	AvailableTo    time.Time
	PreferredTimes []TimeSlot

	// Expertise matching
	MatchScore   float64
	MatchDetails map[string]float64
}

// TimeSlot represents an available time slot
type TimeSlot struct {
	StartTime time.Time
	EndTime   time.Time
	TimeZone  string
}

// rankExpertCandidates ranks expert candidates based on criteria match
func (api *AcademicIntegrationAPI) rankExpertCandidates(candidates []*ExpertCandidate, criteria ExpertCriteria) []*ExpertCandidate {
	// Calculate match scores for each candidate
	for _, candidate := range candidates {
		candidate.MatchScore = api.calculateExpertMatchScore(candidate, criteria)
	}

	// Sort by match score (descending)
	for i := 0; i < len(candidates)-1; i++ {
		for j := i + 1; j < len(candidates); j++ {
			if candidates[i].MatchScore < candidates[j].MatchScore {
				candidates[i], candidates[j] = candidates[j], candidates[i]
			}
		}
	}

	return candidates
}

// calculateExpertMatchScore calculates how well a candidate matches the criteria
func (api *AcademicIntegrationAPI) calculateExpertMatchScore(candidate *ExpertCandidate, criteria ExpertCriteria) float64 {
	var totalScore float64
	var weightSum float64

	// Research area match (40% weight)
	areaScore := api.calculateResearchAreaMatch(candidate.ResearchAreas, criteria.ResearchAreas)
	totalScore += areaScore * 0.4
	weightSum += 0.4

	// Publication count (20% weight)
	pubScore := math.Min(float64(candidate.Publications)/math.Max(float64(criteria.Publications), 1.0), 1.0)
	totalScore += pubScore * 0.2
	weightSum += 0.2

	// H-index (15% weight)
	hIndexScore := math.Min(float64(candidate.HIndex)/20.0, 1.0) // Normalize to 20
	totalScore += hIndexScore * 0.15
	weightSum += 0.15

	// Institution preference (10% weight)
	instScore := 0.0
	if len(criteria.Institution) == 0 {
		instScore = 1.0 // No preference
	} else {
		for _, prefInst := range criteria.Institution {
			if candidate.Institution == prefInst {
				instScore = 1.0
				break
			}
		}
	}
	totalScore += instScore * 0.1
	weightSum += 0.1

	// Expertise level (15% weight)
	levelScore := api.calculateExpertiseLevelMatch(candidate, criteria.ExpertiseLevel)
	totalScore += levelScore * 0.15
	weightSum += 0.15

	if weightSum > 0 {
		return totalScore / weightSum
	}

	return 0.0
}

// calculateResearchAreaMatch calculates research area overlap
func (api *AcademicIntegrationAPI) calculateResearchAreaMatch(candidateAreas, criteriaAreas []string) float64 {
	if len(criteriaAreas) == 0 {
		return 1.0 // No specific requirements
	}

	matches := 0
	for _, criteriaArea := range criteriaAreas {
		for _, candidateArea := range candidateAreas {
			if api.isResearchAreaMatch(candidateArea, criteriaArea) {
				matches++
				break
			}
		}
	}

	return float64(matches) / float64(len(criteriaAreas))
}

// isResearchAreaMatch reports whether two research areas are the same string
// (exact, case-sensitive; no fuzzy matching).
func (api *AcademicIntegrationAPI) isResearchAreaMatch(candidateArea, criteriaArea string) bool {
	return candidateArea == criteriaArea
}

// calculateExpertiseLevelMatch calculates expertise level match
func (api *AcademicIntegrationAPI) calculateExpertiseLevelMatch(candidate *ExpertCandidate, requiredLevel string) float64 {
	// Map expertise levels to numeric values
	levelValues := map[string]int{
		"junior": 1,
		"senior": 2,
		"expert": 3,
	}

	// Infer candidate level from publications and h-index
	candidateLevel := "junior"
	if candidate.Publications > 50 && candidate.HIndex > 15 {
		candidateLevel = "expert"
	} else if candidate.Publications > 20 && candidate.HIndex > 8 {
		candidateLevel = "senior"
	}

	candidateLevelValue := levelValues[candidateLevel]
	requiredLevelValue := levelValues[requiredLevel]

	if candidateLevelValue >= requiredLevelValue {
		return 1.0
	}

	return float64(candidateLevelValue) / float64(requiredLevelValue)
}

// autoInviteExperts builds and "sends" invitations to the given candidates;
// see sendEmail, which is a stub.
func (api *AcademicIntegrationAPI) autoInviteExperts(ctx context.Context, candidates []*ExpertCandidate, criteria ExpertCriteria) error {
	for _, candidate := range candidates {
		// Generate personalized invitation
		invitation := api.generatePersonalizedInvitation(candidate, criteria)

		// Send invitation email
		if err := api.sendExpertInvitation(ctx, candidate, invitation); err != nil {
			fmt.Printf("Failed to send invitation to %s: %v\n", candidate.Email, err)
			continue
		}

		// Track invitation in collaboration platform
		api.collaborationPlatform.TrackInvitation(candidate, invitation)

		// Add delay between invitations to avoid spam detection
		time.Sleep(2 * time.Second)
	}

	return nil
}

// generatePersonalizedInvitation generates a personalized invitation for an expert
func (api *AcademicIntegrationAPI) generatePersonalizedInvitation(candidate *ExpertCandidate, criteria ExpertCriteria) *ExpertInvitation {
	invitation := &ExpertInvitation{
		RecipientName:       candidate.Name,
		RecipientEmail:      candidate.Email,
		Institution:         candidate.Institution,
		PersonalizedMessage: api.generatePersonalizedMessage(candidate, criteria),
		ProjectDetails:      api.generateProjectDetails(criteria),
		TimeCommitment:      "4-6 hours over 1-2 sessions",
		Compensation:        "Academic collaboration and co-authorship opportunity",
		Deadline:            time.Now().Add(14 * 24 * time.Hour),
		ContactInfo:         api.getContactInformation(),
	}

	return invitation
}

// generatePersonalizedMessage creates a personalized message for the expert
func (api *AcademicIntegrationAPI) generatePersonalizedMessage(candidate *ExpertCandidate, criteria ExpertCriteria) string {
	message := fmt.Sprintf(
		"Dear Dr. %s,\n\n"+
			"I am writing to invite you to take part in a blind expert assessment of synthetic genomic sequences. "+
			"Your expertise in %s would help us evaluate them.\n\n"+
			"The sequences are generated by k-mer hidden Markov models trained on 1000 Genomes Project data. "+
			"The assessment asks whether experts can tell them apart from real sequences, and where they differ.",
		candidate.Name,
		formatResearchAreas(candidate.ResearchAreas),
	)

	return message
}

// generateProjectDetails creates detailed project information
func (api *AcademicIntegrationAPI) generateProjectDetails(criteria ExpertCriteria) *ProjectDetails {
	return &ProjectDetails{
		Title:       "Expert Validation of Population-Specific Genomic Sequence Generation",
		Description: "Blind expert assessment of synthetic genomic sequences generated by k-mer HMM models trained on 1000 Genomes Project data",
		Methodology: "Three-phase validation: blind assessment, reference comparison, and expert discussion",
		DataSize:    "Set by the validation materials",
		Innovation:  "Population-specific synthetic sequence generation (under evaluation)",
		Applications: []string{
			"Test data for genomics pipelines",
			"Population genomics method development",
		},
		ExpectedOutcomes: []string{
			"A written report of the assessment results",
			"A list of differences between synthetic and real sequences",
			"Model changes based on expert feedback",
		},
	}
}

// sendExpertInvitation sends an invitation email to an expert
func (api *AcademicIntegrationAPI) sendExpertInvitation(ctx context.Context, candidate *ExpertCandidate, invitation *ExpertInvitation) error {
	// Create email payload
	emailPayload := EmailPayload{
		To:      candidate.Email,
		Subject: "Invitation: blind expert assessment of synthetic genomic sequences",
		Body:    api.formatInvitationEmail(invitation),
		Type:    "html",
	}

	// Send via configured email service
	return api.sendEmail(ctx, emailPayload)
}

// GenerateManuscript assembles a manuscript draft from fixed template text.
// Only the expert count and authenticity score come from validationReport;
// the figures, tables and supplementary paths are placeholders.
func (api *AcademicIntegrationAPI) GenerateManuscript(ctx context.Context, validationReport *OverallValidationReport) (*ManuscriptDraft, error) {
	manuscript := &ManuscriptDraft{
		Title:    api.generateManuscriptTitle(validationReport),
		Abstract: api.generateAbstract(validationReport),
		Authors:  api.generateAuthorList(validationReport),

		// Main sections
		Introduction: api.generateIntroduction(validationReport),
		Methods:      api.generateMethods(validationReport),
		Results:      api.generateResults(validationReport),
		Discussion:   api.generateDiscussion(validationReport),
		Conclusion:   api.generateConclusion(validationReport),

		// Supporting materials
		Figures:       api.generateFigures(validationReport),
		Tables:        api.generateTables(validationReport),
		References:    api.generateReferences(validationReport),
		Supplementary: api.generateSupplementaryMaterials(validationReport),

		// Metadata
		Keywords:        api.generateKeywords(validationReport),
		JournalTargets:  api.selectTargetJournals(validationReport),
		JournalVersions: make(map[string]*FormattedManuscript),
		GeneratedAt:     time.Now(),
	}

	// Auto-format for target journals
	for _, journal := range manuscript.JournalTargets {
		formattedVersion, err := api.manuscriptManager.FormatForJournal(manuscript, journal)
		if err != nil {
			return nil, fmt.Errorf("failed to format manuscript for %s: %w", journal.Name, err)
		}
		manuscript.JournalVersions[journal.Name] = formattedVersion
	}

	return manuscript, nil
}

// SubmitToJournals is a stub: nothing is sent to a journal. Each target gets
// a result with status "submitted" and a made-up tracking URL.
func (api *AcademicIntegrationAPI) SubmitToJournals(ctx context.Context, manuscript *ManuscriptDraft) ([]*SubmissionResult, error) {
	var results []*SubmissionResult

	for _, journal := range manuscript.JournalTargets {
		// Get journal-specific version
		version, exists := manuscript.JournalVersions[journal.Name]
		if !exists {
			continue
		}

		// Submit to journal
		result, err := api.submitToJournal(ctx, version, journal)
		if err != nil {
			result = &SubmissionResult{
				Journal:     journal.Name,
				Status:      "failed",
				Error:       err.Error(),
				SubmittedAt: time.Now(),
			}
		}

		results = append(results, result)

		// Add delay between submissions
		time.Sleep(5 * time.Second)
	}

	return results, nil
}

// TrackSubmissionStatus tracks the status of journal submissions
func (api *AcademicIntegrationAPI) TrackSubmissionStatus(ctx context.Context, submissionID string) (*SubmissionStatus, error) {
	// Check submission status across all journals
	for _, institutionAPI := range api.institutionAPIs {
		status, err := institutionAPI.CheckSubmissionStatus(ctx, submissionID)
		if err == nil && status != nil {
			return status, nil
		}
	}

	return nil, errors.NewResourceError("academic_integration", "track_submission", "submission not found")
}

// SetupCollaborationWorkspace sets up a shared workspace for expert collaboration
func (api *AcademicIntegrationAPI) SetupCollaborationWorkspace(ctx context.Context, experts []*ExpertCandidate, projectDetails *ProjectDetails) (*CollaborationWorkspace, error) {
	workspace := &CollaborationWorkspace{
		Name:        fmt.Sprintf("Expert Validation - %s", projectDetails.Title),
		Description: projectDetails.Description,
		CreatedAt:   time.Now(),
		Members:     make([]*WorkspaceMember, 0),
		Channels:    make([]*Channel, 0),
		SharedFiles: make([]*SharedFile, 0),
	}

	// Create channels for different aspects
	channels := []string{
		"general-discussion",
		"blind-assessment",
		"reference-comparison",
		"results-analysis",
		"manuscript-collaboration",
	}

	for _, channelName := range channels {
		channel, err := api.collaborationPlatform.CreateChannel(ctx, workspace.ID, channelName)
		if err != nil {
			return nil, fmt.Errorf("failed to create channel %s: %w", channelName, err)
		}
		workspace.Channels = append(workspace.Channels, channel)
	}

	// Invite experts to workspace
	for _, expert := range experts {
		member := &WorkspaceMember{
			Name:     expert.Name,
			Email:    expert.Email,
			Role:     "expert_validator",
			JoinedAt: time.Now(),
		}

		if err := api.collaborationPlatform.InviteMember(ctx, workspace.ID, member); err != nil {
			fmt.Printf("Failed to invite %s to workspace: %v\n", expert.Email, err)
			continue
		}

		workspace.Members = append(workspace.Members, member)
	}

	// Upload validation materials
	validationFiles := []string{
		"blind_test_sequences",
		"reference_sequences",
		"evaluation_materials",
	}

	for _, fileName := range validationFiles {
		sharedFile, err := api.collaborationPlatform.UploadFile(ctx, workspace.ID, fileName)
		if err != nil {
			fmt.Printf("Failed to upload %s: %v\n", fileName, err)
			continue
		}
		workspace.SharedFiles = append(workspace.SharedFiles, sharedFile)
	}

	return workspace, nil
}

// Helper functions and utility methods

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func formatResearchAreas(areas []string) string {
	if len(areas) == 0 {
		return "genomics research"
	}
	if len(areas) == 1 {
		return areas[0]
	}
	return fmt.Sprintf("%s, and %s", areas[0], areas[len(areas)-1])
}

// Supporting type definitions

type ExpertInvitation struct {
	RecipientName       string
	RecipientEmail      string
	Institution         string
	PersonalizedMessage string
	ProjectDetails      *ProjectDetails
	TimeCommitment      string
	Compensation        string
	Deadline            time.Time
	ContactInfo         *ContactInfo
}

type ProjectDetails struct {
	Title            string
	Description      string
	Methodology      string
	DataSize         string
	Innovation       string
	Applications     []string
	ExpectedOutcomes []string
}

type ContactInfo struct {
	Name        string
	Email       string
	Phone       string
	Institution string
}

type EmailPayload struct {
	To      string
	Subject string
	Body    string
	Type    string
}

type ManuscriptDraft struct {
	Title           string
	Abstract        string
	Authors         []*Author
	Introduction    string
	Methods         string
	Results         string
	Discussion      string
	Conclusion      string
	Figures         []*Figure
	Tables          []*Table
	References      []*Reference
	Supplementary   []*SupplementaryMaterial
	Keywords        []string
	JournalTargets  []*JournalConfig
	JournalVersions map[string]*FormattedManuscript
	GeneratedAt     time.Time
}

type Author struct {
	Name         string
	Email        string
	Institution  string
	ORCID        string
	Contribution []string
}

type Figure struct {
	ID       string
	Title    string
	Caption  string
	FilePath string
	Type     string
}

type Table struct {
	ID      string
	Title   string
	Caption string
	Headers []string
	Rows    [][]string
}

type Reference struct {
	ID      string
	Authors []string
	Title   string
	Journal string
	Year    int
	DOI     string
}

type SupplementaryMaterial struct {
	ID          string
	Title       string
	Description string
	FilePath    string
	Type        string
}

type FormattedManuscript struct {
	Content     string
	Format      string
	WordCount   int
	FigureCount int
	TableCount  int
}

type SubmissionResult struct {
	Journal         string
	Status          string
	SubmissionID    string
	TrackingURL     string
	Error           string
	SubmittedAt     time.Time
	EstimatedReview time.Duration
}

type SubmissionStatus struct {
	SubmissionID        string
	Journal             string
	CurrentStatus       string
	LastUpdated         time.Time
	ReviewerCount       int
	EstimatedCompletion time.Time
	Messages            []StatusMessage
}

type StatusMessage struct {
	Date    time.Time
	Status  string
	Message string
	Actor   string
}

type CollaborationWorkspace struct {
	ID          string
	Name        string
	Description string
	CreatedAt   time.Time
	Members     []*WorkspaceMember
	Channels    []*Channel
	SharedFiles []*SharedFile
}

type WorkspaceMember struct {
	Name     string
	Email    string
	Role     string
	JoinedAt time.Time
}

type Channel struct {
	ID          string
	Name        string
	Description string
	CreatedAt   time.Time
}

type SharedFile struct {
	ID          string
	Name        string
	Description string
	FilePath    string
	UploadedAt  time.Time
}

// Article represents a PubMed article returned by NCBI E-utilities.
type Article struct {
	PMID    string
	Title   string
	Authors []string
	Journal string
	Year    int
	DOI     string
	URL     string
}

// Researcher represents an ORCID researcher profile.
type Researcher struct {
	ORCID        string
	Name         string
	Email        string
	Institutions []string
	Publications int
	Biography    string
}

// Publication represents a CrossRef publication record.
type Publication struct {
	DOI       string
	Title     string
	Authors   []string
	Journal   string
	Year      int
	Publisher string
	URL       string
	Type      string
	Score     float64
}

// --- private API response types ---

type ncbiESearchResponse struct {
	ESearchResult struct {
		IDList []string `json:"idlist"`
	} `json:"esearchresult"`
}

type ncbiESummaryResponse struct {
	Result map[string]ncbiDocSummary `json:"result"`
}

type ncbiDocSummary struct {
	UID     string `json:"uid"`
	Title   string `json:"title"`
	Source  string `json:"source"`
	PubDate string `json:"pubdate"`
	DOI     string `json:"elocationid"`
	Authors []struct {
		Name string `json:"name"`
	} `json:"authors"`
}

type orcidRecord struct {
	Person struct {
		Name struct {
			GivenNames struct {
				Value string `json:"value"`
			} `json:"given-names"`
			FamilyName struct {
				Value string `json:"value"`
			} `json:"family-name"`
		} `json:"name"`
		Biography struct {
			Content string `json:"content"`
		} `json:"biography"`
		Emails struct {
			Email []struct {
				Email   string `json:"email"`
				Primary bool   `json:"primary"`
			} `json:"email"`
		} `json:"emails"`
	} `json:"person"`
	ActivitiesSummary struct {
		Employments struct {
			AffiliationGroup []struct {
				Summaries []struct {
					EmploymentSummary struct {
						Organization struct {
							Name string `json:"name"`
						} `json:"organization"`
					} `json:"employment-summary"`
				} `json:"summaries"`
			} `json:"affiliation-group"`
		} `json:"employments"`
		Works struct {
			Group []json.RawMessage `json:"group"`
		} `json:"works"`
	} `json:"activities-summary"`
}

type crossRefResponse struct {
	Message crossRefWork `json:"message"`
}

type crossRefWork struct {
	DOI            string   `json:"DOI"`
	Title          []string `json:"title"`
	Score          float64  `json:"score"`
	Type           string   `json:"type"`
	Publisher      string   `json:"publisher"`
	URL            string   `json:"URL"`
	ContainerTitle []string `json:"container-title"`
	Created        struct {
		DateParts [][]int `json:"date-parts"`
	} `json:"created"`
	Author []struct {
		Given  string `json:"given"`
		Family string `json:"family"`
	} `json:"author"`
}

// Client implementations for external API clients

type PubMedClient struct {
	apiKey     string
	httpClient *http.Client
}

func NewPubMedClient(apiKey string, httpClient *http.Client) *PubMedClient {
	return &PubMedClient{apiKey: apiKey, httpClient: httpClient}
}

type ORCIDClient struct {
	apiKey     string
	httpClient *http.Client
}

func NewORCIDClient(apiKey string, httpClient *http.Client) *ORCIDClient {
	return &ORCIDClient{apiKey: apiKey, httpClient: httpClient}
}

type CrossRefClient struct {
	apiKey     string
	httpClient *http.Client
}

func NewCrossRefClient(apiKey string, httpClient *http.Client) *CrossRefClient {
	return &CrossRefClient{apiKey: apiKey, httpClient: httpClient}
}

// SearchArticles queries PubMed using the NCBI E-utilities API (eSearch + eSummary).
func (c *PubMedClient) SearchArticles(ctx context.Context, query string) ([]*Article, error) {
	searchURL := "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=20&retmode=json&term=" + url.QueryEscape(query)
	if c.apiKey != "" {
		searchURL += "&api_key=" + url.QueryEscape(c.apiKey)
	}

	req, err := http.NewRequestWithContext(ctx, http.MethodGet, searchURL, nil)
	if err != nil {
		return nil, fmt.Errorf("pubmed_search: %w", err)
	}
	resp, err := c.httpClient.Do(req)
	if err != nil {
		return nil, fmt.Errorf("pubmed_search: esearch: %w", err)
	}
	defer resp.Body.Close()

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("pubmed_search: read esearch: %w", err)
	}
	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("pubmed_search: esearch status %d", resp.StatusCode)
	}

	var searchResp ncbiESearchResponse
	if err := json.Unmarshal(body, &searchResp); err != nil {
		return nil, fmt.Errorf("pubmed_search: parse esearch: %w", err)
	}
	ids := searchResp.ESearchResult.IDList
	if len(ids) == 0 {
		return nil, nil
	}

	summaryURL := "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=" + url.QueryEscape(strings.Join(ids, ","))
	if c.apiKey != "" {
		summaryURL += "&api_key=" + url.QueryEscape(c.apiKey)
	}

	req2, err := http.NewRequestWithContext(ctx, http.MethodGet, summaryURL, nil)
	if err != nil {
		return nil, fmt.Errorf("pubmed_search: %w", err)
	}
	resp2, err := c.httpClient.Do(req2)
	if err != nil {
		return nil, fmt.Errorf("pubmed_search: esummary: %w", err)
	}
	defer resp2.Body.Close()

	body2, err := io.ReadAll(resp2.Body)
	if err != nil {
		return nil, fmt.Errorf("pubmed_search: read esummary: %w", err)
	}
	if resp2.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("pubmed_search: esummary status %d", resp2.StatusCode)
	}

	var summaryResp ncbiESummaryResponse
	if err := json.Unmarshal(body2, &summaryResp); err != nil {
		return nil, fmt.Errorf("pubmed_search: parse esummary: %w", err)
	}

	articles := make([]*Article, 0, len(ids))
	for _, id := range ids {
		doc, ok := summaryResp.Result[id]
		if !ok {
			continue
		}
		authors := make([]string, 0, len(doc.Authors))
		for _, a := range doc.Authors {
			authors = append(authors, a.Name)
		}
		year := 0
		if len(doc.PubDate) >= 4 {
			fmt.Sscanf(doc.PubDate[:4], "%d", &year)
		}
		articles = append(articles, &Article{
			PMID:    doc.UID,
			Title:   doc.Title,
			Authors: authors,
			Journal: doc.Source,
			Year:    year,
			DOI:     doc.DOI,
			URL:     fmt.Sprintf("https://pubmed.ncbi.nlm.nih.gov/%s/", doc.UID),
		})
	}
	return articles, nil
}

// LookupResearcher retrieves an ORCID researcher profile by ORCID iD via the public API.
func (c *ORCIDClient) LookupResearcher(ctx context.Context, orcidID string) (*Researcher, error) {
	apiURL := "https://pub.orcid.org/v3.0/" + url.PathEscape(orcidID) + "/record"

	req, err := http.NewRequestWithContext(ctx, http.MethodGet, apiURL, nil)
	if err != nil {
		return nil, fmt.Errorf("orcid_lookup: %w", err)
	}
	req.Header.Set("Accept", "application/json")
	if c.apiKey != "" {
		req.Header.Set("Authorization", "Bearer "+c.apiKey)
	}

	resp, err := c.httpClient.Do(req)
	if err != nil {
		return nil, fmt.Errorf("orcid_lookup: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode == http.StatusNotFound {
		return nil, fmt.Errorf("orcid_lookup: researcher %s not found", orcidID)
	}
	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("orcid_lookup: status %d", resp.StatusCode)
	}

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("orcid_lookup: read response: %w", err)
	}

	var record orcidRecord
	if err := json.Unmarshal(body, &record); err != nil {
		return nil, fmt.Errorf("orcid_lookup: parse response: %w", err)
	}

	name := strings.TrimSpace(record.Person.Name.GivenNames.Value + " " + record.Person.Name.FamilyName.Value)

	email := ""
	for _, em := range record.Person.Emails.Email {
		if em.Primary {
			email = em.Email
			break
		}
	}
	if email == "" && len(record.Person.Emails.Email) > 0 {
		email = record.Person.Emails.Email[0].Email
	}

	institutions := make([]string, 0)
	for _, group := range record.ActivitiesSummary.Employments.AffiliationGroup {
		for _, summary := range group.Summaries {
			if summary.EmploymentSummary.Organization.Name != "" {
				institutions = append(institutions, summary.EmploymentSummary.Organization.Name)
			}
		}
	}

	return &Researcher{
		ORCID:        orcidID,
		Name:         name,
		Email:        email,
		Institutions: institutions,
		Publications: len(record.ActivitiesSummary.Works.Group),
		Biography:    record.Person.Biography.Content,
	}, nil
}

// ResolveDOI looks up publication metadata via the CrossRef REST API.
func (c *CrossRefClient) ResolveDOI(ctx context.Context, doi string) (*Publication, error) {
	apiURL := "https://api.crossref.org/works/" + url.PathEscape(doi)

	req, err := http.NewRequestWithContext(ctx, http.MethodGet, apiURL, nil)
	if err != nil {
		return nil, fmt.Errorf("crossref_resolve: %w", err)
	}
	req.Header.Set("User-Agent", "phantom-science/1.0 (mailto:research@phantom-science.com)")
	if c.apiKey != "" {
		req.Header.Set("Crossref-Plus-API-Token", "Bearer "+c.apiKey)
	}

	resp, err := c.httpClient.Do(req)
	if err != nil {
		return nil, fmt.Errorf("crossref_resolve: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode == http.StatusNotFound {
		return nil, fmt.Errorf("crossref_resolve: DOI %s not found", doi)
	}
	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("crossref_resolve: status %d", resp.StatusCode)
	}

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("crossref_resolve: read response: %w", err)
	}

	var result crossRefResponse
	if err := json.Unmarshal(body, &result); err != nil {
		return nil, fmt.Errorf("crossref_resolve: parse response: %w", err)
	}

	w := result.Message
	title := ""
	if len(w.Title) > 0 {
		title = w.Title[0]
	}
	journal := ""
	if len(w.ContainerTitle) > 0 {
		journal = w.ContainerTitle[0]
	}
	year := 0
	if len(w.Created.DateParts) > 0 && len(w.Created.DateParts[0]) > 0 {
		year = w.Created.DateParts[0][0]
	}
	authors := make([]string, 0, len(w.Author))
	for _, a := range w.Author {
		if name := strings.TrimSpace(a.Given + " " + a.Family); name != "" {
			authors = append(authors, name)
		}
	}

	return &Publication{
		DOI:       w.DOI,
		Title:     title,
		Authors:   authors,
		Journal:   journal,
		Year:      year,
		Publisher: w.Publisher,
		URL:       w.URL,
		Type:      w.Type,
		Score:     w.Score,
	}, nil
}

type InstitutionAPI struct {
	config     InstitutionConfig
	httpClient *http.Client
}

func NewInstitutionAPI(config InstitutionConfig, httpClient *http.Client) (*InstitutionAPI, error) {
	return &InstitutionAPI{
		config:     config,
		httpClient: httpClient,
	}, nil
}

// SearchExperts queries the institution's REST API for candidates matching criteria.
// Returns an empty slice (no error) when no APIEndpoint is configured.
func (ia *InstitutionAPI) SearchExperts(ctx context.Context, criteria ExpertCriteria) ([]*ExpertCandidate, error) {
	if ia.config.APIEndpoint == "" {
		return []*ExpertCandidate{}, nil
	}

	params := url.Values{}
	if len(criteria.ResearchAreas) > 0 {
		params.Set("research_areas", strings.Join(criteria.ResearchAreas, ","))
	}
	if criteria.ExpertiseLevel != "" {
		params.Set("expertise_level", criteria.ExpertiseLevel)
	}
	if criteria.Publications > 0 {
		params.Set("min_publications", fmt.Sprintf("%d", criteria.Publications))
	}
	if criteria.ExperienceYears > 0 {
		params.Set("min_experience_years", fmt.Sprintf("%d", criteria.ExperienceYears))
	}

	apiURL := ia.config.APIEndpoint + "/experts?" + params.Encode()
	req, err := http.NewRequestWithContext(ctx, http.MethodGet, apiURL, nil)
	if err != nil {
		return nil, fmt.Errorf("institution_search_experts: %w", err)
	}
	if ia.config.APIKey != "" {
		req.Header.Set("Authorization", "Bearer "+ia.config.APIKey)
	}

	resp, err := ia.httpClient.Do(req)
	if err != nil {
		return nil, fmt.Errorf("institution_search_experts: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("institution_search_experts: status %d", resp.StatusCode)
	}

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("institution_search_experts: read response: %w", err)
	}

	var candidates []*ExpertCandidate
	if err := json.Unmarshal(body, &candidates); err != nil {
		return nil, fmt.Errorf("institution_search_experts: parse response: %w", err)
	}
	return candidates, nil
}

// CheckSubmissionStatus queries the institution's API for a journal submission's current status.
func (ia *InstitutionAPI) CheckSubmissionStatus(ctx context.Context, submissionID string) (*SubmissionStatus, error) {
	if ia.config.APIEndpoint == "" {
		return nil, fmt.Errorf("institution_check_status: APIEndpoint not configured")
	}

	apiURL := ia.config.APIEndpoint + "/submissions/" + url.PathEscape(submissionID)
	req, err := http.NewRequestWithContext(ctx, http.MethodGet, apiURL, nil)
	if err != nil {
		return nil, fmt.Errorf("institution_check_status: %w", err)
	}
	if ia.config.APIKey != "" {
		req.Header.Set("Authorization", "Bearer "+ia.config.APIKey)
	}

	resp, err := ia.httpClient.Do(req)
	if err != nil {
		return nil, fmt.Errorf("institution_check_status: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode == http.StatusNotFound {
		return nil, fmt.Errorf("institution_check_status: submission %s not found", submissionID)
	}
	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("institution_check_status: status %d", resp.StatusCode)
	}

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("institution_check_status: read response: %w", err)
	}

	var status SubmissionStatus
	if err := json.Unmarshal(body, &status); err != nil {
		return nil, fmt.Errorf("institution_check_status: parse response: %w", err)
	}
	return &status, nil
}

type ManuscriptManager struct {
	config ManuscriptConfig
}

func NewManuscriptManager(config ManuscriptConfig) *ManuscriptManager {
	return &ManuscriptManager{config: config}
}

func (mm *ManuscriptManager) FormatForJournal(manuscript *ManuscriptDraft, journal *JournalConfig) (*FormattedManuscript, error) {
	// Stub: journal requirements are not applied; the content is a placeholder.
	return &FormattedManuscript{
		Content:     "Formatted manuscript content",
		Format:      journal.PreferredFormat,
		WordCount:   5000,
		FigureCount: len(manuscript.Figures),
		TableCount:  len(manuscript.Tables),
	}, nil
}

type CollaborationPlatform struct {
	config CollaborationConfig
}

func NewCollaborationPlatform(config CollaborationConfig) *CollaborationPlatform {
	return &CollaborationPlatform{config: config}
}

func (cp *CollaborationPlatform) TrackInvitation(candidate *ExpertCandidate, invitation *ExpertInvitation) {
	// Stub: nothing is tracked.
}

func (cp *CollaborationPlatform) CreateChannel(ctx context.Context, workspaceID, channelName string) (*Channel, error) {
	return &Channel{
		ID:          fmt.Sprintf("channel_%s_%d", channelName, time.Now().Unix()),
		Name:        channelName,
		Description: fmt.Sprintf("Channel for %s discussions", channelName),
		CreatedAt:   time.Now(),
	}, nil
}

func (cp *CollaborationPlatform) InviteMember(ctx context.Context, workspaceID string, member *WorkspaceMember) error {
	// Stub: nobody is invited.
	return nil
}

func (cp *CollaborationPlatform) UploadFile(ctx context.Context, workspaceID, fileName string) (*SharedFile, error) {
	return &SharedFile{
		ID:          fmt.Sprintf("file_%s_%d", fileName, time.Now().Unix()),
		Name:        fileName,
		Description: fmt.Sprintf("Validation material: %s", fileName),
		FilePath:    fmt.Sprintf("/shared/%s", fileName),
		UploadedAt:  time.Now(),
	}, nil
}

// Manuscript text generators. These return fixed template text; only the
// expert count and authenticity score are taken from the report.

func (api *AcademicIntegrationAPI) generateManuscriptTitle(report *OverallValidationReport) string {
	return "Expert Validation of Population-Specific Genomic Sequence Generation Using Hidden Markov Models Trained on 1000 Genomes Project Data"
}

func (api *AcademicIntegrationAPI) generateAbstract(report *OverallValidationReport) string {
	return fmt.Sprintf("Background: Synthetic genomic sequences can stand in for real data when real data cannot be shared. Methods: Sequences were generated with k-mer hidden Markov models trained on 1000 Genomes Project data, and %d genomics specialists assessed them blind. Results: The authenticity score, (1 - expert consensus accuracy) x 100, was %.1f%%; 50%% corresponds to chance-level discrimination. Conclusions: [to be written from the results]", report.TotalExperts, report.AuthenticityScore)
}

func (api *AcademicIntegrationAPI) generateAuthorList(report *OverallValidationReport) []*Author {
	return []*Author{
		{
			Name:         "Research Team",
			Email:        "research@phantom-science.com",
			Institution:  "Phantom Science",
			Contribution: []string{"Conceptualization", "Methodology", "Software", "Validation"},
		},
	}
}

func (api *AcademicIntegrationAPI) generateIntroduction(report *OverallValidationReport) string {
	return "Synthetic genomic sequences are used as test data for pipelines and methods when real data cannot be shared. [Introduction to be written.]"
}

func (api *AcademicIntegrationAPI) generateMethods(report *OverallValidationReport) string {
	return "Sequences were generated by k-mer hidden Markov models trained on variant data from the 1000 Genomes Project, and assessed in three phases: blind assessment, reference comparison and expert discussion. [Methods to be completed.]"
}

func (api *AcademicIntegrationAPI) generateResults(report *OverallValidationReport) string {
	return fmt.Sprintf("%d genomics specialists took part. The authenticity score, (1 - expert consensus accuracy) x 100, was %.1f%%; 50%% corresponds to chance-level discrimination. [Results to be completed.]", report.TotalExperts, report.AuthenticityScore)
}

func (api *AcademicIntegrationAPI) generateDiscussion(report *OverallValidationReport) string {
	return "[Discussion to be written from the results.]"
}

func (api *AcademicIntegrationAPI) generateConclusion(report *OverallValidationReport) string {
	return "[Conclusion to be written from the results.]"
}

func (api *AcademicIntegrationAPI) generateFigures(report *OverallValidationReport) []*Figure {
	return []*Figure{
		{
			ID:       "fig1",
			Title:    "Expert Validation Results",
			Caption:  "Distribution of authenticity scores across expert evaluators",
			FilePath: "/figures/expert_validation_results.png",
			Type:     "png",
		},
	}
}

func (api *AcademicIntegrationAPI) generateTables(report *OverallValidationReport) []*Table {
	return []*Table{
		{
			ID:      "table1",
			Title:   "Population-Specific Validation Results",
			Caption: "Summary of expert validation results by population",
			Headers: []string{"Population", "Sequences", "Mean Score", "Consensus"},
			// Rows are not filled from the report yet.
			Rows: [][]string{},
		},
	}
}

func (api *AcademicIntegrationAPI) generateReferences(report *OverallValidationReport) []*Reference {
	return []*Reference{
		{
			ID:      "ref1",
			Authors: []string{"1000 Genomes Project Consortium"},
			Title:   "A global reference for human genetic variation",
			Journal: "Nature",
			Year:    2015,
			DOI:     "10.1038/nature15393",
		},
	}
}

func (api *AcademicIntegrationAPI) generateSupplementaryMaterials(report *OverallValidationReport) []*SupplementaryMaterial {
	return []*SupplementaryMaterial{
		{
			ID:          "supp1",
			Title:       "Expert Validation Package",
			Description: "Complete validation materials used by expert evaluators",
			FilePath:    "/supplementary/expert_validation_package.zip",
			Type:        "zip",
		},
	}
}

func (api *AcademicIntegrationAPI) generateKeywords(report *OverallValidationReport) []string {
	return []string{
		"population genomics",
		"sequence generation",
		"hidden markov models",
		"expert validation",
		"precision medicine",
		"1000 genomes project",
	}
}

func (api *AcademicIntegrationAPI) selectTargetJournals(report *OverallValidationReport) []*JournalConfig {
	// Select journals based on validation results quality
	var count int
	if report.AuthenticityScore >= 75.0 {
		count = min(3, len(api.config.ManuscriptConfig.TargetJournals))
	} else {
		count = min(2, len(api.config.ManuscriptConfig.TargetJournals))
	}
	result := make([]*JournalConfig, count)
	for i := 0; i < count; i++ {
		j := api.config.ManuscriptConfig.TargetJournals[i]
		result[i] = &j
	}
	return result
}

func (api *AcademicIntegrationAPI) submitToJournal(ctx context.Context, manuscript *FormattedManuscript, journal *JournalConfig) (*SubmissionResult, error) {
	// Stub: nothing is submitted.
	return &SubmissionResult{
		Journal:         journal.Name,
		Status:          "submitted",
		SubmissionID:    fmt.Sprintf("sub_%d", time.Now().Unix()),
		TrackingURL:     fmt.Sprintf("https://%s.com/track/submission", journal.Name),
		SubmittedAt:     time.Now(),
		EstimatedReview: journal.ReviewTimeline,
	}, nil
}

func (api *AcademicIntegrationAPI) sendEmail(ctx context.Context, payload EmailPayload) error {
	// Stub: no email is sent.
	return nil
}

func (api *AcademicIntegrationAPI) formatInvitationEmail(invitation *ExpertInvitation) string {
	return fmt.Sprintf(`
		<html>
		<body>
		<h2>Expert Validation Invitation</h2>
		<p>%s</p>
		<h3>Project Details:</h3>
		<p><strong>Title:</strong> %s</p>
		<p><strong>Description:</strong> %s</p>
		<p><strong>Time Commitment:</strong> %s</p>
		<p><strong>Deadline:</strong> %s</p>
		</body>
		</html>
	`, invitation.PersonalizedMessage, invitation.ProjectDetails.Title, invitation.ProjectDetails.Description, invitation.TimeCommitment, invitation.Deadline.Format("January 2, 2006"))
}

func (api *AcademicIntegrationAPI) getContactInformation() *ContactInfo {
	return &ContactInfo{
		Name:        "Research Team",
		Email:       "research@phantom-science.com",
		Phone:       "+1-555-0123",
		Institution: "Phantom Science",
	}
}
