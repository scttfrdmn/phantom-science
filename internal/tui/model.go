package tui

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"time"

	"github.com/charmbracelet/bubbletea"
	"github.com/charmbracelet/lipgloss"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/timing"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Model represents the main TUI state
type Model struct {
	// Configuration
	sequencerConfig types.SequencerConfig
	phantomConfig   *types.PhantomConfig

	// UI State
	currentSection    int
	sequencerChoice   int
	compression       timing.TimeCompression
	samples           []Sample
	sampleInput       SampleInput
	editingMode       bool
	
	// Fingerprinting configuration
	fingerprintingEnabled bool
	customerID           string
	embeddingStrategy    int // 0=Conservative, 1=Balanced, 2=Aggressive

	// HMM model configuration (mirrors sequencerConfig fields for UI state)
	hmmModelPath string // Path to a serialized .hmm.json model (empty = crypto/rand)
	population   string // 1000 Genomes population code (e.g. YRI, CEU, CHB)

	// Simulation state
	simulator       *timing.RunSimulator
	genomicsEngine  *genomics.Engine
	runInProgress   bool
	progress        timing.RunProgress
	ctx             context.Context
	cancel          context.CancelFunc

	// UI dimensions
	width  int
	height int

	// Styles
	styles Styles
}

// Sample represents a sample in the sample sheet
type Sample struct {
	ID          string
	Lane        int
	Index1      string
	Index2      string
	Project     string
	Description string
}

// SampleInput represents the current sample being edited
type SampleInput struct {
	ID          string
	Lane        int
	Index1      string
	Index2      string
	Project     string
	Description string
	Field       int // Current field being edited
}

// Sections of the TUI
const (
	SectionSequencer = iota
	SectionSamples
	SectionParameters
	SectionFingerprinting
	SectionOutput
	SectionActions
)

// Sequencer choices
const (
	ChoiceNovaSeq = iota
	ChoiceNextSeq
	ChoiceMiSeq
)

// Embedding strategy choices
const (
	StrategyConservative = iota
	StrategyBalanced
	StrategyAggressive
)

// Message types
type (
	tickMsg              struct{}
	progressMsg          timing.RunProgress
	fileGeneratedMsg     struct{ filename string }
	simulationCompleteMsg struct{}
	simulationErrorMsg   struct{ err error }
)

// New creates a new TUI model
func New(sequencerConfig types.SequencerConfig, phantomConfig *types.PhantomConfig) Model {
	// Set default phantom config if none provided
	if phantomConfig == nil {
		phantomConfig = &types.PhantomConfig{
			Enterprise: types.EnterpriseConfig{
				CustomerID:   "DEMO",
				ServiceTier:  "standard",
				Fingerprinting: types.FingerprintConfig{
					Enable: false,
				},
			},
			Validation: types.ValidationConfig{
				Enable: false,
			},
		}
	}

	// Set default samples
	samples := []Sample{
		{
			ID:          "Patient_001",
			Lane:        1,
			Index1:      "ATCGATCG",
			Index2:      "CGTAGATC",
			Project:     "Cancer_Study",
			Description: "Tumor sample",
		},
		{
			ID:          "Patient_002",
			Lane:        1,
			Index1:      "GCTAGCTA",
			Index2:      "TACGACGT",
			Project:     "Cancer_Study",
			Description: "Normal sample",
		},
	}

	return Model{
		sequencerConfig:       sequencerConfig,
		phantomConfig:         phantomConfig,
		currentSection:        SectionSequencer,
		sequencerChoice:       ChoiceMiSeq,
		compression:           timing.Compressed100x,
		samples:               samples,
		sampleInput:           SampleInput{Lane: 1},
		fingerprintingEnabled: phantomConfig.Enterprise.Fingerprinting.Enable,
		customerID:            phantomConfig.Enterprise.CustomerID,
		embeddingStrategy:     StrategyBalanced, // Default to balanced strategy
		hmmModelPath:          sequencerConfig.HMMModelPath,
		population:            sequencerConfig.Population,
		styles:                NewStyles(),
	}
}

// Init implements bubbletea.Model
func (m Model) Init() tea.Cmd {
	return tea.Batch(
		tea.EnterAltScreen,
		tickCmd(),
	)
}

// Update implements bubbletea.Model
func (m Model) Update(msg tea.Msg) (tea.Model, tea.Cmd) {
	switch msg := msg.(type) {
	case tea.WindowSizeMsg:
		m.width = msg.Width
		m.height = msg.Height
		return m, nil

	case tea.KeyMsg:
		return m.handleKeyPress(msg)

	case tickMsg:
		if m.runInProgress && m.simulator != nil {
			// Update simulation progress here
			return m, tickCmd()
		}
		return m, tickCmd()

	case progressMsg:
		m.progress = timing.RunProgress(msg)
		return m, nil

	case fileGeneratedMsg:
		// Update file count and current file name
		m.progress.FilesGenerated++
		m.progress.CurrentFile = msg.filename
		return m, nil

	case simulationCompleteMsg:
		m.runInProgress = false
		m.progress.Stage = timing.StageComplete
		m.progress.Percent = 100.0
		if m.cancel != nil {
			m.cancel()
		}
		return m, nil

	case simulationErrorMsg:
		m.runInProgress = false
		if m.cancel != nil {
			m.cancel()
		}
		// Could show error in UI
		return m, nil
	}

	return m, nil
}

// View implements bubbletea.Model
func (m Model) View() string {
	if m.width == 0 || m.height == 0 {
		return "Loading..."
	}

	// Build the main interface
	header := m.renderHeader()
	topPanels := m.renderTopPanels()
	sampleSheet := m.renderSampleSheet()
	bottomPanels := m.renderBottomPanels()
	preview := m.renderPreview()
	actions := m.renderActions()

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		header,
		"",
		topPanels,
		"",
		sampleSheet,
		"",
		bottomPanels,
		"",
		preview,
		"",
		actions,
	)

	// Add border around everything
	return m.styles.MainBorder.Render(content)
}

// Helper functions
func (m Model) handleKeyPress(msg tea.KeyMsg) (Model, tea.Cmd) {
	switch msg.String() {
	case "ctrl+c", "q":
		if m.runInProgress && m.cancel != nil {
			m.cancel()
		}
		return m, tea.Quit

	case "tab":
		if !m.editingMode {
			m.currentSection = (m.currentSection + 1) % 6
		}
		return m, nil

	case "shift+tab":
		if !m.editingMode {
			m.currentSection = (m.currentSection - 1 + 6) % 6
		}
		return m, nil

	case "enter":
		return m.handleEnter()

	case "esc":
		if m.editingMode {
			m.editingMode = false
		}
		return m, nil

	case "up", "down", "left", "right":
		return m.handleArrowKeys(msg.String())

	case " ":
		return m.handleSpace()

	default:
		if m.editingMode {
			return m.handleTextInput(msg.String())
		}
	}

	return m, nil
}

func (m Model) handleEnter() (Model, tea.Cmd) {
	switch m.currentSection {
	case SectionActions:
		return m.handleActionSelection()
	case SectionSamples:
		if !m.editingMode {
			m.editingMode = true
			m.sampleInput = SampleInput{Lane: 1}
		}
		return m, nil
	}
	return m, nil
}

func (m Model) handleActionSelection() (Model, tea.Cmd) {
	if m.runInProgress {
		// Stop simulation
		if m.cancel != nil {
			m.cancel()
		}
		m.runInProgress = false
		return m, nil
	}

	// Start simulation
	m = m.updateConfigFromUI()
	
	ctx, cancel := context.WithCancel(context.Background())
	m.ctx = ctx
	m.cancel = cancel

	// Create output directory if it doesn't exist
	if err := os.MkdirAll(m.sequencerConfig.OutputDir, 0750); err != nil {
		// Handle error - could show in UI
		return m, nil
	}

	// Initialize genomics engine for actual FASTQ generation
	m.genomicsEngine = genomics.NewEngine(m.sequencerConfig)
	
	simulator := timing.NewRunSimulator(m.sequencerConfig, m.compression)
	
	// Set up callback to actually generate FASTQ files
	simulator.FileCallback = func(filename string, batchNum int, readNum int) {
		fullPath := filepath.Join(m.sequencerConfig.OutputDir, filename)
		
		// Generate the actual FASTQ file
		if err := m.genomicsEngine.GenerateFASTQFile(fullPath, batchNum, readNum); err != nil {
			// Log error but continue simulation
			fmt.Printf("Error generating FASTQ file %s: %v\n", filename, err)
		}
	}
	
	simulator.ProgressCallback = func(stage string, percent float64, eta time.Duration) {
		// Progress updates are handled by the timing simulator
		// and displayed in the UI automatically
	}

	m.simulator = simulator
	m.runInProgress = true
	
	// Initialize progress
	totalFiles, totalBytes := simulator.EstimateTotalOutput()
	m.progress = timing.RunProgress{
		Stage:          timing.StageInitialization,
		Percent:        0.0,
		ETA:            time.Duration(simulator.Timing.TotalRunTimeHours) * time.Hour / time.Duration(m.compression),
		FilesGenerated: 0,
		TotalFiles:     totalFiles,
		BytesGenerated: 0,
		TotalBytes:     totalBytes,
		CurrentFile:    "",
	}

	// Start simulation in background
	return m, startSimulation(ctx, simulator)
}

func (m Model) handleArrowKeys(key string) (Model, tea.Cmd) {
	switch m.currentSection {
	case SectionSequencer:
		switch key {
		case "up":
			if m.sequencerChoice > 0 {
				m.sequencerChoice--
			}
		case "down":
			if m.sequencerChoice < 2 {
				m.sequencerChoice++
			}
		}
	}
	return m, nil
}

func (m Model) handleSpace() (Model, tea.Cmd) {
	switch m.currentSection {
	case SectionParameters:
		// Toggle paired-end, compression, etc.
		return m, nil
	}
	return m, nil
}

func (m Model) handleTextInput(input string) (Model, tea.Cmd) {
	// Handle text input for sample editing
	return m, nil
}

func (m Model) updateConfigFromUI() Model {
	// Update sequencer model based on choice
	switch m.sequencerChoice {
	case ChoiceNovaSeq:
		m.sequencerConfig.SequencerModel = "NovaSeq 6000"
		m.sequencerConfig.FlowCellType = "S2"
	case ChoiceNextSeq:
		m.sequencerConfig.SequencerModel = "NextSeq 2000"
		m.sequencerConfig.FlowCellType = "P3"
	case ChoiceMiSeq:
		m.sequencerConfig.SequencerModel = "MiSeq"
		m.sequencerConfig.FlowCellType = "v3"
	}

	// Update other parameters from UI state
	if len(m.samples) > 0 {
		m.sequencerConfig.SampleName = m.samples[0].ID
	}
	m.sequencerConfig.TotalBatches = len(m.samples)

	// Sync HMM model fields back to config
	m.sequencerConfig.HMMModelPath = m.hmmModelPath
	m.sequencerConfig.Population = m.population

	return m
}

// Tick command for regular updates
func tickCmd() tea.Cmd {
	return tea.Tick(time.Millisecond*100, func(t time.Time) tea.Msg {
		return tickMsg{}
	})
}

// Start simulation command
func startSimulation(ctx context.Context, simulator *timing.RunSimulator) tea.Cmd {
	return tea.Cmd(func() tea.Msg {
		err := simulator.Run(ctx)
		if err != nil {
			return simulationErrorMsg{err: err}
		}
		return simulationCompleteMsg{}
	})
}