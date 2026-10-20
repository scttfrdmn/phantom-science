// Package tui is the interactive terminal UI (phantom-science tui). It paces
// FASTQ generation with the run timeline from package timing and shows run
// progress, the sample being generated and the fingerprinting settings.
package tui

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sync"
	"time"

	"github.com/charmbracelet/bubbletea"
	"github.com/charmbracelet/lipgloss"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/timing"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Model represents the main TUI state
type Model struct {
	version string // shown in the header

	// Configuration
	sequencerConfig types.SequencerConfig
	phantomConfig   *types.PhantomConfig

	// UI State
	currentSection int
	// sequencerChoice is the highlighted preset (ChoiceNovaSeq..ChoiceMiSeq),
	// or ChoiceConfig when the configured sequencer is none of the presets.
	sequencerChoice int
	// sequencerChanged is set once the user picks a preset; until then the
	// configured sequencer_model and flow_cell_type are used unchanged.
	sequencerChanged bool
	compression      timing.TimeCompression
	samples          []Sample

	// Fingerprinting settings the run will use (after ResolveFingerprinting).
	fingerprintingEnabled bool
	customerID            string
	fingerprintLayers     []string

	// Simulation state
	simulator      *timing.RunSimulator
	genomicsEngine *genomics.Engine
	runInProgress  bool
	runID          int       // identifies the current run's messages
	run            *runState // shared with the simulator goroutine; nil before the first run
	progress       timing.RunProgress
	err            error // last run error, shown in View
	ctx            context.Context
	cancel         context.CancelFunc

	// UI dimensions
	width  int
	height int

	// Styles
	styles Styles
}

// Sample is one row of the sample sheet.
type Sample struct {
	ID          string
	Lane        int
	Index1      string
	Index2      string
	Project     string
	Description string
}

// Sections of the TUI
const (
	SectionSequencer = iota
	SectionSamples
	SectionParameters
	SectionFingerprinting
	SectionOutput
	SectionActions
	numSections
)

// Sequencer choices. ChoiceConfig means the configured sequencer is not one
// of the presets and is used as configured.
const (
	ChoiceConfig  = -1
	ChoiceNovaSeq = 0
	ChoiceNextSeq = 1
	ChoiceMiSeq   = 2
)

// sequencerPresets are the sequencers selectable with the arrow keys, with
// the flow cell each one uses when picked.
var sequencerPresets = []struct {
	name     string
	flowCell string
	icon     string
}{
	{"NovaSeq 6000", "S2", "🚀"},
	{"NextSeq 2000", "P3", "⚡"},
	{"MiSeq", "v3", "🔬"},
}

// compressionOptions are the selectable time compressions, in display order.
var compressionOptions = []timing.TimeCompression{
	timing.RealTime,
	timing.Compressed10x,
	timing.Compressed100x,
	timing.Compressed1000x,
}

// stageStopped and stageFailed label runs that did not complete.
const (
	stageStopped timing.RunStage = "Stopped"
	stageFailed  timing.RunStage = "Failed"
)

// Message types
type (
	tickMsg               struct{}
	simulationCompleteMsg struct{ runID int }
	simulationErrorMsg    struct {
		runID int
		err   error
	}
)

// runState is the progress of a run, written by the simulator goroutine and
// read on each tick.
type runState struct {
	mu         sync.Mutex
	progress   timing.RunProgress
	started    time.Time
	fileErrors int
	firstErr   error
}

func (r *runState) setStage(stage string, percent float64, eta time.Duration) {
	r.mu.Lock()
	defer r.mu.Unlock()
	r.progress.Stage = timing.RunStage(stage)
	r.progress.Percent = percent
	r.progress.ETA = eta
}

func (r *runState) fileDone(path string, size int64, err error) {
	r.mu.Lock()
	defer r.mu.Unlock()
	if err != nil {
		r.fileErrors++
		if r.firstErr == nil {
			r.firstErr = fmt.Errorf("generating %s: %w", filepath.Base(path), err)
		}
		return
	}
	r.progress.FilesGenerated++
	r.progress.BytesGenerated += size
	r.progress.CurrentFile = path
}

// snapshot returns the progress and an error summarising failed files.
func (r *runState) snapshot() (timing.RunProgress, time.Duration, error) {
	r.mu.Lock()
	defer r.mu.Unlock()
	var err error
	if r.firstErr != nil {
		err = fmt.Errorf("%d file(s) failed; first: %w", r.fileErrors, r.firstErr)
	}
	return r.progress, time.Since(r.started), err
}

// New creates a new TUI model. version is shown in the header.
func New(sequencerConfig types.SequencerConfig, phantomConfig *types.PhantomConfig, version string) Model {
	resolved := sequencerConfig
	resolved.ResolveFingerprinting(phantomConfig)

	return Model{
		version:               version,
		sequencerConfig:       sequencerConfig,
		phantomConfig:         phantomConfig,
		currentSection:        SectionSequencer,
		sequencerChoice:       presetIndex(sequencerConfig.SequencerModel),
		compression:           timing.Compressed100x,
		samples:               []Sample{{ID: sequencerConfig.SampleName, Lane: sequencerConfig.Lane}},
		fingerprintingEnabled: resolved.Fingerprinting.Enable,
		customerID:            resolved.Enterprise.CustomerID,
		fingerprintLayers:     resolved.Fingerprinting.SteganographyLayers,
		styles:                NewStyles(),
	}
}

// presetIndex returns the preset matching model, or ChoiceConfig.
func presetIndex(model string) int {
	for i, p := range sequencerPresets {
		if p.name == model {
			return i
		}
	}
	return ChoiceConfig
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
		if m.runInProgress {
			m = m.syncProgress()
		}
		return m, tickCmd()

	case simulationCompleteMsg:
		if msg.runID != m.runID || !m.runInProgress {
			return m, nil // a stopped or earlier run
		}
		m = m.syncProgress()
		m.runInProgress = false
		m.progress.Stage = timing.StageComplete
		m.progress.Percent = 100.0
		m.progress.ETA = 0
		if m.cancel != nil {
			m.cancel()
		}
		return m, nil

	case simulationErrorMsg:
		if msg.runID != m.runID || !m.runInProgress {
			return m, nil // a stopped or earlier run
		}
		m = m.syncProgress()
		m.runInProgress = false
		if m.cancel != nil {
			m.cancel()
		}
		if !errors.Is(msg.err, context.Canceled) {
			m.progress.Stage = stageFailed
			m.err = msg.err
		}
		return m, nil
	}

	return m, nil
}

// syncProgress copies the run's progress into the model.
func (m Model) syncProgress() Model {
	if m.run == nil {
		return m
	}
	p, _, err := m.run.snapshot()
	totalFiles, totalBytes := m.progress.TotalFiles, m.progress.TotalBytes
	m.progress = p
	m.progress.TotalFiles, m.progress.TotalBytes = totalFiles, totalBytes
	if err != nil {
		m.err = err
	}
	return m
}

// View implements bubbletea.Model
func (m Model) View() string {
	if m.width == 0 || m.height == 0 {
		return "Loading..."
	}

	parts := []string{m.renderHeader()}
	if m.err != nil {
		parts = append(parts, "", m.renderError())
	}
	parts = append(parts,
		"",
		m.renderTopPanels(),
		"",
		m.renderSampleSheet(),
		"",
		m.renderBottomPanels(),
		"",
		m.renderPreview(),
		"",
		m.renderActions(),
	)

	return m.styles.MainBorder.Render(lipgloss.JoinVertical(lipgloss.Left, parts...))
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
		m.currentSection = (m.currentSection + 1) % numSections
		return m, nil

	case "shift+tab":
		m.currentSection = (m.currentSection - 1 + numSections) % numSections
		return m, nil

	case "enter":
		if m.currentSection == SectionActions {
			return m.handleActionSelection()
		}
		return m, nil

	case "up", "down", "left", "right":
		return m.handleArrowKeys(msg.String())
	}

	return m, nil
}

func (m Model) handleActionSelection() (Model, tea.Cmd) {
	if m.runInProgress {
		// Stop simulation
		if m.cancel != nil {
			m.cancel()
		}
		m = m.syncProgress()
		m.runInProgress = false
		m.progress.Stage = stageStopped
		m.progress.ETA = 0
		return m, nil
	}

	// Start simulation
	m = m.updateConfigFromUI()
	m.err = nil

	if err := os.MkdirAll(m.sequencerConfig.OutputDir, 0750); err != nil {
		m.err = fmt.Errorf("creating output directory: %w", err)
		return m, nil
	}

	ctx, cancel := context.WithCancel(context.Background())
	m.ctx = ctx
	m.cancel = cancel

	// Initialize genomics engine for actual FASTQ generation
	engineConfig := m.sequencerConfig
	engineConfig.ResolveFingerprinting(m.phantomConfig)
	engine := genomics.NewEngine(engineConfig)
	engine.SetToolVersion(m.version)
	m.genomicsEngine = engine

	simulator := timing.NewRunSimulator(m.sequencerConfig, m.compression)
	run := &runState{started: time.Now()}
	outputDir := m.sequencerConfig.OutputDir

	simulator.FileCallback = func(filename string, batchNum int, readNum int) {
		fullPath := filepath.Join(outputDir, filename)
		if err := engine.GenerateFASTQFile(fullPath, batchNum, readNum); err != nil {
			run.fileDone(fullPath, 0, err)
			return
		}
		// The engine may rename the file (gzip output); report the final one.
		final := fullPath
		if files := engine.GeneratedFiles(); len(files) > 0 {
			final = files[len(files)-1]
		}
		var size int64
		if fi, err := os.Stat(final); err == nil {
			size = fi.Size()
		}
		run.fileDone(final, size, nil)
	}
	simulator.ProgressCallback = run.setStage

	m.simulator = simulator
	m.run = run
	m.runID++
	m.runInProgress = true

	totalFiles, totalBytes := simulator.EstimateTotalOutput()
	m.progress = timing.RunProgress{
		Stage:      timing.StageInitialization,
		TotalFiles: totalFiles,
		TotalBytes: totalBytes,
	}

	return m, startSimulation(ctx, simulator, engine, m.runID)
}

func (m Model) handleArrowKeys(key string) (Model, tea.Cmd) {
	if m.runInProgress {
		return m, nil // the running simulation does not pick up changes
	}
	switch m.currentSection {
	case SectionSequencer:
		choice := m.sequencerChoice
		switch key {
		case "up":
			if choice > 0 {
				choice--
			}
		case "down":
			if choice < len(sequencerPresets)-1 {
				choice++
			}
		}
		if choice != m.sequencerChoice {
			m.sequencerChoice = choice
			m.sequencerChanged = true
		}
	case SectionOutput:
		i := compressionIndex(m.compression)
		switch key {
		case "left":
			if i > 0 {
				i--
			}
		case "right":
			if i < len(compressionOptions)-1 {
				i++
			}
		}
		m.compression = compressionOptions[i]
	}
	return m, nil
}

// compressionIndex returns c's position in compressionOptions (0 if absent).
func compressionIndex(c timing.TimeCompression) int {
	for i, o := range compressionOptions {
		if o == c {
			return i
		}
	}
	return 0
}

// updateConfigFromUI applies the sequencer the user picked. Until the user
// picks one, the configured sequencer_model and flow_cell_type are kept.
func (m Model) updateConfigFromUI() Model {
	if m.sequencerChanged && m.sequencerChoice >= 0 && m.sequencerChoice < len(sequencerPresets) {
		p := sequencerPresets[m.sequencerChoice]
		m.sequencerConfig.SequencerModel = p.name
		m.sequencerConfig.FlowCellType = p.flowCell
	}
	return m
}

// Tick command for regular updates
func tickCmd() tea.Cmd {
	return tea.Tick(time.Millisecond*100, func(t time.Time) tea.Msg {
		return tickMsg{}
	})
}

// startSimulation runs the simulator and closes engine when it ends.
func startSimulation(ctx context.Context, simulator *timing.RunSimulator, engine *genomics.Engine, runID int) tea.Cmd {
	return func() tea.Msg {
		err := simulator.Run(ctx)
		if engine != nil {
			engine.Close()
		}
		if err != nil {
			return simulationErrorMsg{runID: runID, err: err}
		}
		return simulationCompleteMsg{runID: runID}
	}
}
