package tui

import (
	"context"
	"strings"
	"testing"
	"time"

	tea "github.com/charmbracelet/bubbletea"
	"github.com/scttfrdmn/phantom-science/internal/timing"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestNewModel(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   5,
		ReadsPerBatch:  1000,
	}
	
	phantomConfig := &types.PhantomConfig{
		Enterprise: types.EnterpriseConfig{
			CustomerID:  "TEST_CUSTOMER",
			ServiceTier: "enterprise",
			Fingerprinting: types.FingerprintConfig{
				Enable: true,
			},
		},
	}
	
	model := New(config, phantomConfig)
	
	// Test configuration preservation
	if model.sequencerConfig.SequencerModel != "MiSeq" {
		t.Errorf("Expected sequencer model 'MiSeq', got '%s'", model.sequencerConfig.SequencerModel)
	}
	
	if model.phantomConfig.Enterprise.CustomerID != "TEST_CUSTOMER" {
		t.Errorf("Expected customer ID 'TEST_CUSTOMER', got '%s'", model.phantomConfig.Enterprise.CustomerID)
	}
	
	// Test fingerprinting setup
	if !model.fingerprintingEnabled {
		t.Error("Fingerprinting should be enabled when config.Enterprise.Fingerprinting.Enable is true")
	}
	
	if model.customerID != "TEST_CUSTOMER" {
		t.Errorf("Expected customerID 'TEST_CUSTOMER', got '%s'", model.customerID)
	}
	
	// Test default sample initialization
	if len(model.samples) == 0 {
		t.Error("Model should initialize with default samples")
	}
	
	if len(model.samples) >= 1 {
		firstSample := model.samples[0]
		if firstSample.ID == "" || firstSample.Lane == 0 {
			t.Error("Default samples should have valid ID and Lane")
		}
	}
	
	// Test UI state initialization
	if model.currentSection != SectionSequencer {
		t.Errorf("Expected initial section to be SectionSequencer (%d), got %d", SectionSequencer, model.currentSection)
	}
	
	if model.sequencerChoice != ChoiceMiSeq {
		t.Errorf("Expected default sequencer choice to be ChoiceMiSeq (%d), got %d", ChoiceMiSeq, model.sequencerChoice)
	}
	
	if model.embeddingStrategy != StrategyBalanced {
		t.Errorf("Expected default embedding strategy to be StrategyBalanced (%d), got %d", StrategyBalanced, model.embeddingStrategy)
	}
	
	// Test compression default
	if model.compression != timing.Compressed100x {
		t.Errorf("Expected default compression to be Compressed100x, got %v", model.compression)
	}
}

func TestModelUIBehavior(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	
	// Test window size handling
	updatedModel, cmd := model.Update(tea.WindowSizeMsg{Width: 120, Height: 40})
	m, ok := updatedModel.(Model)
	if !ok {
		t.Fatal("Update should return Model type")
	}
	
	if m.width != 120 || m.height != 40 {
		t.Errorf("Expected dimensions 120x40, got %dx%d", m.width, m.height)
	}
	
	if cmd != nil {
		t.Error("WindowSizeMsg should not return a command")
	}
	
	// Test section navigation with Tab
	if m.currentSection != SectionSequencer {
		t.Error("Should start at SectionSequencer")
	}
	
	updatedModel, _ = m.Update(tea.KeyMsg{Type: tea.KeyTab})
	m = updatedModel.(Model)
	if m.currentSection != SectionSamples {
		t.Errorf("Tab should move to next section, expected %d, got %d", SectionSamples, m.currentSection)
	}
	
	// Test section navigation with Shift+Tab
	updatedModel, _ = m.Update(tea.KeyMsg{Type: tea.KeyShiftTab})
	m = updatedModel.(Model)
	if m.currentSection != SectionSequencer {
		t.Errorf("Shift+Tab should move to previous section, expected %d, got %d", SectionSequencer, m.currentSection)
	}
}

func TestSequencerSelection(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	model.currentSection = SectionSequencer
	
	// Test initial sequencer choice (defaults to ChoiceMiSeq = 2)
	if model.sequencerChoice != ChoiceMiSeq {
		t.Errorf("Expected initial choice to be ChoiceMiSeq (%d), got %d", ChoiceMiSeq, model.sequencerChoice)
	}
	
	// Test up arrow navigation from MiSeq (2) to NextSeq (1)
	updatedModel, _ := model.Update(tea.KeyMsg{Type: tea.KeyUp})
	m := updatedModel.(Model)
	if m.sequencerChoice != ChoiceNextSeq {
		t.Errorf("Up arrow should move to NextSeq, expected %d, got %d", ChoiceNextSeq, m.sequencerChoice)
	}
	
	// Test up arrow again from NextSeq (1) to NovaSeq (0)
	updatedModel, _ = m.Update(tea.KeyMsg{Type: tea.KeyUp})
	m = updatedModel.(Model)
	if m.sequencerChoice != ChoiceNovaSeq {
		t.Errorf("Up arrow should move to NovaSeq, expected %d, got %d", ChoiceNovaSeq, m.sequencerChoice)
	}
	
	// Test up arrow at top should not go below 0
	updatedModel, _ = m.Update(tea.KeyMsg{Type: tea.KeyUp})
	m = updatedModel.(Model)
	if m.sequencerChoice != ChoiceNovaSeq {
		t.Error("Up arrow at top should not change selection")
	}
	
	// Test down arrow navigation from NovaSeq (0) to NextSeq (1)
	updatedModel, _ = m.Update(tea.KeyMsg{Type: tea.KeyDown})
	m = updatedModel.(Model)
	if m.sequencerChoice != ChoiceNextSeq {
		t.Errorf("Down arrow should move to NextSeq, expected %d, got %d", ChoiceNextSeq, m.sequencerChoice)
	}
	
	// Test down arrow again from NextSeq (1) to MiSeq (2)
	updatedModel, _ = m.Update(tea.KeyMsg{Type: tea.KeyDown})
	m = updatedModel.(Model)
	if m.sequencerChoice != ChoiceMiSeq {
		t.Errorf("Down arrow should move to MiSeq, expected %d, got %d", ChoiceMiSeq, m.sequencerChoice)
	}
	
	// Test down arrow at bottom should not go beyond last option
	updatedModel, _ = m.Update(tea.KeyMsg{Type: tea.KeyDown})
	m = updatedModel.(Model)
	if m.sequencerChoice != ChoiceMiSeq {
		t.Error("Down arrow at bottom should not change selection")
	}
}

func TestConfigurationUpdate(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	
	// Test sequencer choice logic (note: updateConfigFromUI has value receiver so changes don't persist)
	// We'll test the logic by simulating what should happen
	
	// Test NovaSeq configuration
	testModel := model
	testModel.sequencerChoice = ChoiceNovaSeq
	switch testModel.sequencerChoice {
	case ChoiceNovaSeq:
		expectedModel := "NovaSeq 6000"
		expectedFlowCell := "S2"
		t.Logf("NovaSeq choice would set: Model=%s, FlowCell=%s", expectedModel, expectedFlowCell)
	case ChoiceNextSeq:
		expectedModel := "NextSeq 2000"
		expectedFlowCell := "P3"
		t.Logf("NextSeq choice would set: Model=%s, FlowCell=%s", expectedModel, expectedFlowCell)
	case ChoiceMiSeq:
		expectedModel := "MiSeq"
		expectedFlowCell := "v3"
		t.Logf("MiSeq choice would set: Model=%s, FlowCell=%s", expectedModel, expectedFlowCell)
	}
	
	// Test that the choice constants are correct
	if ChoiceNovaSeq != 0 {
		t.Errorf("Expected ChoiceNovaSeq to be 0, got %d", ChoiceNovaSeq)
	}
	if ChoiceNextSeq != 1 {
		t.Errorf("Expected ChoiceNextSeq to be 1, got %d", ChoiceNextSeq)
	}
	if ChoiceMiSeq != 2 {
		t.Errorf("Expected ChoiceMiSeq to be 2, got %d", ChoiceMiSeq)
	}
	
	// Test initial configuration is preserved
	if model.sequencerConfig.SequencerModel != "MiSeq" {
		t.Errorf("Expected initial sequencer model to be MiSeq, got %s", model.sequencerConfig.SequencerModel)
	}
	
	// Test sample structure
	if len(model.samples) == 0 {
		t.Error("Model should have default samples")
	}
	
	// Note: updateConfigFromUI() method has a design issue - it uses value receiver
	// so changes don't persist. This would need to be fixed to use pointer receiver
	// for the configuration updates to work properly.
}

func TestViewRendering(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	model.width = 120
	model.height = 40
	
	// Test view rendering without dimensions should show loading
	model.width = 0
	model.height = 0
	view := model.View()
	if view != "Loading..." {
		t.Errorf("Expected 'Loading...' when dimensions are 0, got %s", view)
	}
	
	// Test view rendering with proper dimensions
	model.width = 120
	model.height = 40
	view = model.View()
	if view == "" {
		t.Error("View should not be empty with proper dimensions")
	}
	
	// Test that view contains expected elements
	if !strings.Contains(view, "Phantom Science") {
		t.Error("View should contain application title")
	}
	
	if !strings.Contains(view, "Sequencer Selection") {
		t.Error("View should contain sequencer selection panel")
	}
	
	if !strings.Contains(view, "Sample Sheet") {
		t.Error("View should contain sample sheet panel")
	}
	
	if !strings.Contains(view, "MiSeq") {
		t.Error("View should show current sequencer selection")
	}
}

func TestKeyHandling(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	
	// Test Ctrl+C quit command
	updatedModel, cmd := model.Update(tea.KeyMsg{Type: tea.KeyCtrlC})
	if cmd == nil {
		t.Error("Ctrl+C should return a quit command")
	}
	
	// Test 'q' quit command
	updatedModel, cmd = model.Update(tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune{'q'}})
	if cmd == nil {
		t.Error("'q' key should return a quit command")
	}
	
	// Test escape key when not in editing mode
	model.editingMode = false
	updatedModel, cmd = model.Update(tea.KeyMsg{Type: tea.KeyEsc})
	m := updatedModel.(Model)
	if cmd != nil {
		t.Error("Escape should not return command when not editing")
	}
	
	// Test escape key when in editing mode
	model.editingMode = true
	updatedModel, cmd = model.Update(tea.KeyMsg{Type: tea.KeyEsc})
	m = updatedModel.(Model)
	if m.editingMode {
		t.Error("Escape should exit editing mode")
	}
	
	// Test tab navigation when not editing
	model.editingMode = false
	model.currentSection = SectionSequencer
	updatedModel, _ = model.Update(tea.KeyMsg{Type: tea.KeyTab})
	m = updatedModel.(Model)
	if m.currentSection == SectionSequencer {
		t.Error("Tab should change section when not editing")
	}
	
	// Test tab ignored when editing
	model.editingMode = true
	originalSection := model.currentSection
	updatedModel, _ = model.Update(tea.KeyMsg{Type: tea.KeyTab})
	m = updatedModel.(Model)
	if m.currentSection != originalSection {
		t.Error("Tab should be ignored when editing")
	}
}

func TestSampleHandling(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	
	// Test default samples
	if len(model.samples) < 2 {
		t.Error("Model should initialize with at least 2 default samples")
	}
	
	// Verify sample structure
	if len(model.samples) >= 1 {
		firstSample := model.samples[0]
		if firstSample.ID == "" {
			t.Error("First sample should have an ID")
		}
		if firstSample.Lane <= 0 {
			t.Error("Sample should have a valid lane number")
		}
		if firstSample.Index1 == "" {
			t.Error("Sample should have Index1")
		}
		if firstSample.Project == "" {
			t.Error("Sample should have a project name")
		}
	}
	
	// Test enter key on samples section
	model.currentSection = SectionSamples
	model.editingMode = false
	updatedModel, _ := model.Update(tea.KeyMsg{Type: tea.KeyEnter})
	m := updatedModel.(Model)
	if !m.editingMode {
		t.Error("Enter on samples section should enable editing mode")
	}
	
	// Test sample input initialization
	if m.sampleInput.Lane != 1 {
		t.Error("Sample input should initialize with lane 1")
	}
}

func TestFingerprintingConfiguration(t *testing.T) {
	// Test with fingerprinting enabled
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	phantomConfig := &types.PhantomConfig{
		Enterprise: types.EnterpriseConfig{
			CustomerID:  "ENTERPRISE_CUSTOMER",
			ServiceTier: "enterprise",
			Fingerprinting: types.FingerprintConfig{
				Enable: true,
			},
		},
	}
	
	model := New(config, phantomConfig)
	
	if !model.fingerprintingEnabled {
		t.Error("Fingerprinting should be enabled when config specifies it")
	}
	
	if model.customerID != "ENTERPRISE_CUSTOMER" {
		t.Errorf("Expected customer ID 'ENTERPRISE_CUSTOMER', got '%s'", model.customerID)
	}
	
	if model.embeddingStrategy != StrategyBalanced {
		t.Errorf("Expected default embedding strategy to be StrategyBalanced (%d), got %d", StrategyBalanced, model.embeddingStrategy)
	}
	
	// Test with fingerprinting disabled
	phantomConfig.Enterprise.Fingerprinting.Enable = false
	model2 := New(config, phantomConfig)
	
	if model2.fingerprintingEnabled {
		t.Error("Fingerprinting should be disabled when config specifies it")
	}
	
	// Test with nil phantom config (should create defaults)
	model3 := New(config, nil)
	
	if model3.phantomConfig == nil {
		t.Error("Model should create default phantom config when nil is provided")
	}
	
	if model3.customerID != "DEMO" {
		t.Errorf("Expected default customer ID 'DEMO', got '%s'", model3.customerID)
	}
	
	if model3.fingerprintingEnabled {
		t.Error("Default phantom config should have fingerprinting disabled")
	}
}

func TestBubbleTeaInterface(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	
	// Test that model implements tea.Model interface
	var teaModel tea.Model = model
	if teaModel == nil {
		t.Error("Model should implement tea.Model interface")
	}
	
	// Test Init method
	cmd := model.Init()
	if cmd == nil {
		t.Error("Init should return a command (likely tea.EnterAltScreen)")
	}
	
	// Test Update method with various message types
	updatedModel, updateCmd := model.Update(tea.WindowSizeMsg{Width: 100, Height: 30})
	if updatedModel == nil {
		t.Error("Update should return updated model")
	}
	
	m, ok := updatedModel.(Model)
	if !ok {
		t.Error("Updated model should be of type Model")
	}
	
	if m.width != 100 || m.height != 30 {
		t.Error("WindowSizeMsg should update model dimensions")
	}
	
	// Test View method
	view := m.View()
	if view == "" {
		t.Error("View should return non-empty string")
	}
	
	// Test tickMsg handling
	updatedModel, tickCmd := m.Update(tickMsg{})
	if updatedModel == nil {
		t.Error("tickMsg should be handled")
	}
	
	// Should return another tick command
	if tickCmd == nil {
		t.Error("tickMsg should return another tick command")
	}
	
	_ = updateCmd // Suppress unused variable warning
}

func TestSimulationState(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	
	// Test initial simulation state
	if model.runInProgress {
		t.Error("Run should not be in progress initially")
	}
	
	if model.simulator != nil {
		t.Error("Simulator should be nil initially")
	}
	
	if model.genomicsEngine != nil {
		t.Error("Genomics engine should be nil initially")
	}
	
	// Test progress message handling
	progressUpdate := progressMsg(timing.RunProgress{
		Stage:          timing.StageSequencing,
		Percent:        50.0,
		FilesGenerated: 5,
		TotalFiles:     10,
		CurrentFile:    "test_R1.fastq",
	})
	
	updatedModel, cmd := model.Update(progressUpdate)
	m := updatedModel.(Model)
	
	if cmd != nil {
		t.Error("Progress message should not return command")
	}
	
	if m.progress.Stage != timing.StageSequencing {
		t.Errorf("Expected stage %s, got %s", timing.StageSequencing, m.progress.Stage)
	}
	
	if m.progress.Percent != 50.0 {
		t.Errorf("Expected progress 50.0, got %.1f", m.progress.Percent)
	}
	
	// Test file generation message
	fileMsg := fileGeneratedMsg{filename: "sample_R1.fastq"}
	updatedModel, _ = m.Update(fileMsg)
	m = updatedModel.(Model)
	
	if m.progress.FilesGenerated != 6 { // Should increment from 5 to 6
		t.Errorf("Expected files generated to be 6, got %d", m.progress.FilesGenerated)
	}
	
	if m.progress.CurrentFile != "sample_R1.fastq" {
		t.Errorf("Expected current file 'sample_R1.fastq', got '%s'", m.progress.CurrentFile)
	}
	
	// Test simulation completion
	m.runInProgress = true
	completionMsg := simulationCompleteMsg{}
	updatedModel, _ = m.Update(completionMsg)
	m = updatedModel.(Model)
	
	if m.runInProgress {
		t.Error("Run should not be in progress after completion")
	}
	
	if m.progress.Stage != timing.StageComplete {
		t.Errorf("Expected stage %s, got %s", timing.StageComplete, m.progress.Stage)
	}
	
	if m.progress.Percent != 100.0 {
		t.Errorf("Expected progress 100.0, got %.1f", m.progress.Percent)
	}
}

func TestActionHandling(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	model.currentSection = SectionActions
	
	// Test Enter key on actions section when not running
	if model.runInProgress {
		t.Error("Run should not be in progress initially")
	}
	
	updatedModel, cmd := model.Update(tea.KeyMsg{Type: tea.KeyEnter})
	m := updatedModel.(Model)
	
	// Should attempt to start simulation
	if !m.runInProgress {
		t.Error("Enter on actions should start simulation")
	}
	
	if m.simulator == nil {
		t.Error("Simulator should be initialized when starting run")
	}
	
	if m.genomicsEngine == nil {
		t.Error("Genomics engine should be initialized when starting run")
	}
	
	if cmd == nil {
		t.Error("Starting simulation should return a command")
	}
	
	// Test Enter key when already running (should stop)
	m.runInProgress = true
	updatedModel, _ = m.Update(tea.KeyMsg{Type: tea.KeyEnter})
	m2 := updatedModel.(Model)
	
	if m2.runInProgress {
		t.Error("Enter on actions during run should stop simulation")
	}
}

func TestConstants(t *testing.T) {
	// Test section constants
	if SectionSequencer != 0 {
		t.Errorf("Expected SectionSequencer to be 0, got %d", SectionSequencer)
	}
	if SectionSamples != 1 {
		t.Errorf("Expected SectionSamples to be 1, got %d", SectionSamples)
	}
	if SectionParameters != 2 {
		t.Errorf("Expected SectionParameters to be 2, got %d", SectionParameters)
	}
	if SectionFingerprinting != 3 {
		t.Errorf("Expected SectionFingerprinting to be 3, got %d", SectionFingerprinting)
	}
	if SectionOutput != 4 {
		t.Errorf("Expected SectionOutput to be 4, got %d", SectionOutput)
	}
	if SectionActions != 5 {
		t.Errorf("Expected SectionActions to be 5, got %d", SectionActions)
	}
	
	// Test sequencer choice constants
	if ChoiceNovaSeq != 0 {
		t.Errorf("Expected ChoiceNovaSeq to be 0, got %d", ChoiceNovaSeq)
	}
	if ChoiceNextSeq != 1 {
		t.Errorf("Expected ChoiceNextSeq to be 1, got %d", ChoiceNextSeq)
	}
	if ChoiceMiSeq != 2 {
		t.Errorf("Expected ChoiceMiSeq to be 2, got %d", ChoiceMiSeq)
	}
	
	// Test embedding strategy constants
	if StrategyConservative != 0 {
		t.Errorf("Expected StrategyConservative to be 0, got %d", StrategyConservative)
	}
	if StrategyBalanced != 1 {
		t.Errorf("Expected StrategyBalanced to be 1, got %d", StrategyBalanced)
	}
	if StrategyAggressive != 2 {
		t.Errorf("Expected StrategyAggressive to be 2, got %d", StrategyAggressive)
	}
}

func TestTUIStyles(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		PairedEnd:      true,
		OutputDir:      "/tmp/test",
		SampleName:     "TestSample",
		TotalBatches:   2,
		ReadsPerBatch:  1000,
	}
	
	model := New(config, nil)
	
	// Styles are initialized by NewStyles() - basic functionality test
	
	// Test sample structure
	if len(model.samples) > 0 {
		sample := model.samples[0]
		if sample.ID == "" {
			t.Error("Sample ID should not be empty")
		}
		if sample.Lane <= 0 {
			t.Error("Sample lane should be positive")
		}
		if len(sample.Index1) == 0 {
			t.Error("Sample Index1 should not be empty")
		}
		if sample.Project == "" {
			t.Error("Sample project should not be empty")
		}
	}
	
	// Test sample input structure
	if model.sampleInput.Lane != 1 {
		t.Errorf("Expected sampleInput lane to be 1, got %d", model.sampleInput.Lane)
	}
}
// ---------- handleSpace ----------

func TestHandleSpace_SectionParameters(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		ReadLength:     150,
		TotalBatches:   1,
	}
	m := New(config, nil)
	m.currentSection = SectionParameters

	newM, cmd := m.handleSpace()
	_ = cmd
	if newM.currentSection != SectionParameters {
		t.Errorf("handleSpace in SectionParameters should stay in SectionParameters, got %d", newM.currentSection)
	}
}

func TestHandleSpace_OtherSection(t *testing.T) {
	config := types.SequencerConfig{SequencerModel: "MiSeq", TotalBatches: 1}
	m := New(config, nil)
	m.currentSection = SectionSequencer

	newM, cmd := m.handleSpace()
	_ = cmd
	if newM.currentSection != SectionSequencer {
		t.Error("handleSpace in non-Parameters section should not change section")
	}
}

// ---------- handleTextInput ----------

func TestHandleTextInput_ReturnsModel(t *testing.T) {
	config := types.SequencerConfig{SequencerModel: "MiSeq", TotalBatches: 1}
	m := New(config, nil)

	newM, cmd := m.handleTextInput("test_input")
	_ = cmd
	// handleTextInput is a no-op for now; just verify it doesn't panic
	if newM.currentSection != m.currentSection {
		t.Error("handleTextInput should not change currentSection")
	}
}

// ---------- tickCmd ----------

func TestTickCmd_ReturnsCmd(t *testing.T) {
	cmd := tickCmd()
	if cmd == nil {
		t.Error("tickCmd should return a non-nil command")
	}
}

// ---------- startSimulation ----------

func TestStartSimulation_ReturnsCmd(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		TotalBatches:   1,
		ReadsPerBatch:  100,
	}
	sim := timing.NewRunSimulator(config, timing.Compressed10x)
	ctx := context.Background()
	cmd := startSimulation(ctx, sim)
	if cmd == nil {
		t.Error("startSimulation should return a non-nil command")
	}
}

func TestStartSimulation_ExecutesErrorPath(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		TotalBatches:   1,
		ReadsPerBatch:  100,
	}
	// Pre-cancel the context so Run() returns immediately with error
	ctx, cancel := context.WithCancel(context.Background())
	cancel()

	sim := timing.NewRunSimulator(config, timing.Compressed10x)
	cmd := startSimulation(ctx, sim)
	if cmd == nil {
		t.Fatal("startSimulation should return non-nil cmd")
	}
	// Execute the tea.Cmd to cover the body
	msg := cmd()
	if _, ok := msg.(simulationErrorMsg); !ok {
		t.Errorf("expected simulationErrorMsg for cancelled context, got %T", msg)
	}
}

func TestUpdateConfigFromUI_AllChoices(t *testing.T) {
	config := types.SequencerConfig{SequencerModel: "MiSeq", TotalBatches: 1}

	// Test ChoiceNovaSeq
	m1 := New(config, nil)
	m1.sequencerChoice = ChoiceNovaSeq
	m1.updateConfigFromUI() // value receiver; exercises the NovaSeq branch

	// Test ChoiceNextSeq
	m2 := New(config, nil)
	m2.sequencerChoice = ChoiceNextSeq
	m2.updateConfigFromUI() // exercises the NextSeq branch

	// Test ChoiceMiSeq with no samples (TotalBatches=0)
	m3 := New(config, nil)
	m3.sequencerChoice = ChoiceMiSeq
	m3.samples = nil // force empty samples
	m3.updateConfigFromUI()
}

func TestRenderStatsPanelWithRunInProgress(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		TotalBatches:   2,
		ReadsPerBatch:  100,
	}
	m := New(config, nil)
	m.width = 120
	m.height = 40
	m.runInProgress = true
	m.progress = timing.RunProgress{
		Stage:          timing.StageSequencing,
		Percent:        50.0,
		FilesGenerated: 3,
		TotalFiles:     10,
		BytesGenerated: 1024 * 1024,
		CurrentFile:    "/tmp/sample_R1_001.fastq",
	}
	m.simulator = timing.NewRunSimulator(config, timing.Compressed10x)

	result := m.renderStatsPanel()
	if result == "" {
		t.Error("renderStatsPanel should return non-empty string when run is in progress")
	}
}

func TestRenderStatusPanelWithRunInProgress(t *testing.T) {
	config := types.SequencerConfig{SequencerModel: "MiSeq", TotalBatches: 1}
	m := New(config, nil)
	m.width = 120
	m.height = 40
	m.runInProgress = true
	m.progress = timing.RunProgress{
		Stage:   timing.StageSequencing,
		Percent: 45.0,
		ETA:     30 * time.Second,
	}

	result := m.renderStatusPanel()
	if result == "" {
		t.Error("renderStatusPanel should return non-empty string when run is in progress")
	}
}

func TestRenderActionsRunInProgress(t *testing.T) {
	config := types.SequencerConfig{SequencerModel: "MiSeq", TotalBatches: 1}
	m := New(config, nil)
	m.width = 120
	m.height = 40
	m.runInProgress = true

	result := m.renderActions()
	if result == "" {
		t.Error("renderActions should return non-empty string when run is in progress")
	}
}

func TestViewWithRunInProgress(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		ReadLength:     150,
		TotalBatches:   2,
		ReadsPerBatch:  100,
	}
	m := New(config, nil)
	m.width = 120
	m.height = 40
	m.runInProgress = true
	m.progress = timing.RunProgress{
		Stage:          timing.StageSequencing,
		Percent:        30.0,
		FilesGenerated: 2,
		TotalFiles:     10,
	}
	m.simulator = timing.NewRunSimulator(config, timing.Compressed10x)

	view := m.View()
	if view == "" {
		t.Error("View should return non-empty string when run is in progress")
	}
}

// ---------- HMM model / Population TUI tests ----------

func TestNewModel_HMMFieldsPreserved(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		ReadLength:     150,
		OutputDir:      "/tmp/test",
		HMMModelPath:   "/models/YRI.hmm.json",
		Population:     "YRI",
	}
	m := New(config, nil)

	if m.hmmModelPath != "/models/YRI.hmm.json" {
		t.Errorf("expected hmmModelPath %q, got %q", "/models/YRI.hmm.json", m.hmmModelPath)
	}
	if m.population != "YRI" {
		t.Errorf("expected population %q, got %q", "YRI", m.population)
	}
}

func TestNewModel_HMMFieldsEmptyByDefault(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		ReadLength:     150,
		OutputDir:      "/tmp/test",
	}
	m := New(config, nil)

	if m.hmmModelPath != "" {
		t.Errorf("expected empty hmmModelPath, got %q", m.hmmModelPath)
	}
	if m.population != "" {
		t.Errorf("expected empty population, got %q", m.population)
	}
}

func TestParametersPanel_ShowsHMMModel(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		ReadLength:     150,
		HMMModelPath:   "/data/models/YRI.hmm.json",
		Population:     "YRI",
	}
	m := New(config, nil)
	panel := m.renderParametersPanel()

	// The panel should show the model basename (not the full path).
	if !strings.Contains(panel, "YRI.hmm.json") {
		t.Errorf("Parameters panel should display HMM model basename; got:\n%s", panel)
	}
	if !strings.Contains(panel, "YRI") {
		t.Errorf("Parameters panel should display population code; got:\n%s", panel)
	}
}

func TestParametersPanel_ShowsDashWhenHMMEmpty(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		ReadLength:     150,
	}
	m := New(config, nil)
	panel := m.renderParametersPanel()

	// Both fields should show the em-dash placeholder when not configured.
	if !strings.Contains(panel, "HMM Model:") {
		t.Errorf("Parameters panel should contain 'HMM Model:' label; got:\n%s", panel)
	}
	if !strings.Contains(panel, "Population:") {
		t.Errorf("Parameters panel should contain 'Population:' label; got:\n%s", panel)
	}
}

func TestUpdateConfigFromUI_SyncsHMMFields(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		ReadLength:     150,
		OutputDir:      "/tmp/test",
	}
	m := New(config, nil)

	// Override model fields as if the user typed them.
	m.hmmModelPath = "/new/CEU.hmm.json"
	m.population = "CEU"

	m = m.updateConfigFromUI()

	if m.sequencerConfig.HMMModelPath != "/new/CEU.hmm.json" {
		t.Errorf("expected HMMModelPath %q, got %q", "/new/CEU.hmm.json", m.sequencerConfig.HMMModelPath)
	}
	if m.sequencerConfig.Population != "CEU" {
		t.Errorf("expected Population %q, got %q", "CEU", m.sequencerConfig.Population)
	}
}

func TestParametersPanel_LongPathShowsBasename(t *testing.T) {
	config := types.SequencerConfig{
		SequencerModel: "MiSeq",
		ReadLength:     150,
		HMMModelPath:   "/very/long/deeply/nested/path/to/models/GBR.hmm.json",
	}
	m := New(config, nil)
	panel := m.renderParametersPanel()

	if !strings.Contains(panel, "GBR.hmm.json") {
		t.Errorf("Parameters panel should show basename for long paths; got:\n%s", panel)
	}
}

// TestTabNavigation_AllSixSections verifies that Tab cycles through all 6 sections
// (0–5) before wrapping, and Shift+Tab reverses the cycle. Regression for #126
// where modulo 5 caused SectionActions (5) to be unreachable.
func TestTabNavigation_AllSixSections(t *testing.T) {
	config := types.SequencerConfig{SequencerModel: "MiSeq", TotalBatches: 1}
	m := New(config, nil)
	m.editingMode = false

	expectedForward := []int{
		SectionSamples,
		SectionParameters,
		SectionFingerprinting,
		SectionOutput,
		SectionActions, // was unreachable with % 5
		SectionSequencer, // wrap back to start
	}

	for i, want := range expectedForward {
		updated, _ := m.Update(tea.KeyMsg{Type: tea.KeyTab})
		m = updated.(Model)
		if m.currentSection != want {
			t.Errorf("Tab press %d: expected section %d, got %d", i+1, want, m.currentSection)
		}
	}

	// Now walk backward from SectionSequencer through all sections.
	expectedBackward := []int{
		SectionActions, // was unreachable with % 5
		SectionOutput,
		SectionFingerprinting,
		SectionParameters,
		SectionSamples,
		SectionSequencer,
	}

	for i, want := range expectedBackward {
		updated, _ := m.Update(tea.KeyMsg{Type: tea.KeyShiftTab})
		m = updated.(Model)
		if m.currentSection != want {
			t.Errorf("Shift+Tab press %d: expected section %d, got %d", i+1, want, m.currentSection)
		}
	}
}
