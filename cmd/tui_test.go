package cmd

import (
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func TestTUICommand(t *testing.T) {
	// Test that TUI command is properly initialized
	if tuiCmd == nil {
		t.Fatal("TUI command should not be nil")
	}

	if tuiCmd.Use != "tui" {
		t.Errorf("Expected Use to be 'tui', got %s", tuiCmd.Use)
	}

	if tuiCmd.Short != "Launch Terminal User Interface for sequencer simulation" {
		t.Errorf("Expected Short description about TUI, got %s", tuiCmd.Short)
	}

	// Test that TUI command has a run function
	if tuiCmd.Run == nil {
		t.Error("TUI command should have a Run function")
	}

	// Test command name validation
	if len(tuiCmd.Use) == 0 {
		t.Error("TUI command should have a non-empty use name")
	}

	// Test that it's a valid command name (no spaces)
	if strings.Contains(tuiCmd.Use, " ") {
		t.Error("TUI command use should not contain spaces")
	}

	// Test aliases (if any)
	if len(tuiCmd.Aliases) > 0 {
		for _, alias := range tuiCmd.Aliases {
			if alias == "" {
				t.Error("TUI command aliases should not be empty")
			}
		}
	}
}

func TestGetDefaultSequencerConfig(t *testing.T) {
	// Test the default sequencer config generation
	config := getDefaultSequencerConfig()

	if config.SequencerModel == "" {
		t.Error("Default config should have a sequencer model")
	}

	if config.SampleName == "" {
		t.Error("Default config should have a sample name")
	}

	if config.ReadLength <= 0 {
		t.Error("Default config should have positive read length")
	}

	if config.ReadsPerBatch <= 0 {
		t.Error("Default config should have positive reads per batch")
	}

	if config.TotalBatches <= 0 {
		t.Error("Default config should have positive total batches")
	}

	if config.BaseQuality <= 0 {
		t.Error("Default config should have positive base quality")
	}

	if config.OutputDir == "" {
		t.Error("Default config should have output directory")
	}

	// Test that the config is suitable for TUI
	if config.FlowCellType == "" {
		t.Error("Default config should have flow cell type for TUI")
	}

	if config.Lane <= 0 {
		t.Error("Default config should have valid lane number")
	}

	if config.BatchInterval == "" {
		t.Error("Default config should have batch interval")
	}

	// Verify this produces a usable configuration
	if !config.PairedEnd {
		// Paired end should be true for comprehensive testing
		// This is not an error, just a note about the default
		t.Log("Default config has paired end disabled - single end mode")
	}
}

func TestGetDefaultSequencerConfigOutputDir(t *testing.T) {
	// Test that default config has the expected output directory
	config := getDefaultSequencerConfig()

	expectedPath := "./sequencer_output"
	if config.OutputDir != expectedPath {
		t.Errorf("Expected output dir %s, got %s", expectedPath, config.OutputDir)
	}
}

func TestGetDefaultSequencerConfigComplete(t *testing.T) {
	// Comprehensive test of default sequencer configuration
	config := getDefaultSequencerConfig()

	// Test specific values for quality assurance
	if config.SequencerModel != "MiSeq" {
		t.Errorf("Expected default sequencer model 'MiSeq', got %s", config.SequencerModel)
	}

	if config.FlowCellType != "v3" {
		t.Errorf("Expected default flow cell type 'v3', got %s", config.FlowCellType)
	}

	if config.SampleName != "Demo_Sample" {
		t.Errorf("Expected default sample name 'Demo_Sample', got %s", config.SampleName)
	}

	if config.Lane != 1 {
		t.Errorf("Expected default lane 1, got %d", config.Lane)
	}

	if config.ReadLength != 150 {
		t.Errorf("Expected default read length 150, got %d", config.ReadLength)
	}

	if config.ReadsPerBatch != 1000 {
		t.Errorf("Expected default reads per batch 1000, got %d", config.ReadsPerBatch)
	}

	if config.BatchInterval != "5m" {
		t.Errorf("Expected default batch interval '5m', got %s", config.BatchInterval)
	}

	if config.TotalBatches != 10 {
		t.Errorf("Expected default total batches 10, got %d", config.TotalBatches)
	}

	if !config.PairedEnd {
		t.Error("Expected default paired end to be true")
	}

	if config.BaseQuality != 30 {
		t.Errorf("Expected default base quality 30, got %d", config.BaseQuality)
	}

	if config.QualityDecay != 0.1 {
		t.Errorf("Expected default quality decay 0.1, got %f", config.QualityDecay)
	}

	if config.Q30Percentage != 85.0 {
		t.Errorf("Expected default Q30 percentage 85.0, got %f", config.Q30Percentage)
	}
}

func TestTUICommandDescription(t *testing.T) {
	// Test TUI command descriptions and documentation
	if tuiCmd.Long == "" {
		t.Error("TUI command should have a long description")
	}

	if !strings.Contains(tuiCmd.Long, "Terminal User Interface") {
		t.Error("TUI description should mention Terminal User Interface")
	}

	if !strings.Contains(tuiCmd.Long, "NovaSeq") {
		t.Error("TUI description should mention NovaSeq sequencer")
	}

	if !strings.Contains(tuiCmd.Long, "NextSeq") {
		t.Error("TUI description should mention NextSeq sequencer")
	}

	if !strings.Contains(tuiCmd.Long, "MiSeq") {
		t.Error("TUI description should mention MiSeq sequencer")
	}

	if !strings.Contains(tuiCmd.Long, "sample sheet") {
		t.Error("TUI description should mention the sample sheet")
	}

	if !strings.Contains(tuiCmd.Long, "time compression") {
		t.Error("TUI description should mention time compression")
	}

	// Test examples in description
	if !strings.Contains(tuiCmd.Long, "phantom-science tui") {
		t.Error("TUI description should contain usage examples")
	}
}

func TestTUICommandParent(t *testing.T) {
	// Test that TUI command is properly attached to root
	if tuiCmd.Parent() != rootCmd {
		t.Error("TUI command should have root command as parent")
	}

	// Test that command can be found
	found := false
	for _, cmd := range rootCmd.Commands() {
		if cmd == tuiCmd {
			found = true
			break
		}
	}
	if !found {
		t.Error("TUI command should be registered with root command")
	}
}

func TestTUICommandExecution(t *testing.T) {
	// Test TUI command execution setup (without actually running TUI)
	if tuiCmd.Run == nil {
		t.Fatal("TUI command should have a Run function")
	}

	if !tuiCmd.Runnable() {
		t.Error("TUI command should be runnable")
	}

	// Test argument validation
	err := tuiCmd.ValidateArgs([]string{})
	if err != nil {
		t.Errorf("TUI command should accept empty args, got error: %v", err)
	}

	err = tuiCmd.ValidateArgs([]string{"extra", "args"})
	if err != nil {
		t.Errorf("TUI command should accept extra args, got error: %v", err)
	}
}

func TestDefaultSequencerConfigTypes(t *testing.T) {
	// Test that default config returns proper types.SequencerConfig
	config := getDefaultSequencerConfig()

	// Verify it's the correct type
	var _ types.SequencerConfig = config

	// Test that all required fields are set for a valid configuration
	if config.SequencerModel == "" {
		t.Error("SequencerModel should not be empty")
	}

	if config.FlowCellType == "" {
		t.Error("FlowCellType should not be empty")
	}

	if config.SampleName == "" {
		t.Error("SampleName should not be empty")
	}

	if config.OutputDir == "" {
		t.Error("OutputDir should not be empty")
	}

	if config.BatchInterval == "" {
		t.Error("BatchInterval should not be empty")
	}

	// Test numerical constraints
	if config.Lane < 1 {
		t.Error("Lane should be >= 1")
	}

	if config.ReadLength < 1 {
		t.Error("ReadLength should be >= 1")
	}

	if config.ReadsPerBatch < 1 {
		t.Error("ReadsPerBatch should be >= 1")
	}

	if config.TotalBatches < 1 {
		t.Error("TotalBatches should be >= 1")
	}

	if config.BaseQuality < 1 {
		t.Error("BaseQuality should be >= 1")
	}

	if config.QualityDecay < 0 {
		t.Error("QualityDecay should be >= 0")
	}

	if config.Q30Percentage < 0 || config.Q30Percentage > 100 {
		t.Error("Q30Percentage should be between 0 and 100")
	}
}

func TestTUICommandUsage(t *testing.T) {
	// Test TUI command usage information
	usage := tuiCmd.UsageString()
	if usage == "" {
		t.Error("TUI command should have usage information")
	}

	if !strings.Contains(usage, "tui") {
		t.Error("Usage should contain command name 'tui'")
	}

	// Test that help includes flag descriptions
	if !strings.Contains(usage, "config") {
		t.Error("Usage should mention config flag")
	}

	if !strings.Contains(usage, "phantom-config") {
		t.Error("Usage should mention phantom-config flag")
	}

	if !strings.Contains(usage, "enterprise") {
		t.Error("Usage should mention enterprise flag")
	}
}

func TestTUICommandFlags(t *testing.T) {
	// Test that TUI command inherits flags from parent
	parentFlags := tuiCmd.InheritedFlags()

	// Should inherit verbose flag
	verboseFlag := parentFlags.Lookup("verbose")
	if verboseFlag == nil {
		t.Error("TUI command should inherit --verbose flag")
	}

	// Should inherit phantom flag (config is redefined locally)
	phantomFlag := parentFlags.Lookup("phantom")
	if phantomFlag == nil {
		t.Error("TUI command should inherit --phantom flag")
	}

	// Test TUI-specific flags
	tuiConfigFlag := tuiCmd.Flags().Lookup("config")
	if tuiConfigFlag == nil {
		t.Error("TUI command should have --config flag")
	}

	phantomConfigFlag := tuiCmd.Flags().Lookup("phantom-config")
	if phantomConfigFlag == nil {
		t.Error("TUI command should have --phantom-config flag")
	}

	enterpriseFlag := tuiCmd.Flags().Lookup("enterprise")
	if enterpriseFlag == nil {
		t.Error("TUI command should have --enterprise flag")
	}

	// Test flag defaults
	if tuiConfigFlag.DefValue != "" {
		t.Error("TUI config flag should have empty default")
	}

	if phantomConfigFlag.DefValue != "" {
		t.Error("TUI phantom-config flag should have empty default")
	}

	if enterpriseFlag.DefValue != "false" {
		t.Error("TUI enterprise flag should default to false")
	}

	// Test flag usage descriptions
	if tuiConfigFlag.Usage == "" {
		t.Error("TUI config flag should have usage description")
	}

	if phantomConfigFlag.Usage == "" {
		t.Error("TUI phantom-config flag should have usage description")
	}

	if enterpriseFlag.Usage == "" {
		t.Error("TUI enterprise flag should have usage description")
	}
}

func TestGetDefaultSequencerConfigConsistency(t *testing.T) {
	// Test that getDefaultSequencerConfig returns consistent results
	config1 := getDefaultSequencerConfig()
	config2 := getDefaultSequencerConfig()

	// Should return identical configurations
	if config1.SequencerModel != config2.SequencerModel {
		t.Error("getDefaultSequencerConfig should return consistent results")
	}

	if config1.FlowCellType != config2.FlowCellType {
		t.Error("getDefaultSequencerConfig should return consistent FlowCellType")
	}

	if config1.SampleName != config2.SampleName {
		t.Error("getDefaultSequencerConfig should return consistent SampleName")
	}

	if config1.ReadLength != config2.ReadLength {
		t.Error("getDefaultSequencerConfig should return consistent ReadLength")
	}

	if config1.OutputDir != config2.OutputDir {
		t.Error("getDefaultSequencerConfig should return consistent OutputDir")
	}
}

func TestTUICommandMetadata(t *testing.T) {
	// Test various command metadata
	if tuiCmd.Name() != "tui" {
		t.Errorf("Expected command name 'tui', got '%s'", tuiCmd.Name())
	}

	if tuiCmd.CommandPath() != "phantom-science tui" {
		t.Errorf("Expected command path 'phantom-science tui', got '%s'", tuiCmd.CommandPath())
	}

	// Test that command is not hidden
	if tuiCmd.Hidden {
		t.Error("TUI command should not be hidden")
	}

	// Test that command is not deprecated
	if tuiCmd.Deprecated != "" {
		t.Error("TUI command should not be deprecated")
	}
}
