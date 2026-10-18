package cmd

import (
	"bytes"
	"os"
	"strings"
	"testing"

	"github.com/spf13/cobra"
)

func TestRootCommand(t *testing.T) {
	// Test that root command is properly initialized
	if rootCmd == nil {
		t.Fatal("Root command should not be nil")
	}

	if rootCmd.Use != "phantom-science" {
		t.Errorf("Expected Use to be 'phantom-science', got %s", rootCmd.Use)
	}

	if rootCmd.Short != "Synthetic sequencing data generator" {
		t.Errorf("Expected Short description, got %s", rootCmd.Short)
	}

	// Test that command has expected persistent flags
	configFlag := rootCmd.PersistentFlags().Lookup("config")
	if configFlag == nil {
		t.Error("Expected --config flag to exist")
	}

	phantomFlag := rootCmd.PersistentFlags().Lookup("phantom")
	if phantomFlag == nil {
		t.Error("Expected --phantom flag to exist")
	}

	verboseFlag := rootCmd.PersistentFlags().Lookup("verbose")
	if verboseFlag == nil {
		t.Error("Expected --verbose flag to exist")
	}

	versionFlag := rootCmd.Flags().Lookup("version")
	if versionFlag == nil {
		t.Error("Expected --version flag to exist")
	}
}

func TestCommandValidation(t *testing.T) {
	// Test that root command has proper validation
	if rootCmd.Run == nil {
		t.Error("Root command should have a Run function")
	}

	// Test version constants
	if Version == "" {
		t.Error("Version should not be empty")
	}

	if BuildDate == "" {
		t.Error("BuildDate should not be empty")
	}
	
	// Test that constants have expected format
	if !strings.Contains(Version, ".") {
		t.Error("Version should contain version format with dots")
	}
	
	// Test Run function behavior (without actually running CLI)
	// We test the function exists and can be called safely
	if rootCmd.Runnable() != true {
		t.Error("Root command should be runnable")
	}
}

func TestExecuteFunction(t *testing.T) {
	// Test that Execute function exists by checking if rootCmd can execute
	// We test the setup rather than actually calling Execute()
	err := rootCmd.ValidateArgs([]string{})
	if err != nil {
		t.Errorf("Root command should accept empty args, got error: %v", err)
	}
}

func TestVersionHandling(t *testing.T) {
	// Test version flag functionality
	cmd := &cobra.Command{
		Use: "test",
		Run: func(cmd *cobra.Command, args []string) {},
	}
	cmd.Flags().Bool("version", false, "version flag")
	
	// Test that version constants are set
	if Version == "" {
		t.Error("Version constant should not be empty")
	}
	
	if BuildDate == "" {
		t.Error("BuildDate constant should not be empty")
	}
	
	// The VERSION file is the source of truth (make version, make
	// sync-version, scripts/check-version.sh); the compiled-in constant must
	// match it.
	raw, err := os.ReadFile("../VERSION")
	if err != nil {
		t.Fatalf("reading VERSION file: %v", err)
	}
	if want := strings.TrimSpace(string(raw)); Version != want {
		t.Errorf("Version = %q, but VERSION file says %q (run make sync-version)", Version, want)
	}
}

func TestRootCommandArgs(t *testing.T) {
	// Test different argument scenarios
	testCases := []struct {
		name string
		args []string
		should string
	}{
		{"empty args", []string{}, "accept empty args"},
		{"single arg", []string{"test"}, "accept single arg"},
		{"multiple args", []string{"test", "args"}, "accept multiple args"},
	}
	
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			err := rootCmd.ValidateArgs(tc.args)
			if err != nil {
				t.Errorf("Root command should %s, got error: %v", tc.should, err)
			}
		})
	}
}

func TestRootCommandStructure(t *testing.T) {
	// Test command structure and relationships
	if len(rootCmd.Commands()) == 0 {
		t.Error("Root command should have subcommands")
	}
	
	// Test that TUI command is added
	found := false
	for _, cmd := range rootCmd.Commands() {
		if cmd.Use == "tui" {
			found = true
			break
		}
	}
	if !found {
		t.Error("TUI command should be added to root command")
	}
	
	// Test command descriptions
	if rootCmd.Long == "" {
		t.Error("Root command should have a long description")
	}
	
	if !strings.Contains(rootCmd.Long, "fingerprint") {
		t.Error("Root command description should mention fingerprinting")
	}
	
	if !strings.Contains(rootCmd.Long, "FASTQ") {
		t.Error("Root command description should mention FASTQ files")
	}
	
	// Test root command metadata
	if rootCmd.HasSubCommands() == false {
		t.Error("Root command should have subcommands")
	}
	
	if rootCmd.IsAvailableCommand() == false {
		t.Error("Root command should be available")
	}
	
	// Test that root command is runnable
	if !rootCmd.Runnable() {
		t.Error("Root command should be runnable")
	}
}

func TestRootCommandFlags(t *testing.T) {
	// Test flag defaults and types
	configFlag := rootCmd.PersistentFlags().Lookup("config")
	if configFlag == nil {
		t.Fatal("Config flag should exist")
	}
	
	if configFlag.DefValue != "" {
		t.Error("Config flag should have empty default value")
	}
	
	if configFlag.Shorthand != "c" {
		t.Errorf("Config flag should have shorthand 'c', got '%s'", configFlag.Shorthand)
	}
	
	// Test phantom flag
	phantomFlag := rootCmd.PersistentFlags().Lookup("phantom")
	if phantomFlag == nil {
		t.Fatal("Phantom flag should exist")
	}
	
	if phantomFlag.Shorthand != "p" {
		t.Errorf("Phantom flag should have shorthand 'p', got '%s'", phantomFlag.Shorthand)
	}
	
	// Test verbose flag
	verboseFlag := rootCmd.PersistentFlags().Lookup("verbose")
	if verboseFlag == nil {
		t.Fatal("Verbose flag should exist")
	}
	
	if verboseFlag.Shorthand != "v" {
		t.Errorf("Verbose flag should have shorthand 'v', got '%s'", verboseFlag.Shorthand)
	}
	
	// Test command flag
	commandFlag := rootCmd.Flags().Lookup("command")
	if commandFlag == nil {
		t.Fatal("Command flag should exist")
	}
	
	if commandFlag.DefValue != "generate" {
		t.Errorf("Command flag should default to 'generate', got '%s'", commandFlag.DefValue)
	}
	
	// Test version flag
	versionFlag := rootCmd.Flags().Lookup("version")
	if versionFlag == nil {
		t.Fatal("Version flag should exist")
	}
	
	if versionFlag.DefValue != "false" {
		t.Error("Version flag should default to false")
	}
	
	// Test flag types
	if configFlag.Value.Type() != "string" {
		t.Error("Config flag should be string type")
	}
	
	if phantomFlag.Value.Type() != "string" {
		t.Error("Phantom flag should be string type")
	}
	
	if verboseFlag.Value.Type() != "bool" {
		t.Error("Verbose flag should be bool type")
	}
	
	if versionFlag.Value.Type() != "bool" {
		t.Error("Version flag should be bool type")
	}
	
	// Test flag usage descriptions
	if configFlag.Usage == "" {
		t.Error("Config flag should have usage description")
	}
	
	if phantomFlag.Usage == "" {
		t.Error("Phantom flag should have usage description")
	}
	
	if verboseFlag.Usage == "" {
		t.Error("Verbose flag should have usage description")
	}
}

func TestRootCommandPreRun(t *testing.T) {
	// Test PreRun function existence and basic structure
	if rootCmd.PreRun == nil {
		t.Error("Root command should have a PreRun function for version handling")
	}
	
	// Test that PreRun doesn't panic with mock command
	mockCmd := &cobra.Command{Use: "test"}
	mockCmd.Flags().Bool("version", false, "version flag")
	
	// This should not panic when version is false
	defer func() {
		if r := recover(); r != nil {
			t.Errorf("PreRun should not panic when version is false, got: %v", r)
		}
	}()
	
	// Call PreRun with version flag set to false (should not exit)
	rootCmd.PreRun(mockCmd, []string{})
	
	// Note: We can't test the actual version output without mocking os.Exit
	// but we can test that the function exists and handles the flag properly
}

func TestCommandUsage(t *testing.T) {
	// Test command usage information
	usage := rootCmd.UsageString()
	if usage == "" {
		t.Error("Root command should have usage information")
	}
	
	if !strings.Contains(usage, "phantom-science") {
		t.Error("Usage should contain command name")
	}
	
	if !strings.Contains(usage, "config") {
		t.Error("Usage should mention config flag")
	}
	
	// Test help text
	helpText := rootCmd.Long
	if !strings.Contains(helpText, "Copyright") {
		t.Error("Help text should contain copyright information")
	}
	
	// Test that help mentions key features
	for _, banned := range []string{"enterprise-grade", "production-ready", "Enterprise"} {
		if strings.Contains(helpText, banned) {
			t.Errorf("Help should not claim %q", banned)
		}
	}
	
	if !strings.Contains(helpText, "FASTQ") {
		t.Error("Help should mention FASTQ generation")
	}
	
	if !strings.Contains(helpText, "sequencer") {
		t.Error("Help should mention sequencer simulation")
	}
	
	// Test short description
	if rootCmd.Short == "" {
		t.Error("Root command should have short description")
	}
	
	// Test command name
	if rootCmd.Name() != "phantom-science" {
		t.Errorf("Expected root command name 'phantom-science', got '%s'", rootCmd.Name())
	}
}

func TestRootCommandOutput(t *testing.T) {
	// Test command output formatting
	// Capture command output for testing
	var buf bytes.Buffer
	rootCmd.SetOut(&buf)
	rootCmd.SetErr(&buf)
	
	// Test help output
	rootCmd.SetArgs([]string{"--help"})
	
	// Reset args after test
	defer rootCmd.SetArgs([]string{})
	
	// The help should be available
	helpCmd := rootCmd.HelpFunc()
	if helpCmd == nil {
		t.Error("Root command should have help function")
	}
	
	// Test usage template
	usageTemplate := rootCmd.UsageTemplate()
	if usageTemplate == "" {
		t.Error("Root command should have usage template")
	}
	
	// Test help template
	helpTemplate := rootCmd.HelpTemplate()
	if helpTemplate == "" {
		t.Error("Root command should have help template")
	}
}