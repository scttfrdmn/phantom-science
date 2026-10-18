// Package cli implements the phantom-science root command's actions: loading
// sequencer and phantom configs, generating FASTQ, validating files with
// package validation, and verifying their fingerprints.
package cli

import (
	"fmt"
	"os"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// App holds the loaded configs for one invocation of the root command.
type App struct {
	config        *types.SequencerConfig
	phantomConfig *types.PhantomConfig
	verbose       bool
}

// NewApp loads the sequencer config at configPath and, when phantomPath is
// set, the phantom config. Without a phantom config, fingerprinting and
// validation are disabled.
func NewApp(configPath string, phantomPath string, verbose bool) (*App, error) {
	loader := NewConfigLoader(verbose)

	config, err := loader.LoadSequencerConfig(configPath)
	if err != nil {
		return nil, fmt.Errorf("failed to load sequencer config: %w", err)
	}

	var phantomConfig *types.PhantomConfig
	if phantomPath != "" {
		phantomConfig, err = loader.LoadPhantomConfig(phantomPath)
		if err != nil {
			return nil, fmt.Errorf("failed to load phantom config: %w", err)
		}
	} else {
		// Use minimal defaults when no phantom config provided
		phantomConfig = &types.PhantomConfig{
			Enterprise: types.EnterpriseConfig{
				Fingerprinting: types.FingerprintConfig{Enable: false},
			},
			Validation: types.ValidationConfig{Enable: false},
		}
	}

	return &App{
		config:        config,
		phantomConfig: phantomConfig,
		verbose:       verbose,
	}, nil
}

// Generate writes the configured FASTQ output and, when the phantom config
// enables validation, validates the files it wrote.
func (a *App) Generate() error {
	fmt.Println("🧬 Starting Phantom Science genomic data generation...")

	// Initialize genomics engine with the effective fingerprinting settings
	config := *a.config
	config.ResolveFingerprinting(a.phantomConfig)
	engine := genomics.NewEngine(config)
	defer engine.Close()

	// Generate synthetic data
	if err := engine.Generate(); err != nil {
		return fmt.Errorf("generation failed: %w", err)
	}

	// Validate the generated files if the phantom config enables validation
	if cfg, ok := EffectiveValidation(a.phantomConfig.Validation); ok {
		lengths, err := ReadLengthValidation(config)
		if err != nil {
			return err
		}
		if err := ValidateFiles(os.Stdout, engine.GeneratedFiles(), cfg, lengths); err != nil {
			return fmt.Errorf("validation failed: %w", err)
		}
	}

	fmt.Println("✅ Generation completed successfully")
	return nil
}

// Validate runs the validation framework on each FASTQ file in paths, using
// the phantom config's validation tests (see EffectiveValidation) and
// DefaultValidationTests when the phantom config does not enable validation.
func (a *App) Validate(paths []string) error {
	fmt.Println("🔍 Running validation framework...")
	cfg := DefaultValidationTests
	if a.phantomConfig != nil {
		if effective, ok := EffectiveValidation(a.phantomConfig.Validation); ok {
			cfg = effective
		}
	}
	return ValidateFiles(os.Stdout, paths, cfg)
}

// Verify checks each FASTQ file in paths for a steganographic fingerprint.
func (a *App) Verify(paths []string, opts VerifyOptions) error {
	fmt.Println("🔐 Running fingerprint verification...")
	return VerifyFiles(os.Stdout, paths, opts)
}

// GetConfig returns the sequencer configuration
func (a *App) GetConfig() *types.SequencerConfig {
	return a.config
}

// GetPhantomConfig returns the phantom configuration
func (a *App) GetPhantomConfig() *types.PhantomConfig {
	return a.phantomConfig
}

// SetConfig sets the sequencer configuration
func (a *App) SetConfig(config types.SequencerConfig) {
	a.config = &config
}

// SetPhantomConfig sets the phantom configuration
func (a *App) SetPhantomConfig(config *types.PhantomConfig) {
	a.phantomConfig = config
}
