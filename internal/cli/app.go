package cli

import (
    "fmt"
    "os"

    "github.com/scttfrdmn/phantom-science/internal/genomics"
    "github.com/scttfrdmn/phantom-science/pkg/types"
)

type App struct {
    config        *types.SequencerConfig
    phantomConfig *types.PhantomConfig
    verbose       bool
}

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
        if err := ValidateFiles(os.Stdout, engine.GeneratedFiles(), cfg); err != nil {
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
func (a *App) Verify(paths []string, expect VerifyExpectation) error {
    fmt.Println("🔐 Running fingerprint verification...")
    return VerifyFiles(os.Stdout, paths, expect)
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

// Legacy functions for backward compatibility - now delegated to ConfigLoader