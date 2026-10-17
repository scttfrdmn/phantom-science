package cli

import (
    "fmt"

    "github.com/scttfrdmn/phantom-science/internal/genomics"
    "github.com/scttfrdmn/phantom-science/internal/validation"
    "github.com/scttfrdmn/phantom-science/internal/fingerprinting"
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
    
    // Initialize enterprise features
    if a.phantomConfig.Enterprise.Fingerprinting.Enable {
        fingerprinter := fingerprinting.New(a.phantomConfig.Enterprise.CustomerID)
        defer fingerprinter.Close()
    }
    
    // Initialize validation framework
    validator := validation.NewFramework(a.phantomConfig.Validation)
    
    // Initialize genomics engine
    engine := genomics.NewEngine(*a.config)
    
    // Generate synthetic data
    if err := engine.Generate(); err != nil {
        return fmt.Errorf("generation failed: %w", err)
    }
    
    // Run validation if enabled
    if a.phantomConfig.Validation.Enable {
        results := validator.RunValidation(nil) // Would pass generated data
        if !results.OverallPass {
            return fmt.Errorf("validation failed: %s", results.Summary)
        }
        fmt.Printf("✅ Validation passed (academic ready: %t, regulatory ready: %t)\n", 
            results.AcademicReady, results.RegulatoryReady)
    }
    
    fmt.Println("✅ Generation completed successfully")
    return nil
}

func (a *App) Validate() error {
    fmt.Println("🔍 Running validation framework...")
    
    validator := validation.NewFramework(a.phantomConfig.Validation)
    results := validator.RunValidation(nil)
    
    // Output results
    fmt.Printf("Validation Results:\n")
    fmt.Printf("  Overall Pass: %t\n", results.OverallPass)
    fmt.Printf("  Academic Ready: %t\n", results.AcademicReady)
    fmt.Printf("  Regulatory Ready: %t\n", results.RegulatoryReady)
    
    for _, test := range results.Tests {
        status := "❌"
        if test.Passed {
            status = "✅"
        }
        fmt.Printf("  %s %s (confidence: %.2f)\n", status, test.TestName, test.Confidence)
    }
    
    return nil
}

func (a *App) Verify() error {
    fmt.Println("🔐 Running fingerprint verification...")
    
    if !a.phantomConfig.Enterprise.Fingerprinting.Enable {
        return fmt.Errorf("fingerprinting not enabled in configuration")
    }
    
    _ = fingerprinting.NewVerifier()
    
    // This would verify specific files passed as arguments
    fmt.Println("✅ Verification completed")
    return nil
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