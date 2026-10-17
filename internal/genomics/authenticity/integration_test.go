package authenticity

import (
	"fmt"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// mockGenomicsEngine implements the genomicsEngine interface for testing.
type mockGenomicsEngine struct {
	generateCalled        bool
	generateFASTQCalled   bool
	closeCalled           bool
	generateErr           error
	generateFASTQErr      error
}

func (m *mockGenomicsEngine) Generate() error {
	m.generateCalled = true
	return m.generateErr
}

func (m *mockGenomicsEngine) GenerateFASTQFile(filename string, batch, read int) error {
	m.generateFASTQCalled = true
	return m.generateFASTQErr
}

func (m *mockGenomicsEngine) Close() {
	m.closeCalled = true
}

// newTestIntegrationManager creates a minimal integration manager for testing.
// It uses the real authenticity engine (un-initialized) so we test the error paths too.
func newTestIntegrationManager() (*IntegrationManager, *Engine) {
	engine := NewEngine(nil)
	cfg := types.SequencerConfig{
		SampleName:    "TestSample",
		Lane:          1,
		ReadLength:    50,
		BaseQuality:   30,
		QualityDecay:  0.1,
		ReadsPerBatch: 2,
		OutputDir:     "/tmp",
	}
	manager := &IntegrationManager{
		authenticityEngine: engine,
		config:             cfg,
	}
	return manager, engine
}

// TestNewIntegrationManager tests the constructor.
func TestNewIntegrationManager(t *testing.T) {
	engine := NewEngine(nil)
	defer engine.Close()

	cfg := types.SequencerConfig{
		SampleName: "TestSample",
		OutputDir:  "/tmp",
	}

	manager := NewIntegrationManager(engine, nil, cfg)
	require.NotNil(t, manager)
}

// TestIntegrationManager_generateAuthenticQualityScores tests quality score generation.
func TestIntegrationManager_generateAuthenticQualityScores(t *testing.T) {
	manager, engine := newTestIntegrationManager()
	defer engine.Close()

	sequence := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")

	qualities, err := manager.generateAuthenticQualityScores(sequence)
	require.NoError(t, err)
	assert.Len(t, qualities, len(sequence))

	// All quality chars should be printable ASCII in the valid range.
	for _, q := range qualities {
		assert.GreaterOrEqual(t, q, byte('!'), "quality char too low: %c (%d)", q, q)
		assert.LessOrEqual(t, q, byte('J'), "quality char too high: %c (%d)", q, q)
	}
}

// TestIntegrationManager_applyBiologicalQualityPatterns tests the quality adjuster.
func TestIntegrationManager_applyBiologicalQualityPatterns(t *testing.T) {
	manager, engine := newTestIntegrationManager()
	defer engine.Close()

	// Build a varied sequence.
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")

	// Test at several positions.
	positions := []int{0, 5, 25, len(seq) - 1}
	for _, pos := range positions {
		if pos < len(seq) {
			quality := manager.applyBiologicalQualityPatterns(30, seq[pos], pos, seq)
			assert.GreaterOrEqual(t, quality, 0)
		}
	}
}

// TestIntegrationManager_applyBiologicalQualityPatterns_GCBias tests GC-biased quality.
func TestIntegrationManager_applyBiologicalQualityPatterns_GCBias(t *testing.T) {
	manager, engine := newTestIntegrationManager()
	defer engine.Close()

	// Create a GC-rich sequence.
	gcRichSeq := []byte("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")

	// Mid-sequence G should get reduced quality due to GC bias.
	qualityG := manager.applyBiologicalQualityPatterns(30, 'G', 20, gcRichSeq)
	// Quality should not exceed base value.
	assert.LessOrEqual(t, qualityG, 30)
}

// TestIntegrationManager_calculateLocalGCBias tests local GC content calculation.
func TestIntegrationManager_calculateLocalGCBias(t *testing.T) {
	manager, engine := newTestIntegrationManager()
	defer engine.Close()

	allGC := []byte("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")
	allAT := []byte("ATATATATATATATATATATATATATATATATATATATATATATATATAT")
	mixed := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")

	gcBiasGC := manager.calculateLocalGCBias(len(allGC)/2, allGC)
	gcBiasAT := manager.calculateLocalGCBias(len(allAT)/2, allAT)
	gcBiasMixed := manager.calculateLocalGCBias(len(mixed)/2, mixed)

	assert.Greater(t, gcBiasGC, 0.5, "GC-rich seq should have high GC bias")
	assert.Less(t, gcBiasAT, 0.1, "AT-rich seq should have low GC bias")
	assert.Greater(t, gcBiasMixed, 0.0)
	assert.Less(t, gcBiasMixed, 1.0)
}

// TestIntegrationManager_calculateLocalGCBias_EdgeCases tests boundary positions.
func TestIntegrationManager_calculateLocalGCBias_EdgeCases(t *testing.T) {
	manager, engine := newTestIntegrationManager()
	defer engine.Close()

	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")

	// Position 0 (beginning).
	bias0 := manager.calculateLocalGCBias(0, seq)
	assert.GreaterOrEqual(t, bias0, 0.0)
	assert.LessOrEqual(t, bias0, 1.0)

	// Position at end.
	biasEnd := manager.calculateLocalGCBias(len(seq)-1, seq)
	assert.GreaterOrEqual(t, biasEnd, 0.0)
	assert.LessOrEqual(t, biasEnd, 1.0)

	// Empty sequence (total==0 branch).
	biasEmpty := manager.calculateLocalGCBias(0, []byte{})
	assert.InDelta(t, 0.0, biasEmpty, 1e-9)
}

// TestIntegrationManager_calculateHomopolymerLength tests homopolymer run length.
func TestIntegrationManager_calculateHomopolymerLength(t *testing.T) {
	manager, engine := newTestIntegrationManager()
	defer engine.Close()

	tests := []struct {
		name     string
		seq      []byte
		pos      int
		expected int
	}{
		{
			name:     "no_run",
			seq:      []byte("ATCGATCG"),
			pos:      0, // A at position 0 - backward nothing, forward nothing
			expected: 1,
		},
		{
			name:     "run_of_4_As",
			seq:      []byte("AAAAAT"),
			pos:      2, // A in the middle of the AAAAA run
			expected: 5,
		},
		{
			name:     "position_beyond_end",
			seq:      []byte("AT"),
			pos:      10, // Beyond end
			expected: 0,
		},
		{
			name:     "single_char",
			seq:      []byte("A"),
			pos:      0,
			expected: 1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			length := manager.calculateHomopolymerLength(tt.pos, tt.seq)
			assert.Equal(t, tt.expected, length)
		})
	}
}

// TestNewEnhancedGenomicsEngine tests the enhanced engine constructor.
func TestNewEnhancedGenomicsEngine(t *testing.T) {
	authEngine := NewEngine(nil)
	defer authEngine.Close()

	mockEngine := &mockGenomicsEngine{}
	cfg := types.SequencerConfig{
		SampleName: "TestSample",
		OutputDir:  "/tmp",
	}

	enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)
	require.NotNil(t, enhanced)

	assert.True(t, enhanced.IsAuthenticityEnabled())
}

// TestEnhancedGenomicsEngine_SetAuthenticityEnabled tests the enabled toggle.
func TestEnhancedGenomicsEngine_SetAuthenticityEnabled(t *testing.T) {
	authEngine := NewEngine(nil)
	defer authEngine.Close()

	mockEngine := &mockGenomicsEngine{}
	cfg := types.SequencerConfig{OutputDir: "/tmp"}

	enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)
	require.NotNil(t, enhanced)

	enhanced.SetAuthenticityEnabled(false)
	assert.False(t, enhanced.IsAuthenticityEnabled())

	enhanced.SetAuthenticityEnabled(true)
	assert.True(t, enhanced.IsAuthenticityEnabled())
}

// TestEnhancedGenomicsEngine_GetAuthenticityStatus tests the status reporter.
func TestEnhancedGenomicsEngine_GetAuthenticityStatus(t *testing.T) {
	t.Run("with_engine", func(t *testing.T) {
		authEngine := NewEngine(nil)
		defer authEngine.Close()

		mockEngine := &mockGenomicsEngine{}
		cfg := types.SequencerConfig{OutputDir: "/tmp"}

		enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)
		status := enhanced.GetAuthenticityStatus()
		assert.True(t, status["available"].(bool))
		assert.True(t, status["enabled"].(bool))
	})

	t.Run("without_engine", func(t *testing.T) {
		enhanced := &EnhancedGenomicsEngine{
			originalEngine:     &mockGenomicsEngine{},
			authenticityEngine: nil,
			useAuthenticity:    false,
		}
		status := enhanced.GetAuthenticityStatus()
		assert.False(t, status["available"].(bool))
	})
}

// TestEnhancedGenomicsEngine_Generate_NoAuth tests generate with authenticity disabled.
func TestEnhancedGenomicsEngine_Generate_NoAuth(t *testing.T) {
	mockEngine := &mockGenomicsEngine{}

	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     mockEngine,
		authenticityEngine: nil,
		useAuthenticity:    false,
	}

	err := enhanced.Generate()
	require.NoError(t, err)
	assert.True(t, mockEngine.generateCalled)
}

// TestEnhancedGenomicsEngine_GenerateFASTQFile_NoAuth tests FASTQ generation without auth.
func TestEnhancedGenomicsEngine_GenerateFASTQFile_NoAuth(t *testing.T) {
	mockEngine := &mockGenomicsEngine{}

	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     mockEngine,
		authenticityEngine: nil,
		useAuthenticity:    false,
	}

	err := enhanced.GenerateFASTQFile("/tmp/test.fastq", 1, 1)
	require.NoError(t, err)
	assert.True(t, mockEngine.generateFASTQCalled)
}

// TestEnhancedGenomicsEngine_Close tests resource cleanup.
func TestEnhancedGenomicsEngine_Close(t *testing.T) {
	authEngine := NewEngine(nil)

	mockEngine := &mockGenomicsEngine{}
	cfg := types.SequencerConfig{OutputDir: "/tmp"}

	enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)
	require.NotNil(t, enhanced)

	enhanced.Close()
	assert.True(t, mockEngine.closeCalled)
}

// TestEnhancedGenomicsEngine_UpdateAuthenticityConfig tests config update.
func TestEnhancedGenomicsEngine_UpdateAuthenticityConfig(t *testing.T) {
	authEngine := NewEngine(nil)
	defer authEngine.Close()

	mockEngine := &mockGenomicsEngine{}
	cfg := types.SequencerConfig{OutputDir: "/tmp"}

	enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)

	err := enhanced.UpdateAuthenticityConfig(types.SequencerConfig{
		SequencerModel: "illumina",
		ReadLength:     100,
	})
	assert.NoError(t, err)
}

// TestEnhancedGenomicsEngine_UpdateAuthenticityConfig_NoEngine tests no-engine error.
func TestEnhancedGenomicsEngine_UpdateAuthenticityConfig_NoEngine(t *testing.T) {
	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     &mockGenomicsEngine{},
		authenticityEngine: nil,
	}

	err := enhanced.UpdateAuthenticityConfig(types.SequencerConfig{})
	assert.Error(t, err)
}

// TestEnhancedGenomicsEngine_ConfigureForPlatform tests platform configuration.
func TestEnhancedGenomicsEngine_ConfigureForPlatform(t *testing.T) {
	platforms := []string{"illumina", "pacbio", "nanopore"}
	for _, platform := range platforms {
		t.Run(platform, func(t *testing.T) {
			authEngine := NewEngine(nil)
			defer authEngine.Close()

			mockEngine := &mockGenomicsEngine{}
			cfg := types.SequencerConfig{OutputDir: "/tmp"}

			enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)
			err := enhanced.ConfigureForPlatform(platform)
			assert.NoError(t, err)
		})
	}
}

// TestEnhancedGenomicsEngine_ConfigureForPlatform_NoEngine tests no-engine error.
func TestEnhancedGenomicsEngine_ConfigureForPlatform_NoEngine(t *testing.T) {
	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     &mockGenomicsEngine{},
		authenticityEngine: nil,
	}

	err := enhanced.ConfigureForPlatform("illumina")
	assert.Error(t, err)
}

// TestEnhancedGenomicsEngine_GetReferenceStatistics_NoEngine tests no-engine error.
func TestEnhancedGenomicsEngine_GetReferenceStatistics_NoEngine(t *testing.T) {
	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     &mockGenomicsEngine{},
		authenticityEngine: nil,
		useAuthenticity:    false,
	}

	_, err := enhanced.GetReferenceStatistics()
	assert.Error(t, err)
}

// TestEnhancedGenomicsEngine_ValidateLastGeneratedSequence_NoEngine tests no-engine error.
func TestEnhancedGenomicsEngine_ValidateLastGeneratedSequence_NoEngine(t *testing.T) {
	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     &mockGenomicsEngine{},
		authenticityEngine: nil,
		useAuthenticity:    false,
	}

	_, err := enhanced.ValidateLastGeneratedSequence([]byte("ATCG"))
	assert.Error(t, err)
}

// TestEnhancedGenomicsEngine_GetAuthenticityScore_NoEngine tests no-engine error.
func TestEnhancedGenomicsEngine_GetAuthenticityScore_NoEngine(t *testing.T) {
	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     &mockGenomicsEngine{},
		authenticityEngine: nil,
		useAuthenticity:    false,
	}

	_, err := enhanced.GetAuthenticityScore([]byte("ATCG"))
	assert.Error(t, err)
}

// TestEnhancedGenomicsEngine_GetIntegrationInfo tests the integration info report.
func TestEnhancedGenomicsEngine_GetIntegrationInfo(t *testing.T) {
	authEngine := NewEngine(nil)
	defer authEngine.Close()

	mockEngine := &mockGenomicsEngine{}
	cfg := types.SequencerConfig{OutputDir: "/tmp"}

	enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)

	info := enhanced.GetIntegrationInfo()
	require.NotNil(t, info)

	assert.True(t, info["authenticity_available"].(bool))
	assert.True(t, info["authenticity_enabled"].(bool))
	assert.Contains(t, info, "supported_platforms")
	assert.Contains(t, info, "features")
}

// TestNewBackwardsCompatibilityWrapper tests the BCW constructor.
func TestNewBackwardsCompatibilityWrapper(t *testing.T) {
	authEngine := NewEngine(nil)
	defer authEngine.Close()

	mockEngine := &mockGenomicsEngine{}
	cfg := types.SequencerConfig{OutputDir: "/tmp"}

	enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)
	bcw := NewBackwardsCompatibilityWrapper(enhanced)
	require.NotNil(t, bcw)
}

// TestBackwardsCompatibilityWrapper_Generate tests the wrapped Generate call.
func TestBackwardsCompatibilityWrapper_Generate(t *testing.T) {
	mockEngine := &mockGenomicsEngine{}

	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     mockEngine,
		authenticityEngine: nil,
		useAuthenticity:    false,
	}

	bcw := NewBackwardsCompatibilityWrapper(enhanced)
	err := bcw.Generate()
	require.NoError(t, err)
	assert.True(t, mockEngine.generateCalled)
}

// TestBackwardsCompatibilityWrapper_GenerateFASTQFile tests wrapped FASTQ generation.
func TestBackwardsCompatibilityWrapper_GenerateFASTQFile(t *testing.T) {
	mockEngine := &mockGenomicsEngine{}

	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     mockEngine,
		authenticityEngine: nil,
		useAuthenticity:    false,
	}

	bcw := NewBackwardsCompatibilityWrapper(enhanced)
	err := bcw.GenerateFASTQFile("/tmp/test.fastq", 1, 1)
	require.NoError(t, err)
	assert.True(t, mockEngine.generateFASTQCalled)
}

// TestBackwardsCompatibilityWrapper_Close tests wrapped Close.
func TestBackwardsCompatibilityWrapper_Close(t *testing.T) {
	mockEngine := &mockGenomicsEngine{}

	enhanced := &EnhancedGenomicsEngine{
		originalEngine:     mockEngine,
		authenticityEngine: nil,
		useAuthenticity:    false,
	}

	bcw := NewBackwardsCompatibilityWrapper(enhanced)
	bcw.Close()
	assert.True(t, mockEngine.closeCalled)
}

// TestBackwardsCompatibilityWrapper_GetEnhancedFeatures tests feature access.
func TestBackwardsCompatibilityWrapper_GetEnhancedFeatures(t *testing.T) {
	authEngine := NewEngine(nil)
	defer authEngine.Close()

	mockEngine := &mockGenomicsEngine{}
	cfg := types.SequencerConfig{OutputDir: "/tmp"}

	enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)
	bcw := NewBackwardsCompatibilityWrapper(enhanced)

	features := bcw.GetEnhancedFeatures()
	assert.NotNil(t, features)
	assert.Equal(t, enhanced, features)
}

// TestEnhancedGenomicsEngine_Generate_WithAuth_InitFailure tests auth init failure fallback.
func TestEnhancedGenomicsEngine_Generate_WithAuth_InitFailure(t *testing.T) {
	// Use an uninitialized engine pointing to a non-existent path.
	// Initialize() will fail, and the engine should fall back to the original.
	authEngine := NewEngine(&Config{
		ReferenceGenomePath: "/nonexistent/path",
		Assembly:            "GRCh38",
		KmerSizes:           []int{4},
		JSThreshold:         0.01,
		ValidationInterval:  1 * time.Minute,
		MaxLoadedChromosomes: 3,
		EnableLazyLoading:   false, // direct load = will fail
	})
	defer authEngine.Close()

	mockEngine := &mockGenomicsEngine{}
	cfg := types.SequencerConfig{OutputDir: "/tmp"}

	enhanced := NewEnhancedGenomicsEngine(mockEngine, authEngine, cfg)
	require.NotNil(t, enhanced)
	assert.True(t, enhanced.IsAuthenticityEnabled())

	// Generate should fall back to original engine when authenticity init fails.
	err := enhanced.Generate()
	// Original mock doesn't error.
	require.NoError(t, err)
	assert.True(t, mockEngine.generateCalled)
}

// Ensure time is imported.
var _ = fmt.Sprintf
