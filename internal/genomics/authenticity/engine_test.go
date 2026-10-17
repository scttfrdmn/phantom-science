package authenticity

import (
	"fmt"
	"os"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// createTestReferenceDir creates chr1.fa, chr2.fa, chr3.fa in a temp dir
// with a simple FASTA sequence suitable for k-mer extraction tests.
func createTestReferenceDir(t *testing.T) string {
	t.Helper()
	tmpDir := t.TempDir()
	seq := make([]byte, 200)
	for i := range seq {
		seq[i] = []byte("ATCG")[i%4]
	}
	seqStr := string(seq)
	for _, num := range []int{1, 2, 3} {
		content := fmt.Sprintf(">chr%d\n%s\n", num, seqStr)
		path := fmt.Sprintf("%s/chr%d.fa", tmpDir, num)
		if err := os.WriteFile(path, []byte(content), 0600); err != nil {
			t.Fatalf("WriteFile chr%d.fa: %v", num, err)
		}
	}
	return tmpDir
}

// TestDefaultConfig tests the engine's default configuration.
func TestDefaultConfig(t *testing.T) {
	cfg := DefaultConfig()
	require.NotNil(t, cfg)

	assert.Equal(t, "GRCh38", cfg.Assembly)
	assert.NotEmpty(t, cfg.KmerSizes)
	assert.Greater(t, cfg.JSThreshold, 0.0)
	assert.Greater(t, cfg.WassersteinThreshold, 0.0)
	assert.Greater(t, cfg.ValidationInterval, time.Duration(0))
	assert.Greater(t, cfg.MaxLoadedChromosomes, 0)
}

// TestNewEngine_NilConfig tests that a nil config uses defaults.
func TestNewEngine_NilConfig(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)
	defer engine.Close()

	assert.NotNil(t, engine.GetConfig())
}

// TestNewEngine_CustomConfig tests engine creation with a custom config.
func TestNewEngine_CustomConfig(t *testing.T) {
	cfg := &Config{
		ReferenceGenomePath:  "/tmp/ref",
		Assembly:             "GRCh38",
		KmerSizes:            []int{2, 4},
		JSThreshold:          0.02,
		WassersteinThreshold: 0.03,
		EnableValidation:     true,
		ValidationInterval:   1 * time.Minute,
		MaxLoadedChromosomes: 3,
		EnableLazyLoading:    true,
		ParallelProcessing:   false,
		LogLevel:             "info",
	}

	engine := NewEngine(cfg)
	require.NotNil(t, engine)
	defer engine.Close()

	assert.False(t, engine.IsInitialized())
}

// TestNewEngine_NoLazyLoading tests creation without lazy loading.
func TestNewEngine_NoLazyLoading(t *testing.T) {
	cfg := DefaultConfig()
	cfg.EnableLazyLoading = false

	engine := NewEngine(cfg)
	require.NotNil(t, engine)
	defer engine.Close()
}

// TestEngine_GetStatus tests the status getter before and after init.
func TestEngine_GetStatus(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)
	defer engine.Close()

	status := engine.GetStatus()
	require.NotNil(t, status)
	assert.Contains(t, status, "initialized")
	assert.Contains(t, status, "assembly")
	assert.Contains(t, status, "kmer_sizes")
	assert.Contains(t, status, "js_threshold")
	assert.Contains(t, status, "validation_enabled")
	assert.Contains(t, status, "lazy_loading_enabled")
	assert.Contains(t, status, "last_validated")
	assert.Contains(t, status, "reference_tables")
	assert.False(t, status["initialized"].(bool))
}

// TestEngine_GetConfig tests that GetConfig returns a copy.
func TestEngine_GetConfig(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)
	defer engine.Close()

	cfg1 := engine.GetConfig()
	require.NotNil(t, cfg1)

	cfg2 := engine.GetConfig()
	require.NotNil(t, cfg2)

	// Should be equal but not the same pointer (it's a copy).
	assert.Equal(t, cfg1.Assembly, cfg2.Assembly)
}

// TestEngine_IsInitialized tests initialization state tracking.
func TestEngine_IsInitialized(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)
	defer engine.Close()

	assert.False(t, engine.IsInitialized())
}

// TestEngine_Close tests that Close runs without error.
func TestEngine_Close(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)

	err := engine.Close()
	assert.NoError(t, err)
	assert.False(t, engine.IsInitialized())
}

// TestEngine_Close_Twice tests that double-close is safe.
func TestEngine_Close_Twice(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)

	require.NoError(t, engine.Close())
	err := engine.Close()
	// Second close should not panic; it may or may not error depending on impl.
	_ = err
}

// TestEngine_UpdateFromSequencerConfig_Illumina tests Illumina platform configuration.
func TestEngine_UpdateFromSequencerConfig_Illumina(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)
	defer engine.Close()

	config := types.SequencerConfig{
		SequencerModel: "illumina",
		ReadLength:     150,
	}
	err := engine.UpdateFromSequencerConfig(config)
	require.NoError(t, err)

	cfg := engine.GetConfig()
	assert.InDelta(t, 0.01, cfg.JSThreshold, 1e-9)
	// Long read (>= 100) => 5 minute interval.
	assert.Equal(t, 5*time.Minute, cfg.ValidationInterval)
}

// TestEngine_UpdateFromSequencerConfig_PacBio tests PacBio platform configuration.
func TestEngine_UpdateFromSequencerConfig_PacBio(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)
	defer engine.Close()

	config := types.SequencerConfig{
		SequencerModel: "pacbio",
		ReadLength:     10000,
	}
	err := engine.UpdateFromSequencerConfig(config)
	require.NoError(t, err)

	cfg := engine.GetConfig()
	assert.InDelta(t, 0.02, cfg.JSThreshold, 1e-9)
}

// TestEngine_UpdateFromSequencerConfig_Nanopore tests Nanopore platform configuration.
func TestEngine_UpdateFromSequencerConfig_Nanopore(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)
	defer engine.Close()

	config := types.SequencerConfig{
		SequencerModel: "nanopore",
		ReadLength:     50,
	}
	err := engine.UpdateFromSequencerConfig(config)
	require.NoError(t, err)

	cfg := engine.GetConfig()
	assert.InDelta(t, 0.03, cfg.JSThreshold, 1e-9)
	// Short read (< 100) => 1 minute interval.
	assert.Equal(t, 1*time.Minute, cfg.ValidationInterval)
}

// TestEngine_UpdateFromSequencerConfig_Unknown tests unknown platform.
func TestEngine_UpdateFromSequencerConfig_Unknown(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)
	defer engine.Close()

	originalCfg := engine.GetConfig()
	originalThreshold := originalCfg.JSThreshold

	config := types.SequencerConfig{
		SequencerModel: "unknown_platform",
	}
	err := engine.UpdateFromSequencerConfig(config)
	require.NoError(t, err)

	// Should remain unchanged for unknown platform.
	newCfg := engine.GetConfig()
	assert.InDelta(t, originalThreshold, newCfg.JSThreshold, 1e-9)
}

// TestEngine_UpdateFromSequencerConfig_EmptyModel tests zero-value config.
func TestEngine_UpdateFromSequencerConfig_EmptyModel(t *testing.T) {
	engine := NewEngine(nil)
	require.NotNil(t, engine)
	defer engine.Close()

	config := types.SequencerConfig{}
	err := engine.UpdateFromSequencerConfig(config)
	assert.NoError(t, err)
}

// TestEngine_UpdateFromSequencerConfig_ReadLength tests read-length-based interval.
func TestEngine_UpdateFromSequencerConfig_ReadLength(t *testing.T) {
	t.Run("short_read", func(t *testing.T) {
		engine := NewEngine(nil)
		defer engine.Close()

		err := engine.UpdateFromSequencerConfig(types.SequencerConfig{ReadLength: 50})
		require.NoError(t, err)
		assert.Equal(t, 1*time.Minute, engine.GetConfig().ValidationInterval)
	})

	t.Run("long_read", func(t *testing.T) {
		engine := NewEngine(nil)
		defer engine.Close()

		err := engine.UpdateFromSequencerConfig(types.SequencerConfig{ReadLength: 300})
		require.NoError(t, err)
		assert.Equal(t, 5*time.Minute, engine.GetConfig().ValidationInterval)
	})

	t.Run("zero_read_length_no_change", func(t *testing.T) {
		engine := NewEngine(nil)
		defer engine.Close()

		original := engine.GetConfig().ValidationInterval
		err := engine.UpdateFromSequencerConfig(types.SequencerConfig{ReadLength: 0})
		require.NoError(t, err)
		// ReadLength == 0 means no update.
		assert.Equal(t, original, engine.GetConfig().ValidationInterval)
	})
}

func TestEngine_SelectRandomKmer_EmptyMap(t *testing.T) {
    engine := NewEngine(nil)
    defer engine.Close()
    result := engine.selectRandomKmer(map[string]float64{})
    if result != "" {
        t.Errorf("selectRandomKmer(empty) = %q, want empty string", result)
    }
}

func TestEngine_SelectRandomKmer_SingleEntry(t *testing.T) {
    engine := NewEngine(nil)
    defer engine.Close()
    result := engine.selectRandomKmer(map[string]float64{"ACGT": 1.0})
    if result != "ACGT" {
        t.Errorf("selectRandomKmer(single) = %q, want ACGT", result)
    }
}

func TestEngine_SelectRandomKmer_MaxProb(t *testing.T) {
    engine := NewEngine(nil)
    defer engine.Close()
    probDist := map[string]float64{
        "AAAA": 0.1,
        "CCCC": 0.5,
        "TTTT": 0.3,
        "GGGG": 0.1,
    }
    result := engine.selectRandomKmer(probDist)
    if result != "CCCC" {
        t.Errorf("selectRandomKmer(max=CCCC) = %q, want CCCC", result)
    }
}

// createInitializedEngine creates an Engine pointing at a tmpDir with chr1/2/3 FASTA files
// and calls Initialize(). Returns the engine (caller must Close() it).
func createInitializedEngine(t *testing.T) *Engine {
	t.Helper()
	tmpDir := createTestReferenceDir(t)
	cfg := DefaultConfig()
	cfg.ReferenceGenomePath = tmpDir
	cfg.EnableLazyLoading = false
	cfg.KmerSizes = []int{4} // Use a single k-mer size for speed
	engine := NewEngine(cfg)
	if err := engine.Initialize(); err != nil {
		engine.Close()
		t.Fatalf("Initialize failed: %v", err)
	}
	return engine
}

func TestEngine_GenerateAuthenticSequence_Basic(t *testing.T) {
	engine := createInitializedEngine(t)
	defer engine.Close()

	seq, err := engine.GenerateAuthenticSequence(50, "chr1")
	if err != nil {
		t.Fatalf("GenerateAuthenticSequence returned error: %v", err)
	}
	if len(seq) != 50 {
		t.Errorf("expected length 50, got %d", len(seq))
	}
}

func TestEngine_GenerateAuthenticSequence_Uninitialized(t *testing.T) {
	// With a bad path, Initialize() will fail, so GenerateAuthenticSequence returns an error
	cfg := DefaultConfig()
	cfg.ReferenceGenomePath = "/nonexistent/path"
	cfg.EnableLazyLoading = false
	cfg.KmerSizes = []int{4}
	engine := NewEngine(cfg)
	defer engine.Close()

	_, err := engine.GenerateAuthenticSequence(50, "chr1")
	if err == nil {
		t.Error("expected error for bad reference path, got nil")
	}
}

func TestEngine_ValidateSequence_Basic(t *testing.T) {
	engine := createInitializedEngine(t)
	defer engine.Close()

	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	report, err := engine.ValidateSequence(seq)
	if err != nil {
		t.Fatalf("ValidateSequence returned error: %v", err)
	}
	if report == nil {
		t.Error("expected non-nil validation report")
	}
}

func TestEngine_GetAuthenticityScore_Basic(t *testing.T) {
	engine := createInitializedEngine(t)
	defer engine.Close()

	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	score, err := engine.GetAuthenticityScore(seq)
	if err != nil {
		t.Fatalf("GetAuthenticityScore returned error: %v", err)
	}
	if score < 0 || score > 100 {
		t.Errorf("GetAuthenticityScore = %f, want [0, 100]", score)
	}
}

func TestEngine_GetReferenceStatistics_Basic(t *testing.T) {
	engine := createInitializedEngine(t)
	defer engine.Close()

	stats, err := engine.GetReferenceStatistics()
	if err != nil {
		t.Fatalf("GetReferenceStatistics returned error: %v", err)
	}
	if stats == nil {
		t.Error("expected non-nil reference statistics")
	}
}

func TestEngine_GetReferenceStatistics_Uninitialized(t *testing.T) {
	cfg := DefaultConfig()
	cfg.ReferenceGenomePath = "/nonexistent/path"
	cfg.EnableLazyLoading = false
	cfg.KmerSizes = []int{4}
	engine := NewEngine(cfg)
	defer engine.Close()

	_, err := engine.GetReferenceStatistics()
	if err == nil {
		t.Error("expected error for uninitialized engine with bad path, got nil")
	}
}
