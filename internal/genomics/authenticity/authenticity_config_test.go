package authenticity

import (
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestDefaultAuthenticityConfig tests the default config constructor.
func TestDefaultAuthenticityConfig(t *testing.T) {
	cfg := DefaultAuthenticityConfig()
	require.NotNil(t, cfg)

	assert.True(t, cfg.Enabled)
	assert.Equal(t, "GRCh38", cfg.ReferenceGenome)
	assert.Equal(t, "EUR", cfg.Population)
	assert.NotEmpty(t, cfg.KmerModelSizes)
	assert.Greater(t, cfg.JSdivergenceThreshold, 0.0)
	assert.LessOrEqual(t, cfg.JSdivergenceThreshold, 1.0)
	assert.Greater(t, cfg.MaxMemoryUsage, int64(0))
	assert.Greater(t, cfg.NumWorkers, 0)
}

// TestEnhancedAuthenticityConfig tests the enhanced config constructor.
func TestEnhancedAuthenticityConfig(t *testing.T) {
	cfg := EnhancedAuthenticityConfig()
	require.NotNil(t, cfg)

	assert.Equal(t, "enhanced", cfg.QualityRealism)
	assert.Equal(t, "comprehensive", cfg.ValidationLevel)
	assert.True(t, cfg.EnablePopulationGenetics)
	assert.True(t, cfg.EnableGenomicFeatures)
	// Enhanced should have at least as much memory as default.
	defaultCfg := DefaultAuthenticityConfig()
	assert.GreaterOrEqual(t, cfg.MaxMemoryUsage, defaultCfg.MaxMemoryUsage)
}

// TestExpertAuthenticityConfig tests the expert config constructor.
func TestExpertAuthenticityConfig(t *testing.T) {
	cfg := ExpertAuthenticityConfig()
	require.NotNil(t, cfg)

	assert.Equal(t, "expert", cfg.QualityRealism)
	assert.Equal(t, "expert", cfg.ValidationLevel)
	assert.Contains(t, cfg.KmerModelSizes, 8)
	// Expert should be stricter than default.
	assert.LessOrEqual(t, cfg.JSdivergenceThreshold, 0.01)
}

// TestAuthenticityConfig_Validate_Valid tests valid configurations pass validation.
func TestAuthenticityConfig_Validate_Valid(t *testing.T) {
	tests := []struct {
		name  string
		setup func(*AuthenticityConfig)
	}{
		{
			name: "default_with_path",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
			},
		},
		{
			name: "GRCh39",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.ReferenceGenome = "GRCh39"
			},
		},
		{
			name: "AFR_population",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.Population = "AFR"
			},
		},
		{
			name: "EAS_population",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.Population = "EAS"
			},
		},
		{
			name: "SAS_population",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.Population = "SAS"
			},
		},
		{
			name: "AMR_population",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.Population = "AMR"
			},
		},
		{
			name: "expert_validation",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.ValidationLevel = "expert"
			},
		},
		{
			name: "comprehensive_validation",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.ValidationLevel = "comprehensive"
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			cfg := DefaultAuthenticityConfig()
			tt.setup(cfg)
			err := cfg.Validate()
			assert.NoError(t, err)
		})
	}
}

// TestAuthenticityConfig_Validate_Disabled tests that disabled config skips validation.
func TestAuthenticityConfig_Validate_Disabled(t *testing.T) {
	cfg := DefaultAuthenticityConfig()
	cfg.Enabled = false
	cfg.ReferenceGenomePath = "" // would be invalid if enabled
	err := cfg.Validate()
	assert.NoError(t, err, "disabled config should not require a path")
}

// TestAuthenticityConfig_Validate_Invalid tests various invalid configurations.
func TestAuthenticityConfig_Validate_Invalid(t *testing.T) {
	tests := []struct {
		name      string
		setup     func(*AuthenticityConfig)
		expectMsg string
	}{
		{
			name: "bad_reference_genome",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.ReferenceGenome = "hg19"
			},
			expectMsg: "unsupported reference genome",
		},
		{
			name: "missing_path",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = ""
			},
			expectMsg: "path must be specified",
		},
		{
			name: "bad_population",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.Population = "UNKNOWN"
			},
			expectMsg: "unsupported population",
		},
		{
			name: "empty_kmer_sizes",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.KmerModelSizes = []int{}
			},
			expectMsg: "at least one k-mer model size",
		},
		{
			name: "kmer_size_too_small",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.KmerModelSizes = []int{1}
			},
			expectMsg: "k-mer size must be between",
		},
		{
			name: "kmer_size_too_large",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.KmerModelSizes = []int{9}
			},
			expectMsg: "k-mer size must be between",
		},
		{
			name: "JS_threshold_zero",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.JSdivergenceThreshold = 0
			},
			expectMsg: "JS divergence threshold",
		},
		{
			name: "JS_threshold_too_large",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.JSdivergenceThreshold = 1.5
			},
			expectMsg: "JS divergence threshold",
		},
		{
			name: "bad_quality_realism",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.QualityRealism = "extreme"
			},
			expectMsg: "unsupported quality realism",
		},
		{
			name: "bad_validation_level",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.ValidationLevel = "legendary"
			},
			expectMsg: "unsupported validation level",
		},
		{
			name: "memory_too_low",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.MaxMemoryUsage = 1024 // 1 KB
			},
			expectMsg: "max memory usage too low",
		},
		{
			name: "workers_too_few",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.NumWorkers = 0
			},
			expectMsg: "number of workers",
		},
		{
			name: "workers_too_many",
			setup: func(c *AuthenticityConfig) {
				c.ReferenceGenomePath = "/tmp/ref.fa"
				c.NumWorkers = 65
			},
			expectMsg: "number of workers",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			cfg := DefaultAuthenticityConfig()
			tt.setup(cfg)
			err := cfg.Validate()
			require.Error(t, err)
			assert.True(t, strings.Contains(err.Error(), tt.expectMsg),
				"expected error to contain %q, got: %v", tt.expectMsg, err)
		})
	}
}

// TestAuthenticityConfig_GetPopulationConfig tests population config lookup for all populations.
func TestAuthenticityConfig_GetPopulationConfig(t *testing.T) {
	populations := []string{"EUR", "AFR", "EAS", "SAS", "AMR"}
	for _, pop := range populations {
		t.Run(pop, func(t *testing.T) {
			cfg := DefaultAuthenticityConfig()
			cfg.Population = pop
			popCfg := cfg.GetPopulationConfig()
			require.NotNil(t, popCfg)
			assert.NotEmpty(t, popCfg.PopulationName)
			assert.Greater(t, popCfg.EffectivePopSize, 0)
			assert.Greater(t, popCfg.MutationRate, 0.0)
			assert.Greater(t, popCfg.RecombinationRate, 0.0)
		})
	}

	// Unknown population should fall back to EUR.
	t.Run("unknown_population", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		cfg.Population = "UNKNOWN"
		popCfg := cfg.GetPopulationConfig()
		require.NotNil(t, popCfg)
		assert.Equal(t, "European", popCfg.PopulationName)
	})
}

// TestAuthenticityConfig_GetQualityConfig tests quality config lookup for known platforms.
func TestAuthenticityConfig_GetQualityConfig(t *testing.T) {
	platforms := []string{"illumina", "pacbio", "nanopore"}
	for _, platform := range platforms {
		t.Run(platform, func(t *testing.T) {
			cfg := DefaultAuthenticityConfig()
			qCfg := cfg.GetQualityConfig(platform)
			require.NotNil(t, qCfg)
			assert.Equal(t, platform, qCfg.Platform)
			assert.Greater(t, qCfg.BaseQuality, 0.0)
			assert.Greater(t, qCfg.ErrorRate, 0.0)
		})
	}

	// Unknown platform falls back to Illumina.
	t.Run("unknown_platform", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		qCfg := cfg.GetQualityConfig("unknown")
		require.NotNil(t, qCfg)
		assert.Equal(t, "illumina", qCfg.Platform)
	})
}

// TestAuthenticityConfig_GetValidationConfig tests validation config lookup for all levels.
func TestAuthenticityConfig_GetValidationConfig(t *testing.T) {
	levels := []string{"basic", "comprehensive", "expert"}
	for _, level := range levels {
		t.Run(level, func(t *testing.T) {
			cfg := DefaultAuthenticityConfig()
			cfg.ValidationLevel = level
			vCfg := cfg.GetValidationConfig()
			require.NotNil(t, vCfg)
			assert.True(t, vCfg.EnableJSdivergence)
			assert.Greater(t, vCfg.SampleSize, 0)
			assert.Greater(t, vCfg.ConfidenceLevel, 0.0)
			assert.LessOrEqual(t, vCfg.ConfidenceLevel, 1.0)
		})
	}

	// Unknown level falls back to basic.
	t.Run("unknown_level", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		cfg.ValidationLevel = "unknown"
		vCfg := cfg.GetValidationConfig()
		require.NotNil(t, vCfg)
	})
}

// TestAuthenticityConfig_CachePaths tests cache path generation helpers.
func TestAuthenticityConfig_CachePaths(t *testing.T) {
	cfg := DefaultAuthenticityConfig()
	cfg.CacheDirectory = "/tmp/cache"
	cfg.ReferenceGenome = "GRCh38"
	cfg.Population = "EUR"

	t.Run("GetCacheFilePath", func(t *testing.T) {
		path := cfg.GetCacheFilePath("test.cache")
		assert.Contains(t, path, "/tmp/cache")
		assert.Contains(t, path, "test.cache")
	})

	t.Run("GetKmerCachePath", func(t *testing.T) {
		for _, k := range []int{2, 4, 6, 8} {
			path := cfg.GetKmerCachePath(k)
			assert.Contains(t, path, "/tmp/cache")
			assert.Contains(t, path, "GRCh38")
			assert.Contains(t, path, "EUR")
		}
	})

	t.Run("GetModelCachePath", func(t *testing.T) {
		path := cfg.GetModelCachePath("hmm")
		assert.Contains(t, path, "/tmp/cache")
		assert.Contains(t, path, "hmm")
		assert.Contains(t, path, "GRCh38")
	})
}

// TestAuthenticityConfig_IsHighMemoryMode tests high-memory mode detection.
func TestAuthenticityConfig_IsHighMemoryMode(t *testing.T) {
	t.Run("not_high_memory", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		cfg.MaxMemoryUsage = 1 * 1024 * 1024 * 1024 // 1 GB
		assert.False(t, cfg.IsHighMemoryMode())
	})

	t.Run("exactly_2GB_not_high", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		cfg.MaxMemoryUsage = 2 * 1024 * 1024 * 1024 // 2 GB (not strictly >)
		assert.False(t, cfg.IsHighMemoryMode())
	})

	t.Run("is_high_memory", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		cfg.MaxMemoryUsage = 3 * 1024 * 1024 * 1024 // 3 GB
		assert.True(t, cfg.IsHighMemoryMode())
	})
}

// TestAuthenticityConfig_IsHighPerformanceMode tests high-performance mode detection.
func TestAuthenticityConfig_IsHighPerformanceMode(t *testing.T) {
	t.Run("parallel_disabled", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		cfg.ParallelProcessing = false
		cfg.NumWorkers = 16
		assert.False(t, cfg.IsHighPerformanceMode())
	})

	t.Run("not_enough_workers", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		cfg.ParallelProcessing = true
		cfg.NumWorkers = 4
		assert.False(t, cfg.IsHighPerformanceMode())
	})

	t.Run("exactly_8_workers", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		cfg.ParallelProcessing = true
		cfg.NumWorkers = 8
		assert.True(t, cfg.IsHighPerformanceMode())
	})

	t.Run("many_workers", func(t *testing.T) {
		cfg := DefaultAuthenticityConfig()
		cfg.ParallelProcessing = true
		cfg.NumWorkers = 16
		assert.True(t, cfg.IsHighPerformanceMode())
	})
}

// TestAuthenticityConfig_KmerSizeMethods tests min/max k-mer size helpers.
func TestAuthenticityConfig_KmerSizeMethods(t *testing.T) {
	cfg := DefaultAuthenticityConfig()
	cfg.KmerModelSizes = []int{4, 2, 8, 6}

	assert.Equal(t, 8, cfg.GetMaxKmerSize())
	assert.Equal(t, 2, cfg.GetMinKmerSize())
}

// TestAuthenticityConfig_GetMinKmerSize_Empty tests min kmer size with no sizes.
func TestAuthenticityConfig_GetMinKmerSize_Empty(t *testing.T) {
	cfg := DefaultAuthenticityConfig()
	cfg.KmerModelSizes = []int{}
	// Should return default of 2.
	assert.Equal(t, 2, cfg.GetMinKmerSize())
}

// TestAuthenticityConfig_GetMaxKmerSize_Empty tests max kmer size with no sizes.
func TestAuthenticityConfig_GetMaxKmerSize_Empty(t *testing.T) {
	cfg := DefaultAuthenticityConfig()
	cfg.KmerModelSizes = []int{}
	assert.Equal(t, 0, cfg.GetMaxKmerSize())
}

// TestAuthenticityConfig_SingleKmerSize tests with a single kmer size.
func TestAuthenticityConfig_SingleKmerSize(t *testing.T) {
	cfg := DefaultAuthenticityConfig()
	cfg.KmerModelSizes = []int{5}

	assert.Equal(t, 5, cfg.GetMaxKmerSize())
	assert.Equal(t, 5, cfg.GetMinKmerSize())
}
