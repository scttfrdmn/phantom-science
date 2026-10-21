//go:build !release

// Package e2e checks the internal/testing E2E framework end to end through
// its public API: reads come from the real generator, every value asserted
// on is computed from those reads, and platforms the generator does not
// model are skipped with testing.ErrPlatformNotModeled. External tools are
// exercised in tests/integration.
package e2e

import (
	"errors"
	"os"
	"path/filepath"
	"reflect"
	"regexp"
	"strconv"
	"strings"
	"sync"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	pt "github.com/scttfrdmn/phantom-science/internal/testing"
)

// Header formats written by the engine.
var (
	illuminaHeader = regexp.MustCompile(`^[^\s:]+:\d+:[^\s:]+:\d+:\d+:\d+:\d+ [12]:[YN]:\d+:\S+`)
	ontHeader      = regexp.MustCompile(`^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12} runid=[0-9a-f]{40} read=\d+ ch=\d+ start_time=\S+ flow_cell_id=\S+ protocol_group_id=\S+ sample_id=\S+`)
)

// generateOrSkip generates reads for config, skipping the test when the
// generator has no model for the platform.
func generateOrSkip(t *testing.T, f *pt.E2EFramework, c *pt.PlatformConfig) *pt.SyntheticData {
	t.Helper()
	data, err := f.GenerateSyntheticData(c)
	if errors.Is(err, pt.ErrPlatformNotModeled) {
		t.Skipf("%s: %v", c.Platform, err)
	}
	require.NoError(t, err)
	return data
}

// TestPlatforms generates reads for every platform profile, writes them as
// FASTQ and checks the file against the profile: read count, header format,
// read lengths and Phred+33 qualities, and that the reported statistics
// equal those of the file read back.
func TestPlatforms(t *testing.T) {
	f := pt.NewE2EFramework()
	for _, platform := range pt.AllPlatforms {
		t.Run(string(platform), func(t *testing.T) {
			c := f.CreatePlatformConfig(platform, "YRI")
			require.NotNil(t, c)
			data := generateOrSkip(t, f, c)

			dir := t.TempDir()
			r1 := filepath.Join(dir, "reads_R1.fastq")
			require.NoError(t, f.WriteFASTQ(data, r1, filepath.Join(dir, "reads_R2.fastq")))
			back, err := pt.LoadFASTQ(r1, 0)
			require.NoError(t, err)
			require.Len(t, back.Sequences, c.Reads)
			assert.Equal(t, data.Stats, back.Stats, "stats of the written FASTQ")

			header := illuminaHeader
			if c.EnginePlatform == "ont" {
				header = ontHeader
			}
			for i, h := range back.Headers {
				if !header.MatchString(h) {
					t.Fatalf("read %d header %q does not match the %s format", i, h, c.EnginePlatform)
				}
			}
			if c.ReadLengthDistribution.Dist == "" {
				assert.Equal(t, c.ReadLength, back.Stats.MinReadLength)
				assert.Equal(t, c.ReadLength, back.Stats.MaxReadLength)
			} else {
				assert.GreaterOrEqual(t, back.Stats.MinReadLength, c.ReadLengthDistribution.Min)
				assert.LessOrEqual(t, back.Stats.MaxReadLength, c.ReadLengthDistribution.Max)
				assert.Greater(t, back.Stats.N50, int(back.Stats.MeanReadLength*0.9),
					"N50 of a right-skewed length distribution should not be far below the mean")
			}
			assert.GreaterOrEqual(t, back.Stats.MeanQuality, 2.0)
			assert.LessOrEqual(t, back.Stats.MeanQuality, 93.0)
			t.Logf("%s: %+v", platform, back.Stats)
		})
	}
}

// TestPositiveControl is the classifier calibration on each modeled
// platform, for every feature set: two independent seeded draws must be
// indistinguishable (holdout accuracy within 0.5 ± 4σ); skewing one draw's
// GC content to 0.7 must be detected (accuracy above 0.9) by the sets with
// sequence features, and the other quality model by the sets with quality
// features.
func TestPositiveControl(t *testing.T) {
	f := pt.NewE2EFramework()
	for _, platform := range pt.AllPlatforms {
		t.Run(string(platform), func(t *testing.T) {
			ca, cb := f.CreatePlatformConfig(platform, "CEU"), f.CreatePlatformConfig(platform, "CEU")
			ca.Seed, cb.Seed = 1001, 1002
			ca.Reads, cb.Reads = 1000, 1000
			a := generateOrSkip(t, f, ca)
			b := generateOrSkip(t, f, cb)
			alt := generateOrSkip(t, f, cb.AlternativeQualityModel())
			gc := pt.GCSkewed(b, 0.7, 3)

			for _, set := range pt.AllFeatureSets {
				null, err := pt.HoldoutAccuracyFor(set, a, b, 0)
				require.NoError(t, err)
				band := pt.ChanceBand(null.Holdout, 4)
				assert.InDelta(t, 0.5, null.Accuracy, band, "%s null accuracy on %d held-out reads", set, null.Holdout)

				gcAcc, err := pt.HoldoutAccuracyFor(set, a, gc, 0)
				require.NoError(t, err)
				qAcc, err := pt.HoldoutAccuracyFor(set, a, alt, 0)
				require.NoError(t, err)
				if set != pt.FeatureQuality {
					assert.Greater(t, gcAcc.Accuracy, 0.9, "%s GC-skewed control accuracy", set)
				}
				if set != pt.FeatureSequence {
					assert.Greater(t, qAcc.Accuracy, 0.9, "%s quality-model control accuracy", set)
				}
				t.Logf("%s %s: null %.4f (band ±%.4f, %d held out), GC control %.4f, quality control %.4f",
					platform, set, null.Accuracy, band, null.Holdout, gcAcc.Accuracy, qAcc.Accuracy)
			}
		})
	}
}

// TestRealReads compares real reads with generated ones when a FASTQ is
// supplied in PHANTOM_REAL_FASTQ (platform from PHANTOM_REAL_PLATFORM,
// default illumina_novaseq_6000; PHANTOM_REAL_READS reads, default 2000).
// Reads are generated from uniform random bases, or from the reference
// FASTA in PHANTOM_REFERENCE (optionally only PHANTOM_REGION of it).
//
// For each feature set the real-versus-generated accuracy is reported, not
// graded. The test fails if alternate halves of the real reads can be told
// apart (the real-versus-real null, outside 0.5 ± 4σ), or if a feature set
// has no power on the real reads: a 0.3 GC shift must be detected by the
// sets with sequence features, a one-Phred quality shift by the sets with
// quality features. (Real reads vary in GC far more from read to read than
// uniform ones, so 0.2, the shift used on generated reads, is not enough:
// on 1000 real NovaSeq reads it gives 0.896.)
func TestRealReads(t *testing.T) {
	path := os.Getenv("PHANTOM_REAL_FASTQ")
	if path == "" {
		t.Skip("PHANTOM_REAL_FASTQ not set: no real reads to compare with")
	}
	platform := pt.SequencingPlatform(os.Getenv("PHANTOM_REAL_PLATFORM"))
	if platform == "" {
		platform = pt.IlluminaNovaSeq6000
	}
	n := 2000
	if s := os.Getenv("PHANTOM_REAL_READS"); s != "" {
		var err error
		n, err = strconv.Atoi(s)
		require.NoError(t, err, "PHANTOM_REAL_READS")
	}
	real, err := pt.LoadFASTQ(path, n)
	require.NoError(t, err)
	require.GreaterOrEqual(t, len(real.Sequences), 8, "need at least eight real reads")

	f := pt.NewE2EFramework()
	c := f.CreatePlatformConfig(platform, "CEU")
	require.NotNil(t, c, "unknown PHANTOM_REAL_PLATFORM %q", platform)
	c.Reads = len(real.Sequences)
	c.ReferencePath, c.ReferenceRegion = os.Getenv("PHANTOM_REFERENCE"), os.Getenv("PHANTOM_REGION")
	if real.Stats.MinReadLength == real.Stats.MaxReadLength && c.ReadLengthDistribution.Dist == "" {
		c.ReadLength = real.Stats.MinReadLength
	}
	synth := generateOrSkip(t, f, c)
	t.Logf("real (%s, %d reads) vs %s generated from %s %s %s", path, len(real.Sequences), platform, c.BaseSource(), c.ReferencePath, c.ReferenceRegion)
	t.Logf("real stats %+v; generated stats %+v", real.Stats, synth.Stats)

	even, odd := pt.SplitAlternate(real)
	gcOdd := pt.GCSkewed(odd, odd.Stats.GCContent+0.3, 1)
	qOdd := pt.QualityShifted(odd, 1)
	for _, set := range pt.AllFeatureSets {
		gen, err := pt.HoldoutAccuracyFor(set, real, synth, 5)
		require.NoError(t, err)
		t.Logf("%s (%s): real vs generated %.4f on %d held-out reads (chance band 0.5 ± %.4f); top features %+v",
			set, gen.Classifier, gen.Accuracy, gen.Holdout, pt.ChanceBand(gen.Holdout, 4), gen.Top)

		null, err := pt.HoldoutAccuracyFor(set, even, odd, 0)
		require.NoError(t, err)
		band := pt.ChanceBand(null.Holdout, 4)
		assert.InDelta(t, 0.5, null.Accuracy, band, "%s: real vs real (alternate halves) on %d held-out reads", set, null.Holdout)

		if set != pt.FeatureQuality {
			power, err := pt.HoldoutAccuracyFor(set, even, gcOdd, 0)
			require.NoError(t, err)
			assert.Greater(t, power.Accuracy, 0.9, "%s: power on real reads (GC +0.3)", set)
		}
		if set != pt.FeatureSequence {
			power, err := pt.HoldoutAccuracyFor(set, even, qOdd, 0)
			require.NoError(t, err)
			assert.Greater(t, power.Accuracy, 0.9, "%s: power on real reads (quality −1)", set)
		}
	}
}

// TestConcurrentGeneration generates every modeled platform concurrently
// from one framework and checks each result equals a sequential run with
// the same seed (run with -race to check the framework's locking).
func TestConcurrentGeneration(t *testing.T) {
	f := pt.NewE2EFramework()
	var platforms []pt.SequencingPlatform
	for _, p := range pt.AllPlatforms {
		if f.Modeled(p) {
			platforms = append(platforms, p)
		}
	}
	require.NotEmpty(t, platforms)

	got := make([]*pt.SyntheticData, len(platforms))
	errs := make([]error, len(platforms))
	var wg sync.WaitGroup
	for i, p := range platforms {
		wg.Add(1)
		go func() {
			defer wg.Done()
			c := f.CreatePlatformConfig(p, "CEU")
			c.Reads = 200
			c.TechnologyParams["worker"] = strings.Repeat("x", i) // must not leak into other configs
			got[i], errs[i] = f.GenerateSyntheticData(c)
		}()
	}
	wg.Wait()

	for i, p := range platforms {
		require.NoError(t, errs[i], "%s", p)
		c := f.CreatePlatformConfig(p, "CEU")
		c.Reads = 200
		_, leaked := c.TechnologyParams["worker"]
		assert.False(t, leaked, "%s: a config mutation leaked into the framework profile", p)
		want, err := f.GenerateSyntheticData(c)
		require.NoError(t, err)
		assert.True(t, reflect.DeepEqual(want, got[i]), "%s: concurrent result differs from a sequential run", p)
	}
}
