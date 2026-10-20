//go:build !release

package population_test

import (
	"reflect"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics/population"
)

// ─────────────────────────────────────────────────────────────────────────────
// PopulationPanel
// ─────────────────────────────────────────────────────────────────────────────

func TestNewPopulationPanel(t *testing.T) {
	panel := population.NewPopulationPanel("YRI", 50)
	require.NotNil(t, panel)
	assert.Equal(t, "YRI", panel.PopulationName)
	assert.Equal(t, 50, panel.SampleSize)
	assert.Equal(t, 12300, panel.Ne)
	assert.Equal(t, population.PresetMutationRate, panel.MutationRate)

	other := population.NewPopulationPanel("not-a-preset", 5)
	assert.Equal(t, 10000, other.Ne)
	assert.Equal(t, 1e-8, other.MutationRate)
}

func TestPopulationPanel_LoadFromSimulation(t *testing.T) {
	panel := population.NewPopulationPanel("CEU", 5)
	panel.SetSeed(7)
	require.NoError(t, panel.LoadFromSimulation(200))
	assert.NotNil(t, panel.GetVariants())
}

// loadPanel builds a panel with an elevated mutation rate so that a short
// region carries variants.
func loadPanel(t *testing.T, name string, seed int64, seedIt bool) *population.PopulationPanel {
	t.Helper()
	p := population.NewPopulationPanel(name, 8)
	p.Ne, p.MutationRate = 100, 1e-4
	if seedIt {
		p.SetSeed(seed)
	}
	require.NoError(t, p.LoadFromSimulation(5000))
	return p
}

// LoadFromSimulation draws its seed from the panel's own random source: the
// default source is fixed (two fresh panels agree), SetSeed selects another,
// and successive loads differ.
func TestPopulationPanel_LoadFromSimulation_Deterministic(t *testing.T) {
	a := loadPanel(t, "XYZ", 0, false)
	b := loadPanel(t, "XYZ", 0, false)
	require.NotEmpty(t, a.GetVariants())
	assert.True(t, reflect.DeepEqual(a.GetVariants(), b.GetVariants()), "default seed must be fixed")

	c := loadPanel(t, "XYZ", 42, true)
	d := loadPanel(t, "XYZ", 42, true)
	e := loadPanel(t, "XYZ", 43, true)
	assert.True(t, reflect.DeepEqual(c.GetVariants(), d.GetVariants()), "same seed")
	assert.False(t, reflect.DeepEqual(c.GetVariants(), e.GetVariants()), "different seed")

	first := c.GetVariants()
	require.NoError(t, c.LoadFromSimulation(5000))
	assert.False(t, reflect.DeepEqual(first, c.GetVariants()), "successive loads draw new seeds")
}

func TestPopulationPanel_GetVariants_ValidFrequencies(t *testing.T) {
	panel := loadPanel(t, "CHB", 111, true)
	for _, v := range panel.GetVariants() {
		assert.Greater(t, v.Frequency, 0.0)
		assert.Less(t, v.Frequency, 1.0)
		assert.NotZero(t, v.ReferenceBase)
		assert.NotZero(t, v.AlternateBase)
	}
}

func TestPopulationPanel_GetAlleleFrequencies_Empty(t *testing.T) {
	panel := population.NewPopulationPanel("empty_pop", 5)
	assert.Empty(t, panel.GetAlleleFrequencies())
}

func TestPopulationPanel_GetAlleleFrequencies_AfterLoad(t *testing.T) {
	panel := loadPanel(t, "JPT", 222, true)
	freqs := panel.GetAlleleFrequencies()
	require.NotEmpty(t, freqs)
	// Raw derived-allele frequencies (not normalised) in [0, 1].
	for site, freq := range freqs {
		assert.GreaterOrEqualf(t, freq, 0.0, "site %s", site)
		assert.LessOrEqualf(t, freq, 1.0, "site %s", site)
	}
}

func TestPopulationPanel_GetAlleleFrequencies_NormalizedFromVariants(t *testing.T) {
	// When AlleleFrequencies is empty but Variants are present,
	// GetAlleleFrequencies normalizes to sum == 1.
	panel := population.NewPopulationPanel("JPT", 6)
	panel.Variants = []population.Variant{
		{Position: 100, Chromosome: "chr1", ReferenceBase: 'A', AlternateBase: 'T', Frequency: 0.3},
		{Position: 200, Chromosome: "chr1", ReferenceBase: 'C', AlternateBase: 'G', Frequency: 0.5},
		{Position: 300, Chromosome: "chr1", ReferenceBase: 'T', AlternateBase: 'A', Frequency: 0.1},
	}

	freqs := panel.GetAlleleFrequencies()
	require.Len(t, freqs, 3)
	total := 0.0
	for _, f := range freqs {
		total += f
	}
	assert.InDelta(t, 1.0, total, 1e-9, "allele frequencies should sum to 1")
}

func TestPopulationPanel_GetDemographics(t *testing.T) {
	panel := population.NewPopulationPanel("AFR", 10)
	assert.NotEmpty(t, panel.GetDemographics())
}

// A demographic event set on the panel changes the simulation: a severe
// persistent bottleneck from 5 generations ago (N=2) collapses the genealogy,
// so far fewer sites segregate than under the constant-size model.
func TestPopulationPanel_SetDemographicEvent(t *testing.T) {
	segregating := func(withEvent bool) int {
		p := population.NewPopulationPanel("XYZ", 10)
		p.Ne, p.MutationRate = 1000, 1e-5
		p.SetSeed(5)
		if withEvent {
			p.SetDemographicEvent(population.DemographicEvent{
				Type: "bottleneck", Time: 5,
				Parameters: map[string]float64{"population_size": 2},
			})
			assert.Equal(t, 1, p.GetPopulationStats()["demographic_events"])
		}
		require.NoError(t, p.LoadFromSimulation(100_000))
		return len(p.GetVariants())
	}
	base, bottleneck := segregating(false), segregating(true)
	t.Logf("segregating sites: constant %d, bottleneck %d", base, bottleneck)
	assert.Less(t, bottleneck*4, base)
}

func TestPopulationPanel_GetPopulationStats(t *testing.T) {
	panel := population.NewPopulationPanel("SAS", 5)
	panel.Ne, panel.MutationRate = 100, 1e-4
	panel.SetSeed(333)
	require.NoError(t, panel.LoadFromSimulation(200))

	stats := panel.GetPopulationStats()
	assert.Equal(t, "SAS", stats["population_name"])
	assert.Equal(t, 5, stats["sample_size"])
	assert.Equal(t, 100, stats["effective_population_size"])
	het := stats["heterozygosity"].(float64)
	assert.GreaterOrEqual(t, het, 0.0)
	assert.LessOrEqual(t, het, 0.5)
}
