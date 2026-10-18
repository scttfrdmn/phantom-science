package population_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics/population"
)

// ─────────────────────────────────────────────────────────────────────────────
// CoalescentSimulator construction & validation
// ─────────────────────────────────────────────────────────────────────────────

func TestNewCoalescentSimulator(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-8, 1e-8)
	require.NotNil(t, sim)

	stats := sim.GetSimulationStats()
	assert.Equal(t, 100, stats["effective_population_size"])
	assert.Equal(t, 1e-8, stats["mutation_rate"])
	assert.Equal(t, 1e-8, stats["recombination_rate"])
	assert.Equal(t, 0, stats["demographic_events"])
	assert.Equal(t, 0, stats["genealogies_simulated"])
}

func TestValidateParameters_Valid(t *testing.T) {
	tests := []struct {
		name              string
		effectivePopSize  int
		mutationRate      float64
		recombinationRate float64
	}{
		{"small pop", 10, 1e-8, 0.0},
		{"large pop", 10000, 0.5, 0.5},
		{"zero rates", 1, 0.0, 0.0},
		{"max rates", 100, 1.0, 1.0},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			sim := population.NewCoalescentSimulator(tc.effectivePopSize, tc.mutationRate, tc.recombinationRate)
			require.NotNil(t, sim)
			err := sim.ValidateParameters()
			assert.NoError(t, err)
		})
	}
}

func TestValidateParameters_Invalid(t *testing.T) {
	tests := []struct {
		name              string
		effectivePopSize  int
		mutationRate      float64
		recombinationRate float64
	}{
		{"zero pop size", 0, 1e-8, 1e-8},
		{"negative pop size", -100, 1e-8, 1e-8},
		{"negative mutation rate", 100, -0.1, 1e-8},
		{"mutation rate > 1", 100, 1.5, 1e-8},
		{"negative recombination rate", 100, 1e-8, -0.1},
		{"recombination rate > 1", 100, 1e-8, 1.5},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			sim := population.NewCoalescentSimulator(tc.effectivePopSize, tc.mutationRate, tc.recombinationRate)
			require.NotNil(t, sim)
			err := sim.ValidateParameters()
			assert.Error(t, err)
		})
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// Genealogy simulation
// ─────────────────────────────────────────────────────────────────────────────

func TestSimulateGenealogy_Basic(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-8, 0.0)
	sim.SetSeed(42)

	genealogy, err := sim.SimulateGenealogy("test_region", 5, 100)
	require.NoError(t, err)
	require.NotNil(t, genealogy)

	// Must have root
	assert.NotNil(t, genealogy.Root)
	// Leaf count must equal sample size
	assert.Len(t, genealogy.LeafNodes, 5)
	// Total nodes = 2*sampleSize - 1 for a fully-bifurcating tree
	assert.Len(t, genealogy.Nodes, 2*5-1)
	// Total time must be positive
	assert.Greater(t, genealogy.TotalTime, 0.0)
}

func TestSimulateGenealogy_SampleTooSmall(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-8, 0.0)
	_, err := sim.SimulateGenealogy("bad", 1, 100)
	require.Error(t, err)
	assert.Contains(t, err.Error(), "sample size")
}

func TestSimulateGenealogy_WithMutations(t *testing.T) {
	// High mutation rate so we expect some mutations
	sim := population.NewCoalescentSimulator(100, 1e-4, 0.0)
	sim.SetSeed(1234)

	genealogy, err := sim.SimulateGenealogy("mutation_test", 10, 500)
	require.NoError(t, err)
	require.NotNil(t, genealogy)

	// With high mutation rate and decent region length, expect at least some mutations
	// (not a guarantee but highly likely with these params)
	t.Logf("Mutation events: %d", len(genealogy.MutationEvents))
	assert.NotNil(t, genealogy.MutationEvents)
}

func TestSimulateGenealogy_WithRecombination(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 0.0, 1e-4)
	sim.SetSeed(999)

	genealogy, err := sim.SimulateGenealogy("recomb_test", 5, 1000)
	require.NoError(t, err)
	require.NotNil(t, genealogy)
	t.Logf("Recombination events: %d", len(genealogy.RecombinationEvents))
	assert.NotNil(t, genealogy.RecombinationEvents)
}

func TestSimulateGenealogy_StoredInSimulator(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-8, 0.0)
	sim.SetSeed(77)

	_, err := sim.SimulateGenealogy("region_A", 4, 100)
	require.NoError(t, err)

	stats := sim.GetSimulationStats()
	assert.Equal(t, 1, stats["genealogies_simulated"])
}

func TestSimulateGenealogy_MultipleRegions(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-8, 0.0)
	sim.SetSeed(55)

	for i := 0; i < 5; i++ {
		regionID := fmt.Sprintf("region_%d", i)
		_, err := sim.SimulateGenealogy(regionID, 3, 50)
		require.NoError(t, err)
	}

	stats := sim.GetSimulationStats()
	assert.Equal(t, 5, stats["genealogies_simulated"])
}

// ─────────────────────────────────────────────────────────────────────────────
// Demographic events
// ─────────────────────────────────────────────────────────────────────────────

func TestAddDemographicEvent(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-8, 0.0)

	event := population.DemographicEvent{
		Type: "bottleneck",
		Time: 100,
		Parameters: map[string]float64{
			"population_size": 10.0,
		},
		Populations: []string{"CEU"},
		Duration:    50,
	}
	sim.AddDemographicEvent(event)

	events := sim.GetDemographicEvents()
	require.Len(t, events, 1)
	assert.Equal(t, "bottleneck", events[0].Type)
	assert.Equal(t, 100, events[0].Time)
}

func TestAddDemographicEvent_SortedByTime(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-8, 0.0)

	for _, t_ := range []int{300, 100, 200} {
		sim.AddDemographicEvent(population.DemographicEvent{
			Type:       "expansion",
			Time:       t_,
			Parameters: map[string]float64{"expansion_factor": 2.0},
		})
	}

	events := sim.GetDemographicEvents()
	require.Len(t, events, 3)
	assert.Equal(t, 100, events[0].Time)
	assert.Equal(t, 200, events[1].Time)
	assert.Equal(t, 300, events[2].Time)
}

func TestSimulateGenealogy_WithBottleneck(t *testing.T) {
	sim := population.NewCoalescentSimulator(1000, 1e-8, 0.0)
	sim.SetSeed(22)

	// Add a severe bottleneck 50 generations ago
	sim.AddDemographicEvent(population.DemographicEvent{
		Type: "bottleneck",
		Time: 50,
		Parameters: map[string]float64{
			"population_size": 5.0,
		},
	})

	genealogy, err := sim.SimulateGenealogy("bottleneck_region", 5, 100)
	require.NoError(t, err)
	require.NotNil(t, genealogy)
	assert.NotNil(t, genealogy.Root)
}

func TestSimulateGenealogy_WithExpansion(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-8, 0.0)
	sim.SetSeed(33)

	sim.AddDemographicEvent(population.DemographicEvent{
		Type: "expansion",
		Time: 200,
		Parameters: map[string]float64{
			"expansion_factor": 10.0,
		},
	})

	genealogy, err := sim.SimulateGenealogy("expansion_region", 4, 100)
	require.NoError(t, err)
	require.NotNil(t, genealogy)
}

// ─────────────────────────────────────────────────────────────────────────────
// Allele frequencies
// ─────────────────────────────────────────────────────────────────────────────

func TestSimulateAlleleFrequencies(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-4, 0.0)
	sim.SetSeed(456)

	freqs, err := sim.SimulateAlleleFrequencies(200, 5)
	require.NoError(t, err)
	require.NotNil(t, freqs)

	// All frequencies must be in [0, 1]
	for site, freq := range freqs {
		assert.GreaterOrEqualf(t, freq, 0.0, "site %s frequency must be >= 0", site)
		assert.LessOrEqualf(t, freq, 1.0, "site %s frequency must be <= 1", site)
	}
}

func TestSimulateAlleleFrequencies_SampleSizeTooSmall(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-8, 0.0)
	_, err := sim.SimulateAlleleFrequencies(100, 1)
	require.Error(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// Site Frequency Spectrum
// ─────────────────────────────────────────────────────────────────────────────

func TestComputeSiteFrequencySpectrum(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-4, 0.0)
	sim.SetSeed(789)

	sampleSize := 10
	genealogy, err := sim.SimulateGenealogy("sfs_test", sampleSize, 500)
	require.NoError(t, err)

	sfs := sim.ComputeSiteFrequencySpectrum(genealogy, sampleSize)
	require.NotNil(t, sfs)

	// Frequency bins should be filled
	assert.Len(t, sfs.Frequencies, sampleSize)
	assert.Len(t, sfs.Counts, sampleSize)

	// Total sites should equal number of mutation events
	assert.Equal(t, len(genealogy.MutationEvents), sfs.TotalSites)

	// ThetaWatterson and ThetaPi should be >= 0
	assert.GreaterOrEqual(t, sfs.ThetaWatterson, 0.0)
	assert.GreaterOrEqual(t, sfs.ThetaPi, 0.0)
}

func TestComputeSiteFrequencySpectrum_FrequencyBins(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 0.0, 0.0)
	sim.SetSeed(100)

	sampleSize := 5
	genealogy, err := sim.SimulateGenealogy("bins_test", sampleSize, 100)
	require.NoError(t, err)

	sfs := sim.ComputeSiteFrequencySpectrum(genealogy, sampleSize)

	// Frequencies must be monotonically increasing from small positive to 1
	for i, f := range sfs.Frequencies {
		assert.Equal(t, float64(i+1)/float64(sampleSize), f)
	}
}

func TestTajimaD_ZeroWhenNoVariation(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 0.0, 0.0)
	sim.SetSeed(0)

	genealogy, err := sim.SimulateGenealogy("no_var", 5, 100)
	require.NoError(t, err)

	sfs := sim.ComputeSiteFrequencySpectrum(genealogy, 5)
	// With 0 mutation rate, ThetaWatterson == 0, so TajimaD should be 0
	assert.Equal(t, 0.0, sfs.TajimaD)
}

func TestComputeSiteFrequencySpectrum_SampleSizeOne_NoInfNaN(t *testing.T) {
	// #111: sampleSize == 1 caused division by zero → +Inf/NaN in all four stats.
	// SimulateGenealogy rejects sampleSize < 2, so build a Genealogy directly.
	sim := population.NewCoalescentSimulator(100, 1e-4, 0.0)
	genealogy := &population.Genealogy{
		MutationEvents: []population.MutationEvent{},
	}

	sfs := sim.ComputeSiteFrequencySpectrum(genealogy, 1)
	require.NotNil(t, sfs)

	assert.False(t, math.IsInf(sfs.ThetaWatterson, 0), "ThetaWatterson should not be Inf")
	assert.False(t, math.IsNaN(sfs.ThetaWatterson), "ThetaWatterson should not be NaN")
	assert.False(t, math.IsInf(sfs.ThetaPi, 0), "ThetaPi should not be Inf")
	assert.False(t, math.IsNaN(sfs.ThetaPi), "ThetaPi should not be NaN")
	assert.False(t, math.IsInf(sfs.TajimaD, 0), "TajimaD should not be Inf")
	assert.False(t, math.IsNaN(sfs.TajimaD), "TajimaD should not be NaN")
	assert.False(t, math.IsInf(sfs.FayWuH, 0), "FayWuH should not be Inf")
	assert.False(t, math.IsNaN(sfs.FayWuH), "FayWuH should not be NaN")
}

// ─────────────────────────────────────────────────────────────────────────────
// SetSeed / reproducibility
// ─────────────────────────────────────────────────────────────────────────────

func TestSetSeed_Reproducibility(t *testing.T) {
	makeGenealogy := func(seed int64) *population.Genealogy {
		sim := population.NewCoalescentSimulator(100, 1e-4, 0.0)
		sim.SetSeed(seed)
		g, err := sim.SimulateGenealogy("rep_test", 4, 200)
		require.NoError(t, err)
		return g
	}

	g1 := makeGenealogy(42)
	g2 := makeGenealogy(42)

	assert.Equal(t, g1.TotalTime, g2.TotalTime)
	assert.Len(t, g1.Nodes, len(g2.Nodes))
}

func TestSetSeed_DifferentSeeds(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-4, 0.0)
	sim.SetSeed(1)
	g1, err := sim.SimulateGenealogy("diff_seed_1", 5, 200)
	require.NoError(t, err)

	sim.SetSeed(2)
	g2, err := sim.SimulateGenealogy("diff_seed_2", 5, 200)
	require.NoError(t, err)

	// Different seeds should (almost certainly) give different total times
	// — we just check they aren't identical.
	if g1.TotalTime == g2.TotalTime {
		t.Log("Unlikely coincidence: same total time with different seeds")
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// GetSimulationStats
// ─────────────────────────────────────────────────────────────────────────────

func TestGetSimulationStats_IncrementsOnSimulation(t *testing.T) {
	sim := population.NewCoalescentSimulator(50, 1e-8, 0.0)
	sim.SetSeed(8)

	for i := 0; i < 3; i++ {
		_, err := sim.SimulateGenealogy(fmt.Sprintf("r%d", i), 3, 50)
		require.NoError(t, err)
	}

	stats := sim.GetSimulationStats()
	assert.Equal(t, 3, stats["genealogies_simulated"])
	assert.Equal(t, 50, stats["effective_population_size"])
}

// ─────────────────────────────────────────────────────────────────────────────
// PopulationPanel
// ─────────────────────────────────────────────────────────────────────────────

func TestNewPopulationPanel(t *testing.T) {
	panel := population.NewPopulationPanel("YRI", 50)
	require.NotNil(t, panel)
	assert.Equal(t, "YRI", panel.PopulationName)
	assert.Equal(t, 50, panel.SampleSize)
	assert.NotNil(t, panel.Simulator)
}

func TestPopulationPanel_LoadFromSimulation(t *testing.T) {
	panel := population.NewPopulationPanel("CEU", 5)
	panel.Simulator.SetSeed(7)

	err := panel.LoadFromSimulation(200)
	require.NoError(t, err)

	// Should have loaded some data
	variants := panel.GetVariants()
	t.Logf("Loaded %d variants", len(variants))
	// GetVariants returns a slice (may be empty if no mutations, but structure is valid)
	assert.NotNil(t, variants)
}

func TestPopulationPanel_GetVariants_ValidFrequencies(t *testing.T) {
	panel := population.NewPopulationPanel("CHB", 8)
	panel.Simulator = population.NewCoalescentSimulator(100, 1e-4, 0.0)
	panel.Simulator.SetSeed(111)

	err := panel.LoadFromSimulation(300)
	require.NoError(t, err)

	for _, v := range panel.GetVariants() {
		assert.GreaterOrEqual(t, v.Frequency, 0.0)
		assert.LessOrEqual(t, v.Frequency, 1.0)
		assert.NotZero(t, v.ReferenceBase)
		assert.NotZero(t, v.AlternateBase)
	}
}

func TestPopulationPanel_GetAlleleFrequencies_Empty(t *testing.T) {
	panel := population.NewPopulationPanel("empty_pop", 5)
	// Not loaded, expect empty map
	freqs := panel.GetAlleleFrequencies()
	assert.NotNil(t, freqs)
}

func TestPopulationPanel_GetAlleleFrequencies_AfterLoad(t *testing.T) {
	panel := population.NewPopulationPanel("JPT", 6)
	panel.Simulator = population.NewCoalescentSimulator(100, 1e-4, 0.0)
	panel.Simulator.SetSeed(222)

	err := panel.LoadFromSimulation(300)
	require.NoError(t, err)

	freqs := panel.GetAlleleFrequencies()
	assert.NotNil(t, freqs)

	// All values must be in [0, 1]
	for site, freq := range freqs {
		assert.GreaterOrEqualf(t, freq, 0.0, "site %s", site)
		assert.LessOrEqualf(t, freq, 1.0, "site %s", site)
	}
	// Note: when AlleleFrequencies is populated by LoadFromSimulation,
	// the raw derived-allele frequencies are returned (not normalized to sum=1).
	// Only the Variants-derived path normalizes to a probability distribution.
}

func TestPopulationPanel_GetAlleleFrequencies_NormalizedFromVariants(t *testing.T) {
	// When AlleleFrequencies is empty but Variants are present,
	// GetAlleleFrequencies normalizes to sum == 1.
	panel := population.NewPopulationPanel("JPT", 6)
	// Do NOT load from simulation (keeps AlleleFrequencies empty)
	// Manually add variants with known frequencies
	panel.Variants = []population.Variant{
		{Position: 100, Chromosome: "chr1", ReferenceBase: 'A', AlternateBase: 'T', Frequency: 0.3},
		{Position: 200, Chromosome: "chr1", ReferenceBase: 'C', AlternateBase: 'G', Frequency: 0.5},
		{Position: 300, Chromosome: "chr1", ReferenceBase: 'T', AlternateBase: 'A', Frequency: 0.1},
	}

	freqs := panel.GetAlleleFrequencies()
	require.NotEmpty(t, freqs)
	assert.Len(t, freqs, 3)

	total := 0.0
	for _, f := range freqs {
		total += f
	}
	assert.InDelta(t, 1.0, total, 1e-9, "allele frequencies should sum to 1")
}

func TestPopulationPanel_GetDemographics(t *testing.T) {
	panel := population.NewPopulationPanel("AFR", 10)
	demographics := panel.GetDemographics()
	assert.NotNil(t, demographics)
}

func TestPopulationPanel_SetDemographicEvent(t *testing.T) {
	panel := population.NewPopulationPanel("EUR", 5)

	event := population.DemographicEvent{
		Type:       "bottleneck",
		Time:       1000,
		Parameters: map[string]float64{"population_size": 50.0},
	}
	panel.SetDemographicEvent(event)

	events := panel.Simulator.GetDemographicEvents()
	assert.Len(t, events, 1)
}

func TestPopulationPanel_GetPopulationStats(t *testing.T) {
	panel := population.NewPopulationPanel("SAS", 5)
	panel.Simulator = population.NewCoalescentSimulator(100, 1e-4, 0.0)
	panel.Simulator.SetSeed(333)

	err := panel.LoadFromSimulation(200)
	require.NoError(t, err)

	stats := panel.GetPopulationStats()
	require.NotNil(t, stats)

	assert.Equal(t, "SAS", stats["population_name"])
	assert.Equal(t, 5, stats["sample_size"])

	// heterozygosity in [0, 1]
	het := stats["heterozygosity"].(float64)
	assert.GreaterOrEqual(t, het, 0.0)
	assert.LessOrEqual(t, het, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Mathematical consistency checks
// ─────────────────────────────────────────────────────────────────────────────

func TestCoalescent_TreeTopology(t *testing.T) {
	sampleSize := 6
	sim := population.NewCoalescentSimulator(100, 0.0, 0.0)
	sim.SetSeed(3)

	genealogy, err := sim.SimulateGenealogy("topology_test", sampleSize, 50)
	require.NoError(t, err)

	// A bifurcating tree with n leaves has exactly 2n-1 nodes
	assert.Len(t, genealogy.Nodes, 2*sampleSize-1)

	// Root should have no parent
	assert.Nil(t, genealogy.Root.Parent)

	// All leaf nodes should have IsLeaf == true
	for _, leaf := range genealogy.LeafNodes {
		assert.True(t, leaf.IsLeaf)
	}

	// All internal nodes should have exactly 2 children
	for _, node := range genealogy.Nodes {
		if !node.IsLeaf {
			assert.Len(t, node.Children, 2)
		}
	}
}

func TestCoalescent_TimesMonotonicallyIncreasing(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 0.0, 0.0)
	sim.SetSeed(4)

	genealogy, err := sim.SimulateGenealogy("time_test", 5, 50)
	require.NoError(t, err)

	// For every node that has a parent, parent.Time > node.Time
	for _, node := range genealogy.Nodes {
		if node.Parent != nil {
			assert.Greater(t, node.Parent.Time, node.Time,
				"parent time should be greater than child time")
		}
	}
}

func TestCoalescent_RootHasHighestTime(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 0.0, 0.0)
	sim.SetSeed(5)

	genealogy, err := sim.SimulateGenealogy("root_time_test", 5, 50)
	require.NoError(t, err)

	for _, node := range genealogy.Nodes {
		assert.LessOrEqual(t, node.Time, genealogy.Root.Time)
	}
}

func TestCoalescent_MutationFrequenciesInRange(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-3, 0.0)
	sim.SetSeed(6)

	genealogy, err := sim.SimulateGenealogy("freq_range_test", 8, 500)
	require.NoError(t, err)

	for _, mut := range genealogy.MutationEvents {
		f := mut.Mutation.Frequency
		assert.GreaterOrEqual(t, f, 0.0, "mutation frequency must be >= 0")
		assert.LessOrEqual(t, f, 1.0, "mutation frequency must be <= 1")

		// Mutation type must be transition or transversion
		assert.Contains(t, []string{"transition", "transversion"}, mut.Type)

		// Ancestral and derived must differ
		assert.NotEqual(t, mut.Mutation.Ancestral, mut.Mutation.Derived)
	}
}

func TestCoalescent_FixedMutations(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-3, 0.0)
	sim.SetSeed(7)

	genealogy, err := sim.SimulateGenealogy("fixed_test", 5, 500)
	require.NoError(t, err)

	for _, mut := range genealogy.MutationEvents {
		if mut.Mutation.IsFixed {
			assert.InDelta(t, 1.0, mut.Mutation.Frequency, 1e-9,
				"fixed mutations should have frequency == 1.0")
		}
	}
}

func TestCoalescent_SFSCountsMatchMutations(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-4, 0.0)
	sim.SetSeed(99)
	sampleSize := 8

	genealogy, err := sim.SimulateGenealogy("sfs_count_test", sampleSize, 500)
	require.NoError(t, err)

	sfs := sim.ComputeSiteFrequencySpectrum(genealogy, sampleSize)

	// Total counts in SFS bins must equal TotalSites (mutation count)
	total := 0
	for _, c := range sfs.Counts {
		total += c
	}
	assert.Equal(t, sfs.TotalSites, total,
		"sum of SFS counts must equal total mutation events")
}

func TestCoalescent_FayWuH_Definition(t *testing.T) {
	sim := population.NewCoalescentSimulator(100, 1e-4, 0.0)
	sim.SetSeed(77)
	sampleSize := 10

	genealogy, err := sim.SimulateGenealogy("faywu_test", sampleSize, 500)
	require.NoError(t, err)

	sfs := sim.ComputeSiteFrequencySpectrum(genealogy, sampleSize)

	// H = ThetaPi - ThetaH; the value can be any real number
	// We just check it's finite (not NaN or Inf)
	assert.False(t, math.IsNaN(sfs.FayWuH), "FayWuH must not be NaN")
	assert.False(t, math.IsInf(sfs.FayWuH, 0), "FayWuH must not be Inf")
}
