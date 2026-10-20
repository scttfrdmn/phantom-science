//go:build !release

package population_test

import (
	"reflect"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	coalescent "github.com/scttfrdmn/phantom-science/internal/genomics/population"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

const (
	testStart  = 42904241 // 0-based start of a 50 kb window
	testLength = 50000
)

func TestSimulateHaplotypes_Shape(t *testing.T) {
	res, err := population.SimulateHaplotypes("CEU", "chr1", testStart, testLength, 6, 1)
	require.NoError(t, err)
	assert.Equal(t, population.Population("CEU"), res.Population)
	assert.Equal(t, "CEU", res.Preset)
	assert.Equal(t, 6, res.NumSamples)
	assert.Equal(t, len(res.Variants), res.NumVariants)
	assert.Equal(t, res.NumVariants, res.Statistics.SegregatingSites)
	require.Len(t, res.HaplotypeMatrix, 6)
	require.NotEmpty(t, res.Variants, "a 50 kb CEU sample of 6 haplotypes has segregating sites")
	assert.InDelta(t, 100*coalescent.PresetRecombinationRate*testLength, res.Statistics.GeneticMapLengthCM, 1e-12)

	prev := -1
	for i, v := range res.Variants {
		assert.GreaterOrEqual(t, v.Position, testStart)
		assert.Less(t, v.Position, testStart+testLength)
		assert.Greater(t, v.Position, prev, "positions strictly increase (infinite sites)")
		prev = v.Position
		require.Len(t, v.Alleles, 2)
		assert.NotEqual(t, v.Alleles[0], v.Alleles[1])
		assert.Equal(t, v.Alleles[0], v.AncestralState)
		carriers := 0
		for h := range res.HaplotypeMatrix {
			x := res.HaplotypeMatrix[h][i]
			require.True(t, x == 0 || x == 1)
			carriers += x
		}
		assert.Greater(t, carriers, 0, "segregating: some haplotype carries ALT")
		assert.Less(t, carriers, 6, "segregating: some haplotype carries REF")
		assert.InDelta(t, float64(carriers)/6, v.Frequency, 1e-12)
	}
}

func TestSimulateHaplotypes_Deterministic(t *testing.T) {
	a, err := population.SimulateHaplotypes("YRI", "chr7", 1000, 20000, 6, 42)
	require.NoError(t, err)
	b, err := population.SimulateHaplotypes("YRI", "chr7", 1000, 20000, 6, 42)
	require.NoError(t, err)
	assert.True(t, reflect.DeepEqual(a, b), "same seed must give identical haplotypes")

	differ := false
	for seed := uint64(43); seed < 48 && !differ; seed++ {
		c, err := population.SimulateHaplotypes("YRI", "chr7", 1000, 20000, 6, seed)
		require.NoError(t, err)
		differ = !reflect.DeepEqual(a.HaplotypeMatrix, c.HaplotypeMatrix) || !reflect.DeepEqual(a.Variants, c.Variants)
	}
	assert.True(t, differ, "different seeds should give different samples")
}

// Averaged over seeds, YRI carries more segregating sites than CEU, as in
// 1000 Genomes and the calibrated presets.
func TestSimulateHaplotypes_YRIMoreDiverseThanCEU(t *testing.T) {
	const n = 150
	var yri, ceu int
	for seed := uint64(1); seed <= n; seed++ {
		y, err := population.SimulateHaplotypes("YRI", "chr1", 0, 20000, 6, seed)
		require.NoError(t, err)
		c, err := population.SimulateHaplotypes("CEU", "chr1", 0, 20000, 6, seed)
		require.NoError(t, err)
		yri += y.NumVariants
		ceu += c.NumVariants
	}
	assert.Greater(t, yri, ceu, "mean S: YRI %.1f, CEU %.1f", float64(yri)/n, float64(ceu)/n)
}

func TestSimulateHaplotypes_Populations(t *testing.T) {
	for _, pop := range []string{"YRI", "CEU", "CHB", "MXL"} {
		res, err := population.SimulateHaplotypes(population.Population(pop), "chr1", 0, 10000, 6, 3)
		require.NoError(t, err, pop)
		assert.Equal(t, pop, res.Preset)
	}
	// 1000 Genomes aliases are simulated with the closest preset.
	res, err := population.SimulateHaplotypes("GBR", "chr1", 0, 10000, 6, 3)
	require.NoError(t, err)
	assert.Equal(t, "CEU", res.Preset)
	assert.Equal(t, population.Population("GBR"), res.Population)

	supported := population.SupportedPopulations()
	for _, want := range []string{"CEU", "CHB", "MXL", "YRI", "GBR→CEU", "JPT→CHB"} {
		assert.Contains(t, supported, want)
	}
}

func TestSimulateHaplotypes_Errors(t *testing.T) {
	_, err := population.SimulateHaplotypes("FIN", "chr1", 0, 10000, 6, 1)
	require.Error(t, err, "FIN has no preset")
	assert.Contains(t, err.Error(), `"FIN"`)
	assert.Contains(t, err.Error(), "CEU")
	assert.Contains(t, err.Error(), "YRI")

	_, err = population.SimulateHaplotypes("", "chr1", 0, 10000, 6, 1)
	assert.Error(t, err)
	_, err = population.SimulateHaplotypes("CEU", "chr1", 0, 0, 6, 1)
	assert.Error(t, err)
	_, err = population.SimulateHaplotypes("CEU", "chr1", 0, 10000, 1, 1)
	assert.Error(t, err)
	_, err = population.SimulateHaplotypes("CEU", "chr1", -1, 10000, 6, 1)
	assert.Error(t, err)
}

// The simulator is pure Go: it runs with nothing on PATH.
func TestSimulateHaplotypes_NoExternalProcess(t *testing.T) {
	t.Setenv("PATH", t.TempDir())
	_, err := population.SimulateHaplotypes("CHB", "chr12", 102815492, testLength, 6, 9)
	require.NoError(t, err)
}
