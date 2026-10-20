//go:build !release

package population_test

import (
	"math"
	"reflect"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics/population"
)

func constantDemography(ne float64) population.Demography {
	return population.Demography{
		Pops:        []population.PopSpec{{Name: "A", Ne: ne}},
		Mu:          1.25e-8,
		LocusLength: 10000,
	}
}

// Constant N: E[T2] = 2N and E[π per locus] = θ = 4NμL (diploid Ne).
func TestStructured_ConstantN_PairwiseCoalescence(t *testing.T) {
	const ne, nLoci = 5000.0, 4000
	d := constantDemography(ne)
	jp, err := population.SimulateJointPanel(d, map[string]int{"A": 2}, nLoci, 1)
	require.NoError(t, err)

	mT, seT := meanSE(jp.LocusTMRCA)
	assert.InDeltaf(t, 2*ne, mT, 4*seT, "E[T2] = %.0f ± %.0f, want %.0f", mT, seT, 2*ne)

	theta := 4 * ne * d.Mu * float64(d.LocusLength)
	mPi, sePi := meanSE(perLocusPi(jp, "A"))
	assert.InDeltaf(t, theta, mPi, 4*sePi, "E[π] = %.3f ± %.3f per locus, want θ = %.3f", mPi, sePi, theta)
}

// Constant N, larger sample: Watterson's θ and π both estimate θ, and
// Tajima's D is ≈ 0. This is also the negative control for the growth test.
func TestStructured_ConstantN_WattersonMatchesPi(t *testing.T) {
	const ne, nLoci, n = 5000.0, 2000, 20
	d := constantDemography(ne)
	jp, err := population.SimulateJointPanel(d, map[string]int{"A": n}, nLoci, 2)
	require.NoError(t, err)
	sum, err := jp.Summary("A")
	require.NoError(t, err)

	a1 := 0.0
	for i := 1; i < n; i++ {
		a1 += 1 / float64(i)
	}
	theta := 4 * ne * d.Mu * float64(d.LocusLength)
	thetaW := make([]float64, nLoci)
	for l, s := range sum.PerLocusS {
		thetaW[l] = float64(s) / a1
	}
	mW, seW := meanSE(thetaW)
	mPi, sePi := meanSE(sum.PerLocusPi)
	assert.InDeltaf(t, theta, mW, 4*seW, "θ_W = %.3f ± %.3f, want %.3f", mW, seW, theta)
	assert.InDeltaf(t, theta, mPi, 4*sePi, "π = %.3f ± %.3f, want %.3f", mPi, sePi, theta)
	mDiff, seDiff := meanSE(diff(thetaW, sum.PerLocusPi))
	assert.InDeltaf(t, 0, mDiff, 4*seDiff, "θ_W − π = %.3f ± %.3f", mDiff, seDiff)

	mD, seD := meanSE(finite(sum.PerLocusD))
	assert.InDeltaf(t, 0, mD, 4*seD, "constant-N Tajima's D = %.3f ± %.3f, want ≈ 0", mD, seD)
}

// Exponential growth produces an excess of rare variants: Tajima's D < 0.
func TestStructured_ExponentialGrowth_NegativeTajimaD(t *testing.T) {
	const nLoci, n = 2000, 20
	d := population.Demography{
		Pops:        []population.PopSpec{{Name: "A", Ne: 5000, GrowthRate: 0.02, GrowthStart: 200}},
		Mu:          1.25e-8,
		LocusLength: 10000,
	}
	jp, err := population.SimulateJointPanel(d, map[string]int{"A": n}, nLoci, 3)
	require.NoError(t, err)
	sum, err := jp.Summary("A")
	require.NoError(t, err)

	mD, seD := meanSE(finite(sum.PerLocusD))
	t.Logf("growth: Tajima's D = %.3f ± %.3f", mD, seD)
	assert.Less(t, mD+4*seD, 0.0, "exponential growth must give Tajima's D < 0")
}

func twoPopSplit(ne float64, split int) population.Demography {
	return population.Demography{
		Pops: []population.PopSpec{
			{Name: "A", Ne: ne},
			{Name: "B", Ne: ne, SplitFrom: "A", SplitTime: split},
		},
		Mu:          1.25e-8,
		LocusLength: 10000,
	}
}

// Drift-only split of two equal-size populations: Hudson's F_ST has
// expectation 1 − T_within/T_between = T/(T + 2N), increases with T, and is
// ≈ 0 at T = 0.
func TestStructured_SplitFST_MatchesTheory(t *testing.T) {
	const ne, nLoci, n = 5000.0, 1500, 40
	prev := -1.0
	for i, split := range []int{0, 500, 2000, 5000} {
		jp, err := population.SimulateJointPanel(twoPopSplit(ne, split), map[string]int{"A": n, "B": n}, nLoci, uint64(10+i))
		require.NoError(t, err)

		hud, se := jackknife(nLoci, 50, func(keep func(int) bool) float64 { return hudsonFST(jp, "A", "B", keep) })
		want := float64(split) / (float64(split) + 2*ne)
		t.Logf("T=%d: Hudson F_ST = %.4f ± %.4f (theory %.4f)", split, hud, se, want)
		assert.InDeltaf(t, want, hud, 4*se, "T=%d: Hudson F_ST %.4f ± %.4f, theory %.4f", split, hud, se, want)
		assert.Greater(t, hud, prev, "F_ST must increase with split time")
		prev = hud

		if split == 0 {
			wc := wcFST(jp, "A", "B", func(int) bool { return true })
			assert.Less(t, wc, 0.01, "T=0 (one panmictic population): WC F_ST must be ≈ 0")
		}
	}
}

func TestSimulateJointPanel_Deterministic(t *testing.T) {
	d := population.PresetDemography()
	samples := map[string]int{"YRI": 10, "CEU": 10}
	a, err := population.SimulateJointPanel(d, samples, 20, 42)
	require.NoError(t, err)
	b, err := population.SimulateJointPanel(d, samples, 20, 42)
	require.NoError(t, err)
	assert.True(t, reflect.DeepEqual(a, b), "same seed must give identical panels")

	c, err := population.SimulateJointPanel(d, samples, 20, 43)
	require.NoError(t, err)
	assert.False(t, reflect.DeepEqual(a.LocusTMRCA, c.LocusTMRCA), "different seeds should differ")
}

func TestSimulateJointPanel_SiteInvariants(t *testing.T) {
	jp, err := population.SimulateJointPanel(population.PresetDemography(),
		map[string]int{"YRI": 20, "CEU": 20, "CHB": 20}, 50, 7)
	require.NoError(t, err)
	require.NotEmpty(t, jp.Sites)
	seen := map[int]bool{}
	for _, s := range jp.Sites {
		assert.False(t, seen[s.Position], "infinite sites: positions are unique")
		seen[s.Position] = true
		assert.Equal(t, s.Locus, s.Position/jp.LocusLength)
		assert.NotEqual(t, s.Ancestral, s.Derived)
		total := 0
		for i, c := range s.DerivedCount {
			total += c
			assert.LessOrEqual(t, c, jp.Haplotypes[i])
			g := 0
			for _, x := range s.Genotypes[i] {
				g += x
			}
			assert.Equal(t, c, g, "genotype dosages must sum to the derived count")
			require.Len(t, s.Carriers[i], jp.Haplotypes[i])
			carried := 0
			for h, x := range s.Carriers[i] {
				require.True(t, x == 0 || x == 1)
				carried += int(x)
				if s.Genotypes[i] != nil && h%2 == 1 {
					assert.Equal(t, s.Genotypes[i][h/2], int(s.Carriers[i][h-1])+int(x),
						"haplotypes 2k and 2k+1 make up genotype k")
				}
			}
			assert.Equal(t, c, carried, "carrier rows must sum to the derived count")
		}
		// A mutation is below the root, so it is never carried by every
		// sampled haplotype, and never by none.
		assert.Greater(t, total, 0)
		assert.Less(t, total, 60)
	}
	for l, tm := range jp.LocusTMRCA {
		assert.Greater(t, jp.LocusTreeLength[l], tm, "tree length exceeds TMRCA for n > 2")
	}
}

// Two sites on the same genealogy are either nested, disjoint or identical
// in their carrier sets (no recombination within a locus), so no pair of
// sites shows all four gametes.
func TestSimulateJointPanel_CarriersFourGameteFree(t *testing.T) {
	jp, err := population.SimulateJointPanel(population.PresetDemography(),
		map[string]int{"CEU": 12}, 3, 11)
	require.NoError(t, err)
	byLocus := map[int][]population.JointSite{}
	for _, s := range jp.Sites {
		byLocus[s.Locus] = append(byLocus[s.Locus], s)
	}
	pairs := 0
	for _, sites := range byLocus {
		for a := 0; a < len(sites); a++ {
			for b := a + 1; b < len(sites); b++ {
				var gametes [2][2]bool
				for h := range sites[a].Carriers[0] {
					gametes[sites[a].Carriers[0][h]][sites[b].Carriers[0][h]] = true
				}
				pairs++
				assert.False(t, gametes[0][0] && gametes[0][1] && gametes[1][0] && gametes[1][1],
					"sites %d and %d show all four gametes", sites[a].Position, sites[b].Position)
			}
		}
	}
	require.Greater(t, pairs, 0)
}

func TestDemographyValidate(t *testing.T) {
	ok := population.PresetDemography()
	require.NoError(t, ok.Validate())

	cases := map[string]population.Demography{
		"no pops":        {Mu: 1e-8},
		"two roots":      {Pops: []population.PopSpec{{Name: "A", Ne: 10}, {Name: "B", Ne: 10}}},
		"unknown parent": {Pops: []population.PopSpec{{Name: "A", Ne: 10}, {Name: "B", Ne: 10, SplitFrom: "Z", SplitTime: 5}}},
		"duplicate":      {Pops: []population.PopSpec{{Name: "A", Ne: 10}, {Name: "A", Ne: 10}}},
		"bad Ne":         {Pops: []population.PopSpec{{Name: "A", Ne: 0}}},
		"child older than parent split": {Pops: []population.PopSpec{
			{Name: "A", Ne: 10},
			{Name: "B", Ne: 10, SplitFrom: "A", SplitTime: 100},
			{Name: "C", Ne: 10, SplitFrom: "B", SplitTime: 200},
		}},
	}
	for name, d := range cases {
		assert.Errorf(t, d.Validate(), "%s should be invalid", name)
	}

	_, err := population.SimulateJointPanel(ok, map[string]int{"XXX": 4}, 1, 1)
	assert.Error(t, err, "unknown sampled population")
	_, err = population.SimulateJointPanel(ok, map[string]int{"YRI": 1}, 1, 1)
	assert.Error(t, err, "fewer than 2 haplotypes")
	_, err = population.SimulateJointPanel(ok, map[string]int{"YRI": 4}, 0, 1)
	assert.Error(t, err, "zero loci")
}

func finite(xs []float64) []float64 {
	out := make([]float64, 0, len(xs))
	for _, x := range xs {
		if !math.IsNaN(x) {
			out = append(out, x)
		}
	}
	return out
}
