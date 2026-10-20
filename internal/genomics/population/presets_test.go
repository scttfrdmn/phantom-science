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

// presetSample simulates YRI, CEU and CHB jointly (50 diploids each).
func presetSample(t *testing.T, nLoci int) *population.JointPanel {
	t.Helper()
	jp, err := population.SimulateJointPanel(population.PresetDemography(),
		map[string]int{"YRI": 100, "CEU": 100, "CHB": 100}, nLoci, 20261018)
	require.NoError(t, err)
	return jp
}

func presetLoci(t *testing.T) int {
	if testing.Short() {
		return 600
	}
	return 2000
}

// 1000 Genomes Phase 3: nucleotide diversity and rare-variant counts are
// highest in African populations; CEU and CHB are similar. Differences are
// tested per locus (paired, the populations share genealogies) at 4σ.
func TestPresets_DiversityOrdering(t *testing.T) {
	nLoci := presetLoci(t)
	jp := presetSample(t, nLoci)
	bp := float64(jp.LocusLength)

	pi := map[string][]float64{}
	for _, p := range []string{"YRI", "CEU", "CHB"} {
		pi[p] = perLocusPi(jp, p)
		m, se := meanSE(pi[p])
		t.Logf("%s: π = %.3e ± %.1e per bp", p, m/bp, se/bp)
	}
	for _, p := range []string{"CEU", "CHB"} {
		m, se := meanSE(diff(pi["YRI"], pi[p]))
		assert.Greaterf(t, m-4*se, 0.0, "π(YRI) − π(%s) = %.3f ± %.3f per locus, want > 0 at 4σ", p, m, se)

		mr, ser := meanSE(diff(perLocusRare(jp, "YRI"), perLocusRare(jp, p)))
		assert.Greaterf(t, mr-4*ser, 0.0, "rare variants YRI − %s = %.3f ± %.3f per locus, want > 0", p, mr, ser)
	}
	// CHB ≈ CEU: the relative difference is small (within 10%).
	mC, _ := meanSE(pi["CEU"])
	m, se := meanSE(diff(pi["CHB"], pi["CEU"]))
	t.Logf("π(CHB) − π(CEU) = %.3f ± %.3f per locus (%.1f%%)", m, se, 100*m/mC)
	assert.Less(t, math.Abs(m)/mC, 0.10, "CEU and CHB diversity should be similar")
}

// Weir–Cockerham F_ST between preset populations must fall in the literature
// ranges recorded in internal/population/fst_analysis.py (removed after v0.9.0), and continental
// (African vs non-African) differentiation must exceed CEU–CHB.
func TestPresets_FSTLiteratureRanges(t *testing.T) {
	nLoci := presetLoci(t)
	jp := presetSample(t, nLoci)

	fst := func(a, b string) (float64, float64) {
		return jackknife(nLoci, 50, func(keep func(int) bool) float64 { return wcFST(jp, a, b, keep) })
	}
	ranges := []struct {
		a, b   string
		lo, hi float64
	}{
		{"YRI", "CEU", 0.10, 0.16},
		{"YRI", "CHB", 0.09, 0.15},
		{"CEU", "CHB", 0.06, 0.12},
	}
	for _, r := range ranges {
		f, se := fst(r.a, r.b)
		t.Logf("F_ST %s–%s = %.4f ± %.4f (literature %.2f–%.2f)", r.a, r.b, f, se, r.lo, r.hi)
		assert.GreaterOrEqualf(t, f, r.lo, "%s–%s F_ST below literature range", r.a, r.b)
		assert.LessOrEqualf(t, f, r.hi, "%s–%s F_ST above literature range", r.a, r.b)
	}

	for _, p := range []string{"CEU", "CHB"} {
		d, se := jackknife(nLoci, 50, func(keep func(int) bool) float64 {
			return wcFST(jp, "YRI", p, keep) - wcFST(jp, "CEU", "CHB", keep)
		})
		assert.Greaterf(t, d-4*se, 0.0, "F_ST(YRI–%s) − F_ST(CEU–CHB) = %.4f ± %.4f, want > 0 at 4σ", p, d, se)
	}
}

// Two samples of the same population must give F_ST ≈ 0. (The pre-fix test
// reported ≈ 0.33 for any pair of panels, including identical populations.)
func TestPresets_SamePopulationNullFST(t *testing.T) {
	d := population.PresetDemography()
	ceu, ok := d.Pop("CEU")
	require.True(t, ok)
	// CEU2 merges into CEU at time 0: one panmictic population, two samples.
	d.Pops = append(d.Pops, population.PopSpec{Name: "CEU2", Ne: ceu.Ne, SplitFrom: "CEU", SplitTime: 0})
	jp, err := population.SimulateJointPanel(d, map[string]int{"CEU": 100, "CEU2": 100}, presetLoci(t), 5)
	require.NoError(t, err)
	f := wcFST(jp, "CEU", "CEU2", func(int) bool { return true })
	t.Logf("same-population WC F_ST = %.4f", f)
	assert.Less(t, f, 0.02)
}

func TestPresets_LookupAndAliases(t *testing.T) {
	assert.Equal(t, []string{"CEU", "CHB", "MXL", "YRI"}, population.PresetNames())
	for name, want := range map[string]string{"YRI": "YRI", "LWK": "YRI", "EUR": "CEU", "GBR": "CEU", "JPT": "CHB", "PUR": "MXL"} {
		got, ok := population.ResolvePreset(name)
		assert.True(t, ok, name)
		assert.Equal(t, want, got, name)
	}
	_, ok := population.ResolvePreset("UNKNOWN")
	assert.False(t, ok)
	_, err := population.LookupPreset("UNKNOWN")
	assert.Error(t, err)

	spec, err := population.LookupPreset("CEU")
	require.NoError(t, err)
	assert.Equal(t, 11418.0, spec.Ne)
	assert.Equal(t, "YRI", spec.SplitFrom)
}

// Issue #5: the population name must select the demographic model.
func TestNewPopulationPanel_UsesPreset(t *testing.T) {
	yri := population.NewPopulationPanel("YRI", 10)
	ceu := population.NewPopulationPanel("CEU", 10)
	assert.True(t, yri.HasPreset())
	assert.Equal(t, "YRI", yri.PresetName())
	assert.Equal(t, 12300, yri.Ne)
	assert.Equal(t, 11418, ceu.Ne)
	assert.Equal(t, "YRI", ceu.GetDemographics()["split_from"])

	lwk := population.NewPopulationPanel("LWK", 10)
	assert.Equal(t, "YRI", lwk.PresetName())
	assert.Equal(t, "LWK", lwk.PopulationName)

	unknown := population.NewPopulationPanel("XYZ", 10)
	assert.False(t, unknown.HasPreset())
	_, err := population.NewPresetPopulationPanel("XYZ", 10)
	assert.Error(t, err)
	p, err := population.NewPresetPopulationPanel("CHB", 10)
	require.NoError(t, err)
	assert.Equal(t, "CHB", p.PresetName())
}

// Separately loaded YRI and CEU panels come from different models: YRI has
// more segregating sites per locus (independent simulations, 4σ).
func TestPopulationPanel_LoadFromSimulation_PopulationSpecific(t *testing.T) {
	const region = 10_000_000 // 1000 loci
	perLocusS := func(name string, seed int64) []float64 {
		p := population.NewPopulationPanel(name, 50)
		p.SetSeed(seed)
		require.NoError(t, p.LoadFromSimulation(region))
		out := make([]float64, region/population.DefaultLocusLength)
		for _, v := range p.GetVariants() {
			out[v.Position/population.DefaultLocusLength]++
		}
		return out
	}
	mY, seY := meanSE(perLocusS("YRI", 1))
	mC, seC := meanSE(perLocusS("CEU", 2))
	t.Logf("S per locus: YRI %.2f ± %.2f, CEU %.2f ± %.2f", mY, seY, mC, seC)
	se := math.Sqrt(seY*seY + seC*seC)
	assert.Greater(t, mY-mC-4*se, 0.0)
}

// Variants and AlleleFrequencies come from one simulation, and genotypes are
// built from the simulated haplotypes.
func TestPopulationPanel_LoadFromSimulation_Consistent(t *testing.T) {
	p := population.NewPopulationPanel("CHB", 20)
	p.SetSeed(9)
	require.NoError(t, p.LoadFromSimulation(200_000))
	variants := p.GetVariants()
	freqs := p.GetAlleleFrequencies()
	require.NotEmpty(t, variants)
	assert.Len(t, freqs, len(variants), "single-population panel: every site segregates")
	for _, v := range variants {
		f, ok := freqs[population.SiteKey(v.Position)]
		require.True(t, ok, "variant %d missing from AlleleFrequencies", v.Position)
		assert.Equal(t, v.Frequency, f)
		require.Len(t, v.Genotypes, 20)
		dosage := 0
		for _, g := range v.Genotypes {
			assert.Contains(t, []int{0, 1, 2}, g)
			dosage += g
		}
		assert.InDelta(t, v.Frequency*40, float64(dosage), 1e-9)
	}

	// Same seed → identical panel; reloading replaces rather than appends.
	q := population.NewPopulationPanel("CHB", 20)
	q.SetSeed(9)
	require.NoError(t, q.LoadFromSimulation(200_000))
	assert.True(t, reflect.DeepEqual(variants, q.GetVariants()))
	require.NoError(t, q.LoadFromSimulation(200_000))
	assert.Len(t, q.GetAlleleFrequencies(), len(q.GetVariants()), "reload must replace, not append")
}

func TestPopulationPanel_LoadFromSimulation_Errors(t *testing.T) {
	p := population.NewPopulationPanel("CEU", 0)
	assert.Error(t, p.LoadFromSimulation(1000))
	p = population.NewPopulationPanel("CEU", 5)
	assert.Error(t, p.LoadFromSimulation(0))
	p = population.NewPopulationPanel("XYZ", 5)
	p.Ne = 0
	assert.Error(t, p.LoadFromSimulation(1000))
}

// Issue #6: panels built jointly share one site set, and the same seed gives
// the same panels.
func TestNewJointPanels(t *testing.T) {
	pops := []string{"YRI", "CEU", "CHB"}
	panels, jp, err := population.NewJointPanels(pops, 10, 30, 77)
	require.NoError(t, err)
	require.Len(t, panels, 3)
	keys := func(m map[string]float64) []string {
		out := make([]string, 0, len(m))
		for k := range m {
			out = append(out, k)
		}
		return out
	}
	ref := panels["YRI"].GetAlleleFrequencies()
	assert.Len(t, ref, len(jp.Sites))
	for _, p := range pops {
		assert.ElementsMatch(t, keys(ref), keys(panels[p].GetAlleleFrequencies()), "%s must share the joint site set", p)
		assert.Equal(t, jp.Frequencies(panels[p].PresetName()), panels[p].GetAlleleFrequencies())
	}

	again, _, err := population.NewJointPanels(pops, 10, 30, 77)
	require.NoError(t, err)
	for _, p := range pops {
		assert.Equal(t, panels[p].GetVariants(), again[p].GetVariants())
		assert.Equal(t, panels[p].GetAlleleFrequencies(), again[p].GetAlleleFrequencies())
	}

	_, _, err = population.NewJointPanels([]string{"CEU", "EUR"}, 10, 5, 1)
	assert.Error(t, err, "two names resolving to the same preset")
	_, _, err = population.NewJointPanels([]string{"XYZ"}, 10, 5, 1)
	assert.Error(t, err)
	_, _, err = population.NewJointPanels([]string{"YRI"}, 0, 5, 1)
	assert.Error(t, err)
}
