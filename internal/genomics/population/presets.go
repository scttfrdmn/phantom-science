package population

import (
	"fmt"
	"sort"
)

// Demographic presets for the 1000 Genomes populations modelled by this
// package.
//
// Per-population parameters (Ne, growth rate, bottleneck start and size,
// mutation rate) are copied from the static table in
// internal/population/bridge.go (GetPopulationParameters). bridge.go,
// demographics.py and fst_analysis.py, cited below as the sources of these
// values, were removed after v0.9.0; see them in the v0.9.0 tag. Note that
// internal/population/demographics.py used different values for MXL
// (Ne 7000, growth 0.008, bottleneck 500 gen / 800) and lists its YRI growth
// from the present back to the bottleneck, which as an msprime model would
// shrink YRI to ~0 by 2000 generations; the values below follow bridge.go.
//
// bridge.go does not give growth onsets, bottleneck durations, or a population
// tree, so these are added here:
//
//   - Growth onset: 204 generations ago for every population (Tennessen et
//     al. 2012, recent growth over the last 5.1 kya; YRI's 1.66%/gen rate in
//     bridge.go is theirs). Gravel et al. 2011 apply the EUR/EAS rates in
//     bridge.go (0.38%, 0.48%) from 920 generations, but starting that growth
//     from the bridge.go Ne (rather than Gravel's ~1000-individual founders)
//     gives CEU/CHB present sizes of 0.4–0.7 million and more segregating
//     sites than YRI, contradicting 1000 Genomes; the shorter window keeps
//     YRI the most diverse in both π and S.
//   - Bottleneck durations (bridge.go gives start and size only): 100
//     generations by default (the bottleneck length used in
//     internal/population/demographics.py), YRI 20, CEU 400 (the
//     out-of-Africa bottleneck, cf. Gravel et al. 2011: ~1100 generations at
//     N≈1900).
//   - Topology (no migration, no admixture): CEU splits from YRI (out of
//     Africa, 2800 generations = bridge.go YRI FoundingTimeGenerations), CHB
//     splits from CEU (1800, calibrated), MXL splits from CHB (600 = bridge.go
//     MXL FoundingTimeGenerations). MXL is modelled as a
//     Native-American-like branch; its recent European/African admixture is
//     not modelled.
//   - Split times and bottleneck durations are calibrated so the simulated
//     Weir–Cockerham F_ST falls inside the literature ranges recorded in
//     internal/population/fst_analysis.py (AFR–EUR 0.10–0.16, AFR–EAS
//     0.09–0.15, EUR–EAS 0.06–0.12) and π follows the 1000 Genomes Phase 3
//     ordering (YRI highest). See presets_test.go for the calibration checks.
//   - CHB's bridge.go bottleneck starts at 1800 generations, the CEU–CHB
//     split, and MXL's (1500) predates its 600-generation founding, so both
//     fall in ancestral branches (where the CEU out-of-Africa bottleneck
//     applies) and have no effect of their own.
//
// Measured with these values (2000 loci × 10 kb, 50 diploids per population;
// see presets_test.go): π ≈ 6.0e-4 (YRI), 5.8e-4 (CHB), 5.6e-4 (CEU) per
// bp; Weir–Cockerham F_ST ≈ 0.14 (YRI–CEU), 0.13–0.14 (YRI–CHB), 0.095
// (CEU–CHB), 0.03 (CHB–MXL; MXL is weakly differentiated because its
// admixture is not modelled).
const (
	// PresetMutationRate is the per-bp per-generation mutation rate (bridge.go).
	PresetMutationRate = 1.25e-8
	// PresetRecombinationRate is the per-bp per-generation recombination rate
	// (bridge.go). Loci are simulated without recombination; the value is
	// recorded on the panel's simulator for reference.
	PresetRecombinationRate = 1.13e-8

	// Calibrated split times (generations ago).
	presetSplitCEUfromYRI = 2800
	presetSplitCHBfromCEU = 1800
	presetSplitMXLfromCHB = 600

	// presetGrowthStart is the onset of recent exponential growth for every
	// preset population (Tennessen et al. 2012: 5.1 kya ≈ 204 generations).
	presetGrowthStart = 204
)

// presetPops returns the preset population tree.
func presetPops() []PopSpec {
	return []PopSpec{
		{
			Name: "YRI", Ne: 12300,
			GrowthRate: 0.0166, GrowthStart: presetGrowthStart,
			BottleneckStart: 2000, BottleneckDuration: 20, BottleneckNe: 1000,
		},
		{
			Name: "CEU", Ne: 11418,
			GrowthRate: 0.0038, GrowthStart: presetGrowthStart,
			BottleneckStart: 1600, BottleneckDuration: 400, BottleneckNe: 2000,
			SplitFrom: "YRI", SplitTime: presetSplitCEUfromYRI,
		},
		{
			Name: "CHB", Ne: 9025,
			GrowthRate: 0.0048, GrowthStart: presetGrowthStart,
			BottleneckStart: 1800, BottleneckDuration: 100, BottleneckNe: 1500,
			SplitFrom: "CEU", SplitTime: presetSplitCHBfromCEU,
		},
		{
			Name: "MXL", Ne: 8500,
			GrowthRate: 0.0045, GrowthStart: presetGrowthStart,
			BottleneckStart: 1500, BottleneckDuration: 100, BottleneckNe: 1200,
			SplitFrom: "CHB", SplitTime: presetSplitMXLfromCHB,
		},
	}
}

// PresetDemography returns the joint YRI/CEU/CHB/MXL demographic model.
func PresetDemography() Demography {
	return Demography{
		Pops:        presetPops(),
		Mu:          PresetMutationRate,
		LocusLength: DefaultLocusLength,
	}
}

// presetAliases maps other 1000 Genomes population and super-population codes
// to the closest modelled preset. These are proxies: e.g. LWK and GWD are
// simulated with the YRI model.
var presetAliases = map[string]string{
	"AFR": "YRI", "LWK": "YRI", "GWD": "YRI",
	"EUR": "CEU", "GBR": "CEU", "TSI": "CEU",
	"EAS": "CHB", "ASN": "CHB", "JPT": "CHB", "CHS": "CHB",
	"AMR": "MXL", "PUR": "MXL",
}

// PresetNames returns the names of the modelled preset populations.
func PresetNames() []string {
	pops := presetPops()
	names := make([]string, len(pops))
	for i, p := range pops {
		names[i] = p.Name
	}
	sort.Strings(names)
	return names
}

// ResolvePreset returns the preset population used to model name: the name
// itself for YRI/CEU/CHB/MXL, or the proxy preset for an alias (see
// presetAliases). ok is false for populations without a preset.
func ResolvePreset(name string) (preset string, ok bool) {
	for _, p := range presetPops() {
		if p.Name == name {
			return name, true
		}
	}
	if p, ok := presetAliases[name]; ok {
		return p, true
	}
	return "", false
}

// LookupPreset returns the preset PopSpec for a population name or alias.
func LookupPreset(name string) (PopSpec, error) {
	preset, ok := ResolvePreset(name)
	if !ok {
		return PopSpec{}, fmt.Errorf("no demographic preset for population %q (modelled: %v)", name, PresetNames())
	}
	spec, _ := PresetDemography().Pop(preset)
	return spec, nil
}
