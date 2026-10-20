// Package population simulates neutral population-genetic variation with the
// coalescent.
//
// SimulateJointPanel runs a structured coalescent over a population tree
// (splits and piecewise-exponential size histories; no migration, admixture
// or within-locus recombination) on unlinked loci, so that every population
// is scored at the same sites. PopulationPanel and NewJointPanels wrap it
// using the demographic presets for YRI, CEU, CHB and MXL in presets.go;
// other population names fall back to a constant-size model.
package population

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sync"
)

// PopulationPanel represents a panel of population genetic data
type PopulationPanel struct {
	PopulationName string
	// SampleSize is the number of diploid individuals; panels are simulated
	// from 2·SampleSize haplotypes.
	SampleSize int
	// Ne is the effective size of this population's own branch (for presets,
	// its pre-growth Ne) and MutationRate the per-site, per-generation
	// mutation rate used by LoadFromSimulation.
	Ne           int
	MutationRate float64
	// AlleleFrequencies maps SiteKey(position) to the derived-allele
	// frequency. For panels built together by NewJointPanels every panel has
	// the same keys (sites monomorphic in a population have frequency 0 or 1).
	AlleleFrequencies map[string]float64
	Demographics      map[string]interface{}
	// Variants lists the sites segregating in this panel, with genotypes
	// formed from the simulated haplotypes (individual k = haplotypes 2k, 2k+1).
	Variants []Variant
	preset   string             // preset population used for the demographic model ("" = none)
	events   []DemographicEvent // extra events on this branch (SetDemographicEvent)
	rng      *rand.Rand         // seeds LoadFromSimulation; fixed default, see SetSeed
	mu       sync.RWMutex
}

// Variant represents a genetic variant in the population panel
type Variant struct {
	Position      int
	Chromosome    string
	ReferenceBase byte
	AlternateBase byte
	Frequency     float64
	Genotypes     []int // 0=homozygous ref, 1=heterozygous, 2=homozygous alt
	Quality       float64
}

// Default single-population model for populations without a preset.
const (
	defaultPanelNe           = 10000
	defaultPanelMutationRate = 1e-8
)

// NewPopulationPanel creates a population panel. Populations with a
// demographic preset (YRI, CEU, CHB, MXL and the aliases accepted by
// ResolvePreset) are simulated under the joint preset model in presets.go.
// Any other name falls back to a constant-size model (Ne 10000, μ 1e-8) with
// no demographic events; use HasPreset to detect this, or
// NewPresetPopulationPanel to get an error instead. The panel's random source
// starts from a fixed seed (as SetSeed(0)), so LoadFromSimulation is
// reproducible; call SetSeed to choose another.
func NewPopulationPanel(populationName string, sampleSize int) *PopulationPanel {
	pp := &PopulationPanel{
		PopulationName:    populationName,
		SampleSize:        sampleSize,
		AlleleFrequencies: make(map[string]float64),
		Demographics:      make(map[string]interface{}),
		Variants:          make([]Variant, 0),
		rng:               rand.New(rand.NewPCG(0, 0)), // #nosec G404 -- simulation, not crypto
	}
	spec, err := LookupPreset(populationName)
	if err != nil {
		pp.Ne, pp.MutationRate = defaultPanelNe, defaultPanelMutationRate
		pp.Demographics["model"] = "default constant-size (no preset)"
		pp.Demographics["effective_population_size"] = float64(defaultPanelNe)
		pp.Demographics["mutation_rate"] = defaultPanelMutationRate
		return pp
	}
	pp.preset = spec.Name
	pp.Ne, pp.MutationRate = int(spec.Ne), PresetMutationRate
	pp.Demographics = presetDemographics(spec)
	return pp
}

// NewPresetPopulationPanel is like NewPopulationPanel but returns an error for
// populations that have no demographic preset.
func NewPresetPopulationPanel(populationName string, sampleSize int) (*PopulationPanel, error) {
	if _, err := LookupPreset(populationName); err != nil {
		return nil, err
	}
	return NewPopulationPanel(populationName, sampleSize), nil
}

// HasPreset reports whether the panel is simulated under a demographic preset.
func (pp *PopulationPanel) HasPreset() bool { return pp.preset != "" }

// PresetName returns the preset population modelling this panel ("" if none).
// For aliases such as "LWK" this is the proxy preset ("YRI").
func (pp *PopulationPanel) PresetName() string { return pp.preset }

func presetDemographics(spec PopSpec) map[string]interface{} {
	return map[string]interface{}{
		"model":                           "preset",
		"preset":                          spec.Name,
		"effective_population_size":       spec.Ne,
		"present_population_size":         spec.presentSize(),
		"growth_rate":                     spec.GrowthRate,
		"growth_start_generations":        spec.GrowthStart,
		"bottleneck_start_generations":    spec.BottleneckStart,
		"bottleneck_duration_generations": spec.BottleneckDuration,
		"bottleneck_size":                 spec.BottleneckNe,
		"split_from":                      spec.SplitFrom,
		"split_time_generations":          spec.SplitTime,
		"mutation_rate":                   PresetMutationRate,
	}
}

// LoadFromSimulation simulates regionLength bp for this population and
// replaces the panel's Variants and AlleleFrequencies. The region is split
// into unlinked loci of DefaultLocusLength bp (one genealogy each), and both
// Variants and AlleleFrequencies come from that single simulation. The
// simulation seed is drawn from the panel's random source (see SetSeed), so
// repeated calls give different panels and the sequence of panels is
// reproducible.
//
// For preset populations the whole preset population tree is used (so the
// panel's lineages pass through its ancestral populations); the panel's Ne,
// MutationRate and extra events apply to this population's branch. Other
// populations use a single constant-size population with those parameters.
func (pp *PopulationPanel) LoadFromSimulation(regionLength int) error {
	pp.mu.Lock()
	defer pp.mu.Unlock()

	if regionLength <= 0 {
		return fmt.Errorf("region length must be positive")
	}
	if pp.SampleSize < 1 {
		return fmt.Errorf("sample size must be at least 1 individual")
	}
	d, label, err := pp.singlePopulationDemography()
	if err != nil {
		return err
	}
	locusLen := DefaultLocusLength
	if regionLength < locusLen {
		locusLen = regionLength
	}
	d.LocusLength = locusLen
	nLoci := (regionLength + locusLen - 1) / locusLen

	jp, err := SimulateJointPanel(d, map[string]int{label: 2 * pp.SampleSize}, nLoci, pp.rng.Uint64())
	if err != nil {
		return fmt.Errorf("failed to simulate population %s: %w", pp.PopulationName, err)
	}
	pp.loadFromJoint(jp, label)
	return nil
}

// singlePopulationDemography returns the demography used to simulate this
// panel on its own and the name of the sampled population in it.
func (pp *PopulationPanel) singlePopulationDemography() (Demography, string, error) {
	ne, mu, events := pp.Ne, pp.MutationRate, append([]DemographicEvent(nil), pp.events...)
	if ne < 1 {
		return Demography{}, "", fmt.Errorf("effective population size must be positive")
	}
	if pp.preset != "" {
		d := PresetDemography()
		for i := range d.Pops {
			if d.Pops[i].Name == pp.preset {
				d.Pops[i].Ne = float64(ne)
				d.Pops[i].Events = events
			}
		}
		d.Mu = mu
		return d, pp.preset, nil
	}
	label := pp.PopulationName
	if label == "" {
		label = "sample"
	}
	return Demography{
		Pops: []PopSpec{{Name: label, Ne: float64(ne), Events: events}},
		Mu:   mu,
	}, label, nil
}

// loadFromJoint fills Variants and AlleleFrequencies from population label of
// a joint panel. The caller must hold pp.mu (or own pp exclusively).
func (pp *PopulationPanel) loadFromJoint(jp *JointPanel, label string) {
	i := jp.PopIndex(label)
	n := jp.Haplotypes[i]
	freqs := make(map[string]float64, len(jp.Sites))
	variants := make([]Variant, 0, len(jp.Sites))
	for _, s := range jp.Sites {
		c := s.DerivedCount[i]
		f := float64(c) / float64(n)
		freqs[SiteKey(s.Position)] = f
		if c == 0 || c == n {
			continue
		}
		variants = append(variants, Variant{
			Position:      s.Position,
			Chromosome:    "chr1", // Default chromosome
			ReferenceBase: s.Ancestral,
			AlternateBase: s.Derived,
			Frequency:     f,
			Genotypes:     s.Genotypes[i],
			Quality:       30.0, // Default quality score
		})
	}
	pp.AlleleFrequencies = freqs
	pp.Variants = variants
}

// NewJointPanels simulates the named populations together under the preset
// demography (one genealogy per locus spanning all of them) and returns one
// panel per name, with sampleSize diploid individuals each, plus the joint
// sample. All panels share the same site keys, which is what cross-population
// statistics such as F_ST need. Names are resolved with ResolvePreset; unknown
// names, or two names resolving to the same preset, are errors.
func NewJointPanels(populations []string, sampleSize, nLoci int, seed uint64) (map[string]*PopulationPanel, *JointPanel, error) {
	if sampleSize < 1 {
		return nil, nil, fmt.Errorf("sample size must be at least 1 individual")
	}
	samples := make(map[string]int, len(populations))
	resolved := make(map[string]string, len(populations))
	for _, name := range populations {
		preset, ok := ResolvePreset(name)
		if !ok {
			return nil, nil, fmt.Errorf("no demographic preset for population %q (modelled: %v)", name, PresetNames())
		}
		if _, dup := samples[preset]; dup {
			return nil, nil, fmt.Errorf("population %q resolves to preset %s more than once", name, preset)
		}
		samples[preset] = 2 * sampleSize
		resolved[name] = preset
	}
	jp, err := SimulateJointPanel(PresetDemography(), samples, nLoci, seed)
	if err != nil {
		return nil, nil, err
	}
	panels := make(map[string]*PopulationPanel, len(populations))
	for _, name := range populations {
		pp := NewPopulationPanel(name, sampleSize)
		pp.loadFromJoint(jp, resolved[name])
		pp.Demographics["joint_populations"] = append([]string(nil), jp.Populations...)
		panels[name] = pp
	}
	return panels, jp, nil
}

// GetVariants returns all variants in the panel
func (pp *PopulationPanel) GetVariants() []Variant {
	pp.mu.RLock()
	defer pp.mu.RUnlock()
	return pp.Variants
}

// GetAlleleFrequencies returns allele frequencies normalized as a probability distribution
func (pp *PopulationPanel) GetAlleleFrequencies() map[string]float64 {
	pp.mu.RLock()
	defer pp.mu.RUnlock()

	// If we have cached normalized frequencies, return them
	if len(pp.AlleleFrequencies) > 0 {
		return pp.AlleleFrequencies
	}

	// Generate normalized frequencies from variants
	if len(pp.Variants) > 0 {
		frequencies := make(map[string]float64)
		totalWeight := 0.0

		for i, variant := range pp.Variants {
			siteID := fmt.Sprintf("site_%d", i)
			// Use minor allele frequency as the weight, but ensure it's positive
			weight := math.Min(variant.Frequency, 1.0-variant.Frequency)
			if weight <= 0 {
				weight = 0.001 // Small positive value for monomorphic sites
			}
			frequencies[siteID] = weight
			totalWeight += weight
		}

		// Normalize to make a proper probability distribution (sum = 1.0)
		if totalWeight > 0 {
			for siteID := range frequencies {
				frequencies[siteID] /= totalWeight
			}
		}

		return frequencies
	}

	return make(map[string]float64)
}

// GetDemographics returns demographic information
func (pp *PopulationPanel) GetDemographics() map[string]interface{} {
	pp.mu.RLock()
	defer pp.mu.RUnlock()
	return pp.Demographics
}

// SetDemographicEvent adds a demographic event to the population
func (pp *PopulationPanel) SetDemographicEvent(event DemographicEvent) {
	pp.mu.Lock()
	defer pp.mu.Unlock()
	pp.events = sortedEvents(append(pp.events, event))
}

// SetSeed reseeds the panel's random source, which seeds LoadFromSimulation.
func (pp *PopulationPanel) SetSeed(seed int64) {
	pp.mu.Lock()
	defer pp.mu.Unlock()
	pp.rng = rand.New(rand.NewPCG(uint64(seed), 0)) // #nosec G115 G404
}

// GetPopulationStats returns population statistics
func (pp *PopulationPanel) GetPopulationStats() map[string]interface{} {
	pp.mu.RLock()
	defer pp.mu.RUnlock()

	// Calculate heterozygosity from variants
	heterozygosity := 0.0
	if len(pp.Variants) > 0 {
		totalHeterozygosity := 0.0
		for _, variant := range pp.Variants {
			// Heterozygosity = 2 * p * (1-p)
			p := variant.Frequency
			totalHeterozygosity += float64(2 * p * (1 - p))
		}
		heterozygosity = totalHeterozygosity / float64(len(pp.Variants))
	}

	stats := map[string]interface{}{
		"population_name": pp.PopulationName,
		"sample_size":     pp.SampleSize,
		"variants":        len(pp.Variants),
		"allele_sites":    len(pp.AlleleFrequencies),
		"heterozygosity":  heterozygosity,
	}

	stats["effective_population_size"] = pp.Ne
	stats["mutation_rate"] = pp.MutationRate
	stats["demographic_events"] = len(pp.events)

	return stats
}
