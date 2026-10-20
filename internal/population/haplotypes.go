// Package population simulates phased haplotypes for a small sample from one
// 1000 Genomes population, for trio generation in clinical scenarios
// (internal/genomics/clinical). It is pure Go: the haplotypes come from the
// structured coalescent and demographic presets in
// internal/genomics/population, with no Python or external process.
package population

import (
	"fmt"
	"sort"
	"strings"

	coalescent "github.com/scttfrdmn/phantom-science/internal/genomics/population"
)

// Population is a 1000 Genomes population or super-population code.
type Population string

// HaplotypeVariant is one segregating site of a simulated sample.
type HaplotypeVariant struct {
	Position       int      `json:"position"` // 0-based genomic coordinate
	Alleles        []string `json:"alleles"`  // [ancestral (REF), derived (ALT)]
	AncestralState string   `json:"ancestral_state"`
	Frequency      float64  `json:"frequency"` // derived-allele frequency in the sample
}

// HaplotypeStatistics summarises a simulated sample.
type HaplotypeStatistics struct {
	SegregatingSites int     `json:"segregating_sites"`
	NumHaplotypes    int     `json:"n_haplotypes"`
	DiversityPi      float64 `json:"diversity_pi"` // mean pairwise differences per bp
	// GeneticMapLengthCM is the genetic length of the simulated region at
	// the preset recombination rate (100 × rate × length).
	GeneticMapLengthCM float64 `json:"genetic_map_length_cm"`
}

// HaplotypeSimulationResult is a phased sample of NumSamples haplotypes over
// one region.
type HaplotypeSimulationResult struct {
	Population     Population          `json:"population"`
	Preset         string              `json:"preset"` // demographic preset that modelled Population
	Chromosome     string              `json:"chromosome"`
	RegionStart    int                 `json:"region_start"` // 0-based start of the simulated region
	SequenceLength int                 `json:"sequence_length"`
	NumVariants    int                 `json:"num_variants"`
	NumSamples     int                 `json:"num_samples"` // haplotypes (rows of HaplotypeMatrix)
	Variants       []HaplotypeVariant  `json:"variants"`
	Statistics     HaplotypeStatistics `json:"statistics"`
	// HaplotypeMatrix[h][i] is 1 when haplotype h carries the derived
	// (ALT) allele of Variants[i], else 0.
	HaplotypeMatrix [][]int `json:"haplotype_matrix"`
}

// thousandGenomesCodes are the 26 1000 Genomes Phase 3 populations, the five
// super-populations and the Phase 1 code ASN.
var thousandGenomesCodes = []string{
	"ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN",
	"GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL",
	"PEL", "PJL", "PUR", "STU", "TSI", "YRI",
	"AFR", "AMR", "EAS", "EUR", "SAS", "ASN",
}

// SupportedPopulations returns the population codes SimulateHaplotypes
// accepts, sorted, each with the preset that models it ("GBR→CEU"). Only
// YRI, CEU, CHB and MXL have their own demographic model; the other codes
// are proxies simulated with the closest preset.
func SupportedPopulations() []string {
	var out []string
	for _, code := range thousandGenomesCodes {
		if preset, ok := coalescent.ResolvePreset(code); ok {
			if preset == code {
				out = append(out, code)
			} else {
				out = append(out, code+"→"+preset)
			}
		}
	}
	sort.Strings(out)
	return out
}

// SimulateHaplotypes simulates numHaplotypes phased haplotypes of pop over a
// region of sequenceLength bp starting at the 0-based coordinate regionStart
// on chromosome. The sample is drawn from the joint YRI/CEU/CHB/MXL
// structured coalescent (coalescent.PresetDemography) as a single
// non-recombining locus covering the whole region, under the infinite-sites
// model; variant positions are regionStart plus the offset within the locus.
// The result is deterministic for a seed. A population without a preset is
// an error listing the supported codes.
func SimulateHaplotypes(pop Population, chromosome string, regionStart, sequenceLength, numHaplotypes int, seed uint64) (*HaplotypeSimulationResult, error) {
	preset, ok := coalescent.ResolvePreset(string(pop))
	if !ok {
		return nil, fmt.Errorf("haplotype simulation: no demographic model for population %q; supported: %s",
			pop, strings.Join(SupportedPopulations(), ", "))
	}
	if sequenceLength < 1 {
		return nil, fmt.Errorf("haplotype simulation: sequence length must be positive, got %d", sequenceLength)
	}
	if numHaplotypes < 2 {
		return nil, fmt.Errorf("haplotype simulation: need at least 2 haplotypes, got %d", numHaplotypes)
	}
	if regionStart < 0 {
		return nil, fmt.Errorf("haplotype simulation: negative region start %d", regionStart)
	}

	d := coalescent.PresetDemography()
	d.LocusLength = sequenceLength
	jp, err := coalescent.SimulateJointPanel(d, map[string]int{preset: numHaplotypes}, 1, seed)
	if err != nil {
		return nil, fmt.Errorf("haplotype simulation: %w", err)
	}
	summary, err := jp.Summary(preset)
	if err != nil {
		return nil, fmt.Errorf("haplotype simulation: %w", err)
	}

	res := &HaplotypeSimulationResult{
		Population:      pop,
		Preset:          preset,
		Chromosome:      chromosome,
		RegionStart:     regionStart,
		SequenceLength:  sequenceLength,
		NumVariants:     len(jp.Sites),
		NumSamples:      numHaplotypes,
		Variants:        make([]HaplotypeVariant, len(jp.Sites)),
		HaplotypeMatrix: make([][]int, numHaplotypes),
		Statistics: HaplotypeStatistics{
			SegregatingSites:   len(jp.Sites),
			NumHaplotypes:      numHaplotypes,
			DiversityPi:        summary.Pi,
			GeneticMapLengthCM: 100 * coalescent.PresetRecombinationRate * float64(sequenceLength),
		},
	}
	for h := range res.HaplotypeMatrix {
		res.HaplotypeMatrix[h] = make([]int, len(jp.Sites))
	}
	for i, s := range jp.Sites {
		anc, der := string(s.Ancestral), string(s.Derived)
		res.Variants[i] = HaplotypeVariant{
			Position:       regionStart + s.Position,
			Alleles:        []string{anc, der},
			AncestralState: anc,
			Frequency:      float64(s.DerivedCount[0]) / float64(numHaplotypes),
		}
		for h, c := range s.Carriers[0] {
			res.HaplotypeMatrix[h][i] = int(c)
		}
	}
	return res, nil
}
