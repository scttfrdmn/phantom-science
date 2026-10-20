//go:build !release

package clinical

import (
	"fmt"
	"math"
	"math/rand"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

// TrioHaplotypes holds the resolved haplotypes for all three family members.
// Each haplotype is a slice indexed by variant index in Variants,
// containing allele values (0=ref, 1=alt).
// DeNovoVariant is NOT in Variants; it is planted on ProbandH1 only.
type TrioHaplotypes struct {
	MotherH1, MotherH2   []int
	FatherH1, FatherH2   []int
	ProbandH1, ProbandH2 []int
	Variants             []population.HaplotypeVariant
	// DeNovoVariant is not in Variants; the generator plants it on ProbandH1 only.
	DeNovoVariant *genomics.VariantRecord
}

// GenerateTrioHaplotypes takes a phased HaplotypeSimulationResult with at
// least 6 haplotypes (population.SimulateHaplotypes), designates rows 0-3 as
// mother_h1, mother_h2, father_h1 and father_h2 (rows 4 and 5 are unused
// background), then simulates Mendelian transmission: proband_h1 is a
// recombinant of the mother's pair and proband_h2 of the father's pair. The
// de novo variant is recorded for ProbandH1 only. Background variants that
// overlap the de novo variant's REF span are dropped, so the planted allele
// is never shadowed by a simulated one.
func GenerateTrioHaplotypes(
	result *population.HaplotypeSimulationResult,
	deNovo *genomics.VariantRecord,
	rng *rand.Rand,
) (*TrioHaplotypes, error) {
	if result == nil {
		return nil, fmt.Errorf("nil simulation result")
	}
	if len(result.HaplotypeMatrix) < 6 {
		return nil, fmt.Errorf("haplotype matrix has %d rows; need at least 6 haplotypes",
			len(result.HaplotypeMatrix))
	}
	nVariants := len(result.Variants)
	// Validate matrix columns match variant count.
	for i, row := range result.HaplotypeMatrix[:6] {
		if len(row) != nVariants {
			return nil, fmt.Errorf("row %d has %d entries; expected %d variants", i, len(row), nVariants)
		}
	}

	keep := make([]int, 0, nVariants)
	for i, v := range result.Variants {
		if !overlapsDeNovo(v, deNovo) {
			keep = append(keep, i)
		}
	}
	variants := make([]population.HaplotypeVariant, len(keep))
	positions := make([]int, len(keep))
	for k, i := range keep {
		variants[k] = result.Variants[i]
		positions[k] = result.Variants[i].Position
	}
	row := func(r int) []int {
		h := make([]int, len(keep))
		for k, i := range keep {
			h[k] = result.HaplotypeMatrix[r][i]
		}
		return h
	}
	motherH1, motherH2 := row(0), row(1)
	fatherH1, fatherH2 := row(2), row(3)

	// Genetic length of the region; 0 (or unset) means no crossovers.
	cm := max(result.Statistics.GeneticMapLengthCM, 0)

	// Simulate recombination for proband from each parent.
	probandH1 := simulateCrossover(motherH1, motherH2, positions, cm, rng)
	probandH2 := simulateCrossover(fatherH1, fatherH2, positions, cm, rng)

	return &TrioHaplotypes{
		MotherH1: motherH1, MotherH2: motherH2,
		FatherH1: fatherH1, FatherH2: fatherH2,
		ProbandH1: probandH1, ProbandH2: probandH2,
		Variants:      variants,
		DeNovoVariant: deNovo,
	}, nil
}

// overlapsDeNovo reports whether background variant v (0-based Position)
// overlaps the REF span of the 1-based de novo variant.
func overlapsDeNovo(v population.HaplotypeVariant, dn *genomics.VariantRecord) bool {
	if dn == nil || dn.Position <= 0 {
		return false
	}
	vLen := 1
	if len(v.Alleles) > 0 && len(v.Alleles[0]) > 0 {
		vLen = len(v.Alleles[0])
	}
	dLen := int64(max(len(dn.Ref), 1))
	vs, ds := int64(v.Position), dn.Position-1
	return vs < ds+dLen && ds < vs+int64(vLen)
}

// simulateCrossover returns a recombinant haplotype from h1 and h2 using a
// Poisson crossover model: the number of crossovers is Poisson with mean
// geneticMapCM / 100 (Morgans), their positions are uniform over the span of
// the variant positions, and the starting strand is h1 or h2 with
// probability 1/2, so without a crossover either parental haplotype is
// transmitted whole with equal probability.
func simulateCrossover(h1, h2 []int, positions []int, geneticMapCM float64, rng *rand.Rand) []int {
	result := make([]int, len(h1))
	if len(h1) == 0 {
		return result
	}

	strand := rng.Intn(2) // 0 = h1, 1 = h2
	numCrossovers := poissonSample(geneticMapCM/100.0, rng)

	lo, hi := positions[0], positions[len(positions)-1]+1
	crossovers := make([]int, numCrossovers)
	for i := range crossovers {
		crossovers[i] = lo + rng.Intn(hi-lo)
	}
	sortInts(crossovers)

	// Switch strands at each crossover.
	xIdx := 0
	for i, pos := range positions {
		for xIdx < len(crossovers) && crossovers[xIdx] <= pos {
			strand = 1 - strand
			xIdx++
		}
		if strand == 0 {
			result[i] = h1[i]
		} else {
			result[i] = h2[i]
		}
	}
	return result
}

// poissonSample draws a sample from Poisson(lambda) using the Knuth algorithm.
func poissonSample(lambda float64, rng *rand.Rand) int {
	if lambda <= 0 {
		return 0
	}
	L := math.Exp(-lambda)
	k := 0
	p := 1.0
	for p > L {
		k++
		p *= rng.Float64()
	}
	return k - 1
}

// sortInts sorts a small slice of ints using insertion sort (avoids sort import for tiny slices).
func sortInts(a []int) {
	for i := 1; i < len(a); i++ {
		key := a[i]
		j := i - 1
		for j >= 0 && a[j] > key {
			a[j+1] = a[j]
			j--
		}
		a[j+1] = key
	}
}
