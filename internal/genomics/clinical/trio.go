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

// GenerateTrioHaplotypes takes a msprime HaplotypeSimulationResult (≥6 haplotypes
// required — request sampleSize≥3 with includeMatrix=true), designates the first 6
// rows as mother_h1, mother_h2, father_h1, father_h2, bg_h1, bg_h2, then simulates
// Mendelian recombination to produce proband_h1 (recombinant of mother's pair) and
// proband_h2 (recombinant of father's pair), and records deNovo on ProbandH1.
func GenerateTrioHaplotypes(
    result *population.HaplotypeSimulationResult,
    deNovo *genomics.VariantRecord,
    rng *rand.Rand,
) (*TrioHaplotypes, error) {
    if result == nil {
        return nil, fmt.Errorf("nil simulation result")
    }
    if len(result.HaplotypeMatrix) < 6 {
        return nil, fmt.Errorf("haplotype matrix has %d rows; need ≥6 (sampleSize≥3 with includeMatrix=true)",
            len(result.HaplotypeMatrix))
    }
    nVariants := len(result.Variants)
    // Validate matrix columns match variant count.
    for i, row := range result.HaplotypeMatrix[:6] {
        if len(row) != nVariants {
            return nil, fmt.Errorf("row %d has %d entries; expected %d variants", i, len(row), nVariants)
        }
    }

    motherH1 := result.HaplotypeMatrix[0]
    motherH2 := result.HaplotypeMatrix[1]
    fatherH1 := result.HaplotypeMatrix[2]
    fatherH2 := result.HaplotypeMatrix[3]
    // rows 4,5 are background haplotypes (unused for trio members)

    // Extract variant positions for crossover simulation.
    positions := make([]int, nVariants)
    for i, v := range result.Variants {
        positions[i] = v.Position
    }

    cm := result.Statistics.GeneticMapLengthCM
    if cm <= 0 {
        cm = 100.0 // default 100 cM (approximately one chromosome arm)
    }

    // Simulate recombination for proband from each parent.
    probandH1 := simulateCrossover(motherH1, motherH2, positions, cm, rng)
    probandH2 := simulateCrossover(fatherH1, fatherH2, positions, cm, rng)

    return &TrioHaplotypes{
        MotherH1: motherH1, MotherH2: motherH2,
        FatherH1: fatherH1, FatherH2: fatherH2,
        ProbandH1: probandH1, ProbandH2: probandH2,
        Variants:      result.Variants,
        DeNovoVariant: deNovo,
    }, nil
}

// simulateCrossover returns a recombinant haplotype from h1 and h2 using a
// Poisson crossover model. Expected crossovers = geneticMapLengthCM / 100.
// Crossover positions are drawn uniformly from the physical positions of variants.
func simulateCrossover(h1, h2 []int, positions []int, geneticMapCM float64, rng *rand.Rand) []int {
    result := make([]int, len(h1))
    if len(h1) == 0 {
        return result
    }

    lambda := geneticMapCM / 100.0
    numCrossovers := poissonSample(lambda, rng)

    if numCrossovers == 0 {
        // No crossover: inherit entire h1 (arbitrary choice of starting strand).
        copy(result, h1)
        return result
    }

    // Draw crossover positions uniformly from [0, maxPos).
    maxPos := positions[len(positions)-1] + 1
    if maxPos <= 0 {
        maxPos = 1000000
    }
    crossovers := make([]int, numCrossovers)
    for i := range crossovers {
        crossovers[i] = rng.Intn(maxPos)
    }
    // Sort crossovers.
    sortInts(crossovers)

    // Build recombinant: start on h1, switch strands at each crossover.
    strand := 0 // 0 = h1, 1 = h2
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
