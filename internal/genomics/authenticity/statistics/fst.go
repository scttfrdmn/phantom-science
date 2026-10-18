//go:build !release

package statistics

import "math"

// PopulationSample holds allele frequency data for a single population at a single locus.
type PopulationSample struct {
	// AlleleFreq is the frequency of the derived (or minor) allele in [0, 1].
	AlleleFreq float64
	// SampleSize is n_j in ComputeFST: the number of sampled allele copies
	// (chromosomes), not individuals.
	SampleSize int
}

// FSTResult holds the computed Weir-Cockerham (1984) F_ST estimator and its
// intermediate quantities.
type FSTResult struct {
	// FST is the estimated fixation index, clamped to [0, 1].
	FST float64
	// MSG is the mean square deviation between populations.
	MSG float64
	// MSW is the mean square deviation within populations.
	MSW float64
	// Nc is the effective sample-size correction factor.
	Nc float64
	// Valid is false when F_ST cannot be computed (e.g. MSG ≤ 0, < 2 populations,
	// or all sample sizes are zero).
	Valid bool
}

// ComputeFST returns the Weir-Cockerham (1984) F_ST (theta) estimator for a
// single biallelic locus, in the form without the observed-heterozygosity
// term (i.e. treating the sample as allele copies). Negative estimates are
// clamped to 0.
//
// Mathematical derivation (r populations, allele frequency p_j in population j,
// sample size n_j, N = Σ n_j):
//
//	p̄  = Σ n_j·p_j / N                      weighted mean allele frequency
//	n_c = (N − Σ n_j²/N) / (r−1)             size-correction factor
//
//	MSG = Σ n_j·(p_j − p̄)² / (r−1)          mean sq. between populations
//	MSW = Σ n_j·p_j·(1−p_j) / (N−r)          mean sq. within populations
//
//	F_ST = (MSG − MSW) / (MSG + (n_c−1)·MSW)
//
// Returns an invalid FSTResult when MSG ≤ 0 (no between-population variance) or
// when the input is degenerate.
//
// Reference: Weir, B.S. & Cockerham, C.C. (1984). Estimating F-statistics for
// the analysis of population structure. Evolution, 38(6), 1358-1370.
func ComputeFST(populations []PopulationSample) FSTResult {
	r := len(populations)
	if r < 2 {
		return FSTResult{}
	}

	// Compute N and weighted mean allele frequency p̄.
	N := 0
	sumNjSq := 0.0
	for _, pop := range populations {
		N += pop.SampleSize
	}
	if N == 0 {
		return FSTResult{}
	}
	Nf := float64(N)

	pBar := 0.0
	for _, pop := range populations {
		pBar += float64(pop.SampleSize) * pop.AlleleFreq
	}
	pBar /= Nf

	// Size-correction factor n_c.
	for _, pop := range populations {
		nj := float64(pop.SampleSize)
		sumNjSq += nj * nj
	}
	nc := (Nf - sumNjSq/Nf) / float64(r-1)

	// MSG: mean square between populations.
	msgSum := 0.0
	for _, pop := range populations {
		diff := pop.AlleleFreq - pBar
		msgSum += float64(pop.SampleSize) * diff * diff
	}
	msg := msgSum / float64(r-1)

	// MSW: mean square within populations.
	msw := 0.0
	if N-r > 0 {
		mswSum := 0.0
		for _, pop := range populations {
			mswSum += float64(pop.SampleSize) * pop.AlleleFreq * (1.0 - pop.AlleleFreq)
		}
		msw = mswSum / float64(N-r)
	}

	// F_ST cannot be estimated when there is no between-population variance.
	if msg <= 0 {
		return FSTResult{FST: 0, MSG: msg, MSW: msw, Nc: nc, Valid: false}
	}

	denom := msg + (nc-1)*msw
	if denom == 0 {
		return FSTResult{FST: 0, MSG: msg, MSW: msw, Nc: nc, Valid: false}
	}

	fst := (msg - msw) / denom

	// Clamp to [0, 1].
	if math.IsNaN(fst) || fst < 0 {
		fst = 0
	}
	if fst > 1 {
		fst = 1
	}

	return FSTResult{FST: fst, MSG: msg, MSW: msw, Nc: nc, Valid: true}
}

// ComputeMultiLocusFST averages Weir-Cockerham F_ST across multiple loci, weighting
// equally. Only valid loci (FSTResult.Valid == true) are included in the average.
//
// Returns an invalid FSTResult when no locus yields a valid estimate.
func ComputeMultiLocusFST(loci [][]PopulationSample) FSTResult {
	if len(loci) == 0 {
		return FSTResult{}
	}

	sumFST := 0.0
	n := 0
	for _, locus := range loci {
		res := ComputeFST(locus)
		if res.Valid {
			sumFST += res.FST
			n++
		}
	}

	if n == 0 {
		return FSTResult{}
	}

	return FSTResult{
		FST:   sumFST / float64(n),
		Valid: true,
	}
}
