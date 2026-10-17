//go:build !release

package statistics

import (
	"fmt"
	"math"
)

// LDResult holds the pairwise linkage disequilibrium statistics computed from
// phased haplotype counts.
type LDResult struct {
	// D is the linkage disequilibrium coefficient: D = p_AB − p_A·p_B.
	D float64
	// DPrime is the normalised LD coefficient D' ∈ [−1, 1].
	DPrime float64
	// RSquared is the squared correlation coefficient r² ∈ [0, 1].
	RSquared float64
	// FreqA is the frequency of allele A at the first locus.
	FreqA float64
	// FreqB is the frequency of allele B at the second locus.
	FreqB float64
	// Valid is false when LD cannot be computed (e.g. monomorphic locus).
	Valid bool
}

// ComputeLD calculates pairwise linkage disequilibrium statistics from phased
// haplotype counts.
//
// Haplotype notation (A/a = alleles at locus 1, B/b = alleles at locus 2):
//
//	nAB  = count of haplotypes carrying A and B
//	nAb  = count of haplotypes carrying A and b
//	naB  = count of haplotypes carrying a and B
//	nab  = count of haplotypes carrying a and b
//
// Derived quantities (N = total haplotypes):
//
//	p_A = (n_AB + n_Ab) / N,   q_a = 1 − p_A
//	p_B = (n_AB + n_aB) / N,   q_b = 1 − p_B
//	p_AB = n_AB / N
//
//	D    = p_AB − p_A·p_B
//	r²   = D² / (p_A·q_a·p_B·q_b)
//
//	D_max = min(p_A·q_b, q_a·p_B)   if D > 0
//	       max(−p_A·p_B, −q_a·q_b)  if D ≤ 0
//	D'   = D / D_max
//
// Returns an error for monomorphic loci (p_A = 0 or 1, or p_B = 0 or 1).
//
// References:
//   - Lewontin, R.C. (1964). The interaction of selection and linkage.
//     Genetics, 49(1), 49-67.
//   - Hill, W.G. & Robertson, A. (1968). Linkage disequilibrium in finite
//     populations. Theoretical and Applied Genetics, 38(6), 226-231.
func ComputeLD(nAB, nAb, naB, nab int) (LDResult, error) {
	if nAB < 0 || nAb < 0 || naB < 0 || nab < 0 {
		return LDResult{}, fmt.Errorf("haplotype counts must be non-negative")
	}

	N := nAB + nAb + naB + nab
	if N == 0 {
		return LDResult{}, fmt.Errorf("total haplotype count N is zero")
	}

	Nf := float64(N)
	pA := (float64(nAB) + float64(nAb)) / Nf
	qA := 1.0 - pA
	pB := (float64(nAB) + float64(naB)) / Nf
	qB := 1.0 - pB

	// Guard against monomorphic loci.
	if pA == 0 || pA == 1 {
		return LDResult{FreqA: pA, FreqB: pB, Valid: false}, nil
	}
	if pB == 0 || pB == 1 {
		return LDResult{FreqA: pA, FreqB: pB, Valid: false}, nil
	}

	pAB := float64(nAB) / Nf
	D := pAB - pA*pB

	// r²: squared Pearson correlation between allele indicator variables.
	denom := pA * qA * pB * qB
	var rSq float64
	if denom > 0 {
		rSq = (D * D) / denom
		// Clamp to [0, 1] due to floating-point rounding.
		if rSq > 1.0 {
			rSq = 1.0
		}
	}

	// D': normalised LD coefficient.
	var dPrime float64
	if D > 0 {
		dMax := math.Min(pA*qB, qA*pB)
		if dMax > 0 {
			dPrime = D / dMax
		}
	} else if D < 0 {
		dMax := math.Max(-pA*pB, -qA*qB)
		if dMax < 0 {
			dPrime = D / dMax
		}
	}
	// D = 0 → D' = 0 (no LD).

	// Clamp D' to [−1, 1] for numerical stability.
	if dPrime > 1.0 {
		dPrime = 1.0
	}
	if dPrime < -1.0 {
		dPrime = -1.0
	}

	return LDResult{
		D:        D,
		DPrime:   dPrime,
		RSquared: rSq,
		FreqA:    pA,
		FreqB:    pB,
		Valid:    true,
	}, nil
}
