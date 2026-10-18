package statistics

import (
	"fmt"
	"math"
)

// HWEResult holds the result of a Hardy-Weinberg Equilibrium chi-square test.
type HWEResult struct {
	// ChiSquared is the χ² test statistic (df = 1).
	ChiSquared float64
	// PValue is the two-tailed p-value derived from the chi-square distribution (df=1).
	// Computed as erfc(sqrt(χ²/2)), which is exact for df=1.
	PValue float64
	// ObservedFreq is the MLE allele frequency p̂ = (2·n_AA + n_Aa) / (2·N).
	ObservedFreq float64
	// PassesHWE is true when PValue > 0.05 (cannot reject HWE at α=0.05).
	PassesHWE bool
}

// TestHWE performs a chi-square test for Hardy-Weinberg Equilibrium.
//
// Given observed genotype counts n_HomRef (AA), n_Het (Aa), n_HomAlt (aa):
//
//	N  = n_AA + n_Aa + n_aa
//	p̂  = (2·n_AA + n_Aa) / (2·N)       MLE allele frequency
//	E_AA = N·p̂²
//	E_Aa = 2·N·p̂·(1−p̂)
//	E_aa = N·(1−p̂)²
//
//	χ² = (n_AA−E_AA)²/E_AA + (n_Aa−E_Aa)²/E_Aa + (n_aa−E_aa)²/E_aa  (df=1)
//	p  = erfc(sqrt(χ²/2))
//
// Returns an error for monomorphic sites (p̂ = 0 or 1) to guard against
// division by zero in the expected-count calculation.
//
// This is the asymptotic chi-square test, not the exact test described by
// Wigginton, Cutler & Abecasis (2005), A note on exact tests of
// Hardy-Weinberg equilibrium, Am. J. Hum. Genet. 76, 887-893, which they
// recommend over chi-square for small samples and rare alleles.
func TestHWE(nHomRef, nHet, nHomAlt int) (HWEResult, error) {
	if nHomRef < 0 || nHet < 0 || nHomAlt < 0 {
		return HWEResult{}, fmt.Errorf("genotype counts must be non-negative")
	}

	N := nHomRef + nHet + nHomAlt
	if N == 0 {
		return HWEResult{}, fmt.Errorf("total count N is zero")
	}

	// MLE allele frequency.
	pHat := (2.0*float64(nHomRef) + float64(nHet)) / (2.0 * float64(N))

	if pHat == 0 || pHat == 1 {
		return HWEResult{}, fmt.Errorf("monomorphic site: allele frequency is %.1f (division by zero in expected counts)", pHat)
	}

	qHat := 1.0 - pHat
	Nf := float64(N)

	// Expected genotype counts under HWE.
	eAA := Nf * pHat * pHat
	eAa := 2.0 * Nf * pHat * qHat
	eaa := Nf * qHat * qHat

	// Chi-square statistic (df = 1 because one allele frequency is estimated).
	chiSq := 0.0
	if eAA > 0 {
		diff := float64(nHomRef) - eAA
		chiSq += diff * diff / eAA
	}
	if eAa > 0 {
		diff := float64(nHet) - eAa
		chiSq += diff * diff / eAa
	}
	if eaa > 0 {
		diff := float64(nHomAlt) - eaa
		chiSq += diff * diff / eaa
	}

	// p-value for df=1: p = erfc(sqrt(χ²/2)).
	pValue := math.Erfc(math.Sqrt(chiSq / 2.0))

	return HWEResult{
		ChiSquared:   chiSq,
		PValue:       pValue,
		ObservedFreq: pHat,
		PassesHWE:    pValue > 0.05,
	}, nil
}
