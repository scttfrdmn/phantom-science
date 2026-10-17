//go:build !release

package statistics

import (
	"math"
	"testing"
)

func TestHWE_PerfectEquilibrium(t *testing.T) {
	// Counts derived directly from HWE with p=0.6, N=1000:
	//   E_AA = 1000*0.36 = 360, E_Aa = 1000*0.48 = 480, E_aa = 1000*0.16 = 160
	res, err := TestHWE(360, 480, 160)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// χ² should be very close to 0 → p-value close to 1.
	if res.ChiSquared > 1e-6 {
		t.Errorf("expected χ²≈0 for perfect HWE input, got %.8f", res.ChiSquared)
	}
	if res.PValue < 0.95 {
		t.Errorf("expected p-value ≥ 0.95 for perfect HWE, got %.4f", res.PValue)
	}
	if !res.PassesHWE {
		t.Error("expected PassesHWE=true for perfect HWE input")
	}
}

func TestHWE_ExcessHomozygosity(t *testing.T) {
	// Strongly excess homozygosity (inbreeding): far from HWE.
	// p≈0.5: expected ~250 AA, 500 Aa, 250 aa, but we supply 450, 100, 450.
	res, err := TestHWE(450, 100, 450)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if res.PValue >= 0.01 {
		t.Errorf("expected p-value < 0.01 for extreme inbreeding, got %.4f", res.PValue)
	}
	if res.PassesHWE {
		t.Error("expected PassesHWE=false for extreme inbreeding")
	}
}

func TestHWE_MonomorphicSite(t *testing.T) {
	// All individuals are homozygous ref → allele freq = 1.0 → monomorphic.
	_, err := TestHWE(100, 0, 0)
	if err == nil {
		t.Error("expected error for monomorphic site (p=1), got nil")
	}

	_, err = TestHWE(0, 0, 100)
	if err == nil {
		t.Error("expected error for monomorphic site (p=0), got nil")
	}
}

func TestHWE_SmallN(t *testing.T) {
	// N = 10: manually compute expected values and verify.
	// Observed: AA=4, Aa=4, aa=2, N=10
	// p̂ = (8+4)/20 = 0.6, q̂ = 0.4
	// E_AA=10*0.36=3.6, E_Aa=10*0.48=4.8, E_aa=10*0.16=1.6
	// χ² = (4-3.6)²/3.6 + (4-4.8)²/4.8 + (2-1.6)²/1.6
	//     = 0.16/3.6 + 0.64/4.8 + 0.16/1.6
	//     ≈ 0.04444 + 0.13333 + 0.10000 = 0.27778
	res, err := TestHWE(4, 4, 2)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expectedChiSq := 0.04444 + 0.13333 + 0.10000
	if math.Abs(res.ChiSquared-expectedChiSq) > 0.001 {
		t.Errorf("χ² mismatch: expected ≈%.5f, got %.5f", expectedChiSq, res.ChiSquared)
	}
	if math.Abs(res.ObservedFreq-0.6) > 1e-9 {
		t.Errorf("allele freq mismatch: expected 0.6, got %.9f", res.ObservedFreq)
	}
}

func TestHWE_ZeroTotal(t *testing.T) {
	_, err := TestHWE(0, 0, 0)
	if err == nil {
		t.Error("expected error for N=0, got nil")
	}
}

func TestHWE_NegativeCounts(t *testing.T) {
	_, err := TestHWE(-1, 5, 3)
	if err == nil {
		t.Error("expected error for negative count, got nil")
	}
}

func TestHWE_PValueRange(t *testing.T) {
	// p-value must always be in [0, 1].
	cases := [][3]int{
		{500, 0, 500},   // extreme excess homozygosity
		{250, 500, 250}, // perfect HWE at p=0.5
		{1, 1, 998},     // near-monomorphic alt
	}
	for _, c := range cases {
		res, err := TestHWE(c[0], c[1], c[2])
		if err != nil {
			// Monomorphic-adjacent cases may return an error; skip.
			continue
		}
		if res.PValue < 0 || res.PValue > 1 {
			t.Errorf("p-value out of [0,1] for (%d,%d,%d): %.6f", c[0], c[1], c[2], res.PValue)
		}
	}
}
