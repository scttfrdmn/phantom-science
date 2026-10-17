//go:build !release

package statistics

import (
	"math"
	"testing"
)

func TestComputeLD_CompleteLD(t *testing.T) {
	// No recombinant haplotypes: n_Ab = n_aB = 0
	// → D = p_AB - p_A*p_B = p_A - p_A*(p_A) = p_A*(1-p_A) = p_A*q_A
	// → r² = D²/(p_A*q_A*p_B*q_B) = (p_A*q_A)²/(p_A*q_A*p_A*q_A) = 1
	// → D' = 1
	res, err := ComputeLD(30, 0, 0, 70)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !res.Valid {
		t.Fatal("expected Valid=true for complete LD")
	}
	if math.Abs(res.RSquared-1.0) > 1e-6 {
		t.Errorf("expected r²=1.0 for complete LD, got %.6f", res.RSquared)
	}
	if math.Abs(res.DPrime-1.0) > 1e-6 {
		t.Errorf("expected D'=1.0 for complete LD, got %.6f", res.DPrime)
	}
}

func TestComputeLD_Independence(t *testing.T) {
	// Alleles at the two loci are in linkage equilibrium.
	// p_A = 0.5, p_B = 0.5 → D = 0 when n_AB = N*p_A*p_B
	// N=100: n_AB=25, n_Ab=25, n_aB=25, n_ab=25
	res, err := ComputeLD(25, 25, 25, 25)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !res.Valid {
		t.Fatal("expected Valid=true for independence")
	}
	if math.Abs(res.D) > 1e-9 {
		t.Errorf("expected D=0 for independence, got %.9f", res.D)
	}
	if math.Abs(res.RSquared) > 1e-9 {
		t.Errorf("expected r²=0 for independence, got %.9f", res.RSquared)
	}
}

func TestComputeLD_MonomorphicFirstLocus(t *testing.T) {
	// All haplotypes carry the same allele at locus 1 → monomorphic.
	res, err := ComputeLD(50, 50, 0, 0)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if res.Valid {
		t.Error("expected Valid=false for monomorphic locus 1")
	}
}

func TestComputeLD_MonomorphicSecondLocus(t *testing.T) {
	res, err := ComputeLD(50, 0, 50, 0)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if res.Valid {
		t.Error("expected Valid=false for monomorphic locus 2")
	}
}

func TestComputeLD_NegativeCounts(t *testing.T) {
	_, err := ComputeLD(-1, 25, 25, 25)
	if err == nil {
		t.Error("expected error for negative count")
	}
}

func TestComputeLD_ZeroN(t *testing.T) {
	_, err := ComputeLD(0, 0, 0, 0)
	if err == nil {
		t.Error("expected error for N=0")
	}
}

func TestComputeLD_AsymmetricLD(t *testing.T) {
	// Partial LD: D' and r² should both be in (0, 1).
	// n_AB=40, n_Ab=10, n_aB=10, n_ab=40; N=100
	// p_A = 0.5, p_B = 0.5
	// p_AB = 0.4, D = 0.4 - 0.25 = 0.15
	// D_max = min(p_A*q_B, q_A*p_B) = min(0.25,0.25) = 0.25
	// D' = 0.15/0.25 = 0.6
	// r² = 0.15²/(0.25*0.25) = 0.0225/0.0625 = 0.36
	res, err := ComputeLD(40, 10, 10, 40)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !res.Valid {
		t.Fatal("expected Valid=true")
	}
	if math.Abs(res.D-0.15) > 1e-9 {
		t.Errorf("D mismatch: expected 0.15, got %.9f", res.D)
	}
	if math.Abs(res.DPrime-0.6) > 1e-9 {
		t.Errorf("D' mismatch: expected 0.6, got %.9f", res.DPrime)
	}
	if math.Abs(res.RSquared-0.36) > 1e-9 {
		t.Errorf("r² mismatch: expected 0.36, got %.9f", res.RSquared)
	}
}

func TestComputeLD_RSquaredBounded(t *testing.T) {
	// r² must always be in [0, 1].
	cases := [][4]int{
		{10, 40, 40, 10},
		{1, 1, 1, 1},
		{100, 1, 1, 100},
	}
	for _, c := range cases {
		res, err := ComputeLD(c[0], c[1], c[2], c[3])
		if err != nil {
			continue
		}
		if !res.Valid {
			continue
		}
		if res.RSquared < 0 || res.RSquared > 1.0+1e-9 {
			t.Errorf("r² out of [0,1] for (%d,%d,%d,%d): %.6f", c[0], c[1], c[2], c[3], res.RSquared)
		}
		if res.DPrime < -1.0-1e-9 || res.DPrime > 1.0+1e-9 {
			t.Errorf("D' out of [-1,1] for (%d,%d,%d,%d): %.6f", c[0], c[1], c[2], c[3], res.DPrime)
		}
	}
}
