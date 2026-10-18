//go:build !release

package population_test

import (
	"math"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/statistics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/population"
)

// meanSE returns the mean and standard error of xs.
func meanSE(xs []float64) (float64, float64) {
	n := float64(len(xs))
	m := 0.0
	for _, x := range xs {
		m += x
	}
	m /= n
	ss := 0.0
	for _, x := range xs {
		ss += (x - m) * (x - m)
	}
	return m, math.Sqrt(ss / (n - 1) / n)
}

// jackknife computes a statistic over loci and its delete-one-block jackknife
// standard error. stat receives a predicate selecting the loci to use.
func jackknife(nLoci, nBlocks int, stat func(keep func(locus int) bool) float64) (float64, float64) {
	full := stat(func(int) bool { return true })
	vals := make([]float64, nBlocks)
	mean := 0.0
	for b := 0; b < nBlocks; b++ {
		bb := b
		vals[b] = stat(func(l int) bool { return l*nBlocks/nLoci != bb })
		mean += vals[b]
	}
	mean /= float64(nBlocks)
	ss := 0.0
	for _, v := range vals {
		ss += (v - mean) * (v - mean)
	}
	g := float64(nBlocks)
	return full, math.Sqrt((g - 1) / g * ss)
}

// wcFST is the Weir–Cockerham F_ST between two sampled populations of a joint
// panel, combined over sites as a ratio of averages (Σ numerator / Σ
// denominator), using statistics.ComputeFST for the per-site components.
func wcFST(jp *population.JointPanel, a, b string, keep func(int) bool) float64 {
	ia, ib := jp.PopIndex(a), jp.PopIndex(b)
	na, nb := jp.Haplotypes[ia], jp.Haplotypes[ib]
	num, den := 0.0, 0.0
	for _, s := range jp.Sites {
		if !keep(s.Locus) {
			continue
		}
		ca, cb := s.DerivedCount[ia], s.DerivedCount[ib]
		if (ca == 0 && cb == 0) || (ca == na && cb == nb) {
			continue
		}
		r := statistics.ComputeFST([]statistics.PopulationSample{
			{AlleleFreq: float64(ca) / float64(na), SampleSize: na},
			{AlleleFreq: float64(cb) / float64(nb), SampleSize: nb},
		})
		num += r.MSG - r.MSW
		den += r.MSG + (r.Nc-1)*r.MSW
	}
	if den == 0 {
		return 0
	}
	return num / den
}

// hudsonFST is Hudson's F_ST (Hudson et al. 1992, as in Bhatia et al. 2013),
// ratio of averages over sites. Its expectation is 1 − T_within/T_between.
func hudsonFST(jp *population.JointPanel, a, b string, keep func(int) bool) float64 {
	ia, ib := jp.PopIndex(a), jp.PopIndex(b)
	n1, n2 := float64(jp.Haplotypes[ia]), float64(jp.Haplotypes[ib])
	num, den := 0.0, 0.0
	for _, s := range jp.Sites {
		if !keep(s.Locus) {
			continue
		}
		p1 := float64(s.DerivedCount[ia]) / n1
		p2 := float64(s.DerivedCount[ib]) / n2
		num += (p1-p2)*(p1-p2) - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
		den += p1*(1-p2) + p2*(1-p1)
	}
	return num / den
}

// perLocusPi returns π per locus (differences per locus) for one population.
func perLocusPi(jp *population.JointPanel, pop string) []float64 {
	s, err := jp.Summary(pop)
	if err != nil {
		panic(err)
	}
	return s.PerLocusPi
}

// perLocusRare counts, per locus, sites with derived count in (0, 5% of n].
func perLocusRare(jp *population.JointPanel, pop string) []float64 {
	i := jp.PopIndex(pop)
	n := jp.Haplotypes[i]
	out := make([]float64, jp.NumLoci)
	for _, s := range jp.Sites {
		c := s.DerivedCount[i]
		if c > 0 && float64(c) <= 0.05*float64(n) {
			out[s.Locus]++
		}
	}
	return out
}

func diff(a, b []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] - b[i]
	}
	return out
}
