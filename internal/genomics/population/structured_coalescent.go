package population

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sort"
)

// DefaultLocusLength is the length (bp) of one simulated locus. Loci are
// treated as unlinked and internally non-recombining: each gets its own
// genealogy.
const DefaultLocusLength = 10000

// PopSpec describes one population (one branch of the population tree) in a
// multi-population demographic model. All times are in generations before the
// present.
type PopSpec struct {
	Name string

	// Ne is the effective size of the population before its recent growth
	// epoch (the long-term branch size).
	Ne float64

	// GrowthRate is the per-generation forward-time exponential growth rate of
	// the recent epoch [0, GrowthStart). The present-day size is
	// Ne·exp(GrowthRate·GrowthStart).
	GrowthRate float64
	// GrowthStart is how many generations ago growth began (0 = no growth).
	GrowthStart int

	// Bottleneck: the size is BottleneckNe during
	// [BottleneckStart, BottleneckStart+BottleneckDuration). Ignored when
	// BottleneckDuration <= 0 or BottleneckNe <= 0. A bottleneck that starts
	// at or after SplitTime falls in the ancestral population's branch and
	// has no effect on this population.
	BottleneckStart    int
	BottleneckDuration int
	BottleneckNe       float64

	// SplitTime is when (generations ago) this population split from
	// SplitFrom. Going backwards in time, at SplitTime every remaining lineage
	// of this population moves into SplitFrom. The root population has
	// SplitFrom == "" and SplitTime is ignored.
	SplitTime int
	SplitFrom string

	// Events are additional DemographicEvents applied to this branch (see
	// size_history.go for their semantics). They only act on [0, SplitTime).
	Events []DemographicEvent
}

// presentSize returns the present-day size implied by Ne and recent growth.
func (p PopSpec) presentSize() float64 {
	if p.GrowthRate != 0 && p.GrowthStart > 0 {
		return p.Ne * math.Exp(p.GrowthRate*float64(p.GrowthStart))
	}
	return p.Ne
}

// epochs returns the size history of this branch.
func (p PopSpec) epochs() []epoch {
	events := append([]DemographicEvent(nil), p.Events...)
	if p.GrowthRate != 0 && p.GrowthStart > 0 {
		events = append(events, DemographicEvent{
			Type:       "expansion",
			Time:       0,
			Duration:   p.GrowthStart,
			Parameters: map[string]float64{"growth_rate": p.GrowthRate},
		})
	}
	if p.BottleneckDuration > 0 && p.BottleneckNe > 0 {
		events = append(events, DemographicEvent{
			Type:       "bottleneck",
			Time:       p.BottleneckStart,
			Duration:   p.BottleneckDuration,
			Parameters: map[string]float64{"population_size": p.BottleneckNe},
		})
	}
	return buildEpochs(p.presentSize(), events)
}

// Demography is a multi-population demographic model: a tree of populations
// connected by splits (no migration, no admixture).
type Demography struct {
	Pops        []PopSpec
	Mu          float64 // per-bp per-generation mutation rate
	LocusLength int     // bp per locus (DefaultLocusLength when 0)
}

// Pop returns the spec of the named population.
func (d Demography) Pop(name string) (PopSpec, bool) {
	for _, p := range d.Pops {
		if p.Name == name {
			return p, true
		}
	}
	return PopSpec{}, false
}

// Validate checks that the demography is a well-formed population tree.
func (d Demography) Validate() error {
	if len(d.Pops) == 0 {
		return fmt.Errorf("demography has no populations")
	}
	if d.Mu < 0 || math.IsNaN(d.Mu) {
		return fmt.Errorf("mutation rate must be non-negative")
	}
	if d.LocusLength < 0 {
		return fmt.Errorf("locus length must be non-negative")
	}
	idx := make(map[string]int, len(d.Pops))
	for i, p := range d.Pops {
		if p.Name == "" {
			return fmt.Errorf("population %d has no name", i)
		}
		if _, dup := idx[p.Name]; dup {
			return fmt.Errorf("duplicate population %q", p.Name)
		}
		if !(p.Ne >= 1) {
			return fmt.Errorf("population %q: Ne must be >= 1", p.Name)
		}
		if p.GrowthStart < 0 || p.BottleneckStart < 0 || p.BottleneckDuration < 0 || p.SplitTime < 0 {
			return fmt.Errorf("population %q: times must be non-negative", p.Name)
		}
		idx[p.Name] = i
	}
	roots := 0
	for _, p := range d.Pops {
		if p.SplitFrom == "" {
			roots++
			continue
		}
		j, ok := idx[p.SplitFrom]
		if !ok {
			return fmt.Errorf("population %q splits from unknown population %q", p.Name, p.SplitFrom)
		}
		parent := d.Pops[j]
		if parent.SplitFrom != "" && parent.SplitTime < p.SplitTime {
			return fmt.Errorf("population %q splits from %q at %d, after %q itself split (at %d)",
				p.Name, p.SplitFrom, p.SplitTime, p.SplitFrom, parent.SplitTime)
		}
	}
	if roots != 1 {
		return fmt.Errorf("demography must have exactly one root population, found %d", roots)
	}
	// Cycle check: follow parents.
	for _, p := range d.Pops {
		seen := map[string]bool{}
		for cur := p; cur.SplitFrom != ""; cur = d.Pops[idx[cur.SplitFrom]] {
			if seen[cur.Name] {
				return fmt.Errorf("population tree has a cycle through %q", cur.Name)
			}
			seen[cur.Name] = true
		}
	}
	return nil
}

// JointSite is one segregating site (one mutation, infinite-sites model) in a
// joint multi-population sample.
type JointSite struct {
	Locus     int
	Position  int // global position: Locus*LocusLength + offset within locus
	Ancestral byte
	Derived   byte
	Age       float64 // generations ago the mutation arose
	// DerivedCount[i] is the number of sampled haplotypes of
	// JointPanel.Populations[i] carrying the derived allele.
	DerivedCount []int
	// Genotypes[i] holds 0/1/2 derived-allele dosages for the diploid
	// individuals of Populations[i] (haplotypes 2k and 2k+1 form individual
	// k). Nil for a population with an odd haplotype count.
	Genotypes [][]int
	// Carriers[i][h] is 1 when haplotype h of Populations[i] carries the
	// derived allele and 0 otherwise: the phased haplotype matrix of the
	// sample at this site (one row per population, one column per
	// haplotype).
	Carriers [][]byte
}

// JointPanel is a sample from several populations simulated on shared
// genealogies, one per unlinked locus, so every population is scored at the
// same set of sites.
type JointPanel struct {
	Populations []string // sampled populations (sorted)
	Haplotypes  []int    // haplotype sample size per population
	NumLoci     int
	LocusLength int
	Mu          float64
	Sites       []JointSite
	// LocusTMRCA and LocusTreeLength are the TMRCA and total branch length
	// (generations) of each locus genealogy.
	LocusTMRCA      []float64
	LocusTreeLength []float64
}

// SiteKey returns the allele-frequency map key for a global position.
func SiteKey(position int) string { return fmt.Sprintf("site_%d", position) }

// PopIndex returns the index of a sampled population or -1.
func (jp *JointPanel) PopIndex(name string) int {
	for i, p := range jp.Populations {
		if p == name {
			return i
		}
	}
	return -1
}

// Frequencies returns the derived-allele frequency of every site in the joint
// panel for one population, keyed by SiteKey. Sites monomorphic in that
// population are included (frequency 0 or 1), so maps from different
// populations of the same JointPanel share exactly the same keys.
func (jp *JointPanel) Frequencies(pop string) map[string]float64 {
	i := jp.PopIndex(pop)
	if i < 0 {
		return nil
	}
	out := make(map[string]float64, len(jp.Sites))
	n := float64(jp.Haplotypes[i])
	for _, s := range jp.Sites {
		out[SiteKey(s.Position)] = float64(s.DerivedCount[i]) / n
	}
	return out
}

// DiversitySummary holds standard single-population summary statistics.
type DiversitySummary struct {
	SegregatingSites int
	Pi               float64 // mean pairwise differences per bp
	ThetaW           float64 // Watterson's θ per bp
	TajimaD          float64 // mean of per-locus Tajima's D over loci with S > 0
	// PerLocusPi and PerLocusS are per-locus π (differences per locus) and S.
	PerLocusPi []float64
	PerLocusS  []int
	// PerLocusD is Tajima's D per locus (NaN where S == 0).
	PerLocusD []float64
}

// Summary computes π, Watterson's θ and Tajima's D for one population.
func (jp *JointPanel) Summary(pop string) (DiversitySummary, error) {
	i := jp.PopIndex(pop)
	if i < 0 {
		return DiversitySummary{}, fmt.Errorf("population %q not sampled", pop)
	}
	n := jp.Haplotypes[i]
	if n < 2 {
		return DiversitySummary{}, fmt.Errorf("population %q: need at least 2 haplotypes", pop)
	}
	sum := DiversitySummary{
		PerLocusPi: make([]float64, jp.NumLoci),
		PerLocusS:  make([]int, jp.NumLoci),
		PerLocusD:  make([]float64, jp.NumLoci),
	}
	nf := float64(n)
	for _, s := range jp.Sites {
		c := s.DerivedCount[i]
		if c == 0 || c == n {
			continue
		}
		sum.PerLocusS[s.Locus]++
		sum.PerLocusPi[s.Locus] += 2 * float64(c) * float64(n-c) / (nf * (nf - 1))
	}
	a1 := harmonic(n - 1)
	totalPi, dSum, dCount := 0.0, 0.0, 0
	for l := 0; l < jp.NumLoci; l++ {
		sum.SegregatingSites += sum.PerLocusS[l]
		totalPi += sum.PerLocusPi[l]
		d := TajimasD(sum.PerLocusPi[l], sum.PerLocusS[l], n)
		sum.PerLocusD[l] = d
		if !math.IsNaN(d) {
			dSum += d
			dCount++
		}
	}
	bp := float64(jp.NumLoci * jp.LocusLength)
	sum.Pi = totalPi / bp
	sum.ThetaW = float64(sum.SegregatingSites) / a1 / bp
	if dCount > 0 {
		sum.TajimaD = dSum / float64(dCount)
	} else {
		sum.TajimaD = math.NaN()
	}
	return sum, nil
}

// harmonic returns Σ_{i=1..m} 1/i.
func harmonic(m int) float64 {
	h := 0.0
	for i := 1; i <= m; i++ {
		h += 1 / float64(i)
	}
	return h
}

// TajimasD returns Tajima's (1989) D for one locus from the mean number of
// pairwise differences pi, the number of segregating sites s and the number of
// sampled haplotypes n. It returns NaN when s == 0 or n < 4 (undefined).
func TajimasD(pi float64, s, n int) float64 {
	if s == 0 || n < 4 {
		return math.NaN()
	}
	nf := float64(n)
	a1 := harmonic(n - 1)
	a2 := 0.0
	for i := 1; i < n; i++ {
		a2 += 1 / float64(i*i)
	}
	b1 := (nf + 1) / (3 * (nf - 1))
	b2 := 2 * (nf*nf + nf + 3) / (9 * nf * (nf - 1))
	c1 := b1 - 1/a1
	c2 := b2 - (nf+2)/(a1*nf) + a2/(a1*a1)
	e1 := c1 / a1
	e2 := c2 / (a1*a1 + a2)
	sf := float64(s)
	v := e1*sf + e2*sf*(sf-1)
	if v <= 0 {
		return math.NaN()
	}
	return (pi - sf/a1) / math.Sqrt(v)
}

// SimulateJointPanel simulates nLoci unlinked loci under the multi-population
// demography d. samples maps population name → number of sampled haplotypes.
// Each locus gets one genealogy spanning all sampled populations, so every
// population is scored at the same sites. The simulation is deterministic for
// a given seed.
//
// Algorithm (backwards in time): lineages are grouped by population; within a
// population with k lineages coalescence occurs at rate k(k−1)/(2N(t)), where
// N(t) is the population's piecewise constant/exponential size history.
// Waiting times are drawn exactly by inverting the cumulative hazard across
// epoch boundaries. Competing populations draw independent waiting times and
// the earliest wins (valid because the processes are independent Poisson
// processes with deterministic intensities). At a population's SplitTime its
// remaining lineages move into SplitFrom. Mutations are placed on each branch
// as a Poisson process with rate μ·LocusLength per generation (infinite
// sites); the derived-allele count of each population is the number of its
// sampled leaves below the mutated branch.
func SimulateJointPanel(d Demography, samples map[string]int, nLoci int, seed uint64) (*JointPanel, error) {
	if err := d.Validate(); err != nil {
		return nil, err
	}
	if nLoci < 1 {
		return nil, fmt.Errorf("nLoci must be at least 1")
	}
	locusLen := d.LocusLength
	if locusLen == 0 {
		locusLen = DefaultLocusLength
	}

	popIdx := make(map[string]int, len(d.Pops))
	for i, p := range d.Pops {
		popIdx[p.Name] = i
	}

	names := make([]string, 0, len(samples))
	total := 0
	for name, n := range samples {
		if _, ok := popIdx[name]; !ok {
			return nil, fmt.Errorf("sampled population %q is not in the demography", name)
		}
		if n < 0 {
			return nil, fmt.Errorf("negative sample size for %q", name)
		}
		if n == 0 {
			continue
		}
		names = append(names, name)
		total += n
	}
	if total < 2 {
		return nil, fmt.Errorf("need at least 2 sampled haplotypes in total, got %d", total)
	}
	sort.Strings(names)

	jp := &JointPanel{
		Populations:     names,
		Haplotypes:      make([]int, len(names)),
		NumLoci:         nLoci,
		LocusLength:     locusLen,
		Mu:              d.Mu,
		LocusTMRCA:      make([]float64, nLoci),
		LocusTreeLength: make([]float64, nLoci),
	}
	for i, name := range names {
		jp.Haplotypes[i] = samples[name]
	}

	// Per-demography-population data.
	nPops := len(d.Pops)
	epochsOf := make([][]epoch, nPops)
	parent := make([]int, nPops)
	splitAt := make([]float64, nPops)
	for i, p := range d.Pops {
		epochsOf[i] = p.epochs()
		if p.SplitFrom == "" {
			parent[i] = -1
			splitAt[i] = math.Inf(1)
		} else {
			parent[i] = popIdx[p.SplitFrom]
			splitAt[i] = float64(p.SplitTime)
		}
	}

	// Leaves: sample population index (into names) and haplotype index within
	// that population.
	leafSample := make([]int, 0, total)
	leafHap := make([]int, 0, total)
	leafDemPop := make([]int, 0, total)
	for si, name := range names {
		for h := 0; h < samples[name]; h++ {
			leafSample = append(leafSample, si)
			leafHap = append(leafHap, h)
			leafDemPop = append(leafDemPop, popIdx[name])
		}
	}

	rng := rand.New(rand.NewPCG(seed, seed^0x9e3779b97f4a7c15)) // #nosec G404 -- simulation, not crypto
	sim := &jointLocusSim{
		nPops:      nPops,
		epochsOf:   epochsOf,
		parent:     parent,
		splitAt:    splitAt,
		leafDemPop: leafDemPop,
		nSampled:   len(names),
		leafSample: leafSample,
		rng:        rng,
	}

	for l := 0; l < nLoci; l++ {
		sim.simulateTree()
		jp.LocusTMRCA[l] = sim.time[len(sim.time)-1]
		treeLen := 0.0
		for v := 0; v < len(sim.time)-1; v++ {
			treeLen += sim.time[sim.par[v]] - sim.time[v]
		}
		jp.LocusTreeLength[l] = treeLen
		sim.addMutations(jp, l, d.Mu, locusLen, treeLen, leafHap)
	}
	return jp, nil
}

// jointLocusSim holds reusable state for simulating one locus genealogy.
type jointLocusSim struct {
	nPops      int
	epochsOf   [][]epoch
	parent     []int
	splitAt    []float64
	leafDemPop []int
	nSampled   int
	leafSample []int
	rng        *rand.Rand

	// Tree (rebuilt per locus). Node v: time[v], par[v] (−1 for root),
	// children kids[v]; leaves are 0..len(leafDemPop)-1; the root is last.
	time []float64
	par  []int
	kids [][2]int
}

func (s *jointLocusSim) simulateTree() {
	nLeaves := len(s.leafDemPop)
	s.time = s.time[:0]
	s.par = s.par[:0]
	s.kids = s.kids[:0]
	for v := 0; v < nLeaves; v++ {
		s.time = append(s.time, 0)
		s.par = append(s.par, -1)
		s.kids = append(s.kids, [2]int{-1, -1})
	}

	lineages := make([][]int, s.nPops)
	for v, p := range s.leafDemPop {
		lineages[p] = append(lineages[p], v)
	}
	active := make([]bool, s.nPops)
	for i := range active {
		active[i] = true
	}

	t := 0.0
	remaining := nLeaves
	for remaining > 1 {
		// Next split among active populations that still hold lineages or
		// whose parents might receive them.
		nextSplit := math.Inf(1)
		for p := 0; p < s.nPops; p++ {
			if active[p] && s.parent[p] >= 0 && s.splitAt[p] < nextSplit {
				nextSplit = s.splitAt[p]
			}
		}

		best, bestPop := math.Inf(1), -1
		if nextSplit > t {
			for p := 0; p < s.nPops; p++ {
				k := len(lineages[p])
				if !active[p] || k < 2 {
					continue
				}
				rate := float64(k*(k-1)) / 4 // k(k−1)/2 pairs, each at 1/(2N)
				tc := nextEventTime(s.epochsOf[p], rate, t, nextSplit, s.rng.ExpFloat64())
				if tc < best {
					best, bestPop = tc, p
				}
			}
		}

		if bestPop >= 0 && best < nextSplit {
			t = best
			ls := lineages[bestPop]
			k := len(ls)
			i := s.rng.IntN(k)
			j := s.rng.IntN(k - 1)
			if j >= i {
				j++
			}
			a, b := ls[i], ls[j]
			v := len(s.time)
			s.time = append(s.time, t)
			s.par = append(s.par, -1)
			s.kids = append(s.kids, [2]int{a, b})
			s.par[a], s.par[b] = v, v
			// Remove i and j (swap-remove, larger index first), add v.
			if i < j {
				i, j = j, i
			}
			ls[i] = ls[k-1]
			ls = ls[:k-1]
			ls[j] = ls[len(ls)-1]
			ls = ls[:len(ls)-1]
			lineages[bestPop] = append(ls, v)
			remaining--
			continue
		}

		if math.IsInf(nextSplit, 1) {
			// Cannot happen for a validated demography: the root population
			// always has an infinite-hazard final epoch.
			panic("structured coalescent: lineages can never coalesce")
		}
		// Apply all splits at nextSplit (children may chain into parents that
		// split at the same time; iterate until stable).
		t = nextSplit
		for changed := true; changed; {
			changed = false
			for p := 0; p < s.nPops; p++ {
				if active[p] && s.parent[p] >= 0 && s.splitAt[p] <= t {
					q := s.parent[p]
					lineages[q] = append(lineages[q], lineages[p]...)
					lineages[p] = nil
					active[p] = false
					changed = true
				}
			}
		}
	}
}

// addMutations drops infinite-sites mutations on the current tree.
func (s *jointLocusSim) addMutations(jp *JointPanel, locus int, mu float64, locusLen int, treeLen float64, leafHap []int) {
	if mu <= 0 || treeLen <= 0 {
		return
	}
	nMut := poissonExact(s.rng, mu*float64(locusLen)*treeLen)
	if nMut == 0 {
		return
	}
	nBranches := len(s.time) - 1 // every node but the root
	cum := make([]float64, nBranches)
	acc := 0.0
	for v := 0; v < nBranches; v++ {
		acc += s.time[s.par[v]] - s.time[v]
		cum[v] = acc
	}
	used := make(map[int]bool, nMut)
	bases := [4]byte{'A', 'C', 'G', 'T'}

	for m := 0; m < nMut; m++ {
		r := s.rng.Float64() * acc
		v := sort.SearchFloat64s(cum, r)
		if v >= nBranches {
			v = nBranches - 1
		}
		// Infinite sites: one mutation per position.
		pos := s.rng.IntN(locusLen)
		for tries := 0; used[pos] && tries < 4*locusLen; tries++ {
			pos = s.rng.IntN(locusLen)
		}
		if used[pos] {
			continue // locus saturated; drop the mutation
		}
		used[pos] = true

		anc := s.rng.IntN(4)
		der := (anc + 1 + s.rng.IntN(3)) % 4
		lo, hi := s.time[v], s.time[s.par[v]]

		site := JointSite{
			Locus:        locus,
			Position:     locus*locusLen + pos,
			Ancestral:    bases[anc],
			Derived:      bases[der],
			Age:          lo + s.rng.Float64()*(hi-lo),
			DerivedCount: make([]int, s.nSampled),
			Genotypes:    make([][]int, s.nSampled),
			Carriers:     make([][]byte, s.nSampled),
		}
		for i, n := range jp.Haplotypes {
			if n%2 == 0 {
				site.Genotypes[i] = make([]int, n/2)
			}
			site.Carriers[i] = make([]byte, n)
		}
		s.forLeavesBelow(v, func(leaf int) {
			si := s.leafSample[leaf]
			site.DerivedCount[si]++
			site.Carriers[si][leafHap[leaf]] = 1
			if g := site.Genotypes[si]; g != nil {
				g[leafHap[leaf]/2]++
			}
		})
		jp.Sites = append(jp.Sites, site)
	}
	// Keep sites ordered by position within the locus.
	start := len(jp.Sites)
	for start > 0 && jp.Sites[start-1].Locus == locus {
		start--
	}
	ls := jp.Sites[start:]
	sort.Slice(ls, func(a, b int) bool { return ls[a].Position < ls[b].Position })
}

// forLeavesBelow calls fn for every leaf in the subtree rooted at v.
func (s *jointLocusSim) forLeavesBelow(v int, fn func(int)) {
	stack := []int{v}
	for len(stack) > 0 {
		u := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		k := s.kids[u]
		if k[0] < 0 {
			fn(u)
			continue
		}
		stack = append(stack, k[0], k[1])
	}
}

// poissonExact draws a Poisson(lambda) variate exactly, as a sum of Poisson
// variates with mean <= 20 drawn by Knuth's multiplication method.
func poissonExact(rng *rand.Rand, lambda float64) int {
	if !(lambda > 0) {
		return 0
	}
	k := 0
	for lambda > 0 {
		l := math.Min(lambda, 20)
		lambda -= l
		limit := math.Exp(-l)
		p := rng.Float64()
		for p > limit {
			k++
			p *= rng.Float64()
		}
	}
	return k
}
