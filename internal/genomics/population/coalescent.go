// Package population simulates neutral population-genetic variation with the
// coalescent.
//
// SimulateJointPanel runs a structured coalescent over a population tree
// (splits and piecewise-exponential size histories; no migration, admixture
// or within-locus recombination) on unlinked loci, so that every population
// is scored at the same sites. PopulationPanel and NewJointPanels wrap it
// using the demographic presets for YRI, CEU, CHB and MXL in presets.go;
// other population names fall back to a constant-size model.
// CoalescentSimulator is the older single-population simulator.
package population

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sort"
	"sync"
	"time"
)

// CoalescentSimulator simulates single-population coalescent genealogies
// under a piecewise population-size history, with mutations placed on
// branches by length. Recombination events are recorded but do not change
// the genealogy, and migration and admixture events are ignored. Its random
// source is seeded from the clock until SetSeed is called.
type CoalescentSimulator struct {
	effectivePopSize   int
	mutationRate       float64
	recombinationRate  float64
	demographicHistory []DemographicEvent
	genealogies        map[string]*Genealogy
	rng                *rand.Rand
	mu                 sync.RWMutex
}

// DemographicEvent represents a demographic event in population history
type DemographicEvent struct {
	Type        string             // "bottleneck", "expansion", "size_change", "migration", "admixture" (see size_history.go)
	Time        int                // Generations ago
	Parameters  map[string]float64 // Event-specific parameters
	Populations []string           // Affected populations
	Duration    int                // Duration in generations
}

// Genealogy represents a genealogical tree for a genomic region
type Genealogy struct {
	Nodes               []*CoalescentNode
	Root                *CoalescentNode
	LeafNodes           []*CoalescentNode
	TotalTime           float64
	MutationEvents      []MutationEvent
	RecombinationEvents []RecombinationEvent
}

// CoalescentNode represents a node in the coalescent tree
type CoalescentNode struct {
	ID           int
	Time         float64
	Parent       *CoalescentNode
	Children     []*CoalescentNode
	Population   string
	Mutations    []Mutation
	IsLeaf       bool
	LineageCount int
}

// MutationEvent represents a mutation event in the genealogy
type MutationEvent struct {
	Position int
	Time     float64
	Lineage  *CoalescentNode
	Mutation Mutation
	Type     string // "transition" or "transversion"
}

// RecombinationEvent represents a recombination event
type RecombinationEvent struct {
	Position int
	Time     float64
	Lineages []*CoalescentNode
}

// Mutation represents a single mutation
type Mutation struct {
	Position  int
	Ancestral byte
	Derived   byte
	Time      float64
	Frequency float64
	IsFixed   bool
}

// AlleleFrequencySpectrum represents the site frequency spectrum
type AlleleFrequencySpectrum struct {
	Frequencies    []float64 // Frequency bins
	Counts         []int     // Count of sites in each bin
	TotalSites     int
	ThetaWatterson float64
	ThetaPi        float64
	TajimaD        float64
	FayWuH         float64
}

// NewCoalescentSimulator creates a new coalescent simulator
func NewCoalescentSimulator(effectivePopSize int, mutationRate, recombinationRate float64) *CoalescentSimulator {
	return &CoalescentSimulator{
		effectivePopSize:   effectivePopSize,
		mutationRate:       mutationRate,
		recombinationRate:  recombinationRate,
		demographicHistory: make([]DemographicEvent, 0),
		genealogies:        make(map[string]*Genealogy),
		rng:                rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
	}
}

// AddDemographicEvent adds a demographic event to the simulation
func (s *CoalescentSimulator) AddDemographicEvent(event DemographicEvent) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.demographicHistory = append(s.demographicHistory, event)

	// Sort events by time (most recent first)
	sort.Slice(s.demographicHistory, func(i, j int) bool {
		return s.demographicHistory[i].Time < s.demographicHistory[j].Time
	})
}

// SimulateGenealogy simulates a coalescent genealogy for a genomic region
func (s *CoalescentSimulator) SimulateGenealogy(regionID string, sampleSize int, regionLength int) (*Genealogy, error) {
	s.mu.Lock()
	defer s.mu.Unlock()

	if sampleSize < 2 {
		return nil, fmt.Errorf("sample size must be at least 2")
	}

	genealogy := &Genealogy{
		Nodes:               make([]*CoalescentNode, 0),
		LeafNodes:           make([]*CoalescentNode, 0),
		MutationEvents:      make([]MutationEvent, 0),
		RecombinationEvents: make([]RecombinationEvent, 0),
	}

	// Create initial sample (leaf nodes)
	activeLineages := make([]*CoalescentNode, sampleSize)
	for i := 0; i < sampleSize; i++ {
		node := &CoalescentNode{
			ID:           i,
			Time:         0.0,
			Children:     make([]*CoalescentNode, 0),
			Population:   "modern",
			Mutations:    make([]Mutation, 0),
			IsLeaf:       true,
			LineageCount: 1,
		}
		activeLineages[i] = node
		genealogy.Nodes = append(genealogy.Nodes, node)
		genealogy.LeafNodes = append(genealogy.LeafNodes, node)
	}

	currentTime := 0.0
	nodeCounter := sampleSize

	// Piecewise size history N(t) from the effective size and demographic
	// events (see size_history.go).
	epochs := buildEpochs(float64(s.effectivePopSize), s.demographicHistory)

	// Simulate coalescent process
	for len(activeLineages) > 1 {
		// With n lineages there are n(n-1)/2 pairs, each coalescing at rate
		// 1/(2N(t)) for a diploid population of N(t) individuals, so E[T2] = 2N.
		// The waiting time is drawn exactly across epoch boundaries and through
		// exponential-growth epochs.
		n := len(activeLineages)
		rate := float64(n*(n-1)) / 4.0
		currentTime = nextEventTime(epochs, rate, currentTime, math.Inf(1), s.rng.ExpFloat64())
		if math.IsInf(currentTime, 1) {
			return nil, fmt.Errorf("demographic history prevents coalescence")
		}

		// Choose two lineages to coalesce
		i := s.rng.IntN(n)
		j := s.rng.IntN(n - 1)
		if j >= i {
			j++
		}

		// Create new internal node
		newNode := &CoalescentNode{
			ID:           nodeCounter,
			Time:         currentTime,
			Children:     []*CoalescentNode{activeLineages[i], activeLineages[j]},
			Population:   s.getPopulationAtTime(currentTime),
			Mutations:    make([]Mutation, 0),
			IsLeaf:       false,
			LineageCount: activeLineages[i].LineageCount + activeLineages[j].LineageCount,
		}

		// Set parent relationships
		activeLineages[i].Parent = newNode
		activeLineages[j].Parent = newNode

		genealogy.Nodes = append(genealogy.Nodes, newNode)
		nodeCounter++

		// Remove coalesced lineages and add new one
		newActiveLineages := make([]*CoalescentNode, 0, n-1)
		for k := 0; k < n; k++ {
			if k != i && k != j {
				newActiveLineages = append(newActiveLineages, activeLineages[k])
			}
		}
		newActiveLineages = append(newActiveLineages, newNode)
		activeLineages = newActiveLineages
	}

	// Set root
	genealogy.Root = activeLineages[0]
	genealogy.TotalTime = currentTime

	// Add mutations to the genealogy
	s.addMutationsToGenealogy(genealogy, regionLength)

	// Add recombination events if needed
	if s.recombinationRate > 0 {
		s.addRecombinationEvents(genealogy, regionLength)
	}

	// Store genealogy
	s.genealogies[regionID] = genealogy

	return genealogy, nil
}

// getAdjustedPopSize returns the effective population size N(t) at time t
// (generations ago) implied by the demographic history. Bottlenecks last for
// their Duration (a Duration of 0 is a persistent size change), expansions
// shrink the population going back in time (exponentially over Duration, or
// instantaneously by expansion_factor), and event types that do not change
// size are ignored. See size_history.go for the full semantics.
func (s *CoalescentSimulator) getAdjustedPopSize(currentTime float64) int {
	n, _ := sizeAt(float64(s.effectivePopSize), sortedEvents(s.demographicHistory), currentTime)
	return int(math.Round(n))
}

// getPopulationAtTime returns the population name for a given time (generations
// ago) by consulting the demographic history.  The most recent (smallest Time)
// event whose range contains t and that names at least one population wins.
// Falls back to "ancestral" when no event matches.
func (s *CoalescentSimulator) getPopulationAtTime(t float64) string {
	bestTime := -1
	bestPop := "ancestral"
	for _, event := range s.demographicHistory {
		if len(event.Populations) == 0 {
			continue
		}
		start := float64(event.Time)
		// A point event (Duration==0) applies from its Time into the past
		// (i.e., t >= start).  A ranged event applies within [start, start+Duration].
		end := start
		if event.Duration > 0 {
			end = start + float64(event.Duration)
		}
		if t >= start && t <= end {
			// Prefer the most recent (smallest Time) event.
			if bestTime < 0 || event.Time < bestTime {
				bestTime = event.Time
				bestPop = event.Populations[0]
			}
		}
	}
	return bestPop
}

// addMutationsToGenealogy adds mutations to the genealogy based on the mutation rate
func (s *CoalescentSimulator) addMutationsToGenealogy(genealogy *Genealogy, regionLength int) {
	totalBranchLength := s.calculateTotalBranchLength(genealogy)
	expectedMutations := s.mutationRate * float64(regionLength) * totalBranchLength

	// Poisson number of mutations
	numMutations := s.poissonRandom(expectedMutations)

	for i := 0; i < numMutations; i++ {
		// Choose random position
		position := s.rng.IntN(regionLength)

		// Choose random branch (weighted by branch length)
		branch := s.chooseBranchWeighted(genealogy)
		if branch == nil {
			continue
		}

		// Choose mutation time uniformly along branch
		branchLength := s.getBranchLength(branch)
		mutationTime := branch.Time - s.rng.Float64()*branchLength

		// Create mutation
		mutation := Mutation{
			Position:  position,
			Ancestral: s.getRandomBase(),
			Derived:   s.getRandomBase(),
			Time:      mutationTime,
			Frequency: 0.0, // Will be calculated later
			IsFixed:   false,
		}

		// Ensure ancestral and derived are different
		for mutation.Derived == mutation.Ancestral {
			mutation.Derived = s.getRandomBase()
		}

		// Add mutation to branch
		branch.Mutations = append(branch.Mutations, mutation)

		// Create mutation event
		mutationEvent := MutationEvent{
			Position: position,
			Time:     mutationTime,
			Lineage:  branch,
			Mutation: mutation,
			Type:     s.getMutationType(mutation.Ancestral, mutation.Derived),
		}

		genealogy.MutationEvents = append(genealogy.MutationEvents, mutationEvent)
	}

	// Calculate allele frequencies
	s.calculateAlleleFrequencies(genealogy)
}

// calculateTotalBranchLength calculates the total branch length of the genealogy
func (s *CoalescentSimulator) calculateTotalBranchLength(genealogy *Genealogy) float64 {
	totalLength := 0.0

	for _, node := range genealogy.Nodes {
		if node.Parent != nil {
			branchLength := node.Parent.Time - node.Time
			totalLength += branchLength
		}
	}

	return totalLength
}

// chooseBranchWeighted chooses a branch weighted by its length
func (s *CoalescentSimulator) chooseBranchWeighted(genealogy *Genealogy) *CoalescentNode {
	// Calculate branch lengths
	branchWeights := make([]float64, 0)
	nodes := make([]*CoalescentNode, 0)

	for _, node := range genealogy.Nodes {
		if node.Parent != nil {
			branchLength := node.Parent.Time - node.Time
			branchWeights = append(branchWeights, branchLength)
			nodes = append(nodes, node)
		}
	}

	if len(branchWeights) == 0 {
		return nil
	}

	// Choose weighted random branch
	totalWeight := 0.0
	for _, weight := range branchWeights {
		totalWeight += weight
	}

	r := s.rng.Float64() * totalWeight
	cumulative := 0.0

	for i, weight := range branchWeights {
		cumulative += weight
		if r <= cumulative {
			return nodes[i]
		}
	}

	return nodes[len(nodes)-1]
}

// getBranchLength returns the length of a branch
func (s *CoalescentSimulator) getBranchLength(node *CoalescentNode) float64 {
	if node.Parent == nil {
		return 0.0
	}
	return node.Parent.Time - node.Time
}

// getRandomBase returns a random DNA base
func (s *CoalescentSimulator) getRandomBase() byte {
	bases := []byte{'A', 'T', 'G', 'C'}
	return bases[s.rng.IntN(4)]
}

// getMutationType determines the type of mutation
func (s *CoalescentSimulator) getMutationType(ancestral, derived byte) string {
	transitions := map[string]bool{
		"A->G": true, "G->A": true,
		"C->T": true, "T->C": true,
	}

	mutKey := fmt.Sprintf("%c->%c", ancestral, derived)
	if transitions[mutKey] {
		return "transition"
	}
	return "transversion"
}

// calculateAlleleFrequencies calculates allele frequencies for all mutations
func (s *CoalescentSimulator) calculateAlleleFrequencies(genealogy *Genealogy) {
	for i := range genealogy.MutationEvents {
		mutation := &genealogy.MutationEvents[i].Mutation
		frequency := s.calculateMutationFrequency(genealogy, mutation)
		mutation.Frequency = frequency
		mutation.IsFixed = frequency == 1.0
	}
}

// calculateMutationFrequency calculates the frequency of a mutation in the sample
func (s *CoalescentSimulator) calculateMutationFrequency(genealogy *Genealogy, mutation *Mutation) float64 {
	// Count how many leaf nodes inherit this mutation
	inheritingLeaves := 0

	for _, leaf := range genealogy.LeafNodes {
		if s.nodeInheritsMutation(leaf, mutation) {
			inheritingLeaves++
		}
	}

	return float64(inheritingLeaves) / float64(len(genealogy.LeafNodes))
}

// nodeInheritsMutation checks if a node inherits a specific mutation
func (s *CoalescentSimulator) nodeInheritsMutation(node *CoalescentNode, mutation *Mutation) bool {
	current := node

	for current != nil {
		for _, nodeMutation := range current.Mutations {
			if nodeMutation.Position == mutation.Position && nodeMutation.Time == mutation.Time {
				return true
			}
		}
		current = current.Parent
	}

	return false
}

// addRecombinationEvents records Poisson-distributed recombination events on
// the genealogy. They are bookkeeping only: the tree, and so the mutations and
// allele frequencies, are the same as without them (no ancestral
// recombination graph is built).
func (s *CoalescentSimulator) addRecombinationEvents(genealogy *Genealogy, regionLength int) {
	totalTime := genealogy.TotalTime
	expectedRecombinations := s.recombinationRate * float64(regionLength) * totalTime

	numRecombinations := s.poissonRandom(expectedRecombinations)

	for i := 0; i < numRecombinations; i++ {
		position := s.rng.IntN(regionLength)
		time := s.rng.Float64() * totalTime

		// Find lineages active at this time
		activeLineages := s.getActiveLineagesAtTime(genealogy, time)

		if len(activeLineages) >= 2 {
			// Choose two lineages for recombination
			lineage1 := activeLineages[s.rng.IntN(len(activeLineages))]
			lineage2 := activeLineages[s.rng.IntN(len(activeLineages))]

			recombEvent := RecombinationEvent{
				Position: position,
				Time:     time,
				Lineages: []*CoalescentNode{lineage1, lineage2},
			}

			genealogy.RecombinationEvents = append(genealogy.RecombinationEvents, recombEvent)
		}
	}
}

// getActiveLineagesAtTime returns lineages active at a given time
func (s *CoalescentSimulator) getActiveLineagesAtTime(genealogy *Genealogy, time float64) []*CoalescentNode {
	activeLineages := make([]*CoalescentNode, 0)

	for _, node := range genealogy.Nodes {
		if node.Time <= time {
			if node.Parent == nil || node.Parent.Time > time {
				activeLineages = append(activeLineages, node)
			}
		}
	}

	return activeLineages
}

// SimulateAlleleFrequencies simulates one genealogy over a region of numSites
// bp and returns the sample frequency of each mutation keyed "site_<pos>".
// Mutations that hit the same position share a key, so the last one wins.
func (s *CoalescentSimulator) SimulateAlleleFrequencies(numSites int, sampleSize int) (map[string]float64, error) {
	frequencies := make(map[string]float64)

	// Simulate genealogy for the region
	genealogy, err := s.SimulateGenealogy("simulation", sampleSize, numSites)
	if err != nil {
		return nil, err
	}

	// Extract allele frequencies from mutations
	for _, mutationEvent := range genealogy.MutationEvents {
		site := fmt.Sprintf("site_%d", mutationEvent.Position)
		frequencies[site] = mutationEvent.Mutation.Frequency
	}

	return frequencies, nil
}

// ComputeSiteFrequencySpectrum computes the site frequency spectrum
func (s *CoalescentSimulator) ComputeSiteFrequencySpectrum(genealogy *Genealogy, sampleSize int) *AlleleFrequencySpectrum {
	spectrum := &AlleleFrequencySpectrum{
		Frequencies: make([]float64, sampleSize),
		Counts:      make([]int, sampleSize),
		TotalSites:  len(genealogy.MutationEvents),
	}

	// Fill frequency bins
	for i := 0; i < sampleSize; i++ {
		spectrum.Frequencies[i] = float64(i+1) / float64(sampleSize)
	}

	// Count mutations in each frequency bin
	for _, mutationEvent := range genealogy.MutationEvents {
		frequency := mutationEvent.Mutation.Frequency
		// Round, not truncate: k/n*n can fall just below k in floating point.
		binIndex := int(math.Round(frequency*float64(sampleSize))) - 1
		if binIndex < 0 {
			binIndex = 0
		}
		if binIndex >= sampleSize {
			binIndex = sampleSize - 1
		}
		spectrum.Counts[binIndex]++
	}

	// Calculate summary statistics (require at least 2 samples; all helpers divide by sampleSize-1)
	if sampleSize < 2 {
		return spectrum // ThetaWatterson/ThetaPi/TajimaD/FayWuH remain 0.0
	}
	spectrum.ThetaWatterson = s.calculateThetaWatterson(spectrum, sampleSize)
	spectrum.ThetaPi = s.calculateThetaPi(spectrum, sampleSize)
	spectrum.TajimaD = s.calculateTajimaD(spectrum.ThetaPi, spectrum.ThetaWatterson, sampleSize)
	spectrum.FayWuH = s.calculateFayWuH(spectrum, sampleSize)

	return spectrum
}

// calculateThetaWatterson calculates Watterson's theta
func (s *CoalescentSimulator) calculateThetaWatterson(spectrum *AlleleFrequencySpectrum, sampleSize int) float64 {
	segregatingSites := 0
	for _, count := range spectrum.Counts {
		segregatingSites += count
	}

	harmonicNumber := 0.0
	for i := 1; i < sampleSize; i++ {
		harmonicNumber += 1.0 / float64(i)
	}

	return float64(segregatingSites) / harmonicNumber
}

// calculateThetaPi returns Tajima's (1983) θ_π, the mean number of pairwise
// differences: Σ_i i(n−i)ξ_i / C(n,2) = Σ 2f(1−f)ξ_f · n/(n−1) with f = i/n.
// Under the constant-size neutral model E[θ_π] = θ.
func (s *CoalescentSimulator) calculateThetaPi(spectrum *AlleleFrequencySpectrum, sampleSize int) float64 {
	pi := 0.0

	for i, count := range spectrum.Counts {
		frequency := spectrum.Frequencies[i]
		pi += 2 * frequency * (1 - frequency) * float64(count)
	}

	return pi * float64(sampleSize) / float64(sampleSize-1)
}

// calculateTajimaD calculates Tajima's D statistic
func (s *CoalescentSimulator) calculateTajimaD(thetaPi, thetaWatterson float64, sampleSize int) float64 {
	if thetaWatterson == 0 {
		return 0.0
	}

	a1 := 0.0
	a2 := 0.0
	for i := 1; i < sampleSize; i++ {
		a1 += 1.0 / float64(i)
		a2 += 1.0 / float64(i*i)
	}

	b1 := float64(sampleSize+1) / (3.0 * float64(sampleSize-1))
	b2 := 2.0 * (float64(sampleSize*sampleSize) + float64(sampleSize) + 3.0) / (9.0 * float64(sampleSize) * float64(sampleSize-1))

	c1 := b1 - 1.0/a1
	c2 := b2 - float64(sampleSize+2)/(a1*float64(sampleSize)) + a2/(a1*a1)

	e1 := c1 / a1
	e2 := c2 / (a1*a1 + a2)

	// Var(π − θ_W) = e1·S + e2·S·(S−1), with S = θ_W·a1 segregating sites.
	segregating := thetaWatterson * a1
	variance := e1*segregating + e2*segregating*(segregating-1)

	if variance <= 0 {
		return 0.0 // undefined when variance is non-positive; guards math.Sqrt(negative) → NaN
	}

	return (thetaPi - thetaWatterson) / math.Sqrt(variance)
}

// calculateFayWuH returns Fay and Wu's (2000) unnormalized H = θ_π − θ_H, with
// θ_H = Σ_i i²ξ_i / C(n,2) = Σ 2f²ξ_f · n/(n−1). E[H] = 0 under neutrality.
func (s *CoalescentSimulator) calculateFayWuH(spectrum *AlleleFrequencySpectrum, sampleSize int) float64 {
	thetaH := 0.0

	for i, count := range spectrum.Counts {
		frequency := spectrum.Frequencies[i]
		thetaH += frequency * frequency * float64(count)
	}

	thetaH *= 2.0 * float64(sampleSize) / float64(sampleSize-1)

	return spectrum.ThetaPi - thetaH
}

// poissonRandom generates a Poisson random number
func (s *CoalescentSimulator) poissonRandom(lambda float64) int {
	if lambda < 30.0 {
		// Use Knuth's algorithm for small lambda
		L := math.Exp(-lambda)
		k := 0
		p := 1.0

		for p > L {
			k++
			p *= s.rng.Float64()
		}

		return k - 1
	} else {
		// Use normal approximation for large lambda; clamp to 0 (Poisson is non-negative)
		normal := s.rng.NormFloat64()
		result := int(lambda + math.Sqrt(lambda)*normal + 0.5)
		if result < 0 {
			return 0
		}
		return result
	}
}

// SetSeed sets the random seed for reproducible simulations
func (s *CoalescentSimulator) SetSeed(seed int64) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.rng = rand.New(rand.NewPCG(uint64(seed), 0)) // #nosec G115 G404
}

// GetSimulationStats returns statistics about the simulation
func (s *CoalescentSimulator) GetSimulationStats() map[string]interface{} {
	s.mu.RLock()
	defer s.mu.RUnlock()

	stats := map[string]interface{}{
		"effective_population_size": s.effectivePopSize,
		"mutation_rate":             s.mutationRate,
		"recombination_rate":        s.recombinationRate,
		"demographic_events":        len(s.demographicHistory),
		"genealogies_simulated":     len(s.genealogies),
	}

	return stats
}

// ValidateParameters validates simulation parameters
func (s *CoalescentSimulator) ValidateParameters() error {
	if s.effectivePopSize <= 0 {
		return fmt.Errorf("effective population size must be positive")
	}

	if s.mutationRate < 0 || s.mutationRate > 1 {
		return fmt.Errorf("mutation rate must be between 0 and 1")
	}

	if s.recombinationRate < 0 || s.recombinationRate > 1 {
		return fmt.Errorf("recombination rate must be between 0 and 1")
	}

	return nil
}

// nextSeed draws a seed for a derived simulation from the simulator's RNG.
func (s *CoalescentSimulator) nextSeed() uint64 {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.rng.Uint64()
}

// branchParameters returns a snapshot of the effective size, mutation rate
// and demographic events.
func (s *CoalescentSimulator) branchParameters() (int, float64, []DemographicEvent) {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.effectivePopSize, s.mutationRate, append([]DemographicEvent(nil), s.demographicHistory...)
}

// GetDemographicEvents returns the demographic events in the simulation
func (s *CoalescentSimulator) GetDemographicEvents() []DemographicEvent {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.demographicHistory
}

// PopulationPanel represents a panel of population genetic data
type PopulationPanel struct {
	PopulationName string
	// SampleSize is the number of diploid individuals; panels are simulated
	// from 2·SampleSize haplotypes.
	SampleSize int
	// Simulator carries this population's branch parameters: its effective
	// size (for presets, the pre-growth Ne of the population's own branch),
	// mutation rate, extra demographic events (SetDemographicEvent) and the
	// RNG used by LoadFromSimulation (seed it with Simulator.SetSeed).
	Simulator *CoalescentSimulator
	// AlleleFrequencies maps SiteKey(position) to the derived-allele
	// frequency. For panels built together by NewJointPanels every panel has
	// the same keys (sites monomorphic in a population have frequency 0 or 1).
	AlleleFrequencies map[string]float64
	Demographics      map[string]interface{}
	// Variants lists the sites segregating in this panel, with genotypes
	// formed from the simulated haplotypes (individual k = haplotypes 2k, 2k+1).
	Variants []Variant
	preset   string // preset population used for the demographic model ("" = none)
	mu       sync.RWMutex
}

// Variant represents a genetic variant in the population panel
type Variant struct {
	Position      int
	Chromosome    string
	ReferenceBase byte
	AlternateBase byte
	Frequency     float64
	Genotypes     []int // 0=homozygous ref, 1=heterozygous, 2=homozygous alt
	Quality       float64
}

// Default single-population model for populations without a preset.
const (
	defaultPanelNe                = 10000
	defaultPanelMutationRate      = 1e-8
	defaultPanelRecombinationRate = 1e-8
)

// NewPopulationPanel creates a population panel. Populations with a
// demographic preset (YRI, CEU, CHB, MXL and the aliases accepted by
// ResolvePreset) are simulated under the joint preset model in presets.go.
// Any other name falls back to a constant-size model (Ne 10000, μ 1e-8) with
// no demographic events; use HasPreset to detect this, or
// NewPresetPopulationPanel to get an error instead.
func NewPopulationPanel(populationName string, sampleSize int) *PopulationPanel {
	pp := &PopulationPanel{
		PopulationName:    populationName,
		SampleSize:        sampleSize,
		AlleleFrequencies: make(map[string]float64),
		Demographics:      make(map[string]interface{}),
		Variants:          make([]Variant, 0),
	}
	spec, err := LookupPreset(populationName)
	if err != nil {
		pp.Simulator = NewCoalescentSimulator(defaultPanelNe, defaultPanelMutationRate, defaultPanelRecombinationRate)
		pp.Demographics["model"] = "default constant-size (no preset)"
		pp.Demographics["effective_population_size"] = float64(defaultPanelNe)
		pp.Demographics["mutation_rate"] = defaultPanelMutationRate
		return pp
	}
	pp.preset = spec.Name
	pp.Simulator = NewCoalescentSimulator(int(spec.Ne), PresetMutationRate, PresetRecombinationRate)
	pp.Demographics = presetDemographics(spec)
	return pp
}

// NewPresetPopulationPanel is like NewPopulationPanel but returns an error for
// populations that have no demographic preset.
func NewPresetPopulationPanel(populationName string, sampleSize int) (*PopulationPanel, error) {
	if _, err := LookupPreset(populationName); err != nil {
		return nil, err
	}
	return NewPopulationPanel(populationName, sampleSize), nil
}

// HasPreset reports whether the panel is simulated under a demographic preset.
func (pp *PopulationPanel) HasPreset() bool { return pp.preset != "" }

// PresetName returns the preset population modelling this panel ("" if none).
// For aliases such as "LWK" this is the proxy preset ("YRI").
func (pp *PopulationPanel) PresetName() string { return pp.preset }

func presetDemographics(spec PopSpec) map[string]interface{} {
	return map[string]interface{}{
		"model":                           "preset",
		"preset":                          spec.Name,
		"effective_population_size":       spec.Ne,
		"present_population_size":         spec.presentSize(),
		"growth_rate":                     spec.GrowthRate,
		"growth_start_generations":        spec.GrowthStart,
		"bottleneck_start_generations":    spec.BottleneckStart,
		"bottleneck_duration_generations": spec.BottleneckDuration,
		"bottleneck_size":                 spec.BottleneckNe,
		"split_from":                      spec.SplitFrom,
		"split_time_generations":          spec.SplitTime,
		"mutation_rate":                   PresetMutationRate,
	}
}

// LoadFromSimulation simulates regionLength bp for this population and
// replaces the panel's Variants and AlleleFrequencies. The region is split
// into unlinked loci of DefaultLocusLength bp (one genealogy each), and both
// Variants and AlleleFrequencies come from that single simulation. The
// simulation is seeded from the panel's Simulator RNG, so Simulator.SetSeed
// makes it reproducible.
//
// For preset populations the whole preset population tree is used (so the
// panel's lineages pass through its ancestral populations); the panel's
// Simulator supplies this population's branch Ne, the mutation rate and any
// extra events. Other populations use a single constant-size population with
// the Simulator's parameters and events.
func (pp *PopulationPanel) LoadFromSimulation(regionLength int) error {
	pp.mu.Lock()
	defer pp.mu.Unlock()

	if regionLength <= 0 {
		return fmt.Errorf("region length must be positive")
	}
	if pp.SampleSize < 1 {
		return fmt.Errorf("sample size must be at least 1 individual")
	}
	d, label, err := pp.singlePopulationDemography()
	if err != nil {
		return err
	}
	locusLen := DefaultLocusLength
	if regionLength < locusLen {
		locusLen = regionLength
	}
	d.LocusLength = locusLen
	nLoci := (regionLength + locusLen - 1) / locusLen

	jp, err := SimulateJointPanel(d, map[string]int{label: 2 * pp.SampleSize}, nLoci, pp.Simulator.nextSeed())
	if err != nil {
		return fmt.Errorf("failed to simulate population %s: %w", pp.PopulationName, err)
	}
	pp.loadFromJoint(jp, label)
	return nil
}

// singlePopulationDemography returns the demography used to simulate this
// panel on its own and the name of the sampled population in it.
func (pp *PopulationPanel) singlePopulationDemography() (Demography, string, error) {
	ne, mu, events := pp.Simulator.branchParameters()
	if ne < 1 {
		return Demography{}, "", fmt.Errorf("effective population size must be positive")
	}
	if pp.preset != "" {
		d := PresetDemography()
		for i := range d.Pops {
			if d.Pops[i].Name == pp.preset {
				d.Pops[i].Ne = float64(ne)
				d.Pops[i].Events = events
			}
		}
		d.Mu = mu
		return d, pp.preset, nil
	}
	label := pp.PopulationName
	if label == "" {
		label = "sample"
	}
	return Demography{
		Pops: []PopSpec{{Name: label, Ne: float64(ne), Events: events}},
		Mu:   mu,
	}, label, nil
}

// loadFromJoint fills Variants and AlleleFrequencies from population label of
// a joint panel. The caller must hold pp.mu (or own pp exclusively).
func (pp *PopulationPanel) loadFromJoint(jp *JointPanel, label string) {
	i := jp.PopIndex(label)
	n := jp.Haplotypes[i]
	freqs := make(map[string]float64, len(jp.Sites))
	variants := make([]Variant, 0, len(jp.Sites))
	for _, s := range jp.Sites {
		c := s.DerivedCount[i]
		f := float64(c) / float64(n)
		freqs[SiteKey(s.Position)] = f
		if c == 0 || c == n {
			continue
		}
		variants = append(variants, Variant{
			Position:      s.Position,
			Chromosome:    "chr1", // Default chromosome
			ReferenceBase: s.Ancestral,
			AlternateBase: s.Derived,
			Frequency:     f,
			Genotypes:     s.Genotypes[i],
			Quality:       30.0, // Default quality score
		})
	}
	pp.AlleleFrequencies = freqs
	pp.Variants = variants
}

// NewJointPanels simulates the named populations together under the preset
// demography (one genealogy per locus spanning all of them) and returns one
// panel per name, with sampleSize diploid individuals each, plus the joint
// sample. All panels share the same site keys, which is what cross-population
// statistics such as F_ST need. Names are resolved with ResolvePreset; unknown
// names, or two names resolving to the same preset, are errors.
func NewJointPanels(populations []string, sampleSize, nLoci int, seed uint64) (map[string]*PopulationPanel, *JointPanel, error) {
	if sampleSize < 1 {
		return nil, nil, fmt.Errorf("sample size must be at least 1 individual")
	}
	samples := make(map[string]int, len(populations))
	resolved := make(map[string]string, len(populations))
	for _, name := range populations {
		preset, ok := ResolvePreset(name)
		if !ok {
			return nil, nil, fmt.Errorf("no demographic preset for population %q (modelled: %v)", name, PresetNames())
		}
		if _, dup := samples[preset]; dup {
			return nil, nil, fmt.Errorf("population %q resolves to preset %s more than once", name, preset)
		}
		samples[preset] = 2 * sampleSize
		resolved[name] = preset
	}
	jp, err := SimulateJointPanel(PresetDemography(), samples, nLoci, seed)
	if err != nil {
		return nil, nil, err
	}
	panels := make(map[string]*PopulationPanel, len(populations))
	for _, name := range populations {
		pp := NewPopulationPanel(name, sampleSize)
		pp.loadFromJoint(jp, resolved[name])
		pp.Demographics["joint_populations"] = append([]string(nil), jp.Populations...)
		panels[name] = pp
	}
	return panels, jp, nil
}

// GetVariants returns all variants in the panel
func (pp *PopulationPanel) GetVariants() []Variant {
	pp.mu.RLock()
	defer pp.mu.RUnlock()
	return pp.Variants
}

// GetAlleleFrequencies returns allele frequencies normalized as a probability distribution
func (pp *PopulationPanel) GetAlleleFrequencies() map[string]float64 {
	pp.mu.RLock()
	defer pp.mu.RUnlock()

	// If we have cached normalized frequencies, return them
	if len(pp.AlleleFrequencies) > 0 {
		return pp.AlleleFrequencies
	}

	// Generate normalized frequencies from variants
	if len(pp.Variants) > 0 {
		frequencies := make(map[string]float64)
		totalWeight := 0.0

		for i, variant := range pp.Variants {
			siteID := fmt.Sprintf("site_%d", i)
			// Use minor allele frequency as the weight, but ensure it's positive
			weight := math.Min(variant.Frequency, 1.0-variant.Frequency)
			if weight <= 0 {
				weight = 0.001 // Small positive value for monomorphic sites
			}
			frequencies[siteID] = weight
			totalWeight += weight
		}

		// Normalize to make a proper probability distribution (sum = 1.0)
		if totalWeight > 0 {
			for siteID := range frequencies {
				frequencies[siteID] /= totalWeight
			}
		}

		return frequencies
	}

	return make(map[string]float64)
}

// GetDemographics returns demographic information
func (pp *PopulationPanel) GetDemographics() map[string]interface{} {
	pp.mu.RLock()
	defer pp.mu.RUnlock()
	return pp.Demographics
}

// SetDemographicEvent adds a demographic event to the population
func (pp *PopulationPanel) SetDemographicEvent(event DemographicEvent) {
	pp.mu.Lock()
	defer pp.mu.Unlock()
	pp.Simulator.AddDemographicEvent(event)
}

// GetPopulationStats returns population statistics
func (pp *PopulationPanel) GetPopulationStats() map[string]interface{} {
	pp.mu.RLock()
	defer pp.mu.RUnlock()

	// Calculate heterozygosity from variants
	heterozygosity := 0.0
	if len(pp.Variants) > 0 {
		totalHeterozygosity := 0.0
		for _, variant := range pp.Variants {
			// Heterozygosity = 2 * p * (1-p)
			p := variant.Frequency
			totalHeterozygosity += 2 * p * (1 - p)
		}
		heterozygosity = totalHeterozygosity / float64(len(pp.Variants))
	}

	stats := map[string]interface{}{
		"population_name": pp.PopulationName,
		"sample_size":     pp.SampleSize,
		"variants":        len(pp.Variants),
		"allele_sites":    len(pp.AlleleFrequencies),
		"heterozygosity":  heterozygosity,
	}

	// Add simulator stats
	simStats := pp.Simulator.GetSimulationStats()
	for k, v := range simStats {
		stats["sim_"+k] = v
	}

	return stats
}
