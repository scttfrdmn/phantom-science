package population

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sort"
	"sync"
	"time"
)

// CoalescentSimulator implements coalescent theory for population genetics simulation
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
	Type          string                 // "bottleneck", "expansion", "migration", "admixture"
	Time          int                    // Generations ago
	Parameters    map[string]float64     // Event-specific parameters
	Populations   []string               // Affected populations
	Duration      int                    // Duration in generations
}

// Genealogy represents a genealogical tree for a genomic region
type Genealogy struct {
	Nodes         []*CoalescentNode
	Root          *CoalescentNode
	LeafNodes     []*CoalescentNode
	TotalTime     float64
	MutationEvents []MutationEvent
	RecombinationEvents []RecombinationEvent
}

// CoalescentNode represents a node in the coalescent tree
type CoalescentNode struct {
	ID            int
	Time          float64
	Parent        *CoalescentNode
	Children      []*CoalescentNode
	Population    string
	Mutations     []Mutation
	IsLeaf        bool
	LineageCount  int
}

// MutationEvent represents a mutation event in the genealogy
type MutationEvent struct {
	Position    int
	Time        float64
	Lineage     *CoalescentNode
	Mutation    Mutation
	Type        string // "transition", "transversion", "indel"
}

// RecombinationEvent represents a recombination event
type RecombinationEvent struct {
	Position    int
	Time        float64
	Lineages    []*CoalescentNode
}

// Mutation represents a single mutation
type Mutation struct {
	Position    int
	Ancestral   byte
	Derived     byte
	Time        float64
	Frequency   float64
	IsFixed     bool
}

// AlleleFrequencySpectrum represents the site frequency spectrum
type AlleleFrequencySpectrum struct {
	Frequencies   []float64  // Frequency bins
	Counts        []int      // Count of sites in each bin
	TotalSites    int
	ThetaWatterson float64
	ThetaPi       float64
	TajimaD       float64
	FayWuH        float64
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

	// Simulate coalescent process
	for len(activeLineages) > 1 {
		// Calculate coalescence rate
		n := len(activeLineages)
		coalescenceRate := float64(n*(n-1)) / (2.0 * float64(s.effectivePopSize))
		
		// Apply demographic events
		adjustedPopSize := s.getAdjustedPopSize(currentTime)
		coalescenceRate = float64(n*(n-1)) / (2.0 * float64(adjustedPopSize))

		// Time to next coalescence event
		timeToCoalescence := s.rng.ExpFloat64() / coalescenceRate
		currentTime += timeToCoalescence

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

// getAdjustedPopSize returns the effective population size adjusted for demographic events
func (s *CoalescentSimulator) getAdjustedPopSize(currentTime float64) int {
	adjustedSize := s.effectivePopSize

	for _, event := range s.demographicHistory {
		eventTime := float64(event.Time)
		
		if currentTime >= eventTime {
			switch event.Type {
			case "bottleneck":
				if bottleneckSize, ok := event.Parameters["population_size"]; ok {
					adjustedSize = int(bottleneckSize)
				}
			case "expansion":
				if expansionFactor, ok := event.Parameters["expansion_factor"]; ok {
					adjustedSize = int(float64(adjustedSize) * expansionFactor)
				}
			}
		}
	}

	return adjustedSize
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

// addRecombinationEvents adds recombination events to the genealogy
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

// SimulateAlleleFrequencies simulates allele frequencies for a number of sites
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
		binIndex := int(frequency * float64(sampleSize)) - 1
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

// calculateThetaPi calculates nucleotide diversity (pi)
func (s *CoalescentSimulator) calculateThetaPi(spectrum *AlleleFrequencySpectrum, sampleSize int) float64 {
	pi := 0.0
	
	for i, count := range spectrum.Counts {
		frequency := spectrum.Frequencies[i]
		pi += frequency * (1 - frequency) * float64(count)
	}
	
	return pi * float64(sampleSize) / float64(sampleSize-1)
}

// calculateTajimaD calculates Tajima's D statistic
func (s *CoalescentSimulator) calculateTajimaD(thetaPi, thetaWatterson float64, sampleSize int) float64 {
	if thetaWatterson == 0 {
		return 0.0
	}
	
	// Simplified calculation - in practice, this would include variance calculations
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
	
	variance := e1*thetaWatterson + e2*thetaWatterson*thetaWatterson
	
	if variance <= 0 {
		return 0.0 // undefined when variance is non-positive; guards math.Sqrt(negative) → NaN
	}

	return (thetaPi - thetaWatterson) / math.Sqrt(variance)
}

// calculateFayWuH calculates Fay and Wu's H statistic
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
		"mutation_rate":            s.mutationRate,
		"recombination_rate":       s.recombinationRate,
		"demographic_events":       len(s.demographicHistory),
		"genealogies_simulated":    len(s.genealogies),
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

// GetDemographicEvents returns the demographic events in the simulation
func (s *CoalescentSimulator) GetDemographicEvents() []DemographicEvent {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.demographicHistory
}

// PopulationPanel represents a panel of population genetic data
type PopulationPanel struct {
	PopulationName    string
	SampleSize        int
	Simulator         *CoalescentSimulator
	AlleleFrequencies map[string]float64
	Demographics      map[string]interface{}
	Variants          []Variant
	mu                sync.RWMutex
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

// NewPopulationPanel creates a new population panel
func NewPopulationPanel(populationName string, sampleSize int) *PopulationPanel {
	simulator := NewCoalescentSimulator(10000, 1e-8, 1e-8) // Default parameters
	
	return &PopulationPanel{
		PopulationName:    populationName,
		SampleSize:        sampleSize,
		Simulator:         simulator,
		AlleleFrequencies: make(map[string]float64),
		Demographics:      make(map[string]interface{}),
		Variants:          make([]Variant, 0),
	}
}

// LoadFromSimulation loads panel data from coalescent simulation
func (pp *PopulationPanel) LoadFromSimulation(regionLength int) error {
	pp.mu.Lock()
	defer pp.mu.Unlock()
	
	// Simulate genealogy for the population
	genealogy, err := pp.Simulator.SimulateGenealogy(pp.PopulationName, pp.SampleSize, regionLength)
	if err != nil {
		return fmt.Errorf("failed to simulate genealogy: %w", err)
	}
	
	// Extract variants from genealogy
	for _, mutationEvent := range genealogy.MutationEvents {
		variant := Variant{
			Position:      mutationEvent.Position,
			Chromosome:    "chr1", // Default chromosome
			ReferenceBase: mutationEvent.Mutation.Ancestral,
			AlternateBase: mutationEvent.Mutation.Derived,
			Frequency:     mutationEvent.Mutation.Frequency,
			Quality:       30.0, // Default quality score
		}
		
		// Generate genotypes based on frequency
		variant.Genotypes = pp.generateGenotypes(variant.Frequency)
		pp.Variants = append(pp.Variants, variant)
	}
	
	// Store allele frequencies
	frequencies, err := pp.Simulator.SimulateAlleleFrequencies(regionLength, pp.SampleSize)
	if err != nil {
		return fmt.Errorf("failed to simulate allele frequencies: %w", err)
	}
	
	pp.AlleleFrequencies = frequencies
	
	return nil
}

// generateGenotypes generates genotypes based on allele frequency
func (pp *PopulationPanel) generateGenotypes(frequency float64) []int {
	genotypes := make([]int, pp.SampleSize)
	
	// Use Hardy-Weinberg equilibrium
	p := frequency
	q := 1.0 - p
	
	for i := 0; i < pp.SampleSize; i++ {
		r := rand.Float64() // #nosec G404
		if r < q*q {
			genotypes[i] = 0 // homozygous reference
		} else if r < q*q + 2*p*q {
			genotypes[i] = 1 // heterozygous
		} else {
			genotypes[i] = 2 // homozygous alternate
		}
	}
	
	return genotypes
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