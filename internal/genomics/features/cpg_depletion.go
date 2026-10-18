package features

import (
	"fmt"
	"math"
	"math/rand/v2"
	"strings"
	"time"
)

// CpGDepletionModel applies a heuristic CpG→TpG depletion to sequences
// outside detected CpG islands and scores CpG observed/expected ratios.
// Islands use the Gardiner-Garden & Frommer (1987) criteria (200 bp window,
// GC ≥ 0.5, CpG o/e ≥ 0.6). The depletion probabilities are fixed parameters,
// not a population-genetic mutation model.
type CpGDepletionModel struct {
	depletionRate      float64           // Rate of CpG -> TpG transitions
	islandLocations    []GenomicInterval // Stored by LoadCpGIslands; not used by ApplyDepletion
	methylationRate    float64           // Methylation rate in non-island regions
	transitionProb     float64           // Not used by any method (reported only)
	gcContentThreshold float64           // GC content threshold for CpG islands
	lengthThreshold    int               // Not used; the island window is fixed at 200 bp
	rng                *rand.Rand        // Random number generator
}

// GenomicInterval represents a genomic interval [start, end)
type GenomicInterval struct {
	Start      int
	End        int
	Chromosome string
	Strand     string
	GCContent  float64
	CpGRatio   float64
}

// CpGIsland represents a CpG island with its properties
type CpGIsland struct {
	Interval    GenomicInterval
	ObservedCpG int
	ExpectedCpG float64
	CpGRatio    float64
	GCContent   float64
	Length      int
	IsProtected bool // Whether protected from methylation
}

// DepletionResult represents the result of applying CpG depletion
type DepletionResult struct {
	OriginalSequence  []byte
	DepletedSequence  []byte
	CpGSitesOriginal  int
	CpGSitesRemaining int
	DepletionRate     float64
	IslandRegions     []CpGIsland
	TransitionEvents  []TransitionEvent
}

// TransitionEvent represents a C->T transition event
type TransitionEvent struct {
	Position         int
	OriginalBase     byte
	NewBase          byte
	InCpGIsland      bool
	MethylationState string // "methylated", "unmethylated", "hemimethylated"
	TransitionType   string // "deamination", "replication_error"
}

// NewCpGDepletionModel creates a CpG depletion model. Its random source is
// seeded from the clock; call SetSeed for reproducible output.
func NewCpGDepletionModel() *CpGDepletionModel {
	return &CpGDepletionModel{
		depletionRate:      0.8,                                                     // per-call probability factor for a methylated non-island CpG
		methylationRate:    0.7,                                                     // fraction of non-island CpG sites treated as methylated
		transitionProb:     0.02,                                                    // unused
		gcContentThreshold: 0.5,                                                     // 50% GC content threshold
		lengthThreshold:    200,                                                     // 200bp minimum length
		rng:                rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
	}
}

// SetDepletionRate sets the CpG depletion rate
func (m *CpGDepletionModel) SetDepletionRate(rate float64) {
	if rate >= 0 && rate <= 1 {
		m.depletionRate = rate
	}
}

// SetMethylationRate sets the methylation rate
func (m *CpGDepletionModel) SetMethylationRate(rate float64) {
	if rate >= 0 && rate <= 1 {
		m.methylationRate = rate
	}
}

// LoadCpGIslands stores known CpG island locations. ApplyDepletion ignores
// them and detects islands from the sequence instead.
func (m *CpGDepletionModel) LoadCpGIslands(islands []CpGIsland) {
	m.islandLocations = make([]GenomicInterval, len(islands))
	for i, island := range islands {
		m.islandLocations[i] = island.Interval
	}
}

// ApplyDepletion turns CpG sites into TpG: outside detected islands with
// probability methylationRate × depletionRate × a context factor, inside
// islands with probability 0.05. Only upper-case CG is considered.
func (m *CpGDepletionModel) ApplyDepletion(sequence []byte) (*DepletionResult, error) {
	if len(sequence) == 0 {
		return nil, fmt.Errorf("sequence cannot be empty")
	}

	result := &DepletionResult{
		OriginalSequence: make([]byte, len(sequence)),
		DepletedSequence: make([]byte, len(sequence)),
		TransitionEvents: make([]TransitionEvent, 0),
	}

	// Copy original sequence
	copy(result.OriginalSequence, sequence)
	copy(result.DepletedSequence, sequence)

	// Identify CpG islands in the sequence
	islands := m.identifyCpGIslands(sequence)
	result.IslandRegions = islands

	// Count original CpG sites
	result.CpGSitesOriginal = m.countCpGSites(sequence)

	// Apply depletion to non-island regions
	for i := 0; i < len(sequence)-1; i++ {
		if sequence[i] == 'C' && sequence[i+1] == 'G' {
			// Check if this CpG site is in an island
			inIsland := m.isInCpGIsland(i, islands)

			if !inIsland {
				// Apply depletion based on methylation rate and sequence context
				if m.shouldDeplete(i, sequence) {
					// C -> T transition (deamination of methylated cytosine)
					result.DepletedSequence[i] = 'T'

					// Record transition event
					event := TransitionEvent{
						Position:         i,
						OriginalBase:     'C',
						NewBase:          'T',
						InCpGIsland:      false,
						MethylationState: "methylated",
						TransitionType:   "deamination",
					}
					result.TransitionEvents = append(result.TransitionEvents, event)
				}
			} else {
				// CpG islands are protected from methylation
				// But can still have some depletion at a much lower rate
				if m.rng.Float64() < 0.05 { // fixed 5% per call inside islands
					result.DepletedSequence[i] = 'T'

					event := TransitionEvent{
						Position:         i,
						OriginalBase:     'C',
						NewBase:          'T',
						InCpGIsland:      true,
						MethylationState: "unmethylated",
						TransitionType:   "replication_error",
					}
					result.TransitionEvents = append(result.TransitionEvents, event)
				}
			}
		}
	}

	// Count remaining CpG sites
	result.CpGSitesRemaining = m.countCpGSites(result.DepletedSequence)

	// Calculate actual depletion rate
	if result.CpGSitesOriginal > 0 {
		result.DepletionRate = float64(result.CpGSitesOriginal-result.CpGSitesRemaining) / float64(result.CpGSitesOriginal)
	}

	return result, nil
}

// identifyCpGIslands identifies CpG islands in a sequence
func (m *CpGDepletionModel) identifyCpGIslands(sequence []byte) []CpGIsland {
	islands := make([]CpGIsland, 0)
	windowSize := 200
	step := 50

	for i := 0; i <= len(sequence)-windowSize; i += step {
		window := sequence[i : i+windowSize]

		// Calculate GC content
		gcContent := m.calculateGCContent(window)

		// Calculate CpG ratio
		cpgRatio := m.calculateCpGRatio(window)

		// Check if this window qualifies as a CpG island
		if gcContent >= m.gcContentThreshold && cpgRatio >= 0.6 {
			// Extend the island
			island := m.extendCpGIsland(sequence, i, windowSize)

			// Check if this island overlaps with existing islands
			if !m.overlapsWithExistingIsland(island, islands) {
				islands = append(islands, island)
			}
		}
	}

	return islands
}

// extendCpGIsland extends a CpG island to its natural boundaries
func (m *CpGDepletionModel) extendCpGIsland(sequence []byte, start, initialLength int) CpGIsland {
	// Extend to the left
	leftBoundary := start
	for leftBoundary > 0 {
		testStart := leftBoundary - 50
		if testStart < 0 {
			testStart = 0
		}

		testWindow := sequence[testStart : leftBoundary+initialLength]
		if m.calculateGCContent(testWindow) < m.gcContentThreshold*0.8 {
			break
		}
		leftBoundary = testStart
	}

	// Extend to the right
	rightBoundary := start + initialLength
	for rightBoundary < len(sequence) {
		testEnd := rightBoundary + 50
		if testEnd > len(sequence) {
			testEnd = len(sequence)
		}

		testWindow := sequence[leftBoundary:testEnd]
		if m.calculateGCContent(testWindow) < m.gcContentThreshold*0.8 {
			break
		}
		rightBoundary = testEnd
	}

	// Create island
	islandSequence := sequence[leftBoundary:rightBoundary]
	island := CpGIsland{
		Interval: GenomicInterval{
			Start: leftBoundary,
			End:   rightBoundary,
		},
		ObservedCpG: m.countCpGSites(islandSequence),
		GCContent:   m.calculateGCContent(islandSequence),
		Length:      rightBoundary - leftBoundary,
		IsProtected: true,
	}

	// Calculate expected CpG sites using separate G and C counts (consistent
	// with calculateCpGRatio): E[CpG] ≈ P(G)×P(C)×N = gCount×cCount/N
	gCount, cCount := 0, 0
	for _, base := range islandSequence {
		if base == 'G' {
			gCount++
		} else if base == 'C' {
			cCount++
		}
	}
	if len(islandSequence) > 0 {
		island.ExpectedCpG = float64(gCount) * float64(cCount) / float64(len(islandSequence))
	}

	// Calculate CpG ratio
	if island.ExpectedCpG > 0 {
		island.CpGRatio = float64(island.ObservedCpG) / island.ExpectedCpG
	}

	return island
}

// overlapsWithExistingIsland checks if a new island overlaps with existing ones
func (m *CpGDepletionModel) overlapsWithExistingIsland(newIsland CpGIsland, existingIslands []CpGIsland) bool {
	for _, existing := range existingIslands {
		// Check for overlap
		if newIsland.Interval.Start < existing.Interval.End && newIsland.Interval.End > existing.Interval.Start {
			return true
		}
	}
	return false
}

// isInCpGIsland checks if a position is within a CpG island
func (m *CpGDepletionModel) isInCpGIsland(position int, islands []CpGIsland) bool {
	for _, island := range islands {
		if position >= island.Interval.Start && position < island.Interval.End {
			return true
		}
	}
	return false
}

// shouldDeplete determines if a CpG site should be depleted
func (m *CpGDepletionModel) shouldDeplete(position int, sequence []byte) bool {
	// Check if the cytosine is methylated
	if m.rng.Float64() > m.methylationRate {
		return false // Not methylated, no depletion
	}

	// Check local sequence context for methylation probability
	contextMethylation := m.getContextMethylationProbability(position, sequence)

	// Apply depletion based on context and evolutionary time
	depleteProb := m.depletionRate * contextMethylation

	return m.rng.Float64() < depleteProb
}

// getContextMethylationProbability returns a context factor for a CpG site.
// The factor is capped at 1.0, so the CCGG and GCGC boosts have no effect;
// only CGCG (0.8) changes the result.
func (m *CpGDepletionModel) getContextMethylationProbability(position int, sequence []byte) float64 {
	baseProbability := 1.0

	// Check for specific sequence contexts that affect methylation
	if position >= 2 && position < len(sequence)-2 {
		context := string(sequence[position-2 : position+3])

		// Some contexts are more likely to be methylated
		switch {
		case strings.Contains(context, "CCGG"):
			baseProbability *= 1.2 // Higher methylation in CCGG context
		case strings.Contains(context, "GCGC"):
			baseProbability *= 1.1 // Higher methylation in GCGC context
		case strings.Contains(context, "CGCG"):
			baseProbability *= 0.8 // Lower methylation in CGCG repeats
		}
	}

	// Cap at 1.0
	if baseProbability > 1.0 {
		baseProbability = 1.0
	}

	return baseProbability
}

// calculateGCContent calculates GC content of a sequence
func (m *CpGDepletionModel) calculateGCContent(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0.0
	}

	gcCount := 0
	totalCount := 0

	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
		if base == 'A' || base == 'T' || base == 'G' || base == 'C' {
			totalCount++
		}
	}

	if totalCount == 0 {
		return 0.0
	}

	return float64(gcCount) / float64(totalCount)
}

// calculateCpGRatio calculates the CpG observed/expected ratio
func (m *CpGDepletionModel) calculateCpGRatio(sequence []byte) float64 {
	if len(sequence) < 2 {
		return 0.0
	}

	// Count CpG dinucleotides
	cpgCount := m.countCpGSites(sequence)

	// Count G and C bases
	gCount := 0
	cCount := 0
	for _, base := range sequence {
		if base == 'G' {
			gCount++
		} else if base == 'C' {
			cCount++
		}
	}

	// Calculate expected CpG frequency
	expectedCpG := float64(gCount) * float64(cCount) / float64(len(sequence))

	if expectedCpG == 0 {
		return 0.0
	}

	return float64(cpgCount) / expectedCpG
}

// countCpGSites counts CpG dinucleotides in a sequence
func (m *CpGDepletionModel) countCpGSites(sequence []byte) int {
	count := 0
	for i := 0; i < len(sequence)-1; i++ {
		if sequence[i] == 'C' && sequence[i+1] == 'G' {
			count++
		}
	}
	return count
}

// ValidateDepletion validates the depletion result
func (m *CpGDepletionModel) ValidateDepletion(result *DepletionResult) error {
	if result == nil {
		return fmt.Errorf("depletion result is nil")
	}

	if len(result.OriginalSequence) != len(result.DepletedSequence) {
		return fmt.Errorf("sequence lengths do not match")
	}

	if result.CpGSitesOriginal < result.CpGSitesRemaining {
		return fmt.Errorf("depleted sequence has more CpG sites than original")
	}

	// Check that only C->T transitions occurred
	for i, event := range result.TransitionEvents {
		if event.OriginalBase != 'C' || event.NewBase != 'T' {
			return fmt.Errorf("invalid transition event at index %d: %c->%c", i, event.OriginalBase, event.NewBase)
		}
	}

	return nil
}

// GetDepletionStats returns statistics about the depletion process
func (m *CpGDepletionModel) GetDepletionStats(result *DepletionResult) map[string]interface{} {
	stats := map[string]interface{}{
		"original_cpg_sites":  result.CpGSitesOriginal,
		"remaining_cpg_sites": result.CpGSitesRemaining,
		"depletion_rate":      result.DepletionRate,
		"transition_events":   len(result.TransitionEvents),
		"cpg_islands_found":   len(result.IslandRegions),
		"sequence_length":     len(result.OriginalSequence),
	}

	// Calculate island statistics
	if len(result.IslandRegions) > 0 {
		totalIslandLength := 0
		totalIslandCpG := 0
		avgGCContent := 0.0

		for _, island := range result.IslandRegions {
			totalIslandLength += island.Length
			totalIslandCpG += island.ObservedCpG
			avgGCContent += island.GCContent
		}

		stats["total_island_length"] = totalIslandLength
		stats["total_island_cpg"] = totalIslandCpG
		stats["avg_island_gc_content"] = avgGCContent / float64(len(result.IslandRegions))
		stats["island_coverage"] = float64(totalIslandLength) / float64(len(result.OriginalSequence))
	}

	// Calculate transition statistics
	transitionsInIslands := 0
	transitionsOutsideIslands := 0

	for _, event := range result.TransitionEvents {
		if event.InCpGIsland {
			transitionsInIslands++
		} else {
			transitionsOutsideIslands++
		}
	}

	stats["transitions_in_islands"] = transitionsInIslands
	stats["transitions_outside_islands"] = transitionsOutsideIslands

	return stats
}

// SimulateEvolutionaryDepletion calls ApplyDepletion generations times,
// feeding each output into the next call. It has no time scale or mutation
// rate; each call is one round of the same fixed probabilities.
func (m *CpGDepletionModel) SimulateEvolutionaryDepletion(sequence []byte, generations int) (*DepletionResult, error) {
	currentSequence := make([]byte, len(sequence))
	copy(currentSequence, sequence)

	allTransitionEvents := make([]TransitionEvent, 0)

	// Apply depletion over multiple generations
	for gen := 0; gen < generations; gen++ {
		result, err := m.ApplyDepletion(currentSequence)
		if err != nil {
			return nil, err
		}

		// Update sequence for next generation
		currentSequence = result.DepletedSequence

		// Accumulate transition events
		allTransitionEvents = append(allTransitionEvents, result.TransitionEvents...)
	}

	// Create final result
	finalResult := &DepletionResult{
		OriginalSequence:  sequence,
		DepletedSequence:  currentSequence,
		CpGSitesOriginal:  m.countCpGSites(sequence),
		CpGSitesRemaining: m.countCpGSites(currentSequence),
		TransitionEvents:  allTransitionEvents,
		IslandRegions:     m.identifyCpGIslands(currentSequence),
	}

	if finalResult.CpGSitesOriginal > 0 {
		finalResult.DepletionRate = float64(finalResult.CpGSitesOriginal-finalResult.CpGSitesRemaining) / float64(finalResult.CpGSitesOriginal)
	}

	return finalResult, nil
}

// SetSeed sets the random seed for reproducible depletion
func (m *CpGDepletionModel) SetSeed(seed int64) {
	m.rng = rand.New(rand.NewPCG(uint64(seed), 0)) // #nosec G115 G404
}

// GetModelParameters returns the current model parameters
func (m *CpGDepletionModel) GetModelParameters() map[string]interface{} {
	return map[string]interface{}{
		"depletion_rate":       m.depletionRate,
		"methylation_rate":     m.methylationRate,
		"transition_prob":      m.transitionProb,
		"gc_content_threshold": m.gcContentThreshold,
		"length_threshold":     m.lengthThreshold,
		"known_islands":        len(m.islandLocations),
	}
}

// CpGDepletionDetails contains detailed CpG depletion analysis results
type CpGDepletionDetails struct {
	CpGRatio          float64 // observed / expected CpG; > 1.0 means CpG-enriched
	DepletionScore    float64 // clamped to [0, 1]: 1.0 = fully depleted, 0.0 = at or above expected
	RawDepletionScore float64 // unclamped 1.0 - CpGRatio; negative for CpG-enriched sequences
}

// AnalyzeSequence analyzes CpG depletion in a DNA sequence
func (m *CpGDepletionModel) AnalyzeSequence(seq []byte) (float64, *CpGDepletionDetails, error) {
	if len(seq) == 0 {
		return 0.0, &CpGDepletionDetails{}, nil
	}
	cpgCount := m.countCpGSites(seq)
	gCount := 0
	cCount := 0
	for _, b := range seq {
		if b == 'G' {
			gCount++
		} else if b == 'C' {
			cCount++
		}
	}
	expectedCpG := float64(gCount) * float64(cCount) / float64(len(seq))
	cpgRatio := 0.0
	if expectedCpG > 0 {
		cpgRatio = float64(cpgCount) / expectedCpG
	}
	// Score: 1.0 = fully depleted (ratio ≈ 0), 0.0 = at or above expected (ratio ≥ 1).
	// Clamped to [0, 1] — callers needing to distinguish CpG-enriched sequences from
	// merely non-depleted ones should inspect RawDepletionScore, which is negative when
	// cpgRatio > 1.0.
	rawScore := 1.0 - cpgRatio
	depletionScore := math.Max(0, rawScore)
	return depletionScore, &CpGDepletionDetails{
		CpGRatio:          cpgRatio,
		DepletionScore:    depletionScore,
		RawDepletionScore: rawScore,
	}, nil
}
