package authenticity

import (
	"fmt"
	"math/rand/v2"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/errors"
)

// AdaptiveHMMModel represents an adaptive HMM model that dynamically adjusts for sequence complexity
type AdaptiveHMMModel struct {
	k                      int
	stateSpace             []string
	transitionMatrix       [][]float64
	emissionMatrix         [][]float64
	startProbabilities     []float64
	
	// Adaptive complexity features
	complexityMemory       []string              // Recent k-mer history for complexity tracking
	complexityWindowSize   int                   // Size of complexity tracking window
	targetComplexity       float64               // Target sequence complexity (0.6-0.8)
	currentComplexity      float64               // Current measured complexity
	adaptiveWeights        map[string]float64    // Dynamic weights based on complexity
	
	// Biological constraints
	gcContentTarget        float64
	gcContentTolerance     float64
	currentGCContent       float64
	homopolymerMaxLength   int
	cpgDepletionRatio      float64
	
	// Training and generation state
	totalKmers             uint64
	trained                bool
	mu                     sync.RWMutex
	rng                    *rand.Rand
	
	// Performance metrics
	generationStats        GenerationStats
}

// GenerationStats tracks generation performance and quality
type GenerationStats struct {
	TotalSequences      int
	TotalBases          int64
	AverageComplexity   float64
	AverageGCContent    float64
	MaxHomopolymer      int
	AverageJSDivergence float64
}

// NewAdaptiveHMMModel creates a new adaptive HMM model
func NewAdaptiveHMMModel(k int) *AdaptiveHMMModel {
	if k < 2 || k > 8 {
		// Use proper error handling instead of panic
		err := errors.NewValidationError("adaptive_hmm_model", "new_model", fmt.Sprintf("k must be between 2 and 8, got %d", k))
		errors.HandleError(err)
		return nil
	}

	model := &AdaptiveHMMModel{
		k:                    k,
		rng:                  rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
		complexityWindowSize: 100,              // Track last 100 k-mers
		targetComplexity:     0.7,              // Target 70% complexity
		gcContentTarget:      0.42,             // Target 42% GC content
		gcContentTolerance:   0.015,            // ±1.5% tolerance (enhanced)
		homopolymerMaxLength: 8,                // Max 8-base homopolymer runs
		cpgDepletionRatio:    0.25,             // CpG depletion ratio
		adaptiveWeights:      make(map[string]float64),
		complexityMemory:     make([]string, 0),
	}

	model.initializeStateSpace()
	model.initializeAdaptiveWeights()
	return model
}

// initializeStateSpace generates all possible (k-1)-mers as states
func (m *AdaptiveHMMModel) initializeStateSpace() {
	bases := []byte{'A', 'T', 'G', 'C'}
	stateLength := m.k - 1

	var generateStates func(current string, remaining int)
	generateStates = func(current string, remaining int) {
		if remaining == 0 {
			m.stateSpace = append(m.stateSpace, current)
			return
		}
		for _, base := range bases {
			generateStates(current+string(base), remaining-1)
		}
	}

	generateStates("", stateLength)

	// Initialize matrices
	numStates := len(m.stateSpace)
	m.transitionMatrix = make([][]float64, numStates)
	m.emissionMatrix = make([][]float64, numStates)
	m.startProbabilities = make([]float64, numStates)

	for i := range m.transitionMatrix {
		m.transitionMatrix[i] = make([]float64, 4) // 4 bases
		m.emissionMatrix[i] = make([]float64, 42)  // Quality scores 0-41
	}
}

// initializeAdaptiveWeights initializes weights for adaptive complexity control
func (m *AdaptiveHMMModel) initializeAdaptiveWeights() {
	// Initialize with uniform weights
	for _, state := range m.stateSpace {
		m.adaptiveWeights[state] = 1.0
	}
}

// TrainFromReferenceGenome trains the adaptive model with complexity awareness
func (m *AdaptiveHMMModel) TrainFromReferenceGenome(refGenome GenomeInterface) error {
	m.mu.Lock()
	defer m.mu.Unlock()

	// Initialize counters
	transitionCounts := make([][]int, len(m.stateSpace))
	emissionCounts := make([][]int, len(m.stateSpace))
	startCounts := make([]int, len(m.stateSpace))

	for i := range transitionCounts {
		transitionCounts[i] = make([]int, 4)
		emissionCounts[i] = make([]int, 42)
	}

	// Process each chromosome with complexity tracking
	chromosomeNames := refGenome.GetChromosomeNames()
	for _, chrName := range chromosomeNames {
		chr, err := refGenome.LoadChromosome(chrName)
		if err != nil {
			continue
		}

		err = m.trainOnChromosomeAdaptive(chr, transitionCounts, emissionCounts, startCounts)
		if err != nil {
			return errors.WrapError(err, errors.BiologicalModelError, "adaptive_hmm_model", "train_from_reference", fmt.Sprintf("failed to train on chromosome %s", chrName))
		}
	}

	// Normalize with adaptive complexity control
	m.normalizeWithComplexityControl(transitionCounts, emissionCounts, startCounts)

	m.trained = true
	return nil
}

// trainOnChromosomeAdaptive trains on a single chromosome with adaptive complexity tracking
func (m *AdaptiveHMMModel) trainOnChromosomeAdaptive(chr ChromosomeInterface, 
	transitionCounts, emissionCounts [][]int, startCounts []int) error {
	
	sequence := chr.GetSequence()
	if len(sequence) < m.k {
		return nil
	}

	// Track complexity during training
	var kmerWindow []string
	
	// Process sequence with adaptive complexity awareness
	for i := 0; i <= len(sequence)-m.k; i++ {
		kmer := string(sequence[i : i+m.k])
		
		// Skip k-mers with non-standard bases
		if !m.isValidKmer(kmer) {
			continue
		}

		// Update complexity window
		kmerWindow = append(kmerWindow, kmer)
		if len(kmerWindow) > m.complexityWindowSize {
			kmerWindow = kmerWindow[1:]
		}

		// Calculate current window complexity
		windowComplexity := m.calculateWindowComplexity(kmerWindow)
		
		// Apply adaptive weighting based on complexity
		weight := m.getAdaptiveWeight(kmer, windowComplexity)
		
		// Apply biological filters
		if !m.passesBiologicalFilters(kmer, kmerWindow) {
			continue
		}

		state := kmer[:m.k-1]
		nextBase := kmer[m.k-1]

		stateIndex := m.getStateIndex(state)
		if stateIndex == -1 {
			continue
		}

		baseIndex := m.getBaseIndex(nextBase)
		if baseIndex == -1 {
			continue
		}

		// Count transitions with adaptive weighting
		if i == 0 {
			startCounts[stateIndex] += int(weight * 100)
		} else {
			transitionCounts[stateIndex][baseIndex] += int(weight * 100)
		}

		// Generate quality score
		qualityScore := m.generateAdaptiveQuality(kmer, windowComplexity)
		if qualityScore >= 0 && qualityScore < 42 {
			emissionCounts[stateIndex][qualityScore] += int(weight * 100)
		}

		m.totalKmers++
	}

	return nil
}

// calculateWindowComplexity calculates the complexity of a k-mer window
func (m *AdaptiveHMMModel) calculateWindowComplexity(window []string) float64 {
	if len(window) < 2 {
		return 0.0
	}

	// Count unique k-mers in window
	uniqueKmers := make(map[string]bool)
	for _, kmer := range window {
		uniqueKmers[kmer] = true
	}

	// Calculate complexity as ratio of unique k-mers to total k-mers
	return float64(len(uniqueKmers)) / float64(len(window))
}

// getAdaptiveWeight calculates adaptive weight based on complexity needs
func (m *AdaptiveHMMModel) getAdaptiveWeight(kmer string, currentComplexity float64) float64 {
	baseWeight := 1.0
	
	// If current complexity is too low, favor rare k-mers
	if currentComplexity < m.targetComplexity {
		// Check if this k-mer would increase complexity
		if m.wouldIncreaseComplexity(kmer) {
			baseWeight *= 2.0 // Favor complexity-increasing k-mers
		} else {
			baseWeight *= 0.5 // Reduce weight for common k-mers
		}
	}
	
	// If current complexity is too high, favor common k-mers
	if currentComplexity > m.targetComplexity + 0.1 {
		if m.wouldIncreaseComplexity(kmer) {
			baseWeight *= 0.5 // Reduce weight for rare k-mers
		} else {
			baseWeight *= 1.5 // Favor common k-mers
		}
	}
	
	return baseWeight
}

// wouldIncreaseComplexity checks if adding this k-mer would increase complexity
func (m *AdaptiveHMMModel) wouldIncreaseComplexity(kmer string) bool {
	// Check if k-mer is in recent memory
	for _, memKmer := range m.complexityMemory {
		if memKmer == kmer {
			return false // Would not increase complexity
		}
	}
	return true // Would increase complexity
}

// passesBiologicalFilters applies biological realism filters
func (m *AdaptiveHMMModel) passesBiologicalFilters(kmer string, window []string) bool {
	// Homopolymer filter
	if !m.passesHomopolymerFilter(kmer) {
		return false
	}
	
	// GC content filter
	if !m.passesGCFilter(kmer, window) {
		return false
	}
	
	// CpG depletion filter
	if !m.passesCpGFilter(kmer) {
		return false
	}
	
	return true
}

// passesHomopolymerFilter checks homopolymer run length
func (m *AdaptiveHMMModel) passesHomopolymerFilter(kmer string) bool {
	currentRun := 1
	for i := 1; i < len(kmer); i++ {
		if kmer[i] == kmer[i-1] {
			currentRun++
			if currentRun > m.homopolymerMaxLength {
				return false
			}
		} else {
			currentRun = 1
		}
	}
	return true
}

// passesGCFilter checks GC content constraints with enhanced control
func (m *AdaptiveHMMModel) passesGCFilter(kmer string, window []string) bool {
	// Calculate current window GC content
	gcCount := 0
	totalBases := 0
	
	for _, windowKmer := range window {
		for _, base := range windowKmer {
			totalBases++
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}
	}
	
	if totalBases == 0 {
		return true
	}
	
	currentGC := float64(gcCount) / float64(totalBases)
	
	// Check if adding this k-mer would improve GC content
	kmerGC := 0
	for _, base := range kmer {
		if base == 'G' || base == 'C' {
			kmerGC++
		}
	}
	
	kmerGCRatio := float64(kmerGC) / float64(len(kmer))
	
	// Enhanced GC content correction with tighter control
	deviationFromTarget := currentGC - m.gcContentTarget
	absDeviation := deviationFromTarget
	if absDeviation < 0 {
		absDeviation = -absDeviation
	}
	
	// If we're significantly outside tolerance, be very aggressive
	if absDeviation > m.gcContentTolerance * 2 {
		if deviationFromTarget > 0 {
			// Too GC-rich, strongly favor AT-rich k-mers
			return kmerGCRatio < 0.3
		} else {
			// Too AT-rich, strongly favor GC-rich k-mers
			return kmerGCRatio > 0.7
		}
	}
	
	// If we're moderately outside tolerance, be moderately aggressive
	if absDeviation > m.gcContentTolerance {
		if deviationFromTarget > 0 {
			// Too GC-rich, favor AT-rich k-mers
			return kmerGCRatio < 0.35
		} else {
			// Too AT-rich, favor GC-rich k-mers
			return kmerGCRatio > 0.65
		}
	}
	
	// If we're within tolerance, apply gentle corrections
	if deviationFromTarget > 0 {
		// Slightly too GC-rich, gently favor AT-rich k-mers
		return kmerGCRatio < 0.5
	} else if deviationFromTarget < 0 {
		// Slightly too AT-rich, gently favor GC-rich k-mers
		return kmerGCRatio > 0.5
	}
	
	return true
}

// passesCpGFilter applies CpG depletion constraints
func (m *AdaptiveHMMModel) passesCpGFilter(kmer string) bool {
	// Count CpG dinucleotides in k-mer
	cpgCount := 0
	for i := 0; i < len(kmer)-1; i++ {
		if kmer[i] == 'C' && kmer[i+1] == 'G' {
			cpgCount++
		}
	}
	
	// Apply CpG depletion
	if cpgCount > 0 {
		return m.rng.Float64() < m.cpgDepletionRatio
	}
	
	return true
}

// generateAdaptiveQuality generates quality scores based on complexity context
func (m *AdaptiveHMMModel) generateAdaptiveQuality(kmer string, complexity float64) int {
	baseQuality := 30
	
	// Higher complexity sequences tend to have slightly higher quality
	if complexity > 0.6 {
		baseQuality = 32
	} else if complexity < 0.3 {
		baseQuality = 28
	}
	
	// GC content affects quality
	gcCount := 0
	for _, base := range kmer {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}
	
	if gcCount > len(kmer)/2 {
		baseQuality -= 2
	}
	
	// Add randomness and clamp
	quality := baseQuality + int(m.rng.NormFloat64()*3)
	if quality < 0 {
		quality = 0
	}
	if quality >= 42 {
		quality = 41
	}
	
	return quality
}

// normalizeWithComplexityControl normalizes probabilities with complexity control
func (m *AdaptiveHMMModel) normalizeWithComplexityControl(
	transitionCounts, emissionCounts [][]int, startCounts []int) {
	
	// Normalize start probabilities
	startTotal := 0
	for _, count := range startCounts {
		startTotal += count
	}
	
	if startTotal > 0 {
		for i, count := range startCounts {
			m.startProbabilities[i] = float64(count) / float64(startTotal)
		}
	} else {
		for i := range m.startProbabilities {
			m.startProbabilities[i] = 1.0 / float64(len(m.startProbabilities))
		}
	}

	// Normalize transition probabilities with complexity awareness
	for i, counts := range transitionCounts {
		total := 0
		for _, count := range counts {
			total += count
		}
		
		if total > 0 {
			for j, count := range counts {
				baseProb := float64(count) / float64(total)
				
				// Apply adaptive complexity weighting
				state := m.stateSpace[i]
				nextBase := m.getBaseFromIndex(j)
				nextKmer := state + string(nextBase)
				
				complexityWeight := m.getComplexityWeight(nextKmer)
				biologicalWeight := m.getBiologicalWeight(state, nextBase)
				
				finalProb := baseProb * complexityWeight * biologicalWeight
				m.transitionMatrix[i][j] = finalProb
			}
			
			// Renormalize
			m.renormalizeTransitionRow(i)
		} else {
			// Uniform distribution
			for j := range m.transitionMatrix[i] {
				m.transitionMatrix[i][j] = 0.25
			}
		}
	}

	// Normalize emission probabilities
	for i, counts := range emissionCounts {
		total := 0
		for _, count := range counts {
			total += count
		}
		
		if total > 0 {
			for j, count := range counts {
				m.emissionMatrix[i][j] = float64(count) / float64(total)
			}
		} else {
			for j := range m.emissionMatrix[i] {
				m.emissionMatrix[i][j] = m.getRealisticQualityProb(j)
			}
		}
	}
}

// getComplexityWeight calculates weight based on complexity needs
func (m *AdaptiveHMMModel) getComplexityWeight(kmer string) float64 {
	// Favor k-mers that would increase complexity
	if weight, exists := m.adaptiveWeights[kmer]; exists {
		return weight
	}
	return 1.0
}

// getBiologicalWeight calculates biological realism weight with enhanced GC control
func (m *AdaptiveHMMModel) getBiologicalWeight(state string, nextBase byte) float64 {
	weight := 1.0
	
	// Enhanced GC content bias with adaptive correction
	// Calculate current local GC content from state
	localGC := 0.0
	if len(state) > 0 {
		gcCount := 0
		for _, base := range state {
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}
		localGC = float64(gcCount) / float64(len(state))
	}
	
	// Calculate deviation from target GC content
	gcDeviation := localGC - m.gcContentTarget
	
	// Apply adaptive GC content weighting
	if nextBase == 'G' || nextBase == 'C' {
		// GC base - reduce weight if we're already GC-rich
		baseGCWeight := m.gcContentTarget / 0.5
		if gcDeviation > m.gcContentTolerance {
			// Too GC-rich, penalize GC bases more heavily
			penalty := 1.0 + (gcDeviation - m.gcContentTolerance) * 10.0
			baseGCWeight /= penalty
		} else if gcDeviation < -m.gcContentTolerance {
			// Too AT-rich, favor GC bases more heavily
			bonus := 1.0 + (-gcDeviation - m.gcContentTolerance) * 5.0
			baseGCWeight *= bonus
		}
		weight *= baseGCWeight
	} else {
		// AT base - reduce weight if we're already AT-rich
		baseATWeight := (1.0 - m.gcContentTarget) / 0.5
		if gcDeviation < -m.gcContentTolerance {
			// Too AT-rich, penalize AT bases more heavily
			penalty := 1.0 + (-gcDeviation - m.gcContentTolerance) * 10.0
			baseATWeight /= penalty
		} else if gcDeviation > m.gcContentTolerance {
			// Too GC-rich, favor AT bases more heavily
			bonus := 1.0 + (gcDeviation - m.gcContentTolerance) * 5.0
			baseATWeight *= bonus
		}
		weight *= baseATWeight
	}
	
	// CpG depletion
	if len(state) > 0 && state[len(state)-1] == 'C' && nextBase == 'G' {
		weight *= m.cpgDepletionRatio
	}
	
	// Homopolymer penalty
	if len(state) > 0 && state[len(state)-1] == nextBase {
		weight *= 0.7
	}
	
	return weight
}

// GenerateSequence generates a sequence with adaptive complexity control
func (m *AdaptiveHMMModel) GenerateSequence(length int) (*GenerationResult, error) {
	// Read trained/k under RLock before entering the generation lock.
	m.mu.RLock()
	trained := m.trained
	k := m.k
	m.mu.RUnlock()

	if !trained {
		return nil, errors.NewValidationError("adaptive_hmm_model", "generate_sequence", "model must be trained before generating sequences")
	}

	if length < k {
		return nil, errors.NewValidationError("adaptive_hmm_model", "generate_sequence", fmt.Sprintf("sequence length must be at least k=%d", k))
	}

	start := time.Now()

	// Generate with adaptive complexity control.
	// generateAdaptiveSequence acquires a full write lock internally because it
	// writes m.complexityMemory during generation.
	sequence, qualities, states := m.generateAdaptiveSequence(length)

	// Calculate k-mer frequencies (local computation, no shared state).
	kmerFreqs := make(map[string]int)
	for i := 0; i <= len(sequence)-k; i++ {
		kmer := string(sequence[i : i+k])
		kmerFreqs[kmer]++
	}

	// Update shared generation stats under a write lock.
	m.mu.Lock()
	m.updateGenerationStats(sequence, qualities)
	m.mu.Unlock()

	return &GenerationResult{
		Sequence:        sequence,
		Qualities:       qualities,
		KmerFrequencies: kmerFreqs,
		GenerationTime:  time.Since(start),
		States:          states,
	}, nil
}

// generateAdaptiveSequence generates sequence with adaptive complexity control.
// Must be called without holding m.mu — acquires a full write lock because it
// writes m.complexityMemory throughout the generation loop.
func (m *AdaptiveHMMModel) generateAdaptiveSequence(length int) ([]byte, []byte, []int) {
	m.mu.Lock()
	defer m.mu.Unlock()

	sequence := make([]byte, length)
	qualities := make([]byte, length)
	states := make([]int, length)

	// Initialize complexity tracking
	m.complexityMemory = make([]string, 0)
	
	// Choose initial state
	currentState := m.sampleFromDistribution(m.startProbabilities)
	states[0] = currentState

	// Generate initial k-mer
	initialKmer := m.stateSpace[currentState]
	for i, base := range initialKmer {
		sequence[i] = byte(base) // #nosec G115
		qualities[i] = byte(m.generateAdaptiveQuality(initialKmer, 0.5)) // #nosec G115
	}

	// Add initial k-mer to memory
	if len(initialKmer) == m.k {
		m.updateComplexityMemory(initialKmer)
	}

	// Generate remaining sequence with adaptive complexity control
	for i := m.k - 1; i < length; i++ {
		// Calculate current sequence complexity
		currentComplexity := m.calculateCurrentComplexity(sequence, i)
		
		// Generate next base with adaptive control
		nextBase := m.generateAdaptiveBase(currentState, sequence, i, currentComplexity)
		nextQuality := m.generateAdaptiveQuality(string(nextBase), currentComplexity)

		sequence[i] = nextBase
		qualities[i] = byte(nextQuality) // #nosec G115

		// Update complexity memory
		if i >= m.k-1 {
			newKmer := string(sequence[i-m.k+1 : i+1])
			m.updateComplexityMemory(newKmer)
		}

		// Update state
		if i < length-1 {
			newState := string(sequence[i-m.k+2 : i+1])
			currentState = m.getStateIndex(newState)
			if currentState == -1 {
				currentState = 0
			}
			states[i] = currentState
		}
	}

	return sequence, qualities, states
}

// calculateCurrentComplexity calculates complexity of current sequence
func (m *AdaptiveHMMModel) calculateCurrentComplexity(sequence []byte, pos int) float64 {
	if pos < m.k {
		return 0.5 // Default complexity
	}
	
	// Look at recent window
	windowStart := pos - m.complexityWindowSize
	if windowStart < 0 {
		windowStart = 0
	}
	
	windowEnd := pos
	if windowEnd > len(sequence) {
		windowEnd = len(sequence)
	}
	
	windowSeq := sequence[windowStart:windowEnd]
	
	// Count unique k-mers in window
	uniqueKmers := make(map[string]bool)
	for i := 0; i <= len(windowSeq)-m.k; i++ {
		kmer := string(windowSeq[i : i+m.k])
		uniqueKmers[kmer] = true
	}
	
	totalPossibleKmers := len(windowSeq) - m.k + 1
	if totalPossibleKmers <= 0 {
		return 0.5
	}
	
	return float64(len(uniqueKmers)) / float64(totalPossibleKmers)
}

// generateAdaptiveBase generates next base with adaptive complexity control
func (m *AdaptiveHMMModel) generateAdaptiveBase(state int, sequence []byte, pos int, complexity float64) byte {
	transitionProbs := make([]float64, 4)
	copy(transitionProbs, m.transitionMatrix[state])

	// Apply adaptive complexity adjustments
	for i := 0; i < 4; i++ {
		nextBase := []byte{'A', 'T', 'G', 'C'}[i]
		
		// Check if this base would create a k-mer that increases complexity
		if pos >= m.k-1 {
			potentialKmer := string(sequence[pos-m.k+1:pos]) + string(nextBase)
			
			// If complexity is too low, favor new k-mers
			if complexity < m.targetComplexity {
				if m.wouldIncreaseComplexity(potentialKmer) {
					transitionProbs[i] *= 1.5
				} else {
					transitionProbs[i] *= 0.7
				}
			}
			
			// If complexity is too high, favor repeated k-mers
			if complexity > m.targetComplexity + 0.1 {
				if m.wouldIncreaseComplexity(potentialKmer) {
					transitionProbs[i] *= 0.7
				} else {
					transitionProbs[i] *= 1.3
				}
			}
		}
		
		// Apply biological constraints
		if pos > 0 {
			prevBase := sequence[pos-1]
			if nextBase == prevBase {
				// Check homopolymer length
				runLength := 1
				for j := pos - 2; j >= 0 && sequence[j] == prevBase; j-- {
					runLength++
				}
				
				if runLength >= 5 {
					transitionProbs[i] *= 0.2
				}
			}
		}
	}

	// Renormalize
	sum := 0.0
	for _, prob := range transitionProbs {
		sum += prob
	}
	if sum > 0 {
		for i := range transitionProbs {
			transitionProbs[i] /= sum
		}
	}

	baseIndex := m.sampleFromDistribution(transitionProbs)
	return []byte{'A', 'T', 'G', 'C'}[baseIndex]
}

// updateComplexityMemory updates the complexity memory window
func (m *AdaptiveHMMModel) updateComplexityMemory(kmer string) {
	m.complexityMemory = append(m.complexityMemory, kmer)
	if len(m.complexityMemory) > m.complexityWindowSize {
		m.complexityMemory = m.complexityMemory[1:]
	}
}

// updateGenerationStats updates generation statistics
func (m *AdaptiveHMMModel) updateGenerationStats(sequence, qualities []byte) {
	m.generationStats.TotalSequences++
	m.generationStats.TotalBases += int64(len(sequence))
	
	// Calculate complexity
	complexity := m.calculateSequenceComplexity(sequence)
	m.generationStats.AverageComplexity = (m.generationStats.AverageComplexity*float64(m.generationStats.TotalSequences-1) + complexity) / float64(m.generationStats.TotalSequences)
	
	// Calculate GC content
	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}
	gcContent := float64(gcCount) / float64(len(sequence))
	m.generationStats.AverageGCContent = (m.generationStats.AverageGCContent*float64(m.generationStats.TotalSequences-1) + gcContent) / float64(m.generationStats.TotalSequences)
	
	// Calculate max homopolymer
	maxHomopolymer := m.findMaxHomopolymerRun(sequence)
	if maxHomopolymer > m.generationStats.MaxHomopolymer {
		m.generationStats.MaxHomopolymer = maxHomopolymer
	}
}

// calculateSequenceComplexity calculates overall sequence complexity
func (m *AdaptiveHMMModel) calculateSequenceComplexity(sequence []byte) float64 {
	if len(sequence) < m.k {
		return 0.0
	}
	
	uniqueKmers := make(map[string]bool)
	for i := 0; i <= len(sequence)-m.k; i++ {
		kmer := string(sequence[i : i+m.k])
		uniqueKmers[kmer] = true
	}
	
	totalKmers := len(sequence) - m.k + 1
	return float64(len(uniqueKmers)) / float64(totalKmers)
}

// Utility methods

func (m *AdaptiveHMMModel) K() int {
	return m.k
}

func (m *AdaptiveHMMModel) IsTrained() bool {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.trained
}

func (m *AdaptiveHMMModel) TotalKmers() uint64 {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.totalKmers
}

func (m *AdaptiveHMMModel) StateSpace() []string {
	return m.stateSpace
}

func (m *AdaptiveHMMModel) SetSeed(seed int64) {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.rng = rand.New(rand.NewPCG(uint64(seed), 0)) // #nosec G115 G404
}

func (m *AdaptiveHMMModel) ValidateModel() error {
	m.mu.RLock()
	defer m.mu.RUnlock()
	
	if !m.trained {
		return errors.NewValidationError("adaptive_hmm_model", "validate_model", "model has not been trained")
	}
	
	return nil
}

func (m *AdaptiveHMMModel) GetGenerationStats() GenerationStats {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.generationStats
}

// Helper methods

func (m *AdaptiveHMMModel) isValidKmer(kmer string) bool {
	for _, base := range kmer {
		if base != 'A' && base != 'T' && base != 'G' && base != 'C' {
			return false
		}
	}
	return true
}

func (m *AdaptiveHMMModel) getStateIndex(state string) int {
	for i, s := range m.stateSpace {
		if s == state {
			return i
		}
	}
	return -1
}

func (m *AdaptiveHMMModel) getBaseIndex(base byte) int {
	switch base {
	case 'A': return 0
	case 'T': return 1
	case 'G': return 2
	case 'C': return 3
	default: return -1
	}
}

func (m *AdaptiveHMMModel) getBaseFromIndex(index int) byte {
	bases := []byte{'A', 'T', 'G', 'C'}
	if index >= 0 && index < len(bases) {
		return bases[index]
	}
	return 'N'
}

func (m *AdaptiveHMMModel) sampleFromDistribution(probs []float64) int {
	r := m.rng.Float64()
	cumulative := 0.0
	
	for i, prob := range probs {
		cumulative += prob
		if r <= cumulative {
			return i
		}
	}
	
	return len(probs) - 1
}

func (m *AdaptiveHMMModel) getRealisticQualityProb(quality int) float64 {
	if quality >= 20 && quality <= 35 {
		return 0.8 / 16.0
	}
	return 0.2 / 26.0
}

func (m *AdaptiveHMMModel) renormalizeTransitionRow(row int) {
	sum := 0.0
	for _, prob := range m.transitionMatrix[row] {
		sum += prob
	}
	
	if sum > 0 {
		for j := range m.transitionMatrix[row] {
			m.transitionMatrix[row][j] /= sum
		}
	}
}

func (m *AdaptiveHMMModel) findMaxHomopolymerRun(sequence []byte) int {
	if len(sequence) == 0 {
		return 0
	}
	
	maxRun := 1
	currentRun := 1
	
	for i := 1; i < len(sequence); i++ {
		if sequence[i] == sequence[i-1] {
			currentRun++
		} else {
			if currentRun > maxRun {
				maxRun = currentRun
			}
			currentRun = 1
		}
	}
	
	if currentRun > maxRun {
		maxRun = currentRun
	}
	
	return maxRun
}

// TrainFromSequence trains the model on a raw sequence
func (m *AdaptiveHMMModel) TrainFromSequence(sequence []byte) error {
	if len(sequence) < m.k {
		return nil
	}
	
	m.mu.Lock()
	defer m.mu.Unlock()
	
	// Process sequence and update model
	for i := 0; i <= len(sequence)-m.k; i++ {
		kmer := string(sequence[i : i+m.k])
		if !m.isValidKmer(kmer) {
			continue
		}
		
		// Update k-mer counts
		m.totalKmers++
		
		// Update complexity memory
		m.updateComplexityMemory(kmer)
		
		// Train transitions if we have a previous k-mer
		if i > 0 {
			prevKmer := string(sequence[i-1 : i-1+m.k])
			if m.isValidKmer(prevKmer) {
				// Update transition probabilities
				prevStateIndex := m.getStateIndex(prevKmer[:m.k-1])
				nextBase := kmer[m.k-1]
				nextBaseIndex := m.getBaseIndex(nextBase)
				
				if prevStateIndex >= 0 && nextBaseIndex >= 0 {
					m.transitionMatrix[prevStateIndex][nextBaseIndex] += 1.0
				}
			}
		}
	}
	
	m.trained = true
	return nil
}