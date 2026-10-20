package authenticity

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/errors"
)

// AdaptiveHMMModel is an order-(k-1) Markov chain like HMMKmerModel (no
// state is hidden despite the name), with heuristic controls at generation
// time.
//
// Training (TrainFromSequence, TrainFromReferenceGenome) accumulates raw
// transition counts: every k-mer of A/C/G/T counts one observation of its
// last base following its first k-1 bases. Training may be split across any
// number of calls; for a sequence cut into chunks, carrying the last k-1
// bases of each chunk into the next gives exactly the counts of one-shot
// training. Finalize turns the counts into probabilities with pseudocount
// PseudoCount (see markov.go): every row sums to 1, unseen contexts are
// uniform, and the start distribution is the smoothed (k-1)-mer frequency.
// Finalize runs automatically before generation and on Save.
//
// Generation samples its first k-1 bases from the start distribution and
// then each base from the transition row of the preceding k-1 bases,
// reweighted toward a target fraction of distinct k-mers in a 100-k-mer
// window and against homopolymer runs of 5 or more, so generated sequences
// are not an unbiased sample of the chain. Qualities come from
// generateAdaptiveQuality, a heuristic, not from training data.
type AdaptiveHMMModel struct {
	k                  int
	stateSpace         []string
	transitionCounts   [][]float64 // raw counts c(s, b) from training
	transitionMatrix   [][]float64 // P(b | s), valid when finalized
	startProbabilities []float64   // π(s), valid when finalized
	finalized          bool        // probabilities are up to date with the counts

	// Generation-time complexity control
	complexityMemory     []string // Recent k-mer history for complexity tracking
	complexityWindowSize int      // Size of complexity tracking window
	targetComplexity     float64  // Target fraction of distinct k-mers in the window

	// Training and generation state
	totalKmers uint64
	trained    bool
	mu         sync.RWMutex
	rng        *rand.Rand

	// Performance metrics
	generationStats GenerationStats
}

// GenerationStats accumulates summary statistics over generated sequences.
// AverageJSDivergence is never computed and stays zero.
type GenerationStats struct {
	TotalSequences      int
	TotalBases          int64
	AverageComplexity   float64
	AverageGCContent    float64
	MaxHomopolymer      int
	AverageJSDivergence float64
}

// NewAdaptiveHMMModel creates a new, untrained adaptive model of order k-1
// (2 <= k <= 8); it returns nil for any other k. Its random source starts
// from DefaultSeed (see SetSeed).
func NewAdaptiveHMMModel(k int) *AdaptiveHMMModel {
	if k < 2 || k > 8 {
		// Use proper error handling instead of panic
		err := errors.NewValidationError("adaptive_hmm_model", "new_model", fmt.Sprintf("k must be between 2 and 8, got %d", k))
		errors.HandleError(err)
		return nil
	}

	model := &AdaptiveHMMModel{
		k:                    k,
		rng:                  newModelRNG(DefaultSeed),
		complexityWindowSize: 100, // Track last 100 k-mers
		targetComplexity:     0.7, // Target fraction of distinct k-mers in the window
		complexityMemory:     make([]string, 0),
	}

	model.initializeStateSpace()
	return model
}

// initializeStateSpace generates all (k-1)-mers as states, in state-index
// order, and allocates the count and probability matrices.
func (m *AdaptiveHMMModel) initializeStateSpace() {
	m.stateSpace = generateStateSpace(m.k)
	numStates := len(m.stateSpace)
	m.transitionCounts = newCountMatrix(numStates, 4)
	m.transitionMatrix = newCountMatrix(numStates, 4)
	m.startProbabilities = make([]float64, numStates)
}

// generateStateSpace returns every (k-1)-mer over A,T,G,C in state-index
// order (first base most significant).
func generateStateSpace(k int) []string {
	states := []string{""}
	for i := 0; i < k-1; i++ {
		next := make([]string, 0, len(states)*4)
		for _, s := range states {
			for _, b := range basesATGC {
				next = append(next, s+string(b))
			}
		}
		states = next
	}
	return states
}

// TrainFromSequence adds the k-mer transitions of sequence to the model's
// counts. k-mers containing anything other than A/C/G/T (case-insensitive)
// are skipped. A sequence shorter than k adds nothing. To train on a long
// sequence in chunks, prefix each chunk with the last k-1 bases of the
// previous one.
func (m *AdaptiveHMMModel) TrainFromSequence(sequence []byte) error {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.addCounts(sequence)
	return nil
}

// addCounts counts sequence's k-mers. The caller must hold m.mu.
func (m *AdaptiveHMMModel) addCounts(sequence []byte) {
	if n := countKmers(sequence, m.k, m.transitionCounts); n > 0 {
		m.totalKmers += n
		m.trained = true
		m.finalized = false
	}
}

// TrainFromReferenceGenome adds the k-mer transitions of every chromosome of
// refGenome to the counts (chromosomes are independent: no k-mer spans two)
// and finalizes the model. Chromosomes that fail to load are skipped.
func (m *AdaptiveHMMModel) TrainFromReferenceGenome(refGenome GenomeInterface) error {
	m.mu.Lock()
	defer m.mu.Unlock()

	for _, chrName := range refGenome.GetChromosomeNames() {
		chr, err := refGenome.LoadChromosome(chrName)
		if err != nil {
			continue
		}
		m.addCounts(chr.GetSequence())
	}
	if m.trained {
		m.finalizeLocked()
	}
	return nil
}

// Finalize estimates the transition and start probabilities from the counts
// accumulated so far (see markov.go). It is idempotent; training after
// Finalize marks the model stale, and the next Finalize, GenerateSequence or
// Save re-estimates from all counts.
func (m *AdaptiveHMMModel) Finalize() {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.finalizeLocked()
}

// finalizeLocked is Finalize with m.mu held for writing.
func (m *AdaptiveHMMModel) finalizeLocked() {
	if m.finalized {
		return
	}
	estimateChain(m.transitionCounts, m.transitionMatrix, m.startProbabilities)
	m.finalized = true
}

// wouldIncreaseComplexity reports whether kmer is absent from the recent
// k-mer window.
func (m *AdaptiveHMMModel) wouldIncreaseComplexity(kmer string) bool {
	for _, memKmer := range m.complexityMemory {
		if memKmer == kmer {
			return false
		}
	}
	return true
}

// generateAdaptiveQuality returns a heuristic quality: 30 (32 above 0.6
// complexity, 28 below 0.3), minus 2 for GC-rich k-mers, plus Gaussian noise
// with SD 3, clamped to [0, 41]. The adjustments are not fitted to data.
func (m *AdaptiveHMMModel) generateAdaptiveQuality(kmer string, complexity float64) int {
	baseQuality := 30

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

// GenerateSequence generates a sequence with adaptive complexity control.
// It finalizes the model first if training has added counts since the last
// Finalize.
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

	// generateAdaptiveSequence acquires a full write lock internally because it
	// may finalize and writes m.complexityMemory during generation.
	sequence, qualities, states := m.generateAdaptiveSequence(length)

	// Calculate k-mer frequencies (local computation, no shared state).
	kmerFreqs := make(map[string]int)
	for i := 0; i <= len(sequence)-k; i++ {
		kmer := string(sequence[i : i+k])
		kmerFreqs[kmer]++
	}

	// Update shared generation stats under a write lock.
	m.mu.Lock()
	m.updateGenerationStats(sequence)
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
// Must be called without holding m.mu. States[j] is the state (context
// x[j:j+k-1]) that emitted base j+k-1.
func (m *AdaptiveHMMModel) generateAdaptiveSequence(length int) ([]byte, []byte, []int) {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.finalizeLocked()

	k := m.k
	sequence := make([]byte, length)
	qualities := make([]byte, length)
	states := make([]int, 0, length-k+1)

	m.complexityMemory = make([]string, 0, m.complexityWindowSize)

	// Initial k-1 bases from the start distribution.
	currentState := sampleIndex(m.rng, m.startProbabilities)
	initialKmer := m.stateSpace[currentState]
	for i := 0; i < len(initialKmer); i++ {
		sequence[i] = initialKmer[i]
		qualities[i] = byte(m.generateAdaptiveQuality(initialKmer, 0.5)) // #nosec G115
	}

	for i := k - 1; i < length; i++ {
		states = append(states, currentState)
		currentComplexity := m.calculateCurrentComplexity(sequence, i)

		nextBase := m.generateAdaptiveBase(currentState, sequence, i, currentComplexity)
		sequence[i] = nextBase
		qualities[i] = byte(m.generateAdaptiveQuality(string(nextBase), currentComplexity)) // #nosec G115

		m.updateComplexityMemory(string(sequence[i-k+1 : i+1]))
		currentState = contextIndex(sequence[i-k+2 : i+1])
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

// generateAdaptiveBase draws the next base from the (smoothed, strictly
// positive) transition row of state, reweighted for complexity and against
// long homopolymers. Every base keeps a positive weight.
func (m *AdaptiveHMMModel) generateAdaptiveBase(state int, sequence []byte, pos int, complexity float64) byte {
	var transitionProbs [4]float64
	copy(transitionProbs[:], m.transitionMatrix[state])

	for i, nextBase := range basesATGC {
		// Check if this base would create a k-mer that increases complexity
		if pos >= m.k-1 {
			potentialKmer := string(sequence[pos-m.k+1:pos]) + string(nextBase)

			// If complexity is too low, favor new k-mers
			if complexity < m.targetComplexity {
				if m.wouldIncreaseComplexity(potentialKmer) {
					transitionProbs[i] = float64(transitionProbs[i] * 1.5)
				} else {
					transitionProbs[i] = float64(transitionProbs[i] * 0.7)
				}
			}

			// If complexity is too high, favor repeated k-mers
			if complexity > m.targetComplexity+0.1 {
				if m.wouldIncreaseComplexity(potentialKmer) {
					transitionProbs[i] = float64(transitionProbs[i] * 0.7)
				} else {
					transitionProbs[i] = float64(transitionProbs[i] * 1.3)
				}
			}
		}

		// Penalise extending a homopolymer run of 5 or more.
		if pos > 0 {
			prevBase := sequence[pos-1]
			if nextBase == prevBase {
				runLength := 1
				for j := pos - 2; j >= 0 && sequence[j] == prevBase; j-- {
					runLength++
				}
				if runLength >= 5 {
					transitionProbs[i] = float64(transitionProbs[i] * 0.2)
				}
			}
		}
	}

	return basesATGC[sampleIndex(m.rng, transitionProbs[:])]
}

// updateComplexityMemory updates the complexity memory window
func (m *AdaptiveHMMModel) updateComplexityMemory(kmer string) {
	m.complexityMemory = append(m.complexityMemory, kmer)
	if len(m.complexityMemory) > m.complexityWindowSize {
		m.complexityMemory = m.complexityMemory[1:]
	}
}

// updateGenerationStats updates generation statistics
func (m *AdaptiveHMMModel) updateGenerationStats(sequence []byte) {
	m.generationStats.TotalSequences++
	m.generationStats.TotalBases += int64(len(sequence))

	// Calculate complexity
	complexity := m.calculateSequenceComplexity(sequence)
	m.generationStats.AverageComplexity = (float64(m.generationStats.AverageComplexity*float64(m.generationStats.TotalSequences-1)) + complexity) / float64(m.generationStats.TotalSequences)

	// Calculate GC content
	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}
	gcContent := float64(gcCount) / float64(len(sequence))
	m.generationStats.AverageGCContent = (float64(m.generationStats.AverageGCContent*float64(m.generationStats.TotalSequences-1)) + gcContent) / float64(m.generationStats.TotalSequences)

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

// K returns the k-mer size (the model has order k-1).
func (m *AdaptiveHMMModel) K() int {
	return m.k
}

// IsTrained reports whether at least one k-mer has been counted.
func (m *AdaptiveHMMModel) IsTrained() bool {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.trained
}

// TotalKmers returns the number of k-mers counted in training.
func (m *AdaptiveHMMModel) TotalKmers() uint64 {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.totalKmers
}

// StateSpace returns the (k-1)-mer states in state-index order.
func (m *AdaptiveHMMModel) StateSpace() []string {
	return m.stateSpace
}

// SetSeed reseeds the model's random source; the same seed gives the same
// generated sequences.
func (m *AdaptiveHMMModel) SetSeed(seed int64) {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.rng = newModelRNG(seed)
}

// ValidateModel checks that the model is trained, its matrices have the
// right shape, and (once finalized) that every transition row and the start
// distribution sum to 1.
func (m *AdaptiveHMMModel) ValidateModel() error {
	m.mu.RLock()
	defer m.mu.RUnlock()

	if !m.trained {
		return errors.NewValidationError("adaptive_hmm_model", "validate_model", "model has not been trained")
	}
	return validateChain("adaptive_hmm_model", m.k, m.stateSpace, m.transitionCounts, m.transitionMatrix, m.startProbabilities, m.finalized)
}

// validateChain checks the shapes of a chain's matrices and, when
// checkProbs is set, that its rows and start distribution sum to 1.
func validateChain(component string, k int, stateSpace []string, counts, trans [][]float64, start []float64, checkProbs bool) error {
	numStates := 1 << (2 * (k - 1))
	if len(stateSpace) != numStates || len(trans) != numStates || len(start) != numStates ||
		(counts != nil && len(counts) != numStates) {
		return errors.NewValidationError(component, "validate_model", fmt.Sprintf("expected %d states for k=%d", numStates, k))
	}
	for i, s := range stateSpace {
		if contextIndex([]byte(s)) != i || len(s) != k-1 {
			return errors.NewValidationError(component, "validate_model", fmt.Sprintf("state %d is %q, not in state-index order", i, s))
		}
	}
	for i := range trans {
		if len(trans[i]) != 4 || (counts != nil && len(counts[i]) != 4) {
			return errors.NewValidationError(component, "validate_model", fmt.Sprintf("state %d does not have 4 transitions", i))
		}
	}
	if !checkProbs {
		return nil
	}
	const tol = 1e-9
	for i, row := range trans {
		sum := 0.0
		for _, p := range row {
			if p < 0 || math.IsNaN(p) {
				return errors.NewValidationError(component, "validate_model", fmt.Sprintf("state %d has an invalid transition probability", i))
			}
			sum += p
		}
		if math.Abs(sum-1) > tol {
			return errors.NewValidationError(component, "validate_model", fmt.Sprintf("transition probabilities for state %d sum to %f, not 1", i, sum))
		}
	}
	sum := 0.0
	for _, p := range start {
		if p < 0 || math.IsNaN(p) {
			return errors.NewValidationError(component, "validate_model", "invalid start probability")
		}
		sum += p
	}
	if math.Abs(sum-1) > tol {
		return errors.NewValidationError(component, "validate_model", fmt.Sprintf("start probabilities sum to %f, not 1", sum))
	}
	return nil
}

// GetGenerationStats returns the statistics accumulated over all generated
// sequences.
func (m *AdaptiveHMMModel) GetGenerationStats() GenerationStats {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.generationStats
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
