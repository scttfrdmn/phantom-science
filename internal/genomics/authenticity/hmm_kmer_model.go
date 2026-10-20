package authenticity

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/errors"
)

// HMMKmerModel is an order-(k-1) Markov chain for sequence generation. Each
// state is the preceding k-1 bases and transitions give the next base, so no
// state is hidden despite the name. Transition and start probabilities are
// estimated from k-mer counts with pseudocount PseudoCount (see markov.go).
// Each state also has a quality distribution, trained from simulateQuality's
// heuristic values rather than real base qualities.
type HMMKmerModel struct {
	k                  int
	stateSpace         []string     // All possible (k-1)-mers
	transitionMatrix   [][]float64  // State -> base probabilities
	emissionMatrix     [][]float64  // State -> quality probabilities
	startProbabilities []float64    // Initial state probabilities
	totalKmers         uint64       // Total k-mers processed during training
	trained            bool         // Whether model has been trained
	mu                 sync.RWMutex // Thread safety
	rng                *rand.Rand   // Random number generator
}

// KmerState represents a state in the HMM (k-1 mer)
type KmerState struct {
	Index       int
	Kmer        string
	Transitions map[byte]float64 // Base -> probability
	Emissions   map[int]float64  // Quality -> probability
}

// GenerationResult represents the result of sequence generation
type GenerationResult struct {
	Sequence        []byte
	Qualities       []byte
	KmerFrequencies map[string]int
	GenerationTime  time.Duration
	States          []int // State sequence for debugging
}

// NewHMMKmerModel creates a new HMM k-mer model
func NewHMMKmerModel(k int) *HMMKmerModel {
	if k < 2 || k > 8 {
		// Use proper error handling instead of panic
		err := errors.NewValidationError("hmm_kmer_model", "new_model", fmt.Sprintf("k must be between 2 and 8, got %d", k))
		errors.HandleError(err)
		return nil
	}

	model := &HMMKmerModel{
		k:   k,
		rng: rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
	}

	// Initialize state space (all possible k-1 mers)
	model.stateSpace = model.generateStateSpace()

	// Initialize matrices
	numStates := len(model.stateSpace)
	model.transitionMatrix = make([][]float64, numStates)
	model.emissionMatrix = make([][]float64, numStates)
	model.startProbabilities = make([]float64, numStates)

	for i := 0; i < numStates; i++ {
		model.transitionMatrix[i] = make([]float64, 4) // 4 bases: A, T, G, C
		model.emissionMatrix[i] = make([]float64, 42)  // Phred scores 0-41
	}

	return model
}

// generateStateSpace generates all possible (k-1)-mers in state-index order.
func (m *HMMKmerModel) generateStateSpace() []string {
	return generateStateSpace(m.k)
}

// getStateIndex returns the index of a (k-1)-mer state, or -1 if it is not one.
func (m *HMMKmerModel) getStateIndex(state string) int {
	if len(state) != m.k-1 {
		return -1
	}
	return contextIndex([]byte(state))
}

// GenomeInterface defines the interface for genome data sources
type GenomeInterface interface {
	GetChromosomeNames() []string
	LoadChromosome(name string) (ChromosomeInterface, error)
}

// ChromosomeInterface defines the interface for chromosome data
type ChromosomeInterface interface {
	GetSequence() []byte
	GetName() string
	GetLength() int
}

// TrainFromReferenceGenome trains the model on every chromosome of refGenome,
// replacing any earlier training. Each k-mer of A/C/G/T (case-insensitive)
// within a chromosome counts one transition; k-mers spanning two chromosomes
// are not formed. Transition and start probabilities are then estimated with
// pseudocount PseudoCount (see markov.go), so every row sums to 1, unseen
// contexts are uniform and the start distribution is the smoothed (k-1)-mer
// frequency. Chromosomes that fail to load are skipped; it is an error if no
// k-mer was counted.
func (m *HMMKmerModel) TrainFromReferenceGenome(refGenome GenomeInterface) error {
	m.mu.Lock()
	defer m.mu.Unlock()

	numStates := len(m.stateSpace)
	transitionCounts := newCountMatrix(numStates, 4)
	emissionCounts := make([][]int, numStates)
	for i := range emissionCounts {
		emissionCounts[i] = make([]int, 42)
	}

	var total uint64
	for _, chrName := range refGenome.GetChromosomeNames() {
		chr, err := refGenome.LoadChromosome(chrName)
		if err != nil {
			continue // Skip problematic chromosomes
		}
		total += m.trainOnChromosome(chr, transitionCounts, emissionCounts)
	}
	if total == 0 {
		return errors.NewValidationError("hmm_kmer_model", "train_from_reference", fmt.Sprintf("no valid %d-mers in the reference", m.k))
	}

	estimateChain(transitionCounts, m.transitionMatrix, m.startProbabilities)
	m.normalizeEmissions(emissionCounts)
	m.totalKmers = total
	m.trained = true
	return nil
}

// trainOnChromosome counts the k-mers of one chromosome and returns how many
// were counted.
func (m *HMMKmerModel) trainOnChromosome(chr ChromosomeInterface, transitionCounts [][]float64, emissionCounts [][]int) uint64 {
	sequence := chr.GetSequence()
	var n uint64
	for i := 0; i+m.k <= len(sequence); i++ {
		kmer := sequence[i : i+m.k]
		stateIndex := contextIndex(kmer[:m.k-1])
		baseIndex := baseIndex(kmer[m.k-1])
		if stateIndex < 0 || baseIndex < 0 {
			continue // k-mer with N or another non-ACGT byte
		}
		transitionCounts[stateIndex][baseIndex]++

		// No real quality data is used: qualities come from simulateQuality.
		emissionCounts[stateIndex][m.simulateQuality(string(kmer))]++
		n++
	}
	return n
}

// normalizeEmissions converts quality counts to per-state distributions with
// a pseudocount of 0.1 per quality value.
func (m *HMMKmerModel) normalizeEmissions(emissionCounts [][]int) {
	const smoothingFactor = 0.1
	for i := range emissionCounts {
		total := 0.0
		for j := 0; j < 42; j++ {
			total += float64(emissionCounts[i][j]) + smoothingFactor
		}
		for j := 0; j < 42; j++ {
			m.emissionMatrix[i][j] = (float64(emissionCounts[i][j]) + smoothingFactor) / total
		}
	}
}

// GenerateSequence samples a sequence of the given length from the trained
// transition probabilities, with a quality per base sampled from the state's
// quality distribution.
func (m *HMMKmerModel) GenerateSequence(length int) (*GenerationResult, error) {
	// Generation advances m.rng, so it needs the write lock.
	m.mu.Lock()
	defer m.mu.Unlock()

	if !m.trained {
		return nil, errors.NewValidationError("hmm_kmer_model", "generate_sequence", "model must be trained before generating sequences")
	}

	if length < m.k {
		return nil, errors.NewValidationError("hmm_kmer_model", "generate_sequence", fmt.Sprintf("sequence length %d must be at least k=%d", length, m.k))
	}

	startTime := time.Now()

	sequence := make([]byte, length)
	qualities := make([]byte, length)
	states := make([]int, length-m.k+1)
	kmerFreqs := make(map[string]int)

	// Choose initial state
	currentState := m.chooseInitialState()
	states[0] = currentState

	// Generate first k-1 bases from initial state
	initialKmer := m.stateSpace[currentState]
	for i := 0; i < len(initialKmer); i++ {
		sequence[i] = initialKmer[i]
		qualities[i] = byte(m.sampleQuality(currentState)) // #nosec G115
	}

	// Generate remaining bases
	for i := m.k - 1; i < length; i++ {
		// Sample next base from transition probabilities
		nextBase := m.sampleNextBase(currentState)
		sequence[i] = nextBase
		qualities[i] = byte(m.sampleQuality(currentState)) // #nosec G115

		// Update k-mer frequency tracking
		if i >= m.k-1 {
			kmer := string(sequence[i-m.k+1 : i+1])
			kmerFreqs[kmer]++
		}

		// Update state for next iteration
		if i < length-1 {
			nextState := m.getNextState(currentState, nextBase)
			currentState = nextState
			stateIndex := i - m.k + 2
			if stateIndex < len(states) {
				states[stateIndex] = currentState
			}
		}
	}

	generationTime := time.Since(startTime)

	return &GenerationResult{
		Sequence:        sequence,
		Qualities:       qualities,
		KmerFrequencies: kmerFreqs,
		GenerationTime:  generationTime,
		States:          states,
	}, nil
}

// chooseInitialState draws the initial state from the start probabilities.
func (m *HMMKmerModel) chooseInitialState() int {
	return sampleIndex(m.rng, m.startProbabilities)
}

// sampleNextBase draws the next base from the state's transition row.
func (m *HMMKmerModel) sampleNextBase(stateIndex int) byte {
	return basesATGC[sampleIndex(m.rng, m.transitionMatrix[stateIndex])]
}

// sampleQuality draws a quality score from the state's emission row.
func (m *HMMKmerModel) sampleQuality(stateIndex int) int {
	return sampleIndex(m.rng, m.emissionMatrix[stateIndex])
}

// getNextState determines the next state given current state and observed base
func (m *HMMKmerModel) getNextState(currentState int, observedBase byte) int {
	if m.k == 2 {
		// For k=2, state is just the observed base
		return m.getBaseIndex(observedBase)
	}

	// For k>2, slide the window
	currentKmer := m.stateSpace[currentState]
	nextKmer := currentKmer[1:] + string(observedBase)

	return m.getStateIndex(nextKmer)
}

// getBaseIndex returns the index of a base (A=0, T=1, G=2, C=3), or -1.
func (m *HMMKmerModel) getBaseIndex(base byte) int {
	return baseIndex(base)
}

// simulateQuality returns a heuristic quality: 30 minus a GC-deviation
// penalty and a homopolymer penalty, clamped to [0, 41]. (An earlier version
// also subtracted 0.01 per base of offset into the training chromosome, which
// is not a read position: past about 3,000 bases every quality was 0.)
func (m *HMMKmerModel) simulateQuality(kmer string) byte {
	baseQuality := 30.0

	// GC bias effect
	gcContent := m.calculateGCContent(kmer)
	gcBias := math.Abs(gcContent-0.5) * 10.0

	// Homopolymer effect
	homopolymerPenalty := m.calculateHomopolymerPenalty(kmer)

	finalQuality := baseQuality - gcBias - homopolymerPenalty

	// Clamp to valid range
	if finalQuality < 0 {
		finalQuality = 0
	}
	if finalQuality > 41 {
		finalQuality = 41
	}

	return byte(finalQuality)
}

// calculateGCContent calculates GC content of a k-mer
func (m *HMMKmerModel) calculateGCContent(kmer string) float64 {
	gcCount := 0
	for _, base := range kmer {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}
	return float64(gcCount) / float64(len(kmer))
}

// calculateHomopolymerPenalty calculates quality penalty for homopolymer runs
func (m *HMMKmerModel) calculateHomopolymerPenalty(kmer string) float64 {
	maxRun := 1
	currentRun := 1

	for i := 1; i < len(kmer); i++ {
		if kmer[i] == kmer[i-1] {
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

	// Penalty increases with run length
	if maxRun >= 3 {
		return float64(maxRun-2) * 2.0
	}

	return 0.0
}

// GetModelStats returns statistics about the trained model
func (m *HMMKmerModel) GetModelStats() map[string]interface{} {
	m.mu.RLock()
	defer m.mu.RUnlock()

	stats := map[string]interface{}{
		"k":           m.k,
		"num_states":  len(m.stateSpace),
		"trained":     m.trained,
		"total_kmers": m.totalKmers,
	}

	if m.trained {
		// Calculate entropy of transition matrix
		transitionEntropy := 0.0
		for i := 0; i < len(m.stateSpace); i++ {
			for j := 0; j < 4; j++ {
				prob := m.transitionMatrix[i][j]
				if prob > 0 {
					transitionEntropy -= prob * math.Log2(prob)
				}
			}
		}
		stats["transition_entropy"] = transitionEntropy

		// Calculate most probable states
		maxStartProb := 0.0
		maxStartState := ""
		for i, prob := range m.startProbabilities {
			if prob > maxStartProb {
				maxStartProb = prob
				maxStartState = m.stateSpace[i]
			}
		}
		stats["most_probable_start_state"] = maxStartState
		stats["max_start_probability"] = maxStartProb
	}

	return stats
}

// K returns the k-mer size
func (m *HMMKmerModel) K() int {
	return m.k
}

// IsTrained returns whether the model has been trained
func (m *HMMKmerModel) IsTrained() bool {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.trained
}

// TotalKmers returns the total number of k-mers processed during training
func (m *HMMKmerModel) TotalKmers() uint64 {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.totalKmers
}

// StateSpace returns the state space of the model
func (m *HMMKmerModel) StateSpace() []string {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.stateSpace
}

// SetSeed sets the random seed for reproducible generation
func (m *HMMKmerModel) SetSeed(seed int64) {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.rng = rand.New(rand.NewPCG(uint64(seed), 0)) // #nosec G115 G404
}

// HMMGenerationResult represents the result of HMM generation
type HMMGenerationResult struct {
	Sequence        []byte
	Qualities       []byte
	KmerFrequencies map[string]int
	GenerationTime  time.Duration
	States          []int
}

// ValidateModel validates the trained model for consistency
func (m *HMMKmerModel) ValidateModel() error {
	m.mu.RLock()
	defer m.mu.RUnlock()

	if !m.trained {
		return errors.NewValidationError("hmm_kmer_model", "validate_model", "model is not trained")
	}

	// Check transition matrix normalization
	for i := 0; i < len(m.stateSpace); i++ {
		sum := 0.0
		for j := 0; j < 4; j++ {
			sum += m.transitionMatrix[i][j]
		}
		if math.Abs(sum-1.0) > 1e-6 {
			return errors.NewValidationError("hmm_kmer_model", "validate_model", fmt.Sprintf("transition probabilities for state %d do not sum to 1.0: %f", i, sum))
		}
	}

	// Check emission matrix normalization
	for i := 0; i < len(m.stateSpace); i++ {
		sum := 0.0
		for j := 0; j < 42; j++ {
			sum += m.emissionMatrix[i][j]
		}
		if math.Abs(sum-1.0) > 1e-6 {
			return errors.NewValidationError("hmm_kmer_model", "validate_model", fmt.Sprintf("emission probabilities for state %d do not sum to 1.0: %f", i, sum))
		}
	}

	// Check start probabilities normalization
	sum := 0.0
	for _, prob := range m.startProbabilities {
		sum += prob
	}
	if math.Abs(sum-1.0) > 1e-6 {
		return errors.NewValidationError("hmm_kmer_model", "validate_model", fmt.Sprintf("start probabilities do not sum to 1.0: %f", sum))
	}

	return nil
}
