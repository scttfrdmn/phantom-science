package authenticity

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/errors"
)

// HMMKmerModel represents a Hidden Markov Model for k-mer based sequence generation
type HMMKmerModel struct {
	k                  int
	stateSpace         []string              // All possible (k-1)-mers
	transitionMatrix   [][]float64           // State -> base probabilities
	emissionMatrix     [][]float64           // State -> quality probabilities
	startProbabilities []float64             // Initial state probabilities
	totalKmers         uint64                // Total k-mers processed during training
	trained            bool                  // Whether model has been trained
	mu                 sync.RWMutex          // Thread safety
	rng                *rand.Rand            // Random number generator
}

// KmerState represents a state in the HMM (k-1 mer)
type KmerState struct {
	Index    int
	Kmer     string
	Transitions map[byte]float64  // Base -> probability
	Emissions   map[int]float64   // Quality -> probability
}

// GenerationResult represents the result of sequence generation
type GenerationResult struct {
	Sequence       []byte
	Qualities      []byte
	KmerFrequencies map[string]int
	GenerationTime time.Duration
	States         []int  // State sequence for debugging
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

// generateStateSpace generates all possible (k-1)-mers
func (m *HMMKmerModel) generateStateSpace() []string {
	bases := []string{"A", "T", "G", "C"}
	states := []string{""}

	for i := 0; i < m.k-1; i++ {
		var newStates []string
		for _, state := range states {
			for _, base := range bases {
				newStates = append(newStates, state+base)
			}
		}
		states = newStates
	}

	return states
}

// getStateIndex returns the index of a state in the state space
func (m *HMMKmerModel) getStateIndex(state string) int {
	for i, s := range m.stateSpace {
		if s == state {
			return i
		}
	}
	return -1 // Should never happen with valid input
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

// TrainFromReferenceGenome trains the HMM model using reference genome k-mer frequencies
func (m *HMMKmerModel) TrainFromReferenceGenome(refGenome GenomeInterface) error {
	m.mu.Lock()
	defer m.mu.Unlock()

	// Initialize counters
	transitionCounts := make([][]int, len(m.stateSpace))
	emissionCounts := make([][]int, len(m.stateSpace))
	startCounts := make([]int, len(m.stateSpace))

	for i := 0; i < len(m.stateSpace); i++ {
		transitionCounts[i] = make([]int, 4)
		emissionCounts[i] = make([]int, 42)
	}

	// Train on each chromosome
	for _, chrName := range refGenome.GetChromosomeNames() {
		chr, err := refGenome.LoadChromosome(chrName)
		if err != nil {
			continue // Skip problematic chromosomes
		}

		err = m.trainOnChromosome(chr, transitionCounts, emissionCounts, startCounts)
		if err != nil {
			return errors.WrapError(err, errors.BiologicalModelError, "hmm_kmer_model", "train_from_reference", fmt.Sprintf("failed to train on chromosome %s", chrName))
		}
	}

	// Convert counts to probabilities
	m.normalizeProbabilities(transitionCounts, emissionCounts, startCounts)
	m.trained = true

	return nil
}

// trainOnChromosome trains the model on a single chromosome
func (m *HMMKmerModel) trainOnChromosome(chr ChromosomeInterface, 
	transitionCounts, emissionCounts [][]int, startCounts []int) error {
	
	sequence := chr.GetSequence()
	if len(sequence) < m.k {
		return nil // Skip very short sequences
	}

	// Process sequence with sliding window
	for i := 0; i <= len(sequence)-m.k; i++ {
		kmer := string(sequence[i : i+m.k])
		
		// Skip k-mers with N bases
		if containsN(kmer) {
			continue
		}

		// Extract state (k-1 mer) and observation (next base)
		state := kmer[:m.k-1]
		nextBase := kmer[m.k-1]

		stateIndex := m.getStateIndex(state)
		if stateIndex == -1 {
			continue // Invalid state
		}

		baseIndex := m.getBaseIndex(nextBase)
		if baseIndex == -1 {
			continue // Invalid base
		}

		// Update transition counts
		transitionCounts[stateIndex][baseIndex]++

		// Update start counts for first k-mer
		if i == 0 {
			startCounts[stateIndex]++
		}

		// Simulate quality scores based on position and context
		// This is a simplified model - in practice, we'd use real quality data
		quality := m.simulateQuality(i, kmer)
		qualityIndex := int(quality)
		if qualityIndex >= 0 && qualityIndex < 42 {
			emissionCounts[stateIndex][qualityIndex]++
		}

		m.totalKmers++
	}

	return nil
}

// normalizeProbabilities converts counts to probabilities with smoothing
func (m *HMMKmerModel) normalizeProbabilities(transitionCounts, emissionCounts [][]int, startCounts []int) {
	const smoothingFactor = 0.1 // Laplace smoothing

	// Normalize transition probabilities
	for i := 0; i < len(m.stateSpace); i++ {
		total := 0.0
		for j := 0; j < 4; j++ {
			total += float64(transitionCounts[i][j]) + smoothingFactor
		}
		
		if total > 0 {
			for j := 0; j < 4; j++ {
				m.transitionMatrix[i][j] = (float64(transitionCounts[i][j]) + smoothingFactor) / total
			}
		} else {
			// Uniform distribution if no data
			for j := 0; j < 4; j++ {
				m.transitionMatrix[i][j] = 0.25
			}
		}
	}

	// Normalize emission probabilities
	for i := 0; i < len(m.stateSpace); i++ {
		total := 0.0
		for j := 0; j < 42; j++ {
			total += float64(emissionCounts[i][j]) + smoothingFactor
		}
		
		if total > 0 {
			for j := 0; j < 42; j++ {
				m.emissionMatrix[i][j] = (float64(emissionCounts[i][j]) + smoothingFactor) / total
			}
		} else {
			// Normal distribution centered at quality 20
			for j := 0; j < 42; j++ {
				m.emissionMatrix[i][j] = m.normalDistribution(float64(j), 20.0, 5.0)
			}
		}
	}

	// Normalize start probabilities
	startTotal := 0.0
	for i := 0; i < len(m.stateSpace); i++ {
		startTotal += float64(startCounts[i]) + smoothingFactor
	}
	
	if startTotal > 0 {
		for i := 0; i < len(m.stateSpace); i++ {
			m.startProbabilities[i] = (float64(startCounts[i]) + smoothingFactor) / startTotal
		}
	} else {
		// Uniform distribution
		prob := 1.0 / float64(len(m.stateSpace))
		for i := 0; i < len(m.stateSpace); i++ {
			m.startProbabilities[i] = prob
		}
	}
}

// GenerateSequence generates a biologically authentic sequence using the trained HMM
func (m *HMMKmerModel) GenerateSequence(length int) (*GenerationResult, error) {
	m.mu.RLock()
	defer m.mu.RUnlock()

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

// chooseInitialState chooses the initial state based on start probabilities
func (m *HMMKmerModel) chooseInitialState() int {
	r := m.rng.Float64()
	cumulative := 0.0
	
	for i, prob := range m.startProbabilities {
		cumulative += prob
		if r <= cumulative {
			return i
		}
	}
	
	// Fallback to last state
	return len(m.stateSpace) - 1
}

// sampleNextBase samples the next base from transition probabilities
func (m *HMMKmerModel) sampleNextBase(stateIndex int) byte {
	r := m.rng.Float64()
	cumulative := 0.0
	bases := []byte{'A', 'T', 'G', 'C'}
	
	for i, prob := range m.transitionMatrix[stateIndex] {
		cumulative += prob
		if r <= cumulative {
			return bases[i]
		}
	}
	
	// Fallback to 'A'
	return 'A'
}

// sampleQuality samples a quality score from emission probabilities
func (m *HMMKmerModel) sampleQuality(stateIndex int) int {
	r := m.rng.Float64()
	cumulative := 0.0
	
	for i, prob := range m.emissionMatrix[stateIndex] {
		cumulative += prob
		if r <= cumulative {
			return i
		}
	}
	
	// Fallback to quality 20
	return 20
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

// getBaseIndex returns the index of a base (A=0, T=1, G=2, C=3)
func (m *HMMKmerModel) getBaseIndex(base byte) int {
	switch base {
	case 'A':
		return 0
	case 'T':
		return 1
	case 'G':
		return 2
	case 'C':
		return 3
	default:
		return -1
	}
}

// simulateQuality simulates quality scores based on position and context
func (m *HMMKmerModel) simulateQuality(position int, kmer string) byte {
	// Simplified quality model - in practice, this would be more sophisticated
	baseQuality := 30.0
	
	// Quality degradation with position
	positionDecay := float64(position) * 0.01
	
	// GC bias effect
	gcContent := m.calculateGCContent(kmer)
	gcBias := math.Abs(gcContent-0.5) * 10.0
	
	// Homopolymer effect
	homopolymerPenalty := m.calculateHomopolymerPenalty(kmer)
	
	finalQuality := baseQuality - positionDecay - gcBias - homopolymerPenalty
	
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

// normalDistribution calculates normal distribution probability
func (m *HMMKmerModel) normalDistribution(x, mean, stddev float64) float64 {
	return math.Exp(-0.5*math.Pow((x-mean)/stddev, 2)) / (stddev * math.Sqrt(2*math.Pi))
}

// containsN checks if a string contains 'N' bases
func containsN(s string) bool {
	for _, c := range s {
		if c == 'N' {
			return true
		}
	}
	return false
}

// GetModelStats returns statistics about the trained model
func (m *HMMKmerModel) GetModelStats() map[string]interface{} {
	m.mu.RLock()
	defer m.mu.RUnlock()

	stats := map[string]interface{}{
		"k":             m.k,
		"num_states":    len(m.stateSpace),
		"trained":       m.trained,
		"total_kmers":   m.totalKmers,
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
	Sequence       []byte
	Qualities      []byte
	KmerFrequencies map[string]int
	GenerationTime time.Duration
	States         []int
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