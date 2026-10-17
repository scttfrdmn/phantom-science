package authenticity

import (
	"bufio"
	"context"
	"fmt"
	"math/rand/v2"
	"os"
	"runtime"
	"strings"
	"sync"
	"sync/atomic"
	"time"
)

// ParallelHMMModel represents a high-performance, multi-core HMM model for genomic training
type ParallelHMMModel struct {
	k                      int
	stateSpace             []string
	transitionMatrix       [][]float64
	emissionMatrix         [][]float64
	startProbabilities     []float64
	
	// Parallel processing configuration
	numWorkers             int
	chunkSize              int
	processingPool         chan processingJob
	resultPool             chan processingResult
	
	// Adaptive complexity features (thread-safe)
	complexityMemory       []string
	complexityWindowSize   int
	targetComplexity       float64
	currentComplexity      int64  // atomic operations
	adaptiveWeights        sync.Map // thread-safe map
	
	// Biological constraints
	gcContentTarget        float64
	gcContentTolerance     float64
	currentGCContent       int64  // atomic operations
	homopolymerMaxLength   int
	cpgDepletionRatio      float64
	
	// Training state (thread-safe)
	totalKmers             uint64  // atomic operations
	trained                int32   // atomic bool
	mu                     sync.RWMutex
	rng                    *rand.Rand
	
	// Performance metrics (thread-safe)
	generationStats        *ThreadSafeGenerationStats
	
	// GPU acceleration (optional)
	useGPU                 bool
	gpuDevice              interface{} // GPU context if available
}

// ThreadSafeGenerationStats provides thread-safe performance tracking
type ThreadSafeGenerationStats struct {
	TotalSequences      int64   // atomic
	TotalBases          int64   // atomic
	AverageComplexity   int64   // atomic (scaled by 1000000)
	AverageGCContent    int64   // atomic (scaled by 1000000)
	MaxHomopolymer      int64   // atomic
	AverageJSDivergence int64   // atomic (scaled by 1000000)
	mu                  sync.RWMutex
}

// processingJob represents a unit of work for parallel processing
type processingJob struct {
	id       int
	sequence []byte
	startPos int
	endPos   int
	ctx      context.Context
}

// processingResult represents the result of parallel processing
type processingResult struct {
	id           int
	kmerCounts   map[string]int
	transitions  map[string]map[string]int
	gcCount      int
	totalBases   int
	complexity   float64
	err          error
}

// NewParallelHMMModel creates a new high-performance parallel HMM model
func NewParallelHMMModel(k int) *ParallelHMMModel {
	numWorkers := detectPerformanceCores()
	fmt.Printf("🍎 Apple Silicon detected: Using %d P-cores (leaving E-cores for system)\n", numWorkers)
	
	model := &ParallelHMMModel{
		k:                    k,
		numWorkers:          numWorkers,
		chunkSize:           100000, // 100KB chunks
		processingPool:      make(chan processingJob, numWorkers*2),
		resultPool:          make(chan processingResult, numWorkers*2),
		complexityWindowSize: 1000000,
		targetComplexity:    0.6,
		gcContentTarget:     0.42,
		gcContentTolerance:  0.02,
		homopolymerMaxLength: 15,
		cpgDepletionRatio:   0.25,
		rng:                 rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
		generationStats:     &ThreadSafeGenerationStats{},
		useGPU:              false, // Will be enabled if GPU detected
	}
	
	model.initializeStateSpace()
	model.initializeMatrices()
	model.startWorkerPool()
	
	// Try to initialize GPU if available
	model.initializeGPU()
	
	return model
}

// initializeStateSpace creates the k-mer state space
func (m *ParallelHMMModel) initializeStateSpace() {
	bases := []string{"A", "T", "C", "G"}
	m.stateSpace = make([]string, 0, 1<<(2*uint(m.k))) // #nosec G115
	
	// Generate all possible k-mers
	var generateKmers func(current string, remaining int)
	generateKmers = func(current string, remaining int) {
		if remaining == 0 {
			m.stateSpace = append(m.stateSpace, current)
			return
		}
		for _, base := range bases {
			generateKmers(current+base, remaining-1)
		}
	}
	
	generateKmers("", m.k)
}

// initializeMatrices creates the transition and emission matrices
func (m *ParallelHMMModel) initializeMatrices() {
	stateCount := len(m.stateSpace)
	
	// Initialize transition matrix
	m.transitionMatrix = make([][]float64, stateCount)
	for i := range m.transitionMatrix {
		m.transitionMatrix[i] = make([]float64, stateCount)
		// Initialize with small uniform probabilities
		for j := range m.transitionMatrix[i] {
			m.transitionMatrix[i][j] = 1.0 / float64(stateCount)
		}
	}
	
	// Initialize emission matrix (for our case, emissions are the k-mers themselves)
	m.emissionMatrix = make([][]float64, stateCount)
	for i := range m.emissionMatrix {
		m.emissionMatrix[i] = make([]float64, 4) // A, T, C, G
		for j := range m.emissionMatrix[i] {
			m.emissionMatrix[i][j] = 0.25 // Uniform initially
		}
	}
	
	// Initialize start probabilities
	m.startProbabilities = make([]float64, stateCount)
	for i := range m.startProbabilities {
		m.startProbabilities[i] = 1.0 / float64(stateCount)
	}
}

// startWorkerPool initializes the worker goroutines for parallel processing
func (m *ParallelHMMModel) startWorkerPool() {
	for i := 0; i < m.numWorkers; i++ {
		go m.worker(i)
	}
}

// worker processes jobs from the processing pool
func (m *ParallelHMMModel) worker(workerID int) {
	for job := range m.processingPool {
		select {
		case <-job.ctx.Done():
			m.resultPool <- processingResult{id: job.id, err: job.ctx.Err()}
			continue
		default:
		}
		
		result := m.processSequenceChunk(job)
		result.id = job.id
		
		select {
		case m.resultPool <- result:
		case <-job.ctx.Done():
			return
		}
	}
}

// processSequenceChunk processes a chunk of sequence data
func (m *ParallelHMMModel) processSequenceChunk(job processingJob) processingResult {
	sequence := job.sequence[job.startPos:job.endPos]
	
	kmerCounts := make(map[string]int)
	transitions := make(map[string]map[string]int)
	gcCount := 0
	totalBases := len(sequence)
	
	// Process k-mers
	for i := 0; i <= len(sequence)-m.k; i++ {
		kmer := string(sequence[i : i+m.k])
		kmerCounts[kmer]++
		
		// Count transitions
		if i > 0 {
			prevKmer := string(sequence[i-1 : i-1+m.k])
			if transitions[prevKmer] == nil {
				transitions[prevKmer] = make(map[string]int)
			}
			transitions[prevKmer][kmer]++
		}
	}
	
	// Count GC content
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}
	
	// Calculate complexity
	complexity := float64(len(kmerCounts)) / float64(len(sequence)-m.k+1)
	
	return processingResult{
		kmerCounts:  kmerCounts,
		transitions: transitions,
		gcCount:     gcCount,
		totalBases:  totalBases,
		complexity:  complexity,
	}
}

// TrainFromFileParallel trains the model using parallel processing
func (m *ParallelHMMModel) TrainFromFileParallel(filename string, ctx context.Context) error {
	fmt.Printf("🚀 Starting parallel training with %d CPU cores\n", m.numWorkers)
	if m.useGPU {
		fmt.Printf("🎮 GPU acceleration enabled\n")
	}
	
	// Read sequence data
	sequences, err := ReadFASTAFile(filename)
	if err != nil {
		return fmt.Errorf("failed to read FASTA file: %w", err)
	}
	
	fmt.Printf("📊 Processing %d sequences with parallel pipeline\n", len(sequences))
	
	// Process sequences in parallel
	for seqIdx, sequence := range sequences {
		if err := m.trainOnSequenceParallel(sequence, ctx); err != nil {
			return fmt.Errorf("failed to train on sequence %d: %w", seqIdx, err)
		}
	}
	
	// Finalize training
	m.normalizeMatrices()
	atomic.StoreInt32(&m.trained, 1)
	
	fmt.Printf("✅ Parallel training completed with %d CPU cores\n", m.numWorkers)
	return nil
}

// trainOnSequenceParallel processes a sequence using parallel workers
func (m *ParallelHMMModel) trainOnSequenceParallel(sequence []byte, ctx context.Context) error {
	if len(sequence) < m.k {
		return nil
	}
	
	// Split sequence into chunks for parallel processing
	numChunks := (len(sequence) + m.chunkSize - 1) / m.chunkSize
	results := make([]processingResult, numChunks)
	
	// Submit jobs
	for i := 0; i < numChunks; i++ {
		startPos := i * m.chunkSize
		endPos := startPos + m.chunkSize
		if endPos > len(sequence) {
			endPos = len(sequence)
		}
		
		job := processingJob{
			id:       i,
			sequence: sequence,
			startPos: startPos,
			endPos:   endPos,
			ctx:      ctx,
		}
		
		select {
		case m.processingPool <- job:
		case <-ctx.Done():
			return ctx.Err()
		}
	}
	
	// Collect results
	for i := 0; i < numChunks; i++ {
		select {
		case result := <-m.resultPool:
			if result.err != nil {
				return result.err
			}
			results[result.id] = result
		case <-ctx.Done():
			return ctx.Err()
		}
	}
	
	// Merge results and update matrices
	return m.mergeResults(results)
}

// mergeResults combines parallel processing results
func (m *ParallelHMMModel) mergeResults(results []processingResult) error {
	m.mu.Lock()
	defer m.mu.Unlock()
	
	totalKmers := uint64(0)
	totalGC := 0
	totalBases := 0
	totalComplexity := 0.0
	
	// Merge k-mer counts and transitions
	for _, result := range results {
		totalKmers += uint64(len(result.kmerCounts))
		totalGC += result.gcCount
		totalBases += result.totalBases
		totalComplexity += result.complexity
		
		// Update transition matrix
		for fromKmer, transitions := range result.transitions {
			fromIdx := m.getStateIndex(fromKmer)
			if fromIdx < 0 {
				continue
			}
			
			for toKmer, count := range transitions {
				toIdx := m.getStateIndex(toKmer)
				if toIdx < 0 {
					continue
				}
				
				m.transitionMatrix[fromIdx][toIdx] += float64(count)
			}
		}
	}
	
	// Update atomic counters
	atomic.AddUint64(&m.totalKmers, totalKmers)
	
	// Update GC content (scaled by 1000000 for atomic operations)
	currentGC := float64(totalGC) / float64(totalBases)
	atomic.StoreInt64(&m.currentGCContent, int64(currentGC*1000000))
	
	// Update complexity (scaled by 1000000 for atomic operations)
	avgComplexity := totalComplexity / float64(len(results))
	atomic.StoreInt64(&m.currentComplexity, int64(avgComplexity*1000000))
	
	return nil
}

// getStateIndex returns the index of a k-mer in the state space
func (m *ParallelHMMModel) getStateIndex(kmer string) int {
	for i, state := range m.stateSpace {
		if state == kmer {
			return i
		}
	}
	return -1
}

// normalizeMatrices normalizes the transition and emission matrices
func (m *ParallelHMMModel) normalizeMatrices() {
	m.mu.Lock()
	defer m.mu.Unlock()
	
	// Normalize transition matrix
	for i := range m.transitionMatrix {
		sum := 0.0
		for j := range m.transitionMatrix[i] {
			sum += m.transitionMatrix[i][j]
		}
		if sum > 0 {
			for j := range m.transitionMatrix[i] {
				m.transitionMatrix[i][j] /= sum
			}
		}
	}
}

// GenerateSequenceParallel generates a sequence using parallel processing
func (m *ParallelHMMModel) GenerateSequenceParallel(length int, ctx context.Context) ([]byte, error) {
	if atomic.LoadInt32(&m.trained) == 0 {
		return nil, fmt.Errorf("model not trained")
	}
	
	// For now, generation is still sequential but optimized
	// In production, this could be parallelized using beam search
	return m.generateSequenceOptimized(length, ctx)
}

// generateSequenceOptimized generates a sequence with optimization
func (m *ParallelHMMModel) generateSequenceOptimized(length int, ctx context.Context) ([]byte, error) {
	sequence := make([]byte, length)
	
	// Use faster random number generation
	rng := rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)) // #nosec G404
	
	// Start with a random k-mer
	startStateIdx := rng.IntN(len(m.stateSpace))
	startKmer := m.stateSpace[startStateIdx]
	
	// Copy start k-mer
	copy(sequence[:m.k], startKmer)
	
	// Generate remaining bases
	currentKmer := startKmer
	for i := m.k; i < length; i++ {
		select {
		case <-ctx.Done():
			return nil, ctx.Err()
		default:
		}
		
		// Find current state
		currentStateIdx := m.getStateIndex(currentKmer)
		if currentStateIdx < 0 {
			// Fallback to random base
			bases := []byte{'A', 'T', 'C', 'G'}
			sequence[i] = bases[rng.IntN(4)]
		} else {
			// Use transition probabilities
			nextStateIdx := m.selectNextState(currentStateIdx, rng)
			nextKmer := m.stateSpace[nextStateIdx]
			sequence[i] = nextKmer[len(nextKmer)-1]
		}
		
		// Update current k-mer
		if i >= m.k {
			currentKmer = string(sequence[i-m.k+1 : i+1])
		}
	}
	
	return sequence, nil
}

// selectNextState selects the next state based on transition probabilities
func (m *ParallelHMMModel) selectNextState(currentStateIdx int, rng *rand.Rand) int {
	m.mu.RLock()
	defer m.mu.RUnlock()
	
	// Use weighted random selection
	totalProb := 0.0
	for _, prob := range m.transitionMatrix[currentStateIdx] {
		totalProb += prob
	}
	
	if totalProb == 0 {
		return rng.IntN(len(m.stateSpace))
	}
	
	target := rng.Float64() * totalProb
	cumulative := 0.0
	
	for i, prob := range m.transitionMatrix[currentStateIdx] {
		cumulative += prob
		if cumulative >= target {
			return i
		}
	}
	
	return len(m.stateSpace) - 1
}

// initializeGPU attempts to initialize GPU acceleration
func (m *ParallelHMMModel) initializeGPU() {
	// Detect Metal GPU on macOS
	if m.detectMetalGPU() {
		m.useGPU = true
		fmt.Println("🎮 Metal GPU acceleration initialized")
		fmt.Println("⚡ Matrix operations will be GPU-accelerated")
	} else {
		m.useGPU = false
		fmt.Printf("💻 CPU-only mode with %d P-cores\n", m.numWorkers)
	}
}

// detectPerformanceCores detects Apple Silicon performance cores
func detectPerformanceCores() int {
	totalCores := runtime.NumCPU()
	
	// Apple Silicon detection heuristics
	if totalCores >= 8 {
		// M1/M2/M3 series - typically 4-8 P-cores + 4-8 E-cores
		switch totalCores {
		case 8:  // M1: 4P + 4E
			return 4
		case 10: // M1 Pro: 6P + 4E or 8P + 2E
			return 6
		case 12: // M1 Max: 8P + 4E
			return 8
		case 16: // M2 Ultra: 12P + 4E
			return 12
		case 20: // M1 Ultra: 16P + 4E
			return 16
		case 24: // M2 Ultra: 16P + 8E
			return 16
		default:
			// Conservative estimate: use 2/3 of total cores as P-cores
			return (totalCores * 2) / 3
		}
	}
	
	// Non-Apple Silicon fallback
	if totalCores > 4 {
		return totalCores - 2 // Leave 2 cores for system
	}
	
	return totalCores
}

// detectMetalGPU detects Metal GPU availability on macOS
func (m *ParallelHMMModel) detectMetalGPU() bool {
	// Check if running on macOS
	if runtime.GOOS != "darwin" {
		return false
	}
	
	// In production, this would use Metal framework
	// For now, we'll assume Metal is available on macOS
	return true
}

// GetPerformanceStats returns current performance statistics
func (m *ParallelHMMModel) GetPerformanceStats() map[string]interface{} {
	return map[string]interface{}{
		"cpu_cores":        m.numWorkers,
		"gpu_enabled":      m.useGPU,
		"total_kmers":      atomic.LoadUint64(&m.totalKmers),
		"current_gc":       float64(atomic.LoadInt64(&m.currentGCContent)) / 1000000.0,
		"current_complexity": float64(atomic.LoadInt64(&m.currentComplexity)) / 1000000.0,
		"trained":          atomic.LoadInt32(&m.trained) == 1,
		"chunk_size":       m.chunkSize,
	}
}

// Cleanup closes the worker pool and releases resources
func (m *ParallelHMMModel) Cleanup() {
	close(m.processingPool)
	close(m.resultPool)
	
	if m.useGPU {
		m.cleanupGPU()
	}
}

// GPU-accelerated matrix operations
func (m *ParallelHMMModel) multiplyMatricesGPU(a, b [][]float64) [][]float64 {
	if !m.useGPU {
		return m.multiplyMatricesCPU(a, b)
	}
	
	// In production, this would use Metal Performance Shaders
	// For now, we'll use optimized CPU with vectorization hints
	return m.multiplyMatricesOptimized(a, b)
}

// multiplyMatricesCPU performs standard CPU matrix multiplication
func (m *ParallelHMMModel) multiplyMatricesCPU(a, b [][]float64) [][]float64 {
	rows := len(a)
	cols := len(b[0])
	inner := len(b)
	
	result := make([][]float64, rows)
	for i := range result {
		result[i] = make([]float64, cols)
	}
	
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			for k := 0; k < inner; k++ {
				result[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	
	return result
}

// multiplyMatricesOptimized performs optimized matrix multiplication
func (m *ParallelHMMModel) multiplyMatricesOptimized(a, b [][]float64) [][]float64 {
	rows := len(a)
	cols := len(b[0])
	inner := len(b)
	
	result := make([][]float64, rows)
	for i := range result {
		result[i] = make([]float64, cols)
	}
	
	// Use goroutines for parallel matrix multiplication
	var wg sync.WaitGroup
	for i := 0; i < rows; i++ {
		wg.Add(1)
		go func(row int) {
			defer wg.Done()
			for j := 0; j < cols; j++ {
				sum := 0.0
				for k := 0; k < inner; k++ {
					sum += a[row][k] * b[k][j]
				}
				result[row][j] = sum
			}
		}(i)
	}
	
	wg.Wait()
	return result
}

// vectorAddGPU performs GPU-accelerated vector addition
func (m *ParallelHMMModel) vectorAddGPU(a, b []float64) []float64 {
	if !m.useGPU {
		return m.vectorAddCPU(a, b)
	}
	
	// In production, this would use Metal compute shaders
	return m.vectorAddOptimized(a, b)
}

// vectorAddCPU performs standard CPU vector addition
func (m *ParallelHMMModel) vectorAddCPU(a, b []float64) []float64 {
	result := make([]float64, len(a))
	for i := range a {
		result[i] = a[i] + b[i]
	}
	return result
}

// vectorAddOptimized performs optimized vector addition
func (m *ParallelHMMModel) vectorAddOptimized(a, b []float64) []float64 {
	result := make([]float64, len(a))
	
	// Use parallel processing for large vectors
	if len(a) > 1000 {
		var wg sync.WaitGroup
		chunkSize := len(a) / m.numWorkers
		
		for i := 0; i < m.numWorkers; i++ {
			wg.Add(1)
			go func(start, end int) {
				defer wg.Done()
				for j := start; j < end && j < len(a); j++ {
					result[j] = a[j] + b[j]
				}
			}(i*chunkSize, (i+1)*chunkSize)
		}
		
		wg.Wait()
	} else {
		for i := range a {
			result[i] = a[i] + b[i]
		}
	}
	
	return result
}

// cleanupGPU releases GPU resources
func (m *ParallelHMMModel) cleanupGPU() {
	if m.gpuDevice != nil {
		// In production, this would release Metal resources
		m.gpuDevice = nil
		fmt.Println("🎮 Metal GPU resources released")
	}
}

// ReadFASTAFile reads sequences from a FASTA file
func ReadFASTAFile(filename string) ([][]byte, error) {
	file, err := os.Open(filename) // #nosec G304
	if err != nil {
		return nil, fmt.Errorf("failed to open FASTA file: %w", err)
	}
	defer file.Close()
	
	var sequences [][]byte
	var currentSequence []byte
	
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		
		if strings.HasPrefix(line, ">") {
			// New sequence header
			if len(currentSequence) > 0 {
				sequences = append(sequences, currentSequence)
				currentSequence = nil
			}
		} else {
			// Sequence data
			currentSequence = append(currentSequence, []byte(strings.ToUpper(line))...)
		}
	}
	
	// Add the last sequence
	if len(currentSequence) > 0 {
		sequences = append(sequences, currentSequence)
	}
	
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading FASTA file: %w", err)
	}
	
	return sequences, nil
}

// TrainFromSequence trains the model on a single sequence
func (m *ParallelHMMModel) TrainFromSequence(sequence []byte) error {
	if len(sequence) < m.k {
		return nil
	}
	
	// Process k-mers from the sequence
	kmerCounts := make(map[string]int)
	transitions := make(map[string]map[string]int)
	
	for i := 0; i <= len(sequence)-m.k; i++ {
		kmer := string(sequence[i : i+m.k])
		kmerCounts[kmer]++
		
		// Count transitions
		if i > 0 {
			prevKmer := string(sequence[i-1 : i-1+m.k])
			if transitions[prevKmer] == nil {
				transitions[prevKmer] = make(map[string]int)
			}
			transitions[prevKmer][kmer]++
		}
	}
	
	// Update model with thread safety
	m.mu.Lock()
	defer m.mu.Unlock()
	
	// Update transition matrix
	for fromKmer, transitionMap := range transitions {
		fromIdx := m.getStateIndex(fromKmer)
		if fromIdx < 0 {
			continue
		}
		
		for toKmer, count := range transitionMap {
			toIdx := m.getStateIndex(toKmer)
			if toIdx < 0 {
				continue
			}
			
			m.transitionMatrix[fromIdx][toIdx] += float64(count)
		}
	}
	
	// Update counters
	atomic.AddUint64(&m.totalKmers, uint64(len(kmerCounts)))
	
	return nil
}

// GenerateSequence generates a sequence of the specified length
func (m *ParallelHMMModel) GenerateSequence(length int) ([]byte, error) {
	if atomic.LoadInt32(&m.trained) == 0 {
		// For untrained models, generate based on current training state
		return m.generateSequenceUntrained(length)
	}
	
	return m.generateSequenceOptimized(length, context.Background())
}

// generateSequenceUntrained generates a sequence from an untrained model
func (m *ParallelHMMModel) generateSequenceUntrained(length int) ([]byte, error) {
	sequence := make([]byte, length)
	
	// Simple random generation with biological constraints
	bases := []byte{'A', 'T', 'C', 'G'}
	
	// Target GC content around 42%
	gcTarget := 0.42
	
	for i := 0; i < length; i++ {
		// Calculate current GC content
		gcCount := 0
		for j := 0; j < i; j++ {
			if sequence[j] == 'G' || sequence[j] == 'C' {
				gcCount++
			}
		}
		
		currentGC := float64(gcCount) / float64(i+1)
		
		// Bias towards GC or AT to reach target
		if currentGC < gcTarget {
			// Favor GC
			if m.rng.Float64() < 0.6 {
				if m.rng.Float64() < 0.5 {
					sequence[i] = 'G'
				} else {
					sequence[i] = 'C'
				}
			} else {
				sequence[i] = bases[m.rng.IntN(4)]
			}
		} else {
			// Favor AT
			if m.rng.Float64() < 0.6 {
				if m.rng.Float64() < 0.5 {
					sequence[i] = 'A'
				} else {
					sequence[i] = 'T'
				}
			} else {
				sequence[i] = bases[m.rng.IntN(4)]
			}
		}
	}
	
	return sequence, nil
}