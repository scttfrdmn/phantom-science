package authenticity

import "math/rand/v2"

// Shared estimation and sampling for the order-(k-1) Markov chains
// (HMMKmerModel and AdaptiveHMMModel).
//
// A state is a (k-1)-mer over the alphabet A,T,G,C (in that order); its
// index is the base-4 number of its bases, first base most significant, which
// is also its position in the model's state space. Training counts, for every
// k-mer x[i:i+k] of the input made of A/C/G/T only (case-insensitive), one
// observation of base x[i+k-1] following context x[i:i+k-1]. Estimation
// adds a pseudocount α = PseudoCount to every cell (Laplace smoothing):
//
//	P(b | s) = (c(s,b) + α) / (n(s) + 4α)        n(s) = Σ_b c(s,b)
//	π(s)     = (n(s) + α)   / (N + S·α)          N = Σ_s n(s), S = 4^(k-1)
//
// so every row is a proper distribution, contexts never seen in training get
// the uniform row, and the start distribution π is the smoothed frequency of
// each (k-1)-mer as a context.

// PseudoCount is the Laplace pseudocount α added to every transition and start
// count when a model's probabilities are estimated.
const PseudoCount = 1.0

// DefaultSeed seeds a model's random source until SetSeed is called, so an
// unseeded model (new or loaded) is still reproducible; models never seed
// from the clock. The genomics engine reseeds every model per output stream
// with SetSeed.
const DefaultSeed int64 = 0

// newModelRNG returns the random source a model uses for seed, as
// PCG(uint64(seed), 0).
func newModelRNG(seed int64) *rand.Rand {
	return rand.New(rand.NewPCG(uint64(seed), 0)) // #nosec G115 G404 -- simulation, not crypto
}

// basesATGC is the base order of state indices and transition columns.
var basesATGC = [4]byte{'A', 'T', 'G', 'C'}

// baseIndex returns the column of base b (A=0, T=1, G=2, C=3; lower case
// accepted), or -1 for any other byte.
func baseIndex(b byte) int {
	switch b {
	case 'A', 'a':
		return 0
	case 'T', 't':
		return 1
	case 'G', 'g':
		return 2
	case 'C', 'c':
		return 3
	}
	return -1
}

// contextIndex returns the state index of the (k-1)-mer ctx, or -1 if it
// contains a byte other than A/C/G/T.
func contextIndex(ctx []byte) int {
	idx := 0
	for _, b := range ctx {
		bi := baseIndex(b)
		if bi < 0 {
			return -1
		}
		idx = idx*4 + bi
	}
	return idx
}

// countKmers adds the transitions of every valid k-mer in seq to counts
// (counts[state][base]) and returns how many k-mers were counted.
func countKmers(seq []byte, k int, counts [][]float64) uint64 {
	var n uint64
	for i := 0; i+k <= len(seq); i++ {
		s := contextIndex(seq[i : i+k-1])
		b := baseIndex(seq[i+k-1])
		if s < 0 || b < 0 {
			continue
		}
		counts[s][b]++
		n++
	}
	return n
}

// estimateChain fills trans (rows of 4) and start from counts using
// PseudoCount, as documented above.
func estimateChain(counts [][]float64, trans [][]float64, start []float64) {
	total := 0.0
	for s, row := range counts {
		n := row[0] + row[1] + row[2] + row[3]
		total += n
		den := n + 4*PseudoCount
		for b := 0; b < 4; b++ {
			trans[s][b] = (row[b] + PseudoCount) / den
		}
		start[s] = n // normalised below
	}
	den := total + float64(float64(len(counts))*PseudoCount)
	for s := range start {
		start[s] = (start[s] + PseudoCount) / den
	}
}

// sampleIndex draws an index from the weights w (which need not sum to 1).
// Only indices with positive weight can be returned; if no weight is
// positive the draw is uniform.
func sampleIndex(rng *rand.Rand, w []float64) int {
	total := 0.0
	last := -1
	for i, p := range w {
		if p > 0 {
			total += p
			last = i
		}
	}
	if last < 0 {
		return rng.IntN(len(w))
	}
	r := rng.Float64() * total
	acc := 0.0
	for i, p := range w {
		if p <= 0 {
			continue
		}
		acc += p
		if r < acc {
			return i
		}
	}
	return last // round-off: r landed at the very top
}

// newCountMatrix returns a zeroed states×cols matrix.
func newCountMatrix(states, cols int) [][]float64 {
	m := make([][]float64, states)
	for i := range m {
		m[i] = make([]float64, cols)
	}
	return m
}
