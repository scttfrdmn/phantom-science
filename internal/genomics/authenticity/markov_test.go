//go:build !release

package authenticity

import (
	"fmt"
	"math"
	"math/rand/v2"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// randomDNA returns n bases drawn uniformly from ACGT with a fixed seed, with
// an N every nEvery bases (0 = none).
func randomDNA(seed uint64, n, nEvery int) []byte {
	r := rand.New(rand.NewPCG(seed, 1))
	out := make([]byte, n)
	for i := range out {
		out[i] = "ACGT"[r.IntN(4)]
		if nEvery > 0 && i%nEvery == nEvery-1 {
			out[i] = 'N'
		}
	}
	return out
}

// trainChunked trains m on seq in chunks of chunk bases, prefixing each chunk
// with the last k-1 bases seen (as cmd/train does).
func trainChunked(t *testing.T, m *AdaptiveHMMModel, seq []byte, chunk int) {
	t.Helper()
	var carry []byte
	for off := 0; off < len(seq); off += chunk {
		end := min(off+chunk, len(seq))
		buf := append(append([]byte(nil), carry...), seq[off:end]...)
		require.NoError(t, m.TrainFromSequence(buf))
		keep := min(m.K()-1, len(buf))
		carry = append(carry[:0], buf[len(buf)-keep:]...)
	}
}

// Hand-computed example, k=2 on "AACAG": transitions A→A, A→C, C→A, A→G, so
// n(A)=3, n(C)=1, n(G)=n(T)=0 and, with α=1,
// P(·|A) = (2,1,2,2)/7 over (A,T,G,C), P(·|C) = (2,1,1,1)/5, P(·|G) = P(·|T) = 1/4,
// π = (n(s)+1)/(4+4) = (4,1,1,2)/8 over states (A,T,G,C).
func TestAdaptiveHMM_TransitionProbabilitiesExact(t *testing.T) {
	m := NewAdaptiveHMMModel(2)
	require.NoError(t, m.TrainFromSequence([]byte("AACAG")))
	m.Finalize()

	A, T, G, C := 0, 1, 2, 3
	assert.Equal(t, []float64{2.0 / 7, 1.0 / 7, 2.0 / 7, 2.0 / 7}, m.transitionMatrix[A])
	assert.Equal(t, []float64{2.0 / 5, 1.0 / 5, 1.0 / 5, 1.0 / 5}, m.transitionMatrix[C])
	assert.Equal(t, []float64{0.25, 0.25, 0.25, 0.25}, m.transitionMatrix[G])
	assert.Equal(t, []float64{0.25, 0.25, 0.25, 0.25}, m.transitionMatrix[T])
	assert.Equal(t, []float64{4.0 / 8, 1.0 / 8, 1.0 / 8, 2.0 / 8}, m.startProbabilities)
	assert.Equal(t, uint64(4), m.TotalKmers())
}

// The context of the k-mer x[i:i+k] is x[i:i+k-1] (not x[i-1:i+k-2]): on
// "ACGT" with k=3 the observations are AC→G and CG→T only.
func TestAdaptiveHMM_ContextIsPrecedingKMinus1Bases(t *testing.T) {
	m := NewAdaptiveHMMModel(3)
	require.NoError(t, m.TrainFromSequence([]byte("ACGT")))
	ac, cg := contextIndex([]byte("AC")), contextIndex([]byte("CG"))
	assert.Equal(t, []float64{0, 0, 1, 0}, m.transitionCounts[ac], "AC→G")
	assert.Equal(t, []float64{0, 1, 0, 0}, m.transitionCounts[cg], "CG→T")
	total := 0.0
	for _, row := range m.transitionCounts {
		for _, c := range row {
			total += c
		}
	}
	assert.Equal(t, 2.0, total)
}

// On a random sequence every count equals a direct k-mer tally, and every
// probability is exactly (c+α)/(n+4α).
func TestAdaptiveHMM_ProbabilitiesMatchEmpiricalCounts(t *testing.T) {
	seq := randomDNA(3, 20000, 97)
	for _, k := range []int{2, 4, 6} {
		t.Run(fmt.Sprintf("k=%d", k), func(t *testing.T) {
			m := NewAdaptiveHMMModel(k)
			require.NoError(t, m.TrainFromSequence(seq))
			m.Finalize()

			tally := map[string]float64{}
			for i := 0; i+k <= len(seq); i++ {
				km := string(seq[i : i+k])
				if contextIndex([]byte(km)) >= 0 {
					tally[km]++
				}
			}
			for s, ctx := range m.stateSpace {
				n := 0.0
				for _, b := range basesATGC {
					n += tally[ctx+string(b)]
				}
				for j, b := range basesATGC {
					c := tally[ctx+string(b)]
					require.Equal(t, c, m.transitionCounts[s][j])
					require.Equal(t, (c+PseudoCount)/(n+4*PseudoCount), m.transitionMatrix[s][j])
				}
			}
		})
	}
}

// Every transition row and the start distribution sum to 1 (1e-12), for
// both models and every k.
func TestHMM_RowsSumToOne(t *testing.T) {
	seq := randomDNA(4, 5000, 0)
	genome := &MockReferenceGenome{chromosomes: map[string]*MockChromosome{
		"chr1": {Name: "chr1", Sequence: seq, Length: len(seq)},
	}}
	check := func(t *testing.T, trans [][]float64, start []float64) {
		for i, row := range trans {
			sum := 0.0
			for _, p := range row {
				require.Greater(t, p, 0.0)
				sum += p
			}
			require.InDelta(t, 1.0, sum, 1e-12, "row %d", i)
		}
		sum := 0.0
		for _, p := range start {
			require.Greater(t, p, 0.0)
			sum += p
		}
		require.InDelta(t, 1.0, sum, 1e-12, "start")
	}
	for k := 2; k <= 8; k++ {
		t.Run(fmt.Sprintf("k=%d", k), func(t *testing.T) {
			a := NewAdaptiveHMMModel(k)
			require.NoError(t, a.TrainFromSequence(seq))
			a.Finalize()
			check(t, a.transitionMatrix, a.startProbabilities)
			require.NoError(t, a.ValidateModel())

			km := NewHMMKmerModel(k)
			require.NoError(t, km.TrainFromReferenceGenome(genome))
			check(t, km.transitionMatrix, km.startProbabilities)
			require.NoError(t, km.ValidateModel())
		})
	}
}

// Training in 7-base chunks with a (k-1)-base carry gives exactly the
// one-shot counts and probabilities.
func TestAdaptiveHMM_ChunkedTrainingEqualsOneShot(t *testing.T) {
	seq := randomDNA(5, 3001, 61)
	for k := 2; k <= 8; k++ {
		t.Run(fmt.Sprintf("k=%d", k), func(t *testing.T) {
			one := NewAdaptiveHMMModel(k)
			require.NoError(t, one.TrainFromSequence(seq))
			chunked := NewAdaptiveHMMModel(k)
			trainChunked(t, chunked, seq, 7)

			assert.Equal(t, one.TotalKmers(), chunked.TotalKmers())
			assert.Equal(t, one.transitionCounts, chunked.transitionCounts)
			one.Finalize()
			chunked.Finalize()
			assert.Equal(t, one.transitionMatrix, chunked.transitionMatrix)
			assert.Equal(t, one.startProbabilities, chunked.startProbabilities)
		})
	}
}

// Finalize is idempotent, and training after Finalize re-finalizes from all
// the counts on the next generation.
func TestAdaptiveHMM_FinalizeIdempotentAndRefinalizes(t *testing.T) {
	m := NewAdaptiveHMMModel(3)
	require.NoError(t, m.TrainFromSequence([]byte("ACGTACGTAC")))
	m.Finalize()
	first := append([]float64(nil), m.startProbabilities...)
	m.Finalize()
	assert.Equal(t, first, m.startProbabilities)

	require.NoError(t, m.TrainFromSequence([]byte("GGGGGGGG")))
	assert.False(t, m.finalized)
	_, err := m.GenerateSequence(50)
	require.NoError(t, err)
	assert.True(t, m.finalized)

	ref := NewAdaptiveHMMModel(3)
	require.NoError(t, ref.TrainFromSequence([]byte("ACGTACGTAC")))
	require.NoError(t, ref.TrainFromSequence([]byte("GGGGGGGG")))
	ref.Finalize()
	assert.Equal(t, ref.transitionMatrix, m.transitionMatrix)
	assert.Equal(t, ref.startProbabilities, m.startProbabilities)
}

// TrainFromReferenceGenome counts each chromosome as TrainFromSequence does.
func TestAdaptiveHMM_TrainFromReferenceGenomeMatchesSequences(t *testing.T) {
	a, b := randomDNA(6, 900, 0), randomDNA(7, 1100, 50)
	genome := &MockReferenceGenome{chromosomes: map[string]*MockChromosome{
		"chr1": {Name: "chr1", Sequence: a, Length: len(a)},
		"chr2": {Name: "chr2", Sequence: b, Length: len(b)},
	}}
	g := NewAdaptiveHMMModel(4)
	require.NoError(t, g.TrainFromReferenceGenome(genome))
	s := NewAdaptiveHMMModel(4)
	require.NoError(t, s.TrainFromSequence(a))
	require.NoError(t, s.TrainFromSequence(b))
	s.Finalize()
	assert.Equal(t, s.transitionCounts, g.transitionCounts)
	assert.Equal(t, s.transitionMatrix, g.transitionMatrix)
}

// isPeriodic reports whether seq repeats with some period p <= maxP.
func isPeriodic(seq []byte, maxP int) bool {
	for p := 1; p <= maxP; p++ {
		periodic := true
		for i := 0; i+p < len(seq); i++ {
			if seq[i] != seq[i+p] {
				periodic = false
				break
			}
		}
		if periodic {
			return true
		}
	}
	return false
}

// baseCounts returns the A,T,G,C counts of seq.
func baseCounts(seq []byte) [4]int {
	var c [4]int
	for _, b := range seq {
		if i := baseIndex(b); i >= 0 {
			c[i]++
		}
	}
	return c
}

// A model trained on 20 bases (most contexts unseen) used to generate all C
// (unseen rows were all zero and sampling fell through to the last base). It
// must now generate a non-degenerate, aperiodic sequence using every base.
func TestAdaptiveHMM_TinyTrainingNotDegenerate(t *testing.T) {
	m := NewAdaptiveHMMModel(4)
	require.NoError(t, m.TrainFromSequence([]byte("ACGTAGCTAGGATCCATGCA")))
	m.SetSeed(1)
	res, err := m.GenerateSequence(10_000)
	require.NoError(t, err)
	c := baseCounts(res.Sequence)
	t.Logf("A,T,G,C counts: %v", c)
	for i, n := range c {
		assert.Greater(t, n, 1000, "base %c under-represented", basesATGC[i])
	}
	assert.False(t, isPeriodic(res.Sequence, 200))
}

// Even with a single context seen in training, the start state is sampled
// (first bases vary with the seed) and unseen rows are uniform, not 'C'.
func TestAdaptiveHMM_StartSampledAndUnseenRowsUniform(t *testing.T) {
	m := NewAdaptiveHMMModel(3)
	require.NoError(t, m.TrainFromSequence([]byte("AAAAAAAAAA")))
	starts := map[string]bool{}
	for seed := int64(0); seed < 200; seed++ {
		m.SetSeed(seed)
		res, err := m.GenerateSequence(3)
		require.NoError(t, err)
		starts[string(res.Sequence[:2])] = true
	}
	assert.Greater(t, len(starts), 8, "start (k-1)-mers seen: %v", starts)
	m.Finalize()
	assert.Equal(t, []float64{0.25, 0.25, 0.25, 0.25}, m.transitionMatrix[contextIndex([]byte("CC"))])
}

// stationaryBaseFreqs returns the stationary base frequencies of an order-(k-1)
// chain by power iteration.
func stationaryBaseFreqs(trans [][]float64) [4]float64 {
	n := len(trans)
	pi := make([]float64, n)
	for i := range pi {
		pi[i] = 1 / float64(n)
	}
	next := make([]float64, n)
	for iter := 0; iter < 2000; iter++ {
		clear(next)
		for s, row := range trans {
			for b, p := range row {
				next[(s*4+b)%n] += pi[s] * p
			}
		}
		pi, next = next, pi
	}
	var f [4]float64
	for s, row := range trans {
		for b, p := range row {
			f[b] += pi[s] * p
		}
	}
	return f
}

// HMMKmerModel generation is an unbiased sample of the smoothed chain: base
// frequencies from a 20-base training input agree with the chain's
// stationary expectation within 4 standard errors (batch means over 50
// batches of 4 kb). The A-rich k=2 case has a strongly non-uniform
// stationary distribution.
func TestHMMKmerModel_TinyTrainingMatchesStationary(t *testing.T) {
	for _, tc := range []struct {
		k   int
		seq string
	}{{4, "ACGTAGCTAGGATCCATGCA"}, {2, "AAAAAAAAAAAACCGTAAAA"}} {
		t.Run(fmt.Sprintf("k=%d", tc.k), func(t *testing.T) {
			genome := &MockReferenceGenome{chromosomes: map[string]*MockChromosome{
				"chr1": {Name: "chr1", Sequence: []byte(tc.seq), Length: len(tc.seq)},
			}}
			m := NewHMMKmerModel(tc.k)
			require.NoError(t, m.TrainFromReferenceGenome(genome))
			m.SetSeed(2)
			const batches, batchLen = 50, 4000
			res, err := m.GenerateSequence(batches * batchLen)
			require.NoError(t, err)
			want := stationaryBaseFreqs(m.transitionMatrix)
			for b := 0; b < 4; b++ {
				means := make([]float64, batches)
				for j := range means {
					means[j] = float64(baseCounts(res.Sequence[j*batchLen : (j+1)*batchLen])[b]) / batchLen
				}
				mean, ss := 0.0, 0.0
				for _, x := range means {
					mean += x
				}
				mean /= batches
				for _, x := range means {
					ss += (x - mean) * (x - mean)
				}
				se := math.Sqrt(ss / (batches - 1) / batches)
				t.Logf("%c: %.4f, stationary %.4f, se %.4f", basesATGC[b], mean, want[b], se)
				assert.LessOrEqual(t, math.Abs(mean-want[b]), 4*se, "base %c", basesATGC[b])
			}
			assert.False(t, isPeriodic(res.Sequence[:10_000], 200))
		})
	}
}

// Same seed, same output; different seeds differ (both models).
func TestHMM_SeedDeterminism(t *testing.T) {
	seq := randomDNA(8, 2000, 0)
	genome := &MockReferenceGenome{chromosomes: map[string]*MockChromosome{
		"chr1": {Name: "chr1", Sequence: seq, Length: len(seq)},
	}}
	type model interface {
		SetSeed(int64)
		GenerateSequence(int) (*GenerationResult, error)
	}
	newModels := func() map[string]model {
		a := NewAdaptiveHMMModel(4)
		require.NoError(t, a.TrainFromSequence(seq))
		k := NewHMMKmerModel(4)
		require.NoError(t, k.TrainFromReferenceGenome(genome))
		return map[string]model{"adaptive": a, "kmer": k}
	}
	gen := func(m model, seed int64) *GenerationResult {
		m.SetSeed(seed)
		r, err := m.GenerateSequence(2000)
		require.NoError(t, err)
		return r
	}
	ms1, ms2 := newModels(), newModels()
	for name := range ms1 {
		t.Run(name, func(t *testing.T) {
			a, b := gen(ms1[name], 17), gen(ms2[name], 17)
			assert.Equal(t, a.Sequence, b.Sequence)
			assert.Equal(t, a.Qualities, b.Qualities)
			c := gen(ms1[name], 18)
			assert.NotEqual(t, a.Sequence, c.Sequence)
			// Reseeding restarts the stream.
			assert.Equal(t, a.Sequence, gen(ms1[name], 17).Sequence)
		})
	}
}

// Load → SetSeed → Generate is byte-identical across two loads of an
// HMMKmerModel file.
func TestHMMKmerModel_LoadDeterministic(t *testing.T) {
	m := NewHMMKmerModel(5)
	require.NoError(t, m.TrainFromReferenceGenome(createSyntheticReferenceGenome(t)))
	path := t.TempDir() + "/k.hmm.json"
	require.NoError(t, m.Save(path))
	gen := func() *GenerationResult {
		l, err := LoadHMMKmerModel(path)
		require.NoError(t, err)
		l.SetSeed(3)
		r, err := l.GenerateSequence(3000)
		require.NoError(t, err)
		return r
	}
	a, b := gen(), gen()
	assert.Equal(t, a.Sequence, b.Sequence)
	assert.Equal(t, a.Qualities, b.Qualities)
}

// sampleIndex never returns a zero-weight index, and is uniform on an
// all-zero vector instead of returning the last index.
func TestSampleIndex(t *testing.T) {
	r := rand.New(rand.NewPCG(1, 1))
	var seen [4]int
	for i := 0; i < 4000; i++ {
		j := sampleIndex(r, []float64{0, 0.5, 0, 0.5})
		require.Contains(t, []int{1, 3}, j)
		seen[sampleIndex(r, []float64{0, 0, 0, 0})]++
	}
	for _, n := range seen {
		assert.Greater(t, n, 800)
	}
}
