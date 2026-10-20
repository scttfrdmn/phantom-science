//go:build !release

package authenticity

import (
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// Models never seed from the clock: two fresh (or freshly loaded) models
// that are never reseeded generate the same output, and that output equals
// SetSeed(DefaultSeed).
func TestHMM_DefaultSeedIsFixed(t *testing.T) {
	seq := randomDNA(9, 3000, 0)
	genome := &MockReferenceGenome{chromosomes: map[string]*MockChromosome{
		"chr1": {Name: "chr1", Sequence: seq, Length: len(seq)},
	}}
	newAdaptive := func() *AdaptiveHMMModel {
		m := NewAdaptiveHMMModel(4)
		require.NoError(t, m.TrainFromSequence(seq))
		return m
	}
	newKmer := func() *HMMKmerModel {
		m := NewHMMKmerModel(4)
		require.NoError(t, m.TrainFromReferenceGenome(genome))
		return m
	}
	gen := func(m SequenceGenerator) []byte {
		r, err := m.GenerateSequence(1000)
		require.NoError(t, err)
		return append(r.Sequence, r.Qualities...)
	}

	a1, a2, a3 := newAdaptive(), newAdaptive(), newAdaptive()
	a3.SetSeed(DefaultSeed)
	assert.Equal(t, gen(a1), gen(a2))
	assert.Equal(t, gen(newAdaptive()), gen(a3))

	k1, k2, k3 := newKmer(), newKmer(), newKmer()
	k3.SetSeed(DefaultSeed)
	assert.Equal(t, gen(k1), gen(k2))
	assert.Equal(t, gen(newKmer()), gen(k3))

	dir := t.TempDir()
	paths := []string{filepath.Join(dir, "a.hmm.json"), filepath.Join(dir, "k.hmm.json")}
	require.NoError(t, a1.Save(paths[0]))
	require.NoError(t, k1.Save(paths[1]))
	for _, p := range paths {
		l1, err := LoadModel(p)
		require.NoError(t, err)
		l2, err := LoadModel(p)
		require.NoError(t, err)
		assert.Equal(t, gen(l1), gen(l2), p)
	}
}
