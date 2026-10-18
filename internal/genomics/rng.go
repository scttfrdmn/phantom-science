package genomics

// Seedable randomness for read generation.
//
// Every output stream draws from its own math/rand/v2 PCG generator derived
// from (run seed, stream ID), so output is reproducible from the seed and does
// not depend on the order (or, later, the concurrency) in which streams are
// generated. crypto/rand is used only to pick a seed when none is configured.

import (
	crand "crypto/rand"
	"encoding/binary"
	"fmt"
	"io"
	"math/rand/v2"
)

// Stream IDs. File output uses one stream per (batch, mate) file, plus, for
// paired-end output, the batch's pair stream (mate pairMate) for the
// per-read values both mate files share; streamed output (GenerateToWriter)
// uses one stream per batch with mate 0. Streams with the high bit set are
// reserved for engine-internal use.
const (
	rawStream uint64 = 1 << 63 // GenerateRawSequence
	pairMate         = 3       // mate slot of a batch's pair stream
)

// flowCellTag separates flow-cell seed derivation from every other use of
// the run seed.
const flowCellTag uint64 = 0x666c6f7763656c6c // "flowcell"

// flowCellSeed returns the run seed of flow cell k (1-based) of a
// multi-flow-cell run; attempt > 0 re-derives it when an earlier seed gave
// a flow-cell ID already in use. The result is never zero.
func flowCellSeed(seed uint64, k, attempt int) uint64 {
	// #nosec G115 -- k and attempt are small and positive
	s := splitmix64(seed ^ splitmix64(flowCellTag+uint64(k)<<16+uint64(attempt)))
	if s == 0 {
		s = 1
	}
	return s
}

// streamID returns the RNG stream for one output batch and mate (0 for an
// interleaved stream holding both mates).
func streamID(batch, mate int) uint64 {
	return uint64(batch)<<2 | uint64(mate&3) // #nosec G115 -- batch and mate are positive
}

// newStream returns an independent, deterministic generator for stream of
// the run seeded with seed.
func newStream(seed, stream uint64) *rand.Rand {
	// #nosec G404 -- simulation randomness, reproducible by design
	return rand.New(rand.NewPCG(seed, splitmix64(stream)))
}

// splitmix64 scrambles stream IDs so that neighbouring IDs give unrelated
// PCG increments.
func splitmix64(x uint64) uint64 {
	x += 0x9e3779b97f4a7c15
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9
	x = (x ^ (x >> 27)) * 0x94d049bb133111eb
	return x ^ (x >> 31)
}

// ResolveSeed returns seed, or a random non-zero seed drawn from crypto/rand
// when seed is zero; drawn reports whether a seed was drawn.
func ResolveSeed(seed uint64) (resolved uint64, drawn bool, err error) {
	if seed != 0 {
		return seed, false, nil
	}
	var b [8]byte
	for resolved == 0 {
		if _, err := io.ReadFull(crand.Reader, b[:]); err != nil {
			return 0, false, fmt.Errorf("drawing a random seed: %w", err)
		}
		resolved = binary.LittleEndian.Uint64(b[:])
	}
	return resolved, true, nil
}

// Seed returns the run seed. Passing it as SequencerConfig.Seed (or --seed)
// reproduces this engine's output.
func (e *Engine) Seed() uint64 { return e.seed }

// Rand returns a deterministic generator for stream of this engine's run.
// Packages that build their own read loops (clinical scenarios) use it so
// their output is reproducible from the seed too. Streams with the high bit
// set are reserved for the engine.
func (e *Engine) Rand(stream uint64) *rand.Rand {
	return newStream(e.seed, stream)
}

// logSeed reports a drawn seed once, so the run can be reproduced.
func (e *Engine) logSeed() {
	if !e.seedDrawn {
		return
	}
	e.seedLogOnce.Do(func() {
		fmt.Fprintf(e.logw, "   [info] seed %d (pass --seed %d to reproduce this run)\n", e.seed, e.seed)
	})
}

// seedable is implemented by HMM models whose generator can be reseeded.
type seedable interface {
	SetSeed(seed int64)
}

// reseed reseeds model (when it supports seeding) from r, so HMM-generated
// reads are reproducible from the run seed as well.
func reseed(model any, r *rand.Rand) {
	if s, ok := model.(seedable); ok {
		s.SetSeed(r.Int64())
	}
}

// randomBases fills a new slice of n bases drawn uniformly from ACGT.
func randomBases(r *rand.Rand, n int) []byte {
	const bases = "ATCG"
	seq := make([]byte, n)
	for i := 0; i < n; {
		v := r.Uint64()
		for k := 0; k < 32 && i < n; k++ {
			seq[i] = bases[v&3]
			v >>= 2
			i++
		}
	}
	return seq
}
