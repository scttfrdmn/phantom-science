package genomics

// Read-length extension point (#10).
//
// Contract:
//
//   - lengthSampler is called once per output stream (one file, or one batch
//     of streamed output) and returns a function the pipeline calls once per
//     read (once per pair for paired-end reads) with that stream's generator.
//     It must draw only from r, so output stays reproducible from the seed.
//     Returned lengths must be >= 1. An unsupported configuration returns an
//     error (wrap ErrNotImplemented for planned features).
//   - fixedLengthReads reports whether reads must come out at exactly the
//     sampled length. When true (Illumina), the pipeline generates slack
//     bases to absorb deletions and trims back after variant injection, so
//     indels never change the read length and no 'N' padding is needed. When
//     false (long reads), the read is generated at the sampled length and an
//     indel changes the read length by exactly the indel size.
//
// The default reproduces fixed-length output: every read is
// SequencerConfig.ReadLength bases.

import (
	"fmt"
	"math/rand/v2"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// lengthSampler returns the per-read length sampler for the configured
// read-length distribution.
func (e *Engine) lengthSampler() (func(r *rand.Rand) int, error) {
	spec := e.config.ReadLengthDistribution
	switch spec.Dist {
	case "", types.ReadLengthFixed:
		n := e.config.ReadLength
		return func(*rand.Rand) int { return n }, nil
	case types.ReadLengthLognormal, types.ReadLengthGamma, types.ReadLengthEmpirical:
		return nil, fmt.Errorf("read-length distribution %q is %w (#10)", spec.Dist, ErrNotImplemented)
	default:
		return nil, fmt.Errorf("unknown read-length distribution %q: want fixed, lognormal, gamma or empirical", spec.Dist)
	}
}

// fixedLengthReads reports whether every read is trimmed to exactly its
// sampled length after variant injection.
func (e *Engine) fixedLengthReads() bool {
	return true
}
