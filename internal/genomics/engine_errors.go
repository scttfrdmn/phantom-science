package genomics

// Quality and sequencing-error extension points (#11).
//
// Contract:
//
//   - qualityModel is called once per output stream and returns a function
//     that produces n Phred scores (0..93, not ASCII) for a read of n bases
//     whose sequence did not come with HMM qualities. It must draw only from
//     r. The default is the linear decay BaseQuality - i*QualityDecay,
//     clamped to [2, 41], and draws nothing.
//   - errorModel is called once per output stream and returns a function
//     applied to every read after variant injection and trimming, before the
//     header is formatted. It receives the read's bases and Phred qualities
//     (same length), may modify them in place or return new slices of any
//     (equal) length, and must draw only from r. A nil function means no
//     sequencing errors, which is the default.
//
// Either factory returns an error for an unsupported configuration (wrap
// ErrNotImplemented for planned features).

import (
	"fmt"
	"math/rand/v2"
)

// Basecaller tiers for ErrorModelSpec.Tier.
const (
	tierSup  = "sup"
	tierHac  = "hac"
	tierFast = "fast"
)

// qualityModel returns the quality generator for reads without HMM
// qualities.
func (e *Engine) qualityModel() (func(r *rand.Rand, n int) []byte, error) {
	base, decay := e.config.BaseQuality, e.config.QualityDecay
	return func(_ *rand.Rand, n int) []byte {
		q := make([]byte, n)
		for i := range q {
			v := base - int(float64(i)*decay)
			if v < 2 {
				v = 2
			}
			if v > maxIlluminaPhred {
				v = maxIlluminaPhred
			}
			q[i] = byte(v) // #nosec G115 -- clamped to [2, 41]
		}
		return q
	}, nil
}

// errorModel returns the sequencing-error model, or nil for none.
func (e *Engine) errorModel() (func(r *rand.Rand, seq, qual []byte) ([]byte, []byte), error) {
	switch tier := e.config.ErrorModel.Tier; tier {
	case "":
		return nil, nil
	case tierSup, tierHac, tierFast:
		return nil, fmt.Errorf("basecaller tier %q error model is %w (#11)", tier, ErrNotImplemented)
	default:
		return nil, fmt.Errorf("unknown basecaller tier %q: want sup, hac or fast", tier)
	}
}
