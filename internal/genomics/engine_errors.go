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
	"math/rand/v2"

	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore/errmodel"
)

// tierModel returns the ONT error model of the configured basecaller tier,
// or nil when no tier is set.
func (e *Engine) tierModel() (*errmodel.Model, error) {
	tier := e.config.ErrorModel.Tier
	if tier == "" {
		return nil, nil
	}
	t, err := errmodel.ParseTier(tier)
	if err != nil {
		return nil, err
	}
	return errmodel.New(t)
}

// qualityModel returns the quality generator for reads without HMM
// qualities: the basecaller tier's distribution when a tier is set, else the
// linear decay.
func (e *Engine) qualityModel() (func(r *rand.Rand, n int) []byte, error) {
	m, err := e.tierModel()
	if err != nil {
		return nil, err
	}
	if m != nil {
		return m.Qualities, nil
	}
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

// errorModel returns the sequencing-error model, or nil for none. With a
// basecaller tier set it applies that tier's substitution, insertion,
// deletion and homopolymer-deletion rates; the read's qualities are kept as
// baselines and penalised at and next to errors (see package errmodel).
func (e *Engine) errorModel() (func(r *rand.Rand, seq, qual []byte) ([]byte, []byte), error) {
	m, err := e.tierModel()
	if err != nil || m == nil {
		return nil, err
	}
	return m.ApplyQualities, nil
}
