package genomics

// Quality and sequencing-error extension points (#11, #27).
//
// Contract:
//
//   - qualityModel is called once per output stream and returns a function
//     that produces n Phred scores (0..93, not ASCII) for one read of mate
//     (1, or 2 for the second read of a pair). It must draw only from r.
//   - errorModel is called once per output stream and returns a function
//     applied to every read after variant injection and trimming, before the
//     header is formatted. It receives the read's bases and Phred qualities
//     (same length), may modify them in place or return new slices of any
//     (equal) length, and must draw only from r. A nil function means no
//     sequencing errors.
//
// Platform models:
//
//   - ONT with a basecaller tier (ErrorModel.Tier): the tier's quality
//     distribution and its substitution, insertion, deletion and
//     homopolymer-deletion rates (package nanopore/errmodel). Without a
//     tier, ONT reads get the linear decay and no errors.
//   - Illumina, quality_model "binned" (the default): RTA3-binned per-cycle
//     qualities (Q2/Q12/Q23/Q37) from package illumina, read 2 lower than
//     read 1, and with errors "quality" (the default for binned) a
//     substitution with probability 10^(−Q/10) per base and Q2 bases
//     written as N. The binned model replaces HMM qualities.
//   - Illumina, quality_model "linear": the earlier model,
//     BaseQuality − i·QualityDecay clamped to [2, 41], which draws nothing;
//     errors default to "none". Single-end linear output is byte-identical
//     to the output before #27 (HMM qualities are kept when the model has
//     them).
//
// A basecaller tier on a non-ONT platform is a configuration error (it used
// to apply the ONT model to Illumina reads).
//
// Either factory returns an error for an unsupported configuration (wrap
// ErrNotImplemented for planned features).

import (
	"fmt"
	"math/rand/v2"

	"github.com/scttfrdmn/phantom-science/internal/genomics/illumina"
	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore/errmodel"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// tierModel returns the ONT error model of the configured basecaller tier,
// or nil when no tier is set.
func (e *Engine) tierModel() (*errmodel.Model, error) {
	tier := e.config.ErrorModel.Tier
	if tier == "" {
		return nil, nil
	}
	if !e.isONT() {
		return nil, fmt.Errorf("error_model tier %q is an Oxford Nanopore basecaller tier: it needs platform ont, not %s (Illumina errors are set with illumina.errors)", tier, e.platform())
	}
	t, err := errmodel.ParseTier(tier)
	if err != nil {
		return nil, err
	}
	return errmodel.New(t)
}

// illuminaSpec returns the Illumina read model with defaults applied.
func (e *Engine) illuminaSpec() types.IlluminaSpec {
	return e.config.Illumina.WithDefaults()
}

// checkIllumina validates the illumina section: it must be valid, and is
// only accepted on the Illumina platform.
func (e *Engine) checkIllumina() error {
	if e.isONT() {
		if e.config.Illumina != (types.IlluminaSpec{}) {
			return fmt.Errorf("the illumina section applies to platform illumina, not ont")
		}
		return nil
	}
	return e.config.Illumina.Validate()
}

// binnedQualities reports whether reads get Illumina binned qualities.
func (e *Engine) binnedQualities() bool {
	return !e.isONT() && e.illuminaSpec().QualityModel == types.IlluminaQualityBinned
}

// qualityModel returns the quality generator for reads without HMM
// qualities (and, for binned Illumina qualities, for every read).
func (e *Engine) qualityModel() (func(r *rand.Rand, n, mate int) []byte, error) {
	m, err := e.tierModel()
	if err != nil {
		return nil, err
	}
	if m != nil {
		return func(r *rand.Rand, n, _ int) []byte { return m.Qualities(r, n) }, nil
	}
	if e.binnedQualities() {
		r1, r2 := illumina.ParamsFor(1), illumina.ParamsFor(2)
		return func(r *rand.Rand, n, mate int) []byte {
			if mate == 2 {
				return r2.Qualities(r, n)
			}
			return r1.Qualities(r, n)
		}, nil
	}
	base, decay := e.config.BaseQuality, e.config.QualityDecay
	return func(_ *rand.Rand, n, _ int) []byte {
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

// errorModel returns the sequencing-error model, or nil for none: the ONT
// tier model when a tier is set (qualities are kept as baselines and
// penalised at and next to errors, see package errmodel), else, on
// Illumina with errors "quality", quality-matched substitutions.
func (e *Engine) errorModel() (func(r *rand.Rand, seq, qual []byte) ([]byte, []byte), error) {
	m, err := e.tierModel()
	if err != nil {
		return nil, err
	}
	if m != nil {
		return m.ApplyQualities, nil
	}
	if !e.isONT() && e.illuminaSpec().Errors == types.IlluminaErrorsQuality {
		return illumina.ApplyErrors, nil
	}
	return nil, nil
}
