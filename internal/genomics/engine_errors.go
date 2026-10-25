package genomics

// Quality and sequencing-error extension points (#11, #27).
//
// Contract:
//
//   - qualityModel is called once per output stream and returns a function
//     that produces n Phred scores (0..93, not ASCII) for one read of mate
//     (1, or 2 for the second read of a pair). seq is the read's bases in
//     sequencing orientation (len(seq) == n) before sequencing errors, or
//     nil where qualities are drawn before the bases (the legacy flows of
//     ONT and linear Illumina qualities). class is the read's quality class
//     drawn for the pair (pairClasses), or -1 for the model to draw it. It
//     must draw only from r.
//   - pairClasses is called once per output stream and returns nil, or a
//     function that draws both mates' quality classes of one pair from r
//     (the pair's shared stream), so that the R1 and R2 files regenerate
//     the same classes.
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
//   - Illumina, quality_model "binned" (the default): qualities in the
//     NovaSeq 6000 bins (Q2/Q11/Q25/Q37) from package illumina, fitted to
//     a real run (illumina.profile, #35; by default
//     novaseq6000-pooled, three runs pooled, #36) and drawn
//     given the read's bases (each cycle's Q11 and
//     Q25 probabilities depend on its 4-base context, #35; the two reads of
//     a pair draw their quality classes together, #35), and with errors
//     "quality" (the default for binned)
//     substitutions at the measured rate of each bin and cycle, with the
//     measured spectrum, and Q2 bases written as N. The binned model
//     replaces HMM qualities.
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
// only accepted on the Illumina platform; a profile must exist and needs
// the binned quality model.
func (e *Engine) checkIllumina() error {
	if e.isONT() {
		if e.config.Illumina != (types.IlluminaSpec{}) {
			return fmt.Errorf("the illumina section applies to platform illumina, not ont")
		}
		return nil
	}
	if err := e.config.Illumina.Validate(); err != nil {
		return err
	}
	if e.config.Illumina.Profile != "" && e.illuminaSpec().QualityModel != types.IlluminaQualityBinned {
		return fmt.Errorf("illumina profile %q needs quality_model binned, not %s", e.config.Illumina.Profile, e.illuminaSpec().QualityModel)
	}
	_, err := e.illuminaProfile()
	return err
}

// illuminaProfile returns the configured Illumina run profile (the default
// when none is set).
func (e *Engine) illuminaProfile() (*illumina.Profile, error) {
	return illumina.LookupProfile(e.config.Illumina.Profile)
}

// insertSizes returns the insert-size distribution of paired-end reads, with
// defaults applied: illumina.insert_size when set (any field); else, with
// illumina.profile set and binned qualities, that profile's measured one;
// else the types default, lognormal 455 ± 100. An unset profile thus gives
// the default profile's qualities (illumina.DefaultProfile) but not its
// run's inserts; see the illumina package's profile documentation.
func (e *Engine) insertSizes() (types.InsertSizeSpec, error) {
	il := e.config.Illumina
	if il.InsertSize != (types.InsertSizeSpec{}) || il.Profile == "" || !e.binnedQualities() {
		return e.illuminaSpec().InsertSize, nil
	}
	prof, err := e.illuminaProfile()
	if err != nil {
		return types.InsertSizeSpec{}, err
	}
	return prof.InsertSize, nil
}

// binnedQualities reports whether reads get Illumina binned qualities.
func (e *Engine) binnedQualities() bool {
	return !e.isONT() && e.illuminaSpec().QualityModel == types.IlluminaQualityBinned
}

// qualityModel returns the quality generator for reads without HMM
// qualities (and, for binned Illumina qualities, for every read).
func (e *Engine) qualityModel() (func(r *rand.Rand, seq []byte, n, mate, class int) []byte, error) {
	m, err := e.tierModel()
	if err != nil {
		return nil, err
	}
	if m != nil {
		return func(r *rand.Rand, _ []byte, n, _, _ int) []byte { return m.Qualities(r, n) }, nil
	}
	if e.binnedQualities() {
		prof, err := e.illuminaProfile()
		if err != nil {
			return nil, err
		}
		r1, r2 := prof.Read1, prof.Read2
		return func(r *rand.Rand, seq []byte, n, mate, class int) []byte {
			q := r1
			if mate == 2 {
				q = r2
			}
			if len(seq) != n {
				return q.Qualities(r, n)
			}
			return q.QualitiesForClass(r, seq, class)
		}, nil
	}
	base, decay := e.config.BaseQuality, e.config.QualityDecay
	return func(_ *rand.Rand, _ []byte, n, _, _ int) []byte {
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

// pairClasses returns the draw of a pair's two quality classes, or nil when
// the mates' qualities are drawn independently (single-end reads, models
// other than binned Illumina qualities). It draws one value from r.
func (e *Engine) pairClasses() func(r *rand.Rand) [2]int {
	if !e.config.PairedEnd || !e.binnedQualities() {
		return nil
	}
	prof, err := e.illuminaProfile()
	if err != nil || !prof.Pairs.Coupled() {
		return nil // an unknown profile fails in qualityModel
	}
	pp := prof.Pairs
	return func(r *rand.Rand) [2]int {
		k1, k2 := pp.Classes(r.Float64())
		return [2]int{k1, k2}
	}
}

// errorModel returns the sequencing-error model, or nil for none: the ONT
// tier model when a tier is set (qualities are kept as baselines and
// penalised at and next to errors, see package errmodel), else, on
// Illumina with errors "quality", quality-matched substitutions (the
// profile's, for binned qualities).
func (e *Engine) errorModel() (func(r *rand.Rand, seq, qual []byte) ([]byte, []byte), error) {
	m, err := e.tierModel()
	if err != nil {
		return nil, err
	}
	if m != nil {
		return m.ApplyQualities, nil
	}
	if !e.isONT() && e.illuminaSpec().Errors == types.IlluminaErrorsQuality {
		if !e.binnedQualities() {
			return illumina.ApplyErrors, nil
		}
		prof, err := e.illuminaProfile()
		if err != nil {
			return nil, err
		}
		return prof.Errors.Apply, nil
	}
	return nil, nil
}
