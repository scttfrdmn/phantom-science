package illumina

// Run profiles (#35). Each profile is the whole model fitted, by the same
// pipeline (scripts/acceptance/illumina_calibrate.sh), to the first 200,000
// read pairs of one cited NovaSeq 6000 run: per-cycle bin curves, read
// classes and their mate coupling, stickiness, context factors, the
// substitution model and the run's insert sizes. The runs differ: the
// quality level, the share of Q11, how often the quality changes between
// cycles and the Q37 error rate are properties of the run, and the insert
// size of its library. Profiles are named for the instrument and the run's
// accession.
//
// # The default
//
// DefaultProfile, used when illumina.profile is not set, is
// novaseq6000-drr394916 (TruSeq DNA PCR-Free; mean Q 35.29 / 34.60): of the
// fitted runs it is the most central against NovaSeq 6000 runs that were
// never fitted (#35). Releases before it defaulted to novaseq6000-drr276626,
// which is still available by name.
//
// The insert size is resolved separately, because it belongs to the
// library, not to the qualities:
//
//   - illumina.insert_size set (any field): it wins, whatever the profile,
//     with types.IlluminaSpec.WithDefaults filling the fields left empty;
//   - else illumina.profile set (any profile, novaseq6000-drr394916
//     included): that profile's measured inserts (Profile.InsertSize;
//     novaseq6000-drr394916 lognormal 560 ± 115);
//   - else (no profile, no insert size): lognormal 455 ± 100
//     (types.DefaultInsertSizeMean and DefaultInsertSizeSD), not the
//     default profile's 560 ± 115, so that default paired output keeps the
//     insert sizes of earlier releases.
//
// So an unset profile and an explicit novaseq6000-drr394916 give the same
// qualities, mate coupling, errors and context factors, and differ only in
// the inserts of paired reads. The quality_model linear does not use a
// profile (its substitution model, with illumina.errors quality, is
// ApplyErrors, unchanged).

import (
	"fmt"
	"slices"
	"strings"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// RunInfo describes the run a profile was fitted to (ENA metadata).
type RunInfo struct {
	Accession, Study string
	Instrument       string
	Center           string
	// Library is the library construction protocol as submitted.
	Library    string
	ReadLength int
	// Pairs is the number of read pairs the profile was fitted to.
	Pairs int
}

// RunStats are headline statistics of the measured reads: per mate (read
// 1, read 2) the mean quality, the fraction of bases at Q30 or more and of
// Q11 and Q25 bases, and the proper pairs' insert sizes.
type RunStats struct {
	MeanQ, FracQ30, FracQ11, FracQ25   [2]float64
	InsertMedian, InsertMean, InsertSD float64
}

// Profile is the Illumina model fitted to one run.
type Profile struct {
	// Name selects the profile (illumina.profile), e.g.
	// "novaseq6000-drr276626".
	Name  string
	Run   RunInfo
	Stats RunStats
	// Read1 and Read2 are the quality models of the two mates (Read1 also
	// of single-end reads); Pairs couples their classes.
	Read1, Read2 QualityParams
	Pairs        PairParams
	Errors       ErrorParams
	// InsertSize is the run's measured insert-size distribution, used when
	// illumina.insert_size is not set.
	InsertSize types.InsertSizeSpec

	// freq is the 4-base context frequency of the measured reads, over
	// which the context factors average to 1.
	freq *[256]float64
}

// DefaultProfile is the profile used when illumina.profile is not set.
// Its qualities, mate coupling, errors and context factors are used; its
// insert sizes are not (see the package's profile documentation above:
// unset, the insert size is lognormal 455 ± 100).
const DefaultProfile = "novaseq6000-drr394916"

// profiles lists the fitted profiles, the default first.
var profiles = []*Profile{&profileDRR394916, &profileDRR276626, &profileERR15600747}

// ProfileNames returns the names of the fitted profiles, the default first.
func ProfileNames() []string {
	out := make([]string, len(profiles))
	for i, p := range profiles {
		out[i] = p.Name
	}
	return out
}

// Profiles returns the fitted profiles, the default first.
func Profiles() []*Profile { return slices.Clone(profiles) }

// LookupProfile returns the profile called name (case-insensitive), or the
// default for "".
func LookupProfile(name string) (*Profile, error) {
	if name == "" {
		name = DefaultProfile
	}
	for _, p := range profiles {
		if strings.EqualFold(p.Name, name) {
			return p, nil
		}
	}
	return nil, fmt.Errorf("unknown illumina profile %q: want one of %s", name, strings.Join(ProfileNames(), ", "))
}

// Params returns the quality model of mate (1 or 2; 1 for single-end
// reads).
func (p *Profile) Params(mate int) QualityParams {
	if mate == 2 {
		return p.Read2
	}
	return p.Read1
}

// The default profile's (DefaultProfile's) models.
var (
	// Read1 is the default quality model of read 1 (of single-end reads
	// too).
	Read1 = profileDRR394916.Read1
	// Read2 is the default quality model of read 2.
	Read2 = profileDRR394916.Read2
	// NovaSeqPairs couples the default models' read classes.
	NovaSeqPairs = profileDRR394916.Pairs
	// NovaSeqErrors is the default profile's substitution model.
	NovaSeqErrors = profileDRR394916.Errors
)

// linearErrors is the substitution model of the linear quality model
// (ApplyErrors): novaseq6000-drr276626's, the default profile when the
// linear model gained quality-matched errors, kept so that linear output
// does not change with the default profile.
var linearErrors = profileDRR276626.Errors

// ParamsFor returns the default quality parameters of mate (1 or 2).
func ParamsFor(mate int) QualityParams {
	if mate == 2 {
		return Read2
	}
	return Read1
}
