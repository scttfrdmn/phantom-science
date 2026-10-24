package illumina

// Run profiles (#35). Each profile is the whole model fitted, by the same
// pipeline (scripts/acceptance/illumina_calibrate.sh), to the first 200,000
// read pairs of one cited NovaSeq 6000 run: per-cycle bin curves, read
// classes and their mate coupling, stickiness, context factors, the
// substitution model and the run's insert sizes. The runs differ: the
// quality level, the share of Q11, how often the quality changes between
// cycles and the Q37 error rate are properties of the run, and the insert
// size of its library. Profiles are named for the instrument and the run's
// accession; DefaultProfile is the model of earlier releases.

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
const DefaultProfile = "novaseq6000-drr276626"

// profiles lists the fitted profiles, the default first.
var profiles = []*Profile{&profileDRR276626, &profileDRR394916, &profileERR15600747}

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

// The default profile's models, as in earlier releases.
var (
	// Read1 is the default quality model of read 1 (of single-end reads
	// too).
	Read1 = profileDRR276626.Read1
	// Read2 is the default quality model of read 2.
	Read2 = profileDRR276626.Read2
	// NovaSeqPairs couples the default models' read classes.
	NovaSeqPairs = profileDRR276626.Pairs
	// NovaSeqErrors is the default substitution model.
	NovaSeqErrors = profileDRR276626.Errors
)

// ParamsFor returns the default quality parameters of mate (1 or 2).
func ParamsFor(mate int) QualityParams {
	if mate == 2 {
		return Read2
	}
	return Read1
}
