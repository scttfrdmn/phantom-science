// Package signal synthesises nanopore raw current traces from k-mer pore
// models, following squigulator's algorithm (Gamaarachchi et al., "Simulation
// of nanopore sequencing signal data with tunable parameters", Genome
// Research 2024; https://github.com/hasindu2008/squigulator), reimplemented
// in Go.
//
// A read is turned into signal by walking its k-mers in translocation order
// (5'→3' for DNA, 3'→5' for direct RNA). Each k-mer dwells for a
// gamma-distributed number of samples (mean SamplingRate/SpeedBPS) and each
// sample is drawn from N(level, noise) plus a slow baseline drift; the pA
// trace is then digitised to int16 as MinKNOW does:
//
//	pA  = (raw + offset) × range / digitisation
//	raw = round(pA × digitisation / range − offset), clamped to int16
//
// Stalls (dwell × 5–20) and skips (zero dwell) are inserted at low rates, an
// adapter/leader segment is prepended, and an optional poly-A (or poly-T)
// plateau can be added for RNA / cDNA.
//
// Synthesis is deterministic: the same generator state gives the same
// signal. SynthesizeRead seeds the generator from (run seed, read ID), so any
// read's signal is reproducible from its UUID alone, independent of the
// order reads are generated in.
//
// The package uses only the standard library. Pore models are not shipped:
// scripts/fetch-pore-model.sh downloads ONT's published tables.
package signal

import (
	"fmt"
	"math"
	"sort"
)

// Preset names.
const (
	PresetDNAR10 = "dna_r10.4.1_e8.2_400bps_5khz"
	PresetRNA004 = "rna004_130bps_4khz"
)

// Profile holds the chemistry and device parameters of a synthetic run.
// Lengths in samples are derived from SamplingRate; levels are in pA.
type Profile struct {
	Name string

	// Device.
	SamplingRate float64 // Hz
	Digitisation float64 // ADC levels (2048 on PromethION, 8192 on MinION)
	Range        float64 // pA spanned by Digitisation levels
	OffsetMean   float64 // per-read ADC offset (raw units), drawn N(OffsetMean, OffsetStd)
	OffsetStd    float64
	// MedianBeforeMean/Std give the per-read open-pore current (pA) before
	// capture, reported as Read.MedianBefore (slow5 median_before).
	MedianBeforeMean float64
	MedianBeforeStd  float64

	// Translocation.
	SpeedBPS float64 // bases per second; mean dwell = SamplingRate/SpeedBPS samples
	DwellCV  float64 // coefficient of variation of the gamma dwell (0 = fixed dwell)
	// Reverse walks the sequence 3'→5' (direct RNA). Model k-mers are always
	// looked up in 5'→3' orientation, as ONT's RNA tables are stored.
	Reverse bool
	// LevelTau is the time constant, in samples, of the transition between
	// consecutive segment levels. 0 gives squigulator's step function: every
	// sample of a k-mer is centred on its model level. With τ > 0 the mean
	// level relaxes exponentially from the level reached at the end of the
	// previous segment, prev, towards the k-mer's level, target: the i-th
	// sample (0-based) of a segment is centred on
	//
	//	target + (prev − target)·exp(−(i+1)/τ)
	//
	// a first-order low-pass of the step function. A segment shorter than a
	// few τ hands an unfinished level to the next one, and a skipped k-mer
	// (no samples) leaves the level unchanged. The first segment of a read
	// starts at its own level. Event.LevelPA still records the k-mer's
	// target level.
	LevelTau float64

	// Noise.
	LevelNoiseScale float64 // multiplies the model's per-k-mer stdv (squigulator --amp-noise)
	WhiteNoisePA    float64 // extra Gaussian per-sample noise, pA
	DriftPA         float64 // stationary SD of the baseline drift (Ornstein-Uhlenbeck), pA
	DriftTau        float64 // drift correlation time, seconds
	ReadShiftSD     float64 // per-read additive level shift SD, pA
	ReadScaleSD     float64 // per-read multiplicative level scale SD (relative)

	// Events.
	StallProb float64 // per-k-mer probability of a stall
	StallMin  float64 // stall dwell multiplier range (uniform)
	StallMax  float64
	SkipProb  float64 // per-k-mer probability of a skip (zero samples)

	// Prefix and poly-A.
	AdapterSeconds float64 // total adapter segment (leader plateau + adapter k-mers); 0 = none
	AdapterLevelPA float64 // leader plateau level
	AdapterNoisePA float64 // leader plateau per-sample noise
	AdapterSeq     string  // adapter sequence translocated after the leader (may be empty)
	PolyANoisePA   float64 // per-sample noise of a poly-A/T plateau (low variance)

	// Pore-model scaling recommended for this chemistry's ONT
	// kmer_models table (standardised units → pA), used by ModelOptions.
	// ModelStdv is the per-sample level stdv, in table units, given to
	// k-mers of a table without a level_stdv column (ONT's tables have
	// none); a table's own stdv column takes precedence.
	ModelK     int
	ModelScale float64
	ModelShift float64
	ModelStdv  float64
}

// Kit 14 ligation adapter top strand as it appears at the start of R10.4.1
// DNA reads (the "LA"/"NA" adapter, as listed by Porechop/Dorado's adapter
// trimmer): poly-T leader, then the Y-adapter.
const kit14Adapter = "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"

// Presets. Device values (digitisation, sampling rate, range, offset and
// median_before distributions) are squigulator v0.4's dna-r10-prom and
// rna004-prom profiles (src/sim.c). The model scaling for ONT's
// standardised kmer_models tables is the fixed pA standardisation in
// Dorado's config.toml: mean 93.6924, stdev 23.5067 for
// dna_r10.4.1_e8.2_400bps_sup@v5.2.0, mean 80, stdev 17 for
// rna004_130bps_sup@v5.3.0.
//
// The DNA dwell, stall, skip and noise values were chosen against real
// R10.4.1 data (FLO-PRO114M, 5 kHz; ont-open-data pgx_as_2025.07) with
// scripts/acceptance/dorado_identity.sh and Dorado 2.1.2 SUP v5.2.0, measuring
// dwell from Dorado's move tables. Real reads: mean dwell 13.1 samples/base,
// CV 1.12 (a long right tail with few short dwells), 98.5% median identity.
// With these settings synthetic reads have mean dwell 13.0, CV 0.88 and
// 90.8% identity, against 91.5% for squigulator 0.5.0 on the same truth.
// Short dwells cost accuracy (gamma CV 0.5: 88.9%, CV 1.0: 67.8%) and skips
// are deletions by construction, so the dwell CV is squigulator's 4/13 and
// the tail comes from stalls; extra white noise and drift lowered identity
// without a measured basis and are off.
//
// The remaining DNA gap was ablated in scripts/acceptance/calibrate
// (results/dna.tsv; 2000 reads per condition, bootstrap 95% CIs; Dorado
// 2.1.2 SUP v5.2.0 on Metal). Real reads: 98.54% [98.52, 98.58]. Synthetic,
// this profile: 90.80% on uniform random reads, 91.75% [91.67, 91.82] on
// GRCh38 pseudo-reads (composition, +0.95). Swapping ONT's 9-mer table for
// squigulator's built-in R10 levels gives 93.81% (+2.06; squigulator with
// ONT's table loses the same 2.1 points), and dropping stalls a further
// +0.15 at the cost of a dwell CV further from real (0.61 vs 1.08). Level
// smoothing (LevelTau 0.5-8 samples: -0.3 to -19 points), noise matched to
// real reads (x1.27: -0.17), a +5 pA level shift (-0.04) and a per-k-mer
// stdv (-0.93) do not help, so LevelTau is 0 and the rest is unchanged.
// The ~4.6 points left are the k-mer level model: real per-base levels
// fit the tables at r = 0.70-0.72 against 0.85 for synthetic signal, and
// CpG k-mers sit 3.8 pA lower in native (methylated) DNA.
//
// The RNA004 values were measured on real RNA004 reads (FLO-PRO004RA,
// 4 kHz; ont-open-data rna-modbase-validation_2025.03 control_rep1, 3000
// reads, Dorado 2.1.2 rna004_130bps_sup@v5.3.0; results/rna.tsv): the
// pre-sequence segment Dorado trims is ~0.7 s at ~62 pA with ~1.5 pA
// noise, per-sample noise in the read is ~3.1 pA, and Dorado's moves give
// 37.5 samples/base (CV 1.15). On 1200 synthetic reads of the same 100-nt
// constructs the old values (1.5 s adapter at 100 pA, 1.3 pA model noise
// plus white noise and drift, skips) left 38% of reads unaligned (63.4%
// median identity of the rest, 0.12 of the target matched per read); these
// values align 72% at 72.9% (0.28 matched), against squigulator 0.5.0's
// 98% at 82.9% (0.39) and real reads' 99.9% at 97.7% (0.84). A per-k-mer
// stdv (squigulator's RNA004 table: 2.3-6.4 pA, 5-95%, through a pore model with a
// level_stdv column) lifts that to 97% at 77.8% (0.62); ONT's table has
// none, so ModelStdv here is one constant (0.176 x 17 = 3.0 pA). Level
// smoothing did not help.
var presets = map[string]Profile{
	PresetDNAR10: {
		Name:         PresetDNAR10,
		SamplingRate: 5000, Digitisation: 2048, Range: 281.345551,
		OffsetMean: -127.5655735, OffsetStd: 19.377283387665,
		MedianBeforeMean: 189.87607393756, MedianBeforeStd: 15.788097978713,
		SpeedBPS: 400, DwellCV: 0.31,
		LevelNoiseScale: 1, WhiteNoisePA: 0,
		DriftPA: 0, DriftTau: 2.0,
		ReadShiftSD: 3.0, ReadScaleSD: 0.03,
		StallProb: 0.004, StallMin: 5, StallMax: 20, SkipProb: 0,
		// 0.2 s matches the 0.1–0.3 s adapter time of nanopore.ChannelModel.
		AdapterSeconds: 0.2, AdapterLevelPA: 110, AdapterNoisePA: 8, AdapterSeq: kit14Adapter,
		PolyANoisePA: 1.0,
		ModelK:       9, ModelScale: 23.506, ModelShift: 93.692, ModelStdv: DefaultModelStdv,
	},
	PresetRNA004: {
		Name:         PresetRNA004,
		SamplingRate: 4000, Digitisation: 2048, Range: 299.432068,
		OffsetMean: -259.421128, OffsetStd: 16.010841823643,
		MedianBeforeMean: 205.63935594369, MedianBeforeStd: 8.3994882799157,
		SpeedBPS: 130, DwellCV: 0.31, Reverse: true,
		LevelNoiseScale: 1, WhiteNoisePA: 0,
		DriftPA: 0, DriftTau: 2.0,
		ReadShiftSD: 3.0, ReadScaleSD: 0.03,
		StallProb: 0.004, StallMin: 5, StallMax: 20, SkipProb: 0,
		// The segment Dorado trims before the read: median 2800 samples
		// (0.7 s) at ~62 pA, ~1.5 pA noise, in real RNA004 reads.
		AdapterSeconds: 0.7, AdapterLevelPA: 62, AdapterNoisePA: 1.5,
		PolyANoisePA: 1.0,
		ModelK:       9, ModelScale: 17, ModelShift: 80, ModelStdv: 0.176,
	},
}

// Preset returns the named profile.
func Preset(name string) (Profile, error) {
	p, ok := presets[name]
	if !ok {
		return Profile{}, fmt.Errorf("unknown signal profile %q (have %v)", name, PresetNames())
	}
	return p, nil
}

// PresetNames lists the preset names in sorted order.
func PresetNames() []string {
	names := make([]string, 0, len(presets))
	for n := range presets {
		names = append(names, n)
	}
	sort.Strings(names)
	return names
}

// tunable maps the parameter names accepted by Set to their fields.
func (p *Profile) tunable() map[string]*float64 {
	return map[string]*float64{
		"speed_bps": &p.SpeedBPS, "dwell_cv": &p.DwellCV,
		"level_noise_scale": &p.LevelNoiseScale, "white_noise_pa": &p.WhiteNoisePA,
		"drift_pa": &p.DriftPA, "drift_tau": &p.DriftTau,
		"read_shift_sd": &p.ReadShiftSD, "read_scale_sd": &p.ReadScaleSD,
		"stall_prob": &p.StallProb, "stall_min": &p.StallMin, "stall_max": &p.StallMax,
		"skip_prob": &p.SkipProb, "level_tau": &p.LevelTau,
		"adapter_seconds": &p.AdapterSeconds, "adapter_level_pa": &p.AdapterLevelPA,
		"adapter_noise_pa": &p.AdapterNoisePA, "polya_noise_pa": &p.PolyANoisePA,
		"offset_mean": &p.OffsetMean, "offset_std": &p.OffsetStd,
		"median_before_mean": &p.MedianBeforeMean, "median_before_std": &p.MedianBeforeStd,
		"model_scale": &p.ModelScale, "model_shift": &p.ModelShift, "model_stdv": &p.ModelStdv,
	}
}

// ParameterNames lists the parameter names accepted by Set, sorted.
func ParameterNames() []string {
	var p Profile
	names := make([]string, 0, len(p.tunable()))
	for n := range p.tunable() {
		names = append(names, n)
	}
	sort.Strings(names)
	return names
}

// Set overrides one modelling parameter by name (see ParameterNames), for
// calibrating a profile against real data. Device geometry (sampling rate,
// digitisation, range) and k cannot be changed. Call Validate afterwards.
func (p *Profile) Set(name string, v float64) error {
	f, ok := p.tunable()[name]
	if !ok {
		return fmt.Errorf("unknown signal profile parameter %q (have %v)", name, ParameterNames())
	}
	*f = v
	return nil
}

// MeanDwell is the mean number of samples per base, SamplingRate/SpeedBPS.
func (p Profile) MeanDwell() float64 { return p.SamplingRate / p.SpeedBPS }

// ModelOptions returns load options that convert this chemistry's ONT
// kmer_models table (standardised units) to pA, give its k-mers the
// profile's ModelStdv, and check its k.
func (p Profile) ModelOptions() []LoadOption {
	opts := []LoadOption{WithScaling(p.ModelScale, p.ModelShift), WithDefaultStdv(p.ModelStdv)}
	if p.ModelK > 0 {
		opts = append(opts, WithK(p.ModelK))
	}
	return opts
}

// Digitiser returns the profile's ADC conversion for a given per-read offset.
func (p Profile) Digitiser(offset float64) Digitiser {
	return Digitiser{Digitisation: p.Digitisation, Offset: offset, Range: p.Range}
}

// Validate reports parameters that would make synthesis meaningless.
func (p Profile) Validate() error {
	pos := []struct {
		name string
		v    float64
	}{{"sampling_rate", p.SamplingRate}, {"digitisation", p.Digitisation}, {"range", p.Range}, {"speed_bps", p.SpeedBPS}}
	for _, f := range pos {
		if !(f.v > 0) || math.IsInf(f.v, 0) {
			return fmt.Errorf("signal profile %s must be positive and finite, got %g", f.name, f.v)
		}
	}
	nonneg := []struct {
		name string
		v    float64
	}{
		{"offset_std", p.OffsetStd}, {"median_before_std", p.MedianBeforeStd}, {"dwell_cv", p.DwellCV},
		{"level_noise_scale", p.LevelNoiseScale}, {"white_noise_pa", p.WhiteNoisePA}, {"level_tau", p.LevelTau},
		{"drift_pa", p.DriftPA}, {"read_shift_sd", p.ReadShiftSD}, {"read_scale_sd", p.ReadScaleSD},
		{"adapter_seconds", p.AdapterSeconds}, {"adapter_noise_pa", p.AdapterNoisePA}, {"polya_noise_pa", p.PolyANoisePA},
		{"model_stdv", p.ModelStdv},
	}
	for _, f := range nonneg {
		if !(f.v >= 0) || math.IsInf(f.v, 0) {
			return fmt.Errorf("signal profile %s must be non-negative and finite, got %g", f.name, f.v)
		}
	}
	if p.DriftPA > 0 && !(p.DriftTau > 0) {
		return fmt.Errorf("signal profile drift_tau must be positive when drift_pa > 0, got %g", p.DriftTau)
	}
	if p.MeanDwell() < 1 {
		return fmt.Errorf("signal profile mean dwell %.3g samples/base < 1 (sampling_rate/speed_bps)", p.MeanDwell())
	}
	for _, f := range []struct {
		name string
		v    float64
	}{{"stall_prob", p.StallProb}, {"skip_prob", p.SkipProb}} {
		if !(f.v >= 0 && f.v <= 1) {
			return fmt.Errorf("signal profile %s must be in [0,1], got %g", f.name, f.v)
		}
	}
	if p.StallProb+p.SkipProb > 1 {
		return fmt.Errorf("signal profile stall_prob + skip_prob = %g > 1", p.StallProb+p.SkipProb)
	}
	if p.StallProb > 0 && !(p.StallMin >= 1 && p.StallMax >= p.StallMin) {
		return fmt.Errorf("signal profile stall multiplier range [%g, %g] must satisfy 1 <= min <= max", p.StallMin, p.StallMax)
	}
	return nil
}
