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
	ModelK     int
	ModelScale float64
	ModelShift float64
}

// Kit 14 ligation adapter top strand as it appears at the start of R10.4.1
// DNA reads (the "LA"/"NA" adapter, as listed by Porechop/Dorado's adapter
// trimmer): poly-T leader, then the Y-adapter.
const kit14Adapter = "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"

// Presets. Device values (digitisation, sampling rate, range, offset and
// median_before distributions) are squigulator v0.4's dna-r10-prom and
// rna004-prom profiles (src/sim.c). The model scaling for ONT's
// standardised kmer_models tables is Dorado's per-model standardisation
// (dna_r10.4.1_e8.2_400bps_*@v5: mean 93.69 pA, stdev 23.51 pA); the
// RNA004 values are approximate. Dwell CV, noise, drift, stall/skip rates
// and the adapter leader are modelling choices, not measured values: they
// are starting points to be calibrated against real pod5 data with
// scripts/acceptance/dorado_identity.sh.
var presets = map[string]Profile{
	PresetDNAR10: {
		Name:         PresetDNAR10,
		SamplingRate: 5000, Digitisation: 2048, Range: 281.345551,
		OffsetMean: -127.5655735, OffsetStd: 19.377283387665,
		MedianBeforeMean: 189.87607393756, MedianBeforeStd: 15.788097978713,
		SpeedBPS: 400, DwellCV: 0.5,
		LevelNoiseScale: 1, WhiteNoisePA: 1.0,
		DriftPA: 1.0, DriftTau: 2.0,
		ReadShiftSD: 3.0, ReadScaleSD: 0.03,
		StallProb: 0.004, StallMin: 5, StallMax: 20, SkipProb: 0.004,
		// 0.2 s matches the 0.1–0.3 s adapter time of nanopore.ChannelModel.
		AdapterSeconds: 0.2, AdapterLevelPA: 110, AdapterNoisePA: 8, AdapterSeq: kit14Adapter,
		PolyANoisePA: 1.0,
		ModelK:       9, ModelScale: 23.506, ModelShift: 93.692,
	},
	PresetRNA004: {
		Name:         PresetRNA004,
		SamplingRate: 4000, Digitisation: 2048, Range: 299.432068,
		OffsetMean: -259.421128, OffsetStd: 16.010841823643,
		MedianBeforeMean: 205.63935594369, MedianBeforeStd: 8.3994882799157,
		SpeedBPS: 130, DwellCV: 0.6, Reverse: true,
		LevelNoiseScale: 1, WhiteNoisePA: 1.0,
		DriftPA: 1.0, DriftTau: 2.0,
		ReadShiftSD: 3.0, ReadScaleSD: 0.03,
		StallProb: 0.004, StallMin: 5, StallMax: 20, SkipProb: 0.004,
		// Direct-RNA reads begin with a long (~1–2 s) adapter/motor segment.
		AdapterSeconds: 1.5, AdapterLevelPA: 100, AdapterNoisePA: 8,
		PolyANoisePA: 1.0,
		ModelK:       9, ModelScale: 16.0, ModelShift: 79.0,
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

// MeanDwell is the mean number of samples per base, SamplingRate/SpeedBPS.
func (p Profile) MeanDwell() float64 { return p.SamplingRate / p.SpeedBPS }

// ModelOptions returns load options that convert this chemistry's ONT
// kmer_models table (standardised units) to pA and check its k.
func (p Profile) ModelOptions() []LoadOption {
	opts := []LoadOption{WithScaling(p.ModelScale, p.ModelShift)}
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
		{"level_noise_scale", p.LevelNoiseScale}, {"white_noise_pa", p.WhiteNoisePA},
		{"drift_pa", p.DriftPA}, {"read_shift_sd", p.ReadShiftSD}, {"read_scale_sd", p.ReadScaleSD},
		{"adapter_seconds", p.AdapterSeconds}, {"adapter_noise_pa", p.AdapterNoisePA}, {"polya_noise_pa", p.PolyANoisePA},
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
