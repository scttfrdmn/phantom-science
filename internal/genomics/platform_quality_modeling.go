package genomics

import (
	"math"
	"math/rand/v2"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/population"
)

// SequencingPlatform represents the supported sequencing platforms
type SequencingPlatform string

const (
	NovaSeq6000  SequencingPlatform = "NovaSeq 6000"
	NextSeq2000  SequencingPlatform = "NextSeq 2000"
	MiSeq        SequencingPlatform = "MiSeq"
	PacBioRevio  SequencingPlatform = "PacBio Revio"
	DNBSeqT20x2  SequencingPlatform = "DNBSEQ-T20x2"
	ElementAVITI SequencingPlatform = "Element AVITI"
	PromethION24 SequencingPlatform = "PromethION 24"
)

// QualityProfile represents platform-specific quality characteristics
type QualityProfile struct {
	Platform                 SequencingPlatform
	Technology               string // Illumina, PacBio, ONT, MGI, Element
	ReadLength               ReadLengthSpec
	ErrorProfile             ErrorProfile
	QualityScoreDistribution QualityDistribution
	PopulationSpecificBias   map[population.Population]PlatformBias
	GCBiasCharacteristics    GCBiasProfile
	HomopolymerHandling      HomopolymerProfile
	SequenceContextEffects   ContextEffects
}

// ReadLengthSpec defines read length characteristics
type ReadLengthSpec struct {
	MinLength    int
	MaxLength    int
	MeanLength   int
	Distribution string // "fixed", "normal", "gamma"
	Parameters   map[string]float64
}

// ErrorProfile defines platform-specific error characteristics
type ErrorProfile struct {
	BaseErrorRate      float64
	IndelErrorRate     float64
	SubstitutionBias   map[string]map[string]float64 // A->T, etc.
	PositionDependency PositionErrorProfile
	QualityScoreModel  QualityScoreModel
}

// QualityDistribution defines quality score patterns
type QualityDistribution struct {
	MinQuality     int
	MaxQuality     int
	MeanQuality    float64
	StdQuality     float64
	DistributionFn string // "normal", "beta", "gamma"
	CycleDecay     CycleDecayProfile
}

// PlatformBias represents population-specific biases for a platform
type PlatformBias struct {
	GCContentBias     float64 // Bias factor for GC-rich regions
	RepeatRegionBias  float64 // Bias for repetitive sequences
	RecombinationBias float64 // Bias related to recombination hotspots
	AllelicBias       float64 // Bias toward reference vs alternate alleles
}

// GCBiasProfile defines GC content bias characteristics
type GCBiasProfile struct {
	BiasStrength float64            // Overall strength of GC bias
	OptimalGC    float64            // Optimal GC content (usually ~50%)
	BiasFunction string             // "quadratic", "sigmoid", "linear"
	RegionSizeBp int                // Size of regions for GC calculation
	Parameters   map[string]float64 // Function-specific parameters
}

// HomopolymerProfile defines homopolymer run characteristics
type HomopolymerProfile struct {
	MaxRunLength   int                // Maximum homopolymer length supported
	ErrorIncrease  map[int]float64    // Error rate increase by run length
	QualityPenalty map[int]int        // Quality score penalty by run length
	BaseSpecific   map[string]float64 // Base-specific homopolymer effects
}

// ContextEffects defines sequence context-dependent effects
type ContextEffects struct {
	MotifEffects    map[string]float64 // Effect of specific sequence motifs
	DimerBias       map[string]float64 // Dinucleotide bias effects
	SecondaryStruct map[string]float64 // Secondary structure effects
	CpGEffects      CpGProfile         // CpG-specific effects
}

// PositionErrorProfile defines position-dependent error characteristics
type PositionErrorProfile struct {
	StartErrors     []float64 // Error rates for read start positions
	EndErrors       []float64 // Error rates for read end positions
	CycleDependency bool      // Whether error depends on sequencing cycle
	DecayFunction   string    // "linear", "exponential", "sigmoid"
}

// QualityScoreModel defines how quality scores are calculated
type QualityScoreModel struct {
	BaseQuality    int     // Base quality score
	ErrorModel     string  // "phred", "solexa", "illumina"
	RandomNoise    float64 // Random noise in quality calculation
	SystematicBias float64 // Systematic bias in quality scores
}

// CycleDecayProfile defines quality decay over sequencing cycles
type CycleDecayProfile struct {
	DecayRate     float64 // Rate of quality decay per cycle
	DecayFunction string  // "linear", "exponential", "logarithmic"
	MinQuality    int     // Minimum quality floor
	PairEndDiff   float64 // Difference between R1 and R2 quality
}

// CpGProfile defines CpG-specific sequencing effects
type CpGProfile struct {
	CpGBias        float64 // Bias in CpG sites
	MethylationSim bool    // Simulate methylation effects
	ConversionRate float64 // Bisulfite conversion rate (if applicable)
}

// MultiPlatformQualityModeler holds per-platform quality profiles and
// per-population parameters, and generates heuristic quality strings from
// them. It is experimental: the FASTQ generation pipeline does not use it
// (only the cmd/experiments programs do). The profile values are hand-set
// approximations without cited sources, and its RNG is seeded from the
// clock, so results are not reproducible.
type MultiPlatformQualityModeler struct {
	platformProfiles map[SequencingPlatform]*QualityProfile
	populationData   map[population.Population]*PopulationQualityData
	rng              *rand.Rand
}

// PopulationQualityData stores per-population parameters. The values are
// placeholders, not derived from 1000 Genomes data.
type PopulationQualityData struct {
	Population     population.Population
	GCContentMean  float64
	GCContentStd   float64
	RepeatContent  float64
	RecombHotspots []int // Positions of recombination hotspots
	AllelicBalance float64
}

// NewMultiPlatformQualityModeler creates a new multi-platform quality modeler
func NewMultiPlatformQualityModeler() *MultiPlatformQualityModeler {
	modeler := &MultiPlatformQualityModeler{
		platformProfiles: make(map[SequencingPlatform]*QualityProfile),
		populationData:   make(map[population.Population]*PopulationQualityData),
		rng:              rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
	}

	// Initialize platform profiles
	modeler.initializePlatformProfiles()

	// Initialize population-specific data
	modeler.initializePopulationData()

	return modeler
}

// initializePlatformProfiles sets up quality profiles for all 7 platforms
func (m *MultiPlatformQualityModeler) initializePlatformProfiles() {
	// Illumina NovaSeq 6000 - High throughput, short reads
	m.platformProfiles[NovaSeq6000] = &QualityProfile{
		Platform:   NovaSeq6000,
		Technology: "Illumina",
		ReadLength: ReadLengthSpec{
			MinLength:    75,
			MaxLength:    150,
			MeanLength:   150,
			Distribution: "fixed",
			Parameters:   map[string]float64{"length": 150},
		},
		ErrorProfile: ErrorProfile{
			BaseErrorRate:  0.001, // 0.1% error rate
			IndelErrorRate: 0.0001,
			SubstitutionBias: map[string]map[string]float64{
				"A": {"T": 0.4, "G": 0.3, "C": 0.3},
				"T": {"A": 0.4, "G": 0.3, "C": 0.3},
				"G": {"C": 0.4, "A": 0.3, "T": 0.3},
				"C": {"G": 0.4, "A": 0.3, "T": 0.3},
			},
			PositionDependency: PositionErrorProfile{
				CycleDependency: true,
				DecayFunction:   "linear",
			},
		},
		QualityScoreDistribution: QualityDistribution{
			MinQuality:     2,
			MaxQuality:     41,
			MeanQuality:    35.0,
			StdQuality:     5.0,
			DistributionFn: "normal",
			CycleDecay: CycleDecayProfile{
				DecayRate:     0.1,
				DecayFunction: "linear",
				MinQuality:    10,
				PairEndDiff:   2.0,
			},
		},
		GCBiasCharacteristics: GCBiasProfile{
			BiasStrength: 0.15,
			OptimalGC:    0.50,
			BiasFunction: "quadratic",
			RegionSizeBp: 100,
			Parameters:   map[string]float64{"curvature": 2.0},
		},
		HomopolymerHandling: HomopolymerProfile{
			MaxRunLength:   20,
			ErrorIncrease:  map[int]float64{3: 1.5, 4: 2.0, 5: 3.0, 6: 4.0},
			QualityPenalty: map[int]int{3: 2, 4: 4, 5: 6, 6: 8},
			BaseSpecific:   map[string]float64{"A": 1.0, "T": 1.0, "G": 1.2, "C": 1.2},
		},
	}

	// Illumina NextSeq 2000 - Mid throughput
	m.platformProfiles[NextSeq2000] = &QualityProfile{
		Platform:   NextSeq2000,
		Technology: "Illumina",
		ReadLength: ReadLengthSpec{
			MinLength:    75,
			MaxLength:    300,
			MeanLength:   150,
			Distribution: "fixed",
			Parameters:   map[string]float64{"length": 150},
		},
		ErrorProfile: ErrorProfile{
			BaseErrorRate:  0.0008, // Slightly better than NovaSeq
			IndelErrorRate: 0.00008,
			SubstitutionBias: map[string]map[string]float64{
				"A": {"T": 0.35, "G": 0.35, "C": 0.30},
				"T": {"A": 0.35, "G": 0.35, "C": 0.30},
				"G": {"C": 0.40, "A": 0.30, "T": 0.30},
				"C": {"G": 0.40, "A": 0.30, "T": 0.30},
			},
		},
		QualityScoreDistribution: QualityDistribution{
			MinQuality:     2,
			MaxQuality:     41,
			MeanQuality:    36.0,
			StdQuality:     4.5,
			DistributionFn: "normal",
		},
		GCBiasCharacteristics: GCBiasProfile{
			BiasStrength: 0.12,
			OptimalGC:    0.50,
			BiasFunction: "quadratic",
			RegionSizeBp: 100,
		},
	}

	// Illumina MiSeq - Lower throughput, higher quality
	m.platformProfiles[MiSeq] = &QualityProfile{
		Platform:   MiSeq,
		Technology: "Illumina",
		ReadLength: ReadLengthSpec{
			MinLength:    75,
			MaxLength:    300,
			MeanLength:   250,
			Distribution: "fixed",
			Parameters:   map[string]float64{"length": 250},
		},
		ErrorProfile: ErrorProfile{
			BaseErrorRate:  0.0005, // Higher quality
			IndelErrorRate: 0.00005,
		},
		QualityScoreDistribution: QualityDistribution{
			MinQuality:     2,
			MaxQuality:     41,
			MeanQuality:    38.0,
			StdQuality:     3.5,
			DistributionFn: "normal",
		},
		GCBiasCharacteristics: GCBiasProfile{
			BiasStrength: 0.08,
			OptimalGC:    0.50,
			BiasFunction: "quadratic",
			RegionSizeBp: 100,
		},
	}

	// PacBio Revio - Long reads, higher error rate
	m.platformProfiles[PacBioRevio] = &QualityProfile{
		Platform:   PacBioRevio,
		Technology: "PacBio",
		ReadLength: ReadLengthSpec{
			MinLength:    10000,
			MaxLength:    100000,
			MeanLength:   25000,
			Distribution: "gamma",
			Parameters:   map[string]float64{"shape": 2.0, "scale": 12500},
		},
		ErrorProfile: ErrorProfile{
			BaseErrorRate:  0.01,  // 1% error rate
			IndelErrorRate: 0.008, // Higher indel rate
			SubstitutionBias: map[string]map[string]float64{
				"A": {"T": 0.33, "G": 0.33, "C": 0.34},
				"T": {"A": 0.33, "G": 0.33, "C": 0.34},
				"G": {"C": 0.33, "A": 0.33, "T": 0.34},
				"C": {"G": 0.33, "A": 0.33, "T": 0.34},
			},
		},
		QualityScoreDistribution: QualityDistribution{
			MinQuality:     2,
			MaxQuality:     93,
			MeanQuality:    25.0,
			StdQuality:     8.0,
			DistributionFn: "gamma",
		},
		HomopolymerHandling: HomopolymerProfile{
			MaxRunLength:   50,
			ErrorIncrease:  map[int]float64{5: 2.0, 10: 4.0, 15: 8.0, 20: 15.0},
			QualityPenalty: map[int]int{5: 5, 10: 10, 15: 15, 20: 20},
		},
	}

	// MGI DNBSEQ-T20x2 - Chinese platform
	m.platformProfiles[DNBSeqT20x2] = &QualityProfile{
		Platform:   DNBSeqT20x2,
		Technology: "MGI",
		ReadLength: ReadLengthSpec{
			MinLength:    50,
			MaxLength:    150,
			MeanLength:   150,
			Distribution: "fixed",
			Parameters:   map[string]float64{"length": 150},
		},
		ErrorProfile: ErrorProfile{
			BaseErrorRate:  0.0012, // Similar to Illumina
			IndelErrorRate: 0.0001,
		},
		QualityScoreDistribution: QualityDistribution{
			MinQuality:     2,
			MaxQuality:     40,
			MeanQuality:    33.0,
			StdQuality:     5.5,
			DistributionFn: "normal",
		},
		GCBiasCharacteristics: GCBiasProfile{
			BiasStrength: 0.18, // Slightly higher GC bias
			OptimalGC:    0.48,
			BiasFunction: "quadratic",
			RegionSizeBp: 100,
		},
	}

	// Element AVITI - Novel platform
	m.platformProfiles[ElementAVITI] = &QualityProfile{
		Platform:   ElementAVITI,
		Technology: "Element",
		ReadLength: ReadLengthSpec{
			MinLength:    75,
			MaxLength:    300,
			MeanLength:   150,
			Distribution: "fixed",
			Parameters:   map[string]float64{"length": 150},
		},
		ErrorProfile: ErrorProfile{
			BaseErrorRate:  0.0015, // New platform, potentially higher error
			IndelErrorRate: 0.0002,
		},
		QualityScoreDistribution: QualityDistribution{
			MinQuality:     2,
			MaxQuality:     40,
			MeanQuality:    30.0,
			StdQuality:     6.0,
			DistributionFn: "normal",
		},
		GCBiasCharacteristics: GCBiasProfile{
			BiasStrength: 0.20,
			OptimalGC:    0.52,
			BiasFunction: "sigmoid",
			RegionSizeBp: 150,
		},
	}

	// Oxford Nanopore PromethION 24 - Long reads, higher error
	m.platformProfiles[PromethION24] = &QualityProfile{
		Platform:   PromethION24,
		Technology: "ONT",
		ReadLength: ReadLengthSpec{
			MinLength:    1000,
			MaxLength:    200000,
			MeanLength:   15000,
			Distribution: "gamma",
			Parameters:   map[string]float64{"shape": 1.5, "scale": 10000},
		},
		ErrorProfile: ErrorProfile{
			BaseErrorRate:  0.03,  // 3% error rate
			IndelErrorRate: 0.025, // Very high indel rate
			SubstitutionBias: map[string]map[string]float64{
				"A": {"T": 0.5, "G": 0.25, "C": 0.25},
				"T": {"A": 0.5, "G": 0.25, "C": 0.25},
				"G": {"C": 0.4, "A": 0.3, "T": 0.3},
				"C": {"G": 0.4, "A": 0.3, "T": 0.3},
			},
		},
		QualityScoreDistribution: QualityDistribution{
			MinQuality:     2,
			MaxQuality:     60,
			MeanQuality:    15.0,
			StdQuality:     7.0,
			DistributionFn: "gamma",
		},
		HomopolymerHandling: HomopolymerProfile{
			MaxRunLength:   100,
			ErrorIncrease:  map[int]float64{3: 3.0, 5: 6.0, 10: 12.0, 15: 25.0},
			QualityPenalty: map[int]int{3: 8, 5: 15, 10: 25, 15: 35},
		},
	}
}

// initializePopulationData sets up population-specific quality modeling data
func (m *MultiPlatformQualityModeler) initializePopulationData() {
	// Placeholder values; not derived from literature or 1000 Genomes data.
	populations := []population.Population{
		population.YRI, population.LWK, population.GWD, // African
		population.CEU, population.GBR, population.TSI, // European
		population.CHB, population.JPT, population.CHS, // East Asian
		population.MXL, population.PUR, // Americas
	}

	for _, pop := range populations {
		m.populationData[pop] = &PopulationQualityData{
			Population:     pop,
			GCContentMean:  m.getPopulationGCContent(pop),
			GCContentStd:   m.getPopulationGCVariation(pop),
			RepeatContent:  m.getPopulationRepeatContent(pop),
			AllelicBalance: m.getPopulationAllelicBalance(pop),
		}
	}
}

// getPopulationGCContent returns a per-population GC content. The small
// differences between populations are placeholders with no empirical basis:
// genome-wide GC content does not differ measurably by ancestry.
func (m *MultiPlatformQualityModeler) getPopulationGCContent(pop population.Population) float64 {
	gcContent := map[population.Population]float64{
		population.YRI: 0.406,
		population.LWK: 0.408,
		population.GWD: 0.407,
		population.CEU: 0.410,
		population.GBR: 0.409,
		population.TSI: 0.411,
		population.CHB: 0.412,
		population.JPT: 0.413,
		population.CHS: 0.412,
		population.MXL: 0.409,
		population.PUR: 0.410,
	}

	if gc, exists := gcContent[pop]; exists {
		return gc
	}
	return 0.410 // Default GC content
}

// getPopulationGCVariation returns a per-population GC content spread.
// The values are placeholders with no empirical basis.
func (m *MultiPlatformQualityModeler) getPopulationGCVariation(pop population.Population) float64 {
	if pop == population.YRI || pop == population.LWK || pop == population.GWD {
		return 0.025
	}
	if pop == population.CHB || pop == population.JPT || pop == population.CHS {
		return 0.018
	}
	return 0.022
}

// getPopulationRepeatContent returns a per-population repeat fraction.
// The values are placeholders with no empirical basis.
func (m *MultiPlatformQualityModeler) getPopulationRepeatContent(pop population.Population) float64 {
	if pop == population.YRI || pop == population.LWK || pop == population.GWD {
		return 0.42
	}
	return 0.45
}

// getPopulationAllelicBalance returns a per-population allelic balance.
// The 0.48 for MXL and PUR is a placeholder with no empirical basis.
func (m *MultiPlatformQualityModeler) getPopulationAllelicBalance(pop population.Population) float64 {
	if pop == population.MXL || pop == population.PUR {
		return 0.48
	}
	return 0.50
}

// GenerateQualityScores returns heuristic Phred scores for sequence: the
// platform's mean quality with cycle decay, a GC factor, a homopolymer
// penalty and Gaussian noise, clamped to the platform's range. Unknown
// platforms get Q30 +/- 5 noise.
func (m *MultiPlatformQualityModeler) GenerateQualityScores(
	platform SequencingPlatform,
	pop population.Population,
	sequence string,
	readPosition int,
) []int {
	profile := m.platformProfiles[platform]
	if profile == nil {
		return m.generateDefaultQualityScores(len(sequence))
	}

	popData := m.populationData[pop]
	qualities := make([]int, len(sequence))

	for i, base := range sequence {
		baseQuality := m.calculateBaseQuality(profile, popData, string(base), i, readPosition, sequence)
		qualities[i] = int(math.Max(float64(profile.QualityScoreDistribution.MinQuality),
			math.Min(float64(profile.QualityScoreDistribution.MaxQuality), baseQuality)))
	}

	return qualities
}

// calculateBaseQuality calculates quality score for a specific base
func (m *MultiPlatformQualityModeler) calculateBaseQuality(
	profile *QualityProfile,
	popData *PopulationQualityData,
	base string,
	position int,
	readPosition int,
	sequence string,
) float64 {
	// Start with base quality
	baseQuality := profile.QualityScoreDistribution.MeanQuality

	// Apply cycle decay if applicable
	if profile.QualityScoreDistribution.CycleDecay.DecayRate > 0 {
		decay := profile.QualityScoreDistribution.CycleDecay.DecayRate * float64(position)
		baseQuality -= decay
	}

	// Apply GC bias
	gcBias := m.calculateGCBias(profile, popData, sequence, position)
	baseQuality *= gcBias

	// Apply homopolymer penalty
	homopolymerPenalty := m.calculateHomopolymerPenalty(profile, sequence, position)
	baseQuality -= float64(homopolymerPenalty)

	// Add random noise
	noise := m.rng.NormFloat64() * profile.QualityScoreDistribution.StdQuality
	baseQuality += noise

	return baseQuality
}

// calculateGCBias calculates GC bias effect on quality
func (m *MultiPlatformQualityModeler) calculateGCBias(
	profile *QualityProfile,
	popData *PopulationQualityData,
	sequence string,
	position int,
) float64 {
	regionSize := profile.GCBiasCharacteristics.RegionSizeBp
	start := int(math.Max(0, float64(position-regionSize/2)))
	end := int(math.Min(float64(len(sequence)), float64(position+regionSize/2)))

	region := sequence[start:end]
	gcContent := m.calculateGCContent(region)

	optimalGC := profile.GCBiasCharacteristics.OptimalGC
	biasStrength := profile.GCBiasCharacteristics.BiasStrength

	// Quadratic bias function
	deviation := gcContent - optimalGC
	bias := 1.0 - biasStrength*deviation*deviation

	return math.Max(0.5, math.Min(1.5, bias))
}

// calculateGCContent calculates GC content of a sequence
func (m *MultiPlatformQualityModeler) calculateGCContent(sequence string) float64 {
	if len(sequence) == 0 {
		return 0.5
	}

	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}

	return float64(gcCount) / float64(len(sequence))
}

// calculateHomopolymerPenalty calculates quality penalty for homopolymer runs
func (m *MultiPlatformQualityModeler) calculateHomopolymerPenalty(
	profile *QualityProfile,
	sequence string,
	position int,
) int {
	if position >= len(sequence) {
		return 0
	}

	base := sequence[position]
	runLength := 1

	// Count backward
	for i := position - 1; i >= 0 && sequence[i] == base; i-- {
		runLength++
	}

	// Count forward
	for i := position + 1; i < len(sequence) && sequence[i] == base; i++ {
		runLength++
	}

	// Apply penalty based on run length
	if penalty, exists := profile.HomopolymerHandling.QualityPenalty[runLength]; exists {
		return penalty
	}

	// Linear penalty for longer runs
	if runLength > 6 {
		return runLength * 2
	}

	return 0
}

// generateDefaultQualityScores generates default quality scores when no profile exists
func (m *MultiPlatformQualityModeler) generateDefaultQualityScores(length int) []int {
	qualities := make([]int, length)
	for i := range qualities {
		qualities[i] = 30 + int(m.rng.NormFloat64()*5) // Q30 ± 5
		if qualities[i] < 2 {
			qualities[i] = 2
		}
		if qualities[i] > 40 {
			qualities[i] = 40
		}
	}
	return qualities
}

// ValidateMultiPlatformQuality checks that the profile table is internally
// consistent (all seven platforms present, Illumina base error rate at most
// 0.2%, long-read mean length at least 10 kb). It does not compare against
// real data, and "quality_generation_functional" is always true.
func (m *MultiPlatformQualityModeler) ValidateMultiPlatformQuality() map[string]bool {
	validation := make(map[string]bool)

	// Test that all 7 platforms have profiles
	expectedPlatforms := []SequencingPlatform{
		NovaSeq6000, NextSeq2000, MiSeq, PacBioRevio,
		DNBSeqT20x2, ElementAVITI, PromethION24,
	}

	platformsConfigured := 0
	for _, platform := range expectedPlatforms {
		if _, exists := m.platformProfiles[platform]; exists {
			platformsConfigured++
		}
	}

	validation["all_platforms_configured"] = platformsConfigured == len(expectedPlatforms)
	validation["population_data_loaded"] = len(m.populationData) >= 10
	validation["quality_generation_functional"] = true // Basic functionality test

	// Test platform-specific characteristics
	illuminaPlatforms := []SequencingPlatform{NovaSeq6000, NextSeq2000, MiSeq}
	longReadPlatforms := []SequencingPlatform{PacBioRevio, PromethION24}

	illuminaProperlyConfigured := true
	for _, platform := range illuminaPlatforms {
		if profile := m.platformProfiles[platform]; profile != nil {
			if profile.Technology != "Illumina" || profile.ErrorProfile.BaseErrorRate > 0.002 {
				illuminaProperlyConfigured = false
				break
			}
		}
	}

	longReadProperlyConfigured := true
	for _, platform := range longReadPlatforms {
		if profile := m.platformProfiles[platform]; profile != nil {
			if profile.ReadLength.MeanLength < 10000 {
				longReadProperlyConfigured = false
				break
			}
		}
	}

	validation["illumina_platforms_configured"] = illuminaProperlyConfigured
	validation["long_read_platforms_configured"] = longReadProperlyConfigured

	return validation
}

// GetSupportedPlatforms returns list of all supported platforms
func (m *MultiPlatformQualityModeler) GetSupportedPlatforms() []SequencingPlatform {
	return []SequencingPlatform{
		NovaSeq6000, NextSeq2000, MiSeq, PacBioRevio,
		DNBSeqT20x2, ElementAVITI, PromethION24,
	}
}

// GetPlatformProfile returns the quality profile for a specific platform
func (m *MultiPlatformQualityModeler) GetPlatformProfile(platform SequencingPlatform) *QualityProfile {
	return m.platformProfiles[platform]
}

// GetPopulationData returns population-specific quality data
func (m *MultiPlatformQualityModeler) GetPopulationData(pop population.Population) *PopulationQualityData {
	return m.populationData[pop]
}
