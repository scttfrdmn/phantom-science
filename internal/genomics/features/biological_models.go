// Package features scores DNA sequences against simple biological heuristics:
// GC content, CpG observed/expected ratio, dinucleotide entropy, codon usage,
// k-mer distance to hard-coded primate profiles, and position-based chromatin,
// mutational-signature and selection annotations.
//
// Most models are heuristics, not implementations of the methods they are
// named after. The chromatin, selection and phylogenetic models use
// annotations synthesised from genomic position (no ENCODE, PhyloP or
// PhastCons data is loaded), and the scores they return are weighted sums of
// ad-hoc terms, not calibrated probabilities. The codon tables for human,
// yeast and E. coli come from the Kazusa Codon Usage Database. No generation
// or CLI code path uses this package.
package features

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sort"

	"github.com/scttfrdmn/phantom-science/internal/errors"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// BiologicalModelingEngine runs the enabled feature models on a sequence and
// combines their scores into one weighted heuristic score.
type BiologicalModelingEngine struct {
	gcBiasModel       *GCBiasModel
	cpgDepletionModel *CpGDepletionModel
	repeatModel       *RepetitiveElementModel
	codonUsageModel   *CodonUsageModel
	phylogeneticModel *PhylogeneticDistanceModel

	// Advanced models
	chromatinModel *ChromatinAccessibilityModel
	mutationModel  *MutationalSignatureModel
	selectionModel *SelectionPressureModel

	// Configuration
	config BiologicalModelConfig
}

// BiologicalModelConfig selects which feature models the engine runs.
// GCWindowSize and MutationalProcesses are passed to their models;
// CpGThreshold, RepeatComplexityThreshold, PhylogeneticThreshold and
// PopulationPanels have no effect on any score.
type BiologicalModelConfig struct {
	EnableGCBias            bool
	EnableCpGDepletion      bool
	EnableRepeatElements    bool
	EnableCodonUsage        bool
	EnablePhylogenetic      bool
	EnableChromatin         bool
	EnableMutationSignature bool
	EnableSelection         bool

	// Model parameters
	GCWindowSize              int
	CpGThreshold              float64
	RepeatComplexityThreshold float64
	PhylogeneticThreshold     float64

	// Population-specific parameters
	PopulationPanels    []string
	MutationalProcesses []MutationalProcess
}

// CodonUsageModel models codon usage bias in protein-coding regions
type CodonUsageModel struct {
	codonFrequencies    map[string]float64
	aminoAcidMapping    map[string][]string
	expectedFrequencies map[string]float64
	rscu                map[string]float64 // Relative Synonymous Codon Usage
	tissueSpecific      bool
	expressionLevel     float64
	organism            string
	codonTable          *GeneticCode
	rng                 *rand.Rand
}

// PhylogeneticDistanceModel compares a sequence's 4-mer frequencies with
// hard-coded primate profiles and a fixed primate tree of approximate
// divergence times. The non-human profiles are the human profile scaled by a
// constant per species, so the comparison carries no real phylogenetic signal.
type PhylogeneticDistanceModel struct {
	speciesDistances  map[string]map[string]float64
	substitutionRates map[string]float64
	branchLengths     map[string]float64
	treeTopology      *PhylogeneticTree
}

// PhylogeneticTree represents phylogenetic relationships
type PhylogeneticTree struct {
	Nodes   map[string]*TreeNode
	Root    *TreeNode
	Species []string
}

// TreeNode represents a node in the phylogenetic tree
type TreeNode struct {
	Name         string
	Children     []*TreeNode
	Parent       *TreeNode
	BranchLength float64
	Bootstrap    float64
}

// ChromatinAccessibilityModel scores a region from a synthetic chromatin
// annotation: states, histone-mark levels and transcription-factor sites are
// assigned from genomic position (chr1–chr10 only) by fixed rules, not loaded
// from ENCODE or Roadmap data.
type ChromatinAccessibilityModel struct {
	accessibilityScores map[string]float64 // genomic region -> accessibility
	histoneMarks        map[string]map[string]float64
	tfBindingSites      map[string][]TranscriptionFactorSite
	chromatinStates     map[string]ChromatinState
}

// TranscriptionFactorSite represents a transcription factor binding site
type TranscriptionFactorSite struct {
	Name       string
	Position   int64
	Chromosome string
	Score      float64
	Motif      string
}

// ChromatinState represents different chromatin states
type ChromatinState int

const (
	ActivePromoter ChromatinState = iota
	WeakPromoter
	PoisedPromoter
	StrongEnhancer
	WeakEnhancer
	Txn_Transition
	Txn_Elongation
	Weak_Txn
	Polycomb_Repressed
	Heterochromatin_low
	Repetitive_CNV
	Inactive_Low
)

// MutationalSignatureModel compares a sequence's trinucleotide composition
// with five hand-written substitution profiles named after COSMIC signatures.
type MutationalSignatureModel struct {
	signatures        map[string]*MutationalSignature
	trinucleotideFreq map[string]float64
	processActivity   map[string]float64
}

// MutationalSignature represents a specific mutational process
type MutationalSignature struct {
	Name                 string
	Process              MutationalProcess
	TrinucleotideProfile map[string]float64
	Etiology             string
	CancerTypes          []string
}

// MutationalProcess represents different mutational processes
type MutationalProcess int

const (
	SpontaneousDeamination MutationalProcess = iota
	UVExposure
	TobaccoSmoking
	DNARepairDeficiency
	ClockLikeMutation
	APOBECActivity
	ChemicalMutagenesis
	AgeRelatedMutation
)

// SelectionPressureModel scores a region from synthetic selection
// annotations: genomic context, selection coefficients, conservation scores
// and per-gene dN/dS ratios are all derived from position or index by fixed
// formulas, not from real annotations or alignments.
type SelectionPressureModel struct {
	purifyingSelection map[string]float64 // region -> selection coefficient
	positiveSelection  map[string]float64
	balancingSelection map[string]float64
	dNdSRatios         map[string]float64 // gene -> dN/dS ratio
	conservationScores map[string]float64 // region -> synthetic conservation score (not PhyloP/PhastCons)
}

// NewBiologicalModelingEngine creates a new biological modeling engine
func NewBiologicalModelingEngine(config BiologicalModelConfig) (*BiologicalModelingEngine, error) {
	engine := &BiologicalModelingEngine{
		config: config,
	}

	// Initialize enabled models
	if config.EnableGCBias {
		var err error
		engine.gcBiasModel, err = NewGCBiasModel(config.GCWindowSize)
		if err != nil {
			return nil, errors.WrapError(err, errors.BiologicalModelError, "biological_modeling_engine", "initialize_gc_bias", "failed to initialize GC bias model")
		}
	}

	if config.EnableCpGDepletion {
		engine.cpgDepletionModel = NewCpGDepletionModel()
	}

	if config.EnableRepeatElements {
		var err error
		engine.repeatModel, err = NewRepetitiveElementModel(config.RepeatComplexityThreshold)
		if err != nil {
			return nil, fmt.Errorf("failed to initialize repeat element model: %w", err)
		}
	}

	if config.EnableCodonUsage {
		var err error
		engine.codonUsageModel, err = NewCodonUsageModel()
		if err != nil {
			return nil, fmt.Errorf("failed to initialize codon usage model: %w", err)
		}
	}

	if config.EnablePhylogenetic {
		var err error
		engine.phylogeneticModel, err = NewPhylogeneticDistanceModel(config.PhylogeneticThreshold)
		if err != nil {
			return nil, fmt.Errorf("failed to initialize phylogenetic model: %w", err)
		}
	}

	if config.EnableChromatin {
		var err error
		engine.chromatinModel, err = NewChromatinAccessibilityModel()
		if err != nil {
			return nil, fmt.Errorf("failed to initialize chromatin model: %w", err)
		}
	}

	if config.EnableMutationSignature {
		var err error
		engine.mutationModel, err = NewMutationalSignatureModel(config.MutationalProcesses)
		if err != nil {
			return nil, fmt.Errorf("failed to initialize mutational signature model: %w", err)
		}
	}

	if config.EnableSelection {
		var err error
		engine.selectionModel, err = NewSelectionPressureModel()
		if err != nil {
			return nil, fmt.Errorf("failed to initialize selection pressure model: %w", err)
		}
	}

	return engine, nil
}

// AnalyzeBiologicalRealism runs each enabled model on the sequence and
// returns the per-model scores and their weighted mean. A model that returns
// an error is skipped. The "gc_bias" score is the raw GC fraction, so it is
// not a measure of realism.
func (bme *BiologicalModelingEngine) AnalyzeBiologicalRealism(sequence []byte, kmerProfile *kmer.PackedFrequencyTable,
	chromosome string, position int64) (*BiologicalRealismResult, error) {

	result := &BiologicalRealismResult{
		Chromosome:     chromosome,
		Position:       position,
		SequenceLength: int64(len(sequence)),
		FeatureScores:  make(map[string]float64),
		ModelResults:   make(map[string]interface{}),
	}

	// GC bias analysis
	if bme.gcBiasModel != nil {
		gcScore, gcDetails, err := bme.gcBiasModel.AnalyzeSequence(sequence)
		if err == nil {
			result.FeatureScores["gc_bias"] = gcScore
			result.ModelResults["gc_bias_details"] = gcDetails
		}
	}

	// CpG depletion analysis
	if bme.cpgDepletionModel != nil {
		cpgScore, cpgDetails, err := bme.cpgDepletionModel.AnalyzeSequence(sequence)
		if err == nil {
			result.FeatureScores["cpg_depletion"] = cpgScore
			result.ModelResults["cpg_depletion_details"] = cpgDetails
		}
	}

	// Repetitive element analysis
	if bme.repeatModel != nil {
		repeatScore, repeatDetails, err := bme.repeatModel.AnalyzeSequence(sequence)
		if err == nil {
			result.FeatureScores["repetitive_elements"] = repeatScore
			result.ModelResults["repeat_details"] = repeatDetails
		}
	}

	// Codon usage analysis (if sequence length allows)
	if bme.codonUsageModel != nil && len(sequence) >= 60 { // minimum for codon analysis
		codonScore, codonDetails, err := bme.codonUsageModel.AnalyzeSequence(sequence)
		if err == nil {
			result.FeatureScores["codon_usage"] = codonScore
			result.ModelResults["codon_details"] = codonDetails
		}
	}

	// Phylogenetic analysis
	if bme.phylogeneticModel != nil {
		phyloScore, phyloDetails, err := bme.phylogeneticModel.AnalyzeSequence(sequence, kmerProfile)
		if err == nil {
			result.FeatureScores["phylogenetic_distance"] = phyloScore
			result.ModelResults["phylogenetic_details"] = phyloDetails
		}
	}

	// Chromatin accessibility analysis
	if bme.chromatinModel != nil {
		chromatinScore, chromatinDetails, err := bme.chromatinModel.AnalyzeRegion(chromosome, position, position+int64(len(sequence)))
		if err == nil {
			result.FeatureScores["chromatin_accessibility"] = chromatinScore
			result.ModelResults["chromatin_details"] = chromatinDetails
		}
	}

	// Mutational signature analysis
	if bme.mutationModel != nil {
		mutationScore, mutationDetails, err := bme.mutationModel.AnalyzeSequence(sequence)
		if err == nil {
			result.FeatureScores["mutational_signature"] = mutationScore
			result.ModelResults["mutation_details"] = mutationDetails
		}
	}

	// Selection pressure analysis
	if bme.selectionModel != nil {
		selectionScore, selectionDetails, err := bme.selectionModel.AnalyzeRegion(chromosome, position, position+int64(len(sequence)))
		if err == nil {
			result.FeatureScores["selection_pressure"] = selectionScore
			result.ModelResults["selection_details"] = selectionDetails
		}
	}

	// Calculate overall biological realism score
	result.OverallScore = bme.calculateOverallScore(result.FeatureScores)
	result.RealismCategory = bme.categorizeRealism(result.OverallScore)

	return result, nil
}

// calculateOverallScore is the weighted mean of the feature scores present.
// The weights are fixed by hand, not fitted to data.
func (bme *BiologicalModelingEngine) calculateOverallScore(featureScores map[string]float64) float64 {
	if len(featureScores) == 0 {
		return 0.0
	}

	// Weighted combination of feature scores
	weights := map[string]float64{
		"gc_bias":                 0.15,
		"cpg_depletion":           0.15,
		"repetitive_elements":     0.10,
		"codon_usage":             0.12,
		"phylogenetic_distance":   0.20,
		"chromatin_accessibility": 0.10,
		"mutational_signature":    0.08,
		"selection_pressure":      0.10,
	}

	weightedSum := 0.0
	totalWeight := 0.0

	for feature, score := range featureScores {
		if weight, exists := weights[feature]; exists {
			weightedSum += score * weight
			totalWeight += weight
		}
	}

	if totalWeight == 0 {
		return 0.0
	}

	return weightedSum / totalWeight
}

// categorizeRealism categorizes the overall biological realism
func (bme *BiologicalModelingEngine) categorizeRealism(overallScore float64) RealismCategory {
	if overallScore >= 0.85 {
		return HighlyRealistic
	} else if overallScore >= 0.70 {
		return ModeratelyRealistic
	} else if overallScore >= 0.50 {
		return SomewhatRealistic
	} else {
		return PoorlyRealistic
	}
}

// BiologicalRealismResult holds the per-model scores and details from
// AnalyzeBiologicalRealism.
type BiologicalRealismResult struct {
	Chromosome      string
	Position        int64
	SequenceLength  int64
	OverallScore    float64
	RealismCategory RealismCategory
	FeatureScores   map[string]float64
	ModelResults    map[string]interface{}
}

// RealismCategory buckets OverallScore at fixed cut-offs (0.5, 0.7, 0.85).
// The cut-offs are arbitrary, not calibrated against real data.
type RealismCategory int

const (
	PoorlyRealistic RealismCategory = iota
	SomewhatRealistic
	ModeratelyRealistic
	HighlyRealistic
)

func (rc RealismCategory) String() string {
	switch rc {
	case PoorlyRealistic:
		return "Poorly Realistic"
	case SomewhatRealistic:
		return "Somewhat Realistic"
	case ModeratelyRealistic:
		return "Moderately Realistic"
	case HighlyRealistic:
		return "Highly Realistic"
	default:
		return "Unknown"
	}
}

// NewCodonUsageModel creates a new codon usage bias model
func NewCodonUsageModel() (*CodonUsageModel, error) {
	model := &CodonUsageModel{
		codonFrequencies:    make(map[string]float64),
		aminoAcidMapping:    make(map[string][]string),
		expectedFrequencies: make(map[string]float64),
		rscu:                make(map[string]float64),
		expressionLevel:     1.0,
		organism:            "human",
		rng:                 rand.New(rand.NewPCG(42, 0)), // #nosec G404
	}

	// Initialize standard genetic code mapping
	model.initializeGeneticCode()

	// Set expected codon frequencies based on human genome data
	model.initializeExpectedFrequencies()

	// Initialize codon table
	model.codonTable = getStandardGeneticCode()

	// Codon usage frequencies for the organism (used by ApplyCodonBias).
	model.loadDefaultCodonUsage()

	return model, nil
}

// initializeGeneticCode sets up the standard genetic code
func (cum *CodonUsageModel) initializeGeneticCode() {
	// Standard genetic code - amino acid to codons mapping
	geneticCode := map[string][]string{
		"A": {"GCT", "GCC", "GCA", "GCG"},               // Alanine
		"R": {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}, // Arginine
		"N": {"AAT", "AAC"},                             // Asparagine
		"D": {"GAT", "GAC"},                             // Aspartic acid
		"C": {"TGT", "TGC"},                             // Cysteine
		"Q": {"CAA", "CAG"},                             // Glutamine
		"E": {"GAA", "GAG"},                             // Glutamic acid
		"G": {"GGT", "GGC", "GGA", "GGG"},               // Glycine
		"H": {"CAT", "CAC"},                             // Histidine
		"I": {"ATT", "ATC", "ATA"},                      // Isoleucine
		"L": {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}, // Leucine
		"K": {"AAA", "AAG"},                             // Lysine
		"M": {"ATG"},                                    // Methionine (start)
		"F": {"TTT", "TTC"},                             // Phenylalanine
		"P": {"CCT", "CCC", "CCA", "CCG"},               // Proline
		"S": {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}, // Serine
		"T": {"ACT", "ACC", "ACA", "ACG"},               // Threonine
		"W": {"TGG"},                                    // Tryptophan
		"Y": {"TAT", "TAC"},                             // Tyrosine
		"V": {"GTT", "GTC", "GTA", "GTG"},               // Valine
		"*": {"TAA", "TAG", "TGA"},                      // Stop codons
	}

	cum.aminoAcidMapping = geneticCode
}

// initializeExpectedFrequencies sets expected codon frequencies from the
// Kazusa Codon Usage Database human table (Nakamura et al. 2000, NAR 28:292).
func (cum *CodonUsageModel) initializeExpectedFrequencies() {
	// Human codon usage frequencies (per thousand codons)
	humanCodonFreq := map[string]float64{
		"TTT": 17.6, "TTC": 20.3, "TTA": 7.7, "TTG": 12.9,
		"TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4,
		"TAT": 12.2, "TAC": 15.3, "TAA": 1.0, "TAG": 0.8,
		"TGT": 10.6, "TGC": 12.6, "TGA": 1.6, "TGG": 13.2,
		"CTT": 13.2, "CTC": 19.6, "CTA": 7.2, "CTG": 39.6,
		"CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
		"CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2,
		"CGT": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4,
		"ATT": 16.0, "ATC": 20.8, "ATA": 7.5, "ATG": 22.0,
		"ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
		"AAT": 17.0, "AAC": 19.1, "AAA": 24.4, "AAG": 31.9,
		"AGT": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0,
		"GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1,
		"GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
		"GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
		"GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
	}

	// Normalize to probabilities
	total := 0.0
	for _, freq := range humanCodonFreq {
		total += freq
	}

	for codon, freq := range humanCodonFreq {
		cum.expectedFrequencies[codon] = freq / total
	}
}

// AnalyzeSequence counts codons in all three forward frames (no ORF
// detection) and scores them against the human table as
// exp(−χ²/k), where χ² is summed over the k observed codons. 1.0 means the
// counts match human codon usage; the score falls towards 0 as they diverge.
func (cum *CodonUsageModel) AnalyzeSequence(sequence []byte) (float64, *CodonUsageDetails, error) {
	if len(sequence) < 3 {
		return 0.0, nil, errors.NewValidationError("codon_usage_model", "analyze_sequence", "sequence too short for codon analysis")
	}

	// Count codons in sequence
	codonCounts := make(map[string]int)
	totalCodons := 0

	// Analyze in all three reading frames
	for frame := 0; frame < 3; frame++ {
		for i := frame; i+2 < len(sequence); i += 3 {
			codon := string(sequence[i : i+3])
			// Convert to uppercase and check for valid bases
			codon = toUpper(codon)
			if isValidCodon(codon) {
				codonCounts[codon]++
				totalCodons++
			}
		}
	}

	if totalCodons == 0 {
		return 0.0, nil, errors.NewValidationError("codon_usage_model", "analyze_sequence", "no valid codons found in sequence")
	}

	// Calculate observed frequencies
	observedFreq := make(map[string]float64)
	for codon, count := range codonCounts {
		observedFreq[codon] = float64(count) / float64(totalCodons)
	}

	// Calculate codon usage bias score
	biasScore := cum.calculateCodonBias(observedFreq)

	// Calculate RSCU values
	rscu := cum.calculateRSCU(observedFreq)

	details := &CodonUsageDetails{
		TotalCodons:     totalCodons,
		ObservedFreq:    observedFreq,
		ExpectedFreq:    cum.expectedFrequencies,
		RSCU:            rscu,
		BiasScore:       biasScore,
		PreferredCodons: cum.identifyPreferredCodons(rscu),
	}

	return biasScore, details, nil
}

// calculateCodonBias calculates the overall codon usage bias score
func (cum *CodonUsageModel) calculateCodonBias(observedFreq map[string]float64) float64 {
	// Calculate chi-square-like statistic for codon usage bias
	chiSquare := 0.0
	degreesOfFreedom := 0

	for codon, observed := range observedFreq {
		if expected, exists := cum.expectedFrequencies[codon]; exists && expected > 0 {
			chiSquare += math.Pow(observed-expected, 2) / expected
			degreesOfFreedom++
		}
	}

	if degreesOfFreedom == 0 {
		return 0.0
	}

	// Normalize to 0-1 scale (lower values indicate more bias)
	normalizedScore := math.Exp(-chiSquare / float64(degreesOfFreedom))

	return normalizedScore
}

// calculateRSCU calculates Relative Synonymous Codon Usage
func (cum *CodonUsageModel) calculateRSCU(observedFreq map[string]float64) map[string]float64 {
	rscu := make(map[string]float64)

	// Calculate RSCU for each amino acid
	for aa, codons := range cum.aminoAcidMapping {
		if aa == "*" { // Skip stop codons
			continue
		}

		// Calculate mean frequency for this amino acid
		totalFreq := 0.0
		validCodons := 0

		for _, codon := range codons {
			if freq, exists := observedFreq[codon]; exists {
				totalFreq += freq
				validCodons++
			}
		}

		if validCodons == 0 || totalFreq == 0 {
			continue
		}

		meanFreq := totalFreq / float64(len(codons))

		// Calculate RSCU for each codon
		for _, codon := range codons {
			if freq, exists := observedFreq[codon]; exists {
				if meanFreq > 0 {
					rscu[codon] = freq / meanFreq
				}
			}
		}
	}

	return rscu
}

// identifyPreferredCodons identifies preferred codons (RSCU > 1.5)
func (cum *CodonUsageModel) identifyPreferredCodons(rscu map[string]float64) []string {
	var preferred []string

	for codon, rscuValue := range rscu {
		if rscuValue > 1.5 {
			preferred = append(preferred, codon)
		}
	}

	sort.Strings(preferred)
	return preferred
}

// CodonUsageDetails contains detailed codon usage analysis results
type CodonUsageDetails struct {
	TotalCodons     int
	ObservedFreq    map[string]float64
	ExpectedFreq    map[string]float64
	RSCU            map[string]float64
	BiasScore       float64
	PreferredCodons []string
}

// isValidCodon checks if a codon contains only valid DNA bases
func isValidCodon(codon string) bool {
	if len(codon) != 3 {
		return false
	}

	validBases := map[byte]bool{'A': true, 'T': true, 'G': true, 'C': true}

	for i := 0; i < 3; i++ {
		if !validBases[codon[i]] {
			return false
		}
	}

	return true
}

// Helper function to convert strings to uppercase
func toUpper(s string) string {
	result := make([]byte, len(s))
	for i, b := range []byte(s) {
		if b >= 'a' && b <= 'z' {
			result[i] = b - 32
		} else {
			result[i] = b
		}
	}
	return string(result)
}
