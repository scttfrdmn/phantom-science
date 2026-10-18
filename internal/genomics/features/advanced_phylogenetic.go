package features

import (
	"fmt"
	"math"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// NewPhylogeneticDistanceModel creates a phylogenetic distance model with the
// built-in primate table. threshold is ignored.
func NewPhylogeneticDistanceModel(threshold float64) (*PhylogeneticDistanceModel, error) {
	model := &PhylogeneticDistanceModel{
		speciesDistances:  make(map[string]map[string]float64),
		substitutionRates: make(map[string]float64),
		branchLengths:     make(map[string]float64),
	}

	// Initialize with primate phylogenetic data
	if err := model.initializePrimatePhylogeny(); err != nil {
		return nil, fmt.Errorf("failed to initialize primate phylogeny: %w", err)
	}

	return model, nil
}

// initializePrimatePhylogeny sets up phylogenetic relationships for primates
func (pdm *PhylogeneticDistanceModel) initializePrimatePhylogeny() error {
	// Species in the phylogeny
	species := []string{"Human", "Chimpanzee", "Bonobo", "Gorilla", "Orangutan", "Macaque", "Marmoset"}

	// Distance matrix of rounded, approximate divergence times in millions of
	// years (MYA); not taken from a specific published tree
	distanceMatrix := map[string]map[string]float64{
		"Human": {
			"Human":      0.0,
			"Chimpanzee": 6.0,
			"Bonobo":     6.0,
			"Gorilla":    8.0,
			"Orangutan":  14.0,
			"Macaque":    25.0,
			"Marmoset":   40.0,
		},
		"Chimpanzee": {
			"Human":      6.0,
			"Chimpanzee": 0.0,
			"Bonobo":     1.0,
			"Gorilla":    8.0,
			"Orangutan":  14.0,
			"Macaque":    25.0,
			"Marmoset":   40.0,
		},
		"Bonobo": {
			"Human":      6.0,
			"Chimpanzee": 1.0,
			"Bonobo":     0.0,
			"Gorilla":    8.0,
			"Orangutan":  14.0,
			"Macaque":    25.0,
			"Marmoset":   40.0,
		},
		"Gorilla": {
			"Human":      8.0,
			"Chimpanzee": 8.0,
			"Bonobo":     8.0,
			"Gorilla":    0.0,
			"Orangutan":  14.0,
			"Macaque":    25.0,
			"Marmoset":   40.0,
		},
		"Orangutan": {
			"Human":      14.0,
			"Chimpanzee": 14.0,
			"Bonobo":     14.0,
			"Gorilla":    14.0,
			"Orangutan":  0.0,
			"Macaque":    25.0,
			"Marmoset":   40.0,
		},
		"Macaque": {
			"Human":      25.0,
			"Chimpanzee": 25.0,
			"Bonobo":     25.0,
			"Gorilla":    25.0,
			"Orangutan":  25.0,
			"Macaque":    0.0,
			"Marmoset":   40.0,
		},
		"Marmoset": {
			"Human":      40.0,
			"Chimpanzee": 40.0,
			"Bonobo":     40.0,
			"Gorilla":    40.0,
			"Orangutan":  40.0,
			"Macaque":    40.0,
			"Marmoset":   0.0,
		},
	}

	pdm.speciesDistances = distanceMatrix

	// Placeholder substitution rate, the same for every species; not used by
	// any score
	for _, sp := range species {
		pdm.substitutionRates[sp] = 1.0e-8 // Base rate
	}

	// Branch lengths: divergence time (MYA) × 1e-8, an arbitrary scale
	for _, sp := range species {
		pdm.branchLengths[sp] = distanceMatrix["Human"][sp] * 1.0e-8
	}

	// Build tree topology
	pdm.treeTopology = pdm.buildPrimateTree()

	return nil
}

// buildPrimateTree constructs the primate phylogenetic tree
func (pdm *PhylogeneticDistanceModel) buildPrimateTree() *PhylogeneticTree {
	tree := &PhylogeneticTree{
		Nodes:   make(map[string]*TreeNode),
		Species: []string{"Human", "Chimpanzee", "Bonobo", "Gorilla", "Orangutan", "Macaque", "Marmoset"},
	}

	// Create leaf nodes for each species
	for _, species := range tree.Species {
		node := &TreeNode{
			Name:         species,
			BranchLength: pdm.branchLengths[species],
			Bootstrap:    100.0, // Fixed value, not computed
		}
		tree.Nodes[species] = node
	}

	// Internal nodes: a fixed topology with hand-set branch lengths and
	// bootstrap values; no tree is reconstructed from data

	// Human-Chimp-Bonobo clade
	humanChimpClade := &TreeNode{
		Name:         "Human-Chimp-Bonobo",
		Children:     []*TreeNode{tree.Nodes["Human"]},
		BranchLength: 6.0e-8,
		Bootstrap:    95.0,
	}

	chimpBonoboClade := &TreeNode{
		Name:         "Chimp-Bonobo",
		Children:     []*TreeNode{tree.Nodes["Chimpanzee"], tree.Nodes["Bonobo"]},
		BranchLength: 1.0e-8,
		Bootstrap:    98.0,
	}

	humanChimpClade.Children = append(humanChimpClade.Children, chimpBonoboClade)

	// Set parent relationships
	tree.Nodes["Human"].Parent = humanChimpClade
	chimpBonoboClade.Parent = humanChimpClade
	tree.Nodes["Chimpanzee"].Parent = chimpBonoboClade
	tree.Nodes["Bonobo"].Parent = chimpBonoboClade

	// Great ape clade (including Gorilla)
	greatApeClade := &TreeNode{
		Name:         "Great-Apes",
		Children:     []*TreeNode{humanChimpClade, tree.Nodes["Gorilla"]},
		BranchLength: 8.0e-8,
		Bootstrap:    92.0,
	}

	humanChimpClade.Parent = greatApeClade
	tree.Nodes["Gorilla"].Parent = greatApeClade

	// Hominidae (including Orangutan)
	hominidaeClade := &TreeNode{
		Name:         "Hominidae",
		Children:     []*TreeNode{greatApeClade, tree.Nodes["Orangutan"]},
		BranchLength: 14.0e-8,
		Bootstrap:    88.0,
	}

	greatApeClade.Parent = hominidaeClade
	tree.Nodes["Orangutan"].Parent = hominidaeClade

	// Catarrhini (including Macaque)
	catarrhiniClade := &TreeNode{
		Name:         "Catarrhini",
		Children:     []*TreeNode{hominidaeClade, tree.Nodes["Macaque"]},
		BranchLength: 25.0e-8,
		Bootstrap:    85.0,
	}

	hominidaeClade.Parent = catarrhiniClade
	tree.Nodes["Macaque"].Parent = catarrhiniClade

	// Root (including Marmoset as outgroup)
	root := &TreeNode{
		Name:         "Primates",
		Children:     []*TreeNode{catarrhiniClade, tree.Nodes["Marmoset"]},
		BranchLength: 0.0,
		Bootstrap:    100.0,
	}

	catarrhiniClade.Parent = root
	tree.Nodes["Marmoset"].Parent = root
	tree.Root = root

	return tree
}

// AnalyzeSequence scores a sequence with heuristics: the RMS distance between
// kmerProfile and each hard-coded species 4-mer profile, the correlation of
// those distances with divergence time, GC deviation from 0.41 (labelled a
// conservation score), the transition/transversion ratio of adjacent bases,
// and a Gaussian fit of distances to divergence times. No phylogenetic
// inference is done. Sequences shorter than 100 bp return an error.
func (pdm *PhylogeneticDistanceModel) AnalyzeSequence(sequence []byte, kmerProfile *kmer.PackedFrequencyTable) (float64, *PhylogeneticAnalysisDetails, error) {
	if len(sequence) < 100 {
		return 0.0, nil, fmt.Errorf("sequence too short for phylogenetic analysis")
	}

	// Calculate evolutionary distances based on k-mer composition
	distances := pdm.calculateEvolutionaryDistances(kmerProfile)

	// Estimate phylogenetic signal strength
	phyloSignal := pdm.calculatePhylogeneticSignal(distances)

	// Calculate conservation score
	conservationScore := pdm.calculateConservationScore(sequence)

	// Identify most likely evolutionary origin
	likelyOrigin := pdm.identifyEvolutionaryOrigin(distances)

	// Calculate substitution pattern consistency
	substitutionConsistency := pdm.analyzeSubstitutionPatterns(sequence)

	details := &PhylogeneticAnalysisDetails{
		EvolutionaryDistances:    distances,
		PhylogeneticSignal:       phyloSignal,
		ConservationScore:        conservationScore,
		LikelyEvolutionaryOrigin: likelyOrigin,
		SubstitutionConsistency:  substitutionConsistency,
		TreeLikelihood:           pdm.calculateTreeLikelihood(distances),
	}

	// Overall phylogenetic realism score
	realismScore := pdm.calculatePhylogeneticRealism(details)

	return realismScore, details, nil
}

// calculateEvolutionaryDistances estimates evolutionary distances based on k-mer composition
func (pdm *PhylogeneticDistanceModel) calculateEvolutionaryDistances(kmerProfile *kmer.PackedFrequencyTable) map[string]float64 {
	distances := make(map[string]float64)

	// Hard-coded k-mer profiles (see getSpeciesKmerProfiles)
	speciesKmerProfiles := pdm.getSpeciesKmerProfiles()

	for species, expectedProfile := range speciesKmerProfiles {
		distance := pdm.calculateKmerDistance(kmerProfile, expectedProfile)
		distances[species] = distance
	}

	return distances
}

// getSpeciesKmerProfiles returns placeholder 4-mer profiles. The human profile
// is 16 made-up values for the AA-prefixed 4-mers, not derived from genomic
// data. Every other species is the human profile times a constant: the
// "variation" terms depend only on len(kmer), which is always 4.
func (pdm *PhylogeneticDistanceModel) getSpeciesKmerProfiles() map[string]map[string]float64 {
	profiles := make(map[string]map[string]float64)

	// Human profile (baseline)
	humanProfile := map[string]float64{
		"AAAA": 0.062, "AAAT": 0.049, "AAAG": 0.058, "AAAC": 0.051,
		"AATA": 0.041, "AATT": 0.042, "AATG": 0.047, "AATC": 0.039,
		"AAGA": 0.055, "AAGT": 0.048, "AAGG": 0.062, "AAGC": 0.053,
		"AACA": 0.047, "AACT": 0.051, "AACG": 0.043, "AACC": 0.058,
	}
	profiles["Human"] = humanProfile

	// Chimpanzee: human profile times a constant
	chimpProfile := make(map[string]float64)
	for kmer, freq := range humanProfile {
		// Constant factor (len(kmer) is always 4), not random
		variation := 1.0 + (0.02 * (math.Sin(float64(len(kmer))) - 0.5))
		chimpProfile[kmer] = freq * variation
	}
	profiles["Chimpanzee"] = chimpProfile

	// Other species with increasing divergence
	for _, species := range []string{"Bonobo", "Gorilla", "Orangutan", "Macaque", "Marmoset"} {
		profile := make(map[string]float64)
		divergence := pdm.speciesDistances["Human"][species] / 40.0 // Normalize by max distance

		for kmer, freq := range humanProfile {
			// Increase variation with evolutionary distance
			variation := 1.0 + (divergence * 0.1 * (math.Cos(float64(len(kmer))) - 0.5))
			profile[kmer] = freq * variation
		}
		profiles[species] = profile
	}

	return profiles
}

// calculateKmerDistance calculates distance between observed and expected k-mer profiles
func (pdm *PhylogeneticDistanceModel) calculateKmerDistance(observed *kmer.PackedFrequencyTable, expected map[string]float64) float64 {
	if observed == nil {
		return math.Inf(1)
	}

	distance := 0.0
	totalComparisons := 0

	// Convert packed array to frequency map for comparison
	observedFreqs := observed.ToFrequencyTable()

	for kmer, expectedFreq := range expected {
		observedFreq := 0.0
		if freq, exists := observedFreqs.Frequencies[kmer]; exists {
			observedFreq = float64(freq)
		}

		// Calculate squared difference
		diff := observedFreq - expectedFreq
		distance += diff * diff
		totalComparisons++
	}

	if totalComparisons == 0 {
		return math.Inf(1)
	}

	// Return root mean square distance
	return math.Sqrt(distance / float64(totalComparisons))
}

// calculatePhylogeneticSignal estimates the strength of phylogenetic signal
func (pdm *PhylogeneticDistanceModel) calculatePhylogeneticSignal(distances map[string]float64) float64 {
	if len(distances) < 2 {
		return 0.0
	}

	// Calculate correlation between observed distances and expected phylogenetic distances
	var observedDist, expectedDist []float64

	for species, obsDist := range distances {
		if expDist, exists := pdm.speciesDistances["Human"][species]; exists {
			observedDist = append(observedDist, obsDist)
			expectedDist = append(expectedDist, expDist)
		}
	}

	if len(observedDist) < 2 {
		return 0.0
	}

	// Calculate Pearson correlation
	correlation := calculateCorrelation(observedDist, expectedDist)

	// Return absolute correlation as signal strength
	if correlation < 0 {
		return -correlation
	}
	return correlation
}

// calculateConservationScore returns exp(−5·|GC − 0.41|). It is a GC-content
// heuristic, not a measure of conservation.
func (pdm *PhylogeneticDistanceModel) calculateConservationScore(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0.0
	}

	// GC content stands in for conservation
	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' || base == 'g' || base == 'c' {
			gcCount++
		}
	}

	gcContent := float64(gcCount) / float64(len(sequence))

	// Human genome average GC content is ~41% (Lander et al. 2001)
	expectedGC := 0.41
	gcDeviation := math.Abs(gcContent - expectedGC)

	// Convert to conservation score (0-1, higher = more conserved)
	conservationScore := math.Exp(-5.0 * gcDeviation)

	return conservationScore
}

// identifyEvolutionaryOrigin identifies the most likely evolutionary origin
func (pdm *PhylogeneticDistanceModel) identifyEvolutionaryOrigin(distances map[string]float64) string {
	minDistance := math.Inf(1)
	closestSpecies := ""

	for species, distance := range distances {
		if distance < minDistance {
			minDistance = distance
			closestSpecies = species
		}
	}

	return closestSpecies
}

// analyzeSubstitutionPatterns scores the ratio of transition-type to
// transversion-type adjacent base pairs against 2.0. Adjacent bases in one
// sequence are not substitutions, so this is a composition heuristic, not a
// Ts/Tv estimate.
func (pdm *PhylogeneticDistanceModel) analyzeSubstitutionPatterns(sequence []byte) float64 {
	if len(sequence) < 10 {
		return 0.0
	}

	// Count transition and transversion substitutions
	transitions := 0   // A<->G, C<->T
	transversions := 0 // A<->C, A<->T, G<->C, G<->T

	// Analyze adjacent base pairs as a proxy for substitution patterns
	for i := 0; i < len(sequence)-1; i++ {
		base1 := sequence[i]
		base2 := sequence[i+1]

		// Check for transition patterns
		if (base1 == 'A' && base2 == 'G') || (base1 == 'G' && base2 == 'A') ||
			(base1 == 'C' && base2 == 'T') || (base1 == 'T' && base2 == 'C') {
			transitions++
		}

		// Check for transversion patterns
		if (base1 == 'A' && (base2 == 'C' || base2 == 'T')) ||
			(base1 == 'G' && (base2 == 'C' || base2 == 'T')) ||
			(base1 == 'C' && (base2 == 'A' || base2 == 'G')) ||
			(base1 == 'T' && (base2 == 'A' || base2 == 'G')) {
			transversions++
		}
	}

	totalSubstitutions := transitions + transversions
	if totalSubstitutions == 0 || transversions == 0 {
		return 0.5 // Neutral score; ratio undefined when no transversions present
	}

	// Calculate transition/transversion ratio
	// Expected ratio in mammalian evolution is ~2:1
	tsRatio := float64(transitions) / float64(transversions)
	expectedRatio := 2.0

	// Score based on how close the ratio is to expected
	ratioScore := math.Exp(-math.Abs(tsRatio-expectedRatio) / expectedRatio)

	return ratioScore
}

// calculateTreeLikelihood scores the k-mer distances under a Gaussian centred
// on each species' divergence time, rescaled to [0, 1]. The two quantities
// have different units, so this is a heuristic, not a phylogenetic likelihood.
func (pdm *PhylogeneticDistanceModel) calculateTreeLikelihood(distances map[string]float64) float64 {
	if len(distances) == 0 {
		return 0.0
	}

	logLikelihood := 0.0
	count := 0

	for species, observedDist := range distances {
		if expectedDist, exists := pdm.speciesDistances["Human"][species]; exists {
			// Assume normal distribution around expected distance
			variance := expectedDist * 0.1 // arbitrary: variance = 0.1 × divergence time
			if variance == 0 {
				variance = 0.01
			}

			// Log-likelihood under normal distribution
			diff := observedDist - expectedDist
			logLikelihood += -0.5 * (diff*diff/variance + math.Log(2*math.Pi*variance))
			count++
		}
	}

	if count == 0 {
		return 0.0
	}

	// Convert to 0-1 scale (higher = better fit)
	avgLogLikelihood := logLikelihood / float64(count)
	likelihood := math.Exp(avgLogLikelihood / 10.0) // Scale factor

	if likelihood > 1.0 {
		likelihood = 1.0
	}

	return likelihood
}

// calculatePhylogeneticRealism calculates overall phylogenetic realism score
func (pdm *PhylogeneticDistanceModel) calculatePhylogeneticRealism(details *PhylogeneticAnalysisDetails) float64 {
	// Weighted combination of different phylogenetic measures
	weights := map[string]float64{
		"phylogenetic_signal":      0.3,
		"conservation_score":       0.2,
		"substitution_consistency": 0.2,
		"tree_likelihood":          0.3,
	}

	score := 0.0
	score += details.PhylogeneticSignal * weights["phylogenetic_signal"]
	score += details.ConservationScore * weights["conservation_score"]
	score += details.SubstitutionConsistency * weights["substitution_consistency"]
	score += details.TreeLikelihood * weights["tree_likelihood"]

	return score
}

// PhylogeneticAnalysisDetails contains detailed phylogenetic analysis results
type PhylogeneticAnalysisDetails struct {
	EvolutionaryDistances    map[string]float64
	PhylogeneticSignal       float64
	ConservationScore        float64
	LikelyEvolutionaryOrigin string
	SubstitutionConsistency  float64
	TreeLikelihood           float64
}

// calculateCorrelation calculates Pearson correlation coefficient
func calculateCorrelation(x, y []float64) float64 {
	if len(x) != len(y) || len(x) < 2 {
		return 0.0
	}

	n := float64(len(x))

	// Calculate means
	meanX, meanY := 0.0, 0.0
	for i := 0; i < len(x); i++ {
		meanX += x[i]
		meanY += y[i]
	}
	meanX /= n
	meanY /= n

	// Calculate correlation
	numerator := 0.0
	sumXSq := 0.0
	sumYSq := 0.0

	for i := 0; i < len(x); i++ {
		dx := x[i] - meanX
		dy := y[i] - meanY
		numerator += dx * dy
		sumXSq += dx * dx
		sumYSq += dy * dy
	}

	denominator := math.Sqrt(sumXSq * sumYSq)
	if denominator == 0 {
		return 0.0
	}

	return numerator / denominator
}
