package features

import (
	"fmt"
	"math"
)

// NewSelectionPressureModel builds the synthetic selection annotation for
// chr1–chr10: conservation scores every 1 kb up to 100 Mb (10^6 entries),
// selection coefficients, and 1000 synthetic genes per gene class with
// formula-generated dN/dS ratios.
func NewSelectionPressureModel() (*SelectionPressureModel, error) {
	model := &SelectionPressureModel{
		purifyingSelection: make(map[string]float64),
		positiveSelection:  make(map[string]float64),
		balancingSelection: make(map[string]float64),
		dNdSRatios:         make(map[string]float64),
		conservationScores: make(map[string]float64),
	}

	// Initialize with genomic selection patterns
	if err := model.initializeSelectionPatterns(); err != nil {
		return nil, fmt.Errorf("failed to initialize selection patterns: %w", err)
	}

	return model, nil
}

// initializeSelectionPatterns sets up selection pressure patterns across the genome
func (spm *SelectionPressureModel) initializeSelectionPatterns() error {
	// Initialize conservation scores based on genomic elements
	spm.initializeConservationScores()

	// Initialize dN/dS ratios for protein-coding genes
	spm.initializeDnDsRatios()

	// Initialize selection coefficients for different genomic regions
	spm.initializeSelectionCoefficients()

	return nil
}

// initializeConservationScores fills synthetic conservation scores from the
// position-based context in determineGenomicContext; no PhyloP or PhastCons
// data is used.
func (spm *SelectionPressureModel) initializeConservationScores() {
	chromosomes := []string{"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"}

	for _, chr := range chromosomes {
		// Generate conservation scores based on genomic position and element type
		for pos := int64(0); pos < 100000000; pos += 1000 {
			regionKey := fmt.Sprintf("%s:%d", chr, pos)

			// Model different conservation levels based on genomic context
			conservationScore := spm.calculatePositionalConservation(chr, pos)
			spm.conservationScores[regionKey] = conservationScore
		}
	}
}

// calculatePositionalConservation calculates position-specific conservation
func (spm *SelectionPressureModel) calculatePositionalConservation(chromosome string, position int64) float64 {
	// Simple model based on genomic features

	// Base conservation (intergenic regions)
	baseConservation := 0.1

	// Determine genomic context
	context := spm.determineGenomicContext(chromosome, position)

	switch context {
	case "exon":
		// Protein-coding exons are highly conserved
		return 0.8 + 0.1*math.Sin(float64(position)/10000.0) // 0.7-0.9

	case "utr5":
		// 5' UTRs are moderately conserved
		return 0.6 + 0.1*math.Sin(float64(position)/8000.0) // 0.5-0.7

	case "utr3":
		// 3' UTRs are moderately conserved
		return 0.5 + 0.1*math.Sin(float64(position)/12000.0) // 0.4-0.6

	case "promoter":
		// Promoters are highly conserved
		return 0.85 + 0.05*math.Sin(float64(position)/5000.0) // 0.8-0.9

	case "enhancer":
		// Enhancers show variable conservation
		return 0.4 + 0.2*math.Sin(float64(position)/15000.0) // 0.2-0.6

	case "intron":
		// Introns show low to moderate conservation
		return 0.2 + 0.1*math.Sin(float64(position)/20000.0) // 0.1-0.3

	case "repetitive":
		// Repetitive elements are poorly conserved
		return 0.05 + 0.05*math.Sin(float64(position)/25000.0) // 0.0-0.1

	default:
		// Intergenic regions
		return baseConservation
	}
}

// determineGenomicContext assigns a context from position mod 100 kb (the
// same promoter/exon/UTR/intron/enhancer/repeat/intergenic layout repeats
// every 100 kb). It is not a gene annotation; chromosome is ignored.
func (spm *SelectionPressureModel) determineGenomicContext(chromosome string, position int64) string {
	positionMod := position % 100000

	if positionMod < 5000 {
		return "promoter"
	} else if positionMod < 15000 {
		return "exon"
	} else if positionMod < 18000 {
		return "utr5"
	} else if positionMod < 22000 {
		return "utr3"
	} else if positionMod < 45000 {
		return "intron"
	} else if positionMod < 50000 {
		return "enhancer"
	} else if positionMod < 60000 {
		return "repetitive"
	} else {
		return "intergenic"
	}
}

// initializeDnDsRatios creates synthetic genes per class with dN/dS ratios
// from a deterministic sine formula around a hand-set class mean; no
// alignments are analysed.
func (spm *SelectionPressureModel) initializeDnDsRatios() {
	// Gene classes with hand-set mean dN/dS ratios
	geneClasses := map[string]struct {
		meanDnDs    float64
		variance    float64
		description string
	}{
		"housekeeping": {
			meanDnDs:    0.15,
			variance:    0.05,
			description: "Essential housekeeping genes under strong purifying selection",
		},
		"metabolic": {
			meanDnDs:    0.25,
			variance:    0.08,
			description: "Metabolic genes with moderate constraint",
		},
		"immune": {
			meanDnDs:    0.8,
			variance:    0.3,
			description: "Immune genes under balancing/positive selection",
		},
		"sensory": {
			meanDnDs:    0.6,
			variance:    0.2,
			description: "Sensory genes with relaxed constraint",
		},
		"reproductive": {
			meanDnDs:    1.2,
			variance:    0.4,
			description: "Reproductive genes under positive selection",
		},
		"transcription": {
			meanDnDs:    0.2,
			variance:    0.06,
			description: "Transcription factors under strong constraint",
		},
		"development": {
			meanDnDs:    0.18,
			variance:    0.05,
			description: "Developmental genes under strong constraint",
		},
	}

	// Generate dN/dS ratios for genes in each class
	for geneClass, params := range geneClasses {
		numGenesInClass := 1000 // 1000 synthetic genes per class

		for i := 0; i < numGenesInClass; i++ {
			geneID := fmt.Sprintf("%s_gene_%d", geneClass, i)

			// Deterministic variation around the class mean (not random)
			noise := params.variance * (math.Sin(float64(i)*3.14159/100.0) - 0.5)
			dnds := params.meanDnDs + noise

			// Ensure dN/dS stays positive
			if dnds < 0.01 {
				dnds = 0.01
			}

			spm.dNdSRatios[geneID] = dnds
		}
	}
}

// initializeSelectionCoefficients sets up selection coefficients
func (spm *SelectionPressureModel) initializeSelectionCoefficients() {
	chromosomes := []string{"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"}

	for _, chr := range chromosomes {
		for pos := int64(0); pos < 100000000; pos += 10000 {
			regionKey := fmt.Sprintf("%s:%d-%d", chr, pos, pos+10000)

			// Calculate selection coefficients based on genomic context
			context := spm.determineGenomicContext(chr, pos)

			purifying := spm.calculatePurifyingSelection(context, pos)
			positive := spm.calculatePositiveSelection(context, pos)
			balancing := spm.calculateBalancingSelection(context, pos)

			spm.purifyingSelection[regionKey] = purifying
			spm.positiveSelection[regionKey] = positive
			spm.balancingSelection[regionKey] = balancing
		}
	}
}

// calculatePurifyingSelection calculates purifying selection strength
func (spm *SelectionPressureModel) calculatePurifyingSelection(context string, position int64) float64 {
	// Selection coefficients for different genomic contexts
	basePurifying := map[string]float64{
		"exon":       0.8,  // Strong purifying selection
		"promoter":   0.75, // Strong selection on regulatory elements
		"utr5":       0.6,  // Moderate selection
		"utr3":       0.5,  // Moderate selection
		"enhancer":   0.4,  // Variable selection
		"intron":     0.2,  // Weak selection
		"intergenic": 0.1,  // Very weak selection
		"repetitive": 0.05, // Almost no selection
	}

	base := basePurifying[context]
	if base == 0 {
		base = 0.1 // Default
	}

	// Add position-dependent variation
	variation := 0.1 * math.Sin(float64(position)/7500.0)
	selection := base + variation

	// Ensure positive values
	if selection < 0 {
		selection = 0.01
	}
	if selection > 1.0 {
		selection = 1.0
	}

	return selection
}

// calculatePositiveSelection calculates positive selection strength
func (spm *SelectionPressureModel) calculatePositiveSelection(context string, position int64) float64 {
	// Positive selection is rare and context-specific
	basePositive := map[string]float64{
		"exon":       0.05, // Rare in most exons
		"promoter":   0.02, // Very rare in promoters
		"utr5":       0.01,
		"utr3":       0.03, // Slightly higher in 3' UTRs
		"enhancer":   0.08, // Higher in regulatory elements
		"intron":     0.02,
		"intergenic": 0.01,
		"repetitive": 0.15, // Higher in repetitive elements
	}

	base := basePositive[context]

	// Synthetic "hotspot": the first 50 kb of every megabase (not real
	// immune gene clusters)
	if (position % 1000000) < 50000 {
		base *= 5.0 // Higher positive selection
	}

	// Deterministic sinusoidal variation (not stochastic)
	variation := 0.02 * math.Sin(float64(position)/12000.0)
	selection := base + variation

	if selection < 0 {
		selection = 0.0
	}
	if selection > 1.0 {
		selection = 1.0
	}

	return selection
}

// calculateBalancingSelection calculates balancing selection strength
func (spm *SelectionPressureModel) calculateBalancingSelection(context string, position int64) float64 {
	// Balancing selection is very rare but important in specific regions
	baseBalancing := map[string]float64{
		"exon":       0.01, // Rare but present (e.g., MHC genes)
		"promoter":   0.005,
		"utr5":       0.002,
		"utr3":       0.005,
		"enhancer":   0.01,
		"intron":     0.002,
		"intergenic": 0.001,
		"repetitive": 0.001,
	}

	base := baseBalancing[context]

	// Synthetic "MHC-like" region: the first 100 kb of every 5 Mb (not the
	// real HLA locus)
	if (position % 5000000) < 100000 {
		base *= 20.0 // Much higher balancing selection
	}

	// Deterministic sinusoidal variation
	variation := 0.005 * math.Sin(float64(position)/15000.0)
	selection := base + variation

	if selection < 0 {
		selection = 0.0
	}
	if selection > 1.0 {
		selection = 1.0
	}

	return selection
}

// AnalyzeRegion looks up the synthetic selection coefficients, conservation
// score and dN/dS ratio for a region and returns a weighted heuristic score.
// It uses only coordinates, never sequence, so it says nothing about whether
// a sequence shows signs of selection.
func (spm *SelectionPressureModel) AnalyzeRegion(chromosome string, start, end int64) (float64, *SelectionAnalysisDetails, error) {
	if start >= end {
		return 0.0, nil, fmt.Errorf("invalid genomic coordinates: %s:%d-%d", chromosome, start, end)
	}

	// Get selection coefficients for the region
	regionKey := fmt.Sprintf("%s:%d-%d", chromosome, (start/10000)*10000, ((end/10000)+1)*10000)

	purifying := spm.getSelectionCoefficient(spm.purifyingSelection, regionKey, start)
	positive := spm.getSelectionCoefficient(spm.positiveSelection, regionKey, start)
	balancing := spm.getSelectionCoefficient(spm.balancingSelection, regionKey, start)

	// Get conservation score
	conservationKey := fmt.Sprintf("%s:%d", chromosome, (start/1000)*1000)
	conservation := spm.getConservationScore(conservationKey, start)

	// Calculate effective selection strength
	effectiveSelection := spm.calculateEffectiveSelection(purifying, positive, balancing)

	// Determine dominant selection type
	dominantSelection := spm.determineDominantSelection(purifying, positive, balancing)

	// Calculate selection consistency across the region
	consistency := spm.calculateSelectionConsistency(chromosome, start, end)

	// Get dN/dS information if in coding region
	dndsInfo := spm.getDnDsInformation(chromosome, start, end)

	details := &SelectionAnalysisDetails{
		PurifyingSelection:   purifying,
		PositiveSelection:    positive,
		BalancingSelection:   balancing,
		ConservationScore:    conservation,
		EffectiveSelection:   effectiveSelection,
		DominantSelection:    dominantSelection,
		SelectionConsistency: consistency,
		DnDsInformation:      dndsInfo,
		RegionSize:           end - start,
	}

	// Calculate overall selection realism score
	realismScore := spm.calculateSelectionRealism(details)

	return realismScore, details, nil
}

// getSelectionCoefficient retrieves selection coefficient with interpolation
func (spm *SelectionPressureModel) getSelectionCoefficient(selectionMap map[string]float64, regionKey string, position int64) float64 {
	if coeff, exists := selectionMap[regionKey]; exists {
		return coeff
	}

	// Default values based on position
	context := spm.determineGenomicContext("chr1", position) // Use chr1 as default

	switch context {
	case "exon":
		return 0.7
	case "promoter":
		return 0.6
	case "enhancer":
		return 0.3
	default:
		return 0.1
	}
}

// getConservationScore retrieves conservation score with interpolation
func (spm *SelectionPressureModel) getConservationScore(conservationKey string, position int64) float64 {
	if score, exists := spm.conservationScores[conservationKey]; exists {
		return score
	}

	// Default conservation based on genomic context
	return spm.calculatePositionalConservation("chr1", position)
}

// calculateEffectiveSelection calculates effective selection strength
func (spm *SelectionPressureModel) calculateEffectiveSelection(purifying, positive, balancing float64) float64 {
	// Effective selection considers all types of selection
	// Positive and balancing selection can counteract purifying selection

	effectiveS := purifying - (positive+balancing)/2.0

	// Ensure reasonable bounds
	if effectiveS < 0 {
		effectiveS = 0
	}
	if effectiveS > 1.0 {
		effectiveS = 1.0
	}

	return effectiveS
}

// determineDominantSelection determines the dominant type of selection
func (spm *SelectionPressureModel) determineDominantSelection(purifying, positive, balancing float64) SelectionType {
	if purifying > positive && purifying > balancing {
		if purifying > 0.5 {
			return StrongPurifying
		} else {
			return WeakPurifying
		}
	} else if positive > balancing {
		if positive > 0.1 {
			return PositiveSelection
		} else {
			return NeutralEvolution
		}
	} else if balancing > 0.05 {
		return BalancingSelection
	} else {
		return NeutralEvolution
	}
}

// calculateSelectionConsistency calculates consistency of selection across region
func (spm *SelectionPressureModel) calculateSelectionConsistency(chromosome string, start, end int64) float64 {
	// Sample selection coefficients across the region
	numSamples := int((end - start) / 5000) // Sample every 5kb
	if numSamples < 2 {
		return 1.0 // Perfect consistency for small regions
	}
	if numSamples > 20 {
		numSamples = 20 // Limit sampling
	}

	var purifyingValues []float64
	step := (end - start) / int64(numSamples)

	for i := 0; i < numSamples; i++ {
		pos := start + int64(i)*step
		regionKey := fmt.Sprintf("%s:%d-%d", chromosome, (pos/10000)*10000, ((pos/10000)+1)*10000)
		purifying := spm.getSelectionCoefficient(spm.purifyingSelection, regionKey, pos)
		purifyingValues = append(purifyingValues, purifying)
	}

	// Calculate coefficient of variation
	if len(purifyingValues) < 2 {
		return 1.0
	}

	sum := 0.0
	for _, val := range purifyingValues {
		sum += val
	}
	mean := sum / float64(len(purifyingValues))

	if mean == 0 {
		return 1.0
	}

	variance := 0.0
	for _, val := range purifyingValues {
		diff := val - mean
		variance += diff * diff
	}
	variance /= float64(len(purifyingValues) - 1) // Bessel's correction; len >= 2 guaranteed above
	stddev := math.Sqrt(variance)

	cv := stddev / mean

	// Convert to consistency score (higher = more consistent)
	consistency := math.Exp(-cv)

	return consistency
}

// getDnDsInformation returns a synthetic gene's dN/dS ratio when the region
// midpoint falls in a synthetic exon. The gene class is picked by iterating
// a Go map, so which class (and ratio) is returned varies between calls.
func (spm *SelectionPressureModel) getDnDsInformation(chromosome string, start, end int64) *DnDsInfo {
	midpoint := (start + end) / 2
	context := spm.determineGenomicContext(chromosome, midpoint)

	if context != "exon" {
		return nil // No dN/dS for non-coding regions
	}

	// Check if we have dN/dS data for this gene
	var dnds float64
	var geneClass string

	// Try to find matching gene class
	for class, _ := range map[string]bool{"housekeeping": true, "metabolic": true, "immune": true, "sensory": true, "reproductive": true, "transcription": true, "development": true} {
		testGeneID := fmt.Sprintf("%s_gene_%d", class, int(midpoint/10000)%1000)
		if ratio, exists := spm.dNdSRatios[testGeneID]; exists {
			dnds = ratio
			geneClass = class
			break
		}
	}

	if geneClass == "" {
		// Default to housekeeping
		dnds = 0.2
		geneClass = "housekeeping"
	}

	return &DnDsInfo{
		DnDsRatio:     dnds,
		GeneClass:     geneClass,
		SelectionType: spm.interpretDnDsRatio(dnds),
	}
}

// interpretDnDsRatio interprets dN/dS ratio in terms of selection
func (spm *SelectionPressureModel) interpretDnDsRatio(dnds float64) SelectionType {
	if dnds < 0.3 {
		return StrongPurifying
	} else if dnds < 0.7 {
		return WeakPurifying
	} else if dnds < 1.2 {
		return NeutralEvolution
	} else {
		return PositiveSelection
	}
}

// calculateSelectionRealism calculates overall selection realism score
func (spm *SelectionPressureModel) calculateSelectionRealism(details *SelectionAnalysisDetails) float64 {
	// Weighted combination of selection features
	weights := map[string]float64{
		"conservation_score":    0.25,
		"effective_selection":   0.2,
		"selection_consistency": 0.2,
		"selection_balance":     0.15,
		"dnds_consistency":      0.1,
		"region_size_factor":    0.1,
	}

	// Calculate selection balance (realistic mix of selection types)
	selectionBalance := spm.calculateSelectionBalance(details)

	// Calculate dN/dS consistency
	dndsConsistency := 0.5 // Default neutral
	if details.DnDsInformation != nil {
		dndsConsistency = spm.calculateDnDsConsistency(details)
	}

	// Region size factor (larger regions should show more consistency)
	regionSizeFactor := math.Min(float64(details.RegionSize)/10000.0, 1.0)

	// Combine scores
	realismScore := 0.0
	realismScore += details.ConservationScore * weights["conservation_score"]
	realismScore += details.EffectiveSelection * weights["effective_selection"]
	realismScore += details.SelectionConsistency * weights["selection_consistency"]
	realismScore += selectionBalance * weights["selection_balance"]
	realismScore += dndsConsistency * weights["dnds_consistency"]
	realismScore += regionSizeFactor * weights["region_size_factor"]

	return realismScore
}

// calculateSelectionBalance calculates how realistic the mix of selection types is
func (spm *SelectionPressureModel) calculateSelectionBalance(details *SelectionAnalysisDetails) float64 {
	// Realistic selection should be dominated by purifying selection
	// with small amounts of positive and balancing selection

	total := details.PurifyingSelection + details.PositiveSelection + details.BalancingSelection
	if total == 0 {
		return 0.5 // Neutral
	}

	purifyingFraction := details.PurifyingSelection / total
	positiveFraction := details.PositiveSelection / total
	balancingFraction := details.BalancingSelection / total

	// Target fractions, hand-set (not measured genomic averages)
	idealPurifying := 0.7
	idealPositive := 0.05
	idealBalancing := 0.02

	// Calculate how close to ideal
	purifyingScore := 1.0 - math.Abs(purifyingFraction-idealPurifying)
	positiveScore := 1.0 - math.Abs(positiveFraction-idealPositive)*10.0    // Higher penalty for deviation
	balancingScore := 1.0 - math.Abs(balancingFraction-idealBalancing)*20.0 // Higher penalty for deviation

	// Ensure scores stay positive
	if purifyingScore < 0 {
		purifyingScore = 0
	}
	if positiveScore < 0 {
		positiveScore = 0
	}
	if balancingScore < 0 {
		balancingScore = 0
	}

	// Weighted average
	balance := 0.6*purifyingScore + 0.25*positiveScore + 0.15*balancingScore

	return balance
}

// expectedDnDsRatio maps each selection type to a representative dN/dS
// ratio. Only the interpretation (ω < 1 purifying, ω = 1 neutral, ω > 1
// positive) comes from Kimura (1983) and Yang & Bielawski (2000); the
// specific values, and 1.5 for balancing selection, are illustrative choices.
var expectedDnDsRatio = map[SelectionType]float64{
	StrongPurifying:    0.1,
	WeakPurifying:      0.3,
	NeutralEvolution:   1.0,
	PositiveSelection:  2.0,
	BalancingSelection: 1.5,
}

// calculateDnDsConsistency returns a continuous score in [0, 1] measuring how
// well the observed dN/dS ratio agrees with the canonical ratio expected for
// the dominant selection type.
//
// Score = exp(−|log(observed / expected)|)
//
//	1.0  → observed ratio matches expectation exactly
//	→ 0  → increasing divergence on a log scale
//
// Returns 0.0 (cannot compute) when DnDsInformation is nil, when the observed
// ratio is non-positive, or when the dominant selection type has no known
// canonical ratio.
func (spm *SelectionPressureModel) calculateDnDsConsistency(details *SelectionAnalysisDetails) float64 {
	if details.DnDsInformation == nil {
		return 0.0
	}

	observed := details.DnDsInformation.DnDsRatio
	if observed <= 0 {
		return 0.0
	}

	expected, ok := expectedDnDsRatio[details.DominantSelection]
	if !ok || expected <= 0 {
		return 0.0
	}

	logRatio := math.Abs(math.Log(observed / expected))
	return math.Exp(-logRatio)
}

// SelectionAnalysisDetails contains detailed selection pressure analysis results
type SelectionAnalysisDetails struct {
	PurifyingSelection   float64
	PositiveSelection    float64
	BalancingSelection   float64
	ConservationScore    float64
	EffectiveSelection   float64
	DominantSelection    SelectionType
	SelectionConsistency float64
	DnDsInformation      *DnDsInfo
	RegionSize           int64
}

// SelectionType represents different types of natural selection
type SelectionType int

const (
	StrongPurifying SelectionType = iota
	WeakPurifying
	NeutralEvolution
	PositiveSelection
	BalancingSelection
)

func (st SelectionType) String() string {
	selectionNames := map[SelectionType]string{
		StrongPurifying:    "Strong Purifying Selection",
		WeakPurifying:      "Weak Purifying Selection",
		NeutralEvolution:   "Neutral Evolution",
		PositiveSelection:  "Positive Selection",
		BalancingSelection: "Balancing Selection",
	}

	if name, exists := selectionNames[st]; exists {
		return name
	}

	return "Unknown"
}

// DnDsInfo contains dN/dS analysis information
type DnDsInfo struct {
	DnDsRatio     float64
	GeneClass     string
	SelectionType SelectionType
}
