package features

import (
	"fmt"
	"math"
)

// NewChromatinAccessibilityModel creates a new chromatin accessibility model
func NewChromatinAccessibilityModel() (*ChromatinAccessibilityModel, error) {
	model := &ChromatinAccessibilityModel{
		accessibilityScores: make(map[string]float64),
		histoneMarks:        make(map[string]map[string]float64),
		tfBindingSites:      make(map[string][]TranscriptionFactorSite),
		chromatinStates:     make(map[string]ChromatinState),
	}
	
	// Initialize with default chromatin landscape data
	if err := model.initializeChromatinLandscape(); err != nil {
		return nil, fmt.Errorf("failed to initialize chromatin landscape: %w", err)
	}
	
	return model, nil
}

// initializeChromatinLandscape sets up default chromatin accessibility patterns
func (cam *ChromatinAccessibilityModel) initializeChromatinLandscape() error {
	// Initialize histone modification patterns for different chromatin states
	cam.initializeHistoneMarks()
	
	// Initialize transcription factor binding sites
	cam.initializeTFBindingSites()
	
	// Initialize chromatin state annotations
	cam.initializeChromatinStates()
	
	return nil
}

// initializeHistoneMarks sets up histone modification patterns
func (cam *ChromatinAccessibilityModel) initializeHistoneMarks() {
	// Different histone marks associated with various chromatin states
	histoneMarks := map[string]map[string]float64{
		// Active promoter signatures
		"active_promoter": {
			"H3K4me3":  0.85, // Strong promoter mark
			"H3K27ac":  0.75, // Active enhancer/promoter
			"H3K9ac":   0.70, // Transcriptional activation
			"H3K4me1":  0.40, // Enhancer mark (lower at promoters)
			"H3K27me3": 0.05, // Repressive mark (low)
			"H3K9me3":  0.02, // Heterochromatin mark (very low)
		},
		
		// Enhancer signatures
		"strong_enhancer": {
			"H3K4me1":  0.80, // Strong enhancer mark
			"H3K27ac":  0.85, // Active enhancer mark
			"H3K4me3":  0.15, // Promoter mark (low)
			"H3K9ac":   0.60, // Transcriptional activation
			"H3K27me3": 0.08, // Repressive mark (low)
			"H3K9me3":  0.03, // Heterochromatin mark (very low)
		},
		
		// Transcribed regions
		"transcribed": {
			"H3K36me3": 0.75, // Transcriptional elongation
			"H3K4me1":  0.30, // Some enhancer activity
			"H3K27ac":  0.40, // Some active marks
			"H3K4me3":  0.20, // Some promoter activity
			"H3K27me3": 0.10, // Low repression
			"H3K9me3":  0.05, // Low heterochromatin
		},
		
		// Polycomb repressed regions
		"polycomb_repressed": {
			"H3K27me3": 0.80, // Strong Polycomb mark
			"H3K4me3":  0.05, // Very low active marks
			"H3K27ac":  0.03,
			"H3K4me1":  0.08,
			"H3K9ac":   0.02,
			"H3K9me3":  0.15, // Some heterochromatin
		},
		
		// Heterochromatin
		"heterochromatin": {
			"H3K9me3":  0.85, // Strong heterochromatin mark
			"H3K27me3": 0.30, // Some Polycomb
			"H3K4me3":  0.01, // Very low active marks
			"H3K27ac":  0.01,
			"H3K4me1":  0.02,
			"H3K9ac":   0.01,
		},
		
		// Inactive regions
		"inactive": {
			"H3K4me3":  0.05, // Low active marks
			"H3K27ac":  0.03,
			"H3K4me1":  0.08,
			"H3K9ac":   0.04,
			"H3K27me3": 0.20, // Moderate repression
			"H3K9me3":  0.25, // Moderate heterochromatin
		},
	}
	
	cam.histoneMarks = histoneMarks
}

// initializeTFBindingSites sets up transcription factor binding patterns
func (cam *ChromatinAccessibilityModel) initializeTFBindingSites() {
	// Common transcription factors and their binding patterns
	commonTFs := []struct {
		name  string
		motif string
		score float64
	}{
		{"CTCF", "CCGCGNGGNGGCAG", 0.85},      // Insulator
		{"SP1", "GGGCGG", 0.75},               // GC-rich promoters
		{"AP1", "TGAGTCA", 0.70},              // Stress response
		{"NFY", "CCAAT", 0.65},                // CCAAT box
		{"ELF1", "AGGAAGT", 0.60},             // ETS family
		{"GABPA", "CGGAAGT", 0.58},            // ETS family
		{"YY1", "GCGCCATNTT", 0.72},           // Polycomb recruitment
		{"REST", "TTCAGCACCATGGACAG", 0.80},   // Transcriptional repressor
	}
	
	// Simulate TF binding sites across different chromosomes
	chromosomes := []string{"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"}
	
	for _, chr := range chromosomes {
		sites := make([]TranscriptionFactorSite, 0)
		
		// Generate binding sites at regular intervals (simplified model)
		for pos := int64(1000); pos < 100000000; pos += 5000 {
			for _, tf := range commonTFs {
				// Probabilistic binding based on TF score
				if (pos*17 + int64(len(tf.name)))%100 < int64(tf.score*100) {
					site := TranscriptionFactorSite{
						Name:       tf.name,
						Position:   pos,
						Chromosome: chr,
						Score:      tf.score + 0.1*(math.Sin(float64(pos)/1000.0)),
						Motif:      tf.motif,
					}
					sites = append(sites, site)
				}
			}
		}
		
		cam.tfBindingSites[chr] = sites
	}
}

// initializeChromatinStates sets up chromatin state annotations
func (cam *ChromatinAccessibilityModel) initializeChromatinStates() {
	// Simulate chromatin states across genomic regions
	chromosomes := []string{"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"}
	
	for _, chr := range chromosomes {
		// Assign chromatin states based on genomic position patterns
		for pos := int64(0); pos < 100000000; pos += 10000 {
			regionKey := fmt.Sprintf("%s:%d-%d", chr, pos, pos+10000)
			
			// Simple model based on position
			positionFactor := float64(pos % 1000000) / 1000000.0
			
			var state ChromatinState
			if positionFactor < 0.05 {
				state = ActivePromoter
			} else if positionFactor < 0.15 {
				state = StrongEnhancer
			} else if positionFactor < 0.25 {
				state = WeakEnhancer
			} else if positionFactor < 0.40 {
				state = Txn_Elongation
			} else if positionFactor < 0.50 {
				state = Weak_Txn
			} else if positionFactor < 0.65 {
				state = Polycomb_Repressed
			} else if positionFactor < 0.80 {
				state = Heterochromatin_low
			} else if positionFactor < 0.90 {
				state = Inactive_Low
			} else {
				state = Repetitive_CNV
			}
			
			cam.chromatinStates[regionKey] = state
		}
	}
}

// AnalyzeRegion analyzes chromatin accessibility for a genomic region
func (cam *ChromatinAccessibilityModel) AnalyzeRegion(chromosome string, start, end int64) (float64, *ChromatinAnalysisDetails, error) {
	if start >= end {
		return 0.0, nil, fmt.Errorf("invalid genomic coordinates: %s:%d-%d", chromosome, start, end)
	}
	
	// Get chromatin state for the region
	chromatinState := cam.getChromatinState(chromosome, start, end)
	
	// Calculate accessibility score based on chromatin state
	accessibilityScore := cam.calculateAccessibilityScore(chromatinState)
	
	// Get histone modification profile
	histoneProfile := cam.getHistoneProfile(chromatinState)
	
	// Find overlapping transcription factor binding sites
	overlappingTFs := cam.findOverlappingTFs(chromosome, start, end)
	
	// Calculate regulatory potential
	regulatoryPotential := cam.calculateRegulatoryPotential(chromatinState, overlappingTFs)
	
	// Calculate DNase I hypersensitivity score (proxy for accessibility)
	dnaseScore := cam.calculateDNaseScore(chromatinState, accessibilityScore)
	
	details := &ChromatinAnalysisDetails{
		ChromatinState:       chromatinState,
		AccessibilityScore:   accessibilityScore,
		HistoneProfile:       histoneProfile,
		OverlappingTFs:       overlappingTFs,
		RegulatoryPotential:  regulatoryPotential,
		DNaseScore:          dnaseScore,
		RegionSize:          end - start,
	}
	
	// Overall chromatin realism score
	realismScore := cam.calculateChromatinRealism(details)
	
	return realismScore, details, nil
}

// getChromatinState determines the predominant chromatin state for a region
func (cam *ChromatinAccessibilityModel) getChromatinState(chromosome string, start, end int64) ChromatinState {
	// Find the region that overlaps with our query
	regionStart := (start / 10000) * 10000
	regionKey := fmt.Sprintf("%s:%d-%d", chromosome, regionStart, regionStart+10000)
	
	if state, exists := cam.chromatinStates[regionKey]; exists {
		return state
	}
	
	// Default to inactive if not found
	return Inactive_Low
}

// calculateAccessibilityScore calculates accessibility based on chromatin state
func (cam *ChromatinAccessibilityModel) calculateAccessibilityScore(state ChromatinState) float64 {
	// Accessibility scores for different chromatin states
	accessibilityMap := map[ChromatinState]float64{
		ActivePromoter:      0.90,
		WeakPromoter:        0.70,
		PoisedPromoter:      0.60,
		StrongEnhancer:      0.85,
		WeakEnhancer:        0.65,
		Txn_Transition:      0.50,
		Txn_Elongation:      0.55,
		Weak_Txn:           0.40,
		Polycomb_Repressed:  0.15,
		Heterochromatin_low: 0.10,
		Repetitive_CNV:      0.08,
		Inactive_Low:        0.20,
	}
	
	if score, exists := accessibilityMap[state]; exists {
		return score
	}
	
	return 0.20 // Default low accessibility
}

// getHistoneProfile gets the histone modification profile for a chromatin state
func (cam *ChromatinAccessibilityModel) getHistoneProfile(state ChromatinState) map[string]float64 {
	stateNames := map[ChromatinState]string{
		ActivePromoter:      "active_promoter",
		StrongEnhancer:      "strong_enhancer",
		WeakEnhancer:        "strong_enhancer", // Use same as strong enhancer
		Txn_Elongation:      "transcribed",
		Weak_Txn:           "transcribed",
		Polycomb_Repressed:  "polycomb_repressed",
		Heterochromatin_low: "heterochromatin",
		Inactive_Low:        "inactive",
	}
	
	if stateName, exists := stateNames[state]; exists {
		if profile, exists := cam.histoneMarks[stateName]; exists {
			return profile
		}
	}
	
	// Return default inactive profile
	return cam.histoneMarks["inactive"]
}

// findOverlappingTFs finds transcription factors binding sites overlapping with the region
func (cam *ChromatinAccessibilityModel) findOverlappingTFs(chromosome string, start, end int64) []TranscriptionFactorSite {
	var overlapping []TranscriptionFactorSite
	
	if sites, exists := cam.tfBindingSites[chromosome]; exists {
		for _, site := range sites {
			if site.Position >= start && site.Position <= end {
				overlapping = append(overlapping, site)
			}
		}
	}
	
	return overlapping
}

// calculateRegulatoryPotential calculates the regulatory potential of the region
func (cam *ChromatinAccessibilityModel) calculateRegulatoryPotential(state ChromatinState, tfSites []TranscriptionFactorSite) float64 {
	// Base regulatory potential from chromatin state
	baseRegulatory := map[ChromatinState]float64{
		ActivePromoter:     0.95,
		WeakPromoter:       0.75,
		PoisedPromoter:     0.60,
		StrongEnhancer:     0.90,
		WeakEnhancer:       0.70,
		Txn_Transition:     0.40,
		Txn_Elongation:     0.30,
		Weak_Txn:          0.25,
		Polycomb_Repressed: 0.05,
		Heterochromatin_low: 0.02,
		Repetitive_CNV:     0.01,
		Inactive_Low:       0.15,
	}
	
	basePotential := baseRegulatory[state]
	
	// Add contribution from transcription factor binding sites
	tfContribution := 0.0
	for _, site := range tfSites {
		tfContribution += site.Score * 0.1 // Each TF contributes up to 0.1
	}
	
	// Cap the total potential at 1.0
	totalPotential := basePotential + tfContribution
	if totalPotential > 1.0 {
		totalPotential = 1.0
	}
	
	return totalPotential
}

// calculateDNaseScore calculates DNase I hypersensitivity score
func (cam *ChromatinAccessibilityModel) calculateDNaseScore(state ChromatinState, accessibilityScore float64) float64 {
	// DNase score is highly correlated with accessibility but with some noise
	baseScore := accessibilityScore
	
	// Add chromatin state-specific modulation
	stateModulation := map[ChromatinState]float64{
		ActivePromoter:     1.1,  // Higher DNase sensitivity
		StrongEnhancer:     1.15, // Highest sensitivity
		WeakEnhancer:       1.05,
		Txn_Elongation:     0.9,  // Lower in gene bodies
		Polycomb_Repressed: 0.7,  // Reduced accessibility
		Heterochromatin_low: 0.5, // Very low accessibility
		Inactive_Low:       0.8,
	}
	
	modulation := stateModulation[state]
	if modulation == 0 {
		modulation = 1.0
	}
	
	dnaseScore := baseScore * modulation
	
	// Add some biological noise (±10%)
	noise := 0.1 * (math.Sin(float64(int(state)*17)) - 0.5)
	dnaseScore += noise
	
	// Ensure score stays in [0, 1] range
	if dnaseScore < 0 {
		dnaseScore = 0
	}
	if dnaseScore > 1 {
		dnaseScore = 1
	}
	
	return dnaseScore
}

// calculateChromatinRealism calculates overall chromatin realism score
func (cam *ChromatinAccessibilityModel) calculateChromatinRealism(details *ChromatinAnalysisDetails) float64 {
	// Weighted combination of chromatin features
	weights := map[string]float64{
		"accessibility":        0.3,
		"regulatory_potential": 0.25,
		"dnase_score":         0.2,
		"tf_density":          0.15,
		"histone_consistency": 0.1,
	}
	
	// Calculate TF density score
	tfDensity := float64(len(details.OverlappingTFs)) / float64(details.RegionSize/1000.0)
	if tfDensity > 1.0 {
		tfDensity = 1.0
	}
	
	// Calculate histone mark consistency
	histoneConsistency := cam.calculateHistoneConsistency(details.HistoneProfile, details.ChromatinState)
	
	// Combine scores
	realismScore := 0.0
	realismScore += details.AccessibilityScore * weights["accessibility"]
	realismScore += details.RegulatoryPotential * weights["regulatory_potential"]
	realismScore += details.DNaseScore * weights["dnase_score"]
	realismScore += tfDensity * weights["tf_density"]
	realismScore += histoneConsistency * weights["histone_consistency"]
	
	return realismScore
}

// calculateHistoneConsistency checks consistency of histone marks with chromatin state
func (cam *ChromatinAccessibilityModel) calculateHistoneConsistency(observedProfile map[string]float64, state ChromatinState) float64 {
	expectedProfile := cam.getHistoneProfile(state)
	
	consistency := 0.0
	markCount := 0
	
	for mark, expectedLevel := range expectedProfile {
		if observedLevel, exists := observedProfile[mark]; exists {
			// Calculate how close observed is to expected
			diff := math.Abs(observedLevel - expectedLevel)
			markConsistency := 1.0 - diff // 1.0 = perfect match, 0.0 = completely off
			consistency += markConsistency
			markCount++
		}
	}
	
	if markCount == 0 {
		return 0.5 // Neutral score if no marks to compare
	}
	
	return consistency / float64(markCount)
}

// ChromatinAnalysisDetails contains detailed chromatin analysis results
type ChromatinAnalysisDetails struct {
	ChromatinState      ChromatinState
	AccessibilityScore  float64
	HistoneProfile      map[string]float64
	OverlappingTFs      []TranscriptionFactorSite
	RegulatoryPotential float64
	DNaseScore         float64
	RegionSize         int64
}

// String returns a string representation of the chromatin state
func (cs ChromatinState) String() string {
	stateNames := map[ChromatinState]string{
		ActivePromoter:      "Active Promoter",
		WeakPromoter:        "Weak Promoter", 
		PoisedPromoter:      "Poised Promoter",
		StrongEnhancer:      "Strong Enhancer",
		WeakEnhancer:        "Weak Enhancer",
		Txn_Transition:      "Transcriptional Transition",
		Txn_Elongation:      "Transcriptional Elongation",
		Weak_Txn:           "Weak Transcription",
		Polycomb_Repressed:  "Polycomb Repressed",
		Heterochromatin_low: "Heterochromatin",
		Repetitive_CNV:      "Repetitive/CNV",
		Inactive_Low:        "Inactive",
	}
	
	if name, exists := stateNames[cs]; exists {
		return name
	}
	
	return "Unknown"
}