package features

import (
	"fmt"
	"math"
)

// NewMutationalSignatureModel creates a new mutational signature model
func NewMutationalSignatureModel(processes []MutationalProcess) (*MutationalSignatureModel, error) {
	model := &MutationalSignatureModel{
		signatures:        make(map[string]*MutationalSignature),
		trinucleotideFreq: make(map[string]float64),
		processActivity:   make(map[string]float64),
	}

	// Initialize known mutational signatures
	if err := model.initializeMutationalSignatures(); err != nil {
		return nil, fmt.Errorf("failed to initialize mutational signatures: %w", err)
	}

	// Set activity levels for specified processes
	model.setProcessActivity(processes)

	return model, nil
}

// initializeMutationalSignatures sets up five profiles named after COSMIC v2
// signatures 1, 2, 4, 7 and 15. The trinucleotide weights are hand-written
// approximations, not the published COSMIC values (for example, signature 1
// here is not dominated by C>T at NpCpG).
func (msm *MutationalSignatureModel) initializeMutationalSignatures() error {

	// Signature 1: Age-related (clock-like)
	sig1 := &MutationalSignature{
		Name:        "Signature_1",
		Process:     ClockLikeMutation,
		Etiology:    "Age-related spontaneous deamination of 5-methylcytosine",
		CancerTypes: []string{"All"},
		TrinucleotideProfile: map[string]float64{
			"ACA>ATA": 0.012, "ACA>AGA": 0.008, "ACA>AAA": 0.004,
			"ACC>ATC": 0.015, "ACC>AGC": 0.010, "ACC>AAC": 0.005,
			"ACG>ATG": 0.008, "ACG>AGG": 0.012, "ACG>AAG": 0.003,
			"ACT>ATT": 0.018, "ACT>AGT": 0.014, "ACT>AAT": 0.007,
			"CCA>CTA": 0.022, "CCA>CGA": 0.008, "CCA>CAA": 0.004,
			"CCC>CTC": 0.025, "CCC>CGC": 0.010, "CCC>CAC": 0.005,
			"CCG>CTG": 0.035, "CCG>CGG": 0.015, "CCG>CAG": 0.008,
			"CCT>CTT": 0.045, "CCT>CGT": 0.018, "CCT>CAT": 0.010,
			"GCA>GTA": 0.008, "GCA>GAA": 0.005, "GCA>GGA": 0.003,
			"GCC>GTC": 0.010, "GCC>GAC": 0.007, "GCC>GGC": 0.004,
			"GCG>GTG": 0.012, "GCG>GAG": 0.008, "GCG>GGG": 0.005,
			"GCT>GTT": 0.015, "GCT>GAT": 0.010, "GCT>GGT": 0.006,
			"TCA>TTA": 0.032, "TCA>TGA": 0.012, "TCA>TAA": 0.008,
			"TCC>TTC": 0.028, "TCC>TGC": 0.010, "TCC>TAC": 0.007,
			"TCG>TTG": 0.025, "TCG>TGG": 0.008, "TCG>TAG": 0.005,
			"TCT>TTT": 0.038, "TCT>TGT": 0.014, "TCT>TAT": 0.009,
		},
	}
	msm.signatures["Signature_1"] = sig1

	// Signature 2: APOBEC cytidine deaminase activity
	sig2 := &MutationalSignature{
		Name:        "Signature_2",
		Process:     APOBECActivity,
		Etiology:    "Activity of APOBEC family of cytidine deaminases",
		CancerTypes: []string{"Bladder", "Cervical", "Breast", "Head_neck", "Lung"},
		TrinucleotideProfile: map[string]float64{
			"ACA>ATA": 0.005, "ACA>AGA": 0.003, "ACA>AAA": 0.002,
			"ACC>ATC": 0.008, "ACC>AGC": 0.005, "ACC>AAC": 0.003,
			"ACG>ATG": 0.012, "ACG>AGG": 0.008, "ACG>AAG": 0.004,
			"ACT>ATT": 0.015, "ACT>AGT": 0.010, "ACT>AAT": 0.005,
			"CCA>CTA": 0.045, "CCA>CGA": 0.015, "CCA>CAA": 0.008,
			"CCC>CTC": 0.038, "CCC>CGC": 0.012, "CCC>CAC": 0.006,
			"CCG>CTG": 0.055, "CCG>CGG": 0.018, "CCG>CAG": 0.010,
			"CCT>CTT": 0.085, "CCT>CGT": 0.028, "CCT>CAT": 0.015,
			"GCA>GTA": 0.012, "GCA>GAA": 0.008, "GCA>GGA": 0.004,
			"GCC>GTC": 0.015, "GCC>GAC": 0.010, "GCC>GGC": 0.005,
			"GCG>GTG": 0.018, "GCG>GAG": 0.012, "GCG>GGG": 0.006,
			"GCT>GTT": 0.025, "GCT>GAT": 0.016, "GCT>GGT": 0.008,
			"TCA>TTA": 0.165, "TCA>TGA": 0.055, "TCA>TAA": 0.028,
			"TCC>TTC": 0.142, "TCC>TGC": 0.048, "TCC>TAC": 0.024,
			"TCG>TTG": 0.125, "TCG>TGG": 0.042, "TCG>TAG": 0.021,
			"TCT>TTT": 0.185, "TCT>TGT": 0.062, "TCT>TAT": 0.031,
		},
	}
	msm.signatures["Signature_2"] = sig2

	// Signature 4: Tobacco smoking
	sig4 := &MutationalSignature{
		Name:        "Signature_4",
		Process:     TobaccoSmoking,
		Etiology:    "Tobacco smoking - exposure to benzo[a]pyrene diol epoxide",
		CancerTypes: []string{"Lung", "Head_neck", "Esophageal", "Bladder"},
		TrinucleotideProfile: map[string]float64{
			"ACA>ATA": 0.008, "ACA>AGA": 0.012, "ACA>AAA": 0.005,
			"ACC>ATC": 0.010, "ACC>AGC": 0.015, "ACC>AAC": 0.006,
			"ACG>ATG": 0.006, "ACG>AGG": 0.025, "ACG>AAG": 0.004,
			"ACT>ATT": 0.012, "ACT>AGT": 0.028, "ACT>AAT": 0.008,
			"CCA>CTA": 0.015, "CCA>CGA": 0.035, "CCA>CAA": 0.012,
			"CCC>CTC": 0.018, "CCC>CGC": 0.042, "CCC>CAC": 0.015,
			"CCG>CTG": 0.025, "CCG>CGG": 0.085, "CCG>CAG": 0.028,
			"CCT>CTT": 0.032, "CCT>CGT": 0.125, "CCT>CAT": 0.045,
			"GCA>GTA": 0.028, "GCA>GAA": 0.015, "GCA>GGA": 0.008,
			"GCC>GTC": 0.035, "GCC>GAC": 0.018, "GCC>GGC": 0.010,
			"GCG>GTG": 0.055, "GCG>GAG": 0.025, "GCG>GGG": 0.012,
			"GCT>GTT": 0.078, "GCT>GAT": 0.035, "GCT>GGT": 0.018,
			"TCA>TTA": 0.045, "TCA>TGA": 0.025, "TCA>TAA": 0.012,
			"TCC>TTC": 0.038, "TCC>TGC": 0.022, "TCC>TAC": 0.010,
			"TCG>TTG": 0.032, "TCG>TGG": 0.018, "TCG>TAG": 0.008,
			"TCT>TTT": 0.055, "TCT>TGT": 0.028, "TCT>TAT": 0.015,
		},
	}
	msm.signatures["Signature_4"] = sig4

	// Signature 7: UV exposure
	sig7 := &MutationalSignature{
		Name:        "Signature_7",
		Process:     UVExposure,
		Etiology:    "UV light exposure - formation of pyrimidine dimers",
		CancerTypes: []string{"Melanoma", "Skin"},
		TrinucleotideProfile: map[string]float64{
			"ACA>ATA": 0.005, "ACA>AGA": 0.003, "ACA>AAA": 0.002,
			"ACC>ATC": 0.008, "ACC>AGC": 0.005, "ACC>AAC": 0.003,
			"ACG>ATG": 0.004, "ACG>AGG": 0.003, "ACG>AAG": 0.001,
			"ACT>ATT": 0.012, "ACT>AGT": 0.008, "ACT>AAT": 0.004,
			"CCA>CTA": 0.025, "CCA>CGA": 0.008, "CCA>CAA": 0.004,
			"CCC>CTC": 0.185, "CCC>CGC": 0.028, "CCC>CAC": 0.015,
			"CCG>CTG": 0.015, "CCG>CGG": 0.005, "CCG>CAG": 0.003,
			"CCT>CTT": 0.285, "CCT>CGT": 0.045, "CCT>CAT": 0.025,
			"GCA>GTA": 0.008, "GCA>GAA": 0.005, "GCA>GGA": 0.003,
			"GCC>GTC": 0.012, "GCC>GAC": 0.008, "GCC>GGC": 0.004,
			"GCG>GTG": 0.006, "GCG>GAG": 0.004, "GCG>GGG": 0.002,
			"GCT>GTT": 0.018, "GCT>GAT": 0.012, "GCT>GGT": 0.006,
			"TCA>TTA": 0.035, "TCA>TGA": 0.012, "TCA>TAA": 0.006,
			"TCC>TTC": 0.125, "TCC>TGC": 0.042, "TCC>TAC": 0.021,
			"TCG>TTG": 0.015, "TCG>TGG": 0.005, "TCG>TAG": 0.003,
			"TCT>TTT": 0.165, "TCT>TGT": 0.055, "TCT>TAT": 0.028,
		},
	}
	msm.signatures["Signature_7"] = sig7

	// Signature 15: DNA mismatch repair deficiency
	sig15 := &MutationalSignature{
		Name:        "Signature_15",
		Process:     DNARepairDeficiency,
		Etiology:    "Defective DNA mismatch repair",
		CancerTypes: []string{"Colorectal", "Endometrial", "Stomach"},
		TrinucleotideProfile: map[string]float64{
			"ACA>ATA": 0.015, "ACA>AGA": 0.010, "ACA>AAA": 0.008,
			"ACC>ATC": 0.018, "ACC>AGC": 0.012, "ACC>AAC": 0.010,
			"ACG>ATG": 0.025, "ACG>AGG": 0.018, "ACG>AAG": 0.012,
			"ACT>ATT": 0.032, "ACT>AGT": 0.025, "ACT>AAT": 0.015,
			"CCA>CTA": 0.028, "CCA>CGA": 0.015, "CCA>CAA": 0.012,
			"CCC>CTC": 0.035, "CCC>CGC": 0.018, "CCC>CAC": 0.015,
			"CCG>CTG": 0.045, "CCG>CGG": 0.025, "CCG>CAG": 0.018,
			"CCT>CTT": 0.055, "CCT>CGT": 0.032, "CCT>CAT": 0.025,
			"GCA>GTA": 0.022, "GCA>GAA": 0.015, "GCA>GGA": 0.010,
			"GCC>GTC": 0.028, "GCC>GAC": 0.018, "GCC>GGC": 0.012,
			"GCG>GTG": 0.035, "GCG>GAG": 0.022, "GCG>GGG": 0.015,
			"GCT>GTT": 0.042, "GCT>GAT": 0.028, "GCT>GGT": 0.018,
			"TCA>TTA": 0.065, "TCA>TGA": 0.035, "TCA>TAA": 0.025,
			"TCC>TTC": 0.078, "TCC>TGC": 0.042, "TCC>TAC": 0.032,
			"TCG>TTG": 0.085, "TCG>TGG": 0.045, "TCG>TAG": 0.035,
			"TCT>TTT": 0.095, "TCT>TGT": 0.055, "TCT>TAT": 0.042,
		},
	}
	msm.signatures["Signature_15"] = sig15

	return nil
}

// setProcessActivity sets activity levels for specified mutational processes
func (msm *MutationalSignatureModel) setProcessActivity(processes []MutationalProcess) {
	// Default activity levels
	defaultActivity := map[MutationalProcess]float64{
		ClockLikeMutation:      0.3, // Always present at moderate level
		SpontaneousDeamination: 0.2, // Background level
		UVExposure:             0.0, // Only in skin cancers
		TobaccoSmoking:         0.0, // Only with exposure
		DNARepairDeficiency:    0.0, // Only when deficient
		APOBECActivity:         0.1, // Low background
		ChemicalMutagenesis:    0.0, // Only with exposure
		AgeRelatedMutation:     0.4, // Increases with age
	}

	// Set activities for all processes
	for process, activity := range defaultActivity {
		msm.processActivity[process.String()] = activity
	}

	// Increase activity for specified processes
	for _, process := range processes {
		msm.processActivity[process.String()] = 0.8 // High activity
	}
}

// AnalyzeSequence scores a sequence's trinucleotide composition against the
// signature profiles. The profiles are keyed by substitution ("ACA>ATA") and
// the observed frequencies by trinucleotide ("ACA"), so the cosine
// similarities, signature contributions and dominant signatures are always
// zero; the score comes from the composition, hotspot and bias terms only.
// No mutations are called, since a single sequence has none to call.
func (msm *MutationalSignatureModel) AnalyzeSequence(sequence []byte) (float64, *MutationalSignatureDetails, error) {
	if len(sequence) < 50 {
		return 0.0, nil, fmt.Errorf("sequence too short for mutational signature analysis")
	}

	// Extract trinucleotide context frequencies
	trinucleotideFreq := msm.extractTrinucleotideFrequencies(sequence)

	// Calculate signature contributions
	signatureContributions := msm.calculateSignatureContributions(trinucleotideFreq)

	// Identify dominant signatures
	dominantSignatures := msm.identifyDominantSignatures(signatureContributions)

	// Calculate mutational load
	mutationalLoad := msm.calculateMutationalLoad(sequence)

	// Assess signature consistency
	signatureConsistency := msm.assessSignatureConsistency(signatureContributions)

	// Calculate trinucleotide bias
	trinucleotideBias := msm.calculateTrinucleotideBias(trinucleotideFreq)

	details := &MutationalSignatureDetails{
		TrinucleotideFrequencies: trinucleotideFreq,
		SignatureContributions:   signatureContributions,
		DominantSignatures:       dominantSignatures,
		MutationalLoad:           mutationalLoad,
		SignatureConsistency:     signatureConsistency,
		TrinucleotideBias:        trinucleotideBias,
		SequenceLength:           int64(len(sequence)),
	}

	// Calculate overall mutational signature realism
	realismScore := msm.calculateMutationalRealism(details)

	return realismScore, details, nil
}

// extractTrinucleotideFrequencies extracts trinucleotide context frequencies
func (msm *MutationalSignatureModel) extractTrinucleotideFrequencies(sequence []byte) map[string]float64 {
	trinucleotideCount := make(map[string]int)
	totalTrinucleotides := 0

	// Extract all trinucleotides
	for i := 0; i <= len(sequence)-3; i++ {
		trinucleotide := string(sequence[i : i+3])
		trinucleotide = toUpper(trinucleotide)

		// Only count valid trinucleotides
		if isValidTrinucleotide(trinucleotide) {
			trinucleotideCount[trinucleotide]++
			totalTrinucleotides++
		}
	}

	// Convert counts to frequencies
	frequencies := make(map[string]float64)
	for trinucleotide, count := range trinucleotideCount {
		frequencies[trinucleotide] = float64(count) / float64(totalTrinucleotides)
	}

	return frequencies
}

// calculateSignatureContributions calculates the contribution of each signature
func (msm *MutationalSignatureModel) calculateSignatureContributions(observedFreq map[string]float64) map[string]float64 {
	contributions := make(map[string]float64)

	for sigName, signature := range msm.signatures {
		// Calculate similarity between observed and signature profile
		similarity := msm.calculateProfileSimilarity(observedFreq, signature.TrinucleotideProfile)

		// Weight by process activity
		processActivity := msm.processActivity[signature.Process.String()]
		weightedContribution := similarity * processActivity

		contributions[sigName] = weightedContribution
	}

	return contributions
}

// calculateProfileSimilarity calculates similarity between two trinucleotide profiles
func (msm *MutationalSignatureModel) calculateProfileSimilarity(profile1, profile2 map[string]float64) float64 {
	// Calculate cosine similarity
	dotProduct := 0.0
	norm1 := 0.0
	norm2 := 0.0

	// Get all trinucleotides present in either profile
	allTrinucleotides := make(map[string]bool)
	for trinuc := range profile1 {
		allTrinucleotides[trinuc] = true
	}
	for trinuc := range profile2 {
		allTrinucleotides[trinuc] = true
	}

	for trinuc := range allTrinucleotides {
		freq1 := profile1[trinuc]
		freq2 := profile2[trinuc]

		dotProduct += freq1 * freq2
		norm1 += freq1 * freq1
		norm2 += freq2 * freq2
	}

	if norm1 == 0 || norm2 == 0 {
		return 0.0
	}

	return dotProduct / (math.Sqrt(norm1) * math.Sqrt(norm2))
}

// identifyDominantSignatures identifies the most prominent signatures
func (msm *MutationalSignatureModel) identifyDominantSignatures(contributions map[string]float64) []string {
	type sigContrib struct {
		signature    string
		contribution float64
	}

	var sortedContribs []sigContrib
	for sig, contrib := range contributions {
		sortedContribs = append(sortedContribs, sigContrib{sig, contrib})
	}

	// Sort by contribution (descending)
	for i := 0; i < len(sortedContribs)-1; i++ {
		for j := i + 1; j < len(sortedContribs); j++ {
			if sortedContribs[i].contribution < sortedContribs[j].contribution {
				sortedContribs[i], sortedContribs[j] = sortedContribs[j], sortedContribs[i]
			}
		}
	}

	// Return top signatures with contribution > threshold
	var dominant []string
	threshold := 0.1

	for _, sc := range sortedContribs {
		if sc.contribution > threshold {
			dominant = append(dominant, sc.signature)
		}
	}

	return dominant
}

// calculateMutationalLoad returns the density of CpG, TT and CC
// dinucleotides (candidate hotspots), not an observed mutation burden.
func (msm *MutationalSignatureModel) calculateMutationalLoad(sequence []byte) float64 {
	if len(sequence) <= 1 {
		return 0.0
	}

	// Simple proxy: count potential mutation hotspots
	hotspots := 0

	for i := 0; i <= len(sequence)-2; i++ {
		dinucleotide := string(sequence[i : i+2])
		dinucleotide = toUpper(dinucleotide)

		// CpG dinucleotides are mutation hotspots
		if dinucleotide == "CG" {
			hotspots++
		}

		// TT and CC (potential UV damage sites)
		if dinucleotide == "TT" || dinucleotide == "CC" {
			hotspots++
		}
	}

	// Normalize by sequence length
	return float64(hotspots) / float64(len(sequence)-1)
}

// assessSignatureConsistency assesses consistency of signature contributions
func (msm *MutationalSignatureModel) assessSignatureConsistency(contributions map[string]float64) float64 {
	if len(contributions) == 0 {
		return 0.0
	}

	// Calculate coefficient of variation
	sum := 0.0
	count := 0

	for _, contrib := range contributions {
		sum += contrib
		count++
	}

	if count == 0 || sum == 0 {
		return 0.0
	}

	// Single signature: CV is undefined; perfect consistency by convention.
	if count < 2 {
		return 1.0
	}

	mean := sum / float64(count)

	// Calculate sample variance (Bessel's correction: ÷N-1, not ÷N).
	variance := 0.0
	for _, contrib := range contributions {
		diff := contrib - mean
		variance += diff * diff
	}
	variance /= float64(count - 1)

	stddev := math.Sqrt(variance)

	// Coefficient of variation
	cv := stddev / mean

	// Convert to consistency score (0-1, higher = more consistent)
	return math.Exp(-cv)
}

// calculateTrinucleotideBias calculates trinucleotide composition bias using a
// chi-square goodness-of-fit test against a uniform null distribution.
//
// All 64 possible trinucleotides are included in the statistic — even those
// absent from trinucleotideFreq — so that sequences dominated by few
// trinucleotides (high real bias) are not artificially deflated.
// The expected proportion per trinucleotide under the null is 1/64.
// df = 63; bias score is 0 for a uniform distribution and approaches 1 for a
// maximally biased distribution.
func (msm *MutationalSignatureModel) calculateTrinucleotideBias(trinucleotideFreq map[string]float64) float64 {
	const nCategories = 64
	expectedPerCategory := 1.0 / float64(nCategories)

	// Sum chi-square over all 64 trinucleotides, including unobserved ones.
	chiSquare := 0.0
	for _, b1 := range []byte("ACGT") {
		for _, b2 := range []byte("ACGT") {
			for _, b3 := range []byte("ACGT") {
				tri := string([]byte{b1, b2, b3})
				obs := trinucleotideFreq[tri] // 0.0 if absent
				diff := obs - expectedPerCategory
				chiSquare += (diff * diff) / expectedPerCategory
			}
		}
	}

	// Normalise by degrees of freedom (df = nCategories - 1 = 63).
	normalizedChiSquare := chiSquare / float64(nCategories-1)

	// Convert to bias score in [0, 1]; higher = more biased.
	return 1.0 - math.Exp(-normalizedChiSquare/10.0)
}

// calculateMutationalRealism calculates overall mutational signature realism
func (msm *MutationalSignatureModel) calculateMutationalRealism(details *MutationalSignatureDetails) float64 {
	// Weighted combination of mutational signature features
	weights := map[string]float64{
		"signature_consistency": 0.3,
		"dominant_signatures":   0.25,
		"mutational_load":       0.2,
		"trinucleotide_bias":    0.15,
		"profile_completeness":  0.1,
	}

	// Calculate score for dominant signatures
	dominantScore := 0.0
	if len(details.DominantSignatures) > 0 {
		dominantScore = math.Min(float64(len(details.DominantSignatures))/3.0, 1.0)
	}

	// Calculate profile completeness
	profileCompleteness := float64(len(details.TrinucleotideFrequencies)) / 64.0

	// Combine scores
	realismScore := 0.0
	realismScore += details.SignatureConsistency * weights["signature_consistency"]
	realismScore += dominantScore * weights["dominant_signatures"]
	realismScore += math.Min(details.MutationalLoad*5.0, 1.0) * weights["mutational_load"]
	realismScore += details.TrinucleotideBias * weights["trinucleotide_bias"]
	realismScore += profileCompleteness * weights["profile_completeness"]

	return realismScore
}

// MutationalSignatureDetails contains detailed mutational signature analysis results
type MutationalSignatureDetails struct {
	TrinucleotideFrequencies map[string]float64
	SignatureContributions   map[string]float64
	DominantSignatures       []string
	MutationalLoad           float64
	SignatureConsistency     float64
	TrinucleotideBias        float64
	SequenceLength           int64
}

// isValidTrinucleotide checks if a trinucleotide contains only valid DNA bases
func isValidTrinucleotide(trinuc string) bool {
	if len(trinuc) != 3 {
		return false
	}

	validBases := map[byte]bool{'A': true, 'T': true, 'G': true, 'C': true}

	for i := 0; i < 3; i++ {
		if !validBases[trinuc[i]] {
			return false
		}
	}

	return true
}

// String returns a string representation of the mutational process
func (mp MutationalProcess) String() string {
	processNames := map[MutationalProcess]string{
		SpontaneousDeamination: "Spontaneous Deamination",
		UVExposure:             "UV Exposure",
		TobaccoSmoking:         "Tobacco Smoking",
		DNARepairDeficiency:    "DNA Repair Deficiency",
		ClockLikeMutation:      "Clock-like Mutation",
		APOBECActivity:         "APOBEC Activity",
		ChemicalMutagenesis:    "Chemical Mutagenesis",
		AgeRelatedMutation:     "Age-related Mutation",
	}

	if name, exists := processNames[mp]; exists {
		return name
	}

	return "Unknown"
}
