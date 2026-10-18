package features

import (
	"fmt"
	"math"
	"math/rand/v2"
)

// OpenReadingFrame represents a protein-coding sequence
type OpenReadingFrame struct {
	Start      int    // Start position (0-based)
	End        int    // End position (0-based)
	Frame      int    // Reading frame (0, 1, or 2)
	Strand     string // "0" (forward) or "1" (reverse complement)
	Sequence   []byte // Nucleotide sequence
	Protein    []byte // Translated protein sequence
	Length     int    // Length in codons
	GCContent  float64
	StartCodon string         // Start codon (ATG, etc.)
	StopCodon  string         // Stop codon (TAA, TAG, TGA)
	CodonCount map[string]int // Codon frequency counts
	CAI        float64        // Codon Adaptation Index
}

// GeneticCode represents the genetic code mapping
type GeneticCode struct {
	CodonToAA   map[string]string   // Codon -> amino acid
	AAToCodon   map[string][]string // Amino acid -> codons
	StartCodons []string            // Start codons
	StopCodons  []string            // Stop codons
}

// CodonUsageTable represents codon usage frequencies for an organism
type CodonUsageTable struct {
	Organism            string
	Tissue              string
	GeneCount           int
	CodonFrequencies    map[string]float64
	RelativeFrequencies map[string]float64 // Relative to synonymous codons
	OptimalCodons       map[string]string  // Most frequent codon per amino acid
}

// CodonBiasResult represents the result of applying codon bias
type CodonBiasResult struct {
	OriginalSequence []byte
	BiasedSequence   []byte
	CodonChanges     []CodonChange
	ExpressionLevel  float64
	BiasStrength     float64
	ORFs             []OpenReadingFrame
}

// CodonChange represents a codon substitution
type CodonChange struct {
	Position      int
	OriginalCodon string
	NewCodon      string
	AminoAcid     string
	ChangeType    string // "optimal", "frequent", "rare"
	tRNAAbundance float64
}

// ORF represents an ORF in the sequence
type ORF struct {
	Start      int
	End        int
	Strand     string
	StartCodon string
	StopCodon  string
	Length     int
	GCContent  float64
	CodonCount map[string]int
	CAI        float64 // Codon Adaptation Index
}

// getStandardGeneticCode returns the standard genetic code
func getStandardGeneticCode() *GeneticCode {
	gc := &GeneticCode{
		CodonToAA: map[string]string{
			"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
			"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
			"TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
			"TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
			"CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
			"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
			"CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
			"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
			"ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
			"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
			"AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
			"AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
			"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
			"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
			"GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
			"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
		},
		AAToCodon:   make(map[string][]string),
		StartCodons: []string{"ATG"},
		StopCodons:  []string{"TAA", "TAG", "TGA"},
	}

	// Build reverse mapping
	for codon, aa := range gc.CodonToAA {
		gc.AAToCodon[aa] = append(gc.AAToCodon[aa], codon)
	}

	return gc
}

// loadDefaultCodonUsage loads default codon usage frequencies for the
// organism. NewCodonUsageModel does not call it, so codonFrequencies (used by
// ApplyCodonBias and GetOptimalCodons) is empty until this or LoadCodonTable
// runs.
func (m *CodonUsageModel) loadDefaultCodonUsage() {
	switch m.organism {
	case "human", "homo_sapiens":
		m.loadHumanCodonUsage()
	case "mouse", "mus_musculus":
		m.loadMouseCodonUsage()
	case "yeast", "saccharomyces_cerevisiae":
		m.loadYeastCodonUsage()
	case "ecoli", "escherichia_coli":
		m.loadEcoliCodonUsage()
	default:
		m.loadUniversalCodonUsage()
	}
}

// loadHumanCodonUsage loads human codon usage frequencies.
// Values are per-thousand-codons from the Kazusa Codon Usage Database
// (Nakamura et al. 2000, NAR 28:292; organism: Homo sapiens, all GenBank
// coding sequences, not only highly expressed genes).
func (m *CodonUsageModel) loadHumanCodonUsage() {
	humanUsage := map[string]float64{
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

	// Normalize to frequencies
	total := 0.0
	for _, freq := range humanUsage {
		total += freq
	}

	for codon, freq := range humanUsage {
		m.codonFrequencies[codon] = freq / total
	}

	// Build amino acid mapping
	m.buildAminoAcidMapping()
}

// loadMouseCodonUsage loads approximate mouse codon usage frequencies. The
// values are the human table with small hand-made adjustments, not a
// published mouse table.
func (m *CodonUsageModel) loadMouseCodonUsage() {
	mouseUsage := map[string]float64{
		"TTT": 17.1, "TTC": 20.8, "TTA": 7.4, "TTG": 12.6,
		"TCT": 15.5, "TCC": 17.4, "TCA": 12.5, "TCG": 4.6,
		"TAT": 12.5, "TAC": 15.0, "TAA": 1.1, "TAG": 0.9,
		"TGT": 10.8, "TGC": 12.3, "TGA": 1.5, "TGG": 13.5,
		"CTT": 13.5, "CTC": 19.3, "CTA": 7.5, "CTG": 39.2,
		"CCT": 17.8, "CCC": 19.5, "CCA": 17.2, "CCG": 7.2,
		"CAT": 11.2, "CAC": 14.8, "CAA": 12.6, "CAG": 33.9,
		"CGT": 4.8, "CGC": 10.1, "CGA": 6.5, "CGG": 11.1,
		"ATT": 16.3, "ATC": 20.5, "ATA": 7.8, "ATG": 21.7,
		"ACT": 13.4, "ACC": 18.6, "ACA": 15.4, "ACG": 6.4,
		"AAT": 17.3, "AAC": 18.8, "AAA": 24.7, "AAG": 31.6,
		"AGT": 12.4, "AGC": 19.2, "AGA": 12.5, "AGG": 11.7,
		"GTT": 11.3, "GTC": 14.2, "GTA": 7.4, "GTG": 27.8,
		"GCT": 18.7, "GCC": 27.4, "GCA": 16.1, "GCG": 7.7,
		"GAT": 22.1, "GAC": 24.8, "GAA": 29.3, "GAG": 39.3,
		"GGT": 11.1, "GGC": 21.9, "GGA": 16.8, "GGG": 16.2,
	}

	// Normalize and store
	total := 0.0
	for _, freq := range mouseUsage {
		total += freq
	}
	for codon, freq := range mouseUsage {
		m.codonFrequencies[codon] = freq / total
	}

	m.buildAminoAcidMapping()
}

// loadUniversalCodonUsage loads universal (balanced) codon usage
func (m *CodonUsageModel) loadUniversalCodonUsage() {
	// Equal frequency for all synonymous codons
	for aa, codons := range m.codonTable.AAToCodon {
		if aa == "*" {
			continue // Skip stop codons
		}
		frequency := 1.0 / float64(len(codons))
		for _, codon := range codons {
			m.codonFrequencies[codon] = frequency / 64.0 // Normalize across all codons
		}
	}

	m.buildAminoAcidMapping()
}

// buildAminoAcidMapping builds the amino acid to codon mapping
func (m *CodonUsageModel) buildAminoAcidMapping() {
	m.aminoAcidMapping = make(map[string][]string)

	for codon, aa := range m.codonTable.CodonToAA {
		if aa != "*" { // Skip stop codons
			m.aminoAcidMapping[aa] = append(m.aminoAcidMapping[aa], codon)
		}
	}
}

// LoadCodonTable loads a custom codon usage table
func (m *CodonUsageModel) LoadCodonTable(table *CodonUsageTable) {
	m.codonFrequencies = make(map[string]float64)
	for codon, freq := range table.CodonFrequencies {
		m.codonFrequencies[codon] = freq
	}
	m.buildAminoAcidMapping()
}

// SetExpressionLevel sets the gene expression level (affects codon bias strength)
func (m *CodonUsageModel) SetExpressionLevel(level float64) {
	if level >= 0 {
		m.expressionLevel = level
	}
}

// ApplyCodonBias finds ORFs (ATG to stop, at least 150 bp) in all six frames
// and resamples each internal codon from its synonymous codons with
// probability proportional to frequency^(1+expressionLevel). Reverse-strand
// ORF coordinates refer to the reverse complement but are applied to the
// forward sequence.
func (m *CodonUsageModel) ApplyCodonBias(sequence []byte) (*CodonBiasResult, error) {
	if len(sequence) == 0 {
		return nil, fmt.Errorf("sequence cannot be empty")
	}

	result := &CodonBiasResult{
		OriginalSequence: make([]byte, len(sequence)),
		BiasedSequence:   make([]byte, len(sequence)),
		CodonChanges:     make([]CodonChange, 0),
		ExpressionLevel:  m.expressionLevel,
	}

	copy(result.OriginalSequence, sequence)
	copy(result.BiasedSequence, sequence)

	// Find open reading frames
	orfs := m.findOpenReadingFrames(sequence)
	result.ORFs = orfs

	// Apply codon bias to each ORF
	for _, orf := range orfs {
		changes := m.applyCodonBiasToORF(result.BiasedSequence, orf)
		result.CodonChanges = append(result.CodonChanges, changes...)
	}

	// Calculate bias strength
	result.BiasStrength = m.calculateBiasStrength(result.CodonChanges)

	return result, nil
}

// findOpenReadingFrames finds ORFs in the sequence
func (m *CodonUsageModel) findOpenReadingFrames(sequence []byte) []OpenReadingFrame {
	orfs := make([]OpenReadingFrame, 0)

	// Search in all three reading frames on both strands
	for strand := 0; strand < 2; strand++ {
		var searchSeq []byte
		if strand == 0 {
			searchSeq = sequence
		} else {
			searchSeq = m.reverseComplement(sequence)
		}

		for frame := 0; frame < 3; frame++ {
			frameORFs := m.findORFsInFrame(searchSeq, frame, strand)
			orfs = append(orfs, frameORFs...)
		}
	}

	return orfs
}

// findORFsInFrame finds ORFs in a specific reading frame
func (m *CodonUsageModel) findORFsInFrame(sequence []byte, frame int, strand int) []OpenReadingFrame {
	orfs := make([]OpenReadingFrame, 0)

	if len(sequence) < frame+3 {
		return orfs
	}

	for i := frame; i < len(sequence)-2; i += 3 {
		codon := string(sequence[i : i+3])

		// Check for start codon
		if m.isStartCodon(codon) {
			orf := m.extendORF(sequence, i, strand)
			if orf.Length >= 150 { // Minimum ORF length (50 codons)
				orfs = append(orfs, orf)
			}
		}
	}

	return orfs
}

// extendORF extends an ORF from a start codon to a stop codon
func (m *CodonUsageModel) extendORF(sequence []byte, start int, strand int) OpenReadingFrame {
	orf := OpenReadingFrame{
		Start:      start,
		Strand:     fmt.Sprintf("%d", strand),
		StartCodon: string(sequence[start : start+3]),
		CodonCount: make(map[string]int),
	}

	// Extend until stop codon or end of sequence
	for i := start; i < len(sequence)-2; i += 3 {
		codon := string(sequence[i : i+3])
		orf.CodonCount[codon]++

		if m.isStopCodon(codon) {
			orf.End = i + 3
			orf.StopCodon = codon
			break
		}
	}

	// If no stop codon found, extend to end of sequence
	if orf.End == 0 {
		orf.End = len(sequence)
	}

	orf.Length = orf.End - orf.Start
	orf.GCContent = m.calculateGCContent(sequence[orf.Start:orf.End])
	orf.CAI = m.calculateCodonAdaptationIndex(orf.CodonCount)

	return orf
}

// applyCodonBiasToORF applies codon bias to a single ORF
func (m *CodonUsageModel) applyCodonBiasToORF(sequence []byte, orf OpenReadingFrame) []CodonChange {
	changes := make([]CodonChange, 0)

	for i := orf.Start; i < orf.End-2; i += 3 {
		codon := string(sequence[i : i+3])

		// Skip start and stop codons
		if m.isStartCodon(codon) || m.isStopCodon(codon) {
			continue
		}

		// Get amino acid
		aa := m.codonTable.CodonToAA[codon]
		if aa == "" {
			continue
		}

		// Choose optimal codon based on usage frequency and expression level
		newCodon := m.chooseOptimalCodon(aa, m.expressionLevel)

		if newCodon != codon {
			// Apply the change
			copy(sequence[i:i+3], []byte(newCodon))

			change := CodonChange{
				Position:      i,
				OriginalCodon: codon,
				NewCodon:      newCodon,
				AminoAcid:     aa,
				ChangeType:    m.getChangeType(codon, newCodon),
				tRNAAbundance: m.gettRNAAbundance(newCodon),
			}
			changes = append(changes, change)
		}
	}

	return changes
}

// chooseOptimalCodon chooses the optimal codon for an amino acid
func (m *CodonUsageModel) chooseOptimalCodon(aa string, expressionLevel float64) string {
	codons := m.aminoAcidMapping[aa]
	if len(codons) == 0 {
		return ""
	}

	// Calculate selection probabilities based on usage frequency and expression level
	probabilities := make([]float64, len(codons))
	totalProb := 0.0

	for i, codon := range codons {
		baseFreq := m.codonFrequencies[codon]

		// Adjust probability based on expression level
		// Higher expression -> stronger bias toward optimal codons
		adjustedFreq := math.Pow(baseFreq, 1.0+expressionLevel)
		probabilities[i] = adjustedFreq
		totalProb += adjustedFreq
	}

	// Normalize probabilities
	for i := range probabilities {
		probabilities[i] /= totalProb
	}

	// Choose codon based on probabilities
	r := m.rng.Float64()
	cumulative := 0.0

	for i, prob := range probabilities {
		cumulative += prob
		if r <= cumulative {
			return codons[i]
		}
	}

	return codons[len(codons)-1]
}

// getChangeType determines the type of codon change
func (m *CodonUsageModel) getChangeType(originalCodon, newCodon string) string {
	originalFreq := m.codonFrequencies[originalCodon]
	newFreq := m.codonFrequencies[newCodon]

	if newFreq > originalFreq*1.5 {
		return "optimal"
	} else if newFreq > originalFreq {
		return "frequent"
	} else {
		return "rare"
	}
}

// gettRNAAbundance returns the codon's usage frequency × 100 as a stand-in
// for tRNA abundance; no tRNA data is used.
func (m *CodonUsageModel) gettRNAAbundance(codon string) float64 {
	return m.codonFrequencies[codon] * 100.0
}

// calculateBiasStrength calculates the overall bias strength
func (m *CodonUsageModel) calculateBiasStrength(changes []CodonChange) float64 {
	if len(changes) == 0 {
		return 0.0
	}

	strengthSum := 0.0
	for _, change := range changes {
		originalFreq := m.codonFrequencies[change.OriginalCodon]
		newFreq := m.codonFrequencies[change.NewCodon]

		if originalFreq > 0 {
			strengthSum += newFreq / originalFreq
		}
	}

	return strengthSum / float64(len(changes))
}

// calculateCodonAdaptationIndex calculates the CAI for an ORF
func (m *CodonUsageModel) calculateCodonAdaptationIndex(codonCount map[string]int) float64 {
	if len(codonCount) == 0 {
		return 0.0
	}

	geometricMean := 1.0
	totalCodons := 0

	for codon, count := range codonCount {
		aa := m.codonTable.CodonToAA[codon]
		if aa == "" || aa == "*" {
			continue
		}

		// Find maximum frequency for this amino acid
		maxFreq := 0.0
		for _, synonymousCodon := range m.aminoAcidMapping[aa] {
			freq := m.codonFrequencies[synonymousCodon]
			if freq > maxFreq {
				maxFreq = freq
			}
		}

		// Calculate relative adaptiveness
		if maxFreq > 0 {
			relativeAdaptiveness := m.codonFrequencies[codon] / maxFreq
			geometricMean *= math.Pow(relativeAdaptiveness, float64(count))
			totalCodons += count
		}
	}

	if totalCodons == 0 {
		return 0.0
	}

	return math.Pow(geometricMean, 1.0/float64(totalCodons))
}

// calculateGCContent calculates GC content of a sequence
func (m *CodonUsageModel) calculateGCContent(sequence []byte) float64 {
	if len(sequence) == 0 {
		return 0.0
	}

	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}

	return float64(gcCount) / float64(len(sequence))
}

// isStartCodon checks if a codon is a start codon
func (m *CodonUsageModel) isStartCodon(codon string) bool {
	for _, start := range m.codonTable.StartCodons {
		if codon == start {
			return true
		}
	}
	return false
}

// isStopCodon checks if a codon is a stop codon
func (m *CodonUsageModel) isStopCodon(codon string) bool {
	for _, stop := range m.codonTable.StopCodons {
		if codon == stop {
			return true
		}
	}
	return false
}

// reverseComplement returns the reverse complement of a sequence
func (m *CodonUsageModel) reverseComplement(sequence []byte) []byte {
	complement := map[byte]byte{
		'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
		'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
		'N': 'N', 'n': 'n',
	}

	result := make([]byte, len(sequence))
	for i, base := range sequence {
		result[len(sequence)-1-i] = complement[base]
	}

	return result
}

// SetSeed sets the random seed for reproducible codon selection
func (m *CodonUsageModel) SetSeed(seed int64) {
	m.rng = rand.New(rand.NewPCG(uint64(seed), 0)) // #nosec G115 G404
}

// GetModelStats returns statistics about the codon usage model
func (m *CodonUsageModel) GetModelStats() map[string]interface{} {
	stats := map[string]interface{}{
		"organism":         m.organism,
		"tissue_specific":  m.tissueSpecific,
		"expression_level": m.expressionLevel,
		"total_codons":     len(m.codonFrequencies),
		"amino_acids":      len(m.aminoAcidMapping),
	}

	// Calculate entropy of codon usage
	entropy := 0.0
	for _, freq := range m.codonFrequencies {
		if freq > 0 {
			entropy -= freq * math.Log2(freq)
		}
	}
	stats["codon_entropy"] = entropy

	// Find most and least frequent codons
	maxFreq := 0.0
	minFreq := 1.0
	var mostFrequentCodon, leastFrequentCodon string

	for codon, freq := range m.codonFrequencies {
		if freq > maxFreq {
			maxFreq = freq
			mostFrequentCodon = codon
		}
		if freq < minFreq {
			minFreq = freq
			leastFrequentCodon = codon
		}
	}

	stats["most_frequent_codon"] = mostFrequentCodon
	stats["least_frequent_codon"] = leastFrequentCodon
	stats["max_frequency"] = maxFreq
	stats["min_frequency"] = minFreq

	return stats
}

// ValidateCodonUsage validates the codon usage model
func (m *CodonUsageModel) ValidateCodonUsage() error {
	// Check that all codons are present
	for codon := range m.codonTable.CodonToAA {
		if _, exists := m.codonFrequencies[codon]; !exists {
			return fmt.Errorf("missing frequency for codon %s", codon)
		}
	}

	// Check that frequencies sum to approximately 1
	total := 0.0
	for _, freq := range m.codonFrequencies {
		total += freq
	}

	if math.Abs(total-1.0) > 0.01 {
		return fmt.Errorf("codon frequencies do not sum to 1.0: %f", total)
	}

	return nil
}

// Additional utility functions for specific codon usage patterns

// loadYeastCodonUsage loads S. cerevisiae codon usage frequencies.
// Values are per-thousand-codons from the Kazusa Codon Usage Database
// (Nakamura et al. 2000, NAR 28:292; organism: Saccharomyces cerevisiae).
// Yeast has very strong codon bias: e.g. GAA (Glu) is 2.4× more common
// than GAG; CTG is only ~10‰ of all codons (vs. 52‰ in E. coli).
func (m *CodonUsageModel) loadYeastCodonUsage() {
	yeastUsage := map[string]float64{
		"TTT": 26.1, "TTC": 18.4, "TTA": 26.2, "TTG": 27.2,
		"CTT": 11.3, "CTC": 5.4, "CTA": 13.4, "CTG": 10.0,
		"ATT": 30.1, "ATC": 17.2, "ATA": 17.7, "ATG": 20.9,
		"GTT": 21.1, "GTC": 11.8, "GTA": 11.8, "GTG": 10.8,
		"TCT": 23.5, "TCC": 14.2, "TCA": 18.7, "TCG": 8.6,
		"CCT": 13.5, "CCC": 6.8, "CCA": 18.3, "CCG": 5.3,
		"ACT": 20.3, "ACC": 12.7, "ACA": 17.9, "ACG": 7.9,
		"GCT": 21.2, "GCC": 12.6, "GCA": 16.5, "GCG": 6.2,
		"TAT": 18.8, "TAC": 14.8, "TAA": 1.1, "TAG": 0.5,
		"CAT": 13.6, "CAC": 7.8, "CAA": 27.3, "CAG": 12.1,
		"AAT": 35.9, "AAC": 24.8, "AAA": 41.9, "AAG": 30.8,
		"GAT": 37.6, "GAC": 20.8, "GAA": 45.5, "GAG": 19.2,
		"TGT": 8.1, "TGC": 4.8, "TGA": 0.7, "TGG": 10.4,
		"CGT": 6.4, "CGC": 2.6, "CGA": 3.0, "CGG": 1.7,
		"AGT": 14.2, "AGC": 9.8, "AGA": 21.4, "AGG": 9.2,
		"GGT": 23.9, "GGC": 9.8, "GGA": 10.9, "GGG": 5.5,
	}

	total := 0.0
	for _, freq := range yeastUsage {
		total += freq
	}
	for codon, freq := range yeastUsage {
		m.codonFrequencies[codon] = freq / total
	}
	m.buildAminoAcidMapping()
}

// loadEcoliCodonUsage loads E. coli K-12 MG1655 codon usage frequencies.
// Values are per-thousand-codons from the Kazusa Codon Usage Database
// (Nakamura et al. 2000, NAR 28:292; organism: Escherichia coli K-12).
// E. coli has strong GC-biased codon preferences in highly expressed genes:
// CTG dominates Leu (~52‰), CCG dominates Pro (~23‰), GCG dominates Ala (~34‰),
// and CGT/CGC dominate Arg — patterns entirely absent from the universal table.
func (m *CodonUsageModel) loadEcoliCodonUsage() {
	ecoliUsage := map[string]float64{
		"TTT": 21.9, "TTC": 16.8, "TTA": 13.7, "TTG": 13.1,
		"CTT": 11.1, "CTC": 10.9, "CTA": 3.9, "CTG": 52.4,
		"ATT": 28.7, "ATC": 25.0, "ATA": 7.1, "ATG": 27.9,
		"GTT": 18.1, "GTC": 15.3, "GTA": 10.8, "GTG": 26.1,
		"TCT": 8.5, "TCC": 8.6, "TCA": 7.1, "TCG": 8.8,
		"CCT": 7.1, "CCC": 5.5, "CCA": 8.4, "CCG": 23.3,
		"ACT": 8.8, "ACC": 23.4, "ACA": 6.7, "ACG": 13.6,
		"GCT": 16.2, "GCC": 25.5, "GCA": 20.4, "GCG": 33.5,
		"TAT": 16.2, "TAC": 12.3, "TAA": 2.0, "TAG": 0.3,
		"CAT": 13.3, "CAC": 9.7, "CAA": 15.1, "CAG": 28.7,
		"AAT": 17.7, "AAC": 22.0, "AAA": 33.5, "AAG": 10.3,
		"GAT": 32.1, "GAC": 19.0, "GAA": 39.1, "GAG": 17.8,
		"TGT": 5.1, "TGC": 6.4, "TGA": 1.0, "TGG": 15.2,
		"CGT": 20.8, "CGC": 22.0, "CGA": 3.6, "CGG": 5.4,
		"AGT": 8.5, "AGC": 16.1, "AGA": 2.1, "AGG": 1.2,
		"GGT": 24.3, "GGC": 29.9, "GGA": 8.0, "GGG": 11.1,
	}

	total := 0.0
	for _, freq := range ecoliUsage {
		total += freq
	}
	for codon, freq := range ecoliUsage {
		m.codonFrequencies[codon] = freq / total
	}
	m.buildAminoAcidMapping()
}

// GetOptimalCodons returns the most optimal codons for each amino acid
func (m *CodonUsageModel) GetOptimalCodons() map[string]string {
	optimal := make(map[string]string)

	for aa, codons := range m.aminoAcidMapping {
		maxFreq := 0.0
		optimalCodon := ""

		for _, codon := range codons {
			freq := m.codonFrequencies[codon]
			if freq > maxFreq {
				maxFreq = freq
				optimalCodon = codon
			}
		}

		optimal[aa] = optimalCodon
	}

	return optimal
}
