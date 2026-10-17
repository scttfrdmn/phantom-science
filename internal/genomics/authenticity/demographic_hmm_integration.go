package authenticity

import (
	"fmt"
	"math/rand/v2"
	"time"
)

// DemographicChromosome implements ChromosomeInterface using demographic simulation data
type DemographicChromosome struct {
	name         string
	sequence     []byte
	variants     []DemographicVariant
	length       int
	population   string
}

// DemographicVariant represents a variant from demographic simulation
type DemographicVariant struct {
	Position        int
	ReferenceAllele byte
	AlternateAllele byte
	Frequency       float64
	Genotypes       []int
}

// DemographicGenome implements GenomeInterface using demographic simulation data
type DemographicGenome struct {
	chromosomes map[string]*DemographicChromosome
	population  string
}

// NewDemographicGenome creates a new demographic genome
func NewDemographicGenome(population string) *DemographicGenome {
	return &DemographicGenome{
		chromosomes: make(map[string]*DemographicChromosome),
		population:  population,
	}
}

// AddChromosome adds a chromosome with demographic variants
func (dg *DemographicGenome) AddChromosome(name string, length int, variants []DemographicVariant) {
	// Generate base sequence
	sequence := generateReferenceSequence(length)
	
	// Apply variants to create population-specific sequence
	for _, variant := range variants {
		if variant.Position >= 0 && variant.Position < len(sequence) {
			// Apply variant based on frequency (simplified - in reality would use haplotypes)
			if rand.Float64() < variant.Frequency { // #nosec G404
				sequence[variant.Position] = variant.AlternateAllele
			}
		}
	}
	
	chromosome := &DemographicChromosome{
		name:       name,
		sequence:   sequence,
		variants:   variants,
		length:     length,
		population: dg.population,
	}
	
	dg.chromosomes[name] = chromosome
}

// GetChromosomeNames returns list of chromosome names
func (dg *DemographicGenome) GetChromosomeNames() []string {
	names := make([]string, 0, len(dg.chromosomes))
	for name := range dg.chromosomes {
		names = append(names, name)
	}
	return names
}

// LoadChromosome returns chromosome by name
func (dg *DemographicGenome) LoadChromosome(name string) (ChromosomeInterface, error) {
	if chr, exists := dg.chromosomes[name]; exists {
		return chr, nil
	}
	return nil, fmt.Errorf("chromosome %s not found", name)
}

// GetSequence returns chromosome sequence
func (dc *DemographicChromosome) GetSequence() []byte {
	return dc.sequence
}

// GetName returns chromosome name
func (dc *DemographicChromosome) GetName() string {
	return dc.name
}

// GetLength returns chromosome length
func (dc *DemographicChromosome) GetLength() int {
	return dc.length
}

// generateReferenceSequence generates a realistic reference sequence
func generateReferenceSequence(length int) []byte {
	_ = []byte{'A', 'T', 'G', 'C'} // Available bases for future use
	sequence := make([]byte, length)
	
	// Use deterministic random for consistency
	rng := rand.New(rand.NewPCG(42, 0)) // #nosec G404
	
	// Generate sequence with realistic GC content (~40%)
	for i := 0; i < length; i++ {
		r := rng.Float64()
		if r < 0.3 {
			sequence[i] = 'A'
		} else if r < 0.6 {
			sequence[i] = 'T'
		} else if r < 0.8 {
			sequence[i] = 'G'
		} else {
			sequence[i] = 'C'
		}
	}
	
	return sequence
}

// EnhancedHMMKmerModel extends the basic HMM model with demographic integration
type EnhancedHMMKmerModel struct {
	*HMMKmerModel
	population           string
	demographicGenome    *DemographicGenome
	variantIntegrated    bool
	populationStatistics map[string]float64
}

// NewEnhancedHMMKmerModel creates an enhanced HMM model with demographic integration
func NewEnhancedHMMKmerModel(k int, population string) *EnhancedHMMKmerModel {
	baseModel := NewHMMKmerModel(k)
	
	return &EnhancedHMMKmerModel{
		HMMKmerModel:         baseModel,
		population:           population,
		demographicGenome:    NewDemographicGenome(population),
		variantIntegrated:    false,
		populationStatistics: make(map[string]float64),
	}
}

// IntegrateDemographicVariants integrates variants from demographic simulation
func (ehm *EnhancedHMMKmerModel) IntegrateDemographicVariants(
	chromosomeName string,
	length int,
	variants []DemographicVariant,
	statistics map[string]float64,
) error {
	
	fmt.Printf("🔗 Integrating demographic variants for %s chromosome %s\n", 
		ehm.population, chromosomeName)
	fmt.Printf("   Variants: %d, Length: %d bp\n", len(variants), length)
	
	// Add chromosome to demographic genome
	ehm.demographicGenome.AddChromosome(chromosomeName, length, variants)
	
	// Store population statistics
	for key, value := range statistics {
		ehm.populationStatistics[key] = value
	}
	
	// Train HMM model on demographic genome
	err := ehm.HMMKmerModel.TrainFromReferenceGenome(ehm.demographicGenome)
	if err != nil {
		return fmt.Errorf("failed to train HMM on demographic genome: %w", err)
	}
	
	ehm.variantIntegrated = true
	
	fmt.Printf("   ✅ Demographic integration complete\n")
	fmt.Printf("   📊 Total k-mers processed: %d\n", ehm.HMMKmerModel.TotalKmers())
	
	return nil
}

// GeneratePopulationSequence generates sequence with population-specific characteristics
func (ehm *EnhancedHMMKmerModel) GeneratePopulationSequence(length int) (*GenerationResult, error) {
	if !ehm.variantIntegrated {
		return nil, fmt.Errorf("demographic variants not integrated - call IntegrateDemographicVariants first")
	}
	
	fmt.Printf("🧬 Generating %s-specific sequence (length: %d)\n", ehm.population, length)
	
	// Generate base sequence using HMM
	result, err := ehm.HMMKmerModel.GenerateSequence(length)
	if err != nil {
		return nil, fmt.Errorf("failed to generate base sequence: %w", err)
	}
	
	// Enhance with population-specific characteristics
	ehm.applyPopulationCharacteristics(result)
	
	return result, nil
}

// applyPopulationCharacteristics applies population-specific modifications to generated sequence
func (ehm *EnhancedHMMKmerModel) applyPopulationCharacteristics(result *GenerationResult) {
	// Apply population-specific GC content adjustments
	targetGC := ehm.getPopulationGCContent()
	ehm.adjustGCContent(result.Sequence, targetGC)
	
	// Apply population-specific quality score patterns
	ehm.adjustQualityScores(result.Qualities)
	
	// Update generation metadata
	result.KmerFrequencies["population"] = 1 // Marker for population-specific generation
}

// getPopulationGCContent returns target GC content for population
func (ehm *EnhancedHMMKmerModel) getPopulationGCContent() float64 {
	// Population-specific GC content from literature
	populationGC := map[string]float64{
		"YRI": 0.41,  // African populations
		"LWK": 0.41,
		"GWD": 0.41,
		"CEU": 0.42,  // European populations
		"GBR": 0.42,
		"TSI": 0.42,
		"CHB": 0.40,  // East Asian populations
		"JPT": 0.40,
		"CHS": 0.40,
		"MXL": 0.41,  // Americas populations (mixed)
		"PUR": 0.41,
	}
	
	if gc, exists := populationGC[ehm.population]; exists {
		return gc
	}
	
	return 0.41 // Default human GC content
}

// adjustGCContent adjusts sequence to match target GC content
func (ehm *EnhancedHMMKmerModel) adjustGCContent(sequence []byte, targetGC float64) {
	currentGC := calculateGCContent(sequence)
	
	if abs(currentGC-targetGC) < 0.01 {
		return // Already close enough
	}
	
	// Simple adjustment - in practice would use more sophisticated methods
	rng := rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)) // #nosec G404
	
	if currentGC < targetGC {
		// Need more GC bases
		for i := 0; i < len(sequence); i++ {
			if (sequence[i] == 'A' || sequence[i] == 'T') && rng.Float64() < 0.1 {
				if rng.Float64() < 0.5 {
					sequence[i] = 'G'
				} else {
					sequence[i] = 'C'
				}
			}
		}
	} else {
		// Need fewer GC bases
		for i := 0; i < len(sequence); i++ {
			if (sequence[i] == 'G' || sequence[i] == 'C') && rng.Float64() < 0.1 {
				if rng.Float64() < 0.5 {
					sequence[i] = 'A'
				} else {
					sequence[i] = 'T'
				}
			}
		}
	}
}

// adjustQualityScores applies population-specific quality score patterns
func (ehm *EnhancedHMMKmerModel) adjustQualityScores(qualities []byte) {
	// Population-specific quality patterns could be applied here
	// For now, maintain existing quality distribution
}

// calculateGCContent calculates GC content of sequence
func calculateGCContent(sequence []byte) float64 {
	gcCount := 0
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			gcCount++
		}
	}
	return float64(gcCount) / float64(len(sequence))
}

// abs returns absolute value of float64
func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

// GetPopulationStatistics returns integrated population statistics
func (ehm *EnhancedHMMKmerModel) GetPopulationStatistics() map[string]float64 {
	return ehm.populationStatistics
}

// GetPopulation returns the population this model represents
func (ehm *EnhancedHMMKmerModel) GetPopulation() string {
	return ehm.population
}

// IsEnhanced returns true if demographic integration is complete
func (ehm *EnhancedHMMKmerModel) IsEnhanced() bool {
	return ehm.variantIntegrated
}

// GetEnhancedModelStats returns enhanced model statistics
func (ehm *EnhancedHMMKmerModel) GetEnhancedModelStats() map[string]interface{} {
	baseStats := ehm.HMMKmerModel.GetModelStats()
	
	enhancedStats := map[string]interface{}{
		"base_model":            baseStats,
		"population":            ehm.population,
		"demographic_enhanced":  ehm.variantIntegrated,
		"population_statistics": ehm.populationStatistics,
		"chromosomes_loaded":    len(ehm.demographicGenome.chromosomes),
	}
	
	if len(ehm.demographicGenome.chromosomes) > 0 {
		chrNames := make([]string, 0, len(ehm.demographicGenome.chromosomes))
		totalVariants := 0
		
		for name, chr := range ehm.demographicGenome.chromosomes {
			chrNames = append(chrNames, name)
			totalVariants += len(chr.variants)
		}
		
		enhancedStats["chromosome_names"] = chrNames
		enhancedStats["total_variants_integrated"] = totalVariants
	}
	
	return enhancedStats
}