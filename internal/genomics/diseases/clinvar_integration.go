package diseases

import (
	"context"
	"fmt"
	"io"
	"math/rand"
	"net/http"
	"net/url"
	"strings"
	"time"
)

// ClinVarIntegrator handles integration with ClinVar database for disease simulation
type ClinVarIntegrator struct {
	baseURL   string
	client    *http.Client
	apiKey    string // For NCBI E-utilities if needed
	rateLimit time.Duration
	cache     map[string]*ClinVarRecord
	rng       *rand.Rand
}

// ClinVarRecord represents a ClinVar variant record
type ClinVarRecord struct {
	VariationID          string                 `json:"variation_id"`
	AlleleID             string                 `json:"allele_id"`
	Type                 string                 `json:"type"`
	Name                 string                 `json:"name"`
	GeneSymbol           string                 `json:"gene_symbol"`
	Chromosome           string                 `json:"chromosome"`
	Position             int64                  `json:"position"`
	ReferenceAllele      string                 `json:"reference_allele"`
	AlternateAllele      string                 `json:"alternate_allele"`
	ClinicalSignificance string                 `json:"clinical_significance"`
	ReviewStatus         string                 `json:"review_status"`
	Conditions           []ClinVarCondition     `json:"conditions"`
	Frequencies          map[string]float64     `json:"frequencies"`
	MolecularConsequence string                 `json:"molecular_consequence"`
	ProteinChange        string                 `json:"protein_change"`
	LastUpdated          time.Time              `json:"last_updated"`
	Citations            []string               `json:"citations"`
	Submissions          []ClinVarSubmission    `json:"submissions"`
}

// ClinVarCondition represents a clinical condition associated with a variant
type ClinVarCondition struct {
	Name          string   `json:"name"`
	ConceptID     string   `json:"concept_id"`
	DatabaseName  string   `json:"database_name"`
	ModeOfInheritance string `json:"mode_of_inheritance"`
	Phenotypes    []string `json:"phenotypes"`
}

// ClinVarSubmission represents submitter information
type ClinVarSubmission struct {
	Submitter             string    `json:"submitter"`
	ClinicalSignificance  string    `json:"clinical_significance"`
	LastEvaluated         time.Time `json:"last_evaluated"`
	ReviewStatus          string    `json:"review_status"`
	SubmissionDate        time.Time `json:"submission_date"`
}

// DiseaseSimulationConfig configures disease-specific simulation parameters
type DiseaseSimulationConfig struct {
	DiseaseName           string                 `json:"disease_name"`
	GeneSymbols          []string               `json:"gene_symbols"`
	PathogenicVariants   []*ClinVarRecord       `json:"pathogenic_variants"`
	CarrierFrequency     float64                `json:"carrier_frequency"`
	PopulationFrequencies map[string]float64    `json:"population_frequencies"`
	PenetranceModel      *PenetranceModel       `json:"penetrance_model"`
	InheritancePattern   string                 `json:"inheritance_pattern"`
	SyntheticCohortSize  int                    `json:"synthetic_cohort_size"`
	AffectedRatio        float64                `json:"affected_ratio"`
}

// PenetranceModel defines age-dependent penetrance curves
type PenetranceModel struct {
	AgeGroups   []AgeGroup `json:"age_groups"`
	Gender      string     `json:"gender"` // "male", "female", "both"
	ModifyingFactors map[string]float64 `json:"modifying_factors"`
}

// AgeGroup defines penetrance for specific age ranges
type AgeGroup struct {
	MinAge     int     `json:"min_age"`
	MaxAge     int     `json:"max_age"`
	Penetrance float64 `json:"penetrance"`
}

// SyntheticDiseaseResult contains generated disease cohort data
type SyntheticDiseaseResult struct {
	DiseaseConfig    *DiseaseSimulationConfig `json:"disease_config"`
	GeneratedSamples []*SyntheticPatient      `json:"generated_samples"`
	Statistics       *CohortStatistics        `json:"statistics"`
	ValidationMetrics *ValidationMetrics      `json:"validation_metrics"`
	GenerationTime   time.Duration            `json:"generation_time"`
}

// SyntheticPatient represents a synthetic patient with disease-relevant variants
type SyntheticPatient struct {
	PatientID           string                    `json:"patient_id"`
	Population          string                    `json:"population"`
	Age                 int                       `json:"age"`
	Gender              string                    `json:"gender"`
	AffectionStatus     string                    `json:"affection_status"` // "affected", "carrier", "unaffected"
	Genotype            map[string]*PatientVariant `json:"genotype"`
	PhenotypeProfile    []string                  `json:"phenotype_profile"`
	FamilyHistory       *FamilyHistory            `json:"family_history"`
	ClinicalAnnotations []string                  `json:"clinical_annotations"`
}

// PatientVariant represents a variant in a synthetic patient
type PatientVariant struct {
	VariantID           string  `json:"variant_id"`
	Chromosome          string  `json:"chromosome"`
	Position            int64   `json:"position"`
	ReferenceAllele     string  `json:"reference_allele"`
	AlternateAllele     string  `json:"alternate_allele"`
	Zygosity            string  `json:"zygosity"` // "homozygous", "heterozygous"
	ClinicalSignificance string `json:"clinical_significance"`
	InheritedFrom       string  `json:"inherited_from"` // "maternal", "paternal", "de_novo"
	PenetranceScore     float64 `json:"penetrance_score"`
}

// FamilyHistory contains relevant family medical history
type FamilyHistory struct {
	ParentalAffection map[string]bool `json:"parental_affection"`
	SiblingAffection  []bool          `json:"sibling_affection"`
	ExtendedFamily    map[string]bool `json:"extended_family"`
}

// CohortStatistics provides statistical summary of generated cohort
type CohortStatistics struct {
	TotalSamples        int                    `json:"total_samples"`
	AffectedCount       int                    `json:"affected_count"`
	CarrierCount        int                    `json:"carrier_count"`
	UnaffectedCount     int                    `json:"unaffected_count"`
	VariantFrequencies  map[string]float64     `json:"variant_frequencies"`
	PopulationBreakdown map[string]int         `json:"population_breakdown"`
	AgeDistribution     map[string]int         `json:"age_distribution"`
	GenderDistribution  map[string]int         `json:"gender_distribution"`
	GenotypeFrequencies map[string]float64     `json:"genotype_frequencies"`
}

// ValidationMetrics compares synthetic cohort against expected disease characteristics
type ValidationMetrics struct {
	CarrierFrequencyAccuracy    float64                `json:"carrier_frequency_accuracy"`
	PopulationFrequencyAccuracy map[string]float64     `json:"population_frequency_accuracy"`
	PenetranceModelAccuracy     float64                `json:"penetrance_model_accuracy"`
	InheritancePatternAccuracy  float64                `json:"inheritance_pattern_accuracy"`
	ClinicalRealism            float64                `json:"clinical_realism"`
	BiologicalConsistency      float64                `json:"biological_consistency"`
	OverallQualityScore        float64                `json:"overall_quality_score"`
}

// Well-characterized diseases for initial implementation
var WellCharacterizedDiseases = map[string]*DiseaseSimulationConfig{
	"cystic_fibrosis": {
		DiseaseName:         "Cystic Fibrosis",
		GeneSymbols:        []string{"CFTR"},
		CarrierFrequency:   0.04,  // ~1 in 25 for European populations
		PopulationFrequencies: map[string]float64{
			"EUR": 0.04,
			"AFR": 0.01,
			"EAS": 0.005,
			"AMR": 0.025,
			"SAS": 0.01,
		},
		InheritancePattern: "autosomal_recessive",
		AffectedRatio:     0.25, // 25% affected in carrier x carrier crosses
	},
	"huntington_disease": {
		DiseaseName:         "Huntington Disease",
		GeneSymbols:        []string{"HTT"},
		CarrierFrequency:   0.0001, // ~1 in 10,000
		PopulationFrequencies: map[string]float64{
			"EUR": 0.0001,
			"AFR": 0.00005,
			"EAS": 0.00002,
			"AMR": 0.00008,
			"SAS": 0.00005,
		},
		InheritancePattern: "autosomal_dominant",
		AffectedRatio:     0.5, // 50% affected in affected x unaffected crosses
		PenetranceModel: &PenetranceModel{
			AgeGroups: []AgeGroup{
				{MinAge: 30, MaxAge: 40, Penetrance: 0.1},
				{MinAge: 40, MaxAge: 50, Penetrance: 0.5},
				{MinAge: 50, MaxAge: 60, Penetrance: 0.8},
				{MinAge: 60, MaxAge: 100, Penetrance: 0.95},
			},
			Gender: "both",
		},
	},
	"brca_hereditary_cancer": {
		DiseaseName:         "BRCA-Associated Hereditary Breast and Ovarian Cancer",
		GeneSymbols:        []string{"BRCA1", "BRCA2"},
		CarrierFrequency:   0.002, // ~1 in 500 for pathogenic variants
		PopulationFrequencies: map[string]float64{
			"EUR": 0.002,
			"AFR": 0.003,  // Higher frequency of VUS, but similar pathogenic
			"EAS": 0.001,
			"AMR": 0.0015,
			"SAS": 0.0018,
			"AJ":  0.025,  // Ashkenazi Jewish founder variants
		},
		InheritancePattern: "autosomal_dominant",
		AffectedRatio:     0.5,
		PenetranceModel: &PenetranceModel{
			AgeGroups: []AgeGroup{
				{MinAge: 30, MaxAge: 40, Penetrance: 0.1},
				{MinAge: 40, MaxAge: 50, Penetrance: 0.3},
				{MinAge: 50, MaxAge: 60, Penetrance: 0.6},
				{MinAge: 60, MaxAge: 70, Penetrance: 0.8},
				{MinAge: 70, MaxAge: 100, Penetrance: 0.85},
			},
			Gender: "female", // Primary penetrance for breast/ovarian cancer
			ModifyingFactors: map[string]float64{
				"family_history": 1.5,
				"reproductive_history": 1.2,
			},
		},
	},
}

// NewClinVarIntegrator creates a new ClinVar integration instance
func NewClinVarIntegrator(apiKey string) *ClinVarIntegrator {
	return &ClinVarIntegrator{
		baseURL:   "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
		client:    &http.Client{Timeout: 30 * time.Second},
		apiKey:    apiKey,
		rateLimit: 100 * time.Millisecond, // NCBI rate limiting
		cache:     make(map[string]*ClinVarRecord),
		rng:       rand.New(rand.NewSource(time.Now().UnixNano())), //nolint:gosec
	}
}

// FetchPathogenicVariants retrieves pathogenic variants for specified genes
func (cv *ClinVarIntegrator) FetchPathogenicVariants(ctx context.Context, geneSymbols []string) ([]*ClinVarRecord, error) {
	var allVariants []*ClinVarRecord
	
	for _, gene := range geneSymbols {
		variants, err := cv.fetchVariantsForGene(ctx, gene)
		if err != nil {
			return nil, fmt.Errorf("failed to fetch variants for gene %s: %w", gene, err)
		}
		allVariants = append(allVariants, variants...)
	}
	
	// Filter for pathogenic/likely pathogenic variants
	pathogenicVariants := cv.filterPathogenicVariants(allVariants)
	
	return pathogenicVariants, nil
}

// fetchVariantsForGene retrieves ClinVar variants for a specific gene
func (cv *ClinVarIntegrator) fetchVariantsForGene(ctx context.Context, geneSymbol string) ([]*ClinVarRecord, error) {
	// Check cache first
	cacheKey := fmt.Sprintf("gene_%s", geneSymbol)
	if cached, exists := cv.cache[cacheKey]; exists {
		return []*ClinVarRecord{cached}, nil
	}
	
	// Rate limiting
	time.Sleep(cv.rateLimit)
	
	// Build search query for ClinVar (URL-encode the full term to produce a valid URL)
	rawQuery := fmt.Sprintf("%s[gene] AND (pathogenic[clin_significance] OR likely pathogenic[clin_significance])", geneSymbol)
	query := url.QueryEscape(rawQuery)

	searchURL := fmt.Sprintf("%sesearch.fcgi?db=clinvar&term=%s&retmode=json", cv.baseURL, query)
	if cv.apiKey != "" {
		searchURL += "&api_key=" + cv.apiKey
	}
	
	req, err := http.NewRequestWithContext(ctx, "GET", searchURL, nil)
	if err != nil {
		return nil, fmt.Errorf("failed to create request: %w", err)
	}
	
	resp, err := cv.client.Do(req) // #nosec G704 -- baseURL is hardcoded to NCBI API endpoint
	if err != nil {
		return nil, fmt.Errorf("failed to fetch ClinVar data: %w", err)
	}
	defer resp.Body.Close()
	
	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("ClinVar API returned status %d", resp.StatusCode)
	}
	
	body, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("failed to read response body: %w", err)
	}
	
	// Parse the response (simplified - actual ClinVar API parsing would be more complex)
	variants, err := cv.parseClinVarResponse(body, geneSymbol)
	if err != nil {
		return nil, fmt.Errorf("failed to parse ClinVar response: %w", err)
	}
	
	return variants, nil
}

// parseClinVarResponse parses ClinVar API response into ClinVarRecord structures
func (cv *ClinVarIntegrator) parseClinVarResponse(responseBody []byte, geneSymbol string) ([]*ClinVarRecord, error) {
	// Simplified parsing - in reality this would handle the complex ClinVar XML/JSON structure
	var variants []*ClinVarRecord
	
	// For demonstration, create synthetic ClinVar records for well-known variants
	switch geneSymbol {
	case "CFTR":
		variants = append(variants, &ClinVarRecord{
			VariationID:          "5",
			AlleleID:             "7",
			Type:                 "single nucleotide variant",
			Name:                 "NM_000492.3(CFTR):c.1521_1523delCTT (p.Phe508del)",
			GeneSymbol:           "CFTR",
			Chromosome:           "7",
			Position:             117199644,
			ReferenceAllele:      "CTT",
			AlternateAllele:      "-",
			ClinicalSignificance: "Pathogenic",
			ReviewStatus:         "criteria provided, multiple submitters, no conflicts",
			Conditions: []ClinVarCondition{
				{
					Name:              "Cystic fibrosis",
					ConceptID:         "C0010674",
					DatabaseName:      "MedGen",
					ModeOfInheritance: "Autosomal recessive inheritance",
				},
			},
			Frequencies: map[string]float64{
				"EUR": 0.024,
				"AFR": 0.005,
			},
			MolecularConsequence: "frameshift variant",
			ProteinChange:        "p.Phe508del",
			LastUpdated:          time.Now(),
		})
		
	case "HTT":
		variants = append(variants, &ClinVarRecord{
			VariationID:          "10",
			AlleleID:             "25",
			Type:                 "microsatellite",
			Name:                 "NM_002111.6(HTT):c.52CAG(40)",
			GeneSymbol:           "HTT",
			Chromosome:           "4",
			Position:             3074877,
			ReferenceAllele:      "(CAG)35",
			AlternateAllele:      "(CAG)40",
			ClinicalSignificance: "Pathogenic",
			ReviewStatus:         "criteria provided, multiple submitters, no conflicts",
			Conditions: []ClinVarCondition{
				{
					Name:              "Huntington disease",
					ConceptID:         "C0020179",
					DatabaseName:      "MedGen",
					ModeOfInheritance: "Autosomal dominant inheritance",
				},
			},
			Frequencies: map[string]float64{
				"EUR": 0.0001,
			},
			MolecularConsequence: "repeat expansion",
			ProteinChange:        "p.(Gln)40",
			LastUpdated:          time.Now(),
		})
		
	case "BRCA1":
		variants = append(variants, &ClinVarRecord{
			VariationID:          "43971",
			AlleleID:             "58076",
			Type:                 "single nucleotide variant",
			Name:                 "NM_007294.3(BRCA1):c.68_69delAG (p.Glu23fs)",
			GeneSymbol:           "BRCA1",
			Chromosome:           "17",
			Position:             43124016,
			ReferenceAllele:      "AG",
			AlternateAllele:      "-",
			ClinicalSignificance: "Pathogenic",
			ReviewStatus:         "criteria provided, multiple submitters, no conflicts",
			Conditions: []ClinVarCondition{
				{
					Name:              "Hereditary breast and ovarian cancer syndrome",
					ConceptID:         "C2676676",
					DatabaseName:      "MedGen",
					ModeOfInheritance: "Autosomal dominant inheritance",
				},
			},
			Frequencies: map[string]float64{
				"EUR": 0.0001,
				"AJ":  0.0125, // Ashkenazi Jewish founder mutation
			},
			MolecularConsequence: "frameshift variant",
			ProteinChange:        "p.Glu23fs",
			LastUpdated:          time.Now(),
		})
	}
	
	return variants, nil
}

// filterPathogenicVariants filters variants to include only pathogenic/likely pathogenic
func (cv *ClinVarIntegrator) filterPathogenicVariants(variants []*ClinVarRecord) []*ClinVarRecord {
	var pathogenic []*ClinVarRecord
	
	pathogenicTerms := []string{
		"Pathogenic",
		"Likely pathogenic",
		"Pathogenic/Likely pathogenic",
	}
	
	for _, variant := range variants {
		for _, term := range pathogenicTerms {
			if strings.Contains(variant.ClinicalSignificance, term) {
				pathogenic = append(pathogenic, variant)
				break
			}
		}
	}
	
	return pathogenic
}

// GenerateDiseaseSimulation creates synthetic disease cohorts
func (cv *ClinVarIntegrator) GenerateDiseaseSimulation(ctx context.Context, config *DiseaseSimulationConfig) (*SyntheticDiseaseResult, error) {
	startTime := time.Now()
	
	// Fetch pathogenic variants if not provided
	if len(config.PathogenicVariants) == 0 {
		variants, err := cv.FetchPathogenicVariants(ctx, config.GeneSymbols)
		if err != nil {
			return nil, fmt.Errorf("failed to fetch pathogenic variants: %w", err)
		}
		config.PathogenicVariants = variants
	}
	
	// Generate synthetic patients
	patients, err := cv.generateSyntheticPatients(config)
	if err != nil {
		return nil, fmt.Errorf("failed to generate synthetic patients: %w", err)
	}
	
	// Calculate cohort statistics
	statistics := cv.calculateCohortStatistics(patients)
	
	// Validate against expected disease characteristics
	validationMetrics := cv.validateCohortRealism(config, patients, statistics)
	
	result := &SyntheticDiseaseResult{
		DiseaseConfig:     config,
		GeneratedSamples:  patients,
		Statistics:        statistics,
		ValidationMetrics: validationMetrics,
		GenerationTime:    time.Since(startTime),
	}
	
	return result, nil
}

// generateSyntheticPatients creates synthetic patients based on disease configuration
func (cv *ClinVarIntegrator) generateSyntheticPatients(config *DiseaseSimulationConfig) ([]*SyntheticPatient, error) {
	var patients []*SyntheticPatient
	
	// Determine number of affected vs unaffected patients
	numAffected := int(float64(config.SyntheticCohortSize) * config.AffectedRatio)
	numUnaffected := config.SyntheticCohortSize - numAffected
	
	// Generate affected patients
	for i := 0; i < numAffected; i++ {
		patient, err := cv.generateAffectedPatient(config, i)
		if err != nil {
			return nil, fmt.Errorf("failed to generate affected patient %d: %w", i, err)
		}
		patients = append(patients, patient)
	}
	
	// Generate unaffected patients (carriers and non-carriers)
	for i := 0; i < numUnaffected; i++ {
		patient, err := cv.generateUnaffectedPatient(config, numAffected+i)
		if err != nil {
			return nil, fmt.Errorf("failed to generate unaffected patient %d: %w", i, err)
		}
		patients = append(patients, patient)
	}
	
	return patients, nil
}

// generateAffectedPatient creates a synthetic affected patient
func (cv *ClinVarIntegrator) generateAffectedPatient(config *DiseaseSimulationConfig, patientIndex int) (*SyntheticPatient, error) {
	patient := &SyntheticPatient{
		PatientID:       fmt.Sprintf("SYN_AFFECTED_%s_%04d", strings.ToUpper(config.DiseaseName), patientIndex),
		AffectionStatus: "affected",
		Genotype:        make(map[string]*PatientVariant),
	}
	
	// Select population (weighted by population frequencies)
	patient.Population = cv.selectPopulation(config.PopulationFrequencies)
	
	// Generate age and gender
	patient.Age = cv.generateAgeForDisease(config)
	patient.Gender = cv.generateGender()
	
	// Generate disease-causing genotype based on inheritance pattern
	err := cv.generateDiseaseGenotype(patient, config, true)
	if err != nil {
		return nil, fmt.Errorf("failed to generate disease genotype: %w", err)
	}
	
	// Generate phenotype profile
	patient.PhenotypeProfile = cv.generatePhenotypeProfile(config, patient)
	
	// Generate family history
	patient.FamilyHistory = cv.generateFamilyHistory(config, true)
	
	return patient, nil
}

// generateUnaffectedPatient creates a synthetic unaffected patient (carrier or non-carrier)
func (cv *ClinVarIntegrator) generateUnaffectedPatient(config *DiseaseSimulationConfig, patientIndex int) (*SyntheticPatient, error) {
	// Determine if patient is a carrier based on carrier frequency
	isCarrier := cv.randomFloat() < config.CarrierFrequency
	
	status := "unaffected"
	if isCarrier {
		status = "carrier"
	}
	
	patient := &SyntheticPatient{
		PatientID:       fmt.Sprintf("SYN_%s_%s_%04d", strings.ToUpper(status), strings.ToUpper(config.DiseaseName), patientIndex),
		AffectionStatus: status,
		Genotype:        make(map[string]*PatientVariant),
	}
	
	// Select population
	patient.Population = cv.selectPopulation(config.PopulationFrequencies)
	
	// Generate age and gender
	patient.Age = cv.generateAgeForGeneral()
	patient.Gender = cv.generateGender()
	
	// Generate genotype
	err := cv.generateDiseaseGenotype(patient, config, false)
	if err != nil {
		return nil, fmt.Errorf("failed to generate genotype: %w", err)
	}
	
	// Generate family history
	patient.FamilyHistory = cv.generateFamilyHistory(config, false)
	
	return patient, nil
}

// Helper methods for patient generation (simplified implementations)

func (cv *ClinVarIntegrator) selectPopulation(frequencies map[string]float64) string {
	// Weighted random selection based on population frequencies
	r := cv.randomFloat()
	cumulative := 0.0
	
	for pop, freq := range frequencies {
		cumulative += freq
		if r <= cumulative {
			return pop
		}
	}
	
	return "EUR" // Default fallback
}

func (cv *ClinVarIntegrator) generateAgeForDisease(config *DiseaseSimulationConfig) int {
	if config.PenetranceModel == nil || len(config.PenetranceModel.AgeGroups) == 0 {
		return cv.randomInt(18, 80)
	}

	// Compute incremental penetrance weights for each age group.
	// AgeGroups is assumed to be ordered with monotonically increasing cumulative
	// penetrance (e.g., 0.1 by age 40, 0.5 by age 50, 0.8 by age 60).
	// The incremental weight for group i is Penetrance[i] - Penetrance[i-1].
	groups := config.PenetranceModel.AgeGroups
	weights := make([]float64, len(groups))
	prev := 0.0
	total := 0.0
	for i, g := range groups {
		w := g.Penetrance - prev
		if w < 0 {
			w = 0 // guard against non-monotone data
		}
		weights[i] = w
		total += w
		prev = g.Penetrance
	}

	// Fall back to uniform if all weights are zero.
	if total <= 0 {
		return cv.randomInt(18, 80)
	}

	// Weighted random selection of an age group.
	r := cv.randomFloat() * total
	cumulative := 0.0
	for i, g := range groups {
		cumulative += weights[i]
		if r <= cumulative {
			return cv.randomInt(g.MinAge, g.MaxAge)
		}
	}
	// Fallback to last group if floating-point rounding slips past.
	g := groups[len(groups)-1]
	return cv.randomInt(g.MinAge, g.MaxAge)
}

func (cv *ClinVarIntegrator) generateAgeForGeneral() int {
	return cv.randomInt(18, 80)
}

func (cv *ClinVarIntegrator) generateGender() string {
	if cv.randomFloat() < 0.5 {
		return "male"
	}
	return "female"
}

func (cv *ClinVarIntegrator) generateDiseaseGenotype(patient *SyntheticPatient, config *DiseaseSimulationConfig, isAffected bool) error {
	// Implementation would be disease-specific based on inheritance pattern
	// This is a simplified version
	
	if len(config.PathogenicVariants) == 0 {
		return fmt.Errorf("no pathogenic variants available for genotype generation")
	}

	// Select a variant weighted by its allele frequency in the patient's population.
	// Falls back to equal weights when no frequency data is available.
	variant := cv.selectVariantByFrequency(config.PathogenicVariants, patient.Population)
	
	patientVariant := &PatientVariant{
		VariantID:           variant.VariationID,
		Chromosome:          variant.Chromosome,
		Position:            variant.Position,
		ReferenceAllele:     variant.ReferenceAllele,
		AlternateAllele:     variant.AlternateAllele,
		ClinicalSignificance: variant.ClinicalSignificance,
	}
	
	// Set zygosity based on inheritance pattern and affection status
	switch config.InheritancePattern {
	case "autosomal_recessive":
		if isAffected {
			patientVariant.Zygosity = "homozygous"
		} else if patient.AffectionStatus == "carrier" {
			patientVariant.Zygosity = "heterozygous"
		}
	case "autosomal_dominant":
		if isAffected || patient.AffectionStatus == "carrier" {
			patientVariant.Zygosity = "heterozygous"
		}
	}
	
	patient.Genotype[variant.GeneSymbol] = patientVariant
	
	return nil
}

// selectVariantByFrequency picks a variant from the slice using population-specific
// allele frequencies as weights. If a variant has no frequency recorded for the
// given population it receives weight 1.0 (equal chance). The total weight is
// normalised before the draw so that the result is a proper probability distribution.
func (cv *ClinVarIntegrator) selectVariantByFrequency(variants []*ClinVarRecord, population string) *ClinVarRecord {
	if len(variants) == 1 {
		return variants[0]
	}

	weights := make([]float64, len(variants))
	total := 0.0
	for i, v := range variants {
		w := 1.0 // default weight when no frequency data
		if v.Frequencies != nil {
			if freq, ok := v.Frequencies[population]; ok && freq > 0 {
				w = freq
			}
		}
		weights[i] = w
		total += w
	}

	r := cv.randomFloat() * total
	cumulative := 0.0
	for i, w := range weights {
		cumulative += w
		if r <= cumulative {
			return variants[i]
		}
	}
	return variants[len(variants)-1]
}

func (cv *ClinVarIntegrator) generatePhenotypeProfile(config *DiseaseSimulationConfig, patient *SyntheticPatient) []string {
	// Generate realistic phenotype based on disease and genotype
	var phenotypes []string
	
	for _, variant := range config.PathogenicVariants {
		for _, condition := range variant.Conditions {
			phenotypes = append(phenotypes, condition.Phenotypes...)
		}
	}
	
	return phenotypes
}

func (cv *ClinVarIntegrator) generateFamilyHistory(config *DiseaseSimulationConfig, hasDisease bool) *FamilyHistory {
	history := &FamilyHistory{
		ParentalAffection: make(map[string]bool),
		ExtendedFamily:    make(map[string]bool),
	}
	
	// Generate family history based on inheritance pattern
	if hasDisease && config.InheritancePattern == "autosomal_dominant" {
		// Higher probability of affected parent
		history.ParentalAffection["mother"] = cv.randomFloat() < 0.5
		history.ParentalAffection["father"] = cv.randomFloat() < 0.5
	}
	
	return history
}

func (cv *ClinVarIntegrator) calculateCohortStatistics(patients []*SyntheticPatient) *CohortStatistics {
	stats := &CohortStatistics{
		TotalSamples:        len(patients),
		VariantFrequencies:  make(map[string]float64),
		PopulationBreakdown: make(map[string]int),
		AgeDistribution:     make(map[string]int),
		GenderDistribution:  make(map[string]int),
		GenotypeFrequencies: make(map[string]float64),
	}
	
	for _, patient := range patients {
		// Count affection status
		switch patient.AffectionStatus {
		case "affected":
			stats.AffectedCount++
		case "carrier":
			stats.CarrierCount++
		case "unaffected":
			stats.UnaffectedCount++
		}
		
		// Count population breakdown
		stats.PopulationBreakdown[patient.Population]++
		
		// Count gender distribution
		stats.GenderDistribution[patient.Gender]++
		
		// Age distribution (simplified)
		ageGroup := fmt.Sprintf("%d-%d", (patient.Age/10)*10, ((patient.Age/10)*10)+9)
		stats.AgeDistribution[ageGroup]++
	}
	
	return stats
}

func (cv *ClinVarIntegrator) validateCohortRealism(config *DiseaseSimulationConfig, patients []*SyntheticPatient, stats *CohortStatistics) *ValidationMetrics {
	metrics := &ValidationMetrics{
		PopulationFrequencyAccuracy: make(map[string]float64),
	}
	if len(patients) == 0 {
		return metrics // cannot compute frequency ratios with no patients (#113)
	}

	// Validate carrier frequency
	expectedCarrierCount := int(float64(len(patients)) * config.CarrierFrequency)
	actualCarrierCount := stats.CarrierCount
	
	if expectedCarrierCount > 0 {
		metrics.CarrierFrequencyAccuracy = 1.0 - float64(abs(expectedCarrierCount-actualCarrierCount))/float64(expectedCarrierCount)
	}
	
	// Validate population frequencies
	for pop, expectedFreq := range config.PopulationFrequencies {
		actualCount := stats.PopulationBreakdown[pop]
		actualFreq := float64(actualCount) / float64(len(patients))
		
		if expectedFreq > 0 {
			metrics.PopulationFrequencyAccuracy[pop] = 1.0 - float64(abs(int(expectedFreq*100)-int(actualFreq*100)))/100.0
		}
	}
	
	// Calculate overall quality score
	scores := []float64{
		metrics.CarrierFrequencyAccuracy,
	}
	
	for _, score := range metrics.PopulationFrequencyAccuracy {
		scores = append(scores, score)
	}
	
	if len(scores) > 0 {
		sum := 0.0
		for _, score := range scores {
			sum += score
		}
		metrics.OverallQualityScore = sum / float64(len(scores))
	}
	
	return metrics
}

// Utility functions
func (cv *ClinVarIntegrator) randomFloat() float64 {
	return cv.rng.Float64()
}

func (cv *ClinVarIntegrator) randomInt(min, max int) int {
	if max <= min {
		return min
	}
	return min + cv.rng.Intn(max-min+1)
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}