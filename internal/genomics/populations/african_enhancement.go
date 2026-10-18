// Package populations builds toy genotype panels for a fixed registry of
// African populations (five 1000 Genomes populations and AWI-Gen), each with
// allele frequencies for a handful of rsIDs (APOE, ACTN3, LCT, CYP2C19).
//
// It does not query any literature or H3Africa resource: the "literature
// search" returns synthetic placeholder entries, and the model's allele
// frequency accuracy compares the registry with a copy of itself. Samples
// carry genotypes for the registry rsIDs only. No generation or CLI code path
// uses this package.
package populations

import (
	"context"
	"fmt"
	"log"
	"math"
	"math/rand/v2"
	"net/http"
	"time"
)

// africanFST contains approximate pairwise F_ST values between each African
// population and the four non-African 1000 Genomes super-populations. They
// are attributed to Auton et al. (2015) Nature 526:68-74, but the numbers have
// not been checked against that paper, and the AWI row (AWI-Gen is an
// H3Africa cohort, not a 1000 Genomes population) cannot come from it.
// Populations not listed fall back to YRI (see africanFSTForPop).
var africanFST = map[string]map[string]float64{
	// Continental African populations
	"YRI": {"EUR": 0.153, "EAS": 0.190, "AMR": 0.123, "SAS": 0.160},
	"LWK": {"EUR": 0.143, "EAS": 0.180, "AMR": 0.113, "SAS": 0.150},
	"GWD": {"EUR": 0.148, "EAS": 0.185, "AMR": 0.118, "SAS": 0.155},
	"MSL": {"EUR": 0.150, "EAS": 0.188, "AMR": 0.120, "SAS": 0.158},
	"ESN": {"EUR": 0.151, "EAS": 0.189, "AMR": 0.121, "SAS": 0.159},
	// Diaspora / admixed populations have lower FST with EUR due to European admixture
	"AWI": {"EUR": 0.107, "EAS": 0.152, "AMR": 0.077, "SAS": 0.120},
	"ASW": {"EUR": 0.107, "EAS": 0.152, "AMR": 0.077, "SAS": 0.120},
	"ACB": {"EUR": 0.115, "EAS": 0.158, "AMR": 0.085, "SAS": 0.128},
}

// africanFSTForPop returns the FST map for a given African population ID,
// falling back to YRI (the largest and most studied sub-population) when the
// ID is not found in the lookup table.
func africanFSTForPop(populationID string) map[string]float64 {
	if fst, ok := africanFST[populationID]; ok {
		return fst
	}
	return africanFST["YRI"]
}

// AfricanPopulationEnhancer builds AfricanPopulationModels from
// H3AfricaPopulationRegistry and synthetic literature entries, and generates
// toy samples from them. Its random source is seeded from the clock.
type AfricanPopulationEnhancer struct {
	literatureDB     *LiteratureDatabase
	h3africaDB       *H3AfricaDatabase
	populationModels map[string]*AfricanPopulationModel
	client           *http.Client
	cache            map[string]*LiteratureEntry
	rng              *rand.Rand
}

// H3AfricaDatabase holds the built-in population registry. Despite the name
// it is not loaded from H3Africa; StudyRegistry and VariantCatalog are empty.
type H3AfricaDatabase struct {
	PopulationPanels map[string]*H3AfricaPopulation `json:"population_panels"`
	StudyRegistry    []*H3AfricaStudy               `json:"study_registry"`
	VariantCatalog   map[string]*AfricanVariant     `json:"variant_catalog"`
	LastUpdated      time.Time                      `json:"last_updated"`
}

// H3AfricaPopulation describes one population in the built-in registry.
type H3AfricaPopulation struct {
	PopulationID            string               `json:"population_id"`
	PopulationName          string               `json:"population_name"`
	Country                 string               `json:"country"`
	Region                  string               `json:"region"`
	LanguageFamily          string               `json:"language_family"`
	EthnicGroup             string               `json:"ethnic_group"`
	SampleSize              int                  `json:"sample_size"`
	StudyPeriod             []int                `json:"study_period"`
	AlleleFrequencies       map[string]float64   `json:"allele_frequencies"`
	LinkagePatterns         map[string]*LDBlock  `json:"linkage_patterns"`
	DemographicParams       *DemographicHistory  `json:"demographic_params"`
	PharmacogenomicVariants []*PGxVariant        `json:"pharmacogenomic_variants"`
	StructuralVariants      []*StructuralVariant `json:"structural_variants"`
	GeographicCoords        []float64            `json:"geographic_coords"`
	Citations               []string             `json:"citations"`
}

// H3AfricaStudy represents a study from the H3Africa consortium
type H3AfricaStudy struct {
	StudyID               string    `json:"study_id"`
	StudyName             string    `json:"study_name"`
	PrincipalInvestigator string    `json:"principal_investigator"`
	Institution           string    `json:"institution"`
	Country               string    `json:"country"`
	StudyType             string    `json:"study_type"`
	Populations           []string  `json:"populations"`
	SampleSize            int       `json:"sample_size"`
	SequencingPlatform    string    `json:"sequencing_platform"`
	DataType              string    `json:"data_type"`
	PublicationDate       time.Time `json:"publication_date"`
	DOI                   string    `json:"doi"`
	DataAvailability      string    `json:"data_availability"`
	EthicsApproval        string    `json:"ethics_approval"`
}

// AfricanVariant represents a population-specific variant
type AfricanVariant struct {
	VariantID         string                 `json:"variant_id"`
	Chromosome        string                 `json:"chromosome"`
	Position          int64                  `json:"position"`
	ReferenceAllele   string                 `json:"reference_allele"`
	AlternateAllele   string                 `json:"alternate_allele"`
	PopulationFreqs   map[string]float64     `json:"population_frequencies"`
	FunctionalClass   string                 `json:"functional_class"`
	GeneSymbol        string                 `json:"gene_symbol"`
	Consequence       string                 `json:"consequence"`
	ClinicalRelevance string                 `json:"clinical_relevance"`
	DiscoveryStudy    string                 `json:"discovery_study"`
	ValidationStatus  string                 `json:"validation_status"`
	Annotations       map[string]interface{} `json:"annotations"`
}

// LDBlock represents linkage disequilibrium patterns
type LDBlock struct {
	BlockID           string             `json:"block_id"`
	Chromosome        string             `json:"chromosome"`
	StartPos          int64              `json:"start_position"`
	EndPos            int64              `json:"end_position"`
	Variants          []string           `json:"variants"`
	HaplotypeFreqs    map[string]float64 `json:"haplotype_frequencies"`
	RecombinationRate float64            `json:"recombination_rate"`
	Population        string             `json:"population"`
}

// DemographicHistory represents population demographic parameters
type DemographicHistory struct {
	EffectivePopulationSize []PopulationEpoch  `json:"effective_population_size"`
	MigrationEvents         []*MigrationEvent  `json:"migration_events"`
	BottleneckEvents        []*BottleneckEvent `json:"bottleneck_events"`
	PopulationSplits        []*PopulationSplit `json:"population_splits"`
	CoalescentParameters    *CoalescentParams  `json:"coalescent_parameters"`
}

// PopulationEpoch represents population size at different time periods
type PopulationEpoch struct {
	TimeAgo            int     `json:"time_ago_generations"`
	PopulationSize     int     `json:"population_size"`
	GrowthRate         float64 `json:"growth_rate"`
	ConfidenceInterval []int   `json:"confidence_interval"`
}

// MigrationEvent represents population migration
type MigrationEvent struct {
	SourcePopulation string  `json:"source_population"`
	TargetPopulation string  `json:"target_population"`
	MigrationRate    float64 `json:"migration_rate"`
	StartTime        int     `json:"start_time_generations"`
	EndTime          int     `json:"end_time_generations"`
	MigrationPulse   bool    `json:"migration_pulse"`
}

// BottleneckEvent represents population bottlenecks
type BottleneckEvent struct {
	PopulationName string  `json:"population_name"`
	BottleneckTime int     `json:"bottleneck_time_generations"`
	OriginalSize   int     `json:"original_size"`
	BottleneckSize int     `json:"bottleneck_size"`
	RecoveryTime   int     `json:"recovery_time"`
	Severity       float64 `json:"severity"`
}

// PopulationSplit represents population divergence events
type PopulationSplit struct {
	AncestralPopulation string   `json:"ancestral_population"`
	DerivedPopulations  []string `json:"derived_populations"`
	SplitTime           int      `json:"split_time_generations"`
	DivergenceIndex     float64  `json:"divergence_index"`
}

// CoalescentParams represents coalescent model parameters
type CoalescentParams struct {
	ThetaWatterson      float64 `json:"theta_watterson"`
	ThetaPi             float64 `json:"theta_pi"`
	TajimasD            float64 `json:"tajimas_d"`
	FuAndLiD            float64 `json:"fu_and_li_d"`
	NucleotideDiversity float64 `json:"nucleotide_diversity"`
	SegregationSites    int     `json:"segregation_sites"`
}

// PGxVariant represents pharmacogenomic variants
type PGxVariant struct {
	VariantID      string             `json:"variant_id"`
	GeneSymbol     string             `json:"gene_symbol"`
	DrugName       string             `json:"drug_name"`
	Phenotype      string             `json:"phenotype"`
	PopulationFreq map[string]float64 `json:"population_frequency"`
	ClinicalImpact string             `json:"clinical_impact"`
	StudyEvidence  []string           `json:"study_evidence"`
}

// StructuralVariant represents large structural variants
type StructuralVariant struct {
	SVID             string             `json:"sv_id"`
	Type             string             `json:"type"` // DEL, DUP, INV, TRA
	Chromosome       string             `json:"chromosome"`
	StartPos         int64              `json:"start_position"`
	EndPos           int64              `json:"end_position"`
	Size             int                `json:"size"`
	PopulationFreq   map[string]float64 `json:"population_frequency"`
	GeneOverlap      []string           `json:"gene_overlap"`
	FunctionalImpact string             `json:"functional_impact"`
}

// AfricanPopulationModel contains enhanced modeling parameters for African populations
type AfricanPopulationModel struct {
	Population          *H3AfricaPopulation   `json:"population"`
	SyntheticParameters *SyntheticGenParams   `json:"synthetic_parameters"`
	ValidationMetrics   *PopValidationMetrics `json:"validation_metrics"`
	GenerationHistory   []*GenerationEvent    `json:"generation_history"`
	QualityAssessment   *QualityMetrics       `json:"quality_assessment"`
}

// SyntheticGenParams contains parameters for synthetic data generation
type SyntheticGenParams struct {
	VariantDensity      float64              `json:"variant_density"`
	TransitionTsv       float64              `json:"transition_transversion_ratio"`
	IndelRatio          float64              `json:"indel_ratio"`
	StructuralVarRate   float64              `json:"structural_variant_rate"`
	RecombinationMap    map[string][]float64 `json:"recombination_map"`
	MutationRates       map[string]float64   `json:"mutation_rates"`
	SelectionCoeff      map[string]float64   `json:"selection_coefficients"`
	PopulationStructure *PopulationStructure `json:"population_structure"`
	// AlleleFreqs holds the allele frequencies (keyed by variant ID) that the
	// model encodes as its generation target, copied from the source population
	// data. Used in validation to compare model parameters against the reference.
	AlleleFreqs map[string]float64 `json:"allele_frequencies"`
}

// PopulationStructure defines population genetic structure
type PopulationStructure struct {
	Fst                  map[string]float64 `json:"fst_values"`
	AdmixtureProportions map[string]float64 `json:"admixture_proportions"`
	IsolationByDistance  bool               `json:"isolation_by_distance"`
	SubpopulationSizes   map[string]int     `json:"subpopulation_sizes"`
	MigrationMatrix      [][]float64        `json:"migration_matrix"`
}

// PopValidationMetrics contains validation metrics for population models
type PopValidationMetrics struct {
	AlleleFreqAccuracy     float64            `json:"allele_freq_accuracy"`
	LinkageAccuracy        float64            `json:"linkage_accuracy"`
	DemographicFit         float64            `json:"demographic_fit"`
	StructuralVarAccuracy  float64            `json:"structural_variant_accuracy"`
	PharmacogenomicFit     float64            `json:"pharmacogenomic_fit"`
	OverallValidationScore float64            `json:"overall_validation_score"`
	ComparisonMetrics      map[string]float64 `json:"comparison_metrics"`
}

// GenerationEvent tracks synthetic data generation events
type GenerationEvent struct {
	Timestamp       time.Time              `json:"timestamp"`
	EventType       string                 `json:"event_type"`
	Parameters      map[string]interface{} `json:"parameters"`
	ResultMetrics   map[string]float64     `json:"result_metrics"`
	DataGenerated   int                    `json:"data_generated"`
	ValidationScore float64                `json:"validation_score"`
}

// QualityMetrics provides quality assessment of synthetic data
type QualityMetrics struct {
	BiologicalRealism     float64 `json:"biological_realism"`
	StatisticalFidelity   float64 `json:"statistical_fidelity"`
	PopulationSpecificity float64 `json:"population_specificity"`
	ClinicalRelevance     float64 `json:"clinical_relevance"`
	OverallQuality        float64 `json:"overall_quality"`
}

// LiteratureDatabase manages genomics literature integration
type LiteratureDatabase struct {
	Entries        map[string]*LiteratureEntry `json:"entries"`
	SearchIndex    map[string][]string         `json:"search_index"`
	PopulationTags map[string][]string         `json:"population_tags"`
	LastUpdated    time.Time                   `json:"last_updated"`
}

// LiteratureEntry represents a genomics literature entry
type LiteratureEntry struct {
	PubMedID           string                 `json:"pubmed_id"`
	DOI                string                 `json:"doi"`
	Title              string                 `json:"title"`
	Authors            []string               `json:"authors"`
	Journal            string                 `json:"journal"`
	PublicationDate    time.Time              `json:"publication_date"`
	Abstract           string                 `json:"abstract"`
	PopulationsFocused []string               `json:"populations_focused"`
	StudyType          string                 `json:"study_type"`
	SampleSize         int                    `json:"sample_size"`
	ExtractedData      *ExtractedGenomicsData `json:"extracted_data"`
	QualityScore       float64                `json:"quality_score"`
	CurationStatus     string                 `json:"curation_status"`
}

// ExtractedGenomicsData contains automatically extracted genomics parameters
type ExtractedGenomicsData struct {
	AlleleFrequencies     map[string]map[string]float64 `json:"allele_frequencies"`
	LinkageData           map[string]*LDBlock           `json:"linkage_data"`
	DemographicParams     *DemographicHistory           `json:"demographic_params"`
	PopulationComparisons map[string]float64            `json:"population_comparisons"`
	StatisticalSummaries  map[string]float64            `json:"statistical_summaries"`
	ClinicalAssociations  []string                      `json:"clinical_associations"`
	Confidence            float64                       `json:"confidence"`
}

// H3AfricaPopulationRegistry is the built-in population table: five 1000
// Genomes African populations (YRI, LWK, GWD, MSL, ESN) and the AWI-Gen
// cohort. Sample sizes and coordinates are descriptive; the allele
// frequencies are approximate and have not been checked against the cited
// sources.
var H3AfricaPopulationRegistry = map[string]*H3AfricaPopulation{
	"YRI": {
		PopulationID:   "YRI",
		PopulationName: "Yoruba in Ibadan, Nigeria",
		Country:        "Nigeria",
		Region:         "West Africa",
		LanguageFamily: "Niger-Congo",
		EthnicGroup:    "Yoruba",
		SampleSize:     108,
		StudyPeriod:    []int{2002, 2005},
		AlleleFrequencies: map[string]float64{
			"rs429358":  0.205, // APOE4 variant
			"rs7412":    0.123, // APOE2 variant
			"rs1815739": 0.000, // ACTN3 R577X (practically absent)
			"rs4988235": 0.001, // LCT lactase persistence (very low)
		},
		GeographicCoords: []float64{7.3775, 3.9470}, // Ibadan coordinates
		Citations: []string{
			"10.1038/nature15393",       // 1000 Genomes
			"10.1038/s41586-020-2859-7", // H3Africa AWI-Gen
		},
	},
	"LWK": {
		PopulationID:   "LWK",
		PopulationName: "Luhya in Webuye, Kenya",
		Country:        "Kenya",
		Region:         "East Africa",
		LanguageFamily: "Niger-Congo",
		EthnicGroup:    "Luhya",
		SampleSize:     99,
		StudyPeriod:    []int{2002, 2005},
		AlleleFrequencies: map[string]float64{
			"rs429358":  0.176,
			"rs7412":    0.141,
			"rs1815739": 0.010,
			"rs4988235": 0.023,
		},
		GeographicCoords: []float64{0.6167, 34.7667}, // Webuye coordinates
		Citations: []string{
			"10.1038/nature15393",
			"10.1186/s13073-018-0555-x", // H3Africa kidney disease study
		},
	},
	"GWD": {
		PopulationID:   "GWD",
		PopulationName: "Gambian in Western Divisions in the Gambia",
		Country:        "Gambia",
		Region:         "West Africa",
		LanguageFamily: "Niger-Congo",
		EthnicGroup:    "Gambian",
		SampleSize:     113,
		StudyPeriod:    []int{2010, 2012},
		AlleleFrequencies: map[string]float64{
			"rs429358":  0.193,
			"rs7412":    0.088,
			"rs1815739": 0.000,
			"rs4988235": 0.002,
		},
		GeographicCoords: []float64{13.4432, -15.3101}, // Banjul coordinates
		PharmacogenomicVariants: []*PGxVariant{
			{
				VariantID:      "rs4244285",
				GeneSymbol:     "CYP2C19",
				DrugName:       "Clopidogrel",
				Phenotype:      "Poor Metabolizer",
				PopulationFreq: map[string]float64{"GWD": 0.155},
				ClinicalImpact: "Reduced efficacy",
			},
		},
	},
	"MSL": {
		PopulationID:   "MSL",
		PopulationName: "Mende in Sierra Leone",
		Country:        "Sierra Leone",
		Region:         "West Africa",
		LanguageFamily: "Niger-Congo",
		EthnicGroup:    "Mende",
		SampleSize:     85,
		StudyPeriod:    []int{2010, 2012},
		AlleleFrequencies: map[string]float64{
			"rs429358":  0.221,
			"rs7412":    0.106,
			"rs1815739": 0.000,
			"rs4988235": 0.000,
		},
		GeographicCoords: []float64{7.9465, -11.1986}, // Bo coordinates
	},
	"ESN": {
		PopulationID:   "ESN",
		PopulationName: "Esan in Nigeria",
		Country:        "Nigeria",
		Region:         "West Africa",
		LanguageFamily: "Niger-Congo",
		EthnicGroup:    "Esan",
		SampleSize:     99,
		StudyPeriod:    []int{2010, 2012},
		AlleleFrequencies: map[string]float64{
			"rs429358":  0.212,
			"rs7412":    0.131,
			"rs1815739": 0.005,
			"rs4988235": 0.000,
		},
		GeographicCoords: []float64{6.7591, 6.1745}, // Edo State coordinates
	},
	"AWI": {
		PopulationID:   "AWI",
		PopulationName: "Africa Wits-INDEPTH Partnership for Genomic Studies",
		Country:        "Multi-country",
		Region:         "Sub-Saharan Africa",
		LanguageFamily: "Mixed",
		EthnicGroup:    "Multi-ethnic",
		SampleSize:     10603,
		StudyPeriod:    []int{2013, 2016},
		AlleleFrequencies: map[string]float64{
			"rs429358":  0.198,
			"rs7412":    0.119,
			"rs1815739": 0.003,
			"rs4988235": 0.012,
		},
		Citations: []string{
			"10.1038/s41586-020-2859-7",  // Main AWI-Gen publication
			"10.1038/s41467-021-22207-5", // AWI-Gen genomic analysis
		},
	},
}

// NewAfricanPopulationEnhancer creates a new African population enhancement system
func NewAfricanPopulationEnhancer() *AfricanPopulationEnhancer {
	enhancer := &AfricanPopulationEnhancer{
		populationModels: make(map[string]*AfricanPopulationModel),
		client:           &http.Client{Timeout: 30 * time.Second},
		cache:            make(map[string]*LiteratureEntry),
		rng:              rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)), // #nosec G404
	}

	// Initialize H3Africa database
	enhancer.h3africaDB = &H3AfricaDatabase{
		PopulationPanels: H3AfricaPopulationRegistry,
		LastUpdated:      time.Now(),
	}

	// Initialize literature database
	enhancer.literatureDB = &LiteratureDatabase{
		Entries:        make(map[string]*LiteratureEntry),
		SearchIndex:    make(map[string][]string),
		PopulationTags: make(map[string][]string),
		LastUpdated:    time.Now(),
	}

	return enhancer
}

// EnhanceAfricanPopulations builds a model for each listed registry
// population from the registry entry and synthetic literature entries (no
// literature is retrieved). It fails if a population is not in the registry.
func (ape *AfricanPopulationEnhancer) EnhanceAfricanPopulations(ctx context.Context, targetPopulations []string) (map[string]*AfricanPopulationModel, error) {
	enhancedModels := make(map[string]*AfricanPopulationModel)

	for _, popID := range targetPopulations {
		if basePop, exists := ape.h3africaDB.PopulationPanels[popID]; exists {
			// Create enhanced model for this population
			model, err := ape.createEnhancedModel(ctx, basePop)
			if err != nil {
				return nil, fmt.Errorf("failed to enhance model for population %s: %w", popID, err)
			}

			enhancedModels[popID] = model
			ape.populationModels[popID] = model
		} else {
			return nil, fmt.Errorf("population %s not found in H3Africa registry", popID)
		}
	}

	return enhancedModels, nil
}

// createEnhancedModel builds a population model from the registry entry and
// the synthetic literature entries.
func (ape *AfricanPopulationEnhancer) createEnhancedModel(ctx context.Context, basePop *H3AfricaPopulation) (*AfricanPopulationModel, error) {
	// Search for relevant literature
	literatureEntries, err := ape.searchRelevantLiterature(ctx, basePop.PopulationID, basePop.Country)
	if err != nil {
		return nil, fmt.Errorf("failed to search literature: %w", err)
	}

	// Extract genomics parameters from literature
	extractedParams, err := ape.extractGenomicsParameters(literatureEntries, basePop)
	if err != nil {
		return nil, fmt.Errorf("failed to extract parameters: %w", err)
	}

	// Create synthetic generation parameters
	syntheticParams := ape.createSyntheticParameters(basePop, extractedParams)

	// Validate enhanced model
	validationMetrics := ape.validateEnhancedModel(basePop, syntheticParams)

	// Create enhanced model
	model := &AfricanPopulationModel{
		Population:          basePop,
		SyntheticParameters: syntheticParams,
		ValidationMetrics:   validationMetrics,
		GenerationHistory:   []*GenerationEvent{},
		QualityAssessment:   ape.assessModelQuality(syntheticParams, validationMetrics),
	}

	return model, nil
}

// searchRelevantLiterature is a stub: it makes no network request and
// returns one synthetic entry per search term.
func (ape *AfricanPopulationEnhancer) searchRelevantLiterature(ctx context.Context, populationID, country string) ([]*LiteratureEntry, error) {
	searchTerms := []string{
		fmt.Sprintf("%s population genetics", populationID),
		fmt.Sprintf("%s genomics", country),
		fmt.Sprintf("African %s allele frequency", populationID),
		"H3Africa consortium",
		fmt.Sprintf("%s pharmacogenomics", country),
	}

	var entries []*LiteratureEntry

	// Synthetic entries with fixed values (not taken from any study)
	for _, term := range searchTerms {
		entry := ape.createSyntheticLiteratureEntry(term, populationID, country)
		entries = append(entries, entry)
	}

	return entries, nil
}

// createSyntheticLiteratureEntry creates a placeholder literature entry. The
// IDs, title, authors and journal are synthetic, and the extracted values are
// the same fixed numbers for every population apart from F_ST.
func (ape *AfricanPopulationEnhancer) createSyntheticLiteratureEntry(searchTerm, populationID, country string) *LiteratureEntry {
	entry := &LiteratureEntry{
		PubMedID:           fmt.Sprintf("PMID_%s_%d", populationID, time.Now().Unix()%10000),
		DOI:                fmt.Sprintf("synthetic/example.%s.%d", populationID, time.Now().Unix()%1000),
		Title:              fmt.Sprintf("Population genetics analysis of %s populations in %s", populationID, country),
		Authors:            []string{"Author A", "Author B"},
		Journal:            "synthetic placeholder",
		PublicationDate:    time.Now().AddDate(-1, 0, 0),
		PopulationsFocused: []string{populationID},
		StudyType:          "Population genetics",
		SampleSize:         1000,
		QualityScore:       0.85,
		CurationStatus:     "synthetic",
	}

	// Extracted data: F_ST from africanFST, other values fixed
	fst := africanFSTForPop(populationID)
	entry.ExtractedData = &ExtractedGenomicsData{
		AlleleFrequencies: make(map[string]map[string]float64),
		LinkageData:       make(map[string]*LDBlock),
		PopulationComparisons: map[string]float64{
			"fst_eur": fst["EUR"],
			"fst_eas": fst["EAS"],
		},
		StatisticalSummaries: map[string]float64{
			"nucleotide_diversity": 0.00125,
			"tajimas_d":            0.32,
			"theta_watterson":      0.00108,
		},
		Confidence: 0.90,
	}

	// Fixed allele frequencies, the same for every population
	entry.ExtractedData.AlleleFrequencies[populationID] = map[string]float64{
		"rs429358":  0.20,  // APOE4
		"rs7412":    0.12,  // APOE2
		"rs1815739": 0.001, // ACTN3
		"rs4988235": 0.005, // LCT
		"rs4244285": 0.15,  // CYP2C19
		"rs776746":  0.60,  // CYP3A5
	}

	return entry
}

// extractGenomicsParameters combines the entries' allele frequencies and
// summaries. Repeated values are merged as (existing + value·weight)/2 and
// then divided by the total weight, which is not a weighted mean. The result
// is passed to createSyntheticParameters, which ignores it.
func (ape *AfricanPopulationEnhancer) extractGenomicsParameters(entries []*LiteratureEntry, basePop *H3AfricaPopulation) (*ExtractedGenomicsData, error) {
	combined := &ExtractedGenomicsData{
		AlleleFrequencies:     make(map[string]map[string]float64),
		LinkageData:           make(map[string]*LDBlock),
		PopulationComparisons: make(map[string]float64),
		StatisticalSummaries:  make(map[string]float64),
		ClinicalAssociations:  []string{},
	}

	// Aggregate data from multiple literature sources
	totalWeight := 0.0

	for _, entry := range entries {
		weight := entry.QualityScore * float64(entry.SampleSize) / 1000.0
		totalWeight += weight

		if entry.ExtractedData != nil {
			// Weighted combination of allele frequencies
			for pop, variants := range entry.ExtractedData.AlleleFrequencies {
				if pop == basePop.PopulationID {
					if combined.AlleleFrequencies[pop] == nil {
						combined.AlleleFrequencies[pop] = make(map[string]float64)
					}

					for variant, freq := range variants {
						if existing, exists := combined.AlleleFrequencies[pop][variant]; exists {
							combined.AlleleFrequencies[pop][variant] = (existing + freq*weight) / 2.0
						} else {
							combined.AlleleFrequencies[pop][variant] = freq * weight
						}
					}
				}
			}

			// Combine statistical summaries
			for metric, value := range entry.ExtractedData.StatisticalSummaries {
				if existing, exists := combined.StatisticalSummaries[metric]; exists {
					combined.StatisticalSummaries[metric] = (existing + value*weight) / 2.0
				} else {
					combined.StatisticalSummaries[metric] = value * weight
				}
			}
		}
	}

	// Normalize by total weight
	if totalWeight > 0 {
		for pop := range combined.AlleleFrequencies {
			for variant := range combined.AlleleFrequencies[pop] {
				combined.AlleleFrequencies[pop][variant] /= totalWeight
			}
		}

		for metric := range combined.StatisticalSummaries {
			combined.StatisticalSummaries[metric] /= totalWeight
		}
	}

	// Confidence = sample-size-weighted average quality score of contributing entries.
	// Returns 0 when no entries are provided or all have zero sample size.
	totalSizeWeight := 0.0
	weightedQualitySum := 0.0
	for _, entry := range entries {
		sizeW := float64(entry.SampleSize) / 1000.0
		totalSizeWeight += sizeW
		weightedQualitySum += entry.QualityScore * sizeW
	}
	if totalSizeWeight > 0 {
		combined.Confidence = weightedQualitySum / totalSizeWeight
	} else if len(entries) > 0 {
		// All entries have zero sample size; fall back to unweighted mean quality.
		sumQ := 0.0
		for _, entry := range entries {
			sumQ += entry.QualityScore
		}
		combined.Confidence = sumQ / float64(len(entries))
	} else {
		combined.Confidence = 0.0
	}

	return combined, nil
}

// createSyntheticParameters sets generation parameters. Apart from the
// allele frequencies (copied from the registry) and F_ST, the values are the
// same fixed numbers for every population; extracted is not used.
func (ape *AfricanPopulationEnhancer) createSyntheticParameters(basePop *H3AfricaPopulation, extracted *ExtractedGenomicsData) *SyntheticGenParams {
	// Copy allele frequencies from the source population so that
	// validateEnhancedModel can compare the model's encoded target against
	// the reference without the circular ±5% noise stub.
	alleleFreqs := make(map[string]float64, len(basePop.AlleleFrequencies))
	for id, freq := range basePop.AlleleFrequencies {
		alleleFreqs[id] = freq
	}

	params := &SyntheticGenParams{
		VariantDensity:    0.00125, // fixed; matches the synthetic entries' nucleotide diversity
		TransitionTsv:     2.1,     // human germline Ts/Tv (~2.0–2.1)
		IndelRatio:        0.12,    // fixed, not population-specific
		StructuralVarRate: 0.0008,  // SVs per kb, fixed
		RecombinationMap:  make(map[string][]float64),
		MutationRates:     make(map[string]float64),
		SelectionCoeff:    make(map[string]float64),
		AlleleFreqs:       alleleFreqs,
	}

	// Mutation rates, the same for every population
	params.MutationRates["snp"] = 1.2e-8   // Per base per generation
	params.MutationRates["indel"] = 1.5e-9 // Per base per generation
	params.MutationRates["sv"] = 8.0e-6    // Per kb per generation

	// Illustrative selection coefficients for three variants (not estimates)
	params.SelectionCoeff["rs4988235"] = -0.02 // LCT lactase persistence (slightly deleterious in traditional African diet)
	params.SelectionCoeff["rs334"] = 0.15      // Sickle cell trait (balancing selection)
	params.SelectionCoeff["rs1815739"] = 0.0   // ACTN3 (neutral in African context)

	// Population structure: F_ST from africanFST; admixture proportions and
	// subpopulation sizes are fixed placeholders
	fst := africanFSTForPop(basePop.PopulationID)
	params.PopulationStructure = &PopulationStructure{
		Fst: map[string]float64{
			"EUR": fst["EUR"],
			"EAS": fst["EAS"],
			"AMR": fst["AMR"],
			"SAS": fst["SAS"],
		},
		AdmixtureProportions: map[string]float64{
			"african":        0.98,
			"european":       0.01,
			"middle_eastern": 0.01,
		},
		IsolationByDistance: true,
		SubpopulationSizes: map[string]int{
			"rural": 60000,
			"urban": 40000,
		},
	}

	return params
}

// validateEnhancedModel computes five coverage-style metrics. Allele
// frequency accuracy compares params.AlleleFreqs with the registry it was
// copied from, so it is 1.0 whenever any frequency is non-zero; the other
// metrics count characterised epochs, LD blocks, SVs and PGx variants and are
// 0 when the registry has none. Nothing is compared with independent data.
func (ape *AfricanPopulationEnhancer) validateEnhancedModel(basePop *H3AfricaPopulation, params *SyntheticGenParams) *PopValidationMetrics {
	metrics := &PopValidationMetrics{
		ComparisonMetrics: make(map[string]float64),
	}

	// Validate allele frequency accuracy: compare params.AlleleFreqs (the model's
	// encoded target frequencies) against basePop.AlleleFrequencies (the reference).
	// Skip zero-frequency variants to avoid division by zero.
	accuracySum := 0.0
	variantCount := 0

	for variantID, expectedFreq := range basePop.AlleleFrequencies {
		if expectedFreq <= 0 {
			continue
		}
		predictedFreq, ok := params.AlleleFreqs[variantID]
		if !ok {
			continue
		}
		accuracy := 1.0 - math.Min(1.0, math.Abs(expectedFreq-predictedFreq)/expectedFreq)
		accuracySum += accuracy
		variantCount++
	}

	if variantCount > 0 {
		metrics.AlleleFreqAccuracy = accuracySum / float64(variantCount)
	}
	metrics.ComparisonMetrics["allele_freq_variants"] = float64(variantCount)

	// Validate demographic fit: requires DemographicParams with population epoch data.
	// Without it, fit cannot be estimated; returns 0.0 so callers know it is unmeasured.
	if basePop.DemographicParams != nil && len(basePop.DemographicParams.EffectivePopulationSize) > 0 {
		// Proxy: score scales with the number of characterised population epochs (capped at 1.0).
		// A more complete demographic model covers more time periods.
		epochCount := float64(len(basePop.DemographicParams.EffectivePopulationSize))
		const minEpochs = 3.0
		metrics.DemographicFit = math.Min(1.0, epochCount/minEpochs)
		metrics.ComparisonMetrics["demographic_epochs"] = epochCount
	} else {
		metrics.DemographicFit = 0.0
		log.Printf("validateEnhancedModel: DemographicFit not computed for %s — no DemographicParams available", basePop.PopulationID)
	}

	// Validate linkage accuracy: requires LinkagePatterns characterising LD blocks.
	if len(basePop.LinkagePatterns) > 0 {
		// Proxy: score scales with the number of characterised LD blocks (capped at 1.0).
		blockCount := float64(len(basePop.LinkagePatterns))
		const minBlocks = 5.0
		metrics.LinkageAccuracy = math.Min(1.0, blockCount/minBlocks)
		metrics.ComparisonMetrics["linkage_blocks"] = blockCount
	} else {
		metrics.LinkageAccuracy = 0.0
		log.Printf("validateEnhancedModel: LinkageAccuracy not computed for %s — no LinkagePatterns available", basePop.PopulationID)
	}

	// Validate structural variant accuracy: compare observed SV density with model parameter.
	if len(basePop.StructuralVariants) > 0 && params.StructuralVarRate > 0 {
		const genomeKb = 3_000_000.0 // approximate human genome size in kb
		observedRate := float64(len(basePop.StructuralVariants)) / genomeKb
		ratio := observedRate / params.StructuralVarRate
		// Score decays from 1.0 as the log-ratio diverges from zero.
		metrics.StructuralVarAccuracy = math.Max(0, 1.0-math.Abs(math.Log(ratio)))
		metrics.ComparisonMetrics["structural_variants"] = float64(len(basePop.StructuralVariants))
	} else {
		metrics.StructuralVarAccuracy = 0.0
		log.Printf("validateEnhancedModel: StructuralVarAccuracy not computed for %s — no StructuralVariants or zero StructuralVarRate", basePop.PopulationID)
	}

	// Validate pharmacogenomic fit: fraction of PGx variants covered in allele frequency data.
	if len(basePop.PharmacogenomicVariants) > 0 {
		covered := 0
		for _, pgx := range basePop.PharmacogenomicVariants {
			if pgx == nil {
				continue
			}
			if _, ok := basePop.AlleleFrequencies[pgx.VariantID]; ok {
				covered++
			}
		}
		metrics.PharmacogenomicFit = float64(covered) / float64(len(basePop.PharmacogenomicVariants))
		metrics.ComparisonMetrics["pgx_variants"] = float64(len(basePop.PharmacogenomicVariants))
		metrics.ComparisonMetrics["pgx_variants_covered"] = float64(covered)
	} else {
		metrics.PharmacogenomicFit = 0.0
		log.Printf("validateEnhancedModel: PharmacogenomicFit not computed for %s — no PharmacogenomicVariants available", basePop.PopulationID)
	}

	// Calculate overall validation score as mean of all five metrics.
	// Metrics that could not be computed contribute 0.0, which is an honest
	// representation of an unmeasured value.
	scores := []float64{
		metrics.AlleleFreqAccuracy,
		metrics.DemographicFit,
		metrics.LinkageAccuracy,
		metrics.StructuralVarAccuracy,
		metrics.PharmacogenomicFit,
	}

	sum := 0.0
	for _, score := range scores {
		sum += score
	}
	metrics.OverallValidationScore = sum / float64(len(scores))

	return metrics
}

// assessModelQuality relabels four validation metrics (e.g. DemographicFit as
// BiologicalRealism) and averages them.
func (ape *AfricanPopulationEnhancer) assessModelQuality(params *SyntheticGenParams, validation *PopValidationMetrics) *QualityMetrics {
	quality := &QualityMetrics{
		BiologicalRealism:     validation.DemographicFit,
		StatisticalFidelity:   validation.AlleleFreqAccuracy,
		PopulationSpecificity: validation.LinkageAccuracy,
		ClinicalRelevance:     validation.PharmacogenomicFit,
	}

	// Calculate overall quality
	scores := []float64{
		quality.BiologicalRealism,
		quality.StatisticalFidelity,
		quality.PopulationSpecificity,
		quality.ClinicalRelevance,
	}

	sum := 0.0
	for _, score := range scores {
		sum += score
	}
	quality.OverallQuality = sum / float64(len(scores))

	return quality
}

// GenerateAfricanSyntheticData generates sampleSize toy samples for a
// population that EnhanceAfricanPopulations has modelled. The samples are
// not validated; ValidationStatus is always "unvalidated".
func (ape *AfricanPopulationEnhancer) GenerateAfricanSyntheticData(ctx context.Context, populationID string, sampleSize int) (*AfricanSyntheticDataset, error) {
	model, exists := ape.populationModels[populationID]
	if !exists {
		return nil, fmt.Errorf("population model not found for %s", populationID)
	}

	startTime := time.Now()

	// Generate synthetic samples
	samples, err := ape.generateSyntheticSamples(model, sampleSize)
	if err != nil {
		return nil, fmt.Errorf("failed to generate synthetic samples: %w", err)
	}

	// Create result dataset
	dataset := &AfricanSyntheticDataset{
		PopulationID:     populationID,
		PopulationModel:  model,
		SampleSize:       sampleSize,
		GeneratedSamples: samples,
		GenerationTime:   time.Since(startTime),
		QualityMetrics:   model.QualityAssessment,
		ValidationStatus: "unvalidated",
	}

	// Record generation event
	event := &GenerationEvent{
		Timestamp:       time.Now(),
		EventType:       "synthetic_data_generation",
		Parameters:      map[string]interface{}{"sample_size": sampleSize},
		DataGenerated:   sampleSize,
		ValidationScore: model.QualityAssessment.OverallQuality,
	}

	model.GenerationHistory = append(model.GenerationHistory, event)

	return dataset, nil
}

// AfricanSyntheticDataset represents generated synthetic data for African populations
type AfricanSyntheticDataset struct {
	PopulationID     string                    `json:"population_id"`
	PopulationModel  *AfricanPopulationModel   `json:"population_model"`
	SampleSize       int                       `json:"sample_size"`
	GeneratedSamples []*SyntheticAfricanSample `json:"generated_samples"`
	GenerationTime   time.Duration             `json:"generation_time"`
	QualityMetrics   *QualityMetrics           `json:"quality_metrics"`
	ValidationStatus string                    `json:"validation_status"`
}

// SyntheticAfricanSample represents a synthetic sample from an African population
type SyntheticAfricanSample struct {
	SampleID         string                    `json:"sample_id"`
	PopulationID     string                    `json:"population_id"`
	GeographicOrigin []float64                 `json:"geographic_origin"`
	Genotype         map[string]*SampleVariant `json:"genotype"`
	PGxProfile       []*PGxVariant             `json:"pharmacogenomic_profile"`
	Ancestry         map[string]float64        `json:"ancestry_proportions"`
	QualityScore     float64                   `json:"quality_score"`
}

// SampleVariant represents a variant in a synthetic sample
type SampleVariant struct {
	VariantID       string  `json:"variant_id"`
	Chromosome      string  `json:"chromosome"`
	Position        int64   `json:"position"`
	ReferenceAllele string  `json:"reference_allele"`
	AlternateAllele string  `json:"alternate_allele"`
	Zygosity        string  `json:"zygosity"`
	QualityScore    float64 `json:"quality_score"`
	PopulationFreq  float64 `json:"population_frequency"`
}

// generateSyntheticSamples includes each registry variant in a sample with
// probability equal to its allele frequency, then draws its zygosity
// independently from Hardy-Weinberg proportions, so an included variant can
// be "0/0".
func (ape *AfricanPopulationEnhancer) generateSyntheticSamples(model *AfricanPopulationModel, sampleSize int) ([]*SyntheticAfricanSample, error) {
	var samples []*SyntheticAfricanSample

	for i := 0; i < sampleSize; i++ {
		sample := &SyntheticAfricanSample{
			SampleID:         fmt.Sprintf("SYN_%s_%04d", model.Population.PopulationID, i+1),
			PopulationID:     model.Population.PopulationID,
			GeographicOrigin: model.Population.GeographicCoords,
			Genotype:         make(map[string]*SampleVariant),
			Ancestry:         make(map[string]float64),
		}

		// Generate genotype based on population allele frequencies
		for variant, freq := range model.Population.AlleleFrequencies {
			if ape.randomFloat() < freq {
				ref, alt := ape.generateAlleles(variant)
				sampleVariant := &SampleVariant{
					VariantID:       variant,
					Chromosome:      ape.getVariantChromosome(variant),
					Position:        ape.getVariantPosition(variant),
					ReferenceAllele: ref,
					AlternateAllele: alt,
					Zygosity:        ape.generateZygosity(freq),
					QualityScore:    0.95,
					PopulationFreq:  freq,
				}
				sample.Genotype[variant] = sampleVariant
			}
		}

		// Add pharmacogenomic variants
		if model.Population.PharmacogenomicVariants != nil {
			sample.PGxProfile = append([]*PGxVariant{}, model.Population.PharmacogenomicVariants...)
		}

		// Set ancestry proportions
		if model.SyntheticParameters.PopulationStructure != nil {
			sample.Ancestry = make(map[string]float64)
			for ancestry, proportion := range model.SyntheticParameters.PopulationStructure.AdmixtureProportions {
				sample.Ancestry[ancestry] = proportion + (ape.randomFloat()-0.5)*0.1 // Add individual variation
			}
		}

		sample.QualityScore = model.QualityAssessment.OverallQuality
		samples = append(samples, sample)
	}

	return samples, nil
}

// Helper methods

// variantAlleles maps known rsIDs to [ref, alt] alleles as usually reported.
// rs4988235 (LCT) and rs776746 (CYP3A5) are given on the gene (minus)
// strand, not the GRCh38 forward strand.
var variantAlleles = map[string][2]string{
	"rs429358":  {"T", "C"}, // APOE ε4  — chr19:44908684
	"rs7412":    {"C", "T"}, // APOE ε2  — chr19:44908822
	"rs1815739": {"C", "T"}, // ACTN3 R577X
	"rs4988235": {"C", "T"}, // LCT upstream (lactase persistence)
	"rs4244285": {"G", "A"}, // CYP2C19*2 splicing variant
	"rs776746":  {"A", "G"}, // CYP3A5*3 splicing variant
}

// generateAlleles returns the REF and ALT allele for a variant.
// Known rsIDs return the alleles in variantAlleles.  Unknown variants use
// Ts/Tv-weighted random selection (ratio ≈ 2.1, the observed human germline value).
func (ape *AfricanPopulationEnhancer) generateAlleles(variantID string) (ref, alt string) {
	if alleles, ok := variantAlleles[variantID]; ok {
		return alleles[0], alleles[1]
	}
	bases := [4]byte{'A', 'T', 'G', 'C'}
	refBase := bases[int(ape.randomFloat()*4)%4]
	altBase := ape.pickAlternateAllele(refBase)
	return string([]byte{refBase}), string([]byte{altBase})
}

// pickAlternateAllele selects an alternate base using the empirical human
// germline Ts/Tv ratio of 2.1 (transitions are 2.1× more likely than each
// transversion class).
func (ape *AfricanPopulationEnhancer) pickAlternateAllele(ref byte) byte {
	const tsTv = 2.1
	type altWeight struct {
		base   byte
		weight float64
	}
	var alts []altWeight
	switch ref {
	case 'A':
		alts = []altWeight{{'G', tsTv}, {'C', 1.0}, {'T', 1.0}} // G = transition
	case 'T':
		alts = []altWeight{{'C', tsTv}, {'A', 1.0}, {'G', 1.0}} // C = transition
	case 'G':
		alts = []altWeight{{'A', tsTv}, {'C', 1.0}, {'T', 1.0}} // A = transition
	default: // 'C'
		alts = []altWeight{{'T', tsTv}, {'A', 1.0}, {'G', 1.0}} // T = transition
	}
	total := 0.0
	for _, a := range alts {
		total += a.weight
	}
	r := ape.randomFloat() * total
	cum := 0.0
	for _, a := range alts {
		cum += a.weight
		if r <= cum {
			return a.base
		}
	}
	return alts[len(alts)-1].base
}

func (ape *AfricanPopulationEnhancer) randomFloat() float64 {
	if ape.rng == nil {
		ape.rng = rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)) // #nosec G404
	}
	return ape.rng.Float64()
}

func (ape *AfricanPopulationEnhancer) getVariantChromosome(variantID string) string {
	// Simplified mapping
	variantChromosomes := map[string]string{
		"rs429358":  "19", // APOE
		"rs7412":    "19", // APOE
		"rs1815739": "11", // ACTN3
		"rs4988235": "2",  // LCT
		"rs4244285": "10", // CYP2C19
		"rs776746":  "7",  // CYP3A5
	}

	if chr, exists := variantChromosomes[variantID]; exists {
		return chr
	}
	return "1" // Default
}

func (ape *AfricanPopulationEnhancer) getVariantPosition(variantID string) int64 {
	// Positions mix assemblies: rs429358, rs7412, rs4988235 and rs776746 are
	// GRCh37; rs1815739 and rs4244285 are GRCh38
	variantPositions := map[string]int64{
		"rs429358":  45411941,  // APOE
		"rs7412":    45412079,  // APOE
		"rs1815739": 66560624,  // ACTN3
		"rs4988235": 136608646, // LCT
		"rs4244285": 94781859,  // CYP2C19
		"rs776746":  99270539,  // CYP3A5
	}

	if pos, exists := variantPositions[variantID]; exists {
		return pos
	}
	return 1000000 // Default
}

func (ape *AfricanPopulationEnhancer) generateZygosity(freq float64) string {
	// Zygosity drawn from Hardy-Weinberg genotype proportions
	r := ape.randomFloat()

	// Assuming biallelic variants
	homRef := (1 - freq) * (1 - freq)
	het := 2 * freq * (1 - freq)

	if r < homRef {
		return "0/0" // Homozygous reference
	} else if r < homRef+het {
		return "0/1" // Heterozygous
	} else {
		return "1/1" // Homozygous alternate
	}
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
