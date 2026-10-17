package authenticity_test

// coverage_boost_test.go exercises code paths not covered by the primary test file.

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ============================================================
// Helper - build a large sequences map to pass minimum sample-size
// requirements (SampleSizeRequired / 2 = 50 for PCA, 100 for FST).
// ============================================================

func buildLargeSeqs(n int) map[genomics.SequencingPlatform]map[population.Population][]string {
	bases := []string{
		"ATCGATCGATCG", "GCTAGCTAGCTA", "TTTTAAAACCCC", "GGGGCCCCATCG",
		"AAATTTGGGCCC", "CCCTTTAAAGGG", "GCGCATATATGC", "ATATGCGCATCG",
		"AACCGGTTAACCGGTT", "GGTTAACCGGTTAACC",
	}
	build := func(count int) []string {
		out := make([]string, count)
		for i := range out {
			out[i] = bases[i%len(bases)]
		}
		return out
	}
	return map[genomics.SequencingPlatform]map[population.Population][]string{
		genomics.NovaSeq6000: {
			population.YRI: build(n),
			population.CEU: build(n),
		},
		genomics.PacBioRevio: {
			population.CHB: build(n),
			population.MXL: build(n),
		},
	}
}

// ============================================================
// IndistinguishabilityTester – larger sample coverage
// ============================================================

func TestIndistinguishabilityTester_LargeSample_AllTestsPassing(t *testing.T) {
	it := authenticity.NewIndistinguishabilityTester()
	seqs := buildLargeSeqs(110) // > 100 = SampleSizeRequired for PCA

	result, err := it.TestCrossPlatformIndistinguishability(seqs, 110)
	require.NoError(t, err)
	require.NotNil(t, result)

	// With 110 samples all classification tests should be able to run.
	assert.True(t, len(result.PassedTests)+len(result.FailedTests) > 0)
	assert.NotEmpty(t, result.StatisticalSummary.MultipleTestingCorrection)
}

func TestIndistinguishabilityTester_LargeSample_FST(t *testing.T) {
	// FST requires SampleSizeRequired/2 = 100 samples
	it := authenticity.NewIndistinguishabilityTester()
	seqs := buildLargeSeqs(210)

	result, err := it.TestCrossPlatformIndistinguishability(seqs, 210)
	require.NoError(t, err)
	require.NotNil(t, result)

	for _, score := range result.PlatformScores {
		assert.GreaterOrEqual(t, score, 0.0)
		assert.LessOrEqual(t, score, 1.0)
	}
}

// ============================================================
// AdversarialDefenseSystem – additional coverage paths
// ============================================================

func TestAdversarialDefenseSystem_RecommendedDefenses_NotNil(t *testing.T) {
	ads := authenticity.NewAdversarialDefenseSystem()
	seqs := makeSeqs()

	result, err := ads.TestAdversarialResistance(seqs, 5)
	require.NoError(t, err)
	assert.NotNil(t, result.RecommendedDefenses)
}

func TestAdversarialDefenseSystem_DefenseEffectiveness_AllStrategies(t *testing.T) {
	ads := authenticity.NewAdversarialDefenseSystem()
	seqs := buildLargeSeqs(10)

	result, err := ads.TestAdversarialResistance(seqs, 5)
	require.NoError(t, err)

	// Defense effectiveness map should have entries for initialised strategies.
	assert.NotEmpty(t, result.DefenseEffectiveness)
}

func TestAdversarialDefenseSystem_ModelSpecificResults_Keys(t *testing.T) {
	ads := authenticity.NewAdversarialDefenseSystem()
	seqs := makeSeqs()

	result, err := ads.TestAdversarialResistance(seqs, 3)
	require.NoError(t, err)

	// All 5 adversarial models should appear.
	assert.True(t, len(result.ModelSpecificResults) == 5, "expected 5 models, got %d", len(result.ModelSpecificResults))
}

func TestAdversarialDefenseSystem_SecurityLevel_AllRanges(t *testing.T) {
	// Run with a tiny sample to exercise a particular score range
	ads := authenticity.NewAdversarialDefenseSystem()
	seqs := makeSeqs()

	result, err := ads.TestAdversarialResistance(seqs, 3)
	require.NoError(t, err)

	validLevels := map[string]bool{
		"Maximum_Security":      true,
		"High_Security":         true,
		"Moderate_Security":     true,
		"Low_Security":          true,
		"Critical_Vulnerability": true,
	}
	assert.True(t, validLevels[result.SecurityLevel])
}

// ============================================================
// FingerprintObfuscationSystem – exercise all strategy types
// ============================================================

func TestFingerprintObfuscation_AllStrategiesExercised(t *testing.T) {
	fos := authenticity.NewFingerprintObfuscationSystem()
	fos.InitializeObfuscationStrategies()
	fos.InitializeFingerprintAnalyzers()
	fos.InitializeCrossPlatformModels()

	seqs := buildLargeSeqs(5)

	results, err := fos.ObfuscatePlatformFingerprints(seqs, 0.8)
	require.NoError(t, err)
	require.NotNil(t, results)

	for _, result := range *results {
		assert.NotNil(t, result.AnalysisResults)
		assert.GreaterOrEqual(t, result.QualityPreservation, 0.0)
		assert.LessOrEqual(t, result.QualityPreservation, 1.0)
		assert.GreaterOrEqual(t, result.PopulationPreservation, 0.0)
		assert.LessOrEqual(t, result.PopulationPreservation, 1.0)
	}
}

func TestFingerprintObfuscation_GetterMethods(t *testing.T) {
	fos := authenticity.NewFingerprintObfuscationSystem()

	// Before init – maps are empty but not nil
	assert.NotNil(t, fos.GetObfuscationStrategies())
	assert.NotNil(t, fos.GetFingerprintAnalyzers())
	assert.NotNil(t, fos.GetCrossPlatformModels())

	// After init – maps are populated
	fos.InitializeObfuscationStrategies()
	fos.InitializeFingerprintAnalyzers()
	fos.InitializeCrossPlatformModels()

	assert.NotEmpty(t, fos.GetObfuscationStrategies())
	assert.NotEmpty(t, fos.GetFingerprintAnalyzers())
	assert.NotEmpty(t, fos.GetCrossPlatformModels())
}

func TestFingerprintObfuscation_MultiPlatformSequences(t *testing.T) {
	fos := authenticity.NewFingerprintObfuscationSystem()
	seqs := map[genomics.SequencingPlatform]map[population.Population][]string{
		genomics.NovaSeq6000: {
			population.YRI: {"ATCGATCGATCG", "GCTAGCTAGCTA"},
			population.CEU: {"GGGGCCCCATCG", "CCCCATCGATCG"},
		},
		genomics.PromethION24: {
			population.CHB: {"AAAAAAAAGGGGGGGG", "TTTTTTTTCCCCCCCC"},
		},
		genomics.PacBioRevio: {
			population.MXL: {"ATCGATCGATCGATCG", "GCTAGCTAGCTAGCTA"},
		},
		genomics.MiSeq: {
			population.YRI: {"ATCG", "GCTA"},
		},
		genomics.NextSeq2000: {
			population.CEU: {"ATCGATCGATCGATCGATCG"},
		},
		genomics.DNBSeqT20x2: {
			population.CHB: {"GCTAGCTAGCTA"},
		},
		genomics.ElementAVITI: {
			population.MXL: {"ATCGATCGATCG"},
		},
	}

	results, err := fos.ObfuscatePlatformFingerprints(seqs, 0.9)
	require.NoError(t, err)
	require.NotNil(t, results)
	assert.Equal(t, 7, len(*results))
}

// ============================================================
// FinalValidationSystem – all protocol levels
// ============================================================

func TestFinalValidationSystem_Advanced_Protocol(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	seqs := buildLargeSeqs(50)

	result, err := fvs.PerformComprehensiveValidation(seqs, "Advanced_Comprehensive_Validation", "Gold_Advanced_Compliance")
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Equal(t, "Advanced", result.ValidationLevel)
}

func TestFinalValidationSystem_Expert_Protocol(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	seqs := buildLargeSeqs(20)

	result, err := fvs.PerformComprehensiveValidation(seqs, "Expert_Research_Grade_Validation", "Platinum_Expert_Compliance")
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Equal(t, "Expert", result.ValidationLevel)
}

func TestFinalValidationSystem_Certification_Protocol(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	seqs := buildLargeSeqs(20)

	result, err := fvs.PerformComprehensiveValidation(seqs, "Certification_Formal_Validation", "Diamond_Ultimate_Compliance")
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Equal(t, "Certification", result.ValidationLevel)
}

func TestFinalValidationSystem_GetterMethods_AfterInit(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	fvs.InitializeValidationProtocols()
	fvs.InitializeBenchmarkDatasets()
	fvs.InitializeCertificationTests()

	protocols := fvs.GetValidationProtocols()
	datasets := fvs.GetBenchmarkDatasets()
	tests := fvs.GetCertificationTests()

	assert.NotEmpty(t, protocols)
	assert.NotEmpty(t, datasets)
	assert.NotEmpty(t, tests)
}

func TestFinalValidationSystem_QualityGrade_ExpertProtocol(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	seqs := buildLargeSeqs(15)

	result, err := fvs.PerformComprehensiveValidation(seqs, "Expert_Research_Grade_Validation", "Diamond_Ultimate_Compliance")
	require.NoError(t, err)
	require.NotNil(t, result)

	validGrades := map[string]bool{
		"A+": true, "A": true, "B+": true, "B": true,
		"C+": true, "C": true, "D": true, "F": true,
	}
	assert.True(t, validGrades[result.QualityAssessment.OverallQualityGrade])
}

func TestFinalValidationSystem_RecommendationsReport_Fields(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	seqs := makeSeqs()

	result, err := fvs.PerformComprehensiveValidation(seqs, "Basic_Statistical_Indistinguishability", "Bronze_Basic_Compliance")
	require.NoError(t, err)
	require.NotNil(t, result)

	rr := result.RecommendationsReport
	// These slices may be empty but must not be nil
	assert.NotNil(t, rr.HighPriorityRecommendations)
	assert.NotNil(t, rr.MediumPriorityRecommendations)
	assert.NotNil(t, rr.LowPriorityRecommendations)
	assert.NotNil(t, rr.TechnicalImprovements)
	assert.NotNil(t, rr.StrategicInitiatives)
	assert.NotNil(t, rr.ComplianceRequirements)
}

func TestFinalValidationSystem_StatisticalSummary_Populated(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	seqs := buildLargeSeqs(10)

	result, err := fvs.PerformComprehensiveValidation(seqs, "Basic_Statistical_Indistinguishability", "Bronze_Basic_Compliance")
	require.NoError(t, err)
	require.NotNil(t, result)

	ss := result.StatisticalSummary
	assertFiniteOrBounded(t, ss.StatisticalIndistinguishability, 0.0, 1.0)
	assertFiniteOrBounded(t, ss.CrossPlatformConsistency, 0.0, 1.0)
	assertFiniteOrBounded(t, ss.BiologicalRealismScore, 0.0, 1.0)
	assert.NotNil(t, ss.PlatformSpecificRates)
	assert.NotNil(t, ss.PopulationSpecificRates)
}

// ============================================================
// AdvancedValidation2025 – empty data path
// ============================================================

func TestAdvancedValidation2025_EmptyData(t *testing.T) {
	av := authenticity.NewAdvancedValidation2025()
	emptyGenerated := map[genomics.SequencingPlatform][]authenticity.GenomicSequence{}
	emptyReal := map[genomics.SequencingPlatform][]authenticity.GenomicSequence{}

	result, err := av.ValidateAgainstStateOfTheArt2025(emptyGenerated, emptyReal, []population.Population{population.YRI})
	require.NoError(t, err)
	require.NotNil(t, result)
	// All sub-validators use placeholder returns, so no error expected even with empty data.
	assert.NotNil(t, result.DPGANValidation)
}

// ============================================================
// AdvancedValidation2025 – multiple populations
// ============================================================

func TestAdvancedValidation2025_MultiplePopulations(t *testing.T) {
	av := authenticity.NewAdvancedValidation2025()

	generated := map[genomics.SequencingPlatform][]authenticity.GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "ATCGATCG", QualityScores: []int{35, 35, 35, 35, 35, 35, 35, 35}, Population: population.YRI},
			{Sequence: "GCTAGCTA", QualityScores: []int{30, 30, 30, 30, 30, 30, 30, 30}, Population: population.CEU},
			{Sequence: "GGGGCCCC", QualityScores: []int{32, 32, 32, 32, 32, 32, 32, 32}, Population: population.CHB},
		},
		genomics.PacBioRevio: {
			{Sequence: "AAAATTTT", QualityScores: []int{25, 25, 25, 25, 25, 25, 25, 25}, Population: population.MXL},
		},
	}
	real := map[genomics.SequencingPlatform][]authenticity.GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "TTTTAAAA", QualityScores: []int{38, 38, 38, 38, 38, 38, 38, 38}, Population: population.YRI},
		},
	}
	targetPops := []population.Population{population.YRI, population.CEU, population.CHB, population.MXL}

	result, err := av.ValidateAgainstStateOfTheArt2025(generated, real, targetPops)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.NotNil(t, result.ThreatModelCompliance)
	assert.NotNil(t, result.ThreatModelCompliance.RecommendedImprovements)
}

// ============================================================
// FinalValidationSystem – different certification targets
// ============================================================

func TestFinalValidationSystem_BronzeCertification(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	seqs := makeSeqs()

	result, err := fvs.PerformComprehensiveValidation(seqs, "Basic_Statistical_Indistinguishability", "Bronze_Basic_Compliance")
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.NotEmpty(t, result.CertificationStatus.CertificationStandard)
	assert.NotEmpty(t, result.CertificationStatus.CertificationValidity)
	assert.NotEmpty(t, result.CertificationStatus.CertificationAuthority)
}

func TestFinalValidationSystem_DiamondCertification(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	seqs := makeSeqs()

	result, err := fvs.PerformComprehensiveValidation(seqs, "Basic_Statistical_Indistinguishability", "Diamond_Ultimate_Compliance")
	require.NoError(t, err)
	require.NotNil(t, result)
	// The score may not meet Diamond requirements, but the system should evaluate it.
	assert.NotNil(t, result.CertificationStatus)
}

func TestFinalValidationSystem_UnknownCertificationTarget(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	seqs := makeSeqs()

	// When target certification is unknown, system falls back to Bronze
	result, err := fvs.PerformComprehensiveValidation(seqs, "Basic_Statistical_Indistinguishability", "Unknown_Cert_Level")
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.NotNil(t, result.CertificationStatus)
}

// ============================================================
// IndistinguishabilityTester – population pair scoring edge cases
// ============================================================

func TestIndistinguishabilityTester_SinglePlatformMultiPopulation(t *testing.T) {
	it := authenticity.NewIndistinguishabilityTester()
	seqs := map[genomics.SequencingPlatform]map[population.Population][]string{
		genomics.NovaSeq6000: {
			population.YRI: {"ATCGATCGATCG", "GCTAGCTAGCTA", "TTTTAAAACCCC"},
			population.CEU: {"GCTAGCTAGCTA", "ATCGATCGATCG", "AAAATTTTGGGG"},
			population.CHB: {"GGGGCCCCATCG", "CCCCATCGATCG", "TTTAAAGGGCCC"},
		},
	}

	result, err := it.TestCrossPlatformIndistinguishability(seqs, 3)
	require.NoError(t, err)
	require.NotNil(t, result)
	// 3 populations -> 3 pairs
	assert.Equal(t, 3, len(result.PopulationPairScores))
}

func TestIndistinguishabilityTester_RecommendationsContent(t *testing.T) {
	it := authenticity.NewIndistinguishabilityTester()
	seqs := makeSeqs()

	result, err := it.TestCrossPlatformIndistinguishability(seqs, 3)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.NotEmpty(t, result.RecommendedActions)
}

// ============================================================
// Comprehensive result data structure checks
// ============================================================

func TestAdvancedValidationResult2025_AllFields(t *testing.T) {
	av := authenticity.NewAdvancedValidation2025()
	generated := map[genomics.SequencingPlatform][]authenticity.GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "ATCGATCG", QualityScores: []int{35, 35, 35, 35, 35, 35, 35, 35}, Population: population.YRI},
		},
	}
	real := map[genomics.SequencingPlatform][]authenticity.GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "GCTAGCTA", QualityScores: []int{30, 30, 30, 30, 30, 30, 30, 30}, Population: population.CEU},
		},
	}

	result, err := av.ValidateAgainstStateOfTheArt2025(generated, real, []population.Population{population.YRI, population.CEU})
	require.NoError(t, err)
	require.NotNil(t, result)

	// Check all sub-result structs are populated
	assert.NotNil(t, result.TabPFNValidation.VulnerabilityReport)
	assert.NotNil(t, result.TabPFNValidation.VulnerabilityReport.Vulnerabilities)
	assert.NotNil(t, result.TabPFNValidation.VulnerabilityReport.Recommendations)

	assert.GreaterOrEqual(t, result.DPGANValidation.MembershipInferenceResistance, 0.0)
	assert.GreaterOrEqual(t, result.HaploNetValidation.StructureInferenceResistance, 0.0)
	assert.GreaterOrEqual(t, result.FederatedValidation.GradientInversionResistance, 0.0)
	assert.GreaterOrEqual(t, result.CrossPlatformValidation.BEADSEffectiveness, 0.0)
	assert.GreaterOrEqual(t, result.EnsembleValidation.StackingResistance, 0.0)
}

// ============================================================
// Obfuscation result field checks
// ============================================================

func TestObfuscationResult_AnalysisResultsPopulated(t *testing.T) {
	fos := authenticity.NewFingerprintObfuscationSystem()
	seqs := map[genomics.SequencingPlatform]map[population.Population][]string{
		genomics.NovaSeq6000: {
			population.YRI: {"ATCGATCGATCG", "GCTAGCTAGCTA"},
		},
	}

	results, err := fos.ObfuscatePlatformFingerprints(seqs, 0.9)
	require.NoError(t, err)
	require.NotNil(t, results)

	for _, result := range *results {
		assert.NotNil(t, result.AnalysisResults)
		// There should be one result per fingerprint analyzer
		assert.NotEmpty(t, result.AnalysisResults)

		for analyzerName, analysisResult := range result.AnalysisResults {
			assert.NotEmpty(t, analyzerName)
			assert.NotNil(t, analysisResult)
			assert.NotEmpty(t, analysisResult.AnalyzerName)
		}
	}
}

// ============================================================
// CrossPlatformModel fields
// ============================================================

func TestCrossPlatformModels_UnifiedProfiles(t *testing.T) {
	fos := authenticity.NewFingerprintObfuscationSystem()
	fos.InitializeCrossPlatformModels()

	models := fos.GetCrossPlatformModels()
	require.NotEmpty(t, models)

	universalModel, exists := models["Universal_Neutral_Model"]
	require.True(t, exists)
	assert.NotEmpty(t, universalModel.ModelName)
	assert.Greater(t, universalModel.IndistinguishabilityScore, 0.0)
	assert.NotEmpty(t, universalModel.PlatformCompatibility)
	assert.NotEmpty(t, universalModel.UnifiedErrorProfile.SubstitutionPattern)
	assert.NotEmpty(t, universalModel.UnifiedQualityProfile.CycleDecayPattern)
}

// ============================================================
// BenchmarkDataset fields
// ============================================================

func TestBenchmarkDataset_Fields(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	fvs.InitializeBenchmarkDatasets()

	datasets := fvs.GetBenchmarkDatasets()
	require.NotEmpty(t, datasets)

	ref, exists := datasets["Reference_Multi_Platform_Dataset"]
	require.True(t, exists)
	assert.Equal(t, "Reference", ref.DatasetType)
	assert.Greater(t, ref.TotalSequences, 0)
	assert.NotEmpty(t, ref.PlatformCoverage)
	assert.NotEmpty(t, ref.PopulationCoverage)
	assert.Greater(t, ref.QualityMetrics.DataCompleteness, 0.0)
}

// ============================================================
// CertificationTest fields
// ============================================================

func TestCertificationTest_Fields(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	fvs.InitializeCertificationTests()

	tests := fvs.GetCertificationTests()

	for level, test := range tests {
		assert.NotEmpty(t, level)
		assert.NotEmpty(t, test.TestName)
		assert.NotEmpty(t, test.CertificationLevel)
		assert.Greater(t, test.RequiredAccuracy, 0.0)
		assert.Greater(t, test.CertificationCriteria.MinimumAccuracy, 0.0)
	}
}

// ============================================================
// ValidationProtocol fields
// ============================================================

func TestValidationProtocol_Fields(t *testing.T) {
	fvs := authenticity.NewFinalValidationSystem()
	fvs.InitializeValidationProtocols()

	protocols := fvs.GetValidationProtocols()

	for name, protocol := range protocols {
		assert.NotEmpty(t, name)
		assert.NotEmpty(t, protocol.ProtocolName)
		assert.NotEmpty(t, protocol.ValidationLevel)
		assert.NotEmpty(t, protocol.TestCategories)
		assert.Greater(t, protocol.RequiredSampleSize, 0)
		assert.Greater(t, protocol.PassingThreshold, 0.0)
		assert.NotEmpty(t, protocol.StatisticalRigor.StatisticalTests)
		assert.NotEmpty(t, protocol.QualityAssurance.DataIntegrityChecks)
	}
}

// ============================================================
// DefenseStrategy fields
// ============================================================

func TestDefenseStrategy_Fields(t *testing.T) {
	ads := authenticity.NewAdversarialDefenseSystem()
	ads.InitializeDefenseStrategies()

	strategies := ads.GetDefenseStrategies()

	for name, strategy := range strategies {
		assert.NotEmpty(t, name)
		assert.NotEmpty(t, strategy.StrategyName)
		assert.NotEmpty(t, strategy.DefenseType)
		assert.NotEmpty(t, strategy.TargetedModels)
		assert.Greater(t, strategy.EffectivenessScore, 0.0)
		assert.NotEmpty(t, strategy.DefenseParameters)
	}
}

// ============================================================
// AdversarialModel fields
// ============================================================

func TestAdversarialModel_Fields(t *testing.T) {
	ads := authenticity.NewAdversarialDefenseSystem()
	ads.InitializeAdversarialModels()

	models := ads.GetAdversarialModels()

	for name, model := range models {
		assert.NotEmpty(t, name)
		assert.NotEmpty(t, model.ModelName)
		assert.NotEmpty(t, model.ModelType)
		assert.NotEmpty(t, model.FeatureSet)
		assert.Greater(t, model.TrainingAccuracy, 0.0)
		assert.Greater(t, model.ValidationAccuracy, 0.0)
		assert.NotEmpty(t, model.FeatureWeights)
		assert.NotEmpty(t, model.ModelParameters)
		assert.Greater(t, model.AttackCapability.MaxFeatureDimensions, 0)
	}
}

// ============================================================
// FingerprintAnalyzer fields
// ============================================================

func TestFingerprintAnalyzer_Fields(t *testing.T) {
	fos := authenticity.NewFingerprintObfuscationSystem()
	fos.InitializeFingerprintAnalyzers()

	analyzers := fos.GetFingerprintAnalyzers()

	for name, analyzer := range analyzers {
		assert.NotEmpty(t, name)
		assert.NotEmpty(t, analyzer.AnalyzerName)
		assert.NotEmpty(t, analyzer.AnalysisType)
		assert.NotEmpty(t, analyzer.DetectionMethod)
		assert.Greater(t, analyzer.Sensitivity, 0.0)
		assert.Greater(t, analyzer.Specificity, 0.0)
		assert.Greater(t, analyzer.AccuracyScore, 0.0)
		assert.NotEmpty(t, analyzer.FeatureWeights)
	}
}

// ============================================================
// ObfuscationStrategy fields
// ============================================================

func TestObfuscationStrategy_Fields(t *testing.T) {
	fos := authenticity.NewFingerprintObfuscationSystem()
	fos.InitializeObfuscationStrategies()

	strategies := fos.GetObfuscationStrategies()

	for name, strategy := range strategies {
		assert.NotEmpty(t, name)
		assert.NotEmpty(t, strategy.StrategyName)
		assert.NotEmpty(t, strategy.ObfuscationType)
		assert.NotEmpty(t, strategy.TargetPlatforms)
		assert.Greater(t, strategy.ObfuscationStrength, 0.0)
		assert.Greater(t, strategy.PreservationScore, 0.0)
		assert.Greater(t, strategy.DetectionResistance, 0.0)
		assert.NotEmpty(t, strategy.ObfuscationParameters)
	}
}

// ============================================================
// IndistinguishabilityResult field completeness
// ============================================================

func TestIndistinguishabilityResult_AllFields(t *testing.T) {
	it := authenticity.NewIndistinguishabilityTester()
	seqs := buildLargeSeqs(110)

	result, err := it.TestCrossPlatformIndistinguishability(seqs, 110)
	require.NoError(t, err)
	require.NotNil(t, result)

	assert.NotNil(t, result.PlatformScores)
	assert.NotNil(t, result.PopulationPairScores)
	assert.NotNil(t, result.PassedTests)
	assert.NotNil(t, result.FailedTests)
	assert.NotNil(t, result.RecommendedActions)

	summary := result.StatisticalSummary
	assert.GreaterOrEqual(t, summary.MeanClassificationRate, 0.0)
	assert.LessOrEqual(t, summary.MeanClassificationRate, 1.0)
	assert.GreaterOrEqual(t, summary.PercentagePassing, 0.0)
	assert.LessOrEqual(t, summary.PercentagePassing, 100.0)
}

// ============================================================
// ClassificationTest fields
// ============================================================

func TestClassificationTest_FieldValues(t *testing.T) {
	it := authenticity.NewIndistinguishabilityTester()
	it.InitializeClassificationTests()

	tests := it.GetClassificationTests()
	require.NotEmpty(t, tests)

	for key, ct := range tests {
		assert.NotEmpty(t, key)
		assert.NotEmpty(t, ct.TestName)
		assert.NotEmpty(t, ct.TestType)
		assert.Greater(t, ct.Sensitivity, 0.0)
		assert.Greater(t, ct.SpecificityThreshold, 0.0)
		assert.Greater(t, ct.PowerAnalysis.SampleSizeRequired, 0)
		assert.Greater(t, ct.PowerAnalysis.PowerAchieved, 0.0)
	}
}
