package genomics_test

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ---------------------------------------------------------------------------
// SequencingPlatform constants
// ---------------------------------------------------------------------------

func TestSequencingPlatform_Constants(t *testing.T) {
	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}
	assert.Len(t, platforms, 7)
	// All platform names should be non-empty strings
	for _, p := range platforms {
		assert.NotEmpty(t, string(p))
	}
}

// ---------------------------------------------------------------------------
// MultiPlatformQualityModeler
// ---------------------------------------------------------------------------

func TestNewMultiPlatformQualityModeler(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	require.NotNil(t, m)
}

func TestMultiPlatformQualityModeler_GetSupportedPlatforms(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	platforms := m.GetSupportedPlatforms()
	assert.Len(t, platforms, 7)
}

func TestMultiPlatformQualityModeler_GetPlatformProfile(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}
	for _, p := range platforms {
		profile := m.GetPlatformProfile(p)
		require.NotNil(t, profile, "profile for %s should not be nil", p)
		assert.Equal(t, p, profile.Platform)
		assert.NotEmpty(t, profile.Technology)
	}
}

func TestMultiPlatformQualityModeler_GetPlatformProfile_Unknown(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	profile := m.GetPlatformProfile(genomics.SequencingPlatform("UnknownPlatform"))
	assert.Nil(t, profile)
}

func TestMultiPlatformQualityModeler_GetPopulationData(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()

	pops := []population.Population{
		population.YRI, population.LWK, population.GWD,
		population.CEU, population.GBR, population.TSI,
		population.CHB, population.JPT, population.CHS,
		population.MXL, population.PUR,
	}
	for _, pop := range pops {
		data := m.GetPopulationData(pop)
		require.NotNil(t, data, "population data for %s should not be nil", pop)
		assert.Equal(t, pop, data.Population)
		assert.Greater(t, data.GCContentMean, 0.0)
		assert.Less(t, data.GCContentMean, 1.0)
		assert.Greater(t, data.GCContentStd, 0.0)
		assert.Greater(t, data.RepeatContent, 0.0)
	}
}

func TestMultiPlatformQualityModeler_GetPopulationData_Unknown(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	data := m.GetPopulationData(population.Population("UNKNOWN_POP"))
	assert.Nil(t, data)
}

func TestMultiPlatformQualityModeler_GenerateQualityScores_AllPlatforms(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	sequence := "ATCGATCGATCGATCGATCGATCGATCGATCG" // 32 bases

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			qualities := m.GenerateQualityScores(platform, population.CEU, sequence, 0)
			require.Len(t, qualities, len(sequence))

			profile := m.GetPlatformProfile(platform)
			require.NotNil(t, profile)

			for _, q := range qualities {
				assert.GreaterOrEqual(t, q, profile.QualityScoreDistribution.MinQuality)
				assert.LessOrEqual(t, q, profile.QualityScoreDistribution.MaxQuality)
			}
		})
	}
}

func TestMultiPlatformQualityModeler_GenerateQualityScores_UnknownPlatform(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	qualities := m.GenerateQualityScores(
		genomics.SequencingPlatform("UnknownPlatform"),
		population.CEU,
		"ATCGATCG",
		0,
	)
	// Falls back to default quality scores
	require.Len(t, qualities, 8)
	for _, q := range qualities {
		assert.GreaterOrEqual(t, q, 2)
		assert.LessOrEqual(t, q, 40)
	}
}

func TestMultiPlatformQualityModeler_GenerateQualityScores_AllPopulations(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	sequence := "GCGCGCGCGCGCGCGCGCGC" // GC-rich

	pops := []population.Population{
		population.YRI, population.LWK, population.GWD,
		population.CEU, population.GBR, population.TSI,
		population.CHB, population.JPT, population.CHS,
		population.MXL, population.PUR,
	}

	for _, pop := range pops {
		t.Run(string(pop), func(t *testing.T) {
			qualities := m.GenerateQualityScores(genomics.NovaSeq6000, pop, sequence, 0)
			assert.Len(t, qualities, len(sequence))
		})
	}
}

func TestMultiPlatformQualityModeler_GenerateQualityScores_ReadPosition(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	sequence := "ATCGATCGATCG"

	// Different read positions should produce different quality distributions
	q0 := m.GenerateQualityScores(genomics.NovaSeq6000, population.CEU, sequence, 0)
	q100 := m.GenerateQualityScores(genomics.NovaSeq6000, population.CEU, sequence, 100)

	assert.Len(t, q0, len(sequence))
	assert.Len(t, q100, len(sequence))
}

func TestMultiPlatformQualityModeler_GenerateQualityScores_HomopolymerSequence(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	// Homopolymer run of A
	sequence := "AAAAAAAAAAAAAAAAAAAAATCGATCGATCG"
	qualities := m.GenerateQualityScores(genomics.NovaSeq6000, population.CEU, sequence, 0)
	assert.Len(t, qualities, len(sequence))
}

func TestMultiPlatformQualityModeler_GenerateQualityScores_ShortSequence(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	qualities := m.GenerateQualityScores(genomics.NovaSeq6000, population.CEU, "A", 0)
	assert.Len(t, qualities, 1)
}

func TestMultiPlatformQualityModeler_ValidateMultiPlatformQuality(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	validation := m.ValidateMultiPlatformQuality()

	assert.True(t, validation["all_platforms_configured"])
	assert.True(t, validation["population_data_loaded"])
	assert.True(t, validation["quality_generation_functional"])
	assert.True(t, validation["illumina_platforms_configured"])
	assert.True(t, validation["long_read_platforms_configured"])
}

// ---------------------------------------------------------------------------
// QualityProfile fields
// ---------------------------------------------------------------------------

func TestQualityProfile_Illumina_Characteristics(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()

	illuminaPlatforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
	}

	for _, p := range illuminaPlatforms {
		profile := m.GetPlatformProfile(p)
		require.NotNil(t, profile)
		assert.Equal(t, "Illumina", profile.Technology)
		assert.Less(t, profile.ErrorProfile.BaseErrorRate, 0.002)
		assert.Less(t, profile.ReadLength.MeanLength, 10000) // Short reads
	}
}

func TestQualityProfile_LongRead_Characteristics(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()

	longReadPlatforms := []genomics.SequencingPlatform{
		genomics.PacBioRevio,
		genomics.PromethION24,
	}

	for _, p := range longReadPlatforms {
		profile := m.GetPlatformProfile(p)
		require.NotNil(t, profile)
		assert.Greater(t, profile.ReadLength.MeanLength, 10000)
		assert.Greater(t, profile.ErrorProfile.BaseErrorRate, 0.001)
		assert.Greater(t, profile.ErrorProfile.IndelErrorRate, 0.0)
	}
}

func TestQualityProfile_ReadLength_Spec(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	profile := m.GetPlatformProfile(genomics.NovaSeq6000)
	require.NotNil(t, profile)

	rl := profile.ReadLength
	assert.LessOrEqual(t, rl.MinLength, rl.MaxLength)
	assert.LessOrEqual(t, rl.MinLength, rl.MeanLength)
	assert.LessOrEqual(t, rl.MeanLength, rl.MaxLength)
	assert.NotEmpty(t, rl.Distribution)
}

func TestQualityProfile_GCBias_Characteristics(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()

	gcBiasPlatforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
	}

	for _, p := range gcBiasPlatforms {
		profile := m.GetPlatformProfile(p)
		require.NotNil(t, profile)
		assert.Greater(t, profile.GCBiasCharacteristics.BiasStrength, 0.0)
		assert.Greater(t, profile.GCBiasCharacteristics.OptimalGC, 0.0)
		assert.Less(t, profile.GCBiasCharacteristics.OptimalGC, 1.0)
	}
}

func TestQualityProfile_HomopolymerHandling(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()

	// These platforms should have homopolymer profiles
	platformsWithHomopolymer := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.PacBioRevio,
		genomics.PromethION24,
	}

	for _, p := range platformsWithHomopolymer {
		profile := m.GetPlatformProfile(p)
		require.NotNil(t, profile)
		assert.Greater(t, profile.HomopolymerHandling.MaxRunLength, 0, "platform %s", p)
	}
}

// ---------------------------------------------------------------------------
// PlatformErrorSignatureManager
// ---------------------------------------------------------------------------

func TestNewPlatformErrorSignatureManager(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	require.NotNil(t, m)
}

func TestPlatformErrorSignatureManager_GetSupportedPlatforms(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	platforms := m.GetSupportedPlatforms()
	assert.Len(t, platforms, 7)
}

func TestPlatformErrorSignatureManager_GetErrorSignature_AllPlatforms(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	expectedPlatforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}

	for _, p := range expectedPlatforms {
		sig := m.GetErrorSignature(p)
		require.NotNil(t, sig, "signature for %s should not be nil", p)
		assert.Equal(t, p, sig.Platform)
		assert.NotEmpty(t, sig.SignatureFingerprint)
		assert.NotEmpty(t, sig.ErrorPatterns)
	}
}

func TestPlatformErrorSignatureManager_GetErrorSignature_Unknown(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	sig := m.GetErrorSignature(genomics.SequencingPlatform("UnknownPlatform"))
	assert.Nil(t, sig)
}

func TestPlatformErrorSignatureManager_GetPopulationErrorData_AllPopulations(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	pops := []population.Population{
		population.YRI, population.LWK, population.GWD,
		population.CEU, population.GBR, population.TSI,
		population.CHB, population.JPT, population.CHS,
		population.MXL, population.PUR,
	}

	for _, pop := range pops {
		data := m.GetPopulationErrorData(pop)
		require.NotNil(t, data, "error data for %s should not be nil", pop)
		assert.Equal(t, pop, data.Population)
		assert.Greater(t, data.BaseErrorModifier, 0.0)
		assert.Greater(t, data.GCBiasModifier, 0.0)
	}
}

func TestPlatformErrorSignatureManager_GetPopulationErrorData_Unknown(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	data := m.GetPopulationErrorData(population.Population("UNKNOWN_POP"))
	assert.Nil(t, data)
}

func TestPlatformErrorSignatureManager_ValidatePlatformSignatures(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	validation := m.ValidatePlatformSignatures()

	assert.True(t, validation["all_platform_signatures_configured"])
	assert.True(t, validation["population_error_data_loaded"])
	assert.True(t, validation["error_application_functional"])
	assert.True(t, validation["signature_fingerprints_unique"])
	assert.True(t, validation["illumina_platforms_represented"])
	assert.True(t, validation["long_read_platforms_represented"])
	assert.True(t, validation["illumina_error_rates_appropriate"])
	assert.True(t, validation["long_read_indel_patterns"])
}

func TestPlatformErrorSignatureManager_ApplyErrorSignature_AllPlatforms(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	sequence := "ATCGATCGATCGATCGATCGATCGATCGATCG" // 32 bases
	qualities := make([]int, len(sequence))
	for i := range qualities {
		qualities[i] = 30
	}

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			modSeq, modQual, result := m.ApplyErrorSignature(
				platform, population.CEU, sequence, qualities, 0,
			)
			assert.True(t, result.Success, "platform %s: %s", platform, result.Message)
			assert.Equal(t, platform, result.Platform)
			assert.Equal(t, population.CEU, result.Population)
			assert.Len(t, modSeq, len(sequence))
			assert.Len(t, modQual, len(qualities))
			assert.GreaterOrEqual(t, result.ErrorsApplied, 0)
		})
	}
}

func TestPlatformErrorSignatureManager_ApplyErrorSignature_UnknownPlatform(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	seq := "ATCG"
	quals := []int{30, 30, 30, 30}

	retSeq, retQuals, result := m.ApplyErrorSignature(
		genomics.SequencingPlatform("UnknownPlatform"),
		population.CEU,
		seq,
		quals,
		0,
	)

	assert.False(t, result.Success)
	assert.Equal(t, seq, retSeq)
	assert.Equal(t, quals, retQuals)
}

func TestPlatformErrorSignatureManager_ApplyErrorSignature_UnknownPopulation(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	sequence := "ATCGATCGATCGATCG"
	qualities := make([]int, len(sequence))
	for i := range qualities {
		qualities[i] = 25
	}

	// Unknown population – should use default PopulationErrorData
	modSeq, modQual, result := m.ApplyErrorSignature(
		genomics.NovaSeq6000,
		population.Population("UNKNOWN"),
		sequence,
		qualities,
		0,
	)

	assert.True(t, result.Success)
	assert.Len(t, modSeq, len(sequence))
	assert.Len(t, modQual, len(qualities))
}

func TestPlatformErrorSignatureManager_ApplyErrorSignature_ShortSequence(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	modSeq, modQual, result := m.ApplyErrorSignature(
		genomics.NovaSeq6000,
		population.CEU,
		"A",
		[]int{30},
		0,
	)
	assert.True(t, result.Success)
	assert.NotEmpty(t, modSeq)
	assert.Len(t, modQual, 1)
}

func TestPlatformErrorSignatureManager_ApplyErrorSignature_AllPopulations(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	sequence := "GCGCGCGCGCGCGCGCGCGC"
	qualities := make([]int, len(sequence))
	for i := range qualities {
		qualities[i] = 35
	}

	pops := []population.Population{
		population.YRI, population.LWK, population.GWD,
		population.CEU, population.GBR, population.TSI,
		population.CHB, population.JPT, population.CHS,
		population.MXL, population.PUR,
	}

	for _, pop := range pops {
		t.Run(string(pop), func(t *testing.T) {
			_, _, result := m.ApplyErrorSignature(
				genomics.MiSeq, pop, sequence, qualities, 0,
			)
			assert.True(t, result.Success)
		})
	}
}

// ---------------------------------------------------------------------------
// Error pattern helper functions
// ---------------------------------------------------------------------------

func TestPlatformErrorSignature_LinearDecayPattern(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	// NovaSeq should have linear decay in substitution pattern
	sig := m.GetErrorSignature(genomics.NovaSeq6000)
	require.NotNil(t, sig)
	subPattern, exists := sig.ErrorPatterns["substitution"]
	assert.True(t, exists)
	assert.Len(t, subPattern.FrequencyPattern, 150)

	// Verify it's a linear decay (start < end)
	first := subPattern.FrequencyPattern[0]
	last := subPattern.FrequencyPattern[149]
	assert.NotEqual(t, first, last) // Should differ
}

func TestPlatformErrorSignature_ConstantRatePattern(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	// PacBio should have constant rates
	sig := m.GetErrorSignature(genomics.PacBioRevio)
	require.NotNil(t, sig)
	subPattern, exists := sig.ErrorPatterns["substitution"]
	assert.True(t, exists)
	assert.Len(t, subPattern.FrequencyPattern, 25000)

	// All values should be equal (constant rate)
	for _, rate := range subPattern.FrequencyPattern {
		assert.InDelta(t, subPattern.FrequencyPattern[0], rate, 1e-10)
	}
}

func TestPlatformErrorSignature_LongReadIndels(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	for _, p := range []genomics.SequencingPlatform{genomics.PacBioRevio, genomics.PromethION24} {
		sig := m.GetErrorSignature(p)
		require.NotNil(t, sig)

		_, hasIns := sig.ErrorPatterns["insertion"]
		assert.True(t, hasIns, "platform %s should have insertion patterns", p)

		_, hasDel := sig.ErrorPatterns["deletion"]
		assert.True(t, hasDel, "platform %s should have deletion patterns", p)
	}
}

func TestPlatformErrorSignature_QualityFingerprints(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.PromethION24,
	}

	for _, p := range platforms {
		sig := m.GetErrorSignature(p)
		require.NotNil(t, sig)
		assert.NotEmpty(t, sig.QualityFingerprints.ScoreDistribution, "platform %s", p)

		// Quality scores should be probabilities that sum to approximately 1
		total := 0.0
		for _, prob := range sig.QualityFingerprints.ScoreDistribution {
			assert.GreaterOrEqual(t, prob, 0.0)
			total += prob
		}
		assert.InDelta(t, 1.0, total, 0.01, "quality distribution for %s should sum to 1", p)
	}
}

func TestPlatformErrorSignature_SystematicBiases(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	for _, p := range []genomics.SequencingPlatform{
		genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq,
	} {
		sig := m.GetErrorSignature(p)
		require.NotNil(t, sig)
		// Illumina should have minimal strand bias
		assert.Less(t, sig.SystematicBiases.StrandBias, 0.05, "platform %s strand bias", p)
		// Reference bias around 0.5
		assert.Greater(t, sig.SystematicBiases.ReferenceBias, 0.48)
		assert.Less(t, sig.SystematicBiases.ReferenceBias, 0.60)
	}

	// Long read platforms should have zero duplication bias
	for _, p := range []genomics.SequencingPlatform{genomics.PacBioRevio, genomics.PromethION24} {
		sig := m.GetErrorSignature(p)
		require.NotNil(t, sig)
		assert.Equal(t, 0.0, sig.SystematicBiases.DuplicationBias, "platform %s", p)
	}
}

func TestPlatformErrorSignature_GCBias_Illumina(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	sig := m.GetErrorSignature(genomics.NovaSeq6000)
	require.NotNil(t, sig)

	gcBias := sig.QualityFingerprints.GCBiasSignature
	assert.Greater(t, gcBias.WindowSize, 0)
	assert.Len(t, gcBias.BiasFunction, 101) // 0-100%
}

func TestPlatformErrorSignature_HomopolymerSignature_LongRead(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	for _, p := range []genomics.SequencingPlatform{genomics.PacBioRevio, genomics.PromethION24} {
		sig := m.GetErrorSignature(p)
		require.NotNil(t, sig)
		assert.NotEmpty(t, sig.QualityFingerprints.HomopolymerSignature.RunLengthErrors)
		assert.NotEmpty(t, sig.QualityFingerprints.HomopolymerSignature.BaseSpecificBias)
	}
}

func TestPlatformErrorSignature_ContextDependency(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	sig := m.GetErrorSignature(genomics.NovaSeq6000)
	require.NotNil(t, sig)

	subPattern := sig.ErrorPatterns["substitution"]
	assert.NotEmpty(t, subPattern.ContextDependency)
	// GC bias context should be present for Illumina
	assert.Contains(t, subPattern.ContextDependency, "GG")
	assert.Contains(t, subPattern.ContextDependency, "CC")
}

// ---------------------------------------------------------------------------
// Quality score model correctness
// ---------------------------------------------------------------------------

func TestMultiPlatformQualityModeler_QualityBoundaries(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()

	// Long homopolymer sequence to test boundary clamping
	sequence := "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" // 58 A's

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.PacBioRevio,
		genomics.PromethION24,
	}

	for _, p := range platforms {
		profile := m.GetPlatformProfile(p)
		require.NotNil(t, profile)
		qualities := m.GenerateQualityScores(p, population.CEU, sequence, 0)
		for i, q := range qualities {
			assert.GreaterOrEqual(t, q, profile.QualityScoreDistribution.MinQuality,
				"position %d, platform %s: quality below minimum", i, p)
			assert.LessOrEqual(t, q, profile.QualityScoreDistribution.MaxQuality,
				"position %d, platform %s: quality above maximum", i, p)
		}
	}
}

func TestMultiPlatformQualityModeler_GCRichSequence(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	// 100% GC sequence
	sequence := "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" // 40 bases
	qualities := m.GenerateQualityScores(genomics.NovaSeq6000, population.CEU, sequence, 0)
	assert.Len(t, qualities, len(sequence))
}

func TestMultiPlatformQualityModeler_ATRichSequence(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	// 100% AT sequence
	sequence := "ATATATATATATATATATATATATATATATATATATATAT"
	qualities := m.GenerateQualityScores(genomics.NovaSeq6000, population.CEU, sequence, 0)
	assert.Len(t, qualities, len(sequence))
}

func TestMultiPlatformQualityModeler_EmptySequence(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	qualities := m.GenerateQualityScores(genomics.NovaSeq6000, population.CEU, "", 0)
	assert.Empty(t, qualities)
}

// ---------------------------------------------------------------------------
// Population-specific quality differences
// ---------------------------------------------------------------------------

func TestMultiPlatformQualityModeler_PopulationGCContent_AfricanLower(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()

	africanPops := []population.Population{population.YRI, population.LWK, population.GWD}
	europeanPops := []population.Population{population.CEU, population.GBR, population.TSI}

	for _, pop := range africanPops {
		data := m.GetPopulationData(pop)
		require.NotNil(t, data)
		// African populations should have slightly lower GC content
		for _, eur := range europeanPops {
			eurData := m.GetPopulationData(eur)
			require.NotNil(t, eurData)
			assert.LessOrEqual(t, data.GCContentMean, eurData.GCContentMean+0.01)
		}
	}
}

func TestMultiPlatformQualityModeler_PopulationRepeatContent_AfricanLower(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()

	africanData := m.GetPopulationData(population.YRI)
	require.NotNil(t, africanData)

	europeanData := m.GetPopulationData(population.CEU)
	require.NotNil(t, europeanData)

	// African populations have lower repeat content
	assert.Less(t, africanData.RepeatContent, europeanData.RepeatContent)
}

// ---------------------------------------------------------------------------
// ErrorApplicationResult
// ---------------------------------------------------------------------------

func TestErrorApplicationResult_Fields(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()
	_, _, result := m.ApplyErrorSignature(
		genomics.NovaSeq6000,
		population.CEU,
		"ATCGATCGATCG",
		[]int{30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
		5,
	)

	assert.Equal(t, genomics.NovaSeq6000, result.Platform)
	assert.Equal(t, population.CEU, result.Population)
	assert.True(t, result.Success)
	assert.NotEmpty(t, result.Message)
	assert.GreaterOrEqual(t, result.ErrorsApplied, 0)
}

// ---------------------------------------------------------------------------
// GC bias function correctness
// ---------------------------------------------------------------------------

func TestPlatformErrorSignature_GCBias_Quadratic(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	sig := m.GetErrorSignature(genomics.NovaSeq6000)
	require.NotNil(t, sig)

	biasFunc := sig.QualityFingerprints.GCBiasSignature.BiasFunction
	require.Len(t, biasFunc, 101)

	// At 50% GC (index 50), bias should be maximal (closest to 1.0)
	// At 0% or 100% GC, bias should be lower
	assert.Greater(t, biasFunc[50], biasFunc[0])
	assert.Greater(t, biasFunc[50], biasFunc[100])
}

func TestPlatformErrorSignature_GCBias_AllValuesPositive(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	sig := m.GetErrorSignature(genomics.NovaSeq6000)
	require.NotNil(t, sig)

	for i, v := range sig.QualityFingerprints.GCBiasSignature.BiasFunction {
		assert.Greater(t, v, 0.0, "bias[%d] should be positive", i)
	}
}

// ---------------------------------------------------------------------------
// Population error data modifiers
// ---------------------------------------------------------------------------

func TestPopulationErrorData_AfricanModifier(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	for _, pop := range []population.Population{population.YRI, population.LWK, population.GWD} {
		data := m.GetPopulationErrorData(pop)
		require.NotNil(t, data)
		// African populations have slightly lower error rate modifier
		assert.InDelta(t, 0.98, data.BaseErrorModifier, 0.01)
	}
}

func TestPopulationErrorData_NonAfricanModifier(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	for _, pop := range []population.Population{population.CEU, population.CHB} {
		data := m.GetPopulationErrorData(pop)
		require.NotNil(t, data)
		assert.InDelta(t, 1.0, data.BaseErrorModifier, 0.01)
	}
}

func TestPopulationErrorData_GCBiasModifier(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	// East Asian should have higher GC bias modifier
	chbData := m.GetPopulationErrorData(population.CHB)
	require.NotNil(t, chbData)
	assert.InDelta(t, 1.05, chbData.GCBiasModifier, 0.01)

	// European should be 1.0
	ceuData := m.GetPopulationErrorData(population.CEU)
	require.NotNil(t, ceuData)
	assert.InDelta(t, 1.0, ceuData.GCBiasModifier, 0.01)
}

func TestPopulationErrorData_RepeatTolerance(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	africanData := m.GetPopulationErrorData(population.YRI)
	require.NotNil(t, africanData)
	assert.InDelta(t, 1.10, africanData.RepeatTolerance, 0.01)

	europeanData := m.GetPopulationErrorData(population.CEU)
	require.NotNil(t, europeanData)
	assert.InDelta(t, 1.0, europeanData.RepeatTolerance, 0.01)
}

func TestPopulationErrorData_QualityBias(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	yriData := m.GetPopulationErrorData(population.YRI)
	require.NotNil(t, yriData)
	assert.InDelta(t, 1.02, yriData.QualityBias, 0.01)

	chbData := m.GetPopulationErrorData(population.CHB)
	require.NotNil(t, chbData)
	assert.InDelta(t, 0.99, chbData.QualityBias, 0.01)
}

// ---------------------------------------------------------------------------
// MGI / Element platform specifics
// ---------------------------------------------------------------------------

func TestPlatformErrorSignature_MGI_DNBSeq(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	sig := m.GetErrorSignature(genomics.DNBSeqT20x2)
	require.NotNil(t, sig)
	assert.Equal(t, "MGI_DNB_SEQUENCING_V1", sig.SignatureFingerprint)

	// Should have higher GC bias than Illumina
	gcBiasNova := m.GetErrorSignature(genomics.NovaSeq6000).ErrorPatterns["substitution"].ContextDependency["GG"]
	gcBiasDNB := sig.ErrorPatterns["substitution"].ContextDependency["GG"]
	assert.Greater(t, gcBiasDNB, gcBiasNova)
}

func TestPlatformErrorSignature_ElementAVITI(t *testing.T) {
	m := genomics.NewPlatformErrorSignatureManager()

	sig := m.GetErrorSignature(genomics.ElementAVITI)
	require.NotNil(t, sig)
	assert.Equal(t, "ELEMENT_AVIDITY_CHEMISTRY", sig.SignatureFingerprint)

	// Element AVITI has CpG bias
	subPattern := sig.ErrorPatterns["substitution"]
	assert.Contains(t, subPattern.ContextDependency, "CpG")
}

// ---------------------------------------------------------------------------
// Quality score cycle decay
// ---------------------------------------------------------------------------

func TestMultiPlatformQualityModeler_CycleDecay_NovaSeq(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	profile := m.GetPlatformProfile(genomics.NovaSeq6000)
	require.NotNil(t, profile)

	cycleDecay := profile.QualityScoreDistribution.CycleDecay
	assert.Greater(t, cycleDecay.DecayRate, 0.0)
	assert.Greater(t, cycleDecay.MinQuality, 0)
	assert.Greater(t, cycleDecay.PairEndDiff, 0.0)
}

// ---------------------------------------------------------------------------
// SubstitutionBias table
// ---------------------------------------------------------------------------

func TestQualityProfile_SubstitutionBias_SumsToOne(t *testing.T) {
	m := genomics.NewMultiPlatformQualityModeler()
	profile := m.GetPlatformProfile(genomics.NovaSeq6000)
	require.NotNil(t, profile)

	for fromBase, transitions := range profile.ErrorProfile.SubstitutionBias {
		sum := 0.0
		for _, prob := range transitions {
			sum += prob
		}
		assert.InDelta(t, 1.0, sum, 0.01, "substitution probabilities from %s should sum to 1", fromBase)
	}
}

// ---------------------------------------------------------------------------
// CrossPlatformBiasCorrector
// ---------------------------------------------------------------------------

func TestNewCrossPlatformBiasCorrector(t *testing.T) {
	corrector := genomics.NewCrossPlatformBiasCorrector()
	require.NotNil(t, corrector)
}

func TestNewBEADSNormalizer(t *testing.T) {
	n := genomics.NewBEADSNormalizer()
	require.NotNil(t, n)
}

func TestNewMANCIEIntegrator(t *testing.T) {
	n := genomics.NewMANCIEIntegrator()
	require.NotNil(t, n)
}

func TestNewQuantileNormalizer(t *testing.T) {
	n := genomics.NewQuantileNormalizer()
	require.NotNil(t, n)
}

func TestCrossPlatformBiasCorrector_CorrectCrossPlatformBias(t *testing.T) {
	corrector := genomics.NewCrossPlatformBiasCorrector()
	data := map[genomics.SequencingPlatform][]float64{
		genomics.NovaSeq6000: {1.0, 2.0, 3.0, 4.0, 5.0},
		genomics.MiSeq:       {1.1, 2.1, 3.1, 4.1, 5.1},
		genomics.PacBioRevio: {0.9, 1.9, 2.9, 3.9, 4.9},
	}

	result, err := corrector.CorrectCrossPlatformBias(data, population.CEU)
	require.NoError(t, err)
	require.NotNil(t, result)

	assert.Equal(t, data, result.OriginalData)
	assert.NotNil(t, result.CorrectedData)
	assert.Contains(t, result.CorrectedData, genomics.NovaSeq6000)
	assert.GreaterOrEqual(t, result.OverallQuality, 0.0)
	assert.NotNil(t, result.CorrectionMetrics)
}

func TestCrossPlatformBiasCorrector_EmptyData(t *testing.T) {
	corrector := genomics.NewCrossPlatformBiasCorrector()
	data := map[genomics.SequencingPlatform][]float64{}

	result, err := corrector.CorrectCrossPlatformBias(data, population.YRI)
	require.NoError(t, err)
	assert.NotNil(t, result)
	assert.Empty(t, result.CorrectedData)
}

func TestCrossPlatformBiasCorrector_SinglePlatform(t *testing.T) {
	corrector := genomics.NewCrossPlatformBiasCorrector()
	data := map[genomics.SequencingPlatform][]float64{
		genomics.NovaSeq6000: {10.0, 20.0, 30.0},
	}

	result, err := corrector.CorrectCrossPlatformBias(data, population.CHB)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Len(t, result.CorrectedData, 1)
}

func TestBEADSNormalizer_ApplyBEADSCorrection_AllPopulations(t *testing.T) {
	n := genomics.NewBEADSNormalizer()
	data := map[genomics.SequencingPlatform][]float64{
		genomics.NovaSeq6000: {1.0, 2.0, 3.0},
		genomics.NextSeq2000: {1.5, 2.5, 3.5},
	}

	pops := []population.Population{
		population.YRI, population.CEU, population.CHB,
		population.MXL, population.GBR,
	}

	for _, pop := range pops {
		result, err := n.ApplyBEADSCorrection(data, pop)
		require.NoError(t, err, "population %s", pop)
		assert.NotNil(t, result)
		assert.Len(t, result.CorrectedData, 2)
		assert.Contains(t, result.BiasFactors, genomics.NovaSeq6000)
	}
}

func TestMANCIEIntegrator_ApplyMANCIEIntegration(t *testing.T) {
	m := genomics.NewMANCIEIntegrator()
	data := map[genomics.SequencingPlatform][]float64{
		genomics.NovaSeq6000: {1.0, 2.0, 3.0, 4.0},
		genomics.MiSeq:       {1.1, 2.1, 3.1, 4.1},
		genomics.PacBioRevio: {0.9, 1.9, 2.9, 3.9},
	}

	result, err := m.ApplyMANCIEIntegration(data)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.NotNil(t, result.IntegratedData)
	assert.Len(t, result.IntegratedData, 3)
	assert.GreaterOrEqual(t, result.TransformationQuality, 0.0)
}

func TestMANCIEIntegrator_ApplyMANCIEIntegration_EmptyData(t *testing.T) {
	m := genomics.NewMANCIEIntegrator()
	result, err := m.ApplyMANCIEIntegration(map[genomics.SequencingPlatform][]float64{})
	require.NoError(t, err)
	assert.NotNil(t, result)
	assert.Empty(t, result.IntegratedData)
}

func TestQuantileNormalizer_ApplyQuantileNormalization(t *testing.T) {
	q := genomics.NewQuantileNormalizer()
	data := map[genomics.SequencingPlatform][]float64{
		genomics.NovaSeq6000: {5.0, 3.0, 1.0, 4.0, 2.0},
		genomics.MiSeq:       {6.0, 4.0, 2.0, 5.0, 3.0},
	}

	result, err := q.ApplyQuantileNormalization(data)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Len(t, result.NormalizedData, 2)
	assert.NotEmpty(t, result.ReferenceDistribution)
	assert.GreaterOrEqual(t, result.NormalizationQuality, 0.0)
}

func TestQuantileNormalizer_ApplyQuantileNormalization_SingleElement(t *testing.T) {
	q := genomics.NewQuantileNormalizer()
	data := map[genomics.SequencingPlatform][]float64{
		genomics.NovaSeq6000: {1.0},
	}

	result, err := q.ApplyQuantileNormalization(data)
	require.NoError(t, err)
	assert.Len(t, result.NormalizedData, 1)
}

func TestCrossPlatformCorrectionResult_Fields(t *testing.T) {
	corrector := genomics.NewCrossPlatformBiasCorrector()
	data := map[genomics.SequencingPlatform][]float64{
		genomics.NovaSeq6000: {1.0, 2.0, 3.0},
	}

	result, err := corrector.CorrectCrossPlatformBias(data, population.CEU)
	require.NoError(t, err)

	assert.Equal(t, data, result.OriginalData)
	assert.Equal(t, 0.0, result.BiasReductionScore)    // default
	assert.Equal(t, 0.0, result.StatisticalIndistinguishability) // default
	assert.Greater(t, result.OverallQuality, 0.0)
}

func TestBEADSNormalizer_LargerDataset(t *testing.T) {
	n := genomics.NewBEADSNormalizer()
	// Create larger dataset to exercise more code paths
	largeData := make([]float64, 100)
	for i := range largeData {
		largeData[i] = float64(i) * 0.1
	}
	data := map[genomics.SequencingPlatform][]float64{
		genomics.NovaSeq6000: largeData,
		genomics.DNBSeqT20x2: largeData,
		genomics.ElementAVITI: largeData,
	}

	result, err := n.ApplyBEADSCorrection(data, population.JPT)
	require.NoError(t, err)
	assert.Len(t, result.CorrectedData, 3)
}

func TestCrossPlatformBiasCorrector_AllSevenPlatforms(t *testing.T) {
	corrector := genomics.NewCrossPlatformBiasCorrector()
	data := make(map[genomics.SequencingPlatform][]float64)
	for i, p := range []genomics.SequencingPlatform{
		genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq,
		genomics.PacBioRevio, genomics.DNBSeqT20x2, genomics.ElementAVITI,
		genomics.PromethION24,
	} {
		vals := make([]float64, 10)
		for j := range vals {
			vals[j] = float64(i*10 + j)
		}
		data[p] = vals
	}

	result, err := corrector.CorrectCrossPlatformBias(data, population.TSI)
	require.NoError(t, err)
	assert.Len(t, result.CorrectedData, 7)
	assert.Len(t, result.CorrectionMetrics, 7)
}
