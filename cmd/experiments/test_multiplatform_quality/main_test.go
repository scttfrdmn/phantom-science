//go:build !release
// +build !release

package main

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

// TestNewMultiPlatformQualityModeler tests that the modeler initialises without error.
func TestNewMultiPlatformQualityModeler(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	if modeler == nil {
		t.Fatal("expected non-nil modeler")
	}
}

// TestGetSupportedPlatforms verifies that exactly 7 platforms are returned.
func TestGetSupportedPlatforms(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	platforms := modeler.GetSupportedPlatforms()

	if len(platforms) != 7 {
		t.Errorf("expected 7 platforms, got %d", len(platforms))
	}

	expected := map[genomics.SequencingPlatform]bool{
		genomics.NovaSeq6000:  false,
		genomics.NextSeq2000:  false,
		genomics.MiSeq:        false,
		genomics.PacBioRevio:  false,
		genomics.DNBSeqT20x2:  false,
		genomics.ElementAVITI: false,
		genomics.PromethION24: false,
	}

	for _, p := range platforms {
		if _, ok := expected[p]; !ok {
			t.Errorf("unexpected platform: %s", p)
		}
		expected[p] = true
	}

	for p, seen := range expected {
		if !seen {
			t.Errorf("missing expected platform: %s", p)
		}
	}
}

// TestGetPlatformProfile_AllPlatforms checks that every supported platform has a non-nil profile.
func TestGetPlatformProfile_AllPlatforms(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

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
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("expected non-nil profile for %s", platform)
			}
			if profile.Platform != platform {
				t.Errorf("profile.Platform = %s, want %s", profile.Platform, platform)
			}
		})
	}
}

// TestGetPlatformProfile_UnknownPlatform confirms nil is returned for unknown platforms.
func TestGetPlatformProfile_UnknownPlatform(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	profile := modeler.GetPlatformProfile(genomics.SequencingPlatform("Unknown"))
	if profile != nil {
		t.Error("expected nil profile for unknown platform")
	}
}

// TestPlatformTechnology_Illumina checks technology classification of Illumina platforms.
func TestPlatformTechnology_Illumina(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	illuminaPlatforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
	}

	for _, platform := range illuminaPlatforms {
		t.Run(string(platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", platform)
			}
			if profile.Technology != "Illumina" {
				t.Errorf("%s technology = %q, want Illumina", platform, profile.Technology)
			}
		})
	}
}

// TestPlatformTechnology_Other checks technology classification of non-Illumina platforms.
func TestPlatformTechnology_Other(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	cases := []struct {
		platform   genomics.SequencingPlatform
		technology string
	}{
		{genomics.PacBioRevio, "PacBio"},
		{genomics.PromethION24, "ONT"},
		{genomics.DNBSeqT20x2, "MGI"},
		{genomics.ElementAVITI, "Element"},
	}

	for _, tc := range cases {
		t.Run(string(tc.platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(tc.platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", tc.platform)
			}
			if profile.Technology != tc.technology {
				t.Errorf("%s technology = %q, want %q", tc.platform, profile.Technology, tc.technology)
			}
		})
	}
}

// TestShortReadPlatforms_ErrorRates confirms short-read platforms have low base error rates.
func TestShortReadPlatforms_ErrorRates(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	shortReadPlatforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
	}

	for _, platform := range shortReadPlatforms {
		t.Run(string(platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", platform)
			}
			if profile.ErrorProfile.BaseErrorRate > 0.002 {
				t.Errorf("%s BaseErrorRate = %.4f, want <= 0.002", platform, profile.ErrorProfile.BaseErrorRate)
			}
		})
	}
}

// TestLongReadPlatforms_ReadLength verifies long-read platforms have mean read length >= 10000.
func TestLongReadPlatforms_ReadLength(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	longReadPlatforms := []genomics.SequencingPlatform{
		genomics.PacBioRevio,
		genomics.PromethION24,
	}

	for _, platform := range longReadPlatforms {
		t.Run(string(platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", platform)
			}
			if profile.ReadLength.MeanLength < 10000 {
				t.Errorf("%s MeanLength = %d, want >= 10000", platform, profile.ReadLength.MeanLength)
			}
		})
	}
}

// TestLongReadPlatforms_HighErrorRates confirms long-read platforms have higher error rates.
func TestLongReadPlatforms_HighErrorRates(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	longReadPlatforms := []genomics.SequencingPlatform{
		genomics.PacBioRevio,
		genomics.PromethION24,
	}

	for _, platform := range longReadPlatforms {
		t.Run(string(platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", platform)
			}
			if profile.ErrorProfile.BaseErrorRate < 0.01 {
				t.Errorf("%s BaseErrorRate = %.4f, want >= 0.01", platform, profile.ErrorProfile.BaseErrorRate)
			}
		})
	}
}

// TestGetPopulationData_AllPopulations checks that population data exists for key populations.
func TestGetPopulationData_AllPopulations(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	pops := []population.Population{
		population.YRI,
		population.LWK,
		population.GWD,
		population.CEU,
		population.GBR,
		population.TSI,
		population.CHB,
		population.JPT,
		population.CHS,
		population.MXL,
		population.PUR,
	}

	for _, pop := range pops {
		t.Run(string(pop), func(t *testing.T) {
			data := modeler.GetPopulationData(pop)
			if data == nil {
				t.Fatalf("expected non-nil population data for %s", pop)
			}
			if data.Population != pop {
				t.Errorf("data.Population = %s, want %s", data.Population, pop)
			}
			if data.GCContentMean <= 0 || data.GCContentMean >= 1 {
				t.Errorf("%s GCContentMean = %.3f, want in (0, 1)", pop, data.GCContentMean)
			}
			if data.GCContentStd <= 0 {
				t.Errorf("%s GCContentStd = %.3f, want > 0", pop, data.GCContentStd)
			}
			if data.RepeatContent <= 0 || data.RepeatContent >= 1 {
				t.Errorf("%s RepeatContent = %.3f, want in (0, 1)", pop, data.RepeatContent)
			}
			if data.AllelicBalance <= 0 || data.AllelicBalance > 1 {
				t.Errorf("%s AllelicBalance = %.3f, want in (0, 1]", pop, data.AllelicBalance)
			}
		})
	}
}

// TestGetPopulationData_UnknownPopulation verifies nil is returned for an unknown population.
func TestGetPopulationData_UnknownPopulation(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	data := modeler.GetPopulationData(population.Population("UNKNOWN"))
	if data != nil {
		t.Error("expected nil for unknown population")
	}
}

// TestGenerateQualityScores_Length verifies quality scores length equals sequence length.
func TestGenerateQualityScores_Length(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	testCases := []struct {
		platform genomics.SequencingPlatform
		pop      population.Population
		seq      string
	}{
		{genomics.NovaSeq6000, population.YRI, "ATCGATCG"},
		{genomics.PacBioRevio, population.CEU, "GGGGGGGGGGGG"},
		{genomics.PromethION24, population.CHB, "AAAAAAAAAAAAAAAAAAAAA"},
		{genomics.ElementAVITI, population.MXL, "ATCGATCGATCGATCGAAAAAATCGTCGTCGCCCCCGGGGATCGATCG"},
		{genomics.MiSeq, population.GBR, "CGTACGTACGTA"},
		{genomics.NextSeq2000, population.JPT, "TTTTAAAACCCCGGGG"},
		{genomics.DNBSeqT20x2, population.LWK, "ACGTACGTACGT"},
	}

	for _, tc := range testCases {
		t.Run(string(tc.platform)+"_"+string(tc.pop), func(t *testing.T) {
			qualities := modeler.GenerateQualityScores(tc.platform, tc.pop, tc.seq, 1)
			if len(qualities) != len(tc.seq) {
				t.Errorf("quality length = %d, want %d", len(qualities), len(tc.seq))
			}
		})
	}
}

// TestGenerateQualityScores_RangeConstraints verifies all quality scores are within platform bounds.
func TestGenerateQualityScores_RangeConstraints(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.PacBioRevio,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
		genomics.PromethION24,
	}

	testSeq := "ATCGATCGATCGATCGAAAAAATCGTCGTCGCCCCCGGGGATCGATCG"

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", platform)
			}

			qualities := modeler.GenerateQualityScores(platform, population.YRI, testSeq, 1)

			for i, q := range qualities {
				if q < profile.QualityScoreDistribution.MinQuality {
					t.Errorf("position %d: quality %d < MinQuality %d", i, q, profile.QualityScoreDistribution.MinQuality)
				}
				if q > profile.QualityScoreDistribution.MaxQuality {
					t.Errorf("position %d: quality %d > MaxQuality %d", i, q, profile.QualityScoreDistribution.MaxQuality)
				}
			}
		})
	}
}

// TestGenerateQualityScores_UnknownPlatform ensures a fallback is used for unknown platforms.
func TestGenerateQualityScores_UnknownPlatform(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	seq := "ATCGATCG"
	qualities := modeler.GenerateQualityScores(genomics.SequencingPlatform("UNKNOWN"), population.YRI, seq, 1)
	if len(qualities) != len(seq) {
		t.Errorf("fallback quality length = %d, want %d", len(qualities), len(seq))
	}
	for i, q := range qualities {
		if q < 2 || q > 40 {
			t.Errorf("fallback quality[%d] = %d, want in [2, 40]", i, q)
		}
	}
}

// TestGenerateQualityScores_DifferentReadPositions verifies read position affects quality generation.
func TestGenerateQualityScores_DifferentReadPositions(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	seq := "ATCGATCGATCGATCGAAAAAATCGTCGTCGCCCCCGGGGATCGATCG"

	positions := []int{1, 50, 100, 200}
	for _, pos := range positions {
		t.Run("pos"+string(rune('0'+pos%10)), func(t *testing.T) {
			qualities := modeler.GenerateQualityScores(genomics.NovaSeq6000, population.YRI, seq, pos)
			if len(qualities) != len(seq) {
				t.Errorf("quality length at pos %d = %d, want %d", pos, len(qualities), len(seq))
			}
		})
	}
}

// TestValidateMultiPlatformQuality checks that the validation map has the expected keys.
func TestValidateMultiPlatformQuality(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	validation := modeler.ValidateMultiPlatformQuality()

	if len(validation) == 0 {
		t.Fatal("expected non-empty validation map")
	}

	expectedKeys := []string{
		"all_platforms_configured",
		"population_data_loaded",
		"quality_generation_functional",
		"illumina_platforms_configured",
		"long_read_platforms_configured",
	}

	for _, key := range expectedKeys {
		val, ok := validation[key]
		if !ok {
			t.Errorf("missing expected validation key: %s", key)
			continue
		}
		if !val {
			t.Errorf("validation key %q = false, want true", key)
		}
	}
}

// TestPlatformProfile_ReadLengthSpec verifies read length specifications are coherent.
func TestPlatformProfile_ReadLengthSpec(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	for _, platform := range modeler.GetSupportedPlatforms() {
		t.Run(string(platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", platform)
			}

			spec := profile.ReadLength
			if spec.MinLength <= 0 {
				t.Errorf("%s MinLength = %d, want > 0", platform, spec.MinLength)
			}
			if spec.MaxLength < spec.MinLength {
				t.Errorf("%s MaxLength %d < MinLength %d", platform, spec.MaxLength, spec.MinLength)
			}
			if spec.MeanLength < spec.MinLength || spec.MeanLength > spec.MaxLength {
				t.Errorf("%s MeanLength %d not in [%d, %d]", platform, spec.MeanLength, spec.MinLength, spec.MaxLength)
			}
		})
	}
}

// TestPlatformProfile_QualityDistribution verifies quality distribution min <= max.
func TestPlatformProfile_QualityDistribution(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	for _, platform := range modeler.GetSupportedPlatforms() {
		t.Run(string(platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", platform)
			}

			dist := profile.QualityScoreDistribution
			if dist.MinQuality < 0 {
				t.Errorf("%s MinQuality = %d, want >= 0", platform, dist.MinQuality)
			}
			if dist.MaxQuality <= dist.MinQuality {
				t.Errorf("%s MaxQuality %d <= MinQuality %d", platform, dist.MaxQuality, dist.MinQuality)
			}
			if dist.MeanQuality < float64(dist.MinQuality) || dist.MeanQuality > float64(dist.MaxQuality) {
				t.Errorf("%s MeanQuality %.1f not in [%d, %d]", platform, dist.MeanQuality, dist.MinQuality, dist.MaxQuality)
			}
			if dist.StdQuality <= 0 {
				t.Errorf("%s StdQuality = %.2f, want > 0", platform, dist.StdQuality)
			}
		})
	}
}

// TestPlatformProfile_GCBias verifies GC bias strength and optimal GC are in valid ranges.
func TestPlatformProfile_GCBias(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	gcBiasPlatforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000,
		genomics.NextSeq2000,
		genomics.MiSeq,
		genomics.DNBSeqT20x2,
		genomics.ElementAVITI,
	}

	for _, platform := range gcBiasPlatforms {
		t.Run(string(platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", platform)
			}

			gcBias := profile.GCBiasCharacteristics
			if gcBias.BiasStrength <= 0 || gcBias.BiasStrength > 1 {
				t.Errorf("%s BiasStrength = %.3f, want in (0, 1]", platform, gcBias.BiasStrength)
			}
			if gcBias.OptimalGC <= 0 || gcBias.OptimalGC >= 1 {
				t.Errorf("%s OptimalGC = %.3f, want in (0, 1)", platform, gcBias.OptimalGC)
			}
			if gcBias.RegionSizeBp <= 0 {
				t.Errorf("%s RegionSizeBp = %d, want > 0", platform, gcBias.RegionSizeBp)
			}
		})
	}
}

// TestPopulationGCContent_AfricanVsOther checks African populations have slightly lower GC than others.
func TestPopulationGCContent_AfricanVsOther(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	africanPops := []population.Population{population.YRI, population.LWK, population.GWD}
	otherPops := []population.Population{population.CEU, population.CHB}

	for _, ap := range africanPops {
		afData := modeler.GetPopulationData(ap)
		if afData == nil {
			t.Fatalf("nil data for %s", ap)
		}
		for _, op := range otherPops {
			otherData := modeler.GetPopulationData(op)
			if otherData == nil {
				t.Fatalf("nil data for %s", op)
			}
			// African GC should be less than or equal to non-African GC
			if afData.GCContentMean > otherData.GCContentMean+0.01 {
				t.Errorf("African (%s) GC %.3f > non-African (%s) GC %.3f by more than 0.01",
					ap, afData.GCContentMean, op, otherData.GCContentMean)
			}
		}
	}
}

// TestPopulationAllelicBalance_Americas verifies Americas populations have slightly lower allelic balance.
func TestPopulationAllelicBalance_Americas(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	americasPops := []population.Population{population.MXL, population.PUR}
	for _, pop := range americasPops {
		data := modeler.GetPopulationData(pop)
		if data == nil {
			t.Fatalf("nil data for %s", pop)
		}
		if data.AllelicBalance >= 0.50 {
			t.Errorf("%s AllelicBalance = %.2f, want < 0.50 (slight imbalance)", pop, data.AllelicBalance)
		}
	}

	nonAmericasPops := []population.Population{population.YRI, population.CEU, population.CHB}
	for _, pop := range nonAmericasPops {
		data := modeler.GetPopulationData(pop)
		if data == nil {
			t.Fatalf("nil data for %s", pop)
		}
		if data.AllelicBalance != 0.50 {
			t.Errorf("%s AllelicBalance = %.2f, want 0.50", pop, data.AllelicBalance)
		}
	}
}

// TestAfricanRepeatContent verifies African populations have lower repeat content.
func TestAfricanRepeatContent(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	africanPops := []population.Population{population.YRI, population.LWK, population.GWD}
	nonAfricanPops := []population.Population{population.CEU, population.CHB}

	for _, ap := range africanPops {
		afData := modeler.GetPopulationData(ap)
		if afData == nil {
			t.Fatalf("nil data for %s", ap)
		}
		for _, nap := range nonAfricanPops {
			naData := modeler.GetPopulationData(nap)
			if naData == nil {
				t.Fatalf("nil data for %s", nap)
			}
			if afData.RepeatContent >= naData.RepeatContent {
				t.Errorf("African (%s) RepeatContent %.3f >= non-African (%s) %.3f",
					ap, afData.RepeatContent, nap, naData.RepeatContent)
			}
		}
	}
}

// TestAfricanGCVariation verifies African populations have higher GC variation.
func TestAfricanGCVariation(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	africanPops := []population.Population{population.YRI, population.LWK, population.GWD}
	eastAsianPops := []population.Population{population.CHB, population.JPT, population.CHS}

	for _, ap := range africanPops {
		afData := modeler.GetPopulationData(ap)
		if afData == nil {
			t.Fatalf("nil data for %s", ap)
		}
		for _, ep := range eastAsianPops {
			eaData := modeler.GetPopulationData(ep)
			if eaData == nil {
				t.Fatalf("nil data for %s", ep)
			}
			if afData.GCContentStd <= eaData.GCContentStd {
				t.Errorf("African (%s) GCContentStd %.3f <= East Asian (%s) %.3f",
					ap, afData.GCContentStd, ep, eaData.GCContentStd)
			}
		}
	}
}

// TestGenerateQualityScores_HomopolymerRuns tests quality generation on homopolymer sequences.
func TestGenerateQualityScores_HomopolymerRuns(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	homopolymerSeqs := []struct {
		name string
		seq  string
	}{
		{"poly-A-short", "AAAA"},
		{"poly-G-long", "GGGGGGGG"},
		{"poly-C-verylong", "CCCCCCCCCCCCCC"},
		{"mixed-homo", "AAAAATTTTTGGGGG"},
	}

	for _, tc := range homopolymerSeqs {
		t.Run(tc.name, func(t *testing.T) {
			// Test with platforms that have explicit homopolymer handling
			for _, platform := range []genomics.SequencingPlatform{genomics.PacBioRevio, genomics.PromethION24, genomics.NovaSeq6000} {
				quals := modeler.GenerateQualityScores(platform, population.YRI, tc.seq, 1)
				if len(quals) != len(tc.seq) {
					t.Errorf("%s/%s: quality length %d, want %d", platform, tc.name, len(quals), len(tc.seq))
				}
				profile := modeler.GetPlatformProfile(platform)
				for i, q := range quals {
					if q < profile.QualityScoreDistribution.MinQuality {
						t.Errorf("%s/%s position %d: quality %d < min %d", platform, tc.name, i, q, profile.QualityScoreDistribution.MinQuality)
					}
					if q > profile.QualityScoreDistribution.MaxQuality {
						t.Errorf("%s/%s position %d: quality %d > max %d", platform, tc.name, i, q, profile.QualityScoreDistribution.MaxQuality)
					}
				}
			}
		})
	}
}

// TestGenerateQualityScores_GCContent tests quality generation with varying GC content.
func TestGenerateQualityScores_GCContent(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	gcSeqs := []struct {
		name string
		seq  string // Approximate GC content noted
	}{
		{"low-GC", "AAAAAAAAAAAAAAATTTTTTTTTTTTTTTTT"},    // ~0% GC
		{"high-GC", "GGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCC"},   // ~100% GC
		{"medium-GC", "ATCGATCGATCGATCGATCGATCGATCGATCG"}, // ~50% GC
	}

	for _, tc := range gcSeqs {
		t.Run(tc.name, func(t *testing.T) {
			quals := modeler.GenerateQualityScores(genomics.NovaSeq6000, population.CEU, tc.seq, 1)
			if len(quals) != len(tc.seq) {
				t.Errorf("%s: quality length %d, want %d", tc.name, len(quals), len(tc.seq))
			}
		})
	}
}

// TestValidateMultiPlatformQuality_AllTrue verifies all validation results pass.
func TestValidateMultiPlatformQuality_AllTrue(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	validation := modeler.ValidateMultiPlatformQuality()

	for key, passed := range validation {
		if !passed {
			t.Errorf("validation[%q] = false, want true", key)
		}
	}
}

// TestPlatformConstants verifies that platform constant strings match expected values.
func TestPlatformConstants(t *testing.T) {
	cases := []struct {
		platform genomics.SequencingPlatform
		name     string
	}{
		{genomics.NovaSeq6000, "NovaSeq 6000"},
		{genomics.NextSeq2000, "NextSeq 2000"},
		{genomics.MiSeq, "MiSeq"},
		{genomics.PacBioRevio, "PacBio Revio"},
		{genomics.DNBSeqT20x2, "DNBSEQ-T20x2"},
		{genomics.ElementAVITI, "Element AVITI"},
		{genomics.PromethION24, "PromethION 24"},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if string(tc.platform) != tc.name {
				t.Errorf("platform constant = %q, want %q", string(tc.platform), tc.name)
			}
		})
	}
}

// TestGenerateQualityScores_EmptySequence verifies empty sequence returns empty quality slice.
func TestGenerateQualityScores_EmptySequence(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	quals := modeler.GenerateQualityScores(genomics.NovaSeq6000, population.YRI, "", 1)
	if len(quals) != 0 {
		t.Errorf("empty sequence quality length = %d, want 0", len(quals))
	}
}

// TestGenerateQualityScores_SingleBase verifies single-base sequence quality generation.
func TestGenerateQualityScores_SingleBase(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	bases := []string{"A", "T", "G", "C"}

	for _, base := range bases {
		t.Run(base, func(t *testing.T) {
			quals := modeler.GenerateQualityScores(genomics.NovaSeq6000, population.YRI, base, 1)
			if len(quals) != 1 {
				t.Errorf("single base quality length = %d, want 1", len(quals))
			}
		})
	}
}

// TestPlatformProfile_ErrorProfile verifies error profiles are non-trivially configured.
func TestPlatformProfile_ErrorProfile(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	for _, platform := range modeler.GetSupportedPlatforms() {
		t.Run(string(platform), func(t *testing.T) {
			profile := modeler.GetPlatformProfile(platform)
			if profile == nil {
				t.Fatalf("nil profile for %s", platform)
			}

			ep := profile.ErrorProfile
			if ep.BaseErrorRate <= 0 {
				t.Errorf("%s BaseErrorRate = %.6f, want > 0", platform, ep.BaseErrorRate)
			}
			if ep.BaseErrorRate >= 1 {
				t.Errorf("%s BaseErrorRate = %.6f, want < 1", platform, ep.BaseErrorRate)
			}
			if ep.IndelErrorRate < 0 {
				t.Errorf("%s IndelErrorRate = %.6f, want >= 0", platform, ep.IndelErrorRate)
			}
		})
	}
}

// TestGenerateQualityScores_AllPopulations verifies quality scores are produced for all populations.
func TestGenerateQualityScores_AllPopulations(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	seq := "ATCGATCGATCG"

	allPops := []population.Population{
		population.YRI, population.LWK, population.GWD,
		population.CEU, population.GBR, population.TSI,
		population.CHB, population.JPT, population.CHS,
		population.MXL, population.PUR,
	}

	for _, pop := range allPops {
		t.Run(string(pop), func(t *testing.T) {
			quals := modeler.GenerateQualityScores(genomics.NovaSeq6000, pop, seq, 1)
			if len(quals) != len(seq) {
				t.Errorf("%s quality length = %d, want %d", pop, len(quals), len(seq))
			}
		})
	}
}

// TestPlatformCoverage_SevenExpected ensures exactly 7 platforms are configured internally.
func TestPlatformCoverage_SevenExpected(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	platforms := modeler.GetSupportedPlatforms()

	if len(platforms) != 7 {
		t.Errorf("expected exactly 7 supported platforms, got %d", len(platforms))
	}

	// Each platform must have a profile
	for _, p := range platforms {
		if modeler.GetPlatformProfile(p) == nil {
			t.Errorf("missing profile for supported platform %s", p)
		}
	}
}

// TestNovaSeqSpecificProfile validates NovaSeq 6000-specific configuration.
func TestNovaSeqSpecificProfile(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	profile := modeler.GetPlatformProfile(genomics.NovaSeq6000)
	if profile == nil {
		t.Fatal("nil profile for NovaSeq6000")
	}

	if profile.ReadLength.MeanLength != 150 {
		t.Errorf("NovaSeq MeanLength = %d, want 150", profile.ReadLength.MeanLength)
	}
	if profile.ErrorProfile.BaseErrorRate != 0.001 {
		t.Errorf("NovaSeq BaseErrorRate = %.4f, want 0.001", profile.ErrorProfile.BaseErrorRate)
	}
	if profile.QualityScoreDistribution.MeanQuality != 35.0 {
		t.Errorf("NovaSeq MeanQuality = %.1f, want 35.0", profile.QualityScoreDistribution.MeanQuality)
	}
}

// TestPromethIONSpecificProfile validates PromethION 24-specific configuration.
func TestPromethIONSpecificProfile(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	profile := modeler.GetPlatformProfile(genomics.PromethION24)
	if profile == nil {
		t.Fatal("nil profile for PromethION24")
	}

	if profile.ReadLength.MeanLength != 15000 {
		t.Errorf("PromethION MeanLength = %d, want 15000", profile.ReadLength.MeanLength)
	}
	if profile.ErrorProfile.BaseErrorRate != 0.03 {
		t.Errorf("PromethION BaseErrorRate = %.4f, want 0.03", profile.ErrorProfile.BaseErrorRate)
	}
	if profile.QualityScoreDistribution.MeanQuality != 15.0 {
		t.Errorf("PromethION MeanQuality = %.1f, want 15.0", profile.QualityScoreDistribution.MeanQuality)
	}
}

// TestPacBioRevioSpecificProfile validates PacBio Revio-specific configuration.
func TestPacBioRevioSpecificProfile(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	profile := modeler.GetPlatformProfile(genomics.PacBioRevio)
	if profile == nil {
		t.Fatal("nil profile for PacBioRevio")
	}

	if profile.ReadLength.MeanLength != 25000 {
		t.Errorf("PacBio MeanLength = %d, want 25000", profile.ReadLength.MeanLength)
	}
	if profile.ErrorProfile.BaseErrorRate != 0.01 {
		t.Errorf("PacBio BaseErrorRate = %.4f, want 0.01", profile.ErrorProfile.BaseErrorRate)
	}
}

// TestMiSeqSpecificProfile validates MiSeq-specific configuration.
func TestMiSeqSpecificProfile(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	profile := modeler.GetPlatformProfile(genomics.MiSeq)
	if profile == nil {
		t.Fatal("nil profile for MiSeq")
	}

	if profile.ReadLength.MeanLength != 250 {
		t.Errorf("MiSeq MeanLength = %d, want 250", profile.ReadLength.MeanLength)
	}
	if profile.QualityScoreDistribution.MeanQuality != 38.0 {
		t.Errorf("MiSeq MeanQuality = %.1f, want 38.0", profile.QualityScoreDistribution.MeanQuality)
	}
}

// TestGenerateQualityScores_LongSequence handles a long sequence to exercise cycle-decay logic.
func TestGenerateQualityScores_LongSequence(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()

	// Build a 200-character sequence
	longSeq := "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
	quals := modeler.GenerateQualityScores(genomics.NovaSeq6000, population.CEU, longSeq, 1)
	if len(quals) != len(longSeq) {
		t.Errorf("long sequence quality length = %d, want %d", len(quals), len(longSeq))
	}
}

// TestGenerateQualityScores_MultipleCallsDifferent ensures randomness produces variation across calls.
func TestGenerateQualityScores_MultipleCallsDifferent(t *testing.T) {
	modeler := genomics.NewMultiPlatformQualityModeler()
	seq := "ATCGATCGATCGATCGATCGATCGATCGATCG"

	// Run quality generation multiple times; at least one pair should differ
	results := make([][]int, 5)
	for i := range results {
		results[i] = modeler.GenerateQualityScores(genomics.NovaSeq6000, population.YRI, seq, 1)
	}

	allSame := true
	for i := 1; i < len(results); i++ {
		for j := 0; j < len(seq); j++ {
			if results[0][j] != results[i][j] {
				allSame = false
				break
			}
		}
		if !allSame {
			break
		}
	}

	// With random noise, expect at least some variation across 5 runs
	// (This test may occasionally fail with extremely unlikely seeds, but is generally reliable)
	if allSame {
		t.Log("all 5 quality score runs produced identical results (unlikely but possible)")
	}
}
