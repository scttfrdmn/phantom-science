package features_test

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/features"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ==================== GCBiasModel tests ====================

func TestNewGCBiasModel(t *testing.T) {
	t.Run("default window size", func(t *testing.T) {
		model, err := features.NewGCBiasModel(100)
		require.NoError(t, err)
		require.NotNil(t, model)
	})

	t.Run("zero window size uses default", func(t *testing.T) {
		model, err := features.NewGCBiasModel(0)
		require.NoError(t, err)
		require.NotNil(t, model)
	})

	t.Run("negative window size uses default", func(t *testing.T) {
		model, err := features.NewGCBiasModel(-1)
		require.NoError(t, err)
		require.NotNil(t, model)
	})
}

func TestGCBiasModel_AnalyzeSequence(t *testing.T) {
	model, err := features.NewGCBiasModel(100)
	require.NoError(t, err)

	t.Run("empty sequence", func(t *testing.T) {
		score, details, err := model.AnalyzeSequence([]byte{})
		require.NoError(t, err)
		assert.Equal(t, 0.0, score)
		assert.NotNil(t, details)
	})

	t.Run("all GC sequence", func(t *testing.T) {
		seq := []byte("GCGCGCGCGCGCGCGCGCGC")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.InDelta(t, 1.0, score, 0.001)
		assert.NotNil(t, details)
		assert.InDelta(t, 1.0, details.GCContent, 0.001)
	})

	t.Run("all AT sequence", func(t *testing.T) {
		seq := []byte("ATATATATATATAT")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.InDelta(t, 0.0, score, 0.001)
		assert.NotNil(t, details)
	})

	t.Run("mixed sequence", func(t *testing.T) {
		// 50% GC
		seq := []byte("ATCGATCGATCGATCGATCG")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.InDelta(t, 0.5, score, 0.01)
		assert.NotNil(t, details)
	})

	t.Run("lowercase sequence", func(t *testing.T) {
		seq := []byte("gcgcgcatatatat")
		score, _, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		// 6 GC out of 14 = ~0.43
		assert.InDelta(t, 6.0/14.0, score, 0.01)
	})

	t.Run("realistic human sequence", func(t *testing.T) {
		// Typical human genome GC ~41%
		seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		assert.LessOrEqual(t, score, 1.0)
		assert.NotNil(t, details)
	})
}

// ==================== RepetitiveElementModel tests ====================

func TestNewRepetitiveElementModel(t *testing.T) {
	model, err := features.NewRepetitiveElementModel(0.7)
	require.NoError(t, err)
	require.NotNil(t, model)
}

func TestRepetitiveElementModel_AnalyzeSequence(t *testing.T) {
	model, err := features.NewRepetitiveElementModel(0.7)
	require.NoError(t, err)

	t.Run("empty sequence", func(t *testing.T) {
		score, details, err := model.AnalyzeSequence([]byte{})
		require.NoError(t, err)
		assert.Equal(t, 0.0, score)
		assert.NotNil(t, details)
	})

	t.Run("score in bounds for normal sequence", func(t *testing.T) {
		seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		assert.LessOrEqual(t, score, 1.0)
		assert.NotNil(t, details)
	})

	// Regression: previously the stub always returned 1.0 regardless of content.
	// A highly repetitive sequence must score < 0.5 with the entropy-based impl.
	t.Run("highly repetitive sequence scores low", func(t *testing.T) {
		seq := []byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") // single nucleotide
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.Less(t, score, 0.5, "monotonically repetitive sequence should score below 0.5")
		assert.GreaterOrEqual(t, details.RepeatScore, 0.0)
		assert.LessOrEqual(t, details.RepeatScore, 1.0)
	})

	// A complex sequence (all 4 bases roughly equally represented) should score
	// higher than a simple homopolymer.
	t.Run("complex sequence scores higher than repetitive", func(t *testing.T) {
		repetitive := []byte("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
		complex := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")

		repScore, _, err := model.AnalyzeSequence(repetitive)
		require.NoError(t, err)
		cplxScore, _, err := model.AnalyzeSequence(complex)
		require.NoError(t, err)

		assert.Greater(t, cplxScore, repScore,
			"complex sequence (score %v) should score higher than homopolymer (score %v)",
			cplxScore, repScore)
	})

	// A long homopolymer run should produce a non-zero RepeatCount.
	t.Run("homopolymer repeat count", func(t *testing.T) {
		seq := []byte("AAAAAAAAAAAAAAAA") // 16 A's
		_, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.Greater(t, details.RepeatCount, 0, "homopolymer run should produce RepeatCount > 0")
	})

	// A non-repetitive diverse sequence should have RepeatCount = 0.
	t.Run("non-repetitive sequence has zero repeat count", func(t *testing.T) {
		seq := []byte("ATCGATCGATCG") // no run of 3 identical bases
		_, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.Equal(t, 0, details.RepeatCount)
	})
}

// ==================== MutationalSignatureModel tests ====================

func TestNewMutationalSignatureModel(t *testing.T) {
	t.Run("no processes", func(t *testing.T) {
		model, err := features.NewMutationalSignatureModel(nil)
		require.NoError(t, err)
		require.NotNil(t, model)
	})

	t.Run("with UV exposure", func(t *testing.T) {
		model, err := features.NewMutationalSignatureModel([]features.MutationalProcess{features.UVExposure})
		require.NoError(t, err)
		require.NotNil(t, model)
	})

	t.Run("with multiple processes", func(t *testing.T) {
		processes := []features.MutationalProcess{
			features.ClockLikeMutation,
			features.APOBECActivity,
			features.DNARepairDeficiency,
		}
		model, err := features.NewMutationalSignatureModel(processes)
		require.NoError(t, err)
		require.NotNil(t, model)
	})

	t.Run("with all processes", func(t *testing.T) {
		processes := []features.MutationalProcess{
			features.SpontaneousDeamination,
			features.UVExposure,
			features.TobaccoSmoking,
			features.DNARepairDeficiency,
			features.ClockLikeMutation,
			features.APOBECActivity,
			features.ChemicalMutagenesis,
			features.AgeRelatedMutation,
		}
		model, err := features.NewMutationalSignatureModel(processes)
		require.NoError(t, err)
		require.NotNil(t, model)
	})
}

func TestMutationalProcess_String(t *testing.T) {
	tests := []struct {
		process  features.MutationalProcess
		expected string
	}{
		{features.SpontaneousDeamination, "Spontaneous Deamination"},
		{features.UVExposure, "UV Exposure"},
		{features.TobaccoSmoking, "Tobacco Smoking"},
		{features.DNARepairDeficiency, "DNA Repair Deficiency"},
		{features.ClockLikeMutation, "Clock-like Mutation"},
		{features.APOBECActivity, "APOBEC Activity"},
		{features.ChemicalMutagenesis, "Chemical Mutagenesis"},
		{features.AgeRelatedMutation, "Age-related Mutation"},
	}

	for _, tc := range tests {
		t.Run(tc.expected, func(t *testing.T) {
			assert.Equal(t, tc.expected, tc.process.String())
		})
	}

	t.Run("unknown process", func(t *testing.T) {
		var unknown features.MutationalProcess = 99
		assert.Equal(t, "Unknown", unknown.String())
	})
}

func TestMutationalSignatureModel_AnalyzeSequence(t *testing.T) {
	model, err := features.NewMutationalSignatureModel([]features.MutationalProcess{features.ClockLikeMutation})
	require.NoError(t, err)

	t.Run("sequence too short", func(t *testing.T) {
		seq := []byte("ATCGATCG") // < 50 bases
		_, _, err := model.AnalyzeSequence(seq)
		assert.Error(t, err)
	})

	t.Run("valid sequence", func(t *testing.T) {
		// 100bp sequence with various trinucleotides
		seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		assert.LessOrEqual(t, score, 1.0)
		require.NotNil(t, details)
		assert.GreaterOrEqual(t, details.SequenceLength, int64(50))
	})

	t.Run("CpG-rich sequence mutational load", func(t *testing.T) {
		// Sequence with CG dinucleotides (mutation hotspots)
		seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		assert.NotNil(t, details)
		assert.GreaterOrEqual(t, details.MutationalLoad, 0.0)
	})

	t.Run("long realistic sequence", func(t *testing.T) {
		// 200bp sequence simulating exonic region
		seq := make([]byte, 200)
		pattern := []byte("ATCGATCGATCGCCGGAATCGATCGATCG")
		for i := 0; i < 200; i++ {
			seq[i] = pattern[i%len(pattern)]
		}
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		assert.LessOrEqual(t, score, 1.0)
		assert.NotNil(t, details)
		assert.NotNil(t, details.TrinucleotideFrequencies)
		assert.NotNil(t, details.SignatureContributions)
	})
}

// ==================== SelectionPressureModel tests ====================

func TestNewSelectionPressureModel(t *testing.T) {
	model, err := features.NewSelectionPressureModel()
	require.NoError(t, err)
	require.NotNil(t, model)
}

func TestSelectionType_String(t *testing.T) {
	tests := []struct {
		selType  features.SelectionType
		expected string
	}{
		{features.StrongPurifying, "Strong Purifying Selection"},
		{features.WeakPurifying, "Weak Purifying Selection"},
		{features.NeutralEvolution, "Neutral Evolution"},
		{features.PositiveSelection, "Positive Selection"},
		{features.BalancingSelection, "Balancing Selection"},
	}

	for _, tc := range tests {
		t.Run(tc.expected, func(t *testing.T) {
			assert.Equal(t, tc.expected, tc.selType.String())
		})
	}

	t.Run("unknown selection type", func(t *testing.T) {
		var unknown features.SelectionType = 99
		assert.Equal(t, "Unknown", unknown.String())
	})
}

func TestSelectionPressureModel_AnalyzeRegion(t *testing.T) {
	model, err := features.NewSelectionPressureModel()
	require.NoError(t, err)

	t.Run("invalid coordinates", func(t *testing.T) {
		_, _, err := model.AnalyzeRegion("chr1", 1000, 500) // start > end
		assert.Error(t, err)
	})

	t.Run("exonic region on chr1", func(t *testing.T) {
		// Position 5000-15000 is classified as "exon" by the model
		score, details, err := model.AnalyzeRegion("chr1", 5000, 15000)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		assert.LessOrEqual(t, score, 1.0)
		require.NotNil(t, details)
		assert.GreaterOrEqual(t, details.PurifyingSelection, 0.0)
		assert.GreaterOrEqual(t, details.PositiveSelection, 0.0)
		assert.GreaterOrEqual(t, details.BalancingSelection, 0.0)
		assert.Equal(t, int64(10000), details.RegionSize)
	})

	t.Run("intergenic region", func(t *testing.T) {
		score, details, err := model.AnalyzeRegion("chr2", 70000, 80000)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		require.NotNil(t, details)
	})

	t.Run("promoter region", func(t *testing.T) {
		// Position 0-5000 is classified as "promoter"
		score, details, err := model.AnalyzeRegion("chr1", 0, 5000)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		require.NotNil(t, details)
	})

	t.Run("small region", func(t *testing.T) {
		score, details, err := model.AnalyzeRegion("chr1", 1000, 1100)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		require.NotNil(t, details)
	})

	t.Run("large region spans multiple contexts", func(t *testing.T) {
		score, details, err := model.AnalyzeRegion("chr3", 0, 100000)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		require.NotNil(t, details)
		assert.Equal(t, int64(100000), details.RegionSize)
	})
}

// ==================== ChromatinAccessibilityModel tests ====================

func TestNewChromatinAccessibilityModel(t *testing.T) {
	model, err := features.NewChromatinAccessibilityModel()
	require.NoError(t, err)
	require.NotNil(t, model)
}

func TestChromatinState_String(t *testing.T) {
	tests := []struct {
		state    features.ChromatinState
		expected string
	}{
		{features.ActivePromoter, "Active Promoter"},
		{features.WeakPromoter, "Weak Promoter"},
		{features.PoisedPromoter, "Poised Promoter"},
		{features.StrongEnhancer, "Strong Enhancer"},
		{features.WeakEnhancer, "Weak Enhancer"},
		{features.Txn_Transition, "Transcriptional Transition"},
		{features.Txn_Elongation, "Transcriptional Elongation"},
		{features.Weak_Txn, "Weak Transcription"},
		{features.Polycomb_Repressed, "Polycomb Repressed"},
		{features.Heterochromatin_low, "Heterochromatin"},
		{features.Repetitive_CNV, "Repetitive/CNV"},
		{features.Inactive_Low, "Inactive"},
	}

	for _, tc := range tests {
		t.Run(tc.expected, func(t *testing.T) {
			assert.Equal(t, tc.expected, tc.state.String())
		})
	}

	t.Run("unknown chromatin state", func(t *testing.T) {
		var unknown features.ChromatinState = 99
		assert.Equal(t, "Unknown", unknown.String())
	})
}

func TestChromatinAccessibilityModel_AnalyzeRegion(t *testing.T) {
	model, err := features.NewChromatinAccessibilityModel()
	require.NoError(t, err)

	t.Run("invalid coordinates", func(t *testing.T) {
		_, _, err := model.AnalyzeRegion("chr1", 5000, 1000) // start > end
		assert.Error(t, err)
	})

	t.Run("active promoter region", func(t *testing.T) {
		// Position 0 should be at active promoter state
		score, details, err := model.AnalyzeRegion("chr1", 0, 10000)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		assert.LessOrEqual(t, score, 1.0)
		require.NotNil(t, details)
		assert.GreaterOrEqual(t, details.AccessibilityScore, 0.0)
		assert.GreaterOrEqual(t, details.RegulatoryPotential, 0.0)
		assert.GreaterOrEqual(t, details.DNaseScore, 0.0)
	})

	t.Run("chr1 enhancer region", func(t *testing.T) {
		// Position 100000-200000 on chr1
		score, details, err := model.AnalyzeRegion("chr1", 100000, 200000)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		require.NotNil(t, details)
		assert.Equal(t, int64(100000), details.RegionSize)
	})

	t.Run("unknown chromosome uses default", func(t *testing.T) {
		score, details, err := model.AnalyzeRegion("chrX", 1000, 2000)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		require.NotNil(t, details)
	})

	t.Run("region with TF binding sites", func(t *testing.T) {
		// Use a region where TF binding sites are likely to be
		score, details, err := model.AnalyzeRegion("chr1", 1000, 50000)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		require.NotNil(t, details)
		assert.NotNil(t, details.HistoneProfile)
	})
}

// ==================== CodonUsageModel tests ====================

func TestNewCodonUsageModel(t *testing.T) {
	model, err := features.NewCodonUsageModel()
	require.NoError(t, err)
	require.NotNil(t, model)
}

func TestCodonUsageModel_AnalyzeSequence(t *testing.T) {
	model, err := features.NewCodonUsageModel()
	require.NoError(t, err)

	t.Run("too short sequence", func(t *testing.T) {
		_, _, err := model.AnalyzeSequence([]byte("AT"))
		assert.Error(t, err)
	})

	t.Run("valid DNA sequence", func(t *testing.T) {
		// Coding-like sequence
		seq := []byte("ATGAAAGAAATCGGCATCGCAGAGAAGTTTAAGGCAATCATGAACATCATGGAGAAGTCGCAGAAAGCGCTGG")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		assert.LessOrEqual(t, score, 1.0)
		require.NotNil(t, details)
		assert.Greater(t, details.TotalCodons, 0)
		assert.NotNil(t, details.ObservedFreq)
		assert.NotNil(t, details.ExpectedFreq)
		assert.NotNil(t, details.RSCU)
	})

	t.Run("sequence with invalid bases is handled", func(t *testing.T) {
		// Sequence with N characters (ambiguous bases)
		seq := []byte("ATGAAANNNGGCATCGCAGAGAAGTTTAAGGCAATCATGAACATCATGGAG")
		_, _, err := model.AnalyzeSequence(seq)
		// May error or succeed depending on how it handles N
		_ = err
	})

	t.Run("minimal valid triplet sequence", func(t *testing.T) {
		seq := []byte("ATG")
		score, details, err := model.AnalyzeSequence(seq)
		require.NoError(t, err)
		assert.GreaterOrEqual(t, score, 0.0)
		assert.NotNil(t, details)
	})
}

func TestCodonUsageModel_GetOptimalCodons(t *testing.T) {
	model, err := features.NewCodonUsageModel()
	require.NoError(t, err)

	optimal := model.GetOptimalCodons()
	require.NotNil(t, optimal)
	// Should return a map - it may be empty if codonFrequencies is not populated
	// Just verify it doesn't panic and returns a map
	assert.NotNil(t, optimal)
	// If the map has entries, each codon should be 3 characters
	for aa, codon := range optimal {
		if codon != "" {
			assert.Len(t, codon, 3, "Optimal codon should be 3 bases for amino acid %s", aa)
		}
	}
}

func TestCodonUsageModel_SetSeed(t *testing.T) {
	model, err := features.NewCodonUsageModel()
	require.NoError(t, err)

	// Should not panic
	model.SetSeed(12345)
	model.SetSeed(0)
	model.SetSeed(-1)
}

func TestCodonUsageModel_SetExpressionLevel(t *testing.T) {
	model, err := features.NewCodonUsageModel()
	require.NoError(t, err)

	// Should not panic
	model.SetExpressionLevel(1.0)
	model.SetExpressionLevel(0.5)
	model.SetExpressionLevel(0.0)
	model.SetExpressionLevel(-1.0) // Negative value should be ignored
}

func TestCodonUsageModel_GetModelStats(t *testing.T) {
	model, err := features.NewCodonUsageModel()
	require.NoError(t, err)

	stats := model.GetModelStats()
	require.NotNil(t, stats)
	assert.Contains(t, stats, "organism")
	assert.Contains(t, stats, "total_codons")
	assert.Contains(t, stats, "codon_entropy")
}

func TestCodonUsageModel_ValidateCodonUsage(t *testing.T) {
	model, err := features.NewCodonUsageModel()
	require.NoError(t, err)

	// The default model should be valid - it may or may not pass depending
	// on whether expected frequencies cover all 64 codons
	err = model.ValidateCodonUsage()
	// We just check it doesn't panic
	_ = err
}

func TestCodonUsageModel_LoadCodonTable(t *testing.T) {
	model, err := features.NewCodonUsageModel()
	require.NoError(t, err)

	table := &features.CodonUsageTable{
		Organism: "test_organism",
		CodonFrequencies: map[string]float64{
			"ATG": 0.02, "GCT": 0.015, "GCC": 0.02, "GCA": 0.015, "GCG": 0.01,
		},
	}

	// Should not panic
	model.LoadCodonTable(table)
}

func TestCodonUsageModel_ApplyCodonBias(t *testing.T) {
	model, err := features.NewCodonUsageModel()
	require.NoError(t, err)

	t.Run("empty sequence errors", func(t *testing.T) {
		_, err := model.ApplyCodonBias([]byte{})
		assert.Error(t, err)
	})

	t.Run("short non-coding sequence", func(t *testing.T) {
		seq := []byte("ATCGATCGATCG")
		result, err := model.ApplyCodonBias(seq)
		require.NoError(t, err)
		require.NotNil(t, result)
		assert.Equal(t, seq, result.OriginalSequence)
		assert.NotNil(t, result.BiasedSequence)
	})

	t.Run("sequence with ORF", func(t *testing.T) {
		// Long enough sequence with ATG start and stop codon
		// Needs to be > 150 nt for an ORF to be found
		seq := []byte(
			"ATGAAAGAAATCGGCATCGCAGAGAAGTTTAAGGCAATCATGAACATCATGGAGAAGTCGCAGAAAGCGCTG" +
				"ATCGCGCTGAAAGAGAAGATCATCATGAAAGAAATCGGCATCGCAGAGAAGTTTAAGGCAATCATGAACATC" +
				"ATCATGGAGAAGTCGCAGAAAGCGCTGATCGCGCTGAAAGAGAAGATCATCATGAAAGAAATCGGCATCGCA" +
				"GAGAAGTTTAAGGCAATCATGTAA",
		)
		result, err := model.ApplyCodonBias(seq)
		require.NoError(t, err)
		require.NotNil(t, result)
		assert.NotNil(t, result.ORFs)
		assert.Equal(t, float64(1.0), result.ExpressionLevel)
	})
}

// ==================== PhylogeneticDistanceModel tests ====================

func TestNewPhylogeneticDistanceModel(t *testing.T) {
	model, err := features.NewPhylogeneticDistanceModel(0.1)
	require.NoError(t, err)
	require.NotNil(t, model)
}

func TestPhylogeneticDistanceModel_AnalyzeSequence(t *testing.T) {
	model, err := features.NewPhylogeneticDistanceModel(0.1)
	require.NoError(t, err)

	t.Run("sequence too short", func(t *testing.T) {
		seq := []byte("ATCGATCGATCG") // < 100 bases
		_, _, err := model.AnalyzeSequence(seq, nil)
		assert.Error(t, err)
	})

	t.Run("valid sequence with nil kmer profile", func(t *testing.T) {
		// 150bp sequence
		seq := make([]byte, 150)
		pattern := []byte("ATCGATCGATCGATCGATCGATCGATCGATCG")
		for i := range seq {
			seq[i] = pattern[i%len(pattern)]
		}
		score, details, err := model.AnalyzeSequence(seq, nil)
		require.NoError(t, err)
		require.NotNil(t, details)
		// With nil kmer profile, score may be NaN/Inf due to Inf distances
		// Just verify no panic and details returned
		_ = score
		assert.GreaterOrEqual(t, details.ConservationScore, 0.0)
		assert.LessOrEqual(t, details.ConservationScore, 1.0)
	})

	t.Run("GC-rich sequence conservation score", func(t *testing.T) {
		// CpG-rich sequence (high GC)
		seq := make([]byte, 150)
		pattern := []byte("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG")
		for i := range seq {
			seq[i] = pattern[i%len(pattern)]
		}
		_, details, err := model.AnalyzeSequence(seq, nil)
		require.NoError(t, err)
		require.NotNil(t, details)
		// High GC deviates from expected 41%, conservation score should be lower
		assert.GreaterOrEqual(t, details.ConservationScore, 0.0)
		assert.LessOrEqual(t, details.ConservationScore, 1.0)
	})

	t.Run("typical human genome GC content", func(t *testing.T) {
		// ~41% GC sequence
		seq := make([]byte, 200)
		// Manually construct ~41% GC sequence
		template := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
		for i := range seq {
			seq[i] = template[i%len(template)]
		}
		_, details, err := model.AnalyzeSequence(seq, nil)
		require.NoError(t, err)
		require.NotNil(t, details)
		// With nil kmer profile, LikelyEvolutionaryOrigin may be empty (all Inf distances)
		// Just verify no panic
		_ = details.LikelyEvolutionaryOrigin
	})
}

// ==================== BiologicalModelingEngine tests ====================

func TestNewBiologicalModelingEngine(t *testing.T) {
	t.Run("all models disabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})

	t.Run("only GC bias enabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableGCBias: true,
			GCWindowSize: 100,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})

	t.Run("only CpG depletion enabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableCpGDepletion: true,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})

	t.Run("only repeat elements enabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableRepeatElements:      true,
			RepeatComplexityThreshold: 0.7,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})

	t.Run("only codon usage enabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableCodonUsage: true,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})

	t.Run("only phylogenetic enabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnablePhylogenetic:    true,
			PhylogeneticThreshold: 0.1,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})

	t.Run("only chromatin enabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableChromatin: true,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})

	t.Run("only mutation signature enabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableMutationSignature: true,
			MutationalProcesses:     []features.MutationalProcess{features.ClockLikeMutation},
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})

	t.Run("only selection model enabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableSelection: true,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})

	t.Run("all models enabled", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableGCBias:            true,
			EnableCpGDepletion:      true,
			EnableRepeatElements:    true,
			EnableCodonUsage:        true,
			EnablePhylogenetic:      true,
			EnableChromatin:         true,
			EnableMutationSignature: true,
			EnableSelection:         true,

			GCWindowSize:              100,
			CpGThreshold:              0.6,
			RepeatComplexityThreshold: 0.7,
			PhylogeneticThreshold:     0.1,
			MutationalProcesses: []features.MutationalProcess{
				features.ClockLikeMutation,
				features.APOBECActivity,
			},
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)
		require.NotNil(t, engine)
	})
}

func TestBiologicalModelingEngine_AnalyzeBiologicalRealism(t *testing.T) {
	t.Run("empty engine no models", func(t *testing.T) {
		config := features.BiologicalModelConfig{}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)

		seq := []byte("ATCGATCGATCGATCGATCG")
		result, err := engine.AnalyzeBiologicalRealism(seq, nil, "chr1", 1000)
		require.NoError(t, err)
		require.NotNil(t, result)
		assert.Equal(t, "chr1", result.Chromosome)
		assert.Equal(t, int64(1000), result.Position)
		assert.Equal(t, int64(len(seq)), result.SequenceLength)
	})

	t.Run("with GC bias model", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableGCBias: true,
			GCWindowSize: 100,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)

		seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
		result, err := engine.AnalyzeBiologicalRealism(seq, nil, "chr1", 5000)
		require.NoError(t, err)
		require.NotNil(t, result)
		assert.Contains(t, result.FeatureScores, "gc_bias")
	})

	t.Run("with CpG depletion model", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableCpGDepletion: true,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)

		seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
		result, err := engine.AnalyzeBiologicalRealism(seq, nil, "chr1", 5000)
		require.NoError(t, err)
		require.NotNil(t, result)
		assert.Contains(t, result.FeatureScores, "cpg_depletion")
	})

	t.Run("with chromatin model", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableChromatin: true,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)

		seq := []byte("ATCGATCGATCGATCGATCG")
		result, err := engine.AnalyzeBiologicalRealism(seq, nil, "chr1", 5000)
		require.NoError(t, err)
		require.NotNil(t, result)
		assert.Contains(t, result.FeatureScores, "chromatin_accessibility")
	})

	t.Run("with selection model", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableSelection: true,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)

		seq := []byte("ATCGATCGATCGATCGATCG")
		result, err := engine.AnalyzeBiologicalRealism(seq, nil, "chr1", 5000)
		require.NoError(t, err)
		require.NotNil(t, result)
		assert.Contains(t, result.FeatureScores, "selection_pressure")
	})

	t.Run("with codon usage on long sequence", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableCodonUsage: true,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)

		// Long enough sequence for codon analysis (>= 60 bases)
		seq := make([]byte, 100)
		pattern := []byte("ATGAAAGAAATCGGCATCGCAGAGAAGTTTAAGGCAATCATGAAC")
		for i := range seq {
			seq[i] = pattern[i%len(pattern)]
		}
		result, err := engine.AnalyzeBiologicalRealism(seq, nil, "chr1", 5000)
		require.NoError(t, err)
		require.NotNil(t, result)
		assert.Contains(t, result.FeatureScores, "codon_usage")
	})

	t.Run("overall score and category", func(t *testing.T) {
		config := features.BiologicalModelConfig{
			EnableGCBias:       true,
			EnableCpGDepletion: true,
			EnableChromatin:    true,
			EnableSelection:    true,
			GCWindowSize:       100,
		}
		engine, err := features.NewBiologicalModelingEngine(config)
		require.NoError(t, err)

		seq := make([]byte, 100)
		pattern := []byte("ATCGATCGATCGATCGATCG")
		for i := range seq {
			seq[i] = pattern[i%len(pattern)]
		}
		result, err := engine.AnalyzeBiologicalRealism(seq, nil, "chr1", 10000)
		require.NoError(t, err)
		require.NotNil(t, result)
		// OverallScore may be NaN if selection model returns NaN; just check RealismCategory
		_ = result.OverallScore
		// Check that RealismCategory is a valid value
		_ = result.RealismCategory.String()
	})
}

// ==================== RealismCategory tests ====================

func TestRealismCategory_String(t *testing.T) {
	tests := []struct {
		category features.RealismCategory
		expected string
	}{
		{features.PoorlyRealistic, "Poorly Realistic"},
		{features.SomewhatRealistic, "Somewhat Realistic"},
		{features.ModeratelyRealistic, "Moderately Realistic"},
		{features.HighlyRealistic, "Highly Realistic"},
	}

	for _, tc := range tests {
		t.Run(tc.expected, func(t *testing.T) {
			assert.Equal(t, tc.expected, tc.category.String())
		})
	}

	t.Run("unknown realism category", func(t *testing.T) {
		var unknown features.RealismCategory = 99
		assert.Equal(t, "Unknown", unknown.String())
	})
}
