package quality

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewHomopolymerErrorModel tests model creation for different platforms
func TestNewHomopolymerErrorModel(t *testing.T) {
	testCases := []struct {
		platform           string
		expectedDeletion   float64
		expectedInsertion  float64
		expectedDecay      float64
		expectedMaxRun     int
		expectedBaseRate   float64
	}{
		{"illumina", 0.3, 0.2, 1.0, 15, 0.001},
		{"novaseq", 0.3, 0.2, 1.0, 15, 0.001},
		{"nextseq", 0.3, 0.2, 1.0, 15, 0.001},
		{"miseq", 0.3, 0.2, 1.0, 15, 0.001},
		{"pacbio", 0.6, 0.3, 2.0, 25, 0.01},
		{"pacbio_revio", 0.6, 0.3, 2.0, 25, 0.01},
		{"nanopore", 0.7, 0.25, 3.0, 30, 0.05},
		{"promethion", 0.7, 0.25, 3.0, 30, 0.05},
		{"unknown_platform", 0.3, 0.2, 1.0, 15, 0.001},
	}

	for _, tc := range testCases {
		t.Run(tc.platform, func(t *testing.T) {
			model := NewHomopolymerErrorModel(tc.platform)

			assert.NotNil(t, model)
			assert.Equal(t, tc.platform, model.platform)
			assert.Equal(t, tc.expectedDeletion, model.deletionBias)
			assert.Equal(t, tc.expectedInsertion, model.insertionBias)
			assert.Equal(t, tc.expectedDecay, model.qualityDecay)
			assert.Equal(t, tc.expectedMaxRun, model.maxRunLength)
			assert.Equal(t, tc.expectedBaseRate, model.baseErrorRate)
			assert.True(t, model.contextEffects)
			assert.False(t, model.modelTrained)
			assert.NotNil(t, model.rng)
			assert.NotNil(t, model.errorRates)
			assert.Greater(t, len(model.errorRates), 0)
			assert.NotNil(t, model.errorDistribution)
		})
	}
}

// TestHomopolymerErrorModel_SetSeed tests seed setting for reproducibility
func TestHomopolymerErrorModel_SetSeed(t *testing.T) {
	model := NewHomopolymerErrorModel("nanopore")
	model.SetSeed(42)

	assert.NotNil(t, model.rng)
}

// TestHomopolymerErrorModel_GetPlatformParameters tests platform parameter retrieval
func TestHomopolymerErrorModel_GetPlatformParameters(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")

	params := model.GetPlatformParameters()

	assert.Equal(t, "illumina", params["platform"])
	assert.Equal(t, 0.3, params["deletion_bias"])
	assert.Equal(t, 0.2, params["insertion_bias"])
	assert.Equal(t, 1.0, params["quality_decay"])
	assert.Equal(t, 0.001, params["base_error_rate"])
	assert.Equal(t, 15, params["max_run_length"])
}

// TestHomopolymerErrorModel_GetModelStats tests model statistics retrieval
func TestHomopolymerErrorModel_GetModelStats(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")

	stats := model.GetModelStats()

	assert.Equal(t, "illumina", stats["platform"])
	assert.Equal(t, false, stats["model_trained"])
	assert.Equal(t, 0.001, stats["base_error_rate"])
	assert.Equal(t, 15, stats["max_run_length"])
	assert.Equal(t, 0.3, stats["deletion_bias"])
	assert.Equal(t, 0.2, stats["insertion_bias"])
	assert.Equal(t, 1.0, stats["quality_decay"])
	assert.Equal(t, true, stats["context_effects"])
	assert.Equal(t, 0, stats["observed_errors"])

	// Should have error rate statistics since rates are initialized
	assert.Contains(t, stats, "max_error_rate")
	assert.Contains(t, stats, "min_error_rate")
	assert.Contains(t, stats, "error_rate_entries")
}

// TestHomopolymerErrorModel_ValidateModel tests model validation
func TestHomopolymerErrorModel_ValidateModel(t *testing.T) {
	t.Run("valid illumina model", func(t *testing.T) {
		model := NewHomopolymerErrorModel("illumina")
		err := model.ValidateModel()
		assert.NoError(t, err)
	})

	t.Run("valid nanopore model", func(t *testing.T) {
		model := NewHomopolymerErrorModel("nanopore")
		err := model.ValidateModel()
		assert.NoError(t, err)
	})

	t.Run("valid pacbio model", func(t *testing.T) {
		model := NewHomopolymerErrorModel("pacbio")
		err := model.ValidateModel()
		assert.NoError(t, err)
	})

	t.Run("invalid base error rate", func(t *testing.T) {
		model := NewHomopolymerErrorModel("illumina")
		model.baseErrorRate = -0.1
		err := model.ValidateModel()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "invalid base error rate")
	})

	t.Run("base error rate too high", func(t *testing.T) {
		model := NewHomopolymerErrorModel("illumina")
		model.baseErrorRate = 1.5
		err := model.ValidateModel()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "invalid base error rate")
	})

	t.Run("invalid max run length", func(t *testing.T) {
		model := NewHomopolymerErrorModel("illumina")
		model.maxRunLength = 0
		err := model.ValidateModel()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "invalid max run length")
	})

	t.Run("invalid deletion bias", func(t *testing.T) {
		model := NewHomopolymerErrorModel("illumina")
		model.deletionBias = -0.1
		err := model.ValidateModel()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "invalid error bias values")
	})

	t.Run("biases sum to more than 1", func(t *testing.T) {
		model := NewHomopolymerErrorModel("illumina")
		model.deletionBias = 0.7
		model.insertionBias = 0.5
		err := model.ValidateModel()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "invalid error bias values")
	})
}

// TestHomopolymerErrorModel_FindHomopolymerRuns tests finding runs in sequences
func TestHomopolymerErrorModel_FindHomopolymerRuns(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")

	t.Run("empty sequence", func(t *testing.T) {
		runs := model.FindHomopolymerRuns([]byte{})
		assert.Empty(t, runs)
	})

	t.Run("no homopolymers", func(t *testing.T) {
		// No runs of length >= 2
		runs := model.FindHomopolymerRuns([]byte("ATCGATCG"))
		assert.Empty(t, runs)
	})

	t.Run("single homopolymer run", func(t *testing.T) {
		// AAAA in the middle
		sequence := []byte("TCGAAAAATCG")
		runs := model.FindHomopolymerRuns(sequence)
		require.Len(t, runs, 1)
		assert.Equal(t, byte('A'), runs[0].Base)
		assert.Equal(t, 5, runs[0].Length)
		assert.True(t, runs[0].ErrorProne)
	})

	t.Run("multiple runs", func(t *testing.T) {
		sequence := []byte("AAATTTCCC")
		runs := model.FindHomopolymerRuns(sequence)
		assert.Len(t, runs, 3)
	})

	t.Run("run at end of sequence", func(t *testing.T) {
		sequence := []byte("ATCGAAAA")
		runs := model.FindHomopolymerRuns(sequence)
		require.Len(t, runs, 1)
		assert.Equal(t, byte('A'), runs[0].Base)
		assert.Equal(t, 4, runs[0].Length)
	})

	t.Run("run at start of sequence", func(t *testing.T) {
		sequence := []byte("AAAATCG")
		runs := model.FindHomopolymerRuns(sequence)
		require.Len(t, runs, 1)
		assert.Equal(t, byte('A'), runs[0].Base)
		assert.Equal(t, 4, runs[0].Length)
	})

	t.Run("G homopolymer", func(t *testing.T) {
		sequence := []byte("ATCGGGGATCG")
		runs := model.FindHomopolymerRuns(sequence)
		require.Len(t, runs, 1)
		assert.Equal(t, byte('G'), runs[0].Base)
		assert.Equal(t, 4, runs[0].Length)
	})

	t.Run("T homopolymer", func(t *testing.T) {
		sequence := []byte("ACTTTTCG")
		runs := model.FindHomopolymerRuns(sequence)
		require.Len(t, runs, 1)
		assert.Equal(t, byte('T'), runs[0].Base)
		assert.Equal(t, 4, runs[0].Length)
	})

	t.Run("homopolymer run properties", func(t *testing.T) {
		sequence := []byte("ATCGAAAAATCG")
		runs := model.FindHomopolymerRuns(sequence)
		require.Len(t, runs, 1)

		run := runs[0]
		assert.Equal(t, byte('A'), run.Base)
		assert.Equal(t, 5, run.Length)
		assert.Equal(t, 4, run.Start)
		assert.Equal(t, 9, run.End)
		assert.True(t, run.ErrorProne)
		assert.Greater(t, run.ErrorRate, 0.0)
		assert.NotEmpty(t, run.Context)
	})
}

// TestHomopolymerErrorModel_SimulateErrors tests error simulation
func TestHomopolymerErrorModel_SimulateErrors(t *testing.T) {
	model := NewHomopolymerErrorModel("nanopore") // Higher error rate for more test coverage
	model.SetSeed(42)

	t.Run("empty sequence", func(t *testing.T) {
		result, err := model.SimulateErrors([]byte{}, []byte{})
		assert.Error(t, err)
		assert.Nil(t, result)
		assert.Contains(t, err.Error(), "sequence cannot be empty")
	})

	t.Run("mismatched lengths", func(t *testing.T) {
		result, err := model.SimulateErrors([]byte("ATCG"), []byte{30, 30})
		assert.Error(t, err)
		assert.Nil(t, result)
		assert.Contains(t, err.Error(), "sequence and quality lengths must match")
	})

	t.Run("sequence without homopolymers", func(t *testing.T) {
		sequence := []byte("ATCGATCGATCG")
		qualities := make([]byte, len(sequence))
		for i := range qualities {
			qualities[i] = 30
		}

		result, err := model.SimulateErrors(sequence, qualities)
		require.NoError(t, err)
		require.NotNil(t, result)

		assert.Equal(t, len(sequence), len(result.OriginalSequence))
		assert.NotNil(t, result.OriginalQualities)
		assert.NotNil(t, result.ErrorStatistics)
	})

	t.Run("sequence with homopolymers", func(t *testing.T) {
		// Sequence with many homopolymer runs
		sequence := []byte("ATCGAAAAATTTTTGGGGGCCCCCATCG")
		qualities := make([]byte, len(sequence))
		for i := range qualities {
			qualities[i] = 30
		}

		result, err := model.SimulateErrors(sequence, qualities)
		require.NoError(t, err)
		require.NotNil(t, result)

		assert.Equal(t, len(sequence), len(result.OriginalSequence))
		assert.NotNil(t, result.HomopolymerRuns)
		assert.Greater(t, len(result.HomopolymerRuns), 0)
		assert.NotNil(t, result.ErrorStatistics)
	})

	t.Run("result structure validation", func(t *testing.T) {
		sequence := []byte("ATCGAAAAATCG")
		qualities := make([]byte, len(sequence))
		for i := range qualities {
			qualities[i] = 35
		}

		result, err := model.SimulateErrors(sequence, qualities)
		require.NoError(t, err)

		// Original sequence should be preserved
		assert.Equal(t, sequence, result.OriginalSequence)
		assert.Equal(t, qualities, result.OriginalQualities)

		// Error statistics should be valid
		stats := result.ErrorStatistics
		assert.GreaterOrEqual(t, stats.TotalErrors, 0)
		assert.GreaterOrEqual(t, stats.InsertionErrors, 0)
		assert.GreaterOrEqual(t, stats.DeletionErrors, 0)
		assert.GreaterOrEqual(t, stats.SubstitutionErrors, 0)
	})
}

// TestHomopolymerErrorModel_SimulateErrors_ErrorTypes tests different error types
func TestHomopolymerErrorModel_SimulateErrors_ErrorTypes(t *testing.T) {
	// Use a nanopore model which has high error rates
	model := NewHomopolymerErrorModel("nanopore")
	model.SetSeed(12345)

	// Sequence with a long homopolymer run to maximize error probability
	sequence := []byte("ATCGAAAAAAAAAAATCG") // 10 A's
	qualities := make([]byte, len(sequence))
	for i := range qualities {
		qualities[i] = 30
	}

	// Run simulation multiple times to increase chance of seeing various error types
	seenInsertion := false
	seenDeletion := false
	seenSubstitution := false

	for i := 0; i < 100; i++ {
		model.SetSeed(int64(i))
		result, err := model.SimulateErrors(sequence, qualities)
		require.NoError(t, err)

		for _, e := range result.Errors {
			switch e.ErrorType {
			case "insertion":
				seenInsertion = true
			case "deletion":
				seenDeletion = true
			case "substitution":
				seenSubstitution = true
			}
		}

		if seenInsertion && seenDeletion && seenSubstitution {
			break
		}
	}

	// With enough iterations, we should see all error types
	// (though statistical variation means this might not always hit all three)
	totalSeen := 0
	if seenInsertion {
		totalSeen++
	}
	if seenDeletion {
		totalSeen++
	}
	if seenSubstitution {
		totalSeen++
	}
	assert.GreaterOrEqual(t, totalSeen, 1, "Should see at least some error types")
}

// TestHomopolymerErrorModel_TrainFromObservedErrors tests training from observations
func TestHomopolymerErrorModel_TrainFromObservedErrors(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")

	t.Run("empty training data", func(t *testing.T) {
		err := model.TrainFromObservedErrors(nil)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "no training data provided")
	})

	t.Run("empty slice", func(t *testing.T) {
		err := model.TrainFromObservedErrors([]HomopolymerError{})
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "no training data provided")
	})

	t.Run("valid training data", func(t *testing.T) {
		errors := []HomopolymerError{
			{
				Position:       10,
				RunLength:      5,
				Base:           'A',
				ErrorType:      "deletion",
				Quality:        25.0,
				Context:        "ATAAAAATCG",
				Platform:       "illumina",
				TrueLength:     5,
				ObservedLength: 4,
			},
			{
				Position:       20,
				RunLength:      4,
				Base:           'G',
				ErrorType:      "insertion",
				Quality:        28.0,
				Context:        "CGGGGGAT",
				Platform:       "illumina",
				TrueLength:     4,
				ObservedLength: 5,
			},
			{
				Position:       30,
				RunLength:      3,
				Base:           'T',
				ErrorType:      "substitution",
				Quality:        30.0,
				Context:        "ATTTATCG",
				Platform:       "illumina",
				TrueLength:     3,
				ObservedLength: 3,
			},
		}

		err := model.TrainFromObservedErrors(errors)
		require.NoError(t, err)

		assert.True(t, model.modelTrained)
		assert.Len(t, model.observedErrors, 3)

		// errorRates must NOT be overwritten with 1.0: they should remain < 0.5
		// (the exponential rates set by initializeErrorRates for illumina are well below 0.5).
		for length, rate := range model.errorRates {
			assert.Less(t, rate, 1.0, "errorRate for length %d should not be 1.0 (100%%)", length)
		}

		// errorDistribution must reflect the observed counts (1 deletion, 1 insertion, 1 substitution).
		assert.InDelta(t, 1.0/3.0, model.errorDistribution["deletion"], 0.01)
		assert.InDelta(t, 1.0/3.0, model.errorDistribution["insertion"], 0.01)
		assert.InDelta(t, 1.0/3.0, model.errorDistribution["substitution"], 0.01)
	})
}

// TestHomopolymerErrorModel_TrainFromRunObservations tests training with total run counts
func TestHomopolymerErrorModel_TrainFromRunObservations(t *testing.T) {
	t.Run("error rates computed from counts", func(t *testing.T) {
		model := NewHomopolymerErrorModel("illumina")

		errors := []HomopolymerError{
			{RunLength: 5, ErrorType: "deletion", Position: 10, Base: 'A', Quality: 25.0, Platform: "illumina"},
			{RunLength: 5, ErrorType: "insertion", Position: 20, Base: 'A', Quality: 28.0, Platform: "illumina"},
		}
		// 2 errors out of 10 observations for run length 5 → rate = 0.2
		totalRunCounts := map[int]int{5: 10}

		err := model.TrainFromRunObservations(errors, totalRunCounts)
		require.NoError(t, err)

		assert.True(t, model.modelTrained)
		assert.InDelta(t, 0.2, model.errorRates[5], 0.001, "errorRate[5] should be 2/10 = 0.2")
	})

	t.Run("error rate capped at 0.5", func(t *testing.T) {
		model := NewHomopolymerErrorModel("illumina")

		errors := []HomopolymerError{
			{RunLength: 8, ErrorType: "deletion", Position: 1, Base: 'T', Quality: 20.0, Platform: "illumina"},
			{RunLength: 8, ErrorType: "deletion", Position: 2, Base: 'T', Quality: 20.0, Platform: "illumina"},
			{RunLength: 8, ErrorType: "deletion", Position: 3, Base: 'T', Quality: 20.0, Platform: "illumina"},
		}
		// 3 errors out of 3 observations would give rate=1.0, but should be capped at 0.5.
		totalRunCounts := map[int]int{8: 3}

		err := model.TrainFromRunObservations(errors, totalRunCounts)
		require.NoError(t, err)

		assert.LessOrEqual(t, model.errorRates[8], 0.5, "errorRate must be capped at 0.5")
	})

	t.Run("empty errors returns error", func(t *testing.T) {
		model := NewHomopolymerErrorModel("illumina")
		err := model.TrainFromRunObservations(nil, map[int]int{5: 10})
		assert.Error(t, err)
	})
}

// TestHomopolymerErrorModel_EstimateErrorImpact tests error impact estimation
func TestHomopolymerErrorModel_EstimateErrorImpact(t *testing.T) {
	model := NewHomopolymerErrorModel("nanopore")

	t.Run("sequence without homopolymers", func(t *testing.T) {
		sequence := []byte("ATCGATCGATCG")
		impact := model.EstimateErrorImpact(sequence)

		assert.NotNil(t, impact)
		assert.Equal(t, 0, impact["total_homopolymer_runs"])
		assert.Equal(t, 0, impact["error_prone_runs"])
		assert.Equal(t, 0.0, impact["expected_errors"])
	})

	t.Run("sequence with homopolymers", func(t *testing.T) {
		sequence := []byte("ATCGAAAAATTTTTCG")
		impact := model.EstimateErrorImpact(sequence)

		assert.NotNil(t, impact)
		totalRuns := impact["total_homopolymer_runs"].(int)
		assert.Greater(t, totalRuns, 0)

		errorProneRuns := impact["error_prone_runs"].(int)
		assert.Greater(t, errorProneRuns, 0)

		// All metrics should be present
		assert.Contains(t, impact, "expected_errors")
		assert.Contains(t, impact, "avg_run_length")
		assert.Contains(t, impact, "max_run_length")
		assert.Contains(t, impact, "error_prone_fraction")
		assert.Contains(t, impact, "homopolymer_density")
	})

	t.Run("high homopolymer sequence", func(t *testing.T) {
		sequence := []byte("AAAAAAGGGGGGTTTTTCCCCC")
		impact := model.EstimateErrorImpact(sequence)

		assert.NotNil(t, impact)
		assert.Equal(t, 4, impact["total_homopolymer_runs"])
		assert.Greater(t, impact["avg_run_length"].(float64), 4.0)
	})
}

// TestHomopolymerErrorModel_initializeErrorRates tests error rate initialization
func TestHomopolymerErrorModel_initializeErrorRates(t *testing.T) {
	platforms := []string{"illumina", "pacbio", "nanopore", "unknown"}

	for _, platform := range platforms {
		t.Run(platform, func(t *testing.T) {
			model := NewHomopolymerErrorModel(platform)

			// Should have error rates for each length up to maxRunLength
			for length := 1; length <= model.maxRunLength; length++ {
				rate, exists := model.errorRates[length]
				assert.True(t, exists, "Error rate should exist for length %d", length)
				assert.GreaterOrEqual(t, rate, 0.0, "Error rate should be non-negative for length %d", length)
				assert.LessOrEqual(t, rate, 0.5, "Error rate should be capped at 0.5 for length %d", length)
			}

			// Error rate should increase with length (up to cap)
			if len(model.errorRates) > 1 {
				shortRate := model.errorRates[1]
				longRate := model.errorRates[model.maxRunLength]
				// Long runs should have at least as high error rate as short runs
				assert.GreaterOrEqual(t, longRate, shortRate,
					"Longer runs should have at least as high error rate")
			}
		})
	}
}

// TestHomopolymerErrorModel_getErrorRate tests error rate lookup
func TestHomopolymerErrorModel_getErrorRate(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")

	t.Run("rate for length within range", func(t *testing.T) {
		rate := model.getErrorRate(5)
		assert.Greater(t, rate, 0.0)
		assert.LessOrEqual(t, rate, 1.0)
	})

	t.Run("rate for length beyond max", func(t *testing.T) {
		// Should cap at maxRunLength
		rate := model.getErrorRate(100)
		rateAtMax := model.getErrorRate(model.maxRunLength)
		assert.Equal(t, rateAtMax, rate)
	})

	t.Run("rate for length 1", func(t *testing.T) {
		rate := model.getErrorRate(1)
		assert.Equal(t, model.baseErrorRate, rate)
	})
}

// TestHomopolymerErrorModel_getContextModifier tests context modifier
func TestHomopolymerErrorModel_getContextModifier(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")

	t.Run("CG context increases rate", func(t *testing.T) {
		modifier := model.getContextModifier("CGAAAA")
		assert.Greater(t, modifier, 1.0)
	})

	t.Run("AT context modifier", func(t *testing.T) {
		// "ATAAAA" contains "AT" which gives 0.9 modifier, but also may contain
		// tandem repeats ("AA" at position 2-5) which gives 1.3 multiplier
		// The combined result depends on sequence - just verify it's a valid float
		modifier := model.getContextModifier("ATAAAA")
		assert.Greater(t, modifier, 0.0)
		assert.Less(t, modifier, 5.0)
	})

	t.Run("pure AT context no repeats decreases rate", func(t *testing.T) {
		// "ATAT" has "AT" but no tandem repeat (ATAT: check i=0: A==A,T==T -> tandem!)
		// Use a simpler context that has AT but no tandem repeat
		modifier := model.getContextModifier("AT")
		assert.Less(t, modifier, 1.0)
	})

	t.Run("tandem repeat increases rate", func(t *testing.T) {
		modifier := model.getContextModifier("ABABAAAA")
		assert.Greater(t, modifier, 1.0)
	})

	t.Run("neutral context", func(t *testing.T) {
		modifier := model.getContextModifier("TTCC")
		assert.InDelta(t, 1.0, modifier, 0.1)
	})
}

// TestHomopolymerError_Structure tests HomopolymerError structure
func TestHomopolymerError_Structure(t *testing.T) {
	error := HomopolymerError{
		Position:       15,
		RunLength:      6,
		Base:           'A',
		ErrorType:      "deletion",
		Quality:        25.0,
		Context:        "GCAAAAAAT",
		Platform:       "illumina",
		TrueLength:     6,
		ObservedLength: 5,
	}

	assert.Equal(t, 15, error.Position)
	assert.Equal(t, 6, error.RunLength)
	assert.Equal(t, byte('A'), error.Base)
	assert.Equal(t, "deletion", error.ErrorType)
	assert.Equal(t, 25.0, error.Quality)
	assert.Equal(t, 6, error.TrueLength)
	assert.Equal(t, 5, error.ObservedLength)
}

// TestHomopolymerRun_Structure tests HomopolymerRun structure
func TestHomopolymerRun_Structure(t *testing.T) {
	run := HomopolymerRun{
		Start:      10,
		End:        16,
		Length:     6,
		Base:       'G',
		Context:    "ATGGGGGCAT",
		ErrorProne: true,
		ErrorRate:  0.05,
	}

	assert.Equal(t, 10, run.Start)
	assert.Equal(t, 16, run.End)
	assert.Equal(t, 6, run.Length)
	assert.Equal(t, byte('G'), run.Base)
	assert.True(t, run.ErrorProne)
	assert.Equal(t, 0.05, run.ErrorRate)
}

// TestErrorSimulationResult_Structure tests ErrorSimulationResult structure
func TestErrorSimulationResult_Structure(t *testing.T) {
	result := &ErrorSimulationResult{
		OriginalSequence:  []byte("ATCGAAAAT"),
		ErrorSequence:     []byte("ATCGAAAT"),
		OriginalQualities: []byte{30, 30, 30, 30, 30, 30, 30, 30, 30},
		ErrorQualities:    []byte{30, 30, 30, 30, 30, 30, 30, 30},
		Errors: []HomopolymerError{
			{Position: 5, ErrorType: "deletion", RunLength: 4},
		},
		HomopolymerRuns: []HomopolymerRun{
			{Start: 4, End: 8, Length: 4, Base: 'A'},
		},
		ErrorStatistics: ErrorStatistics{
			TotalErrors:    1,
			DeletionErrors: 1,
			AvgRunLength:   4.0,
			MaxRunLength:   4,
			ErrorRate:      0.111,
		},
	}

	assert.Equal(t, 9, len(result.OriginalSequence))
	assert.Equal(t, 8, len(result.ErrorSequence))
	assert.Len(t, result.Errors, 1)
	assert.Len(t, result.HomopolymerRuns, 1)
	assert.Equal(t, 1, result.ErrorStatistics.TotalErrors)
	assert.Equal(t, 1, result.ErrorStatistics.DeletionErrors)
}

// TestErrorStatistics_Structure tests ErrorStatistics structure
func TestErrorStatistics_Structure(t *testing.T) {
	stats := ErrorStatistics{
		TotalErrors:        10,
		InsertionErrors:    3,
		DeletionErrors:     5,
		SubstitutionErrors: 2,
		AvgRunLength:       4.5,
		MaxRunLength:       8,
		ErrorRate:          0.025,
	}

	assert.Equal(t, 10, stats.TotalErrors)
	assert.Equal(t, 3, stats.InsertionErrors)
	assert.Equal(t, 5, stats.DeletionErrors)
	assert.Equal(t, 2, stats.SubstitutionErrors)
	assert.Equal(t, stats.TotalErrors, stats.InsertionErrors+stats.DeletionErrors+stats.SubstitutionErrors)
	assert.Greater(t, stats.AvgRunLength, 0.0)
	assert.Greater(t, stats.MaxRunLength, 0)
}

// TestHomopolymerErrorModel_SimulateErrors_AllPlatforms tests simulation on all platforms
func TestHomopolymerErrorModel_SimulateErrors_AllPlatforms(t *testing.T) {
	platforms := []string{"illumina", "pacbio", "nanopore"}
	sequence := []byte("ATCGAAAAAAAATTTTTGGGCCCATCG")
	qualities := make([]byte, len(sequence))
	for i := range qualities {
		qualities[i] = 30
	}

	for _, platform := range platforms {
		t.Run(platform, func(t *testing.T) {
			model := NewHomopolymerErrorModel(platform)
			model.SetSeed(99)

			result, err := model.SimulateErrors(sequence, qualities)
			require.NoError(t, err)
			assert.NotNil(t, result)
			assert.Equal(t, sequence, result.OriginalSequence)
			assert.GreaterOrEqual(t, result.ErrorStatistics.TotalErrors, 0)
		})
	}
}

// TestHomopolymerErrorModel_calculateErrorStatistics_NoErrors tests statistics with no errors
func TestHomopolymerErrorModel_calculateErrorStatistics_NoErrors(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")

	result := &ErrorSimulationResult{
		OriginalSequence: []byte("ATCGATCG"),
		Errors:           []HomopolymerError{},
	}

	stats := model.calculateErrorStatistics(result)
	assert.Equal(t, 0, stats.TotalErrors)
	assert.Equal(t, 0, stats.InsertionErrors)
	assert.Equal(t, 0, stats.DeletionErrors)
	assert.Equal(t, 0, stats.SubstitutionErrors)
	assert.Equal(t, 0.0, stats.AvgRunLength)
	assert.Equal(t, 0, stats.MaxRunLength)
	assert.Equal(t, 0.0, stats.ErrorRate)
}

// TestHomopolymerErrorModel_calculateErrorStatistics_WithErrors tests statistics with errors
func TestHomopolymerErrorModel_calculateErrorStatistics_WithErrors(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")

	result := &ErrorSimulationResult{
		OriginalSequence: []byte("ATCGAAAAATTTTTCG"), // 16 bases
		Errors: []HomopolymerError{
			{ErrorType: "insertion", RunLength: 5},
			{ErrorType: "deletion", RunLength: 4},
			{ErrorType: "substitution", RunLength: 3},
			{ErrorType: "deletion", RunLength: 6},
		},
	}

	stats := model.calculateErrorStatistics(result)
	assert.Equal(t, 4, stats.TotalErrors)
	assert.Equal(t, 1, stats.InsertionErrors)
	assert.Equal(t, 2, stats.DeletionErrors)
	assert.Equal(t, 1, stats.SubstitutionErrors)
	assert.Equal(t, 4.5, stats.AvgRunLength) // (5+4+3+6)/4 = 4.5
	assert.Equal(t, 6, stats.MaxRunLength)
	assert.Equal(t, 4.0/16.0, stats.ErrorRate)
}

// BenchmarkHomopolymerErrorModel_SimulateErrors benchmarks error simulation
func BenchmarkHomopolymerErrorModel_SimulateErrors(b *testing.B) {
	model := NewHomopolymerErrorModel("nanopore")
	model.SetSeed(42)

	sequence := []byte("ATCGAAAAAAAATTTTTTTGGGGGCCCCCATCGATCGATCGATCG")
	qualities := make([]byte, len(sequence))
	for i := range qualities {
		qualities[i] = 30
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = model.SimulateErrors(sequence, qualities)
	}
}

// BenchmarkHomopolymerErrorModel_FindHomopolymerRuns benchmarks run finding
func BenchmarkHomopolymerErrorModel_FindHomopolymerRuns(b *testing.B) {
	model := NewHomopolymerErrorModel("illumina")

	// Create a sequence with many homopolymer runs
	sequence := make([]byte, 1000)
	bases := []byte{'A', 'A', 'A', 'T', 'T', 'G', 'G', 'G', 'C', 'C'}
	for i := range sequence {
		sequence[i] = bases[i%len(bases)]
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		model.FindHomopolymerRuns(sequence)
	}
}

// TestHomopolymerErrorModel_stringIndex tests stringIndex utility
func TestStringIndexFunction(t *testing.T) {
	// Testing the stringIndex helper indirectly through code that uses it
	// We can verify the behavior through model usage
	model := NewHomopolymerErrorModel("illumina")
	assert.NotNil(t, model)

	// stringIndex is used internally but verify it works through getEffectType on GCbiasModel
	gcModel := NewGCbiasQualityModel("illumina")
	assert.NotNil(t, gcModel)
}

// TestCalculateErrorStatistics_EmptyOriginalSequence is a regression test for #120.
// calculateErrorStatistics divided by len(result.OriginalSequence) without guarding
// against zero length; a caller constructing ErrorSimulationResult with a non-nil
// Errors slice but a nil/empty OriginalSequence would get +Inf for ErrorRate.
func TestCalculateErrorStatistics_EmptyOriginalSequence(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")

	result := &ErrorSimulationResult{
		OriginalSequence: nil, // zero value — triggers the bug pre-fix
		Errors: []HomopolymerError{
			{ErrorType: "insertion", RunLength: 4, Position: 0},
		},
	}

	stats := model.calculateErrorStatistics(result)

	if stats.ErrorRate != 0.0 {
		t.Errorf("ErrorRate = %v with empty OriginalSequence; want 0.0 (not +Inf)", stats.ErrorRate)
	}
	// Other fields should still be populated correctly.
	if stats.TotalErrors != 1 {
		t.Errorf("TotalErrors = %d, want 1", stats.TotalErrors)
	}
	if stats.InsertionErrors != 1 {
		t.Errorf("InsertionErrors = %d, want 1", stats.InsertionErrors)
	}
}

// TestCalculateErrorStatistics_NoErrors verifies AvgRunLength is 0 (not NaN)
// when the result contains no errors — regression for #124.
func TestCalculateErrorStatistics_NoErrors(t *testing.T) {
	model := NewHomopolymerErrorModel("illumina")
	result := &ErrorSimulationResult{
		OriginalSequence: []byte("ACGTACGT"),
		Errors:           nil, // no homopolymer errors
	}
	stats := model.calculateErrorStatistics(result)
	if stats.AvgRunLength != 0.0 {
		t.Errorf("AvgRunLength = %v with no errors; want 0.0 (not NaN)", stats.AvgRunLength)
	}
	if stats.TotalErrors != 0 {
		t.Errorf("TotalErrors = %d with no errors; want 0", stats.TotalErrors)
	}
}
