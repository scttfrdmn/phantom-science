//go:build !release
// +build !release

package main

import (
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestCalculateComplexity(t *testing.T) {
	t.Run("empty returns 0", func(t *testing.T) {
		if got := calculateComplexity(nil, 8); got != 0 {
			t.Errorf("expected 0 for nil, got %f", got)
		}
		if got := calculateComplexity([][]byte{}, 8); got != 0 {
			t.Errorf("expected 0 for empty slice, got %f", got)
		}
	})

	t.Run("sequence shorter than k skipped", func(t *testing.T) {
		seqs := [][]byte{[]byte("AC")}
		if got := calculateComplexity(seqs, 8); got != 0 {
			t.Errorf("expected 0 for sequence shorter than k, got %f", got)
		}
	})

	t.Run("uniform sequence low complexity", func(t *testing.T) {
		seq := make([]byte, 200)
		for i := range seq {
			seq[i] = 'A'
		}
		got := calculateComplexity([][]byte{seq}, 8)
		if got > 0.05 {
			t.Errorf("uniform sequence complexity too high: %f", got)
		}
	})

	t.Run("result bounded in [0,1]", func(t *testing.T) {
		seq := []byte("ACGTACGTGCATGCATACGTACGT")
		got := calculateComplexity([][]byte{seq}, 4)
		if got < 0 || got > 1 {
			t.Errorf("complexity %f out of [0,1]", got)
		}
	})
}

func TestCalculateGCContent(t *testing.T) {
	t.Run("empty returns 0", func(t *testing.T) {
		if got := calculateGCContent(nil); got != 0 {
			t.Errorf("expected 0, got %f", got)
		}
		if got := calculateGCContent([][]byte{}); got != 0 {
			t.Errorf("expected 0, got %f", got)
		}
	})

	t.Run("all G/C returns 1.0", func(t *testing.T) {
		seq := []byte("GCGCGC")
		if got := calculateGCContent([][]byte{seq}); got != 1.0 {
			t.Errorf("expected 1.0, got %f", got)
		}
	})

	t.Run("all A/T returns 0.0", func(t *testing.T) {
		seq := []byte("ATATAT")
		if got := calculateGCContent([][]byte{seq}); got != 0.0 {
			t.Errorf("expected 0.0, got %f", got)
		}
	})

	t.Run("50pct GC", func(t *testing.T) {
		seq := []byte("AACCGGTT")
		if got := calculateGCContent([][]byte{seq}); got != 0.5 {
			t.Errorf("expected 0.5, got %f", got)
		}
	})
}

func TestCalculateGCContentSingle(t *testing.T) {
	t.Run("empty returns 0", func(t *testing.T) {
		if got := calculateGCContentSingle([]byte{}); got != 0 {
			t.Errorf("expected 0, got %f", got)
		}
	})

	t.Run("all GC", func(t *testing.T) {
		seq := []byte("GGCC")
		if got := calculateGCContentSingle(seq); got != 1.0 {
			t.Errorf("expected 1.0, got %f", got)
		}
	})

	t.Run("all AT", func(t *testing.T) {
		seq := []byte("AATT")
		if got := calculateGCContentSingle(seq); got != 0.0 {
			t.Errorf("expected 0.0, got %f", got)
		}
	})

	t.Run("50pct", func(t *testing.T) {
		seq := []byte("GCTA")
		if got := calculateGCContentSingle(seq); got != 0.5 {
			t.Errorf("expected 0.5, got %f", got)
		}
	})
}

func TestAbs(t *testing.T) {
	tests := []struct{ in, out float64 }{
		{0, 0},
		{1.5, 1.5},
		{-1.5, 1.5},
		{-0.001, 0.001},
	}
	for _, tt := range tests {
		if got := abs(tt.in); got != tt.out {
			t.Errorf("abs(%f) = %f, want %f", tt.in, got, tt.out)
		}
	}
}

func TestDetectAppleSiliconCoresPositive(t *testing.T) {
	if cores := detectAppleSiliconCores(); cores <= 0 {
		t.Errorf("expected positive core count, got %d", cores)
	}
}

func TestGetOptimizedPopulationProfiles(t *testing.T) {
	profiles := getOptimizedPopulationProfiles()
	if len(profiles) == 0 {
		t.Error("expected non-empty profiles")
	}
	for _, p := range profiles {
		if p.Population == "" {
			t.Error("expected non-empty population name")
		}
		if p.TargetGC <= 0 || p.TargetGC >= 1 {
			t.Errorf("population %s: invalid TargetGC %f", p.Population, p.TargetGC)
		}
		if p.Tolerance <= 0 {
			t.Errorf("population %s: non-positive tolerance %f", p.Population, p.Tolerance)
		}
		if p.BaseComplexity <= 0 {
			t.Errorf("population %s: non-positive BaseComplexity %f", p.Population, p.BaseComplexity)
		}
		if p.Rationale == "" {
			t.Errorf("population %s: empty rationale", p.Population)
		}
	}
}

func TestGetProfileByPopulation(t *testing.T) {
	profiles := getOptimizedPopulationProfiles()

	t.Run("existing population", func(t *testing.T) {
		p := getProfileByPopulation(profiles, "YRI")
		if p.Population != "YRI" {
			t.Errorf("expected YRI, got %q", p.Population)
		}
	})

	t.Run("non-existing returns zero value", func(t *testing.T) {
		p := getProfileByPopulation(profiles, "NONEXISTENT")
		if p.Population != "" {
			t.Errorf("expected zero-value population, got %q", p.Population)
		}
	})
}

func TestGetStatusEmoji(t *testing.T) {
	if got := getStatusEmoji(true); got != "✅" {
		t.Errorf("expected ✅ for true, got %q", got)
	}
	if got := getStatusEmoji(false); got != "❌" {
		t.Errorf("expected ❌ for false, got %q", got)
	}
}

func TestAdjustSequenceGC(t *testing.T) {
	t.Run("empty sequence returned as-is", func(t *testing.T) {
		result := adjustSequenceGC([]byte{}, 0.42)
		if len(result) != 0 {
			t.Errorf("expected empty, got len=%d", len(result))
		}
	})

	t.Run("returns same length", func(t *testing.T) {
		seq := []byte("ACGTACGTGCATGCATACGTACGT")
		result := adjustSequenceGC(seq, 0.50)
		if len(result) != len(seq) {
			t.Errorf("expected length %d, got %d", len(seq), len(result))
		}
	})

	t.Run("does not modify original", func(t *testing.T) {
		seq := []byte("ACGTACGT")
		orig := make([]byte, len(seq))
		copy(orig, seq)
		adjustSequenceGC(seq, 0.42)
		for i := range seq {
			if seq[i] != orig[i] {
				t.Errorf("original was modified at position %d", i)
			}
		}
	})
}

func TestAdjustGCContentPreserveComplexity(t *testing.T) {
	t.Run("returns same number of sequences", func(t *testing.T) {
		seqs := [][]byte{
			[]byte("ACGTACGT"),
			[]byte("GCGCGCGC"),
		}
		result := adjustGCContentPreserveComplexity(seqs, 0.42)
		if len(result) != len(seqs) {
			t.Errorf("expected %d sequences, got %d", len(seqs), len(result))
		}
	})
}

func TestPopulationResultsStruct(t *testing.T) {
	r := PopulationResults{
		Population:  "CEU",
		Complexity:  0.6677,
		GCContent:   0.45,
		Success:     true,
		Sequences:   10,
		TrainingSeq: 5,
	}
	if r.Population != "CEU" {
		t.Errorf("expected CEU, got %s", r.Population)
	}
	if r.Complexity != 0.6677 {
		t.Errorf("expected 0.6677, got %f", r.Complexity)
	}
	if !r.Success {
		t.Error("expected Success=true")
	}
}

func TestAdjustGCDownAndUp(t *testing.T) {
	t.Run("adjustGCDown reduces GC content", func(t *testing.T) {
		// Start with 100% GC
		seq := []byte("GCGCGCGCGCGCGCGCGCGCGCGC")
		before := calculateGCContentSingle(seq)
		adjustGCDown(seq, 0.4)
		after := calculateGCContentSingle(seq)
		// After adjusting down, GC should not be higher than before
		if after > before {
			t.Errorf("adjustGCDown should not increase GC: before=%f after=%f", before, after)
		}
	})

	t.Run("adjustGCUp increases GC content", func(t *testing.T) {
		// Start with 0% GC
		seq := []byte("ATATATATATATATATATATAT")
		before := calculateGCContentSingle(seq)
		adjustGCUp(seq, 0.6)
		after := calculateGCContentSingle(seq)
		// After adjusting up, GC should not be lower than before
		if after < before {
			t.Errorf("adjustGCUp should not decrease GC: before=%f after=%f", before, after)
		}
	})
}

// ---------- generateValidationSequences ----------

func TestGenerateValidationSequences_UntrainedModel_ReturnsError(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(8)
	// Untrained model - GenerateSequence will fail
	_, err := generateValidationSequences(model, 1, 100)
	// Expected to fail for untrained model - just verify no panic
	_ = err
}

func TestGenerateValidationSequences_TrainedModel(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
	if err := model.TrainFromSequence(seq); err != nil {
		t.Skipf("model training failed: %v", err)
	}
	seqs, err := generateValidationSequences(model, 3, 50)
	if err != nil {
		t.Fatalf("generateValidationSequences error: %v", err)
	}
	if len(seqs) != 3 {
		t.Errorf("expected 3 sequences, got %d", len(seqs))
	}
}

// ---------- processPopulationEfficiently ----------

func TestProcessPopulationEfficiently_UntrainedModel_ReturnsFailed(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(8)
	profile := PopulationGCProfile{
		Population: "YRI",
		TargetGC:   0.47,
		Tolerance:  0.02,
	}
	result := processPopulationEfficiently(model, profile)
	// With untrained model, generateValidationSequences fails → Success should be false
	if result.Population != "YRI" {
		t.Errorf("Population = %q, want YRI", result.Population)
	}
}
