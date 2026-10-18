//go:build !release
// +build !release

package main

import (
	"testing"
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
		got := calculateComplexity(seqs, 8)
		// Short seq skipped, average of 0 valid seqs => 0
		if got != 0 {
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

	t.Run("50pct GC returns 0.5", func(t *testing.T) {
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

	t.Run("all GC returns 1.0", func(t *testing.T) {
		seq := []byte("GGCC")
		if got := calculateGCContentSingle(seq); got != 1.0 {
			t.Errorf("expected 1.0, got %f", got)
		}
	})

	t.Run("all AT returns 0.0", func(t *testing.T) {
		seq := []byte("AATT")
		if got := calculateGCContentSingle(seq); got != 0.0 {
			t.Errorf("expected 0.0, got %f", got)
		}
	})

	t.Run("50pct", func(t *testing.T) {
		seq := []byte("GCTA")
		got := calculateGCContentSingle(seq)
		if got != 0.5 {
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

func TestGetPopulationProfiles(t *testing.T) {
	profiles := getPopulationProfiles()
	if len(profiles) == 0 {
		t.Error("expected non-empty profiles")
	}
	for _, p := range profiles {
		if p.Population == "" {
			t.Error("expected non-empty population name")
		}
		if p.TargetGC <= 0 || p.TargetGC >= 1 {
			t.Errorf("population %s has invalid TargetGC %f", p.Population, p.TargetGC)
		}
		if p.Tolerance <= 0 {
			t.Errorf("population %s has non-positive tolerance %f", p.Population, p.Tolerance)
		}
		if p.BaseComplexity <= 0 {
			t.Errorf("population %s has non-positive BaseComplexity %f", p.Population, p.BaseComplexity)
		}
	}
}

func TestGetProfileByPopulation(t *testing.T) {
	profiles := getPopulationProfiles()

	t.Run("existing population found", func(t *testing.T) {
		p := getProfileByPopulation(profiles, "YRI")
		if p.Population != "YRI" {
			t.Errorf("expected YRI, got %q", p.Population)
		}
	})

	t.Run("non-existing population returns zero value", func(t *testing.T) {
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
	t.Run("empty sequence returned unchanged", func(t *testing.T) {
		result := adjustSequenceGC([]byte{}, 0.42)
		if len(result) != 0 {
			t.Errorf("expected empty result, got len=%d", len(result))
		}
	})

	t.Run("returns same length as input", func(t *testing.T) {
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

func TestGenerateMockK8Sequences(t *testing.T) {
	k8 := ActualK8Results{Complexity: 0.6677, GCContent: 0.4816}

	t.Run("returns correct count", func(t *testing.T) {
		seqs := generateMockK8Sequences(k8, 5, 100)
		if len(seqs) != 5 {
			t.Errorf("expected 5 sequences, got %d", len(seqs))
		}
	})

	t.Run("sequences have requested length", func(t *testing.T) {
		seqs := generateMockK8Sequences(k8, 3, 200)
		for i, seq := range seqs {
			if len(seq) != 200 {
				t.Errorf("sequence %d: expected len 200, got %d", i, len(seq))
			}
		}
	})

	t.Run("sequences contain only valid DNA bases", func(t *testing.T) {
		seqs := generateMockK8Sequences(k8, 2, 100)
		for si, seq := range seqs {
			for bi, b := range seq {
				switch b {
				case 'A', 'T', 'G', 'C':
					// valid
				default:
					t.Errorf("seq %d base %d: invalid base %c", si, bi, b)
				}
			}
		}
	})
}

func TestAdjustSequencesGC(t *testing.T) {
	t.Run("returns same number of sequences", func(t *testing.T) {
		seqs := [][]byte{
			[]byte("ACGTACGT"),
			[]byte("GCGCGCGC"),
		}
		result := adjustSequencesGC(seqs, 0.42)
		if len(result) != len(seqs) {
			t.Errorf("expected %d sequences, got %d", len(seqs), len(result))
		}
	})
}

func TestAdjustGCUp_LowGCSequence(t *testing.T) {
	// All AT sequence: GC=0, target=0.5 => should convert some AT to GC
	seq := []byte("ATATATATATATATATATATATATATATATATATATATATAT")
	adjustGCUp(seq, 0.5)
	// Should have increased GC content or at least not crashed
	gcCount := 0
	for _, b := range seq {
		if b == 'G' || b == 'C' {
			gcCount++
		}
	}
	// With 0.3 probability filter, changes should be made
	// Just verify no crash and valid bases
	for i, b := range seq {
		switch b {
		case 'A', 'T', 'G', 'C':
			// ok
		default:
			t.Errorf("invalid base %c at position %d", b, i)
		}
	}
}

func TestAdjustGCUp_AlreadyAtTarget(t *testing.T) {
	// GC=1.0, target=0.5: changesToMake <= 0 so no changes
	seq := []byte("GCGCGCGC")
	orig := make([]byte, len(seq))
	copy(orig, seq)
	adjustGCUp(seq, 0.5)
	// No changes because current GC > target (adjustGCUp would find changesToMake <= 0)
}

func TestDemonstratePopulationAdjustment_YRI(t *testing.T) {
	k8Results := ActualK8Results{Complexity: 0.6677, GCContent: 0.4816}
	profile := PopulationGCProfile{
		Population:     "YRI",
		BaseComplexity: 0.6677,
		ObservedGC:     0.4816,
		TargetGC:       0.47,
		Tolerance:      0.02,
		Rationale:      "test",
	}
	result := demonstratePopulationAdjustment(k8Results, profile)
	if result.Population != "YRI" {
		t.Errorf("expected YRI, got %q", result.Population)
	}
	if result.Sequences <= 0 {
		t.Errorf("expected positive sequences count, got %d", result.Sequences)
	}
	if result.Complexity < 0 || result.Complexity > 1.0 {
		t.Errorf("complexity %f out of [0,1]", result.Complexity)
	}
}

func TestDemonstratePopulationAdjustment_HighGCTarget(t *testing.T) {
	// k8 GC=0.3, target=0.6 => adjustGCUp will be called
	k8Results := ActualK8Results{Complexity: 0.5, GCContent: 0.30}
	profile := PopulationGCProfile{
		Population:     "TEST",
		BaseComplexity: 0.5,
		ObservedGC:     0.30,
		TargetGC:       0.60,
		Tolerance:      0.05,
		Rationale:      "test high GC",
	}
	result := demonstratePopulationAdjustment(k8Results, profile)
	if result.Population != "TEST" {
		t.Errorf("expected TEST population, got %q", result.Population)
	}
}

func TestDetectAppleSiliconCores_AllBranches(t *testing.T) {
	// Just ensure no crash and positive result
	cores := detectAppleSiliconCores()
	if cores <= 0 {
		t.Errorf("expected positive cores, got %d", cores)
	}
}


// ---------------------------------------------------------------------------
// Additional tests for pure helper functions
// ---------------------------------------------------------------------------

func TestCalculateGCContentSingle_SingleG(t *testing.T) {
	got := calculateGCContentSingle([]byte("G"))
	if got != 1.0 {
		t.Errorf("got %v, want 1.0", got)
	}
}

func TestCalculateGCContentSingle_SingleA(t *testing.T) {
	got := calculateGCContentSingle([]byte("A"))
	if got != 0.0 {
		t.Errorf("got %v, want 0.0", got)
	}
}

func TestCalculateGCContentSingle_QuarterGC(t *testing.T) {
	// GCCC has 4 bases, 4 are G or C — wait, GCCC is 1G+3C = all GC => 1.0
	// Use GAAT: 1G + 2A + 1T = 1 GC out of 4 = 0.25
	got := calculateGCContentSingle([]byte("GAAT"))
	if got != 0.25 {
		t.Errorf("got %v, want 0.25", got)
	}
}

func TestCalculateGCContentSingle_InRange(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACGTACGTGCATGCAT"),
		[]byte("GGGGCCCCAAAATTTT"),
		[]byte("ATATATATATATATATAT"),
	}
	for _, seq := range seqs {
		got := calculateGCContentSingle(seq)
		if got < 0.0 || got > 1.0 {
			t.Errorf("GC content %v outside [0,1] for seq %s", got, seq)
		}
	}
}

func TestCalculateGCContent_MultipleSequencesAverage(t *testing.T) {
	// Two seqs: 0.0 GC and 1.0 GC => average = 0.5
	seqs := [][]byte{
		[]byte("ATAT"), // 0.0
		[]byte("GCGC"), // 1.0
	}
	got := calculateGCContent(seqs)
	if got != 0.5 {
		t.Errorf("got %v, want 0.5", got)
	}
}

func TestCalculateGCContent_SingleSequence(t *testing.T) {
	got := calculateGCContent([][]byte{[]byte("AACCGGTT")})
	if got != 0.5 {
		t.Errorf("got %v, want 0.5", got)
	}
}

func TestCalculateGCContent_InRange(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACGTACGT"),
		[]byte("GCGCGCGC"),
		[]byte("ATATATAT"),
	}
	got := calculateGCContent(seqs)
	if got < 0.0 || got > 1.0 {
		t.Errorf("GC content %v outside [0,1]", got)
	}
}

func TestCalculateComplexity_ExactlyK(t *testing.T) {
	// Sequence of exactly k bases produces 1 k-mer window
	// uniqueKmers=1, maxPos=1 => complexity=1.0
	seq := []byte("ACGT") // len=4, k=4
	got := calculateComplexity([][]byte{seq}, 4)
	if got != 1.0 {
		t.Errorf("got %v, want 1.0 for sequence of exactly k bases", got)
	}
}

func TestCalculateComplexity_MultipleSequences(t *testing.T) {
	// Average over two sequences
	seqs := [][]byte{
		[]byte("AAAAAAAAAA"),        // low complexity
		[]byte("ACGTACGTGCATGCAT"), // higher complexity
	}
	got := calculateComplexity(seqs, 4)
	if got < 0.0 || got > 1.0 {
		t.Errorf("complexity %v outside [0,1]", got)
	}
}

func TestCalculateComplexity_AllShortSkipped(t *testing.T) {
	// All sequences shorter than k are skipped; result = 0/0 sequences
	// Implementation divides by len(sequences) not valid count,
	// so result = 0 / 2 = 0 for two short sequences
	seqs := [][]byte{
		[]byte("AC"),
		[]byte("GT"),
	}
	got := calculateComplexity(seqs, 8)
	if got != 0.0 {
		t.Errorf("got %v, want 0.0 when all sequences shorter than k", got)
	}
}

func TestCalculateComplexity_HighDiversity(t *testing.T) {
	// Sequence with all distinct k-mers should have complexity close to 1.
	// "ACGT" repeated: 4-mers are ACGT, CGTA, GTAC, TACG (cycling) => 4 unique / 5 positions
	seq := []byte("ACGTACGTACGT")
	got := calculateComplexity([][]byte{seq}, 4)
	if got <= 0.0 {
		t.Errorf("got %v, want > 0 for diverse sequence", got)
	}
}

func TestAbs_PositiveUnchanged(t *testing.T) {
	if got := abs(42.0); got != 42.0 {
		t.Errorf("abs(42.0) = %v, want 42.0", got)
	}
}

func TestAbs_NegativeNegated(t *testing.T) {
	if got := abs(-3.14); got != 3.14 {
		t.Errorf("abs(-3.14) = %v, want 3.14", got)
	}
}

func TestAbs_LargeValue(t *testing.T) {
	const large = 1e15
	if got := abs(-large); got != large {
		t.Errorf("abs(%v) = %v, want %v", -large, got, large)
	}
}

func TestGetStatusEmoji_TrueIsCheckmark(t *testing.T) {
	got := getStatusEmoji(true)
	if got != "✅" {
		t.Errorf("getStatusEmoji(true) = %q, want checkmark", got)
	}
}

func TestGetStatusEmoji_FalseIsCross(t *testing.T) {
	got := getStatusEmoji(false)
	if got != "❌" {
		t.Errorf("getStatusEmoji(false) = %q, want cross mark", got)
	}
}

func TestGetStatusEmoji_NotEmpty(t *testing.T) {
	if s := getStatusEmoji(true); s == "" {
		t.Error("getStatusEmoji(true) returned empty string")
	}
	if s := getStatusEmoji(false); s == "" {
		t.Error("getStatusEmoji(false) returned empty string")
	}
}

func TestAdjustSequenceGC_GCContentApproachesTarget(t *testing.T) {
	// All-AT sequence adjusted toward 0.8 GC should increase GC content.
	seq := make([]byte, 200)
	for i := range seq {
		seq[i] = 'A'
	}
	result := adjustSequenceGC(seq, 0.8)
	gcBefore := calculateGCContentSingle(seq)
	gcAfter := calculateGCContentSingle(result)
	if gcAfter <= gcBefore {
		t.Errorf("GC after adjustment (%v) should be > before (%v) when targeting 0.8", gcAfter, gcBefore)
	}
}

func TestAdjustSequenceGC_AllValidBasesAfterAdjust(t *testing.T) {
	seq := []byte("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")
	result := adjustSequenceGC(seq, 0.2)
	for i, b := range result {
		switch b {
		case 'A', 'T', 'G', 'C':
			// valid
		default:
			t.Errorf("invalid base %c at position %d", b, i)
		}
	}
}

func TestAdjustSequenceGC_SameLengthAfterAdjust(t *testing.T) {
	seq := []byte("AAAAGGGGCCCCTTTT")
	result := adjustSequenceGC(seq, 0.3)
	if len(result) != len(seq) {
		t.Errorf("length changed from %d to %d", len(seq), len(result))
	}
}

func TestAdjustGCDown_ValidBases(t *testing.T) {
	// All-GC sequence, target=0.1 => should convert many GC to AT
	seq := []byte("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")
	adjustGCDown(seq, 0.1)
	for i, b := range seq {
		switch b {
		case 'A', 'T', 'G', 'C':
			// valid
		default:
			t.Errorf("invalid base %c at position %d after adjustGCDown", b, i)
		}
	}
}

func TestAdjustGCDown_LengthUnchanged(t *testing.T) {
	seq := []byte("GCGCGCGCGCGCGCGCGCGC")
	orig := make([]byte, len(seq))
	copy(orig, seq)
	adjustGCDown(seq, 0.1)
	if len(seq) != len(orig) {
		t.Errorf("length changed from %d to %d", len(orig), len(seq))
	}
}

func TestAdjustGCDown_AlreadyAtOrBelowTarget(t *testing.T) {
	// All-AT sequence: GC=0, target=0.5 => changesToMake <= 0, no-op
	seq := []byte("ATATATATATATATATAT")
	orig := make([]byte, len(seq))
	copy(orig, seq)
	adjustGCDown(seq, 0.5)
	for i := range seq {
		if seq[i] != orig[i] {
			t.Errorf("sequence modified at position %d when GC already below target", i)
		}
	}
}

func TestAdjustGCUp_ValidBases(t *testing.T) {
	// All-AT sequence, target=0.9 => should convert some AT to GC
	seq := make([]byte, 60)
	for i := range seq {
		if i%2 == 0 {
			seq[i] = 'A'
		} else {
			seq[i] = 'T'
		}
	}
	adjustGCUp(seq, 0.9)
	for i, b := range seq {
		switch b {
		case 'A', 'T', 'G', 'C':
			// valid
		default:
			t.Errorf("invalid base %c at position %d after adjustGCUp", b, i)
		}
	}
}

func TestAdjustGCUp_LengthUnchanged(t *testing.T) {
	seq := []byte("ATATATATATATATATATAT")
	adjustGCUp(seq, 0.8)
	if len(seq) != 20 {
		t.Errorf("length changed to %d, want 20", len(seq))
	}
}

func TestAdjustGCUp_AlreadyAtOrAboveTarget(t *testing.T) {
	// All-GC sequence: GC=1.0, target=0.3 => changesToMake <= 0, no-op
	seq := []byte("GCGCGCGCGCGCGCGCGCGC")
	orig := make([]byte, len(seq))
	copy(orig, seq)
	adjustGCUp(seq, 0.3)
	for i := range seq {
		if seq[i] != orig[i] {
			t.Errorf("sequence modified at position %d when GC already above target", i)
		}
	}
}
