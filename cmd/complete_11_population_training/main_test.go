//go:build !release
// +build !release

package main

import (
	"os"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func TestGetStatusEmoji_Success(t *testing.T) {
	if got := getStatusEmoji(true); got != "✅" {
		t.Errorf("getStatusEmoji(true) = %q, want ✅", got)
	}
}

func TestGetStatusEmoji_Failure(t *testing.T) {
	if got := getStatusEmoji(false); got != "❌" {
		t.Errorf("getStatusEmoji(false) = %q, want ❌", got)
	}
}

func TestAbs_Positive(t *testing.T) {
	if got := abs(5.0); got != 5.0 {
		t.Errorf("abs(5.0) = %f, want 5.0", got)
	}
}

func TestAbs_Negative(t *testing.T) {
	if got := abs(-5.0); got != 5.0 {
		t.Errorf("abs(-5.0) = %f, want 5.0", got)
	}
}

func TestAbs_Zero(t *testing.T) {
	if got := abs(0.0); got != 0.0 {
		t.Errorf("abs(0.0) = %f, want 0.0", got)
	}
}

func TestDetectAppleSiliconCores_Positive(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestCalculateGCInWindow_Empty(t *testing.T) {
	got := calculateGCInWindow([]byte{}, 0, 10)
	if got != 0.0 {
		t.Errorf("calculateGCInWindow(empty) = %f, want 0.0", got)
	}
}

func TestCalculateGCInWindow_AllGC(t *testing.T) {
	seq := []byte("GCGCGCGCGCGCGCGC")
	got := calculateGCInWindow(seq, 16, 10)
	if got < 0.99 || got > 1.01 {
		t.Errorf("calculateGCInWindow(all GC) = %f, want ~1.0", got)
	}
}

func TestCalculateGCInWindow_AllAT(t *testing.T) {
	seq := []byte("ATATATATATAT")
	got := calculateGCInWindow(seq, 12, 10)
	if got != 0.0 {
		t.Errorf("calculateGCInWindow(all AT) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	if got := calculateSequenceComplexity([]byte{}, 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(empty) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_TooShort(t *testing.T) {
	if got := calculateSequenceComplexity([]byte("AT"), 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(too short) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_Repetitive(t *testing.T) {
	seq := []byte("AAAAAAAAAAAAAAAA")
	got := calculateSequenceComplexity(seq, 4)
	if got >= 0.5 {
		t.Errorf("calculateSequenceComplexity(all-A) = %f, expected < 0.5", got)
	}
}

func TestCalculateSequenceComplexity_Diverse(t *testing.T) {
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCG")
	got := calculateSequenceComplexity(seq, 4)
	if got <= 0 {
		t.Errorf("calculateSequenceComplexity(diverse) = %f, expected > 0", got)
	}
}

func TestCalculateSequenceGCContent_Empty(t *testing.T) {
	if got := calculateSequenceGCContent([]byte{}); got != 0.0 {
		t.Errorf("calculateSequenceGCContent(empty) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceGCContent_AllGC(t *testing.T) {
	if got := calculateSequenceGCContent([]byte("GCGCGCGC")); got != 1.0 {
		t.Errorf("calculateSequenceGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestCalculateSequenceGCContent_AllAT(t *testing.T) {
	if got := calculateSequenceGCContent([]byte("ATATATAT")); got != 0.0 {
		t.Errorf("calculateSequenceGCContent(all AT) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceGCContent_HalfGC(t *testing.T) {
	got := calculateSequenceGCContent([]byte("GCATGCAT"))
	if got < 0.49 || got > 0.51 {
		t.Errorf("calculateSequenceGCContent(half GC) = %f, want ~0.5", got)
	}
}

func TestHashString_Deterministic(t *testing.T) {
	// Same input should produce same hash
	if hashString("YRI") != hashString("YRI") {
		t.Error("hashString is not deterministic")
	}
}

func TestHashString_DifferentStrings(t *testing.T) {
	// Different strings should (usually) produce different hashes
	if hashString("YRI") == hashString("CEU") {
		t.Error("hashString(YRI) == hashString(CEU), expected different values")
	}
}

func TestHashString_NonNegative(t *testing.T) {
	pops := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"}
	for _, pop := range pops {
		if hashString(pop) < 0 {
			t.Errorf("hashString(%q) = %d, want non-negative", pop, hashString(pop))
		}
	}
}

func TestGetComplete1000GenomesTargets_NotEmpty(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) == 0 {
		t.Error("getComplete1000GenomesTargets() returned empty slice")
	}
}

func TestGetComplete1000GenomesTargets_Has11Populations(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) != 11 {
		t.Errorf("expected 11 populations, got %d", len(targets))
	}
}

func TestGetComplete1000GenomesTargets_ValidGCTargets(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	for _, target := range targets {
		if target.TargetGC < 0 || target.TargetGC > 1 {
			t.Errorf("population %s has invalid GC target: %f", target.Population, target.TargetGC)
		}
	}
}

func TestGeneratePopulationSyntheticSequence_Length(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) == 0 {
		t.Skip("no targets available")
	}
	seq := generatePopulationSyntheticSequence(targets[0], 100)
	if len(seq) != 100 {
		t.Errorf("expected length 100, got %d", len(seq))
	}
}

func TestGeneratePopulationSyntheticSequence_ValidBases(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) == 0 {
		t.Skip("no targets available")
	}
	seq := generatePopulationSyntheticSequence(targets[0], 50)
	for i, base := range seq {
		if base != 'A' && base != 'T' && base != 'G' && base != 'C' {
			t.Errorf("invalid base at position %d: %c", i, base)
		}
	}
}

// ---------- hashString additional coverage (negative hash branch) ----------

func TestHashString_Empty_C11(t *testing.T) {
	got := hashString("")
	if got != 0 {
		t.Errorf("hashString(\"\") = %d, want 0", got)
	}
}

func TestHashString_LongString_NonNegative(t *testing.T) {
	// Long strings can produce large intermediate values; result must be >= 0
	longInputs := []string{
		"AFRIKAANS_POPULATION_LONG_NAME_EXCEEDS_HASH",
		"THIS_IS_A_VERY_LONG_STRING_THAT_MIGHT_OVERFLOW_AND_GO_NEGATIVE",
	}
	for _, s := range longInputs {
		got := hashString(s)
		if got < 0 {
			t.Errorf("hashString(%q) = %d, want >= 0", s, got)
		}
	}
}

func TestHashString_Distinct_C11(t *testing.T) {
	a := hashString("YRI")
	b := hashString("MXL")
	if a == b {
		t.Errorf("hashString(\"YRI\") == hashString(\"MXL\") = %d, expected different", a)
	}
}

// ---------- calculateGCInWindow additional coverage ----------

func TestCalculateGCInWindow_ZeroPos_C11(t *testing.T) {
	// pos=0, any windowSize → range [max(0-w,0)..0) is empty → returns 0.0
	seq := []byte("GCATAT")
	got := calculateGCInWindow(seq, 0, 4)
	if got != 0.0 {
		t.Errorf("calculateGCInWindow(pos=0) = %v, want 0.0", got)
	}
}

// ---------- calculateSequenceComplexity additional coverage ----------

func TestCalculateSequenceComplexity_ExactK(t *testing.T) {
	// Exactly k=6 chars → one k-mer total, 1 unique → complexity = 1.0
	got := calculateSequenceComplexity([]byte("ACGTAC"), 6)
	if got != 1.0 {
		t.Errorf("calculateSequenceComplexity(exact k) = %v, want 1.0", got)
	}
}

// ---------- displayPopulationOverview ----------

func TestDisplayPopulationOverview_NoCrash(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	displayPopulationOverview(targets)
}

func TestDisplayPopulationOverview_Empty(t *testing.T) {
	displayPopulationOverview([]PopulationTarget{})
}

// ---------- analyzePopulationResults ----------

func TestAnalyzePopulationResults_Empty(t *testing.T) {
	analyzePopulationResults(map[string]*PopulationModel{}, []PopulationTarget{})
}

func TestAnalyzePopulationResults_WithUntrained(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	models := map[string]*PopulationModel{
		targets[0].Population: {
			Population:  targets[0].Population,
			Model:       nil,
			TargetGC:    targets[0].TargetGC,
			Region:      targets[0].Region,
			Individuals: targets[0].Individuals,
			Trained:     false,
			Complexity:  0.0,
		},
	}
	analyzePopulationResults(models, targets[:1])
}

func TestAnalyzePopulationResults_WithTrained(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	models := map[string]*PopulationModel{
		targets[0].Population: {
			Population:  targets[0].Population,
			Model:       nil,
			TargetGC:    targets[0].TargetGC,
			Region:      targets[0].Region,
			Individuals: targets[0].Individuals,
			Trained:     true,
			Complexity:  0.75,
		},
	}
	analyzePopulationResults(models, targets[:1])
}

// ---------- analyzeRegionalResults ----------

func TestAnalyzeRegionalResults_Empty(t *testing.T) {
	analyzeRegionalResults(map[string]*PopulationModel{}, []PopulationTarget{})
}

func TestAnalyzeRegionalResults_WithData(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	models := map[string]*PopulationModel{}
	for _, target := range targets {
		models[target.Population] = &PopulationModel{
			Population:  target.Population,
			TargetGC:    target.TargetGC,
			Region:      target.Region,
			Individuals: target.Individuals,
			Trained:     true,
			Complexity:  0.70,
		}
	}
	analyzeRegionalResults(models, targets)
}

// ---------- generateValidationSequences ----------

func TestGenerateValidationSequences_EmptyModels(t *testing.T) {
	// No targets have models, so it should skip all and not crash
	targets := getComplete1000GenomesTargets()
	generateValidationSequences(map[string]*PopulationModel{}, targets)
}

func TestGenerateValidationSequences_AllUntrained(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	models := map[string]*PopulationModel{}
	for _, target := range targets {
		models[target.Population] = &PopulationModel{
			Population: target.Population,
			Trained:    false,
		}
	}
	generateValidationSequences(models, targets)
}

// ---------- writeFASTASequence ----------

func TestWriteFASTASequence_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	path := tmpDir + "/test.fasta"
	model := &PopulationModel{
		Population:  "YRI",
		Region:      "Africa",
		TargetGC:    0.42,
		Complexity:  0.75,
		Trained:     true,
		Individuals: 108,
	}
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCGATCG")
	err := writeFASTASequence(path, "YRI", 32, 1, seq, model)
	if err != nil {
		t.Fatalf("writeFASTASequence returned error: %v", err)
	}
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("expected file to exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("expected non-empty FASTA file")
	}
}

func TestWriteFASTASequence_BadPath(t *testing.T) {
	model := &PopulationModel{Population: "YRI", Region: "Africa"}
	err := writeFASTASequence("/nonexistent/dir/test.fasta", "YRI", 32, 1, []byte("ATCG"), model)
	if err == nil {
		t.Error("expected error for bad path")
	}
}

// ---------- trainModelWithSyntheticData ----------

func TestTrainModelWithSyntheticData_NoCrash(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) == 0 {
		t.Skip("no targets")
	}
	model := authenticity.NewAdaptiveHMMModel(4) // small k for speed
	err := trainModelWithSyntheticData(model, targets[0])
	// May succeed or fail depending on model state, but shouldn't panic
	_ = err
}


// ============================================================
// Additional coverage for pure helper functions
// ============================================================

// ---------- abs (additional edge cases) ----------

func TestAbs_LargePositive_C11(t *testing.T) {
	if got := abs(1e15); got != 1e15 {
		t.Errorf("abs(1e15) = %v, want 1e15", got)
	}
}

func TestAbs_LargeNegative_C11(t *testing.T) {
	if got := abs(-1e15); got != 1e15 {
		t.Errorf("abs(-1e15) = %v, want 1e15", got)
	}
}

func TestAbs_SmallFraction_C11(t *testing.T) {
	if got := abs(-0.0001); got != 0.0001 {
		t.Errorf("abs(-0.0001) = %v, want 0.0001", got)
	}
}

// ---------- getStatusEmoji (additional cases) ----------

func TestGetStatusEmoji_ReturnsDifferentValues_C11(t *testing.T) {
	if getStatusEmoji(true) == getStatusEmoji(false) {
		t.Error("getStatusEmoji(true) and getStatusEmoji(false) must return different strings")
	}
}

func TestGetStatusEmoji_NonEmpty_C11(t *testing.T) {
	if getStatusEmoji(true) == "" {
		t.Error("getStatusEmoji(true) must not be empty")
	}
	if getStatusEmoji(false) == "" {
		t.Error("getStatusEmoji(false) must not be empty")
	}
}

// ---------- hashString (additional cases) ----------

func TestHashString_SingleChar_C11(t *testing.T) {
	// Single-char input: hash = 0*31 + int('A') = int('A')
	got := hashString("A")
	if got != int('A') {
		t.Errorf("hashString(\"A\") = %d, want %d", got, int('A'))
	}
}

func TestHashString_Reproducible_C11(t *testing.T) {
	input := "CEU"
	first := hashString(input)
	second := hashString(input)
	if first != second {
		t.Errorf("hashString(%q) not reproducible: %d vs %d", input, first, second)
	}
}

func TestHashString_AllPopulations_Unique_C11(t *testing.T) {
	pops := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"}
	seen := make(map[int]string)
	for _, pop := range pops {
		h := hashString(pop)
		if prev, exists := seen[h]; exists {
			t.Errorf("hashString(%q) == hashString(%q) = %d (collision)", pop, prev, h)
		}
		seen[h] = pop
	}
}

// ---------- calculateGCInWindow (additional cases) ----------

func TestCalculateGCInWindow_WindowLargerThanSequence_C11(t *testing.T) {
	// windowSize > len: start clamped to 0, entire sequence counted
	seq := []byte("GCGC") // 100% GC
	got := calculateGCInWindow(seq, 4, 100)
	if got < 0.99 || got > 1.01 {
		t.Errorf("calculateGCInWindow(window>len) = %v, want ~1.0", got)
	}
}

func TestCalculateGCInWindow_MixedWindow_C11(t *testing.T) {
	// seq = "ATGCATGC", pos=4, windowSize=4 → window [0,4) = "ATGC" → 2/4 = 0.5
	seq := []byte("ATGCATGC")
	got := calculateGCInWindow(seq, 4, 4)
	if got < 0.49 || got > 0.51 {
		t.Errorf("calculateGCInWindow(mixed) = %v, want ~0.5", got)
	}
}

func TestCalculateGCInWindow_PartialWindowAtStart_C11(t *testing.T) {
	// pos=3, windowSize=10: start = max(3-10, 0) = 0, counts seq[0:3] = "GCA"
	// G=GC, C=GC, A=not → gcCount=2, totalBases=3 → 2/3 ≈ 0.666
	seq := []byte("GCATAT")
	got := calculateGCInWindow(seq, 3, 10)
	if got < 0.65 || got > 0.68 {
		t.Errorf("calculateGCInWindow(partial start) = %v, want ~0.666", got)
	}
}

func TestCalculateGCInWindow_ResultInRange_C11(t *testing.T) {
	seq := []byte("AAATTTGGGCCC")
	for pos := 1; pos <= len(seq); pos++ {
		got := calculateGCInWindow(seq, pos, 4)
		if got < 0 || got > 1 {
			t.Errorf("calculateGCInWindow(pos=%d) = %v, out of [0,1]", pos, got)
		}
	}
}

// ---------- calculateSequenceComplexity (additional cases) ----------

func TestCalculateSequenceComplexity_K1_AllSame_C11(t *testing.T) {
	// k=1, all 'A': 1 unique kmer out of 8 total → 1/8
	seq := []byte("AAAAAAAA")
	got := calculateSequenceComplexity(seq, 1)
	want := 1.0 / float64(len(seq))
	if got < want-1e-9 || got > want+1e-9 {
		t.Errorf("calculateSequenceComplexity(all-A, k=1) = %v, want %v", got, want)
	}
}

func TestCalculateSequenceComplexity_K1_FourBases_C11(t *testing.T) {
	// k=1, 4 distinct bases in 4-char seq → 4/4 = 1.0
	got := calculateSequenceComplexity([]byte("ATGC"), 1)
	if got < 0.99 || got > 1.01 {
		t.Errorf("calculateSequenceComplexity(ATGC, k=1) = %v, want 1.0", got)
	}
}

func TestCalculateSequenceComplexity_BoundedByOne_C11(t *testing.T) {
	seq := []byte("ACGTACGTACGT")
	for k := 1; k <= 6; k++ {
		got := calculateSequenceComplexity(seq, k)
		if got < 0 || got > 1.0 {
			t.Errorf("calculateSequenceComplexity(k=%d) = %v, out of [0,1]", k, got)
		}
	}
}

// ---------- calculateSequenceGCContent (additional cases) ----------

func TestCalculateSequenceGCContent_SingleG_C11(t *testing.T) {
	if got := calculateSequenceGCContent([]byte("G")); got != 1.0 {
		t.Errorf("calculateSequenceGCContent(\"G\") = %v, want 1.0", got)
	}
}

func TestCalculateSequenceGCContent_SingleA_C11(t *testing.T) {
	if got := calculateSequenceGCContent([]byte("A")); got != 0.0 {
		t.Errorf("calculateSequenceGCContent(\"A\") = %v, want 0.0", got)
	}
}

func TestCalculateSequenceGCContent_SingleC_C11(t *testing.T) {
	if got := calculateSequenceGCContent([]byte("C")); got != 1.0 {
		t.Errorf("calculateSequenceGCContent(\"C\") = %v, want 1.0", got)
	}
}

func TestCalculateSequenceGCContent_SingleT_C11(t *testing.T) {
	if got := calculateSequenceGCContent([]byte("T")); got != 0.0 {
		t.Errorf("calculateSequenceGCContent(\"T\") = %v, want 0.0", got)
	}
}

func TestCalculateSequenceGCContent_ResultInRange_C11(t *testing.T) {
	seqs := [][]byte{
		[]byte("ATCG"),
		[]byte("GGGGGGGG"),
		[]byte("AAAAAAAA"),
		[]byte("GCGCGCGCATATATAT"),
	}
	for _, seq := range seqs {
		got := calculateSequenceGCContent(seq)
		if got < 0 || got > 1 {
			t.Errorf("calculateSequenceGCContent(%q) = %v, out of [0,1]", seq, got)
		}
	}
}

// ---------- generatePopulationSyntheticSequence (additional cases) ----------

func TestGeneratePopulationSyntheticSequence_ZeroLength_C11(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) == 0 {
		t.Skip("no targets available")
	}
	seq := generatePopulationSyntheticSequence(targets[0], 0)
	if len(seq) != 0 {
		t.Errorf("expected empty sequence for length=0, got %d bases", len(seq))
	}
}

func TestGeneratePopulationSyntheticSequence_Reproducible_C11(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) == 0 {
		t.Skip("no targets available")
	}
	seq1 := generatePopulationSyntheticSequence(targets[0], 200)
	seq2 := generatePopulationSyntheticSequence(targets[0], 200)
	if string(seq1) != string(seq2) {
		t.Error("generatePopulationSyntheticSequence not reproducible for same population")
	}
}

func TestGeneratePopulationSyntheticSequence_DifferentPopulationsDiffer_C11(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) < 2 {
		t.Skip("need at least 2 targets")
	}
	seq1 := generatePopulationSyntheticSequence(targets[0], 500)
	seq2 := generatePopulationSyntheticSequence(targets[1], 500)
	if string(seq1) == string(seq2) {
		t.Error("different populations produced identical sequences")
	}
}

func TestGeneratePopulationSyntheticSequence_GCApproachesTarget_C11(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) == 0 {
		t.Skip("no targets available")
	}
	for _, target := range targets {
		seq := generatePopulationSyntheticSequence(target, 5000)
		gc := calculateSequenceGCContent(seq)
		tolerance := target.Tolerance
		if tolerance == 0 {
			tolerance = 0.05
		}
		// Allow generous band: target ± 3*tolerance
		lo := target.TargetGC - 3*tolerance
		hi := target.TargetGC + 3*tolerance
		if gc < lo || gc > hi {
			t.Errorf("population %s: GC=%.4f not within [%.4f, %.4f]",
				target.Population, gc, lo, hi)
		}
	}
}

func TestGeneratePopulationSyntheticSequence_AllPopulationsValidBases_C11(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) == 0 {
		t.Skip("no targets available")
	}
	for _, target := range targets {
		seq := generatePopulationSyntheticSequence(target, 300)
		for i, b := range seq {
			if b != 'A' && b != 'T' && b != 'G' && b != 'C' {
				t.Errorf("population %s: invalid base %c at position %d",
					target.Population, b, i)
				break
			}
		}
	}
}

func TestGeneratePopulationSyntheticSequence_LargeLength_C11(t *testing.T) {
	targets := getComplete1000GenomesTargets()
	if len(targets) == 0 {
		t.Skip("no targets available")
	}
	const bigLen = 10000
	seq := generatePopulationSyntheticSequence(targets[0], bigLen)
	if len(seq) != bigLen {
		t.Errorf("expected length %d, got %d", bigLen, len(seq))
	}
}
