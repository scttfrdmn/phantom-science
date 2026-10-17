//go:build !release

package authenticity

import (
	"math"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/population"
)

func newFOS() *FingerprintObfuscationSystem {
	return NewFingerprintObfuscationSystem()
}

func makeSeqMap(seqs []string) map[population.Population][]string {
	return map[population.Population][]string{population.CEU: seqs}
}

// ---------------------------------------------------------------------------
// calculateAverageSubstitutionRate
// ---------------------------------------------------------------------------

func TestCalculateAverageSubstitutionRate_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateAverageSubstitutionRate(makeSeqMap(nil))
	if got <= 0 {
		t.Errorf("empty input: got %v, want > 0 (biological floor)", got)
	}
}

func TestCalculateAverageSubstitutionRate_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ACGTACGTACGT", "GGGGCCCCAAAATTTT"})
	a := fos.calculateAverageSubstitutionRate(seqs)
	b := fos.calculateAverageSubstitutionRate(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// Sequences with extreme GC (all GC) deviate more from 42% than balanced sequences.
func TestCalculateAverageSubstitutionRate_HighGCHigherRate(t *testing.T) {
	fos := newFOS()
	balanced := makeSeqMap([]string{"ACGTACGTACGTACGT"}) // GC ~= 50%
	highGC := makeSeqMap([]string{"GCGCGCGCGCGCGCGC"})  // GC = 100%
	rateBalanced := fos.calculateAverageSubstitutionRate(balanced)
	rateHigh := fos.calculateAverageSubstitutionRate(highGC)
	if rateHigh <= rateBalanced {
		t.Errorf("100%% GC rate=%v should exceed balanced GC rate=%v", rateHigh, rateBalanced)
	}
}

// ---------------------------------------------------------------------------
// calculateAverageIndelRate
// ---------------------------------------------------------------------------

func TestCalculateAverageIndelRate_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateAverageIndelRate(makeSeqMap(nil))
	if got <= 0 {
		t.Errorf("empty input: got %v, want > 0 (biological floor)", got)
	}
}

func TestCalculateAverageIndelRate_NoHomopolymers(t *testing.T) {
	fos := newFOS()
	// Alternating bases — no run ≥ 5.
	seqs := makeSeqMap([]string{"ACGTACGTACGT", "TGCATGCA"})
	got := fos.calculateAverageIndelRate(seqs)
	if got != 0.0 {
		t.Errorf("no homopolymers ≥ 5: got %v, want 0.0", got)
	}
}

func TestCalculateAverageIndelRate_LongHomopolymer(t *testing.T) {
	fos := newFOS()
	// 10-base homopolymer in a 12-base sequence → fraction = 10/12.
	seqs := makeSeqMap([]string{"AAAAAAAAAA" + "GC"})
	got := fos.calculateAverageIndelRate(seqs)
	want := 10.0 / 12.0
	if math.Abs(got-want) > 1e-9 {
		t.Errorf("10/12 homopolymer: got %v, want %v", got, want)
	}
}

func TestCalculateAverageIndelRate_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"AAAAAACGT", "TTTTTGCA"})
	a := fos.calculateAverageIndelRate(seqs)
	b := fos.calculateAverageIndelRate(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculateContextDependency
// ---------------------------------------------------------------------------

func TestCalculateContextDependency_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateContextDependency(makeSeqMap(nil))
	if got <= 0 {
		t.Errorf("empty input: got %v, want > 0", got)
	}
}

// A sequence with no CpG dinucleotides should give O/E = 0.
func TestCalculateContextDependency_NoCpG(t *testing.T) {
	fos := newFOS()
	// Only GC, no CpG (no C followed by G).
	seqs := makeSeqMap([]string{"GGGGGCCCCC"}) // GG, GC, CC patterns — no CG
	got := fos.calculateContextDependency(seqs)
	if got != 0.0 {
		t.Errorf("no CpG: got %v, want 0.0", got)
	}
}

// A purely CG-repeat sequence has very high O/E.
func TestCalculateContextDependency_AllCpG(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"CGCGCGCGCG"}) // all CpG
	got := fos.calculateContextDependency(seqs)
	if got <= 1.0 {
		t.Errorf("all-CpG: got %v, want > 1.0 (enriched)", got)
	}
}

func TestCalculateContextDependency_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ACGTCGACGT", "GCGTATCG"})
	a := fos.calculateContextDependency(seqs)
	b := fos.calculateContextDependency(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculateQualityCorrelation
// ---------------------------------------------------------------------------

func TestCalculateQualityCorrelation_EmptyQualities(t *testing.T) {
	fos := newFOS()
	got := fos.calculateQualityCorrelation(makeSeqMap(nil), nil)
	if got < -1 || got > 1 {
		t.Errorf("empty: got %v outside [-1,1]", got)
	}
}

// Decreasing quality produces a negative correlation.
func TestCalculateQualityCorrelation_DecreasingQuality(t *testing.T) {
	fos := newFOS()
	quals := [][]int{{40, 35, 30, 25, 20, 15, 10, 5}}
	got := fos.calculateQualityCorrelation(makeSeqMap(nil), quals)
	if got >= 0 {
		t.Errorf("decreasing quality: got %v, want < 0", got)
	}
}

// Increasing quality produces a positive correlation.
func TestCalculateQualityCorrelation_IncreasingQuality(t *testing.T) {
	fos := newFOS()
	quals := [][]int{{5, 10, 15, 20, 25, 30, 35, 40}}
	got := fos.calculateQualityCorrelation(makeSeqMap(nil), quals)
	if got <= 0 {
		t.Errorf("increasing quality: got %v, want > 0", got)
	}
}

// Constant quality → correlation is 0 (or undefined → clamped to 0).
func TestCalculateQualityCorrelation_ConstantQuality(t *testing.T) {
	fos := newFOS()
	quals := [][]int{{30, 30, 30, 30, 30}}
	got := fos.calculateQualityCorrelation(makeSeqMap(nil), quals)
	if got != 0.0 {
		t.Errorf("constant quality: got %v, want 0.0", got)
	}
}

func TestCalculateQualityCorrelation_InRange(t *testing.T) {
	fos := newFOS()
	quals := [][]int{{40, 38, 35, 30, 25}, {20, 22, 28, 35, 38}}
	got := fos.calculateQualityCorrelation(makeSeqMap(nil), quals)
	if got < -1 || got > 1 {
		t.Errorf("mixed quality: got %v outside [-1,1]", got)
	}
}

// ---------------------------------------------------------------------------
// calculateHomopolymerPatterns
// ---------------------------------------------------------------------------

func TestCalculateHomopolymerPatterns_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateHomopolymerPatterns(makeSeqMap(nil))
	if got <= 0 {
		t.Errorf("empty: got %v, want > 0 (biological floor)", got)
	}
}

func TestCalculateHomopolymerPatterns_NoLongRuns(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ACGTACGTACGT"})
	got := fos.calculateHomopolymerPatterns(seqs)
	if got != 0.0 {
		t.Errorf("no runs ≥ 4: got %v, want 0.0", got)
	}
}

func TestCalculateHomopolymerPatterns_AllSameBase(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"AAAAAAAAAA"}) // 10-base run ≥ 4
	got := fos.calculateHomopolymerPatterns(seqs)
	if got != 1.0 {
		t.Errorf("all-A: got %v, want 1.0", got)
	}
}

func TestCalculateHomopolymerPatterns_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"AAAACGTTTTT", "GCGCAAAA"})
	a := fos.calculateHomopolymerPatterns(seqs)
	b := fos.calculateHomopolymerPatterns(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

func TestCalculateHomopolymerPatterns_InRange(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"AAACGT", "GGGGACGT", "ACGTTTTT"})
	got := fos.calculateHomopolymerPatterns(seqs)
	if got < 0 || got > 1 {
		t.Errorf("got %v outside [0,1]", got)
	}
}

// ---------------------------------------------------------------------------
// calculateCycleDecay
// ---------------------------------------------------------------------------

func TestCalculateCycleDecay_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateCycleDecay(nil)
	if got < 0 || got > 1 {
		t.Errorf("empty: got %v outside [0,1]", got)
	}
}

func TestCalculateCycleDecay_FlatQuality(t *testing.T) {
	fos := newFOS()
	// Constant quality: no decay.
	quals := [][]int{{30, 30, 30, 30, 30, 30, 30, 30, 30, 30}}
	got := fos.calculateCycleDecay(quals)
	if got != 0.0 {
		t.Errorf("flat quality: got %v, want 0.0", got)
	}
}

func TestCalculateCycleDecay_StrongDecay(t *testing.T) {
	fos := newFOS()
	// Quality drops from 40 down to 1 across 10 cycles.
	quals := [][]int{{40, 36, 32, 28, 24, 20, 16, 12, 8, 4}}
	got := fos.calculateCycleDecay(quals)
	if got <= 0 {
		t.Errorf("strong decay: got %v, want > 0", got)
	}
	if got > 1 {
		t.Errorf("strong decay: got %v, want <= 1", got)
	}
}

func TestCalculateCycleDecay_Deterministic(t *testing.T) {
	fos := newFOS()
	quals := [][]int{{40, 38, 35, 30, 25, 20}}
	a := fos.calculateCycleDecay(quals)
	b := fos.calculateCycleDecay(quals)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculateDistributionShape
// ---------------------------------------------------------------------------

func TestCalculateDistributionShape_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateDistributionShape(nil)
	if got < 0 || got > 1 {
		t.Errorf("empty: got %v outside [0,1]", got)
	}
}

func TestCalculateDistributionShape_Deterministic(t *testing.T) {
	fos := newFOS()
	quals := [][]int{{30, 35, 20, 40, 25}}
	a := fos.calculateDistributionShape(quals)
	b := fos.calculateDistributionShape(quals)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// With constant quality there is no spread so all values are at mean — fraction
// outside [mean-σ, mean+σ] = 0 (σ = 0 so the interval is a point).
func TestCalculateDistributionShape_ConstantQuality(t *testing.T) {
	fos := newFOS()
	quals := [][]int{{30, 30, 30, 30, 30}}
	got := fos.calculateDistributionShape(quals)
	if got != 0.0 {
		t.Errorf("constant quality: got %v, want 0.0", got)
	}
}

// ---------------------------------------------------------------------------
// calculateChemistrySignature
// ---------------------------------------------------------------------------

func TestCalculateChemistrySignature_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateChemistrySignature(makeSeqMap(nil))
	if got < 0 || got > 1 {
		t.Errorf("empty: got %v outside [0,1]", got)
	}
}

func TestCalculateChemistrySignature_NoExtremeGC(t *testing.T) {
	fos := newFOS()
	// GC ≈ 50% — within [30%, 70%].
	seqs := makeSeqMap([]string{"ACGTACGTACGT", "GCATGCAT"})
	got := fos.calculateChemistrySignature(seqs)
	if got != 0.0 {
		t.Errorf("no extreme GC: got %v, want 0.0", got)
	}
}

func TestCalculateChemistrySignature_AllExtremeGC(t *testing.T) {
	fos := newFOS()
	// 100% GC — all reads extreme.
	seqs := makeSeqMap([]string{"GCGCGCGCGCGC", "GCGCGCGCGCGC"})
	got := fos.calculateChemistrySignature(seqs)
	if got != 1.0 {
		t.Errorf("all extreme GC: got %v, want 1.0", got)
	}
}

func TestCalculateChemistrySignature_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ACGTACGT", "GCGCGCGC", "TTTTTTTT"})
	a := fos.calculateChemistrySignature(seqs)
	b := fos.calculateChemistrySignature(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculateReadLengthPattern
// ---------------------------------------------------------------------------

func TestCalculateReadLengthPattern_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateReadLengthPattern(makeSeqMap(nil))
	if got != 0.0 {
		t.Errorf("empty: got %v, want 0.0", got)
	}
}

func TestCalculateReadLengthPattern_UniformLength(t *testing.T) {
	fos := newFOS()
	// All reads same length → CV = 0.
	seqs := makeSeqMap([]string{"ACGT", "ACGT", "ACGT"})
	got := fos.calculateReadLengthPattern(seqs)
	if got != 0.0 {
		t.Errorf("uniform lengths: got %v, want 0.0", got)
	}
}

func TestCalculateReadLengthPattern_VariableLengths(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"A", "ACGTACGTACGT", "ACGT"})
	got := fos.calculateReadLengthPattern(seqs)
	if got <= 0 {
		t.Errorf("variable lengths: got %v, want > 0", got)
	}
}

func TestCalculateReadLengthPattern_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ACGT", "ACGTACGT", "A"})
	a := fos.calculateReadLengthPattern(seqs)
	b := fos.calculateReadLengthPattern(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculateErrorProfile
// ---------------------------------------------------------------------------

func TestCalculateErrorProfile_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateErrorProfile(makeSeqMap(nil))
	if got < 0 || got > 1 {
		t.Errorf("empty: got %v outside [0,1]", got)
	}
}

func TestCalculateErrorProfile_NoShortRuns(t *testing.T) {
	fos := newFOS()
	// Alternating: no run ≥ 3.
	seqs := makeSeqMap([]string{"ACGTACGTACGT"})
	got := fos.calculateErrorProfile(seqs)
	if got != 0.0 {
		t.Errorf("no runs ≥ 3: got %v, want 0.0", got)
	}
}

func TestCalculateErrorProfile_AllHomopolymer(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"AAAAAAAAAA"})
	got := fos.calculateErrorProfile(seqs)
	if got != 1.0 {
		t.Errorf("all homopolymer: got %v, want 1.0", got)
	}
}

func TestCalculateErrorProfile_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"AAAGCGT", "CCCATG"})
	a := fos.calculateErrorProfile(seqs)
	b := fos.calculateErrorProfile(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculatePlatformArtifacts
// ---------------------------------------------------------------------------

func TestCalculatePlatformArtifacts_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculatePlatformArtifacts(makeSeqMap(nil))
	if got < 0 || got > 1 {
		t.Errorf("empty: got %v outside [0,1]", got)
	}
}

func TestCalculatePlatformArtifacts_NormalGC(t *testing.T) {
	fos := newFOS()
	// GC ≈ 50% — within [20%, 80%].
	seqs := makeSeqMap([]string{"ACGTACGTACGT"})
	got := fos.calculatePlatformArtifacts(seqs)
	if got != 0.0 {
		t.Errorf("normal GC: got %v, want 0.0", got)
	}
}

func TestCalculatePlatformArtifacts_AllExtreme(t *testing.T) {
	fos := newFOS()
	// 100% GC > 80% threshold.
	seqs := makeSeqMap([]string{"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"})
	got := fos.calculatePlatformArtifacts(seqs)
	if got != 1.0 {
		t.Errorf("all extreme: got %v, want 1.0", got)
	}
}

func TestCalculatePlatformArtifacts_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ACGT", "GCGC", "TTTT"})
	a := fos.calculatePlatformArtifacts(seqs)
	b := fos.calculatePlatformArtifacts(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculateGCBiasStrength
// ---------------------------------------------------------------------------

func TestCalculateGCBiasStrength_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateGCBiasStrength(makeSeqMap(nil))
	if got < 0 || got > 1 {
		t.Errorf("empty: got %v outside [0,1]", got)
	}
}

func TestCalculateGCBiasStrength_UniformGC(t *testing.T) {
	fos := newFOS()
	// All reads identical GC → std dev = 0.
	seqs := makeSeqMap([]string{"ACGT", "ACGT", "ACGT"})
	got := fos.calculateGCBiasStrength(seqs)
	if got != 0.0 {
		t.Errorf("uniform GC: got %v, want 0.0", got)
	}
}

func TestCalculateGCBiasStrength_HighVariability(t *testing.T) {
	fos := newFOS()
	// One 100% GC and one 0% GC → large std dev.
	seqs := makeSeqMap([]string{"GCGCGCGCGCGC", "ATATATATAT"})
	got := fos.calculateGCBiasStrength(seqs)
	if got <= 0 {
		t.Errorf("high variability: got %v, want > 0", got)
	}
}

func TestCalculateGCBiasStrength_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ACGT", "GCGC", "ATAT"})
	a := fos.calculateGCBiasStrength(seqs)
	b := fos.calculateGCBiasStrength(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculateWindowPatterns
// ---------------------------------------------------------------------------

func TestCalculateWindowPatterns_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateWindowPatterns(makeSeqMap(nil))
	if got < 0 {
		t.Errorf("empty: got %v, want >= 0", got)
	}
}

func TestCalculateWindowPatterns_ShortSeq(t *testing.T) {
	fos := newFOS()
	// Shorter than window size: variance = 0.
	seqs := makeSeqMap([]string{"ACGT"})
	got := fos.calculateWindowPatterns(seqs)
	if got != 0.0 {
		t.Errorf("short seq: got %v, want 0.0", got)
	}
}

func TestCalculateWindowPatterns_UniformGC(t *testing.T) {
	fos := newFOS()
	// All windows have same GC → variance = 0.
	seqs := makeSeqMap([]string{"ACGTACGTACGTACGTACGTACGTACGTACGT"}) // 32 bases, 50% GC throughout
	got := fos.calculateWindowPatterns(seqs)
	if got != 0.0 {
		t.Errorf("uniform GC windows: got %v, want 0.0", got)
	}
}

func TestCalculateWindowPatterns_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ACGTACGTACGTACGTACGTACGTACGTACGT"})
	a := fos.calculateWindowPatterns(seqs)
	b := fos.calculateWindowPatterns(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculatePopulationGCDeviation
// ---------------------------------------------------------------------------

func TestCalculatePopulationGCDeviation_SinglePopulation(t *testing.T) {
	fos := newFOS()
	// Only one population → std dev undefined → 0.
	seqs := makeSeqMap([]string{"ACGT", "ACGT"})
	got := fos.calculatePopulationGCDeviation(seqs)
	if got != 0.0 {
		t.Errorf("single population: got %v, want 0.0", got)
	}
}

func TestCalculatePopulationGCDeviation_TwoIdenticalPops(t *testing.T) {
	fos := newFOS()
	// Two populations with same GC → std dev = 0.
	seqs := map[population.Population][]string{
		population.CEU: {"ACGT", "ACGT"},
		population.YRI: {"ACGT", "ACGT"},
	}
	got := fos.calculatePopulationGCDeviation(seqs)
	if got != 0.0 {
		t.Errorf("identical populations: got %v, want 0.0", got)
	}
}

func TestCalculatePopulationGCDeviation_TwoDifferentPops(t *testing.T) {
	fos := newFOS()
	// One pop all-GC, other all-AT → large std dev.
	seqs := map[population.Population][]string{
		population.CEU: {"GCGCGCGCGCGC"},
		population.YRI: {"ATATATATATAT"},
	}
	got := fos.calculatePopulationGCDeviation(seqs)
	if got <= 0 {
		t.Errorf("different populations: got %v, want > 0", got)
	}
}

func TestCalculatePopulationGCDeviation_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := map[population.Population][]string{
		population.CEU: {"ACGT", "GCGC"},
		population.YRI: {"ATAT", "TTTT"},
	}
	a := fos.calculatePopulationGCDeviation(seqs)
	b := fos.calculatePopulationGCDeviation(seqs)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}

// ---------------------------------------------------------------------------
// calculateEnsembleConsensus
// ---------------------------------------------------------------------------

func TestCalculateEnsembleConsensus_Empty(t *testing.T) {
	fos := newFOS()
	got := fos.calculateEnsembleConsensus(makeSeqMap(nil), nil)
	if got < 0 || got > 1 {
		t.Errorf("empty: got %v outside [0,1]", got)
	}
}

// When GC and quality are positively correlated, result should be > 0.5.
func TestCalculateEnsembleConsensus_PositiveCorrelation(t *testing.T) {
	fos := newFOS()
	// Low-GC seq gets low quality; high-GC seq gets high quality.
	seqs := makeSeqMap([]string{"ATATATATAT", "GCGCGCGCGCGC"}) // GC: 0.0, 1.0
	quals := [][]int{
		{10, 10, 10, 10, 10}, // low Q for low GC
		{40, 40, 40, 40, 40}, // high Q for high GC
	}
	got := fos.calculateEnsembleConsensus(seqs, quals)
	if got <= 0.5 {
		t.Errorf("positive correlation: got %v, want > 0.5", got)
	}
}

// When GC and quality are negatively correlated, result should be < 0.5.
func TestCalculateEnsembleConsensus_NegativeCorrelation(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ATATATATAT", "GCGCGCGCGCGC"}) // GC: 0.0, 1.0
	quals := [][]int{
		{40, 40, 40, 40, 40}, // high Q for low GC
		{10, 10, 10, 10, 10}, // low Q for high GC
	}
	got := fos.calculateEnsembleConsensus(seqs, quals)
	if got >= 0.5 {
		t.Errorf("negative correlation: got %v, want < 0.5", got)
	}
}

func TestCalculateEnsembleConsensus_Deterministic(t *testing.T) {
	fos := newFOS()
	seqs := makeSeqMap([]string{"ACGT", "GCGC"})
	quals := [][]int{{30, 35}, {20, 25}}
	a := fos.calculateEnsembleConsensus(seqs, quals)
	b := fos.calculateEnsembleConsensus(seqs, quals)
	if a != b {
		t.Errorf("not deterministic: %v != %v", a, b)
	}
}
