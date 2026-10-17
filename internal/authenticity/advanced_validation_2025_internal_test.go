//go:build !release

package authenticity

import (
	"math"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
)

// ---------------------------------------------------------------------------
// AncestryClassifier.PredictProbability
// ---------------------------------------------------------------------------

func TestPredictProbability_Untrained(t *testing.T) {
	a := &AncestryClassifier{}
	f := TabularFeature{Population: population.CEU}
	got := a.PredictProbability(f)
	if got != 0.5 {
		t.Errorf("untrained: got %v, want 0.5", got)
	}
}

func TestPredictProbability_OneCentroid(t *testing.T) {
	a := &AncestryClassifier{}
	// Train on a single population — no second centroid, so confidence must be 1.0.
	a.Train([]TabularFeature{
		{Population: population.CEU, GCContent: 0.42, Length: 100, Complexity: 0.8, Quality: 35},
	})
	f := TabularFeature{Population: population.CEU, GCContent: 0.42, Length: 100, Complexity: 0.8, Quality: 35}
	got := a.PredictProbability(f)
	if got != 1.0 {
		t.Errorf("one centroid: got %v, want 1.0", got)
	}
}

func TestPredictProbability_TwoCentroids_HighConfidence(t *testing.T) {
	a := &AncestryClassifier{}
	// Two very distinct populations.
	a.Train([]TabularFeature{
		{Population: population.CEU, GCContent: 0.1, Quality: 10},
		{Population: population.YRI, GCContent: 0.9, Quality: 40},
	})
	// Query right at the CEU centroid → distance to CEU ≈ 0, distance to YRI >> 0.
	f := TabularFeature{Population: population.CEU, GCContent: 0.1, Quality: 10}
	got := a.PredictProbability(f)
	if got <= 0.5 {
		t.Errorf("high-confidence case: got %v, want > 0.5", got)
	}
}

func TestPredictProbability_Range(t *testing.T) {
	a := &AncestryClassifier{}
	a.Train([]TabularFeature{
		{Population: population.CEU, GCContent: 0.3},
		{Population: population.YRI, GCContent: 0.7},
	})
	f := TabularFeature{Population: population.CEU, GCContent: 0.5}
	got := a.PredictProbability(f)
	if got < 0 || got > 1 {
		t.Errorf("probability %v outside [0,1]", got)
	}
}

// ---------------------------------------------------------------------------
// GANDiscriminator.Predict
// ---------------------------------------------------------------------------

func TestGANDiscriminatorPredict_EmptySeq(t *testing.T) {
	g := &GANDiscriminator{}
	seq := GenomicSequence{Sequence: ""}
	got := g.Predict(seq)
	if got != 0.5 {
		t.Errorf("empty seq: got %v, want 0.5", got)
	}
}

func TestGANDiscriminatorPredict_PerfectGC(t *testing.T) {
	g := &GANDiscriminator{}
	// Sequence with ~42% GC and high quality → score > 0.5.
	seq := GenomicSequence{
		Sequence: "ACGTACGTACGTACGTACGTACGT", // GC = 50% (close to 42%)
		QualityScores: []int{
			38, 38, 38, 38, 38, 38, 38, 38,
			38, 38, 38, 38, 38, 38, 38, 38,
			38, 38, 38, 38, 38, 38, 38, 38,
		},
	}
	got := g.Predict(seq)
	if got <= 0.5 {
		t.Errorf("balanced GC + high quality: got %v, want > 0.5", got)
	}
}

func TestGANDiscriminatorPredict_ExtremeGC(t *testing.T) {
	g := &GANDiscriminator{}
	// All-G sequence: GC = 100%, deviation = 0.58 > 0.20 → gcScore = 0.
	// Low quality → qualScore near 0 → total near 0.
	seq := GenomicSequence{
		Sequence:      "GGGGGGGGGGGGGGGGGGGG",
		QualityScores: []int{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
	}
	got := g.Predict(seq)
	if got >= 0.3 {
		t.Errorf("extreme GC + low quality: got %v, want < 0.3", got)
	}
}

func TestGANDiscriminatorPredict_InRange(t *testing.T) {
	g := &GANDiscriminator{}
	seqs := []GenomicSequence{
		{Sequence: "ACGTACGT", QualityScores: []int{30, 30, 30, 30, 30, 30, 30, 30}},
		{Sequence: "GGGGGGGG"},
		{Sequence: "AAAAAAAA", QualityScores: []int{10}},
	}
	for _, seq := range seqs {
		got := g.Predict(seq)
		if got < 0 || got > 1 {
			t.Errorf("score %v outside [0,1] for seq %q", got, seq.Sequence)
		}
	}
}

// ---------------------------------------------------------------------------
// TabularTransformer.CalculateAttention
// ---------------------------------------------------------------------------

func TestCalculateAttention_Empty(t *testing.T) {
	tr := &TabularTransformer{}
	got := tr.CalculateAttention(nil)
	if got != nil {
		t.Errorf("empty input: got %v, want nil", got)
	}
}

func TestCalculateAttention_SingleFeature(t *testing.T) {
	tr := &TabularTransformer{}
	features := []TabularFeature{{GCContent: 0.42, Length: 100, Complexity: 0.8, Quality: 35}}
	got := tr.CalculateAttention(features)
	if len(got) != 1 || len(got[0]) != 1 {
		t.Fatalf("single feature: got %dx? matrix, want 1x1", len(got))
	}
	if math.Abs(got[0][0]-1.0) > 1e-9 {
		t.Errorf("single-feature attention weight = %v, want 1.0", got[0][0])
	}
}

func TestCalculateAttention_RowsSumToOne(t *testing.T) {
	tr := &TabularTransformer{}
	features := []TabularFeature{
		{GCContent: 0.40, Length: 100, Quality: 30},
		{GCContent: 0.50, Length: 150, Quality: 35},
		{GCContent: 0.45, Length: 120, Quality: 32},
	}
	got := tr.CalculateAttention(features)
	if len(got) != 3 {
		t.Fatalf("expected 3 rows, got %d", len(got))
	}
	for i, row := range got {
		if len(row) != 3 {
			t.Errorf("row %d: expected 3 cols, got %d", i, len(row))
			continue
		}
		sum := 0.0
		for _, v := range row {
			sum += v
		}
		if math.Abs(sum-1.0) > 1e-9 {
			t.Errorf("row %d sum = %v, want 1.0", i, sum)
		}
	}
}

func TestCalculateAttention_ValuesInRange(t *testing.T) {
	tr := &TabularTransformer{}
	features := []TabularFeature{
		{GCContent: 0.42}, {GCContent: 0.55}, {GCContent: 0.30},
	}
	for _, row := range tr.CalculateAttention(features) {
		for _, v := range row {
			if v < 0 || v > 1 {
				t.Errorf("attention weight %v outside [0,1]", v)
			}
		}
	}
}

// ---------------------------------------------------------------------------
// SelectionInferenceEngine.DetectSelection
// ---------------------------------------------------------------------------

func TestDetectSelection_Empty(t *testing.T) {
	s := &SelectionInferenceEngine{}
	got := s.DetectSelection(map[genomics.SequencingPlatform][]GenomicSequence{})
	if len(got) != 0 {
		t.Errorf("empty data: got %v entries, want 0", len(got))
	}
}

func TestDetectSelection_IdenticalSequences(t *testing.T) {
	s := &SelectionInferenceEngine{}
	// Identical sequences → GC variance = 0, nucleotide diversity = 0.
	data := map[genomics.SequencingPlatform][]GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "ACGTACGT"},
			{Sequence: "ACGTACGT"},
			{Sequence: "ACGTACGT"},
		},
	}
	got := s.DetectSelection(data)
	entry, ok := got["NovaSeq 6000"]
	if !ok {
		t.Fatalf("NovaSeq 6000 key missing from result")
	}
	stats, ok := entry.(map[string]float64)
	if !ok {
		t.Fatalf("unexpected value type for entry: %T", entry)
	}
	if stats["gc_variance"] != 0.0 {
		t.Errorf("identical seqs: gc_variance = %v, want 0.0", stats["gc_variance"])
	}
	if stats["nucleotide_diversity"] != 0.0 {
		t.Errorf("identical seqs: nucleotide_diversity = %v, want 0.0", stats["nucleotide_diversity"])
	}
}

func TestDetectSelection_DiverseSequences(t *testing.T) {
	s := &SelectionInferenceEngine{}
	// Diverse sequences → nucleotide_diversity > 0.
	data := map[genomics.SequencingPlatform][]GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "AAAAAAAA"},
			{Sequence: "CCCCCCCC"},
			{Sequence: "GGGGGGGG"},
		},
	}
	got := s.DetectSelection(data)
	entry, ok := got["NovaSeq 6000"]
	if !ok {
		t.Fatalf("NovaSeq 6000 key missing")
	}
	stats := entry.(map[string]float64)
	if stats["nucleotide_diversity"] <= 0 {
		t.Errorf("diverse seqs: nucleotide_diversity = %v, want > 0", stats["nucleotide_diversity"])
	}
}

// ---------------------------------------------------------------------------
// DemographicInferenceEngine.InferDemographics
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// SelectionInferenceEngine.DetectSelection — Bessel's correction (#86)
// ---------------------------------------------------------------------------

func TestDetectSelection_SampleVariance_TwoSeqs(t *testing.T) {
	s := &SelectionInferenceEngine{}
	// Two sequences with GC 0.0 and 1.0: sample variance = 0.5, pop variance = 0.25.
	data := map[genomics.SequencingPlatform][]GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "AAAA"}, // GC = 0.0
			{Sequence: "GGGG"}, // GC = 1.0
		},
	}
	got := s.DetectSelection(data)
	entry, ok := got["NovaSeq 6000"]
	if !ok {
		t.Fatalf("expected NovaSeq 6000 key")
	}
	stats := entry.(map[string]float64)
	// sample variance = ((0-0.5)² + (1-0.5)²) / (2-1) = 0.5
	const want = 0.5
	if math.Abs(stats["gc_variance"]-want) > 1e-9 {
		t.Errorf("gc_variance = %v, want %v (Bessel-corrected)", stats["gc_variance"], want)
	}
}

func TestDetectSelection_SingleSeq_ZeroVariance(t *testing.T) {
	s := &SelectionInferenceEngine{}
	// Single sequence: variance undefined → must be 0.0, not NaN.
	data := map[genomics.SequencingPlatform][]GenomicSequence{
		genomics.NovaSeq6000: {{Sequence: "ACGTACGT"}},
	}
	got := s.DetectSelection(data)
	stats := got["NovaSeq 6000"].(map[string]float64)
	v := stats["gc_variance"]
	if v != 0.0 || v != v { // second condition catches NaN
		t.Errorf("single seq: gc_variance = %v, want 0.0", v)
	}
}

// ---------------------------------------------------------------------------
// GANDiscriminator gcMean + FinetuneDicriminator (#89)
// ---------------------------------------------------------------------------

func TestGANDiscriminatorPredict_UsesGcMean(t *testing.T) {
	// A discriminator with gcMean=1.0 scores a 100%-GC sequence more highly
	// than the default discriminator (gcMean=0.42) does.
	defaultDisc := &GANDiscriminator{}
	tunedDisc := &GANDiscriminator{gcMean: 1.0, gcMeanSet: true}
	seq := GenomicSequence{
		Sequence:      "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC", // GC = 100%
		QualityScores: []int{35, 35, 35, 35, 35},
	}
	defaultScore := defaultDisc.Predict(seq)
	tunedScore := tunedDisc.Predict(seq)
	if tunedScore <= defaultScore {
		t.Errorf("tuned(gcMean=1.0) score %v should exceed default(gcMean=0.42) score %v for 100%%-GC seq",
			tunedScore, defaultScore)
	}
}

func TestFinetuneDicriminator_UpdatesGcMean(t *testing.T) {
	ft := &DiscriminatorFinetuning{}
	base := &GANDiscriminator{} // gcMean = 0 (uses default 0.42)

	// All sequences are 100% GC → finetuned gcMean should be ≈ 1.0.
	data := map[genomics.SequencingPlatform][]GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "GCGCGCGCGCGCGCGC"},
			{Sequence: "GGGGGGGGGGGGGGGG"},
		},
	}
	finetuned := ft.FinetuneDicriminator(base, data)
	if math.Abs(finetuned.gcMean-1.0) > 1e-9 {
		t.Errorf("finetuned gcMean = %v, want 1.0 for all-GC data", finetuned.gcMean)
	}
}

func TestFinetuneDicriminator_NoData_ReturnsSame(t *testing.T) {
	ft := &DiscriminatorFinetuning{}
	base := &GANDiscriminator{gcMean: 0.42}
	finetuned := ft.FinetuneDicriminator(base, nil)
	if finetuned != base {
		t.Errorf("no data: expected original discriminator pointer returned unchanged")
	}
}

func TestFinetuneDicriminator_DoesNotMutateOriginal(t *testing.T) {
	ft := &DiscriminatorFinetuning{}
	base := &GANDiscriminator{gcMean: 0.42}
	data := map[genomics.SequencingPlatform][]GenomicSequence{
		genomics.NovaSeq6000: {{Sequence: "GCGCGCGCGCGC"}},
	}
	_ = ft.FinetuneDicriminator(base, data)
	if base.gcMean != 0.42 {
		t.Errorf("original discriminator gcMean mutated to %v; expected 0.42", base.gcMean)
	}
}

func TestFinetuneDicriminator_ChangesAccuracy(t *testing.T) {
	ft := &DiscriminatorFinetuning{}
	disc := &GANDiscriminator{} // default gcMean=0.42

	// Training data: all AT sequences (GC = 0%)
	allAT := map[genomics.SequencingPlatform][]GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "ATATATATATAT"},
			{Sequence: "ATATATATAT"},
		},
	}
	finetuned := ft.FinetuneDicriminator(disc, allAT)

	// finetuned gcMean should be ≈ 0 (all AT sequences)
	if finetuned.gcMean >= 0.1 {
		t.Errorf("finetuned gcMean = %v, want near 0 for all-AT data", finetuned.gcMean)
	}

	// An AT-rich sequence should score higher with the finetuned discriminator.
	atSeq := GenomicSequence{Sequence: "ATATATATATAT", QualityScores: []int{35, 35, 35}}
	defaultScore := disc.Predict(atSeq)
	tunedScore := finetuned.Predict(atSeq)
	if tunedScore <= defaultScore {
		t.Errorf("AT-rich seq: finetuned score %v should exceed default score %v",
			tunedScore, defaultScore)
	}
}

// ---------------------------------------------------------------------------
// DemographicInferenceEngine.InferDemographics
// ---------------------------------------------------------------------------

func TestInferDemographics_Empty(t *testing.T) {
	d := &DemographicInferenceEngine{}
	got := d.InferDemographics(map[genomics.SequencingPlatform][]GenomicSequence{}, nil)
	if len(got) != 0 {
		t.Errorf("empty data: got %d entries, want 0", len(got))
	}
}

func TestInferDemographics_BasicStats(t *testing.T) {
	d := &DemographicInferenceEngine{}
	data := map[genomics.SequencingPlatform][]GenomicSequence{
		genomics.NovaSeq6000: {
			{Sequence: "ACGTACGT"},   // len=8, GC=50%
			{Sequence: "ACGTACGTAC"}, // len=10, GC=50%
		},
	}
	pops := []population.Population{population.CEU, population.YRI}
	got := d.InferDemographics(data, pops)

	entry, ok := got["NovaSeq 6000"]
	if !ok {
		t.Fatalf("NovaSeq 6000 key missing")
	}
	stats, ok := entry.(map[string]interface{})
	if !ok {
		t.Fatalf("unexpected value type %T", entry)
	}

	if stats["sequence_count"] != 2 {
		t.Errorf("sequence_count = %v, want 2", stats["sequence_count"])
	}
	wantMeanLen := float64(8+10) / 2.0
	if stats["mean_length"] != wantMeanLen {
		t.Errorf("mean_length = %v, want %v", stats["mean_length"], wantMeanLen)
	}
	if stats["population_count"] != 2 {
		t.Errorf("population_count = %v, want 2", stats["population_count"])
	}
	meanGC, ok := stats["mean_gc"].(float64)
	if !ok || meanGC < 0 || meanGC > 1 {
		t.Errorf("mean_gc = %v, want float in [0,1]", stats["mean_gc"])
	}
}
