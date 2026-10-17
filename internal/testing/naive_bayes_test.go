//go:build !release

package testing

import (
	"math"
	"math/rand"
	"strings"
	"testing"
)

// generateSeq creates a random DNA sequence of the given length with the
// specified GC probability using the provided RNG.
func generateSeq(rng *rand.Rand, length int, gcProb float64) string {
	const at = "AT"
	const gc = "GC"
	var sb strings.Builder
	sb.Grow(length)
	for i := 0; i < length; i++ {
		if rng.Float64() < gcProb {
			sb.WriteByte(gc[rng.Intn(2)])
		} else {
			sb.WriteByte(at[rng.Intn(2)])
		}
	}
	return sb.String()
}

// ─── Feature extraction ──────────────────────────────────────────────────────

func TestExtract4merFrequencies_CountSum(t *testing.T) {
	seq := "ACGTACGTACGT" // length 12 → 9 4-mers
	counts := Extract4merFrequencies(seq)
	if len(counts) != kmerVocabSize {
		t.Fatalf("expected %d features, got %d", kmerVocabSize, len(counts))
	}
	total := 0.0
	for _, c := range counts {
		total += c
	}
	expected := float64(len(seq) - kmerSize + 1) // 12 - 4 + 1 = 9
	if math.Abs(total-expected) > 1e-9 {
		t.Errorf("counts sum mismatch: expected %.0f, got %.0f", expected, total)
	}
}

func TestExtract4merFrequencies_ShortSequence(t *testing.T) {
	for _, seq := range []string{"", "A", "AC", "ACG"} {
		counts := Extract4merFrequencies(seq)
		for i, c := range counts {
			if c != 0 {
				t.Errorf("short seq %q: expected all-zero, got counts[%d]=%.0f", seq, i, c)
			}
		}
	}
}

func TestExtract4merFrequencies_NormalisedSumsToOne(t *testing.T) {
	seq := "ACGTACGTACGTACGT" // 16 bases → 13 4-mers
	counts := Extract4merFrequencies(seq)
	total := 0.0
	for _, c := range counts {
		total += c
	}
	// Normalise manually.
	norm := make([]float64, len(counts))
	for i, c := range counts {
		norm[i] = c / total
	}
	normSum := 0.0
	for _, v := range norm {
		normSum += v
	}
	if math.Abs(normSum-1.0) > 1e-9 {
		t.Errorf("normalised vector does not sum to 1: %.12f", normSum)
	}
}

// ─── Classifier correctness ──────────────────────────────────────────────────

func TestMultinomialNB_Separable(t *testing.T) {
	// Class 0: GC-rich (70% GC), Class 1: AT-rich (70% AT).
	// These distributions are very different in 4-mer space → CV accuracy > 90%.
	rng := rand.New(rand.NewSource(12345)) // #nosec G404
	const n = 100
	const seqLen = 200

	features := make([][]float64, 2*n)
	labels := make([]int, 2*n)
	for i := 0; i < n; i++ {
		features[i] = Extract4merFrequencies(generateSeq(rng, seqLen, 0.70)) // GC-rich
		labels[i] = 0
	}
	for i := 0; i < n; i++ {
		features[n+i] = Extract4merFrequencies(generateSeq(rng, seqLen, 0.30)) // AT-rich
		labels[n+i] = 1
	}

	nb := NewMultinomialNB(1.0)
	mean, _, err := nb.CrossValidate(features, labels, 5)
	if err != nil {
		t.Fatalf("CrossValidate failed: %v", err)
	}
	if mean < 0.90 {
		t.Errorf("separable test: expected CV accuracy > 0.90, got %.4f", mean)
	}
}

func TestMultinomialNB_Indistinguishable(t *testing.T) {
	// Both classes drawn from the same distribution → classifier at ~chance.
	rng := rand.New(rand.NewSource(99999)) // #nosec G404
	const n = 100
	const seqLen = 200

	features := make([][]float64, 2*n)
	labels := make([]int, 2*n)
	for i := 0; i < n; i++ {
		features[i] = Extract4merFrequencies(generateSeq(rng, seqLen, 0.50))
		labels[i] = 0
	}
	for i := 0; i < n; i++ {
		features[n+i] = Extract4merFrequencies(generateSeq(rng, seqLen, 0.50))
		labels[n+i] = 1
	}

	nb := NewMultinomialNB(1.0)
	mean, _, err := nb.CrossValidate(features, labels, 5)
	if err != nil {
		t.Fatalf("CrossValidate failed: %v", err)
	}
	if mean < 0.40 || mean > 0.60 {
		t.Errorf("indistinguishable test: expected CV accuracy in [0.40, 0.60], got %.4f", mean)
	}
}

func TestMultinomialNB_CVDeterministic(t *testing.T) {
	// With the same seed embedded in CrossValidate (42), results should be
	// deterministic across two calls on the same data.
	rng := rand.New(rand.NewSource(777)) // #nosec G404
	const n = 60

	features := make([][]float64, 2*n)
	labels := make([]int, 2*n)
	for i := 0; i < n; i++ {
		features[i] = Extract4merFrequencies(generateSeq(rng, 100, 0.65))
		labels[i] = 0
	}
	for i := 0; i < n; i++ {
		features[n+i] = Extract4merFrequencies(generateSeq(rng, 100, 0.35))
		labels[n+i] = 1
	}

	nb1 := NewMultinomialNB(1.0)
	m1, s1, err := nb1.CrossValidate(features, labels, 5)
	if err != nil {
		t.Fatalf("first CV failed: %v", err)
	}
	nb2 := NewMultinomialNB(1.0)
	m2, s2, err := nb2.CrossValidate(features, labels, 5)
	if err != nil {
		t.Fatalf("second CV failed: %v", err)
	}
	if m1 != m2 || s1 != s2 {
		t.Errorf("CV is not deterministic: mean (%.4f vs %.4f), std (%.4f vs %.4f)", m1, m2, s1, s2)
	}
}

func TestMultinomialNB_CVStdReasonable(t *testing.T) {
	rng := rand.New(rand.NewSource(42)) // #nosec G404
	const n = 80

	features := make([][]float64, 2*n)
	labels := make([]int, 2*n)
	for i := 0; i < n; i++ {
		features[i] = Extract4merFrequencies(generateSeq(rng, 150, 0.68))
		labels[i] = 0
	}
	for i := 0; i < n; i++ {
		features[n+i] = Extract4merFrequencies(generateSeq(rng, 150, 0.32))
		labels[n+i] = 1
	}

	nb := NewMultinomialNB(1.0)
	_, std, err := nb.CrossValidate(features, labels, 5)
	if err != nil {
		t.Fatalf("CrossValidate failed: %v", err)
	}
	if std > 0.15 {
		t.Errorf("CV std is unreasonably large: %.4f (expected < 0.15)", std)
	}
}

func TestMultinomialNB_FitPredictBasic(t *testing.T) {
	// Train on perfectly separable data; all training samples should be predicted correctly.
	features := [][]float64{
		Extract4merFrequencies(strings.Repeat("GCGCGC", 20)), // GC class 0
		Extract4merFrequencies(strings.Repeat("GCGCGC", 20)),
		Extract4merFrequencies(strings.Repeat("ATATAT", 20)), // AT class 1
		Extract4merFrequencies(strings.Repeat("ATATAT", 20)),
	}
	labels := []int{0, 0, 1, 1}

	nb := NewMultinomialNB(1.0)
	if err := nb.Fit(features, labels); err != nil {
		t.Fatalf("Fit failed: %v", err)
	}

	for i, feat := range features {
		pred := nb.Predict(feat)
		if pred != labels[i] {
			t.Errorf("sample %d: expected label %d, got %d", i, labels[i], pred)
		}
	}
}

func TestMultinomialNB_PredictProbaNormalised(t *testing.T) {
	features := [][]float64{
		Extract4merFrequencies(strings.Repeat("GCGC", 25)),
		Extract4merFrequencies(strings.Repeat("ATAT", 25)),
	}
	nb := NewMultinomialNB(1.0)
	if err := nb.Fit(features, []int{0, 1}); err != nil {
		t.Fatalf("Fit failed: %v", err)
	}
	proba := nb.PredictProba(features[0])
	sum := proba[0] + proba[1]
	if math.Abs(sum-1.0) > 1e-9 {
		t.Errorf("PredictProba does not sum to 1: %.12f", sum)
	}
}

func TestMultinomialNB_EmptyTrainingSet(t *testing.T) {
	nb := NewMultinomialNB(1.0)
	if err := nb.Fit([][]float64{}, []int{}); err == nil {
		t.Error("expected error for empty training set")
	}
}

func TestMultinomialNB_MismatchedLabels(t *testing.T) {
	nb := NewMultinomialNB(1.0)
	features := [][]float64{Extract4merFrequencies("ACGT")}
	if err := nb.Fit(features, []int{0, 1}); err == nil {
		t.Error("expected error for features/labels length mismatch")
	}
}
