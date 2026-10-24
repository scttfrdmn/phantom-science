package testing

import (
	"errors"
	"fmt"
	"math"
	"math/rand/v2"
	"sort"
)

// MLClassifier separates two read sets with a Multinomial Naive Bayes
// classifier on 4-mer count vectors (one vector per read). Class 0 is the
// first set passed to Train and TestAccuracy, class 1 the second.
type MLClassifier struct {
	nb      *MultinomialNB
	trained bool
}

// errEmptySet is returned when a read set has no reads.
var errEmptySet = errors.New("read set has no sequences")

// CreateMLClassifier returns an untrained classifier.
func (f *E2EFramework) CreateMLClassifier() *MLClassifier {
	return &MLClassifier{nb: NewMultinomialNB(1.0)}
}

// Train fits the classifier to a (class 0) and b (class 1).
func (c *MLClassifier) Train(a, b *SyntheticData) error {
	features, labels, err := labelledFeatures(a, b)
	if err != nil {
		return fmt.Errorf("MLClassifier.Train: %w", err)
	}
	if err := c.nb.Fit(features, labels); err != nil {
		return fmt.Errorf("MLClassifier.Train: %w", err)
	}
	c.trained = true
	return nil
}

// TestAccuracy returns the fraction of reads in a and b that the classifier
// assigns to their own class. A trained classifier is scored on a and b as
// held-out data; an untrained one is scored by 5-fold cross-validation on
// them. Accuracy near 0.5 means this classifier cannot tell the sets apart;
// it does not rule out other classifiers or features.
func (c *MLClassifier) TestAccuracy(a, b *SyntheticData) (float64, error) {
	features, labels, err := labelledFeatures(a, b)
	if err != nil {
		return 0, fmt.Errorf("MLClassifier.TestAccuracy: %w", err)
	}
	if !c.trained {
		mean, _, err := NewMultinomialNB(1.0).CrossValidate(features, labels, 5)
		if err != nil {
			return 0, fmt.Errorf("MLClassifier.TestAccuracy: %w", err)
		}
		return mean, nil
	}
	correct := 0
	for i, f := range features {
		if c.nb.Predict(f) == labels[i] {
			correct++
		}
	}
	return float64(correct) / float64(len(features)), nil
}

// labelledFeatures returns 4-mer vectors for a (label 0) and b (label 1).
func labelledFeatures(a, b *SyntheticData) ([][]float64, []int, error) {
	if a == nil || b == nil || len(a.Sequences) == 0 || len(b.Sequences) == 0 {
		return nil, nil, errEmptySet
	}
	features := make([][]float64, 0, len(a.Sequences)+len(b.Sequences))
	labels := make([]int, 0, cap(features))
	for label, d := range []*SyntheticData{a, b} {
		for _, seq := range d.Sequences {
			features = append(features, Extract4merFrequencies(string(seq)))
			labels = append(labels, label)
		}
	}
	return features, labels, nil
}

// split returns the even-indexed reads of d and the odd-indexed ones, with
// their qualities and headers, and for paired reads their mates.
func split(d *SyntheticData) (even, odd *SyntheticData) {
	even = &SyntheticData{Platform: d.Platform, Population: d.Population, Seed: d.Seed}
	odd = &SyntheticData{Platform: d.Platform, Population: d.Population, Seed: d.Seed}
	if d.Mates != nil {
		even.Mates, odd.Mates = split(d.Mates)
	}
	for i, seq := range d.Sequences {
		dst := even
		if i%2 == 1 {
			dst = odd
		}
		dst.Sequences = append(dst.Sequences, seq)
		if i < len(d.QualityScores) {
			dst.QualityScores = append(dst.QualityScores, d.QualityScores[i])
		}
		if i < len(d.Headers) {
			dst.Headers = append(dst.Headers, d.Headers[i])
		}
	}
	return even, odd
}

// SplitAlternate returns the even-indexed reads of d and the odd-indexed
// ones (with qualities and headers), each with its Stats computed. Two
// halves of one read set are the null comparison for that set: a classifier
// separating them measures nothing but chance. Alternating reads, rather
// than the first and second half of the file, keeps both halves drawn from
// the same tiles of a real run.
func SplitAlternate(d *SyntheticData) (even, odd *SyntheticData) {
	even, odd = split(d)
	for _, h := range []*SyntheticData{even, odd} {
		for x := h; x != nil; x = x.Mates {
			x.Stats = ComputeReadStats(x.Sequences, x.QualityScores)
		}
	}
	return even, odd
}

// HoldoutAccuracy trains a classifier on the even-indexed reads of a and b
// and returns its accuracy on the odd-indexed reads, with the number of
// held-out reads. Each set needs at least two reads. It uses the sequence
// feature set (4-mer counts, Multinomial Naive Bayes); HoldoutAccuracyFor
// covers the quality-aware sets.
func HoldoutAccuracy(a, b *SyntheticData) (accuracy float64, holdout int, err error) {
	r, err := HoldoutAccuracyFor(FeatureSequence, a, b, 0)
	if err != nil {
		return 0, 0, err
	}
	return r.Accuracy, r.Holdout, nil
}

// FeatureDiff is how one feature differs between two read sets: its mean
// per read in each set and the standardised mean difference
// (MeanB − MeanA) / pooled SD (Cohen's d).
type FeatureDiff struct {
	Name       string  `json:"name"`
	MeanA      float64 `json:"mean_a"`
	MeanB      float64 `json:"mean_b"`
	EffectSize float64 `json:"effect_size"`
}

// Holdout is the result of HoldoutAccuracyFor.
type Holdout struct {
	FeatureSet FeatureSet `json:"feature_set"`
	Classifier string     `json:"classifier"`
	// Accuracy is the fraction of held-out reads assigned to their own set.
	Accuracy float64 `json:"accuracy"`
	// Holdout is the number of held-out reads (both sets).
	Holdout int `json:"holdout_reads"`
	// Top lists the features with the largest |EffectSize| on the training
	// reads, largest first. 4-mer features are per-read frequencies (count
	// over the read's 4-mers) so that reads of different lengths compare.
	Top []FeatureDiff `json:"top_features,omitempty"`
}

// HoldoutAccuracyFor trains set's classifier (ClassifierName) on the
// even-indexed reads of a (class 0) and b (class 1) and scores it on the
// odd-indexed reads. It also ranks the features by effect size on the
// training reads and keeps the top ones (none when top is 0). Each set
// needs at least two reads, qualities for the quality-aware sets, and
// mates (SyntheticData.Mates) for FeaturePair.
func HoldoutAccuracyFor(set FeatureSet, a, b *SyntheticData, top int) (*Holdout, error) {
	if a == nil || b == nil || len(a.Sequences) < 2 || len(b.Sequences) < 2 {
		return nil, fmt.Errorf("HoldoutAccuracy: each set needs at least two reads")
	}
	trainA, testA := split(a)
	trainB, testB := split(b)
	xTrain, yTrain, err := labelledSetFeatures(set, trainA, trainB)
	if err != nil {
		return nil, fmt.Errorf("HoldoutAccuracy: %w", err)
	}
	xTest, yTest, err := labelledSetFeatures(set, testA, testB)
	if err != nil {
		return nil, fmt.Errorf("HoldoutAccuracy: %w", err)
	}
	c := newClassifier(set)
	if err := c.Fit(xTrain, yTrain); err != nil {
		return nil, fmt.Errorf("HoldoutAccuracy: %w", err)
	}
	correct := 0
	for i, x := range xTest {
		if c.Predict(x) == yTest[i] {
			correct++
		}
	}
	r := &Holdout{FeatureSet: set, Classifier: ClassifierName(set),
		Accuracy: float64(correct) / float64(len(xTest)), Holdout: len(xTest)}
	if top > 0 {
		names, err := FeatureNames(set)
		if err != nil {
			return nil, err
		}
		r.Top = topFeatureDiffs(names, xTrain, yTrain, set == FeatureSequence || set == FeatureCombined, top)
	}
	return r, nil
}

// labelledSetFeatures returns set's feature vectors for a (label 0) and b
// (label 1).
func labelledSetFeatures(set FeatureSet, a, b *SyntheticData) ([][]float64, []int, error) {
	if a == nil || b == nil || len(a.Sequences) == 0 || len(b.Sequences) == 0 {
		return nil, nil, errEmptySet
	}
	features := make([][]float64, 0, len(a.Sequences)+len(b.Sequences))
	labels := make([]int, 0, cap(features))
	for label, d := range []*SyntheticData{a, b} {
		if set == FeaturePair && !d.Paired() {
			return nil, nil, fmt.Errorf("feature set %s: set %d is not paired (one mate per read)", set, label)
		}
		for i, seq := range d.Sequences {
			var x []float64
			var err error
			if set == FeaturePair {
				m := d.Mates
				x, err = PairFeatures(seq, qualAt(d, i), m.Sequences[i], qualAt(m, i))
			} else {
				x, err = ReadFeatures(set, seq, qualAt(d, i))
			}
			if err != nil {
				return nil, nil, fmt.Errorf("read %d of set %d: %w", i, label, err)
			}
			features = append(features, x)
			labels = append(labels, label)
		}
	}
	return features, labels, nil
}

// qualAt returns the qualities of read i of d, or nil.
func qualAt(d *SyntheticData, i int) []byte {
	if i < len(d.QualityScores) {
		return d.QualityScores[i]
	}
	return nil
}

// topFeatureDiffs ranks features by |Cohen's d| between class 0 and 1 and
// returns the top ones. With kmerFreq the first kmerVocabSize features are
// 4-mer counts and are converted to per-read frequencies first.
func topFeatureDiffs(names []string, x [][]float64, y []int, kmerFreq bool, top int) []FeatureDiff {
	d := len(names)
	var sum, sumSq [2][]float64
	var n [2]float64
	for c := range 2 {
		sum[c], sumSq[c] = make([]float64, d), make([]float64, d)
	}
	for i, row := range x {
		c := y[i]
		n[c]++
		total := 0.0
		if kmerFreq {
			for k := 0; k < kmerVocabSize && k < len(row); k++ {
				total += row[k]
			}
		}
		for k, v := range row {
			if kmerFreq && k < kmerVocabSize {
				v = 0
				if total > 0 {
					v = row[k] / total
				}
			}
			sum[c][k] += v
			sumSq[c][k] += v * v
		}
	}
	diffs := make([]FeatureDiff, 0, d)
	for k := range d {
		ma, mb := sum[0][k]/n[0], sum[1][k]/n[1]
		va := math.Max(0, sumSq[0][k]/n[0]-ma*ma)
		vb := math.Max(0, sumSq[1][k]/n[1]-mb*mb)
		pooled := math.Sqrt((n[0]*va + n[1]*vb) / (n[0] + n[1]))
		es := 0.0
		if pooled > 1e-12 {
			es = (mb - ma) / pooled
		}
		diffs = append(diffs, FeatureDiff{Name: names[k], MeanA: ma, MeanB: mb, EffectSize: es})
	}
	sort.SliceStable(diffs, func(i, j int) bool { return math.Abs(diffs[i].EffectSize) > math.Abs(diffs[j].EffectSize) })
	return diffs[:min(top, len(diffs))]
}

// QualityShifted returns a copy of data with every quality lowered by
// delta Phred (floored at Q2), bases unchanged. It is a positive control
// for the quality-aware feature sets on reads whose quality model cannot be
// switched, such as real reads. Mates, when present, are shifted too.
func QualityShifted(data *SyntheticData, delta int) *SyntheticData {
	out := *data
	out.QualityScores = make([][]byte, len(data.QualityScores))
	for i, q := range data.QualityScores {
		s := make([]byte, len(q))
		for j, c := range q {
			s[j] = byte(max(int(c)-delta, 33+2))
		}
		out.QualityScores[i] = s
	}
	out.Stats = ComputeReadStats(out.Sequences, out.QualityScores)
	if data.Mates != nil {
		out.Mates = QualityShifted(data.Mates, delta)
	}
	return &out
}

// ChanceBand returns the half-width of the sigmas-sigma band around 0.5 for
// accuracy measured on n held-out reads: under the null hypothesis that the
// sets cannot be told apart, correct calls are Binomial(n, 0.5).
func ChanceBand(n int, sigmas float64) float64 {
	if n <= 0 {
		return math.Inf(1)
	}
	return sigmas * math.Sqrt(0.25/float64(n))
}

// GCSkewed returns a copy of data whose expected GC content is targetGC: each
// A or T base is replaced by G or C with probability
// q = (targetGC - g) / (1 - g), where g is data's measured GC content. When
// targetGC <= g the sequences are copied unchanged. The input is not
// modified; headers and qualities are shared with it, and Stats is
// recomputed.
//
// It is the positive control for classifier comparisons: a classifier that
// cannot separate reads from their GC-skewed counterparts cannot detect any
// real distributional difference of that size either. The seed makes the
// skew deterministic. Mates, when present, are skewed too (seed+1).
func GCSkewed(data *SyntheticData, targetGC float64, seed uint64) *SyntheticData {
	g := data.Stats.GCContent
	if data.Stats.Reads == 0 {
		g = ComputeReadStats(data.Sequences, nil).GCContent
	}
	q := 0.0
	if targetGC > g && g < 1 {
		q = math.Min(1, (targetGC-g)/(1-g))
	}
	rng := rand.New(rand.NewPCG(seed, 0x6763736b6577)) // #nosec G404 -- deterministic test fixture
	out := *data
	out.Sequences = make([][]byte, len(data.Sequences))
	for i, seq := range data.Sequences {
		s := make([]byte, len(seq))
		for j, b := range seq {
			if (b == 'A' || b == 'T') && rng.Float64() < q {
				b = "GC"[rng.IntN(2)]
			}
			s[j] = b
		}
		out.Sequences[i] = s
	}
	out.Stats = ComputeReadStats(out.Sequences, out.QualityScores)
	if data.Mates != nil {
		out.Mates = GCSkewed(data.Mates, targetGC, seed+1)
	}
	return &out
}
