package testing

import (
	"errors"
	"fmt"
	"math"
	"math/rand/v2"
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

// split returns the even-indexed reads of d and the odd-indexed ones.
func split(d *SyntheticData) (train, test *SyntheticData) {
	train, test = &SyntheticData{}, &SyntheticData{}
	for i, seq := range d.Sequences {
		if i%2 == 0 {
			train.Sequences = append(train.Sequences, seq)
		} else {
			test.Sequences = append(test.Sequences, seq)
		}
	}
	return train, test
}

// HoldoutAccuracy trains a classifier on the even-indexed reads of a and b
// and returns its accuracy on the odd-indexed reads, with the number of
// held-out reads. Each set needs at least two reads.
func HoldoutAccuracy(a, b *SyntheticData) (accuracy float64, holdout int, err error) {
	if a == nil || b == nil || len(a.Sequences) < 2 || len(b.Sequences) < 2 {
		return 0, 0, fmt.Errorf("HoldoutAccuracy: each set needs at least two reads")
	}
	trainA, testA := split(a)
	trainB, testB := split(b)
	c := &MLClassifier{nb: NewMultinomialNB(1.0)}
	if err := c.Train(trainA, trainB); err != nil {
		return 0, 0, err
	}
	acc, err := c.TestAccuracy(testA, testB)
	if err != nil {
		return 0, 0, err
	}
	return acc, len(testA.Sequences) + len(testB.Sequences), nil
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
// skew deterministic.
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
	return &out
}
