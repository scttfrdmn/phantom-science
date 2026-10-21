package testing

import (
	"fmt"
	"math"
	"sync"
)

// LogisticRegression is an L2-regularised two-class logistic regression on
// standardised features, for feature sets with continuous values where
// Multinomial Naive Bayes does not apply.
//
// Fit standardises each feature with the training mean and standard
// deviation (a feature constant in training is dropped), then minimises the
// mean log loss plus (Lambda/2)·|w|² (the bias is not penalised) by
// full-batch gradient descent with Adam steps for Iterations steps from
// w = 0. The loss is convex, so the result does not depend on a random
// start; gradient sums are taken over fixed row chunks in a fixed order, so
// a fit is deterministic for a given input on a given CPU architecture.
type LogisticRegression struct {
	// Lambda is the L2 penalty on the standardised weights.
	Lambda float64
	// Iterations is the number of full-batch Adam steps.
	Iterations int
	// LearningRate is the Adam step size.
	LearningRate float64

	mean, scale []float64 // standardisation; scale 0 drops the feature
	w           []float64
	b           float64
	trained     bool
}

// Default LogisticRegression settings.
const (
	defaultLRLambda       = 1e-3
	defaultLRIterations   = 300
	defaultLRLearningRate = 0.05
	// lrChunk is the number of rows per gradient chunk.
	lrChunk = 2048
)

// NewLogisticRegression returns an untrained classifier with the default
// settings.
func NewLogisticRegression() *LogisticRegression {
	return &LogisticRegression{Lambda: defaultLRLambda, Iterations: defaultLRIterations, LearningRate: defaultLRLearningRate}
}

// Fit trains on features (rows of equal length) and labels (0 or 1, both
// present).
func (lr *LogisticRegression) Fit(features [][]float64, labels []int) error {
	if len(features) == 0 {
		return fmt.Errorf("empty training set")
	}
	if len(features) != len(labels) {
		return fmt.Errorf("features and labels length mismatch: %d vs %d", len(features), len(labels))
	}
	d := len(features[0])
	if d == 0 {
		return fmt.Errorf("feature vectors have zero length")
	}
	var classCount [2]int
	for i, l := range labels {
		if l < 0 || l > 1 {
			return fmt.Errorf("label %d out of range [0,1]", l)
		}
		if len(features[i]) != d {
			return fmt.Errorf("sample %d has %d features, expected %d", i, len(features[i]), d)
		}
		classCount[l]++
	}
	if classCount[0] == 0 || classCount[1] == 0 {
		return fmt.Errorf("both classes need training samples (have %d and %d)", classCount[0], classCount[1])
	}

	n := float64(len(features))
	lr.mean = make([]float64, d)
	lr.scale = make([]float64, d)
	for _, x := range features {
		for k, v := range x {
			lr.mean[k] += v
		}
	}
	for k := range lr.mean {
		lr.mean[k] /= n
	}
	for _, x := range features {
		for k, v := range x {
			dv := v - lr.mean[k]
			lr.scale[k] += dv * dv
		}
	}
	for k := range lr.scale {
		sd := math.Sqrt(lr.scale[k] / n)
		lr.scale[k] = 0
		if sd > 1e-12 {
			lr.scale[k] = 1 / sd
		}
	}
	z := make([][]float64, len(features))
	for i, x := range features {
		z[i] = lr.standardise(x)
	}

	lr.w = make([]float64, d)
	lr.b = 0
	const beta1, beta2, eps = 0.9, 0.999, 1e-8
	m, v := make([]float64, d+1), make([]float64, d+1)
	grad := make([]float64, d+1)
	for t := 1; t <= lr.Iterations; t++ {
		lr.gradient(z, labels, grad)
		for k := 0; k < d; k++ {
			grad[k] += lr.Lambda * lr.w[k]
		}
		c1 := 1 - math.Pow(beta1, float64(t))
		c2 := 1 - math.Pow(beta2, float64(t))
		for k := 0; k <= d; k++ {
			m[k] = beta1*m[k] + (1-beta1)*grad[k]
			v[k] = beta2*v[k] + (1-beta2)*grad[k]*grad[k]
			step := lr.LearningRate * (m[k] / c1) / (math.Sqrt(v[k]/c2) + eps)
			if k < d {
				lr.w[k] -= step
			} else {
				lr.b -= step
			}
		}
	}
	lr.trained = true
	return nil
}

// gradient writes the mean log-loss gradient (weights, then bias) to grad.
func (lr *LogisticRegression) gradient(z [][]float64, labels []int, grad []float64) {
	d := len(lr.w)
	chunks := (len(z) + lrChunk - 1) / lrChunk
	parts := make([][]float64, chunks)
	var wg sync.WaitGroup
	for c := range chunks {
		wg.Add(1)
		go func() {
			defer wg.Done()
			g := make([]float64, d+1)
			for i := c * lrChunk; i < min(len(z), (c+1)*lrChunk); i++ {
				r := sigmoid(lr.score(z[i])) - float64(labels[i])
				for k, x := range z[i] {
					g[k] += r * x
				}
				g[d] += r
			}
			parts[c] = g
		}()
	}
	wg.Wait()
	clear(grad)
	for _, g := range parts { // fixed order: deterministic sums
		for k, x := range g {
			grad[k] += x
		}
	}
	for k := range grad {
		grad[k] /= float64(len(z))
	}
}

func (lr *LogisticRegression) standardise(x []float64) []float64 {
	z := make([]float64, len(x))
	for k, v := range x {
		z[k] = (v - lr.mean[k]) * lr.scale[k]
	}
	return z
}

func (lr *LogisticRegression) score(z []float64) float64 {
	s := lr.b
	for k, x := range z {
		s += lr.w[k] * x
	}
	return s
}

func sigmoid(s float64) float64 {
	if s >= 0 {
		return 1 / (1 + math.Exp(-s))
	}
	e := math.Exp(s)
	return e / (1 + e)
}

// PredictProba returns P(class 1 | x); 0.5 when untrained or when x has
// the wrong length.
func (lr *LogisticRegression) PredictProba(x []float64) float64 {
	if !lr.trained || len(x) != len(lr.w) {
		return 0.5
	}
	return sigmoid(lr.score(lr.standardise(x)))
}

// Predict returns the more probable class (1 on a tie).
func (lr *LogisticRegression) Predict(x []float64) int {
	if lr.PredictProba(x) >= 0.5 {
		return 1
	}
	return 0
}

// stackedFolds is the number of cross-fitting folds StackedClassifier uses
// for the sequence score of its training rows.
const stackedFolds = 5

// StackedClassifier classifies combined feature vectors: the first
// kmerVocabSize elements are 4-mer counts, the rest continuous features.
// The counts are reduced to one score, the log-odds of a Multinomial Naive
// Bayes fitted to them, and a LogisticRegression is fitted to that score
// and the continuous features. Training rows get out-of-fold scores (5-fold
// cross-fitting, rows assigned to folds by index) so the regression learns
// how far to trust the score on reads the Naive Bayes has not seen; at
// prediction the score comes from a Naive Bayes fitted to every training
// row. Fitting the 256 counts in the regression directly instead adds 256
// weakly informative, noisy weights and does worse than either part alone.
type StackedClassifier struct {
	nb *MultinomialNB
	lr *LogisticRegression
}

// NewStackedClassifier returns an untrained StackedClassifier.
func NewStackedClassifier() *StackedClassifier {
	return &StackedClassifier{nb: NewMultinomialNB(1.0), lr: NewLogisticRegression()}
}

// Fit trains on combined vectors (longer than kmerVocabSize) and labels.
func (s *StackedClassifier) Fit(features [][]float64, labels []int) error {
	if len(features) == 0 {
		return fmt.Errorf("empty training set")
	}
	if len(features) != len(labels) {
		return fmt.Errorf("features and labels length mismatch: %d vs %d", len(features), len(labels))
	}
	for i, x := range features {
		if len(x) <= kmerVocabSize {
			return fmt.Errorf("sample %d has %d features; want the %d 4-mer counts followed by other features", i, len(x), kmerVocabSize)
		}
	}
	counts := make([][]float64, len(features))
	for i, x := range features {
		counts[i] = x[:kmerVocabSize]
	}
	score := make([]float64, len(features))
	if len(features) >= 2*stackedFolds {
		for f := range stackedFolds {
			var tx [][]float64
			var ty []int
			for i := range counts {
				if i%stackedFolds != f {
					tx, ty = append(tx, counts[i]), append(ty, labels[i])
				}
			}
			nb := NewMultinomialNB(1.0)
			if err := nb.Fit(tx, ty); err != nil {
				return fmt.Errorf("sequence score fold %d: %w", f, err)
			}
			for i := f; i < len(counts); i += stackedFolds {
				score[i] = nb.LogOdds(counts[i])
			}
		}
	}
	if err := s.nb.Fit(counts, labels); err != nil {
		return err
	}
	if len(features) < 2*stackedFolds { // too few rows to cross-fit
		for i := range counts {
			score[i] = s.nb.LogOdds(counts[i])
		}
	}
	rows := make([][]float64, len(features))
	for i, x := range features {
		rows[i] = append([]float64{score[i]}, x[kmerVocabSize:]...)
	}
	return s.lr.Fit(rows, labels)
}

// Predict returns the predicted class of a combined vector.
func (s *StackedClassifier) Predict(x []float64) int {
	if len(x) <= kmerVocabSize {
		return s.lr.Predict(nil)
	}
	return s.lr.Predict(append([]float64{s.nb.LogOdds(x[:kmerVocabSize])}, x[kmerVocabSize:]...))
}

// Weights returns the fitted weights on the standardised features (0 for a
// feature dropped as constant); nil when untrained.
func (lr *LogisticRegression) Weights() []float64 {
	if !lr.trained {
		return nil
	}
	return append([]float64(nil), lr.w...)
}
