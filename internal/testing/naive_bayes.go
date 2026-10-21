package testing

import (
	"fmt"
	"math"
	"math/rand"
)

const (
	// kmerSize is the k for 4-mer feature extraction (256 possible 4-mers).
	kmerSize = 4
	// kmerVocabSize is 4^kmerSize = 256.
	kmerVocabSize = 256
)

// kmerIndex maps a DNA 4-mer string to an integer index in [0, 255].
// Bases are encoded A=0, C=1, G=2, T=3.
// Unknown bases are mapped to index 0 (treated as AAAA).
func kmerIndex(kmer []byte) int {
	idx := 0
	for _, b := range kmer {
		idx <<= 2
		switch b {
		case 'A', 'a':
			idx |= 0
		case 'C', 'c':
			idx |= 1
		case 'G', 'g':
			idx |= 2
		case 'T', 't':
			idx |= 3
		}
	}
	return idx
}

// Extract4merFrequencies returns a 256-element count vector, one per canonical
// 4-mer, computed by sliding a window of width 4 over seq.
//
// Sequences shorter than 4 bases return a zero vector.
// The counts sum to len(seq)−3 for sequences of length ≥ 4.
func Extract4merFrequencies(seq string) []float64 {
	counts := make([]float64, kmerVocabSize)
	if len(seq) < kmerSize {
		return counts
	}
	data := []byte(seq)
	for i := 0; i <= len(data)-kmerSize; i++ {
		idx := kmerIndex(data[i : i+kmerSize])
		counts[idx]++
	}
	return counts
}

// MultinomialNB is a Multinomial Naive Bayes classifier for two classes
// (label 0 = real, label 1 = synthetic) trained on k-mer count vectors.
//
// Training uses Laplace smoothing with parameter alpha:
//
//	θ_{c,k} = (count_k(c) + α) / (Σ_j count_j(c) + α·V)
//
// where V is the vocabulary size (256 for 4-mers).
//
// Prediction uses log-probability:
//
//	log P(c|x) = log P(c) + Σ_k x_k · log θ_{c,k}
type MultinomialNB struct {
	alpha     float64      // Laplace smoothing parameter
	logTheta  [2][]float64 // log smoothed feature probabilities per class
	logPrior  [2]float64   // log class priors
	nFeatures int          // expected feature vector length
	trained   bool
}

// NewMultinomialNB creates a new MultinomialNB with Laplace smoothing parameter alpha.
// A typical value is alpha=1.0 (add-one smoothing).
func NewMultinomialNB(alpha float64) *MultinomialNB {
	return &MultinomialNB{alpha: alpha}
}

// Fit trains the classifier on the given feature matrix and labels.
// labels must contain only 0 and 1. features[i] must have the same length for all i.
func (nb *MultinomialNB) Fit(features [][]float64, labels []int) error {
	if len(features) == 0 {
		return fmt.Errorf("empty training set")
	}
	if len(features) != len(labels) {
		return fmt.Errorf("features and labels length mismatch: %d vs %d", len(features), len(labels))
	}

	nFeatures := len(features[0])
	if nFeatures == 0 {
		return fmt.Errorf("feature vectors have zero length")
	}

	// Count samples per class.
	classCount := [2]int{}
	for _, l := range labels {
		if l < 0 || l > 1 {
			return fmt.Errorf("label %d out of range [0,1]", l)
		}
		classCount[l]++
	}
	for c := 0; c < 2; c++ {
		if classCount[c] == 0 {
			return fmt.Errorf("class %d has no training samples", c)
		}
	}

	// Accumulate feature counts per class.
	featureSum := [2][]float64{
		make([]float64, nFeatures),
		make([]float64, nFeatures),
	}
	for i, feats := range features {
		if len(feats) != nFeatures {
			return fmt.Errorf("sample %d has %d features, expected %d", i, len(feats), nFeatures)
		}
		c := labels[i]
		for k, v := range feats {
			featureSum[c][k] += v
		}
	}

	// Compute smoothed log-probabilities.
	V := float64(nFeatures)
	for c := 0; c < 2; c++ {
		nb.logTheta[c] = make([]float64, nFeatures)
		totalCount := 0.0
		for _, v := range featureSum[c] {
			totalCount += v
		}
		denom := totalCount + nb.alpha*V
		for k := 0; k < nFeatures; k++ {
			nb.logTheta[c][k] = math.Log((featureSum[c][k] + nb.alpha) / denom)
		}
	}

	// Log class priors.
	total := float64(len(labels))
	for c := 0; c < 2; c++ {
		nb.logPrior[c] = math.Log(float64(classCount[c]) / total)
	}

	nb.nFeatures = nFeatures
	nb.trained = true
	return nil
}

// Predict returns the predicted class label (0 or 1) for a feature vector.
func (nb *MultinomialNB) Predict(feature []float64) int {
	proba := nb.PredictProba(feature)
	if proba[1] >= proba[0] {
		return 1
	}
	return 0
}

// PredictProba returns the posterior probability [P(class=0|x), P(class=1|x)].
// Uses log-sum-exp for numerical stability.
func (nb *MultinomialNB) PredictProba(feature []float64) [2]float64 {
	if !nb.trained {
		return [2]float64{0.5, 0.5}
	}

	logProbs := [2]float64{}
	for c := 0; c < 2; c++ {
		lp := nb.logPrior[c]
		for k, v := range feature {
			if v > 0 {
				lp += v * nb.logTheta[c][k]
			}
		}
		logProbs[c] = lp
	}

	// Softmax via log-sum-exp.
	maxLP := math.Max(logProbs[0], logProbs[1])
	sum := math.Exp(logProbs[0]-maxLP) + math.Exp(logProbs[1]-maxLP)
	return [2]float64{
		math.Exp(logProbs[0]-maxLP) / sum,
		math.Exp(logProbs[1]-maxLP) / sum,
	}
}

// LogOdds returns log P(class=1|x) − log P(class=0|x); 0 when untrained.
func (nb *MultinomialNB) LogOdds(feature []float64) float64 {
	if !nb.trained {
		return 0
	}
	lo := nb.logPrior[1] - nb.logPrior[0]
	for k, v := range feature {
		if v > 0 {
			lo += v * (nb.logTheta[1][k] - nb.logTheta[0][k])
		}
	}
	return lo
}

// CrossValidate performs k-fold stratified cross-validation and returns the
// mean and standard deviation of per-fold accuracy.
//
// The data is shuffled with the provided rng (use rand.New(rand.NewSource(seed))
// for reproducibility). folds must be ≥ 2.
func (nb *MultinomialNB) CrossValidate(features [][]float64, labels []int, folds int) (mean, std float64, err error) {
	if folds < 2 {
		return 0, 0, fmt.Errorf("folds must be ≥ 2, got %d", folds)
	}
	if len(features) < folds {
		return 0, 0, fmt.Errorf("fewer samples (%d) than folds (%d)", len(features), folds)
	}

	// Separate indices by class for stratified splitting.
	class0, class1 := []int{}, []int{}
	for i, l := range labels {
		if l == 0 {
			class0 = append(class0, i)
		} else {
			class1 = append(class1, i)
		}
	}
	if len(class0) == 0 || len(class1) == 0 {
		return 0, 0, fmt.Errorf("both classes must have at least one sample")
	}

	// Use a fixed local RNG so results are reproducible.
	rng := rand.New(rand.NewSource(42)) // #nosec G404
	shuffleInts(rng, class0)
	shuffleInts(rng, class1)

	// Build folds by interleaving class0 and class1 chunks.
	foldIndices := make([][]int, folds)
	assignToFolds(class0, folds, foldIndices)
	assignToFolds(class1, folds, foldIndices)

	accuracies := make([]float64, folds)
	for f := 0; f < folds; f++ {
		testIdx := foldIndices[f]
		trainIdx := []int{}
		for g := 0; g < folds; g++ {
			if g != f {
				trainIdx = append(trainIdx, foldIndices[g]...)
			}
		}

		trainFeats := make([][]float64, len(trainIdx))
		trainLabels := make([]int, len(trainIdx))
		for i, idx := range trainIdx {
			trainFeats[i] = features[idx]
			trainLabels[i] = labels[idx]
		}

		fold := NewMultinomialNB(nb.alpha)
		if e := fold.Fit(trainFeats, trainLabels); e != nil {
			return 0, 0, fmt.Errorf("fold %d training failed: %w", f, e)
		}

		correct := 0
		for _, idx := range testIdx {
			if fold.Predict(features[idx]) == labels[idx] {
				correct++
			}
		}
		accuracies[f] = float64(correct) / float64(len(testIdx))
	}

	// Mean and std of accuracies.
	sum := 0.0
	for _, a := range accuracies {
		sum += a
	}
	mean = sum / float64(folds)

	sumSq := 0.0
	for _, a := range accuracies {
		d := a - mean
		sumSq += d * d
	}
	// Bessel's correction: sample std dev uses N-1. folds >= 2 is enforced above.
	std = math.Sqrt(sumSq / float64(folds-1))
	return mean, std, nil
}

// shuffleInts shuffles a slice of ints in-place using the provided RNG.
func shuffleInts(rng *rand.Rand, s []int) {
	for i := len(s) - 1; i > 0; i-- {
		j := rng.Intn(i + 1)
		s[i], s[j] = s[j], s[i]
	}
}

// assignToFolds distributes indices round-robin across folds.
func assignToFolds(indices []int, folds int, foldIndices [][]int) {
	for i, idx := range indices {
		f := i % folds
		foldIndices[f] = append(foldIndices[f], idx)
	}
}
