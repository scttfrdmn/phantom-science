package testing

import (
	"fmt"
	"math"
)

// FeatureSet names the per-read features a read-set comparison classifies
// on (#29).
type FeatureSet string

const (
	// FeatureSequence is the 256 4-mer counts of the read's bases,
	// classified with Multinomial Naive Bayes. It sees base composition and
	// local sequence context, not qualities.
	FeatureSequence FeatureSet = "sequence"
	// FeatureQuality is the quality profile of the read (QualityFeatureNames),
	// classified with logistic regression. It sees the quality model and the
	// no-calls, and the base composition only through the mean quality of
	// each called base.
	FeatureQuality FeatureSet = "quality"
	// FeatureCombined is the 4-mer counts followed by the quality features,
	// classified with logistic regression.
	FeatureCombined FeatureSet = "combined"
)

// AllFeatureSets lists the feature sets in report order.
var AllFeatureSets = []FeatureSet{FeatureSequence, FeatureQuality, FeatureCombined}

// Quality feature layout. Qualities are Phred+33; values at or above
// qualHistBins-1 share the last histogram bin.
const (
	qualHistBins = 51 // Q0..Q49 and Q>=50
	// qualBlocks is the number of equal read segments (cycle blocks) whose
	// mean quality is a feature.
	qualBlocks = 10
	// lowQualThreshold: a base with Phred below it is low quality. On
	// NovaSeq 6000 binned output this is Q2 and Q12.
	lowQualThreshold = 20
)

// qualityFeatureNames is the name of each element of QualityFeatures.
var qualityFeatureNames = func() []string {
	var n []string
	for q := 0; q < qualHistBins-1; q++ {
		n = append(n, fmt.Sprintf("qual_frac_Q%d", q))
	}
	n = append(n, fmt.Sprintf("qual_frac_Q>=%d", qualHistBins-1))
	for b := 1; b <= qualBlocks; b++ {
		n = append(n, fmt.Sprintf("mean_qual_block_%d/%d", b, qualBlocks))
	}
	n = append(n,
		"mean_qual", "sd_qual",
		fmt.Sprintf("lowq_runs_per_100bp(Q<%d)", lowQualThreshold),
		fmt.Sprintf("lowq_longest_run(Q<%d)", lowQualThreshold),
		fmt.Sprintf("lowq_tail_run(Q<%d)", lowQualThreshold),
		"qual_change_frac",
		"n_frac",
		"mean_qual_A", "mean_qual_C", "mean_qual_G", "mean_qual_T",
	)
	return n
}()

// QualityFeatureNames returns the names of the QualityFeatures elements.
func QualityFeatureNames() []string { return append([]string(nil), qualityFeatureNames...) }

// QualityFeatures returns the quality profile of one read, with qual its
// Phred+33 qualities (len(qual) == len(seq)):
//
//   - the fraction of bases at each Phred value Q0..Q49 and Q>=50;
//   - the mean quality of each of 10 equal cycle blocks (block b holds
//     cycles [b·L/10, (b+1)·L/10); an empty block, in reads shorter than
//     10, takes the read's mean);
//   - the read's mean and standard deviation of quality;
//   - low-quality runs (maximal runs of Phred < 20): runs per 100 bases,
//     longest run, and the length of the run ending at the last cycle;
//   - the fraction of adjacent cycles whose qualities differ;
//   - the fraction of N bases;
//   - the mean quality of the bases called A, C, G and T (the read's mean
//     when the base does not occur).
//
// A read without bases gives a zero vector.
func QualityFeatures(seq, qual []byte) []float64 {
	f := make([]float64, len(qualityFeatureNames))
	n := len(qual)
	if n == 0 {
		return f
	}
	hist := f[:qualHistBins]
	blocks := f[qualHistBins : qualHistBins+qualBlocks]
	rest := f[qualHistBins+qualBlocks:]

	var sum, sumSq float64
	var blockSum [qualBlocks]float64
	var blockN [qualBlocks]int
	var baseSum [4]float64
	var baseN [4]int
	var runs, longest, run, changes, ns int
	for i, c := range qual {
		q := int(c) - 33
		if q < 0 {
			q = 0
		}
		hist[min(q, qualHistBins-1)]++
		sum += float64(q)
		sumSq += float64(q) * float64(q)
		b := i * qualBlocks / n
		blockSum[b] += float64(q)
		blockN[b]++
		if q < lowQualThreshold {
			if run == 0 {
				runs++
			}
			run++
			longest = max(longest, run)
		} else {
			run = 0
		}
		if i > 0 && qual[i] != qual[i-1] {
			changes++
		}
		if i < len(seq) {
			switch seq[i] {
			case 'A', 'a':
				baseSum[0] += float64(q)
				baseN[0]++
			case 'C', 'c':
				baseSum[1] += float64(q)
				baseN[1]++
			case 'G', 'g':
				baseSum[2] += float64(q)
				baseN[2]++
			case 'T', 't':
				baseSum[3] += float64(q)
				baseN[3]++
			case 'N', 'n':
				ns++
			}
		}
	}
	fn := float64(n)
	for i := range hist {
		hist[i] /= fn
	}
	mean := sum / fn
	for b := range blocks {
		blocks[b] = mean
		if blockN[b] > 0 {
			blocks[b] = blockSum[b] / float64(blockN[b])
		}
	}
	rest[0] = mean
	rest[1] = math.Sqrt(math.Max(0, sumSq/fn-mean*mean))
	rest[2] = 100 * float64(runs) / fn
	rest[3] = float64(longest)
	rest[4] = float64(run) // the run still open at the last cycle
	if n > 1 {
		rest[5] = float64(changes) / float64(n-1)
	}
	rest[6] = float64(ns) / fn
	for k := range 4 {
		rest[7+k] = mean
		if baseN[k] > 0 {
			rest[7+k] = baseSum[k] / float64(baseN[k])
		}
	}
	return f
}

// kmerNames returns the 4-mer of each Extract4merFrequencies index.
var kmerNames = func() []string {
	n := make([]string, kmerVocabSize)
	for i := range n {
		b := make([]byte, kmerSize)
		for j := kmerSize - 1; j >= 0; j-- {
			b[j] = "ACGT"[(i>>(2*(kmerSize-1-j)))&3]
		}
		n[i] = string(b)
	}
	return n
}()

// FeatureNames returns the name of each element of set's feature vector.
func FeatureNames(set FeatureSet) ([]string, error) {
	switch set {
	case FeatureSequence:
		return append([]string(nil), kmerNames...), nil
	case FeatureQuality:
		return QualityFeatureNames(), nil
	case FeatureCombined:
		return append(append([]string(nil), kmerNames...), qualityFeatureNames...), nil
	}
	return nil, fmt.Errorf("unknown feature set %q", set)
}

// ReadFeatures returns set's feature vector for one read. The quality sets
// need qualities (len(qual) == len(seq)).
func ReadFeatures(set FeatureSet, seq, qual []byte) ([]float64, error) {
	if set != FeatureSequence && len(qual) != len(seq) {
		return nil, fmt.Errorf("feature set %s needs qualities: read has %d bases and %d qualities", set, len(seq), len(qual))
	}
	switch set {
	case FeatureSequence:
		return Extract4merFrequencies(string(seq)), nil
	case FeatureQuality:
		return QualityFeatures(seq, qual), nil
	case FeatureCombined:
		return append(Extract4merFrequencies(string(seq)), QualityFeatures(seq, qual)...), nil
	}
	return nil, fmt.Errorf("unknown feature set %q", set)
}

// BinaryClassifier is a two-class classifier on feature vectors.
type BinaryClassifier interface {
	Fit(features [][]float64, labels []int) error
	Predict(feature []float64) int
}

// newClassifier returns the classifier used for set: Multinomial Naive
// Bayes for the 4-mer counts, which are multinomial counts, logistic
// regression for the continuous quality features, and for the combined set
// a logistic regression on the quality features and the Naive Bayes
// sequence score (StackedClassifier).
func newClassifier(set FeatureSet) BinaryClassifier {
	switch set {
	case FeatureSequence:
		return NewMultinomialNB(1.0)
	case FeatureCombined:
		return NewStackedClassifier()
	}
	return NewLogisticRegression()
}

// ClassifierName describes the classifier HoldoutAccuracyFor uses for set.
func ClassifierName(set FeatureSet) string {
	switch set {
	case FeatureSequence:
		return "multinomial naive Bayes"
	case FeatureCombined:
		return "stacked NB score + LR"
	}
	return "logistic regression"
}
