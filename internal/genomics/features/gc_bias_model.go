package features

import "math"

// GCBiasModel reports the GC content of a sequence.
type GCBiasModel struct {
	windowSize int
}

// GCBiasDetails holds the GC fraction. WindowScores is never filled.
type GCBiasDetails struct {
	GCContent    float64
	WindowScores []float64
}

// NewGCBiasModel creates a GC content model. windowSize (default 100) is
// stored but not used.
func NewGCBiasModel(windowSize int) (*GCBiasModel, error) {
	if windowSize <= 0 {
		windowSize = 100
	}
	return &GCBiasModel{windowSize: windowSize}, nil
}

// AnalyzeSequence returns the whole-sequence GC fraction (G and C, either
// case, over all bytes) as the score. There is no windowing and no
// comparison with an expected bias.
func (m *GCBiasModel) AnalyzeSequence(seq []byte) (float64, *GCBiasDetails, error) {
	if len(seq) == 0 {
		return 0.0, &GCBiasDetails{}, nil
	}
	gc := 0
	for _, b := range seq {
		if b == 'G' || b == 'C' || b == 'g' || b == 'c' {
			gc++
		}
	}
	score := float64(gc) / float64(len(seq))
	return score, &GCBiasDetails{GCContent: score}, nil
}

// RepetitiveElementModel scores sequence complexity from dinucleotide
// entropy. It does not detect repeat families (no RepeatMasker or Dfam data).
type RepetitiveElementModel struct {
	complexityThreshold float64
}

// RepetitiveElementDetails contains the complexity score and homopolymer run
// count.
type RepetitiveElementDetails struct {
	RepeatScore float64
	RepeatCount int
}

// NewRepetitiveElementModel creates a repetitive element model.
// complexityThreshold is stored but not used.
func NewRepetitiveElementModel(complexityThreshold float64) (*RepetitiveElementModel, error) {
	return &RepetitiveElementModel{complexityThreshold: complexityThreshold}, nil
}

// AnalyzeSequence analyzes repetitive elements in a DNA sequence.
//
// Algorithm: entropy-based complexity score using dinucleotide (2-mer) frequencies.
//
//	entropy = −Σ p_i · log₂(p_i)   over all observed dinucleotides
//	score   = entropy / log₂(min(16, n−1))   (normalised to [0, 1])
//
// score ≈ 1.0 → high complexity (non-repetitive)
// score ≈ 0.0 → low complexity (highly repetitive, e.g. AAAAAAA…)
//
// RepeatCount is the number of homopolymer runs of length ≥ 3.
func (m *RepetitiveElementModel) AnalyzeSequence(seq []byte) (float64, *RepetitiveElementDetails, error) {
	if len(seq) == 0 {
		return 0.0, &RepetitiveElementDetails{}, nil
	}

	// Count dinucleotide (2-mer) frequencies.
	type dinuc [2]byte
	kmerCounts := make(map[dinuc]int, 16)
	for i := 0; i < len(seq)-1; i++ {
		kmerCounts[dinuc{seq[i], seq[i+1]}]++
	}

	total := float64(len(seq) - 1)
	entropy := 0.0
	for _, count := range kmerCounts {
		p := float64(count) / total
		if p > 0 {
			entropy -= p * math.Log2(p)
		}
	}

	// Normalise by the maximum possible entropy.
	// maxDistinct is bounded by 16 (all AA,AT,…,TT) and by the number of available
	// dinucleotide positions (len(seq)-1).
	maxDistinct := math.Min(16, total)
	var repeatScore float64
	if maxDistinct > 1 {
		repeatScore = math.Min(1.0, math.Max(0.0, entropy/math.Log2(maxDistinct)))
	} else {
		// Single-base sequences are treated as fully repetitive.
		repeatScore = 0.0
	}

	// Count homopolymer runs of length ≥ 3 (each run counted once).
	repeatCount := 0
	runLen := 1
	for i := 1; i < len(seq); i++ {
		if seq[i] == seq[i-1] {
			runLen++
			if runLen == 3 {
				repeatCount++
			}
		} else {
			runLen = 1
		}
	}

	return repeatScore, &RepetitiveElementDetails{
		RepeatScore: repeatScore,
		RepeatCount: repeatCount,
	}, nil
}
