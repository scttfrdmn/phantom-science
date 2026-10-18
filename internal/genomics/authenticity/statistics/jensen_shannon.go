// Package statistics implements distribution and population-genetics
// statistics: Jensen-Shannon divergence and related distances between k-mer
// distributions, and (in non-release builds only) Weir-Cockerham F_ST, a
// Hardy-Weinberg chi-square test and pairwise linkage disequilibrium.
package statistics

import (
	"fmt"
	"math"
	"sort"
)

// JensenShannonDivergence calculates the Jensen-Shannon divergence between two probability distributions
type JensenShannonDivergence struct {
	Base      float64 // Base for logarithm (2 for bits, math.E for nats)
	Smoothing float64 // Smoothing parameter for zero probabilities
}

// NewJensenShannonDivergence creates a new Jensen-Shannon divergence calculator
func NewJensenShannonDivergence() *JensenShannonDivergence {
	return &JensenShannonDivergence{
		Base:      2.0,   // Use base 2 by default (bits)
		Smoothing: 1e-10, // Small smoothing to avoid log(0)
	}
}

// Calculate computes the Jensen-Shannon divergence between two probability distributions
func (js *JensenShannonDivergence) Calculate(p, q map[string]float64) (float64, error) {
	if len(p) == 0 || len(q) == 0 {
		return 0.0, fmt.Errorf("probability distributions cannot be empty")
	}

	// Validate that probabilities sum to 1
	if err := js.validateProbabilityDistribution(p); err != nil {
		return 0.0, fmt.Errorf("invalid probability distribution p: %w", err)
	}
	if err := js.validateProbabilityDistribution(q); err != nil {
		return 0.0, fmt.Errorf("invalid probability distribution q: %w", err)
	}

	// Get union of all keys
	allKeys := js.getUnionKeys(p, q)

	// Create normalized distributions with smoothing
	pNorm := js.normalizeWithSmoothing(p, allKeys)
	qNorm := js.normalizeWithSmoothing(q, allKeys)

	// Calculate the average distribution M = (P + Q) / 2
	m := make(map[string]float64)
	for key := range allKeys {
		m[key] = (pNorm[key] + qNorm[key]) / 2.0
	}

	// Calculate KL divergences
	klPM := js.calculateKLDivergence(pNorm, m)
	klQM := js.calculateKLDivergence(qNorm, m)

	// Jensen-Shannon divergence is the average of the two KL divergences
	jsDivergence := (klPM + klQM) / 2.0

	return jsDivergence, nil
}

// CalculateDistance returns the Jensen-Shannon distance (square root of divergence)
func (js *JensenShannonDivergence) CalculateDistance(p, q map[string]float64) (float64, error) {
	divergence, err := js.Calculate(p, q)
	if err != nil {
		return 0.0, err
	}
	return math.Sqrt(divergence), nil
}

// CalculateFromFrequencies calculates JS divergence from frequency maps
func (js *JensenShannonDivergence) CalculateFromFrequencies(freqP, freqQ map[string]int) (float64, error) {
	// Convert frequency maps to probability distributions
	probP := js.frequenciesToProbabilities(freqP)
	probQ := js.frequenciesToProbabilities(freqQ)

	return js.Calculate(probP, probQ)
}

// frequenciesToProbabilities converts frequency map to probability distribution
func (js *JensenShannonDivergence) frequenciesToProbabilities(frequencies map[string]int) map[string]float64 {
	total := 0
	for _, count := range frequencies {
		total += count
	}

	if total == 0 {
		return make(map[string]float64)
	}

	probabilities := make(map[string]float64)
	for key, count := range frequencies {
		probabilities[key] = float64(count) / float64(total)
	}

	return probabilities
}

// validateProbabilityDistribution checks if a map represents a valid probability distribution
func (js *JensenShannonDivergence) validateProbabilityDistribution(dist map[string]float64) error {
	sum := 0.0
	for key, prob := range dist {
		if prob < 0 {
			return fmt.Errorf("negative probability for key %s: %f", key, prob)
		}
		sum += prob
	}

	// Allow small numerical errors
	if math.Abs(sum-1.0) > 1e-6 {
		return fmt.Errorf("probabilities do not sum to 1: sum = %f", sum)
	}

	return nil
}

// getUnionKeys returns all unique keys from both distributions
func (js *JensenShannonDivergence) getUnionKeys(p, q map[string]float64) map[string]bool {
	allKeys := make(map[string]bool)

	for key := range p {
		allKeys[key] = true
	}
	for key := range q {
		allKeys[key] = true
	}

	return allKeys
}

// normalizeWithSmoothing creates a normalized distribution with smoothing for missing keys
func (js *JensenShannonDivergence) normalizeWithSmoothing(dist map[string]float64, allKeys map[string]bool) map[string]float64 {
	normalized := make(map[string]float64)

	// Add smoothing to all keys
	for key := range allKeys {
		if prob, exists := dist[key]; exists {
			normalized[key] = prob + js.Smoothing
		} else {
			normalized[key] = js.Smoothing
		}
	}

	// Renormalize to sum to 1
	sum := 0.0
	for _, prob := range normalized {
		sum += prob
	}

	if sum > 0 {
		for key := range normalized {
			normalized[key] /= sum
		}
	}

	return normalized
}

// calculateKLDivergence computes Kullback-Leibler divergence KL(P || Q)
func (js *JensenShannonDivergence) calculateKLDivergence(p, q map[string]float64) float64 {
	klDiv := 0.0

	for key, pProb := range p {
		qProb := q[key]

		if pProb > 0 && qProb > 0 {
			klDiv += pProb * math.Log(pProb/qProb) / math.Log(js.Base)
		}
	}

	return klDiv
}

// CompareDistributions provides a comprehensive comparison of two probability distributions
func (js *JensenShannonDivergence) CompareDistributions(p, q map[string]float64) (*DistributionComparison, error) {
	// Calculate JS divergence and distance
	jsDivergence, err := js.Calculate(p, q)
	if err != nil {
		return nil, err
	}

	jsDistance := math.Sqrt(jsDivergence)

	// Calculate additional metrics
	comparison := &DistributionComparison{
		JSDivergence:          jsDivergence,
		JSDistance:            jsDistance,
		TotalVariation:        js.calculateTotalVariation(p, q),
		HellingerDistance:     js.calculateHellingerDistance(p, q),
		BhattacharyyaDistance: js.calculateBhattacharyyaDistance(p, q),
		WassersteinDistance:   js.calculateWassersteinDistance(p, q),
		ChiSquared:            js.calculateChiSquared(p, q),
		KeysOnlyInP:           js.getKeysOnlyIn(p, q),
		KeysOnlyInQ:           js.getKeysOnlyIn(q, p),
		CommonKeys:            js.getCommonKeys(p, q),
	}

	return comparison, nil
}

// DistributionComparison contains comprehensive comparison metrics
type DistributionComparison struct {
	JSDivergence          float64
	JSDistance            float64
	TotalVariation        float64
	HellingerDistance     float64
	BhattacharyyaDistance float64
	WassersteinDistance   float64
	ChiSquared            float64
	KeysOnlyInP           []string
	KeysOnlyInQ           []string
	CommonKeys            []string
}

// calculateTotalVariation computes the total variation distance
func (js *JensenShannonDivergence) calculateTotalVariation(p, q map[string]float64) float64 {
	allKeys := js.getUnionKeys(p, q)
	totalVar := 0.0

	for key := range allKeys {
		pProb := p[key]
		qProb := q[key]
		totalVar += math.Abs(pProb - qProb)
	}

	return totalVar / 2.0
}

// calculateHellingerDistance computes the Hellinger distance
func (js *JensenShannonDivergence) calculateHellingerDistance(p, q map[string]float64) float64 {
	allKeys := js.getUnionKeys(p, q)
	sum := 0.0

	for key := range allKeys {
		pProb := p[key]
		qProb := q[key]
		sum += math.Pow(math.Sqrt(pProb)-math.Sqrt(qProb), 2)
	}

	return math.Sqrt(sum / 2.0)
}

// calculateBhattacharyyaDistance computes the Bhattacharyya distance between
// two discrete probability distributions represented as k-mer frequency maps.
//
// B(P,Q) = -ln( Σ_k √(p_k · q_k) )
//
// For completely disjoint distributions (no shared support), the Bhattacharyya
// coefficient is 0 and the distance is +Inf. Callers should check for Inf.
func (js *JensenShannonDivergence) calculateBhattacharyyaDistance(p, q map[string]float64) float64 {
	allKeys := js.getUnionKeys(p, q)
	sum := 0.0

	for key := range allKeys {
		pProb := p[key]
		qProb := q[key]
		if pProb > 0 && qProb > 0 {
			sum += math.Sqrt(pProb * qProb)
		}
	}

	if sum <= 0 {
		// Completely disjoint distributions: Bhattacharyya distance is +Inf.
		return math.Inf(1)
	}
	return -math.Log(sum)
}

// calculateWassersteinDistance computes the 1-Wasserstein distance on the
// lexicographically-ordered k-mer support of two discrete probability distributions.
// It integrates |CDF_P − CDF_Q| over the N sorted support points and normalises by
// N−1 (the maximum possible sum), yielding a result in [0, 1].
//
// Note: this uses uniform transport cost (lexicographic rank). It is NOT the
// earth-mover distance with Hamming cost between k-mer strings. For biological
// comparisons of k-mer distributions, Jensen-Shannon divergence is preferred.
func (js *JensenShannonDivergence) calculateWassersteinDistance(p, q map[string]float64) float64 {
	allKeys := js.getUnionKeys(p, q)
	if len(allKeys) == 0 {
		return 0.0
	}

	keys := make([]string, 0, len(allKeys))
	for k := range allKeys {
		keys = append(keys, k)
	}
	sort.Strings(keys)

	if len(keys) <= 1 {
		return 0.0
	}

	cdfP, cdfQ := 0.0, 0.0
	distance := 0.0
	for _, k := range keys {
		cdfP += p[k]
		cdfQ += q[k]
		distance += math.Abs(cdfP - cdfQ)
	}

	// Normalise by N−1: the maximum sum of |CDF differences| over N ordered points.
	return distance / float64(len(keys)-1)
}

// calculateChiSquared computes the Chi-squared statistic
func (js *JensenShannonDivergence) calculateChiSquared(p, q map[string]float64) float64 {
	allKeys := js.getUnionKeys(p, q)
	chiSq := 0.0

	for key := range allKeys {
		pProb := p[key]
		qProb := q[key]
		if qProb > 0 {
			chiSq += math.Pow(pProb-qProb, 2) / qProb
		}
	}

	return chiSq
}

// getKeysOnlyIn returns keys that exist only in the first distribution
func (js *JensenShannonDivergence) getKeysOnlyIn(first, second map[string]float64) []string {
	var keys []string
	for key := range first {
		if _, exists := second[key]; !exists {
			keys = append(keys, key)
		}
	}
	return keys
}

// getCommonKeys returns keys that exist in both distributions
func (js *JensenShannonDivergence) getCommonKeys(p, q map[string]float64) []string {
	var keys []string
	for key := range p {
		if _, exists := q[key]; exists {
			keys = append(keys, key)
		}
	}
	return keys
}

// BatchComparison performs Jensen-Shannon divergence calculations for multiple distribution pairs
type BatchComparison struct {
	calculator *JensenShannonDivergence
	results    []BatchResult
}

// BatchResult represents the result of a single comparison in a batch
type BatchResult struct {
	ID           string
	JSDivergence float64
	JSDistance   float64
	Error        error
}

// NewBatchComparison creates a new batch comparison instance
func NewBatchComparison() *BatchComparison {
	return &BatchComparison{
		calculator: NewJensenShannonDivergence(),
		results:    make([]BatchResult, 0),
	}
}

// AddComparison adds a comparison to the batch
func (bc *BatchComparison) AddComparison(id string, p, q map[string]float64) {
	jsDivergence, err := bc.calculator.Calculate(p, q)
	var jsDistance float64

	if err == nil {
		jsDistance = math.Sqrt(jsDivergence)
	}

	bc.results = append(bc.results, BatchResult{
		ID:           id,
		JSDivergence: jsDivergence,
		JSDistance:   jsDistance,
		Error:        err,
	})
}

// GetResults returns all batch comparison results
func (bc *BatchComparison) GetResults() []BatchResult {
	return bc.results
}

// GetSuccessfulResults returns only successful comparisons
func (bc *BatchComparison) GetSuccessfulResults() []BatchResult {
	var successful []BatchResult
	for _, result := range bc.results {
		if result.Error == nil {
			successful = append(successful, result)
		}
	}
	return successful
}

// GetAverageJSDivergence returns the average JS divergence of successful comparisons
func (bc *BatchComparison) GetAverageJSDivergence() float64 {
	successful := bc.GetSuccessfulResults()
	if len(successful) == 0 {
		return 0.0
	}

	sum := 0.0
	for _, result := range successful {
		sum += result.JSDivergence
	}

	return sum / float64(len(successful))
}

// GetMaxJSDivergence returns the maximum JS divergence among successful comparisons
func (bc *BatchComparison) GetMaxJSDivergence() float64 {
	successful := bc.GetSuccessfulResults()
	if len(successful) == 0 {
		return 0.0
	}

	max := successful[0].JSDivergence
	for _, result := range successful {
		if result.JSDivergence > max {
			max = result.JSDivergence
		}
	}

	return max
}

// GetMinJSDivergence returns the minimum JS divergence among successful comparisons
func (bc *BatchComparison) GetMinJSDivergence() float64 {
	successful := bc.GetSuccessfulResults()
	if len(successful) == 0 {
		return 0.0
	}

	min := successful[0].JSDivergence
	for _, result := range successful {
		if result.JSDivergence < min {
			min = result.JSDivergence
		}
	}

	return min
}

// IsWithinThreshold checks if all successful comparisons are within the specified threshold
func (bc *BatchComparison) IsWithinThreshold(threshold float64) bool {
	successful := bc.GetSuccessfulResults()
	if len(successful) == 0 {
		return false
	}

	for _, result := range successful {
		if result.JSDivergence > threshold {
			return false
		}
	}

	return true
}

// GetStatistics returns comprehensive statistics about the batch results
func (bc *BatchComparison) GetStatistics() BatchStatistics {
	successful := bc.GetSuccessfulResults()

	stats := BatchStatistics{
		TotalComparisons:      len(bc.results),
		SuccessfulComparisons: len(successful),
		FailedComparisons:     len(bc.results) - len(successful),
	}

	if len(successful) > 0 {
		stats.AverageJSDivergence = bc.GetAverageJSDivergence()
		stats.MaxJSDivergence = bc.GetMaxJSDivergence()
		stats.MinJSDivergence = bc.GetMinJSDivergence()

		// Calculate standard deviation
		mean := stats.AverageJSDivergence
		sumSquares := 0.0
		for _, result := range successful {
			diff := result.JSDivergence - mean
			sumSquares += diff * diff
		}
		n := len(successful)
		if n > 1 {
			stats.StdDevJSDivergence = math.Sqrt(sumSquares / float64(n-1))
		}
	}

	return stats
}

// BatchStatistics contains comprehensive statistics about batch comparison results
type BatchStatistics struct {
	TotalComparisons      int
	SuccessfulComparisons int
	FailedComparisons     int
	AverageJSDivergence   float64
	MaxJSDivergence       float64
	MinJSDivergence       float64
	StdDevJSDivergence    float64
}
