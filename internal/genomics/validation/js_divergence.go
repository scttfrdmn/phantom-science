package validation

import (
	"fmt"
	"math"
	"sort"
)

// JSdivergenceCalculator computes the Jensen-Shannon divergence (base 2,
// bounded [0,1]) between two discrete distributions, each first mixed with a
// small uniform component so zero probabilities stay finite.
type JSdivergenceCalculator struct {
	smoothingFactor float64
	confidenceLevel float64
	minProbability  float64
}

// ValidationResult represents the result of a JS divergence calculation
type ValidationResult struct {
	JSdivergence       float64
	IsSignificant      bool       // divergence above a fixed threshold; see IsSignificant
	ConfidenceInterval [2]float64 // heuristic interval; see calculateConfidenceInterval
	SampleSize         int        // number of distinct categories, not observations
	Smoothed           bool
	Recommendations    []string
}

// NewJSdivergenceCalculator creates a new JS divergence calculator
func NewJSdivergenceCalculator() *JSdivergenceCalculator {
	return &JSdivergenceCalculator{
		smoothingFactor: 0.001, // Small smoothing factor for zero probabilities
		confidenceLevel: 0.95,  // 95% confidence interval
		minProbability:  1e-10, // Minimum probability threshold
	}
}

// SetSmoothingFactor sets the smoothing factor for zero probabilities
func (calc *JSdivergenceCalculator) SetSmoothingFactor(factor float64) {
	calc.smoothingFactor = factor
}

// SetConfidenceLevel sets the confidence level for statistical testing
func (calc *JSdivergenceCalculator) SetConfidenceLevel(level float64) {
	calc.confidenceLevel = level
}

// Calculate computes the Jensen-Shannon divergence between two probability distributions
func (calc *JSdivergenceCalculator) Calculate(p, q map[string]float64) (*ValidationResult, error) {
	if len(p) == 0 || len(q) == 0 {
		return nil, fmt.Errorf("cannot calculate JS divergence for empty distributions")
	}

	// Validate that both distributions are proper probability distributions
	if err := calc.validateDistribution(p); err != nil {
		return nil, fmt.Errorf("invalid distribution P: %w", err)
	}
	if err := calc.validateDistribution(q); err != nil {
		return nil, fmt.Errorf("invalid distribution Q: %w", err)
	}

	// Get all unique keys from both distributions
	allKeys := calc.getAllKeys(p, q)

	// Apply smoothing to handle zero probabilities
	pSmoothed := calc.applySmoothingToDistribution(p, allKeys)
	qSmoothed := calc.applySmoothingToDistribution(q, allKeys)

	// Calculate the average distribution M = 0.5 * (P + Q)
	m := calc.calculateAverageDistribution(pSmoothed, qSmoothed, allKeys)

	// Calculate KL divergences
	klPM := calc.calculateKLDivergence(pSmoothed, m, allKeys)
	klQM := calc.calculateKLDivergence(qSmoothed, m, allKeys)

	// Calculate JS divergence: JS(P||Q) = 0.5 * KL(P||M) + 0.5 * KL(Q||M)
	jsDivergence := 0.5*klPM + 0.5*klQM

	// Determine if smoothing was applied
	smoothed := calc.smoothingFactor > 0

	// Number of categories (the observation count is not known here)
	sampleSize := len(allKeys)

	// Check statistical significance
	isSignificant := calc.IsSignificant(jsDivergence, sampleSize)

	// Heuristic interval, not a bootstrap
	confInterval := calc.calculateConfidenceInterval(jsDivergence, sampleSize)

	// Generate recommendations
	recommendations := calc.generateRecommendations(jsDivergence, isSignificant, smoothed)

	return &ValidationResult{
		JSdivergence:       jsDivergence,
		IsSignificant:      isSignificant,
		ConfidenceInterval: confInterval,
		SampleSize:         sampleSize,
		Smoothed:           smoothed,
		Recommendations:    recommendations,
	}, nil
}

// IsSignificant reports whether divergence exceeds a fixed cutoff: 0.01, or
// 0.015 below 1000 and 0.02 below 100 categories. It is a heuristic cutoff,
// not a hypothesis test, and the thresholds are not taken from a reference.
func (calc *JSdivergenceCalculator) IsSignificant(divergence float64, sampleSize int) bool {
	threshold := 0.01

	// Looser cutoff for fewer categories
	if sampleSize < 100 {
		threshold = 0.02
	} else if sampleSize < 1000 {
		threshold = 0.015
	}

	return divergence > threshold
}

// validateDistribution checks if a distribution is valid (non-negative, sums to ~1)
func (calc *JSdivergenceCalculator) validateDistribution(dist map[string]float64) error {
	if len(dist) == 0 {
		return fmt.Errorf("distribution is empty")
	}

	sum := 0.0
	for key, prob := range dist {
		if prob < 0 {
			return fmt.Errorf("negative probability for key %s: %f", key, prob)
		}
		sum += prob
	}

	// Allow small tolerance for floating point precision
	if math.Abs(sum-1.0) > 1e-6 {
		return fmt.Errorf("distribution does not sum to 1.0: sum = %f", sum)
	}

	return nil
}

// getAllKeys returns all unique keys from both distributions
func (calc *JSdivergenceCalculator) getAllKeys(p, q map[string]float64) []string {
	keySet := make(map[string]bool)

	for key := range p {
		keySet[key] = true
	}
	for key := range q {
		keySet[key] = true
	}

	keys := make([]string, 0, len(keySet))
	for key := range keySet {
		keys = append(keys, key)
	}

	// Sort keys for consistent ordering
	sort.Strings(keys)
	return keys
}

// applySmoothingToDistribution applies smoothing to handle zero probabilities
func (calc *JSdivergenceCalculator) applySmoothingToDistribution(dist map[string]float64, allKeys []string) map[string]float64 {
	smoothed := make(map[string]float64)

	// Calculate uniform distribution for smoothing
	uniformProb := calc.smoothingFactor / float64(len(allKeys))

	// Apply smoothing: P_smooth = (1-α)*P + α*U
	for _, key := range allKeys {
		originalProb := dist[key] // will be 0 if key doesn't exist
		smoothedProb := (1-calc.smoothingFactor)*originalProb + uniformProb
		smoothed[key] = math.Max(smoothedProb, calc.minProbability)
	}

	// Renormalize to ensure sum = 1
	return calc.renormalizeDistribution(smoothed)
}

// calculateAverageDistribution calculates M = 0.5 * (P + Q)
func (calc *JSdivergenceCalculator) calculateAverageDistribution(p, q map[string]float64, allKeys []string) map[string]float64 {
	m := make(map[string]float64)

	for _, key := range allKeys {
		m[key] = 0.5 * (p[key] + q[key])
	}

	return m
}

// calculateKLDivergence calculates KL(P||Q) = sum(P * log(P/Q))
func (calc *JSdivergenceCalculator) calculateKLDivergence(p, q map[string]float64, allKeys []string) float64 {
	klDiv := 0.0

	for _, key := range allKeys {
		pProb := p[key]
		qProb := q[key]

		// Skip if probability is effectively zero
		if pProb < calc.minProbability {
			continue
		}

		// Ensure we don't divide by zero
		if qProb < calc.minProbability {
			qProb = calc.minProbability
		}

		klDiv += pProb * math.Log2(pProb/qProb)
	}

	return klDiv
}

// renormalizeDistribution ensures the distribution sums to 1.0
func (calc *JSdivergenceCalculator) renormalizeDistribution(dist map[string]float64) map[string]float64 {
	sum := 0.0
	for _, prob := range dist {
		sum += prob
	}

	if sum == 0.0 {
		return dist
	}

	normalized := make(map[string]float64)
	for key, prob := range dist {
		normalized[key] = prob / sum
	}

	return normalized
}

// calculateConfidenceInterval calculates a confidence interval for the JS divergence
func (calc *JSdivergenceCalculator) calculateConfidenceInterval(jsDivergence float64, sampleSize int) [2]float64 {
	// Heuristic: SE = sqrt(JS / categories). This is not a derived standard
	// error and the interval has no guaranteed coverage; a bootstrap over the
	// underlying observations would be needed for that.
	standardError := math.Sqrt(jsDivergence / float64(sampleSize))

	// z for the configured level (0.90, 0.99; anything else uses 0.95)
	zScore := 1.96
	if calc.confidenceLevel == 0.99 {
		zScore = 2.576
	} else if calc.confidenceLevel == 0.90 {
		zScore = 1.645
	}

	margin := zScore * standardError

	lowerBound := math.Max(0.0, jsDivergence-margin)
	upperBound := jsDivergence + margin

	return [2]float64{lowerBound, upperBound}
}

// generateRecommendations generates actionable recommendations based on JS divergence
func (calc *JSdivergenceCalculator) generateRecommendations(jsDivergence float64, isSignificant bool, smoothed bool) []string {
	var recommendations []string

	if jsDivergence < 0.001 {
		recommendations = append(recommendations, "Excellent: Distributions are virtually identical")
	} else if jsDivergence < 0.01 {
		recommendations = append(recommendations, "Good: Distributions meet authenticity threshold")
	} else if jsDivergence < 0.05 {
		recommendations = append(recommendations, "Moderate: Consider improving model parameters")
	} else {
		recommendations = append(recommendations, "Poor: Significant differences detected - model revision needed")
	}

	if isSignificant {
		recommendations = append(recommendations, "Divergence exceeds the significance cutoff (heuristic threshold, not a hypothesis test)")
		recommendations = append(recommendations, "Consider increasing model complexity or training data")
	}

	if smoothed {
		recommendations = append(recommendations, "Note: Smoothing applied for zero probabilities")
	}

	if jsDivergence > 0.1 {
		recommendations = append(recommendations, "High divergence: Check for systematic bias in generation")
	}

	return recommendations
}

// CompareMultipleDistributions compares multiple distributions pairwise
func (calc *JSdivergenceCalculator) CompareMultipleDistributions(distributions map[string]map[string]float64) (map[string]map[string]*ValidationResult, error) {
	if len(distributions) < 2 {
		return nil, fmt.Errorf("need at least 2 distributions for comparison")
	}

	results := make(map[string]map[string]*ValidationResult)

	// Get all distribution names
	names := make([]string, 0, len(distributions))
	for name := range distributions {
		names = append(names, name)
	}

	// Compare all pairs
	for i := 0; i < len(names); i++ {
		for j := i + 1; j < len(names); j++ {
			name1, name2 := names[i], names[j]

			result, err := calc.Calculate(distributions[name1], distributions[name2])
			if err != nil {
				return nil, fmt.Errorf("failed to calculate JS divergence between %s and %s: %w", name1, name2, err)
			}

			// Initialize nested map if needed
			if results[name1] == nil {
				results[name1] = make(map[string]*ValidationResult)
			}
			if results[name2] == nil {
				results[name2] = make(map[string]*ValidationResult)
			}

			results[name1][name2] = result
			results[name2][name1] = result // Symmetric
		}
	}

	return results, nil
}

// GetAuthenticityScore maps JS divergence to 0-100 as 100*exp(-100*JS). The
// decay constant is arbitrary; the score is a rescaling, not a probability.
func (calc *JSdivergenceCalculator) GetAuthenticityScore(jsDivergence float64) float64 {
	// Convert JS divergence to a score where 0 divergence = 100 score
	// Using exponential decay: score = 100 * exp(-k * divergence)
	k := 100.0 // Decay constant
	score := 100.0 * math.Exp(-k*jsDivergence)

	// Cap at 100
	if score > 100.0 {
		score = 100.0
	}

	return score
}

// GetQualityGrade returns a quality grade based on JS divergence
func (calc *JSdivergenceCalculator) GetQualityGrade(jsDivergence float64) string {
	if jsDivergence < 0.001 {
		return "Excellent"
	} else if jsDivergence < 0.01 {
		return "Good"
	} else if jsDivergence < 0.05 {
		return "Fair"
	} else if jsDivergence < 0.1 {
		return "Poor"
	} else {
		return "Failed"
	}
}
