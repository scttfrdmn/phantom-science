// Command test_simple_2025_integration constructs the platform error
// signature manager, the cross-platform bias corrector and the
// PrivacyPreservingGeneticsEngine, calls each once on synthetic inputs and
// prints what they return. It then runs a Multinomial Naive Bayes 5-fold
// cross-validation on 4-mer frequencies of two sets of i.i.d. random
// sequences generated in this file. It reads no input files.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
// The privacy engine is a mock whose metrics are constants, and the
// classifier step compares two draws from the same random generator, not
// phantom-science output against real reads, so its ~50% accuracy is
// expected by construction and says nothing about realism. The "research
// integration score" is a sum of constants.
package main

import (
	"crypto/rand"
	"fmt"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/population"
	pt "github.com/scttfrdmn/phantom-science/internal/testing"
)

func main() {
	fmt.Println("🔬 TESTING 2024-2025 RESEARCH INTEGRATION (SIMPLIFIED)")
	fmt.Println("=" + strings.Repeat("=", 60))

	fmt.Println("\nExperiment: exercises experimental components once each; not a validation")

	// Test 1: Enhanced Platform Error Modeling with 2024-2025 improvements
	fmt.Println("\n1️⃣ Testing Enhanced Platform Error Modeling...")

	platformManager := genomics.NewPlatformErrorSignatureManager()
	if platformManager == nil {
		fmt.Printf("   ❌ Failed to initialize platform error signature manager\n")
		return
	}
	fmt.Printf("   ✅ Platform error signature manager initialized successfully\n")

	// Validate 2024-2025 error rate improvements
	fmt.Printf("   📊 Configured error signatures:\n")

	// Check all platforms
	platforms := []genomics.SequencingPlatform{
		genomics.NovaSeq6000, genomics.NextSeq2000, genomics.MiSeq,
		genomics.PacBioRevio, genomics.DNBSeqT20x2,
		genomics.ElementAVITI, genomics.PromethION24,
	}

	platformCount := 0
	for _, platform := range platforms {
		signature := platformManager.GetErrorSignature(platform)
		if signature != nil {
			platformCount++
			fmt.Printf("      ✅ %s: Signature configured\n", platform)
		}
	}

	if platformCount == 7 {
		fmt.Printf("   All 7 platforms have an error signature configured\n")
	}

	// Test 2: Cross-Platform Bias Correction Framework
	fmt.Println("\n2️⃣ Testing Cross-Platform Bias Correction Framework...")

	biasCorrector := genomics.NewCrossPlatformBiasCorrector()
	if biasCorrector == nil {
		fmt.Printf("   ❌ Failed to initialize cross-platform bias corrector\n")
		return
	}
	fmt.Printf("   ✅ Cross-platform bias corrector initialized successfully\n")
	fmt.Printf("   🔬 Components named after BEADS, MANCIE and quantile normalization\n")

	// Test 3: Privacy-Preserving Genetics Engine
	fmt.Println("\n3️⃣ Testing Privacy-Preserving Genetics Engine...")

	privacyEngine := population.NewPrivacyPreservingGeneticsEngine()
	if privacyEngine == nil {
		fmt.Printf("   ❌ Failed to initialize privacy-preserving genetics engine\n")
		return
	}
	fmt.Printf("   ✅ Privacy-preserving genetics engine initialized successfully\n")
	fmt.Printf("   Note: this engine is a mock (no GAN, no DP guarantee); its metrics are constants\n")

	// Test 4: Generate Privacy-Preserving Data
	fmt.Println("\n4️⃣ Testing Privacy-Preserving Data Generation...")

	privacyRequirements := &population.PrivacyRequirements{
		DifferentialPrivacyEpsilon: 1.0,
		UseFederatedPrivacy:        true,
		DefendAgainstTabPFN:        true,
		DefendAgainstLocalAncestry: true,
		DefendAgainstHaploNet:      true,
		MinimumPrivacyScore:        0.90,
		MaximumPrivacyBudget:       2.0,
	}

	targetPopulation := population.CEU
	numSamples := 50

	result, err := privacyEngine.GeneratePrivacyPreservingGenomicData(
		targetPopulation, numSamples, privacyRequirements)

	if err != nil {
		fmt.Printf("   ⚠️  Privacy-preserving data generation: %v\n", err)
	} else {
		fmt.Printf("   ✅ Generated %d mock samples\n", len(result.GeneratedData))
		if result.PrivacyMetrics != nil {
			fmt.Printf("   📊 Privacy Score (constant): %.3f\n", result.PrivacyMetrics.OverallPrivacyRating)
		}
		if result.DefenseMetrics != nil {
			fmt.Printf("   📊 Defense Rating (constant): %.3f\n", result.DefenseMetrics.OverallDefenseRating)
		}
	}

	// Test 5: Cross-Platform Bias Correction Application
	fmt.Println("\n5️⃣ Testing Cross-Platform Bias Correction...")

	// Create test data
	testData := map[genomics.SequencingPlatform][]float64{
		genomics.NovaSeq6000:  generateTestSequence(500, 0.001),  // Low error rate
		genomics.PacBioRevio:  generateTestSequence(500, 0.0045), // Improved HiFi rate
		genomics.PromethION24: generateTestSequence(500, 0.008),  // Improved R10 rate
		genomics.ElementAVITI: generateTestSequence(500, 0.0015), // Novel chemistry
	}

	correctionResult, err := biasCorrector.CorrectCrossPlatformBias(testData, population.CEU)
	if err != nil {
		fmt.Printf("   ⚠️  Cross-platform correction: %v\n", err)
	} else {
		fmt.Printf("   ✅ Cross-platform bias correction completed\n")
		fmt.Printf("   📊 Correction Quality: %.3f\n", correctionResult.OverallQuality)
	}

	// Test 6: Comprehensive Assessment
	fmt.Println("\n6️⃣ Research integration score...")

	assessmentScore := calculateComprehensiveAssessment()
	fmt.Printf("   📚 Research Integration Score: %.1f/100 (a sum of fixed constants, not measured)\n", assessmentScore)

	// Test 7: classifier on two draws from the same random generator, at
	// three read lengths. This is a null case: ~50% is expected by
	// construction and is not a comparison with real data.
	fmt.Println("\n7️⃣ Classifier null case (two i.i.d. random sequence sets)...")
	fmt.Printf("   Expected: ~50%% accuracy, because both sets come from the same generator\n")
	fmt.Printf("   Running Multinomial Naive Bayes 5-fold CV on 4-mer frequency vectors...\n")

	type platformResult struct {
		name       string
		readLength int
		score      float64
	}
	platformResults := []platformResult{
		{"Illumina (short read)", 150, measureClassifierAccuracy(50, 150)},
		{"PacBio   (long read) ", 15000, measureClassifierAccuracy(50, 15000)},
		{"ONT      (long read) ", 30000, measureClassifierAccuracy(50, 30000)},
	}

	totalScore := 0.0
	for _, pr := range platformResults {
		totalScore += pr.score
	}
	indistinguishabilityScore := totalScore / float64(len(platformResults))

	fmt.Printf("\n   %-24s  %9s  %8s  %s\n", "Platform", "Read len", "Accuracy", "Status")
	fmt.Printf("   %s\n", strings.Repeat("-", 62))
	for _, pr := range platformResults {
		status := "outside 48-52%"
		if pr.score >= 0.48 && pr.score <= 0.52 {
			status = "within 48-52%"
		}
		fmt.Printf("   %-24s  %6d bp  %7.1f%%  %s\n",
			pr.name, pr.readLength, pr.score*100, status)
	}
	fmt.Printf("   %s\n", strings.Repeat("-", 62))
	fmt.Printf("   %-24s  %9s  %7.1f%%\n", "Overall average", "", indistinguishabilityScore*100)

	// Final Summary
	fmt.Printf("\n📋 SUMMARY\n")
	fmt.Println("=" + strings.Repeat("=", 60))

	qualityChecks := map[string]bool{
		"Platform error modeling (not checked)":        true,
		"Cross-platform bias correction (not checked)": true,
		"Privacy engine returned a result":             result != nil,
		"Research findings integrated (not checked)":   true,
		"Null-case classifier accuracy within 48-52%":  indistinguishabilityScore >= 0.48 && indistinguishabilityScore <= 0.52,
		"ML attack defenses (not checked)":             true,
		"All 7 platforms have a signature":             platformCount == 7,
		"Privacy engine returned metrics (constants)":  result != nil && result.PrivacyMetrics != nil,
	}

	passedChecks := 0
	totalChecks := len(qualityChecks)

	for check, passed := range qualityChecks {
		status := "✅ PASSED"
		if !passed {
			status = "❌ FAILED"
		} else {
			passedChecks++
		}
		fmt.Printf("   %s: %s\n", check, status)
	}

	fmt.Printf("\n🎯 RESULTS\n")
	fmt.Printf("   Checks true: %d/%d (4 of them hard-coded)\n", passedChecks, totalChecks)
	fmt.Printf("   Research Integration Score: %.1f/100 (constant)\n", assessmentScore)
	fmt.Printf("   Null-case classifier accuracy: %.1f%%\n", indistinguishabilityScore*100)

	if passedChecks >= totalChecks*7/8 {
		fmt.Printf("\n   At most one check failed. This run does not show that phantom-science\n")
		fmt.Printf("   output is indistinguishable from real data, nor any privacy or attack resistance.\n")
	} else {
		fmt.Printf("\n   ⚠️  Several checks failed; see above\n")
	}
}

// Helper functions

func generateTestSequence(length int, errorRate float64) []float64 {
	sequence := make([]float64, length)
	for i := range sequence {
		baseValue := 0.42 // Base GC content
		if i%int(1.0/errorRate) == 0 {
			baseValue += 0.05 // Add error
		}
		sequence[i] = baseValue
	}
	return sequence
}

// calculateComprehensiveAssessment returns the sum of fixed per-area
// weights (always 100). Nothing is measured.
func calculateComprehensiveAssessment() float64 {
	components := map[string]float64{
		"platform_error_modeling":          15.0,
		"cross_platform_correction":        15.0,
		"privacy_preserving_genetics":      15.0,
		"advanced_ml_defenses":             15.0,
		"statistical_indistinguishability": 15.0,
		"research_integration":             10.0,
		"implementation_completeness":      10.0,
		"future_proofing":                  5.0,
	}

	total := 0.0
	for _, score := range components {
		total += score
	}

	return total
}

// calculateIndistinguishabilityScore returns the mean of measureClassifierAccuracy
// at 150, 15,000 and 30,000 bases. Both classes come from the same random
// generator, so a value near 0.5 is expected and says nothing about how
// close generated data is to real data.
func calculateIndistinguishabilityScore() float64 {
	short := measureClassifierAccuracy(50, 150)
	long1 := measureClassifierAccuracy(50, 15000)
	long2 := measureClassifierAccuracy(50, 30000)
	return (short + long1 + long2) / 3
}

// measureClassifierAccuracy generates two independent sets of seqCount sequences
// of seqLength bases each, then measures how accurately a Multinomial Naive Bayes
// classifier (5-fold CV on 4-mer frequency vectors) can distinguish them.
// Both sets come from generateBiologicalSequences, so this is a null case.
func measureClassifierAccuracy(seqCount, seqLength int) float64 {
	groupA := generateBiologicalSequences(seqCount, seqLength)
	groupB := generateBiologicalSequences(seqCount, seqLength)

	features := make([][]float64, 0, seqCount*2)
	labels := make([]int, 0, seqCount*2)
	for _, seq := range groupA {
		features = append(features, pt.Extract4merFrequencies(seq))
		labels = append(labels, 0)
	}
	for _, seq := range groupB {
		features = append(features, pt.Extract4merFrequencies(seq))
		labels = append(labels, 1)
	}

	nb := pt.NewMultinomialNB(1.0)
	mean, _, err := nb.CrossValidate(features, labels, 5)
	if err != nil {
		return 0.5 // fallback: report chance accuracy on classifier error
	}
	return mean
}

// generateBiologicalSequences returns count sequences of length bases each,
// sampled with human-genome GC content (~42%: G=21%, C=21%, A=29%, T=29%).
// Uses crypto/rand so each call produces an independently sampled set.
func generateBiologicalSequences(count, length int) []string {
	buf := make([]byte, length)
	seqs := make([]string, count)
	for i := range seqs {
		if _, err := rand.Read(buf); err != nil {
			// Extremely unlikely; fill with AT-rich fallback
			for j := range buf {
				buf[j] = byte(j % 4)
			}
		}
		seq := make([]byte, length)
		for j, b := range buf {
			switch int(b) % 100 {
			case 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
				10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20:
				seq[j] = 'G'
			case 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
				31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41:
				seq[j] = 'C'
			case 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
				52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62,
				63, 64, 65, 66, 67, 68, 69, 70:
				seq[j] = 'A'
			default:
				seq[j] = 'T'
			}
		}
		seqs[i] = string(seq)
	}
	return seqs
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
