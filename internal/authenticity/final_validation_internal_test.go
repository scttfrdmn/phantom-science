//go:build !release

package authenticity

import (
	"math"
	"testing"
)

// makeTestResults builds a slice of IndividualTestResult with controlled scores and labels.
func makeTestResults(scores []float64, labels []bool) []*IndividualTestResult {
	results := make([]*IndividualTestResult, len(scores))
	for i := range scores {
		results[i] = &IndividualTestResult{
			TestResult:    scores[i],
			PassingStatus: labels[i],
		}
	}
	return results
}

var fvsInstance = NewFinalValidationSystem()

// ---------------------------------------------------------------------------
// generateConfusionMatrix
// ---------------------------------------------------------------------------

// TestGenerateConfusionMatrix_Empty returns all-zero matrix for empty input.
func TestGenerateConfusionMatrix_Empty(t *testing.T) {
	cm := fvsInstance.generateConfusionMatrix(nil)
	if cm.TruePositives+cm.TrueNegatives+cm.FalsePositives+cm.FalseNegatives != 0 {
		t.Errorf("empty input: expected all zero, got TP=%d TN=%d FP=%d FN=%d",
			cm.TruePositives, cm.TrueNegatives, cm.FalsePositives, cm.FalseNegatives)
	}
}

// TestGenerateConfusionMatrix_PerfectPositiveClassifier scores=1.0 for all positives,
// 0.0 for all negatives → TP=N_pos, TN=N_neg, FP=FN=0.
func TestGenerateConfusionMatrix_PerfectPositiveClassifier(t *testing.T) {
	scores := []float64{1.0, 1.0, 0.0, 0.0}
	labels := []bool{true, true, false, false}
	cm := fvsInstance.generateConfusionMatrix(makeTestResults(scores, labels))
	if cm.TruePositives != 2 || cm.TrueNegatives != 2 || cm.FalsePositives != 0 || cm.FalseNegatives != 0 {
		t.Errorf("perfect classifier: TP=%d TN=%d FP=%d FN=%d", cm.TruePositives, cm.TrueNegatives, cm.FalsePositives, cm.FalseNegatives)
	}
}

// TestGenerateConfusionMatrix_PerfectNegativeClassifier all scores=0.0 → FN=N_pos, TN=N_neg.
func TestGenerateConfusionMatrix_PerfectNegativeClassifier(t *testing.T) {
	scores := []float64{0.0, 0.0, 0.0, 0.0}
	labels := []bool{true, true, false, false}
	cm := fvsInstance.generateConfusionMatrix(makeTestResults(scores, labels))
	if cm.FalseNegatives != 2 || cm.TrueNegatives != 2 || cm.TruePositives != 0 || cm.FalsePositives != 0 {
		t.Errorf("all-negative: TP=%d TN=%d FP=%d FN=%d", cm.TruePositives, cm.TrueNegatives, cm.FalsePositives, cm.FalseNegatives)
	}
}

// TestGenerateConfusionMatrix_TotalsMatchInput verifies TP+TN+FP+FN == len(results).
func TestGenerateConfusionMatrix_TotalsMatchInput(t *testing.T) {
	scores := []float64{0.9, 0.8, 0.4, 0.3, 0.6, 0.2}
	labels := []bool{true, false, true, false, true, false}
	cm := fvsInstance.generateConfusionMatrix(makeTestResults(scores, labels))
	total := cm.TruePositives + cm.TrueNegatives + cm.FalsePositives + cm.FalseNegatives
	if total != len(scores) {
		t.Errorf("total=%d, want %d", total, len(scores))
	}
}

// TestGenerateConfusionMatrix_ClassificationMatrixNotNil checks the map is initialised.
func TestGenerateConfusionMatrix_ClassificationMatrixNotNil(t *testing.T) {
	cm := fvsInstance.generateConfusionMatrix(makeTestResults([]float64{0.6}, []bool{true}))
	if cm.ClassificationMatrix == nil {
		t.Error("ClassificationMatrix should not be nil")
	}
}

// ---------------------------------------------------------------------------
// generateROCAnalysis
// ---------------------------------------------------------------------------

// TestGenerateROCAnalysis_Empty returns safe defaults.
func TestGenerateROCAnalysis_Empty(t *testing.T) {
	roc := fvsInstance.generateROCAnalysis(nil)
	if roc.AUC < 0 || roc.AUC > 1 {
		t.Errorf("empty: AUC=%v out of [0,1]", roc.AUC)
	}
	if roc.OptimalThreshold <= 0 {
		t.Errorf("empty: OptimalThreshold=%v, want > 0", roc.OptimalThreshold)
	}
}

// TestGenerateROCAnalysis_PerfectClassifier: score=1 for positives, 0 for negatives → AUC≈1.
func TestGenerateROCAnalysis_PerfectClassifier(t *testing.T) {
	scores := []float64{1.0, 1.0, 0.0, 0.0}
	labels := []bool{true, true, false, false}
	roc := fvsInstance.generateROCAnalysis(makeTestResults(scores, labels))
	if roc.AUC < 0.95 {
		t.Errorf("perfect classifier AUC=%v, want ≥ 0.95", roc.AUC)
	}
	if roc.OptimalThreshold <= 0 || roc.OptimalThreshold > 1 {
		t.Errorf("OptimalThreshold=%v out of (0,1]", roc.OptimalThreshold)
	}
}

// TestGenerateROCAnalysis_AUCInRange ensures AUC is always in [0,1] for mixed data.
func TestGenerateROCAnalysis_AUCInRange(t *testing.T) {
	scores := []float64{0.9, 0.8, 0.7, 0.4, 0.3, 0.2, 0.6, 0.5}
	labels := []bool{true, false, true, false, true, false, true, false}
	roc := fvsInstance.generateROCAnalysis(makeTestResults(scores, labels))
	if roc.AUC < 0 || roc.AUC > 1 {
		t.Errorf("mixed: AUC=%v out of [0,1]", roc.AUC)
	}
	if len(roc.ROCCurvePoints) == 0 {
		t.Error("ROCCurvePoints should not be empty")
	}
}

// TestGenerateROCAnalysis_OptimalThresholdPositive verifies threshold > 0.
func TestGenerateROCAnalysis_OptimalThresholdPositive(t *testing.T) {
	scores := []float64{0.9, 0.1, 0.8, 0.2}
	labels := []bool{true, false, true, false}
	roc := fvsInstance.generateROCAnalysis(makeTestResults(scores, labels))
	if roc.OptimalThreshold <= 0 {
		t.Errorf("OptimalThreshold=%v, want > 0", roc.OptimalThreshold)
	}
}

// ---------------------------------------------------------------------------
// generatePrecisionRecallAnalysis
// ---------------------------------------------------------------------------

// TestGeneratePRAnalysis_Empty returns safe defaults.
func TestGeneratePRAnalysis_Empty(t *testing.T) {
	pr := fvsInstance.generatePrecisionRecallAnalysis(nil)
	if pr.AUPRC < 0 || pr.AUPRC > 1 {
		t.Errorf("empty: AUPRC=%v out of [0,1]", pr.AUPRC)
	}
	if pr.OptimalThreshold <= 0 {
		t.Errorf("empty: OptimalThreshold=%v, want > 0", pr.OptimalThreshold)
	}
}

// TestGeneratePRAnalysis_PerfectClassifier: score=1 for positives → high AUPRC.
func TestGeneratePRAnalysis_PerfectClassifier(t *testing.T) {
	scores := []float64{1.0, 1.0, 0.0, 0.0}
	labels := []bool{true, true, false, false}
	pr := fvsInstance.generatePrecisionRecallAnalysis(makeTestResults(scores, labels))
	if pr.AUPRC < 0.9 {
		t.Errorf("perfect classifier AUPRC=%v, want ≥ 0.9", pr.AUPRC)
	}
}

// TestGeneratePRAnalysis_AUPRCInRange ensures AUPRC is in [0,1].
func TestGeneratePRAnalysis_AUPRCInRange(t *testing.T) {
	scores := []float64{0.9, 0.8, 0.7, 0.4, 0.3, 0.2, 0.6, 0.5}
	labels := []bool{true, false, true, false, true, false, true, false}
	pr := fvsInstance.generatePrecisionRecallAnalysis(makeTestResults(scores, labels))
	if pr.AUPRC < 0 || pr.AUPRC > 1 {
		t.Errorf("mixed: AUPRC=%v out of [0,1]", pr.AUPRC)
	}
	if pr.F1AtOptimal < 0 || pr.F1AtOptimal > 1 {
		t.Errorf("F1AtOptimal=%v out of [0,1]", pr.F1AtOptimal)
	}
	if len(pr.PRCurvePoints) == 0 {
		t.Error("PRCurvePoints should not be empty")
	}
}

// TestGeneratePRAnalysis_OptimalThresholdPositive verifies threshold > 0.
func TestGeneratePRAnalysis_OptimalThresholdPositive(t *testing.T) {
	scores := []float64{0.9, 0.1, 0.8, 0.2}
	labels := []bool{true, false, true, false}
	pr := fvsInstance.generatePrecisionRecallAnalysis(makeTestResults(scores, labels))
	if pr.OptimalThreshold <= 0 {
		t.Errorf("OptimalThreshold=%v, want > 0", pr.OptimalThreshold)
	}
}

// TestGeneratePRAnalysis_F1Consistent verifies F1 = 2PR/(P+R) at optimal.
func TestGeneratePRAnalysis_F1Consistent(t *testing.T) {
	scores := []float64{0.9, 0.8, 0.4, 0.3}
	labels := []bool{true, false, true, false}
	pr := fvsInstance.generatePrecisionRecallAnalysis(makeTestResults(scores, labels))
	p := pr.PrecisionAtOptimal
	r := pr.RecallAtOptimal
	if p+r > 0 {
		expectedF1 := 2 * p * r / (p + r)
		if math.Abs(pr.F1AtOptimal-expectedF1) > 0.01 {
			t.Errorf("F1AtOptimal=%v, expected 2*%.2f*%.2f/(%.2f+%.2f)=%.3f",
				pr.F1AtOptimal, p, r, p, r, expectedF1)
		}
	}
}
