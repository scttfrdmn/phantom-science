//go:build !release
// +build !release

package main

import (
	"runtime"
	"testing"
	"time"
)

func TestGetStatusEmoji_Success(t *testing.T) {
	if got := getStatusEmoji(true); got != "✅" {
		t.Errorf("getStatusEmoji(true) = %q, want ✅", got)
	}
}

func TestGetStatusEmoji_Failure(t *testing.T) {
	if got := getStatusEmoji(false); got != "❌" {
		t.Errorf("getStatusEmoji(false) = %q, want ❌", got)
	}
}

func TestAbs_Positive(t *testing.T) {
	if got := abs(5.0); got != 5.0 {
		t.Errorf("abs(5.0) = %f, want 5.0", got)
	}
}

func TestAbs_Negative(t *testing.T) {
	if got := abs(-5.0); got != 5.0 {
		t.Errorf("abs(-5.0) = %f, want 5.0", got)
	}
}

func TestAbs_Zero(t *testing.T) {
	if got := abs(0.0); got != 0.0 {
		t.Errorf("abs(0.0) = %f, want 0.0", got)
	}
}

func TestDetectAppleSiliconCores_Positive(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

// ---------- trainPopulationModel ----------

func TestTrainPopulationModel_YRI(t *testing.T) {
	target := PopulationTarget{Population: "YRI", TargetGC: 0.47, Tolerance: 0.02}
	model := trainPopulationModel(target)
	if model == nil {
		t.Fatal("trainPopulationModel(YRI) returned nil")
	}
	if model.Population != "YRI" {
		t.Errorf("Population = %q, want \"YRI\"", model.Population)
	}
	if !model.Trained {
		t.Error("expected Trained=true")
	}
	if model.Complexity <= 0 {
		t.Errorf("Complexity = %f, want > 0", model.Complexity)
	}
}

func TestTrainPopulationModel_CEU(t *testing.T) {
	target := PopulationTarget{Population: "CEU", TargetGC: 0.45, Tolerance: 0.02}
	model := trainPopulationModel(target)
	if model == nil {
		t.Fatal("trainPopulationModel(CEU) returned nil")
	}
	if model.Complexity != 0.6844 {
		t.Errorf("CEU Complexity = %f, want 0.6844", model.Complexity)
	}
}

func TestTrainPopulationModel_CHB(t *testing.T) {
	target := PopulationTarget{Population: "CHB", TargetGC: 0.42, Tolerance: 0.02}
	model := trainPopulationModel(target)
	if model == nil {
		t.Fatal("trainPopulationModel(CHB) returned nil")
	}
	if model.Complexity != 0.6900 {
		t.Errorf("CHB Complexity = %f, want 0.6900", model.Complexity)
	}
}

func TestTrainPopulationModel_Global(t *testing.T) {
	target := PopulationTarget{Population: "Global", TargetGC: 0.44, Tolerance: 0.02}
	model := trainPopulationModel(target)
	if model.Complexity != 0.6627 {
		t.Errorf("Global Complexity = %f, want 0.6627", model.Complexity)
	}
}

func TestTrainPopulationModel_Default(t *testing.T) {
	target := PopulationTarget{Population: "UNKNOWN_POP", TargetGC: 0.44, Tolerance: 0.02}
	model := trainPopulationModel(target)
	if model.Complexity != 0.6677 {
		t.Errorf("Default Complexity = %f, want 0.6677", model.Complexity)
	}
}

// ---------- getTotalTrainingTime ----------

func TestGetTotalTrainingTime_Empty(t *testing.T) {
	total := getTotalTrainingTime(map[string]*PopulationModel{})
	if total != 0 {
		t.Errorf("getTotalTrainingTime({}) = %v, want 0", total)
	}
}

func TestGetTotalTrainingTime_WithModels(t *testing.T) {
	models := map[string]*PopulationModel{
		"YRI": {Population: "YRI", TrainTime: time.Second * 5},
		"CEU": {Population: "CEU", TrainTime: time.Second * 3},
	}
	total := getTotalTrainingTime(models)
	if total != time.Second*8 {
		t.Errorf("getTotalTrainingTime = %v, want 8s", total)
	}
}

func TestGetTotalTrainingTime_SingleModel(t *testing.T) {
	models := map[string]*PopulationModel{
		"CHB": {Population: "CHB", TrainTime: time.Millisecond * 500},
	}
	total := getTotalTrainingTime(models)
	if total != time.Millisecond*500 {
		t.Errorf("getTotalTrainingTime = %v, want 500ms", total)
	}
}

func TestGetTotalTrainingTime_ZeroDuration(t *testing.T) {
	models := map[string]*PopulationModel{
		"YRI": {Population: "YRI", TrainTime: 0},
		"CEU": {Population: "CEU", TrainTime: 0},
	}
	total := getTotalTrainingTime(models)
	if total != 0 {
		t.Errorf("getTotalTrainingTime = %v, want 0", total)
	}
}

// ---------- PopulationModel struct ----------

func TestPopulationModel_Fields(t *testing.T) {
	target := PopulationTarget{Population: "YRI", TargetGC: 0.47, Tolerance: 0.02, Rationale: "test"}
	model := trainPopulationModel(target)

	if model.Population != "YRI" {
		t.Errorf("Population = %q, want YRI", model.Population)
	}
	if model.TargetGC != 0.47 {
		t.Errorf("TargetGC = %f, want 0.47", model.TargetGC)
	}
	if model.Model == nil {
		t.Error("Model should not be nil after training")
	}
	if !model.Trained {
		t.Error("Trained should be true after trainPopulationModel")
	}
}

func TestPopulationTarget_Fields(t *testing.T) {
	pt := PopulationTarget{
		Population: "CEU",
		TargetGC:   0.45,
		Tolerance:  0.02,
		Rationale:  "European",
	}
	if pt.Population != "CEU" {
		t.Errorf("Population = %q, want CEU", pt.Population)
	}
	if pt.TargetGC != 0.45 {
		t.Errorf("TargetGC = %f, want 0.45", pt.TargetGC)
	}
	if pt.Tolerance != 0.02 {
		t.Errorf("Tolerance = %f, want 0.02", pt.Tolerance)
	}
	if pt.Rationale != "European" {
		t.Errorf("Rationale = %q, want European", pt.Rationale)
	}
}

// ---------- abs edge cases ----------

func TestAbs_LargePositive(t *testing.T) {
	got := abs(1e9)
	if got != 1e9 {
		t.Errorf("abs(1e9) = %f, want 1e9", got)
	}
}

func TestAbs_LargeNegative(t *testing.T) {
	got := abs(-1e9)
	if got != 1e9 {
		t.Errorf("abs(-1e9) = %f, want 1e9", got)
	}
}

func TestAbs_SmallFraction(t *testing.T) {
	got := abs(-0.001)
	if got != 0.001 {
		t.Errorf("abs(-0.001) = %f, want 0.001", got)
	}
}

// ---------- getStatusEmoji table-driven ----------

func TestGetStatusEmoji_TableDriven(t *testing.T) {
	cases := []struct {
		success bool
		want    string
	}{
		{true, "✅"},
		{false, "❌"},
	}
	for _, tc := range cases {
		got := getStatusEmoji(tc.success)
		if got != tc.want {
			t.Errorf("getStatusEmoji(%v) = %q, want %q", tc.success, got, tc.want)
		}
	}
}

// ---------- detectAppleSiliconCores boundary ----------

func TestDetectAppleSiliconCores_ReturnsPositive(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestDetectAppleSiliconCores_NotExceedTotal(t *testing.T) {
	totalCPU := runtime.NumCPU()
	got := detectAppleSiliconCores()
	if got > totalCPU {
		t.Errorf("detectAppleSiliconCores() = %d exceeds NumCPU() = %d", got, totalCPU)
	}
}

func TestDetectAppleSiliconCores_AtLeastOne(t *testing.T) {
	// The result must always be at least 1 regardless of architecture
	got := detectAppleSiliconCores()
	if got < 1 {
		t.Errorf("detectAppleSiliconCores() = %d, expected at least 1", got)
	}
}

// ---------- trainPopulationModel completeness ----------

func TestTrainPopulationModel_ReturnedModelHasExpectedComplexity(t *testing.T) {
	tests := []struct {
		pop             string
		wantComplexity  float64
	}{
		{"YRI", 0.6772},
		{"CEU", 0.6844},
		{"CHB", 0.6900},
		{"Global", 0.6627},
		{"OTHER", 0.6677},
	}
	for _, tt := range tests {
		t.Run(tt.pop, func(t *testing.T) {
			target := PopulationTarget{Population: tt.pop, TargetGC: 0.45, Tolerance: 0.02}
			model := trainPopulationModel(target)
			if model == nil {
				t.Fatalf("trainPopulationModel(%q) returned nil", tt.pop)
			}
			if model.Complexity != tt.wantComplexity {
				t.Errorf("Complexity = %f, want %f", model.Complexity, tt.wantComplexity)
			}
		})
	}
}

func TestTrainPopulationModel_TrainedFlagAlwaysTrue(t *testing.T) {
	populations := []string{"YRI", "CEU", "CHB", "Global", "UNKNOWN"}
	for _, pop := range populations {
		t.Run(pop, func(t *testing.T) {
			target := PopulationTarget{Population: pop, TargetGC: 0.45, Tolerance: 0.02}
			model := trainPopulationModel(target)
			if model == nil {
				t.Fatalf("trainPopulationModel(%q) returned nil", pop)
			}
			if !model.Trained {
				t.Errorf("population %q: Trained = false, want true", pop)
			}
		})
	}
}

func TestTrainPopulationModel_PopulationFieldPreserved(t *testing.T) {
	target := PopulationTarget{Population: "CHB", TargetGC: 0.46, Tolerance: 0.02}
	model := trainPopulationModel(target)
	if model.Population != "CHB" {
		t.Errorf("Population = %q, want CHB", model.Population)
	}
}

func TestTrainPopulationModel_TargetGCPreserved(t *testing.T) {
	target := PopulationTarget{Population: "Global", TargetGC: 0.42, Tolerance: 0.02}
	model := trainPopulationModel(target)
	if model.TargetGC != 0.42 {
		t.Errorf("TargetGC = %f, want 0.42", model.TargetGC)
	}
}
