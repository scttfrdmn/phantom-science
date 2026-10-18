//go:build !release
// +build !release

package main

import (
	"testing"
)

func TestGenerateTestData_Length(t *testing.T) {
	data := generateTestData(100, 0.05)
	if len(data) != 100 {
		t.Errorf("generateTestData(100, 0.05) len = %d, want 100", len(data))
	}
}

func TestGenerateTestData_ZeroError(t *testing.T) {
	data := generateTestData(50, 0.0)
	for i, v := range data {
		if v != 0.5 {
			t.Errorf("data[%d] = %f, want 0.5 with zero error rate", i, v)
		}
	}
}

func TestCalculateIntegrationScore_Positive(t *testing.T) {
	score := calculateIntegrationScore()
	if score <= 0 {
		t.Errorf("calculateIntegrationScore() = %f, want > 0", score)
	}
}

func TestDetermineCertificationLevel_Diamond(t *testing.T) {
	if got := determineCertificationLevel(97.0); got != "Diamond" {
		t.Errorf("determineCertificationLevel(97.0) = %q, want \"Diamond\"", got)
	}
}

func TestDetermineCertificationLevel_Platinum(t *testing.T) {
	if got := determineCertificationLevel(92.0); got != "Platinum" {
		t.Errorf("determineCertificationLevel(92.0) = %q, want \"Platinum\"", got)
	}
}

func TestDetermineCertificationLevel_Gold(t *testing.T) {
	if got := determineCertificationLevel(87.0); got != "Gold" {
		t.Errorf("determineCertificationLevel(87.0) = %q, want \"Gold\"", got)
	}
}

func TestDetermineCertificationLevel_Silver(t *testing.T) {
	if got := determineCertificationLevel(77.0); got != "Silver" {
		t.Errorf("determineCertificationLevel(77.0) = %q, want \"Silver\"", got)
	}
}

func TestDetermineCertificationLevel_Bronze(t *testing.T) {
	if got := determineCertificationLevel(67.0); got != "Bronze" {
		t.Errorf("determineCertificationLevel(67.0) = %q, want \"Bronze\"", got)
	}
}

func TestDetermineCertificationLevel_None(t *testing.T) {
	if got := determineCertificationLevel(50.0); got != "None" {
		t.Errorf("determineCertificationLevel(50.0) = %q, want \"None\"", got)
	}
}

func TestCalculateThreatCompliance_NotEmpty(t *testing.T) {
	compliance := calculateThreatCompliance()
	if len(compliance) == 0 {
		t.Error("calculateThreatCompliance() returned empty map")
	}
}

func TestCalculateThreatCompliance_ValidRange(t *testing.T) {
	compliance := calculateThreatCompliance()
	for k, v := range compliance {
		if v < 0.0 || v > 1.0 {
			t.Errorf("compliance[%q] = %f, want in [0,1]", k, v)
		}
	}
}
