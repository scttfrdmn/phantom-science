package tui

import (
	"testing"
)

func TestNewStyles(t *testing.T) {
	styles := NewStyles()

	// Test that styles object is created by calling a method
	result := styles.CompressionLabel(100)
	if result == "" {
		t.Error("NewStyles should return a configured Styles object")
	}
}

func TestStylesProgressBarView(t *testing.T) {
	styles := NewStyles()

	// Test progress bar rendering with different percentages
	testCases := []struct {
		percent float64
		width   int
	}{
		{0.0, 20},
		{50.0, 20},
		{100.0, 20},
		{25.5, 30},
	}

	for _, tc := range testCases {
		result := styles.ProgressBarView(tc.percent, tc.width)
		if result == "" {
			t.Errorf("ProgressBarView should return non-empty string for percent=%.1f, width=%d", tc.percent, tc.width)
		}

		// Result should have some reasonable length (not exact due to styling)
		if len(result) < tc.width/2 {
			t.Errorf("ProgressBarView result seems too short: %s", result)
		}
	}
}

func TestStylesCompressionLabel(t *testing.T) {
	styles := NewStyles()

	testCases := []int{1, 10, 100, 1000}

	for _, compression := range testCases {
		result := styles.CompressionLabel(compression)
		if result == "" {
			t.Errorf("CompressionLabel should return non-empty string for compression=%d", compression)
		}

		// Should contain the compression value in some form
		// (We can't check exact format due to styling, but should not be empty)
	}
}

func TestStylesStatusLabel(t *testing.T) {
	styles := NewStyles()

	// Test all documented status values
	testCases := []string{
		"Ready",
		"Running", "Initialization", "Cluster Generation", "Sequencing", "Base Calling",
		"Complete",
		"unknown-default",
	}

	for _, status := range testCases {
		result := styles.StatusLabel(status)
		if result == "" {
			t.Errorf("StatusLabel should return non-empty string for status=%s", status)
		}
	}
}

func TestStylesSequencerIcon(t *testing.T) {
	styles := NewStyles()

	testCases := []string{
		"NovaSeq 6000",
		"NextSeq 2000",
		"MiSeq",
		"Unknown Model",
	}

	for _, model := range testCases {
		result := styles.SequencerIcon(model)
		if result == "" {
			t.Errorf("SequencerIcon should return non-empty string for model=%s", model)
		}

		// Should return some kind of icon or representation
		// The exact format depends on implementation but shouldn't be empty
	}
}

func TestStylesConsistency(t *testing.T) {
	// Test that creating multiple Styles objects produces consistent results
	styles1 := NewStyles()
	styles2 := NewStyles()

	// Same inputs should produce same outputs
	progress1 := styles1.ProgressBarView(50.0, 20)
	progress2 := styles2.ProgressBarView(50.0, 20)

	if progress1 != progress2 {
		t.Error("NewStyles should produce consistent results")
	}

	compression1 := styles1.CompressionLabel(100)
	compression2 := styles2.CompressionLabel(100)

	if compression1 != compression2 {
		t.Error("CompressionLabel should produce consistent results")
	}
}

func TestStylesEdgeCases(t *testing.T) {
	styles := NewStyles()

	// Test edge cases for progress bar
	result := styles.ProgressBarView(-1.0, 10) // Negative percentage
	if result == "" {
		t.Error("ProgressBarView should handle negative percentage gracefully")
	}

	result = styles.ProgressBarView(150.0, 10) // Over 100%
	if result == "" {
		t.Error("ProgressBarView should handle percentage over 100 gracefully")
	}

	result = styles.ProgressBarView(50.0, 0) // Zero width
	if result == "" {
		t.Error("ProgressBarView should handle zero width gracefully")
	}

	// Test empty strings - may return empty, that's okay
	styles.StatusLabel("")

	styles.SequencerIcon("")
}
