package tui

import (
	"os"
	"path/filepath"
	"testing"
)

func TestNewSampleManager(t *testing.T) {
	manager := NewSampleManager()
	if manager == nil {
		t.Fatal("NewSampleManager should not return nil")
	}
	
	if manager.samples == nil {
		t.Error("SampleManager should initialize samples slice")
	}
	
	if len(manager.samples) != 0 {
		t.Error("SampleManager should start with empty samples")
	}
}

func TestSampleManagerAddSample(t *testing.T) {
	manager := NewSampleManager()
	
	sample := Sample{
		ID:   "Test_Sample",
		Lane: 1,
	}
	
	manager.AddSample(sample)
	
	if len(manager.samples) != 1 {
		t.Errorf("Expected 1 sample, got %d", len(manager.samples))
	}
	
	if manager.samples[0].ID != "Test_Sample" {
		t.Errorf("Expected sample ID 'Test_Sample', got %s", manager.samples[0].ID)
	}
	
	if manager.samples[0].Lane != 1 {
		t.Errorf("Expected lane 1, got %d", manager.samples[0].Lane)
	}
}

func TestSampleManagerGetSamples(t *testing.T) {
	manager := NewSampleManager()
	
	samples := manager.GetSamples()
	if len(samples) != 0 {
		t.Error("New manager should have no samples")
	}
	
	// Add some samples
	manager.AddSample(Sample{ID: "Sample1", Lane: 1})
	manager.AddSample(Sample{ID: "Sample2", Lane: 2})
	
	samples = manager.GetSamples()
	if len(samples) != 2 {
		t.Errorf("Expected 2 samples, got %d", len(samples))
	}
}

func TestSampleManagerGetSampleCount(t *testing.T) {
	manager := NewSampleManager()
	
	if manager.GetSampleCount() != 0 {
		t.Error("New manager should have 0 samples")
	}
	
	manager.AddSample(Sample{ID: "Sample1", Lane: 1})
	if manager.GetSampleCount() != 1 {
		t.Errorf("Expected 1 sample, got %d", manager.GetSampleCount())
	}
	
	manager.AddSample(Sample{ID: "Sample2", Lane: 2})
	if manager.GetSampleCount() != 2 {
		t.Errorf("Expected 2 samples, got %d", manager.GetSampleCount())
	}
}

func TestSampleManagerClearAll(t *testing.T) {
	manager := NewSampleManager()
	
	// Add some samples
	manager.AddSample(Sample{ID: "Sample1", Lane: 1})
	manager.AddSample(Sample{ID: "Sample2", Lane: 2})
	
	if manager.GetSampleCount() != 2 {
		t.Errorf("Expected 2 samples before clear, got %d", manager.GetSampleCount())
	}
	
	manager.ClearAll()
	
	if manager.GetSampleCount() != 0 {
		t.Errorf("Expected 0 samples after clear, got %d", manager.GetSampleCount())
	}
}

func TestSampleManagerExportToCSV(t *testing.T) {
	manager := NewSampleManager()
	tmpDir := t.TempDir()
	
	// Add test samples
	manager.AddSample(Sample{ID: "Sample1", Lane: 1})
	manager.AddSample(Sample{ID: "Sample2", Lane: 2})
	
	csvFile := filepath.Join(tmpDir, "test_samples.csv")
	err := manager.ExportToCSV(csvFile)
	
	if err != nil {
		t.Errorf("ExportToCSV should not error: %v", err)
	}
	
	// Check file exists
	if _, err := os.Stat(csvFile); os.IsNotExist(err) {
		t.Error("CSV file should be created")
	}
	
	// Check file content
	content, err := os.ReadFile(csvFile)
	if err != nil {
		t.Errorf("Should be able to read CSV file: %v", err)
	}
	
	contentStr := string(content)
	if !containsString(contentStr, "Sample1") {
		t.Error("CSV should contain Sample1")
	}
	
	if !containsString(contentStr, "Sample2") {
		t.Error("CSV should contain Sample2")
	}
}

func TestSampleManagerImportFromCSV(t *testing.T) {
	manager := NewSampleManager()
	tmpDir := t.TempDir()
	
	// Create test CSV content (6 columns: ID, Lane, Index1, Index2, Project, Description)
	csvContent := `Sample_ID,Lane,Index1,Index2,Project,Description
Test_Sample1,1,ATCG,TACG,Project1,Description1
Test_Sample2,2,GCTA,CTAG,Project2,Description2`
	
	csvFile := filepath.Join(tmpDir, "test_import.csv")
	err := os.WriteFile(csvFile, []byte(csvContent), 0600)
	if err != nil {
		t.Fatalf("Failed to create test CSV: %v", err)
	}
	
	err = manager.ImportFromCSV(csvFile)
	if err != nil {
		t.Errorf("ImportFromCSV should not error: %v", err)
	}
	
	if manager.GetSampleCount() != 2 {
		t.Errorf("Expected 2 imported samples, got %d", manager.GetSampleCount())
	}
	
	samples := manager.GetSamples()
	if len(samples) > 0 && samples[0].ID != "Test_Sample1" {
		t.Errorf("Expected first sample ID 'Test_Sample1', got %s", samples[0].ID)
	}
}

func TestSampleManagerImportFromInvalidCSV(t *testing.T) {
	manager := NewSampleManager()
	tmpDir := t.TempDir()
	
	// Create invalid CSV content
	invalidCSV := `Invalid,CSV,Content
No,Lane,Column`
	
	csvFile := filepath.Join(tmpDir, "invalid.csv")
	err := os.WriteFile(csvFile, []byte(invalidCSV), 0600)
	if err != nil {
		t.Fatalf("Failed to create invalid CSV: %v", err)
	}
	
	err = manager.ImportFromCSV(csvFile)
	// Should handle gracefully - either error or skip invalid rows
	
	// Original sample count should be unchanged if error
	if err != nil && manager.GetSampleCount() != 0 {
		t.Error("Sample count should remain 0 if import fails")
	}
}

func TestSampleManagerImportFromNonexistentFile(t *testing.T) {
	manager := NewSampleManager()
	
	err := manager.ImportFromCSV("/nonexistent/file.csv")
	if err == nil {
		t.Error("ImportFromCSV should error for nonexistent file")
	}
}

func TestSampleManagerRemoveSample(t *testing.T) {
	manager := NewSampleManager()
	manager.AddSample(Sample{ID: "Sample1", Lane: 1})
	manager.AddSample(Sample{ID: "Sample2", Lane: 2})
	manager.AddSample(Sample{ID: "Sample3", Lane: 3})

	// Remove middle element
	manager.RemoveSample(1)

	if manager.GetSampleCount() != 2 {
		t.Errorf("Expected 2 samples after removal, got %d", manager.GetSampleCount())
	}

	samples := manager.GetSamples()
	if samples[0].ID != "Sample1" {
		t.Errorf("Expected first sample 'Sample1', got '%s'", samples[0].ID)
	}
	if samples[1].ID != "Sample3" {
		t.Errorf("Expected second sample 'Sample3', got '%s'", samples[1].ID)
	}

	// Try to remove out-of-bounds index (should not panic)
	manager.RemoveSample(100)
	if manager.GetSampleCount() != 2 {
		t.Error("Count should not change for out-of-bounds removal")
	}

	// Try to remove negative index (should not panic)
	manager.RemoveSample(-1)
	if manager.GetSampleCount() != 2 {
		t.Error("Count should not change for negative index removal")
	}
}

func TestSampleManagerUpdateSample(t *testing.T) {
	manager := NewSampleManager()
	manager.AddSample(Sample{ID: "Sample1", Lane: 1})
	manager.AddSample(Sample{ID: "Sample2", Lane: 2})

	// Update first sample
	updatedSample := Sample{ID: "UpdatedSample1", Lane: 1, Index1: "ATCGATCG"}
	manager.UpdateSample(0, updatedSample)

	samples := manager.GetSamples()
	if samples[0].ID != "UpdatedSample1" {
		t.Errorf("Expected updated ID 'UpdatedSample1', got '%s'", samples[0].ID)
	}
	if samples[0].Index1 != "ATCGATCG" {
		t.Errorf("Expected Index1 'ATCGATCG', got '%s'", samples[0].Index1)
	}

	// Second sample unchanged
	if samples[1].ID != "Sample2" {
		t.Errorf("Expected second sample unchanged as 'Sample2', got '%s'", samples[1].ID)
	}

	// Update out-of-bounds index (should not panic)
	manager.UpdateSample(100, Sample{ID: "NoOp"})
	if manager.GetSampleCount() != 2 {
		t.Error("Count should not change for out-of-bounds update")
	}
}

func TestSampleManagerValidateSamples(t *testing.T) {
	manager := NewSampleManager()

	// No samples - should have no issues
	issues := manager.ValidateSamples()
	if len(issues) != 0 {
		t.Errorf("Empty manager should have no validation issues, got %d", len(issues))
	}

	// Add valid samples
	manager.AddSample(Sample{ID: "Sample1", Lane: 1, Index1: "ATCGATCG"})
	manager.AddSample(Sample{ID: "Sample2", Lane: 1, Index1: "CGTAGATC"})

	issues = manager.ValidateSamples()
	if len(issues) != 0 {
		t.Errorf("Valid samples should have no issues, got %d: %v", len(issues), issues)
	}

	// Add duplicate ID
	manager.AddSample(Sample{ID: "Sample1", Lane: 2, Index1: "GCTAGCTA"})
	issues = manager.ValidateSamples()
	hasDuplicateID := false
	for _, issue := range issues {
		if containsString(issue, "Duplicate Sample ID") {
			hasDuplicateID = true
			break
		}
	}
	if !hasDuplicateID {
		t.Error("Should detect duplicate sample ID")
	}
}

func TestSampleManagerValidateSamplesEmptyID(t *testing.T) {
	manager := NewSampleManager()
	manager.AddSample(Sample{ID: "", Lane: 1, Index1: "ATCGATCG"})

	issues := manager.ValidateSamples()
	hasEmptyID := false
	for _, issue := range issues {
		if containsString(issue, "Empty Sample ID") {
			hasEmptyID = true
			break
		}
	}
	if !hasEmptyID {
		t.Error("Should detect empty sample ID")
	}
}

func TestSampleManagerValidateSamplesMissingIndex(t *testing.T) {
	manager := NewSampleManager()
	manager.AddSample(Sample{ID: "Sample1", Lane: 1}) // No Index1 or Index2

	issues := manager.ValidateSamples()
	hasMissingIndex := false
	for _, issue := range issues {
		if containsString(issue, "Missing index sequences") {
			hasMissingIndex = true
			break
		}
	}
	if !hasMissingIndex {
		t.Error("Should detect missing index sequences")
	}
}

func TestSampleManagerGetUniqueLanes(t *testing.T) {
	manager := NewSampleManager()

	// Empty manager
	lanes := manager.GetUniqueLanes()
	if len(lanes) != 0 {
		t.Errorf("Expected 0 lanes for empty manager, got %d", len(lanes))
	}

	// Add samples in multiple lanes
	manager.AddSample(Sample{ID: "Sample1", Lane: 1})
	manager.AddSample(Sample{ID: "Sample2", Lane: 1}) // Same lane
	manager.AddSample(Sample{ID: "Sample3", Lane: 2})
	manager.AddSample(Sample{ID: "Sample4", Lane: 3})

	lanes = manager.GetUniqueLanes()
	if len(lanes) != 3 {
		t.Errorf("Expected 3 unique lanes, got %d", len(lanes))
	}

	// Verify all expected lanes are present
	laneSet := make(map[int]bool)
	for _, lane := range lanes {
		laneSet[lane] = true
	}
	for _, expectedLane := range []int{1, 2, 3} {
		if !laneSet[expectedLane] {
			t.Errorf("Expected lane %d to be in unique lanes", expectedLane)
		}
	}
}

func TestSampleManagerGetSamplesByLane(t *testing.T) {
	manager := NewSampleManager()

	// Add samples in multiple lanes
	manager.AddSample(Sample{ID: "Sample1", Lane: 1})
	manager.AddSample(Sample{ID: "Sample2", Lane: 1})
	manager.AddSample(Sample{ID: "Sample3", Lane: 2})

	// Get lane 1 samples
	lane1Samples := manager.GetSamplesByLane(1)
	if len(lane1Samples) != 2 {
		t.Errorf("Expected 2 samples in lane 1, got %d", len(lane1Samples))
	}

	// Get lane 2 samples
	lane2Samples := manager.GetSamplesByLane(2)
	if len(lane2Samples) != 1 {
		t.Errorf("Expected 1 sample in lane 2, got %d", len(lane2Samples))
	}

	// Get non-existent lane
	lane3Samples := manager.GetSamplesByLane(3)
	if len(lane3Samples) != 0 {
		t.Errorf("Expected 0 samples in lane 3, got %d", len(lane3Samples))
	}
}

func TestSuggestIndexes(t *testing.T) {
	indexes := SuggestIndexes()

	if len(indexes) == 0 {
		t.Error("SuggestIndexes should return some indexes")
	}

	// All returned indexes should be valid
	for _, idx := range indexes {
		if !ValidateIndexSequence(idx) {
			t.Errorf("SuggestIndexes returned invalid index: %s", idx)
		}
	}
}

func TestValidateIndexSequence(t *testing.T) {
	testCases := []struct {
		sequence string
		valid    bool
	}{
		{"ATCGATCG", true},  // Valid 8-base sequence
		{"ATCGAT", true},    // Valid 6-base sequence
		{"ATCGATCGAT", true}, // Valid 10-base sequence
		{"", false},           // Empty
		{"ATCGATCX", false},  // Invalid character X
		{"atcgatcg", true},   // Lowercase (should be accepted)
		{"ATCG", false},      // Too short (4 bases)
	}

	for _, tc := range testCases {
		result := ValidateIndexSequence(tc.sequence)
		if result != tc.valid {
			t.Errorf("ValidateIndexSequence(%q) = %v, expected %v", tc.sequence, result, tc.valid)
		}
	}
}

// Helper function to check if string contains substring
func containsString(s, substr string) bool {
	return len(s) >= len(substr) && (s == substr || containsStringHelper(s, substr))
}

func containsStringHelper(s, substr string) bool {
	for i := 0; i <= len(s)-len(substr); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}