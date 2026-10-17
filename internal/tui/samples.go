package tui

import (
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/security"
)

// SampleManager handles sample sheet operations
type SampleManager struct {
	samples []Sample
}

// NewSampleManager creates a new sample manager
func NewSampleManager() *SampleManager {
	return &SampleManager{
		samples: []Sample{},
	}
}

// AddSample adds a new sample to the sheet
func (sm *SampleManager) AddSample(sample Sample) {
	sm.samples = append(sm.samples, sample)
}

// RemoveSample removes a sample by index
func (sm *SampleManager) RemoveSample(index int) {
	if index >= 0 && index < len(sm.samples) {
		sm.samples = append(sm.samples[:index], sm.samples[index+1:]...)
	}
}

// UpdateSample updates a sample at the given index
func (sm *SampleManager) UpdateSample(index int, sample Sample) {
	if index >= 0 && index < len(sm.samples) {
		sm.samples[index] = sample
	}
}

// GetSamples returns all samples
func (sm *SampleManager) GetSamples() []Sample {
	return sm.samples
}

// ClearAll removes all samples
func (sm *SampleManager) ClearAll() {
	sm.samples = []Sample{}
}

// ImportFromCSV imports samples from a CSV file
func (sm *SampleManager) ImportFromCSV(filename string) error {
	// Validate file path for security
	if err := security.ValidateFilePath(filename, ""); err != nil {
		return fmt.Errorf("invalid CSV file path: %w", err)
	}

	file, err := os.Open(filename) // #nosec G304 - path validated by security.ValidateFilePath
	if err != nil {
		return fmt.Errorf("failed to open CSV file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		return fmt.Errorf("failed to read CSV: %w", err)
	}

	if len(records) == 0 {
		return fmt.Errorf("CSV file is empty")
	}

	// Skip header row if it exists
	startRow := 0
	if len(records) > 0 && strings.ToLower(records[0][0]) == "sample_id" {
		startRow = 1
	}

	// Clear existing samples
	sm.samples = []Sample{}

	// Parse each row
	for i := startRow; i < len(records); i++ {
		record := records[i]
		if len(record) < 6 {
			continue // Skip malformed rows
		}

		lane, err := strconv.Atoi(record[1])
		if err != nil {
			lane = 1 // Default lane
		}

		sample := Sample{
			ID:          record[0],
			Lane:        lane,
			Index1:      record[2],
			Index2:      record[3],
			Project:     record[4],
			Description: record[5],
		}

		sm.samples = append(sm.samples, sample)
	}

	return nil
}

// ExportToCSV exports samples to a CSV file
func (sm *SampleManager) ExportToCSV(filename string) error {
	// Validate file path for security
	if err := security.ValidateFilePath(filename, ""); err != nil {
		return fmt.Errorf("invalid CSV file path: %w", err)
	}

	file, err := security.SecureFileCreate(filename, "")
	if err != nil {
		return fmt.Errorf("failed to create CSV file: %w", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write header
	header := []string{"Sample_ID", "Lane", "Index_1", "Index_2", "Project", "Description"}
	if err := writer.Write(header); err != nil {
		return fmt.Errorf("failed to write CSV header: %w", err)
	}

	// Write samples
	for _, sample := range sm.samples {
		record := []string{
			sample.ID,
			strconv.Itoa(sample.Lane),
			sample.Index1,
			sample.Index2,
			sample.Project,
			sample.Description,
		}
		if err := writer.Write(record); err != nil {
			return fmt.Errorf("failed to write CSV record: %w", err)
		}
	}

	return nil
}

// ValidateSamples checks for common issues in the sample sheet
func (sm *SampleManager) ValidateSamples() []string {
	var issues []string

	// Check for duplicate sample IDs
	idMap := make(map[string]int)
	for i, sample := range sm.samples {
		if prevIndex, exists := idMap[sample.ID]; exists {
			issues = append(issues, fmt.Sprintf("Duplicate Sample ID '%s' at rows %d and %d", sample.ID, prevIndex+1, i+1))
		}
		idMap[sample.ID] = i
	}

	// Check for duplicate index combinations within the same lane
	indexMap := make(map[string][]int)
	for i, sample := range sm.samples {
		key := fmt.Sprintf("L%d-%s-%s", sample.Lane, sample.Index1, sample.Index2)
		indexMap[key] = append(indexMap[key], i)
	}

	for key, indices := range indexMap {
		if len(indices) > 1 {
			rowNums := make([]string, len(indices))
			for i, idx := range indices {
				rowNums[i] = strconv.Itoa(idx + 1)
			}
			issues = append(issues, fmt.Sprintf("Duplicate index combination '%s' at rows: %s", key, strings.Join(rowNums, ", ")))
		}
	}

	// Check for empty required fields
	for i, sample := range sm.samples {
		if sample.ID == "" {
			issues = append(issues, fmt.Sprintf("Empty Sample ID at row %d", i+1))
		}
		if sample.Index1 == "" && sample.Index2 == "" {
			issues = append(issues, fmt.Sprintf("Missing index sequences at row %d", i+1))
		}
	}

	return issues
}

// GetSampleCount returns the number of samples
func (sm *SampleManager) GetSampleCount() int {
	return len(sm.samples)
}

// GetUniqueLanes returns a list of unique lanes used
func (sm *SampleManager) GetUniqueLanes() []int {
	laneMap := make(map[int]bool)
	for _, sample := range sm.samples {
		laneMap[sample.Lane] = true
	}

	var lanes []int
	for lane := range laneMap {
		lanes = append(lanes, lane)
	}
	return lanes
}

// GetSamplesByLane returns samples for a specific lane
func (sm *SampleManager) GetSamplesByLane(lane int) []Sample {
	var result []Sample
	for _, sample := range sm.samples {
		if sample.Lane == lane {
			result = append(result, sample)
		}
	}
	return result
}

// SuggestIndexes returns commonly used index sequences
func SuggestIndexes() []string {
	return []string{
		"ATCGATCG", "CGTAGATC", "GCTAGCTA", "TACGACGT",
		"AGTCAACA", "CATGCCTA", "GTAGAGAG", "TCCTGAGC",
		"CGATGTTT", "TGACCACT", "AACGTGAT", "GCGGCTAG",
		"ATGCCTAA", "GAGTGGAT", "ACAGTGGT", "CGCTCATT",
	}
}

// ValidateIndexSequence checks if an index sequence is valid
func ValidateIndexSequence(sequence string) bool {
	if len(sequence) == 0 {
		return false
	}
	
	// Check for valid DNA bases
	for _, char := range strings.ToUpper(sequence) {
		if char != 'A' && char != 'T' && char != 'G' && char != 'C' {
			return false
		}
	}
	
	// Common index lengths are 6, 8, or 10 bases
	validLengths := []int{6, 8, 10}
	for _, length := range validLengths {
		if len(sequence) == length {
			return true
		}
	}
	
	return false
}