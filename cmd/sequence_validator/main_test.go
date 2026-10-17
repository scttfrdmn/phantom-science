//go:build !release
// +build !release

package main

import (
	"encoding/json"
	"os"
	"testing"
)

func TestParseFilename_Valid(t *testing.T) {
	pop, length, seqNum, err := parseFilename("YRI_L1000_S01.fasta")
	if err != nil {
		t.Fatalf("Expected no error, got: %v", err)
	}
	if pop != "YRI" {
		t.Errorf("Expected population 'YRI', got '%s'", pop)
	}
	if length != 1000 {
		t.Errorf("Expected length 1000, got %d", length)
	}
	if seqNum != 1 {
		t.Errorf("Expected seqNum 1, got %d", seqNum)
	}
}

func TestParseFilename_Invalid(t *testing.T) {
	_, _, _, err := parseFilename("invalid_name.fasta")
	if err == nil {
		t.Error("Expected error for invalid filename format")
	}
}

func TestParseFilename_InvalidLength(t *testing.T) {
	_, _, _, err := parseFilename("YRI_Lxxx_S01.fasta")
	if err == nil {
		t.Error("Expected error for non-numeric length")
	}
}

func TestParseFilename_InvalidSeqNum(t *testing.T) {
	_, _, _, err := parseFilename("YRI_L1000_Sxxx.fasta")
	if err == nil {
		t.Error("Expected error for non-numeric sequence number")
	}
}

func TestParseHeaderMetadata_WithGC(t *testing.T) {
	header := ">YRI_L1000_S01 population=YRI region=Coding model_target_gc=0.496"
	gc, region := parseHeaderMetadata(header)
	if gc < 0.495 || gc > 0.497 {
		t.Errorf("Expected GC ~0.496, got %f", gc)
	}
	if region != "Coding" {
		t.Errorf("Expected region 'Coding', got '%s'", region)
	}
}

func TestParseHeaderMetadata_NoGC(t *testing.T) {
	header := ">YRI_L1000_S01 population=YRI"
	gc, region := parseHeaderMetadata(header)
	if gc != 0.0 {
		t.Errorf("Expected GC 0.0 when not present, got %f", gc)
	}
	if region != "" {
		t.Errorf("Expected empty region when not present, got '%s'", region)
	}
}

func TestCountValidBases_AllValid(t *testing.T) {
	valid, invalid := countValidBases("ATCGATCG")
	if valid != 8 {
		t.Errorf("Expected 8 valid bases, got %d", valid)
	}
	if invalid != 0 {
		t.Errorf("Expected 0 invalid bases, got %d", invalid)
	}
}

func TestCountValidBases_WithInvalid(t *testing.T) {
	valid, invalid := countValidBases("ATCGNNXX")
	if valid != 4 {
		t.Errorf("Expected 4 valid bases, got %d", valid)
	}
	if invalid != 4 {
		t.Errorf("Expected 4 invalid bases, got %d", invalid)
	}
}

func TestCountValidBases_LowerCase(t *testing.T) {
	valid, invalid := countValidBases("atcgatcg")
	if valid != 8 {
		t.Errorf("Expected 8 valid bases (lowercase), got %d", valid)
	}
	if invalid != 0 {
		t.Errorf("Expected 0 invalid bases, got %d", invalid)
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	gc := calculateGCContent("GCGCGCGC")
	if gc < 0.99 || gc > 1.01 {
		t.Errorf("Expected GC ~1.0 for all-GC, got %f", gc)
	}
}

func TestCalculateGCContent_NoGC(t *testing.T) {
	gc := calculateGCContent("ATATATAT")
	if gc > 0.01 {
		t.Errorf("Expected GC ~0.0 for no-GC, got %f", gc)
	}
}

func TestCalculateGCContent_FiftyPercent(t *testing.T) {
	gc := calculateGCContent("GCATGCAT")
	if gc < 0.49 || gc > 0.51 {
		t.Errorf("Expected GC ~0.5, got %f", gc)
	}
}

func TestCalculateGCContent_Empty(t *testing.T) {
	gc := calculateGCContent("")
	if gc != 0.0 {
		t.Errorf("Expected GC 0.0 for empty sequence, got %f", gc)
	}
}

func TestCalculateComplexity_HighComplexity(t *testing.T) {
	// Diverse sequence should have high complexity
	seq := "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
	complexity := calculateComplexity(seq, 4)
	if complexity <= 0 {
		t.Errorf("Expected positive complexity for diverse sequence, got %f", complexity)
	}
}

func TestCalculateComplexity_LowComplexity(t *testing.T) {
	// Repetitive sequence should have low complexity
	seq := "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	complexity := calculateComplexity(seq, 4)
	if complexity >= 0.5 {
		t.Errorf("Expected low complexity for repetitive sequence, got %f", complexity)
	}
}

func TestCalculateComplexity_TooShort(t *testing.T) {
	complexity := calculateComplexity("AT", 8)
	if complexity != 0.0 {
		t.Errorf("Expected 0.0 for sequence shorter than k, got %f", complexity)
	}
}

func TestFindMaxHomopolymer_Single(t *testing.T) {
	max := findMaxHomopolymer("ATCG")
	if max != 1 {
		t.Errorf("Expected max 1 for no homopolymers, got %d", max)
	}
}

func TestFindMaxHomopolymer_Run(t *testing.T) {
	max := findMaxHomopolymer("ATAAAATCG")
	if max != 4 {
		t.Errorf("Expected max 4 for AAAA, got %d", max)
	}
}

func TestFindMaxHomopolymer_AtEnd(t *testing.T) {
	max := findMaxHomopolymer("ATCGGGGG")
	if max != 5 {
		t.Errorf("Expected max 5 for GGGGG at end, got %d", max)
	}
}

func TestFindMaxHomopolymer_Empty(t *testing.T) {
	max := findMaxHomopolymer("")
	if max != 0 {
		t.Errorf("Expected max 0 for empty sequence, got %d", max)
	}
}

func TestHasStopCodons_Present(t *testing.T) {
	if !hasStopCodons("ATCGTAA") {
		t.Error("Expected stop codon TAA to be detected")
	}
	if !hasStopCodons("ATCGTAG") {
		t.Error("Expected stop codon TAG to be detected")
	}
	if !hasStopCodons("ATCGTGA") {
		t.Error("Expected stop codon TGA to be detected")
	}
}

func TestHasStopCodons_Absent(t *testing.T) {
	if hasStopCodons("ATCGATCG") {
		t.Error("Should not detect stop codon in ATCGATCG")
	}
}

func TestHasStopCodons_CaseInsensitive(t *testing.T) {
	if !hasStopCodons("atcgtaa") {
		t.Error("Stop codon detection should be case insensitive")
	}
}

func TestClassifyComplexity(t *testing.T) {
	tests := []struct {
		complexity float64
		expected   string
	}{
		{0.9, "Very High"},
		{0.7, "High"},
		{0.5, "Medium"},
		{0.3, "Low"},
		{0.1, "Very Low"},
	}

	for _, tt := range tests {
		result := classifyComplexity(tt.complexity)
		if result != tt.expected {
			t.Errorf("classifyComplexity(%f) = %q, expected %q", tt.complexity, result, tt.expected)
		}
	}
}

func TestAbsFunction(t *testing.T) {
	if abs(5.0) != 5.0 {
		t.Errorf("abs(5.0) should be 5.0")
	}
	if abs(-5.0) != 5.0 {
		t.Errorf("abs(-5.0) should be 5.0")
	}
	if abs(0.0) != 0.0 {
		t.Errorf("abs(0.0) should be 0.0")
	}
}

// ---- Tests for file-writing and aggregation functions ----

func TestSaveValidationReport_CreatesFile(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-save-report-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	report := &ValidationReport{
		TotalSequences: 2,
		ValidSequences: 1,
		OverallScore:   0.75,
		ValidationCriteria: "test",
		PopulationResults:  []ValidationResult{},
		PopulationSummary:  []PopulationSummary{},
		RegionalAnalysis:   map[string]float64{"Africa": 0.8},
	}

	outFile := tmpDir + "/report.json"
	if err := saveValidationReport(report, outFile); err != nil {
		t.Fatalf("saveValidationReport returned error: %v", err)
	}

	info, err := os.Stat(outFile)
	if err != nil {
		t.Fatalf("output file not created: %v", err)
	}
	if info.Size() == 0 {
		t.Error("output file is empty")
	}
}

func TestSaveValidationReport_ValidJSON(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-save-json-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	report := &ValidationReport{
		TotalSequences: 5,
		ValidSequences: 4,
		OverallScore:   0.9,
		ValidationCriteria: "k-mer complexity",
		PopulationResults:  []ValidationResult{},
		PopulationSummary:  []PopulationSummary{},
		RegionalAnalysis:   map[string]float64{},
	}

	outFile := tmpDir + "/report.json"
	if err := saveValidationReport(report, outFile); err != nil {
		t.Fatalf("saveValidationReport: %v", err)
	}

	data, err := os.ReadFile(outFile)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}

	// Parse back and verify fields
	var loaded ValidationReport
	if err := json.Unmarshal(data, &loaded); err != nil {
		t.Fatalf("json.Unmarshal: %v", err)
	}
	if loaded.TotalSequences != 5 {
		t.Errorf("TotalSequences = %d, want 5", loaded.TotalSequences)
	}
	if loaded.ValidSequences != 4 {
		t.Errorf("ValidSequences = %d, want 4", loaded.ValidSequences)
	}
	if loaded.OverallScore < 0.89 || loaded.OverallScore > 0.91 {
		t.Errorf("OverallScore = %f, want ~0.9", loaded.OverallScore)
	}
}

func TestSaveValidationReport_InvalidPath(t *testing.T) {
	report := &ValidationReport{}
	err := saveValidationReport(report, "/nonexistent/dir/report.json")
	if err == nil {
		t.Error("expected error writing to nonexistent directory, got nil")
	}
}

func TestCalculateValidationScore_PerfectSequence(t *testing.T) {
	result := &ValidationResult{
		Complexity:   0.8,
		GCDeviation:  0.0,
		InvalidBases: 0,
		Issues:       []string{},
	}
	score := calculateValidationScore(result)
	if score < 0.0 || score > 1.0 {
		t.Errorf("score %f out of [0,1] range", score)
	}
	// High complexity, no deviation, no issues -> should be near top
	if score < 0.9 {
		t.Errorf("expected high score for perfect sequence, got %f", score)
	}
}

func TestCalculateValidationScore_LowComplexity(t *testing.T) {
	result := &ValidationResult{
		Complexity:   0.1,
		GCDeviation:  0.0,
		InvalidBases: 0,
		Issues:       []string{},
	}
	score := calculateValidationScore(result)
	if score < 0.0 || score > 1.0 {
		t.Errorf("score %f out of [0,1] range", score)
	}
	// Low complexity penalizes
	if score >= 0.9 {
		t.Errorf("expected lower score for low-complexity sequence, got %f", score)
	}
}

func TestCalculateValidationScore_ManyIssues(t *testing.T) {
	result := &ValidationResult{
		Complexity:   0.7,
		GCDeviation:  0.0,
		InvalidBases: 0,
		Issues:       []string{"issue1", "issue2", "issue3", "issue4", "issue5"},
	}
	score := calculateValidationScore(result)
	if score < 0.0 || score > 1.0 {
		t.Errorf("score %f out of [0,1] range", score)
	}
}

func TestCalculateValidationScore_LargeGCDeviation(t *testing.T) {
	result := &ValidationResult{
		Complexity:   0.7,
		GCDeviation:  0.5, // large deviation
		InvalidBases: 0,
		Issues:       []string{},
	}
	score := calculateValidationScore(result)
	if score < 0.0 || score > 1.0 {
		t.Errorf("score %f out of [0,1] range", score)
	}
	// Large GC deviation should produce a lower score
	perfectResult := &ValidationResult{
		Complexity:   0.7,
		GCDeviation:  0.0,
		InvalidBases: 0,
		Issues:       []string{},
	}
	perfectScore := calculateValidationScore(perfectResult)
	if score >= perfectScore {
		t.Errorf("large GC deviation score (%f) should be lower than zero-deviation score (%f)", score, perfectScore)
	}
}

func TestCalculateValidationScore_InvalidBases(t *testing.T) {
	result := &ValidationResult{
		Complexity:   0.7,
		GCDeviation:  0.0,
		InvalidBases: 5,
		Issues:       []string{},
	}
	score := calculateValidationScore(result)
	if score < 0.0 || score > 1.0 {
		t.Errorf("score %f out of [0,1] range", score)
	}
}

func TestGenerateValidationReport_Empty(t *testing.T) {
	report := generateValidationReport([]ValidationResult{})
	if report == nil {
		t.Fatal("expected non-nil report")
	}
	if report.TotalSequences != 0 {
		t.Errorf("TotalSequences = %d, want 0", report.TotalSequences)
	}
	if report.ValidSequences != 0 {
		t.Errorf("ValidSequences = %d, want 0", report.ValidSequences)
	}
	if report.OverallScore != 0.0 {
		t.Errorf("OverallScore = %f, want 0.0 for empty input", report.OverallScore)
	}
}

func TestGenerateValidationReport_SingleValid(t *testing.T) {
	results := []ValidationResult{
		{
			Population:      "YRI",
			Region:          "Africa",
			IsValid:         true,
			ValidationScore: 0.9,
			GCContent:       0.42,
			Complexity:      0.75,
			Issues:          []string{},
		},
	}
	report := generateValidationReport(results)
	if report.TotalSequences != 1 {
		t.Errorf("TotalSequences = %d, want 1", report.TotalSequences)
	}
	if report.ValidSequences != 1 {
		t.Errorf("ValidSequences = %d, want 1", report.ValidSequences)
	}
	if report.OverallScore < 0.89 || report.OverallScore > 0.91 {
		t.Errorf("OverallScore = %f, want ~0.9", report.OverallScore)
	}
	if len(report.PopulationSummary) != 1 {
		t.Errorf("PopulationSummary len = %d, want 1", len(report.PopulationSummary))
	}
	if _, exists := report.RegionalAnalysis["Africa"]; !exists {
		t.Error("RegionalAnalysis missing 'Africa' key")
	}
}

func TestGenerateValidationReport_MixedValidity(t *testing.T) {
	results := []ValidationResult{
		{Population: "YRI", Region: "Africa", IsValid: true, ValidationScore: 0.9, Issues: []string{}},
		{Population: "CEU", Region: "Europe", IsValid: false, ValidationScore: 0.4, Issues: []string{"low GC"}},
		{Population: "YRI", Region: "Africa", IsValid: true, ValidationScore: 0.8, Issues: []string{}},
	}
	report := generateValidationReport(results)
	if report.TotalSequences != 3 {
		t.Errorf("TotalSequences = %d, want 3", report.TotalSequences)
	}
	if report.ValidSequences != 2 {
		t.Errorf("ValidSequences = %d, want 2", report.ValidSequences)
	}
	// Should have summaries for YRI and CEU
	if len(report.PopulationSummary) != 2 {
		t.Errorf("PopulationSummary len = %d, want 2", len(report.PopulationSummary))
	}
}

func TestGenerateValidationReport_SortedByScore(t *testing.T) {
	results := []ValidationResult{
		{Population: "LOW", Region: "Region1", IsValid: false, ValidationScore: 0.3, Issues: []string{}},
		{Population: "HIGH", Region: "Region2", IsValid: true, ValidationScore: 0.9, Issues: []string{}},
	}
	report := generateValidationReport(results)
	if len(report.PopulationSummary) < 2 {
		t.Skip("not enough populations to test sort order")
	}
	// First summary should have higher avg score than second
	if report.PopulationSummary[0].AvgScore < report.PopulationSummary[1].AvgScore {
		t.Errorf("PopulationSummary not sorted descending: [0]=%f < [1]=%f",
			report.PopulationSummary[0].AvgScore, report.PopulationSummary[1].AvgScore)
	}
}

func TestValidateAllSequences_EmptyDir(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-validate-empty-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	results, err := validateAllSequences(tmpDir)
	if err != nil {
		t.Fatalf("unexpected error for empty dir: %v", err)
	}
	if len(results) != 0 {
		t.Errorf("expected 0 results for empty dir, got %d", len(results))
	}
}

func TestValidateAllSequences_WithFASTAFiles(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-validate-fasta-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	// Write a valid FASTA file with the expected name format
	fastaContent := ">YRI_L50_S01 population=YRI region=Africa model_target_gc=0.42\n" +
		"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC\n"
	fastaPath := tmpDir + "/YRI_L50_S01.fasta"
	if err := os.WriteFile(fastaPath, []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	results, err := validateAllSequences(tmpDir)
	if err != nil {
		t.Fatalf("validateAllSequences: %v", err)
	}
	if len(results) != 1 {
		t.Fatalf("expected 1 result, got %d", len(results))
	}
	if results[0].Population != "YRI" {
		t.Errorf("Population = %q, want YRI", results[0].Population)
	}
	if results[0].Region != "Africa" {
		t.Errorf("Region = %q, want Africa", results[0].Region)
	}
}

func TestValidateAllSequences_SkipsNonFASTA(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-validate-skip-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	// Write a non-FASTA file (should be skipped)
	if err := os.WriteFile(tmpDir+"/YRI_L50_S01.txt", []byte("not a fasta"), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	// Write a FASTA file with bad filename (should be skipped with warning)
	if err := os.WriteFile(tmpDir+"/bad_name.fasta", []byte(">header\nATCG\n"), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	results, err := validateAllSequences(tmpDir)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// bad_name.fasta fails parseFilename (3-part check) so is skipped
	if len(results) != 0 {
		t.Errorf("expected 0 results (non-FASTA and unparseable skipped), got %d", len(results))
	}
}

func TestValidateAllSequences_InvalidFASTAContent(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-validate-invalid-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	// FASTA with valid name but empty content → validateSequenceFile returns error → skip with warning
	if err := os.WriteFile(tmpDir+"/YRI_L50_S01.fasta", []byte(""), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	results, err := validateAllSequences(tmpDir)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// Empty file fails validateSequenceFile so is skipped
	if len(results) != 0 {
		t.Errorf("expected 0 results for invalid FASTA content, got %d", len(results))
	}
}

func TestValidateAllSequences_RoundTrip(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "test-validate-roundtrip-*")
	if err != nil {
		t.Fatalf("MkdirTemp: %v", err)
	}
	defer os.RemoveAll(tmpDir)

	// Write two FASTA files for two different populations
	files := []struct {
		name    string
		content string
	}{
		{
			"CEU_L100_S01.fasta",
			">CEU_L100_S01 population=CEU region=Europe model_target_gc=0.41\n" +
				"GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT\n",
		},
		{
			"YRI_L100_S02.fasta",
			">YRI_L100_S02 population=YRI region=Africa model_target_gc=0.43\n" +
				"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n",
		},
	}
	for _, f := range files {
		if err := os.WriteFile(tmpDir+"/"+f.name, []byte(f.content), 0600); err != nil {
			t.Fatalf("WriteFile %s: %v", f.name, err)
		}
	}

	results, err := validateAllSequences(tmpDir)
	if err != nil {
		t.Fatalf("validateAllSequences: %v", err)
	}
	if len(results) != 2 {
		t.Fatalf("expected 2 results, got %d", len(results))
	}

	// Now generate a report and save it
	report := generateValidationReport(results)
	outFile := tmpDir + "/out.json"
	if err := saveValidationReport(report, outFile); err != nil {
		t.Fatalf("saveValidationReport: %v", err)
	}

	data, err := os.ReadFile(outFile)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	var loaded ValidationReport
	if err := json.Unmarshal(data, &loaded); err != nil {
		t.Fatalf("json.Unmarshal: %v", err)
	}
	if loaded.TotalSequences != 2 {
		t.Errorf("loaded TotalSequences = %d, want 2", loaded.TotalSequences)
	}
}

// ---------- calculateValidationScore clamp-to-zero branch ----------

func TestCalculateValidationScore_ClampedToZero(t *testing.T) {
	// Many issues + low complexity + large GC deviation + invalid bases → score < 0 → clamped to 0
	result := &ValidationResult{
		Complexity:   0.1,  // < 0.3 → -0.3
		GCDeviation:  0.5,  // gcPenalty = 1.0
		InvalidBases: 1,    // -0.5
		Issues:       []string{"i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8"}, // 8 × -0.1 = -0.8
	}
	score := calculateValidationScore(result)
	if score != 0.0 {
		t.Errorf("expected score clamped to 0.0, got %f", score)
	}
}

// ---------- displayValidationResults ----------

func TestDisplayValidationResults_EmptyReport(t *testing.T) {
	report := &ValidationReport{
		TotalSequences: 10,
		ValidSequences: 8,
		OverallScore:   0.8,
		PopulationSummary: []PopulationSummary{},
		RegionalAnalysis:  map[string]float64{},
	}
	// Should not panic
	displayValidationResults(report)
}

func TestDisplayValidationResults_WithIssuesAndRegions(t *testing.T) {
	report := &ValidationReport{
		TotalSequences: 20,
		ValidSequences: 18,
		OverallScore:   0.9,
		PopulationSummary: []PopulationSummary{
			{
				Population:     "YRI",
				Region:         "Africa",
				TotalSequences: 10,
				ValidSequences: 9,
				AvgGCContent:   0.42,
				AvgComplexity:  0.75,
				AvgScore:       0.88,
				Issues:         map[string]int{},
			},
			{
				Population:     "CEU",
				Region:         "Europe",
				TotalSequences: 10,
				ValidSequences: 9,
				AvgGCContent:   0.41,
				AvgComplexity:  0.72,
				AvgScore:       0.85,
				// More than 2 issues covers the issueCount >= 2 break
				Issues: map[string]int{"low_complexity": 1, "high_gc": 2, "extra_issue": 3},
			},
		},
		RegionalAnalysis: map[string]float64{
			"Africa": 0.88,
			"Europe": 0.85,
		},
	}
	// Should not panic; covers PopulationSummary loop, Issues display, regional analysis
	displayValidationResults(report)
}
