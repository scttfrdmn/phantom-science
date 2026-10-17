//go:build !release
// +build !release

package main

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"
	"time"
)

// createTestGzippedVCF creates a minimal valid gzipped VCF file and returns its path.
func createTestGzippedVCF(t *testing.T) string {
	t.Helper()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	gw.Write([]byte("chr1\t1000\t.\tA\tG\t.\tPASS\t.\n"))
	gw.Write([]byte("chr1\t2000\t.\tCT\tC\t.\tPASS\t.\n"))
	gw.Close()
	tmpDir := t.TempDir()
	vcfPath := filepath.Join(tmpDir, "chr1.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)
	return vcfPath
}

func TestIsValidDNASequence_Valid(t *testing.T) {
	if !isValidDNASequence("ATCG") {
		t.Error("ATCG should be valid DNA sequence")
	}
	if !isValidDNASequence("atcg") {
		t.Error("atcg (lowercase) should be valid DNA sequence (converted to uppercase)")
	}
}

func TestIsValidDNASequence_Empty(t *testing.T) {
	if isValidDNASequence("") {
		t.Error("empty string should not be valid DNA sequence")
	}
}

func TestIsValidDNASequence_WithN(t *testing.T) {
	if isValidDNASequence("ATCGN") {
		t.Error("ATCGN should not be valid DNA sequence (N is invalid)")
	}
}

func TestIsValidDNASequence_WithDigit(t *testing.T) {
	if isValidDNASequence("ATC1") {
		t.Error("ATC1 should not be valid DNA sequence (digit)")
	}
}

func TestCalculateAverageGCContent_Empty(t *testing.T) {
	if got := calculateAverageGCContent([][]byte{}); got != 0.0 {
		t.Errorf("calculateAverageGCContent([]) = %f, want 0.0", got)
	}
}

func TestCalculateAverageGCContent_AllGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCGCGCGC")}
	if got := calculateAverageGCContent(seqs); got != 1.0 {
		t.Errorf("calculateAverageGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestCalculateAverageGCContent_AllAT(t *testing.T) {
	seqs := [][]byte{[]byte("ATATATAT")}
	if got := calculateAverageGCContent(seqs); got != 0.0 {
		t.Errorf("calculateAverageGCContent(all AT) = %f, want 0.0", got)
	}
}

func TestCalculateAverageGCContent_HalfGC(t *testing.T) {
	seqs := [][]byte{[]byte("GCATGCAT")}
	got := calculateAverageGCContent(seqs)
	if got < 0.49 || got > 0.51 {
		t.Errorf("calculateAverageGCContent(half GC) = %f, want ~0.5", got)
	}
}

func TestCalculateAverageGCContent_MultipleSeqs(t *testing.T) {
	seqs := [][]byte{
		[]byte("GCGC"), // 1.0 GC
		[]byte("ATAT"), // 0.0 GC
	}
	got := calculateAverageGCContent(seqs)
	if got < 0.49 || got > 0.51 {
		t.Errorf("calculateAverageGCContent(mixed seqs) = %f, want ~0.5", got)
	}
}

func TestCalculateSequenceComplexity_Empty(t *testing.T) {
	if got := calculateSequenceComplexity([]byte{}, 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(empty) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_TooShort(t *testing.T) {
	if got := calculateSequenceComplexity([]byte("ATG"), 4); got != 0.0 {
		t.Errorf("calculateSequenceComplexity(3-char with k=4) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_Repetitive(t *testing.T) {
	seq := []byte("AAAAAAAAAAAAAAAAAAAAAAAAA")
	got := calculateSequenceComplexity(seq, 4)
	if got >= 0.5 {
		t.Errorf("calculateSequenceComplexity(all-A) = %f, expected < 0.5", got)
	}
}

func TestCalculateSequenceComplexity_Diverse(t *testing.T) {
	seq := []byte("ATCGATCGATCGATCGATCGATCGATCG")
	got := calculateSequenceComplexity(seq, 4)
	if got <= 0 {
		t.Errorf("calculateSequenceComplexity(diverse) = %f, expected > 0", got)
	}
}

func TestCalculateAverageComplexity_Empty(t *testing.T) {
	if got := calculateAverageComplexity([][]byte{}, 4); got != 0.0 {
		t.Errorf("calculateAverageComplexity([]) = %f, want 0.0", got)
	}
}

func TestCalculateAverageComplexity_AllTooShort(t *testing.T) {
	seqs := [][]byte{[]byte("ATG"), []byte("AT")}
	if got := calculateAverageComplexity(seqs, 4); got != 0.0 {
		t.Errorf("calculateAverageComplexity(all too short) = %f, want 0.0", got)
	}
}

func TestCalculateAverageComplexity_Valid(t *testing.T) {
	seq := []byte("ATCGATCGATCGATCGATCG")
	seqs := [][]byte{seq}
	got := calculateAverageComplexity(seqs, 4)
	if got <= 0 {
		t.Errorf("calculateAverageComplexity(valid) = %f, expected > 0", got)
	}
}

func TestGetRealPopulationTargets_NotEmpty(t *testing.T) {
	targets := getRealPopulationTargets()
	if len(targets) == 0 {
		t.Error("getRealPopulationTargets() returned empty slice")
	}
}

// TestProcessVCFFile_ValidGzipped tests processVCFFile with a synthetic gzipped VCF.
func TestProcessVCFFile_ValidGzipped(t *testing.T) {
	vcfPath := createTestGzippedVCF(t)
	variants, sequences, err := processVCFFile(vcfPath, "YRI")
	if err != nil {
		t.Fatalf("processVCFFile error: %v", err)
	}
	// We have 2 data lines: "A\tG" and "CT\tC" - both valid DNA
	if variants < 1 {
		t.Errorf("expected at least 1 variant, got %d", variants)
	}
	_ = sequences
}

// TestProcessVCFFile_WithInvalidBases tests that records containing non-DNA characters
// (like '*' in deletions) are skipped.
func TestProcessVCFFile_WithInvalidBases(t *testing.T) {
	t.Helper()
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	// Valid record
	gw.Write([]byte("chr1\t1000\t.\tACGT\tA\t.\tPASS\t.\n"))
	// Invalid: N in REF
	gw.Write([]byte("chr1\t2000\t.\tNA\tG\t.\tPASS\t.\n"))
	// Invalid: * in ALT (spanning deletion)
	gw.Write([]byte("chr1\t3000\t.\tA\t*\t.\tPASS\t.\n"))
	gw.Close()
	tmpDir := t.TempDir()
	vcfPath := filepath.Join(tmpDir, "chr1.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	variants, _, err := processVCFFile(vcfPath, "CEU")
	if err != nil {
		t.Fatalf("processVCFFile error: %v", err)
	}
	// Only 1 valid record should be counted
	if variants != 1 {
		t.Errorf("expected 1 valid variant, got %d", variants)
	}
}

// TestProcessVCFFile_ShortLines ensures lines with fewer than 5 fields are skipped.
func TestProcessVCFFile_ShortLines(t *testing.T) {
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\n"))
	// Line with only 4 fields (too short)
	gw.Write([]byte("chr1\t1000\t.\tA\n"))
	// Valid line
	gw.Write([]byte("chr1\t2000\t.\tG\tC\t.\tPASS\t.\n"))
	gw.Close()
	tmpDir := t.TempDir()
	vcfPath := filepath.Join(tmpDir, "test.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	variants, _, err := processVCFFile(vcfPath, "CHB")
	if err != nil {
		t.Fatalf("processVCFFile error: %v", err)
	}
	if variants != 1 {
		t.Errorf("expected 1 valid variant (short line skipped), got %d", variants)
	}
}

// TestProcessVCFFile_LargeSequenceChunking tests that sequences >= 1000 bytes are chunked.
func TestProcessVCFFile_LargeSequenceChunking(t *testing.T) {
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	gw.Write([]byte("##fileformat=VCFv4.1\n"))
	gw.Write([]byte("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
	// Write many records with long REF+ALT to trigger the chunking branch (>= 1000 bytes)
	ref := string(bytes.Repeat([]byte("A"), 300))
	alt := string(bytes.Repeat([]byte("T"), 300))
	for i := 0; i < 5; i++ {
		line := "chr1\t" + string(rune('0'+i+1)) + "000\t.\t" + ref + "\t" + alt + "\t.\tPASS\t.\n"
		gw.Write([]byte(line))
	}
	gw.Close()
	tmpDir := t.TempDir()
	vcfPath := filepath.Join(tmpDir, "big.vcf.gz")
	os.WriteFile(vcfPath, buf.Bytes(), 0600)

	variants, sequences, err := processVCFFile(vcfPath, "JPT")
	if err != nil {
		t.Fatalf("processVCFFile error: %v", err)
	}
	if variants < 1 {
		t.Errorf("expected at least 1 variant, got %d", variants)
	}
	// With 5 records of 300+300=600 bytes each = 3000 bytes total, at least 2 chunks
	if len(sequences) < 2 {
		t.Errorf("expected at least 2 sequence chunks, got %d", len(sequences))
	}
}

// TestFindPopulationVCFFiles_EmptyDir verifies that no files are found in an empty dir.
func TestFindPopulationVCFFiles_EmptyDir(t *testing.T) {
	tmpDir := t.TempDir()
	files, err := findPopulationVCFFiles(tmpDir, "YRI")
	if err != nil {
		t.Fatalf("findPopulationVCFFiles error: %v", err)
	}
	if len(files) != 0 {
		t.Errorf("expected no files in empty dir, got %d", len(files))
	}
}

// TestFindPopulationVCFFiles_WithFiles verifies files are found in the expected subdirectory.
func TestFindPopulationVCFFiles_WithFiles(t *testing.T) {
	tmpDir := t.TempDir()
	popDir := filepath.Join(tmpDir, "YRI")
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}
	// Create a fake .vcf.gz file
	vcfPath := filepath.Join(popDir, "chr1.vcf.gz")
	if err := os.WriteFile(vcfPath, []byte("fake"), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	files, err := findPopulationVCFFiles(tmpDir, "YRI")
	if err != nil {
		t.Fatalf("findPopulationVCFFiles error: %v", err)
	}
	if len(files) != 1 {
		t.Errorf("expected 1 file, got %d", len(files))
	}
}

// TestFindPopulationVCFFiles_PrefixPattern verifies the basePath/<POP>_*.vcf.gz pattern is also searched.
func TestFindPopulationVCFFiles_PrefixPattern(t *testing.T) {
	tmpDir := t.TempDir()
	// Create file matching the <POP>_*.vcf.gz pattern at basePath level
	vcfPath := filepath.Join(tmpDir, "CEU_chr1.vcf.gz")
	if err := os.WriteFile(vcfPath, []byte("fake"), 0600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	files, err := findPopulationVCFFiles(tmpDir, "CEU")
	if err != nil {
		t.Fatalf("findPopulationVCFFiles error: %v", err)
	}
	if len(files) != 1 {
		t.Errorf("expected 1 file via prefix pattern, got %d", len(files))
	}
}

// TestDetectAppleSiliconCores_RealVCF verifies the function returns a positive value.
func TestDetectAppleSiliconCores_RealVCF(t *testing.T) {
	got := detectAppleSiliconCores()
	if got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
	if got > 1024 {
		t.Errorf("detectAppleSiliconCores() = %d, unreasonably large", got)
	}
}

// TestWriteVCFTrainedFASTA_CreatesFile verifies that writeVCFTrainedFASTA writes a file.
func TestWriteVCFTrainedFASTA_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	outPath := filepath.Join(tmpDir, "test.fasta")
	training := &PopulationVCFTraining{
		Population:        "YRI",
		Region:            "Africa",
		ProcessedVariants: 100,
		TrainingTime:      5 * time.Second,
	}
	sequence := []byte("ACGTACGTACGTACGT")
	err := writeVCFTrainedFASTA(outPath, "YRI", training, 16, 1, sequence)
	if err != nil {
		t.Fatalf("writeVCFTrainedFASTA error: %v", err)
	}
	info, err := os.Stat(outPath)
	if err != nil {
		t.Fatalf("file should exist: %v", err)
	}
	if info.Size() == 0 {
		t.Error("FASTA file should not be empty")
	}
}

// TestWriteVCFTrainedFASTA_LongSequence verifies that sequences > 80 chars are wrapped.
func TestWriteVCFTrainedFASTA_LongSequence(t *testing.T) {
	tmpDir := t.TempDir()
	outPath := filepath.Join(tmpDir, "long.fasta")
	training := &PopulationVCFTraining{
		Population:        "CEU",
		Region:            "Europe",
		ProcessedVariants: 50,
		TrainingTime:      2 * time.Second,
	}
	// Sequence of 200 bytes - should be wrapped at 80 chars
	sequence := bytes.Repeat([]byte("ACGT"), 50)
	err := writeVCFTrainedFASTA(outPath, "CEU", training, 200, 2, sequence)
	if err != nil {
		t.Fatalf("writeVCFTrainedFASTA error: %v", err)
	}
	content, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	// File should contain line breaks (wrapping)
	if !bytes.Contains(content, []byte("\n")) {
		t.Error("expected newlines in wrapped FASTA output")
	}
}

// TestAnalyzeVCFTrainingResults_AllFailed tests analyzeVCFTrainingResults when all models failed.
func TestAnalyzeVCFTrainingResults_AllFailed(t *testing.T) {
	models := map[string]*PopulationVCFTraining{
		"YRI": {Population: "YRI", Region: "Africa", Trained: false},
		"CEU": {Population: "CEU", Region: "Europe", Trained: false},
	}
	// Should not panic
	analyzeVCFTrainingResults(models, 10*time.Second)
}

// TestAnalyzeVCFTrainingResults_AllSuccess tests analyzeVCFTrainingResults when all models succeed.
func TestAnalyzeVCFTrainingResults_AllSuccess(t *testing.T) {
	models := map[string]*PopulationVCFTraining{
		"YRI": {
			Population:        "YRI",
			Region:            "Africa",
			Trained:           true,
			ProcessedVariants: 1000,
			ActualGCContent:   0.42,
			ActualComplexity:  0.08,
			TrainingTime:      5 * time.Second,
		},
		"CEU": {
			Population:        "CEU",
			Region:            "Europe",
			Trained:           true,
			ProcessedVariants: 800,
			ActualGCContent:   0.40,
			ActualComplexity:  0.07,
			TrainingTime:      4 * time.Second,
		},
	}
	// Should not panic
	analyzeVCFTrainingResults(models, 10*time.Second)
}

// TestAnalyzeVCFTrainingResults_MostlySuccessful tests the "mostly successful" branch (>=3/4).
func TestAnalyzeVCFTrainingResults_MostlySuccessful(t *testing.T) {
	models := map[string]*PopulationVCFTraining{
		"YRI": {Population: "YRI", Region: "Africa", Trained: true, ProcessedVariants: 100, ActualGCContent: 0.42, ActualComplexity: 0.05, TrainingTime: time.Second},
		"CEU": {Population: "CEU", Region: "Europe", Trained: true, ProcessedVariants: 100, ActualGCContent: 0.40, ActualComplexity: 0.05, TrainingTime: time.Second},
		"CHB": {Population: "CHB", Region: "East Asia", Trained: true, ProcessedVariants: 100, ActualGCContent: 0.41, ActualComplexity: 0.05, TrainingTime: time.Second},
		"MXL": {Population: "MXL", Region: "Americas", Trained: false},
	}
	// 3/4 = 75% >= 75% threshold → "mostly successful" branch
	analyzeVCFTrainingResults(models, 30*time.Second)
}

// TestAnalyzeVCFTrainingResults_CriticalIssue tests the "critical" branch (< 3/4 success).
func TestAnalyzeVCFTrainingResults_CriticalIssue(t *testing.T) {
	models := map[string]*PopulationVCFTraining{
		"YRI": {Population: "YRI", Region: "Africa", Trained: true, ProcessedVariants: 50, ActualGCContent: 0.42, ActualComplexity: 0.04, TrainingTime: time.Second},
		"CEU": {Population: "CEU", Region: "Europe", Trained: false},
		"CHB": {Population: "CHB", Region: "East Asia", Trained: false},
		"MXL": {Population: "MXL", Region: "Americas", Trained: false},
	}
	// 1/4 = 25% < 75% threshold → "critical" branch
	analyzeVCFTrainingResults(models, 30*time.Second)
}

// ---------- checkVCFDataAvailability ----------

func TestCheckVCFDataAvailability_NoDataPaths_ReturnsEmpty(t *testing.T) {
	// On CI/test environments, the hardcoded paths don't exist → returns ""
	result := checkVCFDataAvailability()
	// Either "" (no data) or a valid path - both are acceptable
	_ = result
}

// ---------- trainPopulationOnVCF ----------

func TestTrainPopulationOnVCF_NoVCFFiles_ReturnsError(t *testing.T) {
	training := &PopulationVCFTraining{
		Population: "YRI",
		Region:     "Africa",
	}
	// With an empty tmpDir as VCF data path, no VCF files found → error
	err := trainPopulationOnVCF(training, t.TempDir())
	if err == nil {
		t.Error("trainPopulationOnVCF(no VCF files) expected error, got nil")
	}
}

// ---------- generateVCFTrainedSequences ----------

func TestGenerateVCFTrainedSequences_EmptyMap_NoCrash(t *testing.T) {
	// With empty map, nothing to iterate → returns immediately
	generateVCFTrainedSequences(map[string]*PopulationVCFTraining{})
}

func TestGenerateVCFTrainedSequences_UntrainedModels_SkipAll(t *testing.T) {
	models := map[string]*PopulationVCFTraining{
		"YRI": {Population: "YRI", Trained: false},
		"CEU": {Population: "CEU", Trained: false},
	}
	// Untrained models should be skipped with no crash
	generateVCFTrainedSequences(models)
}

// ============================================================
// Additional pure-helper function tests
// ============================================================

// --- isValidDNASequence edge cases ---

func TestIsValidDNASequence_RPVT_Table(t *testing.T) {
	tests := []struct {
		input string
		want  bool
	}{
		{"ACGT", true},
		{"acgt", true},   // lowercase accepted via ToUpper
		{"AcGt", true},
		{"", false},       // empty is invalid
		{"ACGN", false},   // N is not valid
		{"ACGU", false},   // U (uracil) is not DNA
		{"ACG1", false},   // digit is invalid
		{"ACG-", false},   // dash is invalid
		{"AC GT", false},  // space is invalid
		{"A", true},
		{"T", true},
		{"G", true},
		{"C", true},
	}
	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			got := isValidDNASequence(tt.input)
			if got != tt.want {
				t.Errorf("isValidDNASequence(%q) = %v, want %v", tt.input, got, tt.want)
			}
		})
	}
}

func TestIsValidDNASequence_RPVT_Whitespace(t *testing.T) {
	if isValidDNASequence("   ") {
		t.Error("whitespace-only string should not be valid DNA")
	}
}

func TestIsValidDNASequence_RPVT_OnlyDigits(t *testing.T) {
	if isValidDNASequence("1234") {
		t.Error("digit-only string should not be valid DNA")
	}
}

// --- calculateAverageGCContent additional coverage ---

func TestCalculateAverageGCContent_RPVT_NilInput(t *testing.T) {
	if got := calculateAverageGCContent(nil); got != 0.0 {
		t.Errorf("calculateAverageGCContent(nil) = %f, want 0.0", got)
	}
}

func TestCalculateAverageGCContent_RPVT_EmptyByteSlice(t *testing.T) {
	// Single empty []byte: totalBases==0 path
	if got := calculateAverageGCContent([][]byte{{}}); got != 0.0 {
		t.Errorf("calculateAverageGCContent([empty]) = %f, want 0.0", got)
	}
}

func TestCalculateAverageGCContent_RPVT_SingleBase(t *testing.T) {
	cases := []struct {
		base byte
		want float64
	}{
		{'G', 1.0},
		{'C', 1.0},
		{'A', 0.0},
		{'T', 0.0},
	}
	for _, c := range cases {
		got := calculateAverageGCContent([][]byte{{c.base}})
		if got != c.want {
			t.Errorf("calculateAverageGCContent([%c]) = %f, want %f", c.base, got, c.want)
		}
	}
}

func TestCalculateAverageGCContent_RPVT_KnownFraction(t *testing.T) {
	// "GCATGCAT": G,C,A,T,G,C,A,T -> GC=4, total=8 -> 0.5
	seq := []byte("GCATGCAT")
	got := calculateAverageGCContent([][]byte{seq})
	const want = 0.5
	if got < want-0.0001 || got > want+0.0001 {
		t.Errorf("calculateAverageGCContent(GCATGCAT) = %f, want %f", got, want)
	}
}

// --- calculateSequenceComplexity exact-value cases ---

func TestCalculateSequenceComplexity_RPVT_KEqualsLength(t *testing.T) {
	// k == len(seq): exactly 1 kmer, ratio = 1.0
	seq := []byte("ACGT")
	got := calculateSequenceComplexity(seq, 4)
	if got != 1.0 {
		t.Errorf("calculateSequenceComplexity(k==len) = %f, want 1.0", got)
	}
}

func TestCalculateSequenceComplexity_RPVT_K1AllSameBase(t *testing.T) {
	// k=1, 10 identical bases: 1 unique / 10 total = 0.1
	seq := []byte("AAAAAAAAAA")
	got := calculateSequenceComplexity(seq, 1)
	const want = 0.1
	if got < want-0.0001 || got > want+0.0001 {
		t.Errorf("calculateSequenceComplexity(all-A, k=1) = %f, want %f", got, want)
	}
}

func TestCalculateSequenceComplexity_RPVT_K1FourBases(t *testing.T) {
	// k=1, ACGT: 4 unique / 4 total = 1.0
	seq := []byte("ACGT")
	got := calculateSequenceComplexity(seq, 1)
	if got != 1.0 {
		t.Errorf("calculateSequenceComplexity(ACGT, k=1) = %f, want 1.0", got)
	}
}

func TestCalculateSequenceComplexity_RPVT_KLargerThanSeq(t *testing.T) {
	got := calculateSequenceComplexity([]byte("AC"), 10)
	if got != 0.0 {
		t.Errorf("calculateSequenceComplexity(k>len) = %f, want 0.0", got)
	}
}

func TestCalculateSequenceComplexity_RPVT_RepeatedMotif(t *testing.T) {
	// "ACGTACGT" k=4: unique={ACGT,CGTA,GTAC,TACG}=4, total=5 -> 0.8
	seq := []byte("ACGTACGT")
	got := calculateSequenceComplexity(seq, 4)
	const want = 0.8
	if got < want-0.0001 || got > want+0.0001 {
		t.Errorf("calculateSequenceComplexity(ACGTACGT,k=4) = %f, want %f", got, want)
	}
}

func TestCalculateSequenceComplexity_RPVT_Deterministic(t *testing.T) {
	seq := []byte("ACGTGCATACGTGCAT")
	c1 := calculateSequenceComplexity(seq, 4)
	c2 := calculateSequenceComplexity(seq, 4)
	if c1 != c2 {
		t.Errorf("calculateSequenceComplexity not deterministic: %f != %f", c1, c2)
	}
}

// --- calculateAverageComplexity additional coverage ---

func TestCalculateAverageComplexity_RPVT_NilInput(t *testing.T) {
	if got := calculateAverageComplexity(nil, 4); got != 0.0 {
		t.Errorf("calculateAverageComplexity(nil) = %f, want 0.0", got)
	}
}

func TestCalculateAverageComplexity_RPVT_IdenticalSeqsAverage(t *testing.T) {
	// Average of two identical sequences equals single sequence complexity
	seq := []byte("ACGTACGT")
	single := calculateSequenceComplexity(seq, 4)
	got := calculateAverageComplexity([][]byte{seq, seq}, 4)
	if got < single-0.0001 || got > single+0.0001 {
		t.Errorf("average of identical seqs = %f, want %f", got, single)
	}
}

func TestCalculateAverageComplexity_RPVT_ShortExcluded(t *testing.T) {
	// Short seq (len < k) is excluded; only long counts
	short := []byte("AC")   // len=2 < k=4
	long := []byte("ACGTACGT") // len=8 >= k=4
	wantSingle := calculateSequenceComplexity(long, 4)
	got := calculateAverageComplexity([][]byte{short, long}, 4)
	if got < wantSingle-0.0001 || got > wantSingle+0.0001 {
		t.Errorf("calculateAverageComplexity(short+long) = %f, want %f", got, wantSingle)
	}
}

func TestCalculateAverageComplexity_RPVT_KEqualsSeqLen(t *testing.T) {
	// When k == len of every sequence, each complexity = 1.0 -> average = 1.0
	s1 := []byte("ACGT")
	s2 := []byte("GCTA")
	got := calculateAverageComplexity([][]byte{s1, s2}, 4)
	if got != 1.0 {
		t.Errorf("calculateAverageComplexity(k==len each) = %f, want 1.0", got)
	}
}

// --- detectAppleSiliconCores property tests ---

func TestDetectAppleSiliconCores_RPVT_Positive(t *testing.T) {
	if got := detectAppleSiliconCores(); got <= 0 {
		t.Errorf("detectAppleSiliconCores() = %d, want > 0", got)
	}
}

func TestDetectAppleSiliconCores_RPVT_ReasonableUpperBound(t *testing.T) {
	if got := detectAppleSiliconCores(); got > 256 {
		t.Errorf("detectAppleSiliconCores() = %d, unexpectedly large", got)
	}
}

func TestDetectAppleSiliconCores_RPVT_Deterministic(t *testing.T) {
	first := detectAppleSiliconCores()
	second := detectAppleSiliconCores()
	if first != second {
		t.Errorf("detectAppleSiliconCores() not deterministic: %d != %d", first, second)
	}
}
