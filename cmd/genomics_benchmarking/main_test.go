//go:build !release
// +build !release

package main

import (
	"fmt"
	"os"
	"path/filepath"
	"testing"
)

func TestSameRegion_SameAfrica(t *testing.T) {
	gb := &GenomicsBenchmark{}
	if !gb.sameRegion("YRI", "LWK") {
		t.Error("YRI and LWK should be in same region (Africa)")
	}
}

func TestSameRegion_SameEurope(t *testing.T) {
	gb := &GenomicsBenchmark{}
	if !gb.sameRegion("CEU", "GBR") {
		t.Error("CEU and GBR should be in same region (Europe)")
	}
}

func TestSameRegion_DifferentRegions(t *testing.T) {
	gb := &GenomicsBenchmark{}
	if gb.sameRegion("YRI", "CEU") {
		t.Error("YRI (Africa) and CEU (Europe) should not be in same region")
	}
}

func TestSameRegion_SameEastAsia(t *testing.T) {
	gb := &GenomicsBenchmark{}
	if !gb.sameRegion("CHB", "JPT") {
		t.Error("CHB and JPT should be in same region (East Asia)")
	}
}

func TestSameRegion_SameAmericas(t *testing.T) {
	gb := &GenomicsBenchmark{}
	if !gb.sameRegion("MXL", "PUR") {
		t.Error("MXL and PUR should be in same region (Americas)")
	}
}

func TestExtractPopulationCall_WithPopulationLine(t *testing.T) {
	gb := &GenomicsBenchmark{}
	output := "Processing sequence...\nPopulation: YRI\nConfidence: 0.85"
	got := gb.extractPopulationCall(output)
	if got != "YRI" {
		t.Errorf("extractPopulationCall = %q, want \"YRI\"", got)
	}
}

func TestExtractPopulationCall_NoPopulationLine(t *testing.T) {
	gb := &GenomicsBenchmark{}
	output := "No population info here"
	got := gb.extractPopulationCall(output)
	if got != "Unknown" {
		t.Errorf("extractPopulationCall = %q, want \"Unknown\"", got)
	}
}

func TestExtractPopulationCall_Empty(t *testing.T) {
	gb := &GenomicsBenchmark{}
	got := gb.extractPopulationCall("")
	if got != "Unknown" {
		t.Errorf("extractPopulationCall(\"\") = %q, want \"Unknown\"", got)
	}
}

func TestNewGenomicsBenchmark_NotNil(t *testing.T) {
	gb := NewGenomicsBenchmark()
	if gb == nil {
		t.Error("NewGenomicsBenchmark() returned nil")
	}
}

func TestNewGenomicsBenchmark_HasTools(t *testing.T) {
	gb := NewGenomicsBenchmark()
	if len(gb.Tools) == 0 {
		t.Error("NewGenomicsBenchmark() has no tools")
	}
}

func TestNewGenomicsBenchmark_HasPopulations(t *testing.T) {
	gb := NewGenomicsBenchmark()
	if len(gb.Populations) == 0 {
		t.Error("NewGenomicsBenchmark() has no populations")
	}
}

func TestNewGenomicsBenchmark_HasSequenceLengths(t *testing.T) {
	gb := NewGenomicsBenchmark()
	if len(gb.SequenceLengths) == 0 {
		t.Error("NewGenomicsBenchmark() has no sequence lengths")
	}
}

func TestParsePhantomResults_ExactMatch(t *testing.T) {
	gb := NewGenomicsBenchmark()
	output := "Processing...\nPopulation: YRI\nConfidence: 0.9"
	if got := gb.parsePhantomResults(output, "YRI"); got != 1.0 {
		t.Errorf("parsePhantomResults(exact match) = %f, want 1.0", got)
	}
}

func TestParsePhantomResults_RegionalMatch(t *testing.T) {
	gb := NewGenomicsBenchmark()
	output := "Population: LWK"
	// YRI and LWK are both African
	if got := gb.parsePhantomResults(output, "YRI"); got != 0.8 {
		t.Errorf("parsePhantomResults(regional match) = %f, want 0.8", got)
	}
}

func TestParsePhantomResults_NoMatch(t *testing.T) {
	gb := NewGenomicsBenchmark()
	output := "Population: CEU"
	// CEU (European) vs YRI (African) - different regions
	if got := gb.parsePhantomResults(output, "YRI"); got != 0.0 {
		t.Errorf("parsePhantomResults(no match) = %f, want 0.0", got)
	}
}

func TestParsePhantomResults_Empty(t *testing.T) {
	gb := NewGenomicsBenchmark()
	if got := gb.parsePhantomResults("", "YRI"); got != 0.0 {
		t.Errorf("parsePhantomResults(\"\") = %f, want 0.0", got)
	}
}

// ---------- readSequenceFromFasta ----------

func TestReadSequenceFromFasta_ValidFile(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := filepath.Join(tmpDir, "test.fasta")
	content := ">SEQ001 population=YRI\nATCGATCGATCG\nTTTAAACCC\n"
	_ = os.WriteFile(fastaPath, []byte(content), 0600)

	gb := NewGenomicsBenchmark()
	seq, err := gb.readSequenceFromFasta(fastaPath)
	if err != nil {
		t.Fatalf("readSequenceFromFasta error: %v", err)
	}
	if seq != "ATCGATCGATCGTTTAAACCC" {
		t.Errorf("sequence = %q, want \"ATCGATCGATCGTTTAAACCC\"", seq)
	}
}

func TestReadSequenceFromFasta_NotFound(t *testing.T) {
	gb := NewGenomicsBenchmark()
	_, err := gb.readSequenceFromFasta("/nonexistent/seq.fasta")
	if err == nil {
		t.Error("expected error for nonexistent file")
	}
}

func TestReadSequenceFromFasta_HeaderOnly(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := filepath.Join(tmpDir, "header_only.fasta")
	_ = os.WriteFile(fastaPath, []byte(">HEADER_ONLY\n"), 0600)

	gb := NewGenomicsBenchmark()
	seq, err := gb.readSequenceFromFasta(fastaPath)
	if err != nil {
		t.Fatalf("error: %v", err)
	}
	if seq != "" {
		t.Errorf("expected empty sequence for header-only FASTA, got %q", seq)
	}
}

// ---------- CreateBenchmarkDatasets ----------

func TestCreateBenchmarkDatasets_CreatesOutputDir(t *testing.T) {
	tmpDir := t.TempDir()
	gb := &GenomicsBenchmark{
		DataDir:         tmpDir + "/nodata",
		OutputDir:       filepath.Join(tmpDir, "output"),
		Tools:           []string{"ADMIXTURE", "fastSTRUCTURE", "PLINK", "PhantomGenomics"},
		Populations:     []string{"YRI", "CEU"},
		SequenceLengths: []int{1000},
		Results:         []BenchmarkResult{},
	}
	if err := gb.CreateBenchmarkDatasets(); err != nil {
		t.Fatalf("CreateBenchmarkDatasets error: %v", err)
	}
	if _, err := os.Stat(gb.OutputDir); err != nil {
		t.Error("expected output dir to exist")
	}
}

func TestCreateBenchmarkDatasets_BadOutputDir(t *testing.T) {
	gb := &GenomicsBenchmark{
		DataDir:         "/nodata",
		OutputDir:       "/nonexistent/path/xyz/output",
		Tools:           []string{"ADMIXTURE"},
		Populations:     []string{"YRI"},
		SequenceLengths: []int{1000},
		Results:         []BenchmarkResult{},
	}
	if err := gb.CreateBenchmarkDatasets(); err == nil {
		t.Error("expected error for bad output dir")
	}
}

// ---------- GenerateBenchmarkReport ----------

func TestGenerateBenchmarkReport_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	gb := &GenomicsBenchmark{
		DataDir:         tmpDir + "/nodata",
		OutputDir:       tmpDir,
		Tools:           []string{"PhantomGenomics", "ADMIXTURE"},
		Populations:     []string{"YRI", "CEU"},
		SequenceLengths: []int{1000},
		Results:         []BenchmarkResult{},
	}
	if err := gb.GenerateBenchmarkReport(); err != nil {
		t.Fatalf("GenerateBenchmarkReport error: %v", err)
	}
	if _, err := os.Stat(filepath.Join(tmpDir, "GENOMICS_BENCHMARKING_REPORT.md")); err != nil {
		t.Error("expected report file to exist")
	}
}

func TestGenerateBenchmarkReport_WithResults(t *testing.T) {
	tmpDir := t.TempDir()
	gb := &GenomicsBenchmark{
		DataDir:         tmpDir + "/nodata",
		OutputDir:       tmpDir,
		Tools:           []string{"PhantomGenomics"},
		Populations:     []string{"YRI"},
		SequenceLengths: []int{1000},
		Results: []BenchmarkResult{
			{Tool: "PhantomGenomics", Population: "YRI", SequenceLength: 1000, Accuracy: 0.85},
		},
	}
	if err := gb.GenerateBenchmarkReport(); err != nil {
		t.Fatalf("GenerateBenchmarkReport(with results) error: %v", err)
	}
}

func TestGenerateBenchmarkReport_BadOutputDir(t *testing.T) {
	gb := &GenomicsBenchmark{
		OutputDir: "/nonexistent/path/xyz",
		Tools:     []string{"ADMIXTURE"},
		Results:   []BenchmarkResult{},
	}
	if err := gb.GenerateBenchmarkReport(); err == nil {
		t.Error("expected error for bad output dir")
	}
}

// TestRunPhantomGenomicsBenchmark_EmptyPopulations covers function entry/exit with empty loop
func TestRunPhantomGenomicsBenchmark_EmptyPopulations(t *testing.T) {
	gb := &GenomicsBenchmark{
		Populations:     []string{},
		SequenceLengths: []int{},
		Results:         []BenchmarkResult{},
	}
	if err := gb.RunPhantomGenomicsBenchmark(); err != nil {
		t.Errorf("RunPhantomGenomicsBenchmark() error: %v", err)
	}
}

// TestRunPhantomGenomicsBenchmark_WithNonexistentFiles covers the error path (continue) in loop
func TestRunPhantomGenomicsBenchmark_WithNonexistentFiles(t *testing.T) {
	tmpDir := t.TempDir()
	gb := &GenomicsBenchmark{
		DataDir:         filepath.Join(tmpDir, "nonexistent"),
		Populations:     []string{"YRI"},
		SequenceLengths: []int{1000},
		Results:         []BenchmarkResult{},
	}
	// All sequences will fail (nonexistent DataDir), so all iterations hit 'continue'
	if err := gb.RunPhantomGenomicsBenchmark(); err != nil {
		t.Errorf("RunPhantomGenomicsBenchmark() error: %v", err)
	}
	// No results should be added since all iterations were skipped
	if len(gb.Results) != 0 {
		t.Errorf("expected 0 results for nonexistent files, got %d", len(gb.Results))
	}
}

// TestPrepareADMIXTUREDatasets_WithRealFiles covers the success path (lines 115-118)
func TestPrepareADMIXTUREDatasets_WithRealFiles(t *testing.T) {
	tmpDir := t.TempDir()

	// Create FASTA files in DataDir/pop/ directory structure
	pop := "YRI"
	popDir := filepath.Join(tmpDir, "data", pop)
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("failed to create pop dir: %v", err)
	}

	// Create sample FASTA file
	filename := fmt.Sprintf("%s_PARALLEL_VCF_L1000_S01.fasta", pop)
	fastaPath := filepath.Join(popDir, filename)
	content := fmt.Sprintf(">%s\nATCGATCGATCGATCG\n", filename)
	if err := os.WriteFile(fastaPath, []byte(content), 0600); err != nil {
		t.Fatalf("failed to create FASTA file: %v", err)
	}

	// Create the admixture subdirectory
	admixDir := filepath.Join(tmpDir, "output", "admixture")
	if err := os.MkdirAll(admixDir, 0750); err != nil {
		t.Fatalf("failed to create admixture dir: %v", err)
	}

	gb := &GenomicsBenchmark{
		DataDir:         filepath.Join(tmpDir, "data"),
		OutputDir:       filepath.Join(tmpDir, "output"),
		Populations:     []string{pop},
		SequenceLengths: []int{1000},
		Results:         []BenchmarkResult{},
	}

	if err := gb.prepareADMIXTUREDatasets(); err != nil {
		t.Errorf("prepareADMIXTUREDatasets() error: %v", err)
	}
}
