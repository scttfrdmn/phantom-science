//go:build !release
// +build !release

package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

func TestSequenceToGenotype_ABase(t *testing.T) {
	pb := &PLINKBenchmark{}
	// Build a sequence with A at positions 1000, 2000, 3000
	seq := make([]byte, 3001)
	for i := range seq {
		seq[i] = 'T'
	}
	seq[999] = 'A'  // position 1000 (1-based)
	seq[1999] = 'A' // position 2000
	seq[2999] = 'A' // position 3000
	got := pb.sequenceToGenotype(string(seq))
	if got != "A A A A A A" {
		t.Errorf("sequenceToGenotype = %q, want \"A A A A A A\"", got)
	}
}

func TestSequenceToGenotype_TooShort(t *testing.T) {
	pb := &PLINKBenchmark{}
	// Short sequence - all positions will be missing
	got := pb.sequenceToGenotype("ATCG")
	if got != "0 0 0 0 0 0" {
		t.Errorf("sequenceToGenotype(short) = %q, want \"0 0 0 0 0 0\"", got)
	}
}

func TestSequenceToGenotype_MixedBases(t *testing.T) {
	pb := &PLINKBenchmark{}
	seq := make([]byte, 3001)
	for i := range seq {
		seq[i] = 'N' // invalid
	}
	seq[999] = 'G'  // position 1000
	seq[1999] = 'C' // position 2000
	seq[2999] = 'T' // position 3000
	got := pb.sequenceToGenotype(string(seq))
	if got != "G G C C T T" {
		t.Errorf("sequenceToGenotype = %q, want \"G G C C T T\"", got)
	}
}

// ---------- NewPLINKBenchmark ----------

func TestNewPLINKBenchmark_NotNil(t *testing.T) {
	pb := NewPLINKBenchmark()
	if pb == nil {
		t.Error("NewPLINKBenchmark() returned nil")
	}
}

func TestNewPLINKBenchmark_HasPopulations(t *testing.T) {
	pb := NewPLINKBenchmark()
	if len(pb.Populations) == 0 {
		t.Error("NewPLINKBenchmark() has no populations")
	}
}

// ---------- readSequenceFromFasta ----------

func TestReadSequenceFromFasta_ValidFile(t *testing.T) {
	tmpDir := t.TempDir()
	fastaPath := filepath.Join(tmpDir, "test.fasta")
	content := ">SEQ001 population=YRI\natcgatcgatcg\nTTTAAACCC\n"
	_ = os.WriteFile(fastaPath, []byte(content), 0600)

	pb := &PLINKBenchmark{}
	seq, err := pb.readSequenceFromFasta(fastaPath)
	if err != nil {
		t.Fatalf("readSequenceFromFasta error: %v", err)
	}
	// Should be uppercase (ToUpper is called)
	if seq != "ATCGATCGATCGTTTAAACCC" {
		t.Errorf("sequence = %q, want \"ATCGATCGATCGTTTAAACCC\"", seq)
	}
}

func TestReadSequenceFromFasta_NotFound(t *testing.T) {
	pb := &PLINKBenchmark{}
	_, err := pb.readSequenceFromFasta("/nonexistent/seq.fasta")
	if err == nil {
		t.Error("expected error for nonexistent file")
	}
}

// ---------- ConvertSequencesToPED ----------

func TestConvertSequencesToPED_CreatesFiles(t *testing.T) {
	tmpDir := t.TempDir()
	pb := &PLINKBenchmark{
		DataDir:     tmpDir + "/nodata",
		OutputDir:   filepath.Join(tmpDir, "ped_output"),
		Populations: []string{"YRI", "CEU"},
		Results:     make([]PLINKResult, 0),
	}
	if err := pb.ConvertSequencesToPED(); err != nil {
		t.Fatalf("ConvertSequencesToPED error: %v", err)
	}
	if _, err := os.Stat(filepath.Join(pb.OutputDir, "benchmark_sequences.ped")); err != nil {
		t.Error("expected PED file to exist")
	}
	if _, err := os.Stat(filepath.Join(pb.OutputDir, "benchmark_sequences.map")); err != nil {
		t.Error("expected MAP file to exist")
	}
}

func TestConvertSequencesToPED_BadOutputDir(t *testing.T) {
	pb := &PLINKBenchmark{
		DataDir:     "/nodata",
		OutputDir:   "/nonexistent/xyz/output",
		Populations: []string{"YRI"},
		Results:     make([]PLINKResult, 0),
	}
	if err := pb.ConvertSequencesToPED(); err == nil {
		t.Error("expected error for bad output dir")
	}
}

// ---------- RunPLINKAnalysis ----------

func TestRunPLINKAnalysis_NoCrash(t *testing.T) {
	tmpDir := t.TempDir()
	pb := &PLINKBenchmark{
		DataDir:   tmpDir,
		OutputDir: tmpDir,
		PLINKPath: "/nonexistent/plink_binary",
		Results:   make([]PLINKResult, 0),
	}
	// RunPLINKAnalysis returns nil even if PLINK isn't installed (records error in results)
	err := pb.RunPLINKAnalysis()
	if err != nil {
		t.Fatalf("RunPLINKAnalysis returned unexpected error: %v", err)
	}
	if len(pb.Results) != 1 {
		t.Errorf("expected 1 result, got %d", len(pb.Results))
	}
}

// ---------- GenerateReport ----------

func TestGenerateReport_CreatesFile(t *testing.T) {
	tmpDir := t.TempDir()
	pb := &PLINKBenchmark{
		OutputDir: tmpDir,
		Results:   []PLINKResult{},
	}
	if err := pb.GenerateReport(); err != nil {
		t.Fatalf("GenerateReport error: %v", err)
	}
	if _, err := os.Stat(filepath.Join(tmpDir, "PLINK_BENCHMARK_REPORT.md")); err != nil {
		t.Error("expected PLINK benchmark report file to exist")
	}
}

func TestGenerateReport_WithResults(t *testing.T) {
	tmpDir := t.TempDir()
	pb := &PLINKBenchmark{
		OutputDir: tmpDir,
		Results: []PLINKResult{
			{Population: "YRI", ProcessingTime: time.Second, Success: true, OutputFiles: []string{"file1.txt"}},
			{Population: "CEU", ProcessingTime: time.Second * 2, Success: false, ErrorMessage: "PLINK not found"},
		},
	}
	if err := pb.GenerateReport(); err != nil {
		t.Fatalf("GenerateReport(with results) error: %v", err)
	}
}

func TestGenerateReport_BadOutputDir(t *testing.T) {
	pb := &PLINKBenchmark{
		OutputDir: "/nonexistent/xyz/output",
		Results:   []PLINKResult{},
	}
	if err := pb.GenerateReport(); err == nil {
		t.Error("expected error for bad output dir")
	}
}

// Ensure time import is used
var _ = time.Second

// TestRunPLINKAnalysis_WithMockBinary covers the else branch (success path) in RunPLINKAnalysis
func TestRunPLINKAnalysis_WithMockBinary(t *testing.T) {
	tmpDir := t.TempDir()

	// Create a mock PLINK binary (shell script) that exits successfully
	mockPlink := filepath.Join(tmpDir, "plink")
	if err := os.WriteFile(mockPlink, []byte("#!/bin/sh\nexit 0\n"), 0700); err != nil {
		t.Fatalf("failed to create mock plink script: %v", err)
	}

	pb := &PLINKBenchmark{
		DataDir:   tmpDir,
		OutputDir: tmpDir,
		PLINKPath: mockPlink,
		Results:   make([]PLINKResult, 0),
	}

	err := pb.RunPLINKAnalysis()
	if err != nil {
		t.Fatalf("RunPLINKAnalysis returned unexpected error: %v", err)
	}
	if len(pb.Results) != 1 {
		t.Errorf("expected 1 result, got %d", len(pb.Results))
	}
	if !pb.Results[0].Success {
		t.Errorf("expected successful result with mock binary, got: %v", pb.Results[0].ErrorMessage)
	}
}

// TestConvertSequencesToPED_WithRealFASTAFiles covers the success path inside the loop
func TestConvertSequencesToPED_WithRealFASTAFiles(t *testing.T) {
	tmpDir := t.TempDir()

	// Create FASTA files in the expected directory structure: DataDir/pop/pop_PARALLEL_VCF_Llen_Snum.fasta
	pop := "YRI"
	popDir := filepath.Join(tmpDir, "data", pop)
	if err := os.MkdirAll(popDir, 0750); err != nil {
		t.Fatalf("failed to create pop dir: %v", err)
	}

	// Create one sample FASTA file that's long enough for genotype extraction
	seq := strings.Repeat("ATCGATCGATCG", 300) // 3600 chars, well beyond position 3000
	for i := 1; i <= 1; i++ {
		filename := fmt.Sprintf("%s_PARALLEL_VCF_L1000_S%02d.fasta", pop, i)
		fastaPath := filepath.Join(popDir, filename)
		content := fmt.Sprintf(">%s\n%s\n", filename, seq)
		if err := os.WriteFile(fastaPath, []byte(content), 0600); err != nil {
			t.Fatalf("failed to create FASTA file: %v", err)
		}
	}

	pb := &PLINKBenchmark{
		DataDir:     filepath.Join(tmpDir, "data"),
		OutputDir:   filepath.Join(tmpDir, "output"),
		Populations: []string{pop},
		Results:     make([]PLINKResult, 0),
	}

	if err := pb.ConvertSequencesToPED(); err != nil {
		t.Fatalf("ConvertSequencesToPED error: %v", err)
	}

	// Verify PED and MAP files were created
	if _, err := os.Stat(filepath.Join(pb.OutputDir, "benchmark_sequences.ped")); err != nil {
		t.Error("expected PED file to exist")
	}
}
