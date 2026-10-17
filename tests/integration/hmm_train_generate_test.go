//go:build !release

package integration

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// TestHMMTrainGeneratePipeline exercises the complete round-trip:
//
//  1. Write a minimal synthetic VCF (~50 variants) to a temp directory.
//  2. Train an adaptive HMM model from the VCF using the authenticity package
//     (mirrors the cmd/train adaptive path) and save to model.hmm.json.
//  3. Verify model.hmm.json exists and contains a valid ModelEnvelope JSON object
//     with a non-empty model_type field.
//  4. Build a SequencerConfig pointing at the trained model and call engine.Generate().
//  5. Assert: FASTQ file exists; every sequence line contains only ATGCN characters;
//     quality line length equals ReadLength.
func TestHMMTrainGeneratePipeline(t *testing.T) {
	tmpDir := t.TempDir()
	vcfPath := filepath.Join(tmpDir, "test.vcf")
	modelPath := filepath.Join(tmpDir, "model.hmm.json")
	outputDir := filepath.Join(tmpDir, "output")

	// Step 1: write synthetic VCF
	writeSyntheticVCFForHMM(t, vcfPath, 50)

	// Step 2: train adaptive HMM model and save to disk
	if err := trainAdaptiveHMMFromVCF(t, vcfPath, modelPath); err != nil {
		t.Fatalf("train step failed: %v", err)
	}

	// Step 3: verify model file is valid JSON with ModelEnvelope shape
	raw, err := os.ReadFile(modelPath)
	if err != nil {
		t.Fatalf("model file not found after training: %v", err)
	}
	var envelope struct {
		ModelType string `json:"model_type"`
	}
	if err := json.Unmarshal(raw, &envelope); err != nil {
		t.Fatalf("model file is not valid JSON: %v", err)
	}
	if envelope.ModelType == "" {
		t.Fatalf("model file missing 'model_type' field; content prefix: %s", raw[:minLen(200, len(raw))])
	}

	// Step 4: create engine with trained model
	cfg := types.SequencerConfig{
		HMMModelPath:  modelPath,
		SampleName:    "HMMIntegration",
		Lane:          1,
		OutputDir:     outputDir,
		ReadLength:    100,
		ReadsPerBatch: 5,
		TotalBatches:  1,
		PairedEnd:     false,
	}
	engine := genomics.NewEngine(cfg)
	defer engine.Close()

	// Step 5: generate FASTQ
	if err := engine.Generate(); err != nil {
		t.Fatalf("generation step failed: %v", err)
	}

	// Step 6: assert FASTQ file exists and has valid content
	fastqPath := filepath.Join(outputDir, fmt.Sprintf("HMMIntegration_S1_L%03d_%03d.fastq", cfg.Lane, 1))
	if _, err := os.Stat(fastqPath); os.IsNotExist(err) {
		t.Fatalf("expected FASTQ file not found: %s", fastqPath)
	}
	verifyFASTQContent(t, fastqPath, cfg.ReadLength)
}

// writeSyntheticVCFForHMM creates a minimal plain-text VCF file with numVariants records.
func writeSyntheticVCFForHMM(t *testing.T, path string, numVariants int) {
	t.Helper()
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("create VCF: %v", err)
	}
	defer f.Close()

	fmt.Fprintln(f, "##fileformat=VCFv4.1")
	fmt.Fprintln(f, "#CHROM\tPOS\tID\tREF\tALT")

	refs := []string{"ATCG", "GCTA", "TATC", "CGAT", "AATG", "CCAT", "GGCA", "TTAC"}
	alts := []string{"ATCC", "GCTG", "TATA", "CGAC", "AATT", "CCAG", "GGCG", "TTAA"}

	for i := 0; i < numVariants; i++ {
		ref := refs[i%len(refs)]
		alt := alts[i%len(alts)]
		fmt.Fprintf(f, "chr1\t%d\t.\t%s\t%s\n", (i+1)*10, ref, alt)
	}
}

// trainAdaptiveHMMFromVCF replicates the adaptive train path from cmd/train.go using
// the authenticity package directly to avoid subprocess overhead.
func trainAdaptiveHMMFromVCF(t *testing.T, vcfPath, outPath string) error {
	t.Helper()
	model := authenticity.NewAdaptiveHMMModel(4)
	if model == nil {
		return fmt.Errorf("failed to create adaptive HMM model")
	}

	f, err := os.Open(vcfPath)
	if err != nil {
		return fmt.Errorf("open VCF: %w", err)
	}
	defer f.Close()

	const chunkSize = 100 * 1024
	var buf []byte
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.SplitN(line, "\t", 6)
		if len(fields) < 5 {
			continue
		}
		buf = append(buf, []byte(fields[3])...) // REF
		if parts := strings.SplitN(fields[4], ",", 2); len(parts) > 0 {
			buf = append(buf, []byte(parts[0])...) // first ALT
		}
		if len(buf) >= chunkSize {
			if err := model.TrainFromSequence(buf); err != nil {
				return fmt.Errorf("TrainFromSequence: %w", err)
			}
			buf = buf[:0]
		}
	}
	if len(buf) > 0 {
		if err := model.TrainFromSequence(buf); err != nil {
			return fmt.Errorf("TrainFromSequence (final): %w", err)
		}
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("scan VCF: %w", err)
	}
	return model.Save(outPath)
}

var hmmDNARE = regexp.MustCompile(`^[ATGCNatgcn]+$`)

// verifyFASTQContent checks that every record in the FASTQ file has a valid
// sequence (ATGCN only) and that sequence/quality lengths match readLength.
func verifyFASTQContent(t *testing.T, path string, readLength int) {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open FASTQ for verification: %v", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	lineNum := 0
	records := 0
	for scanner.Scan() {
		line := scanner.Text()
		switch lineNum % 4 {
		case 0: // header
			if !strings.HasPrefix(line, "@") {
				t.Errorf("FASTQ line %d: expected '@' header, got %q", lineNum+1, line)
			}
		case 1: // sequence
			if !hmmDNARE.MatchString(line) {
				t.Errorf("FASTQ line %d: sequence contains non-DNA characters: %q", lineNum+1, line)
			}
			if len(line) != readLength {
				t.Errorf("FASTQ line %d: sequence length %d, want %d", lineNum+1, len(line), readLength)
			}
		case 2: // separator
			if line != "+" {
				t.Errorf("FASTQ line %d: expected '+', got %q", lineNum+1, line)
			}
		case 3: // quality
			if len(line) != readLength {
				t.Errorf("FASTQ line %d: quality length %d != read length %d", lineNum+1, len(line), readLength)
			}
			records++
		}
		lineNum++
	}
	if err := scanner.Err(); err != nil {
		t.Fatalf("scan FASTQ: %v", err)
	}
	if records == 0 {
		t.Error("FASTQ file is empty — no records found")
	}
}

func minLen(a, b int) int {
	if a < b {
		return a
	}
	return b
}
