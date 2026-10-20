//go:build !release

package cmd

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

// makeSyntheticVCF writes a minimal VCF file with numVariants variant records
// and returns the file path. All alleles use valid ATGC bases.
func makeSyntheticVCF(t *testing.T, numVariants int) string {
	t.Helper()
	f, err := os.CreateTemp(t.TempDir(), "test_*.vcf")
	if err != nil {
		t.Fatalf("create temp VCF: %v", err)
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
	return f.Name()
}

// --- Command structure tests ---

func TestTrainCommand(t *testing.T) {
	if trainCmd == nil {
		t.Fatal("trainCmd should not be nil")
	}
	if trainCmd.Use != "train" {
		t.Errorf("expected Use 'train', got %q", trainCmd.Use)
	}
	if trainCmd.RunE == nil {
		t.Error("trainCmd should have a RunE function")
	}
	for _, name := range []string{"vcf", "out", "model-type", "kmer-size"} {
		if trainCmd.Flags().Lookup(name) == nil {
			t.Errorf("trainCmd should have --%s flag", name)
		}
	}
}

func TestTrainCommand_InheritsPopulation(t *testing.T) {
	// --population is a persistent flag on rootCmd; trainCmd must inherit it.
	if trainCmd.InheritedFlags().Lookup("population") == nil {
		t.Error("trainCmd should inherit --population from rootCmd")
	}
}

func TestTrainCommand_IsRegistered(t *testing.T) {
	found := false
	for _, cmd := range rootCmd.Commands() {
		if cmd.Use == "train" {
			found = true
			break
		}
	}
	if !found {
		t.Error("trainCmd should be registered with rootCmd")
	}
}

func TestTrainCommand_FlagDefaults(t *testing.T) {
	modelTypeFlag := trainCmd.Flags().Lookup("model-type")
	if modelTypeFlag == nil {
		t.Fatal("--model-type flag should exist")
	}
	if modelTypeFlag.DefValue != "adaptive" {
		t.Errorf("expected default model-type 'adaptive', got %q", modelTypeFlag.DefValue)
	}

	kmerSizeFlag := trainCmd.Flags().Lookup("kmer-size")
	if kmerSizeFlag == nil {
		t.Fatal("--kmer-size flag should exist")
	}
	if kmerSizeFlag.DefValue != "4" {
		t.Errorf("expected default kmer-size '4', got %q", kmerSizeFlag.DefValue)
	}
}

// --- runTrain validation tests ---

func TestRunTrain_NoVCFFiles(t *testing.T) {
	err := runTrain(nil, "out.hmm.json", "adaptive", "", 4)
	if err == nil {
		t.Fatal("expected error when no VCF files provided")
	}
	if !strings.Contains(err.Error(), "--vcf") {
		t.Errorf("expected error to mention --vcf, got: %v", err)
	}
}

func TestRunTrain_InvalidKmerSize(t *testing.T) {
	for _, k := range []int{0, 1, 9, 100} {
		err := runTrain([]string{"x.vcf"}, "out.hmm.json", "adaptive", "", k)
		if err == nil {
			t.Errorf("expected error for k=%d", k)
		}
	}
}

func TestRunTrain_InvalidModelType(t *testing.T) {
	vcf := makeSyntheticVCF(t, 1)
	err := runTrain([]string{vcf}, filepath.Join(t.TempDir(), "out.hmm.json"), "bogus", "", 4)
	if err == nil {
		t.Fatal("expected error for unknown model type")
	}
	if !strings.Contains(err.Error(), "bogus") {
		t.Errorf("expected error to mention the bad model type, got: %v", err)
	}
}

func TestRunTrain_MissingVCFFile(t *testing.T) {
	err := runTrain(
		[]string{"/nonexistent/path/file.vcf"},
		filepath.Join(t.TempDir(), "out.hmm.json"),
		"adaptive", "", 4,
	)
	if err == nil {
		t.Fatal("expected error for missing VCF file")
	}
}

// --- Integration tests ---

func TestRunTrain_AdaptiveModel(t *testing.T) {
	vcf := makeSyntheticVCF(t, 200)
	outPath := filepath.Join(t.TempDir(), "adaptive.hmm.json")

	if err := runTrain([]string{vcf}, outPath, "adaptive", "YRI", 4); err != nil {
		t.Fatalf("runTrain adaptive: %v", err)
	}

	// Output file must exist.
	if _, err := os.Stat(outPath); err != nil {
		t.Fatalf("output file not created: %v", err)
	}

	// Must load back via the unified LoadModel.
	model, err := authenticity.LoadModel(outPath)
	if err != nil {
		t.Fatalf("LoadModel: %v", err)
	}
	if !model.IsTrained() {
		t.Error("loaded model should report IsTrained() == true")
	}

	// Must generate a valid sequence.
	result, err := model.GenerateSequence(100)
	if err != nil {
		t.Fatalf("GenerateSequence: %v", err)
	}
	if len(result.Sequence) != 100 {
		t.Errorf("expected 100-base sequence, got %d", len(result.Sequence))
	}
}

func TestRunTrain_KmerModel(t *testing.T) {
	vcf := makeSyntheticVCF(t, 200)
	outPath := filepath.Join(t.TempDir(), "kmer.hmm.json")

	if err := runTrain([]string{vcf}, outPath, "kmer", "", 4); err != nil {
		t.Fatalf("runTrain kmer: %v", err)
	}

	if _, err := os.Stat(outPath); err != nil {
		t.Fatalf("output file not created: %v", err)
	}

	model, err := authenticity.LoadModel(outPath)
	if err != nil {
		t.Fatalf("LoadModel: %v", err)
	}
	if !model.IsTrained() {
		t.Error("loaded kmer model should report IsTrained() == true")
	}
}

func TestRunTrain_MultipleVCFFiles(t *testing.T) {
	vcf1 := makeSyntheticVCF(t, 100)
	vcf2 := makeSyntheticVCF(t, 100)
	outPath := filepath.Join(t.TempDir(), "multi.hmm.json")

	if err := runTrain([]string{vcf1, vcf2}, outPath, "adaptive", "", 4); err != nil {
		t.Fatalf("runTrain multi-VCF: %v", err)
	}
	if _, err := os.Stat(outPath); err != nil {
		t.Fatalf("output file not created: %v", err)
	}
}

// --- Unit tests for helper functions ---

func TestAppendVCFAlleles(t *testing.T) {
	tests := []struct {
		ref, alt string
		want     string
	}{
		{"ATCG", "ATCC", "ATCGATCC"},
		{"A", ".", "A"},
		{"A", "T,C", "AT"}, // multi-allelic: only first ALT
		{"", "T", "T"},
		{"A", "", "A"},
	}
	for _, tt := range tests {
		var buf strings.Builder
		appendVCFAlleles(&buf, tt.ref, tt.alt)
		if got := buf.String(); got != tt.want {
			t.Errorf("appendVCFAlleles(%q, %q) = %q, want %q", tt.ref, tt.alt, got, tt.want)
		}
	}
}

func TestOpenVCFScanner_PlainFile(t *testing.T) {
	f, err := os.CreateTemp(t.TempDir(), "*.vcf")
	if err != nil {
		t.Fatal(err)
	}
	fmt.Fprintln(f, "##header")
	fmt.Fprintln(f, "chr1\t1\t.\tA\tT")
	f.Close()

	f2, err := os.Open(f.Name())
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()

	sc, cleanup, err := openVCFScanner(f2, f.Name())
	if err != nil {
		t.Fatalf("openVCFScanner: %v", err)
	}
	if cleanup != nil {
		defer cleanup() //nolint:errcheck
	}

	var lines []string
	for sc.Scan() {
		lines = append(lines, sc.Text())
	}
	if len(lines) != 2 {
		t.Errorf("expected 2 lines, got %d", len(lines))
	}
}

func TestVCFChromosome(t *testing.T) {
	chr := &vcfChromosome{name: "chr1", sequence: []byte("ATCGATCG")}
	if chr.GetName() != "chr1" {
		t.Errorf("expected name chr1, got %q", chr.GetName())
	}
	if string(chr.GetSequence()) != "ATCGATCG" {
		t.Errorf("unexpected sequence: %q", chr.GetSequence())
	}
	if chr.GetLength() != 8 {
		t.Errorf("expected length 8, got %d", chr.GetLength())
	}
}

func TestVCFGenome(t *testing.T) {
	genome := &vcfGenome{
		chromosomes: map[string]*vcfChromosome{
			"chr1": {name: "chr1", sequence: []byte("ATCGATCG")},
		},
		names: []string{"chr1"},
	}

	names := genome.GetChromosomeNames()
	if len(names) != 1 || names[0] != "chr1" {
		t.Errorf("expected [chr1], got %v", names)
	}

	chr, err := genome.LoadChromosome("chr1")
	if err != nil {
		t.Fatalf("LoadChromosome: %v", err)
	}
	if chr.GetName() != "chr1" {
		t.Errorf("expected name chr1, got %q", chr.GetName())
	}
	if string(chr.GetSequence()) != "ATCGATCG" {
		t.Errorf("unexpected sequence: %q", chr.GetSequence())
	}
	if chr.GetLength() != 8 {
		t.Errorf("expected length 8, got %d", chr.GetLength())
	}

	_, err = genome.LoadChromosome("chrX")
	if err == nil {
		t.Error("expected error for missing chromosome")
	}
}

// modelData returns the data section of a saved model file (the envelope's
// serialized_at timestamp differs between saves).
func modelData(t *testing.T, path string) string {
	t.Helper()
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	var env authenticity.ModelEnvelope
	if err := json.Unmarshal(raw, &env); err != nil {
		t.Fatal(err)
	}
	return string(env.Data)
}

// Streaming the VCF in 7-base chunks (with the k-1 carry) must give exactly
// the model of one-shot training: no k-mer spanning a chunk boundary is lost.
func TestStreamVCFToAdaptive_ChunkedEqualsOneShot(t *testing.T) {
	vcf := makeSyntheticVCF(t, 300)
	for _, k := range []int{2, 4, 8} {
		t.Run(fmt.Sprintf("k=%d", k), func(t *testing.T) {
			dir := t.TempDir()
			paths := map[int]string{}
			for _, chunk := range []int{7, 1 << 30} {
				m := authenticity.NewAdaptiveHMMModel(k)
				if err := streamVCFToAdaptive(vcf, m, chunk); err != nil {
					t.Fatal(err)
				}
				paths[chunk] = filepath.Join(dir, fmt.Sprintf("%d.hmm.json", chunk))
				if err := m.Save(paths[chunk]); err != nil {
					t.Fatal(err)
				}
			}
			if modelData(t, paths[7]) != modelData(t, paths[1<<30]) {
				t.Error("chunked training differs from one-shot training")
			}
		})
	}
}

// A model trained by `train` on a tiny VCF (16 bases) generates every base
// (it used to generate only C).
func TestRunTrain_TinyVCFNotDegenerate(t *testing.T) {
	for _, modelType := range []string{"adaptive", "kmer"} {
		t.Run(modelType, func(t *testing.T) {
			vcf := makeSyntheticVCF(t, 2)
			out := filepath.Join(t.TempDir(), "m.hmm.json")
			if err := runTrain([]string{vcf}, out, modelType, "", 4); err != nil {
				t.Fatal(err)
			}
			m, err := authenticity.LoadModel(out)
			if err != nil {
				t.Fatal(err)
			}
			m.(interface{ SetSeed(int64) }).SetSeed(1)
			res, err := m.GenerateSequence(10_000)
			if err != nil {
				t.Fatal(err)
			}
			counts := map[byte]int{}
			for _, b := range res.Sequence {
				counts[b]++
			}
			t.Logf("base counts: %v", counts)
			for _, b := range []byte("ACGT") {
				if counts[b] < 1000 {
					t.Errorf("base %c: %d of 10000", b, counts[b])
				}
			}
		})
	}
}
