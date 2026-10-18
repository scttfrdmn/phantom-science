//go:build !release
// +build !release

package testing

import (
	"bytes"
	"errors"
	"math"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
	"testing"
	"time"
)

func TestNewE2EFramework_NotNil(t *testing.T) {
	f := NewE2EFramework()
	if f == nil {
		t.Error("NewE2EFramework() returned nil")
	}
}

func TestNewE2EFramework_HasPlatformConfigs(t *testing.T) {
	f := NewE2EFramework()
	if len(f.platformConfigs) == 0 {
		t.Error("NewE2EFramework() has no platform configs")
	}
}

func TestCreatePlatformConfig_NotNil(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	if config == nil {
		t.Error("CreatePlatformConfig returned nil")
	}
}

func TestCreatePlatformConfig_CorrectPlatform(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(PacBioRevio, "YRI")
	if config.Platform != PacBioRevio {
		t.Errorf("CreatePlatformConfig platform = %v, want %v", config.Platform, PacBioRevio)
	}
}

func TestGenerateSyntheticData_NotNil(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	if data == nil {
		t.Error("GenerateSyntheticData returned nil data")
	}
}

func TestGenerateSyntheticData_HasSequences(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	if len(data.Sequences) == 0 {
		t.Error("GenerateSyntheticData returned data with no sequences")
	}
}

func TestCalculateGCContent_AllGC(t *testing.T) {
	f := NewE2EFramework()
	seqs := [][]byte{[]byte("GCGCGCGC")}
	got := f.calculateGCContent(seqs)
	if got != 1.0 {
		t.Errorf("calculateGCContent(all GC) = %f, want 1.0", got)
	}
}

func TestCalculateGCContent_AllAT(t *testing.T) {
	f := NewE2EFramework()
	seqs := [][]byte{[]byte("ATATAT")}
	got := f.calculateGCContent(seqs)
	if got != 0.0 {
		t.Errorf("calculateGCContent(all AT) = %f, want 0.0", got)
	}
}

func TestCalculateGCContent_Empty(t *testing.T) {
	f := NewE2EFramework()
	got := f.calculateGCContent([][]byte{})
	if got != 0.0 {
		t.Errorf("calculateGCContent([]) = %f, want 0.0", got)
	}
}

func TestAnalyzeAuthenticity_InRange(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	score, err := f.AnalyzeAuthenticity(data)
	if err != nil {
		t.Fatalf("AnalyzeAuthenticity error: %v", err)
	}
	if score < 0.0 || score > 100.0 {
		t.Errorf("AnalyzeAuthenticity score = %f, want in [0,100]", score)
	}
}

func TestCreateMLClassifier_NotNil(t *testing.T) {
	f := NewE2EFramework()
	classifier := f.CreateMLClassifier()
	if classifier == nil {
		t.Error("CreateMLClassifier returned nil")
	}
}

func TestMLClassifier_Train(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	realData, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	syntheticData, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	classifier := f.CreateMLClassifier()
	err = classifier.Train(realData, syntheticData)
	if err != nil {
		t.Errorf("Train error: %v", err)
	}
}

func TestMLClassifier_TestAccuracy_InRange(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	realData, _ := f.GenerateSyntheticData(config)
	syntheticData, _ := f.GenerateSyntheticData(config)
	classifier := f.CreateMLClassifier()
	_ = classifier.Train(realData, syntheticData)
	accuracy := classifier.TestAccuracy(realData, syntheticData)
	if accuracy < 0.0 || accuracy > 1.0 {
		t.Errorf("TestAccuracy = %f, want in [0,1]", accuracy)
	}
}

func TestMLClassifier_TestAccuracy_PacBioRevio(t *testing.T) {
	// PacBio reads are 15,000 bp — verify the classifier handles long reads and
	// produces an accuracy in [0,1]. Two independently-generated sets from the
	// same distribution should yield accuracy near 0.5.
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(PacBioRevio, "CEU")
	realData, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(PacBioRevio): %v", err)
	}
	syntheticData, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(PacBioRevio) second call: %v", err)
	}
	classifier := f.CreateMLClassifier()
	_ = classifier.Train(realData, syntheticData)
	accuracy := classifier.TestAccuracy(realData, syntheticData)
	if accuracy < 0.0 || accuracy > 1.0 {
		t.Errorf("PacBio TestAccuracy = %f, want in [0,1]", accuracy)
	}
}

func TestMLClassifier_TestAccuracy_ONTPromethION24(t *testing.T) {
	// ONT reads are 30,000 bp — verify the classifier handles ultra-long reads.
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(ONTPromethION24, "YRI")
	realData, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(ONTPromethION24): %v", err)
	}
	syntheticData, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(ONTPromethION24) second call: %v", err)
	}
	classifier := f.CreateMLClassifier()
	_ = classifier.Train(realData, syntheticData)
	accuracy := classifier.TestAccuracy(realData, syntheticData)
	if accuracy < 0.0 || accuracy > 1.0 {
		t.Errorf("ONT TestAccuracy = %f, want in [0,1]", accuracy)
	}
}

func TestRunBenchmark_NotNil(t *testing.T) {
	f := NewE2EFramework()
	benchmark := f.RunBenchmark(IlluminaNovaSeq6000)
	if benchmark == nil {
		t.Error("RunBenchmark returned nil")
	}
}

func TestRunBenchmark_CorrectPlatform(t *testing.T) {
	f := NewE2EFramework()
	benchmark := f.RunBenchmark(PacBioRevio)
	if benchmark.Platform != PacBioRevio {
		t.Errorf("RunBenchmark platform = %v, want %v", benchmark.Platform, PacBioRevio)
	}
}

// Regression: CPU usage was read from /cpu/classes/total:cpu-seconds, which is
// CPU time *available* (GOMAXPROCS × wall time), so it reported ~GOMAXPROCS×100%.
func TestRunBenchmark_CPUUsageIsShareOfMachine(t *testing.T) {
	f := NewE2EFramework()
	benchmark := f.RunBenchmark(PacBioRevio)
	if benchmark.CPUUsagePercent < 0 || benchmark.CPUUsagePercent > 100 {
		t.Errorf("CPUUsagePercent = %.1f, want in [0, 100]", benchmark.CPUUsagePercent)
	}
}

func TestCPUSecondsSinceStart_Monotonic(t *testing.T) {
	before := cpuSecondsSinceStart()
	deadline := time.Now().Add(50 * time.Millisecond)
	for time.Now().Before(deadline) {
	}
	if after := cpuSecondsSinceStart(); runtime.GOOS != "windows" && after <= before {
		t.Errorf("CPU time did not advance during a busy loop: before=%.4f after=%.4f", before, after)
	}
}

func TestGenerateBalancedDataset_NotNil(t *testing.T) {
	f := NewE2EFramework()
	realData, synData, err := f.GenerateBalancedDataset(IlluminaNovaSeq6000, 10)
	if err != nil {
		t.Fatalf("GenerateBalancedDataset error: %v", err)
	}
	if realData == nil || synData == nil {
		t.Error("GenerateBalancedDataset returned nil data")
	}
}

// Regression: the "real" set used to be a shallow copy of the synthetic set,
// so every read had an identical twin with the opposite label and classifier
// accuracy fell systematically below chance.
func TestGenerateBalancedDataset_IndependentSets(t *testing.T) {
	f := NewE2EFramework()
	realData, synData, err := f.GenerateBalancedDataset(IlluminaNovaSeq6000, 50)
	if err != nil {
		t.Fatalf("GenerateBalancedDataset error: %v", err)
	}
	if len(realData.Sequences) != 50 || len(synData.Sequences) != 50 {
		t.Fatalf("got %d real / %d synthetic reads, want 50 each",
			len(realData.Sequences), len(synData.Sequences))
	}

	synSeqs := make(map[string]bool, len(synData.Sequences))
	for _, s := range synData.Sequences {
		synSeqs[string(s)] = true
	}
	for i, s := range realData.Sequences {
		if synSeqs[string(s)] {
			t.Fatalf("real read %d also appears in the synthetic set", i)
		}
	}
	if &realData.QualityScores[0][0] == &synData.QualityScores[0][0] {
		t.Error("real and synthetic sets share quality-score storage")
	}
}

func TestGenerateBalancedDataset_RejectsNonPositiveSize(t *testing.T) {
	f := NewE2EFramework()
	if _, _, err := f.GenerateBalancedDataset(IlluminaNovaSeq6000, 0); err == nil {
		t.Error("expected error for numSamples = 0")
	}
}

// A trained classifier must score held-out data with the fitted model: train
// on unskewed vs all-GC reads, then check it labels fresh reads correctly.
func TestMLClassifier_TestAccuracy_UsesTrainedModel(t *testing.T) {
	f := NewE2EFramework()
	allGC := func(d *SyntheticData) *SyntheticData {
		out := *d
		out.Sequences = make([][]byte, len(d.Sequences))
		for i, s := range d.Sequences {
			out.Sequences[i] = bytes.Map(func(r rune) rune {
				switch r {
				case 'A':
					return 'G'
				case 'T':
					return 'C'
				}
				return r
			}, s)
		}
		return &out
	}

	trainReal, trainSyn, err := f.GenerateBalancedDataset(IlluminaNovaSeq6000, 100)
	if err != nil {
		t.Fatalf("GenerateBalancedDataset error: %v", err)
	}
	testReal, testSyn, err := f.GenerateBalancedDataset(IlluminaNovaSeq6000, 100)
	if err != nil {
		t.Fatalf("GenerateBalancedDataset error: %v", err)
	}

	c := f.CreateMLClassifier()
	if err := c.Train(trainReal, allGC(trainSyn)); err != nil {
		t.Fatalf("Train error: %v", err)
	}
	if acc := c.TestAccuracy(testReal, allGC(testSyn)); acc < 0.99 {
		t.Errorf("held-out accuracy = %.3f, want ≥ 0.99", acc)
	}
	// Swapping the labels must invert the result, which cross-validation on
	// the held-out set (ignoring the fitted model) would not do.
	if acc := c.TestAccuracy(allGC(testSyn), testReal); acc > 0.01 {
		t.Errorf("label-swapped held-out accuracy = %.3f, want ≤ 0.01", acc)
	}
}

// Tests for applyPopulationAdjustments (via CreatePlatformConfig)
func TestCreatePlatformConfig_YRI(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "YRI")
	if config == nil {
		t.Fatal("CreatePlatformConfig(YRI) returned nil")
	}
	if config.Population != "YRI" {
		t.Errorf("population = %q, want YRI", config.Population)
	}
}

func TestCreatePlatformConfig_CHB(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CHB")
	if config == nil {
		t.Fatal("CreatePlatformConfig(CHB) returned nil")
	}
	if config.Population != "CHB" {
		t.Errorf("population = %q, want CHB", config.Population)
	}
}

func TestCreatePlatformConfig_JPT(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "JPT")
	if config == nil {
		t.Fatal("CreatePlatformConfig(JPT) returned nil")
	}
}

func TestCreatePlatformConfig_MXL(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "MXL")
	if config == nil {
		t.Fatal("CreatePlatformConfig(MXL) returned nil")
	}
}

func TestCreatePlatformConfig_UnknownPopulation(t *testing.T) {
	f := NewE2EFramework()
	// Unknown populations should still work (no specific adjustments)
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "UNKNOWN")
	if config == nil {
		t.Fatal("CreatePlatformConfig(UNKNOWN) returned nil")
	}
}

// Tests for applyPlatformSpecificModeling (via GenerateSyntheticData with different platforms)
func TestGenerateSyntheticData_PacBioRevio(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(PacBioRevio, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(PacBioRevio) error: %v", err)
	}
	if data == nil {
		t.Error("GenerateSyntheticData(PacBioRevio) returned nil")
	}
}

func TestGenerateSyntheticData_ONTPromethION24(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(ONTPromethION24, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(ONTPromethION24) error: %v", err)
	}
	if data == nil {
		t.Error("GenerateSyntheticData(ONTPromethION24) returned nil")
	}
}

func TestGenerateSyntheticData_IlluminaNextSeq(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNextSeq2000, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(IlluminaNextSeq) error: %v", err)
	}
	if data == nil {
		t.Error("GenerateSyntheticData(IlluminaNextSeq) returned nil")
	}
}

func TestGenerateSyntheticData_ElementAVITI(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(ElementAVITI, "CEU")
	if config == nil {
		// ElementAVITI might not have a config; skip if so
		t.Skip("ElementAVITI not configured")
	}
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(ElementAVITI) error: %v", err)
	}
	if data == nil {
		t.Error("GenerateSyntheticData(ElementAVITI) returned nil")
	}
}

// Tests for applyPopulationSpecificFeatures (via GenerateSyntheticData with different populations)
func TestGenerateSyntheticData_YRI(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "YRI")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(YRI) error: %v", err)
	}
	if data == nil {
		t.Error("GenerateSyntheticData(YRI) returned nil")
	}
}

func TestGenerateSyntheticData_CHB(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CHB")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(CHB) error: %v", err)
	}
	if data == nil {
		t.Error("GenerateSyntheticData(CHB) returned nil")
	}
}

func TestGenerateSyntheticData_JPT(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "JPT")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(JPT) error: %v", err)
	}
	if data == nil {
		t.Error("GenerateSyntheticData(JPT) returned nil")
	}
}

func TestGenerateSyntheticData_MXL(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "MXL")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(MXL) error: %v", err)
	}
	if data == nil {
		t.Error("GenerateSyntheticData(MXL) returned nil")
	}
}

// Test AnalyzeAuthenticity with different platforms to cover analyzePlatformAuthenticity branches
func TestAnalyzeAuthenticity_PacBioRevio(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(PacBioRevio, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	score, err := f.AnalyzeAuthenticity(data)
	if err != nil {
		t.Fatalf("AnalyzeAuthenticity(PacBioRevio) error: %v", err)
	}
	if score < 0.0 || score > 100.0 {
		t.Errorf("AnalyzeAuthenticity score = %f, want in [0,100]", score)
	}
}

func TestAnalyzeAuthenticity_ONT(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(ONTPromethION24, "YRI")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	score, err := f.AnalyzeAuthenticity(data)
	if err != nil {
		t.Fatalf("AnalyzeAuthenticity(ONT) error: %v", err)
	}
	if score < 0.0 || score > 100.0 {
		t.Errorf("AnalyzeAuthenticity score = %f, want in [0,100]", score)
	}
}

func TestAnalyzeAuthenticity_NilData(t *testing.T) {
	f := NewE2EFramework()
	_, err := f.AnalyzeAuthenticity(nil)
	if err == nil {
		t.Error("AnalyzeAuthenticity(nil) should return error")
	}
}

// Test RunBenchmark with different platforms to cover addPlatformSpecificBenchmarks
func TestRunBenchmark_ONT(t *testing.T) {
	f := NewE2EFramework()
	benchmark := f.RunBenchmark(ONTPromethION24)
	if benchmark == nil {
		t.Error("RunBenchmark(ONT) returned nil")
	}
}

func TestRunBenchmark_IlluminaNextSeq(t *testing.T) {
	f := NewE2EFramework()
	benchmark := f.RunBenchmark(IlluminaNextSeq2000)
	if benchmark == nil {
		t.Error("RunBenchmark(IlluminaNextSeq) returned nil")
	}
}

// Test RunWorkflow with different workflow engines
func TestRunWorkflow_CromwellWDL(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	result, err := f.RunWorkflow(CromwellWDL, IlluminaNovaSeq6000, data)
	if err != nil {
		t.Fatalf("RunWorkflow(Cromwell) error: %v", err)
	}
	if result == nil {
		t.Error("RunWorkflow returned nil")
	}
}

func TestRunWorkflow_Nextflow(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, _ := f.GenerateSyntheticData(config)
	result, err := f.RunWorkflow(NextflowNFL, IlluminaNovaSeq6000, data)
	if err != nil {
		t.Fatalf("RunWorkflow(Nextflow) error: %v", err)
	}
	if result.SuccessRate <= 0 {
		t.Errorf("RunWorkflow(Nextflow) SuccessRate = %f, want > 0", result.SuccessRate)
	}
}

func TestRunWorkflow_Snakemake(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, _ := f.GenerateSyntheticData(config)
	result, err := f.RunWorkflow(SnakemakeRule, IlluminaNovaSeq6000, data)
	if err != nil {
		t.Fatalf("RunWorkflow(Snakemake) error: %v", err)
	}
	if result == nil {
		t.Error("RunWorkflow(Snakemake) returned nil")
	}
}

// Test file writing functions
func TestWriteFASTQ(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}

	r1 := t.TempDir() + "/test_R1.fastq"
	r2 := t.TempDir() + "/test_R2.fastq"
	if err := f.WriteFASTQ(data, r1, r2); err != nil {
		t.Errorf("WriteFASTQ error: %v", err)
	}
}

func TestWriteReferenceGenome(t *testing.T) {
	f := NewE2EFramework()
	fname := t.TempDir() + "/ref.fa"
	if err := f.WriteReferenceGenome(fname); err != nil {
		t.Errorf("WriteReferenceGenome error: %v", err)
	}
}

func TestWriteVCF(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	fname := t.TempDir() + "/test.vcf"
	if err := f.WriteVCF(data, fname); err != nil {
		t.Errorf("WriteVCF error: %v", err)
	}
}

func TestWriteBAM(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	fname := t.TempDir() + "/test.bam"
	err = f.WriteBAM(data, fname)
	if _, lookErr := exec.LookPath("samtools"); lookErr != nil {
		if !errors.Is(err, ErrBAMUnsupported) {
			t.Fatalf("WriteBAM without samtools: err = %v, want ErrBAMUnsupported", err)
		}
		t.Skip("samtools not installed; BAM conversion not exercised")
	}
	if err != nil {
		t.Fatalf("WriteBAM error: %v", err)
	}

	// The BAM must hold one record per read; the old header stub held none.
	out, err := exec.Command("samtools", "view", "-c", fname).Output()
	if err != nil {
		t.Fatalf("samtools view -c: %v", err)
	}
	if got := strings.TrimSpace(string(out)); got != strconv.Itoa(len(data.Sequences)) {
		t.Errorf("BAM record count = %s, want %d", got, len(data.Sequences))
	}
	if err := exec.Command("samtools", "quickcheck", fname).Run(); err != nil {
		t.Errorf("samtools quickcheck failed: %v", err)
	}
}

// TestWriteBAM_NoSamtools checks that WriteBAM refuses rather than writing a
// fake file when samtools is missing.
func TestWriteBAM_NoSamtools(t *testing.T) {
	orig := lookPath
	t.Cleanup(func() { lookPath = orig })
	lookPath = func(string) (string, error) { return "", exec.ErrNotFound }

	f := NewE2EFramework()
	data, err := f.GenerateSyntheticData(f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU"))
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	fname := t.TempDir() + "/test.bam"
	if err := f.WriteBAM(data, fname); !errors.Is(err, ErrBAMUnsupported) {
		t.Fatalf("WriteBAM err = %v, want ErrBAMUnsupported", err)
	}
	if _, err := os.Stat(fname); !os.IsNotExist(err) {
		t.Errorf("WriteBAM wrote %s despite missing samtools", fname)
	}
}

// TestWriteSAM checks the header and that every read appears as an unmapped
// record with its sequence and qualities.
func TestWriteSAM(t *testing.T) {
	f := NewE2EFramework()
	data, err := f.GenerateSyntheticData(f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU"))
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	fname := t.TempDir() + "/test.sam"
	if err := f.WriteSAM(data, fname); err != nil {
		t.Fatalf("WriteSAM error: %v", err)
	}
	content, err := os.ReadFile(fname)
	if err != nil {
		t.Fatal(err)
	}

	refLen := 0
	for _, s := range data.Sequences {
		refLen += len(s)
	}
	var header, records []string
	for _, line := range strings.Split(strings.TrimRight(string(content), "\n"), "\n") {
		if strings.HasPrefix(line, "@") {
			header = append(header, line)
		} else {
			records = append(records, line)
		}
	}
	if len(header) < 2 || !strings.HasPrefix(header[0], "@HD\tVN:") {
		t.Fatalf("missing @HD header: %q", header)
	}
	if want := "@SQ\tSN:chr1\tLN:" + strconv.Itoa(refLen); header[1] != want {
		t.Errorf("@SQ line = %q, want %q", header[1], want)
	}
	if want := "@RG\tID:" + samReadGroupID + "\tSM:CEU\tPL:ILLUMINA"; len(header) < 3 || header[2] != want {
		t.Errorf("@RG line = %q, want %q", header[2:], want)
	}
	if len(records) != len(data.Sequences) {
		t.Fatalf("records = %d, want one per read (%d)", len(records), len(data.Sequences))
	}
	for i, rec := range records {
		fields := strings.Split(rec, "\t")
		if len(fields) != 12 {
			t.Fatalf("record %d has %d fields, want 11 mandatory plus RG: %q", i, len(fields), rec)
		}
		if fields[11] != "RG:Z:"+samReadGroupID {
			t.Errorf("record %d tag = %q, want RG:Z:%s", i, fields[11], samReadGroupID)
		}
		if fields[0] != strings.TrimPrefix(data.Headers[i], "@") {
			t.Errorf("record %d QNAME = %q, want %q", i, fields[0], data.Headers[i])
		}
		if got := strings.Join(fields[1:9], " "); got != "4 * 0 0 * * 0 0" {
			t.Errorf("record %d fields 2-9 = %q, want unmapped \"4 * 0 0 * * 0 0\"", i, got)
		}
		if fields[9] != string(data.Sequences[i]) || fields[10] != string(data.QualityScores[i]) {
			t.Errorf("record %d SEQ/QUAL do not match read %d", i, i)
		}
	}

	if err := f.WriteSAM(&SyntheticData{}, fname); err == nil {
		t.Error("WriteSAM with no reads should fail")
	}
}

// Test GeneratePerformanceReport (prints to stdout)
func TestGeneratePerformanceReport_NoError(t *testing.T) {
	f := NewE2EFramework()
	benchmarks := map[SequencingPlatform]*PlatformBenchmark{
		IlluminaNovaSeq6000: f.RunBenchmark(IlluminaNovaSeq6000),
		PacBioRevio:         f.RunBenchmark(PacBioRevio),
	}
	// Should not panic
	f.GeneratePerformanceReport(benchmarks)
}

// Test GenerateBalancedDataset with different platforms
func TestGenerateBalancedDataset_PacBio(t *testing.T) {
	f := NewE2EFramework()
	realData, synData, err := f.GenerateBalancedDataset(PacBioRevio, 5)
	if err != nil {
		t.Fatalf("GenerateBalancedDataset(PacBio) error: %v", err)
	}
	if realData == nil || synData == nil {
		t.Error("GenerateBalancedDataset(PacBio) returned nil data")
	}
}

// Test analyzePopulationAuthenticity branches (via AnalyzeAuthenticity with populated data)
func TestAnalyzeAuthenticity_CEU(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	score, err := f.AnalyzeAuthenticity(data)
	if err != nil {
		t.Fatalf("AnalyzeAuthenticity(CEU) error: %v", err)
	}
	if score < 0.0 || score > 100.0 {
		t.Errorf("AnalyzeAuthenticity score = %f, want in [0,100]", score)
	}
}

func TestAnalyzeAuthenticity_CHB(t *testing.T) {
	f := NewE2EFramework()
	config := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CHB")
	data, err := f.GenerateSyntheticData(config)
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	score, err := f.AnalyzeAuthenticity(data)
	if err != nil {
		t.Fatalf("AnalyzeAuthenticity(CHB) error: %v", err)
	}
	if score < 0.0 || score > 100.0 {
		t.Errorf("AnalyzeAuthenticity score = %f, want in [0,100]", score)
	}
}

// TestGCSkewed checks the positive-control transform: GC content rises by
// about (1-g)*p, the input is untouched, and the same seed reproduces it.
func TestGCSkewed(t *testing.T) {
	f := NewE2EFramework()
	data, err := f.generateSyntheticDataN(f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU"), 500)
	if err != nil {
		t.Fatalf("generateSyntheticDataN: %v", err)
	}
	orig := make([][]byte, len(data.Sequences))
	for i, s := range data.Sequences {
		orig[i] = append([]byte(nil), s...)
	}

	const p = 0.4
	before := f.calculateGCContent(data.Sequences)
	skewed := GCSkewed(data, p, 7)
	after := f.calculateGCContent(skewed.Sequences)

	want := before + (1-before)*p
	if math.Abs(after-want) > 0.01 {
		t.Errorf("GC content after skew = %.4f, want %.4f ± 0.01 (before %.4f)", after, want, before)
	}
	for i := range orig {
		if !bytes.Equal(orig[i], data.Sequences[i]) {
			t.Fatalf("GCSkewed modified input sequence %d", i)
		}
		if len(skewed.Sequences[i]) != len(orig[i]) {
			t.Fatalf("sequence %d length changed: %d -> %d", i, len(orig[i]), len(skewed.Sequences[i]))
		}
	}
	again := GCSkewed(data, p, 7)
	for i := range skewed.Sequences {
		if !bytes.Equal(skewed.Sequences[i], again.Sequences[i]) {
			t.Fatalf("same seed produced different output at sequence %d", i)
		}
	}
	if same := GCSkewed(data, 0, 7); !bytes.Equal(same.Sequences[0], orig[0]) {
		t.Error("p=0 should leave sequences unchanged")
	}
}
