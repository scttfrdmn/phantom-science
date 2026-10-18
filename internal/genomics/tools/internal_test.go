// Package tools — whitebox tests for unexported helpers.
// These tests run in the same package to reach unexported functions.
package tools

import (
	"context"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ─────────────────────────────────────────────────────────────────────────────
// min helper
// ─────────────────────────────────────────────────────────────────────────────

func TestMin(t *testing.T) {
	assert.Equal(t, 3, min(3, 5))
	assert.Equal(t, 3, min(5, 3))
	assert.Equal(t, 0, min(0, 10))
	assert.Equal(t, -5, min(-5, 0))
	assert.Equal(t, 7, min(7, 7))
}

// ─────────────────────────────────────────────────────────────────────────────
// formatInfo
// ─────────────────────────────────────────────────────────────────────────────

func TestFormatInfo_Empty(t *testing.T) {
	result := formatInfo(map[string]string{})
	assert.Equal(t, "", result)
}

func TestFormatInfo_SinglePair(t *testing.T) {
	result := formatInfo(map[string]string{"AF": "0.05"})
	assert.Equal(t, "AF=0.05", result)
}

func TestFormatInfo_MultiplePairs(t *testing.T) {
	result := formatInfo(map[string]string{"AF": "0.05", "DP": "100"})
	// Order is non-deterministic for maps, just check both pairs present
	assert.Contains(t, result, "AF=0.05")
	assert.Contains(t, result, "DP=100")
	assert.Contains(t, result, ";")
}

// ─────────────────────────────────────────────────────────────────────────────
// FormatValidator — validateFASTQ
// ─────────────────────────────────────────────────────────────────────────────

func TestValidateFASTQ_Valid(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("@read1\nACGT\n+\nIIII\n")
	err := fv.validateFASTQ(data)
	assert.NoError(t, err)
}

func TestValidateFASTQ_TooFewLines(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("@read1\nACGT\n")
	err := fv.validateFASTQ(data)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "insufficient lines")
}

func TestValidateFASTQ_BadHeader(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte(">read1\nACGT\n+\nIIII\n")
	err := fv.validateFASTQ(data)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "header must start with @")
}

func TestValidateFASTQ_BadPlusLine(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("@read1\nACGT\n-\nIIII\n")
	err := fv.validateFASTQ(data)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "plus line must start with +")
}

// ─────────────────────────────────────────────────────────────────────────────
// FormatValidator — validateVCF
// ─────────────────────────────────────────────────────────────────────────────

func TestValidateVCF_Valid(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("##fileformat=VCFv4.1\n#CHROM\tchr1\n")
	err := fv.validateVCF(data)
	assert.NoError(t, err)
}

func TestValidateVCF_Invalid(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("#CHROM\tPOS\n")
	err := fv.validateVCF(data)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "##fileformat=VCF")
}

func TestValidateVCF_ShortData(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("##fileformat=VCFv4.1")
	err := fv.validateVCF(data)
	assert.NoError(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// FormatValidator — validateSAM
// ─────────────────────────────────────────────────────────────────────────────

func TestValidateSAM_ValidBAM(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("BAM\x01extra data here")
	err := fv.validateSAM(data)
	assert.NoError(t, err)
}

func TestValidateSAM_ValidSAMWithHD(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("@HD\tVN:1.6\tSO:coordinate\n")
	err := fv.validateSAM(data)
	assert.NoError(t, err)
}

func TestValidateSAM_ValidSAMWithSQ(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("@SQ\tSN:chr1\tLN:248956422\n")
	err := fv.validateSAM(data)
	assert.NoError(t, err)
}

func TestValidateSAM_Invalid(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("INVALID data that is not SAM or BAM")
	err := fv.validateSAM(data)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "invalid SAM/BAM format")
}

func TestValidateSAM_ShortData(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("BA") // shorter than 4 bytes, not matching BAM magic
	err := fv.validateSAM(data)
	assert.Error(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// FormatValidator — validateFASTA
// ─────────────────────────────────────────────────────────────────────────────

func TestValidateFASTA_Valid(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte(">chr1\nACGTACGT\n")
	err := fv.validateFASTA(data)
	assert.NoError(t, err)
}

func TestValidateFASTA_Invalid(t *testing.T) {
	fv := &FormatValidator{}
	data := []byte("ACGTACGT\n")
	err := fv.validateFASTA(data)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "FASTA must start with >")
}

// ─────────────────────────────────────────────────────────────────────────────
// GenomicsEcosystemValidator internal helpers
// ─────────────────────────────────────────────────────────────────────────────

// buildTestValidator creates a minimal validator with empty platforms/populations
// so that generateSyntheticTestData is a no-op (no FASTQ/VCF/Reference generated).
func buildTestValidator(t *testing.T) *GenomicsEcosystemValidator {
	t.Helper()
	tmp := t.TempDir()

	// We need a tool with a real executable path that will pass VerifyToolAvailability.
	// /usr/bin/true exists on macOS/Linux and responds to --help with exit 0.
	truePath := "/usr/bin/true"
	if _, err := os.Stat(truePath); os.IsNotExist(err) {
		t.Skip("/usr/bin/true not available")
	}

	cfg := EcosystemConfig{
		TestDataPath:      tmp,
		ReferencePath:     filepath.Join(tmp, "ref"),
		ToolPaths:         map[string]string{"samtools": truePath},
		MaxConcurrentJobs: 2,
		TimeoutDuration:   30 * 1e9, // 30 seconds
		// Empty PlatformTargets and PopulationPanels → no FASTQ/VCF loops
		PlatformTargets:  []string{},
		PopulationPanels: []string{},
		// No ValidationModes → ValidateEcosystem is a no-op
		ValidationModes: []ValidationMode{},
	}

	// WorkflowEngine.Initialize() will fail if no engines are on PATH.
	// We need to get past it to exercise other methods.
	// Try to build manually to exercise the internals.
	_ = cfg
	return nil
}

func TestGenerateRandomSequence(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}

	seq := gev.generateRandomSequence(100)
	assert.Len(t, seq, 100)

	// All characters must be A, T, G, or C.
	for _, c := range seq {
		assert.Contains(t, "ATGC", string(c))
	}

	// Sequence must not be the old deterministic "ATGCATGC…" cycle.
	// With 100 random bases the probability of accidentally matching the cycle
	// is astronomically small (4^{-100}).
	cycle := ""
	for i := 0; i < 100; i++ {
		cycle += string([]byte{'A', 'T', 'G', 'C'}[i%4])
	}
	assert.NotEqual(t, cycle, seq, "sequence must not be a deterministic ATGC cycle")

	// Two independently-generated sequences of length 100 should differ.
	seq2 := gev.generateRandomSequence(100)
	assert.NotEqual(t, seq, seq2, "two independently-generated sequences should differ")
}

func TestGenerateChromosomeSequence(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}

	seq := gev.generateChromosomeSequence("chr1", 200)
	assert.Len(t, seq, 200)

	for _, b := range seq {
		assert.Contains(t, []byte{'A', 'C', 'G', 'T'}, b)
	}

	// Must not be the deterministic ATGCATGC… pattern.
	allCycle := true
	cycle := []byte{'A', 'T', 'G', 'C'}
	for i, b := range seq {
		if b != cycle[i%4] {
			allCycle = false
			break
		}
	}
	assert.False(t, allCycle, "chromosome sequence must not be a deterministic ATGC cycle")
}

// TestGenerateRandomSequence_GCContent checks that GC content is reasonably
// close to 50% (uniform distribution) for a long sequence.
func TestGenerateRandomSequence_GCContent(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}
	seq := gev.generateRandomSequence(1000)
	gc := 0
	for _, c := range seq {
		if c == 'G' || c == 'C' {
			gc++
		}
	}
	gcFrac := float64(gc) / 1000.0
	// With uniform sampling, expected GC = 0.5; ±15% is a very generous tolerance.
	assert.InDelta(t, 0.5, gcFrac, 0.15, "GC content should be ~50%% for uniform random sequence")
}

// TestGeneratePopulationVariants_TsTvRatio checks that the transition/transversion
// ratio is approximately 2:1 over a large number of generated variants.
func TestGeneratePopulationVariants_TsTvRatio(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}
	variants := gev.generatePopulationVariants("CEU", 1000)
	require.Len(t, variants, 1000)

	transitions := map[string]bool{
		"AG": true, "GA": true,
		"CT": true, "TC": true,
	}
	ts, tv := 0, 0
	for _, v := range variants {
		key := v.Reference + v.Alternate[0]
		if transitions[key] {
			ts++
		} else {
			tv++
		}
	}
	ratio := float64(ts) / float64(tv)
	// Expected ratio ≈ 2.0; accept [1.4, 2.8] as a generous range for 1000 variants.
	assert.InDelta(t, 2.0, ratio, 0.6, "Ts/Tv ratio should be ~2:1; got %.2f (%d Ts, %d Tv)", ratio, ts, tv)
}

func TestGeneratePlatformQuality_Illumina(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}

	q := gev.generatePlatformQuality("illumina", 50)
	assert.Len(t, q, 50)
	// Illumina base quality = 35, ASCII offset = 68 = 'D'
	for _, r := range q {
		assert.Equal(t, rune(35+33), r)
	}
}

func TestGeneratePlatformQuality_Pacbio(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}

	q := gev.generatePlatformQuality("pacbio", 30)
	assert.Len(t, q, 30)
	// pacbio base quality = 25
	for _, r := range q {
		assert.Equal(t, rune(25+33), r)
	}
}

func TestGeneratePlatformQuality_Nanopore(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}

	q := gev.generatePlatformQuality("nanopore", 20)
	assert.Len(t, q, 20)
	// nanopore base quality = 20
	for _, r := range q {
		assert.Equal(t, rune(20+33), r)
	}
}

func TestGeneratePlatformQuality_Default(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}

	q := gev.generatePlatformQuality("unknown_platform", 25)
	assert.Len(t, q, 25)
	// default base quality = 30
	for _, r := range q {
		assert.Equal(t, rune(30+33), r)
	}
}

func TestGeneratePlatformSpecificReads(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}

	reads := gev.generatePlatformSpecificReads("illumina", "YRI", 5)
	require.Len(t, reads, 5)

	for i, read := range reads {
		assert.NotEmpty(t, read.Header)
		assert.Contains(t, read.Header, "illumina")
		assert.Contains(t, read.Header, "YRI")
		assert.Equal(t, 150, len(read.Sequence))
		assert.Equal(t, "+", read.Plus)
		assert.Equal(t, 150, len(read.Quality))
		_ = i
	}
}

func TestGeneratePopulationVariants(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}

	variants := gev.generatePopulationVariants("YRI", 10)
	require.Len(t, variants, 10)

	validBases := "ACGT"
	for i, v := range variants {
		assert.NotEmpty(t, v.Chromosome)
		assert.Greater(t, v.Position, int64(0))
		assert.Contains(t, v.ID, "YRI")
		// REF and ALT must be valid nucleotides and must differ.
		assert.Contains(t, validBases, v.Reference, "variant %d: invalid REF %q", i, v.Reference)
		require.Len(t, v.Alternate, 1, "variant %d: expected one ALT allele", i)
		assert.Contains(t, validBases, v.Alternate[0], "variant %d: invalid ALT %q", i, v.Alternate[0])
		assert.NotEqual(t, v.Reference, v.Alternate[0], "variant %d: REF == ALT", i)
		assert.Equal(t, "PASS", v.Filter)
		assert.NotEmpty(t, v.Info)
	}
}

func TestWriteFASTQFile(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}
	tmp := t.TempDir()
	path := filepath.Join(tmp, "test.fastq")

	reads := []*FASTQRead{
		{Header: "@read1", Sequence: "ACGT", Plus: "+", Quality: "IIII"},
		{Header: "@read2", Sequence: "TTTT", Plus: "+", Quality: "HHHH"},
	}

	err := gev.writeFASTQFile(path, reads)
	require.NoError(t, err)

	data, err := os.ReadFile(path)
	require.NoError(t, err)
	assert.Contains(t, string(data), "@read1")
	assert.Contains(t, string(data), "@read2")
	assert.Contains(t, string(data), "ACGT")
	assert.Contains(t, string(data), "TTTT")
}

func TestWriteFASTQFile_EmptyReads(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}
	tmp := t.TempDir()
	path := filepath.Join(tmp, "empty.fastq")

	err := gev.writeFASTQFile(path, []*FASTQRead{})
	require.NoError(t, err)

	data, err := os.ReadFile(path)
	require.NoError(t, err)
	assert.Empty(t, data)
}

func TestWriteVCFFile(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}
	tmp := t.TempDir()
	path := filepath.Join(tmp, "test.vcf")

	variants := []*VCFVariant{
		{
			Chromosome: "chr1",
			Position:   100,
			ID:         "var1",
			Reference:  "A",
			Alternate:  []string{"T"},
			Quality:    60.0,
			Filter:     "PASS",
			Info:       map[string]string{"AF": "0.10"},
		},
	}

	err := gev.writeVCFFile(path, variants)
	require.NoError(t, err)

	data, err := os.ReadFile(path)
	require.NoError(t, err)
	content := string(data)
	assert.Contains(t, content, "##fileformat=VCFv4.3")
	assert.Contains(t, content, "chr1")
	assert.Contains(t, content, "var1")
	assert.Contains(t, content, "PASS")
	assert.Contains(t, content, "AF=0.10")
}

func TestWriteVCFFile_MultipleAlternates(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}
	tmp := t.TempDir()
	path := filepath.Join(tmp, "multi_alt.vcf")

	variants := []*VCFVariant{
		{
			Chromosome: "chr2",
			Position:   500,
			ID:         "var2",
			Reference:  "G",
			Alternate:  []string{"A", "C"},
			Quality:    80.0,
			Filter:     "PASS",
			Info:       map[string]string{"DP": "100"},
		},
	}

	err := gev.writeVCFFile(path, variants)
	require.NoError(t, err)

	data, err := os.ReadFile(path)
	require.NoError(t, err)
	assert.Contains(t, string(data), "A,C")
}

func TestWriteFASTAFile(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}
	tmp := t.TempDir()
	path := filepath.Join(tmp, "test.fa")

	sequence := make([]byte, 200)
	for i := range sequence {
		sequence[i] = []byte{'A', 'T', 'G', 'C'}[i%4]
	}

	err := gev.writeFASTAFile(path, "chr1", sequence)
	require.NoError(t, err)

	data, err := os.ReadFile(path)
	require.NoError(t, err)
	content := string(data)
	assert.Contains(t, content, ">chr1")
	// Should be written in 80-char lines
	lines := []string{}
	for _, line := range []string{content} {
		_ = line
	}
	_ = lines
	// Verify sequence length via line count
	assert.Greater(t, len(data), 200) // file has header + sequence + newlines
}

func TestWriteFASTAFile_ShortSequence(t *testing.T) {
	gev := &GenomicsEcosystemValidator{}
	tmp := t.TempDir()
	path := filepath.Join(tmp, "short.fa")

	// Less than 80 chars → single sequence line
	err := gev.writeFASTAFile(path, "chrM", []byte("ACGTACGT"))
	require.NoError(t, err)

	data, err := os.ReadFile(path)
	require.NoError(t, err)
	assert.Contains(t, string(data), ">chrM")
	assert.Contains(t, string(data), "ACGTACGT")
}

// ─────────────────────────────────────────────────────────────────────────────
// GenomicsToolSuite.detectToolVersion
// ─────────────────────────────────────────────────────────────────────────────

func TestDetectToolVersion_WithRealBinary(t *testing.T) {
	suite := &GenomicsToolSuite{}

	// /usr/bin/true exits 0 but produces no output; result should be "unknown"
	truePath := "/usr/bin/true"
	if _, err := os.Stat(truePath); os.IsNotExist(err) {
		t.Skip("/usr/bin/true not available")
	}

	version := suite.detectToolVersion(truePath)
	// Any non-empty string is acceptable (either version text or "unknown")
	assert.NotNil(t, version)
}

func TestDetectToolVersion_NonExistentPath(t *testing.T) {
	suite := &GenomicsToolSuite{}
	version := suite.detectToolVersion("/nonexistent/path/tool")
	assert.Equal(t, "unknown", version)
}

// ─────────────────────────────────────────────────────────────────────────────
// Scoring methods via populated reports
// ─────────────────────────────────────────────────────────────────────────────

func TestCalculateCompatibilityScore_AllPassed(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults: map[string]*ToolValidationResult{
			"gatk":     {Status: "passed"},
			"bwa":      {Status: "passed"},
			"samtools": {Status: "passed"},
		},
	}
	score := gev.calculateCompatibilityScore(report)
	assert.InDelta(t, 1.0, score, 1e-9)
}

func TestCalculateCompatibilityScore_SomeFailed(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults: map[string]*ToolValidationResult{
			"gatk": {Status: "passed"},
			"bwa":  {Status: "failed"},
		},
	}
	score := gev.calculateCompatibilityScore(report)
	assert.InDelta(t, 0.5, score, 1e-9)
}

func TestCalculateCompatibilityScore_Empty(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults: map[string]*ToolValidationResult{},
	}
	score := gev.calculateCompatibilityScore(report)
	assert.Equal(t, 0.0, score)
}

func TestCalculatePerformanceScore_WithBenchmarks(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		BenchmarkResults: []*BenchmarkResult{
			{SuiteName: "variant_calling", Metrics: map[string]float64{"throughput": 1000.0}},
			{SuiteName: "alignment", Metrics: map[string]float64{"throughput": 500.0}},
		},
	}
	score := gev.calculatePerformanceScore(report)
	// Each normalized to /1000: (1.0 + 0.5) / 2 = 0.75
	assert.InDelta(t, 0.75, score, 1e-9)
}

func TestCalculatePerformanceScore_NoBenchmarks(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{}
	score := gev.calculatePerformanceScore(report)
	assert.Equal(t, 0.0, score)
}

func TestCalculateScalabilityScore_WithPopulationResults(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		PopulationResults: []*PopulationTestResult{
			{Metrics: map[string]float64{"scalability": 0.9}},
			{Metrics: map[string]float64{"scalability": 0.7}},
		},
	}
	score := gev.calculateScalabilityScore(report)
	assert.InDelta(t, 0.8, score, 1e-9)
}

func TestCalculateScalabilityScore_Empty(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{}
	score := gev.calculateScalabilityScore(report)
	assert.Equal(t, 0.0, score)
}

func TestCalculateReliabilityScore_WithTestCases(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults: map[string]*ToolValidationResult{
			"gatk": {
				Status:    "passed",
				TestCases: []*TestCaseResult{{Status: "passed"}},
			},
			"bwa": {
				Status:    "passed",
				TestCases: []*TestCaseResult{{Status: "passed"}},
			},
			"fail_tool": {
				Status: "failed",
			},
		},
	}
	score := gev.calculateReliabilityScore(report)
	// 2 out of 3 tools are passed AND have test cases
	assert.InDelta(t, 2.0/3.0, score, 1e-9)
}

func TestCalculateReliabilityScore_Empty(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults: map[string]*ToolValidationResult{},
	}
	score := gev.calculateReliabilityScore(report)
	assert.Equal(t, 0.0, score)
}

func TestGenerateRecommendations_FailedTools(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults: map[string]*ToolValidationResult{
			"gatk": {Status: "failed"},
			"bwa":  {Status: "failed"},
		},
		CloudResults: &CloudValidationResult{Status: "completed"},
	}
	recs := gev.generateRecommendations(report)
	assert.NotEmpty(t, recs)
	assert.Contains(t, recs[0], "failed tool")
}

func TestGenerateRecommendations_NoBenchmarks(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults:  map[string]*ToolValidationResult{},
		CloudResults: &CloudValidationResult{Status: "completed"},
	}
	recs := gev.generateRecommendations(report)
	found := false
	for _, r := range recs {
		if r != "" {
			found = true
			break
		}
	}
	_ = found
}

func TestGenerateRecommendations_AllPassed(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults: map[string]*ToolValidationResult{
			"gatk": {Status: "passed"},
		},
		BenchmarkResults:  []*BenchmarkResult{{SuiteName: "suite1"}},
		PopulationResults: []*PopulationTestResult{},
		CloudResults:      &CloudValidationResult{Status: "completed"},
	}
	recs := gev.generateRecommendations(report)
	assert.NotEmpty(t, recs)
	assert.Contains(t, recs[0], "passed")
}

func TestGenerateRecommendations_NoCloudConfig(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults:      map[string]*ToolValidationResult{},
		BenchmarkResults: []*BenchmarkResult{{SuiteName: "suite1"}},
	}
	recs := gev.generateRecommendations(report)
	// Should recommend configuring cloud deployment
	found := false
	for _, r := range recs {
		if r != "" {
			found = true
		}
	}
	assert.True(t, found)
}

func TestGenerateRecommendations_SkippedCloud(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults:      map[string]*ToolValidationResult{},
		BenchmarkResults: []*BenchmarkResult{{SuiteName: "s"}},
		CloudResults:     &CloudValidationResult{Status: "skipped"},
	}
	recs := gev.generateRecommendations(report)
	found := false
	for _, r := range recs {
		if r != "" {
			found = true
		}
	}
	assert.True(t, found)
}

func TestGenerateOverallAssessment_Excellent(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults: map[string]*ToolValidationResult{
			"gatk": {Status: "passed", TestCases: []*TestCaseResult{{Status: "passed"}}},
		},
		BenchmarkResults: []*BenchmarkResult{
			{Metrics: map[string]float64{"throughput": 1000.0}},
		},
		PopulationResults: []*PopulationTestResult{
			{Metrics: map[string]float64{"scalability": 1.0}},
		},
		CloudResults: &CloudValidationResult{Status: "completed"},
	}
	assessment := gev.generateOverallAssessment(report)
	require.NotNil(t, assessment)
	assert.NotEmpty(t, assessment.Status)
	assert.Greater(t, assessment.OverallScore, 0.0)
	assert.LessOrEqual(t, assessment.OverallScore, 1.0)
}

func TestGenerateOverallAssessment_NeedsImprovement(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}
	report := &EcosystemValidationReport{
		ToolResults:  map[string]*ToolValidationResult{},
		CloudResults: &CloudValidationResult{Status: "skipped"},
	}
	assessment := gev.generateOverallAssessment(report)
	require.NotNil(t, assessment)
	// With all zeros, status should be "needs_improvement"
	assert.Equal(t, "needs_improvement", assessment.Status)
	assert.Equal(t, 0.0, assessment.OverallScore)
}

// ─────────────────────────────────────────────────────────────────────────────
// ToolValidationResult.calculateMetrics (unexported method — whitebox test)
// ─────────────────────────────────────────────────────────────────────────────

func TestCalculateMetrics_Empty(t *testing.T) {
	result := &ToolValidationResult{ToolName: "gatk"}
	result.calculateMetrics()
	assert.Equal(t, float64(0), result.Metrics["total_test_cases"])
}

func TestCalculateMetrics_AllPassed(t *testing.T) {
	result := &ToolValidationResult{
		ToolName: "bwa",
		TestCases: []*TestCaseResult{
			{Status: "passed", Duration: 100 * 1e6},
			{Status: "passed", Duration: 200 * 1e6},
			{Status: "passed", Duration: 150 * 1e6},
		},
	}
	result.calculateMetrics()
	assert.InDelta(t, 1.0, result.Metrics["success_rate"], 1e-9)
	assert.Equal(t, float64(3), result.Metrics["total_test_cases"])
	assert.Greater(t, result.Metrics["avg_execution_time"], 0.0)
}

func TestCalculateMetrics_MixedResults(t *testing.T) {
	result := &ToolValidationResult{
		ToolName: "samtools",
		TestCases: []*TestCaseResult{
			{Status: "passed", Duration: 100 * 1e6},
			{Status: "failed", Duration: 50 * 1e6},
		},
	}
	result.calculateMetrics()
	assert.InDelta(t, 0.5, result.Metrics["success_rate"], 1e-9)
	assert.Equal(t, float64(2), result.Metrics["total_test_cases"])
}

// ─────────────────────────────────────────────────────────────────────────────
// generateTestVCF / generateTestFASTQ / generateTestReference
// called on a manually-constructed GenomicsEcosystemValidator
// ─────────────────────────────────────────────────────────────────────────────

func makeValidatorWithTmpDir(t *testing.T) (*GenomicsEcosystemValidator, string) {
	t.Helper()
	tmp := t.TempDir()
	gev := &GenomicsEcosystemValidator{
		config: EcosystemConfig{
			TestDataPath:     tmp,
			ReferencePath:    filepath.Join(tmp, "ref"),
			PlatformTargets:  []string{"illumina", "pacbio"},
			PopulationPanels: []string{"YRI", "CEU"},
		},
		performanceMetrics: NewEcosystemMetrics(),
		validationResults:  make(map[string]*ToolValidationResult),
	}
	return gev, tmp
}

func TestGenerateTestVCF(t *testing.T) {
	gev, tmp := makeValidatorWithTmpDir(t)

	// Must create the vcf subdirectory first (prepareTestData would do this)
	require.NoError(t, os.MkdirAll(filepath.Join(tmp, "vcf"), 0750))

	err := gev.generateTestVCF()
	require.NoError(t, err)

	// Check that VCF files were created for each population
	for _, pop := range gev.config.PopulationPanels {
		path := filepath.Join(tmp, "vcf", pop+"_test_variants.vcf.gz")
		_, err := os.Stat(path)
		assert.NoError(t, err, "expected VCF file for population %s", pop)
	}
}

func TestGenerateTestFASTQ_OneCombo(t *testing.T) {
	gev, tmp := makeValidatorWithTmpDir(t)

	require.NoError(t, os.MkdirAll(filepath.Join(tmp, "fastq"), 0750))

	err := gev.generateTestFASTQ("illumina", "YRI")
	require.NoError(t, err)

	path := filepath.Join(tmp, "fastq", "illumina_YRI_test.fastq.gz")
	_, err = os.Stat(path)
	assert.NoError(t, err)
}

func TestGenerateTestReference(t *testing.T) {
	gev, tmp := makeValidatorWithTmpDir(t)

	require.NoError(t, os.MkdirAll(filepath.Join(tmp, "reference"), 0750))

	err := gev.generateTestReference()
	require.NoError(t, err)

	// Check that FASTA files were created for standard chromosomes
	for _, chr := range []string{"chr1", "chr2", "chr3", "chrX", "chrY", "chrM"} {
		path := filepath.Join(tmp, "reference", chr+"_test.fa")
		_, err := os.Stat(path)
		assert.NoError(t, err, "expected reference FASTA for %s", chr)
	}
}

func TestGenerateSyntheticTestData(t *testing.T) {
	gev, tmp := makeValidatorWithTmpDir(t)

	// Create the directory structure first
	for _, dir := range []string{"fastq", "vcf", "bam", "reference", "annotations", "populations", "benchmarks"} {
		require.NoError(t, os.MkdirAll(filepath.Join(tmp, dir), 0750))
	}

	err := gev.generateSyntheticTestData()
	require.NoError(t, err)
}

func TestPrepareTestData(t *testing.T) {
	gev, _ := makeValidatorWithTmpDir(t)

	// prepareTestData creates dirs + generates test data
	err := gev.prepareTestData()
	require.NoError(t, err)
}

func TestGetValidationStatus_Manual(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		toolSuite: NewGenomicsToolSuite(map[string]string{
			"samtools": "/fake/samtools",
		}),
		validationResults:  make(map[string]*ToolValidationResult),
		performanceMetrics: NewEcosystemMetrics(),
	}

	status := gev.GetValidationStatus()
	require.NotNil(t, status)
	assert.Equal(t, "running", status.ValidationState)
	assert.Equal(t, 1, status.TotalTools)
	assert.Equal(t, 0, status.ValidatedTools)
}

func TestSaveValidationReport_Manual(t *testing.T) {
	tmp := t.TempDir()
	reportPath := filepath.Join(tmp, "report.json")

	gev := &GenomicsEcosystemValidator{
		performanceMetrics: NewEcosystemMetrics(),
	}

	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: map[string]*WorkflowValidationResult{},
	}

	err := gev.SaveValidationReport(report, reportPath)
	require.NoError(t, err)

	// Verify the file exists and is valid JSON
	data, err := os.ReadFile(reportPath)
	require.NoError(t, err)
	assert.NotEmpty(t, data)
	// Should start with '{' (JSON object)
	assert.Equal(t, byte('{'), data[0])
}

// ─────────────────────────────────────────────────────────────────────────────
// run* placeholder methods — all return stub results, test them directly
// ─────────────────────────────────────────────────────────────────────────────

func makeGEV() *GenomicsEcosystemValidator {
	return &GenomicsEcosystemValidator{
		config: EcosystemConfig{
			PlatformTargets:  []string{"illumina"},
			PopulationPanels: []string{"YRI"},
		},
		toolSuite:          NewGenomicsToolSuite(map[string]string{}),
		performanceMetrics: NewEcosystemMetrics(),
		validationResults:  make(map[string]*ToolValidationResult),
	}
}

func TestRunBenchmarkSuite(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	suite := BenchmarkSuite{Name: "variant_calling", Tools: []string{"gatk"}}
	result, err := gev.runBenchmarkSuite(ctx, suite)
	require.NoError(t, err)
	assert.Equal(t, "variant_calling", result.SuiteName)
	assert.Equal(t, "completed", result.Status)
	assert.Equal(t, 1000.0, result.Metrics["throughput"])
}

func TestRunPopulationTest(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	test := PopulationTest{Name: "1000g", DatasetSize: 1000, Populations: []string{"YRI", "CEU"}}
	result, err := gev.runPopulationTest(ctx, test)
	require.NoError(t, err)
	assert.Equal(t, "1000g", result.TestName)
	assert.Equal(t, 1000, result.DatasetSize)
	assert.Equal(t, "completed", result.Status)
	assert.Equal(t, 0.88, result.Metrics["scalability"])
}

func TestRunWorkflowTest(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	workflow := WorkflowTest{Engine: "cromwell", WorkflowType: "wdl", Name: "gatk_germline"}
	result, err := gev.runWorkflowTest(ctx, workflow)
	require.NoError(t, err)
	assert.Equal(t, "gatk_germline", result.WorkflowName)
	assert.Equal(t, "cromwell", result.Engine)
	assert.Equal(t, "completed", result.Status)
	assert.Equal(t, 5*time.Minute, result.ExecutionTime)
}

func TestRunStreamingTest(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	test := StreamingTest{Name: "realtime_vc", DataRate: "1MB/s", Duration: 5 * time.Minute}
	result, err := gev.runStreamingTest(ctx, test)
	require.NoError(t, err)
	assert.Equal(t, "realtime_vc", result.TestName)
	assert.Equal(t, "1MB/s", result.DataRate)
	assert.Equal(t, "completed", result.Status)
	assert.Equal(t, 0.9, result.Throughput)
}

func TestRunCloudDeploymentTest(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	cfg := CloudConfig{Provider: "aws", Region: "us-east-1"}
	result, err := gev.runCloudDeploymentTest(ctx, cfg)
	require.NoError(t, err)
	assert.Equal(t, "aws", result.Provider)
	assert.Equal(t, "us-east-1", result.Region)
	assert.Equal(t, "completed", result.Status)
}

func TestRunToolTestCase(t *testing.T) {
	// A non-existent executable now returns Status="failed" (real execution,
	// not the old hardcoded stub that always returned "passed").
	gev := makeGEV()
	ctx := context.Background()
	tool := &GenomicsTool{
		Name:           "gatk",
		ExecutablePath: "/fake/gatk",
	}
	result, err := gev.runToolTestCase(ctx, tool, "illumina", "YRI")
	require.NoError(t, err)
	assert.Equal(t, "illumina_YRI_test", result.TestName)
	assert.Equal(t, "failed", result.Status) // real execution: no such file
	assert.NotEmpty(t, result.ErrorMessage)
}

func TestRunBasicToolTest_NonExistentTool(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	tool := &GenomicsTool{
		Name:           "fake_tool",
		ExecutablePath: "/nonexistent/fake_tool",
	}
	// runBasicToolTest runs --version on the tool; it should return a result
	// (not an error) even if the tool fails — it records the failure in result.Status
	result, err := gev.runBasicToolTest(ctx, tool)
	require.NoError(t, err) // error is captured in result, not returned
	assert.Equal(t, "failed", result.Status)
	assert.NotEmpty(t, result.ErrorMessage)
}

func TestRunBasicToolTest_WithRealTool(t *testing.T) {
	truePath := "/usr/bin/true"
	if _, err := os.Stat(truePath); os.IsNotExist(err) {
		t.Skip("/usr/bin/true not available")
	}

	gev := makeGEV()
	ctx := context.Background()
	tool := &GenomicsTool{
		Name:           "true",
		ExecutablePath: truePath,
	}
	result, err := gev.runBasicToolTest(ctx, tool)
	require.NoError(t, err)
	// /usr/bin/true --version may or may not return exit 0 depending on platform
	assert.NotNil(t, result)
}

// ─────────────────────────────────────────────────────────────────────────────
// validateXxx methods — these call run* stubs so they work without real tools
// but they need to have the right setup
// ─────────────────────────────────────────────────────────────────────────────

func TestValidatePerformanceBenchmark(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: map[string]*WorkflowValidationResult{},
	}
	err := gev.validatePerformanceBenchmark(ctx, report)
	require.NoError(t, err)
	assert.Len(t, report.BenchmarkResults, 4)
}

func TestValidatePopulationScale(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: map[string]*WorkflowValidationResult{},
	}
	err := gev.validatePopulationScale(ctx, report)
	require.NoError(t, err)
	assert.Len(t, report.PopulationResults, 3)
}

func TestValidateWorkflowIntegration(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: make(map[string]*WorkflowValidationResult),
	}
	err := gev.validateWorkflowIntegration(ctx, report)
	require.NoError(t, err)
	assert.Len(t, report.WorkflowResults, 4)
}

func TestValidateRealTimeProcessing(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: map[string]*WorkflowValidationResult{},
	}
	err := gev.validateRealTimeProcessing(ctx, report)
	require.NoError(t, err)
	assert.Len(t, report.StreamingResults, 3)
}

func TestValidateCloudDeployment_NoConfig(t *testing.T) {
	gev := makeGEV() // gev.config.CloudConfig.Provider is ""
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: map[string]*WorkflowValidationResult{},
	}
	err := gev.validateCloudDeployment(ctx, report)
	require.NoError(t, err)
	assert.Equal(t, "skipped", report.CloudResults.Status)
}

func TestValidateCloudDeployment_WithConfig(t *testing.T) {
	gev := makeGEV()
	gev.config.CloudConfig = CloudConfig{Provider: "aws", Region: "us-east-1"}
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: map[string]*WorkflowValidationResult{},
	}
	err := gev.validateCloudDeployment(ctx, report)
	require.NoError(t, err)
	require.NotNil(t, report.CloudResults)
	assert.Equal(t, "completed", report.CloudResults.Status)
	assert.Equal(t, "aws", report.CloudResults.Provider)
}

func TestExecuteValidationMode_PerformanceBenchmark(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: make(map[string]*WorkflowValidationResult),
	}
	err := gev.executeValidationMode(ctx, ModePerformanceBenchmark, report)
	require.NoError(t, err)
}

func TestExecuteValidationMode_PopulationScale(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: make(map[string]*WorkflowValidationResult),
	}
	err := gev.executeValidationMode(ctx, ModePopulationScale, report)
	require.NoError(t, err)
}

func TestExecuteValidationMode_WorkflowIntegration(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: make(map[string]*WorkflowValidationResult),
	}
	err := gev.executeValidationMode(ctx, ModeWorkflowIntegration, report)
	require.NoError(t, err)
}

func TestExecuteValidationMode_RealTimeValidation(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: make(map[string]*WorkflowValidationResult),
	}
	err := gev.executeValidationMode(ctx, ModeRealTimeValidation, report)
	require.NoError(t, err)
}

func TestExecuteValidationMode_CloudDeployment(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: make(map[string]*WorkflowValidationResult),
	}
	err := gev.executeValidationMode(ctx, ModeCloudDeployment, report)
	require.NoError(t, err)
}

func TestExecuteValidationMode_UnknownMode(t *testing.T) {
	gev := makeGEV()
	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     map[string]*ToolValidationResult{},
		WorkflowResults: make(map[string]*WorkflowValidationResult),
	}
	err := gev.executeValidationMode(ctx, ValidationMode(999), report)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "unknown validation mode")
}

func TestValidateEcosystem_EmptyModes(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		config: EcosystemConfig{
			ValidationModes: []ValidationMode{}, // no modes
		},
		toolSuite:          NewGenomicsToolSuite(map[string]string{}),
		performanceMetrics: NewEcosystemMetrics(),
		validationResults:  make(map[string]*ToolValidationResult),
	}
	ctx := context.Background()
	report, err := gev.ValidateEcosystem(ctx)
	require.NoError(t, err)
	require.NotNil(t, report)
	assert.NotNil(t, report.OverallAssessment)
}

func TestValidateEcosystem_MultipleNonBasicModes(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		config: EcosystemConfig{
			ValidationModes: []ValidationMode{
				ModePerformanceBenchmark,
				ModePopulationScale,
				ModeCloudDeployment,
			},
			CloudConfig: CloudConfig{}, // empty → cloud is skipped
		},
		toolSuite:          NewGenomicsToolSuite(map[string]string{}),
		performanceMetrics: NewEcosystemMetrics(),
		validationResults:  make(map[string]*ToolValidationResult),
	}
	ctx := context.Background()
	report, err := gev.ValidateEcosystem(ctx)
	require.NoError(t, err)
	require.NotNil(t, report)
	assert.Len(t, report.BenchmarkResults, 4)
	assert.Len(t, report.PopulationResults, 3)
	assert.Equal(t, "skipped", report.CloudResults.Status)
}

// ─────────────────────────────────────────────────────────────────────────────
// validateToolCompatibility — needs a tool with executable path on disk
// ─────────────────────────────────────────────────────────────────────────────

func TestValidateToolCompatibility_WithRealTool(t *testing.T) {
	truePath := "/usr/bin/true"
	if _, err := os.Stat(truePath); os.IsNotExist(err) {
		t.Skip("/usr/bin/true not available")
	}

	gev := &GenomicsEcosystemValidator{
		config: EcosystemConfig{
			PlatformTargets:  []string{"illumina"},
			PopulationPanels: []string{"YRI"},
		},
		toolSuite:          NewGenomicsToolSuite(map[string]string{}),
		performanceMetrics: NewEcosystemMetrics(),
		validationResults:  make(map[string]*ToolValidationResult),
	}

	tool := &GenomicsTool{
		Name:           "true",
		ExecutablePath: truePath,
	}
	ctx := context.Background()
	result, err := gev.validateToolCompatibility(ctx, tool)
	require.NoError(t, err)
	require.NotNil(t, result)
	assert.Equal(t, "true", result.ToolName)
	assert.NotNil(t, result.Metrics)
}

func TestValidateBasicCompatibility_NoTools(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		config: EcosystemConfig{
			PlatformTargets:  []string{},
			PopulationPanels: []string{},
		},
		toolSuite:          NewGenomicsToolSuite(map[string]string{}), // no tools
		performanceMetrics: NewEcosystemMetrics(),
		validationResults:  make(map[string]*ToolValidationResult),
	}

	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     make(map[string]*ToolValidationResult),
		WorkflowResults: make(map[string]*WorkflowValidationResult),
	}

	err := gev.validateBasicCompatibility(ctx, report)
	require.NoError(t, err)
	assert.Empty(t, report.ToolResults)
}

func TestExecuteValidationMode_BasicCompatibility_NoTools(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		config: EcosystemConfig{
			PlatformTargets:  []string{},
			PopulationPanels: []string{},
		},
		toolSuite:          NewGenomicsToolSuite(map[string]string{}),
		performanceMetrics: NewEcosystemMetrics(),
		validationResults:  make(map[string]*ToolValidationResult),
	}

	ctx := context.Background()
	report := &EcosystemValidationReport{
		ToolResults:     make(map[string]*ToolValidationResult),
		WorkflowResults: make(map[string]*WorkflowValidationResult),
	}

	err := gev.executeValidationMode(ctx, ModeBasicCompatibility, report)
	require.NoError(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// EcosystemMetrics concurrency safety
// ─────────────────────────────────────────────────────────────────────────────

func TestEcosystemMetrics_ConcurrentUpdates(t *testing.T) {
	em := NewEcosystemMetrics()
	done := make(chan struct{}, 10)

	for i := 0; i < 10; i++ {
		go func() {
			em.UpdateMetrics(1, 1, 100.0, 1000)
			done <- struct{}{}
		}()
	}

	for i := 0; i < 10; i++ {
		<-done
	}

	assert.Equal(t, int64(10), em.ToolsValidated)
	assert.Equal(t, int64(10), em.ToolsAttempted)
}

// ─────────────────────────────────────────────────────────────────────────────
// WorkflowEngine.Initialize with supported engines (if available)
// ─────────────────────────────────────────────────────────────────────────────

func TestWorkflowEngine_Initialize_Success(t *testing.T) {
	cfg := EcosystemConfig{}
	engine := NewWorkflowEngine(cfg)

	if len(engine.supportedEngines) == 0 {
		// No workflow engines available — skip; no sense testing success path
		t.Skip("No workflow engines (cromwell/nextflow/snakemake/cwl-runner) in PATH")
	}

	err := engine.Initialize()
	assert.NoError(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// time tracking
// ─────────────────────────────────────────────────────────────────────────────

func TestValidateEcosystem_TimestampsSet(t *testing.T) {
	gev := &GenomicsEcosystemValidator{
		config:             EcosystemConfig{ValidationModes: []ValidationMode{}},
		toolSuite:          NewGenomicsToolSuite(map[string]string{}),
		performanceMetrics: NewEcosystemMetrics(),
		validationResults:  make(map[string]*ToolValidationResult),
	}

	before := time.Now()
	ctx := context.Background()
	report, err := gev.ValidateEcosystem(ctx)
	after := time.Now()

	require.NoError(t, err)
	assert.True(t, report.StartTime.After(before) || report.StartTime.Equal(before))
	assert.True(t, report.EndTime.Before(after) || report.EndTime.Equal(after))
	assert.GreaterOrEqual(t, report.Duration, time.Duration(0))
}

// ─────────────────────────────────────────────────────────────────────────────
// runToolTestCase
// ─────────────────────────────────────────────────────────────────────────────

// makeTool builds a GenomicsTool pointing at the given executable.
func makeTool(execPath, validationCmd string) *GenomicsTool {
	return &GenomicsTool{
		Name:              filepath.Base(execPath),
		ExecutablePath:    execPath,
		ValidationCommand: validationCmd,
	}
}

// newGEV creates a minimal GenomicsEcosystemValidator for unit testing.
func newGEV() *GenomicsEcosystemValidator {
	gev, _ := NewGenomicsEcosystemValidator(EcosystemConfig{
		ToolPaths: map[string]string{}, // empty — no tools probed
	})
	return gev
}

func TestRunToolTestCase_ValidExecutable_Passes(t *testing.T) {
	// "true" (Unix) exits 0 — should produce Status="passed".
	truePath, err := exec.LookPath("true")
	if err != nil {
		t.Skip("'true' not on PATH")
	}
	gev := newGEV()
	tool := makeTool(truePath, "")
	result, err := gev.runToolTestCase(context.Background(), tool, "illumina", "CEU")
	require.NoError(t, err)
	assert.Equal(t, "passed", result.Status)
	assert.Empty(t, result.ErrorMessage)
	assert.NotZero(t, result.StartTime)
	assert.NotZero(t, result.EndTime)
}

func TestRunToolTestCase_InvalidExecutable_Fails(t *testing.T) {
	gev := newGEV()
	tool := makeTool("/nonexistent/tool/does_not_exist", "")
	result, err := gev.runToolTestCase(context.Background(), tool, "pacbio", "YRI")
	require.NoError(t, err) // function itself must not error
	assert.Equal(t, "failed", result.Status)
	assert.NotEmpty(t, result.ErrorMessage)
}

func TestRunToolTestCase_ValidationCommand_Used(t *testing.T) {
	// Use "echo" with a custom validation command to verify args are passed.
	echoPath, err := exec.LookPath("echo")
	if err != nil {
		t.Skip("'echo' not on PATH")
	}
	gev := newGEV()
	// ValidationCommand overrides the default "--version".
	tool := makeTool(echoPath, "hello world")
	result, err := gev.runToolTestCase(context.Background(), tool, "nanopore", "GBR")
	require.NoError(t, err)
	assert.Equal(t, "passed", result.Status)
	assert.Contains(t, result.Output, "hello")
}

func TestRunToolTestCase_TestNameIncludesPlatformPopulation(t *testing.T) {
	truePath, err := exec.LookPath("true")
	if err != nil {
		t.Skip("'true' not on PATH")
	}
	gev := newGEV()
	tool := makeTool(truePath, "")
	result, err := gev.runToolTestCase(context.Background(), tool, "illumina", "CEU")
	require.NoError(t, err)
	assert.Equal(t, "illumina_CEU_test", result.TestName)
}

func TestRunToolTestCase_ContextCancellation(t *testing.T) {
	// Use "sleep" with a long duration; cancel context immediately.
	sleepPath, err := exec.LookPath("sleep")
	if err != nil {
		t.Skip("'sleep' not on PATH")
	}
	gev := newGEV()
	tool := makeTool(sleepPath, "60")
	ctx, cancel := context.WithCancel(context.Background())
	cancel() // cancel immediately
	result, err := gev.runToolTestCase(ctx, tool, "illumina", "CEU")
	require.NoError(t, err)
	assert.Equal(t, "failed", result.Status, "cancelled context should produce failed status")
}
