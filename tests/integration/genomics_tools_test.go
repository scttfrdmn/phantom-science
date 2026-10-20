//go:build !release

package integration

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	pt "github.com/scttfrdmn/phantom-science/internal/testing"
)

// TestGATKWorkflowIntegration tests integration with GATK (Genome Analysis Toolkit)
func TestGATKWorkflowIntegration(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping GATK integration tests in short mode")
	}

	// The pipeline runs bwa and samtools as well as GATK.
	requireTools(t, "gatk", "bwa", "samtools")

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	platforms := []pt.SequencingPlatform{
		pt.IlluminaNovaSeq6000,
		pt.PacBioRevio,
		pt.ONTPromethION24,
	}

	populations := []string{"YRI", "CEU", "CHB"}

	for _, platform := range platforms {
		t.Run(string(platform), func(t *testing.T) {
			for _, population := range populations {
				t.Run(fmt.Sprintf("Population_%s", population), func(t *testing.T) {
					testGATKPipeline(t, framework, platform, population)
				})
			}
		})
	}
}

// testGATKPipeline tests the complete GATK variant calling pipeline
func testGATKPipeline(t *testing.T, framework *pt.E2EFramework, platform pt.SequencingPlatform, population string) {
	t.Logf("🔬 Testing GATK pipeline with %s platform, %s population", platform, population)

	// Step 1: Generate synthetic genomic data
	syntheticData := generateOrSkip(t, framework, platform, population)

	// Step 2: Create temporary working directory
	workDir, err := os.MkdirTemp("", fmt.Sprintf("gatk_test_%s_%s", platform, population))
	require.NoError(t, err)
	defer os.RemoveAll(workDir)

	// Step 3: Generate FASTQ files from synthetic data
	fastqR1 := filepath.Join(workDir, "synthetic_R1.fastq")
	fastqR2 := filepath.Join(workDir, "synthetic_R2.fastq")
	err = framework.WriteFASTQ(syntheticData, fastqR1, fastqR2)
	require.NoError(t, err)

	// Step 4: Build reference genome from the generated sequences so BWA can
	// map the reads back, validating that the FASTQ output is well-formed.
	refGenome := filepath.Join(workDir, "reference.fasta")
	err = framework.WriteReferenceFromData(syntheticData, refGenome)
	require.NoError(t, err)

	// Step 5: Run GATK preprocessing pipeline
	gatkResults := runGATKPreprocessing(t, workDir, fastqR1, fastqR2, refGenome)

	// Step 6: Validate GATK pipeline results
	validateGATKResults(t, gatkResults, platform, population)

	// Step 7: Run variant calling
	variantResults := runGATKVariantCalling(t, workDir, gatkResults.AlignedBAM, refGenome)

	// Step 8: Validate variant calling results
	validateVariantCallingResults(t, variantResults)

	t.Logf("✅ GATK pipeline completed successfully for %s %s", platform, population)
}

// TestBWAAlignmentIntegration tests integration with BWA aligner
func TestBWAAlignmentIntegration(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping BWA integration tests in short mode")
	}

	requireTools(t, "samtools")

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	// Test BWA-MEM (short reads) and BWA-MEM2 (long reads)
	alignmentTests := []struct {
		platform pt.SequencingPlatform
		aligner  string
		binary   string
	}{
		{pt.IlluminaNovaSeq6000, "bwa-mem", "bwa"},
		{pt.IlluminaNextSeq2000, "bwa-mem", "bwa"},
		{pt.PacBioRevio, "bwa-mem2", "bwa-mem2"},
		{pt.ONTPromethION24, "minimap2", "minimap2"}, // Better for very long reads
	}

	for _, test := range alignmentTests {
		t.Run(fmt.Sprintf("%s_%s", test.platform, test.aligner), func(t *testing.T) {
			requireTools(t, test.binary)
			testBWAAlignment(t, framework, test.platform, test.aligner)
		})
	}
}

// testBWAAlignment tests BWA alignment with synthetic data
func testBWAAlignment(t *testing.T, framework *pt.E2EFramework, platform pt.SequencingPlatform, aligner string) {
	t.Logf("🧬 Testing %s alignment with %s platform", aligner, platform)

	// Generate synthetic data
	syntheticData := generateOrSkip(t, framework, platform, "CEU")

	// Create working directory
	workDir, err := os.MkdirTemp("", fmt.Sprintf("bwa_test_%s", platform))
	require.NoError(t, err)
	defer os.RemoveAll(workDir)

	// Write FASTQ and reference files
	fastqR1 := filepath.Join(workDir, "reads_R1.fastq")
	fastqR2 := filepath.Join(workDir, "reads_R2.fastq")
	refGenome := filepath.Join(workDir, "reference.fasta")

	err = framework.WriteFASTQ(syntheticData, fastqR1, fastqR2)
	require.NoError(t, err)
	// Build reference from generated sequences so reads map back with high efficiency.
	err = framework.WriteReferenceFromData(syntheticData, refGenome)
	require.NoError(t, err)

	// Run alignment
	alignmentResult := runBWAAlignment(t, workDir, fastqR1, fastqR2, refGenome, aligner)

	// Validate alignment results
	validateAlignmentResults(t, alignmentResult, syntheticData, platform)

	t.Logf("✅ %s alignment completed: %.2f%% mapped, %.3f MAPQ",
		aligner, alignmentResult.MappingRate*100, alignmentResult.AverageMAPQ)
}

// TestBCFtoolsIntegration tests integration with BCFtools for variant processing
func TestBCFtoolsIntegration(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping BCFtools integration tests in short mode")
	}

	requireTools(t, "bcftools")

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	t.Run("Variant Processing Workflow", func(t *testing.T) {
		testBCFtoolsVariantProcessing(t, framework)
	})
}

// testBCFtoolsVariantProcessing tests BCFtools variant processing capabilities
func testBCFtoolsVariantProcessing(t *testing.T, framework *pt.E2EFramework) {
	t.Log("🧪 Testing BCFtools variant processing workflow")

	// Two call sets against one reference: WriteVCF's sites come from the
	// reference built from the reads, so both files share their sites.
	syntheticData := generateOrSkip(t, framework, pt.IlluminaNovaSeq6000, "CEU")

	workDir, err := os.MkdirTemp("", "bcftools_test")
	require.NoError(t, err)
	defer os.RemoveAll(workDir)

	var vcfFiles []string
	var written []int
	for i := 1; i <= 2; i++ {
		vcfFile := filepath.Join(workDir, fmt.Sprintf("sample_%d.vcf", i))
		require.NoError(t, framework.WriteVCF(syntheticData, vcfFile))
		vcfFiles = append(vcfFiles, vcfFile)
		stats, err := parseVCFStatistics(vcfFile)
		require.NoError(t, err)
		written = append(written, stats.VariantCount)
	}

	// Run BCFtools operations
	bcftoolsResult := runBCFtoolsOperations(t, workDir, vcfFiles)

	// Validate BCFtools results against what was written
	validateBCFtoolsResults(t, bcftoolsResult, written)

	t.Log("✅ BCFtools variant processing completed successfully")
}

// TestSamtoolsIntegration tests integration with Samtools
func TestSamtoolsIntegration(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping Samtools integration tests in short mode")
	}

	requireTools(t, "samtools")

	framework := pt.NewE2EFramework()

	t.Run("BAM Processing", func(t *testing.T) {
		testSamtoolsBAMProcessing(t, framework)
	})

	t.Run("Statistics Generation", func(t *testing.T) {
		testSamtoolsStatistics(t, framework)
	})
}

// TestPicardIntegration tests integration with Picard tools
func TestPicardIntegration(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping Picard integration tests in short mode")
	}

	requireTools(t, "picard")

	framework := pt.NewE2EFramework()

	t.Run("Quality Metrics", func(t *testing.T) {
		testPicardQualityMetrics(t, framework)
	})
}

// Data structures for results

type GATKResults struct {
	AlignedBAM     string
	ProcessingTime time.Duration
	MappingRate    float64
}

type VariantCallingResults struct {
	VCFFile        string
	VariantCount   int
	SNPCount       int
	IndelCount     int
	QualityScore   float64
	FilterPassRate float64
}

type AlignmentResults struct {
	BAMFile     string
	MappingRate float64
	AverageMAPQ float64
}

type BCFtoolsResults struct {
	MergedVCF      string
	VariantStats   VariantStatistics
	MergedVCFStats VCFStatistics
}

type VariantStatistics struct {
	TotalVariants     int
	SNPCount          int
	IndelCount        int
	TransitionCount   int
	TransversionCount int
}

// Implementation functions

// generateOrSkip generates reads for platform, skipping the test when the
// generator has no read model for it (pt.ErrPlatformNotModeled).
func generateOrSkip(t *testing.T, framework *pt.E2EFramework, platform pt.SequencingPlatform, population string) *pt.SyntheticData {
	t.Helper()
	data, err := framework.GenerateSyntheticData(framework.CreatePlatformConfig(platform, population))
	if errors.Is(err, pt.ErrPlatformNotModeled) {
		t.Skipf("%s: %v", platform, err)
	}
	require.NoError(t, err)
	return data
}

func isToolAvailable(tool string) bool {
	_, err := exec.LookPath(tool)
	return err == nil
}

// requireTools skips the test unless every tool is on PATH. Once a tool is
// present, its failures fail the test; results are never substituted.
func requireTools(t *testing.T, tools ...string) {
	t.Helper()
	for _, tool := range tools {
		if !isToolAvailable(tool) {
			t.Skipf("%s not available - skipping", tool)
		}
	}
}

// runTool runs name with args in dir and returns its stdout, failing the
// test with the tool's stderr when it exits non-zero.
func runTool(t *testing.T, dir, name string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	var stderr strings.Builder
	cmd.Stderr = &stderr
	out, err := cmd.Output()
	require.NoError(t, err, "%s %s: %s", name, strings.Join(args, " "), tailLines(stderr.String(), 20))
	return out
}

func runGATKPreprocessing(t *testing.T, workDir, fastqR1, fastqR2, refGenome string) *GATKResults {
	start := time.Now()

	runTool(t, workDir, "samtools", "faidx", refGenome)
	runTool(t, workDir, "bwa", "index", refGenome)

	// Align with BWA; HaplotypeCaller needs a read group.
	samFile := filepath.Join(workDir, "aligned.sam")
	sam := runTool(t, workDir, "bwa", "mem", "-R", `@RG\tID:synthetic\tSM:synthetic\tPL:ILLUMINA`, refGenome, fastqR1, fastqR2)
	require.NoError(t, os.WriteFile(samFile, sam, 0o600), "Failed to write SAM file")

	bamFile := filepath.Join(workDir, "aligned.bam")
	runTool(t, workDir, "samtools", "view", "-bS", samFile, "-o", bamFile)
	sortedBAM := filepath.Join(workDir, "aligned_sorted.bam")
	runTool(t, workDir, "samtools", "sort", bamFile, "-o", sortedBAM)
	runTool(t, workDir, "samtools", "index", sortedBAM)

	mappingRate, err := parseMappingRate(string(runTool(t, workDir, "samtools", "flagstat", sortedBAM)))
	require.NoError(t, err)

	return &GATKResults{
		AlignedBAM:     sortedBAM,
		ProcessingTime: time.Since(start),
		MappingRate:    mappingRate,
	}
}

func runGATKVariantCalling(t *testing.T, workDir, bamFile, refGenome string) *VariantCallingResults {
	dictFile := strings.TrimSuffix(refGenome, filepath.Ext(refGenome)) + ".dict"
	runTool(t, workDir, "gatk", "CreateSequenceDictionary", "-R", refGenome, "-O", dictFile)

	vcfFile := filepath.Join(workDir, "variants.vcf")
	runTool(t, workDir, "gatk", "HaplotypeCaller", "-R", refGenome, "-I", bamFile, "-O", vcfFile)

	stats, err := parseVCFStatistics(vcfFile)
	require.NoError(t, err)

	return &VariantCallingResults{
		VCFFile:        vcfFile,
		VariantCount:   stats.VariantCount,
		SNPCount:       stats.SNPCount,
		IndelCount:     stats.IndelCount,
		QualityScore:   stats.QualityScore,
		FilterPassRate: stats.FilterPassRate,
	}
}

func runBWAAlignment(t *testing.T, workDir, fastqR1, fastqR2, refGenome, aligner string) *AlignmentResults {
	var sam []byte
	switch aligner {
	case "bwa-mem":
		runTool(t, workDir, "bwa", "index", refGenome)
		sam = runTool(t, workDir, "bwa", "mem", refGenome, fastqR1, fastqR2)
	case "bwa-mem2":
		runTool(t, workDir, "bwa-mem2", "index", refGenome)
		sam = runTool(t, workDir, "bwa-mem2", "mem", refGenome, fastqR1, fastqR2)
	case "minimap2":
		sam = runTool(t, workDir, "minimap2", "-ax", "map-ont", refGenome, fastqR1)
	default:
		t.Fatalf("unknown aligner %q", aligner)
	}

	samFile := filepath.Join(workDir, "aligned.sam")
	require.NoError(t, os.WriteFile(samFile, sam, 0o600))
	bamFile := filepath.Join(workDir, "aligned.bam")
	runTool(t, workDir, "samtools", "view", "-bS", samFile, "-o", bamFile)

	mappingRate, err := parseMappingRate(string(runTool(t, workDir, "samtools", "flagstat", bamFile)))
	require.NoError(t, err)
	// Mean MAPQ over mapped primary alignments (-F 0x904: not unmapped,
	// secondary or supplementary).
	avgMAPQ, err := parseMeanMAPQ(string(runTool(t, workDir, "samtools", "view", "-F", "0x904", bamFile)))
	require.NoError(t, err)

	return &AlignmentResults{
		BAMFile:     bamFile,
		MappingRate: mappingRate,
		AverageMAPQ: avgMAPQ,
	}
}

func runBCFtoolsOperations(t *testing.T, workDir string, vcfFiles []string) *BCFtoolsResults {
	// bcftools merge needs bgzipped, indexed inputs.
	compressed := make([]string, 0, len(vcfFiles))
	for _, vcf := range vcfFiles {
		gz := vcf + ".gz"
		runTool(t, workDir, "bcftools", "view", "-Oz", "-o", gz, vcf)
		runTool(t, workDir, "bcftools", "index", gz)
		compressed = append(compressed, gz)
	}

	mergedVCF := filepath.Join(workDir, "merged.vcf")
	runTool(t, workDir, "bcftools", append([]string{"merge", "-o", mergedVCF}, compressed...)...)

	variantStats, err := parseBCFtoolsStats(string(runTool(t, workDir, "bcftools", "stats", mergedVCF)))
	require.NoError(t, err)
	vcfStats, err := parseVCFStatistics(mergedVCF)
	require.NoError(t, err)

	return &BCFtoolsResults{
		MergedVCF:      mergedVCF,
		VariantStats:   variantStats,
		MergedVCFStats: vcfStats,
	}
}

func testSamtoolsBAMProcessing(t *testing.T, framework *pt.E2EFramework) {
	t.Log("🔧 Testing Samtools BAM processing")

	syntheticData := generateOrSkip(t, framework, pt.IlluminaNovaSeq6000, "CEU")

	workDir, err := os.MkdirTemp("", "samtools_test")
	require.NoError(t, err)
	defer os.RemoveAll(workDir)

	bamFile := filepath.Join(workDir, "test.bam")
	writeBAMOrSkip(t, framework, syntheticData, bamFile)

	testSamtoolsOperations(t, workDir, bamFile, len(syntheticData.Sequences))
}

// writeBAMOrSkip writes syntheticData as BAM, skipping the test when samtools
// (needed for the conversion) is not installed and failing on real errors.
func writeBAMOrSkip(t *testing.T, framework *pt.E2EFramework, syntheticData *pt.SyntheticData, bamFile string) {
	t.Helper()
	err := framework.WriteBAM(syntheticData, bamFile)
	if errors.Is(err, pt.ErrBAMUnsupported) {
		t.Skip("samtools not available - cannot write BAM")
	}
	require.NoError(t, err)
}

// Validation functions

func validateGATKResults(t *testing.T, results *GATKResults, platform pt.SequencingPlatform, population string) {
	assert.NotEmpty(t, results.AlignedBAM, "GATK should produce aligned BAM file")
	assert.Greater(t, results.MappingRate, 0.8, "Mapping rate should be >80%")
	assert.Less(t, results.ProcessingTime, 30*time.Minute, "Processing should complete within 30 minutes")

	// Platform-specific validations
	switch platform {
	case pt.IlluminaNovaSeq6000:
		assert.Greater(t, results.MappingRate, 0.95, "NovaSeq should have >95% mapping rate")
	case pt.PacBioRevio:
		assert.Greater(t, results.MappingRate, 0.90, "PacBio should have >90% mapping rate")
	case pt.ONTPromethION24:
		assert.Greater(t, results.MappingRate, 0.85, "ONT should have >85% mapping rate")
	}

	t.Logf("GATK Results - Mapping: %.1f%%, Time: %v", results.MappingRate*100, results.ProcessingTime)
}

// validateVariantCallingResults checks that HaplotypeCaller wrote a parseable
// VCF. The reference is built from the reads themselves, so the reads carry
// no variants relative to it and no variant count is expected.
func validateVariantCallingResults(t *testing.T, results *VariantCallingResults) {
	assert.NotEmpty(t, results.VCFFile, "Should produce VCF file")
	t.Logf("Variant Calling - Variants: %d, SNPs: %d, Indels: %d, Quality: %.1f",
		results.VariantCount, results.SNPCount, results.IndelCount, results.QualityScore)
}

func validateAlignmentResults(t *testing.T, results *AlignmentResults, syntheticData *pt.SyntheticData, platform pt.SequencingPlatform) {
	assert.NotEmpty(t, results.BAMFile, "Should produce BAM file")
	assert.Greater(t, results.MappingRate, 0.8, "Mapping rate should be >80%")
	assert.Greater(t, results.AverageMAPQ, 20.0, "Average MAPQ should be >20")

	// Platform-specific expectations
	switch platform {
	case pt.IlluminaNovaSeq6000, pt.IlluminaNextSeq2000:
		assert.GreaterOrEqual(t, results.MappingRate, 0.95, "Illumina should have high mapping rate")
	case pt.PacBioRevio:
		assert.Greater(t, results.MappingRate, 0.90, "PacBio should have good mapping rate")
		assert.Greater(t, results.AverageMAPQ, 30.0, "PacBio should have high MAPQ")
	case pt.ONTPromethION24:
		assert.Greater(t, results.MappingRate, 0.85, "ONT should have reasonable mapping rate")
	}
}

// validateBCFtoolsResults checks the merge against the inputs. WriteVCF
// writes one transition SNV per 1000 bp of the reference, and both inputs
// share those sites, so the merged file holds the largest input's sites,
// all of them transitions.
func validateBCFtoolsResults(t *testing.T, results *BCFtoolsResults, written []int) {
	want := 0
	for _, n := range written {
		want = max(want, n)
	}
	require.Greater(t, want, 0, "inputs should contain variants")
	assert.NotEmpty(t, results.MergedVCF, "Should produce merged VCF")
	assert.Equal(t, want, results.VariantStats.TotalVariants, "bcftools stats record count")
	assert.Equal(t, want, results.VariantStats.SNPCount, "every input record is a SNP")
	assert.Equal(t, 0, results.VariantStats.IndelCount)
	assert.Equal(t, want, results.VariantStats.TransitionCount, "every input record is a transition")
	assert.Equal(t, 0, results.VariantStats.TransversionCount)
	assert.Equal(t, want, results.MergedVCFStats.VariantCount, "merged VCF record count")
	assert.Equal(t, 1.0, results.MergedVCFStats.FilterPassRate, "all input records PASS")

	t.Logf("BCFtools Results - Variants: %d, Ts: %d, Tv: %d",
		results.VariantStats.TotalVariants, results.VariantStats.TransitionCount, results.VariantStats.TransversionCount)
}

// Helper functions for parsing tool outputs

// parseMappingRate extracts the mapped fraction from samtools flagstat
// output (the "N + M mapped (98.50% : N/A)" line).
func parseMappingRate(flagstatOutput string) (float64, error) {
	for _, line := range strings.Split(flagstatOutput, "\n") {
		if !strings.Contains(line, " mapped (") || strings.Contains(line, "primary mapped") {
			continue
		}
		start := strings.Index(line, "(")
		end := strings.Index(line, "%")
		if start == -1 || end <= start {
			continue
		}
		percent, err := strconv.ParseFloat(strings.TrimSpace(line[start+1:end]), 64)
		if err != nil {
			return 0, fmt.Errorf("parsing mapped percentage in %q: %w", line, err)
		}
		return percent / 100.0, nil
	}
	return 0, fmt.Errorf("no mapped-percentage line in samtools flagstat output")
}

// parseMeanMAPQ returns the mean MAPQ (column 5) of SAM records.
func parseMeanMAPQ(samRecords string) (float64, error) {
	var sum, n int
	for _, line := range strings.Split(samRecords, "\n") {
		if line == "" || strings.HasPrefix(line, "@") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			return 0, fmt.Errorf("SAM record with %d fields: %q", len(fields), line)
		}
		q, err := strconv.Atoi(fields[4])
		if err != nil {
			return 0, fmt.Errorf("bad MAPQ %q: %w", fields[4], err)
		}
		sum += q
		n++
	}
	if n == 0 {
		return 0, fmt.Errorf("no SAM records")
	}
	return float64(sum) / float64(n), nil
}

// VCFStatistics summarises the records of a VCF file.
type VCFStatistics struct {
	VariantCount   int
	SNPCount       int
	IndelCount     int
	QualityScore   float64 // mean QUAL over records with a numeric QUAL
	FilterPassRate float64
}

// parseVCFStatistics counts the records of vcfFile. It fails if the file
// cannot be read or has no #CHROM header; a VCF with no records gives zero
// counts and rates.
func parseVCFStatistics(vcfFile string) (VCFStatistics, error) {
	f, err := os.Open(vcfFile)
	if err != nil {
		return VCFStatistics{}, err
	}
	defer f.Close()

	var st VCFStatistics
	var passCount, qualCount int
	var qualSum float64
	header := false
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "#CHROM") {
			header = true
			continue
		}
		if strings.HasPrefix(line, "#") || line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 8 {
			return VCFStatistics{}, fmt.Errorf("%s: record with %d fields: %q", vcfFile, len(fields), line)
		}
		st.VariantCount++
		if ref, alt := fields[3], fields[4]; len(ref) == 1 && len(alt) == 1 {
			st.SNPCount++
		} else {
			st.IndelCount++
		}
		if q, err := strconv.ParseFloat(fields[5], 64); err == nil {
			qualSum += q
			qualCount++
		}
		if fields[6] == "PASS" || fields[6] == "." {
			passCount++
		}
	}
	if err := sc.Err(); err != nil {
		return VCFStatistics{}, err
	}
	if !header {
		return VCFStatistics{}, fmt.Errorf("%s: no #CHROM header line", vcfFile)
	}
	if qualCount > 0 {
		st.QualityScore = qualSum / float64(qualCount)
	}
	if st.VariantCount > 0 {
		st.FilterPassRate = float64(passCount) / float64(st.VariantCount)
	}
	return st, nil
}

// parseBCFtoolsStats reads the summary-number (SN) and transition/
// transversion (TSTV) sections of bcftools stats output.
func parseBCFtoolsStats(statsOutput string) (VariantStatistics, error) {
	var st VariantStatistics
	sawRecords, sawTSTV := false, false
	atoi := func(s string) (int, error) { return strconv.Atoi(strings.TrimSpace(s)) }
	for _, line := range strings.Split(statsOutput, "\n") {
		fields := strings.Split(line, "\t")
		switch {
		case len(fields) >= 4 && fields[0] == "SN":
			n, err := atoi(fields[3])
			if err != nil {
				continue
			}
			switch fields[2] {
			case "number of records:":
				st.TotalVariants, sawRecords = n, true
			case "number of SNPs:":
				st.SNPCount = n
			case "number of indels:":
				st.IndelCount = n
			}
		case len(fields) >= 4 && fields[0] == "TSTV" && !sawTSTV:
			ts, err1 := atoi(fields[2])
			tv, err2 := atoi(fields[3])
			if err1 != nil || err2 != nil {
				return VariantStatistics{}, fmt.Errorf("bad TSTV line %q", line)
			}
			st.TransitionCount, st.TransversionCount, sawTSTV = ts, tv, true
		}
	}
	if !sawRecords {
		return VariantStatistics{}, fmt.Errorf("no \"number of records\" line in bcftools stats output")
	}
	if !sawTSTV {
		return VariantStatistics{}, fmt.Errorf("no TSTV line in bcftools stats output")
	}
	return st, nil
}

func testSamtoolsOperations(t *testing.T, workDir, bamFile string, wantReads int) {
	run := func(args ...string) string {
		t.Helper()
		cmd := exec.Command("samtools", args...)
		cmd.Dir = workDir
		out, err := cmd.CombinedOutput()
		require.NoError(t, err, "samtools %s: %s", strings.Join(args, " "), out)
		return string(out)
	}

	count := strings.TrimSpace(run("view", "-c", bamFile))
	assert.Equal(t, strconv.Itoa(wantReads), count, "BAM should hold one record per read")

	sorted := filepath.Join(workDir, "sorted.bam")
	run("sort", "-o", sorted, bamFile)
	run("index", sorted)
	_, err := os.Stat(sorted + ".bai")
	require.NoError(t, err, "samtools index should write a .bai")

	flagstat := run("flagstat", bamFile)
	assert.Contains(t, flagstat, fmt.Sprintf("%d + 0 in total", wantReads))
}

func testPicardQualityMetrics(t *testing.T, framework *pt.E2EFramework) {
	t.Log("📊 Testing Picard quality metrics")

	syntheticData := generateOrSkip(t, framework, pt.IlluminaNovaSeq6000, "CEU")

	workDir, err := os.MkdirTemp("", "picard_test")
	require.NoError(t, err)
	defer os.RemoveAll(workDir)

	bamFile := filepath.Join(workDir, "test.bam")
	writeBAMOrSkip(t, framework, syntheticData, bamFile)

	testPicardOperations(t, workDir, bamFile, len(syntheticData.Sequences))
}

func testPicardOperations(t *testing.T, workDir, bamFile string, wantReads int) {
	run := func(args ...string) {
		t.Helper()
		cmd := exec.Command("picard", args...)
		cmd.Dir = workDir
		out, err := cmd.CombinedOutput()
		require.NoError(t, err, "picard %s: %s", strings.Join(args, " "), tailLines(string(out), 20))
	}

	run("ValidateSamFile", "I="+bamFile, "MODE=SUMMARY")

	// The reads are unaligned, so the summary is a single UNPAIRED row whose
	// TOTAL_READS must equal the reads written.
	metrics := filepath.Join(workDir, "alignment_summary.txt")
	run("CollectAlignmentSummaryMetrics", "I="+bamFile, "O="+metrics)
	content, err := os.ReadFile(metrics)
	require.NoError(t, err)
	totalReads := -1
	var columns []string
	for _, line := range strings.Split(string(content), "\n") {
		fields := strings.Split(line, "\t")
		if len(fields) > 1 && fields[0] == "CATEGORY" {
			columns = fields
			continue
		}
		if columns != nil && len(fields) == len(columns) && (fields[0] == "UNPAIRED" || fields[0] == "PAIR") {
			for i, c := range columns {
				if c == "TOTAL_READS" {
					totalReads, _ = strconv.Atoi(fields[i])
				}
			}
			break
		}
	}
	assert.Equal(t, wantReads, totalReads, "Picard TOTAL_READS should match the reads written")
}

func tailLines(s string, n int) string {
	lines := strings.Split(strings.TrimRight(s, "\n"), "\n")
	if len(lines) > n {
		lines = lines[len(lines)-n:]
	}
	return strings.Join(lines, "\n")
}

// testSamtoolsStatistics tests Samtools statistics generation
func testSamtoolsStatistics(t *testing.T, framework *pt.E2EFramework) {
	t.Log("📊 Testing Samtools statistics generation")

	syntheticData := generateOrSkip(t, framework, pt.IlluminaNovaSeq6000, "CEU")

	workDir, err := os.MkdirTemp("", "samtools_stats_test")
	require.NoError(t, err)
	defer os.RemoveAll(workDir)

	bamFile := filepath.Join(workDir, "test.bam")
	writeBAMOrSkip(t, framework, syntheticData, bamFile)

	cmd := exec.Command("samtools", "stats", bamFile)
	cmd.Dir = workDir
	out, err := cmd.Output()
	require.NoError(t, err, "samtools stats")
	assert.Contains(t, string(out), fmt.Sprintf("SN\traw total sequences:\t%d", len(syntheticData.Sequences)),
		"samtools stats should count every read")
}
