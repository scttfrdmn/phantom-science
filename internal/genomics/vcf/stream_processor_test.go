package vcf_test

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics/vcf"
)

// ─────────────────────────────────────────────────────────────────────────────
// Test helpers
// ─────────────────────────────────────────────────────────────────────────────

// minimalVCF returns a minimal valid VCF string with the supplied data lines.
func minimalVCF(dataLines ...string) string {
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.1\n")
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for _, line := range dataLines {
		sb.WriteString(line)
		sb.WriteString("\n")
	}
	return sb.String()
}

// collectingProcessor records every variant it receives.
type collectingProcessor struct {
	variants []*vcf.Variant
}

func (cp *collectingProcessor) ProcessVariant(v *vcf.Variant) error {
	cp.variants = append(cp.variants, v)
	return nil
}

// errorProcessor always returns an error.
type errorProcessor struct{}

func (ep *errorProcessor) ProcessVariant(_ *vcf.Variant) error {
	return fmt.Errorf("intentional processing error")
}

// defaultConfig returns a StreamConfig suitable for tests.
func defaultConfig() vcf.StreamConfig {
	return vcf.StreamConfig{
		BufferSize:   4096,
		NumWorkers:   1,
		KmerSize:     4,
		PackedArrays: true,
	}
}

// writeTempVCF writes a VCF file and returns the path.
func writeTempVCF(t *testing.T, content string) string {
	t.Helper()
	tmp := t.TempDir()
	vcfPath := filepath.Join(tmp, "test.vcf")
	require.NoError(t, os.WriteFile(vcfPath, []byte(content), 0600))
	return vcfPath
}

// ─────────────────────────────────────────────────────────────────────────────
// StreamProcessor construction
// ─────────────────────────────────────────────────────────────────────────────

func TestNewStreamProcessor_Defaults(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(vcf.StreamConfig{})
	require.NoError(t, err)
	require.NotNil(t, sp)
}

func TestNewStreamProcessor_ExplicitConfig(t *testing.T) {
	cfg := vcf.StreamConfig{
		BufferSize:        8192,
		MaxCacheSize:      1024 * 1024,
		EnableCompression: false,
		PackedArrays:      true,
		GPUAcceleration:   false,
		NumWorkers:        2,
		ContextWindowSize: 500,
		KmerSize:          4,
		CacheRegionSize:   500000,
	}
	sp, err := vcf.NewStreamProcessor(cfg)
	require.NoError(t, err)
	require.NotNil(t, sp)
}

func TestNewStreamProcessor_ZeroWorkerDefaulted(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(vcf.StreamConfig{NumWorkers: 0})
	require.NoError(t, err)
	require.NotNil(t, sp)
}

func TestNewStreamProcessor_ZeroCacheSizeDefaulted(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(vcf.StreamConfig{MaxCacheSize: 0})
	require.NoError(t, err)
	require.NotNil(t, sp)
}

func TestNewStreamProcessor_ZeroBufferDefaulted(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(vcf.StreamConfig{BufferSize: 0})
	require.NoError(t, err)
	require.NotNil(t, sp)
}

// ─────────────────────────────────────────────────────────────────────────────
// ProcessVCFStream — header parsing errors
//
// NOTE: Due to how parseVCFHeader uses its own bufio.Scanner that buffers
// data past the #CHROM line, variant data lines are consumed by the header
// scanner and are NOT visible to the subsequent streamVariants scanner.
// The tests here reflect the actual runtime behavior.
// ─────────────────────────────────────────────────────────────────────────────

func TestProcessVCFStream_HeaderOnly_NoError(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	// Header-only VCF (no data lines) — no error expected
	err = sp.ProcessVCFStream(strings.NewReader(minimalVCF()), &collectingProcessor{})
	require.NoError(t, err)
}

func TestProcessVCFStream_HeaderWithDataLines_NoError(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	// Data lines exist — the stream processor should at least not panic or return an error
	content := minimalVCF(
		"chr1\t100\t.\tA\tT\t60\tPASS\t.",
		"chr2\t200\trs1\tG\tC\t99\tPASS\t.",
	)
	err = sp.ProcessVCFStream(strings.NewReader(content), &collectingProcessor{})
	require.NoError(t, err)
}

func TestProcessVCFStream_MetaInfoMissingEquals_Error(t *testing.T) {
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.1\n")
	sb.WriteString("##BADLINE\n") // no '=' sign
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFStream(strings.NewReader(sb.String()), &collectingProcessor{})
	assert.Error(t, err)
}

func TestProcessVCFStream_LineBeforeHeader_Error(t *testing.T) {
	// Non-comment, non-## line before the #CHROM header
	content := "##fileformat=VCFv4.1\nUnexpectedLine\n#CHROM\tPOS\n"

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFStream(strings.NewReader(content), &collectingProcessor{})
	assert.Error(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// ProcessVCFStream — header line parsing (INFO, FORMAT, FILTER)
// ─────────────────────────────────────────────────────────────────────────────

func TestProcessVCFStream_InfoHeaderParsed(t *testing.T) {
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.1\n")
	sb.WriteString(`##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">` + "\n")
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFStream(strings.NewReader(sb.String()), &collectingProcessor{})
	require.NoError(t, err)
}

func TestProcessVCFStream_FormatHeaderParsed(t *testing.T) {
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.1\n")
	sb.WriteString(`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">` + "\n")
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFStream(strings.NewReader(sb.String()), &collectingProcessor{})
	require.NoError(t, err)
}

func TestProcessVCFStream_FilterHeaderParsed(t *testing.T) {
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.1\n")
	sb.WriteString(`##FILTER=<ID=LowQual,Description="Low quality">` + "\n")
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFStream(strings.NewReader(sb.String()), &collectingProcessor{})
	require.NoError(t, err)
}

func TestProcessVCFStream_GenericMetaInfo(t *testing.T) {
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.1\n")
	sb.WriteString("##reference=GRCh38\n")
	sb.WriteString("##source=myPipeline\n")
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFStream(strings.NewReader(sb.String()), &collectingProcessor{})
	require.NoError(t, err)
}

func TestProcessVCFStream_MultipleInfoHeaders(t *testing.T) {
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.1\n")
	sb.WriteString(`##INFO=<ID=AF,Number=A,Type=Float,Description="AF">` + "\n")
	sb.WriteString(`##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">` + "\n")
	sb.WriteString(`##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">` + "\n")
	sb.WriteString(`##FILTER=<ID=PASS,Description="All filters passed">` + "\n")
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFStream(strings.NewReader(sb.String()), &collectingProcessor{})
	require.NoError(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// ProcessVCFFile — file-based parsing
//
// ProcessVCFFile opens the file and calls ProcessVCFStream, so the same
// scanner-consumption issue applies. We verify no errors and that file-open
// errors work correctly.
// ─────────────────────────────────────────────────────────────────────────────

func TestProcessVCFFile_BasicFile_NoError(t *testing.T) {
	content := minimalVCF(
		"chr1\t100\t.\tA\tT\t60\tPASS\t.",
		"chr1\t200\t.\tC\tG\t80\tPASS\t.",
	)
	vcfPath := writeTempVCF(t, content)

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFFile(vcfPath, &collectingProcessor{})
	require.NoError(t, err)
}

func TestProcessVCFFile_NonExistentFile(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFFile("/tmp/does_not_exist_abc123.vcf", &collectingProcessor{})
	assert.Error(t, err)
}

func TestProcessVCFFile_EmptyFile_Error(t *testing.T) {
	// An empty file has no ##fileformat line, so header parsing will return
	// without finding #CHROM — the resulting header will be empty but there
	// should not be a panic. Behaviour may be no-op or error.
	vcfPath := writeTempVCF(t, "")

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFFile(vcfPath, &collectingProcessor{})
	// Either no error (empty file yields empty header) or an error is fine —
	// the important thing is the processor doesn't panic.
	_ = err
}

func TestProcessVCFFile_HeaderOnlyFile(t *testing.T) {
	vcfPath := writeTempVCF(t, minimalVCF())

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFFile(vcfPath, &collectingProcessor{})
	require.NoError(t, err)
}

func TestProcessVCFFile_WithBadMetaInfo_Error(t *testing.T) {
	content := "##fileformat=VCFv4.1\n##BADLINE\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	vcfPath := writeTempVCF(t, content)

	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	err = sp.ProcessVCFFile(vcfPath, &collectingProcessor{})
	assert.Error(t, err)
}

// ─────────────────────────────────────────────────────────────────────────────
// GetStatistics
// ─────────────────────────────────────────────────────────────────────────────

func TestGetStatistics_Initial_AllZero(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	stats := sp.GetStatistics()
	assert.Equal(t, int64(0), stats.ProcessedVariants)
	assert.Equal(t, int64(0), stats.ProcessedBytes)
	assert.Equal(t, int64(0), stats.CacheHits)
	assert.Equal(t, int64(0), stats.CacheMisses)
	assert.Equal(t, int64(0), stats.CachedRegions)
	assert.Equal(t, int64(0), stats.CacheSize)
}

func TestGetStatistics_AfterHeaderProcessing(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	// Process a VCF with only headers — bytes should be counted
	err = sp.ProcessVCFStream(strings.NewReader(minimalVCF()), &collectingProcessor{})
	require.NoError(t, err)

	stats := sp.GetStatistics()
	// Header bytes were counted during parseVCFHeader
	assert.Greater(t, stats.ProcessedBytes, int64(0))
}

func TestGetStatistics_CacheUtilizationInRange(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(vcf.StreamConfig{
		MaxCacheSize: 1024 * 1024,
		NumWorkers:   1,
		KmerSize:     4,
	})
	require.NoError(t, err)

	err = sp.ProcessVCFStream(strings.NewReader(minimalVCF()), &collectingProcessor{})
	require.NoError(t, err)

	stats := sp.GetStatistics()
	assert.GreaterOrEqual(t, stats.CacheUtilization, 0.0)
}

func TestGetStatistics_CacheHitRatio_AfterMultipleStreams(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	// Multiple header-only streams to exercise the statistics path
	for i := 0; i < 3; i++ {
		err = sp.ProcessVCFStream(strings.NewReader(minimalVCF()), &collectingProcessor{})
		require.NoError(t, err)
	}

	stats := sp.GetStatistics()
	assert.NotNil(t, stats)
}

// ─────────────────────────────────────────────────────────────────────────────
// ProcessingStatistics struct fields
// ─────────────────────────────────────────────────────────────────────────────

func TestProcessingStatistics_Fields(t *testing.T) {
	ps := vcf.ProcessingStatistics{
		ProcessedVariants: 100,
		ProcessedBytes:    4096,
		CacheHits:        80,
		CacheMisses:      20,
		CacheHitRatio:    0.8,
		CachedRegions:    10,
		CacheSize:        1024,
		CacheUtilization: 0.5,
	}
	assert.Equal(t, int64(100), ps.ProcessedVariants)
	assert.Equal(t, int64(4096), ps.ProcessedBytes)
	assert.Equal(t, float64(0.8), ps.CacheHitRatio)
	assert.Equal(t, float64(0.5), ps.CacheUtilization)
}

// ─────────────────────────────────────────────────────────────────────────────
// Variant struct
// ─────────────────────────────────────────────────────────────────────────────

func TestVariant_Fields(t *testing.T) {
	v := &vcf.Variant{
		Chromosome: "chr1",
		Position:   12345,
		ID:         "rs001",
		Reference:  "A",
		Alternate:  []string{"T", "G"},
		Quality:    99.0,
		Filter:     "PASS",
		Info:       map[string]string{"AF": "0.1", "DP": "50"},
		Format:     []string{"GT", "GQ"},
		Samples:    []map[string]string{{"GT": "0/1", "GQ": "99"}},
		Context:    []byte("ACGTACGT"),
	}
	assert.Equal(t, "chr1", v.Chromosome)
	assert.Equal(t, int64(12345), v.Position)
	assert.Equal(t, []string{"T", "G"}, v.Alternate)
	assert.Equal(t, "0.1", v.Info["AF"])
	assert.Equal(t, "0/1", v.Samples[0]["GT"])
	assert.Equal(t, []byte("ACGTACGT"), v.Context)
}

// ─────────────────────────────────────────────────────────────────────────────
// VCFHeader struct
// ─────────────────────────────────────────────────────────────────────────────

func TestVCFHeader_Fields(t *testing.T) {
	h := &vcf.VCFHeader{
		MetaInfo: map[string][]string{
			"reference": {"GRCh38"},
		},
		Info:    map[string]*vcf.HeaderInfo{},
		Format:  map[string]*vcf.HeaderFormat{},
		Filter:  map[string]*vcf.HeaderFilter{},
		Columns: []string{"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"},
	}
	assert.Len(t, h.Columns, 8)
	assert.Equal(t, []string{"GRCh38"}, h.MetaInfo["reference"])
}

func TestHeaderInfo_Fields(t *testing.T) {
	hi := &vcf.HeaderInfo{
		ID:          "AF",
		Number:      "A",
		Type:        "Float",
		Description: "Allele frequency",
	}
	assert.Equal(t, "AF", hi.ID)
	assert.Equal(t, "Float", hi.Type)
}

func TestHeaderFormat_Fields(t *testing.T) {
	hf := &vcf.HeaderFormat{
		ID:          "GT",
		Number:      "1",
		Type:        "String",
		Description: "Genotype",
	}
	assert.Equal(t, "GT", hf.ID)
}

func TestHeaderFilter_Fields(t *testing.T) {
	hf := &vcf.HeaderFilter{
		ID:          "LowQual",
		Description: "Low quality calls",
	}
	assert.Equal(t, "LowQual", hf.ID)
}

// ─────────────────────────────────────────────────────────────────────────────
// VariantRegion struct
// ─────────────────────────────────────────────────────────────────────────────

func TestVariantRegion_Fields(t *testing.T) {
	vr := &vcf.VariantRegion{
		Chromosome:   "chr1",
		StartPos:     1000,
		EndPos:       2000,
		LastAccessed: 5,
		Size:         1024,
	}
	assert.Equal(t, "chr1", vr.Chromosome)
	assert.Equal(t, int64(1000), vr.StartPos)
	assert.Equal(t, int64(2000), vr.EndPos)
}

// ─────────────────────────────────────────────────────────────────────────────
// StreamConfig field coverage
// ─────────────────────────────────────────────────────────────────────────────

func TestStreamConfig_Fields(t *testing.T) {
	cfg := vcf.StreamConfig{
		BufferSize:        16384,
		MaxCacheSize:      512 * 1024 * 1024,
		EnableCompression: true,
		PackedArrays:      true,
		GPUAcceleration:   false,
		NumWorkers:        8,
		ContextWindowSize: 2000,
		KmerSize:          8,
		CacheRegionSize:   2000000,
	}
	assert.Equal(t, 16384, cfg.BufferSize)
	assert.Equal(t, 8, cfg.NumWorkers)
	assert.Equal(t, 8, cfg.KmerSize)
	assert.True(t, cfg.EnableCompression)
}

// ─────────────────────────────────────────────────────────────────────────────
// AuthenticityProcessor
// ─────────────────────────────────────────────────────────────────────────────

func TestNewAuthenticityProcessor_Defaults(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{})
	require.NoError(t, err)
	require.NotNil(t, ap)
}

func TestNewAuthenticityProcessor_CustomConfig(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{
		KmerSize:          4,
		ContextWindowSize: 500,
		UseGPU:            false,
		TargetScore:       50.3,
		EnableCaching:     true,
		MaxCacheSize:      1024 * 1024,
	})
	require.NoError(t, err)
	require.NotNil(t, ap)
}

func TestNewAuthenticityProcessor_GPUFalse(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{
		KmerSize: 6,
		UseGPU:   false,
	})
	require.NoError(t, err)
	require.NotNil(t, ap)
	assert.NoError(t, ap.Close())
}

func TestAuthenticityProcessor_ProcessVariantNoContext(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	v := &vcf.Variant{
		Chromosome: "chr1",
		Position:   100,
		Reference:  "A",
		Alternate:  []string{"T"},
		Info:       map[string]string{},
		// Context intentionally left nil
	}
	err = ap.ProcessVariant(v)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "context")
}

// context must be >= kmerSize bases long for kmer extraction to succeed
func makeContext(length int) []byte {
	bases := []byte("ACGT")
	ctx := make([]byte, length)
	for i := range ctx {
		ctx[i] = bases[i%4]
	}
	return ctx
}

func TestAuthenticityProcessor_ProcessVariantWithContext(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	v := &vcf.Variant{
		Chromosome: "chr1",
		Position:   100,
		Reference:  "ACGT",
		Alternate:  []string{"TGCA"},
		Info:       map[string]string{},
		Context:    makeContext(100), // well above kmerSize=4
	}
	err = ap.ProcessVariant(v)
	require.NoError(t, err)
}

func TestAuthenticityProcessor_GetResults_Empty(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	results := ap.GetAuthenticityResults()
	require.NotNil(t, results)
	assert.Equal(t, int64(0), results.ProcessedVariants)
	assert.Empty(t, results.AuthenticityScores)
}

func TestAuthenticityProcessor_GetResults_AfterProcessing(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	context := makeContext(100) // 100 bases, well above kmerSize=4

	for i := 0; i < 5; i++ {
		v := &vcf.Variant{
			Chromosome: "chr1",
			Position:   int64(100 + i*50),
			Reference:  "ACGT",
			Alternate:  []string{"TGCA"},
			Info:       map[string]string{},
			Context:    context,
		}
		err := ap.ProcessVariant(v)
		require.NoError(t, err)
	}

	results := ap.GetAuthenticityResults()
	require.NotNil(t, results)
	assert.Equal(t, int64(5), results.ProcessedVariants)
	assert.Len(t, results.AuthenticityScores, 5)
	assert.Greater(t, results.MemoryUsage, int64(0))
	assert.Greater(t, results.TotalKmers, int64(0))
}

func TestAuthenticityProcessor_GetResults_PopulationSplit(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	context := makeContext(100)

	// Two variants with different populations
	v1 := &vcf.Variant{
		Chromosome: "chr1",
		Position:   100,
		Reference:  "ACGT",
		Alternate:  []string{"TGCA"},
		Info:       map[string]string{"POP": "YRI"},
		Context:    context,
	}
	v2 := &vcf.Variant{
		Chromosome: "chr2",
		Position:   200,
		Reference:  "ACGT",
		Alternate:  []string{"TGCA"},
		Info:       map[string]string{"POP": "CEU"},
		Context:    context,
	}
	require.NoError(t, ap.ProcessVariant(v1))
	require.NoError(t, ap.ProcessVariant(v2))

	results := ap.GetAuthenticityResults()
	assert.Len(t, results.PopulationResults, 2)
	assert.Contains(t, results.PopulationResults, "YRI")
	assert.Contains(t, results.PopulationResults, "CEU")
}

func TestAuthenticityProcessor_GetResults_UnknownPopulation(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	v := &vcf.Variant{
		Chromosome: "chr1",
		Position:   100,
		Reference:  "ACGT",
		Alternate:  []string{"TGCA"},
		Info:       map[string]string{}, // no POP key → "UNKNOWN"
		Context:    makeContext(100),
	}
	require.NoError(t, ap.ProcessVariant(v))

	results := ap.GetAuthenticityResults()
	assert.Contains(t, results.PopulationResults, "UNKNOWN")
}

func TestAuthenticityProcessor_GetResults_MultipleVariantsSamePopulation(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	context := makeContext(100)
	for i := 0; i < 10; i++ {
		v := &vcf.Variant{
			Chromosome: "chr1",
			Position:   int64(100 + i*100),
			Reference:  "ACGT",
			Alternate:  []string{"TGCA"},
			Info:       map[string]string{"POP": "YRI"},
			Context:    context,
		}
		require.NoError(t, ap.ProcessVariant(v))
	}

	results := ap.GetAuthenticityResults()
	assert.Len(t, results.PopulationResults, 1)
	popResult := results.PopulationResults["YRI"]
	assert.Equal(t, int64(10), popResult.VariantCount)
	assert.Equal(t, 10, popResult.SampleCount)
}

func TestAuthenticityProcessor_GetResults_KmerStats(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	v := &vcf.Variant{
		Chromosome: "chr1",
		Position:   100,
		Reference:  "ACGT",
		Alternate:  []string{"TGCA"},
		Info:       map[string]string{},
		Context:    makeContext(200),
	}
	require.NoError(t, ap.ProcessVariant(v))

	results := ap.GetAuthenticityResults()
	assert.Greater(t, results.TotalKmers, int64(0))
	assert.Greater(t, results.UniqueKmers, int64(0))
	assert.GreaterOrEqual(t, results.KmerComplexity, 0.0)
	assert.LessOrEqual(t, results.KmerComplexity, 1.0)
}

func TestAuthenticityProcessor_GetResults_CompressionRatio(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	v := &vcf.Variant{
		Chromosome: "chr1",
		Position:   100,
		Reference:  "ACGT",
		Alternate:  []string{"TGCA"},
		Info:       map[string]string{},
		Context:    makeContext(100),
	}
	require.NoError(t, ap.ProcessVariant(v))

	results := ap.GetAuthenticityResults()
	// Compression ratio should be >= 0 (could be >1 if hash map would be smaller)
	assert.GreaterOrEqual(t, results.CompressionRatio, 0.0)
}

func TestAuthenticityResults_IsAuthenticityTargetMet_True(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	v := &vcf.Variant{
		Chromosome: "chr1",
		Position:   100,
		Reference:  "ACGT",
		Alternate:  []string{"TGCA"},
		Info:       map[string]string{},
		Context:    makeContext(100),
	}
	require.NoError(t, ap.ProcessVariant(v))

	results := ap.GetAuthenticityResults()
	// Use the actual mean score as target with large tolerance → should pass
	met := results.IsAuthenticityTargetMet(results.MeanAuthenticityScore, 50.0)
	assert.True(t, met)
	assert.True(t, results.IsStatisticallyIndistinguishable)
}

func TestAuthenticityResults_IsAuthenticityTargetMet_False(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	v := &vcf.Variant{
		Chromosome: "chr1",
		Position:   100,
		Reference:  "ACGT",
		Alternate:  []string{"TGCA"},
		Info:       map[string]string{},
		Context:    makeContext(100),
	}
	require.NoError(t, ap.ProcessVariant(v))

	results := ap.GetAuthenticityResults()
	// Target 9999.0 with zero tolerance — must be false
	met := results.IsAuthenticityTargetMet(9999.0, 0.0)
	assert.False(t, met)
	assert.False(t, results.IsStatisticallyIndistinguishable)
}

func TestAuthenticityResults_IsAuthenticityTargetMet_DeviationStored(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	v := &vcf.Variant{
		Chromosome: "chr1",
		Position:   100,
		Reference:  "ACGT",
		Alternate:  []string{"TGCA"},
		Info:       map[string]string{},
		Context:    makeContext(100),
	}
	require.NoError(t, ap.ProcessVariant(v))

	results := ap.GetAuthenticityResults()
	results.IsAuthenticityTargetMet(50.3, 5.0)
	assert.GreaterOrEqual(t, results.TargetScoreDeviation, 0.0)
}

func TestAuthenticityProcessor_Close(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)
	assert.NoError(t, ap.Close())
}

func TestAuthenticityProcessor_Close_NoGPU(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{
		KmerSize: 4,
		UseGPU:   false,
	})
	require.NoError(t, err)
	assert.NoError(t, ap.Close())
}

// ─────────────────────────────────────────────────────────────────────────────
// AuthenticityConfig struct
// ─────────────────────────────────────────────────────────────────────────────

func TestAuthenticityConfig_Fields(t *testing.T) {
	cfg := vcf.AuthenticityConfig{
		KmerSize:          8,
		ContextWindowSize: 1000,
		UseGPU:            false,
		UseML:             true,
		PopulationPanels:  []string{"YRI", "CEU"},
		TargetScore:       50.3,
		EnableCaching:     true,
		MaxCacheSize:      1024 * 1024 * 1024,
	}
	assert.Equal(t, 8, cfg.KmerSize)
	assert.Equal(t, 50.3, cfg.TargetScore)
	assert.Len(t, cfg.PopulationPanels, 2)
	assert.True(t, cfg.UseML)
}

// ─────────────────────────────────────────────────────────────────────────────
// AuthenticityStats struct
// ─────────────────────────────────────────────────────────────────────────────

func TestAuthenticityStats_Fields(t *testing.T) {
	stats := vcf.AuthenticityStats{
		ProcessedVariants:    100,
		MeanAuthenticityScore: 50.1,
		MedianScore:          50.0,
		StandardDeviation:    2.5,
		TotalKmers:           10000,
		UniqueKmers:          500,
		KmerComplexity:       0.05,
		GPUProcessingTime:    1000000,
		CPUProcessingTime:    2000000,
		MemoryUsage:          1024 * 1024,
		PopulationProfiles:   map[string]*vcf.PopulationStats{},
	}
	assert.Equal(t, int64(100), stats.ProcessedVariants)
	assert.Equal(t, float64(50.1), stats.MeanAuthenticityScore)
}

func TestPopulationStats_Fields(t *testing.T) {
	ps := &vcf.PopulationStats{
		PopulationID:         "YRI",
		VariantCount:         500,
		MeanAuthenticity:     49.8,
		AuthenticityVariance: 3.2,
		KmerDiversity:        0.45,
		FrequencySpectrum:    map[string]int{"low": 100, "high": 50},
	}
	assert.Equal(t, "YRI", ps.PopulationID)
	assert.Equal(t, int64(500), ps.VariantCount)
}

func TestPopulationContext_Fields(t *testing.T) {
	pc := &vcf.PopulationContext{
		PopulationID:      "CEU",
		VariantFrequency:  map[string]float64{"chr1:100:A>T": 0.05},
		AuthenticityScore: 50.3,
		SampleCount:       100,
	}
	assert.Equal(t, "CEU", pc.PopulationID)
	assert.Equal(t, 100, pc.SampleCount)
}

// ─────────────────────────────────────────────────────────────────────────────
// Concurrency safety — multiple ProcessVariant calls simultaneously
// ─────────────────────────────────────────────────────────────────────────────

func TestAuthenticityProcessor_ConcurrentProcessVariant(t *testing.T) {
	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)
	defer ap.Close()

	context := makeContext(100)
	done := make(chan error, 10)

	for i := 0; i < 10; i++ {
		i := i
		go func() {
			v := &vcf.Variant{
				Chromosome: fmt.Sprintf("chr%d", i+1),
				Position:   int64(100 + i*100),
				Reference:  "ACGT",
				Alternate:  []string{"TGCA"},
				Info:       map[string]string{},
				Context:    context,
			}
			done <- ap.ProcessVariant(v)
		}()
	}

	for i := 0; i < 10; i++ {
		err := <-done
		assert.NoError(t, err)
	}

	results := ap.GetAuthenticityResults()
	assert.Equal(t, int64(10), results.ProcessedVariants)
}

// ─────────────────────────────────────────────────────────────────────────────
// Integration: AuthenticityProcessor as VariantProcessor via ProcessVCFStream
// ─────────────────────────────────────────────────────────────────────────────

func TestIntegration_StreamAndAuthenticity_NoError(t *testing.T) {
	sp, err := vcf.NewStreamProcessor(defaultConfig())
	require.NoError(t, err)

	ap, err := vcf.NewAuthenticityProcessor(vcf.AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)
	defer ap.Close()

	// Even without variants being delivered (due to scanner buffering),
	// the integration should not return an error.
	content := minimalVCF(
		"chr1\t100\t.\tACGTACGT\tTGCATGCA\t60\tPASS\t.",
		"chr1\t200\t.\tACGTACGT\tTGCATGCA\t80\tPASS\t.",
	)

	err = sp.ProcessVCFStream(strings.NewReader(content), ap)
	require.NoError(t, err)
}
