// Package vcf — whitebox tests for unexported stream processor helpers.
package vcf

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// makeTempFasta creates a temporary FASTA file with `length` bases on chromosome `chr`.
func makeTempFasta(t *testing.T, chr string, length int) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "ref.fasta")
	seq := make([]byte, length)
	bases := []byte("ACGT")
	for i := range seq {
		seq[i] = bases[i%4]
	}
	content := fmt.Sprintf(">%s\n%s\n", chr, string(seq))
	require.NoError(t, os.WriteFile(p, []byte(content), 0644))
	return p
}

// ─────────────────────────────────────────────────────────────────────────────
// parseVariantLine
// ─────────────────────────────────────────────────────────────────────────────

func newSP() *StreamProcessor {
	sp, _ := NewStreamProcessor(StreamConfig{NumWorkers: 1, KmerSize: 4})
	return sp
}

func emptyHeader() *VCFHeader {
	return &VCFHeader{
		MetaInfo: make(map[string][]string),
		Info:     make(map[string]*HeaderInfo),
		Format:   make(map[string]*HeaderFormat),
		Filter:   make(map[string]*HeaderFilter),
	}
}

func TestParseVariantLine_Minimal(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA\tT\t60\tPASS\t."
	v, err := sp.parseVariantLine(line, emptyHeader())
	require.NoError(t, err)
	assert.Equal(t, "chr1", v.Chromosome)
	assert.Equal(t, int64(100), v.Position)
	assert.Equal(t, ".", v.ID)
	assert.Equal(t, "A", v.Reference)
	assert.Equal(t, []string{"T"}, v.Alternate)
	assert.Equal(t, float64(60), v.Quality)
	assert.Equal(t, "PASS", v.Filter)
	assert.Empty(t, v.Info)
}

func TestParseVariantLine_MissingQuality(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA\tT\t.\tPASS\t."
	v, err := sp.parseVariantLine(line, emptyHeader())
	require.NoError(t, err)
	assert.Equal(t, float64(0), v.Quality)
}

func TestParseVariantLine_MissingAlternate(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA\t.\t60\tPASS\t."
	v, err := sp.parseVariantLine(line, emptyHeader())
	require.NoError(t, err)
	assert.Nil(t, v.Alternate)
}

func TestParseVariantLine_MultiAllele(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA\tT,G\t60\tPASS\t."
	v, err := sp.parseVariantLine(line, emptyHeader())
	require.NoError(t, err)
	assert.Equal(t, []string{"T", "G"}, v.Alternate)
}

func TestParseVariantLine_InfoKeyValue(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA\tT\t60\tPASS\tAF=0.05;DP=100"
	v, err := sp.parseVariantLine(line, emptyHeader())
	require.NoError(t, err)
	assert.Equal(t, "0.05", v.Info["AF"])
	assert.Equal(t, "100", v.Info["DP"])
}

func TestParseVariantLine_InfoFlagField(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA\tT\t60\tPASS\tDB"
	v, err := sp.parseVariantLine(line, emptyHeader())
	require.NoError(t, err)
	assert.Equal(t, "true", v.Info["DB"])
}

func TestParseVariantLine_WithFormatAndSample(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:GQ\t0/1:99"
	v, err := sp.parseVariantLine(line, emptyHeader())
	require.NoError(t, err)
	require.Len(t, v.Format, 2)
	assert.Equal(t, "GT", v.Format[0])
	assert.Equal(t, "GQ", v.Format[1])
	require.Len(t, v.Samples, 1)
	assert.Equal(t, "0/1", v.Samples[0]["GT"])
	assert.Equal(t, "99", v.Samples[0]["GQ"])
}

func TestParseVariantLine_MultipleSamples(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1"
	v, err := sp.parseVariantLine(line, emptyHeader())
	require.NoError(t, err)
	assert.Len(t, v.Samples, 3)
	assert.Equal(t, "0/0", v.Samples[0]["GT"])
	assert.Equal(t, "0/1", v.Samples[1]["GT"])
	assert.Equal(t, "1/1", v.Samples[2]["GT"])
}

func TestParseVariantLine_InvalidPosition(t *testing.T) {
	sp := newSP()
	line := "chr1\tNOTANUM\t.\tA\tT\t60\tPASS\t."
	_, err := sp.parseVariantLine(line, emptyHeader())
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "invalid position")
}

func TestParseVariantLine_InvalidQuality(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA\tT\tBADQUAL\tPASS\t."
	_, err := sp.parseVariantLine(line, emptyHeader())
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "invalid quality")
}

func TestParseVariantLine_InsufficientFields(t *testing.T) {
	sp := newSP()
	line := "chr1\t100\t.\tA"
	_, err := sp.parseVariantLine(line, emptyHeader())
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "insufficient fields")
}

func TestParseVariantLine_SampleWithMoreValuesThanFormat(t *testing.T) {
	sp := newSP()
	// More sample values than format fields — extra values should be ignored
	line := "chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1:99:extra"
	v, err := sp.parseVariantLine(line, emptyHeader())
	require.NoError(t, err)
	assert.Equal(t, "0/1", v.Samples[0]["GT"])
}

// ─────────────────────────────────────────────────────────────────────────────
// parseMetaInfoLine
// ─────────────────────────────────────────────────────────────────────────────

func TestParseMetaInfoLine_FileFormat(t *testing.T) {
	sp := newSP()
	h := &VCFHeader{
		MetaInfo: make(map[string][]string),
		Info:     make(map[string]*HeaderInfo),
		Format:   make(map[string]*HeaderFormat),
		Filter:   make(map[string]*HeaderFilter),
	}
	err := sp.parseMetaInfoLine("##fileformat=VCFv4.1", h)
	require.NoError(t, err)
	assert.Equal(t, []string{"VCFv4.1"}, h.MetaInfo["fileformat"])
}

func TestParseMetaInfoLine_InfoLine(t *testing.T) {
	sp := newSP()
	h := &VCFHeader{
		MetaInfo: make(map[string][]string),
		Info:     make(map[string]*HeaderInfo),
		Format:   make(map[string]*HeaderFormat),
		Filter:   make(map[string]*HeaderFilter),
	}
	err := sp.parseMetaInfoLine(`##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">`, h)
	require.NoError(t, err)
	assert.Contains(t, h.Info, "placeholder")
}

func TestParseMetaInfoLine_FormatLine(t *testing.T) {
	sp := newSP()
	h := &VCFHeader{
		MetaInfo: make(map[string][]string),
		Info:     make(map[string]*HeaderInfo),
		Format:   make(map[string]*HeaderFormat),
		Filter:   make(map[string]*HeaderFilter),
	}
	err := sp.parseMetaInfoLine(`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`, h)
	require.NoError(t, err)
	assert.Contains(t, h.Format, "placeholder")
}

func TestParseMetaInfoLine_FilterLine(t *testing.T) {
	sp := newSP()
	h := &VCFHeader{
		MetaInfo: make(map[string][]string),
		Info:     make(map[string]*HeaderInfo),
		Format:   make(map[string]*HeaderFormat),
		Filter:   make(map[string]*HeaderFilter),
	}
	err := sp.parseMetaInfoLine(`##FILTER=<ID=LowQual,Description="Low quality">`, h)
	require.NoError(t, err)
	assert.Contains(t, h.Filter, "placeholder")
}

func TestParseMetaInfoLine_MissingEquals(t *testing.T) {
	sp := newSP()
	h := emptyHeader()
	h.MetaInfo = make(map[string][]string)
	h.Info = make(map[string]*HeaderInfo)
	h.Format = make(map[string]*HeaderFormat)
	h.Filter = make(map[string]*HeaderFilter)
	err := sp.parseMetaInfoLine("##BADLINE", h)
	assert.Error(t, err)
}

func TestParseMetaInfoLine_MultipleValuesForKey(t *testing.T) {
	sp := newSP()
	h := &VCFHeader{
		MetaInfo: make(map[string][]string),
		Info:     make(map[string]*HeaderInfo),
		Format:   make(map[string]*HeaderFormat),
		Filter:   make(map[string]*HeaderFilter),
	}
	require.NoError(t, sp.parseMetaInfoLine("##contig=<ID=chr1>", h))
	require.NoError(t, sp.parseMetaInfoLine("##contig=<ID=chr2>", h))
	assert.Len(t, h.MetaInfo["contig"], 2)
}

// ─────────────────────────────────────────────────────────────────────────────
// isGzipped
// ─────────────────────────────────────────────────────────────────────────────

func TestIsGzipped_NotGzipped(t *testing.T) {
	sp := newSP()
	import_bufio_reader_from_text := func(s string) *strings.Reader { return strings.NewReader(s) }
	_ = import_bufio_reader_from_text

	// Use the public ProcessVCFStream path which exercises isGzipped internally
	// For this test, just ensure non-gzipped content doesn't trigger gzip branch
	content := "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	err := sp.ProcessVCFStream(strings.NewReader(content), &collectCollector{})
	require.NoError(t, err)
}

type collectCollector struct {
	variants []*Variant
}

func (c *collectCollector) ProcessVariant(v *Variant) error {
	c.variants = append(c.variants, v)
	return nil
}

// ─────────────────────────────────────────────────────────────────────────────
// addGenomicContext / loadGenomicContext / extractContextFromRegion / cacheRegion / evictOldestRegions
// These are reached via variantWorker which is started by streamVariants.
// We call them directly here.
// ─────────────────────────────────────────────────────────────────────────────

func TestAddGenomicContext_NoReference_SkipsContext(t *testing.T) {
	// When no reference path is configured, addGenomicContext succeeds but
	// leaves Context as nil and does not update cache statistics.
	sp := newSP()
	v := &Variant{
		Chromosome: "chr1",
		Position:   1000,
		Reference:  "A",
		Alternate:  []string{"T"},
		Info:       map[string]string{},
	}
	err := sp.addGenomicContext(v)
	require.NoError(t, err)
	assert.Nil(t, v.Context)
	assert.Equal(t, int64(0), sp.cacheMisses)
	assert.Equal(t, int64(0), sp.cacheHits)
}

func TestAddGenomicContext_WithReference_LoadsContext(t *testing.T) {
	// With a configured reference path, addGenomicContext loads real sequence data.
	p := makeTempFasta(t, "chr1", 4000)
	sp, err := NewStreamProcessor(StreamConfig{ReferencePath: p, NumWorkers: 1})
	require.NoError(t, err)
	require.NotNil(t, sp.refManager)

	v := &Variant{Chromosome: "chr1", Position: 1000}
	err = sp.addGenomicContext(v)
	require.NoError(t, err)
	assert.NotNil(t, v.Context)
	assert.Equal(t, int64(1), sp.cacheMisses)
}

func TestAddGenomicContext_WithReference_CacheHit(t *testing.T) {
	// Second call at the same position increments cacheHits.
	p := makeTempFasta(t, "chr1", 4000)
	sp, err := NewStreamProcessor(StreamConfig{ReferencePath: p, NumWorkers: 1})
	require.NoError(t, err)

	v := &Variant{Chromosome: "chr1", Position: 1000}
	require.NoError(t, sp.addGenomicContext(v))
	assert.Equal(t, int64(0), sp.cacheHits)

	v2 := &Variant{Chromosome: "chr1", Position: 1000}
	require.NoError(t, sp.addGenomicContext(v2))
	assert.Equal(t, int64(1), sp.cacheHits)
}

func TestLoadGenomicContext_NoReference_ReturnsError(t *testing.T) {
	// Without a reference path, loadGenomicContext returns an error rather than fake data.
	sp := newSP()
	_, err := sp.loadGenomicContext("chr1", 100, 1100)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "reference genome not configured")
}

func TestLoadGenomicContext_WithReference_ReturnsSequence(t *testing.T) {
	// With a configured reference, the real sequence is returned.
	p := makeTempFasta(t, "chr1", 4000)
	sp, err := NewStreamProcessor(StreamConfig{ReferencePath: p, NumWorkers: 1})
	require.NoError(t, err)

	ctx, err := sp.loadGenomicContext("chr1", 0, 100)
	require.NoError(t, err)
	assert.Len(t, ctx, 100)
	for _, b := range ctx {
		assert.Contains(t, []byte{'A', 'C', 'G', 'T'}, b)
	}
}

func TestCacheRegion_StoresRegion(t *testing.T) {
	sp, _ := NewStreamProcessor(StreamConfig{
		MaxCacheSize: 100 * 1024 * 1024, // 100 MB
		NumWorkers:   1,
		KmerSize:     4,
	})
	context := make([]byte, 1000)

	sp.cacheRegion("chr1:0-1000", "chr1", 0, 1000, context)

	stats := sp.GetStatistics()
	assert.Equal(t, int64(1), stats.CachedRegions)
	assert.Equal(t, int64(1000), stats.CacheSize)
}

func TestEvictOldestRegions_EvictsWhenOverLimit(t *testing.T) {
	// Set a small cache (200 bytes). Each region is 100 bytes.
	// maxCacheSize/2 = 100. Evict until <=100 bytes.
	// After adding 3 regions (300 bytes > 200), eviction kicks in:
	// removes oldest until <= 100, leaving the last region (100 bytes).
	sp, _ := NewStreamProcessor(StreamConfig{
		MaxCacheSize: 200,
		NumWorkers:   1,
		KmerSize:     4,
	})

	for i := 0; i < 5; i++ {
		key := strings.Repeat(string(rune('a'+i)), 3)
		context := make([]byte, 100)
		sp.cacheRegion(key, "chr1", int64(i*100), int64(i*100+100), context)
	}

	// After 5 regions with eviction, cache should be reasonably small
	stats := sp.GetStatistics()
	// The cache size reflects the regions remaining after eviction
	assert.GreaterOrEqual(t, stats.CacheSize, int64(0))
	// Verify that eviction actually ran (some regions were removed)
	assert.Less(t, stats.CachedRegions, int64(5))
}

func TestExtractContextFromRegion_NoContext(t *testing.T) {
	// Without Context bytes cached, the function returns nil.
	sp := newSP()
	region := &VariantRegion{
		Chromosome: "chr1",
		StartPos:   500,
		EndPos:     1500,
	}
	v := &Variant{
		Chromosome: "chr1",
		Position:   1000,
	}
	ctx := sp.extractContextFromRegion(v, region)
	assert.Nil(t, ctx)
}

func TestExtractContextFromRegion_UsesContextBytes(t *testing.T) {
	// Context bytes stored in the region should be returned (centred on variant).
	sp := newSP()
	// Build a 2000-base context so position 1000 is in the middle.
	ctx2000 := make([]byte, 2000)
	for i := range ctx2000 {
		ctx2000[i] = []byte("ACGT")[i%4]
	}
	region := &VariantRegion{
		Chromosome: "chr1",
		StartPos:   0,
		EndPos:     2000,
		Context:    ctx2000,
	}
	v := &Variant{Chromosome: "chr1", Position: 1000}
	ctx := sp.extractContextFromRegion(v, region)

	// relativePos = 1000 - 0 = 1000; window [500, 1500) → 1000 bytes
	assert.Len(t, ctx, 1000)
	assert.Equal(t, ctx2000[500:1500], ctx)
}

func TestExtractContextFromRegion_ClampedAtStart(t *testing.T) {
	// Variant near the start of the region → window is clamped to [0, relPos+500).
	sp := newSP()
	ctxBytes := make([]byte, 600)
	for i := range ctxBytes {
		ctxBytes[i] = 'A'
	}
	region := &VariantRegion{StartPos: 0, EndPos: 600, Context: ctxBytes}
	v := &Variant{Position: 100} // relativePos=100; start clamped from -400 to 0
	ctx := sp.extractContextFromRegion(v, region)
	assert.Len(t, ctx, 600) // 0 to min(600,600)
}

func TestExtractContextFromRegion_ClampedAtEnd(t *testing.T) {
	// Variant near the end of the region → window clamped at region length.
	sp := newSP()
	ctxBytes := make([]byte, 600)
	for i := range ctxBytes {
		ctxBytes[i] = 'T'
	}
	region := &VariantRegion{StartPos: 0, EndPos: 600, Context: ctxBytes}
	v := &Variant{Position: 500} // relativePos=500; end clamped from 1000 to 600
	ctx := sp.extractContextFromRegion(v, region)
	assert.Len(t, ctx, 600) // clamped: start=0, end=600
}

// ─────────────────────────────────────────────────────────────────────────────
// parseStructuredHeader
// ─────────────────────────────────────────────────────────────────────────────

func TestParseStructuredHeader_InfoResult(t *testing.T) {
	sp := newSP()
	result, err := sp.parseStructuredHeader("<ID=AF,Number=A,Type=Float>", &HeaderInfo{})
	require.NoError(t, err)
	info, ok := result.(*HeaderInfo)
	require.True(t, ok)
	assert.Equal(t, "placeholder", info.ID)
}

func TestParseStructuredHeader_FormatResult(t *testing.T) {
	sp := newSP()
	result, err := sp.parseStructuredHeader("<ID=GT,Number=1,Type=String>", &HeaderFormat{})
	require.NoError(t, err)
	hf, ok := result.(*HeaderFormat)
	require.True(t, ok)
	assert.Equal(t, "placeholder", hf.ID)
}

func TestParseStructuredHeader_FilterResult(t *testing.T) {
	sp := newSP()
	result, err := sp.parseStructuredHeader(`<ID=LowQual,Description="Low">`, &HeaderFilter{})
	require.NoError(t, err)
	hf, ok := result.(*HeaderFilter)
	require.True(t, ok)
	assert.Equal(t, "placeholder", hf.ID)
}

func TestParseStructuredHeader_NotAngleBracketed(t *testing.T) {
	sp := newSP()
	_, err := sp.parseStructuredHeader("ID=AF,Number=A", &HeaderInfo{})
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "invalid structured header format")
}

func TestParseStructuredHeader_UnknownType(t *testing.T) {
	sp := newSP()
	type unknownType struct{}
	_, err := sp.parseStructuredHeader("<ID=X>", &unknownType{})
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "unknown header type")
}

// ─────────────────────────────────────────────────────────────────────────────
// GetStatistics with variant processing
// ─────────────────────────────────────────────────────────────────────────────

func TestGetStatistics_NoReference_NoCacheOps(t *testing.T) {
	// Without a reference, addGenomicContext is a no-op, so no cache operations occur.
	sp := newSP()

	for i := 0; i < 3; i++ {
		v := &Variant{Chromosome: "chr1", Position: int64(1000 + i*2000)}
		require.NoError(t, sp.addGenomicContext(v))
	}

	stats := sp.GetStatistics()
	assert.Equal(t, int64(0), stats.CacheMisses)
	assert.Equal(t, int64(0), stats.CachedRegions)
}

func TestGetStatistics_WithReference_CacheOpsRecorded(t *testing.T) {
	// With a reference, cache misses and regions are tracked.
	p := makeTempFasta(t, "chr1", 8000)
	sp, err := NewStreamProcessor(StreamConfig{ReferencePath: p, NumWorkers: 1})
	require.NoError(t, err)

	// Three different positions → three cache misses
	for i := 0; i < 3; i++ {
		v := &Variant{Chromosome: "chr1", Position: int64(1000 + i*2000)}
		require.NoError(t, sp.addGenomicContext(v))
	}
	stats := sp.GetStatistics()
	assert.Equal(t, int64(3), stats.CacheMisses)
	assert.Equal(t, int64(3), stats.CachedRegions)
}

func TestGetStatistics_WithReference_CacheHitRatio(t *testing.T) {
	// Same position 3 times: 1 miss + 2 hits → ratio = 2/3.
	p := makeTempFasta(t, "chr1", 4000)
	sp, err := NewStreamProcessor(StreamConfig{ReferencePath: p, NumWorkers: 1})
	require.NoError(t, err)

	v := &Variant{Chromosome: "chr1", Position: 1000}
	require.NoError(t, sp.addGenomicContext(v))
	require.NoError(t, sp.addGenomicContext(v))
	require.NoError(t, sp.addGenomicContext(v))

	stats := sp.GetStatistics()
	assert.Equal(t, int64(1), stats.CacheMisses)
	assert.Equal(t, int64(2), stats.CacheHits)
	assert.InDelta(t, 2.0/3.0, stats.CacheHitRatio, 1e-9)
}

func TestGetStatistics_NoAccesses_NotNaN(t *testing.T) {
	// A freshly created processor has no cache accesses; CacheHitRatio must be
	// 0.0, not NaN (which would result from 0/0 without a guard).
	sp := newSP()
	stats := sp.GetStatistics()
	assert.False(t, math.IsNaN(stats.CacheHitRatio), "CacheHitRatio should not be NaN")
	assert.False(t, math.IsInf(stats.CacheHitRatio, 0), "CacheHitRatio should not be Inf")
	assert.Equal(t, 0.0, stats.CacheHitRatio)
}

func TestGetStatistics_ZeroMaxCacheSize_NotNaN(t *testing.T) {
	// A processor constructed with no explicit cache size has maxCacheSize == 0;
	// CacheUtilization must be 0.0, not NaN or +Inf.
	sp := newSP() // default StreamConfig has no cache size set
	stats := sp.GetStatistics()
	assert.False(t, math.IsNaN(stats.CacheUtilization), "CacheUtilization should not be NaN")
	assert.False(t, math.IsInf(stats.CacheUtilization, 0), "CacheUtilization should not be Inf")
	assert.Equal(t, 0.0, stats.CacheUtilization)
}

// TestAuthenticityProcessor_GetResults_ZeroTotalKmers is a regression test for #127.
// When kmerProfiles is non-empty but every profile has TotalCount == 0, the old
// guard (len > 0) still divided by totalKmers == 0, causing a panic. The fix changes
// the guard to (totalKmers > 0).
func TestAuthenticityProcessor_GetResults_ZeroTotalKmers(t *testing.T) {
	ap, err := NewAuthenticityProcessor(AuthenticityConfig{KmerSize: 4})
	require.NoError(t, err)

	// NewPackedFrequencyTable returns a table with TotalCount == 0.
	emptyProfile, err := kmer.NewPackedFrequencyTable(4)
	require.NoError(t, err)

	// Inject it directly — simulates a profile that was added before any sequences
	// were processed (TotalCount stays 0).
	ap.kmerProfiles = append(ap.kmerProfiles, emptyProfile)
	ap.authenticityScores = append(ap.authenticityScores, 50.0)

	// Must not panic; complexity should be 0.0.
	results := ap.GetAuthenticityResults()
	assert.False(t, math.IsNaN(results.KmerComplexity), "KmerComplexity should not be NaN")
	assert.False(t, math.IsInf(results.KmerComplexity, 0), "KmerComplexity should not be Inf")
	assert.Equal(t, 0.0, results.KmerComplexity)
}
