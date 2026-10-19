package reference_test

import (
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

// smallFasta is a compact FASTA string used across many tests.
const smallFasta = ">chr1\nATCGATCGATCGATCGATCG\n>chr2\nGGCCGGCCGGCC\n>chrX\nNNNNATCGATCGNNNN\n"

// writeFasta writes content to a temp file and returns its path.
func writeFasta(t *testing.T, content string) string {
	t.Helper()
	f, err := os.CreateTemp(t.TempDir(), "test_*.fasta")
	require.NoError(t, err)
	_, err = f.WriteString(content)
	require.NoError(t, err)
	require.NoError(t, f.Close())
	return f.Name()
}

// ---------------------------------------------------------------------------
// FastaParser
// ---------------------------------------------------------------------------

func TestNewFastaParser_Nonexistent(t *testing.T) {
	_, err := reference.NewFastaParser("/no/such/file.fasta")
	assert.Error(t, err)
}

func TestNewFastaParser_ValidFile(t *testing.T) {
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)
	assert.NotNil(t, fp)
}

func TestFastaParser_OpenClose(t *testing.T) {
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)

	// Open twice – second call should be no-op
	require.NoError(t, fp.Open())
	require.NoError(t, fp.Open())

	// Close twice – second call should be no-op
	require.NoError(t, fp.Close())
	require.NoError(t, fp.Close())
}

func TestFastaParser_Open_BadFile(t *testing.T) {
	// Create a parser for a file we will delete
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)
	require.NoError(t, os.Remove(path))
	err = fp.Open()
	assert.Error(t, err)
}

func TestFastaParser_BuildIndex(t *testing.T) {
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)

	index, err := fp.BuildIndex()
	require.NoError(t, err)
	require.NotNil(t, index)

	assert.ElementsMatch(t, []string{"chr1", "chr2", "chrX"}, index.OrderedNames)
	assert.Contains(t, index.Chromosomes, "chr1")
	assert.Contains(t, index.Chromosomes, "chr2")
	assert.Contains(t, index.Chromosomes, "chrX")

	chr1 := index.Chromosomes["chr1"]
	assert.Greater(t, chr1.Length, int64(0))
	assert.Equal(t, "chr1", chr1.Name)
}

func TestFastaParser_BuildIndex_MultiLineSequence(t *testing.T) {
	multiLine := ">chr1\nATCGATCG\nATCGATCG\nATCGATCG\n>chr2\nGGCC\n"
	path := writeFasta(t, multiLine)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)

	index, err := fp.BuildIndex()
	require.NoError(t, err)
	assert.Equal(t, int64(24), index.Chromosomes["chr1"].Length)
	assert.Equal(t, int64(4), index.Chromosomes["chr2"].Length)
}

func TestFastaParser_LoadIndex_NoFaiFile_BuildsFromFasta(t *testing.T) {
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)

	// No .fai exists – LoadIndex should build from FASTA
	index, err := fp.LoadIndex()
	require.NoError(t, err)
	assert.NotNil(t, index)
	assert.Contains(t, index.Chromosomes, "chr1")
}

func TestFastaParser_SaveAndLoadIndex(t *testing.T) {
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)

	// Build and save
	index, err := fp.BuildIndex()
	require.NoError(t, err)
	require.NoError(t, fp.SaveIndex(index))

	// Reload from .fai
	fp2, err := reference.NewFastaParser(path)
	require.NoError(t, err)
	loaded, err := fp2.LoadIndex()
	require.NoError(t, err)

	assert.Len(t, loaded.Chromosomes, 3)
	for name, orig := range index.Chromosomes {
		reloaded, ok := loaded.Chromosomes[name]
		require.True(t, ok, "chromosome %s should be in reloaded index", name)
		assert.Equal(t, orig.Length, reloaded.Length)
	}
}

func TestFastaParser_GetChromosomeNames(t *testing.T) {
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)

	index, err := fp.BuildIndex()
	require.NoError(t, err)

	names := fp.GetChromosomeNames(index)
	assert.ElementsMatch(t, []string{"chr1", "chr2", "chrX"}, names)
}

func TestFastaParser_GetChromosomeLength(t *testing.T) {
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)

	index, err := fp.BuildIndex()
	require.NoError(t, err)

	t.Run("existing chromosome", func(t *testing.T) {
		length, err := fp.GetChromosomeLength(index, "chr1")
		require.NoError(t, err)
		assert.Equal(t, int64(20), length)
	})

	t.Run("missing chromosome", func(t *testing.T) {
		_, err := fp.GetChromosomeLength(index, "chrMissing")
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "not found")
	})
}

func TestFastaParser_GetSequence_InvalidCoords(t *testing.T) {
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)

	index, err := fp.BuildIndex()
	require.NoError(t, err)

	t.Run("negative start", func(t *testing.T) {
		_, err := fp.GetSequence(index, "chr1", -1, 5)
		assert.Error(t, err)
	})

	t.Run("start >= end", func(t *testing.T) {
		_, err := fp.GetSequence(index, "chr1", 10, 5)
		assert.Error(t, err)
	})

	t.Run("end beyond chromosome", func(t *testing.T) {
		_, err := fp.GetSequence(index, "chr1", 0, 9999)
		assert.Error(t, err)
	})

	t.Run("missing chromosome", func(t *testing.T) {
		_, err := fp.GetSequence(index, "chrMissing", 0, 5)
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "not found")
	})
}

func TestFastaParser_GetSequence_ValidRegion(t *testing.T) {
	// Multi-line FASTA with known content so we can verify extracted sequence
	fasta := ">chr1\nATCGATCGATCGATCGATCG\n"
	path := writeFasta(t, fasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)

	index, err := fp.BuildIndex()
	require.NoError(t, err)

	seq, err := fp.GetSequence(index, "chr1", 0, 4)
	require.NoError(t, err)
	assert.Equal(t, "ATCG", string(seq))
}

// ---------------------------------------------------------------------------
// KmerExtractor
// ---------------------------------------------------------------------------

func TestNewKmerExtractor(t *testing.T) {
	ext := reference.NewKmerExtractor(4)
	require.NotNil(t, ext)
}

func TestKmerExtractor_SetWindowSize(t *testing.T) {
	ext := reference.NewKmerExtractor(3)
	ext.SetWindowSize(500) // Should not panic
}

func TestKmerExtractor_SetSkipNRegions(t *testing.T) {
	ext := reference.NewKmerExtractor(3)
	ext.SetSkipNRegions(false)
	ext.SetSkipNRegions(true)
}

func TestKmerExtractor_ExtractKmerFrequencies_TooShort(t *testing.T) {
	ext := reference.NewKmerExtractor(10)
	err := ext.ExtractKmerFrequencies([]byte("ATCG"))
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "less than k")
}

func TestKmerExtractor_ExtractKmerFrequencies_Basic(t *testing.T) {
	ext := reference.NewKmerExtractor(3)
	seq := []byte("ATCGATCGATCGATCGATCG")
	err := ext.ExtractKmerFrequencies(seq)
	require.NoError(t, err)

	table := ext.GetFrequencyTable()
	require.NotNil(t, table)
	assert.Equal(t, 3, table.K)
	assert.Greater(t, table.TotalKmers, uint64(0))
	assert.Greater(t, len(table.Frequencies), 0)

	// Frequencies should sum to ~1.0
	sum := 0.0
	for _, f := range table.Frequencies {
		sum += f
	}
	assert.InDelta(t, 1.0, sum, 0.01)
}

func TestKmerExtractor_ExtractKmerFrequencies_SkipN(t *testing.T) {
	ext := reference.NewKmerExtractor(3)
	ext.SetSkipNRegions(true)
	seq := []byte("ATCGNNN") // last 3 have N
	err := ext.ExtractKmerFrequencies(seq)
	require.NoError(t, err)

	// Without Ns we should have fewer kmers
	table := ext.GetFrequencyTable()
	for kmer := range table.Frequencies {
		assert.NotContains(t, kmer, "N", "no N-containing kmers should appear")
	}
}

func TestKmerExtractor_ExtractKmerFrequencies_NoSkipN(t *testing.T) {
	ext := reference.NewKmerExtractor(3)
	ext.SetSkipNRegions(false)
	// All N - isValidKmer will reject because N is not A/T/G/C
	seq := []byte("NNNATCGATCG")
	err := ext.ExtractKmerFrequencies(seq)
	require.NoError(t, err)
	// NNN is invalid, NNAT etc. are invalid – only pure ATCG kmers counted
	table := ext.GetFrequencyTable()
	for kmer := range table.Frequencies {
		for _, base := range kmer {
			assert.Contains(t, "ATGC", string(base))
		}
	}
}

func TestKmerExtractor_ExtractKmerFrequenciesFromChromosome(t *testing.T) {
	ext := reference.NewKmerExtractor(4)
	chr := &reference.Chromosome{
		Name:     "chr1",
		Sequence: []byte("ATCGATCGATCGATCGATCG"),
		Length:   20,
	}
	err := ext.ExtractKmerFrequenciesFromChromosome(chr)
	require.NoError(t, err)

	table := ext.GetFrequencyTable()
	assert.Greater(t, table.TotalKmers, uint64(0))
}

func TestKmerExtractor_GetObservedKmers(t *testing.T) {
	ext := reference.NewKmerExtractor(2)
	err := ext.ExtractKmerFrequencies([]byte("ATCGATCGATCG"))
	require.NoError(t, err)
	observed := ext.GetObservedKmers()
	assert.Greater(t, observed, 0)
}

func TestKmerExtractor_ComputeExpectedKmers(t *testing.T) {
	ext := reference.NewKmerExtractor(1)
	assert.Equal(t, 4, ext.ComputeExpectedKmers(1))
	assert.Equal(t, 16, ext.ComputeExpectedKmers(2))
	assert.Equal(t, 64, ext.ComputeExpectedKmers(3))
}

func TestKmerExtractor_GenerateAllPossibleKmers(t *testing.T) {
	ext := reference.NewKmerExtractor(1)

	kmers1 := ext.GenerateAllPossibleKmers(1)
	assert.Len(t, kmers1, 4)

	kmers2 := ext.GenerateAllPossibleKmers(2)
	assert.Len(t, kmers2, 16)

	kmers3 := ext.GenerateAllPossibleKmers(3)
	assert.Len(t, kmers3, 64)
}

func TestKmerExtractor_NormalizeFrequencies(t *testing.T) {
	ext := reference.NewKmerExtractor(2)

	freqs := map[string]float64{"AT": 3.0, "CG": 1.0, "GC": 2.0}
	normalized := ext.NormalizeFrequencies(freqs)

	sum := 0.0
	for _, v := range normalized {
		sum += v
	}
	assert.InDelta(t, 1.0, sum, 1e-9)
}

func TestKmerExtractor_NormalizeFrequencies_ZeroSum(t *testing.T) {
	ext := reference.NewKmerExtractor(2)
	freqs := map[string]float64{"AT": 0.0}
	normalized := ext.NormalizeFrequencies(freqs)
	// Should return the same map unchanged
	assert.Equal(t, 0.0, normalized["AT"])
}

func TestKmerExtractor_GetKmerComplexity_Empty(t *testing.T) {
	ext := reference.NewKmerExtractor(3)
	// Before any extraction complexity should be 0
	complexity := ext.GetKmerComplexity()
	assert.Equal(t, 0.0, complexity)
}

func TestKmerExtractor_GetKmerComplexity_AfterExtraction(t *testing.T) {
	ext := reference.NewKmerExtractor(2)
	err := ext.ExtractKmerFrequencies([]byte("ATCGATCGATCGATCGATCG"))
	require.NoError(t, err)

	complexity := ext.GetKmerComplexity()
	assert.GreaterOrEqual(t, complexity, 0.0)
	assert.LessOrEqual(t, complexity, 1.0)
}

func TestKmerExtractor_GetTopKmers(t *testing.T) {
	ext := reference.NewKmerExtractor(2)
	err := ext.ExtractKmerFrequencies([]byte("ATCGATCGATCGATCGATCG"))
	require.NoError(t, err)

	topKmers := ext.GetTopKmers(3)
	assert.LessOrEqual(t, len(topKmers), 3)

	// Verify sorted descending
	for i := 1; i < len(topKmers); i++ {
		assert.GreaterOrEqual(t, topKmers[i-1].Frequency, topKmers[i].Frequency)
	}
}

func TestKmerExtractor_GetTopKmers_MoreThanAvailable(t *testing.T) {
	ext := reference.NewKmerExtractor(2)
	err := ext.ExtractKmerFrequencies([]byte("ATCG"))
	require.NoError(t, err)

	// Requesting 100 but fewer unique kmers exist
	topKmers := ext.GetTopKmers(100)
	assert.LessOrEqual(t, len(topKmers), 100)
}

func TestKmerExtractor_GetKmerStats(t *testing.T) {
	ext := reference.NewKmerExtractor(3)
	err := ext.ExtractKmerFrequencies([]byte("ATCGATCGATCGATCG"))
	require.NoError(t, err)

	stats := ext.GetKmerStats()
	assert.Equal(t, 3, stats["k"])
	assert.Greater(t, stats["total_kmers"], uint64(0))
	assert.GreaterOrEqual(t, stats["unique_kmers"], 0)
	assert.Equal(t, 64, stats["expected_kmers"])
}

func TestKmerExtractor_MergeFrequencyTables_Empty(t *testing.T) {
	ext := reference.NewKmerExtractor(3)
	_, err := ext.MergeFrequencyTables(nil)
	assert.Error(t, err)
}

func TestKmerExtractor_MergeFrequencyTables_DifferentK(t *testing.T) {
	ext := reference.NewKmerExtractor(3)

	t1 := &reference.KmerFrequencyTable{K: 3, Frequencies: map[string]float64{"ATG": 1.0}, TotalKmers: 1}
	t2 := &reference.KmerFrequencyTable{K: 4, Frequencies: map[string]float64{"ATCG": 1.0}, TotalKmers: 1}

	_, err := ext.MergeFrequencyTables([]*reference.KmerFrequencyTable{t1, t2})
	assert.Error(t, err)
}

func TestKmerExtractor_MergeFrequencyTables_Valid(t *testing.T) {
	ext := reference.NewKmerExtractor(3)

	t1 := &reference.KmerFrequencyTable{K: 3, Frequencies: map[string]float64{"ATG": 0.5, "CGA": 0.5}, TotalKmers: 10}
	t2 := &reference.KmerFrequencyTable{K: 3, Frequencies: map[string]float64{"ATG": 0.4, "TCG": 0.6}, TotalKmers: 10}

	merged, err := ext.MergeFrequencyTables([]*reference.KmerFrequencyTable{t1, t2})
	require.NoError(t, err)
	assert.Equal(t, 3, merged.K)
	assert.Greater(t, merged.TotalKmers, uint64(0))
	assert.Contains(t, merged.Frequencies, "ATG")
}

func TestKmerExtractor_CompareFrequencies_DifferentK(t *testing.T) {
	ext := reference.NewKmerExtractor(3)

	t1 := &reference.KmerFrequencyTable{K: 3, Frequencies: map[string]float64{"ATG": 1.0}}
	t2 := &reference.KmerFrequencyTable{K: 4, Frequencies: map[string]float64{"ATCG": 1.0}}

	result := ext.CompareFrequencies(t1, t2)
	assert.Nil(t, result)
}

func TestKmerExtractor_CompareFrequencies_Valid(t *testing.T) {
	ext := reference.NewKmerExtractor(3)

	t1 := &reference.KmerFrequencyTable{K: 3, Frequencies: map[string]float64{"ATG": 0.6, "CGA": 0.4}}
	t2 := &reference.KmerFrequencyTable{K: 3, Frequencies: map[string]float64{"ATG": 0.3, "TCG": 0.7}}

	diffs := ext.CompareFrequencies(t1, t2)
	require.NotNil(t, diffs)
	// ATG diff = 0.3 - 0.6 = -0.3
	assert.InDelta(t, -0.3, diffs["ATG"], 1e-9)
	// TCG diff = 0.7 - 0.0 = 0.7
	assert.InDelta(t, 0.7, diffs["TCG"], 1e-9)
	// CGA diff = 0.0 - 0.4 = -0.4
	assert.InDelta(t, -0.4, diffs["CGA"], 1e-9)
}

// ---------------------------------------------------------------------------
// ReferenceGenomeParser (genome_parser.go – functions not in genome_parser_test.go)
// ---------------------------------------------------------------------------

func TestReferenceGenomeParser_LoadChromosome_WithSequenceDetails(t *testing.T) {
	// Use in-memory FASTA with known GC content
	// chr1: 50% GC, chr2: all GC (100%), chrN: all N
	fasta := ">chr1\nATGCATGCATGC\n>chr2\nGGGGCCCC\n>chrN\nNNNNNNNNNNNN\n"
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "g.fasta")
	require.NoError(t, os.WriteFile(fastaPath, []byte(fasta), 0644))

	parser := reference.NewReferenceGenomeParser(fastaPath)
	require.NoError(t, parser.CreateIndex())
	require.NoError(t, parser.LoadIndex())

	chr1, err := parser.LoadChromosome("chr1")
	require.NoError(t, err)
	assert.InDelta(t, 0.5, chr1.GCContent, 0.01)
	assert.Equal(t, 0, len(chr1.NRegions))

	chr2, err := parser.LoadChromosome("chr2")
	require.NoError(t, err)
	assert.InDelta(t, 1.0, chr2.GCContent, 0.01)

	chrN, err := parser.LoadChromosome("chrN")
	require.NoError(t, err)
	assert.Greater(t, len(chrN.NRegions), 0)
}

// ---------------------------------------------------------------------------
// SimpleReferenceGenome / SimpleChromosome  (simple_genome.go)
// ---------------------------------------------------------------------------

func TestNewSimpleReferenceGenome(t *testing.T) {
	g := reference.NewSimpleReferenceGenome("/some/dir")
	require.NotNil(t, g)
}

func TestNewReferenceGenome_Alias(t *testing.T) {
	g := reference.NewReferenceGenome("/some/dir")
	require.NotNil(t, g)
}

func TestSimpleReferenceGenome_GetChromosomeNames_EmptyDir(t *testing.T) {
	dir := t.TempDir()
	g := reference.NewSimpleReferenceGenome(dir)
	names := g.GetChromosomeNames()
	assert.Empty(t, names)
}

func TestSimpleReferenceGenome_GetChromosomeNames_BadDir(t *testing.T) {
	g := reference.NewSimpleReferenceGenome("/no/such/directory")
	names := g.GetChromosomeNames()
	assert.Nil(t, names)
}

func TestSimpleReferenceGenome_GetChromosomeNames_WithFiles(t *testing.T) {
	dir := t.TempDir()
	// Create some chromosome FASTA files
	for _, name := range []string{"chr1.fasta", "chr2.fasta", "chrX.fasta"} {
		require.NoError(t, os.WriteFile(filepath.Join(dir, name), []byte(">chr\nATCG\n"), 0644))
	}
	// Non-matching files should be ignored
	require.NoError(t, os.WriteFile(filepath.Join(dir, "README.txt"), []byte("ignored"), 0644))
	require.NoError(t, os.WriteFile(filepath.Join(dir, "notchr.fasta"), []byte(">notchr\nATCG\n"), 0644))

	g := reference.NewSimpleReferenceGenome(dir)
	names := g.GetChromosomeNames()
	assert.ElementsMatch(t, []string{"chr1", "chr2", "chrX"}, names)
}

func TestSimpleReferenceGenome_LoadChromosome_MissingFile(t *testing.T) {
	dir := t.TempDir()
	g := reference.NewSimpleReferenceGenome(dir)
	_, err := g.LoadChromosome("chr1")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "not found")
}

func TestSimpleReferenceGenome_LoadChromosome_FromFile(t *testing.T) {
	dir := t.TempDir()
	fasta := ">chr1\nATCGATCGATCG\n"
	require.NoError(t, os.WriteFile(filepath.Join(dir, "chr1.fasta"), []byte(fasta), 0644))

	g := reference.NewSimpleReferenceGenome(dir)
	chr, err := g.LoadChromosome("chr1")
	require.NoError(t, err)
	require.NotNil(t, chr)
	assert.Equal(t, "chr1", chr.GetName())
	assert.Equal(t, 12, chr.GetLength())
	assert.Equal(t, []byte("ATCGATCGATCG"), chr.GetSequence())
}

func TestSimpleReferenceGenome_LoadChromosome_FromGRCh38Fallback(t *testing.T) {
	dir := t.TempDir()
	// No chr22.fasta, but GRCh38.fasta exists and contains chr22
	fasta := ">chr22\nGGCCATCG\n>chr1\nATCG\n"
	require.NoError(t, os.WriteFile(filepath.Join(dir, "GRCh38.fasta"), []byte(fasta), 0644))

	g := reference.NewSimpleReferenceGenome(dir)
	chr, err := g.LoadChromosome("chr22")
	require.NoError(t, err)
	assert.Equal(t, "chr22", chr.GetName())
}

func TestSimpleReferenceGenome_LoadChromosome_MultipleInFile(t *testing.T) {
	dir := t.TempDir()
	fasta := ">chr1\nATCGATCG\n>chr2\nGGCC\n"
	require.NoError(t, os.WriteFile(filepath.Join(dir, "chr1.fasta"), []byte(fasta), 0644))

	g := reference.NewSimpleReferenceGenome(dir)
	// Should stop reading when the next chromosome header is found
	chr, err := g.LoadChromosome("chr1")
	require.NoError(t, err)
	assert.Equal(t, 8, chr.GetLength())
}

// TestSimpleReferenceGenome_LoadChromosome_ExactHeaderMatch: a record is
// selected by the first whitespace-delimited token of its header, so "chr1"
// must not also collect chr10 and chr11.
func TestSimpleReferenceGenome_LoadChromosome_ExactHeaderMatch(t *testing.T) {
	dir := t.TempDir()
	fasta := ">chr10\nGGGG\n>chr1 dna:chromosome chromosome:GRCh38:1\nAC\nGT\n>chr11\nCCCC\n"
	require.NoError(t, os.WriteFile(filepath.Join(dir, "GRCh38.fasta"), []byte(fasta), 0644))

	g := reference.NewSimpleReferenceGenome(dir)
	for name, want := range map[string]string{"chr1": "ACGT", "chr10": "GGGG", "chr11": "CCCC"} {
		chr, err := g.LoadChromosome(name)
		require.NoError(t, err, name)
		assert.Equal(t, want, string(chr.GetSequence()), name)
	}
	_, err := g.LoadChromosome("chr")
	assert.Error(t, err, "a prefix of a record name must not match")
}

func TestSimpleReferenceGenome_LoadChromosome_NotFoundInFile(t *testing.T) {
	dir := t.TempDir()
	fasta := ">chr1\nATCG\n"
	require.NoError(t, os.WriteFile(filepath.Join(dir, "chr99.fasta"), []byte(fasta), 0644))

	g := reference.NewSimpleReferenceGenome(dir)
	_, err := g.LoadChromosome("chr99")
	assert.Error(t, err)
}

// ---------------------------------------------------------------------------
// LazyReferenceLoader
// ---------------------------------------------------------------------------

func makeDummyFastaFile(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "ref.fasta")
	require.NoError(t, os.WriteFile(p, []byte(">chr1\nATCG\n"), 0644))
	return p
}

func TestNewLazyReferenceLoader_EmptyFile(t *testing.T) {
	_, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{})
	assert.Error(t, err)
}

func TestNewLazyReferenceLoader_MissingFile(t *testing.T) {
	_, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile: "/no/such/file.fasta",
	})
	assert.Error(t, err)
}

func TestNewLazyReferenceLoader_ValidFile(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile: p,
	})
	require.NoError(t, err)
	assert.NotNil(t, loader)
}

func TestNewLazyReferenceLoader_Defaults(t *testing.T) {
	p := makeDummyFastaFile(t)
	// Negative values for cache / compression should be corrected to defaults
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile:    p,
		MaxCacheSize:     -1,
		CompressionLevel: -5,
		KmerSize:         0,
	})
	require.NoError(t, err)
	assert.NotNil(t, loader)
}

func TestLazyReferenceLoader_GetSequence_Basic(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile: p,
	})
	require.NoError(t, err)

	seq, err := loader.GetSequence("chr1", 0, 10)
	require.NoError(t, err)
	assert.Len(t, seq, 10)
}

func TestLazyReferenceLoader_GetSequence_InvalidLength(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: p})
	require.NoError(t, err)

	// Zero-length request
	_, err = loader.GetSequence("chr1", 5, 5)
	assert.Error(t, err)
}

func TestLazyReferenceLoader_GetSequence_CacheHit(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: p})
	require.NoError(t, err)

	// First call populates cache
	seq1, err := loader.GetSequence("chr1", 0, 10)
	require.NoError(t, err)

	// Second call should hit cache
	seq2, err := loader.GetSequence("chr1", 0, 10)
	require.NoError(t, err)

	assert.Equal(t, seq1, seq2)

	stats := loader.GetCacheStatistics()
	assert.Greater(t, stats.CacheHits, int64(0))
}

func TestLazyReferenceLoader_GetCacheStatistics(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: p})
	require.NoError(t, err)

	_, _ = loader.GetSequence("chr1", 0, 100)
	_, _ = loader.GetSequence("chr1", 0, 100)
	_, _ = loader.GetSequence("chr1", 100, 200)

	stats := loader.GetCacheStatistics()
	assert.GreaterOrEqual(t, stats.CacheHits, int64(1))
	assert.GreaterOrEqual(t, stats.CacheMisses, int64(1))
	assert.GreaterOrEqual(t, stats.HitRatio, 0.0)
	assert.LessOrEqual(t, stats.HitRatio, 1.0)
	assert.GreaterOrEqual(t, stats.CachedRegions, int64(0))
}

// TestLazyReferenceLoader_GetCacheStatistics_ZeroMisses verifies that
// AvgLoadTime is 0.0 (not +Inf or NaN) when no cache misses have occurred.
// This is the regression test for the divide-by-zero fix (issue #75).
func TestLazyReferenceLoader_GetCacheStatistics_ZeroMisses(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: p})
	require.NoError(t, err)

	// Fresh loader: no requests yet → cacheMisses == 0.
	stats := loader.GetCacheStatistics()
	if math.IsInf(stats.AvgLoadTime, 0) || math.IsNaN(stats.AvgLoadTime) {
		t.Errorf("AvgLoadTime = %v with zero misses; want 0.0 (no +Inf/NaN)", stats.AvgLoadTime)
	}
	assert.Equal(t, 0.0, stats.AvgLoadTime)

	// Warm the cache, then repeat: all subsequent calls are hits, still zero misses.
	regions := []reference.GenomicRegion{{Chromosome: "chr1", Start: 0, End: 50}}
	_ = loader.WarmupCache(regions, 0)
	_, _ = loader.GetSequence("chr1", 0, 50) // hits warmed region

	stats2 := loader.GetCacheStatistics()
	if math.IsInf(stats2.AvgLoadTime, 0) || math.IsNaN(stats2.AvgLoadTime) {
		t.Errorf("AvgLoadTime = %v after cache-only hits; want 0.0", stats2.AvgLoadTime)
	}
}

// TestLazyReferenceLoader_GetCacheStatistics_ZeroMaxSize verifies that
// CacheUtilization is 0.0 (not +Inf/NaN) when MaxCacheSize is zero.
func TestLazyReferenceLoader_GetCacheStatistics_ZeroMaxSize(t *testing.T) {
	p := makeDummyFastaFile(t)
	// LazyLoaderConfig with MaxCacheSize == 0 should produce safe statistics.
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile: p,
		MaxCacheSize:  0,
	})
	require.NoError(t, err)

	stats := loader.GetCacheStatistics()
	if math.IsInf(stats.CacheUtilization, 0) || math.IsNaN(stats.CacheUtilization) {
		t.Errorf("CacheUtilization = %v with zero MaxCacheSize; want 0.0", stats.CacheUtilization)
	}
	assert.Equal(t, 0.0, stats.CacheUtilization)
}

func TestLazyReferenceLoader_ClearCache(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: p})
	require.NoError(t, err)

	_, _ = loader.GetSequence("chr1", 0, 100)
	loader.ClearCache()

	stats := loader.GetCacheStatistics()
	assert.Equal(t, int64(0), stats.CachedRegions)
	assert.Equal(t, int64(0), stats.CacheHits)
}

func TestLazyReferenceLoader_Close(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: p})
	require.NoError(t, err)

	_, _ = loader.GetSequence("chr1", 0, 50)
	assert.NoError(t, loader.Close())
}

func TestLazyReferenceLoader_WarmupCache(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: p})
	require.NoError(t, err)

	regions := []reference.GenomicRegion{
		{Chromosome: "chr1", Start: 0, End: 50},
		{Chromosome: "chr1", Start: 50, End: 100},
	}
	err = loader.WarmupCache(regions, 0)
	assert.NoError(t, err)
}

func TestLazyReferenceLoader_GetKmerProfile(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile:      p,
		EnablePackedArrays: false,
	})
	require.NoError(t, err)

	profile, err := loader.GetKmerProfile("chr1", 0, 100, 4)
	require.NoError(t, err)
	assert.NotNil(t, profile)
}

func TestLazyReferenceLoader_PreloadRegions_Invalid(t *testing.T) {
	p := makeDummyFastaFile(t)
	// Invalid region format
	_, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile:  p,
		PreloadRegions: []string{"bad_region_format"},
		KmerSize:       4,
	})
	assert.Error(t, err)
}

func TestLazyReferenceLoader_CacheEviction(t *testing.T) {
	p := makeDummyFastaFile(t)
	// Very small cache to force eviction (1 byte)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile: p,
		MaxCacheSize:  1,
	})
	require.NoError(t, err)

	// Load multiple distinct regions to trigger eviction
	for i := int64(0); i < 5; i++ {
		start := i * 100
		end := start + 100
		_, _ = loader.GetSequence("chr1", start, end)
	}
	// No panic is the success condition
}

// ---------------------------------------------------------------------------
// ReferenceManager
// ---------------------------------------------------------------------------

func TestNewReferenceManager_NoDefaultFile(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{})
	require.NoError(t, err)
	assert.NotNil(t, rm)
}

func TestNewReferenceManager_WithDefaultFile(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultReferenceFile: p,
		DefaultCacheSize:     1024 * 1024,
	})
	require.NoError(t, err)
	assert.NotNil(t, rm)
}

func TestNewReferenceManager_BadDefaultFile(t *testing.T) {
	_, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultReferenceFile: "/no/such/file.fasta",
	})
	assert.Error(t, err)
}

func TestNewReferenceManager_Defaults(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultCacheSize: -1,
		PrefetchWindow:   -1,
		MaxLoaders:       0,
		CompressionLevel: -5,
	})
	require.NoError(t, err)
	assert.NotNil(t, rm)
}

func TestReferenceManager_GetSequence_NoLoader(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{})
	require.NoError(t, err)

	_, err = rm.GetSequence("chr1", 0, 100)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "no reference loader available")
}

func TestReferenceManager_GetKmerProfile_NoLoader(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{})
	require.NoError(t, err)

	_, err = rm.GetKmerProfile("chr1", 0, 100, 4)
	assert.Error(t, err)
}

func TestReferenceManager_RegisterReferenceGenome(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{MaxLoaders: 5})
	require.NoError(t, err)

	err = rm.RegisterReferenceGenome("hg38", p, 0)
	assert.NoError(t, err)
}

func TestReferenceManager_RegisterReferenceGenome_BadFile(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{MaxLoaders: 5})
	require.NoError(t, err)

	err = rm.RegisterReferenceGenome("bad", "/no/such/file.fasta", 0)
	assert.Error(t, err)
}

func TestReferenceManager_RegisterReferenceGenome_ExceedsMaxLoaders(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{MaxLoaders: 1})
	require.NoError(t, err)

	require.NoError(t, rm.RegisterReferenceGenome("g1", p, 0))

	err = rm.RegisterReferenceGenome("g2", p, 0)
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "maximum number of loaders")
}

func TestReferenceManager_GetSequence_WithCustomLoader(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{MaxLoaders: 5})
	require.NoError(t, err)
	require.NoError(t, rm.RegisterReferenceGenome("hg38", p, 0))

	seq, err := rm.GetSequence("chr1", 0, 50, "hg38")
	require.NoError(t, err)
	assert.Len(t, seq, 50)
}

func TestReferenceManager_GetSequence_WithDefaultLoader(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultReferenceFile: p,
	})
	require.NoError(t, err)

	seq, err := rm.GetSequence("chr1", 0, 50)
	require.NoError(t, err)
	assert.Len(t, seq, 50)
}

func TestReferenceManager_GetSequence_WithPrefetching(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultReferenceFile: p,
		EnablePrefetching:    true,
		PrefetchWindow:       10,
	})
	require.NoError(t, err)

	seq, err := rm.GetSequence("chr1", 100, 150)
	require.NoError(t, err)
	assert.Len(t, seq, 50)
}

func TestReferenceManager_GetCombinedStatistics(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultReferenceFile: p,
	})
	require.NoError(t, err)

	_, _ = rm.GetSequence("chr1", 0, 100)
	_, _ = rm.GetSequence("chr1", 0, 100)

	stats := rm.GetCombinedStatistics()
	assert.GreaterOrEqual(t, stats.TotalRequests, int64(2))
	assert.Contains(t, stats.LoaderStats, "default")
}

func TestReferenceManager_WarmupReference_NoLoader(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{})
	require.NoError(t, err)

	err = rm.WarmupReference([]reference.GenomicRegion{{Chromosome: "chr1", Start: 0, End: 100}}, 0)
	assert.Error(t, err)
}

func TestReferenceManager_WarmupReference_WithLoader(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{DefaultReferenceFile: p})
	require.NoError(t, err)

	regions := []reference.GenomicRegion{{Chromosome: "chr1", Start: 0, End: 100}}
	err = rm.WarmupReference(regions, 0)
	assert.NoError(t, err)
}

func TestReferenceManager_OptimizeCache(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{DefaultReferenceFile: p})
	require.NoError(t, err)

	// Should not panic
	assert.NotPanics(t, func() { rm.OptimizeCache() })
}

func TestReferenceManager_ClearAllCaches(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{DefaultReferenceFile: p})
	require.NoError(t, err)
	_, _ = rm.GetSequence("chr1", 0, 100)

	assert.NotPanics(t, func() { rm.ClearAllCaches() })
}

func TestReferenceManager_ClearAllCaches_NoLoader(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{})
	require.NoError(t, err)
	// Should not panic when there's no default loader
	assert.NotPanics(t, func() { rm.ClearAllCaches() })
}

func TestReferenceManager_Close(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{DefaultReferenceFile: p})
	require.NoError(t, err)

	err = rm.Close()
	assert.NoError(t, err)
}

func TestReferenceManager_GetReferenceInfo(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultReferenceFile: p,
		MaxLoaders:           10,
	})
	require.NoError(t, err)
	require.NoError(t, rm.RegisterReferenceGenome("custom", p, 0))

	info := rm.GetReferenceInfo()
	assert.True(t, info.HasDefault)
	assert.Equal(t, 10, info.MaxLoaders)
	assert.Contains(t, info.CustomReferences, "custom")
	assert.Equal(t, 2, info.TotalLoaders) // default + custom
}

func TestReferenceManager_GetReferenceInfo_NoDefault(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{})
	require.NoError(t, err)
	info := rm.GetReferenceInfo()
	assert.False(t, info.HasDefault)
}

func TestReferenceManager_GetOptimalRegionSize(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{})
	require.NoError(t, err)

	size := rm.GetOptimalRegionSize()
	assert.Greater(t, size, int64(0))
}

func TestReferenceManager_PrefetchRegions_NoLoader(t *testing.T) {
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{})
	require.NoError(t, err)

	err = rm.PrefetchRegions([]reference.GenomicRegion{{Chromosome: "chr1", Start: 0, End: 100}}, 0)
	assert.Error(t, err)
}

func TestReferenceManager_PrefetchRegions_WithLoader(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{DefaultReferenceFile: p})
	require.NoError(t, err)

	regions := []reference.GenomicRegion{{Chromosome: "chr1", Start: 0, End: 100}}
	err = rm.PrefetchRegions(regions, 4)
	assert.NoError(t, err) // Prefetching is async
}

// ---------------------------------------------------------------------------
// Additional coverage for edge cases
// ---------------------------------------------------------------------------

func TestKmerFrequencyTable_ZeroTotalKmers(t *testing.T) {
	// GetFrequencyTable when no extraction has been done yet
	ext := reference.NewKmerExtractor(3)
	table := ext.GetFrequencyTable()
	assert.Equal(t, 3, table.K)
	assert.Equal(t, uint64(0), table.TotalKmers)
}

func TestFastaParser_SaveIndex_BadPath(t *testing.T) {
	// Build an index from a real file
	path := writeFasta(t, smallFasta)
	fp, err := reference.NewFastaParser(path)
	require.NoError(t, err)
	index, err := fp.BuildIndex()
	require.NoError(t, err)

	// Create a parser whose index path points to a non-existent subdirectory
	// so SaveIndex (which calls os.Create) will fail.
	// We give it an existing FASTA file so NewFastaParser succeeds,
	// but the index file path will be in a non-existent dir.
	//
	// The index path is constructed as fastaPath+".fai", so we create a file
	// whose path ends with something in a non-existent dir.
	badDir := t.TempDir()
	// Write a valid FASTA at a path whose ".fai" sibling is in a non-existent subdir.
	// Since the index path = fastaPath+".fai" and fastaPath = badDir/noexist/ref.fasta,
	// we can't easily control this without the file existing.
	// Alternative: build a parser pointing to a real file, then
	// verify SaveIndex fails when the index path's parent dir doesn't exist.
	// We do this by copying the real file to a path inside badDir/sub/,
	// but sub/ doesn't exist – so creating the .fai fails.
	nonExistentSubdir := filepath.Join(badDir, "subdir_that_does_not_exist")
	realFastaInSubdir := filepath.Join(nonExistentSubdir, "ref.fasta")
	// nonExistentSubdir doesn't exist, so NewFastaParser will fail trying to
	// open the file. Instead, create a parser for a file whose .fai resides
	// in a directory we make unwritable only after parser creation.
	// Simplest approach: create a dir, create the FASTA file inside, then
	// make the dir read-only so .fai creation fails.
	subDir := filepath.Join(badDir, "subdir")
	require.NoError(t, os.MkdirAll(subDir, 0755))
	realFasta := filepath.Join(subDir, "ref.fasta")
	require.NoError(t, os.WriteFile(realFasta, []byte(smallFasta), 0644))

	fpBad, err := reference.NewFastaParser(realFasta)
	require.NoError(t, err)

	// Make the subdir read-only so the .fai file can't be created
	require.NoError(t, os.Chmod(subDir, 0555))
	t.Cleanup(func() { _ = os.Chmod(subDir, 0755) })

	err = fpBad.SaveIndex(index)
	// On macOS/Linux, creating a file in a read-only dir should fail
	if err == nil {
		t.Skip("filesystem does not enforce read-only directory (possibly running as root)")
	}
	assert.Error(t, err)
	_ = realFastaInSubdir // suppress unused variable warning
}

func TestReferenceManager_GetKmerProfile_WithCustomLoader(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{MaxLoaders: 5})
	require.NoError(t, err)
	require.NoError(t, rm.RegisterReferenceGenome("hg38", p, 0))

	profile, err := rm.GetKmerProfile("chr1", 0, 100, 4, "hg38")
	require.NoError(t, err)
	assert.NotNil(t, profile)
}

func TestReferenceManager_GetSequence_UnknownCustomLoader(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{DefaultReferenceFile: p})
	require.NoError(t, err)

	// "nonexistent" loader falls back to default
	seq, err := rm.GetSequence("chr1", 0, 50, "nonexistent")
	require.NoError(t, err)
	assert.Len(t, seq, 50)
}

// ---------------------------------------------------------------------------
// LargeCache eviction path
// ---------------------------------------------------------------------------

func TestLazyReferenceLoader_LargeRegionExceedsLimit(t *testing.T) {
	p := makeDummyFastaFile(t)
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: p})
	require.NoError(t, err)

	// Requesting more than 100MB should fail
	_, err = loader.GetSequence("chr1", 0, 200*1024*1024)
	assert.Error(t, err)
}

// ---------------------------------------------------------------------------
// Chromosome structure exported fields
// ---------------------------------------------------------------------------

func TestChromosome_ExportedFields(t *testing.T) {
	chr := &reference.Chromosome{
		Name:      "chrTest",
		Sequence:  []byte("ATCG"),
		Length:    4,
		GCContent: 0.5,
		NRegions:  []reference.GenomicInterval{{Start: 0, End: 2}},
	}
	assert.Equal(t, "chrTest", chr.Name)
	assert.Equal(t, 4, chr.Length)
	assert.Equal(t, 0.5, chr.GCContent)
	assert.Len(t, chr.NRegions, 1)
}

func TestGenomicRegion_Fields(t *testing.T) {
	gr := reference.GenomicRegion{Chromosome: "chr1", Start: 1000, End: 2000}
	assert.Equal(t, "chr1", gr.Chromosome)
	assert.Equal(t, int64(1000), gr.Start)
	assert.Equal(t, int64(2000), gr.End)
}

// ---------------------------------------------------------------------------
// Preload with kmerSize=0 (no kmer profile) path
// ---------------------------------------------------------------------------

func TestLazyReferenceLoader_NoPreloadRegions_KmerZero(t *testing.T) {
	p := makeDummyFastaFile(t)
	// Empty PreloadRegions = no preload attempt, just verify creation succeeds
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile:  p,
		PreloadRegions: []string{}, // no preload
		KmerSize:       0,
	})
	require.NoError(t, err)
	assert.NotNil(t, loader)
}

// ---------------------------------------------------------------------------
// strings helper check
// ---------------------------------------------------------------------------

func TestChromosome_NRegionsAtBoundaries(t *testing.T) {
	dir := t.TempDir()
	// Sequence that starts AND ends with Ns
	fasta := ">chrBound\nNNATCGNN\n"
	fastaPath := filepath.Join(dir, "bound.fasta")
	require.NoError(t, os.WriteFile(fastaPath, []byte(fasta), 0644))

	parser := reference.NewReferenceGenomeParser(fastaPath)
	require.NoError(t, parser.CreateIndex())
	require.NoError(t, parser.LoadIndex())

	chr, err := parser.LoadChromosome("chrBound")
	require.NoError(t, err)
	assert.Equal(t, 2, len(chr.NRegions))

	// Verify boundaries
	assert.Equal(t, 0, chr.NRegions[0].Start)
	assert.Equal(t, 2, chr.NRegions[0].End)
	assert.Equal(t, 6, chr.NRegions[1].Start)
	assert.Equal(t, 8, chr.NRegions[1].End)
}

// ---------------------------------------------------------------------------
// LazyReferenceLoader with kmerSize set but no preload regions
// ---------------------------------------------------------------------------

func TestLazyReferenceLoader_NoPreloadRegions_WithKmer(t *testing.T) {
	p := makeDummyFastaFile(t)
	// No PreloadRegions – kmerSize only matters when preloading
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{
		ReferenceFile:  p,
		PreloadRegions: []string{},
		KmerSize:       4,
	})
	require.NoError(t, err)
	assert.NotNil(t, loader)
}

// ensure strings import is used (avoids import error if inlined above)
var _ = strings.Contains

// ---------------------------------------------------------------------------
// Regression: fasta_parser panics on malformed (empty) FASTA header (#125)
// ---------------------------------------------------------------------------

func TestFastaParser_BuildIndex_EmptyHeader(t *testing.T) {
	// A FASTA file where a header line is exactly ">" with no chromosome name
	// previously caused strings.Fields("")[0] → panic (index out of range).
	fastaContent := ">chr1\nACGTACGT\n>\nGGGG\n>chr2\nTTTTTTTT\n"
	f, err := os.CreateTemp(t.TempDir(), "empty_header_*.fa")
	require.NoError(t, err)
	_, err = f.WriteString(fastaContent)
	require.NoError(t, err)
	f.Close()

	parser, err := reference.NewFastaParser(f.Name())
	require.NoError(t, err)
	require.NoError(t, parser.Open())
	defer parser.Close()

	// Must not panic; malformed header is silently skipped.
	index, err := parser.BuildIndex()
	require.NoError(t, err)
	// chr1 and chr2 must be indexed; the empty header must be absent.
	names := parser.GetChromosomeNames(index)
	assert.Contains(t, names, "chr1")
	assert.Contains(t, names, "chr2")
	_, hasEmpty := index.Chromosomes[""]
	assert.False(t, hasEmpty, "empty header should not be indexed")
}

func TestFastaParser_GetSequence_EmptyHeader(t *testing.T) {
	// Same malformed file — getSequenceLinearScan must not panic when scanning
	// past the empty header line.
	fastaContent := ">\nAAAAAAAA\n>chr1\nACGTACGT\n"
	f, err := os.CreateTemp(t.TempDir(), "empty_header2_*.fa")
	require.NoError(t, err)
	_, err = f.WriteString(fastaContent)
	require.NoError(t, err)
	f.Close()

	parser, err := reference.NewFastaParser(f.Name())
	require.NoError(t, err)
	require.NoError(t, parser.Open())
	defer parser.Close()

	// Build index so GetSequence can use it; must not panic during indexing.
	index, err := parser.BuildIndex()
	require.NoError(t, err)
	// chr1 must be reachable even though the file starts with a malformed header.
	seq, err := parser.GetSequence(index, "chr1", 0, 8)
	require.NoError(t, err)
	assert.Equal(t, "ACGTACGT", string(seq))
}
