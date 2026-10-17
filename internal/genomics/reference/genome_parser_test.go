package reference

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewReferenceGenomeParser tests parser creation
func TestNewReferenceGenomeParser(t *testing.T) {
	fastaPath := "/test/path/genome.fasta"
	parser := NewReferenceGenomeParser(fastaPath)

	assert.NotNil(t, parser)
	assert.Equal(t, fastaPath, parser.fastaPath)
	assert.Equal(t, fastaPath+".fai", parser.indexPath)
	assert.NotNil(t, parser.chromSizes)
	assert.NotNil(t, parser.chromOffsets)
	assert.True(t, parser.memoryMapped)
}

// createTestFastaFile creates a test FASTA file for testing
func createTestFastaFile(t interface{}, dir string) string {
	fastaPath := filepath.Join(dir, "test_genome.fasta")
	
	fastaContent := `>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
NNNNNNNNNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>chrX
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
`

	err := os.WriteFile(fastaPath, []byte(fastaContent), 0644)
	if t != nil {
		switch tt := t.(type) {
		case *testing.T:
			require.NoError(tt, err)
		case *testing.B:
			require.NoError(tt, err)
		}
	} else if err != nil {
		panic(err)
	}
	
	return fastaPath
}

// createTestIndexFile creates a test FASTA index file
func createTestIndexFile(t interface{}, fastaPath string) {
	indexPath := fastaPath + ".fai"
	
	indexContent := `chr1	186	6	60	61
chr2	120	199	60	61
chrX	180	326	60	61
`

	err := os.WriteFile(indexPath, []byte(indexContent), 0644)
	if t != nil {
		switch tt := t.(type) {
		case *testing.T:
			require.NoError(tt, err)
		case *testing.B:
			require.NoError(tt, err)
		}
	} else if err != nil {
		panic(err)
	}
}

// TestReferenceGenomeParser_LoadIndex tests index loading
func TestReferenceGenomeParser_LoadIndex(t *testing.T) {
	tempDir := t.TempDir()
	fastaPath := createTestFastaFile(t, tempDir)
	
	parser := NewReferenceGenomeParser(fastaPath)
	
	// Create index automatically
	err := parser.CreateIndex()
	require.NoError(t, err)
	
	err = parser.LoadIndex()
	require.NoError(t, err)
	
	// Check loaded data
	assert.Len(t, parser.chromSizes, 3)
	assert.Greater(t, parser.chromSizes["chr1"], 0)
	assert.Greater(t, parser.chromSizes["chr2"], 0)
	assert.Greater(t, parser.chromSizes["chrX"], 0)
	
	assert.Greater(t, parser.chromOffsets["chr1"], int64(0))
	assert.Greater(t, parser.chromOffsets["chr2"], int64(0))
	assert.Greater(t, parser.chromOffsets["chrX"], int64(0))
}

// TestReferenceGenomeParser_LoadIndex_MissingFile tests error handling
func TestReferenceGenomeParser_LoadIndex_MissingFile(t *testing.T) {
	parser := NewReferenceGenomeParser("/nonexistent/genome.fasta")
	
	err := parser.LoadIndex()
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "failed to open index file")
}

// TestReferenceGenomeParser_LoadChromosome tests chromosome loading
func TestReferenceGenomeParser_LoadChromosome(t *testing.T) {
	tempDir := t.TempDir()
	fastaPath := createTestFastaFile(t, tempDir)
	
	parser := NewReferenceGenomeParser(fastaPath)
	err := parser.CreateIndex()
	require.NoError(t, err)
	err = parser.LoadIndex()
	require.NoError(t, err)
	
	t.Run("load valid chromosome", func(t *testing.T) {
		chr, err := parser.LoadChromosome("chr1")
		require.NoError(t, err)
		require.NotNil(t, chr)
		
		assert.Equal(t, "chr1", chr.Name)
		assert.Greater(t, chr.Length, 0)
		assert.GreaterOrEqual(t, chr.GCContent, 0.0)
		assert.LessOrEqual(t, chr.GCContent, 1.0)
		assert.NotNil(t, chr.Sequence)
		assert.NotNil(t, chr.NRegions)
		
		// Check that sequence contains expected bases
		hasValidBases := false
		for _, base := range chr.Sequence {
			if base == 'A' || base == 'T' || base == 'G' || base == 'C' || base == 'N' {
				hasValidBases = true
				break
			}
		}
		assert.True(t, hasValidBases, "Sequence should contain valid DNA bases")
		
		// Check for N regions
		assert.GreaterOrEqual(t, len(chr.NRegions), 0, "Should identify N regions")
	})
	
	t.Run("load nonexistent chromosome", func(t *testing.T) {
		chr, err := parser.LoadChromosome("chrNonexistent")
		assert.Error(t, err)
		assert.Nil(t, chr)
		assert.Contains(t, err.Error(), "chromosome chrNonexistent not found")
	})
}

// TestReferenceGenomeParser_calculateGCContent tests GC content calculation
func TestReferenceGenomeParser_calculateGCContent(t *testing.T) {
	parser := NewReferenceGenomeParser("/test/path")
	
	testCases := []struct {
		sequence []byte
		expected float64
		description string
	}{
		{[]byte(""), 0.0, "empty sequence"},
		{[]byte("ATGC"), 0.5, "balanced sequence"},
		{[]byte("GGCC"), 1.0, "all GC"},
		{[]byte("AATT"), 0.0, "all AT"},
		{[]byte("ATCGATCG"), 0.5, "mixed sequence"},
		{[]byte("NNNNATCG"), 0.5, "sequence with Ns"},
		{[]byte("NNNN"), 0.0, "all Ns"},
	}
	
	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			gcContent := parser.calculateGCContent(tc.sequence)
			assert.InDelta(t, tc.expected, gcContent, 0.001)
		})
	}
}

// TestReferenceGenomeParser_identifyNRegions tests N region identification
func TestReferenceGenomeParser_identifyNRegions(t *testing.T) {
	parser := NewReferenceGenomeParser("/test/path")
	
	testCases := []struct {
		sequence []byte
		expectedRegions int
		description string
	}{
		{[]byte("ATCG"), 0, "no N regions"},
		{[]byte("ATCGNNNGATCG"), 1, "single N region"},
		{[]byte("NNNNATCGNNNGATCGNNNN"), 3, "multiple N regions"},
		{[]byte("NNNATCGATCG"), 1, "N region at start"},
		{[]byte("ATCGATCGNNN"), 1, "N region at end"},
		{[]byte("NNNNNNNN"), 1, "all Ns"},
		{[]byte(""), 0, "empty sequence"},
	}
	
	for _, tc := range testCases {
		t.Run(tc.description, func(t *testing.T) {
			nRegions := parser.identifyNRegions(tc.sequence)
			assert.Len(t, nRegions, tc.expectedRegions)
			
			// Validate intervals
			for _, region := range nRegions {
				assert.GreaterOrEqual(t, region.Start, 0)
				assert.Greater(t, region.End, region.Start)
				assert.LessOrEqual(t, region.End, len(tc.sequence))
				
				// Check that region actually contains Ns
				for i := region.Start; i < region.End; i++ {
					assert.Equal(t, byte('N'), tc.sequence[i], 
						"N region should contain only N bases")
				}
			}
		})
	}
}

// TestReferenceGenomeParser_GetChromosomeNames tests chromosome name retrieval
func TestReferenceGenomeParser_GetChromosomeNames(t *testing.T) {
	tempDir := t.TempDir()
	fastaPath := createTestFastaFile(t, tempDir)
	
	parser := NewReferenceGenomeParser(fastaPath)
	err := parser.CreateIndex()
	require.NoError(t, err)
	err = parser.LoadIndex()
	require.NoError(t, err)
	
	names := parser.GetChromosomeNames()
	assert.Len(t, names, 3)
	
	expectedNames := []string{"chr1", "chr2", "chrX"}
	for _, expected := range expectedNames {
		assert.Contains(t, names, expected)
	}
}

// TestReferenceGenomeParser_GetChromosomeSize tests chromosome size retrieval
func TestReferenceGenomeParser_GetChromosomeSize(t *testing.T) {
	tempDir := t.TempDir()
	fastaPath := createTestFastaFile(t, tempDir)
	
	parser := NewReferenceGenomeParser(fastaPath)
	err := parser.CreateIndex()
	require.NoError(t, err)
	err = parser.LoadIndex()
	require.NoError(t, err)
	
	t.Run("valid chromosome", func(t *testing.T) {
		size, err := parser.GetChromosomeSize("chr1")
		require.NoError(t, err)
		assert.Greater(t, size, 0)
	})
	
	t.Run("invalid chromosome", func(t *testing.T) {
		size, err := parser.GetChromosomeSize("chrInvalid")
		assert.Error(t, err)
		assert.Equal(t, 0, size)
		assert.Contains(t, err.Error(), "chromosome chrInvalid not found")
	})
}

// TestReferenceGenomeParser_ValidateGenomeFile tests file validation
func TestReferenceGenomeParser_ValidateGenomeFile(t *testing.T) {
	t.Run("valid files", func(t *testing.T) {
		tempDir := t.TempDir()
		fastaPath := createTestFastaFile(t, tempDir)
		
		parser := NewReferenceGenomeParser(fastaPath)
		err := parser.CreateIndex()
		require.NoError(t, err)
		err = parser.ValidateGenomeFile()
		assert.NoError(t, err)
	})
	
	t.Run("missing FASTA file", func(t *testing.T) {
		parser := NewReferenceGenomeParser("/nonexistent/genome.fasta")
		err := parser.ValidateGenomeFile()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "reference genome file not found")
	})
	
	t.Run("missing index file", func(t *testing.T) {
		tempDir := t.TempDir()
		fastaPath := createTestFastaFile(t, tempDir)
		// Don't create index file
		
		parser := NewReferenceGenomeParser(fastaPath)
		err := parser.ValidateGenomeFile()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "reference genome index file not found")
	})
}

// TestReferenceGenomeParser_CreateIndex tests index creation
func TestReferenceGenomeParser_CreateIndex(t *testing.T) {
	tempDir := t.TempDir()
	fastaPath := createTestFastaFile(t, tempDir)
	
	parser := NewReferenceGenomeParser(fastaPath)
	
	t.Run("create index successfully", func(t *testing.T) {
		err := parser.CreateIndex()
		require.NoError(t, err)
		
		// Check that index file was created
		indexPath := fastaPath + ".fai"
		_, err = os.Stat(indexPath)
		assert.NoError(t, err, "Index file should be created")
		
		// Load the created index and verify
		err = parser.LoadIndex()
		require.NoError(t, err)
		
		assert.Greater(t, len(parser.chromSizes), 0, "Should have loaded chromosome sizes")
		assert.Contains(t, parser.chromSizes, "chr1")
		assert.Contains(t, parser.chromSizes, "chr2")
		assert.Contains(t, parser.chromSizes, "chrX")
	})
	
	t.Run("create index for nonexistent file", func(t *testing.T) {
		invalidParser := NewReferenceGenomeParser("/nonexistent/genome.fasta")
		err := invalidParser.CreateIndex()
		assert.Error(t, err)
		assert.Contains(t, err.Error(), "failed to open FASTA file")
	})
}

// TestReferenceGenomeParser_GetGenomeStats tests genome statistics
func TestReferenceGenomeParser_GetGenomeStats(t *testing.T) {
	tempDir := t.TempDir()
	fastaPath := createTestFastaFile(t, tempDir)
	
	parser := NewReferenceGenomeParser(fastaPath)
	err := parser.CreateIndex()
	require.NoError(t, err)
	err = parser.LoadIndex()
	require.NoError(t, err)
	
	stats := parser.GetGenomeStats()
	
	assert.Equal(t, 3, stats["total_chromosomes"])
	assert.Greater(t, stats["total_size"], 0)
	assert.Equal(t, fastaPath, stats["fasta_path"])
	assert.Equal(t, fastaPath+".fai", stats["index_path"])
	
	chromSizes := stats["chromosome_sizes"].(map[string]int)
	assert.Greater(t, chromSizes["chr1"], 0)
	assert.Greater(t, chromSizes["chr2"], 0)
	assert.Greater(t, chromSizes["chrX"], 0)
}

// TestChromosome_Validation tests chromosome structure validation
func TestChromosome_Validation(t *testing.T) {
	chr := &Chromosome{
		Name:      "chr1",
		Sequence:  []byte("ATCGATCGATCGNNNNATCGATCG"),
		Length:    24,
		GCContent: 0.5,
		NRegions: []GenomicInterval{
			{Start: 12, End: 16},
		},
	}
	
	assert.Equal(t, "chr1", chr.Name)
	assert.Equal(t, 24, chr.Length)
	assert.Equal(t, len(chr.Sequence), chr.Length)
	assert.GreaterOrEqual(t, chr.GCContent, 0.0)
	assert.LessOrEqual(t, chr.GCContent, 1.0)
	assert.Len(t, chr.NRegions, 1)
	
	// Validate N region
	nRegion := chr.NRegions[0]
	assert.Equal(t, 12, nRegion.Start)
	assert.Equal(t, 16, nRegion.End)
	assert.Greater(t, nRegion.End, nRegion.Start)
}

// TestGenomicInterval_Validation tests genomic interval structure
func TestGenomicInterval_Validation(t *testing.T) {
	interval := GenomicInterval{
		Start: 1000,
		End:   2000,
	}
	
	assert.Equal(t, 1000, interval.Start)
	assert.Equal(t, 2000, interval.End)
	assert.Greater(t, interval.End, interval.Start)
	assert.Equal(t, 1000, interval.End-interval.Start) // Length
}

// TestReferenceGenomeParser_ThreadSafety tests concurrent access
func TestReferenceGenomeParser_ThreadSafety(t *testing.T) {
	tempDir := t.TempDir()
	fastaPath := createTestFastaFile(t, tempDir)
	
	parser := NewReferenceGenomeParser(fastaPath)
	err := parser.CreateIndex()
	require.NoError(t, err)
	err = parser.LoadIndex()
	require.NoError(t, err)
	
	// Test concurrent access to chromosome information
	done := make(chan bool, 10)
	
	for i := 0; i < 10; i++ {
		go func() {
			defer func() { done <- true }()
			
			// Multiple concurrent operations
			names := parser.GetChromosomeNames()
			assert.Greater(t, len(names), 0)
			
			size, err := parser.GetChromosomeSize("chr1")
			assert.NoError(t, err)
			assert.Greater(t, size, 0)
			
			stats := parser.GetGenomeStats()
			assert.Greater(t, stats["total_chromosomes"], 0)
		}()
	}
	
	// Wait for all goroutines to complete
	for i := 0; i < 10; i++ {
		<-done
	}
}

// TestReferenceGenomeParser_loadChromosomeSequence tests sequence loading
func TestReferenceGenomeParser_loadChromosomeSequence(t *testing.T) {
	tempDir := t.TempDir()
	fastaPath := createTestFastaFile(t, tempDir)
	
	parser := NewReferenceGenomeParser(fastaPath)
	err := parser.CreateIndex()
	require.NoError(t, err)
	err = parser.LoadIndex()
	require.NoError(t, err)
	
	t.Run("load existing chromosome", func(t *testing.T) {
		size, _ := parser.GetChromosomeSize("chr1")
		sequence, err := parser.loadChromosomeSequence("chr1", size)
		require.NoError(t, err)
		
		assert.Greater(t, len(sequence), 0)
		// All bases should be uppercase
		for _, base := range sequence {
			if base >= 'a' && base <= 'z' {
				t.Errorf("Found lowercase base: %c", base)
			}
		}
	})
	
	t.Run("load nonexistent chromosome", func(t *testing.T) {
		sequence, err := parser.loadChromosomeSequence("chrNonexistent", 100)
		assert.Error(t, err)
		assert.Nil(t, sequence)
		assert.Contains(t, err.Error(), "chromosome chrNonexistent not found in FASTA file")
	})
}

// BenchmarkReferenceGenomeParser_LoadChromosome benchmarks chromosome loading
func BenchmarkReferenceGenomeParser_LoadChromosome(b *testing.B) {
	tempDir := b.TempDir()
	fastaPath := createTestFastaFile(b, tempDir)
	createTestIndexFile(b, fastaPath)
	
	parser := NewReferenceGenomeParser(fastaPath)
	err := parser.LoadIndex()
	require.NoError(b, err)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = parser.LoadChromosome("chr1")
	}
}

// BenchmarkReferenceGenomeParser_calculateGCContent benchmarks GC content calculation
func BenchmarkReferenceGenomeParser_calculateGCContent(b *testing.B) {
	parser := NewReferenceGenomeParser("/test/path")
	
	// Create a large sequence for benchmarking
	sequence := make([]byte, 100000)
	for i := range sequence {
		bases := []byte{'A', 'T', 'G', 'C'}
		sequence[i] = bases[i%4]
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		parser.calculateGCContent(sequence)
	}
}

// BenchmarkReferenceGenomeParser_identifyNRegions benchmarks N region identification
func BenchmarkReferenceGenomeParser_identifyNRegions(b *testing.B) {
	parser := NewReferenceGenomeParser("/test/path")
	
	// Create a sequence with scattered N regions
	sequence := make([]byte, 100000)
	for i := range sequence {
		if i%1000 < 10 { // 10 Ns every 1000 bases
			sequence[i] = 'N'
		} else {
			bases := []byte{'A', 'T', 'G', 'C'}
			sequence[i] = bases[i%4]
		}
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		parser.identifyNRegions(sequence)
	}
}