//go:build !release

package reference_test

import (
	"compress/gzip"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// readFixture builds two same-length-named chromosomes with different
// random sequence, soft-masked (lower-case) stretches and 7-bp lines.
func readFixture(t *testing.T) (fasta string, seqs map[string]string) {
	t.Helper()
	rng := rand.New(rand.NewPCG(5, 0))
	seqs = map[string]string{}
	var b strings.Builder
	for _, name := range []string{"chr1", "chr2"} {
		s := make([]byte, 503)
		for i := range s {
			s[i] = "ACGT"[rng.IntN(4)]
			if i%50 < 5 {
				s[i] += 'a' - 'A'
			}
		}
		seqs[name] = strings.ToUpper(string(s))
		b.WriteString(">" + name + " test record\n")
		for i := 0; i < len(s); i += 7 {
			b.Write(s[i:min(i+7, len(s))])
			b.WriteByte('\n')
		}
	}
	return b.String(), seqs
}

// TestLazyReferenceLoader_ReadsFASTA: regions must be the FASTA's own bases
// (upper-cased), for plain files with and without a .fai index and for gzip
// files, and chromosomes with equal-length names must not share sequence.
func TestLazyReferenceLoader_ReadsFASTA(t *testing.T) {
	fasta, seqs := readFixture(t)
	dir := t.TempDir()

	plain := filepath.Join(dir, "plain.fa")
	require.NoError(t, os.WriteFile(plain, []byte(fasta), 0o644))

	indexed := filepath.Join(dir, "indexed.fa")
	require.NoError(t, os.WriteFile(indexed, []byte(fasta), 0o644))
	fp, err := reference.NewFastaParser(indexed)
	require.NoError(t, err)
	idx, err := fp.BuildIndex()
	require.NoError(t, err)
	require.NoError(t, fp.SaveIndex(idx))
	require.FileExists(t, indexed+".fai")

	gz := filepath.Join(dir, "ref.fa.gz")
	f, err := os.Create(gz)
	require.NoError(t, err)
	zw := gzip.NewWriter(f)
	_, err = zw.Write([]byte(fasta))
	require.NoError(t, err)
	require.NoError(t, zw.Close())
	require.NoError(t, f.Close())

	regions := [][2]int64{{0, 1}, {0, 7}, {3, 18}, {6, 8}, {100, 401}, {496, 503}, {0, 503}}
	for _, path := range []string{plain, indexed, gz} {
		loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: path})
		require.NoError(t, err)
		for _, chr := range []string{"chr1", "chr2"} {
			for _, r := range regions {
				got, err := loader.GetSequence(chr, r[0], r[1])
				require.NoError(t, err, "%s %s:%d-%d", filepath.Base(path), chr, r[0], r[1])
				assert.Equal(t, seqs[chr][r[0]:r[1]], string(got), "%s %s:%d-%d", filepath.Base(path), chr, r[0], r[1])
			}
		}
		a, _ := loader.GetSequence("chr1", 0, 100)
		b, _ := loader.GetSequence("chr2", 0, 100)
		assert.NotEqual(t, string(a), string(b), "chr1 and chr2 must not share sequence")
	}
}

func TestLazyReferenceLoader_ReadsFASTA_ClipsAndRejects(t *testing.T) {
	fasta, seqs := readFixture(t)
	p := filepath.Join(t.TempDir(), "ref.fa")
	require.NoError(t, os.WriteFile(p, []byte(fasta), 0o644))
	loader, err := reference.NewLazyReferenceLoader(reference.LazyLoaderConfig{ReferenceFile: p})
	require.NoError(t, err)

	// Clipped to the chromosome, like samtools faidx.
	got, err := loader.GetSequence("chr1", -20, 10)
	require.NoError(t, err)
	assert.Equal(t, seqs["chr1"][:10], string(got))
	got, err = loader.GetSequence("chr1", 490, 900)
	require.NoError(t, err)
	assert.Equal(t, seqs["chr1"][490:], string(got))

	_, err = loader.GetSequence("chr1", 600, 700)
	assert.Error(t, err, "region entirely past the end")
	_, err = loader.GetSequence("chr3", 0, 10)
	assert.Error(t, err, "chromosome not in the reference")
}
