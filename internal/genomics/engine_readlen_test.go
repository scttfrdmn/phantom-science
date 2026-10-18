//go:build !release

package genomics

// Engine-level tests for read-length distributions (#10): sampled lengths
// follow the configured distribution, output is reproducible from the seed,
// fixed mode is unchanged, and indels change variable-length reads by
// exactly the indel size.

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"math"
	"slices"
	"sort"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/readlen"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func lognormalConfig(dir string) types.SequencerConfig {
	c := goldenConfig(dir)
	c.Seed = 1
	c.ReadsPerBatch = 1500
	c.TotalBatches = 2
	c.ReadLengthDistribution = types.ReadLengthSpec{Dist: types.ReadLengthLognormal, N50: 5000, Median: 2500}
	return c
}

func readLengths(recs []fastqRec) []int {
	ls := make([]int, len(recs))
	for i, r := range recs {
		ls[i] = len(r.seq)
	}
	return ls
}

// ksStat is the KS distance of integer samples against a discrete CDF.
func ksStat(xs []int, cdf func(int) float64) float64 {
	xs = slices.Clone(xs)
	sort.Ints(xs)
	n := float64(len(xs))
	var d float64
	for i := 0; i < len(xs); {
		j := i
		for j < len(xs) && xs[j] == xs[i] {
			j++
		}
		d = math.Max(d, math.Abs(float64(j)/n-cdf(xs[i])))
		d = math.Max(d, math.Abs(float64(i)/n-cdf(xs[i]-1)))
		i = j
	}
	return d
}

func TestReadLen_EngineLognormal(t *testing.T) {
	c := lognormalConfig(t.TempDir())
	out := streamOutput(t, c, GenerateOpts{})
	recs := parseFASTQ(t, string(out))
	if len(recs) != c.ReadsPerBatch*c.TotalBatches {
		t.Fatalf("%d reads, want %d", len(recs), c.ReadsPerBatch*c.TotalBatches)
	}
	for _, r := range recs {
		if len(r.qual) != len(r.seq) || !acgt.MatchString(r.seq) {
			t.Fatalf("bad record: %d bases, %d qualities", len(r.seq), len(r.qual))
		}
	}
	d, err := readlen.New(c.ReadLengthDistribution, c.ReadLength)
	if err != nil {
		t.Fatal(err)
	}
	ls := readLengths(recs)
	n := float64(len(ls))
	if ks, crit := ksStat(ls, d.CDF), 1.63/math.Sqrt(n); ks >= crit {
		t.Errorf("KS D = %.4f >= %.4f: lengths do not follow lognormal(N50 5000, median 2500)", ks, crit)
	}
	var sum float64
	for _, l := range ls {
		sum += float64(l)
	}
	if mean, tol := sum/n, 4*d.SD()/math.Sqrt(n); math.Abs(mean-d.Mean()) > tol {
		t.Errorf("mean length %.0f, want %.0f ± %.0f", mean, d.Mean(), tol)
	}

	// Same seed: byte-identical; different seed: different lengths.
	if again := streamOutput(t, c, GenerateOpts{}); !bytes.Equal(out, again) {
		t.Error("same seed produced different output")
	}
	c.Seed = 2
	if other := readLengths(parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{})))); slices.Equal(other, ls) {
		t.Error("different seeds produced identical lengths")
	}
}

// Paired-end reads: in an interleaved stream both mates of a pair share the
// sampled length; file output (one RNG stream per mate file) is
// reproducible.
func TestReadLen_EnginePaired(t *testing.T) {
	c := lognormalConfig("")
	c.PairedEnd = true
	c.ReadsPerBatch = 50
	recs := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{})))
	ls := readLengths(recs)
	for i := 0; i+1 < len(ls); i += 2 {
		if ls[i] != ls[i+1] {
			t.Fatalf("pair %d: mate lengths %d and %d differ", i/2, ls[i], ls[i+1])
		}
	}
	if slices.Min(ls) == slices.Max(ls) {
		t.Errorf("all reads have length %d", ls[0])
	}
	a, b := fileOutputs(t, c), fileOutputs(t, c)
	if len(a) != 2*c.TotalBatches {
		t.Fatalf("%d files, want %d", len(a), 2*c.TotalBatches)
	}
	for name, data := range a {
		if !bytes.Equal(data, b[name]) {
			t.Errorf("%s differs between runs with the same seed", name)
		}
		if ls := readLengths(parseFASTQ(t, string(data))); slices.Min(ls) == slices.Max(ls) {
			t.Errorf("%s: all reads have length %d", name, ls[0])
		}
	}
}

// Fixed mode is byte-identical to the output before #10 (hashes recorded
// from main at 88fe7e5).
func TestReadLen_FixedModeUnchanged(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Seed = 3
	c.PairedEnd = true
	sum := func(b []byte) string { h := sha256.Sum256(b); return hex.EncodeToString(h[:]) }
	const wantStream = "ca0fde4f46c349a7bb79d03af3eca2862950683d7020238ac0f41ac05532fd62"
	for _, dist := range []string{"", types.ReadLengthFixed} {
		c.ReadLengthDistribution.Dist = dist
		got := sum(streamOutput(t, c, GenerateOpts{SampleID: "s1", Variants: []VariantRecord{{Ref: "AC", Alt: "A", Freq: 0.5}}}))
		if got != wantStream {
			t.Errorf("dist %q: stream output hash %s, want %s", dist, got, wantStream)
		}
	}
	region := &GenomicRegion{Chromosome: "chr1", Start: 1000, End: 5000}
	const wantRegion = "6932a4f95a3a73241ae9b4918000c00af6242788b89ea3228bacc3690a5e74a8"
	got := sum(streamOutput(t, c, GenerateOpts{RefRegion: region, Variants: []VariantRecord{{Chromosome: "chr1", Position: 3000, Ref: "ACGT", Alt: "A", Freq: 0.5}}}))
	if got != wantRegion {
		t.Errorf("region stream output hash %s, want %s", got, wantRegion)
	}
	files := fileOutputs(t, c)
	names := make([]string, 0, len(files))
	for n := range files {
		names = append(names, n)
	}
	sort.Strings(names)
	h := sha256.New()
	for _, n := range names {
		h.Write([]byte(n))
		h.Write(files[n])
	}
	const wantFiles = "cfe66cc418a7100e234f03a1a0c2c0a92d54d4b694669fc66b4d721dbebb3748"
	if got := hex.EncodeToString(h.Sum(nil)); got != wantFiles {
		t.Errorf("file output hash %s, want %s", got, wantFiles)
	}
}

// identity returns variants with ALT replaced by REF: they consume the same
// randomness as variants but leave reads unchanged.
func identity(vs []VariantRecord) []VariantRecord {
	out := slices.Clone(vs)
	for i := range out {
		out[i].Alt = out[i].Ref
	}
	return out
}

// In variable-length mode an indel changes the read length by exactly the
// indel size, and nothing is padded or trimmed: each read equals the
// baseline read (same seed, identity alleles) with the allele spliced in.
func TestReadLen_IndelChangesLengthExactly(t *testing.T) {
	cases := []struct {
		name     string
		variants []VariantRecord
		region   *GenomicRegion
	}{
		{"blind-deletion", []VariantRecord{{Ref: "ACGTA", Alt: "A", Freq: 1}}, nil},
		{"blind-insertion", []VariantRecord{{Ref: "GATC", Alt: "GATTTTTTTC", Freq: 1}}, nil},
		{"positioned-deletion", []VariantRecord{{Chromosome: "chr1", Position: 10001, Ref: "CAGTACGTAC", Alt: "C", Freq: 1}},
			&GenomicRegion{Chromosome: "chr1", Start: 0, End: 20000}},
		{"positioned-insertion", []VariantRecord{{Chromosome: "chr1", Position: 10001, Ref: "T", Alt: "TGGGGGG", Freq: 1}},
			&GenomicRegion{Chromosome: "chr1", Start: 0, End: 20000}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			c := lognormalConfig(t.TempDir())
			c.ReadLengthDistribution = types.ReadLengthSpec{Dist: types.ReadLengthLognormal, N50: 2000, Median: 1000, Max: 15000}
			v := tc.variants[0]
			delta := len(v.Alt) - len(v.Ref)
			base := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{Variants: identity(tc.variants), RefRegion: tc.region})))
			got := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{Variants: tc.variants, RefRegion: tc.region})))
			if len(base) != len(got) {
				t.Fatalf("%d reads vs %d baseline reads", len(got), len(base))
			}
			changed := 0
			for i := range got {
				b, g := base[i], got[i]
				if strings.ContainsRune(g.seq, 'N') || len(g.qual) != len(g.seq) {
					t.Fatalf("read %d: N padding or misaligned qualities", i)
				}
				var off int
				if tc.region == nil {
					off = strings.Index(b.seq, v.Ref)
				} else {
					start := headerReadStart(t, g.header)
					off = int(v.Position - 1 - start)
					if off < 0 || off+len(v.Ref) > len(b.seq) {
						off = -1
					}
				}
				if off < 0 {
					if g.seq != b.seq {
						t.Fatalf("read %d without the variant changed", i)
					}
					continue
				}
				if len(g.seq)-len(b.seq) != delta {
					t.Fatalf("read %d: length %d → %d, want change %d", i, len(b.seq), len(g.seq), delta)
				}
				if want := b.seq[:off] + v.Alt + b.seq[off+len(v.Ref):]; g.seq != want {
					t.Fatalf("read %d: allele not spliced in place", i)
				}
				changed++
			}
			if changed < 20 {
				t.Errorf("only %d of %d reads carried the variant; test is too weak", changed, len(got))
			}
			t.Logf("%d of %d reads changed by %d bases", changed, len(got), delta)
		})
	}
}
