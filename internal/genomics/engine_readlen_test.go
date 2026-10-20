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

// Fixed mode is byte-identical to the output before #10. Since #27 this
// holds for single-end output with quality_model linear (goldenConfig) and
// for ONT: their hashes below were recorded from main at fb54e79, before
// #27, and are unchanged. Paired-end output changed in #27 (each pair is now
// the two ends of one fragment), so its hashes were re-pinned then.
func TestReadLen_FixedModeUnchanged(t *testing.T) {
	sum := func(b []byte) string { h := sha256.Sum256(b); return hex.EncodeToString(h[:]) }
	filesHash := func(files map[string][]byte) string {
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
		return hex.EncodeToString(h.Sum(nil))
	}
	blind := GenerateOpts{SampleID: "s1", Variants: []VariantRecord{{Ref: "AC", Alt: "A", Freq: 0.5}}}
	region := GenerateOpts{RefRegion: &GenomicRegion{Chromosome: "chr1", Start: 1000, End: 5000},
		Variants: []VariantRecord{{Chromosome: "chr1", Position: 3000, Ref: "ACGT", Alt: "A", Freq: 0.5}}}
	for _, tc := range []struct {
		name                          string
		paired                        bool
		wantStream, wantRegion, files string
	}{
		{"single-end linear (pre-#27)", false,
			"80f61e1c8210aad398d38c1298467776d4c293aabd6ab3acbcf51a1652c2622b",
			"219f343a7cd54353cff58649b79e3f6768d6d1d0d3b87cadc811ae780cbe0881",
			"214c36b82f0b8a0a5238ba06373b6f1b8d36a79608e1978e0d23b4961c429096"},
		{"paired linear (re-pinned in #27: fragment model)", true,
			"7846ad0f69f85a4bd223ba04efdaa416864e522bb3eae67194833476d5589f3d",
			"49eca4018346fbd5de77f7f880f4a4e67f4a64a4f3ddb7f5cd92d997acfadc6a",
			"a29ed34b62fc0d77a2aa9bee27aac59910cd3862a523f2bfdd52e19e709c9631"},
	} {
		c := goldenConfig(t.TempDir())
		c.Seed = 3
		c.PairedEnd = tc.paired
		if tc.paired {
			// Pinned with inserts of 350 ± 50, the default before the
			// DRR276626 measurement.
			c.Illumina.InsertSize = types.InsertSizeSpec{Mean: 350, SD: 50}
		}
		for _, dist := range []string{"", types.ReadLengthFixed} {
			c.ReadLengthDistribution.Dist = dist
			if got := sum(streamOutput(t, c, blind)); got != tc.wantStream {
				t.Errorf("%s, dist %q: stream output hash %s, want %s", tc.name, dist, got, tc.wantStream)
			}
		}
		if got := sum(streamOutput(t, c, region)); got != tc.wantRegion {
			t.Errorf("%s: region stream output hash %s, want %s", tc.name, got, tc.wantRegion)
		}
		if got := filesHash(fileOutputs(t, c)); got != tc.files {
			t.Errorf("%s: file output hash %s, want %s", tc.name, got, tc.files)
		}
	}
	o := tierConfig(t, "hac")
	o.ReadsPerBatch = 20
	const wantONT = "d3125bb4bf68bcf2689b0216a9dd5ba3554bab049ca02dc37385c1b7ad8767aa" // pre-#27
	if got := sum(streamOutput(t, o, GenerateOpts{})); got != wantONT {
		t.Errorf("ONT tier stream output hash %s, want %s", got, wantONT)
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
