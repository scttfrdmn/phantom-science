//go:build !release

package genomics

import (
	"bytes"
	"io"
	"math"
	"os"
	"sort"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore/errmodel"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// tierConfig is an ONT run with basecaller tier (#27: tiers apply only to
// platform ont).
func tierConfig(t *testing.T, tier string) types.SequencerConfig {
	c := goldenConfig(t.TempDir())
	c.Platform, c.Illumina = "ont", types.IlluminaSpec{}
	c.ReadLength = 2000
	c.ReadsPerBatch = 100
	c.TotalBatches = 2
	c.Seed = 5
	c.ErrorModel.Tier = tier
	return c
}

// meanDeletionRate is the expected per-base deletion probability of m on
// uniform random sequence (a base lies in a run of length L with
// probability L·(3/4)²·(1/4)^(L−1)).
func meanDeletionRate(m *errmodel.Model) float64 {
	d := 0.0
	for L := 1; L <= 64; L++ {
		d += float64(L) * 9 / 16 * math.Pow(0.25, float64(L-1)) * m.DeletionRate(L)
	}
	return d
}

func phred(qual string) []byte {
	q := []byte(qual)
	for i := range q {
		q[i] -= '!'
	}
	return q
}

// TestErrorModel_TierOnStreamPath runs GenerateToWriter with each tier:
// read lengths shift by the tier's insertion − deletion rate, per-read ONT
// mean qualities peak at the tier mode (ordered sup > hac > fast), the output
// differs from the error-free run with the same seed, and it is reproducible.
func TestErrorModel_TierOnStreamPath(t *testing.T) {
	plain := parseFASTQ(t, string(streamOutput(t, tierConfig(t, ""), GenerateOpts{})))
	for _, r := range plain {
		if len(r.seq) != 2000 || len(r.qual) != 2000 {
			t.Fatalf("error-free read has length %d/%d, want 2000", len(r.seq), len(r.qual))
		}
	}
	var medians []float64
	for _, tier := range errmodel.Tiers {
		t.Run(string(tier), func(t *testing.T) {
			c := tierConfig(t, string(tier))
			out := streamOutput(t, c, GenerateOpts{})
			if !bytes.Equal(out, streamOutput(t, c, GenerateOpts{})) {
				t.Fatal("same seed and tier produced different output")
			}
			recs := parseFASTQ(t, string(out))
			if len(recs) != len(plain) {
				t.Fatalf("%d reads, want %d", len(recs), len(plain))
			}
			m, err := errmodel.New(tier)
			if err != nil {
				t.Fatal(err)
			}
			p := m.Params()
			var sumLen float64
			var meanQs []float64
			changed := 0
			for i, r := range recs {
				if len(r.qual) != len(r.seq) {
					t.Fatalf("read %d: %d qualities for %d bases", i, len(r.qual), len(r.seq))
				}
				if !acgt.MatchString(r.seq) {
					t.Fatalf("read %d has non-ACGT bases", i)
				}
				if r.seq != plain[i].seq {
					changed++
				}
				sumLen += float64(len(r.seq))
				meanQs = append(meanQs, errmodel.MeanQ(phred(r.qual)))
			}
			if changed < len(recs)-1 {
				t.Errorf("only %d of %d reads differ from the error-free run", changed, len(recs))
			}
			// Length = L + insertions − deletions: mean shift L·(Ins − D̄),
			// per-read variance ≈ L·(Ins + D̄).
			L, n := float64(c.ReadLength), float64(len(recs))
			del := meanDeletionRate(m)
			want := L * (1 + p.Ins - del)
			tol := 4 * math.Sqrt(L*(p.Ins+del)/n)
			if got := sumLen / n; math.Abs(got-want) > tol {
				t.Errorf("mean read length %.1f, want %.1f ± %.1f (4σ)", got, want, tol)
			} else {
				t.Logf("mean read length %.1f (want %.1f ± %.1f)", got, want, tol)
			}
			sort.Float64s(meanQs)
			med := meanQs[len(meanQs)/2]
			if math.Abs(med-p.ModalQ) > 1 {
				t.Errorf("median read MeanQ %.2f, want Q%.0f ± 1", med, p.ModalQ)
			}
			medians = append(medians, med)
		})
	}
	if len(medians) == 3 && !(medians[0] > medians[1] && medians[1] > medians[2]) {
		t.Errorf("median MeanQ not ordered sup > hac > fast: %v", medians)
	}
}

// TestErrorModel_TierFileOutput checks the file path (Generate) applies the
// tier too.
func TestErrorModel_TierFileOutput(t *testing.T) {
	c := tierConfig(t, "fast")
	c.ReadLength, c.ReadsPerBatch = 500, 20
	files := fileOutputs(t, c)
	if len(files) == 0 {
		t.Fatal("no files written")
	}
	for name, data := range files {
		varied := false
		for _, r := range parseFASTQ(t, string(data)) {
			if len(r.qual) != len(r.seq) {
				t.Fatalf("%s: %d qualities for %d bases", name, len(r.qual), len(r.seq))
			}
			if len(r.seq) != c.ReadLength {
				varied = true
			}
		}
		if !varied {
			t.Errorf("%s: every read has length %d; no indels applied", name, c.ReadLength)
		}
	}
}

// TestErrorModel_TierNeedsONT: a basecaller tier on Illumina is a
// configuration error (#27; it used to apply the ONT model to Illumina
// reads), on both output paths.
func TestErrorModel_TierNeedsONT(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.ErrorModel.Tier = "hac"
	for _, run := range []func(*Engine) error{
		func(e *Engine) error { return e.GenerateToWriter(io.Discard, GenerateOpts{}) },
		func(e *Engine) error { return e.Generate() },
	} {
		e := NewEngine(c)
		err := run(e)
		e.Close()
		if err == nil || !strings.Contains(err.Error(), "platform ont") {
			t.Errorf("tier on illumina: err = %v, want a platform error", err)
		}
	}
	if entries, _ := os.ReadDir(c.OutputDir); len(entries) != 0 {
		t.Errorf("rejected config wrote %d files", len(entries))
	}
}

// TestErrorModel_NoTierUnchanged pins the linear model: no tier and
// quality_model linear give the linear-decay qualities, no errors and
// fixed-length reads.
func TestErrorModel_NoTierUnchanged(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Seed = 9
	want := decayQuality(c)
	for _, r := range parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{}))) {
		if r.qual != want || len(r.seq) != c.ReadLength {
			t.Fatalf("default output changed: len %d, qual %q", len(r.seq), r.qual)
		}
	}
}
