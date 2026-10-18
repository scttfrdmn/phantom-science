//go:build !release

package validation

import (
	"math/rand/v2"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/readlen"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// fastqWithLengths builds a FASTQ string with one read per length.
func fastqWithLengths(lengths []int) string {
	var b strings.Builder
	for i, l := range lengths {
		b.WriteString("@r")
		b.WriteString(strings.Repeat("x", i%3))
		b.WriteString("\n")
		b.WriteString(strings.Repeat("A", l))
		b.WriteString("\n+\n")
		b.WriteString(strings.Repeat("I", l))
		b.WriteString("\n")
	}
	return b.String()
}

func lognormalSample(t *testing.T, n int, seed uint64) (*readlen.Distribution, []int) {
	t.Helper()
	d, err := readlen.New(types.ReadLengthSpec{Dist: "lognormal", N50: 5000, Median: 2500, Min: 200, Max: 60000}, 0)
	if err != nil {
		t.Fatal(err)
	}
	r := rand.New(rand.NewPCG(seed, 1))
	ls := make([]int, n)
	for i := range ls {
		ls[i] = d.Sample(r)
	}
	return d, ls
}

func distributedExpectation(d *readlen.Distribution) ReadLengthExpectation {
	lo, hi := d.Bounds()
	return ReadLengthExpectation{Mode: ReadLengthsDistributed, Min: lo, Max: hi, Mean: d.Mean(), SD: d.SD(), Name: d.Dist()}
}

func readLengthResult(f *Framework, fastq string) types.ValidationResult {
	return f.runStatisticalTest("read_length_consistency", fastq)
}

func TestReadLengthModes(t *testing.T) {
	cfg := types.ValidationConfig{Enable: true}
	fixed := fastqWithLengths([]int{150, 150, 150, 150})
	d, ls := lognormalSample(t, 2000, 9)
	variable := fastqWithLengths(ls)

	strict := NewFramework(cfg)
	if res := readLengthResult(strict, fixed); !res.Passed {
		t.Errorf("strict, fixed data: %s", res.Message)
	}
	if res := readLengthResult(strict, variable); res.Passed {
		t.Errorf("strict, variable data passed: %s", res.Message)
	}
	exact := NewFramework(cfg, WithReadLengths(ReadLengthExpectation{Length: 150}))
	if res := readLengthResult(exact, fixed); !res.Passed {
		t.Errorf("strict length 150, fixed 150 data: %s", res.Message)
	}
	if res := readLengthResult(NewFramework(cfg, WithReadLengths(ReadLengthExpectation{Length: 100})), fixed); res.Passed {
		t.Errorf("strict length 100 accepted 150-base reads: %s", res.Message)
	}

	anyLen := NewFramework(cfg, WithReadLengths(ReadLengthExpectation{Mode: ReadLengthsAny}))
	for name, data := range map[string]string{"fixed": fixed, "variable": variable} {
		if res := readLengthResult(anyLen, data); !res.Passed {
			t.Errorf("any, %s data: %s", name, res.Message)
		}
	}
	if res := readLengthResult(anyLen, variable); !strings.Contains(res.Message, "N50") {
		t.Errorf("any mode does not report N50: %s", res.Message)
	}

	dist := NewFramework(cfg, WithReadLengths(distributedExpectation(d)))
	res := readLengthResult(dist, variable)
	if !res.Passed {
		t.Errorf("distributed, matching data: %s", res.Message)
	}
	t.Log(res.Message)
}

func TestReadLengthDistributed_Failures(t *testing.T) {
	cfg := types.ValidationConfig{Enable: true}
	d, ls := lognormalSample(t, 2000, 10)
	exp := distributedExpectation(d)
	f := NewFramework(cfg, WithReadLengths(exp))

	for name, mutate := range map[string]func([]int){
		"below-min": func(xs []int) { xs[7] = 199 },
		"above-max": func(xs []int) { xs[3] = 60001 },
		// A 15% scale error: the mean moves by ~0.15·3535 ≈ 530 bases,
		// well beyond 4·SD/√n ≈ 310.
		"mean-shift": func(xs []int) {
			for i := range xs {
				xs[i] = min(max(int(float64(xs[i])*1.15), 200), 60000)
			}
		},
	} {
		xs := append([]int(nil), ls...)
		mutate(xs)
		if res := readLengthResult(f, fastqWithLengths(xs)); res.Passed {
			t.Errorf("%s: passed: %s", name, res.Message)
		}
	}

	// Bounds alone (no expected mean) accept the shifted data.
	exp.Mean, exp.SD = 0, 0
	xs := append([]int(nil), ls...)
	for i := range xs {
		xs[i] = min(int(float64(xs[i])*1.15), 60000)
	}
	if res := readLengthResult(NewFramework(cfg, WithReadLengths(exp)), fastqWithLengths(xs)); !res.Passed {
		t.Errorf("bounds-only check failed: %s", res.Message)
	}
}

// The 4σ mean check must not reject data drawn from the distribution: over
// 200 independent samples a false rejection would have probability ~1e-2.
func TestReadLengthDistributed_NoFalseRejects(t *testing.T) {
	cfg := types.ValidationConfig{Enable: true}
	for seed := range uint64(200) {
		d, ls := lognormalSample(t, 300, 1000+seed)
		f := NewFramework(cfg, WithReadLengths(distributedExpectation(d)))
		if res := readLengthResult(f, fastqWithLengths(ls)); !res.Passed {
			t.Fatalf("seed %d: %s", seed, res.Message)
		}
	}
}

func TestSummarizeLengths(t *testing.T) {
	st := summarizeLengths([]int{4, 2, 6, 3, 5})
	if st.n != 5 || st.min != 2 || st.max != 6 || st.mean != 4 || st.n50 != 5 {
		t.Errorf("summary = %+v, want n 5 min 2 max 6 mean 4 N50 5", st)
	}
	if st := summarizeLengths([]int{100, 1, 1, 1}); st.n50 != 100 {
		t.Errorf("N50 = %d, want 100", st.n50)
	}
}
