//go:build !release

package main

import (
	"bytes"
	"strings"
	"testing"
)

func TestGenerateTestSequences_ShapeAndAlphabet(t *testing.T) {
	seqs := generateTestSequences(7, 123, 1)
	if len(seqs) != 7 {
		t.Fatalf("got %d sequences, want 7", len(seqs))
	}
	for i, s := range seqs {
		if len(s) != 123 {
			t.Errorf("sequence %d has length %d, want 123", i, len(s))
		}
		if strings.Trim(string(s), "ACGT") != "" {
			t.Errorf("sequence %d has a non-ACGT base: %q", i, s)
		}
	}
}

func TestGenerateTestSequences_SeedDeterminism(t *testing.T) {
	a := generateTestSequences(3, 500, 42)
	b := generateTestSequences(3, 500, 42)
	c := generateTestSequences(3, 500, 43)
	for i := range a {
		if !bytes.Equal(a[i], b[i]) {
			t.Fatalf("same seed gave different sequence %d", i)
		}
	}
	if bytes.Equal(a[0], c[0]) {
		t.Fatal("different seeds gave the same sequence")
	}
}

// The old generator produced period-4 sequences with only a handful of
// distinct k-mers, so the benchmark measured a degenerate table. Uniform
// random bases should have each base near 25% and, at k=4, all 256 k-mers.
func TestGenerateTestSequences_UniformComposition(t *testing.T) {
	seqs := generateTestSequences(10, 2000, 7)
	counts := map[byte]int{}
	for _, s := range seqs {
		for _, b := range s {
			counts[b]++
		}
	}
	for _, b := range []byte("ACGT") {
		frac := float64(counts[b]) / 20000
		if frac < 0.23 || frac > 0.27 {
			t.Errorf("base %c fraction %.3f, want about 0.25", b, frac)
		}
	}
	table, _, err := countMap(seqs, 4)
	if err != nil {
		t.Fatal(err)
	}
	if len(table.Frequencies) != 256 {
		t.Errorf("distinct 4-mers = %d, want 256", len(table.Frequencies))
	}
}

func TestCounters_AgreeAndCountEveryWindow(t *testing.T) {
	const n, length, k = 20, 300, 6
	seqs := generateTestSequences(n, length, 3)
	m, _, err := countMap(seqs, k)
	if err != nil {
		t.Fatal(err)
	}
	p, _, err := countPacked(seqs, k)
	if err != nil {
		t.Fatal(err)
	}
	want := n * (length - k + 1)
	if m.TotalCount != want || p.TotalCount != want {
		t.Errorf("total k-mers map=%d packed=%d, want %d", m.TotalCount, p.TotalCount, want)
	}
	if err := verifyResults(m, p, k); err != nil {
		t.Errorf("verifyResults on matching tables: %v", err)
	}
}

func TestVerifyResults_DetectsTamperedCount(t *testing.T) {
	seqs := generateTestSequences(5, 200, 9)
	m, _, err := countMap(seqs, 4)
	if err != nil {
		t.Fatal(err)
	}
	p, _, err := countPacked(seqs, 4)
	if err != nil {
		t.Fatal(err)
	}
	// Move one count from one k-mer to another: totals and distinct counts
	// stay the same, so only the per-k-mer comparison can catch it.
	top := m.GetTopKmers(2)
	a, b := top[0].Kmer, top[1].Kmer
	m.Frequencies[a]--
	m.Frequencies[b]++
	err = verifyResults(m, p, 4)
	if err == nil {
		t.Fatal("verifyResults accepted tables with different per-k-mer counts")
	}
	if !strings.Contains(err.Error(), a) && !strings.Contains(err.Error(), b) {
		t.Errorf("error %q does not name the mismatched k-mer (%s or %s)", err, a, b)
	}
}

func TestVerifyResults_DetectsTotalMismatch(t *testing.T) {
	seqs := generateTestSequences(2, 100, 1)
	m, _, _ := countMap(seqs, 3)
	p, _, _ := countPacked(seqs, 3)
	m.TotalCount++
	if err := verifyResults(m, p, 3); err == nil || !strings.Contains(err.Error(), "total count") {
		t.Errorf("verifyResults = %v, want a total-count mismatch", err)
	}
}

func TestMeasure_ReportsRetainedAndAllocatedHeap(t *testing.T) {
	const size = 16 << 20
	buf, m, err := measure(func() ([]byte, error) {
		b := make([]byte, size)
		for i := range b {
			b[i] = byte(i)
		}
		return b, nil
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(buf) != size {
		t.Fatalf("measure did not return fn's value")
	}
	if m.AllocatedBytes < size {
		t.Errorf("allocated %d bytes, want at least %d", m.AllocatedBytes, size)
	}
	// The result is still live at the second GC, so most of it is retained.
	if m.RetainedBytes < size*9/10 {
		t.Errorf("retained %d bytes, want at least %d", m.RetainedBytes, size*9/10)
	}
	if m.Duration <= 0 {
		t.Errorf("duration %v, want > 0", m.Duration)
	}
}

func TestRatio(t *testing.T) {
	cases := []struct {
		a, b float64
		want string
	}{
		{4, 2, "2.00x"},
		{1, 0, "n/a"},
		{1, -5, "n/a"},
		{-1, 5, "n/a"},
	}
	for _, c := range cases {
		if got := ratio(c.a, c.b); got != c.want {
			t.Errorf("ratio(%v, %v) = %q, want %q", c.a, c.b, got, c.want)
		}
	}
}

func TestRun_FullComparison(t *testing.T) {
	var out bytes.Buffer
	cfg := config{k: 5, sequences: 30, length: 200, seed: 1, accuracy: true, memory: true, verbose: true}
	if err := run(cfg, &out); err != nil {
		t.Fatalf("run: %v\n%s", err, out.String())
	}
	s := out.String()
	for _, want := range []string{"Hash map", "Packed array", "results match: true", "Top 5 k-mers", "Retained heap by k"} {
		if !strings.Contains(s, want) {
			t.Errorf("output missing %q:\n%s", want, s)
		}
	}
	// 30 * (200 - 5 + 1) k-mers, reported by both counters.
	if strings.Count(s, "5880 total") != 2 {
		t.Errorf("want both counters to report 5880 k-mers:\n%s", s)
	}
}

func TestRun_SingleCounter(t *testing.T) {
	var out bytes.Buffer
	if err := run(config{k: 4, sequences: 5, length: 50, seed: 1, packedOnly: true}, &out); err != nil {
		t.Fatal(err)
	}
	if strings.Contains(out.String(), "Hash map") || strings.Contains(out.String(), "Comparison") {
		t.Errorf("-packed-only ran the hash-map counter:\n%s", out.String())
	}
}

func TestRun_RejectsBadConfig(t *testing.T) {
	base := config{k: 8, sequences: 10, length: 100, seed: 1}
	cases := map[string]func(c *config){
		"k too large":       func(c *config) { c.k = maxK + 1 },
		"k zero":            func(c *config) { c.k = 0 },
		"no sequences":      func(c *config) { c.sequences = 0 },
		"length below k":    func(c *config) { c.length = 7 },
		"both only flags":   func(c *config) { c.packedOnly, c.mapOnly = true, true },
		"accuracy one side": func(c *config) { c.accuracy, c.mapOnly = true, true },
	}
	for name, mutate := range cases {
		t.Run(name, func(t *testing.T) {
			c := base
			mutate(&c)
			if err := run(c, &bytes.Buffer{}); err == nil {
				t.Error("run accepted an invalid config")
			}
		})
	}
}
