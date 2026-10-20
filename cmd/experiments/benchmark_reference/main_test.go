//go:build !release

package main

import (
	"bufio"
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// parseFASTA returns record names and sequences in order, and fails on a
// sequence line longer than maxLine.
func parseFASTA(t *testing.T, data []byte, maxLine int) ([]string, []string) {
	t.Helper()
	var names, seqs []string
	sc := bufio.NewScanner(bytes.NewReader(data))
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, ">") {
			names = append(names, line[1:])
			seqs = append(seqs, "")
			continue
		}
		if len(line) > maxLine {
			t.Fatalf("sequence line of %d bases, want at most %d", len(line), maxLine)
		}
		seqs[len(seqs)-1] += line
	}
	return names, seqs
}

func TestWriteTestReference_Format(t *testing.T) {
	var buf bytes.Buffer
	if err := writeTestReference(&buf, 3, 205, 1); err != nil {
		t.Fatal(err)
	}
	names, seqs := parseFASTA(t, buf.Bytes(), 80)
	if strings.Join(names, ",") != "chr1,chr2,chr3" {
		t.Fatalf("records %v, want chr1,chr2,chr3", names)
	}
	for i, s := range seqs {
		if len(s) != 205 {
			t.Errorf("%s has %d bases, want 205", names[i], len(s))
		}
		if strings.Trim(s, "ACGT") != "" {
			t.Errorf("%s has non-ACGT bases", names[i])
		}
	}
	if seqs[0] == seqs[1] {
		t.Error("chromosomes are identical; want independent random sequence")
	}
}

func TestWriteTestReference_SeedDeterminism(t *testing.T) {
	var a, b, c bytes.Buffer
	_ = writeTestReference(&a, 2, 1000, 5)
	_ = writeTestReference(&b, 2, 1000, 5)
	_ = writeTestReference(&c, 2, 1000, 6)
	if !bytes.Equal(a.Bytes(), b.Bytes()) {
		t.Error("same seed gave different references")
	}
	if bytes.Equal(a.Bytes(), c.Bytes()) {
		t.Error("different seeds gave the same reference")
	}
}

func writeFile(t *testing.T, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "ref.fasta")
	if err := os.WriteFile(p, []byte(content), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestReadLayout_UsesTheReferencesOwnNames(t *testing.T) {
	p := writeFile(t, ">contigA desc\nACGTACGTAC\nACG\n>contigB\nTTTT\n")
	chroms, err := readLayout(p)
	if err != nil {
		t.Fatal(err)
	}
	want := []chromInfo{{"contigA", 13}, {"contigB", 4}}
	if len(chroms) != len(want) {
		t.Fatalf("readLayout = %v, want %v", chroms, want)
	}
	for i := range want {
		if chroms[i] != want[i] {
			t.Errorf("chrom %d = %v, want %v", i, chroms[i], want[i])
		}
	}
}

func TestReadLayout_Errors(t *testing.T) {
	if _, err := readLayout(filepath.Join(t.TempDir(), "missing.fasta")); err == nil {
		t.Error("readLayout accepted a missing file")
	}
	if _, err := readLayout(writeFile(t, "")); err == nil {
		t.Error("readLayout accepted a reference with no sequences")
	}
}

func TestPlanRegions_StayInsideChromosomes(t *testing.T) {
	chroms := []chromInfo{{"a", 1000}, {"b", 350}, {"c", 50}}
	regions := planRegions(chroms, 30, 100)
	length := map[string]int64{"a": 1000, "b": 350, "c": 50}
	for i, r := range regions {
		if want := chroms[i%3].Name; r.Chromosome != want {
			t.Errorf("region %d on %s, want %s", i, r.Chromosome, want)
		}
		if r.Start < 0 || r.End > length[r.Chromosome] || r.Start >= r.End {
			t.Errorf("region %d = %s:%d-%d outside chromosome of length %d", i, r.Chromosome, r.Start, r.End, length[r.Chromosome])
		}
	}
	// Chromosome a: 10 adjacent tiles of 100 bases.
	if r := regions[3]; r.Start != 100 || r.End != 200 {
		t.Errorf("second tile on a = %d-%d, want 100-200", r.Start, r.End)
	}
	// Chromosome b has 3 whole tiles; the fourth wraps to the start.
	if r := regions[1+3*3]; r.Start != 0 || r.End != 100 {
		t.Errorf("fourth tile on b = %d-%d, want 0-100 (wrapped)", r.Start, r.End)
	}
	// Chromosome c is shorter than a region: the region is the whole of it.
	if r := regions[2]; r.Start != 0 || r.End != 50 {
		t.Errorf("region on c = %d-%d, want 0-50", r.Start, r.End)
	}
}

func newManager(t *testing.T, path string, readAhead bool) *reference.ReferenceManager {
	t.Helper()
	m, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultReferenceFile: path,
		DefaultCacheSize:     64 << 20,
		EnablePrefetching:    readAhead,
		EnablePackedArrays:   true,
	})
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { m.Close() })
	return m
}

func generatedRegions(t *testing.T, chroms int, length int64, n int, size int64) (string, []reference.GenomicRegion) {
	t.Helper()
	p, err := createTestReference(chroms, length, 1)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { os.Remove(p); os.Remove(p + ".fai") })
	layout, err := readLayout(p)
	if err != nil {
		t.Fatal(err)
	}
	return p, planRegions(layout, n, size)
}

func TestCacheBenchmark_ColdMissesThenWarmHits(t *testing.T) {
	for _, readAhead := range []bool{false, true} {
		p, regions := generatedRegions(t, 3, 20000, 24, 1000)
		phases := runCacheBenchmark(newManager(t, p, readAhead), regions)
		cold, warm := phases[0], phases[1]
		if cold.Errors != 0 || warm.Errors != 0 {
			t.Fatalf("read-ahead=%t: errors cold=%d warm=%d: %v %v", readAhead, cold.Errors, warm.Errors, cold.FirstErr, warm.FirstErr)
		}
		if cold.Misses != 24 || cold.Hits != 0 {
			t.Errorf("read-ahead=%t: cold pass hits=%d misses=%d, want 0 and 24", readAhead, cold.Hits, cold.Misses)
		}
		if warm.Hits != 24 || warm.Misses != 0 {
			t.Errorf("read-ahead=%t: warm pass hits=%d misses=%d, want 24 and 0", readAhead, warm.Hits, warm.Misses)
		}
	}
}

func TestKmerBenchmark_CountsEveryWindowAndCachesProfiles(t *testing.T) {
	const n, size, k = 12, 500, 6
	p, regions := generatedRegions(t, 2, 10000, n, size)
	phases := runKmerBenchmark(newManager(t, p, false), regions, k)
	build, cached := phases[0], phases[1]
	if build.Errors != 0 || cached.Errors != 0 {
		t.Fatalf("errors build=%d cached=%d: %v %v", build.Errors, cached.Errors, build.FirstErr, cached.FirstErr)
	}
	// The generated reference is pure ACGT, so every window is a k-mer.
	if want := int64(n * (size - k + 1)); build.Kmers != want || cached.Kmers != want {
		t.Errorf("k-mers build=%d cached=%d, want %d", build.Kmers, cached.Kmers, want)
	}
	if cached.Hits != n || cached.Misses != 0 {
		t.Errorf("cached pass hits=%d misses=%d, want %d and 0", cached.Hits, cached.Misses, n)
	}
}

func TestPrefetchBenchmark_ReadsEveryRegion(t *testing.T) {
	p, regions := generatedRegions(t, 2, 20000, 10, 1000)
	ph := runPrefetchBenchmark(newManager(t, p, false), regions)
	if ph.Errors != 0 {
		t.Fatalf("errors=%d: %v", ph.Errors, ph.FirstErr)
	}
	// Each region is loaded by the reader or the background prefetch, and
	// read once by the reader, so there are at least 10 lookups.
	if ph.Hits+ph.Misses < 10 {
		t.Errorf("hits+misses = %d, want at least 10", ph.Hits+ph.Misses)
	}
}

func TestRun_GeneratedReferenceIsRemoved(t *testing.T) {
	var out bytes.Buffer
	cfg := config{
		cacheSize: 64 << 20, k: 4, regions: 20, regionSize: 1000, benchmark: "all",
		prefetchWindow: 10000, genChroms: 2, genLength: 30000, seed: 1, verbose: true,
	}
	if err := run(cfg, &out); err != nil {
		t.Fatalf("run: %v\n%s", err, out.String())
	}
	s := out.String()
	for _, want := range []string{"cache: cold read", "cache: warm read", "kmer: build profile", "kmer: cached profile", "prefetch: read during PrefetchRegions"} {
		if !strings.Contains(s, want) {
			t.Errorf("output missing phase %q:\n%s", want, s)
		}
	}
	if strings.Contains(s, "first error") {
		t.Errorf("a phase reported errors:\n%s", s)
	}
	m := regexp.MustCompile(`generated reference: (\S+)`).FindStringSubmatch(s)
	if m == nil {
		t.Fatalf("output does not name the generated reference:\n%s", s)
	}
	if _, err := os.Stat(m[1]); !errors.Is(err, os.ErrNotExist) {
		t.Errorf("generated reference %s still exists after run (stat err %v)", m[1], err)
	}
}

func TestRun_UserReferenceWithOwnContigNames(t *testing.T) {
	var buf bytes.Buffer
	_ = writeTestReference(&buf, 1, 5000, 3)
	content := strings.Replace(buf.String(), ">chr1", ">scaffold_7", 1)
	cfg := config{reference: writeFile(t, content), cacheSize: 8 << 20, k: 4, regions: 10, regionSize: 400, benchmark: "cache"}
	var out bytes.Buffer
	if err := run(cfg, &out); err != nil {
		t.Fatalf("run: %v\n%s", err, out.String())
	}
	if !regexp.MustCompile(`cache: cold read\s+10\s+0\s`).MatchString(out.String()) {
		t.Errorf("want 10 regions and 0 errors in the cold pass:\n%s", out.String())
	}
}

func TestRun_RejectsBadConfig(t *testing.T) {
	base := config{cacheSize: 1 << 20, k: 8, regions: 10, regionSize: 100, benchmark: "all", genChroms: 1, genLength: 1000}
	cases := map[string]func(c *config){
		"unknown benchmark": func(c *config) { c.benchmark = "gpu" },
		"k too large":       func(c *config) { c.k = maxK + 1 },
		"no regions":        func(c *config) { c.regions = 0 },
		"empty regions":     func(c *config) { c.regionSize = 0 },
		"negative cache":    func(c *config) { c.cacheSize = -1 },
		"no generated data": func(c *config) { c.genLength = 0 },
		"missing reference": func(c *config) { c.reference = "/nonexistent/ref.fasta" },
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
