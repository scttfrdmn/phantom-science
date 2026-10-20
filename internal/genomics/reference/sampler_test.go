//go:build !release

package reference_test

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// fixtureContig describes one sequence of a test FASTA.
type fixtureContig struct {
	name   string
	length int
	nRuns  [][2]int // [start, end) runs of N
	soft   [][2]int // [start, end) soft-masked (lower-case) runs
}

// writeFixture writes a multi-line FASTA with the given contigs and returns
// its path and the upper-cased sequence of each contig.
func writeFixture(t testing.TB, dir string, contigs []fixtureContig, width int, eol string, withFai bool) (string, map[string][]byte) {
	t.Helper()
	r := rand.New(rand.NewPCG(7, 11))
	var buf bytes.Buffer
	want := make(map[string][]byte)
	for _, c := range contigs {
		seq := make([]byte, c.length)
		for i := range seq {
			seq[i] = "ACGT"[r.IntN(4)]
		}
		for _, run := range c.nRuns {
			for i := run[0]; i < run[1]; i++ {
				seq[i] = 'N'
			}
		}
		want[c.name] = append([]byte(nil), seq...)
		for _, run := range c.soft {
			for i := run[0]; i < run[1]; i++ {
				seq[i] = bytes.ToLower(seq[i : i+1])[0]
			}
		}
		fmt.Fprintf(&buf, ">%s description text%s", c.name, eol)
		for i := 0; i < len(seq); i += width {
			buf.Write(seq[i:min(i+width, len(seq))])
			buf.WriteString(eol)
		}
	}
	path := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(path, buf.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
	if withFai {
		p, err := reference.NewFastaParser(path)
		if err != nil {
			t.Fatal(err)
		}
		idx, err := p.BuildIndex()
		if err != nil {
			t.Fatal(err)
		}
		if err := p.SaveIndex(idx); err != nil {
			t.Fatal(err)
		}
	}
	return path, want
}

var standardContigs = []fixtureContig{
	{name: "chr1", length: 2000, nRuns: [][2]int{{0, 40}, {900, 1000}}, soft: [][2]int{{200, 450}, {1500, 1530}}},
	{name: "chr2", length: 1000, soft: [][2]int{{0, 1000}}},
	{name: "chrUn", length: 300, nRuns: [][2]int{{50, 200}}},
}

func newSampler(t testing.TB, path, region string) *reference.Sampler {
	t.Helper()
	s, err := reference.NewSampler(path, region)
	if err != nil {
		t.Fatalf("NewSampler(%q): %v", region, err)
	}
	return s
}

func revcomp(b []byte) []byte {
	out := make([]byte, len(b))
	for i, c := range b {
		out[len(b)-1-i] = map[byte]byte{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[c]
	}
	return out
}

type frag struct {
	chrom string
	start int64
	minus bool
	seq   string
}

func draw(t testing.TB, s *reference.Sampler, seed uint64, count, n int) []frag {
	t.Helper()
	r := rand.New(rand.NewPCG(seed, 99))
	out := make([]frag, count)
	for i := range out {
		c, st, m, seq, err := s.Fragment(r, n)
		if err != nil {
			t.Fatalf("Fragment: %v", err)
		}
		out[i] = frag{c, st, m, string(seq)}
	}
	return out
}

func TestSamplerDeterministicWithAndWithoutFai(t *testing.T) {
	noFai, _ := writeFixture(t, t.TempDir(), standardContigs, 60, "\n", false)
	withFai, _ := writeFixture(t, t.TempDir(), standardContigs, 60, "\n", true)
	if _, err := os.Stat(noFai + ".fai"); err == nil {
		t.Fatal("NewSampler must not need to write a .fai")
	}

	a := draw(t, newSampler(t, noFai, ""), 1, 300, 100)
	b := draw(t, newSampler(t, noFai, ""), 1, 300, 100)
	c := draw(t, newSampler(t, withFai, ""), 1, 300, 100)
	d := draw(t, newSampler(t, noFai, ""), 2, 300, 100)
	for i := range a {
		if a[i] != b[i] || a[i] != c[i] {
			t.Fatalf("fragment %d differs for the same seed: %+v / %+v / %+v", i, a[i], b[i], c[i])
		}
	}
	same := 0
	for i := range a {
		if a[i] == d[i] {
			same++
		}
	}
	if same > 3 {
		t.Errorf("seeds 1 and 2 gave %d identical fragments of 300", same)
	}
}

func TestSamplerFragmentsMatchFASTA(t *testing.T) {
	path, want := writeFixture(t, t.TempDir(), standardContigs, 60, "\n", true)
	s := newSampler(t, path, "")
	parser, err := reference.NewFastaParser(path)
	if err != nil {
		t.Fatal(err)
	}
	index, err := parser.LoadIndex()
	if err != nil {
		t.Fatal(err)
	}

	const n = 75
	for i, f := range draw(t, s, 3, 500, n) {
		fwd := want[f.chrom][f.start : f.start+n]
		expect := fwd
		if f.minus {
			expect = revcomp(fwd)
		}
		if f.seq != string(expect) {
			t.Fatalf("fragment %d %s:%d minus=%v\n got %s\nwant %s", i, f.chrom, f.start, f.minus, f.seq, expect)
		}
		gs, err := parser.GetSequence(index, f.chrom, f.start, f.start+n)
		if err != nil {
			t.Fatal(err)
		}
		if !bytes.Equal(bytes.ToUpper(gs), fwd) {
			t.Fatalf("fragment %d disagrees with GetSequence", i)
		}
		bs, err := s.Bases(f.chrom, f.start, f.start+n)
		if err != nil || !bytes.Equal(bs, fwd) {
			t.Fatalf("Bases(%s,%d) = %s, %v; want %s", f.chrom, f.start, bs, err, fwd)
		}
		if strings.Trim(f.seq, "ACGT") != "" {
			t.Fatalf("fragment %d contains a non-ACGT base: %s", i, f.seq)
		}
	}

	// Bases returns soft-masked bases upper-cased and N as N.
	bs, err := s.Bases("chr1", 890, 910)
	if err != nil || string(bs) != string(want["chr1"][890:910]) || !strings.HasSuffix(string(bs), strings.Repeat("N", 10)) {
		t.Errorf("Bases(chr1,890,910) = %s, %v", bs, err)
	}

	contigs := s.Contigs()
	if len(contigs) != 3 || contigs[0] != (reference.Contig{Name: "chr1", Length: 2000}) ||
		contigs[2] != (reference.Contig{Name: "chrUn", Length: 300}) {
		t.Errorf("Contigs() = %+v", contigs)
	}
}

func TestSamplerNeverReturnsN(t *testing.T) {
	// chrUn is half N; only windows of 50 bases fit beside the N run.
	path, _ := writeFixture(t, t.TempDir(), []fixtureContig{
		{name: "chrUn", length: 300, nRuns: [][2]int{{50, 250}}},
	}, 60, "\n", false)
	s := newSampler(t, path, "")
	for i, f := range draw(t, s, 4, 2000, 50) {
		if strings.Trim(f.seq, "ACGT") != "" || (f.start != 0 && f.start != 250) {
			t.Fatalf("fragment %d at %d: %s", i, f.start, f.seq)
		}
	}
	r := rand.New(rand.NewPCG(1, 1))
	if _, _, _, _, err := s.Fragment(r, 51); err == nil || !strings.Contains(err.Error(), "longer than every N-free stretch") {
		t.Errorf("Fragment(51) with no clean window: err = %v", err)
	}

	// One clean window among 600,000: the bounded redraw gives up.
	path, _ = writeFixture(t, t.TempDir(), []fixtureContig{
		{name: "gap", length: 600_000, nRuns: [][2]int{{50, 600_000}}},
	}, 60, "\n", false)
	if _, _, _, _, err := newSampler(t, path, "").Fragment(r, 50); err == nil || !strings.Contains(err.Error(), "free of N") {
		t.Errorf("Fragment with one clean window in 600k: err = %v", err)
	}
}

func TestSamplerStrandBalance(t *testing.T) {
	path, _ := writeFixture(t, t.TempDir(), standardContigs, 60, "\n", false)
	const draws = 4000
	minus := 0
	for _, f := range draw(t, newSampler(t, path, ""), 5, draws, 100) {
		if f.minus {
			minus++
		}
	}
	sigma := math.Sqrt(draws * 0.25)
	if math.Abs(float64(minus)-draws/2) > 4*sigma {
		t.Errorf("minus strand %d of %d; expected %d ± %.0f", minus, draws, draws/2, 4*sigma)
	}
	t.Logf("minus strand: %d / %d (4σ = %.0f)", minus, draws, 4*sigma)
}

func TestSamplerContigWeighting(t *testing.T) {
	path, want := writeFixture(t, t.TempDir(), standardContigs, 60, "\n", false)
	s := newSampler(t, path, "")
	nonN := map[string]float64{}
	var total float64
	for name, seq := range want {
		c := float64(len(seq) - bytes.Count(seq, []byte("N")))
		nonN[name] = c
		total += c
	}
	const draws = 20000
	got := map[string]int{}
	for _, f := range draw(t, s, 6, draws, 10) {
		got[f.chrom]++
	}
	var chi2 float64
	for name, w := range nonN {
		p := w / total
		mean, sd := draws*p, math.Sqrt(draws*p*(1-p))
		if math.Abs(float64(got[name])-mean) > 4*sd {
			t.Errorf("%s drawn %d times; expected %.0f ± %.0f (non-N share %.3f)", name, got[name], mean, 4*sd, p)
		}
		chi2 += (float64(got[name]) - mean) * (float64(got[name]) - mean) / mean
	}
	t.Logf("contig counts %v, non-N %v, chi2 = %.2f (df 2)", got, nonN, chi2)
	if chi2 > 13.8 { // p = 0.001 for 2 degrees of freedom
		t.Errorf("chi-square %.2f > 13.8", chi2)
	}

	// Contigs with no N-free stretch as long as the fragment are excluded,
	// not errors: chrUn's longest is 100 bases, chr1's 1000.
	for _, f := range draw(t, s, 7, 500, 120) {
		if f.chrom == "chrUn" {
			t.Fatalf("120-base fragment from chrUn at %d", f.start)
		}
	}
	for _, f := range draw(t, s, 7, 50, 1000) {
		if f.chrom == "chr1" && f.start != 1000 {
			t.Fatalf("1000-base fragment from chr1 at %d", f.start)
		}
	}
}

func TestSamplerRegion(t *testing.T) {
	path, want := writeFixture(t, t.TempDir(), standardContigs, 60, "\n", true)
	for _, region := range []string{"chr1:1,001-1,800", "1:1001-1800"} {
		s := newSampler(t, path, region)
		for i, f := range draw(t, s, 8, 1000, 100) {
			if f.chrom != "chr1" || f.start < 1000 || f.start+100 > 1800 {
				t.Fatalf("%s: fragment %d at %s:%d outside the region", region, i, f.chrom, f.start)
			}
			fwd := want["chr1"][f.start : f.start+100]
			if f.minus {
				fwd = revcomp(fwd)
			}
			if f.seq != string(fwd) {
				t.Fatalf("%s: fragment %d sequence mismatch", region, i)
			}
		}
		if b, err := s.Bases("chr1", 1000, 1800); err != nil || !bytes.Equal(b, want["chr1"][1000:1800]) {
			t.Errorf("%s: Bases over the whole region: %v", region, err)
		}
		for _, bad := range [][2]int64{{999, 1100}, {1700, 1801}} {
			if _, err := s.Bases("chr1", bad[0], bad[1]); err == nil {
				t.Errorf("%s: Bases(%v) outside the region succeeded", region, bad)
			}
		}
		if _, err := s.Bases("chr2", 0, 10); err == nil {
			t.Errorf("%s: Bases on a contig outside the region succeeded", region)
		}
		if _, _, _, _, err := s.Fragment(rand.New(rand.NewPCG(1, 1)), 801); err == nil {
			t.Errorf("%s: fragment longer than the region succeeded", region)
		}
		if len(s.Contigs()) != 3 {
			t.Errorf("%s: Contigs() should list every sequence, got %d", region, len(s.Contigs()))
		}
	}

	// Whole-contig and open-ended regions.
	s := newSampler(t, path, "chr2")
	if b, err := s.Bases("chr2", 0, 1000); err != nil || !bytes.Equal(b, want["chr2"]) {
		t.Errorf("region chr2: %v", err)
	}
	s = newSampler(t, path, "chr2:901")
	if b, err := s.Bases("chr2", 900, 1000); err != nil || !bytes.Equal(b, want["chr2"][900:]) {
		t.Errorf("region chr2:901: %v", err)
	}
	s = newSampler(t, path, "chr2:1-1")
	if b, err := s.Bases("chr2", 0, 1); err != nil || b[0] != want["chr2"][0] {
		t.Errorf("region chr2:1-1: %v", err)
	}
}

func TestSamplerChrAliasing(t *testing.T) {
	path, want := writeFixture(t, t.TempDir(), []fixtureContig{
		{name: "1", length: 500}, {name: "MT", length: 200}, {name: "X", length: 300},
	}, 70, "\n", false)
	s := newSampler(t, path, "")
	for _, q := range []struct {
		alias, name string
	}{{"chr1", "1"}, {"1", "1"}, {"chrM", "MT"}, {"MT", "MT"}, {"chrX", "X"}} {
		b, err := s.Bases(q.alias, 10, 60)
		if err != nil || !bytes.Equal(b, want[q.name][10:60]) {
			t.Errorf("Bases(%s): %v", q.alias, err)
		}
	}
	if _, err := s.Bases("chr2", 0, 10); err == nil {
		t.Error("Bases(chr2) on a FASTA without it succeeded")
	}
	for _, f := range draw(t, newSampler(t, path, "chrM:11-150"), 9, 100, 20) {
		if f.chrom != "MT" || f.start < 10 || f.start+20 > 150 {
			t.Fatalf("region chrM:11-150 gave %s:%d", f.chrom, f.start)
		}
	}

	pathChr, wantChr := writeFixture(t, t.TempDir(), []fixtureContig{{name: "chr1", length: 500}, {name: "chrM", length: 100}}, 60, "\n", true)
	s = newSampler(t, pathChr, "1:100-200")
	if b, err := s.Bases("1", 99, 200); err != nil || !bytes.Equal(b, wantChr["chr1"][99:200]) {
		t.Errorf("unprefixed query on chr-named FASTA: %v", err)
	}
	if b, err := newSampler(t, pathChr, "").Bases("MT", 0, 100); err != nil || !bytes.Equal(b, wantChr["chrM"]) {
		t.Errorf("Bases(MT) on chrM: %v", err)
	}
}

func TestSamplerFileFormats(t *testing.T) {
	contigs := []fixtureContig{
		{name: "long", length: 200_000, soft: [][2]int{{100, 70_000}}, nRuns: [][2]int{{150_000, 150_500}}},
		{name: "short", length: 77},
	}
	for _, tc := range []struct {
		name  string
		width int
		eol   string
	}{
		{"single-line contigs longer than 64 KiB", 1 << 20, "\n"},
		{"CRLF line endings", 61, "\r\n"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			path, want := writeFixture(t, t.TempDir(), contigs, tc.width, tc.eol, false)
			s := newSampler(t, path, "")
			for name, seq := range want {
				b, err := s.Bases(name, 0, int64(len(seq)))
				if err != nil || !bytes.Equal(b, seq) {
					t.Fatalf("Bases(%s) whole contig: %v", name, err)
				}
			}
			r := newSampler(t, path, "long:140,001-160,000")
			if b, err := r.Bases("long", 140_000, 160_000); err != nil || !bytes.Equal(b, want["long"][140_000:160_000]) {
				t.Fatalf("region load: %v", err)
			}
		})
	}
}

func TestSamplerConcurrent(t *testing.T) {
	path, _ := writeFixture(t, t.TempDir(), standardContigs, 60, "\n", true)
	s := newSampler(t, path, "")
	const workers, per = 8, 500
	want := make([][]frag, workers)
	for w := range want {
		want[w] = draw(t, s, uint64(100+w), per, 120)
	}
	got := make([][]frag, workers)
	errs := make([]error, workers)
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func(w int) {
			defer wg.Done()
			r := rand.New(rand.NewPCG(uint64(100+w), 99))
			for i := 0; i < per; i++ {
				c, st, m, seq, err := s.Fragment(r, 120)
				if err != nil {
					errs[w] = err
					return
				}
				got[w] = append(got[w], frag{c, st, m, string(seq)})
				if _, err := s.Bases(c, st, st+120); err != nil {
					errs[w] = err
					return
				}
			}
		}(w)
	}
	wg.Wait()
	for w := range got {
		if errs[w] != nil {
			t.Fatalf("worker %d: %v", w, errs[w])
		}
		for i := range got[w] {
			if got[w][i] != want[w][i] {
				t.Fatalf("worker %d fragment %d differs from the sequential run", w, i)
			}
		}
	}
}

func TestSamplerErrors(t *testing.T) {
	dir := t.TempDir()
	path, _ := writeFixture(t, dir, standardContigs, 60, "\n", false)

	expectErr := func(what, p, region, substr string) {
		t.Helper()
		_, err := reference.NewSampler(p, region)
		if err == nil || !strings.Contains(err.Error(), substr) {
			t.Errorf("%s: err = %v, want it to mention %q", what, err, substr)
		}
	}
	expectErr("missing file", filepath.Join(dir, "absent.fa"), "", "no such file")

	gz := filepath.Join(dir, "ref.fa.gz")
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	raw, _ := os.ReadFile(path)
	zw.Write(raw)
	zw.Close()
	if err := os.WriteFile(gz, buf.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
	expectErr("gzip", gz, "", "only uncompressed FASTA")

	expectErr("unknown contig", path, "chr9:1-10", "not found")
	expectErr("unknown contig name", path, "chr9", "not found")
	expectErr("end past contig", path, "chr2:1-1001", "past the end")
	expectErr("end before start", path, "chr2:50-10", "before start")
	expectErr("zero start", path, "chr2:0-10", "positive")
	expectErr("all-N region", path, "chr1:901-1000", "no A/C/G/T")

	s := newSampler(t, path, "")
	r := rand.New(rand.NewPCG(1, 2))
	if _, _, _, _, err := s.Fragment(r, 1001); err == nil || !strings.Contains(err.Error(), "longer than") {
		t.Errorf("fragment longer than every contig: %v", err)
	}
	if _, _, _, _, err := s.Fragment(r, 0); err == nil {
		t.Error("Fragment(0) succeeded")
	}
	for _, bad := range [][2]int64{{-1, 10}, {10, 10}, {1990, 2001}} {
		if _, err := s.Bases("chr1", bad[0], bad[1]); err == nil {
			t.Errorf("Bases(chr1, %v) succeeded", bad)
		}
	}

	ragged := filepath.Join(t.TempDir(), "ragged.fa")
	os.WriteFile(ragged, []byte(">a\nACGT\nAC\nACGT\n"), 0o600)
	expectErr("ragged lines", ragged, "", "different lengths")

	// A .fai that does not describe the file (offsets shifted by one).
	stale := filepath.Join(t.TempDir(), "stale.fa")
	os.WriteFile(stale, []byte(">a\nACGTACGT\nACGT\n>b\nGGGG\n"), 0o600)
	os.WriteFile(stale+".fai", []byte("a\t12\t2\t8\t9\nb\t4\t20\t4\t5\n"), 0o600)
	expectErr("stale .fai", stale, "", "does not match")
	os.WriteFile(stale+".fai", []byte("a\t12\t3\t8\t9\nb\t4\t20\t4\t5\n"), 0o600)
	if b, err := newSampler(t, stale, "").Bases("b", 0, 4); err != nil || string(b) != "GGGG" {
		t.Errorf("correct .fai: %s, %v", b, err)
	}
}

// BenchmarkSamplerFragment draws 150-base fragments from a 10 Mb contig.
func BenchmarkSamplerFragment(b *testing.B) {
	path, _ := writeFixture(b, b.TempDir(), []fixtureContig{
		{name: "chr1", length: 10_000_000, nRuns: [][2]int{{0, 10_000}, {5_000_000, 5_050_000}}, soft: [][2]int{{1_000_000, 2_000_000}}},
	}, 60, "\n", false)
	s := newSampler(b, path, "")
	r := rand.New(rand.NewPCG(1, 2))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, _, _, err := s.Fragment(r, 150); err != nil {
			b.Fatal(err)
		}
	}
}
