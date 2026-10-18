//go:build !release

package transcriptome

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"reflect"
	"regexp"
	"slices"
	"sort"
	"strconv"
	"strings"
	"testing"
)

// --- fixtures ---------------------------------------------------------------

// Hand-checked transcripts of testdata/tiny_genome.fna + tiny.gtf.
//
// chr1 = ACGTACGTAC TTTTTTTTTT GATTACAGAT CCCCCCCCCC TGCATGCAAA GGGGGGGGGG
//
//	1-10       11-20      21-30      31-40      41-50      51-60
const (
	// TXA (+): exons 1-10, 21-30, 41-50.
	wantTXA = "ACGTACGTAC" + "GATTACAGAT" + "TGCATGCAAA"
	// TXC (+): exons 1-10, 41-50 (TXA without exon 2).
	wantTXC = "ACGTACGTAC" + "TGCATGCAAA"
	// TXB (−): exons 43-49 (CATGCAA), 22-28 (ATTACAG), 3-8 (GTACGT) in
	// transcript order; each reverse-complemented:
	// rc(CATGCAA)=TTGCATG, rc(ATTACAG)=CTGTAAT, rc(GTACGT)=ACGTAC.
	wantTXB = "TTGCATG" + "CTGTAAT" + "ACGTAC"
	// TXB without its exon 2.
	wantTXBSkip2 = "TTGCATG" + "ACGTAC"
	// TXD (−, single exon chr2:5-64) = rc(GGATCCAAGCTTGAATTCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAA).
	wantTXD = "TTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGAATTCAAGCTTGGATCC"
)

func fixtureTable(t *testing.T) (*Table, *Annotation) {
	t.Helper()
	ann, err := ParseGTF(filepath.Join("testdata", "tiny.gtf"))
	if err != nil {
		t.Fatal(err)
	}
	tab, st, err := BuildFromGenome(filepath.Join("testdata", "tiny_genome.fna"), ann)
	if err != nil {
		t.Fatal(err)
	}
	if st.MissingChrom != 1 {
		t.Errorf("MissingChrom = %d, want 1 (TXE on chr3)", st.MissingChrom)
	}
	return tab, ann
}

// randomTable builds n transcripts of random sequence with lengths in
// [minLen, maxLen]; with exons, each gets 1–6 exons of random lengths on
// alternating strands with consistent genomic coordinates.
func randomTable(t *testing.T, seed uint64, n, minLen, maxLen int, exons bool) *Table {
	t.Helper()
	r := rand.New(rand.NewPCG(seed, 99))
	var ts []*Transcript
	pos := int64(1000)
	for i := range n {
		L := minLen + r.IntN(maxLen-minLen+1)
		seq := make([]byte, L)
		for j := range seq {
			seq[j] = "ACGT"[r.IntN(4)]
		}
		tx := &Transcript{ID: "T" + strconv.Itoa(i), GeneID: "G" + strconv.Itoa(i/2), Seq: seq}
		if exons {
			k := 1 + r.IntN(6)
			k = min(k, L/20)
			// Cut L into k pieces of at least 10 bases.
			cuts := []int{0, L}
			for len(cuts) < k+1 {
				c := 10 + r.IntN(L-19)
				ok := true
				for _, x := range cuts {
					if abs(x-c) < 10 {
						ok = false
					}
				}
				if ok {
					cuts = append(cuts, c)
				}
			}
			sort.Ints(cuts)
			tx.Chrom = "chrS"
			tx.Strand = "+-"[i%2]
			var ex []Exon
			for j := 1; j < len(cuts); j++ {
				l := int64(cuts[j] - cuts[j-1])
				ex = append(ex, Exon{Start: pos, End: pos + l - 1})
				pos += l + 100
			}
			if tx.Strand == '-' {
				slices.Reverse(ex)
			}
			tx.Exons = ex
		}
		ts = append(ts, tx)
	}
	tab, err := NewTable(ts)
	if err != nil {
		t.Fatal(err)
	}
	return tab
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func noErrConfig(reads int, seed uint64) Config {
	c := DefaultConfig(reads, seed)
	c.Tier = TierNone
	c.Run.Channels = 512
	return c
}

func collect(t *testing.T, s *Simulator) []*Read {
	t.Helper()
	var out []*Read
	if err := s.Generate(func(r *Read) error { out = append(out, r); return nil }); err != nil {
		t.Fatal(err)
	}
	return out
}

// --- parsing ---------------------------------------------------------------

func TestParseFASTAHeader(t *testing.T) {
	for _, c := range []struct{ h, id, gene string }{
		{"ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|", "ENST00000456328.2", "ENSG00000290825.1"},
		{"ENST00000632684.1 cdna chromosome:GRCh38:7:142786213:142786624:1 gene:ENSG00000282431.1 gene_biotype:TR_D_gene", "ENST00000632684.1", "ENSG00000282431.1"},
		{"tx1 some description", "tx1", "tx1"},
		{"tx2", "tx2", "tx2"},
	} {
		id, gene := ParseFASTAHeader(c.h)
		if id != c.id || gene != c.gene {
			t.Errorf("ParseFASTAHeader(%q) = %q, %q; want %q, %q", c.h, id, gene, c.id, c.gene)
		}
	}
}

func TestBuildFromGenome_HandChecked(t *testing.T) {
	tab, ann := fixtureTable(t)
	if ann.Skipped != 1 {
		t.Errorf("Skipped = %d, want 1 (TXU has strand '.')", ann.Skipped)
	}
	want := map[string]string{"TXA": wantTXA, "TXC": wantTXC, "TXB": wantTXB, "TXD": wantTXD}
	var ids []string
	for _, tx := range tab.Transcripts {
		ids = append(ids, tx.ID)
	}
	if !reflect.DeepEqual(ids, []string{"TXA", "TXC", "TXB", "TXD"}) {
		t.Fatalf("transcripts %v, want GTF order TXA TXC TXB TXD", ids)
	}
	for id, seq := range want {
		tx := tab.Lookup(id)
		if string(tx.Seq) != seq {
			t.Errorf("%s = %s, want %s", id, tx.Seq, seq)
		}
	}
	b := tab.Lookup("TXB")
	if b.Strand != '-' || b.GeneID != "G2" || !reflect.DeepEqual(b.Exons, []Exon{{43, 49}, {22, 28}, {3, 8}}) {
		t.Errorf("TXB structure = %c %s %v, want - G2 [43-49 22-28 3-8]", b.Strand, b.GeneID, b.Exons)
	}
	a := tab.Lookup("TXA")
	if !reflect.DeepEqual(a.Exons, []Exon{{1, 10}, {21, 30}, {41, 50}}) {
		t.Errorf("TXA exons %v not sorted ascending", a.Exons)
	}
}

func TestBuildFromGenome_ExonPastEnd(t *testing.T) {
	dir := t.TempDir()
	g := filepath.Join(dir, "g.fa")
	os.WriteFile(g, []byte(">c\nACGTACGT\n"), 0o600)
	gtf := filepath.Join(dir, "a.gtf")
	os.WriteFile(gtf, []byte("c\tx\texon\t5\t9\t.\t+\t.\ttranscript_id \"t\";\n"), 0o600)
	ann, err := ParseGTF(gtf)
	if err != nil {
		t.Fatal(err)
	}
	if _, _, err := BuildFromGenome(g, ann); err == nil || !strings.Contains(err.Error(), "past chromosome end") {
		t.Errorf("err = %v, want past chromosome end", err)
	}
}

func TestParseGTF_Errors(t *testing.T) {
	for name, gtf := range map[string]string{
		"overlap":      "c\tx\texon\t1\t10\t.\t+\t.\ttranscript_id \"t\";\nc\tx\texon\t10\t20\t.\t+\t.\ttranscript_id \"t\";\n",
		"mixed strand": "c\tx\texon\t1\t10\t.\t+\t.\ttranscript_id \"t\";\nc\tx\texon\t20\t30\t.\t-\t.\ttranscript_id \"t\";\n",
		"no id":        "c\tx\texon\t1\t10\t.\t+\t.\tgene_id \"g\";\n",
		"columns":      "c\tx\texon\t1\t10\n",
		"coords":       "c\tx\texon\t10\t1\t.\t+\t.\ttranscript_id \"t\";\n",
	} {
		if _, err := parseGTF(strings.NewReader(gtf)); err == nil {
			t.Errorf("%s: no error", name)
		}
	}
}

func TestTranscriptFASTA_GzipAndAnnotation(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "tx.fa.gz")
	var buf bytes.Buffer
	gz := gzip.NewWriter(&buf)
	// TXA split over lines and lower case; TXB with the wrong length; U read as T.
	gz.Write([]byte(">TXA|G1|x\n" + strings.ToLower(wantTXA[:13]) + "\n" + wantTXA[13:] + "\n>TXB\nACGU\n>TXZ\nACGT\n"))
	gz.Close()
	os.WriteFile(p, buf.Bytes(), 0o600)
	tab, err := LoadTranscriptFASTA(p)
	if err != nil {
		t.Fatal(err)
	}
	if got := string(tab.Lookup("TXA").Seq); got != wantTXA {
		t.Errorf("TXA = %s", got)
	}
	if got := string(tab.Lookup("TXB").Seq); got != "ACGT" {
		t.Errorf("TXB = %s, want U read as T", got)
	}
	_, ann := fixtureTable(t)
	att, mis := tab.AttachAnnotation(ann)
	if att != 1 || mis != 1 {
		t.Errorf("attached %d, mismatched %d; want 1, 1", att, mis)
	}
	if a := tab.Lookup("TXA"); a.Strand != '+' || len(a.Exons) != 3 {
		t.Errorf("TXA not annotated: %c %v", a.Strand, a.Exons)
	}
	if b := tab.Lookup("TXB"); b.HasStructure() {
		t.Error("TXB with mismatched length got a structure")
	}

	dup := filepath.Join(dir, "dup.fa")
	os.WriteFile(dup, []byte(">a\nAC\n>a\nGT\n"), 0o600)
	if _, err := LoadTranscriptFASTA(dup); err == nil {
		t.Error("duplicate IDs accepted")
	}
}

func TestParseProfile(t *testing.T) {
	tab := randomTable(t, 1, 4, 50, 60, false) // T0..T3
	// Version suffixes: rename T3 to a versioned ID.
	tab.Transcripts[3].ID = "T3.2"
	delete(tab.byID, "T3")
	tab.byID["T3.2"] = 3
	cases := []struct {
		name, in string
		want     []float64
		unmatch  int
	}{
		{"plain", "T0\t1\nT1\t3\nX\t5\n", []float64{0.25e6, 0.75e6, 0, 0}, 1},
		{"header", "# comment\ntranscript_id\ttpm\nT0\t2\nT2\t2\n", []float64{0.5e6, 0, 0.5e6, 0}, 0},
		{"salmon", "Name\tLength\tEffectiveLength\tTPM\tNumReads\nT1\t50\t40\t10\t999\nT3\t50\t40\t30\t1\n", []float64{0, 0.25e6, 0, 0.75e6}, 0},
		{"bambu", "TXNAME\tGENEID\ts1\ts2\nT0\tG0\t1\t1\nT3.7\tG1\t3\t3\nBambuTx1\tG1\t9\t9\n", []float64{0.25e6, 0, 0, 0.75e6}, 1},
	}
	for _, c := range cases {
		got, st, err := ParseProfile(strings.NewReader(c.in), tab)
		if err != nil {
			t.Errorf("%s: %v", c.name, err)
			continue
		}
		for i := range got {
			if math.Abs(got[i]-c.want[i]) > 1e-6 {
				t.Errorf("%s: tpm = %v, want %v", c.name, got, c.want)
				break
			}
		}
		if st.Unmatched != c.unmatch {
			t.Errorf("%s: unmatched = %d, want %d", c.name, st.Unmatched, c.unmatch)
		}
	}
	for _, bad := range []string{"T0\t-1\n", "T0\t1\nT1\tabc\n", "X\t1\n", "T0\t0\n"} {
		if _, _, err := ParseProfile(strings.NewReader(bad), tab); err == nil {
			t.Errorf("ParseProfile(%q): no error", bad)
		}
	}
}

func TestConfigValidate(t *testing.T) {
	tab := randomTable(t, 1, 3, 100, 200, false)
	tpm := []float64{1, 1, 1}
	for name, mut := range map[string]func(*Config){
		"reads":   func(c *Config) { c.Reads = 0 },
		"seed":    func(c *Config) { c.Seed = 0 },
		"full":    func(c *Config) { c.FullLengthFraction = 1.5 },
		"trunc":   func(c *Config) { c.TruncMean = 0 },
		"minfrag": func(c *Config) { c.MinFragment = 0 },
		"polya":   func(c *Config) { c.PolyAMin, c.PolyAMax = 50, 20 },
		"novel":   func(c *Config) { c.NovelRate = -0.1 },
		"tier":    func(c *Config) { c.Tier = "ultra" },
		"kit":     func(c *Config) { c.Kit.SSP = "ACGN" },
	} {
		c := DefaultConfig(10, 1)
		mut(&c)
		if _, err := NewSimulator(tab, tpm, c); err == nil {
			t.Errorf("%s: no error", name)
		}
	}
	if _, err := NewSimulator(tab, []float64{0, 0, 0}, DefaultConfig(10, 1)); err == nil {
		t.Error("zero expression accepted")
	}
}

// --- novel isoforms --------------------------------------------------------

func TestDesignNovel_Fixture(t *testing.T) {
	tab, _ := fixtureTable(t)
	c := noErrConfig(10, 7)
	c.NovelRate = 0.5
	s, err := NewSimulator(tab, []float64{1, 1, 1, 1}, c)
	if err != nil {
		t.Fatal(err)
	}
	novel := map[string]*Species{}
	for _, sp := range s.Species() {
		if sp.Novel != NovelNone {
			novel[sp.Parent.ID] = sp
		}
	}
	// TXA's only internal exon skip equals annotated TXC, and its last exon
	// (10 nt) is too short for an alternative 3′ end; TXC has 2 exons and a
	// 10 nt last exon.
	if novel["TXA"] != nil || novel["TXC"] != nil {
		t.Errorf("unexpected novel isoforms for TXA/TXC: %+v %+v", novel["TXA"], novel["TXC"])
	}
	b := novel["TXB"]
	if b == nil || b.Novel != NovelSkipExon || string(b.Seq) != wantTXBSkip2 || b.ID != "TXB_skip2" ||
		!reflect.DeepEqual(b.Exons, []Exon{{43, 49}, {3, 8}}) || b.IntronChain() != "9-42" {
		t.Fatalf("TXB novel = %+v, want skip of exon 2 = %s with intron 9-42", b, wantTXBSkip2)
	}
	d := novel["TXD"]
	if d == nil || d.Novel != NovelAlt3Prime {
		t.Fatalf("TXD novel = %+v, want alt_3prime", d)
	}
	n := len(d.Seq)
	if n < 12 || n > 48 || string(d.Seq) != wantTXD[:n] || d.Detail != "end="+strconv.Itoa(n) {
		t.Errorf("TXD alt3 = %s (%s), want a prefix of %s keeping 12..48 bases", d.Seq, d.Detail, wantTXD)
	}
	// Minus strand: shortening the 3′ end raises the exon's genomic start.
	if !reflect.DeepEqual(d.Exons, []Exon{{64 - int64(n) + 1, 64}}) {
		t.Errorf("TXD alt3 exons = %v, want [%d-64]", d.Exons, 64-n+1)
	}
	// TPMs split by the novel rate and still sum to 10^6.
	sum := 0.0
	for _, sp := range s.Species() {
		sum += sp.TPM
	}
	if math.Abs(sum-1e6) > 1e-6 || math.Abs(b.TPM-0.125e6) > 1e-6 {
		t.Errorf("TPM sum %g (want 1e6), TXB_skip2 %g (want 125000)", sum, b.TPM)
	}
}

// senseFragment strips the kit primers and poly(A) from an error-free read,
// returning the transcript part, or fails.
func senseFragment(t *testing.T, s *Simulator, rd *Read) []byte {
	t.Helper()
	seq := rd.Record.Seq
	if rd.Truth.Reverse {
		seq = ReverseComplement(seq)
	}
	k := s.cfg.Kit
	pre := k.Adapter5 + k.SSP
	suf := strings.Repeat("A", rd.Truth.PolyA) + string(ReverseComplement([]byte(k.VNP))) + k.Adapter3
	if !bytes.HasPrefix(seq, []byte(pre)) || !bytes.HasSuffix(seq, []byte(suf)) || len(seq) < len(pre)+len(suf) {
		t.Fatalf("read %s lacks primer/poly(A) structure: %s", rd.Record.Name, seq)
	}
	return seq[len(pre) : len(seq)-len(suf)]
}

// expectedIsoform rebuilds a species' sequence independently from its
// parent and Detail.
func expectedIsoform(t *testing.T, sp *Species) []byte {
	t.Helper()
	p := sp.Parent
	switch sp.Novel {
	case NovelNone:
		return p.Seq
	case NovelSkipExon:
		k, _ := strconv.Atoi(strings.TrimPrefix(sp.Detail, "skipped_exon="))
		var out []byte
		off := 0
		for i, e := range p.Exons {
			if i != k-1 {
				out = append(out, p.Seq[off:off+e.Len()]...)
			}
			off += e.Len()
		}
		if k < 2 || k > len(p.Exons)-1 {
			t.Fatalf("%s skips non-internal exon %d of %d", sp.ID, k, len(p.Exons))
		}
		return out
	case NovelAlt3Prime:
		n, _ := strconv.Atoi(strings.TrimPrefix(sp.Detail, "end="))
		last := len(p.Seq)
		if p.HasStructure() {
			last = p.Exons[len(p.Exons)-1].Len()
		}
		if keep := n - (len(p.Seq) - last); keep < (2*last+9)/10 || keep > last*8/10 {
			t.Fatalf("%s keeps %d of a %d-base last exon", sp.ID, keep, last)
		}
		return p.Seq[:n]
	}
	t.Fatalf("unknown novel type %q", sp.Novel)
	return nil
}

func TestReads_ErrorFreeStructure(t *testing.T) {
	tab := randomTable(t, 3, 60, 150, 900, true)
	c := noErrConfig(4000, 11)
	c.NovelRate = 0.3
	c.Kit.Adapter5, c.Kit.Adapter3 = "GGTTAA", "CCA"
	tpm := DefaultTPM(tab, ZipfModel{Alpha: 0.5, Sigma: 0.3}, c.Seed)
	s, err := NewSimulator(tab, tpm, c)
	if err != nil {
		t.Fatal(err)
	}
	var full, skip, alt3, rev int
	for _, rd := range collect(t, s) {
		tr := rd.Truth
		frag := senseFragment(t, s, rd)
		iso := expectedIsoform(t, tr.Species)
		if !bytes.Equal(iso, tr.Species.Seq) {
			t.Fatalf("%s sequence differs from its definition", tr.Species.ID)
		}
		if tr.End != len(iso) || tr.Start != tr.Trunc5p || !bytes.Equal(frag, iso[tr.Start:]) {
			t.Fatalf("read %s is not the 3′ suffix [%d:] of %s", rd.Record.Name, tr.Start, tr.Species.ID)
		}
		if tr.FullLength != (tr.Start == 0) {
			t.Fatalf("full_length %v with start %d", tr.FullLength, tr.Start)
		}
		if !bytes.Equal(rd.Record.Seq, rd.Construct) || len(rd.Record.Qual) != len(rd.Record.Seq) {
			t.Fatal("tier none changed the construct")
		}
		if tr.FullLength {
			full++
		}
		switch tr.Species.Novel {
		case NovelSkipExon:
			skip++
		case NovelAlt3Prime:
			alt3++
		}
		if tr.Reverse {
			rev++
		}
	}
	if skip == 0 || alt3 == 0 || full == 0 {
		t.Errorf("coverage: %d skip, %d alt3, %d full-length reads", skip, alt3, full)
	}
	if d := math.Abs(float64(rev) - 2000); d > 4*math.Sqrt(4000*0.25) {
		t.Errorf("%d of 4000 reads reverse, want 2000 ± %.0f", rev, 4*math.Sqrt(1000))
	}
}

// Exon-skip reads from the genome-built fixture equal the hand-built
// concatenation.
func TestReads_FixtureMinusStrandSkip(t *testing.T) {
	tab, _ := fixtureTable(t)
	c := noErrConfig(400, 5)
	c.NovelRate = 0.5
	c.FullLengthFraction = 1
	s, err := NewSimulator(tab, []float64{0, 0, 1, 0}, c) // TXB only
	if err != nil {
		t.Fatal(err)
	}
	n := map[string]int{}
	for _, rd := range collect(t, s) {
		f := string(senseFragment(t, s, rd))
		switch rd.Truth.Species.Novel {
		case NovelSkipExon:
			if f != wantTXBSkip2 {
				t.Fatalf("skip read %s, want %s", f, wantTXBSkip2)
			}
		case NovelNone:
			if f != wantTXB {
				t.Fatalf("TXB read %s, want %s", f, wantTXB)
			}
		}
		n[rd.Truth.Species.ID]++
	}
	// Binomial(400, 0.5) within 4σ.
	if d := math.Abs(float64(n["TXB_skip2"]) - 200); d > 4*10 || n["TXB"]+n["TXB_skip2"] != 400 {
		t.Errorf("counts %v, want TXB_skip2 ≈ 200 ± 40 and nothing else", n)
	}
}

// --- statistics ------------------------------------------------------------

func TestCounts_MultinomialInTPM(t *testing.T) {
	tab := randomTable(t, 5, 40, 120, 200, true)
	const N = 20000
	c := noErrConfig(N, 21)
	c.NovelRate = 0.2
	tpm := DefaultTPM(tab, ZipfModel{Alpha: 1, Sigma: 0.5}, c.Seed)
	s, err := NewSimulator(tab, tpm, c)
	if err != nil {
		t.Fatal(err)
	}
	counts := map[*Species]int{}
	for _, rd := range collect(t, s) {
		counts[rd.Truth.Species]++
	}
	worst := 0.0
	for _, sp := range s.Species() {
		p := sp.TPM / 1e6
		mu, sd := N*p, math.Sqrt(N*p*(1-p))
		z := (float64(counts[sp]) - mu) / max(sd, 1e-9)
		worst = max(worst, math.Abs(z))
		if math.Abs(z) > 4 {
			t.Errorf("%s: %d reads, want %.1f ± 4·%.1f", sp.ID, counts[sp], mu, sd)
		}
	}
	t.Logf("%d species, max |z| = %.2f", len(s.Species()), worst)
}

func TestZipf_RankFrequencySlope(t *testing.T) {
	for _, m := range []ZipfModel{{1, 0.5}, {0.8, 1}, {1.2, 0}} {
		const n = 20000
		tpm, rank := m.TPM(rand.New(rand.NewPCG(3, 4)), n)
		sum := 0.0
		var sx, sy, sxx, sxy float64
		for i := range tpm {
			sum += tpm[i]
			x, y := math.Log(float64(rank[i])), math.Log(tpm[i])
			sx, sy, sxx, sxy = sx+x, sy+y, sxx+x*x, sxy+x*y
		}
		if math.Abs(sum-1e6) > 1e-3 {
			t.Errorf("TPM sum %g", sum)
		}
		Sxx := sxx - sx*sx/n
		slope := (sxy - sx*sy/n) / Sxx
		se := m.Sigma / math.Sqrt(Sxx)
		if math.Abs(slope+m.Alpha) > 4*se+1e-9 {
			t.Errorf("%+v: slope %.4f, want %.4f ± 4·%.4f", m, slope, -m.Alpha, se)
		}
		t.Logf("%+v: slope %.4f (SE %.4f)", m, slope, se)
	}
}

// The default model puts half of all molecules in a few hundred
// transcripts of a GENCODE-sized table.
func TestZipf_HalfMassInFewHundred(t *testing.T) {
	m := ZipfModel{DefaultZipfAlpha, DefaultExpressionSigma}
	tpm, _ := m.TPM(rand.New(rand.NewPCG(9, 9)), 250000)
	sort.Sort(sort.Reverse(sort.Float64Slice(tpm)))
	acc, k := 0.0, 0
	for acc < 0.5e6 {
		acc += tpm[k]
		k++
	}
	if k < 100 || k > 1000 {
		t.Errorf("top %d transcripts carry half the TPM, want a few hundred", k)
	}
	t.Logf("half the TPM in the top %d of 250000", k)
}

// truncGeomMoments returns the exact mean and variance of P(k) ∝ q^k on
// 0..K, and P(k < x).
func truncGeomMoments(q float64, K, x int) (mean, variance, below float64) {
	var z, m1, m2, b float64
	w := 1.0
	for k := 0; k <= K; k++ {
		z += w
		m1 += w * float64(k)
		m2 += w * float64(k) * float64(k)
		if k < x {
			b += w
		}
		w *= q
	}
	mean = m1 / z
	return mean, m2/z - mean*mean, b / z
}

func TestFragmentsAndPolyA(t *testing.T) {
	tab := randomTable(t, 8, 1, 3000, 3000, false)
	const N = 20000
	c := noErrConfig(N, 31)
	c.NovelRate = 0
	s, err := NewSimulator(tab, []float64{1}, c)
	if err != nil {
		t.Fatal(err)
	}
	L := 3000
	var full int
	var ks, pas []float64
	for _, rd := range collect(t, s) {
		tr := rd.Truth
		if tr.PolyA < c.PolyAMin || tr.PolyA > c.PolyAMax {
			t.Fatalf("poly(A) %d outside [%d, %d]", tr.PolyA, c.PolyAMin, c.PolyAMax)
		}
		pas = append(pas, float64(tr.PolyA-c.PolyAMin))
		if tr.FullLength {
			full++
			continue
		}
		keep := L - tr.Trunc5p
		if keep < c.MinFragment || keep > L-1 {
			t.Fatalf("truncated fragment keeps %d of %d bases", keep, L)
		}
		ks = append(ks, float64(keep-c.MinFragment))
	}
	// Full-length fraction.
	F := c.FullLengthFraction
	if d := math.Abs(float64(full)/N - F); d > 4*math.Sqrt(F*(1-F)/N) {
		t.Errorf("full-length fraction %.4f, want %.2f ± %.4f", float64(full)/N, F, 4*math.Sqrt(F*(1-F)/N))
	}
	check := func(name string, xs []float64, q float64, K, x int) {
		mean, v, below := truncGeomMoments(q, K, x)
		n := float64(len(xs))
		got := 0.0
		nb := 0
		for _, k := range xs {
			got += k
			if k < float64(x) {
				nb++
			}
		}
		got /= n
		if math.Abs(got-mean) > 4*math.Sqrt(v/n) {
			t.Errorf("%s mean %.2f, want %.2f ± %.2f", name, got, mean, 4*math.Sqrt(v/n))
		}
		pb := float64(nb) / n
		if math.Abs(pb-below) > 4*math.Sqrt(below*(1-below)/n) {
			t.Errorf("%s P(k < %d) = %.4f, want %.4f", name, x, pb, below)
		}
		t.Logf("%s: mean %.2f (model %.2f), P(k<%d) %.4f (model %.4f)", name, got, mean, x, pb, below)
	}
	check("truncation excess", ks, math.Exp(-1/c.TruncMean), L-1-c.MinFragment, 300)
	check("poly(A) excess", pas, geomQ(c.PolyAMean), c.PolyAMax-c.PolyAMin, 25)
}

func TestFragment_ShortTranscriptsFullLength(t *testing.T) {
	tab := randomTable(t, 8, 1, 80, 80, false) // shorter than MinFragment + 1
	c := noErrConfig(200, 2)
	c.FullLengthFraction = 0
	s, err := NewSimulator(tab, []float64{1}, c)
	if err != nil {
		t.Fatal(err)
	}
	for _, rd := range collect(t, s) {
		if !rd.Truth.FullLength || rd.Truth.Start != 0 {
			t.Fatal("transcript shorter than the minimum fragment was truncated")
		}
	}
}

func TestTruncGeom_Exact(t *testing.T) {
	r := rand.New(rand.NewPCG(1, 2))
	q, K := 0.7, 5
	const n = 100000
	var h [6]int
	for range n {
		h[truncGeom(r, q, K)]++
	}
	z := (1 - math.Pow(q, float64(K+1))) / (1 - q)
	for k, c := range h {
		p := math.Pow(q, float64(k)) / z
		if math.Abs(float64(c)/n-p) > 4*math.Sqrt(p*(1-p)/n) {
			t.Errorf("P(%d) = %.4f, want %.4f", k, float64(c)/n, p)
		}
	}
	if truncGeom(r, 0, 10) != 0 || truncGeom(r, 0.5, 0) != 0 {
		t.Error("degenerate truncGeom not 0")
	}
}

func TestReads_ErrorModelApplied(t *testing.T) {
	tab := randomTable(t, 4, 10, 500, 800, false)
	c := DefaultConfig(300, 13)
	c.Run.Channels = 512
	s, err := NewSimulator(tab, DefaultTPM(tab, ZipfModel{1, 0.5}, c.Seed), c)
	if err != nil {
		t.Fatal(err)
	}
	diff := 0
	for _, rd := range collect(t, s) {
		if len(rd.Record.Qual) != len(rd.Record.Seq) {
			t.Fatal("quality length mismatch")
		}
		if !bytes.Equal(rd.Record.Seq, rd.Construct) {
			diff++
		}
	}
	// SUP: ~3% errors per base over ~700 bases; every read should differ.
	if diff < 290 {
		t.Errorf("only %d of 300 SUP reads differ from their construct", diff)
	}
}

// --- output ----------------------------------------------------------------

var ontHeader = regexp.MustCompile(`^@([0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}) runid=([0-9a-f]{40}) read=(\d+) ch=(\d+) start_time=(\S+) flow_cell_id=(\S+) protocol_group_id=(\S+) sample_id=(\S+)$`)

func readTSV(t *testing.T, path string) [][]string {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	var rows [][]string
	for _, l := range strings.Split(strings.TrimSuffix(string(b), "\n"), "\n") {
		rows = append(rows, strings.Split(l, "\t"))
	}
	return rows
}

type fqRec struct{ header, seq, qual string }

func readFASTQFile(t *testing.T, path string) []fqRec {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	var r *bufio.Scanner
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			t.Fatal(err)
		}
		r = bufio.NewScanner(gz)
	} else {
		r = bufio.NewScanner(f)
	}
	r.Buffer(make([]byte, 1<<20), 1<<24)
	var out []fqRec
	var lines []string
	for r.Scan() {
		lines = append(lines, r.Text())
	}
	if len(lines)%4 != 0 {
		t.Fatalf("FASTQ has %d lines", len(lines))
	}
	for i := 0; i < len(lines); i += 4 {
		if lines[i+2] != "+" {
			t.Fatalf("bad FASTQ separator %q", lines[i+2])
		}
		out = append(out, fqRec{lines[i], lines[i+1], lines[i+3]})
	}
	return out
}

func TestWriteRun_FilesConsistent(t *testing.T) {
	tab, _ := fixtureTable(t)
	const N = 500
	c := DefaultConfig(N, 17)
	c.NovelRate = 0.4
	c.Run.Channels = 512
	s, err := NewSimulator(tab, []float64{1, 2, 3, 4}, c)
	if err != nil {
		t.Fatal(err)
	}
	dir := t.TempDir()
	res, err := WriteRun(dir, s, OutputOptions{Gzip: true})
	if err != nil {
		t.Fatal(err)
	}
	if res.FASTQ != filepath.Join(dir, "reads.fastq.gz") || res.Reads != N {
		t.Errorf("result %+v", res)
	}
	fq := readFASTQFile(t, res.FASTQ)
	truth := readTSV(t, filepath.Join(dir, TruthName))
	summary := readTSV(t, filepath.Join(dir, SummaryName))
	counts := readTSV(t, filepath.Join(dir, CountsName))
	if len(fq) != N || len(truth) != N+1 || len(summary) != N+1 {
		t.Fatalf("%d FASTQ records, %d truth rows, %d summary rows; want %d each (+ headers)", len(fq), len(truth)-1, len(summary)-1, N)
	}
	if !reflect.DeepEqual(truth[0], TruthColumns) || !reflect.DeepEqual(counts[0], CountsColumns) {
		t.Errorf("headers %v / %v", truth[0], counts[0])
	}
	info := s.RunInfo()
	ids := map[string]bool{}
	for i, r := range fq {
		m := ontHeader.FindStringSubmatch(r.header)
		if m == nil {
			t.Fatalf("bad ONT header %q", r.header)
		}
		if m[2] != info.RunID || m[6] != info.FlowCellID {
			t.Errorf("header run/flow cell %s %s, want %s %s", m[2], m[6], info.RunID, info.FlowCellID)
		}
		if len(r.qual) != len(r.seq) {
			t.Fatalf("read %d: quality length %d != %d", i, len(r.qual), len(r.seq))
		}
		for _, q := range r.qual {
			if q < '!'+1 || q > '!'+50 {
				t.Fatalf("quality char %q outside Q1..Q50", q)
			}
		}
		if ids[m[1]] {
			t.Fatalf("duplicate read ID %s", m[1])
		}
		ids[m[1]] = true
		tr, sm := truth[i+1], summary[i+1]
		if tr[0] != m[1] || sm[1] != m[1] {
			t.Fatalf("row %d: truth %s, summary %s, FASTQ %s", i, tr[0], sm[1], m[1])
		}
		if tr[12] != strconv.Itoa(len(r.seq)) || sm[6] != tr[12] || sm[0] != "reads.fastq.gz" || sm[3] != m[4] {
			t.Fatalf("row %d: lengths truth %s summary %s FASTQ %d; summary file %s ch %s vs %s", i, tr[12], sm[6], len(r.seq), sm[0], sm[3], m[4])
		}
		start, _ := strconv.Atoi(tr[9])
		end, _ := strconv.Atoi(tr[10])
		isoLen, _ := strconv.Atoi(tr[11])
		if end != isoLen || start < 0 || start >= end || (tr[6] == "true") != (start == 0) || (tr[4] == "") != (tr[3] == tr[1]) {
			t.Fatalf("inconsistent truth row %v", tr)
		}
	}
	// counts.tsv: reads sum to N and match truth; TPM sums to 10^6.
	perIso := map[string]int{}
	for _, tr := range truth[1:] {
		perIso[tr[3]]++
	}
	total, tpm := 0, 0.0
	for _, row := range counts[1:] {
		n, _ := strconv.Atoi(row[7])
		v, _ := strconv.ParseFloat(row[5], 64)
		total += n
		tpm += v
		if perIso[row[0]] != n {
			t.Errorf("counts.tsv %s = %d, truth has %d", row[0], n, perIso[row[0]])
		}
	}
	if total != N || math.Abs(tpm-1e6) > 1e-3 {
		t.Errorf("counts sum %d (want %d), TPM sum %g", total, N, tpm)
	}
	// Novel isoforms: TXB_skip2 and TXD_alt3, with the GTF.
	novel := readTSV(t, filepath.Join(dir, NovelTSVName))
	if len(novel) != 3 || novel[1][0] != "TXB_skip2" || novel[1][8] != "9-42" || novel[2][0] != "TXD_alt3" {
		t.Errorf("novel_isoforms.tsv = %v", novel)
	}
	gtf, err := os.ReadFile(filepath.Join(dir, NovelGTFName))
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(gtf), "chr1\tphantom_science\texon\t3\t8\t.\t-\t.\tgene_id \"G2\"; transcript_id \"TXB_skip2\"") {
		t.Errorf("novel GTF lacks TXB_skip2 exon 3-8:\n%s", gtf)
	}
	// Round trip: the novel GTF parses, and the isoforms rebuild from the
	// genome to the simulated sequences.
	ann, err := ParseGTF(filepath.Join(dir, NovelGTFName))
	if err != nil {
		t.Fatal(err)
	}
	rebuilt, _, err := BuildFromGenome(filepath.Join("testdata", "tiny_genome.fna"), ann)
	if err != nil {
		t.Fatal(err)
	}
	for _, sp := range s.Species() {
		if sp.Novel == NovelNone {
			continue
		}
		if got := rebuilt.Lookup(sp.ID); got == nil || !bytes.Equal(got.Seq, sp.Seq) {
			t.Errorf("%s from novel GTF does not match the simulated isoform", sp.ID)
		}
	}
}

func TestWriteRun_Deterministic(t *testing.T) {
	tab := randomTable(t, 6, 30, 200, 600, true)
	run := func(seed uint64, gz bool) string {
		c := DefaultConfig(300, seed)
		c.Run.Channels = 512
		c.NovelRate = 0.2
		s, err := NewSimulator(tab, DefaultTPM(tab, ZipfModel{1, 0.5}, seed), c)
		if err != nil {
			t.Fatal(err)
		}
		dir := t.TempDir()
		if _, err := WriteRun(dir, s, OutputOptions{Gzip: gz}); err != nil {
			t.Fatal(err)
		}
		return dir
	}
	a, b, other := run(42, true), run(42, true), run(43, true)
	for _, name := range []string{FASTQName + ".gz", TruthName, CountsName, NovelTSVName, NovelGTFName, SummaryName} {
		x, _ := os.ReadFile(filepath.Join(a, name))
		y, _ := os.ReadFile(filepath.Join(b, name))
		if len(x) == 0 || !bytes.Equal(x, y) {
			t.Errorf("%s differs between runs with the same seed (or is empty)", name)
		}
	}
	x, _ := os.ReadFile(filepath.Join(a, FASTQName+".gz"))
	z, _ := os.ReadFile(filepath.Join(other, FASTQName+".gz"))
	if bytes.Equal(x, z) {
		t.Error("different seeds gave identical FASTQ")
	}
	// Plain and gzip output carry the same reads.
	plain := run(42, false)
	p := readFASTQFile(t, filepath.Join(plain, FASTQName))
	g := readFASTQFile(t, filepath.Join(a, FASTQName+".gz"))
	if !reflect.DeepEqual(p, g) {
		t.Error("plain and gzip FASTQ differ")
	}
}

// Two isoforms of a gene that differ only in their middle exon would give
// the same exon-skip isoform; only the first gets it, the second falls back
// to an alternative 3′ end.
func TestDesignNovel_NoDuplicateNovelIsoforms(t *testing.T) {
	a := strings.Repeat("ACGTTGCA", 8) // 64
	b1 := strings.Repeat("GGCATT", 10) // 60
	b2 := strings.Repeat("TTACGG", 10) // 60
	c := strings.Repeat("CATGACCA", 8) // 64
	ex := []Exon{{1, 64}, {101, 160}, {201, 264}}
	tab, err := NewTable([]*Transcript{
		{ID: "T1", GeneID: "G", Seq: []byte(a + b1 + c), Chrom: "c", Strand: '+', Exons: ex},
		{ID: "T2", GeneID: "G", Seq: []byte(a + b2 + c), Chrom: "c", Strand: '+', Exons: ex},
	})
	if err != nil {
		t.Fatal(err)
	}
	nov := designNovel(rand.New(rand.NewPCG(1, 1)), tab)
	if nov[0] == nil || nov[0].Novel != NovelSkipExon || string(nov[0].Seq) != a+c {
		t.Fatalf("T1 novel = %+v, want skip of exon 2", nov[0])
	}
	if nov[1] == nil || nov[1].Novel != NovelAlt3Prime {
		t.Fatalf("T2 novel = %+v, want alt_3prime (its skip duplicates T1's)", nov[1])
	}
}
