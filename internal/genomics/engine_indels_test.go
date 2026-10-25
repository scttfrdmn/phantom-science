//go:build !release

package genomics

// Tests of Illumina indel errors (#30) in the engine: reads keep their
// length and stay their fragment's ends apart from errors, the errors
// settings turn them on and off, truth tables and variant offsets follow
// them, and output stays the same for any number of workers.

import (
	"bytes"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/illumina"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// highIndels is an indel model with rates high enough that most reads have
// an error (background ones anywhere, slips in repeats).
var highIndels = illumina.IndelParams{
	HPIns:      []float64{0.004, 0.006, 0.01, 0.02, 0.05},
	HPDel:      []float64{0.005, 0.008, 0.012, 0.03, 0.06},
	STRIns:     []float64{0.01, 0.02, 0.05},
	STRDel:     []float64{0.01, 0.02, 0.05},
	BgIns:      0.002,
	BgDel:      0.002,
	HPInsSize:  [2][]float64{{0.8, 0.2}, {0.6, 0.4}},
	HPDelSize:  [2][]float64{{0.9, 0.1}, {0.6, 0.4}},
	STRInsSize: []float64{1},
	STRDelSize: []float64{1},
	BgInsSize:  []float64{0.6, 0.3, 0.1},
	BgDelSize:  []float64{0, 0.7, 0.3},
	InsQual:    [4]float64{0, 0.3, 0.1, 0.6},
	DelQual:    [4]float64{0, 0.5, 0.2, 0.3},
}

// prefixDistance is the edit distance between read and the best-matching
// prefix of template (both must start together; the template may run on).
func prefixDistance(read, template string) int {
	n, m := len(read), min(len(template), len(read)+40)
	prev, cur := make([]int, m+1), make([]int, m+1)
	for j := range prev {
		prev[j] = j
	}
	for i := 1; i <= n; i++ {
		cur[0] = i
		for j := 1; j <= m; j++ {
			d := prev[j-1]
			if read[i-1] != template[j-1] {
				d++
			}
			cur[j] = min(d, prev[j]+1, cur[j-1]+1)
		}
		prev, cur = cur, prev
	}
	return minOf(prev)
}

func minOf(x []int) int {
	m := x[0]
	for _, v := range x {
		m = min(m, v)
	}
	return m
}

// globalDistance is the edit distance between a and b.
func globalDistance(a, b string) int {
	prev, cur := make([]int, len(b)+1), make([]int, len(b)+1)
	for j := range prev {
		prev[j] = j
	}
	for i := 1; i <= len(a); i++ {
		cur[0] = i
		for j := 1; j <= len(b); j++ {
			d := prev[j-1]
			if a[i-1] != b[j-1] {
				d++
			}
			cur[j] = min(d, prev[j]+1, cur[j-1]+1)
		}
		prev, cur = cur, prev
	}
	return prev[len(b)]
}

// hamming counts the positions where a and b (of equal length) differ.
func hamming(a, b string) int {
	n := 0
	for i := range a {
		if a[i] != b[i] {
			n++
		}
	}
	return n
}

// repeatBases returns n random bases with a long homopolymer or tandem
// repeat every 60 bases, where indel errors concentrate.
func repeatBases(r *rand.Rand, n int) []byte {
	b := make([]byte, 0, n+40)
	for len(b) < n {
		for range 40 {
			b = append(b, "ACGT"[r.IntN(4)])
		}
		switch r.IntN(3) {
		case 0:
			b = append(b, bytes.Repeat([]byte{"ACGT"[r.IntN(4)]}, 15+r.IntN(10))...)
		case 1:
			b = append(b, bytes.Repeat([]byte("CA"), 8+r.IntN(6))...)
		default:
			b = append(b, bytes.Repeat([]byte("AGAT"), 5+r.IntN(4))...)
		}
	}
	return b[:n]
}

// repeatReference writes a two-contig reference of repeatBases.
func repeatReference(t *testing.T) testReference {
	t.Helper()
	r := rand.New(rand.NewPCG(30, 30)) // #nosec G404 -- test data
	ref := testReference{path: filepath.Join(t.TempDir(), "repeats.fa"), seqs: map[string]string{
		"chr1": string(repeatBases(r, 30000)), "chr2": string(repeatBases(r, 20000))}}
	var buf bytes.Buffer
	for _, name := range []string{"chr1", "chr2"} {
		seq := ref.seqs[name]
		buf.WriteString(">" + name + "\n")
		for i := 0; i < len(seq); i += 60 {
			buf.WriteString(seq[i:min(i+60, len(seq))] + "\n")
		}
	}
	if err := os.WriteFile(ref.path, buf.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
	return ref
}

// applyIndels rebuilds a read of n bases from its template and errors
// (inserted bases taken from the read).
func applyIndels(tmpl, read []byte, edits []illumina.Indel, n int) []byte {
	var out []byte
	i, e := 0, 0
	for len(out) < n {
		for e < len(edits) && edits[e].Pos == i && edits[e].Cycle == len(out) {
			ed := edits[e]
			e++
			if ed.Insertion {
				out = append(out, read[ed.Cycle:ed.Cycle+ed.Len]...)
			} else {
				i += ed.Len
			}
		}
		if len(out) >= n || i >= len(tmpl) {
			break
		}
		out = append(out, tmpl[i])
		i++
	}
	return out
}

// Each mate is its end of the fragment with exactly its indel errors: R1
// reads the fragment from its start, R2 its reverse complement, both on
// into the fragment after deletions, keeping the read length; without
// errors they are frag[:L] and revcomp(frag[F−L:]) as before.
func TestIndels_MateIsFragmentEndWithErrors(t *testing.T) {
	c := binnedConfig(t.TempDir())
	e := NewEngine(c)
	defer e.Close()
	e.indelParams = &highIndels
	p, err := e.newPipeline(GenerateOpts{}, nil)
	if err != nil {
		t.Fatal(err)
	}
	if p.indels == nil {
		t.Fatal("binned qualities with errors quality: no indel model")
	}
	p.applyErrors = nil                 // substitutions apart: the read is the template with its indels
	src := rand.New(rand.NewPCG(30, 4)) // #nosec G404 -- test data
	L := c.ReadLength
	withErrors := 0
	for i := range 3000 {
		frag := repeatBases(src, L+i%200)
		for _, mate := range []int{1, 2} {
			r := rand.New(rand.NewPCG(uint64(i), uint64(mate))) // #nosec G404 -- test stream
			seq, qual, edits := p.mate(r, frag, L, mate, -1)
			if len(seq) != L || len(qual) != L {
				t.Fatalf("fragment %d mate %d: %d bases, %d qualities, want %d", i, mate, len(seq), len(qual), L)
			}
			tmpl := frag
			if mate == 2 {
				tmpl = seqread.ReverseComplement(frag)
			}
			if got := applyIndels(tmpl, seq, edits, L); !bytes.Equal(got, seq) {
				t.Fatalf("fragment %d mate %d: read is not its fragment end with errors %+v", i, mate, edits)
			}
			if used := illumina.TemplateUsed(L, edits); used > len(frag) {
				t.Fatalf("fragment %d mate %d: read %d bases of a %d-base fragment", i, mate, used, len(frag))
			}
			if len(edits) == 0 {
				want := frag[:L]
				if mate == 2 {
					want = seqread.ReverseComplement(frag[len(frag)-L:])
				}
				if !bytes.Equal(seq, want) {
					t.Fatalf("fragment %d mate %d: error-free read is not the fragment end", i, mate)
				}
				continue
			}
			withErrors++
		}
	}
	if withErrors < 1000 {
		t.Fatalf("%d of 6000 reads with indel errors", withErrors)
	}
}

// In engine output (errors quality: substitutions and indels) every read
// keeps the read length and is its end of the fragment apart from a few
// errors; some reads are shifted against the fragment (an indel), which
// errors "substitutions" never does; errors "none" gives the fragment ends
// exactly.
func TestIndels_ErrorsSettings(t *testing.T) {
	run := func(t *testing.T, errs string) ([]fastqRec, [][]byte) {
		c := binnedConfig(t.TempDir())
		c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 100, 3000, 1
		c.Illumina.Errors = errs
		c.Illumina.InsertSize = types.InsertSizeSpec{Mean: 300, SD: 40}
		e := NewEngine(c)
		defer e.Close()
		e.indelParams = &highIndels
		s := &repeatSampler{r: rand.New(rand.NewPCG(1, 2))} // #nosec G404 -- test data
		e.fragments = s.sample
		var buf bytes.Buffer
		if err := e.GenerateToWriter(&buf, GenerateOpts{}); err != nil {
			t.Fatal(err)
		}
		return parseFASTQ(t, buf.String()), s.frags
	}
	const L = 100
	for _, errs := range []string{types.IlluminaErrorsQuality, types.IlluminaErrorsSubstitutions, types.IlluminaErrorsNone} {
		t.Run(errs, func(t *testing.T) {
			recs, frags := run(t, errs)
			shifted := 0
			for i, f := range frags {
				r1, r2 := recs[2*i].seq, recs[2*i+1].seq
				if len(r1) != L || len(r2) != L || len(recs[2*i].qual) != L {
					t.Fatalf("pair %d: reads of %d and %d bases, want %d", i, len(r1), len(r2), L)
				}
				t1, t2 := string(f), string(seqread.ReverseComplement(f))
				if errs == types.IlluminaErrorsNone {
					if r1 != t1[:L] || r2 != t2[:L] {
						t.Fatalf("pair %d: errors none, reads are not the fragment ends", i)
					}
					continue
				}
				for _, m := range [2][2]string{{r1, t1}, {r2, t2}} {
					d, h := prefixDistance(m[0], m[1]), hamming(m[0], m[1][:L])
					if d > 25 {
						t.Fatalf("pair %d: read %d edits from its fragment end", i, d)
					}
					if h > d+8 {
						shifted++
					}
				}
			}
			switch {
			case errs == types.IlluminaErrorsQuality && shifted < 500:
				t.Errorf("errors quality: %d of %d reads shifted by an indel", shifted, 2*len(frags))
			case errs == types.IlluminaErrorsSubstitutions && shifted > 0:
				t.Errorf("errors substitutions: %d reads shifted (indels made)", shifted)
			}
			t.Logf("%s: %d of %d reads shifted", errs, shifted, 2*len(frags))
		})
	}
}

// repeatSampler is a FragmentSampler of repeatBases that records the
// fragments it returns.
type repeatSampler struct {
	r     *rand.Rand
	frags [][]byte
}

func (s *repeatSampler) sample(_ *rand.Rand, n int) (string, int64, bool, []byte, error) {
	f := repeatBases(s.r, n)
	s.frags = append(s.frags, f)
	return "chrT", int64(len(s.frags)) * 1000, false, bytes.Clone(f), nil
}

// Indel errors are made only by the binned model with errors quality: not
// by the linear model (its output is unchanged), nor with errors
// substitutions or none.
func TestIndels_OnlyBinnedQuality(t *testing.T) {
	for _, tc := range []struct {
		model, errs string
		want        bool
	}{
		{"", "", true},
		{types.IlluminaQualityBinned, types.IlluminaErrorsQuality, true},
		{types.IlluminaQualityBinned, types.IlluminaErrorsSubstitutions, false},
		{types.IlluminaQualityBinned, types.IlluminaErrorsNone, false},
		{types.IlluminaQualityLinear, types.IlluminaErrorsQuality, false},
		{types.IlluminaQualityLinear, "", false},
	} {
		c := binnedConfig(t.TempDir())
		c.Illumina.QualityModel, c.Illumina.Errors = tc.model, tc.errs
		e := NewEngine(c)
		m, err := e.indelModel()
		e.Close()
		if err != nil || (m != nil) != tc.want {
			t.Errorf("quality_model %q errors %q: indel model %v (err %v), want %v", tc.model, tc.errs, m != nil, err, tc.want)
		}
	}
	if err := (types.IlluminaSpec{Errors: "indels"}).Validate(); err == nil || !strings.Contains(err.Error(), "substitutions") {
		t.Errorf("errors indels: %v, want an error naming the values", err)
	}
	if err := (types.IlluminaSpec{Errors: types.IlluminaErrorsSubstitutions}).Validate(); err != nil {
		t.Errorf("errors substitutions: %v", err)
	}
}

// Single-end reads placed in a region (no reference) carry their planted
// alleles at the offsets the variants table gives, which follow the indel
// errors before them in the read; some offsets moved.
func TestIndels_SingleEndVariantOffsets(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.PairedEnd = false
	c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 100, 3000, 1
	var vs []VariantRecord
	for pos := int64(1002); pos < 1400; pos += 7 {
		vs = append(vs, VariantRecord{ID: "v", Chromosome: "chr5", Position: pos, Ref: "A", Alt: "G", Freq: 1})
	}
	c.VariantsPath = writeVCF(t, vs)
	c.VariantsRegion = "chr5:1001-1400"
	e := NewEngineLog(c, &bytes.Buffer{})
	defer e.Close()
	e.indelParams = &highIndels
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	r1path := e.outputPath(1, 1)
	reads := parseFASTQ(t, string(mustRead(t, r1path)))
	rows := parsePlanted(t, mustRead(t, variantsPath(r1path)))
	var bad, moved, n int
	for _, a := range rows {
		r := reads[a.index-1]
		if len(r.seq) != c.ReadLength {
			t.Fatalf("read of %d bases", len(r.seq))
		}
		if a.off < 0 || a.off >= len(r.seq) {
			t.Fatalf("allele at offset %d of a %d-base read", a.off, len(r.seq))
		}
		n++
		if r.seq[a.off] != 'G' {
			bad++
		}
		if int64(a.off) != a.pos-1-headerReadStart(t, r.header) {
			moved++
		}
	}
	if n < 10000 || bad > n/50 {
		t.Errorf("%d of %d alleles not at their offset", bad, n)
	}
	if moved < n/20 {
		t.Errorf("only %d of %d offsets moved by indel errors", moved, n)
	}
	t.Logf("%d alleles, %d moved by indels, %d not at their offset (substitutions)", n, moved, bad)
}

// Single-end reads from a reference with indel errors (default model,
// repeat-rich reference): each truth row's span is the reference the read
// covers (minus-strand reads reverse complemented), so the read aligns to
// it end to end with a few errors.
func TestIndels_SingleEndReferenceTruth(t *testing.T) {
	ref := repeatReference(t)
	c := binnedConfig(t.TempDir())
	c.PairedEnd = false
	c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 100, 4000, 1
	c.ReferencePath = ref.path
	e := NewEngineLog(c, &bytes.Buffer{})
	defer e.Close()
	e.indelParams = &highIndels
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	path := e.outputPath(1, 1)
	reads := parseFASTQ(t, string(mustRead(t, path)))
	rows := parseTruth(t, mustRead(t, truthPath(path)))
	changed := 0
	for i, row := range rows {
		read := reads[i].seq
		if len(read) != c.ReadLength {
			t.Fatalf("read of %d bases", len(read))
		}
		span := ref.seqs[row.chrom][row.start : row.start+int64(row.length)]
		fwd := read
		if row.minus {
			fwd = string(seqread.ReverseComplement([]byte(read)))
		}
		if d := globalDistance(fwd, span); d > 20 {
			t.Fatalf("row %d (%s:%d %v, %d bases): read %d edits from its span", i, row.chrom, row.start, row.minus, row.length, d)
		}
		if row.length != c.ReadLength {
			changed++
		}
	}
	if changed < 200 {
		t.Errorf("only %d of %d spans changed by indel errors", changed, len(rows))
	}
}

// With indel errors (default model, repeat-rich reference) paired output is
// the same for 1 and 8 workers, and GenerateFASTQFile writes Generate's
// files.
func TestIndels_WorkersAndGenerateFASTQFileIdentical(t *testing.T) {
	ref := repeatReference(t)
	c := refConfig(t, ref)
	c.Illumina.Errors = ""
	c.ReadsPerBatch, c.TotalBatches = 300, 4
	c.Workers = 1
	one, _ := runTree(t, c)
	c.OutputDir, c.Workers = t.TempDir(), 8
	eight, _ := runTree(t, c)
	if len(one) != len(eight) {
		t.Fatalf("%d files with 1 worker, %d with 8", len(one), len(eight))
	}
	for name, data := range one {
		if !bytes.Equal(data, eight[name]) {
			t.Errorf("%s differs between 1 and 8 workers", name)
		}
	}
	c.Workers, c.OutputDir = 0, t.TempDir()
	files, e := runTree(t, c)
	for mate := 1; mate <= 2; mate++ {
		tui := NewEngine(c)
		path := filepath.Join(c.OutputDir, "tui.fastq")
		if err := tui.GenerateFASTQFile(path, 3, mate); err != nil {
			t.Fatal(err)
		}
		tui.Close()
		if got := mustRead(t, path); !bytes.Equal(got, files[filepath.Base(e.outputPath(3, mate))]) {
			t.Errorf("GenerateFASTQFile(batch 3, mate %d) differs from Generate's", mate)
		}
	}
	// The default model makes indel errors here: some reads do not match
	// their reference span base for base.
	rows := parseTruth(t, files[filepath.Base(truthPath(e.outputPath(1, 1)))])
	r1 := parseFASTQ(t, string(files[filepath.Base(e.outputPath(1, 1))]))
	shifted := 0
	for i, row := range rows {
		fwd := ref.seqs[row.chrom][row.start : row.start+int64(row.length)]
		tmpl := fwd
		if row.minus {
			tmpl = string(seqread.ReverseComplement([]byte(fwd)))
		}
		if hamming(r1[i].seq, tmpl[:c.ReadLength]) > prefixDistance(r1[i].seq, tmpl)+8 {
			shifted++
		}
	}
	if shifted == 0 {
		t.Error("no read 1 shifted by an indel error")
	}
	t.Logf("%d of %d read 1 shifted by indel errors", shifted, len(rows))

	// Regenerate (verify --manifest) reproduces every read, indel errors
	// included.
	m := e.Manifest()
	var want []string
	var runFiles []RunFile
	n := 0
	for _, o := range m.Outputs {
		if o.Role != manifest.RoleFASTQ {
			continue
		}
		for _, r := range parseFASTQ(t, string(files[o.Path])) {
			want = append(want, r.seq)
		}
		runFiles = append(runFiles, RunFile{FlowCell: 1, Batch: n/2 + 1, Mate: n%2 + 1, Reads: int(o.Reads)})
		n++
	}
	cfg := m.Config
	cfg.Seed = m.Seed
	cfg.ReferencePath = ref.path
	var got []string
	if _, err := NewEngineLog(cfg, nil).Regenerate(RegenOptions{Files: runFiles}, func(_ RunFile, seq []byte) error {
		got = append(got, string(seq))
		return nil
	}); err != nil {
		t.Fatal(err)
	}
	if strings.Join(got, "\n") != strings.Join(want, "\n") {
		t.Errorf("Regenerate: %d reads, not the run's %d", len(got), len(want))
	}
}
