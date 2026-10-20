//go:build !release

package genomics

// Tests of reads sampled from a reference genome (#27): with sequencing
// errors off, every read is the reference at its recorded origin.

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// testReference is a FASTA fixture: chr1 (20 kb, an N run at
// [8000, 8600) and a lower-case stretch) and chr2 (6 kb).
type testReference struct {
	path string
	seqs map[string]string // upper-cased forward bases
}

func writeTestReference(t *testing.T) testReference {
	t.Helper()
	r := rand.New(rand.NewPCG(5, 6))
	gen := func(n int) []byte {
		b := make([]byte, n)
		for i := range b {
			b[i] = "ACGT"[r.IntN(4)]
		}
		return b
	}
	chr1 := gen(20000)
	for i := 8000; i < 8600; i++ {
		chr1[i] = 'N'
	}
	chr2 := gen(6000)
	ref := testReference{path: filepath.Join(t.TempDir(), "ref.fa"), seqs: map[string]string{"chr1": string(chr1), "chr2": string(chr2)}}
	var buf bytes.Buffer
	for _, name := range []string{"chr1", "chr2"} {
		seq := []byte(ref.seqs[name])
		if name == "chr1" {
			seq = append([]byte(nil), seq...)
			copy(seq[100:300], bytes.ToLower(seq[100:300])) // soft-masked
		}
		buf.WriteString(">" + name + " test contig\n")
		for i := 0; i < len(seq); i += 60 {
			buf.Write(seq[i:min(i+60, len(seq))])
			buf.WriteByte('\n')
		}
	}
	if err := os.WriteFile(ref.path, buf.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
	return ref
}

// refConfig is binned paired Illumina output from ref with sequencing
// errors off.
func refConfig(t *testing.T, ref testReference) types.SequencerConfig {
	c := binnedConfig(t.TempDir())
	c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 100, 400, 3
	c.Illumina.Errors = types.IlluminaErrorsNone
	c.ReferencePath = ref.path
	return c
}

// truthRow is one line of a truth table.
type truthRow struct {
	name, chrom string
	start       int64
	minus       bool
	length      int
}

func parseTruth(t *testing.T, data []byte) []truthRow {
	t.Helper()
	sc := bufio.NewScanner(bytes.NewReader(data))
	if !sc.Scan() || sc.Text()+"\n" != truthHeader {
		t.Fatalf("truth table does not start with the header: %q", sc.Text())
	}
	var rows []truthRow
	for sc.Scan() {
		f := strings.Split(sc.Text(), "\t")
		if len(f) != 5 || (f[3] != "+" && f[3] != "-") {
			t.Fatalf("bad truth line %q", sc.Text())
		}
		start, err1 := strconv.ParseInt(f[2], 10, 64)
		n, err2 := strconv.Atoi(f[4])
		if err1 != nil || err2 != nil {
			t.Fatalf("bad truth line %q", sc.Text())
		}
		rows = append(rows, truthRow{f[0], f[1], start, f[3] == "-", n})
	}
	return rows
}

// readName is a FASTQ record's read name.
func readName(r fastqRec) string { return strings.TrimPrefix(strings.Fields(r.header)[0], "@") }

// checkPairs checks every pair against its truth row: R1 and R2 are the two
// ends of the reference fragment, which holds no N.
func checkPairs(t *testing.T, ref testReference, rows []truthRow, r1, r2 []fastqRec, L int) (minus int) {
	t.Helper()
	if len(rows) != len(r1) || len(r1) != len(r2) {
		t.Fatalf("%d truth rows, %d R1, %d R2", len(rows), len(r1), len(r2))
	}
	for i, row := range rows {
		seq, ok := ref.seqs[row.chrom]
		if !ok {
			t.Fatalf("row %d: unknown contig %s", i, row.chrom)
		}
		F := row.length
		if F < L || row.start < 0 || row.start+int64(F) > int64(len(seq)) {
			t.Fatalf("row %d: fragment %s:%d+%d out of range (read length %d)", i, row.chrom, row.start, F, L)
		}
		frag := seq[row.start : row.start+int64(F)]
		if strings.ContainsRune(frag, 'N') {
			t.Fatalf("row %d: fragment %s:%d+%d spans N", i, row.chrom, row.start, F)
		}
		if readName(r1[i]) != row.name || readName(r2[i]) != row.name {
			t.Fatalf("row %d: names %s / %s, truth %s", i, readName(r1[i]), readName(r2[i]), row.name)
		}
		want1, want2 := frag[:L], string(seqread.ReverseComplement([]byte(frag[F-L:])))
		if row.minus {
			minus++
			want1, want2 = string(seqread.ReverseComplement([]byte(frag[F-L:]))), frag[:L]
		}
		if r1[i].seq != want1 {
			t.Fatalf("row %d (%s:%d %v): R1 is not the reference fragment start", i, row.chrom, row.start, row.minus)
		}
		if r2[i].seq != want2 {
			t.Fatalf("row %d (%s:%d %v): R2 is not the reverse complement of the fragment end", i, row.chrom, row.start, row.minus)
		}
	}
	return minus
}

// Paired file output: R1 is the reference fragment's first L bases (as R1
// reads it), R2 the reverse complement of its last L, fragments avoid N,
// both strands and both contigs occur, inserts follow the configured
// distribution, and the truth tables are in the manifest with their hashes.
func TestReference_PairedReadsAreReferenceFragments(t *testing.T) {
	ref := writeTestReference(t)
	c := refConfig(t, ref)
	tree, e := runTree(t, c)
	m := e.Manifest()
	if err := m.VerifyDir(c.OutputDir); err != nil {
		t.Fatal(err)
	}
	var truths []manifest.Output
	for _, o := range m.Outputs {
		if o.Role == manifest.RoleTruth {
			truths = append(truths, o)
		}
	}
	if len(truths) != c.TotalBatches {
		t.Fatalf("%d truth outputs in the manifest, want %d", len(truths), c.TotalBatches)
	}
	if m.Totals.Files != 2*c.TotalBatches {
		t.Errorf("manifest totals count %d read files, want %d (truth tables are not read files)", m.Totals.Files, 2*c.TotalBatches)
	}
	var minus, pairs int
	contigs := map[string]int{}
	var sum, sq float64
	for b := 1; b <= c.TotalBatches; b++ {
		r1name, r2name := filepath.Base(e.outputPath(b, 1)), filepath.Base(e.outputPath(b, 2))
		tname := filepath.Base(truthPath(e.outputPath(b, 1)))
		if tname != strings.TrimSuffix(r1name, ".fastq")+".truth.tsv" {
			t.Fatalf("truth table %s for %s", tname, r1name)
		}
		var found bool
		for _, o := range truths {
			if o.Path == tname {
				found = true
				if o.Reads != int64(c.ReadsPerBatch) {
					t.Errorf("%s: manifest reads %d, want %d", tname, o.Reads, c.ReadsPerBatch)
				}
			}
		}
		if !found {
			t.Fatalf("%s not in the manifest", tname)
		}
		rows := parseTruth(t, tree[tname])
		minus += checkPairs(t, ref, rows, parseFASTQ(t, string(tree[r1name])), parseFASTQ(t, string(tree[r2name])), c.ReadLength)
		for _, row := range rows {
			contigs[row.chrom]++
			sum += float64(row.length)
			sq += float64(row.length) * float64(row.length)
			pairs++
		}
	}
	n := float64(pairs)
	if p := float64(minus) / n; math.Abs(p-0.5) > 3*math.Sqrt(0.25/n) {
		t.Errorf("minus-strand fraction %.3f over %d pairs, want 0.5", p, pairs)
	}
	// chr1 has 19400 A/C/G/T bases, chr2 6000: chr2 gets ~23.6%.
	if p, want := float64(contigs["chr2"])/n, 6000.0/25400; math.Abs(p-want) > 3*math.Sqrt(want*(1-want)/n) {
		t.Errorf("chr2 fraction %.3f, want %.3f", p, want)
	}
	mean := sum / n
	sd := math.Sqrt(sq/n - mean*mean)
	if math.Abs(mean-350) > 3*50/math.Sqrt(n) {
		t.Errorf("insert mean %.1f, want 350", mean)
	}
	t.Logf("%d pairs: %d minus strand, contigs %v, insert %.1f ± %.1f", pairs, minus, contigs, mean, sd)
}

// --region restricts every fragment to the region.
func TestReference_RegionRestrictsOrigins(t *testing.T) {
	ref := writeTestReference(t)
	c := refConfig(t, ref)
	c.ReferenceRegion = "chr2:1,001-4000"
	tree, e := runTree(t, c)
	for b := 1; b <= c.TotalBatches; b++ {
		rows := parseTruth(t, tree[filepath.Base(truthPath(e.outputPath(b, 1)))])
		checkPairs(t, ref, rows, parseFASTQ(t, string(tree[filepath.Base(e.outputPath(b, 1))])),
			parseFASTQ(t, string(tree[filepath.Base(e.outputPath(b, 2))])), c.ReadLength)
		for _, row := range rows {
			if row.chrom != "chr2" || row.start < 1000 || row.start+int64(row.length) > 4000 {
				t.Fatalf("fragment %s:%d+%d outside chr2:1001-4000", row.chrom, row.start, row.length)
			}
		}
	}
	if got, err := e.ReferenceBases("2", 1000, 1010); err != nil || string(got) != ref.seqs["chr2"][1000:1010] {
		t.Errorf("ReferenceBases(2, 1000, 1010) = %q, %v", got, err)
	}
	if _, err := e.ReferenceBases("chr2", 0, 10); err == nil {
		t.Error("ReferenceBases outside the region succeeded")
	}
}

// Single-end reads (binned and linear qualities) are reference substrings
// of the read length on the recorded strand.
func TestReference_SingleEnd(t *testing.T) {
	ref := writeTestReference(t)
	for _, model := range []string{types.IlluminaQualityBinned, types.IlluminaQualityLinear} {
		t.Run(model, func(t *testing.T) {
			c := refConfig(t, ref)
			c.PairedEnd = false
			c.Illumina.QualityModel = model
			tree, e := runTree(t, c)
			for b := 1; b <= c.TotalBatches; b++ {
				recs := parseFASTQ(t, string(tree[filepath.Base(e.outputPath(b, 1))]))
				rows := parseTruth(t, tree[filepath.Base(truthPath(e.outputPath(b, 1)))])
				if len(rows) != len(recs) {
					t.Fatalf("%d truth rows for %d reads", len(rows), len(recs))
				}
				for i, row := range rows {
					want := ref.seqs[row.chrom][row.start : row.start+int64(c.ReadLength)]
					if row.minus {
						want = string(seqread.ReverseComplement([]byte(want)))
					}
					if row.length != c.ReadLength || recs[i].seq != want || readName(recs[i]) != row.name {
						t.Fatalf("read %d is not the reference at %s:%d (minus %v)", i, row.chrom, row.start, row.minus)
					}
				}
			}
		})
	}
}

// Output does not depend on the number of workers (files, truth tables and
// the manifest), including when a bases target cuts a batch; GenerateFASTQFile
// writes the same file as Generate.
func TestReference_WorkersAndGenerateFASTQFileIdentical(t *testing.T) {
	ref := writeTestReference(t)
	c := refConfig(t, ref)
	c.TotalBatches = 0
	c.ReadsPerBatch = 150
	c.TotalBases = 150 * 2 * 100 * 5.5 // cuts batch 6
	c.Workers = 1
	one, _ := runTree(t, c)
	c.OutputDir, c.Workers = t.TempDir(), 8
	eight, _ := runTree(t, c)
	if len(one) != len(eight) || len(one) < 1+3*6 {
		t.Fatalf("%d files with 1 worker, %d with 8", len(one), len(eight))
	}
	for name, data := range one {
		if !bytes.Equal(data, eight[name]) {
			t.Errorf("%s differs between 1 and 8 workers", name)
		}
	}

	c.TotalBases, c.TotalBatches, c.Workers = 0, 2, 0
	c.OutputDir = t.TempDir()
	files, e := runTree(t, c)
	for mate := 1; mate <= 2; mate++ {
		tui := NewEngine(c)
		path := filepath.Join(c.OutputDir, "tui.fastq")
		if err := tui.GenerateFASTQFile(path, 2, mate); err != nil {
			t.Fatal(err)
		}
		tui.Close()
		got, err := os.ReadFile(path)
		if err != nil {
			t.Fatal(err)
		}
		if name := filepath.Base(e.outputPath(2, mate)); !bytes.Equal(got, files[name]) {
			t.Errorf("GenerateFASTQFile(batch 2, mate %d) differs from Generate's %s", mate, name)
		}
	}
}

// Streamed (stdout) output samples from the reference as file output does:
// interleaved mates are fragment ends, both reference substrings, R1 and
// revcomp(R2) on one contig within the insert size; and with gzip.
func TestReference_StreamOutput(t *testing.T) {
	ref := writeTestReference(t)
	c := refConfig(t, ref)
	c.TotalBatches = 2
	recs := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{})))
	if len(recs) != 2*c.TotalBatches*c.ReadsPerBatch {
		t.Fatalf("%d records", len(recs))
	}
	for i := 0; i < len(recs); i += 2 {
		r1, r2 := recs[i].seq, recs[i+1].seq
		rc2 := string(seqread.ReverseComplement([]byte(r2)))
		found := false
		for _, s := range ref.seqs {
			for _, pair := range [][2]string{{r1, rc2}, {string(seqread.ReverseComplement([]byte(rc2))), string(seqread.ReverseComplement([]byte(r1)))}} {
				a := strings.Index(s, pair[0])
				b := strings.Index(s, pair[1])
				if a >= 0 && b >= a && b+len(pair[1])-a <= 800 {
					found = true
				}
			}
		}
		if !found {
			t.Fatalf("pair %d: R1 and revcomp(R2) are not the ends of one reference fragment", i/2)
		}
		if readName(recs[i]) != readName(recs[i+1]) {
			t.Fatalf("pair %d: mate names differ", i/2)
		}
	}

	c.Output.Compression = "gzip"
	gz := streamOutput(t, c, GenerateOpts{})
	zr, err := gzip.NewReader(bytes.NewReader(gz))
	if err != nil {
		t.Fatal(err)
	}
	var plain bytes.Buffer
	if _, err := plain.ReadFrom(zr); err != nil {
		t.Fatal(err)
	}
	c.Output.Compression = ""
	if !bytes.Equal(plain.Bytes(), streamOutput(t, c, GenerateOpts{})) {
		t.Error("gzip stream differs from the plain stream")
	}
}

// ONT reads (lognormal lengths, no basecaller tier: no errors) are
// reference substrings, recorded in the truth table. From a region shorter
// than some sampled lengths, those lengths are redrawn, not an error.
func TestReference_ONT(t *testing.T) {
	ref := writeTestReference(t)
	for name, tc := range map[string]struct {
		region  string
		maxLen  int64 // longest N-free stretch; 0: lengths are not refitted
		longest int
	}{
		"whole": {"", 11400, 5000},
		// Lognormal N50 2000 / median 1500 without a max draws lengths
		// above 2500 often; they are redrawn.
		"small region": {"chr2:3001-5500", 2500, 2500},
	} {
		t.Run(name, func(t *testing.T) {
			c := ontConfig(t.TempDir())
			c.ReadLength = 1500
			c.ReadLengthDistribution = types.ReadLengthSpec{Dist: types.ReadLengthLognormal, N50: 2000, Median: 1500}
			if tc.region == "" {
				c.ReadLengthDistribution.Max = 5000
			}
			c.ReferencePath, c.ReferenceRegion = ref.path, tc.region
			tree, e := runTree(t, c)
			if e.refMaxRead != tc.maxLen {
				t.Errorf("longest N-free stretch %d, want %d", e.refMaxRead, tc.maxLen)
			}
			var n, longest int
			for b := 1; b <= c.TotalBatches; b++ {
				recs := parseFASTQ(t, string(tree[filepath.Base(e.outputPath(b, 1))]))
				rows := parseTruth(t, tree[filepath.Base(truthPath(e.outputPath(b, 1)))])
				if len(rows) != len(recs) {
					t.Fatalf("%d truth rows for %d reads", len(rows), len(recs))
				}
				for i, row := range rows {
					want := ref.seqs[row.chrom][row.start : row.start+int64(row.length)]
					if row.minus {
						want = string(seqread.ReverseComplement([]byte(want)))
					}
					if recs[i].seq != want || readName(recs[i]) != row.name || !ontHeader.MatchString(recs[i].header) {
						t.Fatalf("ONT read %d is not the reference at %s:%d+%d", i, row.chrom, row.start, row.length)
					}
					longest = max(longest, row.length)
					n++
				}
			}
			if longest > tc.longest {
				t.Errorf("read of %d bases, want at most %d", longest, tc.longest)
			}
			t.Logf("%d ONT reads checked, longest %d", n, longest)
		})
	}
}

// Positioned variants are applied in reference coordinates: a variant whose
// REF matches is applied (ALT at frequency 1) to every fragment covering it,
// on both strands; a mismatching REF is skipped and counted in one warning,
// or with strict_ref fails before anything is written.
func TestReference_VariantsAndRefMismatch(t *testing.T) {
	ref := writeTestReference(t)
	c := refConfig(t, ref)
	c.ReferenceRegion = "chr2:2001-3000"
	c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 100, 2000, 1
	c.Illumina.InsertSize = types.InsertSizeSpec{Mean: 200, SD: 20}
	next := map[byte]byte{'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
	s := ref.seqs["chr2"]
	const pos, pos2 = 2500, 2700 // 1-based
	alt := next[s[pos-1]]
	bad := next[s[pos2-1]] // a REF that does not match the reference
	alt2 := next[bad]
	variants := []VariantRecord{
		{Chromosome: "chr2", Position: pos, Ref: s[pos-1 : pos], Alt: string(alt), Freq: 1},
		{Chromosome: "2", Position: pos2, Ref: string(bad), Alt: string(alt2), Freq: 1},
		{Ref: "ACGT", Alt: "TTTT", Freq: 1}, // position-blind: not used
	}
	var log bytes.Buffer
	e := NewEngineLog(c, &log)
	var buf bytes.Buffer
	if err := e.GenerateToWriter(&buf, GenerateOpts{Variants: variants}); err != nil {
		t.Fatal(err)
	}
	e.Close()
	if !strings.Contains(log.String(), "1 variants skipped: their REF does not match the reference") {
		t.Errorf("no mismatch warning in %q", log.String())
	}
	if !strings.Contains(log.String(), "1 variants without CHROM/POS are ignored") {
		t.Errorf("no position-blind warning in %q", log.String())
	}
	recs := parseFASTQ(t, buf.String())
	mut := s[:pos-1] + string(alt) + s[pos:]
	var covered, applied2 int
	for i := 0; i < len(recs); i += 2 {
		r1, r2 := recs[i].seq, recs[i+1].seq
		for _, read := range []string{r1, r2} {
			rc := string(seqread.ReverseComplement([]byte(read)))
			fwd := read
			if !strings.Contains(mut, fwd) && !strings.Contains(s, fwd) {
				fwd = rc
			}
			a := strings.Index(mut, fwd)
			if a < 0 {
				t.Fatalf("read %q is not from the (mutated) reference", read)
			}
			if a <= pos-1 && pos-1 < a+len(fwd) {
				covered++
				if strings.Contains(s, fwd) {
					t.Fatalf("read covering %d carries the REF allele", pos)
				}
			}
			if a <= pos2-1 && pos2-1 < a+len(fwd) && fwd[pos2-1-a] == alt2 {
				applied2++
			}
		}
	}
	if covered < 100 {
		t.Errorf("only %d reads cover the variant", covered)
	}
	if applied2 != 0 {
		t.Errorf("the mismatching variant was applied to %d reads", applied2)
	}
	t.Logf("%d reads carry the ALT at %d", covered, pos)

	c.StrictRef = true
	e = NewEngineLog(c, &log)
	defer e.Close()
	buf.Reset()
	err := e.GenerateToWriter(&buf, GenerateOpts{Variants: variants})
	if err == nil || !strings.Contains(err.Error(), "strict_ref") || !strings.Contains(err.Error(), "2:2700") {
		t.Fatalf("strict_ref: err = %v", err)
	}
	if buf.Len() != 0 {
		t.Errorf("strict_ref wrote %d bytes before failing", buf.Len())
	}
}

// A deletion near the fragment end keeps reads at the read length: the
// fragment is extended from the reference by the deleted length.
func TestReference_DeletionKeepsLength(t *testing.T) {
	ref := writeTestReference(t)
	c := refConfig(t, ref)
	c.PairedEnd = false
	c.ReferenceRegion = "chr2:2001-2400"
	c.ReadsPerBatch, c.TotalBatches = 2000, 1
	s := ref.seqs["chr2"]
	const pos = 2200
	v := VariantRecord{Chromosome: "chr2", Position: pos, Ref: s[pos-1 : pos+9], Alt: s[pos-1 : pos], Freq: 1} // 9-base deletion
	recs := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{Variants: []VariantRecord{v}})))
	mut := s[:pos] + s[pos+9:]
	var deleted int
	for _, r := range recs {
		if len(r.seq) != c.ReadLength {
			t.Fatalf("read of %d bases", len(r.seq))
		}
		fwd := r.seq
		if !strings.Contains(s, fwd) && !strings.Contains(mut, fwd) {
			fwd = string(seqread.ReverseComplement([]byte(fwd)))
		}
		switch {
		case strings.Contains(s, fwd):
		case strings.Contains(mut, fwd):
			deleted++
		default:
			t.Fatalf("read is neither reference nor reference with the deletion")
		}
	}
	if deleted < 100 {
		t.Errorf("only %d reads carry the deletion", deleted)
	}
}

// Configuration errors surface from Generate and GenerateToWriter before
// anything is written.
func TestReference_ConfigErrors(t *testing.T) {
	ref := writeTestReference(t)
	gz := filepath.Join(t.TempDir(), "ref.fa.gz")
	{
		f, err := os.Create(gz)
		if err != nil {
			t.Fatal(err)
		}
		zw := gzip.NewWriter(f)
		data, _ := os.ReadFile(ref.path)
		zw.Write(data)
		zw.Close()
		f.Close()
	}
	cases := map[string]struct {
		edit func(*types.SequencerConfig)
		want string
	}{
		"gzip":           {func(c *types.SequencerConfig) { c.ReferencePath = gz }, "gzip"},
		"missing":        {func(c *types.SequencerConfig) { c.ReferencePath = filepath.Join(t.TempDir(), "none.fa") }, "open reference"},
		"bad region":     {func(c *types.SequencerConfig) { c.ReferenceRegion = "chr9:1-100" }, "chr9"},
		"region too far": {func(c *types.SequencerConfig) { c.ReferenceRegion = "chr2:1-7000" }, "past the end"},
		"region only": {func(c *types.SequencerConfig) {
			c.ReferencePath = ""
			c.ReferenceRegion = "chr2"
		}, "needs reference_path"},
		"hmm path":  {func(c *types.SequencerConfig) { c.HMMModelPath = "model.hmm.json" }, "cannot be combined with an HMM model"},
		"model dir": {func(c *types.SequencerConfig) { c.ModelDir, c.Population = t.TempDir(), "YRI" }, "cannot be combined with an HMM model"},
		"fragment longer than the region": {func(c *types.SequencerConfig) {
			c.ReferenceRegion = "chr2:1-150"
		}, "longer than every N-free stretch"},
	}
	for name, tc := range cases {
		t.Run(name, func(t *testing.T) {
			c := refConfig(t, ref)
			tc.edit(&c)
			e := NewEngineLog(c, nil)
			defer e.Close()
			var buf bytes.Buffer
			err := e.GenerateToWriter(&buf, GenerateOpts{})
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Fatalf("GenerateToWriter: err = %v, want %q", err, tc.want)
			}
			if buf.Len() != 0 && !strings.Contains(tc.want, "N-free") {
				t.Errorf("%d bytes written before the error", buf.Len())
			}
			err = e.Generate()
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Fatalf("Generate: err = %v, want %q", err, tc.want)
			}
			if tc.want != "longer than every N-free stretch" {
				if entries, _ := os.ReadDir(c.OutputDir); len(entries) != 0 {
					t.Errorf("Generate left %d entries in the output directory", len(entries))
				}
			}
		})
	}
	// ReferenceBases without a reference.
	c := binnedConfig(t.TempDir())
	e := NewEngine(c)
	defer e.Close()
	if _, err := e.ReferenceBases("chr1", 0, 10); err == nil {
		t.Error("ReferenceBases without a reference succeeded")
	}
}

// Without a reference the pipeline is unchanged: region-mode read names
// still encode the start (the golden and region tests cover the bytes).
func TestReference_NoReferenceKeepsRegionNames(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Seed = 3
	recs := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{RefRegion: &GenomicRegion{Chromosome: "chr1", Start: 100000, End: 101000}})))
	for _, r := range recs {
		f := strings.Split(readName(r), ":")
		if f[6] != "1" {
			t.Fatalf("region read %s: y = %s, want 1 (start/1e5)", readName(r), f[6])
		}
	}
}
