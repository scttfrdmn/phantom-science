//go:build !release

package genomics

// Tests of VCF variants planted in file output (#32): allele frequencies,
// shared mate alleles, indels, the variants table, the manifest record,
// workers byte-identity, GenerateFASTQFile parity and regeneration.

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// writeVCF writes a VCF of the given records (CHROM POS ID REF ALT AF).
func writeVCF(t *testing.T, recs []VariantRecord) string {
	t.Helper()
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for _, v := range recs {
		id := v.ID
		if id == "" {
			id = "."
		}
		fmt.Fprintf(&b, "%s\t%d\t%s\t%s\t%s\t.\tPASS\tAF=%g;DP=30\n", v.Chromosome, v.Position, id, v.Ref, v.Alt, v.Freq)
	}
	path := filepath.Join(t.TempDir(), "calls.vcf")
	if err := os.WriteFile(path, []byte(b.String()), 0o600); err != nil {
		t.Fatal(err)
	}
	return path
}

// plantedRow is one line of a variants table.
type plantedRow struct {
	name     string
	index    int
	chrom    string
	pos      int64
	id       string
	ref, alt string
	isAlt    bool
	off      int
}

func parsePlanted(t *testing.T, data []byte) []plantedRow {
	t.Helper()
	sc := bufio.NewScanner(bytes.NewReader(data))
	if !sc.Scan() || sc.Text()+"\n" != variantsHeader {
		t.Fatalf("variants table does not start with the header: %q", sc.Text())
	}
	var rows []plantedRow
	for sc.Scan() {
		f := strings.Split(sc.Text(), "\t")
		if len(f) != 9 || (f[7] != "REF" && f[7] != "ALT") {
			t.Fatalf("bad variants line %q", sc.Text())
		}
		index, err0 := strconv.Atoi(f[1])
		pos, err1 := strconv.ParseInt(f[3], 10, 64)
		off, err2 := strconv.Atoi(f[8])
		if err0 != nil || err1 != nil || err2 != nil {
			t.Fatalf("bad variants line %q", sc.Text())
		}
		rows = append(rows, plantedRow{f[0], index, f[2], pos, f[4], f[5], f[6], f[7] == "ALT", off})
	}
	return rows
}

// siteVariants are the variants of the reference-mode tests, on chr2 of
// the test reference: an AF=1 SNV, an AF=0 SNV, an AF=0.3 SNV, an AF=1
// insertion and an AF=1 deletion.
func siteVariants(ref testReference) []VariantRecord {
	s := ref.seqs["chr2"]
	next := map[byte]byte{'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
	snv := func(id string, pos int64, af float64) VariantRecord {
		return VariantRecord{ID: id, Chromosome: "chr2", Position: pos, Ref: s[pos-1 : pos], Alt: string(next[s[pos-1]]), Freq: af}
	}
	return []VariantRecord{
		snv("snv_af1", 2300, 1),
		snv("snv_af0", 2600, 0),
		snv("snv_af03", 3000, 0.3),
		{ID: "ins", Chromosome: "chr2", Position: 3400, Ref: s[3399:3400], Alt: s[3399:3400] + "TTTTTGG", Freq: 1},
		{ID: "del", Chromosome: "chr2", Position: 3700, Ref: s[3699:3706], Alt: s[3699:3700], Freq: 1},
	}
}

// expectedFragment rebuilds a fragment from its truth row and the alleles
// the variants table records for it: the reference from the fragment start,
// each allele written over its REF (right to left), cut to the fragment
// length, reverse complemented on the minus strand.
func expectedFragment(t *testing.T, ref testReference, row truthRow, alleles []plantedRow) string {
	t.Helper()
	s := ref.seqs[row.chrom]
	end := min(int64(len(s)), row.start+int64(row.length)+64)
	fwd := s[row.start:end]
	sorted := append([]plantedRow(nil), alleles...)
	sort.Slice(sorted, func(i, j int) bool { return sorted[i].pos > sorted[j].pos })
	for _, a := range sorted {
		off := int(a.pos - 1 - row.start)
		if off < 0 || off+len(a.ref) > len(fwd) || fwd[off:off+len(a.ref)] != a.ref {
			t.Fatalf("%s: allele of %s:%d at offset %d does not sit on its REF", row.name, a.chrom, a.pos, off)
		}
		allele := a.ref
		if a.isAlt {
			allele = a.alt
		}
		fwd = fwd[:off] + allele + fwd[off+len(a.ref):]
	}
	fwd = fwd[:row.length]
	if row.minus {
		return string(seqread.ReverseComplement([]byte(fwd)))
	}
	return fwd
}

// checkSiteCounts checks the per-variant allele counts of the variants
// table against the AFs: AF 1 always ALT, AF 0 never, AF 0.3 within 4σ;
// and that every fragment whose span holds a variant's REF has a line for
// it. A fragment that ends within maxSlack bases of the loaded region's end
// (regionEnd, 0-based exclusive) cannot be extended for a deletion and is
// then used without variants (see applyReferenceVariants): it is exempt
// when it has no lines at all.
func checkSiteCounts(t *testing.T, vs []VariantRecord, rows []truthRow, planted map[string][]plantedRow, regionEnd int64) {
	t.Helper()
	maxSlack := 0
	for _, v := range vs {
		maxSlack += max(0, len(v.Ref)-len(v.Alt))
	}
	for _, v := range vs {
		var n, alt int
		for _, row := range rows {
			end := row.start + int64(row.length)
			// Indels planted to the left move the variant in the fragment.
			var shift int64
			for _, a := range planted[row.name] {
				if a.pos < v.Position && a.isAlt {
					shift += int64(len(a.alt) - len(a.ref))
				}
			}
			covers := row.chrom == v.Chromosome && v.Position-1 >= row.start && v.Position-1+shift+int64(len(v.Ref)) <= end
			if covers && len(planted[row.name]) == 0 && end+int64(maxSlack) > regionEnd {
				covers = false
			}
			var got *plantedRow
			for i, a := range planted[row.name] {
				if a.id == v.ID {
					got = &planted[row.name][i]
				}
			}
			if covers && got == nil {
				t.Fatalf("%s (%+v) covers %s but the variants table has no allele for it (has %+v)", row.name, row, v.ID, planted[row.name])
			}
			if got == nil {
				continue
			}
			n++
			if got.isAlt {
				alt++
			}
		}
		switch v.Freq {
		case 1:
			if alt != n {
				t.Errorf("%s (AF 1): %d of %d fragments ALT", v.ID, alt, n)
			}
		case 0:
			if alt != 0 {
				t.Errorf("%s (AF 0): %d of %d fragments ALT", v.ID, alt, n)
			}
		default:
			mean, sd := float64(n)*v.Freq, math.Sqrt(float64(n)*v.Freq*(1-v.Freq))
			if math.Abs(float64(alt)-mean) > 4*sd {
				t.Errorf("%s (AF %g): %d of %d ALT, want %.1f ± 4×%.1f", v.ID, v.Freq, alt, n, mean, sd)
			}
		}
		if n < 100 {
			t.Errorf("%s: only %d fragments carry an allele", v.ID, n)
		}
		t.Logf("%s AF %g: %d of %d ALT", v.ID, v.Freq, alt, n)
	}
}

// File output with a reference and a VCF: every read is the reference
// fragment with the alleles the variants table records (so the table is
// exact, indels included, and both mates carry the fragment's alleles), the
// per-site counts follow the AFs, and the manifest records the VCF by hash
// and the tables with their hashes.
func TestVariants_FileOutputReference(t *testing.T) {
	for _, paired := range []bool{true, false} {
		t.Run(fmt.Sprintf("paired=%v", paired), func(t *testing.T) {
			ref := writeTestReference(t)
			c := refConfig(t, ref)
			c.PairedEnd = paired
			c.ReferenceRegion = "chr2:2001-4000"
			c.ReadsPerBatch, c.TotalBatches = 1500, 2
			// Inserts of 180 ± 20: most pairs' mates overlap, so both read
			// some sites.
			c.Illumina.InsertSize = types.InsertSizeSpec{Mean: 180, SD: 20}
			vs := siteVariants(ref)
			c.VariantsPath = writeVCF(t, vs)
			tree, e := runTree(t, c)
			m := e.Manifest()
			if err := m.VerifyDir(c.OutputDir); err != nil {
				t.Fatal(err)
			}
			checkVariantsManifest(t, m, c, len(e.mates())*c.TotalBatches)

			var rows []truthRow
			planted := map[string][]plantedRow{}
			var bothMates, mateSame int
			for b := 1; b <= c.TotalBatches; b++ {
				r1name := filepath.Base(e.outputPath(b, 1))
				vname := strings.TrimSuffix(r1name, ".fastq") + ".variants.tsv"
				if _, ok := tree[vname]; !ok {
					t.Fatalf("no %s in the output", vname)
				}
				batchRows := parseTruth(t, tree[filepath.Base(truthPath(e.outputPath(b, 1)))])
				for _, a := range parsePlanted(t, tree[vname]) {
					planted[a.name] = append(planted[a.name], a)
				}
				r1 := parseFASTQ(t, string(tree[r1name]))
				var r2 []fastqRec
				if paired {
					r2 = parseFASTQ(t, string(tree[filepath.Base(e.outputPath(b, 2))]))
				}
				for i, row := range batchRows {
					frag := expectedFragment(t, ref, row, planted[row.name])
					L := c.ReadLength
					if r1[i].seq != frag[:L] {
						t.Fatalf("batch %d read %d (%s): R1 is not the fragment with its recorded alleles", b, i, row.name)
					}
					if !paired {
						continue
					}
					if want := string(seqread.ReverseComplement([]byte(frag[len(frag)-L:]))); r2[i].seq != want {
						t.Fatalf("batch %d read %d (%s): R2 is not the fragment with its recorded alleles", b, i, row.name)
					}
					// Both mates read the AF 0.3 site: they carry one allele.
					for _, a := range planted[row.name] {
						if a.id != "snv_af03" {
							continue
						}
						o := a.off
						if row.minus {
							o = len(frag) - 1 - a.off
						}
						if o < L && o >= len(frag)-L {
							bothMates++
							r2fwd := seqread.ReverseComplement([]byte(r2[i].seq))
							if r1[i].seq[o] == r2fwd[o-(len(frag)-L)] {
								mateSame++
							}
						}
					}
				}
				rows = append(rows, batchRows...)
			}
			checkSiteCounts(t, vs, rows, planted, 4000)
			if paired && (bothMates < 20 || mateSame != bothMates) {
				t.Errorf("AF 0.3 site read by both mates in %d pairs, the same base in %d", bothMates, mateSame)
			}
		})
	}
}

// checkVariantsManifest checks the manifest's VCF record (base name, size,
// SHA-256, no path) and the variants tables' outputs.
func checkVariantsManifest(t *testing.T, m *manifest.Manifest, c types.SequencerConfig, readFiles int) {
	t.Helper()
	n, sum, err := manifest.HashFile(c.VariantsPath)
	if err != nil {
		t.Fatal(err)
	}
	if m.Variants == nil || *m.Variants != (manifest.Input{Name: "calls.vcf", Bytes: n, SHA256: sum}) {
		t.Errorf("manifest variants = %+v, want calls.vcf %d bytes sha256 %s", m.Variants, n, sum)
	}
	data, err := m.Marshal()
	if err != nil {
		t.Fatal(err)
	}
	if strings.Contains(string(data), filepath.Dir(c.VariantsPath)) || m.Config.VariantsPath != "" {
		t.Errorf("the manifest records the VCF's path")
	}
	var back manifest.Manifest
	if err := json.Unmarshal(data, &back); err != nil || back.Variants == nil || *back.Variants != *m.Variants {
		t.Errorf("manifest variants do not round-trip: %+v, %v", back.Variants, err)
	}
	tables := 0
	for _, o := range m.Outputs {
		if strings.HasSuffix(o.Path, ".variants.tsv") {
			tables++
			if o.Role != manifest.RoleTruth || o.SHA256 == "" {
				t.Errorf("variants table output %+v", o)
			}
		}
	}
	if want := c.TotalBatches; tables != want {
		t.Errorf("%d variants tables in the manifest, want %d", tables, want)
	}
	if m.Totals.Files != readFiles {
		t.Errorf("manifest counts %d read files, want %d", m.Totals.Files, readFiles)
	}
}

// Without a reference, positioned variants place file-output reads in the
// span of the variants: each allele the table records is in the read at its
// offset, every single-end read whose span holds an AF 1 variant carries
// it, and the AF 0.3 counts follow the AF.
func TestVariants_FileOutputRegionMode(t *testing.T) {
	for _, paired := range []bool{false, true} {
		t.Run(fmt.Sprintf("paired=%v", paired), func(t *testing.T) {
			c := binnedConfig(t.TempDir())
			c.PairedEnd = paired
			c.Illumina.Errors = types.IlluminaErrorsNone
			c.Illumina.InsertSize = types.InsertSizeSpec{Mean: 300, SD: 30}
			c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 100, 2000, 1
			vs := []VariantRecord{
				{ID: "a", Chromosome: "chr5", Position: 1000, Ref: "A", Alt: "AGGGGGG", Freq: 1},
				{ID: "b", Chromosome: "chr5", Position: 1100, Ref: "C", Alt: "T", Freq: 0.3},
				{ID: "c", Chromosome: "chr5", Position: 1200, Ref: "GTTTT", Alt: "G", Freq: 1},
			}
			c.VariantsPath = writeVCF(t, vs)
			placed := "chr5:900-1304" // the variants' span plus one read length
			if paired {
				// Fragments of up to ~400 bases need a wider region.
				c.VariantsRegion = "chr5:501-1700"
				placed = c.VariantsRegion
			}
			var log bytes.Buffer
			e := NewEngineLog(c, &log)
			if err := e.Generate(); err != nil {
				t.Fatal(err)
			}
			if !strings.Contains(log.String(), "Placing reads in "+placed) {
				t.Errorf("placement not logged: %q", log.String())
			}
			r1path := e.outputPath(1, 1)
			r1 := parseFASTQ(t, string(mustRead(t, r1path)))
			rows := parsePlanted(t, mustRead(t, variantsPath(r1path)))
			byIndex := map[int][]plantedRow{}
			for _, a := range rows {
				byIndex[a.index] = append(byIndex[a.index], a)
			}
			var altB, nB int
			for i, r := range r1 {
				start := headerReadStart(t, r.header)
				for _, a := range byIndex[i+1] {
					if a.name != readName(r) {
						t.Fatalf("read %d is %s, the variants table names %s", i+1, readName(r), a.name)
					}
					allele := a.ref
					if a.isAlt {
						allele = a.alt
					}
					if a.off >= len(r.seq) {
						continue // paired: an allele of the fragment that R1 does not read
					}
					allele = allele[:min(len(allele), len(r.seq)-a.off)]
					if r.seq[a.off:a.off+len(allele)] != allele {
						t.Fatalf("%s: %s allele %q not at offset %d: %s", readName(r), a.id, allele, a.off, r.seq)
					}
					if a.id == "b" {
						nB++
						if a.isAlt {
							altB++
						}
					}
				}
				if paired {
					continue
				}
				for _, v := range vs {
					if v.Freq != 1 || v.Position-1 < start || v.Position-1+int64(len(v.Ref)) > start+int64(len(r.seq)) {
						continue
					}
					var found bool
					for _, a := range byIndex[i+1] {
						found = found || (a.id == v.ID && a.isAlt)
					}
					if !found {
						t.Fatalf("read at %d covers %s (AF 1) but does not carry it", start, v.ID)
					}
				}
			}
			mean, sd := float64(nB)*0.3, math.Sqrt(float64(nB)*0.21)
			if nB < 100 || math.Abs(float64(altB)-mean) > 4*sd {
				t.Errorf("AF 0.3 variant: %d of %d ALT", altB, nB)
			}
			t.Logf("AF 0.3 variant: %d of %d ALT", altB, nB)
		})
	}
}

func mustRead(t *testing.T, path string) []byte {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	return data
}

// With variants, output is byte-identical for 1 and 8 workers (FASTQ,
// truth and variants tables, manifest), GenerateFASTQFile reproduces a
// file of the run, and Regenerate reproduces every read given the VCF and
// refuses without it or with another VCF.
func TestVariants_WorkersParityAndRegenerate(t *testing.T) {
	ref := writeTestReference(t)
	c := refConfig(t, ref)
	c.ReferenceRegion = "chr2:2001-4000"
	c.TotalBatches = 0
	c.ReadsPerBatch = 150
	c.TotalBases = 150 * 2 * 100 * 5.5 // cuts batch 6
	c.VariantsPath = writeVCF(t, siteVariants(ref))
	c.Workers = 1
	one, e1 := runTree(t, c)
	c.OutputDir, c.Workers = t.TempDir(), 8
	eight, _ := runTree(t, c)
	if len(one) != len(eight) || len(one) < 1+4*6 {
		t.Fatalf("%d files with 1 worker, %d with 8", len(one), len(eight))
	}
	for name, data := range one {
		if !bytes.Equal(data, eight[name]) {
			t.Errorf("%s differs between 1 and 8 workers", name)
		}
	}

	for mate := 1; mate <= 2; mate++ {
		tui := NewEngine(c)
		path := filepath.Join(c.OutputDir, "tui.fastq")
		if err := tui.GenerateFASTQFile(path, 2, mate); err != nil {
			t.Fatal(err)
		}
		if name := filepath.Base(e1.outputPath(2, mate)); !bytes.Equal(mustRead(t, path), one[name]) {
			t.Errorf("GenerateFASTQFile(batch 2, mate %d) differs from Generate's %s", mate, name)
		}
	}

	m := e1.Manifest()
	var want []string
	var files []RunFile
	n := 0
	for _, o := range m.Outputs {
		if o.Role != manifest.RoleFASTQ {
			continue
		}
		for _, r := range parseFASTQ(t, string(one[o.Path])) {
			want = append(want, r.seq)
		}
		files = append(files, RunFile{FlowCell: 1, Batch: n/2 + 1, Mate: n%2 + 1, Reads: int(o.Reads)})
		n++
	}
	regen := func(cfg types.SequencerConfig) ([]string, error) {
		var got []string
		_, err := NewEngineLog(cfg, nil).Regenerate(RegenOptions{Files: files, Variants: m.Variants, Reference: m.Reference}, func(_ RunFile, seq []byte) error {
			got = append(got, string(seq))
			return nil
		})
		return got, err
	}
	cfg := m.Config
	cfg.Seed = m.Seed
	cfg.ReferencePath = ref.path // the manifest records it by hash (#37)
	if _, err := regen(cfg); err == nil || !strings.Contains(err.Error(), "calls.vcf") {
		t.Errorf("Regenerate without the VCF: err = %v", err)
	}
	other := siteVariants(ref)[:2]
	cfg.VariantsPath = writeVCF(t, other)
	if _, err := regen(cfg); err == nil || !strings.Contains(err.Error(), "not calls.vcf") {
		t.Errorf("Regenerate with another VCF: err = %v", err)
	}
	cfg.VariantsPath = c.VariantsPath
	got, err := regen(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if strings.Join(got, "\n") != strings.Join(want, "\n") {
		t.Errorf("Regenerate with the VCF: %d reads, want the run's %d", len(got), len(want))
	}
}

// strict_ref fails file output before anything is written; without it the
// mismatching variant is skipped with a warning.
func TestVariants_FileStrictRef(t *testing.T) {
	ref := writeTestReference(t)
	c := refConfig(t, ref)
	c.ReferenceRegion = "chr2:2001-4000"
	vs := siteVariants(ref)
	vs[0].Ref = string(map[byte]byte{'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[vs[0].Ref[0]])
	c.VariantsPath = writeVCF(t, vs)
	c.StrictRef = true
	c.OutputDir = filepath.Join(t.TempDir(), "out")
	e := NewEngineLog(c, nil)
	if err := e.Generate(); err == nil || !strings.Contains(err.Error(), "strict_ref") || !strings.Contains(err.Error(), "chr2:2300") {
		t.Fatalf("strict_ref: err = %v", err)
	}
	if _, err := os.Stat(c.OutputDir); !os.IsNotExist(err) {
		t.Errorf("strict_ref created the output directory: %v", err)
	}
	c.StrictRef = false
	var log bytes.Buffer
	e = NewEngineLog(c, &log)
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(log.String(), "1 variants skipped: their REF does not match the reference") {
		t.Errorf("no mismatch warning in %q", log.String())
	}
	for _, a := range parsePlanted(t, mustRead(t, variantsPath(e.outputPath(1, 1)))) {
		if a.id == "snv_af1" {
			t.Fatalf("the mismatching variant was planted: %+v", a)
		}
	}
}

// Configured variants also drive streamed output, and GenerateOpts
// variants cannot be passed as well.
func TestVariants_ConfiguredStream(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Seed = 32
	c.Illumina.Errors = types.IlluminaErrorsNone
	c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 50, 300, 1
	c.VariantsPath = writeVCF(t, []VariantRecord{{Chromosome: "chr3", Position: 1000, Ref: "A", Alt: "AGGGGGGGG", Freq: 1}})
	got := streamOutput(t, c, GenerateOpts{})
	c2 := c
	c2.VariantsPath = ""
	want := streamOutput(t, c2, GenerateOpts{
		Variants:  []VariantRecord{{Chromosome: "chr3", Position: 1000, Ref: "A", Alt: "AGGGGGGGG", Freq: 1}},
		RefRegion: &GenomicRegion{Chromosome: "chr3", Start: 999 - 50, End: 1000 + 50},
	})
	if !bytes.Equal(got, want) {
		t.Error("configured variants stream differs from the same variants passed in GenerateOpts")
	}
	e := NewEngineLog(c, nil)
	err := e.GenerateToWriter(&bytes.Buffer{}, GenerateOpts{Variants: []VariantRecord{{Ref: "A", Alt: "C", Freq: 1}}})
	if err == nil || !strings.Contains(err.Error(), "variants_path") {
		t.Errorf("GenerateOpts variants with variants_path: err = %v", err)
	}
	c.VariantsRegion, c.ReferencePath = "chr3:1-500", "x.fa"
	if err := NewEngineLog(c, nil).Generate(); err == nil {
		t.Error("variants_region with reference_path should fail")
	}
}

// The allele log shifts recorded offsets by later (further left) indels,
// drops overwritten alleles and those trimmed off.
func TestPlantLog(t *testing.T) {
	v1, v2, v3 := &VariantRecord{ID: "1"}, &VariantRecord{ID: "2"}, &VariantRecord{ID: "3"}
	l := &plantLog{}
	l.splice(v1, true, 20, 1, 1) // SNV at 20
	l.splice(v2, true, 10, 1, 4) // insertion of 3 at 10: v1 moves to 23
	l.splice(v3, true, 22, 3, 1) // overwrites v1 (position-blind)
	got := map[string]int{}
	for _, a := range l.alleles {
		got[a.v.ID] = a.off
	}
	if len(got) != 2 || got["2"] != 10 || got["3"] != 22 {
		t.Errorf("alleles %v, want 2@10 3@22", got)
	}
	l.trim(22)
	if len(l.alleles) != 1 || l.alleles[0].v != v2 {
		t.Errorf("after trim to 22: %+v", l.alleles)
	}
	l.splice(v1, true, 20, 1, 5) // an insertion running past the end is kept
	l.trim(22)
	if len(l.alleles) != 2 {
		t.Errorf("after trim of an insertion at 20 to 22: %+v", l.alleles)
	}
	var nilLog *plantLog
	nilLog.splice(v1, true, 0, 1, 1) // no-op
	nilLog.trim(0)
	nilLog.reset()
}

// LoadVariantsVCF reads IDs, positions, the first ALT and AF, from plain
// and gzip VCFs.
func TestLoadVariantsVCF_IDsAndGzip(t *testing.T) {
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\trs1\tA\tG,T\t.\t.\tDP=3;AF=0.25,0.1\n" +
		"chr1\t200\t.\tC\t.\t.\t.\t.\n" + // no ALT: skipped
		"chr1\tx\t.\tC\tT\t.\t.\t.\n" // bad POS: position-blind
	dir := t.TempDir()
	plain := filepath.Join(dir, "v.vcf")
	if err := os.WriteFile(plain, []byte(body), 0o600); err != nil {
		t.Fatal(err)
	}
	gz := filepath.Join(dir, "v.vcf.gz")
	if err := os.WriteFile(gz, gzipBytes(t, []byte(body)), 0o600); err != nil {
		t.Fatal(err)
	}
	for _, p := range []string{plain, gz} {
		vs, err := LoadVariantsVCF(p)
		if err != nil {
			t.Fatal(err)
		}
		want := []VariantRecord{
			{ID: "rs1", Chromosome: "chr1", Position: 100, Ref: "A", Alt: "G", Freq: 0.25},
			{Chromosome: "chr1", Position: 0, Ref: "C", Alt: "T", Freq: 0.5},
		}
		if fmt.Sprint(vs) != fmt.Sprint(want) {
			t.Errorf("%s: %+v, want %+v", filepath.Base(p), vs, want)
		}
	}
}

func gzipBytes(t *testing.T, data []byte) []byte {
	t.Helper()
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	if _, err := zw.Write(data); err != nil {
		t.Fatal(err)
	}
	if err := zw.Close(); err != nil {
		t.Fatal(err)
	}
	return buf.Bytes()
}
