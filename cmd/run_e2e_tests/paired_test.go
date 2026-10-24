//go:build !release
// +build !release

package main

import (
	"bytes"
	"encoding/csv"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	phantomtesting "github.com/scttfrdmn/phantom-science/internal/testing"
)

// writePairs generates read pairs with config c and writes them to
// dir/name_1.fastq and dir/name_2.fastq.
func writePairs(t *testing.T, f *phantomtesting.E2EFramework, c *phantomtesting.PlatformConfig, dir, name string) (string, string) {
	t.Helper()
	c.PairedEnd = true
	d, err := f.GenerateSyntheticData(c)
	if err != nil {
		t.Fatal(err)
	}
	p1, p2 := filepath.Join(dir, name+"_1.fastq"), filepath.Join(dir, name+"_2.fastq")
	if err := f.WriteFASTQ(d, p1, p2); err != nil {
		t.Fatal(err)
	}
	return p1, p2
}

// TestRunAll_Paired runs the harness on paired stand-in reads generated
// from a reference fixture: the generated reads are pairs, read 2 and the
// pairs are compared and calibrated, the generated reads are written, and
// a second run (three Phred lower, 150 bp) is compared with the first after
// cutting both to 150 bp.
func TestRunAll_Paired(t *testing.T) {
	dir := t.TempDir()
	ref := filepath.Join(dir, "ref.fa")
	if err := phantomtesting.WriteRandomReference(ref, []int{60000}, 4); err != nil {
		t.Fatal(err)
	}
	f := phantomtesting.NewE2EFramework()
	c := f.CreatePlatformConfig(phantomtesting.IlluminaNovaSeq6000, "CEU")
	c.Seed, c.Reads, c.ReadLength, c.ReferencePath = 77, 600, 151, ref
	a1, a2 := writePairs(t, f, c, dir, "a")

	cb := f.CreatePlatformConfig(phantomtesting.IlluminaNovaSeq6000, "CEU")
	cb.Seed, cb.Reads, cb.ReadLength, cb.ReferencePath, cb.PairedEnd = 78, 600, 150, ref, true
	b, err := f.GenerateSyntheticData(cb)
	if err != nil {
		t.Fatal(err)
	}
	b1, b2 := filepath.Join(dir, "b_1.fastq"), filepath.Join(dir, "b_2.fastq")
	if err := f.WriteFASTQ(phantomtesting.QualityShifted(b, 3), b1, b2); err != nil {
		t.Fatal(err)
	}

	args := quickArgs(t)
	args.benchmark, args.reference, args.platforms, args.reads = false, ref, "novaseq", 600
	args.realFASTQ, args.realFASTQ2, args.realFASTQB, args.realFASTQB2 = a1, a2, b1, b2
	args.writeGenerated = filepath.Join(dir, "gen")
	if err := validateArgs(args); err != nil {
		t.Fatal(err)
	}
	r := runAll(f, []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000}, args)
	if r.Failures() != 0 {
		t.Fatalf("failures %d: %+v %s %s", r.Failures(), r.Platforms[0].Errors, r.RealError, r.CrossError)
	}
	p := r.Platforms[0]
	if !r.Paired || !p.Calibration.Paired || len(p.Calibration.Sets) != 4 || p.Calibration.Set(phantomtesting.FeaturePair) == nil {
		t.Fatalf("calibration %+v, want paired sets with the pair set", p.Calibration)
	}
	if s := p.Calibration.Set(phantomtesting.FeaturePair); !s.NullPass || !s.ControlPass {
		t.Errorf("pair calibration %+v, want the null and the quality control passed", s)
	}

	x := p.Realism
	if x == nil || !x.Paired || x.RealReads != 600 || x.GeneratedReads != 600 || x.Read2 == nil || x.Pair == nil {
		t.Fatalf("realism %+v, want 600 real and generated pairs with read 2 and pair comparisons", x)
	}
	if len(x.Read2.RealVsReal) != 3 || len(x.Read2.RealVsGenerated) != 3 || len(x.Pair.RealVsGenerated) != 1 {
		t.Fatalf("read 2 / pair comparisons %d %d %d", len(x.Read2.RealVsReal), len(x.Read2.RealVsGenerated), len(x.Pair.RealVsGenerated))
	}
	if g := x.Pair.RealVsGenerated[0]; g.FeatureSet != phantomtesting.FeaturePair || g.Holdout.Holdout != 600 || x.Pair.RealVsReal[0].Holdout.Holdout != 300 {
		t.Errorf("pair comparison %+v / %+v, want 600 and 300 held out", g, x.Pair.RealVsReal[0])
	}
	// Same model, same reference: the stand-in cannot be told from the
	// generated pairs by any set.
	for _, cs := range [][]*Comparison{x.RealVsGenerated, x.Read2.RealVsGenerated, x.Pair.RealVsGenerated} {
		for _, g := range cs {
			if math.Abs(g.Accuracy-0.5) > g.Band {
				t.Errorf("%s: stand-in vs generated %.4f outside 0.5 ± %.4f", g.FeatureSet, g.Accuracy, g.Band)
			}
		}
	}

	gen, err := phantomtesting.LoadPairedFASTQ(filepath.Join(args.writeGenerated, string(p.Platform)+"_1.fastq"),
		filepath.Join(args.writeGenerated, string(p.Platform)+"_2.fastq"), 0)
	if err != nil || len(gen.Sequences) != 600 || gen.Stats.MaxReadLength != 151 {
		t.Fatalf("written generated pairs: %v", err)
	}

	cr := r.CrossRun
	if cr == nil || cr.ReadLength != 150 || cr.ReadsA != 600 || cr.ReadsB != 600 || len(cr.Read1) != 3 || len(cr.Read2) != 3 || len(cr.Pair) != 1 {
		t.Fatalf("cross run %+v, want both runs cut to 150 bp and every set compared", cr)
	}
	if cr.MeanQualityB > cr.MeanQualityA-2 {
		t.Errorf("cross-run mean qualities %.2f / %.2f, want B 3 lower", cr.MeanQualityA, cr.MeanQualityB)
	}
	for _, cs := range [][]*Comparison{cr.Read1, cr.Read2, cr.Pair} {
		for _, c := range cs {
			seq := c.FeatureSet == phantomtesting.FeatureSequence
			if seq && math.Abs(c.Accuracy-0.5) > c.Band {
				t.Errorf("cross-run %s: %.4f, want chance (same bases model)", c.FeatureSet, c.Accuracy)
			}
			if !seq && c.Accuracy < 0.9 {
				t.Errorf("cross-run %s: %.4f, want > 0.9 (qualities 3 lower)", c.FeatureSet, c.Accuracy)
			}
		}
	}

	var out bytes.Buffer
	printReport(&out, r, args)
	for _, want := range []string{"REAL RUN A VERSUS REAL RUN B", "both cut to 150 bp", "600 real pairs", "  pair ", "read 2: real mean Q"} {
		if !strings.Contains(out.String(), want) {
			t.Errorf("output lacks %q:\n%s", want, out.String())
		}
	}
	files, err := writeReports(r, args.outputDir)
	if err != nil {
		t.Fatal(err)
	}
	fh, err := os.Open(files[1])
	if err != nil {
		t.Fatal(err)
	}
	defer fh.Close()
	rows, err := csv.NewReader(fh).ReadAll()
	if err != nil || len(rows) != 2 || len(rows[0]) != len(rows[1]) {
		t.Fatalf("CSV %v %v", rows, err)
	}
	for i, h := range rows[0] {
		if (h == "real_vs_generated_accuracy_pair" || h == "r2_real_vs_generated_accuracy_quality" || h == "null_accuracy_pair") && rows[1][i] == "" {
			t.Errorf("CSV column %s is empty", h)
		}
	}

	// -real-only keeps the nulls and the cross run, without generating.
	args.realOnly, args.authenticity, args.writeGenerated = true, false, ""
	r = runAll(f, []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000}, args)
	if x := r.Platforms[0].Realism; r.Failures() != 0 || x == nil || x.GeneratedReads != 0 || len(x.RealVsGenerated) != 0 ||
		len(x.Read2.RealVsReal) != 3 || x.Pair.RealVsGenerated != nil || r.CrossRun == nil {
		t.Errorf("-real-only: failures %d, realism %+v", r.Failures(), x)
	}

	// An unpaired run B against paired run A is an error, counted.
	args.realFASTQB2 = ""
	if r := runAll(f, nil, args); r.CrossError == "" || r.Failures() != 1 {
		t.Errorf("paired A vs unpaired B: CrossError %q, failures %d", r.CrossError, r.Failures())
	}
	// Paired real reads with a profile without pairs is a failure.
	args.realFASTQB = ""
	if r := runAll(f, []phantomtesting.SequencingPlatform{phantomtesting.ONTPromethION24}, args); r.Failures() != 1 {
		t.Errorf("paired reads on ONT: failures %d, errors %v", r.Failures(), r.Platforms[0].Errors)
	}
}
