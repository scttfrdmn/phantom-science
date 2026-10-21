//go:build !release
// +build !release

package main

import (
	"bytes"
	"encoding/csv"
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	phantomtesting "github.com/scttfrdmn/phantom-science/internal/testing"
)

func TestParsePlatformList(t *testing.T) {
	if got := parsePlatformList("all"); len(got) != 7 {
		t.Errorf("all: %d platforms, want 7", len(got))
	}
	if got := parsePlatformList("illumina"); len(got) != 3 {
		t.Errorf("illumina: %d platforms, want 3", len(got))
	}
	if got := parsePlatformList("novaseq, illumina,novaseq"); len(got) != 3 || got[0] != phantomtesting.IlluminaNovaSeq6000 {
		t.Errorf("duplicates: %v", got)
	}
	if got := parsePlatformList("nonsense"); len(got) != 0 {
		t.Errorf("nonsense: %v", got)
	}
}

func TestValidateArgs(t *testing.T) {
	ref := filepath.Join(t.TempDir(), "ref.fa")
	if err := phantomtesting.WriteRandomReference(ref, []int{1000}, 1); err != nil {
		t.Fatal(err)
	}
	for _, ok := range []*CommandLineArgs{
		{platforms: "novaseq", outputDir: t.TempDir(), report: true, seed: 1},
		{platforms: "novaseq", seed: 1, reference: ref, region: "ctg1:1-500", reads: 100},
	} {
		if err := validateArgs(ok); err != nil {
			t.Errorf("valid args %+v: %v", ok, err)
		}
	}
	bad := []*CommandLineArgs{
		{platforms: "novaseq", report: true, seed: 1},
		{platforms: "bogus", seed: 1},
		{platforms: "novaseq", seed: 0},
		{platforms: "novaseq", seed: 1, realFASTQ: filepath.Join(t.TempDir(), "missing.fastq")},
		{platforms: "novaseq", seed: 1, region: "chr1"},
		{platforms: "novaseq", seed: 1, reference: filepath.Join(t.TempDir(), "missing.fa")},
		{platforms: "novaseq", seed: 1, reads: 3},
		{platforms: "novaseq", seed: 1, reads: -1},
	}
	for i, a := range bad {
		if err := validateArgs(a); err == nil {
			t.Errorf("case %d: want error for %+v", i, a)
		}
	}
}

func TestGradeSet(t *testing.T) {
	band := 4 * math.Sqrt(0.25/1000)
	cases := []struct {
		set                phantomtesting.FeatureSet
		null, gc, qual     float64
		nullPass, ctrlPass bool
	}{
		{phantomtesting.FeatureSequence, 0.5, 0.99, 0.5, true, true},
		{phantomtesting.FeatureSequence, 0.5 + band*0.99, 0.91, 0.5, true, true},
		{phantomtesting.FeatureSequence, 0.5 - band*1.01, 0.99, 0.5, false, true},
		{phantomtesting.FeatureSequence, 0.5, 0.9, 0.99, true, false},
		{phantomtesting.FeatureQuality, 0.5, 0.5, 0.99, true, true},
		{phantomtesting.FeatureQuality, 0.5, 0.99, 0.9, true, false},
		{phantomtesting.FeatureCombined, 0.5, 0.99, 0.99, true, true},
		{phantomtesting.FeatureCombined, 0.5, 0.99, 0.5, true, false},
		{phantomtesting.FeatureCombined, 0.5, 0.5, 0.99, true, false},
	}
	for _, c := range cases {
		g := gradeSet(c.set, 1000, c.null, c.gc, c.qual)
		if g.NullPass != c.nullPass || g.ControlPass != c.ctrlPass {
			t.Errorf("%s null %.4f gc %.2f qual %.2f: pass %v/%v, want %v/%v", c.set, c.null, c.gc, c.qual, g.NullPass, g.ControlPass, c.nullPass, c.ctrlPass)
		}
		if math.Abs(g.Band-band) > 1e-15 {
			t.Errorf("band %v, want %v", g.Band, band)
		}
	}
	cal := &Calibration{}
	if cal.Passed() {
		t.Error("a calibration without feature sets passed")
	}
	cal.Sets = []*SetCalibration{gradeSet(phantomtesting.FeatureSequence, 1000, 0.5, 0.99, 0.5), gradeSet(phantomtesting.FeatureQuality, 1000, 0.5, 0.5, 0.5)}
	if cal.Passed() {
		t.Error("a calibration with a failed set passed")
	}
}

// quickArgs returns arguments for a fast run.
func quickArgs(t *testing.T) *CommandLineArgs {
	return &CommandLineArgs{outputDir: t.TempDir(), seed: 1, quick: true, benchmark: true, authenticity: true, report: true}
}

// TestRunPlatform_Measured runs a modeled platform end to end and checks
// the reported values come from the generated reads.
func TestRunPlatform_Measured(t *testing.T) {
	f := phantomtesting.NewE2EFramework()
	args := quickArgs(t)
	for _, p := range []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000, phantomtesting.ONTPromethION24} {
		r := runPlatform(f, p, nil, args)
		if r.Status != statusMeasured || len(r.Errors) != 0 {
			t.Fatalf("%s: status %s, errors %v", p, r.Status, r.Errors)
		}
		c := f.CreatePlatformConfig(p, "CEU")
		data, err := f.GenerateSyntheticData(c)
		if err != nil {
			t.Fatal(err)
		}
		if *r.Reads != data.Stats {
			t.Errorf("%s: reported reads %+v, want the seed-%d data's stats %+v", p, *r.Reads, c.Seed, data.Stats)
		}
		if r.Benchmark == nil || r.Benchmark.WallSeconds <= 0 || r.Benchmark.BasesPerSecond <= 0 {
			t.Errorf("%s: benchmark %+v, want measured positive values", p, r.Benchmark)
		}
		if r.Calibration == nil || !r.Calibration.Passed() || len(r.Calibration.Sets) != len(phantomtesting.AllFeatureSets) {
			t.Fatalf("%s: calibration %+v, want every feature set's checks passed", p, r.Calibration)
		}
		for _, s := range r.Calibration.Sets {
			if s.Holdout != quickCalibrationReads {
				t.Errorf("%s %s: holdout %d, want %d (half of each %d-read set)", p, s.FeatureSet, s.Holdout, quickCalibrationReads, quickCalibrationReads)
			}
		}
		if r.Calibration.Bases.Source != phantomtesting.BasesUniform || r.Bases.Source != phantomtesting.BasesUniform {
			t.Errorf("%s: bases %+v, want uniform", p, r.Bases)
		}
		if r.Realism != nil {
			t.Errorf("%s: realism reported without real reads", p)
		}
	}
}

// TestRunPlatform_NotModeled checks that unmodeled platforms are reported
// as such, with nothing measured and no failure.
func TestRunPlatform_NotModeled(t *testing.T) {
	f := phantomtesting.NewE2EFramework()
	for _, p := range []phantomtesting.SequencingPlatform{phantomtesting.PacBioRevio, phantomtesting.MGI_DNBSEQ_T20x2, phantomtesting.ElementAVITI} {
		r := runPlatform(f, p, nil, quickArgs(t))
		if r.Status != statusNotModeled || r.Reads != nil || r.Benchmark != nil || r.Calibration != nil || r.failures() != 0 {
			t.Errorf("%s: %+v, want not_modeled with nothing measured and no failure", p, r)
		}
	}
}

// TestRunAll_RealFASTQ checks the real-versus-generated comparison on a
// supplied FASTQ (here generated reads, GC-skewed and 151 bp long,
// standing in for real ones).
func TestRunAll_RealFASTQ(t *testing.T) {
	f := phantomtesting.NewE2EFramework()
	c := f.CreatePlatformConfig(phantomtesting.IlluminaNovaSeq6000, "CEU")
	c.Seed, c.Reads, c.ReadLength = 99, 300, 151
	stand, err := f.GenerateSyntheticData(c)
	if err != nil {
		t.Fatal(err)
	}
	dir := t.TempDir()
	path := filepath.Join(dir, "real.fastq")
	if err := f.WriteFASTQ(phantomtesting.GCSkewed(stand, 0.7, 1), path, filepath.Join(dir, "r2.fastq")); err != nil {
		t.Fatal(err)
	}

	args := quickArgs(t)
	args.benchmark, args.authenticity, args.realFASTQ = false, false, path
	r := runAll(f, []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000, phantomtesting.PacBioRevio}, args)
	if r.Failures() != 0 {
		t.Fatalf("failures %d: %+v", r.Failures(), r.Platforms[0].Errors)
	}
	x := r.Platforms[0].Realism
	if x == nil || x.RealReads != 300 || x.GeneratedReads != 300 || x.GeneratedReadLength != 151 || x.Bases.Source != phantomtesting.BasesUniform {
		t.Fatalf("realism %+v, want 300 real and generated 151 bp uniform reads", x)
	}
	for _, set := range phantomtesting.AllFeatureSets {
		g, n := find(x.RealVsGenerated, set), find(x.RealVsReal, set)
		if g == nil || n == nil || g.Holdout.Holdout != 300 || n.Holdout.Holdout != 150 {
			t.Fatalf("%s: comparisons %+v / %+v, want 300 and 150 held out", set, g, n)
		}
		if len(g.Top) != topFeatures {
			t.Errorf("%s: %d top features, want %d", set, len(g.Top), topFeatures)
		}
	}
	if acc := find(x.RealVsGenerated, phantomtesting.FeatureSequence).Accuracy; acc < 0.9 {
		t.Errorf("GC-skewed stand-in vs generated reads: sequence accuracy %.3f, want > 0.9", acc)
	}
	if r.Platforms[1].Realism != nil {
		t.Error("unmodeled platform has a realism result")
	}

	args.realFASTQ = filepath.Join(dir, "r2.fastq")
	if err := os.WriteFile(args.realFASTQ, []byte("not fastq\n"), 0600); err != nil {
		t.Fatal(err)
	}
	if r := runAll(f, []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000}, args); r.RealError == "" || r.Failures() != 1 {
		t.Errorf("unreadable real FASTQ: RealError %q, failures %d; want an error counted as a failure", r.RealError, r.Failures())
	}
}

// TestRunAll_Reference runs the harness with -reference and -region on a
// small fixture: reads standing in for real ones are generated from the
// fixture, the generated reads come from the same region, and the report
// records the mode.
func TestRunAll_Reference(t *testing.T) {
	dir := t.TempDir()
	ref := filepath.Join(dir, "ref.fa")
	if err := phantomtesting.WriteRandomReference(ref, []int{4000, 20000}, 3); err != nil {
		t.Fatal(err)
	}
	f := phantomtesting.NewE2EFramework()
	c := f.CreatePlatformConfig(phantomtesting.IlluminaNovaSeq6000, "CEU")
	c.Seed, c.Reads, c.ReferencePath, c.ReferenceRegion = 77, 400, ref, "ctg2"
	stand, err := f.GenerateSyntheticData(c)
	if err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(dir, "real.fastq")
	if err := f.WriteFASTQ(stand, path, filepath.Join(dir, "r2.fastq")); err != nil {
		t.Fatal(err)
	}

	args := quickArgs(t)
	args.benchmark, args.realFASTQ, args.reference, args.region = false, path, ref, "ctg2"
	r := runAll(f, []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000}, args)
	if r.Failures() != 0 {
		t.Fatalf("failures %d: %+v", r.Failures(), r.Platforms[0].Errors)
	}
	want := Bases{Source: phantomtesting.BasesReference, Reference: ref, Region: "ctg2"}
	p := r.Platforms[0]
	if r.Bases != want || p.Bases != want || p.Calibration.Bases != want || p.Realism.Bases != want {
		t.Errorf("bases %+v / %+v / %+v / %+v, want %+v", r.Bases, p.Bases, p.Calibration.Bases, p.Realism.Bases, want)
	}
	// Generated from the same region with the same model, the stand-in
	// cannot be told from the generated reads by any feature set.
	for _, g := range p.Realism.RealVsGenerated {
		if math.Abs(g.Accuracy-0.5) > g.Band {
			t.Errorf("%s: stand-in from the reference vs generated %.4f outside 0.5 ± %.4f", g.FeatureSet, g.Accuracy, g.Band)
		}
	}
	var out bytes.Buffer
	printReport(&out, r, args)
	if !strings.Contains(out.String(), "generated from reference "+ref+" region ctg2") {
		t.Errorf("output does not name the reference mode:\n%s", out.String())
	}
}

func TestReportFailures(t *testing.T) {
	pass := &Calibration{Sets: []*SetCalibration{gradeSet(phantomtesting.FeatureSequence, 1000, 0.5, 0.99, 0.5)}}
	fail := &Calibration{Sets: []*SetCalibration{gradeSet(phantomtesting.FeatureSequence, 1000, 0.7, 0.99, 0.5)}}
	badNull := &Realism{RealVsReal: []*Comparison{{Holdout: &phantomtesting.Holdout{Accuracy: 0.7, Holdout: 1000}, Band: 0.06}}}
	goodNull := &Realism{RealVsReal: []*Comparison{{Holdout: &phantomtesting.Holdout{Accuracy: 0.52, Holdout: 1000}, Band: 0.06}}}
	r := &Report{Platforms: []*PlatformResult{
		{Status: statusNotModeled},
		{Status: statusMeasured, Calibration: pass, Realism: goodNull},
		{Status: statusMeasured, Calibration: fail},
		{Status: statusMeasured, Calibration: pass, Realism: badNull},
		{Status: statusError, Errors: []string{"generate: boom"}},
	}}
	if got := r.Failures(); got != 3 {
		t.Errorf("Failures = %d, want 3 (one failed calibration, one real-vs-real outside the band, one error)", got)
	}
}

// TestPrintAndWriteReports checks that the printed table carries the
// measured values and no score or verdict, and that the report files are
// written and parse.
func TestPrintAndWriteReports(t *testing.T) {
	f := phantomtesting.NewE2EFramework()
	args := quickArgs(t)
	r := runAll(f, []phantomtesting.SequencingPlatform{phantomtesting.IlluminaMiSeq, phantomtesting.ElementAVITI}, args)

	var out bytes.Buffer
	printReport(&out, r, args)
	text := out.String()
	for _, want := range []string{"illumina_miseq", "NOT MODELED", "PASS", "uniform random bases", "sequence", "quality", "combined", "No real reads supplied", "Failures: 0"} {
		if !strings.Contains(text, want) {
			t.Errorf("output lacks %q:\n%s", want, text)
		}
	}
	for _, banned := range []string{"Overall Score", "EXCELLENT", "FAIR", "GOOD:", "POOR", "Authenticity"} {
		if strings.Contains(text, banned) {
			t.Errorf("output contains heuristic %q", banned)
		}
	}

	files, err := writeReports(r, args.outputDir)
	if err != nil {
		t.Fatal(err)
	}
	b, err := os.ReadFile(files[0])
	if err != nil {
		t.Fatal(err)
	}
	var back Report
	if err := json.Unmarshal(b, &back); err != nil {
		t.Fatalf("JSON report does not parse: %v", err)
	}
	if len(back.Platforms) != 2 || back.Platforms[0].Reads == nil || *back.Platforms[0].Reads != *r.Platforms[0].Reads ||
		back.Platforms[0].Calibration == nil || len(back.Platforms[0].Calibration.Sets) != 3 {
		t.Errorf("JSON report platforms %+v do not match the run", back.Platforms)
	}
	fh, err := os.Open(files[1])
	if err != nil {
		t.Fatal(err)
	}
	defer fh.Close()
	rows, err := csv.NewReader(fh).ReadAll()
	if err != nil {
		t.Fatalf("CSV does not parse: %v", err)
	}
	if len(rows) != 3 || rows[1][0] != "illumina_miseq" || rows[2][1] != statusNotModeled || len(rows[0]) != len(rows[1]) {
		t.Fatalf("CSV rows %v", rows)
	}
	col := func(name string) string {
		for i, h := range rows[0] {
			if h == name {
				return rows[1][i]
			}
		}
		t.Fatalf("no CSV column %s", name)
		return ""
	}
	if col("bases_source") != "uniform" || col("calibration_pass") != "true" || col("null_accuracy_quality") == "" || col("real_vs_generated_accuracy_combined") != "" {
		t.Errorf("CSV row %v", rows[1])
	}
}
