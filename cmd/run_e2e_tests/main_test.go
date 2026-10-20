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
	ok := &CommandLineArgs{platforms: "novaseq", outputDir: t.TempDir(), report: true, seed: 1}
	if err := validateArgs(ok); err != nil {
		t.Errorf("valid args: %v", err)
	}
	bad := []*CommandLineArgs{
		{platforms: "novaseq", report: true, seed: 1},
		{platforms: "bogus", seed: 1},
		{platforms: "novaseq", seed: 0},
		{platforms: "novaseq", seed: 1, realFASTQ: filepath.Join(t.TempDir(), "missing.fastq")},
	}
	for i, a := range bad {
		if err := validateArgs(a); err == nil {
			t.Errorf("case %d: want error for %+v", i, a)
		}
	}
}

func TestGradeCalibration(t *testing.T) {
	band := 4 * math.Sqrt(0.25/1000)
	cases := []struct {
		null, control      float64
		nullPass, ctrlPass bool
	}{
		{0.5, 0.99, true, true},
		{0.5 + band*0.99, 0.91, true, true},
		{0.5 - band*1.01, 0.99, false, true},
		{0.5, 0.9, true, false},
	}
	for _, c := range cases {
		g := gradeCalibration(1000, 1000, c.null, c.control)
		if g.NullPass != c.nullPass || g.ControlPass != c.ctrlPass || g.Passed() != (c.nullPass && c.ctrlPass) {
			t.Errorf("null %.4f control %.2f: pass %v/%v, want %v/%v", c.null, c.control, g.NullPass, g.ControlPass, c.nullPass, c.ctrlPass)
		}
		if math.Abs(g.Band-band) > 1e-15 {
			t.Errorf("band %v, want %v", g.Band, band)
		}
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
		if r.Calibration == nil || !r.Calibration.Passed() {
			t.Errorf("%s: calibration %+v, want both checks passed", p, r.Calibration)
		}
		if r.Calibration != nil && r.Calibration.Holdout != quickCalibrationReads {
			t.Errorf("%s: holdout %d, want %d (half of each %d-read set)", p, r.Calibration.Holdout, quickCalibrationReads, quickCalibrationReads)
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
// supplied FASTQ (here generated reads, standing in for real ones).
func TestRunAll_RealFASTQ(t *testing.T) {
	f := phantomtesting.NewE2EFramework()
	c := f.CreatePlatformConfig(phantomtesting.IlluminaNovaSeq6000, "CEU")
	c.Seed, c.Reads = 99, 300
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
	if x == nil || x.RealReads != 300 || x.Holdout != 300 {
		t.Fatalf("realism %+v, want 300 real reads and 300 held out", x)
	}
	if x.Accuracy < 0.9 {
		t.Errorf("GC-skewed stand-in vs generated reads: accuracy %.3f, want > 0.9", x.Accuracy)
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

func TestReportFailures(t *testing.T) {
	r := &Report{Platforms: []*PlatformResult{
		{Status: statusNotModeled},
		{Status: statusMeasured, Calibration: gradeCalibration(1000, 1000, 0.5, 0.99)},
		{Status: statusMeasured, Calibration: gradeCalibration(1000, 1000, 0.7, 0.99)},
		{Status: statusError, Errors: []string{"generate: boom"}},
	}}
	if got := r.Failures(); got != 2 {
		t.Errorf("Failures = %d, want 2 (one failed calibration, one error)", got)
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
	for _, want := range []string{"illumina_miseq", "NOT MODELED", "PASS", "No real reads supplied", "Failures: 0"} {
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
	if len(back.Platforms) != 2 || back.Platforms[0].Reads == nil || *back.Platforms[0].Reads != *r.Platforms[0].Reads {
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
	if len(rows) != 3 || rows[1][0] != "illumina_miseq" || rows[2][1] != statusNotModeled {
		t.Errorf("CSV rows %v", rows)
	}
}
