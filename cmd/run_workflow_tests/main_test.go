//go:build !release
// +build !release

package main

import (
	"bytes"
	"encoding/csv"
	"encoding/json"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"testing"

	phantomtesting "github.com/scttfrdmn/phantom-science/internal/testing"
)

// withLookPath replaces the executable lookup for the duration of the test.
// found lists the program names that resolve; everything else is missing.
func withLookPath(t *testing.T, found ...string) {
	t.Helper()
	orig := lookPath
	t.Cleanup(func() { lookPath = orig })
	lookPath = func(name string) (string, error) {
		for _, f := range found {
			if f == name {
				return "/fake/bin/" + name, nil
			}
		}
		return "", exec.ErrNotFound
	}
	t.Setenv("CROMWELL_JAR", "")
}

// withNoEngines makes every workflow engine unavailable for the test.
func withNoEngines(t *testing.T) { withLookPath(t) }

// withFakeEngine makes engine resolve to a shell script with the given body.
func withFakeEngine(t *testing.T, engine, body string) {
	t.Helper()
	if runtime.GOOS == "windows" {
		t.Skip("fake engines are shell scripts")
	}
	script := filepath.Join(t.TempDir(), engine)
	if err := os.WriteFile(script, []byte("#!/bin/sh\n"+body+"\n"), 0700); err != nil { // #nosec G306 -- test executable
		t.Fatal(err)
	}
	orig := lookPath
	t.Cleanup(func() { lookPath = orig })
	lookPath = func(name string) (string, error) {
		if name == engine {
			return script, nil
		}
		return "", exec.ErrNotFound
	}
	t.Setenv("CROMWELL_JAR", "")
}

// requireRealEngine skips unless engine is genuinely installed.
func requireRealEngine(t *testing.T, engine string) {
	t.Helper()
	if testing.Short() {
		t.Skipf("skipping real %s run in -short mode", engine)
	}
	if !isWorkflowEngineAvailable(engine) {
		t.Skipf("%s not installed", engine)
	}
}

func testArgs(t *testing.T) *CommandLineArgs {
	return &CommandLineArgs{outputDir: t.TempDir(), timeout: 600, skipUnavailable: true}
}

func TestIsWorkflowEngineAvailable_NotOnPath(t *testing.T) {
	withNoEngines(t)
	for _, engine := range []string{"cromwell", "nextflow", "snakemake"} {
		if isWorkflowEngineAvailable(engine) {
			t.Errorf("isWorkflowEngineAvailable(%s) = true with nothing on PATH", engine)
		}
		if _, _, err := engineCommand(engine); !errors.Is(err, errEngineUnavailable) {
			t.Errorf("engineCommand(%s) error = %v, want errEngineUnavailable", engine, err)
		}
	}
}

func TestIsWorkflowEngineAvailable_OnPath(t *testing.T) {
	withLookPath(t, "cromwell", "nextflow", "snakemake")
	for _, engine := range []string{"cromwell", "nextflow", "snakemake"} {
		if !isWorkflowEngineAvailable(engine) {
			t.Errorf("isWorkflowEngineAvailable(%s) = false with %s on PATH", engine, engine)
		}
	}
	withLookPath(t, "unknown")
	if isWorkflowEngineAvailable("unknown") {
		t.Error("unknown engines must not be reported available")
	}
}

func TestEngineCommand_CromwellJar(t *testing.T) {
	jar := filepath.Join(t.TempDir(), "cromwell.jar")
	if err := os.WriteFile(jar, []byte("jar"), 0600); err != nil {
		t.Fatal(err)
	}

	withLookPath(t, "java")
	t.Setenv("CROMWELL_JAR", jar)
	program, lead, err := engineCommand("cromwell")
	if err != nil {
		t.Fatalf("engineCommand(cromwell) with CROMWELL_JAR: %v", err)
	}
	if program != "/fake/bin/java" || len(lead) != 2 || lead[0] != "-jar" || lead[1] != jar {
		t.Errorf("engineCommand(cromwell) = %q %q, want java -jar %s", program, lead, jar)
	}

	t.Setenv("CROMWELL_JAR", filepath.Join(t.TempDir(), "missing.jar"))
	if isWorkflowEngineAvailable("cromwell") {
		t.Error("CROMWELL_JAR pointing at a missing file must not count as available")
	}

	withLookPath(t) // no java either
	t.Setenv("CROMWELL_JAR", jar)
	if isWorkflowEngineAvailable("cromwell") {
		t.Error("CROMWELL_JAR without java on PATH must not count as available")
	}
}

func TestParseLists(t *testing.T) {
	for in, want := range map[string][]string{
		"all":            {"cromwell", "nextflow", "snakemake"},
		"wdl,nf":         {"cromwell", "nextflow"},
		"snake":          {"snakemake"},
		"all,cromwell":   {"cromwell", "nextflow", "snakemake"},
		"unknown_engine": {},
	} {
		if got := parseEngineList(in); strings.Join(got, ",") != strings.Join(want, ",") {
			t.Errorf("parseEngineList(%q) = %v, want %v", in, got, want)
		}
	}

	for in, want := range map[string]int{"all": 7, "major": 3, "illumina": 3, "longread": 2, "all,novaseq": 7, "bogus": 0} {
		if got := parsePlatformList(in); len(got) != want {
			t.Errorf("parsePlatformList(%q) = %d platforms, want %d", in, len(got), want)
		}
	}
	for in, want := range map[string]phantomtesting.SequencingPlatform{
		"novaseq": phantomtesting.IlluminaNovaSeq6000, "nextseq": phantomtesting.IlluminaNextSeq2000,
		"miseq": phantomtesting.IlluminaMiSeq, "pacbio": phantomtesting.PacBioRevio, "ont": phantomtesting.ONTPromethION24,
		"dnbseq": phantomtesting.MGI_DNBSEQ_T20x2, "aviti": phantomtesting.ElementAVITI,
	} {
		if got := parsePlatformList(in); len(got) != 1 || got[0] != want {
			t.Errorf("parsePlatformList(%q) = %v, want [%s]", in, got, want)
		}
	}

	for in, want := range map[string][]string{
		"basic":              {"basic"},
		"gatk":               {"gatk_germline", "gatk_somatic"},
		"all,basic":          {"basic", "gatk", "assembly", "variants", "qc", "population"},
		"my-custom-workflow": {"my-custom-workflow"},
		"nfcore":             {"nfcore_sarek", "nfcore_nanoseq", "nfcore_hlatyping"},
	} {
		if got := parseWorkflowList(in); strings.Join(got, ",") != strings.Join(want, ",") {
			t.Errorf("parseWorkflowList(%q) = %v, want %v", in, got, want)
		}
	}
}

func TestValidateArgs(t *testing.T) {
	valid := CommandLineArgs{outputDir: "/tmp/test", timeout: 120, engines: "all", platforms: "all", workflows: "all"}
	if err := validateArgs(&valid); err != nil {
		t.Errorf("valid args: %v", err)
	}
	for name, mutate := range map[string]func(*CommandLineArgs){
		"empty output": func(a *CommandLineArgs) { a.outputDir = "" },
		"low timeout":  func(a *CommandLineArgs) { a.timeout = 30 },
		"engines":      func(a *CommandLineArgs) { a.engines = "invalid" },
		"platforms":    func(a *CommandLineArgs) { a.platforms = "invalid" },
		"workflows":    func(a *CommandLineArgs) { a.workflows = "," },
	} {
		a := valid
		mutate(&a)
		if err := validateArgs(&a); err == nil {
			t.Errorf("%s: want error", name)
		}
	}
}

func TestWorkflowSources(t *testing.T) {
	if wdl := createBasicWDL("gatk", "illumina"); !strings.Contains(wdl, "workflow Gatk {") || !strings.Contains(wdl, readCountFile) {
		t.Errorf("createBasicWDL:\n%s", wdl)
	}
	if inputs := createWDLInputs("gatk", "illumina", "/data/reads.fastq"); !strings.Contains(inputs, `"Gatk.fastq": "/data/reads.fastq"`) {
		t.Errorf("createWDLInputs should bind the FASTQ path, got:\n%s", inputs)
	}
	if nf := createBasicNextflow("variant-calling", "illumina"); !strings.Contains(nf, "process countReads") {
		t.Errorf("createBasicNextflow:\n%s", nf)
	}
	if sm := createBasicSnakemake("qc", "illumina"); !strings.Contains(sm, "rule count_reads") || !strings.Contains(sm, inputFASTQ) {
		t.Errorf("createBasicSnakemake:\n%s", sm)
	}
	for in, want := range map[string]string{
		"gatk":            "Gatk",
		"gatk_germline":   "Gatk_germline",
		"variant-calling": "Variant_calling",
		"1kg":             "W1kg",
		"-x":              "W_x",
		"":                "Workflow",
	} {
		if got := wdlIdentifier(in); got != want {
			t.Errorf("wdlIdentifier(%q) = %q, want %q", in, got, want)
		}
	}
}

func novaSeqData(t *testing.T) *phantomtesting.SyntheticData {
	t.Helper()
	f := phantomtesting.NewE2EFramework()
	data, err := f.GenerateSyntheticData(f.CreatePlatformConfig(phantomtesting.IlluminaNovaSeq6000, "CEU"))
	if err != nil {
		t.Fatalf("GenerateSyntheticData error: %v", err)
	}
	return data
}

// TestCalculateDataIntegrity_ChecksReadCount verifies integrity reflects
// whether the engine's read count matches the reads written. A former
// implementation scored any existing file (even "mock") at 0.95.
func TestCalculateDataIntegrity_ChecksReadCount(t *testing.T) {
	data := novaSeqData(t)
	n := len(data.Sequences)

	dir := t.TempDir()
	write := func(name, content string) string {
		path := filepath.Join(dir, name)
		if err := os.WriteFile(path, []byte(content), 0600); err != nil {
			t.Fatal(err)
		}
		return path
	}
	good := write("good.txt", strconv.Itoa(n)+"\n")
	wrong := write("wrong.txt", strconv.Itoa(n-1)+"\n")
	mock := write("mock.txt", "mock output")
	missing := filepath.Join(dir, "missing.txt")

	cases := []struct {
		name  string
		files []string
		want  float64
	}{
		{"matching count", []string{good}, 1.0},
		{"wrong count", []string{wrong}, 0.0},
		{"mock content", []string{mock}, 0.0},
		{"missing file", []string{missing}, 0.0},
		{"one of two", []string{good, wrong}, 0.5},
		{"no files", nil, 0.0},
	}
	for _, c := range cases {
		if got := calculateDataIntegrity(data, c.files); got != c.want {
			t.Errorf("%s: calculateDataIntegrity = %f, want %f", c.name, got, c.want)
		}
	}
}

type executeFunc func(string, *phantomtesting.SyntheticData, phantomtesting.SequencingPlatform, string, *CommandLineArgs) ([]string, error)

var executors = map[string]executeFunc{
	"cromwell":  executeCromwellWorkflow,
	"nextflow":  executeNextflowWorkflow,
	"snakemake": executeSnakemakeWorkflow,
}

// TestExecuteWorkflow_EngineUnavailable checks that an unavailable engine is
// an error and that no placeholder output is written.
func TestExecuteWorkflow_EngineUnavailable(t *testing.T) {
	withNoEngines(t)
	data := novaSeqData(t)
	for engine, execute := range executors {
		t.Run(engine, func(t *testing.T) {
			dir := t.TempDir()
			files, err := execute(dir, data, phantomtesting.IlluminaNovaSeq6000, "basic", &CommandLineArgs{outputDir: dir})
			if !errors.Is(err, errEngineUnavailable) {
				t.Errorf("error = %v, want errEngineUnavailable", err)
			}
			if len(files) != 0 {
				t.Errorf("returned output files %v for an engine that did not run", files)
			}
		})
	}
}

// TestExecuteWorkflow_RealEngine runs each installed engine on the
// generated reads and checks the reported read count matches what was
// written.
func TestExecuteWorkflow_RealEngine(t *testing.T) {
	data := novaSeqData(t)
	for _, engine := range []string{"snakemake", "nextflow", "cromwell"} {
		execute := executors[engine]
		t.Run(engine, func(t *testing.T) {
			requireRealEngine(t, engine)
			dir := t.TempDir()
			files, err := execute(dir, data, phantomtesting.IlluminaNovaSeq6000, "basic", &CommandLineArgs{outputDir: dir, timeout: 600})
			if err != nil {
				t.Fatalf("%s run failed: %v", engine, err)
			}
			if got := calculateDataIntegrity(data, files); got != 1.0 {
				t.Errorf("integrity = %f, want 1.0 (files %v)", got, files)
			}
		})
	}
}

func TestRunWorkflowTest_Unavailable(t *testing.T) {
	withNoEngines(t)
	f := phantomtesting.NewE2EFramework()
	for _, engine := range []string{"cromwell", "nextflow", "snakemake"} {
		r := runWorkflowTest(f, engine, phantomtesting.IlluminaNovaSeq6000, "basic", testArgs(t))
		if r.Status != statusSkipped || r.ExitStatus != nil || r.Reason == "" {
			t.Errorf("%s unavailable with skip-unavailable: %+v, want skipped with a reason and no exit status", engine, r)
		}
	}
	args := testArgs(t)
	args.skipUnavailable = false
	if r := runWorkflowTest(f, "nextflow", phantomtesting.IlluminaNovaSeq6000, "basic", args); r.Status != statusFailed || r.Reason == "" {
		t.Errorf("unavailable without skip: %+v, want failed with a reason", r)
	}
	if r := runWorkflowTest(f, "unsupported_engine", phantomtesting.IlluminaNovaSeq6000, "basic", testArgs(t)); r.Status != statusFailed {
		t.Errorf("unsupported engine: status %s, want failed (not skipped)", r.Status)
	}
}

// TestRunWorkflowTest_NotModeled checks that a platform the generator has no
// model for is skipped, not run on fabricated reads.
func TestRunWorkflowTest_NotModeled(t *testing.T) {
	withLookPath(t, "snakemake")
	f := phantomtesting.NewE2EFramework()
	for _, p := range []phantomtesting.SequencingPlatform{phantomtesting.PacBioRevio, phantomtesting.ElementAVITI, phantomtesting.MGI_DNBSEQ_T20x2} {
		r := runWorkflowTest(f, "snakemake", p, "basic", testArgs(t))
		if r.Status != statusSkipped || !strings.Contains(r.Reason, phantomtesting.ErrPlatformNotModeled.Error()) || r.Reads != 0 {
			t.Errorf("%s: %+v, want skipped as not modeled", p, r)
		}
	}
}

// TestRunWorkflowTest_FakeEngineFailure checks that a failing engine's exit
// status is recorded.
func TestRunWorkflowTest_FakeEngineFailure(t *testing.T) {
	withFakeEngine(t, "snakemake", "echo boom >&2\nexit 3")
	r := runWorkflowTest(phantomtesting.NewE2EFramework(), "snakemake", phantomtesting.IlluminaNovaSeq6000, "basic", testArgs(t))
	if r.Status != statusFailed || r.ExitStatus == nil || *r.ExitStatus != 3 || !strings.Contains(r.Reason, "boom") {
		t.Errorf("%+v (exit %v), want failed with exit status 3 and the engine's output", r, r.ExitStatus)
	}
	if r.Reads != 1000 || r.WallSeconds <= 0 {
		t.Errorf("Reads %d, WallSeconds %f; want 1000 and a measured time", r.Reads, r.WallSeconds)
	}
}

// TestRunWorkflowTest_FakeEngineWrongCount checks that an engine exiting 0
// with the wrong read count fails on data integrity.
func TestRunWorkflowTest_FakeEngineWrongCount(t *testing.T) {
	// Arguments: --cores 1 --snakefile <file> --directory <workDir>.
	withFakeEngine(t, "snakemake", `echo 7 > "$6/`+readCountFile+`"`)
	r := runWorkflowTest(phantomtesting.NewE2EFramework(), "snakemake", phantomtesting.IlluminaNovaSeq6000, "basic", testArgs(t))
	if r.Status != statusFailed || r.ExitStatus == nil || *r.ExitStatus != 0 {
		t.Fatalf("%+v, want failed with exit status 0", r)
	}
	if r.DataIntegrity != 0 || len(r.ReportedReads) != 1 || r.ReportedReads[0] != "7" {
		t.Errorf("integrity %f, reported %v; want 0 and [7]", r.DataIntegrity, r.ReportedReads)
	}
}

// TestRunWorkflowTest_RealEngines runs each installed engine through
// runWorkflowTest on Illumina and ONT reads.
func TestRunWorkflowTest_RealEngines(t *testing.T) {
	f := phantomtesting.NewE2EFramework()
	for _, engine := range []string{"snakemake", "nextflow", "cromwell"} {
		for _, p := range []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000, phantomtesting.ONTPromethION24} {
			t.Run(engine+"/"+string(p), func(t *testing.T) {
				requireRealEngine(t, engine)
				r := runWorkflowTest(f, engine, p, "basic", testArgs(t))
				if r.Status != statusPassed {
					t.Fatalf("%s run failed: %s", engine, r.Reason)
				}
				want := strconv.Itoa(f.CreatePlatformConfig(p, "CEU").Reads)
				if r.DataIntegrity != 1.0 || r.ExitStatus == nil || *r.ExitStatus != 0 || len(r.ReportedReads) == 0 || r.ReportedReads[0] != want {
					t.Errorf("%+v, want integrity 1, exit 0 and %s reported reads", r, want)
				}
			})
		}
	}
}

// TestRunWorkflowTests_DryRun is a regression test: with -dry-run nothing
// was recorded, so the summary said no engine was available. Dry-run runs
// must be counted as planned (or skipped where they would be), nothing may
// execute, and no pass/fail verdict may be printed.
func TestRunWorkflowTests_DryRun(t *testing.T) {
	withLookPath(t, "cromwell", "nextflow")
	f := phantomtesting.NewE2EFramework()
	args := testArgs(t)
	args.dryRun = true
	engines := []string{"cromwell", "nextflow", "snakemake"}
	platforms := []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000, phantomtesting.ONTPromethION24, phantomtesting.PacBioRevio}
	workflows := []string{"basic", "gatk_germline"}
	results := runWorkflowTests(f, engines, platforms, workflows, args)

	// cromwell and nextflow × 2 modeled platforms × 2 workflows are planned;
	// snakemake (not installed) and PacBio (not modeled) runs are skipped.
	c := results.Counts()
	if want := (OutcomeCounts{Planned: 8, Skipped: 10}); c != want {
		t.Errorf("counts %+v, want %+v", c, want)
	}
	for _, r := range results.Runs {
		if r.ExitStatus != nil || r.Reads != 0 || len(r.OutputFiles) != 0 {
			t.Errorf("dry run executed %s/%s/%s: %+v", r.Engine, r.Platform, r.Workflow, r)
		}
	}
	if entries, _ := os.ReadDir(args.outputDir); len(entries) != 0 {
		t.Errorf("dry run created workspaces: %v", entries)
	}

	var out bytes.Buffer
	printTestSummary(&out, results, args)
	text := out.String()
	if !strings.Contains(text, "Dry run: 8 planned, 10 would be skipped") {
		t.Errorf("summary lacks the planned count:\n%s", text)
	}
	for _, banned := range []string{"no workflow engine was available", "No workflow was executed", "Passed", "Score", "EXCELLENT", "POOR"} {
		if strings.Contains(text, banned) {
			t.Errorf("dry-run summary contains %q:\n%s", banned, text)
		}
	}
}

// TestRunWorkflowTests_NoEngines checks that with nothing installed every
// run is skipped, nothing passes, and the summary says so.
func TestRunWorkflowTests_NoEngines(t *testing.T) {
	withNoEngines(t)
	args := testArgs(t)
	results := runWorkflowTests(phantomtesting.NewE2EFramework(), []string{"cromwell", "nextflow", "snakemake"},
		[]phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000}, []string{"basic"}, args)
	if c := results.Counts(); c != (OutcomeCounts{Skipped: 3}) {
		t.Errorf("counts %+v, want 3 skipped", c)
	}
	var out bytes.Buffer
	printTestSummary(&out, results, args)
	if !strings.Contains(out.String(), "Passed 0, failed 0, skipped 3.") || !strings.Contains(out.String(), "No workflow was executed") {
		t.Errorf("summary:\n%s", out.String())
	}
}

func TestCounts(t *testing.T) {
	r := &WorkflowTestResults{Runs: []*WorkflowTestResult{
		{Status: statusPassed}, {Status: statusFailed}, {Status: statusPassed},
		{Status: statusSkipped}, {Status: statusPlanned},
	}}
	if c := r.Counts(); c != (OutcomeCounts{Passed: 2, Failed: 1, Skipped: 1, Planned: 1}) {
		t.Errorf("Counts = %+v", c)
	}
}

// TestWriteReports checks that the JSON and CSV reports are written and
// hold every run.
func TestWriteReports(t *testing.T) {
	withFakeEngine(t, "snakemake", "exit 2")
	args := testArgs(t)
	results := runWorkflowTests(phantomtesting.NewE2EFramework(), []string{"snakemake", "nextflow"},
		[]phantomtesting.SequencingPlatform{phantomtesting.IlluminaMiSeq}, []string{"basic"}, args)
	files, err := writeReports(results, args.outputDir)
	if err != nil {
		t.Fatal(err)
	}

	b, err := os.ReadFile(files[0])
	if err != nil {
		t.Fatal(err)
	}
	var back WorkflowTestResults
	if err := json.Unmarshal(b, &back); err != nil {
		t.Fatalf("JSON report does not parse: %v", err)
	}
	if len(back.Runs) != 2 || back.Runs[0].Status != statusFailed || back.Runs[0].ExitStatus == nil || *back.Runs[0].ExitStatus != 2 ||
		back.Runs[1].Status != statusSkipped {
		t.Errorf("JSON runs %+v", back.Runs)
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
	if len(rows) != 3 || rows[1][3] != statusFailed || rows[1][5] != "2" || rows[2][3] != statusSkipped {
		t.Errorf("CSV rows %v", rows)
	}
}
