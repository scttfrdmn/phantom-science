// Command run_workflow_tests runs generated reads through Cromwell (WDL),
// Nextflow and Snakemake and reports, for each run, the engine's exit
// status, the wall time and whether the read count the workflow reports
// matches the reads written. It computes no scores.
package main

import (
	"bufio"
	"context"
	"encoding/csv"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"time"
	"unicode"

	"github.com/scttfrdmn/phantom-science/internal/testing"
)

// CommandLineArgs holds command line arguments for workflow testing.
type CommandLineArgs struct {
	engines         string
	platforms       string
	workflows       string
	outputDir       string
	timeout         int
	verbose         bool
	dryRun          bool
	skipUnavailable bool
	generateReport  bool
	help            bool
}

func main() {
	args := parseCommandLineArgs()
	if args.help {
		printUsage()
		return
	}
	if err := validateArgs(args); err != nil {
		log.Fatalf("invalid arguments: %v", err)
	}
	if err := os.MkdirAll(args.outputDir, 0750); err != nil {
		log.Fatalf("create output directory %s: %v", args.outputDir, err)
	}

	results := runWorkflowTests(testing.NewE2EFramework(),
		parseEngineList(args.engines), parsePlatformList(args.platforms), parseWorkflowList(args.workflows), args)
	printTestSummary(os.Stdout, results, args)

	if args.generateReport {
		files, err := writeReports(results, args.outputDir)
		if err != nil {
			log.Fatalf("write reports: %v", err)
		}
		for _, f := range files {
			fmt.Printf("wrote %s\n", f)
		}
	}
	if c := results.Counts(); c.Failed > 0 {
		log.Fatalf("%d of %d executed workflow runs failed", c.Failed, c.Passed+c.Failed)
	}
}

// parseCommandLineArgs parses command line arguments.
func parseCommandLineArgs() *CommandLineArgs {
	args := &CommandLineArgs{}
	flag.StringVar(&args.engines, "engines", "all",
		"Comma-separated list of workflow engines to test (all, cromwell, nextflow, snakemake)")
	flag.StringVar(&args.platforms, "platforms", "major",
		"Comma-separated list of platforms to test (all, major, illumina, longread, or a single platform)")
	flag.StringVar(&args.workflows, "workflows", "basic",
		"Comma-separated list of workflow names (all, basic, gatk, assembly, variants, qc, population, nfcore, or custom names)")
	flag.StringVar(&args.outputDir, "output", "./workflow_test_results",
		"Directory for workspaces and reports")
	flag.IntVar(&args.timeout, "timeout", 3600, "Timeout for each workflow run in seconds")
	flag.BoolVar(&args.verbose, "verbose", false, "Enable verbose logging")
	flag.BoolVar(&args.dryRun, "dry-run", false, "List the runs that would be made without executing them")
	flag.BoolVar(&args.skipUnavailable, "skip-unavailable", true,
		"Report runs on engines that are not installed as skipped (false: as failed)")
	flag.BoolVar(&args.generateReport, "report", true,
		"Write workflow_results.json and workflow_results.csv to -output")
	flag.BoolVar(&args.help, "help", false, "Show this help message")
	flag.Parse()
	return args
}

// printUsage prints usage information.
func printUsage() {
	fmt.Println(`Usage: run_workflow_tests [OPTIONS]

Runs generated reads through Cromwell (WDL), Nextflow and Snakemake for each
selected engine, platform profile and workflow name.

OPTIONS:`)
	flag.PrintDefaults()
	fmt.Print(`
WORKFLOW NAMES:
  all       - basic, gatk, assembly, variants, qc and population
  gatk      - gatk_germline, gatk_somatic
  assembly  - pacbio_assembly, ont_assembly
  variants  - variant_calling, variant_annotation
  qc        - quality_control, multiqc
  population- population_genomics, admixture
  Every name runs the same read-count pipeline; names only label runs.

PLATFORMS:
  all, major (NovaSeq, PacBio, ONT), illumina, longread, novaseq, nextseq,
  miseq, pacbio, ont, dnbseq, aviti. PacBio, MGI and AVITI have no read
  model in the generator; their runs are reported as skipped.

WHAT IS EXECUTED:
  Each engine found on PATH (cromwell, or java with CROMWELL_JAR set;
  nextflow; snakemake) runs a minimal pipeline that reads the generated FASTQ
  and writes its record count. A run passes when the engine exits 0 and the
  count equals the reads written. The table reports, per run: status, wall
  time (including engine start-up), the engine's exit status, the reads
  written, the count the workflow reported, and data integrity (the fraction
  of count outputs that match).

  Engines that are not installed are reported as SKIPPED (with
  -skip-unavailable=false, as FAILED). -dry-run lists the runs as PLANNED
  and executes nothing. The command exits non-zero if any executed run
  fails.
`)
}

// validateArgs validates command line arguments.
func validateArgs(args *CommandLineArgs) error {
	if args.outputDir == "" {
		return fmt.Errorf("output directory cannot be empty")
	}
	if args.timeout < 60 {
		return fmt.Errorf("timeout must be at least 60 seconds")
	}
	if len(parseEngineList(args.engines)) == 0 {
		return fmt.Errorf("no valid workflow engines specified")
	}
	if len(parsePlatformList(args.platforms)) == 0 {
		return fmt.Errorf("no valid platforms specified")
	}
	if len(parseWorkflowList(args.workflows)) == 0 {
		return fmt.Errorf("no valid workflows specified")
	}
	return nil
}

// dedupe returns items without repeats, in first-seen order.
func dedupe[T comparable](items []T) []T {
	seen := make(map[T]bool)
	out := make([]T, 0, len(items))
	for _, it := range items {
		if !seen[it] {
			seen[it] = true
			out = append(out, it)
		}
	}
	return out
}

// parseEngineList parses the workflow engine list argument.
func parseEngineList(enginesArg string) []string {
	engineMap := map[string][]string{
		"all":       {"cromwell", "nextflow", "snakemake"},
		"cromwell":  {"cromwell"},
		"nextflow":  {"nextflow"},
		"snakemake": {"snakemake"},
		"wdl":       {"cromwell"},
		"nf":        {"nextflow"},
		"snake":     {"snakemake"},
	}
	var engines []string
	for _, part := range strings.Split(strings.ToLower(enginesArg), ",") {
		engines = append(engines, engineMap[strings.TrimSpace(part)]...)
	}
	return dedupe(engines)
}

// parsePlatformList parses the platform list argument.
func parsePlatformList(platformsArg string) []testing.SequencingPlatform {
	platformMap := map[string][]testing.SequencingPlatform{
		"all":      testing.AllPlatforms,
		"major":    {testing.IlluminaNovaSeq6000, testing.PacBioRevio, testing.ONTPromethION24},
		"illumina": {testing.IlluminaNovaSeq6000, testing.IlluminaNextSeq2000, testing.IlluminaMiSeq},
		"longread": {testing.PacBioRevio, testing.ONTPromethION24},
		"novaseq":  {testing.IlluminaNovaSeq6000},
		"nextseq":  {testing.IlluminaNextSeq2000},
		"miseq":    {testing.IlluminaMiSeq},
		"pacbio":   {testing.PacBioRevio},
		"ont":      {testing.ONTPromethION24},
		"dnbseq":   {testing.MGI_DNBSEQ_T20x2},
		"aviti":    {testing.ElementAVITI},
	}
	var platforms []testing.SequencingPlatform
	for _, part := range strings.Split(strings.ToLower(platformsArg), ",") {
		platforms = append(platforms, platformMap[strings.TrimSpace(part)]...)
	}
	return dedupe(platforms)
}

// parseWorkflowList parses the workflow list argument. Unknown names are
// kept as custom workflow names.
func parseWorkflowList(workflowsArg string) []string {
	workflowMap := map[string][]string{
		"all":        {"basic", "gatk", "assembly", "variants", "qc", "population"},
		"basic":      {"basic"},
		"gatk":       {"gatk_germline", "gatk_somatic"},
		"assembly":   {"pacbio_assembly", "ont_assembly"},
		"variants":   {"variant_calling", "variant_annotation"},
		"qc":         {"quality_control", "multiqc"},
		"population": {"population_genomics", "admixture"},
		"nfcore":     {"nfcore_sarek", "nfcore_nanoseq", "nfcore_hlatyping"},
	}
	var workflows []string
	for _, part := range strings.Split(strings.ToLower(workflowsArg), ",") {
		part = strings.TrimSpace(part)
		if list, ok := workflowMap[part]; ok {
			workflows = append(workflows, list...)
		} else if part != "" {
			workflows = append(workflows, part)
		}
	}
	return dedupe(workflows)
}

// Run statuses.
const (
	statusPassed  = "passed"
	statusFailed  = "failed"
	statusSkipped = "skipped"
	statusPlanned = "planned"
)

// WorkflowTestResult holds the outcome of one engine run.
type WorkflowTestResult struct {
	Engine   string                     `json:"engine"`
	Platform testing.SequencingPlatform `json:"platform"`
	Workflow string                     `json:"workflow"`
	Status   string                     `json:"status"`
	// Reason explains a skip, a plan note or a failure.
	Reason string `json:"reason,omitempty"`
	// WallSeconds is the time from generating the reads to the end of the
	// engine run, including engine start-up.
	WallSeconds float64 `json:"wall_seconds"`
	// ExitStatus is the engine's exit status; nil when it was not run or did
	// not exit normally (timeout or signal).
	ExitStatus *int `json:"exit_status,omitempty"`
	// Reads is the number of reads written for the workflow.
	Reads int `json:"reads"`
	// ReportedReads are the counts the workflow wrote, one per output.
	ReportedReads []string `json:"reported_reads,omitempty"`
	// DataIntegrity is the fraction of count outputs equal to Reads.
	DataIntegrity float64  `json:"data_integrity"`
	OutputFiles   []string `json:"output_files,omitempty"`
}

// OutcomeCounts tallies run statuses.
type OutcomeCounts struct {
	Passed, Failed, Skipped, Planned int
}

// WorkflowTestResults holds every run.
type WorkflowTestResults struct {
	Engines     []string                     `json:"engines"`
	Platforms   []testing.SequencingPlatform `json:"platforms"`
	Workflows   []string                     `json:"workflows"`
	DryRun      bool                         `json:"dry_run"`
	StartedAt   time.Time                    `json:"started_at"`
	WallSeconds float64                      `json:"wall_seconds"`
	Runs        []*WorkflowTestResult        `json:"runs"`
}

// Counts tallies the runs by status.
func (r *WorkflowTestResults) Counts() OutcomeCounts {
	var c OutcomeCounts
	for _, run := range r.Runs {
		switch run.Status {
		case statusPassed:
			c.Passed++
		case statusFailed:
			c.Failed++
		case statusSkipped:
			c.Skipped++
		case statusPlanned:
			c.Planned++
		}
	}
	return c
}

// runWorkflowTests runs (or, with -dry-run, plans) every engine × platform ×
// workflow combination once.
func runWorkflowTests(framework *testing.E2EFramework, engines []string, platforms []testing.SequencingPlatform, workflows []string, args *CommandLineArgs) *WorkflowTestResults {
	start := time.Now()
	results := &WorkflowTestResults{Engines: engines, Platforms: platforms, Workflows: workflows, DryRun: args.dryRun, StartedAt: start}
	for _, engine := range engines {
		for _, platform := range platforms {
			for _, workflow := range workflows {
				var r *WorkflowTestResult
				if args.dryRun {
					r = planWorkflowTest(framework, engine, platform, workflow, args)
				} else {
					r = runWorkflowTest(framework, engine, platform, workflow, args)
				}
				if args.verbose {
					log.Printf("%s %s %s: %s %s", engine, platform, workflow, r.Status, r.Reason)
				}
				results.Runs = append(results.Runs, r)
			}
		}
	}
	results.WallSeconds = time.Since(start).Seconds()
	return results
}

// planWorkflowTest records what runWorkflowTest would do, without running
// anything: skipped where it would skip, planned otherwise.
func planWorkflowTest(framework *testing.E2EFramework, engine string, platform testing.SequencingPlatform, workflow string, args *CommandLineArgs) *WorkflowTestResult {
	r := &WorkflowTestResult{Engine: engine, Platform: platform, Workflow: workflow, Status: statusPlanned}
	if !framework.Modeled(platform) {
		r.Status, r.Reason = statusSkipped, fmt.Sprintf("%v: %s", testing.ErrPlatformNotModeled, platform)
		return r
	}
	if _, _, err := engineCommand(engine); err != nil {
		if errors.Is(err, errEngineUnavailable) && args.skipUnavailable {
			r.Status = statusSkipped
		}
		r.Reason = err.Error()
	}
	return r
}

// runWorkflowTest generates reads for platform, runs workflow on engine and
// checks the count it reports.
func runWorkflowTest(framework *testing.E2EFramework, engine string, platform testing.SequencingPlatform, workflow string, args *CommandLineArgs) *WorkflowTestResult {
	start := time.Now()
	r := &WorkflowTestResult{Engine: engine, Platform: platform, Workflow: workflow, Status: statusFailed}
	defer func() { r.WallSeconds = time.Since(start).Seconds() }()

	if _, _, err := engineCommand(engine); err != nil {
		if errors.Is(err, errEngineUnavailable) && args.skipUnavailable {
			r.Status = statusSkipped
		}
		r.Reason = err.Error()
		return r
	}

	data, err := framework.GenerateSyntheticData(framework.CreatePlatformConfig(platform, "CEU"))
	if errors.Is(err, testing.ErrPlatformNotModeled) {
		r.Status, r.Reason = statusSkipped, err.Error()
		return r
	}
	if err != nil {
		r.Reason = fmt.Sprintf("generate reads: %v", err)
		return r
	}
	r.Reads = len(data.Sequences)

	workDir, err := os.MkdirTemp(args.outputDir, fmt.Sprintf("%s_%s_%s", engine, platform, workflow))
	if err != nil {
		r.Reason = fmt.Sprintf("create workspace: %v", err)
		return r
	}
	defer os.RemoveAll(workDir)

	var outputs []string
	switch engine {
	case "cromwell":
		outputs, err = executeCromwellWorkflow(workDir, data, platform, workflow, args)
	case "nextflow":
		outputs, err = executeNextflowWorkflow(workDir, data, platform, workflow, args)
	case "snakemake":
		outputs, err = executeSnakemakeWorkflow(workDir, data, platform, workflow, args)
	}
	var runErr *engineRunError
	switch {
	case errors.As(err, &runErr):
		r.ExitStatus = runErr.exitStatus
	case err == nil || errors.Is(err, errNoOutput):
		zero := 0
		r.ExitStatus = &zero
	}
	if err != nil {
		r.Reason = err.Error()
		return r
	}

	for _, f := range outputs {
		r.ReportedReads = append(r.ReportedReads, readOutput(f))
	}
	r.OutputFiles = outputs
	r.DataIntegrity = calculateDataIntegrity(data, outputs)
	if r.DataIntegrity < 1.0 {
		r.Reason = fmt.Sprintf("reported read count %v does not match the %d reads written", r.ReportedReads, r.Reads)
		return r
	}
	r.Status = statusPassed
	return r
}

// readOutput returns a count output's trimmed content, or a note when it
// cannot be read.
func readOutput(path string) string {
	b, err := os.ReadFile(path) // #nosec G304 -- engine output inside the test workspace
	if err != nil {
		return "unreadable"
	}
	return strings.TrimSpace(string(b))
}

// Workflow execution functions
//
// Each execute function writes the synthetic reads as FASTQ, writes a minimal
// read-counting pipeline for the engine, runs the engine, and returns the
// pipeline's output files. An error is returned if the engine is unavailable,
// exits non-zero, or produces no output.

// errEngineUnavailable reports that a workflow engine is not installed.
var errEngineUnavailable = errors.New("workflow engine not available")

// errNoOutput reports that an engine exited 0 but wrote no count output.
var errNoOutput = errors.New("workflow produced no read count output")

// engineRunError is returned when an engine run fails.
type engineRunError struct {
	engine string
	// exitStatus is the exit status, or nil when the process did not exit
	// normally (timeout, signal, or it could not be started).
	exitStatus *int
	err        error
	output     string
}

func (e *engineRunError) Error() string {
	return fmt.Sprintf("%s run failed: %v\n%s", e.engine, e.err, e.output)
}

func (e *engineRunError) Unwrap() error { return e.err }

// lookPath locates executables; tests replace it to simulate missing engines.
var lookPath = exec.LookPath

// inputFASTQ is the file name of the reads handed to each workflow.
const inputFASTQ = "reads.fastq"

// readCountFile is the output each workflow writes: the number of FASTQ records.
const readCountFile = "read_count.txt"

// engineCommand returns the executable and leading arguments that invoke
// engine. Cromwell is found as a "cromwell" executable on PATH, or as
// "java -jar $CROMWELL_JAR" when CROMWELL_JAR names an existing file.
func engineCommand(engine string) (string, []string, error) {
	switch engine {
	case "cromwell":
		if path, err := lookPath("cromwell"); err == nil {
			return path, nil, nil
		}
		if jar := os.Getenv("CROMWELL_JAR"); jar != "" {
			if _, err := os.Stat(jar); err != nil {
				return "", nil, fmt.Errorf("%w: cromwell not on PATH and CROMWELL_JAR=%s: %v", errEngineUnavailable, jar, err)
			}
			java, err := lookPath("java")
			if err != nil {
				return "", nil, fmt.Errorf("%w: CROMWELL_JAR is set but java is not on PATH", errEngineUnavailable)
			}
			return java, []string{"-jar", jar}, nil
		}
		return "", nil, fmt.Errorf("%w: cromwell not on PATH and CROMWELL_JAR not set", errEngineUnavailable)
	case "nextflow", "snakemake":
		path, err := lookPath(engine)
		if err != nil {
			return "", nil, fmt.Errorf("%w: %s not on PATH", errEngineUnavailable, engine)
		}
		return path, nil, nil
	default:
		return "", nil, fmt.Errorf("unsupported engine: %s", engine)
	}
}

// isWorkflowEngineAvailable reports whether engine can be invoked.
func isWorkflowEngineAvailable(engine string) bool {
	_, _, err := engineCommand(engine)
	return err == nil
}

// runEngine runs engine with extra arguments in workDir, honouring the
// -timeout flag. On failure it returns an *engineRunError carrying the exit
// status and the tail of the engine's output.
func runEngine(engine, workDir string, args *CommandLineArgs, env []string, extra ...string) error {
	program, lead, err := engineCommand(engine)
	if err != nil {
		return err
	}
	ctx := context.Background()
	if args != nil && args.timeout > 0 {
		var cancel context.CancelFunc
		ctx, cancel = context.WithTimeout(ctx, time.Duration(args.timeout)*time.Second)
		defer cancel()
	}
	cmd := exec.CommandContext(ctx, program, append(lead, extra...)...) // #nosec G204 -- program resolved via LookPath, arguments built internally
	cmd.Dir = workDir
	cmd.Env = append(os.Environ(), env...)
	out, err := cmd.CombinedOutput()
	if err != nil {
		runErr := &engineRunError{engine: engine, err: err, output: tail(string(out), 40)}
		if ctx.Err() != nil {
			runErr.err = ctx.Err()
		} else if cmd.ProcessState != nil && cmd.ProcessState.Exited() {
			code := cmd.ProcessState.ExitCode()
			runErr.exitStatus = &code
		}
		return runErr
	}
	if args != nil && args.verbose {
		log.Printf("     %s output:\n%s", engine, tail(string(out), 20))
	}
	return nil
}

// tail returns the last n lines of s.
func tail(s string, n int) string {
	lines := strings.Split(strings.TrimRight(s, "\n"), "\n")
	if len(lines) > n {
		lines = lines[len(lines)-n:]
	}
	return strings.Join(lines, "\n")
}

// writeInputFASTQ writes the synthetic reads to workDir and returns the absolute path.
func writeInputFASTQ(workDir string, data *testing.SyntheticData) (string, error) {
	path, err := filepath.Abs(filepath.Join(workDir, inputFASTQ))
	if err != nil {
		return "", err
	}
	f, err := os.Create(path) // #nosec G304 -- path inside the test workspace
	if err != nil {
		return "", err
	}
	w := bufio.NewWriter(f)
	for i, seq := range data.Sequences {
		qual := []byte(strings.Repeat("I", len(seq)))
		if i < len(data.QualityScores) && len(data.QualityScores[i]) == len(seq) {
			qual = data.QualityScores[i]
		}
		fmt.Fprintf(w, "@read_%d\n%s\n+\n%s\n", i, seq, qual)
	}
	if err := w.Flush(); err != nil {
		f.Close()
		return "", err
	}
	return path, f.Close()
}

// requireOutputs returns paths unchanged if every path exists, or an error
// wrapping errNoOutput naming the first missing one.
func requireOutputs(engine string, paths []string) ([]string, error) {
	if len(paths) == 0 {
		return nil, fmt.Errorf("%w: %s produced no %s", errNoOutput, engine, readCountFile)
	}
	for _, p := range paths {
		if _, err := os.Stat(p); err != nil {
			return nil, fmt.Errorf("%w: %s did not produce %s: %v", errNoOutput, engine, p, err)
		}
	}
	return paths, nil
}

func executeCromwellWorkflow(workDir string, syntheticData *testing.SyntheticData, platform testing.SequencingPlatform, workflow string, args *CommandLineArgs) ([]string, error) {
	fastq, err := writeInputFASTQ(workDir, syntheticData)
	if err != nil {
		return nil, fmt.Errorf("write input FASTQ: %w", err)
	}
	wdlFile := filepath.Join(workDir, "test.wdl")
	if err := os.WriteFile(wdlFile, []byte(createBasicWDL(workflow, string(platform))), 0600); err != nil {
		return nil, fmt.Errorf("write WDL: %w", err)
	}
	inputsFile := filepath.Join(workDir, "inputs.json")
	if err := os.WriteFile(inputsFile, []byte(createWDLInputs(workflow, string(platform), fastq)), 0600); err != nil {
		return nil, fmt.Errorf("write WDL inputs: %w", err)
	}

	if err := runEngine("cromwell", workDir, args, nil, "run", wdlFile, "-i", inputsFile); err != nil {
		return nil, err
	}

	// Cromwell's local backend writes task outputs under
	// cromwell-executions/<workflow>/<run id>/call-<task>/execution/.
	matches, _ := filepath.Glob(filepath.Join(workDir, "cromwell-executions", "*", "*", "call-CountReads", "execution", readCountFile))
	return requireOutputs("cromwell", matches)
}

func executeNextflowWorkflow(workDir string, syntheticData *testing.SyntheticData, platform testing.SequencingPlatform, workflow string, args *CommandLineArgs) ([]string, error) {
	fastq, err := writeInputFASTQ(workDir, syntheticData)
	if err != nil {
		return nil, fmt.Errorf("write input FASTQ: %w", err)
	}
	nfFile := filepath.Join(workDir, "main.nf")
	if err := os.WriteFile(nfFile, []byte(createBasicNextflow(workflow, string(platform))), 0600); err != nil {
		return nil, fmt.Errorf("write Nextflow script: %w", err)
	}
	configFile := filepath.Join(workDir, "nextflow.config")
	if err := os.WriteFile(configFile, []byte("process.executor = 'local'\n"), 0600); err != nil {
		return nil, fmt.Errorf("write Nextflow config: %w", err)
	}

	outDir := filepath.Join(workDir, "output")
	if err := runEngine("nextflow", workDir, args, []string{"NXF_OFFLINE=true", "NXF_ANSI_LOG=false"},
		"run", nfFile, "-c", configFile, "--fastq", fastq, "--outdir", outDir); err != nil {
		return nil, err
	}
	return requireOutputs("nextflow", []string{filepath.Join(outDir, readCountFile)})
}

func executeSnakemakeWorkflow(workDir string, syntheticData *testing.SyntheticData, platform testing.SequencingPlatform, workflow string, args *CommandLineArgs) ([]string, error) {
	if _, err := writeInputFASTQ(workDir, syntheticData); err != nil {
		return nil, fmt.Errorf("write input FASTQ: %w", err)
	}
	snakeFile := filepath.Join(workDir, "Snakefile")
	if err := os.WriteFile(snakeFile, []byte(createBasicSnakemake(workflow, string(platform))), 0600); err != nil {
		return nil, fmt.Errorf("write Snakefile: %w", err)
	}

	if err := runEngine("snakemake", workDir, args, nil, "--cores", "1", "--snakefile", snakeFile, "--directory", workDir); err != nil {
		return nil, err
	}
	return requireOutputs("snakemake", []string{filepath.Join(workDir, readCountFile)})
}

// wdlIdentifier turns a workflow name into a valid WDL identifier.
func wdlIdentifier(workflow string) string {
	var b strings.Builder
	for i, r := range workflow {
		switch {
		case unicode.IsLetter(r) && r < unicode.MaxASCII:
			if i == 0 {
				r = unicode.ToUpper(r)
			}
			b.WriteRune(r)
		case unicode.IsDigit(r) && r < unicode.MaxASCII:
			if i == 0 {
				b.WriteString("W")
			}
			b.WriteRune(r)
		default:
			b.WriteRune('_')
		}
	}
	if b.Len() == 0 {
		return "Workflow"
	}
	if s := b.String(); s[0] == '_' {
		return "W" + s
	}
	return b.String()
}

// createBasicWDL returns a WDL workflow that counts the records in a FASTQ.
func createBasicWDL(workflow, platform string) string {
	return fmt.Sprintf(`version 1.0

# %s workflow on %s platform: counts FASTQ records.
workflow %s {
    input {
        File fastq
        String sample_name = "synthetic_%s"
    }

    call CountReads {
        input: fastq = fastq
    }

    output {
        File read_count = CountReads.read_count
    }
}

task CountReads {
    input {
        File fastq
    }

    command <<<
        awk 'END { print NR / 4 }' ~{fastq} > %s
    >>>

    output {
        File read_count = "%s"
    }
}
`, workflow, platform, wdlIdentifier(workflow), platform, readCountFile, readCountFile)
}

// createWDLInputs returns the Cromwell inputs JSON for createBasicWDL.
func createWDLInputs(workflow, platform, fastqPath string) string {
	name := wdlIdentifier(workflow)
	return fmt.Sprintf(`{
    "%s.fastq": %s,
    "%s.sample_name": "synthetic_%s"
}`, name, strconv.Quote(fastqPath), name, platform)
}

// createBasicNextflow returns a Nextflow DSL2 script that counts the records
// in params.fastq and publishes the count to params.outdir.
func createBasicNextflow(workflow, platform string) string {
	return fmt.Sprintf(`#!/usr/bin/env nextflow
// %s workflow on %s platform: counts FASTQ records.

params.fastq = null
params.outdir = 'output'

process countReads {
    publishDir params.outdir, mode: 'copy'

    input:
    path fastq

    output:
    path "%s"

    script:
    """
    awk 'END { print NR / 4 }' ${fastq} > %s
    """
}

workflow {
    countReads(Channel.fromPath(params.fastq, checkIfExists: true))
}
`, workflow, platform, readCountFile, readCountFile)
}

// createBasicSnakemake returns a Snakefile that counts the records in reads.fastq.
func createBasicSnakemake(workflow, platform string) string {
	return fmt.Sprintf(`# %s workflow on %s platform: counts FASTQ records.
rule all:
    input:
        "%s"

rule count_reads:
    input:
        "%s"
    output:
        "%s"
    shell:
        "awk 'END {{ print NR / 4 }}' {input} > {output}"
`, workflow, platform, readCountFile, inputFASTQ, readCountFile)
}

// calculateDataIntegrity returns the fraction of read-count output files whose
// count equals the number of reads written. Missing, unreadable or
// unparseable files count as mismatches; no files gives 0.
func calculateDataIntegrity(syntheticData *testing.SyntheticData, outputFiles []string) float64 {
	if syntheticData == nil || len(outputFiles) == 0 {
		return 0.0
	}
	want := len(syntheticData.Sequences)
	matched := 0
	for _, file := range outputFiles {
		content, err := os.ReadFile(file) // #nosec G304 -- engine output inside the test workspace
		if err != nil {
			continue
		}
		got, err := strconv.ParseFloat(strings.TrimSpace(string(content)), 64)
		if err == nil && got == float64(want) {
			matched++
		}
	}
	return float64(matched) / float64(len(outputFiles))
}

// printTestSummary prints one row per run and the outcome counts.
func printTestSummary(w io.Writer, results *WorkflowTestResults, args *CommandLineArgs) {
	rule := strings.Repeat("-", 110)
	mode := "executed"
	if results.DryRun {
		mode = "dry run, nothing executed"
	}
	fmt.Fprintf(w, "Workflow runs: %d engine(s) × %d platform(s) × %d workflow(s) (%s), %.1fs\n",
		len(results.Engines), len(results.Platforms), len(results.Workflows), mode, results.WallSeconds)
	fmt.Fprintln(w, rule)
	fmt.Fprintf(w, "%-10s %-22s %-20s %-8s %8s %5s %6s %9s %9s\n",
		"Engine", "Platform", "Workflow", "Status", "Wall(s)", "Exit", "Reads", "Reported", "Integrity")
	fmt.Fprintln(w, rule)
	for _, r := range results.Runs {
		exit, wall, reads, reported, integrity := "-", "-", "-", "-", "-"
		if r.ExitStatus != nil {
			exit = strconv.Itoa(*r.ExitStatus)
		}
		if r.Status == statusPassed || r.Status == statusFailed {
			wall = fmt.Sprintf("%.2f", r.WallSeconds)
		}
		if r.Reads > 0 {
			reads = strconv.Itoa(r.Reads)
		}
		if len(r.ReportedReads) > 0 {
			reported = strings.Join(r.ReportedReads, ",")
		}
		if len(r.OutputFiles) > 0 {
			integrity = fmt.Sprintf("%.2f", r.DataIntegrity)
		}
		fmt.Fprintf(w, "%-10s %-22s %-20s %-8s %8s %5s %6s %9s %9s\n",
			r.Engine, r.Platform, r.Workflow, strings.ToUpper(r.Status), wall, exit, reads, reported, integrity)
		if r.Reason != "" && r.Status != statusPassed {
			fmt.Fprintf(w, "    %s\n", firstLine(r.Reason))
		}
	}
	fmt.Fprintln(w, rule)

	c := results.Counts()
	if results.DryRun {
		fmt.Fprintf(w, "Dry run: %d planned, %d would be skipped, %d would fail before running; nothing was executed.\n",
			c.Planned, c.Skipped, c.Failed)
		return
	}
	fmt.Fprintf(w, "Passed %d, failed %d, skipped %d.\n", c.Passed, c.Failed, c.Skipped)
	if c.Passed+c.Failed == 0 {
		fmt.Fprintln(w, "No workflow was executed: every run was skipped (install cromwell, nextflow or snakemake to run them).")
	}
}

// firstLine returns s up to its first newline.
func firstLine(s string) string {
	if i := strings.IndexByte(s, '\n'); i >= 0 {
		return s[:i]
	}
	return s
}

// writeReports writes workflow_results.json (every run) and
// workflow_results.csv (one row per run) to dir and returns their paths.
func writeReports(results *WorkflowTestResults, dir string) ([]string, error) {
	if err := os.MkdirAll(dir, 0750); err != nil {
		return nil, err
	}
	jsonPath := filepath.Join(dir, "workflow_results.json")
	b, err := json.MarshalIndent(results, "", "  ")
	if err != nil {
		return nil, err
	}
	if err := os.WriteFile(jsonPath, append(b, '\n'), 0600); err != nil {
		return nil, err
	}

	csvPath := filepath.Join(dir, "workflow_results.csv")
	fh, err := os.Create(csvPath) // #nosec G304 -- caller-chosen output directory
	if err != nil {
		return nil, err
	}
	cw := csv.NewWriter(fh)
	_ = cw.Write([]string{"engine", "platform", "workflow", "status", "wall_seconds", "exit_status", "reads",
		"reported_reads", "data_integrity", "reason"})
	for _, r := range results.Runs {
		exit := ""
		if r.ExitStatus != nil {
			exit = strconv.Itoa(*r.ExitStatus)
		}
		_ = cw.Write([]string{r.Engine, string(r.Platform), r.Workflow, r.Status,
			strconv.FormatFloat(r.WallSeconds, 'f', 3, 64), exit, strconv.Itoa(r.Reads),
			strings.Join(r.ReportedReads, ";"), strconv.FormatFloat(r.DataIntegrity, 'f', -1, 64), firstLine(r.Reason)})
	}
	cw.Flush()
	if err := cw.Error(); err != nil {
		fh.Close()
		return nil, err
	}
	if err := fh.Close(); err != nil {
		return nil, err
	}
	return []string{jsonPath, csvPath}, nil
}
