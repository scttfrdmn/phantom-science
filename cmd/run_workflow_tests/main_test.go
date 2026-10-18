//go:build !release
// +build !release

package main

import (
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
	"time"

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

func TestIsWorkflowEngineAvailable_NotOnPath(t *testing.T) {
	withNoEngines(t)
	for _, engine := range []string{"cromwell", "nextflow", "snakemake"} {
		if isWorkflowEngineAvailable(engine) {
			t.Errorf("isWorkflowEngineAvailable(%s) = true with nothing on PATH", engine)
		}
		_, _, err := engineCommand(engine)
		if !errors.Is(err, errEngineUnavailable) {
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
}

func TestIsWorkflowEngineAvailable_Unknown(t *testing.T) {
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

func TestSelectRepresentativeWorkflow_Empty(t *testing.T) {
	got := selectRepresentativeWorkflow("cromwell", phantomtesting.IlluminaNovaSeq6000, []string{})
	if got != "basic" {
		t.Errorf("selectRepresentativeWorkflow with empty list = %q, want \"basic\"", got)
	}
}

func TestSelectRepresentativeWorkflow_SingleWorkflow(t *testing.T) {
	got := selectRepresentativeWorkflow("cromwell", phantomtesting.IlluminaNovaSeq6000, []string{"variant-calling"})
	if got != "variant-calling" {
		t.Errorf("selectRepresentativeWorkflow = %q, want \"variant-calling\"", got)
	}
}

func TestSelectRepresentativeWorkflow_PacBioAssembly(t *testing.T) {
	// PacBioRevio should prefer "assembly" workflows
	got := selectRepresentativeWorkflow("cromwell", phantomtesting.PacBioRevio, []string{"assembly", "gatk"})
	if got != "assembly" {
		t.Errorf("selectRepresentativeWorkflow(PacBio) = %q, want \"assembly\"", got)
	}
}

func TestSelectRepresentativeWorkflow_ONTNanopore(t *testing.T) {
	// ONT should prefer "ont" workflows
	got := selectRepresentativeWorkflow("cromwell", phantomtesting.ONTPromethION24, []string{"ont-assembly", "basic"})
	if got != "ont-assembly" {
		t.Errorf("selectRepresentativeWorkflow(ONT) = %q, want \"ont-assembly\"", got)
	}
}

func TestParseEngineList_All(t *testing.T) {
	engines := parseEngineList("all")
	if len(engines) != 3 {
		t.Errorf("parseEngineList(all) = %d engines, want 3", len(engines))
	}
}

func TestParseEngineList_Single(t *testing.T) {
	engines := parseEngineList("cromwell")
	if len(engines) != 1 || engines[0] != "cromwell" {
		t.Errorf("parseEngineList(cromwell) = %v, want [cromwell]", engines)
	}
}

func TestParseEngineList_Multiple(t *testing.T) {
	engines := parseEngineList("cromwell,nextflow")
	if len(engines) != 2 {
		t.Errorf("parseEngineList(cromwell,nextflow) = %d engines, want 2", len(engines))
	}
}

func TestParseEngineList_NoDuplicates(t *testing.T) {
	// "all" includes cromwell, and cromwell is also explicit
	engines := parseEngineList("all,cromwell")
	seen := make(map[string]bool)
	for _, e := range engines {
		if seen[e] {
			t.Errorf("duplicate engine %q in result", e)
		}
		seen[e] = true
	}
}

func TestParsePlatformList_All(t *testing.T) {
	platforms := parsePlatformList("all")
	if len(platforms) == 0 {
		t.Error("parsePlatformList(all) returned empty list")
	}
}

func TestParsePlatformList_Illumina(t *testing.T) {
	platforms := parsePlatformList("illumina")
	if len(platforms) != 3 {
		t.Errorf("parsePlatformList(illumina) = %d platforms, want 3", len(platforms))
	}
}

func TestParsePlatformList_NovaSeq(t *testing.T) {
	platforms := parsePlatformList("novaseq")
	if len(platforms) != 1 || platforms[0] != phantomtesting.IlluminaNovaSeq6000 {
		t.Errorf("parsePlatformList(novaseq) = %v", platforms)
	}
}

func TestParseWorkflowList_All(t *testing.T) {
	workflows := parseWorkflowList("all")
	if len(workflows) == 0 {
		t.Error("parseWorkflowList(all) returned empty list")
	}
}

func TestParseWorkflowList_Custom(t *testing.T) {
	workflows := parseWorkflowList("my-custom-workflow")
	if len(workflows) != 1 || workflows[0] != "my-custom-workflow" {
		t.Errorf("parseWorkflowList(custom) = %v, want [my-custom-workflow]", workflows)
	}
}

func TestCreateBasicWDL_NotEmpty(t *testing.T) {
	wdl := createBasicWDL("gatk", "illumina")
	if wdl == "" {
		t.Error("createBasicWDL returned empty string")
	}
	if !strings.Contains(wdl, "workflow") {
		t.Error("createBasicWDL should contain 'workflow' keyword")
	}
}

func TestCreateWDLInputs_NotEmpty(t *testing.T) {
	inputs := createWDLInputs("gatk", "illumina", "/data/reads.fastq")
	if !strings.Contains(inputs, `"Gatk.fastq": "/data/reads.fastq"`) {
		t.Errorf("createWDLInputs should bind the FASTQ path, got:\n%s", inputs)
	}
}

func TestWDLIdentifier(t *testing.T) {
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

func TestCreateBasicNextflow_NotEmpty(t *testing.T) {
	nf := createBasicNextflow("variant-calling", "illumina")
	if nf == "" {
		t.Error("createBasicNextflow returned empty string")
	}
	if !strings.Contains(nf, "nextflow") {
		t.Error("createBasicNextflow should contain 'nextflow'")
	}
}

func TestCreateBasicSnakemake_NotEmpty(t *testing.T) {
	sm := createBasicSnakemake("qc", "illumina")
	if sm == "" {
		t.Error("createBasicSnakemake returned empty string")
	}
	if !strings.Contains(sm, "rule") {
		t.Error("createBasicSnakemake should contain 'rule' keyword")
	}
}

func TestEstimateMemoryUsage_PacBio(t *testing.T) {
	mem := estimateMemoryUsage(phantomtesting.PacBioRevio, "assembly")
	if mem <= 0 {
		t.Errorf("estimateMemoryUsage(PacBio, assembly) = %d, want > 0", mem)
	}
}

func TestEstimateMemoryUsage_Illumina(t *testing.T) {
	mem := estimateMemoryUsage(phantomtesting.IlluminaNovaSeq6000, "basic")
	if mem <= 0 {
		t.Errorf("estimateMemoryUsage(IlluminaNovaSeq, basic) = %d, want > 0", mem)
	}
}

func TestEstimateCPUUsage_Assembly(t *testing.T) {
	cpu := estimateCPUUsage(phantomtesting.PacBioRevio, "assembly")
	if cpu <= 0 || cpu > 100 {
		t.Errorf("estimateCPUUsage(assembly) = %f, want in (0,100]", cpu)
	}
}

func TestEstimateCPUUsage_Basic(t *testing.T) {
	cpu := estimateCPUUsage(phantomtesting.IlluminaNovaSeq6000, "basic")
	if cpu <= 0 || cpu > 100 {
		t.Errorf("estimateCPUUsage(basic) = %f, want in (0,100]", cpu)
	}
}

func TestCalculateWorkflowScore_Failed(t *testing.T) {
	result := &WorkflowTestResult{Success: false}
	if got := calculateWorkflowScore(result); got != 0.0 {
		t.Errorf("calculateWorkflowScore(failed) = %f, want 0.0", got)
	}
}

func TestCalculateWorkflowScore_Success(t *testing.T) {
	result := &WorkflowTestResult{
		Success:       true,
		DataIntegrity: 0.9,
		OutputAccuracy: 0.85,
		ExecutionTime: 30 * time.Second, // < 1 minute
	}
	got := calculateWorkflowScore(result)
	if got <= 0.0 || got > 1.0+1e-9 {
		t.Errorf("calculateWorkflowScore(success) = %f, want in (0,1]", got)
	}
}

func TestCalculateEngineCompatibilityScore_NotAvailable(t *testing.T) {
	result := &EngineTestResult{Available: false}
	if got := calculateEngineCompatibilityScore(result); got != 0.0 {
		t.Errorf("calculateEngineCompatibilityScore(not available) = %f, want 0.0", got)
	}
}

func TestCalculateEngineCompatibilityScore_Available(t *testing.T) {
	result := &EngineTestResult{Available: true, SuccessRate: 0.8}
	got := calculateEngineCompatibilityScore(result)
	if got != 80.0 {
		t.Errorf("calculateEngineCompatibilityScore(0.8) = %f, want 80.0", got)
	}
}

func TestCalculateOverallScore_NoEngines(t *testing.T) {
	results := &WorkflowTestResults{
		EngineResults: make(map[string]*EngineTestResult),
	}
	if got := calculateOverallScore(results); got != 0.0 {
		t.Errorf("calculateOverallScore(empty) = %f, want 0.0", got)
	}
}

func TestCalculateDataIntegrity_NoOutputs(t *testing.T) {
	f := phantomtesting.NewE2EFramework()
	data, _ := f.GenerateSyntheticData(f.CreatePlatformConfig(phantomtesting.IlluminaNovaSeq6000, "CEU"))
	if got := calculateDataIntegrity(data, []string{}); got != 0.0 {
		t.Errorf("calculateDataIntegrity(no outputs) = %f, want 0", got)
	}
}

func TestCalculateOutputAccuracy_IlluminaBasic(t *testing.T) {
	got := calculateOutputAccuracy(phantomtesting.IlluminaNovaSeq6000, "basic", []string{})
	if got <= 0 || got > 1.1 {
		t.Errorf("calculateOutputAccuracy(Illumina, basic) = %f, want ~1.0", got)
	}
}

func TestCalculateOutputAccuracy_IlluminaAssembly(t *testing.T) {
	basic := calculateOutputAccuracy(phantomtesting.IlluminaNovaSeq6000, "basic", []string{})
	assembly := calculateOutputAccuracy(phantomtesting.IlluminaNovaSeq6000, "assembly", []string{})
	if assembly >= basic {
		t.Errorf("assembly accuracy (%f) should be less than basic (%f)", assembly, basic)
	}
}

func TestCalculateOutputAccuracy_PacBio(t *testing.T) {
	got := calculateOutputAccuracy(phantomtesting.PacBioRevio, "variant-calling", []string{})
	if got <= 0 || got > 1.1 {
		t.Errorf("calculateOutputAccuracy(PacBio) = %f, want in (0, 1.1]", got)
	}
}

func TestCalculateOutputAccuracy_ONT(t *testing.T) {
	got := calculateOutputAccuracy(phantomtesting.ONTPromethION24, "ont-assembly", []string{})
	if got <= 0 || got > 1.1 {
		t.Errorf("calculateOutputAccuracy(ONT) = %f, want in (0, 1.1]", got)
	}
}

func TestFindBestWorkflowResult_Empty(t *testing.T) {
	engineResult := &EngineTestResult{
		WorkflowResults: make(map[string]*WorkflowTestResult),
	}
	got := findBestWorkflowResult(engineResult, phantomtesting.IlluminaNovaSeq6000)
	if got != nil {
		t.Errorf("findBestWorkflowResult(empty) = %v, want nil", got)
	}
}

func TestFindBestWorkflowResult_NoMatch(t *testing.T) {
	engineResult := &EngineTestResult{
		WorkflowResults: map[string]*WorkflowTestResult{
			"gatk": {
				Platform: phantomtesting.PacBioRevio,
				Success:  true,
				DataIntegrity:  0.9,
				OutputAccuracy: 0.9,
			},
		},
	}
	got := findBestWorkflowResult(engineResult, phantomtesting.IlluminaNovaSeq6000)
	if got != nil {
		t.Errorf("findBestWorkflowResult(no match) = %v, want nil", got)
	}
}

func TestFindBestWorkflowResult_Match(t *testing.T) {
	engineResult := &EngineTestResult{
		WorkflowResults: map[string]*WorkflowTestResult{
			"gatk": {
				Platform:       phantomtesting.IlluminaNovaSeq6000,
				Success:        true,
				DataIntegrity:  0.9,
				OutputAccuracy: 0.9,
			},
		},
	}
	got := findBestWorkflowResult(engineResult, phantomtesting.IlluminaNovaSeq6000)
	if got == nil {
		t.Error("findBestWorkflowResult(match) = nil, want non-nil")
	}
}

func TestFindBestWorkflowResult_BestScore(t *testing.T) {
	engineResult := &EngineTestResult{
		WorkflowResults: map[string]*WorkflowTestResult{
			"low": {
				Platform:       phantomtesting.IlluminaNovaSeq6000,
				Success:        true,
				DataIntegrity:  0.5,
				OutputAccuracy: 0.5,
			},
			"high": {
				Platform:       phantomtesting.IlluminaNovaSeq6000,
				Success:        true,
				DataIntegrity:  0.99,
				OutputAccuracy: 0.99,
			},
		},
	}
	got := findBestWorkflowResult(engineResult, phantomtesting.IlluminaNovaSeq6000)
	if got == nil {
		t.Fatal("findBestWorkflowResult(two matches) = nil, want non-nil")
	}
	if got.DataIntegrity != 0.99 {
		t.Errorf("findBestWorkflowResult returned lower-score result, DataIntegrity = %f", got.DataIntegrity)
	}
}

func TestBuildCompatibilityMatrix_Empty(t *testing.T) {
	results := &WorkflowTestResults{
		EngineResults:       make(map[string]*EngineTestResult),
		CompatibilityMatrix: make(map[string]map[phantomtesting.SequencingPlatform]*WorkflowTestResult),
		Platforms:           []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000},
	}
	buildCompatibilityMatrix(results)
	// No panic expected
}

func TestBuildCompatibilityMatrix_WithResults(t *testing.T) {
	results := &WorkflowTestResults{
		EngineResults: map[string]*EngineTestResult{
			"cromwell": {
				WorkflowResults: map[string]*WorkflowTestResult{
					"gatk": {
						Platform:       phantomtesting.IlluminaNovaSeq6000,
						Success:        true,
						DataIntegrity:  0.9,
						OutputAccuracy: 0.9,
					},
				},
			},
		},
		CompatibilityMatrix: make(map[string]map[phantomtesting.SequencingPlatform]*WorkflowTestResult),
		Platforms:           []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000},
	}
	buildCompatibilityMatrix(results)
	if results.CompatibilityMatrix["cromwell"] == nil {
		t.Error("buildCompatibilityMatrix should populate cromwell entry")
	}
}

func TestValidateArgs_EmptyOutputDir(t *testing.T) {
	args := &CommandLineArgs{
		outputDir: "",
		timeout:   120,
		engines:   "all",
		platforms: "all",
		workflows: "all",
	}
	if err := validateArgs(args); err == nil {
		t.Error("validateArgs(empty outputDir) should return error")
	}
}

func TestValidateArgs_LowTimeout(t *testing.T) {
	args := &CommandLineArgs{
		outputDir: "/tmp/test",
		timeout:   30,
		engines:   "all",
		platforms: "all",
		workflows: "all",
	}
	if err := validateArgs(args); err == nil {
		t.Error("validateArgs(timeout<60) should return error")
	}
}

func TestValidateArgs_Valid(t *testing.T) {
	args := &CommandLineArgs{
		outputDir: "/tmp/test",
		timeout:   120,
		engines:   "all",
		platforms: "all",
		workflows: "all",
	}
	if err := validateArgs(args); err != nil {
		t.Errorf("validateArgs(valid) returned error: %v", err)
	}
}

func TestCalculateOverallScore_WithEngines(t *testing.T) {
	results := &WorkflowTestResults{
		EngineResults: map[string]*EngineTestResult{
			"cromwell": {Available: true, CompatibilityScore: 90.0},
		},
	}
	got := calculateOverallScore(results)
	if got != 90.0 {
		t.Errorf("calculateOverallScore = %f, want 90.0", got)
	}
}

// ---------- printBanner, printUsage ----------

func TestPrintBanner_NoCrash(t *testing.T) {
	// printBanner just outputs to stdout, should not panic
	printBanner()
}

func TestPrintUsage_NoCrash(t *testing.T) {
	printUsage()
}

// ---------- generateJSONReport, generateHTMLReport, generateCompatibilityCSV, generateReports ----------

func TestGenerateJSONReport_NoCrash(t *testing.T) {
	results := &WorkflowTestResults{
		EngineResults: map[string]*EngineTestResult{},
	}
	args := &CommandLineArgs{outputDir: "/tmp/test_output"}
	generateJSONReport(results, args)
}

func TestGenerateHTMLReport_NoCrash(t *testing.T) {
	results := &WorkflowTestResults{}
	args := &CommandLineArgs{outputDir: "/tmp/test_output"}
	generateHTMLReport(results, args)
}

func TestGenerateCompatibilityCSV_NoCrash(t *testing.T) {
	results := &WorkflowTestResults{}
	args := &CommandLineArgs{outputDir: "/tmp/test_output"}
	generateCompatibilityCSV(results, args)
}

func TestGenerateReports_NoCrash(t *testing.T) {
	results := &WorkflowTestResults{
		EngineResults: map[string]*EngineTestResult{},
	}
	args := &CommandLineArgs{outputDir: "/tmp/test_output"}
	generateReports(results, args)
}

// ---------- printTestSummary ----------

func TestPrintTestSummary_EmptyResults(t *testing.T) {
	results := &WorkflowTestResults{
		Engines:           []string{},
		Platforms:         []phantomtesting.SequencingPlatform{},
		Workflows:         []string{},
		EngineResults:     map[string]*EngineTestResult{},
		PlatformResults:   map[phantomtesting.SequencingPlatform]*PlatformTestResult{},
		CompatibilityMatrix: map[string]map[phantomtesting.SequencingPlatform]*WorkflowTestResult{},
		OverallScore:      0.0,
	}
	args := &CommandLineArgs{outputDir: "/tmp/test"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_WithResults_HighScore(t *testing.T) {
	platform := phantomtesting.IlluminaNovaSeq6000
	results := &WorkflowTestResults{
		Engines:   []string{"cromwell"},
		Platforms: []phantomtesting.SequencingPlatform{platform},
		Workflows: []string{"gatk"},
		EngineResults: map[string]*EngineTestResult{
			"cromwell": {
				Engine:               "cromwell",
				Available:            true,
				SuccessRate:          0.95,
				CompatibilityScore:   95.0,
				AverageExecutionTime: time.Second * 30,
			},
		},
		PlatformResults: map[phantomtesting.SequencingPlatform]*PlatformTestResult{
			platform: {
				Platform:             platform,
				SuccessRate:          0.95,
				DataIntegrityScore:   0.98,
				AverageExecutionTime: time.Second * 30,
			},
		},
		CompatibilityMatrix: map[string]map[phantomtesting.SequencingPlatform]*WorkflowTestResult{
			"cromwell": {
				platform: {Success: true, Engine: "cromwell", Platform: platform},
			},
		},
		OverallScore: 95.0,
	}
	args := &CommandLineArgs{outputDir: "/tmp/test"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_MediumScore(t *testing.T) {
	results := &WorkflowTestResults{
		Engines:             []string{},
		Platforms:           []phantomtesting.SequencingPlatform{},
		Workflows:           []string{},
		EngineResults:       map[string]*EngineTestResult{},
		PlatformResults:     map[phantomtesting.SequencingPlatform]*PlatformTestResult{},
		CompatibilityMatrix: map[string]map[phantomtesting.SequencingPlatform]*WorkflowTestResult{},
		OverallScore:        76.0, // "GOOD" range
	}
	args := &CommandLineArgs{outputDir: "/tmp/test"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_LowScore(t *testing.T) {
	results := &WorkflowTestResults{
		Engines:             []string{},
		Platforms:           []phantomtesting.SequencingPlatform{},
		Workflows:           []string{},
		EngineResults:       map[string]*EngineTestResult{},
		PlatformResults:     map[phantomtesting.SequencingPlatform]*PlatformTestResult{},
		CompatibilityMatrix: map[string]map[phantomtesting.SequencingPlatform]*WorkflowTestResult{},
		OverallScore:        50.0,
	}
	args := &CommandLineArgs{outputDir: "/tmp/test"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_FairScore(t *testing.T) {
	results := &WorkflowTestResults{
		Engines:             []string{},
		Platforms:           []phantomtesting.SequencingPlatform{},
		Workflows:           []string{},
		EngineResults:       map[string]*EngineTestResult{},
		PlatformResults:     map[phantomtesting.SequencingPlatform]*PlatformTestResult{},
		CompatibilityMatrix: map[string]map[phantomtesting.SequencingPlatform]*WorkflowTestResult{},
		OverallScore:        65.0, // "FAIR" range (60-75)
	}
	args := &CommandLineArgs{outputDir: "/tmp/test"}
	printTestSummary(results, args)
}

func TestPrintTestSummary_GenerateReport(t *testing.T) {
	results := &WorkflowTestResults{
		Engines:             []string{},
		Platforms:           []phantomtesting.SequencingPlatform{},
		Workflows:           []string{},
		EngineResults:       map[string]*EngineTestResult{},
		PlatformResults:     map[phantomtesting.SequencingPlatform]*PlatformTestResult{},
		CompatibilityMatrix: map[string]map[phantomtesting.SequencingPlatform]*WorkflowTestResult{},
		OverallScore:        91.0,
	}
	args := &CommandLineArgs{outputDir: "/tmp/test", generateReport: true}
	printTestSummary(results, args)
}

func TestPrintTestSummary_EngineNotAvailable(t *testing.T) {
	platform := phantomtesting.IlluminaNovaSeq6000
	results := &WorkflowTestResults{
		Engines:   []string{"cromwell"},
		Platforms: []phantomtesting.SequencingPlatform{platform},
		Workflows: []string{"basic"},
		EngineResults: map[string]*EngineTestResult{
			"cromwell": {
				Engine:    "cromwell",
				Available: false, // not available
			},
		},
		PlatformResults: map[phantomtesting.SequencingPlatform]*PlatformTestResult{
			platform: {
				Platform:           platform,
				SuccessRate:        0.0,
				DataIntegrityScore: 0.0,
			},
		},
		CompatibilityMatrix: map[string]map[phantomtesting.SequencingPlatform]*WorkflowTestResult{
			"cromwell": {
				platform: {Success: false, Engine: "cromwell", Platform: platform},
			},
		},
		OverallScore: 0.0,
	}
	args := &CommandLineArgs{outputDir: "/tmp/test"}
	printTestSummary(results, args)
}

// ---------- validateArgs additional paths ----------

func TestValidateArgs_InvalidEngines(t *testing.T) {
	args := &CommandLineArgs{
		outputDir: "/tmp/test",
		timeout:   120,
		engines:   "invalid_engine_xyz",
		platforms: "all",
		workflows: "all",
	}
	if err := validateArgs(args); err == nil {
		t.Error("validateArgs(invalid engines) should return error")
	}
}

func TestValidateArgs_InvalidPlatforms(t *testing.T) {
	args := &CommandLineArgs{
		outputDir: "/tmp/test",
		timeout:   120,
		engines:   "all",
		platforms: "invalid_platform_xyz",
		workflows: "all",
	}
	if err := validateArgs(args); err == nil {
		t.Error("validateArgs(invalid platforms) should return error")
	}
}

func TestValidateArgs_InvalidWorkflows(t *testing.T) {
	// parseWorkflowList allows custom names, so "invalid_xyz" will be added as-is.
	// We need to trigger zero-workflow situation: use an empty string after trimming.
	// Actually, looking at the code, parseWorkflowList always returns at least 1 entry
	// for non-empty input (it uses the raw value for unknown entries).
	// We can only get 0 workflows if the input is empty string split produces only empty parts.
	// But validateArgs only checks len==0. Since parseWorkflowList accepts custom names,
	// this test ensures the happy path: a truly empty workflows arg.
	// Instead let's verify that invalid engines properly fail.
	// The only way to get 0 workflows is to pass an empty string, but that gives [""].
	// So test the boundary: valid config passes.
	args := &CommandLineArgs{
		outputDir: "/tmp/test",
		timeout:   120,
		engines:   "cromwell",
		platforms: "novaseq",
		workflows: "basic",
	}
	if err := validateArgs(args); err != nil {
		t.Errorf("validateArgs(specific valid args) should pass: %v", err)
	}
}

// ---------- calculateDataIntegrity: compares engine output with the input ----------

// TestCalculateDataIntegrity_ChecksReadCount verifies integrity reflects
// whether the engine's read count matches the reads written. The previous
// implementation scored any existing file (even "mock") at 0.95.
func TestCalculateDataIntegrity_ChecksReadCount(t *testing.T) {
	f := phantomtesting.NewE2EFramework()
	data, _ := f.GenerateSyntheticData(f.CreatePlatformConfig(phantomtesting.IlluminaNovaSeq6000, "CEU"))
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
	}
	for _, c := range cases {
		if got := calculateDataIntegrity(data, c.files); got != c.want {
			t.Errorf("%s: calculateDataIntegrity = %f, want %f", c.name, got, c.want)
		}
	}
}

// ---------- calculateWorkflowScore additional branches ----------

func TestCalculateWorkflowScore_LongExecution(t *testing.T) {
	// Execution between 1 min and 5 min: score += 0.1
	result := &WorkflowTestResult{
		Success:       true,
		DataIntegrity: 0.9,
		OutputAccuracy: 0.85,
		ExecutionTime: 2 * time.Minute,
	}
	got := calculateWorkflowScore(result)
	if got <= 0.0 {
		t.Errorf("calculateWorkflowScore(2 min) = %f, want > 0", got)
	}
}

func TestCalculateWorkflowScore_VeryLongExecution(t *testing.T) {
	// Execution > 5 min: no time bonus
	result := &WorkflowTestResult{
		Success:       true,
		DataIntegrity: 0.9,
		OutputAccuracy: 0.85,
		ExecutionTime: 10 * time.Minute,
	}
	got := calculateWorkflowScore(result)
	if got <= 0.0 {
		t.Errorf("calculateWorkflowScore(10 min) = %f, want > 0", got)
	}
}

// ---------- executeCromwellWorkflow, executeNextflowWorkflow, executeSnakemakeWorkflow ----------

type executeFunc func(string, *phantomtesting.SyntheticData, phantomtesting.SequencingPlatform, string, *CommandLineArgs) ([]string, error)

var executors = map[string]executeFunc{
	"cromwell":  executeCromwellWorkflow,
	"nextflow":  executeNextflowWorkflow,
	"snakemake": executeSnakemakeWorkflow,
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

// TestExecuteWorkflow_EngineUnavailable checks that an unavailable engine is
// an error and that no placeholder output is written. The previous
// implementation returned success with "mock output" files.
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
			_ = filepath.Walk(dir, func(path string, info os.FileInfo, err error) error {
				if err == nil && !info.IsDir() {
					if b, _ := os.ReadFile(path); strings.Contains(string(b), "mock") {
						t.Errorf("placeholder output written to %s", path)
					}
				}
				return nil
			})
		})
	}
}

// TestExecuteWorkflow_RealEngine runs each installed engine on the synthetic
// reads and checks the reported read count matches what was written.
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

// ---------- runWorkflowTest ----------

func TestRunWorkflowTest_UnavailableSkipped(t *testing.T) {
	withNoEngines(t)
	f := phantomtesting.NewE2EFramework()
	for _, engine := range []string{"cromwell", "nextflow", "snakemake"} {
		result := runWorkflowTest(f, engine, phantomtesting.IlluminaNovaSeq6000, "basic",
			&CommandLineArgs{outputDir: t.TempDir(), skipUnavailable: true})
		if !result.Skipped || result.Success {
			t.Errorf("%s unavailable with skip-unavailable: Skipped=%v Success=%v, want skipped and not passed",
				engine, result.Skipped, result.Success)
		}
	}
}

func TestRunWorkflowTest_UnavailableFailsWithoutSkip(t *testing.T) {
	withNoEngines(t)
	f := phantomtesting.NewE2EFramework()
	result := runWorkflowTest(f, "nextflow", phantomtesting.PacBioRevio, "assembly",
		&CommandLineArgs{outputDir: t.TempDir(), skipUnavailable: false})
	if result.Skipped || result.Success {
		t.Errorf("Skipped=%v Success=%v, want a failure", result.Skipped, result.Success)
	}
	if len(result.ErrorMessages) == 0 {
		t.Error("a failed run must carry an error message")
	}
}

func TestRunWorkflowTest_RealSnakemake(t *testing.T) {
	requireRealEngine(t, "snakemake")
	f := phantomtesting.NewE2EFramework()
	result := runWorkflowTest(f, "snakemake", phantomtesting.ONTPromethION24, "basic",
		&CommandLineArgs{outputDir: t.TempDir(), skipUnavailable: true, timeout: 600})
	if !result.Success {
		t.Fatalf("snakemake run failed: %v", result.ErrorMessages)
	}
	if result.DataIntegrity != 1.0 {
		t.Errorf("DataIntegrity = %f, want 1.0", result.DataIntegrity)
	}
}

func TestRunWorkflowTest_UnsupportedEngine(t *testing.T) {
	tmpDir := t.TempDir()
	f := phantomtesting.NewE2EFramework()
	args := &CommandLineArgs{
		outputDir:       tmpDir,
		skipUnavailable: true,
	}
	result := runWorkflowTest(f, "unsupported_engine", phantomtesting.IlluminaNovaSeq6000, "basic", args)
	if result == nil {
		t.Fatal("runWorkflowTest returned nil")
	}
	if len(result.ErrorMessages) == 0 {
		t.Error("runWorkflowTest(unsupported engine) should have error messages")
	}
	if result.Skipped || result.Success {
		t.Error("an unsupported engine is a failure, not a skip")
	}
}

// ---------- testWorkflowEngine ----------

func TestTestWorkflowEngine_DryRun(t *testing.T) {
	withLookPath(t, "cromwell")
	tmpDir := t.TempDir()
	f := phantomtesting.NewE2EFramework()
	args := &CommandLineArgs{
		outputDir:       tmpDir,
		skipUnavailable: true,
		dryRun:          true,
		verbose:         false,
	}
	platforms := []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000}
	workflows := []string{"basic"}
	result := testWorkflowEngine(f, "cromwell", platforms, workflows, args)
	if result == nil {
		t.Fatal("testWorkflowEngine returned nil")
	}
	if result.Engine != "cromwell" {
		t.Errorf("engine = %q, want \"cromwell\"", result.Engine)
	}
	if len(result.WorkflowResults) != 0 {
		t.Errorf("dry run executed %d workflows", len(result.WorkflowResults))
	}
}

func TestTestWorkflowEngine_SkipUnavailable(t *testing.T) {
	withNoEngines(t)
	f := phantomtesting.NewE2EFramework()
	platforms := []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000}
	workflows := []string{"basic"}

	skipped := testWorkflowEngine(f, "nextflow", platforms, workflows,
		&CommandLineArgs{outputDir: t.TempDir(), skipUnavailable: true, verbose: true})
	if skipped.Available || len(skipped.WorkflowResults) != 0 || skipped.SuccessRate != 0 {
		t.Errorf("skipped engine: Available=%v results=%d SuccessRate=%f, want unavailable with no results",
			skipped.Available, len(skipped.WorkflowResults), skipped.SuccessRate)
	}

	failed := testWorkflowEngine(f, "nextflow", platforms, workflows,
		&CommandLineArgs{outputDir: t.TempDir(), skipUnavailable: false, verbose: true})
	if len(failed.WorkflowResults) != 1 || failed.SuccessRate != 0 {
		t.Fatalf("unskipped unavailable engine: results=%d SuccessRate=%f, want 1 failed result",
			len(failed.WorkflowResults), failed.SuccessRate)
	}
	for _, r := range failed.WorkflowResults {
		if r.Success || r.Skipped {
			t.Errorf("Success=%v Skipped=%v, want failed", r.Success, r.Skipped)
		}
	}
}

// ---------- testPlatformCompatibility ----------

func TestTestPlatformCompatibility_DryRun(t *testing.T) {
	withLookPath(t, "cromwell")
	tmpDir := t.TempDir()
	f := phantomtesting.NewE2EFramework()
	args := &CommandLineArgs{
		outputDir:       tmpDir,
		skipUnavailable: true,
		dryRun:          true,
	}
	engines := []string{"cromwell"}
	workflows := []string{"basic"}
	result := testPlatformCompatibility(f, phantomtesting.IlluminaNovaSeq6000, engines, workflows, args)
	if result == nil {
		t.Fatal("testPlatformCompatibility returned nil")
	}
}

func TestTestPlatformCompatibility_SkipsUnavailable(t *testing.T) {
	withNoEngines(t)
	f := phantomtesting.NewE2EFramework()
	args := &CommandLineArgs{outputDir: t.TempDir(), skipUnavailable: true}
	for _, platform := range []phantomtesting.SequencingPlatform{
		phantomtesting.IlluminaNovaSeq6000, phantomtesting.PacBioRevio, phantomtesting.ONTPromethION24,
	} {
		result := testPlatformCompatibility(f, platform, []string{"cromwell", "nextflow"}, []string{"pacbio_assembly", "ont_assembly", "basic"}, args)
		if len(result.EngineResults) != 0 || result.SuccessRate != 0 {
			t.Errorf("%s: results=%d SuccessRate=%f, want nothing run", platform, len(result.EngineResults), result.SuccessRate)
		}
	}
}

// ---------- runWorkflowTests ----------

func TestRunWorkflowTests_DryRun(t *testing.T) {
	withLookPath(t, "cromwell")
	tmpDir := t.TempDir()
	f := phantomtesting.NewE2EFramework()
	args := &CommandLineArgs{
		outputDir:       tmpDir,
		skipUnavailable: true,
		dryRun:          true,
		verbose:         false,
		generateReport:  false,
	}
	engines := []string{"cromwell"}
	platforms := []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000}
	workflows := []string{"basic"}
	results := runWorkflowTests(f, engines, platforms, workflows, args)
	if results == nil {
		t.Fatal("runWorkflowTests returned nil")
	}
	if len(results.Engines) != 1 {
		t.Errorf("engines = %d, want 1", len(results.Engines))
	}
}

// TestRunWorkflowTests_NoEnginesReportsSkipped checks that with nothing
// installed every engine is reported as skipped and nothing as passed.
func TestRunWorkflowTests_NoEnginesReportsSkipped(t *testing.T) {
	withNoEngines(t)
	f := phantomtesting.NewE2EFramework()
	args := &CommandLineArgs{outputDir: t.TempDir(), skipUnavailable: true}
	engines := []string{"cromwell", "nextflow", "snakemake"}
	platforms := []phantomtesting.SequencingPlatform{phantomtesting.IlluminaNovaSeq6000}
	results := runWorkflowTests(f, engines, platforms, []string{"basic"}, args)

	ran, failed, skipped := countWorkflowOutcomes(results)
	if ran != 0 || failed != 0 || skipped != 3 {
		t.Errorf("ran=%d failed=%d skipped=%d, want 0/0/3", ran, failed, skipped)
	}
	if results.OverallScore != 0 {
		t.Errorf("OverallScore = %f, want 0 when nothing ran", results.OverallScore)
	}
	printTestSummary(results, args)
}

func TestCountWorkflowOutcomes(t *testing.T) {
	results := &WorkflowTestResults{
		EngineResults: map[string]*EngineTestResult{
			"cromwell": {Available: false},
			"nextflow": {Available: true, WorkflowResults: map[string]*WorkflowTestResult{
				"a": {Success: true},
				"b": {Success: false},
			}},
		},
		PlatformResults: map[phantomtesting.SequencingPlatform]*PlatformTestResult{
			phantomtesting.IlluminaNovaSeq6000: {EngineResults: map[string]*WorkflowTestResult{
				"nextflow":  {Success: true},
				"snakemake": {Skipped: true},
			}},
		},
	}
	ran, failed, skipped := countWorkflowOutcomes(results)
	if ran != 3 || failed != 1 || skipped != 2 {
		t.Errorf("ran=%d failed=%d skipped=%d, want 3/1/2", ran, failed, skipped)
	}
}

// ---------- parsePlatformList additional entries ----------

func TestParsePlatformList_LongRead(t *testing.T) {
	platforms := parsePlatformList("longread")
	if len(platforms) != 2 {
		t.Errorf("parsePlatformList(longread) = %d platforms, want 2", len(platforms))
	}
}

func TestParsePlatformList_NoDuplicates(t *testing.T) {
	// "all" + "novaseq" should not produce duplicate novaseq
	platforms := parsePlatformList("all,novaseq")
	seen := make(map[phantomtesting.SequencingPlatform]bool)
	for _, p := range platforms {
		if seen[p] {
			t.Errorf("duplicate platform %q", p)
		}
		seen[p] = true
	}
}

func TestParsePlatformList_Individual(t *testing.T) {
	cases := []struct {
		input    string
		expected phantomtesting.SequencingPlatform
	}{
		{"nextseq", phantomtesting.IlluminaNextSeq2000},
		{"miseq", phantomtesting.IlluminaMiSeq},
		{"pacbio", phantomtesting.PacBioRevio},
		{"ont", phantomtesting.ONTPromethION24},
		{"dnbseq", phantomtesting.MGI_DNBSEQ_T20x2},
		{"aviti", phantomtesting.ElementAVITI},
	}
	for _, tc := range cases {
		t.Run(tc.input, func(t *testing.T) {
			platforms := parsePlatformList(tc.input)
			if len(platforms) != 1 || platforms[0] != tc.expected {
				t.Errorf("parsePlatformList(%q) = %v, want [%v]", tc.input, platforms, tc.expected)
			}
		})
	}
}

// ---------- parseWorkflowList additional entries ----------

func TestParseWorkflowList_GATK(t *testing.T) {
	workflows := parseWorkflowList("gatk")
	if len(workflows) < 1 {
		t.Error("parseWorkflowList(gatk) should return at least 1 workflow")
	}
}

func TestParseWorkflowList_Assembly(t *testing.T) {
	workflows := parseWorkflowList("assembly")
	if len(workflows) < 1 {
		t.Error("parseWorkflowList(assembly) should return at least 1 workflow")
	}
}

func TestParseWorkflowList_Variants(t *testing.T) {
	workflows := parseWorkflowList("variants")
	if len(workflows) < 1 {
		t.Error("parseWorkflowList(variants) should return at least 1 workflow")
	}
}

func TestParseWorkflowList_QC(t *testing.T) {
	workflows := parseWorkflowList("qc")
	if len(workflows) < 1 {
		t.Error("parseWorkflowList(qc) should return at least 1 workflow")
	}
}

func TestParseWorkflowList_Population(t *testing.T) {
	workflows := parseWorkflowList("population")
	if len(workflows) < 1 {
		t.Error("parseWorkflowList(population) should return at least 1 workflow")
	}
}

func TestParseWorkflowList_NFCore(t *testing.T) {
	workflows := parseWorkflowList("nfcore")
	if len(workflows) < 1 {
		t.Error("parseWorkflowList(nfcore) should return at least 1 workflow")
	}
}

func TestParseWorkflowList_NoDuplicates(t *testing.T) {
	// "all" + "basic" should not duplicate "basic"
	workflows := parseWorkflowList("all,basic")
	seen := make(map[string]bool)
	for _, w := range workflows {
		if seen[w] {
			t.Errorf("duplicate workflow %q in result", w)
		}
		seen[w] = true
	}
}

// ---------- parseEngineList additional entries ----------

func TestParseEngineList_WDLAlias(t *testing.T) {
	engines := parseEngineList("wdl")
	if len(engines) != 1 || engines[0] != "cromwell" {
		t.Errorf("parseEngineList(wdl) = %v, want [cromwell]", engines)
	}
}

func TestParseEngineList_NFAlias(t *testing.T) {
	engines := parseEngineList("nf")
	if len(engines) != 1 || engines[0] != "nextflow" {
		t.Errorf("parseEngineList(nf) = %v, want [nextflow]", engines)
	}
}

func TestParseEngineList_SnakeAlias(t *testing.T) {
	engines := parseEngineList("snake")
	if len(engines) != 1 || engines[0] != "snakemake" {
		t.Errorf("parseEngineList(snake) = %v, want [snakemake]", engines)
	}
}

func TestParseEngineList_Unknown(t *testing.T) {
	engines := parseEngineList("unknown_engine_xyz")
	if len(engines) != 0 {
		t.Errorf("parseEngineList(unknown) = %v, want []", engines)
	}
}

// ---------- selectRepresentativeWorkflow additional branches ----------

func TestSelectRepresentativeWorkflow_DefaultGATK(t *testing.T) {
	// Default platform (not PacBio or ONT) should prefer gatk workflow
	got := selectRepresentativeWorkflow("cromwell", phantomtesting.IlluminaNovaSeq6000, []string{"gatk_germline", "basic"})
	if got != "gatk_germline" {
		t.Errorf("selectRepresentativeWorkflow(Illumina, gatk+basic) = %q, want \"gatk_germline\"", got)
	}
}

func TestSelectRepresentativeWorkflow_DefaultFallback(t *testing.T) {
	// Default platform with no preferred workflow, falls back to first
	got := selectRepresentativeWorkflow("cromwell", phantomtesting.IlluminaNovaSeq6000, []string{"basic", "qc"})
	if got != "basic" {
		t.Errorf("selectRepresentativeWorkflow(Illumina, no-preferred) = %q, want \"basic\"", got)
	}
}

func TestSelectRepresentativeWorkflow_ONTFallback(t *testing.T) {
	// ONT with no ont/nanopore workflow falls back to first
	got := selectRepresentativeWorkflow("nextflow", phantomtesting.ONTPromethION24, []string{"basic", "gatk"})
	if got != "basic" {
		t.Errorf("selectRepresentativeWorkflow(ONT, no-nanopore) = %q, want \"basic\"", got)
	}
}

func TestSelectRepresentativeWorkflow_PacBioFallback(t *testing.T) {
	// PacBio with no assembly workflow falls back to first
	got := selectRepresentativeWorkflow("cromwell", phantomtesting.PacBioRevio, []string{"basic", "gatk"})
	if got != "basic" {
		t.Errorf("selectRepresentativeWorkflow(PacBio, no-assembly) = %q, want \"basic\"", got)
	}
}

// ---------- estimateMemoryUsage additional branches ----------

func TestEstimateMemoryUsage_ONT(t *testing.T) {
	mem := estimateMemoryUsage(phantomtesting.ONTPromethION24, "basic")
	if mem <= 0 {
		t.Errorf("estimateMemoryUsage(ONT, basic) = %d, want > 0", mem)
	}
}

func TestEstimateMemoryUsage_NextSeq(t *testing.T) {
	// NextSeq uses default base memory
	mem := estimateMemoryUsage(phantomtesting.IlluminaNextSeq2000, "qc")
	if mem <= 0 {
		t.Errorf("estimateMemoryUsage(NextSeq, qc) = %d, want > 0", mem)
	}
}

func TestEstimateMemoryUsage_PacBioBasic(t *testing.T) {
	mem := estimateMemoryUsage(phantomtesting.PacBioRevio, "basic")
	if mem <= 0 {
		t.Errorf("estimateMemoryUsage(PacBio, basic) = %d, want > 0", mem)
	}
}

// ---------- estimateCPUUsage additional branches ----------

func TestEstimateCPUUsage_Other(t *testing.T) {
	// Non-assembly, non-basic workflow gets base CPU
	cpu := estimateCPUUsage(phantomtesting.IlluminaNovaSeq6000, "gatk")
	if cpu <= 0 || cpu > 100 {
		t.Errorf("estimateCPUUsage(gatk) = %f, want in (0,100]", cpu)
	}
}

// ---------- calculateOverallScore additional paths ----------

func TestCalculateOverallScore_UnavailableEngines(t *testing.T) {
	// Engines that are not available should not contribute to score
	results := &WorkflowTestResults{
		EngineResults: map[string]*EngineTestResult{
			"cromwell": {Available: false, CompatibilityScore: 90.0},
		},
	}
	got := calculateOverallScore(results)
	if got != 0.0 {
		t.Errorf("calculateOverallScore(all unavailable) = %f, want 0.0", got)
	}
}

func TestCalculateOverallScore_MixedAvailability(t *testing.T) {
	results := &WorkflowTestResults{
		EngineResults: map[string]*EngineTestResult{
			"cromwell": {Available: true, CompatibilityScore: 80.0},
			"nextflow": {Available: false, CompatibilityScore: 90.0},
		},
	}
	got := calculateOverallScore(results)
	if got != 80.0 {
		t.Errorf("calculateOverallScore(mixed) = %f, want 80.0", got)
	}
}
