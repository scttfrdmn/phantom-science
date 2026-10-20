//go:build !release

package integration

import (
	"context"
	"encoding/json"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	pt "github.com/scttfrdmn/phantom-science/internal/testing"
)

// These tests run real workflow engines on synthetic reads. Each test skips
// when its engine, or a tool its workflow calls, is not installed; once a
// workflow runs, any failure fails the test. nf-core pipelines download
// pipelines and containers, so they only run when PHANTOM_NFCORE_TESTS=1.

// workflowTimeout bounds a single engine run.
const workflowTimeout = 30 * time.Minute

// TestCromwellWDLIntegration tests integration with Cromwell workflow engine
func TestCromwellWDLIntegration(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping Cromwell WDL integration tests in short mode")
	}

	if !isWorkflowEngineAvailable("cromwell") {
		t.Skip("Cromwell not available - skipping WDL integration tests")
	}

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	// Test different genomics workflows
	workflows := []WDLWorkflow{
		{
			Name:          "GATK_Germline_Variants",
			Description:   "GATK best practices germline variant calling",
			WDLFile:       "gatk_germline.wdl",
			Platforms:     []pt.SequencingPlatform{pt.IlluminaNovaSeq6000, pt.IlluminaNextSeq2000},
			RequiredTools: []string{"bwa", "samtools", "picard", "gatk"},
		},
		{
			Name:          "PacBio_HiFi_Assembly",
			Description:   "PacBio HiFi genome assembly workflow",
			WDLFile:       "pacbio_hifi_assembly.wdl",
			Platforms:     []pt.SequencingPlatform{pt.PacBioRevio},
			RequiredTools: []string{"hifiasm", "quast"},
		},
		{
			Name:          "ONT_Long_Read_Variants",
			Description:   "Oxford Nanopore long-read alignment (first stage of variant calling)",
			WDLFile:       "ont_longread_variants.wdl",
			Platforms:     []pt.SequencingPlatform{pt.ONTPromethION24},
			RequiredTools: []string{"minimap2", "samtools"},
		},
	}

	for _, workflow := range workflows {
		t.Run(workflow.Name, func(t *testing.T) {
			testCromwellWorkflow(t, framework, workflow)
		})
	}
}

// TestNextflowIntegration tests integration with Nextflow workflow engine
func TestNextflowIntegration(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping Nextflow integration tests in short mode")
	}

	if !isWorkflowEngineAvailable("nextflow") {
		t.Skip("Nextflow not available - skipping integration tests")
	}

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	workflows := []NextflowWorkflow{
		{
			Name:        "local/read_count",
			Description: "Local pipeline that counts FASTQ records",
			Platforms:   []pt.SequencingPlatform{pt.IlluminaNovaSeq6000, pt.PacBioRevio, pt.ONTPromethION24},
		},
		{
			Name:        "nf-core/sarek",
			Description: "Comprehensive variant calling and annotation pipeline",
			Pipeline:    "nf-core/sarek",
			Version:     "3.4.0",
			Platforms:   []pt.SequencingPlatform{pt.IlluminaNovaSeq6000, pt.IlluminaNextSeq2000},
		},
		{
			Name:        "nf-core/nanoseq",
			Description: "Nanopore sequencing analysis pipeline",
			Pipeline:    "nf-core/nanoseq",
			Version:     "3.1.0",
			Platforms:   []pt.SequencingPlatform{pt.ONTPromethION24},
		},
		{
			Name:        "nf-core/hlatyping",
			Description: "HLA typing from sequencing data",
			Pipeline:    "nf-core/hlatyping",
			Version:     "2.0.0",
			Platforms:   []pt.SequencingPlatform{pt.IlluminaNovaSeq6000},
		},
	}

	for _, workflow := range workflows {
		t.Run(workflow.Name, func(t *testing.T) {
			testNextflowWorkflow(t, framework, workflow)
		})
	}
}

// TestSnakemakeIntegration tests integration with Snakemake workflow engine
func TestSnakemakeIntegration(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping Snakemake integration tests in short mode")
	}

	if !isWorkflowEngineAvailable("snakemake") {
		t.Skip("Snakemake not available - skipping integration tests")
	}

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	snakemakeWorkflows := []SnakemakeWorkflow{
		{
			Name:          "Multi_Platform_QC",
			Description:   "Quality control across multiple sequencing platforms",
			Snakefile:     "multi_platform_qc.smk",
			Platforms:     []pt.SequencingPlatform{pt.IlluminaNovaSeq6000, pt.PacBioRevio, pt.ONTPromethION24},
			RequiredTools: []string{"fastqc", "multiqc"},
		},
		{
			Name:          "Population_Genomics",
			Description:   "VCF quality filtering and summary statistics",
			Snakefile:     "population_genomics.smk",
			Platforms:     []pt.SequencingPlatform{pt.IlluminaNovaSeq6000},
			RequiredTools: []string{"bcftools"},
		},
		{
			Name:          "Structural_Variants",
			Description:   "Long-read alignment (first stage of structural variant detection)",
			Snakefile:     "structural_variants.smk",
			Platforms:     []pt.SequencingPlatform{pt.PacBioRevio, pt.ONTPromethION24},
			RequiredTools: []string{"minimap2", "samtools"},
		},
	}

	for _, workflow := range snakemakeWorkflows {
		t.Run(workflow.Name, func(t *testing.T) {
			testSnakemakeWorkflow(t, framework, workflow)
		})
	}
}

// TestWorkflowCompatibilityMatrix runs a read-count workflow on every
// installed engine for each platform and checks the count matches.
func TestWorkflowCompatibilityMatrix(t *testing.T) {
	if testing.Short() {
		t.Skip("Skipping workflow compatibility matrix tests in short mode")
	}

	framework := pt.NewE2EFramework()
	require.NotNil(t, framework)

	platforms := []pt.SequencingPlatform{
		pt.IlluminaNovaSeq6000,
		pt.PacBioRevio,
		pt.ONTPromethION24,
		pt.MGI_DNBSEQ_T20x2,
		pt.ElementAVITI,
	}

	engines := []string{"cromwell", "nextflow", "snakemake"}

	compatibilityMatrix := make(map[string]map[pt.SequencingPlatform]WorkflowCompatibilityResult)

	for _, engine := range engines {
		if !isWorkflowEngineAvailable(engine) {
			t.Logf("Workflow engine %s not available, skipping", engine)
			continue
		}

		compatibilityMatrix[engine] = make(map[pt.SequencingPlatform]WorkflowCompatibilityResult)

		for _, platform := range platforms {
			t.Run(fmt.Sprintf("%s_%s", engine, platform), func(t *testing.T) {
				result := testWorkflowCompatibility(t, framework, engine, platform)
				compatibilityMatrix[engine][platform] = result
				require.True(t, result.Success, "%s on %s failed: %s",
					engine, platform, strings.Join(result.ErrorMessages, "; "))
				assert.Equal(t, 1.0, result.QualityMetrics.DataIntegrity,
					"%s should report the number of reads written", engine)
			})
		}
	}

	if len(compatibilityMatrix) == 0 {
		t.Skip("no workflow engine installed")
	}
	generateCompatibilityReport(t, compatibilityMatrix)
}

// Data structures for workflow testing

type WDLWorkflow struct {
	Name          string
	Description   string
	WDLFile       string
	Platforms     []pt.SequencingPlatform
	RequiredTools []string
}

type NextflowWorkflow struct {
	Name        string
	Description string
	// Pipeline is a remote nf-core pipeline; empty means the local read-count script.
	Pipeline  string
	Version   string
	Platforms []pt.SequencingPlatform
}

type SnakemakeWorkflow struct {
	Name          string
	Description   string
	Snakefile     string
	Platforms     []pt.SequencingPlatform
	RequiredTools []string
}

type WorkflowCompatibilityResult struct {
	Platform       pt.SequencingPlatform
	Engine         string
	Success        bool
	ExecutionTime  time.Duration
	OutputFiles    []string
	ErrorMessages  []string
	QualityMetrics WorkflowQualityMetrics
}

// WorkflowQualityMetrics holds values measured from a workflow run.
type WorkflowQualityMetrics struct {
	// DataIntegrity is the fraction of expected outputs that were produced
	// and are valid (non-empty, or matching the input read count).
	DataIntegrity float64
	// ProcessingSpeed is input megabytes divided by wall-clock seconds,
	// including engine start-up.
	ProcessingSpeed float64
}

type WorkflowResults struct {
	WorkflowName    string
	Platform        pt.SequencingPlatform
	Engine          string
	Success         bool
	ExecutionTime   time.Duration
	OutputDirectory string
	OutputFiles     []string
	// EngineOutput is the tail of the engine's combined stdout/stderr.
	EngineOutput string
	Metrics      WorkflowQualityMetrics
}

// workflowInputs are the synthetic input files written for a workflow.
type workflowInputs struct {
	R1, R2, Reference, VCF string
	Reads                  int
	Bytes                  int64
}

// Implementation functions

// cromwellCommand returns the program and leading arguments that invoke
// Cromwell: a "cromwell" executable, or "java -jar $CROMWELL_JAR".
func cromwellCommand() (string, []string, bool) {
	if path, err := exec.LookPath("cromwell"); err == nil {
		return path, nil, true
	}
	if jar := os.Getenv("CROMWELL_JAR"); jar != "" {
		if _, err := os.Stat(jar); err == nil {
			if java, err := exec.LookPath("java"); err == nil {
				return java, []string{"-jar", jar}, true
			}
		}
	}
	return "", nil, false
}

func isWorkflowEngineAvailable(engine string) bool {
	switch engine {
	case "cromwell":
		_, _, ok := cromwellCommand()
		return ok
	case "nextflow", "snakemake":
		return isToolAvailable(engine)
	default:
		return false
	}
}

// skipIfMissingTools skips the test when any tool is not on PATH.
func skipIfMissingTools(t *testing.T, tools []string) {
	t.Helper()
	for _, tool := range tools {
		if !isToolAvailable(tool) {
			t.Skipf("%s not installed - required by this workflow", tool)
		}
	}
}

// writeWorkflowInputs writes FASTQ, reference and VCF files for data into dir.
func writeWorkflowInputs(t *testing.T, framework *pt.E2EFramework, data *pt.SyntheticData, dir string) workflowInputs {
	t.Helper()
	require.NoError(t, os.MkdirAll(dir, 0750))
	in := workflowInputs{
		R1:        filepath.Join(dir, "synthetic_R1.fastq"),
		R2:        filepath.Join(dir, "synthetic_R2.fastq"),
		Reference: filepath.Join(dir, "reference.fasta"),
		VCF:       filepath.Join(dir, "variants.vcf"),
		Reads:     len(data.Sequences),
	}
	require.NoError(t, framework.WriteFASTQ(data, in.R1, in.R2))
	require.NoError(t, framework.WriteReferenceFromData(data, in.Reference))
	require.NoError(t, framework.WriteVCF(data, in.VCF))
	if info, err := os.Stat(in.R1); err == nil {
		in.Bytes = info.Size()
	}
	return in
}

// runEngine runs a command in dir and returns the tail of its output.
func runEngine(dir string, env []string, program string, args ...string) (string, error) {
	ctx, cancel := context.WithTimeout(context.Background(), workflowTimeout)
	defer cancel()
	cmd := exec.CommandContext(ctx, program, args...)
	cmd.Dir = dir
	cmd.Env = append(os.Environ(), env...)
	out, err := cmd.CombinedOutput()
	return tailLines(string(out), 40), err
}

// outputIntegrity returns the fraction of files that exist and are non-empty.
func outputIntegrity(files []string) float64 {
	if len(files) == 0 {
		return 0
	}
	ok := 0
	for _, f := range files {
		if info, err := os.Stat(f); err == nil && info.Size() > 0 {
			ok++
		}
	}
	return float64(ok) / float64(len(files))
}

// readCountMatches reports whether file holds the number want.
func readCountMatches(file string, want int) bool {
	content, err := os.ReadFile(file) // #nosec G304 -- workflow output in a test temp dir
	if err != nil {
		return false
	}
	got, err := strconv.ParseFloat(strings.TrimSpace(string(content)), 64)
	return err == nil && got == float64(want)
}

func processingSpeed(bytes int64, elapsed time.Duration) float64 {
	if elapsed <= 0 {
		return 0
	}
	return float64(bytes) / (1024 * 1024) / elapsed.Seconds()
}

func testCromwellWorkflow(t *testing.T, framework *pt.E2EFramework, workflow WDLWorkflow) {
	t.Logf("🔄 Testing Cromwell WDL workflow: %s", workflow.Name)

	for _, platform := range workflow.Platforms {
		t.Run(string(platform), func(t *testing.T) {
			skipIfMissingTools(t, workflow.RequiredTools)

			syntheticData := generateOrSkip(t, framework, platform, "CEU")

			workDir := t.TempDir()
			inputs := writeWorkflowInputs(t, framework, syntheticData, filepath.Join(workDir, "data"))

			wdlFile := filepath.Join(workDir, workflow.WDLFile)
			require.NoError(t, createWDLWorkflow(wdlFile, workflow.Name, platform))

			inputsFile := filepath.Join(workDir, "inputs.json")
			require.NoError(t, createWDLInputs(inputsFile, workflow.Name, inputs, platform))

			result := runCromwellWorkflow(t, workDir, wdlFile, inputsFile, workflow.Name, platform, inputs.Bytes)

			validateWorkflowResults(t, result, workflow.Name, platform)

			t.Logf("✅ Cromwell workflow %s completed for %s in %v",
				workflow.Name, platform, result.ExecutionTime)
		})
	}
}

func testNextflowWorkflow(t *testing.T, framework *pt.E2EFramework, workflow NextflowWorkflow) {
	t.Logf("🌊 Testing Nextflow workflow: %s", workflow.Name)

	for _, platform := range workflow.Platforms {
		t.Run(string(platform), func(t *testing.T) {
			if workflow.Pipeline != "" {
				if os.Getenv("PHANTOM_NFCORE_TESTS") != "1" {
					t.Skip("nf-core pipelines need network access and a container runtime; set PHANTOM_NFCORE_TESTS=1 to run")
				}
				skipIfMissingTools(t, []string{"docker"})
			}

			syntheticData := generateOrSkip(t, framework, platform, "CEU")

			workDir := t.TempDir()
			inputs := writeWorkflowInputs(t, framework, syntheticData, filepath.Join(workDir, "data"))

			require.NoError(t, prepareNextflowInputs(workDir, inputs, workflow, platform))

			result := runNextflowWorkflow(t, workDir, workflow, platform, inputs)

			validateWorkflowResults(t, result, workflow.Name, platform)

			t.Logf("✅ Nextflow workflow %s completed for %s in %v",
				workflow.Name, platform, result.ExecutionTime)
		})
	}
}

func testSnakemakeWorkflow(t *testing.T, framework *pt.E2EFramework, workflow SnakemakeWorkflow) {
	t.Logf("🐍 Testing Snakemake workflow: %s", workflow.Name)

	for _, platform := range workflow.Platforms {
		t.Run(string(platform), func(t *testing.T) {
			skipIfMissingTools(t, workflow.RequiredTools)

			syntheticData := generateOrSkip(t, framework, platform, "CEU")

			workDir := t.TempDir()
			inputs := writeWorkflowInputs(t, framework, syntheticData, filepath.Join(workDir, "data"))

			snakeFile := filepath.Join(workDir, workflow.Snakefile)
			expected, err := createSnakemakeWorkflow(snakeFile, workflow.Name, platform)
			require.NoError(t, err)

			result := runSnakemakeWorkflow(t, workDir, snakeFile, workflow.Name, platform, expected, inputs.Bytes)

			validateWorkflowResults(t, result, workflow.Name, platform)

			t.Logf("✅ Snakemake workflow %s completed for %s in %v",
				workflow.Name, platform, result.ExecutionTime)
		})
	}
}

func testWorkflowCompatibility(t *testing.T, framework *pt.E2EFramework, engine string, platform pt.SequencingPlatform) WorkflowCompatibilityResult {
	t.Logf("🧪 Testing %s compatibility with %s", engine, platform)

	start := time.Now()
	result := WorkflowCompatibilityResult{Platform: platform, Engine: engine}

	syntheticData := generateOrSkip(t, framework, platform, "CEU")

	workDir := t.TempDir()
	inputs := writeWorkflowInputs(t, framework, syntheticData, filepath.Join(workDir, "data"))

	outputs, err := runReadCountWorkflow(workDir, engine, inputs)
	result.ExecutionTime = time.Since(start)
	if err != nil {
		result.ErrorMessages = append(result.ErrorMessages, err.Error())
		return result
	}

	result.OutputFiles = outputs
	matched := 0
	for _, f := range outputs {
		if readCountMatches(f, inputs.Reads) {
			matched++
		}
	}
	if len(outputs) > 0 {
		result.QualityMetrics.DataIntegrity = float64(matched) / float64(len(outputs))
	}
	result.QualityMetrics.ProcessingSpeed = processingSpeed(inputs.Bytes, result.ExecutionTime)
	result.Success = result.QualityMetrics.DataIntegrity == 1.0
	if !result.Success {
		result.ErrorMessages = append(result.ErrorMessages,
			fmt.Sprintf("read count output %v does not match %d reads", outputs, inputs.Reads))
	}
	return result
}

// runReadCountWorkflow runs a minimal read-count workflow on engine and
// returns the read-count output files.
func runReadCountWorkflow(workDir, engine string, inputs workflowInputs) ([]string, error) {
	switch engine {
	case "cromwell":
		wdl := filepath.Join(workDir, "read_count.wdl")
		if err := os.WriteFile(wdl, []byte(readCountWDL), 0600); err != nil {
			return nil, err
		}
		inputsFile := filepath.Join(workDir, "read_count_inputs.json")
		if err := writeJSON(inputsFile, map[string]interface{}{"ReadCount.fastq": inputs.R1}); err != nil {
			return nil, err
		}
		outputs, out, err := cromwellRun(workDir, wdl, inputsFile)
		if err != nil {
			return nil, fmt.Errorf("cromwell run: %v\n%s", err, out)
		}
		return outputs, nil
	case "nextflow":
		script := filepath.Join(workDir, "read_count.nf")
		if err := os.WriteFile(script, []byte(readCountNextflow), 0600); err != nil {
			return nil, err
		}
		outDir := filepath.Join(workDir, "results")
		out, err := runEngine(workDir, nextflowEnv, "nextflow", "run", script, "--fastq", inputs.R1, "--outdir", outDir)
		if err != nil {
			return nil, fmt.Errorf("nextflow run: %v\n%s", err, out)
		}
		return []string{filepath.Join(outDir, "read_count.txt")}, nil
	case "snakemake":
		snakefile := filepath.Join(workDir, "Snakefile")
		content := fmt.Sprintf(readCountSnakemake, inputs.R1)
		if err := os.WriteFile(snakefile, []byte(content), 0600); err != nil {
			return nil, err
		}
		out, err := runEngine(workDir, nil, "snakemake", "-s", snakefile, "--cores", "1", "--directory", workDir)
		if err != nil {
			return nil, fmt.Errorf("snakemake: %v\n%s", err, out)
		}
		return []string{filepath.Join(workDir, "read_count.txt")}, nil
	default:
		return nil, fmt.Errorf("unsupported engine %s", engine)
	}
}

var nextflowEnv = []string{"NXF_ANSI_LOG=false"}

const readCountWDL = `version 1.0

workflow ReadCount {
    input {
        File fastq
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
        awk 'END { print NR / 4 }' ~{fastq} > read_count.txt
    >>>

    output {
        File read_count = "read_count.txt"
    }
}
`

const readCountNextflow = `#!/usr/bin/env nextflow

params.fastq = null
params.outdir = 'results'

process countReads {
    publishDir params.outdir, mode: 'copy'

    input:
    path fastq

    output:
    path "read_count.txt"

    script:
    """
    awk 'END { print NR / 4 }' ${fastq} > read_count.txt
    """
}

workflow {
    countReads(Channel.fromPath(params.fastq, checkIfExists: true))
}
`

// readCountSnakemake takes the FASTQ path via fmt.Sprintf.
const readCountSnakemake = `rule all:
    input:
        "read_count.txt"

rule count_reads:
    input:
        %q
    output:
        "read_count.txt"
    shell:
        "awk 'END {{ print NR / 4 }}' {input} > {output}"
`

func writeJSON(path string, v interface{}) error {
	b, err := json.MarshalIndent(v, "", "    ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, b, 0600)
}

// cromwellRun runs a WDL with Cromwell and returns the output file paths
// listed in its metadata.
func cromwellRun(workDir, wdlFile, inputsFile string) ([]string, string, error) {
	program, lead, ok := cromwellCommand()
	if !ok {
		return nil, "", fmt.Errorf("cromwell not available")
	}
	metadata := filepath.Join(workDir, "metadata.json")
	args := append(lead, "run", wdlFile, "-i", inputsFile, "-m", metadata)
	out, err := runEngine(workDir, nil, program, args...)
	if err != nil {
		return nil, out, err
	}

	raw, err := os.ReadFile(metadata) // #nosec G304 -- cromwell metadata in a test temp dir
	if err != nil {
		return nil, out, fmt.Errorf("read cromwell metadata: %w", err)
	}
	var meta struct {
		Status  string                 `json:"status"`
		Outputs map[string]interface{} `json:"outputs"`
	}
	if err := json.Unmarshal(raw, &meta); err != nil {
		return nil, out, fmt.Errorf("parse cromwell metadata: %w", err)
	}
	if meta.Status != "Succeeded" {
		return nil, out, fmt.Errorf("cromwell workflow status %q", meta.Status)
	}
	var files []string
	var collect func(v interface{})
	collect = func(v interface{}) {
		switch x := v.(type) {
		case string:
			files = append(files, x)
		case []interface{}:
			for _, e := range x {
				collect(e)
			}
		}
	}
	for _, v := range meta.Outputs {
		collect(v)
	}
	return files, out, nil
}

// Workflow creation functions

func createWDLWorkflow(wdlFile, workflowName string, platform pt.SequencingPlatform) error {
	var wdlContent string

	switch workflowName {
	case "GATK_Germline_Variants":
		wdlContent = `version 1.0

workflow GATKGermlineVariants {
    input {
        File fastq_r1
        File fastq_r2
        File reference_genome
        String sample_name
    }

    call AlignReads {
        input:
            fastq_r1 = fastq_r1,
            fastq_r2 = fastq_r2,
            reference = reference_genome,
            sample_name = sample_name
    }

    call MarkDuplicates {
        input:
            input_bam = AlignReads.aligned_bam,
            sample_name = sample_name
    }

    call HaplotypeCaller {
        input:
            input_bam = MarkDuplicates.dedup_bam,
            reference = reference_genome,
            sample_name = sample_name
    }

    output {
        File final_vcf = HaplotypeCaller.output_vcf
        File final_bam = MarkDuplicates.dedup_bam
    }
}

task AlignReads {
    input {
        File fastq_r1
        File fastq_r2
        File reference
        String sample_name
    }

    command <<<
        set -euo pipefail
        cp ~{reference} ref.fasta
        bwa index ref.fasta
        bwa mem -R '@RG\tID:~{sample_name}\tSM:~{sample_name}\tPL:ILLUMINA' \
            ref.fasta ~{fastq_r1} ~{fastq_r2} | \
            samtools sort -o ~{sample_name}.bam
        samtools index ~{sample_name}.bam
    >>>

    output {
        File aligned_bam = "~{sample_name}.bam"
        File aligned_bai = "~{sample_name}.bam.bai"
    }
}

task MarkDuplicates {
    input {
        File input_bam
        String sample_name
    }

    command <<<
        set -euo pipefail
        picard MarkDuplicates \
            I=~{input_bam} \
            O=~{sample_name}_dedup.bam \
            M=~{sample_name}_metrics.txt
    >>>

    output {
        File dedup_bam = "~{sample_name}_dedup.bam"
        File metrics = "~{sample_name}_metrics.txt"
    }
}

task HaplotypeCaller {
    input {
        File input_bam
        File reference
        String sample_name
    }

    command <<<
        set -euo pipefail
        cp ~{reference} ref.fasta
        samtools faidx ref.fasta
        picard CreateSequenceDictionary R=ref.fasta O=ref.dict
        cp ~{input_bam} in.bam
        samtools index in.bam
        gatk HaplotypeCaller \
            -R ref.fasta \
            -I in.bam \
            -O ~{sample_name}.vcf
    >>>

    output {
        File output_vcf = "~{sample_name}.vcf"
    }
}`

	case "PacBio_HiFi_Assembly":
		wdlContent = `version 1.0

workflow PacBioHiFiAssembly {
    input {
        Array[File] hifi_reads
        String sample_name
    }

    call HifiAsmAssembly {
        input:
            reads = hifi_reads,
            sample_name = sample_name
    }

    call AssemblyQC {
        input:
            assembly = HifiAsmAssembly.primary_assembly,
            sample_name = sample_name
    }

    output {
        File primary_assembly = HifiAsmAssembly.primary_assembly
        File assembly_stats = AssemblyQC.stats
    }
}

task HifiAsmAssembly {
    input {
        Array[File] reads
        String sample_name
    }

    command <<<
        set -euo pipefail
        hifiasm -o ~{sample_name} ~{sep=' ' reads}
        awk '/^S/{print ">"$2; print $3}' ~{sample_name}.bp.p_ctg.gfa > ~{sample_name}.bp.p_ctg.fa
    >>>

    output {
        File primary_assembly = "~{sample_name}.bp.p_ctg.fa"
    }
}

task AssemblyQC {
    input {
        File assembly
        String sample_name
    }

    command <<<
        quast ~{assembly} -o ~{sample_name}_qc
    >>>

    output {
        File stats = "~{sample_name}_qc/report.txt"
    }
}`

	case "ONT_Long_Read_Variants":
		wdlContent = `version 1.0

workflow ONTLongReadVariants {
    input {
        File reads
        File reference
        String sample_name
    }

    call AlignLongReads {
        input:
            reads = reads,
            reference = reference,
            sample_name = sample_name
    }

    output {
        File aligned_bam = AlignLongReads.bam
        File flagstat = AlignLongReads.flagstat
    }
}

task AlignLongReads {
    input {
        File reads
        File reference
        String sample_name
    }

    command <<<
        set -euo pipefail
        minimap2 -ax map-ont ~{reference} ~{reads} | samtools sort -o ~{sample_name}.bam
        samtools flagstat ~{sample_name}.bam > ~{sample_name}.flagstat.txt
    >>>

    output {
        File bam = "~{sample_name}.bam"
        File flagstat = "~{sample_name}.flagstat.txt"
    }
}`

	default:
		return fmt.Errorf("no WDL defined for workflow %s", workflowName)
	}

	return os.WriteFile(wdlFile, []byte(wdlContent), 0600)
}

// createWDLInputs writes the Cromwell inputs JSON for workflowName.
func createWDLInputs(inputsFile, workflowName string, in workflowInputs, platform pt.SequencingPlatform) error {
	sample := fmt.Sprintf("synthetic_%s", platform)
	var inputs map[string]interface{}
	switch workflowName {
	case "GATK_Germline_Variants":
		inputs = map[string]interface{}{
			"GATKGermlineVariants.fastq_r1":         in.R1,
			"GATKGermlineVariants.fastq_r2":         in.R2,
			"GATKGermlineVariants.reference_genome": in.Reference,
			"GATKGermlineVariants.sample_name":      sample,
		}
	case "PacBio_HiFi_Assembly":
		inputs = map[string]interface{}{
			"PacBioHiFiAssembly.hifi_reads":  []string{in.R1},
			"PacBioHiFiAssembly.sample_name": sample,
		}
	case "ONT_Long_Read_Variants":
		inputs = map[string]interface{}{
			"ONTLongReadVariants.reads":       in.R1,
			"ONTLongReadVariants.reference":   in.Reference,
			"ONTLongReadVariants.sample_name": sample,
		}
	default:
		return fmt.Errorf("no inputs defined for workflow %s", workflowName)
	}
	return writeJSON(inputsFile, inputs)
}

// createSnakemakeWorkflow writes the Snakefile for workflowName and returns
// the output files (relative to the working directory) it must produce.
func createSnakemakeWorkflow(snakeFile, workflowName string, platform pt.SequencingPlatform) ([]string, error) {
	var snakemakeContent string
	var expected []string

	switch workflowName {
	case "Multi_Platform_QC":
		snakemakeContent = `rule all:
    input:
        "results/qc_report.html"

rule fastqc:
    input:
        "data/{sample}_R1.fastq",
        "data/{sample}_R2.fastq"
    output:
        "qc/{sample}_R1_fastqc.html",
        "qc/{sample}_R2_fastqc.html"
    shell:
        "mkdir -p qc && fastqc {input} -o qc/"

rule multiqc:
    input:
        expand("qc/{sample}_R{read}_fastqc.html", sample=["synthetic"], read=[1,2])
    output:
        "results/qc_report.html"
    shell:
        "multiqc qc/ -o results/ -n qc_report.html"
`
		expected = []string{"results/qc_report.html"}

	case "Population_Genomics":
		snakemakeContent = `rule all:
    input:
        "results/variant_stats.txt"

rule vcf_filter:
    input:
        "data/variants.vcf"
    output:
        "filtered/variants_filtered.vcf"
    shell:
        "bcftools filter -i 'QUAL>20' {input} > {output}"

rule variant_stats:
    input:
        "filtered/variants_filtered.vcf"
    output:
        "results/variant_stats.txt"
    shell:
        "bcftools stats {input} > {output}"
`
		expected = []string{"filtered/variants_filtered.vcf", "results/variant_stats.txt"}

	case "Structural_Variants":
		preset := "map-hifi"
		if platform == pt.ONTPromethION24 {
			preset = "map-ont"
		}
		snakemakeContent = fmt.Sprintf(`rule all:
    input:
        "results/aligned.flagstat.txt"

rule align:
    input:
        reads="data/synthetic_R1.fastq",
        ref="data/reference.fasta"
    output:
        "aligned/aligned.bam"
    shell:
        "minimap2 -ax %s {input.ref} {input.reads} | samtools sort -o {output}"

rule flagstat:
    input:
        "aligned/aligned.bam"
    output:
        "results/aligned.flagstat.txt"
    shell:
        "samtools flagstat {input} > {output}"
`, preset)
		expected = []string{"aligned/aligned.bam", "results/aligned.flagstat.txt"}

	default:
		return nil, fmt.Errorf("no Snakefile defined for workflow %s", workflowName)
	}

	return expected, os.WriteFile(snakeFile, []byte(snakemakeContent), 0600)
}

// Workflow execution functions

func runCromwellWorkflow(t *testing.T, workDir, wdlFile, inputsFile, workflowName string, platform pt.SequencingPlatform, inputBytes int64) *WorkflowResults {
	start := time.Now()

	outputs, out, err := cromwellRun(workDir, wdlFile, inputsFile)
	elapsed := time.Since(start)
	if err != nil {
		t.Logf("Cromwell execution failed: %v", err)
	}

	return &WorkflowResults{
		WorkflowName:    workflowName,
		Platform:        platform,
		Engine:          "cromwell",
		Success:         err == nil,
		ExecutionTime:   elapsed,
		OutputDirectory: workDir,
		OutputFiles:     outputs,
		EngineOutput:    out,
		Metrics: WorkflowQualityMetrics{
			DataIntegrity:   outputIntegrity(outputs),
			ProcessingSpeed: processingSpeed(inputBytes, elapsed),
		},
	}
}

func runNextflowWorkflow(t *testing.T, workDir string, workflow NextflowWorkflow, platform pt.SequencingPlatform, inputs workflowInputs) *WorkflowResults {
	start := time.Now()
	outDir := filepath.Join(workDir, "results")

	var out string
	var err error
	var outputs []string
	var integrity float64

	if workflow.Pipeline == "" {
		outputs, err = runReadCountWorkflow(workDir, "nextflow", inputs)
		if err == nil && len(outputs) > 0 && readCountMatches(outputs[0], inputs.Reads) {
			integrity = 1.0
		}
	} else {
		args := []string{"run", workflow.Pipeline, "-r", workflow.Version, "-profile", "docker",
			"--input", filepath.Join(workDir, "samplesheet.csv"), "--outdir", outDir}
		if workflow.Pipeline == "nf-core/sarek" {
			args = append(args, "--igenomes_ignore", "--fasta", inputs.Reference)
		}
		out, err = runEngine(workDir, nextflowEnv, "nextflow", args...)
		// nf-core pipelines always write pipeline_info; require it plus at
		// least one other published result.
		entries, _ := os.ReadDir(outDir)
		for _, e := range entries {
			outputs = append(outputs, filepath.Join(outDir, e.Name()))
		}
		if len(entries) >= 2 {
			integrity = 1.0
		}
	}
	elapsed := time.Since(start)
	if err != nil {
		t.Logf("Nextflow execution failed: %v", err)
	}

	return &WorkflowResults{
		WorkflowName:    workflow.Name,
		Platform:        platform,
		Engine:          "nextflow",
		Success:         err == nil,
		ExecutionTime:   elapsed,
		OutputDirectory: outDir,
		OutputFiles:     outputs,
		EngineOutput:    out,
		Metrics: WorkflowQualityMetrics{
			DataIntegrity:   integrity,
			ProcessingSpeed: processingSpeed(inputs.Bytes, elapsed),
		},
	}
}

func runSnakemakeWorkflow(t *testing.T, workDir, snakeFile, workflowName string, platform pt.SequencingPlatform, expected []string, inputBytes int64) *WorkflowResults {
	start := time.Now()

	out, err := runEngine(workDir, nil, "snakemake", "-s", snakeFile, "--cores", "2", "--directory", workDir)
	elapsed := time.Since(start)
	if err != nil {
		t.Logf("Snakemake execution failed: %v", err)
	}

	outputs := make([]string, len(expected))
	for i, rel := range expected {
		outputs[i] = filepath.Join(workDir, rel)
	}

	return &WorkflowResults{
		WorkflowName:    workflowName,
		Platform:        platform,
		Engine:          "snakemake",
		Success:         err == nil,
		ExecutionTime:   elapsed,
		OutputDirectory: workDir,
		OutputFiles:     outputs,
		EngineOutput:    out,
		Metrics: WorkflowQualityMetrics{
			DataIntegrity:   outputIntegrity(outputs),
			ProcessingSpeed: processingSpeed(inputBytes, elapsed),
		},
	}
}

// Helper functions

// prepareNextflowInputs writes the samplesheet an nf-core pipeline expects,
// pointing at the synthetic FASTQ files. The local pipeline needs none.
func prepareNextflowInputs(workDir string, in workflowInputs, workflow NextflowWorkflow, platform pt.SequencingPlatform) error {
	var content string
	switch workflow.Pipeline {
	case "":
		return nil
	case "nf-core/sarek":
		content = fmt.Sprintf("patient,sample,lane,fastq_1,fastq_2\nsynthetic,synthetic,L001,%s,%s\n", in.R1, in.R2)
	case "nf-core/nanoseq":
		content = fmt.Sprintf("group,replicate,barcode,input_file,fasta,gtf\nsynthetic,1,,%s,%s,\n", in.R1, in.Reference)
	case "nf-core/hlatyping":
		content = fmt.Sprintf("sample,fastq_1,fastq_2,seq_type\nsynthetic,%s,%s,dna\n", in.R1, in.R2)
	default:
		content = fmt.Sprintf("sample,fastq_1,fastq_2\nsynthetic,%s,%s\n", in.R1, in.R2)
	}
	return os.WriteFile(filepath.Join(workDir, "samplesheet.csv"), []byte(content), 0600)
}

// Validation functions

func validateWorkflowResults(t *testing.T, result *WorkflowResults, workflowName string, platform pt.SequencingPlatform) {
	require.True(t, result.Success, "Workflow %s should complete successfully for %s; engine output:\n%s",
		workflowName, platform, result.EngineOutput)
	assert.Less(t, result.ExecutionTime, 2*time.Hour, "Workflow should complete within 2 hours")
	assert.Equal(t, 1.0, result.Metrics.DataIntegrity,
		"Workflow %s should produce all expected outputs %v", workflowName, result.OutputFiles)

	t.Logf("Workflow validation passed - Integrity: %.2f, Throughput: %.2f MB/s (including engine start-up)",
		result.Metrics.DataIntegrity, result.Metrics.ProcessingSpeed)
}

func generateCompatibilityReport(t *testing.T, matrix map[string]map[pt.SequencingPlatform]WorkflowCompatibilityResult) {
	t.Log("📊 WORKFLOW COMPATIBILITY MATRIX REPORT")
	t.Log(strings.Repeat("=", 80))

	for engine, platformResults := range matrix {
		if len(platformResults) == 0 {
			continue
		}
		t.Logf("\n🔧 Engine: %s", strings.ToUpper(engine))
		t.Log(strings.Repeat("-", 60))

		successCount := 0
		totalTime := time.Duration(0)

		for platform, result := range platformResults {
			status := "❌ FAIL"
			if result.Success {
				status = "✅ PASS"
				successCount++
			}

			totalTime += result.ExecutionTime

			t.Logf("  %-20s %s (%.1fs, %.1f%% integrity)",
				platform, status, result.ExecutionTime.Seconds(),
				result.QualityMetrics.DataIntegrity*100)
		}

		compatibilityRate := float64(successCount) / float64(len(platformResults)) * 100
		avgTime := totalTime / time.Duration(len(platformResults))

		t.Logf("\n  📈 Summary: %.1f%% compatibility, avg time: %.1fs",
			compatibilityRate, avgTime.Seconds())
	}

	t.Log("\n" + strings.Repeat("=", 80))
}
