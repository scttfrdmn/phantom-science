# Phantom Science Cromwell workflow

`phantom_analysis.wdl` is a small WDL demonstration that reads FASTQ files
written by `phantom-science generate` and writes text, HTML and JSON reports.
It shows that the generated files can be passed through Cromwell; it is not an
analysis pipeline.

What each task does:

| Task | What it does |
|---|---|
| `BatchValidation` | checks that each file's line count is a multiple of 4 and reports the N content of the first 1000 reads; quality scores are not checked |
| `FingerprintAnalysis` | counts the quality characters (total and distinct) in the first 1000 reads. It does **not** read fingerprints, which are encrypted and need the stegano library and the key: use `phantom-science verify` |
| `PerformanceAnalysis` | reports each file's size and read count, and a processing time computed from an assumed 50 MB/s (not a measurement) |
| `GenerateReport` | combines the three reports into `phantom_analysis_report.html` and `workflow_summary.json` |

The inputs `customer_id`, `analysis_type`, `min_quality_score`,
`fingerprint_strategy`, `enable_steganography` and `batch_size` are only
echoed into the reports. `memory_limit_gb` and `cpu_cores` set the runtime of
`PerformanceAnalysis`.

## What has been run

On 2026-10-18 the workflow was run with Cromwell 92 (Homebrew), Docker and the
commands below, and finished with status `Succeeded` on a 10,000-read-pair
file pair. `womtool validate phantom_analysis.wdl` passes. `run_workflow.sh`,
the `./cromwell` wrapper it calls (which downloads Cromwell 85), and
`cromwell.conf` were not tested.

## Running it

Requirements: Java 11 or later, Cromwell, and Docker (every task runs in the
`ubuntu:22.04` image). Cromwell fails when the working directory path contains
a space.

From the repository root, generate the files that `inputs.json` names
(`novaseq_s2_output/Sample_001_S1_L001_R{1,2}_001.fastq`):

```bash
phantom-science generate --config configs/sequencer/v1.0/novaseq-6000-s2.yaml \
  --total-reads 10000 --seed 1
```

Then run the workflow from this directory, since `inputs.json` uses paths
relative to it:

```bash
cd workflows/cromwell
cromwell run phantom_analysis.wdl -i inputs.json -o options.json
```

With `options.json`, the reports are copied to
`cromwell_output/cromwell_output/{validation,fingerprinting,performance,reports}/`,
logs to `cromwell_output/logs/` and per-task logs to
`cromwell_output/call_logs/`. To check the fingerprints of the same files
(needs a stegano build and the key they were generated with):

```bash
phantom-science verify --key-file fingerprint.key novaseq_s2_output/*.fastq
```

`novaseq-6000-s2.yaml` does not enable fingerprinting, so `verify` reports
no fingerprint for these files; add
`-p configs/phantom/v1.0/enterprise-settings.yaml` to `generate` to embed one.

### run_workflow.sh (untested)

```bash
cd workflows/cromwell
./run_workflow.sh           # check, validate inputs, run, summarise
./run_workflow.sh check     # check for the ./cromwell wrapper, Java 11+, the WDL and inputs.json
./run_workflow.sh validate  # check that the FASTQ files exist
./run_workflow.sh clean     # remove the output directory
```

The script runs the repository's `./cromwell` wrapper without `cromwell.conf`.
When `jq` is installed and a FASTQ file named in `inputs.json` is missing,
`validate` writes a three-read placeholder FASTQ at `<repository root>/<path>`
instead of failing; with the default `inputs.json` that is two directories
above the repository.

## Files

```
workflows/cromwell/
├── phantom_analysis.wdl    # workflow definition
├── inputs.json             # inputs: the novaseq_s2_output file pair
├── options.json            # output and log directories, call caching
├── cromwell.conf           # Local and Docker backend configuration (not used by the commands above)
├── run_workflow.sh         # wrapper script (untested)
└── README.md
```

## Custom inputs

Edit `inputs.json` to point at other FASTQ files (plain or `.gz`):

```json
{
  "PhantomAnalysis.fastq_files": [
    "/path/to/sample1_R1.fastq",
    "/path/to/sample1_R2.fastq"
  ]
}
```

`run_workflow_tests` (in `cmd/run_workflow_tests`) runs its own minimal
read-count workflows for Cromwell, Nextflow and Snakemake; it does not use this
directory.
