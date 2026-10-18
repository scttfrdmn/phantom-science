# Phantom Science - Genomics Authenticity Engine

[![Version](https://img.shields.io/badge/version-0.7.0-blue.svg)](CHANGELOG.md)
[![CI](https://github.com/scttfrdmn/phantom-science/actions/workflows/ci.yml/badge.svg)](https://github.com/scttfrdmn/phantom-science/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

A synthetic genomic data generator: population-specific FASTQ from HMMs trained on 1000 Genomes data, clinical scenarios with known ground truth, and validation tooling for checking how closely synthetic data resembles real data. The project is pre-1.0 and under active development; see [Status](#status).

## Core Mission: Statistical Indistinguishability

The Phantom Science engine aims to generate synthetic genomic data that passes the ultimate authenticity test: when real and synthetic data are mixed, an ML classifier should achieve only ~50% accuracy (random chance), proving statistical indistinguishability.

The indistinguishability measurement framework is fully implemented (Multinomial Naive Bayes on 4-mer vectors, 5-fold cross-validation). **Measured accuracy on your data requires running the tool against real VCF input** — see [docs/SCIENTIFIC_VALIDATION.md](docs/SCIENTIFIC_VALIDATION.md) for methodology and how to run validation.

## Validation Status

| Metric | Formula | Status |
|--------|---------|--------|
| k-mer JS divergence | JS(P‖Q) = ½KL(P‖M) + ½KL(Q‖M) | ✅ Implemented |
| F_ST (Weir-Cockerham 1984) | (MSG − MSW) / (MSG + (n_c−1)·MSW) | ✅ Implemented |
| Hardy-Weinberg χ² | Σ(O−E)²/E, 1 df; p-value via erfc | ✅ Implemented |
| LD r² | D²/(p_A·q_a·p_B·q_b) | ✅ Implemented |
| ML classifier accuracy | Multinomial NB on 4-mer vectors, 5-fold CV | ✅ Implemented |
| Measured indistinguishability | ~50% target — random-chance discrimination | ⏳ Run with real VCF data |

See [docs/SCIENTIFIC_VALIDATION.md](docs/SCIENTIFIC_VALIDATION.md) for full methodology.

## Platform Support
- ✅ **7 Sequencing Platforms**: Illumina (NovaSeq, NextSeq, MiSeq), PacBio Revio, ONT PromethION, MGI DNBSEQ, Element AVITI
- ✅ **11 Global Populations**: YRI, CEU, CHB, JPT, MXL, PUR, TSI, GBR, CHS, LWK, GWD
- ✅ **Platform-Specific Authenticity**: SBS chemistry, HiFi consensus, nanopore error patterns

## Features

### 🧬 **Advanced Sequencer Platform Support**
- **Illumina Platforms**: NovaSeq 6000, NextSeq 2000, MiSeq with accurate SBS chemistry modeling
- **PacBio Revio**: HiFi circular consensus sequencing with >99.9% accuracy
- **Oxford Nanopore**: PromethION 24 with R10.4.1 chemistry and homopolymer error modeling
- **MGI DNBSEQ**: T20×2 with DNA nanoball technology
- **Element AVITI**: Novel sequencing-by-binding (SBB) chemistry

### ⚡ **Time Compression**
- **1000x**: 36-hour NovaSeq run → 2.2 minutes
- **100x**: 18-hour NextSeq run → 11 minutes
- **10x**: 65-hour MiSeq run → 6.5 hours
- **1x**: Real-time simulation for authentic experience

### 🖥️ **Interactive Terminal Interface (TUI)**
- Intuitive sequencer selection and configuration
- Real-time sample sheet management with CSV import/export
- Live progress monitoring with file generation tracking
- Configurable time compression for rapid testing

### 🔏 **Steganographic Fingerprinting**
- FASTQ quality-score and SAM tag fingerprints via stegano v0.8.0, embedded in every generated FASTQ file when the `quality_scores` layer is enabled. See [Steganographic Fingerprinting](#steganographic-fingerprinting).

## Status

Pre-1.0: interfaces and output formats may change between minor versions (see [CHANGELOG.md](CHANGELOG.md)). Known bugs and planned work, including the nanopore FASTQ and signal roadmap, are tracked in [issues](https://github.com/scttfrdmn/phantom-science/issues) and [milestones](https://github.com/scttfrdmn/phantom-science/milestones).

## Quick Start

### Installation

Requires Go 1.24 or later.

```bash
# Install without steganographic fingerprinting (pure Go, no C toolchain needed)
CGO_ENABLED=0 go install github.com/scttfrdmn/phantom-science@latest

# Or build from source with fingerprinting (needs a C toolchain and Rust/cargo)
git clone https://github.com/scttfrdmn/phantom-science.git
cd phantom-science
make stegano   # builds the stegano library pinned in go.mod into .stegano/
make build     # writes bin/phantom-science
```

With CGO enabled (the default on macOS and Linux), `go install` fails to link unless the stegano library is available; see [Building with fingerprinting](#building-with-fingerprinting). Prebuilt release binaries are not published yet.

### Launch Interactive TUI

```bash
# Start with default configuration
./phantom-science tui

# Load specific sequencer configuration
./phantom-science tui --config configs/sequencer/v1.0/novaseq-6000-s2.yaml

# Enable enterprise features
./phantom-science tui --phantom-config configs/phantom/v1.0/enterprise.yaml
```

### Command Line Usage

```bash
# Generate FASTQ files directly
./phantom-science generate --config sequencer.json [--phantom phantom.json]

# Validate FASTQ files (built-in suite, or the phantom config's validation tests)
./phantom-science validate output/*.fastq [--phantom phantom.json]

# Verify fingerprints, optionally requiring a customer and dataset
./phantom-science verify output/*.fastq --key-file fingerprint.key [--customer ACME --dataset ds-42]
```

`validate` and `verify` exit non-zero if any file fails. `verify` reports each file's customer, dataset and creation time, and needs a build with the stegano library. When a phantom config enables validation, `generate` also validates every file it writes; if no tests are listed, the built-in suite runs.

## Clinical Scenario Generation

Use the `scenario` subcommand to generate reads with a planted pathogenic variant and a machine-readable ground truth manifest:

```bash
# Resolve from HPO terms (auto-looks up SLC2A1 for GLUT1 deficiency)
phantom-science scenario \
    --spec glut1_scenario.yaml \
    --population CEU

# Scenario YAML format
cat glut1_scenario.yaml
```

```yaml
phenotype:
  hpo_terms: [HP:0001250, HP:0001253]   # seizures + elevated lactate → SLC2A1
sequencer:
  platform: NovaSeq6000
  coverage: 30
  read_length: 150
output:
  proband_fastq: proband.fastq
  ground_truth:  manifest.json
```

Or specify the variant explicitly to skip HPO resolution:

```yaml
phenotype:
  planted_variant:
    chrom: chr1
    pos: 43393429
    ref: C
    alt: T
    zygosity: heterozygous
    acmg_class: pathogenic
sequencer:
  coverage: 50
  read_length: 150
output:
  proband_fastq: proband.fastq
  ground_truth:  manifest.json
```

The generated manifest records the planted variant, expected VAF, region, read count, and wall clock time:

```json
{
  "planted_variants": [{
    "gene": "SLC2A1", "chrom": "chr1", "pos": 43393429,
    "ref": "C", "alt": "T", "zygosity": "heterozygous",
    "expected_vaf": 0.5, "acmg_class": "pathogenic"
  }],
  "sequencer": { "coverage": 30, "region_focused": "chr1:43368429-43418429" },
  "compute_metrics": { "reads_generated": 10000, "wall_clock_ns": 1234567890 }
}
```

Trio generation (de novo variant absent from both parents) is enabled with:

```yaml
trio:
  generate: true
  de_novo: true
output:
  proband_fastq: proband.fastq
  mother_fastq:  mother.fastq
  father_fastq:  father.fastq
  ground_truth:  manifest.json
```

Supported catalog entries for HPO-based resolution:

| HPO Terms | Gene | Variant | Consequence |
|-----------|------|---------|-------------|
| HP:0001250 + HP:0001253 | SLC2A1 | chr1:43393429 C>T | missense |
| HP:0001360 | PAH | chr12:102836819 G>A | missense |
| HP:0002837 | CFTR | chr7:117559590 CTT>C | frameshift |

---

## Training a Population-Specific HMM Model

Use the `train` subcommand to build a Hidden Markov Model from VCF variant data:

```bash
# Adaptive mode (default): streams VCF in chunks, low memory
phantom-science train --vcf variants.vcf.gz --out YRI.hmm.json --population YRI

# Kmer mode: accumulates all sequences, larger k-mer contexts
phantom-science train --vcf chr1.vcf.gz --vcf chr2.vcf.gz \
    --out CEU.hmm.json --model-type kmer --kmer-size 8
```

Flags:

| Flag | Default | Description |
|------|---------|-------------|
| `--vcf` | *(required)* | VCF or VCF.gz input; may be repeated |
| `--out` | *(required)* | Output path for `.hmm.json` model file |
| `--model-type` | `adaptive` | `adaptive` (streaming) or `kmer` (in-memory) |
| `--kmer-size` | `4` | K-mer size (2–8) |
| `--population` | | 1000 Genomes population code stored in model envelope |

## Generating with a Trained Model

Configure via `sequencer.yaml`:

```yaml
hmm_model_path: "/models/YRI.hmm.json"   # explicit path, OR:
population: "YRI"                          # resolved from model_dir
model_dir:  "/models"
```

Or via CLI flags:

```bash
# Load model from explicit path
phantom-science --hmm-model-path /models/YRI.hmm.json --config seq.yaml --command generate

# Load model by population code from a directory
phantom-science --population YRI --model-dir /models --config seq.yaml --command generate
```

When the model path is invalid or absent, the engine silently falls back to `crypto/rand`-based
sequence generation, so existing workflows continue to function without modification.

## TUI Interface

The Terminal User Interface provides an intuitive way to configure and run sequencer simulations:

```
┌─────────────────── Phantom Science Sequencer Simulator v0.6.3 ────────────────────┐
│                                                                                    │
│ ┌─── Sequencer Selection ───┐  ┌─── Run Status ────┐  ┌─── Real-time Stats ─────┐ │
│ │ ○ NovaSeq 6000 (S2)       │  │ Status: Running    │  │ Files: 15/40           │ │
│ │ ○ NextSeq 2000 (P3)       │  │ Progress: [████▒▒] │  │ Size: 2.3 GB          │ │
│ │ ● MiSeq v3 (300bp)        │  │ ETA: 8:23          │  │ Rate: 50.2 MB/s       │ │
│ └───────────────────────────┘  └────────────────────┘  │ Current: Sample_R1.fastq│ │
│                                                         └─────────────────────────┘ │
│ ┌─── Sample Sheet ─────────────────────────────────────────────────────────────┐   │
│ │ Sample_ID    Lane  Index_1    Index_2    Project      Description           │   │
│ │ Patient_001    1   ATCGATCG   CGTAGATC   Cancer_Study  Tumor sample        │   │
│ │ Patient_002    1   GCTAGCTA   TACGACGT   Cancer_Study  Normal sample       │   │
│ │ [Add Sample] [Import CSV] [Clear All]                                      │   │
│ └─────────────────────────────────────────────────────────────────────────────┘   │
└────────────────────────────────────────────────────────────────────────────────────┘
```

### Navigation
- **Tab/Shift+Tab**: Navigate between sections
- **Arrow Keys**: Navigate within sections
- **Enter**: Select/Edit
- **Space**: Toggle checkboxes
- **Ctrl+Q**: Quit

## Configuration

### Sequencer Configuration

```yaml
# configs/sequencer/v1.0/novaseq-6000-s2.yaml
sequencer:
  model: "NovaSeq 6000"
  flow_cell_type: "S2"
  
sample:
  name: "Cancer_Study_001"
  lane: 1
  output_dir: "./sequencer_output"
  
reads:
  length: 150
  reads_per_batch: 1000
  total_batches: 20
  paired_end: true
  
quality:
  base_quality: 35
  quality_decay: 0.05
  q30_percentage: 90.0

performance:
  batch_interval: "5m"
```

### Enterprise Configuration

```yaml
# configs/phantom/v1.0/enterprise.yaml
enterprise:
  customer_id: "ACME_GENOMICS_001"
  license_key: "ent-key-2024-..."
  service_tier: "enterprise"
  compliance_level: "hipaa"
  
  fingerprinting:
    enable: true
    customer_fingerprint: true
    dataset_fingerprint: true
    steganography_layers: ["quality_scores", "read_headers"]
    cryptographic_signing: true
    tamper_detection: true

validation:
  enable: true
  statistical_tests: ["hardy_weinberg", "allele_frequency"]
  population_genetics_tests: ["linkage_disequilibrium"]
  pipeline_compatibility_tests: ["gatk", "bcftools", "freebayes"]
```

## File Generation

### Output Structure
```
sequencer_output/
├── Patient_001_S1_L001_R1_001.fastq    # First batch, Read 1
├── Patient_001_S1_L001_R2_001.fastq    # First batch, Read 2
├── Patient_001_S1_L001_R1_002.fastq    # Second batch, Read 1
├── Patient_001_S1_L001_R2_002.fastq    # Second batch, Read 2
└── ...
```

### Realistic Timing Examples

| Sequencer | Flow Cell | Real Duration | Files | Total Size | 100x Compressed |
|-----------|-----------|---------------|-------|------------|-----------------|
| NovaSeq 6000 | S2 | 36 hours | 40 | ~2.5 TB | 22 minutes |
| NextSeq 2000 | P3 | 18 hours | 20 | ~400 GB | 11 minutes |
| MiSeq | v3 | 65 hours | 10 | ~15 GB | 39 minutes |

### FASTQ Format
Generated files contain realistic sequencing data:

```
@Patient_001:1:HXXXXXX:1:1101:1000:2000 1:N:0:1
TGTAATCTGATACGGTGGAGTATAATTCGATTGGGGAAAAGAAGCTGACTACTCAGAAGAGTGGCAAGGAG
+
??????????>>>>>>>>>>==========<<<<<<<<<<;;;;;;;;;;::::::::::999999999988
```

## Use Cases

### 🔬 **Research & Development**
- Test bioinformatics pipelines with realistic data streams
- Validate file handling and parsing logic
- Train new lab members on sequencer operations
- Performance testing with large datasets

### 🏭 **Production Pipeline Testing**
- **File watcher integration**: Monitor progressive file generation
- **Data transfer validation**: Test upload pipelines with realistic timing
- **Quality control pipeline development**: Test with authentic sequencer parameters
- **Resource scaling and load testing**: Use time compression for rapid iteration

### 📚 **Education & Training**
- Interactive genomics workflow education
- Sequencer operation training without hardware
- Bioinformatics course demonstrations
- Pipeline development workshops

## Steganographic Fingerprinting

Phantom Science uses [stegano](https://github.com/scttfrdmn/stegano) v0.8.0 (through its Go module, `github.com/scttfrdmn/stegano/go`) to embed a fingerprint, such as `customer:dataset:stegano-0.8.0:ts-…`, into genomics data and read it back:

- **FASTQ:** the fingerprint goes into quality-score least-significant bits. Headers and sequences are unchanged, and each quality score changes by at most one Phred step. A typical fingerprint needs about 900 quality characters (six 150 bp or eight 100 bp reads); the exact minimum grows with fingerprint length.
- **SAM:** the fingerprint goes into a `zs:Z:` optional tag. The 11 mandatory fields of each record are unchanged.
- **Extraction** fails cleanly on data with no fingerprint or a damaged one; it never returns a wrong fingerprint.

Fingerprints are **encrypted with a key** you supply: stegano derives an encryption key from it with Argon2id and seals the payload with AES-256-GCM. Without the key a fingerprint cannot be read, and in SAM output the tag no longer shows the customer and dataset IDs. Unencrypted fingerprints are possible only with an explicit `plaintext: true`; anyone with the stegano library can read those. Fingerprints written by stegano 0.5.x and earlier use an older wire format and cannot be extracted.

### Fingerprinting generated FASTQ

Add `quality_scores` to `steganography_layers` to fingerprint every FASTQ file a run writes. The settings can come from a phantom config passed with `-p` (which takes precedence, including its `customer_id`) or from the sequencer config's `enterprise.fingerprinting` section:

```json
"fingerprinting": {
  "enable": true,
  "customer_fingerprint": true,
  "dataset_fingerprint": true,
  "steganography_layers": ["quality_scores"],
  "key_file": "/path/to/fingerprint.key"
}
```

Create a key once, keep it private, and use the same key to verify:

```bash
openssl rand -base64 32 > fingerprint.key && chmod 600 fingerprint.key
```

- **Where the key comes from:** `key_file`, or the `PHANTOM_FINGERPRINT_KEY` environment variable when no file is set. Relative `key_file` paths resolve against the working directory.
- **What a valid key looks like:** text, at least 16 bytes, with no NUL bytes. Key files must not be readable by other users (`chmod 600`).
- **When there is no key:** the `quality_scores` layer refuses to run unless `plaintext: true` is set.

- Each file is fingerprinted once, after it is written, with one fingerprint per run (`customer:dataset:stegano-<version>:ts-<unix time>`). Files are processed in a streaming pass, so file size does not affect memory use, and they keep their 0600 permissions.
- Generation fails, rather than writing unmarked output, when the fingerprint cannot be embedded. That happens when the build has no stegano library (`CGO_ENABLED=0`, or a release build), or when a file is too small to hold it (see the capacity figure above). A file that fails is removed.
- `generate --stdout` cannot carry a quality-score fingerprint and refuses to run with this layer. Remove `quality_scores` to stream.
- `customer_fingerprint` and `dataset_fingerprint` control the separate read-ID tags (`:cust:<id>`, `:ds:<id>`) appended to each read header. These tags are plain text in every read; turn them off if customer or dataset IDs must not be visible.
- `phantom-science verify --key-file fingerprint.key FILE...` reads fingerprints back and reports whether each was encrypted.

### Building with fingerprinting

Fingerprinting needs CGO and the stegano shared library, built at the version pinned in `go.mod`.

```bash
make stegano   # fetches and builds libstegano into .stegano/
make build     # make targets set CGO_LDFLAGS to link .stegano/
```

Outside `make`, set `CGO_LDFLAGS="-L$PWD/.stegano/target/release -Wl,-rpath,$PWD/.stegano/target/release"`. Builds with `CGO_ENABLED=0`, and release builds, use a legacy pure-Go engine that does not produce extractable fingerprints.

## Workflow Integration

Phantom Science supports multiple workflow management systems for genomic analysis:

### Supported Workflow Engines

- **🐍 Snakemake**: Python-based workflow with active execution logs
- **🌊 Nextflow**: Groovy-based workflow with comprehensive reporting
- **🏗️ Cromwell**: WDL-based workflow with enterprise-grade execution

### Workflow Capabilities

Each workflow engine provides identical analysis capabilities:

1. **FASTQ Validation**: Quality and format compliance checking
2. **Steganographic Analysis**: Fingerprint detection and capacity analysis
3. **Performance Metrics**: Throughput and resource utilization analysis
4. **Report Generation**: Comprehensive HTML and JSON reports

### Quick Start Examples

```bash
# Snakemake workflow (existing)
snakemake --cores 4 --configfile configs/analysis.yaml

# Nextflow workflow (reports available)
nextflow run workflows/nextflow/phantom_analysis.nf

# Cromwell workflow (new)
cd workflows/cromwell && ./run_workflow.sh
```

## Architecture

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│  TUI Interface  │    │  Timing          │    │  FASTQ          │
│  (bubbletea)    │◄──►│  Simulator       │◄──►│  Generator      │
│                 │    │                  │    │                 │
└─────────────────┘    └──────────────────┘    └─────────────────┘
         │                       │                       │
         ▼                       ▼                       ▼
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│  Configuration  │    │  Steganographic  │    │  Validation     │
│  Management     │    │  Fingerprinting  │    │  Framework      │
│                 │    │  (stegano v0.8.0)│    │                 │
└─────────────────┘    └──────────────────┘    └─────────────────┘
                                │
                                ▼
                    ┌──────────────────┐
                    │  Quality Score   │
                    │  Steganography   │
                    │                  │
                    └──────────────────┘
                                │
                                ▼
         ┌─────────────────────────────────────────────────────────────┐
         │                 Workflow Integration                        │
         │  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────┐ │
         │  │  Snakemake  │  │  Nextflow   │  │  Cromwell (WDL)     │ │
         │  │  (Python)   │  │  (Groovy)   │  │  (Java)             │ │
         │  └─────────────┘  └─────────────┘  └─────────────────────┘ │
         └─────────────────────────────────────────────────────────────┘
```

## Development

### Requirements
- Go 1.24+
- For fingerprinting: a C toolchain and Rust (cargo)
- Make (optional)

### Building
```bash
# Build the stegano library once (CGO builds link against it)
make stegano

# Standard build (make sets CGO_LDFLAGS for .stegano/)
make build

# Pure-Go build without fingerprinting
CGO_ENABLED=0 go build -o phantom-science .
```

### Testing
```bash
# Run all tests
go test ./...

# Run tests with coverage
go test -coverprofile=coverage.out ./...
go tool cover -html=coverage.out

# Run benchmarks
go test -bench=. ./internal/genomics
```

### Code Quality
```bash
# Linting
golangci-lint run

# Coverage check (requires 85% overall, 80% per file)
./scripts/check-coverage.sh

# Pre-commit hooks
pre-commit run --all-files
```

### Performance Benchmarks

```
BenchmarkFASTQGeneration-10           100    12.3 MB/s    1000 reads/file
BenchmarkLongReadGeneration-10         50    8.7 MB/s     300bp reads
BenchmarkHighVolumeGeneration-10      200    45.2 MB/s    10K reads/file
```

## Contributing

### Development Workflow
1. Fork the repository
2. Create feature branch: `git checkout -b feature/amazing-feature`
3. Make changes with tests
4. Ensure code coverage ≥85%
5. Run pre-commit hooks
6. Submit pull request

### Code Standards
- Go code follows `gofmt` and `golangci-lint` standards
- Test coverage must be ≥85% overall, ≥80% per file
- All public functions must have documentation
- Commit messages follow conventional commits

## License

Licensed under the [Apache License, Version 2.0](LICENSE). See [NOTICE](NOTICE).

Copyright 2025-2026 Scott Friedman
