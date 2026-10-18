# Phantom Science

[![CI](https://github.com/scttfrdmn/phantom-science/actions/workflows/ci.yml/badge.svg)](https://github.com/scttfrdmn/phantom-science/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

Phantom Science generates synthetic sequencing data for testing bioinformatics
pipelines: Illumina-style and Oxford Nanopore FASTQ, ONT cDNA reads with
isoform ground truth, and raw nanopore signal in BLOW5. Output is reproducible
from a seed, and FASTQ files from `generate` can carry an encrypted
steganographic fingerprint that identifies the customer and dataset they were
made for.

The project is pre-1.0. Interfaces and output formats can change between minor
versions; see [CHANGELOG.md](CHANGELOG.md). Features marked "since 0.8.0" are
on `main` but not in a tagged release yet.

## What it does and does not do

- **Reads are synthetic sequence, not samples of a reference genome.** Bases
  are uniform random, or drawn from an HMM trained with `train`. Reads do not
  align to GRCh38. Paired-end mates share a length and start position but are
  generated independently (there is no fragment or insert-size model).
- **Illumina-style output** has Illumina headers and file names. Qualities are
  a fixed per-position ramp (`base_quality` minus `quality_decay` per base),
  the same for every read, and no sequencing errors are added. The
  `sequencer_model`, `flow_cell_type` and `q30_percentage` config fields do not
  change the output.
- **ONT output** (since 0.8.0) has read-length distributions, a per-tier
  basecaller error and quality model, MinKNOW-style headers and a
  `sequencing_summary.txt`.
- **PacBio, MGI and Element** templates in `configs/sequencers/` fail to load
  with an explicit error. Only ONT v2 templates are supported.
- **No claim of statistical indistinguishability from real data is made.** The
  repository has a real-vs-synthetic classifier, but no real reads are bundled
  and no measurement against real data has been published. See
  [docs/SCIENTIFIC_VALIDATION.md](docs/SCIENTIFIC_VALIDATION.md).

## Install

Requires Go 1.24.4 or later. `go.mod` pins the go1.26.9 toolchain, which the
`go` command downloads automatically. Fingerprinting also needs a C toolchain
and Rust (`cargo`) to build the [stegano](https://github.com/scttfrdmn/stegano)
library.

```bash
git clone https://github.com/scttfrdmn/phantom-science.git
cd phantom-science
make stegano   # builds the stegano library pinned in go.mod into .stegano/
make build     # writes bin/phantom-science, linked against .stegano/
```

Outside `make`, set the linker flags yourself:

```bash
L="$PWD/.stegano/target/release"
CGO_LDFLAGS="-L$L -Wl,-rpath,$L" go build -o phantom-science .
```

A pure-Go build needs no C or Rust toolchain, but cannot write or read
quality-score fingerprints:

```bash
CGO_ENABLED=0 go build -o phantom-science .
```

`CGO_ENABLED=0 go install github.com/scttfrdmn/phantom-science@latest` installs
the latest tagged release, which does not have the features marked "since
0.8.0". Prebuilt binaries are not published. A container image is described in
[deploy/production/README.md](deploy/production/README.md).

The examples below run from the repository root and call the binary
`phantom-science`; put `bin/` on your `PATH` or adjust the path.

## Quick start

### Illumina-style FASTQ

A sequencer config sets the sample, read length, batch size and output
directory. `configs/sequencer/` has examples in JSON and v1.0 YAML.

```bash
cat > run.json <<'EOF'
{
  "sample_name": "Sample_001",
  "lane": 1,
  "output_dir": "./illumina_out",
  "read_length": 150,
  "reads_per_batch": 1000,
  "total_batches": 2,
  "paired_end": true,
  "base_quality": 30,
  "quality_decay": 0.1
}
EOF
phantom-science generate --config run.json --seed 1
phantom-science validate illumina_out/*.fastq
```

This writes `Sample_001_S1_L001_R1_001.fastq`, `..._R2_001.fastq` and so on
(one R1/R2 pair per batch), plus `phantom-manifest.json`. `validate` runs nine
format, quality and composition checks and exits non-zero if any file fails.

### Oxford Nanopore FASTQ, then raw signal

```bash
# Reads with MinKNOW-style headers, gzipped, plus sequencing_summary.txt.
phantom-science generate --config configs/sequencers/ont_promethion24.yaml \
    --platform ont --total-bases 5M --seed 1 --output reads

# Raw signal needs an ONT k-mer pore model, which is downloaded, not shipped.
model=$(scripts/fetch-pore-model.sh dna_r10.4.1_e8.2_400bps)
phantom-science signal --input reads/*.fastq.gz --out sig \
    --pore-model "$model" --seed 1
```

`sig/` then holds `<flowcell>_<runid>_0.blow5`, `truth.fa` (the sequence each
signal was made from) and `phantom-manifest.json`. To check the files with
slow5tools, Dorado and other external tools, see
[scripts/acceptance/README.md](scripts/acceptance/README.md).

## Commands

Run `phantom-science <command> --help` for every flag.

| Command | What it does |
|---|---|
| `generate` | Writes FASTQ files from a sequencer config, or streams FASTQ to stdout with `--stdout` |
| `validate FILE...` | Runs the validation checks on FASTQ files (plain or `.gz`) |
| `verify FILE...` | Extracts and reports each file's fingerprint (needs the stegano library) |
| `train` | Trains an HMM from VCF REF/ALT alleles and writes a `.hmm.json` model |
| `scenario` | Writes FASTQ with a planted pathogenic variant and a ground-truth manifest |
| `transcriptome` | Since 0.8.0. Writes ONT cDNA reads from transcripts, with expression and isoform truth |
| `signal` | Since 0.8.0. Synthesizes raw nanopore signal (BLOW5) from FASTQ or FASTA reads |
| `bench` | Since 0.8.0. Reports reads/s, bases/s and MB/s of `generate` for several worker counts |
| `tui` | Interactive terminal UI for Illumina-style runs (see [TUI](#tui)) |

The older form `phantom-science --config FILE --command generate` still
works; prefer the subcommands.

### generate

- **Reproducibility:** `--seed N` (or `seed:` in the config) makes output
  byte-identical across runs. Without it a seed is drawn and printed to
  stderr.
- **How much to write:** `--total-reads`, or `--total-bases 140G` (K/M/G/T
  suffixes). `--flow-cells N` writes N independent runs, each in its own
  subdirectory. Without these, `reads_per_batch × total_batches` applies.
- **Parallelism:** `--workers N` writes files in parallel. Output is
  byte-identical for any worker count.
- **Output:** `--output DIR` overrides `output_dir`. `--compress gzip` writes
  `.fastq.gz`. `--output s3://bucket/prefix/` uploads each finished file to S3
  using the default AWS credential chain. The S3 path is tested only against a
  fake uploader, not against real AWS.
- **Manifest:** every run writes `phantom-manifest.json` with the seed, the
  sanitized config, and the size and sha256 of every output. It never contains
  key material.
- **HMM model:** `--hmm-model-path model.hmm.json`, or
  `--population CEU --model-dir models/` to load `models/CEU.hmm.json`. If the model cannot be
  loaded, `generate` warns and uses uniform random bases.
- **Variants:** `--stdout --variants calls.vcf[.gz]` injects VCF variants
  with their allele frequency. `--region chr:start-end` places reads in a
  region, so positioned variants land at their coordinates. `--sample-id`
  adds `RG:Z:<id>` to each header (for `bwa mem -C`).

### Nanopore options (since 0.8.0)

| Option | Flag | Values |
|---|---|---|
| Platform | `--platform` | `illumina` (default) or `ont`. ONT reads are single-end. |
| Read lengths | `--read-length-dist` | `fixed`, `lognormal` (`--read-length-n50` and `--read-length-median`, or `--read-length-mean` and `--read-length-sd`), `gamma` (`shape` and `scale` in the config, or mean and SD), `empirical` (a `length<TAB>count` histogram named by `read_length_distribution.empirical_path`) |
| Basecaller tier | `--basecaller-tier` | `sup`, `hac`, `fast`: substitutions, insertions, deletions and homopolymer deletions, with per-tier qualities |

With `--platform ont`, read IDs are UUID4s and headers carry `runid`, `read`,
`ch`, `start_time`, `flow_cell_id`, `protocol_group_id` and `sample_id` from a
per-channel timing model (2675 channels, 400 b/s by default). A
PycoQC-style `sequencing_summary.txt` is written next to the reads.
`configs/sequencers/ont_promethion24.yaml` is a PromethION template with an
optional `run:` section.

### transcriptome (since 0.8.0)

```bash
phantom-science transcriptome --transcripts transcripts.fa \
    --reads 100000 --out cdna --seed 1
phantom-science transcriptome --genome genome.fa --gtf genes.gtf \
    --reads 100000 --kit PCS114 --tier sup --novel-rate 0.05 --out cdna
```

Inputs are a transcript FASTA, a FASTA plus GTF, or a genome FASTA plus GTF
(minus-strand transcripts are reverse-complemented). Expression is a
Zipf/log-normal model, or a TPM, salmon `quant.sf` or Bambu counts table
(`--expression-profile`). Reads have PCS109/110/111/114 primers, 5′
truncation, a poly(A) tail, random strand, the basecaller error tiers and
MinKNOW headers. Truth files: `truth.tsv` (per read), `counts.tsv` (true TPM,
expected and realised reads) and `novel_isoforms.tsv` (injected exon-skip and
alternative-3′-end isoforms), plus `novel_isoforms.gtf` when the isoforms have
genomic coordinates.

### signal (since 0.8.0)

`signal` reimplements squigulator's algorithm in Go: k-mer pore-model levels,
gamma-distributed dwell, level and white noise, baseline drift, stalls and
skips, an adapter leader, poly-A/T plateaus and int16 digitisation. Profiles
are `dna_r10.4.1_e8.2_400bps_5khz` and `rna004_130bps_4khz`. Output is BLOW5
following the SLOW5 v0.2.0 spec, 4000 reads per file by default
(`--reads-per-file`), with zlib or zstd record compression and svb-zd signal
compression. Read IDs and run metadata from `generate --platform ont` headers
are kept. Output is reproducible per `--seed` and identical for any
`--workers`. `--truth` takes the error-free sequences, so the signal is made
from the true molecule rather than the basecalled read.

Pore models come from
[nanoporetech/kmer_models](https://github.com/nanoporetech/kmer_models) via
`scripts/fetch-pore-model.sh`. They are not part of this repository; review
ONT's licence before you use or share them.

### Acceptance against external tools

These checks need external software and are run by hand, not in CI. Details
are in [scripts/acceptance/README.md](scripts/acceptance/README.md).

| Feature | Script | Status |
|---|---|---|
| BLOW5 writer (#15) | `slow5_validate.sh` | Passed 132/132 checks with slow5tools 1.5.0 and blue-crab 0.6.0, including a pod5 round trip |
| Signal synthesis (#14) | `dorado_identity.sh` | Not run yet: Dorado basecall identity is unmeasured |
| Transcriptome (#13) | `transcriptome_bambu.sh` | Not run yet: Bambu quantification is unmeasured |
| Sequencing summary (#16) | `pycoqc.sh` | Not run yet |

## Fingerprinting

Phantom Science uses [stegano](https://github.com/scttfrdmn/stegano) v0.8.0
(Go module `github.com/scttfrdmn/stegano/go`) to embed one fingerprint per
output file, `customer:dataset:stegano-<version>:ts-<unix time>`, into the
least-significant bits of the FASTQ quality scores. Headers and sequences are
unchanged, and each quality score changes by at most one Phred step
(`TestSteganoEngine_FASTQRoundTrip`). Extraction fails on files with no
fingerprint or a damaged one (`TestSteganoEngine_CleanFASTQHasNoFingerprint`,
`TestSteganoEngine_CorruptedFingerprintIsRejected`).

Fingerprints are encrypted by default: stegano derives a key from your key
with Argon2id and seals the payload with AES-256-GCM. Without the key the
fingerprint cannot be read.

1. Create a key once and keep it private:

   ```bash
   openssl rand -base64 32 > fingerprint.key && chmod 600 fingerprint.key
   ```

2. Enable the `quality_scores` layer in the sequencer config's
   `enterprise.fingerprinting` section, or in a phantom config passed with
   `-p` (which takes precedence):

   ```json
   "enterprise": {
     "customer_id": "ACME",
     "fingerprinting": {
       "enable": true,
       "steganography_layers": ["quality_scores"],
       "key_file": "fingerprint.key"
     }
   }
   ```

3. Generate as usual, then verify with the same key:

   ```bash
   phantom-science verify --key-file fingerprint.key --customer ACME out/*.fastq.gz
   ```

The rules:

- **Key source:** `fingerprinting.key_file`, or the `PHANTOM_FINGERPRINT_KEY`
  environment variable when no file is set. Relative paths resolve against
  the working directory. A key is text of at least 16 bytes with no NUL bytes,
  and a key file must not be readable by other users (`chmod 600`).
- **No key:** generation refuses to start, unless `plaintext: true` is set.
  Anyone with the stegano library can read plaintext fingerprints.
  `configs/dev/config.json` sets `plaintext: true`.
- **Fail closed:** generation fails rather than write unmarked files. That
  happens without the stegano library (`CGO_ENABLED=0` or `-tags release`
  builds, checked before any file is written), or when a file is too small to
  hold the fingerprint (the file is removed). With a short customer ID an
  encrypted fingerprint fits in six 150 bp reads or eight 100 bp reads. The
  minimum grows with the fingerprint's length.
- **gzip:** the fingerprint is embedded before compression. `verify` and
  `validate` read `.fastq.gz`.
- **Streaming:** `generate --stdout` cannot carry a fingerprint and refuses to
  run with the `quality_scores` layer.
- **Read-ID tags:** `customer_fingerprint` and `dataset_fingerprint` append
  plain-text `:cust:<id>` and `:ds:<id>` tags to every read header. They are
  separate from the encrypted fingerprint. Turn them off if the IDs must not
  be visible.
- **verify:** reports customer, dataset, creation time, stegano version and
  whether the fingerprint was encrypted. `--customer` and `--dataset` require
  specific values. It exits non-zero if any file has no fingerprint, a
  mismatched one, or cannot be read.
- **Not covered:** BLOW5 output from `signal` and reads from `transcriptome`
  and `scenario` do not carry a quality-score fingerprint. The library can also put a fingerprint
  in a SAM `zs:Z:` tag, but no command writes SAM. Fingerprints written by
  stegano 0.5.x and earlier cannot be extracted.

## Training an HMM

```bash
# Adaptive (default): streams the VCF in 100 KB chunks.
phantom-science train --vcf variants.vcf.gz --out models/YRI.hmm.json --population YRI

# k-mer model: holds all sequence in memory, k from 2 to 8.
phantom-science train --vcf chr1.vcf.gz --vcf chr2.vcf.gz \
    --out models/CEU.hmm.json --model-type kmer --kmer-size 8
```

Training reads only the REF and first ALT allele of each VCF record and
counts k-mer transitions in them. `--population` is a label stored in the
model file and used to find it with `--model-dir`. It does not change how the
model is trained. A model trained on very little data produces degenerate,
low-complexity sequence.

Separately, `internal/genomics/population` has a structured coalescent with
demographic presets for YRI, CEU, CHB and MXL (other 1000 Genomes codes map to
the closest preset). It is used by the training pipeline library and tests, not
by `generate`. See [docs/SCIENTIFIC_VALIDATION.md](docs/SCIENTIFIC_VALIDATION.md).

## Clinical scenarios

`scenario` writes reads around a pathogenic variant and a JSON manifest of
what was planted. The variant comes from HPO terms or is given explicitly:

```yaml
phenotype:
  hpo_terms: [HP:0001250, HP:0001253]   # resolves to SLC2A1 chr1:43393429 C>T
  # or: planted_variant: {chrom: chr1, pos: 43393429, ref: C, alt: T,
  #                      zygosity: heterozygous, acmg_class: pathogenic}
sequencer:
  platform: NovaSeq6000
  coverage: 30
  read_length: 150
output:
  proband_fastq: proband.fastq
  ground_truth: manifest.json
```

```bash
phantom-science scenario --spec glut1.yaml --population CEU
```

HPO resolution uses a built-in catalogue of three entries:

| HPO terms | Gene | Variant | Consequence |
|---|---|---|---|
| HP:0001250 + HP:0001253 | SLC2A1 | chr1:43393429 C>T | missense |
| HP:0001360 | PAH | chr12:102836819 G>A | missense |
| HP:0002837 | CFTR | chr7:117559590 CTT>C | frameshift |

Reads cover a 50 kb window around the variant at the given coverage, and the
manifest records the variant, expected VAF, region and read count. As with
`generate`, the surrounding sequence is synthetic. Trio generation
(`trio: {generate: true, de_novo: true}` with `mother_fastq` and
`father_fastq`) runs an msprime haplotype simulation in Python. It needs a
virtualenv at `./phantom_science_env` with `msprime`, `tskit` and `numpy`, and
must be run from the repository root. Without it the trio path fails.

## TUI

`phantom-science tui [--config FILE] [--phantom-config FILE]` opens a terminal
UI that picks NovaSeq 6000 S2, NextSeq 2000 P3 or MiSeq v3 timing and runs a
simulated run at 100× time compression. Files are written by the same engine as
`generate` as the simulated run progresses. Keys: Tab and Shift+Tab move
between panels, Up and Down choose the sequencer, Enter on the actions panel
starts or stops a run, `q` or Ctrl+C quits.

Limitations: the sample sheet shows two fixed example samples and cannot be
edited. The time compression setting, the run parameters, the output settings
(including the gzip checkbox) and the fingerprinting strategy are display-only.
Time compression exists only in the TUI. `generate` writes as fast as it can, and
`batch_interval` has no effect.

## Testing tools

- `CGO_ENABLED=0 go run ./cmd/run_e2e_tests -quick` runs an end-to-end
  harness with its own platform models. Its ML calibration trains a
  Multinomial Naive Bayes classifier on two independent draws from one
  generator (it passes if accuracy is within 0.5 ± 4σ) and on a GC-skewed
  positive control (it passes above 0.9). This tests the measurement, not
  realism. The harness also prints an "overall score" and "authenticity"
  numbers that are internal heuristics, not validated metrics.
- `CGO_ENABLED=0 go run ./cmd/run_workflow_tests` runs a minimal
  read-counting pipeline in Cromwell, Nextflow and Snakemake when they are
  installed, and reports missing engines as skipped. `workflows/cromwell/`
  has a WDL workflow ([workflows/cromwell/README.md](workflows/cromwell/README.md)).
- One-off research programs are in `cmd/experiments/`
  (`go build ./cmd/experiments/<name>`). They are not maintained as part of
  the CLI.

## Documentation

- [docs/SCIENTIFIC_VALIDATION.md](docs/SCIENTIFIC_VALIDATION.md): the statistics, the HMM and the coalescent, and what is and is not validated
- [docs/development/DEVELOPMENT.md](docs/development/DEVELOPMENT.md): building, testing and contributing
- [docs/VERSIONING_AND_CI.md](docs/VERSIONING_AND_CI.md): versioning, changelog and CI
- [scripts/acceptance/README.md](scripts/acceptance/README.md): acceptance scripts for external nanopore tools
- [deploy/production/README.md](deploy/production/README.md): container image
- [configs/README.md](configs/README.md): configuration files

## License

Licensed under the [Apache License, Version 2.0](LICENSE). See [NOTICE](NOTICE).

Copyright 2025-2026 Scott Friedman
