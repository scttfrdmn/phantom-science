# Phantom Science

[![CI](https://github.com/scttfrdmn/phantom-science/actions/workflows/ci.yml/badge.svg)](https://github.com/scttfrdmn/phantom-science/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

Phantom Science generates synthetic sequencing data for testing bioinformatics
pipelines: Illumina-style and Oxford Nanopore FASTQ, ONT cDNA reads with
isoform ground truth, and raw nanopore signal in BLOW5. Output is reproducible
from a seed. Runs from `generate` can be fingerprinted for the customer and
dataset they were made for, without changing the data: the seed is derived
from a secret key, and `verify` attributes a file by regenerating the run.

The project is pre-1.0. Interfaces and output formats can change between minor
versions; see [CHANGELOG.md](CHANGELOG.md). Features marked "since 0.8.0" are
on `main` but not in a tagged release yet.

## What it does and does not do

- **Reads come from a reference genome, or are synthetic.** With
  `--reference genome.fa` (since 0.10.0) every read or fragment is sampled
  from the FASTA on both strands, so reads align back to it. Without it,
  bases are uniform random, or drawn from an HMM trained with `train`, and
  reads do not align to any genome. Paired-end mates are the two ends of one
  fragment, with insert sizes lognormal 455 ± 100 by default, or the run's
  own when `illumina.profile` is named, unless `illumina.insert_size` is
  set.
- **Illumina-style output** (since 0.10.0) has Illumina headers and file
  names, and a quality and error model fitted to real NovaSeq 6000 runs
  (200,000 read pairs each, aligned to GRCh38;
  `scripts/acceptance/illumina_calibrate.sh`). The default profile is
  `novaseq6000-pooled` (since 0.15.0), fitted to three runs together; see
  the profile table below.
  - qualities in the four values NovaSeq 6000 writes (Q2, Q11, Q25, Q37),
    from per-read quality classes and per-cycle fractions; generated reads
    match their runs' per-cycle fractions and mean quality (default: R1
    35.35 and R2 34.98 against 35.36 and 34.99 in the pooled reads);
  - qualities that depend on the read's bases (since 0.12.0): each read
    class's low-bin probabilities follow the 4-base context, as in the real
    run (for example, Q25 is 2.5 times as likely after TT as after CA);
  - substitutions at the measured, cycle-dependent rates of each bin (for
    DRR276626: Q37 3.2e-4; Q25 0.07–0.25%; Q11 2.9–8.3%) with the measured substitution
    spectrum (C→A and G→T dominate); Q2 bases are `N`.

  The two reads of a pair share a quality class, as real mates do (default
  profile: mates' below-Q30 counts correlate 0.56, against 0.56 in its
  pooled reads and 0.37–0.62 in single runs). Real
  runs differ a lot from each other, so the model ships run profiles
  (`illumina.profile`, `generate --illumina-profile`), each fitted to one
  run's (or, pooled, three runs') qualities, mate coupling, errors and
  insert sizes:

  | Profile | Run | Library | Mean Q R1/R2 | ≥Q30 R1/R2 | Insert |
  |---|---|---|---|---|---|
  | `novaseq6000-pooled` (default) | DRR394916 + DRR276626 + ERR15600747 (200,000 pairs each) | the three libraries | 35.4 / 35.0 | 90.7 / 89.0% | 465 ± 140 (455 ± 100 when the profile is not named) |
  | `novaseq6000-drr394916` (default in 0.13.0–0.14.0) | DRR394916 | TruSeq PCR-free | 35.3 / 34.6 | 90.3 / 86.9% | 560 ± 115 |
  | `novaseq6000-drr276626` (default before 0.13.0) | DRR276626 | not recorded | 36.2 / 36.1 | 95.4 / 94.7% | 455 ± 100 |
  | `novaseq6000-err15600747` | ERR15600747 | SureSelect QXT | 34.6 / 34.3 | 86.5 / 85.5% | 390 ± 150 |

  A single-run profile reproduces its run (pair classifier 0.55–0.58
  against it), not other runs: against three held-out runs the closest
  single-run profile scores 0.64–0.66, within the 0.58–0.76 spread of real
  runs against those runs (`scripts/acceptance/illumina_heldout.sh`). The
  pooled default is the closest to the held-out runs: pair classifier
  median 0.615, against 0.648 for drr394916, and below the real-run median
  0.675. Each profile comes from one
  tile of one run. Indels are not modelled (measured 5–7e-5 per base on all
  runs). It is not validated beyond 151 cycles
  or for other instruments (MiSeq writes unbinned qualities; NovaSeq X uses other
  values). `illumina.quality_model: linear`
  restores the old ramp (`base_quality` minus `quality_decay` per base, no
  errors). The `sequencer_model`, `flow_cell_type` and `q30_percentage` config
  fields do not change the output.
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

Requires Go 1.24.4 or later. `go.mod` pins the go1.26.10 toolchain, which the
`go` command downloads automatically. Reading quality-score fingerprints
from files made by 0.9.0 and earlier also needs a C toolchain and Rust
(`cargo`) to build the [stegano](https://github.com/scttfrdmn/stegano)
library; everything else, including fingerprinting, builds in pure Go.

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

A pure-Go build needs no C or Rust toolchain. It fingerprints and verifies
runs; it only lacks extraction of legacy (0.9.0 and earlier) quality-score
fingerprints:

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
| `verify FILE...` | Attributes files to a fingerprinted run by regenerating it (`--manifest`, or `--config` with candidates); without those flags, extracts legacy quality-score fingerprints (needs the stegano library) |
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
  key material or local paths: the files a run reads (reference FASTA, VCF,
  HMM model, read-length histogram) are recorded by base name, size and
  sha256 (`reference`, `variants`, `hmm_model`, `read_length_histogram`), not
  by path.
- **Reference genome:** `--reference genome.fa [--region chr:start-end]`
  samples reads from an uncompressed FASTA (a samtools `.fai` is used when
  present; gzip is rejected). Only the region is loaded; a whole GRCh38 takes
  about 3 GB of memory, so use `--region` for small runs. File output writes a
  `<name>.truth.tsv` per batch (read name, contig, 0-based start, strand,
  fragment length), listed in the manifest. Positioned `--variants` are
  applied in reference coordinates; a REF that does not match the reference is
  skipped with a warning, or is an error with `--strict-ref`. A reference
  cannot be combined with an HMM model. Measured on hg38
  chr22:20,000,000–21,000,000 (200,000 pairs, `minimap2 -ax sr`): 100% of reads
  mapped and properly paired, insert 348 ± 48, and per-quality mismatch rates
  equal to the model's. The manifest records the FASTA by name, size and
  sha256; `verify --reference genome.fa` regenerates such runs (also with
  `--region` and `--variants`) and refuses a FASTA whose sha256 differs.
- **HMM model:** `--hmm-model-path model.hmm.json`, or
  `--population CEU --model-dir models/` to load `models/CEU.hmm.json`. If a
  configured model is missing or cannot be loaded, `generate` fails before
  writing anything. Model files from before 0.10.0 must be retrained.
- **Variants:** `--variants calls.vcf[.gz]` (or `variants_path`) injects
  VCF variants with their allele frequency, in file output or with
  `--stdout` (file output since 0.14.0). Each read, or each pair's fragment,
  draws REF or ALT once, so both mates carry the same allele. With
  `--reference`, variants are placed at their reference coordinates and a REF
  that does not match is skipped with a warning (an error with
  `--strict-ref`, before anything is written). Without a reference,
  `--region chr:start-end` (or `variants_region`) places reads around the
  variants; for paired-end output the region must be wider than a fragment.
  File output writes `<name>.variants.tsv` per batch with one line per allele
  drawn per read or pair (read name, read number, CHROM, POS, ID, REF, ALT,
  the allele drawn, and its forward-strand offset in the read or fragment),
  and the manifest records the VCF by name, size and sha256. `verify
  --variants calls.vcf` regenerates such runs. Measured on 60,000 pairs: every
  read equals the reference plus its recorded alleles, and ALT fractions are
  within 0.7σ of the VCF's AF. `--sample-id` adds `RG:Z:<id>` to each header
  (for `bwa mem -C`).

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
| Signal synthesis (#14) | `dorado_identity.sh` | Passed: Dorado SUP median identity 90.79%, within one point of squigulator (91.47%) on the same truth. Real R10.4.1 reads reach 98.52%; the gap is attributed in [docs/SCIENTIFIC_VALIDATION.md](docs/SCIENTIFIC_VALIDATION.md) §5.3 (composition, level table, k-mer level model). RNA004: 72.9% against 97.7% for real reads |
| Transcriptome (#13) | `transcriptome_bambu.sh` | Completed on GRCh38 chr22: Bambu count correlation 0.872 (Pearson, log1p), novel-isoform precision 0.90 and recall 0.55 over detectable isoforms |
| Sequencing summary (#16) | `pycoqc.sh` | Passed: PycoQC 2.5.2 counted all 537 reads |

## Fingerprinting

Fingerprints never alter the data ([docs/CHARTER.md](docs/CHARTER.md),
principle 2). With fingerprinting on, the run's seed is derived from a secret
key and the customer and dataset IDs:

```
seed = first 8 bytes of HMAC-SHA256(key, "phantom-science/fingerprint/v1" ‖ customer ‖ dataset ‖ salt)
```

Everything else (bases, qualities, lengths, inserts, errors) follows from
that seed, so the output is byte-identical to an unfingerprinted run with the
derived seed (`TestKeyedSeed_OutputEqualsUnfingerprintedRun`) and its
qualities stay on the real bins. Nothing is embedded in the files. Without
the key, nobody can produce data that matches another customer's seed, and
the fingerprint cannot be stripped without discarding the data. It needs no
stegano library or CGO.

1. Create a key once and keep it private:

   ```bash
   openssl rand -base64 32 > fingerprint.key && chmod 600 fingerprint.key
   ```

2. Enable fingerprinting in the sequencer config's top-level or
   `enterprise.fingerprinting` section, or in a phantom config passed with
   `-p` (which takes precedence), and set the customer and dataset:

   ```json
   "dataset_id": "STUDY_1",
   "enterprise": {
     "customer_id": "ACME",
     "fingerprinting": {"enable": true, "key_file": "fingerprint.key"}
   }
   ```

   `generate --dataset-id` sets the dataset too. `--seed` is the salt; without
   it a salt is drawn and logged. The manifest records the scheme, customer,
   dataset and salt, never the key or the derived seed.

3. Verify a suspect file by regenerating candidate runs:

   ```bash
   phantom-science verify --key-file fingerprint.key \
       --manifest acme/phantom-manifest.json --manifest beta/phantom-manifest.json \
       suspect.fastq.gz
   ```

   Or name candidates without manifests: `--config run.json --customer C
   --dataset D --salt S`, or `--candidates file.tsv` (customer, dataset, salt
   per line) with `--config` or a `--manifest`.

   A run that read input files is regenerated with them: pass `--reference`,
   `--variants`, `--hmm-model` or `--read-length-histogram` (each
   repeatable). verify uses only a file whose size and sha256 match the
   manifest's record, and refuses without one. Manifests from 0.14.0 and
   earlier record these files by path; verify then uses the given file of the
   same name (or the only one given), else the recorded path, and warns that
   it cannot be checked.

How verification decides:

- **Matching:** suspect reads are matched to the regenerated run by sequence:
  whole reads on either strand, or containment of sampled 31-mers. Reads that
  were trimmed, shuffled, renamed, reverse-complemented or had every quality
  rewritten still match. A file is attributed to a candidate when at least
  three reads match and matches exceed ten times both a uniform-random chance
  bound and a `--decoy` run (the first candidate regenerated with a random
  seed).
- **Measured:** 10,000 of 10,000 reads of each file of a 200,000-pair
  NovaSeq run matched the right customer (exact), and 0 of 10,000 matched the
  wrong one. After rewriting every quality to `F`, trimming 10 bases from each
  end, renaming and shuffling, 10,000 of 10,000 still matched by 31-mer
  containment. Regenerating one candidate's 400,000 reads took about 2.5 s on
  an M4 Pro.
- **Limits:** `--max-reads` (default 10,000 suspect reads per file),
  `--max-files` and `--max-regen-reads` bound the work. The chance bound
  assumes uniform-random sequence, so for `--reference` or HMM runs use
  `--decoy`. Runs made by `signal` and `scenario` cannot be verified this way
  (`scenario` rejects fingerprinting).

The rules:

- **Key source:** `fingerprinting.key_file`, or the `PHANTOM_FINGERPRINT_KEY`
  environment variable when no file is set. Relative paths resolve against
  the working directory. A key is text of at least 16 bytes with no NUL
  bytes, and a key file must not be readable by other users (`chmod 600`).
- **Required:** a key, `enterprise.customer_id` and `dataset_id`; each
  missing one is an error before anything is written.
- **Layers:** an empty `steganography_layers`, or `keyed_seed`, turns on the
  keyed seed. `quality_scores` is rejected, because it edits qualities.
- **Read-ID labels:** `customer_fingerprint` and `dataset_fingerprint` add
  visible SAM-style comment fields `XC:Z:<customer>` and `XD:Z:<dataset>`
  (tab-separated, after the Illumina comment) to every read header. They are
  labels, not fingerprints, and are off unless set; the layer `read_ids`
  alone writes only the labels.
- **Legacy files:** `verify FILE` with no candidate flags extracts the
  quality-score fingerprints written by phantom-science 0.9.0 and earlier
  (with stegano; needs a build with the stegano library). Generation can no
  longer write them.

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
the closest preset). It is a library with tests; `generate` does not use it. See [docs/SCIENTIFIC_VALIDATION.md](docs/SCIENTIFIC_VALIDATION.md).

## Clinical scenarios

`scenario` writes reads around a pathogenic variant and a JSON manifest of
what was planted. The variant comes from HPO terms or is given explicitly:

```yaml
phenotype:
  hpo_terms: [HP:0001250, HP:0011972]   # resolves to SLC2A1 chr1:42929242 C>T
  # or: planted_variant: {chrom: chr1, pos: 42929242, ref: C, alt: T,
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
phantom-science scenario --spec glut1.yaml --population CEU --reference GRCh38.fa
```

HPO resolution uses a built-in catalogue of three entries. Each was checked on
2026-10-19 against HPO (release 2026-09-01) and ClinVar, with GRCh38
coordinates in VCF form:

| HPO terms | Condition | Variant (GRCh38) | ClinVar | Consequence |
|---|---|---|---|---|
| HP:0001250 Seizure + HP:0011972 Hypoglycorrhachia | GLUT1 deficiency | SLC2A1 c.940G>A p.Gly314Ser, chr1:42929242 C>T | VCV000016113, Pathogenic | missense |
| HP:0004923 Hyperphenylalaninemia | PKU | PAH c.1222C>T p.Arg408Trp, chr12:102840493 G>A | VCV000000577, Pathogenic | missense |
| HP:0012236 Elevated sweat chloride | Cystic fibrosis | CFTR c.1521_1523del p.Phe508del, chr7:117559590 ATCT>A | VCV000007105, Pathogenic | in-frame deletion |

Unmatched terms are an error; set `planted_variant` instead. With
`--reference`, each planted REF allele is checked against that FASTA before
any reads are written. The manifest records the build, ClinVar ID and HGVS.

Reads cover a 50 kb window around the variant at the given coverage, and the
manifest records the variant, expected VAF, region and read count. With
`--reference`, reads are sampled from the reference over the window and the
variant is planted at its coordinate; without it, bases are synthetic.
Scenario reads use the `generate` read model: binned Illumina qualities,
quality-matched errors (`illumina.errors: none` turns them off), and with
`sequencer.paired_end: true` read pairs written to `*_fastq` and
`*_fastq_r2`. Measured with hg38 (GLUT1 and a CF trio, 30×, 2×150): all reads
map and pair to the window, and the planted variant is seen at VAF 0.53
(GLUT1) and 0.43 (F508del, proband only). Trio generation
(`trio: {generate: true, de_novo: true}` with `mother_fastq` and
`father_fastq`) simulates the parents' haplotypes in Go with the structured
coalescent (`internal/genomics/population`), seeded from the scenario seed.
It needs no Python. Trio populations are YRI, CEU, CHB and MXL, plus the 1000
Genomes codes those presets model (GBR, TSI, EUR → CEU; JPT, CHS, EAS, ASN →
CHB; LWK, GWD, AFR → YRI; PUR, AMR → MXL); other codes are an error. The
region is one non-recombining locus, so crossovers come only from the genetic
map length of the region. With `--reference`, each member's haplotypes are the
reference window carrying the simulated sites (REF is the reference base;
ALT is the simulated derived base, or the ancestral base if the derived one
equals the reference), and every read comes from one phased haplotype.

## TUI

`phantom-science tui [--config FILE] [--phantom-config FILE]` opens a terminal
UI that runs a simulated sequencer run with the timing of NovaSeq 6000 S2,
NextSeq 2000 P3 or MiSeq v3, or of the model named in `--config`. Files are
written by the same engine as `generate` as the simulated run progresses, and
the run panel shows measured progress, files and bytes written. Errors appear
in the UI. Keys: Tab and Shift+Tab move between panels, Up and Down choose a
sequencer preset, Left and Right in the output panel choose the time speedup
(1×, 10×, 100× or 1000×), Enter on the actions panel starts or stops a run,
`q` or Ctrl+C quits.

The panels show the loaded configuration (platform, read length, paired-end,
reads per batch, batches, seed, compression, fingerprinting); settings and the
sample sheet are not editable in the TUI. Time compression exists only in the
TUI: `generate` writes as fast as it can, and `batch_interval` has no effect.

## Testing tools

- `CGO_ENABLED=0 go run ./cmd/run_e2e_tests -quick` generates reads with the
  engine for each modeled platform (Illumina NovaSeq, NextSeq, MiSeq and
  ONT; others are reported as not modeled) and prints measured read
  statistics and benchmarks. It calibrates classifiers on three feature sets
  (sequence 4-mers, per-read quality profile, combined) with a null (two
  generated draws, within 0.5 ± 4σ), a GC-skew control and a quality-model
  control. With `-real-fastq FILE` (or `PHANTOM_REAL_FASTQ`) it reports how
  well each feature set separates real reads from generated ones, with a
  real-versus-real null and the most separating features; `-reference
  FASTA [-region]` generates the compared reads from a genome. On NovaSeq
  6000 run DRR276626 with its own profile (`novaseq6000-drr276626`), against hg38:
  sequence 0.519, quality 0.524; against three held-out runs, quality
  0.62–0.76, about as far as those runs are from each other (chance
  0.5; see [docs/SCIENTIFIC_VALIDATION.md](docs/SCIENTIFIC_VALIDATION.md)).
  `-real-fastq2` adds read 2 and a `pair` feature set; `-real-fastq-b`
  compares two real runs, the realistic target. `-report` writes
  `e2e_results.json` and `e2e_summary.csv`.
- `CGO_ENABLED=0 go run ./cmd/run_workflow_tests` runs a minimal
  read-counting pipeline in Cromwell, Nextflow and Snakemake when they are
  installed and checks the reported count against the reads written. Missing
  engines are reported as skipped, and `-dry-run` lists the planned runs. `workflows/cromwell/`
  has a WDL workflow ([workflows/cromwell/README.md](workflows/cromwell/README.md)).
- `cmd/experiments/` holds three developer tools, not part of the CLI:
  `benchmark_memory` (k-mer counter memory and speed), `benchmark_reference`
  (reference-loader cache timings) and `sequence_validator` (FASTA GC,
  complexity and homopolymer checks).

## Documentation

- [docs/CHARTER.md](docs/CHARTER.md): the principles that decide design questions (data fidelity, fingerprints that never alter the data, measured claims)
- [docs/SCIENTIFIC_VALIDATION.md](docs/SCIENTIFIC_VALIDATION.md): the statistics, the HMM and the coalescent, and what is and is not validated
- [docs/development/DEVELOPMENT.md](docs/development/DEVELOPMENT.md): building, testing and contributing
- [docs/VERSIONING_AND_CI.md](docs/VERSIONING_AND_CI.md): versioning, changelog and CI
- [scripts/acceptance/README.md](scripts/acceptance/README.md): acceptance scripts for external nanopore tools
- [deploy/production/README.md](deploy/production/README.md): container image
- [configs/README.md](configs/README.md): configuration files

## License

Licensed under the [Apache License, Version 2.0](LICENSE). See [NOTICE](NOTICE).

Copyright 2025-2026 Scott Friedman
