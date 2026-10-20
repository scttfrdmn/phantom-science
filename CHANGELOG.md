# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- **BREAKING:** `sequence_validator` takes input paths as arguments (it read only `sequences/validation`) and exits 2 when any sequence fails (#27).
- `run_e2e_tests` and the e2e test suites generate reads with the real engine and report only measured values: read statistics, wall and CPU time, allocations, bases/s and classifier accuracies (#23).
  - PacBio, MGI and AVITI, which the engine does not model, are reported as not modeled.
  - `-report` writes `e2e_results.json` and `e2e_summary.csv`.
  - `-seed` is added, and `-real-fastq` (or `PHANTOM_REAL_FASTQ`) compares real reads with generated ones.
- `run_workflow_tests` reports, for each run: status, wall time, engine exit status, reported read count and data integrity. Its JSON and CSV reports are written (#23).
- The TUI shows the loaded configuration and measured run progress (#27):
  - platform, read length, paired-end, reads per batch, batches, seed, compression and resolved fingerprinting settings
  - files, bytes and MB/s written
  - the running version in the header
  - Left and Right select the time speedup.

### Removed
- **Subsystems that reported fabricated constants**, none of which any supported command used (#22, #26):
  - `internal/authenticity`, `internal/gpu` and `internal/monitoring`
  - in `internal/genomics`: `acceleration` (Metal/CUDA), `diseases`, `features`, `populations`, `quality`, `tools`, `training`, `validation` and `vcf`
  - the platform error-signature, platform quality-modeling and cross-platform bias-correction managers
  - kmer's `GPUPackedProcessor`, which never attached a GPU
  - the `prometheus/client_golang` dependency
- `scenario`'s ClinVar fallback (#26). It could never match, and its integrator returned made-up variants. Unmatched HPO terms now fail with a message to set `planted_variant`.
- **56 of the 59 programs under `cmd/experiments/`** (#27). Only `benchmark_memory`, `benchmark_reference` and `sequence_validator` remain, as developer tools.
- `Makefile.cuda`, `Makefile.metal` and `Makefile.streaming` (#27).
- Scripts that only fed removed code (#27): `benchmark-gpu.sh`, `compare-performance.py`, `convert_vcf_to_fasta.sh`, `create_population_dataset.sh`, `download_1000_genomes.sh`, `download_1000genomes_data.sh`, `download_full_genome.sh`, `download_validation_data.sh` and `monitor_and_train.sh`.
- `benchmark_reference -compression` (the loader ignored it) and `-warmup` (#27).
- **Constant-reporting parts of the test harnesses** (#23):
  - `AnalyzeAuthenticity`
  - the per-platform and per-population constants
  - the `RunWorkflow` stub
  - fixed benchmark accuracy, cost and throughput figures
  - the "Overall Score" and verdicts
  - the constant workflow accuracy, memory and CPU figures and the score built from them
  - tests that only read back constants or literal configs
- **From the TUI** (#27):
  - the Chemistry, PhiX and Real-time fields
  - the embedding-strategy selector
  - the sample-sheet editing mode, which trapped Tab
  - the demo sample rows
  - the Load Config, Save Config, Help, Quit, Add Sample, Import CSV and Clear All buttons, which did nothing (`q` still quits)

### Fixed
- `CreatePlatformConfig` returned configs that shared one map, so changing one changed the others (#23).
- `run_workflow_tests -dry-run` said no workflow engine was available. It now lists the planned runs (#23).
- The TUI overwrote `sequencer_model`, `flow_cell_type`, `total_batches` and `sample_name` from `--config` (#27).
- The TUI highlighted the wrong time-compression option (#27).
- The TUI printed generation errors over the screen or dropped them; they are now shown in the UI (#27).
- In the TUI, a stopped run's late message could end a newer run (#27).
- `benchmark_memory` (#27):
  - it generated period-4 sequences with only a handful of distinct k-mers
  - its unsigned memory delta could wrap
- `benchmark_reference` (#27):
  - it assumed chromosomes chr1–chr5
  - it read past chromosome ends and silently dropped the errors
  - it wrote its temporary FASTA one unbuffered byte at a time
- `sequence_validator` (#27):
  - it failed every sequence whose header had no `model_target_gc`
  - it counted N bases in GC content
  - it read extra FASTA records as bases
- `Makefile.memory` and `Makefile.reference` test targets passed package paths without `./` and ran tests that do not exist (#27).
- `tests/integration/workflow_engines_test.go` failed `go vet -tags release` (#23).

## [0.9.0] - 2026-10-19

### Added
- `signal --profile-set NAME=VALUE` overrides signal-profile parameters (dwell, stall, skip, noise, adapter, offset and model scaling) for calibration runs. Overrides are validated before anything is written and recorded in the manifest as `signal.profile_overrides` (#14).
- `dorado_identity.sh --squigulator [PATH]` simulates the same truth with squigulator and uses its median identity, minus `--tolerance` points (default 1.0), as the pass threshold. Without it, `--min-identity` applies (default 90.5%) (#14).
- `transcriptome_bambu.sh` reports `injected_skip_isoforms_min_reads` and `novel_recall_any_fraction`. `MIN_GENE_FRACTION` sets the gene-read share an injected isoform needs to count as detectable (default 0.05) (#13).

### Changed
- **BREAKING (output):** the `dna_r10.4.1_e8.2_400bps_5khz` signal profile is calibrated against real R10.4.1 data with Dorado 2.1.2 SUP v5.2.0. The new values are dwell CV 0.31 (was 0.5), skip probability 0 (was 0.004), extra white noise 0 pA (was 1.0) and drift 0 pA (was 1.0); stalls stay at 0.004. Dorado median identity rises from 88.05% to 90.79%, against 91.47% for squigulator on the same truth, and the mean dwell (13.0 samples/base) matches real reads (13.1). Signal for a given seed differs from 0.8.0 (#14).
- `dorado_identity.sh` passes signal relative to squigulator instead of at a fixed 95%, which no k-mer-table synthesizer reaches. Real reads reached 98.52% and squigulator 91.47% with the same model (#14).
- `transcriptome_bambu.sh` defaults (#13):
  - `NOVEL_RATE` is 0.2 (was 0.05). At 0.05 almost no injected isoform reaches the 5% gene share Bambu requires.
  - `NDR` is 1 (Bambu's novel-discovery-rate filter off), since injected isoforms are the only true novel transcripts. `NDR=` keeps Bambu's default.
  - Recall is scored over isoforms with at least 5 reads and 5% of their gene's reads.
- The DNA signal profile's pore-model scaling (mean 93.692, stdev 23.506) was checked against Dorado's `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` config.toml, and the profile comment now cites it (#14).

## [0.8.0] - 2026-10-18

### Added
- Fingerprint encryption (#3).
  - Quality-score and SAM fingerprints are encrypted with a key using stegano's Argon2id plus AES-256-GCM.
  - The key comes from `fingerprinting.key_file` or the `PHANTOM_FINGERPRINT_KEY` environment variable. Key files must be `chmod 600`, and keys must be text of at least 16 bytes.
  - `verify --key-file` decrypts encrypted fingerprints and reports whether each file's fingerprint was encrypted.
  - New `fingerprinting.WithKey`, `ResolveKey`, `LoadKeyFile`, `ErrKeyRequired` and `ErrWrongKey`.
- `phantom-science verify FILE...` extracts and reports each FASTQ file's fingerprint (customer, dataset, creation time, stegano version). `--customer` and `--dataset` require specific values. It exits non-zero if any file has no fingerprint, a mismatched one, or cannot be read (#2).
- `phantom-science validate FILE...` runs the validation framework on FASTQ files: the phantom config's tests, or a built-in suite of nine format, quality and composition checks (#2).
- `fingerprinting.ParseFingerprint` and `FormatFingerprint`.
- `--seed` (or `seed:` in the config) makes generation byte-reproducible; without it a seed is drawn and logged to stderr so any run can be repeated (#17).
- `generate --region chr:start-end`. `--stdout` with a positioned VCF now places reads around the variants automatically (#17).
- Config fields and flags for platform, read-length distribution, basecaller tier, ONT run metadata, output compression and summary file, and output targets and workers. Values not implemented yet fail with "not implemented yet (#N)".
- Multi-population structured coalescent (`SimulateJointPanel`) with per-population demographic presets for YRI, CEU, CHB and MXL, `NewJointPanels` for panels that share one site set, and `NewPresetPopulationPanel` (#5).
- `testing.GCSkewed` positive-control helper (#7). `E2EFramework.WriteSAM` and `testing.ErrBAMUnsupported` (#18).
- Container image `deploy/production/docker/Dockerfile` (Debian bookworm, stegano built at the go.mod version, non-root user 65532) and `.dockerignore` (#4).
- CI jobs for `go vet`, `govulncheck`, the full test suite with and without CGO, a `-tags release` build, and a Docker image build with a smoke test (#9).
- Read-length distributions for generation (#10): `lognormal` (from N50 and median, or mean and SD), `gamma` (shape and scale, or mean and SD), `empirical` (a `length<TAB>count` histogram), each with optional min/max truncation and reproducible from the seed. New package `internal/genomics/readlen`.
- ONT error and quality model (#11): `--basecaller-tier` / `error_model.tier` (`sup`, `hac`, `fast`) applies per-base substitutions, insertions and deletions (sup 0.8% / 0.8% / 1.0%), a homopolymer deletion excess that grows with run length, and per-tier qualities (modal read Q20 / Q15 / Q10, ONT mean-Q definition). New package `internal/genomics/nanopore/errmodel`.
- ONT output (#12): `platform: ont` writes MinKNOW-style headers (UUID4 read IDs, runid, read, ch, start_time, flow_cell_id, protocol_group_id, sample_id) from a per-channel timing model, plus a PycoQC-compatible `sequencing_summary.txt`. Run metadata defaults derive from the seed (2675 PromethION channels, 400 b/s).
- `output.compression: gzip` / `--compress gzip` writes `.fastq.gz`; quality-score fingerprints are embedded before compression, and `verify` and `validate` read `.gz` files (#12).
- `output.file_pattern` with `{sample_name}`, `{flowcell_id}`, `{run_id}`, `{file_number}`, `{barcode}`, `{lane}` and `{mate}` (#12).
- v2 sequencer templates (`configs/sequencers/*.yaml`) load for Oxford Nanopore, with an optional `run:` section for run settings (#12).
- `validation.WithReadLengths` and `cli.ReadLengthValidation`: generated reads are checked against the configured distribution (bounds and a 4σ mean check) (#10).
- `phantom-science transcriptome` (#13): synthetic ONT cDNA reads (PCS109/110/111/114 primers) from a transcript FASTA, FASTA + GTF, or genome FASTA + GTF (minus strand reverse-complemented). Includes a Zipf/log-normal expression model or a loaded TPM/salmon/Bambu profile, 5′ truncation (3′ capture), geometric poly(A), random strand, the basecaller error tiers, MinKNOW headers and summary, and optional gzip.
- Isoform ground truth (#13): `truth.tsv` (per read), `counts.tsv` (true TPM, expected and realised reads), and `novel_isoforms.tsv`/`.gtf` for injected exon-skip and alternative-3′-end isoforms (`--novel-rate`).
- `internal/genomics/signal` (#14): nanopore raw-signal synthesis from k-mer pore models, reimplementing squigulator's algorithm in Go. It has gamma dwell, level and white noise, baseline drift, stalls and skips, an adapter leader, poly-A/T plateaus and int16 digitisation, and is reproducible per read from the run seed and read UUID. Presets: `dna_r10.4.1_e8.2_400bps_5khz` and `rna004_130bps_4khz`.
- `scripts/fetch-pore-model.sh` (#14) downloads ONT pore models, which are not vendored (MPL-2.0; review the licence).
- `internal/genomics/slow5` (#15): SLOW5 ASCII and BLOW5 writer and reader following the SLOW5 v0.2.0 spec. Record compression is none, zlib or zstd; signal compression is svb-zd, byte-compatible with slow5lib. It supports all auxiliary types including the standard ONT fields, and a chunked writer names files `<flowcell>_<runid>_<n>.blow5` with 4000 reads each. Validated against slow5tools 1.5.0 and blue-crab (pod5 round trip).
- Output targets (#16): `--total-bases` (K/M/G/T), `--total-reads` and `--flow-cells`, with `max_output_gb`/`max_output_tb` as per-flow-cell base caps.
- `--workers N` (#16) generates files in parallel, with byte-identical output for any worker count.
- `--output s3://bucket/prefix/` (#16) spools, finalizes and multipart-uploads each file under keys that mirror the local names, with regenerable-from-seed object tags.
- `phantom-manifest.json` for every run (#16): seed, sanitized config, and the size and sha256 of every output. It never contains key material.
- `phantom-science bench` (#16) reports reads/s, bases/s and MB/s per worker count.
- `phantom-science signal` (#14) synthesizes raw nanopore signal from FASTQ or FASTA reads (plain or gzip) using a k-mer pore model and a chemistry profile.
  - Writes BLOW5 chunked per run as `<flowcell>_<runid>_<n>.blow5`. The default zlib + svb-zd encoding is pod5/blue-crab compatible.
  - Keeps the UUID read IDs and MinKNOW run metadata from `generate --platform ont` as read-group attributes and ONT auxiliary fields. Other reads get version-5 UUIDs and run metadata derived from the seed.
  - Writes `truth.fa`, the sequence each signal was synthesized from (optionally error-free via `--truth`), for `scripts/acceptance/dorado_identity.sh`.
  - Writes a manifest with `signal` and `truth` outputs and a `signal` section recording the profile, pore model and inputs (each by name and sha256).
  - Reproducible per `--seed`, with byte-identical output for any `--workers`.
- `scripts/acceptance/README.md`: a guide to the acceptance scripts, covering required tools, commands, pass criteria and the issue each one closes.
- Acceptance scripts for external tools (run by the user):
  - `scripts/acceptance/transcriptome_bambu.sh`: Pychopper → minimap2 → Bambu (#13)
  - `dorado_identity.sh` (#14)
  - `slow5_validate.sh` (#15)
  - `pycoqc.sh` (#16)

### Changed
- **BREAKING:** the `quality_scores` steganography layer requires an encryption key and refuses to run without one, unless `plaintext: true` is set. Generation checks the key before writing anything (#3).
- `Verifier.Verify` returns a `VerifyResult` (`Found`, `Fingerprint`, `Encrypted`) (#3).
- `configs/dev/config.json` sets `plaintext: true`.
- `fingerprinting.NewVerifier` returns `(Verifier, error)`: a stegano-backed verifier, or an error in builds without stegano. `Verifier.ValidateFile` is removed (#2).
- A phantom config that enables validation without listing tests runs the built-in suite (#2).
- Quality-score fingerprinting is now controlled by the `quality_scores` entry in `steganography_layers`; enabling fingerprinting without that layer only tags read IDs (#1).
- When `quality_scores` is requested, generation fails instead of writing unmarked files:
  - without the stegano library (`CGO_ENABLED=0` or release builds), before any file is written
  - when a file is too small to hold the fingerprint; the file is removed
  - for `generate --stdout`, which cannot carry the fingerprint (#1)
- **BREAKING:** read generation uses seeded math/rand/v2 PCG streams instead of crypto/rand; output formats are unchanged (#17).
- File and streamed output share one read pipeline; Illumina headers, file names and qualities are unchanged (#17).
- HMM model messages go to stderr, so `generate --stdout` emits clean FASTQ.
- `--variants` or `--region` without `--stdout` is an error instead of being ignored.
- Go toolchain pinned to go1.26.9 (`toolchain` directive): go1.24.4 had 39 reachable standard-library vulnerabilities; the `go` directive stays 1.24.4 (#9).
- `run_e2e_tests` runs an ML calibration (1000 train and 1000 holdout reads per class): null accuracy with its 4σ band around 0.5, plus a GC-skew positive control that must exceed 0.9. Failed platforms are reported and the command exits non-zero (#7).
- `stream_vcf -target` defaults to 50.0, the neutral point of its k-mer diversity heuristic (#7).
- The e2e benchmark enforces absolute throughput, latency and CPU thresholds only with `PHANTOM_PERF_STRICT=1`, since they depend on the machine and its load (#7).
- `run_workflow_tests` and the workflow integration tests run Cromwell, Nextflow and Snakemake for real and pass a workflow only when its outputs check out; missing engines are reported as skipped (#18).
- `deploy/production/README.md` documents only the container image and real commands (#4).
- **BREAKING (paths):** 59 one-off programs moved from `cmd/<name>` to `cmd/experiments/<name>`; build them with `go build ./cmd/experiments/<name>`. Only the `phantom-science` CLI, `cmd/run_e2e_tests` and `cmd/run_workflow_tests` stay where they were (#19).
- **BREAKING:** a YAML config with no recognised settings is now an error instead of loading as all zeros; flat YAML configs use the JSON key names; non-ONT v2 templates (PacBio, Element, MGI) fail with an explicit error (#12).
- `phantom-science validate` reports variable read lengths (count, min, max, mean, N50) instead of failing on them; after generation, lengths are checked against the configured distribution (#10).
- With a read-length distribution, an indel variant changes the read length by exactly the indel size (#10).
- `ont_promethion24.yaml`: 2675 channels, 400 b/s, and error rates and SUP accuracy (0.974) consistent with the generator's tier table (#12).
- With an HMM model, each output file uses its own model instance, so a file generated alone matches the same file in a full run (#16).
- `Generate` writes through a hidden spool directory in the output directory and finishes by writing the manifest (#16).
- Added the github.com/klauspost/compress v1.19.2 dependency, for zstd in BLOW5 (#15).
- Documentation audit (#20):
  - The README, `docs/`, CLI help, config files, the Cromwell workflow and the Go doc comments now describe what the code does.
  - Compliance, "production-ready", "indistinguishable" and unbenchmarked performance claims are removed. Stubs and heuristics are named as such.
  - Remaining work is tracked in #22–#27.
- Sample configs (#20):
  - `configs/phantom/v1.0/enterprise-settings.yaml` keeps only the keys the code reads, with valid validation test names. It previously listed `hardy_weinberg_equilibrium`, which the framework rejects.
  - `configs/dev/*.json` drop keys the loader ignores.
  - The `sequencers/*.yaml` templates drop the descriptive `steganography`, `workflows` and `validation` sections.
- Output strings that presented fixed values as results are relabelled, for example the benchmark report's fixed per-platform figures and the TUI's display-only embedding strategy (#20).
- `examples/enterprise_fingerprinting.go` demonstrates encrypted, per-file fingerprinting with a key and verification, and its tests cover the failure cases (#20).
- **BREAKING:** codon-usage ORF and `CodonChange` coordinates are forward-strand, and `Strand` is `"+"`/`"-"` (was `"0"`/`"1"`) (#20).
- An unknown chromosome in VCF genomic context is now an error instead of returning generated bases (#20).
- **BREAKING:** validation results never claim regulatory readiness (#20).
  - `regulatory_ready` is always false. When regulatory settings are present, the summary says regulatory readiness was not assessed.
  - `academic_ready` requires every test to pass and `peer_review_ready`.
  - `compliance_level` is `none`, `tests_passed` or `academic`.
- `scripts/release.sh X.Y.Z` prepares a release for this changelog layout (#20).
  - It moves the [Unreleased] entries, updates the compare links, VERSION and `cmd/root.go`, and runs `validate-changelog` and `check-version`.
  - `--dry-run` prints the diff. It commits, tags and pushes only with `--commit`, `--tag` and `--push`.
- Genomics-tool integration tests skip when a tool is missing and fail when an installed tool fails. They no longer substitute mock results (#20).
- `make validate-configs` uses a repository `.yamllint`, counts every error instead of stopping at the first, and loads the shipped sequencer configs (#20).

### Removed
- Unused "business" Prometheus metrics in `internal/monitoring` (`phantom_customers_active`, `phantom_datasets_processed_total`, `phantom_revenue_generated_total`, `phantom_sla_compliance_percent`, `phantom_customer_satisfaction_score`). Nothing recorded values to them.
- Rust pre-commit hooks and Makefile steps for the removed `crypto-engine/` crate.
- `.github/workflows/enterprise-ci.yml`, which could not pass (#9), and `deploy/production/docker/Dockerfile.enterprise`, which could not build (#4).
- The retired "~50.3%" accuracy target and its grading (#7).
- "mock output" workflow results, always-true engine availability, forced `success = true` and hardcoded workflow metrics (#18).
- `docs/TUI_DESIGN.md`, which described an unbuilt design (#20).
- Configs that failed to load or silently loaded as empty: `configs/sequencer/novaseq-6000.json`, `configs/phantom/enterprise-settings.json`, `configs/phantom/volume-testing.json`, `configs/phantom/v1.0/volume-testing.yaml` and `sample_config.json` (#20).
- The unread `-validate` flag of `cmd/experiments/simulate_disease` (#20).

### Fixed
- `--command verify` printed "Verification completed" without checking anything, and `NewVerifier` returned a mock that reported every file as verified (#2).
- Validation never saw any data:
  - `--command validate` and generation's built-in validation passed `nil` to the validation framework. With no tests configured this reported success, and with tests configured every test failed with "invalid data type".
  - `--command validate` also returned success when validation failed.
  - Generation now validates the files it wrote, and running zero tests no longer counts as a pass (#2).
- Errors from subcommands were printed twice.
- Generated FASTQ files contained no steganographic fingerprint (#1):
  - The generator passed each read's quality string to the fingerprinting engine, which rejected it as not FASTQ and far too small, and the error was discarded.
  - Each finished file is now fingerprinted once with stegano's streaming `EmbedFile`, with one fingerprint per run, and the file keeps its permissions.
  - New `fingerprinting.FileFingerprinter` interface (`EmbedInFASTQFile`, `ExtractFromFASTQFile`).
- Fingerprinting settings never reached the generator from a phantom config (`-p`) or from a sequencer config's `enterprise.fingerprinting` section; only a top-level `fingerprinting` section in the sequencer config was read. `SequencerConfig.ResolveFingerprinting` now applies them for the root command, `generate` and the TUI. `generate` previously ignored `-p` entirely (#1).
- FASTQ write and close errors were ignored.
- `make build` failed with "cannot write multiple packages to non-directory"; it now builds `bin/phantom-science` from the main package. `make build-all` now runs `make stegano` and then `make build`.
- The README's install instructions pointed at a release download that does not exist. They now use `CGO_ENABLED=0 go install` for a pure-Go build, or `make stegano && make build` for fingerprinting.
- Variant allele draws are exact: Freq=1.0 gave REF in 1 of 256 reads, and every VAF was quantised to 1/256 (#17).
- Several variants in one read are applied right to left, so an earlier indel no longer shifts later variants; overlapping variants are skipped and counted (#17).
- Reads with deletions are extended from generated bases instead of padded with `N`, in the engine and in clinical haplotype reads (#17).
- `generate --stdout --variants` injects VCF variants at their coordinates instead of at the first REF match.
- Clinical scenarios are reproducible from the seed, and homozygous variants have VAF exactly 1.0 (#17).
- The `generate` subcommand did not validate its output when the phantom config enabled validation (#2).
- `NewPopulationPanel` ignored the population name, so every population was the same model. Panels now use per-population presets, and `LoadFromSimulation` builds variants, frequencies and genotypes from one seeded simulation (#5).
- Coalescent demography: bottlenecks ignored their duration, expansions were applied backwards in time, there was no exponential growth, and coalescence used a haploid rate (E[T2] = N instead of 2N). Tajima's D used the wrong variance (#5).
- The training pipeline uses the presets instead of hardcoded events older than the coalescent's TMRCA (#5).
- The cross-population F_ST test compared panels with no shared sites (≈0.33 for any pair). It now uses joint panels and Weir–Cockerham F_ST, checks literature ranges and a same-population null (#6).
- `WriteBAM` wrote a two-line SAM header stub; it now writes a real BAM via `samtools view -b`, or returns `ErrBAMUnsupported` (#18).
- `run_e2e_tests` silently skipped platforms whose generation or training failed (#7).
- `stream_vcf` printed "NOT ACHIEVED" in its summary regardless of the result (#7).
- `go vet` lock-copy findings in `internal/gpu` (#8).
- `HomopolymerErrorModel` chose error types in random map order, so results were not reproducible from the seed (#11).
- `HomopolymerErrorModel.SimulateErrors` applied each edit at the original position in the already-edited sequence, so after an indel later edits landed in the wrong place and long reads panicked with an index out of range (#21).
- `scripts/release.sh` built the nonexistent `./cmd/phantom-science`; `Makefile.memory`, `Makefile.reference` and `Makefile.streaming` passed `cmd/<name>` without `./`, which `go build` reads as a standard-library path (#19).
- `.gitignore` entries matched directories as well as root binaries.
- `read_length_consistency` counted a trailing CR in CRLF files (#10).
- The Metal GPU backend's Objective-C compiled on Linux CGO builds; it is now macOS-only (#9).
- In paired-end file output with a read-length distribution, R1 and R2 mates had different lengths; they now share the sampled length, and fixed-length output is unchanged (#16).
- Coalescent θ_π was half of Tajima's estimator (Σf(1−f) instead of Σ2f(1−f)). This biased ThetaPi, Tajima's D (mean −2 under neutrality) and Fay–Wu H low. Site-frequency binning also put sites in the wrong bin (truncation instead of rounding) (#20).
- Mutational-signature contributions were always 0, because substitution keys ("ACA>ATA") were compared with context keys ("ACA"). Signatures are now compared by reference-context distribution on the pyrimidine strand (#20).
- Codon bias (#20):
  - The default codon-usage table was never loaded, giving NaN probabilities and a CAI of 0.
  - Reverse-strand ORFs were edited at reverse-complement offsets on the forward sequence. Edits are now synonymous on both strands.
  - Seeded runs are now reproducible.
- `getDnDsInformation` returned a different gene class and dN/dS per call, depending on map iteration order (#20).
- `SimpleReferenceGenome` matched FASTA headers by substring, so `chr1` also collected `chr10`–`chr19` (#20).
- `ReferenceManager.OptimizeCache` could deadlock against a concurrent writer (recursive RLock). The request and cache counters are now race-free (#20).
- `LazyReferenceLoader` returned generated bases instead of reading the reference. Regions are now read from the FASTA (using `.fai` when present), so chromosomes with equal-length names no longer share sequence (#20).
- The Metal `extract_kmers_small_k` kernel declared 256 KB of threadgroup memory, above the 32 KB limit. Its pipeline failed to build, and with it all Metal setup. It now uses 16 KB, and GPU counts match the CPU for k=3–8 (#20).
- `cmd/experiments` (#20):
  - fixed a data race in `population_specific_training`
  - fixed a panic in `quick_validation` with fewer than five chromosomes
  - `sequence_classifier` now returns a deterministic population instead of a random one
- The Cromwell workflow's validation task hit SIGPIPE under `pipefail`, and `continueOnReturnCode 141` hid it, so the report stopped after the first file (#20).
- `generate --stdout` with an inaccessible `hmm_model_path` wrote a warning line into the FASTQ stream. Config-loader diagnostics now go to stderr (#20).
- `make version` did nothing on macOS, because the VERSION file satisfied the target. All Makefile targets are now `.PHONY` (#20).
- `scripts/ci-status-monitor.sh` printed "null" as every workflow name (#20).
- e2e `calculateVariance` returned the standard deviation (#20).

### Security
- Path validation (`ValidateFilePath`, `SecureWriteFile`, `SecureFileCreate`, `SecureMkdirAll`) checks containment by path component on symlink-resolved paths. Base `/data` no longer admits `/data2/x`, and a symlink inside the base cannot escape it (#20).
- Ecosystem validator reports no longer contain cloud credentials: `CloudConfig.Credentials` is never serialized (#20).
- `GenerateAuthenticFASTQFile` fails closed. If the fingerprint cannot be generated or embedded, the error is returned and the file is removed. Quality-score fingerprints are now actually embedded, one per file, when the `quality_scores` layer is configured; previously each read's quality string was rejected as too small (#20).

## [0.7.0] - 2026-10-17

### Changed
- **BREAKING: Go module path is now `github.com/scttfrdmn/phantom-science`** (was `github.com/phantom-science/phantom-science`), matching the repository so that `go install` and `go get` work.
- **Relicensed under the Apache License 2.0** (see `LICENSE` and `NOTICE`); copyright notices now name Scott Friedman.
- **BREAKING: stegano upgraded to v0.8.0**, consumed through the official Go module `github.com/scttfrdmn/stegano/go`. stegano 0.6.0 changed the payload wire format (unchanged since), so fingerprints embedded with earlier versions can no longer be extracted.
- `internal/fingerprinting` now has a single `SteganoEngine`. It replaces the hand-written CGO bridge (which still targeted the stegano v0.1.1 C API) and the six version-labelled engines (`SteganoV0_1_1` … `SteganoV0_5_1`, including a nonexistent "v0.34.0"), along with their config types and the `NewWithConfig` / `NewV0_*` constructors. `New`, `NewLegacy`, `NewVerifier`, `Engine` and `SAMFingerprintingEngine` are unchanged.
- SAM fingerprints use stegano's `auxiliary` strategy (a `zs:Z:` tag) instead of quality scores, so sequences and qualities in SAM output are untouched.
- The vendored `stegano/` source tree (0.3.2) is removed. `make stegano` (`scripts/build-stegano.sh`) builds the library version pinned in `go.mod` into `.stegano/`, and make targets set `CGO_LDFLAGS` to link it. CGO builds previously failed to link.

### Removed
- Business, strategy, progress-report and planning documents, and all of `docs/` except `SCIENTIFIC_VALIDATION.md`, `TUI_DESIGN.md`, `VERSIONING_AND_CI.md` and `development/DEVELOPMENT.md`.
- Committed build outputs and generated data:
  - 23 compiled binaries at the repository root, plus `examples/examples`
  - the `phantom_science_env/` Python virtualenv
  - benchmark, workflow and validation reports
  - the expert-validation packages
- Redistributed third-party ADMIXTURE binaries (`genomics_tools/dist`); `install_genomics_tools.sh` installs them.
- The unused `crypto-engine/` Rust crate, superseded by stegano, and obsolete bootstrap and monitoring scripts.

### Fixed
- `EmbedSAMRecords` leaked a C string on every call.
- `RunBenchmark` reported CPU time *available* to the process (GOMAXPROCS × wall time) as CPU usage, giving values of 600–2,200%. It now reports process CPU time as a share of machine capacity, from 0 to 100%.
- The `VERSION` file said 0.1.0 while the binary reported 0.6.3, so `make sync-version` would have downgraded the version.
- Test suite:
  - The end-to-end indistinguishability test compared synthetic data with a copy of itself. It now runs a null calibration on independent draws plus a positive control.
  - A GPU backend test failed on every non-macOS CI run.
  - A population test was flaky and tested nothing.
  - The stegano engine tests did not compile without CGO.

### Security
- Sample `license_key` values in `configs/phantom/` are replaced with `REPLACE_WITH_LICENSE_KEY` placeholders.

## [0.6.3] - 2026-03-12

### Added
- **`scenario` subcommand** (`phantom-science scenario --spec <yaml> [--population CEU]`): end-to-end clinical scenario generation that maps HPO phenotype terms to a pathogenic variant, plants it at the exact genomic locus with the expected VAF, and writes a machine-readable ground truth manifest
- **`internal/genomics/clinical` package**: full clinical scenario pipeline
  - `ClinicalScenario` YAML types with `LoadScenario`
  - `GroundTruthManifest` JSON types with `WriteManifest`
  - `ResolveHPO` — deterministic catalog (GLUT1/SLC2A1, PKU/PAH, CF/CFTR) with ClinVar fallback
  - `TrioHaplotypes` + `GenerateTrioHaplotypes` — Mendelian recombination via Poisson crossover model
  - `ScenarioGenerator.Run` — orchestrates HPO resolution, 50 kb windowed region, per-member FASTQ generation, and manifest writing
- **Coordinate-aware variant injection** in `GenerateToWriter`: when `GenerateOpts.RefRegion` is set, each read is assigned a random genomic `readStart` and `injectAtPosition` forces the variant locus to REF or ALT (exact VAF = Freq), replacing the legacy `bytes.Index` scatter approach
- **`VariantRecord.Chromosome` / `VariantRecord.Position`** (1-based VCF): preserved through `loadVariantsVCF` and used by the coordinate-aware injection path; legacy position-blind path unchanged
- **`GenomicRegion`** struct: half-open 0-based `[Start, End)` interval; read start positions encoded in FASTQ tile-coordinate fields for downstream recovery
- **`Engine.GenerateRawSequence`**: exported wrapper enabling external packages to build custom per-read loops (e.g. haplotype-aware generation) on top of the HMM/crypto-rand backend
- **Statistical correctness tests** (`vaf_test.go`): position accuracy (VAF ≈ 0.5 ± 0.07), no false-positive injections, de novo trio isolation, Mendelian crossover statistics, binomial VAF distribution across 20 independent trials

## [0.6.2] - 2026-03-07

### Fixed
- **Wrong standard deviation** in `cmd/test_population_structure_validation`: `diversityStd` was computed as `variance²` instead of `sqrt(variance)`, making the coefficient of variation ~125× too small (#129)
- **Division by zero** in `authenticity_processor.go` `GetAuthenticityResults`: k-mer complexity guard checked `len(kmerProfiles) > 0` but not `totalKmers > 0`; profiles with all-zero `TotalCount` caused a panic (#127)
- **Division by zero** in `cmd/validate_authenticity` `calculateBiologicalMetrics`: `GCContent` and `ATContent` divided by `totalBases` with no zero guard; panicked on nil or all-empty sequences (#128)
- **TUI tab navigation** in `internal/tui/model.go`: Tab/Shift+Tab used modulo 5 for 6 sections, making `SectionActions` unreachable via keyboard (#126)
- **FASTA parser panic** in `fasta_parser.go`: `strings.Fields(line[1:])[0]` panicked on empty or whitespace-only FASTA headers; byte-offset tracking also advanced incorrectly on skip (#125)
- **Nucleotide diversity ranges inverted** in `bridge.go`: YRI and CEU ranges were swapped, contradicting population genetics literature (YRI > CHB > CEU per 1000 Genomes Phase 3) (#122)
- **Division by zero** in `indistinguishability_tester.go` `calculateSequenceStatistics`: mean length and GC content divided by zero when sequences slice or total length was empty (#123)
- **Division by zero** in `homopolymer_errors.go` `calculateErrorStatistics`: `AvgRunLength` divided by zero when `Errors` slice was empty (#124)

## [0.6.1] - 2026-03-06

### Added
- **`generate` subcommand** (`phantom-science generate --stdout --variants <vcf> --sample-id <id>`): streams FASTQ to stdout for zero-local-storage pipelines; `--variants` injects controlled variants at specified allele frequencies; `--sample-id` appends `RG:Z:<id>` to every read header for Queryabl BAMQ integration (#115)

### Fixed
- **Integer overflow** in `cpg_depletion.go`: three sites (`extendCpGIsland`, `calculateCpGRatio`, `AnalyzeSequence`) promoted `gCount*cCount` to `float64` before multiplication, preventing silent overflow on 32-bit platforms (#119, #121)
- **Division by zero** in `homopolymer_errors.go` `calculateErrorStatistics`: `ErrorRate` now guarded when `OriginalSequence` is empty (#120)
- **Division by zero** in `advanced_phylogenetic.go` `analyzeSubstitutionPatterns`: pure-transition sequences (no transversions) no longer produce `+Inf` tsRatio (#117)
- **Division by zero** in `platform_error_signatures.go` `generateLinearDecay`: `length==1` now returns `{startRate}` instead of `NaN` (#118)
- **Inf/NaN noise in DPGAN** via `privacy_preserving_genetics.go` `AllocateBudget`: `DifferentialPrivacyEpsilon ≤ 0` now returns an error before dividing sensitivity by epsilon (#116)
- **Context cancellation gap** in `timing/simulator.go` paired-end R2: ticker receive now wrapped in `select` to honour `ctx.Done()` (#114)
- **False validation pass** in `cohort_realism.go` `validateCohortRealism`: empty patient slice no longer returns a passing score (#113)
- **TajimaD/SFS guards**: `ComputeSiteFrequencySpectrum` with `sampleSize < 2` and TajimaD variance guard for small sample sizes (#111, #112)
- **Statistical formula fixes**: sample std dev Bessel correction, quality score standard error formula (#109, #110)
- **Silent JSON marshal errors** in publication workflow; wrong quality SE formula (#107, #108)
- **Validation always-pass**: `AggregateResults` now requires actual evidence before passing (#106)
- **NaN in frequency tables** and negative `poissonRandom` return (#104, #105)
- **Homopolymer `+Inf`** division sites and `AfricanPopulationEnhancer` RNG seed (#100, #101)
- **Double Phred+33 encoding** and `calculateGCContent` divide-by-zero (#98, #99)
- **Credential fields**: additional `json:"-" yaml:"-"` tags across monitoring, expert orchestrator, and training configs (#95–#97)
- **Bessel's correction** applied consistently in `assessSignatureConsistency` and `calculateVariance` (#93–#94)
- **Real implementations** replacing random-noise stubs in adversarial defense, GAN discriminator, fingerprint obfuscation, and advanced validation (#82–#92)

---

## [0.6.0] - 2026-03-01

### Added
- **PubMed/ORCID/CrossRef/InstitutionAPI HTTP clients** (`internal/genomics/validation`): real network calls with httptest-backed test suite
- **Real system metric collection**: platform-specific files for Linux (`/proc/stat`, `/proc/meminfo`, `/proc/loadavg`, `syscall.Statfs`), Darwin (`unix.SysctlRaw` + `vm_stat` + `sysctl vm.loadavg` + `unix.Statfs`), and stub zeros for other platforms
- **SMTP email sending** in `ExpertValidationOrchestrator`
- **Training validation methods**: JS-divergence, timing, memory, and scalability checks
- **Goroutine lifecycle management**: `ctx`/`wg`/`Shutdown()` in `ExpertValidationOrchestrator`

### Fixed
- Credential fields now carry `json:"-" yaml:"-"` tags (`SMTPPassword`, `PubMedAPIKey`, `ORCIDAPIKey`, `LicenseKey`, `SMTPConfig.Password`) to prevent accidental serialization
- `io/ioutil` deprecation: replaced with `os.*`/`io.*` equivalents across validation and diseases packages
- Directory creation uses `security.SecureMkdirAll` to prevent path-traversal in user-supplied paths

---

## [0.5.0] - 2026-02-28

### Added
- **`train` subcommand** (`phantom-science train --vcf <file> --out <model.hmm.json>`): streams VCF in 100 KB adaptive chunks or accumulates sequences for kmer mode
- **Engine HMM integration**: `Engine.hmmModel` field; `resolveModelPath()` (explicit > population+dir > empty); silent `crypto/rand` fallback on load failure; HMM-sourced Phred quality scores (0–41)
- **CLI flags**: `--hmm-model-path`, `--population`, `--model-dir` as persistent root flags overriding config-file values
- **TUI** exposes `hmmModelPath` and `population` fields; `updateConfigFromUI()` returns updated model
- **Training pipeline DI hooks**: `genomeLoader`, `dataLoader`, `paramsLoader` in `HMMTrainer`/`QualityTrainer`/`PopulationTrainer`; `ModelCheckpointer` fully tested

---

## [0.4.0] - 2026-02-27

### Added
- **HMM model serialization**: `Save`/`Load` for `HMMKmerModel` and `AdaptiveHMMModel` via `ModelEnvelope` JSON envelope (`Format`, `Version`, `ModelType`, `SerializedAt`, `Population`, `Data`)
- **`LoadModel(path)`**: unified entry point that peeks `ModelEnvelope.ModelType` to dispatch to the correct loader
- **Pre-flight config validation**: `warnIfHMMModelInvalid` checks accessibility of configured HMM model paths at startup

---

## [0.3.0] - 2026-01-15

### Added
- Validation export suite: `ValidationExporter` writes JSON/CSV results for downstream analysis
- Wasserstein distance authenticity metric added to statistical validation

---

## [0.2.0] - 2025-12-01

### Added
- **Stegano v0.5.1 integration** (PR #6): ChaCha20-Poly1305 encryption, SAM fingerprint embedding via `EmbedInSAMRecords`/`ExtractFromSAMRecords`, FITS and DICOM format support, ROC/AUC validation, WebAssembly bindings
- `SAMFingerprintingEngine` interface with `Xs:Z:` auxiliary-data embedding
- Engine factory updated to try v0.5.1 first, falling back through older versions

---

## [0.1.0] - 2025-07-13

### Added
- Interactive terminal UI (Bubble Tea): sequencer selection (NovaSeq 6000, NextSeq 2000, MiSeq), sample-sheet CSV import/export, live progress monitoring, and time-compression controls (1×, 10×, 100×, 1000×).
- Sequencer run simulation: stage-based timing (initialisation, cluster generation, sequencing, base calling) with progressive FASTQ file generation, Illumina file naming, paired-end R1/R2 output and quality-score decay.
- Sequencer configurations for PacBio Revio, ONT PromethION 24, MGI DNBSEQ-T20×2 and Element AVITI alongside the Illumina models.
- Cobra-based CLI with separate sequencer and enterprise configuration files (YAML or JSON, v1.0 schema).
- Steganographic fingerprinting through the stegano library over CGO, with read-ID tagging and a pure-Go legacy fallback.
- Prometheus metrics and alerting (`internal/monitoring`) and security utilities (`internal/security`).

### Changed
- The CLI moved from the standard `flag` package to cobra.

### Security
- Sequence generation uses `crypto/rand` (CWE-338).
- FFI size conversions are checked for integer overflow (CWE-190).
- File paths are validated against traversal (CWE-22).
- HTTP servers set read-header, write and idle timeouts.
- Created files use mode 0600 and directories 0750.

## [0.1.0-mvp] - 2025-07-02

### Added
- Basic FASTQ file generation with a flag-based command-line interface and sequencer configuration files.

[Unreleased]: https://github.com/scttfrdmn/phantom-science/compare/v0.9.0...HEAD
[0.9.0]: https://github.com/scttfrdmn/phantom-science/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/scttfrdmn/phantom-science/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/scttfrdmn/phantom-science/releases/tag/v0.7.0
