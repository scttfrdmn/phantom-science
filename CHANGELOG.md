# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `phantom-science verify FILE...` extracts and reports each FASTQ file's fingerprint (customer, dataset, creation time, stegano version). `--customer` and `--dataset` require specific values. It exits non-zero if any file has no fingerprint, a mismatched one, or cannot be read (#2).
- `phantom-science validate FILE...` runs the validation framework on FASTQ files: the phantom config's tests, or a built-in suite of nine format, quality and composition checks (#2).
- `fingerprinting.ParseFingerprint` and `FormatFingerprint`.

### Changed
- `fingerprinting.NewVerifier` returns `(Verifier, error)`: a stegano-backed verifier, or an error in builds without stegano. `Verifier.ValidateFile` is removed (#2).
- A phantom config that enables validation without listing tests runs the built-in suite (#2).
- Quality-score fingerprinting is now controlled by the `quality_scores` entry in `steganography_layers`; enabling fingerprinting without that layer only tags read IDs (#1).
- When `quality_scores` is requested, generation fails instead of writing unmarked files:
  - without the stegano library (`CGO_ENABLED=0` or release builds), before any file is written
  - when a file is too small to hold the fingerprint; the file is removed
  - for `generate --stdout`, which cannot carry the fingerprint (#1)

### Removed
- Unused "business" Prometheus metrics in `internal/monitoring` (`phantom_customers_active`, `phantom_datasets_processed_total`, `phantom_revenue_generated_total`, `phantom_sla_compliance_percent`, `phantom_customer_satisfaction_score`). Nothing recorded values to them.
- Rust pre-commit hooks and Makefile steps for the removed `crypto-engine/` crate.

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

[Unreleased]: https://github.com/scttfrdmn/phantom-science/compare/v0.7.0...HEAD
[0.7.0]: https://github.com/scttfrdmn/phantom-science/releases/tag/v0.7.0
