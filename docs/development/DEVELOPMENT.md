# Development

## Requirements

- Go 1.24.4 or later. `go.mod` pins the go1.26.10 toolchain, which the `go`
  command downloads automatically.
- For fingerprinting (CGO builds): a C toolchain and Rust (`cargo`), to build
  the stegano library.
- Optional: `make`, `pre-commit`, `golangci-lint`, `staticcheck`, `gosec`,
  `yamllint` (`make install-tools` installs most of them).

## Building

```bash
make stegano   # builds the stegano version pinned in go.mod into .stegano/
make build     # bin/phantom-science, linked against .stegano/
```

`make` exports `CGO_LDFLAGS` for `.stegano/`. Outside make:

```bash
L="$PWD/.stegano/target/release"
CGO_LDFLAGS="-L$L -Wl,-rpath,$L" go build -o phantom-science .
```

Without the library, build with `CGO_ENABLED=0`. That build cannot embed or
extract quality-score fingerprints, and generation refuses to run when the
`quality_scores` layer is enabled.

`scripts/build-stegano.sh` builds into `.stegano/target`, not a shared
`CARGO_TARGET_DIR`, so another stegano checkout cannot overwrite the library.

### Build tags

| Tag | Effect |
|---|---|
| `cgo` (CGO_ENABLED=1) | stegano fingerprinting engine (`engine_stegano.go`); without it, `engine_nocgo.go` |
| `release` | excludes all `_test.go` files (they are `//go:build !release`) and uses the legacy fingerprinting engine, which cannot write extractable fingerprints |

There is no GPU code: the Metal and CUDA acceleration packages were removed
because nothing used them.

## Testing

CI runs the suite both ways (see [../VERSIONING_AND_CI.md](../VERSIONING_AND_CI.md)):

```bash
CGO_ENABLED=0 go test -tags '!release' ./... -timeout 15m

L="$PWD/.stegano/target/release"
CGO_LDFLAGS="-L$L -Wl,-rpath,$L" go test -tags '!release' ./... -timeout 20m
```

`-tags '!release'` is what CI passes. It sets no tag, so the `!release` test
files are included, as they are without it.

- `go test -short` shortens some slow tests, such as the coalescent presets.
- `tests/e2e` and `tests/integration` hold the cross-package tests.
- Tests that need external programs skip without them. The slow5tools
  interoperability test in `internal/genomics/slow5` runs only when
  `SLOW5TOOLS` names a slow5tools binary.
- `PHANTOM_PERF_STRICT=1` makes the e2e benchmark enforce absolute
  throughput, latency and CPU thresholds, which depend on the machine.
- `scripts/acceptance/` checks output with external nanopore tools and is run
  by hand. See [../../scripts/acceptance/README.md](../../scripts/acceptance/README.md).

### Writing tests

- Put `//go:build !release` at the top of every test file.
- Use `package X` for white-box tests of unexported code and `package X_test`
  for tests of the public API.
- Test behaviour that can catch a real bug: known values from the literature
  or a specification (F_ST in [0, 1], Hardy-Weinberg p-values, SLOW5 byte
  layouts), round trips, and the error paths. Coverage percentage is not a
  goal. `make coverage` and the pre-commit coverage hook use
  `scripts/check-coverage.sh 85 80`, but the current tree does not meet those
  thresholds (75.8% overall from
  `CGO_ENABLED=0 go test -tags '!release' -coverprofile=c.out ./...` and
  `go tool cover -func=c.out`, 2026-10-18) and CI does not enforce them.
- Generation must be reproducible from the seed. Draw randomness from the
  engine's seeded streams, never from `crypto/rand` or the global `math/rand`.

## Pre-commit hooks

`.pre-commit-config.yaml` configures hooks for whitespace, YAML and JSON
syntax, `gofmt`, `goimports`, `go vet`, `go mod tidy`, unit tests,
`golangci-lint`, `gosec`, `staticcheck`, `yamllint`, `detect-secrets` (with
`.secrets.baseline`) and the coverage script. They are optional and are not
run in CI.

```bash
make install-pre-commit
pre-commit run --all-files
```

## Changes and pull requests

- Branch from `main`. Use conventional commit messages (`feat:`, `fix:`,
  `docs:`, `test:`, `chore:`) with the issue number, for example
  `fix: ... (#21)`.
- Add a line under `## [Unreleased]` in `CHANGELOG.md` and run
  `make validate-changelog`.
- Run `go vet` and both test configurations before opening a pull request.
- Document new commands and flags in the command's `--help` text and in the
  README. Every documented command must run, and every number must come from
  a test or command in the repository.

## Repository layout

```
main.go                  CLI entry point
cmd/                     cobra commands (generate, validate, verify, train,
                         scenario, transcriptome, signal, bench, tui)
cmd/run_e2e_tests/       end-to-end harness and ML calibration
cmd/run_workflow_tests/  Cromwell / Nextflow / Snakemake checks
cmd/experiments/         three developer tools (benchmark_memory, benchmark_reference, sequence_validator)
internal/cli/            config loading and command plumbing
internal/fingerprinting/ stegano engine, keys, verifier
internal/genomics/       generation engine and read pipeline
  authenticity/          HMM models and statistics
  nanopore/              ONT headers, run metadata, summary, errmodel/
  readlen/               read-length distributions
  population/            structured coalescent and presets
  signal/ slow5/         raw-signal synthesis and the BLOW5 writer
  transcriptome/         cDNA sampler and isoform truth
  clinical/              clinical scenarios
internal/storage/s3sink/ S3 upload of finished files
internal/tui/ timing/    terminal UI and its run-timing simulator
internal/validation/     validation framework used by `validate`
pkg/types/               configuration types
configs/                 example sequencer and phantom configs
scripts/                 build, release and acceptance scripts
deploy/production/       container image
docs/                    this documentation
```
