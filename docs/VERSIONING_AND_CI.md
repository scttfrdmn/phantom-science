# Versioning, Changelog and CI

## Versioning

The project follows [Semantic Versioning 2.0.0](https://semver.org/spec/v2.0.0.html).
It is pre-1.0, so a minor release can contain breaking changes; the
changelog marks them **BREAKING**.

The version is kept in two places, which must agree:

- `VERSION`, the source of truth
- `Version` in `cmd/root.go`, printed by `phantom-science --version`

```bash
cat VERSION
make check-version   # VERSION is valid SemVer and matches cmd/root.go
make sync-version    # copy VERSION into cmd/root.go
```

`make version` prints the version on case-sensitive file systems. On macOS's
default case-insensitive file system, make treats the `VERSION` file as the
`version` target and prints "is up to date" instead.

## Changelog

`CHANGELOG.md` follows [Keep a Changelog 1.0.0](https://keepachangelog.com/en/1.0.0/).
Changes go under `## [Unreleased]` in Added, Changed, Deprecated, Removed,
Fixed or Security, with the issue number. A release renames that section to
`## [X.Y.Z] - YYYY-MM-DD`.

```bash
make validate-changelog   # runs scripts/validate-changelog.sh
```

The script checks for:

- the header and its references to Keep a Changelog and Semantic Versioning
- valid version headers
- standard change-type sections
- an `[Unreleased]` section
- YYYY-MM-DD dates
- empty sections
- bullet formatting

Formatting problems are reported as warnings.

## Releases

Releases are cut by hand:

1. Set the new version in `VERSION` and run `make sync-version`.
2. In `CHANGELOG.md`, rename `[Unreleased]` to the version and date, then add
   a new empty `[Unreleased]` section.
3. Run `make check-version validate-changelog` and the test suite.
4. Commit, then tag `vX.Y.Z`.

`scripts/release.sh` automates part of this. Do not rely on it yet:

- Its changelog rewrite expects an `[Unreleased] - vX.Y.Z-dev` heading, which
  the changelog does not use.
- It does not update `cmd/root.go`.
- Its cross-compiled binaries are built without the stegano library, so they
  cannot write or read quality-score fingerprints.

## Continuous integration

`.github/workflows/ci.yml` runs on pushes and pull requests to `main`, and on
manual dispatch. It has no schedule. All jobs run on `ubuntu-latest` with the
Go toolchain pinned in `go.mod`.

| Job | Command |
|---|---|
| go vet | `CGO_ENABLED=0 go vet -tags '!release' ./...` |
| govulncheck | `govulncheck ./...` |
| Test (CGO_ENABLED=0) | `CGO_ENABLED=0 go test -tags '!release' ./... -timeout 15m` |
| Test (CGO, stegano) | `scripts/build-stegano.sh`, then `go test -tags '!release' ./... -timeout 20m` with CGO linked against `.stegano/` |
| Release build | `CGO_ENABLED=0 go build -tags release ./...` |
| Docker build | builds `deploy/production/docker/Dockerfile` and runs `--version` in it |

CI does not run coverage thresholds, gosec, staticcheck, golangci-lint,
yamllint, the config validators or the acceptance scripts in
`scripts/acceptance/`.

To run the CI test jobs locally:

```bash
CGO_ENABLED=0 go vet -tags '!release' ./...
CGO_ENABLED=0 go test -tags '!release' ./... -timeout 15m

make stegano
L="$PWD/.stegano/target/release"
CGO_LDFLAGS="-L$L -Wl,-rpath,$L" go test -tags '!release' ./... -timeout 20m
```

`make ci` is a separate local target, not what GitHub runs. It runs `test`,
`coverage` (85% overall and 80% per file), `lint`, `security`,
`validate-configs` and `validate-changelog`, and does not pass on the current
tree: for example, `make validate-configs` fails on yamllint errors in
`configs/`.

## CI status

`scripts/ci-status-monitor.sh` queries the GitHub API (no token needed for the
public repository) and summarises the latest workflow runs. It needs `curl`
and `jq`.

```bash
make ci-status                                   # once
make ci-watch                                    # every 5 minutes
scripts/ci-status-monitor.sh --json              # JSON summary of the last 10 runs
scripts/ci-status-monitor.sh --watch --webhook https://hooks.example/...   # POST an alert on failures
```

Or with the GitHub CLI: `gh run list -R scttfrdmn/phantom-science`.
