#!/bin/bash
# Phantom Science release script.
#
# Prepares a release in the working tree: moves the [Unreleased] changelog
# entries under "## [X.Y.Z] - YYYY-MM-DD", leaves an empty [Unreleased]
# section, updates the compare links at the bottom of CHANGELOG.md, sets
# VERSION and cmd/root.go Version/BuildDate, and runs `make
# validate-changelog` and `make check-version`.
#
# Nothing is committed, tagged or pushed unless asked (--commit, --tag,
# --push). See --help.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${BLUE}[INFO]${NC} $*"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $*"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VERSION_FILE="$PROJECT_ROOT/VERSION"
CHANGELOG="$PROJECT_ROOT/CHANGELOG.md"
ROOT_GO="$PROJECT_ROOT/cmd/root.go"

usage() {
    cat <<EOF
Usage: $0 <X.Y.Z | patch | minor | major> [options]

Prepares release X.Y.Z (or the next patch/minor/major after VERSION):
  - CHANGELOG.md: the [Unreleased] entries become "## [X.Y.Z] - DATE", an
    empty "## [Unreleased]" is left above it, and the [Unreleased] compare
    link (if present) is updated with a new [X.Y.Z] link
  - VERSION and cmd/root.go Version/BuildDate are set
  - make validate-changelog and make check-version are run
By default the changes are left uncommitted in the working tree.

Options:
  --dry-run       Print the diff of every change without writing anything
  --date DATE     Release date (YYYY-MM-DD; default: today)
  --skip-tests    Do not run go test before preparing the release
  --build         Build a native binary into build/release/. The binary links
                  the stegano library through CGO, so it is built for this
                  machine only; there is no cross-compilation.
  --commit        Commit the release changes
  --tag           Create annotated tag vX.Y.Z (implies --commit)
  --push          Push the current branch and the tag (implies --tag)
  --help          Show this help

Examples:
  $0 0.8.0 --dry-run          # preview
  $0 0.8.0                    # prepare; review with git diff, then commit
  $0 minor --tag              # prepare, commit and tag; push by hand
EOF
}

TARGET=""
DRY_RUN=false
SKIP_TESTS=false
BUILD=false
COMMIT=false
TAG=false
PUSH=false
RELEASE_DATE="$(date +%Y-%m-%d)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        patch | minor | major) TARGET="$1" ;;
        [0-9]*.[0-9]*.[0-9]*) TARGET="$1" ;;
        --dry-run) DRY_RUN=true ;;
        --skip-tests) SKIP_TESTS=true ;;
        --build) BUILD=true ;;
        --commit) COMMIT=true ;;
        --tag) TAG=true COMMIT=true ;;
        --push) PUSH=true TAG=true COMMIT=true ;;
        --date)
            [[ $# -ge 2 ]] || { log_error "--date needs a value"; exit 1; }
            RELEASE_DATE="$2"
            shift
            ;;
        --help | -h)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown argument: $1"
            usage >&2
            exit 1
            ;;
    esac
    shift
done

if [[ -z "$TARGET" ]]; then
    log_error "A version (X.Y.Z) or patch/minor/major is required"
    usage >&2
    exit 1
fi
if ! [[ "$RELEASE_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
    log_error "Invalid --date $RELEASE_DATE (expected YYYY-MM-DD)"
    exit 1
fi

cd "$PROJECT_ROOT"

current_version() {
    if [[ -f "$VERSION_FILE" ]]; then
        sed 's/-.*$//' "$VERSION_FILE" | tr -d '[:space:]'
    else
        echo "0.0.0"
    fi
}

next_version() {
    local current="$1" kind="$2" major minor patch
    IFS='.' read -r major minor patch <<<"$current"
    case "$kind" in
        patch) echo "$major.$minor.$((patch + 1))" ;;
        minor) echo "$major.$((minor + 1)).0" ;;
        major) echo "$((major + 1)).0.0" ;;
        *) echo "$kind" ;;
    esac
}

CURRENT="$(current_version)"
NEW="$(next_version "$CURRENT" "$TARGET")"
if ! [[ "$NEW" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    log_error "Invalid version: $NEW (expected X.Y.Z)"
    exit 1
fi
if grep -q "^## \[$NEW\]" "$CHANGELOG"; then
    log_error "CHANGELOG.md already has a [$NEW] section"
    exit 1
fi
if ! grep -q '^## \[Unreleased\]' "$CHANGELOG"; then
    log_error "CHANGELOG.md has no '## [Unreleased]' section"
    exit 1
fi
# The [Unreleased] section must have entries to release.
if ! awk '/^## \[Unreleased\]/ {in_u=1; next} /^## \[/ {in_u=0} in_u && /^- / {found=1} END {exit !found}' "$CHANGELOG"; then
    log_error "The [Unreleased] section of CHANGELOG.md has no entries to release"
    exit 1
fi

log_info "Current version: $CURRENT"
log_info "New version:     $NEW ($RELEASE_DATE)"

# rewrite_changelog IN OUT writes the released changelog.
rewrite_changelog() {
    awk -v ver="$NEW" -v date="$RELEASE_DATE" '
        # Compare-link lines: [Unreleased]: <base>/compare/vOLD...HEAD
        /^\[Unreleased\]: / && !done_link {
            line = $0
            if (match(line, /\/compare\/[^\/]+\.\.\.HEAD$/)) {
                base = substr(line, length("[Unreleased]: ") + 1, RSTART - length("[Unreleased]: ") - 1)
                prev = substr(line, RSTART + length("/compare/"), RLENGTH - length("/compare/") - length("...HEAD"))
                print "[Unreleased]: " base "/compare/v" ver "...HEAD"
                print "[" ver "]: " base "/compare/" prev "...v" ver
                done_link = 1
                next
            }
        }
        /^## \[Unreleased\]/ && !done_head {
            print "## [Unreleased]"
            print ""
            print "## [" ver "] - " date
            done_head = 1
            next
        }
        { print }
    ' "$1" >"$2"
}

# rewrite_version_go IN OUT sets Version and BuildDate in cmd/root.go.
rewrite_version_go() {
    sed -E \
        -e "s/^([[:space:]]*Version[[:space:]]*=[[:space:]]*)\"[^\"]*\"/\\1\"$NEW\"/" \
        -e "s/^([[:space:]]*BuildDate[[:space:]]*=[[:space:]]*)\"[^\"]*\"/\\1\"$RELEASE_DATE\"/" \
        "$1" >"$2"
}

TMPDIR_RELEASE="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_RELEASE"' EXIT

rewrite_changelog "$CHANGELOG" "$TMPDIR_RELEASE/CHANGELOG.md"
rewrite_version_go "$ROOT_GO" "$TMPDIR_RELEASE/root.go"
echo "$NEW" >"$TMPDIR_RELEASE/VERSION"

grep -q "Version[[:space:]]*=[[:space:]]*\"$NEW\"" "$TMPDIR_RELEASE/root.go" ||
    { log_error "Could not set Version in cmd/root.go"; exit 1; }
grep -q "BuildDate[[:space:]]*=[[:space:]]*\"$RELEASE_DATE\"" "$TMPDIR_RELEASE/root.go" ||
    { log_error "Could not set BuildDate in cmd/root.go"; exit 1; }

if [[ "$DRY_RUN" == true ]]; then
    log_warning "DRY RUN: nothing is written"
    for pair in "CHANGELOG.md:$TMPDIR_RELEASE/CHANGELOG.md" "VERSION:$TMPDIR_RELEASE/VERSION" "cmd/root.go:$TMPDIR_RELEASE/root.go"; do
        rel="${pair%%:*}"
        diff -u --label "a/$rel" --label "b/$rel" "$rel" "${pair#*:}" || true
    done
    log_info "Validating the new changelog..."
    "$SCRIPT_DIR/validate-changelog.sh" "$TMPDIR_RELEASE/CHANGELOG.md" >/dev/null ||
        { log_error "The new CHANGELOG.md would fail validate-changelog"; exit 1; }
    log_success "The new CHANGELOG.md passes validate-changelog"
    [[ "$SKIP_TESTS" == true ]] || log_info "Would run: go test -tags '!release' ./..."
    [[ "$BUILD" == true ]] && log_info "Would build a native binary into build/release/"
    [[ "$COMMIT" == true ]] && log_info "Would commit the release changes"
    [[ "$TAG" == true ]] && log_info "Would create tag v$NEW"
    [[ "$PUSH" == true ]] && log_info "Would push the branch and tag v$NEW"
    exit 0
fi

if ! git diff-index --quiet HEAD --; then
    log_error "Working tree has uncommitted changes; commit or stash them first"
    exit 1
fi
if [[ "$TAG" == true ]] && git rev-parse -q --verify "refs/tags/v$NEW" >/dev/null; then
    log_error "Tag v$NEW already exists"
    exit 1
fi

if [[ "$SKIP_TESTS" != true ]]; then
    log_info "Running tests..."
    stegano_lib="$PROJECT_ROOT/.stegano/target/release"
    if [[ -z "${CGO_LDFLAGS:-}" && -d "$stegano_lib" ]]; then
        export CGO_LDFLAGS="-L$stegano_lib -Wl,-rpath,$stegano_lib"
    fi
    if ! go test -tags '!release' ./... -timeout 900s; then
        log_error "Tests failed"
        exit 1
    fi
    log_success "Tests passed"
fi

cp "$TMPDIR_RELEASE/CHANGELOG.md" "$CHANGELOG"
cp "$TMPDIR_RELEASE/VERSION" "$VERSION_FILE"
cp "$TMPDIR_RELEASE/root.go" "$ROOT_GO"
log_success "Updated CHANGELOG.md, VERSION and cmd/root.go"

restore() {
    git checkout -- CHANGELOG.md VERSION cmd/root.go
}
if ! make --no-print-directory validate-changelog; then
    log_error "validate-changelog failed; restoring CHANGELOG.md, VERSION and cmd/root.go"
    restore
    exit 1
fi
check_out="$(make --no-print-directory check-version 2>&1)" || {
    echo "$check_out"
    log_error "check-version failed; restoring CHANGELOG.md, VERSION and cmd/root.go"
    restore
    exit 1
}
echo "$check_out"
if grep -q "doesn't match" <<<"$check_out"; then
    log_error "VERSION and cmd/root.go disagree; restoring CHANGELOG.md, VERSION and cmd/root.go"
    restore
    exit 1
fi

if [[ "$BUILD" == true ]]; then
    goos="$(go env GOOS)"
    goarch="$(go env GOARCH)"
    out="build/release/phantom-science-$goos-$goarch"
    log_info "Building $out (native only: CGO and the stegano library)..."
    mkdir -p build/release
    go build -o "$out" .
    log_success "Built $out"
fi

if [[ "$COMMIT" == true ]]; then
    git add CHANGELOG.md VERSION cmd/root.go
    git commit -m "chore: release v$NEW"
    log_success "Committed release v$NEW"
fi
if [[ "$TAG" == true ]]; then
    git tag -a "v$NEW" -m "Phantom Science v$NEW

See CHANGELOG.md for the changes in this release."
    log_success "Created tag v$NEW"
fi
if [[ "$PUSH" == true ]]; then
    branch="$(git rev-parse --abbrev-ref HEAD)"
    git push origin "$branch"
    git push origin "v$NEW"
    log_success "Pushed $branch and v$NEW"
fi

echo
log_success "Release v$NEW prepared"
if [[ "$COMMIT" != true ]]; then
    log_info "Review with: git diff -- CHANGELOG.md VERSION cmd/root.go"
    log_info "Then commit and tag by hand: git commit -am \"chore: release v$NEW\" && git tag -a v$NEW"
elif [[ "$PUSH" != true ]]; then
    log_info "Push with: git push origin $(git rev-parse --abbrev-ref HEAD) v$NEW"
fi
