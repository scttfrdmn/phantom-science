#!/bin/bash
# Build the stegano shared library that internal/fingerprinting links against
# in CGO builds. The version comes from go.mod, so the library always matches
# the github.com/scttfrdmn/stegano/go module the Go code is compiled against.
#
# Usage: scripts/build-stegano.sh   (or: make stegano)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODULE="github.com/scttfrdmn/stegano/go"
VERSION="${STEGANO_VERSION:-$(cd "$ROOT" && go list -m -f '{{.Version}}' "$MODULE")}"
WORK="$ROOT/.stegano"
SRC="$WORK/src-$VERSION"
LIB_DIR="$WORK/target/release"

if ! command -v cargo >/dev/null 2>&1; then
    echo "❌ Rust is not installed. Install it from https://rustup.rs/" >&2
    exit 1
fi

if [ ! -d "$SRC" ]; then
    echo "📥 Fetching stegano $VERSION..."
    git clone --quiet --depth 1 --branch "$VERSION" https://github.com/scttfrdmn/stegano.git "$SRC" 2>/dev/null
fi

# Build into a private target directory: with a shared CARGO_TARGET_DIR, a
# build of another stegano checkout would overwrite this library.
echo "📦 Building stegano $VERSION..."
CARGO_TARGET_DIR="$WORK/target" cargo build --release --quiet --manifest-path "$SRC/Cargo.toml"

if [ ! -e "$LIB_DIR/libstegano.dylib" ] && [ ! -e "$LIB_DIR/libstegano.so" ]; then
    echo "❌ Build finished but no libstegano shared library in $LIB_DIR" >&2
    exit 1
fi

echo "✅ stegano $VERSION built in $LIB_DIR"
echo ""
echo "make targets set CGO_LDFLAGS automatically. Outside make, use:"
echo "  export CGO_LDFLAGS=\"-L$LIB_DIR -Wl,-rpath,$LIB_DIR\""
