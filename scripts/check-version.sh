#!/bin/bash
# Semantic versioning validation script

set -euo pipefail

VERSION_FILE="VERSION"
GO_VERSION_FILE="cmd/root.go"

if [ ! -f "$VERSION_FILE" ]; then
    echo "❌ VERSION file not found"
    exit 1
fi

VERSION=$(cat "$VERSION_FILE")
echo "Checking version: $VERSION"

# Semantic versioning regex
SEMVER_REGEX='^[0-9]+\.[0-9]+\.[0-9]+(-[a-zA-Z0-9.-]+)?(\+[a-zA-Z0-9.-]+)?$'

if [[ $VERSION =~ $SEMVER_REGEX ]]; then
    echo "✅ Version follows semantic versioning"
    
    # Extract version components
    MAJOR=$(echo "$VERSION" | cut -d'.' -f1)
    MINOR=$(echo "$VERSION" | cut -d'.' -f2)
    PATCH=$(echo "$VERSION" | cut -d'.' -f3 | cut -d'-' -f1)
    
    echo "  Major: $MAJOR"
    echo "  Minor: $MINOR" 
    echo "  Patch: $PATCH"
    
    if [[ $VERSION =~ -.*$ ]]; then
        PRERELEASE=$(echo "$VERSION" | sed 's/^[0-9]*\.[0-9]*\.[0-9]*-//' | cut -d'+' -f1)
        echo "  Pre-release: $PRERELEASE"
    fi
    
    if [[ $VERSION =~ \+.*$ ]]; then
        BUILD=$(echo "$VERSION" | sed 's/^.*+//')
        echo "  Build metadata: $BUILD"
    fi
    
    # Check if Go version matches
    if [ -f "$GO_VERSION_FILE" ]; then
        GO_VERSION=$(grep 'Version.*=' "$GO_VERSION_FILE" | sed 's/.*"\(.*\)".*/\1/')
        if [ "$VERSION" != "$GO_VERSION" ]; then
            echo "⚠️  Warning: VERSION file ($VERSION) doesn't match Go version ($GO_VERSION)"
            echo "   Run: make sync-version to fix"
        else
            echo "✅ Go version matches VERSION file"
        fi
    fi
    
else
    echo "❌ Version does not follow semantic versioning"
    echo "Expected format: MAJOR.MINOR.PATCH[-PRERELEASE][+BUILD]"
    echo "Examples:"
    echo "  1.0.0"
    echo "  1.0.0-alpha"
    echo "  1.0.0-alpha.1"
    echo "  1.0.0+build.1"
    echo "  1.0.0-beta+exp.sha.5114f85"
    exit 1
fi