#!/bin/bash
# Configuration validation script for Phantom Science
# Validates YAML/JSON syntax (yamllint with the repository's .yamllint),
# schema versions and required fields, then loads every shipped sequencer
# config through the Go config loader.

set -euo pipefail

cd "$(dirname "$0")/.."

echo "🧬 Validating configuration files..."

ERRORS=0
fail() {
    echo "❌ $*"
    ERRORS=$((ERRORS + 1))
}

if ! command -v yamllint >/dev/null 2>&1; then
    echo "❌ yamllint not found, install with: make install-tools"
    exit 1
fi

# Validate YAML/JSON syntax
validate_syntax() {
    local file="$1"
    case "${file##*.}" in
        yaml|yml)
            if ! yamllint -c .yamllint "$file" >/dev/null 2>&1; then
                fail "YAML lint errors in $file"
                yamllint -c .yamllint "$file" || true
            else
                echo "✅ $file: Valid YAML syntax"
            fi
            ;;
        json)
            if ! python3 -m json.tool "$file" >/dev/null 2>&1; then
                fail "JSON syntax error in $file"
            else
                echo "✅ $file: Valid JSON syntax"
            fi
            ;;
    esac
}

# Report the schema version of v1 YAML configs
validate_schema() {
    local file="$1" version
    case "${file##*.}" in
        yaml|yml)
            if grep -q "schema_version:" "$file"; then
                version=$(grep "schema_version:" "$file" | head -1 | cut -d'"' -f2 | cut -d"'" -f2 | xargs)
                if [[ "$version" == "1.0" ]]; then
                    echo "✅ $file: Valid schema version $version"
                else
                    echo "⚠️  $file: Unsupported schema version $version"
                fi
            else
                echo "⚠️  $file: Missing schema_version (legacy format)"
            fi
            ;;
    esac
}

# Check the top-level sections of v1 YAML configs
validate_required_fields() {
    local file="$1" type="$2" field
    local -a required
    case "$type" in
        sequencer) required=("sequencer:" "sample:" "reads:" "quality:") ;;
        phantom) required=("enterprise:" "fingerprinting:" "validation:") ;;
        *) return 0 ;;
    esac
    case "$file" in
        *.yaml|*.yml) ;;
        *) return 0 ;;
    esac
    for field in "${required[@]}"; do
        if ! grep -q "$field" "$file"; then
            fail "$file: Missing required field $field"
        fi
    done
}

# validate_dir DIR TYPE validates every config file under DIR.
validate_dir() {
    local dir="$1" type="$2" file
    [ -d "$dir" ] || return 0
    echo "Validating $dir..."
    # Process substitution (not a pipe) so ERRORS updates are kept.
    while read -r file; do
        validate_syntax "$file"
        if [ "$type" != "template" ]; then
            validate_schema "$file"
            validate_required_fields "$file" "$type"
        fi
    done < <(find "$dir" \( -name "*.yaml" -o -name "*.yml" -o -name "*.json" \) | sort)
}

validate_dir configs/sequencer sequencer
validate_dir configs/phantom phantom
validate_dir configs/sequencers template
validate_dir configs/dev template

# Load every shipped sequencer config with the Go config loader.
echo "Testing configuration loading..."
if command -v go >/dev/null 2>&1; then
    if CGO_ENABLED=0 go test -tags '!release' -count=1 -run 'TestLoadSequencerConfig_AllShippedConfigs' ./internal/cli/ >/dev/null 2>&1; then
        echo "✅ Go config loading successful for configs/sequencer/v1.0 and configs/sequencers"
    else
        fail "Go config loading failed; run: CGO_ENABLED=0 go test -tags '!release' -run TestLoadSequencerConfig_AllShippedConfigs -v ./internal/cli/"
    fi
else
    echo "⚠️  go not found; skipping config loading"
fi

# Final result
if [ "$ERRORS" -eq 0 ]; then
    echo "✅ All configuration files validated successfully!"
    exit 0
else
    echo "❌ Configuration validation failed with $ERRORS errors"
    exit 1
fi
