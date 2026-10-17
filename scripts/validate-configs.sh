#!/bin/bash
# Configuration validation script for Phantom Science
# Validates YAML/JSON syntax and schema compliance

set -euo pipefail

echo "🧬 Validating configuration files..."

ERRORS=0

# Function to validate YAML/JSON syntax
validate_syntax() {
    local file="$1"
    local ext="${file##*.}"
    
    case "$ext" in
        "yaml"|"yml")
            if ! yamllint -d '{extends: default, rules: {line-length: {max: 120}}}' "$file" >/dev/null 2>&1; then
                echo "❌ YAML syntax error in $file"
                yamllint "$file" || true
                ((ERRORS++))
            else
                echo "✅ $file: Valid YAML syntax"
            fi
            ;;
        "json")
            if ! python3 -m json.tool "$file" >/dev/null 2>&1; then
                echo "❌ JSON syntax error in $file"
                ((ERRORS++))
            else
                echo "✅ $file: Valid JSON syntax"
            fi
            ;;
    esac
}

# Function to validate schema version
validate_schema() {
    local file="$1"
    local ext="${file##*.}"
    
    case "$ext" in
        "yaml"|"yml")
            if grep -q "schema_version:" "$file"; then
                VERSION=$(grep "schema_version:" "$file" | head -1 | cut -d'"' -f2 | cut -d"'" -f2 | xargs)
                if [[ "$VERSION" == "1.0" ]]; then
                    echo "✅ $file: Valid schema version $VERSION"
                else
                    echo "⚠️  $file: Unsupported schema version $VERSION"
                fi
            else
                echo "⚠️  $file: Missing schema_version (legacy format)"
            fi
            ;;
    esac
}

# Function to validate required fields
validate_required_fields() {
    local file="$1"
    local type="$2"
    
    case "$type" in
        "sequencer")
            # Check for required sequencer fields
            if [[ "$file" == *.yaml ]] || [[ "$file" == *.yml ]]; then
                REQUIRED_FIELDS=("sequencer:" "sample:" "reads:" "quality:")
                for field in "${REQUIRED_FIELDS[@]}"; do
                    if ! grep -q "$field" "$file"; then
                        echo "❌ $file: Missing required field $field"
                        ((ERRORS++))
                    fi
                done
            fi
            ;;
        "phantom")
            # Check for required phantom fields
            if [[ "$file" == *.yaml ]] || [[ "$file" == *.yml ]]; then
                REQUIRED_FIELDS=("enterprise:" "fingerprinting:" "validation:")
                for field in "${REQUIRED_FIELDS[@]}"; do
                    if ! grep -q "$field" "$file"; then
                        echo "❌ $file: Missing required field $field"
                        ((ERRORS++))
                    fi
                done
            fi
            ;;
    esac
}

# Find and validate all config files
echo "Searching for configuration files..."

# Validate sequencer configs
if [ -d "configs/sequencer" ]; then
    echo "Validating sequencer configurations..."
    find configs/sequencer -name "*.yaml" -o -name "*.yml" -o -name "*.json" | while read -r file; do
        validate_syntax "$file"
        validate_schema "$file"
        validate_required_fields "$file" "sequencer"
    done
fi

# Validate phantom configs  
if [ -d "configs/phantom" ]; then
    echo "Validating phantom configurations..."
    find configs/phantom -name "*.yaml" -o -name "*.yml" -o -name "*.json" | while read -r file; do
        validate_syntax "$file"
        validate_schema "$file"
        validate_required_fields "$file" "phantom"
    done
fi

# Validate any other config files in root
for file in *.yaml *.yml *.json 2>/dev/null; do
    if [ -f "$file" ]; then
        validate_syntax "$file"
        validate_schema "$file"
    fi
done

# Test config loading with Go
echo "Testing configuration loading..."
if command -v go >/dev/null 2>&1; then
    # Test loading each config with the actual Go application
    for config in configs/sequencer/v1.0/*.yaml configs/sequencer/*.json 2>/dev/null; do
        if [ -f "$config" ]; then
            echo "Testing Go config loading: $config"
            if ! timeout 10s go run main.go -config "$config" -version >/dev/null 2>&1; then
                echo "❌ Failed to load config with Go: $config"
                ((ERRORS++))
            else
                echo "✅ Go config loading successful: $config"
            fi
        fi
    done
fi

# Final result
if [ $ERRORS -eq 0 ]; then
    echo "✅ All configuration files validated successfully!"
    exit 0
else
    echo "❌ Configuration validation failed with $ERRORS errors"
    exit 1
fi