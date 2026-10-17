#!/bin/bash
# Coverage checking script for Phantom Science
# Usage: check-coverage.sh <overall_threshold> <individual_threshold>

set -euo pipefail

# Default thresholds
OVERALL_THRESHOLD=${1:-85}
INDIVIDUAL_THRESHOLD=${2:-80}

echo "🧬 Checking Go test coverage..."
echo "Required: ${OVERALL_THRESHOLD}% overall, ${INDIVIDUAL_THRESHOLD}% individual files"

# Create temporary coverage file
COVERAGE_FILE=$(mktemp)
trap "rm -f $COVERAGE_FILE" EXIT

# Run tests with coverage
echo "Running tests with coverage..."
go test ./... -coverprofile="$COVERAGE_FILE" -covermode=atomic

if [ ! -s "$COVERAGE_FILE" ]; then
    echo "❌ No coverage data generated"
    exit 1
fi

# Calculate overall coverage
OVERALL_COVERAGE=$(go tool cover -func="$COVERAGE_FILE" | grep "total:" | awk '{print $3}' | sed 's/%//')

echo "Overall coverage: ${OVERALL_COVERAGE}%"

# Check overall threshold
if (( $(echo "$OVERALL_COVERAGE < $OVERALL_THRESHOLD" | bc -l) )); then
    echo "❌ Overall coverage ${OVERALL_COVERAGE}% is below required ${OVERALL_THRESHOLD}%"
    exit 1
fi

# Check individual file coverage
echo "Checking individual file coverage..."
FAILED_FILES=()

while IFS= read -r line; do
    if [[ $line == *".go:"* ]] && [[ $line != *"total:"* ]]; then
        FILE=$(echo "$line" | awk '{print $1}')
        COVERAGE=$(echo "$line" | awk '{print $3}' | sed 's/%//')
        
        if (( $(echo "$COVERAGE < $INDIVIDUAL_THRESHOLD" | bc -l) )); then
            FAILED_FILES+=("$FILE: ${COVERAGE}%")
        fi
    fi
done < <(go tool cover -func="$COVERAGE_FILE")

# Report results
if [ ${#FAILED_FILES[@]} -gt 0 ]; then
    echo "❌ Files below ${INDIVIDUAL_THRESHOLD}% coverage:"
    printf '  %s\n' "${FAILED_FILES[@]}"
    exit 1
fi

echo "✅ All coverage requirements met!"
echo "  Overall: ${OVERALL_COVERAGE}% (required: ${OVERALL_THRESHOLD}%)"
echo "  Individual files: All above ${INDIVIDUAL_THRESHOLD}%"

# Generate HTML report for local development
if [ "${CI:-false}" != "true" ]; then
    echo "Generating HTML coverage report..."
    go tool cover -html="$COVERAGE_FILE" -o coverage.html
    echo "📊 Coverage report saved to coverage.html"
fi