#!/bin/bash
# Sequencer specification validation script
# Validates parameters against known Illumina sequencer specifications

set -euo pipefail

echo "🧬 Validating sequencer specifications against vendor specs..."

ERRORS=0

# Function to validate NovaSeq 6000 parameters
validate_novaseq_6000() {
    local file="$1"
    echo "Validating NovaSeq 6000 parameters in $file..."
    
    # Extract parameters using yq or grep
    if command -v yq >/dev/null 2>&1; then
        READ_LENGTH=$(yq e '.reads.length' "$file" 2>/dev/null || echo "")
        FLOW_CELL=$(yq e '.sequencer.flow_cell_type' "$file" 2>/dev/null || echo "")
        Q30_PCT=$(yq e '.quality.q30_percentage' "$file" 2>/dev/null || echo "")
    else
        READ_LENGTH=$(grep -A1 "length:" "$file" | tail -1 | awk '{print $1}' || echo "")
        FLOW_CELL=$(grep -A1 "flow_cell_type:" "$file" | tail -1 | awk '{print $1}' | tr -d '"' || echo "")
        Q30_PCT=$(grep -A1 "q30_percentage:" "$file" | tail -1 | awk '{print $1}' || echo "")
    fi
    
    # Validate read lengths (50, 100, 150, 250bp supported)
    if [ -n "$READ_LENGTH" ]; then
        case "$READ_LENGTH" in
            50|100|150|250)
                echo "✅ Valid NovaSeq read length: ${READ_LENGTH}bp"
                ;;
            *)
                echo "❌ Invalid NovaSeq read length: ${READ_LENGTH}bp (supported: 50, 100, 150, 250)"
                ((ERRORS++))
                ;;
        esac
    fi
    
    # Validate flow cell types (S1, S2, S4, SP)
    if [ -n "$FLOW_CELL" ]; then
        case "$FLOW_CELL" in
            S1|S2|S4|SP)
                echo "✅ Valid NovaSeq flow cell: $FLOW_CELL"
                ;;
            *)
                echo "❌ Invalid NovaSeq flow cell: $FLOW_CELL (supported: S1, S2, S4, SP)"
                ((ERRORS++))
                ;;
        esac
    fi
    
    # Validate Q30 percentage (should be ≥90% for NovaSeq)
    if [ -n "$Q30_PCT" ]; then
        if (( $(echo "$Q30_PCT >= 90" | bc -l 2>/dev/null || echo "0") )); then
            echo "✅ Valid NovaSeq Q30 percentage: ${Q30_PCT}%"
        else
            echo "⚠️  Low NovaSeq Q30 percentage: ${Q30_PCT}% (typically ≥90%)"
        fi
    fi
}

# Function to validate NextSeq 2000 parameters
validate_nextseq_2000() {
    local file="$1"
    echo "Validating NextSeq 2000 parameters in $file..."
    
    # Extract parameters
    if command -v yq >/dev/null 2>&1; then
        READ_LENGTH=$(yq e '.reads.length' "$file" 2>/dev/null || echo "")
        FLOW_CELL=$(yq e '.sequencer.flow_cell_type' "$file" 2>/dev/null || echo "")
        CHEMISTRY=$(yq e '.sequencer.chemistry_version' "$file" 2>/dev/null || echo "")
    else
        READ_LENGTH=$(grep -A1 "length:" "$file" | tail -1 | awk '{print $1}' || echo "")
        FLOW_CELL=$(grep -A1 "flow_cell_type:" "$file" | tail -1 | awk '{print $1}' | tr -d '"' || echo "")
        CHEMISTRY=$(grep -A1 "chemistry_version:" "$file" | tail -1 | awk '{print $1}' | tr -d '"' || echo "")
    fi
    
    # Validate read lengths (50-300bp in various configurations)
    if [ -n "$READ_LENGTH" ]; then
        if [ "$READ_LENGTH" -ge 50 ] && [ "$READ_LENGTH" -le 300 ]; then
            echo "✅ Valid NextSeq read length: ${READ_LENGTH}bp"
        else
            echo "❌ Invalid NextSeq read length: ${READ_LENGTH}bp (supported: 50-300bp)"
            ((ERRORS++))
        fi
    fi
    
    # Validate flow cell types (P1, P2 for NextSeq 1000; P3, P4 for NextSeq 2000)
    if [ -n "$FLOW_CELL" ]; then
        case "$FLOW_CELL" in
            P1|P2|P3|P4)
                echo "✅ Valid NextSeq flow cell: $FLOW_CELL"
                ;;
            *)
                echo "❌ Invalid NextSeq flow cell: $FLOW_CELL (supported: P1, P2, P3, P4)"
                ((ERRORS++))
                ;;
        esac
    fi
    
    # Validate chemistry (should be XLEAP-SBS for NextSeq 2000)
    if [ -n "$CHEMISTRY" ] && [ "$CHEMISTRY" != "null" ]; then
        if [[ "$CHEMISTRY" == *"XLEAP"* ]]; then
            echo "✅ Valid NextSeq 2000 chemistry: $CHEMISTRY"
        else
            echo "⚠️  Unexpected NextSeq 2000 chemistry: $CHEMISTRY (expected: XLEAP-SBS)"
        fi
    fi
}

# Function to validate MiSeq parameters
validate_miseq() {
    local file="$1"
    echo "Validating MiSeq parameters in $file..."
    
    # Extract parameters
    if command -v yq >/dev/null 2>&1; then
        READ_LENGTH=$(yq e '.reads.length' "$file" 2>/dev/null || echo "")
        FLOW_CELL=$(yq e '.sequencer.flow_cell_type' "$file" 2>/dev/null || echo "")
        MAX_CYCLES=$(yq e '.performance.max_cycles' "$file" 2>/dev/null || echo "")
    else
        READ_LENGTH=$(grep -A1 "length:" "$file" | tail -1 | awk '{print $1}' || echo "")
        FLOW_CELL=$(grep -A1 "flow_cell_type:" "$file" | tail -1 | awk '{print $1}' | tr -d '"' || echo "")
        MAX_CYCLES=$(grep -A1 "max_cycles:" "$file" | tail -1 | awk '{print $1}' || echo "")
    fi
    
    # Validate read lengths (up to 300bp for v3)
    if [ -n "$READ_LENGTH" ]; then
        case "$FLOW_CELL" in
            v3)
                if [ "$READ_LENGTH" -le 300 ]; then
                    echo "✅ Valid MiSeq v3 read length: ${READ_LENGTH}bp"
                else
                    echo "❌ Invalid MiSeq v3 read length: ${READ_LENGTH}bp (max: 300bp)"
                    ((ERRORS++))
                fi
                ;;
            v2)
                if [ "$READ_LENGTH" -le 250 ]; then
                    echo "✅ Valid MiSeq v2 read length: ${READ_LENGTH}bp"
                else
                    echo "❌ Invalid MiSeq v2 read length: ${READ_LENGTH}bp (max: 250bp)"
                    ((ERRORS++))
                fi
                ;;
            *)
                echo "⚠️  Unknown MiSeq flow cell type: $FLOW_CELL"
                ;;
        esac
    fi
    
    # Validate max cycles for v3 (up to 625)
    if [ -n "$MAX_CYCLES" ] && [ "$FLOW_CELL" = "v3" ]; then
        if [ "$MAX_CYCLES" -le 625 ]; then
            echo "✅ Valid MiSeq v3 max cycles: $MAX_CYCLES"
        else
            echo "❌ Invalid MiSeq v3 max cycles: $MAX_CYCLES (max: 625)"
            ((ERRORS++))
        fi
    fi
}

# Main validation function
validate_sequencer_file() {
    local file="$1"
    
    if [ ! -f "$file" ]; then
        echo "❌ File not found: $file"
        ((ERRORS++))
        return
    fi
    
    echo "Validating sequencer file: $file"
    
    # Determine sequencer type from filename or content
    if [[ "$file" == *"novaseq"* ]]; then
        validate_novaseq_6000 "$file"
    elif [[ "$file" == *"nextseq"* ]]; then
        validate_nextseq_2000 "$file"  
    elif [[ "$file" == *"miseq"* ]]; then
        validate_miseq "$file"
    else
        # Try to determine from content
        if grep -q "NovaSeq" "$file"; then
            validate_novaseq_6000 "$file"
        elif grep -q "NextSeq" "$file"; then
            validate_nextseq_2000 "$file"
        elif grep -q "MiSeq" "$file"; then
            validate_miseq "$file"
        else
            echo "⚠️  Unknown sequencer type in $file"
        fi
    fi
}

# Install yq if available and not installed (improves parsing)
if ! command -v yq >/dev/null 2>&1; then
    echo "Note: yq not available, using basic grep parsing"
    echo "Install yq for better YAML parsing: https://github.com/mikefarah/yq"
fi

# Install bc for arithmetic if not available
if ! command -v bc >/dev/null 2>&1; then
    echo "Warning: bc not available, skipping numeric validations"
fi

# Process all sequencer YAML files
find configs/sequencer -name "*.yaml" -o -name "*.yml" 2>/dev/null | while read -r file; do
    validate_sequencer_file "$file"
done

# Final result
if [ $ERRORS -eq 0 ]; then
    echo "✅ All sequencer specifications validated successfully!"
    exit 0
else
    echo "❌ Sequencer specification validation failed with $ERRORS errors"
    exit 1
fi