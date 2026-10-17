#!/bin/bash

# Phantom Science Cromwell Workflow Execution Script
# This script runs the phantom analysis workflow using Cromwell

set -euo pipefail

# Configuration
WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$WORKFLOW_DIR/../.." && pwd)"
CROMWELL_SCRIPT="$PROJECT_ROOT/cromwell"
WDL_FILE="$WORKFLOW_DIR/phantom_analysis.wdl"
INPUTS_FILE="$WORKFLOW_DIR/inputs.json"
OPTIONS_FILE="$WORKFLOW_DIR/options.json"
OUTPUT_DIR="$PROJECT_ROOT/cromwell_output"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check dependencies
check_dependencies() {
    log "Checking dependencies..."
    
    # Check if cromwell script exists
    if [[ ! -f "$CROMWELL_SCRIPT" ]]; then
        error "Cromwell script not found at: $CROMWELL_SCRIPT"
        exit 1
    fi
    
    # Check if cromwell script is executable
    if [[ ! -x "$CROMWELL_SCRIPT" ]]; then
        warn "Making cromwell script executable..."
        chmod +x "$CROMWELL_SCRIPT"
    fi
    
    # Check if WDL file exists
    if [[ ! -f "$WDL_FILE" ]]; then
        error "WDL file not found at: $WDL_FILE"
        exit 1
    fi
    
    # Check if inputs file exists
    if [[ ! -f "$INPUTS_FILE" ]]; then
        error "Inputs file not found at: $INPUTS_FILE"
        exit 1
    fi
    
    # Check if Java is available
    if ! command -v java &> /dev/null; then
        error "Java not found. Please install Java 11 or higher."
        exit 1
    fi
    
    # Check Java version
    java_version=$(java -version 2>&1 | head -1 | cut -d'"' -f2 | cut -d'.' -f1)
    if [[ "$java_version" -lt 11 ]]; then
        error "Java 11 or higher is required. Current version: $java_version"
        exit 1
    fi
    
    success "All dependencies satisfied"
}

# Validate inputs
validate_inputs() {
    log "Validating inputs..."
    
    # Parse inputs JSON to check if input files exist
    if command -v jq &> /dev/null; then
        fastq_files=$(jq -r '.["PhantomAnalysis.fastq_files"][]' "$INPUTS_FILE")
        
        for fastq in $fastq_files; do
            # Check if file exists (handle relative paths)
            if [[ ! -f "$fastq" && ! -f "$PROJECT_ROOT/$fastq" ]]; then
                warn "FASTQ file not found: $fastq"
                log "Creating sample FASTQ file for demonstration..."
                
                # Create sample FASTQ file
                sample_dir="$PROJECT_ROOT/$(dirname "$fastq")"
                mkdir -p "$sample_dir"
                
                cat > "$PROJECT_ROOT/$fastq" << 'EOF'
@Sample_001:1:HXXXXXX:1:1101:1000:2000 1:N:0:1
TGTAATCTGATACGGTGGAGTATAATTCGATTGGGGAAAAGAAGCTGACTACTCAGAAGAGTGGCAAGGAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@Sample_001:1:HXXXXXX:1:1101:1000:2001 1:N:0:1
CGTGAATCTGATACGGTGGAGTATAATTCGATTGGGGAAAAGAAGCTGACTACTCAGAAGAGTGGCAAGGA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@Sample_001:1:HXXXXXX:1:1101:1000:2002 1:N:0:1
AGTGAATCTGATACGGTGGAGTATAATTCGATTGGGGAAAAGAAGCTGACTACTCAGAAGAGTGGCAAGGA
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
EOF
                success "Created sample FASTQ file: $fastq"
            fi
        done
    else
        warn "jq not found. Skipping input file validation."
    fi
    
    success "Input validation completed"
}

# Run Cromwell workflow
run_workflow() {
    log "Starting Cromwell workflow execution..."
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Build Cromwell command
    cromwell_cmd=(
        "$CROMWELL_SCRIPT"
        run
        "$WDL_FILE"
        -i "$INPUTS_FILE"
    )
    
    # Add options file if it exists
    if [[ -f "$OPTIONS_FILE" ]]; then
        cromwell_cmd+=(-o "$OPTIONS_FILE")
    fi
    
    # Add metadata output
    cromwell_cmd+=(-m "$OUTPUT_DIR/metadata.json")
    
    log "Running command: ${cromwell_cmd[*]}"
    
    # Execute Cromwell workflow
    if "${cromwell_cmd[@]}"; then
        success "Workflow completed successfully!"
    else
        error "Workflow execution failed!"
        exit 1
    fi
}

# Generate summary report
generate_summary() {
    log "Generating execution summary..."
    
    summary_file="$OUTPUT_DIR/execution_summary.txt"
    
    cat > "$summary_file" << EOF
# Phantom Science Cromwell Workflow Execution Summary

## Execution Details
- Workflow: phantom_analysis.wdl
- Started: $(date)
- Output Directory: $OUTPUT_DIR
- Options: $OPTIONS_FILE

## Files Generated
EOF
    
    # List generated files
    if [[ -d "$OUTPUT_DIR" ]]; then
        find "$OUTPUT_DIR" -type f -name "*.txt" -o -name "*.html" -o -name "*.json" | sort >> "$summary_file"
    fi
    
    success "Summary report generated: $summary_file"
}

# Main execution
main() {
    log "🧬 Phantom Science Cromwell Workflow Execution"
    log "=============================================="
    
    check_dependencies
    validate_inputs
    run_workflow
    generate_summary
    
    log "=============================================="
    success "Workflow execution completed successfully!"
    log "Check the output directory: $OUTPUT_DIR"
    
    # Show key output files
    if [[ -f "$OUTPUT_DIR/reports/phantom_analysis_report.html" ]]; then
        log "📊 HTML Report: $OUTPUT_DIR/reports/phantom_analysis_report.html"
    fi
    
    if [[ -f "$OUTPUT_DIR/reports/workflow_summary.json" ]]; then
        log "📋 JSON Summary: $OUTPUT_DIR/reports/workflow_summary.json"
    fi
}

# Handle script arguments
case "${1:-}" in
    "check")
        check_dependencies
        exit 0
        ;;
    "validate")
        validate_inputs
        exit 0
        ;;
    "clean")
        log "Cleaning output directory..."
        rm -rf "$OUTPUT_DIR"
        success "Output directory cleaned"
        exit 0
        ;;
    "help"|"-h"|"--help")
        echo "Usage: $0 [command]"
        echo ""
        echo "Commands:"
        echo "  check     - Check dependencies"
        echo "  validate  - Validate inputs"
        echo "  clean     - Clean output directory"
        echo "  help      - Show this help message"
        echo ""
        echo "Run without arguments to execute the workflow"
        exit 0
        ;;
    "")
        main
        ;;
    *)
        error "Unknown command: $1"
        echo "Use '$0 help' for usage information"
        exit 1
        ;;
esac