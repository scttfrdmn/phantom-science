#!/bin/bash
# Phantom Science GPU Performance Benchmark Suite
# Enterprise-grade performance validation

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Configuration
BENCHMARK_DIR="/tmp/phantom-benchmark-$(date +%s)"
LOG_FILE="$BENCHMARK_DIR/benchmark-results.log"
RESULTS_JSON="$BENCHMARK_DIR/results.json"

# Test configurations
declare -A TEST_CONFIGS=(
    ["small"]="1MB,1000"
    ["medium"]="100MB,10000"  
    ["large"]="1GB,100000"
    ["xlarge"]="10GB,1000000"
)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_FILE"
}

success() {
    echo -e "${GREEN}✅ $1${NC}" | tee -a "$LOG_FILE"
}

warning() {
    echo -e "${YELLOW}⚠️  $1${NC}" | tee -a "$LOG_FILE"
}

error() {
    echo -e "${RED}❌ $1${NC}" | tee -a "$LOG_FILE"
}

setup_benchmark() {
    log "Setting up benchmark environment..."
    mkdir -p "$BENCHMARK_DIR"
    
    # Check prerequisites
    if ! command -v nvidia-smi &> /dev/null; then
        error "NVIDIA GPU not found. GPU benchmarks require CUDA-capable hardware."
        exit 1
    fi
    
    # Check if phantom-science binary exists
    if [ ! -f "$PROJECT_ROOT/phantom-science" ]; then
        log "Building phantom-science binary..."
        cd "$PROJECT_ROOT"
        go build -o phantom-science .
    fi
    
    success "Benchmark environment ready"
}

check_gpu_status() {
    log "Checking GPU status..."
    
    # Get GPU info
    GPU_COUNT=$(nvidia-smi --query-gpu=count --format=csv,noheader,nounits | head -1)
    GPU_MEMORY=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
    
    log "GPU Hardware: $GPU_NAME"
    log "GPU Count: $GPU_COUNT"
    log "GPU Memory: ${GPU_MEMORY}MB"
    
    # Check CUDA
    if command -v nvcc &> /dev/null; then
        CUDA_VERSION=$(nvcc --version | grep "release" | sed 's/.*release //' | sed 's/,.*$//')
        log "CUDA Version: $CUDA_VERSION"
    fi
    
    success "GPU status check complete"
}

run_performance_test() {
    local test_name="$1"
    local config="$2"
    local size="${config%,*}"
    local reads="${config#*,}"
    
    log "Running $test_name test: $size data, $reads reads"
    
    # Create test configuration
    local test_config="$BENCHMARK_DIR/config-$test_name.yaml"
    cat > "$test_config" << EOF
schema_version: "1.0"
config_type: "sequencer"

sequencer:
  model: "NovaSeq 6000"
  flow_cell_type: "S2"

sample:
  name: "benchmark-$test_name"
  lane: 1
  output_dir: "$BENCHMARK_DIR/output-$test_name"

reads:
  length: 150
  reads_per_batch: $reads
  total_batches: 1

enterprise:
  customer_id: "benchmark-customer"
  fingerprinting:
    enable: true
    strategy: "balanced"
EOF

    # Run benchmark
    local start_time=$(date +%s.%N)
    
    "$PROJECT_ROOT/phantom-science" \
        --config "$test_config" \
        --command generate \
        --verbose 2>&1 | tee -a "$LOG_FILE"
    
    local end_time=$(date +%s.%N)
    local duration=$(echo "$end_time - $start_time" | bc)
    
    # Calculate performance metrics
    local output_dir="$BENCHMARK_DIR/output-$test_name"
    if [ -d "$output_dir" ]; then
        local total_size=$(du -sb "$output_dir" | cut -f1)
        local throughput=$(echo "scale=2; $total_size / 1024 / 1024 / $duration" | bc)
        
        success "$test_name: ${duration}s, ${throughput}MB/s"
        
        # Store results
        echo "\"$test_name\": {\"duration\": $duration, \"throughput\": $throughput, \"size\": $total_size}," >> "$RESULTS_JSON.tmp"
    else
        error "$test_name: Output directory not found"
    fi
}

run_gpu_memory_test() {
    log "Running GPU memory utilization test..."
    
    # Monitor GPU memory during large dataset processing
    nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader,nounits > "$BENCHMARK_DIR/gpu-memory-before.txt"
    
    run_performance_test "gpu-stress" "5GB,500000"
    
    nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader,nounits > "$BENCHMARK_DIR/gpu-memory-after.txt"
    
    success "GPU memory test complete"
}

generate_report() {
    log "Generating benchmark report..."
    
    # Finalize JSON results
    echo "{" > "$RESULTS_JSON"
    if [ -f "$RESULTS_JSON.tmp" ]; then
        sed '$ s/,$//' "$RESULTS_JSON.tmp" >> "$RESULTS_JSON"
        rm "$RESULTS_JSON.tmp"
    fi
    echo "}" >> "$RESULTS_JSON"
    
    # Generate summary report
    local report_file="$BENCHMARK_DIR/benchmark-report.md"
    cat > "$report_file" << EOF
# Phantom Science GPU Benchmark Report

**Date**: $(date)
**Hardware**: $GPU_NAME ($GPU_COUNT GPUs, ${GPU_MEMORY}MB each)
**Software**: Phantom Science v$(cd "$PROJECT_ROOT" && git describe --tags --always)

## Performance Results

| Test Case | Duration (s) | Throughput (MB/s) | Status |
|-----------|--------------|-------------------|---------|
EOF
    
    # Parse results and add to table
    if [ -f "$RESULTS_JSON" ]; then
        python3 -c "
import json, sys
try:
    with open('$RESULTS_JSON') as f:
        data = json.load(f)
    for test, metrics in data.items():
        if metrics:
            print(f'| {test} | {metrics[\"duration\"]:.2f} | {metrics[\"throughput\"]:.2f} | ✅ |')
except:
    pass
" >> "$report_file"
    fi
    
    echo "" >> "$report_file"
    echo "## Detailed Logs" >> "$report_file"
    echo "" >> "$report_file"
    echo "\`\`\`" >> "$report_file"
    tail -50 "$LOG_FILE" >> "$report_file"
    echo "\`\`\`" >> "$report_file"
    
    success "Benchmark report generated: $report_file"
    success "Detailed results: $RESULTS_JSON"
}

main() {
    log "Starting Phantom Science GPU Benchmark Suite"
    
    setup_benchmark
    check_gpu_status
    
    # Run performance tests
    for test_name in "${!TEST_CONFIGS[@]}"; do
        run_performance_test "$test_name" "${TEST_CONFIGS[$test_name]}"
    done
    
    # Run specialized tests
    run_gpu_memory_test
    
    generate_report
    
    success "Benchmark suite completed successfully!"
    log "Results available in: $BENCHMARK_DIR"
}

# Execute main function
main "$@"