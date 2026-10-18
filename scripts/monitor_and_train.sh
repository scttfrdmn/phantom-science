#!/bin/bash

# Monitor 1000 Genomes Download and Start High-Performance Training
# Automatically detects when download is complete and launches optimized training

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
POPULATION_DIR="/Volumes/Public/genomics_training/populations"
TRAINING_BINARY="$PROJECT_ROOT/high_performance_training"
LOG_FILE="$PROJECT_ROOT/monitor_and_train.log"
LOCK_FILE="$PROJECT_ROOT/.training_in_progress"

# Create log file
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== 1000 GENOMES DOWNLOAD MONITOR & HIGH-PERFORMANCE TRAINING ==="
echo "Started: $(date)"
echo "Apple Silicon P-core + Metal GPU optimization enabled"
echo "Target: Complete 1000 Genomes dataset (100-200GB)"
echo "Expected training time: 4-8 hours with high-performance acceleration"
echo ""

# Function to check download progress
check_download_progress() {
    local current_size=$(du -sh "$POPULATION_DIR" 2>/dev/null | cut -f1 || echo "0")
    local vcf_files=$(find "$POPULATION_DIR" -name "*.vcf.gz" 2>/dev/null | wc -l || echo "0")
    local download_processes=$(ps aux | grep download_1000genomes_data.sh | grep -v grep | wc -l || echo "0")
    
    echo "📊 Download Status:"
    echo "  Current size: $current_size"
    echo "  VCF files: $vcf_files"
    echo "  Download processes: $download_processes"
    echo ""
    
    # Check if download is complete
    if [[ $download_processes -eq 0 ]] && [[ $vcf_files -gt 200 ]]; then
        echo "✅ Download appears complete!"
        echo "📈 Final size: $current_size"
        echo "📁 Total VCF files: $vcf_files"
        return 0
    else
        echo "⏳ Download still in progress..."
        return 1
    fi
}

# Function to start high-performance training
start_training() {
    echo ""
    echo "🚀 STARTING HIGH-PERFORMANCE TRAINING"
    echo "====================================="
    echo "Training method: Apple Silicon P-cores + Metal GPU"
    echo "Dataset: Complete 1000 Genomes (11 populations)"
    echo "Expected duration: 4-8 hours"
    echo "Expected speedup: 12-32x over single-core"
    echo ""
    
    # Create lock file
    touch "$LOCK_FILE"
    
    # Change to project directory
    cd "$PROJECT_ROOT"
    
    # Start high-performance training
    echo "Launching high-performance training..."
    if "$TRAINING_BINARY" > high_performance_training.log 2>&1; then
        echo "✅ HIGH-PERFORMANCE TRAINING COMPLETED SUCCESSFULLY!"
        echo "📊 Results saved in: high_performance_progress/"
        echo "📈 Training log: high_performance_training.log"
        
        # Show final results
        if [[ -f "high_performance_progress/final_results.txt" ]]; then
            echo ""
            echo "🎉 FINAL RESULTS:"
            cat "high_performance_progress/final_results.txt"
        fi
    else
        echo "❌ Training failed. Check high_performance_training.log for details."
        return 1
    fi
    
    # Remove lock file
    rm -f "$LOCK_FILE"
    
    return 0
}

# Function to show system info
show_system_info() {
    echo "🖥️  SYSTEM INFORMATION:"
    echo "Platform: $(uname -s) $(uname -r)"
    echo "Architecture: $(uname -m)"
    echo "Total CPU cores: $(sysctl -n hw.ncpu 2>/dev/null || echo "unknown")"
    echo "Physical cores: $(sysctl -n hw.physicalcpu 2>/dev/null || echo "unknown")"
    echo "Available memory: $(echo "scale=1; $(sysctl -n hw.memsize 2>/dev/null || echo "0") / 1024 / 1024 / 1024" | bc 2>/dev/null || echo "unknown")GB"
    echo "Available disk: $(df -h "$POPULATION_DIR" 2>/dev/null | tail -1 | awk '{print $4}' || echo "unknown")"
    echo ""
}

# Function to estimate completion time
estimate_completion() {
    local current_size_bytes=$(du -b "$POPULATION_DIR" 2>/dev/null | cut -f1 || echo "0")
    local target_size_bytes=$((100 * 1024 * 1024 * 1024))  # 100GB minimum
    
    if [[ $current_size_bytes -gt 0 ]]; then
        local progress_pct=$(echo "scale=1; $current_size_bytes * 100 / $target_size_bytes" | bc 2>/dev/null || echo "0")
        echo "📈 Download progress: ${progress_pct}%"
        
        if [[ $(echo "$progress_pct > 0" | bc 2>/dev/null) -eq 1 ]]; then
            # Estimate time remaining (very rough)
            local hours_remaining=$(echo "scale=1; (100 - $progress_pct) / 10" | bc 2>/dev/null || echo "unknown")
            echo "⏱️  Estimated time remaining: ~${hours_remaining} hours"
        fi
    fi
    echo ""
}

# Check if training is already in progress
if [[ -f "$LOCK_FILE" ]]; then
    echo "⚠️  Training already in progress (lock file exists)"
    echo "If training is not running, remove: $LOCK_FILE"
    exit 1
fi

# Check if training binary exists
if [[ ! -f "$TRAINING_BINARY" ]]; then
    echo "🔨 Building high-performance training binary..."
    cd "$PROJECT_ROOT"
    if go build -o high_performance_training ./cmd/experiments/high_performance_training; then
        echo "✅ Training binary built successfully"
    else
        echo "❌ Failed to build training binary"
        exit 1
    fi
fi

# Show system information
show_system_info

# Main monitoring loop
echo "🔍 Starting download monitoring..."
echo "Press Ctrl+C to stop monitoring"
echo ""

while true; do
    echo "=== $(date) ==="
    
    # Check if download directory exists
    if [[ ! -d "$POPULATION_DIR" ]]; then
        echo "❌ Population directory not found: $POPULATION_DIR"
        echo "Make sure the 1000 Genomes download is running"
        sleep 60
        continue
    fi
    
    # Check download progress
    if check_download_progress; then
        echo ""
        echo "🎯 DOWNLOAD COMPLETE! Starting high-performance training..."
        echo ""
        
        # Start training
        if start_training; then
            echo ""
            echo "🏆 SUCCESS! High-performance population genomics training completed!"
            echo "Results available in: high_performance_progress/"
            break
        else
            echo "❌ Training failed. Monitoring will continue..."
        fi
    else
        # Show progress estimation
        estimate_completion
    fi
    
    echo "Next check in 5 minutes..."
    echo ""
    sleep 300  # Check every 5 minutes
done

echo ""
echo "=== MONITOR AND TRAIN COMPLETE ==="
echo "Completed: $(date)"
echo "Log file: $LOG_FILE"
echo ""
echo "🚀 HIGH-PERFORMANCE GENOMICS TRAINING PIPELINE COMPLETE!"
echo "Ready for expert validation and production deployment"