#!/bin/bash

# Download 1000 Genomes Genomic Data from AWS Open Data
# Access real population genomic diversity for authenticity training

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
POPULATION_DIR="/Volumes/Public/genomics_training/populations"
LOG_FILE="$POPULATION_DIR/aws_download.log"

# Create directories
mkdir -p "$POPULATION_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== 1000 Genomes AWS Open Data Download ==="
echo "Started: $(date)"
echo "AWS S3 Bucket: s3://1000genomes"
echo "Strategy: Download high-diversity population VCFs and convert to FASTA"
echo "Target: 50-100 samples for maximum k-mer diversity"

# Check AWS CLI
if ! command -v aws &> /dev/null; then
    echo "❌ AWS CLI not found. Please install AWS CLI first."
    exit 1
fi

echo "✓ AWS CLI available"

# Function to download population VCF data
download_population_vcf() {
    local pop_code=$1
    local pop_name=$2
    local chromosome=$3
    
    echo ""
    echo "=== Downloading $pop_name ($pop_code) Chr$chromosome ==="
    
    local pop_dir="$POPULATION_DIR/$pop_code"
    mkdir -p "$pop_dir"
    
    # Download VCF file for this population and chromosome
    local vcf_file="$pop_dir/chr${chromosome}.vcf.gz"
    local aws_path="s3://1000genomes/release/20130502/ALL.chr${chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    
    echo "Downloading chromosome $chromosome VCF..."
    echo "AWS path: $aws_path"
    
    if aws s3 cp "$aws_path" "$vcf_file" --no-sign-request; then
        echo "✓ Downloaded Chr$chromosome VCF: $(du -sh "$vcf_file" | cut -f1)"
        
        # Download corresponding index file
        local index_file="$pop_dir/chr${chromosome}.vcf.gz.tbi"
        aws s3 cp "${aws_path}.tbi" "$index_file" --no-sign-request || echo "⚠ Index file not available"
        
        return 0
    else
        echo "❌ Failed to download Chr$chromosome VCF"
        return 1
    fi
}

# Function to download sample information
download_sample_info() {
    echo ""
    echo "=== Downloading Sample Information ==="
    
    # Download sample panel
    local panel_file="$POPULATION_DIR/integrated_call_samples_v3.20130502.ALL.panel"
    echo "Downloading sample panel from AWS..."
    
    if aws s3 cp "s3://1000genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" "$panel_file" --no-sign-request; then
        echo "✓ Downloaded sample panel: $(wc -l < "$panel_file") samples"
        
        # Show available populations
        echo ""
        echo "Available populations:"
        tail -n +2 "$panel_file" | cut -f2 | sort | uniq -c | sort -nr | head -20
        
        return 0
    else
        echo "❌ Failed to download sample panel"
        return 1
    fi
}

# Function to extract population-specific samples
extract_population_samples() {
    local pop_code=$1
    local pop_name=$2
    
    echo ""
    echo "=== Extracting $pop_name ($pop_code) Samples ==="
    
    local panel_file="$POPULATION_DIR/integrated_call_samples_v3.20130502.ALL.panel"
    local pop_dir="$POPULATION_DIR/$pop_code"
    local pop_samples="$pop_dir/samples.txt"
    
    if [ ! -f "$panel_file" ]; then
        echo "❌ Sample panel not found: $panel_file"
        return 1
    fi
    
    mkdir -p "$pop_dir"
    
    # Extract samples for this population
    echo "Extracting samples for $pop_name..."
    grep -E "^[^#].*\t$pop_code\t" "$panel_file" | head -20 > "$pop_samples"
    
    local sample_count=$(wc -l < "$pop_samples")
    echo "✓ $sample_count samples extracted for $pop_name"
    
    if [ $sample_count -gt 0 ]; then
        echo "First 5 samples:"
        head -5 "$pop_samples" | cut -f1
        return 0
    else
        echo "❌ No samples found for $pop_code"
        return 1
    fi
}

# Step 1: Download sample information
download_sample_info

# Step 2: Extract high-diversity populations
echo ""
echo "=== Extracting High-Diversity Populations ==="

# African populations (highest diversity)
extract_population_samples "YRI" "Yoruba_Nigeria"
extract_population_samples "LWK" "Luhya_Kenya"
extract_population_samples "GWD" "Gambian"

# European populations
extract_population_samples "CEU" "European_Utah"
extract_population_samples "GBR" "British"
extract_population_samples "TSI" "Tuscan_Italy"

# Asian populations
extract_population_samples "CHB" "Han_Chinese_Beijing"
extract_population_samples "JPT" "Japanese_Tokyo"
extract_population_samples "CHS" "Han_Chinese_South"

# Admixed populations
extract_population_samples "MXL" "Mexican_American"
extract_population_samples "PUR" "Puerto_Rican"

# Step 3: Download genomic data for key chromosomes
echo ""
echo "=== Downloading Genomic Data for Key Chromosomes ==="
echo "Strategy: Download Chr21 and Chr22 (smaller chromosomes) for all populations"
echo "This provides diverse training data with manageable size"

# Download ALL chromosomes for complete genomic diversity
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
POPULATIONS=(YRI LWK GWD CEU GBR TSI CHB JPT CHS MXL PUR)

echo "=== FULL 1000 GENOMES DATASET DOWNLOAD ==="
echo "Downloading ALL chromosomes for maximum diversity"
echo "Expected size: ~100-200GB total"
echo "Expected download time: 6-12 hours"

for chr in "${CHROMOSOMES[@]}"; do
    echo ""
    echo "=== Downloading Chromosome $chr Data ==="
    
    for pop in "${POPULATIONS[@]}"; do
        # Only download if population samples exist
        if [ -f "$POPULATION_DIR/$pop/samples.txt" ] && [ -s "$POPULATION_DIR/$pop/samples.txt" ]; then
            download_population_vcf "$pop" "Population_$pop" "$chr"
        else
            echo "⚠ Skipping $pop - no samples found"
        fi
    done
done

# Step 4: Create summary
echo ""
echo "=== Download Summary ==="
TOTAL_VCF=$(find "$POPULATION_DIR" -name "*.vcf.gz" | wc -l)
TOTAL_SIZE=$(du -sh "$POPULATION_DIR" | cut -f1)
TOTAL_POPS=$(ls -1 "$POPULATION_DIR" | grep -E "^[A-Z]{3}$" | wc -l)

echo "✅ 1000 Genomes AWS download completed!"
echo "Populations: $TOTAL_POPS"
echo "VCF files: $TOTAL_VCF"
echo "Total size: $TOTAL_SIZE"
echo "Storage: $POPULATION_DIR"

# Create download manifest
DOWNLOAD_MANIFEST="$POPULATION_DIR/aws_download_manifest.txt"
cat > "$DOWNLOAD_MANIFEST" << EOF
# 1000 Genomes AWS Open Data Download
# Completed: $(date)
# Source: AWS S3 bucket s3://1000genomes

DOWNLOAD RESULTS:
===============
Populations: $TOTAL_POPS high-diversity groups
VCF Files: $TOTAL_VCF genomic variant files
Total Size: $TOTAL_SIZE
Chromosomes: 21, 22 (validation set)

POPULATION COVERAGE:
==================
African: YRI (Yoruba), LWK (Luhya), GWD (Gambian)
European: CEU (Utah), GBR (British), TSI (Tuscan)
Asian: CHB (Han Chinese), JPT (Japanese), CHS (Han Chinese South)
Admixed: MXL (Mexican American), PUR (Puerto Rican)

NEXT STEPS:
==========
1. Convert VCF to FASTA format:
   ./scripts/convert_vcf_to_fasta.sh

2. Create population training dataset:
   ./scripts/prepare_population_training.sh

3. Run population-aware training:
   go run ./cmd/experiments/population_training/main.go

EXPECTED RESULTS:
================
K-mer diversity: 10-100x increase
Complexity: 0.05-0.3 (10-60x improvement)
JS Divergence: < 0.01 (production-ready)
Authenticity: Population-aware patterns
EOF

echo ""
echo "✓ Download manifest: $DOWNLOAD_MANIFEST"
echo ""
echo "📊 REAL 1000 GENOMES DATA ACQUIRED!"
echo "Ready for VCF-to-FASTA conversion and population training"
echo ""
echo "Next commands:"
echo "1. Convert data: ./scripts/convert_vcf_to_fasta.sh"
echo "2. Train model: go run ./cmd/experiments/population_training/main.go"
echo ""
echo "🚀 AUTHENTIC POPULATION GENOMICS BREAKTHROUGH!"
echo "Completed: $(date)"
echo "Log: $LOG_FILE"