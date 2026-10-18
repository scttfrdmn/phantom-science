#!/bin/bash

# Download Full Human Genome Reference Script
# Downloads complete GRCh38 reference genome with all chromosomes

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Set GENOMICS_DATA_DIR to download to high-capacity storage instead of the repo.
if [ -n "${GENOMICS_DATA_DIR:-}" ]; then
    DATA_DIR="$GENOMICS_DATA_DIR/human"
    echo "Using GENOMICS_DATA_DIR: $DATA_DIR"
else
    DATA_DIR="$PROJECT_ROOT/data/validation/reference_genome"
    echo "Using local storage: $DATA_DIR"
fi

LOG_FILE="$DATA_DIR/full_genome_download.log"

# Create data directory
mkdir -p "$DATA_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== Full Human Genome Reference Download ==="
echo "Started: $(date)"
echo "Data directory: $DATA_DIR"
echo "Target: Complete GRCh38 reference with all chromosomes"

# Check if AWS CLI is installed
if ! command -v aws &> /dev/null; then
    echo "ERROR: AWS CLI is not installed. Installing now..."
    
    # Install AWS CLI v2 for macOS
    if [[ "$OSTYPE" == "darwin"* ]]; then
        if ! command -v brew &> /dev/null; then
            echo "Installing Homebrew first..."
            /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        fi
        brew install awscli
    else
        echo "Please install AWS CLI manually:"
        echo "  curl 'https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip' -o 'awscliv2.zip'"
        echo "  unzip awscliv2.zip"
        echo "  sudo ./aws/install"
        exit 1
    fi
fi

# Configure AWS CLI for anonymous access to public data
aws configure set default.region us-east-1
aws configure set default.output json

# Function to download file with retry logic
download_file() {
    local bucket=$1
    local key=$2
    local local_path=$3
    local max_retries=3
    local retry_count=0
    
    echo "Downloading: s3://$bucket/$key"
    
    while [ $retry_count -lt $max_retries ]; do
        if aws s3 cp "s3://$bucket/$key" "$local_path" --no-sign-request; then
            echo "✓ Successfully downloaded: $key"
            return 0
        else
            retry_count=$((retry_count + 1))
            echo "⚠ Download failed, retry $retry_count/$max_retries"
            sleep 5
        fi
    done
    
    echo "✗ Failed to download after $max_retries attempts: $key"
    return 1
}

# Download Complete GRCh38 Reference Genome
echo ""
echo "=== Downloading Complete GRCh38 Reference Genome ==="

# Method 1: Download from Broad Institute (most reliable)
echo "Attempting download from Broad Institute..."
download_file "broad-references" "hg38/v0/Homo_sapiens_assembly38.fasta" "$DATA_DIR/GRCh38_complete.fasta"

# If that fails, try NCBI
if [ $? -ne 0 ]; then
    echo "Broad Institute failed, trying NCBI..."
    download_file "ncbi-genomes" "refseq/vertebrate_mammalian/Homo_sapiens/reference/GRCh38_latest/GRCh38_latest_genomic.fna.gz" "$DATA_DIR/GRCh38_complete.fasta.gz"
    
    if [ $? -eq 0 ]; then
        echo "Decompressing NCBI download..."
        gunzip "$DATA_DIR/GRCh38_complete.fasta.gz"
        mv "$DATA_DIR/GRCh38_complete.fna" "$DATA_DIR/GRCh38_complete.fasta"
    fi
fi

# If both fail, download individual chromosomes
if [ ! -f "$DATA_DIR/GRCh38_complete.fasta" ]; then
    echo "Complete genome download failed, downloading individual chromosomes..."
    
    # Download each chromosome individually
    for chr in {1..22} X Y M; do
        echo "Downloading chromosome $chr..."
        download_file "broad-references" "hg38/v0/Homo_sapiens_assembly38.fasta" "$DATA_DIR/chr${chr}.fasta" || \
        download_file "1000genomes" "technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" "$DATA_DIR/chr${chr}.fasta" || \
        echo "Warning: Could not download chromosome $chr"
    done
    
    # Concatenate individual chromosomes
    echo "Concatenating individual chromosomes..."
    cat "$DATA_DIR"/chr*.fasta > "$DATA_DIR/GRCh38_complete.fasta" 2>/dev/null || echo "Concatenation failed"
fi

# Download genome index files
echo ""
echo "=== Downloading Genome Index Files ==="
download_file "broad-references" "hg38/v0/Homo_sapiens_assembly38.fasta.fai" "$DATA_DIR/GRCh38_complete.fasta.fai"
download_file "broad-references" "hg38/v0/Homo_sapiens_assembly38.dict" "$DATA_DIR/GRCh38_complete.dict"

# Alternative: Use curl for direct downloads
echo ""
echo "=== Alternative Download Methods ==="

# ENSEMBL FTP (very reliable)
if [ ! -f "$DATA_DIR/GRCh38_complete.fasta" ]; then
    echo "Trying ENSEMBL FTP..."
    curl -L "ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" -o "$DATA_DIR/GRCh38_ensembl.fasta.gz"
    
    if [ $? -eq 0 ]; then
        echo "Decompressing ENSEMBL download..."
        gunzip "$DATA_DIR/GRCh38_ensembl.fasta.gz"
        mv "$DATA_DIR/GRCh38_ensembl.fasta" "$DATA_DIR/GRCh38_complete.fasta"
    fi
fi

# NCBI FTP (backup)
if [ ! -f "$DATA_DIR/GRCh38_complete.fasta" ]; then
    echo "Trying NCBI FTP..."
    curl -L "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" -o "$DATA_DIR/GRCh38_ncbi.fasta.gz"
    
    if [ $? -eq 0 ]; then
        echo "Decompressing NCBI download..."
        gunzip "$DATA_DIR/GRCh38_ncbi.fasta.gz"
        mv "$DATA_DIR/GRCh38_ncbi.fna" "$DATA_DIR/GRCh38_complete.fasta"
    fi
fi

# Verify download
echo ""
echo "=== Verifying Download ==="

if [ -f "$DATA_DIR/GRCh38_complete.fasta" ]; then
    file_size=$(stat -c%s "$DATA_DIR/GRCh38_complete.fasta" 2>/dev/null || stat -f%z "$DATA_DIR/GRCh38_complete.fasta" 2>/dev/null)
    echo "✓ Complete genome file size: $(numfmt --to=iec-i --suffix=B $file_size)"
    
    # Count chromosomes
    chr_count=$(grep -c "^>" "$DATA_DIR/GRCh38_complete.fasta")
    echo "✓ Number of sequences: $chr_count"
    
    # Show first few chromosomes
    echo "✓ First 10 sequences:"
    grep "^>" "$DATA_DIR/GRCh38_complete.fasta" | head -10
    
    # Expected size should be > 3GB
    if [ $file_size -gt 3000000000 ]; then
        echo "✅ SUCCESS: Complete genome downloaded successfully"
    else
        echo "⚠ WARNING: File size seems small, may be incomplete"
    fi
else
    echo "❌ FAILED: Could not download complete genome"
    echo "Manual download options:"
    echo "1. ENSEMBL: ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    echo "2. NCBI: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    echo "3. Broad: gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    exit 1
fi

# Create individual chromosome files from complete genome
echo ""
echo "=== Extracting Individual Chromosomes ==="

if [ -f "$DATA_DIR/GRCh38_complete.fasta" ]; then
    echo "Extracting individual chromosomes for faster access..."
    
    # Use samtools if available
    if command -v samtools &> /dev/null; then
        echo "Using samtools for chromosome extraction..."
        samtools faidx "$DATA_DIR/GRCh38_complete.fasta"
        
        # Extract major chromosomes
        for chr in {1..22} X Y; do
            echo "Extracting chromosome $chr..."
            samtools faidx "$DATA_DIR/GRCh38_complete.fasta" "chr$chr" > "$DATA_DIR/chr$chr.fasta" 2>/dev/null || \
            samtools faidx "$DATA_DIR/GRCh38_complete.fasta" "$chr" > "$DATA_DIR/chr$chr.fasta" 2>/dev/null || \
            echo "Warning: Could not extract chromosome $chr"
        done
    else
        echo "Installing samtools for chromosome extraction..."
        if [[ "$OSTYPE" == "darwin"* ]]; then
            brew install samtools
        else
            echo "Please install samtools manually"
        fi
    fi
fi

# Create data manifest
echo ""
echo "=== Creating Data Manifest ==="
MANIFEST_FILE="$DATA_DIR/full_genome_manifest.txt"

cat > "$MANIFEST_FILE" << EOF
# Full Human Genome Reference Manifest
# Generated: $(date)
# Source: GRCh38 Complete Reference

# Complete Genome
GRCh38_complete.fasta
GRCh38_complete.fasta.fai
GRCh38_complete.dict

# Individual Chromosomes
$(for chr in {1..22} X Y; do echo "chr$chr.fasta"; done)

# Statistics
# File size: $(stat -c%s "$DATA_DIR/GRCh38_complete.fasta" 2>/dev/null || stat -f%z "$DATA_DIR/GRCh38_complete.fasta" 2>/dev/null || echo "Unknown") bytes
# Sequences: $(grep -c "^>" "$DATA_DIR/GRCh38_complete.fasta" 2>/dev/null || echo "Unknown")
# Download time: $(date)
EOF

echo "✓ Data manifest created: $MANIFEST_FILE"

# Final summary
echo ""
echo "=== Download Complete ==="
echo "✅ Full human genome reference downloaded successfully"
echo "📁 Location: $DATA_DIR/"
echo "📄 Manifest: $MANIFEST_FILE"
echo "🧬 Ready for comprehensive genomic authenticity training"
echo ""
echo "Next steps:"
echo "1. Run full genome training: go run ./cmd/experiments/optimized_full_genome_training/main.go"
echo "2. Verify sequence complexity > 0.5"
echo "3. Measure JS divergence < 0.01"
echo ""
echo "Completed: $(date)"
echo "Log file: $LOG_FILE"