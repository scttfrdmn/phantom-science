# Sequencer Configuration Files

This directory contains configuration files for real Illumina sequencer models with accurate parameters based on vendor specifications.

## Available Configurations

### NovaSeq 6000 (S2 Flow Cell)
- **File**: `novaseq-6000-s2.json`
- **Read Length**: 2×150bp
- **Output**: Up to 3TB per run
- **Run Time**: ~36 hours
- **Q30**: ≥90%
- **Use Case**: High-throughput whole genome sequencing

### NextSeq 2000 (P3 Flow Cell)
- **File**: `nextseq-2000-p3.json`
- **Read Length**: 2×150bp
- **Output**: Up to 120GB per run
- **Run Time**: ~18 hours
- **Q30**: ≥90%
- **Chemistry**: XLEAP-SBS (2-channel)
- **Use Case**: Mid-throughput applications, exome sequencing

### MiSeq v3
- **File**: `miseq-v3.json`
- **Read Length**: 2×300bp (longest available)
- **Output**: Up to 15GB per run
- **Run Time**: ~65 hours
- **Q30**: ≥80%
- **Chemistry**: v3 (4-channel)
- **Use Case**: Amplicon sequencing, targeted panels, small genomes

### Generic Illumina
- **File**: `simple-run.json`
- **Read Length**: 2×150bp
- **Use Case**: Testing and development

## Parameter Accuracy

All configurations are based on official Illumina specifications from 2024:

- **Cluster Densities**: Match vendor specifications for each flow cell type
- **Loading Concentrations**: Realistic pM concentrations used in practice
- **Run Times**: Actual run durations for different read lengths
- **Quality Metrics**: Real Q30 percentages achievable
- **Chemistry Types**: Accurate 2-channel vs 4-channel chemistry flags
- **Output Capacities**: Maximum throughput per flow cell type

## Usage

```bash
# NovaSeq 6000 high-throughput run
./phantom-science -config configs/sequencer/novaseq-6000-s2.json

# NextSeq 2000 mid-throughput run  
./phantom-science -config configs/sequencer/nextseq-2000-p3.json

# MiSeq long-read amplicon sequencing
./phantom-science -config configs/sequencer/miseq-v3.json
```

## Key Differences Between Platforms

| Platform | Read Length | Output | Run Time | Chemistry | Best For |
|----------|-------------|--------|----------|-----------|----------|
| NovaSeq 6000 | 2×150bp | 3TB | 36h | SBS v1.5 | WGS, large cohorts |
| NextSeq 2000 | 2×150bp | 120GB | 18h | XLEAP-SBS | Exomes, targeted |
| MiSeq v3 | 2×300bp | 15GB | 65h | v3 | Amplicons, validation |

## Parameter Sources

- Illumina NovaSeq 6000 System Specifications
- NextSeq 1000/2000 System Specifications  
- MiSeq System Specifications
- Illumina Knowledge Base (2024 updates)