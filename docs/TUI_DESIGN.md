# Phantom Science TUI Design

## Overview

A Terminal User Interface (TUI) for realistic sequencer simulation that mimics the experience of operating real Illumina sequencers.

## Main Interface Layout

```
┌─────────────────── Phantom Science Sequencer Simulator v0.1.0 ────────────────────┐
│                                                                                    │
│ ┌─── Sequencer Selection ───┐  ┌─── Run Status ────┐  ┌─── Real-time Stats ─────┐ │
│ │ ○ NovaSeq 6000 (S2)       │  │ Status: Ready      │  │ Files Generated: 0      │ │
│ │ ○ NextSeq 2000 (P3)       │  │ Progress: [    ]   │  │ Total Size: 0 MB       │ │
│ │ ● MiSeq v3 (300bp)        │  │ ETA: --:--         │  │ Current Rate: 0 MB/s   │ │
│ └───────────────────────────┘  └────────────────────┘  └─────────────────────────┘ │
│                                                                                    │
│ ┌─── Sample Sheet ─────────────────────────────────────────────────────────────┐   │
│ │ Sample_ID    Lane  Index_1    Index_2    Project      Description           │   │
│ │ Patient_001    1   ATCGATCG   CGTAGATC   Cancer_Study  Tumor sample        │   │
│ │ Patient_002    1   GCTAGCTA   TACGACGT   Cancer_Study  Normal sample       │   │
│ │ [Add Sample] [Import CSV] [Clear All]                                      │   │
│ └─────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                    │
│ ┌─── Run Parameters ──────────┐  ┌─── Output Settings ────────────────────────┐   │
│ │ Read Length: [150] bp        │  │ Output Dir: [./sequencer_output]           │   │
│ │ Paired-End: [✓]             │  │ Compression: [✓] gzip                     │   │
│ │ Chemistry: [v1.5]           │  │ Real-time: [✓] Stream files during run   │   │
│ │ PhiX %: [1.0]               │  │ Time Compression: [1x] [10x] [100x]      │   │
│ └─────────────────────────────┘  └────────────────────────────────────────────┘   │
│                                                                                    │
│ ┌─── File Generation Preview ──────────────────────────────────────────────────┐   │
│ │ Expected Files: 20 FASTQ files (~2.5 GB total)                              │   │
│ │ Naming: Patient_001_S1_L001_R1_001.fastq.gz                                │   │
│ │ Duration: ~36 hours (real) / 2.2 minutes (100x compressed)                 │   │
│ └──────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                    │
│ ┌─── Actions ─────────────────────────────────────────────────────────────────┐   │
│ │ [Start Run] [Load Config] [Save Config] [Stop] [Help] [Quit]               │   │
│ └─────────────────────────────────────────────────────────────────────────────┘   │
└────────────────────────────────────────────────────────────────────────────────────┘
```

## Key Features

### 1. Realistic Sequencer Selection
- **NovaSeq 6000**: High-throughput, multiple flow cell options (S1, S2, S4)
- **NextSeq 2000**: Mid-throughput, XLEAP-SBS chemistry
- **MiSeq v3**: Benchtop, longest reads (300bp)

### 2. Interactive Sample Sheet
- **Live Editing**: Add/remove samples in real-time
- **CSV Import**: Load existing sample sheets
- **Validation**: Check for duplicate barcodes, missing indexes
- **Auto-complete**: Suggest common index sequences

### 3. Real-time Simulation
- **Progress Tracking**: Files appear as they would during real sequencing
- **Realistic Timing**: Based on actual sequencer specifications
- **Time Compression**: 1x, 10x, 100x, 1000x modes for testing
- **File Streaming**: Generate files progressively, not all at once

### 4. Advanced Features
- **Live File Monitor**: Watch files appear in output directory
- **Network Simulation**: Mimic sequencer network behavior
- **Error Injection**: Simulate real sequencing errors/issues
- **Resource Monitoring**: CPU/Memory/Disk usage

## Navigation Controls

```
Tab/Shift+Tab:  Navigate between sections
Arrow Keys:     Navigate within sections  
Enter:          Select/Edit
Space:          Toggle checkboxes
Esc:            Cancel/Back
Ctrl+S:         Save configuration
Ctrl+Q:         Quit
F1:             Help
```

## Use Cases

### 1. Research Simulation
- Test downstream analysis pipelines
- Validate data handling workflows
- Train new lab members on sequencer operations

### 2. Software Development
- Test bioinformatics software with realistic data streams
- Validate file handling and parsing
- Performance testing with large datasets

### 3. Education
- Teach sequencing workflow without expensive hardware
- Demonstrate different sequencer characteristics
- Interactive learning tool for genomics courses

## Time Compression Examples

### NovaSeq 6000 (2×150bp, S2 Flow Cell)
- **Real Time**: 36 hours → 20 FASTQ files (100GB)
- **10x Compressed**: 3.6 hours → Same output
- **100x Compressed**: 22 minutes → Same output  
- **1000x Compressed**: 2.2 minutes → Same output

### File Generation Pattern (Real Time)
```
T+0:00    Run starts, cluster generation
T+2:30    First R1_001.fastq appears (Lane 1)
T+2:35    First R2_001.fastq appears (Lane 1)
T+5:00    Second batch files appear
...
T+36:00   Final files complete
```

## Configuration Integration

### Load Existing Configs
```bash
# Load from YAML
phantom-science --tui --config configs/sequencer/v1.0/novaseq-6000-s2.yaml

# Start with defaults
phantom-science --tui
```

### Save Custom Configs
- Export current TUI settings to YAML
- Save sample sheets as CSV
- Store time compression preferences
- Remember last used settings

## Implementation Architecture

```
TUI Layer (bubbletea/lipgloss)
├── Sequencer Models
├── Sample Sheet Manager  
├── Progress Tracker
├── File Generator (background)
└── Configuration Manager
```