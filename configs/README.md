# Configuration System

Phantom Science uses a versioned configuration system supporting both JSON and YAML formats.

## Configuration Versions

### v1.0 (Current)
- **Location**: `configs/sequencer/v1.0/` and `configs/phantom/v1.0/`
- **Format**: YAML (recommended) or JSON
- **Features**: Full sequencer device specifications with vendor-validated parameters

## File Formats

### YAML (Recommended)
- **Extension**: `.yaml` or `.yml`
- **Benefits**: Comments, better readability, industry standard in bioinformatics
- **Schema validation**: Built-in version validation

### JSON (Legacy Support)
- **Extension**: `.json`
- **Benefits**: Compact, programmatic generation
- **Compatibility**: Backward compatible with existing configs

## Configuration Types

### Sequencer Configurations
Located in `configs/sequencer/v1.0/`

#### Available Platforms
- **`novaseq-6000-s2.yaml`**: High-throughput WGS (3TB output, 36h runtime)
- **`nextseq-2000-p3.yaml`**: Mid-throughput applications (120GB output, 18h runtime)  
- **`miseq-v3.yaml`**: Long-read amplicon sequencing (15GB output, 65h runtime, 300bp reads)

#### Schema Structure
```yaml
schema_version: "1.0"
last_updated: "2024-07-03"
validated_against: "Illumina [Platform] Specifications"
config_type: "sequencer"

sequencer:
  model: "NovaSeq 6000"
  flow_cell_type: "S2"
  # ... hardware specifications

sample:
  name: "Sample_001" 
  # ... sample configuration

reads:
  length: 150
  # ... read parameters

quality:
  base_quality: 35
  # ... quality metrics

performance:
  run_time_hours: 36
  # ... performance specifications

technology:
  sbs_version: "v1.5"
  # ... technology features
```

### Phantom Science Configurations  
Located in `configs/phantom/v1.0/`

#### Available Configurations
- **`enterprise-settings.yaml`**: Full enterprise features with compliance
- **`volume-testing.yaml`**: Optimized for generating millions of files

#### Schema Structure
```yaml
schema_version: "1.0"
config_type: "phantom_enterprise"

enterprise:
  customer_id: "CUST_001"
  # ... customer information

fingerprinting:
  enable: true
  # ... steganography settings

audit_logging:
  enable: true
  # ... compliance logging

validation:
  enable: true
  # ... statistical validation

regulatory_compliance:
  fda_compliance: true
  # ... regulatory settings
```

## Parameter Validation

All configurations are validated against vendor specifications:

### NovaSeq 6000 Specifications
- **Read Lengths**: 50, 100, 150, 250bp
- **Q30**: ≥90% guaranteed
- **Run Times**: 13-44 hours depending on read length
- **Output**: Up to 6TB per run
- **Chemistry**: SBS v1.5

### NextSeq 2000 Specifications  
- **Read Lengths**: 1×50 to 2×300bp (14 configurations)
- **Q30**: ≥90% for shorter reads
- **Chemistry**: XLEAP-SBS (2-channel)
- **Output**: Variable by flow cell (P1-P4)

### MiSeq Specifications
- **Read Lengths**: Up to 2×300bp (longest Illumina reads)
- **Q30**: ≥80% typical for 300bp
- **Chemistry**: v3 (4-channel)
- **Output**: Up to 15GB per run

## Usage Examples

### Basic Sequencer Simulation
```bash
# YAML format (recommended)
./phantom-science -config configs/sequencer/v1.0/miseq-v3.yaml

# JSON format (legacy)
./phantom-science -config configs/sequencer/simple-run.json
```

### Enterprise Features
```bash
# Full enterprise configuration
./phantom-science \
  -config configs/sequencer/v1.0/novaseq-6000-s2.yaml \
  -phantom configs/phantom/v1.0/enterprise-settings.yaml

# Volume testing (optimized for speed)
./phantom-science \
  -config configs/sequencer/v1.0/nextseq-2000-p3.yaml \
  -phantom configs/phantom/v1.0/volume-testing.yaml
```

## Migration Guide

### From JSON to YAML
The system automatically detects and loads both formats. To migrate:

1. **Copy existing JSON config**
2. **Convert to YAML format** with schema version
3. **Add comments** documenting parameter sources
4. **Validate** against vendor specifications

### Version Compatibility
- **v1.0 configs**: Full feature support
- **Legacy configs**: Automatically converted to v1.0 internally
- **Mixed usage**: JSON and YAML can be used together

## Adding New Sequencer Models

To add a new sequencer model:

1. **Research vendor specifications**
2. **Create new YAML config** in `configs/sequencer/v1.0/`
3. **Document parameter sources** in comments
4. **Validate against official specs**
5. **Test with actual generation**

Example template:
```yaml
# [Platform] Configuration
# Based on [Vendor] [Platform] Specifications ([Year])
# Spec Reference: [URL]

schema_version: "1.0"
last_updated: "[DATE]"
validated_against: "[VENDOR SPEC]"
config_type: "sequencer"

# Add platform-specific configurations following the schema
```

## Best Practices

### Configuration Comments
- **Document parameter sources**: Reference vendor specification sheets
- **Explain non-obvious values**: Loading concentrations, PhiX percentages
- **Note technology differences**: 2-channel vs 4-channel chemistry

### Version Management
- **Always specify schema_version**: Enables proper validation
- **Update last_updated field**: Track configuration changes
- **Validate against vendor specs**: Ensure parameter accuracy

### File Organization
- **Use descriptive names**: `novaseq-6000-s2.yaml` vs `config1.yaml`
- **Separate by version**: Keep different schema versions in separate directories
- **Group by type**: Sequencer vs Phantom configurations