# Phantom Science Production Deployment

## Enterprise-Grade Deployment Pipeline

### Prerequisites

- **GPU Hardware**: NVIDIA TITAN RTX or better (24GB+ VRAM)
- **CUDA**: Version 12.9+
- **Go**: Version 1.21+
- **Rust**: Version 1.88+ with cargo
- **OS**: Linux (Rocky 9+ recommended for enterprise)

### Quick Deployment

```bash
# 1. Clone and build
git clone git@github.com:scttfrdmn/phantom-science.git
cd phantom-science

# 2. Build stegano library with GPU support
cd stegano
cargo build --release --features gpu,metal
cd ..

# 3. Build main application
go build -o phantom-science .

# 4. Verify installation
./phantom-science --version
```

### Production Configuration

```yaml
# deploy/production/phantom-enterprise.yaml
schema_version: "1.0"
config_type: "enterprise"

enterprise:
  customer_id: "${CUSTOMER_ID}"
  license_key: "${LICENSE_KEY}"
  service_tier: "enterprise"
  compliance_level: "regulatory"
  
gpu_acceleration:
  enable: true
  memory_limit_gb: 20  # Leave 4GB for system
  parallel_streams: 4
  
monitoring:
  telemetry: true
  performance_logging: true
  health_checks: true
  alert_thresholds:
    memory_usage: 85
    gpu_utilization: 90
```

### Security Hardening

```bash
# Set up service account
sudo useradd -r -s /bin/false phantom-science

# Secure file permissions
sudo chown -R phantom-science:phantom-science /opt/phantom-science
sudo chmod 750 /opt/phantom-science
sudo chmod 640 /opt/phantom-science/config/*.yaml

# GPU access
sudo usermod -a -G video phantom-science
```

### Performance Verification

```bash
# Run performance benchmarks
./phantom-science benchmark --gpu --dataset-size 10GB
# Expected: >400MB/s extraction, >176MB/s embedding

# Health check
./phantom-science health --comprehensive
# Expected: All systems operational, GPU available
```

### Enterprise Support

- **24/7 Support**: support@phantom-science.com
- **Security Incidents**: security@phantom-science.com  
- **Performance Issues**: performance@phantom-science.com
- **Documentation**: https://docs.phantom-science.com/enterprise