# Development Guidelines

This document outlines the development standards and practices for the Phantom Science project.

## Version

**Current Version**: 0.1.0

## Development Rules

### 1. Code Coverage Requirements

**Mandatory coverage standards that must be met before any code is merged:**

- **Overall Project Coverage**: Minimum 85%
- **Individual File Coverage**: Minimum 80% per file
- **Critical Components**: 90%+ coverage required for:
  - `internal/genomics/` (core FASTQ generation)
  - `internal/fingerprinting/` (security features)
  - `internal/validation/` (quality assurance)
  - `pkg/types/` (configuration handling)

### 2. Testing Standards

- **Unit Tests**: Required for all public functions
- **Integration Tests**: Required for CLI commands and config loading
- **Performance Tests**: Required for file generation at scale
- **Security Tests**: Required for fingerprinting and encryption features

### 3. Pre-commit Requirements

All commits must pass the following automated checks:

#### Code Quality
- **Go**: `gofmt`, `golint`, `go vet`, `staticcheck`
- **Rust**: `rustfmt`, `clippy`
- **YAML**: `yamllint`
- **General**: `pre-commit-hooks` (trailing whitespace, large files, etc.)

#### Security Checks
- **Go**: `gosec` (security linter)
- **Rust**: `cargo audit` (dependency vulnerabilities)
- **Secrets**: `detect-secrets` (prevent credential commits)

#### Test Coverage
- **Coverage Report**: Generate and validate coverage meets minimums
- **Coverage Diff**: Ensure new code doesn't reduce overall coverage
- **Test Execution**: All tests must pass

### 4. Code Quality Standards

#### Go Code Standards
- Follow effective Go practices
- Use meaningful variable and function names
- Include godoc comments for all exported functions
- Handle all errors explicitly
- Use structured logging with consistent fields

#### Rust Code Standards
- Follow Rust idioms and best practices
- Use `#[derive]` for common traits when appropriate
- Prefer owned types over borrowed in public APIs
- Use `Result<T, E>` for fallible operations

#### Configuration Standards
- Use YAML for new configurations (JSON supported for legacy)
- Include schema version in all config files
- Document all parameters with inline comments
- Validate against vendor specifications

## Development Workflow

### 1. Feature Development
```bash
# 1. Create feature branch
git checkout -b feature/your-feature-name

# 2. Install pre-commit hooks (first time only)
pre-commit install

# 3. Develop with tests
go test ./... -v -cover
cargo test

# 4. Run pre-commit checks
pre-commit run --all-files

# 5. Commit with coverage verification
git commit -m "feat: your feature description"
```

### 2. Pull Request Requirements

Before submitting a PR, ensure:

- [ ] All tests pass: `make test`
- [ ] Coverage meets requirements: `make coverage`
- [ ] Pre-commit hooks pass: `pre-commit run --all-files`
- [ ] Documentation updated for new features
- [ ] CHANGELOG.md updated with changes
- [ ] No security vulnerabilities: `make security-check`

### 3. Code Review Standards

#### Reviewer Checklist
- [ ] Code follows project standards
- [ ] Tests adequately cover new functionality
- [ ] Documentation is clear and complete
- [ ] No performance regressions
- [ ] Security implications considered
- [ ] Configuration changes are backward compatible

## Project Structure Standards

### Directory Organization
```
phantom-science/
├── cmd/                    # Application entry points
├── internal/               # Private application code
│   ├── cli/               # Command-line interface
│   ├── genomics/          # Core FASTQ generation
│   ├── fingerprinting/    # Security features
│   └── validation/        # Quality assurance
├── pkg/                   # Public library code
│   └── types/            # Shared type definitions
├── configs/              # Configuration files
│   ├── sequencer/v1.0/   # Device-specific configs
│   └── phantom/v1.0/     # Enterprise configs
├── crypto-engine/        # Rust cryptographic code
├── docs/                 # Documentation
└── scripts/              # Build and utility scripts
```

### File Naming Conventions
- **Go files**: `snake_case.go`
- **Test files**: `*_test.go`
- **Config files**: `kebab-case.yaml`
- **Documentation**: `UPPERCASE.md` for root-level docs

## Quality Assurance

### 1. Automated Testing

#### Unit Test Requirements
```go
// Example test structure
func TestFunctionName(t *testing.T) {
    // Arrange
    input := "test input"
    expected := "expected output"
    
    // Act
    result := FunctionName(input)
    
    // Assert
    if result != expected {
        t.Errorf("Expected %s, got %s", expected, result)
    }
}
```

#### Coverage Measurement
```bash
# Generate coverage report
go test ./... -coverprofile=coverage.out
go tool cover -html=coverage.out -o coverage.html

# Check coverage thresholds
go tool cover -func=coverage.out | grep "total:" | awk '{print $3}'
```

### 2. Performance Standards

#### FASTQ Generation Performance
- **Throughput**: Minimum 1M reads/minute for volume testing
- **Memory Usage**: Maximum 4GB for standard operations
- **File Size**: Support files up to 10GB without memory issues

#### Benchmarking Requirements
```go
// Example benchmark
func BenchmarkFASTQGeneration(b *testing.B) {
    for i := 0; i < b.N; i++ {
        generateFASTQFile(testConfig)
    }
}
```

### 3. Security Standards

#### Data Protection
- All customer data must be encrypted at rest
- Steganographic fingerprints must be undetectable
- Audit logs must be tamper-evident
- No secrets in configuration files or code

#### Vulnerability Management
- Monthly dependency security scans
- SAST (Static Application Security Testing) in CI/CD
- Security code review for all changes to:
  - Fingerprinting modules
  - Encryption/decryption functions
  - Authentication/authorization
  - Data handling pipelines

## Tools and Dependencies

### Required Development Tools
- **Go**: 1.21+ with modules enabled
- **Rust**: 1.70+ with Cargo
- **pre-commit**: For automated quality checks
- **golangci-lint**: Advanced Go linting
- **yamllint**: YAML file validation

### Installation
```bash
# Install Go tools
go install honnef.co/go/tools/cmd/staticcheck@latest
go install github.com/securecodewarrior/gosec/v2/cmd/gosec@latest

# Install pre-commit
pip install pre-commit
pre-commit install

# Install Rust tools
rustup component add clippy rustfmt
cargo install cargo-audit
```

## Release Process

### Version Numbering
- Follow semantic versioning (SemVer)
- Format: `MAJOR.MINOR.PATCH`
- Current: `0.1.0` (initial release)

### Release Checklist
- [ ] All tests pass with required coverage
- [ ] Security scan shows no critical vulnerabilities
- [ ] Performance benchmarks meet standards
- [ ] Documentation is up to date
- [ ] CHANGELOG.md includes all changes
- [ ] Version number updated in `main.go`
- [ ] Git tag created and pushed

### Branch Strategy
- **main**: Production-ready code
- **develop**: Integration branch for features
- **feature/***: Individual feature development
- **hotfix/***: Critical bug fixes
- **release/***: Release preparation

## Troubleshooting

### Common Issues

#### Coverage Too Low
```bash
# Identify uncovered code
go test ./... -coverprofile=coverage.out
go tool cover -func=coverage.out | grep -v "100.0%"
```

#### Pre-commit Failures
```bash
# Run individual hooks
pre-commit run gofmt --all-files
pre-commit run golint --all-files
pre-commit run yamllint --all-files
```

#### Build Issues
```bash
# Clean and rebuild
go clean -cache
go mod tidy
go build ./...
```

## Contact

For questions about development standards or processes:
- **Technical Lead**: Create GitHub issue with `development` label
- **Security Issues**: Create GitHub issue with `security` label
- **Documentation**: Create GitHub issue with `documentation` label