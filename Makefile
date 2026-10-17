# Phantom Science - Development Makefile
# Provides standardized commands for development, testing, and quality assurance

.PHONY: help build test coverage lint security clean install-tools pre-commit fmt vet stegano

# The stegano library (built by `make stegano`) for CGO builds of internal/fingerprinting.
STEGANO_LIB_DIR := $(CURDIR)/.stegano/target/release
export CGO_LDFLAGS ?= -L$(STEGANO_LIB_DIR) -Wl,-rpath,$(STEGANO_LIB_DIR)

# Default target
help: ## Show this help message
	@echo "Phantom Science - Development Commands"
	@echo "====================================="
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

# Build targets
build: ## Build the application
	@echo "🔨 Building Phantom Science..."
	@mkdir -p bin
	go build -o bin/phantom-science ./...
	@echo "✅ Build complete: bin/phantom-science"

stegano: ## Build the stegano library version pinned in go.mod (needed for CGO builds)
	@scripts/build-stegano.sh

build-all: ## Build Go and Rust components
	@echo "🔨 Building all components..."
	@mkdir -p bin
	go build -o bin/phantom-science ./...
	cd crypto-engine && cargo build --release
	@echo "✅ Build complete: Go and Rust components"

# Testing targets
test: ## Run all tests
	@echo "🧪 Running tests..."
	go test ./... -v

test-short: ## Run short tests only
	@echo "🧪 Running short tests..."
	go test ./... -short

test-integration: ## Run integration tests (no external tool dependencies)
	@echo "🧪 Running integration tests..."
	go test -tags '!release' ./tests/integration/... -v -timeout 120s

test-e2e: ## Run full end-to-end pipeline tests (slow; skipped with -short)
	@echo "🧪 Running end-to-end tests..."
	go test -tags '!release' ./tests/e2e/... -v -timeout 300s

coverage: ## Generate test coverage report (enforces 85%+ overall, 80%+ individual)
	@echo "📊 Generating coverage report..."
	@scripts/check-coverage.sh 85 80

coverage-html: ## Generate HTML coverage report
	@echo "📊 Generating HTML coverage report..."
	go test ./... -coverprofile=coverage.out -covermode=atomic
	go tool cover -html=coverage.out -o coverage.html
	@echo "📊 Coverage report: coverage.html"

# Code quality targets
lint: ## Run all linters
	@echo "🔍 Running linters..."
	@echo "Go linting..."
	@if command -v golangci-lint >/dev/null 2>&1; then golangci-lint run ./...; else echo "golangci-lint not found, install with: make install-tools"; fi
	@echo "YAML linting..."
	@if command -v yamllint >/dev/null 2>&1; then yamllint configs/; else echo "yamllint not found, install with: make install-tools"; fi

fmt: ## Format all code
	@echo "🎨 Formatting code..."
	go fmt ./...
	@if [ -d crypto-engine ]; then cd crypto-engine && cargo fmt; fi
	@echo "✅ Code formatted"

vet: ## Run go vet
	@echo "🔍 Running go vet..."
	go vet ./...

staticcheck: ## Run staticcheck
	@echo "🔍 Running staticcheck..."
	@if command -v staticcheck >/dev/null 2>&1; then staticcheck ./...; else echo "staticcheck not found, install with: make install-tools"; fi

# Security targets  
security: ## Run security checks
	@echo "🔒 Running security checks..."
	@echo "Go security scan..."
	@if command -v gosec >/dev/null 2>&1; then gosec -quiet ./...; else echo "gosec not found, install with: make install-tools"; fi
	@if [ -d crypto-engine ]; then echo "Rust security audit..." && cd crypto-engine && cargo audit; fi
	@echo "✅ Security checks complete"

# Configuration validation
validate-configs: ## Validate all configuration files
	@echo "⚙️  Validating configurations..."
	@scripts/validate-configs.sh

validate-sequencer-specs: ## Validate sequencer specifications
	@echo "🧬 Validating sequencer specifications..."
	@scripts/validate-sequencer-specs.sh

validate-changelog: ## Validate changelog format (Keep a Changelog standard)
	@echo "📝 Validating changelog format..."
	@scripts/validate-changelog.sh

# Pre-commit and development workflow
pre-commit: ## Run pre-commit hooks on all files
	@echo "🔧 Running pre-commit hooks..."
	@if command -v pre-commit >/dev/null 2>&1; then pre-commit run --all-files; else echo "pre-commit not found, install with: make install-tools"; fi

install-tools: ## Install development tools
	@echo "🛠️  Installing development tools..."
	@echo "Installing Go tools..."
	go install honnef.co/go/tools/cmd/staticcheck@latest
	go install github.com/securecodewarrior/gosec/v2/cmd/gosec@latest
	@if command -v golangci-lint >/dev/null 2>&1; then echo "golangci-lint already installed"; else \
		echo "Installing golangci-lint..."; \
		curl -sSfL https://raw.githubusercontent.com/golangci/golangci-lint/master/install.sh | sh -s -- -b $$(go env GOPATH)/bin latest; \
	fi
	@echo "Installing Python tools..."
	@if command -v pip3 >/dev/null 2>&1; then pip3 install pre-commit yamllint detect-secrets; else echo "pip3 not found, install Python first"; fi
	@if [ -d crypto-engine ]; then \
		echo "Installing Rust tools..."; \
		rustup component add clippy rustfmt; \
		cargo install cargo-audit; \
	fi
	@echo "✅ Development tools installed"

install-pre-commit: ## Install pre-commit hooks
	@echo "🔧 Installing pre-commit hooks..."
	@if command -v pre-commit >/dev/null 2>&1; then pre-commit install; else echo "pre-commit not found, run: make install-tools"; fi
	@echo "✅ Pre-commit hooks installed"

# Performance testing
benchmark: ## Run performance benchmarks
	@echo "⚡ Running benchmarks..."
	go test ./... -bench=. -benchmem
	@if [ -d crypto-engine ]; then cd crypto-engine && cargo bench; fi

# Cleanup targets
clean: ## Clean build artifacts and temporary files
	@echo "🧹 Cleaning up..."
	rm -f bin/phantom-science
	rm -f coverage.out coverage.html
	@if [ -d crypto-engine ]; then rm -rf crypto-engine/target/release/; fi
	rm -rf *_output/
	go clean -cache -testcache
	@echo "✅ Cleanup complete"

clean-all: clean ## Clean everything including downloaded dependencies
	@echo "🧹 Deep cleaning..."
	go clean -modcache
	@if [ -d crypto-engine ]; then cd crypto-engine && cargo clean; fi
	@echo "✅ Deep cleanup complete"

# Version and release
version: ## Show current version
	@echo "Phantom Science v$$(cat VERSION)"

check-version: ## Validate semantic versioning compliance
	@echo "🔍 Checking version compliance..."
	@scripts/check-version.sh

sync-version: ## Sync version across all files
	@echo "🔄 Syncing version across files..."
	@VERSION=$$(cat VERSION) && \
	sed -i.bak "s/Version.*=.*/Version   = \"$$VERSION\"/" cmd/root.go && \
	rm -f cmd/root.go.bak
	@echo "✅ Version synced to $$(cat VERSION)"

ci-status: ## Check CI/CD pipeline status
	@echo "🔍 Checking CI/CD status..."
	@scripts/ci-status-monitor.sh

ci-watch: ## Monitor CI/CD status continuously
	@echo "👀 Starting CI/CD monitoring..."
	@scripts/ci-status-monitor.sh --watch

# Development workflow shortcuts
dev: build test lint ## Full development workflow: build, test, lint
	@echo "✅ Development workflow complete"

ci: test coverage lint security validate-configs validate-changelog ## CI pipeline: test, coverage, lint, security, validation
	@echo "✅ CI pipeline complete"

quick: build test-short ## Quick development check: build and short tests
	@echo "✅ Quick check complete"

# Legacy Docker targets (preserved)
docker-build: ## Build Docker image
	@echo "🐳 Building Docker image..."
	@if [ -f deploy/docker/Dockerfile ]; then docker build -f deploy/docker/Dockerfile -t phantom-science:latest .; else echo "Dockerfile not found"; fi

# Configuration management
config-examples: ## Show available configuration examples
	@echo "📝 Available configuration examples:"
	@echo "Sequencer configs:"
	@find configs/sequencer -name "*.yaml" -o -name "*.json" | head -10
	@echo "Phantom configs:"
	@find configs/phantom -name "*.yaml" -o -name "*.json" | head -5

# Development setup
setup: install-tools install-pre-commit ## Complete development setup
	@echo "🚀 Development environment setup complete!"
	@echo ""
	@echo "Next steps:"
	@echo "  make build     # Build the application"
	@echo "  make test      # Run tests"
	@echo "  make dev       # Full development workflow"