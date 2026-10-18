# Makefile for Metal GPU acceleration support
# Handles Metal shader compilation and integration

# Variables
METAL_SRC_DIR = internal/genomics/acceleration/metal
METAL_SHADER = $(METAL_SRC_DIR)/kmer_kernels.metal
METAL_LIB = $(METAL_SRC_DIR)/kmer_kernels.metallib
METAL_AIR = $(METAL_SRC_DIR)/kmer_kernels.air

# Check if we're on macOS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    HAS_METAL = true
    XCRUN = $(shell which xcrun)
else
    HAS_METAL = false
endif

.PHONY: metal-check metal-compile metal-clean metal-install gpu-benchmark help-metal

# Default target
metal: metal-check metal-compile

# Check Metal environment
metal-check:
ifndef HAS_METAL
	@echo "❌ Metal GPU acceleration is only available on macOS"
	@exit 1
endif
ifndef XCRUN
	@echo "❌ xcrun not found - please install Xcode command line tools"
	@echo "   Run: xcode-select --install"
	@exit 1
endif
	@echo "✅ Metal development environment available"
	@xcrun metal --version

# Compile Metal shaders
metal-compile: metal-check
	@echo "🔨 Compiling Metal shaders..."
	@echo "   Source: $(METAL_SHADER)"
	@echo "   Output: $(METAL_LIB)"
	@mkdir -p $(dir $(METAL_LIB))
	@xcrun -sdk macosx metal -c $(METAL_SHADER) -o $(METAL_AIR)
	@xcrun -sdk macosx metallib $(METAL_AIR) -o $(METAL_LIB)
	@rm -f $(METAL_AIR)
	@echo "✅ Metal shaders compiled successfully"
	@ls -lh $(METAL_LIB)

# Clean Metal build artifacts
metal-clean:
	@echo "🧹 Cleaning Metal build artifacts..."
	@rm -f $(METAL_LIB) $(METAL_AIR)
	@echo "✅ Metal artifacts cleaned"

# Install Metal library to application bundle (for distribution)
metal-install: metal-compile
	@echo "📦 Installing Metal library for distribution..."
	@if [ -d "phantom-science.app" ]; then \
		mkdir -p phantom-science.app/Contents/Resources; \
		cp $(METAL_LIB) phantom-science.app/Contents/Resources/; \
		echo "✅ Metal library installed to application bundle"; \
	else \
		echo "⚠️  No application bundle found - Metal library will be loaded from source"; \
	fi

# Build GPU benchmark tool
gpu-benchmark: metal-compile
	@echo "🚀 Building GPU benchmark tool..."
	@go build -o bin/benchmark_gpu ./cmd/experiments/benchmark_gpu
	@echo "✅ GPU benchmark tool built: bin/benchmark_gpu"

# Run comprehensive GPU benchmarks
gpu-bench-run: gpu-benchmark
	@echo "🧪 Running comprehensive GPU benchmarks..."
	@./bin/benchmark_gpu -k 6 -sequences 10000 -length 1000 -verbose
	@echo ""
	@./bin/benchmark_gpu -k 8 -sequences 5000 -length 2000 -verbose
	@echo ""
	@./bin/benchmark_gpu -k 10 -sequences 2000 -length 5000 -verbose

# Quick GPU performance test
gpu-test: gpu-benchmark
	@echo "⚡ Quick GPU performance test..."
	@./bin/benchmark_gpu -k 6 -sequences 1000 -length 1000

# Build with Metal support
build-metal: metal-compile
	@echo "🔨 Building phantom-science with Metal GPU support..."
	@CGO_ENABLED=1 go build -tags metal -ldflags="-framework Metal -framework Foundation" -o phantom-science .
	@echo "✅ Built with Metal GPU acceleration support"

# Test Metal integration
test-metal: metal-compile
	@echo "🧪 Testing Metal GPU integration..."
	@CGO_ENABLED=1 go test -tags metal ./internal/genomics/acceleration/metal/... -v
	@echo "✅ Metal integration tests completed"

# Profile GPU performance
profile-gpu: gpu-benchmark
	@echo "📊 Profiling GPU performance..."
	@echo "Starting Xcode Instruments for GPU profiling..."
	@echo "Run: ./bin/benchmark_gpu -k 8 -sequences 5000 -length 2000"
	@echo "Then attach Instruments to profile GPU usage"

# Help for Metal-specific targets
help-metal:
	@echo "🍎 Metal GPU Acceleration Targets:"
	@echo ""
	@echo "  metal-check     - Check Metal development environment"
	@echo "  metal-compile   - Compile Metal compute shaders"
	@echo "  metal-clean     - Clean compiled Metal artifacts"
	@echo "  metal-install   - Install Metal library to app bundle"
	@echo ""
	@echo "  gpu-benchmark   - Build GPU benchmark tool"
	@echo "  gpu-bench-run   - Run comprehensive GPU benchmarks"
	@echo "  gpu-test        - Quick GPU performance test"
	@echo ""
	@echo "  build-metal     - Build phantom-science with Metal support"
	@echo "  test-metal      - Run Metal integration tests"
	@echo "  profile-gpu     - Instructions for GPU profiling"
	@echo ""
	@echo "Metal GPU acceleration requires:"
	@echo "  - macOS with Metal support (M1/M2/M3 or Intel with discrete GPU)"
	@echo "  - Xcode command line tools (xcrun metal compiler)"
	@echo "  - CGO enabled for Metal framework integration"

# Performance targets for different k-mer sizes
perf-k4: gpu-benchmark
	@./bin/benchmark_gpu -k 4 -sequences 20000 -length 500

perf-k6: gpu-benchmark
	@./bin/benchmark_gpu -k 6 -sequences 10000 -length 1000

perf-k8: gpu-benchmark
	@./bin/benchmark_gpu -k 8 -sequences 5000 -length 2000

perf-k10: gpu-benchmark
	@./bin/benchmark_gpu -k 10 -sequences 2000 -length 5000

perf-k12: gpu-benchmark
	@./bin/benchmark_gpu -k 12 -sequences 1000 -length 10000

# Memory stress test
perf-memory: gpu-benchmark
	@echo "🧠 GPU memory stress test..."
	@./bin/benchmark_gpu -k 8 -sequences 50000 -length 2000 -verbose

# All performance tests
perf-all: perf-k4 perf-k6 perf-k8 perf-k10 perf-k12 perf-memory