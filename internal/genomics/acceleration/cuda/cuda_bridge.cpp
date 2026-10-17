#include "cuda_bridge.h"
#include <cstring>
#include <cstdlib>

// Conditionally include CUDA headers.  When the CUDA Toolkit is not installed
// (e.g., on macOS developer machines or CI systems) we compile stub
// implementations that return safe "not available" values.  This lets the
// package build without the CUDA SDK while still exporting the symbols that
// the CGO Go wrapper expects.
//
// To build with real CUDA support, install the CUDA Toolkit and build with
// the 'cuda' build tag:
//   CGO_ENABLED=1 go build -tags cuda ./...
#ifdef __has_include
#  if __has_include(<cuda_runtime.h>)
#    define CUDA_AVAILABLE 1
#  endif
#endif

#ifdef CUDA_AVAILABLE
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>
#include <memory>
#include <algorithm>

// External kernel declarations
extern "C" {
    void extract_kmers_kernel(const char* sequences, const void* seq_info,
                             unsigned int* kmer_counts, unsigned int k,
                             unsigned int num_sequences, unsigned int max_kmers);

    void extract_kmers_small_k_kernel(const char* sequences, const void* seq_info,
                                     unsigned int* kmer_counts, unsigned int k,
                                     unsigned int num_sequences, unsigned int max_kmers);

    void extract_kmers_cooperative_kernel(const char* sequences, const void* seq_info,
                                         unsigned int* kmer_counts, unsigned int k,
                                         unsigned int num_sequences, unsigned int max_kmers);

    void clear_kmer_counts_kernel(unsigned int* kmer_counts, unsigned int max_kmers);
}

// CUDA error checking macro
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(error) << std::endl; \
            return -1; \
        } \
    } while(0)

#define CUDA_CHECK_BOOL(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(error) << std::endl; \
            return false; \
        } \
    } while(0)

// Structure for sequence metadata (must match CUDA kernel)
struct SequenceInfo {
    unsigned int offset;
    unsigned int length;
    unsigned int seq_id;
    unsigned int padding;
};

// CUDA device context structure
struct CUDADevice {
    int device_id;
    cudaDeviceProp properties;
    bool is_initialized;
    void* stream;  // cudaStream_t

    // Device capabilities
    size_t max_threads_per_block;
    size_t max_blocks_per_grid;
    size_t max_shared_memory;
    size_t global_memory_size;

    // Performance characteristics
    int compute_capability_major;
    int compute_capability_minor;
    int multiprocessor_count;
    int warp_size;
};

// Check if CUDA is available on this system
extern "C" bool isCudaAvailable() {
    int device_count = 0;
    cudaError_t error = cudaGetDeviceCount(&device_count);

    if (error != cudaSuccess) {
        return false;
    }

    return device_count > 0;
}

// Initialize CUDA device and get capabilities
extern "C" CUDADevice* initializeCudaDevice() {
    // Check CUDA availability
    if (!isCudaAvailable()) {
        std::cerr << "CUDA Error: No CUDA devices available" << std::endl;
        return nullptr;
    }

    // Allocate device structure
    auto* device = new CUDADevice();
    memset(device, 0, sizeof(CUDADevice));

    // Get device with highest compute capability
    int device_id = 0;
    int max_compute_capability = 0;
    int device_count = 0;

    CUDA_CHECK_BOOL(cudaGetDeviceCount(&device_count));

    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp props;
        CUDA_CHECK_BOOL(cudaGetDeviceProperties(&props, i));

        int compute_capability = props.major * 10 + props.minor;
        if (compute_capability > max_compute_capability) {
            max_compute_capability = compute_capability;
            device_id = i;
        }
    }

    // Set device
    CUDA_CHECK_BOOL(cudaSetDevice(device_id));

    // Get device properties
    CUDA_CHECK_BOOL(cudaGetDeviceProperties(&device->properties, device_id));

    device->device_id = device_id;
    device->max_threads_per_block = device->properties.maxThreadsPerBlock;
    device->max_blocks_per_grid = device->properties.maxGridSize[0];
    device->max_shared_memory = device->properties.sharedMemPerBlock;
    device->global_memory_size = device->properties.totalGlobalMem;
    device->compute_capability_major = device->properties.major;
    device->compute_capability_minor = device->properties.minor;
    device->multiprocessor_count = device->properties.multiProcessorCount;
    device->warp_size = device->properties.warpSize;

    // Create CUDA stream for asynchronous operations
    cudaStream_t stream;
    CUDA_CHECK_BOOL(cudaStreamCreate(&stream));
    device->stream = static_cast<void*>(stream);

    device->is_initialized = true;

    std::cout << "CUDA GPU: Initialized device " << device_id
              << " (" << device->properties.name << ")" << std::endl;
    std::cout << "CUDA GPU: Compute capability "
              << device->compute_capability_major << "."
              << device->compute_capability_minor << std::endl;
    std::cout << "CUDA GPU: " << device->multiprocessor_count
              << " SMs, " << device->max_threads_per_block
              << " threads/block" << std::endl;
    std::cout << "CUDA GPU: " << (device->global_memory_size / (1024*1024))
              << " MB global memory" << std::endl;

    return device;
}

// Process k-mers on GPU using CUDA
extern "C" int processKmersOnCuda(CUDADevice* device, const char* sequences,
                                  int num_sequences, int* sequence_lengths,
                                  int k, int* results, int max_results) {
    if (!device || !device->is_initialized) {
        std::cerr << "CUDA Error: Device not initialized" << std::endl;
        return -1;
    }

    if (!sequences || !sequence_lengths || !results || k <= 0 || k > 16) {
        std::cerr << "CUDA Error: Invalid parameters" << std::endl;
        return -2;
    }

    cudaStream_t stream = static_cast<cudaStream_t>(device->stream);

    // Calculate total sequence data size
    size_t total_sequence_size = 0;
    for (int i = 0; i < num_sequences; i++) {
        total_sequence_size += sequence_lengths[i] + 1;  // +1 for null terminator
    }

    std::cout << "CUDA GPU: Processing " << num_sequences
              << " sequences (" << total_sequence_size << " bytes) with k=" << k << std::endl;

    // Check memory requirements
    size_t required_memory = total_sequence_size +
                           num_sequences * sizeof(SequenceInfo) +
                           max_results * sizeof(unsigned int) * 2;

    if (required_memory > device->global_memory_size * 0.8) {
        std::cerr << "CUDA Error: Dataset too large for GPU memory ("
                  << (required_memory / (1024*1024)) << " MB required, "
                  << (device->global_memory_size / (1024*1024)) << " MB available)" << std::endl;
        return -3;
    }

    // Create sequence info structures
    auto sequence_infos = std::make_unique<SequenceInfo[]>(num_sequences);
    unsigned int offset = 0;
    for (int i = 0; i < num_sequences; i++) {
        sequence_infos[i].offset = offset;
        sequence_infos[i].length = sequence_lengths[i];
        sequence_infos[i].seq_id = i;
        sequence_infos[i].padding = 0;
        offset += sequence_lengths[i] + 1;
    }

    // Allocate GPU memory
    char* d_sequences = nullptr;
    SequenceInfo* d_seq_info = nullptr;
    unsigned int* d_kmer_counts = nullptr;

    CUDA_CHECK(cudaMalloc(&d_sequences, total_sequence_size));
    CUDA_CHECK(cudaMalloc(&d_seq_info, num_sequences * sizeof(SequenceInfo)));
    CUDA_CHECK(cudaMalloc(&d_kmer_counts, max_results * sizeof(unsigned int)));

    // Copy data to GPU
    CUDA_CHECK(cudaMemcpyAsync(d_sequences, sequences, total_sequence_size,
                              cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_seq_info, sequence_infos.get(),
                              num_sequences * sizeof(SequenceInfo),
                              cudaMemcpyHostToDevice, stream));

    // Clear k-mer counts buffer
    unsigned int clear_threads = std::min(1024U, static_cast<unsigned int>(device->max_threads_per_block));
    unsigned int clear_blocks = (max_results + clear_threads - 1) / clear_threads;

    clear_kmer_counts_kernel<<<clear_blocks, clear_threads, 0, stream>>>(
        d_kmer_counts, max_results);

    CUDA_CHECK(cudaGetLastError());

    // Choose optimal kernel configuration based on dataset characteristics
    unsigned int threads_per_block;
    unsigned int blocks_per_grid;
    size_t shared_memory_size = 0;

    if (k <= 8 && num_sequences <= 10000) {
        threads_per_block = std::min(256U, static_cast<unsigned int>(device->max_threads_per_block));
        blocks_per_grid = (num_sequences + threads_per_block - 1) / threads_per_block;
        shared_memory_size = (1U << (2 * k)) * sizeof(unsigned int);

        if (shared_memory_size > device->max_shared_memory) {
            shared_memory_size = 0;
        }
    }

    if (shared_memory_size == 0) {
        if (total_sequence_size > 100000000 || num_sequences > 50000) {
            threads_per_block = std::min(512U, static_cast<unsigned int>(device->max_threads_per_block));
            blocks_per_grid = std::min(
                static_cast<unsigned int>(device->multiprocessor_count * 4),
                (num_sequences + threads_per_block - 1) / threads_per_block
            );
        } else {
            threads_per_block = std::min(256U, static_cast<unsigned int>(device->max_threads_per_block));
            blocks_per_grid = (num_sequences + threads_per_block - 1) / threads_per_block;
        }
    }

    // Launch appropriate kernel
    if (shared_memory_size > 0) {
        extract_kmers_small_k_kernel<<<blocks_per_grid, threads_per_block,
                                      shared_memory_size, stream>>>(
            d_sequences, d_seq_info, d_kmer_counts, k, num_sequences, max_results);
    } else if (total_sequence_size > 100000000 || num_sequences > 50000) {
        extract_kmers_cooperative_kernel<<<blocks_per_grid, threads_per_block, 0, stream>>>(
            d_sequences, d_seq_info, d_kmer_counts, k, num_sequences, max_results);
    } else {
        extract_kmers_kernel<<<blocks_per_grid, threads_per_block, 0, stream>>>(
            d_sequences, d_seq_info, d_kmer_counts, k, num_sequences, max_results);
    }

    CUDA_CHECK(cudaGetLastError());

    // Copy results back to CPU
    auto results_buffer = std::make_unique<unsigned int[]>(max_results);
    CUDA_CHECK(cudaMemcpyAsync(results_buffer.get(), d_kmer_counts,
                              max_results * sizeof(unsigned int),
                              cudaMemcpyDeviceToHost, stream));

    // Wait for completion
    CUDA_CHECK(cudaStreamSynchronize(stream));

    // Convert results and calculate total
    int total_kmers = 0;
    for (int i = 0; i < max_results; i++) {
        results[i] = static_cast<int>(results_buffer[i]);
        total_kmers += results[i];
    }

    // Clean up GPU memory
    cudaFree(d_sequences);
    cudaFree(d_seq_info);
    cudaFree(d_kmer_counts);

    std::cout << "CUDA GPU: K-mer extraction completed, "
              << total_kmers << " total k-mers processed" << std::endl;

    return total_kmers;
}

// Get CUDA device name
extern "C" const char* getCudaDeviceName(CUDADevice* device) {
    if (!device || !device->is_initialized) {
        return "Unknown";
    }
    return device->properties.name;
}

// Get maximum threads per block
extern "C" int getCudaMaxThreadsPerBlock(CUDADevice* device) {
    if (!device || !device->is_initialized) {
        return 0;
    }
    return static_cast<int>(device->max_threads_per_block);
}

// Get total global memory
extern "C" unsigned long getCudaGlobalMemorySize(CUDADevice* device) {
    if (!device || !device->is_initialized) {
        return 0;
    }
    return static_cast<unsigned long>(device->global_memory_size);
}

// Get compute capability
extern "C" int getCudaComputeCapability(CUDADevice* device) {
    if (!device || !device->is_initialized) {
        return 0;
    }
    return device->compute_capability_major * 10 + device->compute_capability_minor;
}

// Get multiprocessor count
extern "C" int getCudaMultiprocessorCount(CUDADevice* device) {
    if (!device || !device->is_initialized) {
        return 0;
    }
    return device->multiprocessor_count;
}

// Release CUDA device and resources
extern "C" void releaseCudaDevice(CUDADevice* device) {
    if (!device) {
        return;
    }

    std::cout << "CUDA GPU: Releasing device resources" << std::endl;

    if (device->is_initialized) {
        if (device->stream) {
            cudaStream_t stream = static_cast<cudaStream_t>(device->stream);
            cudaStreamSynchronize(stream);
            cudaStreamDestroy(stream);
        }

        cudaDeviceReset();
        device->is_initialized = false;
    }

    delete device;
}

#else // CUDA_AVAILABLE not defined — provide CPU stubs

// Opaque stub for CUDADevice when the CUDA SDK is not present.
// All functions return "not available" values so that the CGO wrapper can
// detect the absence of CUDA at runtime via isCudaAvailable() == false.
struct CUDADevice {
    int device_id;
    bool is_initialized;
};

extern "C" bool isCudaAvailable() { return false; }

extern "C" CUDADevice* initializeCudaDevice() { return nullptr; }

extern "C" int processKmersOnCuda(CUDADevice* /*device*/, const char* /*sequences*/,
                                  int /*num_sequences*/, int* /*sequence_lengths*/,
                                  int /*k*/, int* /*results*/, int /*max_results*/) {
    return -1;
}

extern "C" const char* getCudaDeviceName(CUDADevice* /*device*/) { return "CUDA unavailable"; }
extern "C" int getCudaMaxThreadsPerBlock(CUDADevice* /*device*/) { return 0; }
extern "C" unsigned long getCudaGlobalMemorySize(CUDADevice* /*device*/) { return 0; }
extern "C" int getCudaComputeCapability(CUDADevice* /*device*/) { return 0; }
extern "C" int getCudaMultiprocessorCount(CUDADevice* /*device*/) { return 0; }
extern "C" void releaseCudaDevice(CUDADevice* device) { free(device); }

#endif // CUDA_AVAILABLE
