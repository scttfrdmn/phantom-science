#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cub/cub.cuh>

// CUDA kernels for parallel k-mer extraction and counting with atomic
// increments into a dense count table. No benchmark in the repository
// measures their performance.

// Base encoding lookup table: A=0, T=1, G=2, C=3, invalid=4
__constant__ unsigned char BASE_ENCODING[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0-15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 16-31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 32-47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 48-63
    4, 0, 4, 3, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 64-79  (A=65->0, C=67->3, G=71->2)
    4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 80-95  (T=84->1)
    4, 0, 4, 3, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 96-111 (a=97->0, c=99->3, g=103->2)
    4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 112-127 (t=116->1)
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 128+
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

// Structure for sequence metadata
struct SequenceInfo {
    unsigned int offset;      // Offset in sequence buffer
    unsigned int length;      // Length of this sequence
    unsigned int seq_id;      // Sequence identifier
    unsigned int padding;     // Alignment padding
};

// Standard k-mer extraction kernel for general-purpose use
__global__ void extract_kmers_kernel(
    const char* __restrict__ sequences,
    const SequenceInfo* __restrict__ seq_info,
    unsigned int* __restrict__ kmer_counts,
    const unsigned int k,
    const unsigned int num_sequences,
    const unsigned int max_kmers
) {
    unsigned int seq_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (seq_idx >= num_sequences) {
        return;
    }
    
    SequenceInfo info = seq_info[seq_idx];
    
    // Skip sequences shorter than k
    if (info.length < k) {
        return;
    }
    
    const char* seq = sequences + info.offset;
    unsigned int kmer_mask = (1u << (2 * k)) - 1;
    unsigned int kmer_hash = 0;
    unsigned int valid_bases = 0;
    
    // Initialize sliding window
    for (unsigned int i = 0; i < k && i < info.length; i++) {
        unsigned char base = BASE_ENCODING[(unsigned char)seq[i]];
        if (base < 4) {  // Valid base
            kmer_hash = (kmer_hash << 2) | base;
            valid_bases++;
        } else {  // Invalid base - reset
            kmer_hash = 0;
            valid_bases = 0;
        }
    }
    
    // Count first k-mer if valid
    if (valid_bases >= k && kmer_hash < max_kmers) {
        atomicAdd(&kmer_counts[kmer_hash], 1);
    }
    
    // Slide through remaining sequence
    for (unsigned int i = k; i < info.length; i++) {
        unsigned char base = BASE_ENCODING[(unsigned char)seq[i]];
        
        if (base < 4) {  // Valid base
            kmer_hash = ((kmer_hash << 2) | base) & kmer_mask;
            valid_bases = min(valid_bases + 1, k);
            
            if (valid_bases >= k && kmer_hash < max_kmers) {
                atomicAdd(&kmer_counts[kmer_hash], 1);
            }
        } else {  // Invalid base - reset window
            kmer_hash = 0;
            valid_bases = 0;
        }
    }
}

// Optimized kernel for small k-mers (k <= 8) using shared memory
__global__ void extract_kmers_small_k_kernel(
    const char* __restrict__ sequences,
    const SequenceInfo* __restrict__ seq_info,
    unsigned int* __restrict__ kmer_counts,
    const unsigned int k,
    const unsigned int num_sequences,
    const unsigned int max_kmers
) {
    // Shared memory for local accumulation
    extern __shared__ unsigned int shared_counts[];
    
    unsigned int tid = threadIdx.x;
    unsigned int block_size = blockDim.x;
    unsigned int seq_idx = blockIdx.x * block_size + tid;
    
    // Initialize shared memory
    for (unsigned int i = tid; i < max_kmers; i += block_size) {
        shared_counts[i] = 0;
    }
    __syncthreads();
    
    // Process sequence if within bounds
    if (seq_idx < num_sequences) {
        SequenceInfo info = seq_info[seq_idx];
        
        if (info.length >= k) {
            const char* seq = sequences + info.offset;
            unsigned int kmer_mask = (1u << (2 * k)) - 1;
            unsigned int kmer_hash = 0;
            unsigned int valid_bases = 0;
            
            // Initialize sliding window
            for (unsigned int i = 0; i < k && i < info.length; i++) {
                unsigned char base = BASE_ENCODING[(unsigned char)seq[i]];
                if (base < 4) {
                    kmer_hash = (kmer_hash << 2) | base;
                    valid_bases++;
                } else {
                    kmer_hash = 0;
                    valid_bases = 0;
                }
            }
            
            // Count first k-mer
            if (valid_bases >= k && kmer_hash < max_kmers) {
                atomicAdd(&shared_counts[kmer_hash], 1);
            }
            
            // Slide through sequence
            for (unsigned int i = k; i < info.length; i++) {
                unsigned char base = BASE_ENCODING[(unsigned char)seq[i]];
                
                if (base < 4) {
                    kmer_hash = ((kmer_hash << 2) | base) & kmer_mask;
                    valid_bases = min(valid_bases + 1, k);
                    
                    if (valid_bases >= k && kmer_hash < max_kmers) {
                        atomicAdd(&shared_counts[kmer_hash], 1);
                    }
                } else {
                    kmer_hash = 0;
                    valid_bases = 0;
                }
            }
        }
    }
    
    __syncthreads();
    
    // Write shared memory results to global memory
    for (unsigned int i = tid; i < max_kmers; i += block_size) {
        if (shared_counts[i] > 0) {
            atomicAdd(&kmer_counts[i], shared_counts[i]);
        }
    }
}

// High-throughput kernel for large datasets using cooperative groups
__global__ void extract_kmers_cooperative_kernel(
    const char* __restrict__ sequences,
    const SequenceInfo* __restrict__ seq_info,
    unsigned int* __restrict__ kmer_counts,
    const unsigned int k,
    const unsigned int num_sequences,
    const unsigned int max_kmers
) {
    // Use cooperative groups for better load balancing
    unsigned int global_tid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int grid_size = gridDim.x * blockDim.x;
    
    unsigned int kmer_mask = (1u << (2 * k)) - 1;
    
    // Process multiple sequences per thread for better utilization
    for (unsigned int seq_idx = global_tid; seq_idx < num_sequences; seq_idx += grid_size) {
        SequenceInfo info = seq_info[seq_idx];
        
        if (info.length < k) {
            continue;
        }
        
        const char* seq = sequences + info.offset;
        unsigned int kmer_hash = 0;
        unsigned int valid_bases = 0;
        
        // Process sequence with optimized loop unrolling
        for (unsigned int i = 0; i < info.length; i++) {
            unsigned char base = BASE_ENCODING[(unsigned char)seq[i]];
            
            if (base < 4) {  // Valid base
                kmer_hash = ((kmer_hash << 2) | base) & kmer_mask;
                valid_bases = min(valid_bases + 1, k);
                
                // Count k-mer if we have k consecutive valid bases
                if (valid_bases >= k && kmer_hash < max_kmers) {
                    atomicAdd(&kmer_counts[kmer_hash], 1);
                }
            } else {  // Invalid base - reset
                kmer_hash = 0;
                valid_bases = 0;
            }
        }
    }
}

// Memory-efficient kernel for very large datasets using streaming
__global__ void extract_kmers_streaming_kernel(
    const char* __restrict__ sequences,
    const unsigned int* __restrict__ sequence_offsets,
    const unsigned int* __restrict__ sequence_lengths,
    unsigned int* __restrict__ kmer_counts,
    const unsigned int k,
    const unsigned int total_bases,
    const unsigned int max_kmers,
    const unsigned int batch_offset
) {
    unsigned int global_tid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int grid_size = gridDim.x * blockDim.x;
    
    unsigned int kmer_mask = (1u << (2 * k)) - 1;
    
    // Process bases in parallel across all sequences
    for (unsigned int base_idx = global_tid + batch_offset; 
         base_idx + k <= total_bases; 
         base_idx += grid_size) {
        
        unsigned int kmer_hash = 0;
        bool valid_kmer = true;
        
        // Build k-mer from k consecutive bases
        #pragma unroll
        for (unsigned int i = 0; i < 16 && i < k; i++) {  // Unroll for common k values
            if (i < k) {
                unsigned char base = BASE_ENCODING[(unsigned char)sequences[base_idx + i]];
                if (base >= 4) {
                    valid_kmer = false;
                    break;
                }
                kmer_hash = (kmer_hash << 2) | base;
            }
        }
        
        // Handle k > 16 case
        for (unsigned int i = 16; i < k; i++) {
            unsigned char base = BASE_ENCODING[(unsigned char)sequences[base_idx + i]];
            if (base >= 4) {
                valid_kmer = false;
                break;
            }
            kmer_hash = (kmer_hash << 2) | base;
        }
        
        // Count valid k-mer
        if (valid_kmer && kmer_hash < max_kmers) {
            atomicAdd(&kmer_counts[kmer_hash], 1);
        }
    }
}

// Optimized kernel using CUB for reduction operations
__global__ void extract_kmers_cub_kernel(
    const char* __restrict__ sequences,
    const SequenceInfo* __restrict__ seq_info,
    unsigned int* __restrict__ kmer_counts,
    const unsigned int k,
    const unsigned int num_sequences,
    const unsigned int max_kmers
) {
    // Use CUB block primitives for efficient reductions
    typedef cub::BlockReduce<unsigned int, 256> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    
    unsigned int tid = threadIdx.x;
    unsigned int seq_idx = blockIdx.x * blockDim.x + tid;
    
    unsigned int local_kmer_count = 0;
    
    if (seq_idx < num_sequences) {
        SequenceInfo info = seq_info[seq_idx];
        
        if (info.length >= k) {
            const char* seq = sequences + info.offset;
            unsigned int kmer_mask = (1u << (2 * k)) - 1;
            unsigned int kmer_hash = 0;
            unsigned int valid_bases = 0;
            
            // Count k-mers in this sequence
            for (unsigned int i = 0; i < info.length; i++) {
                unsigned char base = BASE_ENCODING[(unsigned char)seq[i]];
                
                if (base < 4) {
                    kmer_hash = ((kmer_hash << 2) | base) & kmer_mask;
                    valid_bases = min(valid_bases + 1, k);
                    
                    if (valid_bases >= k && kmer_hash < max_kmers) {
                        atomicAdd(&kmer_counts[kmer_hash], 1);
                        local_kmer_count++;
                    }
                } else {
                    kmer_hash = 0;
                    valid_bases = 0;
                }
            }
        }
    }
    
    // Use CUB for efficient block-level reduction (for statistics)
    unsigned int block_total = BlockReduce(temp_storage).Sum(local_kmer_count);
    
    // Block leader writes statistics (could be used for progress monitoring)
    if (tid == 0 && block_total > 0) {
        // Statistics collection could be added here
    }
}

// Utility kernel for clearing k-mer count buffer
__global__ void clear_kmer_counts_kernel(
    unsigned int* __restrict__ kmer_counts,
    const unsigned int max_kmers
) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < max_kmers) {
        kmer_counts[idx] = 0;
    }
}

// Kernel for prefixing sequence data for better memory layout
__global__ void prepare_sequences_kernel(
    const char* __restrict__ input_sequences,
    const unsigned int* __restrict__ input_lengths,
    char* __restrict__ output_sequences,
    SequenceInfo* __restrict__ seq_info,
    const unsigned int num_sequences
) {
    unsigned int seq_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (seq_idx >= num_sequences) {
        return;
    }
    
    // Each thread recomputes its output offset with a serial O(n) sum over
    // preceding lengths; a CUB prefix scan would avoid the quadratic total work.
    unsigned int output_offset = 0;
    for (unsigned int i = 0; i < seq_idx; i++) {
        output_offset += input_lengths[i] + 1;  // +1 for null terminator
    }
    
    // Copy sequence data
    unsigned int input_offset = 0;
    for (unsigned int i = 0; i < seq_idx; i++) {
        input_offset += input_lengths[i] + 1;
    }
    
    unsigned int length = input_lengths[seq_idx];
    
    // Copy sequence with validation and normalization
    for (unsigned int i = 0; i < length; i++) {
        char base = input_sequences[input_offset + i];
        
        // Normalize to uppercase
        if (base >= 'a' && base <= 'z') {
            base -= 32;
        }
        
        // Validate DNA base
        if (base != 'A' && base != 'T' && base != 'G' && base != 'C') {
            base = 'N';  // Replace invalid bases with N
        }
        
        output_sequences[output_offset + i] = base;
    }
    
    // Add null terminator
    output_sequences[output_offset + length] = '\0';
    
    // Fill sequence info
    seq_info[seq_idx].offset = output_offset;
    seq_info[seq_idx].length = length;
    seq_info[seq_idx].seq_id = seq_idx;
    seq_info[seq_idx].padding = 0;
}