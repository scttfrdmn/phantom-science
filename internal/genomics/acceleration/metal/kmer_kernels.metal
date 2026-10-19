#include <metal_stdlib>
using namespace metal;

// Metal kernels for parallel k-mer extraction and counting (one thread per
// sequence, atomic increments into a dense count table). No benchmark in the
// repository measures their performance.

// Base encoding: A=0, T=1, G=2, C=3
constant uint8_t BASE_ENCODING[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0-15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 16-31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 32-47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 48-63
    4, 0, 4, 3, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 64-79  (A=65->0, C=67->3, G=71->2)
    4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 80-95  (T=84->1)
    4, 0, 4, 3, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 96-111 (a=97->0, c=99->3, g=103->2)
    4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 112-127 (t=116->1)
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 128-143
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 144-159
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 160-175
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 176-191
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 192-207
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 208-223
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 224-239
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4   // 240-255
};

// Struct for passing sequence information to GPU
struct SequenceInfo {
    uint32_t offset;     // Offset in sequence buffer
    uint32_t length;     // Length of this sequence
    uint32_t seq_id;     // Sequence identifier
    uint32_t padding;    // Alignment padding
};

// Kernel for extracting k-mers from DNA sequences
kernel void extract_kmers(
    constant char* sequences [[buffer(0)]],                    // All sequences concatenated
    constant SequenceInfo* seq_info [[buffer(1)]],            // Sequence metadata
    device atomic_uint* kmer_counts [[buffer(2)]],            // Output k-mer frequency counts
    constant uint32_t& k [[buffer(3)]],                       // K-mer size
    constant uint32_t& num_sequences [[buffer(4)]],           // Number of sequences
    constant uint32_t& max_kmers [[buffer(5)]],               // Maximum possible k-mers (4^k)
    uint thread_position_in_grid [[thread_position_in_grid]], // Global thread ID
    uint threads_per_threadgroup [[threads_per_threadgroup]], // Threads per group
    uint threadgroup_position_in_grid [[threadgroup_position_in_grid]] // Threadgroup ID
) {
    // Early exit if thread is out of range
    if (thread_position_in_grid >= num_sequences) {
        return;
    }
    
    uint32_t seq_idx = thread_position_in_grid;
    SequenceInfo info = seq_info[seq_idx];
    
    // Skip sequences shorter than k
    if (info.length < k) {
        return;
    }
    
    // Calculate k-mer mask for efficient bit manipulation
    uint32_t kmer_mask = (1u << (2 * k)) - 1;
    uint32_t kmer_hash = 0;
    uint32_t valid_bases = 0;
    
    // Process sequence starting at offset
    constant char* seq = sequences + info.offset;
    uint32_t num_kmers = info.length - k + 1;
    
    // Initialize sliding window
    for (uint32_t i = 0; i < k && i < info.length; i++) {
        uint8_t base = BASE_ENCODING[uint8_t(seq[i])];
        if (base < 4) {  // Valid base (A, T, G, C)
            kmer_hash = (kmer_hash << 2) | base;
            valid_bases++;
        } else {  // Invalid base (N or other)
            kmer_hash = 0;
            valid_bases = 0;
        }
    }
    
    // Extract first k-mer if we have k valid bases
    if (valid_bases >= k && kmer_hash < max_kmers) {
        atomic_fetch_add_explicit(&kmer_counts[kmer_hash], 1, memory_order_relaxed);
    }
    
    // Slide window through remaining sequence
    for (uint32_t i = k; i < info.length; i++) {
        uint8_t base = BASE_ENCODING[uint8_t(seq[i])];
        
        if (base < 4) {  // Valid base
            // Slide window: remove leftmost base, add rightmost base
            kmer_hash = ((kmer_hash << 2) | base) & kmer_mask;
            valid_bases = min(valid_bases + 1, k);
            
            // Count k-mer if we have k consecutive valid bases
            if (valid_bases >= k && kmer_hash < max_kmers) {
                atomic_fetch_add_explicit(&kmer_counts[kmer_hash], 1, memory_order_relaxed);
            }
        } else {  // Invalid base - reset window
            kmer_hash = 0;
            valid_bases = 0;
        }
    }
}

// Threadgroup histogram size for extract_kmers_small_k: 4096 counters
// (16 KB) cover every k-mer for k <= 6 and fit the 32 KB threadgroup memory
// limit of Apple GPUs (MTLDevice.maxThreadgroupMemoryLength). A 4^8-entry
// table (256 KB) cannot be allocated: the pipeline fails to build, and with
// it the whole Metal device setup in metal_bridge_darwin.m.
constant uint SMALL_K_LOCAL_BINS = 4096;

// Count one k-mer, in the threadgroup histogram when the table fits there
// and directly in device memory otherwise (k = 7, 8).
inline void count_small_k(threadgroup atomic_uint* local_counts,
                          device atomic_uint* kmer_counts,
                          bool use_local, uint32_t kmer_hash, uint max_kmers) {
    if (kmer_hash >= max_kmers) {
        return;
    }
    if (use_local) {
        atomic_fetch_add_explicit(&local_counts[kmer_hash], 1, memory_order_relaxed);
    } else {
        atomic_fetch_add_explicit(&kmer_counts[kmer_hash], 1, memory_order_relaxed);
    }
}

// Kernel for small k-mers (k <= 8), one thread per sequence. For k <= 6 the
// counts are accumulated per threadgroup and then added to device memory.
kernel void extract_kmers_small_k(
    constant char* sequences [[buffer(0)]],
    constant SequenceInfo* seq_info [[buffer(1)]],
    device atomic_uint* kmer_counts [[buffer(2)]],
    constant uint32_t& k [[buffer(3)]],
    constant uint32_t& num_sequences [[buffer(4)]],
    constant uint32_t& max_kmers [[buffer(5)]],
    uint thread_position_in_grid [[thread_position_in_grid]],
    uint thread_position_in_threadgroup [[thread_position_in_threadgroup]],
    uint threads_per_threadgroup [[threads_per_threadgroup]]
) {
    threadgroup atomic_uint local_counts[SMALL_K_LOCAL_BINS];

    // Uniform across the threadgroup, so the barriers below are reached by
    // every thread.
    bool use_local = max_kmers <= SMALL_K_LOCAL_BINS;
    uint local_max = use_local ? max_kmers : 0;
    for (uint i = thread_position_in_threadgroup; i < local_max; i += threads_per_threadgroup) {
        atomic_store_explicit(&local_counts[i], 0, memory_order_relaxed);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Process the sequence assigned to this thread
    if (thread_position_in_grid < num_sequences) {
        uint32_t seq_idx = thread_position_in_grid;
        SequenceInfo info = seq_info[seq_idx];

        if (info.length >= k) {
            constant char* seq = sequences + info.offset;
            uint32_t kmer_mask = (1u << (2 * k)) - 1;
            uint32_t kmer_hash = 0;
            uint32_t valid_bases = 0;

            // Initialize sliding window
            for (uint32_t i = 0; i < k && i < info.length; i++) {
                uint8_t base = BASE_ENCODING[uint8_t(seq[i])];
                if (base < 4) {
                    kmer_hash = (kmer_hash << 2) | base;
                    valid_bases++;
                } else {
                    kmer_hash = 0;
                    valid_bases = 0;
                }
            }

            // Count first k-mer
            if (valid_bases >= k) {
                count_small_k(local_counts, kmer_counts, use_local, kmer_hash, max_kmers);
            }

            // Slide through sequence
            for (uint32_t i = k; i < info.length; i++) {
                uint8_t base = BASE_ENCODING[uint8_t(seq[i])];

                if (base < 4) {
                    kmer_hash = ((kmer_hash << 2) | base) & kmer_mask;
                    valid_bases = min(valid_bases + 1, k);

                    if (valid_bases >= k) {
                        count_small_k(local_counts, kmer_counts, use_local, kmer_hash, max_kmers);
                    }
                } else {
                    kmer_hash = 0;
                    valid_bases = 0;
                }
            }
        }
    }

    // Synchronize threadgroup before writing to global memory
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Add the threadgroup's counts to global memory
    for (uint i = thread_position_in_threadgroup; i < local_max; i += threads_per_threadgroup) {
        uint c = atomic_load_explicit(&local_counts[i], memory_order_relaxed);
        if (c > 0) {
            atomic_fetch_add_explicit(&kmer_counts[i], c, memory_order_relaxed);
        }
    }
}

// Kernel for large datasets (selected above 100 MB of sequence) that strides
// threads across positions instead of sequences
kernel void extract_kmers_coalesced(
    constant char* sequences [[buffer(0)]],
    constant SequenceInfo* seq_info [[buffer(1)]],
    device atomic_uint* kmer_counts [[buffer(2)]],
    constant uint32_t& k [[buffer(3)]],
    constant uint32_t& total_bases [[buffer(4)]],
    constant uint32_t& max_kmers [[buffer(5)]],
    uint thread_position_in_grid [[thread_position_in_grid]],
    uint threads_per_grid [[threads_per_grid]]
) {
    // Each thread processes a stride of bases across all sequences
    uint32_t kmer_mask = (1u << (2 * k)) - 1;
    uint32_t stride = threads_per_grid;
    
    for (uint32_t base_idx = thread_position_in_grid; base_idx + k <= total_bases; base_idx += stride) {
        uint32_t kmer_hash = 0;
        bool valid_kmer = true;
        
        // Build k-mer from k consecutive bases
        for (uint32_t i = 0; i < k; i++) {
            uint8_t base = BASE_ENCODING[uint8_t(sequences[base_idx + i])];
            if (base >= 4) {  // Invalid base
                valid_kmer = false;
                break;
            }
            kmer_hash = (kmer_hash << 2) | base;
        }
        
        // Count valid k-mer
        if (valid_kmer && kmer_hash < max_kmers) {
            atomic_fetch_add_explicit(&kmer_counts[kmer_hash], 1, memory_order_relaxed);
        }
    }
}

// Utility kernel for clearing k-mer count buffer
kernel void clear_kmer_counts(
    device atomic_uint* kmer_counts [[buffer(0)]],
    constant uint32_t& max_kmers [[buffer(1)]],
    uint thread_position_in_grid [[thread_position_in_grid]]
) {
    if (thread_position_in_grid < max_kmers) {
        atomic_store_explicit(&kmer_counts[thread_position_in_grid], 0, memory_order_relaxed);
    }
}

// Kernel for copying atomic results to regular buffer for CPU readback
kernel void copy_results(
    device atomic_uint* kmer_counts [[buffer(0)]],
    device uint* output_counts [[buffer(1)]],
    constant uint32_t& max_kmers [[buffer(2)]],
    uint thread_position_in_grid [[thread_position_in_grid]]
) {
    if (thread_position_in_grid < max_kmers) {
        output_counts[thread_position_in_grid] = atomic_load_explicit(
            &kmer_counts[thread_position_in_grid], 
            memory_order_relaxed
        );
    }
}