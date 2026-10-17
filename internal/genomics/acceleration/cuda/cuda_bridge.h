#ifndef CUDA_BRIDGE_H
#define CUDA_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

// Forward declaration of CUDA device structure
typedef struct CUDADevice CUDADevice;

// Check if CUDA is available on this system
bool isCudaAvailable(void);

// Initialize CUDA device and prepare for compute operations
CUDADevice* initializeCudaDevice(void);

// Process k-mers on GPU using CUDA compute kernels
// Returns the total number of k-mers processed, or negative value on error
int processKmersOnCuda(CUDADevice* device, const char* sequences,
                       int numSequences, int* sequenceLengths,
                       int k, int* results, int maxResults);

// Device capability queries
const char* getCudaDeviceName(CUDADevice* device);
int getCudaMaxThreadsPerBlock(CUDADevice* device);
unsigned long getCudaGlobalMemorySize(CUDADevice* device);
int getCudaComputeCapability(CUDADevice* device);
int getCudaMultiprocessorCount(CUDADevice* device);

// Release CUDA device and associated resources
void releaseCudaDevice(CUDADevice* device);

#ifdef __cplusplus
}
#endif

#endif // CUDA_BRIDGE_H