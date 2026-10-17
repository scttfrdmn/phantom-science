#ifndef METAL_BRIDGE_H
#define METAL_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

// Forward declaration of MetalDevice structure
typedef struct MetalDevice MetalDevice;

// Check if Metal is available on this system
bool isMetalAvailable(void);

// Initialize Metal device and load compute shaders
MetalDevice* initializeMetalDevice(void);

// Process k-mers on GPU using Metal compute shaders
// Returns the total number of k-mers processed, or negative value on error
int processKmersOnGPU(MetalDevice* device, const char* sequences, 
                      int numSequences, int* sequenceLengths, 
                      int k, int* results, int maxResults);

// Device capability queries
const char* getMetalDeviceName(MetalDevice* device);
int getMetalMaxThreadsPerGroup(MetalDevice* device);
unsigned long getMetalMaxBufferLength(MetalDevice* device);

// Release Metal device and associated resources
void releaseMetalDevice(MetalDevice* device);

#ifdef __cplusplus
}
#endif

#endif // METAL_BRIDGE_H