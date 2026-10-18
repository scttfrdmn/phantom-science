#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include "metal_bridge.h"

// Error handling macros
#define CHECK_ERROR(error) \
    if (error != nil) { \
        NSLog(@"Metal Error: %@", error.localizedDescription); \
        return NULL; \
    }

#define LOG_INFO(format, ...) \
    NSLog(@"Metal GPU: " format, ##__VA_ARGS__)

// Metal device structure implementation
struct MetalDevice {
    id<MTLDevice> device;
    id<MTLCommandQueue> commandQueue;
    id<MTLLibrary> library;
    id<MTLComputePipelineState> extractKmersPipeline;
    id<MTLComputePipelineState> extractKmersSmallKPipeline;
    id<MTLComputePipelineState> extractKmersCoalescedPipeline;
    id<MTLComputePipelineState> clearCountsPipeline;
    id<MTLComputePipelineState> copyResultsPipeline;
    NSString* deviceName;
    NSUInteger maxThreadsPerGroup;
    NSUInteger maxBufferLength;
    bool isInitialized;
};

// Check if Metal is available on this system
bool isMetalAvailable() {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        return device != nil;
    }
}

// Initialize Metal device and compile compute kernels
MetalDevice* initializeMetalDevice() {
    @autoreleasepool {
        NSError* error = nil;
        
        // Create Metal device
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            NSLog(@"Metal Error: No Metal-compatible device found");
            return NULL;
        }
        
        LOG_INFO("Initializing Metal device: %@", device.name);
        
        // Create command queue
        id<MTLCommandQueue> commandQueue = [device newCommandQueue];
        if (!commandQueue) {
            NSLog(@"Metal Error: Failed to create command queue");
            return NULL;
        }
        
        // Locate and load the Metal shader library
        NSBundle* bundle = [NSBundle mainBundle];
        NSString* libraryPath = [bundle pathForResource:@"kmer_kernels" ofType:@"metallib"];
        
        id<MTLLibrary> library = nil;
        if (libraryPath) {
            // Try to load pre-compiled Metal library
            NSURL* libraryURL = [NSURL fileURLWithPath:libraryPath];
            library = [device newLibraryWithURL:libraryURL error:&error];
        }
        
        if (!library) {
            // Fallback: compile from source (for development)
            NSString* sourcePath = [[NSBundle mainBundle] pathForResource:@"kmer_kernels" ofType:@"metal"];
            if (!sourcePath) {
                // Try relative path for development builds
                NSString* currentDir = [[NSFileManager defaultManager] currentDirectoryPath];
                sourcePath = [currentDir stringByAppendingPathComponent:@"internal/genomics/acceleration/metal/kmer_kernels.metal"];
            }
            
            if (sourcePath && [[NSFileManager defaultManager] fileExistsAtPath:sourcePath]) {
                NSString* source = [NSString stringWithContentsOfFile:sourcePath 
                                                            encoding:NSUTF8StringEncoding 
                                                               error:&error];
                if (source) {
                    LOG_INFO("Compiling Metal shaders from source: %@", sourcePath);
                    library = [device newLibraryWithSource:source options:nil error:&error];
                }
            }
        }
        
        if (!library) {
            NSLog(@"Metal Error: Failed to load Metal library: %@", error ? error.localizedDescription : @"Unknown error");
            return NULL;
        }
        
        // Create compute pipeline states for each kernel
        id<MTLFunction> extractKmersFunction = [library newFunctionWithName:@"extract_kmers"];
        id<MTLFunction> extractKmersSmallKFunction = [library newFunctionWithName:@"extract_kmers_small_k"];
        id<MTLFunction> extractKmersCoalescedFunction = [library newFunctionWithName:@"extract_kmers_coalesced"];
        id<MTLFunction> clearCountsFunction = [library newFunctionWithName:@"clear_kmer_counts"];
        id<MTLFunction> copyResultsFunction = [library newFunctionWithName:@"copy_results"];
        
        if (!extractKmersFunction || !extractKmersSmallKFunction || !extractKmersCoalescedFunction || 
            !clearCountsFunction || !copyResultsFunction) {
            NSLog(@"Metal Error: Failed to load required compute functions from library");
            return NULL;
        }
        
        id<MTLComputePipelineState> extractKmersPipeline = 
            [device newComputePipelineStateWithFunction:extractKmersFunction error:&error];
        CHECK_ERROR(error);
        
        id<MTLComputePipelineState> extractKmersSmallKPipeline = 
            [device newComputePipelineStateWithFunction:extractKmersSmallKFunction error:&error];
        CHECK_ERROR(error);
        
        id<MTLComputePipelineState> extractKmersCoalescedPipeline = 
            [device newComputePipelineStateWithFunction:extractKmersCoalescedFunction error:&error];
        CHECK_ERROR(error);
        
        id<MTLComputePipelineState> clearCountsPipeline = 
            [device newComputePipelineStateWithFunction:clearCountsFunction error:&error];
        CHECK_ERROR(error);
        
        id<MTLComputePipelineState> copyResultsPipeline = 
            [device newComputePipelineStateWithFunction:copyResultsFunction error:&error];
        CHECK_ERROR(error);
        
        // Allocate and initialize MetalDevice structure
        MetalDevice* metalDevice = (MetalDevice*)malloc(sizeof(MetalDevice));
        if (!metalDevice) {
            NSLog(@"Metal Error: Failed to allocate MetalDevice structure");
            return NULL;
        }
        
        metalDevice->device = device;
        metalDevice->commandQueue = commandQueue;
        metalDevice->library = library;
        metalDevice->extractKmersPipeline = extractKmersPipeline;
        metalDevice->extractKmersSmallKPipeline = extractKmersSmallKPipeline;
        metalDevice->extractKmersCoalescedPipeline = extractKmersCoalescedPipeline;
        metalDevice->clearCountsPipeline = clearCountsPipeline;
        metalDevice->copyResultsPipeline = copyResultsPipeline;
        metalDevice->deviceName = [device.name copy];
        metalDevice->maxThreadsPerGroup = extractKmersPipeline.maxTotalThreadsPerThreadgroup;
        metalDevice->maxBufferLength = device.maxBufferLength;
        metalDevice->isInitialized = true;
        
        // Retain Objective-C objects to prevent deallocation
        [device retain];
        [commandQueue retain];
        [library retain];
        [extractKmersPipeline retain];
        [extractKmersSmallKPipeline retain];
        [extractKmersCoalescedPipeline retain];
        [clearCountsPipeline retain];
        [copyResultsPipeline retain];
        [metalDevice->deviceName retain];
        
        LOG_INFO("Metal device initialized successfully");
        LOG_INFO("Device: %@", metalDevice->deviceName);
        LOG_INFO("Max threads per threadgroup: %lu", (unsigned long)metalDevice->maxThreadsPerGroup);
        LOG_INFO("Max buffer length: %lu MB", (unsigned long)(metalDevice->maxBufferLength / (1024 * 1024)));
        
        return metalDevice;
    }
}

// Process k-mers on GPU using Metal compute shaders
int processKmersOnGPU(MetalDevice* metalDevice, const char* sequences, 
                      int numSequences, int* sequenceLengths, 
                      int k, int* results, int maxResults) {
    @autoreleasepool {
        if (!metalDevice || !metalDevice->isInitialized) {
            NSLog(@"Metal Error: Device not initialized");
            return -1;
        }
        
        if (!sequences || !sequenceLengths || !results || k <= 0 || k > 16) {
            NSLog(@"Metal Error: Invalid parameters");
            return -2;
        }
        
        NSError* error = nil;
        id<MTLDevice> device = metalDevice->device;
        id<MTLCommandQueue> commandQueue = metalDevice->commandQueue;
        
        // Calculate total sequence data size
        size_t totalSequenceSize = 0;
        for (int i = 0; i < numSequences; i++) {
            totalSequenceSize += sequenceLengths[i] + 1; // +1 for null terminator
        }
        
        LOG_INFO("Processing %d sequences (%zu bytes) with k=%d", numSequences, totalSequenceSize, k);
        
        // Create sequence info structures
        typedef struct {
            uint32_t offset;
            uint32_t length;
            uint32_t seq_id;
            uint32_t padding;
        } SequenceInfo;
        
        SequenceInfo* sequenceInfos = (SequenceInfo*)malloc(numSequences * sizeof(SequenceInfo));
        if (!sequenceInfos) {
            NSLog(@"Metal Error: Failed to allocate sequence info array");
            return -3;
        }
        
        uint32_t offset = 0;
        for (int i = 0; i < numSequences; i++) {
            sequenceInfos[i].offset = offset;
            sequenceInfos[i].length = sequenceLengths[i];
            sequenceInfos[i].seq_id = i;
            sequenceInfos[i].padding = 0;
            offset += sequenceLengths[i] + 1;
        }
        
        // Create Metal buffers
        id<MTLBuffer> sequenceBuffer = [device newBufferWithBytes:sequences 
                                                           length:totalSequenceSize 
                                                          options:MTLResourceStorageModeShared];
        
        id<MTLBuffer> sequenceInfoBuffer = [device newBufferWithBytes:sequenceInfos 
                                                               length:numSequences * sizeof(SequenceInfo)
                                                              options:MTLResourceStorageModeShared];
        
        id<MTLBuffer> kmerCountsBuffer = [device newBufferWithLength:maxResults * sizeof(uint32_t)
                                                             options:MTLResourceStorageModeShared];
        
        id<MTLBuffer> outputBuffer = [device newBufferWithLength:maxResults * sizeof(uint32_t)
                                                         options:MTLResourceStorageModeShared];
        
        if (!sequenceBuffer || !sequenceInfoBuffer || !kmerCountsBuffer || !outputBuffer) {
            NSLog(@"Metal Error: Failed to create Metal buffers");
            free(sequenceInfos);
            return -4;
        }
        
        // Create command buffer
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        if (!commandBuffer) {
            NSLog(@"Metal Error: Failed to create command buffer");
            free(sequenceInfos);
            return -5;
        }
        
        // Clear k-mer counts buffer
        id<MTLComputeCommandEncoder> clearEncoder = [commandBuffer computeCommandEncoder];
        [clearEncoder setComputePipelineState:metalDevice->clearCountsPipeline];
        [clearEncoder setBuffer:kmerCountsBuffer offset:0 atIndex:0];
        [clearEncoder setBytes:&maxResults length:sizeof(uint32_t) atIndex:1];
        
        NSUInteger clearThreadsPerGroup = MIN(256, metalDevice->maxThreadsPerGroup);
        NSUInteger clearNumGroups = (maxResults + clearThreadsPerGroup - 1) / clearThreadsPerGroup;
        
        [clearEncoder dispatchThreadgroups:MTLSizeMake(clearNumGroups, 1, 1)
                      threadsPerThreadgroup:MTLSizeMake(clearThreadsPerGroup, 1, 1)];
        [clearEncoder endEncoding];
        
        // Choose optimal kernel based on k-mer size and dataset characteristics
        id<MTLComputePipelineState> pipeline;
        NSUInteger threadsPerGroup;
        NSUInteger numGroups;
        
        if (k <= 8 && numSequences <= 10000) {
            // Use small-k optimized kernel with shared memory
            pipeline = metalDevice->extractKmersSmallKPipeline;
            threadsPerGroup = MIN(256, metalDevice->maxThreadsPerGroup);
            numGroups = (numSequences + threadsPerGroup - 1) / threadsPerGroup;
            LOG_INFO("Using small-k optimized kernel");
        } else if (totalSequenceSize > 100000000) {  // >100MB
            // Use coalesced memory access kernel for large datasets
            pipeline = metalDevice->extractKmersCoalescedPipeline;
            threadsPerGroup = MIN(512, metalDevice->maxThreadsPerGroup);
            numGroups = (totalSequenceSize + threadsPerGroup - 1) / threadsPerGroup;
            LOG_INFO("Using coalesced memory access kernel");
        } else {
            // Use standard kernel
            pipeline = metalDevice->extractKmersPipeline;
            threadsPerGroup = MIN(256, metalDevice->maxThreadsPerGroup);
            numGroups = (numSequences + threadsPerGroup - 1) / threadsPerGroup;
            LOG_INFO("Using standard kernel");
        }
        
        // Encode k-mer extraction kernel
        id<MTLComputeCommandEncoder> extractEncoder = [commandBuffer computeCommandEncoder];
        [extractEncoder setComputePipelineState:pipeline];
        [extractEncoder setBuffer:sequenceBuffer offset:0 atIndex:0];
        [extractEncoder setBuffer:sequenceInfoBuffer offset:0 atIndex:1];
        [extractEncoder setBuffer:kmerCountsBuffer offset:0 atIndex:2];
        [extractEncoder setBytes:&k length:sizeof(uint32_t) atIndex:3];
        
        if (pipeline == metalDevice->extractKmersCoalescedPipeline) {
            uint32_t totalBases = (uint32_t)totalSequenceSize;
            [extractEncoder setBytes:&totalBases length:sizeof(uint32_t) atIndex:4];
        } else {
            uint32_t numSeq = (uint32_t)numSequences;
            [extractEncoder setBytes:&numSeq length:sizeof(uint32_t) atIndex:4];
        }
        
        [extractEncoder setBytes:&maxResults length:sizeof(uint32_t) atIndex:5];
        
        [extractEncoder dispatchThreadgroups:MTLSizeMake(numGroups, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(threadsPerGroup, 1, 1)];
        [extractEncoder endEncoding];
        
        // Copy atomic results to output buffer
        id<MTLComputeCommandEncoder> copyEncoder = [commandBuffer computeCommandEncoder];
        [copyEncoder setComputePipelineState:metalDevice->copyResultsPipeline];
        [copyEncoder setBuffer:kmerCountsBuffer offset:0 atIndex:0];
        [copyEncoder setBuffer:outputBuffer offset:0 atIndex:1];
        [copyEncoder setBytes:&maxResults length:sizeof(uint32_t) atIndex:2];
        
        NSUInteger copyThreadsPerGroup = MIN(256, metalDevice->maxThreadsPerGroup);
        NSUInteger copyNumGroups = (maxResults + copyThreadsPerGroup - 1) / copyThreadsPerGroup;
        
        [copyEncoder dispatchThreadgroups:MTLSizeMake(copyNumGroups, 1, 1)
                     threadsPerThreadgroup:MTLSizeMake(copyThreadsPerGroup, 1, 1)];
        [copyEncoder endEncoding];
        
        // Execute command buffer and wait for completion
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];
        
        if (commandBuffer.error) {
            NSLog(@"Metal Error: Command buffer execution failed: %@", commandBuffer.error.localizedDescription);
            free(sequenceInfos);
            return -6;
        }
        
        // Copy results back to CPU
        uint32_t* outputData = (uint32_t*)[outputBuffer contents];
        int totalKmers = 0;
        for (int i = 0; i < maxResults; i++) {
            results[i] = (int)outputData[i];
            totalKmers += results[i];
        }
        
        LOG_INFO("GPU k-mer extraction completed: %d total k-mers processed", totalKmers);
        
        free(sequenceInfos);
        return totalKmers;
    }
}

// Get Metal device name
const char* getMetalDeviceName(MetalDevice* metalDevice) {
    if (!metalDevice || !metalDevice->isInitialized) {
        return "Unknown";
    }
    return [metalDevice->deviceName UTF8String];
}

// Get maximum threads per threadgroup
int getMetalMaxThreadsPerGroup(MetalDevice* metalDevice) {
    if (!metalDevice || !metalDevice->isInitialized) {
        return 0;
    }
    return (int)metalDevice->maxThreadsPerGroup;
}

// Get maximum buffer length
unsigned long getMetalMaxBufferLength(MetalDevice* metalDevice) {
    if (!metalDevice || !metalDevice->isInitialized) {
        return 0;
    }
    return (unsigned long)metalDevice->maxBufferLength;
}

// Release Metal device and associated resources
void releaseMetalDevice(MetalDevice* metalDevice) {
    @autoreleasepool {
        if (!metalDevice) {
            return;
        }
        
        LOG_INFO("Releasing Metal device resources");
        
        if (metalDevice->isInitialized) {
            // Release Objective-C objects
            [metalDevice->device release];
            [metalDevice->commandQueue release];
            [metalDevice->library release];
            [metalDevice->extractKmersPipeline release];
            [metalDevice->extractKmersSmallKPipeline release];
            [metalDevice->extractKmersCoalescedPipeline release];
            [metalDevice->clearCountsPipeline release];
            [metalDevice->copyResultsPipeline release];
            [metalDevice->deviceName release];
            
            metalDevice->isInitialized = false;
        }
        
        free(metalDevice);
    }
}