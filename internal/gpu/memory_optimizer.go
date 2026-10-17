// Package gpu provides GPU memory optimization for large dataset processing
package gpu

import (
	"fmt"
	"log"
	"runtime"
	"sync"
	"time"
)

// MemoryPool manages GPU memory allocation with optimization for large datasets
type MemoryPool struct {
	maxMemoryGB    float64
	currentUsageGB float64
	pools          map[string]*DevicePool
	mu             sync.RWMutex
	config         *OptimizationConfig
	metrics        *MemoryMetrics
}

// DevicePool manages memory for a specific GPU device
type DevicePool struct {
	deviceID       int
	totalMemoryGB  float64
	availableGB    float64
	allocations    map[string]*Allocation
	lastUsed       time.Time
	mu             sync.RWMutex
}

// Allocation represents a GPU memory allocation
type Allocation struct {
	ID           string
	SizeGB       float64
	Purpose      string
	CreatedAt    time.Time
	LastAccessed time.Time
	Priority     Priority
	DeviceID     int
}

// Priority defines allocation priority levels
type Priority int

const (
	PriorityLow Priority = iota
	PriorityNormal
	PriorityHigh
	PriorityCritical
)

// OptimizationConfig defines GPU memory optimization parameters
type OptimizationConfig struct {
	// Memory management
	MaxMemoryUsagePercent float64 // Max % of GPU memory to use
	PreallocationSizeGB   float64 // Memory to preallocate for performance
	GCThresholdPercent    float64 // Trigger garbage collection threshold
	
	// Performance tuning
	StreamingEnabled      bool    // Enable streaming for large datasets
	ChunkSizeGB          float64 // Size of data chunks for streaming
	ConcurrentStreams    int     // Number of concurrent CUDA streams
	
	// Caching
	CacheEnabled         bool    // Enable intelligent caching
	CacheMaxSizeGB       float64 // Maximum cache size
	CacheTTL             time.Duration // Cache time-to-live
	
	// Auto-optimization
	AdaptiveEnabled      bool    // Enable adaptive optimization
	ProfilingEnabled     bool    // Enable performance profiling
	OptimizationInterval time.Duration // How often to optimize
}

// MemoryMetrics tracks GPU memory usage statistics
type MemoryMetrics struct {
	TotalAllocationsGB   float64
	PeakUsageGB         float64
	AllocationCount     int
	EvictionCount       int
	CacheHitRate        float64
	AverageLatency      time.Duration
	LastOptimization    time.Time
	mu                  sync.RWMutex
}

// DefaultOptimizationConfig returns enterprise-grade default configuration
func DefaultOptimizationConfig() *OptimizationConfig {
	return &OptimizationConfig{
		MaxMemoryUsagePercent: 85.0, // Leave 15% for system
		PreallocationSizeGB:   2.0,  // 2GB preallocated pool
		GCThresholdPercent:    75.0, // GC when 75% full
		
		StreamingEnabled:      true,
		ChunkSizeGB:          1.0,   // 1GB chunks for streaming
		ConcurrentStreams:    4,     // 4 concurrent streams
		
		CacheEnabled:         true,
		CacheMaxSizeGB:       4.0,   // 4GB cache
		CacheTTL:             30 * time.Minute,
		
		AdaptiveEnabled:      true,
		ProfilingEnabled:     true,
		OptimizationInterval: 5 * time.Minute,
	}
}

// NewMemoryPool creates a new GPU memory pool manager
func NewMemoryPool(maxMemoryGB float64, config *OptimizationConfig) (*MemoryPool, error) {
	if config == nil {
		config = DefaultOptimizationConfig()
	}
	
	pool := &MemoryPool{
		maxMemoryGB: maxMemoryGB,
		pools:       make(map[string]*DevicePool),
		config:      config,
		metrics:     &MemoryMetrics{},
	}
	
	// Initialize GPU devices
	if err := pool.initializeDevices(); err != nil {
		return nil, fmt.Errorf("failed to initialize GPU devices: %w", err)
	}
	
	// Start background optimization if enabled
	if config.AdaptiveEnabled {
		go pool.optimizationLoop()
	}
	
	return pool, nil
}

// initializeDevices discovers and initializes available GPU devices
func (mp *MemoryPool) initializeDevices() error {
	// In a real implementation, this would query CUDA/OpenCL for available devices
	// For now, simulate detection of available GPUs
	
	// Simulated GPU detection (replace with actual NVIDIA ML/CUDA calls)
	gpuCount := runtime.NumCPU() // Placeholder - would be nvidia-ml-py equivalent
	if gpuCount > 8 {
		gpuCount = 8 // Reasonable limit for enterprise
	}
	
	for i := 0; i < gpuCount; i++ {
		deviceMemoryGB := 24.0 // NVIDIA TITAN RTX has 24GB
		
		device := &DevicePool{
			deviceID:      i,
			totalMemoryGB: deviceMemoryGB,
			availableGB:   deviceMemoryGB * (mp.config.MaxMemoryUsagePercent / 100.0),
			allocations:   make(map[string]*Allocation),
			lastUsed:      time.Now(),
		}
		
		mp.pools[fmt.Sprintf("cuda:%d", i)] = device
		log.Printf("Initialized GPU device %d with %.1fGB available memory", i, device.availableGB)
	}
	
	return nil
}

// AllocateMemory allocates GPU memory with optimization
func (mp *MemoryPool) AllocateMemory(sizeGB float64, purpose string, priority Priority) (*Allocation, error) {
	mp.mu.Lock()
	defer mp.mu.Unlock()
	
	// Find best device for allocation
	deviceKey, device, err := mp.findBestDevice(sizeGB)
	if err != nil {
		// Try to free memory and retry
		if mp.tryFreeMemory(sizeGB) {
			deviceKey, device, err = mp.findBestDevice(sizeGB)
		}
		if err != nil {
			return nil, fmt.Errorf("cannot allocate %.2fGB: %w", sizeGB, err)
		}
	}
	
	// Create allocation
	allocation := &Allocation{
		ID:           generateAllocationID(),
		SizeGB:       sizeGB,
		Purpose:      purpose,
		CreatedAt:    time.Now(),
		LastAccessed: time.Now(),
		Priority:     priority,
		DeviceID:     device.deviceID,
	}
	
	// Update device state
	device.mu.Lock()
	device.allocations[allocation.ID] = allocation
	device.availableGB -= sizeGB
	device.lastUsed = time.Now()
	device.mu.Unlock()
	
	// Update metrics
	mp.updateMetrics(allocation, true)
	
	log.Printf("Allocated %.2fGB on %s for %s (ID: %s)", sizeGB, deviceKey, purpose, allocation.ID)
	return allocation, nil
}

// findBestDevice finds the optimal GPU device for allocation
func (mp *MemoryPool) findBestDevice(sizeGB float64) (string, *DevicePool, error) {
	var bestDevice *DevicePool
	var bestKey string
	var bestScore float64
	
	for key, device := range mp.pools {
		device.mu.RLock()
		
		// Check if device has enough memory
		if device.availableGB < sizeGB {
			device.mu.RUnlock()
			continue
		}
		
		// Calculate allocation score (higher is better)
		utilizationRatio := (device.totalMemoryGB - device.availableGB) / device.totalMemoryGB
		timeSinceLastUse := time.Since(device.lastUsed).Seconds()
		
		// Prefer less utilized devices with recent usage
		score := (1.0 - utilizationRatio) * 0.7 + (1.0/(timeSinceLastUse+1)) * 0.3
		
		device.mu.RUnlock()
		
		if bestDevice == nil || score > bestScore {
			bestDevice = device
			bestKey = key
			bestScore = score
		}
	}
	
	if bestDevice == nil {
		return "", nil, fmt.Errorf("no device found with %.2fGB available memory", sizeGB)
	}
	
	return bestKey, bestDevice, nil
}

// tryFreeMemory attempts to free memory by evicting low-priority allocations
func (mp *MemoryPool) tryFreeMemory(neededGB float64) bool {
	var freedGB float64
	var candidates []*Allocation
	
	// Collect eviction candidates
	for _, device := range mp.pools {
		device.mu.RLock()
		for _, allocation := range device.allocations {
			if allocation.Priority <= PriorityLow {
				candidates = append(candidates, allocation)
			}
		}
		device.mu.RUnlock()
	}
	
	// Sort by priority and age (evict oldest low-priority first)
	for i := 0; i < len(candidates); i++ {
		for j := i + 1; j < len(candidates); j++ {
			if candidates[i].Priority > candidates[j].Priority ||
				(candidates[i].Priority == candidates[j].Priority && 
				 candidates[i].CreatedAt.After(candidates[j].CreatedAt)) {
				candidates[i], candidates[j] = candidates[j], candidates[i]
			}
		}
	}
	
	// Evict allocations until we have enough memory
	for _, allocation := range candidates {
		if freedGB >= neededGB {
			break
		}
		
		if err := mp.FreeMemory(allocation.ID); err == nil {
			freedGB += allocation.SizeGB
			mp.metrics.mu.Lock()
			mp.metrics.EvictionCount++
			mp.metrics.mu.Unlock()
		}
	}
	
	return freedGB >= neededGB
}

// FreeMemory releases allocated GPU memory
func (mp *MemoryPool) FreeMemory(allocationID string) error {
	mp.mu.Lock()
	defer mp.mu.Unlock()
	
	// Find allocation across all devices
	for _, device := range mp.pools {
		device.mu.Lock()
		if allocation, exists := device.allocations[allocationID]; exists {
			delete(device.allocations, allocationID)
			device.availableGB += allocation.SizeGB
			device.mu.Unlock()
			
			mp.updateMetrics(allocation, false)
			log.Printf("Freed %.2fGB (ID: %s)", allocation.SizeGB, allocationID)
			return nil
		}
		device.mu.Unlock()
	}
	
	return fmt.Errorf("allocation %s not found", allocationID)
}

// optimizationLoop runs background optimization tasks
func (mp *MemoryPool) optimizationLoop() {
	ticker := time.NewTicker(mp.config.OptimizationInterval)
	defer ticker.Stop()
	
	for range ticker.C {
		mp.runOptimization()
	}
}

// runOptimization performs adaptive memory optimization
func (mp *MemoryPool) runOptimization() {
	mp.mu.Lock()
	defer mp.mu.Unlock()
	
	// Clean up expired cache entries
	if mp.config.CacheEnabled {
		mp.cleanupCache()
	}
	
	// Rebalance allocations across devices if needed
	mp.rebalanceAllocations()
	
	// Update optimization timestamp
	mp.metrics.mu.Lock()
	mp.metrics.LastOptimization = time.Now()
	mp.metrics.mu.Unlock()
	
	log.Printf("GPU memory optimization completed")
}

// cleanupCache removes expired cache entries
func (mp *MemoryPool) cleanupCache() {
	now := time.Now()
	
	for _, device := range mp.pools {
		device.mu.Lock()
		for id, allocation := range device.allocations {
			if allocation.Purpose == "cache" && 
			   now.Sub(allocation.LastAccessed) > mp.config.CacheTTL {
				delete(device.allocations, id)
				device.availableGB += allocation.SizeGB
			}
		}
		device.mu.Unlock()
	}
}

// rebalanceAllocations redistributes allocations for optimal performance
func (mp *MemoryPool) rebalanceAllocations() {
	// Calculate load imbalance
	var totalUsed, totalAvailable float64
	deviceLoads := make(map[string]float64)
	
	for key, device := range mp.pools {
		device.mu.RLock()
		used := device.totalMemoryGB - device.availableGB
		deviceLoads[key] = used / device.totalMemoryGB
		totalUsed += used
		totalAvailable += device.totalMemoryGB
		device.mu.RUnlock()
	}
	
	avgLoad := totalUsed / totalAvailable
	
	// If load imbalance > 20%, consider rebalancing
	for key, load := range deviceLoads {
		if load > avgLoad*1.2 {
			log.Printf("Device %s overloaded (%.1f%% vs %.1f%% avg)", key, load*100, avgLoad*100)
			// In a real implementation, this would trigger allocation migration
		}
	}
}

// updateMetrics updates memory usage metrics
func (mp *MemoryPool) updateMetrics(allocation *Allocation, isAllocation bool) {
	mp.metrics.mu.Lock()
	defer mp.metrics.mu.Unlock()
	
	if isAllocation {
		mp.metrics.TotalAllocationsGB += allocation.SizeGB
		mp.metrics.AllocationCount++
		
		if mp.metrics.TotalAllocationsGB > mp.metrics.PeakUsageGB {
			mp.metrics.PeakUsageGB = mp.metrics.TotalAllocationsGB
		}
	} else {
		mp.metrics.TotalAllocationsGB -= allocation.SizeGB
	}
}

// GetMetrics returns current memory usage metrics
func (mp *MemoryPool) GetMetrics() MemoryMetrics {
	mp.metrics.mu.RLock()
	defer mp.metrics.mu.RUnlock()
	return *mp.metrics
}

// GetDeviceStatus returns status of all GPU devices
func (mp *MemoryPool) GetDeviceStatus() map[string]DeviceStatus {
	mp.mu.RLock()
	defer mp.mu.RUnlock()
	
	status := make(map[string]DeviceStatus)
	
	for key, device := range mp.pools {
		device.mu.RLock()
		
		utilizationPercent := ((device.totalMemoryGB - device.availableGB) / device.totalMemoryGB) * 100
		
		status[key] = DeviceStatus{
			DeviceID:           device.deviceID,
			TotalMemoryGB:      device.totalMemoryGB,
			AvailableMemoryGB:  device.availableGB,
			UtilizationPercent: utilizationPercent,
			AllocationCount:    len(device.allocations),
			LastUsed:           device.lastUsed,
		}
		
		device.mu.RUnlock()
	}
	
	return status
}

// DeviceStatus represents the status of a GPU device
type DeviceStatus struct {
	DeviceID           int
	TotalMemoryGB      float64
	AvailableMemoryGB  float64
	UtilizationPercent float64
	AllocationCount    int
	LastUsed           time.Time
}

// generateAllocationID generates a unique allocation identifier
func generateAllocationID() string {
	return fmt.Sprintf("alloc-%d", time.Now().UnixNano())
}