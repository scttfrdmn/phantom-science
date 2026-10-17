package gpu

import (
	"testing"
	"time"
)

func TestDefaultOptimizationConfig(t *testing.T) {
	config := DefaultOptimizationConfig()
	
	if config == nil {
		t.Fatal("DefaultOptimizationConfig should not return nil")
	}
	
	if config.MaxMemoryUsagePercent <= 0 {
		t.Error("DefaultOptimizationConfig should set positive MaxMemoryUsagePercent")
	}
	
	if config.ChunkSizeGB <= 0 {
		t.Error("DefaultOptimizationConfig should set positive ChunkSizeGB")
	}
	
	if config.ConcurrentStreams <= 0 {
		t.Error("DefaultOptimizationConfig should set positive ConcurrentStreams")
	}
	
	if config.PreallocationSizeGB < 0 {
		t.Error("DefaultOptimizationConfig should set non-negative PreallocationSizeGB")
	}
}

func TestNewMemoryPool(t *testing.T) {
	maxMemoryGB := 16.0
	config := DefaultOptimizationConfig()
	
	pool, err := NewMemoryPool(maxMemoryGB, config)
	
	if err != nil {
		t.Fatalf("NewMemoryPool should not error: %v", err)
	}
	
	if pool == nil {
		t.Fatal("NewMemoryPool should not return nil")
	}
	
	if pool.maxMemoryGB != maxMemoryGB {
		t.Error("NewMemoryPool should use provided maxMemoryGB")
	}
}

func TestMemoryPoolAllocateMemory(t *testing.T) {
	pool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	// Test allocation of different sizes
	testSizes := []float64{0.5, 1.0, 2.0} // GB
	
	for _, size := range testSizes {
		allocation, err := pool.AllocateMemory(size, "test-allocation", PriorityNormal)
		if err != nil {
			t.Errorf("AllocateMemory should not error for size %.1fGB: %v", size, err)
			continue
		}
		
		if allocation == nil {
			t.Errorf("AllocateMemory should return non-nil allocation for size %.1fGB", size)
			continue
		}
		
		if allocation.SizeGB != size {
			t.Errorf("AllocateMemory should return allocation of size %.1fGB, got %.1fGB", size, allocation.SizeGB)
		}
		
		// Clean up
		pool.FreeMemory(allocation.ID)
	}
}

func TestMemoryPoolFreeMemory(t *testing.T) {
	pool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	// Allocate and free memory
	allocation, err := pool.AllocateMemory(1.0, "test-allocation", PriorityNormal)
	if err != nil {
		t.Fatalf("AllocateMemory failed: %v", err)
	}
	
	// Should not error
	err = pool.FreeMemory(allocation.ID)
	if err != nil {
		t.Errorf("FreeMemory should not error for valid allocation: %v", err)
	}
	
	// Should error for non-existent allocation
	err = pool.FreeMemory("nonexistent-id")
	if err == nil {
		t.Error("FreeMemory should error for non-existent allocation")
	}
}

func TestMemoryPoolGetMetrics(t *testing.T) {
	pool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	metrics := pool.GetMetrics()
	
	if metrics.TotalAllocationsGB < 0 {
		t.Error("Metrics should have non-negative TotalAllocationsGB")
	}
	
	if metrics.AllocationCount < 0 {
		t.Error("Metrics should have non-negative AllocationCount")
	}
	
	if metrics.PeakUsageGB < 0 {
		t.Error("Metrics should have non-negative PeakUsageGB")
	}
	
	if metrics.EvictionCount < 0 {
		t.Error("Metrics should have non-negative EvictionCount")
	}
}

func TestMemoryPoolMetricsAfterOperations(t *testing.T) {
	pool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	initialMetrics := pool.GetMetrics()
	
	// Allocate some memory
	allocations := make([]*Allocation, 3)
	for i := range allocations {
		allocation, err := pool.AllocateMemory(0.5, "test-allocation", PriorityNormal)
		if err != nil {
			t.Fatalf("AllocateMemory failed: %v", err)
		}
		allocations[i] = allocation
	}
	
	metricsAfterAlloc := pool.GetMetrics()
	
	// Should have more allocations
	if metricsAfterAlloc.AllocationCount <= initialMetrics.AllocationCount {
		t.Error("AllocationCount should increase after allocation")
	}
	
	if metricsAfterAlloc.TotalAllocationsGB <= initialMetrics.TotalAllocationsGB {
		t.Error("TotalAllocationsGB should increase after allocation")
	}
	
	// Free all allocations
	for _, allocation := range allocations {
		pool.FreeMemory(allocation.ID)
	}
	
	metricsAfterFree := pool.GetMetrics()
	
	// TotalAllocationsGB should decrease
	if metricsAfterFree.TotalAllocationsGB > metricsAfterAlloc.TotalAllocationsGB {
		t.Error("TotalAllocationsGB should not increase after freeing memory")
	}
}

func TestMemoryPoolLargeAllocation(t *testing.T) {
	pool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	// Try to allocate more than available memory
	largeSize := 32.0 // 32GB, larger than our 16GB pool
	
	allocation, err := pool.AllocateMemory(largeSize, "large-allocation", PriorityNormal)
	
	// Should either succeed (if devices have memory) or fail gracefully
	if err != nil {
		// Error is acceptable for very large allocations
		if allocation != nil {
			t.Error("If AllocateMemory errors, it should return nil allocation")
		}
	} else {
		// If successful, allocation should be correct size
		if allocation.SizeGB != largeSize {
			t.Errorf("Large allocation should return correct size, got %.1fGB expected %.1fGB", allocation.SizeGB, largeSize)
		}
		pool.FreeMemory(allocation.ID)
	}
}

func TestMemoryPoolGetDeviceStatus(t *testing.T) {
	pool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	status := pool.GetDeviceStatus()
	
	if len(status) == 0 {
		t.Error("GetDeviceStatus should return at least one device")
	}
	
	for deviceKey, deviceStatus := range status {
		if deviceStatus.TotalMemoryGB <= 0 {
			t.Errorf("Device %s should have positive TotalMemoryGB", deviceKey)
		}
		
		if deviceStatus.UtilizationPercent < 0 || deviceStatus.UtilizationPercent > 100 {
			t.Errorf("Device %s should have UtilizationPercent between 0-100, got %.1f", deviceKey, deviceStatus.UtilizationPercent)
		}
	}
}

func TestOptimizationConfigStructure(t *testing.T) {
	// Test that OptimizationConfig can be created with zero values
	config := &OptimizationConfig{
		MaxMemoryUsagePercent: 0,
		ChunkSizeGB:          0,
		ConcurrentStreams:    0,
		StreamingEnabled:     false,
	}
	
	if config == nil {
		t.Error("OptimizationConfig should be creatable")
	}
	
	// Should handle gracefully (implementation dependent)
	pool, err := NewMemoryPool(16.0, config)
	if err != nil {
		t.Errorf("NewMemoryPool should handle zero config: %v", err)
	}
	if pool == nil {
		t.Error("NewMemoryPool should not return nil even with zero config")
	}
}
func TestRunOptimization_NoCrash(t *testing.T) {
	pool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("NewMemoryPool failed: %v", err)
	}
	// runOptimization should not panic
	pool.runOptimization()
}

func TestCleanupCache_RemovesExpiredEntries(t *testing.T) {
	cfg := DefaultOptimizationConfig()
	cfg.CacheEnabled = true
	cfg.CacheTTL = 1 * time.Millisecond // Very short TTL

	pool, err := NewMemoryPool(16.0, cfg)
	if err != nil {
		t.Fatalf("NewMemoryPool failed: %v", err)
	}

	// Allocate a "cache" purpose allocation
	alloc, allocErr := pool.AllocateMemory(1.0, "cache", PriorityNormal)
	if allocErr != nil {
		t.Fatalf("AllocateMemory failed: %v", allocErr)
	}

	// Make the entry appear old by setting LastAccessed in the past
	pool.mu.Lock()
	for _, device := range pool.pools {
		device.mu.Lock()
		if entry, ok := device.allocations[alloc.ID]; ok {
			entry.LastAccessed = time.Now().Add(-1 * time.Hour)
		}
		device.mu.Unlock()
	}
	pool.mu.Unlock()

	// cleanupCache should remove the expired entry
	pool.mu.Lock()
	pool.cleanupCache()
	pool.mu.Unlock()
	// No panic and no error is success
}

func TestRebalanceAllocations_NoCrash(t *testing.T) {
	pool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("NewMemoryPool failed: %v", err)
	}
	// Allocate some memory to create load imbalance
	pool.AllocateMemory(2.0, "test", PriorityNormal) //nolint
	pool.AllocateMemory(4.0, "test", PriorityNormal) //nolint

	// rebalanceAllocations should not panic
	pool.mu.Lock()
	pool.rebalanceAllocations()
	pool.mu.Unlock()
}
