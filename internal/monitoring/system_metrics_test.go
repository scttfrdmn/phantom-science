//go:build !release

package monitoring

import (
	"testing"
)

// TestGetCPUUsage verifies that getCPUUsage returns a value in [0, 100].
func TestGetCPUUsage(t *testing.T) {
	v := getCPUUsage()
	if v < 0 || v > 100 {
		t.Errorf("getCPUUsage() = %f, want value in [0, 100]", v)
	}
}

// TestGetMemoryUsage verifies that getMemoryUsage returns a value in [0, 100].
func TestGetMemoryUsage(t *testing.T) {
	v := getMemoryUsage()
	if v < 0 || v > 100 {
		t.Errorf("getMemoryUsage() = %f, want value in [0, 100]", v)
	}
}

// TestGetDiskUsage verifies that getDiskUsage returns a value in [0, 100].
// On any supported OS the root filesystem will have non-zero total size.
func TestGetDiskUsage(t *testing.T) {
	v := getDiskUsage()
	if v < 0 || v > 100 {
		t.Errorf("getDiskUsage() = %f, want value in [0, 100]", v)
	}
	// Root filesystem should always have some usage.
	if v == 0 {
		t.Log("getDiskUsage() returned 0 — may indicate an unsupported platform or root-mount edge case")
	}
}

// TestGetLoadAverage verifies that getLoadAverage returns a non-negative value.
func TestGetLoadAverage(t *testing.T) {
	v := getLoadAverage()
	if v < 0 {
		t.Errorf("getLoadAverage() = %f, want >= 0", v)
	}
}

// TestGetGPUCount verifies that getGPUCount returns a non-negative value.
func TestGetGPUCount(t *testing.T) {
	n := getGPUCount()
	if n < 0 {
		t.Errorf("getGPUCount() = %d, want >= 0", n)
	}
}

// TestGPUMetrics verifies that GPU metric functions return values in plausible ranges.
func TestGPUMetrics(t *testing.T) {
	for _, fn := range []struct {
		name string
		f    func(int) float64
	}{
		{"getGPUUtilization", getGPUUtilization},
		{"getGPUMemoryUsage", getGPUMemoryUsage},
		{"getGPUTemperature", getGPUTemperature},
		{"getGPUPowerConsumption", getGPUPowerConsumption},
	} {
		v := fn.f(0)
		if v < 0 {
			t.Errorf("%s(0) = %f, want >= 0", fn.name, v)
		}
	}
}

// TestCollectMetrics_NoZeroPanic verifies that collectMetrics() can run without
// panicking, even when GPU count is 0.
func TestCollectMetrics_NoZeroPanic(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Errorf("collectMetrics panicked: %v", r)
		}
	}()

	m, err := NewEnterpriseMonitor(DefaultMonitoringConfig())
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor: %v", err)
	}
	m.collectMetrics()
}
