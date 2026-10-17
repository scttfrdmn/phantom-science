//go:build linux

package monitoring

import (
	"fmt"
	"os"
	"strconv"
	"strings"
	"syscall"
	"time"
)

// getCPUUsage reads /proc/stat twice (100 ms apart) and returns the CPU busy
// percentage across all cores.
func getCPUUsage() float64 {
	type snapshot struct {
		total uint64
		idle  uint64
	}
	read := func() (snapshot, error) {
		data, err := os.ReadFile("/proc/stat")
		if err != nil {
			return snapshot{}, err
		}
		for _, line := range strings.SplitN(string(data), "\n", 2) {
			if !strings.HasPrefix(line, "cpu ") {
				continue
			}
			var user, nice, system, idle, iowait, irq, softirq, steal uint64
			if _, err := fmt.Sscanf(line, "cpu %d %d %d %d %d %d %d %d",
				&user, &nice, &system, &idle, &iowait, &irq, &softirq, &steal); err != nil {
				return snapshot{}, err
			}
			total := user + nice + system + idle + iowait + irq + softirq + steal
			return snapshot{total: total, idle: idle + iowait}, nil
		}
		return snapshot{}, fmt.Errorf("cpu line not found in /proc/stat")
	}

	s1, err := read()
	if err != nil {
		return 0.0
	}
	time.Sleep(100 * time.Millisecond)
	s2, err := read()
	if err != nil {
		return 0.0
	}

	dTotal := float64(s2.total - s1.total)
	dIdle := float64(s2.idle - s1.idle)
	if dTotal <= 0 {
		return 0.0
	}
	return (dTotal - dIdle) / dTotal * 100.0
}

// getMemoryUsage parses /proc/meminfo and returns the percentage of physical
// RAM that is in active use (total − free − buffers − cached).
func getMemoryUsage() float64 {
	data, err := os.ReadFile("/proc/meminfo")
	if err != nil {
		return 0.0
	}
	var total, free, buffers, cached uint64
	for _, line := range strings.Split(string(data), "\n") {
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}
		val, _ := strconv.ParseUint(fields[1], 10, 64)
		switch fields[0] {
		case "MemTotal:":
			total = val
		case "MemFree:":
			free = val
		case "Buffers:":
			buffers = val
		case "Cached:":
			cached = val
		}
	}
	if total == 0 {
		return 0.0
	}
	used := total - free - buffers - cached
	return float64(used) / float64(total) * 100.0
}

// getDiskUsage returns the percentage of the root filesystem that is in use.
func getDiskUsage() float64 {
	var stat syscall.Statfs_t
	if err := syscall.Statfs("/", &stat); err != nil {
		return 0.0
	}
	total := stat.Blocks * uint64(stat.Bsize)
	free := stat.Bfree * uint64(stat.Bsize)
	if total == 0 {
		return 0.0
	}
	return float64(total-free) / float64(total) * 100.0
}

// getLoadAverage returns the 1-minute load average from /proc/loadavg.
func getLoadAverage() float64 {
	data, err := os.ReadFile("/proc/loadavg")
	if err != nil {
		return 0.0
	}
	fields := strings.Fields(string(data))
	if len(fields) < 1 {
		return 0.0
	}
	load, err := strconv.ParseFloat(fields[0], 64)
	if err != nil {
		return 0.0
	}
	return load
}

// GPU metrics require NVML/CGO bindings not available in pure Go.
// The CUDA acceleration package exposes GPU metrics when built with -tags cuda.
// These functions return zero so Prometheus gauges remain at 0 instead of a
// hard-coded placeholder value.

func getGPUCount() int        { return 0 }
func getGPUUtilization(_ int) float64    { return 0.0 }
func getGPUMemoryUsage(_ int) float64    { return 0.0 }
func getGPUTemperature(_ int) float64    { return 0.0 }
func getGPUPowerConsumption(_ int) float64 { return 0.0 }
