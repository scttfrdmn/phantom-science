//go:build darwin

package monitoring

import (
	"encoding/binary"
	"fmt"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
	"time"

	"golang.org/x/sys/unix"
)

// getCPUUsage reads kern.cp_time via sysctl twice (100 ms apart) and returns
// the CPU busy percentage.  kern.cp_time holds five uint64 tick counters:
// user, nice, sys, intr, idle (CP_USER..CP_IDLE).
func getCPUUsage() float64 {
	type snapshot struct {
		total uint64
		idle  uint64
	}
	read := func() (snapshot, error) {
		raw, err := unix.SysctlRaw("kern.cp_time")
		if err != nil {
			return snapshot{}, err
		}
		// 5 × 8 bytes (uint64) on 64-bit Darwin
		if len(raw) < 40 {
			return snapshot{}, fmt.Errorf("kern.cp_time: unexpected size %d", len(raw))
		}
		var total uint64
		for i := 0; i < 5; i++ {
			total += binary.LittleEndian.Uint64(raw[i*8:])
		}
		idle := binary.LittleEndian.Uint64(raw[4*8:]) // CP_IDLE is index 4
		return snapshot{total: total, idle: idle}, nil
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

// getMemoryUsage estimates system memory usage.
// Total physical RAM comes from hw.memsize (uint64 sysctl).
// Available memory is approximated from the vm_stat command output
// (free pages + inactive pages × page size).
func getMemoryUsage() float64 {
	// Total physical RAM
	rawMem, err := unix.SysctlRaw("hw.memsize")
	if err != nil || len(rawMem) < 8 {
		// Fall back to process-level heap ratio
		var ms runtime.MemStats
		runtime.ReadMemStats(&ms)
		if ms.Sys > 0 {
			return float64(ms.HeapInuse) / float64(ms.Sys) * 100.0
		}
		return 0.0
	}
	totalMem := binary.LittleEndian.Uint64(rawMem)

	// Available memory via vm_stat
	out, err := exec.Command("vm_stat").Output()
	if err != nil {
		var ms runtime.MemStats
		runtime.ReadMemStats(&ms)
		if totalMem > 0 {
			return float64(ms.Sys) / float64(totalMem) * 100.0
		}
		return 0.0
	}

	var pageSize uint64 = 16384 // default macOS page size (M-series)
	var freePages, inactivePages uint64
	for _, line := range strings.Split(string(out), "\n") {
		if strings.HasPrefix(line, "Mach Virtual Memory Statistics:") {
			var ps uint64
			if n, _ := fmt.Sscanf(line, "Mach Virtual Memory Statistics: (page size of %d bytes)", &ps); n == 1 {
				pageSize = ps
			}
			continue
		}
		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}
		key := strings.TrimSpace(parts[0])
		val, _ := strconv.ParseUint(strings.TrimSpace(strings.TrimRight(parts[1], ".")), 10, 64)
		switch key {
		case "Pages free":
			freePages = val
		case "Pages inactive":
			inactivePages = val
		}
	}

	availMem := (freePages + inactivePages) * pageSize
	if availMem >= totalMem {
		return 0.0
	}
	return float64(totalMem-availMem) / float64(totalMem) * 100.0
}

// getDiskUsage returns the percentage of the root filesystem that is in use.
func getDiskUsage() float64 {
	var stat unix.Statfs_t
	if err := unix.Statfs("/", &stat); err != nil {
		return 0.0
	}
	total := uint64(stat.Blocks) * uint64(stat.Bsize)
	free := uint64(stat.Bfree) * uint64(stat.Bsize)
	if total == 0 {
		return 0.0
	}
	return float64(total-free) / float64(total) * 100.0
}

// getLoadAverage returns the 1-minute load average via `sysctl -n vm.loadavg`.
// The kernel returns a string like "{ 1.23 4.56 7.89 }".
func getLoadAverage() float64 {
	out, err := exec.Command("sysctl", "-n", "vm.loadavg").Output()
	if err != nil {
		return 0.0
	}
	s := strings.Trim(strings.TrimSpace(string(out)), "{} ")
	fields := strings.Fields(s)
	if len(fields) < 1 {
		return 0.0
	}
	load, err := strconv.ParseFloat(fields[0], 64)
	if err != nil {
		return 0.0
	}
	return load
}

// GPU metrics on macOS require Metal/CGO bindings not available in pure Go.
// These return zero so Prometheus gauges reflect actual absence of data.

func getGPUCount() int                    { return 0 }
func getGPUUtilization(_ int) float64     { return 0.0 }
func getGPUMemoryUsage(_ int) float64     { return 0.0 }
func getGPUTemperature(_ int) float64     { return 0.0 }
func getGPUPowerConsumption(_ int) float64 { return 0.0 }
