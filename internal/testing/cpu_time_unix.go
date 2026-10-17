//go:build unix

package testing

import "syscall"

// cpuSecondsSinceStart returns the cumulative CPU time (user + system) in
// seconds consumed by this process, or 0 if it cannot be read.
func cpuSecondsSinceStart() float64 {
	var ru syscall.Rusage
	if err := syscall.Getrusage(syscall.RUSAGE_SELF, &ru); err != nil {
		return 0
	}
	return float64(ru.Utime.Nano()+ru.Stime.Nano()) / 1e9
}
