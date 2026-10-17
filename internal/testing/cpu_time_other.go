//go:build !unix

package testing

// cpuSecondsSinceStart returns 0: process CPU time is not measured on this
// platform.
func cpuSecondsSinceStart() float64 {
	return 0
}
