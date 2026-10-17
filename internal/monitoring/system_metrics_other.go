//go:build !linux && !darwin

// Package monitoring — system metric stubs for unsupported platforms.
// Real implementations live in system_metrics_linux.go and system_metrics_darwin.go.
package monitoring

func getCPUUsage() float64         { return 0.0 }
func getMemoryUsage() float64      { return 0.0 }
func getDiskUsage() float64        { return 0.0 }
func getLoadAverage() float64      { return 0.0 }
func getGPUCount() int             { return 0 }
func getGPUUtilization(_ int) float64      { return 0.0 }
func getGPUMemoryUsage(_ int) float64      { return 0.0 }
func getGPUTemperature(_ int) float64      { return 0.0 }
func getGPUPowerConsumption(_ int) float64 { return 0.0 }
