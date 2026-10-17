package monitoring

import (
	"context"
	"testing"
	"time"
)

func TestDefaultMonitoringConfig(t *testing.T) {
	config := DefaultMonitoringConfig()
	
	if config == nil {
		t.Fatal("DefaultMonitoringConfig should not return nil")
	}
	
	if !config.MetricsEnabled {
		t.Error("Metrics should be enabled by default")
	}
	
	if config.MetricsInterval != 15*time.Second {
		t.Errorf("Expected metrics interval 15s, got %v", config.MetricsInterval)
	}
	
	if config.PrometheusPort != 9090 {
		t.Errorf("Expected Prometheus port 9090, got %d", config.PrometheusPort)
	}
	
	if !config.AlertingEnabled {
		t.Error("Alerting should be enabled by default")
	}
	
	if !config.ComplianceLogging {
		t.Error("Compliance logging should be enabled by default")
	}
	
	if config.Thresholds.CPUUsagePercent != 85.0 {
		t.Errorf("Expected CPU threshold 85.0, got %.1f", config.Thresholds.CPUUsagePercent)
	}
}

func TestNewEnterpriseMonitor(t *testing.T) {
	config := DefaultMonitoringConfig()
	monitor, err := NewEnterpriseMonitor(config)
	
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor should not error: %v", err)
	}
	
	if monitor == nil {
		t.Fatal("NewEnterpriseMonitor should not return nil")
	}
	
	if monitor.config != config {
		t.Error("Monitor should use provided config")
	}
	
	if monitor.metrics == nil {
		t.Error("Monitor should initialize metrics")
	}
	
	if monitor.alerts == nil {
		t.Error("Monitor should initialize alerts when alerting enabled")
	}
	
	if monitor.healthChecks == nil {
		t.Error("Monitor should initialize health checks map")
	}
	
	if monitor.dashboards == nil {
		t.Error("Monitor should initialize dashboards map")
	}
}

func TestNewEnterpriseMonitorNilConfig(t *testing.T) {
	monitor, err := NewEnterpriseMonitor(nil)
	
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor should not error with nil config: %v", err)
	}
	
	if monitor == nil {
		t.Fatal("NewEnterpriseMonitor should not return nil")
	}
	
	if monitor.config == nil {
		t.Error("Monitor should use default config when nil provided")
	}
}

func TestNewEnterpriseMonitorNoAlerting(t *testing.T) {
	config := DefaultMonitoringConfig()
	config.AlertingEnabled = false
	
	monitor, err := NewEnterpriseMonitor(config)
	
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor should not error: %v", err)
	}
	
	if monitor.alerts != nil {
		t.Error("Monitor should not initialize alerts when alerting disabled")
	}
}

func TestEnterpriseMonitorStartStop(t *testing.T) {
	// Use non-standard port to avoid conflicts
	config := DefaultMonitoringConfig()
	config.PrometheusPort = 19090
	
	monitor, err := NewEnterpriseMonitor(config)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	// Test that monitor is not running initially
	if monitor.running {
		t.Error("Monitor should not be running initially")
	}
	
	// Start monitor
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()
	
	err = monitor.Start(ctx)
	if err != nil {
		t.Errorf("Start should not error: %v", err)
	}
	
	if !monitor.running {
		t.Error("Monitor should be running after Start")
	}
	
	// Test starting already running monitor
	err = monitor.Start(ctx)
	if err == nil {
		t.Error("Start should error when already running")
	}
	
	// Stop monitor
	err = monitor.Stop()
	if err != nil {
		t.Errorf("Stop should not error: %v", err)
	}
	
	if monitor.running {
		t.Error("Monitor should not be running after Stop")
	}
	
	// Test stopping already stopped monitor
	err = monitor.Stop()
	if err != nil {
		t.Errorf("Stop should not error when already stopped: %v", err)
	}
}

func TestRegisterHealthCheck(t *testing.T) {
	monitor, err := NewEnterpriseMonitor(nil)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	checkCalled := false
	testCheck := func(ctx context.Context) error {
		checkCalled = true
		return nil
	}
	
	monitor.RegisterHealthCheck("test-check", testCheck)
	
	if len(monitor.healthChecks) != 1 {
		t.Errorf("Expected 1 health check, got %d", len(monitor.healthChecks))
	}
	
	check, exists := monitor.healthChecks["test-check"]
	if !exists {
		t.Error("Health check should be registered")
	}
	
	// Test that the check function works
	err = check(context.Background())
	if err != nil {
		t.Errorf("Health check should not error: %v", err)
	}
	
	if !checkCalled {
		t.Error("Health check should have been called")
	}
}

func TestPerformHealthChecks(t *testing.T) {
	config := DefaultMonitoringConfig()
	config.HealthCheckTimeout = 100 * time.Millisecond
	
	monitor, err := NewEnterpriseMonitor(config)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	successCheckCalled := false
	failureCheckCalled := false
	
	// Register successful health check
	monitor.RegisterHealthCheck("success-check", func(ctx context.Context) error {
		successCheckCalled = true
		return nil
	})
	
	// Register failing health check
	monitor.RegisterHealthCheck("failure-check", func(ctx context.Context) error {
		failureCheckCalled = true
		return context.DeadlineExceeded
	})
	
	// Perform health checks
	ctx := context.Background()
	monitor.performHealthChecks(ctx)
	
	// Give goroutines time to execute
	time.Sleep(200 * time.Millisecond)
	
	if !successCheckCalled {
		t.Error("Success health check should have been called")
	}
	
	if !failureCheckCalled {
		t.Error("Failure health check should have been called")
	}
}

func TestCollectMetrics(t *testing.T) {
	monitor, err := NewEnterpriseMonitor(nil)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	// Should not panic
	monitor.collectMetrics()
	
	// Verify metrics were set (check that they're not zero/default)
	// Note: Using placeholder implementations so values are fixed
	if monitor.metrics.SystemMetrics.CPUUsagePercent == nil {
		t.Error("CPU usage metric should be set")
	}
	
	if monitor.metrics.SystemMetrics.MemoryUsagePercent == nil {
		t.Error("Memory usage metric should be set")
	}
}

func TestSystemMetricsStructure(t *testing.T) {
	monitor, err := NewEnterpriseMonitor(nil)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	metrics := monitor.metrics.SystemMetrics
	
	if metrics.CPUUsagePercent == nil {
		t.Error("CPUUsagePercent should be initialized")
	}
	
	if metrics.MemoryUsagePercent == nil {
		t.Error("MemoryUsagePercent should be initialized")
	}
	
	if metrics.DiskUsagePercent == nil {
		t.Error("DiskUsagePercent should be initialized")
	}
	
	if metrics.NetworkBytesIn == nil {
		t.Error("NetworkBytesIn should be initialized")
	}
	
	if metrics.NetworkBytesOut == nil {
		t.Error("NetworkBytesOut should be initialized")
	}
	
	if metrics.LoadAverage == nil {
		t.Error("LoadAverage should be initialized")
	}
	
	if metrics.OpenFileDescriptors == nil {
		t.Error("OpenFileDescriptors should be initialized")
	}
}

func TestApplicationMetricsStructure(t *testing.T) {
	monitor, err := NewEnterpriseMonitor(nil)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	metrics := monitor.metrics.AppMetrics
	
	if metrics.RequestsTotal == nil {
		t.Error("RequestsTotal should be initialized")
	}
	
	if metrics.RequestDuration == nil {
		t.Error("RequestDuration should be initialized")
	}
	
	if metrics.ResponseSize == nil {
		t.Error("ResponseSize should be initialized")
	}
	
	if metrics.ErrorsTotal == nil {
		t.Error("ErrorsTotal should be initialized")
	}
	
	if metrics.ActiveConnections == nil {
		t.Error("ActiveConnections should be initialized")
	}
	
	if metrics.QueueDepth == nil {
		t.Error("QueueDepth should be initialized")
	}
	
	if metrics.ProcessingLatency == nil {
		t.Error("ProcessingLatency should be initialized")
	}
}

func TestGPUMetricsStructure(t *testing.T) {
	monitor, err := NewEnterpriseMonitor(nil)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	metrics := monitor.metrics.GPUMetrics
	
	// Test that GaugeVec metrics can be used
	metrics.GPUUtilization.WithLabelValues("0").Set(50.0)
	metrics.GPUMemoryUsage.WithLabelValues("0").Set(75.0)
	metrics.GPUTemperature.WithLabelValues("0").Set(65.0)
	metrics.GPUPowerConsumption.WithLabelValues("0").Set(200.0)
	
	if metrics.CUDAErrors == nil {
		t.Error("CUDAErrors should be initialized")
	}
	
	if metrics.ProcessingThroughput == nil {
		t.Error("ProcessingThroughput should be initialized")
	}
}

func TestBusinessMetricsStructure(t *testing.T) {
	monitor, err := NewEnterpriseMonitor(nil)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	metrics := monitor.metrics.BusinessMetrics
	
	if metrics.CustomersActive == nil {
		t.Error("CustomersActive should be initialized")
	}
	
	if metrics.DatasetsProcessed == nil {
		t.Error("DatasetsProcessed should be initialized")
	}
	
	if metrics.RevenueGenerated == nil {
		t.Error("RevenueGenerated should be initialized")
	}
	
	if metrics.SLACompliance == nil {
		t.Error("SLACompliance should be initialized")
	}
	
	if metrics.CustomerSatisfaction == nil {
		t.Error("CustomerSatisfaction should be initialized")
	}
}

func TestComplianceMetricsStructure(t *testing.T) {
	monitor, err := NewEnterpriseMonitor(nil)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	metrics := monitor.metrics.ComplianceMetrics
	
	if metrics.AuditEventsTotal == nil {
		t.Error("AuditEventsTotal should be initialized")
	}
	
	if metrics.ComplianceViolations == nil {
		t.Error("ComplianceViolations should be initialized")
	}
	
	if metrics.DataRetentionDays == nil {
		t.Error("DataRetentionDays should be initialized")
	}
	
	if metrics.SecurityIncidents == nil {
		t.Error("SecurityIncidents should be initialized")
	}
	
	if metrics.ValidationFailures == nil {
		t.Error("ValidationFailures should be initialized")
	}
}

func TestPrometheusMetricsStructure(t *testing.T) {
	monitor, err := NewEnterpriseMonitor(nil)
	if err != nil {
		t.Fatalf("NewEnterpriseMonitor failed: %v", err)
	}
	
	promMetrics := monitor.metrics.PrometheusMetrics
	
	if promMetrics == nil {
		t.Fatal("PrometheusMetrics should be initialized")
	}
	
	if promMetrics.Registry == nil {
		t.Error("Prometheus Registry should be initialized")
	}
}

func TestMonitoringConfigStructure(t *testing.T) {
	config := &MonitoringConfig{
		MetricsEnabled:       true,
		MetricsInterval:      30 * time.Second,
		MetricsRetention:     7 * 24 * time.Hour,
		PrometheusPort:       9091,
		HealthCheckInterval:  60 * time.Second,
		HealthCheckTimeout:   5 * time.Second,
		AlertingEnabled:      true,
		AlertManagerURL:      "http://localhost:9093",
		PagerDutyKey:         "test-key",
		SlackWebhookURL:      "https://hooks.slack.com/test",
		ComplianceLogging:    true,
		AuditTrailEnabled:    true,
		RegulatorySLAEnabled: true,
		Thresholds: AlertThresholds{
			CPUUsagePercent:        80.0,
			MemoryUsagePercent:     85.0,
			GPUUsagePercent:        90.0,
			DiskUsagePercent:       80.0,
			ResponseTimeMs:         500.0,
			ErrorRatePercent:       0.5,
			ThroughputMBps:         50.0,
			QueueDepth:             500,
			ConcurrentConnections:  5000,
		},
	}
	
	if !config.MetricsEnabled {
		t.Error("MetricsEnabled should be settable")
	}
	
	if config.MetricsInterval != 30*time.Second {
		t.Error("MetricsInterval should be settable")
	}
	
	if config.PrometheusPort != 9091 {
		t.Error("PrometheusPort should be settable")
	}
	
	if config.Thresholds.CPUUsagePercent != 80.0 {
		t.Error("Threshold values should be settable")
	}
}

func TestAlertThresholdsStructure(t *testing.T) {
	thresholds := AlertThresholds{
		CPUUsagePercent:        75.0,
		MemoryUsagePercent:     80.0,
		GPUUsagePercent:        85.0,
		DiskUsagePercent:       70.0,
		ResponseTimeMs:         2000.0,
		ErrorRatePercent:       2.0,
		ThroughputMBps:         25.0,
		QueueDepth:             2000,
		ConcurrentConnections:  15000,
	}
	
	if thresholds.CPUUsagePercent != 75.0 {
		t.Error("AlertThresholds fields should be settable")
	}
	
	if thresholds.ConcurrentConnections != 15000 {
		t.Error("AlertThresholds integer fields should be settable")
	}
}

func TestDashboardStructure(t *testing.T) {
	dashboard := &Dashboard{
		Name:        "System Overview",
		Description: "System performance dashboard",
		Panels: []DashboardPanel{
			{
				Title: "CPU Usage",
				Type:  "graph",
				Query: "phantom_cpu_usage_percent",
				Thresholds: map[string]float64{
					"warning":  75.0,
					"critical": 90.0,
				},
			},
		},
		LastUpdate: time.Now(),
	}
	
	if dashboard.Name != "System Overview" {
		t.Error("Dashboard name should be settable")
	}
	
	if len(dashboard.Panels) != 1 {
		t.Error("Dashboard panels should be settable")
	}
	
	panel := dashboard.Panels[0]
	if panel.Title != "CPU Usage" {
		t.Error("Panel title should be settable")
	}
	
	if panel.Thresholds["warning"] != 75.0 {
		t.Error("Panel thresholds should be settable")
	}
}

func TestHelperFunctions(t *testing.T) {
	// Test placeholder helper functions
	cpuUsage := getCPUUsage()
	if cpuUsage < 0 || cpuUsage > 100 {
		t.Errorf("CPU usage should be between 0-100, got %.1f", cpuUsage)
	}
	
	memUsage := getMemoryUsage()
	if memUsage < 0 || memUsage > 100 {
		t.Errorf("Memory usage should be between 0-100, got %.1f", memUsage)
	}
	
	diskUsage := getDiskUsage()
	if diskUsage < 0 || diskUsage > 100 {
		t.Errorf("Disk usage should be between 0-100, got %.1f", diskUsage)
	}
	
	loadAvg := getLoadAverage()
	if loadAvg < 0 {
		t.Errorf("Load average should not be negative, got %.1f", loadAvg)
	}
	
	gpuCount := getGPUCount()
	if gpuCount < 0 {
		t.Errorf("GPU count should not be negative, got %d", gpuCount)
	}
	
	if gpuCount > 0 {
		gpuUtil := getGPUUtilization(0)
		if gpuUtil < 0 || gpuUtil > 100 {
			t.Errorf("GPU utilization should be between 0-100, got %.1f", gpuUtil)
		}
		
		gpuMem := getGPUMemoryUsage(0)
		if gpuMem < 0 || gpuMem > 100 {
			t.Errorf("GPU memory usage should be between 0-100, got %.1f", gpuMem)
		}
		
		gpuTemp := getGPUTemperature(0)
		if gpuTemp < 0 {
			t.Errorf("GPU temperature should not be negative, got %.1f", gpuTemp)
		}
		
		gpuPower := getGPUPowerConsumption(0)
		if gpuPower < 0 {
			t.Errorf("GPU power consumption should not be negative, got %.1f", gpuPower)
		}
	}
}