// Package monitoring provides enterprise-grade monitoring and alerting for Phantom Science
package monitoring

import (
	"context"
	"fmt"
	"log"
	"net/http"
	"sync"
	"time"

	"github.com/prometheus/client_golang/prometheus"
	"github.com/prometheus/client_golang/prometheus/promhttp"
)

// EnterpriseMonitor provides comprehensive monitoring and alerting
type EnterpriseMonitor struct {
	config         *MonitoringConfig
	metrics        *EnterpriseMetrics
	alerts         *AlertManager
	healthChecks   map[string]HealthCheck
	dashboards     map[string]*Dashboard
	mu             sync.RWMutex
	running        bool
	stopChan       chan struct{}
}

// MonitoringConfig defines enterprise monitoring configuration
type MonitoringConfig struct {
	// Metrics collection
	MetricsEnabled       bool          `yaml:"metrics_enabled"`
	MetricsInterval      time.Duration `yaml:"metrics_interval"`
	MetricsRetention     time.Duration `yaml:"metrics_retention"`
	PrometheusPort       int           `yaml:"prometheus_port"`
	
	// Health checks
	HealthCheckInterval  time.Duration `yaml:"health_check_interval"`
	HealthCheckTimeout   time.Duration `yaml:"health_check_timeout"`
	
	// Alerting
	AlertingEnabled      bool          `yaml:"alerting_enabled"`
	AlertManagerURL      string        `json:"-" yaml:"-"`
	PagerDutyKey         string        `json:"-" yaml:"-"`
	SlackWebhookURL      string        `json:"-" yaml:"-"`
	
	// Compliance
	ComplianceLogging    bool          `yaml:"compliance_logging"`
	AuditTrailEnabled    bool          `yaml:"audit_trail_enabled"`
	RegulatorySLAEnabled bool          `yaml:"regulatory_sla_enabled"`
	
	// Performance thresholds
	Thresholds           AlertThresholds `yaml:"thresholds"`
}

// AlertThresholds defines performance and health thresholds
type AlertThresholds struct {
	CPUUsagePercent        float64 `yaml:"cpu_usage_percent"`
	MemoryUsagePercent     float64 `yaml:"memory_usage_percent"`
	GPUUsagePercent        float64 `yaml:"gpu_usage_percent"`
	DiskUsagePercent       float64 `yaml:"disk_usage_percent"`
	ResponseTimeMs         float64 `yaml:"response_time_ms"`
	ErrorRatePercent       float64 `yaml:"error_rate_percent"`
	ThroughputMBps         float64 `yaml:"throughput_mbps"`
	QueueDepth             int     `yaml:"queue_depth"`
	ConcurrentConnections  int     `yaml:"concurrent_connections"`
}

// EnterpriseMetrics contains all monitoring metrics
type EnterpriseMetrics struct {
	// System metrics
	SystemMetrics    *SystemMetrics
	
	// Application metrics  
	AppMetrics       *ApplicationMetrics
	
	// GPU metrics
	GPUMetrics       *GPUMetrics
	
	// Business metrics
	BusinessMetrics  *BusinessMetrics
	
	// Compliance metrics
	ComplianceMetrics *ComplianceMetrics
	
	// Prometheus metrics
	PrometheusMetrics *PrometheusMetrics
}

// SystemMetrics tracks system-level performance
type SystemMetrics struct {
	CPUUsagePercent    prometheus.Gauge
	MemoryUsagePercent prometheus.Gauge
	DiskUsagePercent   prometheus.Gauge
	NetworkBytesIn     prometheus.Counter
	NetworkBytesOut    prometheus.Counter
	LoadAverage        prometheus.Gauge
	OpenFileDescriptors prometheus.Gauge
}

// ApplicationMetrics tracks application-specific metrics
type ApplicationMetrics struct {
	RequestsTotal         prometheus.Counter
	RequestDuration       prometheus.Histogram
	ResponseSize          prometheus.Histogram
	ErrorsTotal           prometheus.Counter
	ActiveConnections     prometheus.Gauge
	QueueDepth           prometheus.Gauge
	ProcessingLatency    prometheus.Histogram
}

// GPUMetrics tracks GPU performance and utilization
type GPUMetrics struct {
	GPUUtilization       prometheus.GaugeVec
	GPUMemoryUsage       prometheus.GaugeVec
	GPUTemperature       prometheus.GaugeVec
	GPUPowerConsumption  prometheus.GaugeVec
	CUDAErrors          prometheus.Counter
	ProcessingThroughput prometheus.Gauge
}

// BusinessMetrics tracks business and operational KPIs
type BusinessMetrics struct {
	CustomersActive      prometheus.Gauge
	DatasetsProcessed    prometheus.Counter
	RevenueGenerated     prometheus.Counter
	SLACompliance        prometheus.Gauge
	CustomerSatisfaction prometheus.Gauge
}

// ComplianceMetrics tracks regulatory compliance metrics
type ComplianceMetrics struct {
	AuditEventsTotal     prometheus.Counter
	ComplianceViolations prometheus.Counter
	DataRetentionDays    prometheus.Gauge
	SecurityIncidents    prometheus.Counter
	ValidationFailures   prometheus.Counter
}

// PrometheusMetrics wraps all Prometheus metric collections
type PrometheusMetrics struct {
	Registry *prometheus.Registry
}

// HealthCheck represents a health check function
type HealthCheck func(ctx context.Context) error

// Dashboard represents a monitoring dashboard
type Dashboard struct {
	Name        string
	Description string
	Panels      []DashboardPanel
	LastUpdate  time.Time
}

// DashboardPanel represents a panel in a monitoring dashboard  
type DashboardPanel struct {
	Title       string
	Type        string
	Query       string
	Thresholds  map[string]float64
}

// DefaultMonitoringConfig returns enterprise-grade monitoring configuration
func DefaultMonitoringConfig() *MonitoringConfig {
	return &MonitoringConfig{
		MetricsEnabled:       true,
		MetricsInterval:      15 * time.Second,
		MetricsRetention:     30 * 24 * time.Hour, // 30 days
		PrometheusPort:       9090,
		
		HealthCheckInterval:  30 * time.Second,
		HealthCheckTimeout:   10 * time.Second,
		
		AlertingEnabled:      true,
		ComplianceLogging:    true,
		AuditTrailEnabled:    true,
		RegulatorySLAEnabled: true,
		
		Thresholds: AlertThresholds{
			CPUUsagePercent:       85.0,
			MemoryUsagePercent:    90.0,
			GPUUsagePercent:       95.0,
			DiskUsagePercent:      85.0,
			ResponseTimeMs:        1000.0,
			ErrorRatePercent:      1.0,
			ThroughputMBps:        100.0,
			QueueDepth:            1000,
			ConcurrentConnections: 10000,
		},
	}
}

// NewEnterpriseMonitor creates a new enterprise monitoring system
func NewEnterpriseMonitor(config *MonitoringConfig) (*EnterpriseMonitor, error) {
	if config == nil {
		config = DefaultMonitoringConfig()
	}
	
	monitor := &EnterpriseMonitor{
		config:       config,
		healthChecks: make(map[string]HealthCheck),
		dashboards:   make(map[string]*Dashboard),
		stopChan:     make(chan struct{}),
	}
	
	// Initialize metrics
	if err := monitor.initializeMetrics(); err != nil {
		return nil, fmt.Errorf("failed to initialize metrics: %w", err)
	}
	
	// Initialize alert manager
	if config.AlertingEnabled {
		monitor.alerts = NewAlertManager(config)
	}
	
	return monitor, nil
}

// initializeMetrics sets up all Prometheus metrics
func (em *EnterpriseMonitor) initializeMetrics() error {
	registry := prometheus.NewRegistry()
	
	// System metrics
	systemMetrics := &SystemMetrics{
		CPUUsagePercent: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_cpu_usage_percent",
			Help: "Current CPU usage percentage",
		}),
		MemoryUsagePercent: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_memory_usage_percent", 
			Help: "Current memory usage percentage",
		}),
		DiskUsagePercent: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_disk_usage_percent",
			Help: "Current disk usage percentage",
		}),
		NetworkBytesIn: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_network_bytes_received_total",
			Help: "Total bytes received over network",
		}),
		NetworkBytesOut: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_network_bytes_sent_total",
			Help: "Total bytes sent over network",
		}),
		LoadAverage: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_load_average",
			Help: "System load average",
		}),
		OpenFileDescriptors: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_open_file_descriptors",
			Help: "Number of open file descriptors",
		}),
	}
	
	// Application metrics
	appMetrics := &ApplicationMetrics{
		RequestsTotal: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_requests_total",
			Help: "Total number of requests processed",
		}),
		RequestDuration: prometheus.NewHistogram(prometheus.HistogramOpts{
			Name:    "phantom_request_duration_seconds",
			Help:    "Request duration in seconds",
			Buckets: prometheus.DefBuckets,
		}),
		ResponseSize: prometheus.NewHistogram(prometheus.HistogramOpts{
			Name:    "phantom_response_size_bytes",
			Help:    "Response size in bytes",
			Buckets: prometheus.ExponentialBuckets(1024, 2, 10),
		}),
		ErrorsTotal: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_errors_total",
			Help: "Total number of errors",
		}),
		ActiveConnections: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_active_connections",
			Help: "Number of active connections",
		}),
		QueueDepth: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_queue_depth",
			Help: "Current processing queue depth",
		}),
		ProcessingLatency: prometheus.NewHistogram(prometheus.HistogramOpts{
			Name:    "phantom_processing_latency_seconds",
			Help:    "Processing latency in seconds",
			Buckets: prometheus.ExponentialBuckets(0.1, 2, 10),
		}),
	}
	
	// GPU metrics
	gpuMetrics := &GPUMetrics{
		GPUUtilization: *prometheus.NewGaugeVec(
			prometheus.GaugeOpts{
				Name: "phantom_gpu_utilization_percent",
				Help: "GPU utilization percentage",
			},
			[]string{"device_id"},
		),
		GPUMemoryUsage: *prometheus.NewGaugeVec(
			prometheus.GaugeOpts{
				Name: "phantom_gpu_memory_usage_percent",
				Help: "GPU memory usage percentage",
			},
			[]string{"device_id"},
		),
		GPUTemperature: *prometheus.NewGaugeVec(
			prometheus.GaugeOpts{
				Name: "phantom_gpu_temperature_celsius",
				Help: "GPU temperature in Celsius",
			},
			[]string{"device_id"},
		),
		GPUPowerConsumption: *prometheus.NewGaugeVec(
			prometheus.GaugeOpts{
				Name: "phantom_gpu_power_consumption_watts",
				Help: "GPU power consumption in watts",
			},
			[]string{"device_id"},
		),
		CUDAErrors: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_cuda_errors_total",
			Help: "Total number of CUDA errors",
		}),
		ProcessingThroughput: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_processing_throughput_mbps",
			Help: "Current processing throughput in MB/s",
		}),
	}
	
	// Business metrics
	businessMetrics := &BusinessMetrics{
		CustomersActive: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_customers_active",
			Help: "Number of active customers",
		}),
		DatasetsProcessed: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_datasets_processed_total",
			Help: "Total number of datasets processed",
		}),
		RevenueGenerated: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_revenue_generated_total",
			Help: "Total revenue generated",
		}),
		SLACompliance: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_sla_compliance_percent",
			Help: "SLA compliance percentage",
		}),
		CustomerSatisfaction: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_customer_satisfaction_score",
			Help: "Customer satisfaction score",
		}),
	}
	
	// Compliance metrics
	complianceMetrics := &ComplianceMetrics{
		AuditEventsTotal: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_audit_events_total",
			Help: "Total number of audit events",
		}),
		ComplianceViolations: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_compliance_violations_total",
			Help: "Total number of compliance violations",
		}),
		DataRetentionDays: prometheus.NewGauge(prometheus.GaugeOpts{
			Name: "phantom_data_retention_days",
			Help: "Data retention period in days",
		}),
		SecurityIncidents: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_security_incidents_total",
			Help: "Total number of security incidents",
		}),
		ValidationFailures: prometheus.NewCounter(prometheus.CounterOpts{
			Name: "phantom_validation_failures_total",
			Help: "Total number of validation failures",
		}),
	}
	
	// Register all metrics
	metrics := []*SystemMetrics{systemMetrics}
	for _, m := range metrics {
		registry.MustRegister(
			m.CPUUsagePercent,
			m.MemoryUsagePercent,
			m.DiskUsagePercent,
			m.NetworkBytesIn,
			m.NetworkBytesOut,
			m.LoadAverage,
			m.OpenFileDescriptors,
		)
	}
	
	registry.MustRegister(
		appMetrics.RequestsTotal,
		appMetrics.RequestDuration,
		appMetrics.ResponseSize,
		appMetrics.ErrorsTotal,
		appMetrics.ActiveConnections,
		appMetrics.QueueDepth,
		appMetrics.ProcessingLatency,
		
		&gpuMetrics.GPUUtilization,
		&gpuMetrics.GPUMemoryUsage,
		&gpuMetrics.GPUTemperature,
		&gpuMetrics.GPUPowerConsumption,
		gpuMetrics.CUDAErrors,
		gpuMetrics.ProcessingThroughput,
		
		businessMetrics.CustomersActive,
		businessMetrics.DatasetsProcessed,
		businessMetrics.RevenueGenerated,
		businessMetrics.SLACompliance,
		businessMetrics.CustomerSatisfaction,
		
		complianceMetrics.AuditEventsTotal,
		complianceMetrics.ComplianceViolations,
		complianceMetrics.DataRetentionDays,
		complianceMetrics.SecurityIncidents,
		complianceMetrics.ValidationFailures,
	)
	
	em.metrics = &EnterpriseMetrics{
		SystemMetrics:     systemMetrics,
		AppMetrics:        appMetrics,
		GPUMetrics:        gpuMetrics,
		BusinessMetrics:   businessMetrics,
		ComplianceMetrics: complianceMetrics,
		PrometheusMetrics: &PrometheusMetrics{Registry: registry},
	}
	
	return nil
}

// Start begins the monitoring system
func (em *EnterpriseMonitor) Start(ctx context.Context) error {
	em.mu.Lock()
	defer em.mu.Unlock()
	
	if em.running {
		return fmt.Errorf("monitor already running")
	}
	
	em.running = true
	
	// Start metrics collection
	go em.metricsCollectionLoop(ctx)
	
	// Start health checks
	go em.healthCheckLoop(ctx)
	
	// Start Prometheus metrics server
	if em.config.MetricsEnabled {
		go em.startPrometheusServer()
	}
	
	log.Printf("Enterprise monitoring started on port %d", em.config.PrometheusPort)
	return nil
}

// Stop shuts down the monitoring system
func (em *EnterpriseMonitor) Stop() error {
	em.mu.Lock()
	defer em.mu.Unlock()
	
	if !em.running {
		return nil
	}
	
	close(em.stopChan)
	em.running = false
	
	log.Printf("Enterprise monitoring stopped")
	return nil
}

// metricsCollectionLoop collects metrics at regular intervals
func (em *EnterpriseMonitor) metricsCollectionLoop(ctx context.Context) {
	ticker := time.NewTicker(em.config.MetricsInterval)
	defer ticker.Stop()
	
	for {
		select {
		case <-ctx.Done():
			return
		case <-em.stopChan:
			return
		case <-ticker.C:
			em.collectMetrics()
		}
	}
}

// collectMetrics gathers current system and application metrics
func (em *EnterpriseMonitor) collectMetrics() {
	// Collect system metrics
	em.metrics.SystemMetrics.CPUUsagePercent.Set(getCPUUsage())
	em.metrics.SystemMetrics.MemoryUsagePercent.Set(getMemoryUsage())
	em.metrics.SystemMetrics.DiskUsagePercent.Set(getDiskUsage())
	em.metrics.SystemMetrics.LoadAverage.Set(getLoadAverage())
	
	// Collect GPU metrics
	for deviceID := 0; deviceID < getGPUCount(); deviceID++ {
		deviceStr := fmt.Sprintf("%d", deviceID)
		em.metrics.GPUMetrics.GPUUtilization.WithLabelValues(deviceStr).Set(getGPUUtilization(deviceID))
		em.metrics.GPUMetrics.GPUMemoryUsage.WithLabelValues(deviceStr).Set(getGPUMemoryUsage(deviceID))
		em.metrics.GPUMetrics.GPUTemperature.WithLabelValues(deviceStr).Set(getGPUTemperature(deviceID))
		em.metrics.GPUMetrics.GPUPowerConsumption.WithLabelValues(deviceStr).Set(getGPUPowerConsumption(deviceID))
	}
}

// healthCheckLoop performs health checks at regular intervals
func (em *EnterpriseMonitor) healthCheckLoop(ctx context.Context) {
	ticker := time.NewTicker(em.config.HealthCheckInterval)
	defer ticker.Stop()
	
	for {
		select {
		case <-ctx.Done():
			return
		case <-em.stopChan:
			return
		case <-ticker.C:
			em.performHealthChecks(ctx)
		}
	}
}

// performHealthChecks executes all registered health checks
func (em *EnterpriseMonitor) performHealthChecks(ctx context.Context) {
	em.mu.RLock()
	checks := make(map[string]HealthCheck)
	for name, check := range em.healthChecks {
		checks[name] = check
	}
	em.mu.RUnlock()
	
	for name, check := range checks {
		go func(checkName string, checkFunc HealthCheck) {
			checkCtx, cancel := context.WithTimeout(ctx, em.config.HealthCheckTimeout)
			defer cancel()
			
			if err := checkFunc(checkCtx); err != nil {
				log.Printf("Health check failed: %s - %v", checkName, err)
				// Trigger alert if alerting is enabled
				if em.alerts != nil {
					em.alerts.TriggerAlert("health_check_failed", fmt.Sprintf("Health check %s failed: %v", checkName, err))
				}
			}
		}(name, check)
	}
}

// startPrometheusServer starts the Prometheus metrics HTTP server
func (em *EnterpriseMonitor) startPrometheusServer() {
	mux := http.NewServeMux()
	mux.Handle("/metrics", promhttp.HandlerFor(em.metrics.PrometheusMetrics.Registry, promhttp.HandlerOpts{}))
	
	server := &http.Server{
		Addr:              fmt.Sprintf(":%d", em.config.PrometheusPort),
		Handler:           mux,
		ReadTimeout:       10 * time.Second,
		WriteTimeout:      10 * time.Second,
		IdleTimeout:       60 * time.Second,
		ReadHeaderTimeout: 5 * time.Second,
	}
	
	log.Printf("Starting Prometheus metrics server on port %d", em.config.PrometheusPort)
	if err := server.ListenAndServe(); err != nil && err != http.ErrServerClosed {
		log.Printf("Prometheus server error: %v", err)
	}
}

// RegisterHealthCheck registers a new health check
func (em *EnterpriseMonitor) RegisterHealthCheck(name string, check HealthCheck) {
	em.mu.Lock()
	defer em.mu.Unlock()
	em.healthChecks[name] = check
}

// getCPUUsage, getMemoryUsage, getDiskUsage, getLoadAverage, getGPUCount,
// getGPUUtilization, getGPUMemoryUsage, getGPUTemperature, getGPUPowerConsumption
// are defined in platform-specific files:
//   system_metrics_linux.go   (//go:build linux)
//   system_metrics_darwin.go  (//go:build darwin)
//   system_metrics_other.go   (//go:build !linux && !darwin)