// Package monitoring provides enterprise alerting capabilities
package monitoring

import (
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"net/http"
	"sync"
	"time"
)

// AlertManager handles enterprise alerting and notifications
type AlertManager struct {
	config         *MonitoringConfig
	activeAlerts   map[string]*Alert
	alertHistory   []*Alert
	notifiers      []Notifier
	mu             sync.RWMutex
	rateLimiter    *RateLimiter
}

// Alert represents an active alert
type Alert struct {
	ID          string            `json:"id"`
	Name        string            `json:"name"`
	Severity    AlertSeverity     `json:"severity"`
	Message     string            `json:"message"`
	Timestamp   time.Time         `json:"timestamp"`
	Status      AlertStatus       `json:"status"`
	Labels      map[string]string `json:"labels"`
	Annotations map[string]string `json:"annotations"`
	ResolvedAt  *time.Time        `json:"resolved_at,omitempty"`
	Count       int               `json:"count"`
	LastSeen    time.Time         `json:"last_seen"`
}

// AlertSeverity defines alert severity levels
type AlertSeverity string

const (
	SeverityInfo     AlertSeverity = "info"
	SeverityWarning  AlertSeverity = "warning"
	SeverityCritical AlertSeverity = "critical"
	SeverityEmergency AlertSeverity = "emergency"
)

// AlertStatus defines alert status
type AlertStatus string

const (
	StatusFiring   AlertStatus = "firing"
	StatusResolved AlertStatus = "resolved"
	StatusSilenced AlertStatus = "silenced"
)

// Notifier interface for different notification channels
type Notifier interface {
	SendAlert(alert *Alert) error
	GetName() string
}

// RateLimiter prevents alert spam
type RateLimiter struct {
	limits map[string]*RateLimit
	mu     sync.RWMutex
}

// RateLimit tracks rate limiting for specific alert types
type RateLimit struct {
	Count     int
	Window    time.Duration
	LastReset time.Time
	MaxCount  int
}

// SlackNotifier sends alerts to Slack
type SlackNotifier struct {
	WebhookURL string
	Channel    string
	Username   string
}

// PagerDutyNotifier sends alerts to PagerDuty
type PagerDutyNotifier struct {
	IntegrationKey string `json:"-" yaml:"-"`
	Client         *http.Client
}

// EmailNotifier sends alerts via email
type EmailNotifier struct {
	SMTPHost     string
	SMTPPort     int
	Username     string
	Password     string `json:"-" yaml:"-"` // #nosec G117
	FromAddress  string
	ToAddresses  []string
}

// NewAlertManager creates a new alert manager
func NewAlertManager(config *MonitoringConfig) *AlertManager {
	am := &AlertManager{
		config:       config,
		activeAlerts: make(map[string]*Alert),
		alertHistory: make([]*Alert, 0),
		notifiers:    make([]Notifier, 0),
		rateLimiter:  NewRateLimiter(),
	}
	
	// Initialize notifiers based on configuration
	if config.SlackWebhookURL != "" {
		am.AddNotifier(&SlackNotifier{
			WebhookURL: config.SlackWebhookURL,
			Channel:    "#phantom-alerts",
			Username:   "phantom-monitor",
		})
	}
	
	if config.PagerDutyKey != "" {
		am.AddNotifier(&PagerDutyNotifier{
			IntegrationKey: config.PagerDutyKey,
			Client:         &http.Client{Timeout: 30 * time.Second},
		})
	}
	
	return am
}

// NewRateLimiter creates a new rate limiter
func NewRateLimiter() *RateLimiter {
	return &RateLimiter{
		limits: make(map[string]*RateLimit),
	}
}

// TriggerAlert triggers a new alert or updates an existing one
func (am *AlertManager) TriggerAlert(name, message string) {
	am.TriggerAlertWithDetails(name, message, SeverityWarning, nil, nil)
}

// TriggerAlertWithDetails triggers an alert with full details
func (am *AlertManager) TriggerAlertWithDetails(name, message string, severity AlertSeverity, 
	labels, annotations map[string]string) {
	
	am.mu.Lock()
	defer am.mu.Unlock()
	
	alertID := generateAlertID(name, labels)
	
	// Check if alert already exists
	if existingAlert, exists := am.activeAlerts[alertID]; exists {
		// Update existing alert
		existingAlert.Count++
		existingAlert.LastSeen = time.Now()
		existingAlert.Message = message // Update message in case it changed
		
		// Check rate limiting
		if am.rateLimiter.ShouldLimit(name) {
			log.Printf("Alert %s rate limited", name)
			return
		}
	} else {
		// Create new alert
		alert := &Alert{
			ID:          alertID,
			Name:        name,
			Severity:    severity,
			Message:     message,
			Timestamp:   time.Now(),
			Status:      StatusFiring,
			Labels:      labels,
			Annotations: annotations,
			Count:       1,
			LastSeen:    time.Now(),
		}
		
		if labels == nil {
			alert.Labels = make(map[string]string)
		}
		if annotations == nil {
			alert.Annotations = make(map[string]string)
		}
		
		// Add default labels
		alert.Labels["service"] = "phantom-science"
		alert.Labels["environment"] = "production"
		
		am.activeAlerts[alertID] = alert
		am.alertHistory = append(am.alertHistory, alert)
		
		// Send notifications
		go am.sendNotifications(alert)
		
		log.Printf("Alert triggered: %s - %s", name, message)
	}
}

// ResolveAlert resolves an active alert
func (am *AlertManager) ResolveAlert(name string, labels map[string]string) {
	am.mu.Lock()
	defer am.mu.Unlock()
	
	alertID := generateAlertID(name, labels)
	
	if alert, exists := am.activeAlerts[alertID]; exists {
		now := time.Now()
		alert.Status = StatusResolved
		alert.ResolvedAt = &now
		
		delete(am.activeAlerts, alertID)
		
		// Send resolution notification
		go am.sendNotifications(alert)
		
		log.Printf("Alert resolved: %s", name)
	}
}

// sendNotifications sends alerts to all configured notifiers.
// It snapshots the notifiers slice under RLock so that concurrent calls to
// AddNotifier cannot race with the iteration.
func (am *AlertManager) sendNotifications(alert *Alert) {
	am.mu.RLock()
	notifiers := make([]Notifier, len(am.notifiers))
	copy(notifiers, am.notifiers)
	am.mu.RUnlock()

	for _, notifier := range notifiers {
		if err := notifier.SendAlert(alert); err != nil {
			log.Printf("Failed to send alert via %s: %v", notifier.GetName(), err)
		}
	}
}

// AddNotifier adds a new notification channel
func (am *AlertManager) AddNotifier(notifier Notifier) {
	am.mu.Lock()
	defer am.mu.Unlock()
	am.notifiers = append(am.notifiers, notifier)
}

// GetActiveAlerts returns all active alerts
func (am *AlertManager) GetActiveAlerts() []*Alert {
	am.mu.RLock()
	defer am.mu.RUnlock()
	
	alerts := make([]*Alert, 0, len(am.activeAlerts))
	for _, alert := range am.activeAlerts {
		alerts = append(alerts, alert)
	}
	return alerts
}

// GetAlertHistory returns alert history
func (am *AlertManager) GetAlertHistory(limit int) []*Alert {
	am.mu.RLock()
	defer am.mu.RUnlock()
	
	if limit <= 0 || limit > len(am.alertHistory) {
		limit = len(am.alertHistory)
	}
	
	// Return most recent alerts
	start := len(am.alertHistory) - limit
	return am.alertHistory[start:]
}

// ShouldLimit checks if an alert should be rate limited
func (rl *RateLimiter) ShouldLimit(alertName string) bool {
	rl.mu.Lock()
	defer rl.mu.Unlock()
	
	now := time.Now()
	
	limit, exists := rl.limits[alertName]
	if !exists {
		// Create new rate limit
		limit = &RateLimit{
			Count:     0,
			Window:    5 * time.Minute, // 5 minute window
			LastReset: now,
			MaxCount:  5, // Max 5 alerts per window
		}
		rl.limits[alertName] = limit
	}
	
	// Reset if window has expired
	if now.Sub(limit.LastReset) > limit.Window {
		limit.Count = 0
		limit.LastReset = now
	}
	
	limit.Count++
	return limit.Count > limit.MaxCount
}

// SlackNotifier implementation

func (sn *SlackNotifier) SendAlert(alert *Alert) error {
	payload := map[string]interface{}{
		"channel":   sn.Channel,
		"username":  sn.Username,
		"icon_emoji": ":warning:",
		"attachments": []map[string]interface{}{
			{
				"color":     sn.getColorForSeverity(alert.Severity),
				"title":     fmt.Sprintf("Alert: %s", alert.Name),
				"text":      alert.Message,
				"timestamp": alert.Timestamp.Unix(),
				"fields": []map[string]interface{}{
					{
						"title": "Severity",
						"value": string(alert.Severity),
						"short": true,
					},
					{
						"title": "Status", 
						"value": string(alert.Status),
						"short": true,
					},
					{
						"title": "Count",
						"value": alert.Count,
						"short": true,
					},
				},
			},
		},
	}
	
	jsonPayload, err := json.Marshal(payload)
	if err != nil {
		return fmt.Errorf("failed to marshal Slack payload: %w", err)
	}
	
	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()
	req, err := http.NewRequestWithContext(ctx, http.MethodPost, sn.WebhookURL, bytes.NewBuffer(jsonPayload))
	if err != nil {
		return fmt.Errorf("failed to create Slack request: %w", err)
	}
	req.Header.Set("Content-Type", "application/json")
	resp, err := http.DefaultClient.Do(req)
	if err != nil {
		return fmt.Errorf("failed to send Slack notification: %w", err)
	}
	defer resp.Body.Close()
	
	if resp.StatusCode != http.StatusOK {
		body, readErr := io.ReadAll(resp.Body)
		if readErr != nil {
			return fmt.Errorf("Slack API returned status %d (failed to read body: %w)", resp.StatusCode, readErr)
		}
		return fmt.Errorf("Slack API returned status %d: %s", resp.StatusCode, string(body))
	}
	
	return nil
}

func (sn *SlackNotifier) GetName() string {
	return "slack"
}

func (sn *SlackNotifier) getColorForSeverity(severity AlertSeverity) string {
	switch severity {
	case SeverityInfo:
		return "good"
	case SeverityWarning:
		return "warning" 
	case SeverityCritical:
		return "danger"
	case SeverityEmergency:
		return "#8B0000" // Dark red
	default:
		return "warning"
	}
}

// PagerDutyNotifier implementation

func (pd *PagerDutyNotifier) SendAlert(alert *Alert) error {
	eventAction := "trigger"
	if alert.Status == StatusResolved {
		eventAction = "resolve"
	}
	
	payload := map[string]interface{}{
		"routing_key":  pd.IntegrationKey,
		"event_action": eventAction,
		"dedup_key":    alert.ID,
		"payload": map[string]interface{}{
			"summary":        fmt.Sprintf("%s: %s", alert.Name, alert.Message),
			"source":         "phantom-science",
			"severity":       string(alert.Severity),
			"timestamp":      alert.Timestamp.Format(time.RFC3339),
			"custom_details": alert.Labels,
		},
	}
	
	jsonPayload, err := json.Marshal(payload)
	if err != nil {
		return fmt.Errorf("failed to marshal PagerDuty payload: %w", err)
	}
	
	resp, err := pd.Client.Post("https://events.pagerduty.com/v2/enqueue", 
		"application/json", bytes.NewBuffer(jsonPayload))
	if err != nil {
		return fmt.Errorf("failed to send PagerDuty notification: %w", err)
	}
	defer resp.Body.Close()
	
	if resp.StatusCode != http.StatusAccepted {
		body, readErr := io.ReadAll(resp.Body)
		if readErr != nil {
			return fmt.Errorf("PagerDuty API returned status %d (failed to read body: %w)", resp.StatusCode, readErr)
		}
		return fmt.Errorf("PagerDuty API returned status %d: %s", resp.StatusCode, string(body))
	}
	
	return nil
}

func (pd *PagerDutyNotifier) GetName() string {
	return "pagerduty"
}

// EmailNotifier implementation

func (en *EmailNotifier) SendAlert(alert *Alert) error {
	// Placeholder for email implementation
	// In a real implementation, this would use SMTP to send emails
	log.Printf("Would send email alert: %s - %s", alert.Name, alert.Message)
	return nil
}

func (en *EmailNotifier) GetName() string {
	return "email"
}

// Helper functions

func generateAlertID(name string, labels map[string]string) string {
	// Simple alert ID generation based on name and labels
	// In practice, this would be more sophisticated
	id := name
	if labels != nil {
		for k, v := range labels {
			id += fmt.Sprintf("-%s-%s", k, v)
		}
	}
	return id
}

// GetSeverityPriority returns numeric priority for severity comparison
func GetSeverityPriority(severity AlertSeverity) int {
	switch severity {
	case SeverityInfo:
		return 1
	case SeverityWarning:
		return 2
	case SeverityCritical:
		return 3
	case SeverityEmergency:
		return 4
	default:
		return 1
	}
}