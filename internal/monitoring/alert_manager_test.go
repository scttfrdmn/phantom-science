package monitoring

import (
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"net/http/httptest"
	"strings"
	"sync"
	"testing"
	"time"
)

func TestNewAlertManager(t *testing.T) {
	config := &MonitoringConfig{
		AlertingEnabled: true,
		SlackWebhookURL: "https://hooks.slack.com/test",
		PagerDutyKey:    "test-key",
	}
	
	am := NewAlertManager(config)
	
	if am == nil {
		t.Fatal("NewAlertManager should not return nil")
	}
	
	if am.config != config {
		t.Error("AlertManager should use provided config")
	}
	
	if len(am.activeAlerts) != 0 {
		t.Error("AlertManager should start with no active alerts")
	}
	
	if len(am.notifiers) != 2 {
		t.Errorf("Expected 2 notifiers (Slack + PagerDuty), got %d", len(am.notifiers))
	}
}

func TestNewAlertManagerNoNotifiers(t *testing.T) {
	config := &MonitoringConfig{
		AlertingEnabled: true,
		// No webhook URLs or keys provided
	}
	
	am := NewAlertManager(config)
	
	if am == nil {
		t.Fatal("NewAlertManager should not return nil")
	}
	
	if len(am.notifiers) != 0 {
		t.Errorf("Expected 0 notifiers, got %d", len(am.notifiers))
	}
}

func TestTriggerAlert(t *testing.T) {
	config := &MonitoringConfig{AlertingEnabled: true}
	am := NewAlertManager(config)
	
	am.TriggerAlert("test-alert", "Test alert message")
	
	activeAlerts := am.GetActiveAlerts()
	if len(activeAlerts) != 1 {
		t.Errorf("Expected 1 active alert, got %d", len(activeAlerts))
	}
	
	alert := activeAlerts[0]
	if alert.Name != "test-alert" {
		t.Errorf("Expected alert name 'test-alert', got %s", alert.Name)
	}
	
	if alert.Message != "Test alert message" {
		t.Errorf("Expected alert message 'Test alert message', got %s", alert.Message)
	}
	
	if alert.Severity != SeverityWarning {
		t.Errorf("Expected default severity warning, got %s", alert.Severity)
	}
	
	if alert.Status != StatusFiring {
		t.Errorf("Expected status firing, got %s", alert.Status)
	}
	
	if alert.Count != 1 {
		t.Errorf("Expected count 1, got %d", alert.Count)
	}
}

func TestTriggerAlertWithDetails(t *testing.T) {
	config := &MonitoringConfig{AlertingEnabled: true}
	am := NewAlertManager(config)
	
	labels := map[string]string{
		"service": "test-service",
		"region":  "us-west-2",
	}
	
	annotations := map[string]string{
		"runbook": "https://example.com/runbook",
	}
	
	am.TriggerAlertWithDetails("critical-alert", "Critical issue detected", 
		SeverityCritical, labels, annotations)
	
	activeAlerts := am.GetActiveAlerts()
	if len(activeAlerts) != 1 {
		t.Errorf("Expected 1 active alert, got %d", len(activeAlerts))
	}
	
	alert := activeAlerts[0]
	if alert.Severity != SeverityCritical {
		t.Errorf("Expected severity critical, got %s", alert.Severity)
	}
	
	// Check that custom labels are preserved (note: default labels are added)
	if alert.Labels["region"] != "us-west-2" {
		t.Errorf("Expected region label 'us-west-2', got %s", alert.Labels["region"])
	}
	
	if alert.Annotations["runbook"] != "https://example.com/runbook" {
		t.Errorf("Expected runbook annotation, got %s", alert.Annotations["runbook"])
	}
}

func TestTriggerDuplicateAlert(t *testing.T) {
	config := &MonitoringConfig{AlertingEnabled: true}
	am := NewAlertManager(config)
	
	// Trigger same alert twice
	am.TriggerAlert("test-alert", "First message")
	am.TriggerAlert("test-alert", "Updated message")
	
	activeAlerts := am.GetActiveAlerts()
	if len(activeAlerts) != 1 {
		t.Errorf("Expected 1 active alert, got %d", len(activeAlerts))
	}
	
	alert := activeAlerts[0]
	if alert.Count != 2 {
		t.Errorf("Expected count 2, got %d", alert.Count)
	}
	
	if alert.Message != "Updated message" {
		t.Errorf("Expected updated message, got %s", alert.Message)
	}
}

func TestResolveAlert(t *testing.T) {
	config := &MonitoringConfig{AlertingEnabled: true}
	am := NewAlertManager(config)
	
	// Trigger alert
	am.TriggerAlert("test-alert", "Test alert message")
	
	if len(am.GetActiveAlerts()) != 1 {
		t.Error("Should have 1 active alert before resolution")
	}
	
	// Resolve alert
	am.ResolveAlert("test-alert", nil)
	
	if len(am.GetActiveAlerts()) != 0 {
		t.Error("Should have 0 active alerts after resolution")
	}
	
	// Check alert history
	history := am.GetAlertHistory(10)
	if len(history) != 1 {
		t.Errorf("Expected 1 alert in history, got %d", len(history))
	}
	
	alert := history[0]
	if alert.Status != StatusResolved {
		t.Errorf("Expected status resolved, got %s", alert.Status)
	}
	
	if alert.ResolvedAt == nil {
		t.Error("ResolvedAt should be set")
	}
}

func TestResolveNonExistentAlert(t *testing.T) {
	config := &MonitoringConfig{AlertingEnabled: true}
	am := NewAlertManager(config)
	
	// Try to resolve non-existent alert (should not panic)
	am.ResolveAlert("non-existent", nil)
	
	if len(am.GetActiveAlerts()) != 0 {
		t.Error("Should have 0 active alerts")
	}
}

func TestGetAlertHistory(t *testing.T) {
	config := &MonitoringConfig{AlertingEnabled: true}
	am := NewAlertManager(config)
	
	// Trigger multiple different alerts (same name creates one alert with count)
	for i := 0; i < 5; i++ {
		alertName := fmt.Sprintf("test-alert-%d", i)
		am.TriggerAlert(alertName, "Test message")
		time.Sleep(1 * time.Millisecond) // Ensure different timestamps
	}
	
	// Test different limits
	history := am.GetAlertHistory(3)
	if len(history) != 3 {
		t.Errorf("Expected 3 alerts in history, got %d", len(history))
	}
	
	history = am.GetAlertHistory(0)
	if len(history) != 5 {
		t.Errorf("Expected all 5 alerts in history, got %d", len(history))
	}
	
	history = am.GetAlertHistory(10)
	if len(history) != 5 {
		t.Errorf("Expected all 5 alerts in history, got %d", len(history))
	}
}

func TestAddNotifier(t *testing.T) {
	config := &MonitoringConfig{AlertingEnabled: true}
	am := NewAlertManager(config)
	
	initialCount := len(am.notifiers)
	
	// Add custom notifier
	am.AddNotifier(&TestNotifier{name: "test"})
	
	if len(am.notifiers) != initialCount+1 {
		t.Errorf("Expected %d notifiers, got %d", initialCount+1, len(am.notifiers))
	}
}

func TestNewRateLimiter(t *testing.T) {
	rl := NewRateLimiter()
	
	if rl == nil {
		t.Fatal("NewRateLimiter should not return nil")
	}
	
	if rl.limits == nil {
		t.Error("RateLimiter should initialize limits map")
	}
}

func TestRateLimiterShouldLimit(t *testing.T) {
	rl := NewRateLimiter()
	
	// Should not limit initially
	if rl.ShouldLimit("test-alert") {
		t.Error("Should not limit first alert")
	}
	
	// Trigger alerts up to limit
	for i := 1; i < 5; i++ {
		if rl.ShouldLimit("test-alert") {
			t.Errorf("Should not limit alert %d", i+1)
		}
	}
	
	// Should limit after max count
	if !rl.ShouldLimit("test-alert") {
		t.Error("Should limit after max count reached")
	}
}

func TestSlackNotifierGetName(t *testing.T) {
	sn := &SlackNotifier{}
	
	if sn.GetName() != "slack" {
		t.Errorf("Expected name 'slack', got %s", sn.GetName())
	}
}

func TestSlackNotifierGetColorForSeverity(t *testing.T) {
	sn := &SlackNotifier{}
	
	testCases := []struct {
		severity AlertSeverity
		expected string
	}{
		{SeverityInfo, "good"},
		{SeverityWarning, "warning"},
		{SeverityCritical, "danger"},
		{SeverityEmergency, "#8B0000"},
	}
	
	for _, tc := range testCases {
		result := sn.getColorForSeverity(tc.severity)
		if result != tc.expected {
			t.Errorf("For severity %s, expected color %s, got %s", tc.severity, tc.expected, result)
		}
	}
}

func TestSlackNotifierSendAlert(t *testing.T) {
	// Create test server
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.Method != "POST" {
			t.Errorf("Expected POST request, got %s", r.Method)
		}
		
		if r.Header.Get("Content-Type") != "application/json" {
			t.Errorf("Expected Content-Type application/json, got %s", r.Header.Get("Content-Type"))
		}
		
		w.WriteHeader(http.StatusOK)
	}))
	defer server.Close()
	
	sn := &SlackNotifier{
		WebhookURL: server.URL,
		Channel:    "#test",
		Username:   "test-bot",
	}
	
	alert := &Alert{
		ID:        "test-alert",
		Name:      "Test Alert",
		Severity:  SeverityWarning,
		Message:   "Test message",
		Timestamp: time.Now(),
		Status:    StatusFiring,
		Count:     1,
	}
	
	err := sn.SendAlert(alert)
	if err != nil {
		t.Errorf("SendAlert should not error: %v", err)
	}
}

func TestPagerDutyNotifierGetName(t *testing.T) {
	pd := &PagerDutyNotifier{}
	
	if pd.GetName() != "pagerduty" {
		t.Errorf("Expected name 'pagerduty', got %s", pd.GetName())
	}
}

func TestEmailNotifierGetName(t *testing.T) {
	en := &EmailNotifier{}
	
	if en.GetName() != "email" {
		t.Errorf("Expected name 'email', got %s", en.GetName())
	}
}

func TestEmailNotifierSendAlert(t *testing.T) {
	en := &EmailNotifier{}
	
	alert := &Alert{
		Name:    "Test Alert",
		Message: "Test message",
	}
	
	// Should not error (placeholder implementation)
	err := en.SendAlert(alert)
	if err != nil {
		t.Errorf("SendAlert should not error: %v", err)
	}
}

func TestGenerateAlertID(t *testing.T) {
	// Test with no labels
	id1 := generateAlertID("test-alert", nil)
	if id1 != "test-alert" {
		t.Errorf("Expected 'test-alert', got %s", id1)
	}

	// Test with a single label (deterministic - no map ordering issue)
	singleLabel := map[string]string{"service": "test"}
	id2 := generateAlertID("test-alert", singleLabel)
	if id2 == "test-alert" {
		t.Error("Alert ID should include labels")
	}
	if !strings.Contains(id2, "service") || !strings.Contains(id2, "test") {
		t.Errorf("Alert ID should contain label key and value, got: %s", id2)
	}

	// Same single-label inputs should produce same ID (deterministic)
	singleLabel2 := map[string]string{"service": "test"}
	id3 := generateAlertID("test-alert", singleLabel2)
	if id2 != id3 {
		t.Errorf("Same inputs should produce same alert ID: %s vs %s", id2, id3)
	}
}

func TestGetSeverityPriority(t *testing.T) {
	testCases := []struct {
		severity AlertSeverity
		expected int
	}{
		{SeverityInfo, 1},
		{SeverityWarning, 2},
		{SeverityCritical, 3},
		{SeverityEmergency, 4},
		{"unknown", 1},
	}
	
	for _, tc := range testCases {
		result := GetSeverityPriority(tc.severity)
		if result != tc.expected {
			t.Errorf("For severity %s, expected priority %d, got %d", tc.severity, tc.expected, result)
		}
	}
}

func TestAlertSeverityValues(t *testing.T) {
	if SeverityInfo != "info" {
		t.Errorf("Expected SeverityInfo to be 'info', got %s", SeverityInfo)
	}
	
	if SeverityWarning != "warning" {
		t.Errorf("Expected SeverityWarning to be 'warning', got %s", SeverityWarning)
	}
	
	if SeverityCritical != "critical" {
		t.Errorf("Expected SeverityCritical to be 'critical', got %s", SeverityCritical)
	}
	
	if SeverityEmergency != "emergency" {
		t.Errorf("Expected SeverityEmergency to be 'emergency', got %s", SeverityEmergency)
	}
}

func TestAlertStatusValues(t *testing.T) {
	if StatusFiring != "firing" {
		t.Errorf("Expected StatusFiring to be 'firing', got %s", StatusFiring)
	}
	
	if StatusResolved != "resolved" {
		t.Errorf("Expected StatusResolved to be 'resolved', got %s", StatusResolved)
	}
	
	if StatusSilenced != "silenced" {
		t.Errorf("Expected StatusSilenced to be 'silenced', got %s", StatusSilenced)
	}
}

// TestNotifier is a test implementation of the Notifier interface
type TestNotifier struct {
	name string
}

func (tn *TestNotifier) SendAlert(alert *Alert) error {
	return nil
}

func (tn *TestNotifier) GetName() string {
	return tn.name
}
// roundTripFunc is a test helper that implements http.RoundTripper.
type roundTripFunc func(req *http.Request) (*http.Response, error)

func (f roundTripFunc) RoundTrip(req *http.Request) (*http.Response, error) {
	return f(req)
}

// ---------- PagerDutyNotifier.SendAlert ----------

func TestPagerDutyNotifier_SendAlert_Trigger(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.WriteHeader(http.StatusAccepted)
		fmt.Fprint(w, `{"status":"success","dedup_key":"test"}`)
	}))
	defer srv.Close()

	pd := &PagerDutyNotifier{
		IntegrationKey: "test-key",
		Client: &http.Client{
			Transport: roundTripFunc(func(req *http.Request) (*http.Response, error) {
				// Redirect the hardcoded PagerDuty URL to our test server
				req.URL.Host = srv.Listener.Addr().String()
				req.URL.Scheme = "http"
				return http.DefaultTransport.RoundTrip(req)
			}),
		},
	}

	alert := &Alert{
		ID:        "test-alert-1",
		Name:      "TestAlert",
		Message:   "Test alert message",
		Severity:  SeverityCritical,
		Status:    StatusFiring,
		Timestamp: time.Now(),
	}

	err := pd.SendAlert(alert)
	if err != nil {
		t.Errorf("PagerDutyNotifier.SendAlert expected no error, got: %v", err)
	}
}

func TestPagerDutyNotifier_SendAlert_Resolve(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.WriteHeader(http.StatusAccepted)
		fmt.Fprint(w, `{"status":"success","dedup_key":"test"}`)
	}))
	defer srv.Close()

	pd := &PagerDutyNotifier{
		IntegrationKey: "test-key",
		Client: &http.Client{
			Transport: roundTripFunc(func(req *http.Request) (*http.Response, error) {
				req.URL.Host = srv.Listener.Addr().String()
				req.URL.Scheme = "http"
				return http.DefaultTransport.RoundTrip(req)
			}),
		},
	}

	alert := &Alert{
		ID:        "test-alert-2",
		Name:      "TestAlert",
		Message:   "Resolved",
		Severity:  SeverityInfo,
		Status:    StatusResolved,
		Timestamp: time.Now(),
	}

	err := pd.SendAlert(alert)
	if err != nil {
		t.Errorf("PagerDutyNotifier.SendAlert(resolve) expected no error, got: %v", err)
	}
}

func TestPagerDutyNotifier_SendAlert_ServerError(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.WriteHeader(http.StatusInternalServerError)
		fmt.Fprint(w, "internal server error")
	}))
	defer srv.Close()

	pd := &PagerDutyNotifier{
		IntegrationKey: "test-key",
		Client: &http.Client{
			Transport: roundTripFunc(func(req *http.Request) (*http.Response, error) {
				req.URL.Host = srv.Listener.Addr().String()
				req.URL.Scheme = "http"
				return http.DefaultTransport.RoundTrip(req)
			}),
		},
	}

	alert := &Alert{
		ID:        "test-alert-3",
		Name:      "TestAlert",
		Message:   "Error alert",
		Severity:  SeverityCritical,
		Status:    StatusFiring,
		Timestamp: time.Now(),
	}

	err := pd.SendAlert(alert)
	if err == nil {
		t.Error("PagerDutyNotifier.SendAlert with server error expected non-nil error")
	}
}

// ---------- PagerDutyNotifier.GetName ----------

func TestPagerDutyNotifier_GetName(t *testing.T) {
	pd := &PagerDutyNotifier{}
	if got := pd.GetName(); got != "pagerduty" {
		t.Errorf("PagerDutyNotifier.GetName() = %q, want pagerduty", got)
	}
}

// ---------- io.ReadAll error-path tests (#43) ----------

// errReader is an io.Reader that always returns an error.
type errReader struct{ msg string }

func (e errReader) Read(_ []byte) (int, error) {
	return 0, fmt.Errorf("%s", e.msg)
}

func TestSlackNotifier_SendAlert_BodyReadError(t *testing.T) {
	// Return a non-200 status with a body that fails to read.
	pd := &SlackNotifier{
		WebhookURL: "http://unused",
		Channel:    "#test",
	}

	// Swap out http.DefaultTransport so we control the response.
	orig := http.DefaultTransport
	t.Cleanup(func() { http.DefaultTransport = orig })
	http.DefaultTransport = roundTripFunc(func(req *http.Request) (*http.Response, error) {
		return &http.Response{
			StatusCode: http.StatusInternalServerError,
			Body:       io.NopCloser(errReader{"body read failed"}),
		}, nil
	})

	alert := &Alert{
		ID: "s1", Name: "S", Message: "m", Severity: SeverityCritical,
		Status: StatusFiring, Timestamp: time.Now(),
	}
	err := pd.SendAlert(alert)
	if err == nil {
		t.Fatal("expected an error, got nil")
	}
	if !strings.Contains(err.Error(), "failed to read body") {
		t.Errorf("error should mention 'failed to read body', got: %v", err)
	}
	if !strings.Contains(err.Error(), "500") {
		t.Errorf("error should contain status code 500, got: %v", err)
	}
}

func TestPagerDutyNotifier_SendAlert_BodyReadError(t *testing.T) {
	pd := &PagerDutyNotifier{
		IntegrationKey: "key",
		Client: &http.Client{
			Transport: roundTripFunc(func(req *http.Request) (*http.Response, error) {
				return &http.Response{
					StatusCode: http.StatusInternalServerError,
					Body:       io.NopCloser(errReader{"body read failed"}),
				}, nil
			}),
		},
	}

	alert := &Alert{
		ID: "p1", Name: "P", Message: "m", Severity: SeverityCritical,
		Status: StatusFiring, Timestamp: time.Now(),
	}
	err := pd.SendAlert(alert)
	if err == nil {
		t.Fatal("expected an error, got nil")
	}
	if !strings.Contains(err.Error(), "failed to read body") {
		t.Errorf("error should mention 'failed to read body', got: %v", err)
	}
	if !strings.Contains(err.Error(), "500") {
		t.Errorf("error should contain status code 500, got: %v", err)
	}
}

// ---------- Credential serialization tests (#42) ----------

func TestPagerDutyNotifier_IntegrationKey_NotInJSON(t *testing.T) {
	pd := &PagerDutyNotifier{IntegrationKey: "secret-pd-key"}
	data, err := json.Marshal(pd)
	if err != nil {
		t.Fatalf("json.Marshal failed: %v", err)
	}
	if strings.Contains(string(data), "secret-pd-key") {
		t.Errorf("PagerDutyNotifier.IntegrationKey must not appear in JSON output; got: %s", data)
	}
	if strings.Contains(string(data), "IntegrationKey") {
		t.Errorf("PagerDutyNotifier.IntegrationKey field name must not appear in JSON output; got: %s", data)
	}
}

func TestEmailNotifier_Password_NotInJSON(t *testing.T) {
	en := &EmailNotifier{
		SMTPHost:    "smtp.example.com",
		SMTPPort:    587,
		Username:    "user",
		Password:    "super-secret-password",
		FromAddress: "from@example.com",
		ToAddresses: []string{"to@example.com"},
	}
	data, err := json.Marshal(en)
	if err != nil {
		t.Fatalf("json.Marshal failed: %v", err)
	}
	if strings.Contains(string(data), "super-secret-password") {
		t.Errorf("EmailNotifier.Password must not appear in JSON output; got: %s", data)
	}
	if strings.Contains(string(data), `"Password"`) {
		t.Errorf("EmailNotifier.Password field name must not appear in JSON output; got: %s", data)
	}
}

// ---------- Race-condition regression test (#40) ----------

// countingNotifier counts how many alerts it receives, thread-safely.
type countingNotifier struct {
	mu    sync.Mutex
	count int
}

func (n *countingNotifier) SendAlert(_ *Alert) error {
	n.mu.Lock()
	n.count++
	n.mu.Unlock()
	return nil
}

func (n *countingNotifier) GetName() string { return "counting" }

func (n *countingNotifier) Count() int {
	n.mu.Lock()
	defer n.mu.Unlock()
	return n.count
}

// TestSendNotifications_NoRaceWithAddNotifier verifies that concurrent calls to
// TriggerAlert and AddNotifier do not race on the notifiers slice.
// Run with: go test -race ./internal/monitoring/...
func TestSendNotifications_NoRaceWithAddNotifier(t *testing.T) {
	am := NewAlertManager(&MonitoringConfig{AlertingEnabled: true})

	// Seed one notifier so sendNotifications has work to do from the start.
	seed := &countingNotifier{}
	am.AddNotifier(seed)

	const goroutines = 20
	var wg sync.WaitGroup

	// Concurrently trigger alerts.
	for i := 0; i < goroutines; i++ {
		wg.Add(1)
		go func(idx int) {
			defer wg.Done()
			am.TriggerAlert(fmt.Sprintf("race-alert-%d", idx), "concurrent trigger")
		}(i)
	}

	// Concurrently add notifiers while alerts are being triggered.
	for i := 0; i < goroutines; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			am.AddNotifier(&countingNotifier{})
		}()
	}

	wg.Wait()

	// Allow in-flight notification goroutines to finish.
	time.Sleep(50 * time.Millisecond)

	// The seed notifier must have received at least one alert (sanity check).
	if seed.Count() == 0 {
		t.Error("seed notifier should have received at least one alert")
	}
}
