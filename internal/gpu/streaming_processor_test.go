package gpu

import (
	"os"
	"testing"
)

func TestNewStreamingProcessor(t *testing.T) {
	// Create a simple memory pool for testing
	memoryPool, err := NewMemoryPool(16.0, nil) // 16GB test pool
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	config := DefaultStreamingConfig()
	processor := NewStreamingProcessor(memoryPool, config)
	
	if processor == nil {
		t.Fatal("NewStreamingProcessor should not return nil")
	}
}

func TestStreamingProcessorGetMetrics(t *testing.T) {
	memoryPool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	processor := NewStreamingProcessor(memoryPool, nil)
	
	metrics := processor.GetMetrics()
	
	if metrics.TotalStreams < 0 {
		t.Error("TotalStreams should not be negative")
	}
	
	if metrics.ActiveStreams < 0 {
		t.Error("ActiveStreams should not be negative")
	}
	
	if metrics.CompletedStreams < 0 {
		t.Error("CompletedStreams should not be negative")
	}
	
	if metrics.FailedStreams < 0 {
		t.Error("FailedStreams should not be negative")
	}
}

func TestDefaultStreamingConfig(t *testing.T) {
	config := DefaultStreamingConfig()
	
	if config == nil {
		t.Fatal("DefaultStreamingConfig should not return nil")
	}
	
	if config.MaxConcurrentStreams <= 0 {
		t.Error("MaxConcurrentStreams should be positive")
	}
	
	if config.ChunkSizeGB <= 0 {
		t.Error("ChunkSizeGB should be positive")
	}
	
	if config.BufferSizeGB <= 0 {
		t.Error("BufferSizeGB should be positive")
	}
	
	if config.RetryAttempts < 0 {
		t.Error("RetryAttempts should not be negative")
	}
}

func TestStreamingProcessorCancelStream(t *testing.T) {
	memoryPool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	processor := NewStreamingProcessor(memoryPool, nil)
	
	// Test cancelling non-existent stream
	err = processor.CancelStream("nonexistent-stream")
	if err == nil {
		t.Error("CancelStream should error for non-existent stream")
	}
}

func TestStreamingProcessorGetStreamStatus(t *testing.T) {
	memoryPool, err := NewMemoryPool(16.0, nil)
	if err != nil {
		t.Fatalf("Failed to create memory pool: %v", err)
	}
	
	processor := NewStreamingProcessor(memoryPool, nil)
	
	// Test getting status of non-existent stream
	_, err = processor.GetStreamStatus("nonexistent-stream")
	if err == nil {
		t.Error("GetStreamStatus should error for non-existent stream")
	}
}

func TestStreamStatusString(t *testing.T) {
	testCases := []struct {
		status   StreamStatus
		expected string
	}{
		{StatusInitializing, "initializing"},
		{StatusProcessing, "processing"},
		{StatusCompleted, "completed"},
		{StatusFailed, "failed"},
		{StatusCancelled, "cancelled"},
	}
	
	for _, tc := range testCases {
		result := tc.status.String()
		if result != tc.expected {
			t.Errorf("Status %d should return '%s', got '%s'", int(tc.status), tc.expected, result)
		}
	}
}

func TestDatasetInputStructure(t *testing.T) {
	// Test that DatasetInput can be created
	input := DatasetInput{
		InputPath:  "/path/to/input.fastq",
		OutputDir:  "/path/to/output",
		SizeGB:     10.5,
		Format:     "FASTQ",
		Parameters: make(map[string]interface{}),
	}
	
	if input.SizeGB != 10.5 {
		t.Error("DatasetInput fields should be settable")
	}
	
	if input.Format != "FASTQ" {
		t.Error("DatasetInput format should be settable")
	}
}

func TestProcessingResultStructure(t *testing.T) {
	// Test that ProcessingResult can be created
	result := &ProcessingResult{
		ChunkID:    "chunk-001",
		SizeGB:     1.5,
		OutputPath: "/path/to/output.fastq",
		Checksum:   "sha256:abc123",
		Metadata:   make(map[string]interface{}),
	}
	
	if result == nil {
		t.Error("ProcessingResult should be creatable")
	}
	
	if result.SizeGB != 1.5 {
		t.Error("ProcessingResult fields should be settable")
	}
	
	if result.ChunkID != "chunk-001" {
		t.Error("ProcessingResult ChunkID should be settable")
	}
}

func TestProcessingStreamStructure(t *testing.T) {
	// Test that ProcessingStream can be created
	stream := &ProcessingStream{
		ID:              "stream-001",
		DatasetSize:     10.0,
		ChunksTotal:     10,
		ChunksProcessed: 5,
		Status:          StatusProcessing,
	}
	
	if stream == nil {
		t.Error("ProcessingStream should be creatable")
	}
	
	if stream.DatasetSize != 10.0 {
		t.Error("ProcessingStream fields should be settable")
	}
	
	if stream.Status != StatusProcessing {
		t.Error("ProcessingStream status should be settable")
	}
}

func TestStreamingMetricsStructure(t *testing.T) {
	// Test that StreamingMetrics can be created
	metrics := &StreamingMetrics{
		TotalStreams:        5,
		ActiveStreams:       2,
		CompletedStreams:    3,
		FailedStreams:       0,
		TotalDataProcessed:  25.5,
		AverageThroughput:   10.2,
		PeakThroughput:      15.8,
	}
	
	if metrics == nil {
		t.Error("StreamingMetrics should be creatable")
	}
	
	if metrics.TotalStreams != 5 {
		t.Error("StreamingMetrics fields should be settable")
	}
	
	if metrics.AverageThroughput != 10.2 {
		t.Error("StreamingMetrics throughput should be settable")
	}
}
// ---- calculateChecksum tests ----

func TestCalculateChecksum_NonExistentFile(t *testing.T) {
	// A file that does not exist must return "" (no panic, no fake hash).
	got := calculateChecksum("/nonexistent/path/that/does/not/exist.fastq")
	if got != "" {
		t.Errorf("non-existent file: expected empty string, got %q", got)
	}
}

func TestCalculateChecksum_RealFile(t *testing.T) {
	// Write a temp file, compute checksum, verify format and determinism.
	f, err := os.CreateTemp(t.TempDir(), "checksum-*.txt")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	content := []byte("ACGTACGTACGT\n")
	if _, err := f.Write(content); err != nil {
		t.Fatalf("failed to write temp file: %v", err)
	}
	f.Close()

	got1 := calculateChecksum(f.Name())
	got2 := calculateChecksum(f.Name())

	if len(got1) == 0 {
		t.Fatal("expected non-empty checksum for a real file")
	}
	if got1[:7] != "sha256:" {
		t.Errorf("checksum should start with 'sha256:', got %q", got1[:7])
	}
	if got1 != got2 {
		t.Errorf("checksum must be deterministic: first=%q second=%q", got1, got2)
	}
	// SHA-256 hex is 64 chars; "sha256:" prefix → total 71 chars.
	const wantLen = 7 + 64
	if len(got1) != wantLen {
		t.Errorf("expected checksum length %d, got %d (%q)", wantLen, len(got1), got1)
	}
}

func TestCalculateChecksum_DifferentContents(t *testing.T) {
	// Two files with different content must have different checksums.
	dir := t.TempDir()

	writeFile := func(name, content string) string {
		path := dir + "/" + name
		if err := os.WriteFile(path, []byte(content), 0600); err != nil {
			t.Fatalf("failed to write %s: %v", name, err)
		}
		return path
	}

	p1 := writeFile("a.txt", "AAAA")
	p2 := writeFile("b.txt", "TTTT")

	c1 := calculateChecksum(p1)
	c2 := calculateChecksum(p2)
	if c1 == c2 {
		t.Errorf("different files should have different checksums, both got %q", c1)
	}
}
