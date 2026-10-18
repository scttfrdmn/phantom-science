package gpu

import (
	"context"
	"crypto/sha256"
	"fmt"
	"io"
	"log"
	"os"
	"sync"
	"time"
)

// StreamingProcessor simulates chunked processing of a dataset against a
// MemoryPool; see the package comment.
type StreamingProcessor struct {
	memoryPool    *MemoryPool
	config        *StreamingConfig
	activeStreams map[string]*ProcessingStream
	mu            sync.RWMutex
	metricsMu     sync.RWMutex // guards metrics; kept off StreamingMetrics so GetMetrics can return a plain copy
	metrics       *StreamingMetrics
}

// StreamingConfig defines configuration for streaming processing
type StreamingConfig struct {
	MaxConcurrentStreams int           // Maximum concurrent processing streams
	ChunkSizeGB          float64       // Size of each processing chunk
	BufferSizeGB         float64       // Size of input/output buffers
	TimeoutDuration      time.Duration // Timeout for processing operations
	RetryAttempts        int           // Number of retry attempts on failure
	CompressionEnabled   bool          // Enable data compression
	ChecksumValidation   bool          // Enable data integrity validation
}

// ProcessingStream represents an active streaming processing operation
type ProcessingStream struct {
	ID              string
	DatasetSize     float64
	ProcessedSize   float64
	ChunksTotal     int
	ChunksProcessed int
	StartTime       time.Time
	Status          StreamStatus
	Error           error
	ResultChan      chan *ProcessingResult
	CancelFunc      context.CancelFunc
	mu              sync.RWMutex
}

// StreamStatus represents the status of a processing stream
type StreamStatus int

const (
	StatusInitializing StreamStatus = iota
	StatusProcessing
	StatusCompleted
	StatusFailed
	StatusCancelled
)

// ProcessingResult contains the result of a chunk processing operation
type ProcessingResult struct {
	ChunkID    string
	SizeGB     float64
	Duration   time.Duration
	OutputPath string
	Checksum   string
	Error      error
	Metadata   map[string]interface{}
}

// StreamingMetrics tracks streaming processing performance
type StreamingMetrics struct {
	TotalStreams       int
	ActiveStreams      int
	CompletedStreams   int
	FailedStreams      int
	TotalDataProcessed float64
	AverageThroughput  float64
	PeakThroughput     float64
}

// DefaultStreamingConfig returns the default streaming configuration.
func DefaultStreamingConfig() *StreamingConfig {
	return &StreamingConfig{
		MaxConcurrentStreams: 4,
		ChunkSizeGB:          1.0,
		BufferSizeGB:         2.0,
		TimeoutDuration:      30 * time.Minute,
		RetryAttempts:        3,
		CompressionEnabled:   true,
		ChecksumValidation:   true,
	}
}

// NewStreamingProcessor creates a new streaming processor
func NewStreamingProcessor(memoryPool *MemoryPool, config *StreamingConfig) *StreamingProcessor {
	if config == nil {
		config = DefaultStreamingConfig()
	}

	return &StreamingProcessor{
		memoryPool:    memoryPool,
		config:        config,
		activeStreams: make(map[string]*ProcessingStream),
		metrics:       &StreamingMetrics{},
	}
}

// ProcessDataset starts a simulated processing stream for input and returns
// it; chunks are "processed" asynchronously by waiting, not by computation.
func (sp *StreamingProcessor) ProcessDataset(ctx context.Context, input DatasetInput) (*ProcessingStream, error) {
	sp.mu.Lock()
	defer sp.mu.Unlock()

	// Check concurrent stream limit
	if len(sp.activeStreams) >= sp.config.MaxConcurrentStreams {
		return nil, fmt.Errorf("maximum concurrent streams (%d) reached", sp.config.MaxConcurrentStreams)
	}

	// Create processing stream
	streamID := generateStreamID()
	chunksTotal := int((input.SizeGB / sp.config.ChunkSizeGB) + 1)

	ctx, cancel := context.WithTimeout(ctx, sp.config.TimeoutDuration)

	stream := &ProcessingStream{
		ID:          streamID,
		DatasetSize: input.SizeGB,
		ChunksTotal: chunksTotal,
		StartTime:   time.Now(),
		Status:      StatusInitializing,
		ResultChan:  make(chan *ProcessingResult, chunksTotal),
		CancelFunc:  cancel,
	}

	sp.activeStreams[streamID] = stream

	// Update metrics
	sp.metricsMu.Lock()
	sp.metrics.TotalStreams++
	sp.metrics.ActiveStreams++
	sp.metricsMu.Unlock()

	// Start processing in background
	go sp.processStreamAsync(ctx, stream, input)

	log.Printf("Started streaming processing for dataset %.2fGB (Stream ID: %s, %d chunks)",
		input.SizeGB, streamID, chunksTotal)

	return stream, nil
}

// processStreamAsync handles asynchronous stream processing
func (sp *StreamingProcessor) processStreamAsync(ctx context.Context, stream *ProcessingStream, input DatasetInput) {
	defer func() {
		close(stream.ResultChan)
		sp.cleanupStream(stream.ID)
	}()

	stream.mu.Lock()
	stream.Status = StatusProcessing
	stream.mu.Unlock()

	// Process dataset in chunks
	for chunkIndex := 0; chunkIndex < stream.ChunksTotal; chunkIndex++ {
		select {
		case <-ctx.Done():
			stream.mu.Lock()
			stream.Status = StatusCancelled
			stream.Error = ctx.Err()
			stream.mu.Unlock()
			return
		default:
		}

		// Process individual chunk
		result := sp.processChunk(ctx, stream, chunkIndex, input)

		// Update stream progress
		stream.mu.Lock()
		stream.ChunksProcessed++
		stream.ProcessedSize += result.SizeGB

		if result.Error != nil {
			stream.Status = StatusFailed
			stream.Error = result.Error
			stream.mu.Unlock()
			stream.ResultChan <- result
			return
		}
		stream.mu.Unlock()

		// Send result
		stream.ResultChan <- result

		// Update metrics
		sp.updateThroughputMetrics(result.SizeGB, result.Duration)
	}

	// Mark stream as completed
	stream.mu.Lock()
	stream.Status = StatusCompleted
	stream.mu.Unlock()

	log.Printf("Completed streaming processing for Stream ID: %s", stream.ID)
}

// processChunk processes an individual data chunk
func (sp *StreamingProcessor) processChunk(ctx context.Context, stream *ProcessingStream, chunkIndex int, input DatasetInput) *ProcessingResult {
	chunkID := fmt.Sprintf("%s-chunk-%d", stream.ID, chunkIndex)
	startTime := time.Now()

	// Calculate chunk parameters
	chunkOffset := float64(chunkIndex) * sp.config.ChunkSizeGB
	chunkSize := sp.config.ChunkSizeGB
	if chunkOffset+chunkSize > input.SizeGB {
		chunkSize = input.SizeGB - chunkOffset
	}

	// Allocate GPU memory for chunk processing
	allocation, err := sp.memoryPool.AllocateMemory(chunkSize*2, // Input + output buffers
		fmt.Sprintf("chunk-processing-%s", chunkID), PriorityNormal)
	if err != nil {
		return &ProcessingResult{
			ChunkID:  chunkID,
			SizeGB:   chunkSize,
			Duration: time.Since(startTime),
			Error:    fmt.Errorf("GPU memory allocation failed: %w", err),
		}
	}

	defer func() {
		if err := sp.memoryPool.FreeMemory(allocation.ID); err != nil {
			log.Printf("Warning: Failed to free GPU memory for chunk %s: %v", chunkID, err)
		}
	}()

	// Process chunk with retry logic
	var result *ProcessingResult
	for attempt := 0; attempt <= sp.config.RetryAttempts; attempt++ {
		result = sp.executeChunkProcessing(ctx, chunkID, chunkSize, chunkOffset, input, allocation)

		if result.Error == nil {
			break
		}

		if attempt < sp.config.RetryAttempts {
			log.Printf("Chunk processing failed (attempt %d/%d): %v",
				attempt+1, sp.config.RetryAttempts+1, result.Error)
			time.Sleep(time.Duration(attempt+1) * time.Second) // Exponential backoff
		}
	}

	result.Duration = time.Since(startTime)
	return result
}

// executeChunkProcessing simulates one chunk: it waits 100 ms per GB and
// returns a result naming an output file that is never written.
func (sp *StreamingProcessor) executeChunkProcessing(ctx context.Context, chunkID string, chunkSize, chunkOffset float64,
	input DatasetInput, allocation *Allocation) *ProcessingResult {

	// Simulate chunk processing (replace with actual GPU processing)
	result := &ProcessingResult{
		ChunkID:    chunkID,
		SizeGB:     chunkSize,
		OutputPath: fmt.Sprintf("%s/%s.fastq", input.OutputDir, chunkID),
		Metadata:   make(map[string]interface{}),
	}

	// Simulate processing time based on chunk size
	processingTime := time.Duration(chunkSize*100) * time.Millisecond
	select {
	case <-ctx.Done():
		result.Error = ctx.Err()
		return result
	case <-time.After(processingTime):
		// Processing completed
	}

	// Calculate checksum if enabled
	if sp.config.ChecksumValidation {
		result.Checksum = calculateChecksum(result.OutputPath)
	}

	// Add processing metadata
	result.Metadata["device_id"] = allocation.DeviceID
	result.Metadata["allocation_id"] = allocation.ID
	result.Metadata["compression_enabled"] = sp.config.CompressionEnabled
	result.Metadata["chunk_offset"] = chunkOffset

	return result
}

// GetStreamStatus returns the current status of a processing stream
func (sp *StreamingProcessor) GetStreamStatus(streamID string) (*ProcessingStream, error) {
	sp.mu.RLock()
	defer sp.mu.RUnlock()

	stream, exists := sp.activeStreams[streamID]
	if !exists {
		return nil, fmt.Errorf("stream %s not found", streamID)
	}

	// Return a copy to avoid race conditions
	stream.mu.RLock()
	defer stream.mu.RUnlock()

	return &ProcessingStream{
		ID:              stream.ID,
		DatasetSize:     stream.DatasetSize,
		ProcessedSize:   stream.ProcessedSize,
		ChunksTotal:     stream.ChunksTotal,
		ChunksProcessed: stream.ChunksProcessed,
		StartTime:       stream.StartTime,
		Status:          stream.Status,
		Error:           stream.Error,
	}, nil
}

// CancelStream cancels an active processing stream
func (sp *StreamingProcessor) CancelStream(streamID string) error {
	sp.mu.RLock()
	stream, exists := sp.activeStreams[streamID]
	sp.mu.RUnlock()

	if !exists {
		return fmt.Errorf("stream %s not found", streamID)
	}

	stream.CancelFunc()
	log.Printf("Cancelled processing stream: %s", streamID)
	return nil
}

// cleanupStream removes completed or failed streams from active list
func (sp *StreamingProcessor) cleanupStream(streamID string) {
	sp.mu.Lock()
	defer sp.mu.Unlock()

	if stream, exists := sp.activeStreams[streamID]; exists {
		delete(sp.activeStreams, streamID)

		// Update metrics
		sp.metricsMu.Lock()
		sp.metrics.ActiveStreams--
		sp.metrics.TotalDataProcessed += stream.ProcessedSize

		switch stream.Status {
		case StatusCompleted:
			sp.metrics.CompletedStreams++
		case StatusFailed, StatusCancelled:
			sp.metrics.FailedStreams++
		}
		sp.metricsMu.Unlock()
	}
}

// updateThroughputMetrics updates throughput statistics
func (sp *StreamingProcessor) updateThroughputMetrics(sizeGB float64, duration time.Duration) {
	throughput := sizeGB / duration.Seconds() * 1024 // MB/s

	sp.metricsMu.Lock()
	defer sp.metricsMu.Unlock()

	// Update average throughput (simple moving average)
	if sp.metrics.AverageThroughput == 0 {
		sp.metrics.AverageThroughput = throughput
	} else {
		sp.metrics.AverageThroughput = (sp.metrics.AverageThroughput + throughput) / 2
	}

	// Update peak throughput
	if throughput > sp.metrics.PeakThroughput {
		sp.metrics.PeakThroughput = throughput
	}
}

// GetMetrics returns current streaming metrics
func (sp *StreamingProcessor) GetMetrics() StreamingMetrics {
	sp.metricsMu.RLock()
	defer sp.metricsMu.RUnlock()
	return *sp.metrics
}

// DatasetInput represents input parameters for dataset processing
type DatasetInput struct {
	InputPath  string
	OutputDir  string
	SizeGB     float64
	Format     string
	Parameters map[string]interface{}
}

// Helper functions

func generateStreamID() string {
	return fmt.Sprintf("stream-%d", time.Now().UnixNano())
}

// calculateChecksum computes the SHA-256 digest of the file at filePath and
// returns it as a hex string prefixed with "sha256:".  Returns an empty
// string if the file cannot be opened or read.
func calculateChecksum(filePath string) string {
	f, err := os.Open(filePath) // #nosec G304 — path originates from internal StreamingResult
	if err != nil {
		return ""
	}
	defer f.Close()

	h := sha256.New()
	if _, err := io.Copy(h, f); err != nil {
		return ""
	}
	return fmt.Sprintf("sha256:%x", h.Sum(nil))
}

// String methods for enums

func (s StreamStatus) String() string {
	switch s {
	case StatusInitializing:
		return "initializing"
	case StatusProcessing:
		return "processing"
	case StatusCompleted:
		return "completed"
	case StatusFailed:
		return "failed"
	case StatusCancelled:
		return "cancelled"
	default:
		return "unknown"
	}
}
