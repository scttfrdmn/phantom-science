package vcf

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"sync"

	"github.com/scttfrdmn/phantom-science/internal/errors"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// StreamProcessor provides memory-efficient streaming processing of VCF files
type StreamProcessor struct {
	bufferSize       int
	maxCacheSize     int64
	enableCompression bool
	packedArrays     bool
	gpuAcceleration  bool
	numWorkers       int

	// Reference genome for genomic context loading (nil if not configured)
	refManager       *reference.ReferenceManager

	// Internal state
	variantCache     map[string]*VariantRegion
	cacheMutex       sync.RWMutex
	currentCacheSize int64

	// Statistics
	processedVariants int64
	processedBytes    int64
	cacheHits        int64
	cacheMisses      int64
}

// VariantRegion contains cached variant information for a genomic region
type VariantRegion struct {
	Chromosome   string
	StartPos     int64
	EndPos       int64
	Variants     []*Variant
	KmerContext  *kmer.PackedFrequencyTable
	LastAccessed int64
	Size         int64
	Context      []byte
}

// Variant represents a single variant from a VCF file
type Variant struct {
	Chromosome  string
	Position    int64
	ID          string
	Reference   string
	Alternate   []string
	Quality     float64
	Filter      string
	Info        map[string]string
	Format      []string
	Samples     []map[string]string
	Context     []byte  // Surrounding genomic context
}

// StreamConfig provides configuration for VCF streaming
type StreamConfig struct {
	BufferSize          int
	MaxCacheSize        int64   // Maximum cache size in bytes
	EnableCompression   bool    // Enable internal compression
	PackedArrays        bool    // Use packed arrays for k-mer storage
	GPUAcceleration     bool    // Enable GPU processing
	NumWorkers          int     // Number of worker goroutines
	ContextWindowSize   int     // Size of genomic context window
	KmerSize            int     // K-mer size for context analysis
	CacheRegionSize     int64   // Size of genomic regions to cache
	ReferencePath       string  // Path to reference FASTA for genomic context loading
}

// NewStreamProcessor creates a new VCF stream processor
func NewStreamProcessor(config StreamConfig) (*StreamProcessor, error) {
	// Set defaults
	if config.BufferSize <= 0 {
		config.BufferSize = 64 * 1024 // 64KB default buffer
	}
	if config.MaxCacheSize <= 0 {
		config.MaxCacheSize = 1024 * 1024 * 1024 // 1GB default cache
	}
	if config.NumWorkers <= 0 {
		config.NumWorkers = 4 // Default 4 workers
	}
	if config.ContextWindowSize <= 0 {
		config.ContextWindowSize = 1000 // 1KB context window
	}
	if config.KmerSize <= 0 {
		config.KmerSize = 8 // Default k=8
	}
	if config.CacheRegionSize <= 0 {
		config.CacheRegionSize = 1000000 // 1MB regions
	}

	processor := &StreamProcessor{
		bufferSize:        config.BufferSize,
		maxCacheSize:      config.MaxCacheSize,
		enableCompression: config.EnableCompression,
		packedArrays:      config.PackedArrays,
		gpuAcceleration:   config.GPUAcceleration,
		numWorkers:        config.NumWorkers,
		variantCache:      make(map[string]*VariantRegion),
	}

	// Initialize reference manager when a path is configured.
	if config.ReferencePath != "" {
		rm, err := reference.NewReferenceManager(reference.ManagerConfig{
			DefaultReferenceFile: config.ReferencePath,
		})
		if err != nil {
			return nil, fmt.Errorf("failed to initialize reference manager: %w", err)
		}
		processor.refManager = rm
	}

	return processor, nil
}

// ProcessVCFStream processes a VCF file in streaming fashion
func (sp *StreamProcessor) ProcessVCFStream(reader io.Reader, processor VariantProcessor) error {
	// Create buffered reader
	bufReader := bufio.NewReaderSize(reader, sp.bufferSize)
	
	// Check if input is gzipped
	var finalReader io.Reader = bufReader
	if sp.isGzipped(bufReader) {
		gzReader, err := gzip.NewReader(bufReader)
		if err != nil {
			return errors.WrapError(err, errors.VCFParseError, "stream_processor", "process_vcf_stream", "failed to create gzip reader")
		}
		defer gzReader.Close()
		finalReader = bufio.NewReaderSize(gzReader, sp.bufferSize)
	}

	// Parse VCF header
	header, err := sp.parseVCFHeader(finalReader)
	if err != nil {
		return errors.WrapError(err, errors.VCFParseError, "stream_processor", "process_vcf_stream", "failed to parse VCF header")
	}

	// Process variants in streaming fashion
	return sp.streamVariants(finalReader, header, processor)
}

// ProcessVCFFile processes a VCF file by filename
func (sp *StreamProcessor) ProcessVCFFile(filename string, processor VariantProcessor) error {
	file, err := os.Open(filename) // #nosec G304
	if err != nil {
		return errors.WrapError(err, errors.VCFParseError, "stream_processor", "process_vcf_file", fmt.Sprintf("failed to open VCF file %s", filename))
	}
	defer file.Close()

	return sp.ProcessVCFStream(file, processor)
}

// isGzipped checks if the input stream is gzipped
func (sp *StreamProcessor) isGzipped(reader *bufio.Reader) bool {
	magic, err := reader.Peek(2)
	if err != nil || len(magic) < 2 {
		return false
	}
	return magic[0] == 0x1f && magic[1] == 0x8b
}

// parseVCFHeader parses the VCF header section
func (sp *StreamProcessor) parseVCFHeader(reader io.Reader) (*VCFHeader, error) {
	scanner := bufio.NewScanner(reader)
	header := &VCFHeader{
		MetaInfo: make(map[string][]string),
		Info:     make(map[string]*HeaderInfo),
		Format:   make(map[string]*HeaderFormat),
		Filter:   make(map[string]*HeaderFilter),
	}

	for scanner.Scan() {
		line := scanner.Text()
		sp.processedBytes += int64(len(line) + 1) // +1 for newline

		if strings.HasPrefix(line, "##") {
			// Meta-information line
			if err := sp.parseMetaInfoLine(line, header); err != nil {
				return nil, fmt.Errorf("failed to parse meta-info line: %w", err)
			}
		} else if strings.HasPrefix(line, "#CHROM") {
			// Header line with column names
			header.Columns = strings.Split(line[1:], "\t") // Remove leading #
			break
		} else if strings.HasPrefix(line, "#") {
			// Comment line, skip
			continue
		} else {
			return nil, fmt.Errorf("unexpected line before header: %s", line)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading VCF header: %w", err)
	}

	return header, nil
}

// streamVariants processes variants in streaming fashion
func (sp *StreamProcessor) streamVariants(reader io.Reader, header *VCFHeader, processor VariantProcessor) error {
	scanner := bufio.NewScanner(reader)
	
	// Create worker pool for parallel processing
	variantChan := make(chan *Variant, sp.numWorkers*2)
	errorChan := make(chan error, sp.numWorkers)
	var wg sync.WaitGroup

	// Start worker goroutines
	for i := 0; i < sp.numWorkers; i++ {
		wg.Add(1)
		go sp.variantWorker(&wg, variantChan, errorChan, processor)
	}

	// Parse and distribute variants
	go func() {
		defer close(variantChan)
		
		for scanner.Scan() {
			line := scanner.Text()
			sp.processedBytes += int64(len(line) + 1)

			// Skip comment lines
			if strings.HasPrefix(line, "#") {
				continue
			}

			// Parse variant
			variant, err := sp.parseVariantLine(line, header)
			if err != nil {
				errorChan <- fmt.Errorf("failed to parse variant line: %w", err)
				return
			}

			// Send variant to worker pool
			select {
			case variantChan <- variant:
				sp.processedVariants++
			case err := <-errorChan:
				errorChan <- err
				return
			}
		}

		if err := scanner.Err(); err != nil {
			errorChan <- fmt.Errorf("error reading VCF variants: %w", err)
		}
	}()

	// Wait for workers to complete
	wg.Wait()
	close(errorChan)

	// Check for errors
	select {
	case err := <-errorChan:
		return err
	default:
		return nil
	}
}

// variantWorker processes variants in parallel
func (sp *StreamProcessor) variantWorker(wg *sync.WaitGroup, variantChan <-chan *Variant, errorChan chan<- error, processor VariantProcessor) {
	defer wg.Done()

	for variant := range variantChan {
		// Add genomic context
		if err := sp.addGenomicContext(variant); err != nil {
			errorChan <- fmt.Errorf("failed to add genomic context: %w", err)
			return
		}

		// Process variant
		if err := processor.ProcessVariant(variant); err != nil {
			errorChan <- fmt.Errorf("failed to process variant: %w", err)
			return
		}
	}
}

// parseMetaInfoLine parses a VCF meta-information line
func (sp *StreamProcessor) parseMetaInfoLine(line string, header *VCFHeader) error {
	// Remove ## prefix
	line = line[2:]

	// Find key-value separator
	eqIndex := strings.Index(line, "=")
	if eqIndex == -1 {
		return fmt.Errorf("invalid meta-info line format: %s", line)
	}

	key := line[:eqIndex]
	value := line[eqIndex+1:]

	// Parse structured meta-info (INFO, FORMAT, FILTER, etc.)
	switch key {
	case "INFO":
		info, err := sp.parseInfoHeader(value)
		if err != nil {
			return err
		}
		header.Info[info.ID] = info
	case "FORMAT":
		format, err := sp.parseFormatHeader(value)
		if err != nil {
			return err
		}
		header.Format[format.ID] = format
	case "FILTER":
		filter, err := sp.parseFilterHeader(value)
		if err != nil {
			return err
		}
		header.Filter[filter.ID] = filter
	default:
		// Generic meta-info
		if header.MetaInfo[key] == nil {
			header.MetaInfo[key] = make([]string, 0)
		}
		header.MetaInfo[key] = append(header.MetaInfo[key], value)
	}

	return nil
}

// parseVariantLine parses a single variant line
func (sp *StreamProcessor) parseVariantLine(line string, header *VCFHeader) (*Variant, error) {
	fields := strings.Split(line, "\t")
	if len(fields) < 8 {
		return nil, fmt.Errorf("insufficient fields in variant line: %d", len(fields))
	}

	// Parse position
	pos, err := strconv.ParseInt(fields[1], 10, 64)
	if err != nil {
		return nil, fmt.Errorf("invalid position: %s", fields[1])
	}

	// Parse quality
	var quality float64
	if fields[5] != "." {
		quality, err = strconv.ParseFloat(fields[5], 64)
		if err != nil {
			return nil, fmt.Errorf("invalid quality: %s", fields[5])
		}
	}

	// Parse alternate alleles
	var alternates []string
	if fields[4] != "." {
		alternates = strings.Split(fields[4], ",")
	}

	// Parse INFO field
	info := make(map[string]string)
	if fields[7] != "." {
		infoPairs := strings.Split(fields[7], ";")
		for _, pair := range infoPairs {
			if strings.Contains(pair, "=") {
				kv := strings.SplitN(pair, "=", 2)
				info[kv[0]] = kv[1]
			} else {
				info[pair] = "true" // Flag fields
			}
		}
	}

	variant := &Variant{
		Chromosome: fields[0],
		Position:   pos,
		ID:         fields[2],
		Reference:  fields[3],
		Alternate:  alternates,
		Quality:    quality,
		Filter:     fields[6],
		Info:       info,
	}

	// Parse FORMAT and sample data if present
	if len(fields) > 8 {
		formatFields := strings.Split(fields[8], ":")
		variant.Format = formatFields
		variant.Samples = make([]map[string]string, len(fields)-9)

		for i, sampleData := range fields[9:] {
			sampleValues := strings.Split(sampleData, ":")
			sampleMap := make(map[string]string)
			
			for j, value := range sampleValues {
				if j < len(formatFields) {
					sampleMap[formatFields[j]] = value
				}
			}
			variant.Samples[i] = sampleMap
		}
	}

	return variant, nil
}

// addGenomicContext adds surrounding genomic context to a variant.
// When no reference genome is configured (StreamConfig.ReferencePath is empty),
// context loading is skipped and variant.Context remains nil.
func (sp *StreamProcessor) addGenomicContext(variant *Variant) error {
	if sp.refManager == nil {
		// No reference configured; skip context loading.
		return nil
	}

	regionKey := fmt.Sprintf("%s:%d-%d", variant.Chromosome,
		variant.Position-500, variant.Position+500)

	// Check cache first
	sp.cacheMutex.RLock()
	_, found := sp.variantCache[regionKey]
	sp.cacheMutex.RUnlock()

	if found {
		sp.cacheHits++
		// Region metadata is cached; load context from reference for this position.
		context, err := sp.loadGenomicContext(variant.Chromosome, variant.Position-500, variant.Position+500)
		if err != nil {
			return fmt.Errorf("failed to load genomic context (cached region): %w", err)
		}
		variant.Context = context
		return nil
	}

	sp.cacheMisses++

	context, err := sp.loadGenomicContext(variant.Chromosome, variant.Position-500, variant.Position+500)
	if err != nil {
		return fmt.Errorf("failed to load genomic context: %w", err)
	}

	variant.Context = context

	// Cache the region metadata for future lookups.
	sp.cacheRegion(regionKey, variant.Chromosome, variant.Position-500, variant.Position+500, context)

	return nil
}

// loadGenomicContext loads the genomic context window from the configured reference genome.
// Returns an error when no reference genome is configured (StreamConfig.ReferencePath is empty).
func (sp *StreamProcessor) loadGenomicContext(chromosome string, start, end int64) ([]byte, error) {
	if sp.refManager == nil {
		return nil, fmt.Errorf("reference genome not configured; set StreamConfig.ReferencePath to enable genomic context loading")
	}
	return sp.refManager.GetSequence(chromosome, start, end)
}

// extractContextFromRegion extracts genomic context from cached region.
// Returns a slice of the cached Context bytes centred on the variant position,
// up to 500 bases on each side. Returns nil if no context was cached.
func (sp *StreamProcessor) extractContextFromRegion(variant *Variant, region *VariantRegion) []byte {
	if len(region.Context) == 0 {
		return nil
	}
	relativePos := int(variant.Position - region.StartPos)
	const radius = 500
	start := relativePos - radius
	end := relativePos + radius
	if start < 0 {
		start = 0
	}
	if end > len(region.Context) {
		end = len(region.Context)
	}
	if start >= end {
		return nil
	}
	return region.Context[start:end]
}

// cacheRegion caches a genomic region
func (sp *StreamProcessor) cacheRegion(regionKey, chromosome string, start, end int64, context []byte) {
	region := &VariantRegion{
		Chromosome:   chromosome,
		StartPos:     start,
		EndPos:       end,
		LastAccessed: sp.processedVariants,
		Size:         int64(len(context)),
		Context:      context,
	}

	sp.cacheMutex.Lock()
	defer sp.cacheMutex.Unlock()

	// Check cache size and evict if necessary
	sp.currentCacheSize += region.Size
	if sp.currentCacheSize > sp.maxCacheSize {
		sp.evictOldestRegions()
	}

	sp.variantCache[regionKey] = region
}

// evictOldestRegions removes oldest cached regions
func (sp *StreamProcessor) evictOldestRegions() {
	// Simple LRU eviction based on LastAccessed
	type regionEntry struct {
		key    string
		region *VariantRegion
	}

	var regions []regionEntry
	for key, region := range sp.variantCache {
		regions = append(regions, regionEntry{key, region})
	}

	// Sort by LastAccessed (ascending)
	for i := 0; i < len(regions)-1; i++ {
		for j := i + 1; j < len(regions); j++ {
			if regions[i].region.LastAccessed > regions[j].region.LastAccessed {
				regions[i], regions[j] = regions[j], regions[i]
			}
		}
	}

	// Remove oldest regions until under cache limit
	for _, entry := range regions {
		if sp.currentCacheSize <= sp.maxCacheSize/2 { // Evict to 50% of limit
			break
		}
		
		delete(sp.variantCache, entry.key)
		sp.currentCacheSize -= entry.region.Size
	}
}

// GetStatistics returns processing statistics
func (sp *StreamProcessor) GetStatistics() ProcessingStatistics {
	sp.cacheMutex.RLock()
	defer sp.cacheMutex.RUnlock()

	cacheHitRatio := 0.0
	if sp.cacheHits+sp.cacheMisses > 0 {
		cacheHitRatio = float64(sp.cacheHits) / float64(sp.cacheHits+sp.cacheMisses)
	}
	cacheUtil := 0.0
	if sp.maxCacheSize > 0 {
		cacheUtil = float64(sp.currentCacheSize) / float64(sp.maxCacheSize)
	}

	return ProcessingStatistics{
		ProcessedVariants: sp.processedVariants,
		ProcessedBytes:    sp.processedBytes,
		CacheHits:        sp.cacheHits,
		CacheMisses:      sp.cacheMisses,
		CacheHitRatio:    cacheHitRatio,
		CachedRegions:    int64(len(sp.variantCache)),
		CacheSize:        sp.currentCacheSize,
		CacheUtilization: cacheUtil,
	}
}

// ProcessingStatistics contains VCF processing statistics
type ProcessingStatistics struct {
	ProcessedVariants int64
	ProcessedBytes    int64
	CacheHits        int64
	CacheMisses      int64
	CacheHitRatio    float64
	CachedRegions    int64
	CacheSize        int64
	CacheUtilization float64
}

// VariantProcessor interface for processing individual variants
type VariantProcessor interface {
	ProcessVariant(variant *Variant) error
}

// VCFHeader represents the header section of a VCF file
type VCFHeader struct {
	MetaInfo map[string][]string
	Info     map[string]*HeaderInfo
	Format   map[string]*HeaderFormat
	Filter   map[string]*HeaderFilter
	Columns  []string
}

// HeaderInfo represents an INFO header field
type HeaderInfo struct {
	ID          string
	Number      string
	Type        string
	Description string
}

// HeaderFormat represents a FORMAT header field
type HeaderFormat struct {
	ID          string
	Number      string
	Type        string
	Description string
}

// HeaderFilter represents a FILTER header field
type HeaderFilter struct {
	ID          string
	Description string
}

// Helper functions for parsing structured headers
func (sp *StreamProcessor) parseInfoHeader(value string) (*HeaderInfo, error) {
	// Parse <ID=...,Number=...,Type=...,Description="...">
	result, err := sp.parseStructuredHeader(value, &HeaderInfo{})
	if err != nil {
		return nil, err
	}
	return result.(*HeaderInfo), nil
}

func (sp *StreamProcessor) parseFormatHeader(value string) (*HeaderFormat, error) {
	// Parse <ID=...,Number=...,Type=...,Description="...">
	result, err := sp.parseStructuredHeader(value, &HeaderFormat{})
	if err != nil {
		return nil, err
	}
	return result.(*HeaderFormat), nil
}

func (sp *StreamProcessor) parseFilterHeader(value string) (*HeaderFilter, error) {
	// Parse <ID=...,Description="...">
	result, err := sp.parseStructuredHeader(value, &HeaderFilter{})
	if err != nil {
		return nil, err
	}
	return result.(*HeaderFilter), nil
}

func (sp *StreamProcessor) parseStructuredHeader(value string, target interface{}) (interface{}, error) {
	// Simple structured header parser (would be more robust in production)
	if !strings.HasPrefix(value, "<") || !strings.HasSuffix(value, ">") {
		return nil, fmt.Errorf("invalid structured header format: %s", value)
	}

	// Remove angle brackets
	content := value[1 : len(value)-1]
	
	// For now, return a basic parsed structure
	// In production, this would use proper VCF header parsing
	switch target.(type) {
	case *HeaderInfo:
		return &HeaderInfo{ID: "placeholder", Type: "String", Description: content}, nil
	case *HeaderFormat:
		return &HeaderFormat{ID: "placeholder", Type: "String", Description: content}, nil
	case *HeaderFilter:
		return &HeaderFilter{ID: "placeholder", Description: content}, nil
	default:
		return nil, fmt.Errorf("unknown header type")
	}
}