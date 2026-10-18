// Package reference reads reference genome FASTA files and counts k-mers.
//
// FastaParser (plain or gzip FASTA), ReferenceGenomeParser and
// SimpleReferenceGenome (plain FASTA) read real sequence. LazyReferenceLoader and ReferenceManager
// are caches whose loader is a stub: they return placeholder sequence
// generated from the coordinates, not the contents of the reference file.
package reference

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
)

// FastaParser reads plain or gzip-compressed FASTA files and builds, loads
// and saves .fai indexes. Plain files support indexed random access;
// compressed files are scanned linearly.
type FastaParser struct {
	filename        string
	compressionType CompressionType
	indexFile       string

	// Cached file handle for streaming access
	file    *os.File
	reader  io.Reader
	scanner *bufio.Scanner
	isOpen  bool
}

// CompressionType represents different compression formats
type CompressionType int

const (
	NoCompression CompressionType = iota
	GzipCompression
	BgzfCompression // BGZF (samtools); never detected: BGZF files are reported as GzipCompression
)

// ChromosomeInfo contains metadata about a chromosome in the reference
type ChromosomeInfo struct {
	Name         string
	Length       int64
	StartOffset  int64 // Byte offset in file where sequence starts
	LineLength   int   // Bases per line
	BytesPerLine int   // Bytes per line (including newline)
}

// FastaIndex represents a FASTA index (.fai file) for rapid random access
type FastaIndex struct {
	Chromosomes  map[string]*ChromosomeInfo
	OrderedNames []string // Maintain chromosome order
}

// NewFastaParser creates a new FASTA parser for the specified file
func NewFastaParser(filename string) (*FastaParser, error) {
	// Detect compression type
	compressionType, err := detectCompressionType(filename)
	if err != nil {
		return nil, fmt.Errorf("failed to detect compression type: %w", err)
	}

	parser := &FastaParser{
		filename:        filename,
		compressionType: compressionType,
		indexFile:       filename + ".fai",
	}

	return parser, nil
}

// detectCompressionType determines the compression format of the file
func detectCompressionType(filename string) (CompressionType, error) {
	file, err := os.Open(filename) // #nosec G304
	if err != nil {
		return NoCompression, err
	}
	defer file.Close()

	// Read first few bytes to check magic numbers
	magic := make([]byte, 4)
	n, err := file.Read(magic)
	if err != nil && err != io.EOF {
		return NoCompression, err
	}

	if n >= 2 {
		// Check for gzip magic number
		if magic[0] == 0x1f && magic[1] == 0x8b {
			return GzipCompression, nil
		}
	}

	// Check file extension as fallback
	if strings.HasSuffix(strings.ToLower(filename), ".gz") {
		return GzipCompression, nil
	}

	return NoCompression, nil
}

// Open opens the FASTA file for reading
func (fp *FastaParser) Open() error {
	if fp.isOpen {
		return nil
	}

	file, err := os.Open(fp.filename)
	if err != nil {
		return fmt.Errorf("failed to open FASTA file %s: %w", fp.filename, err)
	}

	fp.file = file

	// Set up appropriate reader based on compression
	switch fp.compressionType {
	case GzipCompression:
		gzReader, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return fmt.Errorf("failed to create gzip reader: %w", err)
		}
		fp.reader = gzReader
	default:
		fp.reader = file
	}

	fp.scanner = bufio.NewScanner(fp.reader)
	fp.isOpen = true

	return nil
}

// Close closes the FASTA file
func (fp *FastaParser) Close() error {
	if !fp.isOpen {
		return nil
	}

	var err error

	// Close gzip reader if present
	if gzReader, ok := fp.reader.(*gzip.Reader); ok {
		err = gzReader.Close()
	}

	// Close file handle
	if fp.file != nil {
		if fileErr := fp.file.Close(); fileErr != nil && err == nil {
			err = fileErr
		}
	}

	fp.isOpen = false
	fp.file = nil
	fp.reader = nil
	fp.scanner = nil

	return err
}

// BuildIndex creates a FASTA index for rapid random access
func (fp *FastaParser) BuildIndex() (*FastaIndex, error) {
	if err := fp.Open(); err != nil {
		return nil, err
	}
	defer fp.Close()

	index := &FastaIndex{
		Chromosomes:  make(map[string]*ChromosomeInfo),
		OrderedNames: make([]string, 0),
	}

	var currentChrom *ChromosomeInfo
	var lineCount int64
	var byteOffset int64

	for fp.scanner.Scan() {
		line := fp.scanner.Text()
		lineBytes := int64(len(line) + 1) // +1 for newline

		if strings.HasPrefix(line, ">") {
			// Header line - new chromosome
			if currentChrom != nil {
				// Finalize previous chromosome
				index.Chromosomes[currentChrom.Name] = currentChrom
			}

			fields := strings.Fields(line[1:]) // Remove > and split on whitespace
			if len(fields) == 0 {
				// Malformed header with no chromosome name — skip, but still
				// advance byteOffset so subsequent lines have correct offsets.
				byteOffset += lineBytes
				continue
			}
			chromName := fields[0]

			currentChrom = &ChromosomeInfo{
				Name:        chromName,
				StartOffset: byteOffset + lineBytes,
				Length:      0,
			}

			index.OrderedNames = append(index.OrderedNames, chromName)
			lineCount = 0

		} else if currentChrom != nil {
			// Sequence line
			seqLen := int64(len(line))
			currentChrom.Length += seqLen
			lineCount++

			// Set line parameters from first sequence line
			if lineCount == 1 {
				currentChrom.LineLength = len(line)
				currentChrom.BytesPerLine = len(line) + 1 // +1 for newline
			}
		}

		byteOffset += lineBytes
	}

	// Finalize last chromosome
	if currentChrom != nil {
		index.Chromosomes[currentChrom.Name] = currentChrom
	}

	if err := fp.scanner.Err(); err != nil {
		return nil, fmt.Errorf("error scanning FASTA file: %w", err)
	}

	return index, nil
}

// LoadIndex loads a FASTA index from a .fai file
func (fp *FastaParser) LoadIndex() (*FastaIndex, error) {
	indexFile, err := os.Open(fp.indexFile)
	if err != nil {
		// Index file doesn't exist - build it
		return fp.BuildIndex()
	}
	defer indexFile.Close()

	index := &FastaIndex{
		Chromosomes:  make(map[string]*ChromosomeInfo),
		OrderedNames: make([]string, 0),
	}

	scanner := bufio.NewScanner(indexFile)
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")

		if len(fields) < 5 {
			continue
		}

		chromInfo := &ChromosomeInfo{
			Name: fields[0],
		}

		// Parse numeric fields
		if length, err := parseUint64(fields[1]); err == nil {
			chromInfo.Length = int64(length) // #nosec G115
		}
		if offset, err := parseUint64(fields[2]); err == nil {
			chromInfo.StartOffset = int64(offset) // #nosec G115
		}
		if lineLength, err := parseInt(fields[3]); err == nil {
			chromInfo.LineLength = lineLength
		}
		if bytesPerLine, err := parseInt(fields[4]); err == nil {
			chromInfo.BytesPerLine = bytesPerLine
		}

		index.Chromosomes[chromInfo.Name] = chromInfo
		index.OrderedNames = append(index.OrderedNames, chromInfo.Name)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading FASTA index: %w", err)
	}

	return index, nil
}

// SaveIndex saves a FASTA index to a .fai file
func (fp *FastaParser) SaveIndex(index *FastaIndex) error {
	indexFile, err := os.Create(fp.indexFile)
	if err != nil {
		return fmt.Errorf("failed to create index file %s: %w", fp.indexFile, err)
	}
	defer indexFile.Close()

	writer := bufio.NewWriter(indexFile)
	defer writer.Flush()

	for _, chromName := range index.OrderedNames {
		chromInfo := index.Chromosomes[chromName]
		fmt.Fprintf(writer, "%s\t%d\t%d\t%d\t%d\n",
			chromInfo.Name,
			chromInfo.Length,
			chromInfo.StartOffset,
			chromInfo.LineLength,
			chromInfo.BytesPerLine)
	}

	return nil
}

// GetSequence extracts sequence for a specific genomic region using the index
func (fp *FastaParser) GetSequence(index *FastaIndex, chromosome string, start, end int64) ([]byte, error) {
	chromInfo, exists := index.Chromosomes[chromosome]
	if !exists {
		return nil, fmt.Errorf("chromosome %s not found in reference", chromosome)
	}

	// Validate coordinates
	if start < 0 || end > chromInfo.Length || start >= end {
		return nil, fmt.Errorf("invalid coordinates: %s:%d-%d (chr length: %d)",
			chromosome, start, end, chromInfo.Length)
	}

	// Compressed files are decompressed and scanned from the beginning
	// (no BGZF/.gzi random access)
	if fp.compressionType != NoCompression {
		return fp.getSequenceLinearScan(chromosome, start, end)
	}

	// Calculate file positions using index
	seqLength := end - start
	sequence := make([]byte, 0, seqLength)

	file, err := os.Open(fp.filename)
	if err != nil {
		return nil, fmt.Errorf("failed to open reference file: %w", err)
	}
	defer file.Close()

	// Calculate starting position in file
	lineNum := start / int64(chromInfo.LineLength)
	posInLine := start % int64(chromInfo.LineLength)
	fileOffset := chromInfo.StartOffset + lineNum*int64(chromInfo.BytesPerLine) + posInLine

	// Seek to position
	if _, err := file.Seek(fileOffset, io.SeekStart); err != nil {
		return nil, fmt.Errorf("failed to seek to position: %w", err)
	}

	reader := bufio.NewReader(file)
	remaining := seqLength

	for remaining > 0 {
		// Read characters until newline or remaining count
		line, err := reader.ReadString('\n')
		if err != nil && err != io.EOF {
			return nil, fmt.Errorf("error reading sequence: %w", err)
		}

		// Remove newline
		line = strings.TrimSuffix(line, "\n")
		line = strings.TrimSuffix(line, "\r")

		// Take only what we need
		takeCount := int64(len(line))
		if takeCount > remaining {
			takeCount = remaining
		}

		sequence = append(sequence, line[:takeCount]...)
		remaining -= takeCount

		if err == io.EOF {
			break
		}
	}

	return sequence, nil
}

// getSequenceLinearScan extracts sequence by scanning through the file (for compressed files)
func (fp *FastaParser) getSequenceLinearScan(chromosome string, start, end int64) ([]byte, error) {
	if err := fp.Open(); err != nil {
		return nil, err
	}
	defer fp.Close()

	var currentChrom string
	var currentPos int64
	var inTargetChrom bool
	sequence := make([]byte, 0, end-start)

	for fp.scanner.Scan() {
		line := fp.scanner.Text()

		if strings.HasPrefix(line, ">") {
			// Header line
			headerFields := strings.Fields(line[1:])
			if len(headerFields) == 0 {
				continue // skip malformed header with no chromosome name
			}
			currentChrom = headerFields[0]
			currentPos = 0
			inTargetChrom = (currentChrom == chromosome)

			// If we've passed our target chromosome, we're done
			if inTargetChrom && len(sequence) > 0 {
				break
			}

		} else if inTargetChrom {
			// Sequence line in target chromosome
			lineSeq := []byte(line)
			lineLength := int64(len(lineSeq))

			// Check if this line contains our target region
			lineStart := currentPos
			lineEnd := currentPos + lineLength

			if lineEnd <= start {
				// Line is before our region
				currentPos = lineEnd
				continue
			}

			if lineStart >= end {
				// Line is after our region - we're done
				break
			}

			// Line overlaps with our region
			takeStart := int64(0)
			takeEnd := lineLength

			if lineStart < start {
				takeStart = start - lineStart
			}

			if lineEnd > end {
				takeEnd = end - lineStart
			}

			sequence = append(sequence, lineSeq[takeStart:takeEnd]...)
			currentPos = lineEnd

			// Check if we have enough sequence
			if int64(len(sequence)) >= end-start {
				break
			}
		}
	}

	if err := fp.scanner.Err(); err != nil {
		return nil, fmt.Errorf("error scanning FASTA file: %w", err)
	}

	return sequence, nil
}

// GetChromosomeNames returns all chromosome names in order
func (fp *FastaParser) GetChromosomeNames(index *FastaIndex) []string {
	return index.OrderedNames
}

// GetChromosomeLength returns the length of a specific chromosome
func (fp *FastaParser) GetChromosomeLength(index *FastaIndex, chromosome string) (int64, error) {
	chromInfo, exists := index.Chromosomes[chromosome]
	if !exists {
		return 0, fmt.Errorf("chromosome %s not found", chromosome)
	}
	return chromInfo.Length, nil
}

// Helper functions for parsing
func parseUint64(s string) (uint64, error) {
	var result uint64
	for _, c := range s {
		if c < '0' || c > '9' {
			return 0, fmt.Errorf("invalid number: %s", s)
		}
		result = result*10 + uint64(c-'0')
	}
	return result, nil
}

func parseInt(s string) (int, error) {
	val, err := parseUint64(s)
	return int(val), err // #nosec G115
}
