// Package genomics provides the core FASTQ generation engine, including HMM model
// integration, steganographic fingerprinting, and sequencer-realistic output.
package genomics

import (
    "bufio"
    "bytes"
    "crypto/rand"
    "fmt"
    "io"
    "os"
    "path/filepath"
    "sync"
    "github.com/scttfrdmn/phantom-science/internal/fingerprinting"
    "github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
    "github.com/scttfrdmn/phantom-science/internal/security"
    "github.com/scttfrdmn/phantom-science/pkg/types"
)

type Engine struct {
    config               types.SequencerConfig
    fingerprintingEngine fingerprinting.Engine
    hmmModel             authenticity.SequenceGenerator // nil = use crypto/rand

    // Files written by Generate and GenerateFASTQFile, in order.
    filesMu sync.Mutex
    files   []string

    // keyErr records a fingerprint key that could not be loaded; encrypted
    // reports whether quality-score fingerprints are encrypted with a key.
    keyErr    error
    encrypted bool

    // One fingerprint per engine, shared by every file it writes.
    fingerprintOnce sync.Once
    fingerprint     string
    fingerprintErr  error
}

// qualityLayer is the steganography layer that embeds a fingerprint in FASTQ
// quality scores.
const qualityLayer = "quality_scores"

// qualityFingerprintingEnabled reports whether output files must carry a
// quality-score fingerprint.
func (e *Engine) qualityFingerprintingEnabled() bool {
    if !e.config.Fingerprinting.Enable {
        return false
    }
    for _, layer := range e.config.Fingerprinting.SteganographyLayers {
        if layer == qualityLayer {
            return true
        }
    }
    return false
}

// fileFingerprinter returns the engine's file-level fingerprinter, or an error
// when the fingerprinting engine cannot embed into files (CGO disabled, release
// build, or stegano library unavailable).
func (e *Engine) fileFingerprinter() (fingerprinting.FileFingerprinter, error) {
    if ff, ok := e.fingerprintingEngine.(fingerprinting.FileFingerprinter); ok {
        return ff, nil
    }
    return nil, fmt.Errorf("steganography layer %q needs the stegano library: build with CGO after `make stegano`, or remove %q from steganography_layers", qualityLayer, qualityLayer)
}

// checkQualityFingerprinting reports why quality-score fingerprints cannot be
// written, or nil when they can: the stegano library must be available, and
// the fingerprint must be encrypted with a key unless plaintext is allowed.
func (e *Engine) checkQualityFingerprinting() error {
    if _, err := e.fileFingerprinter(); err != nil {
        return err
    }
    if e.keyErr != nil {
        return e.keyErr
    }
    if !e.encrypted && !e.config.Fingerprinting.Plaintext {
        return fmt.Errorf("steganography layer %q needs an encryption key: set fingerprinting.key_file or %s, or set plaintext: true to embed readable fingerprints", qualityLayer, fingerprinting.KeyEnvVar)
    }
    return nil
}

// runFingerprint returns the fingerprint for this engine's run, generating it
// on first use.
func (e *Engine) runFingerprint() (string, error) {
    e.fingerprintOnce.Do(func() {
        e.fingerprint, e.fingerprintErr = e.fingerprintingEngine.GenerateFingerprint(e.config.DatasetID)
    })
    return e.fingerprint, e.fingerprintErr
}

// embedFileFingerprint embeds the run fingerprint into the finished FASTQ file
// at filename. On failure the unmarked file is removed, so output never
// silently lacks a requested fingerprint.
func (e *Engine) embedFileFingerprint(filename string) error {
    ff, err := e.fileFingerprinter()
    if err == nil {
        var fp string
        if fp, err = e.runFingerprint(); err == nil {
            err = ff.EmbedInFASTQFile(filename, filename, fp)
        }
    }
    if err != nil {
        os.Remove(filename) // #nosec G104 -- best effort; the embed error is returned
        return fmt.Errorf("fingerprinting %s: %w", filepath.Base(filename), err)
    }
    return nil
}

func NewEngine(config types.SequencerConfig) *Engine {
    // Initialize fingerprinting engine with customer ID from config
    customerID := "phantom-science" // Default customer ID
    if config.Enterprise.CustomerID != "" {
        customerID = config.Enterprise.CustomerID
    }

    e := &Engine{config: config}
    var fpOpts []fingerprinting.Option
    if config.Fingerprinting.Enable {
        key, err := fingerprinting.ResolveKey(config.Fingerprinting.KeyFile)
        if err != nil {
            e.keyErr = err
        } else if key != nil {
            fpOpts = append(fpOpts, fingerprinting.WithKey(key))
            e.encrypted = true
        }
    }
    e.fingerprintingEngine = fingerprinting.New(customerID, fpOpts...)
    if modelPath := resolveModelPath(config); modelPath != "" {
        m, err := authenticity.LoadModel(modelPath)
        if err != nil {
            fmt.Printf("   [warn] HMM model load failed (%s): %v — using crypto/rand\n",
                modelPath, err)
        } else {
            e.hmmModel = m
            fmt.Printf("   [info] Loaded HMM model from %s\n", modelPath)
        }
    }
    return e
}

// resolveModelPath returns the effective HMM model file path for config.
// Resolution order:
//  1. Explicit HMMModelPath — highest priority.
//  2. ModelDir/<Population>.hmm.json — when both Population and ModelDir are set.
//  3. Empty string — engine uses crypto/rand.
func resolveModelPath(config types.SequencerConfig) string {
    if config.HMMModelPath != "" {
        return config.HMMModelPath
    }
    if config.Population != "" && config.ModelDir != "" {
        return filepath.Join(config.ModelDir, config.Population+".hmm.json")
    }
    return ""
}

func (e *Engine) Generate() error {
    // Validate configuration
    if e.config.ReadLength <= 0 {
        return fmt.Errorf("invalid read length: must be greater than 0, got %d", e.config.ReadLength)
    }
    if e.config.ReadsPerBatch <= 0 {
        return fmt.Errorf("invalid reads per batch: must be greater than 0, got %d", e.config.ReadsPerBatch)
    }
    if e.config.TotalBatches <= 0 {
        return fmt.Errorf("invalid total batches: must be greater than 0, got %d", e.config.TotalBatches)
    }
    if e.config.SampleName == "" {
        return fmt.Errorf("sample name cannot be empty")
    }
    if e.config.OutputDir == "" {
        return fmt.Errorf("output directory cannot be empty")
    }
    if e.qualityFingerprintingEnabled() {
        if err := e.checkQualityFingerprinting(); err != nil {
            return err
        }
    }
    
    fmt.Println("🧬 Initializing genomics simulation engine...")
    fmt.Printf("   Sample: %s\n", e.config.SampleName)
    fmt.Printf("   Reads per batch: %d\n", e.config.ReadsPerBatch)
    fmt.Printf("   Total batches: %d\n", e.config.TotalBatches)
    
    // Note: crypto/rand doesn't need seeding
    
    // Create output directory with security validation
    if err := security.SecureMkdirAll(e.config.OutputDir, ""); err != nil {
        return fmt.Errorf("failed to create output directory: %w", err)
    }
    
    totalReads := e.config.TotalBatches * e.config.ReadsPerBatch
    fmt.Printf("   Total reads to generate: %d\n", totalReads)
    
    // Generate FASTQ files for each batch
    for batch := 1; batch <= e.config.TotalBatches; batch++ {
        fmt.Printf("   Generating batch %d/%d...\n", batch, e.config.TotalBatches)
        
        if e.config.PairedEnd {
            // Generate paired-end files
            r1File := filepath.Join(e.config.OutputDir, fmt.Sprintf("%s_S1_L%03d_R1_%03d.fastq", e.config.SampleName, e.config.Lane, batch))
            r2File := filepath.Join(e.config.OutputDir, fmt.Sprintf("%s_S1_L%03d_R2_%03d.fastq", e.config.SampleName, e.config.Lane, batch))
            
            if err := e.generateFASTQFile(r1File, batch, 1); err != nil {
                return fmt.Errorf("failed to generate R1 file: %w", err)
            }
            if err := e.generateFASTQFile(r2File, batch, 2); err != nil {
                return fmt.Errorf("failed to generate R2 file: %w", err)
            }
        } else {
            // Generate single-end file
            fastqFile := filepath.Join(e.config.OutputDir, fmt.Sprintf("%s_S1_L%03d_%03d.fastq", e.config.SampleName, e.config.Lane, batch))
            if err := e.generateFASTQFile(fastqFile, batch, 1); err != nil {
                return fmt.Errorf("failed to generate FASTQ file: %w", err)
            }
        }
    }
    
    fmt.Printf("   ✅ Generated %d FASTQ files in %s\n", 
        e.config.TotalBatches * (1 + boolToInt(e.config.PairedEnd)), 
        e.config.OutputDir)
    
    return nil
}

// Close cleans up engine resources
func (e *Engine) Close() {
    if e.fingerprintingEngine != nil {
        e.fingerprintingEngine.Close()
    }
}

// GenerateFASTQFile generates a single FASTQ file (public method for TUI)
func (e *Engine) GenerateFASTQFile(filename string, batch, read int) error {
    return e.generateFASTQFile(filename, batch, read)
}

func (e *Engine) generateFASTQFile(filename string, batch, read int) error {
    // Validate file path against the configured output directory
    if err := security.ValidateFilePath(filename, e.config.OutputDir); err != nil {
        return fmt.Errorf("invalid file path: %w", err)
    }
    
    if e.qualityFingerprintingEnabled() {
        if err := e.checkQualityFingerprinting(); err != nil {
            return err
        }
    }

    file, err := security.SecureFileCreate(filename, e.config.OutputDir)
    if err != nil {
        return err
    }
    if err := e.writeFASTQRecords(file, batch, read); err != nil {
        file.Close() // #nosec G104 -- the write error is returned
        os.Remove(filename) // #nosec G104 -- best effort cleanup of a partial file
        return err
    }
    if err := file.Close(); err != nil {
        return fmt.Errorf("closing %s: %w", filepath.Base(filename), err)
    }
    if e.qualityFingerprintingEnabled() {
        if err := e.embedFileFingerprint(filename); err != nil {
            return err
        }
    }
    e.filesMu.Lock()
    e.files = append(e.files, filename)
    e.filesMu.Unlock()
    return nil
}

// GeneratedFiles returns the FASTQ files this engine has written, in order.
func (e *Engine) GeneratedFiles() []string {
    e.filesMu.Lock()
    defer e.filesMu.Unlock()
    return append([]string(nil), e.files...)
}

// writeFASTQRecords writes one batch of FASTQ records to w.
func (e *Engine) writeFASTQRecords(w io.Writer, batch, read int) error {
    bw := bufio.NewWriterSize(w, 256*1024)
    qualityChars := "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
    
    for readNum := 1; readNum <= e.config.ReadsPerBatch; readNum++ {
        // Generate read ID with batch and read info
        baseReadID := fmt.Sprintf("@%s:%d:HXXXXXX:%d:1101:%d:%d %d:N:0:1", 
            e.config.SampleName, e.config.Lane, batch, readNum*1000, readNum*2000, read)
        
        // Apply fingerprinting to read ID if enabled
        readID := baseReadID
        if e.config.Fingerprinting.Enable && e.fingerprintingEngine != nil {
            if e.config.Fingerprinting.CustomerFingerprint {
                readID = e.fingerprintingEngine.EmbedCustomerID(baseReadID)
            }
            if e.config.Fingerprinting.DatasetFingerprint && e.config.DatasetID != "" {
                readID = e.fingerprintingEngine.EmbedDatasetFingerprint(readID, e.config.DatasetID)
            }
        }
        
        // Generate DNA sequence (HMM model if configured, else crypto/rand).
        sequence, hmmQuals, err := e.generateSequence(e.config.ReadLength)
        if err != nil {
            return fmt.Errorf("failed to generate sequence: %w", err)
        }

        // Encode quality scores.
        qualities := make([]byte, e.config.ReadLength)
        if hmmQuals != nil {
            // HMM emission model provides per-base Phred scores (0–41).
            for i, q := range hmmQuals {
                idx := int(q)
                if idx >= len(qualityChars) {
                    idx = len(qualityChars) - 1
                }
                qualities[i] = qualityChars[idx]
            }
        } else {
            // Config-based quality decay (original path).
            for i := 0; i < e.config.ReadLength; i++ {
                qualityValue := e.config.BaseQuality - int(float64(i)*e.config.QualityDecay)
                if qualityValue < 2 {
                    qualityValue = 2
                }
                if qualityValue > len(qualityChars)-1 {
                    qualityValue = len(qualityChars) - 1
                }
                qualities[i] = qualityChars[qualityValue]
            }
        }
        
        if _, err := fmt.Fprintf(bw, "%s\n%s\n+\n%s\n", readID, sequence, qualities); err != nil {
            return fmt.Errorf("writing FASTQ record: %w", err)
        }
    }

    return bw.Flush()
}

// generateSequence returns a DNA sequence and optional per-base Phred quality scores.
// It tries the HMM model first; if unavailable or erroring, falls back to crypto/rand.
// When falling back, quals is nil (signals the caller to use config-based quality decay).
func (e *Engine) generateSequence(length int) ([]byte, []byte, error) {
    if e.hmmModel != nil {
        result, err := e.hmmModel.GenerateSequence(length)
        if err == nil {
            return result.Sequence, result.Qualities, nil
        }
        // fall through to crypto/rand on per-call HMM failure
    }
    seq, err := e.generateSecureSequence(length)
    return seq, nil, err
}

// generateSecureSequence creates a cryptographically secure random DNA sequence
func (e *Engine) generateSecureSequence(length int) ([]byte, error) {
    bases := []byte("ATCG")
    sequence := make([]byte, length)
    
    // Generate cryptographically secure random bytes
    randomBytes := make([]byte, length)
    _, err := rand.Read(randomBytes)
    if err != nil {
        return nil, fmt.Errorf("failed to generate random bytes: %w", err)
    }
    
    // Convert random bytes to DNA bases
    for i := 0; i < length; i++ {
        sequence[i] = bases[randomBytes[i]%4]
    }
    
    return sequence, nil
}

// VariantRecord holds a single variant for controlled injection into synthetic reads.
// Ref and Alt must be non-empty ACGT strings. Freq is the injection probability in [0,1].
type VariantRecord struct {
    Ref        string
    Alt        string
    Freq       float64
    Chromosome string // "" → position-blind (legacy)
    Position   int64  // 0 → position-blind; 1-based VCF when > 0
}

// GenomicRegion is a half-open [Start, End) interval in 0-based coordinates.
// When set in GenerateOpts, reads are positioned within this region and
// coordinate-aware injection is used for VariantRecords with Chromosome+Position set.
type GenomicRegion struct {
    Chromosome string
    Start      int64 // 0-based inclusive
    End        int64 // 0-based exclusive; must be > Start + ReadLength
}

// GenerateOpts holds optional parameters for GenerateToWriter.
type GenerateOpts struct {
    // Variants are injected into each read with probability Freq.
    Variants []VariantRecord
    // SampleID, when non-empty, appends "\tRG:Z:<id>" to every read header.
    // BWA passes FASTQ comments to SAM when invoked with -C.
    SampleID string
    // RefRegion, when non-nil, enables coordinate-aware injection. Reads are
    // assigned random start positions within the region, and variants with
    // Chromosome+Position set are injected at their exact genomic offsets.
    RefRegion *GenomicRegion // nil → legacy position-blind injection
}

// GenerateToWriter streams all synthetic FASTQ reads to w without writing local
// files. Variant injection and per-read RG tagging are controlled via opts.
// Config-level fingerprinting is applied when enabled.
func (e *Engine) GenerateToWriter(w io.Writer, opts GenerateOpts) error {
    if e.config.ReadLength <= 0 {
        return fmt.Errorf("invalid read length: %d", e.config.ReadLength)
    }
    if e.config.ReadsPerBatch <= 0 {
        return fmt.Errorf("invalid reads per batch: %d", e.config.ReadsPerBatch)
    }
    if e.config.TotalBatches <= 0 {
        return fmt.Errorf("invalid total batches: %d", e.config.TotalBatches)
    }
    if e.qualityFingerprintingEnabled() {
        return fmt.Errorf("steganography layer %q needs file output: a fingerprint cannot be embedded in streamed output; write files instead, or remove %q from steganography_layers", qualityLayer, qualityLayer)
    }

    qualityChars := "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
    bw, ok := w.(*bufio.Writer)
    if !ok {
        bw = bufio.NewWriterSize(w, 256*1024)
    }

    for batch := 1; batch <= e.config.TotalBatches; batch++ {
        for readNum := 1; readNum <= e.config.ReadsPerBatch; readNum++ {
            // Compute read start position for coordinate-aware injection.
            var readStart int64
            if opts.RefRegion != nil {
                regionLen := opts.RefRegion.End - opts.RefRegion.Start - int64(e.config.ReadLength)
                if regionLen <= 0 {
                    return fmt.Errorf("region [%s:%d-%d) too small for read length %d",
                        opts.RefRegion.Chromosome, opts.RefRegion.Start, opts.RefRegion.End, e.config.ReadLength)
                }
                readStart = opts.RefRegion.Start + cryptoRandInt63n(regionLen)
            }
            sides := 1
            if e.config.PairedEnd {
                sides = 2
            }
            for side := 1; side <= sides; side++ {
                var readID string
                if opts.RefRegion != nil {
                    readID = fmt.Sprintf("@%s:%d:HXXXXXX:%d:1101:%d:%d %d:N:0:1",
                        e.config.SampleName, e.config.Lane, batch,
                        int(readStart%100000), int((readStart/100000)%100000), side)
                } else {
                    readID = fmt.Sprintf("@%s:%d:HXXXXXX:%d:1101:%d:%d %d:N:0:1",
                        e.config.SampleName, e.config.Lane, batch, readNum*1000, readNum*2000, side)
                }
                if opts.SampleID != "" {
                    readID += "\tRG:Z:" + opts.SampleID
                }

                // Apply fingerprinting to read ID if enabled.
                if e.config.Fingerprinting.Enable && e.fingerprintingEngine != nil {
                    if e.config.Fingerprinting.CustomerFingerprint {
                        readID = e.fingerprintingEngine.EmbedCustomerID(readID)
                    }
                    if e.config.Fingerprinting.DatasetFingerprint && e.config.DatasetID != "" {
                        readID = e.fingerprintingEngine.EmbedDatasetFingerprint(readID, e.config.DatasetID)
                    }
                }

                sequence, hmmQuals, err := e.generateSequence(e.config.ReadLength)
                if err != nil {
                    return fmt.Errorf("failed to generate sequence: %w", err)
                }

                // Inject controlled variants.
                if len(opts.Variants) > 0 {
                    if opts.RefRegion != nil {
                        sequence = injectAtPosition(sequence, readStart, opts.Variants)
                    } else {
                        sequence = injectVariants(sequence, opts.Variants)
                    }
                }

                // Encode quality scores.
                qualities := make([]byte, e.config.ReadLength)
                if hmmQuals != nil {
                    for i, q := range hmmQuals {
                        idx := int(q)
                        if idx >= len(qualityChars) {
                            idx = len(qualityChars) - 1
                        }
                        qualities[i] = qualityChars[idx]
                    }
                } else {
                    for i := 0; i < e.config.ReadLength; i++ {
                        qv := e.config.BaseQuality - int(float64(i)*e.config.QualityDecay)
                        if qv < 2 {
                            qv = 2
                        }
                        if qv > len(qualityChars)-1 {
                            qv = len(qualityChars) - 1
                        }
                        qualities[i] = qualityChars[qv]
                    }
                }

                if _, err := fmt.Fprintf(bw, "%s\n%s\n+\n%s\n", readID, sequence, qualities); err != nil {
                    return fmt.Errorf("writing FASTQ record: %w", err)
                }
            }
        }
    }
    return bw.Flush()
}

// injectVariants probabilistically substitutes REF alleles with ALT in a read.
// For each variant, a single random byte from crypto/rand determines whether to
// inject (byte/255 < Freq). The read is always returned at the original length:
// longer ALTs are truncated, shorter ALTs are padded with the following bases.
func injectVariants(seq []byte, variants []VariantRecord) []byte {
    result := make([]byte, len(seq))
    copy(result, seq)
    rng := make([]byte, 1)
    for _, v := range variants {
        if len(v.Ref) == 0 || len(v.Alt) == 0 {
            continue
        }
        rand.Read(rng) // #nosec G104 -- crypto/rand never errors on Linux/macOS
        if float64(rng[0])/255.0 >= v.Freq {
            continue
        }
        idx := bytes.Index(result, []byte(v.Ref))
        if idx < 0 {
            continue
        }
        // Splice ALT over REF.
        after := make([]byte, 0, len(seq))
        after = append(after, result[:idx]...)
        after = append(after, []byte(v.Alt)...)
        after = append(after, result[idx+len(v.Ref):]...)
        // Trim or pad to original length.
        if len(after) > len(seq) {
            after = after[:len(seq)]
        }
        for len(after) < len(seq) {
            after = append(after, 'N')
        }
        result = after
    }
    return result
}

// cryptoRandInt63n returns a cryptographically random int64 in [0, n).
// Returns 0 when n <= 0.
func cryptoRandInt63n(n int64) int64 {
    if n <= 0 {
        return 0
    }
    b := make([]byte, 8)
    rand.Read(b) // #nosec G104
    v := uint64(b[0])<<56 | uint64(b[1])<<48 | uint64(b[2])<<40 | uint64(b[3])<<32 |
        uint64(b[4])<<24 | uint64(b[5])<<16 | uint64(b[6])<<8 | uint64(b[7])
    v &= 0x7fffffffffffffff // clear sign bit
    return int64(v % uint64(n))
}

// injectAtPosition injects variants at their exact genomic offsets within a read
// that starts at readStart (0-based). For each variant with Chromosome and Position
// set, the base at offset (Position-1 - readStart) is deterministically set to
// either ALT (with probability Freq) or REF (otherwise), ensuring VAF = Freq.
// Variants without Chromosome/Position are skipped.
func injectAtPosition(seq []byte, readStart int64, variants []VariantRecord) []byte {
    result := make([]byte, len(seq))
    copy(result, seq)
    rng := make([]byte, 1)
    for _, v := range variants {
        if v.Chromosome == "" || v.Position == 0 {
            continue
        }
        if len(v.Ref) == 0 || len(v.Alt) == 0 {
            continue
        }
        // v.Position is 1-based; convert to 0-based offset within this read.
        offset := (v.Position - 1) - readStart
        if offset < 0 || int(offset)+len(v.Ref) > len(result) {
            continue
        }
        rand.Read(rng) // #nosec G104
        var inject []byte
        if float64(rng[0])/255.0 < v.Freq {
            inject = []byte(v.Alt)
        } else {
            inject = []byte(v.Ref)
        }
        ioff := int(offset)
        after := make([]byte, 0, len(result))
        after = append(after, result[:ioff]...)
        after = append(after, inject...)
        after = append(after, result[ioff+len(v.Ref):]...)
        if len(after) > len(result) {
            after = after[:len(result)]
        }
        for len(after) < len(result) {
            after = append(after, 'N')
        }
        result = after
    }
    return result
}

// GenerateRawSequence exposes the internal sequence generator for packages that
// build custom read loops (e.g., haplotype-aware clinical scenario generation).
// Returns DNA bytes and optional per-base quality indices (nil when no HMM model).
func (e *Engine) GenerateRawSequence(readLength int) ([]byte, []byte, error) {
    return e.generateSequence(readLength)
}

func boolToInt(b bool) int {
    if b {
        return 1
    }
    return 0
}