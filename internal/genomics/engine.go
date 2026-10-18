// Package genomics provides the core FASTQ generation engine, including HMM model
// integration, steganographic fingerprinting, and sequencer-realistic output.
//
// Generation runs through one shared per-read pipeline (pipeline.go) with
// extension points for read lengths (engine_readlen.go), qualities and
// sequencing errors (engine_errors.go), and headers, file names and output
// handling (engine_output.go). Randomness is seedable (rng.go).
package genomics

import (
    "bufio"
    "fmt"
    "io"
    "math/rand/v2"
    "os"
    "path/filepath"
    "sync"
    "sync/atomic"

    "github.com/scttfrdmn/phantom-science/internal/fingerprinting"
    "github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
    "github.com/scttfrdmn/phantom-science/internal/security"
    "github.com/scttfrdmn/phantom-science/pkg/types"
)

type Engine struct {
    config               types.SequencerConfig
    fingerprintingEngine fingerprinting.Engine
    hmmModel             authenticity.SequenceGenerator // nil = uniform random bases

    // seed is the run seed (SequencerConfig.Seed, or drawn when unset);
    // seedDrawn reports a drawn seed, which is logged once to logw.
    seed        uint64
    seedDrawn   bool
    seedErr     error
    seedLogOnce sync.Once
    // logw receives diagnostics (seed, HMM model loading). It is stderr so
    // that FASTQ streamed to stdout stays clean.
    logw io.Writer

    // rawMu guards rawRand, the stream behind GenerateRawSequence.
    rawMu   sync.Mutex
    rawRand *rand.Rand

    // overlapSkips counts positioned variants skipped because they overlap
    // a variant already applied to the same read.
    overlapSkips atomic.Int64

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

// NewEngine creates an engine for config. It resolves the run seed (drawing
// one when config.Seed is zero) and loads the HMM model when configured.
// Diagnostics go to stderr.
func NewEngine(config types.SequencerConfig) *Engine {
    // Initialize fingerprinting engine with customer ID from config
    customerID := "phantom-science" // Default customer ID
    if config.Enterprise.CustomerID != "" {
        customerID = config.Enterprise.CustomerID
    }

    e := &Engine{config: config, logw: os.Stderr}
    e.seed, e.seedDrawn, e.seedErr = ResolveSeed(config.Seed)
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
            fmt.Fprintf(e.logw, "   [warn] HMM model load failed (%s): %v — using uniform random bases\n",
                modelPath, err)
        } else {
            e.hmmModel = m
            e.reseedModel(newStream(e.seed, rawStream))
            fmt.Fprintf(e.logw, "   [info] Loaded HMM model from %s\n", modelPath)
        }
    }
    return e
}

// resolveModelPath returns the effective HMM model file path for config.
// Resolution order:
//  1. Explicit HMMModelPath — highest priority.
//  2. ModelDir/<Population>.hmm.json — when both Population and ModelDir are set.
//  3. Empty string — engine uses uniform random bases.
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
    if err := e.checkExtensions(); err != nil {
        return err
    }
    if e.qualityFingerprintingEnabled() {
        if err := e.checkQualityFingerprinting(); err != nil {
            return err
        }
    }
    e.logSeed()

    fmt.Println("🧬 Initializing genomics simulation engine...")
    fmt.Printf("   Sample: %s\n", e.config.SampleName)
    fmt.Printf("   Reads per batch: %d\n", e.config.ReadsPerBatch)
    fmt.Printf("   Total batches: %d\n", e.config.TotalBatches)

    // Create output directory with security validation
    if err := security.SecureMkdirAll(e.config.OutputDir, ""); err != nil {
        return fmt.Errorf("failed to create output directory: %w", err)
    }

    totalReads := e.config.TotalBatches * e.config.ReadsPerBatch
    fmt.Printf("   Total reads to generate: %d\n", totalReads)

    obs, err := e.observer()
    if err != nil {
        return err
    }
    mates := []int{1}
    if e.config.PairedEnd {
        mates = []int{1, 2}
    }
    for batch := 1; batch <= e.config.TotalBatches; batch++ {
        fmt.Printf("   Generating batch %d/%d...\n", batch, e.config.TotalBatches)
        for _, mate := range mates {
            if err := e.generateFASTQFile(e.outputPath(batch, mate), batch, mate, obs); err != nil {
                closeObserver(obs) // #nosec G104 -- the generation error is returned
                if e.config.PairedEnd {
                    return fmt.Errorf("failed to generate R%d file: %w", mate, err)
                }
                return fmt.Errorf("failed to generate FASTQ file: %w", err)
            }
        }
    }
    if err := closeObserver(obs); err != nil {
        return err
    }

    fmt.Printf("   ✅ Generated %d FASTQ files in %s\n",
        e.config.TotalBatches*len(mates),
        e.config.OutputDir)

    return nil
}

// closeObserver closes obs when it is not nil.
func closeObserver(obs ReadObserver) error {
    if obs == nil {
        return nil
    }
    if err := obs.Close(); err != nil {
        return fmt.Errorf("closing read observer: %w", err)
    }
    return nil
}

// Close cleans up engine resources
func (e *Engine) Close() {
    if e.fingerprintingEngine != nil {
        e.fingerprintingEngine.Close()
    }
}

// GenerateFASTQFile generates a single FASTQ file (public method for TUI).
// read is the mate number (1, or 2 for the second file of a pair). The
// file's reads come from the RNG stream of (batch, read), so its content
// depends only on the seed, the configuration, batch and read.
func (e *Engine) GenerateFASTQFile(filename string, batch, read int) error {
    if err := e.checkExtensions(); err != nil {
        return err
    }
    obs, err := e.observer()
    if err != nil {
        return err
    }
    if err := e.generateFASTQFile(filename, batch, read, obs); err != nil {
        closeObserver(obs) // #nosec G104 -- the generation error is returned
        return err
    }
    return closeObserver(obs)
}

func (e *Engine) generateFASTQFile(filename string, batch, read int, obs ReadObserver) error {
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
    if err := e.writeFileStream(file, batch, read, filename, obs); err != nil {
        file.Close()        // #nosec G104 -- the write error is returned
        os.Remove(filename) // #nosec G104 -- best effort cleanup of a partial file
        return err
    }
    if err := file.Close(); err != nil {
        return fmt.Errorf("closing %s: %w", filepath.Base(filename), err)
    }
    final, err := e.finalize(filename)
    if err != nil {
        return err
    }
    e.filesMu.Lock()
    e.files = append(e.files, final)
    e.filesMu.Unlock()
    return nil
}

// GeneratedFiles returns the FASTQ files this engine has written, in order.
func (e *Engine) GeneratedFiles() []string {
    e.filesMu.Lock()
    defer e.filesMu.Unlock()
    return append([]string(nil), e.files...)
}

// SkippedVariantOverlaps returns how many positioned variant applications
// were skipped because the variant overlapped another variant already
// applied to the same read.
func (e *Engine) SkippedVariantOverlaps() int64 {
    return e.overlapSkips.Load()
}

// writeFileStream writes one file's batch of FASTQ records (one mate) to w.
func (e *Engine) writeFileStream(w io.Writer, batch, mate int, filename string, obs ReadObserver) error {
    p, err := e.newPipeline(GenerateOpts{}, obs)
    if err != nil {
        return err
    }
    ww, err := e.wrapOutput(w)
    if err != nil {
        return err
    }
    r := newStream(e.seed, streamID(batch, mate))
    e.reseedModel(r)
    bw := bufio.NewWriterSize(ww, 256*1024)
    var serial int64
    if err := p.writeBatch(bw, r, batch, []int{mate}, filename, &serial); err != nil {
        return err
    }
    if err := bw.Flush(); err != nil {
        return fmt.Errorf("writing FASTQ record: %w", err)
    }
    return ww.Close()
}

// reseedModel reseeds the HMM model (when it supports seeding) from r, so
// HMM-generated reads are reproducible from the run seed as well. The model
// is shared by the engine's streams; generating streams concurrently needs a
// model per stream.
func (e *Engine) reseedModel(r *rand.Rand) {
    if s, ok := e.hmmModel.(seedable); ok {
        s.SetSeed(r.Int64())
    }
}

// generateWith returns n bases and their Phred qualities. It uses the HMM
// model when one is loaded, and falls back to uniform random bases drawn from
// r when there is none or it fails; qualities are then nil, and the pipeline
// applies the quality model.
func (e *Engine) generateWith(r *rand.Rand, n int) ([]byte, []byte, error) {
    if e.hmmModel != nil {
        result, err := e.hmmModel.GenerateSequence(n)
        if err == nil {
            return result.Sequence, result.Qualities, nil
        }
        // fall through to uniform bases on per-call HMM failure
    }
    return randomBases(r, n), nil, nil
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
    // Chromosome+Position set (on the region's chromosome) are injected at
    // their exact genomic offsets; position-blind variants are ignored.
    RefRegion *GenomicRegion // nil → legacy position-blind injection
}

// GenerateToWriter streams all synthetic FASTQ reads to w without writing local
// files. Variant injection and per-read RG tagging are controlled via opts.
// Config-level fingerprinting is applied when enabled. Each batch draws from
// its own RNG stream, so output is reproducible from the seed.
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
    if err := e.checkExtensions(); err != nil {
        return err
    }
    e.logSeed()

    obs, err := e.observer()
    if err != nil {
        return err
    }
    before := e.SkippedVariantOverlaps()
    if err := e.writeStreams(w, opts, obs); err != nil {
        closeObserver(obs) // #nosec G104 -- the generation error is returned
        return err
    }
    if n := e.SkippedVariantOverlaps() - before; n > 0 {
        fmt.Fprintf(e.logw, "   [warn] %d variant applications skipped: the variant overlapped another variant in the same read\n", n)
    }
    return closeObserver(obs)
}

// writeStreams writes every batch to w, one RNG stream per batch.
func (e *Engine) writeStreams(w io.Writer, opts GenerateOpts, obs ReadObserver) error {
    ww, err := e.wrapOutput(w)
    if err != nil {
        return err
    }
    // Reuse the caller's buffer when the output is passed through unchanged.
    var bw *bufio.Writer
    if nw, isNop := ww.(nopWriteCloser); isNop {
        bw, _ = nw.Writer.(*bufio.Writer)
    }
    if bw == nil {
        bw = bufio.NewWriterSize(ww, 256*1024)
    }
    mates := []int{1}
    if e.config.PairedEnd {
        mates = []int{1, 2}
    }
    var serial int64
    for batch := 1; batch <= e.config.TotalBatches; batch++ {
        p, err := e.newPipeline(opts, obs)
        if err != nil {
            return err
        }
        r := newStream(e.seed, streamID(batch, 0))
        e.reseedModel(r)
        if err := p.writeBatch(bw, r, batch, mates, "", &serial); err != nil {
            return err
        }
    }
    if err := bw.Flush(); err != nil {
        return fmt.Errorf("writing FASTQ record: %w", err)
    }
    return ww.Close()
}

// GenerateRawSequence exposes the internal sequence generator for packages that
// build custom read loops (e.g., haplotype-aware clinical scenario generation).
// Returns DNA bytes and optional per-base quality indices (nil when no HMM model).
// It draws from an engine-wide stream; GenerateRawSequenceWith takes the
// generator explicitly.
func (e *Engine) GenerateRawSequence(readLength int) ([]byte, []byte, error) {
    e.rawMu.Lock()
    defer e.rawMu.Unlock()
    if e.rawRand == nil {
        e.rawRand = newStream(e.seed, rawStream|1)
    }
    return e.generateWith(e.rawRand, readLength)
}

// GenerateRawSequenceWith is GenerateRawSequence drawing uniform bases from r
// (for example a generator from Engine.Rand).
func (e *Engine) GenerateRawSequenceWith(r *rand.Rand, readLength int) ([]byte, []byte, error) {
    return e.generateWith(r, readLength)
}
