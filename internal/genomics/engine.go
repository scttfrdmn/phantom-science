// Package genomics provides the core FASTQ generation engine, including HMM model
// integration, keyed-seed fingerprinting, and sequencer-realistic output.
//
// Generation runs through one shared per-read pipeline (pipeline.go) with
// extension points for read lengths (engine_readlen.go), qualities and
// sequencing errors (engine_errors.go), and headers, file names and output
// handling (engine_output.go). Randomness is seedable (rng.go).
package genomics

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"io/fs"
	"math/rand/v2"
	"os"
	"path/filepath"
	"sync"
	"sync/atomic"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/scttfrdmn/phantom-science/internal/security"
	"github.com/scttfrdmn/phantom-science/internal/storage/s3sink"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

type Engine struct {
	config   types.SequencerConfig
	hmmModel authenticity.SequenceGenerator // nil = uniform random bases

	// seed is the run seed (SequencerConfig.Seed, or drawn when unset, or,
	// with a keyed-seed fingerprint, derived from the key; see
	// engine_fingerprint.go); seedDrawn reports a drawn seed, which is
	// logged once to logw.
	seed        uint64
	seedDrawn   bool
	seedErr     error
	seedLogOnce sync.Once
	// logw receives diagnostics (seed, HMM model loading, warnings). It is
	// stderr unless set with NewEngineLog, so that FASTQ streamed to stdout
	// stays clean.
	logw io.Writer

	// customerID is the customer written in read-ID labels.
	customerID string
	// fragments, when set, supplies read and fragment bases (see
	// FragmentSampler); nil uses the HMM model or uniform bases.
	fragments FragmentSampler
	// ref is the reference genome reads are sampled from (ReferencePath),
	// shared read-only by every stream; refErr records a reference that
	// could not be used. See engine_reference.go.
	ref    *reference.Sampler
	refErr error
	// refMaxRead, when > 0, is the longest single-end read the reference
	// can supply; longer sampled lengths are redrawn (see fitLengths).
	refMaxRead int64

	// modelErr records a configured HMM model that could not be loaded;
	// generation fails with it rather than silently using uniform bases.
	modelErr error

	// rawMu guards rawRand, the stream behind GenerateRawSequence.
	rawMu   sync.Mutex
	rawRand *rand.Rand

	// overlapSkips counts positioned variants skipped because they overlap
	// a variant already applied to the same read (in the run's engine; see
	// countOverlapSkips).
	overlapSkips atomic.Int64

	// variants are the variants of variants_path, shared read-only by every
	// stream, and variantsInput the VCF as the manifest records it;
	// varRegion is the region reads are placed in without a reference, and
	// refVars the checked positioned variants with one. varErr records
	// variant settings that cannot be used. See engine_variants.go.
	variants      []VariantRecord
	variantsInput *manifest.Input
	varRegion     *GenomicRegion
	refVars       *refVariants
	varErr        error

	// inputs are the run's other input files as the manifest records them,
	// hashed once (inputsOnce); inputsErr records a file that could not be
	// hashed. See engine_inputs.go.
	inputsOnce sync.Once
	inputs     runInputs
	inputsErr  error

	// Files written by Generate and GenerateFASTQFile, in order.
	filesMu sync.Mutex
	files   []string

	// fp is the run's keyed-seed fingerprint; fpErr records fingerprint
	// settings that cannot be used (no key, quality_scores layer, ...).
	fp    keyedSeed
	fpErr error

	// modelPath is the HMM model file hmmModel was loaded from; file output
	// loads a fresh instance per file from it (models are not safe to share
	// between parallel workers, and adaptive models carry state).
	modelPath string
	// parent is the run's engine when this one generates one flow cell of
	// it.
	parent *Engine
	// toolVersion is recorded in the manifest; uploader receives S3 output.
	toolVersion string
	uploader    s3sink.Uploader
	// manifest describes the last completed Generate run.
	manifest *manifest.Manifest
}

// SetToolVersion sets the tool version recorded in run manifests.
func (e *Engine) SetToolVersion(v string) { e.toolVersion = v }

// SetUploader sets the uploader used for s3:// output (by default the AWS
// SDK's, with the default credential chain).
func (e *Engine) SetUploader(u s3sink.Uploader) { e.uploader = u }

// Manifest returns the manifest of the last completed Generate run, or nil.
func (e *Engine) Manifest() *manifest.Manifest { return e.manifest }

// NewEngine creates an engine for config. It resolves the run seed (drawing
// one when config.Seed is zero) and loads the HMM model when configured.
// Diagnostics go to stderr (see NewEngineLog).
func NewEngine(config types.SequencerConfig) *Engine {
	return NewEngineLog(config, os.Stderr)
}

// NewEngineLog is NewEngine with diagnostics (the drawn seed, HMM model
// loading, skipped variants, the output cap) written to logw instead of
// stderr, for callers such as the TUI that own the terminal. A nil logw
// discards them.
func NewEngineLog(config types.SequencerConfig, logw io.Writer) *Engine {
	if logw == nil {
		logw = io.Discard
	}
	// The customer named by read-ID labels.
	customerID := "phantom-science" // Default customer ID
	if config.Enterprise.CustomerID != "" {
		customerID = config.Enterprise.CustomerID
	}

	e := &Engine{config: config, logw: logw, customerID: customerID}
	e.resolveSeed()
	// The reference first: with both a reference and a model configured,
	// the conflict is the error and neither is loaded.
	e.loadReference()
	if modelPath := resolveModelPath(config); modelPath != "" && config.ReferencePath == "" {
		m, err := authenticity.LoadModel(modelPath)
		switch {
		case err != nil && config.HMMModelPath == "" && errors.Is(err, fs.ErrNotExist):
			e.modelErr = fmt.Errorf("no HMM model for population %s in model_dir %s (expected %s)",
				config.Population, config.ModelDir, modelPath)
		case err != nil:
			e.modelErr = fmt.Errorf("HMM model %s could not be loaded: %w", modelPath, err)
		default:
			e.hmmModel = m
			e.modelPath = modelPath
			e.reseedModel(newStream(e.seed, rawStream))
			fmt.Fprintf(e.logw, "   [info] Loaded HMM model from %s\n", modelPath)
		}
	}
	e.loadVariants()
	return e
}

// resolveModelPath returns the effective HMM model file path for config. A
// configured model that cannot be loaded, or a ModelDir with no
// <Population>.hmm.json, is a generation error; with neither HMMModelPath nor
// ModelDir+Population set, generation uses uniform random bases.
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
	if err := e.checkBatches(); err != nil {
		return err
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
	e.logSeed()

	c := e.config
	fmt.Println("🧬 Initializing genomics simulation engine...")
	fmt.Printf("   Sample: %s\n", c.SampleName)
	fmt.Printf("   Reads per batch: %d\n", c.ReadsPerBatch)
	switch {
	case c.TotalBases > 0:
		fmt.Printf("   Target bases: %d\n", c.TotalBases)
	case c.TotalReads > 0:
		fmt.Printf("   Target reads: %d\n", c.TotalReads)
	default:
		fmt.Printf("   Total batches: %d\n", c.TotalBatches)
		fmt.Printf("   Total reads to generate: %d\n", c.TotalBatches*c.ReadsPerBatch)
	}
	if c.FlowCells > 1 {
		fmt.Printf("   Flow cells: %d\n", c.FlowCells)
	}
	if c.Workers > 1 {
		fmt.Printf("   Workers: %d\n", c.Workers)
	}

	before := len(e.GeneratedFiles())
	skips := e.SkippedVariantOverlaps()
	if err := e.generateRun(); err != nil {
		return err
	}
	e.reportOverlapSkips(skips)
	fmt.Printf("   ✅ Generated %d FASTQ files in %s\n", len(e.GeneratedFiles())-before, c.OutputDir)
	return nil
}

// checkBatches validates the batch count, which output-size targets replace.
func (e *Engine) checkBatches() error {
	c := e.config
	if c.TotalBases == 0 && c.TotalReads == 0 && c.TotalBatches <= 0 {
		return fmt.Errorf("invalid total batches: must be greater than 0, got %d", c.TotalBatches)
	}
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

// Close cleans up engine resources. The engine holds none that need
// releasing; Close is kept for callers that defer it.
func (e *Engine) Close() {}

// GenerateFASTQFile generates a single FASTQ file (public method for TUI).
// read is the mate number (1, or 2 for the second file of a pair). The
// file's reads come from the RNG stream of (batch, read) (and, for
// paired-end output, the batch's pair stream), so its content depends only
// on the seed, the configuration, batch and read: it equals that file of a
// Generate run. Output-size targets and flow cells do not apply; the file
// holds ReadsPerBatch reads.
func (e *Engine) GenerateFASTQFile(filename string, batch, read int) error {
	if err := e.checkExtensions(); err != nil {
		return err
	}
	if s3sink.IsURI(e.config.OutputDir) {
		return fmt.Errorf("GenerateFASTQFile writes local files: output_dir %s is an S3 URI", e.config.OutputDir)
	}
	obs, err := e.observer()
	if err != nil {
		return err
	}
	final, err := e.generateFASTQFile(fileJob{path: filename, batch: batch, mate: read, reads: e.config.ReadsPerBatch, obs: obs})
	if err != nil {
		closeObserver(obs) // #nosec G104 -- the generation error is returned
		return err
	}
	e.filesMu.Lock()
	e.files = append(e.files, final)
	e.filesMu.Unlock()
	return closeObserver(obs)
}

// fileJob describes one FASTQ file to write.
type fileJob struct {
	path               string
	batch, mate, reads int
	obs                ReadObserver
	tally              func(bases int) // told each written read's length
	abort              *atomic.Bool    // stops the file with errAborted when set
	// truth, when set, is the truth table written alongside (reference
	// mode, mate 1); variants, when set, the variants table (variants_path,
	// mate 1).
	truth, variants string
	// plain writes the records uncompressed whatever the output
	// compression (regeneration, which discards them).
	plain bool
}

// generateFASTQFile writes and finalizes one file and returns its final path.
func (e *Engine) generateFASTQFile(j fileJob) (string, error) {
	// Validate file path against the configured output directory
	if err := security.ValidateFilePath(j.path, e.config.OutputDir); err != nil {
		return "", fmt.Errorf("invalid file path: %w", err)
	}

	file, err := security.SecureFileCreate(j.path, e.config.OutputDir)
	if err != nil {
		return "", err
	}
	// side are the tables written alongside: the truth and variants
	// tables, when set.
	side := []struct {
		path string
		f    *os.File
	}{{path: j.truth}, {path: j.variants}}
	fail := func(err error) (string, error) {
		file.Close()      // #nosec G104 -- the write error is returned
		os.Remove(j.path) // #nosec G104 -- best effort cleanup of a partial file
		for _, t := range side {
			if t.f != nil {
				t.f.Close() // #nosec G104
			}
			if t.path != "" {
				os.Remove(t.path) // #nosec G104
			}
		}
		return "", err
	}
	for i := range side {
		if side[i].path == "" {
			continue
		}
		if side[i].f, err = security.SecureFileCreate(side[i].path, e.config.OutputDir); err != nil {
			return fail(err)
		}
	}
	if err := e.writeFileStream(file, side[0].f, side[1].f, j); err != nil {
		return fail(err)
	}
	for i := range side {
		if side[i].f == nil {
			continue
		}
		err := side[i].f.Close()
		side[i].f = nil
		if err != nil {
			return fail(fmt.Errorf("closing %s: %w", filepath.Base(side[i].path), err))
		}
	}
	if err := file.Close(); err != nil {
		return "", fmt.Errorf("closing %s: %w", filepath.Base(j.path), err)
	}
	return e.finalize(j.path)
}

// GeneratedFiles returns the FASTQ files this engine has written, in order
// (s3:// URIs for S3 output).
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

// reportOverlapSkips warns about the overlap skips counted since before.
func (e *Engine) reportOverlapSkips(before int64) {
	if n := e.SkippedVariantOverlaps() - before; n > 0 {
		fmt.Fprintf(e.logw, "   [warn] %d variant applications skipped: the variant overlapped another variant in the same read\n", n)
	}
}

// writeFileStream writes one file's batch of FASTQ records (one mate) to w,
// and, when truth (variants) is not nil, the batch's truth (variants) table
// to it. The file gets its own HMM model instance, seeded from its stream,
// and the engine's configured variants.
func (e *Engine) writeFileStream(w io.Writer, truth, variants *os.File, j fileJob) error {
	p, err := e.newPipeline(e.configuredOpts(), j.obs)
	if err != nil {
		return err
	}
	if e.ref != nil {
		p.refVars = e.refVars
	}
	p.tally, p.abort = j.tally, j.abort
	if p.model, err = e.streamModel(); err != nil {
		return err
	}
	var ww io.WriteCloser = nopWriteCloser{w}
	if !j.plain {
		if ww, err = e.wrapOutput(w); err != nil {
			return err
		}
	}
	r := newStream(e.seed, streamID(j.batch, j.mate))
	shared := r
	if e.config.PairedEnd {
		// Both mate files regenerate the batch's fragments from the pair
		// stream, so the model that generates them is seeded from it too.
		shared = newStream(e.seed, streamID(j.batch, pairMate))
	}
	reseed(p.model, shared)
	bw := bufio.NewWriterSize(ww, 256*1024)
	if truth != nil {
		p.truth = bufio.NewWriterSize(truth, 64*1024)
		p.truth.WriteString(truthHeader) // #nosec G104 -- bufio errors are sticky and surface at Flush
	}
	if variants != nil {
		p.planted = &plantLog{}
		p.variants = bufio.NewWriterSize(variants, 64*1024)
		p.variants.WriteString(variantsHeader) // #nosec G104 -- bufio errors are sticky and surface at Flush
	}
	var serial int64
	spec := batchSpec{batch: j.batch, reads: j.reads, mates: []int{j.mate}, file: j.path, shared: shared}
	if _, err := p.writeBatch(bw, r, spec, &serial); err != nil {
		return err
	}
	if err := bw.Flush(); err != nil {
		return fmt.Errorf("writing FASTQ record: %w", err)
	}
	if p.truth != nil {
		if err := p.truth.Flush(); err != nil {
			return fmt.Errorf("writing truth table: %w", err)
		}
	}
	if p.variants != nil {
		if err := p.variants.Flush(); err != nil {
			return fmt.Errorf("writing variants table: %w", err)
		}
	}
	return ww.Close()
}

// streamModel returns the HMM model for one output file: a fresh instance
// loaded from the model path, so files do not depend on each other or on
// scheduling. Without a path it is the engine's shared model (or nil).
func (e *Engine) streamModel() (authenticity.SequenceGenerator, error) {
	if e.hmmModel == nil || e.modelPath == "" {
		return e.hmmModel, nil
	}
	m, err := authenticity.LoadModel(e.modelPath)
	if err != nil {
		return nil, fmt.Errorf("loading HMM model %s: %w", e.modelPath, err)
	}
	return m, nil
}

// reseedModel reseeds the engine's shared HMM model from r. Streamed output
// (GenerateToWriter) uses the shared model, one stream at a time.
func (e *Engine) reseedModel(r *rand.Rand) {
	reseed(e.hmmModel, r)
}

// generateWith returns n bases and their Phred qualities. It uses the HMM
// model when one is loaded, and falls back to uniform random bases drawn from
// r when there is none or it fails; qualities are then nil, and the pipeline
// applies the quality model.
func (e *Engine) generateWith(r *rand.Rand, n int) ([]byte, []byte, error) {
	return generateFrom(e.hmmModel, r, n)
}

// VariantRecord holds a single variant for controlled injection into synthetic reads.
// Ref and Alt must be non-empty ACGT strings. Freq is the injection probability in [0,1].
type VariantRecord struct {
	ID         string // VCF ID; "" for none
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
	// Variants are injected into each read with probability Freq. With
	// Variants and RefRegion both unset, the engine's configured variants
	// (variants_path, variants_region) are used; setting them as well as
	// variants_path is an error.
	Variants []VariantRecord
	// SampleID, when non-empty, appends "\tRG:Z:<id>" to every read header
	// comment (before any read-ID labels).
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
// A keyed-seed fingerprint applies as for file output: it is the seed. Each
// batch draws from its own RNG stream, so output is reproducible from the seed. TotalReads,
// TotalBases and the max-output cap apply as for Generate (one stream, so
// Workers does not apply); flow cells need file output.
func (e *Engine) GenerateToWriter(w io.Writer, opts GenerateOpts) error {
	if e.config.ReadLength <= 0 {
		return fmt.Errorf("invalid read length: %d", e.config.ReadLength)
	}
	if e.config.ReadsPerBatch <= 0 {
		return fmt.Errorf("invalid reads per batch: %d", e.config.ReadsPerBatch)
	}
	if err := e.checkBatches(); err != nil {
		return err
	}
	if e.config.FlowCells > 1 {
		return fmt.Errorf("flow_cells %d needs file output: each flow cell is written to its own directory", e.config.FlowCells)
	}
	if err := e.checkExtensions(); err != nil {
		return err
	}
	e.logSeed()

	obs, err := e.observer()
	if err != nil {
		return err
	}
	configured := opts.Variants == nil && opts.RefRegion == nil
	if !configured && e.config.VariantsPath != "" {
		closeObserver(obs) // #nosec G104 -- the option error is returned
		return fmt.Errorf("variants_path %s is configured: do not also pass GenerateOpts.Variants or RefRegion", e.config.VariantsPath)
	}
	refVars := e.refVars
	if configured {
		opts.Variants, opts.RefRegion = e.variants, e.varRegion
	} else if refVars, err = e.referenceVariants(opts.Variants); err != nil {
		closeObserver(obs) // #nosec G104 -- the variant error is returned
		return err
	}
	before := e.SkippedVariantOverlaps()
	if err := e.writeStreams(w, opts, obs, refVars); err != nil {
		closeObserver(obs) // #nosec G104 -- the generation error is returned
		return err
	}
	e.reportOverlapSkips(before)
	if !configured {
		refVars.report(e) // configured variants are reported when loaded
	}
	return closeObserver(obs)
}

// writeStreams writes every batch to w, one RNG stream per batch. refVars
// are the positioned variants of reference mode.
func (e *Engine) writeStreams(w io.Writer, opts GenerateOpts, obs ReadObserver, refVars *refVariants) error {
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
	mates := e.mates()
	plan := e.outputPlan()
	b := &budget{p: plan}
	var serial int64
	for batch := 1; plan.has(batch); batch++ {
		p, err := e.newPipeline(opts, obs)
		if err != nil {
			return err
		}
		p.model, p.admit = e.hmmModel, b.admit
		if e.ref != nil {
			p.refVars = refVars
		}
		r := newStream(e.seed, streamID(batch, 0))
		e.reseedModel(r)
		spec := batchSpec{batch: batch, reads: plan.reads(batch), mates: mates, shared: r}
		stopped, err := p.writeBatch(bw, r, spec, &serial)
		if err != nil {
			return err
		}
		if stopped || (plan.bases > 0 && b.bases >= plan.bases) {
			break
		}
	}
	if b.capped {
		fmt.Fprintf(e.logw, "   [warn] output stopped at the max output of %d bases (max_output_gb / max_output_tb)\n", plan.cap)
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
	if e.modelErr != nil {
		return nil, nil, e.modelErr
	}
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
	if e.modelErr != nil {
		return nil, nil, e.modelErr
	}
	return e.generateWith(r, readLength)
}
