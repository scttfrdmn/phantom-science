// Package manifest describes one generation run: the tool version, seed and
// configuration that reproduce it, and every output file with its size and
// SHA-256 (#16). The engine writes it as phantom-manifest.json next to the
// output (in the output directory, or under the S3 prefix).
//
// A manifest is a pure function of the run's seed and configuration: it
// holds no wall-clock time, host name or absolute output path, and leaves
// out settings that do not change the output (output_dir, workers), so two
// runs with the same seed produce byte-identical manifests. It never holds
// key material: the fingerprint key file path is removed and only whether
// fingerprints were encrypted is recorded.
package manifest

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"time"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// FileName is the manifest's file name in the output directory or prefix.
const FileName = "phantom-manifest.json"

// Format and Version identify the manifest schema.
const (
	Format  = "phantom-manifest"
	Version = 1
)

// Output roles.
const (
	RoleFASTQ   = "fastq"   // generated reads
	RoleSummary = "summary" // run summary (ONT sequencing_summary.txt)
	RoleTruth   = "truth"   // ground truth (variants, read origins, source sequences)
	RoleSignal  = "signal"  // raw nanopore signal (BLOW5)
)

// Manifest is the phantom-manifest.json document.
type Manifest struct {
	Format  string `json:"format"`
	Version int    `json:"version"`
	Tool    Tool   `json:"tool"`
	// Seed is the run seed: passing it (with Config) reproduces every
	// output byte for byte (except quality-score fingerprints, which embed
	// the time they were made).
	Seed uint64 `json:"seed"`
	// RunStart is the run's acquisition start, derived from the seed (it is
	// the ONT start_time origin); it is not the wall-clock time of the run.
	// A signal run takes it from its input reads' start times.
	RunStart time.Time `json:"run_start"`
	// Destination is the S3 URI the outputs were uploaded under; empty for
	// local output (paths are relative to the manifest's directory).
	Destination string `json:"destination,omitempty"`
	// Config is the sanitized configuration (see SanitizeConfig). A signal
	// run has none and records its parameters in Signal instead.
	Config types.SequencerConfig `json:"config,omitzero"`
	// Signal holds the parameters of a raw-signal run (the signal
	// command); nil for read generation.
	Signal *Signal `json:"signal,omitempty"`
	// Encrypted reports whether quality-score fingerprints were encrypted
	// with a key (the key itself is never recorded).
	Encrypted bool       `json:"fingerprint_encrypted,omitempty"`
	FlowCells []FlowCell `json:"flow_cells"`
	Outputs   []Output   `json:"outputs"`
	Totals    Totals     `json:"totals"`
}

// Signal describes a raw-signal run: with the manifest Seed and the inputs
// it names, it reproduces every BLOW5 byte. Inputs are recorded by base name,
// size and SHA-256, never by path.
type Signal struct {
	Profile string `json:"profile"`
	// ProfileOverrides are the parameters set with signal --profile-set.
	ProfileOverrides  map[string]float64 `json:"profile_overrides,omitempty"`
	PoreModel         Input              `json:"pore_model"`
	ModelUnits        string             `json:"model_units"` // "pa" or "standardised"
	RecordCompression string             `json:"record_compression"`
	SignalCompression string             `json:"signal_compression"`
	ReadsPerFile      int                `json:"reads_per_file"`
	Inputs            []Input            `json:"inputs"`
	Truth             *Input             `json:"truth,omitempty"`
	// Samples is the total number of signal samples written.
	Samples int64 `json:"samples"`
	// Skipped counts input reads shorter than the pore model's k, which
	// have no signal.
	Skipped int64 `json:"skipped,omitempty"`
}

// Input is an input file of a run.
type Input struct {
	Name   string `json:"name"` // base name
	Bytes  int64  `json:"bytes"`
	SHA256 string `json:"sha256"`
}

// Tool identifies the generator.
type Tool struct {
	Name    string `json:"name"`
	Version string `json:"version"`
}

// FlowCell describes one independent run (flow cell) of the manifest.
type FlowCell struct {
	Index      int    `json:"index"` // 1-based
	Seed       uint64 `json:"seed"`
	RunID      string `json:"run_id,omitempty"`
	FlowCellID string `json:"flow_cell_id,omitempty"`
	// Dir is the flow cell's subdirectory ("" when the run has one flow
	// cell and writes to the output root).
	Dir   string `json:"dir,omitempty"`
	Files int    `json:"files"`
	Reads int64  `json:"reads"`
	Bases int64  `json:"bases"`
	// Capped reports that the max_output cap stopped this flow cell early.
	Capped bool `json:"capped,omitempty"`
}

// Output is one output file.
type Output struct {
	// Path is relative to the output root, with forward slashes.
	Path string `json:"path"`
	// Key is the S3 object key for uploaded outputs.
	Key      string `json:"key,omitempty"`
	Role     string `json:"role"`
	FlowCell int    `json:"flow_cell"`
	Bytes    int64  `json:"bytes"`
	SHA256   string `json:"sha256"`
	// Reads and Bases count the records and bases of a fastq output, or the
	// reads and source bases of a signal output.
	Reads int64 `json:"reads,omitempty"`
	Bases int64 `json:"bases,omitempty"`
}

// Totals sums the read outputs: FASTQ files, or BLOW5 files of a signal run
// (Bytes sums every output).
type Totals struct {
	Files int   `json:"files"`
	Reads int64 `json:"reads"`
	Bases int64 `json:"bases"`
	Bytes int64 `json:"bytes"`
}

// New returns an empty manifest for a run.
func New(toolVersion string, seed uint64, runStart time.Time, config types.SequencerConfig) *Manifest {
	if toolVersion == "" {
		toolVersion = "dev"
	}
	return &Manifest{
		Format:    Format,
		Version:   Version,
		Tool:      Tool{Name: "phantom-science", Version: toolVersion},
		Seed:      seed,
		RunStart:  runStart.UTC(),
		Config:    SanitizeConfig(config),
		FlowCells: []FlowCell{},
		Outputs:   []Output{},
	}
}

// SanitizeConfig returns c without key material, local output locations or
// settings that do not affect the output: fingerprint key file paths, the
// licence key, OutputDir and Workers are cleared.
func SanitizeConfig(c types.SequencerConfig) types.SequencerConfig {
	c.OutputDir = ""
	c.Workers = 0
	c.Fingerprinting.KeyFile = ""
	c.Enterprise.Fingerprinting.KeyFile = ""
	c.Enterprise.LicenseKey = ""
	c.Fingerprinting.SteganographyLayers = append([]string(nil), c.Fingerprinting.SteganographyLayers...)
	c.Enterprise.Fingerprinting.SteganographyLayers = append([]string(nil), c.Enterprise.Fingerprinting.SteganographyLayers...)
	return c
}

// Add appends an output and updates the totals.
func (m *Manifest) Add(o Output) {
	m.Outputs = append(m.Outputs, o)
	m.Totals.Bytes += o.Bytes
	if o.Role == RoleFASTQ || o.Role == RoleSignal {
		m.Totals.Files++
		m.Totals.Reads += o.Reads
		m.Totals.Bases += o.Bases
	}
}

// Marshal returns the manifest as indented JSON with a trailing newline.
func (m *Manifest) Marshal() ([]byte, error) {
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(m); err != nil {
		return nil, fmt.Errorf("encoding manifest: %w", err)
	}
	return buf.Bytes(), nil
}

// Read parses the manifest file at path.
func Read(path string) (*Manifest, error) {
	data, err := os.ReadFile(path) // #nosec G304 -- caller-chosen manifest path
	if err != nil {
		return nil, err
	}
	var m Manifest
	if err := json.Unmarshal(data, &m); err != nil {
		return nil, fmt.Errorf("parsing %s: %w", filepath.Base(path), err)
	}
	if m.Format != Format {
		return nil, fmt.Errorf("%s: format %q, want %q", filepath.Base(path), m.Format, Format)
	}
	return &m, nil
}

// HashFile returns the size and hex SHA-256 of the file at path.
func HashFile(path string) (int64, string, error) {
	f, err := os.Open(path) // #nosec G304 -- an output file of this run
	if err != nil {
		return 0, "", err
	}
	defer f.Close()
	h := sha256.New()
	n, err := io.Copy(h, f)
	if err != nil {
		return 0, "", fmt.Errorf("hashing %s: %w", filepath.Base(path), err)
	}
	return n, hex.EncodeToString(h.Sum(nil)), nil
}

// VerifyDir checks that every output exists under root (the manifest's
// directory) with the recorded size and SHA-256.
func (m *Manifest) VerifyDir(root string) error {
	var errs []error
	for _, o := range m.Outputs {
		n, sum, err := HashFile(filepath.Join(root, filepath.FromSlash(o.Path)))
		switch {
		case err != nil:
			errs = append(errs, err)
		case n != o.Bytes || sum != o.SHA256:
			errs = append(errs, fmt.Errorf("%s: %d bytes sha256 %s, manifest has %d bytes sha256 %s", o.Path, n, sum, o.Bytes, o.SHA256))
		}
	}
	return errors.Join(errs...)
}
