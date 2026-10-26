//go:build !release

package clinical

// Keyed-seed fingerprints and the run manifest of a scenario (#31).
//
// A scenario draws all of its randomness from one run seed: the haplotype
// simulation (seed ^ trioSimSeedSalt), the trio transmission (seed), the
// synthetic window (engine stream backgroundStream of the seed) and every
// member's reads (memberSeed of the seed). With fingerprinting enabled in
// the sequencer config (as for generate: fingerprinting.enable, the key,
// enterprise.customer_id and dataset_id) the run seed is
//
//	fingerprinting.DeriveSeed(key, customer, dataset, salt)
//
// with the configured seed (or a drawn, logged value) as the salt, so a
// fingerprinted scenario is byte-identical to the unfingerprinted scenario
// whose seed is the derived seed, and its statistics are those of any
// scenario run.
//
// Each scenario writes one run manifest, phantom-manifest.json, in the
// directory of its ground_truth file: the sanitized sequencer configuration
// (no paths, no key), the seed or, when fingerprinted, the fingerprint
// record (scheme, customer, dataset, salt) and no seed, the scenario
// specification without output paths, the input files by name, size and
// SHA-256 (reference FASTA, HMM model, read-length histogram), and every
// member FASTQ and the ground truth with size and SHA-256. "phantom-science
// verify --manifest" regenerates the scenario from it (Regenerate).

import (
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	runmanifest "github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Members are the scenario samples, in the order their reads are made.
var Members = []string{"proband", "mother", "father"}

// resolveScenarioSeed returns the scenario's run seed from cfg (see the
// file comment).
func resolveScenarioSeed(cfg types.SequencerConfig) (genomics.RunSeed, error) {
	cfg.ResolveFingerprinting(nil)
	return genomics.ResolveRunSeed(cfg, "scenario")
}

// runRecord describes a scenario run for its run manifest.
type runRecord struct {
	seed      genomics.RunSeed
	reference *runmanifest.Input
	files     []memberFile
	m         *runmanifest.Manifest // set by finish
}

// memberFile is one FASTQ file a scenario wrote.
type memberFile struct {
	member, path string
	reads, bases int64
}

func (r *runRecord) addFile(member, path string, reads, bases int64) {
	r.files = append(r.files, memberFile{member, path, reads, bases})
}

// specRecord is the scenario specification recorded in the run manifest:
// everything but the output paths.
type specRecord struct {
	Phenotype PhenotypeSpec `json:"phenotype"`
	Trio      *TrioSpec     `json:"trio,omitempty"`
	Sequencer SequencerSpec `json:"sequencer"`
}

// finish builds the run's manifest, except for the output files, from the
// generator's configuration and the scenario. cfg is the generation
// configuration (the engines' inputs are recorded from it).
func (r *runRecord) finish(g *ScenarioGenerator, scenario *ClinicalScenario, cfg types.SequencerConfig) error {
	spec, err := json.Marshal(specRecord{scenario.Phenotype, scenario.Trio, scenario.Sequencer})
	if err != nil {
		return fmt.Errorf("recording the scenario: %w", err)
	}
	var members []string
	seen := map[string]bool{}
	for _, f := range r.files {
		if !seen[f.member] {
			seen[f.member] = true
			members = append(members, f.member)
		}
	}
	recorded := r.seed.Seed
	if r.seed.Fingerprint != nil {
		recorded = 0
	}
	m := runmanifest.New("", recorded, time.Time{}, *g.config)
	m.Config.Seed = recorded
	if fp := r.seed.Fingerprint; fp != nil {
		m.Fingerprint = fp
		m.Config.Seed = fp.Salt
	}
	m.Scenario = &runmanifest.Scenario{Spec: spec, Members: members}
	// The engines' other inputs (HMM model, read-length histogram, VCF), as
	// an engine configured like the run's reads them.
	ic := cfg
	ic.ReferencePath, ic.ReferenceRegion, ic.StrictRef = "", "", false
	e := genomics.NewEngineLog(ic, io.Discard)
	defer e.Close()
	if err := e.RecordInputs(m); err != nil {
		return err
	}
	m.Reference = r.reference
	r.m = m
	return nil
}

// WriteRunManifest writes the run manifest of the last Run to path (see the
// file comment); groundTruth, when not empty, is the ground-truth file Run's
// manifest was written to, recorded as an output. Output paths are relative
// to the manifest's directory.
func (g *ScenarioGenerator) WriteRunManifest(path, groundTruth, toolVersion string) error {
	r := g.rec
	if r == nil || r.m == nil {
		return fmt.Errorf("no scenario run to record")
	}
	m := *r.m
	m.Outputs = nil
	m.Totals = runmanifest.Totals{}
	if toolVersion != "" {
		m.Tool.Version = toolVersion
	}
	dir, err := filepath.Abs(filepath.Dir(path))
	if err != nil {
		return err
	}
	rel := func(p string) string {
		abs, err := filepath.Abs(p)
		if err != nil {
			return filepath.Base(p)
		}
		if rp, err := filepath.Rel(dir, abs); err == nil {
			return filepath.ToSlash(rp)
		}
		return filepath.Base(p)
	}
	add := func(p, role, member string, reads, bases int64) error {
		n, sum, err := runmanifest.HashFile(p)
		if err != nil {
			return err
		}
		m.Add(runmanifest.Output{Path: rel(p), Role: role, Member: member, Bytes: n, SHA256: sum, Reads: reads, Bases: bases})
		return nil
	}
	for _, f := range r.files {
		if err := add(f.path, runmanifest.RoleFASTQ, f.member, f.reads, f.bases); err != nil {
			return err
		}
	}
	if groundTruth != "" {
		if err := add(groundTruth, runmanifest.RoleTruth, "", 0, 0); err != nil {
			return err
		}
	}
	data, err := m.Marshal()
	if err != nil {
		return err
	}
	if err := os.WriteFile(path, data, 0o600); err != nil {
		return fmt.Errorf("writing run manifest: %w", err)
	}
	return nil
}

// Regenerate reproduces the reads of a scenario run without writing them
// and calls observe with each read's member and bases, in the order Run
// makes them. config is the run manifest's configuration with Seed set to
// the run seed (the derived seed of a fingerprinted run), fingerprinting
// off and the recorded input files configured (reference_path for the
// scenario reference); sc is the manifest's scenario record. It returns the
// number of reads.
func Regenerate(ctx context.Context, config types.SequencerConfig, sc *runmanifest.Scenario, observe func(member string, seq []byte) error) (int, error) {
	if sc == nil {
		return 0, fmt.Errorf("not a scenario run")
	}
	if observe == nil {
		return 0, fmt.Errorf("scenario regeneration needs an observer")
	}
	var spec specRecord
	if err := json.Unmarshal(sc.Spec, &spec); err != nil {
		return 0, fmt.Errorf("scenario record: %w", err)
	}
	if config.Seed == 0 {
		return 0, fmt.Errorf("scenario regeneration needs the run seed")
	}
	scenario := &ClinicalScenario{Phenotype: spec.Phenotype, Trio: spec.Trio, Sequencer: spec.Sequencer}
	for _, member := range sc.Members {
		switch member {
		case "proband":
			scenario.Output.ProbandFASTQ = member
		case "mother":
			scenario.Output.MotherFASTQ = member
		case "father":
			scenario.Output.FatherFASTQ = member
		default:
			return 0, fmt.Errorf("scenario record: unknown member %q", member)
		}
	}
	if scenario.Output.ProbandFASTQ == "" {
		return 0, fmt.Errorf("scenario record: no proband")
	}
	cfg := config
	g := NewScenarioGenerator(&cfg)
	g.observe = observe
	gt, err := g.Run(ctx, scenario)
	if err != nil {
		return 0, err
	}
	return gt.ComputeMetrics.ReadsGenerated, nil
}

// seqObserver is the FASTQ writer of a regenerated member: it passes the
// sequence line of every record to observe.
type seqObserver struct {
	member  string
	observe func(member string, seq []byte) error
	line    []byte
	lines   int
}

func (s *seqObserver) Write(p []byte) (int, error) {
	n := len(p)
	for len(p) > 0 {
		i := bytes.IndexByte(p, '\n')
		if i < 0 {
			s.line = append(s.line, p...)
			break
		}
		s.line = append(s.line, p[:i]...)
		p = p[i+1:]
		if s.lines%4 == 1 {
			if err := s.observe(s.member, bytes.TrimRight(s.line, "\r")); err != nil {
				return 0, err
			}
		}
		s.lines++
		s.line = s.line[:0]
	}
	return n, nil
}

// newEngine returns an engine for cfg, logging to stderr, or nowhere while
// regenerating.
func (g *ScenarioGenerator) newEngine(cfg types.SequencerConfig) *genomics.Engine {
	if g.observe != nil {
		return genomics.NewEngineLog(cfg, io.Discard)
	}
	return genomics.NewEngine(cfg)
}
