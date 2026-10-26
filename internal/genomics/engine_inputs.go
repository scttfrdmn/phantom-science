package genomics

// Input files of a run as the manifest records them (#37).
//
// A manifest carries no local paths: each file a run reads that changes its
// output (the reference FASTA, the variants VCF, the HMM model and the
// read-length histogram) is recorded by base name, size and SHA-256, and
// manifest.SanitizeConfig removes its path from the recorded configuration.
// The VCF is hashed when it is loaded (loadVariants); the other files are
// hashed once, when a run's manifest is written or a run is regenerated.

import (
	"fmt"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
)

// runInputs are the input files of a run; nil for an input not used.
type runInputs struct {
	reference, model, readLengths *manifest.Input
}

// hashInputs hashes the run's input files (other than the VCF), once.
func (e *Engine) hashInputs() (runInputs, error) {
	e.inputsOnce.Do(func() {
		hash := func(what, path string) *manifest.Input {
			if path == "" || e.inputsErr != nil {
				return nil
			}
			in, err := manifest.HashInput(path)
			if err != nil {
				e.inputsErr = fmt.Errorf("%s: %w", what, err)
			}
			return in
		}
		if e.ref != nil {
			e.inputs.reference = hash("reference", e.config.ReferencePath)
		}
		if e.hmmModel != nil {
			e.inputs.model = hash("HMM model", e.modelPath)
		}
		e.inputs.readLengths = hash("read-length histogram", e.config.ReadLengthDistribution.EmpiricalPath)
	})
	return e.inputs, e.inputsErr
}

// recordInputs records the run's input files in m.
func (e *Engine) recordInputs(m *manifest.Manifest) error {
	in, err := e.hashInputs()
	if err != nil {
		return err
	}
	m.Variants, m.Reference, m.HMMModel, m.ReadLengths = e.variantsInput, in.reference, in.model, in.readLengths
	return nil
}

// RecordInputs records the input files the engine reads (variants VCF,
// reference FASTA, HMM model, read-length histogram) in m, by base name,
// size and SHA-256, for commands that write their own run manifest
// (scenario).
func (e *Engine) RecordInputs(m *manifest.Manifest) error { return e.recordInputs(m) }

// checkRegenInputs reports an input file the run recorded in opts that the
// engine is not configured with: without it the run's reads cannot be
// reproduced. Inputs opts leaves nil are not checked.
func (e *Engine) checkRegenInputs(opts RegenOptions) error {
	in, err := e.hashInputs()
	if err != nil {
		return err
	}
	for _, c := range []struct {
		what, configure string
		want, got       *manifest.Input
	}{
		{"planted variants from", "that VCF (variants_path / --variants)", opts.Variants, e.variantsInput},
		{"sampled reads from reference", "that FASTA (reference_path / --reference)", opts.Reference, in.reference},
		{"sampled bases from HMM model", "that model (hmm_model_path / --hmm-model)", opts.HMMModel, in.model},
		{"sampled read lengths from histogram", "that histogram (read_length_distribution.empirical_path / --read-length-histogram)", opts.ReadLengths, in.readLengths},
	} {
		switch {
		case c.want == nil:
		case c.got == nil:
			return fmt.Errorf("the run %s %s (sha256 %s): configure %s to regenerate it", c.what, c.want.Name, c.want.SHA256, c.configure)
		case !c.got.Same(c.want):
			return fmt.Errorf("the run %s %s (sha256 %s), not %s (sha256 %s)", c.what, c.want.Name, c.want.SHA256, c.got.Name, c.got.SHA256)
		}
	}
	return nil
}
