package cli

// Input files of signal and transcriptome candidates (#31).
//
// A signal manifest records the input reads, the pore model and the --truth
// file of the run, and a transcriptome manifest its transcript FASTA,
// genome FASTA, GTF and expression table, each by base name, size and
// SHA-256. As for generate (verify_inputs.go), verify is given the files
// (--signal-input, --pore-model, --signal-truth; --transcripts, --genome,
// --gtf, --expression-profile) and uses, for each candidate, the one whose
// size and SHA-256 match the record; a candidate whose file is not given
// fails.

import (
	"fmt"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/signalrun"
	"github.com/scttfrdmn/phantom-science/internal/genomics/transcriptome"
)

// suppliedFiles hashes the files given with each flag once.
type suppliedFiles map[string][]suppliedFile

func (s suppliedFiles) get(flag string, paths []string) ([]suppliedFile, error) {
	if have, ok := s[flag]; ok {
		return have, nil
	}
	have := []suppliedFile{}
	for _, p := range paths {
		in, err := manifest.HashInput(p)
		if err != nil {
			return nil, fmt.Errorf("%s: %w", flag, err)
		}
		have = append(have, suppliedFile{p, in})
	}
	s[flag] = have
	return have, nil
}

// pick returns the path of the file given with flag that matches want.
func (s suppliedFiles) pick(c *candidate, want *manifest.Input, flag string, paths []string, file, verb string) (string, error) {
	have, err := s.get(flag, paths)
	if err != nil {
		return "", err
	}
	for _, h := range have {
		if h.in.Same(want) {
			return h.path, nil
		}
	}
	return "", fmt.Errorf("%s: the run %s %s (%d bytes, sha256 %s): pass that %s with %s", c.label, verb, want.Name, want.Bytes, want.SHA256, file, flag)
}

// resolveRunFiles points each signal and transcriptome candidate at the
// given files its manifest records.
func resolveRunFiles(cands []*candidate, opts VerifyOptions) error {
	s := suppliedFiles{}
	for _, c := range cands {
		switch c.kind {
		case kindSignal:
			r := c.sig
			spec := signalrun.RegenSpec{Profile: r.Profile, Sets: r.ProfileOverrides, ModelUnits: r.ModelUnits}
			for i := range r.Inputs {
				p, err := s.pick(c, &r.Inputs[i], "--signal-input", opts.SignalInputs, "read file", "synthesised signal from")
				if err != nil {
					return err
				}
				spec.Inputs = append(spec.Inputs, p)
			}
			p, err := s.pick(c, &r.PoreModel, "--pore-model", opts.PoreModels, "pore model", "used the pore model")
			if err != nil {
				return err
			}
			spec.ModelPath = p
			if r.Truth != nil {
				if spec.TruthPath, err = s.pick(c, r.Truth, "--signal-truth", opts.SignalTruths, "truth file", "took source sequences from"); err != nil {
					return err
				}
			}
			c.sigSpec = spec
		case kindTranscriptome:
			r := c.tx
			spec := transcriptome.SpecFromRecord(r)
			for _, in := range []struct {
				want       *manifest.Input
				flag       string
				paths      []string
				file, verb string
				dst        *string
			}{
				{r.Transcripts, "--transcripts", opts.Transcripts, "transcript FASTA", "sampled molecules from", &spec.TranscriptsPath},
				{r.Genome, "--genome", opts.Genomes, "genome FASTA", "spliced transcripts from", &spec.GenomePath},
				{r.GTF, "--gtf", opts.GTFs, "GTF", "took exon structure from", &spec.GTFPath},
				{r.ExpressionProfile, "--expression-profile", opts.ExpressionProfiles, "expression table", "took expression from", &spec.ProfilePath},
			} {
				if in.want == nil {
					continue
				}
				p, err := s.pick(c, in.want, in.flag, in.paths, in.file, in.verb)
				if err != nil {
					return err
				}
				*in.dst = p
			}
			if spec.TranscriptsPath == "" && spec.GenomePath == "" {
				return fmt.Errorf("%s: the transcriptome manifest records neither transcripts nor a genome", c.label)
			}
			c.txSpec = spec
		}
	}
	return nil
}
