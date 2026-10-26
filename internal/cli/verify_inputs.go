package cli

// Input files of candidate runs (#32, #37).
//
// A manifest records the files a run read (VCF, reference FASTA, HMM model,
// read-length histogram) by base name, size and SHA-256, never by path, so
// verify is given the files (--variants, --reference, --hmm-model,
// --read-length-histogram) and uses, for each candidate, the one whose size
// and SHA-256 match the record. A candidate whose file is not given fails,
// rather than being regenerated without it and reported as not matching.
//
// Manifests of phantom-science 0.14.0 and earlier record the reference, HMM
// model and histogram by path in their configuration instead. For those,
// the given file with the recorded base name (or the only file given) is
// used, else the recorded path; either way with a warning, because nothing
// in the manifest can check the file.

import (
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// runInputs are the input files a candidate's manifest records; nil for
// none.
type runInputs struct {
	variants, reference, model, readLengths *manifest.Input
}

// inputKind is one kind of input file.
type inputKind struct {
	flag string // the verify flag that supplies it
	file string // what the file is, for messages
	verb string // what the run did with it
	// recorded returns the candidate's record of the input.
	recorded func(*runInputs) *manifest.Input
	// supplied returns the files given for it.
	supplied func(VerifyOptions) []string
	// use configures c to read the file at path.
	use func(c *types.SequencerConfig, path string)
	// legacyPath returns the path an older manifest's configuration
	// records ("" for none); nil for an input always recorded by hash.
	legacyPath func(c types.SequencerConfig) string
}

var inputKinds = []inputKind{
	{
		flag: "--variants", file: "VCF", verb: "planted variants from",
		recorded: func(r *runInputs) *manifest.Input { return r.variants },
		supplied: func(o VerifyOptions) []string { return o.Variants },
		use:      func(c *types.SequencerConfig, p string) { c.VariantsPath = p },
	},
	{
		flag: "--reference", file: "FASTA", verb: "sampled reads from",
		recorded:   func(r *runInputs) *manifest.Input { return r.reference },
		supplied:   func(o VerifyOptions) []string { return o.References },
		use:        func(c *types.SequencerConfig, p string) { c.ReferencePath = p },
		legacyPath: func(c types.SequencerConfig) string { return c.ReferencePath },
	},
	{
		flag: "--hmm-model", file: "HMM model", verb: "sampled bases from",
		recorded: func(r *runInputs) *manifest.Input { return r.model },
		supplied: func(o VerifyOptions) []string { return o.HMMModels },
		use:      func(c *types.SequencerConfig, p string) { c.HMMModelPath, c.ModelDir = p, "" },
		legacyPath: func(c types.SequencerConfig) string {
			switch {
			case c.HMMModelPath != "":
				return c.HMMModelPath
			case c.ModelDir != "" && c.Population != "":
				return filepath.Join(c.ModelDir, c.Population+".hmm.json")
			}
			return ""
		},
	},
	{
		flag: "--read-length-histogram", file: "read-length histogram", verb: "sampled read lengths from",
		recorded:   func(r *runInputs) *manifest.Input { return r.readLengths },
		supplied:   func(o VerifyOptions) []string { return o.ReadLengths },
		use:        func(c *types.SequencerConfig, p string) { c.ReadLengthDistribution.EmpiricalPath = p },
		legacyPath: func(c types.SequencerConfig) string { return c.ReadLengthDistribution.EmpiricalPath },
	},
}

// suppliedFile is a file given for an input, hashed once.
type suppliedFile struct {
	path string
	in   *manifest.Input
}

// resolveCandidateInputs points each manifest candidate at the given files
// its manifest records (see the file comment); w receives warnings.
func resolveCandidateInputs(w io.Writer, cands []*candidate, opts VerifyOptions) error {
	have := make([][]suppliedFile, len(inputKinds))
	for i, k := range inputKinds {
		for _, p := range k.supplied(opts) {
			in, err := manifest.HashInput(p)
			if err != nil {
				return fmt.Errorf("%s: %w", k.flag, err)
			}
			have[i] = append(have[i], suppliedFile{p, in})
		}
	}
	for _, c := range cands {
		if !c.fromManifest {
			continue // a --config candidate names its own files
		}
		for i, k := range inputKinds {
			if err := k.resolve(w, c, have[i]); err != nil {
				return err
			}
		}
	}
	return resolveRunFiles(cands, opts)
}

// resolve configures c with its input of kind k from have.
func (k inputKind) resolve(w io.Writer, c *candidate, have []suppliedFile) error {
	if want := k.recorded(&c.inputs); want != nil {
		for _, h := range have {
			if h.in.Same(want) {
				k.use(&c.config, h.path)
				return nil
			}
		}
		return fmt.Errorf("%s: the run %s %s (%d bytes, sha256 %s): pass that %s with %s", c.label, k.verb, want.Name, want.Bytes, want.SHA256, k.file, k.flag)
	}
	if k.legacyPath == nil {
		return nil
	}
	old := k.legacyPath(c.config)
	if old == "" {
		return nil
	}
	use := ""
	for _, h := range have {
		if h.in.Name == filepath.Base(old) {
			use = h.path
			break
		}
	}
	switch {
	case use != "":
	case len(have) == 1:
		use = have[0].path
	case len(have) > 1:
		return fmt.Errorf("%s: the manifest (phantom-science 0.14.0 or earlier) records the %s by path only, %s, and none of the %d files given with %s has that name", c.label, k.file, old, len(have), k.flag)
	default:
		if _, err := os.Stat(old); err != nil {
			return fmt.Errorf("%s: the manifest (phantom-science 0.14.0 or earlier) records the %s by path only, %s, which cannot be read (%v): pass it with %s", c.label, k.file, old, err, k.flag)
		}
		use = old
	}
	fmt.Fprintf(w, "   [warn] %s: the manifest (phantom-science 0.14.0 or earlier) records the %s by path, not by sha256: using %s, which cannot be checked against the run\n", c.label, k.file, use)
	k.use(&c.config, use)
	return nil
}
