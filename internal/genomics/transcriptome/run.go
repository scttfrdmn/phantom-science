package transcriptome

// A transcriptome run as a whole: its parameters (RunSpec), the simulator
// they build, its run manifest (phantom-manifest.json, #31) and its
// regeneration for keyed-seed verification.
//
// Every random value of a run comes from the run seed: the default
// expression (stream 1), the novel isoforms (2), the reads (3), the read
// headers (4) and the run metadata (nanopore.RunInfo.WithDefaults). A
// keyed-seed fingerprint is therefore just the seed: the run is the
// unfingerprinted run whose seed is the derived seed.

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// RunSpec are the parameters of a transcriptome run, without its seed.
type RunSpec struct {
	Reads              int
	Tier               string
	Kit                string
	NovelRate          float64
	FullLengthFraction float64
	Zipf               ZipfModel
	SampleID           string
	Gzip               bool
	// Input files: TranscriptsPath, or GenomePath with GTFPath; GTFPath may
	// also annotate TranscriptsPath. ProfilePath is the expression table
	// ("" for the Zipf model).
	TranscriptsPath, GenomePath, GTFPath, ProfilePath string
}

// Simulator loads the inputs and returns the simulator of the run with
// seed; logw receives notes about the inputs (skipped records, profile
// matches).
func (s RunSpec) Simulator(seed uint64, logw io.Writer) (*Simulator, *Table, error) {
	kit, err := LookupKit(s.Kit)
	if err != nil {
		return nil, nil, err
	}
	if err := s.Zipf.Validate(); err != nil {
		return nil, nil, err
	}
	cfg := DefaultConfig(s.Reads, seed)
	cfg.Tier, cfg.Kit, cfg.NovelRate, cfg.FullLengthFraction = s.Tier, kit, s.NovelRate, s.FullLengthFraction
	cfg.Run.SampleID = s.SampleID
	if err := cfg.Validate(); err != nil {
		return nil, nil, err
	}
	if err := cfg.Run.WithDefaults(seed).Validate(); err != nil {
		return nil, nil, err
	}
	table, err := LoadTable(s.TranscriptsPath, s.GenomePath, s.GTFPath, logw)
	if err != nil {
		return nil, nil, err
	}
	var tpm []float64
	if s.ProfilePath != "" {
		var st ProfileStats
		tpm, st, err = LoadProfile(s.ProfilePath, table)
		if err != nil {
			return nil, nil, err
		}
		fmt.Fprintf(logw, "   [info] expression profile %s: %d rows, %d matched, %d not in the transcript table (values: %s)\n",
			s.ProfilePath, st.Rows, st.Matched, st.Unmatched, st.Column)
	} else {
		tpm = DefaultTPM(table, s.Zipf, seed)
	}
	sim, err := NewSimulator(table, tpm, cfg)
	if err != nil {
		return nil, nil, err
	}
	return sim, table, nil
}

// LoadTable loads the transcript table from a transcript FASTA (optionally
// annotated by a GTF) or from a genome FASTA and a GTF, reporting skipped
// records to logw.
func LoadTable(txPath, genomePath, gtfPath string, logw io.Writer) (*Table, error) {
	var ann *Annotation
	if gtfPath != "" {
		var err error
		if ann, err = ParseGTF(gtfPath); err != nil {
			return nil, err
		}
		if ann.Skipped > 0 {
			fmt.Fprintf(logw, "   [warn] %s: skipped %d unstranded transcripts\n", gtfPath, ann.Skipped)
		}
	}
	if genomePath != "" {
		table, st, err := BuildFromGenome(genomePath, ann)
		if err != nil {
			return nil, err
		}
		if st.MissingChrom > 0 {
			fmt.Fprintf(logw, "   [warn] %d annotated transcripts are on chromosomes missing from %s and were skipped\n", st.MissingChrom, genomePath)
		}
		return table, nil
	}
	if txPath == "" {
		return nil, errors.New("a transcript FASTA, or a genome FASTA with a GTF, is required")
	}
	table, err := LoadTranscriptFASTA(txPath)
	if err != nil {
		return nil, err
	}
	if ann != nil {
		att, mis := table.AttachAnnotation(ann)
		fmt.Fprintf(logw, "   [info] exon structure from %s for %d of %d transcripts", gtfPath, att, len(table.Transcripts))
		if mis > 0 {
			fmt.Fprintf(logw, " (%d matched by ID but differ in length)", mis)
		}
		fmt.Fprintln(logw)
	}
	return table, nil
}

// Record returns the manifest record of the run's parameters, hashing its
// input files.
func (s RunSpec) Record() (*manifest.Transcriptome, error) {
	r := &manifest.Transcriptome{
		Reads: s.Reads, Tier: s.Tier, Kit: s.Kit, NovelRate: s.NovelRate, FullLengthFraction: s.FullLengthFraction,
		ZipfAlpha: s.Zipf.Alpha, ExpressionSigma: s.Zipf.Sigma, SampleID: s.SampleID, Compression: "none",
	}
	if s.Gzip {
		r.Compression = "gzip"
	}
	for _, in := range []struct {
		path string
		dst  **manifest.Input
	}{
		{s.TranscriptsPath, &r.Transcripts}, {s.GenomePath, &r.Genome}, {s.GTFPath, &r.GTF}, {s.ProfilePath, &r.ExpressionProfile},
	} {
		if in.path == "" {
			continue
		}
		h, err := manifest.HashInput(in.path)
		if err != nil {
			return nil, err
		}
		*in.dst = h
	}
	return r, nil
}

// SpecFromRecord returns the parameters a run manifest records; the input
// paths are left empty for the caller to set to the recorded files.
func SpecFromRecord(r *manifest.Transcriptome) RunSpec {
	return RunSpec{
		Reads: r.Reads, Tier: r.Tier, Kit: r.Kit, NovelRate: r.NovelRate, FullLengthFraction: r.FullLengthFraction,
		Zipf: ZipfModel{Alpha: r.ZipfAlpha, Sigma: r.ExpressionSigma}, SampleID: r.SampleID, Gzip: r.Compression == "gzip",
	}
}

// WriteManifest writes the run manifest of a run written by WriteRun into
// res.Dir: the tool version, the seed (0 for a fingerprinted run, whose
// record fp is written instead), the run parameters rec, and every output
// with its size and SHA-256.
func WriteManifest(toolVersion string, seed uint64, fp *manifest.Fingerprint, rec *manifest.Transcriptome, sim *Simulator, res *RunResult) error {
	if fp != nil {
		seed = 0
	}
	m := manifest.New(toolVersion, seed, sim.RunInfo().Start, types.SequencerConfig{})
	m.Fingerprint = fp
	m.Transcriptome = rec
	info := sim.RunInfo()
	m.FlowCells = append(m.FlowCells, manifest.FlowCell{
		Index: 1, Seed: seed, RunID: info.RunID, FlowCellID: info.FlowCellID, Files: 1,
		Reads: int64(res.Reads), Bases: res.Bases,
	})
	outs := []struct{ name, role string }{
		{filepath.Base(res.FASTQ), manifest.RoleFASTQ}, {SummaryName, manifest.RoleSummary},
		{TruthName, manifest.RoleTruth}, {CountsName, manifest.RoleTruth}, {NovelTSVName, manifest.RoleTruth}, {NovelGTFName, manifest.RoleTruth},
	}
	for _, o := range outs {
		p := filepath.Join(res.Dir, o.name)
		n, sum, err := manifest.HashFile(p)
		if errors.Is(err, os.ErrNotExist) && o.name == NovelGTFName {
			continue // written only for isoforms with coordinates
		}
		if err != nil {
			return err
		}
		out := manifest.Output{Path: o.name, Role: o.role, FlowCell: 1, Bytes: n, SHA256: sum}
		if o.role == manifest.RoleFASTQ {
			out.Reads, out.Bases = int64(res.Reads), res.Bases
		}
		m.Add(out)
	}
	data, err := m.Marshal()
	if err != nil {
		return err
	}
	if err := os.WriteFile(filepath.Join(res.Dir, manifest.FileName), data, 0o600); err != nil {
		return fmt.Errorf("writing manifest: %w", err)
	}
	return nil
}

// Regenerate reproduces the reads of the run spec with seed, without
// writing them, and calls observe with every read's bases in order. It
// returns the number of reads.
func Regenerate(spec RunSpec, seed uint64, observe func(seq []byte) error) (int, error) {
	sim, _, err := spec.Simulator(seed, io.Discard)
	if err != nil {
		return 0, err
	}
	n := 0
	err = sim.Generate(func(rd *Read) error {
		n++
		return observe(rd.Record.Seq)
	})
	return n, err
}
