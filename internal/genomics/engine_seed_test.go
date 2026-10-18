//go:build !release

package genomics

import (
	"bytes"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

func streamOutput(t *testing.T, c types.SequencerConfig, opts GenerateOpts) []byte {
	t.Helper()
	e := NewEngine(c)
	defer e.Close()
	var buf bytes.Buffer
	if err := e.GenerateToWriter(&buf, opts); err != nil {
		t.Fatal(err)
	}
	return buf.Bytes()
}

func fileOutputs(t *testing.T, c types.SequencerConfig) map[string][]byte {
	t.Helper()
	c.OutputDir = t.TempDir()
	e := NewEngine(c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	out := map[string][]byte{}
	for _, f := range e.GeneratedFiles() {
		data, err := os.ReadFile(f)
		if err != nil {
			t.Fatal(err)
		}
		out[filepath.Base(f)] = data
	}
	return out
}

func TestSeed_StreamByteIdentical(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.PairedEnd = true
	c.Seed = 42
	opts := GenerateOpts{
		SampleID:  "s1",
		RefRegion: &GenomicRegion{Chromosome: "chr2", Start: 1000, End: 1400},
		Variants: []VariantRecord{
			{Chromosome: "chr2", Position: 1100, Ref: "A", Alt: "T", Freq: 0.5},
			{Chromosome: "chr2", Position: 1150, Ref: "ACGT", Alt: "A", Freq: 0.3},
		},
	}
	a := streamOutput(t, c, opts)
	b := streamOutput(t, c, opts)
	if !bytes.Equal(a, b) {
		t.Fatal("same seed produced different streamed output")
	}
	c.Seed = 43
	if bytes.Equal(a, streamOutput(t, c, opts)) {
		t.Error("different seeds produced identical output")
	}
}

func TestSeed_FilesByteIdentical(t *testing.T) {
	c := goldenConfig("")
	c.PairedEnd = true
	c.Seed = 7
	a, b := fileOutputs(t, c), fileOutputs(t, c)
	if len(a) != 4 {
		t.Fatalf("%d files", len(a))
	}
	for name, data := range a {
		if !bytes.Equal(data, b[name]) {
			t.Errorf("%s differs between runs with the same seed", name)
		}
	}
	// R1 and R2 of a batch come from different streams.
	if bytes.Equal(a["Gold_S1_L003_R1_001.fastq"][100:], a["Gold_S1_L003_R2_001.fastq"][100:]) {
		t.Error("R1 and R2 have identical content")
	}
}

// A file's content depends only on (seed, batch, mate), not on what else the
// engine generated before: the basis for parallel generation (#16).
func TestSeed_FileIndependentOfOrder(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Seed = 9
	all := fileOutputs(t, c)

	e := NewEngine(c)
	defer e.Close()
	path := filepath.Join(c.OutputDir, "only2.fastq")
	if err := e.GenerateFASTQFile(path, 2, 1); err != nil {
		t.Fatal(err)
	}
	data, _ := os.ReadFile(path)
	if !bytes.Equal(data, all["Gold_S1_L003_002.fastq"]) {
		t.Error("batch 2 generated alone differs from batch 2 of a full run")
	}
}

func TestSeed_HMMByteIdentical(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(3)
	if err := model.TrainFromSequence(makeEngineTestSequence(5000)); err != nil {
		t.Fatal(err)
	}
	modelPath := filepath.Join(t.TempDir(), "m.hmm.json")
	if err := model.Save(modelPath); err != nil {
		t.Fatal(err)
	}
	c := goldenConfig(t.TempDir())
	c.HMMModelPath = modelPath
	c.Seed = 5
	a := streamOutput(t, c, GenerateOpts{})
	b := streamOutput(t, c, GenerateOpts{})
	if !bytes.Equal(a, b) {
		t.Error("same seed with an HMM model produced different output")
	}
}

func TestSeed_UnsetIsDrawnAndLogged(t *testing.T) {
	c := goldenConfig(t.TempDir())
	e := NewEngine(c)
	defer e.Close()
	var log bytes.Buffer
	e.logw = &log
	if e.Seed() == 0 {
		t.Fatal("seed not drawn")
	}
	if err := e.GenerateToWriter(io.Discard, GenerateOpts{}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(log.String(), "--seed ") {
		t.Errorf("drawn seed not logged: %q", log.String())
	}
	// Reusing the logged seed reproduces the run.
	c.Seed = e.Seed()
	first := streamOutput(t, c, GenerateOpts{})
	var again bytes.Buffer
	if err := e.GenerateToWriter(&again, GenerateOpts{}); err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(first, again.Bytes()) {
		t.Error("explicit seed did not reproduce the drawn-seed run")
	}
}

func TestSeed_ExplicitNotLogged(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Seed = 1
	e := NewEngine(c)
	defer e.Close()
	var log bytes.Buffer
	e.logw = &log
	if err := e.GenerateToWriter(io.Discard, GenerateOpts{}); err != nil {
		t.Fatal(err)
	}
	if log.Len() != 0 {
		t.Errorf("unexpected log: %q", log.String())
	}
}

// NewEngine must not write to stdout: `generate --stdout` streams FASTQ there.
func TestNewEngine_HMMMessagesNotOnStdout(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(3)
	if err := model.TrainFromSequence(makeEngineTestSequence(2000)); err != nil {
		t.Fatal(err)
	}
	modelPath := filepath.Join(t.TempDir(), "m.hmm.json")
	if err := model.Save(modelPath); err != nil {
		t.Fatal(err)
	}
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	old := os.Stdout
	os.Stdout = w
	c := goldenConfig(t.TempDir())
	c.HMMModelPath = modelPath
	e1 := NewEngine(c)
	c.HMMModelPath = filepath.Join(t.TempDir(), "missing.hmm.json")
	e2 := NewEngine(c)
	var buf bytes.Buffer
	err = e1.GenerateToWriter(&buf, GenerateOpts{})
	os.Stdout = old
	w.Close()
	e1.Close()
	e2.Close()
	if err != nil {
		t.Fatal(err)
	}
	out, _ := io.ReadAll(r)
	if len(out) != 0 {
		t.Errorf("NewEngine/GenerateToWriter wrote to stdout: %q", out)
	}
}

func TestNotImplemented_Errors(t *testing.T) {
	cases := []struct {
		name  string
		edit  func(*types.SequencerConfig)
		issue string
	}{
		{"lognormal", func(c *types.SequencerConfig) { c.ReadLengthDistribution.Dist = "lognormal" }, "#10"},
		{"gamma", func(c *types.SequencerConfig) { c.ReadLengthDistribution.Dist = "gamma" }, "#10"},
		{"empirical", func(c *types.SequencerConfig) { c.ReadLengthDistribution.Dist = "empirical" }, "#10"},
		{"ont", func(c *types.SequencerConfig) { c.Platform = "ont" }, "#12"},
		{"gzip", func(c *types.SequencerConfig) { c.Output.Compression = "gzip" }, "#12"},
		{"summary", func(c *types.SequencerConfig) { c.Output.SummaryFile = "sequencing_summary.txt" }, "#12"},
		{"pattern", func(c *types.SequencerConfig) { c.Output.FilePattern = "{sample}.fastq" }, "#16"},
		{"total-bases", func(c *types.SequencerConfig) { c.TotalBases = 1e9 }, "#16"},
		{"total-reads", func(c *types.SequencerConfig) { c.TotalReads = 1000 }, "#16"},
		{"flow-cells", func(c *types.SequencerConfig) { c.FlowCells = 2 }, "#16"},
		{"workers", func(c *types.SequencerConfig) { c.Workers = 8 }, "#16"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			c := goldenConfig(t.TempDir())
			tc.edit(&c)
			e := NewEngine(c)
			defer e.Close()
			for what, err := range map[string]error{
				"GenerateToWriter":  e.GenerateToWriter(io.Discard, GenerateOpts{}),
				"Generate":          e.Generate(),
				"GenerateFASTQFile": e.GenerateFASTQFile(filepath.Join(c.OutputDir, "x.fastq"), 1, 1),
			} {
				if !errors.Is(err, ErrNotImplemented) || !strings.Contains(err.Error(), tc.issue) {
					t.Errorf("%s: err = %v, want ErrNotImplemented mentioning %s", what, err, tc.issue)
				}
			}
			if files, _ := os.ReadDir(c.OutputDir); len(files) != 0 {
				t.Errorf("output written despite the error: %d files", len(files))
			}
		})
	}
}

func TestDefaults_Accepted(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Platform = "illumina"
	c.ReadLengthDistribution.Dist = "fixed"
	c.Output.Compression = "none"
	c.Workers = 1
	want := streamOutput(t, func() types.SequencerConfig { d := goldenConfig(c.OutputDir); d.Seed = 3; return d }(), GenerateOpts{})
	c.Seed = 3
	if got := streamOutput(t, c, GenerateOpts{}); !bytes.Equal(got, want) {
		t.Error("explicit default settings changed the output")
	}
}

func TestUnknownValues_Rejected(t *testing.T) {
	for name, edit := range map[string]func(*types.SequencerConfig){
		"platform": func(c *types.SequencerConfig) { c.Platform = "pacbio-hifi" },
		"dist":     func(c *types.SequencerConfig) { c.ReadLengthDistribution.Dist = "weibull" },
		"tier":     func(c *types.SequencerConfig) { c.ErrorModel.Tier = "ultra" },
		"compress": func(c *types.SequencerConfig) { c.Output.Compression = "zstd" },
	} {
		c := goldenConfig(t.TempDir())
		edit(&c)
		e := NewEngine(c)
		err := e.GenerateToWriter(io.Discard, GenerateOpts{})
		e.Close()
		if err == nil || errors.Is(err, ErrNotImplemented) || !strings.Contains(err.Error(), "unknown") {
			t.Errorf("%s: err = %v, want an 'unknown ...' error", name, err)
		}
	}
}
