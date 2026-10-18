//go:build !release

package cli

import (
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// repoRoot is the absolute repository root (the loader rejects ".." paths).
var repoRoot = func() string {
	p, err := filepath.Abs("../..")
	if err != nil {
		panic(err)
	}
	return p
}()

func writeConfig(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestLoadSequencerConfig_PromethIONTemplate(t *testing.T) {
	c, err := NewConfigLoader(false).LoadSequencerConfig(filepath.Join(repoRoot, "configs/sequencers/ont_promethion24.yaml"))
	if err != nil {
		t.Fatal(err)
	}
	want := types.SequencerConfig{
		Platform:         "ont",
		SequencerModel:   "Oxford Nanopore Technologies PromethION 24",
		ChemistryVersion: "R10.4.1",
		MaxOutputTB:      14,
		SampleName:       TemplateDefaultSampleName,
		OutputDir:        TemplateDefaultOutputDir,
		ReadLength:       25000,
		ReadsPerBatch:    TemplateDefaultReadsPerBatch,
		TotalBatches:     TemplateDefaultTotalBatches,
		BaseQuality:      20,
		ReadLengthDistribution: types.ReadLengthSpec{
			Dist: types.ReadLengthLognormal, N50: 50000, Median: 25000, Min: 100, Max: 4000000,
		},
		ErrorModel: types.ErrorModelSpec{Tier: "sup"},
		ONT:        types.ONTSpec{Channels: 2675, SpeedBPS: 400},
		Output: types.OutputSpec{
			Compression: "gzip",
			FilePattern: "{sample_name}_{barcode}_{flowcell_id}_{file_number}.fastq.gz",
		},
	}
	if !reflect.DeepEqual(*c, want) {
		t.Errorf("loaded\n%+v\nwant\n%+v", *c, want)
	}

	// The template's output settings work with the engine. (Read lengths and
	// the error tier are #10/#11; use fixed short reads here.)
	c.OutputDir = t.TempDir()
	c.ReadLengthDistribution = types.ReadLengthSpec{}
	c.ErrorModel = types.ErrorModelSpec{}
	c.ReadLength, c.ReadsPerBatch, c.TotalBatches, c.Seed = 500, 5, 2, 1
	e := genomics.NewEngine(*c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	files := e.GeneratedFiles()
	if len(files) != 2 || !strings.HasPrefix(filepath.Base(files[1]), "sample_barcode01_PA") || !strings.HasSuffix(files[1], "_1.fastq.gz") {
		t.Errorf("GeneratedFiles = %v", files)
	}
	if _, err := os.Stat(filepath.Join(c.OutputDir, "sequencing_summary.txt")); err != nil {
		t.Errorf("no sequencing summary: %v", err)
	}
}

// Every shipped sequencer config loads to a non-zero config or fails with an
// explicit error; none silently loads as an empty config.
func TestLoadSequencerConfig_AllShippedConfigs(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join(repoRoot, "configs/sequencers/*.yaml"))
	v1, _ := filepath.Glob(filepath.Join(repoRoot, "configs/sequencer/v1.0/*.yaml"))
	if len(files) < 4 || len(v1) == 0 {
		t.Fatalf("found %d v2 and %d v1 configs", len(files), len(v1))
	}
	for _, f := range append(files, v1...) {
		t.Run(filepath.Base(f), func(t *testing.T) {
			c, err := NewConfigLoader(false).LoadSequencerConfig(f)
			switch {
			case strings.Contains(f, "ont_") || strings.Contains(f, "v1.0"):
				if err != nil {
					t.Fatalf("load failed: %v", err)
				}
				if c.ReadLength <= 0 || c.SampleName == "" || c.ReadsPerBatch <= 0 || c.TotalBatches <= 0 {
					t.Errorf("incomplete config %+v", *c)
				}
			default: // short-read and PacBio templates are not generated yet
				if err == nil {
					t.Fatalf("loaded %+v, want an explicit error", *c)
				}
				if !strings.Contains(err.Error(), "only Oxford Nanopore templates") {
					t.Errorf("error %q does not explain why", err)
				}
			}
		})
	}
}

func TestLoadSequencerConfig_TemplateRunSection(t *testing.T) {
	const tmpl = `version: "2.0"
platform: "ont_minion"
vendor: "Oxford Nanopore Technologies"
model: "MinION"
specifications:
  read_length: {min: 50, max: 100000, n50_target: 20000, median: 9000, distribution: log_normal}
  accuracy:
    default_basecaller: hac
    basecaller_modes:
      hac: { accuracy: 0.95, quality: 15 }
  throughput: {active_channels_per_fc: 512}
`
	l := NewConfigLoader(false)
	c, err := l.LoadSequencerConfig(writeConfig(t, "minion.yaml", tmpl))
	if err != nil {
		t.Fatal(err)
	}
	if c.ErrorModel.Tier != "hac" || c.BaseQuality != 15 || c.ONT.Channels != 512 || c.ONT.SpeedBPS != 0 ||
		c.ReadLengthDistribution.N50 != 20000 || c.Output.Compression != "" ||
		c.SampleName != TemplateDefaultSampleName || c.ReadsPerBatch != TemplateDefaultReadsPerBatch {
		t.Errorf("template without run section: %+v", *c)
	}

	c, err = l.LoadSequencerConfig(writeConfig(t, "minion_run.yaml", tmpl+`run:
  sample_name: HG002
  output_dir: /data/out
  reads_per_batch: 100
  total_batches: 7
  seed: 42
  ont: {flow_cell_id: FAQ12345}
  output: {summary_file: summary.txt}
  read_length_distribution: {dist: fixed}
  read_length: 5000
`))
	if err != nil {
		t.Fatal(err)
	}
	if c.SampleName != "HG002" || c.OutputDir != "/data/out" || c.ReadsPerBatch != 100 || c.TotalBatches != 7 || c.Seed != 42 ||
		c.ONT.FlowCellID != "FAQ12345" || c.ONT.Channels != 512 || c.Output.SummaryFile != "summary.txt" ||
		c.ReadLengthDistribution.Dist != "fixed" || c.ReadLength != 5000 || c.ErrorModel.Tier != "hac" {
		t.Errorf("run section not applied over the template: %+v", *c)
	}

	for name, content := range map[string]string{
		"typo in run":      tmpl + "run:\n  sample_nmae: x\n",
		"run not a map":    tmpl + "run: [1, 2]\n",
		"bad version":      strings.Replace(tmpl, `"2.0"`, `"3.0"`, 1),
		"bad dist":         strings.Replace(tmpl, "log_normal", "weibull", 1),
		"lognormal no n50": strings.Replace(tmpl, "n50_target: 20000, ", "", 1),
		"not ONT":          strings.NewReplacer("ont_minion", "acme_x", "Oxford Nanopore Technologies", "Acme").Replace(tmpl),
	} {
		if c, err := l.LoadSequencerConfig(writeConfig(t, "bad.yaml", content)); err == nil {
			t.Errorf("%s: loaded %+v, want an error", name, *c)
		}
	}
}

// Flat YAML uses the JSON field names; before, snake_case keys other than
// the few with yaml tags were silently dropped.
func TestLoadSequencerConfig_FlatYAMLUsesJSONKeys(t *testing.T) {
	c, err := NewConfigLoader(false).LoadSequencerConfig(writeConfig(t, "flat.yaml", `sample_name: S1
output_dir: ./out
read_length: 150
reads_per_batch: 10
total_batches: 2
paired_end: true
platform: ont
fingerprinting: {enable: true, steganography_layers: [quality_scores]}
ont: {channels: 512}
`))
	if err != nil {
		t.Fatal(err)
	}
	if c.SampleName != "S1" || c.ReadLength != 150 || c.ReadsPerBatch != 10 || !c.PairedEnd || c.Platform != "ont" ||
		!c.Fingerprinting.Enable || len(c.Fingerprinting.SteganographyLayers) != 1 || c.ONT.Channels != 512 {
		t.Errorf("flat YAML: %+v", *c)
	}
}

func TestLoadSequencerConfig_EmptyConfigsRejected(t *testing.T) {
	for name, content := range map[string]string{
		"empty.yaml":    "",
		"unknown.yaml":  "specs:\n  foo: 1\nreadlength: 3\n",
		"empty.json":    "{}",
		"unknown.json":  `{"readLength": 150}`,
		"empty_v1.yaml": "schema_version: \"1.0\"\n",
		"noext":         "foo: bar\n",
	} {
		if c, err := NewConfigLoader(false).LoadSequencerConfig(writeConfig(t, name, content)); err == nil {
			t.Errorf("%s: loaded %+v, want an error", name, *c)
		}
	}
}
