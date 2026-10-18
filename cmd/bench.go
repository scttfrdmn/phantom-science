package cmd

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"text/tabwriter"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"github.com/spf13/cobra"
)

var benchCmd = &cobra.Command{
	Use:   "bench",
	Short: "Benchmark generation of a scenario across worker counts",
	Long: `Generate one scenario (platform, basecaller tier, read-length mix, total bases)
once per worker count and report reads/s, bases/s and MB/s for each.

The scenario is a sequencer config (--config) or a built-in one for --platform
(illumina: 150 bp reads; ont: lognormal lengths, N50 15 kb, median 8 kb), adjusted by the
same flags as generate. Every run uses the same seed, so all runs produce
identical data (the report checks that the manifests match). With --out the
data of each run is kept under <out>/workers-<n>/ with its phantom-manifest.json,
and the report is written to <out>/bench.json, so another benchmark can reuse
the same reproducible inputs. Throughput depends on the machine and its load.

  phantom-science bench --platform ont --basecaller-tier hac --total-bases 500M --workers 1,4,8`,
	RunE: func(cmd *cobra.Command, args []string) error {
		return runBench(cmd, cmd.OutOrStdout())
	},
}

func init() {
	rootCmd.AddCommand(benchCmd)
	f := benchCmd.Flags()
	f.String("config", "", "Sequencer configuration file (default: a built-in scenario for --platform)")
	f.String("platform", "", "Sequencing platform: illumina (default) or ont")
	f.String("basecaller-tier", "", "ONT basecaller tier for the error model: sup, hac or fast")
	f.String("read-length-dist", "", "Read-length distribution: fixed, lognormal, gamma or empirical")
	f.Int("read-length-n50", 0, "Read-length N50 for the lognormal distribution")
	f.Int("read-length-median", 0, "Median read length for the lognormal distribution")
	f.Float64("read-length-mean", 0, "Mean read length")
	f.Float64("read-length-sd", 0, "Read-length standard deviation")
	f.String("compress", "", "Output compression: none (default) or gzip")
	f.Uint64("seed", 1, "Random seed (the same for every run)")
	addScaleFlags(f)
	f.Lookup("total-bases").DefValue = "20M" // applied in runBench unless the config sets a target
	f.String("workers", "1,2,4,8", "Comma-separated worker counts to run")
	f.String("out", "", "Keep each run's output and the report (bench.json) under this directory")
	f.String("format", "table", "Report format: table or json")
}

// BenchScenario is the generation scenario of a bench report.
type BenchScenario struct {
	Platform      string               `json:"platform"`
	Tier          string               `json:"basecaller_tier,omitempty"`
	ReadLength    int                  `json:"read_length"`
	ReadLengths   types.ReadLengthSpec `json:"read_length_distribution,omitzero"`
	ReadsPerBatch int                  `json:"reads_per_batch"`
	PairedEnd     bool                 `json:"paired_end,omitempty"`
	Compression   string               `json:"compression,omitempty"`
	TotalBases    int64                `json:"total_bases"`
	FlowCells     int                  `json:"flow_cells,omitempty"`
	Seed          uint64               `json:"seed"`
}

// BenchResult is one run of a bench report.
type BenchResult struct {
	Workers        int     `json:"workers"`
	Seconds        float64 `json:"seconds"`
	Files          int     `json:"files"`
	Reads          int64   `json:"reads"`
	Bases          int64   `json:"bases"`
	Bytes          int64   `json:"bytes"`
	ReadsPerSecond float64 `json:"reads_per_second"`
	BasesPerSecond float64 `json:"bases_per_second"`
	MBPerSecond    float64 `json:"mb_per_second"`
	ManifestSHA256 string  `json:"manifest_sha256"`
	Dir            string  `json:"dir,omitempty"`
}

// BenchReport is the output of the bench command.
type BenchReport struct {
	ToolVersion string        `json:"tool_version"`
	Scenario    BenchScenario `json:"scenario"`
	Results     []BenchResult `json:"results"`
	// Identical reports that every run produced the same manifest (and so
	// the same output bytes).
	Identical bool `json:"identical"`
}

// benchConfig returns the built-in scenario for platform.
func benchConfig(platform string) types.SequencerConfig {
	c := types.SequencerConfig{
		SampleName: "bench", Lane: 1, ReadLength: 150, ReadsPerBatch: 20000, TotalBatches: 1,
		BaseQuality: 36, QualityDecay: 0.05, Platform: platform,
	}
	if strings.EqualFold(platform, "ont") {
		c.ReadLength, c.ReadsPerBatch = 8000, 500
		c.BaseQuality, c.QualityDecay = 20, 0
		c.ReadLengthDistribution = types.ReadLengthSpec{Dist: types.ReadLengthLognormal, N50: 15000, Median: 8000, Max: 200000}
	}
	return c
}

// parseWorkerCounts parses "1,2,4,8".
func parseWorkerCounts(s string) ([]int, error) {
	var out []int
	for _, f := range strings.Split(s, ",") {
		n, err := strconv.Atoi(strings.TrimSpace(f))
		if err != nil || n < 1 {
			return nil, fmt.Errorf("--workers %q: want comma-separated positive integers", s)
		}
		out = append(out, n)
	}
	return out, nil
}

// runBench implements the bench command, writing the report to out.
func runBench(cmd *cobra.Command, out io.Writer) error {
	flags := cmd.Flags()
	workersFlag, _ := flags.GetString("workers")
	counts, err := parseWorkerCounts(workersFlag)
	if err != nil {
		return err
	}
	format, _ := flags.GetString("format")
	if format != "table" && format != "json" {
		return fmt.Errorf("--format %q: want table or json", format)
	}
	keep, _ := flags.GetString("out")

	var config *types.SequencerConfig
	if path, _ := flags.GetString("config"); path != "" {
		if config, err = cli.NewConfigLoader(false).LoadSequencerConfig(path); err != nil {
			return fmt.Errorf("loading config: %w", err)
		}
	} else {
		platform, _ := flags.GetString("platform")
		c := benchConfig(platform)
		config = &c
	}
	if err := applyGenerateFlags(cmd, config); err != nil {
		return err
	}
	if !flags.Changed("total-bases") && config.TotalBases == 0 && config.TotalReads == 0 {
		config.TotalBases = 20e6
	}
	config.Seed, _ = flags.GetUint64("seed")
	if config.Seed == 0 {
		config.Seed = 1
	}
	config.Fingerprinting = types.FingerprintConfig{} // fingerprints embed a timestamp: runs would differ
	config.Enterprise.Fingerprinting = types.FingerprintConfig{}

	report := BenchReport{
		ToolVersion: Version,
		Scenario: BenchScenario{
			Platform: config.Platform, Tier: config.ErrorModel.Tier, ReadLength: config.ReadLength,
			ReadLengths: config.ReadLengthDistribution, ReadsPerBatch: config.ReadsPerBatch, PairedEnd: config.PairedEnd,
			Compression: config.Output.Compression, TotalBases: config.TotalBases, FlowCells: config.FlowCells, Seed: config.Seed,
		},
		Identical: true,
	}
	if report.Scenario.Platform == "" {
		report.Scenario.Platform = "illumina"
	}
	for _, n := range counts {
		res, err := benchRun(*config, n, keep)
		if err != nil {
			return fmt.Errorf("workers %d: %w", n, err)
		}
		if len(report.Results) > 0 && res.ManifestSHA256 != report.Results[0].ManifestSHA256 {
			report.Identical = false
		}
		report.Results = append(report.Results, res)
	}

	data, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		return err
	}
	data = append(data, '\n')
	if keep != "" {
		if err := os.WriteFile(filepath.Join(keep, "bench.json"), data, 0o600); err != nil {
			return fmt.Errorf("writing bench report: %w", err)
		}
	}
	if format == "json" {
		_, err := out.Write(data)
		return err
	}
	return writeBenchTable(out, report)
}

// benchRun generates config once with n workers and measures it.
func benchRun(config types.SequencerConfig, n int, keep string) (BenchResult, error) {
	res := BenchResult{Workers: n}
	if keep != "" {
		res.Dir = filepath.Join(keep, fmt.Sprintf("workers-%d", n))
		if err := os.RemoveAll(res.Dir); err != nil {
			return res, err
		}
		config.OutputDir = res.Dir
	} else {
		dir, err := os.MkdirTemp("", "phantom-bench-")
		if err != nil {
			return res, err
		}
		defer os.RemoveAll(dir)
		config.OutputDir = dir
	}
	config.Workers = n

	engine := genomics.NewEngine(config)
	defer engine.Close()
	engine.SetToolVersion(Version)
	restore := silenceStdout()
	start := time.Now()
	err := engine.Generate()
	elapsed := time.Since(start)
	restore()
	if err != nil {
		return res, err
	}
	m := engine.Manifest()
	raw, err := os.ReadFile(filepath.Join(config.OutputDir, manifest.FileName))
	if err != nil {
		return res, err
	}
	sum := sha256.Sum256(raw)
	res.ManifestSHA256 = hex.EncodeToString(sum[:])
	res.Seconds = elapsed.Seconds()
	res.Files, res.Reads, res.Bases, res.Bytes = m.Totals.Files, m.Totals.Reads, m.Totals.Bases, m.Totals.Bytes
	if s := res.Seconds; s > 0 {
		res.ReadsPerSecond = float64(res.Reads) / s
		res.BasesPerSecond = float64(res.Bases) / s
		res.MBPerSecond = float64(res.Bytes) / 1e6 / s
	}
	return res, nil
}

// silenceStdout redirects the engine's progress output (written to
// os.Stdout) to the null device and returns a function restoring it.
func silenceStdout() func() {
	null, err := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
	if err != nil {
		return func() {}
	}
	old := os.Stdout
	os.Stdout = null
	return func() {
		os.Stdout = old
		null.Close()
	}
}

// writeBenchTable writes the report as a table.
func writeBenchTable(out io.Writer, r BenchReport) error {
	s := r.Scenario
	dist := s.ReadLengths.Dist
	if dist == "" {
		dist = "fixed"
	}
	fmt.Fprintf(out, "scenario: platform=%s tier=%s lengths=%s read_length=%d total_bases=%d seed=%d\n",
		s.Platform, orDash(s.Tier), dist, s.ReadLength, s.TotalBases, s.Seed)
	tw := tabwriter.NewWriter(out, 0, 0, 2, ' ', tabwriter.AlignRight)
	fmt.Fprintln(tw, "workers\tseconds\tfiles\treads\tbases\tMB\treads/s\tbases/s\tMB/s\t")
	for _, x := range r.Results {
		fmt.Fprintf(tw, "%d\t%.2f\t%d\t%d\t%d\t%.1f\t%.0f\t%.0f\t%.1f\t\n",
			x.Workers, x.Seconds, x.Files, x.Reads, x.Bases, float64(x.Bytes)/1e6, x.ReadsPerSecond, x.BasesPerSecond, x.MBPerSecond)
	}
	if err := tw.Flush(); err != nil {
		return err
	}
	if r.Identical {
		fmt.Fprintln(out, "outputs identical across worker counts: yes")
	} else {
		fmt.Fprintln(out, "outputs identical across worker counts: NO")
	}
	return nil
}

func orDash(s string) string {
	if s == "" {
		return "-"
	}
	return s
}
