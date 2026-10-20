// Command run_e2e_tests generates reads for each platform profile of the
// internal/testing E2E framework with the real generator and reports what
// it measures: read statistics, generation benchmarks, an ML classifier
// calibration and, when real reads are supplied, a real-versus-synthetic
// classifier accuracy. It prints no scores or verdicts; it exits non-zero
// only when a step fails or the calibration fails.
package main

import (
	"encoding/csv"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/testing"
)

// CommandLineArgs holds command line arguments.
type CommandLineArgs struct {
	platforms    string
	outputDir    string
	seed         uint64
	realFASTQ    string
	verbose      bool
	quick        bool
	benchmark    bool
	report       bool
	authenticity bool
	help         bool
}

func main() {
	args := parseCommandLineArgs()
	if args.help {
		printUsage()
		return
	}
	if args.realFASTQ == "" {
		args.realFASTQ = os.Getenv("PHANTOM_REAL_FASTQ")
	}
	if err := validateArgs(args); err != nil {
		log.Fatalf("invalid arguments: %v", err)
	}

	report := runAll(testing.NewE2EFramework(), parsePlatformList(args.platforms), args)
	printReport(os.Stdout, report, args)

	if args.report {
		files, err := writeReports(report, args.outputDir)
		if err != nil {
			log.Fatalf("write reports: %v", err)
		}
		for _, f := range files {
			fmt.Printf("wrote %s\n", f)
		}
	}
	if n := report.Failures(); n > 0 {
		log.Fatalf("%d failure(s); see the table above", n)
	}
}

// parseCommandLineArgs parses command line arguments.
func parseCommandLineArgs() *CommandLineArgs {
	args := &CommandLineArgs{}
	flag.StringVar(&args.platforms, "platforms", "all",
		"Comma-separated list of platforms (all, illumina, longread, pacbio, ont, novaseq, nextseq, miseq, revio, promethion, dnbseq, aviti)")
	flag.StringVar(&args.outputDir, "output", "./e2e_results", "Directory for the JSON and CSV reports")
	flag.Uint64Var(&args.seed, "seed", testing.DefaultSeed, "Base generation seed (non-zero); every dataset's seed is derived from it")
	flag.StringVar(&args.realFASTQ, "real-fastq", "",
		"FASTQ of real reads to compare generated reads with (default $PHANTOM_REAL_FASTQ)")
	flag.BoolVar(&args.verbose, "verbose", false, "Enable verbose logging")
	flag.BoolVar(&args.quick, "quick", false, "Use fewer reads for the calibration")
	flag.BoolVar(&args.benchmark, "benchmark", true, "Measure generation time, CPU and allocation")
	flag.BoolVar(&args.report, "report", true, "Write e2e_results.json and e2e_summary.csv to -output")
	flag.BoolVar(&args.authenticity, "authenticity", true,
		"Run the ML classifier calibration (null comparison and GC-skew positive control)")
	flag.BoolVar(&args.help, "help", false, "Show this help message")
	flag.Parse()
	return args
}

// printUsage prints usage information.
func printUsage() {
	fmt.Println(`Usage: run_e2e_tests [OPTIONS]

Generates reads for each selected platform profile with the phantom-science
engine (seeded, so runs are reproducible) and reports measured values only.

OPTIONS:`)
	flag.PrintDefaults()
	fmt.Printf(`
PLATFORMS:
  The engine models Illumina (novaseq, nextseq, miseq) and Oxford Nanopore
  (ont, promethion). PacBio, MGI and Element AVITI profiles are listed as
  NOT MODELED and skipped; they are not failures.

WHAT IS REPORTED:
  Reads    read count, bases, mean/min/max read length, N50, mean Phred
           quality and GC content, computed from the generated FASTQ.
  Bench    (-benchmark) wall time, process CPU seconds and share of the
           machine, heap bytes allocated, and bases per second.
  Calib    (-authenticity) a Multinomial Naive Bayes classifier on 4-mer
           counts, trained on half of each set and scored on the other half:
             null     two independent generator draws; passes when accuracy
                      is within 0.5 +/- %.0f sigma (binomial).
             control  the second draw GC-skewed to %.2f; passes when
                      accuracy exceeds %.2f, showing the classifier has power.
           %d reads per set (%d with -quick).
  Real     (-real-fastq or $PHANTOM_REAL_FASTQ) the same classifier on real
           reads versus as many generated reads. Reported, not graded.

  The command exits non-zero when generation, benchmarking or report writing
  fails, or when a calibration check fails.
`, nullBandSigmas, gcControlTarget, positiveControlMin, calibrationReads, quickCalibrationReads)
}

// validateArgs validates command line arguments.
func validateArgs(args *CommandLineArgs) error {
	if args.report && args.outputDir == "" {
		return fmt.Errorf("output directory cannot be empty")
	}
	if len(parsePlatformList(args.platforms)) == 0 {
		return fmt.Errorf("no valid platforms specified")
	}
	if args.seed == 0 {
		return fmt.Errorf("seed must be non-zero")
	}
	if args.realFASTQ != "" {
		if _, err := os.Stat(args.realFASTQ); err != nil {
			return fmt.Errorf("real FASTQ: %w", err)
		}
	}
	return nil
}

// parsePlatformList parses the platform list argument.
func parsePlatformList(platformsArg string) []testing.SequencingPlatform {
	platformMap := map[string][]testing.SequencingPlatform{
		"all": testing.AllPlatforms,
		"illumina": {
			testing.IlluminaNovaSeq6000,
			testing.IlluminaNextSeq2000,
			testing.IlluminaMiSeq,
		},
		"longread":   {testing.PacBioRevio, testing.ONTPromethION24},
		"pacbio":     {testing.PacBioRevio},
		"ont":        {testing.ONTPromethION24},
		"novaseq":    {testing.IlluminaNovaSeq6000},
		"nextseq":    {testing.IlluminaNextSeq2000},
		"miseq":      {testing.IlluminaMiSeq},
		"revio":      {testing.PacBioRevio},
		"promethion": {testing.ONTPromethION24},
		"dnbseq":     {testing.MGI_DNBSEQ_T20x2},
		"aviti":      {testing.ElementAVITI},
	}

	seen := make(map[testing.SequencingPlatform]bool)
	var platforms []testing.SequencingPlatform
	for _, part := range strings.Split(strings.ToLower(platformsArg), ",") {
		for _, p := range platformMap[strings.TrimSpace(part)] {
			if !seen[p] {
				seen[p] = true
				platforms = append(platforms, p)
			}
		}
	}
	return platforms
}

// Calibration parameters. Accuracy on n held-out reads under the null is
// Binomial(n, 0.5)/n, so its standard deviation is sqrt(0.25/n).
const (
	// calibrationReads is the reads per set; half train, half are held out.
	calibrationReads      = 2000
	quickCalibrationReads = 400
	// nullBandSigmas sets the width of the pass band around 0.5. At 4 sigma
	// a correct pipeline fails by chance with probability below 1e-4.
	nullBandSigmas = 4.0
	// gcControlTarget is the GC content the positive-control set is skewed
	// to (the generator's reads are about 0.5).
	gcControlTarget = 0.7
	// positiveControlMin is the accuracy the positive control must exceed.
	positiveControlMin = 0.9
)

// Status values for PlatformResult.
const (
	statusMeasured   = "measured"
	statusNotModeled = "not_modeled"
	statusError      = "error"
)

// Calibration holds the ML calibration measured on one platform.
type Calibration struct {
	ReadsPerSet     int     `json:"reads_per_set"`
	Holdout         int     `json:"holdout_reads"`
	NullAccuracy    float64 `json:"null_accuracy"`
	Band            float64 `json:"null_band"`
	ControlTargetGC float64 `json:"control_target_gc"`
	ControlAccuracy float64 `json:"control_accuracy"`
	NullPass        bool    `json:"null_pass"`
	ControlPass     bool    `json:"control_pass"`
}

// Passed reports whether both checks passed.
func (c *Calibration) Passed() bool { return c.NullPass && c.ControlPass }

// gradeCalibration grades null and positive-control accuracies measured on
// holdout reads.
func gradeCalibration(readsPerSet, holdout int, nullAccuracy, controlAccuracy float64) *Calibration {
	band := testing.ChanceBand(holdout, nullBandSigmas)
	return &Calibration{
		ReadsPerSet:     readsPerSet,
		Holdout:         holdout,
		NullAccuracy:    nullAccuracy,
		Band:            band,
		ControlTargetGC: gcControlTarget,
		ControlAccuracy: controlAccuracy,
		NullPass:        math.Abs(nullAccuracy-0.5) <= band,
		ControlPass:     controlAccuracy > positiveControlMin,
	}
}

// Benchmark is the measured generation benchmark of one platform.
type Benchmark struct {
	WallSeconds    float64 `json:"wall_seconds"`
	CPUSeconds     float64 `json:"cpu_seconds"`
	CPUPercent     float64 `json:"cpu_percent_of_machine"`
	AllocatedBytes uint64  `json:"allocated_bytes"`
	BasesPerSecond float64 `json:"bases_per_second"`
}

// Realism is the classifier accuracy on real versus generated reads.
type Realism struct {
	RealReads int     `json:"real_reads"`
	Holdout   int     `json:"holdout_reads"`
	Accuracy  float64 `json:"accuracy"`
	// Band is the 4-sigma chance band around 0.5, for reference.
	Band float64 `json:"chance_band"`
}

// PlatformResult holds everything measured for one platform.
type PlatformResult struct {
	Platform    testing.SequencingPlatform `json:"platform"`
	Status      string                     `json:"status"`
	Errors      []string                   `json:"errors,omitempty"`
	Seed        uint64                     `json:"seed"`
	Reads       *testing.ReadStats         `json:"reads,omitempty"`
	Benchmark   *Benchmark                 `json:"benchmark,omitempty"`
	Calibration *Calibration               `json:"calibration,omitempty"`
	Realism     *Realism                   `json:"realism,omitempty"`
}

// failures counts the failed steps and calibration checks.
func (r *PlatformResult) failures() int {
	n := len(r.Errors)
	if r.Calibration != nil && !r.Calibration.Passed() {
		n++
	}
	return n
}

// Report is the result of a run.
type Report struct {
	StartedAt   time.Time         `json:"started_at"`
	WallSeconds float64           `json:"wall_seconds"`
	Seed        uint64            `json:"seed"`
	Quick       bool              `json:"quick"`
	RealFASTQ   string            `json:"real_fastq,omitempty"`
	RealError   string            `json:"real_fastq_error,omitempty"`
	Platforms   []*PlatformResult `json:"platforms"`
}

// Failures counts failures across platforms, plus a real FASTQ that could
// not be read.
func (r *Report) Failures() int {
	n := 0
	if r.RealError != "" {
		n++
	}
	for _, p := range r.Platforms {
		n += p.failures()
	}
	return n
}

// runAll measures every platform.
func runAll(f *testing.E2EFramework, platforms []testing.SequencingPlatform, args *CommandLineArgs) *Report {
	start := time.Now()
	report := &Report{StartedAt: start, Seed: args.seed, Quick: args.quick, RealFASTQ: args.realFASTQ}

	var real *testing.SyntheticData
	if args.realFASTQ != "" {
		var err error
		if real, err = testing.LoadFASTQ(args.realFASTQ, calibrationReads); err != nil {
			report.RealError = err.Error()
		} else if len(real.Sequences) < 2 {
			report.RealError = fmt.Sprintf("%s holds %d reads; need at least 2", args.realFASTQ, len(real.Sequences))
			real = nil
		}
	}

	for _, p := range platforms {
		if args.verbose {
			log.Printf("measuring %s", p)
		}
		report.Platforms = append(report.Platforms, runPlatform(f, p, real, args))
	}
	report.WallSeconds = time.Since(start).Seconds()
	return report
}

// runPlatform measures one platform. Dataset seeds are derived from
// args.seed: base for the read table and benchmark, base+1 and base+2 for
// the calibration sets, base+3 for the realism comparison.
func runPlatform(f *testing.E2EFramework, platform testing.SequencingPlatform, real *testing.SyntheticData, args *CommandLineArgs) *PlatformResult {
	res := &PlatformResult{Platform: platform, Seed: args.seed, Status: statusMeasured}
	config := func(seed uint64, reads int) *testing.PlatformConfig {
		c := f.CreatePlatformConfig(platform, "CEU")
		c.Seed = seed
		if reads > 0 {
			c.Reads = reads
		}
		return c
	}
	fail := func(step string, err error) {
		res.Status = statusError
		res.Errors = append(res.Errors, fmt.Sprintf("%s: %v", step, err))
	}

	data, err := f.GenerateSyntheticData(config(args.seed, 0))
	if errors.Is(err, testing.ErrPlatformNotModeled) {
		res.Status = statusNotModeled
		return res
	}
	if err != nil {
		fail("generate", err)
		return res
	}
	stats := data.Stats
	res.Reads = &stats

	if args.benchmark {
		if b, err := f.RunBenchmark(config(args.seed, 0)); err != nil {
			fail("benchmark", err)
		} else {
			res.Benchmark = &Benchmark{
				WallSeconds:    b.WallTime.Seconds(),
				CPUSeconds:     b.CPUSeconds,
				CPUPercent:     b.CPUPercent,
				AllocatedBytes: b.AllocatedBytes,
				BasesPerSecond: b.BasesPerSecond,
			}
		}
	}

	if args.authenticity {
		n := calibrationReads
		if args.quick {
			n = quickCalibrationReads
		}
		if c, err := calibrate(f, config(args.seed+1, n), config(args.seed+2, n)); err != nil {
			fail("calibration", err)
		} else {
			res.Calibration = c
		}
	}

	if real != nil {
		synth, err := f.GenerateSyntheticData(config(args.seed+3, len(real.Sequences)))
		if err != nil {
			fail("realism: generate", err)
			return res
		}
		acc, holdout, err := testing.HoldoutAccuracy(real, synth)
		if err != nil {
			fail("realism", err)
			return res
		}
		res.Realism = &Realism{RealReads: len(real.Sequences), Holdout: holdout, Accuracy: acc,
			Band: testing.ChanceBand(holdout, nullBandSigmas)}
	}
	return res
}

// calibrate measures the null comparison (two independent draws) and the
// positive control (the second draw GC-skewed to gcControlTarget).
func calibrate(f *testing.E2EFramework, ca, cb *testing.PlatformConfig) (*Calibration, error) {
	a, err := f.GenerateSyntheticData(ca)
	if err != nil {
		return nil, err
	}
	b, err := f.GenerateSyntheticData(cb)
	if err != nil {
		return nil, err
	}
	nullAcc, holdout, err := testing.HoldoutAccuracy(a, b)
	if err != nil {
		return nil, fmt.Errorf("null comparison: %w", err)
	}
	controlAcc, _, err := testing.HoldoutAccuracy(a, testing.GCSkewed(b, gcControlTarget, cb.Seed))
	if err != nil {
		return nil, fmt.Errorf("positive control: %w", err)
	}
	return gradeCalibration(ca.Reads, holdout, nullAcc, controlAcc), nil
}

// printReport writes the measured tables.
func printReport(w io.Writer, r *Report, args *CommandLineArgs) {
	rule := strings.Repeat("-", 100)
	fmt.Fprintf(w, "E2E measurements: %d platform(s), base seed %d, %.1fs\n", len(r.Platforms), r.Seed, r.WallSeconds)

	fmt.Fprintln(w, "\nREADS (computed from the generated FASTQ)")
	fmt.Fprintln(w, rule)
	fmt.Fprintf(w, "%-24s %7s %11s %9s %7s %7s %7s %7s %6s\n",
		"Platform", "Reads", "Bases", "MeanLen", "MinLen", "MaxLen", "N50", "MeanQ", "GC")
	fmt.Fprintln(w, rule)
	for _, p := range r.Platforms {
		switch {
		case p.Status == statusNotModeled:
			fmt.Fprintf(w, "%-24s NOT MODELED: the generator has no read model for this platform\n", p.Platform)
		case p.Reads == nil:
			fmt.Fprintf(w, "%-24s ERROR: %s\n", p.Platform, strings.Join(p.Errors, "; "))
		default:
			s := p.Reads
			fmt.Fprintf(w, "%-24s %7d %11d %9.1f %7d %7d %7d %7.2f %6.4f\n",
				p.Platform, s.Reads, s.Bases, s.MeanReadLength, s.MinReadLength, s.MaxReadLength, s.N50, s.MeanQuality, s.GCContent)
		}
	}

	if args.benchmark {
		fmt.Fprintln(w, "\nBENCHMARK (one generation of the read table's dataset)")
		fmt.Fprintln(w, rule)
		fmt.Fprintf(w, "%-24s %9s %8s %8s %11s %12s\n", "Platform", "Wall(s)", "CPU(s)", "CPU%", "Alloc(MiB)", "Mbases/s")
		fmt.Fprintln(w, rule)
		for _, p := range r.Platforms {
			if b := p.Benchmark; b != nil {
				fmt.Fprintf(w, "%-24s %9.3f %8.3f %8.1f %11.1f %12.2f\n", p.Platform,
					b.WallSeconds, b.CPUSeconds, b.CPUPercent, float64(b.AllocatedBytes)/(1<<20), b.BasesPerSecond/1e6)
			}
		}
	}

	if args.authenticity {
		fmt.Fprintln(w, "\nCLASSIFIER CALIBRATION (4-mer Naive Bayes, half of each set held out)")
		fmt.Fprintf(w, "null: two independent draws, pass within 0.5 ± %.0fσ; control: second draw skewed to GC %.2f, pass above %.2f\n",
			nullBandSigmas, gcControlTarget, positiveControlMin)
		fmt.Fprintln(w, rule)
		fmt.Fprintf(w, "%-24s %9s %8s %18s %5s %9s %5s\n", "Platform", "Reads/set", "Holdout", "Null accuracy", "", "Control", "")
		fmt.Fprintln(w, rule)
		for _, p := range r.Platforms {
			if c := p.Calibration; c != nil {
				fmt.Fprintf(w, "%-24s %9d %8d %8.4f ± %7.4f %5s %9.4f %5s\n", p.Platform,
					c.ReadsPerSet, c.Holdout, c.NullAccuracy, c.Band, passFail(c.NullPass), c.ControlAccuracy, passFail(c.ControlPass))
			}
		}
	}

	if r.RealFASTQ != "" {
		fmt.Fprintf(w, "\nREAL VERSUS GENERATED (%s)\n", r.RealFASTQ)
		fmt.Fprintln(w, rule)
		if r.RealError != "" {
			fmt.Fprintf(w, "ERROR: %s\n", r.RealError)
		}
		for _, p := range r.Platforms {
			if x := p.Realism; x != nil {
				fmt.Fprintf(w, "%-24s %d real reads, %d held out: accuracy %.4f (chance band 0.5 ± %.4f)\n",
					p.Platform, x.RealReads, x.Holdout, x.Accuracy, x.Band)
			}
		}
	} else {
		fmt.Fprintln(w, "\nNo real reads supplied (-real-fastq or $PHANTOM_REAL_FASTQ), so no real-versus-generated comparison was made.")
	}

	var errs []string
	for _, p := range r.Platforms {
		if p.Reads != nil {
			for _, e := range p.Errors {
				errs = append(errs, fmt.Sprintf("%s: %s", p.Platform, e))
			}
		}
	}
	if len(errs) > 0 {
		fmt.Fprintln(w, "\nERRORS")
		for _, e := range errs {
			fmt.Fprintln(w, "  "+e)
		}
	}
	fmt.Fprintf(w, "\nFailures: %d\n", r.Failures())
}

func passFail(ok bool) string {
	if ok {
		return "PASS"
	}
	return "FAIL"
}

// writeReports writes e2e_results.json (the full report) and
// e2e_summary.csv (one row per platform) to dir and returns their paths.
func writeReports(r *Report, dir string) ([]string, error) {
	if err := os.MkdirAll(dir, 0750); err != nil {
		return nil, err
	}
	jsonPath := filepath.Join(dir, "e2e_results.json")
	b, err := json.MarshalIndent(r, "", "  ")
	if err != nil {
		return nil, err
	}
	if err := os.WriteFile(jsonPath, append(b, '\n'), 0600); err != nil {
		return nil, err
	}

	csvPath := filepath.Join(dir, "e2e_summary.csv")
	fh, err := os.Create(csvPath) // #nosec G304 -- caller-chosen output directory
	if err != nil {
		return nil, err
	}
	cw := csv.NewWriter(fh)
	_ = cw.Write([]string{"platform", "status", "seed", "reads", "bases", "mean_read_length", "min_read_length",
		"max_read_length", "n50", "mean_quality", "gc_content", "wall_seconds", "cpu_seconds", "allocated_bytes",
		"bases_per_second", "null_accuracy", "null_band", "control_accuracy", "calibration_pass", "real_accuracy", "errors"})
	ff := func(v float64) string { return strconv.FormatFloat(v, 'f', -1, 64) }
	for _, p := range r.Platforms {
		row := make([]string, 21)
		row[0], row[1], row[2] = string(p.Platform), p.Status, strconv.FormatUint(p.Seed, 10)
		if s := p.Reads; s != nil {
			row[3], row[4] = strconv.Itoa(s.Reads), strconv.FormatInt(s.Bases, 10)
			row[5], row[6], row[7], row[8] = ff(s.MeanReadLength), strconv.Itoa(s.MinReadLength), strconv.Itoa(s.MaxReadLength), strconv.Itoa(s.N50)
			row[9], row[10] = ff(s.MeanQuality), ff(s.GCContent)
		}
		if b := p.Benchmark; b != nil {
			row[11], row[12], row[13], row[14] = ff(b.WallSeconds), ff(b.CPUSeconds), strconv.FormatUint(b.AllocatedBytes, 10), ff(b.BasesPerSecond)
		}
		if c := p.Calibration; c != nil {
			row[15], row[16], row[17], row[18] = ff(c.NullAccuracy), ff(c.Band), ff(c.ControlAccuracy), strconv.FormatBool(c.Passed())
		}
		if x := p.Realism; x != nil {
			row[19] = ff(x.Accuracy)
		}
		row[20] = strings.Join(p.Errors, "; ")
		_ = cw.Write(row)
	}
	cw.Flush()
	if err := cw.Error(); err != nil {
		fh.Close()
		return nil, err
	}
	if err := fh.Close(); err != nil {
		return nil, err
	}
	return []string{jsonPath, csvPath}, nil
}
