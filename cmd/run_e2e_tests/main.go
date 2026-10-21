// Command run_e2e_tests generates reads for each platform profile of the
// internal/testing E2E framework with the real generator and reports what
// it measures: read statistics, generation benchmarks, a classifier
// calibration per feature set and, when real reads are supplied, real-
// versus-real and real-versus-generated classifier accuracies. Bases are
// uniform random, or sampled from a reference FASTA with -reference. It
// prints no scores or verdicts; it exits non-zero only when a step fails or
// a calibration check fails.
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
	reference    string
	region       string
	reads        int
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
	if args.reference == "" {
		args.reference = os.Getenv("PHANTOM_REFERENCE")
	}
	if args.region == "" {
		args.region = os.Getenv("PHANTOM_REGION")
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
	flag.StringVar(&args.reference, "reference", "",
		"Uncompressed reference FASTA to sample generated reads from (default $PHANTOM_REFERENCE; empty: uniform random bases)")
	flag.StringVar(&args.region, "region", "",
		"Reference region to sample from, contig or contig:start-end, 1-based (default $PHANTOM_REGION; needs -reference)")
	flag.IntVar(&args.reads, "reads", 0,
		fmt.Sprintf("Reads per set for the calibration and the real comparison (default %d, %d with -quick)", calibrationReads, quickCalibrationReads))
	flag.BoolVar(&args.verbose, "verbose", false, "Enable verbose logging")
	flag.BoolVar(&args.quick, "quick", false, "Use fewer reads for the calibration")
	flag.BoolVar(&args.benchmark, "benchmark", true, "Measure generation time, CPU and allocation")
	flag.BoolVar(&args.report, "report", true, "Write e2e_results.json and e2e_summary.csv to -output")
	flag.BoolVar(&args.authenticity, "authenticity", true,
		"Run the classifier calibration (null comparison and positive controls) for each feature set")
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

BASES:
  Uniform random by default. With -reference every generated read is a
  fragment sampled from the FASTA (or from -region of it). The mode is
  printed and recorded in the JSON report.

FEATURE SETS (one feature vector per read; classifiers trained on half of
each set and scored on the other half):
  sequence  256 4-mer counts, multinomial naive Bayes.
  quality   quality histogram, mean quality of 10 cycle blocks, low-quality
            (Q<20) runs, quality changes, N fraction, mean quality per
            called base; logistic regression.
  combined  both, logistic regression.

WHAT IS REPORTED:
  Reads    read count, bases, mean/min/max read length, N50, mean Phred
           quality and GC content, computed from the generated FASTQ.
  Bench    (-benchmark) wall time, process CPU seconds and share of the
           machine, heap bytes allocated, and bases per second.
  Calib    (-authenticity) per feature set:
             null     two independent generator draws; passes when accuracy
                      is within 0.5 +/- %.0f sigma (binomial).
             GC       the second draw GC-skewed to %.2f.
             quality  the second draw generated with the other quality model
                      (Illumina binned <-> linear, ONT tier <-> none).
           A set passes its controls when the ones it can see exceed %.2f:
           sequence GC, quality the quality model, combined both.
           %d reads per set (%d with -quick, or -reads).
  Real     (-real-fastq or $PHANTOM_REAL_FASTQ) per feature set:
             real vs real       alternate reads of the real file against
                                each other; must lie within the null band.
             real vs generated  real reads against as many generated ones
                                (fixed-length profiles take the real read
                                length), with the features that separate
                                them most. Reported, not graded.

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
	if args.reads < 0 || (args.reads > 0 && args.reads < 4) {
		return fmt.Errorf("-reads must be at least 4, got %d", args.reads)
	}
	if args.realFASTQ != "" {
		if _, err := os.Stat(args.realFASTQ); err != nil {
			return fmt.Errorf("real FASTQ: %w", err)
		}
	}
	if args.region != "" && args.reference == "" {
		return fmt.Errorf("-region %q needs -reference", args.region)
	}
	if args.reference != "" {
		if _, err := os.Stat(args.reference); err != nil {
			return fmt.Errorf("reference: %w", err)
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

// readsPerSet returns the reads per calibration set and real comparison.
func (a *CommandLineArgs) readsPerSet() int {
	switch {
	case a.reads > 0:
		return a.reads
	case a.quick:
		return quickCalibrationReads
	}
	return calibrationReads
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
	// to (uniform bases are about 0.5, human reference fragments about 0.41).
	gcControlTarget = 0.7
	// positiveControlMin is the accuracy a positive control must exceed.
	positiveControlMin = 0.9
	// topFeatures is the number of separating features kept per real
	// comparison (the printed table shows topFeaturesPrinted).
	topFeatures        = 10
	topFeaturesPrinted = 5
)

// Status values for PlatformResult.
const (
	statusMeasured   = "measured"
	statusNotModeled = "not_modeled"
	statusError      = "error"
)

// Bases describes where generated bases came from.
type Bases struct {
	// Source is "uniform" or "reference".
	Source    string `json:"source"`
	Reference string `json:"reference,omitempty"`
	Region    string `json:"region,omitempty"`
}

func (b Bases) String() string {
	switch {
	case b.Source != testing.BasesReference:
		return "uniform random bases"
	case b.Region != "":
		return fmt.Sprintf("reference %s region %s", b.Reference, b.Region)
	}
	return "reference " + b.Reference
}

// SetCalibration is the calibration of one feature set.
type SetCalibration struct {
	FeatureSet             testing.FeatureSet `json:"feature_set"`
	Classifier             string             `json:"classifier"`
	Holdout                int                `json:"holdout_reads"`
	NullAccuracy           float64            `json:"null_accuracy"`
	Band                   float64            `json:"null_band"`
	GCControlAccuracy      float64            `json:"gc_control_accuracy"`
	QualityControlAccuracy float64            `json:"quality_control_accuracy"`
	NullPass               bool               `json:"null_pass"`
	ControlPass            bool               `json:"control_pass"`
}

// Calibration holds the classifier calibration measured on one platform.
type Calibration struct {
	ReadsPerSet     int     `json:"reads_per_set"`
	Bases           Bases   `json:"bases"`
	ControlTargetGC float64 `json:"control_target_gc"`
	// QualityControl names the alternative quality model of the quality
	// control set.
	QualityControl string            `json:"quality_control"`
	Sets           []*SetCalibration `json:"feature_sets"`
}

// Passed reports whether every check of every feature set passed.
func (c *Calibration) Passed() bool {
	for _, s := range c.Sets {
		if !s.NullPass || !s.ControlPass {
			return false
		}
	}
	return len(c.Sets) > 0
}

// Set returns the calibration of set, or nil.
func (c *Calibration) Set(set testing.FeatureSet) *SetCalibration {
	for _, s := range c.Sets {
		if s.FeatureSet == set {
			return s
		}
	}
	return nil
}

// gradeSet grades one feature set's null and control accuracies measured
// on holdout reads. The sequence set must detect the GC skew, the quality
// set the quality model, the combined set both.
func gradeSet(set testing.FeatureSet, holdout int, null, gc, qual float64) *SetCalibration {
	band := testing.ChanceBand(holdout, nullBandSigmas)
	var control bool
	switch set {
	case testing.FeatureSequence:
		control = gc > positiveControlMin
	case testing.FeatureQuality:
		control = qual > positiveControlMin
	default:
		control = gc > positiveControlMin && qual > positiveControlMin
	}
	return &SetCalibration{
		FeatureSet: set, Classifier: testing.ClassifierName(set), Holdout: holdout,
		NullAccuracy: null, Band: band, GCControlAccuracy: gc, QualityControlAccuracy: qual,
		NullPass: math.Abs(null-0.5) <= band, ControlPass: control,
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

// Comparison is one feature set's holdout accuracy between two read sets.
type Comparison struct {
	*testing.Holdout
	// Band is the 4-sigma chance band around 0.5.
	Band float64 `json:"chance_band"`
}

// Realism holds the comparisons of real reads (class 0) with each other
// and with generated reads (class 1).
type Realism struct {
	RealReads            int     `json:"real_reads"`
	GeneratedReads       int     `json:"generated_reads"`
	GeneratedReadLength  int     `json:"generated_read_length,omitempty"`
	Bases                Bases   `json:"bases"`
	RealMeanQuality      float64 `json:"real_mean_quality"`
	GeneratedMeanQuality float64 `json:"generated_mean_quality"`
	RealGC               float64 `json:"real_gc"`
	GeneratedGC          float64 `json:"generated_gc"`
	// RealVsReal compares alternate reads of the real file: the null.
	RealVsReal []*Comparison `json:"real_vs_real"`
	// RealVsGenerated compares the real reads with generated ones.
	RealVsGenerated []*Comparison `json:"real_vs_generated"`
}

// nullFailures counts real-versus-real comparisons outside the null band.
func (r *Realism) nullFailures() int {
	n := 0
	for _, c := range r.RealVsReal {
		if math.Abs(c.Accuracy-0.5) > c.Band {
			n++
		}
	}
	return n
}

// find returns set's comparison in cs, or nil.
func find(cs []*Comparison, set testing.FeatureSet) *Comparison {
	for _, c := range cs {
		if c.FeatureSet == set {
			return c
		}
	}
	return nil
}

// PlatformResult holds everything measured for one platform.
type PlatformResult struct {
	Platform    testing.SequencingPlatform `json:"platform"`
	Status      string                     `json:"status"`
	Errors      []string                   `json:"errors,omitempty"`
	Seed        uint64                     `json:"seed"`
	Bases       Bases                      `json:"bases"`
	Reads       *testing.ReadStats         `json:"reads,omitempty"`
	Benchmark   *Benchmark                 `json:"benchmark,omitempty"`
	Calibration *Calibration               `json:"calibration,omitempty"`
	Realism     *Realism                   `json:"realism,omitempty"`
}

// failures counts the failed steps and calibration checks, including a
// real-versus-real comparison outside the null band.
func (r *PlatformResult) failures() int {
	n := len(r.Errors)
	if r.Calibration != nil && !r.Calibration.Passed() {
		n++
	}
	if r.Realism != nil && r.Realism.nullFailures() > 0 {
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
	ReadsPerSet int               `json:"reads_per_set"`
	Bases       Bases             `json:"bases"`
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

// basesOf returns the base source of args.
func basesOf(args *CommandLineArgs) Bases {
	if args.reference == "" {
		return Bases{Source: testing.BasesUniform}
	}
	return Bases{Source: testing.BasesReference, Reference: args.reference, Region: args.region}
}

// runAll measures every platform.
func runAll(f *testing.E2EFramework, platforms []testing.SequencingPlatform, args *CommandLineArgs) *Report {
	start := time.Now()
	report := &Report{StartedAt: start, Seed: args.seed, Quick: args.quick, ReadsPerSet: args.readsPerSet(),
		Bases: basesOf(args), RealFASTQ: args.realFASTQ}

	var real *testing.SyntheticData
	if args.realFASTQ != "" {
		var err error
		if real, err = testing.LoadFASTQ(args.realFASTQ, args.readsPerSet()); err != nil {
			report.RealError = err.Error()
		} else if len(real.Sequences) < 4 {
			report.RealError = fmt.Sprintf("%s holds %d reads; need at least 4", args.realFASTQ, len(real.Sequences))
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

// commonLength returns the length shared by every read, or 0.
func commonLength(d *testing.SyntheticData) int {
	if d.Stats.MinReadLength == d.Stats.MaxReadLength {
		return d.Stats.MinReadLength
	}
	return 0
}

// runPlatform measures one platform. Dataset seeds are derived from
// args.seed: base for the read table and benchmark, base+1 and base+2 for
// the calibration sets (base+2 also for the quality control set), base+3
// for the real comparison.
func runPlatform(f *testing.E2EFramework, platform testing.SequencingPlatform, real *testing.SyntheticData, args *CommandLineArgs) *PlatformResult {
	res := &PlatformResult{Platform: platform, Seed: args.seed, Status: statusMeasured, Bases: basesOf(args)}
	config := func(seed uint64, reads int) *testing.PlatformConfig {
		c := f.CreatePlatformConfig(platform, "CEU")
		c.Seed = seed
		c.ReferencePath, c.ReferenceRegion = args.reference, args.region
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
		n := args.readsPerSet()
		if c, err := calibrate(f, config(args.seed+1, n), config(args.seed+2, n)); err != nil {
			fail("calibration", err)
		} else {
			c.Bases = res.Bases
			res.Calibration = c
		}
	}

	if real != nil {
		c := config(args.seed+3, len(real.Sequences))
		if l := commonLength(real); l > 0 && c.ReadLengthDistribution.Dist == "" {
			c.ReadLength = l
		}
		synth, err := f.GenerateSyntheticData(c)
		if err != nil {
			fail("realism: generate", err)
			return res
		}
		x, err := compareReal(real, synth)
		if err != nil {
			fail("realism", err)
			return res
		}
		x.Bases = res.Bases
		if c.ReadLengthDistribution.Dist == "" {
			x.GeneratedReadLength = c.ReadLength
		}
		res.Realism = x
	}
	return res
}

// compareReal runs the real-versus-real null and the real-versus-generated
// comparison for every feature set.
func compareReal(real, synth *testing.SyntheticData) (*Realism, error) {
	x := &Realism{RealReads: len(real.Sequences), GeneratedReads: len(synth.Sequences),
		RealMeanQuality: real.Stats.MeanQuality, GeneratedMeanQuality: synth.Stats.MeanQuality,
		RealGC: real.Stats.GCContent, GeneratedGC: synth.Stats.GCContent}
	even, odd := testing.SplitAlternate(real)
	for _, set := range testing.AllFeatureSets {
		null, err := testing.HoldoutAccuracyFor(set, even, odd, 0)
		if err != nil {
			return nil, fmt.Errorf("real vs real (%s): %w", set, err)
		}
		x.RealVsReal = append(x.RealVsReal, &Comparison{Holdout: null, Band: testing.ChanceBand(null.Holdout, nullBandSigmas)})
		gen, err := testing.HoldoutAccuracyFor(set, real, synth, topFeatures)
		if err != nil {
			return nil, fmt.Errorf("real vs generated (%s): %w", set, err)
		}
		x.RealVsGenerated = append(x.RealVsGenerated, &Comparison{Holdout: gen, Band: testing.ChanceBand(gen.Holdout, nullBandSigmas)})
	}
	return x, nil
}

// calibrate measures, for each feature set, the null comparison (two
// independent draws), the GC control (the second draw GC-skewed to
// gcControlTarget) and the quality control (the second draw's seed with
// the alternative quality model).
func calibrate(f *testing.E2EFramework, ca, cb *testing.PlatformConfig) (*Calibration, error) {
	a, err := f.GenerateSyntheticData(ca)
	if err != nil {
		return nil, err
	}
	b, err := f.GenerateSyntheticData(cb)
	if err != nil {
		return nil, err
	}
	alt := cb.AlternativeQualityModel()
	bq, err := f.GenerateSyntheticData(alt)
	if err != nil {
		return nil, fmt.Errorf("quality control: %w", err)
	}
	bgc := testing.GCSkewed(b, gcControlTarget, cb.Seed)
	cal := &Calibration{ReadsPerSet: ca.Reads, ControlTargetGC: gcControlTarget, QualityControl: describeQualityModel(alt)}
	for _, set := range testing.AllFeatureSets {
		null, err := testing.HoldoutAccuracyFor(set, a, b, 0)
		if err != nil {
			return nil, fmt.Errorf("null comparison (%s): %w", set, err)
		}
		gc, err := testing.HoldoutAccuracyFor(set, a, bgc, 0)
		if err != nil {
			return nil, fmt.Errorf("GC control (%s): %w", set, err)
		}
		q, err := testing.HoldoutAccuracyFor(set, a, bq, 0)
		if err != nil {
			return nil, fmt.Errorf("quality control (%s): %w", set, err)
		}
		cal.Sets = append(cal.Sets, gradeSet(set, null.Holdout, null.Accuracy, gc.Accuracy, q.Accuracy))
	}
	return cal, nil
}

// describeQualityModel names the quality model of c.
func describeQualityModel(c *testing.PlatformConfig) string {
	if c.EnginePlatform == "ont" {
		if c.ErrorTier == "" {
			return "ont linear (no basecaller tier)"
		}
		return "ont " + c.ErrorTier + " tier"
	}
	if c.IlluminaQualityModel == "" {
		return "illumina binned"
	}
	return "illumina " + c.IlluminaQualityModel
}

// printReport writes the measured tables.
func printReport(w io.Writer, r *Report, args *CommandLineArgs) {
	rule := strings.Repeat("-", 100)
	fmt.Fprintf(w, "E2E measurements: %d platform(s), base seed %d, %.1fs\n", len(r.Platforms), r.Seed, r.WallSeconds)
	fmt.Fprintf(w, "Generated bases: %s\n", r.Bases)

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
		fmt.Fprintln(w, "\nCLASSIFIER CALIBRATION (half of each set held out)")
		fmt.Fprintf(w, "null: two independent draws, pass within 0.5 ± %.0fσ; GC: second draw skewed to GC %.2f; quality: second draw with the other quality model\n",
			nullBandSigmas, gcControlTarget)
		fmt.Fprintf(w, "controls pass above %.2f: sequence on GC, quality on quality, combined on both\n", positiveControlMin)
		fmt.Fprintln(w, rule)
		fmt.Fprintf(w, "%-24s %-9s %9s %8s %18s %5s %8s %8s %5s\n", "Platform", "Features", "Reads/set", "Holdout", "Null accuracy", "", "GC", "Quality", "")
		fmt.Fprintln(w, rule)
		for _, p := range r.Platforms {
			c := p.Calibration
			if c == nil {
				continue
			}
			for _, s := range c.Sets {
				fmt.Fprintf(w, "%-24s %-9s %9d %8d %8.4f ± %7.4f %5s %8.4f %8.4f %5s\n", p.Platform, s.FeatureSet,
					c.ReadsPerSet, s.Holdout, s.NullAccuracy, s.Band, passFail(s.NullPass),
					s.GCControlAccuracy, s.QualityControlAccuracy, passFail(s.ControlPass))
			}
		}
	}

	if r.RealFASTQ != "" {
		fmt.Fprintf(w, "\nREAL VERSUS GENERATED (%s; generated from %s)\n", r.RealFASTQ, r.Bases)
		fmt.Fprintln(w, rule)
		if r.RealError != "" {
			fmt.Fprintf(w, "ERROR: %s\n", r.RealError)
		}
		for _, p := range r.Platforms {
			x := p.Realism
			if x == nil {
				continue
			}
			fmt.Fprintf(w, "%s: %d real reads (mean Q %.2f, GC %.4f), %d generated (mean Q %.2f, GC %.4f",
				p.Platform, x.RealReads, x.RealMeanQuality, x.RealGC, x.GeneratedReads, x.GeneratedMeanQuality, x.GeneratedGC)
			if x.GeneratedReadLength > 0 {
				fmt.Fprintf(w, ", %d bp", x.GeneratedReadLength)
			}
			fmt.Fprintln(w, ")")
			fmt.Fprintf(w, "  %-9s %-24s %8s %18s %5s %8s %18s\n", "Features", "Classifier", "Holdout", "Real vs real", "", "Holdout", "Real vs generated")
			for _, set := range testing.AllFeatureSets {
				n, g := find(x.RealVsReal, set), find(x.RealVsGenerated, set)
				if n == nil || g == nil {
					continue
				}
				fmt.Fprintf(w, "  %-9s %-24s %8d %8.4f ± %7.4f %5s %8d %8.4f ± %7.4f\n", set, g.Classifier,
					n.Holdout.Holdout, n.Accuracy, n.Band, passFail(math.Abs(n.Accuracy-0.5) <= n.Band),
					g.Holdout.Holdout, g.Accuracy, g.Band)
			}
			fmt.Fprintln(w, "  Features separating real (a) from generated (b), by effect size (b−a)/SD on the training reads:")
			for _, set := range []testing.FeatureSet{testing.FeatureSequence, testing.FeatureQuality} {
				g := find(x.RealVsGenerated, set)
				if g == nil {
					continue
				}
				for i, d := range g.Top {
					if i == topFeaturesPrinted {
						break
					}
					fmt.Fprintf(w, "    %-9s %-28s real %10.5f  generated %10.5f  d %+7.3f\n", set, d.Name, d.MeanA, d.MeanB, d.EffectSize)
				}
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

// csvHeader returns the e2e_summary.csv columns.
func csvHeader() []string {
	h := []string{"platform", "status", "seed", "bases_source", "reference", "region", "reads", "bases", "mean_read_length",
		"min_read_length", "max_read_length", "n50", "mean_quality", "gc_content", "wall_seconds", "cpu_seconds",
		"allocated_bytes", "bases_per_second", "calibration_pass"}
	for _, set := range testing.AllFeatureSets {
		s := string(set)
		h = append(h, "null_accuracy_"+s, "null_band_"+s, "gc_control_accuracy_"+s, "quality_control_accuracy_"+s,
			"real_vs_real_accuracy_"+s, "real_vs_generated_accuracy_"+s, "real_band_"+s)
	}
	return append(h, "errors")
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
	_ = cw.Write(csvHeader())
	ff := func(v float64) string { return strconv.FormatFloat(v, 'f', -1, 64) }
	for _, p := range r.Platforms {
		row := []string{string(p.Platform), p.Status, strconv.FormatUint(p.Seed, 10), p.Bases.Source, p.Bases.Reference, p.Bases.Region}
		cells := make([]string, 13)
		if s := p.Reads; s != nil {
			cells[0], cells[1] = strconv.Itoa(s.Reads), strconv.FormatInt(s.Bases, 10)
			cells[2], cells[3], cells[4], cells[5] = ff(s.MeanReadLength), strconv.Itoa(s.MinReadLength), strconv.Itoa(s.MaxReadLength), strconv.Itoa(s.N50)
			cells[6], cells[7] = ff(s.MeanQuality), ff(s.GCContent)
		}
		if b := p.Benchmark; b != nil {
			cells[8], cells[9], cells[10], cells[11] = ff(b.WallSeconds), ff(b.CPUSeconds), strconv.FormatUint(b.AllocatedBytes, 10), ff(b.BasesPerSecond)
		}
		if c := p.Calibration; c != nil {
			cells[12] = strconv.FormatBool(c.Passed())
		}
		row = append(row, cells...)
		for _, set := range testing.AllFeatureSets {
			sc := make([]string, 7)
			if c := p.Calibration; c != nil {
				if s := c.Set(set); s != nil {
					sc[0], sc[1], sc[2], sc[3] = ff(s.NullAccuracy), ff(s.Band), ff(s.GCControlAccuracy), ff(s.QualityControlAccuracy)
				}
			}
			if x := p.Realism; x != nil {
				if n := find(x.RealVsReal, set); n != nil {
					sc[4] = ff(n.Accuracy)
				}
				if g := find(x.RealVsGenerated, set); g != nil {
					sc[5], sc[6] = ff(g.Accuracy), ff(g.Band)
				}
			}
			row = append(row, sc...)
		}
		row = append(row, strings.Join(p.Errors, "; "))
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
