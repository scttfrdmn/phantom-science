// Command run_e2e_tests generates reads for each platform profile of the
// internal/testing E2E framework with the real generator and reports what
// it measures: read statistics, generation benchmarks, a classifier
// calibration per feature set and, when real reads are supplied, real-
// versus-real and real-versus-generated classifier accuracies, for read 1,
// read 2 and read pairs when both files of a pair are supplied, and the
// accuracies between two real runs (-real-fastq-b), the target a model can
// reach. Bases are uniform random, or sampled from a reference FASTA with
// -reference. It prints no scores or verdicts; it exits non-zero only when
// a step fails or a calibration check fails.
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

	"github.com/scttfrdmn/phantom-science/internal/genomics/illumina"
	"github.com/scttfrdmn/phantom-science/internal/testing"
)

// CommandLineArgs holds command line arguments.
type CommandLineArgs struct {
	platforms    string
	outputDir    string
	seed         uint64
	realFASTQ    string
	realFASTQ2   string
	realFASTQB   string
	realFASTQB2  string
	reference    string
	region       string
	reads        int
	verbose      bool
	quick        bool
	benchmark    bool
	report       bool
	authenticity bool
	paired       bool
	// realOnly skips the real-versus-generated comparison: only the
	// real-versus-real nulls and the cross-run comparison are run.
	realOnly bool
	// writeGenerated, when set, is the directory the generated reads of the
	// real comparison are written to as FASTQ.
	writeGenerated string
	// illuminaProfile is the Illumina run profile of the generated reads
	// (illumina.profile; empty: the default).
	illuminaProfile string
	help            bool
}

func main() {
	args := parseCommandLineArgs()
	if args.help {
		printUsage()
		return
	}
	for _, e := range []struct {
		v   *string
		env string
	}{
		{&args.realFASTQ, "PHANTOM_REAL_FASTQ"}, {&args.realFASTQ2, "PHANTOM_REAL_FASTQ2"},
		{&args.realFASTQB, "PHANTOM_REAL_FASTQ_B"}, {&args.realFASTQB2, "PHANTOM_REAL_FASTQ2_B"},
	} {
		if *e.v == "" {
			*e.v = os.Getenv(e.env)
		}
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
		"FASTQ of real reads to compare generated reads with (default $PHANTOM_REAL_FASTQ); read 1 with -real-fastq2")
	flag.StringVar(&args.realFASTQ2, "real-fastq2", "",
		"FASTQ of the real reads 2 (mates of -real-fastq, same order): compare read 2 and pairs too, with generated pairs (default $PHANTOM_REAL_FASTQ2)")
	flag.StringVar(&args.realFASTQB, "real-fastq-b", "",
		"FASTQ of a second real run to compare -real-fastq with (default $PHANTOM_REAL_FASTQ_B)")
	flag.StringVar(&args.realFASTQB2, "real-fastq2-b", "",
		"FASTQ of the second run's reads 2; required with -real-fastq-b when -real-fastq2 is set (default $PHANTOM_REAL_FASTQ2_B)")
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
	flag.BoolVar(&args.paired, "paired", false,
		"Generate read pairs for the calibration and add the pair feature set (implied by -real-fastq2)")
	flag.BoolVar(&args.realOnly, "real-only", false,
		"With real reads, skip the comparison with generated reads: only real versus real and the cross-run comparison")
	flag.StringVar(&args.writeGenerated, "write-generated", "",
		"Directory to write the generated reads of the real comparison to (<platform>_1.fastq, and _2.fastq for pairs)")
	flag.StringVar(&args.illuminaProfile, "illumina-profile", "",
		"Illumina run profile of the generated reads: "+strings.Join(illumina.ProfileNames(), ", ")+" (default "+illumina.DefaultProfile+")")
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

FEATURE SETS (one feature vector per read or pair; classifiers trained on
half of each set and scored on the other half):
  sequence  256 4-mer counts, multinomial naive Bayes.
  quality   quality histogram, mean quality of 10 cycle blocks, low-quality
            (Q<20) runs, quality changes, N fraction, mean quality per
            called base; logistic regression.
  combined  both, logistic regression.
  pair      (paired reads) the quality features of read 1 and of read 2,
            the mates' mean-quality difference, the product and minimum of
            their low-quality fractions and the correlation of their block
            means; logistic regression. No bases, no insert size.

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
           sequence GC, quality and pair the quality model, combined both.
           %d reads per set (%d with -quick, or -reads). With -paired or
           -real-fastq2 the sets are read pairs and the pair set is added.
  Real     (-real-fastq or $PHANTOM_REAL_FASTQ) per feature set:
             real vs real       alternate reads of the real file against
                                each other; must lie within the null band.
             real vs generated  real reads against as many generated ones
                                (fixed-length profiles take the real read
                                length), with the features that separate
                                them most. Reported, not graded.
           With -real-fastq2 the generated reads are pairs (default insert
           sizes) and read 2 and the pairs are compared the same way.
  Cross    (-real-fastq-b) the first real run against the second on every
           feature set (and read 2 and pairs when both are paired), reads
           cut to the shorter run's length: how far apart two real runs
           are. Reported, not graded.

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
	for _, f := range []struct{ name, path string }{
		{"real FASTQ", args.realFASTQ}, {"-real-fastq2", args.realFASTQ2},
		{"-real-fastq-b", args.realFASTQB}, {"-real-fastq2-b", args.realFASTQB2},
	} {
		if f.path != "" {
			if _, err := os.Stat(f.path); err != nil {
				return fmt.Errorf("%s: %w", f.name, err)
			}
		}
	}
	switch {
	case args.realFASTQ2 != "" && args.realFASTQ == "":
		return fmt.Errorf("-real-fastq2 needs -real-fastq (read 1)")
	case args.realFASTQB != "" && args.realFASTQ == "":
		return fmt.Errorf("-real-fastq-b needs -real-fastq")
	case args.realFASTQB2 != "" && args.realFASTQB == "":
		return fmt.Errorf("-real-fastq2-b needs -real-fastq-b (read 1)")
	case args.realFASTQB != "" && (args.realFASTQ2 == "") != (args.realFASTQB2 == ""):
		return fmt.Errorf("both real runs must be paired (-real-fastq2 and -real-fastq2-b) or neither")
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

// pairedRun reports whether the calibration and the real comparison use
// read pairs.
func (a *CommandLineArgs) pairedRun() bool { return a.paired || a.realFASTQ2 != "" }

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
	Paired          bool    `json:"paired,omitempty"`
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
// and pair sets the quality model, the combined set both.
func gradeSet(set testing.FeatureSet, holdout int, null, gc, qual float64) *SetCalibration {
	band := testing.ChanceBand(holdout, nullBandSigmas)
	var control bool
	switch set {
	case testing.FeatureSequence:
		control = gc > positiveControlMin
	case testing.FeatureQuality, testing.FeaturePair:
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

// MateRealism holds the comparisons of the real reads 2 (or pairs) with
// each other and with the generated reads 2 (pairs).
type MateRealism struct {
	RealMeanQuality      float64 `json:"real_mean_quality,omitempty"`
	GeneratedMeanQuality float64 `json:"generated_mean_quality,omitempty"`
	RealGC               float64 `json:"real_gc,omitempty"`
	GeneratedGC          float64 `json:"generated_gc,omitempty"`
	// RealVsReal compares alternate reads (pairs) of the real files.
	RealVsReal []*Comparison `json:"real_vs_real"`
	// RealVsGenerated compares the real reads (pairs) with generated ones.
	RealVsGenerated []*Comparison `json:"real_vs_generated,omitempty"`
}

// Realism holds the comparisons of real reads (class 0) with each other
// and with generated reads (class 1). For paired real reads the top-level
// comparisons are of read 1; Read2 and Pair hold the others.
type Realism struct {
	RealReads            int     `json:"real_reads"`
	GeneratedReads       int     `json:"generated_reads"`
	GeneratedReadLength  int     `json:"generated_read_length,omitempty"`
	Paired               bool    `json:"paired,omitempty"`
	Bases                Bases   `json:"bases"`
	RealMeanQuality      float64 `json:"real_mean_quality"`
	GeneratedMeanQuality float64 `json:"generated_mean_quality"`
	RealGC               float64 `json:"real_gc"`
	GeneratedGC          float64 `json:"generated_gc"`
	// RealVsReal compares alternate reads of the real file: the null.
	RealVsReal []*Comparison `json:"real_vs_real"`
	// RealVsGenerated compares the real reads with generated ones (empty
	// with -real-only).
	RealVsGenerated []*Comparison `json:"real_vs_generated"`
	// Read2 compares the reads 2 on the single-read feature sets, Pair the
	// read pairs on the pair feature set (paired real reads only).
	Read2 *MateRealism `json:"read2,omitempty"`
	Pair  *MateRealism `json:"pair,omitempty"`
}

// nullFailures counts real-versus-real comparisons outside the null band,
// for read 1, read 2 and pairs.
func (r *Realism) nullFailures() int {
	nulls := [][]*Comparison{r.RealVsReal}
	for _, m := range []*MateRealism{r.Read2, r.Pair} {
		if m != nil {
			nulls = append(nulls, m.RealVsReal)
		}
	}
	n := 0
	for _, cs := range nulls {
		for _, c := range cs {
			if math.Abs(c.Accuracy-0.5) > c.Band {
				n++
			}
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

// CrossRun compares two real runs (run A class 0, run B class 1) on every
// feature set: how far apart real runs are, the realistic target for
// real versus generated. Not graded.
type CrossRun struct {
	FASTQA string `json:"fastq_a"`
	FASTQB string `json:"fastq_b"`
	ReadsA int    `json:"reads_a"`
	ReadsB int    `json:"reads_b"`
	// ReadLength is the length both runs' reads were cut to when their
	// lengths differ (0: not cut).
	ReadLength   int     `json:"read_length,omitempty"`
	MeanQualityA float64 `json:"mean_quality_a"`
	MeanQualityB float64 `json:"mean_quality_b"`
	GCA          float64 `json:"gc_a"`
	GCB          float64 `json:"gc_b"`
	// Read1, Read2 (paired runs) and Pair (paired runs) hold the
	// comparisons of the reads 1, the reads 2 and the pairs.
	Read1 []*Comparison `json:"read1"`
	Read2 []*Comparison `json:"read2,omitempty"`
	Pair  []*Comparison `json:"pair,omitempty"`
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
	StartedAt   time.Time `json:"started_at"`
	WallSeconds float64   `json:"wall_seconds"`
	Seed        uint64    `json:"seed"`
	Quick       bool      `json:"quick"`
	ReadsPerSet int       `json:"reads_per_set"`
	Paired      bool      `json:"paired,omitempty"`
	Bases       Bases     `json:"bases"`
	// IlluminaProfile is the Illumina run profile of the generated reads
	// (empty: the default).
	IlluminaProfile string `json:"illumina_profile,omitempty"`
	RealFASTQ       string `json:"real_fastq,omitempty"`
	RealFASTQ2      string `json:"real_fastq2,omitempty"`
	RealError       string `json:"real_fastq_error,omitempty"`
	// RealFASTQB, RealFASTQB2 and CrossRun are the second real run and its
	// comparison with the first (-real-fastq-b).
	RealFASTQB  string            `json:"real_fastq_b,omitempty"`
	RealFASTQB2 string            `json:"real_fastq2_b,omitempty"`
	CrossError  string            `json:"cross_run_error,omitempty"`
	CrossRun    *CrossRun         `json:"cross_run,omitempty"`
	Platforms   []*PlatformResult `json:"platforms"`
}

// Failures counts failures across platforms, plus a real FASTQ that could
// not be read and a failed cross-run comparison.
func (r *Report) Failures() int {
	n := 0
	if r.RealError != "" {
		n++
	}
	if r.CrossError != "" {
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

// loadReal reads up to n real reads from path, paired with the reads 2 in
// path2 when it is set.
func loadReal(path, path2 string, n int) (*testing.SyntheticData, error) {
	var d *testing.SyntheticData
	var err error
	if path2 != "" {
		d, err = testing.LoadPairedFASTQ(path, path2, n)
	} else {
		d, err = testing.LoadFASTQ(path, n)
	}
	if err != nil {
		return nil, err
	}
	if len(d.Sequences) < 4 {
		return nil, fmt.Errorf("%s holds %d reads; need at least 4", path, len(d.Sequences))
	}
	return d, nil
}

// runAll measures every platform, and the cross-run comparison.
func runAll(f *testing.E2EFramework, platforms []testing.SequencingPlatform, args *CommandLineArgs) *Report {
	start := time.Now()
	report := &Report{StartedAt: start, Seed: args.seed, Quick: args.quick, ReadsPerSet: args.readsPerSet(),
		Paired: args.pairedRun(), Bases: basesOf(args), IlluminaProfile: args.illuminaProfile, RealFASTQ: args.realFASTQ, RealFASTQ2: args.realFASTQ2,
		RealFASTQB: args.realFASTQB, RealFASTQB2: args.realFASTQB2}

	var real *testing.SyntheticData
	if args.realFASTQ != "" {
		var err error
		if real, err = loadReal(args.realFASTQ, args.realFASTQ2, args.readsPerSet()); err != nil {
			report.RealError = err.Error()
		}
	}
	if real != nil && args.realFASTQB != "" {
		b, err := loadReal(args.realFASTQB, args.realFASTQB2, args.readsPerSet())
		if err == nil {
			report.CrossRun, err = compareRuns(real, b)
		}
		if err != nil {
			report.CrossError = err.Error()
		} else {
			report.CrossRun.FASTQA, report.CrossRun.FASTQB = args.realFASTQ, args.realFASTQB
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
		if c.EnginePlatform == "illumina" {
			c.IlluminaProfile = args.illuminaProfile
		}
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
	// Only the Illumina profiles generate read pairs.
	pairs := args.pairedRun() && config(args.seed, 0).EnginePlatform == "illumina"

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
		ca, cb := config(args.seed+1, n), config(args.seed+2, n)
		ca.PairedEnd, cb.PairedEnd = pairs, pairs
		if c, err := calibrate(f, ca, cb); err != nil {
			fail("calibration", err)
		} else {
			c.Bases = res.Bases
			res.Calibration = c
		}
	}

	if real != nil {
		if real.Paired() && !pairs {
			fail("realism", fmt.Errorf("paired real reads need a paired-end (Illumina) profile"))
			return res
		}
		var synth *testing.SyntheticData
		c := config(args.seed+3, len(real.Sequences))
		if !args.realOnly {
			c.PairedEnd = real.Paired()
			if l := commonLength(real); l > 0 && c.ReadLengthDistribution.Dist == "" {
				c.ReadLength = l
			}
			if synth, err = f.GenerateSyntheticData(c); err != nil {
				fail("realism: generate", err)
				return res
			}
			if args.writeGenerated != "" {
				if err := writeGenerated(f, synth, args.writeGenerated, platform); err != nil {
					fail("write generated reads", err)
				}
			}
		}
		x, err := compareReal(real, synth)
		if err != nil {
			fail("realism", err)
			return res
		}
		x.Bases = res.Bases
		if synth != nil && c.ReadLengthDistribution.Dist == "" {
			x.GeneratedReadLength = c.ReadLength
		}
		res.Realism = x
	}
	return res
}

// writeGenerated writes the generated reads of the real comparison to
// dir/<platform>_1.fastq and, for pairs, the reads 2 to <platform>_2.fastq.
func writeGenerated(f *testing.E2EFramework, d *testing.SyntheticData, dir string, platform testing.SequencingPlatform) error {
	if err := os.MkdirAll(dir, 0750); err != nil {
		return err
	}
	r1 := filepath.Join(dir, string(platform)+"_1.fastq")
	if !d.Paired() {
		return testing.WriteFASTQFile(d, r1)
	}
	return f.WriteFASTQ(d, r1, filepath.Join(dir, string(platform)+"_2.fastq"))
}

// compare returns the holdout comparison of a (class 0) and b (class 1) on
// each of sets, keeping top separating features.
func compare(sets []testing.FeatureSet, a, b *testing.SyntheticData, top int) ([]*Comparison, error) {
	var out []*Comparison
	for _, set := range sets {
		h, err := testing.HoldoutAccuracyFor(set, a, b, top)
		if err != nil {
			return nil, fmt.Errorf("%s: %w", set, err)
		}
		out = append(out, &Comparison{Holdout: h, Band: testing.ChanceBand(h.Holdout, nullBandSigmas)})
	}
	return out, nil
}

// pairSets is the feature-set list of the pair comparisons.
var pairSets = []testing.FeatureSet{testing.FeaturePair}

// compareReal runs the real-versus-real null and, when synth is not nil,
// the real-versus-generated comparison for every feature set: of the reads
// (reads 1), and for paired reads of the reads 2 and the pairs.
func compareReal(real, synth *testing.SyntheticData) (*Realism, error) {
	x := &Realism{RealReads: len(real.Sequences), Paired: real.Paired(),
		RealMeanQuality: real.Stats.MeanQuality, RealGC: real.Stats.GCContent}
	if synth != nil {
		if real.Paired() != synth.Paired() {
			return nil, fmt.Errorf("real reads paired %v, generated paired %v", real.Paired(), synth.Paired())
		}
		x.GeneratedReads = len(synth.Sequences)
		x.GeneratedMeanQuality, x.GeneratedGC = synth.Stats.MeanQuality, synth.Stats.GCContent
	}
	even, odd := testing.SplitAlternate(real)
	// run fills one comparison pair: the null of the real reads and real
	// against generated.
	run := func(what string, sets []testing.FeatureSet, r, e, o, s *testing.SyntheticData) (null, gen []*Comparison, err error) {
		if null, err = compare(sets, e, o, 0); err != nil {
			return nil, nil, fmt.Errorf("%s real vs real: %w", what, err)
		}
		if s != nil {
			if gen, err = compare(sets, r, s, topFeatures); err != nil {
				return nil, nil, fmt.Errorf("%s real vs generated: %w", what, err)
			}
		}
		return null, gen, nil
	}
	var err error
	if x.RealVsReal, x.RealVsGenerated, err = run("read 1", testing.AllFeatureSets, real, even, odd, synth); err != nil {
		return nil, err
	}
	if !real.Paired() {
		return x, nil
	}
	var synth2 *testing.SyntheticData
	if synth != nil {
		synth2 = synth.Mates
	}
	r2 := &MateRealism{RealMeanQuality: real.Mates.Stats.MeanQuality, RealGC: real.Mates.Stats.GCContent}
	if synth2 != nil {
		r2.GeneratedMeanQuality, r2.GeneratedGC = synth2.Stats.MeanQuality, synth2.Stats.GCContent
	}
	if r2.RealVsReal, r2.RealVsGenerated, err = run("read 2", testing.AllFeatureSets, real.Mates, even.Mates, odd.Mates, synth2); err != nil {
		return nil, err
	}
	pair := &MateRealism{}
	if pair.RealVsReal, pair.RealVsGenerated, err = run("pair", pairSets, real, even, odd, synth); err != nil {
		return nil, err
	}
	x.Read2, x.Pair = r2, pair
	return x, nil
}

// compareRuns compares real run a (class 0) with real run b (class 1) on
// every feature set, after cutting both to the shorter run's read length
// when their lengths differ. Both runs must be paired or neither.
func compareRuns(a, b *testing.SyntheticData) (*CrossRun, error) {
	if a.Paired() != b.Paired() {
		return nil, fmt.Errorf("cross-run comparison: run A paired %v, run B paired %v", a.Paired(), b.Paired())
	}
	x := &CrossRun{ReadsA: len(a.Sequences), ReadsB: len(b.Sequences)}
	if la, lb := a.Stats.MaxReadLength, b.Stats.MaxReadLength; la != lb {
		x.ReadLength = min(la, lb)
		a, b = testing.TrimmedTo(a, x.ReadLength), testing.TrimmedTo(b, x.ReadLength)
	}
	x.MeanQualityA, x.MeanQualityB = a.Stats.MeanQuality, b.Stats.MeanQuality
	x.GCA, x.GCB = a.Stats.GCContent, b.Stats.GCContent
	var err error
	if x.Read1, err = compare(testing.AllFeatureSets, a, b, topFeatures); err != nil {
		return nil, fmt.Errorf("cross-run read 1: %w", err)
	}
	if !a.Paired() {
		return x, nil
	}
	if x.Read2, err = compare(testing.AllFeatureSets, a.Mates, b.Mates, topFeatures); err != nil {
		return nil, fmt.Errorf("cross-run read 2: %w", err)
	}
	if x.Pair, err = compare(pairSets, a, b, topFeatures); err != nil {
		return nil, fmt.Errorf("cross-run pair: %w", err)
	}
	return x, nil
}

// calibrate measures, for each feature set, the null comparison (two
// independent draws), the GC control (the second draw GC-skewed to
// gcControlTarget) and the quality control (the second draw's seed with
// the alternative quality model). For paired configs the sets are read
// pairs, compared on the single-read sets by read 1 and on the pair set.
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
	cal := &Calibration{ReadsPerSet: ca.Reads, Paired: ca.PairedEnd, ControlTargetGC: gcControlTarget, QualityControl: describeQualityModel(alt)}
	sets := testing.AllFeatureSets
	if ca.PairedEnd {
		sets = append(append([]testing.FeatureSet(nil), sets...), testing.FeaturePair)
	}
	for _, set := range sets {
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
	if r.IlluminaProfile != "" {
		fmt.Fprintf(w, "Illumina profile: %s\n", r.IlluminaProfile)
	}

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
		real := r.RealFASTQ
		if r.RealFASTQ2 != "" {
			real += " + " + r.RealFASTQ2
		}
		fmt.Fprintf(w, "\nREAL VERSUS GENERATED (%s; generated from %s)\n", real, r.Bases)
		fmt.Fprintln(w, rule)
		if r.RealError != "" {
			fmt.Fprintf(w, "ERROR: %s\n", r.RealError)
		}
		for _, p := range r.Platforms {
			x := p.Realism
			if x == nil {
				continue
			}
			kind := "reads"
			if x.Paired {
				kind = "pairs"
			}
			fmt.Fprintf(w, "%s: %d real %s (mean Q %.2f, GC %.4f), %d generated (mean Q %.2f, GC %.4f",
				p.Platform, x.RealReads, kind, x.RealMeanQuality, x.RealGC, x.GeneratedReads, x.GeneratedMeanQuality, x.GeneratedGC)
			if x.GeneratedReadLength > 0 {
				fmt.Fprintf(w, ", %d bp", x.GeneratedReadLength)
			}
			fmt.Fprintln(w, ")")
			if x.Read2 != nil {
				fmt.Fprintf(w, "  read 2: real mean Q %.2f, GC %.4f; generated mean Q %.2f, GC %.4f\n",
					x.Read2.RealMeanQuality, x.Read2.RealGC, x.Read2.GeneratedMeanQuality, x.Read2.GeneratedGC)
			}
			fmt.Fprintf(w, "  %-4s %-9s %-24s %8s %18s %5s %8s %18s\n", "Read", "Features", "Classifier", "Holdout", "Real vs real", "", "Holdout", "Real vs generated")
			type block struct {
				read      string
				null, gen []*Comparison
			}
			blocks := []block{{"1", x.RealVsReal, x.RealVsGenerated}}
			if x.Read2 != nil {
				blocks = append(blocks, block{"2", x.Read2.RealVsReal, x.Read2.RealVsGenerated})
			}
			if x.Pair != nil {
				blocks = append(blocks, block{"pair", x.Pair.RealVsReal, x.Pair.RealVsGenerated})
			}
			for _, b := range blocks {
				for _, n := range b.null {
					fmt.Fprintf(w, "  %-4s %-9s %-24s %8d %8.4f ± %7.4f %5s", b.read, n.FeatureSet, n.Classifier,
						n.Holdout.Holdout, n.Accuracy, n.Band, passFail(math.Abs(n.Accuracy-0.5) <= n.Band))
					if g := find(b.gen, n.FeatureSet); g != nil {
						fmt.Fprintf(w, " %8d %8.4f ± %7.4f", g.Holdout.Holdout, g.Accuracy, g.Band)
					}
					fmt.Fprintln(w)
				}
			}
			if len(x.RealVsGenerated) == 0 {
				continue
			}
			fmt.Fprintln(w, "  Features separating real (a) from generated (b), by effect size (b−a)/SD on the training reads:")
			for _, b := range blocks {
				for _, set := range []testing.FeatureSet{testing.FeatureSequence, testing.FeatureQuality, testing.FeaturePair} {
					printTop(w, b.read, find(b.gen, set), "real", "generated")
				}
			}
		}
	} else {
		fmt.Fprintln(w, "\nNo real reads supplied (-real-fastq or $PHANTOM_REAL_FASTQ), so no real-versus-generated comparison was made.")
	}

	if r.RealFASTQB != "" {
		fmt.Fprintf(w, "\nREAL RUN A VERSUS REAL RUN B (A %s, B %s)\n", r.RealFASTQ, r.RealFASTQB)
		fmt.Fprintln(w, rule)
		if r.CrossError != "" {
			fmt.Fprintf(w, "ERROR: %s\n", r.CrossError)
		}
		if x := r.CrossRun; x != nil {
			fmt.Fprintf(w, "A: %d reads (mean Q %.2f, GC %.4f); B: %d reads (mean Q %.2f, GC %.4f)",
				x.ReadsA, x.MeanQualityA, x.GCA, x.ReadsB, x.MeanQualityB, x.GCB)
			if x.ReadLength > 0 {
				fmt.Fprintf(w, "; both cut to %d bp", x.ReadLength)
			}
			fmt.Fprintln(w)
			fmt.Fprintf(w, "  %-4s %-9s %-24s %8s %18s\n", "Read", "Features", "Classifier", "Holdout", "A vs B")
			type block struct {
				read string
				cs   []*Comparison
			}
			blocks := []block{{"1", x.Read1}, {"2", x.Read2}, {"pair", x.Pair}}
			for _, b := range blocks {
				for _, c := range b.cs {
					fmt.Fprintf(w, "  %-4s %-9s %-24s %8d %8.4f ± %7.4f\n", b.read, c.FeatureSet, c.Classifier, c.Holdout.Holdout, c.Accuracy, c.Band)
				}
			}
			fmt.Fprintln(w, "  Features separating run A (a) from run B (b), by effect size (b−a)/SD on the training reads:")
			for _, b := range blocks {
				for _, set := range []testing.FeatureSet{testing.FeatureSequence, testing.FeatureQuality, testing.FeaturePair} {
					printTop(w, b.read, find(b.cs, set), "A", "B")
				}
			}
		}
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

// printTop prints the first topFeaturesPrinted separating features of c.
func printTop(w io.Writer, read string, c *Comparison, a, b string) {
	if c == nil {
		return
	}
	for i, d := range c.Top {
		if i == topFeaturesPrinted {
			break
		}
		fmt.Fprintf(w, "    %-4s %-9s %-34s %s %10.5f  %s %10.5f  d %+7.3f\n", read, c.FeatureSet, d.Name, a, d.MeanA, b, d.MeanB, d.EffectSize)
	}
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
	for _, set := range testing.AllFeatureSets {
		h = append(h, "r2_real_vs_real_accuracy_"+string(set), "r2_real_vs_generated_accuracy_"+string(set))
	}
	h = append(h, "null_accuracy_pair", "quality_control_accuracy_pair", "real_vs_real_accuracy_pair", "real_vs_generated_accuracy_pair")
	return append(h, "errors")
}

// mateCells returns the read-2 and pair cells of p's CSV row.
func mateCells(p *PlatformResult, ff func(float64) string) []string {
	acc := func(cs []*Comparison, set testing.FeatureSet) string {
		if c := find(cs, set); c != nil {
			return ff(c.Accuracy)
		}
		return ""
	}
	var cells []string
	x := p.Realism
	for _, set := range testing.AllFeatureSets {
		if x != nil && x.Read2 != nil {
			cells = append(cells, acc(x.Read2.RealVsReal, set), acc(x.Read2.RealVsGenerated, set))
		} else {
			cells = append(cells, "", "")
		}
	}
	pc := make([]string, 4)
	if c := p.Calibration; c != nil {
		if s := c.Set(testing.FeaturePair); s != nil {
			pc[0], pc[1] = ff(s.NullAccuracy), ff(s.QualityControlAccuracy)
		}
	}
	if x != nil && x.Pair != nil {
		pc[2], pc[3] = acc(x.Pair.RealVsReal, testing.FeaturePair), acc(x.Pair.RealVsGenerated, testing.FeaturePair)
	}
	return append(cells, pc...)
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
		row = append(row, mateCells(p, ff)...)
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
