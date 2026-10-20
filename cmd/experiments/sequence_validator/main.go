// Command sequence_validator is a developer tool, not part of the
// phantom-science CLI. It checks single-record FASTA files named
// <POP>_L<length>_S<n>.fasta (or .fa) found under the paths given on the
// command line, and writes per-file results plus per-population and
// per-region summaries as JSON.
//
// For each file it reports the length against the L<length> in the name,
// bases other than ACGT, GC content of the ACGT bases against the header's
// model_target_gc=<x> (when the header has one), distinct-8-mer fraction and
// the longest homopolymer, and flags values outside the thresholds set by
// the flags. These are format and summary-statistic checks; they say nothing
// about biological realism. Files that cannot be parsed are listed in the
// report as skipped.
//
// Usage:
//
//	sequence_validator [flags] PATH...
//
// Each PATH is a FASTA file or a directory searched recursively. Exit status
// is 0 when every sequence passes, 2 when any fails, and 1 on an error or
// when no sequences are found.
package main

import (
	"bufio"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"time"
)

// complexityK is the k-mer size of the complexity measure.
const complexityK = 8

// thresholds are the limits a sequence must meet to pass.
type thresholds struct {
	GCTolerance    float64 `json:"gc_tolerance"`    // max |GC - model_target_gc|
	MinComplexity  float64 `json:"min_complexity"`  // min distinct-8-mer fraction
	MaxHomopolymer int     `json:"max_homopolymer"` // longest allowed single-base run
}

var defaultThresholds = thresholds{GCTolerance: 0.05, MinComplexity: 0.3, MaxHomopolymer: 10}

// ValidationResult is the result for one FASTA file.
type ValidationResult struct {
	Population     string `json:"population"`
	Region         string `json:"region"`
	Length         int    `json:"length"` // expected, from the file name
	SequenceNumber int    `json:"sequence_number"`
	Filename       string `json:"filename"`

	ActualLength int     `json:"actual_length"`
	GCContent    float64 `json:"gc_content"` // of the ACGT bases
	Complexity   float64 `json:"complexity"` // distinct / total ACGT-only 8-mers
	ValidBases   int     `json:"valid_bases"`
	InvalidBases int     `json:"invalid_bases"`

	MaxHomopolymer int `json:"max_homopolymer"`

	HasTargetGC     bool    `json:"has_target_gc"`
	TargetGC        float64 `json:"target_gc,omitempty"`
	GCDeviation     float64 `json:"gc_deviation,omitempty"`
	ComplexityScore string  `json:"complexity_score"`

	IsValid         bool     `json:"is_valid"`
	ValidationScore float64  `json:"validation_score"`
	Issues          []string `json:"issues"`
}

// PopulationSummary summarises the results for one population.
type PopulationSummary struct {
	Population     string         `json:"population"`
	Region         string         `json:"region"`
	TotalSequences int            `json:"total_sequences"`
	ValidSequences int            `json:"valid_sequences"`
	AvgGCContent   float64        `json:"avg_gc_content"`
	AvgComplexity  float64        `json:"avg_complexity"`
	AvgScore       float64        `json:"avg_score"`
	Issues         map[string]int `json:"issues"`
}

// SkippedFile is a FASTA file that could not be validated.
type SkippedFile struct {
	Path   string `json:"path"`
	Reason string `json:"reason"`
}

// ValidationReport is the JSON report.
type ValidationReport struct {
	TotalSequences     int                 `json:"total_sequences"`
	ValidSequences     int                 `json:"valid_sequences"`
	OverallScore       float64             `json:"overall_score"`
	PopulationResults  []ValidationResult  `json:"population_results"`
	PopulationSummary  []PopulationSummary `json:"population_summary"`
	RegionalAnalysis   map[string]float64  `json:"regional_analysis"`
	Skipped            []SkippedFile       `json:"skipped,omitempty"`
	Thresholds         thresholds          `json:"thresholds"`
	ValidationCriteria string              `json:"validation_criteria"`
	Timestamp          string              `json:"timestamp"`
}

const validationCriteria = "Format and summary-statistic checks: length vs file name, non-ACGT bases, " +
	"GC vs model_target_gc, distinct 8-mer fraction, longest homopolymer. Not a test of biological realism."

func main() {
	os.Exit(run(os.Args[1:], os.Stdout, os.Stderr))
}

func run(args []string, stdout, stderr io.Writer) int {
	fl := flag.NewFlagSet("sequence_validator", flag.ContinueOnError)
	fl.SetOutput(stderr)
	out := fl.String("out", "validation_report.json", "JSON report path")
	th := defaultThresholds
	fl.Float64Var(&th.GCTolerance, "gc-tolerance", th.GCTolerance, "max |GC - model_target_gc|")
	fl.Float64Var(&th.MinComplexity, "min-complexity", th.MinComplexity, "min distinct-8-mer fraction")
	fl.IntVar(&th.MaxHomopolymer, "max-homopolymer", th.MaxHomopolymer, "longest allowed single-base run")
	fl.Usage = func() {
		fmt.Fprintf(stderr, "usage: sequence_validator [flags] PATH...\n\n"+
			"Validates <POP>_L<length>_S<n>.fasta files; each PATH is a file or a directory.\n\n")
		fl.PrintDefaults()
	}
	if err := fl.Parse(args); err != nil {
		return 1
	}
	if fl.NArg() == 0 {
		fl.Usage()
		return 1
	}

	results, skipped, err := validatePaths(fl.Args(), th)
	if err != nil {
		fmt.Fprintf(stderr, "sequence_validator: %v\n", err)
		return 1
	}
	for _, s := range skipped {
		fmt.Fprintf(stderr, "skipped %s: %s\n", s.Path, s.Reason)
	}

	report := generateValidationReport(results, th)
	report.Skipped = skipped
	displayValidationResults(stdout, report)

	if err := saveValidationReport(report, *out); err != nil {
		fmt.Fprintf(stderr, "sequence_validator: write report: %v\n", err)
		return 1
	}
	fmt.Fprintf(stdout, "\nreport written to %s\n", *out)

	switch {
	case report.TotalSequences == 0:
		fmt.Fprintf(stderr, "sequence_validator: no sequences found in %s\n", strings.Join(fl.Args(), " "))
		return 1
	case report.ValidSequences < report.TotalSequences:
		return 2
	}
	return 0
}

func isFASTAName(name string) bool {
	ext := strings.ToLower(filepath.Ext(name))
	return ext == ".fasta" || ext == ".fa"
}

// validatePaths validates every FASTA file under paths. Files that cannot be
// parsed are returned as skipped; an unreadable path is an error.
func validatePaths(paths []string, th thresholds) ([]ValidationResult, []SkippedFile, error) {
	var results []ValidationResult
	var skipped []SkippedFile
	for _, root := range paths {
		err := filepath.WalkDir(root, func(path string, d fs.DirEntry, err error) error {
			if err != nil {
				return err
			}
			if d.IsDir() || !isFASTAName(d.Name()) {
				return nil
			}
			r, verr := validateFile(path, th)
			if verr != nil {
				skipped = append(skipped, SkippedFile{Path: path, Reason: verr.Error()})
				return nil
			}
			results = append(results, *r)
			return nil
		})
		if err != nil {
			return nil, nil, err
		}
	}
	return results, skipped, nil
}

func validateFile(path string, th thresholds) (*ValidationResult, error) {
	population, length, seqNum, err := parseFilename(filepath.Base(path))
	if err != nil {
		return nil, err
	}
	return validateSequenceFile(path, population, length, seqNum, th)
}

// parseFilename parses <POP>_L<length>_S<n>.fasta (or .fa). POP may itself
// contain underscores.
func parseFilename(filename string) (population string, length int, seqNum int, err error) {
	name := strings.TrimSuffix(filename, filepath.Ext(filename))
	parts := strings.Split(name, "_")
	if len(parts) < 3 {
		return "", 0, 0, fmt.Errorf("file name %q is not <POP>_L<length>_S<n>.fasta", filename)
	}
	n := len(parts)
	population = strings.Join(parts[:n-2], "_")
	lengthStr, okL := strings.CutPrefix(parts[n-2], "L")
	seqStr, okS := strings.CutPrefix(parts[n-1], "S")
	if population == "" || !okL || !okS {
		return "", 0, 0, fmt.Errorf("file name %q is not <POP>_L<length>_S<n>.fasta", filename)
	}
	length, err = strconv.Atoi(lengthStr)
	if err != nil || length <= 0 {
		return "", 0, 0, fmt.Errorf("invalid length %q in file name %q", parts[n-2], filename)
	}
	seqNum, err = strconv.Atoi(seqStr)
	if err != nil {
		return "", 0, 0, fmt.Errorf("invalid sequence number %q in file name %q", parts[n-1], filename)
	}
	return population, length, seqNum, nil
}

// readSingleFASTA returns the header (without '>') and sequence of a
// one-record FASTA file.
func readSingleFASTA(r io.Reader) (header, sequence string, err error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 1<<20), 1<<30) // sequences may be on one long line
	if !sc.Scan() {
		if err := sc.Err(); err != nil {
			return "", "", err
		}
		return "", "", errors.New("empty file")
	}
	header = sc.Text()
	if !strings.HasPrefix(header, ">") {
		return "", "", errors.New("not FASTA: first line does not start with '>'")
	}
	var b strings.Builder
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if strings.HasPrefix(line, ">") {
			return "", "", errors.New("more than one record; expected one sequence per file")
		}
		b.WriteString(line)
	}
	if err := sc.Err(); err != nil {
		return "", "", err
	}
	if b.Len() == 0 {
		return "", "", errors.New("empty sequence")
	}
	return header[1:], b.String(), nil
}

func validateSequenceFile(path, population string, expectedLength, seqNum int, th thresholds) (*ValidationResult, error) {
	f, err := os.Open(path) // #nosec G304 -- path comes from the user's command line
	if err != nil {
		return nil, err
	}
	defer f.Close()
	header, seq, err := readSingleFASTA(f)
	if err != nil {
		return nil, err
	}
	targetGC, hasTarget, region := parseHeaderMetadata(header)
	r := validateSequence(seq, population, region, expectedLength, seqNum, path, th)
	if hasTarget {
		r.HasTargetGC, r.TargetGC = true, targetGC
		r.GCDeviation = r.GCContent - targetGC
		if abs(r.GCDeviation) > th.GCTolerance {
			r.Issues = append(r.Issues, fmt.Sprintf("GC content %.3f differs from target %.3f by more than %.3f",
				r.GCContent, targetGC, th.GCTolerance))
		}
	}
	finalise(r)
	return r, nil
}

// parseHeaderMetadata reads region=<name> and model_target_gc=<x> from a
// FASTA header such as "POP_L1000_S01 population=POP region=Africa
// model_target_gc=0.496". hasTarget is false when model_target_gc is
// missing or not a number.
func parseHeaderMetadata(header string) (targetGC float64, hasTarget bool, region string) {
	for _, field := range strings.Fields(header) {
		if v, ok := strings.CutPrefix(field, "region="); ok {
			region = v
		} else if v, ok := strings.CutPrefix(field, "model_target_gc="); ok {
			if gc, err := strconv.ParseFloat(v, 64); err == nil {
				targetGC, hasTarget = gc, true
			}
		}
	}
	return targetGC, hasTarget, region
}

// validateSequence computes the sequence statistics and the length, base,
// complexity and homopolymer issues. The GC-target check needs the header
// and is added by the caller; call finalise afterwards.
func validateSequence(sequence, population, region string, expectedLength, seqNum int, filename string, th thresholds) *ValidationResult {
	seq := strings.ToUpper(sequence)
	r := &ValidationResult{
		Population:     population,
		Region:         region,
		Length:         expectedLength,
		SequenceNumber: seqNum,
		Filename:       filepath.Base(filename),
		Issues:         []string{},
	}
	r.ActualLength = len(seq)
	r.ValidBases, r.InvalidBases = countValidBases(seq)
	r.GCContent = calculateGCContent(seq)
	r.Complexity = calculateComplexity(seq, complexityK)
	r.MaxHomopolymer = findMaxHomopolymer(seq)
	r.ComplexityScore = classifyComplexity(r.Complexity)

	if r.ActualLength != expectedLength {
		r.Issues = append(r.Issues, fmt.Sprintf("length %d, file name says %d", r.ActualLength, expectedLength))
	}
	if r.InvalidBases > 0 {
		r.Issues = append(r.Issues, fmt.Sprintf("%d non-ACGT bases", r.InvalidBases))
	}
	if r.Complexity < th.MinComplexity {
		r.Issues = append(r.Issues, fmt.Sprintf("complexity %.4f below %.4f", r.Complexity, th.MinComplexity))
	}
	if r.MaxHomopolymer > th.MaxHomopolymer {
		r.Issues = append(r.Issues, fmt.Sprintf("homopolymer run of %d bases exceeds %d", r.MaxHomopolymer, th.MaxHomopolymer))
	}
	return r
}

// finalise sets the score and pass/fail once all issues are recorded.
func finalise(r *ValidationResult) {
	r.ValidationScore = calculateValidationScore(r)
	r.IsValid = len(r.Issues) == 0
}

func isACGT(b byte) bool {
	return b == 'A' || b == 'C' || b == 'G' || b == 'T' || b == 'a' || b == 'c' || b == 'g' || b == 't'
}

func countValidBases(sequence string) (valid, invalid int) {
	for i := 0; i < len(sequence); i++ {
		if isACGT(sequence[i]) {
			valid++
		} else {
			invalid++
		}
	}
	return valid, invalid
}

// calculateGCContent returns the G+C fraction of the ACGT bases, ignoring
// N and other symbols; 0 when there are none.
func calculateGCContent(sequence string) float64 {
	gc, acgt := 0, 0
	for i := 0; i < len(sequence); i++ {
		switch sequence[i] {
		case 'G', 'C', 'g', 'c':
			gc++
			acgt++
		case 'A', 'T', 'a', 't':
			acgt++
		}
	}
	if acgt == 0 {
		return 0
	}
	return float64(gc) / float64(acgt)
}

// calculateComplexity returns distinct / total over the k-mers made only of
// ACGT bases (case-insensitive); 0 when there are none.
func calculateComplexity(sequence string, k int) float64 {
	seq := strings.ToUpper(sequence)
	seen := make(map[string]struct{})
	total, run := 0, 0 // run: ACGT bases ending at i
	for i := 0; i < len(seq); i++ {
		if !isACGT(seq[i]) {
			run = 0
			continue
		}
		run++
		if run >= k {
			seen[seq[i-k+1:i+1]] = struct{}{}
			total++
		}
	}
	if total == 0 {
		return 0
	}
	return float64(len(seen)) / float64(total)
}

// findMaxHomopolymer returns the longest run of one repeated base,
// case-insensitively.
func findMaxHomopolymer(sequence string) int {
	if len(sequence) == 0 {
		return 0
	}
	seq := strings.ToUpper(sequence)
	maxRun, cur := 1, 1
	for i := 1; i < len(seq); i++ {
		if seq[i] == seq[i-1] {
			cur++
			maxRun = max(maxRun, cur)
		} else {
			cur = 1
		}
	}
	return maxRun
}

func classifyComplexity(complexity float64) string {
	switch {
	case complexity >= 0.8:
		return "Very High"
	case complexity >= 0.6:
		return "High"
	case complexity >= 0.4:
		return "Medium"
	case complexity >= 0.2:
		return "Low"
	default:
		return "Very Low"
	}
}

// calculateValidationScore is a heuristic summary in [0, 1]: start at 1,
// subtract 0.1 per issue, 0.3 below 0.3 complexity, 2x the GC deviation
// (when there is a target) and 0.5 for any non-ACGT base; add 0.2 at or
// above 0.6 complexity.
func calculateValidationScore(r *ValidationResult) float64 {
	score := 1.0 - 0.1*float64(len(r.Issues))
	if r.Complexity >= 0.6 {
		score += 0.2
	} else if r.Complexity < 0.3 {
		score -= 0.3
	}
	if r.HasTargetGC {
		score -= 2 * abs(r.GCDeviation)
	}
	if r.InvalidBases > 0 {
		score -= 0.5
	}
	return min(max(score, 0), 1)
}

func generateValidationReport(results []ValidationResult, th thresholds) *ValidationReport {
	report := &ValidationReport{
		TotalSequences:     len(results),
		PopulationResults:  results,
		PopulationSummary:  []PopulationSummary{},
		RegionalAnalysis:   make(map[string]float64),
		Thresholds:         th,
		ValidationCriteria: validationCriteria,
		Timestamp:          time.Now().UTC().Format(time.RFC3339),
	}
	if report.PopulationResults == nil {
		report.PopulationResults = []ValidationResult{}
	}

	total := 0.0
	popStats := make(map[string]*PopulationSummary)
	regionScores := make(map[string][]float64)
	for _, r := range results {
		if r.IsValid {
			report.ValidSequences++
		}
		total += r.ValidationScore

		s, ok := popStats[r.Population]
		if !ok {
			s = &PopulationSummary{Population: r.Population, Region: r.Region, Issues: make(map[string]int)}
			popStats[r.Population] = s
		}
		s.TotalSequences++
		if r.IsValid {
			s.ValidSequences++
		}
		s.AvgGCContent += r.GCContent
		s.AvgComplexity += r.Complexity
		s.AvgScore += r.ValidationScore
		for _, issue := range r.Issues {
			s.Issues[issue]++
		}
		regionScores[r.Region] = append(regionScores[r.Region], r.ValidationScore)
	}
	if len(results) > 0 {
		report.OverallScore = total / float64(len(results))
	}

	for _, s := range popStats {
		n := float64(s.TotalSequences)
		s.AvgGCContent /= n
		s.AvgComplexity /= n
		s.AvgScore /= n
		report.PopulationSummary = append(report.PopulationSummary, *s)
	}
	sort.Slice(report.PopulationSummary, func(i, j int) bool {
		a, b := report.PopulationSummary[i], report.PopulationSummary[j]
		if a.AvgScore != b.AvgScore {
			return a.AvgScore > b.AvgScore
		}
		return a.Population < b.Population
	})

	for region, scores := range regionScores {
		sum := 0.0
		for _, s := range scores {
			sum += s
		}
		report.RegionalAnalysis[region] = sum / float64(len(scores))
	}
	return report
}

func percent(n, d int) float64 {
	if d == 0 {
		return 0
	}
	return 100 * float64(n) / float64(d)
}

func displayValidationResults(w io.Writer, report *ValidationReport) {
	fmt.Fprintf(w, "Valid sequences: %d/%d (%.1f%%), mean score %.3f\n",
		report.ValidSequences, report.TotalSequences,
		percent(report.ValidSequences, report.TotalSequences), report.OverallScore)
	if len(report.Skipped) > 0 {
		fmt.Fprintf(w, "Skipped files: %d\n", len(report.Skipped))
	}

	if len(report.PopulationSummary) > 0 {
		fmt.Fprintf(w, "\nBy population:\n")
	}
	for _, s := range report.PopulationSummary {
		region := s.Region
		if region == "" {
			region = "no region"
		}
		fmt.Fprintf(w, "  %s (%s): valid %d/%d (%.1f%%), score %.3f, GC %.3f, complexity %.3f\n",
			s.Population, region, s.ValidSequences, s.TotalSequences,
			percent(s.ValidSequences, s.TotalSequences), s.AvgScore, s.AvgGCContent, s.AvgComplexity)
		issues := make([]string, 0, len(s.Issues))
		for issue := range s.Issues {
			issues = append(issues, issue)
		}
		sort.Slice(issues, func(i, j int) bool {
			if s.Issues[issues[i]] != s.Issues[issues[j]] {
				return s.Issues[issues[i]] > s.Issues[issues[j]]
			}
			return issues[i] < issues[j]
		})
		for _, issue := range issues {
			fmt.Fprintf(w, "    %dx %s\n", s.Issues[issue], issue)
		}
	}

	if len(report.RegionalAnalysis) > 0 {
		regions := make([]string, 0, len(report.RegionalAnalysis))
		for r := range report.RegionalAnalysis {
			regions = append(regions, r)
		}
		sort.Strings(regions)
		fmt.Fprintf(w, "\nMean score by region:\n")
		for _, r := range regions {
			name := r
			if name == "" {
				name = "no region"
			}
			fmt.Fprintf(w, "  %s: %.3f\n", name, report.RegionalAnalysis[r])
		}
	}
}

func saveValidationReport(report *ValidationReport, filename string) error {
	data, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(filename, data, 0o600)
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
